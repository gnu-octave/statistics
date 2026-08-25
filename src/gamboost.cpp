/*
Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This file is part of the statistics package for GNU Octave.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

// The boosted-tree engine shared by the GAM learners, beside the spline engine
// in gam.cpp.  A generalized additive model is additive by construction, so
// every tree of the predictor phase splits on one predictor only.  That is
// what makes this small: no multivariate tree is needed, no surrogate splits,
// no pruning, and the whole fitted shape function of a predictor is a step
// function over that predictor's bins, whatever number of trees produced it.
//
// The scheme is Newton boosting on the deviance.  Each round fits one tree per
// predictor to the current gradient and Hessian, adds it at the running step
// size, and recentres: a shape function is held to mean zero and the constant
// it gives up is absorbed into the intercept, which is the usual GAM
// identifiability convention and is why a fitted intercept moves away from the
// value it was seeded with.
//
// The step is not fixed.  When a round fails to reduce the deviance by enough
// it is retried at half the step, and the fit stops when even the reduced step
// cannot improve it.  That is what MATLAB reports as "Unable to improve the
// model fit" against "Terminated after training the requested number of
// trees", measured on R2024a.  The constants below are ours: the tolerance
// MATLAB stops on is not recoverable from anything it reports, so this engine
// documents its own rather than pretending to match one it cannot read.

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>

// The largest number of cut points a predictor is binned at.  MATLAB reports
// 255 edges for any predictor carrying more distinct values than that
// (measured at 600 and at 2000), so a predictor with fewer distinct values is
// split at every midpoint between them and a richer one at 255 quantiles.
static const octave_idx_type GAMB_MAX_EDGES = 255;

// How many rounds are allowed to pass without a meaningful improvement before
// the fit is declared unable to improve, and what counts as meaningful,
// relative to the deviance the window opened at.  This is the patience form
// scikit-learn and the Explainable Boosting Machine both use; it is measured
// on the training deviance rather than on a held-out split, because a random
// split would make the fit irreproducible and MATLAB's is deterministic.  The
// tolerance is relative and not absolute because this deviance spans fifteen
// orders of magnitude over a run, where an absolute one would be meaningless.
static const octave_idx_type GAMB_PATIENCE = 10;
static const double GAMB_REL_TOL = 1e-6;

// How many times a round may halve its step before it gives up on improving at
// all.  MATLAB reduces its step the same way, which is why its learn-rate
// argument is named an *initial* rate; three is ours, its limit not being
// recoverable from anything it reports.
static const int GAMB_MAX_HALVINGS = 3;

// The fewest observations a predictor tree's leaf may hold.  Measured against
// R2024a on 2026-08-25: over a fixture whose only useful split isolates the
// top K points, MATLAB splits there for K of 5 and above and falls back to the
// cut leaving 5 for every K below it, on fitcgam and fitrgam alike, with bin
// edges identical to ours.  It is fitrtree's own default, which is the tree a
// GAM fits internally.  The interaction fitter is not measured and does not
// apply this.
static const octave_idx_type GAMB_MIN_LEAF = 5;

// One predictor reduced to bin indices, with the cut points that produced
// them.  A missing value gets bin -1 and takes no part in any split.
struct BinnedPredictor
{
  RowVector edges;                  // nbins - 1 cut points, ascending
  Array<octave_idx_type> bin;       // one per observation, -1 if missing
  octave_idx_type nbins;
};

// Bin one predictor.  The cut points are the midpoints between consecutive
// distinct finite values, which is every split a tree could make; where that
// would exceed GAMB_MAX_EDGES the midpoints are thinned to that many, taken at
// equally spaced positions through the sorted distinct values so the cuts stay
// where the data is rather than where its range is.
// The detection grid is coarser than the fitting grid and its size is fixed:
// MATLAB reports 7 cut points, 8 equal-frequency bins, at 60, 250, 1000 and
// 4000 observations alike, and they sit on the octiles to every digit.
static const octave_idx_type GAMB_PAIR_EDGES = 7;

static BinnedPredictor
gamb_bin (const ColumnVector& x, octave_idx_type maxedges)
{
  octave_idx_type n = x.numel ();

  std::vector<double> v;
  v.reserve (n);
  for (octave_idx_type i = 0; i < n; i++)
  {
    if (! octave::math::isnan (x(i)))
    {
      v.push_back (x(i));
    }
  }

  std::sort (v.begin (), v.end ());
  std::vector<double> sorted (v);          // the data, still with its repeats
  v.erase (std::unique (v.begin (), v.end ()), v.end ());

  octave_idx_type nd = (octave_idx_type) v.size ();
  octave_idx_type nf = (octave_idx_type) sorted.size ();
  BinnedPredictor B;

  // A constant predictor admits no split at all: one bin, no cut points.
  if (nd < 2)
  {
    B.edges = RowVector (0);
    B.nbins = 1;
    B.bin = Array<octave_idx_type> (dim_vector (n, 1), 0);
    for (octave_idx_type i = 0; i < n; i++)
    {
      B.bin(i) = octave::math::isnan (x(i)) ? -1 : 0;
    }
    return B;
  }

  octave_idx_type ne = nd - 1;
  if (ne > maxedges)
  {
    ne = maxedges;
  }

  B.edges = RowVector (ne);
  if (ne == nd - 1)
  {
    // Below the cap every midpoint between consecutive distinct values is a
    // cut, so the bins are as fine as any split could be.  MATLAB, scikit-learn
    // and the Explainable Boosting Machine all agree on this rule to the last
    // digit, and it is the case almost every predictor falls into.
    for (octave_idx_type k = 0; k < ne; k++)
    {
      B.edges(k) = 0.5 * (v[k] + v[k + 1]);
    }
  }
  else
  {
    // Above the cap the cuts are equally spaced through the observations
    // rather than through the distinct values, so the bins carry equal counts.
    // Each one is placed midway between the two order statistics that bracket
    // its position, which keeps the convention the same as below the cap.
    // MATLAB's grid is equal-frequency too but its interpolation differs:
    // measured on 1000 standard normal draws the two agree to 0.067 at worst
    // against 0.78 for a spread through the distinct values, so this is the
    // right family and not the exact member.
    for (octave_idx_type k = 0; k < ne; k++)
    {
      double r = (double) (k + 1) * nf / (double) (ne + 1);
      octave_idx_type lo = (octave_idx_type) std::floor (r) - 1;
      if (lo < 0)
      {
        lo = 0;
      }
      if (lo > nf - 2)
      {
        lo = nf - 2;
      }
      B.edges(k) = 0.5 * (sorted[lo] + sorted[lo + 1]);
    }
  }

  // Thinning can put two cuts at the same place when the distinct values are
  // very unevenly spaced; drop the repeats so every bin can be reached.
  octave_idx_type keep = 0;
  for (octave_idx_type k = 0; k < ne; k++)
  {
    if (keep == 0 || B.edges(k) > B.edges(keep - 1))
    {
      B.edges(keep++) = B.edges(k);
    }
  }
  if (keep < ne)
  {
    RowVector e (keep);
    for (octave_idx_type k = 0; k < keep; k++)
    {
      e(k) = B.edges(k);
    }
    B.edges = e;
    ne = keep;
  }

  B.nbins = ne + 1;
  B.bin = Array<octave_idx_type> (dim_vector (n, 1), 0);
  for (octave_idx_type i = 0; i < n; i++)
  {
    if (octave::math::isnan (x(i)))
    {
      B.bin(i) = -1;
      continue;
    }
    // The bin is the count of cut points at or below the value: a binary
    // search, since the cut points are ascending.
    octave_idx_type lo = 0;
    octave_idx_type hi = ne;
    while (lo < hi)
    {
      octave_idx_type mid = lo + (hi - lo) / 2;
      if (x(i) > B.edges(mid))
      {
        lo = mid + 1;
      }
      else
      {
        hi = mid;
      }
    }
    B.bin(i) = lo;
  }

  return B;
}

// One contiguous run of bins, as the tree grows it.
struct BinRegion
{
  octave_idx_type lo;               // first bin, inclusive
  octave_idx_type hi;               // last bin, inclusive
  double gain;                      // best gain available from splitting it
  octave_idx_type cut;              // bin the split would end the left side at
};

// The Newton gain of splitting a region whose gradient and Hessian totals are
// (GL, HL) on the left and (GR, HR) on the right.  A leaf's value is -G/H, so
// this is the deviance drop the split buys.
static inline double
gamb_gain (double GL, double HL, double GR, double HR)
{
  const double eps = 1e-12;
  if (HL <= eps || HR <= eps)
  {
    return -1.0;
  }
  return GL * GL / HL + GR * GR / HR - (GL + GR) * (GL + GR) / (HL + HR);
}

// Find the best cut inside a region, given the prefix sums of the gradient, the
// Hessian and the observation count over bins.  Returns the gain and sets CUT
// to the last bin of the left side; a region that cannot be split usefully
// returns a gain of -1.
//
// A cut leaving fewer than GAMB_MIN_LEAF observations on either side is not
// considered, whatever it gains.  The Hessian guard in gamb_gain cannot stand
// in for this: a logistic Hessian is p*(1-p) and says nothing about how many
// observations produced it, so a single well separated point can carry enough
// curvature to look like a leaf worth having.
static double
gamb_best_cut (const std::vector<double>& G, const std::vector<double>& H,
               const std::vector<double>& C, octave_idx_type lo,
               octave_idx_type hi, octave_idx_type& cut)
{
  double best = -1.0;
  cut = -1;

  double Gtot = G[hi + 1] - G[lo];
  double Htot = H[hi + 1] - H[lo];
  double Ctot = C[hi + 1] - C[lo];

  for (octave_idx_type b = lo; b < hi; b++)
  {
    double CL = C[b + 1] - C[lo];
    if (CL < GAMB_MIN_LEAF || Ctot - CL < GAMB_MIN_LEAF)
    {
      continue;
    }
    double GL = G[b + 1] - G[lo];
    double HL = H[b + 1] - H[lo];
    double g = gamb_gain (GL, HL, Gtot - GL, Htot - HL);
    if (g > best)
    {
      best = g;
      cut = b;
    }
  }

  return best;
}

// Fit one tree to a single binned predictor and accumulate its leaf values,
// scaled by the step, into VAL.  The tree is grown best-first: the region
// offering the largest gain is split, up to MAXSPLITS splits in all, which is
// what MATLAB's MaxNumSplitsPerPredictor counts.  With the default of 1 this
// is a stump.
static void
gamb_fit_tree (const BinnedPredictor& B, const ColumnVector& grad,
               const ColumnVector& hess, octave_idx_type maxsplits,
               double step, ColumnVector& val)
{
  octave_idx_type n = grad.numel ();
  octave_idx_type nb = B.nbins;

  // Gradient and Hessian totals per bin, then their prefix sums, so any
  // region's totals are one subtraction.
  std::vector<double> G (nb + 1, 0.0);
  std::vector<double> H (nb + 1, 0.0);
  std::vector<double> C (nb + 1, 0.0);
  for (octave_idx_type i = 0; i < n; i++)
  {
    octave_idx_type b = B.bin(i);
    if (b < 0)
    {
      continue;
    }
    G[b + 1] += grad(i);
    H[b + 1] += hess(i);
    C[b + 1] += 1.0;
  }
  for (octave_idx_type b = 0; b < nb; b++)
  {
    G[b + 1] += G[b];
    H[b + 1] += H[b];
    C[b + 1] += C[b];
  }

  std::vector<BinRegion> leaves;
  BinRegion root;
  root.lo = 0;
  root.hi = nb - 1;
  root.gain = gamb_best_cut (G, H, C, root.lo, root.hi, root.cut);
  leaves.push_back (root);

  for (octave_idx_type s = 0; s < maxsplits; s++)
  {
    // The leaf that buys the most.
    std::size_t pick = 0;
    double best = -1.0;
    for (std::size_t k = 0; k < leaves.size (); k++)
    {
      if (leaves[k].gain > best)
      {
        best = leaves[k].gain;
        pick = k;
      }
    }
    if (best <= 0.0)
    {
      break;
    }

    BinRegion left, right;
    left.lo = leaves[pick].lo;
    left.hi = leaves[pick].cut;
    right.lo = leaves[pick].cut + 1;
    right.hi = leaves[pick].hi;
    left.gain = gamb_best_cut (G, H, C, left.lo, left.hi, left.cut);
    right.gain = gamb_best_cut (G, H, C, right.lo, right.hi, right.cut);

    leaves[pick] = left;
    leaves.push_back (right);
  }

  // Every leaf contributes its Newton step to the bins it covers.
  for (std::size_t k = 0; k < leaves.size (); k++)
  {
    double Gk = G[leaves[k].hi + 1] - G[leaves[k].lo];
    double Hk = H[leaves[k].hi + 1] - H[leaves[k].lo];
    if (Hk <= 1e-12)
    {
      continue;
    }
    double leaf = step * Gk / Hk;
    for (octave_idx_type b = leaves[k].lo; b <= leaves[k].hi; b++)
    {
      val(b) += leaf;
    }
  }
}

// What a boosted fit produces, in the shape the learners store it.  A shape
// function is a step function: VALUE[j](b) is what predictor j contributes for
// any observation falling in its bin b, so prediction is a lookup and does not
// depend on how many trees built it.
struct GamBoostFit
{
  double intercept;
  std::vector<RowVector> edges;
  std::vector<ColumnVector> value;
  octave_idx_type ntrees;
  std::string reason;
  double deviance;
  ColumnVector residuals;
};

// The bin a value falls in, for a predictor already binned at EDGES.  Matches
// the binary search gamb_bin used, so a training row bins the same way at
// predict time; a missing value has no bin and is signalled by -1.
static octave_idx_type
gamb_bin_of (const RowVector& edges, double v)
{
  if (octave::math::isnan (v))
  {
    return -1;
  }
  octave_idx_type lo = 0;
  octave_idx_type hi = edges.numel ();
  while (lo < hi)
  {
    octave_idx_type mid = lo + (hi - lo) / 2;
    if (v > edges(mid))
    {
      lo = mid + 1;
    }
    else
    {
      hi = mid;
    }
  }
  return lo;
}

// The additive prediction of a fitted model at one row of X.  A predictor that
// is missing contributes nothing rather than poisoning the sum, which is what
// a tree does with a missing value it was never given a surrogate for.
static double
gamb_predict_row (const GamBoostFit& F, const Matrix& X, octave_idx_type i)
{
  double f = F.intercept;
  for (std::size_t j = 0; j < F.value.size (); j++)
  {
    octave_idx_type b = gamb_bin_of (F.edges[j], X(i, (octave_idx_type) j));
    if (b >= 0)
    {
      f += F.value[j](b);
    }
  }
  return f;
}

// The deviance of the current additive prediction.  For a classifier that is
// -2 times the Bernoulli log likelihood, which is the quantity MATLAB reports
// and which starts at 2 n log 2 for a balanced response; for a regression it
// is the residual sum of squares, the Gaussian deviance up to the scale.
static double
gamb_deviance (const ColumnVector& Y, const ColumnVector& f, int method)
{
  octave_idx_type n = Y.numel ();
  double dev = 0.0;

  if (method == 1)
  {
    for (octave_idx_type i = 0; i < n; i++)
    {
      double p = 1.0 / (1.0 + std::exp (-f(i)));
      // A saturated fit puts p at 0 or 1 exactly; clamp so the log is finite
      // and the deviance keeps falling instead of turning into a NaN.
      if (p < 1e-300)
      {
        p = 1e-300;
      }
      if (p > 1.0 - 1e-16)
      {
        p = 1.0 - 1e-16;
      }
      dev += Y(i) * std::log (p) + (1.0 - Y(i)) * std::log (1.0 - p);
    }
    dev *= -2.0;
  }
  else
  {
    for (octave_idx_type i = 0; i < n; i++)
    {
      double r = Y(i) - f(i);
      dev += r * r;
    }
  }

  return dev;
}

// Newton boosting of one tree per predictor per round.  METHOD 1 fits the
// logistic deviance, as a classifier is fitted, and METHOD 2 the squared
// error, as a regression is fitted; the two differ only in the seed, the
// gradient and the Hessian, so the loop is shared.
static GamBoostFit
gamb_boost (const Matrix& X, const ColumnVector& Y, int method,
            octave_idx_type maxtrees, double lrate, octave_idx_type maxsplits,
            int verbose, octave_idx_type numprint,
            const ColumnVector *F0 = nullptr)
{
  octave_idx_type n = X.rows ();
  octave_idx_type d = X.columns ();

  std::vector<BinnedPredictor> B (d);
  for (octave_idx_type j = 0; j < d; j++)
  {
    ColumnVector xj (n);
    for (octave_idx_type i = 0; i < n; i++)
    {
      xj(i) = X(i, j);
    }
    B[j] = gamb_bin (xj, GAMB_MAX_EDGES);
  }

  GamBoostFit F;
  F.edges.resize (d);
  F.value.resize (d);
  for (octave_idx_type j = 0; j < d; j++)
  {
    F.edges[j] = B[j].edges;
    F.value[j] = ColumnVector (B[j].nbins, 0.0);
  }

  // The seed.  A classifier starts at the log-odds of the response mean and a
  // regression at the response mean itself; boosting moves the first, because
  // recentring hands it the constants the shape functions give up, and leaves
  // the second where it is, because a squared-error residual already has mean
  // zero and there is nothing to hand over.
  //
  // Given F0 there is no seed to find: the caller is continuing a fit and
  // hands over the prediction it reached.  What comes back is then the
  // increment, shape values to add to the ones already held and an intercept
  // to add to the one already there, which is the contract gamb_boost_inter
  // works to as well.  A round starts at LRATE whatever its number, so the
  // trees this adds are the trees a longer run would have added.
  ColumnVector f (n);

  if (F0 != nullptr)
  {
    F.intercept = 0.0;
    f = *F0;
  }
  else
  {
    double ybar = 0.0;
    for (octave_idx_type i = 0; i < n; i++)
    {
      ybar += Y(i);
    }
    ybar /= (double) n;

    if (method == 1)
    {
      if (ybar <= 0.0 || ybar >= 1.0)
      {
        // A single-class response has an infinite log-odds and no gradient to
        // follow; the fit is the constant and every shape function stays flat.
        F.intercept = (ybar <= 0.0)
                      ? -octave::numeric_limits<double>::Inf ()
                      : octave::numeric_limits<double>::Inf ();
        F.ntrees = 0;
        F.reason = "Unable to improve the model fit.";
        F.deviance = 0.0;
        F.residuals = ColumnVector (n, 0.0);
        return F;
      }
      F.intercept = std::log (ybar / (1.0 - ybar));
    }
    else
    {
      F.intercept = ybar;
    }

    f = ColumnVector (n, F.intercept);
  }
  double dev = gamb_deviance (Y, f, method);

  // The printed trace, in the columns MATLAB prints.  RelTol here is the
  // relative improvement the round bought, which is what the patience test
  // reads.  MATLAB prints a quantity of its own under that heading which is
  // not recoverable from anything it reports, so this one is ours and the
  // documentation says so rather than implying they agree.
  if (verbose > 0)
  {
    static const char *bar =
      "|========================================================|\n";
    static const char *hdr =
      "| Type | NumTrees |  Deviance  |   RelTol   | LearnRate  |\n";
    octave_stdout << bar << hdr << bar
                  << "|    1D|         0|" << std::setw (12) << dev
                  << "|      -     |      -     |\n";
  }

  ColumnVector grad (n), hess (n);
  std::vector<double> devhist;
  devhist.push_back (dev);
  F.ntrees = 0;
  F.reason = "Terminated after training the requested number of trees.";

  for (octave_idx_type t = 0; t < maxtrees; t++)
  {
    // A round is tried at the running step and retried at half of it while it
    // fails to earn its place, which is the shape of MATLAB's own reported
    // learn rate: constant until the fit plateaus, then reduced once or more
    // on the last round it manages to accept.
    bool accepted = false;
    double step = lrate;

    for (int h = 0; h <= GAMB_MAX_HALVINGS; h++)
    {
      std::vector<ColumnVector> trial (d);
      ColumnVector fnew = f;
      double shift = 0.0;

      // The predictors are fitted in sequence and not in parallel: each tree
      // sees the residual the ones before it have already reduced.  Fitting
      // them all to the same gradient and adding them would overshoot by a
      // factor of the predictor count, and it would also drag the intercept
      // off the response mean in a regression, where MATLAB holds it there.
      for (octave_idx_type j = 0; j < d; j++)
      {
        for (octave_idx_type i = 0; i < n; i++)
        {
          if (method == 1)
          {
            double p = 1.0 / (1.0 + std::exp (-fnew(i)));
            grad(i) = Y(i) - p;
            hess(i) = p * (1.0 - p);
            if (hess(i) < 1e-12)
            {
              hess(i) = 1e-12;
            }
          }
          else
          {
            grad(i) = Y(i) - fnew(i);
            hess(i) = 1.0;
          }
        }

        trial[j] = ColumnVector (B[j].nbins, 0.0);
        gamb_fit_tree (B[j], grad, hess, maxsplits, step, trial[j]);

        // The model value the tree just added, before any bookkeeping.
        for (octave_idx_type i = 0; i < n; i++)
        {
          if (B[j].bin(i) >= 0)
          {
            fnew(i) += trial[j](B[j].bin(i));
          }
        }

        // Recentre over the observations, not over the bins: a shape function
        // is held to mean zero on the data it was fitted to, and the constant
        // it gives up goes to the intercept.  This moves value between the two
        // and leaves the prediction alone, which is why a regression intercept
        // stays at the response mean while a classifier's walks away from its
        // seed: a squared-error residual keeps mean zero from round to round
        // and hands over nothing, a logistic one does not.
        double m = 0.0;
        octave_idx_type cnt = 0;
        for (octave_idx_type i = 0; i < n; i++)
        {
          if (B[j].bin(i) >= 0)
          {
            m += trial[j](B[j].bin(i));
            cnt++;
          }
        }
        if (cnt > 0)
        {
          m /= (double) cnt;
          for (octave_idx_type b = 0; b < B[j].nbins; b++)
          {
            trial[j](b) -= m;
          }
          shift += m;
        }
      }

      double devnew = gamb_deviance (Y, fnew, method);
      double devold = dev;

      // A round only has to improve the deviance at all to be kept.  Whether
      // it improved enough to be worth continuing is the patience test below,
      // which looks across a window of rounds rather than at this one.
      if (devnew < dev)
      {
        for (octave_idx_type j = 0; j < d; j++)
        {
          for (octave_idx_type b = 0; b < B[j].nbins; b++)
          {
            F.value[j](b) += trial[j](b);
          }
        }
        F.intercept += shift;
        f = fnew;
        dev = devnew;
        F.ntrees++;
        accepted = true;
        if (verbose > 0 && (F.ntrees == 1 || F.ntrees % numprint == 0))
        {
          double rel = (devold > 0.0) ? (devold - devnew) / devold : 0.0;
          octave_stdout << "|    1D|" << std::setw (10) << F.ntrees << "|"
                        << std::setw (12) << devnew << "|"
                        << std::setw (12) << rel << "|"
                        << std::setw (12) << step << "|\n";
        }
        break;
      }

      step *= 0.5;
    }

    if (! accepted)
    {
      F.reason = "Unable to improve the model fit.";
      break;
    }

    // Patience: stop once a whole window of rounds has failed to better the
    // deviance the window opened at by more than the tolerance.  A single
    // round that buys little is not evidence of convergence, since the gain
    // from round to round is noisy; a window of them is.
    devhist.push_back (dev);
    octave_idx_type h = (octave_idx_type) devhist.size ();
    if (h > GAMB_PATIENCE)
    {
      double ref = devhist[(std::size_t) (h - 1 - GAMB_PATIENCE)];
      double need = ref - GAMB_REL_TOL * std::fabs (ref);
      bool better = false;
      for (octave_idx_type k = h - GAMB_PATIENCE; k < h; k++)
      {
        if (devhist[(std::size_t) k] < need)
        {
          better = true;
          break;
        }
      }
      if (! better)
      {
        F.reason = "Unable to improve the model fit.";
        break;
      }
    }
  }

  F.deviance = dev;
  F.residuals = ColumnVector (n);
  for (octave_idx_type i = 0; i < n; i++)
  {
    if (method == 1)
    {
      F.residuals(i) = Y(i) - 1.0 / (1.0 + std::exp (-f(i)));
    }
    else
    {
      F.residuals(i) = Y(i) - f(i);
    }
  }

  return F;
}

// What the interaction search reports for one candidate pair.  The p-value is
// not computed here: turning F into a probability needs the incomplete beta,
// and the package already ships a MATLAB-verified fcdf, so the caller does it
// there rather than a second implementation being carried in C++.  Selecting
// which pairs to keep is policy and belongs beside it.
struct GamPairStat
{
  octave_idx_type j;
  octave_idx_type k;
  double F;
  double df1;
  double df2;
};

// Test one pair for an interaction the additive terms have not already
// explained.  The residuals of the predictor phase are laid on the coarse
// grid of the two predictors and decomposed as a two-way layout: what the rows
// explain, what the columns explain, and what only the cells explain.  The
// last is the interaction, and the F ratio against the within-cell error is
// what ranks the pair.
//
// The decomposition is sequential (Type I): with unequal cell counts the row
// and column sums of squares are not orthogonal, so this is an approximation
// rather than an identity, and a small negative interaction term is clamped to
// zero.  That is sound for ranking candidates, which is all it is used for.
static GamPairStat
gamb_pair_stat (const BinnedPredictor& Bj, const BinnedPredictor& Bk,
                const ColumnVector& r, octave_idx_type j, octave_idx_type k)
{
  octave_idx_type n = r.numel ();
  octave_idx_type rj = Bj.nbins;
  octave_idx_type rk = Bk.nbins;

  std::vector<double> csum ((std::size_t) (rj * rk), 0.0);
  std::vector<double> ccnt ((std::size_t) (rj * rk), 0.0);
  std::vector<double> rsum ((std::size_t) rj, 0.0);
  std::vector<double> rcnt ((std::size_t) rj, 0.0);
  std::vector<double> ksum ((std::size_t) rk, 0.0);
  std::vector<double> kcnt ((std::size_t) rk, 0.0);

  double tot = 0.0;
  double cnt = 0.0;

  for (octave_idx_type i = 0; i < n; i++)
  {
    octave_idx_type a = Bj.bin(i);
    octave_idx_type b = Bk.bin(i);
    if (a < 0 || b < 0)
    {
      continue;
    }
    std::size_t c = (std::size_t) (a * rk + b);
    csum[c] += r(i);
    ccnt[c] += 1.0;
    rsum[(std::size_t) a] += r(i);
    rcnt[(std::size_t) a] += 1.0;
    ksum[(std::size_t) b] += r(i);
    kcnt[(std::size_t) b] += 1.0;
    tot += r(i);
    cnt += 1.0;
  }

  GamPairStat S;
  S.j = j;
  S.k = k;
  S.F = 0.0;
  S.df1 = 0.0;
  S.df2 = 0.0;

  if (cnt < 4.0)
  {
    return S;
  }

  double gm = tot / cnt;

  double sscell = 0.0;
  octave_idx_type ncell = 0;
  for (std::size_t c = 0; c < csum.size (); c++)
  {
    if (ccnt[c] > 0.0)
    {
      double m = csum[c] / ccnt[c];
      sscell += ccnt[c] * (m - gm) * (m - gm);
      ncell++;
    }
  }

  double ssrow = 0.0;
  octave_idx_type nrow = 0;
  for (std::size_t a = 0; a < rsum.size (); a++)
  {
    if (rcnt[a] > 0.0)
    {
      double m = rsum[a] / rcnt[a];
      ssrow += rcnt[a] * (m - gm) * (m - gm);
      nrow++;
    }
  }

  double sscol = 0.0;
  octave_idx_type ncol = 0;
  for (std::size_t b = 0; b < ksum.size (); b++)
  {
    if (kcnt[b] > 0.0)
    {
      double m = ksum[b] / kcnt[b];
      sscol += kcnt[b] * (m - gm) * (m - gm);
      ncol++;
    }
  }

  double ssint = sscell - ssrow - sscol;
  if (ssint < 0.0)
  {
    ssint = 0.0;
  }

  double sserr = 0.0;
  for (octave_idx_type i = 0; i < n; i++)
  {
    octave_idx_type a = Bj.bin(i);
    octave_idx_type b = Bk.bin(i);
    if (a < 0 || b < 0)
    {
      continue;
    }
    std::size_t c = (std::size_t) (a * rk + b);
    double m = csum[c] / ccnt[c];
    sserr += (r(i) - m) * (r(i) - m);
  }

  double df1 = (double) ((nrow - 1) * (ncol - 1));
  double df2 = cnt - (double) ncell;

  if (df1 <= 0.0 || df2 <= 0.0 || sserr <= 0.0)
  {
    return S;
  }

  S.F = (ssint / df1) / (sserr / df2);
  S.df1 = df1;
  S.df2 = df2;

  return S;
}

// One pair's fitted interaction surface: a value per cell of the two coarse
// grids, accumulated over every tree that was boosted onto it.  The surface is
// held on the DETECTION grid rather than the fitting grid of the main effects.
// MATLAB reports exactly two grids, so there is no third for interactions to
// live on; a tree limited to MaxNumSplitsPerInteraction splits carves at most
// that many regions plus one, so resolution past a handful of bins per axis
// buys nothing; and the Explainable Boosting Machine bins its interactions
// coarser than its main effects for the same reason.  It also keeps the
// surface at 64 doubles instead of the 520 KB a 255 by 255 accumulator would
// need for every pair.
struct GamInterTerm
{
  octave_idx_type j;
  octave_idx_type k;
  Matrix value;                     // nbins(j) x nbins(k)
};

// A rectangle of the two-dimensional bin grid, as the tree grows it.
struct GridRegion
{
  octave_idx_type r0, r1, c0, c1;
  double gain;
  int dim;                          // 0 splits rows, 1 splits columns
  octave_idx_type cut;              // last row or column of the near side
};

// Totals over a rectangle, from an integral image with one row and column of
// leading zeros.
static inline double
gamb_rect (const Matrix& I, octave_idx_type r0, octave_idx_type r1,
           octave_idx_type c0, octave_idx_type c1)
{
  return I(r1 + 1, c1 + 1) - I(r0, c1 + 1) - I(r1 + 1, c0) + I(r0, c0);
}

// The best cut of a rectangle, over both directions.  Returns the gain and
// sets the direction and the position; a rectangle that cannot be cut usefully
// returns -1.
static double
gamb_best_cut2 (const Matrix& IG, const Matrix& IH, GridRegion& R)
{
  double best = -1.0;
  R.dim = 0;
  R.cut = -1;

  double Gt = gamb_rect (IG, R.r0, R.r1, R.c0, R.c1);
  double Ht = gamb_rect (IH, R.r0, R.r1, R.c0, R.c1);

  for (octave_idx_type a = R.r0; a < R.r1; a++)
  {
    double GL = gamb_rect (IG, R.r0, a, R.c0, R.c1);
    double HL = gamb_rect (IH, R.r0, a, R.c0, R.c1);
    double g = gamb_gain (GL, HL, Gt - GL, Ht - HL);
    if (g > best)
    {
      best = g;
      R.dim = 0;
      R.cut = a;
    }
  }

  for (octave_idx_type b = R.c0; b < R.c1; b++)
  {
    double GL = gamb_rect (IG, R.r0, R.r1, R.c0, b);
    double HL = gamb_rect (IH, R.r0, R.r1, R.c0, b);
    double g = gamb_gain (GL, HL, Gt - GL, Ht - HL);
    if (g > best)
    {
      best = g;
      R.dim = 1;
      R.cut = b;
    }
  }

  return best;
}

// Fit one tree over a pair's grid and accumulate its leaf values, scaled by
// the step, into VAL.  Grown best-first like the univariate trees, and split
// in whichever direction buys more, which is what lets a pair term represent
// something neither predictor could alone.
static void
gamb_fit_tree2 (const BinnedPredictor& Bj, const BinnedPredictor& Bk,
                const ColumnVector& grad, const ColumnVector& hess,
                octave_idx_type maxsplits, double step, Matrix& val)
{
  octave_idx_type n = grad.numel ();
  octave_idx_type rj = Bj.nbins;
  octave_idx_type rk = Bk.nbins;

  Matrix IG (rj + 1, rk + 1, 0.0);
  Matrix IH (rj + 1, rk + 1, 0.0);

  for (octave_idx_type i = 0; i < n; i++)
  {
    octave_idx_type a = Bj.bin(i);
    octave_idx_type b = Bk.bin(i);
    if (a < 0 || b < 0)
    {
      continue;
    }
    IG(a + 1, b + 1) += grad(i);
    IH(a + 1, b + 1) += hess(i);
  }

  for (octave_idx_type a = 1; a <= rj; a++)
  {
    for (octave_idx_type b = 1; b <= rk; b++)
    {
      IG(a, b) += IG(a - 1, b) + IG(a, b - 1) - IG(a - 1, b - 1);
      IH(a, b) += IH(a - 1, b) + IH(a, b - 1) - IH(a - 1, b - 1);
    }
  }

  std::vector<GridRegion> leaves;
  GridRegion root;
  root.r0 = 0;
  root.r1 = rj - 1;
  root.c0 = 0;
  root.c1 = rk - 1;
  root.gain = gamb_best_cut2 (IG, IH, root);
  leaves.push_back (root);

  for (octave_idx_type s = 0; s < maxsplits; s++)
  {
    std::size_t pick = 0;
    double best = -1.0;
    for (std::size_t q = 0; q < leaves.size (); q++)
    {
      if (leaves[q].gain > best)
      {
        best = leaves[q].gain;
        pick = q;
      }
    }
    if (best <= 0.0)
    {
      break;
    }

    GridRegion L = leaves[pick];
    GridRegion Rr = leaves[pick];
    if (leaves[pick].dim == 0)
    {
      L.r1 = leaves[pick].cut;
      Rr.r0 = leaves[pick].cut + 1;
    }
    else
    {
      L.c1 = leaves[pick].cut;
      Rr.c0 = leaves[pick].cut + 1;
    }
    L.gain = gamb_best_cut2 (IG, IH, L);
    Rr.gain = gamb_best_cut2 (IG, IH, Rr);

    leaves[pick] = L;
    leaves.push_back (Rr);
  }

  for (std::size_t q = 0; q < leaves.size (); q++)
  {
    double Gq = gamb_rect (IG, leaves[q].r0, leaves[q].r1,
                           leaves[q].c0, leaves[q].c1);
    double Hq = gamb_rect (IH, leaves[q].r0, leaves[q].r1,
                           leaves[q].c0, leaves[q].c1);
    if (Hq <= 1e-12)
    {
      continue;
    }
    double leaf = step * Gq / Hq;
    for (octave_idx_type a = leaves[q].r0; a <= leaves[q].r1; a++)
    {
      for (octave_idx_type b = leaves[q].c0; b <= leaves[q].c1; b++)
      {
        val(a, b) += leaf;
      }
    }
  }
}

// What the interaction phase produces.
struct GamInterFit
{
  std::vector<RowVector> edges;     // the coarse grid, one per predictor
  std::vector<GamInterTerm> term;   // one per selected pair
  double shift;                     // what recentring handed the intercept
  octave_idx_type ntrees;
  std::string reason;
  double deviance;
  ColumnVector residuals;
};

// Boost trees over the selected pairs, continuing from the additive prediction
// the predictor phase left in F0 rather than refitting it.  MATLAB's own trace
// opens its interaction block at the deviance the predictor block ended on, so
// the two phases share a running fit and differ only in what they are allowed
// to split.  The pairs are fitted in sequence within a round, for the same
// reason the predictors are.
static GamInterFit
gamb_boost_inter (const Matrix& X, const ColumnVector& Y,
                  const ColumnVector& F0, int method, const Matrix& pairs,
                  octave_idx_type maxtrees, double lrate,
                  octave_idx_type maxsplits)
{
  octave_idx_type n = X.rows ();
  octave_idx_type d = X.columns ();
  octave_idx_type np = pairs.rows ();

  std::vector<BinnedPredictor> B ((std::size_t) d);
  for (octave_idx_type j = 0; j < d; j++)
  {
    ColumnVector xj (n);
    for (octave_idx_type i = 0; i < n; i++)
    {
      xj(i) = X(i, j);
    }
    B[(std::size_t) j] = gamb_bin (xj, GAMB_PAIR_EDGES);
  }

  GamInterFit F;
  F.edges.resize ((std::size_t) d);
  for (octave_idx_type j = 0; j < d; j++)
  {
    F.edges[(std::size_t) j] = B[(std::size_t) j].edges;
  }

  F.term.resize ((std::size_t) np);
  for (octave_idx_type q = 0; q < np; q++)
  {
    octave_idx_type j = (octave_idx_type) pairs(q, 0) - 1;
    octave_idx_type k = (octave_idx_type) pairs(q, 1) - 1;
    F.term[(std::size_t) q].j = j;
    F.term[(std::size_t) q].k = k;
    F.term[(std::size_t) q].value = Matrix (B[(std::size_t) j].nbins,
                                            B[(std::size_t) k].nbins, 0.0);
  }

  F.shift = 0.0;
  F.ntrees = 0;
  F.reason = "Terminated after training the requested number of trees.";

  ColumnVector f = F0;
  double dev = gamb_deviance (Y, f, method);
  std::vector<double> devhist;
  devhist.push_back (dev);

  ColumnVector grad (n), hess (n);

  for (octave_idx_type t = 0; t < maxtrees; t++)
  {
    bool accepted = false;
    double step = lrate;

    for (int h = 0; h <= GAMB_MAX_HALVINGS; h++)
    {
      std::vector<Matrix> trial ((std::size_t) np);
      ColumnVector fnew = f;
      double shift = 0.0;

      for (octave_idx_type q = 0; q < np; q++)
      {
        octave_idx_type j = F.term[(std::size_t) q].j;
        octave_idx_type k = F.term[(std::size_t) q].k;
        const BinnedPredictor& Bj = B[(std::size_t) j];
        const BinnedPredictor& Bk = B[(std::size_t) k];

        for (octave_idx_type i = 0; i < n; i++)
        {
          if (method == 1)
          {
            double p = 1.0 / (1.0 + std::exp (-fnew(i)));
            grad(i) = Y(i) - p;
            hess(i) = p * (1.0 - p);
            if (hess(i) < 1e-12)
            {
              hess(i) = 1e-12;
            }
          }
          else
          {
            grad(i) = Y(i) - fnew(i);
            hess(i) = 1.0;
          }
        }

        trial[(std::size_t) q] = Matrix (Bj.nbins, Bk.nbins, 0.0);
        gamb_fit_tree2 (Bj, Bk, grad, hess, maxsplits, step,
                        trial[(std::size_t) q]);

        for (octave_idx_type i = 0; i < n; i++)
        {
          octave_idx_type a = Bj.bin(i);
          octave_idx_type b = Bk.bin(i);
          if (a >= 0 && b >= 0)
          {
            fnew(i) += trial[(std::size_t) q](a, b);
          }
        }

        // Recentred like a shape function, and for the same reason: a pair
        // term carries only what the two predictors do together, and whatever
        // constant it picked up belongs to the intercept.
        double m = 0.0;
        octave_idx_type cnt = 0;
        for (octave_idx_type i = 0; i < n; i++)
        {
          octave_idx_type a = Bj.bin(i);
          octave_idx_type b = Bk.bin(i);
          if (a >= 0 && b >= 0)
          {
            m += trial[(std::size_t) q](a, b);
            cnt++;
          }
        }
        if (cnt > 0)
        {
          m /= (double) cnt;
          for (octave_idx_type a = 0; a < Bj.nbins; a++)
          {
            for (octave_idx_type b = 0; b < Bk.nbins; b++)
            {
              trial[(std::size_t) q](a, b) -= m;
            }
          }
          shift += m;
        }
      }

      double devnew = gamb_deviance (Y, fnew, method);

      if (devnew < dev)
      {
        for (octave_idx_type q = 0; q < np; q++)
        {
          octave_idx_type rj = F.term[(std::size_t) q].value.rows ();
          octave_idx_type rk = F.term[(std::size_t) q].value.columns ();
          for (octave_idx_type a = 0; a < rj; a++)
          {
            for (octave_idx_type b = 0; b < rk; b++)
            {
              F.term[(std::size_t) q].value(a, b)
                += trial[(std::size_t) q](a, b);
            }
          }
        }
        F.shift += shift;
        f = fnew;
        dev = devnew;
        F.ntrees++;
        accepted = true;
        break;
      }

      step *= 0.5;
    }

    if (! accepted)
    {
      F.reason = "Unable to improve the model fit.";
      break;
    }

    devhist.push_back (dev);
    octave_idx_type hn = (octave_idx_type) devhist.size ();
    if (hn > GAMB_PATIENCE)
    {
      double ref = devhist[(std::size_t) (hn - 1 - GAMB_PATIENCE)];
      double need = ref - GAMB_REL_TOL * std::fabs (ref);
      bool better = false;
      for (octave_idx_type q = hn - GAMB_PATIENCE; q < hn; q++)
      {
        if (devhist[(std::size_t) q] < need)
        {
          better = true;
          break;
        }
      }
      if (! better)
      {
        F.reason = "Unable to improve the model fit.";
        break;
      }
    }
  }

  F.deviance = dev;
  F.residuals = ColumnVector (n);
  for (octave_idx_type i = 0; i < n; i++)
  {
    if (method == 1)
    {
      F.residuals(i) = Y(i) - 1.0 / (1.0 + std::exp (-f(i)));
    }
    else
    {
      F.residuals(i) = Y(i) - f(i);
    }
  }

  return F;
}
