/*
Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This file is part of the statistics package for GNU Octave.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <numeric>
#include <random>
#include <string>
#include <vector>

// Growing a binary decision tree by recursive partitioning.
//
// The cost of a tree is the split search, and the search is a sort.  Sorting
// each node's rows afresh for every predictor is what a direct reading of the
// algorithm gives, and it repeats the same comparisons at every level: the
// order of a node's rows along a predictor is the order its parent already
// had, with the rows that went the other way removed.
//
// So the rows are sorted once per predictor at the root, and every split
// partitions those orders in place, stably, into the two children.  A node
// then owns one contiguous range of each predictor's index array, already in
// order, and a level costs a pass over the data rather than a sort of it.
//
// Within a node the search is a single sweep per predictor.  Class weights
// accumulate from the left, so both sides of every candidate cut are known
// from one running total, and no cut is evaluated twice.
//
// Ties are broken exactly as a sweep in this order finds them: the first
// predictor, and within it the lowest cut point, that strictly beats what is
// held.  The arithmetic below is written in the order the reference
// implementation performs it, because a tie decided by a last-bit difference
// would move a split.

// Two splits can be equally good and be computed as differing in the last
// bits, because each predictor accumulates the same weights in its own sort
// order.  Left to decide the winner, that noise makes the tree depend on the
// compiler rather than on the data.  Splits within this much of each other,
// relatively, are therefore held to be equal and the earlier one is kept.
//
// Measured on eight thousand rows against MATLAB R2024a: every one of its
// 1396 splits is within 1.2e-14 of our best, so the noise is of that order,
// while the tolerance may be raised to 1e-6 before the grown tree changes at
// all.  1e-12 sits a hundredfold above the noise and far below any gain
// difference that carries meaning.
static const double GAIN_TIE_TOL = 1e-12;

// One node of the tree under construction.
struct Node
{
  octave_idx_type left, right, parent;
  octave_idx_type cutvar;              // 0 when the node is a leaf
  double cutval;
  octave_idx_type start, stop;         // its range in every sorted order
  octave_idx_type nsize;
  double nweight;
  std::vector<double> catleft;         // the levels sent left, sorted, when
  std::vector<double> catright;        // the cut is on a categorical predictor
};

// The tree grows the same way whichever it is fitting.  Only what a node
// accumulates, how impure that makes it and what a leaf then answers with
// differ, so those three are the only things that branch on the type.
enum Criterion { GDI, DEVIANCE, MSE };

// A regression node carries the weight, the weighted sum of the response and
// its weighted sum of squares, which give the mean and the sum of squared
// error about it in closed form.
struct Moments
{
  double w, wy, wyy;

  void clear (void) { w = wy = wyy = 0.0; }

  void add (double wi, double yi)
  {
    w += wi;
    wy += wi * yi;
    wyy += wi * yi * yi;
  }

  double mean (void) const { return (w > 0.0) ? wy / w : 0.0; }

  // Sum of squared error about the node mean.  Never let rounding take it
  // below zero, which a node holding one distinct value can do.
  double sse (void) const
  {
    if (w <= 0.0)
      return 0.0;
    double v = wyy - wy * wy / w;
    return (v > 0.0) ? v : 0.0;
  }
};

static inline double
impurity (const std::vector<double>& cw, double total, octave_idx_type K,
          Criterion crit)
{
  double imp = 0.0;
  if (crit == GDI)
    {
      for (octave_idx_type k = 0; k < K; k++)
        {
          double pr = cw[k] / total;
          imp += pr * pr;
        }
      return 1.0 - imp;
    }

  for (octave_idx_type k = 0; k < K; k++)
    {
      double pr = cw[k] / total;
      if (pr > 0.0)
        imp -= pr * std::log (pr);
    }
  return imp;
}
// How a node with more than two classes splits a categorical predictor: the
// exact search, one of the three heuristics MATLAB names, or its automatic
// choice between them.
enum CatAlg { CAT_AUTO, CAT_EXACT, CAT_PULLLEFT, CAT_PCA, CAT_OVA };

// Everything the caller settles before the engine starts.  A classifier and a
// regression differ here in three fields and nowhere else.
struct TreeOpts
{
  octave_idx_type K, minparent, minleaf, maxsplits;
  bool mergeleaves;
  double qetoler;
  Criterion crit;
  octave_idx_type nvars;               // predictors tried per node, 0 for all
  std::uint32_t seed;                  // seeds the draw when nvars is set
  std::vector<bool> iscat;             // per predictor, empty when none is
  double maxcat;                       // levels an automatic exact search takes
  CatAlg catalg;
};

// One value in 0 .. range-1, every value equally likely.  The raw output of
// std::mt19937 is fixed by the standard where the distributions built on it
// are not, so drawing from it directly keeps a seeded tree the same under
// every C++ library, and rejecting the low values a bare modulus would favour
// keeps the draw unbiased.
static octave_idx_type
draw_below (std::mt19937& rng, octave_idx_type range)
{
  const std::uint32_t r = static_cast<std::uint32_t> (range);
  const std::uint32_t reject = (0u - r) % r;
  std::uint32_t x;
  do
    x = static_cast<std::uint32_t> (rng ());
  while (x < reject);
  return static_cast<octave_idx_type> (x % r);
}

// The eigenvector of the largest eigenvalue of the symmetric K by K matrix A,
// by cyclic Jacobi rotations, signed so that its largest entry in magnitude
// is positive.
static std::vector<double>
leading_eigenvector (std::vector<double> A, octave_idx_type K)
{
  std::vector<double> V (K * K, 0.0);
  for (octave_idx_type i = 0; i < K; i++)
    V[i * K + i] = 1.0;
  for (int sweep = 0; sweep < 100; sweep++)
    {
      double off = 0.0;
      for (octave_idx_type p = 0; p < K; p++)
        for (octave_idx_type q = p + 1; q < K; q++)
          off += A[p * K + q] * A[p * K + q];
      if (off < 1e-30)
        break;
      for (octave_idx_type p = 0; p < K; p++)
        for (octave_idx_type q = p + 1; q < K; q++)
          {
            const double apq = A[p * K + q];
            if (std::fabs (apq) < 1e-300)
              continue;
            const double theta = (A[q * K + q] - A[p * K + p]) / (2.0 * apq);
            const double t = ((theta >= 0.0) ? 1.0 : -1.0)
                             / (std::fabs (theta)
                                + std::sqrt (theta * theta + 1.0));
            const double c = 1.0 / std::sqrt (t * t + 1.0);
            const double sn = t * c;
            for (octave_idx_type k = 0; k < K; k++)
              {
                const double akp = A[k * K + p];
                const double akq = A[k * K + q];
                A[k * K + p] = c * akp - sn * akq;
                A[k * K + q] = sn * akp + c * akq;
              }
            for (octave_idx_type k = 0; k < K; k++)
              {
                const double apk = A[p * K + k];
                const double aqk = A[q * K + k];
                A[p * K + k] = c * apk - sn * aqk;
                A[q * K + k] = sn * apk + c * aqk;
              }
            for (octave_idx_type k = 0; k < K; k++)
              {
                const double vkp = V[k * K + p];
                const double vkq = V[k * K + q];
                V[k * K + p] = c * vkp - sn * vkq;
                V[k * K + q] = sn * vkp + c * vkq;
              }
          }
    }
  octave_idx_type top = 0;
  for (octave_idx_type i = 1; i < K; i++)
    if (A[i * K + i] > A[top * K + top])
      top = i;
  std::vector<double> v (K);
  octave_idx_type big = 0;
  for (octave_idx_type k = 0; k < K; k++)
    {
      v[k] = V[k * K + top];
      if (std::fabs (v[k]) > std::fabs (v[big]))
        big = k;
    }
  if (v[big] < 0.0)
    for (octave_idx_type k = 0; k < K; k++)
      v[k] = -v[k];
  return v;
}

// The best split of a node's rows on a categorical predictor, as its gain per
// unit of the node's weight, with the levels sent each way.  The rows
// [start, have) of the predictor's order hold its known values, sorted, so
// the rows of one level are contiguous.
//
// A regression orders the levels by their mean response and two classes by
// the probability of the first class, and either takes the best of the L - 1
// cuts of that order, which is the exact search.  More classes search every
// partition when the node holds at most MaxNumCategories levels, the lowest
// level always on the left and ties kept by the first partition in the order
// the other levels' membership of the right counts up in binary.  Above that
// the automatic choice keeps the best of OVAbyClass, PCA and PullLeft for up
// to four classes and of PCA and PullLeft for more, the earlier on a tie, as
// measured on MATLAB R2024a; the heuristics follow MATLAB's descriptions of
// them.  A regression keeps the lower means on the left; every classifier
// split from an order puts the lowest level on the left.
static double
categorical_split (const double *xj, const octave_idx_type *oj,
                   octave_idx_type start, octave_idx_type have,
                   const std::vector<octave_idx_type>& y,
                   const std::vector<double>& resp,
                   const std::vector<double>& w, octave_idx_type K,
                   Criterion crit, bool isreg,
                   const std::vector<double>& present,
                   const Moments& presmom, double wp, double presentimp,
                   double total, const TreeOpts& o,
                   std::vector<double>& catleft,
                   std::vector<double>& catright)
{
  catleft.clear ();
  catright.clear ();

  // Each level's value, row count, weight and class weights or moments
  std::vector<double> lev, lw, lcw;
  std::vector<octave_idx_type> cnt;
  std::vector<Moments> lmom;
  for (octave_idx_type i = start; i < have; i++)
    {
      const octave_idx_type r = oj[i];
      if (lev.empty () || xj[r] != lev.back ())
        {
          lev.push_back (xj[r]);
          cnt.push_back (0);
          lw.push_back (0.0);
          if (isreg)
            {
              Moments m;
              m.clear ();
              lmom.push_back (m);
            }
          else
            lcw.insert (lcw.end (), K, 0.0);
        }
      const std::size_t l = lev.size () - 1;
      cnt[l]++;
      lw[l] += w[r];
      if (isreg)
        lmom[l].add (w[r], resp[r]);
      else
        lcw[l * K + y[r]] += w[r];
    }
  const octave_idx_type L = static_cast<octave_idx_type> (lev.size ());
  if (L < 2)
    return 0.0;
  const octave_idx_type ntot = have - start;

  std::vector<double> running (K), other (K);

  // The gain of sending the levels marked in inL left, and whether both
  // sides keep MinLeaf rows and some weight.
  auto gain_of = [&] (const std::vector<char>& inL, bool& ok) -> double
    {
      octave_idx_type nl = 0;
      double wl = 0.0;
      Moments ml;
      ml.clear ();
      std::fill (running.begin (), running.end (), 0.0);
      for (octave_idx_type l = 0; l < L; l++)
        if (inL[l])
          {
            nl += cnt[l];
            wl += lw[l];
            if (isreg)
              {
                ml.w += lmom[l].w;
                ml.wy += lmom[l].wy;
                ml.wyy += lmom[l].wyy;
              }
            else
              for (octave_idx_type k = 0; k < K; k++)
                running[k] += lcw[l * K + k];
          }
      const double wr = wp - wl;
      ok = (nl >= o.minleaf && ntot - nl >= o.minleaf);
      if (! (wl > 0.0 && wr > 0.0))
        {
          ok = false;
          return -std::numeric_limits<double>::infinity ();
        }
      double impl, impr;
      if (isreg)
        {
          Moments rest;
          rest.w = presmom.w - ml.w;
          rest.wy = presmom.wy - ml.wy;
          rest.wyy = presmom.wyy - ml.wyy;
          impl = ml.sse () / wl;
          impr = rest.sse () / wr;
        }
      else
        {
          for (octave_idx_type k = 0; k < K; k++)
            other[k] = present[k] - running[k];
          impl = impurity (running, wl, K, crit);
          impr = impurity (other, wr, K, crit);
        }
      return (wp * presentimp - wl * impl - wr * impr) / total;
    };

  double best = 0.0;
  bool found = false;
  std::vector<char> bestL;
  auto consider = [&] (const std::vector<char>& inL)
    {
      bool ok;
      const double g = gain_of (inL, ok);
      if (ok && (! found || g > best + std::fabs (best) * GAIN_TIE_TOL))
        {
          best = g;
          bestL = inL;
          found = true;
        }
    };

  // The L - 1 cuts of the levels in ascending order of key, ties kept in
  // the order of the levels themselves
  auto scan_order = [&] (const std::vector<double>& key)
    {
      std::vector<octave_idx_type> ord (L);
      std::iota (ord.begin (), ord.end (), 0);
      std::stable_sort (ord.begin (), ord.end (),
                        [&key] (octave_idx_type a, octave_idx_type b)
                        { return key[a] < key[b]; });
      std::vector<char> inL (L, 0);
      for (octave_idx_type m = 0; m < L - 1; m++)
        {
          inL[ord[m]] = 1;
          consider (inL);
        }
    };

  // The probability of class c at level l
  auto prob = [&] (octave_idx_type l, octave_idx_type c) -> double
    {
      return (lw[l] > 0.0) ? lcw[l * K + c] / lw[l] : 0.0;
    };

  bool lowleft = ! isreg;
  if (isreg)
    {
      std::vector<double> key (L);
      for (octave_idx_type l = 0; l < L; l++)
        key[l] = lmom[l].mean ();
      scan_order (key);
    }
  else if (K <= 2)
    {
      std::vector<double> key (L);
      for (octave_idx_type l = 0; l < L; l++)
        key[l] = prob (l, 0);
      scan_order (key);
    }
  else if (o.catalg == CAT_EXACT
           || (o.catalg == CAT_AUTO && static_cast<double> (L) <= o.maxcat))
    {
      if (L > 31)
        error ("treetrain: an exact categorical search cannot take more "
               "than 31 levels.");
      std::vector<char> inL (L, 1);
      const std::uint64_t stop = std::uint64_t (1) << (L - 1);
      for (std::uint64_t mask = 1; mask < stop; mask++)
        {
          for (octave_idx_type l = 1; l < L; l++)
            inL[l] = ! ((mask >> (l - 1)) & 1);
          consider (inL);
        }
    }
  else
    {
      lowleft = false;
      const bool automatic = (o.catalg == CAT_AUTO);
      if (o.catalg == CAT_OVA || (automatic && K <= 4))
        for (octave_idx_type c = 0; c < K; c++)
          {
            std::vector<double> key (L);
            for (octave_idx_type l = 0; l < L; l++)
              key[l] = -prob (l, c);
            scan_order (key);
          }
      if (o.catalg == CAT_PCA || automatic)
        {
          double ws = 0.0;
          std::vector<double> mean (K, 0.0), C (K * K, 0.0);
          for (octave_idx_type l = 0; l < L; l++)
            {
              ws += lw[l];
              for (octave_idx_type k = 0; k < K; k++)
                mean[k] += lw[l] * prob (l, k);
            }
          for (octave_idx_type k = 0; k < K; k++)
            mean[k] /= ws;
          for (octave_idx_type l = 0; l < L; l++)
            for (octave_idx_type a = 0; a < K; a++)
              for (octave_idx_type b = 0; b < K; b++)
                C[a * K + b] += lw[l] * (prob (l, a) - mean[a])
                                * (prob (l, b) - mean[b]) / ws;
          const std::vector<double> v = leading_eigenvector (C, K);
          std::vector<double> key (L, 0.0);
          for (octave_idx_type l = 0; l < L; l++)
            for (octave_idx_type k = 0; k < K; k++)
              key[l] += prob (l, k) * v[k];
          scan_order (key);
        }
      if (o.catalg == CAT_PULLLEFT || automatic)
        {
          std::vector<char> inL (L, 0);
          octave_idx_type nright = L;
          while (nright > 1)
            {
              std::vector<octave_idx_type> cand;
              for (octave_idx_type c = 0; c < K; c++)
                {
                  octave_idx_type top = -1;
                  for (octave_idx_type l = 0; l < L; l++)
                    if (! inL[l] && (top < 0 || prob (l, c) > prob (top, c)))
                      top = l;
                  if (std::find (cand.begin (), cand.end (), top)
                      == cand.end ())
                    cand.push_back (top);
                }
              octave_idx_type pick = -1;
              double pg = 0.0;
              for (octave_idx_type x : cand)
                {
                  inL[x] = 1;
                  bool ok;
                  const double g = gain_of (inL, ok);
                  inL[x] = 0;
                  if (pick < 0 || g > pg)
                    {
                      pick = x;
                      pg = g;
                    }
                }
              inL[pick] = 1;
              nright--;
              consider (inL);
            }
        }
    }

  if (! found)
    return 0.0;
  if (lowleft && ! bestL[0])
    for (octave_idx_type l = 0; l < L; l++)
      bestL[l] = ! bestL[l];
  for (octave_idx_type l = 0; l < L; l++)
    (bestL[l] ? catleft : catright).push_back (lev[l]);
  return best;
}

// Grow a tree and return it as one node table.  Shared by treetrain, which
// fits it, and used by treepredict through the descent below, so that the two
// halves of the same rule cannot drift apart.
static octave_scalar_map
tree_build (const Matrix& X, const ColumnVector& yv, const ColumnVector& wv,
            const TreeOpts& o)
{
  const octave_idx_type n = X.rows ();
  const octave_idx_type p = X.columns ();

  const octave_idx_type K = o.K;
  const octave_idx_type minparent = o.minparent;
  const octave_idx_type minleaf = o.minleaf;
  const octave_idx_type maxsplits = o.maxsplits;
  const bool mergeleaves = o.mergeleaves;
  const Criterion crit = o.crit;
  const double qetoler = o.qetoler;

  const bool isreg = (crit == MSE);

  // Sampling predictors per node grows the trees of a random forest.  With
  // every predictor tried nothing is drawn, so such a tree is exactly the one
  // grown without sampling.
  const octave_idx_type nvars = (o.nvars > 0 && o.nvars < p) ? o.nvars : p;
  const bool sampling = (nvars < p);
  std::mt19937 rng (o.seed);
  std::vector<octave_idx_type> cand (p);
  std::iota (cand.begin (), cand.end (), 0);

  if (! isreg && K < 1)
    error ("treetrain: NumClasses must be a positive integer.");

  // Class index of each row for a classifier, the response itself for a
  // regression, and the weight either way.
  std::vector<octave_idx_type> y (isreg ? 0 : n);
  std::vector<double> resp (isreg ? n : 0);
  std::vector<double> w (n);
  for (octave_idx_type i = 0; i < n; i++)
    {
      if (isreg)
        resp[i] = yv(i);
      else
        {
          octave_idx_type k = static_cast<octave_idx_type> (yv(i)) - 1;
          if (k < 0 || k >= K)
            error ("treetrain: Y holds a class index outside 1:K.");
          y[i] = k;
        }
      w[i] = wv(i);
    }

  // Row order along each predictor, sorted once and thereafter partitioned.
  // ord[j] holds the rows in ascending order of predictor j; a node owns the
  // slice [start, stop) of every one of them.
  // A row missing every predictor carries no information and is dropped, as
  // MATLAB drops it: it is absent from the root's size, not merely unsplit.
  std::vector<octave_idx_type> keeprow;
  keeprow.reserve (n);
  for (octave_idx_type i = 0; i < n; i++)
    {
      bool allnan = true;
      for (octave_idx_type j = 0; j < p && allnan; j++)
        if (! std::isnan (X(i, j)))
          allnan = false;
      if (! allnan)
        keeprow.push_back (i);
    }
  const octave_idx_type ne = static_cast<octave_idx_type> (keeprow.size ());

  // Missing values sort to the end of every order, so a node's non-missing
  // rows are always a prefix of its range and stay one through a partition.
  // A bare less-than is not a strict weak ordering once NaN is in the data,
  // which std::sort is entitled to crash on, so the comparator is explicit.
  std::vector<std::vector<octave_idx_type>> ord (p);
  for (octave_idx_type j = 0; j < p; j++)
    {
      ord[j] = keeprow;
      const double *xj = X.data () + j * n;
      std::stable_sort (ord[j].begin (), ord[j].end (),
                        [xj] (octave_idx_type a, octave_idx_type b)
                        {
                          bool na = std::isnan (xj[a]);
                          bool nb = std::isnan (xj[b]);
                          if (na != nb)
                            return nb;
                          if (na)
                            return false;
                          return xj[a] < xj[b];
                        });
    }

  std::vector<Node> nodes;
  nodes.reserve (2 * ne + 2);
  Node root;
  root.left = root.right = root.cutvar = 0;
  root.parent = 0;
  root.cutval = std::numeric_limits<double>::quiet_NaN ();
  root.start = 0;
  root.stop = ne;
  root.nsize = ne;
  root.nweight = 0.0;
  nodes.push_back (root);

  // Class weight of each node, laid out node major.
  std::vector<double> cweight (K, 0.0);
  std::vector<double> nodecw;
  nodecw.reserve ((2 * ne + 2) * K);
  // The unweighted count of each class, which is not the weight once the
  // weights are not all alike, and which the learner reports as ClassCount.
  std::vector<double> nodecn;
  nodecn.reserve ((2 * ne + 2) * K);
  std::vector<double> nodesse, nodemean;
  double rootsse = 0.0;

  std::vector<double> running (K);      // class weight left of a cut
  std::vector<double> tally (K);        // class weight of the whole node
  std::vector<double> other (K);        // and right of a cut
  std::vector<double> present (K);      // and of the rows a predictor has
  std::vector<double> tallyn (K);       // and their unweighted count
  std::vector<octave_idx_type> side (n);
  std::vector<octave_idx_type> buffer (n);
  std::vector<octave_idx_type> held (n);

  octave_idx_type numsplits = 0;

  // The tree grows a layer at a time.  Every node of a layer is searched
  // first, in index order, and only then are the splits made, so that when
  // the layer holds more splittable nodes than MaxNumSplits has left, the
  // ones left unsplit are those whose splits gain the least, as MATLAB
  // documents under tree depth control, and not simply the last ones.
  // Splits are made in index order, so a split's children still follow
  // every node of the layer, which is what makes the numbering breadth first.
  octave_idx_type lfirst = 0;
  while (lfirst < static_cast<octave_idx_type> (nodes.size ()))
    {
      const octave_idx_type llast
        = static_cast<octave_idx_type> (nodes.size ());
      const octave_idx_type lsize = llast - lfirst;
      std::vector<octave_idx_type> lvar (lsize, 0);
      std::vector<double> lval (lsize, 0.0);
      std::vector<double> lgain (lsize, 0.0);
      std::vector<std::vector<double>> lleft (lsize), lright (lsize);

      for (octave_idx_type idx = lfirst; idx < llast; idx++)
        {
          Node& nd = nodes[idx];
          const octave_idx_type start = nd.start;
          const octave_idx_type stop = nd.stop;

          std::fill (tally.begin (), tally.end (), 0.0);
          std::fill (tallyn.begin (), tallyn.end (), 0.0);
          Moments nodemom;
          nodemom.clear ();
          double total = 0.0;
          for (octave_idx_type i = start; i < stop; i++)
            {
              octave_idx_type r = ord[0][i];
              if (isreg)
                nodemom.add (w[r], resp[r]);
              else
                {
                  tally[y[r]] += w[r];
                  tallyn[y[r]] += 1.0;
                }
              total += w[r];
            }
          nd.nsize = stop - start;
          nd.nweight = total;
          if (isreg)
            {
              nodesse.push_back (nodemom.sse ());
              nodemean.push_back (nodemom.mean ());
            }
          else
            {
              nodecw.insert (nodecw.end (), tally.begin (), tally.end ());
              nodecn.insert (nodecn.end (), tallyn.begin (), tallyn.end ());
            }

          if (idx == 0 && isreg)
            rootsse = nodemom.sse ();

          // A node with nothing left to separate is a leaf.  For a classifier
          // that is one class; for a regression it is one distinct response,
          // and QEToler holds it to a share of the whole tree's error rather
          // than to exactly zero.
          bool nothingtosplit;
          if (isreg)
            nothingtosplit = (nodemom.sse () <= qetoler * rootsse);
          else
            {
              octave_idx_type nonempty = 0;
              for (octave_idx_type k = 0; k < K; k++)
                if (tally[k] > 0.0)
                  nonempty++;
              nothingtosplit = (nonempty < 2);
            }

          if (nd.nsize < minparent || numsplits >= maxsplits || nothingtosplit)
            continue;

          octave_idx_type bestvar = 0;
          double bestval = std::numeric_limits<double>::quiet_NaN ();
          double bestgain = 0.0;
          std::vector<double> bestleft, bestright, catl, catr;

          // A searched node tries a fresh subset of the predictors, drawn by a
          // partial shuffle and put back in index order, so a tie between two
          // predictors still goes to the lower index as it does without
          // sampling. A node whose subset holds no valid split is a leaf; no
          // second draw is made.
          if (sampling)
            {
              for (octave_idx_type k = 0; k < nvars; k++)
                std::swap (cand[k], cand[k + draw_below (rng, p - k)]);
              std::sort (cand.begin (), cand.begin () + nvars);
            }

          for (octave_idx_type c = 0; c < nvars; c++)
            {
              const octave_idx_type j = cand[c];
              const double *xj = X.data () + j * n;
              const octave_idx_type *oj = &ord[j][0];

              // The rows missing this predictor sit at the end of the range, so
              // the ones that can be split are the prefix before them.
              octave_idx_type have = stop;
              while (have > start && std::isnan (xj[oj[have - 1]]))
                have--;
              if (have - start < 2 * minleaf)
                continue;

              std::fill (present.begin (), present.end (), 0.0);
              Moments presmom, runmom;
              presmom.clear ();
              double wp = 0.0;
              for (octave_idx_type i = start; i < have; i++)
                {
                  octave_idx_type r = oj[i];
                  if (isreg)
                    presmom.add (w[r], resp[r]);
                  else
                    present[y[r]] += w[r];
                  wp += w[r];
                }
              const double presentimp = isreg ? presmom.sse () / wp
                                              : impurity (present, wp, K, crit);

              // A categorical predictor is split into two sets of levels
              if (! o.iscat.empty () && o.iscat[j])
                {
                  const double g
                    = categorical_split (xj, oj, start, have, y, resp, w, K,
                                         crit, isreg, present, presmom, wp,
                                         presentimp, total, o, catl, catr);
                  if (! catl.empty ()
                      && g > bestgain + std::fabs (bestgain) * GAIN_TIE_TOL)
                    {
                      bestgain = g;
                      bestvar = j + 1;
                      bestval = std::numeric_limits<double>::quiet_NaN ();
                      bestleft = catl;
                      bestright = catr;
                    }
                  continue;
                }

              std::fill (running.begin (), running.end (), 0.0);
              runmom.clear ();
              double wl = 0.0;

              // Position i splits after the i-th row of this node's order, so
              // it leaves i + 1 rows on the left.  Only positions where the
              // value changes are cuts, and each side must keep MinLeaf rows.
              const octave_idx_type last = have - start - 1;
              for (octave_idx_type i = 0; i < last; i++)
                {
                  octave_idx_type r = oj[start + i];
                  if (isreg)
                    runmom.add (w[r], resp[r]);
                  else
                    running[y[r]] += w[r];
                  wl += w[r];

                  const octave_idx_type nleft = i + 1;
                  if (nleft < minleaf || (have - start) - nleft < minleaf)
                    continue;

                  const double xa = xj[r];
                  const double xb = xj[oj[start + i + 1]];
                  if (! (xa < xb))
                    continue;

                  const double wr = wp - wl;
                  double impl, impr;
                  if (isreg)
                    {
                      Moments rest;
                      rest.w = presmom.w - runmom.w;
                      rest.wy = presmom.wy - runmom.wy;
                      rest.wyy = presmom.wyy - runmom.wyy;
                      impl = runmom.sse () / wl;
                      impr = rest.sse () / wr;
                    }
                  else
                    {
                      for (octave_idx_type k = 0; k < K; k++)
                        other[k] = present[k] - running[k];
                      impl = impurity (running, wl, K, crit);
                      impr = impurity (other, wr, K, crit);
                    }

                  // The rows that are missing this predictor cannot be placed
                  // by it, so the split is judged on the rows it can place and
                  // then scaled by their share of the node.  Without that
                  // scaling a predictor known for a handful of rows would beat
                  // one known for all of them.
                  const double gain = (wp * presentimp - wl * impl - wr * impr)
                                      / total;

                  if (gain > bestgain + std::fabs (bestgain) * GAIN_TIE_TOL)
                    {
                      bestgain = gain;
                      bestvar = j + 1;
                      bestval = (xa + xb) / 2.0;
                      bestleft.clear ();
                      bestright.clear ();
                    }
                }
            }

          if (bestvar == 0 || bestgain <= 0.0)
            continue;

          // The gain is per unit of the node's weight; weighted by it, it is
          // the reduction in the tree's risk, which is what ranks nodes of
          // different sizes against one another.
          lvar[idx - lfirst] = bestvar;
          lval[idx - lfirst] = bestval;
          lgain[idx - lfirst] = bestgain * total;
          lleft[idx - lfirst] = bestleft;
          lright[idx - lfirst] = bestright;
        }

      // Keep the most successful splits the budget allows.  A stable sort
      // leaves the earlier node ahead on a tie.
      std::vector<octave_idx_type> found;
      for (octave_idx_type i = 0; i < lsize; i++)
        if (lvar[i] != 0)
          found.push_back (i);
      const octave_idx_type budget = maxsplits - numsplits;
      if (static_cast<octave_idx_type> (found.size ()) > budget)
        {
          std::stable_sort (found.begin (), found.end (),
                            [&lgain] (octave_idx_type a, octave_idx_type b)
                            { return lgain[a] > lgain[b]; });
          for (std::size_t i = budget; i < found.size (); i++)
            lvar[found[i]] = 0;
        }

      for (octave_idx_type idx = lfirst; idx < llast; idx++)
        {
          if (lvar[idx - lfirst] == 0)
            continue;

          // Mark which side each row takes, then partition every predictor's
          // order in place, keeping the rows of each side in the order they
          // were. A row missing the chosen predictor cannot be placed on either
          // side, so it descends to neither child and is left in the tail of
          // the range, counted at this node and nowhere below it.  Prediction
          // stops such a row at this node for the same reason.
          const octave_idx_type start = nodes[idx].start;
          const octave_idx_type stop = nodes[idx].stop;
          const octave_idx_type bestvar = lvar[idx - lfirst];
          const double bestval = lval[idx - lfirst];
          const std::vector<double>& catleft = lleft[idx - lfirst];
          const bool bycat = ! catleft.empty ();
          const double *xb = X.data () + (bestvar - 1) * n;
          octave_idx_type nleft = 0, nright = 0;
          for (octave_idx_type i = start; i < stop; i++)
            {
              octave_idx_type r = ord[0][i];
              if (std::isnan (xb[r]))
                side[r] = 2;
              else if (bycat ? std::binary_search (catleft.begin (),
                                                   catleft.end (), xb[r])
                             : (xb[r] < bestval))
                {
                  side[r] = 0;
                  nleft++;
                }
              else
                {
                  side[r] = 1;
                  nright++;
                }
            }

          for (octave_idx_type j = 0; j < p; j++)
            {
              octave_idx_type a = 0, b = 0, c = 0;
              for (octave_idx_type i = start; i < stop; i++)
                {
                  octave_idx_type r = ord[j][i];
                  if (side[r] == 0)
                    ord[j][start + a++] = r;
                  else if (side[r] == 1)
                    buffer[b++] = r;
                  else
                    held[c++] = r;
                }
              for (octave_idx_type i = 0; i < b; i++)
                ord[j][start + a + i] = buffer[i];
              for (octave_idx_type i = 0; i < c; i++)
                ord[j][start + a + b + i] = held[i];
            }

          Node kid;
          kid.left = kid.right = kid.cutvar = 0;
          kid.cutval = std::numeric_limits<double>::quiet_NaN ();
          kid.parent = idx + 1;
          kid.nsize = 0;
          kid.nweight = 0.0;

          nodes[idx].cutvar = bestvar;
          nodes[idx].cutval = bestval;
          nodes[idx].catleft = catleft;
          nodes[idx].catright = lright[idx - lfirst];
          nodes[idx].left = static_cast<octave_idx_type> (nodes.size ()) + 1;
          nodes[idx].right = static_cast<octave_idx_type> (nodes.size ()) + 2;

          kid.start = start;
          kid.stop = start + nleft;
          nodes.push_back (kid);
          kid.start = start + nleft;
          kid.stop = start + nleft + nright;
          nodes.push_back (kid);

          numsplits++;
        }

      lfirst = llast;
    }

  octave_idx_type nn = static_cast<octave_idx_type> (nodes.size ());

  // Collapse a pair of leaves back into their parent when splitting it did
  // not lower the misclassification risk, deepest first so that a merge can
  // expose another.  The risk is the weight the node's majority class does
  // not hold, which is the error criterion rather than the split criterion.
  // A node's class weight is the sum of its rows, and a parent's rows are its
  // two children's rows.  Summed in those two different orders the totals
  // differ in the last bits, which is enough to decide the merge test below:
  // on eight thousand rows every pair that ought to merge missed by about
  // 1e-19, against an observation weight of 1.25e-4.  Recomputing a branch's
  // totals from its children, deepest first, makes the two sides of that test
  // exact rather than nearly equal, and no tolerance is needed.
  // Only where a split placed every row is a parent the sum of its children,
  // and only there can its totals be recomputed from them.  A node that held
  // rows back keeps its own, and the merge test below falls back to the form
  // that does not assume the two agree.
  std::vector<bool> whole (nn, false);
  for (octave_idx_type idx = nn - 1; idx >= 0; idx--)
    {
      octave_idx_type l = nodes[idx].left;
      if (l == 0)
        continue;
      octave_idx_type r = nodes[idx].right;
      if (nodes[idx].nsize != nodes[l-1].nsize + nodes[r-1].nsize)
        continue;
      whole[idx] = true;
      if (! isreg)
        for (octave_idx_type k = 0; k < K; k++)
          nodecw[idx * K + k] = nodecw[(l - 1) * K + k]
                                + nodecw[(r - 1) * K + k];
      nodes[idx].nweight = nodes[l-1].nweight + nodes[r-1].nweight;
    }

  std::vector<bool> keep;
  if (mergeleaves)
    {
      bool merged = true;
      while (merged)
        {
          merged = false;
          for (octave_idx_type idx = nn - 1; idx >= 0; idx--)
            {
              octave_idx_type l = nodes[idx].left;
              octave_idx_type r = nodes[idx].right;
              if (l == 0 || nodes[l-1].left != 0 || nodes[r-1].left != 0)
                continue;

              // Splitting is kept only when it lowers the misclassification
              // risk.  Written as risks the two sides share the node weight,
              // and subtracting it loses the very digits the comparison turns
              // on; cancelled algebraically the test is over class weights
              // alone, and since the parent's are the sums of the children's
              // it is exact when it should be.
              bool collapse;
              if (isreg)
                {
                  // Splitting can only lower the squared error, so a
                  // regression merge is the degenerate case where it lowered
                  // it by nothing at all.
                  collapse = (nodesse[l-1] + nodesse[r-1] >= nodesse[idx]);
                }
              else
                {
                  double pbest = 0.0, lbest = 0.0, rbest = 0.0;
                  for (octave_idx_type k = 0; k < K; k++)
                    {
                      pbest = std::max (pbest, nodecw[idx * K + k]);
                      lbest = std::max (lbest, nodecw[(l - 1) * K + k]);
                      rbest = std::max (rbest, nodecw[(r - 1) * K + k]);
                    }
                  if (whole[idx])
                    collapse = (pbest >= lbest + rbest);
                  else
                    collapse = ((nodes[l-1].nweight - lbest)
                                + (nodes[r-1].nweight - rbest)
                                >= nodes[idx].nweight - pbest);
                }

              if (collapse)
                {
                  nodes[idx].left = nodes[idx].right = 0;
                  nodes[idx].cutvar = 0;
                  nodes[idx].cutval = std::numeric_limits<double>::quiet_NaN ();
                  nodes[idx].catleft.clear ();
                  nodes[idx].catright.clear ();
                  merged = true;
                }
            }
        }
    }

  // Drop what a merge orphaned and renumber what is left, a parent keeping
  // its place ahead of its children.
  keep.assign (nn, false);
  keep[0] = true;
  for (octave_idx_type idx = 0; idx < nn; idx++)
    if (keep[idx] && nodes[idx].left != 0)
      {
        keep[nodes[idx].left - 1] = true;
        keep[nodes[idx].right - 1] = true;
      }

  std::vector<octave_idx_type> newidx (nn, 0);
  octave_idx_type out = 0;
  for (octave_idx_type idx = 0; idx < nn; idx++)
    if (keep[idx])
      newidx[idx] = ++out;

  Matrix children (out, 2, 0.0);
  ColumnVector parent (out, 0.0);
  ColumnVector cutvar (out, 0.0);
  ColumnVector cutval (out, 0.0);
  boolMatrix isbranch (out, 1, false);
  ColumnVector nodesize (out, 0.0);
  ColumnVector nodeweight (out, 0.0);
  Matrix classweight (out, isreg ? 0 : K, 0.0);
  Matrix classcount (out, isreg ? 0 : K, 0.0);
  ColumnVector nodemeanout (isreg ? out : 0, 0.0);
  ColumnVector nodeerror (isreg ? out : 0, 0.0);
  Cell cutcats (out, 2);

  for (octave_idx_type idx = 0; idx < nn; idx++)
    {
      if (! keep[idx])
        continue;
      octave_idx_type i = newidx[idx] - 1;
      children(i, 0) = (nodes[idx].left == 0) ? 0 : newidx[nodes[idx].left - 1];
      children(i, 1) = (nodes[idx].right == 0)
                       ? 0 : newidx[nodes[idx].right - 1];
      parent(i) = (nodes[idx].parent == 0) ? 0 : newidx[nodes[idx].parent - 1];
      cutvar(i) = nodes[idx].cutvar;
      cutval(i) = (nodes[idx].cutvar == 0)
                  ? std::numeric_limits<double>::quiet_NaN ()
                  : nodes[idx].cutval;
      isbranch(i, 0) = (nodes[idx].cutvar != 0);
      cutcats(i, 0) = Matrix ();
      cutcats(i, 1) = Matrix ();
      if (nodes[idx].cutvar != 0 && ! nodes[idx].catleft.empty ())
        {
          RowVector a (nodes[idx].catleft.size ());
          RowVector b (nodes[idx].catright.size ());
          for (std::size_t q = 0; q < nodes[idx].catleft.size (); q++)
            a(q) = nodes[idx].catleft[q];
          for (std::size_t q = 0; q < nodes[idx].catright.size (); q++)
            b(q) = nodes[idx].catright[q];
          cutcats(i, 0) = a;
          cutcats(i, 1) = b;
        }
      nodesize(i) = nodes[idx].nsize;
      nodeweight(i) = nodes[idx].nweight;
      if (isreg)
        {
          nodemeanout(i) = nodemean[idx];
          nodeerror(i) = (nodes[idx].nweight > 0.0)
                         ? nodesse[idx] / nodes[idx].nweight : 0.0;
        }
      else
        for (octave_idx_type k = 0; k < K; k++)
          {
            classweight(i, k) = nodecw[idx * K + k];
            classcount(i, k) = nodecn[idx * K + k];
          }
    }

  octave_scalar_map T;
  T.assign ("Children", children);
  T.assign ("Parent", parent);
  T.assign ("CutPredictorIndex", cutvar);
  T.assign ("CutPoint", cutval);
  T.assign ("CutCategories", cutcats);
  T.assign ("IsBranchNode", isbranch);
  T.assign ("NodeSize", nodesize);
  T.assign ("NodeWeight", nodeweight);
  T.assign ("ClassWeight", classweight);
  T.assign ("ClassCount", classcount);
  T.assign ("NodeMean", nodemeanout);
  T.assign ("NodeError", nodeerror);
  T.assign ("NumNodes", octave_value (out));

  return T;
}

// Send each row of X down the tree and report where it came to rest and what
// that node answers with.  A row is stopped by a node whose split predictor it
// is missing, exactly as growth held such a row back rather than sending it to
// a child; the two are the same rule, which is why they live in one file.  A
// categorical cut sends a row by the set its level is in, and stops a row
// whose level is in neither, a level the node never saw, as MATLAB does.
static void
tree_descend (const Matrix& X, const Matrix& children,
              const ColumnVector& cutvar, const ColumnVector& cutpoint,
              const Matrix& value, Matrix& V, ColumnVector& node,
              const Cell *cats = nullptr)
{
  const octave_idx_type n = X.rows ();
  const octave_idx_type m = value.columns ();
  const octave_idx_type nn = children.rows ();

  V.resize (n, m);
  node.resize (n);

  std::vector<std::vector<double>> catl (nn), catr (nn);
  if (cats)
    for (octave_idx_type i = 0; i < nn; i++)
      {
        if ((*cats)(i, 0).isempty ())
          continue;
        const NDArray a = (*cats)(i, 0).array_value ();
        const NDArray b = (*cats)(i, 1).array_value ();
        catl[i].assign (a.data (), a.data () + a.numel ());
        catr[i].assign (b.data (), b.data () + b.numel ());
        std::sort (catl[i].begin (), catl[i].end ());
        std::sort (catr[i].begin (), catr[i].end ());
      }

  for (octave_idx_type i = 0; i < n; i++)
    {
      octave_idx_type at = 1;
      while (children(at - 1, 0) != 0)
        {
          const octave_idx_type v
            = static_cast<octave_idx_type> (cutvar(at - 1));
          const double x = X(i, v - 1);
          if (std::isnan (x))
            break;
          if (! catl[at - 1].empty ())
            {
              const octave_idx_type from = at;
              if (std::binary_search (catl[from - 1].begin (),
                                      catl[from - 1].end (), x))
                at = static_cast<octave_idx_type> (children(from - 1, 0));
              else if (std::binary_search (catr[from - 1].begin (),
                                           catr[from - 1].end (), x))
                at = static_cast<octave_idx_type> (children(from - 1, 1));
              else
                break;
              continue;
            }
          at = static_cast<octave_idx_type> ((x < cutpoint(at - 1))
                                             ? children(at - 1, 0)
                                             : children(at - 1, 1));
        }
      node(i) = at;
      for (octave_idx_type k = 0; k < m; k++)
        V(i, k) = value(at - 1, k);
    }
}
