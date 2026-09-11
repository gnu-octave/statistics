/*
Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This file is part of the statistics package for GNU Octave.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
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
// Everything the caller settles before the engine starts.  A classifier and a
// regression differ here in three fields and nowhere else.
struct TreeOpts
{
  octave_idx_type K, minparent, minleaf, maxsplits;
  bool mergeleaves;
  double qetoler;
  Criterion crit;
};

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

  // Nodes are visited in index order, so a split's children are visited later
  // in the same pass, which is what makes the numbering breadth first.
  for (octave_idx_type idx = 0;
       idx < static_cast<octave_idx_type> (nodes.size ()); idx++)
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
      // that is one class; for a regression it is one distinct response, and
      // QEToler holds it to a share of the whole tree's error rather than to
      // exactly zero.
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

      for (octave_idx_type j = 0; j < p; j++)
        {
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

          std::fill (running.begin (), running.end (), 0.0);
          runmom.clear ();
          double wl = 0.0;

          // Position i splits after the i-th row of this node's order, so it
          // leaves i + 1 rows on the left.  Only positions where the value
          // changes are cuts, and each side must keep MinLeaf rows.
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

              // The rows that are missing this predictor cannot be placed by
              // it, so the split is judged on the rows it can place and then
              // scaled by their share of the node.  Without that scaling a
              // predictor known for a handful of rows would beat one known
              // for all of them.
              const double gain = (wp * presentimp - wl * impl - wr * impr)
                                  / total;

              if (gain > bestgain + std::fabs (bestgain) * GAIN_TIE_TOL)
                {
                  bestgain = gain;
                  bestvar = j + 1;
                  bestval = (xa + xb) / 2.0;
                }
            }
        }

      if (bestvar == 0 || bestgain <= 0.0)
        continue;

      // Mark which side each row takes, then partition every predictor's
      // order in place, keeping the rows of each side in the order they were.
      // A row missing the chosen predictor cannot be placed on either side,
      // so it descends to neither child and is left in the tail of the range,
      // counted at this node and nowhere below it.  Prediction stops such a
      // row at this node for the same reason.
      const double *xb = X.data () + (bestvar - 1) * n;
      octave_idx_type nleft = 0, nright = 0;
      for (octave_idx_type i = start; i < stop; i++)
        {
          octave_idx_type r = ord[0][i];
          if (std::isnan (xb[r]))
            side[r] = 2;
          else if (xb[r] < bestval)
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
// a child; the two are the same rule, which is why they live in one file.
static void
tree_descend (const Matrix& X, const Matrix& children,
              const ColumnVector& cutvar, const ColumnVector& cutpoint,
              const Matrix& value, Matrix& V, ColumnVector& node)
{
  const octave_idx_type n = X.rows ();
  const octave_idx_type m = value.columns ();

  V.resize (n, m);
  node.resize (n);

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
          at = static_cast<octave_idx_type> ((x < cutpoint(at - 1))
                                             ? children(at - 1, 0)
                                             : children(at - 1, 1));
        }
      node(i) = at;
      for (octave_idx_type k = 0; k < m; k++)
        V(i, k) = value(at - 1, k);
    }
}
