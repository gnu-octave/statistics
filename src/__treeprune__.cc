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

#include <cmath>
#include <string>
#include <vector>
#include <octave/oct.h>

// How close two risks must be to count as equal.  Links equal in exact
// arithmetic differ in their last bits once they have been through a
// division, and treating them as distinct would split one pruning level into
// several.  It is a tolerance on risks, not the one the split gains use.
static const double RISK_TIE_TOL = 1e-12;

// Cost complexity pruning, by weakest link.  A branch node's link is the risk
// it would take on as a leaf, less the risk its subtree carries now, spread
// over the leaves the subtree would give up.  The branch with the smallest
// link is pruned, the sequence repeats on what is left, and the level a node
// is pruned at is its place in that sequence.  Leaves never carry a level, so
// a tree of B branches gives B levels and B + 1 alphas, the first of which is
// zero and stands for the unpruned tree.
//
// RISK is what each node would carry as a leaf and HELD what a branch carries
// for the rows that stop there, on whatever scale the caller measures by.
// Taking both rather than deriving them is what lets a tree be sequenced
// after it was grown, and on a scale the engine knows nothing of: a caller
// holding a cost matrix passes the risk that matrix defines.
static void
tree_prune_sequence (const Matrix& children, const ColumnVector& parent,
                     const std::vector<double>& risk,
                     const std::vector<double>& held,
                     ColumnVector& prunelist, ColumnVector& prunealpha)
{
  const octave_idx_type out = children.rows ();
  prunelist = ColumnVector (out, 0.0);
  prunealpha = ColumnVector (0, 0.0);
  if (out <= 0)
    return;

  std::vector<octave_idx_type> kidl (out), kidr (out);
  for (octave_idx_type i = 0; i < out; i++)
    {
      kidl[i] = static_cast<octave_idx_type> (children(i, 0));
      kidr[i] = static_cast<octave_idx_type> (children(i, 1));
    }

  std::vector<double> subrisk (out);
  std::vector<octave_idx_type> subleaves (out);
  std::vector<bool> reach (out);
  // The branches given up at no cost, which are no part of the sequence
  std::vector<bool> freed (out, false);
  std::vector<double> alphas;
  octave_idx_type level = 0;

  while (kidl[0] != 0)
    {
      // A pruned branch takes its whole subtree with it, so only what is
      // still reachable from the root can be a candidate.  Left in, those
      // orphans go on offering links of their own and split one level
      // into several.
      std::fill (reach.begin (), reach.end (), false);
      reach[0] = true;
      for (octave_idx_type i = 0; i < out; i++)
        if (reach[i] && kidl[i] != 0)
          {
            reach[kidl[i]-1] = true;
            reach[kidr[i]-1] = true;
          }

      // Subtree risk and leaf count, deepest first.
      for (octave_idx_type i = out - 1; i >= 0; i--)
        {
          if (kidl[i] == 0)
            {
              subrisk[i] = risk[i];
              subleaves[i] = 1;
            }
          else
            {
              subrisk[i] = subrisk[kidl[i]-1]
                           + subrisk[kidr[i]-1] + held[i];
              subleaves[i] = subleaves[kidl[i]-1] + subleaves[kidr[i]-1];
            }
        }

      bool any = false;
      double weakest = 0.0;
      for (octave_idx_type i = 0; i < out; i++)
        {
          if (kidl[i] == 0 || ! reach[i])
            continue;
          double link = (risk[i] - subrisk[i]) / (subleaves[i] - 1);
          if (! any || link < weakest)
            {
              weakest = link;
              any = true;
            }
        }

      if (! any)
        break;

      // A subtree that costs nothing to give up is no step of the
      // sequence.  It is what merging leaves would have removed, and it
      // survives to be seen here only when the caller asked for a
      // sequence without a merge.  Such a branch is given up all the
      // same, so that the branches above it are priced on what is left,
      // but it opens no level and takes no alpha: measured on R2024a,
      // where an unmerged iris tree of eleven nodes carries the merged
      // tree's five alphas rather than six.
      const double tol = RISK_TIE_TOL * std::max (std::fabs (risk[0]), 1.0);
      const bool free = (weakest <= tol);

      // Every branch whose link is the weakest goes at this level, not
      // just one of them, and links equal in exact arithmetic can differ
      // in their last bits here as split gains do.
      if (! free)
        {
          level++;
          alphas.push_back (weakest);
        }
      const double cut = weakest + std::fabs (weakest) * RISK_TIE_TOL + tol;
      for (octave_idx_type i = 0; i < out; i++)
        {
          if (kidl[i] == 0 || ! reach[i])
            continue;
          if ((risk[i] - subrisk[i]) / (subleaves[i] - 1) <= cut)
            {
              if (free)
                freed[i] = true;
              else
                prunelist(i) = level;
              kidl[i] = kidr[i] = 0;
            }
        }
    }

  // A branch inside a subtree given up at no cost is no part of the
  // sequence either, its ancestor having left it at no level.  A parent
  // always precedes its children, so one forward pass carries the mark
  // down.
  for (octave_idx_type i = 1; i < out; i++)
    {
      octave_idx_type a = static_cast<octave_idx_type> (parent(i));
      if (a > 0 && freed[a - 1])
        freed[i] = true;
    }

  // A branch that lost an ancestor never came up for pruning on its own
  // account, but it stopped being a branch when that ancestor went, and
  // that is the level it carries.
  for (octave_idx_type i = 0; i < out; i++)
    {
      if (children(i, 0) == 0 || prunelist(i) != 0 || freed[i])
        continue;
      octave_idx_type a = static_cast<octave_idx_type> (parent(i));
      while (a > 0 && prunelist(a - 1) == 0)
        a = static_cast<octave_idx_type> (parent(a - 1));
      if (a > 0)
        prunelist(i) = prunelist(a - 1);
    }

  prunealpha.resize (level + 1);
  prunealpha(0) = 0.0;
  for (octave_idx_type i = 0; i < level; i++)
    prunealpha(i + 1) = alphas[i];
}

// The cost complexity pruning sequence of a tree that has already been grown.
//
// treetrain sequences the tree it grows, but a tree outlives its data: prune
// collapses nodes and the tree that is left has a sequence of its own, and a
// model read back from a file has no data at all.  Taking the risk rather
// than the observations is what serves both, and it is also what lets a
// classifier measure risk on a cost matrix this engine knows nothing about.

DEFUN_DLD (__treeprune__, args, ,
           "-*- texinfo -*-\n\
@deftypefn {statistics} {[@var{PruneList}, @var{PruneAlpha}] =} \
__treeprune__ (@var{Children}, @var{Parent}, @var{risk}, @var{held})\n\
\n\
Cost complexity pruning sequence of a grown tree, by weakest link.\n\
Internal; called by @code{ClassificationTree} and @code{RegressionTree} and\n\
not meant to be used directly.\n\
\n\
A branch node's link is the risk it would take on as a leaf, less the risk\n\
its subtree carries now, spread over the leaves the subtree would give up.\n\
The weakest link is pruned, the sequence repeats on what is left, and the\n\
level a node is pruned at is its place in that sequence.  Leaves never carry\n\
a level, so a tree of @math{B} branch nodes gives @math{B} levels at most and\n\
one more alpha than levels, the first of which is zero and stands for the\n\
unpruned tree.\n\
\n\
@var{Children} is the @math{Nx2} table of child node numbers, zero on a leaf,\n\
and @var{Parent} the node number each node hangs from, zero at the root.\n\
\n\
@var{risk} is the risk each node would carry as a leaf, on whatever scale the\n\
caller measures by: the expected misclassification cost for a classifier and\n\
the squared error about the node's mean for a regression.  Nothing here\n\
interprets it, so a cost matrix enters through this argument alone.\n\
\n\
@var{held} is the risk a node carries on account of the observations that\n\
stop there, missing the predictor it cuts on.  Those observations are in\n\
neither child, so a subtree's risk is its children's plus this, and leaving\n\
it out overstates every link above a node that holds rows back.\n\
\n\
A subtree that costs nothing to give up is no step of the sequence.  It is\n\
what merging leaves would have removed, and MATLAB records neither a level\n\
nor an alpha for it.\n\
\n\
@end deftypefn")
{
  if (args.length () != 4)
    print_usage ();

  if (! args(0).isnumeric () || args(0).iscomplex () || args(0).isempty ()
      || args(0).columns () != 2)
    error ("__treeprune__: Children must be a non-empty N-by-2 real matrix.");

  const Matrix children = args(0).matrix_value ();
  const octave_idx_type n = children.rows ();

  if (! args(1).isnumeric () || args(1).iscomplex ()
      || args(1).numel () != n)
    error ("__treeprune__: Parent must hold one node number per node.");
  if (! args(2).isnumeric () || args(2).iscomplex ()
      || args(2).numel () != n)
    error ("__treeprune__: risk must hold one value per node.");
  if (! args(3).isnumeric () || args(3).iscomplex ()
      || args(3).numel () != n)
    error ("__treeprune__: held must hold one value per node.");

  const ColumnVector parent = args(1).column_vector_value ();
  const ColumnVector riskv = args(2).column_vector_value ();
  const ColumnVector heldv = args(3).column_vector_value ();

  // A child number outside the table would index past the end of every
  // vector below, so it is refused here rather than trusted.
  for (octave_idx_type i = 0; i < n; i++)
    {
      const double l = children(i, 0);
      const double r = children(i, 1);
      if (l != std::floor (l) || r != std::floor (r) || l < 0 || r < 0
          || l > n || r > n || (l == 0) != (r == 0))
        error ("__treeprune__: Children must hold a node number in each row, \
or zero in both on a leaf.");
      if (parent(i) != std::floor (parent(i)) || parent(i) < 0
          || parent(i) > n)
        error ("__treeprune__: Parent must hold a node number, or zero at \
the root.");
    }

  std::vector<double> risk (n), held (n);
  for (octave_idx_type i = 0; i < n; i++)
    {
      risk[i] = riskv(i);
      held[i] = heldv(i);
    }

  ColumnVector prunelist, prunealpha;
  tree_prune_sequence (children, parent, risk, held, prunelist, prunealpha);

  return ovl (prunelist, prunealpha);
}

/*
%!test
%! ## The link of a root whose two leaves cost less than it does
%! C = [2, 3; 0, 0; 0, 0];
%! P = [0; 1; 1];
%! [L, A] = __treeprune__ (C, P, [1; 0.4; 0.4], zeros (3, 1));
%! assert_equal (L', [1, 0, 0]);
%! assert_equal (A', [0, 0.2], 1e-14);

%!test
%! ## What a node holds back is part of its subtree's risk
%! ## The same tree, charged 0.1 for rows that reached neither child, gives
%! ## a link of 0.1 where it gave 0.2.
%! C = [2, 3; 0, 0; 0, 0];
%! P = [0; 1; 1];
%! [L, A] = __treeprune__ (C, P, [1; 0.4; 0.4], [0.1; 0; 0]);
%! assert_equal (L', [1, 0, 0]);
%! assert_equal (A', [0, 0.1], 1e-14);

%!test
%! ## A subtree costing nothing takes no level and no alpha
%! C = [2, 3; 0, 0; 0, 0];
%! P = [0; 1; 1];
%! [L, A] = __treeprune__ (C, P, [0.8; 0.4; 0.4], zeros (3, 1));
%! assert_equal (L', [0, 0, 0]);
%! assert_equal (A', 0);

%!test
%! ## Deeper branches go at their own level, weakest first
%! C = [2, 3; 4, 5; 0, 0; 0, 0; 0, 0];
%! P = [0; 1; 1; 2; 2];
%! [L, A] = __treeprune__ (C, P, [1; 0.5; 0.3; 0.2; 0.2], zeros (5, 1));
%! assert_equal (L', [2, 1, 0, 0, 0]);
%! assert_equal (A', [0, 0.1, 0.2], 1e-14);

%!test
%! ## A branch that lost an ancestor carries the ancestor's level
%! ## Node 2 has the stronger link and never comes up on its own account,
%! ## but it stopped being a branch when node 1 went.
%! C = [2, 3; 4, 5; 0, 0; 0, 0; 0, 0];
%! P = [0; 1; 1; 2; 2];
%! [L, A] = __treeprune__ (C, P, [1; 0.9; 0.3; 0.25; 0.25], zeros (5, 1));
%! assert_equal (L', [1, 1, 0, 0, 0]);
%! assert_equal (A', [0, 0.1], 1e-14);

%!test
%! ## The scale is the caller's and nothing here interprets it
%! C = [2, 3; 0, 0; 0, 0];
%! P = [0; 1; 1];
%! [~, A] = __treeprune__ (C, P, [10; 4; 4], zeros (3, 1));
%! assert_equal (A', [0, 2], 1e-13);

%!test
%! ## A tree of one node has no branch to prune
%! [L, A] = __treeprune__ ([0, 0], 0, 1, 0);
%! assert_equal (L, 0);
%! assert_equal (A, 0);

## The sequence of a tree treetrain grew, on the risk that engine measures
## by: misclassified weight for a classifier, squared error for a regression.
## A learner with a cost matrix builds its own risk and does not come here.
%!function [L, A] = seq (T)
%!  if (isempty (T.ClassWeight))
%!    risk = T.NodeWeight .* T.NodeError;
%!  else
%!    risk = T.NodeWeight - max (T.ClassWeight, [], 2);
%!  endif
%!  held = zeros (T.NumNodes, 1);
%!  br = find (T.Children(:,1) > 0);
%!  if (! isempty (br))
%!    kd = T.Children(br,:);
%!    wh = T.NodeWeight(br) - T.NodeWeight(kd(:,1)) - T.NodeWeight(kd(:,2));
%!    held(br) = wh .* risk(br) ./ T.NodeWeight(br);
%!  endif
%!  [L, A] = __treeprune__ (T.Children, T.Parent, risk, held);
%!endfunction

%!test
%! ## Levels and alphas of the merged iris tree.  Measured on R2024a.
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', 'MergeLeaves', true);
%! [L, A] = seq (treetrain (meas, y, ones (150, 1) / 150, o));
%! assert_equal (L', [4, 0, 3, 2, 0, 1, 0, 0, 0]);
%! assert_equal (A', [0, 1/150, 2/150, 44/150, 50/150], 1e-12);

%!test
%! ## A deeper tree, every branch pruned at its own level
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 2, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', 'MergeLeaves', true);
%! L = seq (treetrain (meas, y, ones (150, 1) / 150, o));
%! assert_equal (L', [5, 0, 4, 3, 1, 2, 2, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0]);

%!test
%! ## A subtree that costs nothing to give up opens no level of the sequence
%! ## The pair of leaves that merging would have removed survives here, and
%! ## giving it up costs nothing, so the eleven node tree carries the merged
%! ## tree's five alphas rather than six.  Measured on R2024a.
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', 'MergeLeaves', false);
%! [L, A] = seq (treetrain (meas, y, ones (150, 1) / 150, o));
%! assert_equal (L', [4, 0, 3, 2, 0, 1, 0, 0, 0, 0, 0]);
%! assert_equal (A', [0, 1/150, 2/150, 44/150, 50/150], 1e-12);

%!test
%! ## A node holding rows back pays for them in the pruning sequence
%! ## Twenty rows have no fourth predictor, so the node cutting on it sends
%! ## them to neither child.  A subtree's risk is its children's plus what
%! ## the node holds back; without it the sequence came out a level short.
%! ## Measured on R2024a, the fit fitctree makes on this fixture.
%! load fisheriris
%! y = grp2idx (species);
%! x = meas;
%! x(51:70, 4) = NaN;
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', 'MergeLeaves', true);
%! [L, A] = seq (treetrain (x, y, ones (150, 1) / 150, o));
%! assert_equal (L', [4, 0, 3, 1, 2, 0, 0, 0, 0]);
%! assert_equal (A', [0, 0.00385185185185185, 0.00593939393939394, ...
%!                    0.286666666666666, 0.333333333333333], 1e-14);

%!test
%! ## A held back node keeps the free subtree rule honest
%! ## Charged for what it holds back, such a node gives a positive link
%! ## where an uncharged one gave zero, so the free subtree left by
%! ## MergeLeaves off must still be the one recognised.  Measured on R2024a.
%! load fisheriris
%! y = grp2idx (species);
%! x = meas;
%! x(51:70, 4) = NaN;
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', 'MergeLeaves', false);
%! [L, A] = seq (treetrain (x, y, ones (150, 1) / 150, o));
%! assert_equal (L', [4, 0, 3, 1, 2, 0, 0, 0, 0, 0, 0]);
%! assert_equal (A', [0, 0.00385185185185185, 0.00593939393939394, ...
%!                    0.286666666666666, 0.333333333333333], 1e-14);

%!test
%! ## A regression tree, sequenced on the squared error about each mean
%! g = mod ((1:20)', 3);
%! x = [(1:20)', g];
%! y = [ones(10,1); 5 * ones(10,1)];
%! o = struct ('NumClasses', 1, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 19, 'SplitCriterion', 'mse', ...
%!             'MergeLeaves', true, 'QEToler', 1e-6);
%! L = seq (treetrain (x, y, ones (20, 1) / 20, o));
%! assert_equal (L', [1, 0, 0]);

%!test
%! ## A regression node holding a row back pays for it too
%! ## One carsmall row has no horsepower and node 2 cuts on horsepower, so
%! ## that row stops there and is in neither child.  Charging it the node's
%! ## own error puts this alpha at 5.99 rather than 6.32.  Measured on
%! ## R2024a, the fit fitrtree makes on this fixture.
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! ok = ! isnan (MPG);
%! o = struct ('NumClasses', 1, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 93, 'SplitCriterion', 'mse', ...
%!             'MergeLeaves', true, 'QEToler', 1e-6);
%! [L, A] = seq (treetrain (X(ok, :), MPG(ok), ones (94, 1) / 94, o));
%! assert_equal (L(1:5)', [17, 16, 14, 15, 13]);
%! assert_equal (numel (A), 18);
%! assert_equal (A(17), 5.99325416896717, 1e-12);
%! assert_equal (A(18), 41.4954735525515, 1e-11);

## Test input validation
%!error <Invalid call> __treeprune__ ([0, 0], 0, 1)
%!error <__treeprune__: Children must be a non-empty N-by-2 real matrix.> ...
%! __treeprune__ ([], [], [], [])
%!error <__treeprune__: Children must be a non-empty N-by-2 real matrix.> ...
%! __treeprune__ ([0, 0, 0], 0, 1, 0)
%!error <__treeprune__: Parent must hold one node number per node.> ...
%! __treeprune__ ([0, 0], [0; 0], 1, 0)
%!error <__treeprune__: risk must hold one value per node.> ...
%! __treeprune__ ([0, 0], 0, [1; 1], 0)
%!error <__treeprune__: held must hold one value per node.> ...
%! __treeprune__ ([0, 0], 0, 1, [0; 0])
%!error <__treeprune__: Children must hold a node number in each row, or zero in both on a leaf.> ...
%! __treeprune__ ([2, 0; 0, 0], [0; 1], [1; 0.5], [0; 0])
%!error <__treeprune__: Children must hold a node number in each row, or zero in both on a leaf.> ...
%! __treeprune__ ([2, 9; 0, 0], [0; 1], [1; 0.5], [0; 0])
%!error <__treeprune__: Parent must hold a node number, or zero at the root.> ...
%! __treeprune__ ([2, 3; 0, 0; 0, 0], [0; 9; 1], [1; 0.4; 0.4], zeros (3, 1))
*/
