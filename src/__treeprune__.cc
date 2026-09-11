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
#include "tree.cpp"

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

%!test
%! ## The same sequence treetrain builds for the tree it grows
%! ## One entry point or the other, the answer is one implementation.
%! load fisheriris
%! y = grp2idx (species);
%! x = meas;
%! x(51:70, 4) = NaN;
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'Prune', true);
%! T = treetrain (x, y, ones (150, 1) / 150, o);
%! risk = T.NodeWeight - max (T.ClassWeight, [], 2);
%! held = zeros (T.NumNodes, 1);
%! br = find (T.Children(:,1) > 0);
%! kids = T.Children(br,:);
%! wh = T.NodeWeight(br) - T.NodeWeight(kids(:,1)) - T.NodeWeight(kids(:,2));
%! held(br) = wh .* risk(br) ./ T.NodeWeight(br);
%! [L, A] = __treeprune__ (T.Children, T.Parent, risk, held);
%! assert_equal (L, T.PruneList);
%! assert_equal (A, T.PruneAlpha, 1e-14);

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
