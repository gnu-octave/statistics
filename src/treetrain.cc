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
#include <octave/ov-struct.h>
#include "tree.cpp"

DEFUN_DLD (treetrain, args, ,
           "-*- texinfo -*-\n\
@deftypefn {statistics} {@var{Mdl} =} treetrain (@var{X}, @var{Y}, @var{W}, @\n\
@var{opts})\n\
\n\
\n\
Grow a binary decision tree by recursive partitioning.\n\
\n\
@code{@var{Mdl} = treetrain (@var{X}, @var{Y}, @var{W}, @var{opts})} grows a\n\
tree on the @math{NxP} predictor matrix @var{X} and returns it as a structure\n\
with one row per node.  It is the fitting engine shared by\n\
@code{ClassificationTree} and @code{RegressionTree}, and is not meant to be\n\
called directly.\n\
\n\
@var{Y} is an @math{Nx1} vector of class indices in @code{1:K} when fitting a\n\
classifier and the response itself when fitting a regression.  @var{W} is an\n\
@math{Nx1} vector of observation weights, already adjusted for any prior.\n\
\n\
@var{opts} is a structure carrying @qcode{NumClasses}, @qcode{MinParent},\n\
@qcode{MinLeaf}, @qcode{MaxSplits}, @qcode{SplitCriterion},\n\
@qcode{MergeLeaves}, @qcode{Prune} and, for a regression, @qcode{QEToler}.\n\
@qcode{SplitCriterion} is @qcode{\"gdi\"} or @qcode{\"deviance\"} for a\n\
classifier and @qcode{\"mse\"} for a regression, and it is what selects\n\
between the two.\n\
\n\
The returned structure holds @qcode{Children}, @qcode{Parent},\n\
@qcode{CutPredictorIndex}, @qcode{CutPoint}, @qcode{IsBranchNode},\n\
@qcode{NodeSize}, @qcode{NodeWeight}, @qcode{NumNodes}, @qcode{PruneList} and\n\
@qcode{PruneAlpha}, plus @qcode{ClassWeight} and @qcode{ClassCount} for a\n\
classifier or @qcode{NodeMean} and @qcode{NodeError} for a regression.\n\
Nodes are numbered as they are created, so a parent always precedes its\n\
children.\n\
\n\
@seealso{treepredict}\n\
@end deftypefn")
{
  if (args.length () != 4)
    error ("treetrain: invalid number of input arguments.");

  Matrix X = args(0).matrix_value ();
  ColumnVector yv = args(1).column_vector_value ();
  ColumnVector wv = args(2).column_vector_value ();
  octave_scalar_map opts = args(3).scalar_map_value ();

  if (yv.numel () != X.rows () || wv.numel () != X.rows ())
    error ("treetrain: X, Y and W must have the same number of rows.");

  TreeOpts o;
  o.K = opts.contents ("NumClasses").idx_type_value ();
  o.minparent = opts.contents ("MinParent").idx_type_value ();
  o.minleaf = opts.contents ("MinLeaf").idx_type_value ();
  o.maxsplits = opts.contents ("MaxSplits").idx_type_value ();
  o.mergeleaves = opts.contents ("MergeLeaves").bool_value ();
  o.prune = opts.isfield ("Prune")
            ? opts.contents ("Prune").bool_value () : false;
  o.qetoler = opts.isfield ("QEToler")
              ? opts.contents ("QEToler").double_value () : 0.0;

  std::string critname = opts.contents ("SplitCriterion").string_value ();
  if (critname == "gdi")
    o.crit = GDI;
  else if (critname == "deviance")
    o.crit = DEVIANCE;
  else if (critname == "mse")
    o.crit = MSE;
  else
    error ("treetrain: unsupported SplitCriterion.");

  if (o.minleaf < 1)
    error ("treetrain: MinLeaf must be a positive integer.");
  if (o.minparent < 1)
    error ("treetrain: MinParent must be a positive integer.");
  if (o.crit != MSE && o.K < 1)
    error ("treetrain: NumClasses must be a positive integer.");
  for (octave_idx_type i = 0; i < wv.numel (); i++)
    if (! (wv(i) >= 0.0))
      error ("treetrain: W must hold non-negative weights.");

  return ovl (tree_build (X, yv, wv, o));
}

/*
%!test
%! ## The iris tree, grown and then merged, as MATLAB grows it
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'Prune', false);
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! assert_equal (T.NumNodes, 11);
%! assert_equal (T.NodeSize', [150, 50, 100, 54, 46, 48, 6, 3, 43, 47, 1]);
%! assert_equal (T.CutPredictorIndex', [3, 0, 4, 3, 3, 4, 0, 0, 0, 0, 0]);
%! assert_equal (T.CutPoint([1, 3, 4, 5, 6])', [2.45, 1.75, 4.95, 4.85, 1.65]);

%!test
%! ## MergeLeaves collapses the one pair that lowers no misclassification
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'Prune', false);
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! assert_equal (T.NumNodes, 9);
%! assert_equal (T.NodeSize', [150, 50, 100, 54, 46, 48, 6, 47, 1]);
%! assert_equal (T.Parent', [0, 1, 1, 3, 3, 4, 4, 6, 6]);

%!test
%! ## Cost complexity pruning, levels and alphas
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'Prune', true);
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! assert_equal (T.PruneList', [4, 0, 3, 2, 0, 1, 0, 0, 0]);
%! assert_equal (T.PruneAlpha', [0, 1/150, 2/150, 44/150, 50/150], 1e-12);

%!test
%! ## A deeper tree, every branch pruned at its own level
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 2, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'Prune', true);
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! assert_equal (T.NumNodes, 17);
%! assert_equal (T.NodeSize', [150, 50, 100, 54, 46, 48, 6, 3, 43, 47, 1, ...
%!                             3, 3, 1, 2, 2, 1]);
%! assert_equal (T.PruneList', [5, 0, 4, 3, 1, 2, 2, 1, 0, 0, 0, 0, 2, ...
%!                              0, 0, 0, 0]);

%!test
%! ## A subtree that costs nothing to give up opens no level of the sequence
%! ## The pair of leaves that merging would have removed survives here, and
%! ## giving it up costs nothing, so the eleven node tree carries the merged
%! ## tree's five alphas rather than six.  Measured on R2024a.
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'Prune', true);
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! assert_equal (T.NumNodes, 11);
%! assert_equal (T.PruneList', [4, 0, 3, 2, 0, 1, 0, 0, 0, 0, 0]);
%! assert_equal (T.PruneAlpha', [0, 1/150, 2/150, 44/150, 50/150], 1e-12);

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
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'Prune', true);
%! T = treetrain (x, y, ones (150, 1) / 150, o);
%! assert_equal (T.NumNodes, 9);
%! assert_equal (T.NodeSize', [150, 50, 100, 45, 55, 25, 1, 8, 46]);
%! assert_equal (T.NodeSize(4) - T.NodeSize(6) - T.NodeSize(7), 19);
%! assert_equal (T.PruneList', [4, 0, 3, 1, 2, 0, 0, 0, 0]);
%! assert_equal (T.PruneAlpha', [0, 0.00385185185185185, ...
%!                               0.00593939393939394, 0.286666666666666, ...
%!                               0.333333333333333], 1e-14);

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
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'Prune', true);
%! T = treetrain (x, y, ones (150, 1) / 150, o);
%! assert_equal (T.NumNodes, 11);
%! assert_equal (T.PruneList', [4, 0, 3, 1, 2, 0, 0, 0, 0, 0, 0]);
%! assert_equal (T.PruneAlpha', [0, 0.00385185185185185, ...
%!                               0.00593939393939394, 0.286666666666666, ...
%!                               0.333333333333333], 1e-14);

%!test
%! ## Prune off leaves the pruning sequence empty
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'Prune', false);
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! assert_equal (isempty (T.PruneList), true);
%! assert_equal (isempty (T.PruneAlpha), true);

%!test
%! ## A row missing the split predictor descends to neither child
%! load fisheriris
%! y = grp2idx (species);
%! x = meas;
%! x(51:60, 4) = NaN;
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'Prune', false);
%! T = treetrain (x, y, ones (150, 1) / 150, o);
%! assert_equal (T.NodeSize', [150, 50, 100, 45, 55, 35, 1, 8, 46, 3, 43]);
%! assert_equal (T.NodeSize(4) - T.NodeSize(6) - T.NodeSize(7), 9);

%!test
%! ## A row missing every predictor is dropped outright
%! load fisheriris
%! y = grp2idx (species);
%! x = meas;
%! x(1:5, :) = NaN;
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'Prune', false);
%! T = treetrain (x, y, ones (150, 1) / 150, o);
%! assert_equal (T.NodeSize(1), 145);
%! assert_equal (T.NodeSize', [145, 45, 100, 54, 46, 48, 6, 3, 43, 47, 1]);

%!test
%! ## A predictor that is entirely missing is simply never chosen
%! load fisheriris
%! y = grp2idx (species);
%! x = meas;
%! x(:, 2) = NaN;
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'Prune', false);
%! T = treetrain (x, y, ones (150, 1) / 150, o);
%! assert_equal (T.CutPredictorIndex', [3, 0, 4, 3, 3, 4, 0, 0, 0, 0, 0]);

%!test
%! ## The split is judged on the rows a predictor has, scaled by their share:
%! ## x1 is known for every row and wins, x2 is known for ten and does not
%! y = [1;1;1;1;1;1;1;1;1;2;1;2;2;2;2;2;2;2;2;2];
%! x2 = NaN (20, 1);
%! x2(1:5) = (1:5)';
%! x2(12:16) = (6:10)';
%! o = struct ('NumClasses', 2, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 19, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'Prune', false);
%! T = treetrain ([(1:20)', x2], y, ones (20, 1) / 20, o);
%! assert_equal (T.CutPredictorIndex(1), 1);
%! assert_equal (T.CutPoint(1), 9.5);
%! assert_equal (T.NodeSize', [20, 9, 11, 2, 9]);

%!test
%! ## MinParent holds a node together and MaxSplits stops the tree
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 60, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'Prune', false);
%! assert_equal (treetrain (meas, y, ones (150, 1) / 150, o).NumNodes, 5);
%! o.MinParent = 10;
%! o.MaxSplits = 2;
%! assert_equal (treetrain (meas, y, ones (150, 1) / 150, o).NumNodes, 5);

%!test
%! ## ClassCount counts rows and ClassWeight weighs them, which are the same
%! ## thing only while the weights are
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'Prune', true);
%! count = [50, 50, 0, 0, 0, 0, 0, 0, 0; ...
%!          50, 0, 50, 49, 1, 47, 2, 47, 0; ...
%!          50, 0, 50, 5, 45, 1, 4, 0, 1];
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! assert_equal (T.ClassCount', count);
%! assert_equal (sum (T.ClassCount, 2), T.NodeSize);
%! w = [3 * ones(50, 1); ones(100, 1)];
%! w = w / sum (w);
%! W = treetrain (meas, y, w, o);
%! assert_equal (W.NodeSize', [150, 50, 100, 54, 46, 48, 6, 47, 1]);
%! assert_equal (W.ClassCount', count);
%! prob = W.ClassWeight ./ sum (W.ClassWeight, 2);
%! assert_equal (prob(1,:), [0.6, 0.2, 0.2], 1e-12);

%!test
%! ## A regression carries a leaf value in place of class counts
%! g = zeros (20, 1);
%! x = [(1:20)', g];
%! y = [ones(10,1); 5 * ones(10,1)];
%! o = struct ('NumClasses', 1, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 19, 'SplitCriterion', 'mse', ...
%!             'MergeLeaves', true, 'Prune', true, 'QEToler', 1e-6);
%! T = treetrain (x, y, ones (20, 1) / 20, o);
%! assert_equal (isempty (T.ClassCount), true);
%! assert_equal (isempty (T.ClassWeight), true);

%!test
%! ## A regression tree, its leaf values and its pruning sequence
%! g = mod ((1:20)', 3);
%! x = [(1:20)', g];
%! y = [ones(10,1); 5 * ones(10,1)];
%! o = struct ('NumClasses', 1, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 19, 'SplitCriterion', 'mse', ...
%!             'MergeLeaves', true, 'Prune', true, 'QEToler', 1e-6);
%! T = treetrain (x, y, ones (20, 1) / 20, o);
%! assert_equal (T.NumNodes, 3);
%! assert_equal (T.NodeSize', [20, 10, 10]);
%! assert_equal (T.NodeMean', [3, 1, 5], 1e-12);
%! assert_equal (T.NodeError', [4, 0, 0], 1e-12);
%! assert_equal (T.PruneList', [1, 0, 0]);

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
%!             'MergeLeaves', true, 'Prune', true, 'QEToler', 1e-6);
%! T = treetrain (X(ok, :), MPG(ok), ones (94, 1) / 94, o);
%! assert_equal (T.NumNodes, 37);
%! assert_equal (T.NodeSize(2) - sum (T.NodeSize(T.Children(2,:))), 1);
%! assert_equal (T.PruneList(1:5)', [17, 16, 14, 15, 13]);
%! assert_equal (numel (T.PruneAlpha), 18);
%! assert_equal (T.PruneAlpha(17), 5.99325416896717, 1e-12);
%! assert_equal (T.PruneAlpha(18), 41.4954735525515, 1e-11);

%!test
%! ## QEToler stops a node whose squared error is a small share of the root's
%! g = mod ((1:40)', 7);
%! x = [(1:40)', g];
%! y = (1:40)';
%! o = struct ('NumClasses', 1, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 39, 'SplitCriterion', 'mse', ...
%!             'MergeLeaves', true, 'Prune', false, 'QEToler', 1e-6);
%! loose = o;
%! loose.QEToler = 0.05;
%! assert_equal (treetrain (x, y, ones (40, 1) / 40, loose).NumNodes ...
%!               < treetrain (x, y, ones (40, 1) / 40, o).NumNodes, true);

%!test
%! ## A single class, a single row and identical rows are each one leaf
%! o = struct ('NumClasses', 2, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 19, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'Prune', true);
%! w = ones (20, 1) / 20;
%! assert_equal (treetrain (rand (20, 2), ones (20, 1), w, o).NumNodes, 1);
%! assert_equal (treetrain ([1, 2], 1, 1, o).NumNodes, 1);
%! y = [ones(10,1); 2*ones(10,1)];
%! assert_equal (treetrain (ones (20, 2), y, w, o).NumNodes, 1);

## Test input validation
%!shared o
%! o = struct ('NumClasses', 2, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 19, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'Prune', true);
%!error <treetrain: invalid number of input arguments.> treetrain (1, 2, 3);
%!error <treetrain: invalid number of input arguments.> ...
%! treetrain (1, 2, 3, o, 5);
%!error <treetrain: X, Y and W must have the same number of rows.> ...
%! treetrain (rand (10, 2), ones (9, 1), ones (10, 1), o);
%!error <treetrain: X, Y and W must have the same number of rows.> ...
%! treetrain (rand (10, 2), ones (10, 1), ones (9, 1), o);
%!error <treetrain: unsupported SplitCriterion.> ...
%! bad = setfield (o, 'SplitCriterion', 'twoing'); ...
%! treetrain (rand (10, 2), ones (10, 1), ones (10, 1), bad);
%!error <treetrain: MinLeaf must be a positive integer.> ...
%! bad = setfield (o, 'MinLeaf', 0); ...
%! treetrain (rand (10, 2), ones (10, 1), ones (10, 1), bad);
%!error <treetrain: MinParent must be a positive integer.> ...
%! bad = setfield (o, 'MinParent', 0); ...
%! treetrain (rand (10, 2), ones (10, 1), ones (10, 1), bad);
%!error <treetrain: NumClasses must be a positive integer.> ...
%! bad = setfield (o, 'NumClasses', 0); ...
%! treetrain (rand (10, 2), ones (10, 1), ones (10, 1), bad);
%!error <treetrain: W must hold non-negative weights.> ...
%! treetrain (rand (10, 2), ones (10, 1), [-1; ones(9, 1)], o);
%!error <treetrain: Y holds a class index outside 1:K.> ...
%! treetrain (rand (10, 2), 5 * ones (10, 1), ones (10, 1), o);
%!error <treetrain: Y holds a class index outside 1:K.> ...
%! treetrain (rand (10, 2), zeros (10, 1), ones (10, 1), o);
*/
