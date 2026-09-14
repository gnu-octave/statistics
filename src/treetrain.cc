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
@qcode{MergeLeaves} and, for a regression, @qcode{QEToler}.\n\
@qcode{SplitCriterion} is @qcode{\"gdi\"} or @qcode{\"deviance\"} for a\n\
classifier and @qcode{\"mse\"} for a regression, and it is what selects\n\
between the two.\n\
\n\
@qcode{NumVariablesToSample}, when present, is the number of predictors\n\
each node chooses its split from, drawn afresh at every node, and\n\
@qcode{Seed}, an integer from 0 to @math{2^32-1}, seeds that draw, so the\n\
same seed grows the same tree.  A node whose predictors hold no valid split\n\
is a leaf.  Without @qcode{NumVariablesToSample} every predictor is tried and\n\
nothing is drawn.\n\
\n\
@qcode{CategoricalPredictors}, when present, lists the columns of @var{X}\n\
that hold levels rather than numbers, and such a column is split into two\n\
sets of levels.  A regression orders the levels by mean response and two\n\
classes by the probability of the first class.  More classes search every\n\
partition when the node holds at most @qcode{MaxNumCategories} levels, 10\n\
by default, and otherwise take MATLAB's heuristics, as\n\
@qcode{AlgorithmForCategorical} selects: @qcode{\"auto\"} (default),\n\
@qcode{\"exact\"}, @qcode{\"pullleft\"}, @qcode{\"pca\"} or\n\
@qcode{\"ovabyclass\"}.\n\
\n\
The returned structure holds @qcode{Children}, @qcode{Parent},\n\
@qcode{CutPredictorIndex}, @qcode{CutPoint}, @qcode{CutCategories},\n\
@qcode{IsBranchNode}, @qcode{NodeSize}, @qcode{NodeWeight} and\n\
@qcode{NumNodes}, plus\n\
@qcode{ClassWeight} and @qcode{ClassCount} for a classifier or\n\
@qcode{NodeMean} and @qcode{NodeError} for a regression.\n\
Nodes are numbered as they are created, so a parent always precedes its\n\
children.\n\
\n\
The cost complexity pruning sequence is not built here.  A tree outlives\n\
the data it was grown on, so the sequence is @code{__treeprune__}, which\n\
takes a node table and the risk the caller measures by.\n\
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
  o.qetoler = opts.isfield ("QEToler")
              ? opts.contents ("QEToler").double_value () : 0.0;

  o.nvars = 0;
  if (opts.isfield ("NumVariablesToSample"))
    {
      double v = opts.contents ("NumVariablesToSample").double_value ();
      if (! (v >= 1.0 && v == std::floor (v)))
        error ("treetrain: NumVariablesToSample must be a positive integer.");
      if (v > X.columns ())
        error ("treetrain: NumVariablesToSample must not exceed the number "
               "of predictors.");
      o.nvars = static_cast<octave_idx_type> (v);
    }
  o.seed = 0;
  if (opts.isfield ("Seed"))
    {
      double v = opts.contents ("Seed").double_value ();
      if (! (v >= 0.0 && v <= 4294967295.0 && v == std::floor (v)))
        error ("treetrain: Seed must be an integer from 0 to 2^32-1.");
      o.seed = static_cast<std::uint32_t> (v);
    }

  o.iscat.clear ();
  if (opts.isfield ("CategoricalPredictors"))
    {
      const NDArray c = opts.contents ("CategoricalPredictors").array_value ();
      if (c.numel () > 0)
        o.iscat.assign (X.columns (), false);
      for (octave_idx_type i = 0; i < c.numel (); i++)
        {
          const double v = c(i);
          if (! (v >= 1.0 && v <= X.columns () && v == std::floor (v)))
            error ("treetrain: CategoricalPredictors must hold indices of "
                   "columns of X.");
          o.iscat[static_cast<octave_idx_type> (v) - 1] = true;
        }
    }
  o.maxcat = 10.0;
  if (opts.isfield ("MaxNumCategories"))
    {
      const double v = opts.contents ("MaxNumCategories").double_value ();
      if (! (v >= 0.0 && (std::isinf (v) || v == std::floor (v))))
        error ("treetrain: MaxNumCategories must be a non-negative integer.");
      o.maxcat = v;
    }
  o.catalg = CAT_AUTO;
  if (opts.isfield ("AlgorithmForCategorical"))
    {
      const std::string a
        = opts.contents ("AlgorithmForCategorical").string_value ();
      if (a == "auto")
        o.catalg = CAT_AUTO;
      else if (a == "exact")
        o.catalg = CAT_EXACT;
      else if (a == "pullleft")
        o.catalg = CAT_PULLLEFT;
      else if (a == "pca")
        o.catalg = CAT_PCA;
      else if (a == "ovabyclass")
        o.catalg = CAT_OVA;
      else
        error ("treetrain: unsupported AlgorithmForCategorical.");
    }

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
%!             'MergeLeaves', false);
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
%!             'MergeLeaves', true);
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! assert_equal (T.NumNodes, 9);
%! assert_equal (T.NodeSize', [150, 50, 100, 54, 46, 48, 6, 47, 1]);
%! assert_equal (T.Parent', [0, 1, 1, 3, 3, 4, 4, 6, 6]);

%!test
%! ## A deeper tree, grown until no node can be split further
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 2, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true);
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! assert_equal (T.NumNodes, 17);
%! assert_equal (T.NodeSize', [150, 50, 100, 54, 46, 48, 6, 3, 43, 47, 1, ...
%!                             3, 3, 1, 2, 2, 1]);

%!test
%! ## A row missing the split predictor descends to neither child
%! load fisheriris
%! y = grp2idx (species);
%! x = meas;
%! x(51:60, 4) = NaN;
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false);
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
%!             'MergeLeaves', false);
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
%!             'MergeLeaves', false);
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
%!             'MergeLeaves', false);
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
%!             'MergeLeaves', false);
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
%!             'MergeLeaves', true);
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
%!             'MergeLeaves', true, 'QEToler', 1e-6);
%! T = treetrain (x, y, ones (20, 1) / 20, o);
%! assert_equal (isempty (T.ClassCount), true);
%! assert_equal (isempty (T.ClassWeight), true);

%!test
%! ## A regression tree and its leaf values
%! g = mod ((1:20)', 3);
%! x = [(1:20)', g];
%! y = [ones(10,1); 5 * ones(10,1)];
%! o = struct ('NumClasses', 1, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 19, 'SplitCriterion', 'mse', ...
%!             'MergeLeaves', true, 'QEToler', 1e-6);
%! T = treetrain (x, y, ones (20, 1) / 20, o);
%! assert_equal (T.NumNodes, 3);
%! assert_equal (T.NodeSize', [20, 10, 10]);
%! assert_equal (T.NodeMean', [3, 1, 5], 1e-12);
%! assert_equal (T.NodeError', [4, 0, 0], 1e-12);

%!test
%! ## QEToler stops a node whose squared error is a small share of the root's
%! g = mod ((1:40)', 7);
%! x = [(1:40)', g];
%! y = (1:40)';
%! o = struct ('NumClasses', 1, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 39, 'SplitCriterion', 'mse', ...
%!             'MergeLeaves', true, 'QEToler', 1e-6);
%! loose = o;
%! loose.QEToler = 0.05;
%! assert_equal (treetrain (x, y, ones (40, 1) / 40, loose).NumNodes ...
%!               < treetrain (x, y, ones (40, 1) / 40, o).NumNodes, true);

%!test
%! ## A single class, a single row and identical rows are each one leaf
%! o = struct ('NumClasses', 2, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 19, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true);
%! w = ones (20, 1) / 20;
%! assert_equal (treetrain (rand (20, 2), ones (20, 1), w, o).NumNodes, 1);
%! assert_equal (treetrain ([1, 2], 1, 1, o).NumNodes, 1);
%! y = [ones(10,1); 2*ones(10,1)];
%! assert_equal (treetrain (ones (20, 2), y, w, o).NumNodes, 1);

## Test input validation
%!test
%! ## Sampling every predictor draws nothing and grows the unsampled tree
%! load fisheriris
%! o = struct ('NumClasses', 3, 'MinParent', 2, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false);
%! w = ones (150, 1) / 150;
%! T = treetrain (meas, grp2idx (species), w, o);
%! o.NumVariablesToSample = 4;
%! o.Seed = 12345;
%! assert_equal (isequaln (treetrain (meas, grp2idx (species), w, o), T), true);

%!test
%! ## The same seed grows the same tree
%! load fisheriris
%! o = struct ('NumClasses', 3, 'MinParent', 2, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'NumVariablesToSample', 1, 'Seed', 7);
%! w = ones (150, 1) / 150;
%! A = treetrain (meas, grp2idx (species), w, o);
%! B = treetrain (meas, grp2idx (species), w, o);
%! assert_equal (isequaln (A, B), true);

%!test
%! ## Different seeds grow different trees
%! load fisheriris
%! o = struct ('NumClasses', 3, 'MinParent', 2, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'NumVariablesToSample', 1, 'Seed', 7);
%! w = ones (150, 1) / 150;
%! assert_equal (treetrain (meas, grp2idx (species), w, o).NumNodes, 23);
%! o.Seed = 8;
%! assert_equal (treetrain (meas, grp2idx (species), w, o).NumNodes, 25);

%!test
%! ## A node whose sampled predictor holds no split is a leaf
%! load fisheriris
%! o = struct ('NumClasses', 3, 'MinParent', 2, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'NumVariablesToSample', 1, 'Seed', 0);
%! X = [ones(150, 1), meas(:, 3)];
%! T = treetrain (X, grp2idx (species), ones (150, 1) / 150, o);
%! assert_equal (T.NumNodes, 1);

%!test
%! ## A predictor that cannot split is never cut, however the draw falls
%! load fisheriris
%! o = struct ('NumClasses', 3, 'MinParent', 2, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'NumVariablesToSample', 1, 'Seed', 14);
%! X = [ones(150, 1), meas(:, 3)];
%! T = treetrain (X, grp2idx (species), ones (150, 1) / 150, o);
%! assert_equal (T.NumNodes, 15);
%! cut = T.CutPredictorIndex(logical (T.IsBranchNode));
%! assert_equal (all (cut == 2), true);

%!test
%! ## A sampled regression tree is reproduced by its seed
%! load fisheriris
%! o = struct ('NumClasses', 0, 'MinParent', 10, 'MinLeaf', 5, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'mse', ...
%!             'MergeLeaves', false, 'NumVariablesToSample', 1, 'Seed', 3);
%! w = ones (150, 1) / 150;
%! A = treetrain (meas(:, 2:4), meas(:, 1), w, o);
%! B = treetrain (meas(:, 2:4), meas(:, 1), w, o);
%! assert_equal (isequaln (A, B), true);

%!shared o
%! o = struct ('NumClasses', 2, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 19, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true);
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
%!error <treetrain: NumVariablesToSample must be a positive integer.> ...
%! bad = setfield (o, 'NumVariablesToSample', 0); ...
%! treetrain (rand (10, 2), ones (10, 1), ones (10, 1), bad);
%!error <treetrain: NumVariablesToSample must be a positive integer.> ...
%! bad = setfield (o, 'NumVariablesToSample', 1.5); ...
%! treetrain (rand (10, 2), ones (10, 1), ones (10, 1), bad);
%!error <treetrain: NumVariablesToSample must not exceed the number of predictors.> ...
%! bad = setfield (o, 'NumVariablesToSample', 3); ...
%! treetrain (rand (10, 2), ones (10, 1), ones (10, 1), bad);
%!error <treetrain: Seed must be an integer from 0 to 2\^32-1.> ...
%! bad = setfield (o, 'Seed', -1); ...
%! treetrain (rand (10, 2), ones (10, 1), ones (10, 1), bad);
%!error <treetrain: Seed must be an integer from 0 to 2\^32-1.> ...
%! bad = setfield (o, 'Seed', 2.5); ...
%! treetrain (rand (10, 2), ones (10, 1), ones (10, 1), bad);
%!error <treetrain: Seed must be an integer from 0 to 2\^32-1.> ...
%! bad = setfield (o, 'Seed', 2^32); ...
%! treetrain (rand (10, 2), ones (10, 1), ones (10, 1), bad);
%!test
%! ## MATLAB parity: two classes order the levels and split them in two
%! k = (0:79)';
%! c = mod (k, 4) + 1;
%! j = floor (k / 4);
%! y = (c == 1) | (c == 2 & mod (j, 4) != 0) | (c == 3 & mod (j, 4) == 0);
%! o = struct ('NumClasses', 2, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 79, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'CategoricalPredictors', 1);
%! T = treetrain ([c, mod(k * 7, 10)], double (y) + 1, ones (80, 1) / 80, o);
%! assert_equal (T.NumNodes, 3);
%! assert_equal (T.CutCategories(1,:), {[1, 2], [3, 4]});
%! assert_equal (isnan (T.CutPoint(1)), true);

%!test
%! ## MATLAB parity: a regression orders the levels by their mean response
%! k = (0:99)';
%! c = mod (k, 5) + 1;
%! means = [3, 1, 4, 1.5, 5];
%! y = means(c)' + 0.1 * sin (k);
%! o = struct ('NumClasses', 1, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 99, 'SplitCriterion', 'mse', ...
%!             'MergeLeaves', true, 'QEToler', 1e-6, ...
%!             'CategoricalPredictors', 1);
%! T = treetrain ([c, mod(k * 3, 8)], y, ones (100, 1) / 100, o);
%! assert_equal (T.NumNodes, 39);
%! assert_equal (T.CutCategories(1,:), {[2, 4], [1, 3, 5]});

%!test
%! ## Every named algorithm finds a split of three classes
%! k = (0:119)';
%! c = mod (k, 6) + 1;
%! y = mod (c + floor (k / 12), 3) + 1;
%! for a = {'exact', 'pullleft', 'pca', 'ovabyclass'}
%!   o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!               'MaxSplits', 1, 'SplitCriterion', 'gdi', ...
%!               'MergeLeaves', false, 'CategoricalPredictors', 1, ...
%!               'AlgorithmForCategorical', a{1});
%!   T = treetrain (c, y, ones (120, 1) / 120, o);
%!   assert_equal (sort ([T.CutCategories{1,:}]), 1:6);
%! endfor

%!error <treetrain: CategoricalPredictors must hold indices of columns of X.> ...
%! o = struct ('NumClasses', 2, 'MinParent', 1, 'MinLeaf', 1, ...
%!             'MaxSplits', 1, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'CategoricalPredictors', 3); ...
%! treetrain (ones (4, 2), [1; 1; 2; 2], ones (4, 1), o);
%!error <treetrain: MaxNumCategories must be a non-negative integer.> ...
%! o = struct ('NumClasses', 2, 'MinParent', 1, 'MinLeaf', 1, ...
%!             'MaxSplits', 1, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'MaxNumCategories', -1); ...
%! treetrain (ones (4, 2), [1; 1; 2; 2], ones (4, 1), o);
%!error <treetrain: unsupported AlgorithmForCategorical.> ...
%! o = struct ('NumClasses', 2, 'MinParent', 1, 'MinLeaf', 1, ...
%!             'MaxSplits', 1, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', false, 'AlgorithmForCategorical', 'bogus'); ...
%! treetrain (ones (4, 2), [1; 1; 2; 2], ones (4, 1), o);
*/
