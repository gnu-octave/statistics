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

DEFUN_DLD (treepredict, args, nargout,
           "-*- texinfo -*-\n\
@deftypefn  {statistics} {@var{V} =} treepredict (@var{X}, @var{Children}, @\n\
@var{CutPredictorIndex}, @var{CutPoint}, @var{Value})\n\
@deftypefnx {statistics} {[@var{V}, @var{node}] =} treepredict (@dots{})\n\
\n\
\n\
Evaluate a binary decision tree on new data.\n\
\n\
@code{@var{V} = treepredict (@var{X}, @var{Children},\n\
@var{CutPredictorIndex}, @var{CutPoint}, @var{Value})} sends each row of the\n\
@math{NxP} matrix @var{X} down the tree and returns the value of the node it\n\
comes to rest at.  It is the prediction engine shared by\n\
@code{ClassificationTree} and @code{RegressionTree}, and is not meant to be\n\
called directly.\n\
\n\
@var{Children}, @var{CutPredictorIndex} and @var{CutPoint} are the node table\n\
returned by @code{treetrain}, a zero in @var{Children} marking a leaf.\n\
@var{Value} carries one row per node: the class probabilities of a\n\
classifier, in which case @var{V} is @math{NxK}, or the mean of a\n\
regression, in which case it is @math{Nx1}.\n\
\n\
@code{[@var{V}, @var{node}] = treepredict (@dots{})} also returns the index\n\
of the node each row came to rest at.\n\
\n\
A row is stopped by the first node whose split predictor it is missing, and\n\
takes that node's value, which is how such a row was held back rather than\n\
sent to a child while the tree was grown.\n\
\n\
@seealso{treetrain}\n\
@end deftypefn")
{
  if (args.length () != 5)
    error ("treepredict: invalid number of input arguments.");

  Matrix X = args(0).matrix_value ();
  Matrix children = args(1).matrix_value ();
  ColumnVector cutvar = args(2).column_vector_value ();
  ColumnVector cutpoint = args(3).column_vector_value ();
  Matrix value = args(4).matrix_value ();

  const octave_idx_type nn = children.rows ();
  if (nn < 1)
    error ("treepredict: the tree must hold at least one node.");
  if (children.columns () != 2)
    error ("treepredict: CHILDREN must have two columns.");
  if (cutvar.numel () != nn || cutpoint.numel () != nn
      || value.rows () != nn)
    error ("treepredict: the node table must be of one length throughout.");

  const octave_idx_type p = X.columns ();
  for (octave_idx_type i = 0; i < nn; i++)
    {
      if (children(i, 0) == 0)
        continue;
      octave_idx_type v = static_cast<octave_idx_type> (cutvar(i));
      if (v < 1 || v > p)
        error ("treepredict: a node cuts on a predictor X does not hold.");
      if (children(i, 0) < 1 || children(i, 0) > nn
          || children(i, 1) < 1 || children(i, 1) > nn)
        error ("treepredict: a node names a child outside the tree.");
    }

  Matrix V;
  ColumnVector node;
  tree_descend (X, children, cutvar, cutpoint, value, V, node);

  if (nargout > 1)
    return ovl (V, node);
  return ovl (V);
}

/*
%!test
%! ## Every row of the training data lands where it was counted
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'Prune', true);
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! prob = T.ClassWeight ./ sum (T.ClassWeight, 2);
%! [score, node] = treepredict (meas, T.Children, T.CutPredictorIndex, ...
%!                              T.CutPoint, prob);
%! assert_equal (size (score), [150, 3]);
%! assert_equal (node(1:12)', 2 * ones (1, 12));
%! assert_equal (sum (score(:,1)), 50, 1e-12);
%! leaf = T.Children(:,1) == 0;
%! assert_equal (all (leaf(node)), true);

%!test
%! ## A row is stopped by the first node whose predictor it is missing
%! load fisheriris
%! y = grp2idx (species);
%! o = struct ('NumClasses', 3, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 149, 'SplitCriterion', 'gdi', ...
%!             'MergeLeaves', true, 'Prune', true);
%! T = treetrain (meas, y, ones (150, 1) / 150, o);
%! prob = T.ClassWeight ./ sum (T.ClassWeight, 2);
%! q = meas([1, 51, 101, 71], :);
%! q(1,3) = NaN;
%! q(2,:) = NaN;
%! [score, node] = treepredict (q, T.Children, T.CutPredictorIndex, ...
%!                              T.CutPoint, prob);
%! assert_equal (node', [1, 1, 5, 5]);
%! assert_equal (score(1,:), prob(1,:));
%! assert_equal (score(3,:), [0, 1/46, 45/46], 1e-12);

%!test
%! ## The value of a regression leaf is its mean
%! g = mod ((1:20)', 3);
%! x = [(1:20)', g];
%! y = [ones(10,1); 5 * ones(10,1)];
%! o = struct ('NumClasses', 1, 'MinParent', 10, 'MinLeaf', 1, ...
%!             'MaxSplits', 19, 'SplitCriterion', 'mse', ...
%!             'MergeLeaves', true, 'Prune', true, 'QEToler', 1e-6);
%! T = treetrain (x, y, ones (20, 1) / 20, o);
%! [yfit, node] = treepredict (x, T.Children, T.CutPredictorIndex, ...
%!                             T.CutPoint, T.NodeMean);
%! assert_equal (size (yfit), [20, 1]);
%! assert_equal (yfit, y, 1e-12);
%! assert_equal (node', [2 * ones(1, 10), 3 * ones(1, 10)]);

%!test
%! ## A tree of one node answers with that node for every row
%! v = [0.25, 0.75];
%! [score, node] = treepredict (rand (7, 3), [0, 0], 0, NaN, v);
%! assert_equal (node, ones (7, 1));
%! assert_equal (score, repmat (v, 7, 1));

%!test
%! ## One output returns the value alone
%! v = [0.25, 0.75];
%! got = treepredict (rand (4, 2), [0, 0], 0, NaN, v);
%! assert_equal (got, repmat (v, 4, 1));

%!test
%! ## A value below the cut point goes left
%! kids = [2, 3; 0, 0; 0, 0];
%! score = treepredict ([1; 3], kids, [1; 0; 0], [2; NaN; NaN], [0; 10; 20]);
%! assert_equal (score, [10; 20]);

## Test input validation
%!error <treepredict: invalid number of input arguments.> ...
%! treepredict (1, 2, 3, 4);
%!error <treepredict: invalid number of input arguments.> ...
%! treepredict (1, 2, 3, 4, 5, 6);
%!error <treepredict: the tree must hold at least one node.> ...
%! treepredict (rand (3, 2), zeros (0, 2), [], [], zeros (0, 1));
%!error <treepredict: CHILDREN must have two columns.> ...
%! treepredict (rand (3, 2), [0, 0, 0], 0, NaN, 1);
%!error <treepredict: the node table must be of one length throughout.> ...
%! treepredict (rand (3, 2), [0, 0], [0; 0], NaN, 1);
%!error <treepredict: the node table must be of one length throughout.> ...
%! treepredict (rand (3, 2), [0, 0], 0, NaN, [1; 2]);
%!error <treepredict: a node cuts on a predictor X does not hold.> ...
%! kids = [2, 3; 0, 0; 0, 0]; ...
%! treepredict (rand (3, 2), kids, [5; 0; 0], [1; NaN; NaN], [0; 1; 2]);
%!error <treepredict: a node names a child outside the tree.> ...
%! kids = [2, 9; 0, 0; 0, 0]; ...
%! treepredict (rand (3, 2), kids, [1; 0; 0], [1; NaN; NaN], [0; 1; 2]);
*/
