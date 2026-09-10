## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{Mdl} =} fitrtree (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitrtree (@dots{}, @var{name}, @var{value})
##
## Fit a binary decision tree for regression.
##
## @code{@var{Mdl} = fitrtree (@var{X}, @var{Y})} grows a binary decision tree
## on the predictor data @var{X} and the response @var{Y}, and returns it as a
## @code{RegressionTree} object.
##
## @itemize
## @item
## @var{X} must be a @math{NxP} numeric matrix of predictor data, where rows
## correspond to observations and columns to predictors.
## @item
## @var{Y} must be a @math{Nx1} numeric vector holding the response of each
## observation in @var{X}.
## @end itemize
##
## An observation whose response is missing is dropped, and the rows kept are
## reported in @code{RowsUsed}.  An observation missing some of its predictors
## is kept: it descends the tree as far as the predictors it does carry allow
## and is answered there.
##
## @code{@var{Mdl} = fitrtree (@dots{}, @var{name}, @var{value})} takes the
## options below.
##
## @multitable @columnfractions 0.24 0.74
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'MaxNumSplits'} @tab A nonnegative integer, the largest number
## of branch nodes the tree may take.  The default is one less than the number
## of observations.
##
## @item @qcode{'MergeLeaves'} @tab @qcode{'on'} (default) or @qcode{'off'}.
## When on, a pair of leaves whose parent is no worse than the two of them
## together is merged back into that parent.
##
## @item @qcode{'MinLeafSize'} @tab A positive integer, the fewest
## observations a leaf may hold.  The default is 1.
##
## @item @qcode{'MinParentSize'} @tab A positive integer, the fewest
## observations a node must hold to be split.  The default is 10.  The value
## the fit uses is @code{max (MinParentSize, 2 * MinLeafSize)}.
##
## @item @qcode{'PredictorNames'} @tab A cell array of character vectors
## naming the columns of @var{X}.
##
## @item @qcode{'Prune'} @tab @qcode{'on'} (default) or @qcode{'off'}.  When
## on, the cost complexity pruning sequence is estimated and reported in
## @code{PruneList} and @code{PruneAlpha}.  The tree returned is the unpruned
## one either way; @code{prune} takes a subtree out of the sequence.
##
## @item @qcode{'PruneCriterion'} @tab @qcode{'mse'}, the only criterion a
## regression tree has.
##
## @item @qcode{'QuadraticErrorTolerance'} @tab A positive scalar.  A node
## whose squared error has fallen to this fraction of the root's is not split
## further.  The default is 1e-6.
##
## @item @qcode{'ResponseName'} @tab A character vector naming the response.
## The default is @qcode{'Y'}.
##
## @item @qcode{'ResponseTransform'} @tab A character vector naming a
## transform to apply to the predicted response, or a function handle.  The
## default is @qcode{'none'}.
##
## @item @qcode{'SplitCriterion'} @tab @qcode{'mse'}, the only criterion a
## regression tree has.
##
## @item @qcode{'Weights'} @tab A nonnegative numeric vector with one element
## per observation.  The default is uniform.
##
## @end multitable
##
## Categorical predictors, surrogate splits and predictor subsampling are not
## implemented, and an option asking for one of them is refused rather than
## quietly ignored.
##
## @seealso{RegressionTree, fitctree, treetrain, treepredict}
## @end deftypefn

function Mdl = fitrtree (X, Y, varargin)

  ## Input validation
  if (nargin < 2)
    error ("fitrtree: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitrtree: name-value arguments must be in pairs.");
  endif
  if (rows (X) != numel (Y))
    error ("fitrtree: number of rows in X and Y must be equal.");
  endif

  ## Parse arguments to the class constructor
  Mdl = RegressionTree (X, Y, varargin{:});

endfunction

## Demo
%!demo
%! ## Grow a regression tree on the carsmall data and look at it
%!
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG, 'MinLeafSize', 15);
%!
%! ## The tree as text: a branch names its cut, a leaf names what it fits
%! view (Mdl);
%!
%! ## How much each predictor contributed
%! predictorImportance (Mdl)
%!
%! ## The mean squared error on the data it was fitted to
%! resubLoss (Mdl)

%!demo
%! ## Prune a tree back and watch the error rise as it gets smaller
%!
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG, 'MinLeafSize', 15);
%!
%! levels = 0:numel (Mdl.PruneAlpha) - 1;
%! leaves = zeros (size (levels));
%! err = zeros (size (levels));
%! for ii = 1:numel (levels)
%!   sub = prune (Mdl, 'Level', levels(ii));
%!   leaves(ii) = sum (! sub.IsBranchNode);
%!   err(ii) = resubLoss (sub);
%! endfor
%! [leaves(:), err(:)]

## Tests
%!test  # MATLAB parity: the tree a default fit grows on carsmall
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG);
%! assert_equal (class (Mdl), 'RegressionTree');
%! assert_equal (Mdl.NumNodes, 37);
%! assert_equal (Mdl.NumObservations, 94);
%! assert_equal (sum (Mdl.RowsUsed), 94);
%! assert_equal (Mdl.NodeSize(1:9)', [94, 58, 36, 40, 17, 8, 28, 18, 22]);
%! assert_equal (Mdl.CutPredictorIndex(1:5)', [1, 3, 3, 1, 2]);
%! assert_equal (Mdl.CutPoint(1:5)', [3085.5, 89, 115, 2162, 5], 1e-12);

%!test  # MATLAB parity: the node statistics of the carsmall tree
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG);
%! assert_equal (Mdl.NodeMean(1:5)', [23.7181, 28.7931, 15.5417, 30.9375, ...
%!                                    24.0882], 1e-4);
%! assert_equal (Mdl.NodeError(1:5)', [63.8859, 30.4400, 9.4219, 24.9648, ...
%!                                     10.1834], 1e-4);
%! assert_equal (Mdl.NodeProbability(1:3)', [94, 58, 36] / 94, 1e-14);
%! assert_equal (Mdl.NodeRisk, Mdl.NodeProbability .* Mdl.NodeError, 1e-14);

%!test  # MATLAB parity: predict, its two outputs and the loss
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG);
%! [yFit, node] = predict (Mdl, X([1, 20, 60], :));
%! assert_equal (node', [20, 13, 36]);
%! assert_equal (yFit', [17.25, 12.3333333333333, 29.1], 1e-12);
%! assert_equal (resubLoss (Mdl), 5.5828069902791, 1e-12);
%! assert_equal (loss (Mdl, X, MPG), 5.5828069902791, 1e-12);

%!test  # MATLAB parity: the pruning sequence of the carsmall tree
%! ## One row of carsmall has no horsepower, and node 2 cuts on horsepower,
%! ## so that row stops there and belongs to neither child.  Its share of
%! ## the node's error is part of the subtree's risk, which is what puts
%! ## this alpha at 5.99 rather than 6.32.
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG);
%! assert_equal (Mdl.NodeSize(2) - sum (Mdl.NodeSize(Mdl.Children(2,:))), 1);
%! assert_equal (Mdl.PruneList(1:5)', [17, 16, 14, 15, 13]);
%! assert_equal (numel (Mdl.PruneAlpha), 18);
%! assert_equal (Mdl.PruneAlpha(17), 5.99325416896717, 1e-12);
%! assert_equal (Mdl.PruneAlpha(18), 41.4954735525515, 1e-11);

%!test  # MATLAB parity: predictor importance discounts what a node holds back
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG);
%! assert_equal (predictorImportance (Mdl), ...
%!               [2.5904, 0.1006, 0.5499], 1e-4);
%! Sub = fitrtree (X, MPG, 'MinLeafSize', 15);
%! assert_equal (predictorImportance (Sub), ...
%!               [11.253155105041, 0, 1.49831354224175], 1e-12);

%!test  # MATLAB parity: a smaller tree, its sequence and its text form
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG, 'MinLeafSize', 15);
%! assert_equal (Mdl.NumNodes, 9);
%! assert_equal (Mdl.NodeSize', [94, 58, 36, 40, 17, 15, 21, 18, 22]);
%! assert_equal (Mdl.PruneList', [4, 3, 1, 2, 0, 0, 0, 0, 0]);
%! assert_equal (Mdl.PruneAlpha', [0, 1.56476063829787, 1.95238622931442, ...
%!                                 5.99325416896717, 41.4954735525515], 1e-11);

%!test  # MATLAB parity: the growth options that stop the tree early
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! assert_equal (fitrtree (X, MPG, 'MinParentSize', 20).NumNodes, 15);
%! assert_equal (fitrtree (X, MPG, 'MinLeafSize', 15).NumNodes, 9);
%! assert_equal (fitrtree (X, MPG, 'MaxNumSplits', 3).NumNodes, 7);
%! assert_equal (fitrtree (X, MPG, 'MergeLeaves', 'off').NumNodes, 37);
%! assert_equal (fitrtree (X, MPG, ...
%!                         'QuadraticErrorTolerance', 0.01).NumNodes, 27);

%!test  # MATLAB parity: turning both reductions off leaves no sequence
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG, 'MergeLeaves', 'off', 'Prune', 'off');
%! assert_equal (isempty (Mdl.PruneList), true);
%! assert_equal (isempty (Mdl.PruneAlpha), true);

%!test  # MATLAB parity: the weights a weighted fit reports and weighs by
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG, 'Weights', (1:100)');
%! assert_equal (sum (Mdl.W), 1, 1e-14);
%! assert_equal (Mdl.W(1), 0.000201328769881216, 1e-15);
%! assert_equal (Mdl.NumNodes, 37);
%! assert_equal (Mdl.NodeSize(1:5)', [94, 56, 38, 27, 29]);
%! assert_equal (Mdl.NodeMean(1:3)', [26.469398026978, 30.3411058363586, ...
%!                                    16.4660894660895], 1e-12);
%! assert_equal (resubLoss (Mdl), 6.51484752783409, 1e-12);
%! assert_equal (loss (Mdl, Mdl.X, Mdl.Y), 7.44189886464653, 1e-12);

%!test  # MATLAB parity: the response transform reaches the prediction
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = fitrtree (X, MPG, 'ResponseTransform', @(y) 2 * y);
%! assert_equal (Mdl.ResponseTransform, '@(y) 2 * y');
%! assert_equal (predict (Mdl, X([1, 20, 60], :))', ...
%!               [34.5, 24.6666666666667, 58.2], 1e-12);
%! assert_equal (resubLoss (Mdl), 633.531443696713, 1e-10);

## Test input validation
%!error<fitrtree: too few arguments.> fitrtree ()
%!error<fitrtree: too few arguments.> fitrtree (ones (4, 1))
%!error<fitrtree: name-value arguments must be in pairs.>
%! fitrtree (ones (4, 2), ones (4, 1), 'K')
%!error<fitrtree: number of rows in X and Y must be equal.>
%! fitrtree (ones (4, 2), ones (3, 1))
%!error<fitrtree: number of rows in X and Y must be equal.>
%! fitrtree (ones (4, 2), ones (3, 1), 'K', 2)
