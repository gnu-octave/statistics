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
## @deftypefn  {statistics} {@var{Mdl} =} fitctree (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitctree (@dots{}, @var{name}, @var{value})
##
## Fit a binary decision tree for classification.
##
## @code{@var{Mdl} = fitctree (@var{X}, @var{Y})} grows a binary decision tree
## on the predictor data @var{X} and the class labels @var{Y}, and returns it
## as a @code{ClassificationTree} object.
##
## @itemize
## @item
## @var{X} must be a @math{NxP} numeric matrix of predictor data, where rows
## correspond to observations and columns to predictors.
## @item
## @var{Y} must be a @math{Nx1} numeric or logical vector, a character array
## with one class name per row, or a cell array of character vectors, holding
## the class label of each observation in @var{X}.  The class names come back
## in the type @var{Y} was given in.
## @end itemize
##
## An observation whose class label is missing is dropped, and the rows kept
## are reported in @code{RowsUsed}.  An observation missing some of its
## predictors is kept: it descends the tree as far as the predictors it does
## carry allow and is answered there.
##
## @code{@var{Mdl} = fitctree (@dots{}, @var{name}, @var{value})} takes the
## options below.
##
## @multitable @columnfractions 0.20 0.78
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'ClassNames'} @tab The classes to fit, of the same type as
## @var{Y}.  Observations of any other class are dropped.
##
## @item @qcode{'Cost'} @tab A square matrix with one row and column per
## class, where element @math{(i,j)} is the cost of classifying an
## observation of class @math{i} into class @math{j}, or a structure with
## fields @qcode{ClassNames} and @qcode{ClassificationCosts}.  The default is
## @code{1 - eye (K)}.  A non-default cost changes the shape of the tree, not
## only what it predicts.
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
## @item @qcode{'Prior'} @tab @qcode{'empirical'} (default),
## @qcode{'uniform'}, a numeric vector with one element per class, or a
## structure with fields @qcode{ClassNames} and @qcode{ClassProbs}.
##
## @item @qcode{'Prune'} @tab @qcode{'on'} (default) or @qcode{'off'}.  When
## on, the cost complexity pruning sequence is estimated and reported in
## @code{PruneList} and @code{PruneAlpha}.  The tree returned is the unpruned
## one either way; @code{prune} takes a subtree out of the sequence.
##
## @item @qcode{'PruneCriterion'} @tab @qcode{'error'}, the only criterion
## implemented.
##
## @item @qcode{'ResponseName'} @tab A character vector naming the response.
## The default is @qcode{'Y'}.
##
## @item @qcode{'ScoreTransform'} @tab A character vector naming a transform
## to apply to the scores, or a function handle.  The default is
## @qcode{'none'}.
##
## @item @qcode{'SplitCriterion'} @tab @qcode{'gdi'} (default), the Gini
## diversity index, or @qcode{'deviance'}, the cross entropy.
##
## @item @qcode{'Weights'} @tab A nonnegative numeric vector with one element
## per observation.  The default is uniform.
##
## @end multitable
##
## Categorical predictors, surrogate splits, predictor subsampling and the
## @qcode{'twoing'} split criterion are not implemented, and an option asking
## for one of them is refused rather than quietly ignored.
##
## @seealso{ClassificationTree, treetrain, treepredict}
## @end deftypefn

function Mdl = fitctree (X, Y, varargin)

  ## Input validation
  if (nargin < 2)
    error ("fitctree: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitctree: name-value arguments must be in pairs.");
  endif
  if (rows (X) != rows (Y))
    error ("fitctree: number of rows in X and Y must be equal.");
  endif

  ## Parse arguments to the class constructor
  Mdl = ClassificationTree (X, Y, varargin{:});

endfunction

## Demo
%!demo
%! ## Grow a classification tree on Fisher's iris data and look at it
%!
%! load fisheriris
%! Mdl = fitctree (meas, species);
%!
%! ## The tree as text: a branch names its cut, a leaf names its class
%! view (Mdl);
%!
%! ## How much each predictor contributed
%! predictorImportance (Mdl)
%!
%! ## The resubstitution error of the whole tree and of its subtrees
%! resubLoss (Mdl)

%!demo
%! ## Prune a tree back and watch the error rise as it gets smaller
%!
%! load fisheriris
%! Mdl = fitctree (meas, species);
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
%!test  # MATLAB parity: the tree a default fit grows on fisheriris
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! assert_equal (class (Mdl), 'ClassificationTree');
%! assert_equal (Mdl.NumNodes, 9);
%! assert_equal (Mdl.NodeSize', [150, 50, 100, 54, 46, 48, 6, 47, 1]);
%! assert_equal (Mdl.Parent', [0, 1, 1, 3, 3, 4, 4, 6, 6]);
%! assert_equal (Mdl.CutPredictorIndex', [3, 0, 4, 3, 0, 4, 0, 0, 0]);

%!test  # MATLAB parity: the cut points and the predictors they name
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! assert_equal (Mdl.CutPoint([1, 3, 4, 6])', [2.45, 1.75, 4.95, 1.65]);
%! assert_equal (isnan (Mdl.CutPoint([2, 5, 7, 8, 9]))', true (1, 5));
%! assert_equal (Mdl.CutPredictor', {'x3', '', 'x4', 'x3', '', 'x4', ...
%!                                   '', '', ''});
%! assert_equal (Mdl.CutType', {'continuous', '', 'continuous', ...
%!                              'continuous', '', 'continuous', '', '', ''});

%!test  # MATLAB parity: the class the fit reports for every node
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! assert_equal (Mdl.NodeClass', {'setosa', 'setosa', 'versicolor', ...
%!                                'versicolor', 'virginica', 'versicolor', ...
%!                                'virginica', 'versicolor', 'virginica'});
%! assert_equal (Mdl.ClassCount, [50, 50, 50; 50, 0, 0; 0, 50, 50; ...
%!                                0, 49, 5; 0, 1, 45; 0, 47, 1; ...
%!                                0, 2, 4; 0, 47, 0; 0, 0, 1]);

%!test  # MATLAB parity: predict, its four outputs and its scores
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! [label, score, node, cnum] = predict (Mdl, meas([1, 60, 120], :));
%! assert_equal (label, {'setosa'; 'versicolor'; 'virginica'});
%! assert_equal (node', [2, 8, 7]);
%! assert_equal (cnum', [1, 2, 3]);
%! assert_equal (score, [1, 0, 0; 0, 1, 0; 0, 1/3, 2/3], 1e-14);

%!test  # MATLAB parity: the losses on the training data
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! assert_equal (resubLoss (Mdl), 0.02, 1e-14);
%! assert_equal (loss (Mdl, meas, species), 0.02, 1e-14);
%! assert_equal (resubEdge (Mdl), 0.9384, 1e-4);

%!test  # MATLAB parity: a fit that keeps every split it made
%! load fisheriris
%! Mdl = fitctree (meas, species, 'MergeLeaves', 'off');
%! assert_equal (Mdl.NumNodes, 11);
%! assert_equal (isempty (Mdl.PruneList), false);

%!test  # MATLAB parity: turning both reductions off leaves no sequence
%! load fisheriris
%! Mdl = fitctree (meas, species, 'MergeLeaves', 'off', 'Prune', 'off');
%! assert_equal (Mdl.NumNodes, 11);
%! assert_equal (isempty (Mdl.PruneList), true);
%! assert_equal (isempty (Mdl.PruneAlpha), true);

%!test  # MATLAB parity: a merged tree carries a sequence with Prune off
%! load fisheriris
%! Mdl = fitctree (meas, species, 'Prune', 'off');
%! assert_equal (Mdl.PruneList', [4, 0, 3, 2, 0, 1, 0, 0, 0]);
%! assert_equal (Mdl.PruneAlpha', [0, 1/150, 2/150, 44/150, 50/150], 1e-14);

%!test  # MATLAB parity: MinParentSize is raised to two leaves
%! load fisheriris
%! Mdl = fitctree (meas, species, 'MinParentSize', 3, 'MinLeafSize', 2);
%! assert_equal (Mdl.ModelParameters.MinParent, 4);
%! assert_equal (Mdl.NodeSize', [150, 50, 100, 54, 46, 48, 6, 3, 3]);

%!test  # MATLAB parity: a leaf size that stops the tree early
%! load fisheriris
%! Mdl = fitctree (meas, species, 'MinLeafSize', 20);
%! assert_equal (Mdl.NumNodes, 5);
%! assert_equal (Mdl.NodeSize', [150, 50, 100, 54, 46]);

%!test  # MATLAB parity: a budget of two splits
%! load fisheriris
%! Mdl = fitctree (meas, species, 'MaxNumSplits', 2);
%! assert_equal (Mdl.NumNodes, 5);
%! assert_equal (Mdl.NodeSize', [150, 50, 100, 54, 46]);

%!test  # MATLAB parity: the deviance criterion and the risk it reports
%! load fisheriris
%! Mdl = fitctree (meas, species, 'SplitCriterion', 'deviance');
%! assert_equal (Mdl.NumNodes, 9);
%! assert_equal (Mdl.CutPoint([1, 3, 4, 6])', [2.45, 1.75, 4.95, 1.65]);
%! assert_equal (Mdl.NodeRisk(1), 0.792481250360578, 1e-12);
%! assert_equal (Mdl.NodeRisk(4), 0.0801, 1e-4);

%!test  # MATLAB parity: the prior a weighted fit reports
%! load fisheriris
%! Mdl = fitctree (meas, species, 'Weights', (1:150)');
%! assert_equal (Mdl.Prior, [1275, 3775, 6275] / 11325, 1e-14);
%! assert_equal (sum (Mdl.W), 1, 1e-14);
%! assert_equal (Mdl.NumNodes, 7);
%! assert_equal (Mdl.NodeSize', [150, 95, 55, 50, 45, 44, 1]);

%!test  # MATLAB parity: a uniform prior over an unbalanced sample
%! ## The shape of this tree is not asserted.  A uniform prior over 50, 20
%! ## and 10 observations makes the split that isolates setosa and the one
%! ## that isolates virginica exactly equal, at a gain of 1/3 each, and the
%! ## two engines keep opposite sides of the tie.
%! load fisheriris
%! inds = [1:50, 51:70, 101:110];
%! Mdl = fitctree (meas(inds, :), species(inds), 'Prior', 'uniform');
%! assert_equal (Mdl.Prior, [1/3, 1/3, 1/3], 1e-14);
%! assert_equal (Mdl.W([1, 51, 71])', [1/150, 1/60, 1/30], 1e-14);
%! assert_equal (sum (Mdl.W), 1, 1e-14);
%! assert_equal (Mdl.NodeSize(1), 80);
%! assert_equal (resubLoss (Mdl), 0, 1e-14);

%!test  # MATLAB parity: a cost matrix reshapes the tree
%! load fisheriris
%! Mdl = fitctree (meas, species, 'Cost', [0, 1, 10; 1, 0, 1; 10, 1, 0]);
%! assert_equal (Mdl.NumNodes, 9);
%! assert_equal (Mdl.NodeSize', [150, 50, 100, 45, 55, 44, 1, 9, 46]);
%! assert_equal (Mdl.CutPoint([1, 3, 4, 5])', [2.45, 4.75, 1.65, 1.75]);
%! assert_equal (Mdl.NodeClass{1}, 'versicolor');
%! assert_equal (Mdl.NodeRisk(1), 0.572916666666667, 1e-12);
%! assert_equal (Mdl.PruneAlpha', [0, 1/150, 43/150, 50/150], 1e-12);

%!test  # MATLAB parity: the label type is the response's own
%! load fisheriris
%! y = strcmp (species, 'setosa');
%! Mdl = fitctree (meas, y);
%! assert_equal (class (Mdl.ClassNames), 'logical');
%! assert_equal (class (predict (Mdl, meas(1, :))), 'logical');
%! Mdl = fitctree (meas, grp2idx (species));
%! assert_equal (Mdl.ClassNames, [1; 2; 3]);
%! assert_equal (Mdl.NodeClass{1}, '1');

%!test  # MATLAB parity: a missing response drops its row, a missing X does not
%! load fisheriris
%! x = meas;
%! x(1:10, 4) = NaN;
%! Mdl = fitctree (x, species);
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (Mdl.RowsUsed, []);
%! y = species;
%! y(1:10) = {''};
%! Mdl = fitctree (meas, y);
%! assert_equal (Mdl.NumObservations, 140);
%! assert_equal (sum (Mdl.RowsUsed), 140);

%!test  # MATLAB parity: ClassNames keeps only the classes it names
%! load fisheriris
%! Mdl = fitctree (meas, species, 'ClassNames', {'setosa', 'versicolor'});
%! assert_equal (Mdl.ClassNames, {'setosa'; 'versicolor'});
%! assert_equal (Mdl.NumObservations, 100);
%! assert_equal (Mdl.NodeSize', [100, 50, 50]);
%! assert_equal (Mdl.Prior, [0.5, 0.5], 1e-14);

## Test input validation
%!error<fitctree: too few arguments.> fitctree ()
%!error<fitctree: too few arguments.> fitctree (ones (4, 1))
%!error<fitctree: name-value arguments must be in pairs.>
%! fitctree (ones (4, 2), ones (4, 1), 'K')
%!error<fitctree: number of rows in X and Y must be equal.>
%! fitctree (ones (4, 2), ones (3, 1))
%!error<fitctree: number of rows in X and Y must be equal.>
%! fitctree (ones (4, 2), ones (3, 1), 'K', 2)
