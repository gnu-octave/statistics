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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{Mdl} =} fitrensemble (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitrensemble (@dots{}, @var{name}, @var{value})
##
## Fit an ensemble of regression trees.
##
## @code{@var{Mdl} = fitrensemble (@var{X}, @var{Y})} grows 100 regression
## trees by LSBoost on the @math{NxP} predictor matrix @var{X} and the numeric
## response @var{Y}, and returns a @code{RegressionEnsemble}.  With
## @qcode{'Method'} set to @qcode{'Bag'}, or with LSBoost resampling, it
## returns a @code{RegressionBaggedEnsemble}.  A row missing a predictor or the
## response is left out.
##
## LSBoost starts from a prediction of zero.  Each tree is fitted, with the
## observation weights, to the residual of the trees before it, and the
## prediction grows by the learning rate times that tree's prediction, the
## learning rate being the tree's weight.  The first tree therefore fits the
## response itself.  Bag grows each tree on a sample drawn in proportion to the
## weights and averages them.
##
## LSBoost resamples when @qcode{'Resample'} is @qcode{'on'} or
## @qcode{'FResample'} or @qcode{'Replace'} is given.  Each tree is then fitted
## to the residual of @code{ceil (FResample * N)} rows, drawn with replacement
## in proportion to the weights or without replacement uniformly, while the
## prediction and @code{FitInfo} run over every row, as in MATLAB R2024a.  The
## ensemble is a @code{RegressionBaggedEnsemble}, which records the rows each
## tree drew and estimates the out-of-bag error.
##
## Name-Value arguments:
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Name} @tab @tab @var{Value}
## @item @qcode{'Method'} @tab @tab @qcode{'LSBoost'} (default) or
## @qcode{'Bag'}.
## @item @qcode{'NumLearningCycles'} @tab @tab A positive integer, the number
## of trees to grow.  The default is 100.
## @item @qcode{'Learners'} @tab @tab @qcode{'tree'} (default) or a template
## from @code{templateTree}, whose options override the defaults: for LSBoost
## @code{MaxNumSplits} 10, @code{MinParentSize} 10 and @code{MinLeafSize} 5;
## for Bag unlimited splits, @code{MinParentSize} 10, @code{MinLeafSize} 5
## and @code{NumVariablesToSample} @code{ceil (P / 3)}.
## @item @qcode{'LearnRate'} @tab @tab A number greater than 0 and no greater
## than 1.  The default is 1.  LSBoost only.
## @item @qcode{'FResample'} @tab @tab The share of the observations each
## tree draws, greater than 0 and no greater than 1.  The default is 1.  Given
## with LSBoost, the ensemble resamples.
## @item @qcode{'Replace'} @tab @tab @qcode{'on'} (default) or @qcode{'off'},
## whether the trees draw with replacement.  Given with LSBoost, the ensemble
## resamples.
## @item @qcode{'Resample'} @tab @tab @qcode{'off'} (default) or
## @qcode{'on'}, whether LSBoost resamples.  Bag always does.
## @item @qcode{'NPrint'} @tab @tab @qcode{'off'} (default) or a positive
## integer @var{n}, to print a line after every @var{n} trees.
## @item @qcode{'Weights'} @tab @tab A nonnegative vector with one weight per
## observation.  The default is uniform.
## @item @qcode{'PredictorNames'} @tab @tab A cell array of character vectors
## naming the columns of @var{X}.
## @item @qcode{'ResponseName'} @tab @tab The name of the response variable.
## @item @qcode{'ResponseTransform'} @tab @tab @qcode{'none'} (default),
## @qcode{'exp'}, @qcode{'log'} or a function handle, applied to the
## predictions.  MATLAB R2024a accepts only a function handle here, failing
## on the named transforms.
## @end multitable
##
## @qcode{'CrossVal'} set to @qcode{'on'}, @qcode{'KFold'},
## @qcode{'Holdout'}, @qcode{'Leaveout'} or @qcode{'CVPartition'}, only one of
## them, fits the ensemble and cross-validates it as @code{crossval} does,
## returning a @code{RegressionPartitionedEnsemble}.
##
## Binning and hyperparameter optimization are not implemented, and an option
## asking for one of them is refused.  @qcode{'CategoricalPredictors'}, as
## indices, as a logical vector with one element per predictor, or as
## @qcode{'all'}, is passed on to every tree, which splits those predictors
## into sets of levels as @code{fitrtree} does.  An ensemble is
## regularized and shrunk afterwards with the @code{regularize}, @code{shrink}
## and @code{cvshrink} methods.
##
## @seealso{RegressionEnsemble, RegressionBaggedEnsemble,
## CompactRegressionEnsemble, templateTree, TreeBagger}
## @end deftypefn

function Mdl = fitrensemble (X, Y, varargin)

  if (nargin < 2)
    error ("fitrensemble: too few input arguments.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("fitrensemble: name-value arguments must be in pairs.");
  endif

  ## The cross-validation options make the fit a cross-validated one: the
  ## ensemble is fitted on all the data first, then on each fold.
  cv = {'crossval', 'kfold', 'holdout', 'leaveout', 'cvpartition'};
  iscv = false (size (varargin));
  method = [];
  resampled = false;
  for i = 1:2:numel (varargin)
    if (ischar (varargin{i}) && any (strcmpi (varargin{i}, cv)))
      iscv(i:i+1) = true;
    elseif (ischar (varargin{i}) && strcmpi (varargin{i}, 'Method'))
      method = varargin{i+1};
    elseif (ischar (varargin{i})
            && any (strcmpi (varargin{i}, {'FResample', 'Replace'})))
      resampled = true;
    elseif (ischar (varargin{i}) && strcmpi (varargin{i}, 'Resample'))
      resampled = resampled || (ischar (varargin{i+1})
                                && strcmpi (varargin{i+1}, 'on'));
    endif
  endfor
  cvargs = varargin(iscv);
  varargin = varargin(! iscv);

  ## LSBoost that resamples is a bagged ensemble too, as in MATLAB.
  isbag = ischar (method) && strcmpi (method, 'Bag');
  if (resampled && ! isbag)
    isbag = true;
    if (isempty (method))
      varargin(end+1:end+2) = {'Method', 'LSBoost'};
    endif
  endif

  if (isbag)
    Mdl = RegressionBaggedEnsemble (X, Y, varargin{:});
  else
    Mdl = RegressionEnsemble (X, Y, varargin{:});
  endif

  if (! isempty (cvargs))
    [P, errmsg] = ensemblePartition (cvargs, Mdl.Y, Mdl.NumObservations, false);
    if (! isempty (errmsg))
      error ("fitrensemble: %s", errmsg);
    endif
    if (! isempty (P))
      Mdl = RegressionPartitionedEnsemble (Mdl, P);
    endif
  endif

endfunction

%!demo
%! ## Boost regression trees to predict sepal length from the other three
%! ## measurements, and watch the training error fall as trees are added.
%! load fisheriris
%! X = meas(:,2:4);
%! y = meas(:,1);
%! Mdl = fitrensemble (X, y, 'NumLearningCycles', 50, 'LearnRate', 0.1);
%! plot (loss (Mdl, X, y, 'Mode', 'cumulative'));
%! xlabel ('Number of trees');
%! ylabel ('Training mean squared error');

%!demo
%! ## A bagged ensemble of regression trees predicts by averaging its trees.
%! load fisheriris
%! rng (42);
%! Mdl = fitrensemble (meas(:,2:4), meas(:,1), 'Method', 'Bag', ...
%!                     'NumLearningCycles', 30);
%! yfit = predict (Mdl, meas([1, 51, 101], 2:4))

## Test output
%!shared X, y, S
%! load fisheriris
%! X = meas(:,2:4);
%! y = meas(:,1);
%! S = templateTree ('MaxNumSplits', 1);

%!test  # MATLAB parity: LSBoost on decision stumps
%! Mdl = fitrensemble (X, y, 'NumLearningCycles', 4, 'Learners', S);
%! assert_equal (class (Mdl), 'RegressionEnsemble');
%! assert_equal (Mdl.TrainedWeights, ones (4, 1));
%! assert_equal (Mdl.FitInfo, [0.263279369032793; 0.185334478476685; ...
%!               0.176267398277275; 0.159714716160712], 1e-13);
%! assert_equal (Mdl.Trained{2}.CutPoint(1), 6.05, 1e-12);
%! assert_equal (Mdl.Trained{2}.NodeMean(2:3)', ...
%!               [-0.070535138620251, 1.105050505050500], 1e-12);

%!test  # MATLAB parity: the first tree fits the response itself
%! Mdl = fitrensemble (X, y, 'NumLearningCycles', 4, 'Learners', S);
%! assert_equal (Mdl.Trained{1}.NodeMean(1), mean (y), 1e-12);
%! assert_equal (predict (Mdl, X([1, 51],:)), ...
%!               [5.085793612916077; 6.335321158531990], 1e-13);

%!test  # MATLAB parity: the learning rate shrinks each tree's step
%! Mdl = fitrensemble (X, y, 'NumLearningCycles', 3, 'Learners', S, ...
%!                     'LearnRate', 0.1);
%! assert_equal (Mdl.TrainedWeights, 0.1 * ones (3, 1), 1e-15);
%! assert_equal (Mdl.FitInfo, [0.263279369032793; 0.258079095922867; ...
%!               0.253496179589799], 1e-13);
%! assert_equal (Mdl.Trained{2}.CutPoint(1), 3.95, 1e-12);
%! assert_equal (predict (Mdl, X([1, 51],:)), ...
%!               [1.391715765572467; 1.746807161347196], 1e-13);

%!test  # MATLAB parity: observation weights enter every tree
%! w = [5 * ones(50, 1); ones(100, 1)];
%! Mdl = fitrensemble (X, y, 'NumLearningCycles', 2, 'Learners', S, ...
%!                     'Weights', w);
%! assert_equal (Mdl.W([1, 51])', [5 / 350, 1 / 350], 1e-15);
%! assert_equal (Mdl.FitInfo, [0.185825250456633; 0.147589707296206], 1e-13);
%! assert_equal (predict (Mdl, X([1, 51],:)), ...
%!               [4.988922545877216; 6.342390194075585], 1e-13);

%!test  # MATLAB parity: LSBoost with its default trees
%! Mdl = fitrensemble (X, y, 'NumLearningCycles', 3);
%! assert_equal (Mdl.FitInfo, [0.079108234913235; 0.058774516561726; ...
%!               0.050203817681049], 1e-13);
%! assert_equal (predict (Mdl, X([1, 2, 51, 52, 101, 150],:)), ...
%!               [5.09375; 4.58949494949495; 6.874251054542892; ...
%!                6.212370889253637; 6.607394901394901; ...
%!                6.032614996997739], 1e-12);

%!test  # MATLAB parity: the defaults of a regression ensemble
%! Mdl = fitrensemble (X, y, 'NumLearningCycles', 2);
%! assert_equal (Mdl.Method, 'LSBoost');
%! assert_equal (Mdl.ModelParameters.LearnRate, 1);
%! assert_equal (Mdl.CombineWeights, 'WeightedSum');
%! assert_equal (Mdl.ResponseTransform, 'none');

%!test  # MATLAB parity: 'Bag' returns a bagged ensemble
%! Mdl = fitrensemble (X, y, 'Method', 'Bag', 'NumLearningCycles', 2);
%! assert_equal (class (Mdl), 'RegressionBaggedEnsemble');
%! assert_equal (Mdl.CombineWeights, 'WeightedAverage');

%!test  # MATLAB parity: progress is printed after every NPrint trees
%! out = evalc (["fitrensemble (X, y, 'NumLearningCycles', 4, ", ...
%!               "'Learners', S, 'NPrint', 2);"]);
%! assert_equal (out, sprintf (["Training LSBoost...\n", ...
%!                              "Grown weak learners: 2\n", ...
%!                              "Grown weak learners: 4\n"]));

## Test input validation
%!error<fitrensemble: too few input arguments.> fitrensemble (X)
%!error<fitrensemble: name-value arguments must be in pairs.> ...
%! fitrensemble (X, y, 'Method')

%!test  # MATLAB parity: 'KFold' returns a cross-validated ensemble
%! CV = fitrensemble (X, y, 'NumLearningCycles', 2, 'Learners', S, 'KFold', 3);
%! assert_equal (class (CV), 'RegressionPartitionedEnsemble');
%! assert_equal (CV.KFold, 3);

%!error<fitrensemble: specify only one of 'CrossVal', 'KFold', 'Holdout', 'Leaveout' and 'CVPartition'.> ...
%! fitrensemble (X, y, 'KFold', 5, 'CrossVal', 'on')
%!error<fitrensemble: 'KFold' must be an integer greater than 1.> ...
%! fitrensemble (X, y, 'KFold', 0)

%!test  # MATLAB parity: LSBoost that resamples is a bagged ensemble
%! load fisheriris
%! M = fitrensemble (meas(:,2:4), meas(:,1), 'NumLearningCycles', 3, ...
%!                   'FResample', 0.5);
%! assert_equal (class (M), 'RegressionBaggedEnsemble');
%! assert_equal (M.Method, 'LSBoost');
%! assert_equal (M.CombineWeights, 'WeightedSum');

%!shared X, yr, Xq
%! k = (0:119)';
%! c = mod (k, 4) + 1;
%! j = floor (k / 4);
%! x2 = mod (k * 7, 10);
%! X = [c, x2];
%! yb = (c == 1) | (c == 2 & mod (j, 4) != 0) | (c == 3 & mod (j, 4) == 0) ...
%!      | (x2 > 7);
%! yr = [3; 1; 4; 1.5];
%! yr = yr(c) + 0.1 * sin (k) + 0.2 * x2;
%! Xq = [1, 0; 3, 5; 5, 0; NaN, 2; 2.5, 9];

%!test  # MATLAB parity: LSBoost trees split a categorical predictor
%! Mdl = fitrensemble (X, yr, 'Method', 'LSBoost', 'NumLearningCycles', 4, ...
%!                     'CategoricalPredictors', 1);
%! assert_equal (Mdl.CategoricalPredictors, 1);
%! assert_equal (Mdl.Trained{1}.CutCategories(1,:), {[2, 4], [1, 3]});
%! ## The first two rows pass a node of the fourth tree where two partitions
%! ## gain the same to 5e-15, which rounding decides; the others do not.
%! yhat = predict (Mdl, Xq);
%! assert_equal (yhat(3:5)', [3.2692634, 3.2028144, 3.4104563], 1e-6);
