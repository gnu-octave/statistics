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
## @deftypefn  {statistics} {@var{Mdl} =} fitcensemble (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcensemble (@dots{}, @var{name}, @var{value})
##
## Fit an ensemble of decision trees for classification.
##
## @code{@var{Mdl} = fitcensemble (@var{X}, @var{Y})} grows 100 boosted
## decision trees on the @math{NxP} predictor matrix @var{X} and the class
## labels @var{Y}, LogitBoost for two classes and AdaBoostM2 for more, and
## returns a @code{ClassificationEnsemble}.  With @qcode{'Method'} set to
## @qcode{'Bag'} it returns a @code{ClassificationBaggedEnsemble}.
##
## @var{Y} holds a class label per row of @var{X}, as a numeric or logical
## vector, a categorical, string or character array, or a cell array of
## character vectors.  A row missing a predictor or a class is left out.
##
## The boosting methods, @var{y} being +1 for the first class and -1 for the
## second, @var{d} the observation weights, which start at @code{W} times the
## total cost of misclassifying each observation's class, and @var{eta} the
## learning rate:
##
## @table @asis
## @item @qcode{'AdaBoostM1'}
## Two classes.  Each tree is grown with the weights @var{d}; @var{h} is +1
## where it predicts the first class and -1 elsewhere, @var{e} its weighted
## error, and its weight @code{eta * log ((1 - e) / e) / 2}.  The weights are
## then multiplied by @code{exp (-weight * y .* h)}.
## @item @qcode{'AdaBoostM2'}
## More than two classes.  A weight is kept for each observation and each
## class other than its own, and each tree is grown with their sums; its
## pseudo-loss @var{e} over its class probabilities gives its weight as for
## AdaBoostM1, and the scores are the weighted sums of the class
## probabilities.
## @item @qcode{'GentleBoost'}
## Two classes.  Each regression tree is fitted to @var{y} with the weights
## @var{d}, its prediction @var{h} added to the score times @var{eta}, and the
## weights multiplied by @code{exp (-eta * y .* h)}.
## @item @qcode{'LogitBoost'}
## Two classes.  With @var{p} the probability of the first class, each
## regression tree is fitted to @code{(y01 - p) ./ (p .* (1 - p))} with the
## weights @code{d .* p .* (1 - p)}; its prediction times @code{eta / 2} is
## added to the score @var{f}, and @code{p = 1 ./ (1 + exp (-f))}.
## @end table
##
## A boosting method whose tree classifies the data without error, zero error
## for AdaBoostM1 or zero pseudo-loss for AdaBoostM2, keeps that tree with the
## weight an error of @code{eps} gives and stops; one whose error is greater
## than 0.5 is not kept and stops.  MATLAB discards a perfect tree, so a fit
## whose first tree is perfect is empty there.  GentleBoost applies the
## learning rate to the fit itself, as the other boosting methods do; MATLAB
## scales only the learner weights, growing the same trees whatever the
## rate.
##
## Name-Value arguments:
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Name} @tab @tab @var{Value}
## @item @qcode{'Method'} @tab @tab @qcode{'AdaBoostM1'},
## @qcode{'AdaBoostM2'}, @qcode{'GentleBoost'}, @qcode{'LogitBoost'} or
## @qcode{'Bag'}.  The default is @qcode{'LogitBoost'} for two classes and
## @qcode{'AdaBoostM2'} for more.
## @item @qcode{'NumLearningCycles'} @tab @tab A positive integer, the number
## of trees to grow.  The default is 100.
## @item @qcode{'Learners'} @tab @tab @qcode{'tree'} (default) or a template
## from @code{templateTree}, whose options override the defaults: for
## boosting @code{MaxNumSplits} 10, @code{MinParentSize} 2 and
## @code{MinLeafSize} 1, the regression trees of GentleBoost and LogitBoost
## taking @code{MinParentSize} 10; for Bag unlimited splits,
## @code{MinParentSize} 2, @code{MinLeafSize} 1 and
## @code{NumVariablesToSample} @code{ceil (sqrt (P))}.
## @item @qcode{'LearnRate'} @tab @tab A number greater than 0 and no greater
## than 1.  The default is 1.  Boosting only.
## @item @qcode{'FResample'} @tab @tab The share of the observations each
## bagged tree draws, greater than 0 and no greater than 1.  The default is
## 1.  Bag only.
## @item @qcode{'Replace'} @tab @tab @qcode{'on'} (default) or @qcode{'off'},
## whether the bagged trees draw with replacement.  Bag only.
## @item @qcode{'NPrint'} @tab @tab @qcode{'off'} (default) or a positive
## integer @var{n}, to print a line after every @var{n} trees.
## @item @qcode{'ClassNames'} @tab @tab The classes to fit, in the order their
## scores are laid out.
## @item @qcode{'Cost'} @tab @tab A @math{KxK} matrix of misclassification
## costs.  The default is @code{1 - eye (K)}.
## @item @qcode{'Prior'} @tab @tab @qcode{'empirical'} (default),
## @qcode{'uniform'}, or a vector with one element per class.
## @item @qcode{'Weights'} @tab @tab A nonnegative vector with one weight per
## observation.  The default is uniform.
## @item @qcode{'PredictorNames'} @tab @tab A cell array of character vectors
## naming the columns of @var{X}.
## @item @qcode{'ResponseName'} @tab @tab The name of the response variable.
## @item @qcode{'ScoreTransform'} @tab @tab A transform applied to the
## returned scores.  The default is @qcode{'none'}.
## @end multitable
##
## The methods @qcode{'Subspace'}, @qcode{'LPBoost'}, @qcode{'TotalBoost'},
## @qcode{'CrossVal'} set to @qcode{'on'}, @qcode{'KFold'},
## @qcode{'Holdout'}, @qcode{'Leaveout'} or @qcode{'CVPartition'}, only one of
## them, fits the ensemble and cross-validates it as @code{crossval} does,
## returning a @code{ClassificationPartitionedEnsemble}.
##
## @qcode{'RobustBoost'} and @qcode{'RUSBoost'}, resampling in a boosting
## method, categorical predictors, binning and hyperparameter optimization
## are not implemented, and an option asking for one of them is refused.
##
## @seealso{ClassificationEnsemble, ClassificationBaggedEnsemble,
## CompactClassificationEnsemble, templateTree, TreeBagger}
## @end deftypefn

function Mdl = fitcensemble (X, Y, varargin)

  if (nargin < 2)
    error ("fitcensemble: too few input arguments.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("fitcensemble: name-value arguments must be in pairs.");
  endif

  ## The cross-validation options make the fit a cross-validated one: the
  ## ensemble is fitted on all the data first, then on each fold.
  cv = {'crossval', 'kfold', 'holdout', 'leaveout', 'cvpartition'};
  iscv = false (size (varargin));
  isbag = false;
  for i = 1:2:numel (varargin)
    if (ischar (varargin{i}) && any (strcmpi (varargin{i}, cv)))
      iscv(i:i+1) = true;
    elseif (ischar (varargin{i}) && strcmpi (varargin{i}, 'Method'))
      isbag = ischar (varargin{i+1}) && strcmpi (varargin{i+1}, 'Bag');
    endif
  endfor
  cvargs = varargin(iscv);
  varargin = varargin(! iscv);

  if (isbag)
    Mdl = ClassificationBaggedEnsemble (X, Y, varargin{:});
  else
    Mdl = ClassificationEnsemble (X, Y, varargin{:});
  endif

  if (! isempty (cvargs))
    [P, errmsg] = ensemblePartition (cvargs, Mdl.Y, Mdl.NumObservations, true);
    if (! isempty (errmsg))
      error ("fitcensemble: %s", errmsg);
    endif
    if (! isempty (P))
      Mdl = ClassificationPartitionedEnsemble (Mdl, P);
    endif
  endif

endfunction

%!demo
%! ## Boost decision stumps to tell versicolor from virginica, and see how
%! ## the training error falls as the stumps accumulate.
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! Mdl = fitcensemble (X, Y, 'Method', 'AdaBoostM1', ...
%!                     'NumLearningCycles', 20, ...
%!                     'Learners', templateTree ('MaxNumSplits', 1));
%! err = loss (Mdl, X, Y, 'Mode', 'cumulative');
%! plot (err);
%! xlabel ('Number of stumps');
%! ylabel ('Training error');

%!demo
%! ## A bagged ensemble of trees classifies all three species, and scores
%! ## each flower with the average of its trees' class probabilities.
%! load fisheriris
%! rng (42);
%! Mdl = fitcensemble (meas, species, 'Method', 'Bag', ...
%!                     'NumLearningCycles', 30);
%! [label, scores] = predict (Mdl, [5.0, 3.4, 1.5, 0.2; 6.7, 3.0, 5.2, 2.3])

## Test output
%!shared X2, Y2, S
%! load fisheriris
%! X2 = meas(51:150,:);
%! Y2 = species(51:150);
%! S = templateTree ('MaxNumSplits', 1);

%!test  # MATLAB parity: AdaBoostM1 on decision stumps
%! Mdl = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                     'NumLearningCycles', 5, 'Learners', S);
%! assert_equal (Mdl.TrainedWeights, [1.375767656520975; 0.993534110774411; ...
%!               0.883434373403998; 0.554364192108758; ...
%!               0.268194447432047], 1e-13);
%! assert_equal (Mdl.FitInfo, [0.06; 0.120567375886525; 0.145932163187856; ...
%!               0.248108033048671; 0.369028016476793], 1e-13);
%! cuts = cellfun (@(t) t.CutPoint(1), Mdl.Trained);
%! assert_equal (cuts, [1.75; 4.95; 4.45; 4.95; 5.15], 1e-12);

%!test  # MATLAB parity: the scores of AdaBoostM1 are signed learner weights
%! Mdl = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                     'NumLearningCycles', 5, 'Learners', S);
%! [label, s] = predict (Mdl, X2([1, 21, 51],:));
%! assert_equal (label, {'versicolor'; 'virginica'; 'virginica'});
%! assert_equal (s(:,1), [1.772037138568098; -0.979498174473852; ...
%!                        -2.966566396022673], 1e-13);
%! assert_equal (s(:,2), -s(:,1));

%!test  # MATLAB parity: the learning rate shrinks AdaBoostM1
%! Mdl = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                     'NumLearningCycles', 4, 'Learners', S, ...
%!                     'LearnRate', 0.5);
%! assert_equal (Mdl.TrainedWeights, [0.687883828260487; 0.595450620008897; ...
%!               0.406803001473525; 0.320537560668013], 1e-13);
%! assert_equal (Mdl.FitInfo(2:4), [0.084570915580601; 0.164212646508385; ...
%!               0.217184427004775], 1e-13);

%!test  # MATLAB parity: a uniform prior sets the starting weights
%! load fisheriris
%! Yv = strcmp (species(41:150), 'virginica');
%! Mdl = fitcensemble (meas(41:150,:), Yv, ...
%!                     'Method', 'AdaBoostM1', 'NumLearningCycles', 3, ...
%!                     'Learners', S, 'Prior', 'uniform');
%! assert_equal (Mdl.W([1, 110])', [0.5 / 60, 0.5 / 50], 1e-15);
%! assert_equal (Mdl.TrainedWeights, [1.406116776793511; 1.083890239183778; ...
%!               0.818511561052510], 1e-13);

%!test  # MATLAB parity: a cost multiplies the starting weights
%! Mdl = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                     'NumLearningCycles', 3, 'Learners', S, ...
%!                     'Cost', [0, 1; 5, 0]);
%! assert_equal (Mdl.W([1, 51])', [0.01, 0.01], 1e-15);
%! assert_equal (Mdl.FitInfo, [0.036666666666667; 0.125511167033658; ...
%!               0.225018706839040], 1e-13);
%! [~, s] = predict (Mdl, X2(1,:));
%! assert_equal (s(1), 3.223215778029441, 1e-13);

%!test  # MATLAB parity: AdaBoostM2 on decision stumps
%! load fisheriris
%! Mdl = fitcensemble (meas, species, 'Method', 'AdaBoostM2', ...
%!                     'NumLearningCycles', 4, 'Learners', S);
%! assert_equal (Mdl.TrainedWeights, [0.549306144334059; 0.496958723838212; ...
%!               0.495870843211405; 0.456824643135089], 1e-13);
%! assert_equal (Mdl.FitInfo, [0.25; 0.270139003161550; 0.270568199265314; ...
%!               0.286253661231616], 1e-13);
%! [~, s] = predict (Mdl, meas([1, 101],:));
%! assert_equal (s, [1.393875268391962, 0.557069431281369, ...
%!                   0.048015654845434; 0, 0.600192648565954, ...
%!                   1.398767705952811], 1e-13);

%!test  # MATLAB parity: GentleBoost on decision stumps
%! Mdl = fitcensemble (X2, Y2, 'Method', 'GentleBoost', ...
%!                     'NumLearningCycles', 3, 'Learners', S);
%! assert_equal (class (Mdl.Trained{1}), 'CompactRegressionTree');
%! assert_equal (Mdl.TrainedWeights, [1; 1; 1]);
%! assert_equal (Mdl.FitInfo, [0.220611916264090; 0.327133481554914; ...
%!               0.606925714189402], 1e-13);
%! [~, s] = predict (Mdl, X2([1, 21, 51],:));
%! assert_equal (s(:,1), [2.112004327552702; -0.903893073252190; ...
%!                        -2.537346885979990], 1e-13);

%!test  # MATLAB parity: GentleBoost with its default trees
%! Mdl = fitcensemble (X2, Y2, 'Method', 'GentleBoost', ...
%!                     'NumLearningCycles', 4);
%! assert_equal (Mdl.FitInfo, [0.08; 0.157657738520697; 0.188119628350251; ...
%!               0.290921726710326], 1e-13);
%! assert_equal (resubLoss (Mdl, 'LossFun', 'exponential'), ...
%!               0.054188751279988, 1e-13);
%! assert_equal (resubEdge (Mdl), 7.127687492916048, 1e-12);

%!test  # the learning rate shrinks GentleBoost, where MATLAB only rescales
%! Mdl = fitcensemble (X2, Y2, 'Method', 'GentleBoost', ...
%!                     'NumLearningCycles', 3, 'Learners', S, ...
%!                     'LearnRate', 0.5);
%! assert_equal (Mdl.TrainedWeights, [0.5; 0.5; 0.5]);
%! assert_equal (Mdl.FitInfo(1), 0.220611916264090, 1e-13);
%! assert_equal (abs (Mdl.FitInfo(2) - 0.327133481554914) > 1e-3, true);

%!test  # MATLAB parity: LogitBoost on decision stumps
%! Mdl = fitcensemble (X2, Y2, 'Method', 'LogitBoost', ...
%!                     'NumLearningCycles', 3, 'Learners', S);
%! assert_equal (Mdl.TrainedWeights, [0.5; 0.5; 0.5]);
%! assert_equal (Mdl.FitInfo, [0.882447665056360; 0.883091315814675; ...
%!               1.164564532959708], 1e-13);
%! [~, s] = predict (Mdl, X2([1, 21, 51],:));
%! assert_equal (s(:,1), [1.906010018452977; 0.134673464507728; ...
%!                        -2.170659689353640], 1e-13);

%!test  # MATLAB parity: the learning rate shrinks LogitBoost
%! Mdl = fitcensemble (X2, Y2, 'Method', 'LogitBoost', ...
%!                     'NumLearningCycles', 3, 'Learners', S, ...
%!                     'LearnRate', 0.5);
%! assert_equal (Mdl.TrainedWeights, [0.25; 0.25; 0.25]);
%! assert_equal (Mdl.FitInfo, [0.882447665056360; 0.826287268663811; ...
%!               0.840635802435344], 1e-13);

%!test  # MATLAB parity: LogitBoost with its default trees
%! Mdl = fitcensemble (X2, Y2, 'Method', 'LogitBoost', ...
%!                     'NumLearningCycles', 6);
%! assert_equal (Mdl.FitInfo, [0.32; 0.341307860560587; ...
%!               0.516570570389474; 0.548550758138539; ...
%!               0.890006883416469; 0.341774284666034], 1e-12);
%! [~, s] = predict (Mdl, X2([1, 2, 21, 100],:));
%! assert_equal (s(:,1), [4.582296394974635; 4.239663719918038; ...
%!                        0.504944771289154; -2.839026377180836], 1e-12);

%!test  # MATLAB parity: an error of exactly one half does not stop the fit
%! X = repmat ([0, 0; 0, 1; 1, 0; 1, 1], 10, 1);
%! Y = repmat ([0; 1; 1; 0], 10, 1);
%! Mdl = fitcensemble (X, Y, 'Method', 'AdaBoostM1', ...
%!                     'NumLearningCycles', 3, 'Learners', S);
%! assert_equal (Mdl.NumTrained, 3);
%! assert_equal (Mdl.FitInfo, [0.5; 0.5; 0.5], 1e-14);

%!test  # a perfect learner is kept and ends the fit, where MATLAB drops it
%! load fisheriris
%! Mdl = fitcensemble (meas, strcmp (species, 'setosa'), ...
%!                     'Method', 'AdaBoostM1', 'NumLearningCycles', 3, ...
%!                     'Learners', S);
%! assert_equal (Mdl.NumTrained, 1);
%! assert_equal (Mdl.ReasonForTermination, ...
%!               'Classification error from the last weak learner is zero.');
%! assert_equal (Mdl.TrainedWeights, log ((1 - eps) / eps) / 2, 1e-12);
%! assert_equal (resubLoss (Mdl), 0);

%!test  # a perfect AdaBoostM2 learner is kept too
%! load fisheriris
%! Mdl = fitcensemble (meas, species, 'NumLearningCycles', 5);
%! assert_equal (Mdl.Method, 'AdaBoostM2');
%! assert_equal (Mdl.NumTrained, 1);
%! assert_equal (Mdl.ReasonForTermination, ...
%!               'Pseudo-loss from the last weak learner is zero.');
%! assert_equal (resubLoss (Mdl), 0);

%!test  # MATLAB parity: LogitBoost is the default for two classes
%! Mdl = fitcensemble (X2, Y2, 'NumLearningCycles', 2);
%! assert_equal (class (Mdl), 'ClassificationEnsemble');
%! assert_equal (Mdl.Method, 'LogitBoost');
%! assert_equal (Mdl.CombineWeights, 'WeightedSum');

%!test  # MATLAB parity: 'Bag' returns a bagged ensemble
%! load fisheriris
%! Mdl = fitcensemble (meas, species, 'Method', 'Bag', ...
%!                     'NumLearningCycles', 2);
%! assert_equal (class (Mdl), 'ClassificationBaggedEnsemble');
%! assert_equal (Mdl.CombineWeights, 'WeightedAverage');

%!test  # MATLAB parity: progress is printed after every NPrint learners
%! out = evalc (["fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ", ...
%!               "'NumLearningCycles', 4, 'Learners', S, 'NPrint', 2);"]);
%! assert_equal (out, sprintf (["Training AdaBoostM1...\n", ...
%!                              "Grown weak learners: 2\n", ...
%!                              "Grown weak learners: 4\n"]));

## Test input validation
%!error<fitcensemble: too few input arguments.> fitcensemble (X2)
%!error<fitcensemble: name-value arguments must be in pairs.> ...
%! fitcensemble (X2, Y2, 'Method')

%!test  # MATLAB parity: 'CrossVal' gives ten folds
%! CV = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                    'NumLearningCycles', 2, 'Learners', S, 'CrossVal', 'on');
%! assert_equal (class (CV), 'ClassificationPartitionedEnsemble');
%! assert_equal (CV.KFold, 10);
%! M = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                   'NumLearningCycles', 2, 'Learners', S, 'CrossVal', 'off');
%! assert_equal (class (M), 'ClassificationEnsemble');

%!error<fitcensemble: specify only one of 'CrossVal', 'KFold', 'Holdout', 'Leaveout' and 'CVPartition'.> ...
%! fitcensemble (X2, Y2, 'KFold', 5, 'Holdout', 0.2)
%!error<fitcensemble: 'KFold' must be an integer greater than 1.> ...
%! fitcensemble (X2, Y2, 'KFold', 1)
%!error<fitcensemble: 'Holdout' must be a number between 0 and 1.> ...
%! fitcensemble (X2, Y2, 'Holdout', 1)
%!error<fitcensemble: 'CrossVal' must be 'on' or 'off'.> ...
%! fitcensemble (X2, Y2, 'CrossVal', true)
%!error<fitcensemble: 'Leaveout' must be 'on' or 'off'.> ...
%! fitcensemble (X2, Y2, 'Leaveout', 1)
%!error<fitcensemble: 'CVPartition' must be a 'cvpartition' object.> ...
%! fitcensemble (X2, Y2, 'CVPartition', 5)
%!error<fitcensemble: 'CVPartition' must partition the observations the ensemble was fitted on.> ...
%! fitcensemble (X2, Y2, 'CVPartition', cvpartition (50, 'KFold', 5))
