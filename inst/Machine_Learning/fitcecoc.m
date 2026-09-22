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
## @deftypefn  {statistics} {@var{Mdl} =} fitcecoc (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcecoc (@var{Tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{Mdl} =} fitcecoc (@var{Tbl}, @var{formula})
## @deftypefnx {statistics} {@var{Mdl} =} fitcecoc (@var{Tbl}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcecoc (@dots{}, @var{name}, @var{value})
##
## Fit a multiclass model from binary learners.
##
## @code{@var{Mdl} = fitcecoc (@var{X}, @var{Y})} turns the multiclass
## problem in @var{X} and @var{Y} into a set of two class problems, fits one
## binary learner to each, and returns them as a @code{ClassificationECOC}
## object.  A coding matrix says which classes each learner calls +1, which
## it calls -1, and which sit it out; a new observation is sent to every
## learner and given the class whose row of that matrix its scores match most
## closely.
##
## @itemize
## @item
## @var{X} must be a @math{NxP} numeric matrix of predictor data.
## @item
## @var{Y} must be a @math{Nx1} vector of class labels, of any type
## @code{ClassNames} accepts, with at least two distinct values.
## @end itemize
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{'Learners'} @tab @tab The binary learner, either a name,
## @qcode{'svm'} (default), @qcode{'tree'}, @qcode{'knn'},
## @qcode{'naivebayes'}, @qcode{'discriminant'}, @qcode{'linear'},
## @qcode{'kernel'} or @qcode{'ensemble'}, a LogitBoost ensemble of 100
## trees, or a template from @code{templateSVM} and its siblings,
## @code{templateEnsemble} among them, which also carries the options that
## learner is to be fitted with.
##
## @item @qcode{'Coding'} @tab @tab The coding design, either a name
## @code{designecoc} accepts, @qcode{'onevsone'} by default, or a coding
## matrix given outright, which sets @code{CodingName} to @qcode{'custom'}.
##
## @item @qcode{'BinaryLoss'} @tab @tab The loss the binary scores are read
## with.  The default follows the learner: @qcode{'exponential'} for an
## AdaBoostM1 or GentleBoost ensemble, @qcode{'binodeviance'} for a
## LogitBoost one, and otherwise @qcode{'hinge'} for one scoring on
## @math{(-Inf,+Inf)} and @qcode{'quadratic'} for one scoring on @math{[0,1]},
## as a bagged, random subspace or RUSBoost ensemble does.
##
## @item @qcode{'ClassNames'} @tab @tab The classes to fit, and the order
## the rows of the coding matrix, @code{Prior} and @code{Cost} take them in.
## By default they are sorted.
##
## @item @qcode{'Cost'} @tab @tab A @math{KxK} matrix of misclassification
## costs.  The default is @code{1 - eye (K)}.
##
## @item @qcode{'Prior'} @tab @tab @qcode{'empirical'} (default),
## @qcode{'uniform'}, or a vector with one element per class.
##
## @item @qcode{'Weights'} @tab @tab A nonnegative numeric vector with one
## element per observation.  The default is uniform.
##
## @item @qcode{'CategoricalPredictors'} @tab @tab The predictors whose values
## are levels, as indices, a logical vector or @qcode{'all'}, passed as given
## to every binary learner, which codes them its own way; a nearest neighbour
## learner takes only @qcode{'all'} and a discriminant learner none.
## A predictor may be named rather than indexed, as a character matrix of one
## padded name per row, a string array or a cellstr; a name must match an entry
## of @qcode{'PredictorNames'} exactly, its case included.
##
## @item @qcode{'PredictorNames'} @tab @tab A cellstr of predictor names.
##
## @item @qcode{'ResponseName'} @tab @tab The name of the response variable.
##
## @item @qcode{'ScoreTransform'} @tab @tab A transform applied to the
## returned scores.  The default is @qcode{'none'}.
## @end multitable
##
## A linear or kernel learner carries no training data of its own, so
## @code{fitcecoc} returns a @code{CompactClassificationECOC} for those two
## and a @code{ClassificationECOC} for the other five, which is what R2024a
## does.
##
## @code{'FitPosterior'} is refused rather than quietly ignored: it installs
## a fitted score transform on each binary learner, which needs the posterior
## fitting this package does not have yet.
##
## An ensemble template may name any classification method of
## @code{fitcensemble}; a regression template is refused, as in MATLAB.
##
## A RUSBoost ensemble scores each class with the weighted sum of its trees'
## class probabilities, which no binary loss can read: MATLAB sets the loss
## to @qcode{'unknown'} and then cannot predict.  Here those scores, and the
## binary scores @code{predict} returns, are divided by the total weight of
## the trees, which gives the weighted mean of their class probabilities, and
## are read with @qcode{'quadratic'} as a bagged ensemble's are.  This follows
## R's @code{ebmc} and @code{adabag} packages and scikit-learn, whose
## @code{AdaBoostClassifier} and imbalanced-learn's @code{RUSBoostClassifier}
## built on it scale a boosted ensemble's weighted votes the same way.
##
## A binary learner that cannot be fitted, such as an AdaBoostM2 ensemble, which
## needs three classes, stops the fit with its own error, where MATLAB warns and
## predicts the majority class.  An AdaBoostM1 learner whose first tree
## separates its two classes keeps that tree, as @code{fitcensemble} does;
## MATLAB keeps no tree and never predicts those classes from it.
##
## Each binary learner is fitted with the weight its observations carry, which
## the prior and any @qcode{'Weights'} decide.  The SVM, KNN, naive Bayes and
## discriminant learners of this package take no observation weights, so when
## those weights are unequal they are given instead the prior their two sides
## hold, the prior MATLAB's learners report, and weights that vary within a
## class, which no prior can express, are refused for them.  MATLAB weighs
## each observation of those learners as well.
##
## @seealso{ClassificationECOC, CompactClassificationECOC, designecoc,
## templateSVM, templateTree}
## @end deftypefn

function Mdl = fitcecoc (X, Y, varargin)

  ## Input validation
  if (nargin < 2)
    error ("fitcecoc: too few input arguments.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("fitcecoc: name-value arguments must be in pairs.");
  endif
  if (! istable (X) && rows (X) != rows (Y))
    error ("fitcecoc: number of rows in X and Y must be equal.");
  endif

  ## A cross-validation option makes the fit a cross-validated one: the
  ## model is fitted on all the data first, then partitioned, which is how
  ## every other fitc* here reads these four.
  cvnames = {'kfold', 'holdout', 'leaveout', 'cvpartition'};
  iscv = cellfun (@(a) ischar (a) && isrow (a) ...
                       && any (strcmpi (a, cvnames)), varargin(1:2:end));
  cvargs = {};
  if (any (iscv))
    if (sum (iscv) > 1)
      error (strcat ("fitcecoc: specify only one of 'KFold', 'Holdout',", ...
                     " 'Leaveout' and 'CVPartition'."));
    endif
    at = 2 * find (iscv) - 1;
    cvargs = varargin(at:at+1);
    varargin(at:at+1) = [];
  endif

  Mdl = ClassificationECOC (X, Y, varargin{:});

  if (! isempty (cvargs))
    Mdl = crossval (Mdl, cvargs{:});
    return;
  endif

  ## A learner that keeps no training data leaves nothing for the full class
  ## to hold, so the fit gives back the compact one.  Measured on R2024a.
  if (any (strcmpi (Mdl.ModelParameters.BinaryLearners.Method, ...
                    {'Linear', 'Kernel'})))
    Mdl = compact (Mdl);
  endif

endfunction

## Tests
%!demo
%! ## Fit from a table, and predict on one
%!
%! load fisheriris
%! T = table (meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!            'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! T.Species = categorical (species);
%!
%! ## A column holding levels is a categorical predictor without being named
%! ## one
%! T.Wide = categorical (meas(:,2) > 3, [false true], {'narrow', 'wide'});
%!
%! ## The response is named by its column, and everything else is a predictor
%! Mdl = fitcecoc (T, 'Species');
%! Mdl.PredictorNames
%! Mdl.CategoricalPredictors
%!
%! ## A model formula names them instead, holding main effects only
%! Mdl2 = fitcecoc (T, 'Species ~ PL + PW');
%! Mdl2.PredictorNames
%!
%! ## predict reads a table by the names the model was fitted on, so the
%! ## columns may come in any order
%! label = predict (Mdl, T(1:5, [6, 5, 4, 3, 2, 1]));
%! label'

%!demo
%! ## Score a table
%!
%! load fisheriris
%! T = table (meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!            'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! T.Species = categorical (species);
%! Mdl = fitcecoc (T, 'Species');
%!
%! ## The response is named by its column, or left out, when it is the
%! ## variable the model was fitted on
%! [loss(Mdl, T, 'Species'), loss(Mdl, T)]
%!
%! ## It may also be given beside a table holding the predictors alone
%! loss (Mdl, T(:, 1:4), T.Species)
%!
%! ## A name-value argument does not stand in for the response: an even
%! ## number of arguments after the table is all name-value, an odd one
%! ## names the response first
%! loss (Mdl, T, 'LossFun', 'classiferror')

%!test  # MATLAB parity: the default fit and what it reports
%! load fisheriris
%! Mdl = fitcecoc (meas, species);
%! assert_equal (class (Mdl), 'ClassificationECOC');
%! assert_equal (numel (properties (Mdl)), 22);
%! assert_equal (Mdl.CodingName, 'onevsone');
%! assert_equal (Mdl.BinaryLoss, 'hinge');
%! assert_equal (Mdl.CodingMatrix, [1, 1, 0; -1, 0, 1; 0, -1, -1]);
%! assert_equal (Mdl.LearnerWeights, [2/3, 2/3, 2/3], 1e-12);
%! assert_equal (Mdl.NumObservations, 150);

%!test  # MATLAB parity: the loss of the default fit
%! ## Our binary SVM is LIBSVM where MATLAB's is SMO, so the scores differ in
%! ## the fourth digit on the pair that is not separable; the labels they
%! ## lead to do not.
%! load fisheriris
%! Mdl = fitcecoc (meas, species);
%! assert_equal (resubLoss (Mdl), 0.0066666666666667, 1e-12);

%!test  # MATLAB parity: a tree code reproduces the decoding exactly
%! ## Our trees match MATLAB's on this fixture, so the whole path can be
%! ## compared and not only the labels.  Measured on R2024a.
%! load fisheriris
%! Mdl = fitcecoc (meas, species, 'Learners', 'tree');
%! assert_equal (Mdl.BinaryLoss, 'quadratic');
%! [~, NegLoss] = predict (Mdl, meas([1, 20, 51, 70, 101, 130], :));
%! assert_equal (NegLoss, ...
%!               [0, -1, -2; ...
%!                0, -1, -2; ...
%!                -2, 0, -1; ...
%!                -2, 0, -1; ...
%!                -2, -0.956994328922495, -0.000472589792060491; ...
%!                -2, -0.444444444444445, -0.111111111111111], 1e-12);

%!test  # MATLAB parity: the default binary loss follows the learner
%! load fisheriris
%! for p = {{'svm', 'hinge'}, {'tree', 'quadratic'}, {'knn', 'quadratic'}, ...
%!          {'naivebayes', 'quadratic'}, {'discriminant', 'quadratic'}}
%!   Mdl = fitcecoc (meas, species, 'Learners', p{1}{1});
%!   assert_equal (Mdl.BinaryLoss, p{1}{2});
%! endfor

%!test  # MATLAB parity: a learner keeping no data gives a compact model
%! load fisheriris
%! assert_equal (class (fitcecoc (meas, species, 'Learners', 'linear')), ...
%!               'CompactClassificationECOC');
%! assert_equal (class (fitcecoc (meas, species, 'Learners', 'kernel')), ...
%!               'CompactClassificationECOC');

%!test  # a template carries the options its learner is fitted with
%! load fisheriris
%! Mdl = fitcecoc (meas, species, 'Learners', templateTree ('MaxNumSplits', 2));
%! assert_equal (Mdl.BinaryLearners{1}.ModelParameters.MaxSplits, 2);

%!test  # MATLAB parity: the coding design is taken by name
%! load fisheriris
%! Mdl = fitcecoc (meas, species, 'Coding', 'onevsall');
%! assert_equal (Mdl.CodingName, 'onevsall');
%! assert_equal (Mdl.CodingMatrix, 2 * eye (3) - 1);
%! assert_equal (numel (Mdl.BinaryLearners), 3);

%!test  # MATLAB parity: a coding matrix given outright is 'custom'
%! load fisheriris
%! M = [1, 0, -1; -1, 1, 0; 0, -1, 1];
%! Mdl = fitcecoc (meas, species, 'Coding', M);
%! assert_equal (Mdl.CodingName, 'custom');
%! assert_equal (Mdl.CodingMatrix, M);

%!test  # what an observation was to each learner is its class's row
%! load fisheriris
%! Mdl = fitcecoc (meas, species);
%! assert_equal (Mdl.BinaryY(1,:), [1, 1, 0]);
%! assert_equal (Mdl.BinaryY(51,:), [-1, 0, 1]);
%! assert_equal (Mdl.BinaryY(101,:), [0, -1, -1]);

%!test  # MATLAB parity: a cross-validation option gives a partitioned model
%! load fisheriris
%! CV = fitcecoc (meas, species, 'KFold', 5);
%! assert_equal (class (CV), 'ClassificationPartitionedECOC');
%! assert_equal (CV.KFold, 5);
%! assert_equal (class (fitcecoc (meas, species, 'Holdout', 0.3)), ...
%!               'ClassificationPartitionedECOC');

## Test input validation
%!test  # MATLAB parity: a character matrix response counts a class per row
%! load fisheriris
%! Mdl = fitcecoc (meas, char (species));
%! assert_equal (size (Mdl.CodingMatrix), [3, 3]);
%! assert_equal (predict (Mdl, meas([1, 60, 120], :)), ...
%!               char ({'setosa'; 'versicolor'; 'virginica'}));

%!error<fitcecoc: specify only one of 'KFold', 'Holdout', 'Leaveout' and 'CVPartition'.> ...
%! fitcecoc (ones (8, 2), [1; 2; 1; 2; 1; 2; 1; 2], 'KFold', 2, 'Holdout', 0.3)
%!error<fitcecoc: too few input arguments.> fitcecoc (ones (4, 2))
%!error<fitcecoc: name-value arguments must be in pairs.> ...
%! fitcecoc (ones (4, 2), [1; 2; 1; 2], 'Coding')
%!error<fitcecoc: number of rows in X and Y must be equal.> ...
%! fitcecoc (ones (4, 2), [1; 2; 1])
%!error<ClassificationECOC: 'nosuch' is not a binary learner.> ...
%! fitcecoc (ones (4, 2), [1; 2; 1; 2], 'Learners', 'nosuch')
%!error<ClassificationECOC: 'FitPosterior' is not implemented> ...
%! fitcecoc (ones (4, 2), [1; 2; 1; 2], 'FitPosterior', true)
%!error<ClassificationECOC: invalid parameter name> ...
%! fitcecoc (ones (4, 2), [1; 2; 1; 2], 'NoSuch', 1)

%!test  # MATLAB parity: classes given as text are sorted
%! load fisheriris
%! k = [101:150, 1:50, 51:100];
%! Mdl = fitcecoc (meas(k,:), species(k));
%! assert_equal (Mdl.ClassNames, {'setosa'; 'versicolor'; 'virginica'});

%!test  # MATLAB parity: a given ClassNames order is kept
%! load fisheriris
%! Mdl = fitcecoc (meas(1:150,:), species(1:150), ...
%!             'ClassNames', {'virginica'; 'setosa'; 'versicolor'});
%! assert_equal (Mdl.ClassNames, {'virginica'; 'setosa'; 'versicolor'});

%!error<ClassificationECOC: not all 'ClassNames' are present in Y.> ...
%! load fisheriris
%! fitcecoc (meas, species, 'ClassNames', [3, 1, 2])

%!error<ClassificationECOC: not all 'ClassNames' are present in Y.> ...
%! load fisheriris
%! fitcecoc (meas, species, 'ClassNames', {'setosa'; 'rose'})

%!test  # MATLAB parity: GentleBoost ensembles as binary learners
%! load fisheriris
%! T = templateEnsemble ('GentleBoost', 5, templateTree ('MaxNumSplits', 1));
%! Mdl = fitcecoc (meas, species, 'Learners', T);
%! assert_equal (Mdl.BinaryLoss, 'exponential');
%! assert_equal (resubLoss (Mdl), 1/30, 1e-14);
%! [~, NegLoss] = predict (Mdl, meas([51, 120],:));
%! assert_equal (NegLoss, [-74.2065795512883, -0.017199167383083, ...
%!                         -4.030127061700748; -74.2065795512883, ...
%!                         -0.197045791294161, -0.32160454535994], 1e-12);

%!test  # MATLAB parity: LogitBoost ensembles as binary learners
%! load fisheriris
%! T = templateEnsemble ('LogitBoost', 5, templateTree ('MaxNumSplits', 1));
%! Mdl = fitcecoc (meas, species, 'Learners', T);
%! assert_equal (Mdl.BinaryLoss, 'binodeviance');
%! assert_equal (resubLoss (Mdl), 1/30, 1e-14);
%! [~, NegLoss] = predict (Mdl, meas([51, 120],:));
%! assert_equal (NegLoss, [-4.844912733252615, -0.002477430762383, ...
%!                         -1.86787026045674; -4.844912733252615, ...
%!                         -0.163235262337723, -0.365681737695726], 1e-12);

%!test  # MATLAB parity: the weights of an AdaBoostM1 binary learner
%! load fisheriris
%! T = templateEnsemble ('AdaBoostM1', 5, templateTree ('MaxNumSplits', 1));
%! Mdl = fitcecoc (meas, species, 'Learners', T);
%! assert_equal (Mdl.BinaryLoss, 'exponential');
%! assert_equal (Mdl.BinaryLearners{3}.TrainedWeights', ...
%!               [1.375767656520975, 0.993534110774411, ...
%!                0.883434373403998, 0.554364192108758, ...
%!                0.268194447432047], 1e-12);

%!test  # MATLAB parity: bagged and random subspace ensembles read posteriors
%! load fisheriris
%! Mdl = fitcecoc (meas, species, 'Learners', ...
%!                 templateEnsemble ('Bag', 3, 'tree'));
%! assert_equal (Mdl.BinaryLoss, 'quadratic');
%! assert_equal (class (Mdl.BinaryLearners{1}), 'ClassificationBaggedEnsemble');
%! Mdl = fitcecoc (meas, species, 'Learners', ...
%!                 templateEnsemble ('Subspace', 3, 'knn', ...
%!                                   'NPredToSample', 2));
%! assert_equal (Mdl.BinaryLoss, 'quadratic');

%!test  # MATLAB parity: an ensemble named as the learner is LogitBoost
%! load fisheriris
%! Mdl = fitcecoc (meas, species, 'Learners', 'ensemble');
%! assert_equal (Mdl.ModelParameters.BinaryLearners.Method, 'LogitBoost');
%! assert_equal (Mdl.ModelParameters.BinaryLearners.NLearn, 100);
%! assert_equal (Mdl.BinaryLoss, 'binodeviance');
%! assert_equal (Mdl.BinaryLearners{1}.NumTrained, 100);

%!test  # an ensemble of binary learners cross-validates
%! load fisheriris
%! T = templateEnsemble ('GentleBoost', 3, templateTree ('MaxNumSplits', 1));
%! CV = crossval (fitcecoc (meas, species, 'Learners', T), 'KFold', 3);
%! assert_equal (class (CV), 'ClassificationPartitionedECOC');

%!test  # RUSBoost binary scores are the weighted mean of tree probabilities
%! load fisheriris
%! T = templateEnsemble ('RUSBoost', 5, templateTree ('MaxNumSplits', 1));
%! Mdl = fitcecoc (meas, species, 'Learners', T);
%! assert_equal (Mdl.BinaryLoss, 'quadratic');
%! B = Mdl.BinaryLearners{1};
%! [~, s] = predict (B, meas(51,:));
%! [~, ~, PBScore] = predict (Mdl, meas(51,:));
%! assert_equal (PBScore(1), s(2) / sum (B.TrainedWeights), 1e-14);
%! assert_equal (all (PBScore >= 0 & PBScore <= 1), true);
%! assert_equal (isfinite (resubLoss (Mdl)), true);
%!error<ClassificationECOC: templates of regression type are not supported.> ...
%! load fisheriris
%! fitcecoc (meas, species, 'Learners', templateEnsemble ('LSBoost', 5, 'tree'))

## Table input
%!shared fecT
%! load fisheriris
%! fecT = table (meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!               'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! fecT.Species = categorical (species);
%! fecT.Wide = categorical (meas(:,2) > 3, [false true], {'narrow', 'wide'});

%!test  # the response is named by a column and the rest are predictors
%! Mdl = fitcecoc (fecT, 'Species');
%! assert_equal (Mdl.PredictorNames, {'SL', 'SW', 'PL', 'PW', 'Wide'});
%! assert_equal (Mdl.ResponseName, 'Species');
%! assert_equal (Mdl.CategoricalPredictors, 5);

%!test  # a model formula names the response and the predictors together
%! Mdl = fitcecoc (fecT, 'Species ~ PW + SL');
%! assert_equal (Mdl.PredictorNames, {'PW', 'SL'});
%! assert_equal (isempty (Mdl.CategoricalPredictors), true);

%!test  # the response may be given beside a table of predictors
%! Mdl = fitcecoc (fecT(:,1:4), fecT.Species);
%! assert_equal (Mdl.PredictorNames, {'SL', 'SW', 'PL', 'PW'});
%! assert_equal (Mdl.ResponseName, 'Y');

%!test  # predict takes a table, matched by name and not by position
%! Mdl = fitcecoc (fecT, 'Species');
%! a = predict (Mdl, fecT);
%! assert_equal (class (a), 'categorical');
%! assert_equal (predict (Mdl, fecT(:, [6, 5, 4, 3, 2, 1])), a);

%!test  # the levels travel with the model when it is made compact
%! Mdl = fitcecoc (fecT, 'Species');
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.PredictorLevels, Mdl.PredictorLevels);
%! assert_equal (predict (CMdl, fecT), predict (Mdl, fecT));

%!error<ClassificationECOC: the table holds no variable 'NoSuch'.> ...
%! fitcecoc (fecT, 'NoSuch')

%!error<ClassificationECOC: a model formula holds main effects only, so no products, powers or wildcards.> ...
%! fitcecoc (fecT, 'Species ~ SL*PW')

%!error<ClassificationECOC.predict: the table holds no predictor 'SW'.> ...
%! predict (fitcecoc (fecT, 'Species'), fecT(:, [1, 3, 4, 5, 6]))
