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
## @deftypefn  {statistics} {@var{ScoreSVMModel} =} fitSVMPosterior (@var{SVMModel})
## @deftypefnx {statistics} {@var{ScoreSVMModel} =} fitSVMPosterior (@var{SVMModel}, @var{name}, @var{value})
## @deftypefnx {statistics} {@var{ScoreSVMModel} =} fitSVMPosterior (@var{SVMModel}, @var{X}, @var{Y})
## @deftypefnx {statistics} {@var{ScoreSVMModel} =} fitSVMPosterior (@var{SVMModel}, @var{Tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{ScoreSVMModel} =} fitSVMPosterior (@var{SVMModel}, @var{Tbl}, @var{Y})
## @deftypefnx {statistics} {@var{ScoreSVMModel} =} fitSVMPosterior (@var{SVMModel}, @var{Tbl})
## @deftypefnx {statistics} {[@var{ScoreSVMModel}, @var{ScoreTransform}] =} fitSVMPosterior (@dots{})
##
## Fit the transformation from SVM scores to posterior probabilities.
##
## @code{@var{ScoreSVMModel} = fitSVMPosterior (@var{SVMModel})} returns the
## two-class SVM classifier @var{SVMModel} with the transformation from the
## score of the second class, @code{@var{SVMModel}.ClassNames(2)}, to its
## posterior probability stored in its @qcode{ScoreTransform} property, so
## that @code{predict}, or @code{kfoldPredict} for a cross-validated model,
## returns posterior probabilities, the two columns summing to 1.
## @var{SVMModel} is a @code{ClassificationSVM}, a
## @code{CompactClassificationSVM} or a @code{ClassificationPartitionedModel}
## cross-validating an SVM, and @var{ScoreSVMModel} is of the same class.
##
## The transformation is fitted to scores of observations the model was not
## trained on:
##
## @itemize
## @item For a @code{ClassificationSVM}, scores from cross-validating it,
## 10-fold unless one of the options below says otherwise.
## @item For a @code{ClassificationPartitionedModel}, its own out-of-fold
## scores, from the partition it was built with; no option may be given.
## @item With predictors @var{X} and their classes @var{Y}, the scores the
## model gives @var{X}, without cross-validation; a
## @code{CompactClassificationSVM} holds no training data and needs them.
## @var{X} may also be a table @var{Tbl}, matched to the predictors by
## variable name, with the classes in the variable @var{ResponseVarName}
## names, in @var{Y}, or in the variable the model was fitted on.  Rows whose
## score is @qcode{NaN} are left out.
## @end itemize
##
## @multitable @columnfractions 0.28 0.7
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'KFold'} @tab The number of folds, an integer greater
## than 1.  The default is 10.
##
## @item @qcode{'Holdout'} @tab The fraction of the data held out, a
## scalar in the range @math{(0,1)}.
##
## @item @qcode{'Leaveout'} @tab @qcode{'on'} for leave-one-out
## cross-validation, or @qcode{'off'} (default).
##
## @item @qcode{'CVPartition'} @tab A @qcode{cvpartition} object over the
## observations of @var{SVMModel}.
## @end multitable
##
## Only one of these options may be given, and only with a
## @code{ClassificationSVM} and no data.  The transformation is one of three:
##
## @itemize
## @item A sigmoid @code{1 / (1 + exp (A * s + B))} of the score @math{s}
## when the scores of the two classes overlap, with @math{A} and @math{B}
## fitted by maximum likelihood, as a logistic regression of the class on the
## score.  Neither the observation weights nor the prior enter the fit, as in
## MATLAB.
## @item A step function when the scores separate the classes: 0 below the
## largest score of the first class, 1 above the smallest score of the second,
## and @code{@var{SVMModel}.Prior(2)} between, with a warning.
## @item A constant when one class has prior probability 0, posterior 1 for
## the other class.
## @end itemize
##
## A @qcode{ScoreTransform} already set on @var{SVMModel} is reset to
## @qcode{'none'} first, with a warning.  The property holds the
## transformation as @qcode{'@@(S)sigmoid(S,A,B)'},
## @qcode{'@@(S)step(S,LowerBound,UpperBound,PositiveClassProbability)'} or
## @qcode{'@@(S)constant(S,K)'}, @math{K} being the index of the predicted
## class.  Each parameter is written exactly, in the fewest significant
## digits that read back to it, where MATLAB shows seven, so a model saved and
## reloaded predicts as before.  The folds of a cross-validated model keep
## @qcode{'none'}, the transformation being applied to the scores they give
## together.
##
## @code{[@var{ScoreSVMModel}, @var{ScoreTransform}] = fitSVMPosterior
## (@dots{})} also returns the parameters in a structure whose field
## @qcode{Type} is @qcode{'sigmoid'}, with fields @qcode{Slope} and
## @qcode{Intercept}, @qcode{'step'}, with fields @qcode{LowerBound},
## @qcode{UpperBound} and @qcode{PositiveClassProbability}, or
## @qcode{'constant'}, with the field @qcode{PredictedClass}.
##
## Reference: J. Platt.  Probabilistic outputs for support vector machines and
## comparisons to regularized likelihood methods.  In @emph{Advances in Large
## Margin Classifiers}, MIT Press, 2000.
##
## @seealso{ClassificationSVM.fitPosterior,
## CompactClassificationSVM.fitPosterior, fitcsvm, cvpartition}
## @end deftypefn

function [ScoreSVMModel, ScoreTransform] = fitSVMPosterior (SVMModel, varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ispart = isa (SVMModel, 'ClassificationPartitionedModel');
  isfull = isa (SVMModel, 'ClassificationSVM');
  if (! (isfull || isa (SVMModel, 'CompactClassificationSVM')
         || (ispart && strcmp (SVMModel.CrossValidatedModel, 'SVM'))))
    error (strcat ("fitSVMPosterior: SVMMODEL must be an SVM classifier", ...
                   " or a cross-validated SVM classifier."));
  endif
  if (classCount (SVMModel.ClassNames) != 2)
    error ("fitSVMPosterior: SVMMODEL must have exactly two classes.");
  endif

  ## Data given after the model, rather than Name-Value options
  hasdata = (numel (varargin) > 0 && (isnumeric (varargin{1})
                                      || islogical (varargin{1})
                                      || istable (varargin{1})));
  if (ispart && numel (varargin) > 0)
    error (strcat ("fitSVMPosterior: no optional argument may be given", ...
                   " for a cross-validated model."));
  endif
  if (! (isfull || ispart || hasdata))
    error (strcat ("fitSVMPosterior: a compact model needs the", ...
                   " predictors X and their classes Y."));
  endif

  if (hasdata)
    X = varargin{1};
    if (istable (X))
      if (numel (varargin) > 2)
        error ("fitSVMPosterior: too many input arguments.");
      elseif (numel (varargin) == 2 && (ischar (varargin{2})
                                        || (isstring (varargin{2})
                                            && isscalar (varargin{2}))))
        name = char (varargin{2});
        if (! any (strcmp (X.Properties.VariableNames, name)))
          error ("fitSVMPosterior: '%s' is not a variable of TBL.", name);
        endif
        Y = X.(name);
      elseif (numel (varargin) == 2)
        Y = varargin{2};
      elseif (any (strcmp (X.Properties.VariableNames, SVMModel.ResponseName)))
        Y = X.(SVMModel.ResponseName);
      else
        error (strcat ("fitSVMPosterior: TBL must hold the response,", ...
                       " or Y or ResponseVarName must be given."));
      endif
    else
      if (numel (varargin) < 2)
        error (strcat ("fitSVMPosterior: the predictors X must be given", ...
                       " with their classes Y."));
      elseif (numel (varargin) > 2)
        error ("fitSVMPosterior: too many input arguments.");
      endif
      Y = varargin{2};
    endif
    if (rows (X) != rows (Y) && rows (X) != numel (Y))
      error ("fitSVMPosterior: Y must have one class for each row of X.");
    endif
  endif

  ## Only the cross-validation options, which cvPartitionOf builds the
  ## partition from; it takes 'CrossVal' too, which is not one of them
  if (isfull && ! hasdata)
    if (any (strcmpi (varargin(1:2:end), 'CrossVal')))
      error ("fitSVMPosterior: unknown parameter name 'CrossVal'.");
    endif
    [part, rest] = cvPartitionOf (varargin, SVMModel.Y, ...
                                  SVMModel.NumObservations, "fitSVMPosterior");
    if (! isempty (rest))
      if (ischar (rest{1}))
        error ("fitSVMPosterior: unknown parameter name '%s'.", rest{1});
      endif
      error (strcat ("fitSVMPosterior: optional arguments must be given", ...
                     " in Name-Value pairs."));
    endif
  endif

  if (! strcmp (SVMModel.ScoreTransform, 'none'))
    warning (strcat ("fitSVMPosterior: resetting 'ScoreTransform' to", ...
                     " 'none' before fitting the transformation."));
    SVMModel.ScoreTransform = 'none';
  endif

  ## Scores of observations the model was not trained on, which a prior of
  ## zero makes unnecessary
  s = [];
  k = [];
  if (all (SVMModel.Prior > 0))
    if (hasdata)
      [~, scores] = predict (SVMModel, X);
      [k, errmsg] = labelIndices (SVMModel.ClassNames, Y);
      if (! isempty (errmsg))
        error ("fitSVMPosterior: %s", errmsg);
      endif
    elseif (ispart)
      [~, scores] = kfoldPredict (SVMModel);
      k = labelIndices (SVMModel.ClassNames, SVMModel.Y);
    else
      CVMdl = crossval (SVMModel, 'CVPartition', part);
      [~, scores] = kfoldPredict (CVMdl);
      k = labelIndices (SVMModel.ClassNames, SVMModel.Y);
    endif
    s = scores(:,2);
  endif

  [st, ScoreTransform, separated] = svmPosterior (s, k, SVMModel.Prior, ...
                                                  SVMModel.ClassNames);
  if (separated)
    warning (strcat ("fitSVMPosterior: the classes are perfectly", ...
                     " separated, so the transformation is a step function."));
  endif
  SVMModel.ScoreTransform = st;
  ScoreSVMModel = SVMModel;

endfunction

%!demo
%! ## Posterior probabilities for two overlapping species
%! load fisheriris
%! X = meas(51:150, 3:4);
%! Y = species(51:150);
%! SVMModel = fitcsvm (X, Y, 'KernelFunction', 'linear');
%! [ScoreSVMModel, ScoreTransform] = fitSVMPosterior (SVMModel)
%! [label, posterior] = predict (ScoreSVMModel, X([1, 2, 60, 61], :))

%!demo
%! ## Two species the scores separate get a step function
%! load fisheriris
%! SVMModel = fitcsvm (meas(1:100, 3:4), species(1:100));
%! [~, ScoreTransform] = fitSVMPosterior (SVMModel)

## Test output
## Values checked against MATLAB R2024a, whose SMO solver gives scores up to
## 1.3e-3 away from LIBSVM's out of fold
%!test
%! load fisheriris
%! Mdl = fitcsvm (meas(51:150, 3:4), species(51:150), ...
%!                'KernelFunction', 'linear');
%! [M, ST] = fitSVMPosterior (Mdl, 'Leaveout', 'on');
%! assert_equal (class (M), 'ClassificationSVM');
%! assert_equal (ST.Type, 'sigmoid');
%! assert_equal (ST.Slope, -3.4011027773863747, 1e-3);
%! assert_equal (ST.Intercept, -0.092575940016464561, 1e-4);
%! [~, ST2] = fitPosterior (Mdl, 'Leaveout', 'on');
%! assert_equal (ST, ST2);
%!test
%! load fisheriris
%! X = meas(51:150, 3:4);
%! Y = species(51:150);
%! Mdl = fitcsvm (X, Y, 'KernelFunction', 'linear');
%! [M, ST] = fitSVMPosterior (Mdl, X, Y);
%! assert_equal (class (M), 'ClassificationSVM');
%! assert_equal (ST.Slope, -3.3591642562674378, 1e-4);
%! assert_equal (ST.Intercept, -0.1511873140962528, 1e-4);
%!test
%! load fisheriris
%! X = meas(51:150, 3:4);
%! Y = species(51:150);
%! C = compact (fitcsvm (X, Y, 'KernelFunction', 'linear'));
%! [M, ST] = fitSVMPosterior (C, X, Y);
%! assert_equal (class (M), 'CompactClassificationSVM');
%! assert_equal (ST.Slope, -3.3591642562674378, 1e-4);
%! [~, ST2] = fitPosterior (C, X, Y);
%! assert_equal (ST, ST2);
%!test
%! load fisheriris
%! Mdl = fitcsvm (meas(51:150, 3:4), species(51:150), ...
%!                'KernelFunction', 'linear');
%! CVMdl = crossval (Mdl, 'Leaveout', 'on');
%! [M, ST] = fitSVMPosterior (CVMdl);
%! assert_equal (class (M), 'ClassificationPartitionedModel');
%! assert_equal (ST.Slope, -3.4011027773863747, 1e-3);
%! tok = regexp (M.ScoreTransform, '^@\(S\)sigmoid\(S,([^,]+),([^)]+)\)$', ...
%!               'tokens', 'once');
%! assert_equal (str2double (tok)', [ST.Slope, ST.Intercept]);
%! assert_equal (M.Trained{1}.ScoreTransform, 'none');
%! [~, p] = kfoldPredict (M);
%! assert_equal (p(1:3,:), [0.9326612658405814, 0.067338734159418587; ...
%!                          0.9825175955481934, 0.017482404451806603; ...
%!                          0.69049950170233787, 0.30950049829766219], 1e-3);
%!test
%! load fisheriris
%! T = table (meas(51:150, 3), meas(51:150, 4), species(51:150), ...
%!            'VariableNames', {'PL', 'PW', 'Species'});
%! C = compact (fitcsvm (T, 'Species', 'KernelFunction', 'linear'));
%! [~, ST] = fitSVMPosterior (C, meas(51:150, 3:4), species(51:150));
%! [~, ST1] = fitSVMPosterior (C, T, 'Species');
%! [~, ST2] = fitSVMPosterior (C, T);
%! [~, ST3] = fitSVMPosterior (C, T(:,1:2), species(51:150));
%! assert_equal (ST1, ST);
%! assert_equal (ST2, ST);
%! assert_equal (ST3, ST);
%!warning <fitSVMPosterior: the classes are perfectly separated, so the transformation is a step function.> ...
%! load fisheriris
%! Mdl = fitcsvm (meas(1:100, 3:4), species(1:100), 'KernelFunction', 'linear');
%! [~, ST] = fitSVMPosterior (Mdl, 'Leaveout', 'on');
%! assert_equal (ST.Type, 'step');
%! assert_equal (ST.LowerBound, -0.8446461538461536, 1e-6);
%! assert_equal (ST.UpperBound, 0.68965517241379271, 1e-6);
%! assert_equal (ST.PositiveClassProbability, 0.5);
%!warning <fitSVMPosterior: resetting 'ScoreTransform' to 'none' before fitting the transformation.> ...
%! load fisheriris
%! X = meas(51:150, 3:4);
%! Y = species(51:150);
%! Mdl = fitcsvm (X, Y, 'KernelFunction', 'linear');
%! Mdl.ScoreTransform = 'logit';
%! [~, ST] = fitSVMPosterior (Mdl, X, Y);
%! assert_equal (ST.Slope, -3.3591642562674378, 1e-4);
%!test
%! load fisheriris
%! Mdl = fitcsvm (meas(51:150, 3:4), species(51:150), 'Prior', [0, 1]);
%! [M, ST] = fitSVMPosterior (Mdl);
%! assert_equal (ST.Type, 'constant');
%! assert_equal (ST.PredictedClass, {'virginica'});
%! assert_equal (M.ScoreTransform, '@(S)constant(S,2)');

## Test input validation
%!error <Invalid call to fitSVMPosterior> fitSVMPosterior ()
%!error <fitSVMPosterior: SVMMODEL must be an SVM classifier or a cross-validated SVM classifier.> ...
%! load fisheriris
%! fitSVMPosterior (fitcknn (meas(51:150, 3:4), species(51:150)))
%!error <fitSVMPosterior: SVMMODEL must be an SVM classifier or a cross-validated SVM classifier.> ...
%! load fisheriris
%! fitSVMPosterior (crossval (fitcknn (meas(51:150, 3:4), species(51:150))))
%!error <fitSVMPosterior: SVMMODEL must have exactly two classes.> ...
%! load fisheriris
%! fitSVMPosterior (fitcsvm (meas(51:150, 3:4), ones (100, 1)))
%!error <fitSVMPosterior: no optional argument may be given for a cross-validated model.> ...
%! load fisheriris
%! CVMdl = crossval (fitcsvm (meas(51:150, 3:4), species(51:150)));
%! fitSVMPosterior (CVMdl, 'KFold', 5)
%!error <fitSVMPosterior: a compact model needs the predictors X and their classes Y.> ...
%! load fisheriris
%! fitSVMPosterior (compact (fitcsvm (meas(51:150, 3:4), species(51:150))))
%!error <fitSVMPosterior: the predictors X must be given with their classes Y.> ...
%! load fisheriris
%! C = compact (fitcsvm (meas(51:150, 3:4), species(51:150)));
%! fitSVMPosterior (C, meas(51:150, 3:4))
%!error <fitSVMPosterior: too many input arguments.> ...
%! load fisheriris
%! C = compact (fitcsvm (meas(51:150, 3:4), species(51:150)));
%! fitSVMPosterior (C, meas(51:150, 3:4), species(51:150), 'Leaveout', 'on')
%!error <fitSVMPosterior: Y must have one class for each row of X.> ...
%! load fisheriris
%! C = compact (fitcsvm (meas(51:150, 3:4), species(51:150)));
%! fitSVMPosterior (C, meas(51:150, 3:4), species(51:100))
%!error <fitSVMPosterior: Y must hold only classes the model was trained on.> ...
%! load fisheriris
%! C = compact (fitcsvm (meas(51:150, 3:4), species(51:150)));
%! fitSVMPosterior (C, meas(1:100, 3:4), species(1:100))
%!error <fitSVMPosterior: 'Kind' is not a variable of TBL.> ...
%! load fisheriris
%! T = table (meas(51:150, 3), meas(51:150, 4), species(51:150), ...
%!            'VariableNames', {'PL', 'PW', 'Species'});
%! C = compact (fitcsvm (T, 'Species'));
%! fitSVMPosterior (C, T, 'Kind')
%!error <fitSVMPosterior: TBL must hold the response, or Y or ResponseVarName must be given.> ...
%! load fisheriris
%! T = table (meas(51:150, 3), meas(51:150, 4), species(51:150), ...
%!            'VariableNames', {'PL', 'PW', 'Species'});
%! C = compact (fitcsvm (T, 'Species'));
%! fitSVMPosterior (C, T(:,1:2))
%!error <fitSVMPosterior: unknown parameter name 'Foo'.> ...
%! load fisheriris
%! fitSVMPosterior (fitcsvm (meas(51:150, 3:4), species(51:150)), 'Foo', 1)
%!error <fitSVMPosterior: unknown parameter name 'CrossVal'.> ...
%! load fisheriris
%! fitSVMPosterior (fitcsvm (meas(51:150, 3:4), species(51:150)), ...
%!                  'CrossVal', 'on')
%!error <fitSVMPosterior: you can use only one of 'KFold', 'Holdout', 'Leaveout', or 'CVPartition' options.> ...
%! load fisheriris
%! fitSVMPosterior (fitcsvm (meas(51:150, 3:4), species(51:150)), ...
%!                  'KFold', 5, 'Leaveout', 'on')
%!error <fitSVMPosterior: 'Holdout' must be a numeric value between 0 and 1.> ...
%! load fisheriris
%! fitSVMPosterior (fitcsvm (meas(51:150, 3:4), species(51:150)), ...
%!                  'Holdout', 1)
