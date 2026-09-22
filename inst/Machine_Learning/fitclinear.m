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
## @deftypefn  {statistics} {@var{Mdl} =} fitclinear (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitclinear (@var{Tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{Mdl} =} fitclinear (@var{Tbl}, @var{formula})
## @deftypefnx {statistics} {@var{Mdl} =} fitclinear (@var{Tbl}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitclinear (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{Mdl}, @var{FitInfo}] =} fitclinear (@dots{})
##
## Fit a linear binary classifier.
##
## @code{@var{Mdl} = fitclinear (@var{X}, @var{Y})} returns a
## @qcode{ClassificationLinear} object fitted to the predictor data @var{X}
## and the two class response @var{Y}, where @var{X} is an @math{NxP}
## numeric matrix and @var{Y} has as many rows as @var{X}.
##
## @code{@var{Mdl} = fitclinear (@dots{}, @var{name}, @var{value})} passes
## the given @qcode{Name-Value} pairs to the model.  They are documented
## under @code{ClassificationLinear}, and the ones most often wanted are
## @qcode{'Learner'}, @qcode{'Regularization'}, @qcode{'Lambda'},
## @qcode{'Solver'} and @qcode{'ObservationsIn'}.
##
## @code{[@var{Mdl}, @var{FitInfo}] = fitclinear (@dots{})} also returns a
## structure describing the optimization: what it converged to, how far it
## got, and which tolerance stopped it.  Its fields follow the solver, so a
## dual fit reports the dual variables and a mini-batch fit the batch it
## stopped on.
##
## @code{@var{Mdl} = fitclinear (@dots{}, @var{cvopt}, @var{value})} returns a
## @code{ClassificationPartitionedLinear}
## instead when one of @qcode{'CrossVal'}, @qcode{'KFold'},
## @qcode{'Holdout'}, @qcode{'Leaveout'} and @qcode{'CVPartition'} is
## given.  A cross-validated model describes no single fit, so
## @var{FitInfo} is not available beside it.
##
## @seealso{ClassificationLinear, ClassificationKernel, fitckernel}
## @end deftypefn

function [Mdl, FitInfo] = fitclinear (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitclinear: too few input arguments.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("fitclinear: name-value arguments must be in pairs.");
  endif

  ## A cross-validation option asks for a partitioned model rather than a
  ## fitted one, and the two are different classes with different methods.
  ## MATLAB refuses a second output there, having no single fit to describe,
  ## and so does this.
  ## 'CrossVal' is the one that carries a value saying whether to cross
  ## validate at all; the other four ask for it by being present.
  cvNames = {'kfold', 'holdout', 'leaveout', 'cvpartition'};
  crossval = false;
  for k = 1:2:numel (varargin)
    if (! ischar (varargin{k}))
      continue;
    endif
    if (any (strcmpi (varargin{k}, cvNames)))
      crossval = true;
    elseif (strcmpi (varargin{k}, 'crossval'))
      val = varargin{k+1};
      if (! (ischar (val) && any (strcmpi (val, {'on', 'off'}))))
        error ("%s: 'CrossVal' must be either 'on' or 'off'.", 'fitclinear');
      endif
      crossval = crossval || strcmpi (val, 'on');
    endif
  endfor

  if (crossval)
    if (nargout > 1)
      error (strcat ("fitclinear: a cross validated model has no", ...
                     " FitInfo to return; ask for the model alone."));
    endif
    Mdl = ClassificationPartitionedLinear (X, Y, varargin{:});
    return;
  endif

  ## 'CrossVal', 'off' has said its piece and the learner does not take it.
  keep = true (1, numel (varargin));
  for k = 1:2:numel (varargin)
    if (ischar (varargin{k}) && strcmpi (varargin{k}, 'crossval'))
      keep(k:k+1) = false;
    endif
  endfor
  varargin = varargin(keep);

  Mdl = ClassificationLinear (X, Y, varargin{:});

  if (nargout > 1)
    FitInfo = fitInfo_ (Mdl);
  endif

endfunction

%!demo
%! ## Fit a linear classifier to the two overlapping iris species and read
%! ## what the optimization did.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! [Mdl, FitInfo] = fitclinear (X, Y, 'Learner', 'logistic')

%!demo
%! ## Fit from a table, and predict on one
%!
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! X = meas(inds,:);
%! T = table (X(:,1), X(:,2), X(:,3), X(:,4), ...
%!            'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! T.Species = species(inds);
%!
%! ## The response is named by its column, and everything else is a predictor
%! Mdl = fitclinear (T, 'Species');
%! Mdl.PredictorNames
%! Mdl.ResponseName
%!
%! ## predict reads a table by the names the model was fitted on, so the
%! ## columns may come in any order
%! label = predict (Mdl, T(1:5, [5, 4, 3, 2, 1]));
%! label'

%!demo
%! ## Score a table
%!
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! X = meas(inds,:);
%! T = table (X(:,1), X(:,2), X(:,3), X(:,4), ...
%!            'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! T.Species = categorical (species(inds));
%! Mdl = fitclinear (T, 'Species');
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

%!test
%! ## The driver returns what the class constructor returns
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! M1 = fitclinear (X, Y);
%! M2 = ClassificationLinear (X, Y);
%! assert_equal (class (M1), 'ClassificationLinear');
%! assert_equal (M1.Beta, M2.Beta);
%! assert_equal (M1.Bias, M2.Bias);

%!test
%! ## The options reach the model
%! load fisheriris
%! Mdl = fitclinear (meas(51:end,:), species(51:end), ...
%!                   'Learner', 'logistic', 'Lambda', 0.05, ...
%!                   'ResponseName', 'species');
%! assert_equal (Mdl.Learner, 'logistic');
%! assert_equal (Mdl.Lambda, 0.05);
%! assert_equal (Mdl.ResponseName, 'species');

%!test
%! ## The second output describes the optimization
%! load fisheriris
%! [~, FitInfo] = fitclinear (meas(51:end,:), species(51:end));
%! assert_equal (FitInfo.Lambda, 0.01);
%! assert_equal (FitInfo.Solver, {'bfgs'});
%! assert_equal (FitInfo.BetaTolerance, 1e-4);
%! assert_equal (FitInfo.GradientTolerance, 1e-6);
%! assert_equal (isfield (FitInfo, 'TerminationStatus'), true);

%!test
%! ## A dual fit reports the passes it took and the dual variables it left
%! load fisheriris
%! [~, FitInfo] = fitclinear (meas(51:end,:), species(51:end), ...
%!                            'Solver', 'dual', 'PassLimit', 20);
%! assert_equal (isfield (FitInfo, 'Alpha'), true);
%! assert_equal (isfield (FitInfo, 'NumPasses'), true);
%! assert_equal (isfield (FitInfo, 'IterationLimit'), false);
%! assert_equal (FitInfo.GradientTolerance, 0);
%! assert_equal (isnan (FitInfo.GradientNorm), true);

%!test
%! ## A mini-batch fit reports the batch it stopped on and the rate it used
%! load fisheriris
%! [~, FitInfo] = fitclinear (meas(51:end,:), species(51:end), ...
%!                            'Solver', 'sgd', 'PassLimit', 5);
%! assert_equal (isfield (FitInfo, 'BatchIndex'), true);
%! assert_equal (isfield (FitInfo, 'OptimalLearnRate'), true);

%!test
%! ## A cross-validation option returns a partitioned model instead
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = fitclinear (X, Y, 'KFold', 5);
%! assert_equal (class (CVMdl), 'ClassificationPartitionedLinear');
%! assert_equal (CVMdl.KFold, 5);
%! assert_equal (numel (CVMdl.Trained), 5);

%!test
%! ## 'CrossVal' on gives the ten folds it defaults to, and 'off' the model
%! ## itself
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = fitclinear (X, Y, 'CrossVal', 'on');
%! assert_equal (CVMdl.KFold, 10);
%! assert_equal (class (fitclinear (X, Y, 'CrossVal', 'off')), ...
%!                     'ClassificationLinear');

%!error<fitclinear: a cross validated model has no FitInfo to return; ask for the model alone.> ...
%! [Mdl, FitInfo] = fitclinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'KFold', 3);

## Test input validation
%!error<fitclinear: too few input arguments.> fitclinear (ones (5, 2))
%!error<fitclinear: name-value arguments must be in pairs.> ...
%! fitclinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Learner')
%!error<ClassificationLinear: 'Learner' must be either 'svm' or 'logistic'.> ...
%! fitclinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Learner', 'tree')

%!test  # MATLAB parity: classes given as text are sorted
%! load fisheriris
%! k = [101:150, 51:100];
%! Mdl = fitclinear (meas(k,:), species(k));
%! assert_equal (Mdl.ClassNames, {'versicolor'; 'virginica'});

%!test  # MATLAB parity: a given ClassNames order is kept
%! load fisheriris
%! Mdl = fitclinear (meas(51:150,:), species(51:150), ...
%!             'ClassNames', {'virginica'; 'versicolor'});
%! assert_equal (Mdl.ClassNames, {'virginica'; 'versicolor'});

%!error<ClassificationLinear: not all 'ClassNames' are present in Y.> ...
%! load fisheriris
%! fitclinear (meas(51:150,:), species(51:150), 'ClassNames', [3, 2])

%!error<ClassificationLinear: not all 'ClassNames' are present in Y.> ...
%! load fisheriris
%! fitclinear (meas(51:150,:), species(51:150), ...
%!             'ClassNames', {'virginica'; 'rose'})

## Table input
%!shared fclT
%! load fisheriris
%! fclI = ! strcmp (species, 'setosa');
%! fclX = meas(fclI,:);
%! fclT = table (fclX(:,1), fclX(:,2), fclX(:,3), fclX(:,4), ...
%!               'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! fclT.Species = species(fclI);

%!test  # the response is named by a column and the rest are predictors
%! Mdl = fitclinear (fclT, 'Species');
%! assert_equal (Mdl.PredictorNames, {'SL', 'SW', 'PL', 'PW'});
%! assert_equal (Mdl.ResponseName, 'Species');

%!test  # a model formula names the response and the predictors together
%! Mdl = fitclinear (fclT, 'Species ~ PW + SL');
%! assert_equal (Mdl.PredictorNames, {'PW', 'SL'});

%!test  # predict takes a table, matched by name and not by position
%! Mdl = fitclinear (fclT, 'Species');
%! a = predict (Mdl, fclT);
%! assert_equal (predict (Mdl, fclT(:, [5, 4, 3, 2, 1])), a);

%!test  # a cross-validated fit takes a table too
%! Mdl = fitclinear (fclT, 'Species', 'KFold', 3);
%! assert_equal (class (Mdl), 'ClassificationPartitionedLinear');
%! assert_equal (Mdl.KFold, 3);

%!error<ClassificationLinear: the table holds no variable 'NoSuch'.> ...
%! fitclinear (fclT, 'NoSuch')

%!error<ClassificationLinear: a model formula holds main effects only, so no products, powers or wildcards.> ...
%! fitclinear (fclT, 'Species ~ SL*PW')
