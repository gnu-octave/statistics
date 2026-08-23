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
## @deftypefn  {statistics} {@var{Mdl} =} fitrnet (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitrnet (@dots{}, @var{name}, @var{value})
##
## Fit a neural network regression model.
##
## @code{@var{Mdl} = fitrnet (@var{X}, @var{Y})} returns a neural network
## regression model, @var{Mdl}, with @var{X} being the predictor data and
## @var{Y} the continuous response of the observations in @var{X}.
##
## @itemize
## @item
## @var{X} must be an @math{NxP} numeric matrix of predictor data, where rows
## correspond to observations and columns to features or variables.
## @item
## @var{Y} must be an @math{Nx1} numeric vector holding the response of the
## corresponding predictor data in @var{X}.  @var{Y} must have the same number
## of rows as @var{X}.
## @end itemize
##
## The network is trained against the mean squared error and its output layer
## applies the identity, so a prediction is an unrestricted real number.  Use
## @code{fitcnet} where the response names a class rather than a quantity.
##
## @code{@var{Mdl} = fitrnet (@dots{}, @var{name}, @var{value})} returns a
## neural network regression model with additional options specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @subheading Model Parameters
##
## @multitable @columnfractions 0.32 0.68
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Standardize'} @tab A logical scalar indicating whether the
## data in @var{X} should be centred and scaled before training.  The same
## transformation is applied by @code{predict}.  The default is @qcode{false}.
##
## @item @qcode{'PredictorNames'} @tab A cell array of character vectors
## specifying the predictor variable names, in the order they appear in
## @var{X}.
##
## @item @qcode{'ResponseName'} @tab A character vector specifying the name of
## the response variable.  The default is @qcode{'Y'}.
##
## @item @qcode{'ResponseTransform'} @tab A character vector naming one of
## @qcode{'none'}, @qcode{'identity'}, @qcode{'exp'} or @qcode{'log'}, or a
## function handle of one argument, applied to the predicted response by
## @code{predict} and @code{resubPredict}.  The default is @qcode{'none'}.
##
## @item @qcode{'LayerSizes'} @tab A vector of positive integers defining the
## number of units in each fully connected hidden layer.  The default value is
## 10, a single hidden layer of ten units.
##
## @item @qcode{'LearningRate'} @tab A positive scalar value that defines the
## learning rate during the gradient descent.  Default value is 0.003.  A
## larger rate can drive every unit of a hidden layer negative, after which a
## rectifier passes no gradient and the network stops training.
## Applies only when @qcode{'Solver'} is @qcode{'sgd'}.
##
## @item @qcode{'Solver'} @tab A character vector naming the solver that
## trains the network, either @qcode{'sgd'} or @qcode{'lbfgs'}.  The
## default is @qcode{'sgd'}, which visits the samples one at a time and
## steps down the gradient of each, running for @qcode{'IterationLimit'}
## epochs.  @qcode{'lbfgs'} minimizes the loss over the whole training
## set at once by limited-memory BFGS, which is the solver MATLAB uses.
## It takes no learning rate, stops on the three tolerances below, and
## reaches a lower training loss in fewer passes over the data, though
## each of its iterations costs several passes where an epoch costs one.
##
## @item @qcode{'GradientTolerance'} @tab A nonnegative scalar.  Training
## stops once the gradient's infinity norm falls to or below it.  The
## default is @qcode{1e-6}.  Applies only when @qcode{'Solver'} is
## @qcode{'lbfgs'}.
##
## @item @qcode{'StepTolerance'} @tab A nonnegative scalar.  Training
## stops once the step's infinity norm falls to or below it.  The
## default is @qcode{1e-6}.  Applies only when @qcode{'Solver'} is
## @qcode{'lbfgs'}.
##
## @item @qcode{'LossTolerance'} @tab A real scalar.  Training stops once
## the training loss falls to or below it.  The test is on the loss
## itself and not on its change, matching MATLAB; pass @code{-Inf} to
## switch it off.  The default is @qcode{1e-6}.  Applies only when
## @qcode{'Solver'} is @qcode{'lbfgs'}.
##
## @item @qcode{'Activations'} @tab A character vector or a cellstr vector
## specifying the activation functions for the hidden layers of the neural
## network, excluding the output layer.  The available activation functions
## are @qcode{'linear'}, @qcode{'none'}, @qcode{'sigmoid'}, @qcode{'relu'},
## @qcode{'tanh'}, @qcode{'lrelu'}, @qcode{'prelu'}, @qcode{'elu'} and
## @qcode{'gelu'}.  The default value is @qcode{'relu'}.
##
## @item @qcode{'OutputLayerActivation'} @tab A character vector specifying
## the activation function for the output layer.  The available functions are
## the same as for @qcode{'Activations'}.  The default value is
## @qcode{'none'}, the identity, which is what a regression output calls for;
## anything else bounds the prediction to that function's range.
##
## @item @qcode{'IterationLimit'} @tab A positive integer scalar specifying
## the maximum number of training iterations.  The default value is 1000.
## Under @qcode{'sgd'} this counts epochs, under
## @qcode{'lbfgs'} solver iterations.
##
## @item @qcode{'DisplayInfo'} @tab A logical scalar indicating whether to
## print information during training.  Default is @qcode{false}.
## @end multitable
##
## @seealso{RegressionNeuralNetwork, fitcnet, fcnntrain, fcnnpredict}
## @end deftypefn

function obj = fitrnet (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitrnet: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitrnet: Name-Value arguments must be in pairs.");
  endif

  ## Check predictor data and response have equal rows
  if (rows (X) != rows (Y))
    error ("fitrnet: number of rows in X and Y must be equal.");
  endif

  ## Parse arguments to classdef constructor
  obj = RegressionNeuralNetwork (X, Y, varargin{:});

endfunction

%!demo
%! ## 1. Predict fuel economy from engine power and weight
%!
%! load carsmall
%! X = [Horsepower, Weight];
%! Mdl = fitrnet (X, MPG, 'Standardize', true, 'IterationLimit', 500);
%!
%! ## Rows carrying a missing value were dropped, so ask about the ones used
%! used = Mdl.RowsUsed;
%! yFit = predict (Mdl, X(used,:));
%! plot (MPG(used), yFit, 'o', [5, 45], [5, 45], 'k-');
%! axis equal;
%! xlabel ('Observed MPG');
%! ylabel ('Predicted MPG');
%! title (sprintf ('Neural network fit, RMSE %.2f', sqrt (resubLoss (Mdl))));

%!demo
%! ## 2. Watching the fit converge
%!
%! load carsmall
%! Mdl = fitrnet ([Horsepower, Weight], MPG, 'Standardize', true, ...
%!                'IterationLimit', 400);
%!
%! ## TrainingHistory records the mean squared error at every iteration
%! h = Mdl.TrainingHistory;
%! semilogy (h.Iteration, h.TrainingLoss, 'linewidth', 1.5);
%! xlabel ('Iteration');
%! ylabel ('Training MSE');
%! title ('The loss recorded is the network''s own, not a running average');

%!demo
%! ## 3. Standardizing matters when the predictors differ in scale
%!
%! ## Horsepower runs to a few hundred and Weight to a few thousand, so the
%! ## heavier column dominates the first layer until both are put on one scale.
%! load carsmall
%! X = [Horsepower, Weight];
%! raw = fitrnet (X, MPG, 'IterationLimit', 400);
%! std_ = fitrnet (X, MPG, 'Standardize', true, 'IterationLimit', 400);
%! printf ('RMSE, raw predictors : %.2f\n', sqrt (resubLoss (raw)));
%! printf ('RMSE, standardized   : %.2f\n', sqrt (resubLoss (std_)));

%!demo
%! ## 4. A network recovers a curve a straight line cannot
%!
%! rand ('seed', 7);
%! randn ('seed', 7);
%! x = linspace (-3, 3, 120)';
%! y = sin (x) + randn (120, 1) * 0.1;
%! Mdl = fitrnet (x, y, 'LayerSizes', [16, 16], 'IterationLimit', 800);
%! plot (x, y, 'o', 'markersize', 4);
%! hold on;
%! plot (x, predict (Mdl, x), 'r-', 'linewidth', 2);
%! plot (x, [ones(120,1), x] * ([ones(120,1), x] \ y), 'k--', 'linewidth', 1.5);
%! hold off;
%! legend ({'data', 'neural network', 'least squares line'});
%! title ('Two hidden layers of sixteen units');

## Test constructor
%!test
%! rand ('seed', 42);
%! X = linspace (-1, 1, 40)';
%! Y = 2 * X + 0.5;
%! Mdl = fitrnet (X, Y, 'IterationLimit', 50);
%! assert_equal (class (Mdl), 'RegressionNeuralNetwork');
%! assert_equal (numel (Mdl.ModelParameters.LayerWeights), 2);
%! assert_equal (size (Mdl.ModelParameters.LayerWeights{1}), [10, 2]);
%! assert_equal (size (Mdl.ModelParameters.LayerWeights{2}), [1, 11]);

## The driver passes its Name-Value pairs straight through
%!test
%! rand ('seed', 42);
%! X = linspace (0, 1, 30)';
%! Mdl = fitrnet (X, 3 * X, 'LayerSizes', [5, 5], 'LearningRate', 0.01, ...
%!                'IterationLimit', 40, 'ResponseName', 'speed');
%! assert_equal (Mdl.LayerSizes, [5, 5]);
%! assert_equal (Mdl.LearningRate, 0.01);
%! assert_equal (Mdl.IterationLimit, 40);
%! assert_equal (Mdl.ResponseName, 'speed');

## A model on real data with missing values trains on the rows that remain
%!test
%! rand ('seed', 42);
%! load carsmall
%! X = [Horsepower, Weight];
%! Mdl = fitrnet (X, MPG, 'Standardize', true, 'IterationLimit', 200);
%! keep = ! isnan (MPG);
%! assert_equal (Mdl.NumObservations, sum (keep));
%! assert_equal (Mdl.RowsUsed, keep);
%! assert_equal (numel (resubPredict (Mdl)), sum (keep));

## Test input validation
%!error<fitrnet: too few arguments.> fitrnet ()
%!error<fitrnet: too few arguments.> fitrnet (ones (4, 1))
%!error<fitrnet: Name-Value arguments must be in pairs.> ...
%! fitrnet (ones (4, 2), ones (4, 1), 'LayerSizes')
%!error<fitrnet: number of rows in X and Y must be equal.> ...
%! fitrnet (ones (4, 2), ones (3, 1))
%!error<fitrnet: number of rows in X and Y must be equal.> ...
%! fitrnet (ones (4, 2), ones (3, 1), 'LayerSizes', 2)
