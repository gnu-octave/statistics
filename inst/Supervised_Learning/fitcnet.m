## Copyright (C) 2024 Pallav Purbia <pallavpurbia@gmail.com>
## Copyright (C) 2024-2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{Mdl} =} fitcnet (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcnet (@dots{}, @var{name}, @var{value})
##
## Fit a Neural Network classification model.
##
## @code{@var{Mdl} = fitcnet (@var{X}, @var{Y})} returns a Neural Network
## classification model, @var{Mdl}, with @var{X} being the predictor data, and
## @var{Y} the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{N*P} numeric matrix of predictor data where rows
## correspond to observations and columns correspond to features or variables.
## @item
## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain any type of
## categorical data. @var{Y} must have same numbers of rows as @var{X}.
## @end itemize
##
## @code{@var{Mdl} = fitcnet (@dots{}, @var{name}, @var{value})} returns a
## Neural Network classification model with additional options specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @subheading Model Parameters
##
## @multitable @columnfractions 0.32 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Standardize'} @tab A boolean flag indicating whether
## the data in @var{X} should be standardized prior to training.
##
## @item @qcode{'PredictorNames'} @tab A cell array of character vectors
## specifying the predictor variable names.  The variable names are assumed to
## be in the same order as they appear in the training data @var{X}.
##
## @item @qcode{'ResponseName'} @tab A character vector specifying the name
## of the response variable.
##
## @item @qcode{'ClassNames'} @tab Names of the classes in the class
## labels, @var{Y}, used for fitting the Neural Network model.
## @qcode{ClassNames} are of the same type as the class labels in @var{Y}.
##
## @item @qcode{'Prior'} @tab A numeric vector specifying the prior
## probabilities for each class.  The order of the elements in @qcode{Prior}
## corresponds to the order of the classes in @qcode{ClassNames}.
##
## @item @qcode{'LayerSizes'} @tab A vector of positive integers that
## defines the sizes of the fully connected layers in the neural network model.
## Each element in LayerSizes corresponds to the number of outputs for the
## respective fully connected layer in the neural network model.
## The default value is 10.
##
## @item @qcode{'LearningRate'} @tab A positive scalar value that defines
## the learning rate during the gradient descent.  Default value is 0.003.
## A larger rate can drive every unit of a hidden layer negative, after which
## a rectifier passes no gradient and the network stops training.
## Applies only when @qcode{'Solver'} is @qcode{'sgd'}.
##
## @item @qcode{'Solver'} @tab A character vector naming the solver that
## trains the network, either @qcode{'lbfgs'} or @qcode{'sgd'}.  The
## default is @qcode{'lbfgs'}, which minimizes the loss over the whole
## training set at once by limited-memory BFGS, as MATLAB does.  It takes
## no learning rate, stops on the three tolerances below, and reaches a
## lower training loss in fewer passes over the data, though each of its
## iterations costs several passes where an epoch costs one.
## @qcode{'sgd'} visits the samples one at a time and steps down the
## gradient of each, running for @qcode{'IterationLimit'} epochs; it was
## the default before version 1.9.0.
##
## @item @qcode{'GradientTolerance'} @tab A nonnegative scalar.  Training
## stops once the gradient's infinity norm falls to or below it, which is
## the quantity MATLAB tests too.  The default is @qcode{1e-6}.  Applies
## only when @qcode{'Solver'} is @qcode{'lbfgs'}.
##
## @item @qcode{'StepTolerance'} @tab A nonnegative scalar.  Training
## stops once the step's infinity norm falls to or below it, which is the
## quantity MATLAB tests too.  The default is @qcode{1e-6}.  Applies only
## when @qcode{'Solver'} is @qcode{'lbfgs'}.
##
## @item @qcode{'LossTolerance'} @tab A real scalar.  Training stops once
## the training loss falls to or below it.  The test is on the loss
## itself and not on its change, matching MATLAB; pass @code{-Inf} to
## switch it off.  The default is @qcode{1e-6}.  Applies only when
## @qcode{'Solver'} is @qcode{'lbfgs'}.
##
## @item @qcode{'Activations'} @tab A character vector or a cellstr vector
## specifying the activation functions for the hidden layers of the neural
## network (excluding the output layer).  The available activation functions
## are @qcode{'linear'}, @qcode{'sigmoid'}, @qcode{'relu'}, @qcode{'tanh'},
## @qcode{'softmax'}, @qcode{'lrelu'}, @qcode{'prelu'}, @qcode{'elu'},
## @qcode{'gelu'}, and @qcode{'none'}.  The default value is @qcode{'relu'}.
##
## @item @qcode{'OutputLayerActivation'} @tab A character vector specifying
## the activation function for the output layer of the neural network.  The
## available activation functions are the same as for @qcode{'Activations'}.
## The default value is @qcode{'softmax'}, which makes the returned scores a
## probability over the classes and trains the network against cross entropy;
## any other value trains it against the mean squared error.
##
## @item @qcode{'IterationLimit'} @tab A positive integer scalar that
## specifies the maximum number of training iterations.  The default value is
## 1000.
## Under @qcode{'sgd'} this counts epochs, under
## @qcode{'lbfgs'} solver iterations.
##
## @item @qcode{'DisplayInfo'} @tab A boolean flag indicating whether to
## print information during training.  Default is @qcode{false}.
##
## @item @qcode{'ScoreTransform'} @tab A character vector defining one of
## the following functions or a user defined function handle, which is used
## for transforming the prediction scores returned by the @code{predict} and
## @code{resubPredict} methods.  Default value is @qcode{'none'}.
## @end multitable
##
## @multitable @columnfractions 0.3 0.75
## @headitem @var{Value} @tab @var{Description}
## @item @qcode{'doublelogit'} @tab @math{1 ./ (1 + exp (-2 * x))}
## @item @qcode{'invlogit'} @tab @math{log (x ./ (1 - x))}
## @item @qcode{'ismax'} @tab Sets the score for the class with the largest
## score to 1, and sets the scores for all other classes to 0
## @item @qcode{'logit'} @tab @math{1 ./ (1 + exp (-x))}
## @item @qcode{'none'} @tab @math{x} (no transformation)
## @item @qcode{'identity'} @tab @math{x} (no transformation)
## @item @qcode{'sign'} @tab @math{-1 for x < 0, 0 for x = 0, 1 for x > 0}
## @item @qcode{'symmetric'} @tab @math{2 * x - 1}
## @item @qcode{'symmetricismax'} @tab Sets the score for the class with
## the largest score to 1, and sets the scores for all other classes to -1
## @item @qcode{'symmetriclogit'} @tab @math{2 ./ (1 + exp (-x)) - 1}
## @end multitable
##
##
## The weights of each layer are drawn from a uniform range whose half-width
## is set by that layer's activation, and the scheme cannot be chosen: a
## rectifying activation (@qcode{'relu'}, @qcode{'lrelu'}, @qcode{'prelu'},
## @qcode{'elu'}, @qcode{'gelu'}) takes the He range
## @math{sqrt (6 / fan_in)}, because it passes only half of its input, and
## the remaining activations take the Glorot range
## @math{sqrt (6 / (fan_in + fan_out))}, which accounts for the backward pass
## as well.  A network whose layers do not share an activation is therefore
## built with both schemes.  What each layer was given is reported by the
## @qcode{LayerWeightsInitializers} field of the fitted model's
## @qcode{ModelParameters}.
##
## @seealso{ClassificationNeuralNetwork}
## @end deftypefn

function obj = fitcnet (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitcnet: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitcnet: Name-Value arguments must be in pairs.");
  endif

  ## Check predictor data and labels have equal rows
  if (rows (X) != rows (Y))
    error ("fitcnet: number of rows in X and Y must be equal.");
  endif

  ## Parse arguments to classdef constructor
  obj = ClassificationNeuralNetwork (X, Y, varargin{:});

endfunction

%!demo
%! ## 1. Train a network on Fisher's iris data and see what it got right
%!
%! load fisheriris
%! Mdl = fitcnet (meas, species);
%! pred_species = resubPredict (Mdl);
%! confusionchart (species, pred_species, 'Title', ...
%!                 'Neural network classification of Fisher''s iris data');

%!demo
%! ## 2. Watching the fit converge
%!
%! load fisheriris
%! Mdl = fitcnet (meas, species, 'IterationLimit', 400);
%!
%! ## TrainingHistory records what the solver converges on.  The default
%! ## solver is lbfgs, so that is the loss and the gradient norm; under
%! ## 'sgd' it is the loss and the accuracy instead.
%! h = Mdl.TrainingHistory;
%! plotyy (h.Iteration, h.TrainingLoss, h.Iteration, h.Gradient);
%! xlabel ('Iteration');
%! title ('Training loss, left, and gradient norm, right');

%!demo
%! ## 3. Rectified hidden layers train faster than sigmoid ones
%!
%! load fisheriris
%! iters = [5, 10, 25, 50, 100, 200, 400];
%! L = zeros (2, numel (iters));
%! for k = 1:numel (iters)
%!   for a = 1:2
%!     act = {'relu', 'sigmoid'}{a};
%!     m = fitcnet (meas, species, 'Activations', act, ...
%!                  'IterationLimit', iters(k));
%!     L(a,k) = loss (m, meas, species, 'LossFun', 'classiferror');
%!   endfor
%! endfor
%! semilogx (iters, L(1,:), 'o-', iters, L(2,:), 's-', 'linewidth', 1.5);
%! xlabel ('Iteration limit');
%! ylabel ('Misclassification rate');
%! legend ({'relu', 'sigmoid'});
%! title ('A sigmoid shrinks the gradient at every layer');

%!demo
%! ## 4. What the network learned, over two predictors
%!
%! load fisheriris
%! X = meas(:,3:4);
%! Mdl = fitcnet (X, species, 'LayerSizes', [12, 12], 'IterationLimit', 400);
%!
%! ## Ask about a grid and paint each point by the answer
%! [gx, gy] = meshgrid (linspace (0.5, 7.5, 120), linspace (0, 3, 120));
%! [~, ~, region] = unique (predict (Mdl, [gx(:), gy(:)]));
%! contourf (gx, gy, reshape (region, size (gx)), [1 2 3]);
%! colormap (summer);
%! hold on;
%! gscatter (X(:,1), X(:,2), species, 'krb', 'ox+');
%! hold off;
%! xlabel ('Petal length');
%! ylabel ('Petal width');
%! title ('Decision regions of a two-layer network');

## Test constructor
%!test
%! load fisheriris
%! x = meas;
%! y = grp2idx (species);
%! Mdl = fitcnet (x, y, 'IterationLimit', 50);
%! assert_equal (class (Mdl), "ClassificationNeuralNetwork");
%! assert_equal (numel (Mdl.LayerWeights), 2);
%! assert_equal (size (Mdl.LayerWeights{1}), [10, 4]);
%! assert_equal (size (Mdl.LayerBiases{1}), [10, 1]);
%! assert_equal (size (Mdl.LayerWeights{2}), [3, 10]);
%! assert_equal (size (Mdl.LayerBiases{2}), [3, 1]);

## Test input validation
%!error<fitcnet: too few arguments.> fitcnet ()
%!error<fitcnet: too few arguments.> fitcnet (ones (4,1))
%!error<fitcnet: Name-Value arguments must be in pairs.>
%! fitcnet (ones (4,2), ones (4, 1), 'LayerSizes')
%!error<fitcnet: number of rows in X and Y must be equal.>
%! fitcnet (ones (4,2), ones (3, 1))
%!error<fitcnet: number of rows in X and Y must be equal.>
%! fitcnet (ones (4,2), ones (3, 1), 'LayerSizes', 2)
