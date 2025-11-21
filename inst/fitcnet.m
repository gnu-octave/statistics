## Copyright (C) 2024 Pallav Purbia <pallavpurbia@gmail.com>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @code{X} must be a @math{NxP} numeric matrix of predictor data where rows
## correspond to observations and columns correspond to features or variables.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
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
## @multitable @columnfractions 0.32 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"Standardize"} @tab @tab A boolean flag indicating whether
## the data in @var{X} should be standardized prior to training.
##
## @item @qcode{"PredictorNames"} @tab @tab A cell array of character vectors
## specifying the predictor variable names.  The variable names are assumed to
## be in the same order as they appear in the training data @var{X}.
##
## @item @qcode{"ResponseName"} @tab @tab A character vector specifying the name
## of the response variable.
##
## @item @qcode{"ClassNames"} @tab @tab Names of the classes in the class
## labels, @var{Y}, used for fitting the Neural Network model.
## @qcode{ClassNames} are of the same type as the class labels in @var{Y}.
##
## @item @qcode{"Prior"} @tab @tab A numeric vector specifying the prior
## probabilities for each class.  The order of the elements in @qcode{Prior}
## corresponds to the order of the classes in @qcode{ClassNames}.
##
## @item @qcode{"LayerSizes"} @tab @tab A vector of positive integers that
## defines the sizes of the fully connected layers in the neural network model.
## Each element in LayerSizes corresponds to the number of outputs for the
## respective fully connected layer in the neural network model.
## The default value is 10.
##
## @item @qcode{"LearningRate"} @tab @tab A positive scalar value that defines
## the learning rate during the gradient descent.  Default value is 0.01.
##
## @item @qcode{"Activations"} @tab @tab A character vector or a cellstr vector
## specifying the activation functions for the hidden layers of the neural
## network (excluding the output layer).  The available activation functions are
## @qcode{'linear'}, @qcode{'sigmoid'}, @qcode{'tanh'}, @qcode{'sigmoid'}, and
## @qcode{'none'}. The default value is @qcode{'sigmoid'}.
##
## @item @qcode{"OutputLayerActivation"} @tab @tab A character vector specifying
## the activation function for the output layer of the neural network.  The
## available activation functions are @qcode{'linear'}, @qcode{'sigmoid'},
## @qcode{'tanh'}, @qcode{'sigmoid'}, and @qcode{'none'}. The default value is
## @qcode{'sigmoid'}.
##
## @item @qcode{"IterationLimit"} @tab @tab A positive integer scalar that
## specifies the maximum number of training iterations.  The default value is
## 1000.
##
## @item @qcode{"DisplayInfo"} @tab @tab A boolean flag indicating whether to
## print information during training.  Default is @qcode{false}.
##
## @item @qcode{"ScoreTransform"} @tab @tab A character vector defining one of
## the following functions or a user defined function handle, which is used
## for transforming the prediction scores returned by the @code{predict} and
## @code{resubPredict} methods.  Default value is @qcode{'none'}.
## @end multitable
##
## @multitable @columnfractions 0.05 0.3 0.75
## @headitem @tab @var{Value} @tab @var{Description}
## @item @tab @qcode{"doublelogit"} @tab @math{1 ./ (1 + exp (-2 * x))}
## @item @tab @qcode{"invlogit"} @tab @math{log (x ./ (1 - x))}
## @item @tab @qcode{"ismax"} @tab Sets the score for the class with the largest
## score to 1, and sets the scores for all other classes to 0
## @item @tab @qcode{"logit"} @tab @math{1 ./ (1 + exp (-x))}
## @item @tab @qcode{"none"} @tab @math{x} (no transformation)
## @item @tab @qcode{"identity"} @tab @math{x} (no transformation)
## @item @tab @qcode{"sign"} @tab @math{-1 for x < 0, 0 for x = 0, 1 for x > 0}
## @item @tab @qcode{"symmetric"} @tab @math{2 * x - 1}
## @item @tab @qcode{"symmetricismax"} @tab Sets the score for the class with
## the largest score to 1, and sets the scores for all other classes to -1
## @item @tab @qcode{"symmetriclogit"} @tab @math{2 ./ (1 + exp (-x)) - 1}
## @end multitable
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
%! ## Train a Neural Network on the Fisher's Iris data set and display
%! ## a confusion chart with the classification results.
%!
%! load fisheriris
%! Mdl = fitcnet (meas, species);
%! pred_species = resubPredict (Mdl);
%! confusionchart (species, pred_species, "Title", ...
%! "Fully Connected Neural Network classification on Fisher's Iris dataset");

## Test constructor
%!test
%! load fisheriris
%! x = meas;
%! y = grp2idx (species);
%! Mdl = fitcnet (x, y, "IterationLimit", 50);
%! assert (class (Mdl), "ClassificationNeuralNetwork");
%! assert (numel (Mdl.ModelParameters.LayerWeights), 2);
%! assert (size (Mdl.ModelParameters.LayerWeights{1}), [10, 5]);
%! assert (size (Mdl.ModelParameters.LayerWeights{2}), [3, 11]);

## Test input validation
%!error<fitcnet: too few arguments.> fitcnet ()
%!error<fitcnet: too few arguments.> fitcnet (ones (4,1))
%!error<fitcnet: Name-Value arguments must be in pairs.>
%! fitcnet (ones (4,2), ones (4, 1), 'LayerSizes')
%!error<fitcnet: number of rows in X and Y must be equal.>
%! fitcnet (ones (4,2), ones (3, 1))
%!error<fitcnet: number of rows in X and Y must be equal.>
%! fitcnet (ones (4,2), ones (3, 1), 'LayerSizes', 2)
