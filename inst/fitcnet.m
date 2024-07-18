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
## @multitable @columnfractions 0.18 0.02 0.8
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
## @item @qcode{"LayerSizes"} @tab @tab A vector of positive integers that
## defines the sizes of the fully connected layers in the neural network model.
## Each element in LayerSizes corresponds to the number of outputs for the
## respective fully connected layer in the neural network model.
## The default value is 10.
##
## @item @qcode{"Activations"} @tab @tab A character vector specifying the
## activation function for the fully connected layers of the neural network
## model.  The available Activation functions are @qcode{'relu'}, which is the
## default, or @qcode{'tanh'}, @qcode{'sigmoid'}, and @qcode{'none'}.
##
## @item @qcode{"LayerWeightsInitializer"} @tab @tab A character vector
## specifying the function to initialize fully connected layer weights. The
## available Layer Weights Initializer are @qcode{'glorot'}, which is the
## default, and @qcode{'he'}.
##
## @item @qcode{"LayerBiasesInitializer"} @tab @tab A character vector
## specifying the type of initial fully connected layer biases. The available
## Layer Biases Initializer are @qcode{'zeros'}, which is the default, and
## @qcode{'ones'}.

##
## @item @qcode{"InitialStepSize"} @tab @tab A positive integer greater than 1
## which specifies the value of k (number of folds).


##
## @item @qcode{"IterationLimit"} @tab @tab A positive integer scalar that
## specifies the maximum number of training iterations. The default value is
## 1e3.
##
## @item @qcode{"GradientTolerance"} @tab @tab A nonnegative scalar that
## specifies the relative gradient tolerance as a termination criterion. The
## default value is 1e-6.
##
## @item @qcode{"LossTolerance"} @tab @tab A nonnegative scalar. If the function
## loss at some iteration is smaller than LossTolerance, then the training
## process terminates. The default value is 1e-6.
##
## @item @qcode{"StepTolerance"} @tab @tab A nonnegative scalar. If the step
## size at some iteration is smaller than StepTolerance, then the training
## process terminates. The default value is 1e-6.
##
## @item @qcode{"Weights"} @tab @tab Observation weights, specified as a
## nonnegative numeric vector. This weights each observation in X with the
## corresponding value in Weights. The length of Weights must equal the number
## of observations in X. By default, Weights is ones(n,1), where n is the number
## of observations in X.
##
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
## No demo for now.

## Test constructor
%!test
## No test for now.

## Test input validation
%!error<fitcnet: too few arguments.> fitcnet ()
%!error<fitcnet: too few arguments.> fitcnet (ones (4,1))
%!error<fitcnet: Name-Value arguments must be in pairs.>
%! fitcnet (ones (4,2), ones (4, 1), 'LayerSizes')
%!error<fitcnet: number of rows in X and Y must be equal.>
%! fitcnet (ones (4,2), ones (3, 1))
%!error<fitcnet: number of rows in X and Y must be equal.>
%! fitcnet (ones (4,2), ones (3, 1), 'LayerSizes', 2)
