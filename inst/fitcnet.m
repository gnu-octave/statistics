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
## @deftypefn  {statistics} {@var{obj} =} fitcnet (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} fitcnet (@dots{}, @var{name}, @var{value})
##
## Fit a Neural Network classification model.
##
## @code{@var{obj} = fitcnet (@var{X}, @var{Y})} returns a Neural Network
## classification model, @var{obj}, with @var{X} being the predictor data,
## and @var{Y} the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the Classification Neural Network model.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain any type of
## categorical data. @var{Y} must have same numbers of Rows as @var{X}.
## @end itemize
##
## @code{@var{obj} = fitcnet (@dots{}, @var{name}, @var{value})} returns a
## Neural Network model with additional options specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.05 0.4 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"LayerSizes"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds).
##
## @item @tab @qcode{"Activations"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds).
##
## @item @tab @qcode{"LayerWeightsInitializer"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds).
##
## @item @tab @qcode{"LayerBiasesInitializer"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds).
##
## @item @tab @qcode{"InitialStepSize"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds).
##
## @item @tab @qcode{"IterationLimit"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds).
##
## @item @tab @qcode{"GradientTolerance"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds).
##
## @item @tab @qcode{"LossTolerance"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds).
##
## @item @tab @qcode{"StepTolerance"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds).
##
## @item @tab @qcode{"Weights"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds).
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
  ## Parse arguments to class def function
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
