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
## @deftypefn  {statistics} {@var{Mdl} =} fitcnb (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcnb (@dots{}, @var{name}, @var{value})
##
## Fit a naive Bayes classification model.
##
## @code{@var{Mdl} = fitcnb (@var{X}, @var{Y})} returns a naive Bayes
## classification model, @var{Mdl}, with @var{X} being the predictor data and
## @var{Y} the class labels of the observations in @var{X}.
##
## @itemize
## @item
## @var{X} must be a @math{N*P} numeric matrix of predictor data where rows
## correspond to observations and columns correspond to features or variables.
## @item
## @var{Y} is an @math{N*1} matrix or cell matrix containing the class labels
## of the corresponding predictor data in @var{X}.  @var{Y} can be numeric,
## logical, a character array or a cell array of character vectors.  @var{Y}
## must have the same number of rows as @var{X}.
## @end itemize
##
## A naive Bayes model fits one univariate density to each predictor within
## each class, and treats the predictors as conditionally independent given the
## class.  An observation's likelihood under a class is therefore the product
## of its per-predictor densities, and its posterior follows by Bayes' rule
## from the class prior.
##
## @code{@var{Mdl} = fitcnb (@dots{}, @var{name}, @var{value})} returns a naive
## Bayes model with additional options specified by @qcode{Name-Value} pair
## arguments listed below.
##
## @subheading Model Parameters
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'PredictorNames'} @tab A cell array of character vectors
## specifying the names of the predictors.  The length of this array must match
## the number of columns in @var{X}.
##
## @item @qcode{'ResponseName'} @tab A character vector specifying the name of
## the response variable.
##
## @item @qcode{'ClassNames'} @tab Names of the classes in the class labels,
## @var{Y}, used for fitting the model.  @qcode{ClassNames} are of the same
## type as the class labels in @var{Y}.  Naming a subset of the classes keeps
## only the observations belonging to them.
##
## @item @qcode{'Prior'} @tab A numeric vector specifying the prior probability
## of each class, in the order of @qcode{ClassNames}, or the character vector
## @qcode{'empirical'} (default) to take the class frequencies, or
## @qcode{'uniform'} to give every class the same probability.
##
## @item @qcode{'Cost'} @tab A square numeric matrix of misclassification
## costs, where @code{Cost(i,j)} is the cost of classifying an observation of
## class @math{i} into class @math{j}.  The default is one off the diagonal and
## zero on it.
##
## @item @qcode{'ScoreTransform'} @tab A character vector naming a transform
## applied to the posterior returned by @code{predict}, or a function handle
## taking and returning a matrix of the same size.  The default is
## @qcode{'none'}.
##
## @item @qcode{'DistributionNames'} @tab A character vector naming the
## distribution fitted to every predictor, or a cell array of character vectors
## naming one per predictor.  Supported are @qcode{'normal'} (default),
## @qcode{'kernel'}, @qcode{'mvmn'} for a categorical predictor, and
## @qcode{'mn'} for token counts.  @qcode{'mn'} describes the whole predictor
## vector at once and so cannot be named for only some predictors.
##
## @item @qcode{'Kernel'} @tab The smoothing kernel of the predictors fitted
## with a kernel density, one of @qcode{'normal'} (default), @qcode{'box'},
## @qcode{'epanechnikov'} or @qcode{'triangle'}, given once for every predictor
## or once per predictor.
##
## @item @qcode{'Support'} @tab The support of the kernel densities, either
## @qcode{'unbounded'} (default), @qcode{'positive'}, or a two element numeric
## vector giving finite bounds.
##
## @item @qcode{'Width'} @tab The bandwidth of the kernel densities, given as a
## scalar, as one value per predictor, as one per class, or as a matrix of one
## per class and predictor.  By default each density chooses its own.
##
## @end multitable
##
## A predictor that takes one value throughout a class has no normal density
## to fit, and that combination of class and predictor is refused rather than
## answered.  Only the combination is refused, not the model: giving that
## predictor a @qcode{'kernel'} or a @qcode{'mvmn'} distribution fits the same
## data, and leaves the other predictors normal.
##
## @seealso{ClassificationNaiveBayes}
## @end deftypefn

function Mdl = fitcnb (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitcnb: too few arguments.");
  endif

  if (mod (nargin, 2) != 0)
    error ("fitcnb: name-value arguments must be in pairs.");
  endif

  ## Check predictor data and labels have equal rows
  if (rows (X) != rows (Y))
    error ("fitcnb: number of rows in X and Y must be equal.");
  endif

  ## Parse arguments to class def function
  Mdl = ClassificationNaiveBayes (X, Y, varargin{:});

endfunction

%!demo
%! ## Fit a naive Bayes classifier to Fisher's iris data and see how often it
%! ## classifies a training observation into its own species.
%!
%! load fisheriris
%! Mdl = fitcnb (meas, species)
%! printf ("resubstitution loss: %g\n", resubLoss (Mdl));

%!demo
%! ## The petal measurements separate the species far better than the sepal
%! ## ones, and a kernel density follows a skewed predictor where a normal
%! ## one cannot.
%!
%! load fisheriris
%! normalMdl = fitcnb (meas, species);
%! kernelMdl = fitcnb (meas, species, 'DistributionNames', 'kernel');
%! printf ("normal  : %g\n", resubLoss (normalMdl));
%! printf ("kernel  : %g\n", resubLoss (kernelMdl));

## Tests
%!test  # the driver returns what the constructor returns
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! assert_equal (class (Mdl), 'ClassificationNaiveBayes');
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (Mdl.ClassNames, unique (species));

%!test  # MATLAB parity: the fitted parameters and the resubstitution loss
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! assert_equal (Mdl.DistributionParameters{1,1}, ...
%!               [5.005999999999998; 0.352489687213451], 1e-13);
%! assert_equal (resubLoss (Mdl), 0.04, 1e-14);

%!test  # name-value arguments reach the constructor
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'Prior', 'uniform', 'ResponseName', 'flower');
%! assert_equal (Mdl.Prior, [1/3, 1/3, 1/3], 1e-15);
%! assert_equal (Mdl.ResponseName, 'flower');

## Test input validation
%!error<fitcnb: too few arguments.> fitcnb ()
%!error<fitcnb: too few arguments.> fitcnb (ones (4, 1))
%!error<fitcnb: name-value arguments must be in pairs.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], ones (4, 1), 'Prior')
%!error<fitcnb: number of rows in X and Y must be equal.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], ones (3, 1))
%!error<fitcnb: number of rows in X and Y must be equal.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], ones (3, 1), 'Prior', 'uniform')
%!error<ClassificationNaiveBayes: a normal distribution cannot be fit for the combination of class 2 and predictor x1. The data has zero variance.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 10, 20], [1; 1; 1; 2])
