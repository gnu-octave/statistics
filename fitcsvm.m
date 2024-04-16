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
## @deftypefn  {statistics} {@var{obj} =} fitcsvm (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} fitcsvm (@dots{}, @var{name}, @var{value})
##
## Fit a Support Vector Machine classification model.
##
## @code{@var{obj} = fitcsvm (@var{X}, @var{Y})} returns a Support Vector Machine
## classification model, @var{obj}, with @var{X} being the predictor data,
## and @var{Y} the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the SVM model.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain any type of
## categorical data. @var{Y} must have same numbers of Rows as @var{X}.
## @item
## @end itemize
##
## @code{@var{obj} = fitcsvm (@dots{}, @var{name}, @var{value})} returns a
## Support Vector Machine model with additional options specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @qcode{"Alpha"} @tab A vector of non negative elements used as
## initial estimates of the alpha coefficients. Each element in the vector
## corresponds to a row in the input data @var(X). The default value of Alpha is:
## The default value of Alpha is:
## - For two-class learning: zeros(size(X,1), 1)
## - For one-class learning: 0.5 * ones(size(X,1), 1)
##
## @item @qcode{"BoxConstraint"} @tab A positive scalar that specifies the
## upper bound of Lagrange multipliers ie C in [0,C]. It determines the trade-off
## between maximizing the margin and minimizing the classification error. The
## default value of BoxConstraint is 1.
##
## @item @qcode{"CacheSize"} @tab Specifies the cache size. It can be:
## @itemize
## @item A positive scalar that specifies the cache size in megabytes (MB).
## @item A string "maximal" which will result in cache large enough to hold the
## entire Gram matrix of size @math{NxN} where N is the number of rows in X.
## The default value is 1000.
## @end itemize
##
## @item @qcode{"CategoricalPredictors"} @tab @tab
##
## @item @qcode{"ClassNames"} @tab @tab
##
## @item @qcode{"ClipAlphas"} @tab @tab
##
## @item @qcode{"Cost"} @tab @tab
##
## @item @qcode{"CrossVal"} @tab @tab
##
## @item @qcode{"CVPartition"} @tab @tab
##
## @item @qcode{"Holdout"} @tab @tab
##
## @item @qcode{"KFold"} @tab @tab
##
## @item @qcode{"Leaveout"} @tab @tab
##
## @item @qcode{"GapTolerance"} @tab @tab
##
## @item @qcode{"DeltaGradientTolerance"} @tab @tab
##
## @item @qcode{"KKTTolerance"} @tab @tab
##
## @item @qcode{"IterationLimit"} @tab @tab
##
## @item @qcode{"KernelFunction"} @tab @tab Specifies the method for computing
## elements of the Gram matrix. It accepts the following options:
## @itemize
## @item 'linear': Computes the linear kernel, which is simply the dot product
## of the input vectors.
## @item 'gaussian' or 'rbf': Computes the Gaussian kernel, also known as the
## radial basis function (RBF) kernel. It measures the similarity between two
## vectors in a high-dimensional space.
## @item 'polynomial': Computes the polynomial kernel, which raises the
## dot product of the input vectors to a specified power.
## @item You can also specify the name of a custom kernel function. It must be of
## the form: function G = KernelFunc(U, V)
## This custom function must take two input matrices, U and V, and return a
## matrix G of size M-by-N, where M and N are the number of rows in U and V.
## @end itemize
##
## @item @qcode{"KernelScale"} @tab @tab
##
## @item @qcode{"KernelOffset"} @tab @tab
##
## @item @qcode{"OptimizeHyperparameters"} @tab @tab
##
## @item @qcode{"PolynomialOrder"} @tab @tab
##
## @item @qcode{"Nu"} @tab @tab
##
## @item @qcode{"NumPrint"} @tab @tab
##
## @item @qcode{"OutlierFraction"} @tab @tab
##
## @item @qcode{"PredictorNames"} @tab @tab
##
## @item @qcode{"Prior"} @tab @tab
##
## @item @qcode{"RemoveDuplicates"} @tab @tab
##
## @item @qcode{"ResponseName"} @tab @tab Response Variable Name, specified as
## a string. If omitted, the default value is @qcode{"Y"}.
##
## @item @qcode{"ScoreTransform"} @tab @tab
##
## @item @qcode{"Solver"} @tab @tab
##
## @item @qcode{"ShrinkagePeriod"} @tab @tab
##
## @item @qcode{"Standardize"} @tab @tab
##
## @item @qcode{"Verbose"} @tab @tab
##
## @item @qcode{"Weights"} @tab @tab
##
## @item @qcode{"DeltaGradientTolerance"} @tab @tab
##
## @end multitable
##
## @seealso{ClassificationSVM, svmtrain, svmpredict}
## @end deftypefn

function obj = fitcsvm (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitcsvm: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitcsvm: Name-Value arguments must be in pairs.");
  endif

  ## Check predictor data and labels have equal rows
  if (rows (X) != rows (Y))
    error ("fitcsvm: number of rows in X and Y must be equal.");
  endif
  ## Parse arguments to class def function
  obj = ClassificationSVM (X, Y, varargin{:});

endfunction

%!demo
## No demo for now.

## Test constructor
%!test
## No test for now.

## Test input validation
%!error<fitcsvm: too few arguments.> fitcsvm ()
%!error<fitcsvm: too few arguments.> fitcsvm (ones (4,1))
%!error<fitcsvm: Name-Value arguments must be in pairs.>
%! fitcsvm (ones (4,2), ones (4, 1), 'Prior')
%!error<fitcsvm: number of rows in X and Y must be equal.>
%! fitcsvm (ones (4,2), ones (3, 1))
%!error<fitcsvm: number of rows in X and Y must be equal.>
%! fitcsvm (ones (4,2), ones (3, 1), 'KFold', 2)
