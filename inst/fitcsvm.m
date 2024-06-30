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
## @deftypefn  {statistics} {@var{Mdl} =} fitcsvm (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcsvm (@dots{}, @var{name}, @var{value})
##
## Fit a Support Vector Machine classification model.
##
## @code{@var{Mdl} = fitcsvm (@var{X}, @var{Y})} returns a Support Vector Machine
## classification model, @var{Mdl}, with @var{X} being the predictor data,
## and @var{Y} the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the SVM model.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain numerical type
## of data. @var{Y} must have same numbers of Rows as @var{X}.
##
## @end itemize
##
## @code{@var{Mdl} = fitcsvm (@dots{}, @var{name}, @var{value})} returns a
## Support Vector Machine model with additional options specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.05 0.4 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"SVMtype"} @tab Specifies the type of SVM. It accepts the
## following options: (Default is 'C_SVC')
##
## @itemize
##
## @item 'C_SVC': C-SVC is the standard SVM formulation for classification tasks.
## It aims to find the optimal hyperplane that separates different classes by
## maximizing the margin between them while allowing some misclassifications.
## The parameter C controls the trade-off between maximizing the margin and
## minimizing the classification error.
##
## @item 'nu_SVC': ν-SVC is a variation of the standard SVM that introduces
## a parameter ν (nu) as an upper bound on the fraction of margin errors and
## a lower bound on the fraction of support vectors. This formulation provides
## more control over the number of support vectors and the margin errors,
## making it useful for specific classification scenarios.
##
## @item 'one_class_SVM': One-Class SVM is used for anomaly detection and
## novelty detection tasks. It aims to separate the data points of a single
## class from the origin in a high-dimensional feature space. This method is
## particularly useful for identifying outliers or unusual patterns in the data.
## For one-class SVM, @var{Y} has no effect and can be any number.
##
## @end itemize
##
## @item @tab @qcode{"KernelFunction"} @tab Specifies the method for computing
## elements of the Gram matrix. It accepts the following options:
## (default is 'gaussian' or 'rbf')
##
## @itemize
##
## @item 'linear': @math{u'*v} Computes the linear kernel, which is simply the
## dot product of the input vectors.
##
## @item 'polynomial': @math{(gamma*u'*v + coef0)^degree} Computes the
## polynomial kernel, which raises the dot product of the input vectors to a
## specified power.
##
## @item 'rbf': @math{exp(-gamma*|u-v|^2)} Computes the radial basis function
## (RBF) kernel. It measures the similarity between two vectors in a
## high-dimensional space.
##
## @item 'sigmoid': @math{tanh(gamma*u'*v + coef0)} Computes the sigmoid kernel,
## which is inspired by the activation function used in neural networks. The
## sigmoid kernel maps the input vectors into a hyperbolic tangent function
## space.
##
## @item 'precomputed': You can also specify precomputed kernel
## (kernel values in training_set_file)
##
## @end itemize
##
## @item @tab @qcode{"PolynomialOrder"} @tab A positive integer that specifies
## the order of polynomial in kernel function. The default value is 3.
##
## @item @tab @qcode{"Gamma"} @tab Specifies the gamma in kernel function.
## The default value is @math{gamma = 1/(number of features)}.
##
## @item @tab @qcode{"KernelOffset"} @tab A nonnegative scalar that specifies
## the @math{coef0} in kernel function. For the polynomial kernel, it influences
## the polynomial's shift, and for the sigmoid kernel, it affects the hyperbolic
## tangent's shift. The default value is 0.
##
## @item @tab @qcode{"BoxConstraint"} @tab A positive scalar that specifies the
## upper bound of Lagrange multipliers ie C in [0,C] i.e, the parameter C of
## class i to weight*C, for C-SVC. It determines the trade-off between
## maximizing the margin and minimizing the classification error. The
## default value of BoxConstraint is 1.
##
## @item @tab @qcode{"Nu"} @tab A positive scalar, in the range (0,1] that
## specifies the parameter nu of nu-SVC, one-class SVM. The default value is
## 0.5.
##
## @item @tab @qcode{"CacheSize"} @tab A positive scalar that specifies the
## cache size. The default value is 100.
##
## @item @tab @qcode{"Tolerance"} @tab A nonnegative scalar that specifies
## the tolerance of termination criterion. The default value is 1e-3.
##
## @item @tab @qcode{"Shrinking"} @tab Specifies whether to use shrinking
## heuristics. It accepts either 0 or 1. The default value is 1.
##
## @item @tab @qcode{"ProbabilityEstimates"} @tab Specifies whether to train
## the model for probability estimates. It accepts either 0 or 1. The default
## value is 0.
##
## @item @tab @qcode{"Weight"} @tab Provided as a structure with two fields:
## Class labels and the corresponding weight which is a scalar which specifies
## the parameter C of class i to weight*C, for C-SVC. The default value is
## 1 for all classes.
##
## @item @tab @qcode{"PredictorNames"} @tab A cell array of character vectors
## specifying the predictor variable names.  The variable names are assumed to
## be in the same order as they appear in the training data @var{X}.
##
## @item @tab @qcode{"ResponseName"} @tab A character vector specifying the name
## of the response variable.
##
## @item @tab @qcode{"KFold"} @tab A positive integer greater than 1 which
## specifies the value of k (number of folds). The dataset is divided into k
## equal-sized subsets (folds). The model is trained k times, each time using
## k-1 subsets for training and the remaining subset for validation. This
## process helps in assessing the model's generalization capability by averaging
## the cross validation accuracy across all folds. The default value is 10.
##
## @end multitable
##
## @seealso{ClassificationSVM, svmtrain, svmpredict}
## @end deftypefn

function Mdl = fitcsvm (X, Y, varargin)

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
  Mdl = ClassificationSVM (X, Y, varargin{:});

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
%! fitcsvm (ones (4,2), ones (4, 1), 'KFold')
%!error<fitcsvm: number of rows in X and Y must be equal.>
%! fitcsvm (ones (4,2), ones (3, 1))
%!error<fitcsvm: number of rows in X and Y must be equal.>
%! fitcsvm (ones (4,2), ones (3, 1), 'KFold', 2)
