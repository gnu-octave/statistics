## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{obj} =} fitrgam (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} fitrgam (@var{X}, @var{Y}, @var{name}, @var{value})
##
## Fit a Generalized Additive Model (GAM) for regression.
##
## @code{@var{obj} = fitrgam (@var{X}, @var{Y})} returns an object of
## class RegressionGAM, with matrix @var{X} containing the predictor data and
## vector @var{Y} containing the continuous response data.
##
## @itemize
## @item
## @var{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the GAM model.
## @item
## @var{Y} must be @math{Nx1} numeric vector containing the response data
## corresponding to the predictor data in @var{X}. @var{Y} must have same
## number of rows as @var{X}.
## @end itemize
##
## @code{@var{obj} = fitrgam (@dots{}, @var{name}, @var{value})} returns
## an object of class RegressionGAM with additional properties specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"predictors"} @tab Predictor Variable names, specified as
## a row vector cell of strings with the same length as the columns in @var{X}.
## If omitted, the program will generate default variable names
## @qcode{(x1, x2, ..., xn)} for each column in @var{X}.
##
## @item @tab @qcode{"responsename"} @tab Response Variable Name, specified as
## a string.  If omitted, the default value is @qcode{"Y"}.
##
## @item @tab @qcode{"formula"} @tab a model specification given as a string in
## the form @qcode{"Y ~ terms"} where @qcode{Y} represents the reponse variable
## and @qcode{terms} the predictor variables.  The formula can be used to
## specify a subset of variables for training model.  For example:
## @qcode{"Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3"} specifies four linear terms
## for the first four columns of for predictor data, and @qcode{x1:x2} and
## @qcode{x2:x3} specify the two interaction terms for 1st-2nd and 3rd-4th
## columns respectively.  Only these terms will be used for training the model,
## but @var{X} must have at least as many columns as referenced in the formula.
## If Predictor Variable names have been defined, then the terms in the formula
## must reference to those.  When @qcode{"formula"} is specified, all terms used
## for training the model are referenced in the @qcode{IntMatrix} field of the
## @var{obj} class object as a matrix containing the column indexes for each
## term including both the predictors and the interactions used.
##
## @item @tab @qcode{"interactions"} @tab a logical matrix, a positive integer
## scalar, or the string @qcode{"all"} for defining the interactions between
## predictor variables.  When given a logical matrix, it must have the same
## number of columns as @var{X} and each row corresponds to a different
## interaction term combining the predictors indexed as @qcode{true}.  Each
## interaction term is appended as a column vector after the available predictor
## column in @var{X}.  When @qcode{"all"} is defined, then all possible
## combinations of interactions are appended in @var{X} before training.  At the
## moment, parsing a positive integer has the same effect as the @qcode{"all"}
## option.  When @qcode{"interactions"} is specified, only the interaction terms
## appended to @var{X} are referenced in the @qcode{IntMatrix} field of the
## @var{obj} class object.
##
## @item @tab @qcode{"knots"} @tab a scalar or a row vector with the same
## columns as @var{X}.  It defines the knots for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 5, hence expanded to @qcode{ones (1, columns (X)) * 5}.  You can
## parse a row vector with different number of knots for each predictor
## variable to be fitted with, although not recommended.
##
## @item @tab @qcode{"order"} @tab a scalar or a row vector with the same
## columns as @var{X}.  It defines the order of the polynomial when training the
## GAM.  As a scalar, it is expanded to a row vector.  The default values is 3,
## hence expanded to @qcode{ones (1, columns (X)) * 3}.  You can parse a row
## vector with different number of polynomial order for each predictor variable
## to be fitted with, although not recommended.
##
## @item @tab @qcode{"dof"} @tab a scalar or a row vector with the same columns
## as @var{X}.  It defines the degrees of freedom for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 8, hence expanded to @qcode{ones (1, columns (X)) * 8}.  You can
## parse a row vector with different degrees of freedom for each predictor
## variable to be fitted with, although not recommended.
##
## @item @tab @qcode{"tol"} @tab a positive scalar to set the tolerance for
## covergence during training. By defaul, it is set to @qcode{1e-3}.
## @end multitable
##
## You can parse either a @qcode{"formula"} or an @qcode{"interactions"}
## optional parameter.  Parsing both parameters will result an error.
## Accordingly, you can only pass up to two parameters among @qcode{"knots"},
## @qcode{"order"}, and @qcode{"dof"} to define the required polynomial for
## training the GAM model.
##
## @seealso{RegressionGAM, regress, regress_gp}
## @end deftypefn

function obj = fitrgam (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitrgam: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitrgam: Name-Value arguments must be in pairs.");
  endif

  ## Check predictor data and labels have equal rows
  if (rows (X) != rows (Y))
    error ("fitrgam: number of rows in X and Y must be equal.");
  endif
  ## Parse arguments to class def function
  obj = RegressionGAM (X, Y, varargin{:});

endfunction

%!demo
%! # Train a RegressionGAM Model for synthetic values
%!
%! f1 = @(x) cos (3 *x);
%! f2 = @(x) x .^ 3;
%!
%! # generate x1 and x2 for f1 and f2
%! x1 = 2 * rand (50, 1) - 1;
%! x2 = 2 * rand (50, 1) - 1;
%!
%! # calculate y
%! y = f1(x1) + f2(x2);
%!
%! # add noise
%! y = y + y .* 0.2 .* rand (50,1);
%! X = [x1, x2];
%!
%! # create an object
%! a = fitrgam (X, y, "tol", 1e-3)


## Test constructor
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [1; 2; 3; 4];
%! a = fitrgam (x, y);
%! assert ({a.X, a.Y}, {x, y})
%! assert ({a.BaseModel.Intercept}, {2.5000})
%! assert ({a.Knots, a.Order, a.DoF}, {[5, 5, 5], [3, 3, 3], [8, 8, 8]})
%! assert ({a.NumObservations, a.NumPredictors}, {4, 3})
%! assert ({a.ResponseName, a.PredictorNames}, {"Y", {"x1", "x2", "x3"}})
%! assert ({a.Formula}, {[]})
%!test
%! x = [1, 2, 3, 4; 4, 5, 6, 7; 7, 8, 9, 1; 3, 2, 1, 2];
%! y = [1; 2; 3; 4];
%! pnames = {"A", "B", "C", "D"};
%! formula = "Y ~ A + B + C + D + A:C";
%! intMat = logical ([1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1;1,0,1,0]);
%! a = fitrgam (x, y, "predictors", pnames, "formula", formula);
%! assert ({a.IntMatrix}, {intMat})
%! assert ({a.ResponseName, a.PredictorNames}, {"Y", pnames})
%! assert ({a.Formula}, {formula})

## Test input validation
%!error<fitrgam: too few arguments.> fitrgam ()
%!error<fitrgam: too few arguments.> fitrgam (ones(10,2))
%!error<fitrgam: Name-Value arguments must be in pairs.>
%! fitrgam (ones (4,2), ones (4, 1), "K")
%!error<fitrgam: number of rows in X and Y must be equal.>
%! fitrgam (ones (4,2), ones (3, 1))
%!error<fitrgam: number of rows in X and Y must be equal.>
%! fitrgam (ones (4,2), ones (3, 1), "K", 2)
