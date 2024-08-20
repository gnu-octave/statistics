## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
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
## @deftypefn  {statistics} {@var{Mdl} =} fitcgam (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcgam (@dots{}, @var{name}, @var{value})
##
## Fit a Generalized Additive Model (GAM) for binary classification.
##
## @code{@var{Mdl} = fitcgam (@var{X}, @var{Y})} returns a a GAM classification
## model, @var{Mdl}, with @var{X} being the predictor data, and @var{Y} the
## binary class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of predictor data where rows
## correspond to observations and columns correspond to features or variables.
## @item
## @code{Y} is @math{Nx1} numeric vector containing binary class labels, 
## typically 0 or 1.
## @end itemize
##
## @code{@var{Mdl} = fitcgam (@dots{}, @var{name}, @var{value})} returns a
## GAM classification model with additional options specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @subheading Model Parameters
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"PredictorNames"} @tab @tab A cell array of character vectors
## specifying the names of the predictors. The length of this array must match
## the number of columns in @var{X}.
##
## @item @qcode{"ResponseName"} @tab @tab A character vector specifying the
## name of the response variable.
##
## @item @qcode{"ClassNames"} @tab @tab Names of the classes in the class
## labels, @var{Y}, used for fitting the Discriminant model. @qcode{ClassNames}
## are of the same type as the class labels in @var{Y}.
##
## @item @qcode{"Cost"} @tab @tab A @math{NxR} numeric matrix containing
## misclassification cost for the corresponding instances in @var{X} where
## @math{R} is the number of unique categories in @var{Y}.  If an instance is
## correctly classified into its category the cost is calculated to be 1,
## otherwise 0. cost matrix can be altered use @code{@var{Mdl.cost} = somecost}.
## default value @qcode{@var{cost} = ones(rows(X),numel(unique(Y)))}.
##
## @item @qcode{"Formula"} @tab @tab A model specification given as a string in
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
## @item @qcode{"Interactions"} @tab @tab A logical matrix, a positive integer
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
## @item @qcode{"Knots"} @tab @tab A scalar or a row vector with the same
## columns as @var{X}.  It defines the knots for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 5, hence expanded to @qcode{ones (1, columns (X)) * 5}.  You can
## parse a row vector with different number of knots for each predictor
## variable to be fitted with, although not recommended.
##
## @item @qcode{"Order"} @tab @tab A scalar or a row vector with the same
## columns as @var{X}.  It defines the order of the polynomial when training the
## GAM.  As a scalar, it is expanded to a row vector.  The default values is 3,
## hence expanded to @qcode{ones (1, columns (X)) * 3}.  You can parse a row
## vector with different number of polynomial order for each predictor variable
## to be fitted with, although not recommended.
##
## @item @qcode{"DoF"} @tab @tab A scalar or a row vector with the same columns
## as @var{X}.  It defines the degrees of freedom for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 8, hence expanded to @qcode{ones (1, columns (X)) * 8}.  You can
## parse a row vector with different degrees of freedom for each predictor
## variable to be fitted with, although not recommended.
##
## @end multitable
## You can parse either a @qcode{"Formula"} or an @qcode{"Interactions"}
## optional parameter.  Parsing both parameters will result an error.
## Accordingly, you can only pass up to two parameters among @qcode{"Knots"},
## @qcode{"Order"}, and @qcode{"DoF"} to define the required polynomial for
## training the GAM model.
##
## @seealso{ClassificationGAM}
## @end deftypefn

function obj = fitcgam (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitcgam: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitcgam: name-value arguments must be in pairs.");
  endif

  ## Check predictor data and labels have equal rows
  if (rows (X) != rows (Y))
    error ("fitcgam: number of rows in X and Y must be equal.");
  endif

  ## Parse arguments to class def function
  obj = ClassificationGAM (X, Y, varargin{:});

endfunction

## Demo
%!demo
%! ## Train a GAM classifier for binary classification
%! ## using specific data and plot the decision boundaries.
%!
%! ## Define specific data
%! X = [1, 2; 2, 3; 3, 3; 4, 5; 5, 5; ...
%!     6, 7; 7, 8; 8, 8; 9, 9; 10, 10];
%! Y = [0; 0; 0; 0; 0; ...
%!     1; 1; 1; 1; 1];
%!
%! ## Train the GAM model
%! obj = fitcgam (X, Y, "Interactions", "all");
%!
%! ## Create a grid of values for prediction
%! x1 = [min(X(:,1)):0.1:max(X(:,1))];
%! x2 = [min(X(:,2)):0.1:max(X(:,2))];
%! [x1G, x2G] = meshgrid (x1, x2);
%! XGrid = [x1G(:), x2G(:)];
%! pred = predict (obj, XGrid);
%!
%! ## Plot decision boundaries and data points
%! predNumeric = str2double (pred);
%! gidx = predNumeric > 0.5;
%!
%! figure
%! scatter(XGrid(gidx,1), XGrid(gidx,2), "markerfacecolor", "magenta");
%! hold on
%! scatter(XGrid(!gidx,1), XGrid(!gidx,2), "markerfacecolor", "red");
%! plot(X(Y == 0, 1), X(Y == 0, 2), "ko", X(Y == 1, 1), X(Y == 1, 2), "kx");
%! xlabel("Feature 1");
%! ylabel("Feature 2");
%! title("Generalized Additive Model (GAM) Decision Boundary");
%! legend({"Class 1 Region", "Class 0 Region", ...
%!       "Class 1 Samples", "Class 0 Samples"}, ...
%!       "location", "northwest")
%! axis tight
%! hold off

## Tests
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [0; 0; 1; 1];
%! PredictorNames = {'Feature1', 'Feature2', 'Feature3'};
%! a = fitcgam (x, y, "PredictorNames", PredictorNames);
%! assert (class (a), "ClassificationGAM");
%! assert ({a.X, a.Y, a.NumObservations}, {x, y, 4})
%! assert ({a.NumPredictors, a.ResponseName}, {3, "Y"})
%! assert (a.ClassNames, {'0'; '1'})
%! assert (a.PredictorNames, PredictorNames)
%! assert (a.BaseModel.Intercept, 0)
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10];
%! y = [1; 0; 1; 0; 1];
%! a = fitcgam (x, y, "interactions", "all");
%! assert (class (a), "ClassificationGAM");
%! assert ({a.X, a.Y, a.NumObservations}, {x, y, 5})
%! assert ({a.NumPredictors, a.ResponseName}, {2, "Y"})
%! assert (a.ClassNames, {'1'; '0'})
%! assert (a.PredictorNames, {'x1', 'x2'})
%! assert (a.ModelwInt.Intercept, 0.4055, 1e-1)
%!test
%! load fisheriris
%! inds = strcmp (species,'versicolor') | strcmp (species,'virginica');
%! X = meas(inds, :);
%! Y = species(inds, :)';
%! Y = strcmp (Y, 'virginica')';
%! a = fitcgam (X, Y, 'Formula', 'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3');
%! assert (class (a), "ClassificationGAM");
%! assert ({a.X, a.Y, a.NumObservations}, {X, Y, 100})
%! assert ({a.NumPredictors, a.ResponseName}, {4, "Y"})
%! assert (a.ClassNames, {'0'; '1'})
%! assert (a.Formula, 'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3')
%! assert (a.PredictorNames, {'x1', 'x2', 'x3', 'x4'})
%! assert (a.ModelwInt.Intercept, 0)

## Test input validation
%!error<fitcgam: too few arguments.> fitcgam ()
%!error<fitcgam: too few arguments.> fitcgam (ones (4,1))
%!error<fitcgam: name-value arguments must be in pairs.>
%! fitcgam (ones (4,2), ones (4, 1), "K")
%!error<fitcgam: number of rows in X and Y must be equal.>
%! fitcgam (ones (4,2), ones (3, 1))
%!error<fitcgam: number of rows in X and Y must be equal.>
%! fitcgam (ones (4,2), ones (3, 1), "K", 2)
