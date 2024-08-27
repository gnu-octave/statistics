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
## @deftypefn  {statistics} {@var{Mdl} =} fitcdiscr (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcdiscr (@dots{}, @var{name}, @var{value})
##
## Fit a Linear Discriminant Analysis classification model.
##
## @code{@var{Mdl} = fitcdiscr (@var{X}, @var{Y})} returns a Linear Discriminant
## Analysis (LDA) classification model, @var{Mdl}, with @var{X} being the
## predictor data, and @var{Y} the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of predictor data where rows
## correspond to observations and columns correspond to features or variables.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can be numerical, logical,
## char array or cell array of character vectors. @var{Y} must have same number
## of rows as @var{X}.
## @end itemize
##
## @code{@var{Mdl} = fitcdiscr (@dots{}, @var{name}, @var{value})} returns a
## Linear Discriminant Analysis model with additional options specified by
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
## @item @qcode{"Prior"} @tab @tab A numeric vector specifying the prior
## probabilities for each class.  The order of the elements in @qcode{Prior}
## corresponds to the order of the classes in @qcode{ClassNames}.
## Alternatively, you can specify @qcode{"empirical"} to use the empirical
## class probabilities or @qcode{"uniform"} to assume equal class probabilities.
##
## @item @qcode{"Cost"} @tab @tab A @math{NxR} numeric matrix containing
## misclassification cost for the corresponding instances in @var{X} where
## @math{R} is the number of unique categories in @var{Y}.  If an instance is
## correctly classified into its category the cost is calculated to be 1,
## otherwise 0. cost matrix can be altered use @code{@var{Mdl.cost} = somecost}.
## default value @qcode{@var{cost} = ones(rows(X),numel(unique(Y)))}.
##
## @item @qcode{"DiscrimType"} @tab @tab A character vector or string scalar
## specifying the type of discriminant analysis to perform. The only supported
## value is @qcode{"linear"}.
##
## @item @qcode{"FillCoeffs"} @tab @tab A character vector or string scalar
## with values @qcode{"on"} or @qcode{"off"} specifying whether to fill the
## coefficients after fitting. If set to @qcode{"on"}, the coefficients are
## computed during model fitting, which can be useful for prediction.
##
## @item @qcode{"Gamma"} @tab @tab A numeric scalar specifying the
## regularization parameter for the covariance matrix. It adjusts the linear
## discriminant analysis to make the model more stable in the presence of
## multicollinearity or small sample sizes. A value of 0 corresponds to no
## regularization, while a value of 1 corresponds to
## a completely regularized model.
##
## @end multitable
## @seealso{ClassificationDiscriminant}
## @end deftypefn

function obj = fitcdiscr (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitcdiscr: too few arguments.");
  endif

  if (mod (nargin, 2) != 0)
    error ("fitcdiscr: name-value arguments must be in pairs.");
  endif

  ## Check predictor data and labels have equal rows
  if (rows (X) != rows (Y))
    error ("fitcdiscr: number of rows in X and Y must be equal.");
  endif

  ## Parse arguments to class def function
  obj = ClassificationDiscriminant (X, Y, varargin{:});

endfunction

%!demo
%! ## Train a linear discriminant classifier for Gamma = 0.5
%! ## and plot the decision boundaries.
%!
%! load fisheriris
%! idx = ! strcmp (species, "setosa");
%! X = meas(idx,3:4);
%! Y = cast (strcmpi (species(idx), "virginica"), "double");
%! obj = fitcdiscr (X, Y, "Gamma", 0.5)
%! x1 = [min(X(:,1)):0.03:max(X(:,1))];
%! x2 = [min(X(:,2)):0.02:max(X(:,2))];
%! [x1G, x2G] = meshgrid (x1, x2);
%! XGrid = [x1G(:), x2G(:)];
%! pred = predict (obj, XGrid);
%! gidx = logical (str2num (cell2mat (pred)));
%!
%! figure
%! scatter (XGrid(gidx,1), XGrid(gidx,2), "markerfacecolor", "magenta");
%! hold on
%! scatter (XGrid(!gidx,1), XGrid(!gidx,2), "markerfacecolor", "red");
%! plot (X(Y == 0, 1), X(Y == 0, 2), "ko", X(Y == 1, 1), X(Y == 1, 2), "kx");
%! xlabel ("Petal length (cm)");
%! ylabel ("Petal width (cm)");
%! title ("Linear Discriminant Analysis Decision Boundary");
%! legend ({"Versicolor Region", "Virginica Region", ...
%!         "Sampled Versicolor", "Sampled Virginica"}, ...
%!         "location", "northwest")
%! axis tight
%! hold off

## Tests
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, "Gamma", 0.5);
%! [label, score, cost] = predict (Mdl, [2, 2, 2, 2]);
%! assert (label, {'versicolor'})
%! assert (score, [0, 0.9999, 0.0001], 1e-4)
%! assert (cost, [1, 0.0001, 0.9999], 1e-4)
%! [label, score, cost] = predict (Mdl, [2.5, 2.5, 2.5, 2.5]);
%! assert (label, {'versicolor'})
%! assert (score, [0, 0.6368, 0.3632], 1e-4)
%! assert (cost, [1, 0.3632, 0.6368], 1e-4)
%! assert (class (Mdl), "ClassificationDiscriminant");
%! assert ({Mdl.X, Mdl.Y, Mdl.NumObservations}, {meas, species, 150})
%! assert ({Mdl.DiscrimType, Mdl.ResponseName}, {"linear", "Y"})
%! assert ({Mdl.Gamma, Mdl.MinGamma}, {0.5, 0})
%! assert (Mdl.ClassNames, unique (species))
%! sigma = [0.265008, 0.046361, 0.083757, 0.019201; ...
%!          0.046361, 0.115388, 0.027622, 0.016355; ...
%!          0.083757, 0.027622, 0.185188, 0.021333; ...
%!          0.019201, 0.016355, 0.021333, 0.041882];
%! assert (Mdl.Sigma, sigma, 1e-6)
%! mu = [5.0060, 3.4280, 1.4620, 0.2460; ...
%!       5.9360, 2.7700, 4.2600, 1.3260; ...
%!       6.5880, 2.9740, 5.5520, 2.0260];
%! assert (Mdl.Mu, mu, 1e-14)
%! assert (Mdl.LogDetSigma, -8.6884, 1e-4)

## Test input validation
%!error<fitcdiscr: too few arguments.> fitcdiscr ()
%!error<fitcdiscr: too few arguments.> fitcdiscr (ones (4,1))
%!error<fitcdiscr: name-value arguments must be in pairs.>
%! fitcdiscr (ones (4,2), ones (4, 1), "K")
%!error<fitcdiscr: number of rows in X and Y must be equal.>
%! fitcdiscr (ones (4,2), ones (3, 1))
%!error<fitcdiscr: number of rows in X and Y must be equal.>
%! fitcdiscr (ones (4,2), ones (3, 1), "K", 2)
