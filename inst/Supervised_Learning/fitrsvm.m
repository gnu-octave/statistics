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
## @deftypefn  {statistics} {@var{Mdl} =} fitrsvm (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitrsvm (@dots{}, @var{name}, @var{value})
##
## Fit a support vector machine regression model.
##
## @code{@var{Mdl} = fitrsvm (@var{X}, @var{Y})} returns a support vector
## regression model, @var{Mdl}, with @var{X} being the predictor data and
## @var{Y} the continuous response of the observations in @var{X}.
##
## @itemize
## @item
## @var{X} must be an @math{NxP} numeric matrix of predictor data, where rows
## correspond to observations and columns to features or variables.
## @item
## @var{Y} must be an @math{Nx1} numeric vector holding the response of the
## corresponding predictor data in @var{X}.  @var{Y} must have the same number
## of rows as @var{X}.
## @end itemize
##
## The model is fitted by @math{epsilon}-insensitive regression: an error
## smaller than @qcode{Epsilon} costs nothing, so only the observations
## outside that tube become support vectors.  Use @code{fitcsvm} where the
## response names a class rather than a quantity.
##
## @code{@var{Mdl} = fitrsvm (@dots{}, @var{name}, @var{value})} returns a
## model with additional options specified by @qcode{Name-Value} pair
## arguments listed below.
##
## @subheading Model Parameters
##
## @multitable @columnfractions 0.32 0.68
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Standardize'} @tab A logical scalar indicating whether the
## data in @var{X} should be centred and scaled before training.  The same
## transformation is applied by @code{predict}.  The default is @qcode{false}.
##
## @item @qcode{'PredictorNames'} @tab A cell array of character vectors
## specifying the predictor variable names, in the order they appear in
## @var{X}.
##
## @item @qcode{'ResponseName'} @tab A character vector specifying the name of
## the response variable.  The default is @qcode{'Y'}.
##
## @item @qcode{'ResponseTransform'} @tab A character vector naming one of
## @qcode{'none'}, @qcode{'identity'}, @qcode{'exp'} or @qcode{'log'}, or a
## function handle of one argument, applied to the predicted response.  The
## default is @qcode{'none'}.
##
## @item @qcode{'Epsilon'} @tab A non-negative scalar, the half-width of the
## insensitive tube.  The default is @code{iqr (@var{Y}) / 13.49}, a robust
## estimate of a tenth of the response's standard deviation, which is what
## MATLAB uses; where that is zero it falls back to @math{0.1}.
##
## @item @qcode{'BoxConstraint'} @tab A positive scalar bounding the dual
## coefficients, the cost of an error outside the tube.  The default is 1.
##
## @item @qcode{'KernelFunction'} @tab A character vector naming the kernel,
## one of @qcode{'linear'}, the default, @qcode{'rbf'}, @qcode{'gaussian'},
## @qcode{'polynomial'} or @qcode{'sigmoid'}.
##
## @item @qcode{'PolynomialOrder'} @tab A positive integer, the order of the
## polynomial kernel.  The default is 3.  It is ignored by every other kernel.
##
## @item @qcode{'KernelScale'} @tab A positive scalar dividing the predictors
## before the kernel is applied.  The default is 1.
##
## @item @qcode{'KernelOffset'} @tab A non-negative scalar added to the kernel
## value.  The default is 0.
##
## @item @qcode{'SVMtype'} @tab A character vector selecting the formulation,
## either @qcode{'eps_svr'}, the default, or @qcode{'nu_svr'}.  MATLAB fits
## only the @math{epsilon} form; @qcode{'nu_svr'} is an Octave extension.
##
## @item @qcode{'Nu'} @tab A scalar in @math{(0, 1]} used by
## @qcode{'nu_svr'}, bounding the fraction of support vectors.  The default
## is 0.5.
##
## @item @qcode{'CacheSize'} @tab A positive scalar, the kernel cache in
## megabytes.  The default is 1000.
##
## @item @qcode{'Tolerance'} @tab A non-negative scalar, the tolerance of the
## termination criterion.  The default is @math{1e-6}.
##
## @item @qcode{'Shrinking'} @tab Either 0 or 1, whether to use the shrinking
## heuristic.  The default is 1.
## @end multitable
##
## @seealso{RegressionSVM, fitcsvm, fitrnet, svmtrain, svmpredict}
## @end deftypefn

function obj = fitrsvm (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitrsvm: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitrsvm: Name-Value arguments must be in pairs.");
  endif

  ## Check predictor data and response have equal rows
  if (rows (X) != rows (Y))
    error ("fitrsvm: number of rows in X and Y must be equal.");
  endif

  ## Parse arguments to classdef constructor
  obj = RegressionSVM (X, Y, varargin{:});

endfunction

%!demo
%! ## 1. Predict fuel economy from engine power and weight
%!
%! load carsmall
%! X = [Horsepower, Weight];
%! Mdl = fitrsvm (X, MPG, 'Standardize', true);
%!
%! ## Rows carrying a missing value were dropped, so ask about the ones used
%! used = Mdl.RowsUsed;
%! yFit = predict (Mdl, X(used,:));
%! plot (MPG(used), yFit, 'o', [5, 45], [5, 45], 'k-');
%! axis equal;
%! xlabel ('Observed MPG');
%! ylabel ('Predicted MPG');
%! title (sprintf ('Linear SVR, RMSE %.2f', sqrt (resubLoss (Mdl))));

%!demo
%! ## 2. The insensitive tube decides who becomes a support vector
%!
%! ## Errors smaller than Epsilon cost nothing, so a wider tube is fitted by
%! ## fewer observations and a narrower one by almost all of them.
%! load carsmall
%! X = [Horsepower, Weight];
%! eps_ = [0.1, 0.5, 1, 2, 4, 8];
%! nsv = zeros (size (eps_));
%! for k = 1:numel (eps_)
%!   m = fitrsvm (X, MPG, 'Standardize', true, 'Epsilon', eps_(k));
%!   nsv(k) = sum (m.IsSupportVector);
%! endfor
%! plot (eps_, nsv, 'o-', 'linewidth', 1.5);
%! xlabel ('Epsilon');
%! ylabel ('Number of support vectors');
%! title ('A wider tube needs fewer support vectors');

%!demo
%! ## 3. A radial kernel fits a curve a linear one cannot
%!
%! rng (42);
%! x = linspace (-3, 3, 120)';
%! y = sin (x) + randn (120, 1) * 0.1;
%! lin = fitrsvm (x, y);
%! rbf = fitrsvm (x, y, 'KernelFunction', 'rbf', 'BoxConstraint', 10);
%! plot (x, y, 'o', 'markersize', 4);
%! hold on;
%! plot (x, predict (lin, x), 'k--', 'linewidth', 1.5);
%! plot (x, predict (rbf, x), 'r-', 'linewidth', 2);
%! hold off;
%! legend ({'data', 'linear kernel', 'rbf kernel'});
%! title ('Support vector regression');

%!demo
%! ## 4. With a linear kernel the model is a plain linear function
%!
%! load carsmall
%! X = [Horsepower, Weight];
%! Mdl = fitrsvm (X, MPG, 'Standardize', true);
%! used = Mdl.RowsUsed;
%! Xs = (X(used,:) - Mdl.Mu) ./ Mdl.Sigma;
%! printf ('max |X*Beta + Bias - predict| = %g\n', ...
%!         max (abs (Xs * Mdl.Beta + Mdl.Bias - resubPredict (Mdl))));

## Test constructor
%!test
%! load carsmall
%! X = [Horsepower, Weight];
%! Mdl = fitrsvm (X, MPG, 'Standardize', true);
%! assert_equal (class (Mdl), 'RegressionSVM');
%! assert_equal (Mdl.NumPredictors, 2);
%! assert_equal (Mdl.ModelParameters.KernelFunction, 'linear');
%! assert_equal (Mdl.ModelParameters.SVMtype, 'eps_svr');
%! assert_equal (Mdl.Epsilon, iqr (MPG(Mdl.RowsUsed)) / 13.49, 1e-12);

## The driver passes its Name-Value pairs straight through
%!test
%! X = [linspace(0, 1, 30)', linspace(1, 2, 30)'];
%! Y = 3 * X(:,1) + 1;
%! Mdl = fitrsvm (X, Y, 'KernelFunction', 'rbf', 'BoxConstraint', 5, ...
%!                'Epsilon', 0.05, 'ResponseName', 'speed');
%! assert_equal (Mdl.ModelParameters.KernelFunction, 'rbf');
%! assert_equal (Mdl.ModelParameters.BoxConstraint, 5);
%! assert_equal (Mdl.Epsilon, 0.05);
%! assert_equal (Mdl.ResponseName, 'speed');
%! assert_equal (isempty (Mdl.Beta), true);

## Test input validation
%!error<fitrsvm: too few arguments.> fitrsvm ()
%!error<fitrsvm: too few arguments.> fitrsvm (ones (4, 1))
%!error<fitrsvm: Name-Value arguments must be in pairs.> ...
%! fitrsvm (ones (4, 2), ones (4, 1), 'KernelFunction')
%!error<fitrsvm: number of rows in X and Y must be equal.> ...
%! fitrsvm (ones (4, 2), ones (3, 1))
%!error<fitrsvm: number of rows in X and Y must be equal.> ...
%! fitrsvm (ones (4, 2), ones (3, 1), 'Epsilon', 1)
