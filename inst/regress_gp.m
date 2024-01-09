## Copyright (c) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
## Copyright (C) 2023-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{Yfit}, @var{Yint}, @var{m}, @var{K}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit})
## @deftypefnx {statistics} {[@var{Yfit}, @var{Yint}, @var{m}, @var{K}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{"linear"})
## @deftypefnx {statistics} {[@var{Yfit}, @var{Yint}, @var{Ysd}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{"rbf"})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{"linear"}, @var{Sp})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @var{Sp})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{"rbf"}, @var{theta})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{"rbf"}, @var{theta}, @var{g})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{"rbf"}, @var{theta}, @var{g}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @var{theta})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @var{theta}, @var{g})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @var{theta}, @var{g}, @var{alpha})
##
## Regression using Gaussian Processes.
##
## @code{[@var{Yfit}, @var{Yint}, @var{m}, @var{K}] = regress_gp (@var{X},
## @var{Y}, @var{Xfit})} will estimate a linear Gaussian Process model @var{m}
## in the form @qcode{@var{Y} = @var{X}' * @var{m}}, where @var{X} is an
## @math{NxP} matrix with @math{N} observations in @math{P} dimensional space
## and @var{Y} is an @math{Nx1} column vector as the dependent variable.  The
## information about errors of the predictions (interpolation/extrapolation) is
## given by the covarianve matrix @var{K}.
## By default, the linear model defines the prior covariance of @var{m} as
## @code{@var{Sp} = 100 * eye (size (@var{X}, 2) + 1)}.  A custom prior
## covariance matrix can be passed as @var{Sp}, which must be a @math{P+1xP+1}
## positive definite matrix.  The model is evaluated for input @var{Xfit}, which
## must have the same columns as @var{X}, and the estimates are returned in
## @var{Yfit} along with the estimated variation in @var{Yint}.
## @qcode{@var{Yint}(:,1)} contains the upper boundary estimate and
## @qcode{@var{Yint}(:,1)} contains the upper boundary estimate with respect to
## @var{Yfit}.
##
## @code{[@var{Yfit}, @var{Yint}, @var{Ysd}, @var{K}] = regress_gp (@var{X},
## @var{Y}, @var{Xfit}, @qcode{"rbf"})} will estimate a Gaussian Process model
## with a Radial Basis Function (RBF) kernel with default parameters
## @qcode{@var{theta} = 5}, which corresponds to the characteristic lengthscale,
## and @qcode{@var{g} = 0.01}, which corresponds to the nugget effect, and
## @qcode{@var{alpha} = 0.05} which defines the confidence level for the
## estimated intervals returned in @var{Yint}.  The function also returns the
## predictive covariance matrix in @var{Ysd}.  For multidimensional predictors
## @var{X} the function will automatically normalize each column to a zero mean
## and a standard deviation to one.
##
## Run @code{demo regress_gp} to see examples.
##
## @seealso{regress, regression_ftest, regression_ttest}
## @end deftypefn

function [Yfit, Yint, varargout] = regress_gp (X, Y, Xfit, varargin)

  ## Check input arguments
  if (nargin < 3)
    print_usage;
  endif
  if (ndims (X) != 2)
    error ("regress_gp: X must be a 2-D matrix.");
  endif
  if (! isvector (Y) || size (Y, 2) != 1)
    error ("regress_gp: Y must be a column vector.");
  endif
  if (size (X, 1) != length (Y))
    error ("regress_gp: rows in X must equal the length of Y.");
  endif
  if (size (X, 2) != size (Xfit, 2))
    error ("regress_gp: X and XI must have the same number of columns.");
  endif

  ## Add defauts
  kernel = "linear";
  Sp = 100 * eye (size (X, 2) + 1);
  theta = 5;
  g = 0.01;
  alpha = 0.05;

  ## Parse extra arguments
  if (nargin > 3)
    tmp = varargin{1};
    if (ischar (tmp) && strcmpi (tmp, "linear"))
      kernel = "linear";
      sinput = true;
    elseif (ischar (tmp) && strcmpi (tmp, "rbf"))
      kernel = "rbf";
      sinput = true;
    elseif (isnumeric (tmp) && ! isscalar (tmp))
      kernel = "linear";
      sinput = false;
      Sp = tmp;
    elseif (isnumeric (tmp) && isscalar (tmp))
      kernel = "rbf";
      sinput = false;
      theta = tmp;
    else
      error ("regress_gp: invalid 4th argument.");
    endif
  endif
  if (nargin > 4)
    tmp = varargin{2};
    if (sinput)
      if (isnumeric (tmp) && ! isscalar (tmp))
        if (strcmpi (kernel, "rbf"))
          error ("regress_gp: theta must be a scalar when using RBF kernel.");
        endif
        Sp = tmp;
        if (! isequal (size (Sp), (size (X, 2) + 1) * [1, 1]))
          error ("regress_gp: wrong size for prior covariance matrix Sp.");
        endif
      elseif (isnumeric (tmp) && isscalar (tmp))
        if (strcmpi (kernel, "linear"))
          error ("regress_gp: wrong size for prior covariance matrix Sp.");
        endif
        theta = tmp;
      else
        error ("regress_gp: invalid 5th argument.");
      endif
    else
      if (strcmpi (kernel, "linear"))
        error ("regress_gp: invalid 5th argument.");
      endif
      g = tmp;
    endif
  endif
  if (nargin > 5)
    tmp = varargin{3};
    if (isnumeric (tmp) && isscalar (tmp) && sinput)
      g = tmp;
    elseif (isnumeric (tmp) && isscalar (tmp) && ! sinput)
      alpha = tmp;
    else
      error ("regress_gp: invalid 6th argument.");
    endif
  endif
  if (nargin > 6)
    tmp = varargin{4};
    if (isnumeric (tmp) && isscalar (tmp) && sinput)
      alpha = tmp;
    else
      error ("regress_gp: invalid 7th argument.");
    endif
  endif

  ## User linear kernel
  if (strcmpi (kernel, "linear"))
    ## Add constant vector
    X = [ones(1,size(X,1)); X'];

    ## Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
    ## Note that in the book the equation (below 2.11) for the A reads
    ## A  = (1/sy^2)*X*X' + inv (Vp);
    ## where sy is the scalar variance of the of the residuals (i.e Y = X' * w + epsilon)
    ## and epsilon is drawn from N(0,sy^2). Vp is the variance of the parameters w.
    ## Note that
    ## (sy^2 * A)^{-1} = (1/sy^2)*A^{-1} =  (X*X' + sy^2 * inv(Vp))^{-1};
    ## and that the formula for the w mean is
    ## (1/sy^2)*A^{-1}*X*Y
    ## Then one obtains
    ## inv(X*X' + sy^2 * inv(Vp))*X*Y
    ## Looking at the formula bloew we see that Sp = (1/sy^2)*Vp
    ## making the regression depend on only one parameter, Sp, and not two.

    ## Xsq = sum (X' .^ 2);
    ## [n, d] = size (X);
    ## ￼sigma = 1/sqrt(2);
    ## Ks = exp (-(Xsq' * ones (1, n) -ones (n, 1) * Xsq + 2 * X * X') / (2 * sigma ^ 2));

    A  = X * X' + inv (Sp);
    K  = inv (A);
    wm = K * X * Y;

    ## Add constant vector
    Xfit = [ones(size(Xfit,1),1), Xfit];

    ## Compute predictions
    Yfit = Xfit*wm;
    Ysd = Xfit * K * Xfit';
    dy = diag (Ysd);
    Yint = [Yfit+dy, Yfit-dy];
    if (nargout > 2)
      varargout{1} = wm;
    endif
    if (nargout > 3)
      varargout{2} = K;
    endif

  endif

  ## User RBF kernel
  if (strcmpi (kernel, "rbf"))
    ## Normalize predictors
    if (size (X, 2) > 1)
      [X, MU, SIGMA] = zscore (X);
      Xfit = (Xfit - MU) ./SIGMA;
    endif
    ## Get number of training samples
    n = size (X, 1);

    ## Calculate squared distance matrix of training input
    D = squareform (pdist (X) .^2);

    ## Compute kernel covariance for training quantities
    S = exp (-D / theta) + g * eye (n);

    ## Compute kernel covariance for testing quantities
    Dxi = squareform (pdist (Xfit) .^ 2);
    Sxi = exp (-Dxi / theta) + g * eye (size (Dxi, 1));

    ## Compute kernel covariance for prediction
    Dx = pdist2 (Xfit, X) .^ 2;
    Sx = exp (-Dx / theta);

    ## Caculate predictive covariance
    K = inv (S);

    ## Calculate response output
    Yfit = Sx * K * Y;

    ## Estimate scale parameter for predictive variance
    scale = (Y' * K * Y) / size (Y, 1);

    ## Calculate standard deviation of the response output
    Ysd = scale * (Sxi - Sx * K * Sx');
    ysd1 = sqrt (diag (Ysd));

    ## Calculate prediction intervals
    Yint = norminv (alpha, 0, ysd1);
    Yint = [Yfit+Yint, Yfit-Yint];

    if (nargout > 2)
      varargout{1} = Ysd;
    endif
  endif

endfunction

%!demo
%! ## Linear fitting of 1D Data
%! rand ("seed", 125);
%! X = 2 * rand (5, 1) - 1;
%! randn ("seed", 25);
%! Y = 2 * X - 1 + 0.3 * randn (5, 1);
%!
%! ## Points for interpolation/extrapolation
%! Xfit = linspace (-2, 2, 10)';
%!
%! ## Fit regression model
%! [Yfit, Yint, m] = regress_gp (X, Y, Xfit);
%!
%! ## Plot fitted data
%! plot (X, Y, "xk", Xfit, Yfit, "r-", Xfit, Yint, "b-");
%! title ("Gaussian process regression with linear kernel");

%!demo
%! ## Linear fitting of 2D Data
%! rand ("seed", 135);
%! X = 2 * rand (4, 2) - 1;
%! randn ("seed", 35);
%! Y = 2 * X(:,1) - 3 * X(:,2) - 1 + 1 * randn (4, 1);
%!
%! ## Mesh for interpolation/extrapolation
%! [x1, x2] = meshgrid (linspace (-1, 1, 10));
%! Xfit = [x1(:), x2(:)];
%!
%! ## Fit regression model
%! [Ypred, Yint, Ysd] = regress_gp (X, Y, Xfit);
%! Ypred = reshape (Ypred, 10, 10);
%! YintU = reshape (Yint(:,1), 10, 10);
%! YintL = reshape (Yint(:,2), 10, 10);
%!
%! ## Plot fitted data
%! plot3 (X(:,1), X(:,2), Y, ".k", "markersize", 16);
%! hold on;
%! h = mesh (x1, x2, Ypred, zeros (10, 10));
%! set (h, "facecolor", "none", "edgecolor", "yellow");
%! h = mesh (x1, x2, YintU, ones (10, 10));
%! set (h, "facecolor", "none", "edgecolor", "cyan");
%! h = mesh (x1, x2, YintL, ones (10, 10));
%! set (h, "facecolor", "none", "edgecolor", "cyan");
%! hold off
%! axis tight
%! view (75, 25)
%! title ("Gaussian process regression with linear kernel");

%!demo
%! ## Projection over basis function with linear kernel
%! pp = [2, 2, 0.3, 1];
%! n = 10;
%! rand ("seed", 145);
%! X = 2 * rand (n, 1) - 1;
%! randn ("seed", 45);
%! Y = polyval (pp, X) + 0.3 * randn (n, 1);
%!
%! ## Powers
%! px = [sqrt(abs(X)), X, X.^2, X.^3];
%!
%! ## Points for interpolation/extrapolation
%! Xfit = linspace (-1, 1, 100)';
%! pxi = [sqrt(abs(Xfit)), Xfit, Xfit.^2, Xfit.^3];
%!
%! ## Define a prior covariance assuming that the sqrt component is not present
%! Sp = 100 * eye (size (px, 2) + 1);
%! Sp(2,2) = 1; # We don't believe the sqrt(abs(X)) is present
%!
%! ## Fit regression model
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, Sp);
%!
%! ## Plot fitted data
%! plot (X, Y, "xk;Data;", Xfit, Yfit, "r-;Estimation;", ...
%!                         Xfit, polyval (pp, Xfit), "g-;True;");
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), "m-;Upper bound;", Xfit, Yint(:,2), "b-;Lower bound;");
%! hold off
%! title ("Linear kernel over basis function with prior covariance");

%!demo
%! ## Projection over basis function with linear kernel
%! pp = [2, 2, 0.3, 1];
%! n = 10;
%! rand ("seed", 145);
%! X = 2 * rand (n, 1) - 1;
%! randn ("seed", 45);
%! Y = polyval (pp, X) + 0.3 * randn (n, 1);
%!
%! ## Powers
%! px = [sqrt(abs(X)), X, X.^2, X.^3];
%!
%! ## Points for interpolation/extrapolation
%! Xfit = linspace (-1, 1, 100)';
%! pxi = [sqrt(abs(Xfit)), Xfit, Xfit.^2, Xfit.^3];
%!
%! ## Fit regression model without any assumption on prior covariance
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi);
%!
%! ## Plot fitted data
%! plot (X, Y, "xk;Data;", Xfit, Yfit, "r-;Estimation;", ...
%!                         Xfit, polyval (pp, Xfit), "g-;True;");
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), "m-;Upper bound;", Xfit, Yint(:,2), "b-;Lower bound;");
%! hold off
%! title ("Linear kernel over basis function without prior covariance");

%!demo
%! ## Projection over basis function with rbf kernel
%! pp = [2, 2, 0.3, 1];
%! n = 10;
%! rand ("seed", 145);
%! X = 2 * rand (n, 1) - 1;
%! randn ("seed", 45);
%! Y = polyval (pp, X) + 0.3 * randn (n, 1);
%!
%! ## Powers
%! px = [sqrt(abs(X)), X, X.^2, X.^3];
%!
%! ## Points for interpolation/extrapolation
%! Xfit = linspace (-1, 1, 100)';
%! pxi = [sqrt(abs(Xfit)), Xfit, Xfit.^2, Xfit.^3];
%!
%! ## Fit regression model with RBF kernel (standard parameters)
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, "rbf");
%!
%! ## Plot fitted data
%! plot (X, Y, "xk;Data;", Xfit, Yfit, "r-;Estimation;", ...
%!                         Xfit, polyval (pp, Xfit), "g-;True;");
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), "m-;Upper bound;", Xfit, Yint(:,2), "b-;Lower bound;");
%! hold off
%! title ("RBF kernel over basis function with standard parameters");
%! text (-0.5, 4, "theta = 5\n g = 0.01");

%!demo
%! ## Projection over basis function with rbf kernel
%! pp = [2, 2, 0.3, 1];
%! n = 10;
%! rand ("seed", 145);
%! X = 2 * rand (n, 1) - 1;
%! randn ("seed", 45);
%! Y = polyval (pp, X) + 0.3 * randn (n, 1);
%!
%! ## Powers
%! px = [sqrt(abs(X)), X, X.^2, X.^3];
%!
%! ## Points for interpolation/extrapolation
%! Xfit = linspace (-1, 1, 100)';
%! pxi = [sqrt(abs(Xfit)), Xfit, Xfit.^2, Xfit.^3];
%!
%! ## Fit regression model with RBF kernel with different parameters
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, "rbf", 10, 0.01);
%!
%! ## Plot fitted data
%! plot (X, Y, "xk;Data;", Xfit, Yfit, "r-;Estimation;", ...
%!                         Xfit, polyval (pp, Xfit), "g-;True;");
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), "m-;Upper bound;", Xfit, Yint(:,2), "b-;Lower bound;");
%! hold off
%! title ("GP regression with RBF kernel and non default parameters");
%! text (-0.5, 4, "theta = 10\n g = 0.01");
%!
%! ## Fit regression model with RBF kernel with different parameters
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, "rbf", 50, 0.01);
%!
%! ## Plot fitted data
%! figure
%! plot (X, Y, "xk;Data;", Xfit, Yfit, "r-;Estimation;", ...
%!                         Xfit, polyval (pp, Xfit), "g-;True;");
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), "m-;Upper bound;", Xfit, Yint(:,2), "b-;Lower bound;");
%! hold off
%! title ("GP regression with RBF kernel and non default parameters");
%! text (-0.5, 4, "theta = 50\n g = 0.01");
%!
%! ## Fit regression model with RBF kernel with different parameters
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, "rbf", 50, 0.001);
%!
%! ## Plot fitted data
%! figure
%! plot (X, Y, "xk;Data;", Xfit, Yfit, "r-;Estimation;", ...
%!                         Xfit, polyval (pp, Xfit), "g-;True;");
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), "m-;Upper bound;", Xfit, Yint(:,2), "b-;Lower bound;");
%! hold off
%! title ("GP regression with RBF kernel and non default parameters");
%! text (-0.5, 4, "theta = 50\n g = 0.001");
%!
%! ## Fit regression model with RBF kernel with different parameters
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, "rbf", 50, 0.05);
%!
%! ## Plot fitted data
%! figure
%! plot (X, Y, "xk;Data;", Xfit, Yfit, "r-;Estimation;", ...
%!                         Xfit, polyval (pp, Xfit), "g-;True;");
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), "m-;Upper bound;", Xfit, Yint(:,2), "b-;Lower bound;");
%! hold off
%! title ("GP regression with RBF kernel and non default parameters");
%! text (-0.5, 4, "theta = 50\n g = 0.05");

%!demo
%! ## RBF fitting on noiseless 1D Data
%! x = [0:2*pi/7:2*pi]';
%! y = 5 * sin (x);
%!
%! ## Predictive grid of 500 equally spaced locations
%! xi = [-0.5:(2*pi+1)/499:2*pi+0.5]';
%!
%! ## Fit regression model with RBF kernel
%! [Yfit, Yint, Ysd] = regress_gp (x, y, xi, "rbf");
%!
%! ## Plot fitted data
%! r = mvnrnd (Yfit, diag (Ysd)', 50);
%! plot (xi, r', "c-");
%! hold on
%! plot (xi, Yfit, "r-", xi, Yint, "b-");
%! plot (x, y, ".k", "markersize", 20)
%! plot (xi, 5 * sin (xi), "-y");
%! xlim ([-0.5,2*pi+0.5]);
%! ylim ([-10,10]);
%! hold off
%! title ("GP regression with RBF kernel on noiseless 1D data");
%! text (-0.5, 4, "theta = 5\n g = 0.01");

%!demo
%! ## RBF fitting on noisy 1D Data
%! x = [0:2*pi/7:2*pi]';
%! x = [x; x];
%! y = 5 * sin (x) + randn (size (x));
%!
%! ## Predictive grid of 500 equally spaced locations
%! xi = [-0.5:(2*pi+1)/499:2*pi+0.5]';
%!
%! ## Fit regression model with RBF kernel
%! [Yfit, Yint, Ysd] = regress_gp (x, y, xi, "rbf");
%!
%! ## Plot fitted data
%! r = mvnrnd (Yfit, diag (Ysd)', 50);
%! plot (xi, r', "c-");
%! hold on
%! plot (xi, Yfit, "r-", xi, Yint, "b-");
%! plot (x, y, ".k", "markersize", 20)
%! plot (xi, 5 * sin (xi), "-y");
%! xlim ([-0.5,2*pi+0.5]);
%! ylim ([-10,10]);
%! hold off
%! title ("GP regression with RBF kernel on noisy 1D data");
%! text (-0.5, 4, "theta = 5\n g = 0.01");

## Test input validation
%!error<Invalid call to regress_gp.> regress_gp (ones (20, 2))
%!error<Invalid call to regress_gp.> regress_gp (ones (20, 2), ones (20, 1))
%!error<regress_gp: X must be a 2-D matrix.> ...
%! regress_gp (ones (20, 2, 3), ones (20, 1), ones (20, 2))
%!error<regress_gp: Y must be a column vector.> ...
%! regress_gp (ones (20, 2), ones (20, 2), ones (20, 2))
%!error<regress_gp: rows in X must equal the length of Y.> ...
%! regress_gp (ones (20, 2), ones (15, 1), ones (20, 2))
%!error<regress_gp: X and XI must have the same number of columns.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (20, 3))
%!error<regress_gp: invalid 4th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), {[3]})
%!error<regress_gp: invalid 4th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "kernel")
%!error<regress_gp: theta must be a scalar when using RBF kernel.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "rbf", ones (4))
%!error<regress_gp: wrong size for prior covariance matrix Sp.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "linear", 1)
%!error<regress_gp: invalid 5th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "rbf", "value")
%!error<regress_gp: invalid 5th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "rbf", {5})
%!error<regress_gp: invalid 5th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), ones (3), 5)
%!error<regress_gp: wrong size for prior covariance matrix Sp.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "linear", 5)
%!error<regress_gp: invalid 6th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "rbf", 5, {5})
%!error<regress_gp: invalid 6th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "rbf", 5, ones (2))
%!error<regress_gp: invalid 6th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 5, 0.01, [1, 1])
%!error<regress_gp: invalid 6th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 5, 0.01, "f")
%!error<regress_gp: invalid 6th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 5, 0.01, "f")
%!error<regress_gp: invalid 7th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "rbf", 5, 0.01, "f")
%!error<regress_gp: invalid 7th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "rbf", 5, 0.01, [1, 1])
%!error<regress_gp: wrong size for prior covariance matrix Sp.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), "linear", 1)
