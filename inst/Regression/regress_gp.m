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
## @deftypefnx {statistics} {[@var{Yfit}, @var{Yint}, @var{m}, @var{K}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{'linear'})
## @deftypefnx {statistics} {[@var{Yfit}, @var{Yint}, @var{Ysd}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{'rbf'})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{'linear'}, @var{Sp})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @var{Sp})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{'rbf'}, @var{theta})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{'rbf'}, @var{theta}, @var{g})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @qcode{'rbf'}, @var{theta}, @var{g}, @var{alpha})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @var{theta})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @var{theta}, @var{g})
## @deftypefnx {statistics} {[@dots{}] =} regress_gp (@var{X}, @var{Y}, @var{Xfit}, @var{theta}, @var{g}, @var{alpha})
##
## Regression using Gaussian Processes.
##
## @code{[@var{Yfit}, @var{Yint}, @var{m}, @var{K}] = regress_gp (@var{X},
## @var{Y}, @var{Xfit})} will estimate a linear Gaussian Process model @var{m}
## in the form @qcode{@var{Y} = @var{X}' * @var{m}}, where @var{X} is an
## @math{N*P} matrix with @math{N} observations in @math{P} dimensional space
## and @var{Y} is an @math{N*1} column vector as the dependent variable.  The
## information about errors of the predictions (interpolation/extrapolation) is
## given by the covariance matrix @var{K}.
## By default, the linear model defines the prior covariance of @var{m} as
## @code{@var{Sp} = 100 * eye (size (@var{X}, 2) + 1)}.  A custom prior
## covariance matrix can be passed as @var{Sp}, which must be a @math{P+1*P+1}
## positive definite matrix.  The model is evaluated for input @var{Xfit}, which
## must have the same columns as @var{X}, and the estimates are returned in
## @var{Yfit} along with the estimated variation in @var{Yint}.
## @qcode{@var{Yint}(:,1)} contains the lower boundary and
## @qcode{@var{Yint}(:,2)} the upper boundary of the interval about
## @var{Yfit}, at the confidence level @qcode{1 - @var{alpha}}.
##
## @code{[@var{Yfit}, @var{Yint}, @var{Ysd}] = regress_gp (@var{X},
## @var{Y}, @var{Xfit}, @qcode{'rbf'})} will estimate a Gaussian Process model
## with a Radial Basis Function (RBF) kernel with default parameters
## @qcode{@var{theta} = 5} and @qcode{@var{g} = 0.01}, which corresponds to the
## nugget effect, and
## @qcode{@var{alpha} = 0.05} which defines the confidence level for the
## estimated intervals returned in @var{Yint}.  The function also returns the
## predictive covariance matrix in @var{Ysd}.  For multidimensional predictors
## @var{X} the function will automatically normalize each column to a zero mean
## and a standard deviation to one.
##
## Four things about the RBF kernel are worth stating, because they decide
## what the numbers mean.
##
## @itemize
## @item
## @var{theta} is not the characteristic lengthscale.  The kernel is
## @code{exp (-d^2 / @var{theta})} with @math{d} the distance between two
## points, so @var{theta} is twice the square of a lengthscale @math{l}, and
## @code{@var{theta} = 2 * l^2}.
##
## @item
## The intervals in @var{Yint} are prediction intervals for a @emph{new
## observation}, not confidence intervals on the mean: the nugget is carried
## in the predictive variance, so the noise of an observation is included.
##
## @item
## The predictors are centred and scaled only when there is more than one of
## them, so @var{theta} is measured in the units of @var{X} for a single
## predictor and in standard deviations for several.
##
## @item
## A nugget of zero is accepted and gives an interpolating process, one that
## reproduces @var{Y} at the training points and reports almost no uncertainty
## there.  It also leaves the kernel matrix rank deficient whenever two inputs
## are close, so in that case the covariance is applied through its
## pseudoinverse and the fit is the minimum-norm solution.  This is a genuine
## answer rather than a refusal, but a small nugget is the better way to ask
## for a smooth fit.
## @end itemize
##
## Run @code{demo regress_gp} to see examples.
##
## @seealso{fitrgp, RegressionGP, regress, regression_ftest,
## regression_ttest}
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

  ## Add defaults
  kernel = 'linear';
  Sp = 100 * eye (size (X, 2) + 1);
  theta = 5;
  g = 0.01;
  alpha = 0.05;

  ## Parse extra arguments
  if (nargin > 3)
    tmp = varargin{1};
    if (ischar (tmp) && strcmpi (tmp, 'linear'))
      kernel = 'linear';
      sinput = true;
    elseif (ischar (tmp) && strcmpi (tmp, 'rbf'))
      kernel = 'rbf';
      sinput = true;
    elseif (isnumeric (tmp) && ! isscalar (tmp))
      kernel = 'linear';
      sinput = false;
      Sp = checkSp (tmp, size (X, 2));
    elseif (isnumeric (tmp) && isscalar (tmp))
      kernel = 'rbf';
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
        if (strcmpi (kernel, 'rbf'))
          error ("regress_gp: theta must be a scalar when using RBF kernel.");
        endif
        Sp = checkSp (tmp, size (X, 2));
      elseif (isnumeric (tmp) && isscalar (tmp))
        if (strcmpi (kernel, 'linear'))
          error ("regress_gp: wrong size for prior covariance matrix Sp.");
        endif
        theta = tmp;
      else
        error ("regress_gp: invalid 5th argument.");
      endif
    else
      if (strcmpi (kernel, 'linear'))
        error ("regress_gp: invalid 5th argument.");
      endif
      ## Every other argument in this function is validated before it is
      ## used.  This one was not, so a non-numeric nugget reached the kernel
      ## and failed inside a matrix product with a message about
      ## nonconformant arguments, naming neither g nor this function.
      if (! (isnumeric (tmp) && isscalar (tmp)))
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
  if (strcmpi (kernel, 'linear'))
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
    ## Looking at the formula below we see that Sp = (1/sy^2)*Vp
    ## making the regression depend on only one parameter, Sp, and not two.

    ## Xsq = sum (X' .^ 2);
    ## [n, d] = size (X);
    ## sigma = 1/sqrt(2);
    ## Ks = exp (-(Xsq' * ones (1, n) -ones (n, 1) * Xsq + 2 * X * X') / (2 * sigma ^ 2));

    ## Sp and A are both positive definite, so their inverses are taken
    ## through a Cholesky factor rather than by inv.  A is built from a Gram
    ## matrix and inherits the square of the predictors' conditioning, which is
    ## exactly where a general inverse loses digits and a triangular solve does
    ## not.  K is still formed because it is the fourth output.
    Rp = chol (Sp);
    A  = X * X' + Rp \ (Rp' \ eye (rows (Sp)));
    Ra = chol (A);
    K  = Ra \ (Ra' \ eye (rows (A)));
    wm = Ra \ (Ra' \ (X * Y));

    ## Add constant vector
    Xfit = [ones(size(Xfit,1),1), Xfit];

    ## Compute predictions
    Yfit = Xfit*wm;
    ## Only the diagonal of Xfit * K * Xfit' is ever read here, so the full
    ## m-by-m matrix is never formed: with A = Ra' * Ra the product is U' * U
    ## for U = Ra' \ Xfit', whose diagonal is the squared column norms of U.
    ## Exact, not an approximation, and O(m) in place of O(m^2).
    U = Ra' \ Xfit';
    Ysd = sum (U .^ 2, 1)';
    ## The diagonal of a covariance is a variance, so it is the square root
    ## that scales the interval, and the interval is that many standard
    ## deviations wide for the confidence level asked for.  This branch used
    ## the variance itself, applied no level at all, and returned its columns
    ## in the opposite order to the other branch.
    dy = norminv (1 - alpha / 2) * sqrt (Ysd);
    Yint = [Yfit-dy, Yfit+dy];
    if (nargout > 2)
      varargout{1} = wm;
    endif
    if (nargout > 3)
      varargout{2} = K;
    endif

  endif

  ## User RBF kernel
  if (strcmpi (kernel, 'rbf'))
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

    ## Compute kernel covariance for prediction
    Dx = pdist2 (Xfit, X) .^ 2;
    Sx = exp (-Dx / theta);

    ## S is symmetric, and positive definite whenever the nugget is, so every
    ## quantity below goes through a solve rather than through inv (S).  The
    ## inverse was not merely slower: at g = 0 it warns "matrix singular to
    ## machine precision" with an rcond of 1e-18 and returns what that implies,
    ## on an argument the function accepts.
    SinvY = Y;
    SinvSx = Sx';
    if (rcond (S) > eps)
      SinvY = S \ Y;
      SinvSx = S \ Sx';
    else
      ## A nugget of zero on repeated or near-repeated inputs leaves S rank
      ## deficient.  The minimum-norm solution is the defensible answer here
      ## and pinv gives it without the warning inv emitted.
      Sinv = pinv (S);
      SinvY = Sinv * Y;
      SinvSx = Sinv * Sx';
    endif

    ## Calculate response output
    Yfit = Sx * SinvY;

    ## Estimate scale parameter for predictive variance
    scale = (Y' * SinvY) / size (Y, 1);

    ## The interval needs only the diagonal of the predictive covariance, and
    ## the test-test prior contributes exactly 1 + g to it, since a point is at
    ## zero distance from itself.  So neither the m-by-m distance matrix nor
    ## the full covariance is built unless the caller asks for the covariance
    ## itself as the third output.
    ysd1 = sqrt (scale * ((1 + g) - sum (Sx' .* SinvSx, 1)'));
    if (nargout > 2)
      Dxi = squareform (pdist (Xfit) .^ 2);
      Sxi = exp (-Dxi / theta) + g * eye (rows (Xfit));
      Ysd = scale * (Sxi - Sx * SinvSx);
    endif

    ## Calculate prediction intervals.  A two sided interval puts alpha/2 in
    ## each tail, so the quantile is of 1 - alpha/2; taking it at alpha
    ## returned a 100*(1-2*alpha) per cent interval under a 100*(1-alpha) per
    ## cent name, which at the default made every interval a 90 per cent one.
    dy = norminv (1 - alpha / 2) * ysd1;
    Yint = [Yfit-dy, Yfit+dy];

    if (nargout > 2)
      varargout{1} = Ysd;
    endif
  endif

endfunction

## Validate the prior covariance of the weights.  Both argument positions that
## accept Sp check it here, because they did not check the same things: the
## fourth-argument form validated nothing at all, and the fifth validated the
## size but not the definiteness the documentation requires.  A singular Sp is
## inverted to infinities and every prediction comes back NaN, warning about a
## matrix the caller never handed to inv.
function Sp = checkSp (Sp, p)

  if (! isequal (size (Sp), (p + 1) * [1, 1]))
    error ("regress_gp: wrong size for prior covariance matrix Sp.");
  endif
  [~, spd] = chol (Sp);
  if (spd != 0)
    error (strcat ("regress_gp: prior covariance matrix Sp must be", ...
                   " positive definite."));
  endif

endfunction

%!demo
%! ## Linear fitting of 1D Data
%! rng (42);
%! X = 2 * rand (5, 1) - 1;
%! Y = 2 * X - 1 + 0.3 * randn (5, 1);
%!
%! ## Points for interpolation/extrapolation
%! Xfit = linspace (-2, 2, 10)';
%!
%! ## Fit regression model
%! [Yfit, Yint, m] = regress_gp (X, Y, Xfit);
%!
%! ## Plot fitted data
%! plot (X, Y, 'xk', Xfit, Yfit, 'r-', Xfit, Yint, 'b-');
%! title ('Gaussian process regression with linear kernel');

%!demo
%! ## Linear fitting of 2D Data
%! rng (42);
%! X = 2 * rand (4, 2) - 1;
%! Y = 2 * X(:,1) - 3 * X(:,2) - 1 + 1 * randn (4, 1);
%!
%! ## Mesh for interpolation/extrapolation
%! [x1, x2] = meshgrid (linspace (-1, 1, 10));
%! Xfit = [x1(:), x2(:)];
%!
%! ## Fit regression model
%! [Ypred, Yint, m] = regress_gp (X, Y, Xfit);
%! Ypred = reshape (Ypred, 10, 10);
%! YintL = reshape (Yint(:,1), 10, 10);
%! YintU = reshape (Yint(:,2), 10, 10);
%!
%! ## Plot fitted data
%! plot3 (X(:,1), X(:,2), Y, '.k', 'markersize', 16);
%! hold on;
%! h = mesh (x1, x2, Ypred, zeros (10, 10));
%! set (h, 'facecolor', 'none', 'edgecolor', 'yellow');
%! h = mesh (x1, x2, YintU, ones (10, 10));
%! set (h, 'facecolor', 'none', 'edgecolor', 'cyan');
%! h = mesh (x1, x2, YintL, ones (10, 10));
%! set (h, 'facecolor', 'none', 'edgecolor', 'cyan');
%! hold off
%! axis tight
%! view (75, 25)
%! title ('Gaussian process regression with linear kernel');

%!demo
%! ## Projection over basis function with linear kernel
%! rng (42);
%! pp = [2, 2, 0.3, 1];
%! n = 10;
%! X = 2 * rand (n, 1) - 1;
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
%! [Yfit, Yint, m] = regress_gp (px, Y, pxi, Sp);
%!
%! ## Plot fitted data
%! plot (X, Y, 'xk;Data;', Xfit, Yfit, 'r-;Estimation;', ...
%!                         Xfit, polyval (pp, Xfit), 'g-;True;');
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), 'b-;Lower bound;', ...
%!       Xfit, Yint(:,2), 'm-;Upper bound;');
%! hold off
%! title ('Linear kernel over basis function with prior covariance');

%!demo
%! ## Projection over basis function with linear kernel
%! rng (42);
%! pp = [2, 2, 0.3, 1];
%! n = 10;
%! X = 2 * rand (n, 1) - 1;
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
%! [Yfit, Yint, m] = regress_gp (px, Y, pxi);
%!
%! ## Plot fitted data
%! plot (X, Y, 'xk;Data;', Xfit, Yfit, 'r-;Estimation;', ...
%!                         Xfit, polyval (pp, Xfit), 'g-;True;');
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), 'b-;Lower bound;', ...
%!       Xfit, Yint(:,2), 'm-;Upper bound;');
%! hold off
%! title ('Linear kernel over basis function without prior covariance');

%!demo
%! ## Projection over basis function with rbf kernel
%! rng (42);
%! pp = [2, 2, 0.3, 1];
%! n = 10;
%! X = 2 * rand (n, 1) - 1;
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
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, 'rbf');
%!
%! ## Plot fitted data
%! plot (X, Y, 'xk;Data;', Xfit, Yfit, 'r-;Estimation;', ...
%!                         Xfit, polyval (pp, Xfit), 'g-;True;');
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), 'b-;Lower bound;', ...
%!       Xfit, Yint(:,2), 'm-;Upper bound;');
%! hold off
%! title ('RBF kernel over basis function with standard parameters');
%! text (-0.5, 4, "theta = 5\n g = 0.01");

%!demo
%! ## Projection over basis function with rbf kernel
%! rng (42);
%! pp = [2, 2, 0.3, 1];
%! n = 10;
%! X = 2 * rand (n, 1) - 1;
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
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, 'rbf', 10, 0.01);
%!
%! ## Plot fitted data
%! plot (X, Y, 'xk;Data;', Xfit, Yfit, 'r-;Estimation;', ...
%!                         Xfit, polyval (pp, Xfit), 'g-;True;');
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), 'b-;Lower bound;', ...
%!       Xfit, Yint(:,2), 'm-;Upper bound;');
%! hold off
%! title ('GP regression with RBF kernel and non default parameters');
%! text (-0.5, 4, "theta = 10\n g = 0.01");
%!
%! ## Fit regression model with RBF kernel with different parameters
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, 'rbf', 50, 0.01);
%!
%! ## Plot fitted data
%! figure
%! plot (X, Y, 'xk;Data;', Xfit, Yfit, 'r-;Estimation;', ...
%!                         Xfit, polyval (pp, Xfit), 'g-;True;');
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), 'b-;Lower bound;', ...
%!       Xfit, Yint(:,2), 'm-;Upper bound;');
%! hold off
%! title ('GP regression with RBF kernel and non default parameters');
%! text (-0.5, 4, "theta = 50\n g = 0.01");
%!
%! ## Fit regression model with RBF kernel with different parameters
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, 'rbf', 50, 0.001);
%!
%! ## Plot fitted data
%! figure
%! plot (X, Y, 'xk;Data;', Xfit, Yfit, 'r-;Estimation;', ...
%!                         Xfit, polyval (pp, Xfit), 'g-;True;');
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), 'b-;Lower bound;', ...
%!       Xfit, Yint(:,2), 'm-;Upper bound;');
%! hold off
%! title ('GP regression with RBF kernel and non default parameters');
%! text (-0.5, 4, "theta = 50\n g = 0.001");
%!
%! ## Fit regression model with RBF kernel with different parameters
%! [Yfit, Yint, Ysd] = regress_gp (px, Y, pxi, 'rbf', 50, 0.05);
%!
%! ## Plot fitted data
%! figure
%! plot (X, Y, 'xk;Data;', Xfit, Yfit, 'r-;Estimation;', ...
%!                         Xfit, polyval (pp, Xfit), 'g-;True;');
%! axis tight
%! axis manual
%! hold on
%! plot (Xfit, Yint(:,1), 'b-;Lower bound;', ...
%!       Xfit, Yint(:,2), 'm-;Upper bound;');
%! hold off
%! title ('GP regression with RBF kernel and non default parameters');
%! text (-0.5, 4, "theta = 50\n g = 0.05");

%!demo
%! ## RBF fitting on noiseless 1D Data
%! rng (42);
%! x = [0:2*pi/7:2*pi]';
%! y = 5 * sin (x);
%!
%! ## Predictive grid of 500 equally spaced locations
%! xi = [-0.5:(2*pi+1)/499:2*pi+0.5]';
%!
%! ## Fit regression model with RBF kernel
%! [Yfit, Yint, Ysd] = regress_gp (x, y, xi, 'rbf');
%!
%! ## Plot fitted data
%! r = mvnrnd (Yfit, diag (Ysd)', 50);
%! plot (xi, r', 'c-');
%! hold on
%! plot (xi, Yfit, 'r-;Estimation;', xi, Yint, 'b-;Confidence interval;');
%! plot (x, y, '.k;Predictor points;', 'markersize', 20)
%! plot (xi, 5 * sin (xi), '-y;True Function;');
%! xlim ([-0.5,2*pi+0.5]);
%! ylim ([-10,10]);
%! hold off
%! title ('GP regression with RBF kernel on noiseless 1D data');
%! text (0, -7, "theta = 5\n g = 0.01");

%!demo
%! ## RBF fitting on noisy 1D Data
%! rng (42);
%! x = [0:2*pi/7:2*pi]';
%! x = [x; x];
%! y = 5 * sin (x) + randn (size (x));
%!
%! ## Predictive grid of 500 equally spaced locations
%! xi = [-0.5:(2*pi+1)/499:2*pi+0.5]';
%!
%! ## Fit regression model with RBF kernel
%! [Yfit, Yint, Ysd] = regress_gp (x, y, xi, 'rbf');
%!
%! ## Plot fitted data
%! r = mvnrnd (Yfit, diag (Ysd)', 50);
%! plot (xi, r', 'c-');
%! hold on
%! plot (xi, Yfit, 'r-;Estimation;', xi, Yint, 'b-;Confidence interval;');
%! plot (x, y, '.k;Predictor points;', 'markersize', 20)
%! plot (xi, 5 * sin (xi), '-y;True Function;');
%! xlim ([-0.5,2*pi+0.5]);
%! ylim ([-10,10]);
%! hold off
%! title ('GP regression with RBF kernel on noisy 1D data');
%! text (0, -7, "theta = 5\n g = 0.01");

## Test input validation
## The value tests this function never had.  Its two branches are checked
## against closed forms computed here, and the RBF branch additionally against
## RegressionGP, whose arithmetic is verified against MATLAB.

%!test
%! ## The RBF mean is the Gaussian process posterior mean, which RegressionGP
%! ## computes too.  Matching the parameterisations: this function writes the
%! ## kernel as exp (-d2 / theta) with a noise variance g added, where
%! ## RegressionGP writes it as SigmaF^2 * exp (-0.5 * d2 / SigmaL^2) with a
%! ## noise standard deviation, so SigmaL = sqrt (theta/2), SigmaF = 1 and
%! ## Sigma = sqrt (g).
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x);
%! xq = [0.15; 0.45; 0.75];
%! theta = 0.5;
%! g = 0.01;
%! yf = regress_gp (x, y, xq, 'rbf', theta, g);
%! Mdl = RegressionGP (x, y, 'KernelFunction', 'squaredexponential', ...
%!                     'KernelParameters', [sqrt(theta/2); 1], ...
%!                     'Sigma', sqrt (g), 'BasisFunction', 'none', ...
%!                     'FitMethod', 'none');
%! assert_equal (yf, predict (Mdl, xq), 1e-10);

%!test
%! ## The RBF interval is the normal quantile of the level times the standard
%! ## deviation, on both sides of the fit
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x);
%! xq = [0.2; 0.5; 0.8];
%! [yf, yi, ys] = regress_gp (x, y, xq, 'rbf');
%! sd = sqrt (diag (ys));
%! assert_equal (yi(:,1), yf - norminv (0.975) * sd, 1e-12);
%! assert_equal (yi(:,2), yf + norminv (0.975) * sd, 1e-12);

%!test
%! ## A looser level gives a narrower interval, and the level reaches the
%! ## interval at all, which it did not before
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x);
%! xq = [0.2; 0.5; 0.8];
%! [~, yi95] = regress_gp (x, y, xq, 'rbf', 5, 0.01, 0.05);
%! [~, yi80] = regress_gp (x, y, xq, 'rbf', 5, 0.01, 0.20);
%! assert (all (yi80(:,2) - yi80(:,1) < yi95(:,2) - yi95(:,1)));

%!test
%! ## The linear mean is the Bayesian linear regression posterior mean
%! x = linspace (0, 1, 15)';
%! y = 2 * x + 0.5;
%! xq = [0.25; 0.75];
%! Sp = 100 * eye (2);
%! yf = regress_gp (x, y, xq);
%! Xd = [ones(1, 15); x'];
%! wm = inv (Xd * Xd' + inv (Sp)) * Xd * y;
%! assert_equal (yf, [ones(2, 1), xq] * wm, 1e-12);

%!test
%! ## The linear interval is built from a standard deviation and carries the
%! ## confidence level, where it used to be a bare variance
%! x = linspace (0, 1, 15)';
%! y = 2 * x + 0.5;
%! xq = [0.25; 0.75];
%! [yf, yi, wm, K] = regress_gp (x, y, xq);
%! sd = sqrt (diag ([ones(2, 1), xq] * K * [ones(2, 1), xq]'));
%! assert_equal (yi(:,1), yf - norminv (0.975) * sd, 1e-12);
%! assert_equal (yi(:,2), yf + norminv (0.975) * sd, 1e-12);

%!test
%! ## Both branches order the interval the same way, lower bound first, as
%! ## every other prediction interval in the package does
%! x = linspace (0, 1, 15)';
%! y = sin (3*x);
%! xq = [0.3; 0.6];
%! [~, yiL] = regress_gp (x, y, xq, 'linear');
%! [~, yiR] = regress_gp (x, y, xq, 'rbf');
%! assert (all (yiL(:,1) < yiL(:,2)));
%! assert (all (yiR(:,1) < yiR(:,2)));

## The Cholesky path and the predictive covariance must agree: the interval is
## built from the diagonal computed directly, the third output from the full
## matrix, and the two are formed by different routes.
%!test
%! rand ('seed', 91);
%! x = 2 * rand (12, 1) - 1;
%! y = sin (3 * x);
%! xq = linspace (-1, 1, 6)';
%! [yf, yi, S] = regress_gp (x, y, xq, 'rbf');
%! assert_equal (yi(:,2) - yf, norminv (0.975) * sqrt (diag (S)), 1e-12);
%! assert_equal (yf - yi(:,1), norminv (0.975) * sqrt (diag (S)), 1e-12);

## Asking for the covariance must not change the fit or the interval.
%!test
%! rand ('seed', 91);
%! x = 2 * rand (12, 1) - 1;
%! y = sin (3 * x);
%! xq = linspace (-1, 1, 6)';
%! [yf2, yi2] = regress_gp (x, y, xq, 'rbf');
%! [yf3, yi3, S] = regress_gp (x, y, xq, 'rbf');
%! assert_equal (yf2, yf3);
%! assert_equal (yi2, yi3);

## Without a nugget the process interpolates its training data.  This used to
## be computed through inv () on a kernel matrix with an rcond of 1e-18, which
## warned and returned what that implies.
%!test
%! rand ('seed', 3);
%! x = sort (2 * rand (10, 1) - 1);
%! y = sin (3 * x) + 0.1;
%! [yf, yi] = regress_gp (x, y, x, 'rbf', 1, 0);
%! assert_equal (yf, y, 1e-4);
%! assert (max (abs (yi(:,2) - yf)) < 1e-2);

## The linear branch reports the posterior covariance of the weights as its
## fourth output, and the interval must follow that same matrix.
%!test
%! rand ('seed', 17);
%! x = 2 * rand (10, 2) - 1;
%! y = x(:,1) - 2 * x(:,2) + 0.1;
%! xq = 2 * rand (4, 2) - 1;
%! [yf, yi, wm, K] = regress_gp (x, y, xq);
%! xa = [ones(4, 1), xq];
%! assert_equal (yi(:,2) - yf, ...
%!               norminv (0.975) * sqrt (diag (xa * K * xa')), 1e-12);

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
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'kernel')
%!error<regress_gp: theta must be a scalar when using RBF kernel.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'rbf', ones (4))
%!error<regress_gp: invalid 5th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (20, 2), 5, 'junk')
%!error<regress_gp: prior covariance matrix Sp must be positive definite.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (20, 2), zeros (3, 3))
%!error<regress_gp: prior covariance matrix Sp must be positive definite.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (20, 2), 'linear', zeros (3, 3))
%!error<regress_gp: wrong size for prior covariance matrix Sp.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (20, 2), eye (4))
%!error<regress_gp: wrong size for prior covariance matrix Sp.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'linear', 1)
%!error<regress_gp: invalid 5th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'rbf', 'value')
%!error<regress_gp: invalid 5th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'rbf', {5})
%!error<regress_gp: invalid 5th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), eye (3), 5)
%!error<regress_gp: wrong size for prior covariance matrix Sp.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'linear', 5)
%!error<regress_gp: invalid 6th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'rbf', 5, {5})
%!error<regress_gp: invalid 6th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'rbf', 5, ones (2))
%!error<regress_gp: invalid 6th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 5, 0.01, [1, 1])
%!error<regress_gp: invalid 6th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 5, 0.01, 'f')
%!error<regress_gp: invalid 6th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 5, 0.01, 'f')
%!error<regress_gp: invalid 7th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'rbf', 5, 0.01, 'f')
%!error<regress_gp: invalid 7th argument.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'rbf', 5, 0.01, [1, 1])
%!error<regress_gp: wrong size for prior covariance matrix Sp.> ...
%! regress_gp (ones (20, 2), ones (20, 1), ones (10, 2), 'linear', 1)
