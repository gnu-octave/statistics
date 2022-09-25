## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {[@var{intercept}, @var{slope}, @var{dev}, @var{dl}, @var{d2l}, @var{P}, @var{stats}] =} logistic_regression (@var{y}, @var{x}, @var{print}, @var{intercept}, @var{slope})
##
## Perform ordinal logistic regression.
##
## Suppose @var{y} takes values in k ordered categories, and let
## @code{P_i (@var{x})} be the cumulative probability that @var{y}
## falls in one of the first i categories given the covariate
## @var{x}.  Then
##
## @example
## [@var{intercept}, @var{slope}] = logistic_regression (@var{y}, @var{x})
## @end example
##
## @noindent
## fits the model
##
## @example
## logit (P_i (@var{x})) = @var{x} * @var{slope} + @var{intercept}_i,   i = 1 @dots{} k-1
## @end example
##
## The number of ordinal categories, k, is taken to be the number
## of distinct values of @code{round (@var{y})}.  If k equals 2,
## @var{y} is binary and the model is ordinary logistic regression.  The
## matrix @var{x} is assumed to have full column rank.
##
## Given @var{y} only, @code{@var{intercept} = logistic_regression (@var{y})}
## fits the model with baseline logit odds only.
##
## The full form is
##
## @example
## @group
## [@var{intercept}, @var{slope}, @var{dev}, @var{dl}, @var{d2l}, @var{P}, @var{stats}]
##    = logistic_regression (@var{y}, @var{x}, @var{print}, @var{intercept}, @var{slope})
## @end group
## @end example
##
## @noindent
## in which all output arguments and all input arguments except @var{y}
## are optional.
##
## Setting @var{print} to 1 requests summary information about the fitted
## model to be displayed.  Setting @var{print} to 2 requests information
## about convergence at each iteration.  Other values request no
## information to be displayed.  The input arguments @var{intercept} and
## @var{slope} give initial estimates for @var{intercept} and @var{slope}.
##
## The returned value @var{dev} holds minus twice the log-likelihood.
##
## The returned values @var{dl} and @var{d2l} are the vector of first
## and the matrix of second derivatives of the log-likelihood with
## respect to @var{intercept} and @var{slope}.
##
## @var{P} holds estimates for the conditional distribution of @var{y}
## given @var{x}.
##
## @var{stats} returns a structure that contains the following fields:
## @itemize
## @item
##   "intercept": intercept coefficients
## @item
##    "slope": slope coefficients
## @item
##    "dfe": degrees of freedom for error
## @item
##    "coeff": regression coefficients (intercepts and slops)
## @item
##    "covb": estimated covariance matrix for coefficients (coeff)
## @item
##    "coeffcorr": correlation matrix for coeff
## @item
##    "se": standard errors of the coeff
## @item
##    "s": theoretical dispersion parameter
## @item
##    "z": z statistics for coeff
## @item
##    "pval": p-values for coeff
## @item
##    "resid": raw residuals
## @end itemize
## @end deftypefn

## Original for MATLAB written by Gordon K Smyth <gks@maths.uq.oz.au>,
## U of Queensland, Australia, on Nov 19, 1990.  Last revision Aug 3,
## 1992.

## Author: Gordon K Smyth <gks@maths.uq.oz.au>,
## Adapted-By: KH <Kurt.Hornik@wu-wien.ac.at> and AP <A.C.Penn@sussex.ac.uk>
## Description: Ordinal logistic regression

## Uses the auxiliary functions logistic_regression_derivatives and
## logistic_regression_likelihood.

function [intercept, slope, dev, dl, d2l, P, stats] = logistic_regression (y, x, print, intercept, slope)

  ## check input
  y = round (vec (y));
  missing = (isnan (y) | any (isnan (x), 2));
  y(missing) = [];
  x(missing,:) = [];
  [my, ny] = size (y);
  if (nargin < 2)
    x = zeros (my, 0);
  endif;
  [mx, nx] = size (x);
  if (mx != my)
    error ("logistic_regression: X and Y must have the same number of observations");
  endif

  ## initial calculations
  tol = 1e-12; incr = 10; decr = 2;
  ymin = min (y); ymax = max (y); yrange = ymax - ymin;
  z  = (y * ones (1, yrange)) == ((y * 0 + 1) * (ymin : (ymax - 1)));
  z1 = (y * ones (1, yrange)) == ((y * 0 + 1) * ((ymin + 1) : ymax));
  z  = z(:, any (z));
  z1 = z1(:, any(z1));
  [mz, nz] = size (z);

  ## starting values
  if (nargin < 3)
    print = 0;
  endif;
  if (nargin < 4)
    g = cumsum (sum (z))' ./ my;
    intercept = log (g ./ (1 - g));
  endif;
  if (nargin < 5)
    slope = zeros (nx, 1);
  endif;
  tb = [intercept; slope];

  ## likelihood and derivatives at starting values
  [g, g1, p, dev] = logistic_regression_likelihood (y, x, tb, z, z1);
  [dl, d2l] = logistic_regression_derivatives (x, z, z1, g, g1, p);
  epsilon = std (vec (d2l)) / 1000;

  ## maximize likelihood using Levenberg modified Newton's method
  iter = 0;
  while (abs (dl' * (d2l \ dl) / length (dl)) > tol)
    iter += 1;
    tbold = tb;
    devold = dev;
    tb = tbold - d2l \ dl;
    [g, g1, p, dev] = logistic_regression_likelihood (y, x, tb, z, z1);
    if ((dev - devold) / (dl' * (tb - tbold)) < 0)
      epsilon /= decr;
    else
      while ((dev - devold) / (dl' * (tb - tbold)) > 0)
        epsilon *= incr;
         if (epsilon > 1e+15)
           error ("logistic_regression: epsilon too large");
         endif
         tb = tbold - (d2l - epsilon * eye (size (d2l))) \ dl;
         [g, g1, p, dev] = logistic_regression_likelihood (y, x, tb, z, z1);
         disp ("epsilon"); disp (epsilon);
      endwhile
    endif
    [dl, d2l] = logistic_regression_derivatives (x, z, z1, g, g1, p);
    if (print == 2)
      disp ("Iteration"); disp (iter);
      disp ("Deviance"); disp (dev);
      disp ("First derivative"); disp (dl');
      disp ("Eigenvalues of second derivative"); disp (eig (d2l)');
    endif
  endwhile

  ## tidy up output

  intercept = tb(1 : nz, 1);
  slope  = tb((nz + 1) : (nz + nx), 1);
  cov = inv (-d2l);
  se = sqrt (diag (cov));

  if (nargout > 5)
    ## Compute predicted probabilities (P)
    if (nx > 0)
      e = ((x * slope) * ones (1, nz)) + ((y * 0 + 1) * intercept');
    else
      e = (y * 0 + 1) * intercept';
    endif
    P = diff ([(y * 0), (exp (e) ./ (1 + exp (e))), (y * 0 + 1)]')';
  endif

  if (nargout > 6)
    ## Create stats structure
    dfe = mx - nx - 1;
    zstat = tb ./ se;
    coeffcorr = cov2corr (cov);
    resid = y - P(:,2);
    stats = struct ("intercept", intercept, ...
                    "slope", slope, ...
                    "dfe", dfe, ...
                    "coeff", tb, ...
                    "cov", cov, ...
                    "coeffcorr", coeffcorr, ...
                    "se", se, ...
                    "s", 1, ...
                    "z", zstat, ...
                    "pval", 2 * normcdf (-abs (zstat)), ...
                    "resid", resid);
  endif

  if (print >= 1)
    printf ("\n");
    printf ("Logistic Regression Results:\n");
    printf ("\n");
    printf ("Number of Iterations: %d\n", iter);
    printf ("Deviance:             %f\n", dev);
    printf ("Parameter Estimates:\n");
    printf ("    Intercept      S.E.\n");
    for i = 1 : nz
      printf ("    %8.4f    %8.4f\n", tb (i), se (i));
    endfor
    if (nx > 0)
      printf ("      Slope        S.E.\n");
      for i = (nz + 1) : (nz + nx)
        printf ("    %8.4f    %8.4f\n", tb (i), se (i));
      endfor
    endif
  endif

endfunction

function [g, g1, p, dev] = logistic_regression_likelihood (y, x, slope, z, z1)

  ## Calculate the likelihood for the ordinal logistic regression model.

  e = exp ([z, x] * slope); e1 = exp ([z1, x] * slope);
  g = e ./ (1 + e); g1 = e1 ./ (1 + e1);
  g = max (y == max (y), g); g1 = min (y > min (y), g1);

  p = g - g1;
  dev = -2 * sum (log (p));

endfunction

function [dl, d2l] = logistic_regression_derivatives (x, z, z1, g, g1, p)

  ## Calculate derivatives of the log-likelihood for ordinal logistic regression

  ## first derivative
  v = g .* (1 - g) ./ p; v1 = g1 .* (1 - g1) ./ p;
  dlogp = [(diag (v) * z - diag (v1) * z1), (diag (v - v1) * x)];
  dl = sum (dlogp)';

  ## second derivative
  w = v .* (1 - 2 * g); w1 = v1 .* (1 - 2 * g1);
  d2l = [z, x]' * diag (w) * [z, x] - [z1, x]' * diag (w1) * [z1, x] ...
      - dlogp' * dlogp;

endfunction

function R = cov2corr (vcov)

   ## Convert covariance matrix to correlation matrix
   sed = sqrt (diag (vcov));
   R = vcov ./ (sed * sed');
   R = (R + R') / 2; # This step ensures that the matrix is positive definite

endfunction

%!test
%! # Output compared to following MATLAB commands
%! # [B, DEV, STATS] = mnrfit(X,Y+1,'model','ordinal');
%! # P = mnrval(B,X)
%! X = [1.489381332449196, 1.1534152241851305; ...
%!      1.8110085304863965, 0.9449666896938425; ...
%!      -0.04453299665130296, 0.34278203449678646; ...
%!      -0.36616019468850347, 1.130254275908322; ...
%!       0.15339143291005095, -0.7921044310668951; ...
%!      -1.6031878794469698, -1.8343471035233376; ...
%!      -0.14349521143198166, -0.6762996896828459; ...
%!      -0.4403818557740143, -0.7921044310668951; ...
%!      -0.7372685001160434, -0.027793137932169563; ...
%!      -0.11875465773681024, 0.5512305689880763];
%! Y = [1,1,1,1,1,0,0,0,0,0]';
%! [INTERCEPT, SLOPE, DEV, DL, D2L, P] = logistic_regression (Y, X, false);
#%! assert (DEV, 5.680728861124, 1e-05);
#%! assert (INTERCEPT(1), -1.10999599948243, 1e-05);
#%! assert (SLOPE(1), -9.12480634225699, 1e-05);
#%! assert (SLOPE(2), -2.18746124517476, 1e-05);
#%! assert (corr(P(:,1),Y), -0.786673288976468, 1e-05);

%!test
%! # Output compared to following MATLAB commands
%! # [B, DEV, STATS] = mnrfit(X,Y+1,'model','ordinal');
%! load carbig
%! X = [Acceleration Displacement Horsepower Weight];
%! miles = [1,1,1,1,1,1,1,1,1,1,NaN,NaN,NaN,NaN,NaN,1,1,NaN,1,1,2,2,1,2,2,2, ...
%!          2,2,2,2,2,1,1,1,1,2,2,2,2,NaN,2,1,1,2,1,1,1,1,1,1,1,1,1,2,2,1,2, ...
%!          2,3,3,3,3,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,2,2,2,2,2, ...
%!          2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,1,1,1,1,1,2,2,2,1,2,2, ...
%!          2,1,1,3,2,2,2,1,2,2,1,2,2,2,1,3,2,3,2,1,1,1,1,1,1,1,1,3,2,2,3,3, ...
%!          2,2,2,2,2,3,2,1,1,1,1,1,1,1,1,1,1,1,2,2,1,3,2,2,2,2,2,2,1,3,2,2, ...
%!          2,2,2,3,2,2,2,2,2,1,1,1,1,2,2,2,2,3,2,3,3,2,1,1,1,3,3,2,2,2,1,2, ...
%!          2,1,1,1,1,1,3,3,3,2,3,1,1,1,1,1,2,2,1,1,1,1,1,3,2,2,2,3,3,3,3,2, ...
%!          2,2,4,3,3,4,3,2,2,2,2,2,2,2,2,2,2,2,1,1,2,1,1,1,3,2,2,3,2,2,2,2, ...
%!          2,1,2,1,3,3,2,2,2,2,2,1,1,1,1,1,1,2,1,3,3,3,2,2,2,2,2,3,3,3,3,2, ...
%!          2,2,3,4,3,3,3,2,2,2,2,3,3,3,3,3,4,2,4,4,4,3,3,4,4,3,3,3,2,3,2,3, ...
%!          2,2,2,2,3,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,NaN,3,2,2,2,2,2,1,2, ...
%!          2,3,3,3,2,2,2,3,3,3,3,3,3,3,3,3,3,3,2,3,2,2,3,3,2,2,4,3,2,3]';
%! [INTERCEPT, SLOPE, DEV, DL, D2L, P] = logistic_regression (miles, X, false);
%! assert (DEV, 433.197174495549, 1e-05);
%! assert (INTERCEPT(1), -16.6895155618903, 1e-05);
%! assert (INTERCEPT(2), -11.7207818178493, 1e-05);
%! assert (INTERCEPT(3), -8.0605768506075, 1e-05);
%! assert (SLOPE(1), 0.104762463756714, 1e-05);
%! assert (SLOPE(2), 0.0103357623191891, 1e-05);
%! assert (SLOPE(3), 0.0645199313242276, 1e-05);
%! assert (SLOPE(4), 0.00166377028388103, 1e-05);