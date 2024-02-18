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
## @deftypefn  {statistics} {@var{x} =} invginv (@var{p}, @var{mu}, @var{lambda})
##
## Inverse of the inverse Gaussian cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the inverse Gaussian distribution with scale parameter @var{mu} and shape
## parameter @var{lambda}.  The size of @var{x} is the common size of @var{p},
## @var{mu}, and @var{lambda}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## The inverse Gaussian CDF is only defined for @qcode{@var{mu} > 0} and
## @qcode{@var{lambda} > 0}.
##
## Further information about the inverse Gaussian distribution can be found at
## @url{https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}
##
## @seealso{invgcdf, invgpdf, invgrnd, invgfit, invglike}
## @end deftypefn

function x = invginv (p, mu, lambda)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("invginv: function called with too few input arguments.");
  endif

  ## Check for common size of P, MU, and LAMBDA
  if (! isscalar (p) || ! isscalar (mu) || ! isscalar(lambda))
    [retval, p, mu, lambda] = common_size (p, mu, lambda);
    vec = true;
    if (retval > 0)
      error ("invginv: P, MU, and LAMBDA must be of common size or scalars.");
    endif
  else
    vec = false;
  endif

  ## Check for X, MU, and LAMBDA being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (lambda))
    error ("invginv: P, MU, and LAMBDA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (mu, "single") || isa (lambda, "single"));
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Return NaNs for out of range values of MU and LAMBDA parameters
  mu(mu <= 0) = NaN;
  lambda(lambda <= 0) = NaN;

  ## Find valid parameters and p-values (handle edges cases below)
  validmulambda = (mu > 0) & (lambda > 0) & (lambda < Inf);
  validp_values = (validmulambda & (p > 0) & (p < 1));
  valid_all = all (validp_values(:));
  valid_any = any (validp_values(:));

  ## Handle edges cases here
  if (! valid_all)
    x(p == 0 & validmulambda) = 0;
    x(p == 1 & validmulambda) = Inf;

    ## Keep valid cases (if any left)
    if (valid_any)
      if (vec)
        p = p(validp_values);
        mu = mu(validp_values);
        lambda = lambda(validp_values);
      endif
    else
      return;
    endif
  endif

  ## Apply Newton's Method to find a root of p = invgcdf(x,mu,lambda)
  ## Choose a starting guess for x0.  Use quantiles from a lognormal
  ## distribution with the same mean (==1) and variance (==lambda0)
  lambda0 = lambda ./ mu;
  lognorm = log (1 ./ lambda0 + 1);
  mulnorm = -0.5 .* lognorm;
  x0 = exp (mulnorm - sqrt (2 .* lognorm) .* erfcinv (2 * p));

  ## Set maximum iterations and tolerance for Newton's Method
  mit = 500;
  tol = eps (class (x0)) .^ (3/4);

  ## Get quantiles
  F = invgcdf (x0, 1, lambda0);
  dF = F - p;
  for it = 1:mit
    ## Compute the Newton step
    f = invgpdf (x0, 1, lambda0);
    h = dF ./ f;
    x0_1 = max (x0/10, min (10 * x0, x0 - h));

    ## Check if tolerance is reached
    complete = (abs (h) <= tol * x0);
    if (all (complete(:)))
      x0 = x0_1;
      break
    endif

    ## Check for increasing error unless tolerance is reached
    dFold = dF;
    for j = 1:25
      F = invgcdf (x0_1, 1, lambda0);
      dF = F - p;
      worse = (abs (dF) > abs (dFold)) & ! complete;
      if (! any (worse(:)))
        break
      endif
      x0_1(worse) = (x0(worse) + x0_1(worse)) / 2;
    endfor

    ## Update for next step
    x0 = x0_1;
  endfor

  ## Issue a warning for exceeding iterations or not converging to tolerance
  notconv = (abs(dF./F) > tol.^(2/3));
  if (it > mit || any (notconv(:)))
    warning (strcat (["invginv: Newton's Method did not converge"], ...
                     [" or exceeded maximum iterations."]));
  endif

  ## Apply the scale factor
  if (valid_all)
    x = x0 .* mu;
  else
    x(validp_values) = x0 .* mu;
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the inverse Gaussian distribution
%! p = 0.001:0.001:0.999;
%! x1 = invginv (p, 1, 0.2);
%! x2 = invginv (p, 1, 1);
%! x3 = invginv (p, 1, 3);
%! x4 = invginv (p, 3, 0.2);
%! x5 = invginv (p, 3, 1);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c", p, x5, "-y")
%! grid on
%! ylim ([0, 3])
%! legend ({"μ = 1, σ = 0.2", "μ = 1, σ = 1", "μ = 1, σ = 3", ...
%!          "μ = 3, σ = 0.2", "μ = 3, σ = 1"}, "location", "northwest")
%! title ("Inverse Gaussian iCDF")
%! xlabel ("probability")
%! ylabel ("x")

## Test output
%!shared p, x
%! p = [0, 0.3829, 0.6827, 1];
%! x = [0, 0.5207, 1.0376, Inf];
%!assert (invginv (p, 1, 1), x, 1e-4);
%!assert (invginv (p, 1, ones (1,4)), x, 1e-4);
%!assert (invginv (p, 1, [-1, 0, 1, 1]), [NaN, NaN, x(3:4)], 1e-4)
%!assert (invginv (p, [-1, 0, 1, 1], 1), [NaN, NaN, x(3:4)], 1e-4)

## Test class of input preserved
%!assert (class (invginv (single ([p, NaN]), 0, 1)), "single")
%!assert (class (invginv ([p, NaN], single (0), 1)), "single")
%!assert (class (invginv ([p, NaN], 0, single (1))), "single")

## Test input validation
%!error<invginv: function called with too few input arguments.> invginv (1)
%!error<invginv: function called with too few input arguments.> invginv (1, 2)
%!error<invginv: P, MU, and LAMBDA must be of common size or scalars.> ...
%! invginv (1, ones (2), ones (3))
%!error<invginv: P, MU, and LAMBDA must be of common size or scalars.> ...
%! invginv (ones (2), 1, ones (3))
%!error<invginv: P, MU, and LAMBDA must be of common size or scalars.> ...
%! invginv (ones (2), ones (3), 1)
%!error<invginv: P, MU, and LAMBDA must not be complex.> invginv (i, 2, 3)
%!error<invginv: P, MU, and LAMBDA must not be complex.> invginv (1, i, 3)
%!error<invginv: P, MU, and LAMBDA must not be complex.> invginv (1, 2, i)
