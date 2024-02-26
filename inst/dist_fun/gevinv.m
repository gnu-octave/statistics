## Copyright (C) 2012 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## @deftypefn  {statistics} {@var{x} =} gevinv (@var{p}, @var{k}, @var{sigma}, @var{mu})
##
## Inverse of the generalized extreme value (GEV) cumulative distribution
## function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the GEV distribution with shape parameter @var{k}, scale parameter
## @var{sigma}, and location parameter @var{mu}.  The size of @var{p} is the
## common size of @var{x}, @var{k}, @var{sigma}, and @var{mu}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## When @qcode{@var{k} < 0}, the GEV is the type III extreme value distribution.
## When @qcode{@var{k} > 0}, the GEV distribution is the type II, or Frechet,
## extreme value distribution.  If @var{W} has a Weibull distribution as
## computed by the @code{wblcdf} function, then @qcode{-@var{W}} has a type III
## extreme value distribution and @qcode{1/@var{W}} has a type II extreme value
## distribution.  In the limit as @var{k} approaches @qcode{0}, the GEV is the
## mirror image of the type I extreme value distribution as computed by the
## @code{evcdf} function.
##
## The mean of the GEV distribution is not finite when @qcode{@var{k} >= 1}, and
## the variance is not finite when @qcode{@var{k} >= 1/2}.  The GEV distribution
## has positive density only for values of @var{x} such that
## @qcode{@var{k} * (@var{x} - @var{mu}) / @var{sigma} > -1}.
##
## Further information about the generalized extreme value distribution can be
## found at
## @url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}
##
## @seealso{gevcdf, gevpdf, gevrnd, gevfit, gevlike, gevstat}
## @end deftypefn

function x = gevinv (p, k, sigma, mu)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("gevinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, K, SIGMA, and MU
  [retval, p, k, sigma, mu] = common_size (p, k, sigma, mu);
  if (retval > 0)
    error ("gevinv: P, K, SIGMA, and MU must be of common size or scalars.");
  endif

  ## Check for P, K, SIGMA, and MU being reals
  if (iscomplex (p) || iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gevinv: P, K, SIGMA, and MU must not be complex.");
  endif

  is_neginf = p == 0;
  is_posinf = p == 1;
  is_nan = p < 0 | p > 1 | isnan (p);
  x = p;

  llP = log (-log (p));
  kllP = k .* llP;

  ## Use the Taylor series expansion of the exponential to
  ## avoid roundoff error or dividing by zero when k is small
  ii = (abs(kllP) < 1E-4);
  x(ii) = mu(ii) - sigma(ii) .* llP(ii) .* (1 - kllP(ii) .* (1 - kllP(ii)));
  x(~ii) = mu(~ii) + (sigma(~ii) ./ k(~ii)) .* (exp(-kllP(~ii)) - 1);
  x(is_neginf) = -Inf;
  x(is_posinf) = Inf;
  x(is_nan) = NaN;
endfunction

%!demo
%! ## Plot various iCDFs from the generalized extreme value distribution
%! p = 0.001:0.001:0.999;
%! x1 = gevinv (p, 1, 1, 1);
%! x2 = gevinv (p, 0.5, 1, 1);
%! x3 = gevinv (p, 1, 1, 5);
%! x4 = gevinv (p, 1, 2, 5);
%! x5 = gevinv (p, 1, 5, 5);
%! x6 = gevinv (p, 1, 0.5, 5);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", ...
%!       p, x4, "-c", p, x5, "-m", p, x6, "-k")
%! grid on
%! ylim ([-1, 10])
%! legend ({"k = 1, σ = 1, μ = 1", "k = 0.5, σ = 1, μ = 1", ...
%!          "k = 1, σ = 1, μ = 5", "k = 1, σ = 2, μ = 5", ...
%!          "k = 1, σ = 5, μ = 5", "k = 1, σ = 0.5, μ = 5"}, ...
%!         "location", "northwest")
%! title ("Generalized extreme value iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!test
%! p = 0.1:0.1:0.9;
%! k = 0;
%! sigma = 1;
%! mu = 0;
%! x = gevinv (p, k, sigma, mu);
%! c = gevcdf(x, k, sigma, mu);
%! assert (c, p, 0.001);
%!test
%! p = 0.1:0.1:0.9;
%! k = 1;
%! sigma = 1;
%! mu = 0;
%! x = gevinv (p, k, sigma, mu);
%! c = gevcdf(x, k, sigma, mu);
%! assert (c, p, 0.001);
%!test
%! p = 0.1:0.1:0.9;
%! k = 0.3;
%! sigma = 1;
%! mu = 0;
%! x = gevinv (p, k, sigma, mu);
%! c = gevcdf(x, k, sigma, mu);
%! assert (c, p, 0.001);

## Test input validation
%!error<gevinv: function called with too few input arguments.> gevinv ()
%!error<gevinv: function called with too few input arguments.> gevinv (1)
%!error<gevinv: function called with too few input arguments.> gevinv (1, 2)
%!error<gevinv: function called with too few input arguments.> gevinv (1, 2, 3)
%!error<gevinv: P, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevinv (ones (3), ones (2), ones(2), ones(2))
%!error<gevinv: P, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevinv (ones (2), ones (3), ones(2), ones(2))
%!error<gevinv: P, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevinv (ones (2), ones (2), ones(3), ones(2))
%!error<gevinv: P, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevinv (ones (2), ones (2), ones(2), ones(3))
%!error<gevinv: P, K, SIGMA, and MU must not be complex.> gevinv (i, 2, 3, 4)
%!error<gevinv: P, K, SIGMA, and MU must not be complex.> gevinv (1, i, 3, 4)
%!error<gevinv: P, K, SIGMA, and MU must not be complex.> gevinv (1, 2, i, 4)
%!error<gevinv: P, K, SIGMA, and MU must not be complex.> gevinv (1, 2, 3, i)
