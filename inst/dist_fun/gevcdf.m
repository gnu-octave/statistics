## Copyright (C) 2012 Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{p} =} gevcdf (@var{x}, @var{k}, @var{sigma}, @var{mu})
## @deftypefnx {statistics} {@var{p} =} gevcdf (@var{x}, @var{k}, @var{sigma}, @var{mu}, @qcode{"upper"})
##
## Generalized extreme value (GEV) cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the GEV distribution with shape parameter @var{k}, scale parameter
## @var{sigma}, and location parameter @var{mu}.  The size of @var{p} is the
## common size of @var{x}, @var{k}, @var{sigma}, and @var{mu}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## @code{[@dots{}] = gevcdf (@var{x}, @var{k}, @var{sigma}, @var{mu}, "upper")}
## computes the upper tail probability of the GEV distribution with parameters
## @var{k}, @var{sigma}, and @var{mu}, at the values in @var{x}.
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
## @subheading References
## @enumerate
## @item
## Rolf-Dieter Reiss and Michael Thomas. @cite{Statistical Analysis of Extreme
## Values with Applications to Insurance, Finance, Hydrology and Other Fields}.
## Chapter 1, pages 16-17, Springer, 2007.
## @end enumerate
##
## @seealso{gevinv, gevpdf, gevrnd, gevfit, gevlike, gevstat}
## @end deftypefn

function p = gevcdf (x, k, sigma, mu, uflag)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("gevcdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 4)
    if (! strcmpi (uflag, "upper"))
      error ("gevcdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, K, SIGMA, and MU
  if (! isscalar (x) || ! isscalar (k) || ! isscalar (sigma) || ! isscalar (mu))
    [err, x, k, sigma, mu] = common_size (x, k, sigma, mu);
    if (err > 0)
      error ("gevcdf: X, K, SIGMA, and MU must be of common size or scalars.");
    endif
  endif

  ## Check for X, K, SIGMA, and MU being reals
  if (iscomplex (x) || iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gevcdf: X, K, SIGMA, and MU must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (k, "single") ...
                        || isa (sigma, "single") || isa (mu, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Prepare output
  p = zeros (size (x), is_class);

  ## Return NaN for out of range parameter SIGMA.
  sigma(sigma <= 0) = NaN;

  ## Calculate z
  z = (x - mu) ./ sigma;

  ## Process k == 0
  k_0 = (abs(k) < eps);
  if (uflag)
    p(k_0) = -expm1 (-exp (-z(k_0)));
  else
    p(k_0) = exp (-exp (-z(k_0)));
  endif

  ## Process k != 0
  k_0 = ! k_0;
  t = z .* k;
  if (uflag)
    p(k_0) = -expm1 (-exp (-(1 ./ k(k_0)) .* log1p (t(k_0))));
  else
    p(k_0) = exp (-exp (-(1 ./ k(k_0)) .* log1p (t(k_0))));
  endif

  ## Return 0 or 1 for 1 + k.*(x-mu)/sigma > 0
  k_1 = k_0 & (t<=-1);
  t(k_1) = 0;
  if uflag == true
    p(k_1) = (k(k_1) >= 0);
  else
    p(k_1) = (k(k_1) < 0);
  endif

endfunction

%!demo
%! ## Plot various CDFs from the generalized extreme value distribution
%! x = -1:0.001:10;
%! p1 = gevcdf (x, 1, 1, 1);
%! p2 = gevcdf (x, 0.5, 1, 1);
%! p3 = gevcdf (x, 1, 1, 5);
%! p4 = gevcdf (x, 1, 2, 5);
%! p5 = gevcdf (x, 1, 5, 5);
%! p6 = gevcdf (x, 1, 0.5, 5);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", ...
%!       x, p4, "-c", x, p5, "-m", x, p6, "-k")
%! grid on
%! xlim ([-1, 10])
%! legend ({"ξ = 1, σ = 1, μ = 1", "ξ = 0.5, σ = 1, μ = 1", ...
%!          "ξ = 1, σ = 1, μ = 5", "ξ = 1, σ = 2, μ = 5", ...
%!          "ξ = 1, σ = 5, μ = 5", "ξ = 1, σ = 0.5, μ = 5"}, ...
%!         "location", "southeast")
%! title ("Generalized extreme value CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!test
%! x = 0:0.5:2.5;
%! sigma = 1:6;
%! k = 1;
%! mu = 0;
%! p = gevcdf (x, k, sigma, mu);
%! expected_p = [0.36788, 0.44933, 0.47237, 0.48323, 0.48954, 0.49367];
%! assert (p, expected_p, 0.001);
%!test
%! x = -0.5:0.5:2.5;
%! sigma = 0.5;
%! k = 1;
%! mu = 0;
%! p = gevcdf (x, k, sigma, mu);
%! expected_p = [0, 0.36788, 0.60653, 0.71653, 0.77880, 0.81873, 0.84648];
%! assert (p, expected_p, 0.001);
%!test # check for continuity for k near 0
%! x = 1;
%! sigma = 0.5;
%! k = -0.03:0.01:0.03;
%! mu = 0;
%! p = gevcdf (x, k, sigma, mu);
%! expected_p = [0.88062, 0.87820, 0.87580, 0.87342, 0.87107, 0.86874, 0.86643];
%! assert (p, expected_p, 0.001);

## Test input validation
%!error<gevcdf: function called with too few input arguments.> gevcdf ()
%!error<gevcdf: function called with too few input arguments.> gevcdf (1)
%!error<gevcdf: function called with too few input arguments.> gevcdf (1, 2)
%!error<gevcdf: function called with too few input arguments.> gevcdf (1, 2, 3)
%!error<gevcdf: function called with too many inputs> ...
%! gevcdf (1, 2, 3, 4, 5, 6)
%!error<gevcdf: invalid argument for upper tail.> gevcdf (1, 2, 3, 4, "tail")
%!error<gevcdf: invalid argument for upper tail.> gevcdf (1, 2, 3, 4, 5)
%!error<gevcdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevcdf (ones (3), ones (2), ones(2), ones(2))
%!error<gevcdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevcdf (ones (2), ones (3), ones(2), ones(2))
%!error<gevcdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevcdf (ones (2), ones (2), ones(3), ones(2))
%!error<gevcdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevcdf (ones (2), ones (2), ones(2), ones(3))
%!error<gevcdf: X, K, SIGMA, and MU must not be complex.> gevcdf (i, 2, 3, 4)
%!error<gevcdf: X, K, SIGMA, and MU must not be complex.> gevcdf (1, i, 3, 4)
%!error<gevcdf: X, K, SIGMA, and MU must not be complex.> gevcdf (1, 2, i, 4)
%!error<gevcdf: X, K, SIGMA, and MU must not be complex.> gevcdf (1, 2, 3, i)
