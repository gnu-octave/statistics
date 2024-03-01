## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{x} =} gaminv (@var{p}, @var{a}, @var{b})
##
## Inverse of the Gamma cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the Gamma distribution with shape parameter @var{a} and scale parameter
## @var{b}.  The size of @var{x} is the common size of @var{p}, @var{a},
## and @var{b}.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## OCTAVE/MATLAB use the alternative parameterization given by the pair
## @math{α, β}, i.e. shape @var{a} and scale @var{b}.  In Wikipedia, the two
## common parameterizations use the pairs @math{k, θ}, as shape and scale, and
## @math{α, β}, as shape and rate, respectively.  The parameter names @var{a}
## and @var{b} used here (for MATLAB compatibility) correspond to the parameter
## notation @math{k, θ} instead of the @math{α, β} as reported in Wikipedia.
##
## Further information about the Gamma distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gamma_distribution}
##
## @seealso{gamcdf, gampdf, gamrnd, gamfit, gamlike, gamstat}
## @end deftypefn

function x = gaminv (p, a, b)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("gaminv: function called with too few input arguments.");
  endif

  ## Check for common size of P, Α, and Β
  if (! isscalar (p) || ! isscalar (a) || ! isscalar (b))
    [retval, p, a, b] = common_size (p, a, b);
    if (retval > 0)
      error ("gaminv: P, Α, and Β must be of common size or scalars.");
    endif
  endif

  ## Check for P, Α, and Β being reals
  if (iscomplex (p) || iscomplex (a) || iscomplex (b))
    error ("gaminv: P, Α, and Β must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (a, "single") || isa (b, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  ## Force NaNs for out of range parameters
  is_nan = ((p < 0) | (p > 1) | isnan (p) ...
        | ! (a > 0) | ! (a < Inf) | ! (b > 0) | ! (b < Inf));
  x(is_nan) = NaN;

  ## Handle edge cases
  is_inf = (p == 1) & (a > 0) & (a < Inf) & (b > 0) & (b < Inf);
  x(is_inf) = Inf;

  ## Handle all other valid cases
  is_valid = find ((p > 0) & (p < 1) & (a > 0) & ...
                   (a < Inf) & (b > 0) & (b < Inf));
  if (! isempty (is_valid))
    if (! isscalar (a) || ! isscalar (b))
      a = a(is_valid);
      b = b(is_valid);
      y = a .* b;
    else
      y = a * b * ones (size (is_valid));
    endif
    p = p(is_valid);

    ## Call GAMMAINCINV to find a root of GAMMAINC
    q = gammaincinv (p, a);
    tol = sqrt (eps (ones (1, 1, class(q))));
    check_cdf = ((abs (gammainc (q, a) - p) ./ p) > tol);
    ## Check for any cdf being far off from tolerance
    if (any (check_cdf(:)))
      warning ("gaminv: calculation failed to converge for some values.");
    endif
    x(is_valid) = q .* b;
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the Gamma distribution
%! p = 0.001:0.001:0.999;
%! x1 = gaminv (p, 1, 2);
%! x2 = gaminv (p, 2, 2);
%! x3 = gaminv (p, 3, 2);
%! x4 = gaminv (p, 5, 1);
%! x5 = gaminv (p, 9, 0.5);
%! x6 = gaminv (p, 7.5, 1);
%! x7 = gaminv (p, 0.5, 1);
%! plot (p, x1, "-r", p, x2, "-g", p, x3, "-y", p, x4, "-m", ...
%!       p, x5, "-k", p, x6, "-b", p, x7, "-c")
%! ylim ([0, 20])
%! grid on
%! legend ({"α = 1, β = 2", "α = 2, β = 2", "α = 3, β = 2", ...
%!          "α = 5, β = 1", "α = 9, β = 0.5", "α = 7.5, β = 1", ...
%!          "α = 0.5, β = 1"}, "location", "northwest")
%! title ("Gamma iCDF")
%! xlabel ("probability")
%! ylabel ("x")

## Test output
%!shared p
%! p = [-1 0 0.63212055882855778 1 2];
%!assert (gaminv (p, ones (1,5), ones (1,5)), [NaN 0 1 Inf NaN], eps)
%!assert (gaminv (p, 1, ones (1,5)), [NaN 0 1 Inf NaN], eps)
%!assert (gaminv (p, ones (1,5), 1), [NaN 0 1 Inf NaN], eps)
%!assert (gaminv (p, [1 -Inf NaN Inf 1], 1), [NaN NaN NaN NaN NaN])
%!assert (gaminv (p, 1, [1 -Inf NaN Inf 1]), [NaN NaN NaN NaN NaN])
%!assert (gaminv ([p(1:2) NaN p(4:5)], 1, 1), [NaN 0 NaN Inf NaN])
%!assert (gaminv ([p(1:2) NaN p(4:5)], 1, 1), [NaN 0 NaN Inf NaN])

## Test for accuracy when p is small. Results compared to Matlab
%!assert (gaminv (1e-16, 1, 1), 1e-16, eps)
%!assert (gaminv (1e-16, 1, 2), 2e-16, eps)
%!assert (gaminv (1e-20, 3, 5), 1.957434012161815e-06, eps)
%!assert (gaminv (1e-15, 1, 1), 1e-15, eps)
%!assert (gaminv (1e-35, 1, 1), 1e-35, eps)

## Test class of input preserved
%!assert (gaminv ([p, NaN], 1, 1), [NaN 0 1 Inf NaN NaN], eps)
%!assert (gaminv (single ([p, NaN]), 1, 1), single ([NaN 0 1 Inf NaN NaN]), ...
%! eps ("single"))
%!assert (gaminv ([p, NaN], single (1), 1), single ([NaN 0 1 Inf NaN NaN]), ...
%! eps ("single"))
%!assert (gaminv ([p, NaN], 1, single (1)), single ([NaN 0 1 Inf NaN NaN]), ...
%! eps ("single"))

## Test input validation
%!error<gaminv: function called with too few input arguments.> gaminv ()
%!error<gaminv: function called with too few input arguments.> gaminv (1)
%!error<gaminv: function called with too few input arguments.> gaminv (1,2)
%!error<gaminv: P, Α, and Β must be of common size or scalars.> ...
%! gaminv (ones (3), ones (2), ones (2))
%!error<gaminv: P, Α, and Β must be of common size or scalars.> ...
%! gaminv (ones (2), ones (3), ones (2))
%!error<gaminv: P, Α, and Β must be of common size or scalars.> ...
%! gaminv (ones (2), ones (2), ones (3))
%!error<gaminv: P, Α, and Β must not be complex.> gaminv (i, 2, 2)
%!error<gaminv: P, Α, and Β must not be complex.> gaminv (2, i, 2)
%!error<gaminv: P, Α, and Β must not be complex.> gaminv (2, 2, i)
