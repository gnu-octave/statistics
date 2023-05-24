## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
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
## @deftypefn  {statistics} {@var{x} =} gaminv (@var{p}, @var{k}, @var{theta})
##
## Inverse of the Gamma cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the Gamma distribution with shape parameter @var{k} and scale parameter
## @var{theta}.  The size of @var{x} is the common size of @var{p}, @var{k},
## and @var{theta}.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## There are two equivalent parameterizations in common use:
## @enumerate
## @item With a shape parameter @math{k} and a scale parameter @math{θ}, which
## is used by @code{gaminv}.
## @item With a shape parameter @math{α = k} and an inverse scale parameter
## @math{β = 1 / θ}, called a rate parameter.
## @end enumerate
##
## Further information about the Gamma distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gamma_distribution}
##
## @seealso{gamcdf, gampdf, gamrnd, gamfit, gamlike, gamstat}
## @end deftypefn

function x = gaminv (p, k, theta)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("gaminv: function called with too few input arguments.");
  endif

  ## Check for common size of P, K, and THETA
  if (! isscalar (p) || ! isscalar (k) || ! isscalar (theta))
    [retval, p, k, theta] = common_size (p, k, theta);
    if (retval > 0)
      error ("gaminv: P, K, and THETA must be of common size or scalars.");
    endif
  endif

  ## Check for P, K, and THETA being reals
  if (iscomplex (p) || iscomplex (k) || iscomplex (theta))
    error ("gaminv: P, K, and THETA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (k, "single") || isa (theta, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  ## Force NaNs for out of range parameters
  is_nan = ((p < 0) | (p > 1) | isnan (p) ...
        | ! (k > 0) | ! (k < Inf) | ! (theta > 0) | ! (theta < Inf));
  x(is_nan) = NaN;

  ## Handle edge cases
  is_inf = (p == 1) & (k > 0) & (k < Inf) & (theta > 0) & (theta < Inf);
  x(is_inf) = Inf;

  ## Handle all other valid cases
  is_valid = find ((p > 0) & (p < 1) & (k > 0) & ...
                   (k < Inf) & (theta > 0) & (theta < Inf));
  if (! isempty (is_valid))
    if (! isscalar (k) || ! isscalar (theta))
      k = k(is_valid);
      theta = theta(is_valid);
      y = k .* theta;
    else
      y = k * theta * ones (size (is_valid));
    endif
    p = p(is_valid);

    ## Call GAMMAINCINV to find k root of GAMMAINC
    q = gammaincinv (p, k);
    tol = sqrt (eps (ones (1, 1, class(q))));
    check_cdf = ((abs (gammainc (q, k) - p) ./ p) > tol);
    ## Check for any cdf being far off from tolerance
    if (any (check_cdf(:)))
      warning ("gaminv: calculation failed to converge for some values.");
    endif
    x(is_valid) = q .* theta;
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
%! legend ({"α = 1, θ = 2", "α = 2, θ = 2", "α = 3, θ = 2", ...
%!          "α = 5, θ = 1", "α = 9, θ = 0.5", "α = 7.5, θ = 1", ...
%!          "α = 0.5, θ = 1"}, "location", "northwest")
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
%!error<gaminv: P, K, and THETA must be of common size or scalars.> ...
%! gaminv (ones (3), ones (2), ones (2))
%!error<gaminv: P, K, and THETA must be of common size or scalars.> ...
%! gaminv (ones (2), ones (3), ones (2))
%!error<gaminv: P, K, and THETA must be of common size or scalars.> ...
%! gaminv (ones (2), ones (2), ones (3))
%!error<gaminv: P, K, and THETA must not be complex.> gaminv (i, 2, 2)
%!error<gaminv: P, K, and THETA must not be complex.> gaminv (2, i, 2)
%!error<gaminv: P, K, and THETA must not be complex.> gaminv (2, 2, i)
