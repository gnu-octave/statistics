## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{y} =} gampdf (@var{x}, @var{k}, @var{theta})
##
## Gamma probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Gamma distribution with shape parameter @var{k} and scale parameter
## @var{theta}.  The size of @var{y} is the common size of @var{x}, @var{k} and
## @var{theta}.  A scalar input functions as a constant matrix of the same size
## as the other inputs.
##
## There are two equivalent parameterizations in common use:
## @enumerate
## @item With a shape parameter @math{k} and a scale parameter @math{θ}, which
## is used by @code{gampdf}.
## @item With a shape parameter @math{α = k} and an inverse scale parameter
## @math{β = 1 / θ}, called a rate parameter.
## @end enumerate
##
## Further information about the Gamma distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gamma_distribution}
##
## @seealso{gamcdf, gaminv, gamrnd, gamfit, gamlike, gamstat}
## @end deftypefn

function y = gampdf (x, k, theta)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("gampdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, K, and THETA
  if (! isscalar (k) || ! isscalar (theta))
    [retval, x, k, theta] = common_size (x, k, theta);
    if (retval > 0)
      error ("gampdf: X, K, and THETA must be of common size or scalars.");
    endif
  endif

  ## Check for X, K, and THETA being reals
  if (iscomplex (x) || iscomplex (k) || iscomplex (theta))
    error ("gampdf: X, K, and THETA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (k, "single") || isa (theta, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Force NaNs for out of range parameters
  is_nan = ! (k > 0) | ! (theta > 0) | isnan (x);
  y(is_nan) = NaN;

  ## Handle all other valid cases
  v = (x >= 0) & (k > 0) & (k <= 1) & (theta > 0);
  if (isscalar (k) && isscalar (theta))
    y(v) = (x(v) .^ (k - 1)) ...
              .* exp (- x(v) / theta) / gamma (k) / (theta ^ k);
  else
    y(v) = (x(v) .^ (k(v) - 1)) ...
              .* exp (- x(v) ./ theta(v)) ./ gamma (k(v)) ./ (theta(v) .^ k(v));
  endif

  v = (x >= 0) & (k > 1) & (theta > 0);
  if (isscalar (k) && isscalar (theta))
    y(v) = exp (- k * log (theta) + (k-1) * log (x(v))
                  - x(v) / theta - gammaln (k));
  else
    y(v) = exp (- k(v) .* log (theta(v)) + (k(v)-1) .* log (x(v))
                  - x(v) ./ theta(v) - gammaln (k(v)));
  endif

endfunction

%!demo
%! ## Plot various PDFs from the Gamma distribution
%! x = 0:0.01:20;
%! y1 = gampdf (x, 1, 2);
%! y2 = gampdf (x, 2, 2);
%! y3 = gampdf (x, 3, 2);
%! y4 = gampdf (x, 5, 1);
%! y5 = gampdf (x, 9, 0.5);
%! y6 = gampdf (x, 7.5, 1);
%! y7 = gampdf (x, 0.5, 1);
%! plot (x, y1, "-r", x, y2, "-g", x, y3, "-y", x, y4, "-m", ...
%!       x, y5, "-k", x, y6, "-b", x, y7, "-c")
%! grid on
%! ylim ([0,0.5])
%! legend ({"α = 1, θ = 2", "α = 2, θ = 2", "α = 3, θ = 2", ...
%!          "α = 5, θ = 1", "α = 9, θ = 0.5", "α = 7.5, θ = 1", ...
%!          "α = 0.5, θ = 1"}, "location", "northeast")
%! title ("Gamma PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 0.5 1 Inf];
%! y = [0 exp(-x(2:end))];
%!assert (gampdf (x, ones (1,5), ones (1,5)), y)
%!assert (gampdf (x, 1, ones (1,5)), y)
%!assert (gampdf (x, ones (1,5), 1), y)
%!assert (gampdf (x, [0 -Inf NaN Inf 1], 1), [NaN NaN NaN NaN y(5)])
%!assert (gampdf (x, 1, [0 -Inf NaN Inf 1]), [NaN NaN NaN 0 y(5)])
%!assert (gampdf ([x, NaN], 1, 1), [y, NaN])

## Test class of input preserved
%!assert (gampdf (single ([x, NaN]), 1, 1), single ([y, NaN]))
%!assert (gampdf ([x, NaN], single (1), 1), single ([y, NaN]))
%!assert (gampdf ([x, NaN], 1, single (1)), single ([y, NaN]))

## Test input validation
%!error<gampdf: function called with too few input arguments.> gampdf ()
%!error<gampdf: function called with too few input arguments.> gampdf (1)
%!error<gampdf: function called with too few input arguments.> gampdf (1,2)
%!error<gampdf: X, K, and THETA must be of common size or scalars.> ...
%! gampdf (ones (3), ones (2), ones (2))
%!error<gampdf: X, K, and THETA must be of common size or scalars.> ...
%! gampdf (ones (2), ones (3), ones (2))
%!error<gampdf: X, K, and THETA must be of common size or scalars.> ...
%! gampdf (ones (2), ones (2), ones (3))
%!error<gampdf: X, K, and THETA must not be complex.> gampdf (i, 2, 2)
%!error<gampdf: X, K, and THETA must not be complex.> gampdf (2, i, 2)
%!error<gampdf: X, K, and THETA must not be complex.> gampdf (2, 2, i)
