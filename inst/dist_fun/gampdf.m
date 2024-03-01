## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{y} =} gampdf (@var{x}, @var{a}, @var{b})
##
## Gamma probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Gamma distribution with shape parameter @var{a} and scale parameter
## @var{b}.  The size of @var{y} is the common size of @var{x}, @var{a} and
## @var{b}.  A scalar input functions as a constant matrix of the same size
## as the other inputs.
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
## @seealso{gamcdf, gaminv, gamrnd, gamfit, gamlike, gamstat}
## @end deftypefn

function y = gampdf (x, a, b)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("gampdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, A, and B
  if (! isscalar (a) || ! isscalar (b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ("gampdf: X, A, and B must be of common size or scalars.");
    endif
  endif

  ## Check for X, A, and B being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b))
    error ("gampdf: X, A, and B must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (a, "single") || isa (b, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Force NaNs for out of range parameters
  is_nan = ! (a > 0) | ! (b > 0) | isnan (x);
  y(is_nan) = NaN;

  ## Handle all other valid cases
  v = (x >= 0) & (a > 0) & (a <= 1) & (b > 0);
  if (isscalar (a) && isscalar (b))
    y(v) = (x(v) .^ (a - 1)) ...
              .* exp (- x(v) / b) / gamma (a) / (b ^ a);
  else
    y(v) = (x(v) .^ (a(v) - 1)) ...
              .* exp (- x(v) ./ b(v)) ./ gamma (a(v)) ./ (b(v) .^ a(v));
  endif

  v = (x >= 0) & (a > 1) & (b > 0);
  if (isscalar (a) && isscalar (b))
    y(v) = exp (- a * log (b) + (a-1) * log (x(v))
                  - x(v) / b - gammaln (a));
  else
    y(v) = exp (- a(v) .* log (b(v)) + (a(v)-1) .* log (x(v))
                  - x(v) ./ b(v) - gammaln (a(v)));
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
%! legend ({"α = 1, β = 2", "α = 2, β = 2", "α = 3, β = 2", ...
%!          "α = 5, β = 1", "α = 9, β = 0.5", "α = 7.5, β = 1", ...
%!          "α = 0.5, β = 1"}, "location", "northeast")
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
%!error<gampdf: X, A, and B must be of common size or scalars.> ...
%! gampdf (ones (3), ones (2), ones (2))
%!error<gampdf: X, A, and B must be of common size or scalars.> ...
%! gampdf (ones (2), ones (3), ones (2))
%!error<gampdf: X, A, and B must be of common size or scalars.> ...
%! gampdf (ones (2), ones (2), ones (3))
%!error<gampdf: X, A, and B must not be complex.> gampdf (i, 2, 2)
%!error<gampdf: X, A, and B must not be complex.> gampdf (2, i, 2)
%!error<gampdf: X, A, and B must not be complex.> gampdf (2, 2, i)
