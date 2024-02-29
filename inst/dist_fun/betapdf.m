## Copyright (C) 2010 Christos Dimitrakakis
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
## @deftypefn  {statistics} {@var{y} =} betapdf (@var{x}, @var{a}, @var{b})
##
## Beta probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Beta distribution with shape parameters @var{a} and @var{b}.  The size
## of @var{y} is the common size of @var{x}, @var{a}, and @var{b}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## Further information about the Beta distribution can be found at
## @url{https://en.wikipedia.org/wiki/Beta_distribution}
##
## @seealso{betacdf, betainv, betarnd, betafit, betalike, betastat}
## @end deftypefn

function y = betapdf (x, a, b)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("betapdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, A, and B
  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ("betapdf: X, A, and B must be of common size or scalars.");
    endif
  endif

  ## Check for X, A, and B being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b))
    error ("betapdf: X, A, and B must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (a, "single") || isa (b, "single"));
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  k = !(a > 0) | !(b > 0) | isnan (x);
  y(k) = NaN;

  k = (x > 0) & (x < 1) & (a > 0) & (b > 0) & ((a != 1) | (b != 1));
  if (isscalar (a) && isscalar (b))
    y(k) = exp ((a - 1) * log (x(k))
                  + (b - 1) * log (1 - x(k))
                  + gammaln (a + b) - gammaln (a) - gammaln (b));
  else
    y(k) = exp ((a(k) - 1) .* log (x(k))
                  + (b(k) - 1) .* log (1 - x(k))
                  + gammaln (a(k) + b(k)) - gammaln (a(k)) - gammaln (b(k)));
  endif

  ## Most important special cases when the density is finite.
  k = (x == 0) & (a == 1) & (b > 0) & (b != 1);
  if (isscalar (a) && isscalar (b))
    y(k) = exp (gammaln (a + b) - gammaln (a) - gammaln (b));
  else
    y(k) = exp (gammaln (a(k) + b(k)) - gammaln (a(k)) - gammaln (b(k)));
  endif

  k = (x == 1) & (b == 1) & (a > 0) & (a != 1);
  if (isscalar (a) && isscalar (b))
    y(k) = exp (gammaln (a + b) - gammaln (a) - gammaln (b));
  else
    y(k) = exp (gammaln (a(k) + b(k)) - gammaln (a(k)) - gammaln (b(k)));
  endif

  k = (x >= 0) & (x <= 1) & (a == 1) & (b == 1);
  y(k) = 1;

  ## Other special case when the density at the boundary is infinite.
  k = (x == 0) & (a < 1);
  y(k) = Inf;

  k = (x == 1) & (b < 1);
  y(k) = Inf;

endfunction

%!demo
%! ## Plot various PDFs from the Beta distribution
%! x = 0.001:0.001:0.999;
%! y1 = betapdf (x, 0.5, 0.5);
%! y2 = betapdf (x, 5, 1);
%! y3 = betapdf (x, 1, 3);
%! y4 = betapdf (x, 2, 2);
%! y5 = betapdf (x, 2, 5);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c", x, y5, "-m")
%! grid on
%! ylim ([0, 2.5])
%! legend ({"α = β = 0.5", "α = 5, β = 1", "α = 1, β = 3", ...
%!          "α = 2, β = 2", "α = 2, β = 5"}, "location", "north")
%! title ("Beta PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 0.5 1 2];
%! y = [0 2 1 0 0];
%!assert (betapdf (x, ones (1, 5), 2 * ones (1, 5)), y)
%!assert (betapdf (x, 1, 2 * ones (1, 5)), y)
%!assert (betapdf (x, ones (1, 5), 2), y)
%!assert (betapdf (x, [0 NaN 1 1 1], 2), [NaN NaN y(3:5)])
%!assert (betapdf (x, 1, 2 * [0 NaN 1 1 1]), [NaN NaN y(3:5)])
%!assert (betapdf ([x, NaN], 1, 2), [y, NaN])

## Test class of input preserved
%!assert (betapdf (single ([x, NaN]), 1, 2), single ([y, NaN]))
%!assert (betapdf ([x, NaN], single (1), 2), single ([y, NaN]))
%!assert (betapdf ([x, NaN], 1, single (2)), single ([y, NaN]))

## Beta (1/2,1/2) == arcsine distribution
%!test
%! x = rand (10,1);
%! y = 1 ./ (pi * sqrt (x .* (1 - x)));
%! assert (betapdf (x, 1/2, 1/2), y, 1e-12);

## Test large input values to betapdf
%!assert (betapdf (0.5, 1000, 1000), 35.678, 1e-3)

## Test input validation
%!error<betapdf: function called with too few input arguments.> betapdf ()
%!error<betapdf: function called with too few input arguments.> betapdf (1)
%!error<betapdf: function called with too few input arguments.> betapdf (1,2)
%!error<betapdf: function called with too many inputs> betapdf (1,2,3,4)
%!error<betapdf: X, A, and B must be of common size or scalars.> ...
%! betapdf (ones (3), ones (2), ones (2))
%!error<betapdf: X, A, and B must be of common size or scalars.> ...
%! betapdf (ones (2), ones (3), ones (2))
%!error<betapdf: X, A, and B must be of common size or scalars.> ...
%! betapdf (ones (2), ones (2), ones (3))
%!error<betapdf: X, A, and B must not be complex.> betapdf (i, 2, 2)
%!error<betapdf: X, A, and B must not be complex.> betapdf (2, i, 2)
%!error<betapdf: X, A, and B must not be complex.> betapdf (2, 2, i)
