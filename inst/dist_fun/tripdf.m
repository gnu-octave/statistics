## Copyright (C) 1997-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
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
## @deftypefn  {statistics} {@var{y} =} tripdf (@var{x}, @var{a}, @var{b}, @var{c})
##
## Triangular probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the triangular distribution with lower limit parameter @var{a}, peak
## location (mode) parameter @var{b}, and upper limit parameter @var{c}.  The
## size of @var{y} is the common size of the input arguments.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## Note that the order of the parameter input arguments has been changed after
## statistics version 1.6.3 in order to be MATLAB compatible with the parameters
## used in the TriangularDistribution probability distribution object.  More
## specifically, the positions of the parameters @var{b} and @var{c} have been
## swapped.  As a result, the naming conventions no longer coincide with those
## used in Wikipedia, in which @math{b} denotes the upper limit and @math{c}
## denotes the mode or peak parameter.
##
## Further information about the triangular distribution can be found at
## @url{https://en.wikipedia.org/wiki/Triangular_distribution}
##
## @seealso{tricdf, triinv, trirnd, tristat}
## @end deftypefn

function y = tripdf (x, a, b, c)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("tripdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, A, B, and C
  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b) || ! isscalar (c))
    [retval, x, a, b, c] = common_size (x, a, b, c);
    if (retval > 0)
      error ("tripdf: X, A, B, and C must be of common size or scalars.");
    endif
  endif

  ## Check for X, A, B, and C being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b) || iscomplex (c))
    error ("tripdf: X, A, B, and C must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (a, "single") || isa (b, "single") ...
                        || isa (c, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Force NaNs for out of range parameters.
  k = isnan (x) | ! (a < c) | ! (b >= a) | ! (b <= c);
  y(k) = NaN;

  k = (x >= a) & (x <= c) & (a < c) & (a <= b) & (b <= c);
  h = 2 ./ (c - a);

  j = k & (a <= x) & (x < b);
  y(j) = h(j) .* (x(j) - a(j)) ./ (b(j) - a(j));
  j = k & (x == b);
  y(j) = h(j);
  j = k & (b < x) & (x <= c);
  y(j) = h(j) .* (c(j) - x(j)) ./ (c(j) - b(j));

endfunction

%!demo
%! ## Plot various CDFs from the triangular distribution
%! x = 0.001:0.001:10;
%! y1 = tripdf (x, 3, 4, 6);
%! y2 = tripdf (x, 1, 2, 5);
%! y3 = tripdf (x, 2, 3, 9);
%! y4 = tripdf (x, 2, 5, 9);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c")
%! grid on
%! xlim ([0, 10])
%! legend ({"a = 3, b = 4, c = 6", "a = 1, b = 2, c = 5", ...
%!          "a = 2, b = 3, c = 9", "a = 2, b = 5, c = 9"}, ...
%!         "location", "northeast")
%! title ("Triangular CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y, deps
%! x = [-1, 0, 0.1, 0.5, 0.9, 1, 2] + 1;
%! y = [0, 0, 0.4, 2, 0.4, 0, 0];
%! deps = 2*eps;
%!assert (tripdf (x, ones (1,7), 1.5*ones (1,7), 2*ones (1,7)), y, deps)
%!assert (tripdf (x, 1*ones (1,7), 1.5, 2), y, deps)
%!assert (tripdf (x, 1, 1.5, 2*ones (1,7)), y, deps)
%!assert (tripdf (x, 1, 1.5*ones (1,7), 2), y, deps)
%!assert (tripdf (x, 1, 1.5, 2), y, deps)
%!assert (tripdf (x, [1, 1, NaN, 1, 1, 1, 1], 1.5, 2), [y(1:2), NaN, y(4:7)], deps)
%!assert (tripdf (x, 1, 1.5, 2*[1, 1, NaN, 1, 1, 1, 1]), [y(1:2), NaN, y(4:7)], deps)
%!assert (tripdf (x, 1, 1.5*[1, 1, NaN, 1, 1, 1, 1], 2), [y(1:2), NaN, y(4:7)], deps)
%!assert (tripdf ([x, NaN], 1, 1.5, 2), [y, NaN], deps)

## Test class of input preserved
%!assert (tripdf (single ([x, NaN]), 1, 1.5, 2), single ([y, NaN]), eps("single"))
%!assert (tripdf ([x, NaN], single (1), 1.5, 2), single ([y, NaN]), eps("single"))
%!assert (tripdf ([x, NaN], 1, 1.5, single (2)), single ([y, NaN]), eps("single"))
%!assert (tripdf ([x, NaN], 1, single (1.5), 2), single ([y, NaN]), eps("single"))

## Test input validation
%!error<tripdf: function called with too few input arguments.> tripdf ()
%!error<tripdf: function called with too few input arguments.> tripdf (1)
%!error<tripdf: function called with too few input arguments.> tripdf (1, 2)
%!error<tripdf: function called with too few input arguments.> tripdf (1, 2, 3)
%!error<tripdf: function called with too many inputs> ...
%! tripdf (1, 2, 3, 4, 5)
%!error<tripdf: X, A, B, and C must be of common size or scalars.> ...
%! tripdf (ones (3), ones (2), ones(2), ones(2))
%!error<tripdf: X, A, B, and C must be of common size or scalars.> ...
%! tripdf (ones (2), ones (3), ones(2), ones(2))
%!error<tripdf: X, A, B, and C must be of common size or scalars.> ...
%! tripdf (ones (2), ones (2), ones(3), ones(2))
%!error<tripdf: X, A, B, and C must be of common size or scalars.> ...
%! tripdf (ones (2), ones (2), ones(2), ones(3))
%!error<tripdf: X, A, B, and C must not be complex.> tripdf (i, 2, 3, 4)
%!error<tripdf: X, A, B, and C must not be complex.> tripdf (1, i, 3, 4)
%!error<tripdf: X, A, B, and C must not be complex.> tripdf (1, 2, i, 4)
%!error<tripdf: X, A, B, and C must not be complex.> tripdf (1, 2, 3, i)
