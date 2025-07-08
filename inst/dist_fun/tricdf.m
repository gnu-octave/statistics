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
## @deftypefn  {statistics} {@var{p} =} tricdf (@var{x}, @var{a}, @var{b}, @var{c})
## @deftypefnx {statistics} {@var{p} =} tricdf (@var{x}, @var{a}, @var{b}, @var{c}, @qcode{"upper"})
##
## Triangular cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the triangular distribution with lower limit parameter @var{a}, peak
## location (mode) parameter @var{b}, and upper limit parameter @var{c}.  The
## size of @var{p} is the common size of the input arguments.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## @code{@var{p} = tricdf (@var{x}, @var{a}, @var{b}, @var{c}, "upper")}
## computes the upper tail probability of the triangular distribution with
## parameters @var{a}, @var{b}, and @var{c}, at the values in @var{x}.
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
## @seealso{triinv, tripdf, trirnd, tristat}
## @end deftypefn

function p = tricdf (x, a, b, c, uflag)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("tricdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 4)
    if (! strcmpi (uflag, "upper"))
      error ("tricdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of A, B, and C
  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b) || ! isscalar (c))
    [retval, x, a, b, c] = common_size (x, a, b, c);
    if (retval > 0)
      error ("tricdf: X, A, B, and C must be of common size or scalars.");
    endif
  endif

  ## Check for X, BETA, and GAMMA being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b) || iscomplex (c))
    error ("tricdf: X, A, B, and C must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (a, "single") || isa (b, "single") ...
                        || isa (c, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Force NaNs for out of range parameters.
  k = isnan (x) | ! (a < c) | ! (b >= a) | ! (b <= c);
  p(k) = NaN;

  ## Find valid values in parameters and data
  k = (a < c) & (a <= b) & (b <= c);
  k1 = (x <= a) & k;
  k2 = (x > a) & (x <= b) & k;
  k3 = (x > b) & (x < c) & k;
  k4 = (x >= c) & k;

  ## Compute triangular CDF
  if (uflag)
    p(k1) = 1;
    p(k2) = 1 - ((x(k2) - a(k2)) .^ 2) ./ ((c(k2) - a(k2)) .* (b(k2) - a(k2)));
    p(k3) = ((c(k3) - x(k3)) .^ 2) ./ ((c(k3) - a(k3)) .* (c(k3) - b(k3)));
  else
    p(k2) = ((x(k2) - a(k2)) .^ 2) ./ ((c(k2) - a(k2)) .* (b(k2) - a(k2)));
    p(k3) = 1 - ((c(k3) - x(k3)) .^ 2) ./ ((c(k3) - a(k3)) .* (c(k3) - b(k3)));
    p(k4) = 1;
  endif

endfunction

%!demo
%! ## Plot various CDFs from the triangular distribution
%! x = 0.001:0.001:10;
%! p1 = tricdf (x, 3, 4, 6);
%! p2 = tricdf (x, 1, 2, 5);
%! p3 = tricdf (x, 2, 3, 9);
%! p4 = tricdf (x, 2, 5, 9);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c")
%! grid on
%! xlim ([0, 10])
%! legend ({"a = 3, b = 4, c = 6", "a = 1, b = 2, c = 5", ...
%!          "a = 2, b = 3, c = 9", "a = 2, b = 5, c = 9"}, ...
%!         "location", "southeast")
%! title ("Triangular CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-1, 0, 0.1, 0.5, 0.9, 1, 2] + 1;
%! y = [0, 0, 0.02, 0.5, 0.98, 1 1];
%!assert (tricdf (x, ones (1,7), 1.5 * ones (1, 7), 2 * ones (1, 7)), y, eps)
%!assert (tricdf (x, 1 * ones (1, 7), 1.5, 2), y, eps)
%!assert (tricdf (x, 1 * ones (1, 7), 1.5, 2, "upper"), 1 - y, eps)
%!assert (tricdf (x, 1, 1.5, 2 * ones (1, 7)), y, eps)
%!assert (tricdf (x, 1, 1.5 * ones (1, 7), 2), y, eps)
%!assert (tricdf (x, 1, 1.5, 2), y, eps)
%!assert (tricdf (x, [1, 1, NaN, 1, 1, 1, 1], 1.5, 2), ...
%! [y(1:2), NaN, y(4:7)], eps)
%!assert (tricdf (x, 1, 1.5, 2*[1, 1, NaN, 1, 1, 1, 1]), ...
%! [y(1:2), NaN, y(4:7)], eps)
%!assert (tricdf (x, 1, 1.5, 2*[1, 1, NaN, 1, 1, 1, 1]), ...
%! [y(1:2), NaN, y(4:7)], eps)
%!assert (tricdf ([x, NaN], 1, 1.5, 2), [y, NaN], eps)

## Test class of input preserved
%!assert (tricdf (single ([x, NaN]), 1, 1.5, 2), ...
%! single ([y, NaN]), eps("single"))
%!assert (tricdf ([x, NaN], single (1), 1.5, 2), ...
%! single ([y, NaN]), eps("single"))
%!assert (tricdf ([x, NaN], 1, single (1.5), 2), ...
%! single ([y, NaN]), eps("single"))
%!assert (tricdf ([x, NaN], 1, 1.5, single (2)), ...
%! single ([y, NaN]), eps("single"))

## Test input validation
%!error<tricdf: function called with too few input arguments.> tricdf ()
%!error<tricdf: function called with too few input arguments.> tricdf (1)
%!error<tricdf: function called with too few input arguments.> tricdf (1, 2)
%!error<tricdf: function called with too few input arguments.> tricdf (1, 2, 3)
%!error<tricdf: function called with too many inputs> ...
%! tricdf (1, 2, 3, 4, 5, 6)
%!error<tricdf: invalid argument for upper tail.> tricdf (1, 2, 3, 4, "tail")
%!error<tricdf: invalid argument for upper tail.> tricdf (1, 2, 3, 4, 5)
%!error<tricdf: X, A, B, and C must be of common size or scalars.> ...
%! tricdf (ones (3), ones (2), ones(2), ones(2))
%!error<tricdf: X, A, B, and C must be of common size or scalars.> ...
%! tricdf (ones (2), ones (3), ones(2), ones(2))
%!error<tricdf: X, A, B, and C must be of common size or scalars.> ...
%! tricdf (ones (2), ones (2), ones(3), ones(2))
%!error<tricdf: X, A, B, and C must be of common size or scalars.> ...
%! tricdf (ones (2), ones (2), ones(2), ones(3))
%!error<tricdf: X, A, B, and C must not be complex.> tricdf (i, 2, 3, 4)
%!error<tricdf: X, A, B, and C must not be complex.> tricdf (1, i, 3, 4)
%!error<tricdf: X, A, B, and C must not be complex.> tricdf (1, 2, i, 4)
%!error<tricdf: X, A, B, and C must not be complex.> tricdf (1, 2, 3, i)
