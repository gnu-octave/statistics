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
## @deftypefn  {statistics} {@var{p} =} unifcdf (@var{x}, @var{a}, @var{b})
## @deftypefnx {statistics} {@var{p} =} unifcdf (@var{x}, @var{a}, @var{b}, @qcode{"upper"})
##
## Continuous uniform cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the continuous uniform distribution with parameters @var{a} and
## @var{b}, which define the lower and upper bounds of the interval
## @qcode{[@var{a}, @var{b}]}.  The size of @var{p} is the common size of
## @var{x}, @var{a}, and @var{b}.  A scalar input functions as a constant matrix
## of the same size as the other inputs.
##
## @code{[@dots{}] = unifcdf (@var{x}, @var{a}, @var{b}, "upper")} computes the
## upper tail probability of the continuous uniform distribution with parameters
## @var{a}, and @var{b}, at the values in @var{x}.
##
## Further information about the continuous uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Continuous_uniform_distribution}
##
## @seealso{unifinv, unifpdf, unifrnd, unifit, unifstat}
## @end deftypefn

function p = unifcdf (x, a, b, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("unifcdf: function called with too few input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 3 && strcmpi (uflag, "upper"))
    uflag = true;
  elseif (nargin > 3  && ! strcmpi (uflag, "upper"))
    error ("unifcdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Check for common size of X, A, and B
  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ("unifcdf: X, A, and B must be of common size or scalars.");
    endif
  endif

  ## Check for X, A, and B being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b))
    error ("unifcdf: X, A, and B must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (a, "single") || isa (b, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Calculate continuous uniform CDF for valid parameter and data range
  k = find(x > a & x < b & a < b);
  if (uflag)
    p(x <= a & a < b) = 1;
    p(x >= b & a < b) = 0;
    if any(k)
      p(k) = (b(k)- x(k)) ./ (b(k) - a(k));
    endif
  else
    p(x <= a & a < b) = 0;
    p(x >= b & a < b) = 1;
    if any(k)
      p(k) = (x(k) - a(k)) ./ (b(k) - a(k));
    endif
  endif

  ## Continue argument check
  p(a >= b) = NaN;
  p(isnan(x) | isnan(a) | isnan(b)) = NaN;

endfunction

%!demo
%! ## Plot various CDFs from the continuous uniform distribution
%! x = 0:0.1:10;
%! p1 = unifcdf (x, 2, 5);
%! p2 = unifcdf (x, 3, 9);
%! plot (x, p1, "-b", x, p2, "-g")
%! grid on
%! xlim ([0, 10])
%! ylim ([0, 1])
%! legend ({"a = 2, b = 5", "a = 3, b = 9"}, "location", "southeast")
%! title ("Continuous uniform CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-1 0 0.5 1 2] + 1;
%! y = [0 0 0.5 1 1];
%!assert (unifcdf (x, ones (1,5), 2*ones (1,5)), y)
%!assert (unifcdf (x, ones (1,5), 2*ones (1,5), "upper"), 1 - y)
%!assert (unifcdf (x, 1, 2*ones (1,5)), y)
%!assert (unifcdf (x, 1, 2*ones (1,5), "upper"), 1 - y)
%!assert (unifcdf (x, ones (1,5), 2), y)
%!assert (unifcdf (x, ones (1,5), 2, "upper"), 1 - y)
%!assert (unifcdf (x, [2 1 NaN 1 1], 2), [NaN 0 NaN 1 1])
%!assert (unifcdf (x, [2 1 NaN 1 1], 2, "upper"), 1 - [NaN 0 NaN 1 1])
%!assert (unifcdf (x, 1, 2*[0 1 NaN 1 1]), [NaN 0 NaN 1 1])
%!assert (unifcdf (x, 1, 2*[0 1 NaN 1 1], "upper"), 1 - [NaN 0 NaN 1 1])
%!assert (unifcdf ([x(1:2) NaN x(4:5)], 1, 2), [y(1:2) NaN y(4:5)])
%!assert (unifcdf ([x(1:2) NaN x(4:5)], 1, 2, "upper"), 1 - [y(1:2) NaN y(4:5)])

## Test class of input preserved
%!assert (unifcdf ([x, NaN], 1, 2), [y, NaN])
%!assert (unifcdf (single ([x, NaN]), 1, 2), single ([y, NaN]))
%!assert (unifcdf ([x, NaN], single (1), 2), single ([y, NaN]))
%!assert (unifcdf ([x, NaN], 1, single (2)), single ([y, NaN]))

## Test input validation
%!error<unifcdf: function called with too few input arguments.> unifcdf ()
%!error<unifcdf: function called with too few input arguments.> unifcdf (1)
%!error<unifcdf: function called with too few input arguments.> unifcdf (1, 2)
%!error<unifcdf: invalid argument for upper tail.> unifcdf (1, 2, 3, 4)
%!error<unifcdf: invalid argument for upper tail.> unifcdf (1, 2, 3, "tail")
%!error<unifcdf: X, A, and B must be of common size or scalars.> ...
%! unifcdf (ones (3), ones (2), ones (2))
%!error<unifcdf: X, A, and B must be of common size or scalars.> ...
%! unifcdf (ones (2), ones (3), ones (2))
%!error<unifcdf: X, A, and B must be of common size or scalars.> ...
%! unifcdf (ones (2), ones (2), ones (3))
%!error<unifcdf: X, A, and B must not be complex.> unifcdf (i, 2, 2)
%!error<unifcdf: X, A, and B must not be complex.> unifcdf (2, i, 2)
%!error<unifcdf: X, A, and B must not be complex.> unifcdf (2, 2, i)
