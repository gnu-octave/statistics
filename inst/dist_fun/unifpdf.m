## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{y} =} unifpdf (@var{x})
## @deftypefnx {statistics} {@var{y} =} unifpdf (@var{x}, @var{a})
## @deftypefnx {statistics} {@var{y} =} unifpdf (@var{x}, @var{a}, @var{b})
##
## Continuous uniform probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the continuous uniform distribution on the interval @qcode{[@var{a},
## @var{b}]}.  The size of @var{y} is the common size of @var{x}, @var{a}, and
## @var{b}.  A scalar input functions as a constant matrix of the same size as
## the other inputs.
##
## Further information about the continuous uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Continuous_uniform_distribution}
##
## @seealso{unifcdf, unifinv, unifrnd, unifit, unifstat}
## @end deftypefn

function y = unifpdf (x, a, b)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("unifpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, A, and B
  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ("unifpdf: X, A, and B must be of common size or scalars.");
    endif
  endif

  ## Check for X, A, and B being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b))
    error ("unifpdf: X, A, and B must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (a, "single") || isa (b, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Calculate continuous uniform PDF for valid parameter and data range
  k = isnan (x) | ! (a < b);
  y(k) = NaN;

  k = (x >= a) & (x <= b) & (a < b);
  y(k) = 1 ./ (b(k) - a(k));

endfunction

%!demo
%! ## Plot various PDFs from the continuous uniform distribution
%! x = 0:0.001:10;
%! y1 = unifpdf (x, 2, 5);
%! y2 = unifpdf (x, 3, 9);
%! plot (x, y1, "-b", x, y2, "-g")
%! grid on
%! xlim ([0, 10])
%! ylim ([0, 0.4])
%! legend ({"a = 2, b = 5", "a = 3, b = 9"}, "location", "northeast")
%! title ("Continuous uniform PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 0.5 1 2] + 1;
%! y = [0 1 1 1 0];
%!assert (unifpdf (x, ones (1,5), 2*ones (1,5)), y)
%!assert (unifpdf (x, 1, 2*ones (1,5)), y)
%!assert (unifpdf (x, ones (1,5), 2), y)
%!assert (unifpdf (x, [2 NaN 1 1 1], 2), [NaN NaN y(3:5)])
%!assert (unifpdf (x, 1, 2*[0 NaN 1 1 1]), [NaN NaN y(3:5)])
%!assert (unifpdf ([x, NaN], 1, 2), [y, NaN])
%!assert (unifpdf (x, 0, 1), [1 1 0 0 0])

## Test class of input preserved
%!assert (unifpdf (single ([x, NaN]), 1, 2), single ([y, NaN]))
%!assert (unifpdf (single ([x, NaN]), single (1), 2), single ([y, NaN]))
%!assert (unifpdf ([x, NaN], 1, single (2)), single ([y, NaN]))

## Test input validation
%!error<unifpdf: function called with too few input arguments.> unifpdf ()
%!error<unifpdf: function called with too few input arguments.> unifpdf (1)
%!error<unifpdf: function called with too few input arguments.> unifpdf (1, 2)
%!error<unifpdf: X, A, and B must be of common size or scalars.> ...
%! unifpdf (ones (3), ones (2), ones (2))
%!error<unifpdf: X, A, and B must be of common size or scalars.> ...
%! unifpdf (ones (2), ones (3), ones (2))
%!error<unifpdf: X, A, and B must be of common size or scalars.> ...
%! unifpdf (ones (2), ones (2), ones (3))
%!error<unifpdf: X, A, and B must not be complex.> unifpdf (i, 2, 2)
%!error<unifpdf: X, A, and B must not be complex.> unifpdf (2, i, 2)
%!error<unifpdf: X, A, and B must not be complex.> unifpdf (2, 2, i)
