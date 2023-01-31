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
## @deftypefn  {statistics} @var{y} = unifpdf (@var{x})
## @deftypefnx {statistics} @var{y} = unifpdf (@var{x}, @var{a})
## @deftypefnx {statistics} @var{y} = unifpdf (@var{x}, @var{a}, @var{b})
##
## Uniform probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the uniform distribution on the interval [@var{a}, @var{b}].
## The size of @var{y} is the common size of the input arguments.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## Default values are @var{a} = 0, @var{b} = 1.
##
## @seealso{unifcdf, unifinv, unifrnd, unifstat}
## @end deftypefn

function y = unifpdf (x, a = 0, b = 1)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ("unifpdf: X, A, and B must be of common size or scalars.");
    endif
  endif

  if (iscomplex (x) || iscomplex (a) || iscomplex (b))
    error ("unifpdf: X, A, and B must not be complex.");
  endif

  if (isa (x, "single") || isa (a, "single") || isa (b, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  k = isnan (x) | !(a < b);
  y(k) = NaN;

  k = (x >= a) & (x <= b) & (a < b);
  y(k) = 1 ./ (b(k) - a(k));

endfunction


%!shared x,y
%! x = [-1 0 0.5 1 2] + 1;
%! y = [0 1 1 1 0];
%!assert (unifpdf (x, ones (1,5), 2*ones (1,5)), y)
%!assert (unifpdf (x, 1, 2*ones (1,5)), y)
%!assert (unifpdf (x, ones (1,5), 2), y)
%!assert (unifpdf (x, [2 NaN 1 1 1], 2), [NaN NaN y(3:5)])
%!assert (unifpdf (x, 1, 2*[0 NaN 1 1 1]), [NaN NaN y(3:5)])
%!assert (unifpdf ([x, NaN], 1, 2), [y, NaN])
%!assert (unifpdf (x), [1 1 0 0 0])
%!assert (unifpdf (x, 0), [1 1 0 0 0])
%!assert (unifpdf (x, 0, 1), [1 1 0 0 0])

## Test class of input preserved
%!assert (unifpdf (single ([x, NaN]), 1, 2), single ([y, NaN]))
%!assert (unifpdf (single ([x, NaN]), single (1), 2), single ([y, NaN]))
%!assert (unifpdf ([x, NaN], 1, single (2)), single ([y, NaN]))

## Test input validation
%!error unifpdf ()
%!error unifpdf (1,2,3,4)
%!error unifpdf (ones (3), ones (2), ones (2))
%!error unifpdf (ones (2), ones (3), ones (2))
%!error unifpdf (ones (2), ones (2), ones (3))
%!error unifpdf (i, 2, 2)
%!error unifpdf (2, i, 2)
%!error unifpdf (2, 2, i)
