## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1997-2015 Kurt Hornik
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {} tricdf (@var{x}, @var{a}, @var{b}, @var{c})
## Compute the cumulative distribution function (CDF) at @var{x} of the
## triangular distribution with parameters @var{a}, @var{b}, and @var{c}
## on the interval [@var{a}, @var{b}].
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: CDF of the triangle distribution

function cdf = tricdf (x, a, b, c)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (a) || ! isscalar (b) || ! isscalar (c))
    [retval, x, a, b, c] = common_size (x, a, b, c);
    if (retval > 0)
      error ("tricdf: X, A, B, and C must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (a) || iscomplex (b) || iscomplex (c))
    error ("tricdf: X, A, B, and C must not be complex");
  endif

  if (isa (x, "single") || isa (a, "single")
      || isa (b, "single") || isa (c, "single"))
    cdf = zeros (size (x), "single");
  else
    cdf = zeros (size (x));
  endif

  k = isnan (x) | !(a < b) | !(c >= a) | !(c <= b) ;
  cdf(k) = NaN;

  k = (x > a) & (a < b) & (a <= c) & (c <= b);
  if (isscalar (a) && isscalar (b) && isscalar (c))
    h = 2 / (b-a);

    k_temp = k & (c <= x);
    full_area = (c-a) * h / 2;
    cdf(k_temp) += full_area;

    k_temp = k & (a < x) & (x < c);
    area = (x(k_temp) - a).^2 * h;
    cdf(k_temp) += area;

    k_temp = k & (b <= x);
    full_area = (b-c) * h / 2;
    cdf(k_temp) += full_area;

    k_temp = k & (c < x) & (x < b);
    area = (b-x(k_temp)).^2 * h;
    cdf(k_temp) += full_area - area;
  else
    h = 2 ./ (b-a);

    k_temp = k & (c <= x);
    full_area = (c(k_temp)-a(k_temp)) .* h(k_temp) / 2;
    cdf(k_temp) += full_area;

    k_temp = k & (a <= x) & (x < c);
    area = (x(k_temp) - a(k_temp)).^2 .* h(k_temp);
    cdf(k_temp) += area;

    k_temp = k & (b <= x);
    full_area = (b(k_temp)-c(k_temp)) .* h(k_temp) / 2;
    cdf(k_temp) += full_area;

    k_temp = k & (c <= x) & (x < b);
    area = (b(k_temp)-x(k_temp)).^2 .* h(k_temp);
    cdf(k_temp) += full_area - area;
  endif

endfunction


%!shared x,y
%! x = [-1, 0, 0.1, 0.5, 0.9, 1, 2] + 1;
%! y = [0, 0, 0.02, 0.5, 0.98, 1 1];
%!assert (tricdf (x, ones (1,7), 2*ones (1,7), 1.5*ones (1,7)), y, eps)
%!assert (tricdf (x, 1*ones (1,7), 2, 1.5), y, eps)
%!assert (tricdf (x, 1, 2*ones (1,7), 1.5), y, eps)
%!assert (tricdf (x, 1, 2, 1.5*ones (1,7)), y, eps)
%!assert (tricdf (x, 1, 2, 1.5), y, eps)
%!assert (tricdf (x, [1, 1, NaN, 1, 1, 1, 1], 2, 1.5), [y(1:2), NaN, y(4:7)], eps)
%!assert (tricdf (x, 1, 2*[1, 1, NaN, 1, 1, 1, 1], 1.5), [y(1:2), NaN, y(4:7)], eps)
%!assert (tricdf (x, 1, 2, 1.5*[1, 1, NaN, 1, 1, 1, 1]), [y(1:2), NaN, y(4:7)], eps)
%!assert (tricdf ([x, NaN], 1, 2, 1.5), [y, NaN], eps)

## Test class of input preserved
%!assert (tricdf (single ([x, NaN]), 1, 2, 1.5), single ([y, NaN]), eps('single'))
%!assert (tricdf ([x, NaN], single (1), 2, 1.5), single ([y, NaN]), eps('single'))
%!assert (tricdf ([x, NaN], 1, single (2), 1.5), single ([y, NaN]), eps('single'))
%!assert (tricdf ([x, NaN], 1, 2, single (1.5)), single ([y, NaN]), eps('single'))

## Test input validation
%!error tricdf ()
%!error tricdf (1)
%!error tricdf (1,2)
%!error tricdf (1,2,3)
%!error tricdf (1,2,3,4,5)
%!error tricdf (1, ones (3), ones (2), ones (2))
%!error tricdf (1, ones (2), ones (3), ones (2))
%!error tricdf (1, ones (2), ones (2), ones (3))
%!error tricdf (i, 2, 2, 2)
%!error tricdf (2, i, 2, 2)
%!error tricdf (2, 2, i, 2)
%!error tricdf (2, 2, 2, i)

