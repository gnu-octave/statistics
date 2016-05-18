## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1995-2015 Kurt Hornik
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
## @deftypefn {} {} triinv (@var{x}, @var{a}, @var{b}, @var{c})
## For each element of @var{x}, compute the quantile (the inverse of the CDF)
## at @var{x} of the triangular distribution with parameters
## @var{a}, @var{b}, and @var{c} on the interval [@var{a}, @var{b}].
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: Quantile function of the triangular distribution

function inv = triinv (x, a, b, c)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (a) || ! isscalar (b) || ! isscalar (c))
    [retval, x, a, b, c] = common_size (x, a, b, c);
    if (retval > 0)
      error ("triinv: X, A, B, and C must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (a) || iscomplex (b) || iscomplex (c))
    error ("triinv: X, A, B, and C must not be complex");
  endif

  if (isa (x, "single") || isa (a, "single") || isa (b, "single"))
    inv = NaN (size (x), "single");
  else
    inv = NaN (size (x));
  endif

  k = (x >= 0) & (x <= 1) & (a < b) & (a <= c) & (c <= b);
  inv(k) = 0;
  
  if (isscalar (a) && isscalar (b) && isscalar(c))
    h = 2 / (b-a);
    w = c-a;
    area1 = h * w / 2;
    j = k & (x <= area1);
    inv(j) += (x(j) * (h/2) * w).^0.5 + a;

    w = b-c;
    j = k & (area1 < x) & (x < 1);
    inv(j) += b - ((1-x(j)) * (h/2) * w).^0.5;

    j = k & (x == 1);
    inv(j) = b;
  else
    h = 2 ./ (b-a);
    w = c-a;
    area1 = h .* w / 2;
    j = k & (x <= area1);
    inv(j) += (x(j) .* (h(j)/2) .* w(j)).^0.5 + a(j);

    w = b-c;
    j = k & (area1 < x) & (x < 1);
    inv(j) += b(j) - ((1-x(j)) .* (h(j)/2) .* w(j)).^0.5;

    j = k & (x == 1);
    inv(j) = b(j);
  endif

endfunction


%!shared x,y
%! x = [-1, 0, 0.02, 0.5, 0.98, 1, 2];
%! y = [NaN, 0, 0.1, 0.5, 0.9, 1, NaN] + 1;
%!assert (triinv (x, ones (1,7), 2*ones (1,7), 1.5*ones (1,7)), y, eps)
%!assert (triinv (x, 1*ones (1,7), 2, 1.5), y, eps)
%!assert (triinv (x, 1, 2*ones (1,7), 1.5), y, eps)
%!assert (triinv (x, 1, 2, 1.5*ones (1,7)), y, eps)
%!assert (triinv (x, 1, 2, 1.5), y, eps)
%!assert (triinv (x, [1, 1, NaN, 1, 1, 1, 1], 2, 1.5), [y(1:2), NaN, y(4:7)], eps)
%!assert (triinv (x, 1, 2*[1, 1, NaN, 1, 1, 1, 1], 1.5), [y(1:2), NaN, y(4:7)], eps)
%!assert (triinv (x, 1, 2, 1.5*[1, 1, NaN, 1, 1, 1, 1]), [y(1:2), NaN, y(4:7)], eps)
%!assert (triinv ([x, NaN], 1, 2, 1.5), [y, NaN], eps)

## Test class of input preserved
%!assert (triinv (single ([x, NaN]), 1, 2, 1.5), single ([y, NaN]), eps('single'))
%!assert (triinv ([x, NaN], single (1), 2, 1.5), single ([y, NaN]), eps('single'))
%!assert (triinv ([x, NaN], 1, single (2), 1.5), single ([y, NaN]), eps('single'))
%!assert (triinv ([x, NaN], 1, 2, single (1.5)), single ([y, NaN]), eps('single'))

## Test input validation
%!error triinv ()
%!error triinv (1)
%!error triinv (1,2)
%!error triinv (1,2,3)
%!error triinv (1,2,3,4,5)
%!error triinv (1, ones (3), ones (2), ones (2))
%!error triinv (1, ones (2), ones (3), ones (2))
%!error triinv (1, ones (2), ones (2), ones (3))
%!error triinv (i, 2, 2, 2)
%!error triinv (2, i, 2, 2)
%!error triinv (2, 2, i, 2)
%!error triinv (2, 2, 2, i)

