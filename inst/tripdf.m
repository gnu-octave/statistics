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
## @deftypefn {} {} tripdf (@var{x}, @var{a}, @var{b}, @var{c})
## Compute the probability density function (PDF) at @var{x} of the triangular
## distribution with parameters @var{a}, @var{b}, and @var{c} on the interval
## [@var{a}, @var{b}].
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: PDF of the triangular distribution

function pdf = tripdf (x, a, b, c)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (a) || ! isscalar (b) || ! isscalar (c))
    [retval, x, a, b, c] = common_size (x, a, b, c);
    if (retval > 0)
      error ("tripdf: X, A, B, and C must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (a) || iscomplex (b) || iscomplex (c))
    error ("tripdf: X, A, B, and C must not be complex");
  endif

  if (isa (x, "single") || isa (a, "single") ...
      || isa (b, "single") || isa (c, "single"))
    pdf = zeros (size (x), "single");
  else
    pdf = zeros (size (x));
  endif

  k = isnan (x) | !(a < b) | !(c >= a) | !(c <= b) ;
  pdf(k) = NaN;

  k = (x >= a) & (x <= b) & (a < b) & (a <= c) & (c <= b);
  h = 2 ./ (b-a);
  if (isscalar (a) && isscalar (b) && isscalar (c))
    j = k & (a <= x) & (x < c);
    pdf(j) = h * (x(j)-a) / (c-a);
    j = k & (x == c);
    pdf(j) = h;
    j = k & (c < x) & (x <= b);
    pdf(j) = h * (b-x(j)) / (b-c);
  else
    j = k & (a <= x) & (x < c);
    pdf(j) = h(j) .* (x(j)-a(j)) ./ (c(j)-a(j));
    j = k & (x == c);
    pdf(j) = h(j);
    j = k & (c < x) & (x <= b);
    pdf(j) = h(j) .* (b(j)-x(j)) ./ (b(j)-c(j));
  endif

endfunction


%!shared x,y,deps
%! x = [-1, 0, 0.1, 0.5, 0.9, 1, 2] + 1;
%! y = [0, 0, 0.4, 2, 0.4, 0, 0];
%! deps = 2*eps;
%!assert (tripdf (x, ones (1,7), 2*ones (1,7), 1.5*ones (1,7)), y, deps)
%!assert (tripdf (x, 1*ones (1,7), 2, 1.5), y, deps)
%!assert (tripdf (x, 1, 2*ones (1,7), 1.5), y, deps)
%!assert (tripdf (x, 1, 2, 1.5*ones (1,7)), y, deps)
%!assert (tripdf (x, 1, 2, 1.5), y, deps)
%!assert (tripdf (x, [1, 1, NaN, 1, 1, 1, 1], 2, 1.5), [y(1:2), NaN, y(4:7)], deps)
%!assert (tripdf (x, 1, 2*[1, 1, NaN, 1, 1, 1, 1], 1.5), [y(1:2), NaN, y(4:7)], deps)
%!assert (tripdf (x, 1, 2, 1.5*[1, 1, NaN, 1, 1, 1, 1]), [y(1:2), NaN, y(4:7)], deps)
%!assert (tripdf ([x, NaN], 1, 2, 1.5), [y, NaN], deps)

## Test class of input preserved
%!assert (tripdf (single ([x, NaN]), 1, 2, 1.5), single ([y, NaN]), eps('single'))
%!assert (tripdf ([x, NaN], single (1), 2, 1.5), single ([y, NaN]), eps('single'))
%!assert (tripdf ([x, NaN], 1, single (2), 1.5), single ([y, NaN]), eps('single'))
%!assert (tripdf ([x, NaN], 1, 2, single (1.5)), single ([y, NaN]), eps('single'))

## Test input validation
%!error tripdf ()
%!error tripdf (1)
%!error tripdf (1,2)
%!error tripdf (1,2,3)
%!error tripdf (1,2,3,4,5)
%!error tripdf (1, ones (3), ones (2), ones (2))
%!error tripdf (1, ones (2), ones (3), ones (2))
%!error tripdf (1, ones (2), ones (2), ones (3))
%!error tripdf (i, 2, 2, 2)
%!error tripdf (2, i, 2, 2)
%!error tripdf (2, 2, i, 2)
%!error tripdf (2, 2, 2, i)

