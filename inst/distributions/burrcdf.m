## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} @var{p} = burrcdf (@var{x}, @var{c}, @var{k})
##
## Burr type XII cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the Burr type XII distribution with parameters @var{a},
## @var{c}, and @var{k}.  The size of @var{x} is the common size of @var{p},
## @var{a}, @var{c}, and @var{k}.  A scalar input functions as a constant matrix
## of the same size as the other inputs.
##
## @seealso{burrinv, burrpdf, burrrnd}
## @end deftypefn

function p = burrcdf (x, a, c, k)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (x) || ! isscalar (a) || ! isscalar (c) || ! isscalar (k))
    [retval, x, a, c, k] = common_size (x, a, c, k);
    if (retval > 0)
      error ("burrcdf: X, A, C, AND K must be of common size or scalars.");
    endif
  endif

  if (iscomplex (x) || iscomplex(a) || iscomplex (c) || iscomplex (k))
    error ("burrcdf: X, A, C, AND K must not be complex.");
  endif

  if (isa (x, "single") || isa (a, "single") || isa (c, "single") ...
      || isa (k, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  j = isnan (x) | ! (a > 0) | ! (c > 0) | ! (k > 0);
  p(j) = NaN;

  j = (x > 0) & (0 < a) & (a < Inf) & (0 < c) & (c < Inf) ...
    & (0 < k) & (k < Inf);
  if (isscalar (a) && isscalar(c) && isscalar(k))
    p(j) = 1 - (1 + (x(j) / a).^c).^(-k);
  else
    p(j) = 1 - (1 + (x(j) ./ a(j)).^c(j)).^(-k(j));
  endif

endfunction


%!shared x,y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 1/2, 2/3, 1];
%!assert (burrcdf (x, ones(1,5), ones (1,5), ones (1,5)), y, eps)
%!assert (burrcdf (x, 1, 1, 1), y, eps)
%!assert (burrcdf (x, [1, 1, NaN, 1, 1], 1, 1), [y(1:2), NaN, y(4:5)], eps)
%!assert (burrcdf (x, 1, [1, 1, NaN, 1, 1], 1), [y(1:2), NaN, y(4:5)], eps)
%!assert (burrcdf (x, 1, 1, [1, 1, NaN, 1, 1]), [y(1:2), NaN, y(4:5)], eps)
%!assert (burrcdf ([x, NaN], 1, 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (burrcdf (single ([x, NaN]), 1, 1, 1), single ([y, NaN]), eps('single'))
%!assert (burrcdf ([x, NaN], single (1), 1, 1), single ([y, NaN]), eps('single'))
%!assert (burrcdf ([x, NaN], 1, single (1), 1), single ([y, NaN]), eps('single'))
%!assert (burrcdf ([x, NaN], 1, 1, single (1)), single ([y, NaN]), eps('single'))

## Test input validation
%!error burrcdf ()
%!error burrcdf (1)
%!error burrcdf (1,2)
%!error burrcdf (1,2,3)
%!error burrcdf (1,2,3,4,5)
%!error burrcdf (ones (3), ones (2), ones(2), ones(2))
%!error burrcdf (ones (2), ones (3), ones(2), ones(2))
%!error burrcdf (ones (2), ones (2), ones(3), ones(2))
%!error burrcdf (ones (2), ones (2), ones(2), ones(3))
%!error burrcdf (i, 2, 2, 2)
%!error burrcdf (2, i, 2, 2)
%!error burrcdf (2, 2, i, 2)
%!error burrcdf (2, 2, 2, i)

