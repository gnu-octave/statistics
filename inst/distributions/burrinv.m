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
## @deftypefn {} {} burrinv (@var{x}, @var{alpha}, @var{c}, @var{k})
## For each element of @var{x}, compute the quantile (the inverse of the CDF)
## at @var{x} of the Burr distribution with scale parameter @var{alpha} and
## shape parameters @var{c} and @var{k}.
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: Quantile function of the Burr distribution

function inv = burrinv (x, alpha, c, k)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (alpha) || ! isscalar (c) || ! isscalar (k) )
    [retval, x, alpha, c, k] = common_size (x, alpha, c, k);
    if (retval > 0)
      error ("burrinv: X, ALPHA, C AND K must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex(alpha) || iscomplex (c) || iscomplex (k))
    error ("burrinv: X, ALPHA, C AND K must not be complex");
  endif

  if (isa (x, "single") || isa (alpha, "single") || isa (c, "single") ...
      || isa (k, "single"))
    inv = zeros (size (x), "single");
  else
    inv = zeros (size (x));
  endif

  j = isnan (x) | (x < 0) | (x > 1) | ! (alpha > 0) | ! (c > 0) | ! (k > 0);
  inv(j) = NaN;

  j = (x == 1) & (0 < alpha) & (alpha < Inf) & (0 < c) & (c < Inf) ...
      & (0 < k) & (k < Inf);
  inv(j) = Inf;

  j = (0 < x) & (x < 1) & (0 < alpha) & (alpha < Inf) & (0 < c) & (c < Inf) ...
      & (0 < k) & (k < Inf);
  if (isscalar (alpha) && isscalar(c) && isscalar(k))
    inv(j) = ((1 - x(j) / alpha).^(-1 / k) - 1).^(1 / c) ;
  else
    inv(j) = ((1 - x(j) ./ alpha(j)).^(-1 ./ k(j)) - 1).^(1 ./ c(j)) ;
  endif

endfunction


%!shared x,y
%! x = [-Inf, -1, 0, 1/2, 1, 2, Inf];
%! y = [NaN, NaN, 0, 1 , Inf, NaN, NaN];
%!assert (burrinv (x, ones (1,7), ones (1,7), ones(1,7)), y, eps)
%!assert (burrinv (x, 1, 1, 1), y, eps)
%!assert (burrinv (x, [1, 1, 1, NaN, 1, 1, 1], 1, 1), [y(1:3), NaN, y(5:7)], eps)
%!assert (burrinv (x, 1, [1, 1, 1, NaN, 1, 1, 1], 1), [y(1:3), NaN, y(5:7)], eps)
%!assert (burrinv (x, 1, 1, [1, 1, 1, NaN, 1, 1, 1]), [y(1:3), NaN, y(5:7)], eps)
%!assert (burrinv ([x, NaN], 1, 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (burrinv (single ([x, NaN]), 1, 1, 1), single ([y, NaN]), eps('single'))
%!assert (burrinv ([x, NaN], single (1), 1, 1), single ([y, NaN]), eps('single'))
%!assert (burrinv ([x, NaN], 1, single (1), 1), single ([y, NaN]), eps('single'))
%!assert (burrinv ([x, NaN], 1, 1, single (1)), single ([y, NaN]), eps('single'))

## Test input validation
%!error burrinv ()
%!error burrinv (1)
%!error burrinv (1,2)
%!error burrinv (1,2,3)
%!error burrinv (1,2,3,4,5)
%!error burrinv (ones (3), ones (2), ones(2), ones(2))
%!error burrinv (ones (2), ones (3), ones(2), ones(2))
%!error burrinv (ones (2), ones (2), ones(3), ones(2))
%!error burrinv (ones (2), ones (2), ones(2), ones(3))
%!error burrinv (i, 2, 2, 2)
%!error burrinv (2, i, 2, 2)
%!error burrinv (2, 2, i, 2)
%!error burrinv (2, 2, 2, i)

