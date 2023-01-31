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
## @deftypefn  {statistics} @var{x} = burrinv (@var{p}, @var{a}, @var{c}, @var{k})
##
## Inverse of the Burr type XII cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Burr type XII distribution with parameters @var{a},
## @var{c}, and @var{k}.  The size of @var{x} is the common size of
## @var{p}, @var{a}, @var{c}, and @var{k}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## @seealso{burrcdf, burrpdf, burrrnd}
## @end deftypefn

function x = burrinv (p, a, c, k)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (p) ||! isscalar (a) || ! isscalar (c) || ! isscalar (k))
    [retval, p, a, c, k] = common_size (p, a, c, k);
    if (retval > 0)
      error ("burrinv: P, A, C, AND K must be of common size or scalars.");
    endif
  endif

  if (iscomplex (p) || iscomplex(a) || iscomplex (c) || iscomplex (k))
    error ("burrinv: P, A, C, AND K must not be complex.");
  endif

  if (isa (p, "single") || isa (a, "single") || isa (c, "single") ...
      || isa (k, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  j = isnan (p) | (p < 0) | (p > 1) | ! (a > 0) | ! (c > 0) | ! (k > 0);
  x(j) = NaN;

  j = (p == 1) & (0 < a) & (a < Inf) & (0 < c) & (c < Inf) ...
      & (0 < k) & (k < Inf);
  x(j) = Inf;

  j = (0 < p) & (p < 1) & (0 < a) & (a < Inf) & (0 < c) & (c < Inf) ...
      & (0 < k) & (k < Inf);
  if (isscalar (a) && isscalar(c) && isscalar(k))
    x(j) = ((1 - p(j) / a).^(-1 / k) - 1).^(1 / c) ;
  else
    x(j) = ((1 - p(j) ./ a(j)).^(-1 ./ k(j)) - 1).^(1 ./ c(j)) ;
  endif

endfunction


%!shared p,y
%! p = [-Inf, -1, 0, 1/2, 1, 2, Inf];
%! y = [NaN, NaN, 0, 1 , Inf, NaN, NaN];
%!assert (burrinv (p, ones (1,7), ones (1,7), ones(1,7)), y, eps)
%!assert (burrinv (p, 1, 1, 1), y, eps)
%!assert (burrinv (p, [1, 1, 1, NaN, 1, 1, 1], 1, 1), [y(1:3), NaN, y(5:7)], eps)
%!assert (burrinv (p, 1, [1, 1, 1, NaN, 1, 1, 1], 1), [y(1:3), NaN, y(5:7)], eps)
%!assert (burrinv (p, 1, 1, [1, 1, 1, NaN, 1, 1, 1]), [y(1:3), NaN, y(5:7)], eps)
%!assert (burrinv ([p, NaN], 1, 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (burrinv (single ([p, NaN]), 1, 1, 1), single ([y, NaN]), eps('single'))
%!assert (burrinv ([p, NaN], single (1), 1, 1), single ([y, NaN]), eps('single'))
%!assert (burrinv ([p, NaN], 1, single (1), 1), single ([y, NaN]), eps('single'))
%!assert (burrinv ([p, NaN], 1, 1, single (1)), single ([y, NaN]), eps('single'))

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

