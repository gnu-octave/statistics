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
## @deftypefn {} {} nakainv (@var{x}, @var{m}, @var{w})
## For each element of @var{x}, compute the quantile (the inverse of the CDF)
## at @var{x} of the Nakagami distribution with shape parameter @var{m} and
## scale parameter @var{w}.
##
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: Quantile function of the Nakagami distribution

function inv = nakainv (x, m, w)

  if (nargin != 3)
    print_usage ();
  endif

  if (! isscalar (m) || ! isscalar (w))
    [retval, x, m, w] = common_size (x, m, w);
    if (retval > 0)
      error ("nakainv: X, M and W must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (m) || iscomplex (w))
    error ("nakainv: X, M, and W must not be complex");
  endif

  if (isa (x, "single") || isa (m, "single") || isa (w, "single"))
    inv = zeros (size (x), "single");
  else
    inv = zeros (size (x));
  endif

  k = isnan (x) | ! (0 <= x) | ! (x <= 1) | ! (-Inf < m) | ! (m < Inf) ...
    | ! (0 < w) | ! (w < Inf);
  inv(k) = NaN;

  k = (x == 1) & (-Inf < m) & (m < Inf) & (0 < w) & (w < Inf);
  inv(k) = Inf;

  k = (0 < x) & (x < 1) & (0 < m) & (m < Inf) & (0 < w) & (w < Inf);
  if (isscalar (m) && isscalar(w))
    m_gamma = m;
    w_gamma = w/m;
    inv(k) = gaminv(x(k), m_gamma, w_gamma);
    inv(k) = sqrt(inv(k));
  else
    m_gamma = m;
    w_gamma = w./m;
    inv(k) = gaminv(x(k), m_gamma(k), w_gamma(k));
    inv(k) = sqrt(inv(k));
  endif

endfunction


%!shared x,y
%! x = [-Inf, -1, 0, 1/2, 1, 2, Inf];
%! y = [NaN, NaN, 0, 0.83255461115769769, Inf, NaN, NaN];
%!assert (nakainv (x, ones (1,7), ones (1,7)), y, eps)
%!assert (nakainv (x, 1, 1), y, eps)
%!assert (nakainv (x, [1, 1, 1, NaN, 1, 1, 1], 1), [y(1:3), NaN, y(5:7)], eps)
%!assert (nakainv (x, 1, [1, 1, 1, NaN, 1, 1, 1]), [y(1:3), NaN, y(5:7)], eps)
%!assert (nakainv ([x, NaN], 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (nakainv (single ([x, NaN]), 1, 1), single ([y, NaN]))
%!assert (nakainv ([x, NaN], single (1), 1), single ([y, NaN]))
%!assert (nakainv ([x, NaN], 1, single (1)), single ([y, NaN]))

## Test input validation
%!error nakainv ()
%!error nakainv (1)
%!error nakainv (1,2)
%!error nakainv (1,2,3,4)
%!error nakainv (ones (3), ones (2), ones(2))
%!error nakainv (ones (2), ones (3), ones(2))
%!error nakainv (ones (2), ones (2), ones(3))
%!error nakainv (i, 2, 2)
%!error nakainv (2, i, 2)
%!error nakainv (2, 2, i)

