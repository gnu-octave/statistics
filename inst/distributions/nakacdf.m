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
## @deftypefn {} {} nakacdf (@var{x}, @var{m}, @var{w})
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the Nakagami distribution with shape parameter @var{m}
## and scale parameter @var{w}.
##
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: CDF of the Nakagami distribution

function cdf = nakacdf (x, m, w)

  if (nargin != 3)
    print_usage ();
  endif

  if (! isscalar (m) || ! isscalar (w))
    [retval, x, m, w] = common_size (x, m, w);
    if (retval > 0)
      error ("nakacdf: X, M and W must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (m) || iscomplex (w))
    error ("nakacdf: X, M and W must not be complex");
  endif

  if (isa (x, "single") || isa (m, "single") || isa (w, "single"))
    inv = zeros (size (x), "single");
  else
    inv = zeros (size (x));
  endif

  k = isnan (x) | ! (m > 0) | ! (w > 0);
  cdf(k) = NaN;

  k = (x == Inf) & (0 < m) & (m < Inf) & (0 < w) & (w < Inf);
  cdf(k) = 1;
  
  k = (0 < x) & (x < Inf) & (0 < m) & (m < Inf) & (0 < w) & (w < Inf);
  if (isscalar(x) && isscalar (m) && isscalar(w))
    left = m;
    right = (m/w) * x^2;
    cdf(k) = gammainc(right, left);
  elseif (isscalar (m) && isscalar(w))
    left = m * ones(size(x));
    right = (m/w) * x.^2;
    cdf(k) = gammainc(right(k), left(k));
  else
    left = m .* ones(size(x));
    right = (m./w) .* x.^2;
    cdf(k) = gammainc(right(k), left(k));
  endif

endfunction


%!shared x,y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 0.63212055882855778, 0.98168436111126578, 1];
%!assert (nakacdf (x, ones (1,5), ones (1,5)), y, eps)
%!assert (nakacdf (x, 1, 1), y, eps)
%!assert (nakacdf (x, [1, 1, NaN, 1, 1], 1), [y(1:2), NaN, y(4:5)])
%!assert (nakacdf (x, 1, [1, 1, NaN, 1, 1]), [y(1:2), NaN, y(4:5)])
%!assert (nakacdf ([x, NaN], 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (nakacdf (single ([x, NaN]), 1, 1), single ([y, NaN]), eps('single'))
%!assert (nakacdf ([x, NaN], single (1), 1), single ([y, NaN]), eps('single'))
%!assert (nakacdf ([x, NaN], 1, single (1)), single ([y, NaN]), eps('single'))

## Test input validation
%!error nakacdf ()
%!error nakacdf (1)
%!error nakacdf (1,2)
%!error nakacdf (1,2,3,4)
%!error nakacdf (ones (3), ones (2), ones(2))
%!error nakacdf (ones (2), ones (3), ones(2))
%!error nakacdf (ones (2), ones (2), ones(3))
%!error nakacdf (i, 2, 2)
%!error nakacdf (2, i, 2)
%!error nakacdf (2, 2, i)

