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
## @deftypefn {} {} bbsinv (@var{x}, @var{location}, @var{scale}, @var{shape})
## For each element of @var{x}, compute the quantile (the inverse of the CDF)
## at @var{x} of the Birnbaum-Saunders distribution with parameters
## @var{location}, @var{scale}, and @var{shape}.
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: Quantile function of the Birnbaum-Saunders distribution

function inv = bbsinv (x, location, scale, shape)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (location) || ! isscalar (scale) || ! isscalar(shape))
    [retval, x, location, scale, shape] = ...
        common_size (x, location, scale, shape);
    if (retval > 0)
      error ("bbsinv: X, LOCATION, SCALE and SHAPE must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (location) ...
      || iscomplex (scale) || iscomplex(shape))
    error ("bbsinv: X, LOCATION, SCALE and SHAPE must not be complex");
  endif

  if (isa (x, "single") || isa (location, "single") ...
      || isa (scale, "single") || isa (shape, "single"))
    inv = zeros (size (x), "single");
  else
    inv = zeros (size (x));
  endif

  k = isnan (x) | (x < 0) | (x > 1) | ! (-Inf < location) | ! (location < Inf) ...
      | ! (scale > 0) | ! (scale < Inf) | ! (shape > 0) | ! (shape < Inf);
  inv(k) = NaN;

  k = (x <= 0) & (-Inf < location) & (location < Inf) ...
      & (scale > 0) & (scale < Inf) & (shape > 0) & (shape < Inf);
  inv(k) = 0;

  k = (x == 1) & (-Inf < location) & (location < Inf) ...
      & (scale > 0) & (scale < Inf) & (shape > 0) & (shape < Inf);
  inv(k) = Inf;

  k = (0 < x) & (x < 1) & (location < Inf) & (0 < scale) & (scale < Inf) ...
    & (0 < shape) & (shape < Inf);
  if (isscalar (location) && isscalar(scale) && isscalar(shape))
    y = shape * norminv (x(k));
    inv(k) = location + scale * (y + sqrt (4 + y.^2)).^2 / 4;
  else
    y = shape(k) .* norminv (x(k));
    inv(k) = location(k) + scale(k) .* (y + sqrt (4 + y.^2)).^2 ./ 4;
  endif

endfunction


%!shared x,y,f
%! f = @(x,a,b,c) (a + b * (c * norminv (x) + sqrt (4 + (c * norminv(x))^2))^2) / 4;
%! x = [-1, 0, 1/4, 1/2, 1, 2];
%! y = [0, 0, f(1/4, 0, 1, 1), 1, Inf, NaN];
%!assert (bbsinv (x, zeros (1,6), ones (1,6), ones (1,6)), y)
%!assert (bbsinv (x, zeros (1,6), 1, 1), y)
%!assert (bbsinv (x, 0, ones (1,6), 1), y)
%!assert (bbsinv (x, 0, 1, ones (1,6)), y)
%!assert (bbsinv (x, 0, 1, 1), y)
%!assert (bbsinv (x, [0, 0, 0, NaN, 0, 0], 1, 1), [y(1:3), NaN, y(5:6)])
%!assert (bbsinv (x, 0, [1, 1, 1, NaN, 1, 1], 1), [y(1:3), NaN, y(5:6)])
%!assert (bbsinv (x, 0, 1, [1, 1, 1, NaN, 1, 1]), [y(1:3), NaN, y(5:6)])
%!assert (bbsinv ([x, NaN], 0, 1, 1), [y, NaN])

## Test class of input preserved
%!assert (bbsinv (single ([x, NaN]), 0, 1, 1), single ([y, NaN]))
%!assert (bbsinv ([x, NaN], single (0), 1, 1), single ([y, NaN]))
%!assert (bbsinv ([x, NaN], 0, single (1), 1), single ([y, NaN]))
%!assert (bbsinv ([x, NaN], 0, 1, single (1)), single ([y, NaN]))

## Test input validation
%!error bbsinv ()
%!error bbsinv (1)
%!error bbsinv (1,2,3)
%!error bbsinv (1,2,3,4,5)
%!error bbsinv (ones (3), ones (2), ones(2), ones(2))
%!error bbsinv (ones (2), ones (3), ones(2), ones(2))
%!error bbsinv (ones (2), ones (2), ones(3), ones(2))
%!error bbsinv (ones (2), ones (2), ones(2), ones(3))
%!error bbsinv (i, 2, 3, 4)
%!error bbsinv (1, i, 3, 4)
%!error bbsinv (1, 2, i, 4)
%!error bbsinv (1, 2, 3, i)

