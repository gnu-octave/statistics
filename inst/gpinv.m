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
## @deftypefn {} {} gpinv (@var{x}, @var{location}, @var{scale}, @var{shape})
## For each element of @var{x}, compute the quantile (the inverse of the CDF)
## at @var{x} of the generalized Pareto distribution with parameters
## @var{location}, @var{scale}, and @var{shape}.
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: Quantile function of the generalized Pareto distribution

function inv = gpinv (x, location, scale, shape)
  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (location) || ! isscalar (scale) || ! isscalar (shape))
    [retval, x, location, scale, shape] = ...
        common_size (x, location, scale, shape);
    if (retval > 0)
      error ("gpinv: X, LOCATION, SCALE and SHAPE must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (location) ...
      || iscomplex (scale) || iscomplex (shape))
    error ("gpinv: X, LOCATION, SCALE and SHAPE must not be complex");
  endif

  if (isa (x, "single") || isa (location, "single") ...
      || isa (scale, "single") || isa (shape, "single"))
    inv = zeros (size (x), "single");
  else
    inv = zeros (size (x));
  endif

  k = isnan (x) | ! (0 <= x) | ! (x <= 1) ...
      | ! (-Inf < location) | ! (location < Inf) ...
      | ! (scale > 0) | ! (scale < Inf) ...
      | ! (-Inf < shape) | ! (shape < Inf);
  inv(k) = NaN;

  k = (0 <= x) & (x <= 1) & (-Inf < location) & (location < Inf) ...
      & (scale > 0) & (scale < Inf) & (-Inf < shape) & (shape < Inf);
  if (isscalar (location) && isscalar (scale) && isscalar (shape))
    if (shape == 0)
      inv(k) = -log(1 - x(k));
      inv(k) = scale * inv(k) + location;
    elseif (shape > 0)
      inv(k) = (1 - x(k)).^(-shape) - 1;
      inv(k) = (scale / shape) * inv(k) + location;
    elseif (shape < 0)
      inv(k) = (1 - x(k)).^(-shape) - 1;
      inv(k) = (scale / shape) * inv(k)  + location;
    end
  else
    j = k & (shape == 0);
    if (any (j))
      inv(j) = -log (1 - x(j));
      inv(j) = scale(j) .* inv(j) + location(j);
    endif
    
    j = k & (shape > 0);
    if (any (j))
      inv(j) = (1 - x(j)).^(-shape(j)) - 1;
      inv(j) = (scale(j) ./ shape(j)) .* inv(j) + location(j);
    endif
    
    j = k & (shape < 0);
    if (any (j))
      inv(j) = (1 - x(j)).^(-shape(j)) - 1;
      inv(j) = (scale(j) ./ shape(j)) .* inv(j) + location(j);
    endif
  endif
endfunction


%!shared x,y1,y2,y3
%! x = [-1, 0, 1/2, 1, 2];
%! y1 = [NaN, 0, 0.6931471805599453, Inf, NaN];
%! y2 = [NaN, 0, 1, Inf, NaN];
%! y3 = [NaN, 0, 1/2, 1, NaN];
%!assert (gpinv (x, zeros (1,5), ones (1,5), zeros (1,5)), y1)
%!assert (gpinv (x, zeros (1,5), 1, 0), y1)
%!assert (gpinv (x, 0, ones (1,5), 0), y1)
%!assert (gpinv (x, 0, 1, zeros (1,5)), y1)
%!assert (gpinv (x, 0, 1, 0), y1)
%!assert (gpinv (x, [0, 0, NaN, 0, 0], 1, 0), [y1(1:2), NaN, y1(4:5)])
%!assert (gpinv (x, 0, [1, 1, NaN, 1, 1], 0), [y1(1:2), NaN, y1(4:5)])
%!assert (gpinv (x, 0, 1, [0, 0, NaN, 0, 0]), [y1(1:2), NaN, y1(4:5)])
%!assert (gpinv ([x(1:2), NaN, x(4:5)], 0, 1, 0), [y1(1:2), NaN, y1(4:5)])

%!assert (gpinv (x, zeros (1,5), ones (1,5), ones (1,5)), y2)
%!assert (gpinv (x, zeros (1,5), 1, 1), y2)
%!assert (gpinv (x, 0, ones (1,5), 1), y2)
%!assert (gpinv (x, 0, 1, ones (1,5)), y2)
%!assert (gpinv (x, 0, 1, 1), y2)
%!assert (gpinv (x, [0, 0, NaN, 0, 0], 1, 1), [y2(1:2), NaN, y2(4:5)])
%!assert (gpinv (x, 0, [1, 1, NaN, 1, 1], 1), [y2(1:2), NaN, y2(4:5)])
%!assert (gpinv (x, 0, 1, [1, 1, NaN, 1, 1]), [y2(1:2), NaN, y2(4:5)])
%!assert (gpinv ([x(1:2), NaN, x(4:5)], 0, 1, 1), [y2(1:2), NaN, y2(4:5)])

%!assert (gpinv (x, zeros (1,5), ones (1,5), -ones (1,5)), y3)
%!assert (gpinv (x, zeros (1,5), 1, -1), y3)
%!assert (gpinv (x, 0, ones (1,5), -1), y3)
%!assert (gpinv (x, 0, 1, -ones (1,5)), y3)
%!assert (gpinv (x, 0, 1, -1), y3)
%!assert (gpinv (x, [0, 0, NaN, 0, 0], 1, -1), [y3(1:2), NaN, y3(4:5)])
%!assert (gpinv (x, 0, [1, 1, NaN, 1, 1], -1), [y3(1:2), NaN, y3(4:5)])
%!assert (gpinv (x, 0, 1, -[1, 1, NaN, 1, 1]), [y3(1:2), NaN, y3(4:5)])
%!assert (gpinv ([x(1:2), NaN, x(4:5)], 0, 1, -1), [y3(1:2), NaN, y3(4:5)])

## Test class of input preserved
%!assert (gpinv (single ([x, NaN]), 0, 1, 0), single ([y1, NaN]))
%!assert (gpinv ([x, NaN], single (0), 1, 0), single ([y1, NaN]))
%!assert (gpinv ([x, NaN], 0, single (1), 0), single ([y1, NaN]))
%!assert (gpinv ([x, NaN], 0, 1, single (0)), single ([y1, NaN]))

%!assert (gpinv (single ([x, NaN]), 0, 1, 1), single ([y2, NaN]))
%!assert (gpinv ([x, NaN], single (0), 1, 1), single ([y2, NaN]))
%!assert (gpinv ([x, NaN], 0, single (1), 1), single ([y2, NaN]))
%!assert (gpinv ([x, NaN], 0, 1, single (1)), single ([y2, NaN]))

%!assert (gpinv (single ([x, NaN]), 0, 1, -1), single ([y3, NaN]))
%!assert (gpinv ([x, NaN], single (0), 1, -1), single ([y3, NaN]))
%!assert (gpinv ([x, NaN], 0, single (1), -1), single ([y3, NaN]))
%!assert (gpinv ([x, NaN], 0, 1, single (-1)), single ([y3, NaN]))

## Test input validation
%!error gpinv ()
%!error gpinv (1)
%!error gpinv (1,2)
%!error gpinv (1,2,3)
%!error gpinv (1,2,3,4,5)
%!error gpinv (ones (3), ones (2), ones (2), ones (2))
%!error gpinv (ones (2), ones (3), ones (2), ones (2))
%!error gpinv (ones (2), ones (2), ones (3), ones (2))
%!error gpinv (ones (2), ones (2), ones (2), ones (3))
%!error gpinv (i, 2, 2, 2)
%!error gpinv (2, i, 2, 2)
%!error gpinv (2, 2, i, 2)
%!error gpinv (2, 2, 2, i)

