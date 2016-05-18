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
## @deftypefn {} {} gppdf (@var{x}, @var{location}, @var{scale}, @var{shape})
## Compute the probability density function (PDF) at @var{x} of the
## generalized Pareto distribution with parameters @var{location}, @var{scale},
## and @var{shape}.
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: PDF of the generalized Pareto distribution

function pdf = gppdf (x, location, scale, shape)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (location) || ! isscalar (scale) || ! isscalar (shape))
    [retval, x, location, scale, shape] = ...
      common_size (x, location, scale, shape);
    if (retval > 0)
      error ("gppdf: X, LOCATION, SCALE and SHAPE must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (location) || iscomplex (scale) ...
      || iscomplex (shape))
    error ("gppdf: X, LOCATION, SCALE and SHAPE must not be complex");
  endif

  if (isa (x, "single") || isa (location, "single") || isa (scale, "single") ...
      || isa (shape, "single"))
    pdf = zeros (size (x), "single");
  else
    pdf = zeros (size (x));
  endif

  k = isnan (x) | ! (-Inf < location) | ! (location < Inf) | ...
      ! (scale > 0) | ! (scale < Inf) | ! (-Inf < shape) | ! (shape < Inf);
  pdf(k) = NaN;

  k = (-Inf < x) & (x < Inf) & (-Inf < location) & (location < Inf) & ...
        (scale > 0) & (scale < Inf) & (-Inf < shape) & (shape < Inf);
  if (isscalar (location) && isscalar (scale) && isscalar (shape))
    z = (x - location) / scale;

    j = k & (shape == 0) & (z >= 0);
    if (any (j))
      pdf(j) = exp (-z(j));
    endif

    j = k & (shape > 0) & (z >= 0);
    if (any (j))
      pdf(j) = (shape * z(j) + 1).^(-(shape + 1) / shape);
    endif

    if (shape < 0)
      j = k & (shape < 0) & (0 <= z) & (z <= -1. / shape);
      if (any (j))
        pdf(j) = (shape * z(j) + 1).^(-(shape + 1) / shape);
      endif
    endif
  else
    z = (x - location) ./ scale;

    j = k & (shape == 0) & (z >= 0);
    if (any (j))
      pdf(j) = exp( -z(j));
    endif

    j = k & (shape > 0) & (z >= 0);
    if (any (j))
      pdf(j) = (shape(j) .* z(j) + 1).^(-(shape(j) + 1) ./ shape(j));
    endif

    if (any (shape < 0))
      j = k & (shape < 0) & (0 <= z) & (z <= -1 ./ shape);
      if (any (j))
        pdf(j) = (shape(j) .* z(j) + 1).^(-(shape(j) + 1) ./ shape(j));
      endif
    endif
  endif

endfunction


%!shared x,y1,y2,y3
%! x = [-Inf, -1, 0, 1/2, 1, Inf];
%! y1 = [0, 0, 1, 0.6065306597126334, 0.36787944117144233, 0];
%! y2 = [0, 0, 1, 4/9, 1/4, 0];
%! y3 = [0, 0, 1, 1, 1, 0];
%!assert (gppdf (x, zeros (1,6), ones (1,6), zeros (1,6)), y1, eps)
%!assert (gppdf (x, zeros (1,6), 1, 0), y1, eps)
%!assert (gppdf (x, 0, ones (1,6), 0), y1, eps)
%!assert (gppdf (x, 0, 1, zeros (1,6)), y1, eps)
%!assert (gppdf (x, 0, 1, 0), y1, eps)
%!assert (gppdf (x, [0, 0, 0, NaN, 0, 0], 1, 0), [y1(1:3), NaN, y1(5:6)])
%!assert (gppdf (x, 0, [1, 1, 1, NaN, 1, 1], 0), [y1(1:3), NaN, y1(5:6)])
%!assert (gppdf (x, 0, 1, [0, 0, 0, NaN, 0, 0]), [y1(1:3), NaN, y1(5:6)])
%!assert (gppdf ([x(1:3), NaN, x(5:6)], 0, 1, 0), [y1(1:3), NaN, y1(5:6)])

%!assert (gppdf (x, zeros (1,6), ones (1,6), ones (1,6)), y2, eps)
%!assert (gppdf (x, zeros (1,6), 1, 1), y2, eps)
%!assert (gppdf (x, 0, ones (1,6), 1), y2, eps)
%!assert (gppdf (x, 0, 1, ones (1,6)), y2, eps)
%!assert (gppdf (x, 0, 1, 1), y2, eps)
%!assert (gppdf (x, [0, 0, 0, NaN, 0, 0], 1, 1), [y2(1:3), NaN, y2(5:6)])
%!assert (gppdf (x, 0, [1, 1, 1, NaN, 1, 1], 1), [y2(1:3), NaN, y2(5:6)])
%!assert (gppdf (x, 0, 1, [1, 1, 1, NaN, 1, 1]), [y2(1:3), NaN, y2(5:6)])
%!assert (gppdf ([x(1:3), NaN, x(5:6)], 0, 1, 1), [y2(1:3), NaN, y2(5:6)])

%!assert (gppdf (x, zeros (1,6), ones (1,6), -ones (1,6)), y3, eps)
%!assert (gppdf (x, zeros (1,6), 1, -1), y3, eps)
%!assert (gppdf (x, 0, ones (1,6), -1), y3, eps)
%!assert (gppdf (x, 0, 1, -ones (1,6)), y3, eps)
%!assert (gppdf (x, 0, 1, -1), y3, eps)
%!assert (gppdf (x, [0, 0, 0, NaN, 0, 0], 1, -1), [y3(1:3), NaN, y3(5:6)])
%!assert (gppdf (x, 0, [1, 1, 1, NaN, 1, 1], -1), [y3(1:3), NaN, y3(5:6)])
%!assert (gppdf (x, 0, 1, [-1, -1, -1, NaN, -1, -1]), [y3(1:3), NaN, y3(5:6)])
%!assert (gppdf ([x(1:3), NaN, x(5:6)], 0, 1, -1), [y3(1:3), NaN, y3(5:6)])

## Test class of input preserved
%!assert (gppdf (single ([x, NaN]), 0, 1, 0), single ([y1, NaN]))
%!assert (gppdf ([x, NaN], single (0), 1, 0), single ([y1, NaN]))
%!assert (gppdf ([x, NaN], 0, single (1), 0), single ([y1, NaN]))
%!assert (gppdf ([x, NaN], 0, 1, single (0)), single ([y1, NaN]))

%!assert (gppdf (single ([x, NaN]), 0, 1, 1), single ([y2, NaN]))
%!assert (gppdf ([x, NaN], single (0), 1, 1), single ([y2, NaN]))
%!assert (gppdf ([x, NaN], 0, single (1), 1), single ([y2, NaN]))
%!assert (gppdf ([x, NaN], 0, 1, single (1)), single ([y2, NaN]))

%!assert (gppdf (single ([x, NaN]), 0, 1, -1), single ([y3, NaN]))
%!assert (gppdf ([x, NaN], single (0), 1, -1), single ([y3, NaN]))
%!assert (gppdf ([x, NaN], 0, single (1), -1), single ([y3, NaN]))
%!assert (gppdf ([x, NaN], 0, 1, single (-1)), single ([y3, NaN]))

## Test input validation
%!error gppdf ()
%!error gppdf (1)
%!error gppdf (1,2)
%!error gppdf (1,2,3)
%!error gppdf (1,2,3,4,5)
%!error gppdf (1, ones (3), ones (2), ones (2))
%!error gppdf (1, ones (2), ones (3), ones (2))
%!error gppdf (1, ones (2), ones (2), ones (3))
%!error gppdf (i, 2, 2, 2)
%!error gppdf (2, i, 2, 2)
%!error gppdf (2, 2, i, 2)
%!error gppdf (2, 2, 2, i)

