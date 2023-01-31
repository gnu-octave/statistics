## Copyright (C) 1997-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2018 John Donoghue
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
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
## @deftypefn  {statistics} @var{x} = gpinv (@var{p}, @var{shape}, @var{scale}, @var{location})
##
## Inverse of the generalized Pareto cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the generalized Pareto distribution with parameters
## @var{location}, @var{scale}, and @var{shape}.  The size of @var{x} is the
## common size of the input arguments.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## @seealso{gpcdf, gppdf, gprnd, gpfit, gplike, gpstat}
## @end deftypefn

function x = gpinv (p, shape, scale, location)
  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (location) || ! isscalar (scale) || ! isscalar (shape))
    [retval, p, location, scale, shape] = ...
        common_size (p, location, scale, shape);
    if (retval > 0)
      error (strcat (["gpinv: X, SHAPE, SCALE, and LOCATION must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  if (iscomplex (p) || iscomplex (location) ...
      || iscomplex (scale) || iscomplex (shape))
    error ("gpinv: X, SHAPE, SCALE, and LOCATION must not be complex.");
  endif

  if (isa (p, "single") || isa (location, "single") ...
      || isa (scale, "single") || isa (shape, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  k = isnan (p) | ! (0 <= p) | ! (p <= 1) ...
      | ! (-Inf < location) | ! (location < Inf) ...
      | ! (scale > 0) | ! (scale < Inf) ...
      | ! (-Inf < shape) | ! (shape < Inf);
  x(k) = NaN;

  k = (0 <= p) & (p <= 1) & (-Inf < location) & (location < Inf) ...
      & (scale > 0) & (scale < Inf) & (-Inf < shape) & (shape < Inf);
  if (isscalar (location) && isscalar (scale) && isscalar (shape))
    if (shape == 0)
      x(k) = -log(1 - p(k));
      x(k) = scale * x(k) + location;
    elseif (shape > 0)
      x(k) = (1 - p(k)).^(-shape) - 1;
      x(k) = (scale / shape) * x(k) + location;
    elseif (shape < 0)
      x(k) = (1 - p(k)).^(-shape) - 1;
      x(k) = (scale / shape) * x(k)  + location;
    end
  else
    j = k & (shape == 0);
    if (any (j))
      x(j) = -log (1 - p(j));
      x(j) = scale(j) .* x(j) + location(j);
    endif

    j = k & (shape > 0);
    if (any (j))
      x(j) = (1 - p(j)).^(-shape(j)) - 1;
      x(j) = (scale(j) ./ shape(j)) .* x(j) + location(j);
    endif

    j = k & (shape < 0);
    if (any (j))
      x(j) = (1 - p(j)).^(-shape(j)) - 1;
      x(j) = (scale(j) ./ shape(j)) .* x(j) + location(j);
    endif
  endif
endfunction


%!shared p,y1,y2,y3
%! p = [-1, 0, 1/2, 1, 2];
%! y1 = [NaN, 0, 0.6931471805599453, Inf, NaN];
%! y2 = [NaN, 0, 1, Inf, NaN];
%! y3 = [NaN, 0, 1/2, 1, NaN];
%!assert (gpinv (p, zeros (1,5), ones (1,5), zeros (1,5)), y1)
%!assert (gpinv (p, 0, 1, zeros (1,5)), y1)
%!assert (gpinv (p, 0, ones (1,5), 0), y1)
%!assert (gpinv (p, zeros (1,5), 1, 0), y1)
%!assert (gpinv (p, 0, 1, 0), y1)
%!assert (gpinv (p, 0, 1, [0, 0, NaN, 0, 0]), [y1(1:2), NaN, y1(4:5)])
%!assert (gpinv (p, 0, [1, 1, NaN, 1, 1], 0), [y1(1:2), NaN, y1(4:5)])
%!assert (gpinv (p, [0, 0, NaN, 0, 0], 1, 0), [y1(1:2), NaN, y1(4:5)])
%!assert (gpinv ([p(1:2), NaN, p(4:5)], 0, 1, 0), [y1(1:2), NaN, y1(4:5)])

%!assert (gpinv (p, ones (1,5), ones (1,5), zeros (1,5)), y2)
%!assert (gpinv (p, 1, 1, zeros (1,5)), y2)
%!assert (gpinv (p, 1, ones (1,5), 0), y2)
%!assert (gpinv (p, ones (1,5), 1, 0), y2)
%!assert (gpinv (p, 1, 1, 0), y2)
%!assert (gpinv (p, 1, 1, [0, 0, NaN, 0, 0]), [y2(1:2), NaN, y2(4:5)])
%!assert (gpinv (p, 1, [1, 1, NaN, 1, 1], 0), [y2(1:2), NaN, y2(4:5)])
%!assert (gpinv (p, [1, 1, NaN, 1, 1], 1, 0), [y2(1:2), NaN, y2(4:5)])
%!assert (gpinv ([p(1:2), NaN, p(4:5)], 1, 1, 0), [y2(1:2), NaN, y2(4:5)])

%!assert (gpinv (p, -ones (1,5), ones (1,5), zeros (1,5)), y3)
%!assert (gpinv (p, -1, 1, zeros (1,5)), y3)
%!assert (gpinv (p, -1, ones (1,5), 0), y3)
%!assert (gpinv (p, -ones (1,5), 1, 0), y3)
%!assert (gpinv (p, -1, 1, 0), y3)
%!assert (gpinv (p, -1, 1, [0, 0, NaN, 0, 0]), [y3(1:2), NaN, y3(4:5)])
%!assert (gpinv (p, -1, [1, 1, NaN, 1, 1], 0), [y3(1:2), NaN, y3(4:5)])
%!assert (gpinv (p, -[1, 1, NaN, 1, 1], 1, 0), [y3(1:2), NaN, y3(4:5)])
%!assert (gpinv ([p(1:2), NaN, p(4:5)], -1, 1, 0), [y3(1:2), NaN, y3(4:5)])

## Test class of input preserved
%!assert (gpinv (single ([p, NaN]), 0, 1, 0), single ([y1, NaN]))
%!assert (gpinv ([p, NaN], 0, 1, single (0)), single ([y1, NaN]))
%!assert (gpinv ([p, NaN], 0, single (1), 0), single ([y1, NaN]))
%!assert (gpinv ([p, NaN], single (0), 1, 0), single ([y1, NaN]))

%!assert (gpinv (single ([p, NaN]), 1, 1, 0), single ([y2, NaN]))
%!assert (gpinv ([p, NaN], 1, 1, single (0)), single ([y2, NaN]))
%!assert (gpinv ([p, NaN], 1, single (1), 0), single ([y2, NaN]))
%!assert (gpinv ([p, NaN], single (1), 1, 0), single ([y2, NaN]))

%!assert (gpinv (single ([p, NaN]), -1, 1, 0), single ([y3, NaN]))
%!assert (gpinv ([p, NaN], -1, 1, single (0)), single ([y3, NaN]))
%!assert (gpinv ([p, NaN], -1, single (1), 0), single ([y3, NaN]))
%!assert (gpinv ([p, NaN], single (-1), 1, 0), single ([y3, NaN]))

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

