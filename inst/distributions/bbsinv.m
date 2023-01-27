## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2018 John Donoghue
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
## @deftypefn  {statistics} @var{x} = bbsinv (@var{p}, @var{shape}, @var{scale}, @var{location})
##
## Inverse of the Birnbaum-Saunders cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Birnbaum-Saunders distribution with parameters @var{shape},
## @var{scale}, and @var{location}.  The size of @var{x} is the common size of
## @var{p}, @var{shape}, @var{scale}, and @var{location}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## @seealso{bbscdf, bbspdf, bbsrnd}
## @end deftypefn

function x = bbsinv (p, shape, scale, location)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (location) || ! isscalar (scale) || ! isscalar(shape))
    [retval, p, location, scale, shape] = ...
        common_size (p, location, scale, shape);
    if (retval > 0)
      error (strcat (["bbsinv: P, SHAPE, SCALE, and LOCATION must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  if (iscomplex (p) || iscomplex (location) ...
      || iscomplex (scale) || iscomplex(shape))
    error ("bbsinv: P, SHAPE, SCALE, and LOCATION must not be complex.");
  endif

  if (isa (p, "single") || isa (location, "single") ...
      || isa (scale, "single") || isa (shape, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  k = isnan (p) | (p < 0) | (p > 1) | ! (-Inf < location) | ! (location < Inf) ...
      | ! (scale > 0) | ! (scale < Inf) | ! (shape > 0) | ! (shape < Inf);
  x(k) = NaN;

  k = (p <= 0) & (-Inf < location) & (location < Inf) ...
      & (scale > 0) & (scale < Inf) & (shape > 0) & (shape < Inf);
  x(k) = 0;

  k = (p == 1) & (-Inf < location) & (location < Inf) ...
      & (scale > 0) & (scale < Inf) & (shape > 0) & (shape < Inf);
  x(k) = Inf;

  k = (0 < p) & (p < 1) & (location < Inf) & (0 < scale) & (scale < Inf) ...
    & (0 < shape) & (shape < Inf);
  if (isscalar (location) && isscalar(scale) && isscalar(shape))
    y = shape * norminv (p(k));
    x(k) = location + scale * (y + sqrt (4 + y.^2)).^2 / 4;
  else
    y = shape(k) .* norminv (p(k));
    x(k) = location(k) + scale(k) .* (y + sqrt (4 + y.^2)).^2 ./ 4;
  endif

endfunction


%!shared p,y,f
%! f = @(p,a,b,c) (a + b * (c * norminv (p) + sqrt (4 + (c * norminv(p))^2))^2) / 4;
%! p = [-1, 0, 1/4, 1/2, 1, 2];
%! y = [0, 0, f(1/4, 0, 1, 1), 1, Inf, NaN];
%!assert (bbsinv (p, ones (1,6), ones (1,6), zeros (1,6)), y)
%!assert (bbsinv (p, 1, 1, zeros (1,6)), y)
%!assert (bbsinv (p, 1, ones (1,6), 0), y)
%!assert (bbsinv (p, ones (1,6), 1, 0), y)
%!assert (bbsinv (p, 1, 1, 0), y)
%!assert (bbsinv (p, 1, 1, [0, 0, 0, NaN, 0, 0]), [y(1:3), NaN, y(5:6)])
%!assert (bbsinv (p, 1, [1, 1, 1, NaN, 1, 1], 0), [y(1:3), NaN, y(5:6)])
%!assert (bbsinv (p, [1, 1, 1, NaN, 1, 1], 1, 0), [y(1:3), NaN, y(5:6)])
%!assert (bbsinv ([p, NaN], 1, 1, 0), [y, NaN])

## Test class of input preserved
%!assert (bbsinv (single ([p, NaN]), 1, 1, 0), single ([y, NaN]))
%!assert (bbsinv ([p, NaN], 1, 1, single (0)), single ([y, NaN]))
%!assert (bbsinv ([p, NaN], 1, single (1), 0), single ([y, NaN]))
%!assert (bbsinv ([p, NaN], single (1), 1, 0), single ([y, NaN]))

## Test input validation
%!error bbsinv ()
%!error bbsinv (1)
%!error bbsinv (1,2,3)
%!error bbsinv (1,2,3,4,5)
%!error bbsinv (ones (3), ones (2), ones(2), ones(2))
%!error bbsinv (ones (2), ones (3), ones(2), ones(2))
%!error bbsinv (ones (2), ones (2), ones(3), ones(2))
%!error bbsinv (ones (2), ones (2), ones(2), ones(3))
%!error bbsinv (i, 4, 3, 2)
%!error bbsinv (1, i, 3, 2)
%!error bbsinv (1, 4, i, 2)
%!error bbsinv (1, 4, 3, i)

