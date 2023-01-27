## Copyright (C) 2018 John Donoghue
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1995-2015 Kurt Hornik
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
## @deftypefn  {statistics} @var{y} = bbspdf (@var{x}, @var{shape}, @var{scale}, @var{location})
##
## Birnbaum-Saunders probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Birnbaum-Saunders distribution with parameters @var{shape},
## @var{scale}, and @var{location}.  The size of @var{y} is the common size of
## @var{x}, @var{shape}, @var{scale}, and @var{location}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## @seealso{bbscdf, bbsinv, bbsrnd}
## @end deftypefn

function y = bbspdf (x, shape, scale, location)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (location) || ! isscalar (scale) || ! isscalar(shape))
    [retval, x, location, scale, shape] = ...
        common_size (x, location, scale, shape);
    if (retval > 0)
      error (strcat (["bbspdf: X, SHAPE, SCALE and LOCATION must be of"], ...
                     [" common size or scalars"]));
    endif
  endif

  if (iscomplex (x) || iscomplex (location) ...
      || iscomplex (scale) || iscomplex(shape))
    error ("bbspdf: X, SHAPE, SCALE and LOCATION must not be complex");
  endif

  if (isa (x, "single") || isa (location, "single") || isa (scale, "single") ...
      || isa (shape, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  k = isnan (x) | ! (-Inf < location) | ! (location < Inf) ...
      | ! (scale > 0) | ! (scale < Inf) ...
      | ! (shape > 0) | ! (shape < Inf);

  y(k) = NaN;

  k = (x > location) & (x < Inf) & (-Inf < location) ...
      & (location < Inf) & (0 < scale) & (scale < Inf) ...
      & (0 < shape) & (shape < Inf);

  if (isscalar (location) && isscalar(scale) && isscalar(shape))
    a = x(k) - location;
    b = sqrt(a ./ scale);
    y(k) = ((b + b.^-1) ./ (2 * shape * a)) ...
        .* normpdf ((b - b.^-1) / shape);
  else
    a = x(k) - location(k);
    b = sqrt(a ./ scale(k));
    y(k) = ((b + b.^-1) ./ (2 * shape(k).* a)) ...
        .* normpdf ((b - b.^-1) ./ shape(k));
  endif
endfunction


%!shared x,y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 0.3989422804014327, 0.1647717335503959, 0];
%!assert (bbspdf (x, ones (1,5), ones (1,5), zeros (1,5)), y, eps)
%!assert (bbspdf (x, 1, 1, zeros (1,5)), y, eps)
%!assert (bbspdf (x, 1, ones (1,5), 0), y, eps)
%!assert (bbspdf (x, ones (1,5), 1, 0), y, eps)
%!assert (bbspdf (x, 1, 1, 0), y, eps)
%!assert (bbspdf (x, 1, 1, [0, 0, NaN, 0, 0]), [y(1:2), NaN, y(4:5)], eps)
%!assert (bbspdf (x, 1, [1, 1, NaN, 1, 1], 0), [y(1:2), NaN, y(4:5)], eps)
%!assert (bbspdf (x, [1, 1, NaN, 1, 1], 1, 0), [y(1:2), NaN, y(4:5)], eps)
%!assert (bbspdf ([x, NaN], 1, 1, 0), [y, NaN], eps)

## Test class of input preserved
%!assert (bbspdf (single ([x, NaN]), 1, 1, 0), single ([y, NaN]), eps('single'))
%!assert (bbspdf ([x, NaN], 1, 1, single (0)), single ([y, NaN]), eps('single'))
%!assert (bbspdf ([x, NaN], 1, single (1), 0), single ([y, NaN]), eps('single'))
%!assert (bbspdf ([x, NaN], single (1), 1, 0), single ([y, NaN]), eps('single'))

## Test input validation
%!error bbspdf ()
%!error bbspdf (1)
%!error bbspdf (1,2,3)
%!error bbspdf (1,2,3,4,5)
%!error bbspdf (ones (3), ones (2), ones(2), ones(2))
%!error bbspdf (ones (2), ones (3), ones(2), ones(2))
%!error bbspdf (ones (2), ones (2), ones(3), ones(2))
%!error bbspdf (ones (2), ones (2), ones(2), ones(3))
%!error bbspdf (i, 4, 3, 2)
%!error bbspdf (1, i, 3, 2)
%!error bbspdf (1, 4, i, 2)
%!error bbspdf (1, 4, 3, i)

