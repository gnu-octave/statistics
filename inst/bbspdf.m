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
## @deftypefn {} {} bbspdf (@var{x}, @var{location}, @var{scale}, @var{shape})
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Birnbaum-Saunders distribution with parameters
## @var{location}, @var{scale} and @var{shape}.
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: PDF of the Birnbaum-Saunders distribution

function pdf = bbspdf (x, location, scale, shape)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (location) || ! isscalar (scale) || ! isscalar(shape))
    [retval, x, location, scale, shape] = ...
        common_size (x, location, scale, shape);
    if (retval > 0)
      error ("bbspdf: X, LOCATION, SCALE and SHAPE must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (location) ...
      || iscomplex (scale) || iscomplex(shape))
    error ("bbspdf: X, LOCATION, SCALE and SHAPE must not be complex");
  endif

  if (isa (x, "single") || isa (location, "single") || isa (scale, "single") ...
      || isa (shape, "single"))
    pdf = zeros (size (x), "single");
  else
    pdf = zeros (size (x));
  endif

  k = isnan (x) | ! (-Inf < location) | ! (location < Inf) ...
      | ! (scale > 0) | ! (scale < Inf) ...
      | ! (shape > 0) | ! (shape < Inf);

  pdf(k) = NaN;

  k = (x > location) & (x < Inf) & (-Inf < location) ...
      & (location < Inf) & (0 < scale) & (scale < Inf) ...
      & (0 < shape) & (shape < Inf);

  if (isscalar (location) && isscalar(scale) && isscalar(shape))
    a = x(k) - location;
    b = sqrt(a ./ scale);
    pdf(k) = ((b + b.^-1) ./ (2 * shape * a)) ...
        .* normpdf ((b - b.^-1) / shape);
  else
    a = x(k) - location(k);
    b = sqrt(a ./ scale(k));
    pdf(k) = ((b + b.^-1) ./ (2 * shape(k).* a)) ...
        .* normpdf ((b - b.^-1) ./ shape(k));
  endif
endfunction


%!shared x,y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 0.3989422804014327, 0.1647717335503959, 0];
%!assert (bbspdf (x, zeros (1,5), ones (1,5), ones (1,5)), y, eps)
%!assert (bbspdf (x, zeros (1,5), 1, 1), y, eps)
%!assert (bbspdf (x, 0, ones (1,5), 1), y, eps)
%!assert (bbspdf (x, 0, 1, ones (1,5)), y, eps)
%!assert (bbspdf (x, 0, 1, 1), y, eps)
%!assert (bbspdf (x, [0, 0, NaN, 0, 0], 1, 1), [y(1:2), NaN, y(4:5)], eps)
%!assert (bbspdf (x, 0, [1, 1, NaN, 1, 1], 1), [y(1:2), NaN, y(4:5)], eps)
%!assert (bbspdf (x, 0, 1, [1, 1, NaN, 1, 1]), [y(1:2), NaN, y(4:5)], eps)
%!assert (bbspdf ([x, NaN], 0, 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (bbspdf (single ([x, NaN]), 0, 1, 1), single ([y, NaN]), eps('single'))
%!assert (bbspdf ([x, NaN], single (0), 1, 1), single ([y, NaN]), eps('single'))
%!assert (bbspdf ([x, NaN], 0, single (1), 1), single ([y, NaN]), eps('single'))
%!assert (bbspdf ([x, NaN], 0, 1, single (1)), single ([y, NaN]), eps('single'))

## Test input validation
%!error bbspdf ()
%!error bbspdf (1)
%!error bbspdf (1,2,3)
%!error bbspdf (1,2,3,4,5)
%!error bbspdf (ones (3), ones (2), ones(2), ones(2))
%!error bbspdf (ones (2), ones (3), ones(2), ones(2))
%!error bbspdf (ones (2), ones (2), ones(3), ones(2))
%!error bbspdf (ones (2), ones (2), ones(2), ones(3))
%!error bbspdf (i, 2, 3, 4)
%!error bbspdf (1, i, 3, 4)
%!error bbspdf (1, 2, i, 4)
%!error bbspdf (1, 2, 3, i)

