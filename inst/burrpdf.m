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
## @deftypefn {} {} burrpdf (@var{x}, @var{alpha}, @var{c}, @var{k})
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Burr distribution with scale parameter @var{alpha} and 
## shape parameters @var{c} and @var{k}.
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: PDF of the Burr distribution

function pdf = burrpdf (x, alpha, c, k)

  if (nargin != 4)
    print_usage ();
  endif

  if (! isscalar (alpha) || ! isscalar (c) || ! isscalar (k) )
    [retval, x, alpha, c, k] = common_size (x, alpha, c, k);
    if (retval > 0)
      error ("burrpdf: X, ALPHA, C AND K must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex(alpha) || iscomplex (c) || iscomplex (k))
    error ("burrpdf: X, ALPHA, C AND K must not be complex");
  endif

  if (isa (x, "single") || isa (alpha, "single") ...
      || isa (c, "single") || isa (k, "single"))
    pdf = zeros (size (x), "single");
  else
    pdf = zeros (size (x));
  endif

  j = isnan (x) | ! (alpha > 0) | ! (c > 0) | ! (k > 0);
  pdf(j) = NaN;

  j = (x > 0) & (0 < alpha) & (alpha < Inf) & (0 < c) & (c < Inf) ...
      & (0 < k) & (k < Inf);
  if (isscalar (alpha) && isscalar (c) && isscalar(k))
    pdf(j) = (c * k / alpha) .* (x(j) / alpha).^(c-1) ./ ...
        (1 + (x(j) / alpha).^c).^(k + 1);
  else
    pdf(j) = (c(j) .* k(j) ./ alpha(j) ).* x(j).^(c(j)-1) ./ ...
        (1 + (x(j) ./ alpha(j) ).^c(j) ).^(k(j) + 1);
  endif

endfunction


%!shared x,y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 1/4, 1/9, 0];
%!assert (burrpdf (x, ones(1,5), ones (1,5), ones (1,5)), y)
%!assert (burrpdf (x, 1, 1, 1), y)
%!assert (burrpdf (x, [1, 1, NaN, 1, 1], 1, 1), [y(1:2), NaN, y(4:5)])
%!assert (burrpdf (x, 1, [1, 1, NaN, 1, 1], 1), [y(1:2), NaN, y(4:5)])
%!assert (burrpdf (x, 1, 1, [1, 1, NaN, 1, 1]), [y(1:2), NaN, y(4:5)])
%!assert (burrpdf ([x, NaN], 1, 1, 1), [y, NaN])

## Test class of input preserved
%!assert (burrpdf (single ([x, NaN]), 1, 1, 1), single ([y, NaN]))
%!assert (burrpdf ([x, NaN], single (1), 1, 1), single ([y, NaN]))
%!assert (burrpdf ([x, NaN], 1, single (1), 1), single ([y, NaN]))
%!assert (burrpdf ([x, NaN], 1, 1, single (1)), single ([y, NaN]))

## Test input validation
%!error burrpdf ()
%!error burrpdf (1)
%!error burrpdf (1,2)
%!error burrpdf (1,2,3)
%!error burrpdf (1,2,3,4,5)
%!error burrpdf (ones (3), ones (2), ones(2), ones(2))
%!error burrpdf (ones (2), ones (3), ones(2), ones(2))
%!error burrpdf (ones (2), ones (2), ones(3), ones(2))
%!error burrpdf (ones (2), ones (2), ones(2), ones(3))
%!error burrpdf (i, 2, 2, 2)
%!error burrpdf (2, i, 2, 2)
%!error burrpdf (2, 2, i, 2)
%!error burrpdf (2, 2, 2, i)

