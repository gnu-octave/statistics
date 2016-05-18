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
## @deftypefn {} {} nakapdf (@var{x}, @var{m}, @var{w})
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Nakagami distribution with shape parameter @var{m} and
## scale parameter @var{w}.
## @end deftypefn

## Author: Dag Lyberg <daglyberg80@gmail.com>
## Description: PDF of the Nakagami distribution

function pdf = nakapdf (x, m, w)

  if (nargin != 3)
    print_usage ();
  endif

  if (! isscalar (m) || ! isscalar (w))
    [retval, x, m, w] = common_size (x, m, w);
    if (retval > 0)
      error ("nakapdf: X, M and W must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (m) || iscomplex (w))
    error ("nakapdf: X, M and W must not be complex");
  endif

  if (isa (x, "single") || isa (m, "single") || isa (w, "single"))
    pdf = zeros (size (x), "single");
  else
    pdf = zeros (size (x));
  endif

  k = isnan (x) | ! (m > 0.5) | ! (w > 0);
  pdf(k) = NaN;

  k = (0 < x) & (x < Inf) & (0 < m) & (m < Inf) & (0 < w) & (w < Inf);
  if (isscalar (m) && isscalar(w))
    pdf(k) = exp (log (2) + m*log (m) - log (gamma (m)) - ...
               m*log (w) + (2*m-1) * ...
               log (x(k)) - (m/w) * x(k).^2);
  else
    pdf(k) = exp(log(2) + m(k).*log (m(k)) - log (gamma (m(k))) - ...
               m(k).*log (w(k)) + (2*m(k)-1) ...
               .* log (x(k)) - (m(k)./w(k)) .* x(k).^2);
  endif

endfunction


%!shared x,y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 0.73575888234288467, 0.073262555554936715, 0];
%!assert (nakapdf (x, ones (1,5), ones (1,5)), y, eps)
%!assert (nakapdf (x, 1, 1), y, eps)
%!assert (nakapdf (x, [1, 1, NaN, 1, 1], 1), [y(1:2), NaN, y(4:5)], eps)
%!assert (nakapdf (x, 1, [1, 1, NaN, 1, 1]), [y(1:2), NaN, y(4:5)], eps)
%!assert (nakapdf ([x, NaN], 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (nakapdf (single ([x, NaN]), 1, 1), single ([y, NaN]))
%!assert (nakapdf ([x, NaN], single (1), 1), single ([y, NaN]))
%!assert (nakapdf ([x, NaN], 1, single (1)), single ([y, NaN]))

## Test input validation
%!error nakapdf ()
%!error nakapdf (1)
%!error nakapdf (1,2)
%!error nakapdf (1,2,3,4)
%!error nakapdf (ones (3), ones (2), ones(2))
%!error nakapdf (ones (2), ones (3), ones(2))
%!error nakapdf (ones (2), ones (2), ones(3))
%!error nakapdf (i, 2, 2)
%!error nakapdf (2, i, 2)
%!error nakapdf (2, 2, i)

