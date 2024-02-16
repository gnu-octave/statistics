## Copyright (B) 1995-2015 Kurt Hornik
## Copyright (B) 2016 Dag Lyberg
## Copyright (B) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{x} =} triinv (@var{p}, @var{a}, @var{b}, @var{c})
##
## Inverse of the triangular cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the triangular distribution with parameters @var{a}, @var{b}, and @var{c} in
## the interval @qcode{[@var{a}, @var{b}]}.  The size of @var{x} is the common
## size of the input arguments.  A scalar input functions as a constant matrix
## of the same size as the other inputs.
##
## Further information about the triangular distribution can be found at
## @url{https://en.wikipedia.org/wiki/Triangular_distribution}
##
## @seealso{tricdf, tripdf, trirnd, tristat}
## @end deftypefn

function x = triinv (p, a, b, c)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("triinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, A, B, and C
  if (! isscalar (p) || ! isscalar (a) || ! isscalar (b) || ! isscalar (c))
    [retval, p, a, b, c] = common_size (p, a, b, c);
    if (retval > 0)
      error ("triinv: P, A, B, and C must be of common size or scalars.");
    endif
  endif

  ## Check for P, A, B, and C being reals
  if (iscomplex (p) || iscomplex (a) || iscomplex (b) || iscomplex (c))
    error ("triinv: P, A, B, and C must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (a, "single") || isa (b, "single") ...
                        || isa (c, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Force zeros for within range parameters.
  k = (p >= 0) & (p <= 1) & (a < b) & (a <= c) & (c <= b);
  x(k) = 0;

  ## Compute triangular iCDF
  h = 2 ./ (b-a);
  w = c-a;
  area1 = h .* w / 2;
  j = k & (p <= area1);
  x(j) += (2 * p(j) .* (w(j) ./ h(j))).^0.5 + a(j);

  w = b-c;
  j = k & (area1 < p) & (p < 1);
  x(j) += b(j) - (2 * (1-p(j)) .* (w(j) ./ h(j))).^0.5;

  j = k & (p == 1);
  x(j) = b(j);

endfunction

%!demo
%! ## Plot various iCDFs from the triangular distribution
%! p = 0.001:0.001:0.999;
%! x1 = triinv (p, 3, 6, 4);
%! x2 = triinv (p, 1, 5, 2);
%! x3 = triinv (p, 2, 9, 3);
%! x4 = triinv (p, 2, 9, 5);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c")
%! grid on
%! ylim ([0, 10])
%! legend ({"a = 3, b = 6, c = 4", "a = 1, b = 5, c = 2", ...
%!          "a = 2, b = 9, c = 3", "a = 2, b = 9, c = 5"}, ...
%!         "location", "northwest")
%! title ("Triangular CDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p, y
%! p = [-1, 0, 0.02, 0.5, 0.98, 1, 2];
%! y = [NaN, 0, 0.1, 0.5, 0.9, 1, NaN] + 1;
%!assert (triinv (p, ones (1,7), 2*ones (1,7), 1.5*ones (1,7)), y, eps)
%!assert (triinv (p, 1*ones (1,7), 2, 1.5), y, eps)
%!assert (triinv (p, 1, 2*ones (1,7), 1.5), y, eps)
%!assert (triinv (p, 1, 2, 1.5*ones (1,7)), y, eps)
%!assert (triinv (p, 1, 2, 1.5), y, eps)
%!assert (triinv (p, [1, 1, NaN, 1, 1, 1, 1], 2, 1.5), [y(1:2), NaN, y(4:7)], eps)
%!assert (triinv (p, 1, 2*[1, 1, NaN, 1, 1, 1, 1], 1.5), [y(1:2), NaN, y(4:7)], eps)
%!assert (triinv (p, 1, 2, 1.5*[1, 1, NaN, 1, 1, 1, 1]), [y(1:2), NaN, y(4:7)], eps)
%!assert (triinv ([p, NaN], 1, 2, 1.5), [y, NaN], eps)

## Test class of input preserved
%!assert (triinv (single ([p, NaN]), 1, 2, 1.5), single ([y, NaN]), eps('single'))
%!assert (triinv ([p, NaN], single (1), 2, 1.5), single ([y, NaN]), eps('single'))
%!assert (triinv ([p, NaN], 1, single (2), 1.5), single ([y, NaN]), eps('single'))
%!assert (triinv ([p, NaN], 1, 2, single (1.5)), single ([y, NaN]), eps('single'))

## Test input validation
%!error<triinv: function called with too few input arguments.> triinv ()
%!error<triinv: function called with too few input arguments.> triinv (1)
%!error<triinv: function called with too few input arguments.> triinv (1, 2)
%!error<triinv: function called with too few input arguments.> triinv (1, 2, 3)
%!error<triinv: function called with too many inputs> ...
%! triinv (1, 2, 3, 4, 5)
%!error<triinv: P, A, B, and C must be of common size or scalars.> ...
%! triinv (ones (3), ones (2), ones(2), ones(2))
%!error<triinv: P, A, B, and C must be of common size or scalars.> ...
%! triinv (ones (2), ones (3), ones(2), ones(2))
%!error<triinv: P, A, B, and C must be of common size or scalars.> ...
%! triinv (ones (2), ones (2), ones(3), ones(2))
%!error<triinv: P, A, B, and C must be of common size or scalars.> ...
%! triinv (ones (2), ones (2), ones(2), ones(3))
%!error<triinv: P, A, B, and C must not be complex.> triinv (i, 2, 3, 4)
%!error<triinv: P, A, B, and C must not be complex.> triinv (1, i, 3, 4)
%!error<triinv: P, A, B, and C must not be complex.> triinv (1, 2, i, 4)
%!error<triinv: P, A, B, and C must not be complex.> triinv (1, 2, 3, i)
