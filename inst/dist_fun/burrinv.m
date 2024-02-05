## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
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
## @deftypefn  {statistics} {@var{x} =} burrinv (@var{p}, @var{lambda}, @var{c}, @var{k})
##
## Inverse of the Burr type XII cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the Burr type XII distribution with scale parameter @var{lambda}, first shape
## parameter @var{c}, and second shape parameter @var{k}.  The size of @var{x}
## is the common size of @var{p}, @var{lambda}, @var{c}, and @var{k}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## Further information about the Burr distribution can be found at
## @url{https://en.wikipedia.org/wiki/Burr_distribution}
##
## @seealso{burrcdf, burrpdf, burrrnd, burrfit, burrlike}
## @end deftypefn

function x = burrinv (p, lambda, c, k)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("burrinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, LANBDA, C and K
  if (! isscalar (p) || ! isscalar (lambda) || ! isscalar (c) || ! isscalar (k))
    [retval, p, lambda, c, k] = common_size (p, lambda, c, k);
    if (retval > 0)
      error ("burrinv: P, LAMBDA, C, and K must be of common size or scalars.");
    endif
  endif

  ## Check for P, LANBDA, C and K being reals
  if (iscomplex (p) || iscomplex (lambda) || iscomplex (c) || iscomplex (k))
    error ("burrinv: P, LAMBDA, C, and K must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (lambda, "single") || isa (c, "single") ...
                        || isa (k, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  ## Force NaNs for out of range parameters
  j = isnan (p) | (p < 0) | (p > 1) | ! (lambda > 0) | ! (c > 0) | ! (k > 0);
  x(j) = NaN;

  ## Handle edge cases
  j = (p == 1) & (lambda > 0) & (lambda < Inf) & (c > 0) & (c < Inf) ...
               & (k > 0) & (k < Inf);
  x(j) = Inf;

  ## Handle all other valid cases
  j = (0 < p) & (p < 1) & (0 < lambda) & (lambda < Inf) & (0 < c) & (c < Inf) ...
      & (0 < k) & (k < Inf);
  if (isscalar (lambda) && isscalar(c) && isscalar(k))
    x(j) = ((1 - p(j) / lambda).^(-1 / k) - 1).^(1 / c) ;
  else
    x(j) = ((1 - p(j) ./ lambda(j)).^(-1 ./ k(j)) - 1).^(1 ./ c(j)) ;
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the Burr type XII distribution
%! p = 0.001:0.001:0.999;
%! x1 = burrinv (p, 1, 1, 1);
%! x2 = burrinv (p, 1, 1, 2);
%! x3 = burrinv (p, 1, 1, 3);
%! x4 = burrinv (p, 1, 2, 1);
%! x5 = burrinv (p, 1, 3, 1);
%! x6 = burrinv (p, 1, 0.5, 2);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", ...
%!       p, x4, "-c", p, x5, "-m", p, x6, "-k")
%! grid on
%! ylim ([0, 5])
%! legend ({"λ = 1, c = 1, k = 1", "λ = 1, c = 1, k = 2", ...
%!          "λ = 1, c = 1, k = 3", "λ = 1, c = 2, k = 1", ...
%!          "λ = 1, c = 3, k = 1", "λ = 1, c = 0.5, k = 2"}, ...
%!         "location", "northwest")
%! title ("Burr type XII iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p, y
%! p = [-Inf, -1, 0, 1/2, 1, 2, Inf];
%! y = [NaN, NaN, 0, 1 , Inf, NaN, NaN];
%!assert (burrinv (p, ones (1,7), ones (1,7), ones(1,7)), y, eps)
%!assert (burrinv (p, 1, 1, 1), y, eps)
%!assert (burrinv (p, [1, 1, 1, NaN, 1, 1, 1], 1, 1), [y(1:3), NaN, y(5:7)], eps)
%!assert (burrinv (p, 1, [1, 1, 1, NaN, 1, 1, 1], 1), [y(1:3), NaN, y(5:7)], eps)
%!assert (burrinv (p, 1, 1, [1, 1, 1, NaN, 1, 1, 1]), [y(1:3), NaN, y(5:7)], eps)
%!assert (burrinv ([p, NaN], 1, 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (burrinv (single ([p, NaN]), 1, 1, 1), single ([y, NaN]), eps("single"))
%!assert (burrinv ([p, NaN], single (1), 1, 1), single ([y, NaN]), eps("single"))
%!assert (burrinv ([p, NaN], 1, single (1), 1), single ([y, NaN]), eps("single"))
%!assert (burrinv ([p, NaN], 1, 1, single (1)), single ([y, NaN]), eps("single"))

## Test input validation
%!error<burrinv: function called with too few input arguments.> burrinv ()
%!error<burrinv: function called with too few input arguments.> burrinv (1)
%!error<burrinv: function called with too few input arguments.> burrinv (1, 2)
%!error<burrinv: function called with too few input arguments.> burrinv (1, 2, 3)
%!error<burrinv: function called with too many inputs> ...
%! burrinv (1, 2, 3, 4, 5)
%!error<burrinv: P, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrinv (ones (3), ones (2), ones(2), ones(2))
%!error<burrinv: P, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrinv (ones (2), ones (3), ones(2), ones(2))
%!error<burrinv: P, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrinv (ones (2), ones (2), ones(3), ones(2))
%!error<burrinv: P, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrinv (ones (2), ones (2), ones(2), ones(3))
%!error<burrinv: P, LAMBDA, C, and K must not be complex.> burrinv (i, 2, 3, 4)
%!error<burrinv: P, LAMBDA, C, and K must not be complex.> burrinv (1, i, 3, 4)
%!error<burrinv: P, LAMBDA, C, and K must not be complex.> burrinv (1, 2, i, 4)
%!error<burrinv: P, LAMBDA, C, and K must not be complex.> burrinv (1, 2, 3, i)
