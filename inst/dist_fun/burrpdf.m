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
## @deftypefn  {statistics} {@var{y} =} burrpdf (@var{x}, @var{lambda}, @var{c}, @var{k})
##
## Burr type XII probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Burr type XII distribution with with scale parameter @var{lambda},
## first shape parameter @var{c}, and second shape parameter @var{k}.  The size
## of @var{y} is the common size of @var{x}, @var{lambda}, @var{c}, and @var{k}.
##  A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## Further information about the Burr distribution can be found at
## @url{https://en.wikipedia.org/wiki/Burr_distribution}
##
## @seealso{burrcdf, burrinv, burrrnd, burrfit, burrlike}
## @end deftypefn

function y = burrpdf (x, lambda, c, k)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("burrpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, LANBDA, C and K
  if (! isscalar (x) || ! isscalar (lambda) || ! isscalar (c) || ! isscalar (k))
    [retval, x, lambda, c, k] = common_size (x, lambda, c, k);
    if (retval > 0)
      error ("burrpdf: X, LAMBDA, C, and K must be of common size or scalars.");
    endif
  endif

  ## Check for X, LANBDA, C and K being reals
  if (iscomplex (x) || iscomplex (lambda) || iscomplex (c) || iscomplex (k))
    error ("burrpdf: X, LAMBDA, C, and K must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (lambda, "single") ...
                        || isa (c, "single") || isa (k, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Force NaNs for out of range parameters
  j = isnan (x) | ! (lambda > 0) | ! (c > 0) | ! (k > 0);
  y(j) = NaN;

  ## Find valid values in parameters and data
  j = (x > 0) & (0 < lambda) & (lambda < Inf) & (0 < c) & (c < Inf) ...
              & (0 < k) & (k < Inf);

  ## Compute Burr PDF
  if (isscalar (lambda) && isscalar (c) && isscalar(k))
    y(j) = (c * k / lambda) .* (x(j) / lambda) .^ (c - 1) ./ ...
           (1 + (x(j) / lambda) .^ c) .^ (k + 1);
  else
    y(j) = (c(j) .* k(j) ./ lambda(j) ) .* x(j).^(c(j) - 1) ./ ...
           (1 + (x(j) ./ lambda(j) ) .^ c(j)) .^ (k(j) + 1);
  endif

endfunction

%!demo
%! ## Plot various PDFs from the Burr type XII distribution
%! x = 0.001:0.001:3;
%! y1 = burrpdf (x, 1, 1, 1);
%! y2 = burrpdf (x, 1, 1, 2);
%! y3 = burrpdf (x, 1, 1, 3);
%! y4 = burrpdf (x, 1, 2, 1);
%! y5 = burrpdf (x, 1, 3, 1);
%! y6 = burrpdf (x, 1, 0.5, 2);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", ...
%!       x, y4, "-c", x, y5, "-m", x, y6, "-k")
%! grid on
%! ylim ([0, 2])
%! legend ({"λ = 1, c = 1, k = 1", "λ = 1, c = 1, k = 2", ...
%!          "λ = 1, c = 1, k = 3", "λ = 1, c = 2, k = 1", ...
%!          "λ = 1, c = 3, k = 1", "λ = 1, c = 0.5, k = 2"}, ...
%!         "location", "northeast")
%! title ("Burr type XII PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
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
%!error<burrpdf: function called with too few input arguments.> burrpdf ()
%!error<burrpdf: function called with too few input arguments.> burrpdf (1)
%!error<burrpdf: function called with too few input arguments.> burrpdf (1, 2)
%!error<burrpdf: function called with too few input arguments.> burrpdf (1, 2, 3)
%!error<burrpdf: function called with too many inputs> ...
%! burrpdf (1, 2, 3, 4, 5)
%!error<burrpdf: X, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrpdf (ones (3), ones (2), ones(2), ones(2))
%!error<burrpdf: X, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrpdf (ones (2), ones (3), ones(2), ones(2))
%!error<burrpdf: X, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrpdf (ones (2), ones (2), ones(3), ones(2))
%!error<burrpdf: X, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrpdf (ones (2), ones (2), ones(2), ones(3))
%!error<burrpdf: X, LAMBDA, C, and K must not be complex.> burrpdf (i, 2, 3, 4)
%!error<burrpdf: X, LAMBDA, C, and K must not be complex.> burrpdf (1, i, 3, 4)
%!error<burrpdf: X, LAMBDA, C, and K must not be complex.> burrpdf (1, 2, i, 4)
%!error<burrpdf: X, LAMBDA, C, and K must not be complex.> burrpdf (1, 2, 3, i)

