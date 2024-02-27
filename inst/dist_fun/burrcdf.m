## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{p} =} burrcdf (@var{x}, @var{lambda}, @var{c}, @var{k})
## @deftypefnx {statistics} {@var{p} =} burrcdf (@var{x}, @var{lambda}, @var{c}, @var{k}, @qcode{"upper"})
##
## Burr type XII cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the Burr type XII distribution with scale parameter @var{lambda},
## first shape parameter @var{c}, and second shape parameter @var{k}.  The size
## of @var{p} is the common size of @var{x}, @var{lambda}, @var{c}, and @var{k}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## @code{@var{p} = burrcdf (@var{x}, @var{beta}, @var{gamma}, "upper")}
## computes the upper tail probability of the Birnbaum-Saunders distribution
## with parameters @var{beta} and @var{gamma}, at the values in @var{x}.
##
## Further information about the Burr distribution can be found at
## @url{https://en.wikipedia.org/wiki/Burr_distribution}
##
## @seealso{burrinv, burrpdf, burrrnd, burrfit, burrlike, burrstat}
## @end deftypefn

function p = burrcdf (x, lambda, c, k, uflag)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("burrcdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 4)
    if (! strcmpi (uflag, "upper"))
      error ("burrcdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, LANBDA, C and K
  if (! isscalar (x) || ! isscalar (lambda) || ! isscalar (c) || ! isscalar (k))
    [retval, x, lambda, c, k] = common_size (x, lambda, c, k);
    if (retval > 0)
      error ("burrcdf: X, LAMBDA, C, and K must be of common size or scalars.");
    endif
  endif

  ## Check for X, LANBDA, C and K being reals
  if (iscomplex (x) || iscomplex (lambda) || iscomplex (c) || iscomplex (k))
    error ("burrcdf: X, LAMBDA, C, and K must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (lambda, "single") || isa (c, "single") ...
                        || isa (k, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Force NaNs for out of range parameters
  j = isnan (x) | ! (lambda > 0) | ! (c > 0) | ! (k > 0);
  p(j) = NaN;

  ## Find valid values in parameters and data
  j = (x > 0) & (lambda > 0) & (lambda < Inf) & (c > 0) & (c < Inf) ...
              & (k > 0) & (k < Inf);

  ## Compute Burr CDF
  if (isscalar (lambda) && isscalar(c) && isscalar(k))
    if (uflag)
      p(j) = (1 + (x(j) / lambda) .^ c) .^ (-k);
    else
      p(j) = 1 - (1 + (x(j) / lambda) .^ c) .^ (-k);
    endif
  else
    if (uflag)
      p(j) = (1 + (x(j) ./ lambda(j)) .^ c(j)) .^ (-k(j));
    else
      p(j) = 1 - (1 + (x(j) ./ lambda(j)) .^ c(j)) .^ (-k(j));
    endif
  endif

endfunction

%!demo
%! ## Plot various CDFs from the Burr type XII distribution
%! x = 0.001:0.001:5;
%! p1 = burrcdf (x, 1, 1, 1);
%! p2 = burrcdf (x, 1, 1, 2);
%! p3 = burrcdf (x, 1, 1, 3);
%! p4 = burrcdf (x, 1, 2, 1);
%! p5 = burrcdf (x, 1, 3, 1);
%! p6 = burrcdf (x, 1, 0.5, 2);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", ...
%!       x, p4, "-c", x, p5, "-m", x, p6, "-k")
%! grid on
%! legend ({"λ = 1, c = 1, k = 1", "λ = 1, c = 1, k = 2", ...
%!          "λ = 1, c = 1, k = 3", "λ = 1, c = 2, k = 1", ...
%!          "λ = 1, c = 3, k = 1", "λ = 1, c = 0.5, k = 2"}, ...
%!         "location", "southeast")
%! title ("Burr type XII CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 1/2, 2/3, 1];
%!assert (burrcdf (x, ones(1,5), ones (1,5), ones (1,5)), y, eps)
%!assert (burrcdf (x, 1, 1, 1), y, eps)
%!assert (burrcdf (x, [1, 1, NaN, 1, 1], 1, 1), [y(1:2), NaN, y(4:5)], eps)
%!assert (burrcdf (x, 1, [1, 1, NaN, 1, 1], 1), [y(1:2), NaN, y(4:5)], eps)
%!assert (burrcdf (x, 1, 1, [1, 1, NaN, 1, 1]), [y(1:2), NaN, y(4:5)], eps)
%!assert (burrcdf ([x, NaN], 1, 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (burrcdf (single ([x, NaN]), 1, 1, 1), single ([y, NaN]), eps("single"))
%!assert (burrcdf ([x, NaN], single (1), 1, 1), single ([y, NaN]), eps("single"))
%!assert (burrcdf ([x, NaN], 1, single (1), 1), single ([y, NaN]), eps("single"))
%!assert (burrcdf ([x, NaN], 1, 1, single (1)), single ([y, NaN]), eps("single"))

## Test input validation
%!error<burrcdf: function called with too few input arguments.> burrcdf ()
%!error<burrcdf: function called with too few input arguments.> burrcdf (1)
%!error<burrcdf: function called with too few input arguments.> burrcdf (1, 2)
%!error<burrcdf: function called with too few input arguments.> burrcdf (1, 2, 3)
%!error<burrcdf: function called with too many inputs> ...
%! burrcdf (1, 2, 3, 4, 5, 6)
%!error<burrcdf: invalid argument for upper tail.> burrcdf (1, 2, 3, 4, "tail")
%!error<burrcdf: invalid argument for upper tail.> burrcdf (1, 2, 3, 4, 5)
%!error<burrcdf: X, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrcdf (ones (3), ones (2), ones(2), ones(2))
%!error<burrcdf: X, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrcdf (ones (2), ones (3), ones(2), ones(2))
%!error<burrcdf: X, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrcdf (ones (2), ones (2), ones(3), ones(2))
%!error<burrcdf: X, LAMBDA, C, and K must be of common size or scalars.> ...
%! burrcdf (ones (2), ones (2), ones(2), ones(3))
%!error<burrcdf: X, LAMBDA, C, and K must not be complex.> burrcdf (i, 2, 3, 4)
%!error<burrcdf: X, LAMBDA, C, and K must not be complex.> burrcdf (1, i, 3, 4)
%!error<burrcdf: X, LAMBDA, C, and K must not be complex.> burrcdf (1, 2, i, 4)
%!error<burrcdf: X, LAMBDA, C, and K must not be complex.> burrcdf (1, 2, 3, i)
