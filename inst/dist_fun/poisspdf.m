## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn {statistics} {@var{y} =} poisspdf (@var{x}, @var{lambda})
##
## Poisson probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Poisson distribution with rate parameter @var{lambda}.  The size of
## @var{y} is the common size of @var{x} and @var{lambda}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## Further information about the Poisson distribution can be found at
## @url{https://en.wikipedia.org/wiki/Poisson_distribution}
##
## @seealso{poisscdf, poissinv, poissrnd, poissfit, poisslike, poisstat}
## @end deftypefn

function y = poisspdf (x, lambda)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("poisspdf: function called with too few input arguments.");
  endif

  ## Check for common size of X and LAMBDA
  if (! isscalar (x) || ! isscalar (lambda))
    [retval, x, lambda] = common_size (x, lambda);
    if (retval > 0)
      error ("poisspdf: X and LAMBDA must be of common size or scalars.");
    endif
  endif

  ## Check for X and LAMBDA being reals
  if (iscomplex (x) || iscomplex (lambda))
    error ("poisspdf: X and LAMBDA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (lambda, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Force NaN for out of range parameters or missing data NaN
  k = isnan (x) | ! (lambda > 0);
  y(k) = NaN;

  k = (x >= 0) & (x < Inf) & (x == fix (x)) & (lambda > 0);
  if (isscalar (lambda))
    y(k) = exp (x(k) * log (lambda) - lambda - gammaln (x(k) + 1));
  else
    y(k) = exp (x(k) .* log (lambda(k)) - lambda(k) - gammaln (x(k) + 1));
  endif

endfunction

%!demo
%! ## Plot various PDFs from the Poisson distribution
%! x = 0:20;
%! y1 = poisspdf (x, 1);
%! y2 = poisspdf (x, 4);
%! y3 = poisspdf (x, 10);
%! plot (x, y1, "*b", x, y2, "*g", x, y3, "*r")
%! grid on
%! ylim ([0, 0.4])
%! legend ({"λ = 1", "λ = 4", "λ = 10"}, "location", "northeast")
%! title ("Poisson PDF")
%! xlabel ("values in x (number of occurences)")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 1 2 Inf];
%! y = [0, exp(-1)*[1 1 0.5], 0];
%!assert (poisspdf (x, ones (1,5)), y, eps)
%!assert (poisspdf (x, 1), y, eps)
%!assert (poisspdf (x, [1 0 NaN 1 1]), [y(1) NaN NaN y(4:5)], eps)
%!assert (poisspdf ([x, NaN], 1), [y, NaN], eps)

## Test class of input preserved
%!assert (poisspdf (single ([x, NaN]), 1), single ([y, NaN]), eps ("single"))
%!assert (poisspdf ([x, NaN], single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<poisspdf: function called with too few input arguments.> poisspdf ()
%!error<poisspdf: function called with too few input arguments.> poisspdf (1)
%!error<poisspdf: X and LAMBDA must be of common size or scalars.> ...
%! poisspdf (ones (3), ones (2))
%!error<poisspdf: X and LAMBDA must be of common size or scalars.> ...
%! poisspdf (ones (2), ones (3))
%!error<poisspdf: X and LAMBDA must not be complex.> poisspdf (i, 2)
%!error<poisspdf: X and LAMBDA must not be complex.> poisspdf (2, i)
