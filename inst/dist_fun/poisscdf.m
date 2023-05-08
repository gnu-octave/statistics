## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{p} =} poisscdf (@var{x}, @var{lambda})
## @deftypefnx {statistics} {@var{p} =} poisscdf (@var{x}, @var{lambda}, @qcode{"upper"})
##
## Poisson cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the Poisson distribution with rate parameter @var{lambda}. The
## size of @var{p} is the common size of @var{x} and @var{lambda}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## @code{@var{p} = poisscdf (@var{x}, @var{x}, @var{lambda}, "upper")} computes
## the upper tail probability of the lognormal distribution.
##
## Further information about the Poisson distribution can be found at
## @url{https://en.wikipedia.org/wiki/Poisson_distribution}
##
## @seealso{poissinv, poisspdf, poissrnd, poissfit, poisslike, poisstat}
## @end deftypefn

function p = poisscdf (x, lambda, uflag)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("poisscdf: function called with too few input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin == 3 && strcmpi (uflag, "upper"))
    uflag = true;
  elseif (nargin == 3  && ! strcmpi (uflag, "upper"))
    error ("poisscdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Check for common size of X and LAMBDA
  if (! isscalar (x) || ! isscalar (lambda))
    [retval, x, lambda] = common_size (x, lambda);
    if (retval > 0)
      error ("poisscdf: X and LAMBDA must be of common size or scalars.");
    endif
  endif

  ## Check for X and LAMBDA being reals
  if (iscomplex (x) || iscomplex (lambda))
    error ("poisscdf: X and LAMBDA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (lambda, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Force NaN for out of range parameters or missing data NaN
  is_nan = isnan (x) | isnan (lambda) | (lambda < 0) ...
                     | (isinf (x) & isinf (lambda));
  p(is_nan) = NaN;

  ## Compute P for X >= 0
  x = floor (x);
  k = x >= 0 & ! is_nan & isfinite (lambda);

  ## Return 1 for positive infinite values of X, unless "upper" is given: p = 0
  k1 = isinf (x) & lambda > 0 & isfinite (lambda);
  if (any (k1))
    if (uflag)
      p(k1) = 0;
    else
      p(k1) = 1;
    endif
  endif

  ## Return 1 when X < 0 and "upper" is given
  k1 = x < 0 & lambda > 0 & isfinite (lambda);
  if (any (k1))
    if (uflag)
      p(k1) = 1;
    endif
  endif

  ## Compute Poisson CDF for remaining cases
  x = x(k);
  lambda = lambda(k);
  if (uflag)
    p(k) = gammainc (lambda, x + 1);
  else
    p(k) = gammainc (lambda, x + 1, "upper");
  endif

endfunction

%!demo
%! ## Plot various CDFs from the Poisson distribution
%! x = 0:20;
%! p1 = poisscdf (x, 1);
%! p2 = poisscdf (x, 4);
%! p3 = poisscdf (x, 10);
%! plot (x, p1, "*b", x, p2, "*g", x, p3, "*r")
%! grid on
%! ylim ([0, 1])
%! legend ({"λ = 1", "λ = 4", "λ = 10"}, "location", "southeast")
%! title ("Poisson CDF")
%! xlabel ("values in x (number of occurences)")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-1 0 1 2 Inf];
%! y = [0, gammainc(1, (x(2:4) +1), "upper"), 1];
%!assert (poisscdf (x, ones (1,5)), y)
%!assert (poisscdf (x, 1), y)
%!assert (poisscdf (x, [1 0 NaN 1 1]), [y(1) 1 NaN y(4:5)])
%!assert (poisscdf ([x(1:2) NaN Inf x(5)], 1), [y(1:2) NaN 1 y(5)])

## Test class of input preserved
%!assert (poisscdf ([x, NaN], 1), [y, NaN])
%!assert (poisscdf (single ([x, NaN]), 1), single ([y, NaN]), eps ("single"))
%!assert (poisscdf ([x, NaN], single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<poisscdf: function called with too few input arguments.> poisscdf ()
%!error<poisscdf: function called with too few input arguments.> poisscdf (1)
%!error<poisscdf: invalid argument for upper tail.> poisscdf (1, 2, 3)
%!error<poisscdf: invalid argument for upper tail.> poisscdf (1, 2, "tail")
%!error<poisscdf: X and LAMBDA must be of common size or scalars.> ...
%! poisscdf (ones (3), ones (2))
%!error<poisscdf: X and LAMBDA must be of common size or scalars.> ...
%! poisscdf (ones (2), ones (3))
%!error<poisscdf: X and LAMBDA must not be complex.> poisscdf (i, 2)
%!error<poisscdf: X and LAMBDA must not be complex.> poisscdf (2, i)
