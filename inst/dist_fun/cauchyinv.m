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
## @deftypefn  {statistics} {@var{x} =} cauchyinv (@var{p}, @var{x0}, @var{gamma})
##
## Inverse of the Cauchy cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the Cauchy distribution with location parameter @var{x0} and scale parameter
## @var{gamma}.  The size of @var{x} is the common size of @var{p}, @var{x0},
## and @var{gamma}.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## Further information about the Cauchy distribution can be found at
## @url{https://en.wikipedia.org/wiki/Cauchy_distribution}
##
## @seealso{cauchycdf, cauchypdf, cauchyrnd}
## @end deftypefn

function x = cauchyinv (p, x0, gamma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("cauchyinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, X0, and GAMMA
  if (! isscalar (p) || ! isscalar (x0) || ! isscalar (gamma))
    [retval, p, x0, gamma] = common_size (p, x0, gamma);
    if (retval > 0)
      error (strcat ("cauchyinv: P, X0, and GAMMA must be of", ...
                     " common size or scalars."));
    endif
  endif

  ## Check for P, X0, and GAMMA being reals
  if (iscomplex (p) || iscomplex (x0) || iscomplex (gamma))
    error ("cauchyinv: P, X0, and GAMMA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (x0, "single") || isa (gamma, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Find valid values in parameters
  ok = ! isinf (x0) & (gamma > 0) & (gamma < Inf);

  ## Handle edge cases
  k0 = (p == 0) & ok;
  x(k0) = -Inf;
  k1 = (p == 1) & ok;
  x(k1) = Inf;

  ## Handle all other valid cases
  k = (p > 0) & (p < 1) & ok;

  if (isscalar (x0) && isscalar (gamma))
    x(k) = x0 - gamma * cot (pi * p(k));
  else
    x(k) = x0(k) - gamma(k) .* cot (pi * p(k));
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the Cauchy distribution
%! p = 0.001:0.001:0.999;
%! x1 = cauchyinv (p, 0, 0.5);
%! x2 = cauchyinv (p, 0, 1);
%! x3 = cauchyinv (p, 0, 2);
%! x4 = cauchyinv (p, -2, 1);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c")
%! grid on
%! ylim ([-5, 5])
%! legend ({"x0 = 0, γ = 0.5", "x0 = 0, γ = 1", ...
%!          "x0 = 0, γ = 2", "x0 = -2, γ = 1"}, "location", "northwest")
%! title ("Cauchy iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (cauchyinv (p, ones (1,5), 2 * ones (1,5)), [NaN -Inf 1 Inf NaN], eps)
%!assert (cauchyinv (p, 1, 2 * ones (1,5)), [NaN -Inf 1 Inf NaN], eps)
%!assert (cauchyinv (p, ones (1,5), 2), [NaN -Inf 1 Inf NaN], eps)
%!assert (cauchyinv (p, [1 -Inf NaN Inf 1], 2), [NaN NaN NaN NaN NaN])
%!assert (cauchyinv (p, 1, 2 * [1 0 NaN Inf 1]), [NaN NaN NaN NaN NaN])
%!assert (cauchyinv ([p(1:2) NaN p(4:5)], 1, 2), [NaN -Inf NaN Inf NaN])
%!assert (cauchyinv ([p, NaN], 1, 2), [NaN -Inf 1 Inf NaN NaN], eps)

## Test class of input preserved
%!assert (cauchyinv (single ([p, NaN]), 1, 2), ...
%! single ([NaN -Inf 1 Inf NaN NaN]), eps ("single"))
%!assert (cauchyinv ([p, NaN], single (1), 2), ...
%! single ([NaN -Inf 1 Inf NaN NaN]), eps ("single"))
%!assert (cauchyinv ([p, NaN], 1, single (2)), ...
%! single ([NaN -Inf 1 Inf NaN NaN]), eps ("single"))

## Test input validation
%!error<cauchyinv: function called with too few input arguments.> cauchyinv ()
%!error<cauchyinv: function called with too few input arguments.> cauchyinv (1)
%!error<cauchyinv: function called with too few input arguments.> ...
%! cauchyinv (1, 2)
%!error<cauchyinv: function called with too many inputs> cauchyinv (1, 2, 3, 4)
%!error<cauchyinv: P, X0, and GAMMA must be of common size or scalars.> ...
%! cauchyinv (ones (3), ones (2), ones(2))
%!error<cauchyinv: P, X0, and GAMMA must be of common size or scalars.> ...
%! cauchyinv (ones (2), ones (3), ones(2))
%!error<cauchyinv: P, X0, and GAMMA must be of common size or scalars.> ...
%! cauchyinv (ones (2), ones (2), ones(3))
%!error<cauchyinv: P, X0, and GAMMA must not be complex.> cauchyinv (i, 4, 3)
%!error<cauchyinv: P, X0, and GAMMA must not be complex.> cauchyinv (1, i, 3)
%!error<cauchyinv: P, X0, and GAMMA must not be complex.> cauchyinv (1, 4, i)
