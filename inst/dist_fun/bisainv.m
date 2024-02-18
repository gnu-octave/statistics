## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2018 John Donoghue
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
## @deftypefn  {statistics} {@var{x} =} bisainv (@var{p}, @var{beta}, @var{gamma})
##
## Inverse of the Birnbaum-Saunders cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the Birnbaum-Saunders distribution with scale parameter @var{beta} and shape
## parameter @var{gamma}.  The size of @var{x} is the common size of @var{p},
## @var{beta}, and @var{gamma}.  A scalar input functions as a constant matrix
## of the same size as the other inputs.
##
## Further information about the Birnbaum-Saunders distribution can be found at
## @url{https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution}
##
## @seealso{bisainv, bisapdf, bisarnd, bisafit, bisalike, bisastat}
## @end deftypefn

function x = bisainv (p, beta, gamma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("bisainv: function called with too few input arguments.");
  endif

  ## Check for common size of X, BETA, and GAMMA
  if (! isscalar (p) || ! isscalar (beta) || ! isscalar (gamma))
    [retval, p, beta, gamma] = common_size (p, beta, gamma);
    if (retval > 0)
      error (strcat (["bisainv: P, BETA, and GAMMA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, BETA, and GAMMA being reals
  if (iscomplex (p) || iscomplex (beta) || iscomplex (gamma))
    error ("bisainv: P, BETA, and GAMMA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (beta, "single") || isa (gamma, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  ## Force NaNs for out of range parameters
  kn = isnan (p) | (p < 0) | (p > 1) | ! (beta > 0) | ! (beta < Inf) ...
                                     | ! (gamma > 0) | ! (gamma < Inf);
  x(kn) = NaN;

  ## Find valid values in parameters
  kv = (beta > 0) & (beta < Inf) & (gamma > 0) & (gamma < Inf);

  ## Handle edge cases
  k0 = (p == 0) & kv;
  x(k0) = 0;
  k1 = (p == 1) & kv;
  x(k1) = Inf;

  ## Handle all other valid cases
  k = (p > 0) & (p < 1) & kv;

  if (isscalar (beta) && isscalar (gamma))
    z = -sqrt (2) .* erfcinv (2 .* p(k)) .* gamma;
    x(k) = 0.25 .* beta .* (z + sqrt (4 + z .^ 2)) .^ 2;
  else
    z = -sqrt (2) .* erfcinv (2 .* p(k)) .* gamma(k);
    x(k) = 0.25 .* beta(k) .* (z + sqrt (4 + z .^ 2)) .^ 2;
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the Birnbaum-Saunders distribution
%! p = 0.001:0.001:0.999;
%! x1 = bisainv (p, 1, 0.5);
%! x2 = bisainv (p, 1, 1);
%! x3 = bisainv (p, 1, 2);
%! x4 = bisainv (p, 1, 5);
%! x5 = bisainv (p, 1, 10);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c", p, x5, "-m")
%! grid on
%! ylim ([0, 10])
%! legend ({"β = 1, γ = 0.5", "β = 1, γ = 1", "β = 1, γ = 2", ...
%!          "β = 1, γ = 5", "β = 1, γ = 10"}, "location", "northwest")
%! title ("Birnbaum-Saunders iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

%!demo
%! ## Plot various iCDFs from the Birnbaum-Saunders distribution
%! p = 0.001:0.001:0.999;
%! x1 = bisainv (p, 1, 0.3);
%! x2 = bisainv (p, 2, 0.3);
%! x3 = bisainv (p, 1, 0.5);
%! x4 = bisainv (p, 3, 0.5);
%! x5 = bisainv (p, 5, 0.5);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c", p, x5, "-m")
%! grid on
%! ylim ([0, 10])
%! legend ({"β = 1, γ = 0.3", "β = 2, γ = 0.3", "β = 1, γ = 0.5", ...
%!          "β = 3, γ = 0.5", "β = 5, γ = 0.5"}, "location", "northwest")
%! title ("Birnbaum-Saunders iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p, y, f
%! f = @(p,b,c) (b * (c * norminv (p) + sqrt (4 + (c * norminv(p))^2))^2) / 4;
%! p = [-1, 0, 1/4, 1/2, 1, 2];
%! y = [NaN, 0, f(1/4, 1, 1), 1, Inf, NaN];
%!assert (bisainv (p, ones (1,6), ones (1,6)), y)
%!assert (bisainv (p, 1, ones (1,6)), y)
%!assert (bisainv (p, ones (1,6), 1), y)
%!assert (bisainv (p, 1, 1), y)
%!assert (bisainv (p, 1, [1, 1, 1, NaN, 1, 1]), [y(1:3), NaN, y(5:6)])
%!assert (bisainv (p, [1, 1, 1, NaN, 1, 1], 1), [y(1:3), NaN, y(5:6)])
%!assert (bisainv ([p, NaN], 1, 1), [y, NaN])

## Test class of input preserved
%!assert (bisainv (single ([p, NaN]), 1, 1), single ([y, NaN]), eps ("single"))
%!assert (bisainv ([p, NaN], 1, single (1)), single ([y, NaN]), eps ("single"))
%!assert (bisainv ([p, NaN], single (1), 1), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<bisainv: function called with too few input arguments.> bisainv ()
%!error<bisainv: function called with too few input arguments.> bisainv (1)
%!error<bisainv: function called with too few input arguments.> bisainv (1, 2)
%!error<bisainv: function called with too many inputs> bisainv (1, 2, 3, 4)
%!error<bisainv: P, BETA, and GAMMA must be of common size or scalars.> ...
%! bisainv (ones (3), ones (2), ones(2))
%!error<bisainv: P, BETA, and GAMMA must be of common size or scalars.> ...
%! bisainv (ones (2), ones (3), ones(2))
%!error<bisainv: P, BETA, and GAMMA must be of common size or scalars.> ...
%! bisainv (ones (2), ones (2), ones(3))
%!error<bisainv: P, BETA, and GAMMA must not be complex.> bisainv (i, 4, 3)
%!error<bisainv: P, BETA, and GAMMA must not be complex.> bisainv (1, i, 3)
%!error<bisainv: P, BETA, and GAMMA must not be complex.> bisainv (1, 4, i)
