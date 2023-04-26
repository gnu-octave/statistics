## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2018 John Donoghue
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
## @deftypefn  {statistics} {@var{x} =} bisainv (@var{p}, @var{a}, @var{b}, @var{mu})
##
## Inverse of the Birnbaum-Saunders cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Birnbaum-Saunders distribution with shape parameter
## @var{a}, scale parameter @var{b}, and location parameter @var{mu}.  The size
## of @var{x} is the common size of @var{p}, @var{a}, @var{b}, and @var{mu}.  A
## scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## Further information about the Birnbaum-Saunders distribution can be found at
## @url{https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution}
##
## @seealso{bisacdf, bisapdf, bisarnd, bisafit, bisalike, bisastat}
## @end deftypefn

function x = bisainv (p, a, b, mu)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("bisainv: function called with too few input arguments.");
  endif

  ## Check for common size of X, A, B, and MU
  if (! isscalar (p) || ! isscalar (a) || ! isscalar (b) ...
                     || ! isscalar(mu))
    [retval, p, a, b, mu] = common_size (p, a, b, mu);
    if (retval > 0)
      error (strcat (["bisainv: P, A, B, and MU must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, A, B, and MU being reals
  if (iscomplex (p) || iscomplex (a) || iscomplex (b) ...
                    || iscomplex (mu))
    error ("bisainv: P, A, B, and MU must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (a, "single") || isa (b, "single") ...
                        || isa (mu, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  ## Force NaNs for out of range parameters.
  kn = isnan (p) | (p < 0) | (p > 1) | ! (mu > -Inf) | ! (mu < Inf) ...
                 | ! (b > 0) | ! (b < Inf) | ! (a > 0) | ! (a < Inf);
  x(kn) = NaN;

  ## Find valid values in parameters and data
  kv = (mu > -Inf) & (mu < Inf) & (b > 0) & (b < Inf) & (a > 0) & (a < Inf);

  ## Handle edge cases
  k0 = (p == 0) & kv;
  x(k0) = 0;
  k1 = (p == 1) & kv;
  x(k1) = Inf;

  ## Handle all other valid cases
  k = (p > 0) & (p < 1) & kv;

  if (isscalar (a) && isscalar (b) && isscalar (mu))
    z = -sqrt (2) .* erfcinv (2 .* p(k)) .* a;
    x(k) = mu + 0.25 .* b .* (z + sqrt (4 + z .^ 2)) .^ 2;
  else
    z = -sqrt (2) .* erfcinv (2 .* p(k)) .* a(k);
    x(k) = mu(k) + 0.25 .* b(k) .* (z + sqrt (4 + z .^ 2)) .^ 2;
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the Birnbaum-Saunders distribution
%! p = 0.001:0.001:0.999;
%! x1 = bisainv (p, 0.5, 1, 0);
%! x2 = bisainv (p, 1, 1, 0);
%! x3 = bisainv (p, 2, 1, 0);
%! x4 = bisainv (p, 5, 1, 0);
%! x5 = bisainv (p, 10, 1, 0);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c", p, x5, "-m")
%! grid on
%! ylim ([0, 10])
%! legend ({"α = 0.5, β = 1, μ = 0", "α = 1,    β = 1, μ = 0", ...
%!          "α = 2,    β = 1, μ = 0", "α = 5,    β = 1, μ = 0", ...
%!          "α = 10,  β = 1, μ = 0"}, "location", "northwest")
%! title ("Birnbaum-Saunders iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

%!demo
%! ## Plot various iCDFs from the Birnbaum-Saunders distribution
%! p = 0.001:0.001:0.999;
%! x1 = bisainv (p, 0.3, 1, 0);
%! x2 = bisainv (p, 0.3, 2, 0);
%! x3 = bisainv (p, 0.5, 1, 0);
%! x4 = bisainv (p, 0.5, 3, 0);
%! x5 = bisainv (p, 0.5, 5, 0);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c", p, x5, "-m")
%! grid on
%! ylim ([0, 10])
%! legend ({"α = 0.3, β = 1, μ = 0", "α = 0.3, β = 2, μ = 0", ...
%!          "α = 0.5, β = 1, μ = 0", "α = 0.5, β = 3, μ = 0", ...
%!          "α = 0.5, β = 5, μ = 0"}, "location", "northwest")
%! title ("Birnbaum-Saunders iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test results
%!shared p, y, f
%! f = @(p,a,b,c) (a + b * (c * norminv (p) + sqrt (4 + (c * norminv(p))^2))^2) / 4;
%! p = [-1, 0, 1/4, 1/2, 1, 2];
%! y = [NaN, 0, f(1/4, 0, 1, 1), 1, Inf, NaN];
%!assert (bisainv (p, ones (1,6), ones (1,6), zeros (1,6)), y)
%!assert (bisainv (p, 1, 1, zeros (1,6)), y)
%!assert (bisainv (p, 1, ones (1,6), 0), y)
%!assert (bisainv (p, ones (1,6), 1, 0), y)
%!assert (bisainv (p, 1, 1, 0), y)
%!assert (bisainv (p, 1, 1, [0, 0, 0, NaN, 0, 0]), [y(1:3), NaN, y(5:6)])
%!assert (bisainv (p, 1, [1, 1, 1, NaN, 1, 1], 0), [y(1:3), NaN, y(5:6)])
%!assert (bisainv (p, [1, 1, 1, NaN, 1, 1], 1, 0), [y(1:3), NaN, y(5:6)])
%!assert (bisainv ([p, NaN], 1, 1, 0), [y, NaN])

## Test class of input preserved
%!assert (bisainv (single ([p, NaN]), 1, 1, 0), single ([y, NaN]))
%!assert (bisainv ([p, NaN], 1, 1, single (0)), single ([y, NaN]))
%!assert (bisainv ([p, NaN], 1, single (1), 0), single ([y, NaN]))
%!assert (bisainv ([p, NaN], single (1), 1, 0), single ([y, NaN]))

## Test input validation
%!error<bisainv: function called with too few input arguments.> bisainv ()
%!error<bisainv: function called with too few input arguments.> bisainv (1)
%!error<bisainv: function called with too few input arguments.> bisainv (1, 2)
%!error<bisainv: function called with too few input arguments.> ...
%! bisainv (1, 2, 3)
%!error<bisainv: function called with too many inputs> bisainv (1, 2, 3, 4, 5)
%!error<bisainv: P, A, B, and MU must be of common size or scalars.> ...
%! bisainv (ones (3), ones (2), ones(2), ones(2))
%!error<bisainv: P, A, B, and MU must be of common size or scalars.> ...
%! bisainv (ones (2), ones (3), ones(2), ones(2))
%!error<bisainv: P, A, B, and MU must be of common size or scalars.> ...
%! bisainv (ones (2), ones (2), ones(3), ones(2))
%!error<bisainv: P, A, B, and MU must be of common size or scalars.> ...
%! bisainv (ones (2), ones (2), ones(2), ones(3))
%!error<bisainv: P, A, B, and MU must not be complex.> bisainv (i, 4, 3, 2)
%!error<bisainv: P, A, B, and MU must not be complex.> bisainv (1, i, 3, 2)
%!error<bisainv: P, A, B, and MU must not be complex.> bisainv (1, 4, i, 2)
%!error<bisainv: P, A, B, and MU must not be complex.> bisainv (1, 4, 3, i)

