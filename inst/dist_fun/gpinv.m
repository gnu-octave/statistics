## Copyright (C) 1997-2015 Kurt Hornik
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
## @deftypefn  {statistics} {@var{x} =} gpinv (@var{p}, @var{k}, @var{sigma}, @var{mu})
##
## Inverse of the generalized Pareto cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the generalized Pareto distribution with shape parameter @var{k}, scale
## parameter @var{sigma}, and location parameter @var{mu}.  The size of @var{x}
## is the common size of @var{p}, @var{k}, @var{sigma}, and @var{mu}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## When @qcode{@var{k} = 0} and @qcode{@var{mu} = 0}, the Generalized Pareto CDF
## is equivalent to the exponential distribution.  When @qcode{@var{k} > 0} and
## @code{@var{mu} = @var{k} / @var{k}} the Generalized Pareto is equivalent to
## the Pareto distribution.  The mean of the Generalized Pareto is not finite
## when @qcode{@var{k} >= 1} and the variance is not finite when
## @qcode{@var{k} >= 1/2}.  When @qcode{@var{k} >= 0}, the Generalized Pareto
## has positive density for @qcode{@var{x} > @var{mu}}, or, when
## @qcode{@var{mu} < 0}, for
## @qcode{0 <= (@var{x} - @var{mu}) / @var{sigma} <= -1 / @var{k}}.
##
## Further information about the generalized Pareto distribution can be found at
## @url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
##
## @seealso{gpcdf, gppdf, gprnd, gpfit, gplike, gpstat}
## @end deftypefn

function x = gpinv (p, k, sigma, mu)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("gpinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, K, SIGMA, and MU
  [retval, p, k, sigma, mu] = common_size (p, k, sigma, mu);
  if (retval > 0)
    error ("gpinv: P, K, SIGMA, and MU must be of common size or scalars.");
  endif

  ## Check for P, K, SIGMA, and MU being reals
  if (iscomplex (p) || iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gpinv: P, K, SIGMA, and MU must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (mu, "single") ...
      || isa (sigma, "single") || isa (k, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  ## Return NaNs for out of range values of sigma parameter
  kx = isnan (p) | ! (0 <= p) | ! (p <= 1) ...
                 | ! (-Inf < mu) | ! (mu < Inf) ...
                 | ! (sigma > 0) | ! (sigma < Inf) ...
                 | ! (-Inf < k) | ! (k < Inf);
  x(kx) = NaN;

  kx = (0 <= p) & (p <= 1) & (-Inf < mu) & (mu < Inf) ...
                & (sigma > 0) & (sigma < Inf) & (-Inf < k) & (k < Inf);
  if (isscalar (mu) && isscalar (sigma) && isscalar (k))
    if (k == 0)
      x(kx) = -log(1 - p(kx));
      x(kx) = sigma * x(kx) + mu;
    elseif (k > 0)
      x(kx) = (1 - p(kx)).^(-k) - 1;
      x(kx) = (sigma / k) * x(kx) + mu;
    elseif (k < 0)
      x(kx) = (1 - p(kx)).^(-k) - 1;
      x(kx) = (sigma / k) * x(kx)  + mu;
    end
  else
    j = kx & (k == 0);
    if (any (j))
      x(j) = -log (1 - p(j));
      x(j) = sigma(j) .* x(j) + mu(j);
    endif

    j = kx & (k > 0);
    if (any (j))
      x(j) = (1 - p(j)).^(-k(j)) - 1;
      x(j) = (sigma(j) ./ k(j)) .* x(j) + mu(j);
    endif

    j = kx & (k < 0);
    if (any (j))
      x(j) = (1 - p(j)).^(-k(j)) - 1;
      x(j) = (sigma(j) ./ k(j)) .* x(j) + mu(j);
    endif
  endif
endfunction

%!demo
%! ## Plot various iCDFs from the generalized Pareto distribution
%! p = 0.001:0.001:0.999;
%! x1 = gpinv (p, 1, 1, 0);
%! x2 = gpinv (p, 5, 1, 0);
%! x3 = gpinv (p, 20, 1, 0);
%! x4 = gpinv (p, 1, 2, 0);
%! x5 = gpinv (p, 5, 2, 0);
%! x6 = gpinv (p, 20, 2, 0);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", ...
%!       p, x4, "-c", p, x5, "-m", p, x6, "-k")
%! grid on
%! ylim ([0, 5])
%! legend ({"ξ = 1, σ = 1, μ = 0", "ξ = 5, σ = 1, μ = 0", ...
%!          "ξ = 20, σ = 1, μ = 0", "ξ = 1, σ = 2, μ = 0", ...
%!          "ξ = 5, σ = 2, μ = 0", "ξ = 20, σ = 2, μ = 0"}, ...
%!         "location", "southeast")
%! title ("Generalized Pareto iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p, y1, y2, y3
%! p = [-1, 0, 1/2, 1, 2];
%! y1 = [NaN, 0, 0.6931471805599453, Inf, NaN];
%! y2 = [NaN, 0, 1, Inf, NaN];
%! y3 = [NaN, 0, 1/2, 1, NaN];
%!assert (gpinv (p, zeros (1,5), ones (1,5), zeros (1,5)), y1)
%!assert (gpinv (p, 0, 1, zeros (1,5)), y1)
%!assert (gpinv (p, 0, ones (1,5), 0), y1)
%!assert (gpinv (p, zeros (1,5), 1, 0), y1)
%!assert (gpinv (p, 0, 1, 0), y1)
%!assert (gpinv (p, 0, 1, [0, 0, NaN, 0, 0]), [y1(1:2), NaN, y1(4:5)])
%!assert (gpinv (p, 0, [1, 1, NaN, 1, 1], 0), [y1(1:2), NaN, y1(4:5)])
%!assert (gpinv (p, [0, 0, NaN, 0, 0], 1, 0), [y1(1:2), NaN, y1(4:5)])
%!assert (gpinv ([p(1:2), NaN, p(4:5)], 0, 1, 0), [y1(1:2), NaN, y1(4:5)])
%!assert (gpinv (p, ones (1,5), ones (1,5), zeros (1,5)), y2)
%!assert (gpinv (p, 1, 1, zeros (1,5)), y2)
%!assert (gpinv (p, 1, ones (1,5), 0), y2)
%!assert (gpinv (p, ones (1,5), 1, 0), y2)
%!assert (gpinv (p, 1, 1, 0), y2)
%!assert (gpinv (p, 1, 1, [0, 0, NaN, 0, 0]), [y2(1:2), NaN, y2(4:5)])
%!assert (gpinv (p, 1, [1, 1, NaN, 1, 1], 0), [y2(1:2), NaN, y2(4:5)])
%!assert (gpinv (p, [1, 1, NaN, 1, 1], 1, 0), [y2(1:2), NaN, y2(4:5)])
%!assert (gpinv ([p(1:2), NaN, p(4:5)], 1, 1, 0), [y2(1:2), NaN, y2(4:5)])
%!assert (gpinv (p, -ones (1,5), ones (1,5), zeros (1,5)), y3)
%!assert (gpinv (p, -1, 1, zeros (1,5)), y3)
%!assert (gpinv (p, -1, ones (1,5), 0), y3)
%!assert (gpinv (p, -ones (1,5), 1, 0), y3)
%!assert (gpinv (p, -1, 1, 0), y3)
%!assert (gpinv (p, -1, 1, [0, 0, NaN, 0, 0]), [y3(1:2), NaN, y3(4:5)])
%!assert (gpinv (p, -1, [1, 1, NaN, 1, 1], 0), [y3(1:2), NaN, y3(4:5)])
%!assert (gpinv (p, -[1, 1, NaN, 1, 1], 1, 0), [y3(1:2), NaN, y3(4:5)])
%!assert (gpinv ([p(1:2), NaN, p(4:5)], -1, 1, 0), [y3(1:2), NaN, y3(4:5)])

## Test class of input preserved
%!assert (gpinv (single ([p, NaN]), 0, 1, 0), single ([y1, NaN]))
%!assert (gpinv ([p, NaN], 0, 1, single (0)), single ([y1, NaN]))
%!assert (gpinv ([p, NaN], 0, single (1), 0), single ([y1, NaN]))
%!assert (gpinv ([p, NaN], single (0), 1, 0), single ([y1, NaN]))
%!assert (gpinv (single ([p, NaN]), 1, 1, 0), single ([y2, NaN]))
%!assert (gpinv ([p, NaN], 1, 1, single (0)), single ([y2, NaN]))
%!assert (gpinv ([p, NaN], 1, single (1), 0), single ([y2, NaN]))
%!assert (gpinv ([p, NaN], single (1), 1, 0), single ([y2, NaN]))
%!assert (gpinv (single ([p, NaN]), -1, 1, 0), single ([y3, NaN]))
%!assert (gpinv ([p, NaN], -1, 1, single (0)), single ([y3, NaN]))
%!assert (gpinv ([p, NaN], -1, single (1), 0), single ([y3, NaN]))
%!assert (gpinv ([p, NaN], single (-1), 1, 0), single ([y3, NaN]))

## Test input validation
%!error<gpinv: function called with too few input arguments.> gpinv ()
%!error<gpinv: function called with too few input arguments.> gpinv (1)
%!error<gpinv: function called with too few input arguments.> gpinv (1, 2)
%!error<gpinv: function called with too few input arguments.> gpinv (1, 2, 3)
%!error<gpinv: P, K, SIGMA, and MU must be of common size or scalars.> ...
%! gpinv (ones (3), ones (2), ones(2), ones(2))
%!error<gpinv: P, K, SIGMA, and MU must be of common size or scalars.> ...
%! gpinv (ones (2), ones (3), ones(2), ones(2))
%!error<gpinv: P, K, SIGMA, and MU must be of common size or scalars.> ...
%! gpinv (ones (2), ones (2), ones(3), ones(2))
%!error<gpinv: P, K, SIGMA, and MU must be of common size or scalars.> ...
%! gpinv (ones (2), ones (2), ones(2), ones(3))
%!error<gpinv: P, K, SIGMA, and MU must not be complex.> gpinv (i, 2, 3, 4)
%!error<gpinv: P, K, SIGMA, and MU must not be complex.> gpinv (1, i, 3, 4)
%!error<gpinv: P, K, SIGMA, and MU must not be complex.> gpinv (1, 2, i, 4)
%!error<gpinv: P, K, SIGMA, and MU must not be complex.> gpinv (1, 2, 3, i)
