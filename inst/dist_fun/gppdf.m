## Copyright (C) 2018 John Donoghue
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1997-2015 Kurt Hornik
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
## @deftypefn  {statistics} {@var{y} =} gppdf (@var{x}, @var{k}, @var{sigma}, @var{mu})
##
## Generalized Pareto probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the generalized Pareto distribution with shape parameter @var{k}, scale
## parameter @var{sigma}, and location parameter @var{mu}.  The size of @var{y}
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
## @seealso{gpcdf, gpinv, gprnd, gpfit, gplike, gpstat}
## @end deftypefn

function y = gppdf (x, k, sigma, mu)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("gppdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, K, SIGMA, and MU
  if (! isscalar (x) || ! isscalar (k) || ! isscalar (sigma) || ! isscalar (mu))
    [err, x, k, sigma, mu] = common_size (x, k, sigma, mu);
    if (err > 0)
      error ("gppdf: X, K, SIGMA, and MU must be of common size or scalars.");
    endif
  endif

  ## Check for X, K, SIGMA, and MU being reals
  if (iscomplex (x) || iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gppdf: X, K, SIGMA, and MU must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single") ...
      || isa (k, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Return NaNs for out of range values of sigma parameter
  ky = isnan (x) | ! (-Inf < mu) | ! (mu < Inf) | ...
                   ! (sigma > 0) | ! (sigma < Inf) | ...
                   ! (-Inf < k) | ! (k < Inf);
  y(ky) = NaN;

  ky = (-Inf < x) & (x < Inf) & (-Inf < mu) & (mu < Inf) & ...
        (sigma > 0) & (sigma < Inf) & (-Inf < k) & (k < Inf);
  if (isscalar (mu) && isscalar (sigma) && isscalar (k))
    z = (x - mu) / sigma;

    j = ky & (k == 0) & (z >= 0);
    if (any (j))
      y(j) = exp (-z(j));
    endif

    j = ky & (k > 0) & (z >= 0);
    if (any (j))
      y(j) = (k * z(j) + 1) .^ (-(k + 1) / k) ./ sigma;
    endif

    if (k < 0)
      j = ky & (k < 0) & (0 <= z) & (z <= -1. / k);
      if (any (j))
        y(j) = (k * z(j) + 1) .^ (-(k + 1) / k) ./ sigma;
      endif
    endif
  else
    z = (x - mu) ./ sigma;

    j = ky & (k == 0) & (z >= 0);
    if (any (j))
      y(j) = exp( -z(j));
    endif

    j = ky & (k > 0) & (z >= 0);
    if (any (j))
      y(j) = (k(j) .* z(j) + 1) .^ (-(k(j) + 1) ./ k(j)) ...
                                      ./ sigma(j);
    endif

    if (any (k < 0))
      j = ky & (k < 0) & (0 <= z) & (z <= -1 ./ k);
      if (any (j))
        y(j) = (k(j) .* z(j) + 1) .^ (-(k(j) + 1) ./ k(j)) ...
                                        ./ sigma(j);
      endif
    endif
  endif

endfunction

%!demo
%! ## Plot various PDFs from the generalized Pareto distribution
%! x = 0:0.001:5;
%! y1 = gppdf (x, 1, 1, 0);
%! y2 = gppdf (x, 5, 1, 0);
%! y3 = gppdf (x, 20, 1, 0);
%! y4 = gppdf (x, 1, 2, 0);
%! y5 = gppdf (x, 5, 2, 0);
%! y6 = gppdf (x, 20, 2, 0);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", ...
%!       x, y4, "-c", x, y5, "-m", x, y6, "-k")
%! grid on
%! xlim ([0, 5])
%! ylim ([0, 1])
%! legend ({"ξ = 1, σ = 1, μ = 0", "ξ = 5, σ = 1, μ = 0", ...
%!          "ξ = 20, σ = 1, μ = 0", "ξ = 1, σ = 2, μ = 0", ...
%!          "ξ = 5, σ = 2, μ = 0", "ξ = 20, σ = 2, μ = 0"}, ...
%!         "location", "northeast")
%! title ("Generalized Pareto PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y1, y2, y3
%! x = [-Inf, -1, 0, 1/2, 1, Inf];
%! y1 = [0, 0, 1, 0.6065306597126334, 0.36787944117144233, 0];
%! y2 = [0, 0, 1, 4/9, 1/4, 0];
%! y3 = [0, 0, 1, 1, 1, 0];
%!assert (gppdf (x, zeros (1,6), ones (1,6), zeros (1,6)), y1, eps)
%!assert (gppdf (x, 0, 1, zeros (1,6)), y1, eps)
%!assert (gppdf (x, 0, ones (1,6), 0), y1, eps)
%!assert (gppdf (x, zeros (1,6), 1, 0), y1, eps)
%!assert (gppdf (x, 0, 1, 0), y1, eps)
%!assert (gppdf (x, 0, 1, [0, 0, 0, NaN, 0, 0]), [y1(1:3), NaN, y1(5:6)])
%!assert (gppdf (x, 0, [1, 1, 1, NaN, 1, 1], 0), [y1(1:3), NaN, y1(5:6)])
%!assert (gppdf (x, [0, 0, 0, NaN, 0, 0], 1, 0), [y1(1:3), NaN, y1(5:6)])
%!assert (gppdf ([x(1:3), NaN, x(5:6)], 0, 1, 0), [y1(1:3), NaN, y1(5:6)])
%!assert (gppdf (x, ones (1,6), ones (1,6), zeros (1,6)), y2, eps)
%!assert (gppdf (x, 1, 1, zeros (1,6)), y2, eps)
%!assert (gppdf (x, 1, ones (1,6), 0), y2, eps)
%!assert (gppdf (x, ones (1,6), 1, 0), y2, eps)
%!assert (gppdf (x, 1, 1, 0), y2, eps)
%!assert (gppdf (x, 1, 1, [0, 0, 0, NaN, 0, 0]), [y2(1:3), NaN, y2(5:6)])
%!assert (gppdf (x, 1, [1, 1, 1, NaN, 1, 1], 0), [y2(1:3), NaN, y2(5:6)])
%!assert (gppdf (x, [1, 1, 1, NaN, 1, 1], 1, 0), [y2(1:3), NaN, y2(5:6)])
%!assert (gppdf ([x(1:3), NaN, x(5:6)], 1, 1, 0), [y2(1:3), NaN, y2(5:6)])
%!assert (gppdf (x, -ones (1,6), ones (1,6), zeros (1,6)), y3, eps)
%!assert (gppdf (x, -1, 1, zeros (1,6)), y3, eps)
%!assert (gppdf (x, -1, ones (1,6), 0), y3, eps)
%!assert (gppdf (x, -ones (1,6), 1, 0), y3, eps)
%!assert (gppdf (x, -1, 1, 0), y3, eps)
%!assert (gppdf (x, -1, 1, [0, 0, 0, NaN, 0, 0]), [y3(1:3), NaN, y3(5:6)])
%!assert (gppdf (x, -1, [1, 1, 1, NaN, 1, 1], 0), [y3(1:3), NaN, y3(5:6)])
%!assert (gppdf (x, [-1, -1, -1, NaN, -1, -1], 1, 0), [y3(1:3), NaN, y3(5:6)])
%!assert (gppdf ([x(1:3), NaN, x(5:6)], -1, 1, 0), [y3(1:3), NaN, y3(5:6)])

## Test class of input preserved
%!assert (gppdf (single ([x, NaN]), 0, 1, 0), single ([y1, NaN]))
%!assert (gppdf ([x, NaN], 0, 1, single (0)), single ([y1, NaN]))
%!assert (gppdf ([x, NaN], 0, single (1), 0), single ([y1, NaN]))
%!assert (gppdf ([x, NaN], single (0), 1, 0), single ([y1, NaN]))
%!assert (gppdf (single ([x, NaN]), 1, 1, 0), single ([y2, NaN]))
%!assert (gppdf ([x, NaN], 1, 1, single (0)), single ([y2, NaN]))
%!assert (gppdf ([x, NaN], 1, single (1), 0), single ([y2, NaN]))
%!assert (gppdf ([x, NaN], single (1), 1, 0), single ([y2, NaN]))
%!assert (gppdf (single ([x, NaN]), -1, 1, 0), single ([y3, NaN]))
%!assert (gppdf ([x, NaN], -1, 1, single (0)), single ([y3, NaN]))
%!assert (gppdf ([x, NaN], -1, single (1), 0), single ([y3, NaN]))
%!assert (gppdf ([x, NaN], single (-1), 1, 0), single ([y3, NaN]))

## Test input validation
%!error<gpcdf: function called with too few input arguments.> gpcdf ()
%!error<gpcdf: function called with too few input arguments.> gpcdf (1)
%!error<gpcdf: function called with too few input arguments.> gpcdf (1, 2)
%!error<gpcdf: function called with too few input arguments.> gpcdf (1, 2, 3)
%!error<gpcdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gpcdf (ones (3), ones (2), ones(2), ones(2))
%!error<gpcdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gpcdf (ones (2), ones (3), ones(2), ones(2))
%!error<gpcdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gpcdf (ones (2), ones (2), ones(3), ones(2))
%!error<gpcdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gpcdf (ones (2), ones (2), ones(2), ones(3))
%!error<gpcdf: X, K, SIGMA, and MU must not be complex.> gpcdf (i, 2, 3, 4)
%!error<gpcdf: X, K, SIGMA, and MU must not be complex.> gpcdf (1, i, 3, 4)
%!error<gpcdf: X, K, SIGMA, and MU must not be complex.> gpcdf (1, 2, i, 4)
%!error<gpcdf: X, K, SIGMA, and MU must not be complex.> gpcdf (1, 2, 3, i)
