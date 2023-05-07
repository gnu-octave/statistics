## Copyright (C) 1997-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2018 John Donoghue
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{p} =} gpcdf (@var{x}, @var{k}, @var{sigma}, @var{mu})
## @deftypefnx {statistics} {@var{p} =} gpcdf (@var{x}, @var{k}, @var{sigma}, @var{mu}, @qcode{"upper"})
##
## Generalized Pareto cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the generalized Pareto distribution with shape parameter @var{k},
## scale parameter @var{sigma}, and location parameter @var{mu}.  The size of
## @var{p} is the common size of @var{x}, @var{k}, @var{sigma}, and @var{mu}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## @code{[@dots{}] = gpcdf(@dots{}, "upper")} computes the upper tail
## probability of the generalized Pareto distribution.
##
## When @qcode{@var{k} = 0} and @qcode{@var{mu} = 0}, the Generalized Pareto CDF
## is equivalent to the exponential distribution.  When @qcode{@var{k} > 0} and
## @code{@var{mu} = @var{k} / @var{k}} the Generalized Pareto is equivalent to
## the Pareto distribution.  The mean of the Generalized Pareto is not finite
## when @qcode{@var{k} >= 1} and the variance is not finite when
## @qcode{@var{k} >= 1/2}.  When @qcode{@var{k} >= 0}, the Generalized Pareto
## has positive density for @qcode{@var{x} > @var{mu}}, or, when
## @qcode{@var{mu} < 0},for
## @qcode{0 <= (@var{x} - @var{mu}) / @var{sigma} <= -1 / @var{k}}.
##
## Further information about the generalized Pareto distribution can be found at
## @url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
##
## @seealso{gpinv, gppdf, gprnd, gpfit, gplike, gpstat}
## @end deftypefn

function p = gpcdf (x, k, sigma, mu, uflag)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("gpcdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 4)
    if (! strcmpi (uflag, "upper"))
      error ("gpcdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, K, SIGMA, and MU
  if (! isscalar (x) || ! isscalar (k) || ! isscalar (sigma) || ! isscalar (mu))
    [err, x, k, sigma, mu] = common_size (x, k, sigma, mu);
    if (err > 0)
      error ("gpcdf: X, K, SIGMA, and MU must be of common size or scalars.");
    endif
  endif

  ## Check for X, K, SIGMA, and MU being reals
  if (iscomplex (x) || iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gpcdf: X, K, SIGMA, and MU must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (k, "single") ...
                        || isa (sigma, "single") || isa (mu, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Prepare output
  p = zeros (size (x), is_class);

  ## Return NaNs for out of range values of sigma parameter
  sigma(sigma <= 0) = NaN;

  ## Calculate (x-mu)/sigma => 0 and force zero below that
  z = (x - mu) ./ sigma;
  z(z < 0) = 0;

  ## Compute cases for SHAPE == 0
  kz = (abs (k) < eps (is_class));
  if (uflag)
    p(kz) = exp (-z(kz));
  else
    p(kz) = -expm1 (-z(kz));
  endif

  ## For SHAPE < 0, calculate 0 <= x/sigma <= -1/k and force zero below that
  t = z .* k;
  kt = (t <= -1 & k < -eps (is_class));
  t(kt) = 0;

  ## Compute cases for SHAPE != 0
  kz = ! kz;
  if (uflag)
    p(kz) = exp ((-1 ./ k(kz)) .* log1p (t(kz)));
  else
    p(kz) = -expm1 ((-1 ./ k(kz)) .* log1p (t(kz)));
  endif
  if (uflag)
    p(kt) = 0;
  else
    p(kt) = 1;
  endif

  ## For SHAPE == NaN force p = NaN
  p(isnan (k)) = NaN;

endfunction

%!demo
%! ## Plot various CDFs from the generalized Pareto distribution
%! x = 0:0.001:5;
%! p1 = gpcdf (x, 1, 1, 0);
%! p2 = gpcdf (x, 5, 1, 0);
%! p3 = gpcdf (x, 20, 1, 0);
%! p4 = gpcdf (x, 1, 2, 0);
%! p5 = gpcdf (x, 5, 2, 0);
%! p6 = gpcdf (x, 20, 2, 0);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", ...
%!       x, p4, "-c", x, p5, "-m", x, p6, "-k")
%! grid on
%! xlim ([0, 5])
%! legend ({"ξ = 1, σ = 1, μ = 0", "ξ = 5, σ = 1, μ = 0", ...
%!          "ξ = 20, σ = 1, μ = 0", "ξ = 1, σ = 2, μ = 0", ...
%!          "ξ = 5, σ = 2, μ = 0", "ξ = 20, σ = 2, μ = 0"}, ...
%!         "location", "northwest")
%! title ("Generalized Pareto CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y1, y1u, y2, y2u, y3, y3u
%! x = [-Inf, -1, 0, 1/2, 1, Inf];
%! y1 = [0, 0, 0, 0.3934693402873666, 0.6321205588285577, 1];
%! y1u = [1, 1, 1, 0.6065306597126334, 0.3678794411714423, 0];
%! y2 = [0, 0, 0, 1/3, 1/2, 1];
%! y2u = [1, 1, 1, 2/3, 1/2, 0];
%! y3 = [0, 0, 0, 1/2, 1, 1];
%! y3u = [1, 1, 1, 1/2, 0, 0];
%!assert (gpcdf (x, zeros (1,6), ones (1,6), zeros (1,6)), y1, eps)
%!assert (gpcdf (x, 0, 1, zeros (1,6)), y1, eps)
%!assert (gpcdf (x, 0, ones (1,6), 0), y1, eps)
%!assert (gpcdf (x, zeros (1,6), 1, 0), y1, eps)
%!assert (gpcdf (x, 0, 1, 0), y1, eps)
%!assert (gpcdf (x, 0, 1, [0, 0, 0, NaN, 0, 0]), [y1(1:3), NaN, y1(5:6)], eps)
%!assert (gpcdf (x, 0, [1, 1, 1, NaN, 1, 1], 0), [y1(1:3), NaN, y1(5:6)], eps)
%!assert (gpcdf (x, [0, 0, 0, NaN, 0, 0], 1, 0), [y1(1:3), NaN, y1(5:6)], eps)
%!assert (gpcdf ([x(1:3), NaN, x(5:6)], 0, 1, 0), [y1(1:3), NaN, y1(5:6)], eps)
%!assert (gpcdf (x, zeros (1,6), ones (1,6), zeros (1,6), "upper"), y1u, eps)
%!assert (gpcdf (x, 0, 1, zeros (1,6), "upper"), y1u, eps)
%!assert (gpcdf (x, 0, ones (1,6), 0, "upper"), y1u, eps)
%!assert (gpcdf (x, zeros (1,6), 1, 0, "upper"), y1u, eps)
%!assert (gpcdf (x, 0, 1, 0, "upper"), y1u, eps)
%!assert (gpcdf (x, ones (1,6), ones (1,6), zeros (1,6)), y2, eps)
%!assert (gpcdf (x, 1, 1, zeros (1,6)), y2, eps)
%!assert (gpcdf (x, 1, ones (1,6), 0), y2, eps)
%!assert (gpcdf (x, ones (1,6), 1, 0), y2, eps)
%!assert (gpcdf (x, 1, 1, 0), y2, eps)
%!assert (gpcdf (x, 1, 1, [0, 0, 0, NaN, 0, 0]), [y2(1:3), NaN, y2(5:6)], eps)
%!assert (gpcdf (x, 1, [1, 1, 1, NaN, 1, 1], 0), [y2(1:3), NaN, y2(5:6)], eps)
%!assert (gpcdf (x, [1, 1, 1, NaN, 1, 1], 1, 0), [y2(1:3), NaN, y2(5:6)], eps)
%!assert (gpcdf ([x(1:3), NaN, x(5:6)], 1, 1, 0), [y2(1:3), NaN, y2(5:6)], eps)
%!assert (gpcdf (x, ones (1,6), ones (1,6), zeros (1,6), "upper"), y2u, eps)
%!assert (gpcdf (x, 1, 1, zeros (1,6), "upper"), y2u, eps)
%!assert (gpcdf (x, 1, ones (1,6), 0, "upper"), y2u, eps)
%!assert (gpcdf (x, ones (1,6), 1, 0, "upper"), y2u, eps)
%!assert (gpcdf (x, 1, 1, 0, "upper"), y2u, eps)
%!assert (gpcdf (x, 1, 1, [0, 0, 0, NaN, 0, 0], "upper"), ...
%!                        [y2u(1:3), NaN, y2u(5:6)], eps)
%!assert (gpcdf (x, 1, [1, 1, 1, NaN, 1, 1], 0, "upper"), ...
%!                        [y2u(1:3), NaN, y2u(5:6)], eps)
%!assert (gpcdf (x, [1, 1, 1, NaN, 1, 1], 1, 0, "upper"), ...
%!                        [y2u(1:3), NaN, y2u(5:6)], eps)
%!assert (gpcdf ([x(1:3), NaN, x(5:6)], 1, 1, 0, "upper"), ...
%!                        [y2u(1:3), NaN, y2u(5:6)], eps)
%!assert (gpcdf (x, -ones (1,6), ones (1,6), zeros (1,6)), y3, eps)
%!assert (gpcdf (x, -1, 1, zeros (1,6)), y3, eps)
%!assert (gpcdf (x, -1, ones (1,6), 0), y3, eps)
%!assert (gpcdf (x, -ones (1,6), 1, 0), y3, eps)
%!assert (gpcdf (x, -1, 1, 0), y3, eps)
%!assert (gpcdf (x, -1, 1, [0, 0, 0, NaN, 0, 0]), [y3(1:3), NaN, y3(5:6)], eps)
%!assert (gpcdf (x, -1, [1, 1, 1, NaN, 1, 1], 0), [y3(1:3), NaN, y3(5:6)], eps)
%!assert (gpcdf (x, [-1, -1, -1, NaN, -1, -1], 1, 0), [y3(1:3), NaN, y3(5:6)], eps)
%!assert (gpcdf ([x(1:3), NaN, x(5:6)], -1, 1, 0), [y3(1:3), NaN, y3(5:6)], eps)
%!assert (gpcdf (x, -ones (1,6), ones (1,6), zeros (1,6), "upper"), y3u, eps)
%!assert (gpcdf (x, -1, 1, zeros (1,6), "upper"), y3u, eps)
%!assert (gpcdf (x, -1, ones (1,6), 0, "upper"), y3u, eps)
%!assert (gpcdf (x, -ones (1,6), 1, 0, "upper"), y3u, eps)
%!assert (gpcdf (x, -1, 1, 0, "upper"), y3u, eps)
%!assert (gpcdf (x, -1, 1, [0, 0, 0, NaN, 0, 0], "upper"), ...
%!                         [y3u(1:3), NaN, y3u(5:6)], eps)
%!assert (gpcdf (x, -1, [1, 1, 1, NaN, 1, 1], 0, "upper"), ...
%!                          [y3u(1:3), NaN, y3u(5:6)], eps)
%!assert (gpcdf (x, [-1, -1, -1, NaN, -1, -1], 1, 0, "upper"), ...
%!                          [y3u(1:3), NaN, y3u(5:6)], eps)
%!assert (gpcdf ([x(1:3), NaN, x(5:6)], -1, 1, 0, "upper"), ...
%!                          [y3u(1:3), NaN, y3u(5:6)], eps)

## Test class of input preserved
%!assert (gpcdf (single ([x, NaN]), 0, 1, 0), single ([y1, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], 0, 1, single (0)), single ([y1, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], 0, single (1), 0), single ([y1, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], single (0), 1, 0), single ([y1, NaN]), eps("single"))
%!assert (gpcdf (single ([x, NaN]), 1, 1, 0), single ([y2, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], 1, 1, single (0)), single ([y2, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], 1, single (1), 0), single ([y2, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], single (1), 1, 0), single ([y2, NaN]), eps("single"))
%!assert (gpcdf (single ([x, NaN]), -1, 1, 0), single ([y3, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], -1, 1, single (0)), single ([y3, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], -1, single (1), 0), single ([y3, NaN]), eps("single"))
%!assert (gpcdf ([x, NaN], single (-1), 1, 0), single ([y3, NaN]), eps("single"))

## Test input validation
%!error<gpcdf: function called with too few input arguments.> gpcdf ()
%!error<gpcdf: function called with too few input arguments.> gpcdf (1)
%!error<gpcdf: function called with too few input arguments.> gpcdf (1, 2)
%!error<gpcdf: function called with too few input arguments.> gpcdf (1, 2, 3)
%!error<gpcdf: invalid argument for upper tail.> gpcdf (1, 2, 3, 4, "tail")
%!error<gpcdf: invalid argument for upper tail.> gpcdf (1, 2, 3, 4, 5)
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
