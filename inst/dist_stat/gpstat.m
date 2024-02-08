## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} gpstat (@var{k}, @var{sigma}, @var{mu})
##
## Compute statistics of the generalized Pareto distribution.
##
## @code{[@var{m}, @var{v}] = gpstat (@var{k}, @var{sigma}, @var{mu})}
## returns the mean and variance of the generalized Pareto distribution with
## shape parameter @var{k}, scale parameter @var{sigma}, and location parameter
## @var{mu}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## When @var{k} = 0 and @var{mu} = 0, the generalized Pareto
## distribution is equivalent to the exponential distribution.  When
## @code{@var{k} > 0} and @code{@var{mu} = @var{sigma} / @var{k}},
## the generalized Pareto distribution is equivalent to the Pareto distribution.
## The mean of the generalized Pareto distribution is not finite when
## @code{@var{k} >= 1}, and the variance is not finite when
## @code{@var{k} >= 1/2}.  When @code{@var{k} >= 0}, the generalized
## Pareto distribution has positive density for @code{@var{x} > @var{mu}},
## or, when @code{@var{k} < 0}, for
## @code{0 <= (@var{x} -  @var{mu}) / @var{sigma} <= -1 / @var{k}}.
##
## Further information about the generalized Pareto distribution can be found at
## @url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
##
## @seealso{gpcdf, gpinv, gppdf, gprnd, gpfit, gplike}
## @end deftypefn

function [m, v] = gpstat (k, sigma, mu)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("gpstat: function called with too few input arguments.");
  endif

  ## Check for K, SIGMA, and MU being numeric
  if (! (isnumeric (k) && isnumeric (sigma) && isnumeric (mu)))
    error ("gpstat: K, SIGMA, and MU must be numeric.");
  endif

  ## Check for K, SIGMA, and MU being real
  if (iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gpstat: K, SIGMA, and MU must not be complex.");
  endif

  ## Check for common size of K, SIGMA, and MU
  if (! isscalar (k) || ! isscalar (sigma) || ! isscalar (mu))
    [retval, k, sigma, mu] = common_size (k, sigma, mu);
    if (retval > 0)
      error ("gpstat: K, SIGMA, and MU must be of common size or scalars.");
    endif
  endif

  ## Return NaNs for out of range SCALE parameters.
  sigma(sigma <= 0) = NaN;

  ## Check for appropriate class
  if (isa (k, "single") || isa (sigma, "single") || isa (mu, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Prepare output
  m = NaN (size (k), is_class);
  v = NaN (size (k), is_class);

  ## Compute cases for SHAPE == 0
  knot0 = (abs (k) < eps (is_class));
  m(knot0) = 1;
  v(knot0) = 1;

  ## Compute cases for SHAPE != 0
  knot0 = ! knot0;
  ## SHAPE < 1
  kless = knot0 & (k < 1);
  m(kless) = 1 ./ (1 - k(kless));
  ## SHAPE > 1
  m(k >= 1) = Inf;
  ## SHAPE < 1/2
  ## Find the k~=0 cases and fill in the variance.
  kless = knot0 & (k < 1/2);
  v(kless) = 1 ./ ((1-k(kless)).^2 .* (1-2.*k(kless)));
  ## SHAPE > 1/2
  v(k >= 1/2) = Inf;

  ## Compute mean and variance
  m = mu + sigma .* m;
  v = sigma .^ 2 .* v;

endfunction

## Input validation tests
%!error<gpstat: function called with too few input arguments.> gpstat ()
%!error<gpstat: function called with too few input arguments.> gpstat (1)
%!error<gpstat: function called with too few input arguments.> gpstat (1, 2)
%!error<gpstat: K, SIGMA, and MU must be numeric.> gpstat ({}, 2, 3)
%!error<gpstat: K, SIGMA, and MU must be numeric.> gpstat (1, "", 3)
%!error<gpstat: K, SIGMA, and MU must be numeric.> gpstat (1, 2, "")
%!error<gpstat: K, SIGMA, and MU must not be complex.> gpstat (i, 2, 3)
%!error<gpstat: K, SIGMA, and MU must not be complex.> gpstat (1, i, 3)
%!error<gpstat: K, SIGMA, and MU must not be complex.> gpstat (1, 2, i)
%!error<gpstat: K, SIGMA, and MU must be of common size or scalars.> ...
%! gpstat (ones (3), ones (2), 3)
%!error<gpstat: K, SIGMA, and MU must be of common size or scalars.> ...
%! gpstat (ones (2), 2, ones (3))
%!error<gpstat: K, SIGMA, and MU must be of common size or scalars.> ...
%! gpstat (1, ones (2), ones (3))

## Output validation tests
%!shared x, y
%! x = [-Inf, -1, 0, 1/2, 1, Inf];
%! y = [0, 0.5, 1, 2, Inf, Inf];
%!assert (gpstat (x, ones (1,6), zeros (1,6)), y, eps)

## Test class of input preserved
%!assert (gpstat (single (x), 1, 0), single (y), eps("single"))
%!assert (gpstat (x, single (1), 0), single (y), eps("single"))
%!assert (gpstat (x, 1, single (0)), single (y), eps("single"))
%!assert (gpstat (single ([x, NaN]), 1, 0), single ([y, NaN]), eps("single"))
%!assert (gpstat ([x, NaN], single (1), 0), single ([y, NaN]), eps("single"))
%!assert (gpstat ([x, NaN], 1, single (0)), single ([y, NaN]), eps("single"))
