## Copyright (C) 2012 Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Copyright (C) 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} gevstat (@var{k}, @var{sigma}, @var{mu})
##
## Compute statistics of the generalized extreme value distribution.
##
## @code{[@var{m}, @var{v}] = gevstat (@var{k}, @var{sigma}, @var{mu})} returns
## the mean and variance of the generalized extreme value distribution with
## shape parameter @var{k}, scale parameter @var{sigma}, and location parameter
## @var{mu}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## The mean of the GEV distribution is not finite when @qcode{@var{k} >= 1}, and
## the variance is not finite when @qcode{@var{k} >= 1/2}.  The GEV distribution
## has positive density only for values of @var{x} such that
## @qcode{@var{k} * (@var{x} - @var{mu}) / @var{sigma} > -1}.
##
## Further information about the generalized extreme value distribution can be
## found at
## @url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}
##
## @seealso{gevcdf, gevinv, gevpdf, gevrnd, gevfit, gevlike}
## @end deftypefn

function [m, v] = gevstat (k, sigma, mu)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("gevstat: function called with too few input arguments.");
  endif

  ## Check for K, SIGMA, and MU being numeric
  if (! (isnumeric (k) && isnumeric (sigma) && isnumeric (mu)))
    error ("gevstat: K, SIGMA, and MU must be numeric.");
  endif

  ## Check for K, SIGMA, and MU being real
  if (iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gevstat: K, SIGMA, and MU must not be complex.");
  endif

  ## Check for common size of K, SIGMA, and MU
  if (! isscalar (k) || ! isscalar (sigma) || ! isscalar (mu))
    [retval, k, sigma, mu] = common_size (k, sigma, mu);
    if (retval > 0)
      error ("gevstat: K, SIGMA, and MU must be of common size or scalars.");
    endif
  endif

  ## Euler-Mascheroni constant
  eg = 0.57721566490153286;

  m = v = k;

  ## Find the mean
  m(k >= 1) = Inf;
  m(k == 0) = mu(k == 0) + eg*sigma(k == 0);
  m(k < 1 & k != 0) = mu(k < 1 & k != 0) + sigma(k < 1 & k != 0) .* ...
                        (gamma(1-k(k < 1 & k != 0)) - 1) ./ k(k < 1 & k != 0);

  ## Find the variance
  v(k >= 0.5) = Inf;
  v(k == 0) = (pi^2 / 6) * sigma(k == 0) .^ 2;
  v(k < 0.5 & k != 0) = (gamma(1-2*k(k < 0.5 & k != 0)) - ...
                         gamma(1-k(k < 0.5 & k != 0)).^2) .* ...
                        (sigma(k < 0.5 & k != 0) ./ k(k < 0.5 & k != 0)) .^ 2;

endfunction

## Input validation tests
%!error<gevstat: function called with too few input arguments.> gevstat ()
%!error<gevstat: function called with too few input arguments.> gevstat (1)
%!error<gevstat: function called with too few input arguments.> gevstat (1, 2)
%!error<gevstat: K, SIGMA, and MU must be numeric.> gevstat ({}, 2, 3)
%!error<gevstat: K, SIGMA, and MU must be numeric.> gevstat (1, "", 3)
%!error<gevstat: K, SIGMA, and MU must be numeric.> gevstat (1, 2, "")
%!error<gevstat: K, SIGMA, and MU must not be complex.> gevstat (i, 2, 3)
%!error<gevstat: K, SIGMA, and MU must not be complex.> gevstat (1, i, 3)
%!error<gevstat: K, SIGMA, and MU must not be complex.> gevstat (1, 2, i)
%!error<gevstat: K, SIGMA, and MU must be of common size or scalars.> ...
%! gevstat (ones (3), ones (2), 3)
%!error<gevstat: K, SIGMA, and MU must be of common size or scalars.> ...
%! gevstat (ones (2), 2, ones (3))
%!error<gevstat: K, SIGMA, and MU must be of common size or scalars.> ...
%! gevstat (1, ones (2), ones (3))

## Output validation tests
%!test
%! k = [-1, -0.5, 0, 0.2, 0.4, 0.5, 1];
%! sigma = 2;
%! mu = 1;
%! [m, v] = gevstat (k, sigma, mu);
%! expected_m = [1, 1.4551, 2.1544, 2.6423, 3.4460, 4.0898, Inf];
%! expected_v = [4, 3.4336, 6.5797, 13.3761, 59.3288, Inf, Inf];
%! assert (m, expected_m, -0.001);
%! assert (v, expected_v, -0.001);
