## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} ricestat (@var{s}, @var{sigma})
##
## Compute statistics of the Rician distribution.
##
## @code{[@var{m}, @var{v}] = ricestat (@var{s}, @var{sigma})} returns the mean
## and variance of the Rician distribution with non-centrality (distance)
## parameter @var{s} and scale parameter @var{sigma}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the Rician distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rice_distribution}
##
## @seealso{ricecdf, riceinv, ricepdf, ricernd, ricefit, ricelike}
## @end deftypefn

function [m, v] = ricestat (s, sigma)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("ricestat: function called with too few input arguments.");
  endif

  ## Check for S and SIGMA being numeric
  if (! (isnumeric (s) && isnumeric (sigma)))
    error ("ricestat: S and SIGMA must be numeric.");
  endif

  ## Check for S and SIGMA being real
  if (iscomplex (s) || iscomplex (sigma))
    error ("ricestat: S and SIGMA must not be complex.");
  endif

  ## Check for common size of S and SIGMA
  if (! isscalar (s) || ! isscalar (sigma))
    [retval, s, sigma] = common_size (s, sigma);
    if (retval > 0)
      error ("ricestat: S and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Initialize mean and variance
  if (isa (s, "single") || isa (sigma, "single"))
    m = NaN (size (s), "single");
    v = m;
  else
    m = NaN (size (s));
    v = m;
  endif

  ## Compute mean and variance for valid parameter values.
  k = (s >= 0 & sigma > 0);
  if (any (k(:)))
    thetasq = (s(k) .^ 2) ./ (sigma(k) .^ 2);
    L = Laguerre_half (-0.5 .* thetasq);
    m(k) = sigma(k) .* sqrt (pi / 2) .* L;
    v(k) = 2 * (sigma(k) .^ 2) + s(k) .^ 2 - ...
           (0.5 .* pi .* sigma(k) .^ 2) .* L .^ 2;
  endif

endfunction

function L = Laguerre_half(x)
  L = exp (x ./ 2) .* ((1 - x) .* besseli (0, -x./2) - x .* besseli (1, -x./2));
endfunction

## Input validation tests
%!error<ricestat: function called with too few input arguments.> ricestat ()
%!error<ricestat: function called with too few input arguments.> ricestat (1)
%!error<ricestat: S and SIGMA must be numeric.> ricestat ({}, 2)
%!error<ricestat: S and SIGMA must be numeric.> ricestat (1, "")
%!error<ricestat: S and SIGMA must not be complex.> ricestat (i, 2)
%!error<ricestat: S and SIGMA must not be complex.> ricestat (1, i)
%!error<ricestat: S and SIGMA must be of common size or scalars.> ...
%! ricestat (ones (3), ones (2))
%!error<ricestat: S and SIGMA must be of common size or scalars.> ...
%! ricestat (ones (2), ones (3))

## Output validation tests
%!shared s, sigma
%! s = [2, 0, -1, 1, 4];
%! sigma = [1, NaN, 3, -1, 2];
%!assert (ricestat (s, sigma), [2.2724, NaN, NaN, NaN, 4.5448], 1e-4);
%!assert (ricestat ([s(1:2), s(4:5)], 1), [2.2724, 1.2533, 1.5486, 4.1272], 1e-4);
%!assert (ricestat ([s(1:2), s(4:5)], 3), [4.1665, 3.7599, 3.8637, 5.2695], 1e-4);
%!assert (ricestat ([s(1:2), s(4:5)], 2), [3.0971, 2.5066, 2.6609, 4.5448], 1e-4);
%!assert (ricestat (2, [sigma(1), sigma(3:5)]), [2.2724, 4.1665, NaN, 3.0971], 1e-4);
%!assert (ricestat (0, [sigma(1), sigma(3:5)]), [1.2533, 3.7599, NaN, 2.5066], 1e-4);
%!assert (ricestat (1, [sigma(1), sigma(3:5)]), [1.5486, 3.8637, NaN, 2.6609], 1e-4);
%!assert (ricestat (4, [sigma(1), sigma(3:5)]), [4.1272, 5.2695, NaN, 4.5448], 1e-4);
