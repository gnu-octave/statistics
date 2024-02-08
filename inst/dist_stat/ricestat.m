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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} ricestat (@var{nu}, @var{sigma})
##
## Compute statistics for the Rician distribution.
##
## @code{[@var{m}, @var{v}] = ricestat (@var{nu}, @var{sigma})} returns the mean
## and variance of the Rician distribution with non-centrality (distance)
## parameter @var{nu} and scale parameter @var{sigma}.
##
## The size of @var{m} and @var{v} is the common size of the input arguments.
## Scalar input arguments @var{nu} and @var{sigma} are regarded as constant
## matrices of the same size as the other input.
##
## Further information about the Rician distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rice_distribution}
##
## @seealso{ricecdf, riceinv, ricepdf, ricernd, ricefit, ricelike}
## @end deftypefn

function [m, v] = ricestat (nu, sigma)

  ## Check for valid input arguments
  if (nargin <  2)
    error ("ricestat: too few input arguments.");
  endif

  ## Check and fix size of input arguments
  [err, nu, sigma] = common_size (nu, sigma);
  if (err > 0)
    error ("ricestat: input size mismatch.");
  endif

  ## Initialize mean and variance
  if (isa (nu, "single") || isa (sigma, "single"))
    m = NaN (size (nu), "single");
    v = m;
  else
    m = NaN (size (nu));
    v = m;
  endif

  ## Compute mean and variance for valid parameter values.
  k = (nu >= 0 & sigma > 0);
  if (any (k(:)))
    thetasq = (nu(k) .^ 2) ./ (sigma(k) .^ 2);
    L = Laguerre_half (-0.5 .* thetasq);
    m(k) = sigma(k) .* sqrt (pi / 2) .* L;
    v(k) = 2 * (sigma(k) .^ 2) + nu(k) .^ 2 - (0.5 .* pi .* sigma(k) .^ 2) .* L;
  endif

endfunction

function L = Laguerre_half(x)
  L = exp (x ./ 2) .* ((1 - x) .* besseli (0, -x./2) - x .* besseli (1, -x./2));
endfunction

## Input validation tests
%!error<ricestat: too few input arguments.> p = ricestat ();
%!error<ricestat: too few input arguments.> p = ricestat (1);
%!error<ricestat: input size mismatch.> p = ricestat ([4, 3], [3, 4, 5]);

## Output validation tests
%!shared nu, sigma
%! nu = [2, 0, -1, 1, 4];
%! sigma = [1, NaN, 3, -1, 2];
%!assert (ricestat (nu, sigma), [2.2724, NaN, NaN, NaN, 4.5448], 1e-4);
%!assert (ricestat ([nu(1:2), nu(4:5)], 1), [2.2724, NaN, 1.5486, 4.1272], 1e-4);
%!assert (ricestat ([nu(1:2), nu(4:5)], 3), [5, NaN, 4, 7]);
%!assert (ricestat ([nu(1:2), nu(4:5)], 2), [4, NaN, 3, 6]);
%!assert (ricestat (2, [sigma(1), sigma(3:5)]), [3, 5, NaN, 4]);
%!assert (ricestat (0, [sigma(1), sigma(3:5)]), [NaN, NaN, NaN, NaN]);
%!assert (ricestat (1, [sigma(1), sigma(3:5)]), [2, 4, NaN, 3]);
%!assert (ricestat (4, [sigma(1), sigma(3:5)]), [5, 7, NaN, 6]);
