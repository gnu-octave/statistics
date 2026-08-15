## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {statistics} {@var{c} =} kernelcdf (@var{z}, @var{kernel})
##
## Cumulative distribution function of the unit-variance smoothing kernel,
## evaluated at the standardized distances @var{z}.  For a named kernel this is
## the canonical cdf evaluated at @math{s*z}; @var{kernel} may instead be a
## function handle, whose cumulative is computed numerically over a fine grid.
## Internal helper shared by @code{ksdensity} and @code{mvksdensity}; not
## intended to be called directly.
##
## @end deftypefn

function c = kernelcdf (z, kernel)
  if (is_function_handle (kernel))
    ## Numeric cumulative for a custom kernel over a fine grid.
    g = linspace (-1e3, 1e3, 200001)';
    C = cumtrapz (g, kernel (g));
    C = C / C(end);
    c = reshape (interp1 (g, C, z(:), 'linear', 'extrap'), size (z));
    return;
  endif
  z = kernelscale (kernel) * z;
  switch (kernel)
    case 'normal'
      c = 0.5 * erfc (-z / sqrt (2));
    case 'box'
      c = min (max ((z + 1) / 2, 0), 1);
    case 'triangle'
      zc = max (min (z, 1), -1);
      c = (zc < 0) .* (0.5 * (zc + 1) .^ 2) ...
          + (zc >= 0) .* (0.5 + zc - 0.5 * zc .^ 2);
    case 'epanechnikov'
      zc = max (min (z, 1), -1);
      c = 0.75 * zc - 0.25 * zc .^ 3 + 0.5;
  endswitch
endfunction
