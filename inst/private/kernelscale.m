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
## @deftypefn {statistics} {@var{s} =} kernelscale (@var{kernel})
##
## Standardized-distance scale of a named smoothing kernel, i.e. the standard
## deviation of its canonical form.  The compact kernels are rescaled by this
## factor so that each has unit variance and the bandwidth is comparable across
## kernels (as in MATLAB).  Internal helper shared by @code{ksdensity} and
## @code{mvksdensity}; not intended to be called directly.
##
## @end deftypefn

function s = kernelscale (kernel)
  switch (kernel)
    case 'box'
      s = 1 / sqrt (3);
    case 'triangle'
      s = 1 / sqrt (6);
    case 'epanechnikov'
      s = 1 / sqrt (5);
    otherwise
      s = 1;
  endswitch
endfunction
