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
## @deftypefn {statistics} {@var{k} =} kernelpdf (@var{z}, @var{kernel})
##
## Unit-variance smoothing kernel density @math{K(z)} evaluated at the
## standardized distances @var{z}.  For a named compact kernel @math{K(z) =
## s * K0(s*z)} with @math{K0} the canonical form and @math{s} its standard
## deviation; @var{kernel} may instead be a function handle, which is applied as
## supplied.  Internal helper shared by @code{ksdensity} and @code{mvksdensity};
## not intended to be called directly.
##
## @end deftypefn

function k = kernelpdf (z, kernel)
  if (is_function_handle (kernel))
    k = kernel (z);
    return;
  endif
  s = kernelscale (kernel);
  z = s * z;
  switch (kernel)
    case 'normal'
      k = exp (-0.5 * z .^ 2) / sqrt (2 * pi);
    case 'box'
      k = 0.5 * (abs (z) <= 1);
    case 'triangle'
      k = max (1 - abs (z), 0);
    case 'epanechnikov'
      k = 0.75 * max (1 - z .^ 2, 0);
  endswitch
  k = s * k;
endfunction
