## Copyright (C) 2006 Arno Onken
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{m}, @var{v}] =} raylstat (@var{sigma})
## Returns mean and variance of the Rayleigh distribution.
##
## Arguments are
##
## @itemize
## @item
## @var{sigma} is the parameter of the Rayleigh distribution. The elements
## of @var{sigma} must be positive.
## @end itemize
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the Rayleigh distribution.
##
## @item
## @var{v} is the variance of the Rayleigh distribution.
## @end itemize
##
## Example:
##
## @example
## sigma = 1:6;
## [m, v] = raylstat (sigma)
## @end example
##
## References:
##
## @enumerate
## @item
## W. L. Martinez and A. R. Martinez. @cite{Computational Statistics
## Handbook with MATLAB.} Chapman & Hall/CRC, pages 547-557, 2001.
##
## @item
## Wikipedia contributors. Rayleigh distribution. @cite{Wikipedia, The Free
## Encyclopedia.}
## @uref{http://en.wikipedia.org/w/index.php?title=Rayleigh_distribution&oldid=69294908},
## August 2006.
## @end enumerate
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>
## Description: Moments of the Rayleigh distribution

function [m, v] = raylstat (sigma)

  # Check arguments
  if (nargin != 1)
    usage ("[m, v] = raylstat (sigma)");
  endif

  if (! isempty (sigma) && ! ismatrix (sigma))
    error ("raylstat: sigma must be a numeric matrix");
  endif

  # Calculate moments
  m = sigma .* sqrt (pi ./ 2);
  v = (2 - pi ./ 2) .* sigma .^ 2;

  # Continue argument check
  k = find (! (sigma > 0));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! sigma = 1:6;
%! [m, v] = raylstat (sigma);
%! expected_m = [1.2533, 2.5066, 3.7599, 5.0133, 6.2666, 7.5199];
%! expected_v = [0.4292, 1.7168, 3.8628, 6.8673, 10.7301, 15.4513];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
