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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} geostat (@var{p})
## Returns mean and variance of the geometric distribution
##
## Arguments are
##
## @itemize
## @item
## @var{p} is the rate parameter of the geometric distribution. The
## elements of @var{p} must be probabilities
## @end itemize
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the geometric distribution
## @item
## @var{v} is the variance of the geometric distribution
## @end itemize
##
## Example:
##
## @example
## p = 1 ./ (1:6);
## [m, v] = geostat (p)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/Geometric_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [m, v] = geostat (p)

  # Check arguments
  if (nargin != 1)
    usage ("[m, v] = geostat (p)");
  endif

  if (! isempty (p) && ! ismatrix (p))
    error ("geostat: p must be a numeric matrix");
  endif

  # Calculate moments
  q = 1 - p;
  m = q ./ p;
  v = q ./ (p .^ 2);

  # Continue argument check
  k = find (! (p >= 0) | ! (p <= 1));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! p = 1 ./ (1:6);
%! [m, v] = geostat (p);
%! assert (m, [0, 1, 2, 3, 4, 5], 0.001);
%! assert (v, [0, 2, 6, 12, 20, 30], 0.001);
