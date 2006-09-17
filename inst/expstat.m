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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} expstat (@var{l})
## Returns mean and variance of the exponential distribution
##
## Arguments are
##
## @itemize
## @item
## @var{l} is the rate parameter of the exponential distribution. The
## elements of @var{l} must be positive
## @end itemize
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the exponential distribution
## @item
## @var{v} is the variance of the exponential distribution
## @end itemize
##
## Example:
##
## @example
## l = 1 ./ (1:6);
## [m, v] = expstat (l)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/Exponential_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [m, v] = expstat (l)

  # Check arguments
  if (nargin != 1)
    usage ("[m, v] = expstat (l)");
  endif

  if (! isempty (l) && ! ismatrix (l))
    error ("expstat: l must be a numeric matrix");
  endif

  # Calculate moments
  m = 1 ./ l;
  v = m .^ 2;

  # Continue argument check
  k = find (! (l > 0) | ! (l < Inf));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! l = 1 ./ (1:6);
%! [m, v] = expstat (l);
%! assert (m, [1, 2, 3, 4, 5, 6], 0.001);
%! assert (v, [1, 4, 9, 16, 25, 36], 0.001);
