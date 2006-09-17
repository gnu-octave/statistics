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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} unidstat (@var{n})
## Returns mean and variance of the discrete uniform distribution
##
## Arguments are
##
## @itemize
## @item
## @var{n} is the parameter of the discrete uniform distribution. The elements
## of @var{n} must be positive natural numbers
## @end itemize
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the discrete uniform distribution
## @item
## @var{v} is the variance of the discrete uniform distribution
## @end itemize
##
## Example:
##
## @example
## n = 1:6;
## [m, v] = unidstat (n)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/Uniform_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [m, v] = unidstat (n)

  # Check arguments
  if (nargin != 1)
    usage ("[m, v] = unidstat (n)");
  endif

  if (! isempty (n) && ! ismatrix (n))
    error ("unidstat: n must be a numeric matrix");
  endif

  # Calculate moments
  m = (n + 1) ./ 2;
  v = ((n .^ 2) - 1) ./ 12;

  # Continue argument check
  k = find (! (n > 0) | ! (n < Inf) | ! (n == round (n)));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! n = 1:6;
%! [m, v] = unidstat (n);
%! expected_m = [1.0000, 1.5000, 2.0000, 2.5000, 3.0000, 3.5000];
%! expected_v = [0.0000, 0.2500, 0.6667, 1.2500, 2.0000, 2.9167];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
