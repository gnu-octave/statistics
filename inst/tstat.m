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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} tstat (@var{n})
## Returns mean and variance of the t (Student) distribution
##
## Arguments are
##
## @itemize
## @item
## @var{n} is the parameter of the t (Student) distribution. The elements
## of @var{n} must be positive
## @end itemize
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the t (Student) distribution
## @item
## @var{v} is the variance of the t (Student) distribution
## @end itemize
##
## Example:
##
## @example
## n = 3:8;
## [m, v] = tstat (n)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/T_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [m, v] = tstat (n)

  # Check arguments
  if (nargin != 1)
    usage ("[m, v] = tstat (n)");
  endif

  if (! isempty (n) && ! ismatrix (n))
    error ("tstat: n must be a numeric matrix");
  endif

  # Calculate moments
  m = zeros (size (n));
  v = n ./ (n - 2);

  # Continue argument check
  k = find (! (n > 1) | ! (n < Inf));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif
  k = find (! (n > 2) & (n < Inf));
  if (any (k))
    v (k) = Inf;
  endif

endfunction

%!test
%! n = 3:8;
%! [m, v] = tstat (n);
%! expected_m = [0, 0, 0, 0, 0, 0];
%! expected_v = [3.0000, 2.0000, 1.6667, 1.5000, 1.4000, 1.3333];
%! assert (m, expected_m);
%! assert (v, expected_v, 0.001);
