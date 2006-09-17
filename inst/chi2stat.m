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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} chi2stat (@var{n})
## Returns mean and variance of the chi-square distribution
##
## Arguments are
##
## @itemize
## @item
## @var{n} is the parameter of the chi-square distribution. The elements
## of @var{n} must be positive
## @end itemize
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the chi-square distribution
## @item
## @var{v} is the variance of the chi-square distribution
## @end itemize
##
## Example:
##
## @example
## n = 1:6;
## [m, v] = chi2stat (n)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/Chi-square}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [m, v] = chi2stat (n)

  # Check arguments
  if (nargin != 1)
    usage ("[m, v] = chi2stat (n)");
  endif

  if (! isempty (n) && ! ismatrix (n))
    error ("chi2stat: n must be a numeric matrix");
  endif

  # Calculate moments
  m = n;
  v = 2 .* n;

  # Continue argument check
  k = find (! (n > 0) | ! (n < Inf));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! n = 1:6;
%! [m, v] = chi2stat (n);
%! assert (m, n);
%! assert (v, [2, 4, 6, 8, 10, 12], 0.001);
