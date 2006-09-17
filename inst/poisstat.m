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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} poisstat (@var{lambda})
## Returns mean and variance of the poisson distribution
##
## Arguments are
##
## @itemize
## @item
## @var{lambda} is the parameter of the poisson distribution. The
## elements of @var{lambda} must be positive
## @end itemize
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the poisson distribution
## @item
## @var{v} is the variance of the poisson distribution
## @end itemize
##
## Example:
##
## @example
## lambda = 1 ./ (1:6);
## [m, v] = poisstat (lambda)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/Poisson_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [m, v] = poisstat (lambda)

  # Check arguments
  if (nargin != 1)
    usage ("[m, v] = poisstat (lambda)");
  endif

  if (! isempty (lambda) && ! ismatrix (lambda))
    error ("poisstat: lambda must be a numeric matrix");
  endif

  # Set moments
  m = lambda;
  v = lambda;

  # Continue argument check
  k = find (! (lambda > 0) | ! (lambda < Inf));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! lambda = 1 ./ (1:6);
%! [m, v] = poisstat (lambda);
%! assert (m, lambda);
%! assert (v, lambda);
