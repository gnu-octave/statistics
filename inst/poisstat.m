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
## Returns mean and variance of the Poisson distribution
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{lambda} is the parameter of the Poisson distribution. The
## elements of @var{lambda} must be positive
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{m} is the mean of the Poisson distribution
##
## @item
## @var{v} is the variance of the Poisson distribution
## @end itemize
##
## @subheading Example
##
## @example
## @group
## lambda = 1 ./ (1:6);
## [m, v] = poisstat (lambda)
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Wendy L. Martinez and Angel R. Martinez. @cite{Computational Statistics
## Handbook with MATLAB}. Appendix E, pages 547-557, Chapman & Hall/CRC,
## 2001.
##
## @item
## Athanasios Papoulis. @cite{Probability, Random Variables, and Stochastic
## Processes}. McGraw-Hill, New York, second edition, 1984.
## @end enumerate
## @end deftypefn

## Author: Arno Onken <whyly@whyly.org>
## Description: Moments of the Poisson distribution

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
