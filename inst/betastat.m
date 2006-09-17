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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} betastat (@var{a}, @var{b})
## Returns mean and variance of the beta distribution
##
## Arguments are
##
## @itemize
## @item
## @var{a} is the first parameter of the beta distribution. @var{a} must be
## positive
## @item
## @var{b} is the second parameter of the beta distribution. @var{b} must be
## positive
## @end itemize
## @var{a} and @var{b} must be of common size or one of them must be scalar
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the beta distribution
## @item
## @var{v} is the variance of the beta distribution
## @end itemize
##
## Examples:
##
## @example
## a = 1:6;
## b = 1:0.2:2;
## [m, v] = betastat (a, b)
##
## [m, v] = betastat (a, 1.5)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/Beta_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [m, v] = betastat (a, b)

  # Check arguments
  if (nargin != 2)
    usage ("[m, v] = betastat (a, b)");
  endif

  if (! isempty (a) && ! ismatrix (a))
    error ("betastat: a must be a numeric matrix");
  endif
  if (! isempty (b) && ! ismatrix (b))
    error ("betastat: b must be a numeric matrix");
  endif

  if (! isscalar (a) || ! isscalar (b))
    [retval, a, b] = common_size (a, b);
    if (retval > 0)
      error ("betastat: a and b must be of common size or scalar");
    endif
  endif

  # Calculate moments
  m = a ./ (a + b);
  v = (a .* b) ./ (((a + b) .^ 2) .* (a + b + 1));

  # Continue argument check
  k = find (! (a > 0) | ! (a < Inf) | ! (b > 0) | ! (b < Inf));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! a = 1:6;
%! b = 1:0.2:2;
%! [m, v] = betastat (a, b);
%! expected_m = [0.5000, 0.6250, 0.6818, 0.7143, 0.7353, 0.7500];
%! expected_v = [0.0833, 0.0558, 0.0402, 0.0309, 0.0250, 0.0208];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);

%!test
%! a = 1:6;
%! [m, v] = betastat (a, 1.5);
%! expected_m = [0.4000, 0.5714, 0.6667, 0.7273, 0.7692, 0.8000];
%! expected_v = [0.0686, 0.0544, 0.0404, 0.0305, 0.0237, 0.0188];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
