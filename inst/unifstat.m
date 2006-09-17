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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} unifstat (@var{a}, @var{b})
## Returns mean and variance of the continuous uniform distribution
##
## Arguments are
##
## @itemize
## @item
## @var{a} is the first parameter of the continuous uniform distribution
## @item
## @var{b} is the second parameter of the continuous uniform distribution
## @end itemize
## @var{a} and @var{b} must be of common size or one of them must be scalar
## and @var{a} must be less than @var{b}
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the continuous uniform distribution
## @item
## @var{v} is the variance of the continuous uniform distribution
## @end itemize
##
## Examples:
##
## @example
## a = 1:6;
## b = 2:2:12;
## [m, v] = unifstat (a, b)
##
## [m, v] = unifstat (a, 10)
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

function [m, v] = unifstat (a, b)

  # Check arguments
  if (nargin != 2)
    usage ("[m, v] = unifstat (a, b)");
  endif

  if (! isempty (a) && ! ismatrix (a))
    error ("unifstat: a must be a numeric matrix");
  endif
  if (! isempty (b) && ! ismatrix (b))
    error ("unifstat: b must be a numeric matrix");
  endif

  if (! isscalar (a) || ! isscalar (b))
    [retval, a, b] = common_size (a, b);
    if (retval > 0)
      error ("unifstat: a and b must be of common size or scalar");
    endif
  endif

  # Calculate moments
  m = (a + b) ./ 2;
  v = ((b - a) .^ 2) ./ 12;

  # Continue argument check
  k = find (! (-Inf < a) | ! (a < b) | ! (b < Inf));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! a = 1:6;
%! b = 2:2:12;
%! [m, v] = unifstat (a, b);
%! expected_m = [1.5000, 3.0000, 4.5000, 6.0000, 7.5000, 9.0000];
%! expected_v = [0.0833, 0.3333, 0.7500, 1.3333, 2.0833, 3.0000];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);

%!test
%! a = 1:6;
%! [m, v] = unifstat (a, 10);
%! expected_m = [5.5000, 6.0000, 6.5000, 7.0000, 7.5000, 8.0000];
%! expected_v = [6.7500, 5.3333, 4.0833, 3.0000, 2.0833, 1.3333];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
