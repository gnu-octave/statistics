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
## @deftypefn {Function File} {[@var{mn}, @var{v}] =} fstat (@var{m}, @var{n})
## Returns mean and variance of the F distribution
##
## Arguments are
##
## @itemize
## @item
## @var{m} is the first parameter of the F distribution. The elements
## of @var{m} must be positive
## @item
## @var{n} is the second parameter of the F distribution. The
## elements of @var{n} must be positive
## @end itemize
## @var{m} and @var{n} must be of common size or one of them must be scalar
##
## Return values are
##
## @itemize
## @item
## @var{mn} is the mean of the F distribution. The mean is undefined for
## @var{n} not greater than 2
## @item
## @var{v} is the variance of the F distribution. The variance is undefined
## for @var{n} not greater than 4
## @end itemize
##
## Examples:
##
## @example
## m = 1:6;
## n = 5:10;
## [mn, v] = fstat (m, n)
##
## [mn, v] = fstat (m, 5)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/F_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [mn, v] = fstat (m, n)

  # Check arguments
  if (nargin != 2)
    usage ("[mn, v] = fstat (m, n)");
  endif

  if (! isempty (m) && ! ismatrix (m))
    error ("fstat: m must be a numeric matrix");
  endif
  if (! isempty (n) && ! ismatrix (n))
    error ("fstat: n must be a numeric matrix");
  endif

  if (! isscalar (m) || ! isscalar (n))
    [retval, m, n] = common_size (m, n);
    if (retval > 0)
      error ("fstat: m and n must be of common size or scalar");
    endif
  endif

  # Calculate moments
  mn = n ./ (n - 2);
  v = (2 .* (n .^ 2) .* (m + n - 2)) ./ (m .* ((n - 2) .^ 2) .* (n - 4));

  # Continue argument check
  k = find (! (m > 0) | ! (m < Inf) | ! (n > 2) | ! (n < Inf));
  if (any (k))
    mn (k) = NaN;
    v (k) = NaN;
  endif

  k = find (! (n > 4));
  if (any (k))
    v (k) = NaN;
  endif

endfunction

%!test
%! m = 1:6;
%! n = 5:10;
%! [mn, v] = fstat (m, n);
%! expected_mn = [1.6667, 1.5000, 1.4000, 1.3333, 1.2857, 1.2500];
%! expected_v = [22.2222, 6.7500, 3.4844, 2.2222, 1.5869, 1.2153];
%! assert (mn, expected_mn, 0.001);
%! assert (v, expected_v, 0.001);

%!test
%! m = 1:6;
%! [mn, v] = fstat (m, 5);
%! expected_mn = [1.6667, 1.6667, 1.6667, 1.6667, 1.6667, 1.6667];
%! expected_v = [22.2222, 13.8889, 11.1111, 9.7222, 8.8889, 8.3333];
%! assert (mn, expected_mn, 0.001);
%! assert (v, expected_v, 0.001);
