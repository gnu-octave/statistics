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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} binostat (@var{n}, @var{p})
## Returns mean and variance of the binomial distribution
##
## Arguments are
##
## @itemize
## @item
## @var{n} is the first parameter of the binomial distribution. The elements
## of @var{n} must be natural numbers
## @item
## @var{p} is the second parameter of the binomial distribution. The
## elements of @var{p} must be probabilities
## @end itemize
## @var{n} and @var{p} must be of common size or one of them must be scalar
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the binomial distribution
## @item
## @var{v} is the variance of the binomial distribution
## @end itemize
##
## Examples:
##
## @example
## n = 1:6;
## p = 0:0.2:1;
## [m, v] = binostat (n, p)
##
## [m, v] = binostat (n, 0.5)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/Binomial_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [m, v] = binostat (n, p)

  # Check arguments
  if (nargin != 2)
    usage ("[m, v] = binostat (n, p)");
  endif

  if (! isempty (n) && ! ismatrix (n))
    error ("binostat: n must be a numeric matrix");
  endif
  if (! isempty (p) && ! ismatrix (p))
    error ("binostat: p must be a numeric matrix");
  endif

  if (! isscalar (n) || ! isscalar (p))
    [retval, n, p] = common_size (n, p);
    if (retval > 0)
      error ("binostat: n and p must be of common size or scalar");
    endif
  endif

  # Calculate moments
  m = n .* p;
  v = n .* p .* (1 - p);

  # Continue argument check
  k = find (! (n > 0) | ! (n < Inf) | ! (n == round (n)) | ! (p >= 0) | ! (p <= 1));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! n = 1:6;
%! p = 0:0.2:1;
%! [m, v] = binostat (n, p);
%! expected_m = [0.00, 0.40, 1.20, 2.40, 4.00, 6.00];
%! expected_v = [0.00, 0.32, 0.72, 0.96, 0.80, 0.00];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);

%!test
%! n = 1:6;
%! [m, v] = binostat (n, 0.5);
%! expected_m = [0.50, 1.00, 1.50, 2.00, 2.50, 3.00];
%! expected_v = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
