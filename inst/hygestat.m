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
## @deftypefn {Function File} {[@var{mn}, @var{v}] =} hygestat (@var{m}, @var{t}, @var{n})
## Returns mean and variance of the hypergeometric distribution
##
## Arguments are
##
## @itemize
## @item
## @var{m} is the number of marked items of the hypergeometric distribution.
## The elements of @var{n} must be natural numbers
## @item
## @var{t} is the total size of the population of the hypergeometric
## distribution. The elements of @var{p} must be positive natural numbers
## @item
## @var{n} is the size of the drawn sample of the hypergeometric
## distribution. The elements of @var{p} must be positive natural numbers
## @end itemize
## @var{m}, @var{t} and @var{n} must be of common size or scalar
##
## Return values are
##
## @itemize
## @item
## @var{mn} is the mean of the hypergeometric distribution
## @item
## @var{v} is the variance of the hypergeometric distribution
## @end itemize
##
## Examples:
##
## @example
## m = 0:5;
## t = 4:9;
## n = 1:6;
## [mn, v] = hygestat (m, t, n)
##
## [mn, v] = hygestat (m, t, 2)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/Hypergeometric_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [mn, v] = hygestat (m, t, n)

  # Check arguments
  if (nargin != 3)
    usage ("[m, v] = hygestat (m, t, n)");
  endif

  if (! isempty (m) && ! ismatrix (m))
    error ("hygestat: m must be a numeric matrix");
  endif
  if (! isempty (t) && ! ismatrix (t))
    error ("hygestat: t must be a numeric matrix");
  endif
  if (! isempty (n) && ! ismatrix (n))
    error ("hygestat: n must be a numeric matrix");
  endif

  if (! isscalar (m) || ! isscalar (t) || ! isscalar (n))
    [retval, m, t, n] = common_size (m, t, n);
    if (retval > 0)
      error ("hygestat: m, t and n must be of common size or scalar");
    endif
  endif

  # Calculate moments
  mn = (n .* m) ./ t;
  v = (n .* (m ./ t) .* (1 - m ./ t) .* (t - n)) ./ (t - 1);

  # Continue argument check
  k = find (! (m >= 0) | ! (t >= 0) | ! (n > 0) | ! (m == round (m)) | ! (t == round (t)) | ! (n == round (n)) | ! (m <= t) | ! (n <= t));
  if (any (k))
    mn (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! m = 0:5;
%! t = 4:9;
%! n = 1:6;
%! [mn, v] = hygestat (m, t, n);
%! expected_mn = [0.0000, 0.4000, 1.0000, 1.7143, 2.5000, 3.3333];
%! expected_v = [0.0000, 0.2400, 0.4000, 0.4898, 0.5357, 0.5556];
%! assert (mn, expected_mn, 0.001);
%! assert (v, expected_v, 0.001);

%!test
%! m = 0:5;
%! t = 4:9;
%! [mn, v] = hygestat (m, t, 2);
%! expected_mn = [0.0000, 0.4000, 0.6667, 0.8571, 1.0000, 1.1111];
%! expected_v = [0.0000, 0.2400, 0.3556, 0.4082, 0.4286, 0.4321];
%! assert (mn, expected_mn, 0.001);
%! assert (v, expected_v, 0.001);
