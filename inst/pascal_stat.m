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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} pascal_stat (@var{n}, @var{p})
## Returns mean and variance of the negative binomial distribution
##
## Arguments are
##
## @itemize
## @item
## @var{n} is the first parameter of the negative binomial distribution. The elements
## of @var{n} must be natural numbers
## @item
## @var{p} is the second parameter of the negative binomial distribution. The
## elements of @var{p} must be probabilities
## @end itemize
## @var{n} and @var{p} must be of common size or one of them must be scalar
##
## Return values are
##
## @itemize
## @item
## @var{m} is the mean of the negative binomial distribution
## @item
## @var{v} is the variance of the negative binomial distribution
## @end itemize
##
## Examples:
##
## @example
## n = 1:4;
## p = 0.2:0.2:0.8;
## [m, v] = pascal_stat (n, p)
##
## [m, v] = pascal_stat (n, 0.5)
## @end example
##
## References:
##
## @itemize
## @item
## @cite{Matlab 7.0 documentation (pdf)}
## @item
## @uref{http://en.wikipedia.org/wiki/Pascal_distribution}
## @end itemize
##
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>

function [m, v] = pascal_stat (n, p)

  # Check arguments
  if (nargin != 2)
    usage ("[m, v] = pascal_stat (n, p)");
  endif

  if (! isempty (n) && ! ismatrix (n))
    error ("pascal_stat: n must be a numeric matrix");
  endif
  if (! isempty (p) && ! ismatrix (p))
    error ("pascal_stat: p must be a numeric matrix");
  endif

  if (! isscalar (n) || ! isscalar (p))
    [retval, n, p] = common_size (n, p);
    if (retval > 0)
      error ("pascal_stat: n and p must be of common size or scalar");
    endif
  endif

  # Calculate moments
  q = 1 - p;
  m = n .* q ./ p;
  v = n .* q ./ (p .^ 2);

  # Continue argument check
  k = find (! (n > 0) | ! (n < Inf) | ! (p > 0) | ! (p < 1));
  if (any (k))
    m (k) = NaN;
    v (k) = NaN;
  endif

endfunction

%!test
%! n = 1:4;
%! p = 0.2:0.2:0.8;
%! [m, v] = pascal_stat (n, p);
%! expected_m = [ 4.0000, 3.0000, 2.0000, 1.0000];
%! expected_v = [20.0000, 7.5000, 3.3333, 1.2500];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);

%!test
%! n = 1:4;
%! [m, v] = pascal_stat (n, 0.5);
%! expected_m = [1, 2, 3, 4];
%! expected_v = [2, 4, 6, 8];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
