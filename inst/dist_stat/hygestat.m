## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{mn}, @var{v}] =} hygestat (@var{m}, @var{k}, @var{n})
##
## Compute statistics of the hypergeometric distribution.
##
## @code{[@var{mn}, @var{v}] = hygestat (@var{m}, @var{k}, @var{n})} returns the
## mean and variance of the hypergeometric distribution parameters @var{m},
## @var{k}, and @var{n}.
##
## @itemize
## @item
## @var{m} is the total size of the population of the hypergeometric
## distribution. The elements of @var{m} must be positive natural numbers.
##
## @item
## @var{k} is the number of marked items of the hypergeometric distribution.
## The elements of @var{k} must be natural numbers.
##
## @item
## @var{n} is the size of the drawn sample of the hypergeometric
## distribution. The elements of @var{n} must be positive natural numbers.
## @end itemize
##
## The size of @var{mn} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the hypergeometric distribution can be found at
## @url{https://en.wikipedia.org/wiki/Hypergeometric_distribution}
##
## @seealso{hygecdf, hygeinv, hygepdf, hygernd}
## @end deftypefn

function [mn, v] = hygestat (m, k, n)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("hygestat: function called with too few input arguments.");
  endif

  ## Check for M, K, and N being numeric
  if (! (isnumeric (m) && isnumeric (k) && isnumeric (n)))
    error ("hygestat: M, K, and N must be numeric.");
  endif

  ## Check for M, K, and N being real
  if (iscomplex (m) || iscomplex (k) || iscomplex (n))
    error ("hygestat: M, K, and N must not be complex.");
  endif

  ## Check for common size of M, K, and N
  if (! isscalar (m) || ! isscalar (k) || ! isscalar (n))
    [retval, m, k, n] = common_size (m, k, n);
    if (retval > 0)
      error ("hygestat: M, K, and N must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  mn = (n .* k) ./ m;
  v = (n .* (k ./ m) .* (1 - k ./ m) .* (m - n)) ./ (m - 1);

  ## Continue argument check
  is_nan = find (! (m >= 0) | ! (k >= 0) | ! (n > 0) | ! (m == round (m)) | ...
                 ! (k == round (k)) | ! (n == round (n)) | ! (k <= m) | ...
                 ! (n <= m));
  if (any (is_nan))
    mn(is_nan) = NaN;
    v(is_nan) = NaN;
  endif

endfunction

## Input validation tests
%!error<hygestat: function called with too few input arguments.> hygestat ()
%!error<hygestat: function called with too few input arguments.> hygestat (1)
%!error<hygestat: function called with too few input arguments.> hygestat (1, 2)
%!error<hygestat: M, K, and N must be numeric.> hygestat ({}, 2, 3)
%!error<hygestat: M, K, and N must be numeric.> hygestat (1, "", 3)
%!error<hygestat: M, K, and N must be numeric.> hygestat (1, 2, "")
%!error<hygestat: M, K, and N must not be complex.> hygestat (i, 2, 3)
%!error<hygestat: M, K, and N must not be complex.> hygestat (1, i, 3)
%!error<hygestat: M, K, and N must not be complex.> hygestat (1, 2, i)
%!error<hygestat: M, K, and N must be of common size or scalars.> ...
%! hygestat (ones (3), ones (2), 3)
%!error<hygestat: M, K, and N must be of common size or scalars.> ...
%! hygestat (ones (2), 2, ones (3))
%!error<hygestat: M, K, and N must be of common size or scalars.> ...
%! hygestat (1, ones (2), ones (3))

## Output validation tests
%!test
%! m = 4:9;
%! k = 0:5;
%! n = 1:6;
%! [mn, v] = hygestat (m, k, n);
%! expected_mn = [0.0000, 0.4000, 1.0000, 1.7143, 2.5000, 3.3333];
%! expected_v = [0.0000, 0.2400, 0.4000, 0.4898, 0.5357, 0.5556];
%! assert (mn, expected_mn, 0.001);
%! assert (v, expected_v, 0.001);
%!test
%! m = 4:9;
%! k = 0:5;
%! [mn, v] = hygestat (m, k, 2);
%! expected_mn = [0.0000, 0.4000, 0.6667, 0.8571, 1.0000, 1.1111];
%! expected_v = [0.0000, 0.2400, 0.3556, 0.4082, 0.4286, 0.4321];
%! assert (mn, expected_mn, 0.001);
%! assert (v, expected_v, 0.001);
