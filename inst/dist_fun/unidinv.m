## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 2007-2016 David Bateman
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{x} =} unidinv (@var{p}, @var{N})
##
## Inverse of the discrete uniform cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the discrete uniform distribution with parameter @var{N}, which corresponds
## to the maximum observable value.  @code{unidinv} assumes the integer values
## in the range @math{[1,N]} with equal probability.  The size of @var{x} is the
## common size of @var{p} and @var{N}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## The maximum observable values in @var{N} must be positive integers, otherwise
## @qcode{NaN} is returned.
##
## Warning: The underlying implementation uses the double class and will only
## be accurate for @var{N} < @code{flintmax} (@w{@math{2^{53}}} on
## IEEE 754 compatible systems).
##
## Further information about the discrete uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Discrete_uniform_distribution}
##
## @seealso{unidcdf, unidpdf, unidrnd, unidfit, unidstat}
## @end deftypefn

function x = unidinv (p, N)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("unidinv: function called with too few input arguments.");
  endif

  ## Check for common size of P and N
  if (! isscalar (p) || ! isscalar (N))
    [retval, p, N] = common_size (p, N);
    if (retval > 0)
      error ("unidinv: P and N must be of common size or scalars.");
    endif
  endif

  ## Check for P and N being reals
  if (iscomplex (p) || iscomplex (N))
    error ("unidinv: P and N must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (N, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## For Matlab compatibility, unidinv(0) = NaN
  k = (p > 0) & (p <= 1) & (N > 0 & N == fix (N));
  x(k) = floor (p(k) .* N(k));

endfunction

%!demo
%! ## Plot various iCDFs from the discrete uniform distribution
%! p = 0.001:0.001:0.999;
%! x1 = unidinv (p, 5);
%! x2 = unidinv (p, 9);
%! plot (p, x1, "-b", p, x2, "-g")
%! grid on
%! xlim ([0, 1])
%! ylim ([0, 10])
%! legend ({"N = 5", "N = 9"}, "location", "northwest")
%! title ("Discrete uniform iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (unidinv (p, 10*ones (1,5)), [NaN NaN 5 10 NaN], eps)
%!assert (unidinv (p, 10), [NaN NaN 5 10 NaN], eps)
%!assert (unidinv (p, 10*[0 1 NaN 1 1]), [NaN NaN NaN 10 NaN], eps)
%!assert (unidinv ([p(1:2) NaN p(4:5)], 10), [NaN NaN NaN 10 NaN], eps)

## Test class of input preserved
%!assert (unidinv ([p, NaN], 10), [NaN NaN 5 10 NaN NaN], eps)
%!assert (unidinv (single ([p, NaN]), 10), single ([NaN NaN 5 10 NaN NaN]), eps)
%!assert (unidinv ([p, NaN], single (10)), single ([NaN NaN 5 10 NaN NaN]), eps)

## Test input validation
%!error<unidinv: function called with too few input arguments.> unidinv ()
%!error<unidinv: function called with too few input arguments.> unidinv (1)
%!error<unidinv: P and N must be of common size or scalars.> ...
%! unidinv (ones (3), ones (2))
%!error<unidinv: P and N must be of common size or scalars.> ...
%! unidinv (ones (2), ones (3))
%!error<unidinv: P and N must not be complex.> unidinv (i, 2)
%!error<unidinv: P and N must not be complex.> unidinv (2, i)
