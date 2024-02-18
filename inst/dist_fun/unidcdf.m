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
## @deftypefn  {statistics} {@var{p} =} unidcdf (@var{x}, @var{N})
## @deftypefnx {statistics} {@var{p} =} unidcdf (@var{x}, @var{N}, @qcode{"upper"})
##
## Discrete uniform cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of a discrete uniform distribution with parameter @var{N}, which
## corresponds to the maximum observable value.  @code{unidcdf} assumes the
## integer values in the range @math{[1,N]} with equal probability.  The size of
## @var{p} is the common size of @var{x} and @var{N}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## The maximum observable values in @var{N} must be positive integers, otherwise
## @qcode{NaN} is returned.
##
## @code{[@dots{}] = unidcdf (@var{x}, @var{N}, "upper")} computes the upper
## tail probability of the discrete uniform distribution with maximum observable
## value @var{N}, at the values in @var{x}.
##
## Warning: The underlying implementation uses the double class and will only
## be accurate for @var{N} < @code{flintmax} (@w{@math{2^{53}}} on
## IEEE 754 compatible systems).
##
## Further information about the discrete uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Discrete_uniform_distribution}
##
## @seealso{unidinv, unidpdf, unidrnd, unidfit, unidstat}
## @end deftypefn

function p = unidcdf (x, N, uflag)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("unidcdf: function called with too few input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 2 && strcmpi (uflag, "upper"))
    uflag = true;
  elseif (nargin > 2  && ! strcmpi (uflag, "upper"))
    error ("unidcdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Check for common size of X and N
  if (! isscalar (x) || ! isscalar (N))
    [retval, x, N] = common_size (x, N);
    if (retval > 0)
      error ("unidcdf: X and N must be of common size or scalars.");
    endif
  endif

  ## Check for X and N being reals
  if (iscomplex (x) || iscomplex (N))
    error ("unidcdf: X and N must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (N, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Return 1 for X >= N
  p(x >= N) = 1;

  ## Floor X
  xf = floor (x);

  ## Compute uniform discrete CDF
  k = find (xf >= 1 & xf <= N);
  if any(k)
    p(k) = xf(k) ./ N(k);
  endif

  ## Check for NaNs or floored N <= 0
  is_nan = isnan (x) | ! (N > 0 & N == fix (N));
  if (any (is_nan(:)))
    p(is_nan) = NaN;
  endif

  p(N < 1 | round(N) != N) = NaN;

  if (uflag)  # Compute upper tail
    p = 1 - unidcdf (x, N);
  endif

endfunction

%!demo
%! ## Plot various CDFs from the discrete uniform distribution
%! x = 0:10;
%! p1 = unidcdf (x, 5);
%! p2 = unidcdf (x, 9);
%! plot (x, p1, "*b", x, p2, "*g")
%! grid on
%! xlim ([0, 10])
%! ylim ([0, 1])
%! legend ({"N = 5", "N = 9"}, "location", "southeast")
%! title ("Discrete uniform CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [0 1 2.5 10 11];
%! y = [0, 0.1 0.2 1.0 1.0];
%!assert (unidcdf (x, 10*ones (1,5)), y)
%!assert (unidcdf (x, 10*ones (1,5), "upper"), 1 - y)
%!assert (unidcdf (x, 10), y)
%!assert (unidcdf (x, 10, "upper"), 1 - y)
%!assert (unidcdf (x, 10*[0 1 NaN 1 1]), [NaN 0.1 NaN y(4:5)])
%!assert (unidcdf ([x(1:2) NaN Inf x(5)], 10), [y(1:2) NaN 1 y(5)])

## Test class of input preserved
%!assert (unidcdf ([x, NaN], 10), [y, NaN])
%!assert (unidcdf (single ([x, NaN]), 10), single ([y, NaN]))
%!assert (unidcdf ([x, NaN], single (10)), single ([y, NaN]))

## Test input validation
%!error<unidcdf: function called with too few input arguments.> unidcdf ()
%!error<unidcdf: function called with too few input arguments.> unidcdf (1)
%!error<unidcdf: invalid argument for upper tail.> unidcdf (1, 2, 3)
%!error<unidcdf: invalid argument for upper tail.> unidcdf (1, 2, "tail")
%!error<unidcdf: X and N must be of common size or scalars.> ...
%! unidcdf (ones (3), ones (2))
%!error<unidcdf: X and N must be of common size or scalars.> ...
%! unidcdf (ones (2), ones (3))
%!error<unidcdf: X and N must not be complex.> unidcdf (i, 2)
%!error<unidcdf: X and N must not be complex.> unidcdf (2, i)
