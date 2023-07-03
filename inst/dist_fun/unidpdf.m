## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 2007-2016 David Bateman
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{y} =} unidpdf (@var{x}, @var{N})
##
## Discrete uniform probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the discrete uniform distribution with parameter @var{N}, which
## corresponds to the maximum observable value.  @code{unidpdf} assumes the
## integer values in the range @math{[1,N]} with equal probability.  The size of
## @var{x} is the common size of @var{p} and @var{N}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
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
## @seealso{unidcdf, unidinv, unidrnd, unidfit, unidstat}
## @end deftypefn

function y = unidpdf (x, N)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("unidpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X and N
  if (! isscalar (x) || ! isscalar (N))
    [retval, x, N] = common_size (x, N);
    if (retval > 0)
      error ("unidpdf: X and N must be of common size or scalars.");
    endif
  endif

  ## Check for X and N being reals
  if (iscomplex (x) || iscomplex (N))
    error ("unidpdf: X and N must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (N, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  k = isnan (x) | ! (N > 0 & N == fix (N));
  y(k) = NaN;

  k = ! k & (x >= 1) & (x <= N) & (x == fix (x));
  y(k) = 1 ./ N(k);

endfunction

%!demo
%! ## Plot various PDFs from the discrete uniform distribution
%! x = 0:10;
%! y1 = unidpdf (x, 5);
%! y2 = unidpdf (x, 9);
%! plot (x, y1, "*b", x, y2, "*g")
%! grid on
%! xlim ([0, 10])
%! ylim ([0, 0.25])
%! legend ({"N = 5", "N = 9"}, "location", "northeast")
%! title ("Descrete uniform PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 1 2 10 11];
%! y = [0 0 0.1 0.1 0.1 0];
%!assert (unidpdf (x, 10*ones (1,6)), y)
%!assert (unidpdf (x, 10), y)
%!assert (unidpdf (x, 10*[0 NaN 1 1 1 1]), [NaN NaN y(3:6)])
%!assert (unidpdf ([x, NaN], 10), [y, NaN])

## Test class of input preserved
%!assert (unidpdf (single ([x, NaN]), 10), single ([y, NaN]))
%!assert (unidpdf ([x, NaN], single (10)), single ([y, NaN]))

## Test input validation
%!error<unidpdf: function called with too few input arguments.> unidpdf ()
%!error<unidpdf: function called with too few input arguments.> unidpdf (1)
%!error<unidpdf: X and N must be of common size or scalars.> ...
%! unidpdf (ones (3), ones (2))
%!error<unidpdf: X and N must be of common size or scalars.> ...
%! unidpdf (ones (2), ones (3))
%!error<unidpdf: X and N must not be complex.> unidpdf (i, 2)
%!error<unidpdf: X and N must not be complex.> unidpdf (2, i)
