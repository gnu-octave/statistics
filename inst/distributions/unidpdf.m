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
## @deftypefn  {statistics} @var{y} = unidpdf (@var{x}, @var{df})
##
## Discrete uniform probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of a discrete uniform distribution which assumes the integer
## values 1--@var{df} with equal probability.  The size of @var{p} is the common
## size of @var{x} and @var{df}.  A scalar input functions as a constant matrix
## of the same size as the other inputs.
##
## Warning: The underlying implementation uses the double class and will only
## be accurate for @var{df} < @code{flintmax} (@w{@math{2^{53}}} on
## IEEE 754 compatible systems).
##
## @seealso{unidcdf, unidinv, unidrnd, unidstat}
## @end deftypefn

function y = unidpdf (x, df)

  if (nargin != 2)
    print_usage ();
  endif

  if (! isscalar (x) || ! isscalar (df))
    [retval, x, df] = common_size (x, df);
    if (retval > 0)
      error ("unidpdf: X and N must be of common size or scalars");
    endif
  endif

  if (iscomplex (x) || iscomplex (df))
    error ("unidpdf: X and N must not be complex");
  endif

  if (isa (x, "single") || isa (df, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  k = isnan (x) | ! (df > 0 & df == fix (df));
  y(k) = NaN;

  k = ! k & (x >= 1) & (x <= df) & (x == fix (x));
  y(k) = 1 ./ df(k);

endfunction


%!shared x,y
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
%!error unidpdf ()
%!error unidpdf (1)
%!error unidpdf (1,2,3)
%!error unidpdf (ones (3), ones (2))
%!error unidpdf (ones (2), ones (3))
%!error unidpdf (i, 2)
%!error unidpdf (2, i)
