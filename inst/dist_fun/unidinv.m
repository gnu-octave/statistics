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
## @deftypefn  {statistics} {@var{x} =} unidinv (@var{p}, @var{df})
##
## Inverse of the discrete uniform cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the discrete uniform distribution which assumes the integer
## values 1--@var{df} with equal probability.  The size of @var{p} is the common
## size of @var{x} and @var{df}.  A scalar input functions as a constant matrix
## of the same size as the other inputs.
##
## @seealso{unidcdf, unidpdf, unidrnd, unidstat}
## @end deftypefn

function x = unidinv (p, df)

  ## Check for valid number of input arguments
  if (nargin != 2)
    print_usage ();
  endif

  ## Check for common size of X and DF
  if (! isscalar (p) || ! isscalar (df))
    [retval, p, df] = common_size (p, df);
    if (retval > 0)
      error ("unidcdf: P and DF must be of common size or scalars.");
    endif
  endif

  ## Check for X and DF being reals
  if (iscomplex (p) || iscomplex (df))
    error ("unidinv: P and DF must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (df, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## For Matlab compatibility, unidinv(0) = NaN
  k = (p > 0) & (p <= 1) & (df > 0 & df == fix (df));
  x(k) = floor (p(k) .* df(k));

endfunction


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
%!error unidinv ()
%!error unidinv (1)
%!error unidinv (1,2,3)
%!error unidinv (ones (3), ones (2))
%!error unidinv (ones (2), ones (3))
%!error unidinv (i, 2)
%!error unidinv (2, i)
