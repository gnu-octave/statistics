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
## @deftypefn  {statistics} {@var{p} =} unidcdf (@var{x}, @var{df})
## @deftypefnx {statistics} {@var{p} =} unidcdf (@var{x}, @var{df}, "upper")
##
## Discrete uniform cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of a discrete uniform distribution which assumes the integer
## values 1--@var{df} with equal probability.  The size of @var{p} is the common
## size of @var{x} and @var{df}.  A scalar input functions as a constant matrix
## of the same size as the other inputs.
##
## @code{[@dots{}] = logncdf (@dots{}, "upper")} computes the upper tail
## probability of the lognormal distribution.
##
## @seealso{unidinv, unidpdf, unidrnd, unidstat}
## @end deftypefn

function p = unidcdf (x, df, uflag)

  ## Check for valid number of input arguments
  if (nargin < 2 || nargin > 3)
    error ("unidcdf: invalid number of input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 2 && strcmpi (uflag, "upper"))
    uflag = true;
  elseif (nargin > 2  && ! strcmpi (uflag, "upper"))
    error ("unidcdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Check for common size of X and DF
  if (! isscalar (x) || ! isscalar (df))
    [retval, x, df] = common_size (x, df);
    if (retval > 0)
      error ("unidcdf: X and DF must be of common size or scalars.");
    endif
  endif

  ## Check for X and DF being reals
  if (iscomplex (x) || iscomplex (df))
    error ("unidcdf: X and DF must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (df, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Return 1 for X >= DF
  p(x >= df) = 1;

  ## Floor X
  xf = floor (x);

  ## Compute uniform discrete CDF
  k = find (xf >= 1 & xf <= df);
  if any(k)
    p(k) = xf(k) ./ df(k);
  endif

  ## Check for NaNs or floored DF <= 0
  is_nan = isnan (x) | ! (df > 0 & df == fix (df));
  if (any (is_nan(:)))
    p(is_nan) = NaN;
  endif

  p(df < 1 | round(df) != df) = NaN;

  if (uflag)  # Compute upper tail
    p = 1 - unidcdf (x, df);
  endif

endfunction


%!shared x,y
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
%!error unidcdf ()
%!error unidcdf (1)
%!error unidcdf (1, 2, 3)
%!error unidcdf (1, 2, "upper", 4)
%!error unidcdf (ones (3), ones (2))
%!error unidcdf (ones (2), ones (3))
%!error unidcdf (i, 2)
%!error unidcdf (2, i)
