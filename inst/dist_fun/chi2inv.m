## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{x} =} chi2inv (@var{p}, @var{df})
##
## Inverse of the Chi-square cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the chi-square distribution with @var{df} degrees of freedom.
## The size of @var{x} is the common size of @var{p} and @var{df}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## @seealso{chi2cdf, chi2pdf, chi2rnd, chi2stat}
## @end deftypefn

function x = chi2inv (p, df)

  if (nargin != 2)
    print_usage ();
  endif

  if (! isscalar (p) || ! isscalar (df))
    [retval, p, df] = common_size (p, df);
    if (retval > 0)
      error ("chi2inv: P and DF must be of common size or scalars.");
    endif
  endif

  if (iscomplex (p) || iscomplex (df))
    error ("chi2inv: P and DF must not be complex.");
  endif

  x = gaminv (p, df/2, 2);

endfunction


%!shared p
%! p = [-1 0 0.3934693402873666 1 2];
%!assert (chi2inv (p, 2*ones (1,5)), [NaN 0 1 Inf NaN], 5*eps)
%!assert (chi2inv (p, 2), [NaN 0 1 Inf NaN], 5*eps)
%!assert (chi2inv (p, 2*[0 1 NaN 1 1]), [NaN 0 NaN Inf NaN], 5*eps)
%!assert (chi2inv ([p(1:2) NaN p(4:5)], 2), [NaN 0 NaN Inf NaN], 5*eps)

## Test class of input preserved
%!assert (chi2inv ([p, NaN], 2), [NaN 0 1 Inf NaN NaN], 5*eps)
%!assert (chi2inv (single ([p, NaN]), 2), single ([NaN 0 1 Inf NaN NaN]), 5*eps ("single"))
%!assert (chi2inv ([p, NaN], single (2)), single ([NaN 0 1 Inf NaN NaN]), 5*eps ("single"))

## Test input validation
%!error chi2inv ()
%!error chi2inv (1)
%!error chi2inv (1,2,3)
%!error chi2inv (ones (3), ones (2))
%!error chi2inv (ones (2), ones (3))
%!error chi2inv (i, 2)
%!error chi2inv (2, i)
