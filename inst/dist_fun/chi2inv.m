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
## Inverse of the chi-squared cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the chi-squared distribution with @var{df} degrees of freedom.
## The size of @var{x} is the common size of @var{p} and @var{df}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## Further information about the chi-squared distribution can be found at
## @url{https://en.wikipedia.org/wiki/Chi-squared_distribution}
##
## @seealso{chi2cdf, chi2pdf, chi2rnd, chi2stat}
## @end deftypefn

function x = chi2inv (p, df)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("chi2inv: function called with too few input arguments.");
  endif

  ## Check for common size of P and DF
  if (! isscalar (p) || ! isscalar (df))
    [retval, p, df] = common_size (p, df);
    if (retval > 0)
      error ("chi2inv: P and DF must be of common size or scalars.");
    endif
  endif

  ## Check for X and DF being reals
  if (iscomplex (p) || iscomplex (df))
    error ("chi2inv: P and DF must not be complex.");
  endif

  ## Compute chi-squared iCDF
  x = gaminv (p, df/2, 2);

endfunction

%!demo
%! ## Plot various iCDFs from the chi-squared distribution
%! p = 0.001:0.001:0.999;
%! x1 = chi2inv (p, 1);
%! x2 = chi2inv (p, 2);
%! x3 = chi2inv (p, 3);
%! x4 = chi2inv (p, 4);
%! x5 = chi2inv (p, 6);
%! x6 = chi2inv (p, 9);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", ...
%!       p, x4, "-c", p, x5, "-m", p, x6, "-y")
%! grid on
%! ylim ([0, 8])
%! legend ({"df = 1", "df = 2", "df = 3", ...
%!          "df = 4", "df = 6", "df = 9"}, "location", "northwest")
%! title ("Chi-squared iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
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
%!error<chi2inv: function called with too few input arguments.> chi2inv ()
%!error<chi2inv: function called with too few input arguments.> chi2inv (1)
%!error<chi2inv: function called with too many inputs> chi2inv (1,2,3)
%!error<chi2inv: P and DF must be of common size or scalars.> ...
%! chi2inv (ones (3), ones (2))
%!error<chi2inv: P and DF must be of common size or scalars.> ...
%! chi2inv (ones (2), ones (3))
%!error<chi2inv: P and DF must not be complex.> chi2inv (i, 2)
%!error<chi2inv: P and DF must not be complex.> chi2inv (2, i)
