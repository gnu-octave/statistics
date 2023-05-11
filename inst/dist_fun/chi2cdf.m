## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{p} =} chi2cdf (@var{x}, @var{df})
## @deftypefnx {statistics} {@var{p} =} chi2cdf (@var{x}, @var{df}, @qcode{"upper"})
##
## Chi-squared cumulative distribution function.
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the chi-squared distribution with @var{df} degrees of freedom.  The
## chi-squared density function with @var{df} degrees of freedom is the same as
## a gamma density function with parameters @qcode{@var{df}/2} and @qcode{2}.
##
## The size of @var{p} is the common size of @var{x} and @var{df}. A scalar
## input functions as a constant matrix of the same size as the other input.
##
## @code{@var{p} = chi2cdf (@var{x}, @var{df}, "upper")} computes the upper tail
## probability of the chi-squared distribution with @var{df} degrees of freedom,
## at the values in @var{x}.
##
## Further information about the chi-squared distribution can be found at
## @url{https://en.wikipedia.org/wiki/Chi-squared_distribution}
##
## @seealso{chi2inv, chi2pdf, chi2rnd, chi2stat}
## @end deftypefn

function p = chi2cdf (x, df, uflag)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("chi2cdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 2 && ! strcmpi (uflag, "upper"))
    error ("chi2cdf: invalid argument for upper tail.");
  else
    uflag = [];
  endif

  ## Check for common size of X and DF
  if (! isscalar (x) || ! isscalar (df))
    [err, x, df] = common_size (x, df);
    if (err > 0)
      error ("chi2cdf: X and DF must be of common size or scalars.");
    endif
  endif

  ## Check for X and DF being reals
  if (iscomplex (x) || iscomplex (df))
    error ("chi2cdf: X and DF must not be complex.");
  endif

  ## Compute chi-squared CDF
  p = gamcdf (x, df/2, 2, uflag);

endfunction

%!demo
%! ## Plot various CDFs from the chi-squared distribution
%! x = 0:0.01:8;
%! p1 = chi2cdf (x, 1);
%! p2 = chi2cdf (x, 2);
%! p3 = chi2cdf (x, 3);
%! p4 = chi2cdf (x, 4);
%! p5 = chi2cdf (x, 6);
%! p6 = chi2cdf (x, 9);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", ...
%!       x, p4, "-c", x, p5, "-m", x, p6, "-y")
%! grid on
%! xlim ([0, 8])
%! legend ({"df = 1", "df = 2", "df = 3", ...
%!          "df = 4", "df = 6", "df = 9"}, "location", "southeast")
%! title ("Chi-squared CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, p, u
%! x = [-1, 0, 0.5, 1, 2];
%! p = [0, (1 - exp (-x(2:end) / 2))];
%! u = [1, 0, NaN, 0.3934693402873666, 0.6321205588285577];
%!assert (chi2cdf (x, 2 * ones (1,5)), p, eps)
%!assert (chi2cdf (x, 2), p, eps)
%!assert (chi2cdf (x, 2 * [1, 0, NaN, 1, 1]), [p(1), 1, NaN, p(4:5)], eps)
%!assert (chi2cdf (x, 2 * [1, 0, NaN, 1, 1], "upper"), ...
%!                        [p(1), 1, NaN, u(4:5)], eps)
%!assert (chi2cdf ([x(1:2), NaN, x(4:5)], 2), [p(1:2), NaN, p(4:5)], eps)

## Test class of input preserved
%!assert (chi2cdf ([x, NaN], 2), [p, NaN], eps)
%!assert (chi2cdf (single ([x, NaN]), 2), single ([p, NaN]), eps ("single"))
%!assert (chi2cdf ([x, NaN], single (2)), single ([p, NaN]), eps ("single"))

## Test input validation
%!error<chi2cdf: function called with too few input arguments.> chi2cdf ()
%!error<chi2cdf: function called with too few input arguments.> chi2cdf (1)
%!error<chi2cdf: function called with too many inputs> chi2cdf (1, 2, 3, 4)
%!error<chi2cdf: invalid argument for upper tail.> chi2cdf (1, 2, 3)
%!error<chi2cdf: invalid argument for upper tail.> chi2cdf (1, 2, "uper")
%!error<chi2cdf: X and DF must be of common size or scalars.> ...
%! chi2cdf (ones (3), ones (2))
%!error<chi2cdf: X and DF must be of common size or scalars.> ...
%! chi2cdf (ones (2), ones (3))
%!error<chi2cdf: X and DF must not be complex.> chi2cdf (i, 2)
%!error<chi2cdf: X and DF must not be complex.> chi2cdf (2, i)
