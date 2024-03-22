## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{y} =} fpdf (@var{x}, @var{df1}, @var{df2})
##
## @math{F}-probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the @math{F}-distribution with @var{df1} and @var{df2} degrees of freedom.
## The size of @var{y} is the common size of @var{x}, @var{df1}, and @var{df2}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## Further information about the @math{F}-distribution can be found at
## @url{https://en.wikipedia.org/wiki/F-distribution}
##
## @seealso{fcdf, finv, frnd, fstat}
## @end deftypefn

function y = fpdf (x, df1, df2)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("fpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, DF1, and DF2
  if (! isscalar (x) ||! isscalar (df1) || ! isscalar (df2))
    [retval, x, df1, df2] = common_size (x, df1, df2);
    if (retval > 0)
      error ("fpdf: X, DF1, and DF2 must be of common size or scalars.");
    endif
  endif

  ## Check for X, DF1, and DF2 being reals
  if (iscomplex (x) || iscomplex (df1) || iscomplex (df2))
    error ("fpdf: X, DF1, and DF2 must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (df1, "single") || isa (df2, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  k = isnan (x) | !(df1 > 0) | !(df1 < Inf) | !(df2 > 0) | !(df2 < Inf);
  y(k) = NaN;

  k = (x > 0) & (x < Inf) & (df1 > 0) & (df1 < Inf) & (df2 > 0) & (df2 < Inf);
  if (isscalar (df1) && isscalar (df2))
    tmp = df1 / df2 * x(k);
    y(k) = (exp ((df1/2 - 1) * log (tmp) ...
                   - ((df1 + df2) / 2) * log (1 + tmp)) ...
              * (df1 / df2) ./ beta (df1/2, df2/2));
  else
    tmp = df1(k) .* x(k) ./ df2(k);
    y(k) = (exp ((df1(k)/2 - 1) .* log (tmp) ...
                   - ((df1(k) + df2(k)) / 2) .* log (1 + tmp)) ...
              .* (df1(k) ./ df2(k)) ./ beta (df1(k)/2, df2(k)/2));
  endif

endfunction

%!demo
%! ## Plot various PDFs from the F distribution
%! x = 0.01:0.01:4;
%! y1 = fpdf (x, 1, 1);
%! y2 = fpdf (x, 2, 1);
%! y3 = fpdf (x, 5, 2);
%! y4 = fpdf (x, 10, 1);
%! y5 = fpdf (x, 100, 100);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c", x, y5, "-m")
%! grid on
%! ylim ([0, 2.5])
%! legend ({"df1 = 1, df2 = 2", "df1 = 2, df2 = 1", ...
%!          "df1 = 5, df2 = 2", "df1 = 10, df2 = 1", ...
%!          "df1 = 100, df2 = 100"}, "location", "northeast")
%! title ("F PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 0.5 1 2];
%! y = [0 0 4/9 1/4 1/9];
%!assert (fpdf (x, 2*ones (1,5), 2*ones (1,5)), y, eps)
%!assert (fpdf (x, 2, 2*ones (1,5)), y, eps)
%!assert (fpdf (x, 2*ones (1,5), 2), y, eps)
%!assert (fpdf (x, [0 NaN Inf 2 2], 2), [NaN NaN NaN y(4:5)], eps)
%!assert (fpdf (x, 2, [0 NaN Inf 2 2]), [NaN NaN NaN y(4:5)], eps)
%!assert (fpdf ([x, NaN], 2, 2), [y, NaN], eps)
%!test #F (x, 1, df1) == T distribution (sqrt (x), df1) / sqrt (x)
%! rand ("seed", 1234);    # for reproducibility
%! xr = rand (10,1);
%! xr = xr(x > 0.1 & x < 0.9);
%! yr = tpdf (sqrt (xr), 2) ./ sqrt (xr);
%! assert (fpdf (xr, 1, 2), yr, 5*eps);

## Test class of input preserved
%!assert (fpdf (single ([x, NaN]), 2, 2), single ([y, NaN]), eps ("single"))
%!assert (fpdf ([x, NaN], single (2), 2), single ([y, NaN]), eps ("single"))
%!assert (fpdf ([x, NaN], 2, single (2)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<fpdf: function called with too few input arguments.> fpdf ()
%!error<fpdf: function called with too few input arguments.> fpdf (1)
%!error<fpdf: function called with too few input arguments.> fpdf (1,2)
%!error<fpdf: X, DF1, and DF2 must be of common size or scalars.> ...
%! fpdf (ones (3), ones (2), ones (2))
%!error<fpdf: X, DF1, and DF2 must be of common size or scalars.> ...
%! fpdf (ones (2), ones (3), ones (2))
%!error<fpdf: X, DF1, and DF2 must be of common size or scalars.> ...
%! fpdf (ones (2), ones (2), ones (3))
%!error<fpdf: X, DF1, and DF2 must not be complex.> fpdf (i, 2, 2)
%!error<fpdf: X, DF1, and DF2 must not be complex.> fpdf (2, i, 2)
%!error<fpdf: X, DF1, and DF2 must not be complex.> fpdf (2, 2, i)
