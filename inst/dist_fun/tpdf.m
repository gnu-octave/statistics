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
## @deftypefn {statistics} {@var{p} =} tpdf (@var{x}, @var{df})
##
## Student's T probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Student's T distribution with @var{df} degrees of freedom.  The size
## of @var{y} is the common size of @var{x} and @var{df}.  A scalar input
## functions as a constant matrix of the same size as the other input.
##
## Further information about the Student's T distribution can be found at
## @url{https://en.wikipedia.org/wiki/Student%27s_t-distribution}
##
## @seealso{tcdf, tpdf, trnd, tstat}
## @end deftypefn

function y = tpdf (x, df)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("tpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X and DF
  if (! isscalar (x) || ! isscalar (df))
    [retval, x, df] = common_size (x, df);
    if (retval > 0)
      error ("tpdf: X and DF must be of common size or scalars.");
    endif
  endif

  ## Check for X and DF being reals
  if (iscomplex (x) || iscomplex (df))
    error ("tpdf: X and DF must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (df, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  k = isnan (x) | ! (df > 0);
  y(k) = NaN;

  k = isfinite (x) & (df > 0) & (df < Inf);
  kinf = isfinite (x) & isinf (df);
  if (any (k))
    y(k) = exp (- (df(k) + 1) .* log (1 + x(k) .^ 2 ./ df(k)) / 2) ./ ...
           (sqrt (df(k)) .* beta (df(k)/2, 1/2));
  endif
  if (any (kinf))
    y(kinf) = normpdf (x(kinf));
  endif

endfunction

%!demo
%! ## Plot various PDFs from the Student's T distribution
%! x = -5:0.01:5;
%! y1 = tpdf (x, 1);
%! y2 = tpdf (x, 2);
%! y3 = tpdf (x, 5);
%! y4 = tpdf (x, Inf);
%! plot (x, y1, "-b", x, y2, "g", x, y3, "-r", x, y4, "-m")
%! grid on
%! xlim ([-5, 5])
%! ylim ([0, 0.41])
%! legend ({"df = 1", "df = 2", ...
%!          "df = 5", 'df = \infty'}, "location", "northeast")
%! title ("Student's T PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!test
%! x = rand (10,1);
%! y = 1./(pi * (1 + x.^2));
%! assert (tpdf (x, 1), y, 5*eps);
%!shared x, y
%! x = [-Inf 0 0.5 1 Inf];
%! y = 1./(pi * (1 + x.^2));
%!assert (tpdf (x, ones (1,5)), y, eps)
%!assert (tpdf (x, 1), y, eps)
%!assert (tpdf (x, [0 NaN 1 1 1]), [NaN NaN y(3:5)], eps)
%!assert (tpdf (x, Inf), normpdf (x))

## Test class of input preserved
%!assert (tpdf ([x, NaN], 1), [y, NaN], eps)
%!assert (tpdf (single ([x, NaN]), 1), single ([y, NaN]), eps ("single"))
%!assert (tpdf ([x, NaN], single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<tpdf: function called with too few input arguments.> tpdf ()
%!error<tpdf: function called with too few input arguments.> tpdf (1)
%!error<tpdf: X and DF must be of common size or scalars.> ...
%! tpdf (ones (3), ones (2))
%!error<tpdf: X and DF must be of common size or scalars.> ...
%! tpdf (ones (2), ones (3))
%!error<tpdf: X and DF must not be complex.> tpdf (i, 2)
%!error<tpdf: X and DF must not be complex.> tpdf (2, i)
