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
## @deftypefn  {statistics} @var{y} = fpdf (@var{x}, @var{df1}, @var{df2})
##
## F probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the F distribution with @var{df1} and @var{df2} degrees of
## freedom.  The size of @var{y} is the common size of @var{x}, @var{df1}, and
## @var{df2}.  A scalar input functions as a constant matrix of the same size as
## the other inputs.
##
## @seealso{fcdf, finv, frnd, fstat}
## @end deftypefn

function y = fpdf (x, df1, df2)

  if (nargin != 3)
    print_usage ();
  endif

  if (! isscalar (df1) || ! isscalar (df2))
    [retval, x, df1, df2] = common_size (x, df1, df2);
    if (retval > 0)
      error ("fpdf: X, DF1, and DF2 must be of common size or scalars.");
    endif
  endif

  if (iscomplex (x) || iscomplex (df1) || iscomplex (df2))
    error ("fpdf: X, DF1, and DF2 must not be complex.");
  endif

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


## F (x, 1, df1) == T distribution (sqrt (x), df1) / sqrt (x)
%!test
%! x = rand (10,1);
%! x = x(x > 0.1 & x < 0.9);
%! y = tpdf (sqrt (x), 2) ./ sqrt (x);
%! assert (fpdf (x, 1, 2), y, 5*eps);

%!shared x,y
%! x = [-1 0 0.5 1 2];
%! y = [0 0 4/9 1/4 1/9];
%!assert (fpdf (x, 2*ones (1,5), 2*ones (1,5)), y, eps)
%!assert (fpdf (x, 2, 2*ones (1,5)), y, eps)
%!assert (fpdf (x, 2*ones (1,5), 2), y, eps)
%!assert (fpdf (x, [0 NaN Inf 2 2], 2), [NaN NaN NaN y(4:5)], eps)
%!assert (fpdf (x, 2, [0 NaN Inf 2 2]), [NaN NaN NaN y(4:5)], eps)
%!assert (fpdf ([x, NaN], 2, 2), [y, NaN], eps)

## Test class of input preserved
%!assert (fpdf (single ([x, NaN]), 2, 2), single ([y, NaN]), eps ("single"))
%!assert (fpdf ([x, NaN], single (2), 2), single ([y, NaN]), eps ("single"))
%!assert (fpdf ([x, NaN], 2, single (2)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error fpdf ()
%!error fpdf (1)
%!error fpdf (1,2)
%!error fpdf (1,2,3,4)
%!error fpdf (ones (3), ones (2), ones (2))
%!error fpdf (ones (2), ones (3), ones (2))
%!error fpdf (ones (2), ones (2), ones (3))
%!error fpdf (i, 2, 2)
%!error fpdf (2, i, 2)
%!error fpdf (2, 2, i)
