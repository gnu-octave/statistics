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
## @deftypefn {statistics} @var{p} = tpdf (@var{x}, @var{df})
##
## Student's T probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Student's T distribution with @var{df} degrees of freedom.
##
## The size of @var{y} is the common size of @var{x} and @var{df}. A scalar
## input functions as a constant matrix of the same size as the other input.
##
## @seealso{tcdf, tpdf, trnd, tstat}
## @end deftypefn

function y = tpdf (x, df)

  if (nargin != 2)
    print_usage ();
  endif

  if (! isscalar (df))
    [retval, x, df] = common_size (x, df);
    if (retval > 0)
      error ("tpdf: X and DF must be of common size or scalars.");
    endif
  endif

  if (iscomplex (x) || iscomplex (df))
    error ("tpdf: X and DF must not be complex");
  endif

  if (isa (x, "single") || isa (df, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  k = isnan (x) | !(df > 0) | !(df < Inf);
  y(k) = NaN;

  k = isfinite (x) & (df > 0) & (df < Inf);
  if (isscalar (df))
    y(k) = (exp (- (df + 1) * log (1 + x(k) .^ 2 / df)/2)
               / (sqrt (df) * beta (df/2, 1/2)));
  else
    y(k) = (exp (- (df(k) + 1) .* log (1 + x(k) .^ 2 ./ df(k))/2)
              ./ (sqrt (df(k)) .* beta (df(k)/2, 1/2)));
  endif

endfunction


%!test
%! x = rand (10,1);
%! y = 1./(pi * (1 + x.^2));
%! assert (tpdf (x, 1), y, 5*eps);

%!shared x,y
%! x = [-Inf 0 0.5 1 Inf];
%! y = 1./(pi * (1 + x.^2));
%!assert (tpdf (x, ones (1,5)), y, eps)
%!assert (tpdf (x, 1), y, eps)
%!assert (tpdf (x, [0 NaN 1 1 1]), [NaN NaN y(3:5)], eps)

## Test class of input preserved
%!assert (tpdf ([x, NaN], 1), [y, NaN], eps)
%!assert (tpdf (single ([x, NaN]), 1), single ([y, NaN]), eps ("single"))
%!assert (tpdf ([x, NaN], single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error tpdf ()
%!error tpdf (1)
%!error tpdf (1,2,3)
%!error tpdf (ones (3), ones (2))
%!error tpdf (ones (2), ones (3))
%!error tpdf (i, 2)
%!error tpdf (2, i)
