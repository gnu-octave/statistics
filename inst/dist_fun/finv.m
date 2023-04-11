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
## @deftypefn  {statistics} {@var{x} =} finv (@var{p}, @var{df1}, @var{df2})
##
## Inverse of the F cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the F distribution with @var{df1} and @var{df2} degrees of
## freedom.  The size of @var{x} is the common size of @var{p}, @var{df1}, and
## @var{df2}.  A scalar input functions as a constant matrix of the same size as
## the other inputs.
##
## @seealso{fcdf, fpdf, frnd, fstat}
## @end deftypefn

function x = finv (p, df1, df2)

  if (nargin != 3)
    print_usage ();
  endif

  if (! isscalar (p) || ! isscalar (df1) || ! isscalar (df2))
    [retval, p, df1, df2] = common_size (p, df1, df2);
    if (retval > 0)
      error ("finv: P, DF1, and DF2 must be of common size or scalars.");
    endif
  endif

  if (iscomplex (p) || iscomplex (df1) || iscomplex (df2))
    error ("finv: P, DF1, and DF2 must not be complex.");
  endif

  if (isa (p, "single") || isa (df1, "single") || isa (df2, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  k = (p == 1) & (df1 > 0) & (df1 < Inf) & (df2 > 0) & (df2 < Inf);
  x(k) = Inf;

  k = (p >= 0) & (p < 1) & (df1 > 0) & (df1 < Inf) & (df2 > 0) & (df2 < Inf);
  if (isscalar (df1) && isscalar (df2))
    x(k) = ((1 ./ betainv (1 - p(k), df2/2, df1/2) - 1) * df2 / df1);
  else
    x(k) = ((1 ./ betainv (1 - p(k), df2(k)/2, df1(k)/2) - 1)
              .* df2(k) ./ df1(k));
  endif

endfunction


%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (finv (p, 2*ones (1,5), 2*ones (1,5)), [NaN 0 1 Inf NaN])
%!assert (finv (p, 2, 2*ones (1,5)), [NaN 0 1 Inf NaN])
%!assert (finv (p, 2*ones (1,5), 2), [NaN 0 1 Inf NaN])
%!assert (finv (p, [2 -Inf NaN Inf 2], 2), [NaN NaN NaN NaN NaN])
%!assert (finv (p, 2, [2 -Inf NaN Inf 2]), [NaN NaN NaN NaN NaN])
%!assert (finv ([p(1:2) NaN p(4:5)], 2, 2), [NaN 0 NaN Inf NaN])

## Test class of input preserved
%!assert (finv ([p, NaN], 2, 2), [NaN 0 1 Inf NaN NaN])
%!assert (finv (single ([p, NaN]), 2, 2), single ([NaN 0 1 Inf NaN NaN]))
%!assert (finv ([p, NaN], single (2), 2), single ([NaN 0 1 Inf NaN NaN]))
%!assert (finv ([p, NaN], 2, single (2)), single ([NaN 0 1 Inf NaN NaN]))

## Test input validation
%!error finv ()
%!error finv (1)
%!error finv (1,2)
%!error finv (1,2,3,4)
%!error finv (ones (3), ones (2), ones (2))
%!error finv (ones (2), ones (3), ones (2))
%!error finv (ones (2), ones (2), ones (3))
%!error finv (i, 2, 2)
%!error finv (2, i, 2)
%!error finv (2, 2, i)
