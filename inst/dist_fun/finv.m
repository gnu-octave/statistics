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
## @deftypefn  {statistics} {@var{x} =} finv (@var{p}, @var{df1}, @var{df2})
##
## Inverse of the @math{F}-cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the @math{F}-distribution with @var{df1} and @var{df2} degrees of freedom.
## The size of @var{x} is the common size of @var{p}, @var{df1}, and @var{df2}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## Further information about the @math{F}-distribution can be found at
## @url{https://en.wikipedia.org/wiki/F-distribution}
##
## @seealso{fcdf, fpdf, frnd, fstat}
## @end deftypefn

function x = finv (p, df1, df2)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("finv: function called with too few input arguments.");
  endif

  ## Check for common size of P, DF1, and DF2
  if (! isscalar (p) || ! isscalar (df1) || ! isscalar (df2))
    [retval, p, df1, df2] = common_size (p, df1, df2);
    if (retval > 0)
      error ("finv: P, DF1, and DF2 must be of common size or scalars.");
    endif
  endif

  ## Check for P, DF1, and DF2 being reals
  if (iscomplex (p) || iscomplex (df1) || iscomplex (df2))
    error ("finv: P, DF1, and DF2 must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (df1, "single") || isa (df2, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  k = (p == 1) & (df1 > 0) & (df1 < Inf) & (df2 > 0) & (df2 < Inf);
  x(k) = Inf;

  ## Limit df2 to 1e6 unless it is Inf
  k = (df2 > 1e6) & (df2 < Inf);
  df2(k) = 1e6;

  k = (p >= 0) & (p < 1) & (df1 > 0) & (df1 < Inf) & (df2 > 0) & (df2 < Inf);
  if (isscalar (df1) && isscalar (df2))
    x(k) = ((1 ./ betainv (1 - p(k), df2/2, df1/2) - 1) * df2 / df1);
  else
    x(k) = ((1 ./ betainv (1 - p(k), df2(k)/2, df1(k)/2) - 1)
              .* df2(k) ./ df1(k));
  endif

  ## Handle case when DF2 is infinite
  k = (p >= 0) & (p < 1) & (df1 > 0) & (df1 < Inf) & (df2 == Inf);
  x(k) = chi2inv (p(k), df1(k)) ./ df1(k);

endfunction

%!demo
%! ## Plot various iCDFs from the F distribution
%! p = 0.001:0.001:0.999;
%! x1 = finv (p, 1, 1);
%! x2 = finv (p, 2, 1);
%! x3 = finv (p, 5, 2);
%! x4 = finv (p, 10, 1);
%! x5 = finv (p, 100, 100);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c", p, x5, "-m")
%! grid on
%! ylim ([0, 4])
%! legend ({"df1 = 1, df2 = 2", "df1 = 2, df2 = 1", ...
%!          "df1 = 5, df2 = 2", "df1 = 10, df2 = 1", ...
%!          "df1 = 100, df2 = 100"}, "location", "northwest")
%! title ("F iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (finv (p, 2*ones (1,5), 2*ones (1,5)), [NaN 0 1 Inf NaN])
%!assert (finv (p, 2, 2*ones (1,5)), [NaN 0 1 Inf NaN])
%!assert (finv (p, 2*ones (1,5), 2), [NaN 0 1 Inf NaN])
%!assert (finv (p, [2 -Inf NaN Inf 2], 2), [NaN NaN NaN NaN NaN])
%!assert (finv (p, 2, [2 -Inf NaN Inf 2]), [NaN NaN NaN NaN NaN])
%!assert (finv ([p(1:2) NaN p(4:5)], 2, 2), [NaN 0 NaN Inf NaN])

## Test for bug #66034 (savannah)
%!assert (finv (0.025, 10, 1e6), 0.3247, 1e-4)
%!assert (finv (0.025, 10, 1e7), 0.3247, 1e-4)
%!assert (finv (0.025, 10, 1e10), 0.3247, 1e-4)
%!assert (finv (0.025, 10, 1e255), 0.3247, 1e-4)
%!assert (finv (0.025, 10, Inf), 0.3247, 1e-4)

## Test class of input preserved
%!assert (finv ([p, NaN], 2, 2), [NaN 0 1 Inf NaN NaN])
%!assert (finv (single ([p, NaN]), 2, 2), single ([NaN 0 1 Inf NaN NaN]))
%!assert (finv ([p, NaN], single (2), 2), single ([NaN 0 1 Inf NaN NaN]))
%!assert (finv ([p, NaN], 2, single (2)), single ([NaN 0 1 Inf NaN NaN]))

## Test input validation
%!error<finv: function called with too few input arguments.> finv ()
%!error<finv: function called with too few input arguments.> finv (1)
%!error<finv: function called with too few input arguments.> finv (1,2)
%!error<finv: P, DF1, and DF2 must be of common size or scalars.> ...
%! finv (ones (3), ones (2), ones (2))
%!error<finv: P, DF1, and DF2 must be of common size or scalars.> ...
%! finv (ones (2), ones (3), ones (2))
%!error<finv: P, DF1, and DF2 must be of common size or scalars.> ...
%! finv (ones (2), ones (2), ones (3))
%!error<finv: P, DF1, and DF2 must not be complex.> finv (i, 2, 2)
%!error<finv: P, DF1, and DF2 must not be complex.> finv (2, i, 2)
%!error<finv: P, DF1, and DF2 must not be complex.> finv (2, 2, i)
