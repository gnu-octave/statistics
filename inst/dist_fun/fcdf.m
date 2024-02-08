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
## @deftypefn  {statistics} {@var{p} =} fcdf (@var{x}, @var{df1}, @var{df2})
## @deftypefnx {statistics} {@var{p} =} fcdf (@var{x}, @var{df1}, @var{df2}, @qcode{"upper"})
##
## @math{F}-cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the @math{F}-distribution with @var{df1} and @var{df2} degrees of
## freedom.  The size of @var{p} is the common size of @var{x}, @var{df1}, and
## @var{df2}.  A scalar input functions as a constant matrix of the same size as
## the other inputs.
##
## @code{@var{p} = fcdf (@var{x}, @var{df1}, @var{df2}, "upper")} computes the
## upper tail probability of the @math{F}-distribution with @var{df1} and
## @var{df2} degrees of freedom, at the values in @var{x}.
##
## Further information about the @math{F}-distribution can be found at
## @url{https://en.wikipedia.org/wiki/F-distribution}
##
## @seealso{finv, fpdf, frnd, fstat}
## @end deftypefn

function p = fcdf (x, df1, df2, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("fcdf: function called with too few input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 3 && strcmpi (uflag, "upper"))
    notnan = ! isnan (x);
    x(notnan) = 1 ./ max (0, x(notnan));
    tmp=df1;
    df1=df2;
    df2=tmp;
  elseif (nargin > 3  && ! strcmpi (uflag, "upper"))
    error ("fcdf: invalid argument for upper tail.");
  endif

  ## Check for common size of X, DF1, and DF2
  if (! isscalar (x) || ! isscalar (df1) || ! isscalar (df2))
    [err, x, df1, df2] = common_size (x, df1, df2);
    if (err > 0)
      error ("fcdf: X, DF1, and DF2 must be of common size or scalars.");
    endif
  endif

  ## Check for X, DF1, and DF2 being reals
  if (iscomplex (x) || iscomplex (df1) || iscomplex (df2))
    error ("fcdf: X, DF1, and DF2 must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (df1, "single") || isa (df2, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Check X for NaNs while DFs <= 0 and make P = NaNs
  make_nan = (df1 <= 0 | df2 <= 0 | isnan(x) | isnan(df1) | isnan(df2));
  p(make_nan) = NaN;
  ## Check remaining valid X for Inf values and make P = 1
  is_inf = (x == Inf) & ! make_nan;
  if any (is_inf(:))
    p(is_inf) = 1;
    make_nan = (make_nan | is_inf);
  endif

  ## Compute P when X > 0.
  k = find(x > 0 & ! make_nan & isfinite(df1) & isfinite(df2));
  if (any (k))
    k1 = (df2(k) <= x(k) .* df1(k));
    if (any (k1))
      kk = k(k1);
      xx = df2(kk) ./ (df2(kk) + x(kk) .* df1(kk));
      p(kk) = betainc (xx, df2(kk)/2, df1(kk)/2, "upper");
    end
    if (any (! k1))
      kk = k(! k1);
      num = df1(kk) .* x(kk);
      xx = num ./ (num + df2(kk));
      p(kk) = betainc (xx, df1(kk)/2, df2(kk)/2, "lower");
    endif
  endif

  if any(~isfinite(df1(:)) | ~isfinite(df2(:)))
    k = find (x > 0 & ! make_nan & isfinite (df1) & ! isfinite (df2) & df2 > 0);
    if (any (k))
      p(k) = gammainc (df1(k) .* x(k) ./ 2, df1(k) ./ 2, "lower");
    end
    k = find (x > 0 & ! make_nan & ! isfinite (df1) & df1 > 0 & isfinite (df2));
    if (any (k))
      p(k) = gammainc (df2(k) ./ x(k) ./ 2, df2(k) ./ 2, "upper");
    end
    k = find (x > 0 & ! make_nan & ! isfinite (df1) & df1 > 0 & ...
                                   ! isfinite (df2) & df2 > 0);
    if (any (k))
      if (nargin >= 4 && x(k) == 1)
        p(k) = 0;
      else
        p(k) = (x(k)>=1);
      end
    endif
  endif

endfunction

%!demo
%! ## Plot various CDFs from the F distribution
%! x = 0.01:0.01:4;
%! p1 = fcdf (x, 1, 2);
%! p2 = fcdf (x, 2, 1);
%! p3 = fcdf (x, 5, 2);
%! p4 = fcdf (x, 10, 1);
%! p5 = fcdf (x, 100, 100);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c", x, p5, "-m")
%! grid on
%! legend ({"df1 = 1, df2 = 2", "df1 = 2, df2 = 1", ...
%!          "df1 = 5, df2 = 2", "df1 = 10, df2 = 1", ...
%!          "df1 = 100, df2 = 100"}, "location", "southeast")
%! title ("F CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-1, 0, 0.5, 1, 2, Inf];
%! y = [0, 0, 1/3, 1/2, 2/3, 1];
%!assert (fcdf (x, 2*ones (1,6), 2*ones (1,6)), y, eps)
%!assert (fcdf (x, 2, 2*ones (1,6)), y, eps)
%!assert (fcdf (x, 2*ones (1,6), 2), y, eps)
%!assert (fcdf (x, [0 NaN Inf 2 2 2], 2), [NaN NaN 0.1353352832366127 y(4:6)], eps)
%!assert (fcdf (x, 2, [0 NaN Inf 2 2 2]), [NaN NaN 0.3934693402873666 y(4:6)], eps)
%!assert (fcdf ([x(1:2) NaN x(4:6)], 2, 2), [y(1:2) NaN y(4:6)], eps)

## Test class of input preserved
%!assert (fcdf ([x, NaN], 2, 2), [y, NaN], eps)
%!assert (fcdf (single ([x, NaN]), 2, 2), single ([y, NaN]), eps ("single"))
%!assert (fcdf ([x, NaN], single (2), 2), single ([y, NaN]), eps ("single"))
%!assert (fcdf ([x, NaN], 2, single (2)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<fcdf: function called with too few input arguments.> fcdf ()
%!error<fcdf: function called with too few input arguments.> fcdf (1)
%!error<fcdf: function called with too few input arguments.> fcdf (1, 2)
%!error<fcdf: invalid argument for upper tail.> fcdf (1, 2, 3, 4)
%!error<fcdf: invalid argument for upper tail.> fcdf (1, 2, 3, "tail")
%!error<fcdf: X, DF1, and DF2 must be of common size or scalars.> ...
%! fcdf (ones (3), ones (2), ones (2))
%!error<fcdf: X, DF1, and DF2 must be of common size or scalars.> ...
%! fcdf (ones (2), ones (3), ones (2))
%!error<fcdf: X, DF1, and DF2 must be of common size or scalars.> ...
%! fcdf (ones (2), ones (2), ones (3))
%!error<fcdf: X, DF1, and DF2 must not be complex.> fcdf (i, 2, 2)
%!error<fcdf: X, DF1, and DF2 must not be complex.> fcdf (2, i, 2)
%!error<fcdf: X, DF1, and DF2 must not be complex.> fcdf (2, 2, i)
