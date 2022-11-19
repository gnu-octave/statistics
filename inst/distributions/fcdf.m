## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{p} = tcdf (@var{x}, @var{v1}, @var{v2})
## @deftypefn {Function File} @var{p} = tcdf (@var{x}, @var{v1}, @var{v2}, "upper")
##
## F cumulative distribution function.
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the F distribution with @var{v1} and @var{v2} degrees of
## freedom.
##
## The size of @var{p} is the common size of @var{x}, @var{v1}, and @var{v2}. A
## scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## @code{@var{p} = fcdf (@var{x}, @var{v1}, @var{v2}, "upper")} computes the
## upper tail probability of the F distribution with @var{v1} and @var{v2}
## degrees of freedom at the values in @var{x}.
##
## @seealso{finv, fpdf, frnd, fstat}
## @end deftypefn

function p = fcdf (x, v1, v2, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("fcdf: too few input arguments.");
  endif

  ## Check for 'upper' flag
  if (nargin > 3 && strcmpi (uflag, "upper"))
    notnan = ! isnan (x);
    x(notnan) = 1 ./ max (0, x(notnan));
    tmp=v1;
    v1=v2;
    v2=tmp;
  elseif (nargin > 3  && ! strcmpi (uflag, "upper"))
    error ("fcdf: invalid argument for upper tail.");
  endif

  ## Check for common size of X, V1 and V2
  if (! isscalar (x) || ! isscalar (v1) || ! isscalar (v2))
    [err, x, v1, v2] = common_size (x, v1, v2);
    if (err > 0)
      error ("fcdf: X, V1, and V2 must be of common size or scalars.");
    endif
  endif

  ## Check for X, V1 and V2 being reals
  if (iscomplex (x) || iscomplex (v1) || iscomplex (v2))
    error ("tcdf: X, V1, and V2 must not be complex.");
  endif

  ## Initialize P according to appropriate class type
  if (isa (x, "single") || isa (v1, "single") || isa (v2, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Check X for NaNs while DFs <= 0 and make P = NaNs
  make_nan = (v1 <= 0 | v2 <= 0 | isnan(x) | isnan(v1) | isnan(v2));
  p(make_nan) = NaN;
  ## Check remaining valid X for Inf values and make P = 1
  is_inf = (x == Inf) & ! make_nan;
  if any (is_inf(:))
    p(is_inf) = 1;
    make_nan = (make_nan | is_inf);
  endif

  ## Compute P when X > 0.
  k = find(x > 0 & ! make_nan & isfinite(v1) & isfinite(v2));
  if (any (k))
    k1 = (v2(k) <= x(k) .* v1(k));
    if (any (k1))
      kk = k(k1);
      xx = v2(kk) ./ (v2(kk) + x(kk) .* v1(kk));
      p(kk) = betainc (xx, v2(kk)/2, v1(kk)/2, "upper");
    end
    if (any (! k1))
      kk = k(! k1);
      num = v1(kk) .* x(kk);
      xx = num ./ (num + v2(kk));
      p(kk) = betainc (xx, v1(kk)/2, v2(kk)/2, "lower");
    endif
  endif

  if any(~isfinite(v1(:)) | ~isfinite(v2(:)))
    k = find (x > 0 & ! make_nan & isfinite (v1) & ! isfinite (v2) & v2 > 0);
    if (any (k))
      p(k) = gammainc (v1(k) .* x(k) ./ 2, v1(k) ./ 2, 'lower');
    end
    k = find (x > 0 & ! make_nan & ! isfinite (v1) & v1 > 0 & isfinite (v2));
    if (any (k))
      p(k) = gammainc (v2(k) ./ x(k) ./ 2, v2(k) ./ 2, 'upper');
    end
    k = find (x > 0 & ! make_nan & ! isfinite (v1) & v1 > 0 & ...
                                   ! isfinite (v2) & v2 > 0);
    if (any (k))
      if (nargin >= 4 && x(k) == 1)
        p(k) = 0;
      else
        p(k) = (x(k)>=1);
      end
    endif
  endif

endfunction

%!shared x,y
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
%!error<fcdf: too few input arguments.> fcdf ()
%!error<fcdf: too few input arguments.> fcdf (1)
%!error<fcdf: too few input arguments.> fcdf (1, 2)
%!error<fcdf: invalid argument for upper tail.> fcdf (1, 2, 3, 4)
%!error<fcdf: invalid argument for upper tail.> fcdf (1, 2, 3, "tail")
%!error fcdf (ones (3), ones (2), ones (2))
%!error fcdf (ones (2), ones (3), ones (2))
%!error fcdf (ones (2), ones (2), ones (3))
%!error fcdf (i, 2, 2)
%!error fcdf (2, i, 2)
%!error fcdf (2, 2, i)
