## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 2013-2017 Julien Bect
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
## @deftypefn {Function File} @var{p} = tcdf (@var{x}, @var{v})
## @deftypefnx {Function File} @var{p} = tcdf (@var{x}, @var{v}, "upper")
##
## Student's T cumulative distribution function (cdf).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the t (Student) distribution with @var{v} degrees of
## freedom.
##
## The size of @var{p} is the common size of @var{x} and @var{v}. A scalar input
## functions as a constant matrix of the same size as the other input.
##
## @code{@var{p} = tcdf (@var{x}, @var{v}, "upper")} computes the upper tail
## probability of the Student's T distribution with @var{v} degrees of freedom,
## at the values in @var{x}.
##
## @seealso{tinv, tpdf, trnd, tstat}
## @end deftypefn

function p = tcdf (x, v, uflag)

  ## Check for valid number of input arguments
  if (nargin < 2 || nargin > 3)
    error ("tcdf: invalid number of input arguments.");
  endif

  ## Check for 'upper' flag
  if (nargin > 2 && strcmpi (uflag, "upper"))
    x = -x;
  elseif (nargin > 2  && ! strcmpi (uflag, "upper"))
    error ("tcdf: invalid argument for upper tail.");
  endif

  ## Check for common size of x and v
  if (! isscalar (x) || ! isscalar (v))
    [err, x, v] = common_size (x, v);
    if (err > 0)
      error ("tcdf: X and V must be of common size or scalars.");
    endif
  endif

  ## Check for X and V being reals
  if (iscomplex (x) || iscomplex (v))
    error ("tcdf: X and V must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (v, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Check for NaNs or df <= 0
  is_nan = isnan (x) | ! (v > 0);
  p(is_nan) = NaN;

  ## Check for Inf where df > 0 (for -Inf is already 0)
  k = (x == Inf) & (v > 0);
  p(k) = 1;

  ## Find finite values in X where df > 0
  k = isfinite (x) & (v > 0);

  ## Process more efficiently small positive integer DF up to 1e4
  ks = k & (fix (v) == v) & (v <= 1e4);
  if any (ks(:))
    if (isscalar (v))
      p(ks) = tcdf_integer_df (x(ks), v);
      return
    else
      vu = unique (v(ks));
      for i = 1:numel (vu)
        ki = k & (v == vu(i));
        p(ki) = tcdf_integer_df (x(ki), vu(i));
      endfor
    endif
  endif

  ## Proccess remaining values for DF (non-integers and > 1e4)
  k &= ! ks;

  ## Distinguish between small and big abs(x)
  xx = x .^ 2;
  x_big_abs = (xx > v);

  ## Deal with the case "abs(x) big"
  kk = k & x_big_abs;
  if (isscalar (v))
    p(kk) = betainc (v ./ (v + xx(kk)), v/2, 1/2) / 2;
  else
    p(kk) = betainc (v(kk) ./ (v(kk) + xx(kk)), v(kk)/2, 1/2) / 2;
  endif

  ## Deal with the case "abs(x) small"
  kk = k & ! x_big_abs;
  if (isscalar (v))
    p(kk) = 0.5 * (1 - betainc (xx(kk) ./ (v + xx(kk)), 1/2, v/2, "lower"));
  else
    p(kk) = 0.5 * (1 - betainc (xx(kk) ./ (v(kk) + xx(kk)), ...
                                  1/2, v(kk)/2, "upper"));
  endif

  ## For x > 0, F(x) = 1 - F(-|x|).
  k &= (x > 0);
  if (any (k(:)))
    p(k) = 1 - p(k);
  endif

  ## Special case for Cauchy distribution
  ## Use acot(-x) instead of the usual (atan x)/pi + 0.5 to avoid roundoff error
  xpos = (x > 0);
  c = (v == 1);
  p(c) = xpos(c) + acot (-x(c)) / pi;

  ## Make the result exact for the median
  p(x == 0 & ! is_nan) = 0.5;

endfunction

## Compute the t distribution CDF efficiently (without calling betainc)
## for small positive integer DF up to 1e4
function p = tcdf_integer_df (x, v)

  if (v == 1)
    p = 0.5 + atan(x)/pi;
  elseif (v == 2)
    p = 0.5 + x ./ (2 * sqrt(2 + x .^ 2));
  else
    xs = x ./ sqrt(v);
    xxf = 1 ./ (1 + xs .^ 2);
    u = s = 1;
    if mod (v, 2)  ## odd DF
      m = (v - 1) / 2;
      for i = 2:m
        u .*= (1 - 1/(2*i - 1)) .* xxf;
        s += u;
      endfor
      p = 0.5 + (xs .* xxf .* s + atan(xs)) / pi;
    else           ## even DF
      m = v / 2;
      for i = 1:(m - 1)
        u .*= (1 - 1/(2*i)) .* xxf;
        s += u;
      endfor
      p = 0.5 + (xs .* sqrt(xxf) .* s) / 2;
    endif
  endif
endfunction

%!shared x,y
%! x = [-Inf 0 1 Inf];
%! y = [0 1/2 3/4 1];
%!assert (tcdf (x, ones (1,4)), y, eps)
%!assert (tcdf (x, 1), y, eps)
%!assert (tcdf (x, [0 1 NaN 1]), [NaN 1/2 NaN 1], eps)
%!assert (tcdf ([x(1:2) NaN x(4)], 1), [y(1:2) NaN y(4)], eps)
%!assert (tcdf (2, 3, "upper"), 0.0697, 1e-4)
%!assert (tcdf (205, 5, "upper"), 2.6206e-11, 1e-14)

## Test class of input preserved
%!assert (tcdf ([x, NaN], 1), [y, NaN], eps)
%!assert (tcdf (single ([x, NaN]), 1), single ([y, NaN]), eps ("single"))
%!assert (tcdf ([x, NaN], single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<tcdf: invalid number of input arguments.> tcdf ()
%!error<tcdf: invalid number of input arguments.> tcdf (1)
%!error tcdf (1,2,3,4)
%!error<tcdf: invalid argument for upper tail.> tcdf (1, 2, "uper")
%!error<tcdf: invalid argument for upper tail.> tcdf (1, 2, 3)
%!error<tcdf: X and V must be of common size or scalars.> ...
%! tcdf (ones (3), ones (2))
%!error<tcdf: X and V must be of common size or scalars.> ...
%! tcdf (ones (3), ones (2))
%!error<tcdf: X and V must be of common size or scalars.> ...
%! tcdf (ones (3), ones (2), "upper")
%!error<tcdf: X and V must not be complex.> tcdf (i, 2)
%!error<tcdf: X and V must not be complex.> tcdf (2, i)

## Check some reference values

%!shared tol_rel
%! tol_rel = 10 * eps;

## check accuracy for small positive values
%!assert (tcdf (10^(-10), 2.5), 0.50000000003618087, -tol_rel)
%!assert (tcdf (10^(-11), 2.5), 0.50000000000361809, -tol_rel)
%!assert (tcdf (10^(-12), 2.5), 0.50000000000036181, -tol_rel)
%!assert (tcdf (10^(-13), 2.5), 0.50000000000003618, -tol_rel)
%!assert (tcdf (10^(-14), 2.5), 0.50000000000000362, -tol_rel)
%!assert (tcdf (10^(-15), 2.5), 0.50000000000000036, -tol_rel)
%!assert (tcdf (10^(-16), 2.5), 0.50000000000000004, -tol_rel)

## check accuracy for large negative values
%!assert (tcdf (-10^1, 2.5), 2.2207478836537124e-03, -tol_rel)
%!assert (tcdf (-10^2, 2.5), 7.1916492116661878e-06, -tol_rel)
%!assert (tcdf (-10^3, 2.5), 2.2747463948307452e-08, -tol_rel)
%!assert (tcdf (-10^4, 2.5), 7.1933970159922115e-11, -tol_rel)
%!assert (tcdf (-10^5, 2.5), 2.2747519231756221e-13, -tol_rel)

## # Reference values obtained using Python 2.7.4 and mpmath 0.17
##
## from mpmath import *
##
## mp.dps = 100
##
## def F(x_in, nu_in):
##     x = mpf(x_in);
##     nu = mpf(nu_in);
##     t = nu / (nu + x*x)
##     a = nu / 2
##     b = mpf(0.5)
##     F = betainc(a, b, 0, t, regularized=True) / 2
##     if (x > 0):
##         F = 1 - F
##     return F
##
## nu = 2.5
##
## for i in range(1, 6):
##     x = - power(mpf(10), mpf(i))
##     print "%%!assert (tcdf (-10^%d, 2.5), %s, -eps)" \
##         % (i, nstr(F(x, nu), 17))
##
## for i in range(10, 17):
##     x = power(mpf(10), -mpf(i))
##     print "%%!assert (tcdf (10^(-%d), 2.5), %s, -eps)" \
##         % (i, nstr(F(x, nu), 17))
