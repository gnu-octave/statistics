## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {statistics} {@var{x} =} stdrinv (@var{p}, @var{k}, @var{df})
##
## Inverse of the studentized range cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## of the studentized range distribution for @var{k} groups and @var{df}
## degrees of freedom.  The size of @var{x} is the common size of @var{p},
## @var{k} and @var{df}.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## @code{stdrinv (1 - @var{alpha}, @var{k}, @var{df})} is the critical value
## of Tukey's honestly significant difference at level @var{alpha}.  @var{k}
## must be an integer of at least 2 and @var{df} positive, @code{Inf}
## included; otherwise @var{x} is @code{NaN}, as it is for @var{p} outside
## @math{[0, 1]}.  The quantile is found to about 1e-13 relative accuracy.
##
## MATLAB has no public counterpart of this function.
##
## Further information about the studentized range distribution can be found
## at @url{https://en.wikipedia.org/wiki/Studentized_range_distribution}
##
## Input arguments must be @qcode{double} or @qcode{single}; integer, logical,
## and character arrays are rejected.
##
## @seealso{stdrcdf, stdrpdf, stdrrnd, stdrstat, multcompare}
## @end deftypefn

function x = stdrinv (p, k, df)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("stdrinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, K, and DF
  if (! isscalar (p) || ! isscalar (k) || ! isscalar (df))
    [err, p, k, df] = common_size (p, k, df);
    if (err > 0)
      error ("stdrinv: P, K, and DF must be of common size or scalars.");
    endif
  endif

  ## Check for P, K, and DF being double or single
  if (! (isfloat (p) && isfloat (k) && isfloat (df)))
    error ("stdrinv: P, K, and DF must be double or single.");
  endif

  ## Check for P, K, and DF being reals
  if (iscomplex (p) || iscomplex (k) || iscomplex (df))
    error ("stdrinv: P, K, and DF must not be complex.");
  endif

  ## Check for class type
  if (isa (p, 'single') || isa (k, 'single') || isa (df, 'single'))
    cls = 'single';
  else
    cls = 'double';
  endif
  sz = size (p);
  p = double (p(:));
  k = double (k(:)) .* ones (numel (p), 1);
  df = double (df(:)) .* ones (numel (p), 1);

  valid = (k >= 2) & (k == fix (k)) & isfinite (k) & (df > 0);
  x = NaN (numel (p), 1);
  x(valid & p == 0) = 0;
  x(valid & p == 1) = Inf;

  for i = find (valid & p > 0 & p < 1)'
    x(i) = invert (p(i), k(i), df(i));
  endfor
  x = cast (reshape (x, sz), cls);

endfunction

## Newton's method on u = log (x), matching the log of whichever tail P
## leaves smaller so that either extreme keeps its digits, inside a bracket
## that is halved whenever a step would leave it.
function q = invert (p, k, df)

  lower = (p <= 0.5);
  if (lower)
    target = log (p);
  else
    target = log1p (- p);
  endif
  ## Start from the Bonferroni bound over the k (k - 1) / 2 pairs
  u = log (sqrt (2) * tinv (1 - min (p, 1 - p) / (k * (k - 1)), df));
  if (lower)
    u = log (sqrt (2) * tinv ((1 + p) / 2, df));
  endif
  [h, dh] = gap (u, k, df, lower, target);
  a = -Inf;
  b = Inf;
  for iter = 1:100
    if (h < 0)
      a = u;
    else
      b = u;
    endif
    step = - h / dh;
    un = u + step;
    if (! isfinite (un) || un <= a || un >= b)
      if (isinf (a))
        un = b - 2;
      elseif (isinf (b))
        un = a + 2;
      else
        un = (a + b) / 2;
      endif
    endif
    if (abs (un - u) <= 1e-14 * max (1, abs (u)))
      u = un;
      break;
    endif
    u = un;
    [h, dh] = gap (u, k, df, lower, target);
    if (h == 0)
      break;
    endif
  endfor
  q = exp (u);

endfunction

## The log tail at exp (U) less its target, increasing in U, and its
## derivative in U.
function [h, dh] = gap (u, k, df, lower, target)

  q = exp (u);
  f = stdrpdf (q, k, df);
  if (lower)
    T = stdrcdf (q, k, df);
    h = log (T) - target;
    dh = f * q / T;
  else
    T = stdrcdf (q, k, df, 'upper');
    h = target - log (T);
    dh = f * q / T;
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the studentized range distribution
%! p = 0.01:0.01:0.99;
%! x1 = stdrinv (p, 2, 5);
%! x2 = stdrinv (p, 3, 5);
%! x3 = stdrinv (p, 5, Inf);
%! plot (p, x1, '-b', p, x2, '-g', p, x3, '-r')
%! grid on
%! legend ({'k = 2, df = 5', 'k = 3, df = 5', 'k = 5, df = \infty'}, ...
%!         'location', 'northwest')
%! title ('Studentized range iCDF')
%! xlabel ('probability')
%! ylabel ('values in x')

## Two groups are sqrt (2) times the absolute value of a t variable
%!shared p
%! p = [0.1, 0.5, 0.95, 0.999];
%!assert_equal (stdrinv (p, 2, 5), sqrt (2) * tinv ((1 + p) / 2, 5), -1e-12)
%!assert_equal (stdrinv (p, 2, 1), sqrt (2) * tan (pi * p / 2), -1e-12)

## Quantiles from R 4.5.0 qtukey, accurate to about 1e-7 for DF of 10 or more
%!assert_equal (stdrinv ([0.5, 0.9, 0.95, 0.99], 3, 10), ...
%!              [1.6446889006146324, 3.2703084031559371, ...
%!               3.8767767491915595, 5.2701615370332799], -1e-6)
%!assert_equal (stdrinv (0.95, 20, Inf), 5.0116887946867648, -1e-6)

## Quantiles from MATLAB R2024a, through the Tukey-Kramer limits of
## multcompare
%!assert_equal (stdrinv (0.95, 5, 30), 4.1020790196264514, -1e-4)
%!assert_equal (stdrinv (0.95, 3, 10), 3.8767767552295549, -1e-4)

## The quantile returns the probability it was found for, in either tail
%!assert_equal (stdrcdf (stdrinv (0.975, 4, 7), 4, 7), 0.975, -1e-13)
%!assert_equal (stdrcdf (stdrinv (1e-8, 4, 7), 4, 7), 1e-8, -1e-10)
%!assert_equal (stdrcdf (stdrinv (1 - 1e-10, 3, 20), 3, 20, 'upper'), ...
%!              1e-10, -1e-6)

## Edge values and invalid parameters
%!assert_equal (stdrinv ([0, 1, -1, 2, NaN], 3, 10), [0, Inf, NaN, NaN, NaN])
%!assert_equal (stdrinv (0.5, [1, 2.5, NaN, Inf], 10), NaN (1, 4))
%!assert_equal (stdrinv (0.5, 3, [0, -1, NaN]), NaN (1, 3))

## Test class of input preserved
%!assert_equal (class (stdrinv (single (0.5), 3, 10)), 'single')
%!assert_equal (class (stdrinv (0.5, 3, single (10))), 'single')

## Test input validation
%!error<stdrinv: function called with too few input arguments.> stdrinv ()
%!error<stdrinv: function called with too few input arguments.> stdrinv (1, 2)
%!error<stdrinv: P, K, and DF must be of common size or scalars.> ...
%! stdrinv (ones (3), ones (2), 3)
%!error<stdrinv: P, K, and DF must be double or single.> ...
%! stdrinv (int32 (1), 3, 10)
%!error<stdrinv: P, K, and DF must be double or single.> ...
%! stdrinv (true, 3, 10)
%!error<stdrinv: P, K, and DF must be double or single.> ...
%! stdrinv ('a', 3, 10)
%!error<stdrinv: P, K, and DF must not be complex.> stdrinv (i, 3, 10)
%!error<stdrinv: P, K, and DF must not be complex.> stdrinv (0.5, i, 10)
%!error<stdrinv: P, K, and DF must not be complex.> stdrinv (0.5, 3, i)
