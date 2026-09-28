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
## @deftypefn  {statistics} {@var{p} =} stdrcdf (@var{x}, @var{k}, @var{df})
## @deftypefnx {statistics} {@var{p} =} stdrcdf (@var{x}, @var{k}, @var{df}, @qcode{'upper'})
##
## Studentized range cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the studentized range distribution for @var{k} groups and @var{df}
## degrees of freedom.  The size of @var{p} is the common size of @var{x},
## @var{k} and @var{df}.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## The studentized range is the range of @var{k} independent standard normal
## variables divided by an independent @math{sqrt (chi2 / df)}, where
## @math{chi2} has a chi-square distribution with @var{df} degrees of
## freedom.  It is the distribution behind Tukey's honestly significant
## difference and the Tukey-Kramer multiple comparison procedure.  @var{k}
## must be an integer of at least 2 and @var{df} positive, @code{Inf}
## included, which gives the range of @var{k} standard normals; otherwise
## @var{p} is @code{NaN}.
##
## @code{@var{p} = stdrcdf (@var{x}, @var{k}, @var{df}, "upper")} computes the
## upper tail probability of the studentized range distribution at the values
## in @var{x}.  Each tail is computed directly, so a small upper tail
## probability keeps its relative accuracy.
##
## The distribution is evaluated by numerical integration, accurate to about
## 1e-13 relative to the tail probability for @var{df} up to 1e4, and to
## about 1e-8 for @var{df} of 1e8.
##
## MATLAB has no public counterpart of this function.
##
## Further information about the studentized range distribution can be found
## at @url{https://en.wikipedia.org/wiki/Studentized_range_distribution}
##
## Input arguments must be @qcode{double} or @qcode{single}; integer, logical,
## and character arrays are rejected.
##
## @seealso{stdrinv, stdrpdf, stdrrnd, stdrstat, multcompare}
## @end deftypefn

function p = stdrcdf (x, k, df, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("stdrcdf: function called with too few input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 3 && strcmpi (uflag, 'upper'))
    utail = true;
  elseif (nargin > 3 && ! strcmpi (uflag, 'upper'))
    error ("stdrcdf: invalid argument for upper tail.");
  else
    utail = false;
  endif

  ## Check for common size of X, K, and DF
  if (! isscalar (x) || ! isscalar (k) || ! isscalar (df))
    [err, x, k, df] = common_size (x, k, df);
    if (err > 0)
      error ("stdrcdf: X, K, and DF must be of common size or scalars.");
    endif
  endif

  ## Check for X, K, and DF being double or single
  if (! (isfloat (x) && isfloat (k) && isfloat (df)))
    error ("stdrcdf: X, K, and DF must be double or single.");
  endif

  ## Check for X, K, and DF being reals
  if (iscomplex (x) || iscomplex (k) || iscomplex (df))
    error ("stdrcdf: X, K, and DF must not be complex.");
  endif

  ## Check for class type
  if (isa (x, 'single') || isa (k, 'single') || isa (df, 'single'))
    cls = 'single';
  else
    cls = 'double';
  endif
  sz = size (x);
  x = double (x(:));
  k = double (k(:)) .* ones (numel (x), 1);
  df = double (df(:)) .* ones (numel (x), 1);

  ## Invalid parameters or X give NaN
  valid = (k >= 2) & (k == fix (k)) & isfinite (k) & (df > 0);
  lower = zeros (numel (x), 1);
  lower(! valid | isnan (x)) = NaN;
  lower(valid & x == Inf) = 1;
  uppr = 1 - lower;

  ## Evaluate the positive finite values, one parameter pair at a time
  todo = valid & (x > 0) & isfinite (x);
  if (any (todo))
    pairs = unique ([k(todo), df(todo)], 'rows');
    for i = 1:rows (pairs)
      idx = todo & (k == pairs(i,1)) & (df == pairs(i,2));
      [F, U] = __stdr__ (x(idx), pairs(i,1), pairs(i,2), 'cdf');
      lower(idx) = F;
      uppr(idx) = U;
    endfor
  endif

  if (utail)
    p = uppr;
  else
    p = lower;
  endif
  p = cast (reshape (p, sz), cls);

endfunction

%!demo
%! ## Plot various CDFs from the studentized range distribution
%! x = 0:0.01:8;
%! p1 = stdrcdf (x, 2, 5);
%! p2 = stdrcdf (x, 3, 5);
%! p3 = stdrcdf (x, 5, 5);
%! p4 = stdrcdf (x, 5, Inf);
%! plot (x, p1, '-b', x, p2, '-g', x, p3, '-r', x, p4, '-m')
%! grid on
%! ylim ([0, 1])
%! legend ({'k = 2, df = 5', 'k = 3, df = 5', 'k = 5, df = 5', ...
%!          'k = 5, df = \infty'}, 'location', 'southeast')
%! title ('Studentized range CDF')
%! xlabel ('values in x')
%! ylabel ('probability')

## Two groups are sqrt (2) times the absolute value of a t variable
%!shared x
%! x = [0.5, 2, 4, 8];
%!assert_equal (stdrcdf (x, 2, 5), 1 - 2 * tcdf (-x / sqrt (2), 5), -1e-13)
%!assert_equal (stdrcdf (x, 2, 5, 'upper'), 2 * tcdf (-x / sqrt (2), 5), -1e-13)
%!assert_equal (stdrcdf (x, 2, Inf), erf (x / 2), -1e-13)
%!assert_equal (stdrcdf (x, 2, Inf, 'upper'), erfc (x / 2), -1e-13)

## A far upper tail keeps its relative accuracy
%!assert_equal (stdrcdf (20, 2, 30, 'upper'), ...
%!              betainc (200 / 230, 1/2, 15, 'upper'), -1e-12)

## Upper tail probabilities from R 4.5.0 ptukey, accurate to about 1e-6 for
## DF of 10 or more
%!assert_equal (stdrcdf ([2, 4, 8], 3, 10, 'upper'), ...
%!              [0.37054467503555799, 0.043349349760161804, ...
%!               0.00055885879221251322], -1e-5)
%!assert_equal (stdrcdf ([2, 4, 8], 10, 30, 'upper'), ...
%!              [0.91315280898350459, 0.17214641499562477, ...
%!               0.00013979917001072373], -1e-5)
%!assert_equal (stdrcdf ([2, 4, 8], 20, Inf, 'upper'), ...
%!              [0.9976642515512063, 0.3360234394147692, ...
%!               2.8901613927656555e-06], -1e-5)

## Upper tail probabilities from MATLAB R2024a, through the Tukey-Kramer
## p-values of multcompare, accurate to about 1e-3 above 1e-3
%!assert_equal (stdrcdf ([0.5, 1, 1.5], 3, 10, 'upper'), ...
%!              [0.93386416513356696, 0.76489080347916205, ...
%!               0.55804816568706039], -1e-3)
%!assert_equal (stdrcdf ([0.5, 1.5, 3, 5], 5, 30, 'upper'), ...
%!              [0.99646380721257077, 0.82481287273578452, ...
%!               0.23767360358955125, 0.010890009004856521], -1e-3)

## Edge values and invalid parameters
%!assert_equal (stdrcdf ([-1, 0, Inf, NaN], 3, 10), [0, 0, 1, NaN])
%!assert_equal (stdrcdf ([-1, 0, Inf, NaN], 3, 10, 'upper'), [1, 1, 0, NaN])
%!assert_equal (stdrcdf (2, [1, 2.5, NaN, Inf], 10), NaN (1, 4))
%!assert_equal (stdrcdf (2, 3, [0, -1, NaN]), NaN (1, 3))

## Test class of input preserved
%!assert_equal (class (stdrcdf (single (2), 3, 10)), 'single')
%!assert_equal (class (stdrcdf (2, single (3), 10)), 'single')
%!assert_equal (class (stdrcdf (2, 3, single (10))), 'single')

## Test input validation
%!error<stdrcdf: function called with too few input arguments.> stdrcdf ()
%!error<stdrcdf: function called with too few input arguments.> stdrcdf (1, 2)
%!error<stdrcdf: invalid argument for upper tail.> stdrcdf (1, 2, 3, 'uper')
%!error<stdrcdf: invalid argument for upper tail.> stdrcdf (1, 2, 3, 4)
%!error<stdrcdf: X, K, and DF must be of common size or scalars.> ...
%! stdrcdf (ones (3), ones (2), 3)
%!error<stdrcdf: X, K, and DF must be double or single.> ...
%! stdrcdf (int32 (2), 3, 10)
%!error<stdrcdf: X, K, and DF must be double or single.> ...
%! stdrcdf (true, 3, 10)
%!error<stdrcdf: X, K, and DF must be double or single.> ...
%! stdrcdf ('a', 3, 10)
%!error<stdrcdf: X, K, and DF must not be complex.> stdrcdf (i, 3, 10)
%!error<stdrcdf: X, K, and DF must not be complex.> stdrcdf (2, i, 10)
%!error<stdrcdf: X, K, and DF must not be complex.> stdrcdf (2, 3, i)
