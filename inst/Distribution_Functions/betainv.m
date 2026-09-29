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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{x} =} betainv (@var{p}, @var{a}, @var{b})
##
## Inverse of the Beta distribution (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## of the Beta distribution with shape parameters @var{a} and @var{b}.  The size
## of @var{x} is the common size of @var{x}, @var{a}, and @var{b}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## Further information about the Beta distribution can be found at
## @url{https://en.wikipedia.org/wiki/Beta_distribution}
##
## Input arguments must be @qcode{double} or @qcode{single}; integer, logical,
## and character arrays are rejected.  MATLAB accepts a character array and
## evaluates it at the character codes, which Octave deliberately does not,
## since a character array is an integer type and integers are refused too.
##
## @seealso{betacdf, betapdf, betarnd, betafit, betalike, betastat}
## @end deftypefn

function x = betainv (p, a, b)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("betainv: function called with too few input arguments.");
  endif

  ## Check for common size of P, A, and B
  if (! isscalar (p) || ! isscalar (a) || ! isscalar (b))
    [retval, p, a, b] = common_size (p, a, b);
    if (retval > 0)
      error ("betainv: P, A, and B must be of common size or scalars.");
    endif
  endif

  ## Check for P, A, and B being double or single
  if (! (isfloat (p) && isfloat (a) && isfloat (b)))
    error ("betainv: P, A, and B must be double or single.");
  endif

  ## Check for P, A, and B being reals
  if (iscomplex (p) || iscomplex (a) || iscomplex (b))
    error ("betainv: P, A, and B must not be complex.");
  endif

  ## Check for class type
  if (isa (p, 'single') || isa (a, 'single') || isa (b, 'single'))
    x = zeros (size (p), 'single');
  else
    x = zeros (size (p));
  endif

  k = (p < 0) | (p > 1) | ! (a > 0) | ! (b > 0) | isnan (p);
  x(k) = NaN;

  k = (p == 1) & (a > 0) & (b > 0);
  x(k) = 1;

  k = find ((p > 0) & (p < 1) & (a > 0) & (b > 0));
  if (! isempty (k))
    if (isscalar (a))
      a = a * ones (size (k));
    else
      a = a(k);
    endif
    if (isscalar (b))
      b = b * ones (size (k));
    else
      b = b(k);
    endif
    x(k) = bInverse (double (p(k)(:)), double (a(:)), double (b(:)));
  endif

endfunction

## The quantiles, solved by Newton's method on betainc kept inside a bracket.
## Below the median the lower tail is matched to P, above it the upper tail
## to 1 - P, which is exact there, so neither tail loses its digits to a
## subtraction from 1.  The steps stop on a relative tolerance, so a
## quantile of 1e-20 is resolved as closely as one of 0.5.
function x = bInverse (p, a, b)

  up = p > 0.5;
  q = p;
  q(up) = 1 - p(up);
  lbeta = betaln (a, b);

  ## Start from the leading term of the tail: x^a / (a B) below the median,
  ## (1 - x)^b / (b B) above it, and from the mean where that falls outside
  x = exp ((log (q) + log (a) + lbeta) ./ a);
  x(up) = -expm1 ((log (q(up)) + log (b(up)) + lbeta(up)) ./ b(up));
  bad = ! (x > 0 & x < 1);
  x(bad) = a(bad) ./ (a(bad) + b(bad));

  lo = zeros (size (x));
  hi = ones (size (x));
  todo = true (size (x));
  for iter = 1:200
    i = find (todo);
    if (isempty (i))
      break;
    endif
    xi = x(i);
    ## The tail matched, signed so that it grows with x
    g = betainc (xi, a(i), b(i)) - q(i);
    u = up(i);
    g(u) = q(i)(u) - betainc (xi(u), a(i)(u), b(i)(u), 'upper');
    ## The root lies below where the tail is too large
    over = g > 0;
    hi(i(over)) = xi(over);
    lo(i(! over)) = xi(! over);
    ## Newton's step, and bisection where it leaves the bracket; geometric
    ## bisection where the bracket spans more than a factor of 4
    dens = exp ((a(i) - 1) .* log (xi) + (b(i) - 1) .* log1p (-xi) - lbeta(i));
    xn = xi - g ./ dens;
    xn(g == 0) = xi(g == 0);
    out = ! (xn >= lo(i) & xn <= hi(i));
    geo = out & lo(i) > 0 & hi(i) > 4 * lo(i);
    xn(geo) = sqrt (lo(i)(geo) .* hi(i)(geo));
    fromzero = out & lo(i) == 0;
    xn(fromzero) = hi(i)(fromzero) / 16;
    mid = out & ! geo & ! fromzero;
    xn(mid) = (lo(i)(mid) + hi(i)(mid)) / 2;
    x(i) = xn;
    ## Near a root betainc is noisy by some tens of units in the last place
    ## of its value, below which no step can improve the answer
    done = abs (g) <= 64 * eps (q(i)) | abs (xn - xi) <= 8 * eps (xn) ...
           | (hi(i) - lo(i)) <= 64 * eps (lo(i));
    todo(i(done)) = false;
  endfor

  ## Values still stepping have reached betainc's own noise; only a residual
  ## well above it means the solution was not found
  if (any (todo))
    j = find (todo);
    if (any (abs (bResidual (x(j), q(j), a(j), b(j), up(j))) > 1e-10 * q(j)))
      warning ("betainv: calculation failed to converge for some values.");
    endif
  endif

  ## Step to a neighbouring double for as long as it matches its tail better
  i = (1:numel (x))';
  for step = 1:64
    gx = abs (bResidual (x(i), q(i), a(i), b(i), up(i)));
    moved = false (size (i));
    for d = [-1, 1]
      y = x(i) + d * eps (x(i));
      ok = y >= 0 & y <= 1;
      y = min (max (y, 0), 1);
      gy = abs (bResidual (y, q(i), a(i), b(i), up(i)));
      better = ok & gy < gx;
      x(i(better)) = y(better);
      gx(better) = gy(better);
      moved |= better;
    endfor
    i = i(moved);
    if (isempty (i))
      break;
    endif
  endfor

endfunction

## How far the tail matched at X is from its target.
function g = bResidual (x, q, a, b, up)
  g = betainc (x, a, b) - q;
  g(up) = q(up) - betainc (x(up), a(up), b(up), 'upper');
endfunction


%!demo
%! ## Plot various iCDFs from the Beta distribution
%! p = 0.001:0.001:0.999;
%! x1 = betainv (p, 0.5, 0.5);
%! x2 = betainv (p, 5, 1);
%! x3 = betainv (p, 1, 3);
%! x4 = betainv (p, 2, 2);
%! x5 = betainv (p, 2, 5);
%! plot (p, x1, '-b', p, x2, '-g', p, x3, '-r', p, x4, '-c', p, x5, '-m')
%! grid on
%! legend ({'α = β = 0.5', 'α = 5, β = 1', 'α = 1, β = 3', ...
%!          'α = 2, β = 2', 'α = 2, β = 5'}, 'location', 'southeast')
%! title ('Beta iCDF')
%! xlabel ('probability')
%! ylabel ('values in x')

## Test output
%!shared p
%! p = [-1 0 0.75 1 2];
%!assert_equal (betainv (p, ones (1,5), 2*ones (1,5)), [NaN 0 0.5 1 NaN], eps)
%!assert_equal (betainv (p, 1, 2*ones (1,5)), [NaN 0 0.5 1 NaN], eps)
%!assert_equal (betainv (p, ones (1,5), 2), [NaN 0 0.5 1 NaN], eps)
%!assert_equal (betainv (p, [1 0 NaN 1 1], 2), [NaN NaN NaN 1 NaN])
%!assert_equal (betainv (p, 1, 2*[1 0 NaN 1 1]), [NaN NaN NaN 1 NaN])
%!assert_equal (betainv ([p(1:2) NaN p(4:5)], 1, 2), [NaN 0 NaN 1 NaN])

## Closed forms, deep into both tails
%!test
%! pp = [1e-300, 1e-20, 0.025, 0.5, 1 - 1e-12];
%! assert_equal (betainv (pp, 1, 1000), -expm1 (log1p (-pp) / 1000), -1e-14);
%!assert_equal (betainv ([1e-300, 1e-20, 0.3], 30, 1), ...
%!              [1e-300, 1e-20, 0.3] .^ (1/30), -1e-14)
%!assert_equal (betainv (2e-10, 0.5, 0.5), sin (pi * 1e-10) ^ 2, -1e-14)
## Expected values from MATLAB R2024a, among them the upper tail of a small
## first shape, where core's betaincinv misses
%!assert_equal (betainv (1e-10, 0.5, 50), 1.578669844608891e-22, -1e-12)
%!assert_equal (betainv (0.999, 0.5, 50), 0.103102263418713, -1e-13)
%!assert_equal (betainv (1 - 1e-9, 0.5, 1000), 0.018493954158234, -1e-13)
%!assert_equal (betainv (0.3, 2.5, 7.3), 0.170751332194546, -1e-14)
%!assert_equal (betainv (1e-6, 30, 200), 0.048760009765902, -1e-14)
%!assert_equal (betainv (1e-20, 2, 3), 4.082482904749744e-11, -1e-14)
%!assert_equal (betainv (0.975, 0.5, 0.5), 0.998458666866564, -1e-14)
%!assert_equal (size (betainv (0.3 * ones (2, 3), 2, 5)), [2, 3])

## Test class of input preserved
%!assert_equal (betainv ([p, NaN], 1, 2), [NaN 0 0.5 1 NaN NaN], eps)
%!assert_equal (betainv (single ([p, NaN]), 1, 2), single ([NaN 0 0.5 1 NaN NaN]))
%!assert_equal (betainv ([p, NaN], single (1), 2), single ([NaN 0 0.5 1 NaN NaN]), eps ('single'))
%!assert_equal (betainv ([p, NaN], 1, single (2)), single ([NaN 0 0.5 1 NaN NaN]), eps ('single'))

## Test input validation
%!error<betainv: function called with too few input arguments.> betainv ()
%!error<betainv: function called with too few input arguments.> betainv (1)
%!error<betainv: function called with too few input arguments.> betainv (1,2)
%!error<betainv: function called with too many inputs> betainv (1,2,3,4)
%!error<betainv: P, A, and B must be of common size or scalars.> ...
%! betainv (ones (3), ones (2), ones (2))
%!error<betainv: P, A, and B must be of common size or scalars.> ...
%! betainv (ones (2), ones (3), ones (2))
%!error<betainv: P, A, and B must be of common size or scalars.> ...
%! betainv (ones (2), ones (2), ones (3))
%!error<betainv: P, A, and B must be double or single.> betainv (int32 (2), 2, 2)
%!error<betainv: P, A, and B must be double or single.> betainv (true, 2, 2)
%!error<betainv: P, A, and B must be double or single.> betainv ('a', 2, 2)
%!error<betainv: P, A, and B must not be complex.> betainv (i, 2, 2)
%!error<betainv: P, A, and B must not be complex.> betainv (2, i, 2)
%!error<betainv: P, A, and B must not be complex.> betainv (2, 2, i)
