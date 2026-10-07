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
## Where one degree of freedom is so large against the other that the
## @math{F}-distribution cannot be told from its limit in double precision,
## the quantile is taken from the limit, @code{chi2inv (@var{p}, @var{df1}) /
## @var{df1}} as @var{df2} grows and @code{@var{df2} / chi2inv (1 - @var{p},
## @var{df2})} as @var{df1} grows.  Once one degree of freedom reaches
## @math{10^8} and the other is small, MATLAB returns quantiles whose
## probability is off by about 3%: its @code{finv (0.025, 10, 1e8)} is
## 0.32227, whose probability is 0.0243, where the quantile is 0.32470.
##
## Input arguments must be @qcode{double} or @qcode{single}; integer, logical,
## and character arrays are rejected.  MATLAB accepts a character array and
## evaluates it at the character codes, which Octave deliberately does not,
## since a character array is an integer type and integers are refused too.
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

  ## Check for P, DF1, and DF2 being double or single
  if (! (isfloat (p) && isfloat (df1) && isfloat (df2)))
    error ("finv: P, DF1, and DF2 must be double or single.");
  endif

  ## Check for P, DF1, and DF2 being reals
  if (iscomplex (p) || iscomplex (df1) || iscomplex (df2))
    error ("finv: P, DF1, and DF2 must not be complex.");
  endif

  ## Check for class type
  if (isa (p, 'single') || isa (df1, 'single') || isa (df2, 'single'))
    x = NaN (size (p), 'single');
  else
    x = NaN (size (p));
  endif

  ## Handle both DFs being INF
  kz = df1 == Inf & df2 == Inf;

  k = p == 1 & df1 > 0 & df2 > 0;
  x(k) = Inf;

  ## A DF so large against the other that the Beta form loses more to
  ## rounding than the limit loses to the finite DF is taken as infinite; the
  ## two errors cross where the larger DF squared is 3e15 times the smaller
  big = max (df1, df2) .^ 2 > 3e15 * max (min (df1, df2), 1);
  inf1 = (df1 == Inf | (big & df1 > df2)) & df2 < Inf;
  inf2 = (df2 == Inf | (big & df2 > df1)) & df1 < Inf;

  ## Solve on the tail nearer the answer: the lower one as Beta (DF1/2, DF2/2)
  ## and the upper one as Beta (DF2/2, DF1/2), so that neither cancels
  k = (p >= 0) & (p < 1) & (df1 > 0) & (df2 > 0) & ! inf1 & ! inf2 & ! kz;
  kl = k & (p <= 0.5);
  z = betainv (p(kl), df1(kl) / 2, df2(kl) / 2);
  x(kl) = z ./ (1 - z) .* df2(kl) ./ df1(kl);
  ku = k & (p > 0.5);
  w = betainv (1 - p(ku), df2(ku) / 2, df1(ku) / 2);
  x(ku) = (1 - w) ./ w .* df2(ku) ./ df1(ku);

  ## Limits: DF1 * X is chi-square with DF1 as DF2 grows, and DF2 / X is
  ## chi-square with DF2 as DF1 grows
  k = (p >= 0) & (p < 1) & (df1 > 0) & inf2;
  x(k) = chi2inv (p(k), df1(k)) ./ df1(k);
  k = (p >= 0) & (p < 1) & (df2 > 0) & inf1;
  x(k) = df2(k) ./ (2 * gammaincinv (p(k), df2(k) / 2, 'upper'));

  ## Force instances with df1 = df2 = INF to 0 for p = 0 and to 1 for 0 < p <= 1
  x(kz & p > 0 & p <= 1) = 1;
  x(kz & p == 0) = 0;

endfunction

%!demo
%! ## Plot various iCDFs from the F distribution
%! p = 0.001:0.001:0.999;
%! x1 = finv (p, 1, 1);
%! x2 = finv (p, 2, 1);
%! x3 = finv (p, 5, 2);
%! x4 = finv (p, 10, 1);
%! x5 = finv (p, 100, 100);
%! plot (p, x1, '-b', p, x2, '-g', p, x3, '-r', p, x4, '-c', p, x5, '-m')
%! grid on
%! ylim ([0, 4])
%! legend ({'df1 = 1, df2 = 2', 'df1 = 2, df2 = 1', ...
%!          'df1 = 5, df2 = 2', 'df1 = 10, df2 = 1', ...
%!          'df1 = 100, df2 = 100'}, 'location', 'northwest')
%! title ('F iCDF')
%! xlabel ('probability')
%! ylabel ('values in x')

## Test output
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert_equal (finv (p, 2*ones (1,5), 2*ones (1,5)), [NaN 0 1 Inf NaN])
%!assert_equal (finv (p, 2, 2*ones (1,5)), [NaN 0 1 Inf NaN])
%!assert_equal (finv (p, 2*ones (1,5), 2), [NaN 0 1 Inf NaN])
%!assert_equal (finv (p, [2 -Inf NaN Inf 2], 2), [NaN NaN NaN Inf NaN])
%!assert_equal (finv (p, 2, [2 -Inf NaN Inf 2]), [NaN NaN NaN Inf NaN])
%!assert_equal (finv ([p(1:2) NaN p(4:5)], 2, 2), [NaN 0 NaN Inf NaN])

## Test for bug #66034 (savannah)
%!assert_equal (finv (0.025, 10, 1e6), 0.3247, 1e-4)
%!assert_equal (finv (0.025, 10, 1e7), 0.3247, 1e-4)
%!assert_equal (finv (0.025, 10, 1e10), 0.3247, 1e-4)
%!assert_equal (finv (0.025, 10, 1e255), 0.3247, 1e-4)
%!assert_equal (finv (0.025, 10, Inf), 0.3247, 1e-4)

## Test for issue #203 (Github)
%!test
%! x = finv (0.35, Inf, 4);
%! assert_equal (x, 0.9014, 1e-4)
%!test
%! x = finv (0, Inf, 4);
%! assert_equal (x, 0)
%!test
%! x = finv (1, Inf, 4);
%! assert_equal (x, Inf)
%!test
%! x = finv (0.35, 4, Inf);
%! assert_equal (x, 0.6175, 1e-4)
%!test
%! x = finv (0, 4, Inf);
%! assert_equal (x, 0)
%!test
%! x = finv (1, 4, Inf);
%! assert_equal (x, Inf)
%!test
%! x = finv ([0, 0.000001, 0.35, 1, 1.2], Inf, Inf);
%! assert_equal (x, [0, 1, 1, 1, NaN]);

## Values from MATLAB R2024a
%!assert_equal (finv (1e-8, 2, 2), 1.00000001e-08, -1e-14)
%!assert_equal (finv (1e-12, 2, 2), 1.000000000001e-12, -1e-14)
%!assert_equal (finv (1e-10, 1, 1), 2.46740110027235e-20, -1e-14)
%!assert_equal (finv (1e-6, 10, 20), 0.0285838655483384, -1e-14)
%!assert_equal (finv (0.975, 1e6, 1e6), 1.00392762317843, -1e-10)
%!assert_equal (finv (0.975, 2e6, 2e6), 1.00277565345099, -1e-10)
%!assert_equal (finv (0.975, 1e7, 1e7), 1.00124035874565, -1e-10)
%!assert_equal (finv (0.975, 1e8, 1e8), 1.00039206963839, -1e-10)
%!assert_equal (finv (0.975, 1e8, 1e10), 1.00027858274271, -1e-9)
%!assert_equal (finv (0.35, Inf, 4), 0.901370019458443, -1e-14)
%!assert_equal (finv (0.35, 4, Inf), 0.617521846868827, -1e-14)

## The limits, where MATLAB R2024a gives 0.322269543148685, 8.90372005876577
## and 2.56334042013798e-05
%!assert_equal (finv (0.025, 10, 1e12), 0.32469727802368437, -1e-10)
%!assert_equal (finv (0.975, 1e12, 4), 8.2573219821426864, -1e-10)
%!assert_equal (finv (1e-10, 10, 1e12), 0.0052331065631905406, -1e-10)

## Test class of input preserved
%!assert_equal (finv ([p, NaN], 2, 2), [NaN 0 1 Inf NaN NaN])
%!assert_equal (finv (single ([p, NaN]), 2, 2), single ([NaN 0 1 Inf NaN NaN]))
%!assert_equal (finv ([p, NaN], single (2), 2), single ([NaN 0 1 Inf NaN NaN]))
%!assert_equal (finv ([p, NaN], 2, single (2)), single ([NaN 0 1 Inf NaN NaN]))

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
%!error<finv: P, DF1, and DF2 must be double or single.> finv (int32 (2), 2, 2)
%!error<finv: P, DF1, and DF2 must be double or single.> finv (true, 2, 2)
%!error<finv: P, DF1, and DF2 must be double or single.> finv ('a', 2, 2)
%!error<finv: P, DF1, and DF2 must not be complex.> finv (i, 2, 2)
%!error<finv: P, DF1, and DF2 must not be complex.> finv (2, i, 2)
%!error<finv: P, DF1, and DF2 must not be complex.> finv (2, 2, i)
