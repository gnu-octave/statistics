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
## @deftypefn  {statistics} {@var{h} =} swtest (@var{x})
## @deftypefnx {statistics} {@var{h} =} swtest (@var{x}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{h}, @var{p}] =} swtest (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{swstat}, @var{critval}] =} swtest (@dots{})
##
## Shapiro-Wilk hypothesis test of composite normality.
##
## @code{@var{h} = swtest (@var{x})} performs the Shapiro-Wilk test of the null
## hypothesis that the sample in the vector @var{x} comes from a normal
## distribution with unknown mean and variance, against the alternative that it
## does not come from a normal distribution.  The result @var{h} is 1 if the
## test rejects the null hypothesis at the 5% significance level, and 0
## otherwise.  @var{x} must be a vector of real values; @qcode{NaN} values are
## treated as missing and removed.
##
## The test statistic is
## @tex
## $$ W = \frac{\left( \sum_{i=1}^n a_i x_{(i)} \right)^2}
##             {\sum_{i=1}^n (x_i - \bar{x})^2}, $$
##
## @end tex
## @ifnottex
## @code{W = (sum (a .* sort (x)))^2 / sum ((x - mean (x)).^2)},
## @end ifnottex
## where the weights @math{a}, applied to the sorted sample, are derived from
## the expected values of the order statistics of a standard normal sample of
## size @math{n}.  @math{W} lies in
## @math{(0,1]}, and small values of @math{W} are evidence against normality.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.2 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Alpha'} @tab The significance level, a scalar in the range
## @math{(0,1)}.  The default is @math{0.05}.
##
## @item @qcode{'Method'} @tab The test to perform:
## @qcode{'shapiro-wilk'} (default) or @qcode{'shapiro-francia'}.
## @end multitable
##
## @qcode{'shapiro-wilk'} computes the weights and the p-value with Royston's
## algorithm AS R94, which approximates the weights of Shapiro and Wilk and
## transforms @math{W} to a normal deviate.  It is the algorithm used by R's
## @code{shapiro.test}, and it accepts samples of 3 to 5000 values.  For
## @math{n = 3} the p-value is exact.
##
## @qcode{'shapiro-francia'} performs the Shapiro-Francia test instead, whose
## statistic @math{W'} is the squared correlation between the sorted sample and
## the approximate expected normal order statistics
## @code{norminv (((1:n) - 3/8) / (n + 1/4))}.  Its p-value is Royston's
## normal approximation for @math{W'}, as in @code{sf.test} of R's
## @code{nortest} package, and it accepts samples of 5 to 5000 values.  The
## Shapiro-Francia test is known to be more powerful than the Shapiro-Wilk test
## against leptokurtic alternatives.  The test is never selected automatically;
## the method used is the one requested.
##
## @code{[@var{h}, @var{p}] = swtest (@dots{})} also returns the p-value
## @var{p} of the test, the probability of observing a statistic as small as
## @var{swstat} under the null hypothesis.
##
## @code{[@var{h}, @var{p}, @var{swstat}, @var{critval}] = swtest (@dots{})}
## also returns the test statistic @var{swstat}, @math{W} or @math{W'}, and the
## critical value @var{critval} at significance level @var{alpha}, obtained by
## inverting the same approximation.  The null hypothesis is rejected when
## @code{@var{swstat} < @var{critval}}, which is the same as
## @code{@var{p} < @var{alpha}}.
##
## A sample whose values are all equal has no defined statistic and is refused,
## as is a sample holding an infinite value.  Samples of more than 5000 values
## are refused, since Royston's approximations are calibrated up to that size;
## use @code{adtest} or @code{jbtest} for larger samples.
##
## MATLAB has no Shapiro-Wilk test, so @code{swtest} is specific to Octave.
##
## References:
## @enumerate
## @item
## S. S. Shapiro and M. B. Wilk.  An analysis of variance test for normality
## (complete samples).  @emph{Biometrika}, 52(3-4):591--611, 1965.
## @item
## S. S. Shapiro and R. S. Francia.  An approximate analysis of variance test
## for normality.  @emph{Journal of the American Statistical Association},
## 67(337):215--216, 1972.
## @item
## P. Royston.  A pocket-calculator algorithm for the Shapiro-Francia test for
## non-normality: an application to medicine.  @emph{Statistics in Medicine},
## 12(2):181--184, 1993.
## @item
## P. Royston.  Remark AS R94: a remark on algorithm AS 181: the W-test for
## normality.  @emph{Applied Statistics}, 44(4):547--551, 1995.
## @end enumerate
##
## @seealso{adtest, jbtest, kstest, lillietest}
## @end deftypefn

function [h, p, swstat, critval] = swtest (x, varargin)

  ## Input validation
  if (nargin < 1)
    print_usage ();
  endif
  if (! (isnumeric (x) && isreal (x) && isvector (x)))
    error ("swtest: X must be a vector of real values.");
  endif

  ## Parse optional Name-Value paired arguments
  optNames = {'Alpha', 'Method'};
  dfValues = {0.05, 'shapiro-wilk'};
  [alpha, method, args] = parsePairedArguments (optNames, dfValues, ...
                                                varargin(:));
  if (! isempty (args))
    error (strcat ("swtest: optional arguments must be 'Alpha' or", ...
                   " 'Method' Name-Value pairs."));
  endif
  if (! (isnumeric (alpha) && isscalar (alpha) && isreal (alpha) ...
         && alpha > 0 && alpha < 1))
    error ("swtest: 'Alpha' must be a scalar in the range (0,1).");
  endif
  if (isstring (method) && isscalar (method))
    method = char (method);
  endif
  if (! (ischar (method) && isrow (method) ...
         && any (strcmpi (method, {'shapiro-wilk', 'shapiro-francia'}))))
    error (strcat ("swtest: 'Method' must be 'shapiro-wilk' or", ...
                   " 'shapiro-francia'."));
  endif
  francia = strcmpi (method, 'shapiro-francia');

  ## Remove missing values and check the sample
  x = sort (double (x(! isnan (x))(:)));
  n = numel (x);
  if (any (isinf (x)))
    error ("swtest: X must not contain infinite values.");
  endif
  if (francia && n < 5)
    error (strcat ("swtest: X must contain at least five non-missing", ...
                   " values for the Shapiro-Francia test."));
  elseif (n < 3)
    error ("swtest: X must contain at least three non-missing values.");
  endif
  if (n > 5000)
    error ("swtest: X must contain at most 5000 non-missing values.");
  endif
  if (x(n) == x(1))
    error ("swtest: X must not be constant.");
  endif

  ## Scale by the range and centre the sample
  x = x / (x(n) - x(1));
  x = x - mean (x);

  if (francia)
    ## Approximate expected normal order statistics (Blom's scores)
    a = norminv (((1:n)' - 0.375) / (n + 0.25));
  else
    a = swtest_weights_ (n);
  endif

  ## 1 - W as the squared sine of the angle between A and X, which keeps its
  ## precision when W is close to 1.  It is clamped at 0 because rounding can
  ## take it just below.
  ssa = sumsq (a);
  ssx = sumsq (x);
  sax = a' * x;
  ssassx = sqrt (ssa * ssx);
  w1 = max ((ssassx - sax) * (ssassx + sax) / (ssa * ssx), 0);
  swstat = 1 - w1;

  if (francia)
    ## Royston (1993): log (1 - W') is approximately normal
    u = log (n);
    v = log (u);
    mu = -1.2725 + 1.0521 * (v - u);
    sigma = 1.0308 - 0.26758 * (v + 2 / u);
    p = normcdf ((log (w1) - mu) / sigma, 'upper');
    critval = -expm1 (mu - sigma * norminv (alpha));
  elseif (n == 3)
    ## Exact null distribution for three values
    p = min (max (6 / pi * (asin (sqrt (swstat)) - pi / 3), 0), 1);
    critval = sin (pi / 3 + alpha * pi / 6) ^ 2;
  elseif (n <= 11)
    ## Royston (1995): -log (gamma - log (1 - W)) is approximately normal.
    ## GMA exceeds log (1 - W) for every attainable W.
    gma = polyval ([0.459, -2.273], n);
    mu = polyval ([-6.714e-4, 0.025054, -0.39978, 0.544], n);
    sigma = exp (polyval ([-0.0020322, 0.062767, -0.77857, 1.3822], n));
    p = normcdf ((-log (gma - log (w1)) - mu) / sigma, 'upper');
    critval = -expm1 (gma - exp (-(mu - sigma * norminv (alpha))));
  else
    ## Royston (1995): log (1 - W) is approximately normal
    u = log (n);
    mu = polyval ([0.0038915, -0.083751, -0.31082, -1.5861], u);
    sigma = exp (polyval ([0.0030302, -0.082676, -0.4803], u));
    p = normcdf ((log (w1) - mu) / sigma, 'upper');
    critval = -expm1 (mu - sigma * norminv (alpha));
  endif

  h = double (swstat < critval);

endfunction

## Shapiro-Wilk weights for a sample of size N by Royston's approximation
## (algorithm AS R94), returned antisymmetric and ordered for sorted data.
function a = swtest_weights_ (n)
  nh = floor (n / 2);
  if (n == 3)
    ah = sqrt (0.5);
  else
    m = norminv (((1:nh)' - 0.375) / (n + 0.25));
    summ2 = 2 * sumsq (m);
    ssumm2 = sqrt (summ2);
    rsn = 1 / sqrt (n);
    a1 = polyval ([-2.706056, 4.434685, -2.07119, -0.147981, 0.221157, 0], ...
                  rsn) - m(1) / ssumm2;
    if (n > 5)
      a2 = polyval ([-3.582633, 5.682633, -1.752461, -0.293762, 0.042981, ...
                     0], rsn) - m(2) / ssumm2;
      fac = sqrt ((summ2 - 2 * m(1) ^ 2 - 2 * m(2) ^ 2) ...
                  / (1 - 2 * a1 ^ 2 - 2 * a2 ^ 2));
      ah = [a1; a2; -m(3:nh) / fac];
    else
      fac = sqrt ((summ2 - 2 * m(1) ^ 2) / (1 - 2 * a1 ^ 2));
      ah = [a1; -m(2:nh) / fac];
    endif
  endif
  a = zeros (n, 1);
  a(1:nh) = -ah;
  a(n:-1:n-nh+1) = ah;
endfunction

%!demo
%! ## Test whether a sample departs from normality
%! x = [148 154 158 160 161 162 166 170 182 195 236];
%! [h, p, W] = swtest (x)

%!demo
%! ## The Shapiro-Francia test on the same sample
%! x = [148 154 158 160 161 162 166 170 182 195 236];
%! [h, p, W] = swtest (x, 'Method', 'shapiro-francia')

## Test output
## W and P checked against shapiro.test in R 4.3.3 (2024-02-29)
%!test  # Shapiro and Wilk's example, n <= 11
%! x = [148 154 158 160 161 162 166 170 182 195 236];
%! [h, p, W, c] = swtest (x);
%! assert_equal (h, 1);
%! assert_equal (p, 0.00670381405650293, 1e-12);
%! assert_equal (W, 0.788814694835387, 1e-12);
%! assert_equal (c, 0.855278601529214, 1e-12);
%!test  # n > 11
%! [h, p, W, c] = swtest (log (1:25));
%! assert_equal (h, 1);
%! assert_equal (p, 0.00741545379169184, 1e-12);
%! assert_equal (W, 0.881481586050452, 1e-12);
%! assert_equal (c, 0.91953513814811, 1e-12);
%!test  # n = 5, two approximated weights
%! [h, p, W, c] = swtest ([2 3 5 8 13]);
%! assert_equal (h, 0);
%! assert_equal (p, 0.534654754257205, 1e-12);
%! assert_equal (W, 0.920729192441235, 1e-12);
%! assert_equal (c, 0.775099786396077, 1e-12);
%!test  # n = 4
%! [h, p, W, c] = swtest ([1 2 4 8]);
%! assert_equal (h, 0);
%! assert_equal (p, 0.53808377727497, 1e-12);
%! assert_equal (W, 0.92020267879194, 1e-12);
%! assert_equal (c, 0.762289320516799, 1e-12);
%!test  # n = 3, exact null distribution
%! [h, p, W, c] = swtest ([1 2 4]);
%! assert_equal (h, 0);
%! assert_equal (W, 27 / 28, 1e-14);
%! assert_equal (p, 6 / pi * (asin (sqrt (27 / 28)) - pi / 3), 1e-14);
%! assert_equal (c, sin (pi / 3 + pi / 120) ^ 2, 1e-14);
## R 4.3.3 gives P = 0.99999999999999334 here, its 6/pi and pi/3 being
## written to 15 digits
%!test  # three equally spaced values fit a normal sample exactly
%! [h, p, W] = swtest ([1 2 3]);
%! assert_equal (h, 0);
%! assert_equal (p, 1);
%! assert_equal (W, 1);
## W' and P checked in R 4.3.3 (2024-02-29) against the formula of sf.test
## from the nortest package, which was not installed
%!test  # Shapiro-Francia, n <= 11
%! x = [148 154 158 160 161 162 166 170 182 195 236];
%! [h, p, W, c] = swtest (x, 'Method', 'shapiro-francia');
%! assert_equal (h, 1);
%! assert_equal (p, 0.00734764001456067, 1e-12);
%! assert_equal (W, 0.771381939646386, 1e-12);
%! assert_equal (c, 0.855095971555222, 1e-12);
%!test  # Shapiro-Francia, n > 11
%! [h, p, W, c] = swtest (log (1:25), 'Method', 'shapiro-francia');
%! assert_equal (h, 1);
%! assert_equal (p, 0.00991400561293846, 1e-12);
%! assert_equal (W, 0.882794840134263, 1e-12);
%! assert_equal (c, 0.919670595475726, 1e-12);
%!test  # Shapiro-Francia, n = 5
%! [h, p, W, c] = swtest ([2 3 5 8 13], 'Method', 'shapiro-francia');
%! assert_equal (h, 0);
%! assert_equal (p, 0.592971063346381, 1e-12);
%! assert_equal (W, 0.925681475433048, 1e-12);
%! assert_equal (c, 0.782592834233546, 1e-12);
%!test  # Shapiro-Francia, n = 200, measured with R's nortest sf.test
%! [h, p, W] = swtest (sqrt (1:200), 'Method', 'shapiro-francia');
%! assert_equal (h, 1);
%! assert_equal (p, 8.0650594412092e-06, -1e-10);
%! assert_equal (W, 0.950552730992928, -1e-10);
%!test  # Shapiro-Francia on a strong skew, measured with R's nortest sf.test
%! [h, p, W] = swtest (exp ((1:20) / 5), 'Method', 'shapiro-francia');
%! assert_equal (h, 1);
%! assert_equal (p, 0.00305794323924647, -1e-10);
%! assert_equal (W, 0.824238460116288, -1e-10);
%!test  # Shapiro-Francia on bounded data, measured with R's nortest sf.test
%! [h, p, W] = swtest (sin (1:30), 'Method', 'shapiro-francia');
%! assert_equal (h, 1);
%! assert_equal (p, 0.018296975230196, -1e-10);
%! assert_equal (W, 0.911248253555032, -1e-10);
%!test  # Alpha sets the critical value and leaves P unchanged
%! [h, p, W, c] = swtest (log (1:25), 'Alpha', 0.01);
%! assert_equal (h, 1);
%! assert_equal (p, 0.00741545379169184, 1e-12);
%! assert_equal (c, 0.887697842838192, 1e-12);
%!test  # Alpha decides H
%! assert_equal (swtest ([2 3 5 8 13], 'Alpha', 0.6), 1);
%!test  # Method is case insensitive
%! [~, p] = swtest (log (1:25), 'Method', 'Shapiro-Francia');
%! assert_equal (p, 0.00991400561293846, 1e-12);
%!test  # Method as a string scalar
%! [~, p] = swtest (log (1:25), 'Method', string ('shapiro-francia'));
%! assert_equal (p, 0.00991400561293846, 1e-12);
%!test  # NaNs are removed
%! [h, p, W, c] = swtest ([NaN, log(1:25), NaN]);
%! assert_equal ([h, p, W, c], ...
%!               [1, 0.00741545379169184, 0.881481586050452, ...
%!                0.91953513814811], 1e-12);
%!test  # a column vector gives the same result
%! [h, p, W, c] = swtest (log (1:25)');
%! assert_equal ([h, p, W, c], ...
%!               [1, 0.00741545379169184, 0.881481586050452, ...
%!                0.91953513814811], 1e-12);
%!test  # unsorted data gives the same result
%! [h, p, W, c] = swtest (fliplr (log (1:25)));
%! assert_equal ([h, p, W, c], ...
%!               [1, 0.00741545379169184, 0.881481586050452, ...
%!                0.91953513814811], 1e-12);
%!test  # integer data is tested as double
%! [~, p, W] = swtest (int16 ([148 154 158 160 161 162 166 170 182 195 236]));
%! assert_equal (p, 0.00670381405650293, 1e-12);
%! assert_equal (W, 0.788814694835387, 1e-12);
%!test  # single data is tested as double
%! [~, p, W] = swtest (single ([148 154 158 160 161 162 166 170 182 195 236]));
%! assert_equal (p, 0.00670381405650293, 1e-12);
%! assert_equal (W, 0.788814694835387, 1e-12);
%!test  # a sample of 5000 values is accepted
%! assert_equal (swtest (1:5000), 1);

## Test input validation
%!error <Invalid call to swtest> swtest ()
%!error <swtest: X must be a vector of real values.> swtest (ones (3, 3))
%!error <swtest: X must be a vector of real values.> swtest ({1, 2, 3})
%!error <swtest: X must be a vector of real values.> swtest ([1 2 3i])
%!error <swtest: X must be a vector of real values.> swtest ('abcde')
%!error <swtest: optional arguments must be 'Alpha' or 'Method' Name-Value pairs.> ...
%! swtest (1:10, 0.01)
%!error <swtest: optional arguments must be 'Alpha' or 'Method' Name-Value pairs.> ...
%! swtest (1:10, 'Tail', 'left')
%!error <swtest: 'Alpha' must be a scalar in the range .0,1..> ...
%! swtest (1:10, 'Alpha', 0)
%!error <swtest: 'Alpha' must be a scalar in the range .0,1..> ...
%! swtest (1:10, 'Alpha', 1)
%!error <swtest: 'Alpha' must be a scalar in the range .0,1..> ...
%! swtest (1:10, 'Alpha', [0.01 0.05])
%!error <swtest: 'Alpha' must be a scalar in the range .0,1..> ...
%! swtest (1:10, 'Alpha', '0.05')
%!error <swtest: 'Method' must be 'shapiro-wilk' or 'shapiro-francia'.> ...
%! swtest (1:10, 'Method', 'sf')
%!error <swtest: 'Method' must be 'shapiro-wilk' or 'shapiro-francia'.> ...
%! swtest (1:10, 'Method', 5)
%!error <swtest: X must not contain infinite values.> swtest ([1 2 3 Inf])
%!error <swtest: X must contain at least five non-missing values for the Shapiro-Francia test.> ...
%! swtest (1:4, 'Method', 'shapiro-francia')
%!error <swtest: X must contain at least three non-missing values.> ...
%! swtest ([1 2 NaN])
%!error <swtest: X must contain at most 5000 non-missing values.> ...
%! swtest (1:5001)
%!error <swtest: X must not be constant.> swtest (ones (1, 10))
