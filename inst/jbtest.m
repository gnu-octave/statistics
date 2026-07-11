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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{h} =} jbtest (@var{x})
## @deftypefnx {statistics} {@var{h} =} jbtest (@var{x}, @var{alpha})
## @deftypefnx {statistics} {@var{h} =} jbtest (@var{x}, @var{alpha}, @var{mctol})
## @deftypefnx {statistics} {[@var{h}, @var{p}] =} jbtest (@dots{})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{jbstat}, @var{critval}] =} jbtest (@dots{})
##
## Jarque-Bera hypothesis test of composite normality.
##
## @code{@var{h} = jbtest (@var{x})} performs the Jarque-Bera test of the null
## hypothesis that the sample in the vector @var{x} comes from a normal
## distribution with unknown mean and variance, against the alternative that it
## does not come from a normal distribution.  The result @var{h} is 1 if the
## test rejects the null hypothesis at the 5% significance level, and 0
## otherwise.  @var{x} must be a vector of real values; @qcode{NaN} values are
## treated as missing and removed.
##
## The Jarque-Bera test statistic is
## @tex
## $ JB = \frac{n}{6} \left( s^2 + \frac{(k-3)^2}{4} \right) $,
## @end tex
## @ifnottex
## @code{JB = (n / 6) * (s^2 + (k - 3)^2 / 4)},
## @end ifnottex
## where @math{n} is the sample size, @math{s} is the sample skewness, and
## @math{k} is the sample kurtosis.  Under the null hypothesis it is
## asymptotically chi-square distributed with two degrees of freedom.
##
## @code{@var{h} = jbtest (@var{x}, @var{alpha})} performs the test at the
## significance level @var{alpha}, a scalar in the range @math{(0,1)}.  The
## default is @math{0.05}.
##
## @code{@var{h} = jbtest (@var{x}, @var{alpha}, @var{mctol})} computes a
## Monte-Carlo approximation of the p-value instead of interpolating the
## embedded table.  @var{mctol} is the maximum Monte-Carlo standard
## error allowed for the p-value; the number of simulated samples is chosen
## accordingly.  Use this for small samples, where the chi-square approximation
## is inaccurate, or for significance levels outside @math{[0.001, 0.5]}.
##
## @code{[@var{h}, @var{p}] = jbtest (@dots{})} also returns the p-value
## @var{p} of the test.  @var{p} is clamped to the tabulated range
## @math{[0.001, 0.5]} and a warning is issued when the value lies outside it,
## matching MATLAB.
##
## @code{[@var{h}, @var{p}, @var{jbstat}, @var{critval}] = jbtest (@dots{})} also
## returns the test statistic @var{jbstat} and the critical value @var{critval}
## at significance level @var{alpha}.  The null hypothesis is rejected when
## @code{@var{jbstat} > @var{critval}}.
##
## Note: for @math{n \le 2000} the p-value and critical value are obtained by
## interpolating an embedded critical-value table (the same approach MATLAB
## uses); for larger samples the large-sample chi-square approximation with two
## degrees of freedom is used instead.  The embedded table was generated here by
## Monte-Carlo simulation, so it is itself an estimate of the true null
## quantiles.  MATLAB's table is likewise a Monte-Carlo estimate but from a
## different simulation, so the two tables agree only to about two decimal
## places.  As a result the reported p-value and critical value, and — in a
## narrow band of statistic values around the critical value — the test decision
## @var{h}, can differ slightly from MATLAB in edge cases.  These differences are
## an unavoidable consequence of the Monte-Carlo origin of both tables, not a
## difference in method.  Supply @var{mctol} for a direct Monte-Carlo p-value.
##
## @seealso{kstest, adtest, lillietest}
## @end deftypefn

function [h, p, jbstat, critval] = jbtest (x, alpha, mctol)

  ## Check input arguments
  if (nargin < 1)
    print_usage ();
  endif
  if (! (isnumeric (x) && isreal (x) && isvector (x)))
    error ("jbtest: X must be a vector of real values.");
  endif

  ## Remove missing values and check sample size
  x = x(! isnan (x));
  x = x(:);
  n = numel (x);
  if (n < 2)
    error ("jbtest: X must contain at least two non-missing values.");
  endif

  if (nargin < 2 || isempty (alpha))
    alpha = 0.05;
  endif
  if (! (isnumeric (alpha) && isscalar (alpha) && isreal (alpha) ...
         && alpha > 0 && alpha < 1))
    error ("jbtest: ALPHA must be a scalar in the range (0,1).");
  endif

  domc = (nargin > 2 && ! isempty (mctol));
  if (domc && ! (isnumeric (mctol) && isscalar (mctol) && isreal (mctol) ...
                 && mctol > 0))
    error ("jbtest: MCTOL must be a positive scalar.");
  endif

  ## Jarque-Bera test statistic (biased sample skewness and kurtosis)
  s = skewness (x);
  k = kurtosis (x);
  jbstat = (n / 6) * (s ^ 2 + (k - 3) ^ 2 / 4);

  if (! domc)
    [tsizes, talphas, tcv] = jbtest_table_ ();
    if (n > tsizes(end))
      ## Beyond the tabulated sample sizes: large-sample chi-square with two
      ## degrees of freedom
      p = 1 - chi2cdf (jbstat, 2);
      critval = chi2inv (1 - alpha, 2);
      if (p < 0.001)
        warning ("jbtest:pTooSmall", ...
                 "jbtest: P is less than the smallest tabulated value; returning 0.001.");
        p = 0.001;
      elseif (p > 0.5)
        warning ("jbtest:pTooBig", ...
                 "jbtest: P is greater than the largest tabulated value; returning 0.5.");
        p = 0.5;
      endif
    else
      ## Interpolate the embedded Monte-Carlo critical-value table.  Critical
      ## values as a function of the significance level at this sample size:
      nn = max (n, tsizes(1));
      cvn = interp1 (fliplr (1 ./ tsizes), flipud (tcv), 1 / nn, "linear");
      ## Critical value at the requested significance level (interpolated over
      ## the log of the tabulated levels, as in adtest)
      pp = pchip (log (talphas), cvn);
      aclamp = min (max (alpha, talphas(1)), talphas(end));
      critval = ppval (pp, log (aclamp));
      if (alpha < talphas(1) || alpha > talphas(end))
        warning ("jbtest:alphaRange", ...
                 "jbtest: ALPHA is outside the tabulated range [0.001, 0.5].");
      endif
      ## p-value by inverse interpolation into the table
      if (jbstat > cvn(1))
        warning ("jbtest:pTooSmall", ...
                 "jbtest: P is less than the smallest tabulated value; returning 0.001.");
        p = talphas(1);
      elseif (jbstat < cvn(end))
        warning ("jbtest:pTooBig", ...
                 "jbtest: P is greater than the largest tabulated value; returning 0.5.");
        p = talphas(end);
      else
        i = find (jbstat > cvn, 1, "first");
        logp = fzero (@(x) ppval (pp, x) - jbstat, log (talphas([i-1, i])));
        p = exp (logp);
      endif
    endif
  else
    ## Monte-Carlo approximation of the null distribution of the statistic
    reps = max (1000, ceil (0.25 / mctol ^ 2));
    jbsim = jbtest_simulate_ (n, reps);
    p = (1 + sum (jbsim >= jbstat)) / (reps + 1);
    critval = quantile (jbsim, 1 - alpha);
  endif

  h = double (jbstat > critval);

endfunction

## Simulate REPS values of the Jarque-Bera statistic for standard normal samples
## of size N, evaluated in chunks to bound memory use.
function jbsim = jbtest_simulate_ (n, reps)
  jbsim = zeros (reps, 1);
  chunk = max (1, floor (1e6 / n));
  done = 0;
  while (done < reps)
    m = min (chunk, reps - done);
    z = randn (n, m);
    s = skewness (z, 1, 1);          ## column-wise, biased
    k = kurtosis (z, 1, 1);
    jbsim(done+1:done+m) = (n / 6) .* (s(:) .^ 2 + (k(:) - 3) .^ 2 / 4);
    done += m;
  endwhile
endfunction

## Embedded Jarque-Bera critical-value table, generated by Monte-Carlo
## simulation (1e6 standard-normal samples per sample size).  Row i, column j is
## the upper @code{1 - ALPHAS(j)} quantile of the statistic for a sample of size
## @code{SIZES(i)}.  Values agree with MATLAB's table to about two decimals;
## both are Monte-Carlo estimates of the same quantiles.
function [sizes, alphas, cv] = jbtest_table_ ()
  sizes = [4 5 6 7 8 9 10 11 12 13 14 15 16 18 20 22 25 28 32 37 43 50 ...
           60 75 90 110 140 180 250 350 500 750 1000 2000];
  alphas = [0.001 0.0025 0.005 0.01 0.025 0.05 0.075 0.1 0.15 0.2 0.25 ...
            0.3 0.4 0.5];
  cv = [ ...
   0.9605  0.9570  0.9509  0.9393  0.9056  0.8522  0.8015  0.7554  0.6729  0.6307  0.5950  0.5630  0.5102  0.4739; ...
   1.8291  1.7788  1.7186  1.6283  1.4350  1.2179  1.0626  0.9425  0.7939  0.7298  0.6877  0.6512  0.5896  0.5289; ...
   3.1838  2.9786  2.7665  2.4815  2.0002  1.5549  1.2841  1.1007  0.9195  0.8301  0.7709  0.7229  0.6423  0.5735; ...
   4.9350  4.4225  3.9482  3.3862  2.5284  1.8493  1.4757  1.2531  1.0229  0.9146  0.8430  0.7859  0.6924  0.6091; ...
   6.9241  5.9457  5.1255  4.2348  2.9843  2.0890  1.6446  1.3906  1.1179  0.9911  0.9074  0.8423  0.7352  0.6418; ...
   8.8975  7.4734  6.2462  4.9863  3.3891  2.3208  1.8140  1.5169  1.2052  1.0616  0.9675  0.8934  0.7733  0.6699; ...
  10.9383  8.9031  7.3105  5.6927  3.7457  2.5276  1.9602  1.6254  1.2823  1.1240  1.0198  0.9389  0.8076  0.6943; ...
  12.9458 10.2945  8.2830  6.3396  4.0810  2.7105  2.0797  1.7233  1.3544  1.1818  1.0692  0.9812  0.8393  0.7181; ...
  14.7456 11.5915  9.1227  6.8706  4.3397  2.8724  2.2000  1.8123  1.4194  1.2347  1.1137  1.0193  0.8678  0.7392; ...
  16.6856 12.6696  9.8759  7.3578  4.6172  3.0267  2.3061  1.8960  1.4816  1.2858  1.1563  1.0558  0.8942  0.7574; ...
  18.1583 13.7946 10.6398  7.8586  4.8579  3.1663  2.4152  1.9802  1.5424  1.3333  1.1966  1.0897  0.9196  0.7751; ...
  19.3397 14.5252 11.1616  8.1792  5.0433  3.2907  2.5011  2.0505  1.5970  1.3771  1.2321  1.1202  0.9407  0.7903; ...
  21.1047 15.4426 11.8050  8.5872  5.2432  3.3999  2.5948  2.1230  1.6496  1.4209  1.2700  1.1521  0.9647  0.8075; ...
  22.8485 16.7542 12.7010  9.1686  5.6093  3.6350  2.7427  2.2420  1.7393  1.4936  1.3306  1.2031  1.0012  0.8325; ...
  25.0858 18.2343 13.6105  9.7801  5.9031  3.7947  2.8745  2.3471  1.8213  1.5622  1.3877  1.2519  1.0360  0.8567; ...
  26.4913 18.9880 14.0979 10.2355  6.1599  3.9611  3.0039  2.4524  1.8990  1.6255  1.4404  1.2959  1.0674  0.8788; ...
  28.4147 20.0407 14.9295 10.7391  6.4411  4.1481  3.1440  2.5749  1.9973  1.7047  1.5056  1.3506  1.1061  0.9045; ...
  29.8195 21.0006 15.4959 11.0930  6.6884  4.2933  3.2540  2.6754  2.0865  1.7778  1.5680  1.4044  1.1450  0.9309; ...
  31.2255 21.8267 16.1070 11.5171  6.9558  4.4934  3.4104  2.8046  2.1884  1.8618  1.6391  1.4628  1.1851  0.9590; ...
  32.0166 22.4555 16.5364 11.9241  7.1674  4.6567  3.5449  2.9289  2.2966  1.9503  1.7104  1.5225  1.2280  0.9883; ...
  33.2009 23.1281 17.0640 12.1954  7.3963  4.8309  3.6920  3.0593  2.4063  2.0435  1.7879  1.5880  1.2719  1.0165; ...
  33.4237 23.2982 17.1731 12.3563  7.5793  4.9742  3.8217  3.1810  2.5093  2.1293  1.8594  1.6467  1.3141  1.0447; ...
  33.5896 23.4273 17.3344 12.5402  7.7329  5.1294  3.9730  3.3281  2.6337  2.2310  1.9437  1.7175  1.3627  1.0782; ...
  33.1832 23.4014 17.3724 12.6873  7.8909  5.2810  4.1324  3.4803  2.7681  2.3450  2.0372  1.7951  1.4141  1.1114; ...
  32.1951 22.8349 17.1964 12.5965  7.9188  5.3719  4.2382  3.5961  2.8783  2.4343  2.1145  1.8597  1.4594  1.1399; ...
  30.7519 22.1102 16.7325 12.3753  7.9505  5.4573  4.3608  3.7217  2.9824  2.5237  2.1883  1.9204  1.5003  1.1689; ...
  30.1986 21.4958 16.4958 12.2626  7.9444  5.5527  4.4764  3.8553  3.1041  2.6213  2.2712  1.9881  1.5454  1.1974; ...
  27.6248 20.4179 15.6934 11.9179  7.9354  5.6402  4.6136  3.9895  3.2225  2.7224  2.3508  2.0582  1.5936  1.2290; ...
  26.3970 19.3670 15.1261 11.5818  7.8954  5.7581  4.7527  4.1311  3.3428  2.8256  2.4425  2.1320  1.6437  1.2631; ...
  23.9042 17.9448 14.3083 11.1508  7.7814  5.7994  4.8474  4.2331  3.4422  2.9143  2.5134  2.1895  1.6831  1.2874; ...
  21.3995 16.6461 13.4998 10.7451  7.6847  5.8609  4.9439  4.3331  3.5378  2.9931  2.5803  2.2471  1.7199  1.3113; ...
  19.8917 15.5804 12.8247 10.3682  7.6270  5.9036  5.0108  4.4173  3.6127  3.0576  2.6380  2.2946  1.7544  1.3325; ...
  18.8596 15.0900 12.4719 10.1555  7.5849  5.9373  5.0527  4.4619  3.6487  3.0892  2.6623  2.3156  1.7676  1.3407; ...
  16.3822 13.5505 11.5058  9.6830  7.4839  5.9564  5.1088  4.5241  3.7189  3.1529  2.7147  2.3586  1.7953  1.3615];
endfunction

%!demo
%! ## Test whether a sample departs from normality
%! x = [1 2 3 4 5 6 7 8 9 100];   # last value is an outlier
%! [h, p, jbstat] = jbtest (x)

## Test output against known values (chi-square approximation)
%!test
%! warning ("off", "jbtest:pTooBig", "local");
%! x = [1 2 3 4 5 6 7 8 9 10];
%! [h, p, jbstat, cv] = jbtest (x);
%! assert_equal (h, 0);
%! assert_equal (jbstat, 0.624487, 1e-5);    # skewness 0, kurtosis 1.7758
%! assert_equal (cv, 2.5276, 1e-4);          # tabulated critical value, n=10, a=0.05
%! assert_equal (p, 0.5, 1e-12);             # clamped to the tabulated maximum

## A very normal-looking (platykurtic) sample warns that p exceeds 0.5
%!warning <jbtest: P is greater than the largest tabulated value; returning 0.5.> ...
%! jbtest (1:10);

%!test  # a strongly non-normal sample is rejected
%! x = [zeros(1, 20), 100];
%! h = jbtest (x);
%! assert_equal (h, 1);

%!test  # NaNs are removed
%! assert_equal (jbtest ([1 2 3 4 5 6 7 8 9 10, NaN]), jbtest (1:10));

%!test  # alpha controls the critical value
%! x = randn (1, 50);
%! [~, ~, ~, cv1] = jbtest (x, 0.05);
%! [~, ~, ~, cv2] = jbtest (x, 0.01);
%! assert_equal (cv2 > cv1, true);

%!test  # Monte-Carlo p-value runs and lies in (0,1]
%! x = [1 2 3 4 5 6 7 8 9 10];
%! [h, p] = jbtest (x, 0.05, 0.05);
%! assert_equal (p > 0 && p <= 1, true);

## Test input validation
%!error <Invalid call to jbtest> jbtest ()
%!error <jbtest: X must be a vector of real values.> jbtest (ones (3, 3))
%!error <jbtest: X must be a vector of real values.> jbtest ({1, 2, 3})
%!error <jbtest: X must be a vector of real values.> jbtest ([1 2 3i])
%!error <jbtest: X must contain at least two non-missing values.> jbtest (5)
%!error <jbtest: ALPHA must be a scalar in the range .0,1..> jbtest (1:10, 0)
%!error <jbtest: ALPHA must be a scalar in the range .0,1..> jbtest (1:10, 1)
%!error <jbtest: ALPHA must be a scalar in the range .0,1..> jbtest (1:10, [0.1 0.2])
%!error <jbtest: MCTOL must be a positive scalar.> jbtest (1:10, 0.05, -1)
