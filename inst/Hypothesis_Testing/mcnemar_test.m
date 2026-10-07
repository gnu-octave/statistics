## Copyright (C) 1996-2017 Kurt Hornik
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
## @deftypefn  {statistics} {@var{h} =} mcnemar_test (@var{x})
## @deftypefnx {statistics} {@var{h} =} mcnemar_test (@var{x}, @var{testtype})
## @deftypefnx {statistics} {@var{h} =} mcnemar_test (@dots{}, @qcode{'Alpha'}, @var{alpha})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{stats}] =} mcnemar_test (@dots{})
##
## Perform a McNemar's test on paired nominal data.
##
## @nospell{McNemar's} test is applied to a @math{2*2} contingency table @var{x}
## with a dichotomous trait, with matched pairs of subjects, of data
## cross-classified on the row and column variables to testing the null
## hypothesis of symmetry of the classification probabilities.  More formally,
## the null hypothesis of marginal homogeneity states that the two marginal
## probabilities for each outcome are the same.  The test rests on the two
## discordant counts, @math{b = @var{x}(1,2)} and @math{c = @var{x}(2,1)}.
##
## Under the null, with a sufficiently large number of discordants
## (@math{b + c >= 25}), the test statistic follows a chi-squared distribution
## with 1 degree of freedom.  When the number of discordants is less than 25,
## then the mid-P exact McNemar test is used.
##
## @var{testtype} will force @code{mcnemar_test} to apply a particular method
## for testing the null hypothesis independently of the number of discordants.
## Valid options for @var{testtype}:
## @itemize
## @item @qcode{'asymptotic'} Original McNemar test statistic
## @item @qcode{'corrected'} Edwards' version with continuity correction
## @item @qcode{'exact'} An exact binomial test
## @item @qcode{'mid-p'} The mid-P McNemar test (mid-p binomial test)
## @end itemize
##
## The test decision is returned in @var{h}, which is 1 when the null hypothesis
## is rejected at the significance level @var{alpha} and 0 otherwise, and the
## p-value in @var{p}.  @qcode{'Alpha'} sets @var{alpha}, 0.05 by default, and
## the level @math{100 (1 - alpha)}% of the confidence intervals.  @var{stats}
## is a structure with the following fields:
##
## @multitable @columnfractions 0.2 0.75
## @item @qcode{chi2stat} @tab the chi-squared statistic, for the
## @qcode{'asymptotic'} and @qcode{'corrected'} tests only.
## @item @qcode{df} @tab its degrees of freedom, 1, for those tests only.
## @item @qcode{OddsRatio} @tab the odds ratio of the discordant pairs,
## @math{b / c}.
## @item @qcode{OddsRatioCI} @tab its confidence interval.
## @item @qcode{CohensG} @tab Cohen's @math{g}, @math{b / (b + c) - 0.5}, the
## departure of the discordant split from one half.
## @item @qcode{CohensGCI} @tab its confidence interval.
## @end multitable
##
## Both intervals come from the exact Clopper-Pearson interval of the
## proportion @math{b / (b + c)}, as does the exact test.  With no discordant
## pairs both effect sizes and their intervals are @qcode{NaN}.
##
## Further information about the McNemar's test can be found at
## @url{https://en.wikipedia.org/wiki/McNemar%27s_test}
##
## @seealso{crosstab, chi2test, fishertest}
## @end deftypefn

function [h, p, stats] = mcnemar_test (x, varargin)

  ## Check contingency table
  if (! isequal (size (x), [2, 2]))
    error ("mcnemar_test: X must be a 2x2 matrix.");
  elseif (! (all ((x(:) >= 0)) && all (x(:) == fix (x(:)))))
    error ("mcnemar_test: all entries of X must be non-negative integers.");
  endif

  ## Add defaults
  alpha = 0.05;
  b = x(1,2);
  c = x(2,1);
  if (b + c < 25)
    testtype = 'mid-p';
  else
    testtype = 'asymptotic';
  endif

  ## Parse optional arguments: a test type, and 'Alpha'
  i = 1;
  while (i <= numel (varargin))
    arg = varargin{i};
    if (ischar (arg) && strcmpi (arg, 'alpha'))
      if (i == numel (varargin))
        error ("mcnemar_test: optional arguments must be in pairs.");
      endif
      alpha = varargin{i+1};
      if (! (isnumeric (alpha) && isscalar (alpha) && isreal (alpha)
             && alpha > 0 && alpha < 1))
        error ("mcnemar_test: invalid value for alpha.");
      endif
      i += 2;
    elseif (ischar (arg))
      if (! any (strcmpi (arg, {'exact', 'asymptotic', 'mid-p', 'corrected'})))
        error ("mcnemar_test: invalid value for TESTTYPE.");
      endif
      testtype = lower (arg);
      i += 1;
    else
      error ("mcnemar_test: invalid optional argument.");
    endif
  endwhile

  ## Calculate test.  The exact tests are two-sided and symmetric in the
  ## discordant counts, so they read the tail of the smaller one.
  n = b + c;
  k = min (b, c);
  stats = struct ();
  switch (testtype)
    case 'asymptotic'
      stats.chi2stat = (b - c) .^2 / n;
      stats.df = 1;
      p = chi2cdf (stats.chi2stat, 1, 'upper');
    case 'corrected'
      stats.chi2stat = (abs (b - c) - 1) .^2 / n;
      stats.df = 1;
      p = chi2cdf (stats.chi2stat, 1, 'upper');
    case 'exact'
      p = min (1, 2 * binocdf (k, n, 0.5));
    case 'mid-p'
      p = min (1, 2 * binocdf (k, n, 0.5) - binopdf (k, n, 0.5));
  endswitch
  h = double (p < alpha);

  ## Effect sizes, with intervals from the Clopper-Pearson interval of the
  ## proportion of discordant pairs that fall in b
  if (n > 0)
    if (b == 0)
      lo = 0;
    else
      lo = betainv (alpha / 2, b, c + 1);
    endif
    if (c == 0)
      hi = 1;
    else
      hi = betainv (1 - alpha / 2, b + 1, c);
    endif
    stats.OddsRatio = b / c;
    stats.OddsRatioCI = [lo / (1 - lo), hi / (1 - hi)];
    stats.CohensG = b / n - 0.5;
    stats.CohensGCI = [lo, hi] - 0.5;
  else
    stats.OddsRatio = NaN;
    stats.OddsRatioCI = [NaN, NaN];
    stats.CohensG = NaN;
    stats.CohensGCI = [NaN, NaN];
  endif

endfunction

%!test
%! [h, p, st] = mcnemar_test ([101,121;59,33]);
%! assert_equal (h, 1);
%! assert_equal (p, 3.8151e-06, 1e-10);
%! assert_equal (st.chi2stat, 21.356, 1e-3);
%!test
%! [h, p, st] = mcnemar_test ([59,6;16,80]);
%! assert_equal (h, 1);
%! assert_equal (p, 0.034690, 1e-6);
%! assert_equal (isfield (st, 'chi2stat'), false);
%!test
%! [h, p] = mcnemar_test ([59,6;16,80], 'Alpha', 0.01);
%! assert_equal (h, 0);
%! assert_equal (p, 0.034690, 1e-6);
%!test
%! [h, p] = mcnemar_test ([59,6;16,80], 'mid-p');
%! assert_equal (h, 1);
%! assert_equal (p, 0.034690, 1e-6);
%!test
%! [h, p, st] = mcnemar_test ([59,6;16,80], 'asymptotic');
%! assert_equal (h, 1);
%! assert_equal (p, 0.033006, 1e-6);
%! assert_equal (st.chi2stat, 4.5455, 1e-4);
%! assert_equal (st.df, 1);
%!test
%! [h, p, st] = mcnemar_test ([59,6;16,80], 'exact');
%! assert_equal (h, 0);
%! assert_equal (p, 0.052479, 1e-6);
%! assert_equal (isfield (st, 'chi2stat'), false);
%!test
%! [h, p, st] = mcnemar_test ([59,6;16,80], 'corrected');
%! assert_equal (h, 0);
%! assert_equal (p, 0.055009, 1e-6);
%! assert_equal (st.chi2stat, 3.6818, 1e-4);
%!test
%! [h, p] = mcnemar_test ([59,6;16,80], 'corrected', 'Alpha', 0.1);
%! assert_equal (h, 1);
%! assert_equal (p, 0.055009, 1e-6);
%!test
%! ## Below the resolution of 1 - chi2cdf
%! [~, p, st] = mcnemar_test ([100, 200; 0, 100]);
%! assert_equal (p, erfc (sqrt (st.chi2stat / 2)), -1e-12);
%!test
%! ## The exact tests are symmetric in the discordant counts, values from R
%! [~, p] = mcnemar_test ([59,16;6,80], 'exact');
%! assert_equal (p, 0.052479, 1e-6);
%!test
%! [~, p] = mcnemar_test ([59,16;6,80], 'mid-p');
%! assert_equal (p, 0.034690, 1e-6);
%!test
%! [~, p] = mcnemar_test ([5,7;7,5], 'exact');
%! assert_equal (p, 1);
%!test
%! [~, ~, st] = mcnemar_test ([59,6;16,80]);
%! assert_equal (st.OddsRatio, 6 / 16);
%!test
%! ## Clopper-Pearson interval of 6/22 from R's binom.test
%! [~, ~, st] = mcnemar_test ([59,6;16,80]);
%! assert_equal (st.OddsRatioCI, [0.12018366326892038, 1.0089244510694295], ...
%!               -1e-12);
%!test
%! [~, ~, st] = mcnemar_test ([59,6;16,80]);
%! assert_equal (st.CohensG, 6 / 22 - 0.5, -1e-14);
%!test
%! [~, ~, st] = mcnemar_test ([59,6;16,80]);
%! assert_equal (st.CohensGCI, ...
%!               [0.10728924837039702, 0.50222120126634895] - 0.5, -1e-12);
%!test
%! [~, ~, st] = mcnemar_test ([101,121;59,33], 'Alpha', 0.01);
%! assert_equal (st.OddsRatioCI, [1.3562696986628573, 3.1577641330443376], ...
%!               -1e-12);
%!test
%! ## With no discordant pairs there is no effect to measure
%! [~, ~, st] = mcnemar_test ([5,0;0,5]);
%! assert_equal ([st.OddsRatio, st.CohensG], [NaN, NaN]);
%!test
%! [~, ~, st] = mcnemar_test ([5,3;0,5]);
%! assert_equal ([st.OddsRatio, st.OddsRatioCI(2)], [Inf, Inf]);

%!error<mcnemar_test: X must be a 2x2 matrix.> mcnemar_test (59, 6, 16, 80)
%!error<mcnemar_test: X must be a 2x2 matrix.> mcnemar_test (ones (3, 3))
%!error<mcnemar_test: all entries of X must be non-negative integers.> ...
%! mcnemar_test ([59,6;16,-80])
%!error<mcnemar_test: all entries of X must be non-negative integers.> ...
%! mcnemar_test ([59,6;16,4.5])
%!error<mcnemar_test: invalid optional argument.> ...
%! mcnemar_test ([59,6;16,80], {''})
%!error<mcnemar_test: invalid optional argument.> ...
%! mcnemar_test ([59,6;16,80], 0.05)
%!error<mcnemar_test: optional arguments must be in pairs.> ...
%! mcnemar_test ([59,6;16,80], 'Alpha')
%!error<mcnemar_test: invalid value for alpha.> ...
%! mcnemar_test ([59,6;16,80], 'Alpha', -0.2)
%!error<mcnemar_test: invalid value for alpha.> ...
%! mcnemar_test ([59,6;16,80], 'Alpha', [0.05, 0.1])
%!error<mcnemar_test: invalid value for alpha.> ...
%! mcnemar_test ([59,6;16,80], 'Alpha', 1)
%!error<mcnemar_test: invalid value for TESTTYPE.> ...
%! mcnemar_test ([59,6;16,80], '')
