## Copyright (C) 2016 Andreas Stahel<Andreas.Stahel@bfh.ch>
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
## @deftypefn  {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} binotest (@var{pos}, @var{N}, @var{p0})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} binotest (@var{pos}, @var{N}, @var{p0}, @var{Name}, @var{Value})
##
## Test for probability @var{p} of a binomial sample
##
## Perform a test of the null hypothesis @var{p} ==  @var{p0} for a sample
## of size @var{N} with @var{pos} positive results.
##
##
## Name-Value pair arguments can be used to set various options.
## @qcode{'alpha'} can be used to specify the significance level
## of the test (the default value is 0.05). The option @qcode{'tail'},
## can be used to select the desired alternative hypotheses.  If the
## value is @qcode{'both'} (default) the null is tested against the two-sided
## alternative @code{@var{p} != @var{p0}}. The value of @var{pval} is
## determined by adding the probabilities of all events less or equally
## likely than the observed number @var{pos} of positive events, equally
## likely to within a relative 1e-7, as R's @code{binom.test} compares them.
## If the value of @qcode{'tail'} is @qcode{'right'}
## the one-sided alternative @code{@var{p} > @var{p0}} is considered.
## Similarly for @qcode{'left'}, the one-sided alternative
## @code{@var{p} < @var{p0}} is considered.
##
## If @var{h} is 0 the null hypothesis is accepted, if it is 1 the null
## hypothesis is rejected. The p-value of the test is returned in @var{pval}.
## A 100(1-alpha)% Clopper-Pearson confidence interval for @var{p} is returned
## in @var{ci}, one-sided for a one-sided test.  @var{stats} is a structure
## with the following fields:
##
## @multitable @columnfractions 0.2 0.75
## @item @qcode{phat} @tab the estimated probability, @code{@var{pos} / @var{N}}.
## @item @qcode{CohensH} @tab Cohen's @math{h}, the difference between the
## arcsine transforms @math{2 asin (sqrt (phat)) - 2 asin (sqrt (p0))}.
## @item @qcode{CohensHCI} @tab its confidence interval, @var{ci} under the same
## transform, which is monotone and so keeps its coverage.
## @end multitable
##
## @end deftypefn

function [h, p, ci, stats] = binotest (pos, n, p0, varargin)

  ## Set default arguments
  alpha = 0.05;
  tail  = 'both';

  i = 1;
  while (i <= length (varargin))
    switch (lower (varargin{i}))
      case 'alpha'
        i = i + 1;
        alpha = varargin{i};
      case 'tail'
        i = i + 1;
        tail = varargin{i};
      otherwise
        error ("binotest: Invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile

  if (! (isnumeric (alpha) && isscalar (alpha) && isreal (alpha)
         && alpha > 0 && alpha < 1))
    error ("binotest: invalid value for alpha.");
  endif
  if (! isa (tail, 'char'))
    error ("binotest: tail must be a string.");
  endif

  if (n <= 0)
    error ("binotest: required n > 0.");
  endif
  if (p0 < 0) || (p0 > 1)
    error ("binotest: required 0 <= p0 <= 1.");
  endif
  if (pos < 0) || (pos > n)
    error ("binotest: required 0 <= pos <= n.");
  endif

  ## Based on the "tail" argument determine the P-value, the critical values,
  ## and the confidence interval.
  switch lower (tail)
    case 'both'
      A_low = binoinv (alpha / 2, n, p0) / n;
      A_high = binoinv (1 - alpha / 2, n, p0) / n;
      p_pos = binopdf (pos, n, p0);
      p_all = binopdf ([0:n], n, p0);
      ## Outcomes as likely as the one observed tie to within rounding
      ind = find (p_all <= p_pos * (1 + 1e-7));
      p = min (1, sum (p_all(ind)));
      if (pos == 0)
        p_low = 0;
      else
        p_low = fzero (@(pl) 1 - binocdf (pos - 1, n, pl) - alpha / 2, [0, 1]);
      endif
      if (pos == n)
        p_high = 1;
      else
        p_high = fzero (@(ph) binocdf (pos, n, ph) - alpha / 2, [0, 1]);
      endif
      ci = [p_low, p_high];
    case 'left'
      p = binocdf (pos, n, p0);
      if (pos == n)
        p_high = 1;
      else
        p_high = fzero (@(ph) binocdf (pos, n, ph) - alpha, [0, 1]);
      endif
      ci = [0, p_high];
    case 'right'
      p = binocdf (pos - 1, n, p0, 'upper');
      if (pos == 0)
        p_low = 0;
      else
        p_low = fzero (@(pl) 1 - binocdf (pos - 1, n, pl) - alpha, [0, 1]);
      endif
      ci = [p_low 1];
    otherwise
      error ("binotest: invalid fifth (tail) argument to binotest.");
  endswitch

  ## Determine the test outcome
  ## MATLAB returns this a double instead of a logical array
  h = double (p < alpha);

  ## Estimate, and the effect size with its interval
  stats.phat = pos / n;
  stats.CohensH = 2 * asin (sqrt (stats.phat)) - 2 * asin (sqrt (p0));
  stats.CohensHCI = 2 * asin (sqrt (ci)) - 2 * asin (sqrt (p0));

endfunction

%!demo
%! % flip a coin 1000 times, showing 475 heads
%! % Hypothesis: coin is fair, i.e. p=1/2
%! [h,p_val,ci] = binotest (475,1000,0.5)
%! % Result: h = 0 : null hypothesis not rejected, coin could be fair
%! %         P value 0.12, i.e. hypothesis not rejected for alpha up to 12%
%! %         0.444 <= p <= 0.506 with 95% confidence

%!demo
%! % flip a coin 100 times, showing 65 heads
%! % Alternative: coin shows more heads than tails, i.e. p>1/2
%! [h,p_val,ci] = binotest (65,100,0.5,'tail','right','alpha',0.01)
%! % Result: h = 1 : null hypothesis is rejected, i.e. coin shows more heads than tails
%! %         P value 0.0018, i.e. hypothesis not rejected for alpha up to 0.18%
%! %         0.53 <= p <= 1 with 99% confidence

%!test #example from https://en.wikipedia.org/wiki/Binomial_test
%! [h,p_val,ci] = binotest (51,235,1/6);
%! assert_equal (p_val, 0.0437, 0.00005)
%! [h,p_val,ci] = binotest (51,235,1/6,'tail','right');
%! assert_equal (p_val, 0.027, 0.0005)
%!test
%! [~, p] = binotest (51, 235, 1/6, 'tail', 'right');
%! assert_equal (p, sum (binopdf (51:235, 235, 1/6)), -1e-12);
%!test
%! [~, p] = binotest (51, 235, 1/6, 'tail', 'left');
%! assert_equal (p, sum (binopdf (0:51, 235, 1/6)), -1e-12);
%!test
%! ## 95 of 100 is no evidence that p < 0.2
%! [h, p] = binotest (95, 100, 0.2, 'tail', 'left');
%! assert_equal ([h, p], [0, 1]);
%!test
%! ## Below the resolution of 1 - binocdf
%! [~, p] = binotest (95, 100, 0.2, 'tail', 'right');
%! assert_equal (p, sum (binopdf (95:100, 100, 0.2)), -1e-12);
%!test
%! [~, ~, ci] = binotest (51, 235, 1/6, 'tail', 'right');
%! assert_equal (ci(2), 1);
%!test
%! [~, ~, ci] = binotest (51, 235, 1/6, 'tail', 'left');
%! assert_equal (ci(1), 0);
%!test
%! ## Outcomes as likely as the one observed are counted, ties to within
%! ## rounding included; values from R's binom.test
%! [~, p] = binotest (1, 10, 0.5);
%! assert_equal (p, 0.021484375, -1e-12);
%!test
%! [~, p] = binotest (4, 10, 0.5);
%! assert_equal (p, 0.75390625, -1e-12);
%!test
%! [~, p] = binotest (11, 50, 0.5);
%! assert_equal (p, 9.021490107130641e-05, -1e-10);
%!test
%! [~, ~, ~, st] = binotest (51, 235, 1/6);
%! assert_equal (st.phat, 51 / 235);
%!test
%! [~, ~, ~, st] = binotest (51, 235, 1/6);
%! assert_equal (st.CohensH, ...
%!               2 * asin (sqrt (51 / 235)) - 2 * asin (sqrt (1/6)), -1e-14);
%!test
%! [~, ~, ci, st] = binotest (51, 235, 1/6);
%! assert_equal (st.CohensHCI, 2 * asin (sqrt (ci)) - 2 * asin (sqrt (1/6)), ...
%!               -1e-14);
%!test
%! ## A one-sided interval maps to a one-sided interval
%! [~, ~, ~, st] = binotest (51, 235, 1/6, 'tail', 'right');
%! assert_equal (st.CohensHCI(2), pi - 2 * asin (sqrt (1/6)), -1e-14);

%!error<binotest: Invalid Name argument.> binotest (5, 10, 0.5, 'size', 1)
%!error<binotest: invalid value for alpha.> binotest (5, 10, 0.5, 'alpha', 0)
%!error<binotest: invalid value for alpha.> binotest (5, 10, 0.5, 'alpha', 1)
%!error<binotest: invalid value for alpha.> ...
%! binotest (5, 10, 0.5, 'alpha', [0.05, 0.1])
%!error<binotest: tail must be a string.> binotest (5, 10, 0.5, 'tail', 1)
%!error<binotest: required n > 0.> binotest (0, 0, 0.5)
%!error<binotest: required 0 <= p0 <= 1.> binotest (5, 10, 1.5)
%!error<binotest: required 0 <= pos <= n.> binotest (11, 10, 0.5)
%!error<binotest: invalid fifth \(tail\) argument to binotest.> ...
%! binotest (5, 10, 0.5, 'tail', 'up')
