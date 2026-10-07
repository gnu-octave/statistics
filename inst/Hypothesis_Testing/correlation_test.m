## Copyright (C) 1995-2017 Kurt Hornik
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
## @deftypefn  {statistics} {@var{h} =} correlation_test (@var{x}, @var{y})
## @deftypefnx {statistics} {[@var{h}, @var{p}] =} correlation_test (@var{y}, @var{x})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{stats}] =} correlation_test (@var{y}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} correlation_test (@var{y}, @var{x}, @var{Name}, @var{Value})
##
## Perform a correlation coefficient test to determine whether two samples
## @var{x} and @var{y} come from uncorrelated populations.
##
## @code{@var{h} = correlation_test (@var{y}, @var{x})} tests the null
## hypothesis that the two samples @var{x} and @var{y} come from uncorrelated
## populations.  The result is @var{h} = 0 if the null hypothesis cannot be
## rejected at the significance level @var{alpha}, or @var{h} = 1 if it can.
## @var{y} and @var{x} must be vectors of equal length with finite real
## numbers.
##
## The p-value of the test is returned in @var{p}.  @var{stats} is a
## structure with the following fields:
## @multitable @columnfractions 0.2 0.70
## @headitem Field @tab Value
## @item @qcode{method} @tab the type of correlation coefficient used
## for the test
## @item @qcode{tstat} @tab the t statistic, for Pearson's and Spearman's
## coefficients
## @item @qcode{df} @tab its degrees of freedom, @math{n - 2}
## @item @qcode{zval} @tab the z statistic, for Kendall's coefficient
## @item @qcode{dist} @tab the respective distribution for the test
## @item @qcode{alt} @tab the alternative hypothesis for the test
## @item @qcode{CorrCoef} @tab the correlation coefficient, the effect size
## @item @qcode{CorrCoefCI} @tab its @math{100 (1 - alpha)}% confidence
## interval, one-sided for a one-sided test
## @end multitable
##
## Pearson's and Spearman's coefficients are tested with
## @math{t = r sqrt ((n - 2) / (1 - r^2))} on @math{n - 2} degrees of freedom,
## as R's @code{cor.test} does without an exact test; Kendall's with its
## normal approximation, whose variance allows for ties.  The confidence
## interval is Fisher's @math{z} transform of the coefficient, with standard
## error @math{1 / sqrt (n - 3)} for Pearson's,
## @math{sqrt ((1 + r^2 / 2) / (n - 3))} for Spearman's and
## @math{sqrt (0.437 / (n - 4))} for Kendall's, after Bonett and Wright (2000).
##
##
## @code{[@dots{}] = correlation_test (@dots{}, @var{name}, @var{value})}
## specifies one or more of the following name/value pairs:
##
## @multitable @columnfractions 0.2 0.75
## @headitem Name @tab Value
## @item @qcode{'alpha'} @tab the significance level of @var{h} and the level
## of the confidence interval. Default is 0.05.
##
## @item @qcode{'tail'} @tab a string specifying the alternative hypothesis
## @end multitable
## @multitable @columnfractions 0.25 0.65
## @item @qcode{'both'} @tab @math{corrcoef} is not 0 (two-tailed, default)
## @item @qcode{'left'} @tab @math{corrcoef} is less than 0 (left-tailed)
## @item @qcode{'right'} @tab @math{corrcoef} is greater than 0
## (right-tailed)
## @end multitable
##
## @multitable @columnfractions 0.2 0.75
## @item @qcode{'method'} @tab a string specifying the correlation
## coefficient used for the test
## @end multitable
## @multitable @columnfractions 0.25 0.65
## @item @qcode{'pearson'} @tab Pearson's product moment correlation
## (Default)
## @item @qcode{'kendall'} @tab Kendall's rank correlation tau
## @item @qcode{'spearman'} @tab Spearman's rank correlation rho
## @end multitable
##
## @seealso{regression_ftest, regression_ttest}
## @end deftypefn

function [h, p, stats] = correlation_test (x, y, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  if (! isvector (x) || ! isvector (y) || length (x) != length (y))
    error ("correlation_test: X and Y must be vectors of equal length.");
  endif

  ## Force to column vectors
  x = x(:);
  y = y(:);

  ## Check for finite real numbers in X and Y
  if (! all (isfinite (x)) || ! isreal (x))
    error ("correlation_test: X must contain finite real numbers.");
  endif
  if (! all (isfinite (y(:))) || ! isreal (y))
    error ("correlation_test: Y must contain finite real numbers.");
  endif

  ## Set default arguments
  alpha = 0.05;
  tail = 'both';
  method = 'pearson';

  ## Check additional options
  i = 1;
  while (i <= length (varargin))
    switch lower (varargin{i})
      case 'alpha'
        i = i + 1;
        alpha = varargin{i};
        ## Check for valid alpha
        if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                    alpha <= 0 || alpha >= 1)
          error ("correlation_test: invalid value for alpha.");
        endif
      case 'tail'
        i = i + 1;
        tail = varargin{i};
        if (! any (strcmpi (tail, {'both', 'left', 'right'})))
          error ("correlation_test: invalid value for tail.");
        endif
      case 'method'
        i = i + 1;
        method = varargin{i};
        if (! any (strcmpi (method, {'pearson', 'kendall', 'spearman'})))
          error ("correlation_test: invalid value for method.");
        endif
      otherwise
        error ("correlation_test: invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile

  n = length (x);
  method = lower (method);
  tail = lower (tail);

  if (strcmp (method, 'pearson'))
    r = corr (x, y);
    stats.method = 'Pearson''s product moment correlation';
  elseif (strcmp (method, 'kendall'))
    r = kendall (x, y);
    stats.method = 'Kendall''s rank correlation tau';
  else  # spearman
    r = spearman (x, y);
    stats.method = 'Spearman''s rank correlation rho';
  endif

  if (strcmp (method, 'kendall'))
    ## The concordance score and its variance, allowing for ties in either
    ## sample (Kendall, 1970)
    S = 0;
    for i = 1:n-1
      S += sum (sign (x(i) - x(i+1:n)) .* sign (y(i) - y(i+1:n)));
    endfor
    [~, ~, jx] = unique (x);
    [~, ~, jy] = unique (y);
    t = accumarray (jx, 1);
    u = accumarray (jy, 1);
    v = (n * (n - 1) * (2 * n + 5) - sum (t .* (t - 1) .* (2 * t + 5)) ...
         - sum (u .* (u - 1) .* (2 * u + 5))) / 18 ...
        + sum (t .* (t - 1) .* (t - 2)) * sum (u .* (u - 1) .* (u - 2)) ...
          / (9 * n * (n - 1) * (n - 2)) ...
        + sum (t .* (t - 1)) * sum (u .* (u - 1)) / (2 * n * (n - 1));
    stats.zval = S / sqrt (v);
    stats.dist = 'standard normal';
    lower_tail = normcdf (stats.zval);
    upper_tail = normcdf (-stats.zval);
    se = sqrt (0.437 / (n - 4));
  else
    stats.tstat = r * sqrt ((n - 2) / (1 - r ^ 2));
    stats.df = n - 2;
    stats.dist = 'Student''s t';
    lower_tail = tcdf (stats.tstat, stats.df);
    upper_tail = tcdf (stats.tstat, stats.df, 'upper');
    if (strcmp (method, 'pearson'))
      se = 1 / sqrt (n - 3);
    else
      se = sqrt ((1 + r ^ 2 / 2) / (n - 3));
    endif
  endif

  ## Based on the "tail" argument determine the P-value and the quantile of
  ## the confidence interval
  switch (tail)
    case 'both'
      p = 2 * min (lower_tail, upper_tail);
      z = norminv (1 - alpha / 2);
    case 'right'
      p = upper_tail;
      z = norminv (1 - alpha);
    case 'left'
      p = lower_tail;
      z = norminv (1 - alpha);
  endswitch

  stats.alt = tail;

  ## The coefficient is the effect size; its interval by Fisher's z, one-sided
  ## for a one-sided test
  stats.CorrCoef = r;
  if (isreal (se) && se > 0)
    stats.CorrCoefCI = tanh (atanh (r) + [-1, 1] * z * se);
  else
    stats.CorrCoefCI = [NaN, NaN];
  endif
  if (strcmp (tail, 'right'))
    stats.CorrCoefCI(2) = 1;
  elseif (strcmp (tail, 'left'))
    stats.CorrCoefCI(1) = -1;
  endif

  ## Determine the test outcome
  h = double (p < alpha);

endfunction

%!test
%! x = [6 7 7 9 10 12 13 14 15 17];
%! y = [19 22 27 25 30 28 30 29 25 32];
%! [h, p, stats] = correlation_test (x, y);
%! assert_equal (stats.CorrCoef, corr (x', y'), 1e-14);
%! assert_equal (p, 0.0223, 1e-4);
%!test
%! x = [6 7 7 9 10 12 13 14 15 17]';
%! y = [19 22 27 25 30 28 30 29 25 32]';
%! [h, p, stats] = correlation_test (x, y);
%! assert_equal (stats.CorrCoef, corr (x, y), 1e-14);
%! assert_equal (p, 0.0223, 1e-4);
%!test
%! ## Below the resolution of 1 - tcdf, the p-value of corr in MATLAB R2024a
%! x = (1:30)';
%! [~, p] = correlation_test (x, x + 0.01 * sin (x));
%! assert_equal (p, 6.68262089195668e-88, -1e-7);
%!test
%! ## Values from R's cor.test without exact tests, on tied data
%! x = [1 2 2 3 4 4 4 5 6 7 7 8]';
%! y = [2 1 3 3 5 4 6 6 8 7 9 9]';
%! [~, ~, stats] = correlation_test (x, y);
%! assert_equal (stats.CorrCoefCI, ...
%!               [0.8116031451285135, 0.98487109489995917], -1e-12);
%!test
%! x = [1 2 2 3 4 4 4 5 6 7 7 8]';
%! y = [2 1 3 3 5 4 6 6 8 7 9 9]';
%! [~, p, stats] = correlation_test (x, y, 'tail', 'right');
%! assert_equal ([p, stats.CorrCoefCI], ...
%!               [1.7690247076465537e-06, 0.84452474586233817, 1], -1e-12);
%!test
%! ## Kendall's variance allows for ties
%! x = [1 2 2 3 4 4 4 5 6 7 7 8]';
%! y = [2 1 3 3 5 4 6 6 8 7 9 9]';
%! [~, p, stats] = correlation_test (x, y, 'method', 'kendall');
%! assert_equal ([stats.zval, p], ...
%!               [3.7786519487436312, 0.00015767962747712322], -1e-12);
%!test
%! x = [1 2 2 3 4 4 4 5 6 7 7 8]';
%! y = [2 1 3 3 5 4 6 6 8 7 9 9]';
%! [~, p] = correlation_test (x, y, 'method', 'spearman');
%! assert_equal (p, 1.2599089155461711e-06, -1e-12);
%!test
%! ## Spearman's test is symmetric in the sign of the coefficient
%! x = (1:12)';
%! y = [3 1 4 2 6 5 8 9 7 12 10 11]';
%! [~, p] = correlation_test (x, -y, 'method', 'spearman');
%! assert_equal (p, 2.8428045348547355e-05, -1e-12);
%!test
%! x = (1:12)';
%! y = [3 1 4 2 6 5 8 9 7 12 10 11]';
%! [~, ~, stats] = correlation_test (x, y, 'method', 'spearman');
%! r = stats.CorrCoef;
%! assert_equal (stats.CorrCoefCI, ...
%!               tanh (atanh (r) + [-1, 1] * norminv (0.975) ...
%!                     * sqrt ((1 + r ^ 2 / 2) / 9)), -1e-14);
%!test
%! x = (1:12)';
%! y = [3 1 4 2 6 5 8 9 7 12 10 11]';
%! [~, ~, stats] = correlation_test (x, y, 'method', 'kendall', 'alpha', 0.01);
%! r = stats.CorrCoef;
%! assert_equal (stats.CorrCoefCI, ...
%!               tanh (atanh (r) + [-1, 1] * norminv (0.995) ...
%!                     * sqrt (0.437 / 8)), -1e-14);
%!test
%! [~, ~, stats] = correlation_test ((1:12)', (12:-1:1)' + sin (1:12)', ...
%!                                   'tail', 'left');
%! assert_equal (stats.CorrCoefCI(1), -1);
%!test
%! [~, ~, stats] = correlation_test ((1:12)', sin (1:12)', 'method', 'kendall');
%! assert_equal ([isfield(stats, 'zval'), isfield(stats, 'tstat')], ...
%!               [true, false]);

## Test input validation
%!error<Invalid call to correlation_test.  Correct usage> correlation_test ();
%!error<Invalid call to correlation_test.  Correct usage> correlation_test (1);
%!error<correlation_test: X must contain finite real numbers.> ...
%! correlation_test ([1 2 NaN]', [2 3 4]');
%!error<correlation_test: X must contain finite real numbers.> ...
%! correlation_test ([1 2 Inf]', [2 3 4]');
%!error<correlation_test: X must contain finite real numbers.> ...
%! correlation_test ([1 2 3+i]', [2 3 4]');
%!error<correlation_test: Y must contain finite real numbers.> ...
%! correlation_test ([1 2 3]', [2 3 NaN]');
%!error<correlation_test: Y must contain finite real numbers.> ...
%! correlation_test ([1 2 3]', [2 3 Inf]');
%!error<correlation_test: Y must contain finite real numbers.> ...
%! correlation_test ([1 2 3]', [3 4 3+i]');
%!error<correlation_test: X and Y must be vectors of equal length.> ...
%! correlation_test ([1 2 3]', [3 4 4 5]');
%!error<correlation_test: invalid value for alpha.> ...
%! correlation_test ([1 2 3]', [2 3 4]', 'alpha', 0);
%!error<correlation_test: invalid value for alpha.> ...
%! correlation_test ([1 2 3]', [2 3 4]', 'alpha', 1.2);
%!error<correlation_test: invalid value for alpha.> ...
%! correlation_test ([1 2 3]', [2 3 4]', 'alpha', [.02 .1]);
%!error<correlation_test: invalid value for alpha.> ...
%! correlation_test ([1 2 3]', [2 3 4]', 'alpha', 'a');
%!error<correlation_test: invalid Name argument.> ...
%! correlation_test ([1 2 3]', [2 3 4]', 'some', 0.05);
%!error<correlation_test: invalid value for tail.>  ...
%! correlation_test ([1 2 3]', [2 3 4]', 'tail', 'val');
%!error<correlation_test: invalid value for tail.>  ...
%! correlation_test ([1 2 3]', [2 3 4]', 'alpha', 0.01, 'tail', 'val');
%!error<correlation_test: invalid value for method.>  ...
%! correlation_test ([1 2 3]', [2 3 4]', 'method', 0.01);
%!error<correlation_test: invalid value for method.>  ...
%! correlation_test ([1 2 3]', [2 3 4]', 'method', 'some');
