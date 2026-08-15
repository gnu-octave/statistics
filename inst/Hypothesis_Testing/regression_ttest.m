## Copyright (C) 1995-2017 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{h} =} regression_ttest (@var{y}, @var{x})
## @deftypefnx {statistics} {[@var{h}, @var{pval}] =} regression_ttest (@var{y}, @var{x})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{ci}] =} regression_ttest (@var{y}, @var{x})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} regression_ttest (@var{y}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} regression_ttest (@var{y}, @var{x}, @var{Name}, @var{Value})
##
## Perform a linear regression t-test.
##
## @code{@var{h} = regression_ttest (@var{y}, @var{x})} tests the null
## hypothesis that the slope @math{beta1} of a simple linear regression equals
## 0.  The result is @var{h} = 0 if the null hypothesis cannot be rejected at
## the 5% significance level, or @var{h} = 1 if the null hypothesis can be
## rejected at the 5% level.  @var{y} and @var{x} must be vectors of equal
## length with finite real numbers.
##
## The p-value of the test is returned in @var{pval}.  A @math{100(1-alpha)%}
## confidence interval for @math{beta1} is returned in @var{ci}.  @var{stats} is
## a structure containing the value of the test statistic (@qcode{tstat}),
## the degrees of freedom (@qcode{df}), the slope coefficient (@qcode{beta1}),
## and the intercept (@qcode{beta0}).  Under the null, the test statistic
## @var{stats}.@qcode{tstat} follows a @math{T}-distribution with
## @var{stats}.@qcode{df} degrees of freedom.
##
## @code{[@dots{}] = regression_ttest (@dots{}, @var{name}, @var{value})}
## specifies one or more of the following name/value pairs:
##
## @multitable @columnfractions 0.2 0.75
## @headitem Name @tab Value
## @item @qcode{'alpha'} @tab the significance level. Default is 0.05.
##
## @item @qcode{'tail'} @tab a string specifying the alternative hypothesis
## @end multitable
## @multitable @columnfractions 0.25 0.65
## @item @qcode{'both'} @tab @math{beta1} is not 0 (two-tailed, default)
## @item @qcode{'left'} @tab @math{beta1} is less than 0 (left-tailed)
## @item @qcode{'right'} @tab @math{beta1} is greater than 0 (right-tailed)
## @end multitable
##
## @seealso{regression_ftest, regress, regress_gp}
## @end deftypefn

function [h, pval, ci, stats] = regression_ttest (y, x, varargin)

  ## Check for valid input
  if (nargin < 2)
    print_usage ();
  endif

  ## Check for finite real numbers in Y, X
  if (! all (isfinite (y)) || ! isreal (y))
    error ("regression_ttest: Y must contain finite real numbers.");
  endif
  if (! all (isfinite (x(:))) || ! isreal (x))
    error ("regression_ttest: X must contain finite real numbers.");
  endif

  # Get number of observations
  n = length (y);

  ## Check Y and X have the same number of observations
  if (! isvector (y) || ! isvector (x) || length (x) != n)
    error ("regression_ttest: Y and X must be vectors of equal length.");
  endif


  ## Set default arguments
  alpha = 0.05;
  tail = 'both';

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
          error ("regression_ttest: invalid value for alpha.");
        endif
      case 'tail'
        i = i + 1;
        tail = varargin{i};
        if (! any (strcmpi (tail, {'both', 'left', 'right'})))
          error ("regression_ttest: invalid value for tail.");
        endif
      otherwise
        error ("regression_ttest: invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile

  y_bar = mean (y);
  x_bar = mean (x);

  ## Ordinary least squares slope, taken from the centred data.  Going through
  ## cov (x, y) tied this to whether cov returns the scalar covariance or the
  ## 2*2 matrix; it returns the matrix, which left every call form of this
  ## function dying on a nonconformant operand.
  x_ctr = x - x_bar;
  Sxx = sum (x_ctr .^ 2);
  stats.beta1 = sum (x_ctr .* (y - y_bar)) / Sxx;
  stats.beta0 = y_bar - stats.beta1 * x_bar;

  ## Fitted values, and from them the residual sum of squares.  Evaluating the
  ## fit at x_bar rather than at x returned y_bar for every observation, since
  ## the least squares line passes through the means, so SSE was really the
  ## total sum of squares and the test statistic came out far too small.
  y_hat = stats.beta0 + stats.beta1 * x;
  SSE = sum ((y - y_hat) .^ 2);
  stats.df = n - 2;
  SE = sqrt (SSE / stats.df);

  term = SE / sqrt (Sxx);
  stats.tstat = stats.beta1 / term;

  ## Based on the "tail" argument determine the P-value, the critical values,
  ## and the confidence interval.
  switch lower (tail)
    case 'both'
      pval = 2 * (1 - tcdf (abs (stats.tstat), stats.df));
      tcrit = - tinv (alpha / 2, stats.df);
      ci = [stats.beta1 - tcrit * term; stats.beta1 + tcrit * term];
    case 'left'
      pval = tcdf (stats.tstat, stats.df);
      tcrit = - tinv (alpha, stats.df);
      ci = [-inf; stats.beta1 + tcrit * term];
    case 'right'
      pval = 1 - tcdf (stats.tstat, stats.df);
      tcrit = - tinv (alpha, stats.df);
      ci = [stats.beta1 - tcrit * term; inf];
  endswitch

  ## Determine the test outcome
  h = double (pval < alpha);
  h(isnan (pval)) = NaN;

endfunction

## The function had no test that called it successfully, which is how every
## documented call form came to die on a nonconformant operand.  Check the
## returned quantities against an ordinary least squares fit done here.
%!test
%! x = [1 2 3 4 5 6 7 8 9 10]';
%! y = 2 + 3 * x + [0.3 -0.2 0.5 -0.4 0.1 0.2 -0.3 0.4 -0.1 0.2]';
%! [h, pval, ci, stats] = regression_ttest (y, x);
%! n = numel (x);
%! b = [ones(n,1), x] \ y;
%! resid = y - [ones(n,1), x] * b;
%! se = sqrt (sum (resid .^ 2) / (n - 2) / sum ((x - mean (x)) .^ 2));
%! assert_equal (stats.beta1, b(2), 1e-12);
%! assert_equal (stats.beta0, b(1), 1e-12);
%! assert_equal (stats.df, n - 2);
%! assert_equal (stats.tstat, b(2) / se, 1e-9);
%! assert_equal (pval, 2 * (1 - tcdf (abs (b(2) / se), n - 2)), 1e-12);
%! assert_equal (ci(:)', [b(2) - tinv(0.975, n-2) * se, ...
%!                        b(2) + tinv(0.975, n-2) * se], 1e-9);
%! assert_equal (h, 1);

## The fitted values must come from x, not from its mean.  Evaluating the fit
## at x_bar returns y_bar for every observation, because the least squares
## line passes through the means, which made SSE the total sum of squares.
%!test
%! x = [1 2 3 4 5 6 7 8 9 10]';
%! y = 2 + 3 * x + [0.3 -0.2 0.5 -0.4 0.1 0.2 -0.3 0.4 -0.1 0.2]';
%! [~, ~, ~, stats] = regression_ttest (y, x);
%! SSE = sum ((y - (stats.beta0 + stats.beta1 * x)) .^ 2);
%! SST = sum ((y - mean (y)) .^ 2);
%! assert_equal (stats.tstat, ...
%!               stats.beta1 / sqrt (SSE / stats.df / sum ((x - mean (x)) .^ 2)), ...
%!               1e-9);
%! assert_equal (SSE < SST / 100, true);

## A slope of zero must not be rejected.
%!test
%! x = (1:40)';
%! y = [0.5 -1.2 0.3 0.8 -0.4 1.1 -0.7 0.2 -0.9 0.6 ...
%!      1.3 -0.1 0.4 -1.1 0.7 0.9 -0.5 0.1 -0.8 1.0 ...
%!      -0.3 0.2 1.2 -0.6 0.5 -1.0 0.8 0.3 -0.2 0.9 ...
%!      -1.3 0.6 0.1 -0.7 1.1 -0.4 0.7 0.2 -0.9 0.4]';
%! [h, pval] = regression_ttest (y, x);
%! assert_equal (h, 0);
%! assert_equal (pval > 0.05, true);

## The one-sided tails must bracket the two-sided p-value.
%!test
%! x = [1 2 3 4 5 6 7 8 9 10]';
%! y = 2 + 3 * x + [0.3 -0.2 0.5 -0.4 0.1 0.2 -0.3 0.4 -0.1 0.2]';
%! [~, pb] = regression_ttest (y, x, 'tail', 'both');
%! [~, pr] = regression_ttest (y, x, 'tail', 'right');
%! [~, pl] = regression_ttest (y, x, 'tail', 'left');
%! assert_equal (pr, pb / 2, 1e-12);
%! assert_equal (pl, 1 - pr, 1e-12);

## A row-vector input must give the same answer as a column.
%!test
%! x = [1 2 3 4 5 6 7 8 9 10]';
%! y = 2 + 3 * x + [0.3 -0.2 0.5 -0.4 0.1 0.2 -0.3 0.4 -0.1 0.2]';
%! [~, ~, ~, sc] = regression_ttest (y, x);
%! [~, ~, ~, sr] = regression_ttest (y', x');
%! assert_equal (sr.beta1, sc.beta1, 1e-12);
%! assert_equal (sr.tstat, sc.tstat, 1e-9);

## Test input validation
%!error<Invalid call to regression_ttest.  Correct usage> regression_ttest ();
%!error<Invalid call to regression_ttest.  Correct usage> regression_ttest (1);
%!error<regression_ttest: Y must contain finite real numbers.> ...
%! regression_ttest ([1 2 NaN]', [2 3 4]');
%!error<regression_ttest: Y must contain finite real numbers.> ...
%! regression_ttest ([1 2 Inf]', [2 3 4]');
%!error<regression_ttest: Y must contain finite real numbers.> ...
%! regression_ttest ([1 2 3+i]', [2 3 4]');
%!error<regression_ttest: X must contain finite real numbers.> ...
%! regression_ttest ([1 2 3]', [2 3 NaN]');
%!error<regression_ttest: X must contain finite real numbers.> ...
%! regression_ttest ([1 2 3]', [2 3 Inf]');
%!error<regression_ttest: X must contain finite real numbers.> ...
%! regression_ttest ([1 2 3]', [3 4 3+i]');
%!error<regression_ttest: Y and X must be vectors of equal length.> ...
%! regression_ttest ([1 2 3]', [3 4 4 5]');
%!error<regression_ttest: invalid value for alpha.> ...
%! regression_ttest ([1 2 3]', [2 3 4]', 'alpha', 0);
%!error<regression_ttest: invalid value for alpha.> ...
%! regression_ttest ([1 2 3]', [2 3 4]', 'alpha', 1.2);
%!error<regression_ttest: invalid value for alpha.> ...
%! regression_ttest ([1 2 3]', [2 3 4]', 'alpha', [.02 .1]);
%!error<regression_ttest: invalid value for alpha.> ...
%! regression_ttest ([1 2 3]', [2 3 4]', 'alpha', 'a');
%!error<regression_ttest: invalid Name argument.> ...
%! regression_ttest ([1 2 3]', [2 3 4]', 'some', 0.05);
%!error<regression_ttest: invalid value for tail.>  ...
%! regression_ttest ([1 2 3]', [2 3 4]', 'tail', 'val');
%!error<regression_ttest: invalid value for tail.>  ...
%! regression_ttest ([1 2 3]', [2 3 4]', 'alpha', 0.01, 'tail', 'val');
