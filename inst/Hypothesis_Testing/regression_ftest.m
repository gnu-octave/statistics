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
## @deftypefn  {statistics} {@var{h} =} regression_ftest (@var{y}, @var{x})
## @deftypefnx {statistics} {@var{h} =} regression_ftest (@var{y}, @var{x}, @var{keep})
## @deftypefnx {statistics} {@var{h} =} regression_ftest (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{h}, @var{p}, @var{stats}] =} regression_ftest (@dots{})
##
## F-test for nested linear regression models.
##
## Fit the full linear model of the response @var{y} on every column of
## @var{x}, and a reduced model on the columns of @var{x} indexed by
## @var{keep}, both by least squares, and test the null hypothesis that the
## columns the reduced model leaves out add nothing to it.  @var{keep} holds
## column indices or a logical mask; left out or empty, the reduced model has
## no predictors and the test is of every column of @var{x} together.  Both
## models include a constant term unless @qcode{'Intercept'} is false.
## @var{y} must be a vector with one element per row of @var{x}, and neither
## may contain missing values.
##
## The test decision is returned in @var{h}, which is 1 when the null
## hypothesis is rejected at the significance level @var{alpha} and 0
## otherwise, and the p-value in @var{p}.  @var{stats} is a structure with the
## following fields:
##
## @multitable @columnfractions 0.2 0.75
## @item @qcode{fstat} @tab the F statistic.
## @item @qcode{df1} @tab its numerator degrees of freedom, the number of
## coefficients the reduced model leaves out.
## @item @qcode{df2} @tab its denominator degrees of freedom, the residual
## degrees of freedom of the full model.
## @item @qcode{betafull} @tab the coefficients of the full model, the constant
## first where there is one.
## @item @qcode{betareduced} @tab those of the reduced model.
## @item @qcode{ssefull} @tab the residual sum of squares of the full model.
## @item @qcode{ssereduced} @tab that of the reduced model.
## @item @qcode{CohensF2} @tab Cohen's @math{f^2}, the effect size of the
## columns tested.
## @item @qcode{CohensF2CI} @tab its @math{100 (1 - alpha)}% confidence
## interval.
## @end multitable
##
## The degrees of freedom come from the ranks of the two designs.  Cohen's
## @math{f^2} is taken, as Cohen defines it, from the noncentrality of the F
## statistic, @math{lambda = f^2 (df1 + df2 + 1)}: the point estimate from
## @math{max (F df1 (df2 - 2) / df2 - df1, 0)}, unbiased for the noncentrality
## before its truncation at zero, and the interval by inverting the noncentral
## F distribution.  Where the interval would need a noncentrality above
## @math{10^5}, where @code{ncfcdf} loses accuracy, it is @qcode{NaN}.
##
## The following @var{Name}/@var{Value} pairs are accepted:
##
## @multitable @columnfractions 0.2 0.75
## @item @qcode{'Alpha'} @tab the significance level of @var{h} and the level
## of the confidence interval.  Default is 0.05.
## @item @qcode{'Intercept'} @tab whether both models include a constant term.
## Default is true.
## @end multitable
##
## The same tests, and any linear hypothesis on the coefficients, are available
## from a fitted @code{LinearModel} through its @code{coefTest} method.
##
## @seealso{regression_ttest, regress, fitlm}
## @end deftypefn

function [h, p, stats] = regression_ftest (y, x, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  ## Check for finite real numbers in Y, X
  if (! all (isfinite (y(:))) || ! isreal (y))
    error ("regression_ftest: Y must contain finite real numbers.");
  endif
  if (! all (isfinite (x(:))) || ! isreal (x))
    error ("regression_ftest: X must contain finite real numbers.");
  endif
  [n, v] = size (x);
  if (! (isvector (y) && numel (y) == n))
    error ("regression_ftest: Y must be a vector of length 'rows (X)'.");
  endif
  y = y(:);

  ## The columns of the reduced model, then the options
  keep = [];
  if (! isempty (varargin) && ! ischar (varargin{1}))
    keep = varargin{1};
    varargin(1) = [];
  endif
  alpha = 0.05;
  intercept = true;
  if (mod (numel (varargin), 2) != 0)
    error ("regression_ftest: optional arguments must be in pairs.");
  endif
  for i = 1:2:numel (varargin)
    name = varargin{i};
    value = varargin{i+1};
    if (! ischar (name))
      error ("regression_ftest: invalid Name argument.");
    endif
    switch (lower (name))
      case 'alpha'
        if (! (isnumeric (value) && isscalar (value) && isreal (value)
               && value > 0 && value < 1))
          error ("regression_ftest: invalid value for alpha.");
        endif
        alpha = value;
      case 'intercept'
        if (! (isscalar (value) && (islogical (value) || isnumeric (value))))
          error ("regression_ftest: invalid value for Intercept.");
        endif
        intercept = logical (value);
      otherwise
        error ("regression_ftest: invalid Name argument.");
    endswitch
  endfor

  ## Resolve KEEP to column indices of a proper subset of X
  if (islogical (keep))
    if (numel (keep) != v)
      error (strcat ("regression_ftest: a logical KEEP must have one", ...
                     " element per column of X."));
    endif
    keep = find (keep);
  elseif (! (isempty (keep) || (isnumeric (keep) && isvector (keep)
             && all (keep == fix (keep)) && all (keep >= 1 & keep <= v)
             && numel (unique (keep)) == numel (keep))))
    error (strcat ("regression_ftest: KEEP must hold distinct column", ...
                   " indices of X."));
  endif
  if (numel (keep) >= v)
    error ("regression_ftest: KEEP must leave out at least one column of X.");
  endif

  ## Fit both models by least squares
  if (intercept)
    Xf = [ones(n, 1), x];
    Xr = [ones(n, 1), x(:,keep)];
  else
    Xf = x;
    Xr = x(:,keep);
  endif
  bf = Xf \ y;
  sse_f = sumsq (y - Xf * bf);
  if (columns (Xr) > 0)
    br = Xr \ y;
    sse_r = sumsq (y - Xr * br);
  else
    br = zeros (0, 1);
    sse_r = sumsq (y);
  endif
  df1 = rank (Xf) - rank (Xr);
  df2 = n - rank (Xf);
  if (df2 < 1)
    error ("regression_ftest: too few observations for the full model.");
  endif
  if (df1 < 1)
    error (strcat ("regression_ftest: the columns left out add nothing", ...
                   " to the reduced model's design."));
  endif

  stats.fstat = ((sse_r - sse_f) / df1) / (sse_f / df2);
  stats.df1 = df1;
  stats.df2 = df2;
  stats.betafull = bf;
  stats.betareduced = br;
  stats.ssefull = sse_f;
  stats.ssereduced = sse_r;
  p = fcdf (stats.fstat, df1, df2, 'upper');
  h = double (p < alpha);

  ## Cohen's f^2 from the noncentrality, lambda = f^2 (df1 + df2 + 1)
  [lambda, lambdahat] = __ncfbounds__ (stats.fstat, df1, df2, alpha);
  stats.CohensF2 = lambdahat / (df1 + df2 + 1);
  stats.CohensF2CI = lambda / (df1 + df2 + 1);

endfunction

%!shared X, y
%! X = [1 1; 2 1; 3 2; 4 2; 5 3; 6 3; 7 4; 8 4; 9 5; 10 5; 11 6; 12 6];
%! y = 2 + 3 * X(:,1) - 1.5 * X(:,2) + ...
%!     [0.2 -0.3 0.1 0.4 -0.2 0.3 -0.1 0.2 -0.4 0.1 0.3 -0.2]';
%!test
%! ## Every column together, as R's anova (lm (y ~ 1), lm (y ~ x1 + x2))
%! [h, p, st] = regression_ftest (y, X);
%! assert_equal ([h, st.df1, st.df2], [1, 2, 9]);
%! assert_equal ([st.fstat, p], ...
%!               [4553.7218634686296, 2.9845354551931477e-14], -1e-10);
%!test
%! ## Dropping one column, as R's anova (lm (y ~ x1), lm (y ~ x1 + x2))
%! [~, p, st] = regression_ftest (y, X, 1);
%! assert_equal ([st.fstat, p], [27.05295073929755, 0.00056311880681268267], ...
%!               -1e-10);
%!test
%! ## Without a constant, as R's anova (lm (y ~ 0 + x1), lm (y ~ 0 + x1 + x2))
%! [~, p, st] = regression_ftest (y, X, 1, 'Intercept', false);
%! assert_equal ([st.fstat, p], ...
%!               [0.0076424950016706168, 0.93206240862541445], -1e-10);
%!test
%! [~, ~, st] = regression_ftest (y, X, 1, 'Intercept', false);
%! assert_equal ([st.df1, st.df2], [1, 10]);
%!test
%! [~, pa] = regression_ftest (y, X, [true, false]);
%! [~, pb] = regression_ftest (y, X, 1);
%! assert_equal (pa, pb);
%!test
%! [~, ~, st] = regression_ftest (y, X, 1);
%! assert_equal (st.betafull, [ones(12, 1), X] \ y, -1e-12);
%!test
%! [~, ~, st] = regression_ftest (y, X, 1);
%! assert_equal (st.betareduced, [ones(12, 1), X(:,1)] \ y, -1e-12);
%!test
%! [~, ~, st] = regression_ftest (y, X, 1);
%! assert_equal (st.ssefull, sumsq (y - [ones(12, 1), X] * st.betafull), ...
%!               -1e-12);
%!test
%! ## Cohen's f^2 with lambda = f^2 (df1 + df2 + 1); the bounds are R's pf
%! ## with ncp inverted by uniroot
%! [~, ~, st] = regression_ftest (y, X, 1);
%! assert_equal (st.CohensF2, 1.8219258098493329, -1e-12);
%!test
%! [~, ~, st] = regression_ftest (y, X, 1);
%! assert_equal (st.CohensF2CI, [0.39131195887888376, 6.1366401225022624], ...
%!               -1e-8);
%!test
%! ## A column that adds nothing has an interval from zero
%! [h, ~, st] = regression_ftest (X(:,1) + sin (5 * (1:12)'), X, 1);
%! assert_equal ([h, st.CohensF2CI(1)], [0, 0]);
%!test
%! ## 'Alpha' moves the decision but never the p-value
%! [h1, p1] = regression_ftest (y, X, 1, 'Alpha', 1e-4);
%! [h2, p2] = regression_ftest (y, X, 1, 'Alpha', 0.01);
%! assert_equal ([h1, h2, p1], [0, 1, p2]);
%!test
%! [~, pr] = regression_ftest (y', X, 1);
%! [~, pc] = regression_ftest (y, X, 1);
%! assert_equal (pr, pc);
%!test
%! ## Degrees of freedom follow the ranks of the designs
%! [~, ~, st] = regression_ftest (y, [X, X(:,1)]);
%! assert_equal (st.df1, 2);

## Test input validation
%!error<Invalid call to regression_ftest.  Correct usage> regression_ftest ();
%!error<Invalid call to regression_ftest.  Correct usage> ...
%! regression_ftest ([1 2 3]');
%!error<regression_ftest: Y must contain finite real numbers.> ...
%! regression_ftest ([1 2 NaN]', [2 3 4; 3 4 5]');
%!error<regression_ftest: Y must contain finite real numbers.> ...
%! regression_ftest ([1 2 Inf]', [2 3 4; 3 4 5]');
%!error<regression_ftest: Y must contain finite real numbers.> ...
%! regression_ftest ([1 2 3+i]', [2 3 4; 3 4 5]');
%!error<regression_ftest: X must contain finite real numbers.> ...
%! regression_ftest ([1 2 3]', [2 3 NaN; 3 4 5]');
%!error<regression_ftest: X must contain finite real numbers.> ...
%! regression_ftest ([1 2 3]', [2 3 Inf; 3 4 5]');
%!error<regression_ftest: X must contain finite real numbers.> ...
%! regression_ftest ([1 2 3]', [2 3 4; 3 4 3+i]');
%!error<regression_ftest: Y must be a vector of length 'rows \(X\)'.> ...
%! regression_ftest ([1 2 3]', [2 3; 3 4]');
%!error<regression_ftest: Y must be a vector of length 'rows \(X\)'.> ...
%! regression_ftest ([1 2; 3 4]', [2 3; 3 4]');
%!error<regression_ftest: optional arguments must be in pairs.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', 1, 'Alpha');
%!error<regression_ftest: invalid value for alpha.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', 1, 'Alpha', 0);
%!error<regression_ftest: invalid value for alpha.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', 1, 'Alpha', 1.2);
%!error<regression_ftest: invalid value for alpha.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', 1, 'Alpha', [0.02, 0.1]);
%!error<regression_ftest: invalid value for alpha.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', 1, 'Alpha', 'a');
%!error<regression_ftest: invalid value for Intercept.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', 1, 'Intercept', 'yes');
%!error<regression_ftest: invalid Name argument.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', 1, 'some', 0.05);
%!error<regression_ftest: invalid Name argument.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', 1, 3, 0.05);
%!error<regression_ftest: KEEP must hold distinct column indices of X.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', 3);
%!error<regression_ftest: KEEP must hold distinct column indices of X.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', 1.5);
%!error<regression_ftest: KEEP must hold distinct column indices of X.> ...
%! regression_ftest ((1:6)', [1:6; 2 1 4 3 6 5; 1 1 2 2 3 4]', [1, 1]);
%!error<regression_ftest: a logical KEEP must have one element per column of X.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', true);
%!error<regression_ftest: KEEP must leave out at least one column of X.> ...
%! regression_ftest ((1:5)', [1:5; 2 1 4 3 5]', [1, 2]);
%!error<regression_ftest: too few observations for the full model.> ...
%! regression_ftest ([1 2 3]', [1 2 3; 2 1 3]');
%!error<regression_ftest: the columns left out add nothing to the reduced model's design.> ...
%! regression_ftest ((1:6)', [1:6; 2 1 4 3 6 5; 1:6]', [1, 2]);
