## Copyright (C) 2014 Tony Richardson <richardson.tony@gmail.com>
## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} ttest2 (@var{x}, @var{y})
## @deftypefnx {statistics} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} ttest2 (@var{x}, @var{y}, @var{Name}, @var{Value})
##
## Perform a t-test to compare the means of two groups of data under the null
## hypothesis that the groups are drawn from distributions with the same mean.
##
## @var{x} and @var{y} can be vectors or matrices. For matrices, @qcode{ttest2}
## performs separate t-tests along each column, and returns a vector of results.
## @var{x} and @var{y} must have the same number of columns. The Type I error
## rate of the resulting vector of @var{pval} can be controlled by entering
## @var{pval} as input to the function @qcode{multcompare}.
##
## @qcode{ttest2} treats NaNs as missing values, and ignores them.
##
## For a nested t-test, use @qcode{anova2}.
##
## The argument @qcode{'alpha'} can be used to specify the significance level
## of the test (the default value is 0.05).  The string argument @qcode{'tail'},
## can be used to select the desired alternative hypotheses.  If @qcode{'tail'}
## is @qcode{'both'} (default) the null is tested against the two-sided
## alternative @code{mean (@var{x}) != @var{m}}.  If @qcode{'tail'} is
## @qcode{'right'} the one-sided alternative @code{mean (@var{x}) > @var{m}} is
## considered.  Similarly for @qcode{'left'}, the one-sided alternative
## @code{mean (@var{x}) < @var{m}} is considered.
##
## When @qcode{'vartype'} is @qcode{'equal'} the variances are assumed to be
## equal (this is the default).  When @qcode{'vartype'} is @qcode{'unequal'} the
## variances are not assumed equal.
##
## When argument @var{x} and @var{y} are matrices the @qcode{'dim'} argument can
## be used to select the dimension over which to perform the test.
## (The default is the first non-singleton dimension.)
##
## If @var{h} is 0 the null hypothesis is accepted, if it is 1 the null
## hypothesis is rejected. The p-value of the test is returned in @var{pval}.
## A 100(1-alpha)% confidence interval is returned in @var{ci}. @var{stats}
## is a structure containing the value of the test statistic (@var{tstat}),
## the degrees of freedom (@var{df}) and the sample standard deviation
## (@var{sd}).
##
## @seealso{hotelling_t2test, anova1, hotelling_t2test2, ttest}
## @end deftypefn

function [h, p, ci, stats] = ttest2 (x, y, varargin)

  ## Set defaults
  alpha = 0.05;
  tail = 'both';
  vartype = 'equal';
  ## Find the first non-singleton dimension of x
  dim = min (find (size (x) != 1));
  if (isempty (dim))
    dim = 1;
  endif

  ## Evaluate optional input arguments
  i = 1;
  while ( i <= length (varargin) )
    switch lower (varargin{i})
      case 'alpha'
        i = i + 1;
        alpha = varargin{i};
        if (! (isscalar (alpha) && isnumeric (alpha) && isreal (alpha)
               && alpha > 0 && alpha < 1))
          error ("ttest2: ALPHA must be a scalar between 0 and 1.");
        endif
      case 'tail'
        i = i + 1;
        tail = varargin{i};
      case 'vartype'
        i = i + 1;
        vartype = varargin{i};
      case 'dim'
        i = i + 1;
        dim = varargin{i};
      otherwise
        error ("ttest2: Invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile

  ## Error checking
  if (! isa (tail, 'char'))
    error ("ttest2: tail argument must be a string.");
  endif
  if (size (x, abs (dim - 3)) != size (y, abs (dim - 3)))
    error ("ttest2: The data in a 2-sample t-test must be commensurate")
  endif

  ## Calculate mean, variance and size of each sample
  m = sum (! isnan (x), dim);
  n = sum (! isnan (y), dim);
  x_bar = mean (x, dim, 'omitnan') - mean (y, dim, 'omitnan');
  s1_var = var (x, 0, dim, 'omitnan');
  s2_var = var (y, 0, dim, 'omitnan');

  ## Perform test-specific calculations
  switch lower (vartype)
    case 'equal'
      stats.tstat = [];
      stats.df = (m + n - 2);
      sp_var = ((m - 1) .* s1_var + (n - 1) .* s2_var) ./ stats.df;
      stats.sd = sqrt (sp_var);
      x_bar_std = sqrt (sp_var .* (1 ./ m + 1 ./ n));
      n_sd = 1;
    case 'unequal'
      stats.tstat = [];
      se1 = sqrt (s1_var ./ m);
      se2 = sqrt (s2_var ./ n);
      sp_var = s1_var ./ m + s2_var ./ n;
      stats.df = ((se1 .^ 2 + se2 .^ 2) .^ 2 ./ ...
                  (se1 .^ 4 ./ (m - 1) + se2 .^ 4 ./ (n - 1)));
      stats.sd = [sqrt(s1_var); sqrt(s2_var)];
      x_bar_std = sqrt (sp_var);
      n_sd = 2;
    otherwise
      error ("ttest2: Invalid value for vartype argument.");
  endswitch
  stats.tstat = x_bar ./ x_bar_std;

  ## Based on the "tail" argument determine the P-value, the critical values,
  ## and the confidence interval.
  switch lower (tail)
    case 'both'
      p = 2 * (1 - tcdf (abs (stats.tstat), stats.df));
      tcrit = - tinv (alpha / 2, stats.df);
      ci = [x_bar-tcrit.*x_bar_std; x_bar+tcrit.*x_bar_std];
    case 'left'
      p = tcdf (stats.tstat, stats.df);
      tcrit = - tinv (alpha, stats.df);
      ci = [-inf*ones(size(x_bar)); x_bar+tcrit.*x_bar_std];
    case 'right'
      p = 1 - tcdf (stats.tstat, stats.df);
      tcrit = - tinv (alpha, stats.df);
      ci = [x_bar-tcrit.*x_bar_std; inf*ones(size(x_bar))];
    otherwise
      error ("ttest2: Invalid value for tail argument.");
  endswitch

  ## Reshape the ci array to match MATLAB shaping
  if (isscalar (x_bar) && dim == 2)
    ci = ci(:)';
    stats.sd = stats.sd(:)';
  elseif (size (x_bar, 2) < size (x_bar, 1))
    ci = reshape (ci(:), length (x_bar), 2);
    stats.sd = reshape (stats.sd(:), length (x_bar), n_sd);
  endif

  ## Determine the test outcome
  ## MATLAB returns this a double instead of a logical array
  h = double (p < alpha);

endfunction

%!test
%! a = 1:5;
%! b = 6:10;
%! b(5) = NaN;
%! [h,p,ci,stats] = ttest2 (a,b);
%! assert_equal (h, 1);
%! assert_equal (p, 0.002535996080258229, 1e-14);
%! assert_equal (ci, [-6.822014919225481, -2.17798508077452], 1e-14);
%! assert_equal (stats.tstat, -4.582575694955839, 1e-14);
%! assert_equal (stats.df, 7);
%! assert_equal (stats.sd, 1.4638501094228, 1e-13);
%!error ttest2 ([8:0.1:12], [8:0.1:12], 'tail', 'invalid');
%!error ttest2 ([8:0.1:12], [8:0.1:12], 'tail', 25);

## Reference values from MATLAB R2024a (probe run 2026-08-02).  The single
## existing block covered only the two-sample default, leaving both 'Vartype'
## values, both one-sided tails, 'alpha' and 'dim' unchecked.  'Vartype'
## matters most: 'unequal' is Welch's test, with its own denominator and a
## non-integer degrees of freedom.
%!shared x, y, xm
%! x = [10.2 9.7 11.1 10.5 9.9 10.8 10.1 9.6 10.4 10.7];
%! y = [11.4 10.9 12.6 11.1 13.2 10.4 12.8 11.7];
%! xm = [10.2 11.4; 9.7 10.9; 11.1 12.6; 10.5 11.1; 9.9 13.2];
%!test
%! [h, p, ci, stats] = ttest2 (x, y, 'Vartype', 'equal');
%! assert_equal (h, 1);
%! assert_equal (p, 8.895270853487265e-04, 1e-15);
%! assert_equal (ci, [-2.224122032184765, -0.700877967815235], 1e-13);
%! assert_equal (stats.tstat, -4.070735048482627, 1e-13);
%! assert_equal (stats.df, 16);
%! assert_equal (stats.sd, 0.757411298436985, 1e-14);
%!test
%! ## Welch: a pooled denominator would give the 'equal' answer instead
%! [h, p, ci, stats] = ttest2 (x, y, 'Vartype', 'unequal');
%! assert_equal (h, 1);
%! assert_equal (p, 0.003801117252046, 1e-14);
%! assert_equal (ci, [-2.327640087870322, -0.597359912129679], 1e-13);
%! assert_equal (stats.tstat, -3.784559445642580, 1e-13);
%! assert_equal (stats.df, 9.661941319801253, 1e-13);
%! assert_equal (stats.sd(:)', [0.489897948556636, 1.001338390070295], 1e-14);
%!test
%! ## the default is the equal-variance test
%! [~, pd] = ttest2 (x, y);
%! [~, pe] = ttest2 (x, y, 'Vartype', 'equal');
%! [~, pu] = ttest2 (x, y, 'Vartype', 'unequal');
%! assert_equal (pd, pe, 0);
%! assert_equal (isequal (pd, pu), false);
%!test
%! [h, p, ci] = ttest2 (x, y, 'Tail', 'left');
%! assert_equal (h, 1);
%! assert_equal (p, 4.447635426743633e-04, 1e-15);
%! assert_equal (ci, [-Inf, -0.835253361212791], 1e-13);
%!test
%! [h, p, ci] = ttest2 (x, y, 'Tail', 'right');
%! assert_equal (h, 0);
%! assert_equal (p, 0.999555236457326, 1e-14);
%! assert_equal (ci, [-2.089746638787210, Inf], 1e-13);
%!test
%! ## a non-default alpha must widen the interval and leave the p-value alone
%! [~, p1, ci1] = ttest2 (x, y);
%! [~, p2, ci2] = ttest2 (x, y, 'Alpha', 0.01);
%! assert_equal (ci2, [-2.511854249766015, -0.413145750233985], 1e-13);
%! assert_equal (diff (ci2) > diff (ci1), true);
%! assert_equal (p1, p2, 0);
%!test
%! [h, p] = ttest2 (xm, xm + 1, 'Dim', 1);
%! assert_equal (h, [1, 0]);
%! assert_equal (p, [0.020600636177229, 0.154833732538475], 1e-13);

## ALPHA was not validated at all, so a negative one silently produced a
## [NaN NaN] confidence interval and "do not reject", and a vector one
## quietly vectorised into an undocumented signature.  MATLAB refuses each of
## these with "ALPHA must be a scalar between 0 and 1."
%!error <ttest2: ALPHA must be a scalar between 0 and 1.> ...
%! ttest2 ([8:0.1:12], [9:0.1:13], 'Alpha', -0.05);
%!error <ttest2: ALPHA must be a scalar between 0 and 1.> ...
%! ttest2 ([8:0.1:12], [9:0.1:13], 'Alpha', 0);
%!error <ttest2: ALPHA must be a scalar between 0 and 1.> ...
%! ttest2 ([8:0.1:12], [9:0.1:13], 'Alpha', 1);
%!error <ttest2: ALPHA must be a scalar between 0 and 1.> ...
%! ttest2 ([8:0.1:12], [9:0.1:13], 'Alpha', 1.5);
%!error <ttest2: ALPHA must be a scalar between 0 and 1.> ...
%! ttest2 ([8:0.1:12], [9:0.1:13], 'Alpha', [0.01, 0.05]);
%!error <ttest2: ALPHA must be a scalar between 0 and 1.> ...
%! ttest2 ([8:0.1:12], [9:0.1:13], 'Alpha', 'a');
%!error <ttest2: ALPHA must be a scalar between 0 and 1.> ...
%! ttest2 ([8:0.1:12], [9:0.1:13], 'Alpha', NaN);
%!error <ttest2: ALPHA must be a scalar between 0 and 1.> ...
%! ttest2 ([8:0.1:12], [9:0.1:13], 'Alpha', 2 + 1i);
