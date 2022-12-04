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
## @deftypefn  {Function File} [@var{h}, @var{pval}, @var{ci}, @var{stats}] = ttest2 (@var{x}, @var{y})
## @deftypefnx {Function File} [@var{h}, @var{pval}, @var{ci}, @var{stats}] = ttest2 (@var{x}, @var{y}, @var{Name}, @var{Value})
##
## Perform a t-test to compare the means of two groups of data under the null
## hypothesis that the groups are drawn from distributions with the same mean.
##
## @var{x} and @var{y} can be vectors or matrices. For matrices, ttest2 performs
## separate t-tests along each column, and returns a vector of results. @var{x}
## and @var{y} must have the same number of columns. The Type I error rate of
## the resulting vector of p-values can be controlled by using the p-values as
## input to the function @qcode{multcompare}.
##
## For a nested t-test, use @qcode{anova2}.
##
## The argument @qcode{"alpha"} can be used to specify the significance level
## of the test (the default value is 0.05).  The string
## argument @qcode{"tail"}, can be used to select the desired alternative
## hypotheses.  If @qcode{"tail"} is @qcode{"both"} (default) the null is
## tested against the two-sided alternative @code{mean (@var{x}) != @var{m}}.
## If @qcode{"tail"} is @qcode{"right"} the one-sided
## alternative @code{mean (@var{x}) > @var{m}} is considered.
## Similarly for @qcode{"left"}, the one-sided alternative @code{mean
## (@var{x}) < @var{m}} is considered.  When @qcode{"vartype"} is @qcode{"equal"}
## the variances are assumed to be equal (this is the default).  When
## @qcode{"vartype"} is @qcode{"unequal"} the variances are not assumed equal.
## When argument @var{x} and @var{y} are matrices the @qcode{"dim"} argument can
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
## @seealso{anova2, multcompare}
## @end deftypefn

function [h, p, ci, stats] = ttest2(x, y, varargin)

  ## Set defaults 
  alpha = 0.05;
  tail = "both";
  vartype = "equal";
  ## Find the first non-singleton dimension of x
  dim = min (find (size (x) != 1));
  if (isempty (dim))
    dim = 1;
  endif

  ## Evaluate optional input arguments
  i = 1;
  while ( i <= length(varargin) )
    switch lower(varargin{i})
      case "alpha"
        i = i + 1;
        alpha = varargin{i};
      case "tail"
        i = i + 1;
        tail = varargin{i};
      case "vartype"
        i = i + 1;
        vartype = varargin{i};
      case "dim"
        i = i + 1;
        dim = varargin{i};
      otherwise
        error ("ttest2: Invalid Name argument.");
    endswitch
    i = i + 1;
  endwhile

  ## Error checking
  if (! isa (tail, "char"))
    error ("ttest2: tail argument must be a string.");
  endif
  if (size (x, abs (dim - 3)) != size (y, abs (dim - 3)))
    error ("ttest2: The data in a 2-sample t-test must be commensurate")
  endif

  ## Calculate mean, variance and size of each sample
  m = sum (!isnan (x), dim);
  n = sum (!isnan (y), dim);
  x_bar = mean (x, dim, "omitnan") - mean (y, dim, "omitnan");
  s1_var = nanvar (x, 0, dim);
  s2_var = nanvar (y, 0, dim);

  ## Perform test-specific calculations
  switch lower (vartype)
    case "equal"
      stats.tstat = [];
      stats.df = (m + n - 2);
      sp_var = ((m - 1) .* s1_var + (n - 1) .* s2_var) ./ stats.df;
      stats.sd = sqrt (sp_var);
      x_bar_std = sqrt (sp_var .* (1 ./ m + 1 ./ n));
      n_sd = 1;
    case "unequal"
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
  end
  stats.tstat = x_bar ./ x_bar_std;

  ## Based on the "tail" argument determine the P-value, the critical values,
  ## and the confidence interval.
  switch lower(tail)
    case "both"
      p = 2 * (1 - tcdf (abs (stats.tstat), stats.df));
      tcrit = - tinv (alpha / 2, stats.df);
      ci = [x_bar-tcrit.*x_bar_std; x_bar+tcrit.*x_bar_std];
    case "left"
      p = tcdf (stats.tstat, stats.df);
      tcrit = - tinv (alpha, stats.df);
      ci = [-inf*ones(size(x_bar)); x_bar+tcrit.*x_bar_std];
    case "right"
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
%! [h,p,ci,stats] = ttest2(a,b);
%! assert (h, 1);
%! assert (p, 0.002535996080258229, 1e-14);
%! assert (ci, [-6.822014919225481, -2.17798508077452], 1e-14);
%! assert (stats.tstat, -4.582575694955839, 1e-14);
%! assert (stats.df, 7);
%! assert (stats.sd, 1.4638501094228, 1e-13);
%!error ttest2 ([8:0.1:12], [8:0.1:12], "tail", "invalid");
%!error ttest2 ([8:0.1:12], [8:0.1:12], "tail", 25);
