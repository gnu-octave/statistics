## Copyright (C) 2014 Tony Richardson
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
## @deftypefn  {Function File} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} ttest2 (@var{x}, @var{y})
## @deftypefnx {Function File} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} ttest2 (@var{x}, @var{y}, @var{Name}, @var{Value})
## Test for mean of a normal sample with known variance.
##
## Perform a T-test of the null hypothesis @code{mean (@var{x}) ==
## @var{m}} for a sample @var{x} from a normal distribution with unknown
## mean and unknown std deviation.  Under the null, the test statistic
## @var{t} has a Student's t distribution.
##
## If the second argument @var{y} is a vector, a paired-t test of the
## hypothesis @code{mean (@var{x}) = mean (@var{y})} is performed.
##
## The argument @qcode{"alpha"} can be used to specify the significance level
## of the test (the default value is 0.05).  The string
## argument @qcode{"tail"}, can be used to select the desired alternative
## hypotheses.  If @qcode{"alt"} is @qcode{"both"} (default) the null is 
## tested against the two-sided alternative @code{mean (@var{x}) != @var{m}}.
## If @qcode{"alt"} is @qcode{"right"} the one-sided 
## alternative @code{mean (@var{x}) > @var{m}} is considered.
## Similarly for @qcode{"left"}, the one-sided alternative @code{mean
## (@var{x}) < @var{m}} is considered.  When @qcode{"vartype"} is @qcode{"equal"}
## the variances are assumed to be equal (this is the default).  When
## @qcode{"vartype"} is @qcode{"unequal"} the variances are not assumed equal.
## When argument @var{x} is a matrix the @qcode{"dim"} argument can be 
## used to selection the dimension over which to perform the test.
## (The default is the first non-singleton dimension.)
##
## If @var{h} is 0 the null hypothesis is accepted, if it is 1 the null
## hypothesis is rejected. The p-value of the test is returned in @var{pval}.
## A 100(1-alpha)% confidence interval is returned in @var{ci}. @var{stats}
## is a structure containing the value of the test statistic (@var{tstat}),
## the degrees of freedom (@var{df}) and the sample standard deviation
## (@var{sd}).
##
## @end deftypefn

## Author: Tony Richardson <richardson.tony@gmail.com>

function [h, p, ci, stats] = ttest2(x, y, varargin)
  
  alpha = 0.05;
  tail  = 'both';
  vartype = 'equal';

  % Find the first non-singleton dimension of x
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end

  i = 1;
  while ( i <= length(varargin) )
    switch lower(varargin{i})
      case 'alpha'
        i = i + 1;
        alpha = varargin{i};
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
        error('Invalid Name argument.',[]);
    end
    i = i + 1;
  end
  
  if ~isa(tail, 'char')
    error('Tail argument to ttest2 must be a string\n',[]);
  end
  
  m = size(x, dim);
  n = size(y, dim);
  x_bar = mean(x,dim)-mean(y,dim);
  s1_var = var(x, 0, dim);
  s2_var = var(y, 0, dim);

  switch lower(vartype)
    case 'equal'
      stats.tstat = 0;
      stats.df = (m + n - 2)*ones(size(x_bar));
      sp_var = ((m-1)*s1_var + (n-1)*s2_var)./stats.df;
      stats.sd = sqrt(sp_var);
      x_bar_std = sqrt(sp_var*(1/m+1/n));
    case 'unequal'
      stats.tstat = 0;
      se1 = sqrt(s1_var/m);
      se2 = sqrt(s2_var/n);
      sp_var = s1_var/m + s2_var/n;
      stats.df = ((se1.^2+se2.^2).^2 ./ (se1.^4/(m-1) + se2.^4/(n-1)));
      stats.sd = [sqrt(s1_var); sqrt(s2_var)];
      x_bar_std = sqrt(sp_var);
    otherwise
      error('Invalid fifth (vartype) argument to ttest2\n',[]);
  end

  stats.tstat = x_bar./x_bar_std;

  % Based on the "tail" argument determine the P-value, the critical values,
  % and the confidence interval.
  switch lower(tail)
    case 'both'
      p = 2*(1 - tcdf(abs(stats.tstat),stats.df));
      tcrit = -tinv(alpha/2,stats.df);
      %ci = [x_bar-tcrit*stats.sd; x_bar+tcrit*stats.sd];
      ci = [x_bar-tcrit.*x_bar_std; x_bar+tcrit.*x_bar_std];
    case 'left'
      p = tcdf(stats.tstat,stats.df);
      tcrit = -tinv(alpha,stats.df);
      ci = [-inf*ones(size(x_bar)); x_bar+tcrit.*x_bar_std];
    case 'right'
      p = 1 - tcdf(stats.tstat,stats.df);
      tcrit = -tinv(alpha,stats.df);
      ci = [x_bar-tcrit.*x_bar_std; inf*ones(size(x_bar))];
    otherwise
      error('Invalid fourth (tail) argument to ttest2\n',[]);
  end

  % Reshape the ci array to match MATLAB shaping
  if and(isscalar(x_bar), dim==2)
    ci = ci(:)';
    stats.sd = stats.sd(:)';
  elseif size(x_bar,2)<size(x_bar,1)
    ci = reshape(ci(:),length(x_bar),2);
    stats.sd = reshape(stats.sd(:),length(x_bar),2);
  end

  % Determine the test outcome
  % MATLAB returns this a double instead of a logical array
  h = double(p < alpha);  
end
