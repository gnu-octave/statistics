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
## @deftypefn  {Function File} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} vartest (@var{x}, @var{y})
## @deftypefnx {Function File} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} vartest (@var{x}, @var{y}, @var{Name}, @var{Value})
## Perform a F-test for equal variances.
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
## Description: Test for mean of a normal sample with known variance

function [h, p, ci, stats] = vartest(x, v, varargin)
  
  % Set default arguments
  alpha = 0.05;
  tail  = 'both';
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
      case 'dim'
        i = i + 1;
        dim = varargin{i};
      otherwise
        error('Invalid Name argument.',[]);
    end
    i = i + 1;
  end
  
  if ~isa(tail, 'char')
    error('tail argument to vartest must be a string\n',[]);
  end
    
  s_var = var(x, 0, dim);
  
  df = size(x, dim) - 1;
  stats.chisqstat = df*s_var/v;

  % Based on the "tail" argument determine the P-value, the critical values,
  % and the confidence interval.
  switch lower(tail)
    case 'both'
      p = 2*min(chi2cdf(stats.chisqstat,df),1-chi2cdf(stats.chisqstat,df));
      ci = [df*s_var ./ (chi2inv(1-alpha/2,df)); df*s_var ./ (chi2inv(alpha/2,df))];
    case 'left'
      p = chi2cdf(stats.chisqstat,df);
      chi2crit = chi2inv(alpha,df);
      ci = [zeros(size(stats.chisqstat)); df*s_var ./ (chi2inv(alpha,df))];
    case 'right'
      p = 1 - chi2cdf(stats.chisqstat,df);
      chi2crit = chi2inv(1-alpha,df);
      ci = [df*s_var ./ (chi2inv(1-alpha,df)); inf*ones(size(stats.chisqstat))];
    otherwise
      error('Invalid fourth (tail) argument to vartest\n',[]);
  end

  % Reshape the ci array to match MATLAB shaping
  if and(isscalar(stats.chisqstat), dim==2)
    ci = ci(:)';
  elseif size(stats.chisqstat,2)<size(stats.chisqstat,1)
    ci = reshape(ci(:),length(stats.chisqstat),2);
  end

  stats.df = df*ones(size(stats.chisqstat));
  
  % Determine the test outcome
  % MATLAB returns this a double instead of a logical array
  h = double(p < alpha);  
end
