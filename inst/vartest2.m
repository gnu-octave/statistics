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
## @deftypefn  {Function File} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} vartest2 (@var{x}, @var{y})
## @deftypefnx {Function File} {[@var{h}, @var{pval}, @var{ci}, @var{stats}] =} vartest2 (@var{x}, @var{y}, @var{Name}, @var{Value})
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

function [h, p, ci, stats] = vartest2(x, y, varargin)
  
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
    error('tail argument to vartest2 must be a string\n',[]);
  end
    
  s1_var = var(x, 0, dim);
  s2_var = var(y, 0, dim);
  
  stats.fstat = s1_var ./ s2_var;
  df1= size(x, dim) - 1;
  df2 = size(y, dim) - 1;

  % Based on the "tail" argument determine the P-value, the critical values,
  % and the confidence interval.
  switch lower(tail)
    case 'both'
      p = 2*min(fcdf(stats.fstat,df1,df2),1 - fcdf(stats.fstat,df1,df2));
      fcrit = finv(1-alpha/2,df1,df2);
      ci = [s1_var ./ (fcrit*s2_var); fcrit*s1_var ./ s2_var];
    case 'left'
      p = fcdf(stats.fstat,df1,df2);
      fcrit = finv(alpha,df1,df2);
      ci = [zeros(size(stats.fstat)); s1_var ./ (fcrit*s2_var)];
    case 'right'
      p = 1 - fcdf(stats.fstat,df1,df2);
      fcrit = finv(1-alpha,df1,df2);
      ci = [s1_var ./ (fcrit*s2_var); inf*ones(size(stats.fstat))];
    otherwise
      error('Invalid fourth (tail) argument to vartest2\n',[]);
  end

  % Reshape the ci array to match MATLAB shaping
  if and(isscalar(stats.fstat), dim==2)
    ci = ci(:)';
  elseif size(stats.fstat,2)<size(stats.fstat,1)
    ci = reshape(ci(:),length(stats.fstat),2);
  end

  stats.df1 = df1*ones(size(stats.fstat));
  stats.df2 = df2*ones(size(stats.fstat));
  
  % Determine the test outcome
  % MATLAB returns this a double instead of a logical array
  h = double(p < alpha);  
end
