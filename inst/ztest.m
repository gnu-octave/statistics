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
## @deftypefn  {Function File} {[@var{h}, @var{pval}, @var{ci}, @var{z}, @var{zcrit}] =} ztest (@var{x}, @var{m}, @var{s})
## @deftypefnx {Function File} {[@var{h}, @var{pval}, @var{ci}, @var{z}, @var{zcrit}] =} ztest (@var{x}, @var{m}, @var{s}, @var{Name}, @var{Value})
## Test for mean of a normal sample with known variance.
##
## Perform a Z-test of the null hypothesis @code{mean (@var{x}) == @var{m}}
## for a sample @var{x} from a normal distribution with unknown
## mean and known std deviation @var{s}.  Under the null, the test statistic
## @var{z} follows a standard normal distribution.
##
## Name-Value pair arguments can be used to set various options.
## @qcode{"alpha"} can be used to specify the significance level
## of the test (the default value is 0.05).  @qcode{"tail"}, can be used
## to select the desired alternative hypotheses.  If the value is
## @qcode{"both"} (default) the null is tested against the two-sided 
## alternative @code{mean (@var{x}) != @var{m}}.
## If it is @qcode{"right"} the one-sided alternative @code{mean (@var{x})
## > @var{m}} is considered.  Similarly for @qcode{"left"}, the one-sided 
## alternative @code{mean (@var{x}) < @var{m}} is considered.  
## When argument @var{x} is a matrix, @qcode{"dim"} can be used to selection
## the dimension over which to perform the test.  (The default is the 
## first non-singleton dimension.)
##
## If @var{h} is 0 the null hypothesis is accepted, if it is 1 the null
## hypothesis is rejected. The p-value of the test is returned in @var{pval}.
## A 100(1-alpha)% confidence interval is returned in @var{ci}.  The test statistic
## value is returned in @var{z} and the z critical value in @var{zcrit}.
##
## @end deftypefn

## Author: Tony Richardson <richardson.tony@gmail.com>

function [h, p, ci, zval, zcrit] = ztest(x, m, sigma, varargin)
  
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
    error('tail argument to ztest must be a string\n',[]);
  end
  
  % Calculate the test statistic value (zval)
  n = size(x, dim);
  x_bar = mean(x, dim);
  x_bar_std = sigma/sqrt(n);
  zval = (x_bar - m)./x_bar_std;
  
  % Based on the "tail" argument determine the P-value, the critical values,
  % and the confidence interval.
  switch lower(tail)
    case 'both'
      p = 2*(1 - normcdf(abs(zval)));
      zcrit = -norminv(alpha/2);
      ci = [x_bar-zcrit*x_bar_std; x_bar+zcrit*x_bar_std];
    case 'left'
      p = normcdf(zval);
      zcrit = -norminv(alpha);
      ci = [-inf*ones(size(x_bar)); x_bar+zcrit*x_bar_std];
    case 'right'
      p = 1 - normcdf(zval);
      zcrit = -norminv(alpha);
      ci = [x_bar-zcrit*x_bar_std; inf*ones(size(x_bar))];
    otherwise
      error('Invalid fifth (tail) argument to ztest\n',[]);
  end
  
  % Reshape the ci array to match MATLAB shaping
  if and(isscalar(x_bar), dim==2)
    ci = ci(:)';
  elseif size(x_bar,2)<size(x_bar,1)
    ci = reshape(ci(:),length(x_bar),2);
  end

  % Determine the test outcome
  % MATLAB returns this a double instead of a logical array
  h = double(p < alpha);
  
end
