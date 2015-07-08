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
## @deftypefn  {Function File} {[@var{pval}, @var{h}, @var{stats}] =} signtest (@var{x})
## @deftypefnx {Function File} {[@var{pval}, @var{h}, @var{stats}] =} signtest (@var{x}, @var{m})
## @deftypefnx {Function File} {[@var{pval}, @var{h}, @var{stats}] =} signtest (@var{x}, @var{y})
## @deftypefnx {Function File} {[@var{pval}, @var{h}, @var{stats}] =} signtest (@var{x}, @var{y}, @var{Name}, @var{Value})
## Test for median.
##
## Perform a signtest of the null hypothesis that @var{x} is from a distribution
## that has a zero median.
##
## If the second argument @var{m} is a scalar, the null hypothesis is that
## X has median m.
##
## If the second argument @var{y} is a vector, the null hypothesis is that
## the distribution of @code{@var{x} - @var{y}} has zero median.
##
## The argument @qcode{"alpha"} can be used to specify the significance level
## of the test (the default value is 0.05).  The string
## argument @qcode{"tail"}, can be used to select the desired alternative
## hypotheses.  If @qcode{"alt"} is @qcode{"both"} (default) the null is 
## tested against the two-sided alternative @code{median (@var{x}) != @var{m}}.
## If @qcode{"alt"} is @qcode{"right"} the one-sided 
## alternative @code{median (@var{x}) > @var{m}} is considered.
## Similarly for @qcode{"left"}, the one-sided alternative @code{median
## (@var{x}) < @var{m}} is considered.  When @qcode{"method"} is @qcode{"exact"}
## the p-value is computed using an exact method (this is the default).  When
## @qcode{"method"} is @qcode{"approximate"} a normal approximation is used for the
## test statistic.
##
## The p-value of the test is returned in @var{pval}. If @var{h} is 0 the 
## null hypothesis is accepted, if it is 1 the null hypothesis is rejected. 
## @var{stats} is a structure containing the value of the test statistic
## (@var{sign}) and the value of the z statistic (@var{zval}) (only computed
## when the 'method' is 'approximate'.
##
## @end deftypefn

## Author: Tony Richardson <richardson.tony@gmail.com>

function [p, h, stats] = signtest(x, my, varargin)
  
  my_default = 0;
  alpha = 0.05;
  tail  = 'both';
  method = 'exact';

  % Find the first non-singleton dimension of x
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end

  if (nargin == 1)
    my = my_default;
  end 
  
  i = 1;
  while ( i <= length(varargin) )
    switch lower(varargin{i})
      case 'alpha'
        i = i + 1;
        alpha = varargin{i};
      case 'tail'
        i = i + 1;
        tail = varargin{i};
      case 'method'
        i = i + 1;
        method = varargin{i};
      case 'dim'
        i = i + 1;
        dim = varargin{i};
      otherwise
        error('Invalid Name argument.',[]);
    end
    i = i + 1;
  end
  
  if ~isa(tail, 'char')
    error('tail argument to signtest must be a string\n',[]);
  end

  if ~isa(method, 'char')
    error('method argument to signtest must be a string\n',[]);
  end

  % Set default values if arguments are present but empty
  if isempty(my)
    my = my_default;
  end

  % This adjustment allows everything else to remain the
  % same for both the one-sample t test and paired tests.
  % If second argument is a vector
  if ~isscalar(my)
    x = x - my;
    my = my_default;
  end  

  n = size(x, dim);

  switch lower(method)
    case 'exact'
      stats.zval = nan;
      switch lower(tail)
        case 'both'
          w = min(sum(x<my),sum(x>my));
          pl = binocdf(w, n, 0.5);
          p = 2*min(pl,1-pl);
        case 'left'
          w = sum(x<my);
          p = binocdf(w, n, 0.5);
        case 'right'
          w = sum(x>my);
          p = 1 - binocdf(w, n, 0.5);
        otherwise
          error('Invalid tail argument to signtest\n',[]);
      end
    case 'approximate'
      switch lower(tail)
        case 'both'
          npos = sum(x>my);
          nneg = sum(x<my);
          w = min(npos,nneg);
          stats.zval = (w - 0.5*n - 0.5*sign(npos-nneg))/sqrt(0.25*n);
          pl = normcdf(stats.zval);
          p = 2*min(pl,1-pl);
        case 'left'
          npos = sum(x>my);
          nneg = sum(x<my);
          w = sum(x<my);
          stats.zval = (w - 0.5*n - 0.5*sign(npos-nneg))/sqrt(0.25*n);
          p = normcdf(stats.zval);
        case 'right'
          npos = sum(x>my);
          nneg = sum(x<my);
          w = sum(x>my);
          stats.zval = (w - 0.5*n - 0.5*sign(npos-nneg))/sqrt(0.25*n);
          p = 1-normcdf(stats.zval);
        otherwise
          error('Invalid tail argument to signtest\n',[]);
      end
    otherwise
      error('Invalid method argument to signtest\n',[]);
  end

  stats.sign = w;
  
  h = double(p < alpha);  
  
end
