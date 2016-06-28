## Copyright (C) 2016 Andreas Stahel
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
## @deftypefn  {Function File} {[@var{h}, @var{pval}, @var{ci}] =} binotest (@var{pos},@var{N},@var{p0})
## @deftypefnx  {Function File} {[@var{h}, @var{pval}, @var{ci}] =} binotest (@var{pos},@var{N},@var{p0},@var{Name},@var{Value})
## Test for probability @var{p} of a binomial sample
##
## Perform a test of the null hypothesis @var{p} ==  @var{p0} for a sample
## of size @var{N} with @var{pos} positive results
##
##
## Name-Value pair arguments can be used to set various options.
## @qcode{"alpha"} can be used to specify the significance level
## of the test (the default value is 0.05). The option @qcode{"tail"},
## can be used to select the desired alternative hypotheses.  If the
## value is @qcode{"both"} (default) the null is tested against the two-sided 
## alternative @code{@var{p} != @var{p0}}. The value of @var{pval} is
## determined by adding the probabilities of all event less or equally
## likely than the observed number @var{pos} of positive events.
## If the value of @qcode{"tail"} is @qcode{"right"}
## the one-sided alternative @code{@var{p} > @var{p0}} is considered.
## Similarly for @qcode{"left"}, the one-sided alternative
## @code{@var{p} < @var{p0}} is considered.
##
## If @var{h} is 0 the null hypothesis is accepted, if it is 1 the null
## hypothesis is rejected. The p-value of the test is returned in @var{pval}.
## A 100(1-alpha)% confidence interval is returned in @var{ci}.
##
## @end deftypefn


## Author: Andreas Stahel <Andreas.Stahel@bfh.ch.
## based on the code ttest.m by Tony Richardson <richardson.tony@gmail.com>

function [h, p, ci] = binotest(pos,n,p0,varargin)
  
  % Set default arguments
  alpha = 0.05;
  tail  = 'both';
  
  i = 1;
  while ( i <= length(varargin) )
    switch lower(varargin{i})
      case 'alpha'
        i = i + 1;
        alpha = varargin{i};
      case 'tail'
        i = i + 1;
        tail = varargin{i};
      otherwise
        error('Invalid Name argument.',[]);
    end
    i = i + 1;
  end
  
  if ~isa(tail,'char')
    error('tail argument to vartest must be a string\n',[]);
  end

  if (n<=0)
    error('binotest: required n>0\n',[]);
  end
  if (p0<0)|(p0>1)
    error('binotest: required 0<= p0 <= 1\n',[]);
  end
  if (pos<0)|(pos>n)
    error('binotest: required 0<= pos <= n\n',[]);
  end

  % Based on the "tail" argument determine the P-value, the critical values,
  % and the confidence interval.
  switch lower(tail)
    case 'both'
      A_low = binoinv(alpha/2,n,p0)/n;
      A_high = binoinv(1-alpha/2,n,p0)/n;
      p_pos = binopdf(pos,n,p0);
      p_all = binopdf([0:n],n,p0);
      ind = find(p_all <=p_pos);
%      p = min(1,sum(p_all(ind)));
      p = sum(p_all(ind));
      if pos==0 p_low = 0;
	 else   p_low = fzero(@(pl)1-binocdf(pos-1,n,pl)-alpha/2,[0 1]);
      endif
      if pos==n p_high = 1;
         else   p_high = fzero(@(ph)  binocdf(pos,n,ph)  -alpha/2,[0,1]);
      endif
      ci = [p_low,p_high];
    case 'left'
      p = 1-binocdf(pos-1,n,p0);
      if pos==n p_high = 1;
         else   p_high = fzero(@(ph)  binocdf(pos,n,ph)  -alpha,[0,1]);
      endif
      ci = [0, p_high];
    case 'right'
      p = binocdf(pos,n,p0);
      if pos==0 p_low = 0;
	 else   p_low = fzero(@(pl)1-binocdf(pos-1,n,pl)-alpha,[0 1]);
      endif
      ci = [p_low 1];
    otherwise
      error('Invalid fifth (tail) argument to binotest\n',[]);
  end

  % Determine the test outcome
  % MATLAB returns this a double instead of a logical array
  h = double(p < alpha);
end

%!demo
%! % flip a coin 1000 times, showing 475 heads
%! % Hypothesis: coin is fair, i.e. p=1/2
%! [h,p_val,ci] = binotest(475,1000,0.5)
%! % Result: h = 0 : null hypothesis not rejected, coin could be fair
%! %         P value 0.12, i.e. hypothesis not rejected for alpha up to 12%
%! %         0.444 <= p <= 0.506 with 95% confidence

%!demo
%! % flip a coin 100 times, showing 65 heads
%! % Hypothesis: coin shows less than 50% heads, i.e. p<=1/2
%! [h,p_val,ci] = binotest(65,100,0.5,'tail','left','alpha',0.01)
%! % Result: h = 1 : null hypothesis is rejected, i.e. coin shows more heads than tails
%! %         P value 0.0018, i.e. hypothesis not rejected for alpha up to 0.18%
%! %         0 <= p <= 0.76 with 99% confidence

%!test #example from https://en.wikipedia.org/wiki/Binomial_test
%! [h,p_val,ci] = binotest (51,235,1/6);
%! assert (p_val, 0.0437, 0.00005)
%! [h,p_val,ci] = binotest (51,235,1/6,'tail','left');
%! assert (p_val, 0.027, 0.0005)
