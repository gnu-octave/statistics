## [p,A] = anderson_darling_test(x,uniform|normal|exponential)
##
## Test the hypothesis that x is selected from the given distribution
## using the Anderson-Darling test.
##
## The Anderson-Darling test statistic measures
##
##    A = -n - \sum_{i=1}^n (2i-1)/n log(z_i (1-z_{n-i+1}))
##
## where z_i is the ordered position of the x's in the CDF of the 
## distribution.  Unlike the Kolmogorov-Smirnov statistic, the 
## Anderson-Darling statistic is sensitive to the tails of the 
## distribution.
##
## uniform (0,1):   f(x, 'uniform')
## normal:          f((x-m)/sqrt(v),'normal')
## lognormal:       f(log (x - m)/sqrt(v),'normal')
## exponential:     f(lambda*x, 'exponential')
## weibull:         f((x/scale)^shape, 'exponential')
##
## The uses Marsaglia's approximation to the Anderson-Darling CDF
## to compute the p-values.
##
## These values differ significantly from those published elsewhere,
## but I do not yet understand why.  I have included published values 
## herein (the Acrit values), and the corrections to apply for small n.
## Uncomment the last line if you prefer to use them.
##
## Critical values and small sample adjustments for normal and
## exponential distributions comes from Charles Annis, P.E.
##     http://www.statisticalengineering.com/goodness.htm
## but identical information is available elsewhere on the web.

## Author: Paul Kienzle
## This program is granted to the public domain.
function [p,A] = anderson_darling_test(x,dist)

  if size(x,1) == 1, x=x(:); end
  x = sort(x);
  n = size(x,1);
  switch dist
    case 'normal',
      adj = 1 + .75/n + 2.25/n/n;
      pvals = [ 1, 0.1, 0.05, 0.025, 0.01 ];
      Acrit = [-Inf, 0.631, 0.752, 0.873, 1.035]/adj;
      x = (1 + erf(x/sqrt(2))) / 2;
    case 'uniform',
      ## Put invalid data at the limits of the distribution
      x(x<0) = 0;
      x(x>1) = 1;
      pvals = [ 1, 0.1, 0.05, 0.01 ];
      Acrit = [-Inf, 1.933, 2.492, 3.878];
    case 'exponential',
      ## Put invalid data at the limits of the distribution
      x(x<0) = 0;
      adj = 1 + 0.2/sqrt(n);
      pvals = [ 1, 0.1, 0.05, 0.025, 0.01 ];
      Acrit = [-Inf, 0.637, 0.757, 0.877, 1.038]/adj;
      x = 1 - exp(-x);
    otherwise
      error('use normal, exponential or uniform distribution');
  end
  if any(x<0 | x>1)
    error('Anderson-Darling test requires data in CDF form');
  endif

  i = [1:n]'*ones(1,size(x,2));
  A = -n - sum( (2*i-1) .* (log(x) + log(1-x(n:-1:1,:))) )/n;
  p = anderson_darling_cdf(A, n);
  ##p = pvals(lookup(Acrit,A));
