## [q,A] = anderson_darling_test(x,uniform|normal|exponential)
##
## Test the hypothesis that x is selected from the given distribution
## using the Anderson-Darling test.  Returns a small q value if the data
## are not distributed as expected.
##
## The Anderson-Darling A statistic is calculated as follows:
##
##    A = -n - \sum_{i=1}^n (2i-1)/n log(z_i (1-z_{n-i+1}))
##
## where z_i is the ordered position of the x's in the CDF of the 
## distribution.  Unlike the Kolmogorov-Smirnov statistic, the 
## Anderson-Darling statistic is sensitive to the tails of the 
## distribution.
##
## For 'normal' and 'exponential' distributions, estimate the 
## distribution parameters from the data, convert the values 
## to CDF values, and compare the result to tabluated critical 
## values.  This includes an correction for small n which 
## works well enough for n >= 8, but less so from smaller n.
##
## For 'uniform', assume the values are uniformly distributed
## in (0,1), compute A and return the corresponding p value from
## 1-anderson_darling_cdf(A,n).
## 
## If you are selecting from a known distribution, convert your 
## values into CDF values for the distribution and use 'uniform'.
## Do not use 'uniform' if the distribution parameters are estimated 
## from the data itself, as this sharply biases the A statistic
## toward smaller values.
##
## [1] Stephens, MA; (1986), "Tests based on EDF statistics", in
## D'Agostino, RB; Stephens, MA; (eds.) Goodness-of-fit Techinques.
## New York: Dekker.

## Author: Paul Kienzle
## This program is granted to the public domain.
function [q,A] = anderson_darling_test(x,dist)

  if size(x,1) == 1, x=x(:); end
  x = sort(x);
  n = size(x,1);
  use_cdf = 0;
  switch dist
    case 'normal',
      adj = 1 + (.75 + 2.25/n)/n;
      qvals = [ 1, 0.1, 0.05, 0.025, 0.01 ];
      Acrit = [-Inf, 0.631, 0.752, 0.873, 1.035]/adj;
      x = stdnormal_cdf(zscore(x));

    case 'uniform',
      ## Put invalid data at the limits of the distribution
      ## This will drive the statistic to infinity.
      x(x<0) = 0;
      x(x>1) = 1;
      #qvals = [ 1, 0.1, 0.05, 0.025, 0.01 ];
      #Acrit = [-Inf, 1.933, 2.492, 3.077, 3.878];
      use_cdf = 1;

    case 'XXXweibull',
      adj = 1 + 0.2/sqrt(n);
      Acrit = [-Inf, 0.637, 0.757, 0.877, 1.038]*adj;
      ## XXX FIXME XXX how to fit alpha and sigma?
      x = weibull_cdf(x,ones(n,1)*alpha,ones(n,1)*sigma);

    case 'exponential',
      qvals = [ 1, 0.1, 0.05, 0.025, 0.01 ];
      adj = 1 + 0.6/n;
      Acrit = [-Inf, 1.062, 1.321, 1.591, 1.959]/adj;

      lambda = 1./mean(x);  # exponential parameter estimation
      x = exponential_cdf(x,ones(n,1)*lambda);

    otherwise
      error("Anderson-Darling test for %s not implemented", dist);
  end

  if any(x<0 | x>1)
    error('Anderson-Darling test requires data in CDF form');
  endif

  i = [1:n]'*ones(1,size(x,2));
  A = -n - sum( (2*i-1) .* (log(x) + log(1-x(n:-1:1,:))) )/n;

  if use_cdf
    q = 1-anderson_darling_cdf(A, n);
  else
    q = qvals(lookup(Acrit,A));
  endif


%!demo
%! c = anderson_darling_test(10*rande(12,10000),'exponential');
%! tabulate(100*c,100*[unique(c),1]);
%! % The Fc column should report 100, 250, 500, 1000, 10000 more or less.

%!demo
%! c = anderson_darling_test(randn(12,10000),'normal');
%! tabulate(100*c,100*[unique(c),1]);
%! % The Fc column should report 100, 250, 500, 1000, 10000 more or less.

%!demo
%! c = anderson_darling_test(rand(12,10000),'uniform');
%! hist(100*c,1:2:99);
%! % The histogram should be flat more or less.
