## p = anderson_darling_cdf(A,n)
##
## Return the CDF for the given Anderson-Darling coefficient A
## computed from n values sampled from a distribution. For a 
## vector of random variables x of length n, compute the CDF
## of the values from the distribution from which they are drawn.
## You can uses these values to compute A as follows:
##
##    A = -n - sum( (2*i-1) .* (log(x) + log(1-x(n:-1:1,:))) )/n;
## 
## From the value A, anderson_darling_cdf returns the probability
## that A could be returned from a set of samples.
##
## The algorithm given in [1] claims to be an approximation for the
## Anderson-Darling CDF accurate to 6 decimal points.
##
## Demonstrate using:
##    n=300; reps=10000;
##    z=randn(n,reps);
##    x=sort((1+erf(z/sqrt(2)))/2);
##    i = [1:n]'*ones(1,size(x,2));
##    A = -n - sum( (2*i-1) .* (log(x) + log(1-x(n:-1:1,:))) )/n;
##    p=anderson_darling_cdf(A,n);
##    hist(100*p,[1:100]-0.5);
## You will see that the histogram is basically flat, which is to
## say that the probabilities returned by the Anderson-Darling CDF 
## are distributed uniformly.
##
## You can easily determine the extreme values of p:
##    [junk,idx]=sort(p);
## The histograms of various p aren't  very informative:
##    histfit(z(:,idx(1)),linspace(-3,3,15));
##    histfit(z(:,idx(end/2)),linspace(-3,3,15));
##    histfit(z(:,idx(end)),linspace(-3,3,15));
## More telling is the qqplot:
##    qqplot(z(:,idx(1))); hold on; plot([-3,3],[-3,3],';;'); hold off;
##    qqplot(z(:,idx(end/2))); hold on; plot([-3,3],[-3,3],';;'); hold off;
##    qqplot(z(:,idx(end))); hold on; plot([-3,3],[-3,3],';;'); hold off;
##
## Try a similarly analysis for z uniform:
##    z=rand(n,reps); x=sort(z);
## and for z exponential:
##    z=rande(n,reps); x=sort(1-exp(-z));
##
## [1] Marsaglia, G; Marsaglia JCW; (2004) "Evaluating the Anderson Darling
## distribution", Journal of Statistical Software, 9(2).

## Author: Paul Kienzle
## This code is granted to the public domain.
  
function y = anderson_darling_cdf(z,n)
  y = ADinf(z);
  y += ADerrfix(y,n);
end

function y = ADinf(z)
  y = zeros(size(z));

  idx = (z < 2);
  if any(idx(:))
    p = [.00168691, -.0116720, .0347962, -.0649821, .247105, 2.00012];
    z1 = z(idx);
    y(idx) = exp(-1.2337141./z1)./sqrt(z1).*polyval(p,z1);
  end

  idx = (z >= 2);
  if any(idx(:))
    p = [-.0003146, +.008056, -.082433, +.43424, -2.30695, 1.0776]; 
    y(idx) = exp(-exp(polyval(p,z(idx))));
  end
end
    
function y = ADerrfix(x,n)
  if isscalar(n), n = n*ones(size(x));
  elseif isscalar(x), x = x*ones(size(n));
  end
  y = zeros(size(x));
  c = .01265 + .1757./n;

  idx = (x >= 0.8);
  if any(idx(:))
    p = [255.7844, -1116.360, 1950.646, -1705.091, 745.2337, -130.2137];
    g3 = polyval(p,x(idx));
    y(idx) = g3./n(idx);
  end

  idx = (x < 0.8 & x > c);
  if any(idx(:))
    p = [1.91864, -8.259, 14.458, -14.6538, 6.54034, -.00022633];
    n1 = 1./n(idx);
    c1 = c(idx);
    g2 = polyval(p,(x(idx)-c1)./(.8-c1));
    y(idx) = (.04213 + .01365*n1).*n1 .* g2;
  end

  idx = (x <= c);
  if any(idx(:))
    x1 = x(idx)./c(idx);
    n1 = 1./n(idx);
    g1 = sqrt(x1).*(1-x1).*(49*x1-102);
    y(idx) = ((.0037*n1+.00078).*n1+.00006).*n1 .* g1;
  end
end
