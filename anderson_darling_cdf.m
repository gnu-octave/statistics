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
## Test using:
##    n=300;
##    z=randn(n,1000);
##    x=(1+erf(z/sqrt(2))/2;
##    x=sort(x);
##    i = [1:n]'*ones(1,size(x,2));
##    A = -n - sum( (2*i-1) .* (log(x) + log(1-x(n:-1:1,:))) )/n;
##    p=anderson_darling_cdf(A,n);
##    hist(10*p,[1:10]-0.5);
## You will see that the histogram is basically flat.
## Similarly for uniform:
##    x=rand(n,1000);
## And for exponential:
##    x=1-exp(-rande(n,1000));
##
## Examine some of the more extreme values of p:
##    [junk,idx]=sort(p);
##    hist(z(:,idx(1)),[-4:4]);
##    hist(z(:,idx(2)),[-4:4]);
##    hist(z(:,idx(end)),[-4:4]);
##    hist(z(:,idx(end-1),[-4:4]);
##
## Repeat the experiment for x=z=rand(n,1000) and x=rande(n,1000)
## z=1-exp(-x).
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
  endif

  idx = (x < 0.8 & x > c);
  if any(idx(:))
    p = [1.91864, -8.259, 14.458, -14.6538, 6.54034, -.00022633];
    n1 = 1./n(idx);
    c1 = c(idx);
    g2 = polyval(p,(x(idx)-c1)./(.8-c1));
    y(idx) = (.04213 + .01365*n1).*n1 .* g2;
  endif

  idx = (x <= c);
  if any(idx(:))
    x1 = x(idx)./c(idx);
    n1 = 1./n(idx);
    g1 = sqrt(x1).*(1-x1).*(49*x1-102);
    y(idx) = ((.0037*n1+.00078).*n1+.00006).*n1 .* g1;
  endif
