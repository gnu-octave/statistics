## [y,dy] = polyconf(p,x,s[,1-alpha])
## Produce confidence values for the fitted y. The vector p and structure s are
## returned from polyfit or wpolyfit.  The x values are where you want to compute
## the confidence interval.  1-alpha is the width of the confidence interval.  The
## default is .05 for a 95% confidence interval.  Standard error bars have a
## 1-alpha confidence value of erfc(1/sqrt(2)).
##
## Example:
##  [p,s] = polyfit(x,y,1);
##  xf = linspace(x(1),x(end),150);
##  [yf,dyf] = polyconf(p,xf,s,erfc(1/sqrt(2)));
##  plot(xf,yf,'g-;fit;',xf,yf+dyf,'g.;;',xf,yf-dyf,'g.;;',x,y,'xr;data;');
##  plot(x,y-polyval(p,x),';residuals;',xf,dyf,'g-;;',xf,-dyf,'g-;;');
function [y,dy] = polyconf(p,x,S,alpha)
  if nargin < 4, alpha = 0.05; end
  if !struct_contains(S,'sig'), S.sig = S.normr / S.df; end
  n=length(p)-1;
  k=length(x(:));
  ## Confidence interval for linear system are given by:
  ##    x' p +/- sqrt( Finv(a,1,df) var(x' p) )
  ## where
  ##    var(x' p) = sigma^2 x' inv(A'A) x
  ## Rather than A'A we have R from the QR decomposition of A, but
  ## R'R equals A'A.  Note that R is not upper triangular since we
  ## have already multiplied it by the permutation matrix, but it
  ## is invertible.  Rather than forming the product R'R which is
  ## ill-conditioned, we can rewrite x' inv(A'A) x as the equivalent
  ##    x' inv(R) inv(R') x = t t', for t = x' inv(R)
  ## Since x is a vector, t t' is the inner product sumsq(t).
  ## Note that LAPACK allows us to do this simultaneously for many
  ## different x using sqrt(sumsq(X/R,2)), with each x on a different row.
  ##
  ## For a polynomial fit, x is the set of powers ( x^n ; ... ; 1 ).
  s = sqrt(f_inv(1-alpha,1,S.df))*S.sig;
  A = (x(:) * ones (1, n+1)) .^ (ones (k, 1) * (n:-1:0));
  y=dy=x;
  y(:) = A*p(:);
  dy(:) = s*sqrt(sumsq(A/S.R,2));
