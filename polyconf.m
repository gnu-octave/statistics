## [y,dy] = polyconf(p,x,s)
##
##   Produce prediction intervals for the fitted y. The vector p 
##   and structure s are returned from polyfit or wpolyfit. The 
##   x values are where you want to compute the prediction interval.
##
## polyconf(...,['ci'|'pi'])
##
##   Produce a confidence interval (range of likely values for the
##   mean at x) or a prediction interval (range of likely values 
##   seen when measuring at x).  The prediction interval tells
##   you the width of the distribution at x.  This should be the same
##   regardless of the number of measurements you have for the value
##   at x.  The confidence interval tells you how well you know the
##   mean at x.  It should get smaller as you increase the number of
##   measurements.  Error bars in the physical sciences usually show 
##   a 1-alpha confidence value of erfc(1/sqrt(2)), representing
##   one standandard deviation of uncertainty in the mean.
##
## polyconf(...,1-alpha)
##
##   The width of the prediction interval. The default is .05 for
##   a 95% prediction interval, or erfc(1/sqrt(2)) for a one
##   standard deviation confidence interval.
##
## Example:
##  [p,s] = polyfit(x,y,1);
##  xf = linspace(x(1),x(end),150);
##  [yf,dyf] = polyconf(p,xf,s,'ci');
##  plot(xf,yf,'g-;fit;',xf,yf+dyf,'g.;;',xf,yf-dyf,'g.;;',x,y,'xr;data;');
##  plot(x,y-polyval(p,x),';residuals;',xf,dyf,'g-;;',xf,-dyf,'g-;;');
function [y,dy] = polyconf(p,x,varargin)
  alpha = S = [];
  default_alpha = 0.05;
  pred = 1;
  for i=1:length(varargin)
    v = varargin{i};
    if isstruct(v), S = v;
    elseif isstr(v),
      switch v, 
	case 'ci', pred = 0; default_alpha=erfc(1/sqrt(2));
	case 'pi', pred = 1; default_alpha=0.05;
	otherwise, error("polyconf(...,['pi'|'ci'])");
      end
    elseif isscalar(v), alpha = v;
    else S = [];
    end
  end
  if isempty(alpha), alpha = default_alpha; end
  if (nargout>1 && (isempty(S)||nargin<3)) || nargin < 2
    usage("[y,dy] = polyconf(p,x,s,alpha,['conf'|'pred'])");
  end

  ## Confidence interval for linear system are given by:
  ##    x' p +/- sqrt( Finv(1-a,1,df) var(x' p) )
  ## where for confidence intervals,
  ##    var(x' p) = sigma^2 (x' inv(A'A) x)
  ## and for prediction intervals,
  ##    var(x' p) = sigma^2 (1 + x' inv(A'A) x)
  ##
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
  ## Note: sqrt(F(1-a;1,df)) = T(1-a/2;df)
  ##
  ## For a polynomial fit, x is the set of powers ( x^n ; ... ; 1 ).
  y=dy=x;
  s = t_inv(1-alpha/2,S.df)*S.normr/sqrt(S.df); 
  n=length(p)-1;
  k=length(x(:));
  A = (x(:) * ones (1, n+1)) .^ (ones (k, 1) * (n:-1:0));
  y(:) = A*p(:);
  dy(:) = s*sqrt(pred+sumsq(A/S.R,2));
