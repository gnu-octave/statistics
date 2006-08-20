% y = mvnpdf(x,mu,Sigma)
% Compute multivariate normal pdf for x given mean mu and covariance matrix 
% sigma.  The dimension of x is d x p, mu is 1 x p and sigma is p x p.
%
% 1/y^2 = (2 pi)^p |\Sigma| exp { (x-\mu)' inv(\Sigma) (x-\mu) }
%
% Ref: NIST Engineering Statistics Handbook 6.5.4.2
% http://www.itl.nist.gov/div898/handbook/pmc/section5/pmc542.htm

% Algorithm
%
% Using Cholesky factorization on the positive definite covariance matrix:
%
%    r = chol(sigma);
%
% where r'*r = sigma. Being upper triangular, the determinant of r is 
% trivially the product of the diagonal, and the determinant of sigma 
% is the square of this:
%
%    det = prod(diag(r))^2;
%
% The formula asks for the square root of the determinant, so no need to 
% square it.
%
% The exponential argument A = x' inv(sigma) x
%	
%    A = x' inv(sigma) x = x' inv(r' * r) x = x' inv(r) inv(r') x
%
% Given that inv(r') == inv(r)', at least in theory if not numerically,
%
%    A  = (x' / r) * (x'/r)' = sumsq(x'/r)
%
% The interface takes the parameters to the mvn in columns rather than 
% rows, so we are actually dealing with the transpose:
%
%    A = sumsq(x/r)
%
% and the final result is:
%
%    r = chol(Sigma);
%    y = (2*pi)^(-p/2) * exp(-sumsq((x-mu)/r,2)/2) / prod(diag(r));


% This program is public domain
% Author: Paul Kienzle
function pdf = mvnpdf(x,mu,sigma)
  [d,p] = size(x);
  % mu can be a scalar, a 1xp vector or a nxp matrix
  if nargin == 1, mu = 0; end
  if all(size(mu) == [1,p]), mu = repmat(mu,[d,1]); end
  if nargin < 3
    pdf = (2*pi)^(-p/2) * exp(-sumsq(x-mu,2)/2);
  else
    r = chol(sigma);
    pdf = (2*pi)^(-p/2) * exp(-sumsq((x-mu)/r,2)/2) / prod(diag(r));
  end

%!test
%! mu = [1,-1];
%! sigma = [.9 .4; .4 .3];
%! x = [ 0.5 -1.2; -0.5 -1.4; 0 -1.5];
%! p = [   0.41680003660313; 0.10278162359708; 0.27187267524566 ];
%! q = mvnpdf(x,mu,sigma);
%! assert(p,q,10*eps);
