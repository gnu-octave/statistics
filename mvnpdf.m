% y = mvnpdf(x,mu,sigma)
% Compute multivariate normal pdf for x given mean mu and covariance matrix sigma.
% x is dimension d x p, is 1 x p and sigma is p x p.

% This program is public domain
% Author: Paul Kienzle
function pdf = mvnpdf(x,mu,sigma)
  [d,p] = size(x);
  % mu can be a scalar, a 1xp vector or a nxp matrix
  if nargin == 1, mu = 0; end
  if all(size(mu) == [1,p]), mu = repmat(mu,[d,1]); end
  if nargin < 3
    pdf = (2*pi)^(-p/2) * exp(-0.5 * sumsq(x-mu,2));
  else
    r = chol(sigma);
    pdf = (2*pi)^(-p/2) * exp(-0.5 * sum(((x-mu)/r).^2,2)) / prod(diag(r));
  end

%!test
%! mu = [1,-1];
%! sigma = [.9 .4; .4 .3];
%! x = [ 0.5 -1.2; -0.5 -1.4; 0 -1.5];
%! p = [   0.41680003660313; 0.10278162359708; 0.27187267524566 ];
%! q = mvnpdf(x,mu,sigma);
%! assert(p,q,10*eps);
