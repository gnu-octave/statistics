## Compute principal components of X
## [pc,z,w,Tsq] = princomp(X)
##   pc  the principal components
##   z   the transformed data
##   w   the eigenvalues of the covariance matrix
##   Tsq Hotelling's T^2 statistic for the transformed data
function [pc,z,w,Tsq] = princomp(X)
  C = cov(X);
  [U,D,pc] = svd(C);
  if nargout>1, z = center(X)*pc; end
  if nargout>2, w = diag(D); end
  if nargout>3, Tsq = sumsq(zscore(z),2); end
