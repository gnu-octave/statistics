## Copyright (C) 2003 Iain Murray
##
## This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} @var{s} = mvnrnd (@var{mu}, @var{Sigma})
## @deftypefnx{Function File} @var{s} = mvnrnd (@var{mu}, @var{Sigma}, @var{n})
## @deftypefnx{Function File} @var{s} = mvnrnd (@dots{}, @var{tol})
## Draw @var{n} random @var{d}-dimensional vectors from a multivariate Gaussian distribution with mean @var{mu}(@var{n}x@var{d}) and covariance matrix
## @var{Sigma}(@var{d}x@var{d}).
##
## @var{mu} must be @var{n}-by-@var{d} (or 1-by-@var{d} if @var{n} is given) or a scalar.
##
## If the argument @var{tol} is given the eigenvalues of @var{Sigma} are checked for positivity against -100*tol. The default value of tol is @code{eps*norm (Sigma, "fro")}.
##
## @end deftypefn

function s = mvnrnd (mu, Sigma, K, tol=eps*norm (Sigma, "fro"))

  % Iain Murray 2003 -- I got sick of this simple thing not being in Octave and locking up a stats-toolbox license in Matlab for no good reason.
  % May 2004 take a third arg, cases. Makes it more compatible with Matlab's.

  % Paul Kienzle <pkienzle@users.sf.net>
  % * Add GPL notice.
  % * Add docs for argument K

  % 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
  % * Uses Octave 3.6.2 broadcast.
  % * Stabilizes chol by perturbing Sigma with a epsilon multiple of the identity.
  %   The effect on the generated samples is to add additional independent noise of variance epsilon. Ref: GPML Rasmussen & Williams. 2006. pp 200-201
  % * Improved doc.
  % * Added tolerance to the positive definite check
  % * Used chol with option 'upper'.

  % 2014 Nir Krakauer <nkrakauer@ccny.cuny.edu>
  % * Add tests.
  % * Allow mu to be scalar, in which case it's assumed that all elements share this mean.
  
  
  %perform some input checking
  if ~issquare (Sigma)
    error ('Sigma must be a square covariance matrix.');
  end
    
  d = size(Sigma, 1);

  % If mu is column vector and Sigma not a scalar then assume user didn't read help but let them off and flip mu. Don't be more liberal than this or it will encourage errors (eg what should you do if mu is square?).
  if (size (mu, 2) == 1) && (d != 1)
    mu = mu';
  end

  if nargin >= 3
    n = K;
  else
    n = size(mu, 1); %1 if mu is scalar
  end

  if (~isscalar (mu)) && any(size (mu) != [1,d]) && any(size (mu) != [n,d])
    error ('mu must be nxd, 1xd, or scalar, where Sigma has dimensions dxd.');
  end
  
  warning ("off", "Octave:broadcast","local");

  try
    U = chol (Sigma + tol*eye (d),"upper");
  catch
    [E , Lambda] = eig (Sigma);

    if min (diag (Lambda)) < -100*tol
      error('Sigma must be positive semi-definite. Lowest eigenvalue %g', ...
            min (diag (Lambda)));
    else
      Lambda(Lambda<0) = 0;
    end
    warning ("mvnrnd:InvalidInput","Cholesky factorization failed. Using diagonalized matrix.")
    U = sqrt (Lambda) * E';
  end

  s = randn(n,d)*U + mu;

  warning ("on", "Octave:broadcast");
endfunction

% {{{ END OF CODE --- Guess I should provide an explanation:
%
% We can draw from axis aligned unit Gaussians with randn(d)
%   x ~ A*exp(-0.5*x'*x)
% We can then rotate this distribution using
%   y = U'*x
% Note that
%   x = inv(U')*y
% Our new variable y is distributed according to:
%   y ~ B*exp(-0.5*y'*inv(U'*U)*y)
% or
%   y ~ N(0,Sigma)
% where
%   Sigma = U'*U
% For a given Sigma we can use the chol function to find the corresponding U,
% draw x and find y. We can adjust for a non-zero mean by just adding it on.
%
% But the Cholsky decomposition function doesn't always work...
% Consider Sigma=[1 1;1 1]. Now inv(Sigma) doesn't actually exist, but Matlab's
% mvnrnd provides samples with this covariance st x(1)~N(0,1) x(2)=x(1). The
% fast way to deal with this would do something similar to chol but be clever
% when the rows aren't linearly independent. However, I can't be bothered, so
% another way of doing the decomposition is by diagonalising Sigma (which is
% slower but works).
% if
%   [E,Lambda]=eig(Sigma)
% then
%   Sigma = E*Lambda*E'
% so
%   U = sqrt(Lambda)*E'
% If any Lambdas are negative then Sigma just isn't even positive semi-definite
% so we can give up.
%
% Paul Kienzle adds:
%   Where it exists, chol(Sigma) is numerically well behaved.  chol(hilb(12)) for doubles and for 100 digit floating point differ in the last digit.
%   Where chol(Sigma) doesn't exist, X*sqrt(Lambda)*E' will be somewhat accurate.  For example, the elements of sqrt(Lambda)*E' for hilb(12), hilb(55) and hilb(120) are accurate to around 1e-8 or better.  This was tested using the TNT+JAMA for eig and chol templates, and qlib for 100 digit precision.
% }}}

%!shared m, n, C, rho
%! m = 10; n = 3; rho = 0.4; C = rho*ones(n, n) + (1 - rho)*eye(n);
%!assert(size(mvnrnd(0, C, m)), [m n])
%!assert(size(mvnrnd(zeros(1, n), C, m)), [m n])
%!assert(size(mvnrnd(zeros(n, 1), C, m)), [m n])
%!assert(size(mvnrnd(zeros(m, n), C, m)), [m n])
%!assert(size(mvnrnd(zeros(m, n), C)), [m n])
%!assert(size(mvnrnd(zeros(1, n), C)), [1 n])
%!assert(size(mvnrnd(zeros(n, 1), C)), [1 n])
%!error(mvnrnd(zeros(m+1, n), C, m))
%!error(mvnrnd(zeros(1, n+1), C, m))
%!error(mvnrnd(zeros(n+1, 1), C, m))
%!error(mvnrnd(zeros(m, n), eye(n+1), m))
%!error(mvnrnd(zeros(m, n), eye(n+1, n), m))

