## Copyright (C) 2003 Iain Murray
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} @var{s} = mvnrnd (@var{mu}, @var{Sigma})
## @deftypefnx{Function File} @var{s} = mvnrnd (@var{mu}, @var{Sigma}, @var{n})
## @deftypefnx{Function File} @var{s} = mvnrnd (@dots{}, @var{tol})
## Draw @var{n} random @var{d}-dimensional vectors from a multivariate Gaussian
## distribution with mean @var{mu}(@var{n}x@var{d}) and covariance matrix
## @var{Sigma}(@var{d}x@var{d}).
##
## If @var{n} is given then @var{mu} can be a 1-by-@var{d} vector. If it is not
## given @var{mu} must be @var{n}-by-@var{d}.
##
## If the argument @var{tol} is given the eigenvalues of @var{Sigma} are checked
## for positivity against -100*tol. the default value of tol is @code{eps*norm (Sigma, "fro")}.
##
## @end deftypefn

function s = mvnrnd (mu, Sigma, K, tol=eps*norm (Sigma, "fro"))

  % Iain Murray 2003 -- I got sick of this simple thing not being in Octave and
  %                     locking up a stats-toolbox license in Matlab for no good
  %                     reason.
  % May 2004 take a third arg, cases. Makes it more compatible with Matlab's.

  % Paul Kienzle <pkienzle@users.sf.net>
  % * Add GPL notice.
  % * Add docs for argument K

  % 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
  % * Uses Octave 3.6.2 boradcast.
  % * Stabilizes chol by perturbing Sigma with a epsilon multiple of the identity.
  %   The effect on the generated samples is to add additional independent noise
  %   of variance epsilon. Ref: GPML Rasmussen & Williams. 2006. pp 200-201
  % * Improved doc.
  % * Added tolerance to the positive definite check
  % * Used chol with option 'upper'.

  warning ("off", "Octave:broadcast");

  % If mu is column vector and Sigma not a scalar then assume user didn't read
  % help but let them off and flip mu. Don't be more liberal than this or it will
  % encourage errors (eg what should you do if mu is square?).
  if (size (mu, 2) == 1) && (size (Sigma) != [1,1])
    mu = mu';
  end

  [n,d] = size (mu);

  if nargin >= 3
     n = K;
  end

  if size (Sigma) != [d,d]
    error ('Sigma must have dimensions dxd where mu is nxd.');
  end

  try
    U = chol (Sigma + tol*eye (d),"upper");
  catch
    [E , Lambda] = eig (Sigma);

    if min (diag (Lambda)) < -100*tol
      error('Sigma must be positive semi-definite. Lowest eigenvaluen %g', ...
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
%   Where it exists, chol(Sigma) is numerically well behaved.  chol(hilb(12))
%   for doubles and for 100 digit floating point differ in the last digit.
%   Where chol(Sigma) doesn't exist, X*sqrt(Lambda)*E' will be somewhat
%   accurate.  For example, the elements of sqrt(Lambda)*E' for hilb(12),
%   hilb(55) and hilb(120) are accurate to around 1e-8 or better.  This was
%   tested using the TNT+JAMA for eig and chol templates, and qlib for
%   100 digit precision.
% }}}
