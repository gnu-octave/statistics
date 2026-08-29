## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
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
## @deftypefn {private} {@var{M} =} __lmefit__ (@var{X}, @var{y}, @var{Z}, @var{G}, @var{method})
##
## Estimation core for linear mixed-effects models.  Undocumented internal
## function for @code{fitlmematrix}/@code{fitlme}/@code{LinearMixedModel}.
##
## Fits @code{y = X*beta + Z*b + e} with @code{b ~ N(0, Psi)} and
## @code{e ~ N(0, sigma2*I)} by maximum likelihood (@var{method} = @qcode{"ML"})
## or restricted maximum likelihood (@qcode{"REML"}).  @var{Z} and @var{G} are
## cell arrays, one entry per grouping term: @var{Z}@{k@} is the @var{n}-by-q_k
## random-effects design and @var{G}@{k@} the @var{n}-by-1 grouping index of
## term k.  Each grouping term has a full (unstructured) q_k-by-q_k covariance,
## shared across the levels of that term.
##
## Returns a struct with the fitted quantities: @code{beta}, @code{Psi} (cell,
## one covariance per term), @code{sigma2}, @code{b} (BLUPs), @code{loglik},
## @code{covbeta}, the profiled covariance parameters @code{theta}, and the
## bookkeeping needed to rebuild the random design (@code{Zx}, @code{qk},
## @code{nlev}, @code{levels}).
##
## @end deftypefn

function M = __lmefit__ (X, y, Z, G, method)

  if (nargin < 5 || isempty (method))
    method = "ML";
  endif
  method = upper (method);
  is_reml = strcmp (method, "REML");

  n = rows (X);
  p = columns (X);
  nt = numel (Z);

  ## Build the expanded random design Zx (n-by-N, N = sum_k q_k*nlev_k) and the
  ## per-term bookkeeping.  Column block for (term k, level l) holds the rows of
  ## Z{k} that belong to level l and zeros elsewhere.
  Zx = [];
  qk = zeros (1, nt);
  nlev = zeros (1, nt);
  levels = cell (1, nt);
  gidx = cell (1, nt);
  for k = 1:nt
    qk(k) = columns (Z{k});
    [lev, ~, gi] = unique (G{k}(:));
    nlev(k) = numel (lev);
    levels{k} = lev;
    gidx{k} = gi;
    blk = zeros (n, qk(k) * nlev(k));
    for l = 1:nlev(k)
      rows_l = (gi == l);
      cols_l = (l-1)*qk(k) + (1:qk(k));
      blk(rows_l, cols_l) = Z{k}(rows_l, :);
    endfor
    Zx = [Zx, blk];
  endfor

  ## theta layout: lower-triangular Cholesky entries of the *relative*
  ## covariance D_k = L_k*L_k' (Psi_k = sigma2*D_k) per term, concatenated.
  theta0 = init_theta (qk);

  obj = @(th) profiled_deviance (th, X, y, Zx, qk, nlev, n, p, is_reml);

  ## The profiled deviance is differentiable in theta and the gradient is
  ## closed form, so give it to the optimiser rather than let it difference.
  ## A finite difference takes a step near 1e-8 on an objective whose
  ## curvature in the variance components is around 1e-13, which makes the
  ## search path, and the estimate it stops at, sensitive to the last bits of
  ## the linear algebra and so to the BLAS in use.
  opts = optimset ("TolX", 1e-10, "TolFun", 1e-10, "MaxFunEvals", 2000, ...
                   "MaxIter", 1000, "Display", "off", "GradObj", "on");
  [theta, dev] = fminunc (obj, theta0, opts);

  ## Recover everything at the optimum.
  Dfull = build_Dfull (theta, qk, nlev);
  ## Zx * Dfull * Zx' is symmetric in exact arithmetic but not bitwise, since
  ## the two products are separate BLAS calls.  'chol' reads one triangle, so
  ## which one it reads decides the factor and a threaded or blocked BLAS can
  ## disagree with a reference one.  Symmetrise before every factorisation.
  Mrel = eye (n) + Zx * Dfull * Zx';
  Mrel = (Mrel + Mrel') / 2;
  Rc = chol (Mrel);
  ## Whitened once, as in profiled_deviance above.  Mir is still formed in
  ## full because the BLUPs need inv (Mrel)*r and not just its norm.
  W = Rc' \ X;
  z = Rc' \ y;
  XtMiX = W' * W;
  beta = W \ z;
  r = y - X * beta;
  u = Rc' \ r;
  Mir = Rc \ u;
  rMir = u' * u;
  if (is_reml)
    sigma2 = rMir / (n - p);
  else
    sigma2 = rMir / n;
  endif

  ## BLUPs and covariance of the fixed effects.  The covariance is genuinely
  ## an inverse, but XtMiX is symmetric positive definite and W'*W already
  ## has its factor, so take it from there rather than from a general inverse.
  b = Dfull * (Zx' * Mir);
  Rx = chol (XtMiX);
  covbeta = sigma2 * (Rx \ (Rx' \ eye (p)));
  loglik = -0.5 * dev;

  ## Per-term covariance matrices Psi_k = sigma2 * D_k.
  Psi = cell (1, nt);
  off = 0;
  for k = 1:nt
    m = qk(k)*(qk(k)+1)/2;
    Lk = tril_from_theta (theta(off+(1:m)), qk(k));
    Psi{k} = sigma2 * (Lk * Lk');
    off += m;
  endfor

  M = struct ();
  M.beta = beta;
  M.Psi = Psi;
  M.sigma2 = sigma2;
  M.b = b;
  M.loglik = loglik;
  M.covbeta = covbeta;
  M.theta = theta;
  M.method = method;
  M.Zx = Zx;
  M.qk = qk;
  M.nlev = nlev;
  M.levels = levels;
  M.gidx = gidx;
  M.fitted = X * beta + Zx * b;      ## conditional fitted values
  M.fitted_marg = X * beta;          ## marginal (fixed-effects-only) fit
  M.resid = y - M.fitted;            ## raw conditional residuals
  M.n = n;
  M.p = p;
  M.dfe = n - p;

endfunction

## Profiled -2*log-likelihood at the covariance parameters theta, and its
## gradient.  With M = I + Zx*D*Zx' and nu = n or n-p, the deviance is
## nu*log (rMir/nu) + logdet (M) [+ logdet (X'*inv (M)*X)] + const, and each
## term differentiates through dM/dtheta = Zx*(dD/dtheta)*Zx'.  Because beta
## is the GLS minimiser, its own dependence on theta drops out.  Every
## derivative then contracts to the random effects dimension: with
## V = Rc'\Zx, the three quantities needed are V'*V, V'*u and V'*W, and since
## D = L*L' gives dD = E_ij*L' + L*E_ij', each gradient entry is 2*(S*L)(i,j)
## for the matching accumulated S.
function [dev, grad] = profiled_deviance (theta, X, y, Zx, qk, nlev, n, p, ...
                                          is_reml)
  Dfull = build_Dfull (theta, qk, nlev);
  Mrel = eye (n) + Zx * Dfull * Zx';
  Mrel = (Mrel + Mrel') / 2;
  [Rc, flag] = chol (Mrel);
  if (flag != 0)
    dev = Inf;
    grad = zeros (size (theta));
    return;
  endif
  ## Whiten once rather than solving twice.  Mrel is Rc'*Rc, so
  ## X'*inv (Mrel)*X is (Rc'\X)'*(Rc'\X): three triangular solves here where
  ## six were needed, and XtMiX comes out symmetric to the last bit by
  ## construction rather than needing to be symmetrised.  beta is then the
  ## least squares solution of the whitened system, which avoids forming the
  ## normal equations and squaring the condition number with them.
  W = Rc' \ X;
  z = Rc' \ y;
  XtMiX = W' * W;
  beta = W \ z;
  r = y - X * beta;
  u = Rc' \ r;
  rMir = u' * u;
  logdetM = 2 * sum (log (diag (Rc)));
  if (is_reml)
    nu = n - p;
    Rx = chol (XtMiX);
    ldXMiX = 2 * sum (log (diag (Rx)));
    dev = nu*log (rMir/nu) + logdetM + rMir/(rMir/nu) + ldXMiX + nu*log (2*pi);
  else
    nu = n;
    dev = nu*log (rMir/nu) + logdetM + rMir/(rMir/nu) + nu*log (2*pi);
  endif
  if (nargout < 2)
    return;
  endif

  ## One further triangular solve carries the whole gradient: V'*V is
  ## Zx'*inv (M)*Zx, V'*u is Zx'*inv (M)*r, and V'*W is Zx'*inv (M)*X, so
  ## nothing works in the observation dimension beyond this point.
  V = Rc' \ Zx;
  c = V' * u;
  B = V' * W;
  if (is_reml)
    T = Rx' \ B';
  endif
  grad = zeros (size (theta));
  off = 0;
  col = 0;
  for k = 1:numel (qk)
    q = qk(k);
    m = q*(q+1)/2;
    Lk = tril_from_theta (theta(off+(1:m)), q);
    Sa = zeros (q);
    Sc = zeros (q);
    Sp = zeros (q);
    for l = 1:nlev(k)
      S = col + (1:q);
      col += q;
      Vs = V(:,S);
      Sa += Vs' * Vs;
      cs = c(S);
      Sc += cs * cs';
      if (is_reml)
        Ts = T(:,S);
        Sp += Ts' * Ts;
      endif
    endfor
    ## dD = E_ij*L' + L*E_ij', so tr (dD*S) and v'*dD*v are both 2*(S*L)(i,j)
    ## for symmetric S.  logdet (M) rises with it, rMir and logdet (X'MiX)
    ## fall.
    Ga = 2 * (Sa * Lk);
    Gc = 2 * (Sc * Lk);
    if (is_reml)
      Gp = 2 * (Sp * Lk);
    else
      Gp = zeros (q);
    endif
    idx = 0;
    for j = 1:q
      for i = j:q
        idx += 1;
        grad(off+idx) = -nu * Gc(i,j) / rMir + Ga(i,j) - Gp(i,j);
      endfor
    endfor
    off += m;
  endfor
endfunction

## Block-diagonal relative covariance: D_k repeated over the nlev_k levels.
function Dfull = build_Dfull (theta, qk, nlev)
  blocks = {};
  off = 0;
  for k = 1:numel (qk)
    m = qk(k)*(qk(k)+1)/2;
    Lk = tril_from_theta (theta(off+(1:m)), qk(k));
    Dk = Lk * Lk';
    for l = 1:nlev(k)
      blocks{end+1} = Dk;
    endfor
    off += m;
  endfor
  Dfull = blkdiag (blocks{:});
endfunction

## Lower-triangular q-by-q factor from its q*(q+1)/2 entries (column-major
## lower triangle).
function L = tril_from_theta (th, q)
  L = zeros (q, q);
  idx = 0;
  for j = 1:q
    for i = j:q
      idx += 1;
      L(i, j) = th(idx);
    endfor
  endfor
endfunction

## Starting covariance parameters: identity relative factors (D_k = I).
function theta0 = init_theta (qk)
  theta0 = [];
  for k = 1:numel (qk)
    L = eye (qk(k));
    th = [];
    for j = 1:qk(k)
      for i = j:qk(k)
        th(end+1) = L(i, j);
      endfor
    endfor
    theta0 = [theta0, th];
  endfor
  theta0 = theta0(:);
endfunction
