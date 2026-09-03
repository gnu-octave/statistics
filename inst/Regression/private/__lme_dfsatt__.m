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
## @deftypefn {Private Function} {@var{df} =} __lme_dfsatt__ (@var{X}, @var{y}, @var{Zx}, @var{qk}, @var{nlev}, @var{Psi}, @var{sigma2}, @var{method}, @var{L})
##
## Satterthwaite denominator degrees of freedom for the single-row contrasts in
## the rows of @var{L}.  Undocumented internal helper for
## @code{LinearMixedModel}.
##
## Uses the delta method: for a contrast @var{l}, the denominator variance is
## @code{l*C*l'} with @code{C} the fixed-effect covariance; its sampling
## variance is @code{g'*Veta*g} where @code{g} is the gradient of @code{l*C*l'}
## with respect to the covariance parameters @var{eta} (the lower-triangle
## entries of each @var{Psi} block followed by @var{sigma2}) and @code{Veta} is
## the asymptotic covariance of @var{eta} (twice the inverse of the numerical
## Hessian of the -2 log-(RE)likelihood).  The df is
## @code{2*(l*C*l')^2 / (g'*Veta*g)}.
##
## @end deftypefn

function df = __lme_dfsatt__ (X, y, Zx, qk, nlev, Psi, sigma2, method, L)

  is_reml = strcmp (upper (method), "REML");
  n = rows (X);
  p = columns (X);

  ## Pack the natural covariance parameters: vech(Psi_k) per term, then sigma2.
  eta = [];
  for k = 1:numel (qk)
    eta = [eta; vech_lower(Psi{k})];
  endfor
  eta = [eta; sigma2];
  ne = numel (eta);

  ## Every cross product below is free of eta, and the finite differences ask
  ## for O(ne^2) deviances, so forming them once is what keeps the whole
  ## Satterthwaite calculation out of the observation dimension.
  CP.ZtZ = Zx' * Zx;
  CP.ZtX = Zx' * X;
  CP.Zty = Zx' * y;
  CP.XtX = X' * X;
  CP.Xty = X' * y;
  CP.yty = y' * y;
  dev = @(e) reml_deviance (e, qk, nlev, CP, n, p, is_reml);
  cov = @(e) fixed_cov (e, qk, nlev, CP);

  ## Asymptotic covariance of eta: 2 * inv (Hessian of the -2 log-likelihood).
  hs = 1e-4 * max (abs (eta), 1);
  H = zeros (ne);
  for i = 1:ne
    for j = i:ne
      ei = zeros (ne, 1); ei(i) = hs(i);
      ej = zeros (ne, 1); ej(j) = hs(j);
      H(i,j) = (dev (eta+ei+ej) - dev (eta+ei-ej) ...
                - dev (eta-ei+ej) + dev (eta-ei-ej)) / (4 * hs(i) * hs(j));
      H(j,i) = H(i,j);
    endfor
  endfor
  Veta = 2 * inv (H);

  C0 = cov (eta);
  df = zeros (rows (L), 1);
  for r = 1:rows (L)
    l = L(r,:);
    lcl = l * C0 * l';
    g = zeros (ne, 1);
    for i = 1:ne
      ei = zeros (ne, 1); ei(i) = hs(i);
      g(i) = (l * cov (eta+ei) * l' - l * cov (eta-ei) * l') / (2 * hs(i));
    endfor
    df(r) = 2 * lcl^2 / (g' * Veta * g);
  endfor

endfunction

## Column-major lower-triangle entries of a symmetric matrix.
function v = vech_lower (A)
  q = rows (A);
  v = [];
  for j = 1:q
    for i = j:q
      v(end+1, 1) = A(i, j);
    endfor
  endfor
endfunction

## The block-diagonal random-effects covariance P and the residual variance,
## read out of eta.  V is sigma2*I + Zx*P*Zx', which is never formed: with
## G = Zx'*Zx the identity
##
##   inv (V) = I/s - Zx*P*inv (I + G*P/s)*Zx'/s^2
##
## and its determinant counterpart det (V) = s^n * det (I + P*G/s) take every
## quantity into the random effects dimension.  Note that this parameterises
## by the covariance rather than by a factor of it, and a finite difference
## drives a small variance component negative, so P is routinely indefinite
## here and no Cholesky of it exists.  The eigenvalues of P*G are real, P and
## G being symmetric with G positive semidefinite, and they carry both the
## log-determinant and the test that replaces the old chol of V: V is positive
## definite exactly when s > 0 and every 1 + lambda/s is.
function [P, s] = build_P (eta, qk, nlev)
  blocks = {};
  off = 0;
  for k = 1:numel (qk)
    q = qk(k);
    m = q*(q+1)/2;
    Pk = zeros (q, q);
    idx = 0;
    for j = 1:q
      for i = j:q
        idx += 1;
        Pk(i, j) = eta(off+idx);
        Pk(j, i) = eta(off+idx);
      endfor
    endfor
    for l = 1:nlev(k)
      blocks{end+1} = Pk;
    endfor
    off += m;
  endfor
  P = blkdiag (blocks{:});
  s = eta(end);
endfunction

## X'*inv (V)*X, X'*inv (V)*y, y'*inv (V)*y and logdet (V), all in the random
## effects dimension.  ok is false where the old code's chol of V would have
## set the flag.
function [XtViX, XtViy, ytViy, logdetV, ok] = vquad (eta, qk, nlev, CP, n)
  [P, s] = build_P (eta, qk, nlev);
  XtViX = XtViy = ytViy = logdetV = [];
  ok = false;
  if (s <= 0)
    return;
  endif
  PG = P * CP.ZtZ;
  lam = real (eig (PG));
  if (any (1 + lam / s <= 0))
    return;
  endif
  ok = true;
  logdetV = n * log (s) + sum (log (1 + lam / s));
  ## G*P is the transpose of P*G, both factors being symmetric.
  Wq = eye (columns (CP.ZtZ)) + PG' / s;
  Ax = P * (Wq \ CP.ZtX);
  Ay = P * (Wq \ CP.Zty);
  XtViX = CP.XtX / s - (CP.ZtX' * Ax) / s^2;
  XtViX = (XtViX + XtViX') / 2;
  XtViy = CP.Xty / s - (CP.ZtX' * Ay) / s^2;
  ytViy = CP.yty / s - (CP.Zty' * Ay) / s^2;
endfunction

function d = reml_deviance (eta, qk, nlev, CP, n, p, is_reml)
  [XtViX, XtViy, ytViy, logdetV, ok] = vquad (eta, qk, nlev, CP, n);
  if (! ok)
    d = Inf;
    return;
  endif
  beta = XtViX \ XtViy;
  rVir = ytViy - XtViy' * beta;
  if (is_reml)
    d = logdetV + rVir + 2*sum (log (diag (chol (XtViX)))) + (n-p)*log (2*pi);
  else
    d = logdetV + rVir + n*log (2*pi);
  endif
endfunction

function C = fixed_cov (eta, qk, nlev, CP)
  XtViX = vquad (eta, qk, nlev, CP, 1);
  C = inv (XtViX);
endfunction
