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
## @deftypefn {private} {@var{df} =} __lme_dfsatt__ (@var{X}, @var{y}, @var{Zx}, @var{qk}, @var{nlev}, @var{Psi}, @var{sigma2}, @var{method}, @var{L})
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

  dev = @(e) reml_deviance (e, qk, nlev, X, y, Zx, n, p, is_reml);
  cov = @(e) fixed_cov (e, qk, nlev, X, Zx, n);

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

## Rebuild the marginal covariance V = Zx*Dabs*Zx' + sigma2*I from eta.
function V = build_V (eta, qk, nlev, Zx, n)
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
  sigma2 = eta(end);
  V = Zx * blkdiag (blocks{:}) * Zx' + sigma2 * eye (n);
endfunction

function d = reml_deviance (eta, qk, nlev, X, y, Zx, n, p, is_reml)
  V = build_V (eta, qk, nlev, Zx, n);
  [Rc, flag] = chol (V);
  if (flag != 0)
    d = Inf;
    return;
  endif
  ViX = Rc \ (Rc' \ X);
  Viy = Rc \ (Rc' \ y);
  XtViX = X' * ViX;
  beta = XtViX \ (X' * Viy);
  r = y - X * beta;
  rVir = r' * (Rc \ (Rc' \ r));
  logdetV = 2 * sum (log (diag (Rc)));
  if (is_reml)
    d = logdetV + rVir + 2*sum (log (diag (chol (XtViX)))) + (n-p)*log (2*pi);
  else
    d = logdetV + rVir + n*log (2*pi);
  endif
endfunction

function C = fixed_cov (eta, qk, nlev, X, Zx, n)
  V = build_V (eta, qk, nlev, Zx, n);
  C = inv (X' * (V \ X));
endfunction
