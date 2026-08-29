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
## @deftypefn {private} {@var{M} =} __glmefit__ (@var{X}, @var{y}, @var{Z}, @var{G}, @var{distr}, @var{link}, @var{method})
##
## Estimation core for generalized linear mixed-effects models.  Undocumented
## internal function for @code{fitglme}/@code{GeneralizedLinearMixedModel}.
##
## Fits the model by penalized quasi-likelihood (PQL): it alternates between
## forming the generalized-linear working response and iterative weights, and
## fitting a weighted linear mixed model to that pseudo-response, until
## convergence.  @var{method} @qcode{"mpl"} uses a maximum-likelihood inner fit,
## @qcode{"rempl"} a restricted-maximum-likelihood inner fit; @qcode{"laplace"}
## and @qcode{"approximatelaplace"} use the same point estimates but report the
## Laplace-approximated marginal log-likelihood.
##
## @end deftypefn

function M = __glmefit__ (X, y, Z, G, distr, link, method)

  n = rows (X);
  p = columns (X);
  method = lower (method);
  distr = lower (distr);
  link = lower (link);

  ## Build the expanded random design + per-term bookkeeping (cf. __lmefit__).
  [Zx, qk, nlev, levels, gidx] = expand_random (Z, G, n);

  ## Fixed dispersion for binomial/poisson; estimated for normal/gamma/invgauss.
  fixed_disp = any (strcmp (distr, {"binomial", "poisson"}));
  isreml = strcmp (method, "rempl");

  ## --- initialise beta from a plain GLM (no random effects) ---
  beta = glm_irls (X, y, distr, link);
  b = zeros (columns (Zx), 1);
  theta = init_theta (qk);
  disp0 = 1;

  ## --- PQL outer iterations ---
  for outer = 1:100
    eta = X * beta + Zx * b;
    [mu, dmu] = inv_link (eta, link);
    v = var_fun (mu, distr);
    w = dmu .^ 2 ./ v;                          ## iterative weights (disp = 1)
    zt = eta + (y - mu) ./ dmu;                 ## working response

    ## inner weighted (RE)ML fit of the pseudo-response zt, Var(e)=disp/w
    [theta, disp0] = inner_fit (theta, X, zt, Zx, qk, nlev, w, isreml, ...
                                fixed_disp);
    [beta_n, b, Psi, covbeta, V] = inner_solve (theta, X, zt, Zx, qk, nlev, ...
                                                 w, disp0);
    if (max (abs (beta_n - beta)) < 1e-8)
      beta = beta_n;
      break;
    endif
    beta = beta_n;
  endfor

  ## --- final quantities ---
  eta = X * beta + Zx * b;
  [mu, dmu] = inv_link (eta, link);
  v = var_fun (mu, distr);
  w = dmu .^ 2 ./ v;
  zt = eta + (y - mu) ./ dmu;

  ## pseudo (RE)ML log-likelihood of the final weighted linear mixed model
  ## (REML uses n - p in place of n in the constant term).
  if (isreml), ndf = n - p; else, ndf = n; endif
  pll = -0.5 * (ndf * log (2*pi) ...
                + weighted_dev (theta, X, zt, Zx, qk, nlev, w, isreml, disp0));
  if (any (strcmp (method, {"laplace", "approximatelaplace"})))
    loglik = laplace_loglik (X, y, beta, Zx, qk, nlev, gidx, Psi, distr, link);
  else
    loglik = pll;
  endif

  M = struct ();
  M.beta = beta;
  M.b = b;
  M.Psi = Psi;
  M.dispersion = disp0;
  M.covbeta = covbeta;
  M.loglik = loglik;
  M.mu = mu;
  M.fitted = mu;
  M.resid_raw = y - mu;
  M.resid_pearson = (y - mu) ./ sqrt (var_fun (mu, distr));
  M.Zx = Zx;  M.qk = qk;  M.nlev = nlev;
  M.levels = levels;  M.gidx = gidx;
  M.n = n;  M.p = p;  M.dfe = n - p;
  M.distr = distr;  M.link = link;  M.method = method;

endfunction

## ---- link functions (canonical) ----
function [mu, dmu] = inv_link (eta, link)
  switch (link)
    case "logit"
      mu = 1 ./ (1 + exp (-eta));
      dmu = mu .* (1 - mu);
    case "log"
      mu = exp (eta);
      dmu = mu;
    case "identity"
      mu = eta;
      dmu = ones (size (eta));
    otherwise
      error ("__glmefit__: unsupported link '%s'.", link);
  endswitch
endfunction

function v = var_fun (mu, distr)
  switch (distr)
    case "binomial"
      v = mu .* (1 - mu);
    case "poisson"
      v = mu;
    case "normal"
      v = ones (size (mu));
    otherwise
      error ("__glmefit__: unsupported distribution '%s'.", distr);
  endswitch
endfunction

## plain GLM by iteratively reweighted least squares (no random effects)
function beta = glm_irls (X, y, distr, link)
  beta = zeros (columns (X), 1);
  for it = 1:100
    eta = X * beta;
    [mu, dmu] = inv_link (eta, link);
    w = dmu .^ 2 ./ var_fun (mu, distr);
    z = eta + (y - mu) ./ dmu;
    bn = (X' * (w .* X)) \ (X' * (w .* z));
    if (max (abs (bn - beta)) < 1e-10)
      beta = bn;
      break;
    endif
    beta = bn;
  endfor
endfunction

## ---- expanded random design (shared with the LMM engine's convention) ----
function [Zx, qk, nlev, levels, gidx] = expand_random (Z, G, n)
  nt = numel (Z);
  Zx = [];  qk = zeros (1, nt);  nlev = zeros (1, nt);
  levels = cell (1, nt);  gidx = cell (1, nt);
  for k = 1:nt
    qk(k) = columns (Z{k});
    [lev, ~, gi] = unique (G{k}(:));
    nlev(k) = numel (lev);  levels{k} = lev;  gidx{k} = gi;
    blk = zeros (n, qk(k) * nlev(k));
    for l = 1:nlev(k)
      blk (gi == l, (l-1)*qk(k) + (1:qk(k))) = Z{k}(gi == l, :);
    endfor
    Zx = [Zx, blk];
  endfor
endfunction

## ---- weighted mixed-model inner solver / objective ----
function [theta, disp0] = inner_fit (theta0, X, z, Zx, qk, nlev, w, ...
                                     isreml, fixed_disp)
  obj = @(th) weighted_dev (th, X, z, Zx, qk, nlev, w, isreml, 1);
  ## The weighted deviance is differentiable in theta and its gradient is
  ## closed form, so hand it over rather than let the optimiser difference an
  ## objective whose curvature in the variance components is far below the
  ## step a finite difference has to take.  See __lmefit__ for the same
  ## treatment of the linear case.
  opts = optimset ("TolX", 1e-10, "TolFun", 1e-10, "MaxFunEvals", 2000, ...
                   "Display", "off", "GradObj", "on");
  theta = fminunc (obj, theta0, opts);
  if (fixed_disp)
    disp0 = 1;
  else
    disp0 = profile_dispersion (theta, X, z, Zx, qk, nlev, w, isreml);
  endif
endfunction

## Weighted -2*log-likelihood at the covariance parameters theta, and its
## gradient.  Only Zx*D*Zx' in V depends on theta, so dV/dtheta is
## Zx*(dD/dtheta)*Zx'; beta drops out because it is the GLS minimiser, and
## D = L*L' gives dD = E_ij*L' + L*E_ij', so every entry is 2*(S*L)(i,j) for
## the matching accumulated S.  The dispersion is held at one here and
## profiled afterwards, so there is no scale term to differentiate through.
function [dev, grad] = weighted_dev (theta, X, z, Zx, qk, nlev, w, isreml, ...
                                     disp0)
  n = rows (X);  p = columns (X);
  V = build_V (theta, qk, nlev, Zx, w, disp0);
  [Rc, flag] = chol (V);
  if (flag != 0)
    dev = Inf;
    grad = zeros (size (theta));
    return;
  endif
  ## V is Rc'*Rc, so X'*inv (V)*X is (Rc'\X)'*(Rc'\X): whitening once halves
  ## the triangular solves and makes XtViX symmetric by construction.
  W = Rc' \ X;
  zz = Rc' \ z;
  XtViX = W' * W;
  beta = W \ zz;
  r = z - X * beta;
  u = Rc' \ r;
  dev = 2*sum (log (diag (Rc))) + u' * u;
  if (isreml)
    Rx = chol (XtViX);
    dev += 2 * sum (log (diag (Rx)));
  endif
  if (nargout < 2)
    return;
  endif

  ## One further solve carries the gradient, all of it in the random effects
  ## dimension: Vv'*Vv is Zx'*inv (V)*Zx, Vv'*u is Zx'*inv (V)*r, Vv'*W is
  ## Zx'*inv (V)*X.
  Vv = Rc' \ Zx;
  c = Vv' * u;
  B = Vv' * W;
  if (isreml)
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
      sel = col + (1:q);
      col += q;
      Vs = Vv(:,sel);
      Sa += Vs' * Vs;
      cs = c(sel);
      Sc += cs * cs';
      if (isreml)
        Ts = T(:,sel);
        Sp += Ts' * Ts;
      endif
    endfor
    Ga = 2 * (Sa * Lk);
    Gc = 2 * (Sc * Lk);
    if (isreml)
      Gp = 2 * (Sp * Lk);
    else
      Gp = zeros (q);
    endif
    idx = 0;
    for j = 1:q
      for i = j:q
        idx += 1;
        grad(off+idx) = Ga(i,j) - Gc(i,j) - Gp(i,j);
      endfor
    endfor
    off += m;
  endfor
endfunction

## The dispersion at the fitted covariance parameters: the weighted residual
## sum of squares over its degrees of freedom.  Whitened, as weighted_dev is,
## rather than through an explicit inverse of V.
function d = profile_dispersion (theta, X, z, Zx, qk, nlev, w, isreml)
  n = rows (X);  p = columns (X);
  V = build_V (theta, qk, nlev, Zx, w, 1);
  Rc = chol (V);
  W = Rc' \ X;
  zz = Rc' \ z;
  beta = W \ zz;
  r = z - X * beta;
  u = Rc' \ r;
  if (isreml)
    d = (u' * u) / (n - p);
  else
    d = (u' * u) / n;
  endif
endfunction

function [beta, b, Psi, covbeta, V] = inner_solve (theta, X, z, Zx, qk, ...
                                                   nlev, w, disp0)
  V = build_V (theta, qk, nlev, Zx, w, disp0);
  Vi = inv (V);
  XtViX = X' * Vi * X;
  beta = XtViX \ (X' * Vi * z);
  Dfull = build_Dfull (theta, qk, nlev);
  b = Dfull * (Zx' * (Vi * (z - X * beta)));
  covbeta = inv (XtViX);
  ## per-term covariance matrices
  Psi = cell (1, numel (qk));  off = 0;
  for k = 1:numel (qk)
    m = qk(k)*(qk(k)+1)/2;
    Lk = tril_from_theta (theta(off+(1:m)), qk(k));
    Psi{k} = Lk * Lk';
    off += m;
  endfor
endfunction

## V = Zx * blkdiag(Psi) * Zx' + disp * diag(1/w)
function V = build_V (theta, qk, nlev, Zx, w, disp0)
  n = rows (Zx);
  Dfull = build_Dfull (theta, qk, nlev);
  ## Zx * Dfull * Zx' is symmetric in exact arithmetic but not bitwise, since
  ## the two products are separate BLAS calls.  'chol' reads one triangle, so
  ## a threaded or blocked BLAS can disagree with a reference one on which
  ## factor comes out.  Symmetrise before it is handed to 'chol'.
  V = Zx * Dfull * Zx' + disp0 * diag (1 ./ w);
  V = (V + V') / 2;
endfunction

function Dfull = build_Dfull (theta, qk, nlev)
  blocks = {};  off = 0;
  for k = 1:numel (qk)
    m = qk(k)*(qk(k)+1)/2;
    Lk = tril_from_theta (theta(off+(1:m)), qk(k));
    Dk = Lk * Lk';
    for l = 1:nlev(k), blocks{end+1} = Dk; endfor
    off += m;
  endfor
  Dfull = blkdiag (blocks{:});
endfunction

function L = tril_from_theta (th, q)
  L = zeros (q, q);  idx = 0;
  for j = 1:q
    for i = j:q
      idx += 1;  L(i, j) = th(idx);
    endfor
  endfor
endfunction

function theta0 = init_theta (qk)
  theta0 = [];
  for k = 1:numel (qk)
    for j = 1:qk(k)
      for i = j:qk(k)
        theta0(end+1) = 0.1 * (i == j);
      endfor
    endfor
  endfor
  theta0 = theta0(:);
endfunction

## ---- Laplace-approximated marginal log-likelihood (single grouping term) ----
## Per group: logL_g = log p(y_g | b_hat) - 0.5*b_hat'*inv(P)*b_hat
##                     - 0.5*log|I + P*(Z_g'*W*Z_g)|,  with W the GLM weights at
## the mode b_hat.  The last two terms are combined into a single log-det that
## stays finite as the random-effect variance goes to zero.
function ll = laplace_loglik (X, y, beta, Zx, qk, nlev, gidx, Psi, distr, link)
  if (numel (qk) != 1)
    ll = NaN;                                   ## only single-term supported
    return;
  endif
  q = qk(1);
  P = Psi{1};
  ## A negligible random-effect variance pins the group modes at zero, so the
  ## marginal likelihood is the plain GLM likelihood (as MATLAB reports).
  degenerate = (max (abs (diag (P))) < 1e-8);
  if (! degenerate), Pi = inv (P); endif
  gi = gidx{1};
  ll = 0;
  for g = 1:nlev(1)
    idx = (gi == g);
    Xg = X(idx, :);
    Zg = Zx(idx, (g-1)*q + (1:q));
    yg = y(idx);
    bg = zeros (q, 1);
    if (! degenerate)
      for it = 1:100                            ## Newton for the group mode
        eta = Xg * beta + Zg * bg;
        [mu, dmu] = inv_link (eta, link);
        W = dmu .^ 2 ./ var_fun (mu, distr);
        gr = Zg' * (yg - mu) - Pi * bg;         ## canonical-link score
        H = -(Zg' * (W .* Zg)) - Pi;
        step = H \ gr;
        bg = bg - step;
        if (max (abs (step)) < 1e-10), break; endif
      endfor
    endif
    eta = Xg * beta + Zg * bg;
    [mu, dmu] = inv_link (eta, link);
    W = dmu .^ 2 ./ var_fun (mu, distr);
    if (degenerate)
      pen = 0;
    else
      pen = 0.5 * (bg' * (P \ bg));
    endif
    ll += log_pmf (yg, mu, distr) - pen ...
          - 0.5 * logdet_spd (eye (q) + P * (Zg' * (W .* Zg)));
  endfor
endfunction

function lp = log_pmf (y, mu, distr)
  switch (distr)
    case "binomial"
      lp = sum (y .* log (mu) + (1 - y) .* log (1 - mu));
    case "poisson"
      lp = sum (y .* log (mu) - mu - gammaln (y + 1));
    case "normal"
      lp = sum (-0.5 * (y - mu) .^ 2 - 0.5 * log (2*pi));
    otherwise
      lp = NaN;
  endswitch
endfunction

function ld = logdet_spd (A)
  ## A is built from products and is not bitwise symmetric; see build_V.
  A = (A + A') / 2;
  ld = 2 * sum (log (diag (chol (A))));
endfunction
