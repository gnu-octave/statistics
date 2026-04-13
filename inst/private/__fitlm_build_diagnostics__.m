## Copyright (C) 2026 Jayant Chauhan <0001jayant@gmail.com>
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
## @deftypefn {Private Function} {[@var{Fitted}, @var{Residuals}, @var{Diagnostics}, @
##   @var{MFVNM}] =} __fitlm_build_diagnostics__ (@var{s}, @var{X_full}, @
##   @var{y_full}, @var{incl_mask}, @var{excl_mask})
##
## Expand lm_fit_engine output from n_used rows to full n rows.
##
## Key behaviors confirmed from MATLAB Phase 4 cross-validation:
##   - Excluded rows: Fitted computed (X*beta), Residuals = NaN
##   - Missing rows (NaN in X or y): Fitted = NaN, Residuals = NaN
##   - Diagnostics: Leverage=0, CooksD/Dffits/S2_i/CovRatio=NaN for non-included rows
##
## @var{s} is the struct returned by @code{lm_fit_engine}.
## @var{X_full} is the full n x p_tot design matrix.
## @var{y_full} is the full n x 1 response (NaN for missing rows).
## @var{incl_mask} is n x 1 logical: rows used in fit.
## @var{excl_mask} is n x 1 logical: rows explicitly excluded.
## @end deftypefn

function [Fitted, Residuals, Diagnostics, MFVNM] = ...
    __fitlm_build_diagnostics__ (s, X_full, y_full, incl_mask, excl_mask)

  n     = rows (X_full);
  p_est = s.p_est;
  mse   = s.mse;
  dfe   = s.dfe;

  ## ── Fitted ──────────────────────────────────────────────────────────────
  Fitted = NaN (n, 1);

  ## Included rows
  Fitted(incl_mask) = s.yhat;

  ## Excluded-but-not-missing rows: compute X * beta
  missing_mask = ! incl_mask & ! excl_mask;   ## true = missing (NaN in data)
  excl_not_missing = excl_mask & ! missing_mask;
  if (any (excl_not_missing))
    Fitted(excl_not_missing) = X_full(excl_not_missing, :) * s.beta;
  endif
  ## Missing rows: Fitted stays NaN

  ## ── Residuals (n x 4 struct) ────────────────────────────────────────────
  raw_r  = NaN (n, 1);
  pear_r = NaN (n, 1);
  stud_r = NaN (n, 1);
  std_r  = NaN (n, 1);

  ## Leverage for included rows
  h_inc = s.h;   ## n_used x 1

  ## Pearson = raw / sqrt(MSE)
  ## Standardized = raw / sqrt(MSE * (1 - h))
  ## Studentized (externally) = raw / sqrt(S2_i * (1 - h))
  ##   where S2_i = (dfe * MSE - raw^2 / (1-h)) / (dfe - 1)

  raw_inc = s.resid;

  raw_r(incl_mask)  = raw_inc;
  pear_r(incl_mask) = raw_inc / sqrt (mse);
  std_denom         = sqrt (mse * (1 - h_inc));
  std_r(incl_mask)  = raw_inc ./ std_denom;

  ## Studentized (leave-one-out MSE)
  s2_i_inc = max (0, (dfe * mse - raw_inc .^ 2 ./ (1 - h_inc)) / max (1, dfe - 1));
  stud_denom = sqrt (s2_i_inc .* (1 - h_inc));
  stud_r_inc = raw_inc ./ stud_denom;
  stud_r_inc(stud_denom == 0) = 0;
  stud_r(incl_mask) = stud_r_inc;

  Residuals.Raw          = raw_r;
  Residuals.Pearson      = pear_r;
  Residuals.Studentized  = stud_r;
  Residuals.Standardized = std_r;

  ## ── Diagnostics (n x 7 struct) ──────────────────────────────────────────
  lev        = zeros (n, 1);
  cooks      = NaN   (n, 1);
  dffits_v   = NaN   (n, 1);
  s2_i_v     = NaN   (n, 1);
  covratio_v = NaN   (n, 1);
  dfbetas_v  = zeros (n, p_est);
  hatmat     = zeros (n, n);

  ## Build hat matrix block for included rows
  Q_inc   = s.Q;             ## n_used x p_est
  H_inc   = Q_inc * Q_inc';  ## n_used x n_used
  h_diag  = h_inc;

  lev(incl_mask) = h_diag;

  ## Hat matrix: fill in the included rows/cols
  inc_idx = find (incl_mask);
  for ii = 1:numel (inc_idx)
    for jj = 1:numel (inc_idx)
      hatmat(inc_idx(ii), inc_idx(jj)) = H_inc(ii, jj);
    endfor
  endfor

  ## Cook's distance: D_i = (raw_i^2 * h_i) / (p_est * MSE * (1-h_i)^2)
  h_inc_safe = max (h_inc, eps);
  cooks_inc = (raw_inc .^ 2 .* h_inc) ./ (p_est * mse * (1 - h_inc) .^ 2);
  cooks(incl_mask) = cooks_inc;

  ## DFFITS: raw / sqrt(S2_i * h)
  dffits_inc = raw_inc ./ sqrt (s2_i_inc .* h_inc_safe);
  dffits_inc(h_inc == 0) = 0;
  dffits_v(incl_mask) = dffits_inc;

  ## S2_i (leave-one-out variance)
  s2_i_v(incl_mask) = s2_i_inc;

  ## CovRatio: (1 / (1-h)) ^ p_est * (dfe / (dfe - 1) + raw^2 / (mse * (dfe-1) * (1-h)))^(-p_est)
  ## Simpler: (S2_i / MSE)^p_est / (1-h)
  s2_ratio = s2_i_inc / max (mse, eps);
  covrat_inc = (s2_ratio .^ p_est) ./ (1 - h_inc);
  covrat_inc(h_inc >= 1) = NaN;
  covratio_v(incl_mask) = covrat_inc;

  ## DFBETAS: (beta_full - beta_(-i)) / sqrt(S2_i * diag(XtX_inv))
  ## Approximation via s.Q: dfbetas_i = (Q(i,:)' * resid_i) / (1-h_i) / sqrt(S2_i * diag(vcov))
  vcov_diag = diag (s.vcov);           ## p_est x 1  (from engine vcov sized p_tot x p_tot)
  if (numel (vcov_diag) > p_est)
    vcov_diag = vcov_diag(1:p_est);
  endif
  se_beta = sqrt (max (vcov_diag, 0));  ## p_est x 1

  for ii = 1:numel (inc_idx)
    ri    = raw_inc(ii);
    hi    = h_inc(ii);
    s2i   = s2_i_inc(ii);
    qi    = Q_inc(ii, :)';           ## p_est x 1
    denom = (1 - hi) .* sqrt (s2i) .* se_beta;
    denom(denom == 0) = NaN;
    dfb   = (ri ./ (1 - hi)) * qi ./ sqrt (s2i .* max (vcov_diag, eps));
    dfbetas_v(inc_idx(ii), 1:p_est) = dfb';
  endfor

  Diagnostics.Leverage     = lev;
  Diagnostics.CooksDistance = cooks;
  Diagnostics.Dffits       = dffits_v;
  Diagnostics.S2_i         = s2_i_v;
  Diagnostics.CovRatio     = covratio_v;
  Diagnostics.Dfbetas      = dfbetas_v;
  Diagnostics.HatMatrix    = hatmat;

  ## ── ModelFitVsNullModel ─────────────────────────────────────────────────
  ## Check if model has intercept (intercept row = zeros(1,p) in terms_pred)
  ## We detect via lp_str starting with '1'
  has_int = (s.p_est > 0) && (s.ssr > 0 || s.sse > 0);

  ## If model has intercept: compute F-stat and p-value
  ## If no intercept: all NaN (confirmed Block 10-P4)
  ##
  ## We check via DFR = p_est - 1 if intercept, else p_est
  ## However we don't have has_intercept here. Use heuristic:
  %% detect from sst — if sst == sse+ssr and dfe = n-p_est, we can check
  %% whether the design matrix has an intercept by seeing if beta(1) was ~unconstrained
  %% Simplest: pass intercept flag when calling this function.
  %% Since we don't have it here, let the caller override MFVNM.

  MFVNM.Fstat     = NaN;
  MFVNM.Pvalue    = NaN;
  MFVNM.NullModel = NaN;

endfunction
