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
## @deftypefn {Private Function} {@var{s} =} lm_fit_engine (@var{X}, @var{y})
## @deftypefnx {Private Function} {@var{s} =} lm_fit_engine (@var{X}, @var{y}, @var{weights})
##
## QR-based ordinary least squares regression engine.
##
## This is the single source of truth for all linear regression math in the
## statistics package.  It is called by @code{fitlm} (Phase 4) and by
## @code{LinearModel} methods such as @code{addTerms} and @code{removeTerms}
## (Phase 3).
##
## @var{X} is the @math{n x p} design matrix (including intercept column if
## desired).  @var{y} is the @math{n x 1} response vector.  @var{weights} is
## an optional @math{n x 1} vector of observation weights; pass @code{[]} for
## unweighted OLS.
##
## Returns a struct @var{s} with the following fields:
##
## @table @code
## @item beta
## Coefficient estimates (p x 1).
## @item se
## Standard errors (p x 1).
## @item tstat
## t-statistics (p x 1).
## @item pval
## p-values for two-sided t-test (p x 1).
## @item vcov
## Coefficient covariance matrix (p x p).
## @item yhat
## Fitted values (n x 1).
## @item resid
## Residuals (n x 1).
## @item dfe
## Error degrees of freedom.
## @item sse
## Sum of squared errors (residual sum of squares).
## @item ssr
## Regression sum of squares.
## @item sst
## Total sum of squares (= SSE + SSR).
## @item mse
## Mean squared error (SSE / DFE).
## @item rmse
## Root mean squared error.
## @item logL
## Log-likelihood under normality assumption.
## @item p_est
## Number of estimated (non-rank-deficient) coefficients.
## @item p_tot
## Total number of coefficients in the model.
## @item rank_def
## Logical vector indicating rank-deficient columns.
## @item rsq
## Struct with fields @code{Ordinary} and @code{Adjusted}.
## @item crit
## Struct with fields @code{AIC}, @code{AICc}, @code{BIC}, @code{CAIC}.
## @item n
## Number of observations.
## @item Q
## Economy QR factor for non-deficient columns (n x p_est).  Used by
## @code{LinearModel} to compute the hat matrix without a second QR.
## @item h
## Leverage vector (n x 1), equal to @code{diag(Q*Q')} = row sums of
## squares of @var{Q}.
## @end table
##
## @seealso{CompactLinearModel}
## @end deftypefn

function s = lm_fit_engine (X, y, weights)

  if (nargin < 2)
    print_usage ();
  endif
  if (nargin < 3)
    weights = [];
  endif

  n = rows (X);
  p = columns (X);

  ## Weighted or unweighted QR
  if (isempty (weights))
    [Q, R] = qr (X, 0);
    Qty    = Q' * y;
  else
    sqw    = sqrt (weights(:));
    [Q, R] = qr (bsxfun (@times, sqw, X), 0);
    Qty    = Q' * (sqw .* y);
  endif

  ## Rank deficiency detection
  ## Detect before solving to avoid spurious singular-matrix warnings
  tol      = max ([n, p]) * eps (norm (R, "inf"));
  rank_def = abs (diag (R)) < tol;

  ## Solve for beta
  if (any (rank_def))
    warning ("CompactLinearModel:rankDeficient", ...
             "CompactLinearModel: rank deficient, rank = %d.", ...
             sum (! rank_def));
    ## Solve only the non-deficient subsystem
    idx_good = find (! rank_def);
    beta = zeros (p, 1);
    R_good = R(idx_good, idx_good);
    beta(idx_good) = R_good \ Qty(idx_good);
  else
    beta = R \ Qty;
  endif

  ## Residuals and sums of squares
  yhat  = X * beta;
  resid = y - yhat;
  p_est = sum (! rank_def);           ## NumEstimatedCoefficients
  dfe   = n - p_est;
  sse   = sum (resid .^ 2);
  ybar  = mean (y);
  ssr   = sum ((yhat - ybar) .^ 2);
  sst   = sse + ssr;                  ## Always SSE + SSR (not sum((y-ybar)^2))
  mse   = sse / dfe;
  rmse  = sqrt (mse);

  ## Coefficient covariance
  ## Use pinv-style via R to handle rank deficiency
  Rinv  = zeros (p, p);
  idx   = find (! rank_def);
  R_est = R(idx, idx);
  Rinv(idx, idx) = R_est \ eye (numel (idx));
  vcov  = (Rinv * Rinv') * mse;

  se        = sqrt (diag (vcov));
  tstat     = beta ./ se;
  tstat(rank_def) = NaN;
  se(rank_def)    = 0;
  pval      = 2 * (1 - tcdf (abs (tstat), dfe));
  pval(rank_def)  = NaN;

  ## Log-likelihood and information criteria
  ## Under normality, MLE variance estimate is sigma2 = SSE/n (not MSE = SSE/DFE)
  ## This matches MATLAB's computation of LogLikelihood.
  sigma2_mle = sse / n;
  logL = -n/2 * log (2 * pi * sigma2_mle) - sse / (2 * sigma2_mle);

  ## m = NumEstimatedCoefficients (no +1 for sigma^2)
  ## Verified: for n=150, logL=-34.787, p_est=5 -> AIC = -2*(-34.787)+2*5 = 79.574
  m    = p_est;
  aic  = -2 * logL + 2 * m;
  aicc = aic + (2 * m * (m + 1)) / (n - m - 1);
  bic  = -2 * logL + m * log (n);
  caic = -2 * logL + m * (log (n) + 1);

  ## Pack output struct
  s.beta     = beta;
  s.se       = se;
  s.tstat    = tstat;
  s.pval     = pval;
  s.vcov     = vcov;
  s.yhat     = yhat;
  s.resid    = resid;
  s.n        = n;
  s.dfe      = dfe;
  s.sse      = sse;
  s.ssr      = ssr;
  s.sst      = sst;
  s.mse      = mse;
  s.rmse     = rmse;
  s.logL     = logL;
  s.p_est    = p_est;
  s.p_tot    = p;
  s.rank_def = rank_def;
  s.rsq      = struct ("Ordinary",  ssr / sst, ...
                        "Adjusted", 1 - (sse / dfe) / (sst / (n - 1)));
  s.crit     = struct ("AIC", aic, "AICc", aicc, "BIC", bic, "CAIC", caic);

  ## QR factor and leverage for LinearModel diagnostics
  Q_est = Q(:, ! rank_def);
  s.Q   = Q_est;
  s.h   = sum (Q_est .^ 2, 2);

endfunction

## Simple 2-predictor OLS — linearly independent columns
%!test
%! X = [ones(5,1), [1;2;3;4;5], [2;1;4;3;6]];
%! y = X * [1; 2; -0.5];
%! s = lm_fit_engine (X, y);
%! ## Perfect fit (noise-free)
%! assert (s.beta, [1; 2; -0.5], 1e-10);
%! assert (s.resid, zeros (5, 1), 1e-10);
%! assert (s.p_tot, 3);
%! assert (s.n, 5);

## Known regression — y = 2 + 3*x with noise-free data
%!test
%! x = (1:10)';
%! y = 2 + 3 * x;
%! X = [ones(10,1), x];
%! s = lm_fit_engine (X, y);
%! assert (s.beta, [2; 3], 1e-10);
%! assert (s.sse, 0, 1e-10);
%! assert (s.rsq.Ordinary, 1, 1e-10);

## Rank-deficient design matrix
%!test
%! warning ("off", "CompactLinearModel:rankDeficient");
%! X = [1, 2, 4; 1, 3, 6; 1, 4, 8; 1, 5, 10];  ## col3 = 2*col2
%! y = [1; 2; 3; 4];
%! s = lm_fit_engine (X, y);
%! assert (s.p_est, 2);
%! assert (s.p_tot, 3);
%! assert (any (s.rank_def), true);
%! ## The rank-deficient coefficient should be zeroed
%! assert (s.beta(s.rank_def), 0);
%! assert (s.se(s.rank_def), 0);
%! assert (isnan (s.tstat(s.rank_def)), true);
%! assert (isnan (s.pval(s.rank_def)), true);
%! warning ("on", "CompactLinearModel:rankDeficient");

## SST = SSE + SSR identity
%!test
%! X = [ones(20,1), randn(20,2)];
%! y = X * [1; 2; -1] + 0.5 * randn (20, 1);
%! s = lm_fit_engine (X, y);
%! assert (s.sst, s.sse + s.ssr, 1e-10);
%! assert (s.dfe, 20 - 3);

## ModelCriterion — m = p_est (no +1)
## Verify AIC = -2*logL + 2*p_est
%!test
%! X = [ones(20,1), (1:20)'];
%! y = 3 + 2 * (1:20)' + randn (20, 1) * 0.1;
%! s = lm_fit_engine (X, y);
%! assert (s.crit.AIC, -2 * s.logL + 2 * s.p_est, 1e-10);
%! assert (s.crit.BIC, -2 * s.logL + s.p_est * log (s.n), 1e-10);

## Weighted regression
%!test
%! X = [ones(5,1), [1;2;3;4;5]];
%! y = [2; 4; 5; 4; 5];
%! w = [1; 1; 1; 1; 1];
%! s1 = lm_fit_engine (X, y, []);
%! s2 = lm_fit_engine (X, y, w);
%! ## Equal weights should give same result as unweighted
%! assert (s1.beta, s2.beta, 1e-10);
%! assert (s1.sse, s2.sse, 1e-10);
