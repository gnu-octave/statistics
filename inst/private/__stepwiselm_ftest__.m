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
## @deftypefn {Private Function} {[@var{fstat}, @var{pval}, @var{sse_new}, @
##   @var{mse_new}, @var{df_term}, @var{crit_new}] =} @
##   __stepwiselm_ftest__ (@var{raw_X}, @var{raw_y}, @var{weights}, @
##   @var{incl_mask}, @var{current_terms}, @var{candidate}, @
##   @var{cat_flags}, @var{var_names}, @var{pred_names}, @
##   @var{raw_cols}, @var{has_intercept}, @var{criterion})
##
## Compute the partial F-test statistic for a candidate term action.
##
## This function builds the design matrix for the candidate model using
## @code{__fitlm_build_design__} and calls @code{lm_fit_engine} directly
## (not @code{fitlm}) for maximum performance inside the stepwise loop.
##
## @strong{F-test formula (confirmed §7)}
##
## For @strong{adding} a term:
## @example
## F = (SSE_current - SSE_augmented) / df_term / MSE_augmented
## p = 1 - fcdf(F, df_term, DFE_augmented)
## @end example
##
## For @strong{removing} a term:
## @example
## F = (SSE_reduced - SSE_current) / abs(df_term) / MSE_current
## p = 1 - fcdf(F, abs(df_term), DFE_current)
## @end example
##
## Key rule: the denominator always uses MSE of the @strong{larger}
## (more complex) model.
##
## @seealso{stepwiselm, __stepwiselm_candidates__, lm_fit_engine}
## @end deftypefn

function [fstat, pval, sse_new, mse_new, df_term, crit_new] = ...
    __stepwiselm_ftest__ (raw_X, raw_y, weights, incl_mask, ...
                          current_terms, candidate, cat_flags, ...
                          var_names, pred_names, raw_cols, ...
                          has_intercept, criterion)

  n_pred = numel (pred_names);
  excl_mask = ! incl_mask;   ## __fitlm_build_design__ expects excl_mask

  ## ── Build design matrices for both current and candidate models ──────────

  if (strcmp (candidate.action, "add"))
    ## Adding term: current model is smaller, candidate model is larger
    new_terms = [current_terms; candidate.term_row];

    ## Sort by term order (number of active predictors) for clean design
    term_order = sum (new_terms(:, 1:n_pred), 2);
    [~, sort_idx] = sort (term_order);
    new_terms = new_terms(sort_idx, :);

    ## Build design for current model (smaller)
    [~, ~, X_curr, ~, ~, ~, ~, ~] = ...
      __fitlm_build_design__ (raw_X, raw_y, current_terms, var_names, ...
                              cat_flags, weights, excl_mask, pred_names, ...
                              raw_cols, has_intercept);

    ## Build design for augmented model (larger)
    [~, ~, X_new, ~, ~, ~, ~, ~] = ...
      __fitlm_build_design__ (raw_X, raw_y, new_terms, var_names, ...
                              cat_flags, weights, excl_mask, pred_names, ...
                              raw_cols, has_intercept);

    ## Extract included rows
    X_curr_incl = X_curr(incl_mask, :);
    X_new_incl  = X_new(incl_mask, :);
    y_incl      = raw_y(incl_mask);
    w_incl      = weights(incl_mask);

    ## Fit both models
    warning ("off", "CompactLinearModel:rankDeficient");
    s_curr = lm_fit_engine (X_curr_incl, y_incl, w_incl);
    s_new  = lm_fit_engine (X_new_incl, y_incl, w_incl);
    warning ("on", "CompactLinearModel:rankDeficient");

    ## df_term = difference in design matrix columns (correct for categoricals)
    df_term = columns (X_new_incl) - columns (X_curr_incl);
    if (df_term <= 0)
      df_term = 1;  ## safety fallback
    endif

    ## F = (SSE_curr - SSE_new) / df / MSE_new
    ## Denominator is MSE of the LARGER (augmented) model
    if (s_new.mse > 0 && df_term > 0)
      fstat = (s_curr.sse - s_new.sse) / df_term / s_new.mse;
      pval  = 1 - fcdf (fstat, df_term, s_new.dfe);
    else
      fstat = 0;
      pval  = 1;
    endif

    sse_new = s_new.sse;
    mse_new = s_new.mse;
    crit_new = __compute_criterion__ (s_new, criterion);

  else
    ## Removing term: current model is larger, candidate model is smaller
    ## Build reduced model (without the term)
    red_terms = __remove_term_row__ (current_terms, candidate.term_row);

    ## Build design for current model (larger)
    [~, ~, X_curr, ~, ~, ~, ~, ~] = ...
      __fitlm_build_design__ (raw_X, raw_y, current_terms, var_names, ...
                              cat_flags, weights, excl_mask, pred_names, ...
                              raw_cols, has_intercept);

    ## Build design for reduced model (smaller)
    [~, ~, X_red, ~, ~, ~, ~, ~] = ...
      __fitlm_build_design__ (raw_X, raw_y, red_terms, var_names, ...
                              cat_flags, weights, excl_mask, pred_names, ...
                              raw_cols, has_intercept);

    ## Extract included rows
    X_curr_incl = X_curr(incl_mask, :);
    X_red_incl  = X_red(incl_mask, :);
    y_incl      = raw_y(incl_mask);
    w_incl      = weights(incl_mask);

    ## Fit both models
    warning ("off", "CompactLinearModel:rankDeficient");
    s_curr = lm_fit_engine (X_curr_incl, y_incl, w_incl);
    s_red  = lm_fit_engine (X_red_incl, y_incl, w_incl);
    warning ("on", "CompactLinearModel:rankDeficient");

    ## df_term = difference in design matrix columns
    df_term = columns (X_curr_incl) - columns (X_red_incl);
    if (df_term <= 0)
      df_term = 1;
    endif

    ## F = (SSE_red - SSE_curr) / df / MSE_curr
    ## Denominator is MSE of the LARGER (current) model
    if (s_curr.mse > 0 && df_term > 0)
      fstat = (s_red.sse - s_curr.sse) / df_term / s_curr.mse;
      pval  = 1 - fcdf (fstat, df_term, s_curr.dfe);
    else
      fstat = 0;
      pval  = 1;
    endif

    sse_new = s_red.sse;
    mse_new = s_red.mse;
    crit_new = __compute_criterion__ (s_red, criterion);
  endif

endfunction


## ── Helper: remove a term row from a terms matrix ───────────────────────────

function new_terms = __remove_term_row__ (terms, term_row)
  ## Remove all rows matching term_row from terms matrix.
  mask = true (rows (terms), 1);
  for j = 1:rows (terms)
    if (isequal (terms(j, :), term_row))
      mask(j) = false;
      break;   ## remove only the first match
    endif
  endfor
  new_terms = terms(mask, :);
endfunction


## ── Helper: compute criterion value from lm_fit_engine output ───────────────

function val = __compute_criterion__ (s, criterion)
  ## Compute model selection criterion from lm_fit_engine struct.
  switch (criterion)
    case "sse"
      val = s.sse;
    case "aic"
      val = s.crit.AIC;
    case "bic"
      val = s.crit.BIC;
    case "rsquared"
      val = -s.rsq.Ordinary;   ## negate so lower = better
    case "adjrsquared"
      val = -s.rsq.Adjusted;
    otherwise
      val = s.sse;
  endswitch
endfunction
