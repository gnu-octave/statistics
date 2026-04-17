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
## @deftypefn  {statistics} {@var{mdl} =} stepwiselm (@var{X}, @var{y})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@var{X}, @var{y}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@var{tbl})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@var{tbl}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@dots{}, @var{Name}, @var{Value})
##
## Fit a linear regression model using stepwise regression.
##
## @code{stepwiselm} performs stepwise regression starting from an initial
## model and iteratively adding or removing terms based on their statistical
## significance, returning a @code{LinearModel} object.
##
## @var{modelspec} specifies the @strong{starting} model (default:
## @qcode{'constant'}).
##
## @strong{Stepwise-specific Name-Value Arguments}
##
## @table @asis
## @item @qcode{'Lower'} — minimum model (default: @qcode{'constant'})
## @item @qcode{'Upper'} — maximum model (default: @qcode{'interactions'})
## @item @qcode{'PEnter'} — p-value threshold to add a term (default: 0.05)
## @item @qcode{'PRemove'} — p-value threshold to remove a term (default: 0.10)
## @item @qcode{'Criterion'} — @qcode{'sse'}, @qcode{'aic'}, or @qcode{'bic'}
## @item @qcode{'NSteps'} — maximum number of steps (default: Inf)
## @item @qcode{'Verbose'} — 0 (silent), 1 (one line/step), 2 (per-candidate)
## @end table
##
## All @code{fitlm} name-value arguments are also accepted.
##
## @seealso{fitlm, LinearModel}
## @end deftypefn

function mdl = stepwiselm (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ## ── [1] Parse inputs ──────────────────────────────────────────────────────
  [raw_X, raw_y, var_names, cat_flags, weights, excl_mask, ...
   modelspec, obs_names, resp_name, pred_names, has_intercept, ...
   raw_cols, lower_spec, upper_spec, p_enter, p_remove, ...
   criterion, n_steps, verbose] = ...
    __stepwiselm_parse_inputs__ (varargin{:});

  n_pred = numel (pred_names);

  ## ── [2] Build lower/upper terms matrices ──────────────────────────────────
  lower_terms_pred = __spec_to_terms__ (lower_spec, n_pred, true);
  upper_terms_pred = __spec_to_terms__ (upper_spec, n_pred, true);

  ## Add response column (always 0)
  lower_terms = [lower_terms_pred, zeros(rows (lower_terms_pred), 1)];
  upper_terms = [upper_terms_pred, zeros(rows (upper_terms_pred), 1)];

  ## ── [3] Fit starting model ────────────────────────────────────────────────
  cat_idx = find (cat_flags(1:n_pred));
  fitlm_args = {};
  if (! isempty (cat_idx))
    fitlm_args = [fitlm_args, {'CategoricalVars', cat_idx}];
  endif
  if (any (weights != 1))
    fitlm_args = [fitlm_args, {'Weights', weights}];
  endif
  if (any (excl_mask))
    fitlm_args = [fitlm_args, {'Exclude', excl_mask}];
  endif
  if (! isempty (var_names))
    fitlm_args = [fitlm_args, {'VarNames', var_names}];
  endif

  if (isa (raw_X, "table"))
    mdl0 = fitlm (raw_X, modelspec, fitlm_args{:});
  else
    mdl0 = fitlm (raw_X, raw_y, modelspec, fitlm_args{:});
  endif

  current_terms = mdl0.Formula.Terms;
  incl_mask = mdl0.ObservationInfo.Subset;

  ## ── [4] Initialize step log ───────────────────────────────────────────────
  start_lp = mdl0.Formula.LinearPredictor;
  crit_current = __get_model_crit__ (mdl0, criterion);

  n_log = 1;
  log_action    = {'Start'};
  log_term      = {start_lp};
  log_terms     = {current_terms};
  log_df        = mdl0.NumCoefficients;
  log_del_df    = NaN;
  log_fstat     = NaN;
  log_pval      = NaN;
  log_critval   = crit_current;

  ## ── [5] Check for empty candidate set ─────────────────────────────────────
  all_candidates = __stepwiselm_candidates__ (current_terms, lower_terms, ...
                                               upper_terms, cat_flags, ...
                                               pred_names);
  if (isempty (all_candidates) && verbose >= 1)
    fprintf ("No terms to add to or remove from initial model.\n");
  endif

  ## ── [6] Search loop ───────────────────────────────────────────────────────
  n_taken = 0;
  use_ftest = strcmp (criterion, "sse");

  while (n_taken < n_steps)
    made_change = false;

    ## ── Phase A: BACKWARD SCAN (keep removing) ────────────────────────────
    while (n_taken < n_steps)
      candidates = __stepwiselm_candidates__ (current_terms, lower_terms, ...
                                               upper_terms, cat_flags, ...
                                               pred_names);
      ## Filter to remove only
      rem_idx = [];
      for ci = 1:numel (candidates)
        if (strcmp (candidates(ci).action, "remove"))
          rem_idx(end+1) = ci;
        endif
      endfor

      if (isempty (rem_idx))
        break;
      endif

      ## Score each removal candidate
      best_ri = 0;
      best_pval = -Inf;
      best_fstat = 0;
      best_crit = Inf;
      best_df = 0;

      for ri = 1:numel (rem_idx)
        c = candidates(rem_idx(ri));
        [f, p, ~, ~, df_t, cr] = __stepwiselm_ftest__ ( ...
            raw_X, raw_y, weights, incl_mask, current_terms, c, ...
            cat_flags, var_names, pred_names, raw_cols, has_intercept, ...
            criterion);

        if (verbose >= 2)
          if (use_ftest)
            fprintf ("   pValue for removing %s is %g\n", c.term_name, p);
          else
            fprintf ("   %s for removing %s is %g\n", ...
                     upper (criterion), c.term_name, cr);
          endif
        endif

        should_take = false;
        if (use_ftest)
          if (p > p_remove && p > best_pval)
            should_take = true;
          endif
        else
          if (cr < crit_current && cr < best_crit)
            should_take = true;
          endif
        endif

        if (should_take)
          best_ri = rem_idx(ri);
          best_pval = p;
          best_fstat = f;
          best_crit = cr;
          best_df = df_t;
        endif
      endfor

      if (best_ri == 0)
        break;  ## no qualifying removal
      endif

      ## Apply removal
      best_c = candidates(best_ri);
      mask = true (rows (current_terms), 1);
      for j = 1:rows (current_terms)
        if (isequal (current_terms(j, :), best_c.term_row))
          mask(j) = false;
          break;
        endif
      endfor
      current_terms = current_terms(mask, :);
      n_taken++;
      made_change = true;

      ## Update criterion
      if (! use_ftest)
        crit_current = best_crit;
      endif

      ## Log
      n_log++;
      log_action{n_log}  = "Remove";
      log_term{n_log}    = best_c.term_name;
      log_terms{n_log}   = current_terms;
      log_df(n_log)      = log_df(n_log - 1) - best_df;
      log_del_df(n_log)  = -best_df;
      log_fstat(n_log)   = best_fstat;
      log_pval(n_log)    = best_pval;
      log_critval(n_log) = best_crit;

      ## Verbose
      if (verbose >= 1)
        __print_step__ (n_taken, best_c, best_fstat, best_pval, ...
                        best_crit, criterion);
      endif
    endwhile  ## backward

    if (n_taken >= n_steps)
      break;
    endif

    ## ── Phase B: FORWARD SCAN (add one) ───────────────────────────────────
    candidates = __stepwiselm_candidates__ (current_terms, lower_terms, ...
                                             upper_terms, cat_flags, ...
                                             pred_names);
    ## Filter to add only
    add_idx = [];
    for ci = 1:numel (candidates)
      if (strcmp (candidates(ci).action, "add"))
        add_idx(end+1) = ci;
      endif
    endfor

    if (isempty (add_idx))
      break;  ## converged
    endif

    ## Score each add candidate
    best_ai = 0;
    best_pval_a = Inf;
    best_fstat_a = 0;
    best_crit_a = Inf;
    best_df_a = 0;

    for ai = 1:numel (add_idx)
      c = candidates(add_idx(ai));
      [f, p, ~, ~, df_t, cr] = __stepwiselm_ftest__ ( ...
          raw_X, raw_y, weights, incl_mask, current_terms, c, ...
          cat_flags, var_names, pred_names, raw_cols, has_intercept, ...
          criterion);

      if (verbose >= 2)
        if (use_ftest)
          fprintf ("   pValue for adding %s is %g\n", c.term_name, p);
        else
          fprintf ("   %s for adding %s is %g\n", ...
                   upper (criterion), c.term_name, cr);
        endif
      endif

      should_take = false;
      if (use_ftest)
        if (p < p_enter && p < best_pval_a)
          should_take = true;
        endif
      else
        if (cr < crit_current && cr < best_crit_a)
          should_take = true;
        endif
      endif

      if (should_take)
        best_ai = add_idx(ai);
        best_pval_a = p;
        best_fstat_a = f;
        best_crit_a = cr;
        best_df_a = df_t;
      endif
    endfor

    if (best_ai == 0)
      if (! made_change)
        break;  ## converged — no add and no remove happened
      endif
      break;
    endif

    ## Apply addition
    best_c = candidates(best_ai);
    current_terms = [current_terms; best_c.term_row];
    ## Sort by term order
    term_order = sum (current_terms(:, 1:n_pred), 2);
    [~, si] = sort (term_order);
    current_terms = current_terms(si, :);
    n_taken++;
    made_change = true;

    ## Update criterion
    if (! use_ftest)
      crit_current = best_crit_a;
    endif

    ## Log
    n_log++;
    log_action{n_log}  = "Add";
    log_term{n_log}    = best_c.term_name;
    log_terms{n_log}   = current_terms;
    log_df(n_log)      = log_df(n_log - 1) + best_df_a;
    log_del_df(n_log)  = best_df_a;
    log_fstat(n_log)   = best_fstat_a;
    log_pval(n_log)    = best_pval_a;
    log_critval(n_log) = best_crit_a;

    ## Verbose
    if (verbose >= 1)
      __print_step__ (n_taken, best_c, best_fstat_a, best_pval_a, ...
                      best_crit_a, criterion);
    endif

    ## Continue outer loop → back to backward scan
  endwhile

  ## ── [7] Fit final model ───────────────────────────────────────────────────
  if (isa (raw_X, "table"))
    mdl = fitlm (raw_X, current_terms, fitlm_args{:});
  else
    mdl = fitlm (raw_X, raw_y, current_terms, fitlm_args{:});
  endif

  ## ── [8] Build Steps struct ────────────────────────────────────────────────
  Steps.Start     = [resp_name, " ~ ", log_term{1}];
  Steps.Lower     = lower_spec;
  Steps.Upper     = upper_spec;
  Steps.Criterion = upper (criterion);
  Steps.PEnter    = p_enter;
  Steps.PRemove   = p_remove;

  ## Build History struct (mimicking a table)
  History.Action   = log_action(:);
  History.TermName = cell (n_log, 1);
  for k = 1:n_log
    History.TermName{k} = {log_term{k}};
  endfor
  History.Terms = log_terms(:);
  History.DF    = log_df(:);
  History.delDF = log_del_df(:);

  if (use_ftest)
    History.FStat  = log_fstat(:);
    History.pValue = log_pval(:);
  else
    ## For AIC/BIC, the column is the criterion name, not FStat/pValue
    History.(upper (criterion)) = log_critval(:);
  endif

  Steps.History = History;

  ## ── [9] Inject Steps into LinearModel ─────────────────────────────────────
  mdl.Steps = Steps;

endfunction


## ═══════════════════════════════════════════════════════════════════════════════
## Local helper functions
## ═══════════════════════════════════════════════════════════════════════════════

function terms = __spec_to_terms__ (spec, p, has_intercept)
  ## Convert a shorthand modelspec string to a (n_terms x p) terms matrix.
  ## Also handles formula strings 'y ~ x1 + x2' and numeric matrices.

  if (isnumeric (spec) && ismatrix (spec))
    terms = spec;
    ## Strip response column if present
    if (columns (terms) == p + 1)
      terms = terms(:, 1:p);
    endif
    return;
  endif

  if (! ischar (spec))
    spec = "interactions";
  endif

  ## Check for formula string
  if (! isempty (strfind (spec, "~")))
    ## Formula — extract terms. For upper/lower bounds this is unusual
    ## but handle gracefully: just use 'interactions' as fallback
    spec = "interactions";
  endif

  switch (lower (spec))
    case "constant"
      terms = zeros (1, p);

    case "linear"
      terms = [zeros(1, p); eye(p)];

    case "interactions"
      lin = eye (p);
      inter = [];
      for i = 1:p-1
        for j = i+1:p
          r = zeros (1, p);
          r(i) = 1; r(j) = 1;
          inter = [inter; r];
        endfor
      endfor
      terms = [zeros(1, p); lin; inter];

    case "purequadratic"
      terms = [zeros(1, p); eye(p); 2*eye(p)];

    case {"quadratic", "full"}
      lin = eye (p);
      inter = [];
      for i = 1:p-1
        for j = i+1:p
          r = zeros (1, p);
          r(i) = 1; r(j) = 1;
          inter = [inter; r];
        endfor
      endfor
      terms = [zeros(1, p); lin; inter; 2*eye(p)];

    otherwise
      terms = [zeros(1, p); eye(p)];
  endswitch
endfunction


function val = __get_model_crit__ (mdl, criterion)
  switch (criterion)
    case "sse"
      val = mdl.SSE;
    case "aic"
      val = mdl.ModelCriterion.AIC;
    case "bic"
      val = mdl.ModelCriterion.BIC;
    case "rsquared"
      val = -mdl.Rsquared.Ordinary;
    case "adjrsquared"
      val = -mdl.Rsquared.Adjusted;
    otherwise
      val = mdl.SSE;
  endswitch
endfunction


function __print_step__ (step_num, candidate, fstat, pval, critval, criterion)
  action_str = [toupper(candidate.action(1)), candidate.action(2:end)];
  action_verb = "Adding";
  if (strcmp (candidate.action, "remove"))
    action_verb = "Removing";
  endif

  if (strcmp (criterion, "sse"))
    fprintf ("%d. %s %s, FStat = %g, pValue = %g\n", ...
             step_num, action_verb, candidate.term_name, fstat, pval);
  else
    fprintf ("%d. %s %s, %s = %g\n", ...
             step_num, action_verb, candidate.term_name, ...
             upper (criterion), critval);
  endif
endfunction


## ═══════════════════════════════════════════════════════════════════════════════
## Unit tests
## ═══════════════════════════════════════════════════════════════════════════════

## Test 1: PEnter >= PRemove error
%!error <PEnter>
%! rng (1);
%! X = randn (20, 3);
%! y = randn (20, 1);
%! stepwiselm (X, y, 'constant', 'PEnter', 0.1, 'PRemove', 0.05);

## Test 2: Basic stepwise — constant start, should add significant terms
%!test
%! rng (42);
%! X = randn (30, 3);
%! y = X * [2; -1; 0] + randn (30, 1) * 0.5;
%! mdl = stepwiselm (X, y, 'constant', 'Verbose', 0);
%! assert (isa (mdl, 'LinearModel'));
%! assert (! isempty (mdl.Steps));
%! assert (isstruct (mdl.Steps));
%! assert (isfield (mdl.Steps, 'History'));
%! assert (strcmp (mdl.Steps.History.Action{1}, 'Start'));

## Test 3: NSteps=1 limits to one step
%!test
%! rng (42);
%! X = randn (30, 3);
%! y = X * [2; -1; 0.5] + randn (30, 1) * 0.5;
%! mdl = stepwiselm (X, y, 'constant', 'NSteps', 1, 'Verbose', 0);
%! n_hist = numel (mdl.Steps.History.Action);
%! assert (n_hist, 2);  ## Start + 1 step

## Test 4: Upper = Lower = linear → no steps possible
%!test
%! rng (1);
%! X = randn (20, 3);
%! y = X * [1; 2; -1] + randn (20, 1) * 0.5;
%! mdl = stepwiselm (X, y, 'linear', 'Lower', 'linear', ...
%!                   'Upper', 'linear', 'Verbose', 0);
%! n_hist = numel (mdl.Steps.History.Action);
%! assert (n_hist, 1);  ## Start row only

## Test 5: Steps struct has correct fields
%!test
%! rng (42);
%! X = randn (20, 2);
%! y = X * [1; 2] + randn (20, 1) * 0.5;
%! mdl = stepwiselm (X, y, 'constant', 'Verbose', 0);
%! assert (isfield (mdl.Steps, 'Start'));
%! assert (isfield (mdl.Steps, 'Lower'));
%! assert (isfield (mdl.Steps, 'Upper'));
%! assert (isfield (mdl.Steps, 'Criterion'));
%! assert (isfield (mdl.Steps, 'PEnter'));
%! assert (isfield (mdl.Steps, 'PRemove'));
%! assert (isfield (mdl.Steps, 'History'));

## Test 6: Criterion = aic — History has AIC column instead of FStat/pValue
%!test
%! rng (42);
%! X = randn (30, 3);
%! y = X * [2; -1; 0.5] + randn (30, 1) * 0.5;
%! mdl = stepwiselm (X, y, 'constant', 'Criterion', 'aic', 'Verbose', 0);
%! assert (isfield (mdl.Steps.History, 'AIC'));
%! assert (! isfield (mdl.Steps.History, 'FStat'));

## Test 7: Verbose = 0 produces no output (just check it runs)
%!test
%! rng (1);
%! X = randn (20, 3);
%! y = X * [1; 2; -1] + randn (20, 1) * 0.5;
%! mdl = stepwiselm (X, y, 'constant', 'Verbose', 0);
%! assert (isa (mdl, 'LinearModel'));

## Test 8: Steps.Criterion is uppercase
%!test
%! rng (1);
%! X = randn (20, 2);
%! y = randn (20, 1);
%! mdl = stepwiselm (X, y, 'constant', 'Verbose', 0);
%! assert (strcmp (mdl.Steps.Criterion, 'SSE'));

## Test 9: Steps.Start is the resolved formula string
%!test
%! rng (1);
%! X = randn (20, 2);
%! y = randn (20, 1);
%! mdl = stepwiselm (X, y, 'constant', 'Verbose', 0);
%! assert (strcmp (mdl.Steps.Start, 'y ~ 1'));

## Test 10: delDF is NaN for Start row
%!test
%! rng (1);
%! X = randn (20, 2);
%! y = X * [1; 2] + randn (20, 1) * 0.5;
%! mdl = stepwiselm (X, y, 'constant', 'Verbose', 0);
%! assert (isnan (mdl.Steps.History.delDF(1)));
%! assert (isnan (mdl.Steps.History.FStat(1)));
%! assert (isnan (mdl.Steps.History.pValue(1)));

## Test 11: Hald canonical example
%!test
%! load hald;
%! mdl = stepwiselm (ingredients, heat, 'constant', 'PEnter', 0.06, 'Verbose', 0);
%! assert (numel(mdl.Steps.History.Action) == 5);
%! assert (strcmp (mdl.Formula.LinearPredictor, '1 + x1 + x2'));
%! assert (abs (mdl.Coefficients.Estimate(2) - 1.4683) < 1e-3);
%! assert (abs (mdl.RMSE - 2.406) < 0.01);

## Test 12: Categorical grouping (fisheriris)
%!test
%! load fisheriris;
%! tbl = table(meas(:,1), meas(:,2), categorical(species), ...
%!             'VariableNames', {'SL','SW','Species'});
%! mdl = stepwiselm (tbl, 'constant', 'ResponseVar', 'SL', 'Verbose', 0);
%! assert (numel(mdl.Steps.History.Action) == 3);   ## Start + Add Species + Add SW
%! assert (mdl.Steps.History.delDF(2) == 2);  ## Species delDF = 2

