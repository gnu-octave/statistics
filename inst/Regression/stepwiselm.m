## Copyright (C) 2026 Avanish Salunke <avanishsalunke16@gmail.com>
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
## @deftypefn  {statistics} {@var{mdl} =} stepwiselm (@var{tbl})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@var{tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@var{tbl}, @var{y})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@var{X}, @var{y})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@dots{}, @var{InitialModel})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@dots{}, @var{Name}, @var{Value}, @dots{})
##
## Fit a linear regression model using stepwise regression and return a
## @code{LinearModel} object.
##
## @code{stepwiselm} starts from an initial model and repeatedly searches for
## a term to add to, or remove from, the current model, based on the value of
## the @qcode{'Criterion'} option, until no single addition or removal
## improves the model any further.
##
## @subheading Basic Syntax
##
## @itemize
## @item
## @code{@var{mdl} = stepwiselm (@var{tbl})} fits a stepwise model using the
## variables in the table (or dataset) @var{tbl}, starting from a constant
## model.  By default, the last variable in @var{tbl} is used as the response
## and all other variables are candidate predictors.  Variables that are
## @code{categorical} arrays, cell arrays of character vectors, or logical
## arrays are automatically treated as categorical predictors.
## @item
## @code{@var{mdl} = stepwiselm (@var{tbl}, @var{ResponseVarName})} uses the
## variable named @var{ResponseVarName} in @var{tbl} as the response, and all
## remaining variables in @var{tbl} as candidate predictors.
## @item
## @code{@var{mdl} = stepwiselm (@var{tbl}, @var{y})} uses the variables in
## @var{tbl} as candidate predictors and the external numeric vector @var{y}
## as the response.
## @item
## @code{@var{mdl} = stepwiselm (@var{X}, @var{y})} fits a stepwise model of
## the response @var{y} to the predictor data @var{X}, an @math{N*P} numeric
## or logical matrix.  By default, the predictors are named @qcode{'x1'},
## @qcode{'x2'}, @dots{}, @qcode{'xP'} and the response is named
## @qcode{'y'}.
## @end itemize
##
## @subheading Initial Model, and Lower/Upper Bounds
##
## @code{@var{mdl} = stepwiselm (@dots{}, @var{InitialModel})} additionally
## specifies the model to start the stepwise search from, using any of the
## input combinations shown above.  @var{InitialModel} can be any of the
## following, and the same set of values can also be used for the
## @qcode{'Lower'} and @qcode{'Upper'} options below, which bound the
## smallest and largest set of terms @code{stepwiselm} is allowed to reach.
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Value} @tab @var{Description}
##
## @item @qcode{'constant'} @tab Model contains only an intercept term.  This
## is the default @var{InitialModel} and default @qcode{'Lower'} bound.
##
## @item @qcode{'linear'} @tab Model contains an intercept and one term
## for each predictor variable.
##
## @item @qcode{'interactions'} @tab Model contains an intercept, all
## linear terms, and all pairwise products of distinct predictor variables
## (no squared terms).  This is the default @qcode{'Upper'} bound.
##
## @item @qcode{'purequadratic'} @tab Model contains an intercept, all
## linear terms, and all squared terms.
##
## @item @qcode{'quadratic'} @tab Model contains an intercept, all linear
## terms, all pairwise products of distinct predictor variables, and all
## squared terms.
##
## @item @qcode{'polyijk'} @tab Model is a polynomial with maximum degree
## @math{i} in the first predictor, @math{j} in the second, and so on, given
## as a run of single-digit numerals, one per predictor (e.g.
## @qcode{'poly21'} for two predictors).  The model contains interaction
## terms, but the degree of each interaction term never exceeds the largest
## of the specified per-predictor degrees.
##
## @item terms matrix @tab A @math{T*P} or @math{T*(P+1)} numeric matrix,
## where @math{T} is the number of terms and @math{P} is the number of
## predictor variables, following the same convention as @code{fitlm}'s
## terms matrix.  When @var{InitialModel} is given as a terms matrix, the
## @qcode{'PredictorVars'} option may not also be used.
##
## @item Wilkinson formula @tab A character vector of the form
## @qcode{'y ~ terms'}.  When a formula is combined with @qcode{'ResponseVar'}
## or @qcode{'PredictorVars'}, the formula's response and predictor terms
## must agree with those options, or @code{stepwiselm} errors.
## @end multitable
##
## @subheading Options
##
## @code{@var{mdl} = stepwiselm (@dots{}, @var{Name}, @var{Value}, @dots{})}
## specifies additional options using one or more @qcode{Name-Value} pair
## arguments.
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Criterion'} @tab Criterion used to decide whether a term is
## added or removed at each step.  One of @qcode{'sse'} (default),
## @qcode{'aic'}, @qcode{'bic'}, @qcode{'rsquared'}, or @qcode{'adjrsquared'}.
## For @qcode{'sse'}, the p-value of an F-test comparing the model with and
## without the candidate term is used; for the others, the raw change in the
## named quantity is used directly against @qcode{'PEnter'}/@qcode{'PRemove'}.
##
## @item @qcode{'PEnter'} @tab Threshold to add a term.  Defaults depend on
## @qcode{'Criterion'}: @qcode{0.05} for @qcode{'sse'}, @qcode{0} for
## @qcode{'aic'}/@qcode{'bic'}, @qcode{0.1} for @qcode{'rsquared'}, @qcode{0}
## for @qcode{'adjrsquared'}.
##
## @item @qcode{'PRemove'} @tab Threshold to remove a term.  Defaults depend
## on @qcode{'Criterion'}: @qcode{0.10} for @qcode{'sse'}, @qcode{0.01} for
## @qcode{'aic'}/@qcode{'bic'}, @qcode{0.05} for @qcode{'rsquared'},
## @qcode{-0.05} for @qcode{'adjrsquared'}.
##
## @item @qcode{'NSteps'} @tab Maximum number of add/remove steps to take,
## given as a nonnegative integer.  Default is unlimited.  @qcode{'NSteps'}
## set to @qcode{0} returns the initial model unchanged.
##
## @item @qcode{'Lower'} @tab Model specification (in the same form as
## @var{InitialModel}, above) describing terms that may never be removed from
## the model.  Terms in @qcode{'Lower'} are protected from removal, but are
## not automatically added if absent from @var{InitialModel}.  Default is
## @qcode{'constant'}.
##
## @item @qcode{'Upper'} @tab Model specification (in the same form as
## @var{InitialModel}, above) describing the largest set of terms
## @code{stepwiselm} may add.  Default is @qcode{'interactions'}.
##
## @item @qcode{'Verbose'} @tab Controls how much progress information is
## printed while stepping.  @qcode{0} suppresses all output, @qcode{1}
## (default) prints the action taken at each step, @qcode{2} additionally
## prints the p-value or criterion value considered for every candidate term
## examined at each step.
##
## @item @qcode{'Intercept'} @tab A logical scalar indicating whether the
## initial model includes a constant (intercept) term.  Only applies when
## @var{InitialModel} is a character vector model name (or omitted); ignored
## when @var{InitialModel} is a terms matrix or formula.  Default is
## @qcode{true}.
##
## @item @qcode{'Weights'} @tab A numeric vector of nonnegative observation
## weights, one per observation.  Default is a vector of ones.
##
## @item @qcode{'Exclude'} @tab A numeric or logical vector specifying
## observations to exclude from the fit.
##
## @item @qcode{'CategoricalVars'} @tab Specifies which predictor variables
## are treated as categorical, given as a vector of column indices, a logical
## vector, or a cell array of variable names (table input only).  Each
## categorical predictor with @math{L} distinct categories is expanded into
## @math{L-1} indicator variables, and @code{stepwiselm} always adds or
## removes that entire group of indicator variables together, in a single
## step, treating the categorical predictor as one term.
##
## @item @qcode{'VarNames'} @tab A cell array of character vectors naming the
## predictor and response variables, response last.  Only applies when
## @var{X} and @var{y} are supplied directly, not table input.
##
## @item @qcode{'ResponseVar'} @tab A character vector naming the response
## variable, overriding the response that would otherwise be inferred (the
## last table variable, or @qcode{'y'} for matrix input).
##
## @item @qcode{'PredictorVars'} @tab A cell array of character vectors
## naming which variables in @var{tbl} to consider as candidate predictors.
## By default, all variables in @var{tbl} other than the response are used.
## May not be combined with a terms-matrix @var{InitialModel}, and must
## agree with any formula-based @var{InitialModel}.
## @end multitable
##
## @subheading Algorithm
##
## At each step, @code{stepwiselm} examines every term not currently in the
## model but within the @qcode{'Upper'} bound, and every term currently in
## the model but not protected by the @qcode{'Lower'} bound.  If any term
## outside the model would improve it by at least @qcode{'PEnter'}, the best
## such term is added; otherwise, if any term inside the model falls short of
## @qcode{'PRemove'}, the worst such term is removed.  The process repeats
## until neither an addition nor a removal improves the model, or until
## @qcode{'NSteps'} steps have been taken.
##
## @code{stepwiselm} never adds a higher-order term unless all of its
## lower-order marginal terms are already in the model (e.g. it will not add
## @code{x1:x2^2} unless both @code{x1} and @code{x2^2} are already present),
## and correspondingly never removes a lower-order term that a higher-order
## term still in the model depends on.  At every step, if a term in the
## current model is found to be exactly redundant (linearly dependent on the
## other terms already in the model), it is removed immediately regardless of
## the @qcode{'Criterion'} value.
##
## Because the final model depends on the initial model and the order in
## which terms are considered, @code{stepwiselm} finds a locally, but not
## necessarily globally, optimal model.
##
## Robust fitting cannot be combined with stepwise regression; do not pass
## @qcode{'RobustOpts'} to @code{stepwiselm}.
##
## @var{mdl} is returned as a @code{LinearModel} object.  See also the
## @code{step} method of @code{LinearModel}, which performs a single bounded
## round of stepwise search starting from an already-fitted model.
##
## @seealso{LinearModel, fitlm}
## @end deftypefn
function mdl = stepwiselm (varargin)

  if (nargin < 1)
    error ("stepwiselm: Not enough input arguments.");
  endif

  sw_keys = {'criterion','penter','premove','nsteps','verbose','lower','upper'};
  is_sw   = @(s) (ischar (s) || isstring (s)) && any (strcmpi (char (s), sw_keys));

  criterion = 'sse';
  penter    = [];
  premove   = [];
  nsteps    = Inf;
  verbose   = 1;
  lower_sp  = 'constant';
  upper_sp  = 'interactions';

  fit_args = {};
  i = 1;
  n = numel (varargin);
  while (i <= n)
    a = varargin{i};
    if (is_sw (a))
      if (i == n)
        error ("stepwiselm: Name-Value arguments must be in pairs.");
      endif
      key = lower (char (a));
      val = varargin{i+1};
      if (strcmp (key, 'criterion'))
        if (! (ischar (val) || isstring (val)) || ...
            ! any (strcmpi (char (val), {'sse','aic','bic','rsquared','adjrsquared'})))
          error (["stepwiselm: '" char(val) "' is not a valid value for the" ...
            " 'Criterion' argument.  Valid values are: 'AIC', 'BIC'," ...
            " 'Rsquared', 'AdjRsquared', 'SSE'."]);
        endif
        criterion = lower (char (val));
      elseif (strcmp (key, 'penter'))
        penter = double (val);
      elseif (strcmp (key, 'premove'))
        premove = double (val);
      elseif (strcmp (key, 'nsteps'))
        nsteps = double (val);
      elseif (strcmp (key, 'verbose'))
        verbose = double (val);
      elseif (strcmp (key, 'lower'))
        lower_sp = val;
      elseif (strcmp (key, 'upper'))
        upper_sp = val;
      endif
      i += 2;
    else
      fit_args{end+1} = a;
      i += 1;
    endif
  endwhile

  if (isempty (penter))
    pe_tbl = struct ('sse', 0.05, 'aic', 0, 'bic', 0, 'rsquared', 0.1, 'adjrsquared', 0);
    penter = pe_tbl.(criterion);
  endif
  if (isempty (premove))
    pr_tbl = struct ('sse', 0.10, 'aic', 0.01, 'bic', 0.01, 'rsquared', 0.05, 'adjrsquared', -0.05);
    premove = pr_tbl.(criterion);
  endif

  fit_args = sw_ensure_modelspec (fit_args);
  sw_predictorvars_matrix_conflict (fit_args);
  [positional, first_nv] = sw_split_positional (fit_args);
  init_spec = positional{end};
  nv_part   = fit_args(first_nv:end);
  resp_name = sw_response_name (positional, nv_part);

  probe_positional = positional;
  probe_positional{end} = 'constant';
  probe_nv   = sw_set_responsevar (nv_part, resp_name);
  probe_args = [probe_positional, probe_nv];

  mdl0 = fitlm (probe_args{:});
  if (! isempty (mdl0.Robust))
    error ("stepwiselm: Robust fitting cannot be combined with stepwise regression.");
  endif
  info = LinearModel.sw_extract (mdl0);

  cat_mask = false (1, info.p_raw);
  if (! isempty (info.cat_info.names))
    for j = 1:info.p_raw
      cat_mask(j) = any (strcmp (info.cat_info.names, info.pred_names{j}));
    endfor
  endif

  T_init_raw  = sw_resolve_bound (init_spec, info.pred_names, cat_mask, info);
  T_lower_raw = sw_resolve_bound (lower_sp,  info.pred_names, cat_mask, info);
  T_upper_raw = sw_resolve_bound (upper_sp,  info.pred_names, cat_mask, info);

  terms_cur   = sw_raw_terms_to_encoded (T_init_raw,  info.pred_names, info.cat_info, info.enc_names);
  T_lower_enc = sw_raw_terms_to_encoded (T_lower_raw, info.pred_names, info.cat_info, info.enc_names);
  T_upper_enc = sw_raw_terms_to_encoded (T_upper_raw, info.pred_names, info.cat_info, info.enc_names);

  groups_upper = sw_group_encoded_terms (T_upper_enc, info.pred_names, info.cat_info, info.enc_names);
  groups_lower = sw_group_encoded_terms (T_lower_enc, info.pred_names, info.cat_info, info.enc_names);

  cur_hi  = any (all (terms_cur(:,1:end-1) == 0, 2));
  cur_fit = LinearModel.lm_fit (build_design (terms_cur, info.X_enc), info.y, info.w, false);

  step_no  = 0;
  any_step = false;

  while (step_no < nsteps)

    groups_cur = sw_group_encoded_terms (terms_cur, info.pred_names, info.cat_info, info.enc_names);
    cur_raws   = sw_group_raws (groups_cur);

    add_cands = sw_add_candidates (groups_upper, cur_raws);

    best_add       = [];
    best_add_score = [];
    for c = 1:numel (add_cands)
      g           = add_cands{c};
      trial_terms = [terms_cur; g.term_rows];
      trial_hi    = any (all (trial_terms(:,1:end-1) == 0, 2));
      trial_fit   = LinearModel.lm_fit (build_design (trial_terms, info.X_enc), info.y, info.w, false);
      score       = sw_criterion_score (criterion, info.n_obs, cur_hi, trial_hi, cur_fit, trial_fit);
      if (verbose >= 2)
        sw_print_candidate (criterion, 'add', g.name, score);
      endif
      if (isempty (best_add) || sw_prefer (criterion, 'add', score, best_add_score))
        best_add       = g;
        best_add.fit   = trial_fit;
        best_add.hi    = trial_hi;
        best_add_score = score;
      endif
    endfor

    if (isempty (add_cands) && verbose >= 2)
      printf ("No candidate terms to add\n");
    endif

    if (! isempty (best_add) && sw_threshold_ok (criterion, 'add', best_add_score, penter))
      terms_cur = [terms_cur; best_add.term_rows];
      cur_fit   = best_add.fit;
      cur_hi    = best_add.hi;
      step_no  += 1;
      any_step  = true;
      if (verbose >= 1)
        sw_print_action (criterion, step_no, 'Adding', best_add.name, best_add_score);
      endif
      continue;
    endif

    groups_cur = sw_group_encoded_terms (terms_cur, info.pred_names, info.cat_info, info.enc_names);
    cur_raws   = sw_group_raws (groups_cur);
    lower_raws = sw_group_raws (groups_lower);

    rem_cands = sw_remove_candidates (groups_cur, lower_raws, cur_raws);

    redundant_g = [];
    for c = 1:numel (rem_cands)
      g    = rem_cands{c};
      keep = true (rows (terms_cur), 1);
      keep (g.idx) = false;
      trial_terms = terms_cur (keep, :);
      trial_fit   = LinearModel.lm_fit (build_design (trial_terms, info.X_enc), info.y, info.w, false);
      if (trial_fit.rank_X == cur_fit.rank_X)
        redundant_g     = g;
        redundant_g.fit = trial_fit;
        redundant_g.hi  = any (all (trial_terms(:,1:end-1) == 0, 2));
        break;
      endif
    endfor

    if (! isempty (redundant_g))
      keep = true (rows (terms_cur), 1);
      keep (redundant_g.idx) = false;
      terms_cur = terms_cur (keep, :);
      cur_fit   = redundant_g.fit;
      cur_hi    = redundant_g.hi;
      step_no  += 1;
      any_step  = true;
      if (verbose >= 1)
        printf ("%d. Removing %s, FStat = Inf, pValue = NaN\n", step_no, redundant_g.name);
      endif
      continue;
    endif

    best_rem       = [];
    best_rem_score = [];
    for c = 1:numel (rem_cands)
      g    = rem_cands{c};
      keep = true (rows (terms_cur), 1);
      keep (g.idx) = false;
      trial_terms = terms_cur (keep, :);
      trial_hi    = any (all (trial_terms(:,1:end-1) == 0, 2));
      trial_fit   = LinearModel.lm_fit (build_design (trial_terms, info.X_enc), info.y, info.w, false);
      score       = sw_criterion_score (criterion, info.n_obs, trial_hi, cur_hi, trial_fit, cur_fit);
      if (verbose >= 2)
        sw_print_candidate (criterion, 'remove', g.name, score);
      endif
      if (isempty (best_rem) || sw_prefer (criterion, 'remove', score, best_rem_score))
        best_rem       = g;
        best_rem.fit   = trial_fit;
        best_rem.hi    = trial_hi;
        best_rem_score = score;
      endif
    endfor

    if (isempty (rem_cands) && verbose >= 2)
      printf ("   No candidate terms to remove\n");
    endif

    if (! isempty (best_rem) && sw_threshold_ok (criterion, 'remove', best_rem_score, premove))
      keep = true (rows (terms_cur), 1);
      keep (best_rem.idx) = false;
      terms_cur = terms_cur (keep, :);
      cur_fit   = best_rem.fit;
      cur_hi    = best_rem.hi;
      step_no  += 1;
      any_step  = true;
      if (verbose >= 1)
        sw_print_action (criterion, step_no, 'Removing', best_rem.name, best_rem_score);
      endif
      continue;
    endif

    break;

  endwhile

  if (! any_step && verbose >= 1)
    printf ("No terms to add to or remove from initial model.\n");
  endif

  int_mask  = all (terms_cur(:, 1:end-1) == 0, 2);
  body      = terms_cur(! int_mask, :);
  n_nonzero = sum (body(:, 1:end-1) != 0, 2);
  degree    = sum (body(:, 1:end-1), 2);
  tier      = zeros (rows (body), 1);
  tier(n_nonzero == 1 & degree == 1) = 1;
  tier(n_nonzero == 2)               = 2;
  tier(n_nonzero == 1 & degree > 1)  = 3;
  bitmask = zeros (rows (body), 1);
  for i = 1:rows (body)
    bitmask(i) = sum (2 .^ (find (body(i, 1:end-1)) - 1));
  endfor
  [~, order] = sortrows ([tier, bitmask]);
  terms_cur  = [terms_cur(int_mask, :); body(order, :)];

  groups_final = sw_group_encoded_terms (terms_cur, info.pred_names, info.cat_info, info.enc_names);
  raws_final   = sw_group_raws (groups_final);
  if (isempty (raws_final))
    used_mask = false (1, numel (info.pred_names));
  else
    used_mask = any (raws_final != 0, 1);
  endif
  used_pred_names = info.pred_names(used_mask);
  used_cols = any (terms_cur(:, 1:end-1) != 0, 1);
  terms_cur = terms_cur(:, [used_cols, true]);

  nv_list = {'Intercept', cur_hi, 'PredictorVars', used_pred_names};
  if (! isempty (info.orig_opts.Weights))
    nv_list = [nv_list, {'Weights', info.orig_opts.Weights}];
  endif
  if (! isempty (info.orig_opts.Exclude))
    nv_list = [nv_list, {'Exclude', info.orig_opts.Exclude}];
  endif
  if (! isempty (info.cat_info.names))
    nv_list = [nv_list, {'CategoricalVars', info.cat_info.names}];
  endif

  mdl = fitlm (info.variables, info.response_name, terms_cur, nv_list{:});

endfunction

function [positional, first_nv] = sw_split_positional (fit_args)
  nv_keys_fitlm = {'varnames','intercept','responsevar','predictorvars', ...
                    'categoricalvars','exclude','weights','robustopts'};
  is_nv = @(s) (ischar (s) || isstring (s)) && any (strcmpi (char (s), nv_keys_fitlm));
  n = numel (fit_args);
  first_nv = n + 1;
  for i = 1:n
    if (is_nv (fit_args{i}))
      first_nv = i;
      break;
    endif
  endfor
  positional = fit_args(1:first_nv-1);
endfunction

function fit_args = sw_ensure_modelspec (fit_args)
  [positional, first_nv] = sw_split_positional (fit_args);
  if (isempty (positional))
    return;
  endif
  if (istable (positional{1}))
    tbl       = positional{1};
    col_names = tbl.Properties.VariableNames;
    if (numel (positional) >= 2)
      arg2 = positional{2};
      if ((ischar (arg2) || isstring (arg2)) && ! any (char (arg2) == '~') ...
          && any (strcmp (char (arg2), col_names)))
        fit_args = [{tbl}, fit_args(3:first_nv-1), ...
                     {'ResponseVar', char(arg2)}, fit_args(first_nv:end)];
        [positional, first_nv] = sw_split_positional (fit_args);
      elseif (isnumeric (arg2) || islogical (arg2))
        [nr2, nc2] = size (arg2);
        if (nc2 != width (tbl) && nc2 == 1 && nr2 == height (tbl))
          ## arg2 is an external y vector, not a terms matrix
          n_cols   = width (tbl);
          col_data = cell (1, n_cols);
          for k = 1:n_cols
            col_data{k} = tbl.(col_names{k});
          endfor
          tbl_ext  = table (col_data{:}, double (arg2(:)), ...
                             'VariableNames', [col_names, {'y'}]);
          fit_args = [{tbl_ext}, fit_args(3:first_nv-1), ...
                       {'ResponseVar', 'y'}, fit_args(first_nv:end)];
          [positional, first_nv] = sw_split_positional (fit_args);
        endif
      endif
    endif
    min_pos = 1;
  else
    min_pos = 2;
  endif
  if (numel (positional) <= min_pos)
    fit_args = [fit_args(1:first_nv-1), {'constant'}, fit_args(first_nv:end)];
  endif
endfunction

function resp_name = sw_response_name (positional, nv_part)
  resp_name = '';
  for i = 1:2:numel (nv_part)-1
    if (strcmpi (char (nv_part{i}), 'ResponseVar'))
      resp_name = char (nv_part{i+1});
    endif
  endfor
  if (! isempty (resp_name))
    return;
  endif
  spec = positional{end};
  if ((ischar (spec) || isstring (spec)) && any (char (spec) == '~'))
    tparts    = strsplit (char (spec), '~');
    resp_name = strtrim (tparts{1});
    return;
  endif
  if (istable (positional{1}))
    resp_name = positional{1}.Properties.VariableNames{end};
  endif
endfunction

function nv_out = sw_set_responsevar (nv_part, resp_name)
  nv_out = nv_part;
  if (isempty (resp_name))
    return;
  endif
  found = false;
  for i = 1:2:numel (nv_out)-1
    if (strcmpi (char (nv_out{i}), 'ResponseVar'))
      nv_out{i+1} = resp_name;
      found = true;
    endif
  endfor
  if (! found)
    nv_out = [nv_out, {'ResponseVar', resp_name}];
  endif
endfunction

function T_raw = sw_reduce_to_raw_padded (mdl_bound, pred_names_full)
  groups_bound = sw_group_encoded_terms (mdl_bound.TermsMatrix, ...
    mdl_bound.PredictorNames, mdl_bound.CatLevelInfo, mdl_bound.EncPredictorNames);
  raws_bound = sw_group_raws (groups_bound);

  n_full = numel (pred_names_full);
  T_raw  = zeros (rows (raws_bound), n_full + 1);
  for r = 1:rows (raws_bound)
    for k = 1:numel (mdl_bound.PredictorNames)
      e = raws_bound (r, k);
      if (e == 0)
        continue;
      endif
      j = find (strcmp (pred_names_full, mdl_bound.PredictorNames{k}), 1);
      T_raw (r, j) = e;
    endfor
  endfor
endfunction

function sw_predictorvars_matrix_conflict (fit_args)
  [positional, first_nv] = sw_split_positional (fit_args);
  n = numel (fit_args);
  pv_value = {};
  rv_value = '';
  for i = first_nv:2:n-1
    if ((ischar (fit_args{i}) || isstring (fit_args{i})) ...
        && strcmpi (char (fit_args{i}), 'PredictorVars'))
      pv_value = fit_args{i+1};
    endif
    if ((ischar (fit_args{i}) || isstring (fit_args{i})) ...
        && strcmpi (char (fit_args{i}), 'ResponseVar'))
      rv_value = char (fit_args{i+1});
    endif
  endfor
  if (isempty (pv_value) || isempty (positional))
    return;
  endif
  if (istable (positional{1}))
    min_pos = 1;
  else
    min_pos = 2;
  endif
  if (numel (positional) <= min_pos)
    return;
  endif
  last_extra = positional{end};
  if (isnumeric (last_extra) || islogical (last_extra))
    error ("stepwiselm: You may not specify PredictorVars with a terms matrix.");
  elseif ((ischar (last_extra) || isstring (last_extra)) ...
          && any (char (last_extra) == '~'))
    result = parseWilkinsonFormula (char (last_extra), 'expand');
    conflict_msg = strcat ("stepwiselm: 'PredictorVars' or 'ResponseVar'", ...
                           " values conflict with the formula character", ...
                           " vector or string scalar.");
    formula_resp = '';
    if (! isempty (result.response) && ! isempty (result.response{1}))
      formula_resp = result.response{1}{1};
    endif
    if (! isempty (rv_value) && ! isempty (formula_resp) ...
        && ! strcmp (rv_value, formula_resp))
      error (conflict_msg);
    endif
    vars_in_formula = {};
    for t = 1:numel (result.model)
      vars_in_formula = [vars_in_formula, result.model{t}];
    endfor
    vars_in_formula = unique (vars_in_formula);
    if (! all (ismember (vars_in_formula, pv_value)))
      error (conflict_msg);
    endif
  endif
endfunction

function T_raw = sw_resolve_bound_formula (raw_spec, info)
  nv_list = {};
  if (! isempty (info.orig_opts.Weights))
    nv_list = [nv_list, {'Weights', info.orig_opts.Weights}];
  endif
  if (! isempty (info.orig_opts.Exclude))
    nv_list = [nv_list, {'Exclude', info.orig_opts.Exclude}];
  endif
  if (! isempty (info.cat_info.names))
    nv_list = [nv_list, {'CategoricalVars', info.cat_info.names}];
  endif
  mdl_bound = fitlm (info.variables, info.response_name, raw_spec, nv_list{:});
  T_raw     = sw_reduce_to_raw_padded (mdl_bound, info.pred_names);
endfunction

function T_raw = sw_resolve_bound (spec, pred_names, cat_mask, info)
  p = numel (pred_names);
  if (isempty (spec))
    spec = 'constant';
  endif
  if (isnumeric (spec) || islogical (spec))
    T = double (spec);
    if (columns (T) == p)
      T = [T, zeros(rows (T), 1)];
    elseif (columns (T) != p + 1)
      error ("stepwiselm: Lower/Upper terms matrix must have %d or %d columns.", p, p+1);
    endif
    T_raw = T;
    return;
  endif
  if (! (ischar (spec) || isstring (spec)))
    error ("stepwiselm: Lower/Upper must be a model name, terms matrix, or polynomial specification.");
  endif
  raw_spec = char (spec);
  if (any (raw_spec == '~'))
    T_raw = sw_resolve_bound_formula (raw_spec, info);
    return;
  endif
  s = lower (strtrim (raw_spec));
  poly_tok = regexp (s, '^poly([0-9]+)$', 'tokens');
  if (! isempty (poly_tok))
    dstr = poly_tok{1}{1};
    if (numel (dstr) != p)
      error (strcat ("stepwiselm: the number of digits in a 'polyijk'", ...
        " specification must match the number of predictors."));
    endif
    digits = dstr - '0';
    T_raw  = sw_poly_terms (digits, p, cat_mask);
    return;
  endif
  switch (s)
    case 'constant'
      T_raw = zeros (1, p + 1);
    case 'linear'
      T_raw = [zeros(1,p+1); [eye(p), zeros(p,1)]];
    case {'interactions', 'purequadratic', 'quadratic', 'full'}
      T_raw = sw_keyword_terms (s, p, cat_mask);
    otherwise
      error ("stepwiselm: '%s' is not a valid model specification.", s);
  endswitch
endfunction

function T = sw_keyword_terms (keyword, p, cat_mask)
  lin = [eye(p), zeros(p,1)];
  T   = [zeros(1,p+1); lin];
  if (strcmp (keyword, 'linear'))
    return;
  endif
  if (any (strcmp (keyword, {'interactions','quadratic','full'})))
    for i = 1:p
      for j = i+1:p
        row = zeros (1, p+1);
        row(i) = 1; row(j) = 1;
        T = [T; row];
      endfor
    endfor
  endif
  if (any (strcmp (keyword, {'purequadratic','quadratic'})))
    for j = 1:p
      if (! cat_mask(j))
        row    = zeros (1, p+1);
        row(j) = 2;
        T      = [T; row];
      endif
    endfor
  endif
  if (strcmp (keyword, 'full'))
    T = zeros (1, p+1);
    for k = 1:p
      idx_mat = nchoosek (1:p, k);
      for r = 1:rows (idx_mat)
        row = zeros (1, p+1);
        row(idx_mat(r,:)) = 1;
        T   = [T; row];
      endfor
    endfor
  endif
endfunction

function T = sw_poly_terms (digits, p, cat_mask)
  ranges = cell (1, p);
  for k = 1:p
    if (cat_mask(k))
      ranges{k} = 0:min (digits(k), 1);
    else
      ranges{k} = 0:digits(k);
    endif
  endfor
  grids = cell (1, p);
  [grids{:}] = ndgrid (ranges{:});
  n = numel (grids{1});
  combos = zeros (n, p);
  for k = 1:p
    gk = grids{k};
    combos(:,k) = gk(:);
  endfor
  keep = false (n, 1);
  for r = 1:n
    e  = combos (r, :);
    nz = find (e != 0);
    if (isempty (nz))
      continue;
    endif
    total_deg = sum (e(nz));
    cap       = max (digits(nz));
    if (total_deg <= cap)
      keep(r) = true;
    endif
  endfor
  T = [zeros(1,p); combos(keep, :)];
  T = [T, zeros(rows (T), 1)];
endfunction

function T_enc = sw_raw_terms_to_encoded (T_raw, pred_names, cat_info, enc_names)
  p_enc = numel (enc_names);
  T_enc = zeros (0, p_enc + 1);
  for r = 1:rows (T_raw)
    G = sw_expand_conceptual_term (T_raw(r, 1:end-1), pred_names, cat_info, enc_names);
    G = [G, zeros(rows (G), 1)];
    T_enc = [T_enc; G];
  endfor
endfunction

function Genc = sw_expand_conceptual_term (term_row, pred_names, cat_info, enc_names)
  p_enc = numel (enc_names);
  Genc  = zeros (1, p_enc);
  for j = find (term_row != 0)
    e  = term_row(j);
    ci = [];
    if (! isempty (cat_info.names))
      ci = find (strcmp (cat_info.names, pred_names{j}));
    endif
    if (isempty (ci))
      col = find (strcmp (enc_names, pred_names{j}), 1);
      Genc (:, col) = e;
    else
      levels_j = cat_info.levels{ci};
      cols_j   = [];
      for L = 2:numel (levels_j)
        nm  = sprintf ("%s_%s", pred_names{j}, char (levels_j{L}));
        col = find (strcmp (enc_names, nm), 1);
        cols_j(end+1) = col;
      endfor
      new_rows = zeros (numel (cols_j) * rows (Genc), p_enc);
      r_out = 0;
      for existing = 1:rows (Genc)
        for c = 1:numel (cols_j)
          r_out += 1;
          row = Genc(existing, :);
          row (cols_j(c)) = 1;
          new_rows (r_out, :) = row;
        endfor
      endfor
      Genc = new_rows;
    endif
  endfor
endfunction

function groups = sw_group_encoded_terms (terms_enc, pred_names, cat_info, enc_names)
  n_pred = numel (pred_names);
  n_rows = rows (terms_enc);
  raw_of = zeros (n_rows, n_pred);
  for r = 1:n_rows
    row_raw = zeros (1, n_pred);
    for j = 1:n_pred
      ci = [];
      if (! isempty (cat_info.names))
        ci = find (strcmp (cat_info.names, pred_names{j}));
      endif
      if (isempty (ci))
        col = find (strcmp (enc_names, pred_names{j}), 1);
        if (! isempty (col))
          row_raw(j) = terms_enc (r, col);
        endif
      else
        levels_j = cat_info.levels{ci};
        present  = false;
        for L = 2:numel (levels_j)
          nm  = sprintf ("%s_%s", pred_names{j}, char (levels_j{L}));
          col = find (strcmp (enc_names, nm), 1);
          if (! isempty (col) && terms_enc (r, col) != 0)
            present = true;
          endif
        endfor
        if (present)
          row_raw(j) = 1;
        endif
      endif
    endfor
    raw_of (r, :) = row_raw;
  endfor

  groups = {};
  used   = false (n_rows, 1);
  for r = 1:n_rows
    if (used(r))
      continue;
    endif
    match = false (n_rows, 1);
    for r2 = r:n_rows
      if (! used(r2) && isequal (raw_of(r2,:), raw_of(r,:)))
        match(r2) = true;
      endif
    endfor
    used = used | match;
    g.raw       = raw_of(r, :);
    g.idx       = find (match)';
    g.term_rows = terms_enc (g.idx, :);
    g.name      = sw_term_name (g.raw, pred_names);
    groups{end+1} = g;
  endfor
endfunction

function nm = sw_term_name (raw_row, pred_names)
  if (all (raw_row == 0))
    nm = '(Intercept)';
    return;
  endif
  parts = {};
  for j = find (raw_row != 0)
    if (raw_row(j) == 1)
      parts{end+1} = pred_names{j};
    else
      parts{end+1} = sprintf ("%s^%d", pred_names{j}, raw_row(j));
    endif
  endfor
  nm = strjoin (parts, ':');
endfunction

function raws = sw_group_raws (groups)
  if (isempty (groups))
    raws = [];
    return;
  endif
  raws = zeros (numel (groups), numel (groups{1}.raw));
  for k = 1:numel (groups)
    raws (k, :) = groups{k}.raw;
  endfor
endfunction

function tf = sw_row_present (row, mat)
  tf = false;
  for k = 1:rows (mat)
    if (isequal (mat(k,:), row))
      tf = true;
      return;
    endif
  endfor
endfunction

function ok = sw_hierarchy_ok_add (term_raw, cur_raws)
  nz = find (term_raw != 0);
  if (numel (nz) == 0)
    ok = true;
    return;
  endif
  if (numel (nz) == 1)
    j = nz(1);
    e = term_raw(j);
    if (e <= 1)
      ok = true;
      return;
    endif
    pred_row    = zeros (size (term_raw));
    pred_row(j) = e - 1;
    ok = sw_row_present (pred_row, cur_raws);
    return;
  endif
  ok = true;
  for k = 1:numel (nz)
    j = nz(k);
    reduced    = term_raw;
    reduced(j) = 0;
    if (! sw_row_present (reduced, cur_raws))
      ok = false;
      return;
    endif
  endfor
endfunction

function ok = sw_hierarchy_ok_remove (term_raw, cur_raws)
  ok = true;
  for k = 1:rows (cur_raws)
    s = cur_raws(k, :);
    if (isequal (s, term_raw))
      continue;
    endif
    if (sw_is_predecessor (term_raw, s))
      ok = false;
      return;
    endif
  endfor
endfunction

function tf = sw_is_predecessor (pred_row, term_raw)
  nz = find (term_raw != 0);
  tf = false;
  if (numel (nz) == 0)
    return;
  endif
  if (numel (nz) == 1)
    j = nz(1);
    e = term_raw(j);
    if (e <= 1)
      return;
    endif
    want    = zeros (size (term_raw));
    want(j) = e - 1;
    tf = isequal (want, pred_row);
    return;
  endif
  for k = 1:numel (nz)
    j = nz(k);
    reduced    = term_raw;
    reduced(j) = 0;
    if (isequal (reduced, pred_row))
      tf = true;
      return;
    endif
  endfor
endfunction

function cand = sw_add_candidates (groups_upper, cur_raws)
  cand = {};
  for u = 1:numel (groups_upper)
    ug = groups_upper{u};
    if (sw_row_present (ug.raw, cur_raws))
      continue;
    endif
    if (sw_hierarchy_ok_add (ug.raw, cur_raws))
      cand{end+1} = ug;
    endif
  endfor
endfunction

function cand = sw_remove_candidates (groups_cur, lower_raws, cur_raws)
  cand = {};
  for c = 1:numel (groups_cur)
    cg = groups_cur{c};
    if (sw_row_present (cg.raw, lower_raws))
      continue;
    endif
    if (sw_hierarchy_ok_remove (cg.raw, cur_raws))
      cand{end+1} = cg;
    endif
  endfor
endfunction

function [Fstat, pval, df1] = sw_ftest (fit_small, fit_big)
  df1 = fit_big.rank_X - fit_small.rank_X;
  if (df1 <= 0 || fit_big.DFE <= 0)
    Fstat = NaN;
    pval  = NaN;
    return;
  endif
  Fstat = ((fit_small.SSE - fit_big.SSE) / df1) / max (fit_big.MSE, eps);
  if (Fstat < 0 || isnan (Fstat))
    pval = NaN;
  else
    pval = betainc (fit_big.DFE / (fit_big.DFE + df1 * Fstat), fit_big.DFE / 2, df1 / 2);
  endif
endfunction

function score = sw_criterion_score (criterion, n_obs, hi_without, hi_with, fit_without, fit_with)
  score = struct ();
  if (strcmp (criterion, 'sse'))
    [Fs, p, ~]   = sw_ftest (fit_without, fit_with);
    score.pvalue = p;
    score.Fstat  = Fs;
    return;
  endif
  crit_without = LinearModel.lm_criteria (fit_without, n_obs, hi_without);
  crit_with    = LinearModel.lm_criteria (fit_with,    n_obs, hi_with);
  switch (criterion)
    case 'aic'
      score.benefit     = crit_with.AIC - crit_without.AIC;
      score.abs_with    = crit_with.AIC;
      score.abs_without = crit_without.AIC;
    case 'bic'
      score.benefit     = crit_with.BIC - crit_without.BIC;
      score.abs_with    = crit_with.BIC;
      score.abs_without = crit_without.BIC;
    case 'rsquared'
      score.benefit     = crit_with.Rsquared - crit_without.Rsquared;
      score.abs_with    = crit_with.Rsquared;
      score.abs_without = crit_without.Rsquared;
    case 'adjrsquared'
      score.benefit     = crit_with.AdjRsquared - crit_without.AdjRsquared;
      score.abs_with    = crit_with.AdjRsquared;
      score.abs_without = crit_without.AdjRsquared;
  endswitch
endfunction

function tf = sw_threshold_ok (criterion, mode, score, thresh)
  if (strcmp (criterion, 'sse'))
    if (strcmp (mode, 'add'))
      tf = ! isnan (score.pvalue) && score.pvalue < thresh;
    else
      tf = ! isnan (score.pvalue) && score.pvalue > thresh;
    endif
    return;
  endif
  if (isnan (score.benefit))
    tf = false;
    return;
  endif
  is_aic_bic = any (strcmp (criterion, {'aic','bic'}));
  if (strcmp (mode, 'add'))
    if (is_aic_bic)
      tf = score.benefit < thresh;
    else
      tf = score.benefit > thresh;
    endif
  else
    if (is_aic_bic)
      tf = score.benefit > thresh;
    else
      tf = score.benefit < thresh;
    endif
  endif
endfunction

function tf = sw_prefer (criterion, mode, score_a, score_b)
  if (strcmp (criterion, 'sse'))
    if (strcmp (mode, 'add'))
      tf = score_a.pvalue < score_b.pvalue;
    else
      tf = score_a.pvalue > score_b.pvalue;
    endif
    return;
  endif
  is_aic_bic = any (strcmp (criterion, {'aic','bic'}));
  if (strcmp (mode, 'add'))
    if (is_aic_bic)
      tf = score_a.benefit < score_b.benefit;
    else
      tf = score_a.benefit > score_b.benefit;
    endif
  else
    if (is_aic_bic)
      tf = score_a.benefit > score_b.benefit;
    else
      tf = score_a.benefit < score_b.benefit;
    endif
  endif
endfunction

function sw_print_candidate (criterion, mode, name, score)
  if (strcmp (mode, 'add'))
    verb = 'adding';
  else
    verb = 'removing';
  endif
  if (strcmp (criterion, 'sse'))
    printf ("   pValue for %s %s is %.6g\n", verb, name, score.pvalue);
    return;
  endif
  labels = struct ('aic','AIC','bic','BIC','rsquared','Rsquared', 'adjrsquared','AdjRsquared');
  lbl = labels.(criterion);
  if (strcmp (mode, 'add'))
    val = score.benefit;
  else
    val = -score.benefit;
  endif
  printf ("   Change in %s for %s %s is %.6g\n", lbl, verb, name, val);
endfunction

function sw_print_action (criterion, step_no, verb, name, score)
  if (strcmp (criterion, 'sse'))
    printf ("%d. %s %s, FStat = %.6g, pValue = %.7g\n", step_no, verb, name, score.Fstat, score.pvalue);
    return;
  endif
  labels = struct ('aic','AIC','bic','BIC','rsquared','Rsquared', 'adjrsquared','AdjRsquared');
  lbl = labels.(criterion);
  if (strcmp (verb, 'Adding'))
    val = score.abs_with;
  else
    val = score.abs_without;
  endif
  printf ("%d. %s %s, %s = %.6g\n", step_no, verb, name, lbl, val);
endfunction

%!demo
%!
%! ## Stepwise search from a constant model, with a custom entry threshold.
%! ## Twenty apartments are described by their size, floor number, and
%! ## distance from the city center, along with their monthly rent.
%! ## Starting from an empty model, `stepwiselm` adds one predictor at a
%! ## time as long as it improves the fit by at least `PEnter`, then checks
%! ## whether anything already in the model should come back out.
%! Size     = [45 50 38 62 70 55 48 40 65 58 42 72 35 68 52 60 46 66 39 54]';
%! Floor    = [3 5 2 8 10 4 1 6 9 7 3 12 2 11 5 6 4 10 1 5]';
%! Distance = [12 8 15 5 3 9 18 10 4 7 14 2 20 3 8 6 11 4 17 9]';
%! Rent = 200 + 12*Size - 15*Distance + 2*Floor + 5*sin ((1:20)'/2);
%! X = [Size, Floor, Distance];
%!
%! ## Fit with a looser entry threshold, printing every candidate examined.
%! mdl = stepwiselm (X, Rent, 'PEnter', 0.06, 'Verbose', 2)

%!demo
%!
%! ## Starting from a formula, with the response, predictors, and
%! ## categorical variable all named explicitly.
%! ## Eighteen coffee shops report their weekly ad spend, staff count, and
%! ## sales season, along with weekly sales. The search starts already
%! ## containing `AdSpend`, with `Employees` and the categorical `Season`
%! ## left to consider adding.
%! AdSpend   = [200 350 500 220 370 520 240 390 540 260 410 560 280 430 580 300 450 600]';
%! Employees = [4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6]';
%! Season    = {'Low';'Mid';'Peak';'Low';'Mid';'Peak';'Low';'Mid';'Peak'; ...
%!              'Low';'Mid';'Peak';'Low';'Mid';'Peak';'Low';'Mid';'Peak'};
%! SeasonEffect = [0;150;400;0;150;400;0;150;400;0;150;400;0;150;400;0;150;400];
%! Sales = 1000 + 0.8*AdSpend + 5*Employees + SeasonEffect + 6*sin ((1:18)');
%! T = table (AdSpend, Employees, Season, Sales, ...
%!   'VariableNames', {'AdSpend','Employees','Season','Sales'});
%!
%! ## Fit, treating Season as categorical and naming the response explicitly.
%! mdl = stepwiselm (T, 'Sales ~ 1 + AdSpend', 'ResponseVar', 'Sales', ...
%!   'PredictorVars', {'AdSpend','Employees','Season'}, ...
%!   'CategoricalVars', {'Season'}, 'Verbose', 1)

%!demo
%!
%! ## Bounding the search with terms matrices instead of model-name
%! ## keywords.
%! ## Twelve potted plants are given varying hours of sunlight and amounts
%! ## of water, and their growth is measured. `T_initial` is a constant
%! ## model and `T_upper` allows the two main effects plus their
%! ## interaction, using the same terms-matrix convention as `fitlm`.
%! Sunlight = [2 4 6 8 3 5 7 9 2.5 4.5 6.5 8.5]';
%! Water    = [100 150 200 250 120 170 220 270 110 160 210 260]';
%! Growth = 5 + 0.3*Sunlight + 0.06*Water + 0.4*sin ((1:12)');
%! X = [Sunlight, Water];
%! T_initial = [0 0 0];
%! T_upper   = [0 0 0; 1 0 0; 0 1 0; 1 1 0];
%!
%! ## Fit, printing the p-value considered for every candidate term.
%! mdl = stepwiselm (X, Growth, T_initial, 'Upper', T_upper, 'Verbose', 2)

%!demo
%!
%! ## A categorical predictor is added or removed as one indicator group,
%! ## never one indicator column at a time.
%! ## Eighteen plots receive varying fertilizer amounts across three
%! ## regions, and crop yield is recorded. The upper bound `poly21` allows
%! ## a constant, `Fertilizer`, `Fertilizer^2`, `Region`, and their
%! ## interaction; `stepwiselm` treats the two indicator columns generated
%! ## by the three-level `Region` as a single term throughout.
%! Fertilizer = repmat ([10;20;30;40;50;60], 3, 1);
%! Region = [repmat({'A'},6,1); repmat({'B'},6,1); repmat({'C'},6,1)];
%! RegionEffect = [zeros(6,1); 15*ones(6,1); 35*ones(6,1)];
%! Yield = 20 + 0.5*Fertilizer + RegionEffect + 1.5*sin ((1:18)');
%! T = table (Fertilizer, Region, Yield);
%!
%! ## Fit, printing every candidate considered at each step.
%! mdl = stepwiselm (T, 'Yield ~ Fertilizer', 'Upper', 'poly21', 'Verbose', 2)

%!demo
%!
%! ## Capping the search with `NSteps` before it converges on its own.
%! ## Twenty observations depend on three predictors, all with genuine
%! ## effects. Left alone, `stepwiselm` would add all three; limiting
%! ## `NSteps` to 1 stops the search after only the single best addition.
%! x1 = (1:20)';
%! x2 = mod ((0:19)', 5);
%! x3 = mod ((0:19)', 3);
%! y  = 5 + 2*x1 + 1.5*x2 + 3*x3 + 2*sin ((1:20)'/1.7);
%! X = [x1, x2, x3];
%!
%! ## Fit, but stop after the very first step.
%! mdl = stepwiselm (X, y, 'Upper', 'linear', 'NSteps', 1, 'Verbose', 1)

%!demo
%!
%! ## Raw matrix input with default `x1, x2, ...` naming, and custom entry
%! ## and removal thresholds.
%! ## Fifteen observations depend mainly on `x1`. `PEnter` and `PRemove`
%! ## are set explicitly rather than relying on the defaults for the SSE
%! ## criterion.
%! x1 = (1:15)';
%! x2 = mod ((0:14)', 4);
%! y  = 8 + 0.5*x1 + 0.3*x2 + 2.5*sin ((1:15)'/1.2);
%! X = [x1, x2];
%!
%! ## Fit with a wider entry threshold and a stricter removal threshold.
%! mdl = stepwiselm (X, y, 'PEnter', 0.2, 'PRemove', 0.3, 'Verbose', 1)

%!demo
%!
%! ## The `hald` cement data, mirroring MATLAB's own reference example for
%! ## `stepwiselm`.
%! ## Four chemical percentages in cement (`ingredients`) are used to
%! ## predict heat given off while hardening (`heat`). Starting from a
%! ## constant model, `stepwiselm` adds three terms and then removes one
%! ## of them again once it becomes redundant.
%! load hald
%!
%! ## Fit with a looser entry threshold than the default.
%! mdl = stepwiselm (ingredients, heat, 'PEnter', 0.06, 'Verbose', 1)

%!demo
%!
%! ## The `carsmall` data, mirroring MATLAB's own reference example for
%! ## `stepwiselm` with a terms-matrix bound.
%! ## Fuel economy (`MPG`) is predicted from acceleration and weight.
%! ## `T_initial` is a constant model and `T_upper` allows both main
%! ## effects plus their interaction.
%! load carsmall
%! X = [Acceleration, Weight];
%! T_initial = [0 0 0];
%! T_upper   = [0 0 0; 1 0 0; 0 1 0; 1 1 0];
%!
%! ## Fit, printing the p-value considered for every candidate term.
%! mdl = stepwiselm (X, MPG, T_initial, 'Upper', T_upper, 'Verbose', 2)

%!shared x1, x2, x3, y, X, tbl
%! n = 24;
%! x1 = (1:n)'/n; x2 = sin((1:n)'/4); x3 = cos((1:n)'/5);
%! y = 4 + 2*x1 - x2 + 0.5*x3 + 0.15*sin((1:n)'*0.8);
%! X = [x1, x2, x3];
%! tbl = table (x1, x2, x3, y, 'VariableNames', {'x1','x2','x3','y'});

%!test
%! mdl = stepwiselm (tbl, 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 24);
%! assert_equal (mdl.NumCoefficients, 4);
%! assert_equal (mdl.NumVariables, 4);
%! assert_equal (mdl.NumPredictors, 3);
%! assert_equal (mdl.DFE, 20);
%! assert_equal (mdl.SSE, 0.237882441508, 1e-9);
%! assert_equal (mdl.SSR, 25.3955253054, 1e-7);
%! assert_equal (mdl.SST, 25.6334077469, 1e-7);
%! assert_equal (mdl.MSE, 0.0118941220754, 1e-10);
%! assert_equal (mdl.RMSE, 0.109060176395, 1e-9);
%! assert_equal (mdl.Rsquared.Ordinary, 0.990719827662, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.989327801811, 1e-9);
%! assert_equal (mdl.LogLikelihood, 21.3138652143, 1e-6);
%! assert_equal (mdl.ModelCriterion.AIC, -34.6277304285, 1e-6);
%! assert_equal (mdl.ModelCriterion.AICc, -32.5224672707, 1e-6);
%! assert_equal (mdl.ModelCriterion.BIC, -29.9155151072, 1e-6);
%! assert_equal (mdl.ModelCriterion.CAIC, -25.9155151072, 1e-6);
%! assert_equal (mdl.ModelFitVsNullModel.Fstat, 711.710796991, 1e-5);
%! assert_equal (mdl.ModelFitVsNullModel.NullModel, 'constant');
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.11199897343; 1.79245385284; -1.06961941012; 0.51029628797], 1e-9);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.0761375845531; 0.144191944451; 0.0575954799856; 0.0489829266327], 1e-9);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [54.007478666; 12.4310262939; -18.5712387566; 10.4178399097], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [3.79220676881e-23; 7.26932473407e-11; 4.42354835071e-14; 1.58256682968e-09], 1e-9);
%! assert_equal (mdl.ResponseName, 'y');
%! assert_equal (mdl.PredictorNames, {'x1','x2','x3'});
%! assert_equal (mdl.VariableNames, {'x1','x2','x3','y'});
%! assert_equal (mdl.Formula.HasIntercept, true);
%! assert_equal (mdl.Formula.LinearPredictor, '1 + x1 + x2 + x3');
%! assert_equal (mdl.Formula.NTerms, 4);

%!test
%! mdl = stepwiselm (tbl, 'y', 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.SSE, 0.237882441508, 1e-9);
%! assert_equal (mdl.ResponseName, 'y');
%! assert_equal (mdl.PredictorNames, {'x1','x2','x3'});

%!test
%! tt = table (x1, x2, x3, 'VariableNames', {'x1','x2','x3'});
%! mdl = stepwiselm (tt, y, 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.SSE, 0.237882441508, 1e-9);
%! assert_equal (mdl.ResponseName, 'y');
%! assert_equal (mdl.VariableNames, {'x1','x2','x3','y'});

%!test
%! mdl = stepwiselm (X, y, 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.11199897343; 1.79245385284; -1.06961941012; 0.51029628797], 1e-9);
%! assert_equal (mdl.ResponseName, 'y');
%! assert_equal (mdl.VariableNames, {'x1','x2','x3','y'});

%!test
%! mdl = stepwiselm (tbl, 'y ~ x1', 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.SSE, 0.237882441508, 1e-9);
%! assert_equal (mdl.Formula.LinearPredictor, '1 + x1 + x2 + x3');

%!test
%! mdl = stepwiselm (X, y, 'constant', 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.SSE, 0.237882441508, 1e-9);

%!test
%! mdl = stepwiselm (X, y, 'linear', 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.SSE, 0.237882441508, 1e-9);

%!test
%! mdl = stepwiselm (X, y, 'constant', 'Upper', 'purequadratic', 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 24);
%! assert_equal (mdl.NumCoefficients, 4);
%! assert_equal (mdl.NumVariables, 4);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 20);
%! assert_equal (mdl.SSE, 0.22345826866, 1e-9);
%! assert_equal (mdl.SSR, 25.4099494782, 1e-7);
%! assert_equal (mdl.SST, 25.6334077469, 1e-7);
%! assert_equal (mdl.MSE, 0.011172913433, 1e-10);
%! assert_equal (mdl.RMSE, 0.105702002975, 1e-9);
%! assert_equal (mdl.Rsquared.Ordinary, 0.991282537583, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.98997491822, 1e-9);
%! assert_equal (mdl.LogLikelihood, 22.0644883645, 1e-6);
%! assert_equal (mdl.ModelCriterion.AIC, -36.1289767289, 1e-6);
%! assert_equal (mdl.ModelCriterion.AICc, -34.023713571, 1e-6);
%! assert_equal (mdl.ModelCriterion.BIC, -31.4167614075, 1e-6);
%! assert_equal (mdl.ModelCriterion.CAIC, -27.4167614075, 1e-6);
%! assert_equal (mdl.ModelFitVsNullModel.Fstat, 758.081874544, 1e-5);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x1^2'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.79002491626; -1.87038904608; -0.881529179216; 3.14370698956], 1e-9);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.0903999262734; 0.328540613357; 0.0534940080998; 0.290849609965], 1e-9);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [52.9870445002; -5.6930223237; -16.4790265402; 10.8087027861], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [5.53960001242e-23; 1.4288675776e-05; 4.19790256549e-13; 8.42439951369e-10], 1e-9);
%! assert_equal (mdl.PredictorNames, {'x1','x2'});
%! assert_equal (mdl.Formula.LinearPredictor, '1 + x1 + x2 + x1^2');

%!test
%! mdl = stepwiselm (X, y, 'constant', 'Upper', 'quadratic', 'Verbose', 0);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x1^2'});
%! assert_equal (mdl.SSE, 0.22345826866, 1e-9);
%! assert_equal (mdl.PredictorNames, {'x1','x2'});

%!test
%! T0 = [0 0 0];
%! T1 = [0 0 0; 1 0 0; 0 1 0; 1 1 0];
%! mdl = stepwiselm (X(:,1:2), y, T0, 'Upper', T1, 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 24);
%! assert_equal (mdl.NumCoefficients, 4);
%! assert_equal (mdl.NumVariables, 3);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 20);
%! assert_equal (mdl.SSE, 1.11546087527, 1e-8);
%! assert_equal (mdl.SSR, 24.5179468716, 1e-7);
%! assert_equal (mdl.SST, 25.6334077469, 1e-7);
%! assert_equal (mdl.MSE, 0.0557730437634, 1e-10);
%! assert_equal (mdl.RMSE, 0.236163171903, 1e-9);
%! assert_equal (mdl.Rsquared.Ordinary, 0.956484097383, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.94995671199, 1e-9);
%! assert_equal (mdl.LogLikelihood, 2.77090924057, 1e-6);
%! assert_equal (mdl.ModelCriterion.AIC, 2.45818151886, 1e-6);
%! assert_equal (mdl.ModelCriterion.AICc, 4.56344467675, 1e-6);
%! assert_equal (mdl.ModelCriterion.BIC, 7.17039684025, 1e-6);
%! assert_equal (mdl.ModelCriterion.CAIC, 11.1703968402, 1e-6);
%! assert_equal (mdl.ModelFitVsNullModel.Fstat, 146.534031599, 1e-5);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x1:x2'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.00100246651; 1.35895695997; -0.23206474908; -1.2780495686], 1e-9);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.181530230599; 0.299744817513; 0.272053305331; 0.469486828642], 1e-9);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [22.0404196772; 4.53371294703; -0.853012055108; -2.72222667524], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [1.67738549624e-15; 0.000202251230338; 0.403753045528; 0.0131231039551], 1e-9);
%! assert_equal (mdl.PredictorNames, {'x1','x2'});
%! assert_equal (mdl.VariableNames, {'x1','x2','y'});

%!test
%! mdl = stepwiselm (X, y, 'Criterion', 'aic', 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.ModelCriterion.AIC, -34.6277304285, 1e-6);

%!test
%! mdl = stepwiselm (X, y, 'Criterion', 'bic', 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.ModelCriterion.BIC, -29.9155151072, 1e-6);

%!test
%! mdl = stepwiselm (X, y, 'Criterion', 'rsquared', 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 2);
%! assert_equal (mdl.NumPredictors, 1);
%! assert_equal (mdl.DFE, 22);
%! assert_equal (mdl.SSE, 2.69619594281, 1e-8);
%! assert_equal (mdl.Rsquared.Ordinary, 0.894817108617, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x2'});
%! assert_equal (mdl.Coefficients.Estimate, [4.93053729606; -1.35113880263], 1e-9);
%! assert_equal (mdl.Coefficients.SE, [0.0714593428648; 0.0987629466731], 1e-9);
%! assert_equal (mdl.Coefficients.tStat, [68.9977978861; -13.6806246486], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, [3.2882259984e-27; 3.08475274476e-12], 1e-9);
%! assert_equal (mdl.PredictorNames, {'x2'});

%!test
%! mdl = stepwiselm (X, y, 'Criterion', 'adjrsquared', 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 5);
%! assert_equal (mdl.DFE, 19);
%! assert_equal (mdl.SSE, 0.224234467489, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.989410626691, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3','x1:x3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.58266504435; 1.22125846816; -1.38808827902; 0.143991325796; 1.04450806346], 1e-8);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.444198853188; 0.55023668447; 0.301652972275; 0.344106973415; 0.971297090441], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [10.3166971537; 2.21951480632; -4.6016065035; 0.418449310593; 1.07537443872], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [3.1776147586e-09; 0.0388206335538; 0.000194754077847; 0.680310511095; 0.295673895373], 1e-8);
%! assert_equal (mdl.PredictorNames, {'x1','x2','x3'});

%!test
%! p1 = (1:24)'/24; p2 = sin ((1:24)'/3.5);
%! pr = 3 + 0.9*p1 - 0.2*p2 + 0.1*sin ((1:24)'*1.1);
%! mdl = stepwiselm (p1, pr, 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 24);
%! assert_equal (mdl.NumCoefficients, 2);
%! assert_equal (mdl.DFE, 22);
%! assert_equal (mdl.SSE, 0.422994127869, 1e-9);
%! assert_equal (mdl.SSR, 2.71741844039, 1e-8);
%! assert_equal (mdl.SST, 3.14041256826, 1e-8);
%! assert_equal (mdl.Rsquared.Ordinary, 0.865306191886, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.859183746062, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1'});
%! assert_equal (mdl.Coefficients.Estimate, [2.85858604186; 1.16664998725], 1e-9);
%! assert_equal (mdl.Coefficients.SE, [0.0584250816201; 0.0981336947312], 1e-9);
%! assert_equal (mdl.Coefficients.tStat, [48.9273778074; 11.8883732081], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, [6.03142199019e-24; 4.75617587395e-11], 1e-9);

%!test
%! p1 = (1:24)'/24; p2 = sin ((1:24)'/3.5);
%! pr = 3 + 0.9*p1 - 0.2*p2 + 0.1*sin ((1:24)'*1.1);
%! mdl = stepwiselm (p1, pr, 'PEnter', 0.5, 'PRemove', 0.6, 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1'});
%! assert_equal (mdl.SSE, 0.422994127869, 1e-9);

%!test
%! mdl = stepwiselm (X, y, 'Upper', 'interactions', 'NSteps', 1, 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 2);
%! assert_equal (mdl.DFE, 22);
%! assert_equal (mdl.SSE, 2.69619594281, 1e-8);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x2'});
%! assert_equal (mdl.Coefficients.Estimate, [4.93053729606; -1.35113880263], 1e-9);
%! assert_equal (mdl.PredictorNames, {'x2'});

%!test
%! mdl = stepwiselm (X, y, 'NSteps', 0, 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 1);
%! assert_equal (mdl.NumPredictors, 0);
%! assert_equal (mdl.DFE, 23);
%! assert_equal (mdl.SSE, 25.6334077469, 1e-7);
%! assert_equal (mdl.SSR, 0, 1e-12);
%! assert_equal (mdl.Rsquared.Ordinary, 0, 1e-12);
%! assert_equal (isnan (mdl.ModelFitVsNullModel.Fstat), true);
%! assert_equal (isnan (mdl.ModelFitVsNullModel.NullModel), true);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)'});
%! assert_equal (mdl.Coefficients.Estimate, 4.92948000444, 1e-9);
%! assert_equal (mdl.Coefficients.SE, 0.215493231622, 1e-9);
%! assert_equal (mdl.Coefficients.tStat, 22.8753356537, 1e-6);
%! assert_equal (mdl.Coefficients.pValue, 2.53817914154e-17, 1e-9);
%! assert_equal (isempty (mdl.PredictorNames), true);
%! assert_equal (mdl.Formula.LinearPredictor, '1');

%!test
%! w = 1400 + (1:45)'*9;
%! g = {'A','B','C'}(mod ((0:44), 3) + 1)';
%! m = 55 - 0.005*w + [0 5 -3](mod ((0:44), 3) + 1)' + 0.3*sin ((1:45)'/6);
%! t = table (w, g, m, 'VariableNames', {'Weight','Group','MPG'});
%! mdl = stepwiselm (t, 'MPG ~ Weight', 'CategoricalVars', {'Group'}, 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 45);
%! assert_equal (mdl.NumCoefficients, 4);
%! assert_equal (mdl.NumVariables, 3);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 41);
%! assert_equal (mdl.SSE, 1.7056787721, 1e-8);
%! assert_equal (mdl.SSR, 512.932602391, 1e-6);
%! assert_equal (mdl.SST, 514.638281163, 1e-6);
%! assert_equal (mdl.MSE, 0.0416019212713, 1e-10);
%! assert_equal (mdl.RMSE, 0.203965490393, 1e-9);
%! assert_equal (mdl.Rsquared.Ordinary, 0.996685674513, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.996443162891, 1e-9);
%! assert_equal (mdl.LogLikelihood, 9.78350141322, 1e-6);
%! assert_equal (mdl.ModelCriterion.AIC, -11.5670028264, 1e-6);
%! assert_equal (mdl.ModelCriterion.AICc, -10.5670028264, 1e-6);
%! assert_equal (mdl.ModelCriterion.BIC, -4.34035286736, 1e-6);
%! assert_equal (mdl.ModelCriterion.CAIC, -0.340352867356, 1e-6);
%! assert_equal (mdl.ModelFitVsNullModel.Fstat, 4109.84706734, 1e-4);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','Weight','Group_B','Group_C'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [56.0078644136; -0.00561620967245; 5.01185752618; -2.97710174861], 1e-8);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.41983060818; 0.000260647333731; 0.074514600823; 0.0746252935319], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [133.40586256; -21.5471594973; 67.2600734732; -39.8940038652], 1e-5);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [1.00643907454e-55; 5.62403610622e-24; 1.37550401987e-43; 1.97956717929e-34], 1e-8);
%! assert_equal (mdl.ResponseName, 'MPG');
%! assert_equal (mdl.PredictorNames, {'Weight','Group'});
%! assert_equal (mdl.Formula.LinearPredictor, '1 + Weight + Group');

%!test
%! w = 1400 + (1:45)'*9;
%! code = mod ((0:44), 3)' + 1;
%! m = 55 - 0.005*w + [0 5 -3](mod ((0:44), 3) + 1)' + 0.3*sin ((1:45)'/6);
%! mdl = stepwiselm ([w, code], m, 'CategoricalVars', 2, 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 41);
%! assert_equal (mdl.SSE, 1.7056787721, 1e-8);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2_2','x2_3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [56.007864413564; -0.00561620967245463; 5.01185752617851; -2.97710174860745], 1e-8);
%! assert_equal (mdl.PredictorNames, {'x1','x2'});

%!test
%! ex = false (24, 1);
%! ex([2 9 15]) = true;
%! mdl = stepwiselm (X, y, 'Upper', 'interactions', 'Exclude', ex, 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 21);
%! assert_equal (mdl.DFE, 17);
%! assert_equal (mdl.SSE, 0.198175777, 1e-9);
%! assert_equal (mdl.SSR, 23.9528827198, 1e-7);
%! assert_equal (mdl.SST, 24.1510584968, 1e-7);
%! assert_equal (mdl.Rsquared.Ordinary, 0.991794323341, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.990346262754, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.10558649487; 1.78751811424; -1.08307909721; 0.499781080674], 1e-9);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.0832423204478; 0.15196995208; 0.0630609308312; 0.0532045604032], 1e-9);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [49.320903992; 11.762312811; -17.1751206798; 9.39357598083], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [8.57561115348e-20; 1.36568780583e-09; 3.54709582656e-12; 3.84225548219e-08], 1e-9);
%! assert_equal (mdl.PredictorNames, {'x1','x2','x3'});

%!test
%! mdl = stepwiselm (X, y, 'Upper', 'interactions', 'Exclude', [2 9 15], 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 21);
%! assert_equal (mdl.SSE, 0.198175777, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});

%!test
%! wt = 0.6 + mod ((1:24)', 4)*0.25;
%! mdl = stepwiselm (X, y, 'Upper', 'interactions', 'Weights', wt, 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 24);
%! assert_equal (mdl.DFE, 20);
%! assert_equal (mdl.SSE, 0.238797389485, 1e-9);
%! assert_equal (mdl.SSR, 24.6325434536, 1e-7);
%! assert_equal (mdl.SST, 24.8713408431, 1e-7);
%! assert_equal (mdl.Rsquared.Ordinary, 0.990398692576, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.988958496462, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.13717568894; 1.74127146188; -1.08433885985; 0.50630139226], 1e-9);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.0789083883313; 0.150888905488; 0.0591696809619; 0.0500332545132], 1e-9);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [52.4301126462; 11.5400894204; -18.3259203401; 10.1192975989], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [6.83329714397e-23; 2.70346109195e-10; 5.69051367149e-14; 2.59082042824e-09], 1e-9);

%!test
%! mdl = stepwiselm (X, y, 'constant', 'Intercept', false, 'Upper', 'linear', 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.SSE, 0.237882441508, 1e-9);
%! assert_equal (mdl.Formula.HasIntercept, true);

%!test
%! mdl = stepwiselm (tbl, 'ResponseVar', 'y', 'PredictorVars', {'x1','x2'}, 'Verbose', 0);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.SSE, 1.11546087527, 1e-8);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x1:x2'});
%! assert_equal (mdl.PredictorNames, {'x1','x2'});

%!test
%! mdl = stepwiselm (tbl, 'ResponseVar', 'y', 'PredictorVars', [1 2], 'Verbose', 0);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.SSE, 1.11546087527, 1e-8);
%! assert_equal (mdl.PredictorNames, {'x1','x2'});

%!test
%! mdl = stepwiselm (tbl, 'ResponseVar', 'y', 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.NumPredictors, 3);
%! assert_equal (mdl.SSE, 0.237882441508, 1e-9);
%! assert_equal (mdl.ResponseName, 'y');

%!test
%! mdl = stepwiselm (X, y, 'VarNames', {'Alpha','Beta','Gamma','Score'}, 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','Alpha','Beta','Gamma'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.11199897343; 1.79245385284; -1.06961941012; 0.51029628797], 1e-9);
%! assert_equal (mdl.ResponseName, 'Score');
%! assert_equal (mdl.PredictorNames, {'Alpha','Beta','Gamma'});
%! assert_equal (mdl.VariableNames, {'Alpha','Beta','Gamma','Score'});
%! assert_equal (mdl.Formula.LinearPredictor, '1 + Alpha + Beta + Gamma');

%!test
%! mdl = stepwiselm (tbl, 'y ~ x1 + x2', 'Lower', 'y ~ x2', 'Upper', 'y ~ x1+x2+x3', 'PRemove', 0.99, 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.11199897343; 1.79245385284; -1.06961941012; 0.51029628797], 1e-9);

%!test
%! xu = (1:50)'/50;
%! idx = mod ((1:50), 2);
%! gu = {'Hi','Lo'}(idx + 1)';
%! sl = [7 1](idx + 1)';
%! yu = 2 + sl.*xu + 0.15*sin ((1:50)'/3);
%! t = table (xu, gu, yu, 'VariableNames', {'X','G','Y'});
%! mdl = stepwiselm (t, 'Y ~ X + G', 'Upper', 'quadratic', 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 50);
%! assert_equal (mdl.NumCoefficients, 4);
%! assert_equal (mdl.NumVariables, 3);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 46);
%! assert_equal (mdl.SSE, 0.540477312442, 1e-9);
%! assert_equal (mdl.SSR, 225.563187153, 1e-6);
%! assert_equal (mdl.SST, 226.103664466, 1e-6);
%! assert_equal (mdl.Rsquared.Ordinary, 0.997609603923, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.997453708527, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','X','G_Hi','X:G_Hi'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [2.02876253697; 0.971349141035; 0.0058812055616; 5.98354193726], 1e-8);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.0433841054684; 0.0751585081172; 0.0622864091272; 0.106290181507], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [46.7628066792; 12.9240077453; 0.0944219717272; 56.2943994677], 1e-5);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [1.96152617542e-40; 6.42479379426e-17; 0.92518406634; 4.45945899483e-44], 1e-8);
%! assert_equal (mdl.PredictorNames, {'X','G'});

%!test
%! idx3 = mod ((1:48), 3); idx2 = mod ((1:48), 2);
%! g1 = {'A','B','C'}(idx3 + 1)';
%! g2 = {'P','Q'}(idx2 + 1)';
%! base = [0 4 -3](idx3 + 1)'; qeff = [5 0](idx2 + 1)';
%! cross = 6*((idx3 == 1) & (idx2 == 0))';
%! y2 = 10 + base + qeff + cross + 0.2*sin ((1:48)'/4);
%! t = table (g1, g2, y2, 'VariableNames', {'G1','G2','Y'});
%! mdl = stepwiselm (t, 'Y ~ G1 + G2', 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 48);
%! assert_equal (mdl.NumCoefficients, 6);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 42);
%! assert_equal (mdl.SSE, 1.00031269326, 1e-8);
%! assert_equal (mdl.SSR, 1526.57402706, 1e-6);
%! assert_equal (mdl.SST, 1527.57433975, 1e-6);
%! assert_equal (mdl.Rsquared.Ordinary, 0.999345162676, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.999267205851, 1e-9);
%! assert_equal (mdl.CoefficientNames, ...
%!   {'(Intercept)','G1_C','G1_A','G2_P','G1_C:G2_P','G1_A:G2_P'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [14.0072303291; -7.00943486401; -4.00436689813; 10.9931096893; -5.98569634506; -6.00058514115], 1e-8);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.0545630013001; 0.0771637364424; 0.0771637364424; 0.0771637364424; 0.1091260026; 0.1091260026], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [256.716639395; -90.8384584155; -51.8944141736; 142.464714594; -54.8512380408; -54.9876747811], 1e-5);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [9.40196273084e-69; 7.64513283056e-50; 1.00685653504e-39; 5.03307771166e-58; 1.01525097736e-40; 9.15939827666e-41], 1e-8);
%! assert_equal (mdl.PredictorNames, {'G1','G2'});

%!test
%! idxb = mod ((1:48), 2);
%! ga = {'M','F'}(idxb + 1)';
%! gb = ga;
%! y3 = 50 - 10*strcmp (ga, 'F') + 0.3*sin ((1:48)'/5);
%! t = table (ga, gb, y3, 'VariableNames', {'GA','GB','Y'});
%! mdl = stepwiselm (t, 'Y ~ 1 + GB', 'ResponseVar', 'Y', ...
%!         'PredictorVars', {'GA','GB'}, 'CategoricalVars', {'GA','GB'}, 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 48);
%! assert_equal (mdl.NumCoefficients, 2);
%! assert_equal (mdl.NumPredictors, 1);
%! assert_equal (mdl.DFE, 46);
%! assert_equal (mdl.SSE, 1.94300386071, 1e-8);
%! assert_equal (mdl.SSR, 1199.4398757, 1e-6);
%! assert_equal (mdl.SST, 1201.38287956, 1e-6);
%! assert_equal (mdl.Rsquared.Ordinary, 0.998382693899, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.998347535071, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','GB_M'});
%! assert_equal (mdl.Coefficients.Estimate, [40.0624369099; 9.99766587633], 1e-8);
%! assert_equal (mdl.Coefficients.SE, [0.041951963781; 0.0593290361473], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, [954.959751562; 168.512191088], 1e-5);
%! assert_equal (mdl.Coefficients.pValue, [1.70539389023e-100; 7.42606982902e-66], 1e-8);
%! assert_equal (mdl.PredictorNames, {'GB'});
%! assert_equal (mdl.Formula.LinearPredictor, '1 + GB');

%!test
%! e1 = (1:30)'/30; e2 = sin ((1:30)'/5); e3 = e1 + e2;
%! ye = 3 + 1.5*e1 - e2 + 0.3*sin ((1:30)'/4);
%! mdl = stepwiselm ([e1,e2,e3], ye, 'linear', 'Upper', 'linear', 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 30);
%! assert_equal (mdl.NumCoefficients, 3);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 27);
%! assert_equal (mdl.SSE, 0.736469147619, 1e-9);
%! assert_equal (mdl.SSR, 30.2842015322, 1e-7);
%! assert_equal (mdl.SST, 31.0206706798, 1e-7);
%! assert_equal (mdl.Rsquared.Ordinary, 0.976258761288, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.974500151013, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x2','x3'});
%! assert_equal (mdl.Coefficients.Estimate, [2.83853512391; -2.5770987124; 1.87079571697], 1e-8);
%! assert_equal (mdl.Coefficients.SE, [0.101176877319; 0.131763778382; 0.186693438852], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, [28.0551762332; -19.5584761158; 10.0206827218], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, [1.65333811724e-21; 1.79264031168e-17; 1.35825364979e-10], 1e-9);
%! assert_equal (mdl.PredictorNames, {'x2','x3'});

%!test
%! h1 = [7;1;11;11;7;11;3;1;2;21;1;11;10];
%! h2 = [26;29;56;31;52;55;71;31;54;47;40;66;68];
%! h3 = [6;15;8;8;6;9;17;22;18;4;23;9;8];
%! h4 = [60;52;20;47;33;22;6;44;22;26;34;12;12];
%! yh = [78.5;74.3;104.3;87.6;95.9;109.2;102.7;72.5;93.1;115.9;83.8;113.3;109.4];
%! mdl = stepwiselm ([h1,h2,h3,h4], yh, 'PEnter', 0.06, 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 13);
%! assert_equal (mdl.NumCoefficients, 3);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 10);
%! assert_equal (mdl.SSE, 57.9044831761, 1e-8);
%! assert_equal (mdl.SSR, 2657.85859375, 1e-5);
%! assert_equal (mdl.SST, 2715.76307692, 1e-5);
%! assert_equal (mdl.Rsquared.Ordinary, 0.978678374536, 1e-9);
%! assert_equal (mdl.Rsquared.Adjusted, 0.974414049443, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2'});
%! assert_equal (mdl.Coefficients.Estimate, [52.5773488821; 1.46830574222; 0.662250491275], 1e-8);
%! assert_equal (mdl.Coefficients.SE, [2.2861743345; 0.121300923606; 0.0458547214685], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, [22.9979613053; 12.1046542645; 14.4423620963], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, [5.45657090149e-10; 2.69221217969e-07; 5.02896031564e-08], 1e-9);
%! assert_equal (mdl.PredictorNames, {'x1','x2'});
%! assert_equal (mdl.VariableNames, {'x1','x2','x3','x4','y'});
%! assert_equal (mdl.Formula.LinearPredictor, '1 + x1 + x2');

%!test
%! w = 1400 + (1:45)'*9;
%! g = {'A','B','C'}(mod ((0:44), 3) + 1)';
%! m = 55 - 0.005*w + [0 5 -3](mod ((0:44), 3) + 1)' + 0.3*sin ((1:45)'/6);
%! t = table (w, g, m, 'VariableNames', {'Weight','Group','MPG'});
%! mdl0 = fitlm (t, 'MPG ~ Weight');
%! mdl = step (mdl0, 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','Weight','Group_B','Group_C'});
%! assert_equal (mdl.SSE, 1.7056787721, 1e-8);
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [56.0078644136; -0.00561620967245; 5.01185752618; -2.97710174861], 1e-8);

%!test
%! w = 1400 + (1:45)'*9;
%! g = {'A','B','C'}(mod ((0:44), 3) + 1)';
%! m = 55 - 0.005*w + [0 5 -3](mod ((0:44), 3) + 1)' + 0.3*sin ((1:45)'/6);
%! t = table (w, g, m, 'VariableNames', {'Weight','Group','MPG'});
%! mdl0 = fitlm (t, 'MPG ~ Weight');
%! mdl = step (mdl0, 'Upper', 'quadratic', 'NSteps', 5, 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 5);
%! assert_equal (mdl.DFE, 40);
%! assert_equal (mdl.SSE, 0.987020309951, 1e-8);
%! assert_equal (mdl.Rsquared.Ordinary, 0.998082108646, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','Weight','Group_B','Group_C','Weight^2'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [82.5932822956; -0.0388795480431; 5.01269583682; -2.97710174861; 1.03495141166e-05], 1e-6);

%!test
%! w = 1400 + (1:45)'*9;
%! code = mod ((0:44), 3)' + 1;
%! m = 55 - 0.005*w + [0 5 -3](mod ((0:44), 3) + 1)' + 0.3*sin ((1:45)'/6);
%! mdl0 = fitlm ([w, code], m);
%! mdl = step (mdl0, 'Criterion', 'bic', 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 2);
%! assert_equal (mdl.DFE, 43);
%! assert_equal (mdl.SSE, 443.573715939, 1e-6);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x2'});
%! assert_equal (mdl.Coefficients.Estimate, [50.739060918498; -1.53909676135581], 1e-6);

%!test
%! tbl_noy = table (x1, x2, x3, 'VariableNames', {'x1','x2','x3'});
%! mdl = stepwiselm (tbl_noy, y, 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 24);
%! assert_equal (mdl.NumCoefficients, 4);
%! assert_equal (mdl.NumPredictors, 3);
%! assert_equal (mdl.DFE, 20);
%! assert_equal (mdl.SSE, 0.237882441508, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.1119989734334; 1.79245385284216; -1.06961941011697; 0.510296287970133], 1e-9);
%! assert_equal (mdl.ResponseName, 'y');
%! assert_equal (mdl.PredictorNames, {'x1','x2','x3'});

%!test
%! T_initial = [0 0 0 0];
%! T_upper   = [0 0 0 0; 1 0 0 0; 0 1 0 0; 1 1 0 0];
%! mdl = stepwiselm (X, y, T_initial, 'Upper', T_upper, 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 4);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 20);
%! assert_equal (mdl.SSE, 1.11546087527, 1e-8);
%! assert_equal (mdl.SSR, 24.5179468716, 1e-6);
%! assert_equal (mdl.Rsquared.Ordinary, 0.956484097383, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x1:x2'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.00100246650923; 1.35895695996517; -0.232064749079705; -1.27804956860349], 1e-8);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.181530230599361; 0.299744817513483; 0.272053305331473; 0.469486828642072], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [22.0404196772024; 4.53371294702715; -0.853012055107928; -2.72222667524045], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [1.67738549623508e-15; 0.000202251230338429; 0.403753045528205; 0.0131231039551409], 1e-8);
%! assert_equal (mdl.PredictorNames, {'x1','x2'});

%!test
%! mdl = stepwiselm (X, y, 'Upper', 'poly110', 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 3);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 21);
%! assert_equal (mdl.SSE, 1.52876802397, 1e-8);
%! assert_equal (mdl.SSR, 24.1046397229, 1e-6);
%! assert_equal (mdl.Rsquared.Ordinary, 0.940360328246, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.21600780425081; 1.37121494151013; -0.897420541182461], 1e-8);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.186735917268724; 0.342414105749106; 0.136495702980697], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [22.577380216489; 4.00455156048637; -6.57471643125185], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [3.27879166781824e-16; 0.000642688076385025; 1.64115221424634e-06], 1e-8);
%! assert_equal (mdl.PredictorNames, {'x1','x2'});

%!test
%! p1   = (1:15)'/15;
%! code = repmat ([1;2;3], 5, 1);
%! resp = [5.18127588719375; 9.35081376514746; 3.49974949866041; 5.6242630760159; ...
%!   9.72651388107706; 3.81411200080599; 5.89825501056437; 9.99098641713587; ...
%!   4.10224698823349; 6.23744090586702; 10.3961126341096; 4.57205845018011; ...
%!   6.75484533214212; 10.9323653265385; 5.09379999767747];
%! mdl = stepwiselm ([p1, code], resp, 'CategoricalVars', logical ([0 1]), ...
%!   'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 15);
%! assert_equal (mdl.NumCoefficients, 4);
%! assert_equal (mdl.NumPredictors, 2);
%! assert_equal (mdl.DFE, 11);
%! assert_equal (mdl.SSE, 0.0670192081399, 1e-8);
%! assert_equal (mdl.Rsquared.Ordinary, 0.999296834962, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2_2','x2_3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [5.04173312790335; 1.92317767382845; 4.01193051752317; -1.97924634508893], 1e-8);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.0482103221773782; 0.0712545629266561; 0.0495946318075163; 0.050272494208676], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [104.577876691085; 26.990238868043; 80.8944510989422; -39.3703629836483], 1e-5);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [7.63821769167572e-18; 2.10244709730975e-11; 1.28303158624727e-16; 3.44025340004292e-13], 1e-8);
%! assert_equal (mdl.PredictorNames, {'x1','x2'});

%!test
%! x1n     = x1; yn = y;
%! x1n(5)  = NaN;
%! yn(12)  = NaN;
%! mdl = stepwiselm ([x1n, x2, x3], yn, 'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.NumObservations, 22);
%! assert_equal (mdl.NumCoefficients, 4);
%! assert_equal (mdl.NumPredictors, 3);
%! assert_equal (mdl.DFE, 18);
%! assert_equal (mdl.SSE, 0.220518784715, 1e-8);
%! assert_equal (mdl.SSR, 23.3502845138, 1e-6);
%! assert_equal (mdl.Rsquared.Ordinary, 0.990644409445, 1e-9);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)','x1','x2','x3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [4.11761577459912; 1.79579438452674; -1.05927400080053; 0.512812488395983], 1e-8);
%! assert_equal (mdl.Coefficients.SE, ...
%!   [0.0775488744294778; 0.146527417121933; 0.0594402526920882; 0.0517068469441139], 1e-8);
%! assert_equal (mdl.Coefficients.tStat, ...
%!   [53.0970411226747; 12.2556885243692; -17.8208192735615; 9.91769018424668], 1e-6);
%! assert_equal (mdl.Coefficients.pValue, ...
%!   [3.09620365929824e-21; 3.59135935529609e-10; 6.98703068441474e-13; 1.01407825601795e-08], 1e-8);
%! assert_equal (mdl.PredictorNames, {'x1','x2','x3'});

%!error <Not enough input arguments> stepwiselm ()
%!error <Name-Value arguments must be in pairs> stepwiselm (X, y, 'Verbose')
%!error <is not a valid value for the 'Criterion' argument> ...
%! stepwiselm (X, y, 'Criterion', 'bogus')
%!error <Robust fitting cannot be combined with stepwise regression> ...
%! stepwiselm (X, y, 'RobustOpts', 'on')
%!error <You may not specify PredictorVars with a terms matrix> ...
%! stepwiselm (X, y, [0 0], 'PredictorVars', {'x1', 'x2'})
%!error <'PredictorVars' or 'ResponseVar' values conflict with the formula> ...
%! stepwiselm (tbl, 'y ~ x1', 'PredictorVars', {'x2'})
%!error <'PredictorVars' or 'ResponseVar' values conflict with the formula> ...
%! stepwiselm (tbl, 'y ~ x1', 'ResponseVar', 'x2', 'PredictorVars', {'x1'})
%!error <Lower/Upper terms matrix must have> stepwiselm (X, y, 'Upper', [0 0 0 0 0])
%!error <Lower/Upper must be a model name, terms matrix, or polynomial specification> ...
%! stepwiselm (X, y, 'Upper', {1, 2})
%!error <the number of digits in a 'polyijk'> stepwiselm (X, y, 'Upper', 'poly1')
%!error <is not a valid model specification> stepwiselm (X, y, 'Upper', 'bogusmodel')
