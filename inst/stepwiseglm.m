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
## @deftypefn  {statistics} {@var{mdl} =} stepwiseglm (@var{X}, @var{y})
## @deftypefnx {statistics} {@var{mdl} =} stepwiseglm (@var{X}, @var{y}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} stepwiseglm (@var{tbl})
## @deftypefnx {statistics} {@var{mdl} =} stepwiseglm (@var{tbl}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} stepwiseglm (@dots{}, @var{Name}, @var{Value})
##
## Fit a generalized linear regression model by stepwise term selection.
##
## @code{@var{mdl} = stepwiseglm (@var{X}, @var{y})} starts from a model given
## by @var{modelspec} and repeatedly adds or removes terms, one at a time, until
## no further move improves the selection criterion.  It returns the fitted
## @code{GeneralizedLinearModel} object @var{mdl}, whose @code{Steps} property
## records the term-selection trace.  @var{X} is an @math{n}-by-@math{p} numeric
## predictor matrix and @var{y} the response; @code{@var{mdl} = stepwiseglm
## (@var{tbl})} instead takes the predictors and response from the table
## @var{tbl} (the last column is the response unless overridden).
##
## @var{modelspec} is the @emph{starting} model.  It is a Wilkinson formula
## string (e.g.@: @qcode{'y ~ x1 + x2'}), a keyword (@qcode{'constant'}
## (default), @qcode{'linear'}, @qcode{'interactions'}, @qcode{'purequadratic'},
## @qcode{'quadratic'}, or @qcode{'full'}), or a terms matrix.  The candidate
## terms available to the search are bounded below by @qcode{'Lower'} and above
## by @qcode{'Upper'}.
##
## The following @var{Name}/@var{Value} pairs control the stepwise search:
##
## @multitable @columnfractions 0.2 0.75
## @headitem Name @tab Value
## @item @qcode{'Lower'} @tab the smallest model considered (terms in it are
## never removed).  Defaults to @qcode{'constant'}.
## @item @qcode{'Upper'} @tab the largest model considered (the candidate term
## universe).  Defaults to @qcode{'interactions'}.
## @item @qcode{'Criterion'} @tab the selection criterion: @qcode{'Deviance'}
## (default), @qcode{'sse'}, @qcode{'aic'}, or @qcode{'bic'}.  Under
## @qcode{'Deviance'} and @qcode{'sse'}, terms enter or leave by a chi-squared
## or @math{F} test on the change in deviance; under @qcode{'aic'}/@qcode{'bic'}
## the move that most reduces the information criterion is taken.
## @item @qcode{'PEnter'} @tab the @math{p}-value (or criterion margin) below
## which a term is added.  Defaults to @math{0.05} for @qcode{'Deviance'} and
## @qcode{'sse'}, and @math{0} for @qcode{'aic'}/@qcode{'bic'}.
## @item @qcode{'PRemove'} @tab the @math{p}-value (or criterion margin) above
## which a term is removed.  Defaults to @math{0.10} for @qcode{'Deviance'} and
## @qcode{'sse'}, and @math{0} for @qcode{'aic'}/@qcode{'bic'}.
## @item @qcode{'NSteps'} @tab the maximum number of steps.  Defaults to
## @code{Inf} (run to convergence).
## @item @qcode{'Verbose'} @tab @math{0} to run silently, or @math{1} (default)
## to print each accepted step.
## @end multitable
##
## All @var{Name}/@var{Value} pairs accepted by @code{fitglm} (such as
## @qcode{'Distribution'}, @qcode{'Link'}, @qcode{'Weights'}, @qcode{'Offset'},
## @qcode{'BinomialSize'}, @qcode{'Intercept'}, @qcode{'DispersionFlag'}, and
## @qcode{'Exclude'}) are also accepted and forwarded to the fit.  Categorical
## predictors are not supported by the stepwise search; use @code{fitglm} for a
## fixed model with categorical predictors.
##
## @seealso{GeneralizedLinearModel, fitglm, stepwisefit, glmfit, glmval}
## @end deftypefn

function mdl = stepwiseglm (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ## ------------------------------------------------------------------------ ##
  ## Separate the data/response, the (optional) starting model, and the pairs.
  ## ------------------------------------------------------------------------ ##
  arg1 = varargin{1};
  if (istable (arg1))
    data = arg1;
    resp = [];
    rest = varargin(2:end);
  else
    if (nargin < 2)
      print_usage ();
    endif
    data = arg1;
    resp = varargin{2};
    rest = varargin(3:end);
  endif

  start_spec = 'constant';
  if (! isempty (rest))
    a = rest{1};
    if ((ischar (a) && ! is_param_name (a)) || isnumeric (a))
      start_spec = a;
      rest       = rest(2:end);
    endif
  endif

  ## Split the pairs into stepwise controls and forwarded fitglm options.
  [sw, opts, glmnv] = parse_pairs (rest);

  ## ------------------------------------------------------------------------ ##
  ## Intake: numeric predictor matrix, response, names, and fitting subset.
  ## ------------------------------------------------------------------------ ##
  start_is_formula = ischar (start_spec) && any (start_spec == '~');
  if (istable (data))
    col_names = data.Properties.VariableNames;
    if (ischar (resp) && ! isempty (resp))
      resp_name = resp;
    elseif (! isempty (opts.ResponseVar))
      resp_name = opts.ResponseVar;
    elseif (start_is_formula)
      resp_name = strtrim (strsplit (start_spec, '~'){1});
    else
      resp_name = col_names{end};
    endif
    if (! isempty (opts.PredictorVars))
      pred_names = opts.PredictorVars;
    else
      pred_names = col_names(! strcmp (col_names, resp_name));
    endif
    p = numel (pred_names);
    X = zeros (height (data), p);
    for j = 1:p
      col = data.(pred_names{j});
      if (iscell (col) || isa (col, 'categorical') || ! isnumeric (col))
        error (strcat ("stepwiseglm: categorical predictors are not", ...
                       " supported; use fitglm for a fixed model."));
      endif
      X(:, j) = double (col(:));
    endfor
    y = double (data.(resp_name)(:));
  else
    if (! (isnumeric (data) && isreal (data) && ismatrix (data)))
      error ("stepwiseglm: X must be a real matrix.");
    endif
    if (! (isnumeric (resp) && isreal (resp) && isvector (resp)))
      error ("stepwiseglm: Y must be a real vector.");
    endif
    X = double (data);
    y = double (resp(:));
    p = columns (X);
    if (rows (X) != numel (y))
      error (strcat ("stepwiseglm: X and Y must have the same number of", ...
                     " observations."));
    endif
    if (! isempty (opts.VarNames))
      pred_names = opts.VarNames(1:p)(:)';
      resp_name  = opts.VarNames{end};
    else
      pred_names = arrayfun (@(k) sprintf ("x%d", k), 1:p, ...
                             'UniformOutput', false);
      resp_name  = 'y';
    endif
    if (! isempty (opts.ResponseVar))
      resp_name = opts.ResponseVar;
    endif
  endif

  if (! isempty (opts.CategoricalVars))
    error (strcat ("stepwiseglm: categorical predictors are not supported;", ...
                   " use fitglm for a fixed model."));
  endif

  ## Distribution and criterion validation.
  known_distr = {'normal', 'binomial', 'poisson', 'gamma', ...
                 'inverse gaussian'};
  distr = opts.Distribution;
  if (! any (strcmp (distr, known_distr)))
    error ("stepwiseglm: unknown distribution '%s'.", distr);
  endif
  crit = lower (sw.Criterion);
  if (! any (strcmp (crit, {'deviance', 'sse', 'aic', 'bic'})))
    error ("stepwiseglm: unknown criterion '%s'.", sw.Criterion);
  endif

  ## Fitting subset (drop missing and excluded rows).
  missing_mask  = any (isnan (X), 2) | isnan (y);
  excluded_mask = false (rows (X), 1);
  if (! isempty (opts.Exclude))
    ex = opts.Exclude(:);
    if (islogical (ex))
      excluded_mask(1:numel (ex)) = ex;
    else
      excluded_mask(ex) = true;
    endif
  endif
  subset = ! missing_mask & ! excluded_mask;
  if (! any (subset))
    error (strcat ("stepwiseglm: no observations remain after removing", ...
                   " missing/excluded rows."));
  endif
  Xs = X(subset, :);
  ys = y(subset);
  n  = rows (Xs);

  ## Weights, offset, and binomial trial counts, subset to the fitting rows.
  w_sub = [];
  if (! isempty (opts.Weights))
    w_sub = double (opts.Weights(:))(subset);
  endif
  off_sub = [];
  if (! isempty (opts.Offset))
    off_sub = double (opts.Offset(:))(subset);
  endif
  N_sub = [];
  if (strcmp (distr, 'binomial') && ! isempty (opts.BinomialSize))
    N = opts.BinomialSize(:);
    if (isscalar (N))
      N = N * ones (rows (X), 1);
    endif
    N_sub = N(subset);
  endif
  w_ll = w_sub;
  if (isempty (w_ll))
    w_ll = ones (n, 1);
  endif

  ## Link function and glmfit argument list (shared by every candidate fit).
  if (! isempty (opts.Link))
    linkspec = opts.Link;
  else
    linkspec = default_link (distr);
  endif
  [~, ~, ilink] = getlinkfunctions (linkspec);
  gargs = {'link', linkspec, 'constant', 'off'};
  if (! isempty (w_sub))
    gargs = [gargs, {'weights', w_sub}];
  endif
  if (! isempty (off_sub))
    gargs = [gargs, {'offset', off_sub}];
  endif
  if (! isempty (opts.DispersionFlag))
    if (opts.DispersionFlag)
      gargs = [gargs, {'estdisp', 'on'}];
    else
      gargs = [gargs, {'estdisp', 'off'}];
    endif
  endif
  yfit = ys;
  if (strcmp (distr, 'binomial') && ! isempty (N_sub))
    yfit = [ys .* N_sub, N_sub];
  endif

  ## Whether dispersion is estimated (chooses the F- vs chi-squared test and
  ## the internal criterion label).
  if (! isempty (opts.DispersionFlag))
    estdisp = logical (opts.DispersionFlag);
  else
    estdisp = any (strcmp (distr, {'normal', 'gamma', 'inverse gaussian'}));
  endif

  ## Package the fitting context passed to the inner fitter.
  ctx = struct ('X', Xs, 'yfit', yfit, 'y', ys, 'distr', distr, ...
                'gargs', {gargs}, 'ilink', ilink, 'off', off_sub, ...
                'N', N_sub, 'w', w_ll, 'n', n);

  ## ------------------------------------------------------------------------ ##
  ## Resolve the starting, lower, and upper models to terms matrices.
  ## ------------------------------------------------------------------------ ##
  T_start = resolve_terms (start_spec, pred_names, p, opts.Intercept);
  T_lower = resolve_terms (sw.Lower, pred_names, p, opts.Intercept);
  T_upper = resolve_terms (sw.Upper, pred_names, p, opts.Intercept);

  ## The search operates within [Lower, Upper]; the start must contain Lower.
  T = union_rows (T_start, T_lower);

  ## Entry/removal thresholds.
  if (any (strcmp (crit, {'deviance', 'sse'})))
    penter  = ternary (isempty (sw.PEnter),  0.05, sw.PEnter);
    premove = ternary (isempty (sw.PRemove), 0.10, sw.PRemove);
  else
    penter  = ternary (isempty (sw.PEnter),  0, sw.PEnter);
    premove = ternary (isempty (sw.PRemove), 0, sw.PRemove);
  endif

  ## ------------------------------------------------------------------------ ##
  ## Stepwise search.
  ## ------------------------------------------------------------------------ ##
  ws = warning ();                # silence per-candidate glmfit warnings
  warning ("off", "all");
  unwind_protect

    f0      = fit_terms (T, ctx);
    hist    = start_history (T, f0.dev, p);
    stepnum = 0;
    do
      stepnum += 1;
      [T, hstep] = one_step (T, T_lower, T_upper, ctx, p, crit, estdisp, ...
                             penter, premove, pred_names);
      if (isempty (hstep))
        break;
      endif
      hist(end+1) = hstep;
      if (sw.Verbose)
        print_step (stepnum, hstep, crit, estdisp, pred_names, p);
      endif
    until (stepnum >= sw.NSteps)

  unwind_protect_cleanup
    warning (ws);
  end_unwind_protect

  ## ------------------------------------------------------------------------ ##
  ## Build the final object (canonical term order) and attach the trace.
  ## ------------------------------------------------------------------------ ##
  T_final = canonical_sort (T, p);
  mdl     = GeneralizedLinearModel (data, resp, T_final, glmnv{:});

  steps = struct ();
  steps.Start     = terms_formula (T_start, resp_name, pred_names, p);
  steps.Lower     = terms_formula (T_lower, resp_name, pred_names, p);
  steps.Upper     = terms_formula (T_upper, resp_name, pred_names, p);
  steps.Criterion = criterion_label (crit, estdisp);
  steps.PEnter    = sw.PEnter;
  steps.PRemove   = sw.PRemove;
  steps.History   = history_table (hist, crit, estdisp, pred_names, p);
  mdl = setSteps (mdl, steps);

endfunction

## ========================================================================== ##
## Argument parsing
## ========================================================================== ##

## True if S names a stepwiseglm or forwarded fitglm parameter.
function tf = is_param_name (s)
  tf = ischar (s) && any (strcmpi (s, {'Lower', 'Upper', 'Criterion', ...
       'PEnter', 'PRemove', 'NSteps', 'Verbose', 'Distribution', 'Link', ...
       'Weights', 'Offset', 'BinomialSize', 'Intercept', 'DispersionFlag', ...
       'CategoricalVars', 'Exclude', 'VarNames', 'PredictorVars', ...
       'ResponseVar'}));
endfunction

## Split Name/Value pairs into stepwise controls (SW), the fit options this
## function needs (OPTS), and the raw pairs forwarded to GeneralizedLinearModel
## (GLMNV).
function [sw, opts, glmnv] = parse_pairs (pairs)
  sw = struct ('Lower', 'constant', 'Upper', 'interactions', ...
               'Criterion', 'Deviance', 'PEnter', [], 'PRemove', [], ...
               'NSteps', Inf, 'Verbose', 1);
  opts = struct ('Distribution', 'normal', 'Link', [], 'Weights', [], ...
                 'Offset', [], 'BinomialSize', [], 'DispersionFlag', [], ...
                 'CategoricalVars', [], 'Exclude', [], 'Intercept', true, ...
                 'VarNames', {{}}, 'PredictorVars', [], 'ResponseVar', []);
  glmnv = {};
  if (mod (numel (pairs), 2) != 0)
    error ("stepwiseglm: Name/Value arguments must come in pairs.");
  endif
  for k = 1:2:numel (pairs)
    name = pairs{k};
    val  = pairs{k+1};
    if (! ischar (name))
      error ("stepwiseglm: parameter names must be character vectors.");
    endif
    switch (lower (name))
      case 'lower'
        sw.Lower = val;
      case 'upper'
        sw.Upper = val;
      case 'criterion'
        sw.Criterion = val;
      case 'penter'
        sw.PEnter = val;
      case 'premove'
        sw.PRemove = val;
      case 'nsteps'
        sw.NSteps = val;
      case 'verbose'
        sw.Verbose = val;
      otherwise
        ## Fit option: record what we need and forward the pair unchanged.
        switch (lower (name))
          case 'distribution';    opts.Distribution   = lower (val);
          case 'link';            opts.Link           = val;
          case 'weights';         opts.Weights        = val;
          case 'offset';          opts.Offset         = val;
          case 'binomialsize';    opts.BinomialSize   = val;
          case 'dispersionflag';  opts.DispersionFlag = val;
          case 'categoricalvars'; opts.CategoricalVars = val;
          case 'exclude';         opts.Exclude        = val;
          case 'intercept';       opts.Intercept      = val;
          case 'varnames';        opts.VarNames       = val;
          case 'predictorvars';   opts.PredictorVars  = val;
          case 'responsevar';     opts.ResponseVar    = val;
          otherwise
            error ("stepwiseglm: unknown parameter name '%s'.", name);
        endswitch
        glmnv = [glmnv, {name, val}];
    endswitch
  endfor
endfunction

## ========================================================================== ##
## Term-selection engine
## ========================================================================== ##

## Perform one step: try to add, then to remove.  Returns the updated terms T
## and a history record HSTEP ([] if no move was accepted).
function [T, hstep] = one_step (T, T_lower, T_upper, ctx, p, crit, estdisp, ...
                                penter, premove, pred_names)
  hstep = [];
  fC    = fit_terms (T, ctx);

  if (any (strcmp (crit, {'deviance', 'sse'})))
    ## --- Add phase: pick the addable term with the smallest p-value. -------
    addable = candidates_add (T, T_upper, p);
    best_p = Inf;  best = [];
    for i = 1:rows (addable)
      r  = addable(i, :);
      fN = fit_terms ([T; r], ctx);
      [stat, pval, ddf] = nested_test (fC.dev, fC.dfe, fN.dev, fN.dfe, estdisp);
      if (pval < best_p)
        best_p = pval;
        best   = struct ('row', r, 'stat', stat, 'pval', pval, ...
                         'ddf', ddf, 'dev', fN.dev, 'df', fN.ncoef);
      endif
    endfor
    if (! isempty (best) && best_p < penter)
      T     = [T; best.row];
      hstep = mkhist ('Add', best.row, T, best.df, best.ddf, best.dev, ...
                      best.stat, best.pval, p);
      return;
    endif

    ## --- Remove phase: pick the removable term with the largest p-value. ---
    removable = candidates_remove (T, T_lower, p);
    worst_p = -Inf;  worst = [];
    for i = 1:rows (removable)
      r  = removable(i, :);
      Tn = remove_row (T, r);
      fN = fit_terms (Tn, ctx);
      [stat, pval, ddf] = nested_test (fN.dev, fN.dfe, fC.dev, fC.dfe, estdisp);
      if (pval > worst_p)
        worst_p = pval;
        worst   = struct ('row', r, 'stat', stat, 'pval', pval, ...
                          'ddf', ddf, 'dev', fN.dev, 'df', fN.ncoef, 'Tn', Tn);
      endif
    endfor
    if (! isempty (worst) && worst_p > premove)
      T     = worst.Tn;
      hstep = mkhist ('Remove', worst.row, T, worst.df, worst.ddf, ...
                      worst.dev, worst.stat, worst.pval, p);
    endif

  else
    ## --- Information criterion: take the single best neighbouring move. ----
    curval  = crit_value (fC, crit, ctx);
    best_val = curval;  best = [];
    addable = candidates_add (T, T_upper, p);
    for i = 1:rows (addable)
      r  = addable(i, :);
      fN = fit_terms ([T; r], ctx);
      v  = crit_value (fN, crit, ctx);
      if (v < best_val - 1e-12)
        best_val = v;
        best = struct ('row', r, 'act', 'Add', 'T', [T; r], ...
                       'dev', fN.dev, 'df', fN.ncoef, 'val', v);
      endif
    endfor
    removable = candidates_remove (T, T_lower, p);
    for i = 1:rows (removable)
      r  = removable(i, :);
      Tn = remove_row (T, r);
      fN = fit_terms (Tn, ctx);
      v  = crit_value (fN, crit, ctx);
      if (v < best_val - 1e-12)
        best_val = v;
        best = struct ('row', r, 'act', 'Remove', 'T', Tn, ...
                       'dev', fN.dev, 'df', fN.ncoef, 'val', v);
      endif
    endfor
    if (! isempty (best))
      T     = best.T;
      hstep = mkhist (best.act, best.row, T, best.df, NaN, best.dev, ...
                      best.val, NaN, p);
    endif
  endif
endfunction

## Fit a candidate model given by terms matrix T.
function f = fit_terms (T, ctx)
  Xd = build_design (T, ctx.X);
  [b, dev, stats] = glmfit (Xd, ctx.yfit, ctx.distr, ctx.gargs{:});
  f = struct ('b', b, 'dev', dev, 'dfe', stats.dfe, 's', stats.s, ...
              'ncoef', rows (T), 'X', Xd);
endfunction

## Nested-model test between the smaller (A) and larger (B) model.  Returns the
## chi-squared or F statistic, its p-value, and the degrees of freedom change.
function [stat, pval, ddf] = nested_test (devA, dfeA, devB, dfeB, estdisp)
  ddf  = dfeA - dfeB;             # number of coefficients that differ (>= 1)
  drop = devA - devB;            # deviance reduction from A to B
  if (estdisp)
    stat = (drop / ddf) / (devB / dfeB);
    pval = 1 - fcdf (stat, ddf, dfeB);
  else
    stat = drop;
    pval = 1 - chi2cdf (stat, ddf);
  endif
endfunction

## Information-criterion value (AIC or BIC) of a fitted candidate.
function v = crit_value (f, crit, ctx)
  eta = f.X * f.b;
  if (! isempty (ctx.off))
    eta = eta + ctx.off;
  endif
  mu = ctx.ilink (eta);
  LL = glm_loglik (ctx.distr, ctx.y, mu, ctx.N, ctx.w, f.s);
  k  = f.ncoef;
  switch (crit)
    case 'aic'
      v = -2 * LL + 2 * k;
    case 'bic'
      v = -2 * LL + k * log (ctx.n);
  endswitch
endfunction

## ========================================================================== ##
## Candidate generation and hierarchy
## ========================================================================== ##

## Rows of the upper model that may be added to T (not present and with every
## one-step-down parent present, so the model hierarchy is preserved).
function C = candidates_add (T, T_upper, p)
  C = zeros (0, columns (T));
  for i = 1:rows (T_upper)
    r = T_upper(i, :);
    if (row_member (r, T))
      continue;
    endif
    if (parents_present (r, T, p))
      C = [C; r];
    endif
  endfor
endfunction

## Rows of T that may be removed (not the intercept, not in the lower model,
## and not a required parent of any higher-order term present in T).
function C = candidates_remove (T, T_lower, p)
  C = zeros (0, columns (T));
  for i = 1:rows (T)
    r = T(i, :);
    if (all (r(1:p) == 0))         # intercept
      continue;
    endif
    if (row_member (r, T_lower))
      continue;
    endif
    if (has_superset (r, T, p))
      continue;
    endif
    C = [C; r];
  endfor
endfunction

## True if every one-step-down parent of term R is present in T.
function tf = parents_present (r, T, p)
  tf = true;
  for j = find (r(1:p) != 0)
    par     = r;
    par(j) -= 1;
    if (any (par(1:p) != 0) && ! row_member (par, T))
      tf = false;
      return;
    endif
  endfor
endfunction

## True if some other row of T is a strict superset of term R.
function tf = has_superset (r, T, p)
  tf = false;
  for i = 1:rows (T)
    s = T(i, :);
    if (isequal (s, r))
      continue;
    endif
    if (all (s(1:p) >= r(1:p)) && any (s(1:p) > r(1:p)))
      tf = true;
      return;
    endif
  endfor
endfunction

function tf = row_member (r, T)
  tf = ! isempty (T) && any (all (T == r, 2));
endfunction

function T = remove_row (T, r)
  T(all (T == r, 2), :) = [];
endfunction

function U = union_rows (A, B)
  U = A;
  for i = 1:rows (B)
    if (! row_member (B(i, :), U))
      U = [U; B(i, :)];
    endif
  endfor
endfunction

## ========================================================================== ##
## Spec resolution and term naming
## ========================================================================== ##

## Resolve a model specification (keyword, terms matrix, or formula) to a terms
## matrix over the P predictors, with the trailing response column.
function T = resolve_terms (spec, pred_names, p, intercept)
  if (ischar (spec) && any (spec == '~'))
    res    = parseWilkinsonFormula (spec, 'expand');
    model  = res.model;
    T      = zeros (0, p + 1);
    has_int = false;
    for i = 1:numel (model)
      if (isempty (model{i}))
        has_int = true;
      else
        row = term_row (model{i}, pred_names, p);
        T   = [T; row];
      endif
    endfor
    if (has_int && intercept)
      T = [zeros(1, p + 1); T];
    endif
  else
    [T, ~, ~, emsg] = parse_modelspec (spec, pred_names, p, intercept);
    if (! isempty (emsg))
      error ("stepwiseglm: %s", emsg);
    endif
  endif
endfunction

## Build a terms-matrix row from a cell of factor strings (e.g. {'x1','x2^2'}).
function row = term_row (factors, pred_names, p)
  row = zeros (1, p + 1);
  for k = 1:numel (factors)
    f   = factors{k};
    pw  = 1;
    hat = strfind (f, '^');
    if (! isempty (hat))
      pw = str2double (f(hat(1)+1:end));
      f  = f(1:hat(1)-1);
    endif
    j = find (strcmp (pred_names, f));
    if (isempty (j))
      error ("stepwiseglm: unknown predictor '%s' in formula.", f);
    endif
    row(j) += pw;
  endfor
endfunction

## Display name of a terms-matrix row (e.g. 'x1:x2', 'x3^2', '(Intercept)').
function nm = term_name (r, pred_names, p)
  if (all (r(1:p) == 0))
    nm = '(Intercept)';
    return;
  endif
  parts = {};
  for j = find (r(1:p) != 0)
    if (r(j) == 1)
      parts{end+1} = pred_names{j};
    else
      parts{end+1} = sprintf ("%s^%d", pred_names{j}, r(j));
    endif
  endfor
  nm = strjoin (parts, ':');
endfunction

## Wilkinson-style formula string for a terms matrix (used for Steps.Start etc).
function s = terms_formula (T, resp_name, pred_names, p)
  has_int = false;
  parts   = {};
  for i = 1:rows (T)
    if (all (T(i, 1:p) == 0))
      has_int = true;
    else
      parts{end+1} = term_name (T(i, :), pred_names, p);
    endif
  endfor
  if (has_int)
    rhs = strjoin ([{'1'}, parts], ' + ');
  elseif (isempty (parts))
    rhs = '1';
  else
    rhs = strjoin (parts, ' + ');
  endif
  s = sprintf ("%s ~ %s", resp_name, rhs);
endfunction

## Sort terms into canonical order: intercept first, then by degree, then by
## the sorted list of predictor indices (matching fitglm's coefficient order).
function T = canonical_sort (T, p)
  m   = rows (T);
  deg = zeros (m, 1);
  idx = zeros (m, p);            # padded index lists
  for i = 1:m
    lst = [];
    for j = find (T(i, 1:p) != 0)
      lst = [lst, repmat(j, 1, T(i, j))];
    endfor
    deg(i)             = numel (lst);
    idx(i, 1:numel (lst)) = lst;
  endfor
  [~, ord] = sortrows ([deg, idx]);
  T = T(ord, :);
endfunction

## ========================================================================== ##
## History
## ========================================================================== ##

function h = start_history (T, dev, p)
  h = struct ('Action', 'Start', 'TermName', '1', 'Terms', T, ...
              'DF', rows (T), 'delDF', NaN, 'Deviance', dev, ...
              'Stat', NaN, 'PValue', NaN);
endfunction

function h = mkhist (action, r, T, df, ddf, dev, stat, pval, p)
  ## pred_names not needed here: the term name is resolved by the caller
  h = struct ('Action', action, 'TermName', r, 'Terms', T, 'DF', df, ...
              'delDF', ddf, 'Deviance', dev, 'Stat', stat, 'PValue', pval);
endfunction

## Assemble the History table, naming the statistic column after the criterion.
function tbl = history_table (hist, crit, estdisp, pred_names, p)
  m        = numel (hist);
  Action   = cell (m, 1);
  TermName = cell (m, 1);
  Terms    = cell (m, 1);
  DF       = zeros (m, 1);
  delDF    = zeros (m, 1);
  Deviance = zeros (m, 1);
  Stat     = zeros (m, 1);
  PValue   = zeros (m, 1);
  for i = 1:m
    Action{i} = hist(i).Action;
    if (ischar (hist(i).TermName))
      TermName{i} = hist(i).TermName;
    else
      TermName{i} = term_name (hist(i).TermName, pred_names, p);
    endif
    Terms{i}    = hist(i).Terms;
    DF(i)       = hist(i).DF;
    delDF(i)    = hist(i).delDF;
    Deviance(i) = hist(i).Deviance;
    Stat(i)     = hist(i).Stat;
    PValue(i)   = hist(i).PValue;
  endfor
  switch (crit)
    case 'deviance'
      statname = ternary (estdisp, 'FStat', 'Chi2Stat');
    case 'sse'
      statname = 'FStat';
    otherwise
      statname = upper (crit);
  endswitch
  tbl = table (Action, TermName, Terms, DF, delDF, Deviance, Stat, PValue, ...
               'VariableNames', {'Action', 'TermName', 'Terms', 'DF', ...
               'delDF', 'Deviance', statname, 'PValue'});
endfunction

## Print an accepted step to the console (Verbose = 1).
function print_step (stepnum, h, crit, estdisp, pred_names, p)
  nm = term_name (h.TermName, pred_names, p);
  if (strcmp (h.Action, 'Add'))
    verb = 'Adding';
  else
    verb = 'Removing';
  endif
  switch (crit)
    case 'deviance'
      if (estdisp)
        printf ("%d. %s %s, Deviance = %g, FStat = %g, PValue = %g\n", ...
                stepnum, verb, nm, h.Deviance, h.Stat, h.PValue);
      else
        printf ("%d. %s %s, Deviance = %g, Chi2Stat = %g, PValue = %g\n", ...
                stepnum, verb, nm, h.Deviance, h.Stat, h.PValue);
      endif
    case 'sse'
      printf ("%d. %s %s, FStat = %g, pValue = %g\n", ...
              stepnum, verb, nm, h.Stat, h.PValue);
    case 'aic'
      printf ("%d. %s %s, AIC = %g\n", stepnum, verb, nm, h.Stat);
    case 'bic'
      printf ("%d. %s %s, BIC = %g\n", stepnum, verb, nm, h.Stat);
  endswitch
endfunction

## ========================================================================== ##
## Small helpers ported from GeneralizedLinearModel
## ========================================================================== ##

function s = criterion_label (crit, estdisp)
  switch (crit)
    case 'deviance'
      s = ternary (estdisp, 'deviance_F', 'deviance_chi2');
    otherwise
      s = crit;
  endswitch
endfunction

function spec = default_link (distr)
  switch (distr)
    case 'normal';            spec = 'identity';
    case 'binomial';          spec = 'logit';
    case 'poisson';           spec = 'log';
    case 'gamma';             spec = 'reciprocal';
    case 'inverse gaussian';  spec = -2;
  endswitch
endfunction

function ll = glm_loglik (distr, y, mu, N, w, phi)
  rmin = realmin;
  switch (distr)
    case 'normal'
      ne = sum (w);
      s2 = sum (w .* (y - mu) .^ 2) / ne;
      ll = -0.5 * ne * (log (2 * pi * s2) + 1);
    case 'poisson'
      ll = sum (w .* (y .* log (max (mu, rmin)) - mu - gammaln (y + 1)));
    case 'binomial'
      if (isempty (N))
        N = ones (size (y));
      endif
      yc = y .* N;
      ll = sum (w .* (gammaln (N + 1) - gammaln (yc + 1) ...
           - gammaln (N - yc + 1) + yc .* log (max (mu, rmin)) ...
           + (N - yc) .* log (max (1 - mu, rmin))));
    case 'gamma'
      a  = 1 ./ phi;
      ll = sum (w .* (a .* log (a) - a .* log (mu) + (a - 1) .* log (y) ...
           - a .* y ./ mu - gammaln (a)));
    case 'inverse gaussian'
      ll = sum (w .* (-0.5 * (log (2 * pi * phi .* y .^ 3) ...
           + (y - mu) .^ 2 ./ (phi .* mu .^ 2 .* y))));
  endswitch
endfunction

function out = ternary (cond, a, b)
  if (cond)
    out = a;
  else
    out = b;
  endif
endfunction

%!demo
%! ## Stepwise Poisson regression: start from a constant model and let the
%! ## search add the predictors that matter.
%! X = [0.83, -0.68; 0.22, 0.93; -0.12, 0.72; 0.55, -2.55; 1.89, 1.39; ...
%!      -1.46, -1.18; 1.06, -0.75; -0.89, 0.85; 0.19, -0.71; -0.43, -0.57];
%! y = [2; 2; 2; 1; 9; 0; 2; 1; 1; 1];
%! mdl = stepwiseglm (X, y, 'constant', 'Distribution', 'poisson', ...
%!                    'Upper', 'linear')

%!demo
%! ## Stepwise logistic regression selected by AIC, reported from a table.
%! X = [0.83, 1.02; 0.22, 0.29; -0.12, 0.09; 0.55, 0.56; 1.89, 1.96; ...
%!      -1.46, -0.46; 1.06, -0.87; -0.89, 1.18; 0.19, -1.00; -0.43, -0.04];
%! y = [1; 1; 0; 1; 1; 0; 1; 0; 1; 0];
%! tbl = array2table ([X, y], 'VariableNames', {'x1', 'x2', 'y'});
%! mdl = stepwiseglm (tbl, 'constant', 'Distribution', 'binomial', ...
%!                    'Criterion', 'aic', 'Verbose', 0)

## Shared test data (identical to the MATLAB stepwiseglm verification probes).
%!shared X, yb, yp, yn
%! X = [ 0.83, -0.68,  1.02,  0.02;  0.22,  0.93,  0.29, -0.57; ...
%!      -0.12,  0.72,  0.09, -1.03;  0.55, -2.55,  0.56,  0.97; ...
%!       1.89,  1.39,  1.96,  0.16; -1.46, -1.18, -0.46,  0.33; ...
%!       1.06, -0.75, -0.87,  0.08; -0.89,  0.85,  1.18, -1.04; ...
%!       0.19, -0.71, -1.00,  1.85; -0.43, -0.57, -0.04,  0.75; ...
%!      -0.90, -1.75, -0.61, -0.68;  1.52, -0.10,  0.43, -1.16; ...
%!       0.58, -1.63,  0.08, -0.73;  0.11,  0.87, -0.40, -0.15; ...
%!       1.26, -0.42, -0.95, -1.07; -0.02,  1.09,  1.03,  1.54; ...
%!       0.80,  0.97,  0.53,  0.62; -0.40, -1.18,  2.83, -1.36; ...
%!      -0.61, -0.44, -1.85, -1.57;  1.22, -0.94, -0.26, -0.75; ...
%!      -0.84,  0.11, -0.55, -0.42;  1.66, -1.27, -1.05, -0.98; ...
%!       0.29,  0.59, -0.36, -0.74; -1.10, -1.29,  0.40,  0.35; ...
%!       0.08, -0.90, -0.74, -0.51;  1.48, -0.04, -0.11, -1.47; ...
%!       0.06,  0.05, -0.09, -1.55;  0.84,  0.13, -1.11, -1.24; ...
%!      -0.72, -0.02,  0.75, -0.39;  1.11,  0.34, -1.32, -0.18; ...
%!       1.48, -1.55, -0.06,  1.28; -0.16, -0.24,  1.34, -1.92; ...
%!      -0.67,  0.80,  1.01, -0.07; -0.26, -1.10, -0.33, -0.88; ...
%!      -0.26, -2.42,  0.77, -1.44; -0.62, -0.38, -0.38, -1.77; ...
%!       0.41,  0.90,  0.85,  0.72; -1.16,  0.99,  0.74, -2.00; ...
%!      -0.44, -1.22,  1.06, -1.14;  1.85,  0.21,  0.65, -0.11];
%! yb = [1 1 0 1 1 0 1 0 1 0 0 1 0 1 1 1 1 1 0 1 0 0 0 0 0 1 0 1 0 1 ...
%!       0 1 0 0 1 0 1 0 0 1]';
%! yp = [2 2 2 1 9 0 2 1 1 1 0 3 1 2 2 2 3 1 1 2 1 2 2 0 1 3 1 2 1 3 ...
%!       2 1 1 1 0 1 3 1 0 5]';
%! yn = [3.17 1.66 0.99 3.87 6.22 1.34 4.8 0.58 4.05 2.06 2.81 3.97 3.84 ...
%!       1.11 4.26 2.57 3.07 -1.07 1.7 4.58 0.47 6.41 1.55 1.02 3.29 3.81 ...
%!       1.58 2.65 0.81 2.94 5.93 1.09 1.06 2.55 1.44 0.97 2.8 -0.47 0.81 ...
%!       4.7]';

## Test results (values verified against MATLAB's stepwiseglm)
%!test
%! ## Binomial: the Deviance criterion selects x1 and x3.
%! mdl = stepwiseglm (X, yb, 'constant', 'Distribution', 'binomial', ...
%!                    'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)', 'x1', 'x3'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [-0.567847008659897; 2.24005553968806; 1.10257647012], 1e-8);
%! assert_equal (mdl.Deviance, 35.424118776089266, 1e-9);
%! assert_equal (mdl.LogLikelihood, -17.712059388044633, 1e-9);

%!test
%! ## Poisson: selects x1 and x2 (interaction not significant).
%! mdl = stepwiseglm (X, yp, 'constant', 'Distribution', 'poisson', ...
%!                    'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)', 'x1', 'x2'});
%! assert_equal (mdl.Coefficients.Estimate, ...
%!   [0.242293607853268; 0.672082935896644; 0.459103523205052], 1e-8);
%! assert_equal (mdl.Deviance, 5.62441316952502, 1e-9);

%!test
%! ## Normal (estimated dispersion, F-test): mains plus x1:x4 and x2:x3.
%! mdl = stepwiseglm (X, yn, 'constant', 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, ...
%!   {'(Intercept)', 'x1', 'x2', 'x3', 'x4', 'x1:x4', 'x2:x3'});
%! assert_equal (mdl.NumCoefficients, 7);
%! assert_equal (any (strcmp (mdl.CoefficientNames, 'x2:x3')), true);

%!test
%! ## The AIC criterion is greedier and also brings in x3:x4.
%! mdl = stepwiseglm (X, yn, 'linear', 'Upper', 'interactions', ...
%!                    'Criterion', 'aic', 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 8);
%! assert_equal (any (strcmp (mdl.CoefficientNames, 'x3:x4')), true);

%!test
%! ## BIC reaches the same eight-coefficient model here.
%! mdl = stepwiseglm (X, yn, 'constant', 'Criterion', 'bic', 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 8);
%! assert_equal (any (strcmp (mdl.CoefficientNames, 'x3:x4')), true);

%!test
%! ## A formula-valued Upper bounds the candidate universe.
%! mdl = stepwiseglm (X, yp, 'y ~ x1', 'Distribution', 'poisson', ...
%!                    'Upper', 'y ~ x1 + x2 + x3', 'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)', 'x1', 'x2'});

%!test
%! ## Lower/Upper keywords: search within [linear, quadratic] adds interactions.
%! mdl = stepwiseglm (X, yn, 'linear', 'Lower', 'linear', ...
%!                    'Upper', 'quadratic', 'Verbose', 0);
%! assert_equal (mdl.NumCoefficients, 7);
%! assert_equal (any (strcmp (mdl.CoefficientNames, 'x2:x3')), true);
%! assert_equal (any (strcmp (mdl.CoefficientNames, 'x1:x4')), true);

%!test
%! ## The Steps property records the trace; History matches the deviance test.
%! mdl = stepwiseglm (X, yb, 'constant', 'Distribution', 'binomial', ...
%!                    'Upper', 'interactions', 'Verbose', 0);
%! assert_equal (mdl.Steps.Criterion, 'deviance_chi2');
%! assert_equal (size (mdl.Steps.History, 1), 3);
%! assert_equal (mdl.Steps.History.Action{1}, 'Start');
%! assert_equal (mdl.Steps.History.Chi2Stat(2), 14.8398950335268, 1e-6);

%!test
%! ## Table input, response taken from the last column by default.
%! tbl = array2table ([X(:,1:2), yp], 'VariableNames', {'a', 'b', 'y'});
%! mdl = stepwiseglm (tbl, 'constant', 'Distribution', 'poisson', ...
%!                    'Verbose', 0);
%! assert_equal (mdl.CoefficientNames, {'(Intercept)', 'a', 'b'});

## Test input validation
%!error<Invalid call> stepwiseglm ()
%!error<stepwiseglm: X must be a real matrix.> stepwiseglm ("a", [1;2])
%!error<stepwiseglm: unknown distribution 'wibble'.> ...
%! stepwiseglm ([1, 2; 3, 4], [1; 0], 'Distribution', 'wibble')
%!error<stepwiseglm: unknown criterion 'nope'.> ...
%! stepwiseglm ([1, 2; 3, 4], [1; 0], 'Criterion', 'nope')
%!error<stepwiseglm: unknown parameter name 'foo'.> ...
%! stepwiseglm ([1, 2; 3, 4], [1; 0], 'linear', 'foo', 1)
%!error<categorical predictors are not supported> ...
%! stepwiseglm ([1, 2; 3, 4], [1; 0], 'CategoricalVars', 1)
