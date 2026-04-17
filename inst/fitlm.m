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
## @deftypefn  {statistics} {@var{mdl} =} fitlm (@var{X}, @var{y})
## @deftypefnx {statistics} {@var{mdl} =} fitlm (@var{X}, @var{y}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} fitlm (@var{X}, @var{y}, @var{modelspec}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {@var{mdl} =} fitlm (@var{tbl})
## @deftypefnx {statistics} {@var{mdl} =} fitlm (@var{tbl}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} fitlm (@var{tbl}, @var{modelspec}, @var{Name}, @var{Value})
##
## Fit a linear regression model.
##
## @code{fitlm} fits a linear regression model and returns a
## @code{LinearModel} object.
##
## @strong{Inputs}
##
## @var{X} is an @math{n x p} numeric matrix of predictor values.
## @var{y} is an @math{n x 1} numeric response vector.
## @var{tbl} is a @code{table} object containing all variables.
##
## @var{modelspec} can be:
## @itemize
## @item @qcode{'linear'} (default) — intercept + all linear terms
## @item @qcode{'constant'} — intercept only
## @item @qcode{'interactions'} — linear + all pairwise interactions
## @item @qcode{'purequadratic'} — linear + all squared terms
## @item @qcode{'quadratic'} — linear + interactions + squared terms
## @item A Wilkinson formula string, e.g. @qcode{'y ~ x1 + x2 + x1:x2'}
## @end itemize
##
## @strong{Name-Value Arguments}
##
## @table @asis
## @item @qcode{'Weights'} — n×1 double, per-observation weights
## @item @qcode{'Exclude'} — integer indices or logical vector of rows to exclude
## @item @qcode{'Intercept'} — logical (default: @code{true})
## @item @qcode{'CategoricalVars'} — integer indices or cell of names
## @item @qcode{'VarNames'} — cell of char, override variable names (matrix input)
## @item @qcode{'ResponseVar'} — char, override response variable (table input)
## @item @qcode{'PredictorVars'} — cell or indices, specify predictor subset (table input)
## @end table
##
## @strong{Output}
##
## @var{mdl} is a @code{LinearModel} object with all properties populated.
##
## @seealso{LinearModel, CompactLinearModel, stepwiselm}
## @end deftypefn

function mdl = fitlm (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ## ── [Step 1] Parse inputs ────────────────────────────────────────────────
  [raw_X, raw_y, var_names, cat_flags, weights, excl_mask, ...
   modelspec, obs_names, resp_name, pred_names, has_intercept, raw_cols, ...
   dummy_coding] = ...
    __fitlm_parse_inputs__ (varargin{:});

  n = rows (raw_y);
  p = numel (pred_names);

  ## ── [Step 2] Build design matrix ─────────────────────────────────────────
  [terms_mat, coef_names, X_full, y_full, incl_mask, p_tot, has_intercept, lp_str] = ...
    __fitlm_build_design__ (raw_X, raw_y, modelspec, var_names, cat_flags, ...
                            weights, excl_mask, pred_names, raw_cols, ...
                            has_intercept, dummy_coding);

  ## ── [Step 3] Call math engine on included rows only ──────────────────────
  X_fit = X_full(incl_mask, :);
  y_fit = y_full(incl_mask);
  w_fit = weights(incl_mask);

  warning ("off", "CompactLinearModel:rankDeficient");
  s = lm_fit_engine (X_fit, y_fit, w_fit);
  warning ("on",  "CompactLinearModel:rankDeficient");

  n_used = s.n;

  ## ── [Step 4] Build Formula struct ────────────────────────────────────────
  ## fromFit expects formula.CoefficientNames to be set
  Formula = __fitlm_build_formula__ (terms_mat, var_names, pred_names, ...
                                     resp_name, has_intercept, lp_str, coef_names);
  Formula.CoefficientNames = coef_names;   ## required by fromFit (line 288)

  ## ── [Step 5] Build VariableInfo struct ────────────────────────────────────
  VariableInfo = __fitlm_build_varinfo__ (var_names, cat_flags, raw_cols, ...
                                          Formula.InModel);

  ## ── [Step 6] Build ObservationInfo ───────────────────────────────────────
  ObservationInfo = __fitlm_build_obsinfo__ (weights, excl_mask, raw_X, raw_y, obs_names);

  ## ── [Step 7] Build Variables struct ──────────────────────────────────────
  Variables = __fitlm_build_variables_table__ (raw_cols, var_names, obs_names);

  ## ── [Step 8] Compute missing mask ────────────────────────────────────────
  missing_mask = ObservationInfo.Missing;

  ## ── [Step 9] Call LinearModel.fromFit ────────────────────────────────────
  ##
  ## fromFit signature:
  ##   fromFit(fit_s, X, y, formula, pred_names, resp_name,
  ##           var_names, var_info, obs_names, obs_info,
  ##           weights, excluded, missing_mask, steps, variables_data)
  ##
  ## X and y here are the USED (included) subsets — fromFit computes diagnostics
  mdl = LinearModel.fromFit ( ...
    s, ...                   ## lm_fit_engine output
    X_fit, ...               ## design matrix (used rows only)
    y_fit, ...               ## response (used rows only)
    Formula, ...             ## Formula struct (with CoefficientNames field)
    pred_names, ...          ## predictor variable names
    resp_name, ...           ## response variable name
    var_names, ...           ## all variable names
    VariableInfo, ...        ## VariableInfo struct
    obs_names, ...           ## ObservationNames
    ObservationInfo, ...     ## ObservationInfo struct
    weights, ...             ## full n weights
    excl_mask, ...           ## full n exclusion mask
    missing_mask, ...        ## full n missing mask
    [], ...                  ## Steps = [] for fitlm
    Variables, ...           ## Variables struct
    X_full);                 ## full n design matrix for excluded-row Fitted

endfunction


## ── Private helper — build Variables struct ───────────────────────────────────

function Vars = __fitlm_build_variables_table__ (raw_cols, var_names, obs_names)
  ## Build Variables as a struct where each field = one column (all n rows).
  Vars = struct ();
  for k = 1:numel (var_names)
    vn = var_names{k};
    ## Replace invalid fieldname chars
    fn = regexprep (vn, '[^a-zA-Z0-9_]', '_');
    if (isempty (fn) || ! isempty (regexp (fn(1), '[^a-zA-Z_]')))
      fn = ['v_', fn];
    endif
    Vars.(fn) = raw_cols{k};
  endfor
endfunction


## ─── UNIT TESTS ───────────────────────────────────────────────────────────────

## Test 1: basic matrix input returns LinearModel
%!test
%! rng (42);
%! X = randn (20, 3);
%! y = X * [1; 2; -1] + randn (20, 1) * 0.5;
%! mdl = fitlm (X, y);
%! assert (isa (mdl, 'LinearModel'));
%! assert (mdl.NumPredictors, 3);
%! assert (mdl.NumObservations, 20);
%! assert (mdl.NumCoefficients, 4);
%! assert (isempty (mdl.Steps));
%! assert (mdl.DFE, 16);

## Test 2: VariableNames and ResponseName for matrix input
%!test
%! rng (1);
%! X = randn (10, 2);
%! y = X * [2; -1] + randn (10, 1) * 0.3;
%! mdl = fitlm (X, y);
%! assert (isequal (mdl.VariableNames, {'x1', 'x2', 'y'}));
%! assert (strcmp (mdl.ResponseName, 'y'));
%! assert (isequal (mdl.PredictorNames, {'x1', 'x2'}));
%! assert (mdl.NumVariables, 3);

## Test 3: VarNames override
%!test
%! rng (1);
%! X = randn (15, 2);
%! y = X * [3; -1] + randn (15, 1) * 0.5;
%! mdl = fitlm (X, y, 'VarNames', {'Height', 'Weight', 'BP'});
%! assert (strcmp (mdl.ResponseName, 'BP'));
%! assert (isequal (mdl.PredictorNames, {'Height', 'Weight'}));
%! assert (isequal (mdl.CoefficientNames, {'(Intercept)', 'Height', 'Weight'}));

## Test 4: Formula struct fields
%!test
%! rng (1);
%! X = randn (10, 2);
%! y = randn (10, 1);
%! mdl = fitlm (X, y);
%! assert (isfield (mdl.Formula, 'LinearPredictor'));
%! assert (isfield (mdl.Formula, 'ResponseName'));
%! assert (isfield (mdl.Formula, 'Terms'));
%! assert (isfield (mdl.Formula, 'HasIntercept'));
%! assert (isfield (mdl.Formula, 'InModel'));
%! assert (isfield (mdl.Formula, 'TermNames'));
%! assert (isfield (mdl.Formula, 'PredictorNames'));
%! assert (strcmp (mdl.Formula.LinearPredictor, '1 + x1 + x2'));
%! assert (mdl.Formula.HasIntercept, true);
%! assert (isequal (mdl.Formula.Terms, [0 0 0; 1 0 0; 0 1 0]));

## Test 5: modelspec 'constant'
%!test
%! rng (1);
%! X = randn (10, 3);
%! y = randn (10, 1);
%! mdl = fitlm (X, y, 'constant');
%! assert (mdl.NumCoefficients, 1);
%! assert (isequal (mdl.CoefficientNames, {'(Intercept)'}));

## Test 6: modelspec 'interactions'
%!test
%! rng (1);
%! X = randn (15, 2);
%! y = randn (15, 1);
%! mdl = fitlm (X, y, 'interactions');
%! assert (mdl.NumCoefficients, 4);   ## 1 + x1 + x2 + x1:x2
%! assert (any (strcmp (mdl.CoefficientNames, 'x1:x2')));

## Test 7: Intercept=false
%!test
%! rng (42);
%! X = randn (20, 2);
%! y = X * [3; -1] + randn (20, 1) * 0.3;
%! mdl = fitlm (X, y, 'Intercept', false);
%! assert (mdl.NumCoefficients, 2);
%! assert (mdl.DFE, 18);
%! assert (mdl.Formula.HasIntercept, false);
%! assert (! any (strcmp (mdl.CoefficientNames, '(Intercept)')));

## Test 8: Exclude integer indices
%!test
%! rng (1);
%! X = randn (10, 2);
%! y = X * [2; -1] + randn (10, 1) * 0.3;
%! mdl = fitlm (X, y, 'Exclude', [3, 7]);
%! assert (mdl.NumObservations, 8);
%! assert (mdl.ObservationInfo.Excluded(3), true);
%! assert (mdl.ObservationInfo.Excluded(7), true);
%! assert (mdl.ObservationInfo.Subset(3),   false);
%! assert (isnan (mdl.Residuals.Raw(3)));
%! assert (isnan (mdl.Residuals.Raw(7)));

## Test 9: Exclude logical vector gives same coefficients as integer
%!test
%! rng (1);
%! X = randn (10, 2);
%! y = X * [2; -1] + randn (10, 1) * 0.3;
%! mdl_int = fitlm (X, y, 'Exclude', [3, 7]);
%! excl_log = false (10, 1); excl_log([3,7]) = true;
%! mdl_log = fitlm (X, y, 'Exclude', excl_log);
%! assert (mdl_int.Coefficients.Estimate, mdl_log.Coefficients.Estimate, 1e-12);

## Test 10: NaN in X → Missing flag, Fitted = NaN
%!test
%! rng (1);
%! X = randn (10, 2);
%! y = X * [2; -1] + randn (10, 1) * 0.3;
%! X(5, 1) = NaN;
%! mdl = fitlm (X, y);
%! assert (mdl.NumObservations, 9);
%! assert (mdl.ObservationInfo.Missing(5), true);
%! assert (isnan (mdl.Fitted(5)));
%! assert (isnan (mdl.Residuals.Raw(5)));

## Test 11: Weighted regression
%!test
%! rng (42);
%! X = randn (15, 2);
%! y = X * [3; -1] + randn (15, 1) * 0.5;
%! w = abs (randn (15, 1)) + 0.5;
%! mdl_w = fitlm (X, y, 'Weights', w);
%! mdl_u = fitlm (X, y);
%! assert (mdl_w.DFE, mdl_u.DFE);
%! assert (mdl_w.NumObservations, 15);
%! assert (mdl_w.ObservationInfo.Weights, w, 1e-14);

## Test 12: CategoricalVars on numeric matrix
%!test
%! rng (42);
%! X = [randn(30, 1), repmat([1;2;3], 10, 1), randn(30, 1)];
%! y = X(:,1)*2 + (X(:,2)==2)*1.5 + randn(30,1)*0.3;
%! mdl = fitlm (X, y, 'CategoricalVars', 2);
%! assert (mdl.VariableInfo.IsCategorical(2), true);
%! assert (mdl.NumCoefficients, 5);   ## intercept + x1 + x2_2 + x2_3 + x3
%! assert (mdl.NumPredictors, 3);
%! assert (any (strcmp (mdl.CoefficientNames, 'x2_2')));
%! assert (any (strcmp (mdl.CoefficientNames, 'x2_3')));

## Test 13: compact() preserves key fields
%!test
%! rng (1);
%! X = randn (15, 2);
%! y = X * [2; -1] + randn (15, 1) * 0.3;
%! mdl = fitlm (X, y);
%! cmdl = compact (mdl);
%! assert (isa (cmdl, 'CompactLinearModel'));
%! assert (! isa (cmdl, 'LinearModel'));
%! assert (cmdl.MSE,  mdl.MSE,  1e-10);
%! assert (cmdl.RMSE, mdl.RMSE, 1e-10);
%! assert (cmdl.DFE,  mdl.DFE);

## Test 14: rank-deficient — NumEstimatedCoefficients < NumCoefficients
%!test
%! warning ('off', 'CompactLinearModel:rankDeficient');
%! rng (1);
%! X_base = randn (20, 2);
%! X_rd = [X_base, X_base(:,1) + X_base(:,2)];
%! y = X_base * [2; -1] + randn (20, 1) * 0.3;
%! mdl = fitlm (X_rd, y);
%! assert (mdl.NumCoefficients, 4);
%! assert (mdl.NumEstimatedCoefficients, 3);
%! warning ('on', 'CompactLinearModel:rankDeficient');

## Test 15: Steps is [] for fitlm models
%!test
%! rng (1);
%! X = randn (10, 2);
%! y = randn (10, 1);
%! mdl = fitlm (X, y);
%! assert (isempty (mdl.Steps));

## Test 16: ObservationInfo fields present
%!test
%! rng (1);
%! X = randn (8, 2);
%! y = randn (8, 1);
%! mdl = fitlm (X, y, 'Exclude', 3);
%! oi = mdl.ObservationInfo;
%! assert (isfield (oi, 'Weights'));
%! assert (isfield (oi, 'Excluded'));
%! assert (isfield (oi, 'Missing'));
%! assert (isfield (oi, 'Subset'));

## Test 17: Wilkinson formula string input
%!test
%! rng (1);
%! X = randn (20, 2);
%! y = X * [1; 2] + randn (20, 1) * 0.5;
%! mdl = fitlm (X, y, 'y ~ x1 + x2');
%! assert (isa (mdl, 'LinearModel'));
%! assert (mdl.Formula.HasIntercept, true);
%! assert (strcmp (mdl.Formula.ResponseName, 'y'));

## Test 18: anova() still works on fitlm output
%!test
%! rng (1);
%! X = randn (20, 2);
%! y = X * [1; 2] + randn (20, 1) * 0.5;
%! mdl = fitlm (X, y);
%! T = anova (mdl);
%! assert (isstruct (T));
