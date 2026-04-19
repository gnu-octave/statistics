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

classdef LinearModel < CompactLinearModel
## -*- texinfo -*-
## @deftp {statistics} LinearModel
##
## Linear regression model.
##
## A @code{LinearModel} object encapsulates the results of fitting a linear
## regression model.  It inherits from @code{CompactLinearModel} and adds
## the source training data, observation-level residuals and diagnostics,
## and model modification methods.
##
## A @code{LinearModel} object is created by @code{fitlm} or
## @code{stepwiselm}.  Users do not call the constructor directly.
##
## @strong{Additional Properties (beyond CompactLinearModel)}
##
## @multitable @columnfractions 0.25 0.75
## @headitem Property @tab Description
## @item @code{Variables} @tab Original input data (all rows, including excluded).
## @item @code{ObservationNames} @tab Cell array of row names.
## @item @code{ObservationInfo} @tab Table with columns Weights, Excluded, Missing, Subset.
## @item @code{Fitted} @tab Fitted values for all observations.
## @item @code{Residuals} @tab Table with columns Raw, Pearson, Studentized, Standardized.
## @item @code{Diagnostics} @tab Table with columns Leverage, CooksDistance, Dffits, S2_i, CovRatio, Dfbetas, HatMatrix.
## @item @code{ModelFitVsNullModel} @tab Struct with fields Fstat, Pvalue, NullModel.
## @item @code{Steps} @tab Stepwise fitting history (empty for fitlm).
## @end multitable
##
## @strong{Methods}
##
## @multitable @columnfractions 0.25 0.75
## @headitem Method @tab Description
## @item @code{compact} @tab Strip to CompactLinearModel.
## @item @code{addTerms} @tab Add terms and refit.
## @item @code{removeTerms} @tab Remove terms and refit.
## @item @code{step} @tab Single optimal stepwise step.
## @item @code{anova} @tab ANOVA table (partial SS).
## @end multitable
##
## @seealso{fitlm, CompactLinearModel}
## @end deftp

  properties (SetAccess = protected)
    ## Source data (not on CompactLinearModel)
    Variables = [];
    ObservationNames = {};
    ObservationInfo = [];

    ## Diagnostics (not on CompactLinearModel)
    Fitted = [];
    Residuals = [];
    Diagnostics = [];
    ModelFitVsNullModel = [];

    ## Fitting method
    Steps = [];
  endproperties


  methods (Access = protected)

    function obj = LinearModel (s)
      ## LinearModel  Constructor (internal use only).
      ##
      ## obj = LinearModel (s)
      ##
      ## S is a struct with two groups of fields:
      ##   - All CompactLinearModel fields (passed to parent constructor)
      ##   - LinearModel-only fields: Variables, ObservationNames,
      ##     ObservationInfo, Fitted, Residuals, Diagnostics,
      ##     ModelFitVsNullModel, Steps

      if (nargin < 1)
        ## Allow zero-arg construction for Octave classdef internals
        return;
  endif

      if (!isstruct(s))
          error("LinearModel: constructor argument must be a struct.");
  endif

      ## Call parent constructor with the CompactLinearModel fields
      obj = obj@CompactLinearModel (s);

      ## Assign LinearModel-only properties
      if (isfield (s, "Variables"))
        obj.Variables = s.Variables;
      endif
      if (isfield (s, "ObservationNames"))
        obj.ObservationNames = s.ObservationNames;
      endif
      if (isfield (s, "ObservationInfo"))
        obj.ObservationInfo = s.ObservationInfo;
      endif
      if (isfield (s, "Fitted"))
        obj.Fitted = s.Fitted;
      endif
      if (isfield (s, "Residuals"))
        obj.Residuals = s.Residuals;
      endif
      if (isfield (s, "Diagnostics"))
        obj.Diagnostics = s.Diagnostics;
      endif
      if (isfield (s, "ModelFitVsNullModel"))
        obj.ModelFitVsNullModel = s.ModelFitVsNullModel;
      endif
      if (isfield (s, "Steps"))
        obj.Steps = s.Steps;
      endif

    endfunction

  endmethods

  methods (Static, Access = public)

    function obj = fromFit (fit_s, X, y, formula, pred_names, resp_name, ...
                              var_names, var_info, obs_names, obs_info, ...
                              weights, excluded, missing_mask, steps, ...
                              variables_data, X_full)
      ## fromFit  Build a LinearModel from lm_fit_engine output.
      ##
      ## This static factory method takes the raw lm_fit_engine struct and
      ## all metadata, computes observation-level diagnostics, and returns
      ## a fully-populated LinearModel object.
      ##
      ## Parameters:
      ##   fit_s          - struct from lm_fit_engine(X_used, y_used, w_used)
      ##   X              - design matrix for USED (subset) observations only
      ##   y              - response vector for USED observations only
      ##   formula         - Formula struct
      ##   pred_names      - cell of predictor names
      ##   resp_name       - response variable name (char)
      ##   var_names       - cell of all variable names
      ##   var_info        - variable info struct/table
      ##   obs_names       - cell of observation names (or {})
      ##   obs_info        - struct with Weights, Excluded, Missing, Subset
      ##   weights         - original weights vector (n_total x 1, or [])
      ##   excluded        - logical n_total x 1 excluded mask
      ##   missing_mask    - logical n_total x 1 missing mask
      ##   steps           - Steps struct ([] for fitlm)
      ##   variables_data  - Variables struct/table (source data)

      n_total = numel (excluded);
      n_used = fit_s.n;
      subset = obs_info.Subset;

      ## Compute Fitted for ALL observations (including excluded/missing)
      fitted_all = NaN (n_total, 1);
      fitted_all(subset) = fit_s.yhat;
      ## Excluded-but-not-missing rows: compute yhat = X_full * beta (Phase 4)
      ## MATLAB Block 7: excluded rows DO get Fitted values (X * beta extrapolation)
      if (nargin >= 16 && ! isempty (X_full))
        excl_not_missing = excluded & ! missing_mask;
        if (any (excl_not_missing))
          fitted_all(excl_not_missing) = X_full(excl_not_missing, :) * fit_s.beta;
        endif
      endif

      ## Compute diagnostics for used observations
      h     = fit_s.h;
      raw   = fit_s.resid;
      rmse  = fit_s.rmse;
      mse   = fit_s.mse;
      sse   = fit_s.sse;
      dfe   = fit_s.dfe;
      p_est = fit_s.p_est;

      ## Residual types (for used observations only)
      pearson      = raw / rmse;
      standardized = raw ./ (rmse * sqrt (1 - h));
      S2_i_used    = (sse - raw .^ 2 ./ (1 - h)) / (dfe - 1);
      studentized  = raw ./ (sqrt (S2_i_used) .* sqrt (1 - h));

      ## Expand residuals to all observations (NaN for excluded/missing)
      raw_all          = NaN (n_total, 1);
      pearson_all      = NaN (n_total, 1);
      studentized_all  = NaN (n_total, 1);
      standardized_all = NaN (n_total, 1);

      raw_all(subset)          = raw;
      pearson_all(subset)      = pearson;
      studentized_all(subset)  = studentized;
      standardized_all(subset) = standardized;

      ## Build Residuals table (MATLAB returns n x 4 table)
      if (! isempty (obs_names))
        resid_tbl = table (raw_all, pearson_all, studentized_all, ...
                           standardized_all, ...
                           "VariableNames", {"Raw", "Pearson", ...
                                             "Studentized", "Standardized"}, ...
                           "RowNames", obs_names(:)');
      else
        resid_tbl = table (raw_all, pearson_all, studentized_all, ...
                           standardized_all, ...
                           "VariableNames", {"Raw", "Pearson", ...
                                             "Studentized", "Standardized"});
      endif

      ## Diagnostics for used observations
      CooksDistance_used = (1 / p_est) * (h ./ (1 - h)) .* standardized .^ 2;
      Dffits_used = studentized .* sqrt (h ./ (1 - h));
      ## CovRatio: Candidate 3 formula, verified to 3.55e-14 against MATLAB
      CovRatio_used = (S2_i_used / mse) .^ p_est ./ (1 - h);

      ## Hat matrix from QR factor (n_used x n_used)
      HatMatrix_full = NaN (n_total, n_total);
      Q_est = fit_s.Q;
      H_used = Q_est * Q_est';
      ## Place into full matrix
      idx_used = find (subset);
      HatMatrix_full(idx_used, idx_used) = H_used;
      ## Excluded rows: set to 0
      idx_excl = find (! subset);
      HatMatrix_full(idx_excl, :) = 0;
      HatMatrix_full(:, idx_excl) = 0;

      ## Leverage for all observations
      leverage_all = zeros (n_total, 1);
      leverage_all(subset) = h;

      ## CooksDistance, Dffits, S2_i, CovRatio for all obs
      CooksDistance_all = NaN (n_total, 1);
      Dffits_all        = NaN (n_total, 1);
      S2_i_all          = NaN (n_total, 1);
      CovRatio_all      = NaN (n_total, 1);

      CooksDistance_all(subset) = CooksDistance_used;
      Dffits_all(subset)        = Dffits_used;
      S2_i_all(subset)          = S2_i_used;
      CovRatio_all(subset)      = CovRatio_used;

      ## Dfbetas (n_total x p_tot) via update formula
      p_tot = fit_s.p_tot;
      Dfbetas_all = NaN (n_total, p_tot);
      vcov = fit_s.vcov;
      se_coef = sqrt (diag (vcov));
      for i = 1:n_used
        ii = idx_used(i);
        x_i = X(i, :)';
        e_i = raw(i);
        h_i = h(i);
        s2_i = S2_i_used(i);
        if (s2_i > 0 && (1 - h_i) > 0)
          Dfbetas_all(ii, :) = ...
            ((vcov * x_i * e_i) / ((1 - h_i) * sqrt (s2_i))) ./ se_coef;
        endif
      endfor

      ## Build Diagnostics table (MATLAB returns n x 7 table)
      if (! isempty (obs_names))
        diag_tbl = table (leverage_all, CooksDistance_all, Dffits_all, ...
                          S2_i_all, CovRatio_all, Dfbetas_all, ...
                          HatMatrix_full, ...
                          "VariableNames", {"Leverage", "CooksDistance", ...
                                            "Dffits", "S2_i", "CovRatio", ...
                                            "Dfbetas", "HatMatrix"}, ...
                          "RowNames", obs_names(:)');
      else
        diag_tbl = table (leverage_all, CooksDistance_all, Dffits_all, ...
                          S2_i_all, CovRatio_all, Dfbetas_all, ...
                          HatMatrix_full, ...
                          "VariableNames", {"Leverage", "CooksDistance", ...
                                            "Dffits", "S2_i", "CovRatio", ...
                                            "Dfbetas", "HatMatrix"});
      endif

      ## ModelFitVsNullModel
      if (p_est > 1 && mse > 0)
        df_reg = p_est - 1;
        F_null = (fit_s.ssr / df_reg) / mse;
        p_null = 1 - fcdf (F_null, df_reg, dfe);
      else
        F_null = NaN;
        p_null = NaN;
      endif
      mfvnm = struct ("Fstat", F_null, "Pvalue", p_null, ...
                       "NullModel", "constant");

      ## Build struct for parent + child constructor
      s = struct ();

      ## CompactLinearModel fields
      s.Coefficients             = table (fit_s.beta, fit_s.se, ...
                                           fit_s.tstat, fit_s.pval, ...
                                           "VariableNames", {"Estimate", "SE", "tStat", "pValue"}, ...
                                           "RowNames", formula.CoefficientNames(:)');
      s.CoefficientNames         = formula.CoefficientNames;
      s.CoefficientCovariance    = fit_s.vcov;
      s.NumCoefficients          = fit_s.p_tot;
      s.NumEstimatedCoefficients = fit_s.p_est;
      s.DFE                      = fit_s.dfe;
      s.SSE                      = fit_s.sse;
      s.SSR                      = fit_s.ssr;
      s.SST                      = fit_s.sst;
      s.MSE                      = fit_s.mse;
      s.RMSE                     = fit_s.rmse;
      s.Rsquared                 = fit_s.rsq;
      s.LogLikelihood            = fit_s.logL;
      s.ModelCriterion           = fit_s.crit;
      s.Robust                   = [];
      s.Formula                  = LinearFormula (formula);
      s.NumObservations          = n_used;
      s.NumPredictors            = numel (pred_names);
      s.NumVariables             = numel (var_names);
      s.PredictorNames           = pred_names;
      s.ResponseName             = resp_name;
      ## Wrap VariableInfo as table if it is still a raw struct
      if (isstruct (var_info) && ! isempty (var_info))
        s.VariableInfo = table (var_info.Class, var_info.Range, ...
                                var_info.InModel, var_info.IsCategorical, ...
                                "VariableNames", {"Class", "Range", ...
                                                  "InModel", "IsCategorical"}, ...
                                "RowNames", var_names(:)');
      else
        s.VariableInfo = var_info;
      endif
      s.VariableNames            = var_names;

      ## LinearModel-only fields
      if (nargin >= 15 && ! isempty (variables_data))
        s.Variables = variables_data;
      else
        s.Variables = [];
      endif
      s.ObservationNames  = obs_names;
      ## Wrap ObservationInfo as table if it is still a raw struct
      if (isstruct (obs_info) && ! isempty (obs_info))
        if (! isempty (obs_names))
          s.ObservationInfo = table (obs_info.Weights, obs_info.Excluded, ...
                                     obs_info.Missing, obs_info.Subset, ...
                                     "VariableNames", {"Weights", "Excluded", ...
                                                       "Missing", "Subset"}, ...
                                     "RowNames", obs_names(:)');
        else
          s.ObservationInfo = table (obs_info.Weights, obs_info.Excluded, ...
                                     obs_info.Missing, obs_info.Subset, ...
                                     "VariableNames", {"Weights", "Excluded", ...
                                                       "Missing", "Subset"});
        endif
      else
        s.ObservationInfo = obs_info;
      endif
      s.Fitted            = fitted_all;
      s.Residuals         = resid_tbl;
      s.Diagnostics       = diag_tbl;
      s.ModelFitVsNullModel = mfvnm;
      s.Steps             = steps;

      obj = LinearModel (s);

    endfunction

  endmethods

  ## Public methods
  methods (Access = public)

    function cmdl = compact (obj)
      ## compact  Strip to CompactLinearModel.
      ##
      ## cmdl = compact (mdl)
      ##
      ## Returns a CompactLinearModel with all 23 CompactLinearModel
      ## properties copied.  LinearModel-only properties (Variables,
      ## Residuals, Diagnostics, Fitted, ObservationInfo, ObservationNames,
      ## ModelFitVsNullModel, Steps) are NOT included.

      s = struct ();
      s.Coefficients             = obj.Coefficients;
      s.CoefficientNames         = obj.CoefficientNames;
      s.CoefficientCovariance    = obj.CoefficientCovariance;
      s.NumCoefficients          = obj.NumCoefficients;
      s.NumEstimatedCoefficients = obj.NumEstimatedCoefficients;
      s.DFE                      = obj.DFE;
      s.SSE                      = obj.SSE;
      s.SSR                      = obj.SSR;
      s.SST                      = obj.SST;
      s.MSE                      = obj.MSE;
      s.RMSE                     = obj.RMSE;
      s.Rsquared                 = obj.Rsquared;
      s.LogLikelihood            = obj.LogLikelihood;
      s.ModelCriterion           = obj.ModelCriterion;
      s.Robust                   = obj.Robust;
      s.Formula                  = obj.Formula;
      s.NumObservations          = obj.NumObservations;
      s.NumPredictors            = obj.NumPredictors;
      s.NumVariables             = obj.NumVariables;
      s.PredictorNames           = obj.PredictorNames;
      s.ResponseName             = obj.ResponseName;
      s.VariableInfo             = obj.VariableInfo;
      s.VariableNames            = obj.VariableNames;
      cmdl = CompactLinearModel (s);
    endfunction

    function NewMdl = addTerms (obj, terms)
      ## addTerms  Add terms to the model and refit.
      ##
      ## NewMdl = addTerms (mdl, terms)
      ##
      ## TERMS can be:
      ##   - A string giving the Wilkinson notation for terms to add
      ##     (e.g. 'x2' or 'x1:x2')
      ##   - A terms matrix row(s) to add (binary matrix)
      ##
      ## Returns a NEW LinearModel; does not modify the original.

      if (nargin < 2)
        error ("LinearModel.addTerms: TERMS argument is required.");
      endif

      ## Get current formula info
      current_terms = obj.Formula.Terms;
      current_formula = obj.Formula;
      p_names = obj.PredictorNames;

      if (ischar (terms) || isstring (terms))
        ## Parse the Wilkinson term string into term rows
        terms = char (terms);
        new_rows = __parse_term_string__ (terms, p_names);
      elseif (isnumeric (terms))
        new_rows = terms;
      else
        error ("LinearModel.addTerms: TERMS must be a string or numeric matrix.");
      endif

      ## Pad with response column(s) if Terms matrix is wider than new_rows
      n_cols_terms = columns (current_terms);
      if (columns (new_rows) < n_cols_terms)
        new_rows = [new_rows, zeros(rows (new_rows), ...
                     n_cols_terms - columns (new_rows))];
      endif

      ## Merge: add new rows to existing terms (skip duplicates)
      merged_terms = current_terms;
      for i = 1:rows (new_rows)
        is_dup = false;
        for j = 1:rows (merged_terms)
          if (isequal (new_rows(i, :), merged_terms(j, :)))
            is_dup = true;
            break;
          endif
        endfor
        if (! is_dup)
          merged_terms = [merged_terms; new_rows(i, :)];
        endif
      endfor

      ## Sort terms: intercept first, then by order (sum of row)
      term_order = sum (merged_terms, 2);
      [~, sort_idx] = sort (term_order);
      merged_terms = merged_terms(sort_idx, :);

      ## Rebuild design matrix
      y_vec = __get_y__ (obj);
      w_vec = __get_weights__ (obj);

      ## Build new design matrix (Phase 3B: returns coef names for categoricals)
      [X_new, new_coef_names] = __build_X_from_terms__ (obj, merged_terms);

      ## Refit
      fit_s = lm_fit_engine (X_new, y_vec, w_vec);

      ## Update formula (convert to struct for modification; fromFit re-wraps)
      if (isa (current_formula, "LinearFormula"))
        new_formula = struct ("ResponseName", current_formula.ResponseName, ...
                              "LinearPredictor", current_formula.LinearPredictor, ...
                              "Terms", current_formula.Terms, ...
                              "HasIntercept", current_formula.HasIntercept, ...
                              "InModel", current_formula.InModel, ...
                              "CoefficientNames", {current_formula.CoefficientNames}, ...
                              "PredictorNames", {current_formula.PredictorNames}, ...
                              "VariableNames", {current_formula.VariableNames});
      else
        new_formula = current_formula;
      endif
      new_formula.Terms = merged_terms;
      new_formula.CoefficientNames = new_coef_names;
      new_formula.LinearPredictor = __terms_to_formula_str__ ( ...
                                      merged_terms, p_names, ...
                                      new_formula.HasIntercept);

      ## Build new LinearModel
      subset = obj.ObservationInfo.Subset;
      NewMdl = LinearModel.fromFit (fit_s, X_new, y_vec, new_formula, ...
                 p_names, obj.ResponseName, obj.VariableNames, ...
                 obj.VariableInfo, obj.ObservationNames, ...
                 obj.ObservationInfo, w_vec, ...
                 obj.ObservationInfo.Excluded, ...
                 obj.ObservationInfo.Missing, obj.Steps, ...
                 obj.Variables);

    endfunction

    function NewMdl = removeTerms (obj, terms)
      ## removeTerms  Remove terms from the model and refit.
      ##
      ## NewMdl = removeTerms (mdl, terms)
      ##
      ## TERMS can be:
      ##   - A string giving the Wilkinson notation for terms to remove
      ##   - A terms matrix row(s) to remove (binary matrix)
      ##
      ## Returns a NEW LinearModel; does not modify the original.

      if (nargin < 2)
        error ("LinearModel.removeTerms: TERMS argument is required.");
      endif

      ## Get current formula info
      current_terms = obj.Formula.Terms;
      current_formula = obj.Formula;
      p_names = obj.PredictorNames;

      if (ischar (terms) || isstring (terms))
        terms = char (terms);
        remove_rows = __parse_term_string__ (terms, p_names);
      elseif (isnumeric (terms))
        remove_rows = terms;
      else
        error ("LinearModel.removeTerms: TERMS must be a string or numeric.");
      endif

      ## Pad with response column(s) if Terms matrix is wider than remove_rows
      n_cols_terms = columns (current_terms);
      if (columns (remove_rows) < n_cols_terms)
        remove_rows = [remove_rows, zeros(rows (remove_rows), ...
                        n_cols_terms - columns (remove_rows))];
      endif

      ## Remove matching rows from current terms
      keep_mask = true (rows (current_terms), 1);
      for i = 1:rows (remove_rows)
        for j = 1:rows (current_terms)
          if (isequal (remove_rows(i, :), current_terms(j, :)))
            keep_mask(j) = false;
            break;
          endif
        endfor
      endfor

      if (sum (keep_mask) == rows (current_terms))
        warning ("LinearModel:removeTerms", ...
                 "LinearModel.removeTerms: no matching terms found.");
      endif

      merged_terms = current_terms(keep_mask, :);

      ## Rebuild design matrix
      y_vec = __get_y__ (obj);
      w_vec = __get_weights__ (obj);
      [X_new, new_coef_names] = __build_X_from_terms__ (obj, merged_terms);

      ## Refit
      fit_s = lm_fit_engine (X_new, y_vec, w_vec);

      ## Update formula (convert to struct for modification; fromFit re-wraps)
      if (isa (current_formula, "LinearFormula"))
        new_formula = struct ("ResponseName", current_formula.ResponseName, ...
                              "LinearPredictor", current_formula.LinearPredictor, ...
                              "Terms", current_formula.Terms, ...
                              "HasIntercept", current_formula.HasIntercept, ...
                              "InModel", current_formula.InModel, ...
                              "CoefficientNames", {current_formula.CoefficientNames}, ...
                              "PredictorNames", {current_formula.PredictorNames}, ...
                              "VariableNames", {current_formula.VariableNames});
      else
        new_formula = current_formula;
      endif
      new_formula.Terms = merged_terms;
      new_formula.CoefficientNames = new_coef_names;
      new_formula.LinearPredictor = __terms_to_formula_str__ ( ...
                                      merged_terms, p_names, ...
                                      new_formula.HasIntercept);

      ## Build new LinearModel
      NewMdl = LinearModel.fromFit (fit_s, X_new, y_vec, new_formula, ...
                 p_names, obj.ResponseName, obj.VariableNames, ...
                 obj.VariableInfo, obj.ObservationNames, ...
                 obj.ObservationInfo, w_vec, ...
                 obj.ObservationInfo.Excluded, ...
                 obj.ObservationInfo.Missing, obj.Steps, ...
                 obj.Variables);

    endfunction

    function NewMdl = step (obj, varargin)
      ## step  Take one or more stepwise steps.
      ##
      ## NewMdl = step (mdl)
      ## NewMdl = step (mdl, Name, Value)
      ##
      ## Name-Value pairs:
      ##   'Criterion'  - 'sse' (default), 'aic', 'bic'
      ##   'PEnter'     - p-value threshold for adding (default 0.05)
      ##   'PRemove'    - p-value threshold for removing (default 0.10)
      ##   'NSteps'     - maximum number of steps (default 1)
      ##   'Verbose'    - 0 (silent), 1 (default, print steps)
      ##   'Upper'      - maximum model (default 'interactions')
      ##   'Lower'      - minimum model (default 'constant')

      ## Parse name-value pairs
      criterion  = "sse";
      p_enter    = 0.05;
      p_remove   = 0.10;
      n_steps    = 1;
      verbose    = 1;
      upper_spec = "interactions";
      lower_spec = "constant";

      i = 1;
      while (i <= numel (varargin))
        if (! ischar (varargin{i}))
          error ("LinearModel.step: expected parameter name string.");
        endif
        switch (lower (varargin{i}))
          case "criterion"
            criterion = lower (varargin{i+1});
            i = i + 2;
          case "penter"
            p_enter = varargin{i+1};
            i = i + 2;
          case "premove"
            p_remove = varargin{i+1};
            i = i + 2;
          case "nsteps"
            n_steps = varargin{i+1};
            i = i + 2;
          case "verbose"
            verbose = varargin{i+1};
            i = i + 2;
          case "upper"
            upper_spec = varargin{i+1};
            i = i + 2;
          case "lower"
            lower_spec = varargin{i+1};
            i = i + 2;
          otherwise
            error ("LinearModel.step: unknown parameter '%s'.", varargin{i});
        endswitch
      endwhile

      p_names = obj.PredictorNames;
      n_pred = numel (p_names);

      ## Build upper/lower terms matrices
      upper_terms = [__step_spec_to_terms__(upper_spec, n_pred), ...
                     zeros(rows (__step_spec_to_terms__(upper_spec, n_pred)), 1)];
      lower_terms = [__step_spec_to_terms__(lower_spec, n_pred), ...
                     zeros(rows (__step_spec_to_terms__(lower_spec, n_pred)), 1)];

      ## Build cat_flags from VariableInfo
      cat_flags = false (1, n_pred + 1);
      if ((isstruct (obj.VariableInfo) && isfield (obj.VariableInfo, "IsCategorical")) || ...
          isa (obj.VariableInfo, "table"))
        for k = 1:n_pred
          vi_idx = strcmp (obj.VariableNames, p_names{k});
          if (any (vi_idx))
            cat_flags(k) = obj.VariableInfo.IsCategorical(vi_idx);
          endif
        endfor
      endif

      current_terms = obj.Formula.Terms;
      use_ftest = strcmp (criterion, "sse");

      ## Initialize or extend Steps.History
      if (isempty (obj.Steps))
        start_lp = obj.Formula.LinearPredictor;
        log_action  = {'Start'};
        log_term    = {start_lp};
        log_terms   = {current_terms};
        log_df      = obj.NumCoefficients;
        log_del_df  = NaN;
        log_fstat   = NaN;
        log_pval    = NaN;
        log_critval = __get_criterion__ (obj, criterion);
        n_log = 1;
      else
        prev_h = obj.Steps.History;
        n_log = numel (prev_h.Action);
        log_action  = prev_h.Action;
        log_term    = cell (n_log, 1);
        for k = 1:n_log
          if (iscell (prev_h.TermName{k}))
            log_term{k} = prev_h.TermName{k}{1};
          else
            log_term{k} = prev_h.TermName{k};
          endif
        endfor
        log_terms   = prev_h.Terms;
        log_df      = prev_h.DF;
        log_del_df  = prev_h.delDF;
        if (use_ftest && isfield (prev_h, "FStat"))
          log_fstat = prev_h.FStat;
          log_pval  = prev_h.pValue;
        else
          log_fstat = NaN (n_log, 1);
          log_pval  = NaN (n_log, 1);
        endif
        log_critval = NaN (n_log, 1);
      endif

      ## Check for empty candidate set
      all_cands = __stepwiselm_candidates__ (current_terms, lower_terms, ...
                                              upper_terms, cat_flags, p_names);
      if (isempty (all_cands) && verbose >= 1)
        fprintf ("No terms to add to or remove from initial model.\n");
      endif

      ## Search loop
      n_taken = 0;
      current_mdl = obj;

      while (n_taken < n_steps)
        made_change = false;

        ## Phase A: BACKWARD SCAN (remove until nothing qualifies)
        while (n_taken < n_steps)
          candidates = __stepwiselm_candidates__ ( ...
              current_terms, lower_terms, upper_terms, cat_flags, p_names);
          rem_idx = [];
          for ci = 1:numel (candidates)
            if (strcmp (candidates(ci).action, "remove"))
              rem_idx(end+1) = ci;
            endif
          endfor
          if (isempty (rem_idx)); break; endif

          best_ri = 0; best_pval_r = -Inf; best_fstat_r = 0; best_df_r = 0;
          for ri = 1:numel (rem_idx)
            c = candidates(rem_idx(ri));
            f_stat = 0; f_pval = 1;
            try
              test_mdl = removeTerms (current_mdl, c.term_row(1:n_pred));
              ss_d = test_mdl.SSE - current_mdl.SSE;
              df_d = test_mdl.DFE - current_mdl.DFE;
              if (df_d > 0 && current_mdl.MSE > 0)
                f_stat = (ss_d / df_d) / current_mdl.MSE;
                f_pval = 1 - fcdf (f_stat, df_d, current_mdl.DFE);
              endif
            catch; continue; end_try_catch
            if (f_pval > p_remove && f_pval > best_pval_r)
              best_ri = rem_idx(ri); best_pval_r = f_pval;
              best_fstat_r = f_stat; best_df_r = df_d;
            endif
          endfor
          if (best_ri == 0); break; endif

          best_c = candidates(best_ri);
          current_mdl = removeTerms (current_mdl, best_c.term_row(1:n_pred));
          mask = true (rows (current_terms), 1);
          for j = 1:rows (current_terms)
            if (isequal (current_terms(j,:), best_c.term_row))
              mask(j) = false; break;
            endif
          endfor
          current_terms = current_terms(mask, :);
          n_taken++; made_change = true;
          n_log++;
          log_action{n_log} = "Remove"; log_term{n_log} = best_c.term_name;
          log_terms{n_log} = current_terms;
          log_df(n_log) = current_mdl.NumCoefficients;
          log_del_df(n_log) = -best_df_r;
          log_fstat(n_log) = best_fstat_r; log_pval(n_log) = best_pval_r;
          log_critval(n_log) = __get_criterion__ (current_mdl, criterion);
          if (verbose >= 1)
            fprintf ("%d. Removing %s, FStat = %g, pValue = %g\n", ...
                     n_taken, best_c.term_name, best_fstat_r, best_pval_r);
          endif
        endwhile
        if (n_taken >= n_steps); break; endif

        ## Phase B: FORWARD SCAN (add one term)
        candidates = __stepwiselm_candidates__ ( ...
            current_terms, lower_terms, upper_terms, cat_flags, p_names);
        add_idx = [];
        for ci = 1:numel (candidates)
          if (strcmp (candidates(ci).action, "add"))
            add_idx(end+1) = ci;
          endif
        endfor
        if (isempty (add_idx)); break; endif

        best_ai = 0; best_pval_a = Inf; best_fstat_a = 0;
        best_df_a = 0; best_add_mdl = [];
        for ai = 1:numel (add_idx)
          c = candidates(add_idx(ai));
          f_stat = 0; f_pval = 1;
          try
            test_mdl = addTerms (current_mdl, c.term_row(1:n_pred));
            ss_d = current_mdl.SSE - test_mdl.SSE;
            df_d = current_mdl.DFE - test_mdl.DFE;
            if (df_d > 0 && test_mdl.MSE > 0)
              f_stat = (ss_d / df_d) / test_mdl.MSE;
              f_pval = 1 - fcdf (f_stat, df_d, test_mdl.DFE);
            endif
          catch; continue; end_try_catch
          if (f_pval < p_enter && f_pval < best_pval_a)
            best_ai = add_idx(ai); best_pval_a = f_pval;
            best_fstat_a = f_stat; best_df_a = df_d;
            best_add_mdl = test_mdl;
          endif
        endfor
        if (best_ai == 0)
          if (! made_change); break; endif
          break;
        endif

        best_c = candidates(best_ai);
        current_mdl = best_add_mdl;
        current_terms = [current_terms; best_c.term_row];
        term_order = sum (current_terms(:, 1:n_pred), 2);
        [~, si] = sort (term_order);
        current_terms = current_terms(si, :);
        n_taken++; made_change = true;
        n_log++;
        log_action{n_log} = "Add"; log_term{n_log} = best_c.term_name;
        log_terms{n_log} = current_terms;
        log_df(n_log) = current_mdl.NumCoefficients;
        log_del_df(n_log) = best_df_a;
        log_fstat(n_log) = best_fstat_a; log_pval(n_log) = best_pval_a;
        log_critval(n_log) = __get_criterion__ (current_mdl, criterion);
        if (verbose >= 1)
          fprintf ("%d. Adding %s, FStat = %g, pValue = %g\n", ...
                   n_taken, best_c.term_name, best_fstat_a, best_pval_a);
        endif
      endwhile

      ## Build Steps struct
      resp_n = obj.ResponseName;
      if (isempty (obj.Steps))
        start_formula = [resp_n, " ~ ", log_term{1}];
      else
        start_formula = obj.Steps.Start;
      endif
      Steps.Start = start_formula; Steps.Lower = lower_spec;
      Steps.Upper = upper_spec; Steps.Criterion = upper (criterion);
      Steps.PEnter = p_enter; Steps.PRemove = p_remove;
      History.Action = log_action(:);
      History.TermName = cell (n_log, 1);
      for k = 1:n_log
        History.TermName{k} = {log_term{k}};
      endfor
      History.Terms = log_terms(:);
      History.DF = log_df(:);  History.delDF = log_del_df(:);
      if (use_ftest)
        History.FStat = log_fstat(:); History.pValue = log_pval(:);
      else
        History.(upper (criterion)) = log_critval(:);
      endif
      Steps.History = History;
      current_mdl.Steps = Steps;
      NewMdl = current_mdl;

    endfunction

    function tbl = anova (obj, varargin)
      ## anova  ANOVA table for LinearModel (Partial SS / Type III).
      ##
      ## tbl = anova (mdl)
      ## tbl = anova (mdl, 'component')
      ##
      ## Returns a struct with fields:
      ##   SumSq  - sum of squares for each term + Error
      ##   DF     - degrees of freedom
      ##   MeanSq - mean squares
      ##   F      - F-statistic (NaN for Error row)
      ##   pValue - p-value (NaN for Error row)
      ##   RowNames - cell of term names + 'Error'

      ## Parse optional anova type argument (only 'component' for now)

      current_terms = obj.Formula.Terms;
      p_names = obj.PredictorNames;

      ## Identify non-intercept terms
      is_intercept = all (current_terms == 0, 2);
      term_idx = find (~is_intercept);
      n_terms = numel (term_idx);

      ## Compute partial SS for each term (Type III):
      ## SS(term_k) = SSE(model without term_k) - SSE(full model)
      SumSq  = zeros (n_terms + 1, 1);
      DF     = zeros (n_terms + 1, 1);
      MeanSq = zeros (n_terms + 1, 1);
      F_val  = NaN (n_terms + 1, 1);
      p_val  = NaN (n_terms + 1, 1);
      term_names = cell (n_terms + 1, 1);

      y_vec = __get_y__ (obj);
      w_vec = __get_weights__ (obj);

      for k = 1:n_terms
        ## Remove term k and refit
        reduced_terms = current_terms;
        reduced_terms(term_idx(k), :) = [];

        [X_reduced, ~] = __build_X_from_terms__ (obj, reduced_terms);
        fit_reduced = lm_fit_engine (X_reduced, y_vec, w_vec);

        ## Partial SS
        ss_k = fit_reduced.sse - obj.SSE;
        df_k = fit_reduced.dfe - obj.DFE;

        SumSq(k) = ss_k;
        DF(k) = df_k;
        MeanSq(k) = ss_k / df_k;
        F_val(k) = MeanSq(k) / obj.MSE;
        p_val(k) = 1 - fcdf (F_val(k), df_k, obj.DFE);

        ## Term name
        idx = find (current_terms(term_idx(k), :));
        parts = {};
        for j = 1:numel (idx)
          parts{j} = p_names{idx(j)};
        endfor
        ## Handle power terms
        row = current_terms(term_idx(k), :);
        power_parts = {};
        for j = 1:numel (idx)
          if (row(idx(j)) == 1)
            power_parts{end+1} = p_names{idx(j)};
          else
            power_parts{end+1} = sprintf ("%s^%d", p_names{idx(j)}, ...
                                          row(idx(j)));
          endif
        endfor
        term_names{k} = strjoin (power_parts, ":");
      endfor

      ## Error row
      SumSq(end)  = obj.SSE;
      DF(end)     = obj.DFE;
      MeanSq(end) = obj.MSE;
      F_val(end)  = NaN;
      p_val(end)  = NaN;
      term_names{end} = "Error";

      tbl = struct ("SumSq", SumSq, "DF", DF, "MeanSq", MeanSq, ...
                    "F", F_val, "pValue", p_val, "RowNames", {term_names});

    endfunction

  endmethods

endclassdef

## ============================================================
## Private helper functions
## ============================================================

function rows_out = __parse_term_string__ (term_str, pred_names)
  ## Parse a Wilkinson-style term string into binary term matrix rows.
  ## e.g. 'x2' -> [0 1 0], 'x1:x2' -> [1 1 0]

  n_pred = numel (pred_names);

  ## Split by '+' to get individual terms
  terms = strtrim (strsplit (term_str, "+"));
  rows_out = zeros (numel (terms), n_pred);

  for i = 1:numel (terms)
    ## Split by ':' for interaction terms
    parts = strtrim (strsplit (terms{i}, ":"));
    for j = 1:numel (parts)
      ## Handle power notation (e.g., 'x1^2')
      p = parts{j};
      caret = strfind (p, "^");
      if (! isempty (caret))
        base = strtrim (p(1:caret-1));
        pwr = str2double (p(caret+1:end));
      else
        base = p;
        pwr = 1;
      endif
      ## Find the predictor index
      idx = find (strcmp (pred_names, base));
      if (isempty (idx))
        error ("LinearModel: unrecognized predictor name '%s'.", base);
      endif
      rows_out(i, idx(1)) = pwr;
    endfor
  endfor
endfunction

function [X, coef_names] = __rebuild_X__ (obj)
  ## Rebuild the design matrix from stored Formula + predictor data.
  [X, coef_names] = __build_X_from_terms__ (obj, obj.Formula.Terms);
endfunction

function [X, coef_names] = __build_X_from_terms__ (obj, terms)
  ## Build design matrix from terms matrix and stored predictor data.
  ##
  ## Phase 3B: handles categorical variables via reference (corner-point)
  ## dummy coding.  For categorical predictors with K levels, K-1 dummy
  ## columns are generated (dropping the first level when intercept is
  ## present).  Interactions between categorical and numeric variables
  ## produce a Cartesian product of dummy columns x numeric columns.
  ##
  ## Returns:
  ##   X          - n_obs x p design matrix
  ##   coef_names - 1 x p cell of coefficient name strings

  ## Get predictor data for subset observations
  subset = obj.ObservationInfo.Subset;
  p_names = obj.PredictorNames;
  n_pred = numel (p_names);

  ## Detect intercept
  has_intercept = any (all (terms == 0, 2));

  ## Extract predictor columns and detect categorical status
  n_obs = sum (subset);
  pred_raw = cell (1, n_pred);
  is_categorical = false (1, n_pred);
  cat_levels = cell (1, n_pred);
  cat_indices = cell (1, n_pred);

  for k = 1:n_pred
    pname = p_names{k};
    if (! isempty (obj.Variables) && isa (obj.Variables, "table"))
      col_full = obj.Variables.(pname);
    elseif (! isempty (obj.Variables) && isstruct (obj.Variables))
      if (isfield (obj.Variables, pname))
        col_full = obj.Variables.(pname);
      else
        error ("LinearModel: predictor '%s' not found in Variables.", pname);
      endif
    else
      error (strcat ("LinearModel: cannot rebuild design matrix", ...
                     " without stored Variables data."));
    endif

    ## Extract subset rows
    if (iscell (col_full) || iscellstr (col_full) || isstring (col_full))
      col = col_full(subset);
    elseif (isa (col_full, "categorical"))
      col = col_full(subset);
    else
      col = col_full(subset, :);
    endif

    pred_raw{k} = col;

    ## Detect categorical: use VariableInfo.IsCategorical as primary source
    ## (Fix 4, Phase 3B), with type-sniffing as a secondary fallback.
    vi_is_cat = false;
    if (! isempty (obj.VariableInfo) && ...
        ((isstruct (obj.VariableInfo) && isfield (obj.VariableInfo, "IsCategorical")) || ...
         isa (obj.VariableInfo, "table")))
      vi_idx = strcmp (obj.VariableNames, pname);
      if (any (vi_idx))
        vi_is_cat = obj.VariableInfo.IsCategorical(vi_idx);
      endif
    endif

    ## Type-sniff fallback (cellstr, string, categorical)
    type_is_cat = iscellstr (col) || isstring (col) || isa (col, "categorical");

    if (vi_is_cat || type_is_cat)
      is_categorical(k) = true;
      if (isa (col, "categorical"))
        levs = categories (col);
        cat_levels{k} = cellstr (levs);
        [~, idx] = ismember (col, levs);
        cat_indices{k} = idx;
      elseif (iscellstr (col))
        [u, ~, idx] = unique (col);
        cat_levels{k} = u;
        cat_indices{k} = idx;
      elseif (isstring (col))
        col_c = cellstr (col);
        [u, ~, idx] = unique (col_c);
        cat_levels{k} = u;
        cat_indices{k} = idx;
      else
        ## Numeric column flagged categorical via VariableInfo.IsCategorical
        ## Treat unique numeric values as levels (sorted)
        [u, ~, idx] = unique (col);
        cat_levels{k} = cellfun (@num2str, num2cell (u), "UniformOutput", false);
        cat_indices{k} = idx;
      endif
    else
      is_categorical(k) = false;
    endif
  endfor

  ## Build X column by column from terms
  X = [];
  coef_names = {};

  for i = 1:rows (terms)
    row = terms(i, :);

    if (all (row == 0))
      ## Intercept
      X = [X, ones(n_obs, 1)];
      coef_names{end+1} = "(Intercept)";
      continue;
    endif

    ## Find active predictors for this term
    active_idx = find (row);

    ## Build the term's columns via Cartesian product approach.
    ## Start with a single column of ones, then for each active predictor:
    ##   - If numeric: multiply current block by (data.^power)
    ##   - If categorical: Cartesian product with dummy columns
    current_block = ones (n_obs, 1);
    current_names = {''};

    for ai = 1:numel (active_idx)
      v = active_idx(ai);
      pname = p_names{v};
      pwr = row(v);

      if (is_categorical(v))
        ## Build dummy columns (reference coding)
        levs = cat_levels{v};
        idx = cat_indices{v};
        n_lev = numel (levs);

        if (has_intercept)
          start_lev = 2;
        else
          start_lev = 1;
        endif

        n_dum = n_lev - start_lev + 1;
        dummies = zeros (n_obs, n_dum);
        dum_names = cell (1, n_dum);

        for L = start_lev:n_lev
          col_idx = L - start_lev + 1;
          dummies(:, col_idx) = double (idx == L);
          dum_names{col_idx} = sprintf ("%s_%s", pname, levs{L});
        endfor

        ## Cartesian product of current block and dummies
        next_block = [];
        next_names = {};
        for c1 = 1:columns (current_block)
          for c2 = 1:columns (dummies)
            next_block = [next_block, current_block(:, c1) .* dummies(:, c2)];
            n1 = current_names{c1};
            n2 = dum_names{c2};
            if (isempty (n1))
              next_names{end+1} = n2;
            else
              next_names{end+1} = [n1, ":", n2];
            endif
          endfor
        endfor
        current_block = next_block;
        current_names = next_names;

      else
        ## Numeric predictor: multiply current block by data.^power
        col_data = pred_raw{v};
        if (pwr > 1)
          col_data = col_data .^ pwr;
          disp_name = sprintf ("%s^%d", pname, pwr);
        else
          disp_name = pname;
        endif

        current_block = bsxfun (@times, current_block, col_data);
        for cn = 1:numel (current_names)
          if (isempty (current_names{cn}))
            current_names{cn} = disp_name;
          else
            current_names{cn} = [current_names{cn}, ":", disp_name];
          endif
        endfor
      endif
    endfor

    X = [X, current_block];
    coef_names = [coef_names, current_names];
  endfor
endfunction

function y = __get_y__ (obj)
  ## Get response vector for subset observations.
  subset = obj.ObservationInfo.Subset;
  resp = obj.ResponseName;

  if (! isempty (obj.Variables) && isa (obj.Variables, "table"))
    y_full = obj.Variables.(resp);
    y = y_full(subset);
  elseif (! isempty (obj.Variables) && isstruct (obj.Variables))
    if (isfield (obj.Variables, resp))
      y_full = obj.Variables.(resp);
      y = y_full(subset);
    else
      error ("LinearModel: response '%s' not found in Variables.", resp);
    endif
  else
    ## Reconstruct from Fitted + Residuals
    y = obj.Fitted(subset) + obj.Residuals.Raw(subset);
  endif
endfunction

function w = __get_weights__ (obj)
  ## Get weights vector for subset observations, or [] if unweighted.
  subset = obj.ObservationInfo.Subset;
  w_all = obj.ObservationInfo.Weights;
  if (all (w_all == 1))
    w = [];
  else
    w = w_all(subset);
  endif
endfunction

function score = __get_criterion__ (mdl, criterion)
  ## Get model selection criterion value.
  switch (criterion)
    case "sse"
      score = mdl.SSE;
    case "aic"
      score = mdl.ModelCriterion.AIC;
    case "bic"
      score = mdl.ModelCriterion.BIC;
    case "rsquared"
      score = -mdl.Rsquared.Ordinary;   ## negate so lower = better
    case "adjrsquared"
      score = -mdl.Rsquared.Adjusted;
    otherwise
      score = mdl.SSE;
  endswitch
endfunction

function names = __terms_to_coef_names__ (terms, pred_names, has_intercept)
  ## Convert terms matrix to coefficient names.
  names = {};
  for i = 1:rows (terms)
    row = terms(i, :);
    if (all (row == 0))
      names{end+1} = "(Intercept)";
    else
      parts = {};
      idx = find (row);
      for j = 1:numel (idx)
        if (row(idx(j)) == 1)
          parts{end+1} = pred_names{idx(j)};
        else
          parts{end+1} = sprintf ("%s^%d", pred_names{idx(j)}, row(idx(j)));
        endif
      endfor
      names{end+1} = strjoin (parts, ":");
    endif
  endfor
endfunction

function fstr = __terms_to_formula_str__ (terms, pred_names, has_intercept)
  ## Convert terms matrix to Wilkinson formula string (RHS only).
  parts = {};
  for i = 1:rows (terms)
    row = terms(i, :);
    if (all (row == 0))
      parts{end+1} = "1";
    else
      tparts = {};
      idx = find (row);
      for j = 1:numel (idx)
        if (row(idx(j)) == 1)
          tparts{end+1} = pred_names{idx(j)};
        else
          tparts{end+1} = sprintf ("%s^%d", pred_names{idx(j)}, row(idx(j)));
        endif
      endfor
      parts{end+1} = strjoin (tparts, ":");
    endif
  endfor
  fstr = strjoin (parts, " + ");
endfunction

## ============================================================
## Tests
## ============================================================

## Test: LinearModel class exists and is parseable
%!test
%! assert (exist ("LinearModel", "class") > 0 || ...
%!         exist ("LinearModel") > 0);

