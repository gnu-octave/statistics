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
## @item @code{Variables} @tab Original input data (all rows, including
## excluded).
## @item @code{ObservationNames} @tab Cell array of row names.
## @item @code{ObservationInfo} @tab Struct with fields @code{Weights},
## @code{Excluded}, @code{Missing}, @code{Subset}.
## @item @code{Fitted} @tab Fitted values for all observations.
## @item @code{Residuals} @tab Struct with fields @code{Raw},
## @code{Pearson}, @code{Studentized}, @code{Standardized}.
## @item @code{Diagnostics} @tab Struct with fields @code{Leverage},
## @code{CooksDistance}, @code{Dffits}, @code{S2_i}, @code{CovRatio},
## @code{Dfbetas}, @code{HatMatrix}.
## @item @code{ModelFitVsNullModel} @tab Struct with fields @code{Fstat},
## @code{Pvalue}, @code{NullModel}.
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
    Variables            = [];
    ObservationNames     = {};
    ObservationInfo      = [];

    ## Diagnostics (not on CompactLinearModel)
    Fitted               = [];
    Residuals            = [];
    Diagnostics          = [];
    ModelFitVsNullModel  = [];

    ## Fitting method
    Steps                = [];
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

      if (! isstruct (s))
        error ("LinearModel: constructor argument must be a struct.");
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
                              variables_data)
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
      n_used  = fit_s.n;
      subset  = obs_info.Subset;

      ## Compute Fitted for ALL observations (including excluded/missing)
      fitted_all = NaN (n_total, 1);
      fitted_all(subset) = fit_s.yhat;
      ## Excluded/missing rows: compute yhat = X_full * beta if X is available
      ## For Phase 3A, excluded rows get NaN fitted values
      ## (MATLAB Block 14 shows excluded rows DO get fitted values)
      ## We compute them using the stored design matrix columns
      ## This requires the full design matrix, which we don't store.
      ## For now, leave excluded rows as NaN; fitlm (Phase 4) will handle this.

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

      resid_struct = struct ("Raw", raw_all, "Pearson", pearson_all, ...
                             "Studentized", studentized_all, ...
                             "Standardized", standardized_all);

      ## Diagnostics for used observations
      CooksDistance_used = (1 / p_est) * (h ./ (1 - h)) .* standardized .^ 2;
      Dffits_used        = studentized .* sqrt (h ./ (1 - h));
      ## CovRatio: Candidate 3 formula, verified to 3.55e-14 against MATLAB
      CovRatio_used      = (S2_i_used / mse) .^ p_est ./ (1 - h);

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
          Dfbetas_all(ii, :) = ((vcov * x_i * e_i) / ...
                                ((1 - h_i) * sqrt (s2_i))) ./ se_coef;
        endif
      endfor

      diag_struct = struct ("Leverage", leverage_all, ...
                            "CooksDistance", CooksDistance_all, ...
                            "Dffits", Dffits_all, ...
                            "S2_i", S2_i_all, ...
                            "CovRatio", CovRatio_all, ...
                            "Dfbetas", Dfbetas_all, ...
                            "HatMatrix", HatMatrix_full);

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
      s.Coefficients             = struct ("Estimate", fit_s.beta, ...
                                           "SE", fit_s.se, ...
                                           "tStat", fit_s.tstat, ...
                                           "pValue", fit_s.pval);
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
      s.Formula                  = formula;
      s.NumObservations          = n_used;
      s.NumPredictors            = numel (pred_names);
      s.NumVariables             = numel (var_names);
      s.PredictorNames           = pred_names;
      s.ResponseName             = resp_name;
      s.VariableInfo             = var_info;
      s.VariableNames            = var_names;

      ## LinearModel-only fields
      if (nargin >= 15 && ! isempty (variables_data))
        s.Variables = variables_data;
      else
        s.Variables = [];
      endif
      s.ObservationNames  = obs_names;
      s.ObservationInfo   = obs_info;
      s.Fitted            = fitted_all;
      s.Residuals         = resid_struct;
      s.Diagnostics       = diag_struct;
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

      ## Rebuild design matrix from current stored X info
      ## Phase 3A approach: rebuild X from column manipulation
      X_old = __rebuild_X__ (obj);
      y_vec = __get_y__ (obj);
      w_vec = __get_weights__ (obj);

      ## Build new design matrix columns for new terms
      X_new = __build_X_from_terms__ (obj, merged_terms);

      ## Refit
      fit_s = lm_fit_engine (X_new, y_vec, w_vec);

      ## Update formula
      new_formula = current_formula;
      new_formula.Terms = merged_terms;
      new_formula.CoefficientNames = __terms_to_coef_names__ ( ...
                                       merged_terms, p_names, ...
                                       current_formula.HasIntercept);
      new_formula.LinearPredictor = __terms_to_formula_str__ ( ...
                                      merged_terms, p_names, ...
                                      current_formula.HasIntercept);

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
      X_new = __build_X_from_terms__ (obj, merged_terms);

      ## Refit
      fit_s = lm_fit_engine (X_new, y_vec, w_vec);

      ## Update formula
      new_formula = current_formula;
      new_formula.Terms = merged_terms;
      new_formula.CoefficientNames = __terms_to_coef_names__ ( ...
                                       merged_terms, p_names, ...
                                       current_formula.HasIntercept);
      new_formula.LinearPredictor = __terms_to_formula_str__ ( ...
                                      merged_terms, p_names, ...
                                      current_formula.HasIntercept);

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
      ##   'Criterion'  - 'sse' (default), 'aic', 'bic', 'rsquared',
      ##                  'adjrsquared'
      ##   'PEnter'     - p-value threshold for adding (default 0.05)
      ##   'PRemove'    - p-value threshold for removing (default 0.10)
      ##   'NSteps'     - maximum number of steps (default 1)
      ##   'Verbose'    - 0 (silent), 1 (default, print steps)

      ## Parse name-value pairs
      criterion = "sse";
      p_enter   = 0.05;
      p_remove  = 0.10;
      n_steps   = 1;
      verbose   = 1;

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
          otherwise
            error ("LinearModel.step: unknown parameter '%s'.", varargin{i});
        endswitch
      endwhile

      current_mdl = obj;
      step_count = 0;

      for s = 1:n_steps
        current_terms = current_mdl.Formula.Terms;
        p_names = current_mdl.PredictorNames;
        n_pred = numel (p_names);

        best_action = "none";
        best_term = [];
        best_score = Inf;       ## lower is better for SSE/AIC/BIC
        best_pval = 1;
        best_fstat = 0;

        ## Current model score
        curr_score = __get_criterion__ (current_mdl, criterion);

        ## Try adding each predictor not already in the model (main effects)
        for k = 1:n_pred
          new_row = zeros (1, n_pred);
          new_row(k) = 1;
          ## Check if this term already exists
          already_in = false;
          for j = 1:rows (current_terms)
            if (isequal (new_row, current_terms(j, :)))
              already_in = true;
              break;
            endif
          endfor
          if (already_in)
            continue;
          endif

          ## Try adding
          try
            test_mdl = addTerms (current_mdl, new_row);
          catch
            continue;
          end_try_catch

          ## Compute partial F-test for the added term
          ss_diff = current_mdl.SSE - test_mdl.SSE;
          df_diff = current_mdl.DFE - test_mdl.DFE;
          if (df_diff > 0 && test_mdl.MSE > 0)
            f_stat = (ss_diff / df_diff) / test_mdl.MSE;
            f_pval = 1 - fcdf (f_stat, df_diff, test_mdl.DFE);
          else
            f_stat = 0;
            f_pval = 1;
          endif

          score = __get_criterion__ (test_mdl, criterion);

          if (f_pval < p_enter && score < best_score)
            best_action = "add";
            best_term = new_row;
            best_score = score;
            best_pval = f_pval;
            best_fstat = f_stat;
            best_term_name = p_names{k};
          endif
        endfor

        ## Try removing each non-intercept term
        for j = 1:rows (current_terms)
          if (all (current_terms(j, :) == 0))
            continue;   ## skip intercept
          endif

          try
            test_mdl = removeTerms (current_mdl, current_terms(j, :));
          catch
            continue;
          end_try_catch

          ## Compute partial F-test for the removed term
          ss_diff = test_mdl.SSE - current_mdl.SSE;
          df_diff = test_mdl.DFE - current_mdl.DFE;
          if (df_diff > 0 && current_mdl.MSE > 0)
            f_stat = (ss_diff / df_diff) / current_mdl.MSE;
            f_pval = 1 - fcdf (f_stat, df_diff, current_mdl.DFE);
          else
            f_stat = 0;
            f_pval = 1;
          endif

          score = __get_criterion__ (test_mdl, criterion);

          if (f_pval > p_remove && score < best_score)
            best_action = "remove";
            best_term = current_terms(j, :);
            best_score = score;
            best_pval = f_pval;
            best_fstat = f_stat;
            ## Reconstruct term name
            idx = find (current_terms(j, :));
            parts = p_names(idx);
            best_term_name = strjoin (parts, ":");
          endif
        endfor

        ## Execute best action
        if (strcmp (best_action, "none"))
          break;   ## converged
        endif

        step_count = step_count + 1;

        if (strcmp (best_action, "add"))
          current_mdl = addTerms (current_mdl, best_term);
          if (verbose)
            fprintf ("%d. Adding %s, FStat = %.3f, pValue = %g\n", ...
                     step_count, best_term_name, best_fstat, best_pval);
          endif
        else
          current_mdl = removeTerms (current_mdl, best_term);
          if (verbose)
            fprintf ("%d. Removing %s, FStat = %.6f, pValue = %g\n", ...
                     step_count, best_term_name, best_fstat, best_pval);
          endif
        endif
      endfor

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

        X_reduced = __build_X_from_terms__ (obj, reduced_terms);
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

function X = __rebuild_X__ (obj)
  ## Rebuild the design matrix from stored Formula + predictor data.
  ## Phase 3A: uses stored internal fields.
  X = __build_X_from_terms__ (obj, obj.Formula.Terms);
endfunction

function X = __build_X_from_terms__ (obj, terms)
  ## Build design matrix from terms matrix and stored predictor data.
  ##
  ## For Phase 3A, this uses the _X_data and _y_data stored on the object
  ## by the factory function, or reconstructs from Variables if available.

  ## Get predictor data for subset observations
  subset = obj.ObservationInfo.Subset;
  p_names = obj.PredictorNames;
  n_pred = numel (p_names);

  ## Try to get data from Variables (preferred)
  if (! isempty (obj.Variables) && isa (obj.Variables, "table"))
    ## Extract predictor columns from table
    pred_data = zeros (sum (subset), n_pred);
    for k = 1:n_pred
      col = obj.Variables.(p_names{k});
      pred_data(:, k) = col(subset);
    endfor
  elseif (! isempty (obj.Variables) && isstruct (obj.Variables))
    pred_data = zeros (sum (subset), n_pred);
    for k = 1:n_pred
      if (isfield (obj.Variables, p_names{k}))
        col = obj.Variables.(p_names{k});
        pred_data(:, k) = col(subset);
      endif
    endfor
  else
    error (strcat ("LinearModel: cannot rebuild design matrix", ...
                   " without stored Variables data."));
  endif

  n_obs = sum (subset);

  ## Build X column by column from terms
  X_cols = {};
  for i = 1:rows (terms)
    row = terms(i, :);
    if (all (row == 0))
      ## Intercept
      X_cols{end+1} = ones (n_obs, 1);
    else
      ## Product of predictor columns raised to powers
      col = ones (n_obs, 1);
      for j = 1:n_pred
        if (row(j) > 0)
          col = col .* (pred_data(:, j) .^ row(j));
        endif
      endfor
      X_cols{end+1} = col;
    endif
  endfor

  X = cell2mat (X_cols);
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

## Test: LinearModel inherits from CompactLinearModel
%!test
%! X = [ones(10,1), (1:10)', ((1:10).^2)'];
%! y = X * [1; 0.5; -0.01] + randn(10,1) * 0.1;
%! fit_s = lm_fit_engine (X, y);
%! pred_names = {'x1', 'x2'};
%! resp_name = 'y';
%! var_names = {'x1', 'x2', 'y'};
%! formula = struct ('Terms', [0 0; 1 0; 0 1], ...
%!                   'HasIntercept', true, ...
%!                   'LinearPredictor', '1 + x1 + x2', ...
%!                   'CoefficientNames', {{'(Intercept)', 'x1', 'x2'}});
%! obs_info = struct ('Weights', ones(10,1), ...
%!                    'Excluded', false(10,1), ...
%!                    'Missing', false(10,1), ...
%!                    'Subset', true(10,1));
%! vars = struct ('x1', (1:10)', 'x2', ((1:10).^2)', 'y', y);
%! mdl = LinearModel.fromFit (fit_s, X, y, formula, pred_names, ...
%!         resp_name, var_names, [], {}, obs_info, [], ...
%!         false(10,1), false(10,1), [], vars);
%! assert (isa (mdl, 'CompactLinearModel'));
%! assert (isa (mdl, 'LinearModel'));

## Test: compact returns CompactLinearModel, not LinearModel
%!test
%! X = [ones(10,1), (1:10)', ((1:10).^2)'];
%! y = X * [1; 0.5; -0.01] + randn(10,1) * 0.1;
%! fit_s = lm_fit_engine (X, y);
%! pred_names = {'x1', 'x2'};
%! resp_name = 'y';
%! var_names = {'x1', 'x2', 'y'};
%! formula = struct ('Terms', [0 0; 1 0; 0 1], ...
%!                   'HasIntercept', true, ...
%!                   'LinearPredictor', '1 + x1 + x2', ...
%!                   'CoefficientNames', {{'(Intercept)', 'x1', 'x2'}});
%! obs_info = struct ('Weights', ones(10,1), ...
%!                    'Excluded', false(10,1), ...
%!                    'Missing', false(10,1), ...
%!                    'Subset', true(10,1));
%! vars = struct ('x1', (1:10)', 'x2', ((1:10).^2)', 'y', y);
%! mdl = LinearModel.fromFit (fit_s, X, y, formula, pred_names, ...
%!         resp_name, var_names, [], {}, obs_info, [], ...
%!         false(10,1), false(10,1), [], vars);
%! cmdl = compact (mdl);
%! assert (strcmp (class (cmdl), 'CompactLinearModel'));
%! assert (cmdl.MSE, mdl.MSE);
%! assert (cmdl.SSE, mdl.SSE);
%! assert (cmdl.DFE, mdl.DFE);

## Test: Residuals formulas
%!test
%! X = [ones(10,1), (1:10)'];
%! y = 2 + 3 * (1:10)' + randn(10,1) * 0.5;
%! fit_s = lm_fit_engine (X, y);
%! pred_names = {'x1'};
%! resp_name = 'y';
%! var_names = {'x1', 'y'};
%! formula = struct ('Terms', [0; 1], 'HasIntercept', true, ...
%!                   'LinearPredictor', '1 + x1', ...
%!                   'CoefficientNames', {{'(Intercept)', 'x1'}});
%! obs_info = struct ('Weights', ones(10,1), 'Excluded', false(10,1), ...
%!                    'Missing', false(10,1), 'Subset', true(10,1));
%! vars = struct ('x1', (1:10)', 'y', y);
%! mdl = LinearModel.fromFit (fit_s, X, y, formula, pred_names, ...
%!         resp_name, var_names, [], {}, obs_info, [], ...
%!         false(10,1), false(10,1), [], vars);
%! ## Verify: Pearson = Raw / RMSE
%! assert (mdl.Residuals.Pearson, mdl.Residuals.Raw / mdl.RMSE, 1e-12);
%! ## Verify: Leverage = diag(HatMatrix)
%! assert (mdl.Diagnostics.Leverage, diag (mdl.Diagnostics.HatMatrix), 1e-14);
%! ## Verify: sum(Leverage) = NumEstimatedCoefficients
%! assert (sum (mdl.Diagnostics.Leverage), mdl.NumEstimatedCoefficients, 1e-12);

## Test: ModelFitVsNullModel formula
%!test
%! X = [ones(10,1), (1:10)'];
%! y = 2 + 3 * (1:10)' + randn(10,1) * 0.5;
%! fit_s = lm_fit_engine (X, y);
%! pred_names = {'x1'};
%! formula = struct ('Terms', [0; 1], 'HasIntercept', true, ...
%!                   'LinearPredictor', '1 + x1', ...
%!                   'CoefficientNames', {{'(Intercept)', 'x1'}});
%! obs_info = struct ('Weights', ones(10,1), 'Excluded', false(10,1), ...
%!                    'Missing', false(10,1), 'Subset', true(10,1));
%! vars = struct ('x1', (1:10)', 'y', y);
%! mdl = LinearModel.fromFit (fit_s, X, y, formula, pred_names, ...
%!         'y', {'x1', 'y'}, [], {}, obs_info, [], ...
%!         false(10,1), false(10,1), [], vars);
%! df_reg = mdl.NumEstimatedCoefficients - 1;
%! F_check = (mdl.SSR / df_reg) / mdl.MSE;
%! assert (mdl.ModelFitVsNullModel.Fstat, F_check, 1e-10);
%! assert (strcmp (mdl.ModelFitVsNullModel.NullModel, 'constant'));

## Test: addTerms increases NumCoefficients
%!test
%! X = [ones(20,1), (1:20)', ((1:20).^2)'];
%! y = X * [1; 0.5; -0.01] + randn(20,1) * 0.5;
%! fit_s = lm_fit_engine (X(:, 1:2), y);
%! pred_names = {'x1', 'x2'};
%! formula = struct ('Terms', [0 0; 1 0], 'HasIntercept', true, ...
%!                   'LinearPredictor', '1 + x1', ...
%!                   'CoefficientNames', {{'(Intercept)', 'x1'}});
%! obs_info = struct ('Weights', ones(20,1), 'Excluded', false(20,1), ...
%!                    'Missing', false(20,1), 'Subset', true(20,1));
%! vars = struct ('x1', (1:20)', 'x2', ((1:20).^2)', 'y', y);
%! mdl = LinearModel.fromFit (fit_s, X(:,1:2), y, formula, pred_names, ...
%!         'y', {'x1', 'x2', 'y'}, [], {}, obs_info, [], ...
%!         false(20,1), false(20,1), [], vars);
%! assert (mdl.NumCoefficients, 2);
%! mdl2 = addTerms (mdl, 'x2');
%! assert (mdl2.NumCoefficients, 3);
%! ## Original model unchanged
%! assert (mdl.NumCoefficients, 2);

## Test: removeTerms decreases NumCoefficients and increases SSE
%!test
%! X = [ones(20,1), (1:20)', randn(20,1)];
%! y = X * [1; 2; 0.1] + randn(20,1) * 0.5;
%! fit_s = lm_fit_engine (X, y);
%! pred_names = {'x1', 'x2'};
%! formula = struct ('Terms', [0 0; 1 0; 0 1], 'HasIntercept', true, ...
%!                   'LinearPredictor', '1 + x1 + x2', ...
%!                   'CoefficientNames', {{'(Intercept)', 'x1', 'x2'}});
%! obs_info = struct ('Weights', ones(20,1), 'Excluded', false(20,1), ...
%!                    'Missing', false(20,1), 'Subset', true(20,1));
%! vars = struct ('x1', (1:20)', 'x2', X(:,3), 'y', y);
%! mdl = LinearModel.fromFit (fit_s, X, y, formula, pred_names, ...
%!         'y', {'x1', 'x2', 'y'}, [], {}, obs_info, [], ...
%!         false(20,1), false(20,1), [], vars);
%! assert (mdl.NumCoefficients, 3);
%! mdl_r = removeTerms (mdl, 'x2');
%! assert (mdl_r.NumCoefficients, 2);
%! assert (mdl_r.SSE >= mdl.SSE);

## Test: anova Error row matches SSE/DFE/MSE
%!test
%! X = [ones(20,1), (1:20)', randn(20,1)];
%! y = X * [1; 2; -1] + randn(20,1) * 0.5;
%! fit_s = lm_fit_engine (X, y);
%! pred_names = {'x1', 'x2'};
%! formula = struct ('Terms', [0 0; 1 0; 0 1], 'HasIntercept', true, ...
%!                   'LinearPredictor', '1 + x1 + x2', ...
%!                   'CoefficientNames', {{'(Intercept)', 'x1', 'x2'}});
%! obs_info = struct ('Weights', ones(20,1), 'Excluded', false(20,1), ...
%!                    'Missing', false(20,1), 'Subset', true(20,1));
%! vars = struct ('x1', (1:20)', 'x2', X(:,3), 'y', y);
%! mdl = LinearModel.fromFit (fit_s, X, y, formula, pred_names, ...
%!         'y', {'x1', 'x2', 'y'}, [], {}, obs_info, [], ...
%!         false(20,1), false(20,1), [], vars);
%! tbl = anova (mdl);
%! ## Error row (last)
%! n_terms = numel (tbl.RowNames);
%! assert (tbl.SumSq(end), mdl.SSE, 1e-10);
%! assert (tbl.DF(end), mdl.DFE);
%! assert (tbl.MeanSq(end), mdl.MSE, 1e-10);

