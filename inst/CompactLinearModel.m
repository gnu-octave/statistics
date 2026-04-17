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

classdef CompactLinearModel < handle
## -*- texinfo -*-
## @deftp {statistics} CompactLinearModel
##
## Compact linear regression model.
##
## A @code{CompactLinearModel} object is a compact representation of a linear
## regression model.  It contains the regression results (coefficient
## estimates, goodness-of-fit statistics, and model information) but does not
## contain the training data, residuals, or observation-level diagnostics.
##
## A @code{CompactLinearModel} object is created by calling @code{compact} on
## a @code{LinearModel} object.  Users cannot construct a
## @code{CompactLinearModel} directly.
##
## @strong{Properties}
##
## All properties are read-only.
##
## @multitable @columnfractions 0.30 0.70
## @headitem Property @tab Description
## @item @code{Coefficients} @tab Table of coefficient estimates, SE, tStat,
## and pValue.
## @item @code{CoefficientCovariance} @tab Variance-covariance matrix of the
## coefficient estimates.
## @item @code{CoefficientNames} @tab Cell array of coefficient names.
## @item @code{DFE} @tab Error degrees of freedom.
## @item @code{Formula} @tab Struct describing the model formula.
## @item @code{LogLikelihood} @tab Log-likelihood of the fitted model.
## @item @code{ModelCriterion} @tab Struct with AIC, AICc, BIC, CAIC.
## @item @code{MSE} @tab Mean squared error.
## @item @code{NumCoefficients} @tab Total number of coefficients.
## @item @code{NumEstimatedCoefficients} @tab Number of estimated (non-rank-
## deficient) coefficients.
## @item @code{NumObservations} @tab Number of observations used in fitting.
## @item @code{NumPredictors} @tab Number of predictor variables.
## @item @code{NumVariables} @tab Total number of variables.
## @item @code{PredictorNames} @tab Cell array of predictor variable names.
## @item @code{ResponseName} @tab Name of the response variable.
## @item @code{RMSE} @tab Root mean squared error.
## @item @code{Robust} @tab Robust fitting information (empty for OLS).
## @item @code{Rsquared} @tab Struct with Ordinary and Adjusted R-squared.
## @item @code{SSE} @tab Sum of squared errors.
## @item @code{SSR} @tab Regression sum of squares.
## @item @code{SST} @tab Total sum of squares.
## @item @code{VariableInfo} @tab Table of variable metadata.
## @item @code{VariableNames} @tab Cell array of all variable names.
## @end multitable
##
## @strong{Methods}
##
## @multitable @columnfractions 0.30 0.70
## @headitem Method @tab Description
## @item @code{predict} @tab Predict response for new data.
## @item @code{coefCI} @tab Confidence intervals for coefficients.
## @item @code{coefTest} @tab Linear hypothesis test on coefficients.
## @end multitable
##
## @seealso{fitlm}
## @end deftp

  properties (SetAccess = protected)

    ## Coefficient estimates
    Coefficients              = [];
    CoefficientNames          = {};
    CoefficientCovariance     = [];
    NumCoefficients           = 0;
    NumEstimatedCoefficients  = 0;

    ## Summary statistics
    DFE                       = 0;
    SSE                       = 0;
    SSR                       = 0;
    SST                       = 0;
    MSE                       = 0;
    RMSE                      = 0;
    Rsquared                  = [];
    LogLikelihood             = 0;
    ModelCriterion            = [];

    ## Fitting method
    Robust                    = [];

    ## Input data info (survives compact())
    Formula                   = [];
    NumObservations           = 0;
    NumPredictors             = 0;
    NumVariables              = 0;
    PredictorNames            = {};
    ResponseName              = '';
    VariableInfo              = [];
    VariableNames             = {};

  endproperties

  ## Constructor — Access = protected so only LinearModel (subclass) or
  ## internal factory methods can create CompactLinearModel objects.
  ## When Octave supports (Access = ?LinearModel), switch to that.
  methods (Access = protected)

    function obj = CompactLinearModel (s)
      ## CompactLinearModel  Constructor (internal use only).
      ##
      ## obj = CompactLinearModel (s)
      ##
      ## S is a struct containing all fields needed to populate the object.
      ## This constructor is called by LinearModel.compact() or by fitlm's
      ## internal builder.  Users should not call it directly.

      if (nargin < 1)
        ## Allow zero-arg construction for Octave classdef internals
        return;
      endif

      if (! isstruct (s))
        error ("CompactLinearModel: constructor argument must be a struct.");
      endif

      ## Coefficient results
      obj.Coefficients             = s.Coefficients;
      obj.CoefficientNames         = s.CoefficientNames;
      obj.CoefficientCovariance    = s.CoefficientCovariance;
      obj.NumCoefficients          = s.NumCoefficients;
      obj.NumEstimatedCoefficients = s.NumEstimatedCoefficients;

      ## Summary statistics
      obj.DFE            = s.DFE;
      obj.SSE            = s.SSE;
      obj.SSR            = s.SSR;
      obj.SST            = s.SST;
      obj.MSE            = s.MSE;
      obj.RMSE           = s.RMSE;
      obj.Rsquared       = s.Rsquared;
      obj.LogLikelihood  = s.LogLikelihood;
      obj.ModelCriterion = s.ModelCriterion;

      ## Fitting method
      obj.Robust = s.Robust;

      ## Formula and variable info
      obj.Formula         = s.Formula;
      obj.NumObservations  = s.NumObservations;
      obj.NumPredictors    = s.NumPredictors;
      obj.NumVariables     = s.NumVariables;
      obj.PredictorNames   = s.PredictorNames;
      obj.ResponseName     = s.ResponseName;
      obj.VariableInfo     = s.VariableInfo;
      obj.VariableNames    = s.VariableNames;

    endfunction

  endmethods

  ## Hidden methods: display, disp, subsref
  methods (Hidden)

    ## Custom display — called when variable name is typed at prompt
    function display (obj)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ("%s =\n", in_name);
      endif
      disp (obj);
    endfunction

    ## Custom disp — replicates MATLAB's LinearModel display format
    function disp (obj)

      fprintf ("\n");

      ## Show "Compact" prefix for CompactLinearModel, not for LinearModel
      if (strcmp (class (obj), "CompactLinearModel"))
        fprintf ("Compact linear regression model:\n");
      else
        fprintf ("Linear regression model:\n");
      endif

      ## Formula line: ResponseName ~ LinearPredictor
      lp_str = "";
      if (isa (obj.Formula, "LinearFormula"))
        lp_str = obj.Formula.LinearPredictor;
      elseif (isstruct (obj.Formula) && isfield (obj.Formula, "LinearPredictor"))
        lp_str = obj.Formula.LinearPredictor;
      endif
      if (! isempty (lp_str))
        if (! isempty (obj.ResponseName))
          fprintf ("    %s ~ %s\n", obj.ResponseName, lp_str);
        else
          fprintf ("    %s\n", lp_str);
        endif
      endif

      fprintf ("\n");
      fprintf ("Estimated Coefficients:\n");

      ## --- Build coefficient table display ---
      coef_names = obj.CoefficientNames;
      p = numel (coef_names);

      ## Extract numeric vectors from Coefficients
      if (isstruct (obj.Coefficients))
        est  = obj.Coefficients.Estimate;
        se   = obj.Coefficients.SE;
        ts   = obj.Coefficients.tStat;
        pv   = obj.Coefficients.pValue;
      elseif (isa (obj.Coefficients, "table"))
        est  = obj.Coefficients.Estimate;
        se   = obj.Coefficients.SE;
        ts   = obj.Coefficients.tStat;
        pv   = obj.Coefficients.pValue;
      else
        est = []; se = []; ts = []; pv = [];
      endif

      if (! isempty (est))
        ## Column headers
        hdr_fmt = "    %-20s %12s %12s %12s %12s\n";
        fprintf (hdr_fmt, "", "Estimate", "SE", "tStat", "pValue");

        ## Separator line
        sep_fmt = "    %-20s %12s %12s %12s %12s\n";
        fprintf (sep_fmt, "", "________", "________", "________", "__________");
        fprintf ("\n");

        ## Data rows
        for i = 1:p
          row_fmt = "    %-20s %12.6g %12.6g %12.6g %12.4g\n";
          fprintf (row_fmt, coef_names{i}, est(i), se(i), ts(i), pv(i));
        endfor
      endif

      fprintf ("\n");

      ## Summary footer
      fprintf ("Number of observations: %d, Error degrees of freedom: %d\n", ...
               obj.NumObservations, obj.DFE);
      fprintf ("Root Mean Squared Error: %.4g\n", obj.RMSE);

      if (isstruct (obj.Rsquared))
        fprintf ("R-squared: %.4g,  Adjusted R-Squared: %.4g\n", ...
                 obj.Rsquared.Ordinary, obj.Rsquared.Adjusted);
      endif

      ## F-statistic vs constant model
      ## F = (SSR / (NumEstimatedCoefficients - 1)) / MSE
      if (obj.NumEstimatedCoefficients > 1 && obj.MSE > 0)
        df1 = obj.NumEstimatedCoefficients - 1;
        F_val = (obj.SSR / df1) / obj.MSE;
        F_pval = 1 - fcdf (F_val, df1, obj.DFE);
        fprintf ("F-statistic vs. constant model: %.4g, p-value = %.4g\n", ...
                 F_val, F_pval);
      endif

      fprintf ("\n");

    endfunction

    ## Custom subsref — supports chained dot indexing for properties
    function varargout = subsref (obj, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case "()"
          error (strcat ("Invalid () indexing for referencing values", ...
                         " in a CompactLinearModel object."));
        case "{}"
          error (strcat ("Invalid {} indexing for referencing values", ...
                         " in a CompactLinearModel object."));
        case "."
          if (! ischar (s.subs))
            error (strcat ("CompactLinearModel.subsref: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          ## Check if it's a method call
          switch (s.subs)
            case "predict"
              if (! isempty (chain_s) && strcmp (chain_s(1).type, "()"))
                args = chain_s(1).subs;
                chain_s = chain_s(2:end);
                [varargout{1:nargout}] = predict (obj, args{:});
              else
                error ("CompactLinearModel: predict requires arguments.");
              endif
              if (! isempty (chain_s))
                [varargout{1:nargout}] = subsref (varargout{1}, chain_s);
              endif
              return;
            case "coefCI"
              if (! isempty (chain_s) && strcmp (chain_s(1).type, "()"))
                args = chain_s(1).subs;
                chain_s = chain_s(2:end);
                [varargout{1:nargout}] = coefCI (obj, args{:});
              else
                [varargout{1:nargout}] = coefCI (obj);
              endif
              if (! isempty (chain_s))
                [varargout{1:nargout}] = subsref (varargout{1}, chain_s);
              endif
              return;
            case "coefTest"
              if (! isempty (chain_s) && strcmp (chain_s(1).type, "()"))
                args = chain_s(1).subs;
                chain_s = chain_s(2:end);
                [varargout{1:nargout}] = coefTest (obj, args{:});
              else
                [varargout{1:nargout}] = coefTest (obj);
              endif
              if (! isempty (chain_s))
                [varargout{1:nargout}] = subsref (varargout{1}, chain_s);
              endif
              return;
            case "compact"
              [varargout{1:nargout}] = compact (obj);
              if (! isempty (chain_s))
                [varargout{1:nargout}] = subsref (varargout{1}, chain_s);
              endif
              return;
            case "addTerms"
              if (! isempty (chain_s) && strcmp (chain_s(1).type, "()"))
                args = chain_s(1).subs;
                chain_s = chain_s(2:end);
                [varargout{1:nargout}] = addTerms (obj, args{:});
              else
                error ("CompactLinearModel: addTerms requires arguments.");
              endif
              if (! isempty (chain_s))
                [varargout{1:nargout}] = subsref (varargout{1}, chain_s);
              endif
              return;
            case "removeTerms"
              if (! isempty (chain_s) && strcmp (chain_s(1).type, "()"))
                args = chain_s(1).subs;
                chain_s = chain_s(2:end);
                [varargout{1:nargout}] = removeTerms (obj, args{:});
              else
                error ("CompactLinearModel: removeTerms requires arguments.");
              endif
              if (! isempty (chain_s))
                [varargout{1:nargout}] = subsref (varargout{1}, chain_s);
              endif
              return;
            case "step"
              if (! isempty (chain_s) && strcmp (chain_s(1).type, "()"))
                args = chain_s(1).subs;
                chain_s = chain_s(2:end);
                [varargout{1:nargout}] = step (obj, args{:});
              else
                [varargout{1:nargout}] = step (obj);
              endif
              if (! isempty (chain_s))
                [varargout{1:nargout}] = subsref (varargout{1}, chain_s);
              endif
              return;
            case "anova"
              if (! isempty (chain_s) && strcmp (chain_s(1).type, "()"))
                args = chain_s(1).subs;
                chain_s = chain_s(2:end);
                [varargout{1:nargout}] = anova (obj, args{:});
              else
                [varargout{1:nargout}] = anova (obj);
              endif
              if (! isempty (chain_s))
                [varargout{1:nargout}] = subsref (varargout{1}, chain_s);
              endif
              return;
            case "disp"
              disp (obj);
              return;
          endswitch
          ## Property access
          try
            out = obj.(s.subs);
          catch
            error ("CompactLinearModel.subsref: unrecognized property: '%s'.", ...
                   s.subs);
          end_try_catch
      endswitch
      ## Chained references
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

    ## Prevent direct property assignment by users
    function obj = subsasgn (obj, s, val)
      if (numel (s) > 1)
        error (strcat ("CompactLinearModel.subsasgn:", ...
                       " chained subscripts not allowed."));
      endif
      switch (s.type)
        case "()"
          error (strcat ("Invalid () indexing for assigning values", ...
                         " to a CompactLinearModel object."));
        case "{}"
          error (strcat ("Invalid {} indexing for assigning values", ...
                         " to a CompactLinearModel object."));
        case "."
          error (strcat ("CompactLinearModel.subsasgn:", ...
                         " all properties are read-only: '%s'."), s.subs);
      endswitch
    endfunction

  endmethods

  ## Public methods: predict, coefCI, coefTest
  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactLinearModel} {@var{ypred} =} predict (@var{obj}, @var{Xnew})
    ## @deftypefnx {CompactLinearModel} {[@var{ypred}, @var{yci}] =} predict (@var{obj}, @var{Xnew})
    ## @deftypefnx {CompactLinearModel} {[@var{ypred}, @var{yci}] =} predict (@dots{}, @var{Name}, @var{Value})
    ##
    ## Predict response for new data using a compact linear regression model.
    ##
    ## @code{@var{ypred} = predict (@var{obj}, @var{Xnew})} returns the predicted
    ## response values for the predictor data in @var{Xnew} using the fitted
    ## model @var{obj}.  @var{Xnew} must be a numeric matrix with the same
    ## number of predictor columns as used to fit the model (excluding
    ## intercept — it is added automatically if the model has one).
    ##
    ## @code{[@var{ypred}, @var{yci}] = predict (@dots{})} also returns
    ## confidence or prediction intervals in the @math{n x 2} matrix @var{yci}.
    ##
    ## Additional options can be specified as @qcode{Name-Value} pairs:
    ##
    ## @table @code
    ## @item Alpha
    ## Significance level for the confidence/prediction intervals.
    ## Default is 0.05 (95% intervals).
    ##
    ## @item Prediction
    ## Type of interval: @qcode{'curve'} for confidence intervals on the
    ## mean response (default), or @qcode{'observation'} for prediction
    ## intervals on individual new observations.
    ##
    ## @item Simultaneous
    ## If @code{true}, compute simultaneous Scheffé confidence bands.
    ## Default is @code{false} (pointwise intervals).
    ## @end table
    ##
    ## @seealso{CompactLinearModel, coefCI, coefTest}
    ## @end deftypefn

    function [ypred, yci] = predict (obj, Xnew, varargin)

      if (nargin < 2)
        error ("CompactLinearModel.predict: Xnew is required.");
      endif

      ## Parse name-value pairs
      alpha    = 0.05;
      pred_obs = false;   ## false = curve (CI on mean), true = observation (PI)
      simult   = false;

      i = 1;
      while (i <= numel (varargin))
        if (! ischar (varargin{i}))
          error ("CompactLinearModel.predict: expected parameter name string.");
        endif
        switch (lower (varargin{i}))
          case "alpha"
            alpha = varargin{i+1};
            i = i + 2;
          case "prediction"
            pred_obs = strcmpi (varargin{i+1}, "observation");
            i = i + 2;
          case "simultaneous"
            simult = varargin{i+1};
            i = i + 2;
          otherwise
            error ("CompactLinearModel.predict: unknown parameter '%s'.", ...
                   varargin{i});
        endswitch
      endwhile

      ## Build design matrix: prepend intercept column if model has one
      has_int = false;
      if (isa (obj.Formula, "LinearFormula"))
        has_int = obj.Formula.HasIntercept;
      elseif (isstruct (obj.Formula) && isfield (obj.Formula, "HasIntercept"))
        has_int = obj.Formula.HasIntercept;
      endif
      if (has_int)
        Xmat = [ones(rows (Xnew), 1), Xnew];
      else
        Xmat = Xnew;
      endif

      ## Check dimensions
      if (columns (Xmat) != obj.NumCoefficients)
        error (strcat ("CompactLinearModel.predict: Xnew has wrong", ...
                       " number of columns. Expected %d predictor columns."), ...
               obj.NumCoefficients - double (obj.Formula.HasIntercept));
      endif

      ## Extract coefficient estimates
      if (isa (obj.Coefficients, "table"))
        beta = obj.Coefficients.Estimate;
      elseif (isstruct (obj.Coefficients))
        beta = obj.Coefficients.Estimate;
      else
        error ("CompactLinearModel.predict: invalid Coefficients format.");
      endif

      ## Point predictions
      ypred = Xmat * beta;

      ## Compute intervals if requested
      if (nargout > 1)
        se_fit = sqrt (sum ((Xmat * obj.CoefficientCovariance) .* Xmat, 2));

        if (pred_obs)
          se_pred = sqrt (se_fit .^ 2 + obj.MSE);
        else
          se_pred = se_fit;
        endif

        if (simult)
          p_eff = obj.NumEstimatedCoefficients;
          F_crit = finv (1 - alpha, p_eff, obj.DFE);
          hw = sqrt (p_eff * F_crit) .* se_pred;
        else
          hw = tinv (1 - alpha / 2, obj.DFE) .* se_pred;
        endif

        yci = [ypred - hw, ypred + hw];
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactLinearModel} {@var{ci} =} coefCI (@var{obj})
    ## @deftypefnx {CompactLinearModel} {@var{ci} =} coefCI (@var{obj}, @var{alpha})
    ##
    ## Confidence intervals for coefficient estimates.
    ##
    ## @code{@var{ci} = coefCI (@var{obj})} returns a @math{p x 2} matrix
    ## of 95% confidence intervals for the coefficient estimates, where
    ## @math{p} is the number of coefficients.
    ##
    ## @code{@var{ci} = coefCI (@var{obj}, @var{alpha})} specifies the
    ## significance level.  For example, @code{alpha = 0.01} gives 99%
    ## confidence intervals.
    ##
    ## @seealso{CompactLinearModel, predict, coefTest}
    ## @end deftypefn

    function ci = coefCI (obj, alpha)

      if (nargin < 2)
        alpha = 0.05;
      endif

      t_crit = tinv (1 - alpha / 2, obj.DFE);

      ## Extract estimates and SEs
      if (isa (obj.Coefficients, "table"))
        est = obj.Coefficients.Estimate;
        se  = obj.Coefficients.SE;
      elseif (isstruct (obj.Coefficients))
        est = obj.Coefficients.Estimate;
        se  = obj.Coefficients.SE;
      else
        error ("CompactLinearModel.coefCI: invalid Coefficients format.");
      endif

      ci = [est - t_crit * se, est + t_crit * se];

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactLinearModel} {[@var{pval}, @var{F}, @var{r}] =} coefTest (@var{obj})
    ## @deftypefnx {CompactLinearModel} {[@var{pval}, @var{F}, @var{r}] =} coefTest (@var{obj}, @var{H})
    ## @deftypefnx {CompactLinearModel} {[@var{pval}, @var{F}, @var{r}] =} coefTest (@var{obj}, @var{H}, @var{c})
    ##
    ## Linear hypothesis test on coefficients.
    ##
    ## @code{[@var{pval}, @var{F}, @var{r}] = coefTest (@var{obj})} performs
    ## the overall F-test: tests whether all non-intercept coefficients are
    ## simultaneously zero.  The F-statistic matches the "F-statistic vs.
    ## constant model" shown by @code{disp}.
    ##
    ## @code{[@var{pval}, @var{F}, @var{r}] = coefTest (@var{obj}, @var{H})}
    ## tests the hypothesis @code{@var{H} * beta = 0}, where @var{H} is a
    ## contrast matrix.
    ##
    ## @code{[@var{pval}, @var{F}, @var{r}] = coefTest (@var{obj}, @var{H},
    ## @var{c})} tests the hypothesis @code{@var{H} * beta = @var{c}}.
    ##
    ## The output @var{r} is the numerator degrees of freedom for the
    ## F-test (equal to the number of rows in @var{H}).  The F-statistic
    ## has @var{r} degrees of freedom in the numerator and @code{obj.DFE}
    ## degrees of freedom in the denominator.
    ##
    ## @seealso{CompactLinearModel, predict, coefCI}
    ## @end deftypefn

    function [pval, F, r] = coefTest (obj, H, c)

      ## Extract beta
      if (isa (obj.Coefficients, "table"))
        beta = obj.Coefficients.Estimate;
      elseif (isstruct (obj.Coefficients))
        beta = obj.Coefficients.Estimate;
      else
        error ("CompactLinearModel.coefTest: invalid Coefficients format.");
      endif

      p = obj.NumCoefficients;

      ## Default H: overall F-test (all non-intercept = 0)
      if (nargin < 2 || isempty (H))
        ## H = [zeros(p-1, 1), eye(p-1)]  — tests all columns except intercept
        H = [zeros(p - 1, 1), eye(p - 1)];
      endif

      ## Default c: zero vector
      if (nargin < 3 || isempty (c))
        c = zeros (rows (H), 1);
      endif

      ## Contrast vector (used internally for F computation)
      contrast = H * beta - c;

      ## Wald F-test
      V   = H * obj.CoefficientCovariance * H';
      df1 = rows (H);
      F   = (contrast' * (V \ contrast)) / df1;
      pval = 1 - fcdf (F, df1, obj.DFE);

      ## r = numerator degrees of freedom (= rows of H)
      r = df1;

    endfunction

  endmethods

endclassdef

## Test: CompactLinearModel cannot be constructed directly (protected access)
## This should error since the constructor is protected.
%!error <CompactLinearModel> CompactLinearModel (struct ())

## Test: verify class definition parses without error
## (This test just confirms the file is syntactically valid)
%!test
%! ## Check that the class exists and is loadable
%! assert (exist ("CompactLinearModel", "class") > 0 || ...
%!         exist ("CompactLinearModel") > 0);
