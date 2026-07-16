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
## @deftypefn {statistics} {@var{mdl} =} NonLinearModel (@dots{})
##
## Nonlinear regression model class.
##
## A @code{NonLinearModel} object holds a nonlinear regression fitted by
## @code{fitnlm}, together with its coefficients, fit statistics, and methods
## for inference, prediction, and diagnostics.  Construct one with
## @code{fitnlm}, which documents the accepted inputs and
## @var{Name}/@var{Value} pairs.
##
## The estimated coefficients and their statistics are in the
## @code{Coefficients} table; @code{Rsquared}, @code{ModelCriterion},
## @code{LogLikelihood}, @code{RMSE}, @code{SSE}, @code{SST}, and @code{SSR}
## summarize the fit.  The methods @code{predict}, @code{feval}, @code{random},
## @code{coefCI}, @code{coefTest}, @code{plotResiduals}, @code{plotDiagnostics},
## and @code{plotSlice} operate on the fitted model.
##
## @subheading Fit statistics
##
## The fit statistics follow MATLAB's conventions.  @code{SSE} is the residual
## sum of squares, @code{SST} the total sum of squares of the response about its
## (weighted) mean, and @code{SSR} the regression sum of squares of the fitted
## values about that mean; because the model is nonlinear, @code{SST} does
## @emph{not} in general equal @code{SSR + SSE}.  @code{Rsquared.Ordinary} is
## @code{1 - @var{SSE} / @var{SST}} and @code{Rsquared.Adjusted} corrects for
## the error degrees of freedom.  @code{RMSE} is @code{sqrt (@var{MSE})}, and
## the Gaussian @code{LogLikelihood} uses the maximum-likelihood error variance
## @code{@var{SSE} / n}.  The information criteria in @code{ModelCriterion}
## (@code{AIC}, @code{AICc}, @code{BIC}, @code{CAIC}) count the @math{p}
## coefficients as the only parameters -- the error variance is @emph{not}
## counted.  @code{coefTest} is a Wald test: for a contrast matrix @var{H} it
## forms @code{(@var{H}*b)' * inv (@var{H}*@var{V}*@var{H}') * (@var{H}*b) / r}
## with @var{V} the coefficient covariance and @math{r} the number of rows of
## @var{H}, referred to an @math{F} distribution on @math{r} and @var{DFE}
## degrees of freedom.  The summary printed by @code{disp} instead reports an
## @math{F} statistic versus the zero model, formed from the uncorrected
## regression sum of squares (the sum of the squared fitted values).
##
## @seealso{fitnlm, nlinfit, nlparci, nlpredci, LinearModel,
## GeneralizedLinearModel}
## @end deftypefn

classdef NonLinearModel

  properties (GetAccess = public, SetAccess = protected)
    Coefficients = [];
    CoefficientNames = {};
    CoefficientCovariance = [];
    NumCoefficients = [];
    NumEstimatedCoefficients = [];
    NumPredictors = [];
    NumObservations = [];
    DFE = [];
    MSE = [];
    RMSE = [];
    SSE = [];
    SST = [];
    SSR = [];
    LogLikelihood = [];
    ModelCriterion = [];
    Rsquared = [];
    Residuals = [];
    Fitted = [];
    ErrorModelInfo = [];
    Robust = [];
    Formula = '';
    ResponseName = 'y';
    PredictorNames = {};
    VariableNames = {};
  endproperties

  properties (Access = private, Hidden)
    modelfun_  = [];   # model function handle @(b, X)
    beta_      = [];   # fitted coefficient vector
    X_         = [];   # numeric predictor matrix used for the fit
    y_         = [];   # response vector used for the fit
    w_         = [];   # observation weights (effective)
    J_         = [];   # Jacobian at the solution
    R_         = [];   # raw residuals at the solution
    leverage_  = [];   # hat-matrix diagonal
    istable_   = false;   # true when the model was fit from a table
  endproperties

  methods (Access = public)

    ## Custom display of the object with its variable name.
    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ("%s =\n", in_name);
      endif
      disp (this);
    endfunction

    ## Custom display of the model summary.
    function disp (this)
      fprintf ("\n  Nonlinear regression model:\n");
      if (! isempty (this.Formula))
        fprintf ("      %s\n", this.Formula);
      endif
      if (! isempty (this.Coefficients))
        fprintf ("\n  Estimated Coefficients:\n\n");
        disp (this.Coefficients);
      endif
      fprintf ("\n");
      if (! isempty (this.NumObservations) && ! isempty (this.DFE))
        fprintf (strcat ("Number of observations: %d,", ...
                         " Error degrees of freedom: %d\n"), ...
                 this.NumObservations, this.DFE);
      endif
      if (! isempty (this.RMSE))
        fprintf ("Root Mean Squared Error: %.3g\n", this.RMSE);
      endif
      if (! isempty (this.Rsquared))
        fprintf ("R-Squared: %.3g,  Adjusted R-Squared %.3g\n", ...
                 this.Rsquared.Ordinary, this.Rsquared.Adjusted);
      endif
      if (! isempty (this.beta_) && ! isempty (this.MSE))
        ## F versus the zero model uses the uncorrected regression sum of
        ## squares (a nonlinear model has no guaranteed intercept term).
        yhat  = this.Fitted;
        Fstat = (sum (yhat .^ 2) / this.NumCoefficients) / this.MSE;
        pval  = 1 - fcdf (Fstat, this.NumCoefficients, this.DFE);
        fprintf ("F-statistic vs. zero model: %.3g, p-value = %.3g\n", ...
                 Fstat, pval);
      endif
    endfunction

    ## Class specific subscripted reference.
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("NonLinearModel: () indexing is not supported.", ...
                         "  Use dot notation for properties."));
        case '{}'
          error (strcat ("NonLinearModel: {} indexing is not supported.", ...
                         "  Use dot notation for properties."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("NonLinearModel.subsref: property name must be", ...
                           " a character vector."));
          endif
          if (ismethod (this, s.subs))
            [varargout{1:nargout}] = builtin ('subsref', this, [s, chain_s]);
            return;
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("NonLinearModel.subsref: unknown property", ...
                           " '%s'."), s.subs);
          end_try_catch
      endswitch
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {NonLinearModel} {@var{mdl} =} NonLinearModel (@var{data}, @var{resp}, @var{modelfun}, @var{beta0})
    ## @deftypefnx {NonLinearModel} {@var{mdl} =} NonLinearModel (@dots{}, @var{Name}, @var{Value})
    ##
    ## Fit a nonlinear regression model.  Prefer the @code{fitnlm} function,
    ## which documents the accepted inputs and @var{Name}/@var{Value} pairs.
    ##
    ## @end deftypefn
    function this = NonLinearModel (data, resp, modelfun, beta0, varargin)

      if (nargin == 0)
        return;   # empty object
      endif
      if (nargin < 4)
        error ("NonLinearModel: DATA, RESP, MODELFUN, and BETA0 are required.");
      endif
      if (! is_function_handle (modelfun))
        error ("NonLinearModel: MODELFUN must be a function handle.");
      endif
      if (! (isnumeric (beta0) && isvector (beta0) && isreal (beta0)))
        error ("NonLinearModel: BETA0 must be a real numeric vector.");
      endif

      opts = nlm_parse_nv (varargin);
      beta0 = beta0(:);
      p = numel (beta0);

      ## ---------------------------------------------------------------- ##
      ## Intake: resolve the predictor matrix, response, and names.
      ## ---------------------------------------------------------------- ##
      if (istable (data))
        this.istable_ = true;
        col_names = data.Properties.VariableNames;
        if (! isempty (opts.ResponseVar))
          resp_name = opts.ResponseVar;
        else
          resp_name = col_names{end};
        endif
        if (! isempty (opts.PredictorVars))
          pred_names = opts.PredictorVars;
        else
          pred_names = col_names(! strcmp (col_names, resp_name));
        endif
        y_full = double (data.(resp_name)(:));
        X_full = zeros (numel (y_full), numel (pred_names));
        for j = 1:numel (pred_names)
          X_full(:,j) = double (data.(pred_names{j})(:));
        endfor
      else
        if (! (isnumeric (data) && isreal (data) && ismatrix (data)))
          error ("NonLinearModel: X must be a real matrix.");
        endif
        if (! (isnumeric (resp) && isreal (resp) && isvector (resp)))
          error ("NonLinearModel: Y must be a real vector.");
        endif
        X_full = double (data);
        y_full = double (resp(:));
        if (rows (X_full) != numel (y_full))
          error (strcat ("NonLinearModel: X and Y must have the same", ...
                         " number of observations."));
        endif
        p_raw = columns (X_full);
        if (! isempty (opts.VarNames))
          if (numel (opts.VarNames) != p_raw + 1)
            error ("NonLinearModel: VarNames must have %d elements.", ...
                   p_raw + 1);
          endif
          pred_names = opts.VarNames(1:p_raw)(:)';
          resp_name  = opts.VarNames{end};
        else
          pred_names = arrayfun (@(k) sprintf ("x%d", k), 1:p_raw, ...
                                 'UniformOutput', false);
          resp_name  = 'y';
        endif
        if (! isempty (opts.ResponseVar))
          resp_name = opts.ResponseVar;
        endif
      endif

      ## Exclude and drop missing rows.
      n_total = numel (y_full);
      excluded = false (n_total, 1);
      if (! isempty (opts.Exclude))
        ex = opts.Exclude(:);
        if (islogical (ex))
          excluded(1:numel (ex)) = ex;
        else
          excluded(ex) = true;
        endif
      endif
      missing = any (isnan (X_full), 2) | isnan (y_full);
      keep = ! excluded & ! missing;
      X = X_full(keep,:);
      y = y_full(keep);
      n = numel (y);
      if (n <= p)
        error (strcat ("NonLinearModel: not enough observations to fit", ...
                       " %d coefficients."), p);
      endif

      ## Coefficient names.
      if (! isempty (opts.CoefficientNames))
        if (numel (opts.CoefficientNames) != p)
          error ("NonLinearModel: CoefficientNames must have %d elements.", p);
        endif
        coef_names = opts.CoefficientNames(:)';
      else
        coef_names = arrayfun (@(k) sprintf ("b%d", k), 1:p, ...
                               'UniformOutput', false);
      endif

      ## ---------------------------------------------------------------- ##
      ## Fit via nlinfit, forwarding the fitting options.
      ## ---------------------------------------------------------------- ##
      nvfit = {};
      if (! isempty (opts.Weights))
        nvfit(end+1:end+2) = {"Weights", opts.Weights(keep)};
      endif
      if (! isempty (opts.ErrorModel))
        nvfit(end+1:end+2) = {"ErrorModel", opts.ErrorModel};
      endif
      if (! isempty (opts.RobustWgtFun))
        nvfit(end+1:end+2) = {"RobustWgtFun", opts.RobustWgtFun};
      endif
      if (! isempty (opts.Tune))
        nvfit(end+1:end+2) = {"Tune", opts.Tune};
      endif
      if (! isempty (opts.Options))
        nvfit(end+1:end+2) = {"Options", opts.Options};
      endif

      [beta, R, J, CovB, MSE, EMI] = nlinfit (X, y, modelfun, beta0, nvfit{:});

      ## ---------------------------------------------------------------- ##
      ## Populate the fit statistics.
      ## ---------------------------------------------------------------- ##
      dfe  = n - p;
      yhat = y - R;
      if (isempty (opts.Weights))
        w = ones (n, 1);
      else
        w = double (opts.Weights(keep)(:));
      endif

      SSE = sum (w .* R .^ 2);
      ybar = sum (w .* y) / sum (w);
      SST = sum (w .* (y - ybar) .^ 2);
      SSR = sum (w .* (yhat - ybar) .^ 2);
      se  = sqrt (diag (CovB));
      tstat = beta ./ se;
      pval  = 2 * tcdf (-abs (tstat), dfe);

      ## Hat-matrix diagonal from the weighted Jacobian.
      Jw = sqrt (w) .* J;
      H  = Jw * pinv (Jw' * Jw) * Jw';
      lev = diag (H);

      LL  = -0.5 * n * (log (2 * pi) + log (SSE / n) + 1);
      AIC  = -2 * LL + 2 * p;
      AICc = AIC + 2 * p * (p + 1) / max (n - p - 1, 1);
      BIC  = -2 * LL + p * log (n);
      CAIC = -2 * LL + p * (log (n) + 1);

      ## Residuals table (raw, Pearson, standardized, studentized).
      rmse = sqrt (MSE);
      raw  = R;
      pear = R .* sqrt (w);
      stnd = R ./ (rmse .* sqrt (max (1 - lev, eps)));
      s_i  = sqrt (max ((dfe * MSE - (R .^ 2 .* w) ./ max (1 - lev, eps)) ...
                        / max (dfe - 1, 1), 0));
      stud = R ./ (s_i .* sqrt (max (1 - lev, eps)) + eps);

      ## Assemble the properties.
      this.modelfun_ = modelfun;
      this.beta_     = beta;
      this.X_        = X;
      this.y_        = y;
      this.w_        = w;
      this.J_        = J;
      this.R_        = R;
      this.leverage_ = lev;

      this.Coefficients = table (beta, se, tstat, pval, ...
        'VariableNames', {'Estimate', 'SE', 'tStat', 'pValue'}, ...
        'RowNames', coef_names);
      this.CoefficientNames = coef_names;
      this.CoefficientCovariance = CovB;
      this.NumCoefficients = p;
      this.NumEstimatedCoefficients = p;
      this.NumPredictors = columns (X);
      this.NumObservations = n;
      this.DFE = dfe;
      this.MSE = MSE;
      this.RMSE = rmse;
      this.SSE = SSE;
      this.SST = SST;
      this.SSR = SSR;
      this.LogLikelihood = LL;
      this.ModelCriterion = struct ('AIC', AIC, 'AICc', AICc, ...
                                    'BIC', BIC, 'CAIC', CAIC);
      this.Rsquared = struct ('Ordinary', 1 - SSE / SST, ...
                              'Adjusted', 1 - (SSE / dfe) / (SST / (n - 1)));
      this.Residuals = table (raw, pear, stnd, stud, 'VariableNames', ...
        {'Raw', 'Pearson', 'Standardized', 'Studentized'});
      this.Fitted = yhat;
      this.ErrorModelInfo = EMI;
      this.Robust = ternary (isempty (opts.RobustWgtFun), [], ...
                             struct ('RobustWgtFun', opts.RobustWgtFun));
      this.ResponseName = resp_name;
      this.PredictorNames = pred_names;
      this.VariableNames = [pred_names, {resp_name}];
      this.Formula = build_formula_string (modelfun, coef_names, resp_name);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {NonLinearModel} {@var{yhat} =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {NonLinearModel} {[@var{yhat}, @var{yci}] =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {NonLinearModel} {[@dots{}] =} predict (@dots{}, @var{Name}, @var{Value})
    ##
    ## Predict responses of the nonlinear model @var{mdl} at the new predictor
    ## values @var{Xnew} (a numeric matrix or a table).  With two outputs it
    ## also returns the confidence intervals @var{yci}.  Accepts @qcode{'Alpha'}
    ## (default 0.05), @qcode{'Prediction'} (@qcode{'curve'} or
    ## @qcode{'observation'}), and @qcode{'Simultaneous'} (a logical).
    ##
    ## @end deftypefn
    function [yhat, yci] = predict (mdl, Xnew, varargin)
      alpha = 0.05;  predtype = 'curve';  simul = false;
      for k = 1:2:numel (varargin)
        switch (lower (varargin{k}))
          case 'alpha'
            alpha = varargin{k+1};
          case 'prediction'
            predtype = lower (varargin{k+1});
          case 'simultaneous'
            simul = varargin{k+1};
          otherwise
            error ("NonLinearModel.predict: unknown parameter '%s'.", ...
                   varargin{k});
        endswitch
      endfor
      Xq = mdl.resolve_predictors (Xnew);
      yhat = mdl.modelfun_ (mdl.beta_, Xq);
      yhat = yhat(:);
      if (nargout > 1)
        simopt = 'off';
        if (simul)
          simopt = 'on';
        endif
        [~, delta] = nlpredci (mdl.modelfun_, Xq, mdl.beta_, mdl.R_, ...
                               'Jacobian', mdl.J_, 'MSE', mdl.MSE, ...
                               'PredOpt', predtype, 'SimOpt', simopt, ...
                               'Alpha', alpha);
        yci = [yhat - delta, yhat + delta];
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {NonLinearModel} {@var{yhat} =} feval (@var{mdl}, @var{X})
    ##
    ## Evaluate the fitted model at the predictor values @var{X}, given either
    ## as a single matrix/table or as separate column arguments (one per
    ## predictor).
    ##
    ## @end deftypefn
    function yhat = feval (mdl, varargin)
      if (numel (varargin) == 1)
        Xq = mdl.resolve_predictors (varargin{1});
      else
        Xq = cell2mat (cellfun (@(v) v(:), varargin, 'UniformOutput', false));
      endif
      yhat = mdl.modelfun_ (mdl.beta_, Xq);
      yhat = yhat(:);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {NonLinearModel} {@var{ysim} =} random (@var{mdl}, @var{Xnew})
    ##
    ## Simulate responses from the fitted model at @var{Xnew} (default: the
    ## training predictors), adding Gaussian noise with the model's error
    ## standard deviation.
    ##
    ## @end deftypefn
    function ysim = random (mdl, Xnew)
      if (nargin < 2)
        Xq = mdl.X_;
      else
        Xq = mdl.resolve_predictors (Xnew);
      endif
      yhat = mdl.modelfun_ (mdl.beta_, Xq);
      yhat = yhat(:);
      ysim = yhat + sqrt (mdl.MSE) * randn (numel (yhat), 1);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {NonLinearModel} {@var{ci} =} coefCI (@var{mdl})
    ## @deftypefnx {NonLinearModel} {@var{ci} =} coefCI (@var{mdl}, @var{alpha})
    ##
    ## Confidence intervals for the coefficients at level
    ## @math{100 (1 - @var{alpha})%} (default @var{alpha} = 0.05).
    ##
    ## @end deftypefn
    function ci = coefCI (mdl, alpha)
      if (nargin < 2)
        alpha = 0.05;
      endif
      se = mdl.Coefficients.SE;
      b  = mdl.Coefficients.Estimate;
      t  = tinv (1 - alpha / 2, mdl.DFE);
      ci = [b - t .* se, b + t .* se];
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {NonLinearModel} {[@var{p}, @var{F}, @var{df}] =} coefTest (@var{mdl})
    ## @deftypefnx {NonLinearModel} {[@dots{}] =} coefTest (@var{mdl}, @var{H})
    ##
    ## Wald test of a linear hypothesis on the coefficients.  With no @var{H} it
    ## tests that all coefficients are zero (the model versus the zero model)
    ## and returns the @math{p}-value @var{p}, the @math{F} statistic @var{F},
    ## and its numerator degrees of freedom @var{df}.  @var{H} is an
    ## @math{r}-by-@math{p} contrast matrix testing @code{@var{H} * @var{beta}
    ## = 0}.
    ##
    ## @end deftypefn
    function [p, F, df] = coefTest (mdl, H)
      b = mdl.beta_;
      if (nargin < 2)
        H = eye (numel (b));
      endif
      df = rows (H);
      M  = H * mdl.CoefficientCovariance * H';
      F  = ((H * b)' * pinv (M) * (H * b)) / df;
      p  = 1 - fcdf (F, df, mdl.DFE);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {NonLinearModel} {@var{h} =} plotResiduals (@var{mdl})
    ## @deftypefnx {NonLinearModel} {@var{h} =} plotResiduals (@var{mdl}, @var{plottype})
    ##
    ## Plot the model residuals.  @var{plottype} is @qcode{'histogram'}
    ## (default),
    ## @qcode{'fitted'}, @qcode{'caseorder'}, or @qcode{'probability'}.
    ##
    ## @end deftypefn
    function h = plotResiduals (mdl, plottype, varargin)
      if (nargin < 2)
        plottype = 'histogram';
      endif
      r = mdl.Residuals.Raw;
      switch (lower (plottype))
        case 'histogram'
          h = hist (r);
          xlabel ("Residuals");  ylabel ("Frequency");
        case 'fitted'
          h = plot (mdl.Fitted, r, 'o');
          xlabel ("Fitted values");  ylabel ("Residuals");
          hold on;  plot (xlim (), [0, 0], 'k:');  hold off;
        case 'caseorder'
          h = plot (1:numel (r), r, 'o-');
          xlabel ("Case order");  ylabel ("Residuals");
        case 'probability'
          h = plot (sort (r), ((1:numel (r)) - 0.5) / numel (r), 'o');
          xlabel ("Residuals");  ylabel ("Probability");
        otherwise
          error ("NonLinearModel.plotResiduals: unknown plot type '%s'.", ...
                 plottype);
      endswitch
      title ("Residuals");
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {NonLinearModel} {@var{h} =} plotDiagnostics (@var{mdl})
    ## @deftypefnx {NonLinearModel} {@var{h} =} plotDiagnostics (@var{mdl}, @var{plottype})
    ##
    ## Plot fit diagnostics.  @var{plottype} is @qcode{'leverage'} (default) or
    ## @qcode{'cookd'} (Cook's distance).
    ##
    ## @end deftypefn
    function h = plotDiagnostics (mdl, plottype)
      if (nargin < 2)
        plottype = 'leverage';
      endif
      switch (lower (plottype))
        case 'leverage'
          h = stem (mdl.leverage_);
          xlabel ("Observation");  ylabel ("Leverage");
        case 'cookd'
          lev = mdl.leverage_;
          cd = (mdl.R_ .^ 2 ./ (mdl.NumCoefficients * mdl.MSE)) ...
               .* (lev ./ max (1 - lev, eps) .^ 2);
          h = stem (cd);
          xlabel ("Observation");  ylabel ("Cook's distance");
        otherwise
          error ("NonLinearModel.plotDiagnostics: unknown plot type '%s'.", ...
                 plottype);
      endswitch
      title ("Diagnostics");
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {NonLinearModel} {@var{h} =} plotSlice (@var{mdl})
    ##
    ## Plot the fitted response as each predictor is varied over its observed
    ## range with the others held at their means.
    ##
    ## @end deftypefn
    function h = plotSlice (mdl)
      pn = mdl.NumPredictors;
      xm = mean (mdl.X_, 1);
      h = zeros (pn, 1);
      for j = 1:pn
        subplot (1, pn, j);
        xj = linspace (min (mdl.X_(:,j)), max (mdl.X_(:,j)), 50)';
        Xg = repmat (xm, numel (xj), 1);
        Xg(:,j) = xj;
        yg = mdl.modelfun_ (mdl.beta_, Xg);
        h(j) = plot (xj, yg(:));
        xlabel (mdl.PredictorNames{j});  ylabel (mdl.ResponseName);
      endfor
    endfunction

  endmethods

  methods (Access = private)

    ## Resolve new predictor values from a matrix or a table to a numeric matrix
    ## matching the columns used at fit time.
    function Xq = resolve_predictors (mdl, Xnew)
      if (istable (Xnew))
        Xq = zeros (height (Xnew), numel (mdl.PredictorNames));
        for j = 1:numel (mdl.PredictorNames)
          Xq(:,j) = double (Xnew.(mdl.PredictorNames{j})(:));
        endfor
      else
        Xq = double (Xnew);
      endif
    endfunction

  endmethods

endclassdef

## ---------------------------------------------------------------------------
## Parse the fitnlm/NonLinearModel Name/Value options.
function opts = nlm_parse_nv (nv)

  opts = struct ("CoefficientNames", [], "Weights", [], "ErrorModel", [], ...
                 "RobustWgtFun", [], "Tune", [], "Options", [], ...
                 "PredictorVars", [], "ResponseVar", [], "VarNames", [], ...
                 "Exclude", []);
  if (mod (numel (nv), 2) != 0)
    error ("NonLinearModel: Name/Value arguments must come in pairs.");
  endif
  for k = 1:2:numel (nv)
    name = nv{k};
    if (! ischar (name))
      error ("NonLinearModel: parameter names must be character vectors.");
    endif
    switch (lower (name))
      case 'coefficientnames'
        opts.CoefficientNames = nv{k+1};
      case 'weights'
        opts.Weights = nv{k+1};
      case 'errormodel'
        opts.ErrorModel = nv{k+1};
      case 'robustwgtfun'
        opts.RobustWgtFun = nv{k+1};
      case 'tune'
        opts.Tune = nv{k+1};
      case 'options'
        opts.Options = nv{k+1};
      case 'predictorvars'
        opts.PredictorVars = nv{k+1};
      case 'responsevar'
        opts.ResponseVar = nv{k+1};
      case 'varnames'
        opts.VarNames = nv{k+1};
      case 'exclude'
        opts.Exclude = nv{k+1};
      otherwise
        error ("NonLinearModel: unknown parameter name '%s'.", name);
    endswitch
  endfor

endfunction

## ---------------------------------------------------------------------------
## Build the display formula "resp ~ body" from an anonymous model function,
## substituting b(k) -> coefname and stripping element-wise dots.
function s = build_formula_string (modelfun, coef_names, resp_name)

  fstr = func2str (modelfun);
  ## Extract the argument names and the body of "@(b, x) body".
  tok = regexp (fstr, '^@\(([^)]*)\)(.*)$', 'tokens');
  if (isempty (tok))
    s = sprintf ("%s ~ f(b, X)", resp_name);
    return;
  endif
  args = strtrim (strsplit (tok{1}{1}, ','));
  ## Strip whitespace first: func2str may render "b (1)" with a space.
  body = strrep (tok{1}{2}, " ", "");
  bname = args{1};
  ## Replace b(k) with the k-th coefficient name (longest index first).
  for k = numel (coef_names):-1:1
    body = strrep (body, sprintf ("%s(%d)", bname, k), coef_names{k});
  endfor
  ## Drop element-wise operator dots for a cleaner display.
  body = strrep (body, ".*", "*");
  body = strrep (body, "./", "/");
  body = strrep (body, ".^", "^");
  s = sprintf ("%s ~ %s", resp_name, body);

endfunction

## ---------------------------------------------------------------------------
function out = ternary (cond, a, b)
  if (cond)
    out = a;
  else
    out = b;
  endif
endfunction

%!demo
%! ## Fit an exponential growth model y = b1 * exp (b2 * x) and inspect it.
%! x = (1:10)';
%! y = [2.1; 2.9; 4.2; 5.3; 7.1; 9.4; 12.8; 16.5; 22.1; 29.8];
%! modelfun = @(b, x) b(1) .* exp (b(2) .* x);
%! mdl = fitnlm (x, y, modelfun, [1; 0.3]);
%! disp (mdl.Coefficients)
%! printf ("RMSE = %g,  R^2 = %g\n", mdl.RMSE, mdl.Rsquared.Ordinary);

## Comprehensive property and method coverage
%!shared X, y, modelfun, beta0
%! X = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%! y = [2.1; 2.9; 4.2; 5.3; 7.1; 9.4; 12.8; 16.5; 22.1; 29.8];
%! modelfun = @(b, x) b(1) .* exp (b(2) .* x);
%! beta0 = [1; 0.3];

%!test  # coefficient table (estimate, SE, tStat) verified against MATLAB
%! mdl = fitnlm (X, y, modelfun, beta0);
%! assert_equal (mdl.Coefficients.Estimate, [1.683747025; 0.286911087], 1e-6);
%! assert_equal (mdl.Coefficients.SE, [0.035194899; 0.002350913], 1e-6);
%! assert_equal (mdl.Coefficients.tStat, [47.8406555; 122.042406], -1e-4);
%! assert_equal (mdl.Coefficients.tStat, ...
%!               mdl.Coefficients.Estimate ./ mdl.Coefficients.SE, 1e-8);

%!test  # sums of squares and their internal relationships
%! mdl = fitnlm (X, y, modelfun, beta0);
%! bhat = mdl.Coefficients.Estimate;  fit = modelfun (bhat, X);
%! assert_equal (mdl.Fitted, fit, 1e-8);
%! assert_equal (mdl.SSE, sum ((y - fit) .^ 2), 1e-8);
%! assert_equal (mdl.SST, sum ((y - mean (y)) .^ 2), 1e-6);
%! assert_equal (mdl.SSR, sum ((fit - mean (y)) .^ 2), 1e-6);
%! assert_equal (mdl.SSE, 0.233771954, 1e-7);
%! assert_equal (mdl.SST, 750.976, 1e-3);

%!test  # MSE/RMSE/DFE and the coefficient of determination
%! mdl = fitnlm (X, y, modelfun, beta0);
%! assert_equal (mdl.DFE, 8);
%! assert_equal (mdl.MSE, mdl.SSE / mdl.DFE, 1e-12);
%! assert_equal (mdl.RMSE, sqrt (mdl.MSE), 1e-12);
%! assert_equal (mdl.RMSE, 0.170942956, 1e-7);
%! assert_equal (mdl.Rsquared.Ordinary, 1 - mdl.SSE / mdl.SST, 1e-12);
%! assert_equal (mdl.Rsquared.Ordinary, 0.999688709, 1e-8);
%! assert_equal (mdl.Rsquared.Adjusted, 0.999649798, 1e-8);

%!test  # log-likelihood and information criteria (values and identities)
%! mdl = fitnlm (X, y, modelfun, beta0);
%! ll = mdl.LogLikelihood;  k = mdl.NumEstimatedCoefficients;  n = 10;
%! assert_equal (ll, 4.590586096, 1e-6);
%! assert_equal (mdl.ModelCriterion.AIC, -5.181172193, 1e-6);
%! assert_equal (mdl.ModelCriterion.BIC, -4.576002007, 1e-6);
%! assert_equal (mdl.ModelCriterion.AIC, -2 * ll + 2 * k, 1e-9);
%! assert_equal (mdl.ModelCriterion.BIC, -2 * ll + k * log (n), 1e-9);
%! assert_equal (mdl.ModelCriterion.AICc, ...
%!               -2 * ll + 2 * k + 2 * k * (k + 1) / (n - k - 1), 1e-9);

%!test  # count/size properties and default names
%! mdl = fitnlm (X, y, modelfun, beta0);
%! assert_equal (mdl.NumCoefficients, 2);
%! assert_equal (mdl.NumEstimatedCoefficients, 2);
%! assert_equal (mdl.NumPredictors, 1);
%! assert_equal (mdl.NumObservations, 10);
%! assert_equal (mdl.CoefficientNames, {'b1', 'b2'});
%! assert_equal (mdl.ResponseName, "y");

%!test  # the coefficient covariance is symmetric with SE^2 on the diagonal
%! mdl = fitnlm (X, y, modelfun, beta0);
%! C = mdl.CoefficientCovariance;
%! assert_equal (size (C), [2, 2]);
%! assert_equal (C, C', 1e-14);
%! assert_equal (diag (C), mdl.Coefficients.SE .^ 2, 1e-12);

%!test  # raw residuals are response minus fit
%! mdl = fitnlm (X, y, modelfun, beta0);
%! assert_equal (class (mdl.Residuals), "table");
%! assert_equal (mdl.Residuals.Raw, y - mdl.Fitted, 1e-10);

%!test  # predict returns fitted values (verified against MATLAB) with CIs
%! mdl = fitnlm (X, y, modelfun, beta0);
%! [yhat, yci] = predict (mdl, [2.5; 5.5; 8.5]);
%! assert_equal (yhat, [3.449741842; 8.158274281; 19.293455074], 1e-6);
%! assert_equal (yci(:,1), [3.329126146; 7.997938483; 19.121921613], 1e-5);
%! assert_equal (yci(:,2), [3.570357538; 8.318610079; 19.464988535], 1e-5);
%! assert_equal (all (yci(:,1) <= yhat & yhat <= yci(:,2)), true);

%!test  # predict at the training data reproduces the fitted response
%! mdl = fitnlm (X, y, modelfun, beta0);
%! assert_equal (predict (mdl, X), mdl.Fitted, 1e-8);

%!test  # feval agrees with predict; random draws match the response size
%! mdl = fitnlm (X, y, modelfun, beta0);
%! assert_equal (feval (mdl, [2.5; 5.5]), predict (mdl, [2.5; 5.5]), 1e-12);
%! ysim = random (mdl);
%! assert_equal (size (ysim), [10, 1]);

%!test  # coefCI matches beta +/- t * SE and honours a custom alpha
%! mdl = fitnlm (X, y, modelfun, beta0);
%! b = mdl.Coefficients.Estimate;  se = mdl.Coefficients.SE;
%! t95 = tinv (0.975, mdl.DFE);
%! assert_equal (coefCI (mdl), [b - t95 * se, b + t95 * se], 1e-12);
%! t90 = tinv (0.95, mdl.DFE);
%! assert_equal (coefCI (mdl, 0.10), [b - t90 * se, b + t90 * se], 1e-12);

%!test  # coefTest reports a Wald F statistic versus the zero model
%! mdl = fitnlm (X, y, modelfun, beta0);
%! [p, F, df] = coefTest (mdl);
%! assert_equal (df, 2);
%! assert_equal (F > 1e5, true);
%! assert_equal (p < 1e-10, true);

%!test  # table input gives the same fit as matrix input
%! tbl = table (X, y, "VariableNames", {'x', 'y'});
%! mdl = fitnlm (tbl, modelfun, beta0);
%! assert_equal (mdl.Coefficients.Estimate, [1.683747025; 0.286911087], 1e-6);
%! assert_equal (mdl.CoefficientNames, {'b1', 'b2'});

%!test  # custom coefficient names are stored and used
%! mdl = fitnlm (X, y, modelfun, beta0, "CoefficientNames", {'A', 'k'});
%! assert_equal (mdl.CoefficientNames, {'A', 'k'});

%!test  # disp prints the model header and the coefficient table
%! mdl = fitnlm (X, y, modelfun, beta0);
%! s = evalc ("disp (mdl)");
%! assert_equal (isempty (strfind (s, "Nonlinear regression model")), false);
%! assert_equal (isempty (strfind (s, "Estimate")), false);

%!test  # chained subsref reaches property -> table column -> element
%! mdl = fitnlm (X, y, modelfun, beta0);
%! assert_equal (numel (mdl.Coefficients.Estimate), 2);
%! assert_equal (mdl.Coefficients.Estimate(1), 1.683747025, 1e-6);

%!test  # the residual and slice plots run without error
%! mdl = fitnlm (X, y, modelfun, beta0);
%! hf = figure ("visible", "off");
%! unwind_protect
%!   plotResiduals (mdl);
%!   plotResiduals (mdl, "fitted");
%!   plotDiagnostics (mdl);
%!   plotSlice (mdl);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error<DATA, RESP, MODELFUN, and BETA0 are required> NonLinearModel (1)
%!error<MODELFUN must be a function handle.> ...
%! NonLinearModel ([1; 2], [1; 2], "bad", [1])
%!error<NonLinearModel: \(\) indexing is not supported> ...
%! mdl = fitnlm ([1;2;3;4], [1;2;3;4], @(b, x) b(1) * x, 1); mdl(1);
