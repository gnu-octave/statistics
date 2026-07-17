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
## @deftypefn {statistics} {} LinearMixedModel
##
## Linear mixed-effects model fitted to data.
##
## A @code{LinearMixedModel} object represents a fitted linear mixed-effects
## model
## @tex
## $y = X\beta + Zb + \varepsilon$,
## @end tex
## @ifnottex
## @code{y = X*beta + Z*b + e},
## @end ifnottex
## with fixed effects @var{beta}, random effects @var{b} distributed as
## @code{N(0, Psi)}, and independent errors @code{N(0, sigma2)}.  Objects are
## created with @code{fitlmematrix} (from design matrices).
##
## The estimated fixed effects and their statistics are available through the
## @code{Coefficients} table; the covariance parameters through
## @code{covarianceParameters}; the random-effect BLUPs through
## @code{randomEffects}; and predictions, residuals, and hypothesis tests
## through the @code{predict}, @code{residuals}, @code{anova}, @code{coefTest},
## and @code{coefCI} methods.
##
## @seealso{fitlmematrix, fitlm}
## @end deftypefn

classdef LinearMixedModel

  properties (GetAccess = public, SetAccess = protected)

    ## Estimation method, "ML" or "REML".
    FitMethod = "";

    ## Number of observations used in the fit.
    NumObservations = [];

    ## Number of fixed-effects coefficients.
    NumCoefficients = [];

    ## Number of estimated fixed-effects coefficients (= NumCoefficients).
    NumEstimatedCoefficients = [];

    ## Table of the fixed-effects estimates with columns Estimate, SE, tStat,
    ## DF, pValue, Lower, and Upper, one row per coefficient.
    Coefficients = [];

    ## Covariance matrix of the fixed-effects estimates.
    CoefficientCovariance = [];

    ## Names of the fixed-effects coefficients.
    CoefficientNames = {};

    ## Maximised log-likelihood (restricted log-likelihood when FitMethod is
    ## "REML").
    LogLikelihood = [];

    ## Struct of information criteria: AIC, BIC, LogLikelihood, and Deviance.
    ModelCriterion = [];

    ## Struct with the Ordinary and Adjusted R-squared statistics.
    Rsquared = [];

    ## Error (residual) sum of squares, sum ((y - conditional fit) .^ 2).
    SSE = [];

    ## Regression sum of squares, sum ((conditional fit - mean (y)) .^ 2).
    SSR = [];

    ## Total sum of squares, SSE + SSR.
    SST = [];

    ## Residual variance estimate (sigma2).
    MSE = [];

    ## Residual (error) degrees of freedom, NumObservations - NumCoefficients.
    DFE = [];

    ## Model formula as a character vector (empty for matrix-input fits).
    Formula = "";

    ## Response variable name (empty for matrix-input fits).
    ResponseName = "";

  endproperties

  properties (Access = private, Hidden)
    X_ = [];  y_ = [];  Zcell_ = {};  Gcell_ = {};
    Zx_ = [];  qk_ = [];  nlev_ = [];  levels_ = {};  gidx_ = {};
    beta_ = [];  covbeta_ = [];  Psi_ = {};  sigma2_ = [];  b_ = [];
    theta_ = [];  GroupNames_ = {};  REPred_ = {};
    fitted_ = [];  fitted_marg_ = [];  resid_ = [];
  endproperties

  methods (Hidden)

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ("%s =\n", in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ("\n  Linear mixed-effects model fit by %s\n", this.FitMethod);
      if (! isempty (this.Formula))
        fprintf ("\n  Formula:\n      %s\n", this.Formula);
      endif
      if (! isempty (this.Coefficients))
        fprintf ("\n  Fixed effects coefficients:\n\n");
        disp (this.Coefficients);
      endif
      if (! isempty (this.Psi_))
        fprintf ("\n  Random effects covariance parameters:\n");
        for k = 1:numel (this.Psi_)
          fprintf ("    Group: %s (%d levels)\n", this.GroupNames_{k}, ...
                   this.nlev_(k));
          disp (this.Psi_{k});
        endfor
        fprintf ("    Residual variance (sigma2): %g\n", this.sigma2_);
      endif
      fprintf ("\n");
      if (! isempty (this.NumObservations))
        fprintf ("Number of observations: %d, Error DF: %d\n", ...
                 this.NumObservations, this.DFE);
      endif
      if (! isempty (this.LogLikelihood))
        fprintf ("Log-likelihood: %g\n", this.LogLikelihood);
      endif
      if (! isempty (this.Rsquared) && isstruct (this.Rsquared))
        fprintf ("R-squared: %g,  Adjusted R-Squared: %g\n", ...
                 this.Rsquared.Ordinary, this.Rsquared.Adjusted);
      endif
    endfunction

    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case "()"
          error (strcat ("LinearMixedModel: () indexing is not supported.", ...
                         "  Use dot notation to access properties."));
        case "{}"
          error (strcat ("LinearMixedModel: {} indexing is not supported.", ...
                         "  Use dot notation to access properties."));
        case "."
          if (! ischar (s.subs))
            error ("LinearMixedModel.subsref: property name must be a string.");
          endif
          ## Method calls (with or without arguments) go straight to builtin.
          if (ismethod (this, s.subs))
            [varargout{1:nargout}] = builtin ("subsref", this, [s, chain_s]);
            return;
          endif
          try
            out = this.(s.subs);
          catch
            error ("LinearMixedModel.subsref: unknown property '%s'.", s.subs);
          end_try_catch
      endswitch
      ## Delegate any further indexing to the property's own subsref
      ## (e.g. lme.Coefficients.Estimate on the returned table).
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {LinearMixedModel} {@var{lme} =} LinearMixedModel (@var{info})
    ##
    ## Construct a @code{LinearMixedModel} from a fitted-model info struct.
    ## This constructor is used internally by @code{fitlmematrix}; call that
    ## function rather than the constructor directly.
    ##
    ## @end deftypefn
    function this = LinearMixedModel (info)

      if (nargin == 0)
        return;
      endif

      n = info.n;
      p = numel (info.beta);

      this.X_ = info.X;         this.y_ = info.y;
      this.Zcell_ = info.Zcell; this.Gcell_ = info.Gcell;
      this.Zx_ = info.Zx;       this.qk_ = info.qk;
      this.nlev_ = info.nlev;   this.levels_ = info.levels;
      this.gidx_ = info.gidx;
      this.beta_ = info.beta;   this.covbeta_ = info.covbeta;
      this.Psi_ = info.Psi;     this.sigma2_ = info.sigma2;
      this.b_ = info.b;         this.theta_ = info.theta;
      this.GroupNames_ = info.GroupNames;
      this.REPred_ = info.REPred;
      this.fitted_ = info.fitted;
      this.fitted_marg_ = info.fitted_marg;
      this.resid_ = info.resid;

      if (isfield (info, "Formula")), this.Formula = info.Formula; endif
      if (isfield (info, "ResponseName"))
        this.ResponseName = info.ResponseName;
      endif

      this.FitMethod = info.method;
      this.NumObservations = n;
      this.NumCoefficients = p;
      this.NumEstimatedCoefficients = p;
      this.CoefficientNames = info.CoefficientNames;
      this.CoefficientCovariance = info.covbeta;
      this.LogLikelihood = info.loglik;
      this.MSE = info.sigma2;
      this.DFE = n - p;

      ## Fixed-effects coefficient statistics (Residual df = n - p).
      se = sqrt (diag (info.covbeta));
      tstat = info.beta ./ se;
      dfe = n - p;
      pval = 2 * (1 - tcdf (abs (tstat), dfe));
      tcrit = tinv (0.975, dfe);
      lower = info.beta - tcrit * se;
      upper = info.beta + tcrit * se;
      this.Coefficients = table (info.beta(:), se(:), tstat(:), ...
        repmat (dfe, p, 1), pval(:), lower(:), upper(:), ...
        "VariableNames", {"Estimate", "SE", "tStat", "DF", "pValue", ...
                          "Lower", "Upper"}, ...
        "RowNames", info.CoefficientNames(:));

      ## Sums of squares and R-squared (MATLAB's mixed-model convention:
      ## SST = SSE + SSR, both from the conditional fit).
      ybar = mean (info.y);
      this.SSE = sum ((info.y - info.fitted) .^ 2);
      this.SSR = sum ((info.fitted - ybar) .^ 2);
      this.SST = this.SSE + this.SSR;
      rsq.Ordinary = 1 - this.SSE / this.SST;
      rsq.Adjusted = 1 - (this.SSE / dfe) / (this.SST / (n - 1));
      this.Rsquared = rsq;

      ## Information criteria (total parameter count = fixed + covariance).
      ncov = sum (arrayfun (@(q) q*(q+1)/2, info.qk)) + 1;
      kpar = p + ncov;
      dev = -2 * info.loglik;
      mc.AIC = dev + 2 * kpar;
      mc.BIC = dev + kpar * log (n);
      mc.LogLikelihood = info.loglik;
      mc.Deviance = dev;
      this.ModelCriterion = mc;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearMixedModel} {@var{beta} =} fixedEffects (@var{lme})
    ## @deftypefnx {LinearMixedModel} {[@var{beta}, @var{names}] =} fixedEffects (@var{lme})
    ##
    ## Return the estimated fixed-effects coefficients @var{beta} and,
    ## optionally, their names.
    ##
    ## @end deftypefn
    function [beta, names] = fixedEffects (this)
      beta = this.beta_;
      names = this.CoefficientNames;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearMixedModel} {@var{b} =} randomEffects (@var{lme})
    ## @deftypefnx {LinearMixedModel} {[@var{b}, @var{names}] =} randomEffects (@var{lme})
    ##
    ## Return the best linear unbiased predictors (BLUPs) of the random effects
    ## @var{b} and, optionally, a cell array of @code{group:level:predictor}
    ## labels.
    ##
    ## @end deftypefn
    function [b, names] = randomEffects (this)
      b = this.b_;
      if (nargout > 1)
        names = {};
        for k = 1:numel (this.qk_)
          pn = this.REPred_{k};
          for l = 1:this.nlev_(k)
            for j = 1:this.qk_(k)
              names{end+1, 1} = sprintf ("%s:%s:%s", this.GroupNames_{k}, ...
                num2str (this.levels_{k}(l)), pn{j});
            endfor
          endfor
        endfor
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {LinearMixedModel} {[@var{psi}, @var{mse}] =} covarianceParameters (@var{lme})
    ##
    ## Return the estimated random-effects covariance matrices @var{psi} (a cell
    ## array, one per grouping term) and the residual variance @var{mse}.
    ##
    ## @end deftypefn
    function [psi, mse] = covarianceParameters (this)
      psi = this.Psi_;
      mse = this.sigma2_;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearMixedModel} {@var{yf} =} fitted (@var{lme})
    ## @deftypefnx {LinearMixedModel} {@var{yf} =} fitted (@var{lme}, @qcode{"Conditional"}, @var{tf})
    ##
    ## Return the fitted values.  With @qcode{"Conditional"} true (the default)
    ## the fit includes the random effects (@code{X*beta + Z*b}); with false it
    ## is the marginal fit (@code{X*beta}).
    ##
    ## @end deftypefn
    function yf = fitted (this, varargin)
      cond = true;
      if (numel (varargin) >= 2 && strcmpi (varargin{1}, "Conditional"))
        cond = logical (varargin{2});
      endif
      if (cond)
        yf = this.fitted_;
      else
        yf = this.fitted_marg_;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearMixedModel} {@var{r} =} residuals (@var{lme})
    ## @deftypefnx {LinearMixedModel} {@var{r} =} residuals (@var{lme}, @qcode{"ResidualType"}, @var{type})
    ##
    ## Return the conditional residuals @code{y - (X*beta + Z*b)}.  @var{type}
    ## is @qcode{"Raw"} (default), @qcode{"Pearson"} (raw divided by
    ## @code{sqrt (sigma2)}), or @qcode{"Standardized"} (raw divided by the
    ## square root of its estimated variance).
    ##
    ## @end deftypefn
    function r = residuals (this, varargin)
      type = "raw";
      if (numel (varargin) >= 2 && strcmpi (varargin{1}, "ResidualType"))
        type = lower (varargin{2});
      endif
      raw = this.resid_;
      switch (type)
        case "raw"
          r = raw;
        case "pearson"
          r = raw / sqrt (this.sigma2_);
        case "standardized"
          r = raw ./ sqrt (cond_resid_var (this));
        otherwise
          error ("LinearMixedModel: unknown ResidualType '%s'.", type);
      endswitch
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearMixedModel} {@var{ypred} =} predict (@var{lme}, @var{Xnew}, @var{Znew}, @var{Gnew})
    ## @deftypefnx {LinearMixedModel} {[@var{ypred}, @var{yci}] =} predict (@dots{})
    ## @deftypefnx {LinearMixedModel} {[@dots{}] =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict the response at new fixed-effects design @var{Xnew},
    ## random-effects design @var{Znew}, and grouping @var{Gnew}.  By default
    ## the prediction is conditional on the estimated random effects (levels of
    ## @var{Gnew} not seen in the fit fall back to the marginal prediction).
    ## With
    ## @qcode{"Conditional"} false the marginal prediction @code{Xnew*beta} is
    ## returned.  The second output @var{yci} gives 95% (or @qcode{"Alpha"})
    ## confidence intervals for the marginal mean.
    ##
    ## @end deftypefn
    function [ypred, yci] = predict (this, Xnew, Znew, Gnew, varargin)

      if (nargin < 2)
        error ("LinearMixedModel: predict requires Xnew.");
      endif
      if (nargin < 3), Znew = []; endif
      if (nargin < 4), Gnew = []; endif
      cond = true;
      alpha = 0.05;
      for i = 1:2:numel (varargin)
        switch (lower (varargin{i}))
          case "conditional"
            cond = logical (varargin{i+1});
          case "alpha"
            alpha = varargin{i+1};
        endswitch
      endfor

      ypred = Xnew * this.beta_;
      if (cond && ! isempty (Znew) && ! isempty (Gnew))
        Gnew = Gnew(:);
        for k = 1:numel (this.qk_)
          q = this.qk_(k);
          lev = this.levels_{k};
          ## offset of term k inside the stacked BLUP vector
          off = 0;
          for kk = 1:k-1
            off += this.qk_(kk) * this.nlev_(kk);
          endfor
          for i = 1:rows (Xnew)
            li = find (lev == Gnew(i), 1);
            if (! isempty (li))
              bk = this.b_(off + (li-1)*q + (1:q));
              ypred(i) += Znew(i, :) * bk;
            endif
          endfor
        endfor
      endif

      if (nargout > 1)
        v = sum ((Xnew * this.covbeta_) .* Xnew, 2);   ## marginal mean variance
        tcrit = tinv (1 - alpha/2, this.DFE);
        half = tcrit * sqrt (v);
        ymarg = Xnew * this.beta_;
        yci = [ymarg - half, ymarg + half];
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearMixedModel} {@var{tbl} =} anova (@var{lme})
    ## @deftypefnx {LinearMixedModel} {@var{tbl} =} anova (@var{lme}, @qcode{"DFMethod"}, @var{method})
    ##
    ## Analysis-of-variance table of F-tests for the fixed-effects terms.  Each
    ## row tests one coefficient.  @var{method} selects the denominator degrees
    ## of freedom: @qcode{"Residual"} (default, @code{n - p}) or
    ## @qcode{"Satterthwaite"}.
    ##
    ## @end deftypefn
    function tbl = anova (this, varargin)
      dfmethod = "residual";
      if (numel (varargin) >= 2 && strcmpi (varargin{1}, "DFMethod"))
        dfmethod = lower (varargin{2});
      endif
      p = this.NumCoefficients;
      beta = this.beta_;
      se = sqrt (diag (this.covbeta_));
      Fstat = (beta ./ se) .^ 2;
      DF1 = ones (p, 1);
      if (strcmp (dfmethod, "satterthwaite"))
        DF2 = __lme_dfsatt__ (this.X_, this.y_, this.Zx_, this.qk_, ...
                this.nlev_, this.Psi_, this.sigma2_, this.FitMethod, eye (p));
      elseif (strcmp (dfmethod, "residual"))
        DF2 = repmat (this.DFE, p, 1);
      else
        error ("LinearMixedModel: unknown DFMethod '%s'.", dfmethod);
      endif
      pValue = 1 - fcdf (Fstat, DF1, DF2);
      tbl = table (Fstat(:), DF1(:), DF2(:), pValue(:), ...
        "VariableNames", {"FStat", "DF1", "DF2", "pValue"}, ...
        "RowNames", this.CoefficientNames(:));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearMixedModel} {@var{p} =} coefTest (@var{lme})
    ## @deftypefnx {LinearMixedModel} {@var{p} =} coefTest (@var{lme}, @var{H})
    ## @deftypefnx {LinearMixedModel} {[@var{p}, @var{F}, @var{df1}, @var{df2}] =} coefTest (@dots{})
    ##
    ## Test the linear hypothesis @code{H*beta = 0} with an F-test.  With no
    ## @var{H}, tests that all non-intercept coefficients are zero (intercept
    ## is taken to be the first coefficient).  Returns the p-value and,
    ## optionally, the F-statistic and its numerator and denominator degrees of
    ## freedom (denominator = @code{n - p}).
    ##
    ## @end deftypefn
    function [pval, F, df1, df2] = coefTest (this, H)
      p = this.NumCoefficients;
      if (nargin < 2)
        H = [zeros(p-1, 1), eye(p-1)];    ## all but the (first) intercept
      endif
      if (columns (H) != p)
        error ("LinearMixedModel: H must have one column per coefficient.");
      endif
      Hb = H * this.beta_;
      df1 = rank (H);
      df2 = this.DFE;
      F = (Hb' * ((H * this.covbeta_ * H') \ Hb)) / df1;
      pval = 1 - fcdf (F, df1, df2);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {LinearMixedModel} {@var{ci} =} coefCI (@var{lme})
    ## @deftypefnx {LinearMixedModel} {@var{ci} =} coefCI (@var{lme}, @var{alpha})
    ##
    ## Confidence intervals for the fixed-effects coefficients at level
    ## @code{1 - @var{alpha}} (default @var{alpha} = 0.05).  Row @math{j} holds
    ## the lower and upper bounds for coefficient @math{j}.
    ##
    ## @end deftypefn
    function ci = coefCI (this, alpha)
      if (nargin < 2), alpha = 0.05; endif
      se = sqrt (diag (this.covbeta_));
      tcrit = tinv (1 - alpha/2, this.DFE);
      ci = [this.beta_ - tcrit * se, this.beta_ + tcrit * se];
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {LinearMixedModel} {@var{D} =} designMatrix (@var{lme}, @var{type})
    ##
    ## Return the fixed-effects design matrix (@var{type} = @qcode{"Fixed"},
    ## default) or the expanded random-effects design matrix
    ## (@var{type} = @qcode{"Random"}).
    ##
    ## @end deftypefn
    function D = designMatrix (this, type)
      if (nargin < 2), type = "Fixed"; endif
      switch (lower (type))
        case "fixed"
          D = this.X_;
        case "random"
          D = this.Zx_;
        otherwise
          error ("LinearMixedModel: type must be 'Fixed' or 'Random'.");
      endswitch
    endfunction

  endmethods

  methods (Access = private)

    ## Diagonal of the covariance of the conditional residuals, for standardized
    ## residuals: r_c = Mc*y with Mc = (I - Hx) - Zx*Dabs*Zx'*P.
    function v = cond_resid_var (this)
      n = this.NumObservations;
      ## rebuild the marginal covariance V and the absolute RE covariance Dabs
      blocks = {};
      for k = 1:numel (this.qk_)
        for l = 1:this.nlev_(k)
          blocks{end+1} = this.Psi_{k};
        endfor
      endfor
      Dabs = blkdiag (blocks{:});
      V = this.Zx_ * Dabs * this.Zx_' + this.sigma2_ * eye (n);
      Vi = inv (V);
      X = this.X_;
      XtViX = X' * Vi * X;
      Hx = X * (XtViX \ (X' * Vi));
      P = Vi - Vi * X * (XtViX \ (X' * Vi));
      Mc = (eye (n) - Hx) - this.Zx_ * Dabs * this.Zx_' * P;
      v = diag (Mc * V * Mc');
    endfunction

  endmethods

endclassdef

## Shared MATLAB-verified fixture (fitlmematrix R2026a): random-intercept REML
## model y ~ x + x2 + (1|g) on 42 observations, 6 groups.
%!shared X, yL, grp, xL, x2, lme
%! xL = [0.032760004 0.70410822 -0.8646718 -0.28869454 0.51276678 -1.4975462 ...
%!  -1.4527871 -0.80013541 -1.644209 1.5137701 0.72905543 0.20880758 1.0856145 ...
%!  0.62862577 -0.87409978 1.9178276 0.09748204 0.50697633 1.0247569 ...
%!  -0.92789896 -0.88921018 -0.98322849 -0.031378913 0.86875961 -0.91481141 ...
%!  0.034324163 -0.25025257 -1.0575644 -0.86131607 -0.35355444 0.82950729 ...
%!  -0.36874363 0.061580868 0.55803564 -0.1763803 1.0482413 1.0137831 ...
%!  -0.94876976 -0.010703972 -0.35149845 -1.6828735 -1.0493301]';
%! x2 = [0.68979276 0.0074354814 -0.45697437 -0.5636481 1.4567202 -0.97829955 ...
%!  -1.12922 -0.030542479 1.5847779 -0.87837755 0.24121762 0.68747601 ...
%!  -0.56728765 0.98895053 -0.39350661 0.85326015 0.36524343 0.15824977 ...
%!  -1.7665212 0.59808246 -0.55763708 -1.1982294 -2.1473319 0.22521416 ...
%!  0.37034398 -1.880586 0.052941033 -0.70016994 0.2174853 -1.7797082 ...
%!  0.51971317 -0.35551286 1.9845963 -1.3498848 -0.63514097 -0.78794714 ...
%!  1.3681179 1.4423152 -0.51233905 0.30238864 2.0458136 0.17326323]';
%! yL = [3.5635971 -0.36763498 1.4192715 3.4471073 2.2760668 3.6420287 ...
%!  4.4481522 1.5292777 5.437158 0.44016873 0.7259756 2.9771749 1.2941347 ...
%!  0.26407508 1.9673466 0.77589984 1.2176078 0.8528384 -0.86522876 2.8651388 ...
%!  2.0127931 3.1242841 -0.49164361 1.3192295 5.3782287 -0.40680701 2.3989882 ...
%!  3.606108 1.9193197 1.6713765 2.4383124 1.9314844 3.5112706 0.91644553 ...
%!  0.034923706 -0.15510033 2.1765206 2.5327412 3.1421219 3.3472574 5.343425 ...
%!  3.8775899]';
%! grp = [1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 ...
%!  5 6 1 2 3 4 5 6]';
%! X = [ones(42,1), xL, x2];
%! lme = fitlmematrix (X, yL, ones (42, 1), grp, "FitMethod", "REML", ...
%!         "FixedEffectPredictors", {"(Intercept)", "x", "x2"});

%!test  # returns a LinearMixedModel and basic properties
%! assert (isa (lme, "LinearMixedModel"));
%! assert (lme.NumObservations, 42);
%! assert (lme.NumCoefficients, 3);
%! assert (lme.DFE, 39);
%! assert (lme.FitMethod, "REML");

%!test  # Coefficients table -- SE / tStat / DF / pValue / CI vs MATLAB
%! C = lme.Coefficients;
%! assert (C.Estimate, [1.96839; -1.37926; 0.811747], 1e-4);
%! assert (C.SE, [0.376509; 0.113264; 0.0966571], 1e-5);
%! assert (C.tStat, [5.22801; -12.1773; 8.39822], 1e-4);
%! assert (C.DF, [39; 39; 39]);
%! assert (C.pValue, [6.08512e-06; 7.30179e-15; 2.8079e-10], 1e-8);
%! assert (C.Lower, [1.20683; -1.60835; 0.61624], 1e-4);
%! assert (C.Upper, [2.72996; -1.15016; 1.00725], 1e-4);

%!test  # anova with Residual DF (default): F = t^2, DF2 = n - p
%! a = anova (lme);
%! assert (a.FStat, [27.3321; 148.287; 70.5301], 1e-3);
%! assert (a.DF1, [1; 1; 1]);
%! assert (a.DF2, [39; 39; 39]);

%!test  # anova with Satterthwaite DF vs MATLAB
%! a = anova (lme, "DFMethod", "Satterthwaite");
%! assert (a.DF2, [4.95056; 34.5037; 34.3302], 1e-3);
%! assert (a.pValue, [0.00348739; 4.79374e-14; 7.71731e-10], 1e-6);

%!test  # coefTest: joint test of the non-intercept coefficients
%! [p, F, df1, df2] = coefTest (lme);
%! assert (F, 106.847, 1e-2);
%! assert (df1, 2);  assert (df2, 39);
%! assert (p < 1e-12);   # ~1e-16, below the precision of 1 - fcdf

%!test  # coefTest with an explicit single-row hypothesis (the x coefficient)
%! [p, F, df1] = coefTest (lme, [0 1 0]);
%! assert (F, 148.287, 1e-2);
%! assert (df1, 1);

%!test  # coefCI matches the Coefficients table bounds
%! ci = coefCI (lme);
%! assert (ci(:,1), lme.Coefficients.Lower, 1e-12);
%! assert (ci(:,2), lme.Coefficients.Upper, 1e-12);

%!test  # residuals: raw, Pearson (= raw/sqrt(mse)), Standardized vs MATLAB
%! [~, mse] = covarianceParameters (lme);
%! raw = residuals (lme);
%! assert (raw, yL - fitted (lme), 1e-12);
%! assert (residuals (lme, "ResidualType", "Pearson"), raw / sqrt (mse), 1e-12);
%! rs = residuals (lme, "ResidualType", "Standardized");
%! assert (rs(1:3), [0.1824776; -0.4452186; -2.057201], 1e-5);

%!test  # fitted: conditional includes random effects, marginal does not
%! fc = fitted (lme);
%! fm = fitted (lme, "Conditional", false);
%! assert (fm, X * fixedEffects (lme), 1e-12);
%! assert (any (abs (fc - fm) > 1e-3));

%!test  # predict: conditional (unseen group falls back to marginal) + marginal
%! Xn = [1 0.5 0; 1 -0.5 1; 1 1 -1];
%! yc = predict (lme, Xn, ones (3, 1), [1; 2; 7]);
%! ym = predict (lme, Xn, ones (3, 1), [1; 2; 7], "Conditional", false);
%! assert (yc, [2.254803; 2.351467; -0.2226094], 1e-4);
%! assert (ym, [1.278766; 3.469768; -0.2226094], 1e-4);
%! assert (yc(3), ym(3), 1e-10);   # group 7 unseen -> conditional == marginal

%!test  # predict confidence intervals bracket the marginal mean
%! Xn = [1 0.5 0; 1 -0.5 1];
%! [yp, ci] = predict (lme, Xn, ones (2, 1), [1; 2], "Conditional", false);
%! assert (all (ci(:,1) < yp & yp < ci(:,2)));

%!test  # effect extraction and design matrices
%! [beta, names] = fixedEffects (lme);
%! assert (beta, lme.Coefficients.Estimate, 1e-12);
%! assert (names, {"(Intercept)", "x", "x2"});
%! b = randomEffects (lme);
%! assert (numel (b), 6);                 # 6 intercept BLUPs
%! assert (size (designMatrix (lme, "Fixed")), [42, 3]);
%! assert (size (designMatrix (lme, "Random")), [42, 6]);

%!test  # R-squared and sums of squares (MATLAB SST = SSE + SSR convention)
%! assert (lme.Rsquared.Ordinary, 0.8758549, 1e-6);
%! assert (lme.Rsquared.Adjusted, 0.8694885, 1e-6);
%! assert (lme.SST, lme.SSE + lme.SSR, 1e-10);
%! assert (lme.ModelCriterion.Deviance, -2 * lme.LogLikelihood, 1e-10);

## Error handling
%!error <unknown ResidualType> residuals (lme, "ResidualType", "xxx")
%!error <unknown DFMethod> anova (lme, "DFMethod", "xxx")
%!error <type must be 'Fixed' or 'Random'> designMatrix (lme, "bogus")
%!error <H must have one column per coefficient> coefTest (lme, [1 0])
%!error <indexing is not supported> lme(1)
