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
## @deftypefn {statistics} {} GeneralizedLinearMixedModel
##
## Generalized linear mixed-effects model fitted to data.
##
## A @code{GeneralizedLinearMixedModel} object represents a fitted generalized
## linear mixed-effects model: a generalized linear model whose linear predictor
## @code{X*beta + Z*b} includes normally distributed random effects
## @code{b ~ N(0, Psi)}.  Objects are created with @code{fitglme}.
##
## The model is fitted by penalized quasi-likelihood.  The fixed-effects
## estimates and their statistics are available through the @code{Coefficients}
## table, the covariance parameters through @code{covarianceParameters}, and
## predictions, residuals, and hypothesis tests through the @code{predict},
## @code{residuals}, @code{anova}, @code{coefTest}, and @code{coefCI} methods.
##
## @seealso{fitglme, fitlme, GeneralizedLinearModel}
## @end deftypefn

classdef GeneralizedLinearMixedModel

  properties (GetAccess = public, SetAccess = protected)

    ## Response distribution ("binomial", "poisson", or "normal").
    Distribution = "";

    ## Link function ("logit", "log", or "identity").
    Link = "";

    ## Fitting method ("MPL", "REMPL", "Laplace", or "ApproximateLaplace").
    FitMethod = "";

    ## Estimated dispersion parameter (1 for binomial and poisson).
    Dispersion = [];

    ## Number of observations used in the fit.
    NumObservations = [];

    ## Number of fixed-effects coefficients.
    NumCoefficients = [];

    ## Table of the fixed-effects estimates (Estimate, SE, tStat, DF, pValue,
    ## Lower, Upper), one row per coefficient.
    Coefficients = [];

    ## Covariance matrix of the fixed-effects estimates.
    CoefficientCovariance = [];

    ## Names of the fixed-effects coefficients.
    CoefficientNames = {};

    ## Maximised (pseudo- or Laplace-approximated) log-likelihood.
    LogLikelihood = [];

    ## Struct of information criteria: AIC, BIC, LogLikelihood, and Deviance.
    ModelCriterion = [];

    ## Residual (error) degrees of freedom, NumObservations - NumCoefficients.
    DFE = [];

    ## Model formula (empty for matrix-input fits).
    Formula = "";

    ## Response variable name (empty for matrix-input fits).
    ResponseName = "";

  endproperties

  properties (Access = private, Hidden)
    X_ = [];  y_ = [];  beta_ = [];  covbeta_ = [];  Psi_ = {};  b_ = [];
    mu_ = [];  resraw_ = [];  respear_ = [];
    Zx_ = [];  qk_ = [];  nlev_ = [];  levels_ = {};  gidx_ = {};
    GroupNames_ = {};  REPred_ = {};  distr_ = "";  link_ = "";
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
      fprintf ("\n  Generalized linear mixed-effects model fit by %s\n", ...
               this.FitMethod);
      fprintf ("    Distribution: %s, Link: %s\n", ...
               this.Distribution, this.Link);
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
      endif
      fprintf ("\n");
      if (! isempty (this.NumObservations))
        fprintf ("Number of observations: %d, Error DF: %d\n", ...
                 this.NumObservations, this.DFE);
      endif
      if (! isempty (this.LogLikelihood))
        fprintf ("Log-likelihood: %g\n", this.LogLikelihood);
      endif
    endfunction

    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case "()"
          error (strcat ("GeneralizedLinearMixedModel: () indexing is not", ...
                         " supported; use dot notation."));
        case "{}"
          error (strcat ("GeneralizedLinearMixedModel: {} indexing is not", ...
                         " supported; use dot notation."));
        case "."
          if (! ischar (s.subs))
            error ("GeneralizedLinearMixedModel.subsref: invalid index.");
          endif
          if (ismethod (this, s.subs))
            [varargout{1:nargout}] = builtin ("subsref", this, [s, chain_s]);
            return;
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("GeneralizedLinearMixedModel.subsref: ", ...
                           "unknown property '%s'."), s.subs);
          end_try_catch
      endswitch
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearMixedModel} {@var{glme} =} GeneralizedLinearMixedModel (@var{info})
    ##
    ## Construct from a fitted-model info struct.  Used internally by
    ## @code{fitglme}; call that function rather than the constructor directly.
    ##
    ## @end deftypefn
    function this = GeneralizedLinearMixedModel (info)
      if (nargin == 0)
        return;
      endif
      n = info.n;
      p = info.p;
      this.X_ = info.X;  this.y_ = info.y;
      this.beta_ = info.beta;  this.covbeta_ = info.covbeta;
      this.Psi_ = info.Psi;  this.b_ = info.b;  this.mu_ = info.mu;
      this.resraw_ = info.resid_raw;  this.respear_ = info.resid_pearson;
      this.Zx_ = info.Zx;  this.qk_ = info.qk;  this.nlev_ = info.nlev;
      this.levels_ = info.levels;  this.gidx_ = info.gidx;
      this.GroupNames_ = info.GroupNames;  this.REPred_ = info.REPred;
      this.distr_ = info.distr;  this.link_ = info.link;

      this.Distribution = info.distr;
      this.Link = info.link;
      this.FitMethod = info.FitMethod;
      this.Dispersion = info.dispersion;
      this.NumObservations = n;
      this.NumCoefficients = p;
      this.CoefficientNames = info.CoefficientNames;
      this.CoefficientCovariance = info.covbeta;
      this.LogLikelihood = info.loglik;
      this.DFE = n - p;
      if (isfield (info, "Formula")), this.Formula = info.Formula; endif
      if (isfield (info, "ResponseName"))
        this.ResponseName = info.ResponseName;
      endif

      se = sqrt (diag (info.covbeta));
      tstat = info.beta ./ se;
      dfe = n - p;
      pval = 2 * (1 - tcdf (abs (tstat), dfe));
      tcrit = tinv (0.975, dfe);
      this.Coefficients = table (info.beta(:), se(:), tstat(:), ...
        repmat (dfe, p, 1), pval(:), info.beta(:) - tcrit*se(:), ...
        info.beta(:) + tcrit*se(:), ...
        "VariableNames", {"Estimate", "SE", "tStat", "DF", "pValue", ...
                          "Lower", "Upper"}, ...
        "RowNames", info.CoefficientNames(:));

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
    ## @deftypefn {GeneralizedLinearMixedModel} {[@var{beta}, @var{names}] =} fixedEffects (@var{glme})
    ## Return the fixed-effects coefficients and, optionally, their names.
    ## @end deftypefn
    function [beta, names] = fixedEffects (this)
      beta = this.beta_;
      names = this.CoefficientNames;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearMixedModel} {@var{b} =} randomEffects (@var{glme})
    ## Return the estimated random-effects (the conditional modes).
    ## @end deftypefn
    function b = randomEffects (this)
      b = this.b_;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearMixedModel} {[@var{psi}, @var{disp}] =} covarianceParameters (@var{glme})
    ## Return the random-effects covariance matrices and the dispersion.
    ## @end deftypefn
    function [psi, dispn] = covarianceParameters (this)
      psi = this.Psi_;
      dispn = this.Dispersion;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearMixedModel} {@var{yf} =} fitted (@var{glme})
    ## Return the fitted mean response (conditional on the random effects).
    ## @end deftypefn
    function yf = fitted (this)
      yf = this.mu_;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearMixedModel} {@var{r} =} residuals (@var{glme})
    ## @deftypefnx {GeneralizedLinearMixedModel} {@var{r} =} residuals (@var{glme}, @qcode{"ResidualType"}, @var{type})
    ## Return @qcode{"Raw"} (default) or @qcode{"Pearson"} residuals.
    ## @end deftypefn
    function r = residuals (this, varargin)
      type = "raw";
      if (numel (varargin) >= 2 && strcmpi (varargin{1}, "ResidualType"))
        type = lower (varargin{2});
      endif
      switch (type)
        case "raw"
          r = this.resraw_;
        case "pearson"
          r = this.respear_;
        otherwise
          error (strcat ("GeneralizedLinearMixedModel: unknown", ...
                         " ResidualType '%s'."), type);
      endswitch
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearMixedModel} {@var{ypred} =} predict (@var{glme}, @var{Xnew}, @var{Znew}, @var{Gnew})
    ## Predict the mean response at new data.  With @qcode{"Conditional"} true
    ## (default) the random effects of known grouping levels are added;
    ## unknown levels fall back to the marginal (fixed-effects) prediction.
    ## @end deftypefn
    function ypred = predict (this, Xnew, Znew, Gnew, varargin)
      if (nargin < 4), Znew = []; Gnew = []; endif
      cond = true;
      for i = 1:2:numel (varargin)
        if (strcmpi (varargin{i}, "Conditional"))
          cond = logical (varargin{i+1});
        endif
      endfor
      eta = Xnew * this.beta_;
      if (cond && ! isempty (Znew) && ! isempty (Gnew))
        Gnew = Gnew(:);
        q = this.qk_(1);
        lev = this.levels_{1};
        for i = 1:rows (Xnew)
          li = find (lev == Gnew(i), 1);
          if (! isempty (li))
            eta(i) += Znew(i, :) * this.b_((li-1)*q + (1:q));
          endif
        endfor
      endif
      ypred = link_inv (eta, this.link_);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearMixedModel} {@var{tbl} =} anova (@var{glme})
    ## Analysis-of-deviance table of F-tests for the fixed-effects terms, using
    ## residual denominator degrees of freedom.
    ## @end deftypefn
    function tbl = anova (this)
      p = this.NumCoefficients;
      se = sqrt (diag (this.covbeta_));
      Fstat = (this.beta_ ./ se) .^ 2;
      DF1 = ones (p, 1);
      DF2 = repmat (this.DFE, p, 1);
      pValue = 1 - fcdf (Fstat, DF1, DF2);
      tbl = table (Fstat(:), DF1(:), DF2(:), pValue(:), ...
        "VariableNames", {"FStat", "DF1", "DF2", "pValue"}, ...
        "RowNames", this.CoefficientNames(:));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearMixedModel} {[@var{p}, @var{F}, @var{df1}, @var{df2}] =} coefTest (@var{glme}, @var{H})
    ## F-test of the linear hypothesis @code{H*beta = 0}.
    ## @end deftypefn
    function [pval, F, df1, df2] = coefTest (this, H)
      p = this.NumCoefficients;
      if (nargin < 2)
        H = [zeros(p-1, 1), eye(p-1)];
      endif
      Hb = H * this.beta_;
      df1 = rank (H);
      df2 = this.DFE;
      F = (Hb' * ((H * this.covbeta_ * H') \ Hb)) / df1;
      pval = 1 - fcdf (F, df1, df2);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearMixedModel} {@var{ci} =} coefCI (@var{glme}, @var{alpha})
    ## Confidence intervals for the fixed-effects coefficients.
    ## @end deftypefn
    function ci = coefCI (this, alpha)
      if (nargin < 2), alpha = 0.05; endif
      se = sqrt (diag (this.covbeta_));
      tcrit = tinv (1 - alpha/2, this.DFE);
      ci = [this.beta_ - tcrit * se, this.beta_ + tcrit * se];
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {GeneralizedLinearMixedModel} {@var{D} =} designMatrix (@var{glme}, @var{type})
    ## Return the fixed (@qcode{"Fixed"}, default) or random (@qcode{"Random"})
    ## design matrix.
    ## @end deftypefn
    function D = designMatrix (this, type)
      if (nargin < 2), type = "Fixed"; endif
      switch (lower (type))
        case "fixed"
          D = this.X_;
        case "random"
          D = this.Zx_;
        otherwise
          error (strcat ("GeneralizedLinearMixedModel: type must be", ...
                         " 'Fixed' or 'Random'."));
      endswitch
    endfunction

  endmethods

endclassdef

## Inverse link (canonical), shared with the fitting engine's convention.
function mu = link_inv (eta, link)
  switch (lower (link))
    case "logit"
      mu = 1 ./ (1 + exp (-eta));
    case "log"
      mu = exp (eta);
    case "identity"
      mu = eta;
    otherwise
      error ("GeneralizedLinearMixedModel: unsupported link '%s'.", link);
  endswitch
endfunction

## Shared MATLAB-verified fixture (fitglme R2026a): poisson random intercept.
%!shared glme, tbl
%! xL = [0.032760004 0.70410822 -0.8646718 -0.28869454 0.51276678 -1.4975462 ...
%!  -1.4527871 -0.80013541 -1.644209 1.5137701 0.72905543 0.20880758 1.0856145 ...
%!  0.62862577 -0.87409978 1.9178276 0.09748204 0.50697633 1.0247569 ...
%!  -0.92789896 -0.88921018 -0.98322849 -0.031378913 0.86875961 -0.91481141 ...
%!  0.034324163 -0.25025257 -1.0575644 -0.86131607 -0.35355444 0.82950729 ...
%!  -0.36874363 0.061580868 0.55803564 -0.1763803 1.0482413 1.0137831 ...
%!  -0.94876976 -0.010703972 -0.35149845 -1.6828735 -1.0493301]';
%! yPois = [3 3 1 1 1 1 1 2 1 5 2 0 5 0 1 5 0 2 0 0 1 0 0 5 3 0 1 0 0 1 1 1 2 2 ...
%!  1 1 4 0 1 0 0 1]';
%! g = [1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 ...
%!  6 1 2 3 4 5 6]';
%! tbl = table (yPois, xL, g);
%! glme = fitglme (tbl, "yPois ~ xL + (1 | g)", "Distribution", "poisson", ...
%!                 "FitMethod", "REMPL");

%!test  # object type and basic properties
%! assert (isa (glme, "GeneralizedLinearMixedModel"));
%! assert (glme.NumObservations, 42);
%! assert (glme.NumCoefficients, 2);
%! assert (glme.DFE, 40);
%! assert (glme.Dispersion, 1);

%!test  # Coefficients table
%! C = glme.Coefficients;
%! assert (C.Estimate, [0.23092; 0.67809], 1e-3);
%! assert (C.SE, [0.15395; 0.15220], 1e-3);
%! assert (C.DF, [40; 40]);
%! assert (C.tStat, C.Estimate ./ C.SE, 1e-10);

%!test  # effect extraction and covariance parameters
%! [beta, names] = fixedEffects (glme);
%! assert (beta, glme.Coefficients.Estimate, 1e-12);
%! assert (names(:), {"(Intercept)"; "xL"});
%! b = randomEffects (glme);
%! assert (numel (b), 6);
%! [psi, dispn] = covarianceParameters (glme);
%! assert (psi{1}, 0.015918, 1e-3);
%! assert (dispn, 1);

%!test  # anova F-tests: F = tStat^2 with residual DF
%! a = anova (glme);
%! assert (a.FStat, (glme.Coefficients.tStat) .^ 2, 1e-8);
%! assert (a.DF2, [40; 40]);

%!test  # coefTest and coefCI
%! [p, F, df1, df2] = coefTest (glme, [0 1]);
%! assert (F, glme.Coefficients.tStat(2) ^ 2, 1e-6);
%! assert (df1, 1);  assert (df2, 40);
%! ci = coefCI (glme);
%! assert (ci(:,1), glme.Coefficients.Lower, 1e-12);

%!test  # fitted values are positive counts (log link) and residuals sum sensibly
%! mu = fitted (glme);
%! assert (all (mu > 0));
%! assert (residuals (glme), tbl.yPois - mu, 1e-12);
%! assert (residuals (glme, "ResidualType", "Pearson"), (tbl.yPois - mu) ./ sqrt (mu), 1e-10);

%!test  # predict: conditional (known group) and marginal (unseen group)
%! ym = predict (glme, [1 0.5], [], [], "Conditional", false);
%! yc = predict (glme, [1 0.5], 1, 1);
%! yu = predict (glme, [1 0.5], 1, 99);   # unseen group -> marginal
%! assert (ym, exp (glme.Coefficients.Estimate' * [1; 0.5]), 1e-10);
%! assert (yu, ym, 1e-12);

%!test  # designMatrix and ModelCriterion
%! assert (size (designMatrix (glme, "Fixed")), [42, 2]);
%! assert (size (designMatrix (glme, "Random")), [42, 6]);
%! assert (glme.ModelCriterion.Deviance, -2 * glme.LogLikelihood, 1e-10);

## Error handling
%!error <unknown ResidualType> residuals (glme, "ResidualType", "xxx")
%!error <type must be 'Fixed' or 'Random'> designMatrix (glme, "bogus")
%!error <indexing is not supported> glme(1)
