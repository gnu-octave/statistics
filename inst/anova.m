## Copyright (C) 2026 Aman Behera <aman.behera.systesms@gmail.com>
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

classdef anova < handle
  ## -*- texinfo -*-
  ## @deftp {statistics} anova
  ##
  ## Object-oriented interface for analysis of variance.
  ##
  ## The @code{anova} class provides a unified object-oriented wrapper around
  ## the procedural functions @code{anova1}, @code{anova2}, and @code{anovan}.
  ## It stores the response data, grouping variables, model specification, and
  ## fitted results in one object.  The class chooses the narrowest compatible
  ## backend, delegates the numeric computation to the existing ANOVA
  ## functions, and exposes common follow-up operations such as
  ## @code{summary}, @code{multcompare}, @code{plotDiagnostics},
  ## @code{predict}, and @code{getEffectSizes}.
  ##
  ## Models are fitted lazily.  Methods that need fitted results call
  ## @code{fit} internally when necessary, so users may construct an object and
  ## immediately call inspection or post-hoc methods.
  ##
  ## @seealso{anova1, anova2, anovan, multcompare, fitlm, LinearModel}
  ## @end deftp

  properties (GetAccess = public, SetAccess = private)
    ## Data
    Y                                       ## Response (vector or matrix)
    GROUP                                   ## Grouping variables

    ## Results (populated after fit; empty until then)
    Coefficients   = [];
    AnovaTable     = {};
    Residuals      = [];
    FittedValues   = [];
    DFE            = [];
    MSE            = [];
    DesignMatrix   = [];
    Stats          = struct ();
  endproperties

  properties (Access = public)
    ## Model specification (mutations will trigger lazy refit in Week 2+)
    ModelType   = 'linear';                 ## model name or terms matrix
    SSType      = 3;                        ## 1, 2, or 3
    VarNames    = {};
    Contrasts   = {};
    Alpha       = 0.05;
    Continuous  = [];                       ## indices of continuous predictors
    Random      = [];                       ## indices of random factors
    Weights     = [];
    Display     = 'off';                    ## 'on' or 'off'
  endproperties

  properties (Access = private)
    fitted_     = false;
    dirty_      = true;
    nFactors_   = 0;
    backend_    = '';                       ## 'anova1' | 'anova2' | 'anovan'
    reps_       = [];                       ## replicate count for anova2 backend
    sourceModel_ = [];                       ## LinearModel object, when supplied
  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {@var{obj} =} anova (@var{Y})
    ## @deftypefnx {anova} {@var{obj} =} anova (@var{Y}, @var{GROUP})
    ## @deftypefnx {anova} {@var{obj} =} anova (@dots{}, @var{name}, @var{value})
    ## @deftypefnx {anova} {@var{obj} =} anova (@var{mdl})
    ##
    ## Create an object-oriented analysis of variance model.
    ##
    ## @var{Y} is a non-empty numeric response vector or matrix.  @var{GROUP}
    ## contains grouping variables for vector responses and may be a grouping
    ## vector, a matrix of grouping variables, or a cell array of grouping
    ## vectors.  If @var{GROUP} is omitted and @var{Y} is a matrix, columns of
    ## @var{Y} are treated as groups following @code{anova1} matrix syntax.
    ##
    ## The constructor accepts name-value arguments matching the underlying
    ## ANOVA backends, including @qcode{'Model'}, @qcode{'SSType'},
    ## @qcode{'VarNames'}, @qcode{'Contrasts'}, @qcode{'Alpha'},
    ## @qcode{'Continuous'}, @qcode{'Random'}, @qcode{'Weights'},
    ## @qcode{'Display'}, and @qcode{'Reps'}.  Passing @qcode{'Reps'} selects
    ## the balanced two-way @code{anova2} backend when @var{Y} is a non-vector
    ## matrix.
    ##
    ## @var{mdl} may be a @code{LinearModel} object, in which case the ANOVA
    ## object is populated from the fitted linear model's public properties.
    ##
    ## @seealso{fit, summary, anova1, anova2, anovan, LinearModel}
    ## @end deftypefn
    function obj = anova (Y, GROUP, varargin)

      if (nargin < 1)
        error ("anova: too few input arguments.");
      endif
      if (nargin < 2)
        GROUP = [];
      endif

      if (isa (Y, 'LinearModel'))
        lm_args = varargin;
        if (! isempty (GROUP))
          lm_args = [{GROUP}, lm_args];
        endif
        obj.initLinearModel_ (Y, lm_args{:});
        return;
      endif

      if (! isnumeric (Y) || isempty (Y))
        error ("anova: Y must be a non-empty numeric array.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error ("anova: name-value pairs must come in pairs.");
      endif

      obj.Y     = Y;
      obj.GROUP = GROUP;

      ## Parse name-value pairs (mirrors anovan.m's loop style)
      for idx = 1:2:numel (varargin)
        name  = varargin{idx};
        value = varargin{idx + 1};
        if (! ischar (name))
          error ("anova: parameter name must be a character vector.");
        endif
        switch (lower (name))
          case 'model'
            obj.ModelType = value;
          case 'sstype'
            obj.SSType = value;
          case 'varnames'
            obj.VarNames = value;
          case 'contrasts'
            obj.Contrasts = value;
          case 'alpha'
            obj.Alpha = value;
          case 'continuous'
            obj.Continuous = value;
          case 'random'
            obj.Random = value;
          case 'weights'
            obj.Weights = value;
          case {'display', 'displayopt'}
            obj.Display = value;
          case 'reps'
            obj.reps_ = value;
          otherwise
            error ("anova: parameter '%s' is not supported.", name);
        endswitch
      endfor

      obj.nFactors_ = obj.countFactors_ ();
      obj.validateSpec_ ();
      obj.validateData_ ();
      obj.selectBackend_ ();
      obj.dirty_ = true;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {anova} {} fit (@var{obj})
    ##
    ## Fit the ANOVA model if it has not already been fitted.
    ##
    ## This method is optional for most workflows because methods such as
    ## @code{summary}, @code{multcompare}, @code{plotDiagnostics},
    ## @code{predict}, and @code{getEffectSizes} call @code{fit} lazily when
    ## fitted results are needed.
    ##
    ## @end deftypefn
    function fit (obj)
      obj.ensureFit_ ();
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {anova} {} summary (@var{obj})
    ##
    ## Display the fitted ANOVA table and basic fit statistics.
    ##
    ## The model is fitted first if needed.  The printed table uses the
    ## backend table returned by @code{anova1}, @code{anova2}, or
    ## @code{anovan}, followed by the mean squared error, error degrees of
    ## freedom, and significance level.
    ##
    ## @seealso{disp, fit}
    ## @end deftypefn
    function summary (obj)
      obj.ensureFit_ ();
      atab = obj.AnovaTable;
      if (isempty (atab))
        fprintf ("  anova: no results to display.\n");
        return;
      endif
      sstype_char = obj.sstypeLabel_ ();
      fprintf ("\nANOVA TABLE (Type %s sums-of-squares, backend = %s):\n\n", ...
               sstype_char, obj.backend_);
      obj.printAtab_ (atab);
      if (! isempty (obj.MSE))
        fprintf ("\nMSE: %g    DFE: %g    Alpha: %g\n", ...
                 obj.MSE, obj.DFE, obj.Alpha);
      endif
      fprintf ("\n");
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {anova} {} disp (@var{obj})
    ##
    ## Display a compact description of an @code{anova} object.
    ##
    ## The display includes the selected backend, fit state, number of factors,
    ## sum-of-squares type, and significance level.
    ##
    ## @seealso{summary}
    ## @end deftypefn
    function disp (obj)
      fprintf ("\n  anova object\n");
      fprintf ("    backend  : %s\n", obj.backend_);
      fprintf ("    fitted   : %d\n", obj.fitted_);
      fprintf ("    nFactors : %d\n", obj.nFactors_);
      fprintf ("    SSType   : %d\n", obj.SSType);
      fprintf ("    Alpha    : %g\n\n", obj.Alpha);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {@var{C} =} multcompare (@var{obj})
    ## @deftypefnx {anova} {[@var{C}, @var{M}, @var{H}, @var{GNAMES}] =} multcompare (@var{obj}, @dots{})
    ##
    ## Perform post-hoc multiple comparisons for a fitted ANOVA object.
    ##
    ## The method fits the object if needed and delegates to the package
    ## function @code{multcompare} using the backend @code{Stats} structure.
    ## Additional arguments are passed through unchanged.
    ##
    ## @seealso{multcompare}
    ## @end deftypefn
    function varargout = multcompare (obj, varargin)
      obj.ensureFit_ ();
      if (strcmp (obj.backend_, 'linearmodel'))
        error ("anova.multcompare: LinearModel-backed ANOVA is not supported.");
      endif
      if (isempty (fieldnames (obj.Stats)))
        error ("anova.multcompare: model has no stats to compare.");
      endif
      varargout = cell (1, max (nargout, 1));
      [varargout{:}] = multcompare (obj.Stats, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {} plotDiagnostics (@var{obj})
    ## @deftypefnx {anova} {@var{h} =} plotDiagnostics (@var{obj})
    ## @deftypefnx {anova} {@var{h} =} plotDiagnostics (@var{obj}, @var{name}, @var{value})
    ##
    ## Plot residual diagnostics for an ANOVA object.
    ##
    ## The method creates a four-panel figure containing a Normal Q-Q plot, a
    ## Spread-Location plot, a Residual-Leverage plot, and a Cook's distance
    ## plot.  Diagnostic plots require an @code{anovan}-backed fit or a
    ## @code{LinearModel}-backed object because the fast @code{anova1} and
    ## @code{anova2} backends do not expose residuals and design-matrix
    ## diagnostics.
    ##
    ## Supported name-value arguments are @qcode{'FigureName'} and
    ## @qcode{'Visible'}.
    ##
    ## @seealso{summary}
    ## @end deftypefn
    function h = plotDiagnostics (obj, varargin)
      obj.ensureFit_ ();
      if (! isempty (obj.sourceModel_))
        mdl = obj.sourceModel_;
        h = obj.plotDiagnostics_ (mdl.Residuals.Raw, mdl.Fitted, ...
                                  mdl.Diagnostics.Leverage, ...
                                  mdl.Diagnostics.CooksDistance, ...
                                  mdl.DFE, varargin{:});
        return;
      endif
      if (isempty (obj.Residuals) || isempty (obj.FittedValues) ...
          || isempty (obj.DesignMatrix))
        error (strcat ("anova.plotDiagnostics: diagnostic plots require", ...
                       " an anovan-backed fit."));
      endif

      leverage = obj.leverage_ ();
      if (isfield (obj.Stats, 'CooksD'))
        cooksd = obj.Stats.CooksD;
      else
        cooksd = obj.cooksDistance_ (leverage);
      endif
      h = obj.plotDiagnostics_ (obj.Residuals, obj.FittedValues, ...
                                leverage, cooksd, obj.DFE, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {@var{ypred} =} predict (@var{obj})
    ## @deftypefnx {anova} {@var{ypred} =} predict (@var{obj}, @var{Xnew})
    ##
    ## Predict fitted responses from an ANOVA object.
    ##
    ## With no @var{Xnew}, return fitted values for the training data.  For an
    ## @code{anovan}-backed object, @var{Xnew} must be a numeric design matrix
    ## with one column per coefficient.  For a @code{LinearModel}-backed
    ## object, prediction is delegated to @code{predict} on the stored
    ## @code{LinearModel}, so that table input and categorical encoding remain
    ## owned by @code{LinearModel}.
    ##
    ## @seealso{fitlm, LinearModel}
    ## @end deftypefn
    function ypred = predict (obj, Xnew, varargin)
      obj.ensureFit_ ();
      if (! isempty (obj.sourceModel_))
        if (nargin < 2)
          ypred = predict (obj.sourceModel_);
        else
          ypred = predict (obj.sourceModel_, Xnew, varargin{:});
        endif
        return;
      endif
      if (nargin < 2 || isempty (Xnew))
        ypred = obj.FittedValues;
        return;
      endif
      if (isempty (obj.Coefficients))
        error ("anova.predict: coefficients are unavailable for this backend.");
      endif
      beta = obj.Coefficients(:, 1);
      if (! isnumeric (Xnew) || columns (Xnew) != numel (beta))
        error (strcat ("anova.predict: Xnew must be a numeric design", ...
                       " matrix with %d columns."), numel (beta));
      endif
      ypred = Xnew * beta;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {anova} {@var{ES} =} getEffectSizes (@var{obj})
    ##
    ## Return ANOVA effect-size estimates.
    ##
    ## The returned structure contains fields @code{Source},
    ## @code{EtaSquared}, @code{PartialEtaSquared}, and @code{OmegaSquared}.
    ## For @code{anovan}-backed fits, values are derived from the fitted ANOVA
    ## table.  For @code{LinearModel}-backed fits, model-level values are
    ## computed from @code{SSR}, @code{SSE}, @code{SST}, and @code{MSE}.
    ##
    ## @seealso{summary}
    ## @end deftypefn
    function es = getEffectSizes (obj)
      obj.ensureFit_ ();
      if (strcmp (obj.backend_, 'linearmodel'))
        es = obj.linearModelEffectSizes_ ();
      else
        es = obj.effectSizesFromAtab_ ();
      endif
    endfunction

  endmethods

  methods (Access = private)

    function initLinearModel_ (obj, mdl, varargin)
      if (mod (numel (varargin), 2) != 0)
        error ("anova: name-value pairs must come in pairs.");
      endif
      for idx = 1:2:numel (varargin)
        name = varargin{idx};
        value = varargin{idx + 1};
        if (! ischar (name))
          error ("anova: parameter name must be a character vector.");
        endif
        switch (lower (name))
          case 'alpha'
            obj.Alpha = value;
          case {'display', 'displayopt'}
            obj.Display = value;
          otherwise
            error ("anova: parameter '%s' is not supported.", name);
        endswitch
      endfor
      obj.validateSpec_ ();

      obj.sourceModel_ = mdl;
      obj.backend_ = 'linearmodel';
      obj.nFactors_ = mdl.NumPredictors;
      obj.ModelType = mdl.Formula.Terms;
      obj.VarNames = mdl.VariableNames;
      obj.Y = mdl.Variables{:, mdl.ResponseName};
      obj.GROUP = mdl.Variables(:, mdl.PredictorNames);
      obj.Coefficients = obj.linearModelCoefficients_ (mdl);
      obj.AnovaTable = obj.linearModelAtab_ (mdl);
      obj.Residuals = mdl.Residuals.Raw;
      obj.FittedValues = mdl.Fitted;
      obj.DFE = mdl.DFE;
      obj.MSE = mdl.MSE;
      obj.DesignMatrix = [];
      obj.Stats = obj.linearModelStats_ (mdl);
      obj.fitted_ = true;
      obj.dirty_ = false;
    endfunction

    function validateSpec_ (obj)
      if (! (isnumeric (obj.SSType) && isscalar (obj.SSType) ...
             && any (obj.SSType == [1, 2, 3])))
        error ("anova: SSType must be 1, 2, or 3.");
      endif
      if (! (isnumeric (obj.Alpha) && isscalar (obj.Alpha) ...
             && obj.Alpha > 0 && obj.Alpha < 1))
        error ("anova: Alpha must be a numeric scalar in (0, 1).");
      endif
      if (! ischar (obj.Display) ...
          || ! any (strcmpi (obj.Display, {'on', 'off'})))
        error ("anova: Display must be 'on' or 'off'.");
      endif
      if (ischar (obj.ModelType))
        if (! any (strcmpi (obj.ModelType, ...
                            {'linear', 'interaction', 'full'})))
          error (strcat ("anova: ModelType must be 'linear',", ...
                         " 'interaction', 'full', or a terms matrix."));
        endif
      elseif (! isnumeric (obj.ModelType))
        error ("anova: ModelType must be a string or a numeric terms matrix.");
      endif
      if (! isempty (obj.Continuous) && ! isnumeric (obj.Continuous))
        error ("anova: Continuous must be a numeric index vector.");
      endif
      if (! isempty (obj.Continuous) ...
          && (! isvector (obj.Continuous) ...
              || any (obj.Continuous != fix (obj.Continuous)) ...
              || any (obj.Continuous < 1)))
        error ("anova: Continuous must contain positive integer indices.");
      endif
      if (! isempty (obj.Random) && ! isnumeric (obj.Random))
        error ("anova: Random must be a numeric index vector.");
      endif
      if (! isempty (obj.Random) ...
          && (! isvector (obj.Random) ...
              || any (obj.Random != fix (obj.Random)) ...
              || any (obj.Random < 1)))
        error ("anova: Random must contain positive integer indices.");
      endif
      if (! isempty (obj.Weights) && ! isnumeric (obj.Weights))
        error ("anova: Weights must be numeric.");
      endif
      if (! isempty (obj.reps_) ...
          && ! (isnumeric (obj.reps_) && isscalar (obj.reps_) ...
                && obj.reps_ > 0 && obj.reps_ == fix (obj.reps_)))
        error ("anova: Reps must be a positive integer scalar.");
      endif
    endfunction

    function validateData_ (obj)
      nobs = numel (obj.Y);
      if (! isempty (obj.GROUP))
        obj.validateGroupLength_ (obj.GROUP, nobs);
      endif
      if (! isempty (obj.Weights) && numel (obj.Weights) != nobs)
        error ("anova: Weights must have one value per observation.");
      endif
      if (! isempty (obj.Continuous) && any (obj.Continuous > obj.nFactors_))
        error ("anova: Continuous indices exceed the number of factors.");
      endif
      if (! isempty (obj.Random) && any (obj.Random > obj.nFactors_))
        error ("anova: Random indices exceed the number of factors.");
      endif
      if (isnumeric (obj.ModelType) && ! isempty (obj.ModelType) ...
          && columns (obj.ModelType) != obj.nFactors_)
        error ("anova: terms matrix must have one column per factor.");
      endif
    endfunction

    function validateGroupLength_ (obj, group, nobs)
      if (iscell (group) ...
          && all (cellfun (@(c) isvector (c) || ischar (c), group(:))) ...
          && size (group, 1) == 1)
        for k = 1:numel (group)
          if (numel (group{k}) != nobs)
            error ("anova: GROUP variables must match the number of observations.");
          endif
        endfor
      elseif (isvector (group))
        if (numel (group) != nobs)
          error ("anova: GROUP must match the number of observations.");
        endif
      elseif (rows (group) != nobs)
        error ("anova: GROUP must have one row per observation.");
      endif
    endfunction

    ## Infer the number of factors implied by GROUP / Y.
    ## - Empty GROUP + matrix Y  -> 1 (anova1 matrix form)
    ## - Empty GROUP + vector Y  -> 0 (intercept-only, falls to anovan)
    ## - GROUP is a cell of vectors -> numel (GROUP)
    ## - GROUP is a 2-D numeric / cell matrix -> size (GROUP, 2)
    function nf = countFactors_ (obj)
      if (isempty (obj.GROUP))
        if (ismatrix (obj.Y) && ! isvector (obj.Y))
          nf = 1;
        else
          nf = 0;
        endif
        return;
      endif
      if (iscell (obj.GROUP) ...
          && all (cellfun (@(c) isvector (c) || ischar (c), obj.GROUP(:))) ...
          && size (obj.GROUP, 1) == 1)
        nf = numel (obj.GROUP);
      else
        nf = size (obj.GROUP, 2);
      endif
    endfunction

    ## Backend heuristic.
    ##   anova2  : user passed 'reps' AND Y is a non-vector matrix
    ##             (Y carries the factor structure; reps is required)
    ##   anova1  : 1 factor, no continuous, no weights, SSType == 3
    ##   anovan  : everything else (full generality)
    function selectBackend_ (obj)
      if (! isempty (obj.reps_) && ismatrix (obj.Y) ...
          && ! isvector (obj.Y) && ! any (isnan (obj.Y(:))) ...
          && isempty (obj.Continuous) && isempty (obj.Weights))
        obj.backend_ = 'anova2';
      elseif (obj.nFactors_ == 1 && isempty (obj.Continuous) ...
              && isempty (obj.Weights) && obj.SSType == 3)
        obj.backend_ = 'anova1';
      else
        obj.backend_ = 'anovan';
      endif
    endfunction

    ## Lazy refit guard: only fits when never-fit or spec changed.
    function ensureFit_ (obj)
      if (! obj.fitted_ || obj.dirty_)
        obj.selectBackend_ ();
        obj.fit_ ();
      endif
    endfunction

    ## Dispatch to the selected backend; populate result properties.
    function fit_ (obj)
      switch (obj.backend_)
        case 'anova1'
          obj.fitAnova1_ ();
        case 'anova2'
          obj.fitAnova2_ ();
        case 'anovan'
          obj.fitAnovan_ ();
        case 'linearmodel'
          ## Already populated from the supplied LinearModel object.
      endswitch
      obj.fitted_ = true;
      obj.dirty_  = false;
    endfunction

    function fitAnova1_ (obj)
      if (isvector (obj.Y))
        [~, atab, stats] = anova1 (obj.Y, obj.GROUP, obj.Display);
      else
        [~, atab, stats] = anova1 (obj.Y, [], obj.Display);
      endif
      obj.AnovaTable = atab;
      obj.Stats      = stats;
      obj.DFE        = stats.df;
      obj.MSE        = stats.s ^ 2;     ## anova1 reports sqrt(MSE) as s
      ## Coefficients / Residuals / DesignMatrix / FittedValues are not
      ## exposed by anova1's stats; they remain at their empty defaults.
    endfunction

    function fitAnova2_ (obj)
      modelarg = 'interaction';
      if (ischar (obj.ModelType))
        modelarg = obj.ModelType;
      endif
      [~, atab, stats] = anova2 (obj.Y, obj.reps_, obj.Display, modelarg);
      obj.AnovaTable = atab;
      obj.Stats      = stats;
      obj.DFE        = stats.df;
      obj.MSE        = stats.sigmasq;
      ## anova2 does not return coeffs / resid / X.
    endfunction

    function fitAnovan_ (obj)
      ## Unroll a matrix Y with no GROUP into a vector + synthetic
      ## column index — needed for the anova1-matrix-form fallback.
      if (isempty (obj.GROUP) && ! isvector (obj.Y))
        [n, m] = size (obj.Y);
        y_vec  = obj.Y(:);
        g_vec  = reshape (repmat ((1:m), n, 1), [], 1);
        group_arg = {g_vec};
      else
        y_vec = obj.Y(:);
        group_arg = obj.GROUP;
        if (isempty (group_arg))
          group_arg = {};
        endif
      endif
      [~, atab, stats] = anovan (y_vec, group_arg, ...
                                 obj.buildAnovanArgs_(){:});
      obj.AnovaTable = atab;
      obj.Stats      = stats;
      if (isfield (stats, 'coeffs'))
        obj.Coefficients = stats.coeffs;
      endif
      if (isfield (stats, 'resid'))
        obj.Residuals = stats.resid;
      endif
      if (isfield (stats, 'X'))
        obj.DesignMatrix = stats.X;
      endif
      if (isfield (stats, 'dfe'))
        obj.DFE = stats.dfe;
      endif
      if (isfield (stats, 'mse'))
        obj.MSE = stats.mse;
      endif
      if (! isempty (obj.DesignMatrix) && ! isempty (obj.Coefficients))
        obj.FittedValues = full (obj.DesignMatrix) * obj.Coefficients(:, 1);
      endif
    endfunction

    function s = sstypeLabel_ (obj)
      switch (obj.SSType)
        case 1; s = "I";
        case 2; s = "II";
        otherwise; s = "III";
      endswitch
    endfunction

    function printAtab_ (obj, atab)
      [nrows, ncols] = size (atab);
      col_w = max (12, ceil (80 / max (ncols, 1)));
      for j = 1:ncols
        fprintf ("%-*s", col_w, char (atab{1, j}));
      endfor
      fprintf ("\n%s\n", repmat ("-", 1, col_w * ncols));
      for i = 2:nrows
        for j = 1:ncols
          v = atab{i, j};
          if (ischar (v))
            fprintf ("%-*s", col_w, v);
          elseif (isnumeric (v) && ! isempty (v) && isscalar (v))
            if (isnan (v))
              fprintf ("%-*s", col_w, "NaN");
            elseif (v == fix (v) && abs (v) < 1e6)
              fprintf ("%-*d", col_w, v);
            else
              fprintf ("%-*.*g", col_w, 5, v);
            endif
          else
            fprintf ("%-*s", col_w, "");
          endif
        endfor
        fprintf ("\n");
      endfor
    endfunction

    function nv = buildAnovanArgs_ (obj)
      nv = {'display', obj.Display, 'sstype', obj.SSType, 'alpha', obj.Alpha};
      if (ischar (obj.ModelType) || isnumeric (obj.ModelType))
        if (! (isnumeric (obj.ModelType) && isempty (obj.ModelType)))
          nv = [nv, {'model', obj.ModelType}];
        endif
      endif
      if (! isempty (obj.VarNames))
        nv = [nv, {'varnames', obj.VarNames}];
      endif
      if (! isempty (obj.Continuous))
        nv = [nv, {'continuous', obj.Continuous}];
      endif
      if (! isempty (obj.Random))
        nv = [nv, {'random', obj.Random}];
      endif
      if (! isempty (obj.Weights))
        nv = [nv, {'weights', obj.Weights}];
      endif
      if (! isempty (obj.Contrasts))
        nv = [nv, {'contrasts', obj.Contrasts}];
      endif
    endfunction

    function h = leverage_ (obj)
      X = full (obj.DesignMatrix);
      Q = qr (X, 0);
      h = sum (Q .^ 2, 2);
    endfunction

    function D = cooksDistance_ (obj, leverage)
      p = max (columns (obj.DesignMatrix), 1);
      D = (obj.Residuals .^ 2 ./ max (p * obj.MSE, eps)) ...
          .* leverage ./ max ((1 - leverage) .^ 2, eps);
    endfunction

    function h = plotDiagnostics_ (obj, residuals, fitted, leverage, ...
                                   cooksd, dfe, varargin)
      if (isempty (residuals) || isempty (fitted) || isempty (leverage) ...
          || isempty (cooksd))
        error ("anova.plotDiagnostics: diagnostic inputs must be non-empty.");
      endif

      residuals = residuals(:);
      fitted    = fitted(:);
      leverage  = leverage(:);
      cooksd    = cooksd(:);
      n = numel (residuals);
      if (numel (fitted) != n || numel (leverage) != n ...
          || numel (cooksd) != n)
        error ("anova.plotDiagnostics: diagnostic inputs must match.");
      endif

      fig_name = "Diagnostic Plots: Model Residuals";
      visible = "on";
      if (mod (numel (varargin), 2) != 0)
        error ("anova.plotDiagnostics: name-value pairs must come in pairs.");
      endif
      for k = 1:2:numel (varargin)
        switch (lower (varargin{k}))
          case "figurename"
            fig_name = varargin{k + 1};
          case "visible"
            visible = varargin{k + 1};
          otherwise
            error ("anova.plotDiagnostics: unknown option '%s'.", varargin{k});
        endswitch
      endfor

      mse = sum (residuals .^ 2) / max (dfe, 1);
      t = residuals ./ sqrt (mse * max (1 - leverage, eps));
      [~, DI] = sort (cooksd, "descend");
      nk = min (4, n);

      h = figure ("Name", fig_name, "Visible", visible);

      subplot (2, 2, 1);
      x = ((1:n)' - 0.5) / n;
      [ts, I] = sort (t);
      q = norminv (x);
      plot (q, ts, "ok", "markersize", 3);
      box off;
      grid on;
      xlabel ("Theoretical quantiles");
      ylabel ("Studentized residuals");
      title ("Normal Q-Q Plot");
      arrayfun (@(i) text (q(I == DI(i)), t(DI(i)), ...
                           sprintf ("  %u", DI(i))), 1:nk);
      iqr = [0.25; 0.75];
      yl = quantile (t, iqr, 1, 6);
      xl = norminv (iqr);
      slope = diff (yl) / diff (xl);
      int = yl(1) - slope * xl(1);
      ax1_xlim = get (gca, "XLim");
      hold on;
      plot (ax1_xlim, slope * ax1_xlim + int, "k-");
      hold off;
      set (gca, "Xlim", ax1_xlim);

      subplot (2, 2, 2);
      plot (fitted, sqrt (abs (t)), "ko", "markersize", 3);
      box off;
      xlabel ("Fitted values");
      ylabel ("sqrt ( | Studentized residuals | )");
      title ("Spread-Location Plot");
      ax2_xlim = get (gca, "XLim");
      hold on;
      plot (ax2_xlim, ones (1, 2) * sqrt (2), "k:");
      plot (ax2_xlim, ones (1, 2) * sqrt (3), "k-.");
      plot (ax2_xlim, ones (1, 2) * sqrt (4), "k--");
      hold off;
      arrayfun (@(i) text (fitted(DI(i)), sqrt (abs (t(DI(i)))), ...
                           sprintf ("  %u", DI(i))), 1:nk);
      xlim (ax2_xlim);

      subplot (2, 2, 3);
      plot (leverage, t, "ko", "markersize", 3);
      box off;
      xlabel ("Leverage");
      ylabel ("Studentized residuals");
      title ("Residual-Leverage Plot");
      ax3_xlim = get (gca, "XLim");
      ax3_ylim = get (gca, "YLim");
      hold on;
      plot (ax3_xlim, zeros (1, 2), "k-");
      hold off;
      arrayfun (@(i) text (leverage(DI(i)), t(DI(i)), ...
                           sprintf ("  %u", DI(i))), 1:nk);
      set (gca, "ygrid", "on");
      xlim (ax3_xlim);
      ylim (ax3_ylim);

      subplot (2, 2, 4);
      stem (cooksd, "ko", "markersize", 3);
      box off;
      xlabel ("Obs. number");
      ylabel ("Cook's distance");
      title ("Cook's Distance Stem Plot");
      xlim ([0, n]);
      ax4_xlim = get (gca, "XLim");
      ax4_ylim = get (gca, "YLim");
      hold on;
      plot (ax4_xlim, ones (1, 2) * 4 / max (dfe, eps), "k:");
      plot (ax4_xlim, ones (1, 2) * 0.5, "k-.");
      plot (ax4_xlim, ones (1, 2), "k--");
      hold off;
      arrayfun (@(i) text (DI(i), cooksd(DI(i)), ...
                           sprintf ("  %u", DI(i))), 1:nk);
      xlim (ax4_xlim);
      ylim (ax4_ylim);

      set (findall (gcf, "-property", "FontSize"), "FontSize", 7);
    endfunction

    function coeffs = linearModelCoefficients_ (obj, mdl)
      est = mdl.Coefficients.Estimate;
      se = mdl.Coefficients.SE;
      tstat = mdl.Coefficients.tStat;
      pval = mdl.Coefficients.pValue;
      crit = tinv (1 - obj.Alpha / 2, mdl.DFE);
      coeffs = [est, se, est - crit * se, est + crit * se, tstat, pval];
    endfunction

    function atab = linearModelAtab_ (obj, mdl)
      has_intercept = isfield (mdl.Formula, 'HasIntercept') ...
                      && mdl.Formula.HasIntercept;
      df_model = mdl.NumEstimatedCoefficients - double (has_intercept);
      df_model = max (df_model, 0);
      if (df_model > 0)
        ms_model = mdl.SSR / df_model;
      else
        ms_model = NaN;
      endif
      atab = {'Source', 'SS', 'df', 'MS', 'F', 'Prob>F'; ...
              'Model', mdl.SSR, df_model, ms_model, ...
              mdl.ModelFitVsNullModel.Fstat, mdl.ModelFitVsNullModel.Pvalue; ...
              'Error', mdl.SSE, mdl.DFE, mdl.MSE, '', ''; ...
              'Total', mdl.SST, mdl.NumObservations - double (has_intercept), ...
              '', '', ''};
    endfunction

    function stats = linearModelStats_ (obj, mdl)
      stats = struct ("source", "linearmodel", ...
                      "resid", mdl.Residuals.Raw, ...
                      "coeffs", obj.Coefficients, ...
                      "dfe", mdl.DFE, ...
                      "mse", mdl.MSE, ...
                      "vcov", mdl.CoefficientCovariance, ...
                      "CooksD", mdl.Diagnostics.CooksDistance, ...
                      "alpha", obj.Alpha, ...
                      "varnames", {mdl.VariableNames}, ...
                      "coeffnames", {mdl.CoefficientNames});
    endfunction

    function es = linearModelEffectSizes_ (obj)
      mdl = obj.sourceModel_;
      eta = mdl.SSR / max (mdl.SST, eps);
      partial_eta = mdl.SSR / max (mdl.SSR + mdl.SSE, eps);
      has_intercept = isfield (mdl.Formula, 'HasIntercept') ...
                      && mdl.Formula.HasIntercept;
      df_model = max (mdl.NumEstimatedCoefficients - double (has_intercept), 0);
      omega = (mdl.SSR - df_model * mdl.MSE) / max (mdl.SST + mdl.MSE, eps);
      es = struct ();
      es.Source = {'Model'};
      es.EtaSquared = eta;
      es.PartialEtaSquared = partial_eta;
      es.OmegaSquared = omega;
    endfunction

    function es = effectSizesFromAtab_ (obj)
      atab = obj.AnovaTable;
      if (isempty (atab))
        error ("anova.getEffectSizes: model has no ANOVA table.");
      endif

      source_col = obj.findAtabColumn_ (atab, {'Source'});
      ss_col = obj.findAtabColumn_ (atab, {'SS', 'Sum Sq.', 'Sum Sq'});
      df_col = obj.findAtabColumn_ (atab, {'df', 'd.f.'});
      sources = {};
      ss = [];
      df = [];
      sse = [];
      sst = [];
      for r = 2:rows (atab)
        name = atab{r, source_col};
        if (! ischar (name))
          continue;
        endif
        val_ss = atab{r, ss_col};
        val_df = atab{r, df_col};
        if (strcmpi (name, 'Error'))
          sse = val_ss;
        elseif (strcmpi (name, 'Total'))
          sst = val_ss;
        elseif (isnumeric (val_ss) && isscalar (val_ss))
          sources{end + 1} = name;
          ss(end + 1, 1) = val_ss;
          df(end + 1, 1) = val_df;
        endif
      endfor
      if (isempty (sst))
        sst = sum (ss) + ifelse (isempty (sse), 0, sse);
      endif
      if (isempty (sse))
        sse = max (sst - sum (ss), 0);
      endif
      eta = ss ./ max (sst, eps);
      partial_eta = ss ./ max (ss + sse, eps);
      omega = (ss - df .* obj.MSE) ./ max (sst + obj.MSE, eps);
      es = struct ();
      es.Source = sources;
      es.EtaSquared = eta;
      es.PartialEtaSquared = partial_eta;
      es.OmegaSquared = omega;
    endfunction

    function idx = findAtabColumn_ (obj, atab, names)
      idx = [];
      for k = 1:numel (names)
        hit = find (strcmpi (atab(1, :), names{k}), 1);
        if (! isempty (hit))
          idx = hit;
          return;
        endif
      endfor
      error ("anova: ANOVA table is missing the '%s' column.", names{1});
    endfunction

  endmethods

endclassdef

%!demo
%! ## One-way ANOVA with a formatted summary
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (y, g, "SSType", 2);
%! summary (a);

%!demo
%! ## Post-hoc multiple comparisons
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (y, g, "SSType", 2);
%! C = multcompare (a, "display", "off")

%!demo
%! ## Diagnostic plots for an anovan-backed fit
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (y, g, "SSType", 2);
%! plotDiagnostics (a);

## --- BISTs ---------------------------------------------------------------

## Basic construction: vector Y + single grouping vector
%!test
%! y = [1; 2; 3; 4; 5; 6];
%! g = [1; 1; 2; 2; 3; 3];
%! a = anova (y, g);
%! assert (class (a), "anova");
%! assert (a.Y, y);
%! assert (a.GROUP, g);

## Basic construction: matrix Y, no GROUP
%!test
%! y = magic (4);
%! a = anova (y);
%! assert (class (a), "anova");
%! assert (a.Y, y);
%! assert (isempty (a.GROUP));

## Basic construction: cell of two grouping vectors
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = repmat ([1;1;2;2], 3, 1);
%! a = anova (y, {g1, g2});
%! assert (class (a), "anova");
%! assert (! isempty (strfind (evalc ("disp (a)"), "nFactors : 2")));

## Property defaults
%!test
%! a = anova ([1;2;3;4], [1;1;2;2]);
%! assert (a.ModelType, 'linear');
%! assert (a.SSType, 3);
%! assert (a.Alpha, 0.05);
%! assert (a.Display, 'off');
%! assert (isequal (a.VarNames, {}));
%! assert (isequal (a.Contrasts, {}));
%! assert (isempty (a.Continuous));
%! assert (isempty (a.Random));
%! assert (isempty (a.Weights));

## Result properties stay empty before fit
%!test
%! a = anova ([1;2;3;4], [1;1;2;2]);
%! assert (isempty (a.Coefficients));
%! assert (isempty (a.Residuals));
%! assert (isempty (a.FittedValues));
%! assert (isempty (a.DFE));
%! assert (isempty (a.MSE));
%! assert (isempty (a.DesignMatrix));

## Name-value parsing: SSType, ModelType, Alpha, Display
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = repmat ([1;1;2;2], 3, 1);
%! a = anova (y, {g1, g2}, 'SSType', 2, 'Model', 'full', ...
%!            'Alpha', 0.01, 'Display', 'on');
%! assert (a.SSType, 2);
%! assert (a.ModelType, 'full');
%! assert (a.Alpha, 0.01);
%! assert (a.Display, 'on');

## Name-value parsing: case-insensitive names, displayopt alias
%!test
%! a = anova ([1;2;3;4], [1;1;2;2], 'sstype', 1, 'displayopt', 'on');
%! assert (a.SSType, 1);
%! assert (a.Display, 'on');

## Backend selection: one-way default -> anova1
%!test
%! a = anova ([1;2;3;4;5;6], [1;1;2;2;3;3]);
%! assert (! isempty (strfind (evalc ("disp (a)"), "backend  : anova1")));

## Backend selection: one-way matrix-Y form -> anova1
%!test
%! a = anova (magic (4));
%! str = evalc ("disp (a)");
%! assert (! isempty (strfind (str, "nFactors : 1")));
%! assert (! isempty (strfind (str, "backend  : anova1")));

## Backend selection: matrix Y + explicit 'reps' -> anova2
%!test
%! y = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!      6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (y, [], 'reps', 3);
%! assert (! isempty (strfind (evalc ("disp (a)"), "backend  : anova2")));

## Backend selection: two-factor cell groups without reps -> anovan
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = repmat ([1;1;2;2], 3, 1);
%! a = anova (y, {g1, g2});
%! str = evalc ("disp (a)");
%! assert (! isempty (strfind (str, "nFactors : 2")));
%! assert (! isempty (strfind (str, "backend  : anovan")));

## Backend selection: three factors -> anovan
%!test
%! y = (1:24)';
%! g1 = repmat ([1;2], 12, 1);
%! g2 = repmat ([1;1;2;2], 6, 1);
%! g3 = repmat ([1;1;1;1;2;2;2;2], 3, 1);
%! a = anova (y, {g1, g2, g3});
%! str = evalc ("disp (a)");
%! assert (! isempty (strfind (str, "nFactors : 3")));
%! assert (! isempty (strfind (str, "backend  : anovan")));

## Backend selection: SSType != 3 with 1 factor falls through to anovan
%!test
%! a = anova ([1;2;3;4;5;6], [1;1;2;2;3;3], 'SSType', 2);
%! assert (! isempty (strfind (evalc ("disp (a)"), "backend  : anovan")));

## Backend selection: continuous predictors force anovan
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = (1:12)';                              ## continuous
%! a = anova (y, {g1, g2}, 'Continuous', 2);
%! assert (! isempty (strfind (evalc ("disp (a)"), "backend  : anovan")));

## Backend selection: NaN in a vector response falls through to anovan
%!test
%! y = [1; 2; 3; NaN; 5; 6; 7; 8];
%! g1 = [1; 2; 1; 2; 1; 2; 1; 2];
%! g2 = [1; 1; 2; 2; 1; 1; 2; 2];
%! a = anova (y, {g1, g2});
%! assert (! isempty (strfind (evalc ("disp (a)"), "backend  : anovan")));

## Backend selection: weights force anovan
%!test
%! a = anova ([1;2;3;4;5;6], [1;1;2;2;3;3], 'Weights', ones (6, 1));
%! assert (! isempty (strfind (evalc ("disp (a)"), "backend  : anovan")));

## Invalid input: no args
%!error <anova: too few input arguments.> anova ()

## Invalid input: non-numeric Y
%!error <anova: Y must be a non-empty numeric array.> anova ('abc')

## Invalid input: empty Y
%!error <anova: Y must be a non-empty numeric array.> anova ([])

## Invalid input: odd number of name-value arguments
%!error <anova: name-value pairs must come in pairs.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'SSType')

## Invalid input: unknown parameter name
%!error <anova: parameter 'bogus' is not supported.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'bogus', 1)

## Invalid input: SSType out of range
%!error <anova: SSType must be 1, 2, or 3.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'SSType', 5)

## Invalid input: Alpha out of range
%!error <anova: Alpha must be a numeric scalar in .0, 1..> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'Alpha', 2)

## Invalid input: bad Display value
%!error <anova: Display must be 'on' or 'off'.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'Display', 'maybe')

## Invalid input: parameter name must be a string
%!error <anova: parameter name must be a character vector.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 1, 2)

## Invalid input: bad Reps value
%!error <anova: Reps must be a positive integer scalar.> ...
%!  anova (magic (4), [], 'reps', -2)

## Invalid input: non-integer Reps value
%!error <anova: Reps must be a positive integer scalar.> ...
%!  anova (magic (4), [], 'reps', 1.5)

## Invalid input: bad ModelType
%!error <anova: ModelType must be> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'Model', 'quadratic')

## Invalid input: terms matrix width must match factor count
%!error <anova: terms matrix must have one column per factor.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'Model', [1 0])

## Invalid input: GROUP length mismatch
%!error <anova: GROUP must match the number of observations.> ...
%!  anova ([1;2;3;4], [1;1;2])

## Invalid input: cell GROUP variable length mismatch
%!error <anova: GROUP variables must match the number of observations.> ...
%!  anova ([1;2;3;4], {[1;1;2], [1;2;1;2]})

## Invalid input: Weights length mismatch
%!error <anova: Weights must have one value per observation.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'Weights', [1;1;1])

## Invalid input: Continuous indices must be positive integers
%!error <anova: Continuous must contain positive integer indices.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'Continuous', 1.5)

## Invalid input: Continuous index out of range
%!error <anova: Continuous indices exceed the number of factors.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'Continuous', 2)

## Invalid input: Random indices must be positive integers
%!error <anova: Random must contain positive integer indices.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'Random', 0)

## Invalid input: Random index out of range
%!error <anova: Random indices exceed the number of factors.> ...
%!  anova ([1;2;3;4], [1;1;2;2], 'Random', 2)

## --- Week 2: fit delegation smoke tests --------------------------------

## fit_(): one-way fixture (anova1 backend; falls back to anovan on
## Octave 10.3 where iscategorical is missing — both paths populate
## the unified result surface)
%!test
%! y = [1; 2; 3; 4; 5; 6];
%! g = [1; 1; 2; 2; 3; 3];
%! a = anova (y, g);
%! a.fit ();
%! assert (! isempty (a.AnovaTable));
%! assert (isstruct (a.Stats));
%! assert (! isempty (a.DFE));
%! assert (! isempty (a.MSE));

## fit_(): two-way balanced fixture (anova2 backend, popcorn data)
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (popcorn, [], 'reps', 3);
%! assert (! isempty (strfind (evalc ("disp (a)"), "backend  : anova2")));
%! a.fit ();
%! assert (! isempty (a.AnovaTable));
%! assert (isfield (a.Stats, 'sigmasq'));
%! assert (a.MSE, a.Stats.sigmasq);

## fit_(): N-way fixture (anovan backend, three factors)
%!test
%! y = (1:24)';
%! g1 = repmat ([1;2], 12, 1);
%! g2 = repmat ([1;1;2;2], 6, 1);
%! g3 = repmat ([1;1;1;1;2;2;2;2], 3, 1);
%! a = anova (y, {g1, g2, g3});
%! assert (! isempty (strfind (evalc ("disp (a)"), "backend  : anovan")));
%! a.fit ();
%! assert (! isempty (a.AnovaTable));
%! assert (! isempty (a.Coefficients));
%! assert (! isempty (a.Residuals));
%! assert (! isempty (a.DesignMatrix));
%! assert (rows (a.Residuals), numel (y));

## fit_(): FittedValues = DesignMatrix * Coefficients(:,1) for anovan
%!test
%! y  = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g  = [1;1;1;2;2;2;3;3;3];
%! a  = anova (y, g, 'SSType', 2);
%! a.fit ();
%! assert (numel (a.FittedValues), numel (y));
%! assert (a.FittedValues + a.Residuals, y, 1e-9);

## ensureFit_(): fit() is idempotent (second call does nothing)
%!test
%! a = anova ((1:9)', [1;1;1;2;2;2;3;3;3]);
%! a.fit ();
%! first_table = a.AnovaTable;
%! a.fit ();
%! assert (isequal (a.AnovaTable, first_table));

## buildAnovanArgs_(): SSType and Alpha are forwarded to anovan
%!test
%! y = (1:12)';
%! g = repmat ([1;2;3], 4, 1);
%! a = anova (y, {g}, 'SSType', 2, 'Alpha', 0.10);
%! a.fit ();
%! assert (! isempty (strfind (evalc ("disp (a)"), "backend  : anovan")));
%! assert (a.Stats.alpha, 0.10, 1e-12);

## --- Numeric references -------------------------------------------------

## One-way ANOVA: reference values match R's aov(y ~ factor(g)).
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (y, g, 'SSType', 2);
%! a.fit ();
%! T = a.AnovaTable;
%! assert (T{2, 2}, 126, 1e-12);       # group SS
%! assert (T{3, 2}, 6, 1e-12);         # error SS
%! assert (T{4, 2}, 132, 1e-12);       # total SS
%! assert (T{2, 3}, 2);
%! assert (T{3, 3}, 6);
%! assert (T{2, 4}, 63, 1e-12);
%! assert (T{3, 4}, 1, 1e-12);
%! assert (T{2, 6}, 63, 1e-12);
%! assert (T{2, 7}, 9.3914e-05, 1e-9);
%! assert (a.MSE, 1, 1e-12);
%! assert (a.DFE, 6);

## Fitted values and residuals: reference values from one-way cell means.
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (y, g, 'SSType', 2);
%! assert (predict (a), [2; 2; 2; 5; 5; 5; 11; 11; 11], 1e-12);
%! assert (a.Residuals, [-1; 0; 1; -1; 0; 1; -1; 0; 1], 1e-12);

## Effect sizes: reference formulas eta2=SS/SST and omega2=(SS-df*MSE)/(SST+MSE).
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! es = getEffectSizes (anova (y, g, 'SSType', 2));
%! assert (es.Source, {'X1'});
%! assert (es.EtaSquared, 126 / 132, 1e-12);
%! assert (es.PartialEtaSquared, 126 / (126 + 6), 1e-12);
%! assert (es.OmegaSquared, 124 / 133, 1e-12);

## Two-way balanced ANOVA: popcorn values match the anova2 documentation example.
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (popcorn, [], 'reps', 3);
%! a.fit ();
%! T = a.AnovaTable;
%! assert (T{2, 2}, 15.75, 1e-12);
%! assert (T{3, 2}, 4.5, 1e-12);
%! assert (T{4, 2}, 1.75, 1e-12);
%! assert (T{5, 2}, 22, 1e-12);
%! assert (T{2, 5}, 63, 1e-12);
%! assert (T{3, 5}, 36, 1e-12);
%! assert (a.MSE, 0.125, 1e-12);

## Post-hoc comparisons: reference values checked against R's TukeyHSD output.
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (y, g, 'SSType', 2);
%! C = multcompare (a, 'display', 'off');
%! assert (C(:, 1:2), [1 2; 1 3; 2 3]);
%! assert (C(:, 4), [-3; -9; -6], 1e-12);
%! assert (C(1, 6), 0.010402, 1e-6);
%! assert (C(2, 6), 9.9474e-05, 1e-9);
%! assert (C(3, 6), 0.00064995, 1e-8);

## --- Week 3: summary / disp ---------------------------------------------

## summary(): runs ensureFit_ and prints a table
%!test
%! a = anova ((1:9)', [1;1;1;2;2;2;3;3;3], 'SSType', 2);
%! str = evalc ('summary (a)');
%! assert (! isempty (strfind (str, 'ANOVA TABLE')));
%! assert (! isempty (strfind (str, 'backend = anovan')));

## summary(): includes MSE / DFE / Alpha line
%!test
%! a = anova ((1:9)', [1;1;1;2;2;2;3;3;3], 'SSType', 2, 'Alpha', 0.10);
%! str = evalc ('summary (a)');
%! assert (! isempty (strfind (str, 'Alpha: 0.1')));

## disp(): one-line overview of key fields
%!test
%! a = anova ((1:9)', [1;1;1;2;2;2;3;3;3], 'SSType', 2);
%! str = evalc ('disp (a)');
%! assert (! isempty (strfind (str, 'anova object')));
%! assert (! isempty (strfind (str, 'backend')));
%! assert (! isempty (strfind (str, 'SSType')));

## summary(): SSType label appears in the header (Type I / II / III)
%!test
%! a1 = anova ((1:9)', [1;1;1;2;2;2;3;3;3], 'SSType', 1);
%! a2 = anova ((1:9)', [1;1;1;2;2;2;3;3;3], 'SSType', 2);
%! assert (! isempty (strfind (evalc ('summary (a1)'), 'Type I sums')));
%! assert (! isempty (strfind (evalc ('summary (a2)'), 'Type II sums')));

## --- Week 4: multcompare pass-through ----------------------------------

## multcompare(): anovan backend returns a pairwise comparison matrix
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (y, g, 'SSType', 2);
%! C = multcompare (a, 'display', 'off');
%! assert (! isempty (C));
%! assert (size (C, 1), 3);              ## 3 pairwise comparisons for 3 groups
%! assert (size (C, 2) >= 6);            ## at least i,j,diff,lo,hi,p

## multcompare(): two-way balanced via anova2 backend
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (popcorn, [], 'reps', 3);
%! C = multcompare (a, 'display', 'off', 'estimate', 'column');
%! assert (! isempty (C));
%! assert (size (C, 2) >= 6);

## multcompare(): runs after a fresh construction (triggers ensureFit_)
%!test
%! a = anova ((1:6)', [1;1;2;2;3;3], 'SSType', 2);
%! C = multcompare (a, 'display', 'off');
%! assert (! isempty (C));

## --- Week 5: diagnostic plots ------------------------------------------

## plotDiagnostics(): anovan-backed fit creates the four-panel figure
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%!   g = [1;1;1;2;2;2;3;3;3];
%!   a = anova (y, g, 'SSType', 2);
%!   h = plotDiagnostics (a, 'Visible', 'off');
%!   assert (ishghandle (h));
%!   assert (numel (findall (h, 'type', 'axes')), 4);
%! unwind_protect_cleanup
%!   close (h);
%!   close (hf);
%! end_unwind_protect

## plotDiagnostics(): fast-path backends report the missing diagnostics
%!error <diagnostic plots require>
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! plotDiagnostics (anova (popcorn, [], 'reps', 3));

## plotDiagnostics(): rejects odd name-value arguments
%!error <name-value pairs must come in pairs>
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! plotDiagnostics (anova (y, g, 'SSType', 2), 'Visible');

## plotDiagnostics(): rejects unknown options
%!error <unknown option>
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! plotDiagnostics (anova (y, g, 'SSType', 2), 'BadOption', true);

## --- Week 6: predict / effect sizes / LinearModel bridge ----------------

## predict(): no Xnew returns fitted values
%!test
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! a = anova (y, g, 'SSType', 2);
%! assert (predict (a), a.FittedValues, 1e-9);

## predict(): accepts an explicit design matrix for anovan-backed fits
%!test
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! a = anova (y, g, 'SSType', 2);
%! a.fit ();
%! assert (predict (a, full (a.DesignMatrix)), a.FittedValues, 1e-9);

## predict(): rejects design matrices with the wrong width
%!error <Xnew must be a numeric design matrix>
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! predict (anova (y, g, 'SSType', 2), ones (2, 2));

## getEffectSizes(): anovan-backed fits report effect-size vectors
%!test
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! a = anova (y, g, 'SSType', 2);
%! es = getEffectSizes (a);
%! assert (iscell (es.Source));
%! assert (numel (es.EtaSquared), numel (es.Source));
%! assert (all (isfinite (es.PartialEtaSquared)));

## anova(LinearModel): builds from Avanish's LinearModel surface
%!test
%! X = [1; 2; 3; 4; 5];
%! y = [2; 4; 5; 4; 5];
%! mdl = fitlm (X, y);
%! a = anova (mdl);
%! assert (! isempty (strfind (evalc ("disp (a)"), "backend  : linearmodel")));
%! assert (! isempty (a.AnovaTable));
%! assert (predict (a), mdl.Fitted, 1e-9);
%! es = getEffectSizes (a);
%! assert (es.Source, {'Model'});
%! assert (es.EtaSquared, mdl.SSR / mdl.SST, 1e-9);

## anova(LinearModel): plotDiagnostics uses LinearModel diagnostics
%!test
%! X = [1; 2; 3; 4; 5];
%! y = [2; 4; 5; 4; 5];
%! mdl = fitlm (X, y);
%! a = anova (mdl);
%! h = plotDiagnostics (a, 'Visible', 'off');
%! unwind_protect
%!   assert (ishghandle (h));
%!   assert (numel (findall (h, 'type', 'axes')), 4);
%! unwind_protect_cleanup
%!   close (h);
%! end_unwind_protect
