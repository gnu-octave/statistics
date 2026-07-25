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

classdef anova
  ## -*- texinfo -*-
  ## @deftp {statistics} anova
  ##
  ## Object-oriented interface for analysis of variance.
  ##
  ## The @code{anova} class provides a MATLAB-compatible object interface for
  ## analysis of variance.  It stores factors, response data, model
  ## specification, and fitted results in one object.  The class chooses the
  ## narrowest compatible backend, delegates the numeric computation to the
  ## existing ANOVA functions, and exposes common follow-up operations such as
  ## @code{stats}, @code{groupmeans}, @code{boxchart},
  ## @code{plotComparisons}, @code{varianceComponent}, and
  ## @code{multcompare}.
  ##
  ## Models are fitted lazily.  Methods that need fitted results call
  ## @code{fit} internally when necessary, so users may construct an object and
  ## immediately call inspection or post-hoc methods.
  ##
  ## @seealso{anova1, anova2, anovan, multcompare, fitlm, LinearModel}
  ## @end deftp

  properties (GetAccess = public, SetAccess = private)
    ## Data

    ## -*- texinfo -*-
    ## @deftp {anova} {property} Y
    ##
    ## Response data
    ##
    ## Numeric response vector (or matrix, for the one-way column form) used to
    ## fit the ANOVA model.  This property is read-only.
    ##
    ## @end deftp
    Y

    ## -*- texinfo -*-
    ## @deftp {anova} {property} Factors
    ##
    ## Factor data
    ##
    ## Table whose variables are the factors used to fit the ANOVA model, with
    ## one column per factor named after @code{FactorNames} and one row per
    ## observation.  This property is read-only.
    ##
    ## @end deftp
    Factors         = [];

    ## -*- texinfo -*-
    ## @deftp {anova} {property} Formula
    ##
    ## Model formula
    ##
    ## Character vector describing the fitted ANOVA model formula.  (MATLAB
    ## returns a formula object; Octave returns the equivalent character
    ## vector.)  This property is read-only.
    ##
    ## @end deftp
    Formula         = '';

    ## -*- texinfo -*-
    ## @deftp {anova} {property} FactorNames
    ##
    ## Factor names
    ##
    ## Cell array of character vectors used as factor names in the fitted
    ## ANOVA table.  This property is read-only.
    ##
    ## @end deftp
    FactorNames     = {};

    ## -*- texinfo -*-
    ## @deftp {anova} {property} ExpandedFactorNames
    ##
    ## Coefficient names
    ##
    ## Cell array of character vectors naming the model coefficients when the
    ## selected backend exposes them, otherwise an empty cell array.  This
    ## property is read-only.
    ##
    ## @end deftp
    ExpandedFactorNames = {};

    ## -*- texinfo -*-
    ## @deftp {anova} {property} SumOfSquaresType
    ##
    ## Sum-of-squares type
    ##
    ## Character vector selecting @qcode{'one'}, @qcode{'two'},
    ## @qcode{'three'}, or @qcode{'hierarchical'} sums of squares.  This
    ## property is read-only.
    ##
    ## @end deftp
    SumOfSquaresType = 'three';

    ## -*- texinfo -*-
    ## @deftp {anova} {property} RandomFactors
    ##
    ## Random factors
    ##
    ## Positive integer indices of random factors.  This property is
    ## read-only.
    ##
    ## @end deftp
    RandomFactors   = [];

    ## -*- texinfo -*-
    ## @deftp {anova} {property} CategoricalFactors
    ##
    ## Categorical factors
    ##
    ## Character vector @qcode{'all'} or positive integer indices of factors
    ## treated as categorical.  This property is read-only.
    ##
    ## @end deftp
    CategoricalFactors = 'all';

    ## -*- texinfo -*-
    ## @deftp {anova} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## Character vector used as the response name in formula display.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName    = 'Y';

    ## -*- texinfo -*-
    ## @deftp {anova} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## Scalar number of response observations used by the model.  This
    ## property is read-only.
    ##
    ## @end deftp
    NumObservations = 0;

    ## Results (populated after fit; empty until then)

    ## -*- texinfo -*-
    ## @deftp {anova} {property} Coefficients
    ##
    ## Model coefficient estimates
    ##
    ## Numeric coefficient table returned by the fitted @code{anovan} backend,
    ## or derived from a @code{LinearModel} object.  This property is read-only.
    ##
    ## @end deftp
    Coefficients   = [];

    ## -*- texinfo -*-
    ## @deftp {anova} {property} Residuals
    ##
    ## Model residuals
    ##
    ## Table with variables @code{Raw} (observed minus fitted values) and
    ## @code{Pearson} (raw residuals scaled by the root mean squared error)
    ## when the selected backend exposes residuals, otherwise empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    Residuals      = [];

    ## -*- texinfo -*-
    ## @deftp {anova} {property} Metrics
    ##
    ## Model fit metrics
    ##
    ## Table with variables @code{MSE}, @code{RMSE}, @code{SSE}, @code{SSR},
    ## @code{SST}, @code{RSquared}, and @code{AdjustedRSquared} summarising the
    ## fitted model.  This property is read-only.
    ##
    ## @end deftp
    Metrics        = [];
  endproperties

  properties (GetAccess = public, SetAccess = private, Hidden)
    ## Backend aliases and Octave-specific extensions retained for internal
    ## use and power users; kept off the documented MATLAB property surface.
    GROUP                                   ## raw factor data (fit workhorse)
    ModelType    = 'linear';
    ModelSpecification = 'linear';
    SSType       = 3;
    VarNames     = {};
    Continuous   = [];
    Random       = [];
    NumFactors   = 0;
    AnovaTable   = {};
    FittedValues = [];
    DFE          = [];
    MSE          = [];
    DesignMatrix = [];
    Stats        = struct ();
  endproperties

  properties (Access = private)
    Alpha       = 0.05;
    Contrasts   = {};
    ## Accepted and validated for compatibility with the anova1 / anova2 /
    ## anovan calling convention, but never forwarded to a backend: results are
    ## presented through summary, disp, and plotDiagnostics instead.
    Display     = 'off';
    Weights     = [];
  endproperties

  properties (Access = private)
    fitted_       = false;
    dirty_        = true;
    nFactors_     = 0;
    backend_      = '';                     ## 'anova1' | 'anova2' | 'anovan'
    reps_         = [];                     ## anova2 replicate count
    sourceModel_  = [];                     ## LinearModel object, when supplied
    rawResiduals_ = [];                     ## raw residual vector for internals
  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {@var{obj} =} anova (@var{Y})
    ## @deftypefnx {anova} {@var{obj} =} anova (@var{factors}, @var{Y})
    ## @deftypefnx {anova} {@var{obj} =} anova (@dots{}, @var{name}, @var{value})
    ## @deftypefnx {anova} {@var{obj} =} anova (@var{mdl})
    ##
    ## Create an object-oriented analysis of variance model.
    ##
    ## @var{Y} is a non-empty numeric response vector or matrix.
    ## @var{factors} contains grouping variables for vector responses and may
    ## be a grouping vector, a matrix of grouping variables, or a cell array of
    ## grouping vectors.  If @var{factors} is omitted and @var{Y} is a matrix,
    ## columns of @var{Y} are treated as groups following @code{anova1} matrix
    ## syntax.
    ##
    ## Supported name-value arguments include @qcode{'ModelSpecification'},
    ## @qcode{'SumOfSquaresType'}, @qcode{'FactorNames'},
    ## @qcode{'CategoricalFactors'}, @qcode{'RandomFactors'},
    ## @qcode{'ResponseName'}, @qcode{'Alpha'}, and @qcode{'Display'}.
    ## Passing @qcode{'Reps'} selects the balanced two-way @code{anova2}
    ## backend when @var{Y} is a non-vector matrix.
    ##
    ## @var{mdl} may be a @code{LinearModel} object, in which case the ANOVA
    ## object is populated from the fitted linear model's public properties.
    ##
    ## @end deftypefn
    function obj = anova (factors, Y, varargin)

      if (nargin < 1)
        error ("anova: too few input arguments.");
      endif
      if (nargin < 2)
        Y = [];
      endif

      if (isa (factors, 'LinearModel'))
        lm_args = varargin;
        if (! isempty (Y))
          lm_args = [{Y}, lm_args];
        endif
        obj = obj.initLinearModel_ (factors, lm_args{:});
        return;
      endif

      if (nargin < 2)
        Y = factors;
        factors = [];
      elseif (isempty (Y) && isnumeric (factors))
        Y = factors;
        factors = [];
      elseif (obj.isName_ (Y))
        varargin = [{Y}, varargin];
        Y = factors;
        factors = [];
      endif

      if (! isnumeric (Y) || isempty (Y))
        error ("anova: Y must be a non-empty numeric array.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error ("anova: name-value pairs must come in pairs.");
      endif

      obj.Y     = Y;
      obj.GROUP = factors;

      ## Parse name-value pairs (mirrors anovan.m's loop style)
      for idx = 1:2:numel (varargin)
        name  = varargin{idx};
        value = varargin{idx + 1};
        if (! ischar (name))
          error ("anova: parameter name must be a character vector.");
        endif
        switch (lower (name))
          case {'model', 'modelspecification'}
            obj = obj.setModelSpecification_ (value);
          case {'sstype', 'sumofsquarestype'}
            obj = obj.setSumOfSquaresType_ (value);
          case {'varnames', 'factornames'}
            obj = obj.setFactorNames_ (value);
          case 'contrasts'
            obj.Contrasts = value;
          case 'alpha'
            obj.Alpha = value;
          case {'categoricalfactors'}
            obj.CategoricalFactors = value;
          case {'continuous'}
            obj.Continuous = value;
          case {'random', 'randomfactors'}
            obj.RandomFactors = value;
            obj.Random = value;
          case 'weights'
            obj.Weights = value;
          case {'display', 'displayopt'}
            obj.Display = value;
          case 'responsename'
            obj.ResponseName = value;
          case 'reps'
            obj.reps_ = value;
          otherwise
            error ("anova: parameter '%s' is not supported.", name);
        endswitch
      endfor

      obj.nFactors_ = obj.countFactors_ ();
      obj.NumFactors = obj.nFactors_;
      obj.NumObservations = numel (obj.Y);
      obj = obj.syncCategoricalFactors_ ();
      obj = obj.setFormula_ ();
      obj.validateSpec_ ();
      obj.validateData_ ();
      obj.Factors = obj.buildFactorsTable_ ();
      obj = obj.selectBackend_ ();
      obj.dirty_ = true;
      ## Fit eagerly at construction so result properties are populated on the
      ## returned value object (value-class semantics; matches MATLAB's anova).
      obj = obj.ensureFit_ ();
    endfunction

  endmethods

  methods (Hidden)

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

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {anova} {} disp (@var{obj})
    ##
    ## Display a compact description of an @code{anova} object.
    ##
    ## The display includes the selected backend, fit state, number of factors,
    ## sum-of-squares type, and significance level.
    ##
    ## @end deftypefn
    function disp (obj)
      obj.ensureFit_ ();
      constrained = 'constrained';
      fprintf ("\n  %d-way anova, %s (Type %s) sums of squares.\n\n", ...
               obj.NumFactors, constrained, obj.sstypeLabel_ ());
      if (! isempty (obj.Formula))
        fprintf ("  %s\n\n", obj.Formula);
      endif
      obj.printAtab_ (obj.AnovaTable);
      fprintf ("\n  Properties, Methods\n\n");
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {@var{s} =} stats (@var{obj})
    ## @deftypefnx {anova} {@var{s} =} stats (@var{obj}, @var{type})
    ## @deftypefnx {anova} {@var{s} =} stats (@var{obj}, @qcode{"Component"}, @var{sstype})
    ## @deftypefnx {anova} {[@var{s}, @var{ems}] =} stats (@dots{})
    ##
    ## Return component or summary ANOVA statistics as a table.
    ##
    ## With no @var{type}, or with @qcode{"component"}, return statistics for
    ## each model term, error, and total.  With @qcode{"summary"}, group terms
    ## into linear, nonlinear, and regression rows.  Replicated continuous
    ## designs also report lack-of-fit and pure-error statistics.
    ##
    ## The @qcode{"Component"} form computes the component table using
    ## @var{sstype}, which must be @qcode{"one"}, @qcode{"two"},
    ## @qcode{"three"}, or @qcode{"hierarchical"}.  This request does not
    ## change the object's read-only @code{SumOfSquaresType} property.
    ## The second output @var{ems} contains expected mean-square information
    ## for each model term and the error term.  Its variables are
    ## @code{Type}, @code{ExpectedMeanSquares},
    ## @code{MeanSquaresDenominator}, @code{DFDenominator}, and
    ## @code{FDenominator}.
    ##
    ## @end deftypefn
    function [s, ems] = stats (obj, varargin)
      type = "component";
      sstype = [];
      if (numel (varargin) == 1)
        type = varargin{1};
      elseif (numel (varargin) == 2 && obj.isName_ (varargin{1}) ...
              && strcmpi (varargin{1}, "component"))
        sstype = obj.parseSSType_ (varargin{2});
      elseif (! isempty (varargin))
        error ("anova.stats: invalid input arguments.");
      endif
      if (! obj.isName_ (type) ...
          || ! any (strcmpi (type, {"component", "summary"})))
        error ("anova.stats: type must be 'component' or 'summary'.");
      endif
      if (strcmpi (type, "summary"))
        atab = obj.summaryStats_ ();
      elseif (! isempty (sstype))
        atab = obj.componentStats_ (sstype);
      else
        obj.ensureFit_ ();
        atab = obj.AnovaTable;
      endif
      s = obj.statsTable_ (atab);
      if (nargout > 1)
        if (isempty (sstype))
          sstype = obj.SSType;
        endif
        ems = obj.expectedMeanSquares_ (sstype);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {@var{means} =} groupmeans (@var{obj})
    ## @deftypefnx {anova} {@var{means} =} groupmeans (@var{obj}, @var{factors})
    ##
    ## Return mean response estimates by factor level.
    ##
    ## The returned value is a @code{table} with one row per factor-level
    ## combination and columns for the level, mean, standard error, and
    ## confidence bounds.
    ##
    ## @end deftypefn
    function means = groupmeans (obj, factors, varargin)
      if (nargin < 2 || isempty (factors))
        factors = 1:obj.NumFactors;
      endif
      alpha = obj.parseAlpha_ (varargin{:});
      obj.ensureFit_ ();
      [G, names] = obj.selectedFactorMatrix_ (factors);
      y = obj.Y(:);
      if (isempty (G))
        error ("anova.groupmeans: factors are required for group means.");
      endif
      [levels, ~, gid] = unique (G, 'rows');
      n = accumarray (gid, 1);
      mu = accumarray (gid, y, [], @mean);
      se = sqrt (obj.MSE ./ n);
      crit = tinv (1 - alpha / 2, max (obj.DFE, 1));
      lo = mu - crit * se;
      hi = mu + crit * se;
      values = num2cell (levels, 1);
      values = [values, {mu, se, lo, hi}];
      vnames = [names, {'Mean', 'SE', 'MeanLower', 'MeanUpper'}];
      means = table (values{:}, 'VariableNames', vnames);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {} boxchart (@var{obj})
    ## @deftypefnx {anova} {@var{h} =} boxchart (@var{obj}, @dots{})
    ##
    ## Plot a one-way box chart of response values by factor.
    ##
    ## This method uses @code{boxplot} as the graphics backend in Octave.
    ##
    ## @end deftypefn
    function h = boxchart (obj, varargin)
      obj.ensureFit_ ();
      [G, ~] = obj.selectedFactorMatrix_ (1);
      if (isempty (G))
        error ("anova.boxchart: a one-way factor is required.");
      endif
      h = boxplot (obj.Y(:), G(:, 1), varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {} plotComparisons (@var{obj})
    ## @deftypefnx {anova} {@var{h} =} plotComparisons (@var{obj}, @dots{})
    ##
    ## Plot multiple-comparison intervals for group means.
    ##
    ## The method delegates interval computation to @code{multcompare}.
    ##
    ## @end deftypefn
    function h = plotComparisons (obj, varargin)
      obj.ensureFit_ ();
      [~, ~, h] = multcompare (obj, 'display', 'on', varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {@var{v} =} varianceComponent (@var{obj})
    ## @deftypefnx {anova} {@var{v} =} varianceComponent (@var{obj}, @dots{})
    ##
    ## Return variance component estimates for the error term.
    ##
    ## Random-factor variance component estimates are not yet implemented.
    ##
    ## @end deftypefn
    function v = varianceComponent (obj, varargin)
      alpha = obj.parseAlpha_ (varargin{:});
      obj.ensureFit_ ();
      if (! isempty (obj.RandomFactors))
        error ("anova.varianceComponent: random factors are not implemented.");
      endif
      lo = obj.DFE * obj.MSE / chi2inv (1 - alpha / 2, obj.DFE);
      hi = obj.DFE * obj.MSE / chi2inv (alpha / 2, obj.DFE);
      v = table ({'Error'}, obj.MSE, lo, hi, 'VariableNames', ...
                 {'Source', 'VarianceComponent', ...
                  'VarianceComponentLower', 'VarianceComponentUpper'});
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

  endmethods

  methods (Hidden)

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
      if (isempty (obj.rawResiduals_) || isempty (obj.FittedValues) ...
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
      h = obj.plotDiagnostics_ (obj.rawResiduals_, obj.FittedValues, ...
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

    function tf = isName_ (obj, value)
      tf = ischar (value) || (isstring (value) && isscalar (value));
    endfunction

    function obj = setModelSpecification_ (obj, value)
      if (isstring (value) && isscalar (value))
        value = char (value);
      endif
      obj.ModelSpecification = value;
      if (ischar (value) && strcmpi (value, 'interactions'))
        obj.ModelType = 'interaction';
      else
        obj.ModelType = value;
      endif
    endfunction

    function obj = setSumOfSquaresType_ (obj, value)
      if (isnumeric (value) && isscalar (value))
        names = {'one', 'two', 'three'};
        if (! any (value == [1, 2, 3]))
          error (strcat ("anova: SumOfSquaresType must be 'one',", ...
                         " 'two', 'three', or 'hierarchical'."));
        endif
        obj.SumOfSquaresType = names{value};
        obj.SSType = value;
        return;
      endif
      if (isstring (value) && isscalar (value))
        value = char (value);
      endif
      if (! ischar (value))
        error ("anova: SumOfSquaresType must be a character vector.");
      endif
      switch (lower (value))
        case {'one', 'typei', 'i'}
          obj.SumOfSquaresType = 'one';
          obj.SSType = 1;
        case {'two', 'typeii', 'ii', 'hierarchical'}
          obj.SumOfSquaresType = lower (value);
          obj.SSType = 2;
        case {'three', 'typeiii', 'iii'}
          obj.SumOfSquaresType = 'three';
          obj.SSType = 3;
        otherwise
          error (strcat ("anova: SumOfSquaresType must be 'one',", ...
                         " 'two', 'three', or 'hierarchical'."));
      endswitch
    endfunction

    function value = parseSSType_ (obj, value)
      if (isstring (value) && isscalar (value))
        value = char (value);
      endif
      if (! ischar (value))
        error ("anova.stats: component type must be a character vector.");
      endif
      switch (lower (value))
        case "one"
          value = 1;
        case "two"
          value = 2;
        case "three"
          value = 3;
        case "hierarchical"
          value = "h";
        otherwise
          error (strcat ("anova.stats: component type must be 'one',", ...
                         " 'two', 'three', or 'hierarchical'."));
      endswitch
    endfunction

    function obj = setFactorNames_ (obj, value)
      if (isstring (value))
        value = cellstr (value(:))';
      elseif (ischar (value))
        value = {value};
      endif
      obj.FactorNames = value;
      obj.VarNames = value;
    endfunction

    function obj = syncCategoricalFactors_ (obj)
      if (obj.isName_ (obj.CategoricalFactors) ...
          && strcmpi (char (obj.CategoricalFactors), 'all'))
        obj.Continuous = [];
        return;
      endif
      if (! isnumeric (obj.CategoricalFactors))
        error (strcat ("anova: CategoricalFactors must be 'all' or a", ...
                       " numeric index vector."));
      endif
      cats = obj.CategoricalFactors(:)';
      if (! isempty (cats) ...
          && (! isvector (cats) || any (cats != fix (cats)) ...
              || any (cats < 1) || any (cats > obj.nFactors_)))
        error ("anova: CategoricalFactors must contain valid factor indices.");
      endif
      obj.Continuous = setdiff (1:obj.nFactors_, cats);
    endfunction

    function obj = setFormula_ (obj)
      if (isempty (obj.FactorNames))
        obj.FactorNames = arrayfun (@(k) sprintf ("Factor%d", k), ...
                                    1:obj.nFactors_, ...
                                    'UniformOutput', false);
        obj.VarNames = obj.FactorNames;
      elseif (numel (obj.FactorNames) != obj.nFactors_)
        error ("anova: FactorNames must contain one name per factor.");
      endif
      if (obj.nFactors_ == 0)
        obj.Formula = sprintf ("%s ~ 1", obj.ResponseName);
      else
        terms = obj.FactorNames(1:obj.nFactors_);
        if (ischar (obj.ModelType) && obj.nFactors_ > 1)
          if (strcmpi (obj.ModelType, 'interaction'))
            for a = 1:obj.nFactors_ - 1
              for b = a + 1:obj.nFactors_
                terms{end + 1} = sprintf ("%s:%s", ...
                                           obj.FactorNames{a}, ...
                                           obj.FactorNames{b});
              endfor
            endfor
          elseif (strcmpi (obj.ModelType, 'full'))
            for mask = 1:(2 ^ obj.nFactors_ - 1)
              idx = find (bitget (mask, 1:obj.nFactors_));
              if (numel (idx) > 1)
                terms{end + 1} = strjoin (obj.FactorNames(idx), ':');
              endif
            endfor
          endif
        endif
        rhs = strjoin (terms, ' + ');
        obj.Formula = sprintf ("%s ~ 1 + %s", obj.ResponseName, rhs);
      endif
    endfunction

    function alpha = parseAlpha_ (obj, varargin)
      alpha = obj.Alpha;
      if (isempty (varargin))
        return;
      endif
      if (mod (numel (varargin), 2) != 0)
        error ("anova: name-value pairs must come in pairs.");
      endif
      for k = 1:2:numel (varargin)
        if (! obj.isName_ (varargin{k}))
          error ("anova: parameter name must be a character vector.");
        endif
        switch (lower (char (varargin{k})))
          case 'alpha'
            alpha = varargin{k + 1};
          otherwise
            error ("anova: parameter '%s' is not supported.", varargin{k});
        endswitch
      endfor
      if (! (isnumeric (alpha) && isscalar (alpha) ...
             && alpha > 0 && alpha < 1))
        error ("anova: Alpha must be a numeric scalar in (0, 1).");
      endif
    endfunction

    function [G, names] = selectedFactorMatrix_ (obj, factors)
      if (ischar (factors))
        idx = find (strcmp (obj.FactorNames, factors));
      elseif (iscellstr (factors))
        idx = cellfun (@(s) find (strcmp (obj.FactorNames, s), 1), factors);
      else
        idx = factors;
      endif
      if (isempty (idx))
        G = [];
        names = {};
        return;
      endif
      if (any (idx < 1) || any (idx > obj.NumFactors))
        error ("anova: factor index exceeds the number of factors.");
      endif
      if (isempty (obj.GROUP) && ! isvector (obj.Y))
        [n, m] = size (obj.Y);
        group_arg = reshape (repmat ((1:m), n, 1), [], 1);
        G = group_arg(:, idx);
      elseif (iscell (obj.GROUP) ...
              && all (cellfun (@(c) isvector (c) || ischar (c), ...
                               obj.GROUP(:))) ...
              && size (obj.GROUP, 1) == 1)
        G = cell2mat (cellfun (@(c) c(:), obj.GROUP(idx), ...
                               'UniformOutput', false));
      else
        G = obj.GROUP(:, idx);
      endif
      names = obj.FactorNames(idx);
    endfunction

    function name = factorName_ (obj, idx)
      if (idx >= 1 && idx <= numel (obj.FactorNames))
        name = obj.FactorNames{idx};
      else
        name = sprintf ("Factor%d", idx);
      endif
    endfunction

    ## Build the public Factors table (one named column per factor) from the
    ## raw GROUP data, leaving the internal GROUP alias untouched.
    function tbl = buildFactorsTable_ (obj)
      if (istable (obj.GROUP))
        tbl = obj.GROUP;
        return;
      endif
      if (isempty (obj.GROUP))
        if (isvector (obj.Y))
          tbl = table ();                   ## intercept-only: no factors
          return;
        endif
        [nr, nc] = size (obj.Y);
        if (isempty (obj.reps_))
          ## one-way column form: single synthetic factor = column index
          cols = {reshape(repmat ((1:nc), nr, 1), [], 1)};
        else
          ## balanced two-way (anova2) form: row-block and column factors
          rowfac = repmat (ceil ((1:nr)' / obj.reps_), nc, 1);
          colfac = reshape (repmat ((1:nc), nr, 1), [], 1);
          cols = {rowfac, colfac};
        endif
      elseif (iscell (obj.GROUP) ...
              && all (cellfun (@(c) isvector (c) || ischar (c), ...
                               obj.GROUP(:))) ...
              && size (obj.GROUP, 1) == 1)
        cols = cellfun (@(c) c(:), obj.GROUP, 'UniformOutput', false);
      elseif (isvector (obj.GROUP))
        cols = {obj.GROUP(:)};
      else
        cols = num2cell (obj.GROUP, 1);
      endif
      tbl = table (cols{:}, 'VariableNames', obj.FactorNames);
    endfunction

    ## Build the public Residuals table (Raw and Pearson) from a raw residual
    ## vector.  Pearson residuals scale the raw residuals by the RMSE, matching
    ## the LinearModel definition.
    function tbl = residualsTable_ (obj, raw)
      raw = raw(:);
      pearson = raw ./ sqrt (max (obj.MSE, eps));
      tbl = table (raw, pearson, 'VariableNames', {'Raw', 'Pearson'});
    endfunction

    ## Build the public Metrics table from the fitted ANOVA table and the
    ## error variance, deriving SSE/SSR/SST and the R-squared measures.
    function tbl = metricsTable_ (obj)
      mse = obj.MSE;
      if (isempty (mse))
        mse = NaN;
      endif
      dfe = obj.DFE;
      if (isempty (dfe))
        dfe = NaN;
      endif
      sse = NaN;
      sst = NaN;
      atab = obj.AnovaTable;
      if (! isempty (atab))
        source_col = obj.findAtabColumn_ (atab, {'Source'});
        ss_col = obj.findAtabColumn_ (atab, {'SS', 'Sum Sq.', 'Sum Sq'});
        for r = 2:rows (atab)
          name = atab{r, source_col};
          if (! ischar (name))
            continue;
          endif
          if (strcmpi (name, 'Error'))
            sse = atab{r, ss_col};
          elseif (strcmpi (name, 'Total'))
            sst = atab{r, ss_col};
          endif
        endfor
      endif
      ssr = sst - sse;
      rsq = ssr / max (sst, eps);
      adj = 1 - (1 - rsq) * (obj.NumObservations - 1) / max (dfe, 1);
      tbl = table (mse, sqrt (max (mse, 0)), sse, ssr, sst, rsq, adj, ...
                   'VariableNames', {'MSE', 'RMSE', 'SSE', 'SSR', 'SST', ...
                                     'RSquared', 'AdjustedRSquared'});
    endfunction

    function out = statsTable_ (obj, atab)
      if (isempty (atab))
        out = table ([], [], [], [], [], "VariableNames", ...
                     {"SumOfSquares", "DF", "MeanSquares", "F", "pValue"});
        return;
      endif
      source_col = obj.findAtabColumn_ (atab, {"Source"});
      columns = {{"SumOfSquares", "Sum Sq.", "Sum Sq", "SS"}, ...
                 {"DF", "d.f.", "df"}, ...
                 {"MeanSquares", "Mean Sq.", "Mean Sq", "MS"}, ...
                 {"F"}, {"pValue", "Prob>F"}};
      names = {"SumOfSquares", "DF", "MeanSquares", "F", "pValue"};
      n = rows (atab) - 1;
      values = cell (1, numel (columns));
      for j = 1:numel (columns)
        col = obj.findAtabColumn_ (atab, columns{j});
        values{j} = NaN (n, 1);
        for i = 1:n
          value = atab{i + 1, col};
          if (isnumeric (value) && isscalar (value) && ! isempty (value))
            values{j}(i) = value;
          endif
        endfor
      endfor
      row_names = obj.publicSourceNames_ (atab(2:end, source_col));
      out = table (values{:}, "VariableNames", names, "RowNames", row_names);
    endfunction

    function ems = expectedMeanSquares_ (obj, sstype)
      obj.ensureFit_ ();
      if (strcmp (obj.backend_, "linearmodel"))
        error (strcat ("anova.stats: expected mean squares are unavailable", ...
                       " for LinearModel input."));
      endif

      atab = obj.componentStats_ (sstype);
      sources = obj.publicSourceNames_ (atab(2:end-1, 1));
      nterms = numel (sources) - 1;
      term_names = sources(1:nterms);
      is_random = false (nterms, 1);
      q_coeff = ones (nterms, 1);
      variance_coeff = zeros (nterms, nterms);

      if (strcmp (obj.backend_, "anovan") && isfield (obj.Stats, "terms"))
        terms = obj.Stats.terms;
        if (! isempty (obj.RandomFactors))
          is_random = any (terms(:, obj.RandomFactors), 2);
        endif
        [q_coeff, variance_coeff] = obj.emsCoefficients_ (terms, sstype, ...
                                                          is_random);
      elseif (strcmp (obj.backend_, "anova1"))
        [q_coeff(1), variance_coeff(1, 1)] = obj.oneWayEmsCoefficient_ ();
        is_random(1) = any (obj.RandomFactors == 1);
      endif

      formulas = cell (nterms + 1, 1);
      types = repmat ({"fixed"}, nterms + 1, 1);
      for i = 1:nterms
        pieces = {};
        if (is_random(i))
          types{i} = "random";
        else
          term = obj.emsTerm_ (q_coeff(i), "Q", term_names{i});
          if (! isempty (term))
            pieces{end + 1} = term;
          endif
        endif
        for j = find (is_random(:))'
          if (abs (variance_coeff(i, j)) > 1e-10)
            pieces{end + 1} = obj.emsTerm_ (variance_coeff(i, j), ...
                                            "V", term_names{j});
          endif
        endfor
        pieces{end + 1} = "V(Error)";
        formulas{i} = strjoin (pieces, "+");
      endfor
      types{end} = "random";
      formulas{end} = "V(Error)";

      denominator_ms = [repmat(obj.MSE, nterms, 1); NaN];
      denominator_df = [repmat(obj.DFE, nterms, 1); NaN];
      denominator_formula = [repmat({"MS(Error)"}, nterms, 1); {""}];
      row_names = [term_names; {"Error"}];
      ems = table (string (types), string (formulas), denominator_ms, ...
                   denominator_df, string (denominator_formula), ...
                   "VariableNames", {"Type", "ExpectedMeanSquares", ...
                   "MeanSquaresDenominator", "DFDenominator", ...
                   "FDenominator"}, "RowNames", row_names);
    endfunction

    function names = publicSourceNames_ (obj, names)
      names = regexprep (names, "'$", "");
      names = regexprep (names, "\\*", ":");
      if (strcmp (obj.backend_, "anova1") && ! isempty (names))
        idx = find (strcmpi (names, "Groups"), 1);
        if (! isempty (idx))
          names{idx} = obj.FactorNames{1};
        endif
      elseif (strcmp (obj.backend_, "anova2"))
        idx = find (strcmpi (names, "Columns"), 1);
        if (! isempty (idx))
          names{idx} = obj.FactorNames{1};
        endif
        idx = find (strcmpi (names, "Rows"), 1);
        if (! isempty (idx))
          names{idx} = obj.FactorNames{2};
        endif
        idx = find (strcmpi (names, "Interaction"), 1);
        if (! isempty (idx))
          names{idx} = sprintf ("%s:%s", obj.FactorNames{:});
        endif
      endif
    endfunction

    function [q_coeff, variance_coeff] = emsCoefficients_ (obj, terms, ...
                                                           sstype, is_random)
      X = full (obj.Stats.X);
      dfs = obj.Stats.df(:);
      nterms = rows (terms);
      blocks = cell (nterms, 1);
      first = 2;
      for i = 1:nterms
        blocks{i} = first:first + dfs(i) - 1;
        first += dfs(i);
      endfor
      bases = cell (nterms, 1);
      for i = 1:nterms
        factors = find (terms(i, :));
        if (all (! ismember (factors, obj.Continuous)))
          [~, ~, ids] = unique (obj.Stats.grps(:, factors), "rows");
          bases{i} = sparse (1:numel (ids), ids, 1);
        else
          bases{i} = X(:, blocks{i});
        endif
      endfor

      q_coeff = zeros (nterms, 1);
      variance_coeff = zeros (nterms, nterms);
      for i = 1:nterms
        [included, excluded] = obj.emsModels_ (X, blocks, terms, i, sstype);
        q_coeff(i) = obj.projectionDifference_ (included, excluded, ...
                                                bases{i}) / dfs(i);
        for j = find (is_random(:))'
          variance_coeff(i, j) = obj.projectionDifference_ (...
                                  included, excluded, bases{j}) / dfs(i);
        endfor
      endfor
    endfunction

    function [included, excluded] = emsModels_ (obj, X, blocks, terms, ...
                                                term, sstype)
      nterms = numel (blocks);
      switch (sstype)
        case 1
          included_terms = 1:term;
          excluded_terms = 1:term-1;
        case {2, "h"}
          factors = find (terms(term, :));
          excluded_terms = find (any (! terms(:, factors), 2))';
          included_terms = [term, excluded_terms];
        otherwise
          included_terms = 1:nterms;
          excluded_terms = setdiff (included_terms, term, "stable");
      endswitch
      included = X(:, [1, blocks{included_terms}]);
      excluded = X(:, [1, blocks{excluded_terms}]);
    endfunction

    function value = projectionDifference_ (obj, included, excluded, basis)
      value = obj.projectionEnergy_ (included, basis) ...
              - obj.projectionEnergy_ (excluded, basis);
      value = max (value, 0);
    endfunction

    function value = projectionEnergy_ (obj, design, basis)
      Q = orth (full (design));
      value = sum (sumsq (Q' * basis));
    endfunction

    function [q_coeff, variance_coeff] = oneWayEmsCoefficient_ (obj)
      if (isvector (obj.Y))
        group = obj.GROUP(:);
      else
        [n, groups] = size (obj.Y);
        group = reshape (repmat (1:groups, n, 1), [], 1);
      endif
      [~, ~, ids] = unique (group, "stable");
      Z = sparse (1:numel (ids), ids, 1);
      intercept = ones (rows (Z), 1);
      df = columns (Z) - 1;
      if (df == 0)
        q_coeff = 0;
      else
        q_coeff = obj.projectionDifference_ (Z, intercept, Z) / df;
      endif
      variance_coeff = q_coeff;
    endfunction

    function value = emsTerm_ (obj, coefficient, symbol, name)
      if (! isfinite (coefficient) || abs (coefficient) <= 1e-10)
        value = "";
      elseif (abs (coefficient - 1) <= 1e-10)
        value = sprintf ("%s(%s)", symbol, name);
      else
        value = sprintf ("%.6g*%s(%s)", coefficient, symbol, name);
      endif
    endfunction

    function obj = initLinearModel_ (obj, mdl, varargin)
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
      obj.NumFactors = obj.nFactors_;
      obj.ModelSpecification = mdl.Formula.Terms;
      obj.ModelType = mdl.Formula.Terms;
      obj.FactorNames = mdl.PredictorNames;
      obj.VarNames = mdl.PredictorNames;
      obj.ResponseName = mdl.ResponseName;
      obj.Formula = mdl.Formula.LinearPredictor;
      obj.Y = mdl.Variables{:, mdl.ResponseName};
      obj.GROUP = mdl.Variables(:, mdl.PredictorNames);
      obj.Factors = obj.GROUP;
      obj.ExpandedFactorNames = mdl.CoefficientNames;
      obj.NumObservations = numel (obj.Y);
      obj.Coefficients = obj.linearModelCoefficients_ (mdl);
      obj.AnovaTable = obj.linearModelAtab_ (mdl);
      obj.FittedValues = mdl.Fitted;
      obj.DFE = mdl.DFE;
      obj.MSE = mdl.MSE;
      obj.rawResiduals_ = mdl.Residuals.Raw;
      obj.Residuals = obj.residualsTable_ (obj.rawResiduals_);
      obj.DesignMatrix = [];
      obj.Stats = obj.linearModelStats_ (mdl);
      obj.Metrics = obj.metricsTable_ ();
      obj.fitted_ = true;
      obj.dirty_ = false;
    endfunction

    function validateSpec_ (obj)
      if (! (isnumeric (obj.SSType) && isscalar (obj.SSType) ...
             && any (obj.SSType == [1, 2, 3])))
        error (strcat ("anova: SumOfSquaresType must be 'one',", ...
                       " 'two', 'three', or 'hierarchical'."));
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
          error (strcat ("anova: ModelSpecification must be 'linear',", ...
                         " 'interactions', 'full', or a terms matrix."));
        endif
      elseif (! isnumeric (obj.ModelType))
        error (strcat ("anova: ModelSpecification must be a string", ...
                       " or a numeric terms matrix."));
      endif
      if (! iscellstr (obj.FactorNames))
        error ("anova: FactorNames must be a character vector or cellstr.");
      endif
      if (! obj.isName_ (obj.ResponseName))
        error ("anova: ResponseName must be a character vector.");
      endif
      if (! isempty (obj.Continuous) && ! isnumeric (obj.Continuous))
        error (strcat ("anova: CategoricalFactors must be 'all' or a", ...
                       " numeric index vector."));
      endif
      if (! isempty (obj.Continuous) ...
          && (! isvector (obj.Continuous) ...
              || any (obj.Continuous != fix (obj.Continuous)) ...
              || any (obj.Continuous < 1)))
        error ("anova: CategoricalFactors must contain valid factor indices.");
      endif
      if (! isempty (obj.RandomFactors) && ! isnumeric (obj.RandomFactors))
        error ("anova: RandomFactors must be a numeric index vector.");
      endif
      if (! isempty (obj.RandomFactors) ...
          && (! isvector (obj.RandomFactors) ...
              || any (obj.RandomFactors != fix (obj.RandomFactors)) ...
              || any (obj.RandomFactors < 1)))
        error ("anova: RandomFactors must contain positive integer indices.");
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
        error ("anova: CategoricalFactors must contain valid factor indices.");
      endif
      if (! isempty (obj.RandomFactors) ...
          && any (obj.RandomFactors > obj.nFactors_))
        error ("anova: RandomFactors indices exceed the number of factors.");
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
            error (strcat ("anova: GROUP variables must match the number", ...
                           " of observations."));
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
          if (isempty (obj.reps_))
            nf = 1;
          else
            nf = 2;
          endif
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
    function obj = selectBackend_ (obj)
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
    function obj = ensureFit_ (obj)
      if (! obj.fitted_ || obj.dirty_)
        obj = obj.selectBackend_ ();
        obj = obj.fit_ ();
      endif
    endfunction

    ## Dispatch to the selected backend; populate result properties.
    function obj = fit_ (obj)
      switch (obj.backend_)
        case 'anova1'
          obj = obj.fitAnova1_ ();
        case 'anova2'
          obj = obj.fitAnova2_ ();
        case 'anovan'
          obj = obj.fitAnovan_ ();
        case 'linearmodel'
          ## Already populated from the supplied LinearModel object.
      endswitch
      if (! strcmp (obj.backend_, 'linearmodel'))
        obj.Metrics = obj.metricsTable_ ();
      endif
      obj.fitted_ = true;
      obj.dirty_  = false;
    endfunction

    function obj = fitAnova1_ (obj)
      ## Always run the anova1 backend silently; anova presents results through
      ## its own summary, disp, and plotDiagnostics methods, so the backend must
      ## not print a table or open its boxplot figure during a fit.
      if (isvector (obj.Y))
        [~, atab, stats] = anova1 (obj.Y, obj.GROUP, 'off');
      else
        [~, atab, stats] = anova1 (obj.Y, [], 'off');
      endif
      obj.AnovaTable = atab;
      obj.Stats      = stats;
      obj.DFE        = stats.df;
      obj.MSE        = stats.s ^ 2;     ## anova1 reports sqrt(MSE) as s
      ## Coefficients / Residuals / DesignMatrix / FittedValues are not
      ## exposed by anova1's stats; they remain at their empty defaults.
    endfunction

    function obj = fitAnova2_ (obj)
      modelarg = 'interaction';
      if (ischar (obj.ModelType))
        modelarg = obj.ModelType;
      endif
      ## Always run the anova2 backend silently; see fitAnova1_.
      [~, atab, stats] = anova2 (obj.Y, obj.reps_, 'off', modelarg);
      obj.AnovaTable = atab;
      obj.Stats      = stats;
      obj.DFE        = stats.df;
      obj.MSE        = stats.sigmasq;
      ## anova2 does not return coeffs / resid / X.
    endfunction

    function obj = fitAnovan_ (obj)
      [y_vec, group_arg] = obj.anovanData_ ();
      [~, atab, stats] = anovan (y_vec, group_arg, ...
                                 obj.buildAnovanArgs_(){:});
      obj.AnovaTable = atab;
      obj.Stats      = stats;
      if (isfield (stats, 'coeffs'))
        obj.Coefficients = stats.coeffs;
      endif
      if (isfield (stats, 'coeffnames'))
        obj.ExpandedFactorNames = cellstr (stats.coeffnames);
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
        ## anovan drops observations with a missing (NaN) response or numeric
        ## grouping value, so the design matrix (and fitted values) cover only
        ## the retained rows.  Align the response with those rows before
        ## differencing; expose raw (observed minus fitted) residuals so that
        ## FittedValues + Residuals == Y over the retained observations.
        keep = ! isnan (y_vec);
        if (iscell (group_arg))
          for k = 1:numel (group_arg)
            if (isnumeric (group_arg{k}))
              keep = keep & ! isnan (group_arg{k}(:));
            endif
          endfor
        elseif (isnumeric (group_arg) && ! isempty (group_arg))
          keep = keep & ! any (isnan (group_arg), 2);
        endif
        if (numel (obj.FittedValues) == sum (keep))
          obj.rawResiduals_ = y_vec(keep) - obj.FittedValues;
          obj.Residuals = obj.residualsTable_ (obj.rawResiduals_);
        endif
      endif
    endfunction

    function atab = componentStats_ (obj, sstype)
      if (strcmp (obj.backend_, "linearmodel"))
        error ("anova.stats: component type is unavailable for LinearModel input.");
      endif
      if (strcmp (obj.backend_, "anova2"))
        obj.ensureFit_ ();
        atab = obj.AnovaTable;
        return;
      endif
      [y_vec, group_arg] = obj.anovanData_ ();
      args = obj.buildAnovanArgs_ (sstype);
      [~, atab] = anovan (y_vec, group_arg, args{:});
    endfunction

    function atab = summaryStats_ (obj)
      obj.ensureFit_ ();
      component = obj.componentStats_ (1);
      sources = component(2:end, 1);
      error_idx = find (strcmpi (sources, "Error"), 1);
      total_idx = find (strcmpi (sources, "Total"), 1);
      term_idx = setdiff (1:numel (sources), [error_idx, total_idx], "stable");
      term_ss = cell2mat (component(term_idx + 1, 2));
      term_df = cell2mat (component(term_idx + 1, 3));

      if (strcmp (obj.backend_, "anova2"))
        nonlinear = strcmpi (sources(term_idx), "Interaction");
      elseif (isfield (obj.Stats, "terms") ...
              && rows (obj.Stats.terms) == numel (term_idx))
        nonlinear = (sum (obj.Stats.terms, 2) > 1);
      else
        nonlinear = ! cellfun (@isempty, regexp (sources(term_idx), "[:*]"));
      endif
      linear = ! nonlinear;

      labels = {};
      sum_sq = [];
      df = [];
      if (any (linear))
        labels{end + 1, 1} = "Linear";
        sum_sq(end + 1, 1) = sum (term_ss(linear));
        df(end + 1, 1) = sum (term_df(linear));
      endif
      if (any (nonlinear))
        labels{end + 1, 1} = "NonLinear";
        sum_sq(end + 1, 1) = sum (term_ss(nonlinear));
        df(end + 1, 1) = sum (term_df(nonlinear));
      endif
      labels{end + 1, 1} = "Regression";
      sum_sq(end + 1, 1) = sum (term_ss);
      df(end + 1, 1) = sum (term_df);
      tested_rows = numel (labels);

      error_ss = component{error_idx + 1, 2};
      error_df = component{error_idx + 1, 3};
      labels{end + 1, 1} = "Error";
      sum_sq(end + 1, 1) = error_ss;
      df(end + 1, 1) = error_df;
      error_row = numel (labels);

      lack_row = [];
      if (strcmp (obj.backend_, "anovan") && isempty (obj.Weights) ...
          && isfield (obj.Stats, "grps") && isfield (obj.Stats, "Y"))
        [~, ~, group_id] = unique (obj.Stats.grps, "rows");
        group_mean = accumarray (group_id, obj.Stats.Y, [], @mean);
        pure_ss = sum ((obj.Stats.Y - group_mean(group_id)) .^ 2);
        pure_df = numel (obj.Stats.Y) - max (group_id);
        lack_df = error_df - pure_df;
        if (pure_df > 0 && lack_df > 0)
          labels(end + 1:end + 2, 1) = {"LackOfFit"; "PureError"};
          sum_sq(end + 1:end + 2, 1) = [max(error_ss - pure_ss, 0); pure_ss];
          df(end + 1:end + 2, 1) = [lack_df; pure_df];
          lack_row = numel (labels) - 1;
        endif
      endif

      labels{end + 1, 1} = "Total";
      sum_sq(end + 1, 1) = component{total_idx + 1, 2};
      df(end + 1, 1) = component{total_idx + 1, 3};
      mean_sq = sum_sq ./ df;
      mse = mean_sq(error_row);
      f_stat = mean_sq(1:tested_rows) ./ mse;
      p_value = fcdf (f_stat, df(1:tested_rows), error_df, "upper");

      atab = cell (numel (labels) + 1, 6);
      atab(1, :) = {"Source", "Sum Sq.", "d.f.", "Mean Sq.", "F", "Prob>F"};
      atab(2:end, 1) = labels;
      atab(2:end, 2:4) = num2cell ([sum_sq, df, mean_sq]);
      atab(2:tested_rows + 1, 5:6) = num2cell ([f_stat, p_value]);
      if (! isempty (lack_row))
        lack_f = mean_sq(lack_row) / mean_sq(lack_row + 1);
        lack_p = fcdf (lack_f, df(lack_row), df(lack_row + 1), "upper");
        atab(lack_row + 1, 5:6) = {lack_f, lack_p};
      endif
    endfunction

    function [y_vec, group_arg] = anovanData_ (obj)
      if (isempty (obj.GROUP) && ! isvector (obj.Y))
        [n, m] = size (obj.Y);
        y_vec = obj.Y(:);
        group = reshape (repmat (1:m, n, 1), [], 1);
        group_arg = {group};
      else
        y_vec = obj.Y(:);
        group_arg = obj.GROUP;
        if (isempty (group_arg))
          group_arg = {};
        endif
      endif
    endfunction

    function s = sstypeLabel_ (obj)
      switch (obj.SSType)
        case 1; s = 'I';
        case 2; s = 'II';
        otherwise; s = 'III';
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

    function nv = buildAnovanArgs_ (obj, sstype)
      if (nargin < 2)
        sstype = obj.SSType;
      endif
      ## Run the backend silently; this class owns display and plotting.
      nv = {'display', 'off', 'sstype', sstype, 'alpha', obj.Alpha};
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
      D = (obj.rawResiduals_ .^ 2 ./ max (p * obj.MSE, eps)) ...
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

      fig_name = 'Diagnostic Plots: Model Residuals';
      visible = 'on';
      if (mod (numel (varargin), 2) != 0)
        error ("anova.plotDiagnostics: name-value pairs must come in pairs.");
      endif
      for k = 1:2:numel (varargin)
        switch (lower (varargin{k}))
          case 'figurename'
            fig_name = varargin{k + 1};
          case 'visible'
            visible = varargin{k + 1};
          otherwise
            error ("anova.plotDiagnostics: unknown option '%s'.", varargin{k});
        endswitch
      endfor

      mse = sum (residuals .^ 2) / max (dfe, 1);
      t = residuals ./ sqrt (mse * max (1 - leverage, eps));
      [~, DI] = sort (cooksd, 'descend');
      nk = min (4, n);

      h = figure ('Name', fig_name, 'Visible', visible);

      subplot (2, 2, 1);
      x = ((1:n)' - 0.5) / n;
      [ts, I] = sort (t);
      q = norminv (x);
      plot (q, ts, 'ok', 'markersize', 3);
      box off;
      grid on;
      xlabel ('Theoretical quantiles');
      ylabel ('Studentized residuals');
      title ('Normal Q-Q Plot');
      arrayfun (@(i) text (q(I == DI(i)), t(DI(i)), ...
                           sprintf ("  %u", DI(i))), 1:nk);
      iqr = [0.25; 0.75];
      yl = quantile (t, iqr, 1, 6);
      xl = norminv (iqr);
      slope = diff (yl) / diff (xl);
      int = yl(1) - slope * xl(1);
      ax1_xlim = get (gca, 'XLim');
      hold on;
      plot (ax1_xlim, slope * ax1_xlim + int, 'k-');
      hold off;
      set (gca, 'Xlim', ax1_xlim);

      subplot (2, 2, 2);
      plot (fitted, sqrt (abs (t)), 'ko', 'markersize', 3);
      box off;
      xlabel ('Fitted values');
      ylabel ('sqrt ( | Studentized residuals | )');
      title ('Spread-Location Plot');
      ax2_xlim = get (gca, 'XLim');
      hold on;
      plot (ax2_xlim, ones (1, 2) * sqrt (2), 'k:');
      plot (ax2_xlim, ones (1, 2) * sqrt (3), 'k-.');
      plot (ax2_xlim, ones (1, 2) * sqrt (4), 'k--');
      hold off;
      arrayfun (@(i) text (fitted(DI(i)), sqrt (abs (t(DI(i)))), ...
                           sprintf ("  %u", DI(i))), 1:nk);
      xlim (ax2_xlim);

      subplot (2, 2, 3);
      plot (leverage, t, 'ko', 'markersize', 3);
      box off;
      xlabel ('Leverage');
      ylabel ('Studentized residuals');
      title ('Residual-Leverage Plot');
      ax3_xlim = get (gca, 'XLim');
      ax3_ylim = get (gca, 'YLim');
      hold on;
      plot (ax3_xlim, zeros (1, 2), 'k-');
      hold off;
      arrayfun (@(i) text (leverage(DI(i)), t(DI(i)), ...
                           sprintf ("  %u", DI(i))), 1:nk);
      set (gca, 'ygrid', 'on');
      xlim (ax3_xlim);
      ylim (ax3_ylim);

      subplot (2, 2, 4);
      stem (cooksd, 'ko', 'markersize', 3);
      box off;
      xlabel ('Obs. number');
      ylabel ('Cook''s distance');
      title ('Cook''s Distance Stem Plot');
      xlim ([0, n]);
      ax4_xlim = get (gca, 'XLim');
      ax4_ylim = get (gca, 'YLim');
      hold on;
      plot (ax4_xlim, ones (1, 2) * 4 / max (dfe, eps), 'k:');
      plot (ax4_xlim, ones (1, 2) * 0.5, 'k-.');
      plot (ax4_xlim, ones (1, 2), 'k--');
      hold off;
      arrayfun (@(i) text (DI(i), cooksd(DI(i)), ...
                           sprintf ("  %u", DI(i))), 1:nk);
      xlim (ax4_xlim);
      ylim (ax4_ylim);

      set (findall (gcf, '-property', 'FontSize'), 'FontSize', 7);
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
      has_intercept = isa (mdl.Formula, 'LinearFormula') ...
                      && mdl.Formula.HasIntercept;
      intercept_df = 0;
      if (has_intercept)
        intercept_df = 1;
      endif
      df_model = mdl.NumEstimatedCoefficients - intercept_df;
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
              'Total', mdl.SST, mdl.NumObservations - intercept_df, ...
              '', '', ''};
    endfunction

    function stats = linearModelStats_ (obj, mdl)
      stats = struct ('source', 'linearmodel', ...
                      'resid', mdl.Residuals.Raw, ...
                      'coeffs', obj.Coefficients, ...
                      'dfe', mdl.DFE, ...
                      'mse', mdl.MSE, ...
                      'vcov', mdl.CoefficientCovariance, ...
                      'CooksD', mdl.Diagnostics.CooksDistance, ...
                      'alpha', obj.Alpha, ...
                      'varnames', {mdl.VariableNames}, ...
                      'coeffnames', {mdl.CoefficientNames});
    endfunction

    function es = linearModelEffectSizes_ (obj)
      mdl = obj.sourceModel_;
      eta = mdl.SSR / max (mdl.SST, eps);
      partial_eta = mdl.SSR / max (mdl.SSR + mdl.SSE, eps);
      has_intercept = isa (mdl.Formula, 'LinearFormula') ...
                      && mdl.Formula.HasIntercept;
      intercept_df = 0;
      if (has_intercept)
        intercept_df = 1;
      endif
      df_model = max (mdl.NumEstimatedCoefficients - intercept_df, 0);
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
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! summary (a);

%!demo
%! ## Post-hoc multiple comparisons
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! C = multcompare (a, 'display', 'off')

%!demo
%! ## Diagnostic plots for an anovan-backed fit
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! plotDiagnostics (a);

## --- BISTs ---------------------------------------------------------------

## Basic construction: vector Y + single grouping vector
%!test
%! y = [1; 2; 3; 4; 5; 6];
%! g = [1; 1; 2; 2; 3; 3];
%! a = anova (g, y);
%! assert_equal (class (a), 'anova');
%! assert_equal (a.Y, y);
%! assert_equal (a.GROUP, g);

## Basic construction: matrix Y, no GROUP
%!test
%! y = magic (4);
%! a = anova (y);
%! assert_equal (class (a), 'anova');
%! assert_equal (a.Y, y);
%! assert_equal (a.GROUP, []);

## Basic construction: cell of two grouping vectors
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = repmat ([1;1;2;2], 3, 1);
%! a = anova ({g1, g2}, y);
%! assert_equal (class (a), 'anova');
%! assert_equal (a.NumFactors, 2);

## Property defaults
%!test
%! a = anova ([1;1;2;2], [1;2;3;4]);
%! assert_equal (a.ModelSpecification, 'linear');
%! assert_equal (a.SumOfSquaresType, 'three');
%! assert_equal (a.ResponseName, 'Y');
%! assert_equal (a.FactorNames, {'Factor1'});
%! assert_equal (a.RandomFactors, []);
%! assert_equal (a.CategoricalFactors, 'all');
%! assert_equal (a.NumFactors, 1);
%! assert_equal (a.NumObservations, 4);

## Eager fit at construction populates the table statistics; the anova1
## backend does not expose coefficients, residuals, fitted values, or a
## design matrix, so those stay at their empty defaults
%!test
%! a = anova ([1;1;2;2], [1;2;3;4]);
%! assert_equal (isempty (a.AnovaTable), false);
%! assert_equal (isempty (a.DFE), false);
%! assert_equal (isempty (a.MSE), false);
%! assert_equal (a.Coefficients, []);
%! assert_equal (a.Residuals, []);
%! assert_equal (a.FittedValues, []);
%! assert_equal (a.DesignMatrix, []);

## Name-value parsing: MATLAB-compatible names
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = repmat ([1;1;2;2], 3, 1);
%! a = anova ({g1, g2}, y, 'SumOfSquaresType', 'two', ...
%!            'ModelSpecification', 'full', 'FactorNames', {'A', 'B'});
%! assert_equal (a.SumOfSquaresType, 'two');
%! assert_equal (a.ModelSpecification, 'full');
%! assert_equal (a.FactorNames, {'A', 'B'});

## Name-value parsing: case-insensitive names, displayopt alias
%!test
%! a = anova ([1;1;2;2], [1;2;3;4], 'SumOfSquaresType', 'one', ...
%!            'displayopt', 'on');
%! assert_equal (a.SumOfSquaresType, 'one');

## Display 'on' never reaches a backend: no fit prints a table or opens a
## figure.  The three constructions select the anova1, anova2, and anovan
## backends respectively.
%!test
%! y = [1; 2; 3; 4];
%! g = [1; 1; 2; 2];
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! figs = get (0, 'children');
%! cmd = strcat ("a1 = anova (g, y, 'Display', 'on');", ...
%!               "a2 = anova (popcorn, [], 'reps', 3, 'Display', 'on');", ...
%!               "a3 = anova (g, y, 'SumOfSquaresType', 'one',", ...
%!               " 'Display', 'on');");
%! str = evalc (cmd);
%! assert_equal (isempty (str), true);
%! assert_equal (isempty (setdiff (get (0, 'children'), figs)), true);

## Backend selection: one-way default -> anova1
%!test
%! a = anova ([1;1;2;2;3;3], [1;2;3;4;5;6]);
%! assert_equal (isempty (strfind (evalc ('disp (a)'), '1-way anova')), false);

## Backend selection: one-way matrix-Y form -> anova1
%!test
%! a = anova (magic (4));
%! str = evalc ('disp (a)');
%! assert_equal (isempty (strfind (str, '1-way anova')), false);

## Backend selection: matrix Y + explicit 'reps' -> anova2
%!test
%! y = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!      6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (y, [], 'reps', 3);
%! assert_equal (isempty (strfind (evalc ('disp (a)'), '2-way anova')), false);

## Backend selection: two-factor cell groups without reps -> anovan
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = repmat ([1;1;2;2], 3, 1);
%! a = anova ({g1, g2}, y);
%! str = evalc ('disp (a)');
%! assert_equal (isempty (strfind (str, '2-way anova')), false);

## Backend selection: three factors -> anovan
%!test
%! y = (1:24)';
%! g1 = repmat ([1;2], 12, 1);
%! g2 = repmat ([1;1;2;2], 6, 1);
%! g3 = repmat ([1;1;1;1;2;2;2;2], 3, 1);
%! a = anova ({g1, g2, g3}, y);
%! str = evalc ('disp (a)');
%! assert_equal (isempty (strfind (str, '3-way anova')), false);

## Backend selection: SSType != 3 with 1 factor falls through to anovan
%!test
%! a = anova ([1;1;2;2;3;3], [1;2;3;4;5;6], 'SumOfSquaresType', 'two');
%! assert_equal (isempty (strfind (evalc ('disp (a)'), 'Type II')), false);

## Backend selection: continuous predictors force anovan
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = (1:12)';                              ## continuous
%! a = anova ({g1, g2}, y, 'CategoricalFactors', 1);
%! assert_equal (a.CategoricalFactors, 1);

## Backend selection: NaN in a vector response falls through to anovan
%!test
%! y = [1; 2; 3; NaN; 5; 6; 7; 8];
%! g1 = [1; 2; 1; 2; 1; 2; 1; 2];
%! g2 = [1; 1; 2; 2; 1; 1; 2; 2];
%! a = anova ({g1, g2}, y);
%! assert_equal (a.NumFactors, 2);

## Backend selection: weights force anovan
%!test
%! a = anova ([1;1;2;2;3;3], [1;2;3;4;5;6], 'Weights', ones (6, 1));
%! assert_equal (isempty (stats (a)), false);

## --- Week 2: fit delegation smoke tests --------------------------------

## fit_(): one-way fixture populates the unified result surface
%!test
%! y = [1; 2; 3; 4; 5; 6];
%! g = [1; 1; 2; 2; 3; 3];
%! a = anova (g, y);
%! a.fit ();
%! assert_equal (isempty (a.AnovaTable), false);
%! assert_equal (class (a.Stats), 'struct');
%! assert_equal (isempty (a.DFE), false);
%! assert_equal (isempty (a.MSE), false);

## fit_(): two-way balanced fixture (anova2 backend, popcorn data)
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (popcorn, [], 'reps', 3);
%! assert_equal (! isempty (strfind (evalc ('disp (a)'), '2-way anova')), true);
%! a.fit ();
%! assert_equal (! isempty (a.AnovaTable), true);
%! assert_equal (isfield (a.Stats, 'sigmasq'), true);
%! assert_equal (a.MSE, a.Stats.sigmasq);

## fit_(): N-way fixture (anovan backend, three factors)
%!test
%! y = (1:24)';
%! g1 = repmat ([1;2], 12, 1);
%! g2 = repmat ([1;1;2;2], 6, 1);
%! g3 = repmat ([1;1;1;1;2;2;2;2], 3, 1);
%! a = anova ({g1, g2, g3}, y);
%! assert_equal (a.NumFactors, 3);
%! a.fit ();
%! assert_equal (! isempty (a.AnovaTable), true);
%! assert_equal (! isempty (a.Coefficients), true);
%! assert_equal (! isempty (a.Residuals), true);
%! assert_equal (! isempty (a.DesignMatrix), true);
%! assert_equal (rows (a.Residuals), numel (y));

## fit_(): FittedValues = DesignMatrix * Coefficients(:,1) for anovan
%!test
%! y  = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g  = [1;1;1;2;2;2;3;3;3];
%! a  = anova (g, y, 'SumOfSquaresType', 'two');
%! a.fit ();
%! assert_equal (numel (a.FittedValues), numel (y));
%! assert_equal (a.FittedValues + a.Residuals.Raw, y, 1e-9);

## ensureFit_(): fit() is idempotent (second call does nothing)
%!test
%! a = anova ([1;1;1;2;2;2;3;3;3], (1:9)');
%! a.fit ();
%! first_table = a.AnovaTable;
%! a.fit ();
%! assert_equal (a.AnovaTable, first_table);

## buildAnovanArgs_(): SSType and Alpha are forwarded to anovan
%!test
%! y = (1:12)';
%! g = repmat ([1;2;3], 4, 1);
%! a = anova ({g}, y, 'SumOfSquaresType', 'two', 'Alpha', 0.10);
%! a.fit ();
%! assert_equal (a.SumOfSquaresType, 'two');
%! assert_equal (a.Stats.alpha, 0.10, 1e-12);

## --- Numeric references -------------------------------------------------

## One-way ANOVA: reference values match R's aov(y ~ factor(g)).
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! a.fit ();
%! T = a.AnovaTable;
%! assert_equal (T{2, 2}, 126, 1e-12);
%! assert_equal (T{3, 2}, 6, 1e-12);
%! assert_equal (T{4, 2}, 132, 1e-12);
%! assert_equal (T{2, 3}, 2);
%! assert_equal (T{3, 3}, 6);
%! assert_equal (T{2, 4}, 63, 1e-12);
%! assert_equal (T{3, 4}, 1, 1e-12);
%! assert_equal (T{2, 6}, 63, 1e-12);
%! assert_equal (T{2, 7}, 9.3914e-05, 1e-9);
%! assert_equal (a.MSE, 1, 1e-12);
%! assert_equal (a.DFE, 6);

## Fitted values and residuals: reference values from one-way cell means.
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! assert_equal (predict (a), [2; 2; 2; 5; 5; 5; 11; 11; 11], 1e-12);
%! assert_equal (a.Residuals.Raw, [-1; 0; 1; -1; 0; 1; -1; 0; 1], 1e-12);

## Effect sizes: reference eta2 = SS/SST, omega2 = (SS-df*MSE)/(SST+MSE).
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! es = getEffectSizes (anova (g, y, 'SumOfSquaresType', 'two'));
%! assert_equal (es.Source, {'Factor1'});
%! assert_equal (es.EtaSquared, 126 / 132, 1e-12);
%! assert_equal (es.PartialEtaSquared, 126 / (126 + 6), 1e-12);
%! assert_equal (es.OmegaSquared, 124 / 133, 1e-12);

## Two-way balanced ANOVA: popcorn values match the anova2 doc example.
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (popcorn, [], 'reps', 3);
%! a.fit ();
%! T = a.AnovaTable;
%! assert_equal (T{2, 2}, 15.75, 1e-12);
%! assert_equal (T{3, 2}, 4.5, 1e-12);
%! assert_equal (T{4, 2}, 1.75, 1e-12);
%! assert_equal (T{5, 2}, 22, 1e-12);
%! assert_equal (T{2, 5}, 63, 1e-12);
%! assert_equal (T{3, 5}, 36, 1e-12);
%! assert_equal (a.MSE, 0.125, 1e-12);

## Post-hoc comparisons: reference values checked against R's TukeyHSD output.
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! C = multcompare (a, 'display', 'off');
%! assert_equal (C(:, 1:2), [1 2; 1 3; 2 3]);
%! assert_equal (C(:, 4), [-3; -9; -6], 1e-12);
%! assert_equal (C(1, 6), 0.010402, 1e-6);
%! assert_equal (C(2, 6), 9.9474e-05, 1e-9);
%! assert_equal (C(3, 6), 0.00064995, 1e-8);

## --- Week 3: summary / disp ---------------------------------------------

## summary(): runs ensureFit_ and prints a table
%!test
%! a = anova ([1;1;1;2;2;2;3;3;3], (1:9)', 'SumOfSquaresType', 'two');
%! str = evalc ('summary (a)');
%! assert_equal (! isempty (strfind (str, 'ANOVA TABLE')), true);
%! assert_equal (! isempty (strfind (str, 'backend = anovan')), true);

## summary(): includes MSE / DFE / Alpha line
%!test
%! a = anova ([1;1;1;2;2;2;3;3;3], (1:9)', 'SumOfSquaresType', 'two', ...
%!            'Alpha', 0.10);
%! str = evalc ('summary (a)');
%! assert_equal (! isempty (strfind (str, 'Alpha: 0.1')), true);

## disp(): one-line overview of key fields
%!test
%! a = anova ([1;1;1;2;2;2;3;3;3], (1:9)', 'SumOfSquaresType', 'two');
%! str = evalc ('disp (a)');
%! assert_equal (! isempty (strfind (str, '1-way anova')), true);
%! assert_equal (! isempty (strfind (str, 'Type II')), true);
%! assert_equal (! isempty (strfind (str, 'Properties, Methods')), true);

## summary(): SSType label appears in the header (Type I / II / III)
%!test
%! a1 = anova ([1;1;1;2;2;2;3;3;3], (1:9)', 'SumOfSquaresType', 'one');
%! a2 = anova ([1;1;1;2;2;2;3;3;3], (1:9)', 'SumOfSquaresType', 'two');
%! assert_equal (! isempty (strfind (evalc ('summary (a1)'), ...
%!                                   'Type I sums')), true);
%! assert_equal (! isempty (strfind (evalc ('summary (a2)'), ...
%!                                   'Type II sums')), true);

## --- Week 4: multcompare pass-through ----------------------------------

## multcompare(): anovan backend returns a pairwise comparison matrix
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! C = multcompare (a, 'display', 'off');
%! assert_equal (! isempty (C), true);
%! assert_equal (size (C, 1), 3);        ## 3 pairwise comparisons for 3 groups
%! assert_equal (size (C, 2) >= 6, true);  ## at least i,j,diff,lo,hi,p

## multcompare(): two-way balanced via anova2 backend
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (popcorn, [], 'reps', 3);
%! C = multcompare (a, 'display', 'off', 'estimate', 'column');
%! assert_equal (! isempty (C), true);
%! assert_equal (size (C, 2) >= 6, true);

## multcompare(): runs after a fresh construction (triggers ensureFit_)
%!test
%! a = anova ([1;1;2;2;3;3], (1:6)', 'SumOfSquaresType', 'two');
%! C = multcompare (a, 'display', 'off');
%! assert_equal (! isempty (C), true);

## MATLAB-compatible public methods and properties
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'FactorNames', {'Brand'}, ...
%!            'ResponseName', 'Yield', 'SumOfSquaresType', 'two');
%! assert_equal (a.Formula, 'Yield ~ 1 + Brand');
%! assert_equal (a.Factors.Brand, g);
%! assert_equal (a.Y, y);
%! assert_equal (a.FactorNames, {'Brand'});
%! assert_equal (a.SumOfSquaresType, 'two');
%! assert_equal (! any (strcmp (methods ('anova'), 'predict'), 'all'), true);
%! T = stats (a);
%! M = groupmeans (a);
%! V = varianceComponent (a);
%! assert_equal (T.SumOfSquares(1), 126, 1e-12);
%! assert_equal (M.Mean, [2; 5; 11], 1e-12);
%! assert_equal (V.VarianceComponent, 1, 1e-12);

## Public property surface matches MATLAB's anova object (13 read-only names)
%!test
%! a = anova ([1;1;1;2;2;2;3;3;3], (1:9)');
%! assert_equal (sort (properties (a)), sort ({'Y'; 'Factors'; 'Formula'; ...
%!   'FactorNames'; 'ExpandedFactorNames'; 'SumOfSquaresType'; ...
%!   'RandomFactors'; 'CategoricalFactors'; 'ResponseName'; ...
%!   'NumObservations'; 'Coefficients'; 'Residuals'; 'Metrics'}));

## Factors is a table with one named column per factor
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = repmat ([1;1;2;2], 3, 1);
%! a = anova ({g1, g2}, y, 'FactorNames', {'A', 'B'});
%! assert_equal (istable (a.Factors), true);
%! assert_equal (a.Factors.A, g1);
%! assert_equal (a.Factors.B, g2);

## Residuals is a Raw/Pearson table (Pearson = Raw ./ sqrt (MSE))
%!test
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! a.fit ();
%! assert_equal (istable (a.Residuals), true);
%! assert_equal (a.Residuals.Properties.VariableNames, {'Raw', 'Pearson'});
%! assert_equal (a.Residuals.Pearson, a.Residuals.Raw ./ sqrt (a.MSE), 1e-12);

## Metrics table exposes the fit summary with a correct R-squared
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! a.fit ();
%! M = a.Metrics;
%! assert_equal (M.Properties.VariableNames, {'MSE', 'RMSE', 'SSE', 'SSR', ...
%!   'SST', 'RSquared', 'AdjustedRSquared'});
%! assert_equal (M.SSE, 6, 1e-12);
%! assert_equal (M.SST, 132, 1e-12);
%! assert_equal (M.RSquared, 126 / 132, 1e-12);

## ExpandedFactorNames is populated for the anovan backend
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! a.fit ();
%! assert_equal (iscellstr (a.ExpandedFactorNames), true);
%! assert_equal (! isempty (a.ExpandedFactorNames), true);

## groupmeans(): confidence bounds bracket the group means
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! M = groupmeans (anova (g, y, 'SumOfSquaresType', 'two'));
%! assert_equal (all (M.MeanLower <= M.Mean, 'all'), true);
%! assert_equal (all (M.Mean <= M.MeanUpper, 'all'), true);

## boxchart(): returns a non-empty graphics result
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%!   g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%!   h = boxchart (anova (g, y));
%!   assert_equal (! isempty (h), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Weighted anovan fit exposes raw residuals (FittedValues + Residuals == Y)
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! w = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'Weights', w);
%! a.fit ();
%! assert_equal (a.FittedValues + a.Residuals.Raw, y, 1e-9);

## --- Week 5: diagnostic plots ------------------------------------------

## plotDiagnostics(): anovan-backed fit creates the four-panel figure
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%!   g = [1;1;1;2;2;2;3;3;3];
%!   a = anova (g, y, 'SumOfSquaresType', 'two');
%!   h = plotDiagnostics (a, 'Visible', 'off');
%!   assert_equal (all (ishghandle (h), 'all'), true);
%!   assert_equal (numel (findall (h, 'type', 'axes')), 4);
%! unwind_protect_cleanup
%!   close (h);
%!   close (hf);
%! end_unwind_protect

## --- Week 6: predict / effect sizes / LinearModel bridge ----------------

## predict(): no Xnew returns fitted values
%!test
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! assert_equal (predict (a), a.FittedValues, 1e-9);

## predict(): accepts an explicit design matrix for anovan-backed fits
%!test
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! a.fit ();
%! assert_equal (predict (a, full (a.DesignMatrix)), a.FittedValues, 1e-9);

## getEffectSizes(): anovan-backed fits report effect-size vectors
%!test
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! es = getEffectSizes (a);
%! assert_equal (iscell (es.Source), true);
%! assert_equal (numel (es.EtaSquared), numel (es.Source));
%! assert_equal (all (isfinite (es.PartialEtaSquared), 'all'), true);

## anova(LinearModel): builds from Avanish's LinearModel surface
%!test
%! X = [1; 2; 3; 4; 5];
%! y = [2; 4; 5; 4; 5];
%! mdl = fitlm (X, y);
%! a = anova (mdl);
%! assert_equal (! isempty (strfind (evalc ('disp (a)'), 'anova')), true);
%! assert_equal (! isempty (a.AnovaTable), true);
%! assert_equal (predict (a), mdl.Fitted, 1e-9);
%! es = getEffectSizes (a);
%! assert_equal (es.Source, {'Model'});
%! assert_equal (es.EtaSquared, mdl.SSR / mdl.SST, 1e-9);

## anova(LinearModel): plotDiagnostics uses LinearModel diagnostics
%!test
%! X = [1; 2; 3; 4; 5];
%! y = [2; 4; 5; 4; 5];
%! mdl = fitlm (X, y);
%! a = anova (mdl);
%! h = plotDiagnostics (a, 'Visible', 'off');
%! unwind_protect
%!   assert_equal (all (ishghandle (h), 'all'), true);
%!   assert_equal (numel (findall (h, 'type', 'axes')), 4);
%! unwind_protect_cleanup
%!   close (h);
%! end_unwind_protect

## --- Input validation ---------------------------------------------------

%!error <anova: too few input arguments.> anova ()
%!error <anova: Y must be a non-empty numeric array.> anova ('abc')
%!error <anova: Y must be a non-empty numeric array.> anova ([])
%!error <anova: name-value pairs must come in pairs.> ...
%! anova ([1;1;2;2], [1;2;3;4], 'SumOfSquaresType')
%!error <anova: parameter 'bogus' is not supported.> ...
%! anova ([1;1;2;2], [1;2;3;4], 'bogus', 1)
%!error <anova: SumOfSquaresType must be> ...
%! anova ([1;1;2;2], [1;2;3;4], 'SumOfSquaresType', 5)
%!error <anova: Alpha must be a numeric scalar in .0, 1..> ...
%! anova ([1;1;2;2], [1;2;3;4], 'Alpha', 2)
%!error <anova: Display must be 'on' or 'off'.> ...
%! anova ([1;1;2;2], [1;2;3;4], 'Display', 'maybe')
%!error <anova: parameter name must be a character vector.> ...
%! anova ([1;1;2;2], [1;2;3;4], 1, 2)
%!error <anova: Reps must be a positive integer scalar.> ...
%! anova (magic (4), [], 'reps', -2)
%!error <anova: Reps must be a positive integer scalar.> ...
%! anova (magic (4), [], 'reps', 1.5)
%!error <anova: ModelSpecification must be> ...
%! anova ([1;1;2;2], [1;2;3;4], 'ModelSpecification', 'quadratic')
%!error <anova: terms matrix must have one column per factor.> ...
%! anova ([1;1;2;2], [1;2;3;4], 'Model', [1 0])
%!error <anova: GROUP must match the number of observations.> ...
%! anova ([1;1;2], [1;2;3;4])
%!error <anova: GROUP variables must match the number of observations.> ...
%! anova ({[1;1;2], [1;2;1;2]}, [1;2;3;4])
%!error <anova: Weights must have one value per observation.> ...
%! anova ([1;1;2;2], [1;2;3;4], 'Weights', [1;1;1])
%!error <anova: CategoricalFactors must contain valid factor indices.> ...
%! anova ([1;1;2;2], [1;2;3;4], 'CategoricalFactors', 1.5)
%!error <anova: CategoricalFactors must contain valid factor indices.> ...
%! anova ([1;1;2;2], [1;2;3;4], 'CategoricalFactors', 2)

## --- Week 9: edge cases ------------------------------------------------

%!test
%! a = anova (ones (5, 1), (1:5)');
%! [T, ems] = stats (a);
%! assert_equal (T.DF(1), 0);
%! assert_equal (T.pValue(1), NaN);
%! assert_equal (cellstr (ems.ExpectedMeanSquares), ...
%!               {'V(Error)'; 'V(Error)'});
%! a = anova (ones (5, 1), (1:5)', 'SumOfSquaresType', 'two');
%! T = stats (a);
%! assert_equal (T.DF(1), 0);
%! assert_equal (T.pValue(1), NaN);

%!test
%! a = anova (1, 7);
%! T = stats (a);
%! assert_equal (a.NumObservations, 1);
%! assert_equal (T.pValue(1), NaN);

%!test
%! a = anova ([1; 1; 2; 2], ones (4, 1));
%! T = stats (a);
%! assert_equal (T.SumOfSquares(1), 0);
%! assert_equal (T.pValue(1), 1);

%!test
%! g = categorical ([3; 3; 1; 1], [3, 2, 1]);
%! a = anova (g, (1:4)');
%! T = stats (a);
%! assert_equal (isfinite (T.pValue(1)), true);
%! assert_equal (a.Stats.n, [2, 2]);
%! assert_equal (a.Stats.gnames, {'3'; '1'});
%! assert_equal (a.Stats.means, [1.5, 3.5]);
%! a = anova (g, (1:4)', 'SumOfSquaresType', 'two');
%! T = stats (a);
%! assert_equal (T.F(1), 8, 1e-12);
%! assert_equal (a.Stats.grpnames{1}, {'3'; '1'});

%!test
%! g = kron ((1:120)', ones (2, 1));
%! a = anova (g, (1:240)');
%! T = stats (a);
%! assert_equal (T.DF(1), 119);
%! assert_equal (T.DF(2), 120);

%!test
%! n = 128;
%! group = cell (1, 6);
%! for k = 1:6
%!   group{k} = mod (floor ((0:n-1)' / 2^(k-1)), 2) + 1;
%! endfor
%! a = anova (group, (1:n)', 'ModelSpecification', 'full');
%! stats (a);
%! assert_equal (rows (a.Stats.terms), 63);
%! assert_equal (columns (a.DesignMatrix), 64);

## --- Week 10: configuration integration -------------------------------

%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'FactorNames', {'Treatment'});
%! [component, ems] = stats (a);
%! assert_equal (istable (component), true);
%! assert_equal (component.Properties.VariableNames, ...
%!               {'SumOfSquares', 'DF', 'MeanSquares', 'F', 'pValue'});
%! assert_equal (component.Properties.RowNames, ...
%!               {'Treatment'; 'Error'; 'Total'});
%! assert_equal (component.SumOfSquares, [126; 6; 132], 1e-12);
%! assert_equal (istable (ems), true);
%! assert_equal (ems.Properties.VariableNames, ...
%!               {'Type', 'ExpectedMeanSquares', 'MeanSquaresDenominator', ...
%!                'DFDenominator', 'FDenominator'});
%! assert_equal (ems.Properties.RowNames, {'Treatment'; 'Error'});
%! assert_equal (cellstr (ems.Type), {'fixed'; 'random'});
%! assert_equal (cellstr (ems.ExpectedMeanSquares), ...
%!               {'3*Q(Treatment)+V(Error)'; 'V(Error)'});
%! assert_equal (ems.MeanSquaresDenominator, [1; NaN]);
%! assert_equal (ems.DFDenominator, [6; NaN]);
%! assert_equal (cellstr (ems.FDenominator), {'MS(Error)'; ''});

%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'FactorNames', {'Treatment'}, 'RandomFactors', 1);
%! [~, ems] = stats (a);
%! assert_equal (cellstr (ems.Type), {'random'; 'random'});
%! assert_equal (cellstr (ems.ExpectedMeanSquares), ...
%!               {'3*V(Treatment)+V(Error)'; 'V(Error)'});

%!test
%! y = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!      6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (y, [], 'Reps', 3, 'ModelSpecification', 'interactions', ...
%!            'FactorNames', {'Brand', 'PopperType'});
%! component = stats (a);
%! assert_equal (component.Properties.RowNames, ...
%!               {'Brand'; 'PopperType'; 'Brand:PopperType'; 'Error'; 'Total'});

%!function values = __anova_values__ (tbl)
%!  if (istable (tbl))
%!    values = [tbl.SumOfSquares, tbl.DF, tbl.MeanSquares, tbl.F, tbl.pValue];
%!    return;
%!  endif
%!  if (columns (tbl) == 6)
%!    columns_ = [2, 3, 4, 5, 6];
%!  else
%!    columns_ = [2, 3, 4, 6, 7];
%!  endif
%!  values = NaN (rows (tbl) - 1, numel (columns_));
%!  for i = 2:rows (tbl)
%!    for j = 1:numel (columns_)
%!      value = tbl{i, columns_(j)};
%!      if (isnumeric (value) && isscalar (value))
%!        values(i - 1, j) = value;
%!      endif
%!    endfor
%!  endfor
%!endfunction

%!test
%! y = [24; 26; 25; 24; 15; 17; 20; 16; 25; 29; 27; 19; 18; 21; 20];
%! gender = [1; 1; 1; 1; 1; 1; 1; 1; 2; 2; 2; 2; 2; 2; 2];
%! degree = [1; 1; 1; 1; 0; 0; 0; 0; 1; 1; 1; 0; 0; 0; 0];
%! a = anova ({gender, degree}, y, 'ModelSpecification', 'full');
%! component = stats (a, 'component', 'one');
%! [~, expected] = anovan (y, {gender, degree}, 'model', 'full', ...
%!                         'sstype', 1, 'display', 'off');
%! assert_equal (component.Properties.RowNames, ...
%!               {'Factor1'; 'Factor2'; 'Factor1:Factor2'; 'Error'; 'Total'});
%! assert_equal (__anova_values__ (component), ...
%!               __anova_values__ (expected), 1e-10);
%! assert_equal (a.SumOfSquaresType, 'three');
%! assert_equal (istable (stats (a)), true);

%!test
%! y = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!      6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (y, [], 'Reps', 3, 'ModelSpecification', 'interactions');
%! component = stats (a, 'Component', 'one');
%! assert_equal (__anova_values__ (component), ...
%!               __anova_values__ (stats (a)), 1e-10);
%! summary_table = stats (a, 'summary');
%! assert_equal (summary_table.Properties.RowNames, ...
%!               {'Linear'; 'NonLinear'; 'Regression'; 'Error'; 'Total'});
%! assert_equal (summary_table.SumOfSquares(1:3), ...
%!               [20.25; 1 / 12; 61 / 3], 1e-10);
%! assert_equal (summary_table.DF(1:3), [3; 2; 5]);
%! assert_equal (summary_table.F(1:3), [48.6; 0.3; 29.28], 1e-10);

%!test
%! x = kron ((0:3)', ones (2, 1));
%! y = [0; 0.2; 1; 1.2; 4; 4.2; 9; 9.2];
%! summary_table = stats (anova (x, y, 'CategoricalFactors', []), 'summary');
%! assert_equal (summary_table.Properties.RowNames, ...
%!               {'Linear'; 'Regression'; 'Error'; 'LackOfFit'; ...
%!                'PureError'; 'Total'});
%! assert_equal (summary_table.SumOfSquares(3:5), [8.08; 8; 0.08], 1e-12);
%! assert_equal (summary_table.DF(3:5), [6; 2; 4]);
%! assert_equal (summary_table.F(4), 200, 1e-10);
%! assert_equal (summary_table.pValue(4), 9.802960494567e-05, 1e-12);
%! g = [ones(10, 1); 2 * ones(10, 1)];
%! y = [(0:9)'; (1000:1009)'];
%! summary_table = stats (anova (g, y), 'summary');
%! assert_equal (summary_table.pValue(1) > 0, true);

%!test
%! [demo_code, demo_idx] = test ('anovan', 'grabdemo');
%! assert_equal (numel (demo_idx) - 1, 13);
%! for k = 1:numel (demo_idx) - 1
%!   code = demo_code(demo_idx(k):demo_idx(k + 1) - 1);
%!   code = strrep (code, "'display', 'on'", "'display', 'off'");
%!   code = strrep (code, "anovan (y, g, 'weights', v.^-1)", ...
%!                         "anovan (y, g, 'weights', v.^-1, 'display', 'off')");
%!   code = regexprep (code, '^ *(figure|plot|xlabel) [^\n]*$', '', 'lineanchors');
%!   eval (code);
%!   switch (k)
%!     case 1
%!       a = anova (gender, score, 'FactorNames', {'gender'});
%!     case 2
%!       a = anova ({treatment(:), subject(:)}, score(:), ...
%!                  'ModelSpecification', 'full', 'RandomFactors', 2, ...
%!                  'SumOfSquaresType', 'two', ...
%!                  'FactorNames', {'treatment', 'subject'});
%!     case 3
%!       a = anova (alloy, strength, 'FactorNames', {'alloy'});
%!     case 4
%!       a = anova ({seconds(:), subject(:)}, words(:), ...
%!                  'ModelSpecification', 'full', 'RandomFactors', 2, ...
%!                  'SumOfSquaresType', 'two', ...
%!                  'FactorNames', {'seconds', 'subject'});
%!     case 5
%!       a = anova ({brands(:), popper(:)}, popcorn(:), ...
%!                  'ModelSpecification', 'full', ...
%!                  'FactorNames', {'brands', 'popper'});
%!     case 6
%!       a = anova ({gender, degree}, salary, 'ModelSpecification', 'full', ...
%!                  'FactorNames', {'gender', 'degree'});
%!     case 7
%!       a = anova ({sugar, milk}, babble, 'ModelSpecification', 'full', ...
%!                  'FactorNames', {'sugar', 'milk'});
%!     case 8
%!       a = anova ({drug(:), feedback(:), diet(:)}, BP(:), ...
%!                  'ModelSpecification', 'full', ...
%!                  'FactorNames', {'drug', 'feedback', 'diet'});
%!     case 9
%!       a = anova ({strain, treatment, block}, measurement / 10, ...
%!                  'ModelSpecification', 'full', 'RandomFactors', 3, ...
%!                  'SumOfSquaresType', 'two', ...
%!                  'FactorNames', {'strain', 'treatment', 'block'});
%!     case 10
%!       a = anova ({species, temp}, pulse, 'CategoricalFactors', 1, ...
%!                  'SumOfSquaresType', 'hierarchical', ...
%!                  'FactorNames', {'species', 'temp'});
%!     case 11
%!       model = [1 0 0; 0 1 0; 0 0 1; 1 1 0];
%!       a = anova ({treatment, exercise, age}, score, ...
%!                  'ModelSpecification', model, 'CategoricalFactors', [1, 2], ...
%!                  'SumOfSquaresType', 'hierarchical', ...
%!                  'FactorNames', {'treatment', 'exercise', 'age'});
%!     case 12
%!       a = anova (g, dv, 'FactorNames', {'score'});
%!     case 13
%!       a = anova (g, y, 'Weights', v .^ -1);
%!   endswitch
%!   actual = stats (a);
%!   actual_sources = actual.Properties.RowNames;
%!   expected_sources = regexprep (ATAB(2:end, 1), "'$", "");
%!   expected_sources = regexprep (expected_sources, "\\*", ":");
%!   expected_sources = regexprep (expected_sources, "X([0-9]+)", "Factor$1");
%!   if (strcmp (a.Stats.source, 'anova1'))
%!     actual_sources{1} = expected_sources{1};
%!   endif
%!   assert_equal (actual_sources, expected_sources);
%!   assert_equal (__anova_values__ (actual), ...
%!                 __anova_values__ (ATAB), 1e-8);
%!   if (! any (k == [1, 3, 12]))
%!     assert_equal (a.Residuals, STATS.resid, 1e-8);
%!     assert_equal (size (a.DesignMatrix), size (STATS.X));
%!   endif
%! endfor

%!error <anova: RandomFactors must contain positive integer indices.> ...
%! anova ([1;1;2;2], [1;2;3;4], 'RandomFactors', 0)
%!error <anova: RandomFactors indices exceed the number of factors.> ...
%! anova ([1;1;2;2], [1;2;3;4], 'RandomFactors', 2)
%!error <anova.stats: type must be a character vector.> ...
%! stats (anova ([1;1;1;2;2;2;3;3;3], (1:9)'), 5)
%!error <anova.groupmeans: factors are required for group means.> ...
%! groupmeans (anova ([1;2;3;4;5;6]))
%!error <anova.varianceComponent: random factors are not implemented.> ...
%! varianceComponent (anova ([1;1;1;2;2;2;3;3;3], (1:9)', 'RandomFactors', 1))
%!error <diagnostic plots require>
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! plotDiagnostics (anova (popcorn, [], 'reps', 3));
%!error <name-value pairs must come in pairs>
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! plotDiagnostics (anova (g, y, 'SumOfSquaresType', 'two'), 'Visible');
%!error <unknown option>
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! plotDiagnostics (anova (g, y, 'SumOfSquaresType', 'two'), 'BadOption', true);
%!error <Xnew must be a numeric design matrix>
%! y = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g = [1;1;1;2;2;2;3;3;3];
%! predict (anova (g, y, 'SumOfSquaresType', 'two'), ones (2, 2));

%!error <component type must be a character vector>
%! a = anova ([1; 1; 2; 2], (1:4)');
%! stats (a, 'Component', 1);

%!error <component type must be 'one'>
%! a = anova ([1; 1; 2; 2], (1:4)');
%! stats (a, 'Component', 'typei');

%!error <type must be 'component' or 'summary'>
%! a = anova ([1; 1; 2; 2], (1:4)');
%! stats (a, 'details');

%!error <invalid input arguments>
%! a = anova ([1; 1; 2; 2], (1:4)');
%! stats (a, 'summary', 'one');
