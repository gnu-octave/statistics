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
  ## @seealso{anova1, anova2, anovan, multcompare}
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
    ## Table containing one variable for each factor used to fit the ANOVA
    ## model.  This property is read-only.
    ##
    ## @end deftp
    Factors

    ## -*- texinfo -*-
    ## @deftp {anova} {property} Formula
    ##
    ## Model formula
    ##
    ## Read-only structural formula value with response, predictor, term,
    ## nesting, and linear-predictor fields matching MATLAB's formula object.
    ##
    ## @end deftp
    Formula         = struct ();

    ## -*- texinfo -*-
    ## @deftp {anova} {property} FactorNames
    ##
    ## Factor names
    ##
    ## String row vector containing the factor names used by the fitted ANOVA
    ## model.  This property is read-only.
    ##
    ## @end deftp
    FactorNames     = [];

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
    ## String scalar selecting @qcode{"one"}, @qcode{"two"},
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
    ## Positive integer indices of factors treated as categorical.  This
    ## property is read-only.
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
    ## Numeric vector of fitted coefficient estimates.  This property is
    ## read-only.
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
    Response     = [];
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
    fitted_     = false;
    dirty_      = true;
    nFactors_   = 0;
    backend_    = '';                       ## 'anova1' | 'anova2' | 'anovan'
    reps_       = [];                       ## replicate count for anova2 backend
    continuousSpecified_ = false;
    coefficientStats_ = [];
    rawResiduals_ = [];
    nesting_ = [];
    formulaText_ = '';
    formulaTerms_ = [];
  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {@var{obj} =} anova (@var{Y})
    ## @deftypefnx {anova} {@var{obj} =} anova (@var{factors}, @var{Y})
    ## @deftypefnx {anova} {@var{obj} =} anova (@var{tbl}, @var{Y})
    ## @deftypefnx {anova} {@var{obj} =} anova (@var{tbl}, @var{responseVarName})
    ## @deftypefnx {anova} {@var{obj} =} anova (@var{tbl}, @var{formula})
    ## @deftypefnx {anova} {@var{obj} =} anova (@dots{}, @var{name}, @var{value})
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
    ## @var{tbl} is a table whose variables contain factors.  The response may
    ## be supplied separately, selected by variable name, or specified with a
    ## Wilkinson formula.  @qcode{'FactorNames'} can select a subset of table
    ## variables when a formula is not supplied.
    ##
    ## Supported name-value arguments include @qcode{'ModelSpecification'},
    ## @qcode{'SumOfSquaresType'}, @qcode{'FactorNames'},
    ## @qcode{'CategoricalFactors'}, @qcode{'RandomFactors'},
    ## @qcode{'ResponseName'}, @qcode{'Alpha'}, and @qcode{'Display'}.
    ## Passing @qcode{'Reps'} selects the balanced two-way @code{anova2}
    ## backend when @var{Y} is a non-vector matrix.
    ##
    ## @end deftypefn
    function obj = anova (factors, Y, varargin)

      if (nargin < 1)
        error ("anova: too few input arguments.");
      endif
      if (nargin < 2)
        Y = [];
      endif

      table_factors = [];
      formula_input = "";
      response_locked = false;

      if (isa (factors, "table"))
        if (nargin < 2)
          error (strcat ("anova: table input requires response data or", ...
                         " a response name."));
        endif
        if (mod (numel (varargin), 2) != 0)
          error ("anova: name-value pairs must come in pairs.");
        endif
        response_locked = obj.isName_ (Y);
        [factors, Y, table_factors, table_names, response_name, ...
         formula_input, formula_model, formula_nesting] = ...
          obj.parseTableInput_ (factors, Y, varargin{:});
        obj.VarNames = table_names;
        obj.FactorNames = string (table_names);
        obj.ResponseName = response_name;
        if (! isempty (formula_input))
          obj.formulaText_ = formula_input;
          obj.formulaTerms_ = formula_model;
          obj.ModelSpecification = formula_input;
          obj.ModelType = formula_model;
          obj.nesting_ = formula_nesting;
        endif
      elseif (nargin < 2)
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

      obj.Response = Y;
      obj.Y        = Y(:);
      obj.GROUP    = factors;

      ## Parse name-value pairs (mirrors anovan.m's loop style)
      for idx = 1:2:numel (varargin)
        name  = varargin{idx};
        value = varargin{idx + 1};
        if (! ischar (name))
          error ("anova: parameter name must be a character vector.");
        endif
        switch (lower (name))
          case {'model', 'modelspecification'}
            if (isempty (formula_input))
              obj = obj.setModelSpecification_ (value);
            endif
          case {'sstype', 'sumofsquarestype'}
            obj = obj.setSumOfSquaresType_ (value);
          case {'varnames', 'factornames'}
            if (isempty (formula_input))
              obj = obj.setFactorNames_ (value);
            endif
          case 'contrasts'
            obj.Contrasts = value;
          case 'alpha'
            obj.Alpha = value;
          case {'categoricalfactors'}
            obj.CategoricalFactors = value;
          case {'continuous'}
            obj.Continuous = value;
            obj.continuousSpecified_ = true;
          case {'random', 'randomfactors'}
            obj.RandomFactors = value;
            obj.Random = value;
          case 'weights'
            obj.Weights = value;
          case {'display', 'displayopt'}
            obj.Display = value;
          case 'responsename'
            if (! response_locked)
              obj.ResponseName = value;
            endif
          case 'reps'
            obj.reps_ = value;
          otherwise
            error ("anova: parameter '%s' is not supported.", name);
        endswitch
      endfor

      [obj.Response, obj.GROUP, obj.Weights, table_factors] = ...
        obj.removeMissing_ (obj.Response, obj.GROUP, obj.Weights, ...
                            table_factors);
      obj.GROUP = obj.normalizeGroups_ (obj.GROUP);
      obj.Y = obj.Response(:);
      obj.nFactors_ = obj.countFactors_ ();
      obj.NumFactors = obj.nFactors_;
      obj.NumObservations = numel (obj.Response);
      if (isnumeric (obj.ModelType) && isscalar (obj.ModelType))
        obj.ModelType = min (obj.ModelType, obj.nFactors_);
      endif
      if (isempty (formula_input))
        obj = obj.setFormula_ ();
      endif
      obj.FactorNames = string (obj.VarNames);
      obj.SumOfSquaresType = string (obj.SumOfSquaresType);
      obj = obj.updateFormula_ ();
      obj = obj.syncFactorSelectors_ ();
      obj.validateSpec_ ();
      obj.validateData_ ();
      if (isempty (table_factors))
        obj.Factors = obj.factorTable_ ();
      else
        obj.Factors = table_factors;
      endif
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
      if (! isempty (obj.formulaText_))
        fprintf ("  %s\n\n", obj.formulaText_);
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
        if (obj.NumFactors != 1)
          error (strcat ("anova.groupmeans: factors must be specified for", ...
                         " a multi-factor ANOVA."));
        endif
        factors = obj.VarNames;
      endif
      alpha = obj.parseAlpha_ (varargin{:});
      obj.ensureFit_ ();
      idx = obj.factorIndices_ (factors);
      names = obj.VarNames(idx);
      y = obj.Y(:);
      if (isempty (idx))
        error ("anova.groupmeans: factors are required for group means.");
      endif
      [gid, levels] = findgroups (obj.Factors(:, names));
      valid = isfinite (y) & isfinite (gid) & gid > 0;
      used = unique (gid(valid), "stable");
      [~, compact] = ismember (gid(valid), used);
      levels = levels(used, :);
      n = accumarray (compact, 1);
      mu = accumarray (compact, y(valid), [], @mean);
      se = sqrt (obj.MSE ./ n);
      crit = tinv (1 - alpha / 2, max (obj.DFE, 1));
      lower = mu - crit * se;
      upper = mu + crit * se;
      intervals = table (mu, se, lower, upper, "VariableNames", ...
                         {"Mean", "SE", "MeanLower", "MeanUpper"});
      means = [levels, intervals];
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {} boxchart (@var{obj})
    ## @deftypefnx {anova} {@var{h} =} boxchart (@var{obj}, @dots{})
    ##
    ## Plot response values grouped by up to two categorical factors.
    ##
    ## This method uses @code{boxplot} as the graphics backend in Octave and
    ## returns the native box graphics handles.  A target axes may be supplied
    ## as the first optional argument.
    ##
    ## @end deftypefn
    function h = boxchart (obj, varargin)
      obj.ensureFit_ ();
      [target, varargin] = obj.extractAxes_ (varargin);
      factors = {};
      if (! isempty (varargin))
        candidate = varargin{1};
        if (isstring (candidate))
          candidate = cellstr (candidate(:))';
        elseif (ischar (candidate))
          candidate = {candidate};
        endif
        if (iscellstr (candidate) ...
            && all (ismember (candidate, obj.VarNames)))
          factors = candidate;
          varargin(1) = [];
        endif
      endif
      if (isempty (factors))
        if (obj.NumFactors != 1)
          error (strcat ("anova.boxchart: factors must be specified for", ...
                         " a multi-factor ANOVA."));
        endif
        factors = obj.VarNames;
      endif
      idx = obj.factorIndices_ (factors);
      if (numel (idx) > 2)
        error ("anova.boxchart: at most two factors can be plotted.");
      endif
      if (! all (ismember (idx, obj.CategoricalFactors)))
        error ("anova.boxchart: factors must be categorical.");
      endif
      grouping = cell (1, numel (idx));
      for k = 1:numel (idx)
        grouping{k} = obj.Factors{:, obj.VarNames{idx(k)}};
      endfor
      if (numel (grouping) == 1)
        grouping = grouping{1};
      endif
      if (! isempty (target))
        axes (target);
      endif
      [~, handles] = boxplot (obj.Y(:), grouping, "notch", "on", ...
                              varargin{:});
      h = handles.box;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {} plotComparisons (@var{obj})
    ## @deftypefnx {anova} {@var{h} =} plotComparisons (@var{obj}, @dots{})
    ##
    ## Plot multiple-comparison intervals for model-adjusted group means.
    ##
    ## Clicking a group highlights it and distinguishes groups whose adjusted
    ## comparison is significant at the requested alpha level.  A target axes
    ## may be supplied as the first optional argument.
    ##
    ## @end deftypefn
    function h = plotComparisons (obj, varargin)
      obj.ensureFit_ ();
      [target, varargin] = obj.extractAxes_ (varargin);
      [stats_, args, sidak] = obj.comparisonArguments_ (varargin{:});
      args = obj.replaceOption_ (args, "display", "off");
      alpha = obj.optionValue_ (args, "alpha", obj.Alpha);
      [comparisons, means, ~, names] = multcompare (stats_, args{:});
      if (! isempty (sidak) && ! isempty (comparisons))
        per_comparison = 1 - (1 - sidak) ^ (1 / rows (comparisons));
        critical = tinv (1 - per_comparison / 2, comparisons(1, 8));
        half_width = means(:, 2) * critical / sqrt (2);
        means(:, 3:4) = [means(:, 1) - half_width, ...
                         means(:, 1) + half_width];
        comparisons(:, 6) = 1 - (1 - comparisons(:, 6)) ...
                                .^ rows (comparisons);
      endif
      if (isempty (target))
        target = gca ();
      endif
      h = ancestor (target, "figure");
      cla (target);
      hold (target, "on");
      count = rows (means);
      intervals = zeros (count, 1);
      markers = zeros (count, 1);
      for group = 1:count
        intervals(group) = line (target, means(group, 3:4), ...
                                 [group, group], "color", [0.35, 0.35, 0.35], ...
                                 "linewidth", 1.5);
        markers(group) = line (target, means(group, 1), group, ...
                               "linestyle", "none", "marker", "o", ...
                               "markerfacecolor", [0.00, 0.45, 0.74], ...
                               "markeredgecolor", [0.00, 0.45, 0.74]);
      endfor
      hold (target, "off");
      set (target, "ydir", "reverse", "ytick", 1:count, ...
                   "yticklabel", names);
      ylim (target, [0.5, count + 0.5]);
      xlabel (target, sprintf ("Group mean with %g%% comparison interval", ...
                               100 * (1 - alpha)));
      title (target, "Multiple comparisons");
      for group = 1:count
        set (markers(group), "buttondownfcn", ...
             @(~, ~) obj.highlightComparison_ (group, markers, intervals, ...
                                                comparisons, alpha));
      endfor
      if (count > 0)
        obj.highlightComparison_ (1, markers, intervals, comparisons, alpha);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {@var{v} =} varianceComponent (@var{obj})
    ## @deftypefnx {anova} {@var{v} =} varianceComponent (@var{obj}, @dots{})
    ##
    ## Return variance component estimates for random model terms and error.
    ##
    ## @end deftypefn
    function v = varianceComponent (obj, varargin)
      alpha = obj.parseAlpha_ (varargin{:});
      obj.ensureFit_ ();
      [coefficients, mean_squares, dfs, names] = ...
        obj.varianceSystem_ (obj.Stats, obj.AnovaTable, obj.SSType);
      estimates = coefficients \ mean_squares;
      lower_ms = dfs .* mean_squares ./ chi2inv (1 - alpha / 2, dfs);
      upper_ms = dfs .* mean_squares ./ chi2inv (alpha / 2, dfs);
      inverse = pinv (coefficients);
      lower = sum (max (inverse, 0) .* lower_ms' ...
                   + min (inverse, 0) .* upper_ms', 2);
      upper = sum (max (inverse, 0) .* upper_ms' ...
                   + min (inverse, 0) .* lower_ms', 2);
      ## A variance cannot be negative, so the propagated lower bound is
      ## truncated at zero.  A negative estimate has no interval at all.
      lower = max (lower, 0);
      undefined = (estimates < 0);
      lower(undefined) = NaN;
      upper(undefined) = NaN;
      v = table (estimates, lower, upper, "VariableNames", ...
                 {"VarianceComponent", "VarianceComponentLower", ...
                  "VarianceComponentUpper"}, "RowNames", names);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {anova} {@var{m} =} multcompare (@var{obj})
    ## @deftypefnx {anova} {@var{m} =} multcompare (@var{obj}, @var{factors})
    ## @deftypefnx {anova} {@var{m} =} multcompare (@dots{}, @var{name}, @var{value})
    ##
    ## Perform post-hoc multiple comparisons for a fitted ANOVA object.
    ##
    ## The returned table contains the compared groups, estimated mean
    ## difference, confidence limits, and p-value.  The default critical value
    ## type is @qcode{"tukey-kramer"}.
    ##
    ## @end deftypefn
    function m = multcompare (obj, varargin)
      obj.ensureFit_ ();
      if (isempty (fieldnames (obj.Stats)))
        error ("anova.multcompare: model has no stats to compare.");
      endif
      [stats_, args, sidak, dims] = obj.comparisonArguments_ (varargin{:});
      gid = findgroups (obj.Factors(:, obj.VarNames(dims)));
      if (numel (unique (gid(isfinite (gid) & gid > 0))) < 2)
        m = obj.comparisonTable_ (zeros (0, 6), cell (0, 1), dims);
        return;
      endif
      [C, ~, ~, group_names] = multcompare (stats_, args{:});
      if (! isempty (sidak))
        per_comparison = 1 - (1 - sidak) ^ (1 / rows (C));
        old_critical = tinv (1 - sidak / 2, C(:, 8));
        standard_error = (C(:, 5) - C(:, 3)) ./ (2 * old_critical);
        critical = tinv (1 - per_comparison / 2, C(:, 8));
        C(:, 3) = C(:, 4) - critical .* standard_error;
        C(:, 5) = C(:, 4) + critical .* standard_error;
        C(:, 6) = 1 - (1 - C(:, 6)) .^ rows (C);
      endif
      if (strcmp (stats_.source, "anovan"))
        group_names = obj.comparisonGroups_ (stats_, dims);
      endif
      m = obj.comparisonTable_ (C, group_names, dims);
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
    ## plot.  Diagnostic plots require an @code{anovan}-backed fit because the
    ## fast @code{anova1} and @code{anova2} backends do not expose residuals
    ## and design-matrix diagnostics.
    ##
    ## Supported name-value arguments are @qcode{'FigureName'} and
    ## @qcode{'Visible'}.
    ##
    ## @end deftypefn
    function h = plotDiagnostics (obj, varargin)
      obj.ensureFit_ ();
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
    ## with one column per coefficient.
    ##
    ## @end deftypefn
    function ypred = predict (obj, Xnew, varargin)
      obj.ensureFit_ ();
      if (nargin < 2 || isempty (Xnew))
        ypred = obj.FittedValues;
        return;
      endif
      if (isempty (obj.Coefficients))
        error ("anova.predict: coefficients are unavailable for this backend.");
      endif
      if (! isempty (obj.coefficientStats_))
        beta = obj.coefficientStats_(:, 1);
      else
        beta = obj.Coefficients;
      endif
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
    ## Values are derived from the fitted ANOVA table.
    ##
    ## @end deftypefn
    function es = getEffectSizes (obj)
      obj.ensureFit_ ();
      es = obj.effectSizesFromAtab_ ();
    endfunction

  endmethods

  methods (Access = private)

    function tf = isName_ (obj, value)
      tf = ischar (value) || (isstring (value) && isscalar (value));
    endfunction

    function [target, args] = extractAxes_ (obj, args)
      target = [];
      if (! isempty (args) && isscalar (args{1}) ...
          && ishghandle (args{1}) ...
          && strcmp (get (args{1}, "type"), "axes"))
        target = args{1};
        args(1) = [];
      endif
    endfunction

    function highlightComparison_ (obj, selected, markers, intervals, ...
                                   comparisons, alpha)
      neutral = [0.45, 0.45, 0.45];
      different = [0.85, 0.20, 0.16];
      chosen = [0.00, 0.45, 0.74];
      for group = 1:numel (markers)
        color = neutral;
        if (group == selected)
          color = chosen;
        else
          row = find ((comparisons(:, 1) == selected ...
                       & comparisons(:, 2) == group) ...
                      | (comparisons(:, 2) == selected ...
                         & comparisons(:, 1) == group), 1);
          if (! isempty (row) && comparisons(row, 6) < alpha)
            color = different;
          endif
        endif
        set (markers(group), "markerfacecolor", color, ...
             "markeredgecolor", color);
        set (intervals(group), "color", color);
      endfor
    endfunction

    function [groups, response, factor_table, factor_names, response_name, ...
              formula, model, nesting] = parseTableInput_ (obj, tbl, ...
                                                            response_arg, ...
                                                            varargin)
      variable_names = tbl.Properties.VariableNames;
      formula = "";
      model = [];
      nesting = [];
      if (isnumeric (response_arg))
        response = response_arg;
        response_name = "Y";
        factor_names = variable_names;
      elseif (obj.isName_ (response_arg))
        value = char (response_arg);
        if (! isempty (strfind (value, "~")))
          formula = value;
          schema = parseWilkinsonFormula (formula, "matrix");
          if (! isscalar (schema.ResponseIdx))
            error ("anova: formula must specify one response variable.");
          endif
          response_name = schema.VariableNames{schema.ResponseIdx};
          [factor_names, model] = obj.schemaModel_ (schema);
          ordered = variable_names(ismember (variable_names, factor_names));
          [~, order] = ismember (ordered, factor_names);
          factor_names = ordered;
          model = model(:, order);
          [model, nesting] = obj.nestedSpecification_ (formula, ...
                                                       factor_names, model);
        else
          response_name = value;
          factor_names = variable_names(! strcmp (variable_names, value));
        endif
        if (! any (strcmp (variable_names, response_name)))
          error ("anova: response variable '%s' is not in the table.", ...
                 response_name);
        endif
        response = tbl{:, response_name};
      else
        error (strcat ("anova: table response must be numeric, a response", ...
                       " variable name, or a formula."));
      endif

      if (! isempty (formula))
        missing = ! ismember (factor_names, variable_names);
        if (any (missing))
          error ("anova: formula variable '%s' is not in the table.", ...
                 factor_names{find (missing, 1)});
        endif
      else
        selected = obj.optionValue_ (varargin, "factornames", []);
        if (! isempty (selected))
          if (isstring (selected))
            selected = cellstr (selected(:))';
          elseif (ischar (selected))
            selected = {selected};
          endif
          if (! iscellstr (selected) ...
              || ! all (ismember (selected, factor_names)))
            error (strcat ("anova: FactorNames must identify factor", ...
                           " variables in the table."));
          endif
          factor_names = selected;
        endif
      endif
      if (! isnumeric (response) || isempty (response) || ! isvector (response))
        error (strcat ("anova: table response data must be a non-empty", ...
                       " numeric vector."));
      endif
      if (numel (response) != height (tbl))
        error (strcat ("anova: table and response must contain the same", ...
                       " number of rows."));
      endif

      factor_table = tbl(:, factor_names);
      groups = cell (1, numel (factor_names));
      for k = 1:numel (factor_names)
        groups{k} = tbl{:, factor_names{k}};
        if (ischar (groups{k}) && rows (groups{k}) > 1)
          groups{k} = cellstr (groups{k});
        endif
      endfor
      if (numel (groups) == 1)
        groups = groups{1};
      endif
    endfunction

    function [factor_names, model] = schemaModel_ (obj, schema)
      predictor_idx = setdiff (1:numel (schema.VariableNames), ...
                               schema.ResponseIdx, "stable");
      used = any (schema.Terms(:, predictor_idx) != 0, 1);
      predictor_idx = predictor_idx(used);
      factor_names = {};
      exponents = zeros (1, numel (predictor_idx));
      factor_idx = zeros (1, numel (predictor_idx));
      for k = 1:numel (predictor_idx)
        name = schema.VariableNames{predictor_idx(k)};
        token = regexp (name, '^(.*)\^(\d+)$', 'tokens', 'once');
        if (isempty (token))
          base = name;
          exponents(k) = 1;
        else
          base = token{1};
          exponents(k) = str2double (token{2});
        endif
        idx = find (strcmp (factor_names, base), 1);
        if (isempty (idx))
          factor_names{end + 1} = base;
          idx = numel (factor_names);
        endif
        factor_idx(k) = idx;
      endfor
      model = zeros (rows (schema.Terms), numel (factor_names));
      for k = 1:numel (predictor_idx)
        rows_ = schema.Terms(:, predictor_idx(k)) != 0;
        model(rows_, factor_idx(k)) = exponents(k);
      endfor
      model(! any (model, 2), :) = [];
      model = unique (model, "rows", "stable");
    endfunction

    function [model, nesting] = nestedSpecification_ (obj, formula, ...
                                                       factor_names, model)
      nesting = false (numel (factor_names));
      definitions = regexp (formula, ...
                            '([A-Za-z]\w*)\s*\(([^()]*)\)', ...
                            'tokens');
      for k = 1:numel (definitions)
        child_name = definitions{k}{1};
        parent_names = strtrim (strsplit (definitions{k}{2}, ','));
        child = find (strcmp (factor_names, child_name), 1);
        [found, parents] = ismember (parent_names, factor_names);
        if (isempty (child) || ! all (found))
          error ("anova: nested formula contains an unknown factor.");
        endif
        nesting(child, parents) = true;
        members = [child, parents];
        rows_ = all (model(:, members) > 0, 2) ...
                & sum (model > 0, 2) == numel (members);
        replacement = zeros (1, columns (model));
        replacement(child) = 1;
        model(rows_, :) = repmat (replacement, sum (rows_), 1);
      endfor
      model = unique (model, "rows", "stable");
    endfunction

    function value = optionValue_ (obj, args, option, default)
      value = default;
      for k = 1:2:numel (args)
        if (obj.isName_ (args{k}) && strcmpi (char (args{k}), option))
          value = args{k + 1};
        endif
      endfor
    endfunction

    function factors = factorTable_ (obj)
      if (obj.nFactors_ == 0)
        factors = table ();
        return;
      endif
      if (isempty (obj.GROUP) && ! isvector (obj.Response))
        if (isempty (obj.reps_))
          [n, groups] = size (obj.Response);
          group_values = reshape (repmat (1:groups, n, 1), [], 1);
          columns_ = {group_values};
        else
          [~, columns_] = obj.anova2Data_ ();
        endif
      elseif (iscell (obj.GROUP) ...
              && size (obj.GROUP, 1) == 1 ...
              && numel (obj.GROUP) == obj.nFactors_)
        columns_ = cellfun (@(x) x(:), obj.GROUP, "UniformOutput", false);
      elseif (obj.nFactors_ == 1)
        columns_ = {obj.GROUP(:)};
      else
        columns_ = cell (1, obj.nFactors_);
        for k = 1:obj.nFactors_
          columns_{k} = obj.GROUP(:, k);
        endfor
      endif
      factors = table (columns_{:}, "VariableNames", obj.VarNames);
    endfunction

    function [response, groups, weights, factor_table] = removeMissing_ (...
                                      obj, response, groups, weights, ...
                                      factor_table)
      if (! isvector (response))
        if (! isempty (obj.reps_))
          return;
        elseif (isempty (groups) && any (! isfinite (response(:))))
          [n, levels] = size (response);
          groups = reshape (repmat (1:levels, n, 1), [], 1);
          response = response(:);
        else
          return;
        endif
      endif
      valid = isfinite (response(:));
      columns_ = obj.groupColumns_ (groups);
      for k = 1:numel (columns_)
        group_id = grp2idx (columns_{k});
        if (numel (group_id) != numel (valid))
          return;
        endif
        valid &= isfinite (group_id);
      endfor
      if (all (valid))
        return;
      endif
      if (! any (valid))
        error (strcat ("anova: no complete observations remain after", ...
                       " removing missing data."));
      endif
      response = response(valid);
      if (iscell (groups) && size (groups, 1) == 1 ...
          && numel (groups) == numel (columns_))
        for k = 1:numel (groups)
          value = groups{k};
          if (ischar (value) && rows (value) > 1)
            groups{k} = value(valid, :);
          else
            groups{k} = value(valid);
          endif
        endfor
      elseif (! isempty (groups))
        groups = groups(valid, :);
      endif
      if (! isempty (weights) && numel (weights) != numel (valid))
        return;
      elseif (! isempty (weights))
        weights = weights(valid);
      endif
      if (isa (factor_table, "table") && width (factor_table) > 0)
        factor_table = factor_table(valid, :);
      endif
    endfunction

    function columns_ = groupColumns_ (obj, groups)
      if (isempty (groups))
        columns_ = {};
      elseif (iscell (groups) && size (groups, 1) == 1 ...
              && all (cellfun (@(x) isvector (x) || ischar (x), groups)))
        columns_ = groups;
      elseif (isvector (groups) || ischar (groups))
        columns_ = {groups};
      else
        columns_ = cell (1, columns (groups));
        for k = 1:columns (groups)
          columns_{k} = groups(:, k);
        endfor
      endif
    endfunction

    function groups = normalizeGroups_ (obj, groups)
      if (ischar (groups) && rows (groups) > 1)
        groups = cellstr (groups);
      elseif (iscell (groups) && size (groups, 1) == 1)
        for k = 1:numel (groups)
          if (ischar (groups{k}) && rows (groups{k}) > 1)
            groups{k} = cellstr (groups{k});
          endif
        endfor
      endif
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
        obj.SumOfSquaresType = string (names{value});
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
        case {"one", "typei", "i"}
          obj.SumOfSquaresType = string ("one");
          obj.SSType = 1;
        case {"two", "typeii", "ii"}
          obj.SumOfSquaresType = string ("two");
          obj.SSType = 2;
        case "hierarchical"
          obj.SumOfSquaresType = string ("hierarchical");
          obj.SSType = "h";
        case {"three", "typeiii", "iii"}
          obj.SumOfSquaresType = string ("three");
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
      obj.VarNames = value;
      obj.FactorNames = string (value);
    endfunction

    function obj = syncFactorSelectors_ (obj)
      if (obj.continuousSpecified_)
        continuous = obj.normalizeFactorSelector_ (obj.Continuous, ...
                                                    "CategoricalFactors", ...
                                                    false);
        obj.Continuous = continuous;
        obj.CategoricalFactors = setdiff (1:obj.nFactors_, continuous);
      else
        obj.CategoricalFactors = obj.normalizeFactorSelector_ (...
                                   obj.CategoricalFactors, ...
                                   "CategoricalFactors", true);
        obj.Continuous = setdiff (1:obj.nFactors_, ...
                                  obj.CategoricalFactors);
      endif
      obj.RandomFactors = obj.normalizeFactorSelector_ (obj.RandomFactors, ...
                                                        "RandomFactors", ...
                                                        true);
      obj.Random = obj.RandomFactors;
    endfunction

    function idx = normalizeFactorSelector_ (obj, value, option, allow_all)
      if (isempty (value))
        idx = [];
        return;
      endif
      if (obj.isName_ (value) && strcmpi (char (value), "all"))
        if (! allow_all)
          error ("anova: %s must identify valid factors.", option);
        endif
        idx = 1:obj.nFactors_;
        return;
      endif
      if (islogical (value))
        if (! isvector (value) || numel (value) != obj.nFactors_)
          error ("anova: %s logical input must have one value per factor.", ...
                 option);
        endif
        idx = find (value(:))';
        return;
      endif
      if (isstring (value))
        value = cellstr (value(:))';
      elseif (ischar (value))
        value = {value};
      endif
      if (iscellstr (value))
        [found, idx] = ismember (value, obj.VarNames);
        if (! all (found))
          error ("anova: %s contains an unknown factor name.", option);
        endif
        idx = idx(:)';
        return;
      endif
      if (! isnumeric (value) || ! isvector (value) ...
          || any (value != fix (value)) || any (value < 1) ...
          || any (value > obj.nFactors_))
        error ("anova: %s must contain valid factor indices.", option);
      endif
      idx = unique (value(:)', "stable");
    endfunction

    function obj = setFormula_ (obj)
      if (isempty (obj.VarNames))
        obj.VarNames = arrayfun (@(k) sprintf ("Factor%d", k), ...
                                 1:obj.nFactors_, "UniformOutput", false);
        obj.FactorNames = string (obj.VarNames);
      elseif (numel (obj.VarNames) != obj.nFactors_)
        error ("anova: FactorNames must contain one name per factor.");
      endif
      if (ischar (obj.ModelType) ...
          && ! isempty (strfind (obj.ModelType, "~")))
        schema = parseWilkinsonFormula (obj.ModelType, "matrix");
        if (! isscalar (schema.ResponseIdx))
          error ("anova: ModelSpecification must contain one response.");
        endif
        response = schema.VariableNames{schema.ResponseIdx};
        if (! strcmp (response, obj.ResponseName))
          error (strcat ("anova: formula response must match", ...
                         " ResponseName."));
        endif
        [schema_names, schema_model] = obj.schemaModel_ (schema);
        [found, columns_] = ismember (obj.VarNames, schema_names);
        if (! all (found))
          error ("anova: formula must use names in FactorNames.");
        endif
        model = schema_model(:, columns_);
        [model, obj.nesting_] = obj.nestedSpecification_ (...
                                  obj.ModelType, obj.VarNames, model);
        obj.ModelType = model;
        obj.formulaTerms_ = model;
        obj.formulaText_ = obj.ModelSpecification;
        return;
      endif
      if (obj.nFactors_ == 0)
        obj.formulaTerms_ = zeros (0, 0);
        obj.formulaText_ = sprintf ("%s ~ 1", obj.ResponseName);
      else
        if (ischar (obj.ModelType) ...
            && obj.isPolynomialModel_ (obj.ModelType))
          model = obj.polynomialTerms_ (obj.ModelType);
          obj.ModelType = model;
        else
          model = obj.modelTerms_ (obj.ModelType);
        endif
        obj.formulaTerms_ = model;
        terms = obj.formulaTermNames_ (model);
        rhs = strjoin (terms, " + ");
        obj.formulaText_ = sprintf ("%s ~ 1 + %s", obj.ResponseName, rhs);
      endif
    endfunction

    function model = modelTerms_ (obj, specification)
      if (isnumeric (specification) && ! isscalar (specification))
        model = specification;
        return;
      endif
      if (isnumeric (specification))
        max_order = min (specification, obj.nFactors_);
      elseif (strcmpi (specification, "linear"))
        max_order = 1;
      elseif (strcmpi (specification, "interaction"))
        max_order = min (2, obj.nFactors_);
      else
        max_order = obj.nFactors_;
      endif
      model = zeros (0, obj.nFactors_);
      for order = 1:max_order
        combinations = nchoosek (1:obj.nFactors_, order);
        block = zeros (rows (combinations), obj.nFactors_);
        for row = 1:rows (combinations)
          block(row, combinations(row, :)) = 1;
        endfor
        model = [model; block];
      endfor
    endfunction

    function names = formulaTermNames_ (obj, model)
      names = obj.termNames_ (model);
      if (isempty (obj.nesting_))
        return;
      endif
      for factor = find (any (obj.nesting_, 2))'
        parents = find (obj.nesting_(factor, :));
        nested = sprintf ("%s(%s)", obj.VarNames{factor}, ...
                          strjoin (obj.VarNames(parents), ","));
        for term = find (model(:, factor) > 0)'
          names{term} = strrep (names{term}, obj.VarNames{factor}, nested);
        endfor
      endfor
    endfunction

    function obj = updateFormula_ (obj)
      separator = strfind (obj.formulaText_, "~");
      if (isempty (separator))
        linear_predictor = obj.formulaText_;
      else
        linear_predictor = strtrim (obj.formulaText_(separator(1) + 1:end));
      endif
      term_names = obj.formulaTermNames_ (obj.formulaTerms_);
      nesting = obj.nesting_;
      if (isempty (nesting))
        nesting = false (obj.nFactors_);
      endif
      obj.Formula = struct (...
        "ResponseName", string (obj.ResponseName), ...
        "PredictorNames", string (obj.VarNames), ...
        "Terms", obj.formulaTerms_, ...
        "TermNames", string (term_names(:)), ...
        "HasIntercept", true, ...
        "Nesting", logical (nesting), ...
        "LinearPredictor", string (linear_predictor), ...
        "Text", string (obj.formulaText_));
    endfunction

    function tf = isPolynomialModel_ (obj, model)
      tf = any (strcmpi (model, {"purequadratic", "quadratic"})) ...
           || ! isempty (regexp (lower (model), '^poly[0-9]+$', 'once'));
    endfunction

    function model = polynomialTerms_ (obj, specification)
      n = obj.nFactors_;
      if (strcmpi (specification, "purequadratic"))
        model = [eye(n); 2 * eye(n)];
        return;
      elseif (strcmpi (specification, "quadratic"))
        interactions = zeros (0, n);
        if (n > 1)
          pairs = nchoosek (1:n, 2);
          interactions = zeros (rows (pairs), n);
          for i = 1:rows (pairs)
            interactions(i, pairs(i, :)) = 1;
          endfor
        endif
        model = [eye(n); interactions; 2 * eye(n)];
        return;
      endif

      degrees = specification(5:end) - '0';
      if (numel (degrees) != n)
        error (strcat ("anova: polyIJK must specify one degree per", ...
                       " factor."));
      endif
      model = [];
      for factor = 1:n
        for exponent = 1:degrees(factor)
          row = zeros (1, n);
          row(factor) = exponent;
          model(end + 1, :) = row;
        endfor
      endfor
      grids = cell (1, n);
      ranges = arrayfun (@(d) 0:d, degrees, "UniformOutput", false);
      [grids{:}] = ndgrid (ranges{:});
      combinations = cell2mat (cellfun (@(g) g(:), grids, ...
                                        "UniformOutput", false));
      interactions = combinations(sum (combinations > 0, 2) > 1 ...
                                  & sum (combinations, 2) <= max (degrees), :);
      model = [model; interactions];
    endfunction

    function names = termNames_ (obj, model)
      names = cell (rows (model), 1);
      for row = 1:rows (model)
        pieces = {};
        for column = find (model(row, :) != 0)
          exponent = model(row, column);
          if (exponent == 1)
            pieces{end + 1} = obj.VarNames{column};
          else
            pieces{end + 1} = sprintf ("%s^%d", ...
                                       obj.VarNames{column}, exponent);
          endif
        endfor
        names{row} = strjoin (pieces, ":");
      endfor
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

    function [stats_, args, sidak_alpha, dims] = comparisonArguments_ (...
                                                        obj, varargin)
      factors = {};
      option_names = {"alpha", "criticalvaluetype", "approximate", ...
                      "controlgroup", "display", "displayopt", ...
                      "estimate"};
      if (! isempty (varargin))
        first = varargin{1};
        is_option = obj.isName_ (first) ...
                    && any (strcmpi (char (first), option_names));
        if ((iscellstr (first) || isstring (first) || ischar (first)) ...
            && ! is_option)
          factors = first;
          varargin(1) = [];
        endif
      endif
      if (isempty (factors))
        if (obj.NumFactors != 1)
          error (strcat ("anova.multcompare: factors must be specified for", ...
                         " a multi-factor ANOVA."));
        endif
        dims = 1;
      else
        dims = obj.factorIndices_ (factors);
      endif
      if (mod (numel (varargin), 2) != 0)
        error ("anova.multcompare: name-value pairs must come in pairs.");
      endif

      alpha = 0.05;
      ctype = "tukey-kramer";
      has_approximate = false;
      has_control = false;
      args = {"display", "off"};
      for k = 1:2:numel (varargin)
        if (! obj.isName_ (varargin{k}))
          error ("anova.multcompare: parameter name must be text.");
        endif
        name = lower (char (varargin{k}));
        value = varargin{k + 1};
        switch (name)
          case "alpha"
            alpha = value;
          case "criticalvaluetype"
            if (! obj.isName_ (value))
              error ("anova.multcompare: CriticalValueType must be text.");
            endif
            ctype = lower (char (value));
          case "approximate"
            if (! ((islogical (value) || isnumeric (value)) ...
                   && isscalar (value) && any (value == [0, 1])))
              error (strcat ("anova.multcompare: Approximate must be a", ...
                             " logical scalar."));
            endif
            has_approximate = true;
          case "controlgroup"
            has_control = true;
            args = [args, {"controlgroup", value}];
          case {"display", "displayopt", "estimate"}
            args = [args, {name, value}];
          otherwise
            error ("anova.multcompare: parameter '%s' is not supported.", ...
                   varargin{k});
        endswitch
      endfor
      if (! (isnumeric (alpha) && isscalar (alpha) ...
             && alpha > 0 && alpha < 1))
        error ("anova.multcompare: Alpha must be a numeric scalar in (0, 1).");
      endif
      supported = {"tukey-kramer", "hsd", "dunn-sidak", ...
                   "bonferroni", "scheffe", "dunnett", "lsd"};
      if (! any (strcmp (ctype, supported)))
        error ("anova.multcompare: unsupported CriticalValueType '%s'.", ...
               ctype);
      endif
      if ((has_approximate || has_control) && ! strcmp (ctype, "dunnett"))
        error (strcat ("anova.multcompare: Approximate and ControlGroup", ...
                       " require CriticalValueType 'dunnett'."));
      endif
      if (strcmp (ctype, "dunnett") && numel (dims) != 1)
        error ("anova.multcompare: Dunnett's test requires one factor.");
      endif
      sidak_alpha = [];
      if (strcmp (ctype, "dunn-sidak"))
        sidak_alpha = alpha;
        ctype = "lsd";
      endif
      args = [args, {"alpha", alpha, "criticalvaluetype", ctype}];

      stats_ = obj.Stats;
      if (strcmp (obj.backend_, "anova2"))
        if (numel (dims) == 1)
          estimate = "column";
          if (dims == 2)
            estimate = "row";
          endif
          args = obj.replaceOption_ (args, "estimate", estimate);
        else
          [y, groups] = obj.anova2Data_ ();
          [~, ~, stats_] = anovan (y, groups, obj.buildAnovanArgs_ (){:});
          args = [args, {"dim", dims}];
        endif
      elseif (strcmp (obj.backend_, "anovan"))
        args = [args, {"dim", dims}];
      endif
    endfunction

    function args = replaceOption_ (obj, args, option, value)
      names = args(1:2:end);
      idx = find (cellfun (@(x) obj.isName_ (x) ...
                          && strcmpi (char (x), option), names), 1, "last");
      if (isempty (idx))
        args = [args, {option, value}];
      else
        args{2 * idx} = value;
      endif
    endfunction

    function idx = factorIndices_ (obj, factors)
      if (isstring (factors))
        factors = cellstr (factors(:))';
      elseif (ischar (factors))
        factors = {factors};
      endif
      if (! iscellstr (factors))
        error ("anova.multcompare: factors must be factor names.");
      endif
      [found, idx] = ismember (factors, obj.VarNames);
      if (! all (found) || isempty (idx))
        error ("anova.multcompare: factors contains an unknown factor name.");
      endif
      idx = idx(:)';
    endfunction

    function [y, groups] = anova2Data_ (obj)
      [nr, nc] = size (obj.Response);
      nlevels = nr / obj.reps_;
      y = obj.Response(:);
      column_factor = kron ((1:nc)', ones (nr, 1));
      row_factor = repmat (kron ((1:nlevels)', ones (obj.reps_, 1)), ...
                           nc, 1);
      groups = {column_factor, row_factor};
    endfunction

    function out = comparisonTable_ (obj, C, group_names, dims)
      if (istable (group_names))
        if (numel (dims) == 1)
          group1 = group_names{C(:, 1), 1};
          group2 = group_names{C(:, 2), 1};
          if (ischar (group1) && rows (group1) > 1)
            group1 = cellstr (group1);
            group2 = cellstr (group2);
          endif
        else
          group1 = group_names(C(:, 1), :);
          group2 = group_names(C(:, 2), :);
        endif
      else
        group1 = obj.groupIdentifiers_ (group_names(C(:, 1)), dims);
        group2 = obj.groupIdentifiers_ (group_names(C(:, 2)), dims);
      endif
      out = table (group1, group2, C(:, 4), C(:, 3), C(:, 5), C(:, 6), ...
                   "VariableNames", {"Group1", "Group2", ...
                   "MeanDifference", "MeanDifferenceLower", ...
                   "MeanDifferenceUpper", "pValue"});
    endfunction

    function groups = comparisonGroups_ (obj, stats_, dims)
      df = stats_.df;
      offsets = 1 + cumsum (df);
      terms = find (sum (stats_.terms(:,dims) > 0, 2) ...
                    == sum (stats_.terms > 0, 2));
      design = zeros (rows (stats_.X), columns (stats_.X));
      design(:,1) = 1;
      for term = terms'
        columns_ = offsets(term) - df(term) + 1:offsets(term);
        design(:,columns_) = stats_.X(:,columns_);
      endfor
      unique_design = unique (design, "rows", "stable");
      rows_ = zeros (rows (unique_design), 1);
      for i = 1:rows (unique_design)
        rows_(i) = find (all (design == unique_design(i,:), 2), 1);
      endfor
      selected = obj.Factors(rows_, obj.VarNames(dims));
      groups = selected;
    endfunction

    function groups = groupIdentifiers_ (obj, labels, dims)
      labels = cellstr (labels);
      if (isempty (labels))
        groups = obj.convertGroupValues_ (cell (0, 1), dims(1));
        return;
      endif
      parts = cellfun (@(x) strsplit (x, ", "), labels, ...
                       "UniformOutput", false);
      widths = cellfun (@numel, parts);
      has_names = ! any (widths != widths(1)) ...
                  && all (cellfun (@(x) ! isempty (strfind (x, "=")), ...
                                   [parts{:}]));
      if (! has_names)
        parts = cellfun (@(x) {x}, labels, "UniformOutput", false);
        widths(:) = 1;
      endif
      values = cell (numel (labels), widths(1));
      for i = 1:numel (labels)
        for j = 1:widths(1)
          if (has_names)
            separator = strfind (parts{i}{j}, "=");
            value = parts{i}{j}(separator(1) + 1:end);
          else
            value = parts{i}{j};
          endif
          values{i, j} = strtrim (value);
        endfor
      endfor
      if (widths(1) == 1)
        groups = obj.convertGroupValues_ (values(:, 1), dims(1));
      else
        columns_ = cell (1, widths(1));
        for j = 1:widths(1)
          columns_{j} = obj.convertGroupValues_ (values(:, j), dims(j));
        endfor
        groups = table (columns_{:}, ...
                        "VariableNames", obj.VarNames(dims));
      endif
    endfunction

    function values = convertGroupValues_ (obj, text_values, dim)
      original = obj.Factors{:, obj.VarNames{dim}};
      if (islogical (original))
        values = strcmpi (text_values, "true") ...
                 | strcmp (text_values, "1");
      elseif (isnumeric (original))
        values = cellfun (@str2double, text_values);
      elseif (iscategorical (original))
        values = categorical (text_values, categories (original));
      elseif (isstring (original))
        values = string (text_values);
      else
        values = text_values;
      endif
    endfunction

    function [G, names] = selectedFactorMatrix_ (obj, factors)
      if (ischar (factors))
        idx = find (strcmp (obj.VarNames, factors));
      elseif (iscellstr (factors))
        idx = cellfun (@(s) find (strcmp (obj.VarNames, s), 1), factors);
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
      if (isempty (obj.GROUP) && ! isvector (obj.Response))
        [n, m] = size (obj.Response);
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
      names = obj.VarNames(idx);
    endfunction

    function name = factorName_ (obj, idx)
      if (idx >= 1 && idx <= numel (obj.VarNames))
        name = obj.VarNames{idx};
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
    ## vector.  Pearson residuals scale the raw residuals by the RMSE.
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
      [atab, stats_] = obj.componentFit_ (sstype);
      sources = obj.publicSourceNames_ (atab(2:end-1, 1));
      nterms = numel (sources) - 1;
      term_names = sources(1:nterms);
      is_random = false (nterms, 1);
      q_coeff = ones (nterms, 1);
      variance_coeff = zeros (nterms, nterms);

      if (strcmp (obj.backend_, "anovan") && isfield (stats_, "terms"))
        terms = stats_.terms;
        if (! isempty (obj.RandomFactors))
          is_random = obj.randomTermMask_ (terms, stats_);
        endif
        [q_coeff, variance_coeff] = obj.emsCoefficients_ (stats_, terms, ...
                                                          sstype, is_random);
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
          ## A fixed term that contains this one contributes a Q component of
          ## its own, so the expected mean square names it too.
          for j = find (! is_random(:))'
            if (j != i && abs (variance_coeff(i, j)) > 1e-10)
              pieces{end + 1} = obj.emsTerm_ (variance_coeff(i, j), "Q", ...
                                              term_names{j});
            endif
          endfor
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

      [term_ms, term_df, error_ms, error_df] = obj.meanSquareData_ (atab);
      error_ms = obj.MSE;
      error_df = obj.DFE;
      [denominator_ms, denominator_df, denominator_formula] = ...
        obj.denominatorData_ (variance_coeff, is_random, term_ms, ...
                              term_df, error_ms, error_df, term_names);
      denominator_ms(end + 1, 1) = NaN;
      denominator_df(end + 1, 1) = NaN;
      denominator_formula{end + 1, 1} = "";
      row_names = [term_names; {"Error"}];
      ems = table (string (types), string (formulas), denominator_ms, ...
                   denominator_df, string (denominator_formula), ...
                   "VariableNames", {"Type", "ExpectedMeanSquares", ...
                   "MeanSquaresDenominator", "DFDenominator", ...
                   "FDenominator"}, "RowNames", row_names);
    endfunction

    function names = publicSourceNames_ (obj, names)
      if (strcmp (obj.backend_, "anova1") && ! isempty (names))
        idx = find (strcmpi (names, "Groups"), 1);
        if (! isempty (idx))
          names{idx} = obj.VarNames{1};
        endif
      elseif (strcmp (obj.backend_, "anova2"))
        idx = find (strcmpi (names, "Columns"), 1);
        if (! isempty (idx))
          names{idx} = obj.VarNames{1};
        endif
        idx = find (strcmpi (names, "Rows"), 1);
        if (! isempty (idx))
          names{idx} = obj.VarNames{2};
        endif
        idx = find (strcmpi (names, "Interaction"), 1);
        if (! isempty (idx))
          names{idx} = sprintf ("%s:%s", obj.VarNames{:});
        endif
      endif
    endfunction

    function [q_coeff, variance_coeff] = emsCoefficients_ (obj, stats_, ...
                                            terms, sstype, is_random)
      X = full (stats_.X);
      dfs = stats_.df(:);
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
        if (isfield (stats_, "vnested") && ! isempty (stats_.vnested))
          nested = factors(any (stats_.vnested(factors,:), 2));
          if (! isempty (nested))
            factors = union (factors, find (any (stats_.vnested(nested,:), 1)));
          endif
        endif
        if (all (! ismember (factors, obj.Continuous)))
          [~, ~, ids] = unique (stats_.grps(:, factors), "rows");
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
        ## Every term that contains this one contributes to its expected mean
        ## square, whether that term is fixed or random, so the coefficient is
        ## taken for all of them; which letter names it is decided later.
        for j = 1:nterms
          variance_coeff(i, j) = obj.projectionDifference_ (...
                                  included, excluded, bases{j}) / dfs(i);
        endfor
      endfor
    endfunction

    function [term_ms, term_df, error_ms, error_df] = ...
              meanSquareData_ (obj, atab)
      source_col = obj.findAtabColumn_ (atab, {"Source"});
      ms_col = obj.findAtabColumn_ (atab, ...
                                    {"MeanSquares", "Mean Sq.", ...
                                     "Mean Sq", "MS"});
      df_col = obj.findAtabColumn_ (atab, {"DF", "d.f.", "df"});
      sources = atab(2:end, source_col);
      error_row = find (strcmpi (sources, "Error"), 1);
      total_row = find (strcmpi (sources, "Total"), 1);
      term_rows = setdiff (1:numel (sources), [error_row, total_row], ...
                           "stable");
      term_ms = cell2mat (atab(term_rows + 1, ms_col));
      term_df = cell2mat (atab(term_rows + 1, df_col));
      error_ms = atab{error_row + 1, ms_col};
      error_df = atab{error_row + 1, df_col};
    endfunction

    function [denominator_ms, denominator_df, formulas] = ...
              denominatorData_ (obj, variance_coeff, is_random, term_ms, ...
                                term_df, error_ms, error_df, term_names)
      nterms = numel (term_ms);
      random_terms = find (is_random(:));
      if (isempty (random_terms))
        denominator_ms = repmat (error_ms, nterms, 1);
        denominator_df = repmat (error_df, nterms, 1);
        formulas = repmat ({"MS(Error)"}, nterms, 1);
        return;
      endif
      random_ms = [term_ms(random_terms); error_ms];
      random_df = [term_df(random_terms); error_df];
      coefficients = [variance_coeff(random_terms, random_terms), ...
                      ones(numel (random_terms), 1);
                      zeros(1, numel (random_terms)), 1];
      names = [term_names(random_terms); {"Error"}];
      denominator_ms = zeros (nterms, 1);
      denominator_df = zeros (nterms, 1);
      formulas = cell (nterms, 1);
      for i = 1:nterms
        target = [variance_coeff(i, random_terms), 1];
        own = find (random_terms == i, 1);
        if (! isempty (own))
          target(own) = 0;
        endif
        weights = pinv (coefficients') * target';
        ## pinv leaves rounding dust on the mean squares a denominator does
        ## not use.  Drop it, so a denominator that is one mean square is
        ## exactly that mean square, and so the formula and the number agree.
        weights(abs (weights) <= 1e-10) = 0;
        denominator_ms(i) = weights' * random_ms;
        ## Satterthwaite's approximation, over the mean squares the
        ## denominator actually uses.  A mean square carrying no degrees of
        ## freedom has no sampling distribution to combine, so a denominator
        ## resting on one has no degrees of freedom either.
        used = (weights != 0);
        if (any (used & (random_df(:) <= 0)))
          denominator_df(i) = 0;
        else
          parts = (weights(used) .* random_ms(used)) .^ 2 ./ random_df(used);
          denominator_df(i) = denominator_ms(i) ^ 2 / sum (parts);
        endif
        formulas{i} = obj.meanSquareFormula_ (weights, names);
      endfor
    endfunction

    function formula = meanSquareFormula_ (obj, weights, names)
      pieces = {};
      for i = 1:numel (weights)
        coefficient = weights(i);
        if (abs (coefficient) <= 1e-10)
          continue;
        elseif (abs (coefficient - 1) <= 1e-10)
          term = sprintf ("MS(%s)", names{i});
        elseif (abs (coefficient + 1) <= 1e-10)
          term = sprintf ("-MS(%s)", names{i});
        else
          term = sprintf ("%.6g*MS(%s)", coefficient, names{i});
        endif
        if (! isempty (pieces) && coefficient > 0)
          term = ["+", term];
        endif
        pieces{end + 1} = term;
      endfor
      formula = strjoin (pieces, "");
    endfunction

    function [coefficients, mean_squares, dfs, names] = ...
              varianceSystem_ (obj, stats_, atab, sstype)
      [term_ms, term_df, error_ms, error_df] = obj.meanSquareData_ (atab);
      sources = obj.publicSourceNames_ (atab(2:end-1, 1));
      term_names = sources(1:numel (term_ms));
      is_random = false (numel (term_ms), 1);
      variance_coeff = zeros (numel (term_ms));
      if (! isempty (obj.RandomFactors))
        if (strcmp (obj.backend_, "anovan") && isfield (stats_, "terms"))
          terms = stats_.terms;
          is_random = obj.randomTermMask_ (terms, stats_);
          [~, variance_coeff] = obj.emsCoefficients_ (stats_, terms, ...
                                                       sstype, is_random);
        elseif (strcmp (obj.backend_, "anova1"))
          is_random(1) = true;
          [~, variance_coeff(1, 1)] = obj.oneWayEmsCoefficient_ ();
        endif
      endif
      random_terms = find (is_random);
      coefficients = [variance_coeff(random_terms, random_terms), ...
                      ones(numel (random_terms), 1);
                      zeros(1, numel (random_terms)), 1];
      mean_squares = [term_ms(random_terms); error_ms];
      dfs = [term_df(random_terms); error_df];
      names = [term_names(random_terms); {"Error"}];
    endfunction

    function atab = applyRandomInference_ (obj, atab, stats_, sstype)
      if (isempty (obj.RandomFactors) || ! isfield (stats_, "terms"))
        return;
      endif
      terms = stats_.terms;
      is_random = obj.randomTermMask_ (terms, stats_);
      [~, variance_coeff] = obj.emsCoefficients_ (stats_, terms, ...
                                                   sstype, is_random);
      [term_ms, term_df, error_ms, error_df] = obj.meanSquareData_ (atab);
      names = obj.publicSourceNames_ (atab(2:rows (terms) + 1, 1));
      [denominator_ms, denominator_df] = obj.denominatorData_ (...
          variance_coeff, is_random, term_ms, term_df, error_ms, ...
          error_df, names);
      f_col = obj.findAtabColumn_ (atab, {"F"});
      p_col = obj.findAtabColumn_ (atab, {"pValue", "Prob>F"});
      f_stat = term_ms ./ denominator_ms;
      p_value = fcdf (f_stat, term_df, denominator_df, "upper");
      atab(2:numel (term_ms) + 1, f_col) = num2cell (f_stat);
      atab(2:numel (term_ms) + 1, p_col) = num2cell (p_value);
    endfunction

    function is_random = randomTermMask_ (obj, terms, stats_)
      membership = terms > 0;
      if (isfield (stats_, "vnested") && ! isempty (stats_.vnested))
        for factor = 1:columns (membership)
          parents = find (stats_.vnested(factor,:));
          if (! isempty (parents))
            membership(:, parents) |= membership(:, factor);
          endif
        endfor
      endif
      is_random = any (membership(:, obj.RandomFactors), 2);
    endfunction

    function [included, excluded] = emsModels_ (obj, X, blocks, terms, ...
                                                term, sstype)
      nterms = numel (blocks);
      switch (sstype)
        case 1
          included_terms = 1:term;
          excluded_terms = 1:term-1;
        case 2
          factors = find (terms(term, :));
          excluded_terms = find (any (terms(:, factors) ...
                                      != terms(term, factors), 2))';
          included_terms = [term, excluded_terms];
        case "h"
          factors = find (terms(term, :));
          excluded_terms = find (any (terms(:, factors) ...
                                      < terms(term, factors), 2))';
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
      if (isvector (obj.Response))
        group = obj.GROUP(:);
      else
        [n, groups] = size (obj.Response);
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


    function validateSpec_ (obj)
      if (! ((isnumeric (obj.SSType) && isscalar (obj.SSType) ...
              && any (obj.SSType == [1, 2, 3])) ...
             || (ischar (obj.SSType) && strcmp (obj.SSType, "h"))))
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
                         " 'interactions', 'purequadratic', 'quadratic',", ...
                         " 'polyIJK', 'full', or a terms matrix."));
        endif
      elseif (! isnumeric (obj.ModelType))
        error (strcat ("anova: ModelSpecification must be a string", ...
                       " or a numeric terms matrix."));
      endif
      if (isnumeric (obj.ModelType) && isscalar (obj.ModelType) ...
          && (obj.ModelType != fix (obj.ModelType) || obj.ModelType < 1))
        error (strcat ("anova: integer ModelSpecification must be a", ...
                       " positive integer."));
      endif
      if (isnumeric (obj.ModelType) && ! isscalar (obj.ModelType) ...
          && (any (! isfinite (obj.ModelType(:))) ...
              || any (obj.ModelType(:) < 0) ...
              || any (obj.ModelType(:) != fix (obj.ModelType(:)))))
        error (strcat ("anova: terms matrix entries must be nonnegative", ...
                       " integers."));
      endif
      if (! iscellstr (obj.VarNames))
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
          && ! isscalar (obj.ModelType) ...
          && columns (obj.ModelType) != obj.nFactors_)
        error ("anova: terms matrix must have one column per factor.");
      endif
      if (isnumeric (obj.ModelType) && ! isscalar (obj.ModelType))
        powered = find (any (obj.ModelType > 1, 1));
        if (any (ismember (powered, obj.CategoricalFactors)))
          error (strcat ("anova: polynomial powers require continuous", ...
                         " factors."));
        endif
      endif
    endfunction

    function validateGroupLength_ (obj, group, nobs)
      if (iscell (group) ...
          && all (cellfun (@(c) isvector (c) || ischar (c), group(:))) ...
          && size (group, 1) == 1)
        for k = 1:numel (group)
          if (ischar (group{k}) && rows (group{k}) > 1)
            count = rows (group{k});
          else
            count = numel (group{k});
          endif
          if (count != nobs)
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
        if (ismatrix (obj.Response) && ! isvector (obj.Response))
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
      anova2_model = ischar (obj.ModelType) ...
                     && any (strcmpi (obj.ModelType, ...
                                      {"linear", "interaction", "full"}));
      ## anova2 requires a balanced design and rejects missing observations,
      ## so data holding any NaN is fitted through anovan, which omits them.
      if (! isempty (obj.reps_) && ismatrix (obj.Response) ...
          && ! isvector (obj.Response) ...
          && ! any (isnan (obj.Response(:))) ...
          && isempty (obj.Continuous) && isempty (obj.Weights) ...
          && isempty (obj.RandomFactors) ...
          && (isempty (obj.nesting_) || ! any (obj.nesting_(:))) ...
          && anova2_model)
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
      endswitch
      obj = obj.updatePublicResults_ ();
      obj.fitted_ = true;
      obj.dirty_  = false;
    endfunction

    function obj = fitAnova1_ (obj)
      ## Run the backend silently; this class owns display and plotting.
      if (isvector (obj.Response))
        [~, atab, stats] = anova1 (obj.Response, obj.GROUP, 'off');
      else
        [~, atab, stats] = anova1 (obj.Response, [], 'off');
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
      [~, atab, stats] = anova2 (obj.Response, obj.reps_, 'off', ...
                                 modelarg);
      obj.AnovaTable = atab;
      obj.Stats      = stats;
      obj.DFE = stats.df;
      obj.MSE = stats.sigmasq;
      ## Balanced anova2 does not expose coefficients or residuals.
    endfunction

    function obj = fitAnovan_ (obj)
      [y_vec, group_arg] = obj.anovanData_ ();
      [~, atab, stats] = anovan (y_vec, group_arg, ...
                                 obj.buildAnovanArgs_(){:});
      obj.AnovaTable = obj.applyRandomInference_ (atab, stats, obj.SSType);
      obj.Stats      = stats;
      if (isfield (stats, 'coeffs'))
        obj.coefficientStats_ = stats.coeffs;
        obj.Coefficients = stats.coeffs(:, 1);
      endif
      if (isfield (stats, 'resid'))
        obj.rawResiduals_ = stats.resid;
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
      ## anovan omits any row holding a missing value, so adopt the
      ## observations it actually fitted whenever it dropped some.
      if (isfield (stats, 'Y') && numel (stats.Y) != numel (obj.Y))
        obj.Y = stats.Y(:);
        obj.Response = obj.Y;
        obj.NumObservations = numel (obj.Y);
        if (isfield (stats, 'grps') ...
            && columns (stats.grps) == numel (obj.VarNames) ...
            && rows (stats.grps) == obj.NumObservations)
          obj.GROUP = num2cell (stats.grps, 1);
          obj.Factors = table (obj.GROUP{:}, 'VariableNames', obj.VarNames);
        endif
      endif
      if (! isempty (obj.DesignMatrix) && ! isempty (obj.Coefficients))
        obj.FittedValues = full (obj.DesignMatrix) * obj.Coefficients;
        obj.rawResiduals_ = obj.Y - obj.FittedValues;
      endif
      if (! isempty (obj.rawResiduals_))
        obj.Residuals = obj.residualTable_ (obj.rawResiduals_, obj.MSE);
      endif
    endfunction

    function atab = componentStats_ (obj, sstype)
      [atab] = obj.componentFit_ (sstype);
    endfunction

    function [atab, stats_] = componentFit_ (obj, sstype)
      if (strcmp (obj.backend_, "anova2"))
        obj.ensureFit_ ();
        atab = obj.AnovaTable;
        stats_ = obj.Stats;
        return;
      endif
      [y_vec, group_arg] = obj.anovanData_ ();
      args = obj.buildAnovanArgs_ (sstype);
      [~, atab, stats_] = anovan (y_vec, group_arg, args{:});
      atab = obj.applyRandomInference_ (atab, stats_, sstype);
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
      if (! isempty (obj.reps_) && isempty (obj.GROUP) ...
          && ! isvector (obj.Response))
        [y_vec, group_arg] = obj.anova2Data_ ();
      elseif (isempty (obj.GROUP) && ! isvector (obj.Response))
        [n, m] = size (obj.Response);
        y_vec = obj.Response(:);
        group = reshape (repmat (1:m, n, 1), [], 1);
        group_arg = {group};
      else
        y_vec = obj.Response(:);
        group_arg = obj.GROUP;
        if (isempty (group_arg))
          group_arg = {};
        endif
      endif
    endfunction

    function s = sstypeLabel_ (obj)
      switch (obj.SSType)
        case 1; s = "I";
        case 2; s = "II";
        case 3; s = "III";
        otherwise; s = "hierarchical";
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
      if (! isempty (obj.nesting_) && any (obj.nesting_(:)))
        nv = [nv, {'nested', obj.nesting_}];
      endif
      ## Random terms are retained so their expected mean squares can provide
      ## the correct F denominators and variance component estimates.
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

    function residuals = residualTable_ (obj, raw, mse)
      pearson = raw ./ sqrt (mse);
      residuals = table (raw(:), pearson(:), "VariableNames", ...
                         {"Raw", "Pearson"});
    endfunction

    function obj = updatePublicResults_ (obj)
      if (strcmp (obj.backend_, "anovan") ...
          && ! isempty (obj.coefficientStats_))
        [obj.Coefficients, names] = obj.expandedCoefficients_ ();
        obj.ExpandedFactorNames = string (names);
      elseif (strcmp (obj.backend_, "anova1"))
        obj = obj.populateOneWayResults_ ();
      elseif (isfield (obj.Stats, "coeffnames") ...
              && numel (obj.Stats.coeffnames) == numel (obj.Coefficients))
        obj.ExpandedFactorNames = string (obj.Stats.coeffnames(:));
      elseif (! isempty (obj.Coefficients))
        names = arrayfun (@(k) sprintf ("Coefficient%d", k), ...
                          (1:numel (obj.Coefficients))', ...
                          "UniformOutput", false);
        obj.ExpandedFactorNames = string (names);
      endif
      if (isempty (obj.AnovaTable) || isempty (obj.MSE))
        return;
      endif
      source = obj.findAtabColumn_ (obj.AnovaTable, {"Source"});
      ss = obj.findAtabColumn_ (obj.AnovaTable, ...
                                {"SumOfSquares", "Sum Sq.", "Sum Sq", ...
                                 "SS"});
      df = obj.findAtabColumn_ (obj.AnovaTable, {"DF", "d.f.", "df"});
      names = obj.AnovaTable(2:end, source);
      error_row = find (strcmpi (names, "Error"), 1) + 1;
      total_row = find (strcmpi (names, "Total"), 1) + 1;
      if (isempty (error_row) || isempty (total_row))
        return;
      endif
      sse = obj.AnovaTable{error_row, ss};
      sst = obj.AnovaTable{total_row, ss};
      total_df = obj.AnovaTable{total_row, df};
      if (isnumeric (total_df) && isscalar (total_df) && isfinite (total_df))
        obj.NumObservations = total_df + 1;
      endif
      ssr = sst - sse;
      if (sst == 0)
        rsquared = NaN;
        adjusted = NaN;
      else
        rsquared = ssr / sst;
        adjusted = 1 - (obj.NumObservations - 1) * sse ...
                       / (obj.DFE * sst);
      endif
      obj.Metrics = table (obj.MSE, sqrt (obj.MSE), sse, ssr, sst, ...
                           rsquared, adjusted, "VariableNames", ...
                           {"MSE", "RMSE", "SSE", "SSR", "SST", ...
                            "RSquared", "AdjustedRSquared"});
    endfunction

    function [coefficients, names] = expandedCoefficients_ (obj)
      coefficients = obj.coefficientStats_(1, 1);
      names = {"(Intercept)"};
      first = 2;
      for term = 1:rows (obj.Stats.terms)
        factors = find (obj.Stats.terms(term, :));
        mapping = 1;
        expanded_names = {""};
        for factor = factors
          if (any (factor == obj.Continuous))
            contrast = 1;
            factor_names = {obj.VarNames{factor}};
          elseif (isfield (obj.Stats, "vnested") ...
                  && any (obj.Stats.vnested(factor,:)))
            [contrast, factor_names] = ...
              obj.nestedCoefficientMap_ (factor);
          else
            contrast = obj.Stats.contrasts{factor};
            levels = obj.Stats.grpnames{factor};
            factor_names = cellfun (@(x) sprintf ("(%s==%s)", ...
                                      obj.VarNames{factor}, ...
                                      obj.levelText_ (x)), levels, ...
                                     "UniformOutput", false);
          endif
          mapping = kron (mapping, contrast);
          combined = cell (numel (expanded_names) * numel (factor_names), 1);
          cursor = 1;
          for left = 1:numel (expanded_names)
            for right = 1:numel (factor_names)
              if (isempty (expanded_names{left}))
                combined{cursor} = factor_names{right};
              else
                combined{cursor} = sprintf ("%s:%s", ...
                                             expanded_names{left}, ...
                                             factor_names{right});
              endif
              cursor += 1;
            endfor
          endfor
          expanded_names = combined;
        endfor
        width = obj.Stats.df(term);
        compact = obj.coefficientStats_(first:first + width - 1, 1);
        coefficients = [coefficients; mapping * compact];
        names = [names; expanded_names];
        first += width;
      endfor
    endfunction

    function [mapping, names] = nestedCoefficientMap_ (obj, factor)
      parents = find (obj.Stats.vnested(factor,:));
      [parent_codes, ~, parent_id] = unique (...
        obj.Stats.grps(:, parents), "rows", "stable");
      mapping = zeros (0, 0);
      names = {};
      for parent = 1:rows (parent_codes)
        rows_ = parent_id == parent;
        child_codes = unique (obj.Stats.grps(rows_, factor), "stable");
        nlevels = numel (child_codes);
        contrast = [zeros(1, nlevels - 1); eye(nlevels - 1)] ...
                   - 1 / nlevels;
        mapping = blkdiag (mapping, contrast);
        for child = child_codes(:)'
          pieces = cell (1, numel (parents) + 1);
          pieces{1} = sprintf ("(%s==%s)", obj.VarNames{factor}, ...
                               obj.levelText_ (...
                                 obj.Stats.grpnames{factor}{child}));
          for k = 1:numel (parents)
            code = parent_codes(parent, k);
            pieces{k + 1} = sprintf ("(%s==%s)", ...
                                     obj.VarNames{parents(k)}, ...
                                     obj.levelText_ (...
                                       obj.Stats.grpnames{parents(k)}{code}));
          endfor
          names{end + 1, 1} = strjoin (pieces, ":");
        endfor
      endfor
    endfunction

    function obj = populateOneWayResults_ (obj)
      means = obj.Stats.means(:);
      intercept = mean (means);
      obj.Coefficients = [intercept; means - intercept];
      level_names = cellfun (@(x) sprintf ("(%s==%s)", ...
                              obj.VarNames{1}, obj.levelText_ (x)), ...
                             obj.Stats.gnames(:), "UniformOutput", false);
      obj.ExpandedFactorNames = string ([{"(Intercept)"}; level_names]);
      if (isvector (obj.Response))
        group_id = grp2idx (obj.GROUP);
        fitted = NaN (size (group_id));
        grouped = isfinite (group_id) & group_id > 0;
        fitted(grouped) = means(group_id(grouped));
      else
        fitted = repmat (means', rows (obj.Response), 1);
        fitted = fitted(:);
      endif
      valid = isfinite (obj.Y) & isfinite (fitted);
      obj.FittedValues = NaN (size (obj.Y));
      obj.FittedValues(valid) = fitted(valid);
      obj.rawResiduals_ = obj.Y - obj.FittedValues;
      obj.Residuals = obj.residualTable_ (obj.rawResiduals_, obj.MSE);
    endfunction

    function text = levelText_ (obj, value)
      if (isnumeric (value) || islogical (value))
        text = num2str (value);
      elseif (isstring (value))
        text = char (value);
      else
        text = value;
      endif
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
%! ## Fit an ANOVA object and inspect component and summary statistics
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! aov = anova (g, y, 'FactorNames', {'Treatment'});
%! component = stats (aov)
%! summary_table = stats (aov, 'summary')

%!demo
%! ## Estimate group means and perform post-hoc comparisons
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! aov = anova (g, y, 'SumOfSquaresType', 'two');
%! means = groupmeans (aov)
%! comparisons = multcompare (aov, 'display', 'off')

%!demo
%! ## Plot multiple-comparison intervals
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! aov = anova (g, y, 'SumOfSquaresType', 'two');
%! plotComparisons (aov);

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
%! assert_equal (class (a), "anova");
%! assert_equal (a.Y, y(:));
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
%! assert_equal (char (a.SumOfSquaresType), 'three');
%! assert_equal (a.ResponseName, 'Y');
%! assert_equal (cellstr (a.FactorNames), {'Factor1'});
%! assert_equal (a.RandomFactors, []);
%! assert_equal (a.CategoricalFactors, 1);
%! assert_equal (a.NumFactors, 1);
%! assert_equal (a.NumObservations, 4);

## Eager one-way fits populate the derived public result properties.
%!test
%! a = anova ([1;1;2;2], [1;2;3;4]);
%! assert_equal (a.Coefficients, [2.5; -1; 1], 1e-12);
%! assert_equal (cellstr (a.ExpandedFactorNames), ...
%!               {'(Intercept)'; '(Factor1==1)'; '(Factor1==2)'});
%! assert_equal (a.FittedValues, [1.5; 1.5; 3.5; 3.5], 1e-12);
%! assert_equal (a.Residuals.Raw, [-0.5; 0.5; -0.5; 0.5], 1e-12);
%! assert_equal (a.Residuals.Pearson, ...
%!               [-1; 1; -1; 1] / sqrt (2), 1e-12);
%! assert_equal (a.DesignMatrix, []);

## Name-value parsing: MATLAB-compatible names
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = repmat ([1;1;2;2], 3, 1);
%! a = anova ({g1, g2}, y, 'SumOfSquaresType', 'two', ...
%!            'ModelSpecification', 'full', 'FactorNames', {'A', 'B'});
%! assert_equal (char (a.SumOfSquaresType), 'two');
%! assert_equal (a.ModelSpecification, 'full');
%! assert_equal (cellstr (a.FactorNames), {'A', 'B'});

## Name-value parsing: case-insensitive names, displayopt alias
%!test
%! a = anova ([1;1;2;2], [1;2;3;4], 'SumOfSquaresType', 'one', 'displayopt', 'on');
%! assert_equal (char (a.SumOfSquaresType), 'one');

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
%! assert_equal (a.NumObservations, 7);
%! assert_equal (a.Stats.source, 'anovan');

## Backend selection: weights force anovan
%!test
%! a = anova ([1;1;2;2;3;3], [1;2;3;4;5;6], 'Weights', ones (6, 1));
%! assert_equal (a.Stats.source, 'anovan');
%! assert_equal (stats (a).Properties.RowNames, ...
%!               {'Factor1'; 'Error'; 'Total'});

## --- Week 2: fit delegation smoke tests --------------------------------

## fit_(): one-way fixture populates the unified result surface
%!test
%! y = [1; 2; 3; 4; 5; 6];
%! g = [1; 1; 2; 2; 3; 3];
%! a = anova (g, y);
%! a.fit ();
%! assert_equal (size (a.AnovaTable), [4, 6]);
%! assert_equal (a.Stats.source, 'anova1');
%! assert_equal (a.DFE, 3);
%! assert_equal (a.MSE, 0.5, 1e-12);

## fit_(): two-way balanced fixture (anova2 backend, popcorn data)
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (popcorn, [], 'reps', 3);
%! assert_equal (! isempty (strfind (evalc ('disp (a)'), '2-way anova')), true);
%! a.fit ();
%! assert_equal (size (a.AnovaTable), [5, 6]);
%! assert_equal (a.MSE, 0.125, 1e-12);
%! assert_equal (a.Stats.sigmasq, 0.125, 1e-12);

## fit_(): N-way fixture (anovan backend, three factors)
%!test
%! y = (1:24)';
%! g1 = repmat ([1;2], 12, 1);
%! g2 = repmat ([1;1;2;2], 6, 1);
%! g3 = repmat ([1;1;1;1;2;2;2;2], 3, 1);
%! a = anova ({g1, g2, g3}, y);
%! assert_equal (a.NumFactors, 3);
%! a.fit ();
%! assert_equal (a.Stats.source, 'anovan');
%! assert_equal (size (a.AnovaTable), [6, 9]);
%! assert_equal (size (a.Coefficients), [7, 1]);
%! assert_equal (size (a.Residuals), [24, 2]);
%! assert_equal (size (a.DesignMatrix), [24, 4]);

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
%! assert_equal (char (a.SumOfSquaresType), 'two');
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
%! assert_equal (T{2, 5}, 63, 1e-12);
%! assert_equal (T{3, 5}, 1, 1e-12);
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
%! assert_equal (a.Residuals.Raw, ...
%!               [-1; 0; 1; -1; 0; 1; -1; 0; 1], 1e-12);

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
%! C = multcompare (a, 'CriticalValueType', 'bonferroni');
%! assert_equal (istable (C), true);
%! assert_equal (C.Group1, [1; 1; 2]);
%! assert_equal (C.Group2, [2; 3; 3]);
%! assert_equal (C.MeanDifference, [-3; -9; -6], 1e-12);
%! assert_equal (all (C.pValue >= 0 & C.pValue <= 1), true);

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

## multcompare(): anovan backend returns MATLAB's comparison table
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! C = multcompare (a);
%! assert_equal (size (C), [3, 6]);
%! assert_equal (C.Properties.VariableNames, ...
%!               {'Group1', 'Group2', 'MeanDifference', ...
%!                'MeanDifferenceLower', 'MeanDifferenceUpper', 'pValue'});

## multcompare(): two-way balanced via anova2 backend
%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (popcorn, [], 'reps', 3);
%! C = multcompare (a, 'Factor1');
%! assert_equal (size (C), [3, 6]);

## multcompare(): runs after a fresh construction (triggers ensureFit_)
%!test
%! a = anova ([1;1;2;2;3;3], (1:6)', 'SumOfSquaresType', 'two');
%! C = multcompare (a);
%! assert_equal (size (C), [3, 6]);

## MATLAB-compatible public methods and properties
%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y, "FactorNames", {"Brand"}, ...
%!            "ResponseName", "Yield", "SumOfSquaresType", "two");
%! assert_equal (char (a.Formula.Text), "Yield ~ 1 + Brand");
%! assert_equal (cellstr (a.Formula.PredictorNames), {"Brand"});
%! assert_equal (a.Formula.Terms, 1);
%! assert_equal (istable (a.Factors), true);
%! assert_equal (a.Factors.Brand, g);
%! assert_equal (a.Response, y);
%! assert_equal (isstring (a.FactorNames), true);
%! assert_equal (cellstr (a.FactorNames), {"Brand"});
%! assert_equal (isstring (a.SumOfSquaresType), true);
%! assert_equal (char (a.SumOfSquaresType), "two");
%! assert_equal (isstruct (a.Formula), true);
%! assert_equal (! any (strcmp (methods ("anova"), "predict")), true);
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
%! assert_equal (isstring (a.ExpandedFactorNames), true);
%! assert_equal (cellstr (a.ExpandedFactorNames), ...
%!               {'(Intercept)'; '(Factor1==1)'; '(Factor1==2)'; ...
%!                '(Factor1==3)'});

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
%!   assert_equal (all (ishghandle (h)), true);
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

## boxchart(): returns graphics handles in the requested axes
%!test
%! hfig = figure ("visible", "off");
%! unwind_protect
%!   ax = axes ("parent", hfig);
%!   a = anova ([1;1;1;2;2;2], [1;2;3;4;5;6], ...
%!              "SumOfSquaresType", "two");
%!   h = boxchart (a, ax);
%!   assert_equal (all (ishghandle (h)), true);
%!   assert_equal (all (arrayfun (@(x) ancestor (x, "axes") == ax, h)), true);
%! unwind_protect_cleanup
%!   close (hfig);
%! end_unwind_protect

## plotComparisons(): plots adjusted means in the requested figure
%!test
%! hfig = figure ("visible", "off");
%! unwind_protect
%!   ax = axes ("parent", hfig);
%!   a = anova ([1;1;1;2;2;2;3;3;3], (1:9)', ...
%!              "SumOfSquaresType", "two");
%!   h = plotComparisons (a, ax);
%!   assert_equal (h, hfig);
%!   assert_equal (numel (findall (ax, "type", "line")) >= 6, true);
%! unwind_protect_cleanup
%!   close (hfig);
%! end_unwind_protect

## plotComparisons(): Dunn-Sidak intervals are wider than unadjusted LSD
%!test
%! hfig = figure ("visible", "off");
%! unwind_protect
%!   ax1 = subplot (1, 2, 1, "parent", hfig);
%!   ax2 = subplot (1, 2, 2, "parent", hfig);
%!   a = anova (kron ((1:4)', ones (3, 1)), (1:12)', ...
%!              "SumOfSquaresType", "two");
%!   plotComparisons (a, ax1, "CriticalValueType", "lsd");
%!   plotComparisons (a, ax2, "CriticalValueType", "dunn-sidak");
%!   lsd = findobj (ax1, "type", "line", "linestyle", "-");
%!   sidak = findobj (ax2, "type", "line", "linestyle", "-");
%!   lsd_width = diff (get (lsd(1), "xdata"));
%!   sidak_width = diff (get (sidak(1), "xdata"));
%!   assert_equal (sidak_width > lsd_width, true);
%! unwind_protect_cleanup
%!   close (hfig);
%! end_unwind_protect

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

## --- Week 6: predict / effect sizes -------------------------------------

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
%!  anova ([1;1;2;2], [1;2;3;4], 'ModelSpecification', 'cubic')

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
%! assert_equal (T.DF, [1; 2; 3]);
%! assert_equal (T.F(1), NaN);
%! assert_equal (T.pValue(1), NaN);

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
%!    columns_ = [2, 3, 4, 5, 6];     ## anova1 and anova2 layout
%!  else
%!    columns_ = [2, 3, 5, 6, 7];     ## anovan layout, Mean Sq. after Singular?
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
%! assert_equal (char (a.SumOfSquaresType), 'three');
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
%!   switch (k)
%!     case 2
%!       [~, ATAB, STATS] = anovan (score(:), ...
%!           {treatment(:), subject(:)}, 'model', 'full', 'sstype', 2, ...
%!           'varnames', {'treatment', 'subject'}, 'display', 'off');
%!     case 4
%!       [~, ATAB, STATS] = anovan (words(:), ...
%!           {seconds(:), subject(:)}, 'model', 'full', 'sstype', 2, ...
%!           'varnames', {'seconds', 'subject'}, 'display', 'off');
%!     case 9
%!       [~, ATAB, STATS] = anovan (measurement / 10, ...
%!           {strain, treatment, block}, 'model', 'full', 'sstype', 2, ...
%!           'varnames', {'strain', 'treatment', 'block'}, ...
%!           'display', 'off');
%!   endswitch
%!   actual = stats (a);
%!   actual_sources = actual.Properties.RowNames;
%!   expected_sources = ATAB(2:end, 1);
%!   expected_sources = regexprep (expected_sources, "X([0-9]+)", "Factor$1");
%!   if (strcmp (a.Stats.source, 'anova1'))
%!     actual_sources{1} = expected_sources{1};
%!   endif
%!   assert_equal (actual_sources, expected_sources);
%!   actual_values = __anova_values__ (actual);
%!   expected_values = __anova_values__ (ATAB);
%!   if (any (k == [2, 4, 9]))
%!     assert_equal (actual_values(:, 1:3), expected_values(:, 1:3), 1e-8);
%!   else
%!     assert_equal (actual_values, expected_values, 1e-8);
%!   endif
%!   if (! any (k == [1, 3, 12]))
%!     if (k == 13)
%!       expected_residuals = STATS.Y - full (STATS.X) * STATS.coeffs(:, 1);
%!     else
%!       expected_residuals = STATS.resid;
%!     endif
%!     assert_equal (a.Residuals.Raw, expected_residuals, 1e-8);
%!     assert_equal (size (a.DesignMatrix), size (STATS.X));
%!   endif
%! endfor

## MATLAB-compatible table constructors and public result schemas.
%!test
%! dose = [1; 1; 2; 2; 1; 2];
%! site = {'A'; 'B'; 'A'; 'B'; 'A'; 'B'};
%! yield = [1; 2; 4; 5; 2; 6];
%! tbl = table (dose, site, yield, ...
%!              'VariableNames', {'Dose', 'Site', 'Yield'});
%! a = anova (tbl, 'Yield');
%! assert_equal (cellstr (a.FactorNames), {'Dose', 'Site'});
%! assert_equal (a.ResponseName, 'Yield');
%! assert_equal (istable (a.Factors), true);
%! assert_equal (a.Factors.Properties.VariableNames, {'Dose', 'Site'});
%! assert_equal (a.Y, yield);
%! T = stats (a);
%! direct = stats (anova ({dose, site}, yield, ...
%!                        'FactorNames', {'Dose', 'Site'}));
%! assert_equal (__anova_values__ (T), __anova_values__ (direct), 1e-12);

%!test
%! dose = [1; 1; 2; 2; 1; 2];
%! site = {'A'; 'B'; 'A'; 'B'; 'A'; 'B'};
%! yield = [1; 2; 4; 5; 2; 6];
%! tbl = table (dose, site, yield, ...
%!              'VariableNames', {'Dose', 'Site', 'Yield'});
%! a = anova (tbl, 'Yield ~ Dose + Site + Dose:Site');
%! T = stats (a);
%! direct = stats (anova ({dose, site}, yield, ...
%!                        'FactorNames', {'Dose', 'Site'}, ...
%!                        'ModelSpecification', 'full'));
%! assert_equal (char (a.Formula.Text), ...
%!               'Yield ~ Dose + Site + Dose:Site');
%! assert_equal (a.ModelSpecification, char (a.Formula.Text));
%! assert_equal (__anova_values__ (T), __anova_values__ (direct), 1e-12);

%!test
%! dose = [1; 1; 2; 2];
%! site = [1; 2; 1; 2];
%! y = [1; 2; 4; 5];
%! tbl = table (dose, site, 'VariableNames', {'Dose', 'Site'});
%! a = anova (tbl, y, 'FactorNames', {'Site'});
%! assert_equal (cellstr (a.FactorNames), {'Site'});
%! assert_equal (a.Factors.Site, site);
%! assert_equal (stats (a).DF, [1; 2; 3]);

%!test
%! g1 = [1; 1; 1; 1; 2; 2; 2; 2];
%! g2 = [1; 1; 2; 2; 1; 1; 2; 2];
%! y = (1:8)';
%! a = anova ({g1, g2}, y, 'FactorNames', {'A', 'B'}, ...
%!            'CategoricalFactors', [true, false], 'RandomFactors', 'A');
%! assert_equal (a.CategoricalFactors, 1);
%! assert_equal (a.RandomFactors, 1);
%! b = anova ({g1, g2}, y, 'FactorNames', {'A', 'B'}, ...
%!            'CategoricalFactors', {'A', 'B'}, 'RandomFactors', 'all');
%! assert_equal (b.CategoricalFactors, [1, 2]);
%! assert_equal (b.RandomFactors, [1, 2]);

%!test
%! y = [1; 2; 3; 4; 5; 6];
%! g = [1; 1; 1; 2; 2; 2];
%! a = anova (g, y, 'SumOfSquaresType', 'two');
%! stats (a);
%! assert_equal (isvector (a.Coefficients), true);
%! assert_equal (isstring (a.ExpandedFactorNames), true);
%! assert_equal (istable (a.Residuals), true);
%! assert_equal (a.Residuals.Properties.VariableNames, {'Raw', 'Pearson'});
%! assert_equal (a.Residuals.Pearson, ...
%!               a.Residuals.Raw / sqrt (a.Metrics.MSE), 1e-12);
%! assert_equal (a.Metrics.Properties.VariableNames, ...
%!               {'MSE', 'RMSE', 'SSE', 'SSR', 'SST', ...
%!                'RSquared', 'AdjustedRSquared'});
%! assert_equal (a.Metrics.SST, a.Metrics.SSE + a.Metrics.SSR, 1e-12);
%! assert_equal (a.Coefficients, [3.5; -1.5; 1.5], 1e-12);
%! assert_equal (cellstr (a.ExpandedFactorNames), ...
%!               {'(Intercept)'; '(Factor1==1)'; '(Factor1==2)'});

%!test
%! y = [1; 2; 3; 4; 5; 6; 10; 11; 12];
%! g = [1; 1; 1; 2; 2; 2; 3; 3; 3];
%! a = anova (g, y);
%! m = multcompare (a);
%! assert_equal (istable (m), true);
%! assert_equal (m.Properties.VariableNames, ...
%!               {'Group1', 'Group2', 'MeanDifference', ...
%!                'MeanDifferenceLower', 'MeanDifferenceUpper', 'pValue'});
%! assert_equal (m.MeanDifference, [-3; -9; -6], 1e-12);
%! v = varianceComponent (a);
%! assert_equal (v.Properties.VariableNames, ...
%!               {'VarianceComponent', 'VarianceComponentLower', ...
%!                'VarianceComponentUpper'});
%! assert_equal (v.Properties.RowNames, {'Error'});
%! one_group = multcompare (anova (ones (4, 1), (1:4)'));
%! assert_equal (istable (one_group), true);
%! assert_equal (size (one_group), [0, 6]);

%!test
%! g1 = [1; 1; 1; 1; 2; 2; 2; 2];
%! g2 = {'A'; 'A'; 'B'; 'B'; 'A'; 'A'; 'B'; 'B'};
%! y = [1; 2; 3; 4; 5; 6; 8; 9];
%! a = anova ({g1, g2}, y, 'FactorNames', {'Dose', 'Site'}, ...
%!            'ModelSpecification', 'full');
%! m = multcompare (a, {'Dose', 'Site'}, ...
%!                  'CriticalValueType', 'bonferroni');
%! assert_equal (istable (m.Group1), true);
%! assert_equal (istable (m.Group2), true);
%! assert_equal (m.Group1.Properties.VariableNames, {'Dose', 'Site'});
%! assert_equal (rows (m), 6);

%!test
%! g1 = {"x, y"; "x, y"; "x, y"; "x, y"; "z=1"; "z=1"; "z=1"; "z=1"};
%! g2 = {"a=1"; "a=1"; "b,2"; "b,2"; "a=1"; "a=1"; "b,2"; "b,2"};
%! y = (1:8)';
%! a = anova ({g1, g2}, y, "FactorNames", {"Dose", "Site"}, ...
%!            "ModelSpecification", "full");
%! m = multcompare (a, {"Dose", "Site"});
%! values = [m.Group1.Dose; m.Group2.Dose; m.Group1.Site; m.Group2.Site];
%! assert_equal (any (strcmp (values, "x, y")), true);
%! assert_equal (any (strcmp (values, "z=1")), true);
%! assert_equal (any (strcmp (values, "a=1")), true);
%! assert_equal (any (strcmp (values, "b,2")), true);

%!test
%! g1 = [1; 1; 1; 1; 2; 2; 2; 2];
%! g2 = [1; 1; 2; 2; 1; 1; 2; 2];
%! y = (1:8)';
%! a = anova ({g1, g2}, y, 'FactorNames', {'A', 'B'}, ...
%!            'ResponseName', 'R', ...
%!            'ModelSpecification', 'R ~ A + B + A:B');
%! assert_equal (char (a.Formula.Text), 'R ~ A + B + A:B');
%! assert_equal (stats (a).Properties.RowNames, ...
%!               {'A'; 'B'; 'A:B'; 'Error'; 'Total'});
%! b = anova ({g1, g2}, y, 'FactorNames', {'A', 'B'}, ...
%!            'ModelSpecification', 3);
%! assert_equal (b.ModelSpecification, 3);
%! assert_equal (char (b.Formula.Text), 'Y ~ 1 + A + B + A:B');
%! assert_equal (__anova_values__ (stats (a)), ...
%!               __anova_values__ (stats (b)), 1e-12);

%!test
%! tbl = table ([1; 1; 2; 2], [1; 2; 3; 4], ...
%!              'VariableNames', {'Group', 'Yield'});
%! a = anova (tbl, 'Yield ~ Group', 'FactorNames', {'Yield'}, ...
%!            'ResponseName', 'Ignored', 'ModelSpecification', 'full');
%! assert_equal (cellstr (a.FactorNames), {'Group'});
%! assert_equal (a.ResponseName, 'Yield');
%! assert_equal (a.ModelSpecification, 'Yield ~ Group');

%!test
%! a = anova ([1; 1; 2; 2], [1; NaN; 3; 4]);
%! stats (a);
%! assert_equal (a.NumObservations, 3);
%! assert_equal (a.Y, [1; 3; 4]);
%! assert_equal (height (a.Factors), 3);
%! assert_equal (height (a.Residuals), 3);

## Random-factor variance components match MATLAB's documented carsmall example.
%!test
%! load carsmall
%! a = anova ({Origin, Model_Year}, MPG, "RandomFactors", [1, 2], ...
%!            "FactorNames", {"Origin", "Year"});
%! v = varianceComponent (a);
%! assert_equal (v.Properties.RowNames, {"Origin"; "Year"; "Error"});
%! assert_equal (v.VarianceComponent, [21.337; 44.031; 20.198], 5e-3);
%! assert_equal (v.VarianceComponentLower, [6.1257; 11.176; 15.298], 5e-3);
%! assert_equal (v.VarianceComponentUpper, [139.94; 1765.7; 27.909], 5e-1);

## Interactions containing a random factor remain in the model and provide
## the denominator for fixed-effect tests.
%!test
%! A = kron ([1; 2], ones (12, 1));
%! B = repmat (kron ([1; 2; 3], ones (4, 1)), 2, 1);
%! y = 10 + 2 * (A == 2) + ...
%!     [0;1;-1;0; 1;0;-1;0; -1;0;1;0; 0;1;0;-1; 2;1;0;1; -1;0;1;0];
%! a = anova ({A, B}, y, "FactorNames", {"A", "B"}, ...
%!            "ModelSpecification", "full", "RandomFactors", 2);
%! [s, ems] = stats (a);
%! assert_equal (s.Properties.RowNames, {"A"; "B"; "A:B"; "Error"; "Total"});
%! assert_equal (s.F(1:3), [49; 1; 1], 1e-12);
%! assert_equal (cellstr (ems.Type), {"fixed"; "random"; "random"; "random"});
%! assert_equal (cellstr (ems.FDenominator), ...
%!               {"MS(A:B)"; "MS(A:B)"; "MS(Error)"; ""});
%! v = varianceComponent (a);
%! assert_equal (v.Properties.RowNames, {"B"; "A:B"; "Error"});

## A negative variance-component estimate is reported without an interval.
%!test
%! A = kron ([1; 2], ones (12, 1));
%! B = repmat (kron ([1; 2; 3; 4], ones (3, 1)), 2, 1);
%! y = [5;6;4; 12;13;11; 20;19;21; 8;7;9; ...
%!      7;8;6; 15;14;16; 23;22;24; 10;11;9];
%! a = anova ({A, B}, y, 'FactorNames', {'A', 'B'}, ...
%!            'ModelSpecification', 'full', 'RandomFactors', 2);
%! v = varianceComponent (a);
%! assert_equal (v.Properties.RowNames, {'B'; 'A:B'; 'Error'});
%! assert_equal (v.VarianceComponent, ...
%!               [45.4166666666666; -0.166666666666666; 1], 1e-9);
%! assert_equal (v.VarianceComponentLower, ...
%!               [13.4429180926305; NaN; 0.554682109897800], 1e-9);
%! assert_equal (v.VarianceComponentUpper, ...
%!               [632.517205323886; NaN; 2.31626772541430], 1e-8);

## A non-negative estimate has its lower bound truncated at zero.
%!test
%! A = kron ([1; 2], ones (12, 1));
%! B = repmat (kron ([1; 2; 3; 4], ones (3, 1)), 2, 1);
%! y = [5;6;4; 12;13;11; 20;19;21; 8;7;9; ...
%!      7;8;6; 15;14;16; 23;22;24; 10;11;9];
%! y(1:12) += [0;0;0; 1;1;1; -1;-1;-1; 0.4;0.4;0.4];
%! a = anova ({A, B}, y, 'FactorNames', {'A', 'B'}, ...
%!            'ModelSpecification', 'full', 'RandomFactors', 2);
%! v = varianceComponent (a);
%! assert_equal (v.VarianceComponent(2), 0.253333333333400, 1e-9);
%! assert_equal (v.VarianceComponentLower(2), 0);
%! assert_equal (v.VarianceComponentUpper(2), 7.97098397237600, 1e-9);

## Named polynomial models preserve MATLAB's term definitions.
%!test
%! x = (-2:2)';
%! y = 2 + 3 * x + 4 * x .^ 2;
%! a = anova (x, y, "CategoricalFactors", [], ...
%!            "ModelSpecification", "purequadratic", ...
%!            "FactorNames", {"x"});
%! s = stats (a);
%! assert_equal (a.Stats.terms, [1; 2]);
%! assert_equal (s.Properties.RowNames, {"x"; "x^2"; "Error"; "Total"});
%! assert_equal (a.FittedValues, y, 1e-12);

## Hierarchical sums preserve polynomial power containment.
%!test
%! x = [0; 1; 2; 4; 7; 8; 10; 15];
%! y = [1; 2; 3; 6; 12; 15; 21; 40];
%! a = anova (x, y, "CategoricalFactors", [], ...
%!            "ModelSpecification", [1; 2], ...
%!            "SumOfSquaresType", "hierarchical");
%! s = stats (a);
%! assert_equal (char (a.SumOfSquaresType), "hierarchical");
%! assert_equal (a.SSType, "h");
%! assert_equal (s.SumOfSquares(1:2), ...
%!               [1149.540669856459; 60.32250042332054], 1e-10);

## Accepted aliases are canonicalized in the public property.
%!test
%! a = anova ([1; 1; 2; 2], (1:4)', "SumOfSquaresType", "typeii");
%! assert_equal (char (a.SumOfSquaresType), "two");

%!test
%! [x1, x2] = ndgrid ((-2:2)', (-1:1)');
%! x1 = x1(:);
%! x2 = x2(:);
%! y = 1 + 2*x1 + 3*x2 + 4*x1.*x2 + 5*x1.^2 + 6*x2.^2;
%! a = anova ({x1, x2}, y, "CategoricalFactors", [], ...
%!            "ModelSpecification", "quadratic");
%! stats (a);
%! assert_equal (a.Stats.terms, [1 0; 0 1; 1 1; 2 0; 0 2]);
%! assert_equal (a.FittedValues, y, 1e-10);

%!test
%! [x1, x2] = ndgrid ([-1; 1], (-2:2)');
%! x1 = x1(:);
%! x2 = x2(:);
%! y = x1 + x2.^2 + x2.^3 + x1.*x2 + x1.*x2.^2;
%! a = anova ({x1, x2}, y, "CategoricalFactors", [], ...
%!            "ModelSpecification", "poly13");
%! stats (a);
%! assert_equal (a.Stats.terms, [1 0; 0 1; 0 2; 0 3; 1 1; 1 2]);

%!test
%! x = (-2:2)';
%! tbl = table (x, 2 + 3*x + 4*x.^2, "VariableNames", {"x", "y"});
%! a = anova (tbl, "y ~ x^2", "CategoricalFactors", []);
%! stats (a);
%! assert_equal (a.Stats.terms, [1; 2]);
%! assert_equal (a.FittedValues, tbl.y, 1e-12);

## Nested factors use separate within-parent contrasts.
%!test
%! A = kron ([1; 2], ones (6, 1));
%! B = repmat (kron ([1; 2], ones (3, 1)), 2, 1);
%! y = 10*A + 2*B + repmat ([-1; 0; 1], 4, 1);
%! tbl = table (A, B, y, "VariableNames", {"A", "B", "Y"});
%! a = anova (tbl, "Y ~ A + B(A)");
%! s = stats (a);
%! assert_equal (a.Stats.terms, [1 0; 0 1]);
%! assert_equal (a.Stats.vnested, logical ([0 0; 1 0]));
%! assert_equal (s.Properties.RowNames, {"A"; "B(A)"; "Error"; "Total"});
%! assert_equal (s.DF, [1; 2; 8; 11]);

%!test
%! [A, C, B, R] = ndgrid ([1; 2], [1; 2], [1; 2], [1; 2]);
%! A = A(:); C = C(:); B = B(:); R = R(:);
%! y = A + 2*C + 3*B + 0.1*R;
%! tbl = table (A, B, C, y, "VariableNames", {"A", "B", "C", "Y"});
%! a = anova (tbl, "Y ~ A*C + B(A,C)");
%! s = stats (a);
%! assert_equal (a.Stats.vnested, logical ([0 0 0; 1 0 1; 0 0 0]));
%! assert_equal (s.Properties.RowNames, ...
%!               {"A"; "C"; "A:C"; "B(A,C)"; "Error"; "Total"});
%! assert_equal (s.DF(4), 4);

## Unbalanced NaN cleanup preserves the anovan fallback results.
%!test
%! x = [1, 5; 2, 6; 3, 7; NaN, 8];
%! warning ("off", "all", "local");
%! a = anova (x, [], "Reps", 2, "ModelSpecification", "interactions");
%! stats (a);
%! assert_equal (a.Stats.source, "anovan");
%! assert_equal (a.NumObservations, 7);
%! assert_equal (height (a.Residuals), 7);

## Numeric replicated models bypass anova2 when it cannot represent the terms.
%!test
%! x = [1, 5; 2, 6; 3, 7; 4, 8];
%! a = anova (x, [], "Reps", 2, ...
%!            "ModelSpecification", [1 0; 0 1], ...
%!            "FactorNames", {"A", "B"});
%! s = stats (a);
%! assert_equal (s.Properties.RowNames, {"A"; "B"; "Error"; "Total"});

## Raw residuals are unweighted observation-minus-fit differences.
%!test
%! g = [1; 1; 2; 2; 3; 3];
%! y = [1; 3; 4; 8; 9; 15];
%! a = anova (g, y, "Weights", [1; 2; 1; 3; 1; 4]);
%! stats (a);
%! assert_equal (a.Residuals.Raw, a.Y - a.FittedValues, 1e-12);

%!error <polynomial powers require continuous> ...
%! anova ([1; 1; 2; 2], (1:4)', "ModelSpecification", "purequadratic")

%!error <anova: RandomFactors must contain valid factor indices.> ...
%!  anova ([1;1;2;2], [1;2;3;4], "RandomFactors", 0)

%!error <anova: RandomFactors must contain valid factor indices.> ...
%!  anova ([1;1;2;2], [1;2;3;4], "RandomFactors", 2)

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

## A saturated mixed model still yields every F ratio and p-value.
%!test
%! y = [444 614 423 625 408 856 447 719 ...
%!      764 831 586 782 609 1002 606 766]' / 10;
%! X1 = {'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola', ...
%!       'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola'}';
%! X2 = {'C','T','C','T','C','T','C','T','C','T','C','T','C','T','C','T'}';
%! X3 = [1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2];
%! a = anova ({X1, X2, X3}, y, 'ModelSpecification', 'full', ...
%!            'RandomFactors', 3, 'FactorNames', {'X1', 'X2', 'X3'});
%! s = stats (a);
%! assert_equal (s.DF', [3, 1, 1, 3, 3, 1, 3, 0, 15]);
%! assert_equal (s.SumOfSquares(8), 0);
%! assert_equal (s.pValue(1:6), ...
%!               [0.288811428913179; 0.091455278902114; 0.042134889025806; ...
%!                0.010944863181481; 0.061814763198376; 0.066584161056625], ...
%!               1e-12);
%! assert_equal (s.pValue(7), NaN);

## The denominator of each F ratio comes from the expected mean squares.
%!test
%! y = [444 614 423 625 408 856 447 719 ...
%!      764 831 586 782 609 1002 606 766]' / 10;
%! X1 = {'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola', ...
%!       'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola'}';
%! X2 = {'C','T','C','T','C','T','C','T','C','T','C','T','C','T','C','T'}';
%! X3 = [1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2];
%! a = anova ({X1, X2, X3}, y, 'ModelSpecification', 'full', ...
%!            'RandomFactors', 3, 'FactorNames', {'X1', 'X2', 'X3'});
%! [~, ems] = stats (a);
%! assert_equal (ems.MeanSquaresDenominator(1:7), ...
%!               [47.1575; 47.61; 88.7925; 5.975; 5.975; 5.975; 0], 1e-9);
%! assert_equal (ems.DFDenominator(1:7), ...
%!               [3; 1; 2.610727841363640; 3; 3; 3; 0], 1e-12);

## A term contained in another contributes to that term's expected mean square,
## whether the containing term is fixed or random.
%!test
%! y = [444 614 423 625 408 856 447 719 ...
%!      764 831 586 782 609 1002 606 766]' / 10;
%! X1 = {'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola', ...
%!       'NIH','NIH','BALB/C','BALB/C','A/J','A/J','129/Ola','129/Ola'}';
%! X2 = {'C','T','C','T','C','T','C','T','C','T','C','T','C','T','C','T'}';
%! X3 = [1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2];
%! a = anova ({X1, X2, X3}, y, 'ModelSpecification', 'full', ...
%!            'RandomFactors', 3, 'FactorNames', {'X1', 'X2', 'X3'});
%! [~, ems] = stats (a);
%! assert_equal (cellstr (ems.ExpectedMeanSquares), ...
%!   {'4*Q(X1)+2*Q(X1:X2)+2*V(X1:X3)+V(X1:X2:X3)+V(Error)'; ...
%!    '8*Q(X2)+2*Q(X1:X2)+4*V(X2:X3)+V(X1:X2:X3)+V(Error)'; ...
%!    '8*V(X3)+2*V(X1:X3)+4*V(X2:X3)+V(X1:X2:X3)+V(Error)'; ...
%!    '2*Q(X1:X2)+V(X1:X2:X3)+V(Error)'; ...
%!    '2*V(X1:X3)+V(X1:X2:X3)+V(Error)'; ...
%!    '4*V(X2:X3)+V(X1:X2:X3)+V(Error)'; ...
%!    'V(X1:X2:X3)+V(Error)'; 'V(Error)'});

## A saturated model reports no residual rather than an infinite one.
%!test
%! y = [3; 4; 7; 8];
%! [~, tbl] = anovan (y, {[1;1;2;2], [1;2;1;2]}, 'model', 'full', ...
%!                    'display', 'off');
%! assert_equal (tbl{end-1, 2}, 0);
%! assert_equal (tbl{end-1, 3}, 0);
%! assert_equal (tbl{end-1, 5}, 0);
