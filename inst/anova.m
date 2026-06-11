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
  ## The @code{anova} class provides a unified, stateful, object-oriented
  ## wrapper around the procedural functions @code{anova1}, @code{anova2}, and
  ## @code{anovan}.  The constructor stores the data and model specification,
  ## selects the appropriate backend, and (in later stages) delegates the
  ## actual computation to the chosen procedural function.
  ##
  ## This is the Week 1 skeleton: construction, validation, name-value
  ## parsing, and backend selection only.  Fit delegation and result
  ## population are introduced in subsequent stages.
  ##
  ## @seealso{anova1, anova2, anovan, multcompare}
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
    ModelType   = 'linear';                 ## 'linear', 'interaction', 'full', or terms matrix
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
  endproperties

  methods (Access = public)

    function obj = anova (Y, GROUP, varargin)

      if (nargin < 1)
        error ("anova: too few input arguments.");
      endif
      if (nargin < 2)
        GROUP = [];
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

      obj.validateSpec_ ();
      obj.nFactors_ = obj.countFactors_ ();
      obj.selectBackend_ ();
      obj.dirty_ = true;
    endfunction

    function fit (obj)
      obj.ensureFit_ ();
    endfunction

  endmethods

  methods (Access = private)

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
      if (! isempty (obj.Random) && ! isnumeric (obj.Random))
        error ("anova: Random must be a numeric index vector.");
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

  endmethods

  ## Test-only accessors for the private backend selection state.
  ## Kept hidden to discourage external use while still letting BISTs
  ## verify Week 1 behaviour without piercing object privacy by hand.
  methods (Access = public, Hidden)
    function b = getBackend (obj)
      b = obj.backend_;
    endfunction
    function n = getNumFactors (obj)
      n = obj.nFactors_;
    endfunction
  endmethods

endclassdef

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
%! assert (a.getNumFactors (), 2);

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
%! assert (a.getBackend (), 'anova1');

## Backend selection: one-way matrix-Y form -> anova1
%!test
%! a = anova (magic (4));
%! assert (a.getNumFactors (), 1);
%! assert (a.getBackend (), 'anova1');

## Backend selection: matrix Y + explicit 'reps' -> anova2
%!test
%! y = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!      6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! a = anova (y, [], 'reps', 3);
%! assert (a.getBackend (), 'anova2');

## Backend selection: two-factor cell groups without reps -> anovan
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = repmat ([1;1;2;2], 3, 1);
%! a = anova (y, {g1, g2});
%! assert (a.getNumFactors (), 2);
%! assert (a.getBackend (), 'anovan');

## Backend selection: three factors -> anovan
%!test
%! y = (1:24)';
%! g1 = repmat ([1;2], 12, 1);
%! g2 = repmat ([1;1;2;2], 6, 1);
%! g3 = repmat ([1;1;1;1;2;2;2;2], 3, 1);
%! a = anova (y, {g1, g2, g3});
%! assert (a.getNumFactors (), 3);
%! assert (a.getBackend (), 'anovan');

## Backend selection: SSType != 3 with 1 factor falls through to anovan
%!test
%! a = anova ([1;2;3;4;5;6], [1;1;2;2;3;3], 'SSType', 2);
%! assert (a.getBackend (), 'anovan');

## Backend selection: continuous predictors force anovan
%!test
%! y = (1:12)';
%! g1 = repmat ([1;2;3], 4, 1);
%! g2 = (1:12)';                              ## continuous
%! a = anova (y, {g1, g2}, 'Continuous', 2);
%! assert (a.getBackend (), 'anovan');

## Backend selection: NaN in Y demotes 2-way to anovan
%!test
%! y = [1, 2; 3, NaN; 5, 6; 7, 8];
%! g1 = [1;2;1;2];
%! g2 = [1;1;2;2];
%! a = anova (y, {g1, g2});
%! assert (a.getBackend (), 'anovan');

## Backend selection: weights force anovan
%!test
%! a = anova ([1;2;3;4;5;6], [1;1;2;2;3;3], 'Weights', ones (6, 1));
%! assert (a.getBackend (), 'anovan');

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
%! assert (a.getBackend (), 'anova2');
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
%! assert (a.getBackend (), 'anovan');
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
%! assert (a.getBackend (), 'anovan');
%! assert (a.Stats.alpha, 0.10, 1e-12);
