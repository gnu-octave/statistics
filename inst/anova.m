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
          otherwise
            error ("anova: parameter '%s' is not supported.", name);
        endswitch
      endfor

      obj.validateSpec_ ();
      obj.nFactors_ = obj.countFactors_ ();
      obj.selectBackend_ ();
      obj.dirty_ = true;
    endfunction

  endmethods

  methods (Access = private)

    ## Validate model-spec property values that the constructor accepts.
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

    ## Backend heuristic (verbatim from the proposal):
    ##   anova1  : 1 factor, no continuous, no weights, SSType == 3
    ##   anova2  : 2 factors, no continuous, Y is a matrix, no NaN, no weights
    ##   anovan  : everything else
    function selectBackend_ (obj)
      if (obj.nFactors_ == 1 && isempty (obj.Continuous) ...
          && isempty (obj.Weights) && obj.SSType == 3)
        obj.backend_ = 'anova1';
      elseif (obj.nFactors_ == 2 && isempty (obj.Continuous) ...
              && ismatrix (obj.Y) && ! any (isnan (obj.Y(:))) ...
              && isempty (obj.Weights))
        obj.backend_ = 'anova2';
      else
        obj.backend_ = 'anovan';
      endif
    endfunction

    ## Stubs for Week 2 — declared so the structure is visible but they
    ## intentionally do not synthesize results in Week 1.
    function fit_ (obj)
      error ("anova.fit_: not implemented yet (scheduled for Week 2).");
    endfunction

    function ensureFit_ (obj)
      if (! obj.fitted_ || obj.dirty_)
        obj.selectBackend_ ();
        obj.fit_ ();
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

## Backend selection: two-way balanced matrix -> anova2
%!test
%! y = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!      6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! g1 = repmat ([1;2], 3, 1);
%! g2 = [1;1;1;2;2;2];
%! a = anova (y, {g1, g2});
%! assert (a.getNumFactors (), 2);
%! assert (a.getBackend (), 'anova2');

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
