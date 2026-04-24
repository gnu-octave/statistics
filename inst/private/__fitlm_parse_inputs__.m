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

## -*- texinfo -*-
## @deftypefn {Private Function} {[@var{raw_X}, @var{raw_y}, @var{var_names}, @
##   @var{cat_flags}, @var{weights}, @var{excl_mask}, @var{modelspec}, @
##   @var{obs_names}, @var{resp_name}, @var{pred_names}, @var{has_intercept}] =} @
##   __fitlm_parse_inputs__ (@var{varargin})
##
## Parse fitlm inputs into clean typed outputs.
##
## Returns:
##   raw_X        - n x p numeric matrix (NaN in missing rows preserved)
##   raw_y        - n x 1 numeric response vector
##   var_names    - 1 x (p+1) cell of char: {pred_names..., resp_name}
##   cat_flags    - 1 x (p+1) logical: IsCategorical for each variable
##   weights      - n x 1 double (1.0 if unweighted)
##   excl_mask    - n x 1 logical (true = explicitly excluded)
##   modelspec    - char: Wilkinson formula string OR shorthand string
##   obs_names    - cell of char (row names or {})
##   resp_name    - char: response variable name
##   pred_names   - 1 x p cell of char: predictor variable names
##   has_intercept - logical: true unless explicitly suppressed
##   raw_cols     - cell of original columns (for VariableInfo building)
## @end deftypefn

function [raw_X, raw_y, var_names, cat_flags, weights, excl_mask, ...
          modelspec, obs_names, resp_name, pred_names, has_intercept, ...
          raw_cols, dummy_coding] = ...
    __fitlm_parse_inputs__ (varargin)

  if (numel (varargin) < 1)
    error ("fitlm: at least one input argument required.");
  endif

  ## --- Detect first argument: table or matrix ---
  is_table_input = isa (varargin{1}, "table");

  ## Initialize name-value defaults
  weights_arg     = [];
  excl_arg        = [];
  intercept_arg   = true;
  cat_vars_arg    = [];     ## integer indices or {}
  var_names_arg   = {};     ## VarNames override
  resp_var_arg    = "";     ## ResponseVar override
  pred_vars_arg   = {};     ## PredictorVars override

  if (is_table_input)
    ## ── TABLE INPUT ──────────────────────────────────────────────────────────
    tbl = varargin{1};
    varargin(1) = [];

    ## Extract table metadata
    all_col_names = tbl.Properties.VariableNames;
    n_vars_tbl    = numel (all_col_names);

    ## ObservationNames from row names
    obs_names = tbl.Properties.RowNames;     ## cell or {} if unset

    ## Parse optional modelspec (2nd positional arg if char/not a NV key,
    ## or a numeric terms matrix)
    modelspec = "linear";                    ## default shorthand
    if (! isempty (varargin) && isnumeric (varargin{1}) && ...
        ismatrix (varargin{1}) && ndims (varargin{1}) == 2)
      ## Numeric terms matrix passed directly
      modelspec = varargin{1};
      varargin(1) = [];
    elseif (! isempty (varargin) && ischar (varargin{1}) && ...
        ! __is_nv_key__ (varargin{1}))
      modelspec = varargin{1};
      varargin(1) = [];
    endif

    ## Parse name-value pairs
    [weights_arg, excl_arg, intercept_arg, cat_vars_arg, var_names_arg, ...
     resp_var_arg, pred_vars_arg, dummy_coding_arg] = __parse_nv__ (varargin{:});

    ## Determine response and predictor variables
    if (! isempty (resp_var_arg))
      resp_name = resp_var_arg;
    elseif (! isempty (strfind (modelspec, "~")))
      ## Formula string — response is LHS
      tilde_pos = strfind (modelspec, "~");
      resp_name = strtrim (modelspec(1:tilde_pos(1)-1));
    else
      ## Default: last column of table
      resp_name = all_col_names{end};
    endif

    ## Predictor names
    if (! isempty (pred_vars_arg))
      if (iscell (pred_vars_arg))
        pred_names = pred_vars_arg;
      else
        pred_names = all_col_names(pred_vars_arg);
      endif
    else
      ## All columns except response
      pred_names = all_col_names(! strcmp (all_col_names, resp_name));
    endif

    ## var_names: predictor order then response (matches MATLAB Block 4-P4)
    var_names = [pred_names, {resp_name}];

    ## Extract raw columns
    n = rows (tbl);
    p = numel (pred_names);
    raw_X = zeros (n, p);
    raw_cols = cell (1, p + 1);

    cat_flags = false (1, p + 1);

    ## Handle CategoricalVars — may be names or indices into pred_names
    cat_idx_set = __resolve_cat_vars__ (cat_vars_arg, pred_names);

    for k = 1:p
      col = tbl.(pred_names{k});
      raw_cols{k} = col;
      if (isa (col, "categorical"))
        ## Store numeric codes — raw_X stays as codes for design builder
        cat_flags(k) = true;
        raw_X(:, k) = double (col);   ## integer codes (1-based after sort)
      elseif (iscellstr (col) || isstring (col))
        cat_flags(k) = true;
        raw_X(:, k) = NaN;            ## placeholder; design builder uses raw_cols
      elseif (ismember (k, cat_idx_set))
        cat_flags(k) = true;
        raw_X(:, k) = col;
      else
        raw_X(:, k) = col;
      endif
    endfor

    ## Response column
    raw_y = tbl.(resp_name);
    raw_y = raw_y(:);
    raw_cols{p+1} = raw_y;

    ## cat_flags for response (always false)
    cat_flags(p+1) = false;

    ## No VarNames override for table input

  else
    ## ── MATRIX INPUT ─────────────────────────────────────────────────────────
    X_in = varargin{1};
    if (nargin < 2 || numel (varargin) < 2)
      error ("fitlm: matrix input requires at least X and y arguments.");
    endif
    y_in = varargin{2}(:);
    varargin(1:2) = [];

    n = rows (X_in);
    p = columns (X_in);
    obs_names = {};

    ## Parse optional modelspec (char string, shorthand, or numeric terms matrix)
    modelspec = "linear";
    if (! isempty (varargin) && isnumeric (varargin{1}) && ...
        ismatrix (varargin{1}) && ndims (varargin{1}) == 2)
      ## Numeric terms matrix passed directly
      modelspec = varargin{1};
      varargin(1) = [];
    elseif (! isempty (varargin) && ischar (varargin{1}) && ...
        ! __is_nv_key__ (varargin{1}))
      modelspec = varargin{1};
      varargin(1) = [];
    endif

    ## Parse name-value pairs
    [weights_arg, excl_arg, intercept_arg, cat_vars_arg, var_names_arg, ...
     resp_var_arg, pred_vars_arg, dummy_coding_arg] = __parse_nv__ (varargin{:});

    ## Variable names
    if (! isempty (var_names_arg))
      ## VarNames override: last = response, rest = predictors
      if (numel (var_names_arg) != p + 1)
        error ("fitlm: VarNames must have %d elements (p predictors + 1 response).", p+1);
      endif
      pred_names = var_names_arg(1:p);
      resp_name  = var_names_arg{p+1};
    else
      ## Auto-generate: x1...xp and y
      pred_names = arrayfun (@(j) sprintf ("x%d", j), 1:p, "UniformOutput", false);
      resp_name  = "y";
    endif

    var_names = [pred_names, {resp_name}];

    raw_X = X_in;
    raw_y = y_in;

    ## cat_flags — start all false, then mark CategoricalVars
    cat_flags = false (1, p + 1);
    cat_idx_set = __resolve_cat_vars__ (cat_vars_arg, pred_names);
    if (! isempty (cat_idx_set))
      cat_flags(cat_idx_set) = true;
    endif

    ## raw_cols for VariableInfo building
    raw_cols = cell (1, p + 1);
    for k = 1:p
      raw_cols{k} = raw_X(:, k);
    endfor
    raw_cols{p+1} = raw_y;

  endif  ## is_table_input

  ## --- Weights ---
  if (isempty (weights_arg))
    weights = ones (n, 1);
  else
    weights = weights_arg(:);
    if (numel (weights) != n)
      error ("fitlm: Weights must have %d elements.", n);
    endif
  endif

  ## --- Exclusion mask ---
  excl_mask = false (n, 1);
  if (! isempty (excl_arg))
    if (islogical (excl_arg))
      excl_mask = logical (excl_arg(:));
    else
      idx = round (excl_arg(:));
      excl_mask(idx) = true;
    endif
  endif

  ## --- Intercept ---
  has_intercept = logical (intercept_arg);

  ## If modelspec is a formula string, extract has_intercept from it
  if (ischar (modelspec) && ! isempty (strfind (modelspec, "~")))
    has_intercept = isempty (regexp (modelspec, '-\s*1'));
    if (! has_intercept)
      ## '-1' suppresses intercept; strip it for later parsing
    endif
  endif


  ## DummyVarCoding
  dummy_coding = dummy_coding_arg;

  ## If Intercept=false was explicitly given, override
  if (! intercept_arg)
    has_intercept = false;
  endif

endfunction

## ── Internal helpers ─────────────────────────────────────────────────────────

function result = __is_nv_key__ (s)
  ## Return true if s is a known name-value argument name.
  ## Everything that is NOT a known key and contains '~' or is a shorthand
  ## is passed as modelspec.
  known = {"weights", "exclude", "intercept", "categoricalvars", ...
           "varnames", "robustopts", "responsevar", "predictorvars", ...
           "dummyvarcoding"};
  result = any (strcmpi (s, known));
endfunction

function [weights_arg, excl_arg, intercept_arg, cat_vars_arg, var_names_arg, ...
          resp_var_arg, pred_vars_arg, dummy_coding_arg] = __parse_nv__ (varargin)
  weights_arg     = [];
  excl_arg        = [];
  intercept_arg   = true;
  cat_vars_arg    = [];
  var_names_arg   = {};
  resp_var_arg    = "";
  pred_vars_arg   = {};
  dummy_coding_arg = "reference";

  i = 1;
  while (i <= numel (varargin))
    if (! ischar (varargin{i}))
      error ("fitlm: expected name-value argument name at position %d.", i);
    endif
    key = lower (varargin{i});
    if (i + 1 > numel (varargin))
      error ("fitlm: name-value argument '%s' has no value.", varargin{i});
    endif
    val = varargin{i+1};
    switch (key)
      case "weights"
        weights_arg = val;
      case "exclude"
        excl_arg = val;
      case "intercept"
        intercept_arg = logical (val);
      case "categoricalvars"
        cat_vars_arg = val;
      case "varnames"
        var_names_arg = val;
      case "robustopts"
        ## Accepted but not implemented in Phase 4
      case "dummyvarcoding"
        dummy_coding_arg = lower (val);
      case "responsevar"
        resp_var_arg = val;
      case "predictorvars"
        pred_vars_arg = val;
      otherwise
        warning ("fitlm: ignoring unknown name-value argument '%s'.", varargin{i});
    endswitch
    i = i + 2;
  endwhile
endfunction

function idx_set = __resolve_cat_vars__ (cat_vars_arg, pred_names)
  ## Return integer indices (into pred_names) of categorical predictors.
  idx_set = [];
  if (isempty (cat_vars_arg))
    return;
  endif
  p = numel (pred_names);
  if (isnumeric (cat_vars_arg))
    idx_set = round (cat_vars_arg(:))';
  elseif (iscell (cat_vars_arg))
    for k = 1:numel (cat_vars_arg)
      m = find (strcmp (pred_names, cat_vars_arg{k}));
      if (! isempty (m))
        idx_set(end+1) = m;
      endif
    endfor
  elseif (ischar (cat_vars_arg))
    m = find (strcmp (pred_names, cat_vars_arg));
    if (! isempty (m))
      idx_set = m;
    endif
  endif
endfunction
