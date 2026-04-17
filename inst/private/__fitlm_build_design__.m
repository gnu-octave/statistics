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
## @deftypefn {Private Function} {[@var{terms_mat}, @var{coef_names}, @var{X_full}, @
##   @var{y_full}, @var{incl_mask}, @var{p_tot}, @var{has_intercept}, @
##   @var{lp_str}] =} __fitlm_build_design__ (@var{raw_X}, @var{raw_y}, @
##   @var{modelspec}, @var{var_names}, @var{cat_flags}, @var{weights}, @
##   @var{excl_mask}, @var{pred_names}, @var{raw_cols}, @var{has_intercept})
##
## Build design matrix and terms matrix from parsed fitlm inputs.
##
## Returns:
##   terms_mat    - (n_terms x n_all_vars) binary terms matrix (last col = response = 0)
##   coef_names   - 1 x p_tot cell of char coefficient names
##   X_full       - n x p_tot design matrix (ALL rows, NaN for missing rows)
##   y_full       - n x 1 response (NaN for missing rows)
##   incl_mask    - n x 1 logical: rows used in fit
##   p_tot        - total columns in design matrix
##   has_intercept - logical
##   lp_str       - LinearPredictor RHS string
## @end deftypefn

function [terms_mat, coef_names, X_full, y_full, incl_mask, p_tot, has_intercept, lp_str] = ...
    __fitlm_build_design__ (raw_X, raw_y, modelspec, var_names, cat_flags, ...
                            weights, excl_mask, pred_names, raw_cols, has_intercept)

  n = rows (raw_X);
  p = columns (raw_X);           ## number of predictor VARIABLES
  n_vars = p + 1;                ## predictors + response

  ## ── Step 1: Resolve shorthand → Terms matrix ────────────────────────────
  if (isnumeric (modelspec) && ismatrix (modelspec) && ndims (modelspec) == 2)
    ## Numeric terms matrix passed directly (e.g., from stepwiselm)
    terms_pred = modelspec;
    ## If width = p+1 (includes response column), strip last column
    if (columns (terms_pred) == n_vars)
      terms_pred = terms_pred(:, 1:p);
    endif
    ## Detect intercept from the terms matrix
    has_intercept = any (all (terms_pred == 0, 2));
  else
    is_formula = ! isempty (strfind (modelspec, "~"));

    if (is_formula)
      ## Formula string — parse response name and RHS
      tilde = strfind (modelspec, "~");
      rhs = strtrim (modelspec(tilde(1)+1:end));
      ## Check for unsupported C() notation
      if (! isempty (regexp (rhs, 'C\s*\(')))
        error ("Unable to understand the character vector or string scalar '%s'.", modelspec);
      endif
      ## Check for intercept suppression
      if (! isempty (regexp (rhs, '(^|\s|[+])\s*-\s*1\b')))
        has_intercept = false;
        rhs = regexprep (rhs, '\s*[-]\s*1\b', '');
        rhs = strtrim (rhs);
      endif
      ## Parse the RHS terms
      terms_pred = __parse_rhs_terms__ (rhs, pred_names, has_intercept);
    else
      ## Shorthand string
      terms_pred = __shorthand_to_terms__ (lower (modelspec), p, has_intercept);
    endif
  endif


  ## terms_pred: (n_terms x p) — predictor columns only
  ## Expand to include response column (always 0)
  terms_mat = [terms_pred, zeros(rows(terms_pred), 1)];

  ## ── Step 2: Build design matrix column by column ─────────────────────────
  ##
  ## We iterate over term rows.  For each term:
  ##   - intercept row (all 0): add ones column
  ##   - linear/squared numeric: add X^power column
  ##   - categorical main effect: add k-1 dummy columns
  ##   - interaction: element-wise product of factor block columns

  ## Pre-compute per-predictor "base columns" (for numerics) and
  ## "dummy block" (for categoricals)
  base_cols   = cell (1, p);    ## each cell = n x something
  base_cnames = cell (1, p);    ## cell of cell of char
  cat_levels  = cell (1, p);    ## sorted unique levels for each cat predictor

  for k = 1:p
    if (cat_flags(k))
      ## Categorical predictor — build dummy columns (reference = first level)
      col_raw = raw_cols{k};
      if (isa (col_raw, "categorical"))
        levs = cellstr (categories (col_raw));
        [~, cat_idx] = ismember (cellstr (col_raw), levs);
      elseif (iscellstr (col_raw) || isstring (col_raw))
        if (isstring (col_raw))
          col_raw = cellstr (col_raw);
        endif
        levs = unique (col_raw);
        levs = sort (levs);
        [~, ~, cat_idx] = unique (col_raw);
        ## Re-sort so idx matches sorted levs
        [levs2, ~, cat_idx] = unique (col_raw);
        levs = levs2;
      else
        ## Numeric forced categorical
        vals = raw_X(:, k);
        levs_num = unique (vals(! isnan (vals)));
        levs_num = sort (levs_num);
        levs = arrayfun (@num2str, levs_num, "UniformOutput", false);
        cat_idx = zeros (n, 1);
        for li = 1:numel (levs_num)
          cat_idx(vals == levs_num(li)) = li;
        endfor
      endif
      cat_levels{k} = levs;
      n_lev = numel (levs);
      ## Dummy: drop first level (reference = level 1, lowest alphabetically)
      n_dum = n_lev - 1;
      D = zeros (n, n_dum);
      dum_names = cell (1, n_dum);
      for L = 2:n_lev
        D(:, L-1) = double (cat_idx == L);
        if (isnumeric (levs))
          lev_str = num2str (levs(L));
        else
          lev_str = levs{L};
        endif
        dum_names{L-1} = sprintf ("%s_%s", pred_names{k}, lev_str);
      endfor
      base_cols{k}   = D;
      base_cnames{k} = dum_names;
    else
      ## Numeric predictor — single column; power applied per-term
      base_cols{k}   = raw_X(:, k);
      base_cnames{k} = {pred_names{k}};
    endif
  endfor

  ## Now build X_full column by column from terms_pred
  X_full     = [];
  coef_names = {};

  for i = 1:rows (terms_pred)
    row = terms_pred(i, :);

    if (all (row == 0))
      ## Intercept
      X_full     = [X_full, ones(n, 1)];
      coef_names{end+1} = "(Intercept)";
      continue;
    endif

    ## Find active predictors in this term
    active = find (row);
    ## Start with a block of ones × empty name
    blk   = ones (n, 1);
    blk_n = {""};

    for ai = 1:numel (active)
      v   = active(ai);
      pwr = row(v);

      if (cat_flags(v))
        ## Categorical: Cartesian product with dummy block
        D = base_cols{v};      ## n x (n_lev-1)
        dn = base_cnames{v};   ## 1 x (n_lev-1) cell
        new_blk   = [];
        new_blk_n = {};
        for c1 = 1:columns (blk)
          for c2 = 1:columns (D)
            new_blk = [new_blk, blk(:, c1) .* D(:, c2)];
            n1 = blk_n{c1};
            n2 = dn{c2};
            if (isempty (n1))
              new_blk_n{end+1} = n2;
            else
              new_blk_n{end+1} = [n1, ":", n2];
            endif
          endfor
        endfor
        blk   = new_blk;
        blk_n = new_blk_n;
      else
        ## Numeric: multiply by column^power
        col = base_cols{v};
        if (pwr > 1)
          col = col .^ pwr;
          cname = sprintf ("%s^%d", pred_names{v}, pwr);
        else
          cname = pred_names{v};
        endif
        blk = bsxfun (@times, blk, col);
        for c1 = 1:numel (blk_n)
          if (isempty (blk_n{c1}))
            blk_n{c1} = cname;
          else
            blk_n{c1} = [blk_n{c1}, ":", cname];
          endif
        endfor
      endif
    endfor

    X_full     = [X_full, blk];
    coef_names = [coef_names, blk_n];
  endfor

  ## ── Step 3: Compute incl_mask ───────────────────────────────────────────
  ## Missing: any NaN in the NUMERIC data (raw_X numeric cols or raw_y)
  nan_in_X   = any (isnan (raw_X), 2);
  nan_in_y   = isnan (raw_y);
  ## For categorical columns stored as codes, NaN means missing
  missing_mask = nan_in_X | nan_in_y;

  incl_mask = (! excl_mask) & (! missing_mask);
  y_full    = raw_y;

  p_tot = columns (X_full);

  ## ── Step 4: Build LinearPredictor string ────────────────────────────────
  lp_str = __build_lp_str__ (terms_pred, pred_names, has_intercept);

endfunction


## ── Shorthand → terms_pred (n_terms x p) ────────────────────────────────────

function terms = __shorthand_to_terms__ (spec, p, has_intercept)
  ## Returns (n_terms x p) matrix (no response column yet)
  switch (spec)
    case "constant"
      if (has_intercept)
        terms = zeros (1, p);
      else
        terms = zeros (0, p);
      endif

    case "linear"
      if (has_intercept)
        terms = [zeros(1, p); eye(p)];
      else
        terms = eye (p);
      endif

    case "interactions"
      lin_terms = eye (p);
      inter_rows = [];
      for i = 1:p-1
        for j = i+1:p
          r = zeros (1, p);
          r(i) = 1; r(j) = 1;
          inter_rows = [inter_rows; r];
        endfor
      endfor
      if (has_intercept)
        terms = [zeros(1,p); lin_terms; inter_rows];
      else
        terms = [lin_terms; inter_rows];
      endif

    case "purequadratic"
      lin_terms = eye (p);
      quad_rows = 2 * eye (p);
      if (has_intercept)
        terms = [zeros(1,p); lin_terms; quad_rows];
      else
        terms = [lin_terms; quad_rows];
      endif

    case {"quadratic", "full"}
      lin_terms = eye (p);
      inter_rows = [];
      for i = 1:p-1
        for j = i+1:p
          r = zeros (1, p);
          r(i) = 1; r(j) = 1;
          inter_rows = [inter_rows; r];
        endfor
      endfor
      quad_rows = 2 * eye (p);
      if (has_intercept)
        terms = [zeros(1,p); lin_terms; inter_rows; quad_rows];
      else
        terms = [lin_terms; inter_rows; quad_rows];
      endif

    otherwise
      ## Check for polyNM pattern (e.g., 'poly22', 'poly33')
      m = regexp (spec, '^poly(\d+)$', 'tokens');
      if (! isempty (m))
        digits_str = m{1}{1};
        if (numel (digits_str) != p)
          error ("fitlm: 'poly%s' requires %d digits for %d predictors.", ...
                 digits_str, p, p);
        endif
        max_powers = arrayfun (@(c) str2double (c), digits_str);
        terms = __poly_to_terms__ (max_powers, p, has_intercept);
      else
        ## Unknown — default to linear
        if (has_intercept)
          terms = [zeros(1, p); eye(p)];
        else
          terms = eye (p);
        endif
      endif
  endswitch
endfunction


## ── Formula RHS → terms_pred ─────────────────────────────────────────────────

function terms = __parse_rhs_terms__ (rhs, pred_names, has_intercept)
  ## Parse a Wilkinson RHS string with full operator support:
  ## +, -, *, /, :, ^, (group)^N
  ## Returns (n_terms x p) matrix
  p = numel (pred_names);

  ## Split into positive and negative terms (respecting parentheses)
  [add_toks, sub_toks, has_int_in_rhs, remove_int] = __split_additive__ (rhs);

  ## Expand each positive token into term rows
  add_rows = {};
  for i = 1:numel (add_toks)
    expanded = __expand_token__ (add_toks{i}, pred_names, p);
    add_rows = [add_rows, expanded];
  endfor

  ## Remove subtracted terms
  for i = 1:numel (sub_toks)
    expanded = __expand_token__ (sub_toks{i}, pred_names, p);
    for j = 1:numel (expanded)
      add_rows = __remove_row__ (add_rows, expanded{j});
    endfor
  endfor

  ## Deduplicate
  add_rows = __dedup_rows__ (add_rows);

  ## Assemble terms matrix
  n_terms = numel (add_rows);
  terms_body = zeros (n_terms, p);
  for i = 1:n_terms
    terms_body(i, :) = add_rows{i};
  endfor

  ## Handle intercept
  if (remove_int)
    has_intercept = false;
  endif
  if (has_intercept || has_int_in_rhs)
    terms = [zeros(1, p); terms_body];
  else
    terms = terms_body;
  endif
endfunction


## ── LinearPredictor string builder ───────────────────────────────────────────

function s = __build_lp_str__ (terms_pred, pred_names, has_intercept)
  parts = {};
  for i = 1:rows (terms_pred)
    row = terms_pred(i, :);
    if (all (row == 0))
      parts{end+1} = "1";
    else
      active = find (row);
      tparts = {};
      for j = 1:numel (active)
        v = active(j);
        pwr = row(v);
        if (pwr == 1)
          tparts{end+1} = pred_names{v};
        else
          tparts{end+1} = sprintf ("%s^%d", pred_names{v}, pwr);
        endif
      endfor
      parts{end+1} = strjoin (tparts, ":");
    endif
  endfor
  s = strjoin (parts, " + ");
endfunction


## ── Split RHS on + and - (respecting parentheses) ───────────────────────────

function [add_toks, sub_toks, has_int, remove_int] = __split_additive__ (rhs)
  add_toks = {};
  sub_toks = {};
  has_int = false;
  remove_int = false;
  rhs = strtrim (rhs);
  if (isempty (rhs))
    return;
  endif
  depth = 0;
  current = '';
  sign = 1;
  for i = 1:numel (rhs)
    ch = rhs(i);
    if (ch == '(')
      depth++;
      current = [current, ch];
    elseif (ch == ')')
      depth--;
      current = [current, ch];
    elseif (depth == 0 && (ch == '+' || ch == '-'))
      tok = strtrim (current);
      if (! isempty (tok))
        if (strcmp (tok, '1'))
          if (sign == -1), remove_int = true; else, has_int = true; endif
        else
          if (sign == -1), sub_toks{end+1} = tok; else, add_toks{end+1} = tok; endif
        endif
      endif
      current = '';
      if (ch == '+'), sign = 1; else, sign = -1; endif
    else
      current = [current, ch];
    endif
  endfor
  tok = strtrim (current);
  if (! isempty (tok))
    if (strcmp (tok, '1'))
      if (sign == -1), remove_int = true; else, has_int = true; endif
    else
      if (sign == -1), sub_toks{end+1} = tok; else, add_toks{end+1} = tok; endif
    endif
  endif
endfunction


## ── Expand a single token into term rows ─────────────────────────────────────

function rows = __expand_token__ (tok, pred_names, p)
  tok = strtrim (tok);
  rows = {};

  ## Case 1: (group)^N
  m = regexp (tok, '^\((.+)\)\^(\d+)$', 'tokens');
  if (! isempty (m))
    inner_str = m{1}{1};
    power = str2double (m{1}{2});
    inner_toks = strtrim (strsplit (inner_str, '+'));
    inner_vars = [];
    for i = 1:numel (inner_toks)
      idx = find (strcmp (pred_names, strtrim (inner_toks{i})));
      if (! isempty (idx))
        inner_vars(end+1) = idx;
      endif
    endfor
    if (! isempty (inner_vars))
      rows = __group_power_expand__ (inner_vars, power, p);
    endif
    return;
  endif

  ## Case 2: bare (group)
  m = regexp (tok, '^\((.+)\)$', 'tokens');
  if (! isempty (m))
    inner_str = m{1}{1};
    [add_t, sub_t, ~, ~] = __split_additive__ (inner_str);
    for i = 1:numel (add_t)
      tmp_exp = __expand_token__ (add_t{i}, pred_names, p);
      rows = [rows, tmp_exp];
    endfor
    for i = 1:numel (sub_t)
      sub_r = __expand_token__ (sub_t{i}, pred_names, p);
      for j = 1:numel (sub_r)
        rows = __remove_row__ (rows, sub_r{j});
      endfor
    endfor
    rows = __dedup_rows__ (rows);
    return;
  endif

  ## Case 3: crossing (*)
  if (__has_at_toplevel__ (tok, '*'))
    parts = __split_toplevel__ (tok, '*');
    result = __expand_token__ (parts{1}, pred_names, p);
    for i = 2:numel (parts)
      right = __expand_token__ (parts{i}, pred_names, p);
      result = __crossing_rows__ (result, right);
    endfor
    rows = result;
    return;
  endif

  ## Case 4: nesting (/)
  if (__has_at_toplevel__ (tok, '/'))
    parts = __split_toplevel__ (tok, '/');
    left = __expand_token__ (parts{1}, pred_names, p);
    for i = 2:numel (parts)
      right = __expand_token__ (parts{i}, pred_names, p);
      cross = {};
      for li = 1:numel (left)
        for ri = 1:numel (right)
          cross{end+1} = left{li} + right{ri};
        endfor
      endfor
      left = [left, cross];
    endfor
    rows = __dedup_rows__ (left);
    return;
  endif

  ## Case 5: power ladder (var^N)
  m = regexp (tok, '^([a-zA-Z_][a-zA-Z0-9_]*)\^(\d+)$', 'tokens');
  if (! isempty (m))
    varname = m{1}{1};
    power = str2double (m{1}{2});
    idx = find (strcmp (pred_names, varname));
    if (! isempty (idx))
      for pwr = 1:power
        row = zeros (1, p);
        row(idx) = pwr;
        rows{end+1} = row;
      endfor
      return;
    endif
  endif

  ## Case 6: interaction (A:B)
  if (any (tok == ':'))
    parts = strtrim (strsplit (tok, ':'));
    row = zeros (1, p);
    valid = true;
    for i = 1:numel (parts)
      part = parts{i};
      caret = strfind (part, '^');
      if (! isempty (caret))
        base = strtrim (part(1:caret(1)-1));
        pwr = str2double (part(caret(1)+1:end));
      else
        base = part;
        pwr = 1;
      endif
      idx = find (strcmp (pred_names, base));
      if (isempty (idx))
        valid = false;
        break;
      endif
      row(idx(1)) = pwr;
    endfor
    if (valid)
      rows = {row};
    endif
    return;
  endif

  ## Case 7: simple variable
  idx = find (strcmp (pred_names, tok));
  if (! isempty (idx))
    row = zeros (1, p);
    row(idx) = 1;
    rows = {row};
  endif
endfunction


## ── Group power expansion ────────────────────────────────────────────────────

function rows = __group_power_expand__ (var_indices, max_power, p)
  nv = numel (var_indices);
  all_exps = __enum_exponents__ (nv, max_power);
  rows = {};
  for i = 1:size (all_exps, 1)
    row = zeros (1, p);
    for j = 1:nv
      row(var_indices(j)) = all_exps(i, j);
    endfor
    rows{end+1} = row;
  endfor
endfunction


## ── Enumerate exponent vectors ───────────────────────────────────────────────

function E = __enum_exponents__ (nv, max_deg)
  E = [];
  total = (max_deg + 1) ^ nv;
  for i = 0:total-1
    vec = zeros (1, nv);
    val = i;
    for j = 1:nv
      vec(j) = mod (val, max_deg + 1);
      val = floor (val / (max_deg + 1));
    endfor
    s = sum (vec);
    if (s >= 1 && s <= max_deg)
      E = [E; vec];
    endif
  endfor
endfunction


## ── Crossing rows (A*B expansion) ────────────────────────────────────────────

function rows = __crossing_rows__ (left, right)
  cross = {};
  for li = 1:numel (left)
    for ri = 1:numel (right)
      cross{end+1} = left{li} + right{ri};
    endfor
  endfor
  rows = [left, right, cross];
  rows = __dedup_rows__ (rows);
endfunction


## ── Remove first matching row ────────────────────────────────────────────────

function rows = __remove_row__ (rows, target)
  for i = 1:numel (rows)
    if (isequal (rows{i}, target))
      rows(i) = [];
      return;
    endif
  endfor
endfunction


## ── Deduplicate rows ─────────────────────────────────────────────────────────

function rows = __dedup_rows__ (rows)
  keep = true (numel (rows), 1);
  for i = 1:numel (rows)
    if (! keep(i)), continue; endif
    for j = i+1:numel (rows)
      if (keep(j) && isequal (rows{i}, rows{j}))
        keep(j) = false;
      endif
    endfor
  endfor
  rows = rows(keep);
endfunction


## ── Check for char at top level (outside parens) ─────────────────────────────

function result = __has_at_toplevel__ (tok, ch)
  result = false;
  depth = 0;
  for i = 1:numel (tok)
    if (tok(i) == '('), depth++;
    elseif (tok(i) == ')'), depth--;
    elseif (depth == 0 && tok(i) == ch)
      result = true;
      return;
    endif
  endfor
endfunction


## ── Split on char at top level ───────────────────────────────────────────────

function parts = __split_toplevel__ (tok, ch)
  parts = {};
  depth = 0;
  current = '';
  for i = 1:numel (tok)
    if (tok(i) == '('), depth++; current = [current, tok(i)];
    elseif (tok(i) == ')'), depth--; current = [current, tok(i)];
    elseif (depth == 0 && tok(i) == ch)
      parts{end+1} = strtrim (current);
      current = '';
    else
      current = [current, tok(i)];
    endif
  endfor
  if (! isempty (current))
    parts{end+1} = strtrim (current);
  endif
endfunction


## ── Polynomial shorthand to terms ────────────────────────────────────────────

function terms = __poly_to_terms__ (max_powers, p, has_intercept)
  max_deg = max (max_powers);
  exps = cell (1, p);
  for j = 1:p
    exps{j} = 0:max_powers(j);
  endfor
  grids = cell (1, p);
  [grids{:}] = ndgrid (exps{:});
  n_combos = numel (grids{1});
  all_rows = zeros (n_combos, p);
  for j = 1:p
    all_rows(:, j) = grids{j}(:);
  endfor
  total_deg = sum (all_rows, 2);
  keep = (total_deg >= 1) & (total_deg <= max_deg);
  terms_body = all_rows(keep, :);
  [~, ord] = sortrows ([sum(terms_body, 2), terms_body]);
  terms_body = terms_body(ord, :);
  if (has_intercept)
    terms = [zeros(1, p); terms_body];
  else
    terms = terms_body;
  endif
endfunction
