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
  is_formula = ! isempty (strfind (modelspec, "~"));

  if (is_formula)
    ## Formula string — parse response name and RHS
    tilde = strfind (modelspec, "~");
    rhs = strtrim (modelspec(tilde(1)+1:end));
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
        [~, cat_idx] = ismember (col_raw, levs);
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
      ## Unknown — default to linear
      if (has_intercept)
        terms = [zeros(1, p); eye(p)];
      else
        terms = eye (p);
      endif
  endswitch
endfunction


## ── Formula RHS → terms_pred ─────────────────────────────────────────────────

function terms = __parse_rhs_terms__ (rhs, pred_names, has_intercept)
  ## Parse a Wilkinson RHS string like '1 + x1 + x2 + x1:x2'
  ## Returns (n_terms x p) binary matrix
  p = numel (pred_names);

  ## Split on '+' to get individual term tokens
  tokens = strtrim (strsplit (rhs, "+"));

  rows_list = {};
  has_int_in_rhs = false;

  for ti = 1:numel (tokens)
    tok = strtrim (tokens{ti});
    if (isempty (tok))
      continue;
    endif
    if (strcmp (tok, "1"))
      has_int_in_rhs = true;
      continue;
    endif
    ## Split on ':' for interaction
    parts = strtrim (strsplit (tok, ":"));
    row = zeros (1, p);
    valid = true;
    for pi = 1:numel (parts)
      part = parts{pi};
      ## Handle power: x1^2
      caret = strfind (part, "^");
      if (! isempty (caret))
        base = strtrim (part(1:caret(1)-1));
        pwr  = str2double (part(caret(1)+1:end));
      else
        base = part;
        pwr  = 1;
      endif
      idx = find (strcmp (pred_names, base));
      if (isempty (idx))
        ## Could be response name or unknown — skip
        valid = false;
        break;
      endif
      row(idx(1)) = pwr;
    endfor
    if (valid)
      rows_list{end+1} = row;
    endif
  endfor

  ## Assemble
  n_terms = numel (rows_list);
  terms_body = zeros (n_terms, p);
  for i = 1:n_terms
    terms_body(i, :) = rows_list{i};
  endfor

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
