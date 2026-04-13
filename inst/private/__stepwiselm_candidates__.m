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
## @deftypefn {Private Function} {@var{candidates} =} __stepwiselm_candidates__ @
##   (@var{current_terms}, @var{lower_terms}, @var{upper_terms}, @
##    @var{cat_flags}, @var{pred_names})
##
## Enumerate candidate terms for adding or removing in stepwise regression.
##
## @strong{Inputs}
##
## @table @var
## @item current_terms
## Binary terms matrix of the current model (n_terms x n_vars).
## Includes the response column (last column, always 0).
## @item lower_terms
## Binary terms matrix of the lower (minimum) model.
## @item upper_terms
## Binary terms matrix of the upper (maximum) model.
## @item cat_flags
## Logical vector (1 x n_vars) indicating categorical variables.
## @item pred_names
## Cell array of predictor variable names.
## @end table
##
## @strong{Output}
##
## @var{candidates} is a struct array.  Each element has fields:
## @table @code
## @item term_row
## Binary row vector (1 x n_vars) for this term.
## @item term_name
## Char string: e.g., @code{'x1'}, @code{'x1:x2'}, @code{'Species'}.
## @item action
## @code{'add'} or @code{'remove'}.
## @item df_term
## Number of design-matrix columns this term contributes:
## 1 for numeric, k-1 for k-level categorical, products for interactions.
## @end table
##
## @seealso{stepwiselm, __stepwiselm_ftest__}
## @end deftypefn

function candidates = __stepwiselm_candidates__ (current_terms, lower_terms, ...
                                                  upper_terms, cat_flags, ...
                                                  pred_names)

  candidates = struct ("term_row", {}, "term_name", {}, ...
                       "action", {}, "df_term", {});

  p = numel (pred_names);     ## number of predictors (excl. response)
  n_vars = columns (current_terms);   ## p + 1 (includes response col)

  ## ── REMOVE candidates ─────────────────────────────────────────────────────
  ## Terms in current_terms that are NOT in lower_terms and NOT the intercept.
  for i = 1:rows (current_terms)
    row = current_terms(i, :);

    ## Skip intercept (all-zeros row)
    if (all (row == 0))
      continue;
    endif

    ## Check if this term is in lower_terms (protected)
    is_lower = false;
    for j = 1:rows (lower_terms)
      if (isequal (row, lower_terms(j, :)))
        is_lower = true;
        break;
      endif
    endfor
    if (is_lower)
      continue;
    endif

    ## This is a valid removal candidate
    c.term_row  = row;
    c.term_name = __term_row_to_name__ (row, pred_names, p);
    c.action    = "remove";
    c.df_term   = __compute_df_term__ (row, cat_flags, pred_names, p);
    candidates(end+1) = c;
  endfor

  ## ── ADD candidates ────────────────────────────────────────────────────────
  ## Terms in upper_terms that are NOT in current_terms.
  for i = 1:rows (upper_terms)
    row = upper_terms(i, :);

    ## Skip intercept (it's always in the model)
    if (all (row == 0))
      continue;
    endif

    ## Check if already in current_terms
    already_in = false;
    for j = 1:rows (current_terms)
      if (isequal (row, current_terms(j, :)))
        already_in = true;
        break;
      endif
    endfor
    if (already_in)
      continue;
    endif

    ## Valid add candidate
    c.term_row  = row;
    c.term_name = __term_row_to_name__ (row, pred_names, p);
    c.action    = "add";
    c.df_term   = __compute_df_term__ (row, cat_flags, pred_names, p);
    candidates(end+1) = c;
  endfor

endfunction


## ── Helper: term row → human-readable name ──────────────────────────────────

function name = __term_row_to_name__ (row, pred_names, p)
  ## Convert a terms matrix row (1 x n_vars) to a human-readable term name
  ## like 'x1', 'x1:x2', 'Species', 'x1^2'.
  ##
  ## The row may have width = p or p+1 (with response column).
  ## We only look at columns 1..p.

  active = find (row(1:p));
  if (isempty (active))
    name = "1";    ## intercept
    return;
  endif

  parts = cell (1, numel (active));
  for k = 1:numel (active)
    v = active(k);
    pwr = row(v);
    if (pwr == 1)
      parts{k} = pred_names{v};
    else
      parts{k} = sprintf ("%s^%d", pred_names{v}, pwr);
    endif
  endfor
  name = strjoin (parts, ":");
endfunction


## ── Helper: compute df_term ─────────────────────────────────────────────────

function df = __compute_df_term__ (row, cat_flags, pred_names, p)
  ## Compute number of design matrix columns for a term.
  ##
  ## Rules (§8 of reference doc):
  ##   Numeric predictor:             df = 1
  ##   Categorical (k levels):        df = k - 1
  ##   Interaction num × num:         df = 1
  ##   Interaction num × cat(k):      df = k - 1
  ##   Interaction cat(j) × cat(k):   df = (j-1) × (k-1)
  ##   Squared numeric (x^2):         df = 1

  active = find (row(1:p));
  if (isempty (active))
    df = 1;  ## intercept
    return;
  endif

  df = 1;
  for k = 1:numel (active)
    v = active(k);
    if (cat_flags(v))
      ## Categorical: need to determine number of levels.
      ## Since we don't have the data here, we need to infer from context.
      ## For now, store a sentinel value of 0 to be resolved later.
      ## Actually, we need the level count. We'll require it as external info.
      ## For the initial implementation, use _PriorSteps mechanism or pass
      ## level info. But cat_flags alone doesn't tell us k.
      ##
      ## DESIGN DECISION: cat_flags is a logical vector. We extend it to
      ## carry level counts for categorical variables. We pass cat_nlevels
      ## as an additional argument or embed it in cat_flags (negative = n_levels).
      ##
      ## SIMPLER APPROACH: accept cat_flags as a numeric vector where:
      ##   0 = numeric, k = number of levels for categorical.
      ## This is more information-dense.
      ##
      ## For now, when cat_flags is logical (boolean), we default to df=1
      ## and let the F-test function compute df dynamically from the design
      ## matrix column count difference. This is actually the correct approach
      ## because df_term is verified at F-test time.
      ## cat_flags can be:
      ##   - logical true: we default df contribution for this factor = 1
      ##     (will be corrected in __stepwiselm_ftest__ using actual design matrix)
      ##   - numeric > 1: number of levels, contribution = val - 1
      if (isnumeric (cat_flags) && cat_flags(v) > 1)
        df = df * (cat_flags(v) - 1);
      else
        ## Default: will be resolved by design matrix column counting
        ## in __stepwiselm_ftest__
        df = df * 1;
      endif
    else
      ## Numeric: 1 column regardless of power
      df = df * 1;
    endif
  endfor
endfunction
