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
## @deftypefn {Private Function} {@var{Formula} =} __fitlm_build_formula__ (@
##   @var{terms_mat}, @var{var_names}, @var{pred_names}, @var{resp_name}, @
##   @var{has_intercept}, @var{lp_str}, @var{coef_names})
##
## Build the Formula struct for a LinearModel.
##
## @var{terms_mat} is the (n_terms x n_all_vars) binary matrix including
## the response column (always 0).
##
## Returns a struct with the 14 confirmed fields from MATLAB Block 1-P4.
## @end deftypefn

function Formula = __fitlm_build_formula__ (terms_mat, var_names, pred_names, ...
                                             resp_name, has_intercept, lp_str, coef_names)

  n_terms = rows (terms_mat);
  n_vars  = columns (terms_mat);     ## includes response column
  n_pred  = numel (pred_names);

  ## Build TermNames from terms_mat (using pred columns only)
  terms_pred = terms_mat(:, 1:n_pred);
  term_names = cell (n_terms, 1);
  for i = 1:n_terms
    row = terms_pred(i, :);
    if (all (row == 0))
      term_names{i} = "(Intercept)";
    else
      active = find (row);
      parts = {};
      for j = 1:numel (active)
        v = active(j);
        pwr = row(v);
        if (pwr == 1)
          parts{end+1} = pred_names{v};
        else
          parts{end+1} = sprintf ("%s^%d", pred_names{v}, pwr);
        endif
      endfor
      term_names{i} = strjoin (parts, ":");
    endif
  endfor

  ## InModel: logical vector over var_names (1 = predictor in model, 0 = response)
  ## Length = n_vars.  Response column is always 0.  A predictor is InModel if
  ## any non-intercept term has a non-zero entry in that predictor's column.
  in_model = false (1, n_vars);
  for j = 1:n_pred
    if (any (terms_pred(:, j) != 0))
      in_model(j) = true;
    endif
  endfor
  ## Response column (last) = 0 (already false)

  Formula.LinearPredictor  = lp_str;
  Formula.ResponseName     = resp_name;
  Formula.Terms            = terms_mat;
  Formula.VariableNames    = var_names;     ## e.g. {'x1','x2','x3','y'}
  Formula.HasIntercept     = logical (has_intercept);
  Formula.InModel          = in_model;
  Formula.TermNames        = term_names;
  Formula.PredictorNames   = pred_names;
  Formula.NTerms           = n_terms;
  Formula.NVars            = n_vars;
  Formula.NPredictors      = n_pred;
  Formula.Link             = "identity";
  Formula.ModelFun         = @(b, X) X * b;
  Formula.FunctionCalls    = {};

endfunction
