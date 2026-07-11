## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Private Function} {[@var{X_enc}, @var{enc_names}, @var{cat_info}] =} encode_categorical (@var{X_num}, @var{cat_cols}, @var{pred_names}, @var{cat_levels})
##
## Encode categorical predictors as reference-coded indicator (dummy) variables.
##
## @var{X_num} is the @math{n}-by-@math{p} numeric predictor matrix (categorical
## columns hold 1-based level codes), @var{cat_cols} is a logical vector marking
## the categorical columns, @var{pred_names} is the cell array of predictor
## names, and @var{cat_levels} is a cell array whose @math{j}-th entry lists the
## level labels of predictor @math{j} (empty to infer them from the data).
##
## Numeric predictors are copied through unchanged.  Each categorical predictor
## with @math{k} levels expands to @math{k - 1} indicator columns (the first
## level is the omitted reference), named @qcode{@var{name}_@var{level}}.
##
## @var{X_enc} is the encoded design matrix, @var{enc_names} the cell array of
## its column names, and @var{cat_info} a structure with fields @qcode{names}
## and @qcode{levels} recording the categorical predictors and the level labels
## used, so the same encoding can be reproduced for new data.
##
## This helper is shared by the @code{LinearModel} and
## @code{GeneralizedLinearModel} classes.
##
## @end deftypefn

function [X_enc, enc_names, cat_info] = encode_categorical ( ...
    X_num, cat_cols, pred_names, cat_levels)

  X_enc     = zeros (rows (X_num), 0);
  enc_names = {};
  cat_info.names  = {};
  cat_info.levels = {};

  for j = 1:numel (pred_names)
    if (! cat_cols(j))
      X_enc     = [X_enc, X_num(:, j)];
      enc_names = [enc_names, pred_names{j}];
    else
      levels_j = cat_levels{j};
      if (isempty (levels_j))
        uvals    = sort (unique (X_num(isfinite (X_num(:,j)), j)));
        levels_j = cellstr (num2str (uvals(:)));
      endif
      n_lev = numel (levels_j);
      for L = 2:n_lev
        dummy     = double (X_num(:, j) == L);
        X_enc     = [X_enc, dummy];
        enc_names = [enc_names, [pred_names{j}, '_', char(levels_j{L})]];
      endfor
      cat_info.names{end+1}  = pred_names{j};
      cat_info.levels{end+1} = levels_j;
    endif
  endfor
endfunction
