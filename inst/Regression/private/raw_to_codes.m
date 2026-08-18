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
## @deftypefn {Private Function} {[@var{X_num}, @var{cat_levels}] =} raw_to_codes (@var{data}, @var{X_raw}, @var{tbl}, @var{pred_names}, @var{cat_logical}, @var{n_total})
##
## Convert raw predictor data to numeric level codes for categorical columns.
##
## @var{data} is the predictor data as the caller received it, and decides which
## branch is taken: a table is read variable by variable through @var{tbl},
## anything else column by column through the numeric matrix @var{X_raw}.
## @var{pred_names} names the predictors, @var{cat_logical} marks the ones to
## treat as categorical, and @var{n_total} is the number of observations.
##
## @var{X_num} holds the numeric predictors unchanged and 1-based level codes
## for the categorical ones, and @var{cat_levels} lists the level labels of each
## predictor, empty for a numeric one.  The output feeds
## @code{encode_categorical}, which turns the codes into indicator columns.
##
## This helper is shared by the @code{CoxModel} and
## @code{GeneralizedLinearModel} classes.
##
## @end deftypefn

function [X_num, cat_levels] = raw_to_codes (data, X_raw, tbl, pred_names, ...
                                             cat_logical, n_total)
  p = numel (pred_names);
  X_num      = zeros (n_total, p);
  cat_levels = cell (1, p);
  for j = 1:p
    if (istable (data))
      col = tbl.(pred_names{j});
      if (iscell (col))
        ## Appearance order, so that the omitted reference level is the one the
        ## data shows first; see parseWilkinsonFormula for the formula path.
        [cat_levels{j}, ~, ic] = unique (col, 'stable');
        X_num(:, j) = ic;
      elseif (isa (col, 'categorical'))
        cat_levels{j} = categories (col);
        [~, ic] = ismember (cellstr (col), cat_levels{j});
        X_num(:, j) = ic;
      else
        X_num(:, j) = double (col(:));
        cat_levels{j} = {};
      endif
    else
      if (cat_logical(j))
        uvals = sort (unique (X_raw(isfinite (X_raw(:,j)), j)));
        cat_levels{j} = strtrim (cellstr (num2str (uvals(:))));
        [~, ic] = ismember (X_raw(:,j), uvals);
        X_num(:, j) = ic;
      else
        X_num(:, j) = X_raw(:, j);
        cat_levels{j} = {};
      endif
    endif
  endfor
endfunction
