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
## @deftypefn {Private Function} {[@var{terms}, @var{cat_info}, @var{col_names}] =} terms_from_coefnames (@var{coef_names}, @var{pred_names}, @var{cat_logical}, @var{data}, @var{tbl_sub})
##
## Recover a terms matrix and the categorical level information from the
## coefficient names a model formula produced.
##
## @var{coef_names} lists the coefficient names, @var{pred_names} the model's
## predictors, @var{cat_logical} marks the categorical ones, @var{data} is the
## data as given, and @var{tbl_sub} the table of rows kept for the fit.
##
## @var{terms} is the terms matrix over the encoded design columns,
## @var{cat_info} a structure with fields @qcode{names} and @qcode{levels}, and
## @var{col_names} the name of each column of @var{terms}.  @var{col_names} is
## @emph{not} @var{coef_names}: a factor appearing only inside an interaction or
## a power gets a column of its own without ever being a coefficient, so the two
## lists diverge and anything indexed by the columns of @var{terms} must be
## built from @var{col_names}.  The
## levels of a character or string grouping column are taken in the order the
## data presents them, matching the order the design matrix was built in; a
## @code{categorical} column carries its own category order.
##
## This helper is shared by the @code{LinearModel} and
## @code{GeneralizedLinearModel} classes.
##
## @end deftypefn

function [terms, cat_info, col_names] = terms_from_coefnames (coef_names, pred_names, ...
                                                   cat_logical, data, tbl_sub)
  cat_info.names  = {};
  cat_info.levels = {};
  if (istable (data))
    for j = 1:numel (pred_names)
      if (cat_logical(j))
        col = tbl_sub.(pred_names{j});
        if (iscell (col) || isa (col, 'string'))
          ## In the order the data presents them, as the design matrix was
          ## built: these levels re-encode new data at prediction time, so a
          ## sorted list here would silently map it onto the wrong indicators.
          levels_j = cellstr (unique (col(:), 'stable'));
        elseif (isa (col, 'categorical'))
          levels_j = categories (col);
        elseif (islogical (col) || isnumeric (col))
          ## Logical, or numeric and declared categorical: coded by ascending
          ## value, so the reference level is the smallest.
          vals = double (col(:));
          vals = sort (unique (vals(! isnan (vals))));
          levels_j = arrayfun (@(x) strtrim (num2str (x)), vals, ...
                               'UniformOutput', false);
        else
          levels_j = {};
        endif
        cat_info.names{end+1}  = pred_names{j};
        cat_info.levels{end+1} = levels_j;
      endif
    endfor
  endif

  nc = numel (coef_names);
  atomic = {};
  for t = 1:nc
    if (strcmp (coef_names{t}, '(Intercept)'))
      continue;
    endif
    for f = strsplit (coef_names{t}, ':')
      if (! any (strcmp (atomic, f{1})))
        atomic{end+1} = f{1};
      endif
    endfor
  endfor
  terms = zeros (nc, numel (atomic) + 1);
  for t = 1:nc
    if (strcmp (coef_names{t}, '(Intercept)'))
      continue;
    endif
    for f = strsplit (coef_names{t}, ':')
      terms(t, strcmp (atomic, f{1})) = 1;
    endfor
  endfor
  col_names = [atomic, {''}];
endfunction
