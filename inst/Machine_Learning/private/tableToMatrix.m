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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {private} {@var{Z} =} tableToMatrix (@var{pnames}, @var{lev}, @var{T}, @var{caller})
##
## Map a table onto a named set of predictors, by name and not by position.
##
## @var{pnames} names the predictors wanted and @var{lev} holds, one cell
## each, the levels a predictor read from text was coded through, empty for
## one read from numbers.  A variable @var{pnames} does not name is passed
## over, one it names and the table does not hold is named in the error, and
## a value holding a level is coded as that level was coded before, so a
## table carrying only some of the levels still codes them the same way.
##
## A matrix passes through, its columns taken in the order they come.
##
## A fitted model reaches this through @code{PredictiveModel.tableColumns}.
## @code{shapley} and @code{lime} reach it directly: each holds a model
## without being one, so neither inherits that method, and the mapping is
## written once rather than three times.
##
## @end deftypefn

function Z = tableToMatrix (pnames, lev, T, caller)

  if (! istable (T))
    Z = T;
    return;
  endif
  if (ischar (pnames))
    pnames = cellstr (pnames);
  endif
  have = T.Properties.VariableNames;
  Z = zeros (height (T), numel (pnames));
  for k = 1:numel (pnames)
    j = find (strcmp (have, pnames{k}), 1);
    if (isempty (j))
      error ("%s: the table holds no predictor '%s'.", caller, pnames{k});
    endif
    col = T.(pnames{k});
    if (numel (lev) >= k && ! isempty (lev{k}))
      Z(:,k) = ttmLevelCodes (col, lev{k}, caller, pnames{k});
    else
      if (! (isnumeric (col) && isreal (col)) && ! islogical (col))
        error (strcat ("%s: the table variable '%s' no longer holds what", ...
                       " it held when the model was fitted."), caller, ...
               pnames{k});
      endif
      Z(:,k) = double (col);
    endif
  endfor

endfunction

## The codes a column of levels carried when the model was fitted, so that a
## table holding only some of them still codes them the same way.
function v = ttmLevelCodes (col, lev, caller, name)

  if (isa (col, 'categorical'))
    col = cellstr (col);
  elseif (ischar (col) || isa (col, 'string'))
    col = cellstr (col);
  endif
  if (! iscellstr (col))
    error (strcat ("%s: the table variable '%s' no longer holds what it", ...
                   " held when the model was fitted."), caller, name);
  endif
  v = zeros (numel (col), 1);
  for k = 1:numel (lev)
    v(strcmp (col(:), lev{k})) = k;
  endfor
  if (any (v == 0))
    error (strcat ("%s: the table variable '%s' holds a level the model", ...
                   " was not fitted on."), caller, name);
  endif

endfunction
