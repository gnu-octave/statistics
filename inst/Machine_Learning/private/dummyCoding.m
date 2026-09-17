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
## @deftypefn  {Private Function} {[@var{C}, @var{errmsg}] =} dummyCoding (@var{X}, @var{spec}, @var{names})
## @deftypefnx {Private Function} {@var{XD} =} dummyCoding (@var{X}, @var{C})
##
## Dummy code the categorical predictors of a model.
##
## The first form learns the coding from the training predictors @var{X}.
## @var{spec} names the categorical predictors as indices, as a logical
## vector with one element per predictor, as @qcode{'all'}, or by name as a
## character matrix of one padded name per row, a string array or a cell
## array of character vectors, and @var{names} holds the predictor names.
## A name is matched against @var{names} exactly, so its case must agree;
## the padding blanks of a character matrix are stripped before matching
## and nothing else is.  @qcode{'all'} is read as every predictor even where
## a predictor is named @qcode{'all'}.
##
## Three readings of a name depart from MATLAB R2024a, each deliberately.
## Repeated names give one index where MATLAB repeats it; the indices are a
## row for every form of @var{spec}, where MATLAB returns a column for a
## character matrix and a row for the rest; and a name is read against the
## default names as readily as against given ones, where MATLAB requires
## @qcode{'PredictorNames'} before it will read any name at all.
##
## @var{C} is a structure with the
## fields @qcode{Index}, the sorted indices of the categorical predictors;
## @qcode{Levels}, a cell array with the sorted levels of each predictor,
## empty for a numeric one; @qcode{NumPredictors}; @qcode{ExpandedNames},
## the names of the coded columns; and @qcode{Dummy}, a logical row marking
## the coded columns that hold a level.  A problem with
## @var{spec} leaves @var{C} empty and returns the message in @var{errmsg}.
##
## The second form codes @var{X} with @var{C}.  Every predictor keeps its
## place: a numeric one is copied, and a categorical one becomes a column of
## zeros and ones per level, named @qcode{"x1 == 2"}, as MATLAB R2024a codes
## them.  A value that is not one of the levels, a missing value included,
## makes every column of that predictor @code{NaN}.
##
## @end deftypefn

function [C, errmsg] = dummyCoding (X, spec, names)

  if (isstruct (spec))
    C = code (X, spec);
    return;
  endif

  C = [];
  errmsg = '';
  p = columns (X);
  if (ischar (spec) && strcmpi (spec, 'all'))
    idx = 1:p;
  elseif (islogical (spec) && (isvector (spec) || isempty (spec)))
    if (numel (spec) != p)
      errmsg = strcat ("a logical 'CategoricalPredictors' must have one", ...
                       " element per predictor.");
      return;
    endif
    idx = find (spec(:)');
  elseif (isnumeric (spec) && (isvector (spec) || isempty (spec))
          && all (spec == fix (spec)) && all (spec > 0))
    if (any (spec > p))
      errmsg = strcat ("'CategoricalPredictors' indices must not exceed", ...
                       " the number of predictors.");
      return;
    endif
    idx = unique (double (spec(:)'));
  elseif (ischar (spec) || iscellstr (spec) || isa (spec, 'string'))
    ## A character matrix gives one name per row, its padding blanks stripped;
    ## a cell array or string array gives one name per element, kept whole.
    if (ischar (spec))
      wanted = cellstr (spec)';
    elseif (iscellstr (spec))
      wanted = spec(:)';
    else
      wanted = cellstr (spec(:)');
    endif
    if (isempty (wanted))
      idx = [];
    elseif (numel (names) != p)
      errmsg = strcat ("'CategoricalPredictors' can name a predictor only", ...
                       " where the predictor names are known.");
      return;
    else
      idx = zeros (1, 0);
      for k = 1:numel (wanted)
        j = find (strcmp (names, wanted{k}));
        if (isempty (j))
          errmsg = sprintf (strcat ("'CategoricalPredictors' does not name", ...
                                    " a predictor: '%s'"), wanted{k});
          return;
        endif
        idx(end+1) = j(1);
      endfor
      idx = unique (idx);
    endif
  else
    errmsg = strcat ("'CategoricalPredictors' must be a vector of positive", ...
                     " integers, a logical vector, a character matrix, a", ...
                     " string array, a cell array of character vectors or", ...
                     " 'all'.");
    return;
  endif

  C = struct ('Index', idx, 'NumPredictors', p);
  C.Levels = cell (1, p);
  C.ExpandedNames = {};
  C.Dummy = false (1, 0);
  for j = 1:p
    if (any (idx == j))
      v = X(:,j);
      C.Levels{j} = unique (v(! isnan (v)))';
      C.Dummy = [C.Dummy, true(1, numel (C.Levels{j}))];
      C.ExpandedNames = [C.ExpandedNames, arrayfun(@(l) sprintf ('%s == %s', ...
                                            names{j}, num2str (l)), ...
                                            C.Levels{j}, ...
                                            'UniformOutput', false)];
    else
      C.ExpandedNames{end+1} = names{j};
      C.Dummy(end+1) = false;
    endif
  endfor

endfunction

function XD = code (X, C)

  if (isempty (C.Index))
    XD = X;
    return;
  endif
  XD = zeros (rows (X), 0);
  for j = 1:C.NumPredictors
    lev = C.Levels{j};
    if (isempty (lev) && ! any (C.Index == j))
      XD = [XD, X(:,j)];
    else
      D = double (X(:,j) == lev);
      D(! any (D, 2),:) = NaN;
      XD = [XD, D];
    endif
  endfor

endfunction
