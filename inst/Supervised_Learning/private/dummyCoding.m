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
## vector with one element per predictor, or as @qcode{'all'}, and
## @var{names} holds the predictor names.  @var{C} is a structure with the
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
  else
    errmsg = strcat ("'CategoricalPredictors' must be a vector of positive", ...
                     " integers, a logical vector or 'all'.");
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
