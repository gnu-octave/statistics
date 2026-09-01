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
## @deftypefn {Private Function} {@var{W} =} edgeWeights (@var{args}, @var{Y}, @var{ClassNames}, @var{Prior}, @var{classname}, @var{caller})
##
## Normalized observation weights of an @code{edge} call.
##
## @var{args} is the @code{varargin} of the caller, which accepts a single
## @qcode{'Weights'} Name-Value pair; without one the weights are uniform.
##
## The weights are then normalized @strong{within each class to that class's
## prior}, which is what MATLAB does and is not the same as dividing by their
## total.  It matters whenever the weights are not uniform: on a weighted iris
## edge the two readings are 0.9269539697 and 0.9438468986, and the latter is
## the oracle's.  Normalizing this way keeps a class's total influence equal to
## its prior however the weights inside it are distributed, so reweighting
## observations cannot quietly reweight the classes.
##
## @var{classname} and @var{caller} name the class and the method in the error
## messages, so a class reports its own name rather than this helper's.
##
## @end deftypefn

function W = edgeWeights (args, Y, ClassNames, Prior, classname, caller)

  n = rows (Y);
  W = ones (n, 1);
  for i = 1:2:numel (args)
    if (! (ischar (args{i}) && isrow (args{i})))
      error ("%s.%s: parameter name must be a character vector.", ...
             classname, caller);
    endif
    if (strcmpi (args{i}, 'weights'))
      W = args{i+1};
      if (! (isnumeric (W) && isvector (W)))
        error ("%s.%s: 'Weights' must be a numeric vector.", ...
               classname, caller);
      endif
      if (numel (W) != n)
        error (strcat ("%s.%s: size of 'Weights' must equal the number", ...
                       " of rows in X."), classname, caller);
      endif
    else
      error (strcat ("%s.%s: invalid parameter name in optional paired", ...
                     " arguments."), classname, caller);
    endif
  endfor
  W = double (W(:));

  ## Which class each observation belongs to.  Matching against ClassNames is
  ## the general route; a response already in some other coding, such as the
  ## +1/-1 a support vector machine accepts, falls back to its own grouping.
  gY = classindex (Y, ClassNames);
  if (isempty (gY))
    gY = grp2idx (Y);
    Prior = ones (1, max (gY)) / max (gY);
  endif

  for k = 1:numel (Prior)
    ck = (gY == k);
    if (any (ck) && sum (W(ck)) > 0)
      W(ck) = W(ck) / sum (W(ck)) * Prior(k);
    endif
  endfor

endfunction

## Class index of every element of Y against ClassNames, empty when any
## element fails to match.
function gY = classindex (Y, ClassNames)

  gY = [];
  n = rows (Y);
  idx = cell (n, 1);
  if (iscellstr (ClassNames))
    if (! (iscellstr (Y) || ischar (Y)))
      return;
    endif
    Yc = cellstr (Y);
    for i = 1:n
      idx{i} = find (strcmp (Yc{i}, ClassNames), 1);
    endfor
  elseif (ischar (ClassNames))
    CN = cellstr (ClassNames);
    if (! (iscellstr (Y) || ischar (Y)))
      return;
    endif
    Yc = cellstr (Y);
    for i = 1:n
      idx{i} = find (strcmp (Yc{i}, CN), 1);
    endfor
  else
    if (! (isnumeric (Y) || islogical (Y)))
      return;
    endif
    for i = 1:n
      idx{i} = find (Y(i) == ClassNames(:), 1);
    endfor
  endif
  if (any (cellfun (@isempty, idx)))
    return;
  endif
  gY = cell2mat (idx);

endfunction
