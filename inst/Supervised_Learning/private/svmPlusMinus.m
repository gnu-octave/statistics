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
## @deftypefn {statistics} {@var{Ypm} =} svmPlusMinus (@var{Y}, @var{ClassNames})
##
## Response of a two-class support vector machine in @math{+1/-1} coding.
##
## @var{Y} may hold the class labels, in whatever type @var{ClassNames} holds
## them, which is what @code{margin} and @code{loss} document and what MATLAB
## accepts.  It may also already be numeric @math{+1/-1}, which is the coding
## the underlying solver works in, and is then passed through unchanged.
##
## The first class in @var{ClassNames} maps to @math{+1} and the second to
## @math{-1}, matching the sign convention the decision values carry.
##
## @end deftypefn

function Ypm = svmPlusMinus (Y, ClassNames)

  ## Already the solver's own coding, so nothing to map.  This is checked
  ## first and without consulting ClassNames on purpose: where the classes
  ## themselves are -1 and +1 the two readings coincide, and mapping through
  ## the sorted names would invert the sign of every margin.
  if (isnumeric (Y) && ! isempty (Y) && all (ismember (Y(:), [-1; 1])))
    Ypm = double (Y(:));
    return;
  endif

  if (iscellstr (ClassNames))
    idx = cellfun (@(v) find (strcmp (v, ClassNames), 1), cellstr (Y), ...
                   'UniformOutput', false);
  elseif (ischar (ClassNames))
    idx = arrayfun (@(k) find (strcmp (Y(k,:), cellstr (ClassNames)), 1), ...
                    (1:rows (Y))', 'UniformOutput', false);
  else
    idx = arrayfun (@(v) find (v == ClassNames(:), 1), Y(:), ...
                    'UniformOutput', false);
  endif

  if (any (cellfun (@isempty, idx)))
    error ("svmPlusMinus: Y contains a label that is not in ClassNames.");
  endif
  idx = cell2mat (idx(:));
  Ypm = ones (numel (idx), 1);
  Ypm(idx == 2) = -1;

endfunction
