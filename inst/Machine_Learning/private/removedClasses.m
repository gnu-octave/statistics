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
## @deftypefn {Private Function} {[@var{drop}, @var{Prior}, @var{msg}] =} removedClasses (@var{Prior}, @var{Cost}, @var{ClassNames}, @var{gY})
##
## The classes a classifier removes from its data before fitting.
##
## A class whose prior probability is 0, or whose row of misclassification
## costs is all zeros, can never be predicted and contributes nothing to the
## fit, its observations weighing nothing, so its observations are removed,
## as MATLAB removes them.  @var{Prior} is the prior as a vector in the order
## of @var{ClassNames}, and @var{Cost} the cost matrix, or empty.  @var{drop}
## is a logical row, true for each class removed; @var{Prior} is returned
## with those classes at 0, normalized to sum to 1.  @var{gY} holds the class
## index of each observation; @var{msg} names the removed classes that hold
## any, the text of the warning the caller raises under its own name, and is
## empty when none does, a class absent from the data having nothing to
## remove.
## @end deftypefn

function [drop, Prior, msg] = removedClasses (Prior, Cost, ClassNames, gY)

  drop = (Prior(:)' == 0);
  ## A row of zero costs means a class free to misclassify, which needs
  ## another class to be misclassified as: a model of one class has none
  if (! isempty (Cost) && rows (Cost) > 1)
    drop |= all (Cost == 0, 2)';
  endif
  msg = '';
  if (! any (drop))
    return;
  endif
  Prior(drop) = 0;
  if (sum (Prior) > 0)
    Prior = Prior / sum (Prior);
  endif

  ## The removed classes named as text, whatever the type of ClassNames
  if (iscellstr (ClassNames))
    names = ClassNames(:);
  elseif (ischar (ClassNames) || iscategorical (ClassNames)
          || isstring (ClassNames))
    names = cellstr (ClassNames);
  else
    names = arrayfun (@num2str, ClassNames(:), 'UniformOutput', false);
  endif
  held = false (size (drop));
  held(unique (gY(isfinite (gY) & gY > 0))) = true;
  names = names(drop & held(1:numel (drop)));
  if (isempty (names))
    return;
  endif
  listed = strjoin (strcat ("'", names(:)', "'"), ', ');
  if (numel (names) == 1)
    msg = sprintf (strcat ("the observations of class %s are removed, its", ...
                           " prior probability being 0."), listed);
  else
    msg = sprintf (strcat ("the observations of classes %s are removed,", ...
                           " their prior probabilities being 0."), listed);
  endif

endfunction
