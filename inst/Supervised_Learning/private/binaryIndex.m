## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{gY} =} binaryIndex (@var{ClassNames}, @var{Y}, @var{classname}, @var{caller})
##
## The index into @var{ClassNames} of every element of @var{Y}.
##
## @var{Y} may be given in any of the types a class list takes, and a label
## the model was not trained on is an error rather than a silent zero.
## @var{classname} and @var{caller} name the class and the method in that
## error, so a class reports its own name rather than this helper's.
##
## @end deftypefn

function gY = binaryIndex (ClassNames, Y, classname, caller)

  ## A character matrix names one class per row, so it is compared row by
  ## row rather than element by element, whichever side of the comparison it
  ## is on.
  if (ischar (ClassNames))
    ClassNames = cellstr (ClassNames);
  endif
  if (iscellstr (ClassNames) && ischar (Y))
    Y = cellstr (Y);
  endif
  if (iscellstr (ClassNames) && ! iscellstr (Y))
    error (strcat ("%s.%s: Y must be of the same type as the class", ...
                   " names."), classname, caller);
  endif

  [tf, gY] = ismember (Y(:), ClassNames);
  if (! all (tf))
    error (strcat ("%s.%s: Y must hold only classes the model was", ...
                   " trained on."), classname, caller);
  endif

endfunction
