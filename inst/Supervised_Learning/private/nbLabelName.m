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
##

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{s} =} nbLabelName (@var{C}, @var{i})
##
## One class's name, as a message should print it.
##
## @var{C} is the model's @qcode{ClassNames} in any of the forms a response
## may take, and @var{i} indexes it.  A textual class prints as its name and
## every other kind as its value, which is what MATLAB's own messages do:
## a numeric response names @qcode{class 2} where a cellstr one names
## @qcode{class beta}.
##
## @seealso{ClassificationNaiveBayes}
## @end deftypefn

function s = nbLabelName (C, i)

  if (iscellstr (C))
    s = C{i};
  elseif (ischar (C))
    s = strtrim (C(i,:));
  elseif (islogical (C))
    s = mat2str (C(i));
  else
    s = num2str (C(i));
  endif

endfunction
