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
## @deftypefn {Private Function} {@var{L} =} labelsFromIndex (@var{C}, @var{idx})
## Pick class labels out of a @code{ClassNames} by index, in its own type.
##
## @var{idx} numbers the classes, and @var{L} holds the label of each, in the
## type @var{C} is in.  A character matrix is indexed by row, every other
## accepted type by element, which is the same thing for a column.
## @end deftypefn

function L = labelsFromIndex (C, idx)

  if (ischar (C))
    L = C(idx, :);
  else
    L = C(idx);
  endif

endfunction
