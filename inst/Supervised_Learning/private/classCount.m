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
## @deftypefn {Private Function} {@var{k} =} classCount (@var{C})
## How many classes a @code{ClassNames} holds, whatever type it is in.
##
## A cell array of character vectors, a numeric column and a logical column
## each hold one class per element, but a character matrix holds one per
## @emph{row} and as many columns as the longest name, so @code{numel} counts
## its characters rather than its classes.
## @end deftypefn

function k = classCount (C)

  if (ischar (C))
    k = rows (C);
  else
    k = numel (C);
  endif

endfunction
