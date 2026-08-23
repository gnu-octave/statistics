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
##
## -*- texinfo -*-
## @deftypefn {Private Function} {@var{idx} =} labelIndex (@var{C}, @var{Y}, @var{i})
## Which class the @var{i}th observation of a response belongs to.
##
## @var{idx} indexes @var{C}, and is empty where the label is not one of the
## classes, which is what the callers test for.
##
## A response naming its classes in the rows of a character matrix has to be
## read a row at a time and matched as a whole name: indexing it by element
## takes a single character, and @code{ismember} between two character
## matrices compares them character by character.  The padding such a matrix
## carries is stripped from both sides so that a padded label still matches
## the class name it names.
## @end deftypefn

function idx = labelIndex (C, Y, i)

  if (ischar (C) || ischar (Y))
    if (ischar (Y))
      y = strtrim (Y(i,:));
    else
      y = char (Y(i));
    endif
    idx = find (strcmp (cellstr (C), y));
  else
    idx = find (ismember (C, Y(i)));
  endif

endfunction
