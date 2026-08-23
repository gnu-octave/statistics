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
## @deftypefn {Private Function} {[@var{gY}, @var{errmsg}] =} labelIndices (@var{C}, @var{Y})
## The class each observation of a response belongs to, for the whole response.
##
## @var{gY} indexes @var{C} once per observation, and is zero where a label is
## not one of the classes.  @var{errmsg} is the body of the message the caller
## should raise, or empty when there is nothing wrong; the caller emits it
## under its own @code{class.method} name, as the package's shared validation
## helpers do.
##
## A response naming its classes in the rows of a character matrix is matched
## as whole names, with the padding such a matrix carries stripped from both
## sides.  Resolving every observation in one comparison is what this is for:
## asking per observation inside a loop costs a pass over the classes each
## time.
## @end deftypefn

function [gY, errmsg] = labelIndices (C, Y)

  errmsg = "";

  ## A character matrix names one class per row, so it is compared row by row
  ## rather than element by element, whichever side of the comparison it is on.
  if (ischar (C))
    C = cellstr (C);
  endif
  if (iscellstr (C) && ischar (Y))
    Y = cellstr (Y);
  endif

  if (iscellstr (C) && ! iscellstr (Y))
    gY = zeros (rows (Y), 1);
    errmsg = "Y must be of the same type as the class names.";
    return;
  endif

  [tf, gY] = ismember (Y(:), C(:));
  gY = gY(:);
  if (! all (tf))
    errmsg = "Y must hold only classes the model was trained on.";
  endif

endfunction
