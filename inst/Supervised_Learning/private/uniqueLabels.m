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
## @deftypefn {Private Function} {[@var{C}, @var{ia}, @var{ic}] =} uniqueLabels (@var{Y})
## The distinct class labels of a response, whatever type it is given in.
##
## @code{unique} compares a character matrix element by element, so a response
## naming its classes in the rows of a character matrix comes back as the
## distinct @emph{letters} of those names rather than the names themselves.
## Comparing by rows is what that response needs, and it is what every other
## accepted type already gets: a cell array of character vectors, a numeric
## column and a logical column each hold one label per row already.
##
## The three outputs carry @code{unique}'s own meanings, with rows in place of
## elements where the response is a character matrix.
## @end deftypefn

function [C, ia, ic] = uniqueLabels (Y)

  if (ischar (Y))
    [C, ia, ic] = unique (Y, "rows");
  else
    [C, ia, ic] = unique (Y);
  endif

endfunction
