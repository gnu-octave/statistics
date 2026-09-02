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
## @deftypefn {Private Function} {@var{tf} =} labelsKnown (@var{Y}, @var{C})
## Whether every distinct label in a response is one of the known classes.
##
## @code{ismember} compares two character matrices character by character, so
## a response naming its classes in the rows of one is answered against the
## letters of the class names rather than the names.  Anything textual is
## compared as whole names here; every other accepted type compares by value.
## @end deftypefn

function tf = labelsKnown (Y, C)

  ## Delegated so the matching rules live in one place: labelIndices reports
  ## a zero for any label that is not one of the classes.
  tf = all (labelIndices (C, Y) > 0);

endfunction
