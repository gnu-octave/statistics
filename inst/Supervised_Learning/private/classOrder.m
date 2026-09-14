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
## @deftypefn {Private Function} {[@var{C}, @var{gY}] =} classOrder (@var{Y}, @var{ClassNames})
## The classes of a response in the order a classifier lays them out, and the
## class of each observation.
##
## @var{C} holds the classes present in @var{Y}, in the type @var{Y} holds
## them: sorted, or in the order of @var{ClassNames} when that names any, a
## class it leaves out being left out of @var{C}.  Names match the classes as
## @code{namedClasses} matches them.  @var{gY} is the index into
## @var{C} of each observation, zero for one whose class was left out.
## @end deftypefn

function [C, gY] = classOrder (Y, ClassNames)

  C = uniqueLabels (Y);
  if (! isempty (ClassNames))
    pos = namedClasses (C, ClassNames);
    pos = pos(pos > 0);
    [~, first] = unique (pos, 'first');
    C = C(pos(sort (first)),:);
  endif
  gY = labelIndices (C, Y);

endfunction
