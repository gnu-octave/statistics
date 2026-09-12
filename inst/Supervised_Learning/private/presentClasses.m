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
## @deftypefn {Private Function} {[@var{gY}, @var{gnY}, @var{glY}] =} presentClasses (@var{gY}, @var{gnY}, @var{glY})
##
## Keep only the groups of @code{grp2idx} that hold an observation.
##
## @code{grp2idx} gives a categorical response one group per category, used or
## not, so a response that has lost a class along with its rows, a subset of a
## categorical array included, would still count it.  The unused groups are
## removed and the indices in @var{gY} renumbered to match.  The class labels
## in @var{glY} keep the categories of the response.
##
## @end deftypefn

function [gY, gnY, glY] = presentClasses (gY, gnY, glY)

  used = false (numel (gnY), 1);
  used(gY(! isnan (gY))) = true;
  if (all (used))
    return;
  endif
  remap = cumsum (used);
  have = ! isnan (gY);
  gY(have) = remap(gY(have));
  gnY = gnY(used);
  glY = glY(used,:);

endfunction
