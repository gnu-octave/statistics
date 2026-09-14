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
## @deftypefn {Private Function} {@var{P} =} bagProximity (@var{M}, @var{X}, @var{trees}, @var{NumPrint})
##
## Proximity of the observations under a bagged ensemble.
##
## @var{M} is a @code{TreeBagger} or @code{CompactTreeBagger} object, @var{X}
## the predictors and @var{trees} the indices of the trees to use.  @var{P} is
## the symmetric @math{NxN} matrix whose element @math{(i,j)} is the share of
## those trees that bring observations @math{i} and @math{j} to the same node.
## When @var{NumPrint} is given and positive, a line saying how many trees are
## done is printed after every that many trees.
##
## @end deftypefn

function P = bagProximity (M, X, trees, NumPrint = 0)

  N = rows (X);
  T = numel (trees);
  P = zeros (N);
  for j = 1:T
    L = bagLeaves (M, X, trees(j));
    P += L == L';
    if (NumPrint > 0 && mod (j, NumPrint) == 0)
      printf ("Tree %d done.\n", j);
    endif
  endfor
  P /= T;

endfunction
