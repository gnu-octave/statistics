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
## @deftypefn {Private Function} {@var{P} =} bagTreeOutputs (@var{M}, @var{X}, @var{trees})
##
## What each tree of a bagged ensemble says about each observation.
##
## @var{M} is a @code{TreeBagger} or @code{CompactTreeBagger} object, @var{X}
## the predictors and @var{trees} the indices of the trees to ask.  For a
## classification ensemble @var{P} is @math{NxKxT}, the scores of each tree
## laid out in the ensemble's class order, a class absent from a tree's bag
## scoring zero.  For a regression ensemble @var{P} is @math{Nx1xT}, the
## predicted responses.
##
## @end deftypefn

function P = bagTreeOutputs (M, X, trees)

  N = rows (X);
  T = numel (trees);
  if (strcmp (M.Method, 'classification'))
    P = zeros (N, classCount (M.ClassNames), T);
    for j = 1:T
      [~, s] = predict (M.Trees{trees(j)}, X);
      P(:, M.TreeClassIdx{trees(j)}, j) = s;
    endfor
  else
    P = zeros (N, 1, T);
    for j = 1:T
      P(:, 1, j) = predict (M.Trees{trees(j)}, X);
    endfor
  endif

endfunction
