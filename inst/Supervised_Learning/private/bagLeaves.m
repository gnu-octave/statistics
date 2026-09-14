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
## @deftypefn {Private Function} {@var{L} =} bagLeaves (@var{M}, @var{X}, @var{trees})
##
## Node each observation comes to rest at in each tree of a bagged ensemble.
##
## @var{M} is a @code{TreeBagger} or @code{CompactTreeBagger} object, @var{X}
## the predictors and @var{trees} the indices of the trees to ask.  @var{L} is
## an @math{NxT} matrix holding, for each observation and tree, the number of
## the node the tree's @code{predict} brings it to.
##
## @end deftypefn

function L = bagLeaves (M, X, trees)

  L = zeros (rows (X), numel (trees));
  isclass = strcmp (M.Method, 'classification');
  for j = 1:numel (trees)
    if (isclass)
      [~, ~, L(:,j)] = predict (M.Trees{trees(j)}, X);
    else
      [~, L(:,j)] = predict (M.Trees{trees(j)}, X);
    endif
  endfor

endfunction
