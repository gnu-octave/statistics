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
## @deftypefn {Private Function} {[@var{dcrit}, @var{nsplit}] =} bagSplitStats (@var{Trees}, @var{p})
##
## Split criterion contributions and split counts of a bagged ensemble.
##
## @var{Trees} is the cell array of compact trees and @var{p} the number of
## predictors.  @var{dcrit} is the mean over the trees of each tree's
## @code{predictorImportance}.  @var{nsplit} is the sum over the trees of the
## share of each tree's branch nodes that split on each predictor; a tree
## without branch nodes adds nothing to either.  Both are @math{1xP}.
##
## @end deftypefn

function [dcrit, nsplit] = bagSplitStats (Trees, p)

  dcrit = zeros (1, p);
  nsplit = zeros (1, p);
  T = numel (Trees);
  for t = 1:T
    tr = Trees{t};
    branch = find (tr.Children(:,1) > 0);
    if (isempty (branch))
      continue;
    endif
    dcrit += predictorImportance (tr);
    nsplit += accumarray (tr.CutPredictorIndex(branch)(:), 1, [p, 1])' ...
              / numel (branch);
  endfor
  if (T > 0)
    dcrit /= T;
  endif

endfunction
