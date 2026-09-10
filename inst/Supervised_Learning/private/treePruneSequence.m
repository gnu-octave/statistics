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
## @deftypefn {Private Function} {[@var{PruneList}, @var{PruneAlpha}] =} treePruneSequence (@var{Children}, @var{Parent}, @var{risk}, @var{held})
##
## The cost complexity pruning sequence of a tree, by weakest link.
##
## A branch node's link is the risk it would take on as a leaf, less the risk
## its subtree carries now, spread over the leaves the subtree would give up.
## The weakest link is pruned, the sequence repeats on what is left, and the
## level a node is pruned at is its place in that sequence.  Leaves never
## carry a level, so a tree of @math{B} branch nodes gives @math{B} levels at
## most and one more alpha than levels, the first of which is zero and stands
## for the unpruned tree.
##
## @var{risk} is the risk each node would carry as a leaf, on whatever scale
## the learner measures it: the expected misclassification cost for a
## classifier and the squared error about the node's mean for a regression.
##
## @var{held} is the risk a node carries on account of the observations that
## stop there, missing the predictor it cuts on.  Those observations are in
## neither child, so a subtree's risk is its children's plus this, and
## leaving it out overstates every link above a node that holds rows back.
## Measured on R2024a, where one row of @code{carsmall} missing its
## horsepower moves an alpha from 5.99325416896717 to 6.31708345021804.
##
## A subtree that costs nothing to give up is no step of the sequence.  It is
## what merging leaves would have removed, and MATLAB records neither a level
## nor an alpha for it: measured on R2024a with MergeLeaves off, where such a
## subtree survives to be seen and the tree carries eleven nodes and the
## merged tree's five alphas, not six.
##
## @seealso{ClassificationTree, RegressionTree}
## @end deftypefn

function [PruneList, PruneAlpha] = treePruneSequence (Children, Parent, ...
                                                      risk, held)

  n = rows (Children);
  PruneList = zeros (n, 1);
  alphas = [];
  kidl = Children(:,1);
  kidr = Children(:,2);
  level = 0;
  ## The branches given up at no cost, which are no part of the sequence
  free = false (n, 1);

  while (kidl(1) != 0)

    ## Only what is still reachable from the root can be a candidate: a
    ## pruned branch takes its whole subtree with it, and orphans left in
    ## would go on offering links and split one level into several.
    reach = false (n, 1);
    reach(1) = true;
    for ii = 1:n
      if (reach(ii) && kidl(ii) != 0)
        reach(kidl(ii)) = true;
        reach(kidr(ii)) = true;
      endif
    endfor

    ## Subtree risk and leaf count, deepest first
    subrisk = risk;
    subleaf = ones (n, 1);
    for ii = n:-1:1
      if (kidl(ii) != 0)
        subrisk(ii) = subrisk(kidl(ii)) + subrisk(kidr(ii)) + held(ii);
        subleaf(ii) = subleaf(kidl(ii)) + subleaf(kidr(ii));
      endif
    endfor

    cand = find (reach & kidl != 0);
    if (isempty (cand))
      break;
    endif
    link = (risk(cand) - subrisk(cand)) ./ (subleaf(cand) - 1);
    weakest = min (link);

    ## Every branch whose link is the weakest goes at this level, not just
    ## one of them, and links equal in exact arithmetic can differ in their
    ## last bits as split gains do.
    tol = treeRiskTol () * max (abs (risk(1)), 1);
    cut = weakest + abs (weakest) * treeRiskTol () + tol;
    gone = cand(link <= cut);
    if (weakest > tol)
      level++;
      alphas(end+1, 1) = weakest;
      PruneList(gone) = level;
    else
      free(gone) = true;
    endif
    kidl(gone) = 0;
    kidr(gone) = 0;

  endwhile

  ## A branch inside a subtree given up at no cost is not part of the
  ## sequence either, its ancestor having left it at no level.
  for ii = 2:n
    free(ii) = free(ii) || free(Parent(ii));
  endfor

  ## A branch that lost an ancestor never came up for pruning on its own
  ## account, but it stopped being a branch when that ancestor went, and that
  ## is the level it carries.
  for ii = 1:n
    if (Children(ii,1) == 0 || PruneList(ii) != 0 || free(ii))
      continue;
    endif
    a = Parent(ii);
    while (a > 0 && PruneList(a) == 0)
      a = Parent(a);
    endwhile
    if (a > 0)
      PruneList(ii) = PruneList(a);
    endif
  endfor

  PruneAlpha = [0; alphas];

endfunction

## How close two risks must be to count as equal.  Links that are equal in
## exact arithmetic differ in their last bits once they have been through a
## division, and treating them as distinct would split one pruning level into
## several.
function t = treeRiskTol ()

  t = 1e-12;

endfunction
