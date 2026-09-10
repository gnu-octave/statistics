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
## @deftypefn {Private Function} {@var{S} =} treeCollapse (@var{Children}, @var{Parent}, @var{CutPredictorIndex}, @var{CutPoint}, @var{nodes})
##
## Turn the named branch nodes of a tree into leaves and renumber what is
## left.
##
## Everything below a collapsed node is discarded.  The surviving nodes keep
## their order, so a parent still carries a lower number than its children,
## which is what lets one forward pass find what is reachable.
##
## @var{S} carries the rewritten @qcode{Children}, @qcode{Parent},
## @qcode{CutPredictorIndex} and @qcode{CutPoint}, along with @qcode{keep},
## the indices of the surviving nodes in the old numbering.  A caller subsets
## whatever else it holds per node by @qcode{keep}.
##
## @seealso{ClassificationTree, RegressionTree}
## @end deftypefn

function S = treeCollapse (Children, Parent, CutPredictorIndex, CutPoint, nodes)

  n = rows (Children);
  kid = Children;
  nodes = nodes(kid(nodes,1) > 0);
  kid(nodes,:) = 0;

  ## What is still reachable from the root
  keep = false (n, 1);
  keep(1) = true;
  for ii = 1:n
    if (keep(ii) && kid(ii,1) > 0)
      keep(kid(ii,1)) = true;
      keep(kid(ii,2)) = true;
    endif
  endfor

  idx = find (keep);
  ## Offset by one so that a zero, meaning no such node, renumbers to zero
  renum = zeros (n + 1, 1);
  renum(idx + 1) = 1:numel (idx);
  cutvar = CutPredictorIndex;
  cutval = CutPoint;
  cutvar(nodes) = 0;
  cutval(nodes) = NaN;

  S.keep = idx;
  S.Children = reshape (renum(kid(idx,:) + 1), numel (idx), 2);
  S.Parent = renum(Parent(idx) + 1);
  S.CutPredictorIndex = cutvar(idx);
  S.CutPoint = cutval(idx);

endfunction
