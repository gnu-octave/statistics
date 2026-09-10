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
## @deftypefn {Private Function} {@var{S} =} treeNodeStats (@var{ClassShare}, @var{Prior}, @var{Cost}, @var{crit}, @var{ClassNames})
##
## The node statistics of a classification tree that depend on the prior or
## the cost.
##
## @var{ClassShare} holds, per node and per class, the share of that class's
## total weight which reached the node, and is free of both the prior and the
## cost: multiplying a class's column by a constant leaves it unchanged.
## Everything below is derived from it, so the same table serves whatever
## @var{Prior} and @var{Cost} are assigned afterwards.
##
## @var{crit} is the split criterion the tree was grown under, @qcode{'gdi'}
## or @qcode{'deviance'}, which is the impurity the risk is measured by.
##
## @var{S} is a structure carrying @qcode{NodeProbability},
## @qcode{ClassProbability}, @qcode{NodeError}, @qcode{NodeClass},
## @qcode{NodeRisk} and @qcode{NodeProbAdj}, the last being the probability
## of reaching each node measured on the cost-adjusted weights, which is the
## scale @qcode{NodeRisk} lives on.
##
## The cost enters the risk and nothing else.  Measured on R2024a: a
## classification tree reports every other node statistic on the unadjusted
## weights and measures the risk on the distribution each class's weight is
## scaled by the total cost of misclassifying it, which is the classical way
## a cost matrix enters CART.  Under the default cost every class is scaled
## alike and the adjustment falls away.
##
## @seealso{ClassificationTree, CompactClassificationTree}
## @end deftypefn

function S = treeNodeStats (ClassShare, Prior, Cost, crit, ClassNames)

  cw = ClassShare .* Prior(:)';
  nw = sum (cw, 2);
  CP = zeros (size (cw));
  nz = nw > 0;
  CP(nz,:) = cw(nz,:) ./ nw(nz);
  S.NodeProbability = nw;
  S.ClassProbability = CP;

  ## The class of least expected misclassification cost, the first of the
  ## class names keeping a tie, which is what min returns
  [err, k] = min (CP * Cost, [], 2);
  S.NodeError = err;
  names = classText (ClassNames);
  S.NodeClass = names(k);

  ## The impurity of the cost-adjusted distribution, weighted by the adjusted
  ## probability of reaching the node
  aw = cw .* sum (Cost, 2)';
  anw = sum (aw, 2);
  AP = zeros (size (aw));
  nz = anw > 0;
  AP(nz,:) = aw(nz,:) ./ anw(nz);
  if (anw(1) > 0)
    S.NodeProbAdj = anw / anw(1);
  else
    S.NodeProbAdj = zeros (rows (aw), 1);
  endif
  S.NodeRisk = S.NodeProbAdj .* nodeImpurity (AP, crit);

endfunction

## The class names as text, whatever type they are carried in, which is what
## NodeClass reports and what view prints.
function s = classText (C)

  if (iscellstr (C))
    s = C(:);
  elseif (ischar (C))
    s = cellstr (C);
  else
    s = arrayfun (@(v) num2str (v), C(:), 'UniformOutput', false);
  endif

endfunction

## The impurity of a distribution, by the criterion the tree was grown under.
## The deviance is the entropy in bits halved, which is the scale MATLAB
## reports NodeRisk on: measured on R2024a, an iris tree grown under
## 'deviance' reports a root risk of 0.792481250360577, which is log2 (3) / 2.
function imp = nodeImpurity (P, crit)

  if (strcmp (crit, 'gdi'))
    imp = 1 - sum (P .^ 2, 2);
  else
    L = P;
    L(L <= 0) = 1;      # a class of no weight contributes nothing
    imp = -sum (P .* log2 (L), 2) / 2;
  endif

endfunction
