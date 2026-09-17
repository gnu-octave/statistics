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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{v} =} pdTraverse (@var{T}, @var{vars}, @var{q}, @var{leaf})
##
## What a decision tree answers on average, with some predictors held fixed.
##
## This is Friedman's weighted traversal, and it is what partial dependence
## means for a tree: the average of the tree's answer over the distribution it
## was fitted on, with the predictors named by @var{vars} held at @var{q}.  It
## is not an average over any set of observations and takes none, which is why
## it answers for a compact model, and why MATLAB R2024a gives the same value
## for a tree whatever @qcode{Data} it is handed.
##
## A node that cuts one of the held predictors is followed to the child that
## value reaches.  A node that cuts any other is followed **both** ways, each
## child taking the share of the training weight that reached it.  A held
## value that matches no level of a categorical cut stops there, as it stops
## an observation in @code{predict}.
##
## @var{T} is the tree, @var{vars} the held predictors and @var{q} their
## values in the same order.  @var{leaf} holds what the tree answers at each
## node, one row per node: a column of means for a regression tree and a
## column per class for a classifier.
##
## @var{v} is a row, one element per column of @var{leaf}.
##
## @end deftypefn

function v = pdTraverse (T, vars, q, leaf)

  v = zeros (1, columns (leaf));
  node = 1;
  w = 1;
  while (! isempty (node))
    n = node(1);
    wn = w(1);
    node(1) = [];
    w(1) = [];
    j = T.CutPredictorIndex(n);
    if (j == 0)
      v += wn * leaf(n,:);
      continue;
    endif
    L = T.Children(n,1);
    R = T.Children(n,2);
    k = find (vars == j, 1);
    if (isempty (k))
      ## Not held: both ways, in the proportion the training weight took
      pn = T.NodeProbability(n);
      if (pn > 0)
        pL = T.NodeProbability(L) / pn;
      else
        pL = 0.5;
      endif
      node(end+1:end+2) = [L, R];
      w(end+1:end+2) = [wn * pL, wn * (1 - pL)];
    else
      [child, stop] = pdChild (T, n, q(k), L, R);
      if (stop)
        v += wn * leaf(n,:);
      else
        node(end+1) = child;
        w(end+1) = wn;
      endif
    endif
  endwhile

endfunction

## Which child a held value reaches, and whether it reaches neither.
function [child, stop] = pdChild (T, n, val, L, R)

  stop = false;
  child = 0;
  if (! isempty (T.CutCategories) && ! isempty (T.CutCategories{n,1}))
    if (any (T.CutCategories{n,1} == val))
      child = L;
    elseif (any (T.CutCategories{n,2} == val))
      child = R;
    else
      stop = true;
    endif
  elseif (isnan (val))
    stop = true;
  elseif (val < T.CutPoint(n))
    child = L;
  else
    child = R;
  endif

endfunction
