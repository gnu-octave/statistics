## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{folds} =} __stratified_group_folds__ (@var{classes}, @var{gidx}, @var{k})
##
## Assign whole groups to @var{k} folds while keeping the class balance even.
##
## @var{classes} numbers each observation by its class and @var{gidx} by its
## group; @var{folds} returns the fold each observation is assigned to.  Every
## observation of a group lands in the same fold, so no group is split, and the
## folds are chosen to keep each class spread as evenly as the groups allow.
##
## The two demands conflict: once a group is atomic it carries whatever class
## mix it has, so the class balance can only be approximate.  Groups are placed
## largest first, each into the fold whose class counts it disturbs least,
## which is the greedy assignment scikit-learn's @code{StratifiedGroupKFold}
## uses.
##
## @end deftypefn

function folds = __stratified_group_folds__ (classes, gidx, k)

  nclasses = max (classes);
  ngroups = max (gidx);

  ## Class counts of every group
  gcount = zeros (ngroups, nclasses);
  for g = 1:ngroups
    for c = 1:nclasses
      gcount(g,c) = sum (gidx == g & classes == c);
    endfor
  endfor

  ## Largest groups first: they constrain the balance most, so they are placed
  ## while every fold is still empty enough to take them
  [~, order] = sort (sum (gcount, 2), 'descend');

  fcount = zeros (k, nclasses);   # class counts accumulated per fold
  gfold = zeros (ngroups, 1);
  for i = 1:ngroups
    g = order(i);
    best = 1;
    best_cost = Inf;
    for f = 1:k
      trial = fcount;
      trial(f,:) += gcount(g,:);
      ## How unevenly each class is spread over the folds, averaged
      cost = mean (std (trial, 0, 1));
      ## Break ties toward the emptiest fold, so the folds stay similar in size
      if (cost < best_cost
          || (cost == best_cost && sum (trial(f,:)) < sum (fcount(best,:))))
        best = f;
        best_cost = cost;
      endif
    endfor
    fcount(best,:) += gcount(g,:);
    gfold(g) = best;
  endfor

  folds = gfold(gidx);

endfunction
