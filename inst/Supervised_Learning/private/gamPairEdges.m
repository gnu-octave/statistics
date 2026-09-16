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
## @deftypefn {Private Function} {[@var{PE}, @var{PM}] =} gamPairEdges (@var{TreeModel}, @var{DetectionEdges})
##
## The cut points and missing-value entries of a boosted-tree GAM's pairs.
##
## @var{PE} has one element per row of @code{@var{TreeModel}.Pairs}, each a
## @math{1x2} cell of the cut points on the pair's two predictors, and @var{PM}
## the matching @math{1x3} cells of what a row missing the first predictor,
## the second or both takes, as @code{gamboostinter} returns them in
## @qcode{PairEdges} and @qcode{PairMissing}.  A model saved before pair terms
## had grids of their own holds its surfaces on the detection grid
## @var{DetectionEdges}, one row vector per predictor, which is then read for
## both predictors of each pair; a model saved before pairs held missing
## values gets zeros for them, so such a pair still contributes nothing.
##
## @end deftypefn

function [PE, PM] = gamPairEdges (TreeModel, DetectionEdges)

  pairs = TreeModel.Pairs;
  np = rows (pairs);
  if (isfield (TreeModel, 'PairEdges') && numel (TreeModel.PairEdges) == np)
    PE = TreeModel.PairEdges;
  else
    PE = cell (1, np);
    for q = 1:np
      PE{q} = {DetectionEdges{pairs(q,1)}(:)', DetectionEdges{pairs(q,2)}(:)'};
    endfor
  endif
  if (isfield (TreeModel, 'PairMissing')
      && numel (TreeModel.PairMissing) == np)
    PM = TreeModel.PairMissing;
  else
    PM = cell (1, np);
    for q = 1:np
      [nj, nk] = size (TreeModel.PairValues{q});
      PM{q} = {zeros(1, nk), zeros(1, nj), 0};
    endfor
  endif

endfunction
