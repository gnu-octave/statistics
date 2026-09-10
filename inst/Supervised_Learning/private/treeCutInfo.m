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
## @deftypefn {Private Function} {@var{S} =} treeCutInfo (@var{CutPredictorIndex}, @var{PredictorNames})
##
## The descriptions of a tree's cuts, derived from its node table.
##
## @var{S} is a structure carrying @qcode{IsBranchNode}, @qcode{CutPredictor},
## @qcode{CutType} and @qcode{CutCategories}, along with
## @qcode{CategoricalSplit} and the six @qcode{Surrogate} properties, which
## are empty in the shapes MATLAB reports for a tree with neither categorical
## predictors nor surrogate splits.
##
## Deriving them rather than storing them is what keeps them from falling out
## of step with the node table, which pruning rewrites.
##
## @seealso{ClassificationTree, CompactClassificationTree}
## @end deftypefn

function S = treeCutInfo (CutPredictorIndex, PredictorNames)

  n = numel (CutPredictorIndex);
  br = CutPredictorIndex(:) > 0;
  cutname = repmat ({''}, n, 1);
  cuttype = repmat ({''}, n, 1);
  if (any (br))
    cutname(br) = PredictorNames(CutPredictorIndex(br));
    cuttype(br) = {'continuous'};
  endif

  S.IsBranchNode = br;
  S.CutPredictor = cutname;
  S.CutType = cuttype;
  S.CutCategories = repmat ({zeros(0, 0)}, n, 2);

  ## Categorical predictors and surrogate splits are not implemented, and
  ## these are the shapes MATLAB reports for a tree that has neither.
  S.CategoricalSplit = cell (0, 0);
  S.SurrogateCutCategories = cell (0, 0);
  S.SurrogateCutFlip = cell (0, 0);
  S.SurrogateCutPoint = cell (0, 0);
  S.SurrogateCutType = cell (0, 0);
  S.SurrogateCutPredictor = cell (0, 1);
  S.SurrogatePredictorAssociation = cell (0, 0);

endfunction
