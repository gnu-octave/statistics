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
## @deftypefn {Private Function} {[@var{XC}, @var{E}, @var{cat}, @var{Levels}] =} gamCatCode (@var{TreeModel}, @var{X}, @var{BinEdges})
##
## Code the categorical predictors of a boosted-tree GAM for its engine.
##
## @code{@var{TreeModel}.CategoricalLevels} has one element per column of
## @var{X}: empty for a continuous predictor, else the sorted levels a
## categorical predictor held in training.  In @var{XC} each categorical column
## holds the index of its level, and @qcode{NaN} for a missing or unseen level,
## which the engine treats alike.  @var{E} is @var{BinEdges} with the cut
## points between those indices in place of each categorical predictor's
## entry, @var{cat} flags the categorical columns and @var{Levels} is the cell
## of levels.  A model without the field, fitted before categorical predictors
## were accepted or with none, leaves @var{X} and @var{BinEdges} unchanged.
##
## @end deftypefn

function [XC, E, cat, Levels] = gamCatCode (TreeModel, X, BinEdges)

  XC = X;
  E = {};
  if (nargin > 2)
    E = BinEdges;
  endif
  Levels = cell (1, columns (X));
  if (isstruct (TreeModel) && isfield (TreeModel, 'CategoricalLevels'))
    Levels = TreeModel.CategoricalLevels;
  endif
  cat = ! cellfun ('isempty', Levels);
  for j = find (cat)
    [~, loc] = ismember (X(:,j), Levels{j});
    code = loc;
    code(loc == 0) = NaN;
    XC(:,j) = code;
    if (nargin > 2)
      E{j} = (1:numel (Levels{j}) - 1) + 0.5;
    endif
  endfor

endfunction
