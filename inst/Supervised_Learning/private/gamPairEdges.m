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
## @deftypefn {Private Function} {@var{PE} =} gamPairEdges (@var{TreeModel}, @var{DetectionEdges})
##
## The cut points each interaction pair of a boosted-tree GAM is held on.
##
## @var{PE} has one element per row of @code{@var{TreeModel}.Pairs}, each a
## @math{1x2} cell of the cut points on the pair's two predictors, as
## @code{gamboostinter} returns them in @qcode{PairEdges}.  A model saved
## before pair terms had grids of their own carries none; its surfaces are
## held on the detection grid @var{DetectionEdges}, one row vector per
## predictor, which is then read for both predictors of each pair.
##
## @end deftypefn

function PE = gamPairEdges (TreeModel, DetectionEdges)

  pairs = TreeModel.Pairs;
  if (isfield (TreeModel, 'PairEdges')
      && numel (TreeModel.PairEdges) == rows (pairs))
    PE = TreeModel.PairEdges;
    return;
  endif
  PE = cell (1, rows (pairs));
  for q = 1:rows (pairs)
    PE{q} = {DetectionEdges{pairs(q,1)}(:)', DetectionEdges{pairs(q,2)}(:)'};
  endfor

endfunction
