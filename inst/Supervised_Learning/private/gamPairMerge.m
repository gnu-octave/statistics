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
## @deftypefn {Private Function} {[@var{PE}, @var{PV}, @var{PM}] =} gamPairMerge (@var{PE1}, @var{PV1}, @var{PM1}, @var{PE2}, @var{PV2}, @var{PM2})
##
## Add two sets of interaction surfaces held on different grids.
##
## Each pair's two surfaces are read on the union of their cut points, at a
## point strictly inside every cell, and summed, so the result predicts what
## the two sets predict together.  What a row missing one predictor takes is
## read and summed the same way on the union of the other predictor's cut
## points, and what a row missing both takes is summed.  @var{PE1} and
## @var{PE2} are cells of per-pair cut points and @var{PM1} and @var{PM2} of
## missing-value entries, as @code{gamPairEdges} returns them, and @var{PV1}
## and @var{PV2} the matching cells of surfaces.
##
## @end deftypefn

function [PE, PV, PM] = gamPairMerge (PE1, PV1, PM1, PE2, PV2, PM2)

  np = numel (PV1);
  PE = cell (1, np);
  PV = cell (1, np);
  PM = cell (1, np);
  for q = 1:np
    ej = unique ([PE1{q}{1}(:); PE2{q}{1}(:)])';
    ek = unique ([PE1{q}{2}(:); PE2{q}{2}(:)])';
    xj = cellPoints (ej);
    xk = cellPoints (ek);
    PV{q} = surfaceAt (PE1{q}, PV1{q}, xj, xk) ...
            + surfaceAt (PE2{q}, PV2{q}, xj, xk);
    PE{q} = {ej, ek};
    ## A row missing the first predictor is read by cell of the second, and
    ## the other way round.
    mj = lineAt (PE1{q}{2}, PM1{q}{1}, xk) + lineAt (PE2{q}{2}, PM2{q}{1}, xk);
    mk = lineAt (PE1{q}{1}, PM1{q}{2}, xj) + lineAt (PE2{q}{1}, PM2{q}{2}, xj);
    PM{q} = {mj, mk, PM1{q}{3} + PM2{q}{3}};
  endfor

endfunction

## A point strictly inside each cell the cut points E define.
function x = cellPoints (e)

  if (isempty (e))
    x = 0;
  else
    x = [e(1) - 1, (e(1:end-1) + e(2:end)) / 2, e(end) + 1];
  endif

endfunction

## The surface S held on the cut points E, read at every pair of XJ and XK.
function V = surfaceAt (E, S, xj, xk)

  a = 1 + sum (xj(:) > E{1}(:)', 2);
  b = 1 + sum (xk(:) > E{2}(:)', 2);
  V = S(a, b);

endfunction

## The values L held on the cut points E, read at every point of X.
function v = lineAt (E, L, x)

  idx = 1 + sum (x(:) > E(:)', 2);
  v = L(idx);
  v = v(:)';

endfunction
