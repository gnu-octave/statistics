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
## @deftypefn {Private Function} {[@var{Yfit}, @var{S}, @var{SD}] =} bagPredict (@var{M}, @var{X}, @var{args}, @var{classname}, @var{fixedUse})
##
## Predictions of a bagged ensemble, for @code{predict} and @code{oobPredict}.
##
## @var{M} is a @code{TreeBagger} or @code{CompactTreeBagger} object and
## @var{args} the Name-Value pairs of the call.  @var{classname} names the
## class and method in the error messages.  @var{fixedUse}, when not empty, is
## the out-of-bag matrix that decides which tree answers for which
## observation, and only @qcode{'Trees'} may then be given.
##
## For classification @var{Yfit} holds the labels, @var{S} the scores and
## @var{SD} their standard deviations over the trees.  For regression
## @var{Yfit} holds the responses and @var{S} their standard deviations, and
## @var{SD} is empty.  An observation no tree may answer for takes
## @code{DefaultYfit}, with the prior as its scores and no deviation.
##
## @end deftypefn

function [Yfit, S, SD] = bagPredict (M, X, args, classname, fixedUse)

  if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
    error ("%s: X must be a real numeric matrix.", classname);
  endif
  if (columns (X) != numel (M.PredictorNames))
    error ("%s: X must have one column per predictor.", classname);
  endif

  if (isempty (fixedUse))
    allowed = {'Trees', 'TreeWeights', 'UseInstanceForTree'};
  else
    allowed = {'Trees'};
  endif
  [o, errmsg] = bagArgs (args, M.NumTrees, rows (X), allowed);
  if (! isempty (errmsg))
    error ("%s: %s", classname, errmsg);
  endif
  if (! isempty (fixedUse))
    o.use = fixedUse(:, o.trees);
  endif

  P = bagTreeOutputs (M, X, o.trees);
  [A, SD, none] = bagCombine (P, o.use, o.tw, 'ensemble');

  if (strcmp (M.Method, 'classification'))
    A(none,:) = repmat (M.DefaultScore, sum (none), 1);
    [~, idx] = max (A, [], 2);
    idx(none) = M.DefaultIndex;
    Yfit = missingLabels (M.ClassNames, rows (X));
    have = idx > 0;
    Yfit(have,:) = labelsFromIndex (M.ClassNames, idx(have));
    S = A;
  else
    A(none) = M.DefaultYfit;
    Yfit = A;
    S = SD;
    SD = [];
  endif

endfunction
