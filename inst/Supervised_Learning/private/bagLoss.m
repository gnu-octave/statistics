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
## @deftypefn {Private Function} {@var{out} =} bagLoss (@var{what}, @var{M}, @var{X}, @var{Y}, @var{args}, @var{classname}, @var{fixedUse}, @var{W0})
##
## Error, margin or mean margin of a bagged ensemble.
##
## @var{what} is @qcode{'error'}, @qcode{'margin'} or @qcode{'meanMargin'}.
## @var{M} is a @code{TreeBagger} or @code{CompactTreeBagger} object, @var{X}
## and @var{Y} the data, @var{args} the Name-Value pairs of the call and
## @var{classname} the class and method named in the error messages.
## @var{fixedUse}, when not empty, is the out-of-bag matrix, and
## @qcode{'UseInstanceForTree'} may not then be given.  @var{W0} weighs the
## observations when no @qcode{'Weights'} are given, and is empty for uniform
## weights.
##
## The error is the weighted share of misclassified observations, or the
## weighted mean squared error, over the observations that have a prediction.
## An observation no tree may answer for takes @code{DefaultYfit} with the
## prior as its scores; when @code{DefaultYfit} is the missing label it has no
## prediction and is left out.  The error is a column with one element per
## step in @qcode{'cumulative'} and @qcode{'individual'} mode and a scalar in
## @qcode{'ensemble'} mode; the margin is @math{NxL} and the mean margin
## @math{1xL}.
##
## @end deftypefn

function out = bagLoss (what, M, X, Y, args, classname, fixedUse, W0)

  isclass = strcmp (M.Method, 'classification');
  if (! isclass && ! strcmp (what, 'error'))
    error ("%s: margins are defined only for classification.", classname);
  endif
  if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
    error ("%s: X must be a real numeric matrix.", classname);
  endif
  if (columns (X) != numel (M.PredictorNames))
    error ("%s: X must have one column per predictor.", classname);
  endif
  N = rows (X);
  if (rows (Y) != N)
    error ("%s: X and Y must have the same number of rows.", classname);
  endif

  allowed = {'Mode', 'Trees', 'TreeWeights'};
  if (isempty (fixedUse))
    allowed{end+1} = 'UseInstanceForTree';
  endif
  if (! strcmp (what, 'margin'))
    allowed{end+1} = 'Weights';
  endif
  [o, errmsg] = bagArgs (args, M.NumTrees, N, allowed);
  if (! isempty (errmsg))
    error ("%s: %s", classname, errmsg);
  endif
  if (! isempty (fixedUse))
    o.use = fixedUse(:, o.trees);
  endif
  if (! isempty (o.w))
    w = o.w;
  elseif (! isempty (W0))
    w = W0(:);
  else
    w = ones (N, 1);
  endif

  P = bagTreeOutputs (M, X, o.trees);
  [A, ~, none] = bagCombine (P, o.use, o.tw, o.mode);
  L = columns (none);

  if (isclass)
    [gY, errmsg] = labelIndices (M.ClassNames, Y);
    if (! isempty (errmsg))
      error ("%s: %s", classname, errmsg);
    endif
    gY = gY(:);
    K = columns (A);
    for l = 1:L
      A(none(:,l),:,l) = repmat (M.DefaultScore, sum (none(:,l)), 1);
    endfor
    if (strcmp (what, 'error'))
      [~, idx] = max (A, [], 2);
      idx = reshape (idx, N, L);
      idx(none) = M.DefaultIndex;
      valid = idx > 0;
      miss = double (idx != gY);
      out = sum (w .* miss .* valid, 1) ./ sum (w .* valid, 1);
      out = out(:);
    else
      m = marginsOf (A, gY, L);
      if (M.DefaultIndex == 0)
        m(none) = NaN;
      endif
      if (strcmp (what, 'margin'))
        out = m;
      else
        have = ! isnan (m);
        m(! have) = 0;
        out = sum (w .* m .* have, 1) ./ sum (w .* have, 1);
      endif
    endif
  else
    if (! (isnumeric (Y) && isreal (Y) && isvector (Y)))
      error ("%s: Y must be a real numeric vector.", classname);
    endif
    y = double (Y(:));
    Yh = reshape (A, N, L);
    Yh(none) = M.DefaultYfit;
    valid = ! (isnan (Yh) | isnan (y));
    d = Yh - y;
    d(! valid) = 0;
    out = sum (w .* d .^ 2 .* valid, 1) ./ sum (w .* valid, 1);
    out = out(:);
  endif

endfunction
