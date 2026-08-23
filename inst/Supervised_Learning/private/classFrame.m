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
## @deftypefn {Private Function} {@var{F} =} classFrame (@var{X}, @var{Y}, @var{ClassNames}, @var{Prior}, @var{Cost}, @var{Weights}, @var{classname})
##
## Resolve the classes, the prior, the cost and the observation weights of a
## binary classifier.
##
## The four linear and kernel classifiers all begin the same way: group the
## response, keep only the named classes, drop the rows that are not
## complete, and turn what survives into a signed response and a weight per
## observation.  @var{F} is a structure carrying the results, so that the
## classifier and its cross-validated counterpart cannot drift apart.
##
## Fields: @qcode{X} and @qcode{Y}, the retained data; @qcode{RowsUsed}, a
## logical over the rows as supplied; @qcode{gY}, the class index of each
## retained observation; @qcode{ClassNames}, in the type of @var{Y};
## @qcode{Prior} and @qcode{Cost}, as the model reports them;
## @qcode{Weights}, the retained weights before normalization; @qcode{W},
## normalized within each class to that class's cost-adjusted prior and
## summing to one; @qcode{y}, @math{+1} for the second class and @math{-1}
## for the first; and @qcode{n} and @qcode{p}.
##
## @var{classname} names the caller in every error message.
##
## @end deftypefn

function F = classFrame (X, Y, ClassNames, Prior, Cost, Weights, classname)

  if (! (isnumeric (X) && isreal (X) && ismatrix (X) && ndims (X) == 2))
    error ("%s: invalid values in X.", classname);
  endif
  if (isempty (X))
    error ("%s: X is empty.", classname);
  endif
  if (rows (X) != rows (Y))
    error ("%s: number of rows in X and Y must be equal.", classname);
  endif

  [gY, gnY, glY] = grp2idx (Y);

  ## Keep only the named classes, if any were named.  Names given as text
  ## are matched against grp2idx's own names, which are always a cell array
  ## of character vectors.  A character matrix must be turned into one
  ## first: ismember on two character matrices compares them character by
  ## character and answers a question nobody asked.
  if (! isempty (ClassNames))
    if (iscellstr (ClassNames) || ischar (ClassNames))
      drop = find (! ismember (gnY, cellstr (ClassNames)));
    else
      drop = find (! ismember (glY, ClassNames));
    endif
    for i = 1:numel (drop)
      gY(gY == drop(i)) = NaN;
    endfor
  endif

  ## Weights are validated against the data as supplied, then follow it
  ## through the rows that are kept.
  if (isempty (Weights))
    Weights = ones (rows (X), 1);
  else
    Weights = Weights(:);
    if (numel (Weights) != rows (X))
      error ("%s: 'Weights' must have one element per observation.", ...
             classname);
    endif
  endif

  RowsUsed = ! (isnan (gY) | any (isnan (X), 2));
  ## Index the rows, not the elements: a character matrix of class names has
  ## one row per observation and several columns, and a linear index would
  ## flatten it into single characters.
  F = struct ();
  F.RowsUsed = RowsUsed;
  F.X = X(RowsUsed, :);
  F.Y = Y(RowsUsed, :);
  F.Weights = Weights(RowsUsed);
  if (isempty (F.Y))
    error ("%s: no complete observations in the data.", classname);
  endif

  [gY, gnY, glY] = grp2idx (F.Y);
  nclasses = numel (gnY);
  F.gY = gY;
  F.ClassNames = glY;
  [F.n, F.p] = size (F.X);

  if (nclasses != 2)
    error (strcat ("%s: Y must name exactly two classes, this being a", ...
                   " binary model; use fitcecoc for more than two."), ...
           classname);
  endif

  ## Prior defaults to the weighted frequencies of the training data and
  ## Cost to zero on the diagonal and one elsewhere.
  if (isstruct (Prior))
    Prior = priorFromStruct (Prior, F.ClassNames, classname);
  endif
  freq = accumarray (gY(:), F.Weights(:), [nclasses, 1])' / sum (F.Weights);
  if (isempty (Prior) || (ischar (Prior) && strcmpi (Prior, 'empirical')))
    Prior = freq;
  elseif (ischar (Prior) && strcmpi (Prior, 'uniform'))
    Prior = ones (1, nclasses) / nclasses;
  else
    if (numel (Prior) != nclasses)
      error ("%s: 'Prior' must have one entry per class.", classname);
    endif
    Prior = Prior(:)' / sum (Prior);
  endif
  F.Prior = Prior;

  if (isempty (Cost))
    Cost = ones (nclasses) - eye (nclasses);
  elseif (rows (Cost) != nclasses || columns (Cost) != nclasses)
    error (strcat ("%s: the number of rows and columns in 'Cost' must", ...
                   " correspond to the classes in Y."), classname);
  endif
  F.Cost = Cost;

  ## MATLAB folds the cost matrix into the prior before it weights the
  ## observations: a class that is costlier to misclassify carries more
  ## weight, in proportion to the total cost of getting it wrong.  The Prior
  ## the model reports is the one it was given, not this one.  Measured on
  ## R2024a: a fit with Cost [0 4; 1 0] and one with Prior [0.8 0.2] return
  ## the same coefficients to every digit, and giving both together lands
  ## back on the empirical fit, since [0.2 0.8] scaled by [4 1] is uniform
  ## again.
  adjPrior = Prior .* sum (Cost, 2)';
  if (sum (adjPrior) > 0)
    adjPrior = adjPrior / sum (adjPrior);
  else
    adjPrior = Prior;
  endif

  W = zeros (F.n, 1);
  for k = 1:nclasses
    idx = (gY == k);
    tot = sum (F.Weights(idx));
    if (tot > 0)
      W(idx) = F.Weights(idx) / tot * adjPrior(k);
    endif
  endfor
  if (sum (W) > 0)
    W = W / sum (W);
  endif
  F.W = W;

  ## The positive class is the second of ClassNames, so a positive score
  ## belongs to it.
  F.y = -ones (F.n, 1);
  F.y(gY == 2) = 1;

endfunction
