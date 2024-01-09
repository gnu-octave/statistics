## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{label} =} predict (@var{obj}, @var{XC})
## @deftypefnx {statistics} {[@var{label}, @var{score}, @var{cost}] =} predict (@var{obj}, @var{XC})
##
## Classify new data points into categories using the kNN algorithm from a
## k-Nearest Neighbor classification model.
##
## @code{@var{label} = predict (@var{obj}, @var{XC}} returns the matrix of
## labels predicted for the corresponding instances in @var{XC}, using the
## predictor data in @code{X} and corresponding labels, @code{Y}, stored in the
## k-Nearest Neighbor classification model, @var{obj}.  @var{XC} must be a
## @math{MxP} numeric matrix with the same number of features @math{P} as the
## corresponding predictors of the kNN model in @var{obj}.
##
## @itemize
## @item
## @var{obj} must be a @qcode{ClassificationKNN} object.
## @end itemize
##
## @code{[@var{label}, @var{score}, @var{cost}] = predict (@var{obj}, @var{XC}}
## also returns @var{score}, which contains the predicted class scores or
## posterior probabilities for each instance of the corresponding unique
## classes, and @var{cost}, which is a matrix containing the expected cost of
## the classifications.
##
## @seealso{fitcknn, @@ClassificationKNN/ClassificationKNN}
## @end deftypefn

function [label, score, cost] = predict (obj, XC)

  ## Check for sufficient input arguments
  if (nargin < 2)
    error ("@ClassificationKNN/predict: too few input arguments.");
  endif

  ## Check for valid XC
  if (isempty (XC))
    error ("@ClassificationKNN/predict: XC is empty.");
  elseif (columns (obj.X) != columns (XC))
    error (strcat (["@ClassificationKNN/predict: XC must have the same"], ...
                   [" number of features (columns) as in the kNN model."]));
  endif

  ## Get training data and labels
  X = obj.X(logical (obj.RowsUsed),:);
  Y = obj.Y(logical (obj.RowsUsed),:);

  ## Standardize (if necessary)
  if (obj.Standardize)
    X = (X - obj.Mu) ./ obj.Sigma;
    XC = (XC - obj.Mu) ./ obj.Sigma;
  endif

  ## Train kNN
  if (strcmpi (obj.Distance, "seuclidean"))
    [idx, dist] = knnsearch (X, XC, "k", obj.NumNeighbors, ...
                  "NSMethod", obj.NSMethod, "Distance", "seuclidean", ...
                  "Scale", obj.DistParameter, "sortindices", true, ...
                  "includeties", obj.IncludeTies, "bucketsize", obj.BucketSize);

  elseif (strcmpi (obj.Distance, "mahalanobis"))
    [idx, dist] = knnsearch (X, XC, "k", obj.NumNeighbors, ...
                  "NSMethod", obj.NSMethod, "Distance", "mahalanobis", ...
                  "cov", obj.DistParameter, "sortindices", true, ...
                  "includeties", obj.IncludeTies, "bucketsize", obj.BucketSize);

  elseif (strcmpi (obj.Distance, "minkowski"))
    [idx, dist] = knnsearch (X, XC, "k", obj.NumNeighbors, ...
                  "NSMethod", obj.NSMethod, "Distance", "minkowski", ...
                  "P", obj.DistParameter, "sortindices", true, ...
                  "includeties",obj.IncludeTies, "bucketsize", obj.BucketSize);

  else
    [idx, dist] = knnsearch (X, XC, "k", obj.NumNeighbors, ...
                  "NSMethod", obj.NSMethod, "Distance", obj.Distance, ...
                  "sortindices", true, "includeties", obj.IncludeTies, ...
                  "bucketsize", obj.BucketSize);
  endif

  ## Make prediction
  label = {};
  score = [];
  cost  = [];

  ## Get IDs and labels for each point in training data
  [gY, gnY, glY] = grp2idx (Y);

  ## Evaluate the K nearest neighbours for each new point
  for i = 1:rows (idx)

    ## Get K nearest neighbours
    if (obj.IncludeTies)
      NN_idx = idx{i};
      NNdist = dist{i};
    else
      NN_idx = idx(i,:);
      NNdist = dist(i,:);
    endif
    k = numel (NN_idx);
    kNNgY = gY(NN_idx);

    ## Count frequency for each class
    for c = 1:numel (obj.ClassNames)
      freq(c) = sum (kNNgY == c) / k;
    endfor

    ## Get label according to BreakTies
    if (strcmpi (obj.BreakTies, "smallest"))
      [~, idl] = max (freq);
    else
      idl = find (freq == max (freq));
      tgn = numel (idl);
      if (tgn > 1)
        if (strcmpi (obj.BreakTies, "nearest"))
          for t = 1:tgn
            tgs(t) = find (gY(NN_idx) == idl(t));
          endfor
          [~, idm] = min (tgs);
          idl = idl(idm);
        else      # "random"
          idl = idl(randperm (numel (idl))(1));
        endif
      endif
    endif
    label = [label; gnY{idl}];

    ## Calculate score and cost
    score = [score; freq];
    cost = [cost; 1-freq];

  endfor

endfunction


## Test output
%!shared x, y
%! load fisheriris
%! x = meas;
%! y = species;
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "versicolor"; "virginica"})
%! assert (s, [1, 0, 0; 0, 1, 0; 0, 0, 1])
%! assert (c, [0, 1, 1; 1, 0, 1; 1, 1, 0])
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "mahalanobis");
%! [l, s, c] = predict (obj, xc);
%! assert (s, [0.3, 0.7, 0; 0, 0.9, 0.1; 0.2, 0.2, 0.6], 1e-4)
%! assert (c, [0.7, 0.3, 1; 1, 0.1, 0.9; 0.8, 0.8, 0.4], 1e-4)
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "cosine");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "versicolor"; "virginica"})
%! assert (s, [1, 0, 0; 0, 1, 0; 0, 0.3, 0.7], 1e-4)
%! assert (c, [0, 1, 1; 1, 0, 1; 1, 0.7, 0.3], 1e-4)
%!test
%! xc = [5.2, 4.1, 1.5,	0.1; 5.1,	3.8, 1.9,	0.4; ...
%!         5.1, 3.8, 1.5, 0.3; 4.9, 3.6, 1.4, 0.1];
%! obj = fitcknn (x, y, "NumNeighbors", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "setosa"; "setosa"; "setosa"})
%! assert (s, [1, 0, 0; 1, 0, 0; 1, 0, 0; 1, 0, 0])
%! assert (c, [0, 1, 1; 0, 1, 1; 0, 1, 1; 0, 1, 1])
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.6, 0.4], 1e-4)
%! assert (c, [1, 0.4, 0.6], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "minkowski", "Exponent", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.5, 0.5], 1e-4)
%! assert (c, [1, 0.5, 0.5], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "jaccard");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.9, 0.1, 0], 1e-4)
%! assert (c, [0.1, 0.9, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "mahalanobis");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0.1000, 0.5000, 0.4000], 1e-4)
%! assert (c, [0.9000, 0.5000, 0.6000], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 5, "distance", "jaccard");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.8, 0.2, 0], 1e-4)
%! assert (c, [0.2, 0.8, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 5, "distance", "seuclidean");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 1, 0], 1e-4)
%! assert (c, [1, 0, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "chebychev");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.7, 0.3], 1e-4)
%! assert (c, [1, 0.3, 0.7], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "cityblock");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.6, 0.4], 1e-4)
%! assert (c, [1, 0.4, 0.6], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "cosine");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"virginica"})
%! assert (s, [0, 0.1, 0.9], 1e-4)
%! assert (c, [1, 0.9, 0.1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "correlation");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"virginica"})
%! assert (s, [0, 0.1, 0.9], 1e-4)
%! assert (c, [1, 0.9, 0.1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 30, "distance", "spearman");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 1, 0], 1e-4)
%! assert (c, [1, 0, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 30, "distance", "hamming");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.4333, 0.3333, 0.2333], 1e-4)
%! assert (c, [0.5667, 0.6667, 0.7667], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 5, "distance", "hamming");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.8, 0.2, 0], 1e-4)
%! assert (c, [0.2, 0.8, 1], 1e-4)
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "correlation");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "versicolor"; "virginica"})
%! assert (s, [1, 0, 0; 0, 1, 0; 0, 0.4, 0.6], 1e-4)
%! assert (c, [0, 1, 1; 1, 0, 1; 1, 0.6, 0.4], 1e-4)
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "hamming");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa";"setosa";"setosa"})
%! assert (s, [0.9, 0.1, 0; 1, 0, 0; 0.5, 0, 0.5], 1e-4)
%! assert (c, [0.1, 0.9, 1; 0, 1, 1; 0.5, 1, 0.5], 1e-4)

## Test input validation
%!error<@ClassificationKNN/predict: too few input arguments.> ...
%! predict (ClassificationKNN (ones (4,2), ones (4,1)))
%!error<@ClassificationKNN/predict: XC is empty.> ...
%! predict (ClassificationKNN (ones (4,2), ones (4,1)), [])
%!error<@ClassificationKNN/predict: XC must have the same number of features> ...
%! predict (ClassificationKNN (ones (4,2), ones (4,1)), 1)
