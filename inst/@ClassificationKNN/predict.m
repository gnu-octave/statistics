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
## k-Nearest Neighbor classification model, @var{obj}.
##
## @itemize
## @item
## @code{obj} must be a @qcode{ClassificationKNN} class object.
## @end itemize
##
## @code{@var{label} = predict (@var{obj}, @var{XC}} returns the matrix of
## labels predicted for the corresponding instances in @var{XC}, using the
## predictor data in @code{XC} and corresponding labels, @code{Y}, stored in the
## k-Nearest Neighbor classification model, @var{obj}.  @var{XC} must be a
## @math{MxP} numeric matrix with the same number of features @math{P} as the
## corresponding kNN model in @var{obj}.
##
## @code{[@var{label}, @var{score}, @var{cost}] = predict (@var{obj}, @var{XC}}
## also returns @var{score}, which contains the predicted class scores or
## posterior probabilities for each instance of the corresponding unique
## classes, and @var{cost}, which is a matrix containing the expected cost of
## the classifications.
##
## @seealso{fitcknn, knnpredict}
## @end deftypefn

function [label, score, cost] = predict (obj, XC)

  ## Check for obj being a ClasifficationKNN object
  if (! strcmp (class (obj), "ClasifficationKNN"))
    error (strcat (["@ClassificationKNN/predict: OBJ"], ...
                   [" must a ClasifficationKNN class object."]));
  endif
  ## Check for valid XC
  if (isempty (XC))
    error ("@ClassificationKNN/predict: X is empty.");
  elseif (columns (obj.X) != columns (XC))
    error (strcat (["@ClassificationKNN/predict: XC must have the"], ...
                   [" same number of features as in the kNN model."]));
  endif

  ## Check cost
  if (isempty (obj.cost))
    ## if empty assign all cost = 1
    obj.cost = ones (rows (obj.X), obj.NosClasses);
  endif

  if (isempty (obj.X))
    ## No data in X
    label     = repmat(classNames(1,:),0,1);
    posterior = NaN(0,classNos);
    cost      = NaN(0,classNos);
  else
    ## Calculate the NNs using knnsearch
    [idx, dist] = knnsearch (obj.X, XC, "k", obj.k, ...
                 "NSMethod", obj.NSmethod, "Distance", obj.distance, ...
                 "P", obj.P, "Scale", obj.Scale, "cov", obj.cov, ...
                 "bucketsize", obj.bucketsize, "sortindices", true, ...
                 "includeties",obj.Includeties);

    [label, score, cost_temp] = predictlabel (obj, idx);
    cost  = obj.cost(rows(cost_temp),columns(cost_temp)) .* cost_temp;

    ## Store predicted in the object variables
    obj.NN    = idx;
    obj.label = label;
    obj.score = score;
    obj.cost  = cost;
  endif

endfunction

## Helper function to predict labels
function [labels, score, cost_temp] = predictlabel (obj, idx)
  ## Assign intial values
  freq = [];
  wsum  = sum (obj.weights);
  score = [];
  labels = [];
  cost_temp = [];

  for i = 1:rows (idx)
    if (!isempty (obj.weights))
      ## Weighted kNN
      for id = 1:numel (obj.classNames)
        new_freq = sum (strcmpi (obj.classNames(id,1), obj.Y(idx(i,:)))  .* ...
                        obj.weights) / wsum;
        freq = [freq; new_freq];
        score_tmp = (freq ./ obj.k)';
      endfor
    else
      ## Non-weighted kNN
      for id = 1:size (obj.classNames,1) #u{iu(:),2}
        freq(id,1) = (sum (strcmpi (obj.classNames(id,1), obj.Y(idx(i,:)))));
      endfor
      ## Score calculation
      score_tmp = (freq ./ obj.k)';
      cost_temp  = [cost_temp; ones(1,obj.NosClasses) - score_tmp];
      score = [score; score_tmp];
    endif

    [val, iu] = max (freq);

    ## Set labels for the index idx
    labels = [labels; obj.classNames(iu,1)];

  endfor

endfunction

%!demo
%! ## find 10 nearest neighbour of a point using different distance metrics
%! ## and compare the results by plotting
%!
%! load fisheriris
%! x = meas;
%! y = species;
%! xc = [5, 3, 5, 1.45];
%!
%! ## Create an object
%! a = fitcknn (x, y, "k", 5)
%!
%! ## Predict labels for points in xc
%! predict (a, xc)
%!
%! ## Change properties keeping training data and predict again
%! a.distance = "hamming";
%! a.k = 10;
%! predict (a, xc)

## Test output
%!shared x, y
%! load fisheriris
%! x = meas;
%! y = species;
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj  = ClassificationKNN (x, y, "k", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "versicolor"; "virginica"})
%! assert (s, [1, 0, 0; 0, 1, 0; 0, 0, 1])
%! assert (c, [0, 1, 1; 1, 0, 1; 1, 1, 0])
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "mahalanobis");
%! [l, s, c] = predict (obj, xc);
%! assert (s, [0.3, 0.7, 0; 0, 0.9, 0.1; 0.2, 0.2, 0.6], 1e-4)
%! assert (c, [0.7, 0.3, 1; 1, 0.1, 0.9; 0.8, 0.8, 0.4], 1e-4)
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "cosine");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "versicolor"; "virginica"})
%! assert (s, [1, 0, 0; 0, 1, 0; 0, 0.3, 0.7], 1e-4)
%! assert (c, [0, 1, 1; 1, 0, 1; 1, 0.7, 0.3], 1e-4)
%!test
%! xc = [5.2, 4.1, 1.5,	0.1; 5.1,	3.8, 1.9,	0.4; ...
%!         5.1, 3.8, 1.5, 0.3; 4.9, 3.6, 1.4, 0.1];
%! obj  = ClassificationKNN (x, y, "k", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "setosa"; "setosa"; "setosa"})
%! assert (s, [1, 0, 0; 1, 0, 0; 1, 0, 0; 1, 0, 0])
%! assert (c, [0, 1, 1; 0, 1, 1; 0, 1, 1; 0, 1, 1])
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.6, 0.4], 1e-4)
%! assert (c, [1, 0.4, 0.6], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "minkowski", "P", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.5, 0.5], 1e-4)
%! assert (c, [1, 0.5, 0.5], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "jaccard");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.9, 0.1, 0], 1e-4)
%! assert (c, [0.1, 0.9, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "mahalanobis");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0.1000, 0.5000, 0.4000], 1e-4)
%! assert (c, [0.9000, 0.5000, 0.6000], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 5, "distance", "jaccard");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.8, 0.2, 0], 1e-4)
%! assert (c, [0.2, 0.8, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 5, "distance", "seuclidean");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 1, 0], 1e-4)
%! assert (c, [1, 0, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "chebychev");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.7, 0.3], 1e-4)
%! assert (c, [1, 0.3, 0.7], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "cityblock");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.6, 0.4], 1e-4)
%! assert (c, [1, 0.4, 0.6], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "manhattan");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.6, 0.4], 1e-4)
%! assert (c, [1, 0.4, 0.6], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "cosine");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"virginica"})
%! assert (s, [0, 0.1, 0.9], 1e-4)
%! assert (c, [1, 0.9, 0.1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "correlation");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"virginica"})
%! assert (s, [0, 0.1, 0.9], 1e-4)
%! assert (c, [1, 0.9, 0.1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 30, "distance", "spearman");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 1, 0], 1e-4)
%! assert (c, [1, 0, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 30, "distance", "hamming");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.4333, 0.3333, 0.2333], 1e-4)
%! assert (c, [0.5667, 0.6667, 0.7667], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj  = ClassificationKNN (x, y, "k", 5, "distance", "hamming");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.8, 0.2, 0], 1e-4)
%! assert (c, [0.2, 0.8, 1], 1e-4)
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "correlation");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "versicolor"; "virginica"})
%! assert (s, [1, 0, 0; 0, 1, 0; 0, 0.4, 0.6], 1e-4)
%! assert (c, [0, 1, 1; 1, 0, 1; 1, 0.6, 0.4], 1e-4)
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj  = ClassificationKNN (x, y, "k", 10, "distance", "hamming");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa";"setosa";"setosa"})
%! assert (s, [0.9, 0.1, 0; 1, 0, 0; 0.5, 0, 0.5], 1e-4)
%! assert (c, [0.1, 0.9, 1; 0, 1, 1; 0.5, 1, 0.5], 1e-4)
