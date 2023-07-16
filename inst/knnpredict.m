## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
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
## @deftypefn  {statistics} {@var{label} =} knnpredict (@var{X}, @var{y}, @var{Xclass})
## @deftypefnx {statistics} {@var{label} =} knnpredict (@var{x}, @var{y}, @var{Xclass}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{label}, @var{score}, @var{cost}] =} knnpredict (@var{x}, @var{y}, @var{Xclass}, @var{name}, @var{value})
##
## Classify new data points into categories using kNN algorithm
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the kNN model.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain any type of
## categorical data. @var{Y} must have same numbers of Rows as @var{X}.
## @item
## @code{Xclass} must be a @math{MxP} numeric matrix of query/new points that
## are to be classified into the labels.
## @var{Xclass} must have same numbers of columns as @var{X}.
##
## @emph{ additional parameters that can be passed as name value pairs :}
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"K"} @tab is the number of nearest neighbors to be found
## in the kNN search.  It must be a positive integer value and by default it is
## 1.
## @item @tab @qcode{"weights"} @tab is a @math{Nx1} numeric non-negative matrix
##          of the observational weights, each row in @var{weights} corresponds
##          to the row in @var{Y} and indicates the relative importance or
##          weight to be considered in calculating the Nearest-neighbour,
##          negative values are removed before calculations if weights are
##          specified. default value @qcode{@var{weight} = ones(rows(Y),1)}.
##
## @item @tab @qcode{"P"} @tab is the Minkowski distance exponent and it must be
##          a positive scalar.  This argument is only valid when the selected
##          distance metric is @qcode{"minkowski"}.  By default it is 2.
##
## @item @tab @qcode{"scale"} @tab is the scale parameter for the standardized
##          Euclidean distance and it must be a nonnegative numeric vector of
##          equal length to the number of columns in @var{X}.  This argument is
##          only valid when the selected distance metric is @qcode{"seuclidean"}
##          , in which case each coordinate of @var{X} is scaled by the
##          corresponding element of @qcode{"scale"}, as is each query point in
##          @var{Y}.  By default, the scale parameter is the standard deviation
##          of each coordinate in @var{X}.
##
## @item @tab @qcode{"cov"} @tab is the covariance matrix for computing the
##          mahalanobis distance and it must be a positive definite matrix
##          matching the the number of columns in @var{X}.  This argument is
##          only valid when the selected distance metric is
##          @qcode{"mahalanobis"}.
## @item @tab @qcode{"cost"} @tab is a @math{NxR} numeric matrix containing
##          misclassification cost for the corresponding instances in @var{X}
##          where @var{R} is the number of unique categories in @var{Y}.
##          If an instance is correctly classified into its category the
##          cost is calculated to be 1, If not then 0. default value
##          @qcode{@var{cost} = ones(rows(X),numel(unique(Y)))}.
##
## @item @tab @qcode{"BucketSize"} @tab is the maximum number of data points in
##          the leaf node of the Kd-tree and it must be a positive integer.
##          This argument is only valid when the selected search
##          method is @qcode{"kdtree"}.
##
## @item @tab @qcode{"Distance"} @tab is the distance metric used by
## @code{knnsearch} as specified below:
## @end multitable
##
## @multitable @columnfractions 0.1 0.25 0.65
## @item @tab @qcode{"euclidean"} @tab Euclidean distance.
## @item @tab @qcode{"seuclidean"} @tab standardized Euclidean distance.  Each
## coordinate difference between the rows in @var{X} and the query matrix
## @var{Y} is scaled by dividing by the corresponding element of the standard
## deviation computed from @var{X}.  To specify a different scaling, use the
## @qcode{"scale"} name-value argument.
## @item @tab @qcode{"cityblock"} @tab City block distance.
## @item @tab @qcode{"chebychev"} @tab Chebychev distance (maximum coordinate
## difference).
## @item @tab @qcode{"minkowski"} @tab Minkowski distance.  The default exponent
## is 2.  To specify a different exponent, use the @qcode{"P"} name-value
## argument.
## @item @tab @qcode{"mahalanobis"} @tab Mahalanobis distance, computed using a
## positive definite covariance matrix.  To change the value of the covariance
## matrix, use the @qcode{"cov"} name-value argument.
## @item @tab @qcode{"cosine"} @tab Cosine distance.
## @item @tab @qcode{"correlation"} @tab One minus the sample linear correlation
## between observations (treated as sequences of values).
## @item @tab @qcode{"spearman"} @tab One minus the sample Spearman's rank
## correlation between observations (treated as sequences of values).
## @item @tab @qcode{"hamming"} @tab Hamming distance, which is the percentage
## of coordinates that differ.
## @item @tab @qcode{"jaccard"} @tab One minus the Jaccard coefficient, which is
## the percentage of nonzero coordinates that differ.
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"NSMethod"} @tab is the nearest neighbor search method used
## by @code{knnsearch} as specified below.
## @end multitable
##
## @multitable @columnfractions 0.1 0.25 0.65
## @item @tab @qcode{"kdtree"} @tab Creates and uses a Kd-tree to find nearest
## neighbors.  @qcode{"kdtree"} is the default value when the number of columns
## in @var{X} is less than or equal to 10, @var{X} is not sparse, and the
## distance metric is @qcode{"euclidean"}, @qcode{"cityblock"},
## @qcode{"chebychev"}, or @qcode{"minkowski"}.  Otherwise, the default value is
## @qcode{"exhaustive"}.  This argument is only valid when the distance metric
## is one of the four aforementioned metrics.
## @item @tab @qcode{"exhaustive"} @tab Uses the exhaustive search algorithm by
## computing the distance values from all the points in @var{X} to each point in
## @var{Y}.
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"standardize"} @tab is the flag to indicate if kNN should
##               be calculated standerdizing @var{X}.
## @end multitable
##
## @code{@var{label} = knnpredict (@var{X}, @var{Y}, @var{Xclass})}
## returns the matrix of labels predicted for the corresponding instances
## in @code{Xclass}, using the predictor data in @code{X} and corresponding
## categorical data in @code{Y}. @var{X} is used to train the kNN model and
## values in @var{Xclass} are classified into classes in @var{Y}.
##
## @code{@var{label} = knnpredict (@var{X}, @var{Y}, @var{Xclass}, @var{Name}, @var{Value})}
## returns a matrix @var{label} containing the predicted labels for instances in
## @var{Xclass} with additional parameters given as Name-Value pairs. for example
## value of nearest-neighbours can be set 10 by specifying "K" as 10, similiarly
## other optional parameters can be given.
##
## @code{[@var{label}, @var{score}, @var{cost}] = knnpredict (@dots{})}
## returns the matrix of labels predicted for the corresponding instances
## in @code{Xclass}, using the predictor data in @code{X} and @code{Y}.
##
## In addition it also returns
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"score"} @tab contains predicted class scores or posterior
##             Probabilities for each instances for corresponding unique
##             classes in @var{Y}.
## @item @tab @qcode{"cost"} @tab is a matrix containing expected cost of the
##             classifications. each row of @var{cost} matrix contains the
##             expected cost of classification of observations in @var{Xclass}
##             into each class of unique classes in @var{Y}.
##@end multitable
##
## @end itemize
##
## for demo use demo knnpredict
##
## @seealso{knnsearch}
## @end deftypefn

function [label, score, cost] = knnpredict (X, Y, Xclass, varargin)


  ## check positional parameters
  if (nargin < 3)
	  error ("knnpredict: too few input arguments.");
  endif

  if (rows (X) != rows (Y))
    error ("knnpredict: number of rows in X and Y must be equal.");
  endif

  if (columns (X) != columns (Xclass))
    error ("knnpredict: number of columns in Xclass must be equal to X.");
  endif

  ##check X
  ## remove this do this after removing NaNs
  if (! isnumeric (X) || ! isfinite (X))
    error ("knnpredict: Invalid values in X.");
  endif


  ## process optional parameters

  ## adding default values
    k = 1;
    cost  = [];
    S = [];
    C   = [];
    weights  = [];
    P = 2;
    BS  = 50;
    distance    = "euclidean";
    NSmethod    = "exhaustive";
    standardize = false;


  ## Parse additional parameters in Name/Value pairs
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))
      case "k"
        k = varargin{2};
      case "weights"
        weights = varargin{2};
      case "p"
        P = varargin{2};
      case "cost"
        cost = varargin{2};
      case "scale"
        S = varargin{2};
      case "cov"
        C = varargin{2};
      case "bucketsize"
        BS = varargin{2};
      case "distance"
        distance = varargin{2};
      case "nsmethod"
        NSmethod = varargin{2};
      case "standardize"
        standardize = varargin{2};
      otherwise
        error ("knnsearch: invalid NAME in optional pairs of arguments.");
    endswitch
    varargin (1:2) = [];
  endwhile


  ##------checking optional parameters------##
  ## check k
  if (! isscalar (k) || !isnumeric (k) || k < 1 || k != round (k))
    error ("knnpredict: Invalid value of k.");
  endif

  ## check weights
  if (! isempty (weights))
    if (! isvector (weights) || rows (X) != length (weights) ...
                             || sum(weights < 0) != 0)
        error("knnpredict: Weights has invalid observarions.");
    endif
  endif

  ## check minkowski distance exponent
  if (! isscalar (P) || ! isnumeric (P) || P < 0)
    error ("knnpredict: Invalid value of Minkowski Exponent.");
  endif

  ## check cost
  if (! isempty (cost))
    if (! isscalar (cost) || ! isnumeric(cost) || ...
        columns != numel (unique(Y)) || rows (cost) != rows (X) )
        error ("knnpredict: Invalid value in cost or the size of cost.");
    endif
  else
    ## empty cost
    cost = ones (rows (Xclass),numel (unique (Y)));
  endif

  ## check scale
  if (! isempty (S))
    if (! isscalar (S) || any (S) < 0 || numel (S) != rows (X))
      error ("knnpredict: Invalid value in Scale or the size of scale.");
    endif
  endif

  ## check C
  if (! isempty (C))
    if (! strcmp (distance,"mahalanobis") || ! ismatrix (C) ...
                                          || ! isnumeric (C))
      error (strcat (["knnpredict: Invalid value in cov, cov can only be"], ...
                    [" given for mahalanobis distance."]));
    endif
  endif

  ## check BS
  if (! isscalar (BS) || BS < 0)
    error ("knnpredict: Invalid value of bucketsize.");
  endif

  ##standardize
  if (! islogical(standardize) && standardize != 1 && standardize != 0)
    error ("knnpredict: value of standardize must be boolean.");
  endif


  ##------checking optional params end------##

  if (standardize)
    stdX = std (X,0,1);
    mu   = mean(X);
    X    = (X - mu) ./ stdX;
  endif


  ## calculate distance for each point in Xclass and pass to function
  ## to predict label for the data point
  classNames = unique (Y);
  classNos   = numel (classNames);


  ## main kNN algorithm

  if (isempty (X))
    ## no data in X
    label     = repmat (classNames(1,:),0,1);
    posterior = NaN (0,classNos);
    cost      = NaN (0,classNos);
  else
    ## calculate the NNs using knnsearch
    [idx, dist] = knnsearch (X, Xclass, "k", k, "NSMethod", NSmethod, ...
                  "Distance", distance, "P", P, "Scale", S, ...
                  "cov", C, "bucketsize", BS, "sortindices", true);

    [label, score, cost_temp] = predictlabel (X, Y, idx, k, weights);

    cost  = cost(rows(cost_temp),columns(cost_temp)) .* cost_temp;
  endif

endfunction


##-------function labels-------##
function [labels, score, cost_temp] = predictlabel (x, y, idx, k, weights)

  ## assign intial values
  u = unique (y);
  freq = [];
  wsum  = sum (weights);
  score = [];
  labels = [];
  cost_temp = [];
  classNos = size (u,1);

  # [iu, id] = meshgrid (1:size (u, 1), 1:size (idx,1));
  for i = 1:rows (idx)

    if (!isempty (weights))

      ## weighted kNN
      for id = 1:numel (u)
        freq = [freq; (sum (strcmpi (u(id,1), y(idx(i,:))) .* weights)) / wsum;]
        score_temp = (freq ./ k)';
      endfor
    else

      ## Non-weighted kNN
      for id = 1:size (u,1) #u{iu(:),2}
        freq(id,1) = (sum (strcmpi (u(id,1), y(idx(i,:)))));
      endfor

      ## score calculation
      score_temp = (freq ./ k)';
      cost_temp  = [cost_temp; ones(1,classNos) - score_temp];
      score = [score; score_temp];

    endif

    [val, iu] = max (freq);

    ## set label for the index idx
    labels = [labels; u(iu,1)];
  endfor
endfunction
## ------labels------- ##

%!demo
%! ## find 10 nearest neighbour of a point using different distance metrics
%! ## and compare the results by plotting
%! load fisheriris
%! x = meas(:,3:4);
%! y = species;
%!
%! point = [5, 1.45];
%!
%! ## calculate 10 nearest-neighbours by minkowski distance
%! [id, d] = knnsearch (x, point, "K", 10);
%! ld = knnpredict (x, y, point, "K", 10);
%!
%! ## calculate 10 nearest-neighbours by minkowski distance
%! [idm, dm] = knnsearch (x, point, "K", 10, "distance", "minkowski", "p", 5);
%! lm = knnpredict (x, y, point, "K", 10, "distance", "minkowski", "p", 5);
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idc, dc] = knnsearch (x, point, "K", 10, "distance", "chebychev");
%! lc = knnpredict (x, y, point, "K", 10, "distance", "chebychev");
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [ids, ds] = knnsearch (x, point, "K", 10, "distance", "hamming");
%! ls = knnpredict (x, y, point, "K", 10, "distance", "hamming");
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idb, db] = knnsearch (x, point, "K", 10, "distance", "cityblock");
%! lb = knnpredict (x, y, point, "K", 10, "distance", "cityblock");
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idh, dh] = knnsearch (x, point, "K", 10, "distance", "manhattan");
%! lh = knnpredict (x, y, point, "K", 10, "distance", "manhattan");
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idn, dn] = knnsearch (x, point, "K", 10, "distance", "cosine");
%! ln = knnpredict (x, y, point, "K", 10, "distance", "cosine");
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idj, dj] = knnsearch (x, point, "K", 10, "distance", "jaccard");
%! lj = knnpredict (x, y, point, "K", 10, "distance", "jaccard");
%!
%! fprintf (strcat ([" Labels predicted by different distance metrics  "], ...
%!                 [" for 10 Nearest Neighbours of [5, 1.45] \n"]));
%! fprintf (strcat (["\n Euclidean : "],[ld{1}]));
%! fprintf (strcat (["\n Minkowski : "],[lm{1}]));
%! fprintf (strcat (["\n Chebychev : "],[lc{1}]));
%! fprintf (strcat (["\n Hamming   : "],[ls{1}]));
%! fprintf (strcat (["\n cityblock : "],[lb{1}]));
%! fprintf (strcat (["\n Manhattan : "],[lh{1}]));
%! fprintf (strcat (["\n Cosine    : "],[ln{1}]));
%! fprintf (strcat (["\n Jaccard   : "],[lj{1}]));
%! warning("off");
%! hold on
%! s1 = subplot (2, 2, 1);
%! gscatter(x(:,1), x(:,2), species,[.75 .75 0; 0 .75 .75; .75 0 .75], '.',20);
%! set( s1, 'title', 'euclidean and minkowski distances' );
%! xlabel('Petal length (cm)');
%! ylabel('Petal width (cm)');
%! line(point(1), point(2), 'marker', 'x', 'color', 'k', 'linewidth', 2, 'displayname', 'query point')
%! line (x(id,1), x(id,2), 'color', [0.5 0.5 0.5],'marker', 'o', 'linestyle', 'none', 'markersize', 10, "displayname", "eulcidean")
%! line (x(idm,1), x(idm,2), 'color', [0.5 0.5 0.5],'marker', 'd', 'linestyle', 'none', 'markersize', 10, "displayname", "Minkowski")
%! xlim ([4.5 5.5]);
%! ylim ([1 2]);
%! axis square;
%!
%! s2 = subplot (2, 2, 2);
%! gscatter(x(:,1), x(:,2), species,[.75 .75 0; 0 .75 .75; .75 0 .75], '.',20);
%! set( s2, 'title', 'chebychev  and hamming distance' );
%! xlabel('Petal length (cm)');
%! ylabel('Petal width (cm)');
%! line(point(1), point(2), 'marker', 'x', 'color', 'k', 'linewidth', 2, 'displayname', 'query point')
%! line (x(idc,1), x(idc,2), 'color', [0.5 0.5 0.5],'marker', 'p', 'linestyle', 'none', 'markersize', 10, "displayname", "chebychev")
%! line (x(ids,1), x(ids,2), 'color', [0.5 0.5 0.5],'marker', 's', 'linestyle', 'none', 'markersize', 10, "displayname", "hamming")
%! xlim ([4.5 5.5]);
%! ylim ([1 2]);
%! axis square;
%!
%! s3 = subplot (2, 2, 3);
%! gscatter(x(:,1), x(:,2), species,[.75 .75 0; 0 .75 .75; .75 0 .75], '.',20);
%! set( s3, 'title', 'cityblock and manhattan distance' );
%! xlabel('Petal length (cm)');
%! ylabel('Petal width (cm)');
%! line(point(1), point(2), 'marker', 'x', 'color', 'k', 'linewidth', 2, 'displayname', 'query point')
%! line (x(idb,1), x(idb,2), 'color', [0.5 0.5 0.5],'marker', '>', 'linestyle', 'none', 'markersize', 10, "displayname", "cityblock")
%! line (x(idh,1), x(idh,2), 'color', [0.5 0.5 0.5],'marker', 'd', 'linestyle', 'none', 'markersize', 10, "displayname", "manhattan")
%! xlim ([4.5 5.5]);
%! ylim ([1 2]);
%! axis square;
%!
%! s4 = subplot (2, 2, 4);
%! gscatter(x(:,1), x(:,2), species,[.75 .75 0; 0 .75 .75; .75 0 .75], '.',20);
%! set( s4, 'title', 'cosine and jaccard distance' );
%! xlabel('Petal length (cm)');
%! ylabel('Petal width (cm)');
%! line(point(1), point(2), 'marker', 'x', 'color', 'k', 'linewidth', 2, 'displayname', 'query point')
%! line (x(idn,1), x(idn,2), 'color', [0.5 0.5 0.5],'marker', '<', 'linestyle', 'none', 'markersize', 10, "displayname", "cosine")
%! line (x(idj,1), x(idj,2), 'color', [0.5 0.5 0.5],'marker', 'h', 'linestyle', 'none', 'markersize', 10, "displayname", "jaccard")
%! xlim ([4.5 5.5]);
%! ylim ([1 2]);
%! axis square;
%! hold off
%!
%!
## Test output
%!shared x, y
%! load fisheriris
%! x = meas;
%! y = species;
%!test
%! xnew = [min(x);mean(x);max(x)];
%! [l, s, c] = knnpredict (x, y, xnew, "K", 5);
%! assert (l, {"setosa";"versicolor";"virginica"})
%! assert (s, [1, 0, 0;0, 1, 0;0, 0, 1])
%! assert (c, [0, 1, 1;1, 0, 1;1, 1, 0])
%!test
%! xnew = [min(x);mean(x);max(x)];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance", "mahalanobis");
%! assert (l, {"versicolor";"versicolor";"virginica"})
%! assert (s, [0.3000, 0.7000, 0;0, 0.9000, 0.1000;0.2000, 0.2000, 0.6000], 1e-4)
%! assert (c, [0.7000, 0.3000, 1.0000;1.0000, 0.1000, 0.9000;0.8000, 0.8000, 0.4000], 1e-4)
%!test
%! xnew = [min(x);mean(x);max(x)];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance", "cosine");
%! assert (l, {"setosa";"versicolor";"virginica"})
%! assert (s, [1.0000, 0, 0;0, 1.0000, 0;0, 0.3000, 0.7000], 1e-4)
%! assert (c, [0, 1.0000, 1.0000;1.0000, 0, 1.0000;1.0000, 0.7000, 0.3000], 1e-4)
%!test
%! xnew = [5.2, 4.1, 1.5,	0.1;5.1,	3.8,	1.9,	0.4;5.1,	3.8, 1.5,	0.3;4.9,	3.6,	1.4,	0.1];
%! [l, s, c] = knnpredict (x, y, xnew, "K", 5);
%! assert (l, {"setosa";"setosa";"setosa";"setosa"})
%! assert (s, [1, 0, 0;1, 0, 0;1, 0, 0;1, 0, 0])
%! assert (c, [0, 1, 1;0, 1, 1;0, 1, 1;0, 1, 1])
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 5);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.6000, 0.4000], 1e-4)
%! assert (c, [1.0000, 0.4000, 0.6000], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance", "minkowski", "P", 5);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.5000, 0.5000], 1e-4)
%! assert (c, [1.0000, 0.5000, 0.5000])
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance", "jaccard");
%! assert (l, {"setosa"})
%! assert (s, [0.9000, 0.1000, 0], 1e-4)
%! assert (c, [0.1000, 0.9000, 1.0000], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance", "mahalanobis");
%! assert (l, {"versicolor"})
%! assert (s, [0.1000, 0.5000, 0.4000], 1e-4)
%! assert (c, [0.9000, 0.5000, 0.6000], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 5, "distance","jaccard");
%! assert (l, {"setosa"})
%! assert (s, [0.8000, 0.2000, 0], 1e-4)
%! assert (c, [0.2000, 0.8000, 1.000], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 5, "distance","seuclidean");
%! assert (l, {"versicolor"})
%! assert (s, [0, 1, 0], 1e-4)
%! assert (c, [1, 0, 1], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance","chebychev");
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.7000, 0.3000], 1e-4)
%! assert (c, [1.000, 0.3000, 0.7000], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance","cityblock");
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.6000, 0.4000], 1e-4)
%! assert (c, [1.000, 0.4000, 0.6000], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance","manhattan");
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.6000, 0.4000], 1e-4)
%! assert (c, [1.000, 0.4000, 0.6000], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance","cosine");
%! assert (l, {"virginica"})
%! assert (s, [0, 0.1000, 0.9000], 1e-4)
%! assert (c, [1.000, 0.9000, 0.1000], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance","correlation");
%! assert (l, {"virginica"})
%! assert (s, [0, 0.1000, 0.9000], 1e-4)
%! assert (c, [1.000, 0.9000, 0.1000], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 30, "distance","spearman");
%! assert (l, {"versicolor"})
%! assert (s, [0, 1, 0], 1e-4)
%! assert (c, [1, 0, 1], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 30, "distance","hamming");
%! assert (l, {"setosa"})
%! assert (s, [0.4333, 0.3333, 0.2333], 1e-4)
%! assert (c, [0.5667, 0.6667, 0.7667], 1e-4)
%!test
%! xnew = [5, 3, 5, 1.45];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 5, "distance","hamming");
%! assert (l, {"setosa"})
%! assert (s, [0.8000, 0.2000, 0], 1e-4)
%! assert (c, [0.2000, 0.8000, 1.0000], 1e-4)
%!test
%! xnew = [min(x);mean(x);max(x)];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance", "correlation");
%! assert (l, {"setosa";"versicolor";"virginica"})
%! assert (s, [1.0000, 0, 0;0, 1.0000, 0;0, 0.4000, 0.6000], 1e-4)
%! assert (c, [0, 1.0000, 1.0000;1.0000, 0, 1.0000;1.0000, 0.6000, 0.4000], 1e-4)
%!test
%! xnew = [min(x);mean(x);max(x)];
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance", "hamming");
%! assert (l, {"setosa";"setosa";"setosa"})
%! assert (s, [0.9000, 0.1000, 0;1.000, 0, 0;0.5000, 0, 0.5000], 1e-4)
%! assert (c, [0.1000, 0.9000, 1.0000;0, 1.0000, 1.0000;0.5000, 1.0000, 0.5000], 1e-4)


## Test input validation
%!error<knnpredict: too few input arguments.>knnpredict (1)
%!error<knnpredict: too few input arguments.>knnpredict (1, 2)
%!error<knnpredict: number of rows in X and Y must be equal.> ...
%! knnpredict (ones(4,5), ones(5,4),ones(1,1))
%!error<knnpredict: number of columns in Xclass must be equal to X.> ...
%! knnpredict (ones(4,5), ones(4,1),ones(1,1))
%!error<knnpredict: Invalid values in X.> ...
%! knnpredict ([1,2,3;"a",5,6], ones(2,3),ones(1,3))
%!error<knnpredict: Invalid value of k.> ...
%! knnpredict (ones(4,5),ones(4,1),ones(1,5), "k", -5)
%!error<knnpredict: Weights has invalid observarions.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "weights", [1;2;-3;4])
%!error<knnpredict: Weights has invalid observarions.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "weights", ones(2,1))
%!error<knnpredict: Invalid value of Minkowski Exponent.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "P", -3)
%!error<knnpredict: Invalid value in cost or the size of cost.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "cost", [1,2,3;4,5,6;7,8,9])
%!error<knnpredict: Invalid value in cost or the size of cost.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "cost", ones(4,2))
%!error<knnpredict: Invalid value in Scale or the size of scale.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "Scale", ones(4,2))
%!error<knnpredict: Invalid value in Scale or the size of scale.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "scale", [1;2;3;-4])
%!error<knnpredict: Invalid value in cov, cov can only be given for mahalanobis distance.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "cov", ones(4,2),"distance","euclidean")
%!error<knnpredict: Invalid value of bucketsize.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "bucketsize",-50)
%!error<knnpredict: value of standardize must be boolean.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "standardize","some")
