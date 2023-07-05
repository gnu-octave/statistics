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
## Classify new data points into categries using kNN algorithm
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows correspond
## to observations and columns correspond to features or variables.@var{X} will
## be used to train the kNN model.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain any type of
## categorical data. @var{Y} must have same numbers of Rows as @var{X}.
## @item
## @code{Xclass} must be a @math{MxP} numeric matrix of query/new points that are
## to be classified into the labels. @var{Xclass} must have same numbers of columns
## as @var{X}.
##
## @emph{ additional parameters that can be passed as name value pairs :}
## @item
## @code{K} is the number of nearest neighbours to be considered in kNN search
##          Default value of @qcode{@var{k} = 1}.
## @item
## @code{weights} is a @math{Nx1} numeric non-negative matrix of the observational
##          weights, each row in @var{weights} corresponds to the row in @var{Y}
##          and indicates the relative importance or weight to be considered in
##          calculating the Nearest-neighbour, negative values are removed before
##          calculations if weights are specified.
##          default value @qcode{@var{weight} = ones(rows(Y),1)}.
## @item
## @code{exponent} is the minkowski distance exponent. Default is @qcode{@var{P} = 2}.
## @item
## @code{cost} is a @math{NxR} numeric matrix containing misclassification cost
##             for the corresponding instances in @var{X} where @var{R} is the
##             number of unique categories in @var{Y}. If an instance is correctly
##             classified into its category the cost is calculated to be 1, If not
##             then 0. default value @qcode{@var{cost} = ones(rows(X),numel(unique(Y)))}.
## @item
## @code{scale} is scale for calculating standardized euclidean distance.
##             Default is @qcode{@var{scale} = []}.
## @item
## @code{cov}   is the covariance matrix for computing mahalanobis distance.
## @item
## @code{bucketsize} is maximum number of data points per leaf node of kd-tree.
##              if NSmethod is 'kdtree', Default is @qcode{@var{bucketsize} = 50}.
## @item
## @code{distance}   is the distance metric to be used in calculating the distance
##                   between points. the choice for distance metric are :
## @itemize @minus
## @item
## @strong{'euclidean'}   - Euclidean distance. this is default.
## @item
## @strong{'sqeuclidean'} - squared euclidean distance.
## @item
## @strong{'cityblock'}   - City Block distance.
## @item
## @strong{'chebyshev'}   - Chebyshev distance.
## @item
## @strong{'minkowski'}   - Minkowski distance with exponent @code{exponent}
##                          default is @qcode{@var{P} = 2}.
## @item
## @strong{'mahalanobis'} - Mahalanobis distance calculated using covariance matrix @code{cov}.
## @item
## @strong{'cosine'}      - Cosine distance.
## @item
## @strong{'correlation'} - Correlation distance
## @item
## @strong{'spearman'}    - spearman distance
## @item
## @strong{'jaccard'}     - jaccard distance.
## @item
## @strong{'hamming'}     - Hamming distance.
## @end itemize
##
## @item
## @code{NSMethod}   is nearest neighbour search method. Choices for @var{NSMethod} are:
## @itemize @minus
## @item
## @strong{'exhaustive'} - function computes the distance of each point in @var{Xclass}
##                         from the predictor data @var{X} to find nearest neighbours.
## @item
## @strong{'kdtree'}     - function builds a kdtree of predictor data @var{x}
##                         and then searches for query points in @var{Xclass}.
##                         kdtree search method works best when @var{X} is not sparse
##                         has less dimensions or columns in @var{X} and large
##                         number of observations.
## @end itemize
## default value of @qcode{@var{NSMethod} = 'exhaustive'}.
##
## @item
## @code{standerdize} is the flag to indicate if kNN should be calculated by
##                    standerdizing @var{X}.
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
## In addition it also returns
## @itemize @minus
## @item
## @var{score} contains predicted class scores or posterior Probabilities
##             for each instances for corresponding unique classes in @var{Y}.
## @item
## @var{cost}  is a matrix containing expected cost of the classifications.
##             each row of @var{cost} matrix contains the expected cost of
##             classification of observations in @var{Xclass} into each class of
##             unique classes in @var{Y}.
##
##
## @end itemize
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

  if (rows(X) != rows(Y))
    error ("knnpredict: number of rows in X and Y must be equal.");
  endif

  if (columns(X) != columns(Xclass))
    error ("knnpredict: number of columns in Xclass must be equal to X.");
  endif

  ##check X
  ## remove this do this after removing NaNs
  if ( !isnumeric(X) || !isfinite(X))
    error ("knnpredict: Invalid values in X.");
  endif


  ## process optional parameters

  ## adding default values
    k = 1;
    cost  = [];
    Scale = [];
    cov   = [];
    weights  = [];
    exponent = 2;
    bucketsize  = 50;
    distance    = "euclidean";
    NSmethod    = "exhaustive";
    standerdize = false;


  ## Parse additional parameters in Name/Value pairs
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))
      case "k"
        k = varargin{2};
      case "weights"
        weights = varargin{2};
      case "exponent"
        exponent = varargin{2};
      case "cost"
        cost = varargin{2};
      case "scale"
        Scale = varargin{2};
      case "cov"
        cov = varargin{2};
      case "bucketsize"
        bucketsize = varargin{2};
      case "distance"
        distance = varargin{2};
      case "nsmethod"
        NSmethod = varargin{2};
      case "standerdize"
        standerdize = varargin{2};
      otherwise
        error ("knnsearch: invalid NAME in optional pairs of arguments.");
    endswitch
    varargin(1:2) = [];
  endwhile


  ##------checking optional parameters------##
  ## check k
  if ( !isscalar(k) || !isnumeric(k) || k < 1 || k != round (k))
    error ("knnpredict: Invalid value of k.");
  endif

  ## check weights
  if ( !isempty(weights))
    if ( !isvector(weights) || rows(X) != length(weights) ...
            || sum(weights < 0) != 0)
        error("knnpredict: Weights has invalid observarions.");
    endif
  endif

  ## check minkowski distance exponent
  if ( !isscalar(exponent) || !isnumeric(exponent) || exponent < 0)
    error ("knnpredict: Invalid value of Minkowski Exponent.");
  endif

  ## check cost
  if ( !isempty(cost))
    if ( !isscalar(cost) || !isnumeric(cost) || ...
        columns != numel(unique(Y)) || rows(cost) != rows(X) )
        error ("knnpredict: Invalid value in cost or the size of cost.");
    endif
  else
    ## empty cost
    cost = ones(rows(X),numel(unique(Y)));
  endif

  ## check scale
  if ( !isempty(Scale))
    if ( !isscalar(Scale) || any(Scale) < 0 || numel(Scale) != rows(X))
      error ("knnpredict: Invalid value in Scale or the size of scale.");
    endif
  endif

  ## check cov
  if ( !isempty(cov))
    if ( !strcmp(distance,"mahalanobis") || !ismatrix(cov) || !isnumeric(cov))
      error ("knnpredict: Invalid value in cov, cov can only be given for mahalanobis distance.");
    endif
  endif

  ## check bucketsize
  if ( !isscalar(bucketsize) || bucketsize < 0)
    error ("knnpredict: Invalid value of bucketsize.");
  endif

  ##standerdize
  if ( !islogical(standerdize) && standerdize != 1 && standerdize != 0)
    error ("knnpredict: value of standerdize must be boolean.");
  endif


  ##------checking optional params end------##

  if (standerdize)
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
    label     = repmat(classNames(1,:),0,1);
    posterior = NaN(0,classNos);
    cost      = NaN(0,classNos);
  else
    ## calculate the NNs using knnsearch
    [idx, dist] = knnsearch (X, Xclass, "k", k, "NSMethod", NSmethod, ...
                  "Distance", distance, "P", exponent, "Scale", Scale, "cov", ...
                  cov, "bucketsize", bucketsize, "sortindices", true);

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
%!
%! ## calculate 10 nearest-neighbours by minkowski distance
%! [idm, dm] = knnsearch (x, point, "K", 10, "distance", "minkowski", "p", 5);
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idc, dc] = knnsearch (x, point, "K", 10, "distance", "chebychev");
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [ids, ds] = knnsearch (x, point, "K", 10, "distance", "hamming");
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idb, db] = knnsearch (x, point, "K", 10, "distance", "cityblock");
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idh, dh] = knnsearch (x, point, "K", 10, "distance", "manhattan");
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idn, dn] = knnsearch (x, point, "K", 10, "distance", "cosine");
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idj, dj] = knnsearch (x, point, "K", 10, "distance", "jaccard");
%!
%! hold on
%! s1 = subplot (2, 2, 1)
%! gscatter(x(:,1), x(:,2), species,[.75 .75 0; 0 .75 .75; .75 0 .75], '.',20)
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
%! s2 = subplot (2, 2, 2)
%! gscatter(x(:,1), x(:,2), species,[.75 .75 0; 0 .75 .75; .75 0 .75], '.',20)
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
%! s3 = subplot (2, 2, 3)
%! gscatter(x(:,1), x(:,2), species,[.75 .75 0; 0 .75 .75; .75 0 .75], '.',20)
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
%! s4 = subplot (2, 2, 4)
%! gscatter(x(:,1), x(:,2), species,[.75 .75 0; 0 .75 .75; .75 0 .75], '.',20)
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
%! [l, s, c] = knnpredict (x, y, xnew, "k", 10, "distance", "minkowski", "exponent", 5);
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
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "exponent", -3)
%!error<knnpredict: Invalid value in cost or the size of cost.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "cost", [1,2,3;4,5,6;7,8,9])
%!error<knnpredict: Invalid value in cost or the size of cost.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "cost", ones(4,2))
%!error<knnpredict: Invalid value in Scale or the size of scale.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "scale", ones(4,2))
%!error<knnpredict: Invalid value in Scale or the size of scale.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "scale", [1;2;3;-4])
%!error<knnpredict: Invalid value in cov, cov can only be given for mahalanobis distance.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "cov", ones(4,2),"distance","euclidean")
%!error<knnpredict: Invalid value of bucketsize.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "bucketsize",-50)
%!error<knnpredict: value of standerdize must be boolean.> ...
%! knnpredict (ones(4,3),ones(4,1),ones(1,3), "standerdize","some")
