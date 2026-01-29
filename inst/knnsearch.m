## Copyright (C) 2023-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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
## @deftypefn  {statistics} {@var{idx} =} knnsearch (@var{X}, @var{Y})
## @deftypefnx {statistics} {[@var{idx}, @var{D}] =} knnsearch (@var{X}, @var{Y})
## @deftypefnx {statistics} {[@dots{}] =} knnsearch (@dots{}, @var{name}, @var{value})
##
## Find k-nearest neighbors from input data.
##
## @code{@var{idx} = knnsearch (@var{X}, @var{Y})} finds @math{K} nearest
## neighbors in @var{X} for @var{Y}. It returns @var{idx} which contains indices
## of @math{K} nearest neighbors of each row of @var{Y}, If not specified,
## @qcode{@var{K} = 1}.  @var{X} must be an @math{NxP} numeric matrix of input
## data, where rows correspond to observations and columns correspond to
## features or variables.  @var{Y} is an @math{MxP} numeric matrix with query
## points, which must have the same numbers of column as @var{X}.
##
## @code{[@var{idx}, @var{D}] = knnsearch (@var{X}, @var{Y})} also returns the
## the distances, @var{D}, which correspond to the @math{K} nearest neighbour in
## @var{X} for each @var{Y}
##
## Additional parameters can be specified by @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{"K"} @tab @tab is the number of nearest neighbors to be found
## in the kNN search.  It must be a positive integer value and by default it is
## 1.
##
## @item @qcode{"P"} @tab @tab is the Minkowski distance exponent and it must be
## a positive scalar.  This argument is only valid when the selected distance
## metric is @qcode{"minkowski"}.  By default it is 2.
##
## @item @qcode{"Scale"} @tab @tab is the scale parameter for the standardized
## Euclidean distance and it must be a nonnegative numeric vector of equal
## length to the number of columns in @var{X}.  This argument is only valid when
## the selected distance metric is @qcode{"seuclidean"}, in which case each
## coordinate of @var{X} is scaled by the corresponding element of
## @qcode{"scale"}, as is each query point in @var{Y}.  By default, the scale
## parameter is the standard deviation of each coordinate in @var{X}.
##
## @item @qcode{"Cov"} @tab @tab is the covariance matrix for computing the
## mahalanobis distance and it must be a positive definite matrix matching the
## the number of columns in @var{X}.  This argument is only valid when the
## selected distance metric is @qcode{"mahalanobis"}.
##
## @item @qcode{"BucketSize"} @tab @tab is the maximum number of data points in
## the leaf node of the Kd-tree and it must be a positive integer.  This
## argument is only valid when the selected search method is @qcode{"kdtree"}.
##
## @item @qcode{"SortIndices"} @tab @tab is a boolean flag to sort the returned
## indices in ascending order by distance and it is @qcode{true} by default.
## When the selected search method is @qcode{"exhaustive"} or the
## @qcode{"IncludeTies"} flag is true, @code{knnsearch} always sorts the
## returned indices.
##
## @item @qcode{"Distance"} @tab @tab is the distance metric used by
## @code{knnsearch} as specified below:
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"euclidean"} @tab Euclidean distance.
## @item @tab @qcode{"seuclidean"} @tab standardized Euclidean distance.  Each
## coordinate difference between the rows in @var{X} and the query matrix
## @var{Y} is scaled by dividing by the corresponding element of the standard
## deviation computed from @var{X}.  To specify a different scaling, use the
## @qcode{"Scale"} name-value argument.
## @item @tab @qcode{"cityblock"} @tab City block distance.
## @item @tab @qcode{"chebychev"} @tab Chebychev distance (maximum coordinate
## difference).
## @item @tab @qcode{"minkowski"} @tab Minkowski distance.  The default exponent
## is 2.  To specify a different exponent, use the @qcode{"P"} name-value
## argument.
## @item @tab @qcode{"mahalanobis"} @tab Mahalanobis distance, computed using a
## positive definite covariance matrix.  To change the value of the covariance
## matrix, use the @qcode{"Cov"} name-value argument.
## @item @tab @qcode{"cosine"} @tab Cosine distance.
## @item @tab @qcode{"correlation"} @tab One minus the sample linear correlation
## between observations (treated as sequences of values).
## @item @tab @qcode{"spearman"} @tab One minus the sample Spearman's rank
## correlation between observations (treated as sequences of values).
## @item @tab @qcode{"hamming"} @tab Hamming distance, which is the percentage
## of coordinates that differ.
## @item @tab @qcode{"jaccard"} @tab One minus the Jaccard coefficient, which is
## the percentage of nonzero coordinates that differ.
## @item @tab @var{@@distfun} @tab Custom distance function handle.  A distance
## function of the form @code{function @var{D2} = distfun (@var{XI}, @var{YI})},
## where @var{XI} is a @math{1xP} vector containing a single observation in
## @math{P}-dimensional space, @var{YI} is an @math{NxP} matrix containing an
## arbitrary number of observations in the same @math{P}-dimensional space, and
## @var{D2} is an @math{NxP} vector of distances, where @qcode{(@var{D2}k)} is
## the distance between observations @var{XI} and @qcode{(@var{YI}k,:)}.
## @end multitable
##
## @multitable @columnfractions 0.18 0.02 0.8
## @item @qcode{"NSMethod"} @tab @tab is the nearest neighbor search method used
## by @code{knnsearch} as specified below.
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"kdtree"} @tab Creates and uses a Kd-tree to find nearest
## neighbors.  @qcode{"kdtree"} is the default value when the number of columns
## in @var{X} is less than or equal to 10, @var{X} is not sparse, and the
## distance metric is @qcode{"euclidean"}, @qcode{"cityblock"},
## @qcode{"manhattan"}, @qcode{"chebychev"}, or @qcode{"minkowski"}.  Otherwise,
## the default value is @qcode{"exhaustive"}.  This argument is only valid when
## the distance metric is one of the four aforementioned metrics.
## @item @tab @qcode{"exhaustive"} @tab Uses the exhaustive search algorithm by
## computing the distance values from all the points in @var{X} to each point in
## @var{Y}.
## @end multitable
##
## @multitable @columnfractions 0.18 0.02 0.8
## @item @qcode{"IncludeTies"} @tab @tab is a boolean flag to indicate if the
## returned values should contain the indices that have same distance as the
## @math{K^th} neighbor.  When @qcode{false}, @code{knnsearch} chooses the
## observation with the smallest index among the observations that have the same
## distance from a query point.  When @qcode{true}, @code{knnsearch} includes
## all nearest neighbors whose distances are equal to the @math{K^th} smallest
## distance in the output arguments.  To specify @math{K}, use the @qcode{"K"}
## name-value pair argument.
## @end multitable
##
## @seealso{rangesearch, pdist2, fitcknn}
## @end deftypefn

function [idx, dist] = knnsearch (X, Y, varargin)

  ## Check input data
  if (nargin < 2)
	  error ("knnsearch: too few input arguments.");
  endif

  if (size (X, 2) != size (Y, 2))
	  error ("knnsearch: number of columns in X and Y must match.");
  endif

  ## Add default values
  K = 1;                    # Number of nearest neighbors
  P = 2;                    # Exponent for Minkowski distance
  S = [];                   # Scale for the standardized Euclidean distance
  C = [];                   # Covariance matrix for Mahalanobis distance
  BS = 50;                  # Maximum number of points per leaf node for Kd-tree
  SI = true;                # Sort returned indices according to distance
  Distance = "euclidean";   # Distance metric to be used
  NSMethod = [];            # Nearest neighbor search method
  InclTies = false;         # Include ties for distance with kth neighbor
  DistParameter = [];       # Distance parameter for pdist2

  ## Parse additional parameters in Name/Value pairs
  PSC = 0;
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))
      case "k"
        K = varargin{2};
      case "p"
        P = varargin{2};
        PSC += 1;
      case "scale"
        S = varargin{2};
        PSC += 1;
      case "cov"
        C = varargin{2};
        PSC += 1;
      case "bucketsize"
        BS = varargin{2};
      case "sortindices"
        SI = varargin{2};
      case "distance"
        Distance = varargin{2};
      case "nsmethod"
        NSMethod = varargin{2};
      case "includeties"
        InclTies = varargin{2};
      otherwise
        error ("knnsearch: invalid NAME in optional pairs of arguments.");
    endswitch
    varargin(1:2) = [];
  endwhile

  ## Check input parameters
  if (PSC > 1)
    error ("knnsearch: only a single distance parameter can be defined.");
  endif
  if (! isscalar (K) || ! isnumeric (K) || K < 1 || K != round (K))
    error ("knnsearch: invalid value of K.");
  endif
  if (! isscalar (P) || ! isnumeric (P) || P <= 0)
    error ("knnsearch: invalid value of Minkowski Exponent.");
  endif
  if (! isempty (S))
    if (any (S) < 0 || numel (S) != columns (X)
                    || ! strcmpi (Distance, "seuclidean"))
      error ("knnsearch: invalid value in Scale or the size of Scale.");
    endif
  endif
  if (! isempty (C))
    if (! strcmp (Distance, "mahalanobis") || ! ismatrix (C) || ! isnumeric (C))
      error (strcat ("knnsearch: invalid value in Cov, Cov can only", ...
                     " be given for mahalanobis distance."));
    endif
  endif
  if (! isscalar (BS) || BS < 0)
    error ("knnsearch: invalid value of bucketsize.");
  endif

  ## Select the appropriate distance parameter
  if (strcmpi (Distance, "minkowski"))
    DistParameter = P;
  elseif (strcmpi (Distance, "seuclidean"))
    DistParameter = S;
  elseif (strcmpi (Distance, "mahalanobis"))
    DistParameter = C;
  endif

  ## Check NSMethod and set kdtree as default if the conditions match
  if (isempty (NSMethod))
    ## Set default method 'kdtree' if conditions are satisfied;
    if (! issparse (X) && (columns (X) <= 10) &&
       (strcmpi (Distance, "euclidean") || strcmpi (Distance, "cityblock")
     || strcmpi (Distance, "manhattan") || strcmpi (Distance, "minkowski")
     || strcmpi (Distance, "chebychev")))
      NSMethod = "kdtree";
    else
      NSMethod = "exhaustive";
    endif
  else
    ## Check if kdtree can be used
    if (strcmpi (NSMethod, "kdtree") && ! (strcmpi (Distance, "euclidean")
     || strcmpi (Distance, "cityblock") || strcmpi (Distance, "minkowski")
     || strcmpi (Distance, "chebychev")))
      error (strcat ("knnsearch: 'kdtree' cannot be used with", ...
                     " the given distance metric."));
    endif
  endif

  ## Check for NSMethod
  if (strcmpi (NSMethod, "kdtree"))
    ## Build kdtree and search the query point
    kdtree = __build_kdtree__ (1:size(X,1), 0, X, BS);

    ## Check for ties and sortindices
    if (! InclTies)
      ## Only return k neighbors
      dist = zeros (rows (Y), K);
      idx = zeros (rows (Y), K);
      for i = 1:rows (Y)
        [temp_idx, temp_D] = __search_kdtree__ (kdtree, Y(i,:), K, X, ...
																								Distance, DistParameter, ...
																								false);
        if (SI)
          [sorted_D, sort_idx] = sort (temp_D);
          idx(i,:) = temp_idx(sort_idx);
          dist(i,:) = sorted_D;
        else
          idx(i,:) = temp_idx;
          dist(i,:) = temp_D;
        endif
      endfor
    else
      ## Return all neighbors as cell
      dist = cell (rows (Y), 1);
      idx = cell (rows (Y), 1);
      for i = 1:rows (Y)
        [temp_idx, temp_D] = __search_kdtree__ (kdtree, Y(i,:), K, ...
																								X, Distance, DistParameter, ...
																								false);
        r = temp_D(end) + 1e-10; # Add small epsilon to capture ties
        [idx{i}, dist{i}] = __search_kdtree__ (kdtree, Y(i,:), Inf, X, ...
																							 Distance, DistParameter, ...
																							 true, r);
        if (SI)
          [sorted_D, sort_idx] = sort (dist{i});
          dist{i} = sorted_D;
          idx{i} = idx{i}(sort_idx);
        endif
      endfor
    endif
  else

    ## Calculate all distances
    if (K == 1)
      D = pdist2 (X, Y, Distance, DistParameter);
      D = reshape (D', size (Y, 1), size (X, 1));
      [dist, idx] = min (D, [], 2);
    else            # always sort indices in this case
      if (InclTies)
        dist = cell (rows (Y), 1);
        idx = cell (rows (Y), 1);
        for i = 1:rows (Y)
          D = pdist2 (X, Y(i,:), Distance, DistParameter);
          [dt, id] = sort (D);
          kth_dist = dt (K);
          tied_idx = (dt <= kth_dist);
          dist {i} = dt(tied_idx, :)';
          idx {i} = id(tied_idx, :)';
        endfor
      else
        ## No ties included
        D = pdist2 (X, Y, Distance, DistParameter);
        D = reshape (D', size (Y, 1), size (X, 1));
        [dist, idx] = sort (D, 2);
        dist = dist(:,1:K);
        idx = idx(:,1:K);
      endif
    endif
  endif

endfunction

## buildkdtree
function node = __build_kdtree__ (indices, depth, X, bucket_size)
  if (length (indices) <= bucket_size)
    node = struct ('indices', indices);
  else
    k = size (X, 2);
    axis = mod (depth, k) + 1;
    values = X(indices, axis);
    [sorted_values, sort_idx] = sort (values);
    sorted_indices = indices(sort_idx);
    median_idx = floor ((length (indices) + 1) / 2);
    split_value = sorted_values(median_idx);
    left_indices = indices(values <= split_value);
    right_indices = indices(values > split_value);
    if (isempty (left_indices) || isempty (right_indices))
      node = struct ('indices', indices);
      return;
    endif
    left_node = __build_kdtree__ (left_indices, depth + 1, X, bucket_size);
    right_node = __build_kdtree__ (right_indices, depth + 1, X, bucket_size);
    node = struct ('axis', axis, 'split_value', split_value, ...
                   'left', left_node, 'right', right_node);
  endif
endfunction

## Search KD-tree
function [indices, distances] = __search_kdtree__ (node, query, k, X, dist, ...
																									 distparam, is_range, r)
  if (nargin < 8)
    r = Inf;
  endif
  if (strcmpi (dist, "minkowski"))
    if (! (isscalar (distparam) && isnumeric (distparam) ...
                                && distparam > 0 && isfinite (distparam)))
      error (strcat("knnsearch.__search_kdtree__:", ...
                    " distparam must be a positive finite", ...
                    " scalar for minkowski."));
    endif
  else
    if (! isempty (distparam))
      error (strcat("knnsearch.__search_kdtree__:", ...
                    " distparam must be empty for non-minkowski metrics."));
    endif
  endif
  indices = [];
  distances = [];
  search (node, 0);

  function search (node, depth)
    if (isempty (node))
      return;
    endif

    if (isfield (node, 'indices'))
      leaf_indices = node.indices;
      if (strcmpi (dist, "minkowski"))
        dists = pdist2 (X(leaf_indices,:), query, dist, distparam);
      else
        dists = pdist2 (X(leaf_indices,:), query, dist);
      endif
      if (is_range)
        mask = dists <= r;
        indices = [indices; leaf_indices(mask)'];
        distances = [distances; dists(mask)];
      else
        indices = [indices; leaf_indices'];
        distances = [distances; dists];
        if (length (distances) > k)
          [distances, sort_idx] = sort (distances);
          indices = indices(sort_idx);
          distances = distances(1:k);
          indices = indices(1:k);
        endif
      endif
    else
      axis = node.axis;
      split_value = node.split_value;
      if (query(axis) <= split_value)
        nearer = node.left;
        further = node.right;
      else
        nearer = node.right;
        further = node.left;
      endif

      search (nearer, depth + 1);

      plane_dist = abs (query(axis) - split_value);
      if (is_range)
        max_dist = r;
        if (plane_dist <= max_dist)
          search (further, depth + 1);
        endif
      else
        if (length (distances) < k || plane_dist < distances(end))
          search (further, depth + 1);
        endif
      endif
    endif
  endfunction
endfunction

%!demo
%! ## find 10 nearest neighbour of a point using different distance metrics
%! ## and compare the results by plotting
%! load fisheriris
%! X = meas(:,3:4);
%! Y = species;
%! point = [5, 1.45];
%!
%! ## calculate 10 nearest-neighbours by minkowski distance
%! [id, d] = knnsearch (X, point, "K", 10);
%!
%! ## calculate 10 nearest-neighbours by minkowski distance
%! [idm, dm] = knnsearch (X, point, "K", 10, "distance", "minkowski", "p", 5);
%!
%! ## calculate 10 nearest-neighbours by chebychev distance
%! [idc, dc] = knnsearch (X, point, "K", 10, "distance", "chebychev");
%!
%! ## plotting the results
%! gscatter (X(:,1), X(:,2), species, [.75 .75 0; 0 .75 .75; .75 0 .75], ".", 20);
%! title ("Fisher's Iris Data - Nearest Neighbors with different types of distance metrics");
%! xlabel("Petal length (cm)");
%! ylabel("Petal width (cm)");
%!
%! line (point(1), point(2), "marker", "X", "color", "k", ...
%!       "linewidth", 2, "displayname", "query point")
%! line (X(id,1), X(id,2), "color", [0.5 0.5 0.5], "marker", "o", ...
%!       "linestyle", "none", "markersize", 10, "displayname", "euclidean")
%! line (X(idm,1), X(idm,2), "color", [0.5 0.5 0.5], "marker", "d", ...
%!       "linestyle", "none", "markersize", 10, "displayname", "Minkowski")
%! line (X(idc,1), X(idc,2), "color", [0.5 0.5 0.5], "marker", "p", ...
%!       "linestyle", "none", "markersize", 10, "displayname", "chebychev")
%! xlim ([4.5 5.5]);
%! ylim ([1 2]);
%! axis square;

%!demo
%! ## knnsearch on iris dataset using kdtree method
%! load fisheriris
%! X = meas(:,3:4);
%! gscatter (X(:,1), X(:,2), species, [.75 .75 0; 0 .75 .75; .75 0 .75], ".", 20);
%! title ("Fisher's iris dataset : Nearest Neighbors with kdtree search");
%!
%! ## new point to be predicted
%! point = [5 1.45];
%!
%! line (point(1), point(2), "marker", "X", "color", "k", ...
%!       "linewidth", 2, "displayname", "query point")
%!
%! ## knnsearch using kdtree method
%! [idx, d] = knnsearch (X, point, "K", 10, "NSMethod", "kdtree");
%!
%! ## plotting predicted neighbours
%! line (X(idx,1), X(idx,2), "color", [0.5 0.5 0.5], "marker", "o", ...
%!       "linestyle", "none", "markersize", 10, ...
%!       "displayname", "nearest neighbour")
%! xlim ([4 6])
%! ylim ([1 3])
%! axis square
%! ## details of predicted labels
%! tabulate (species(idx))
%!
%! ctr = point - d(end);
%! diameter = 2 * d(end);
%! ##  Draw a circle around the 10 nearest neighbors.
%! h = rectangle ("position", [ctr, diameter, diameter], "curvature", [1 1]);
%!
%! ## here only 8 neighbours are plotted instead of 10 since the dataset
%! ## contains duplicate values


## Test output
%!shared X, Y
%! X = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
%! Y = [1, 2, 2, 3; 2, 3, 3, 4];
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "euclidean");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! eucldist = @(v,m) sqrt(sumsq(repmat(v,rows(m),1)-m,2));
%! [idx, D] = knnsearch (X, Y, "Distance", eucldist);
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "euclidean", "includeties", true);
%! assert (iscell (idx), true);
%! assert (iscell (D), true)
%! assert (idx {1}, [1]);
%! assert (idx {2}', [1, 2]);
%! assert (D{1}, ones (1, 1) * sqrt (2));
%! assert (D{2}', ones (1, 2) * sqrt (2));
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "euclidean", "k", 2);
%! assert (idx, [1, 2; 1, 2]);
%! assert (D, [sqrt(2), 3.162277660168380; sqrt(2), sqrt(2)], 1e-14);
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "seuclidean");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "seuclidean", "k", 2);
%! assert (idx, [1, 2; 1, 2]);
%! assert (D, [sqrt(2), 3.162277660168380; sqrt(2), sqrt(2)], 1e-14);
%!test
%! xx = [1, 2; 1, 3; 2, 4; 3, 6];
%! yy = [2, 4; 2, 6];
%! [idx, D] = knnsearch (xx, yy, "Distance", "mahalanobis");
%! assert (idx, [3; 2]);
%! assert (D, [0; 3.162277660168377], 1e-14);
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "minkowski");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "minkowski", "p", 3);
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * 1.259921049894873, 1e-14);
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "cityblock");
%! assert (idx, [1; 1]);
%! assert (D, [2; 2]);
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "chebychev");
%! assert (idx, [1; 1]);
%! assert (D, [1; 1]);
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "cosine");
%! assert (idx, [2; 3]);
%! assert (D, [0.005674536395645; 0.002911214328620], 1e-14);
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "correlation");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * 0.051316701949486, 1e-14);
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "spearman");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * 0.051316701949486, 1e-14);
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "hamming");
%! assert (idx, [1; 1]);
%! assert (D, [0.5; 0.5]);
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "jaccard");
%! assert (idx, [1; 1]);
%! assert (D, [0.5; 0.5]);
%!test
%! [idx, D] = knnsearch (X, Y, "Distance", "jaccard", "k", 2);
%! assert (idx, [1, 2; 1, 2]);
%! assert (D, [0.5, 1; 0.5, 0.5]);
%!test
%! a = [1, 5; 1, 2; 2, 2; 1.5, 1.5; 5, 1; 2 -1.34; 1, -3; 4, -4; -3, 1; 8, 9];
%! b = [1, 1];
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree", "includeties", true);
%! assert (iscell (idx), true);
%! assert (iscell (D), true)
%! assert (cell2mat (idx)', [4, 2, 3, 6, 1, 5, 7, 9]);
%! assert (cell2mat (D)', [0.7071, 1.0000, 1.4142, 2.5447, 4.0000, 4.0000, 4.0000, 4.0000], 1e-4);
%!test
%! a = [1, 5; 1, 2; 2, 2; 1.5, 1.5; 5, 1; 2 -1.34; 1, -3; 4, -4; -3, 1; 8, 9];
%! b = [1, 1];
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "exhaustive", "includeties", true);
%! assert (iscell (idx), true);
%! assert (iscell (D), true)
%! assert (cell2mat (idx), [4, 2, 3, 6, 1, 5, 7, 9]);
%! assert (cell2mat (D), [0.7071, 1.0000, 1.4142, 2.5447, 4.0000, 4.0000, 4.0000, 4.0000], 1e-4);
%!test
%! a = [1, 5; 1, 2; 2, 2; 1.5, 1.5; 5, 1; 2 -1.34; 1, -3; 4, -4; -3, 1; 8, 9];
%! b = [1, 1];
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree", "includeties", false);
%! assert (iscell (idx), false);
%! assert (iscell (D), false)
%! assert (idx, [4, 2, 3, 6, 1]);
%! assert (D, [0.7071, 1.0000, 1.4142, 2.5447, 4.0000], 1e-4);
%!test
%! a = [1, 5; 1, 2; 2, 2; 1.5, 1.5; 5, 1; 2 -1.34; 1, -3; 4, -4; -3, 1; 8, 9];
%! b = [1, 1];
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "exhaustive", "includeties", false);
%! assert (iscell (idx), false);
%! assert (iscell (D), false)
%! assert (idx, [4, 2, 3, 6, 1]);
%! assert (D, [0.7071, 1.0000, 1.4142, 2.5447, 4.0000], 1e-4);
%!test
%! load fisheriris
%! a = meas;
%! b = min(meas);
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree");
%! assert (idx, [42, 9, 14, 39, 13]);
%! assert (D, [0.5099, 0.9950, 1.0050, 1.0536, 1.1874], 1e-4);
%!test
%! load fisheriris
%! a = meas;
%! b = mean(meas);
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree");
%! assert (idx, [65, 83, 89, 72, 100]);
%! assert (D, [0.3451, 0.3869, 0.4354, 0.4481, 0.4625], 1e-4);
%!test
%! load fisheriris
%! a = meas;
%! b = max(meas);
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree");
%! assert (idx, [118, 132, 110, 106, 136]);
%! assert (D, [0.7280, 0.9274, 1.3304, 1.5166, 1.6371], 1e-4);
%!
%!test
%! load fisheriris
%! a = meas;
%! b = max(meas);
%! [idx, D] = knnsearch (a, b, "K", 5, "includeties", true);
%! assert (iscell (idx), true);
%! assert (iscell (D), true);
%! assert (cell2mat (idx)', [118, 132, 110, 106, 136]);
%! assert (cell2mat (D)', [0.7280, 0.9274, 1.3304, 1.5166, 1.6371], 1e-4);

## Test input validation
%!error<knnsearch: too few input arguments.> knnsearch (1)
%!error<knnsearch: number of columns in X and Y must match.> ...
%! knnsearch (ones (4, 5), ones (4))
%!error<knnsearch: invalid NAME in optional pairs of arguments.> ...
%! knnsearch (ones (4, 2), ones (3, 2), "Distance", "euclidean", "some", "some")
%!error<knnsearch: only a single distance parameter can be defined.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "scale", ones (1, 5), "P", 3)
%!error<knnsearch: invalid value of K.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "K", 0)
%!error<knnsearch: invalid value of Minkowski Exponent.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "P", -2)
%!error<knnsearch: invalid value in Scale or the size of Scale.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "scale", ones(4,5), "distance", "euclidean")
%!error<knnsearch: invalid value in Cov, Cov can only be given for mahalanobis distance.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "cov", ["some" "some"])
%!error<knnsearch: invalid value in Cov, Cov can only be given for mahalanobis distance.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "cov", ones(4,5), "distance", "euclidean")
%!error<knnsearch: invalid value of bucketsize.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "bucketsize", -1)
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "NSmethod", "kdtree", "distance", "cosine")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "NSmethod", "kdtree", "distance", "mahalanobis")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "NSmethod", "kdtree", "distance", "correlation")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "NSmethod", "kdtree", "distance", "seuclidean")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "NSmethod", "kdtree", "distance", "spearman")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "NSmethod", "kdtree", "distance", "hamming")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch (ones (4, 5), ones (1, 5), "NSmethod", "kdtree", "distance", "jaccard")
%!test
%! X = ones (10, 2);
%! Y = X(1,:);
%! [idx, D] = knnsearch (X, Y, "NSMethod", "kdtree", "BucketSize", 1);
%! assert (numel (idx), 1);
%! assert (idx(1), 1);
