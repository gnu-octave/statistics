## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{idx} =} knnsearch (@var{x}, @var{y})
## @deftypefnx {statistics} {[@var{idx}, @var{D}] =} knnsearch (@var{x}, @var{y})
## @deftypefnx {statistics} {[@dots{}] =} knnsearch (@var{x}, @var{y}, @var{name}, @var{value})
##
## Find k-nearest neighbors using input data
##
## @code{@var{idx} = knnsearch (@var{x}, @var{y})} finds @math{K} nearest
## neighbors in @var{x} for @var{y}. It returns @var{idx} which contains indices
## of @math{K} nearest neighbors of each row of @var{y}, If not specified,
## @qcode{@var{K} = 1}.  @var{x} must be an @math{NxP} numeric matrix of input
## data, where rows correspond to observations and columns correspond to
## features or variables.  @var{y} is an @math{MxP} numeric matrix with query
## points, which must have the same numbers of column as @var{x}.
##
## @code{[@var{idx}, @var{dist}] = knnsearch (@var{x}, @var{y})} returns the
## @math{K} nearest neighbour in @var{x} for each @var{y} with distances
## returned in @var{dist}.
##
## Additional input arguments can be given as name-value pairs.
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"K"} @tab is the number of nearest neighbors to be found
## in the kNN search.  It must be a positive integer value and by default it is
## 1.
##
## @item @tab @qcode{"P"} @tab is the Minkowski distance exponent and it must be
## a positive scalar.  This argument is only valid when the selected distance
## metric is @qcode{"minkowski"}.  By default it is 2.
##
## @item @tab @qcode{"scale"} @tab is the scale parameter for the standardized
## Euclidean distance and it must be a nonnegative numeric vector of equal
## length to the number of columns in @var{X}.  This argument is only valid when
## the selected distance metric is @qcode{"seuclidean"}, in which case each
## coordinate of @var{X} is scaled by the corresponding element of
## @qcode{"scale"}, as is each query point in @var{Y}.  By default, the scale
## parameter is the standard deviation of each coordinate in @var{X}.
##
## @item @tab @qcode{"cov"} @tab is the covariance matrix for computing the
## mahalanobis distance and it must be a positive definite matrix matching the
## the number of columns in @var{X}.  This argument is only valid when the
## selected distance metric is @qcode{"mahalanobis"}.
##
## @item @tab @qcode{"BucketSize"} @tab is the maximum number of data points in
## the leaf node of the Kd-tree and it must be a positive integer.  This
## argument is only valid when the selected search method is @qcode{"kdtree"}.
##
## @item @tab @qcode{"SortIndices"} @tab is a boolean flag to sort the returned
## indices in ascending order by distance and it is @qcode{true} by default.
## When the selected search method is @qcode{"exhaustive"} or the
## @qcode{"IncludeTies"} flag is true, @code{knnsearch} always sorts the
## returned indices.
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
## @item @tab @qcode{"IncludeTies"} @tab is a boolean flag to indicate if the
## returned values should contain the indices that have same distance as the
## @math{K^th} neighbor.  When @qcode{false}, @code{knnsearch} chooses the
## observation with the smallest index among the observations that have the same
## distance from a query point.  When @qcode{true}, @code{knnsearch} includes
## all nearest neighbors whose distances are equal to the @math{K^th} smallest
## distance in the output arguments.  To specify @math{K}, use the @qcode{"K"}
## name-value pair argument.
## @end multitable
##
## @seealso{rangesearch}
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

  ## Parse additional parameters in Name/Value pairs
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))
      case "k"
        K = varargin{2};
      case "p"
        P = varargin{2};
      case "scale"
        S = varargin{2};
      case "cov"
        C = varargin{2};
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
  if (! isscalar (K) || ! isnumeric (K) || K < 1 || K != round (K))
    error ("knnsearch: Invalid value of K.");
  endif
  if (! isscalar (P) || ! isnumeric (P) || P < 0)
    error ("knnsearch: Invalid value of Minkowski Exponent.");
  endif
  if (! isempty (C))
    if (! strcmp (Distance, "mahalanobis") || ! ismatrix (C) || ! isnumeric (C))
      error (strcat (["knnsearch: Invalid value in cov, cov can only"], ...
                     [" be given for mahalanobis distance."]));
    endif
  endif
  if (! isempty (S))
    if (! isscalar (S) || any (S) < 0 || numel (S) != rows (X) ...
                       || ! strcmpi (Distance, "seuclidean"))
      error ("knnsearch: Invalid value in Scale or the size of scale.");
    endif
  endif
  if (! isscalar (BS) || BS < 0)
    error ("knnsearch: Invalid value of bucketsize.");
  endif

  ## check NSMethod and set kdtree as default if the conditions match
  if ( isempty (NSMethod))
    ## set default method 'kdtree' if condintions are satistfied;
    if (! issparse (X) && (columns (X) <= 10) && ...
       ( strcmpi (Distance, "euclidean") || strcmpi (Distance, "cityblock")...
      || strcmpi (Distance, "chebychev") || strcmpi (Distance, "minkowski")))
      NSMethod = "kdtree";
    else
      NSMethod = "exhaustive";
    endif
  else
    ## not empty then check if is exhaustive or kdtree
    if ( strcmpi (NSMethod,"kdtree") && ! ( strcmpi (Distance, "euclidean") ...
      || strcmpi (Distance, "cityblock") || strcmpi (Distance, "chebychev") ...
      || strcmpi (Distance, "minkowski")))
      error ("knnsearch: 'kdtree' cannot be used with the given distance metric.");
    endif
  endif

  ## Check for NSMethod
  if (strcmpi (NSMethod, "kdtree"))
    ## Build kdtree and search the query point
    ret  = buildkdtree (X, BS);

    ## Check for ties and sortindices
    if (! InclTies)
      ## only return k neighbors
      ##  no need for returning cell
      dist = [];
      idx  = [];
      for i = 1:rows (Y)
        NN = findkdtree (ret, Y (i, :), K, Distance, P, S, C);
        D = calc_dist (X(NN,:), Y (i,:), Distance, P, S, C);
        sorted_D = sortrows ([NN, D], 2);
        dist = [dist; sorted_D(1:K, 2)'];
        idx  = [idx;  sorted_D(1:K, 1)'];
      endfor
      if (SI)
        ## rows are already sorted by distance
        dist = (dist);
        idx  = (idx);
      else
        dist = (dist);
        idx  = (idx);
      endif
    else
      ## return all neighbors as cell
      dist = cell (rows (Y), 1);
      idx  = cell (rows (Y), 1);
      for i = 1:rows (Y)
        NN = findkdtree (ret, Y (i, :), K, Distance, P, S, C);
        D = calc_dist (X(NN,:), Y (i,:), Distance, P, S, C);
        sorted_D = sortrows ([NN, D], 2);
        kth_dist = sorted_D (K, 2);
        tied_idx = (sorted_D (:, 2) <= kth_dist);
        dist {i} = sorted_D (tied_idx, 2)';
        idx  {i} = sorted_D (tied_idx, 1)';
      endfor
    endif
  else

    ## Calculate distance and search by exhaustive
    if (K == 1)
      D = calc_dist (X, Y, Distance, P, S, C);
      D = reshape (D, size (Y, 1), size (X, 1));
      [dist, idx] = min (D, [], 2);
    else            # always sort indices in this case
      if (InclTies)
        ## this part needs fixing so that idx is a cell array (column vector
        dist = cell (rows (Y), 1);
        idx  = cell (rows (Y), 1);
        for i = 1:rows (Y)
          D = calc_dist (X, Y(i,:), Distance, P, S, C);
          [dt, id] = sort (D);
          kth_dist = dt (K);
          tied_idx = (dt <= kth_dist);
          dist {i} = dt (tied_idx, :);
          idx  {i} = id (tied_idx, :);
        endfor
      else
        ## no ties included
        D = calc_dist (X, Y, Distance, P, S, C);
        D = reshape (D, size (Y, 1), size (X, 1));
        [dist, idx] = sort (D, 2);
        dist = dist(:,1:K);
        idx  = idx(:,1:K);
      endif
    endif
  endif

endfunction

function D = calc_dist (X, Y, Distance, P, S, C)
  [ix, iy] = meshgrid (1:size (X, 1), 1:size (Y, 1));
  if strcmpi (Distance, "euclidean")
    D = sqrt (sum ((X(ix(:),:) - Y(iy(:),:)) .^ 2, 2));

  elseif strcmpi (Distance, "seuclidean")
    if (isempty (S))
      S = std (X, [], 1);
    endif
    mu = mean (X, 1);
    sx = (X - mu) ./ S;
    sy = (Y - mu) ./ S;
    D = sqrt (sum ((sx(ix(:),:) - sy(iy(:),:)) .^ 2, 2));

  elseif (strcmpi (Distance, "mahalanobis"))
    if isempty(C)
      C = cov (X(! any (isnan (X), 2),:));
    endif
    dxy = X(ix(:),:) - Y(iy(:),:);
    D   = sqrt (sum ((dxy  * inv (C)) .* dxy, 2));

  elseif (strcmpi (Distance, "minkowski"))
    D = sum (abs (X(ix(:),:) - Y(iy(:),:)) .^ P, 2) .^ (1 / P);

  elseif (strcmpi (Distance, "cityblock") || strcmpi (Distance, "manhattan"))
    D = sum (abs (X(ix(:),:) - Y(iy(:),:)), 2);

  elseif (strcmpi (Distance, "chebychev"))
    D = max (abs (X(ix(:),:) - Y(iy(:),:)), [], 2);

  elseif (strcmpi (Distance, "cosine"))
    sx = sum (X .^ 2, 2) .^ (-1 / 2);
    sy = sum (Y .^ 2, 2) .^ (-1 / 2);
    D  = 1 - sum (X(ix(:),:) .* Y(iy(:),:), 2) .* sx(ix(:)) .* sy(iy(:));

  elseif (strcmp (Distance, "correlation"))
    mX = mean (X(ix(:),:), 2);
    mY = mean (Y(iy(:),:), 2);
    xy = sum ((X(ix(:),:) - mX) .* (Y(iy(:),:) - mY), 2);
    xx = sqrt (sum ((X(ix(:),:) - mX) .* (X(ix(:),:) - mX), 2));
    yy = sqrt (sum ((Y(iy(:),:) - mY) .* (Y(iy(:),:) - mY), 2));
    D = 1 - (xy ./ (xx .* yy));

  elseif (strcmpi (Distance, "spearman"))
    for i = 1:size (X, 1)
      rX(i,:) = tiedrank (X(i,:));
    endfor
    for i = 1:size (Y, 1)
      rY(i,:) = tiedrank (Y(i,:));
    endfor
    rM = (size (X, 2) + 1) / 2;
    xy = sum ((rX(ix(:),:) - rM) .* (rY(iy(:),:) - rM), 2);
    xx = sqrt (sum ((rX(ix(:),:) - rM) .* (rX(ix(:),:) - rM), 2));
    yy = sqrt (sum ((rY(iy(:),:) - rM) .* (rY(iy(:),:) - rM), 2));
    D = 1 - (xy ./ (xx .* yy));

  elseif (strcmpi (Distance, "hamming"))
    D = mean (abs (X(ix(:),:) != Y(iy(:),:)), 2);

  elseif (strcmpi (Distance, "jaccard"))
    xy0 = (X(ix(:),:) != 0 | Y(iy(:),:) != 0);
    D = sum ((X(ix(:),:) != Y(iy(:),:)) & xy0, 2) ./ sum (xy0, 2);
  endif

endfunction

## buildkdtree
function ret = buildkdtree_recur (x, r, d, BS)
  count = length (r);
  dimen = size (x, 2);
  if (count == 1)
    ret = struct ("point", r(1), "dimen", d);
  else
    mid = ceil (count / 2);
    ret = struct ("point", r(mid), "dimen", d);
    d = mod (d, dimen) + 1;
    ## Build left sub tree
    if (mid > 1)
      left = r(1:mid-1);
      left_points = x(left,d);
      [val, left_idx] = sort (left_points);
      leftr = left(left_idx);
      ret.left = buildkdtree_recur (x, leftr, d);
    endif
    ## Build right sub tree
    if (count > mid)
      right = r(mid+1:count);
      right_points = x(right,d);
      [val, right_idx] = sort (right_points);
      rightr = right(right_idx);
      ret.right = buildkdtree_recur (x, rightr, d);
    endif
  endif
endfunction

## wrapper function for buildkdtree_recur
function ret = buildkdtree (x, BS)
  [val, r] = sort (x(:,1));
  ret = struct ("data", x, "root", buildkdtree_recur (x, r, 1, BS));
endfunction

function farthest = kdtree_cand_farthest (x, p, cand, distance, P, S, C)
  [val, index] = max (calc_dist (x, p, distance, P, S, C)(cand));
  farthest = cand (index);
endfunction

## function to insert into NN list
function inserted = kdtree_cand_insert (x, p, cand, k, point, distance, P, S, C)
  if (length (cand) < k)
    inserted = [cand; point];
  else
    farthest = kdtree_cand_farthest (x, p, cand, distance, P, S, C);
    if (calc_dist (cand(find(cand == farthest),:), point, distance, P, S, C))
      inserted = [cand; point];
    else
      farthest = kdtree_cand_farthest (x, p, cand, distance, P, S, C);
      cand (find (cand == farthest)) = point;
      inserted = cand;
    endif
  endif
endfunction

## function to search in a kd tree
function neighbours = findkdtree_recur (x, node, p, neighbours, ...
                                        k, distance, P, S, C)
  point = node.point;
  d = node.dimen;
  if (x(point,d) > p(d))
    ## Search in left sub tree
    if (isfield (node, "left"))
      neighbours = findkdtree_recur (x, node.left, p, neighbours, ...
                                     k, distance, P, S, C);
    endif
    ## Add current point if neccessary
    farthest = kdtree_cand_farthest (x, p, neighbours, distance, P, S, C);
    if (length(neighbours) < k || calc_dist (x(point,:), p, distance, P, S, C)
        <= calc_dist (x(farthest,:), p, distance, P, S, C))
      neighbours = kdtree_cand_insert (x, p, neighbours, ...
                                       k, point, distance, P, S, C);
    endif
    ## Search in right sub tree if neccessary
    farthest = kdtree_cand_farthest (x, p, neighbours, distance, P, S, C);
    radius = calc_dist (x(farthest,:), p, distance, P, S, C);
    if (isfield (node, "right") &&
        (length(neighbours) < k || p(d) + radius > x(point,d)))
      neighbours = findkdtree_recur (x, node.right, p, neighbours, ...
                                     k, distance, P, S, C);
    endif
  else
    ## Search in right sub tree
    if (isfield (node, "right"))
      neighbours = findkdtree_recur (x, node.right, p, neighbours, ...
                                     k, distance, P, S, C);
    endif
    ## Add current point if neccessary
    farthest = kdtree_cand_farthest (x, p, neighbours, distance, P, S, C);
    if (length (neighbours) < k || calc_dist (x(point,:), p, distance, P, S, C)
        <= calc_dist (x(farthest,:), p, distance, P, S, C))
      neighbours = kdtree_cand_insert (x, p, neighbours, ...
                                       k, point, distance, P, S, C);
    endif
    ## Search in left sub tree if neccessary
    farthest = kdtree_cand_farthest (x, p, neighbours, distance, P, S, C);
    radius = calc_dist (x(farthest,:), p, distance, P, S, C);
    if (isfield (node, "left") &&
        (length (neighbours) < k || p(d) - radius <= x(point,d)))
      neighbours = findkdtree_recur (x, node.left, p, neighbours, ...
                                     k, distance, P, S, C);
    endif
  endif
endfunction

## wrapper function for findkdtree_recur
function neighbours = findkdtree (tree, p, k, distance, P, S, C)
    x = tree.data;
    root = tree.root;
    neighbours = findkdtree_recur (x, root, p, [], k, distance, P, S, C);
endfunction

%!demo
%! ## find 10 nearest neighbour of a point using different distance metrics
%! ## and compare the results by plotting
%! load fisheriris
%! x = meas(:,3:4);
%! y = species;
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
%! ## plotting the results
%! gscatter (x(:,1), x(:,2), species, [.75 .75 0; 0 .75 .75; .75 0 .75], ".", 20)
%! title ("Fisher's Iris Data - Nearest Neighbors with different types of distance metrics");
%! xlabel("Petal length (cm)");
%! ylabel("Petal width (cm)");
%!
%! line (point(1), point(2), "marker", "x", "color", "k", ...
%!       "linewidth", 2, "displayname", "query point")
%! line (x(id,1), x(id,2), "color", [0.5 0.5 0.5], "marker", "o", ...
%!       "linestyle", "none", "markersize", 10, "displayname", "eulcidean")
%! line (x(idm,1), x(idm,2), "color", [0.5 0.5 0.5], "marker", "d", ...
%!       "linestyle", "none", "markersize", 10, "displayname", "Minkowski")
%! line (x(idc,1), x(idc,2), "color", [0.5 0.5 0.5], "marker", "p", ...
%!       "linestyle", "none", "markersize", 10, "displayname", "chebychev")
%! xlim ([4.5 5.5]);
%! ylim ([1 2]);
%! axis square;
%!
%!demo
%! ## knnsearch on iris dataset using kdtree method
%! load fisheriris
%! x = meas(:,3:4);
%! gscatter (x(:,1), x(:,2), species, [.75 .75 0; 0 .75 .75; .75 0 .75], ".", 20)
%! title ("Fisher's iris dataset : Nearest Neighbors with kdtree search");
%!
%! ## new point to be predicted
%! point = [5 1.45];
%!
%! line (point(1), point(2), "marker", "x", "color", "k", ...
%!       "linewidth", 2, "displayname", "query point")
%!
%! ## knnsearch using kdtree method
%! [idx, d] = knnsearch (x, point, "K", 10, "NSMethod", "kdtree");
%!
%! ## plotting predicted neighbours
%! line (x(idx,1), x(idx,2), "color", [0.5 0.5 0.5], "marker", "o", ...
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
%! contains duplicate values


## Test output
%!shared x, y
%! x = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
%! y = [1, 2, 2, 3; 2, 3, 3, 4];
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "euclidean");
%! assert (idx, [1; 2]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "euclidean", "includeties", true);
%! assert ( iscell (idx), true);
%! assert ( iscell (D), true)
%! assert (idx {1}, [1]);
%! assert (idx {2}, [2, 1]);
%! assert (D {1}, ones (1, 1) * sqrt (2));
%! assert (D {2}, ones (1, 2) * sqrt (2));
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "euclidean", "k", 2);
%! assert (idx, [1, 2; 2, 1]);
%! assert (D, [sqrt(2), 3.162277660168380; sqrt(2), sqrt(2)], 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "seuclidean");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "seuclidean", "k", 2);
%! assert (idx, [1, 2; 1, 2]);
%! assert (D, [sqrt(2), 3.162277660168380; sqrt(2), sqrt(2)], 1e-14);
%!test
%! xx = [1, 2; 1, 3; 2, 4; 3, 6];
%! yy = [2, 4; 2, 6];
%! [idx, D] = knnsearch (xx, yy, "Distance", "mahalanobis");
%! assert (idx, [3; 2]);
%! assert (D, [0; 3.162277660168377], 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "minkowski");
%! assert (idx, [1; 2]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "minkowski", "p", 3);
%! assert (idx, [1; 2]);
%! assert (D, ones (2, 1) * 1.259921049894873, 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "cityblock");
%! assert (idx, [1; 2]);
%! assert (D, [2; 2]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "manhattan");
%! assert (idx, [1; 1]);
%! assert (D, [2; 2]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "chebychev");
%! assert (idx, [1; 2]);
%! assert (D, [1; 1]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "cosine");
%! assert (idx, [2; 3]);
%! assert (D, [0.005674536395645; 0.002911214328620], 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "correlation");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * 0.051316701949486, 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "spearman");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * 0.051316701949486, 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "hamming");
%! assert (idx, [1; 1]);
%! assert (D, [0.5; 0.5]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "jaccard");
%! assert (idx, [1; 1]);
%! assert (D, [0.5; 0.5]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "jaccard", "k", 2);
%! assert (idx, [1, 2; 1, 2]);
%! assert (D, [0.5, 1; 0.5, 0.5]);
%!test
%! a = [1, 5; 1, 2; 2, 2; 1.5, 1.5; 5, 1; 2 -1.34; 1, -3; 4, -4; -3, 1; 8, 9];
%! b = [1, 1];
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree","includeties",true);
%! assert ( iscell (idx), true);
%! assert ( iscell (D), true)
%! assert (cell2mat (idx), [4, 2, 3, 6, 1, 9, 7, 5]);
%! assert (cell2mat (D), [0.7071, 1.0000, 1.4142, 2.5447, 4.0000, 4.0000, 4.0000, 4.0000],1e-4);
%!test
%! a = [1, 5;1,	2;2,	2;1.5,	1.5;5,	1;2	-1.34;1,	-3;4,	-4;-3,	1;8,	9];
%! b = [1, 1];
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree","includeties",false);
%! assert (idx, [4, 2, 3, 6, 1]);
%! assert (D, [0.7071, 1.0000, 1.4142, 2.5447, 4.0000],1e-4);
%!test
%! a = [1, 5;1,	2;2,	2;1.5,	1.5;5,	1;2	-1.34;1,	-3;4,	-4;-3,	1;8,	9];
%! b = [1, 1];
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "exhaustive", "includeties", false);
%! assert (idx, [4, 2, 3, 6, 1]);
%! assert (D, [0.7071, 1.0000, 1.4142, 2.5447, 4.0000],1e-4);
%!test
%! load fisheriris
%! a = meas;
%! b = min(meas);
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree");
%! assert (idx, [42, 9, 14, 39, 13]);
%! assert (D, [0.5099, 0.9950, 1.0050, 1.0536, 1.1874],1e-4);
%!test
%! load fisheriris
%! a = meas;
%! b = mean(meas);
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree");
%! assert (idx, [65, 83, 89, 72, 100]);
%! assert (D, [0.3451, 0.3869, 0.4354, 0.4481, 0.4625],1e-4);
%!test
%! load fisheriris
%! a = meas;
%! b = max(meas);
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree");
%! assert (idx, [118, 132, 110, 106, 136]);
%! assert (D, [0.7280, 0.9274, 1.3304, 1.5166, 1.6371],1e-4);
%!
%!test
%! load fisheriris
%! a = meas;
%! b = max(meas);
%! [idx, D] = knnsearch (a, b, "K", 5, "includeties", true);
%! assert ( iscell (idx), true);
%! assert ( iscell (D), true);
%! assert (cell2mat (idx), [118, 132, 110, 106, 136]);
%! assert (cell2mat (D), [0.7280, 0.9274, 1.3304, 1.5166, 1.6371],1e-4);

## Test input validation
%!error<knnsearch: too few input arguments.> knnsearch (1)
%!error<knnsearch: number of columns in X and Y must match.> ...
%! knnsearch (ones (4, 5), ones (4))
%!error<knnsearch: invalid NAME in optional pairs of arguments.> ...
%! knnsearch (ones (4, 2), ones (3, 2), "Distance", "euclidean", "some", "some")
%!error<knnsearch: Invalid value of K.> ...
%! knnsearch(ones (4, 5), ones (1,5), "K" ,0)
%!error<knnsearch: Invalid value of Minkowski Exponent.> ...
%! knnsearch(ones (4, 5), ones (1,5),"P",-2)
%!error<knnsearch: Invalid value in cov, cov can only be given for mahalanobis distance.> ...
%! knnsearch(ones (4, 5), ones (1, 5), "cov", ["some" "some"])
%!error<knnsearch: Invalid value in cov, cov can only be given for mahalanobis distance.> ...
%! knnsearch(ones (4, 5), ones (1, 5), "cov", ones(4,5), "distance", "euclidean")
%!error<knnsearch: Invalid value in Scale or the size of scale.> ...
%! knnsearch(ones (4, 5), ones (1, 5), "scale", ones(4,5), "distance", "euclidean")
%!error<knnsearch: Invalid value of bucketsize.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"bucketsize",-1)
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"NSmethod", "kdtree", "distance","cosine")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"NSmethod", "kdtree", "distance","mahalanobis")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"NSmethod", "kdtree", "distance","correlation")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"NSmethod", "kdtree", "distance","seuclidean")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"NSmethod", "kdtree", "distance","spearman")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"NSmethod", "kdtree", "distance","hamming")
%!error<knnsearch: 'kdtree' cannot be used with the given distance metric.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"NSmethod", "kdtree", "distance","jaccard")
