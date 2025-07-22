## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{idx} =} rangesearch (@var{X}, @var{Y}, @var{r})
## @deftypefnx {statistics} {[@var{idx}, @var{D}] =} rangesearch (@var{X}, @var{Y}, @var{r})
## @deftypefnx {statistics} {[@dots{}] =} rangesearch (@dots{}, @var{name}, @var{value})
##
## Find all neighbors within specified distance from input data.
##
## @code{@var{idx} = rangesearch (@var{X}, @var{Y}, @var{r})} returns all the
## points in @var{X} that are within distance @var{r} from the points in @var{Y}.
## @var{X} must be an @math{NxP} numeric matrix of input data, where rows
## correspond to observations and columns correspond to features or variables.
## @var{Y} is an @math{MxP} numeric matrix with query points, which must have
## the same numbers of column as @var{X}.  @var{r} must be a nonnegative scalar
## value.  @var{idx} is an @math{Mx1} cell array, where @math{M} is the number
## of observations in @var{Y}.  The vector @qcode{@var{Idx}@{j@}} contains the
## indices of observations (rows) in @var{X} whose distances to
## @qcode{@var{Y}(j,:)} are not greater than @var{r}.
##
## @code{[@var{idx}, @var{D}] = rangesearch (@var{X}, @var{Y}, @var{r})} also
## returns the distances, @var{D}, which correspond to the points in @var{X}
## that are within distance @var{r} from the points in @var{Y}.  @var{D} is an
## @math{Mx1} cell array, where @math{M} is the number of observations in
## @var{Y}.  The vector @qcode{@var{D}@{j@}} contains the distances of
## observations (rows) in @var{X} whose distances to @qcode{@var{Y}(j,:)} are
## not greater than @var{r}.
##
## Additional parameters can be specified by @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
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
## @qcode{"IncludeTies"} flag is true, @code{rangesearch} always sorts the
## returned indices.
##
## @item @qcode{"Distance"} @tab @tab is the distance metric used by
## @code{rangesearch} as specified below:
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
## by @code{rangesearch} as specified below.
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
## @seealso{knnsearch, pdist2}
## @end deftypefn

function [idx, dist] = rangesearch (X, Y, r, varargin)

  ## Check input data
  if (nargin < 3)
    error ("rangesearch: too few input arguments.");
  endif

  if (size (X, 2) != size (Y, 2))
	  error ("rangesearch: number of columns in X and Y must match.");
  endif

  if (! isscalar (r) || ! isnumeric (r) || r < 0)
    error ("rangesearch: radius r must be a nonnegative scalar.");
  endif

  ## Add default values
  P = 2;                    # Exponent for Minkowski distance
  S = [];                   # Scale for the standardized Euclidean distance
  C = [];                   # Covariance matrix for Mahalanobis distance
  BS = 50;                  # Maximum number of points per leaf node for Kd-tree
  SI = true;                # Sort returned indices according to distance
  Distance = "euclidean";   # Distance metric to be used
  NSMethod = [];            # Nearest neighbor search method
  DistParameter = [];       # Distance parameter for pdist2

  ## Parse additional parameters in Name/Value pairs
  PSC = 0;
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))
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
      otherwise
        error ("rangesearch: invalid NAME in optional pairs of arguments.");
    endswitch
    varargin(1:2) = [];
  endwhile

  ## Check input parameters
  if (PSC > 1)
    error ("rangesearch: only a single distance parameter can be defined.");
  endif
  if (! isscalar (P) || ! isnumeric (P) || P <= 0)
    error ("rangesearch: invalid value of Minkowski Exponent.");
  endif
  if (! isempty (S))
    if (any (S) < 0 || numel (S) != columns (X)
                    || ! strcmpi (Distance, "seuclidean"))
      error ("rangesearch: invalid value in Scale or the size of Scale.");
    endif
  endif
  if (! isempty (C))
    if (! strcmp (Distance, "mahalanobis") || ! ismatrix (C) || ! isnumeric (C))
      error (strcat (["rangesearch: invalid value in Cov, Cov can only"], ...
                     [" be given for mahalanobis distance."]));
    endif
  endif
  if (! isscalar (BS) || BS < 0)
    error ("rangesearch: invalid value of bucketsize.");
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
     || strcmpi (Distance, "minkowski") || strcmpi (Distance, "chebychev")))
      NSMethod = "kdtree";
    else
      NSMethod = "exhaustive";
    endif
  else
    ## Check if kdtree can be used
    if (strcmpi (NSMethod, "kdtree") && ! (strcmpi (Distance, "euclidean")
     || strcmpi (Distance, "cityblock") || strcmpi (Distance, "minkowski")
     || strcmpi (Distance, "chebychev")))
      error (strcat (["rangesearch: 'kdtree' cannot be used with"], ...
                     [" the given distance metric."]));
    endif
  endif

  ## Check for NSMethod
  if (strcmpi (NSMethod, "kdtree"))
    ## Build kdtree and search the query point
    kdtree = __build_kdtree__ (1:size(X,1), 0, X, BS);
    ## Return all neighbors as cell
    dist = cell (rows (Y), 1);
    idx = cell (rows (Y), 1);
    for i = 1:rows (Y)
      [temp_idx, temp_D] = __search_kdtree__ (kdtree, Y(i,:), Inf, X, ...
																							Distance, DistParameter, true, r);
      if (SI)
        [sorted_D, sort_idx] = sort (temp_D);
        idx{i} = temp_idx(sort_idx);
        dist{i} = sorted_D;
      else
        idx{i} = temp_idx;
        dist{i} = temp_D;
      endif
    endfor
  else
    ## Calculate all distances
    dist = cell (rows (Y), 1);
    idx = cell (rows (Y), 1);
    for i = 1:rows (Y)
      D = pdist2 (X, Y(i,:), Distance, DistParameter);
      Didx_row = find (D <= r)';
      Dist_row = D(Didx_row)';
      if (SI)
        [S, I] = sort (Dist_row);
        Dist_row = Dist_row(I);
        Didx_row = Didx_row(I);
      endif
      dist{i} = Dist_row;
      idx{i} = Didx_row;
    endfor
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
    if (! (isscalar (distparam) && isnumeric (distparam)
                                && distparam > 0 && isfinite (distparam)))
      error (strcat("rangesearch.__search_kdtree__:", ...
                    " distparam must be a positive finite", ...
                    " scalar for minkowski."));
    endif
  else
    if (! isempty (distparam))
      error (strcat("rangesearch.__search_kdtree__:", ...
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
%! ## Generate 1000 random 2D points from each of five distinct multivariate
%! ## normal distributions that form five separate classes
%! N = 1000;
%! d = 10;
%! randn ("seed", 5);
%! X1 = mvnrnd (d * [0, 0], eye (2), 1000);
%! randn ("seed", 6);
%! X2 = mvnrnd (d * [1, 1], eye (2), 1000);
%! randn ("seed", 7);
%! X3 = mvnrnd (d * [-1, -1], eye (2), 1000);
%! randn ("seed", 8);
%! X4 = mvnrnd (d * [1, -1], eye (2), 1000);
%! randn ("seed", 8);
%! X5 = mvnrnd (d * [-1, 1], eye (2), 1000);
%! X = [X1; X2; X3; X4; X5];
%!
%! ## For each point in X, find the points in X that are within a radius d
%! ## away from the points in X.
%! Idx = rangesearch (X, X, d, "NSMethod", "exhaustive");
%!
%! ## Select the first point in X (corresponding to the first class) and find
%! ## its nearest neighbors within the radius d.  Display these points in
%! ## one color and the remaining points in a different color.
%! x = X(1,:);
%! nearestPoints = X (Idx{1},:);
%! nonNearestIdx = true (size (X, 1), 1);
%! nonNearestIdx(Idx{1}) = false;
%!
%! scatter (X(nonNearestIdx,1), X(nonNearestIdx,2))
%! hold on
%! scatter (nearestPoints(:,1),nearestPoints(:,2))
%! scatter (x(1), x(2), "black", "filled")
%! hold off
%!
%! ## Select the last point in X (corresponding to the fifth class) and find
%! ## its nearest neighbors within the radius d.  Display these points in
%! ## one color and the remaining points in a different color.
%! x = X(end,:);
%! nearestPoints = X (Idx{1},:);
%! nonNearestIdx = true (size (X, 1), 1);
%! nonNearestIdx(Idx{1}) = false;
%!
%! figure
%! scatter (X(nonNearestIdx,1), X(nonNearestIdx,2))
%! hold on
%! scatter (nearestPoints(:,1),nearestPoints(:,2))
%! scatter (x(1), x(2), "black", "filled")
%! hold off


## Test output
%!shared x, y, X, Y
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [2, 3, 4; 1, 4, 3];
%! X = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
%! Y = [1, 2, 2, 3; 2, 3, 3, 4];
%!test
%! [idx, D] = rangesearch (x, y, 4);
%! assert (idx, {[1; 4; 2]; [1; 4]});
%! assert (D, {[1.7321; 3.3166; 3.4641]; [2; 3.4641]}, 1e-4);
%!test
%! [idx, D] = rangesearch (x, y, 4, "NSMethod", "exhaustive");
%! assert (idx, {[1, 4, 2]; [1, 4]});
%! assert (D, {[1.7321, 3.3166, 3.4641]; [2, 3.4641]}, 1e-4);
%!test
%! [idx, D] = rangesearch (x, y, 4, "NSMethod", "kdtree");
%! assert (idx, {[1; 4; 2]; [1; 4]});
%! assert (D, {[1.7321; 3.3166; 3.4641]; [2; 3.4641]}, 1e-4);
%!test
%! [idx, D] = rangesearch (x, y, 4, "SortIndices", true);
%! assert (idx, {[1; 4; 2]; [1; 4]});
%! assert (D, {[1.7321; 3.3166; 3.4641]; [2; 3.4641]}, 1e-4);
%!test
%! [idx, D] = rangesearch (x, y, 4, "SortIndices", false);
%! assert (idx, {[1; 2; 4]; [1; 4]});
%! assert (D, {[1.7321; 3.4641; 3.3166]; [2; 3.4641]}, 1e-4);
%!test
%! [idx, D] = rangesearch (x, y, 4, "NSMethod", "exhaustive", ...
%!                         "SortIndices", false);
%! assert (idx, {[1, 2, 4]; [1, 4]});
%! assert (D, {[1.7321, 3.4641, 3.3166]; [2, 3.4641]}, 1e-4);
%!test
%! eucldist = @(v,m) sqrt(sumsq(repmat(v,rows(m),1)-m,2));
%! [idx, D] = rangesearch (x, y, 4, "Distance", eucldist);
%! assert (idx, {[1, 4, 2]; [1, 4]});
%! assert (D, {[1.7321, 3.3166, 3.4641]; [2, 3.4641]}, 1e-4);
%!test
%! eucldist = @(v,m) sqrt(sumsq(repmat(v,rows(m),1)-m,2));
%! [idx, D] = rangesearch (x, y, 4, "Distance", eucldist, ...
%!                         "NSMethod", "exhaustive");
%! assert (idx, {[1, 4, 2]; [1, 4]});
%! assert (D, {[1.7321, 3.3166, 3.4641]; [2, 3.4641]}, 1e-4);
%!test
%! [idx, D] = rangesearch (x, y, 1.5, "Distance", "seuclidean", ...
%!                         "NSMethod", "exhaustive");
%! assert (idx, {[1, 4, 2]; [1, 4]});
%! assert (D, {[0.6024, 1.0079, 1.2047]; [0.6963, 1.2047]}, 1e-4);
%!test
%! [idx, D] = rangesearch (x, y, 1.5, "Distance", "seuclidean", ...
%!                         "NSMethod", "exhaustive", "SortIndices", false);
%! assert (idx, {[1, 2, 4]; [1, 4]});
%! assert (D, {[0.6024, 1.2047, 1.0079]; [0.6963, 1.2047]}, 1e-4);
%!test
%! [idx, D] = rangesearch (X, Y, 4);
%! assert (idx, {[1; 2]; [1; 2; 3]});
%! assert (D, {[1.4142; 3.1623]; [1.4142; 1.4142; 3.1623]}, 1e-4);
%!test
%! [idx, D] = rangesearch (X, Y, 2);
%! assert (idx, {[1]; [1; 2]});
%! assert (D, {[1.4142]; [1.4142; 1.4142]}, 1e-4);
%!test
%! eucldist = @(v,m) sqrt(sumsq(repmat(v,rows(m),1)-m,2));
%! [idx, D] = rangesearch (X, Y, 4, "Distance", eucldist);
%! assert (idx, {[1, 2]; [1, 2, 3]});
%! assert (D, {[1.4142, 3.1623]; [1.4142, 1.4142, 3.1623]}, 1e-4);
%!test
%! [idx, D] = rangesearch (X, Y, 4, "SortIndices", false);
%! assert (idx, {[1; 2]; [1; 2; 3]});
%! assert (D, {[1.4142; 3.1623]; [1.4142; 1.4142; 3.1623]}, 1e-4);
%!test
%! [idx, D] = rangesearch (X, Y, 4, "Distance", "seuclidean", ...
%!                         "NSMethod", "exhaustive");
%! assert (idx, {[1, 2]; [1, 2, 3]});
%! assert (D, {[1.4142, 3.1623]; [1.4142, 1.4142, 3.1623]}, 1e-4);

## Test input validation
%!error<rangesearch: too few input arguments.> rangesearch (1)
%!error<rangesearch: too few input arguments.> rangesearch (ones (4, 5))
%!error<rangesearch: too few input arguments.> ...
%! rangesearch (ones (4, 5), ones (4))
%!error<rangesearch: number of columns in X and Y must match.> ...
%! rangesearch (ones (4, 5), ones (4), 1)
%!error<rangesearch: invalid NAME in optional pairs of arguments.> ...
%! rangesearch (ones (4, 2), ones (3, 2), 1, "Distance", "euclidean", "some", "some")
%!error<rangesearch: only a single distance parameter can be defined.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "scale", ones (1, 5), "P", 3)
%!error<rangesearch: invalid value of Minkowski Exponent.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "P", -2)
%!error<rangesearch: invalid value in Scale or the size of Scale.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "scale", ones(4,5), "distance", "euclidean")
%!error<rangesearch: invalid value in Cov, Cov can only be given for mahalanobis distance.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "cov", ["some" "some"])
%!error<rangesearch: invalid value in Cov, Cov can only be given for mahalanobis distance.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "cov", ones(4,5), "distance", "euclidean")
%!error<rangesearch: invalid value of bucketsize.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "bucketsize", -1)
%!error<rangesearch: 'kdtree' cannot be used with the given distance metric.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "NSmethod", "kdtree", "distance", "cosine")
%!error<rangesearch: 'kdtree' cannot be used with the given distance metric.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "NSmethod", "kdtree", "distance", "mahalanobis")
%!error<rangesearch: 'kdtree' cannot be used with the given distance metric.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "NSmethod", "kdtree", "distance", "correlation")
%!error<rangesearch: 'kdtree' cannot be used with the given distance metric.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "NSmethod", "kdtree", "distance", "seuclidean")
%!error<rangesearch: 'kdtree' cannot be used with the given distance metric.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "NSmethod", "kdtree", "distance", "spearman")
%!error<rangesearch: 'kdtree' cannot be used with the given distance metric.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "NSmethod", "kdtree", "distance", "hamming")
%!error<rangesearch: 'kdtree' cannot be used with the given distance metric.> ...
%! rangesearch (ones (4, 5), ones (1, 5), 1, "NSmethod", "kdtree", "distance", "jaccard")
