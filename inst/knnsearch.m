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
## @deftypefn  {statistics} {@var{idx} =} knnsearch (@var{x}, @var{y})
## @deftypefnx {statistics} {@var{idx} =} knnsearch (@var{x}, @var{y}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{idx}, @var{D}] =} knnsearch (@dots{})
##
## Find k-nearest neighbors using input data
## @itemize
## @item
## @code{@var{idx} = knnsearch (@var{x}, @var{y})} finds "K" nearest neighbors
## in @var{x} for @var{y}. It returns @var{idx} which contains indices of "K"
## nearest neighbors of each row of @var{y}, If not specified @qcode{@var{K} = 1}.
## @item
## @code{[@var{idx}, @var{dist}] = knnsearch (@var{x}, @var{y})} returns the
## "K" nearest neighbour in @var{x} for each @var{y} with distances returned
## in @var{dist}.
## @end itemize
##
## additional arguments can be given as name-value pairs. such as value for @var{K}
## and distance metric can be specified to be used in search.
##
## @itemize
## @item
## @code{x} must be an @math{NxP} numeric matrix of input data where rows correspond
## to observations and columns correspond to features or variables.
## @item
## @code{y} is an @math{MxP} numeric matrix with query points. @var{y} must have
## same numbers of column as @var{x}.
## @emph{ additional parameters that can be passed as name value pairs :}
## @item
## @code{K} is the number of nearest neighbours to be considered in kNN search
##          Default value of @qcode{@var{k} = 1}.
## @item
## @code{exponent} is the minkowski distance exponent. Default is @qcode{@var{P} = 2}.
##
## @item
## @code{scale} is scale for standardized euclidean distance. Default is @qcode{@var{scale} = []}.
## @item
## @code{cov}   is the cov matrix for computing mahalanobis distance. Default is @qcode{@var{cov} = []}.
## @item
## @code{bucketsize} is maximum number of data points per leaf node of kd-tree.
## if NSmethod is 'kdtree', Default is @qcode{@var{bucketsize} = 50}.
## @item
## @code{sortindices}  is the flag to indicate wheather the returned indices are
## to be sorted by distance. Defaul is @qcode{@var{sortindices} = true}.
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
## default is @qcode{@var{P} = 2}.
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
## @item
## @code{NSMethod}   is nearest neighbour search method. Default is
## @qcode{@var{NSMethod} = 'exhaustive'}. if specified as @qcode{@var{NSMethod} = "kdtree"}
## function builds a kdtree of @var{x} and then searches for query points in @var{y}.
## @item
## @code{@var{includeties}} is a flag to indicate if the returned values should
## contain the indices that have same distance as k'th neighbour. Default is
## @qcode{@var{includeties} = true}.
## @end itemize
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
  NSMethod = "exhaustive";  # Nearest neighbor search method
  incl_ties = false;        # Include ties for distance with kth neighbor

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
        incl_ties = varargin{2};
      otherwise
        error ("knnsearch: invalid NAME in optional pairs of arguments.");
    endswitch
    varargin(1:2) = [];
  endwhile

  ## check input parameters
  if ( !isscalar(K) || !isnumeric(K) || K < 1 || K != round (K))
    error ("knnsearch: Invalid value of k.");
  endif
  if ( !isscalar(P) || !isnumeric(P) || P < 0)
    error ("knnsearch: Invalid value of Minkowski Exponent.");
  endif
  if ( !isempty(C))
    if ( !strcmp(Distance,"mahalanobis") || !ismatrix(C) || !isnumeric(C))
      error ("knnsearch: Invalid value in cov, cov can only be given for mahalanobis distance.");
    endif
  endif
  if ( !isempty(S))
    if ( !isscalar(S) || any(S) < 0 || numel(S) != rows(X) || !strcmpi (Distance,"seuclidean"))
      error ("knnsearch: Invalid value in Scale or the size of scale.");
    endif
  endif
  if ( !isscalar(BS) || BS < 0)
    error ("knnsearch: Invalid value of bucketsize.");
  endif

  ## check for NSMethod
  if (strcmpi(NSMethod,"kdtree"))
    ## build kdtree and search the query point
    ret = buildkdtree (X, BS);
    NN  = findkdtree (ret, Y, K, Distance, P, S, C);
    D   = calc_dist(X(NN,:), Y, Distance, P, S, C);
    ## check for ties and sortindices
    if (incl_ties)
      if (SI)
        sorted_D = sortrows ([NN, D],2);
        dist = sorted_D(:,2)';
        idx  = sorted_D(:,1)';
      else
        sorted_D = sortrows ([NN, D],2);
        dist = sorted_D(:,2)';
        idx  = sorted_D(:,1)';
      endif
    else
      sorted_D = sortrows ([NN, D],2);
      dist = sorted_D(1:K,2)';
      idx  = sorted_D(1:K,1)';
    endif
  else
    ## calculate distance and search by exhaustive
    D = calc_dist(X, Y, Distance, P, S, C);
    D = reshape (D, size (Y, 1), size (X, 1));
    if (K == 1)
      [dist, idx] = min (D, [], 2);
    else
      if (incl_ties)
        [dist, idx] = sort (D, 2);
        kth_dist = dist(K);
        tied_idx = (dist <= kth_dist);
        dist = dist(:,tied_idx);
        idx  = idx(:,tied_idx);
      else
        [dist, idx] = sort (D, 2);
        dist = dist(:,1:K);
        idx = idx(:,1:K);
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
        ret = struct ('point', r(1), 'dimen', d);
    else
        mid = ceil (count / 2);
        ret = struct ('point', r(mid), 'dimen', d);
        d = mod (d,dimen)+1;
        ## Build left sub tree
        if (mid > 1)
            left = r (1:mid-1);
            left_points = x (left,d);
            [val, left_idx] = sort (left_points);
            leftr = left (left_idx);
            ret.left = buildkdtree_recur (x, leftr, d);
        endif
        ## Build right sub tree
        if (count > mid)
            right = r (mid+1:count);
            right_points = x (right,d);
            [val, right_idx] = sort (right_points);
            rightr = right (right_idx);
            ret.right = buildkdtree_recur (x, rightr, d);
        endif
    endif
endfunction

## wrapper function for buildkdtree_recur
function ret = buildkdtree (x, BS)
    [val, r] = sort (x(:,1));
    ret = struct ('data',x,'root', buildkdtree_recur (x, r, 1, BS));
endfunction

function farthest = kdtree_cand_farthest (x, p, cand, distance, P, S, C)
    [val, index] = max(calc_dist (x, p, distance, P, S, C)(cand));
    farthest = cand (index);
endfunction

## function to insert into NN list
function inserted = kdtree_cand_insert (x, p, cand, k, point, distance, P, S, C)
    if (length (cand) < k)
        inserted = [cand; point];
    else
        farthest = kdtree_cand_farthest(x, p, cand, distance, P, S, C);
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
function neighbours = findkdtree_recur (x, node, p, neighbours, k, distance, P, S, C)
    point = node.point;
    d = node.dimen;
    if (x(point,d) > p(d))
        ## Search in left sub tree
        if (isfield(node, 'left'))
            neighbours = findkdtree_recur (x, node.left, p, neighbours, k, distance, P, S, C);
        endif
        ## Add current point if neccessary
        farthest = kdtree_cand_farthest (x, p, neighbours, distance, P, S, C);
        if (length(neighbours) < k || calc_dist (x(point,:), p, distance, P, S, C) <= calc_dist (x(farthest,:), p, distance, P, S, C))
            neighbours = kdtree_cand_insert (x, p, neighbours, k, point, distance, P, S, C);
        endif
        ## Search in right sub tree if neccessary
        farthest = kdtree_cand_farthest (x, p, neighbours, distance, P, S, C);
        radius = calc_dist (x(farthest,:), p, distance, P, S, C);
        if (isfield (node, 'right') && (length(neighbours) < k || p(d) + radius > x(point,d)))
            neighbours = findkdtree_recur (x, node.right, p, neighbours, k, distance, P, S, C);
        endif
    else
        ## Search in right sub tree
        if (isfield (node, 'right'))
            neighbours = findkdtree_recur (x, node.right, p, neighbours, k, distance, P, S, C);
        endif
        ## Add current point if neccessary
        farthest = kdtree_cand_farthest (x, p, neighbours, distance, P, S, C);
        if (length (neighbours) < k || calc_dist (x(point,:), p, distance, P, S, C) <= calc_dist (x(farthest,:), p, distance, P, S, C))
            neighbours = kdtree_cand_insert (x, p, neighbours, k, point, distance, P, S, C);
        endif
        ## Search in left sub tree if neccessary
        farthest = kdtree_cand_farthest (x, p, neighbours, distance, P, S, C);
        radius = calc_dist (x(farthest,:), p, distance, P, S, C);
        if (isfield (node, 'left') && (length (neighbours) < k || p(d) - radius <= x(point,d)))
            neighbours = findkdtree_recur (x, node.left, p, neighbours, k, distance, P, S, C);
        endif
    endif
endfunction

## wrapper function for findkdtree_recur
function neighbours = findkdtree (tree, p, k, distance, P, S, C)
    x = tree.data;
    root = tree.root;
    neighbours = findkdtree_recur (x, root, p, [], k, distance, P, S, C);
endfunction

## Test output
%!shared x, y
%! x = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
%! y = [1, 2, 2, 3; 2, 3, 3, 4];
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "euclidean");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "euclidean", "k", 2);
%! assert (idx, [1, 2; 1, 2]);
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
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "minkowski", "p", 3);
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * 1.259921049894873, 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "cityblock");
%! assert (idx, [1; 1]);
%! assert (D, [2; 2]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "manhattan");
%! assert (idx, [1; 1]);
%! assert (D, [2; 2]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "chebychev");
%! assert (idx, [1; 1]);
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
%! a = [1, 5;1,	2;2,	2;1.5,	1.5;5,	1;2	-1.34;1,	-3;4,	-4;-3,	1;8,	9];
%! b = [1, 1];
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree","includeties",true);
%! assert (idx, [4, 2, 3, 6, 1, 9, 7, 5]);
%! assert (D, [0.7071, 1.0000, 1.4142, 2.5447, 4.0000, 4.0000, 4.0000, 4.0000],1e-4);
%!test
%! a = [1, 5;1,	2;2,	2;1.5,	1.5;5,	1;2	-1.34;1,	-3;4,	-4;-3,	1;8,	9];
%! b = [1, 1];
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree","includeties",false);
%! assert (idx, [4, 2, 3, 6, 1]);
%! assert (D, [0.7071, 1.0000, 1.4142, 2.5447, 4.0000],1e-4);
%!test
%! a = [1, 5;1,	2;2,	2;1.5,	1.5;5,	1;2	-1.34;1,	-3;4,	-4;-3,	1;8,	9];
%! b = [1, 1];
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "exhaustive","includeties",false);
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
%! load fisheriris
%! a = meas;
%! b = max(meas);
%! [idx, D] = knnsearch (a, b, "K", 5, "NSMethod", "kdtree");
%! assert (idx, [118, 132, 110, 106, 136]);
%! assert (D, [0.7280, 0.9274, 1.3304, 1.5166, 1.6371],1e-4);


## Test input validation
%!error<knnsearch: too few input arguments.> knnsearch (1)
%!error<knnsearch: number of columns in X and Y must match.> ...
%! knnsearch (ones (4, 5), ones (4))
%!error<knnsearch: invalid NAME in optional pairs of arguments.> ...
%! knnsearch (ones (4, 2), ones (3, 2), "Distance", "euclidean", "some", "some")
%!error<knnsearch: Invalid value of k.> knnsearch(ones (4, 5), ones (1,5), "K" ,0)
%!error<knnsearch: Invalid value of Minkowski Exponent.> knnsearch(ones (4, 5), ones (1,5),"P",-2)
%!error<knnsearch: Invalid value in cov, cov can only be given for mahalanobis distance.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"cov",["some" "some"])
%!error<knnsearch: Invalid value in cov, cov can only be given for mahalanobis distance.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"cov",ones(4,5),"distance","euclidean")
%!error<knnsearch: Invalid value in Scale or the size of scale.> ...
%! knnsearch(ones (4, 5), ones (1, 5),"scale",ones(4,5),"distance","euclidean")
%!error<knnsearch: Invalid value of bucketsize.> knnsearch(ones (4, 5), ones (1, 5),"bucketsize",-1)
