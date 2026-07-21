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
## @deftypefn  {statistics} {@var{idx} =} kmedoids (@var{X}, @var{k})
## @deftypefnx {statistics} {[@var{idx}, @var{C}] =} kmedoids (@var{X}, @var{k})
## @deftypefnx {statistics} {[@var{idx}, @var{C}, @var{sumd}] =} kmedoids (@var{X}, @var{k})
## @deftypefnx {statistics} {[@var{idx}, @var{C}, @var{sumd}, @var{D}] =} kmedoids (@var{X}, @var{k})
## @deftypefnx {statistics} {[@var{idx}, @var{C}, @var{sumd}, @var{D}, @var{midx}] =} kmedoids (@var{X}, @var{k})
## @deftypefnx {statistics} {[@var{idx}, @var{C}, @var{sumd}, @var{D}, @var{midx}, @var{info}] =} kmedoids (@var{X}, @var{k})
## @deftypefnx {statistics} {[@dots{}] =} kmedoids (@dots{}, @var{name}, @var{value})
##
## Partition observations into @var{k} clusters using the k-medoids algorithm.
##
## @code{@var{idx} = kmedoids (@var{X}, @var{k})} partitions the @math{N*P}
## numeric matrix @var{X} into @var{k} clusters, each represented by one of the
## observations (its @emph{medoid}), and returns the @math{N*1} vector
## @var{idx} of cluster indices.  Rows of @var{X} correspond to observations and
## columns correspond to features or variables.  Unlike @code{kmeans}, whose
## centroids are the mean of each cluster, a medoid is an actual data point,
## which makes k-medoids more robust to outliers and applicable to any distance
## metric.
##
## @code{[@var{idx}, @var{C}, @var{sumd}, @var{D}, @var{midx}, @var{info}] =
## kmedoids (@dots{})} returns additional results:
##
## @multitable @columnfractions 0.18 0.8
## @item @var{C} @tab a @math{k*P} matrix with the coordinates of the @var{k}
## medoids, one per row (@code{@var{C} = @var{X}(@var{midx},:)}).
##
## @item @var{sumd} @tab a @math{k*1} vector with the within-cluster sum of the
## distances from each point to its cluster medoid, measured with the selected
## metric.
##
## @item @var{D} @tab an @math{N*k} matrix with the distance from every
## observation to every medoid.
##
## @item @var{midx} @tab a @math{k*1} vector with the row indices into @var{X}
## of the @var{k} medoids.
##
## @item @var{info} @tab a scalar structure with the fields @qcode{'algorithm'},
## @qcode{'start'}, @qcode{'distance'}, @qcode{'iterations'}, and
## @qcode{'bestReplicate'} describing the chosen run.
## @end multitable
##
## Additional parameters can be specified by @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Distance'} @tab the distance metric, one of
## @qcode{'sqeuclidean'} (default), @qcode{'euclidean'}, @qcode{'seuclidean'},
## @qcode{'cityblock'}, @qcode{'minkowski'}, @qcode{'chebychev'},
## @qcode{'cosine'}, @qcode{'correlation'}, @qcode{'hamming'},
## @qcode{'jaccard'}, @qcode{'spearman'}, @qcode{'mahalanobis'}, or a custom
## distance function handle accepted by @code{pdist2}.
##
## @item @qcode{'Algorithm'} @tab the optimization algorithm, either
## @qcode{'pam'} (default) for Partitioning Around Medoids, which searches over
## all medoid/non-medoid swaps, or @qcode{'small'} for the faster Voronoi
## iteration that reassigns points and re-selects each cluster medoid until
## convergence.
##
## @item @qcode{'Start'} @tab the method used to choose the initial medoids:
## @qcode{'plus'} (default, k-means++), @qcode{'sample'} (a random subset of the
## observations), @qcode{'cluster'} (a preliminary pass on a subsample), or a
## @math{k*P} numeric matrix of starting medoid locations, each snapped to the
## nearest observation.  A @math{k*P*R} array supplies a separate start for each
## of @var{R} replicates.
##
## @item @qcode{'Replicates'} @tab a positive integer number of times to repeat
## the clustering, each with a new set of initial medoids; the solution with the
## lowest total sum of distances is returned.  The default is 1, or the size of
## the third dimension of a numeric @qcode{'Start'}.
##
## @item @qcode{'Options'} @tab a structure, as created by @code{statset}, whose
## @qcode{'MaxIter'} field caps the number of iterations (default 100).
## @end multitable
##
## @seealso{kmeans, linkage, pdist2, dbscan}
## @end deftypefn

function [idx, C, sumd, D, midx, info] = kmedoids (X, k, varargin)

  ## Check number of input arguments
  if (nargin < 2)
    error ("kmedoids: too few input arguments.");
  endif

  ## Validate X
  if (! isnumeric (X) || ! isreal (X) || ndims (X) != 2 || isempty (X))
    error ("kmedoids: X must be a nonempty real numeric matrix.");
  endif
  [N, P] = size (X);

  ## Validate k
  if (! isscalar (k) || ! isnumeric (k) || ! isreal (k) || k < 1
                     || fix (k) != k)
    error ("kmedoids: K must be a positive integer scalar.");
  endif
  if (k > N)
    error ("kmedoids: K cannot exceed the number of observations in X.");
  endif

  ## Defaults
  distance   = "sqeuclidean";
  algorithm  = "pam";
  start      = "plus";
  replicates = [];
  maxiter    = 100;

  ## Parse Name-Value pairs
  if (mod (numel (varargin), 2) != 0)
    error ("kmedoids: each NAME must be followed by a VALUE.");
  endif
  while (numel (varargin) > 0)
    name = varargin{1};
    val  = varargin{2};
    if (! ischar (name))
      error ("kmedoids: optional argument names must be strings.");
    endif
    switch (tolower (name))
      case "distance"
        distance = val;
      case "algorithm"
        algorithm = tolower (val);
      case "start"
        start = val;
      case "replicates"
        replicates = val;
      case "options"
        if (isstruct (val) && isfield (val, "MaxIter")
                           && ! isempty (val.MaxIter))
          maxiter = val.MaxIter;
        endif
      case "onlinephase"
        ## Accepted for compatibility; the medoid update is already exact.
      otherwise
        error ("kmedoids: unknown parameter name '%s'.", name);
    endswitch
    varargin(1:2) = [];
  endwhile

  ## Resolve the distance metric and map it onto a pdist2 metric
  metrics = {"sqeuclidean", "euclidean", "seuclidean", "cityblock", ...
             "minkowski", "chebychev", "cosine", "correlation", "hamming", ...
             "jaccard", "spearman", "mahalanobis"};
  if (ischar (distance))
    dname = tolower (distance);
    if (! any (strcmp (dname, metrics)))
      error ("kmedoids: unsupported distance metric '%s'.", distance);
    endif
    if (strcmp (dname, "sqeuclidean"))
      pmetric = "squaredeuclidean";
    else
      pmetric = dname;
    endif
  elseif (is_function_handle (distance))
    dname   = distance;
    pmetric = distance;
  else
    error ("kmedoids: DISTANCE must be a metric name or a function handle.");
  endif

  ## Validate the algorithm
  switch (algorithm)
    case {"pam", "small"}
      ## supported
    case {"clara", "large"}
      error ("kmedoids: the '%s' algorithm is not implemented.", algorithm);
    otherwise
      error ("kmedoids: unknown algorithm '%s'.", algorithm);
  endswitch

  ## Resolve the start method and the number of replicates
  if (isnumeric (start))
    if (columns (start) != P)
      error ("kmedoids: numeric START must have the same columns as X.");
    endif
    if (rows (start) != k)
      error ("kmedoids: numeric START must have K rows.");
    endif
    startmode = "numeric";
    rep_start = size (start, 3);
    if (isempty (replicates))
      replicates = rep_start;
    elseif (replicates != rep_start)
      error ("kmedoids: REPLICATES must equal the pages of a numeric START.");
    endif
  else
    startmode = tolower (start);
    if (! any (strcmp (startmode, {"plus", "sample", "cluster"})))
      error ("kmedoids: unknown START '%s'.", start);
    endif
    if (isempty (replicates))
      replicates = 1;
    endif
  endif
  if (! isscalar (replicates) || ! isnumeric (replicates) || replicates < 1
                              || fix (replicates) != replicates)
    error ("kmedoids: REPLICATES must be a positive integer scalar.");
  endif

  ## Precompute the full pairwise distance matrix in the chosen metric.  For
  ## 'sqeuclidean' this holds squared distances, so sums match MATLAB.
  Dall = pdist2 (X, X, pmetric);

  ## Run the requested number of replicates and keep the cheapest solution.
  best_cost = Inf;
  best_M    = [];
  best_iter = 0;
  best_rep  = 1;
  for rep = 1:replicates
    switch (startmode)
      case "numeric"
        M = snap_to_data (X, start(:,:,rep));
      case "sample"
        M = randperm (N, k)(:);
      case "plus"
        M = kpp_init (Dall, k);
      case "cluster"
        sub  = randperm (N, max (k, ceil (N / 10)));
        Msub = voronoi_step (Dall(sub, sub), randperm (numel (sub), k)(:), ...
                             maxiter);
        M    = sub(Msub)(:);
    endswitch

    if (strcmp (algorithm, "pam"))
      [M, iter] = pam_swap (Dall, M, maxiter);
    else
      [M, iter] = voronoi_step (Dall, M, maxiter);
    endif

    cost = sum (min (Dall(:, M), [], 2));
    if (cost < best_cost)
      best_cost = cost;
      best_M    = M;
      best_iter = iter;
      best_rep  = rep;
    endif
  endfor
  M = best_M;

  ## Assemble the outputs from the winning medoid set.
  D            = Dall(:, M);
  [mind, idx]  = min (D, [], 2);
  midx         = M(:);
  C            = X(midx, :);
  sumd         = accumarray (idx, mind, [k, 1]);
  info         = struct ("algorithm", algorithm, "start", startmode, ...
                         "distance", dname, "iterations", best_iter, ...
                         "bestReplicate", best_rep);

endfunction

## Snap each starting location to the nearest observation (Euclidean)
function M = snap_to_data (X, coords)
  k = rows (coords);
  M = zeros (k, 1);
  for j = 1:k
    [~, M(j)] = min (sum ((X - coords(j,:)) .^ 2, 2));
  endfor
endfunction

## k-means++ style seeding on the precomputed distance matrix
function M = kpp_init (Dall, k)
  N    = rows (Dall);
  M    = zeros (k, 1);
  M(1) = randi (N);
  d    = Dall(:, M(1));
  for i = 2:k
    total = sum (d);
    if (total <= 0)
      M(i) = randi (N);           # all remaining points coincide
    else
      M(i) = find (cumsum (d) >= rand * total, 1);
    endif
    d = min (d, Dall(:, M(i)));
  endfor
endfunction

## Partitioning Around Medoids: keep the best cost-reducing swap until none help
function [M, iter] = pam_swap (Dall, M, maxiter)
  N       = rows (Dall);
  k       = numel (M);
  curcost = sum (min (Dall(:, M), [], 2));
  iter    = 0;
  do
    iter    += 1;
    bestcost = curcost;
    bj = 0;
    bh = 0;
    ismed     = false (N, 1);
    ismed(M)  = true;
    for j = 1:k
      for h = 1:N
        if (ismed(h))
          continue;
        endif
        Mtry    = M;
        Mtry(j) = h;
        c = sum (min (Dall(:, Mtry), [], 2));
        if (c < bestcost)
          bestcost = c;
          bj = j;
          bh = h;
        endif
      endfor
    endfor
    improved = (bj > 0);
    if (improved)
      M(bj)   = bh;
      curcost = bestcost;
    endif
  until (! improved || iter >= maxiter)
endfunction

## Voronoi iteration ('small'): reassign, then re-select each cluster's medoid
function [M, iter] = voronoi_step (Dall, M, maxiter)
  k    = numel (M);
  iter = 0;
  do
    iter += 1;
    oldM     = M;
    [~, lab] = min (Dall(:, M), [], 2);
    for j = 1:k
      members = find (lab == j);
      if (isempty (members))
        continue;                 # keep the medoid of an empty cluster
      endif
      [~, loc] = min (sum (Dall(members, members), 1));
      M(j)     = members(loc);
    endfor
  until (isequal (M, oldM) || iter >= maxiter)
endfunction

%!demo
%! ## Cluster three noisy blobs and mark the medoids.
%! X = [randn(20,2)*0.4 + 3; randn(20,2)*0.4; randn(20,2)*0.4 + [3 -3]];
%! [idx, C] = kmedoids (X, 3);
%! gscatter (X(:,1), X(:,2), idx);
%! hold on;
%! plot (C(:,1), C(:,2), "kp", "MarkerSize", 14, "MarkerFaceColor", "y");
%! hold off;
%! title ("kmedoids: three clusters with their medoids");

## Exact anchor: numeric Start is RNG-free, verified against MATLAB R2024a
%!test
%! X = [1 1; 1.2 0.8; 0.8 1.1; 10 10; 10.2 9.8; 9.8 10.1; 1 10; 1.1 10.2; ...
%!      0.9 9.8];
%! S = [1 1; 10 10; 1 10];
%! [idx, C, sumd, D, midx] = kmedoids (X, 3, "Start", S);
%! assert_equal (idx, [1; 1; 1; 2; 2; 2; 3; 3; 3]);
%! assert_equal (midx, [1; 4; 7]);
%! assert_equal (C, [1 1; 10 10; 1 10]);
%! assert_equal (sumd, [0.13; 0.13; 0.10], 1e-12);
%! assert_equal (C, X(midx,:));
%! assert_equal (D(6,:), [160.25 0.05 77.45], 1e-12);
%! assert_equal (size (D), [9 3]);

## Default path (sqeuclidean, k-means++): partition and medoid set are unique
%!test
%! X = [1 1; 1.2 0.8; 0.8 1.1; 10 10; 10.2 9.8; 9.8 10.1; 1 10; 1.1 10.2; ...
%!      0.9 9.8];
%! [idx, C, sumd, D, midx, info] = kmedoids (X, 3);
%! assert_equal (sort (midx), [1; 4; 7]);
%! assert_equal (sort (sumd), [0.10; 0.13; 0.13], 1e-12);
%! assert_equal (arrayfun (@(c) numel (unique (idx(idx == c))), 1:3), [1 1 1]);
%! assert_equal (info.algorithm, "pam");
%! assert_equal (info.start, "plus");
%! assert_equal (info.distance, "sqeuclidean");

## Euclidean metric: sumd is the un-squared distance (numeric Start, exact)
%!test
%! X = [1 1; 1.2 0.8; 0.8 1.1; 10 10; 10.2 9.8; 9.8 10.1; 1 10; 1.1 10.2; ...
%!      0.9 9.8];
%! S = [1 1; 10 10; 1 10];
%! [idx, C, sumd, D, midx] = kmedoids (X, 3, "Start", S, "Distance", ...
%!                                     "euclidean");
%! assert_equal (midx, [1; 4; 7]);
%! assert_equal (sumd, [sqrt(0.08)+sqrt(0.05); sqrt(0.08)+sqrt(0.05); ...
%!                      2*sqrt(0.05)], 1e-12);

## Cityblock metric: sumd is the L1 distance (numeric Start, exact)
%!test
%! X = [1 1; 1.2 0.8; 0.8 1.1; 10 10; 10.2 9.8; 9.8 10.1; 1 10; 1.1 10.2; ...
%!      0.9 9.8];
%! S = [1 1; 10 10; 1 10];
%! [idx, C, sumd, D, midx] = kmedoids (X, 3, "Start", S, "Distance", ...
%!                                     "cityblock");
%! assert_equal (midx, [1; 4; 7]);
%! assert_equal (sumd, [0.7; 0.7; 0.6], 1e-12);
%! assert_equal (D(1,:), [0 18 9], 1e-12);

## The 'small' Voronoi algorithm finds the same partition on separated data
%!test
%! X = [1 1; 1.2 0.8; 0.8 1.1; 10 10; 10.2 9.8; 9.8 10.1; 1 10; 1.1 10.2; ...
%!      0.9 9.8];
%! S = [1 1; 10 10; 1 10];
%! [idx, C, sumd, D, midx, info] = kmedoids (X, 3, "Start", S, "Algorithm", ...
%!                                           "small");
%! assert_equal (midx, [1; 4; 7]);
%! assert_equal (sumd, [0.13; 0.13; 0.10], 1e-12);
%! assert_equal (info.algorithm, "small");

## Two clusters, a custom number of replicates runs and returns a valid result
%!test
%! X = [randn(10,2) - 5; randn(10,2) + 5];
%! [idx, C, sumd, D, midx] = kmedoids (X, 2, "Replicates", 3);
%! assert_equal (numel (idx), 20);
%! assert (all (idx >= 1 & idx <= 2));
%! assert_equal (C, X(midx,:));
%! assert_equal (size (D), [20 2]);

## k == number of observations places one medoid on every point
%!test
%! X = [0; 5; 9];
%! [idx, C, sumd, D, midx] = kmedoids (X, 3, "Start", X);
%! assert_equal (sort (midx), [1; 2; 3]);
%! assert_equal (sumd, [0; 0; 0], 1e-12);

## Test input validation
%!error <kmedoids: too few input arguments.> kmedoids (1)
%!error <kmedoids: X must be a nonempty real numeric matrix.> kmedoids ([], 2)
%!error <kmedoids: X must be a nonempty real numeric matrix.> kmedoids ("a", 2)
%!error <kmedoids: K must be a positive integer scalar.> kmedoids (ones (4,2), 0)
%!error <kmedoids: K must be a positive integer scalar.> kmedoids (ones (4,2), 1.5)
%!error <kmedoids: K cannot exceed the number of observations in X.> ...
%! kmedoids (ones (3,2), 4)
%!error <kmedoids: each NAME must be followed by a VALUE.> ...
%! kmedoids (ones (4,2), 2, "Distance")
%!error <kmedoids: unknown parameter name 'foo'.> ...
%! kmedoids (ones (4,2), 2, "foo", "bar")
%!error <kmedoids: unsupported distance metric 'taxicab'.> ...
%! kmedoids (ones (4,2), 2, "Distance", "taxicab")
%!error <kmedoids: the 'clara' algorithm is not implemented.> ...
%! kmedoids (ones (4,2), 2, "Algorithm", "clara")
%!error <kmedoids: unknown algorithm 'foo'.> ...
%! kmedoids (ones (4,2), 2, "Algorithm", "foo")
%!error <kmedoids: unknown START 'middle'.> ...
%! kmedoids (ones (4,2), 2, "Start", "middle")
%!error <kmedoids: numeric START must have K rows.> ...
%! kmedoids (ones (4,2), 2, "Start", [1 1; 2 2; 3 3])
%!error <kmedoids: numeric START must have the same columns as X.> ...
%! kmedoids (ones (4,2), 2, "Start", [1 1 1; 2 2 2])
%!error <kmedoids: REPLICATES must be a positive integer scalar.> ...
%! kmedoids (ones (4,2), 2, "Replicates", 0)
