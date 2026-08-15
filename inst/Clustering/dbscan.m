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
## @deftypefn  {statistics} {@var{idx} =} dbscan (@var{X}, @var{epsilon}, @var{minpts})
## @deftypefnx {statistics} {@var{idx} =} dbscan (@var{D}, @var{epsilon}, @var{minpts}, @qcode{'Distance'}, @qcode{'precomputed'})
## @deftypefnx {statistics} {@var{idx} =} dbscan (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{idx}, @var{corepts}] =} dbscan (@dots{})
##
## Density-Based Spatial Clustering of Applications with Noise (DBSCAN).
##
## @code{@var{idx} = dbscan (@var{X}, @var{epsilon}, @var{minpts})} partitions
## the observations in the @math{N*P} numeric matrix @var{X} into clusters
## using the DBSCAN algorithm with neighborhood radius @var{epsilon} and
## minimum number of neighbors @var{minpts}.  Rows of @var{X} correspond to
## observations and columns correspond to features or variables.  @var{epsilon}
## must be a nonnegative scalar and @var{minpts} a positive integer scalar.
## @var{idx} is an @math{N*1} vector of cluster indices, numbered @math{1} to
## the number of clusters found; observations flagged as noise are assigned the
## value @math{-1}.
##
## A point is a @emph{core point} when at least @var{minpts} observations
## (@strong{including the point itself}) lie within distance @var{epsilon} of
## it.  Clusters grow from core points to every observation that is
## density-reachable from them; a non-core observation that lies within
## @var{epsilon} of a core point becomes a @emph{border point} and joins that
## point's cluster, while an observation that is neither core nor within reach
## of a core point is labelled noise.  A border point that is reachable from
## more than one cluster is assigned to the first cluster that reaches it,
## following the order of the observations in @var{X}.
##
## @code{@var{idx} = dbscan (@var{D}, @var{epsilon}, @var{minpts},
## @qcode{'Distance'}, @qcode{'precomputed'})} treats the @math{N*N} matrix
## @var{D} as a precomputed matrix of pairwise distances between observations,
## such as the output of @code{pdist2}; @qcode{@var{D}(i,j)} is the distance
## between observations @math{i} and @math{j}.
##
## @code{[@var{idx}, @var{corepts}] = dbscan (@dots{})} also returns an
## @math{N*1} logical vector @var{corepts} that is @qcode{true} for each
## observation that is a core point.
##
## Additional parameters can be specified by @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Distance'} @tab is the distance metric used to find neighbors,
## specified as one of the metrics accepted by @code{rangesearch}
## (@qcode{'euclidean'} by default, and also @qcode{'seuclidean'},
## @qcode{'cityblock'}, @qcode{'chebychev'}, @qcode{'minkowski'},
## @qcode{'mahalanobis'}, @qcode{'cosine'}, @qcode{'correlation'},
## @qcode{'spearman'}, @qcode{'hamming'}, @qcode{'jaccard'}, or a custom
## distance function handle), or the string @qcode{'precomputed'} to interpret
## the first input as a matrix of pairwise distances.
##
## @item @qcode{'P'} @tab is the Minkowski distance exponent, a positive scalar.
## This argument is only valid when the selected distance metric is
## @qcode{'minkowski'}.  By default it is 2.
##
## @item @qcode{'Scale'} @tab is the scale parameter for the standardized
## Euclidean distance, a nonnegative numeric vector of length equal to the
## number of columns in @var{X}.  This argument is only valid when the selected
## distance metric is @qcode{'seuclidean'}.
##
## @item @qcode{'Cov'} @tab is the covariance matrix for the mahalanobis
## distance, a positive definite matrix matching the number of columns in
## @var{X}.  This argument is only valid when the selected distance metric is
## @qcode{'mahalanobis'}.
## @end multitable
##
## @seealso{kmeans, rangesearch, pdist2, knnsearch}
## @end deftypefn

function [idx, corepts] = dbscan (X, epsilon, minpts, varargin)

  ## Check number of input arguments
  if (nargin < 3)
    error ("dbscan: too few input arguments.");
  endif

  ## Validate X
  if (! isnumeric (X) || ! isreal (X) || ndims (X) != 2 || isempty (X))
    error ("dbscan: X must be a nonempty real numeric matrix.");
  endif

  ## Validate epsilon
  if (! isscalar (epsilon) || ! isnumeric (epsilon) || ! isreal (epsilon)
                           || epsilon < 0)
    error ("dbscan: EPSILON must be a nonnegative scalar.");
  endif

  ## Validate minpts
  if (! isscalar (minpts) || ! isnumeric (minpts) || ! isreal (minpts)
                          || minpts < 1 || fix (minpts) != minpts)
    error ("dbscan: MINPTS must be a positive integer scalar.");
  endif

  ## Detect a precomputed distance matrix among the optional arguments
  precomputed = false;
  for i = 1:2:numel (varargin)
    if (ischar (varargin{i}) && strcmpi (varargin{i}, 'Distance')
        && i < numel (varargin) && ischar (varargin{i+1})
        && strcmpi (varargin{i+1}, 'precomputed'))
      precomputed = true;
    endif
  endfor

  ## Build the epsilon-neighborhood (including the point itself) of each
  ## observation.  For a precomputed distance matrix simply threshold each row;
  ## otherwise delegate the metric handling to rangesearch, which already
  ## includes the query point and applies a "distance <= epsilon" cutoff.
  if (precomputed)
    if (numel (varargin) != 2)
      error (strcat ("dbscan: no other Name-Value arguments are allowed", ...
                     " with the 'precomputed' distance."));
    endif
    if (rows (X) != columns (X))
      error ("dbscan: X must be a square distance matrix for 'precomputed'.");
    endif
    N = rows (X);
    neigh = cell (N, 1);
    for i = 1:N
      neigh{i} = find (X(i,:) <= epsilon)(:);
    endfor
  else
    N = rows (X);
    neigh = rangesearch (X, X, epsilon, varargin{:});
    neigh = cellfun (@(c) c(:), neigh, "UniformOutput", false);
  endif

  ## A point is a core point when its neighborhood (self included) holds at
  ## least minpts observations.
  corepts = cellfun (@numel, neigh) >= minpts;

  ## Canonical DBSCAN scan.  Labels: 0 = unvisited, -1 = noise, k = cluster k.
  labels = zeros (N, 1);
  C = 0;
  for i = 1:N
    if (labels(i) != 0)
      continue;                     # already assigned to a cluster or as noise
    endif
    if (! corepts(i))
      labels(i) = -1;               # tentatively noise (may become a border)
      continue;
    endif
    ## Start a new cluster and expand it over the growing seed queue.
    C += 1;
    labels(i) = C;
    seeds = neigh{i};
    k = 1;
    while (k <= numel (seeds))
      q = seeds(k);
      k += 1;
      if (labels(q) == -1)
        labels(q) = C;              # border point reclaimed from noise
      elseif (labels(q) == 0)
        labels(q) = C;
        if (corepts(q))
          seeds = [seeds; neigh{q}];   # density-reachable: keep expanding
        endif
      endif
    endwhile
  endfor

  idx = labels;

endfunction

%!demo
%! ## Cluster a set of points with two dense blobs and scattered noise.
%! X = [randn(30,2)*0.3 + 2; randn(30,2)*0.3 - 2; 5*(rand(6,2)-0.5)];
%! idx = dbscan (X, 0.6, 4);
%! gscatter (X(:,1), X(:,2), idx);
%! title ("dbscan: clusters (>=0) and noise (-1)");

## Two well-separated clusters and one noise point (euclidean)
%!test
%! X = [0 0; 0 1; 1 0; 1 1; 10 10; 10 11; 11 10; 5 5];
%! [idx, cp] = dbscan (X, 1.5, 3);
%! assert_equal (idx, [1; 1; 1; 1; 2; 2; 2; -1]);
%! assert_equal (cp, logical ([1; 1; 1; 1; 1; 1; 1; 0]));

## Border point reachable from two clusters joins the first one reached
%!test
%! X = [0; 0.3; 0.6; 0.9; 2.5; 2.8; 3.1; 3.4; 1.7];
%! [idx, cp] = dbscan (X, 1, 4);
%! assert_equal (idx, [1; 1; 1; 1; 2; 2; 2; 2; 1]);
%! assert_equal (cp, logical ([1; 1; 1; 1; 1; 1; 1; 1; 0]));

## Cluster order follows the order of the observations, not geometry
%!test
%! X = [2.5; 2.8; 3.1; 3.4; 0; 0.3; 0.6; 0.9; 1.7];
%! idx = dbscan (X, 1, 4);
%! assert_equal (idx, [1; 1; 1; 1; 2; 2; 2; 2; 1]);

## minpts counts the point itself: two neighbors + self meets minpts = 3
%!test
%! X = [0; 1; 2];
%! [idx, cp] = dbscan (X, 1, 2);
%! assert_equal (idx, [1; 1; 1]);
%! assert_equal (cp, logical ([1; 1; 1]));

## minpts = 1 makes every point a core point and yields no noise
%!test
%! X = [0; 0; 5; 5; 10];
%! [idx, cp] = dbscan (X, 0.5, 1);
%! assert_equal (idx, [1; 1; 2; 2; 3]);
%! assert_equal (cp, logical ([1; 1; 1; 1; 1]));

## Precomputed distance matrix matches the euclidean result
%!test
%! X = [0 0; 0 1; 1 0; 1 1; 10 10; 10 11; 11 10; 5 5];
%! D = pdist2 (X, X);
%! [idx, cp] = dbscan (D, 1.5, 3, "Distance", "precomputed");
%! assert_equal (idx, [1; 1; 1; 1; 2; 2; 2; -1]);
%! assert_equal (cp, logical ([1; 1; 1; 1; 1; 1; 1; 0]));

## A non-euclidean metric is passed through to rangesearch
%!test
%! X = [0 0; 0 1; 1 0; 1 1; 10 10; 10 11; 11 10; 5 5];
%! idx = dbscan (X, 2, 3, "Distance", "cityblock");
%! assert_equal (idx, [1; 1; 1; 1; 2; 2; 2; -1]);

## Every point is noise when no neighborhood reaches minpts
%!test
%! X = [0; 10; 20; 30];
%! [idx, cp] = dbscan (X, 1, 2);
%! assert_equal (idx, [-1; -1; -1; -1]);
%! assert_equal (cp, logical ([0; 0; 0; 0]));

## Test input validation
%!error <dbscan: too few input arguments.> dbscan (1)
%!error <dbscan: too few input arguments.> dbscan (1, 1)
%!error <dbscan: X must be a nonempty real numeric matrix.> dbscan ([], 1, 1)
%!error <dbscan: X must be a nonempty real numeric matrix.> dbscan ("a", 1, 1)
%!error <dbscan: X must be a nonempty real numeric matrix.> dbscan (i, 1, 1)
%!error <dbscan: EPSILON must be a nonnegative scalar.> dbscan (ones (3,2), [1 2], 1)
%!error <dbscan: EPSILON must be a nonnegative scalar.> dbscan (ones (3,2), -1, 1)
%!error <dbscan: EPSILON must be a nonnegative scalar.> dbscan (ones (3,2), "a", 1)
%!error <dbscan: MINPTS must be a positive integer scalar.> dbscan (ones (3,2), 1, 0)
%!error <dbscan: MINPTS must be a positive integer scalar.> dbscan (ones (3,2), 1, 1.5)
%!error <dbscan: MINPTS must be a positive integer scalar.> dbscan (ones (3,2), 1, [1 2])
%!error <dbscan: X must be a square distance matrix for 'precomputed'.> ...
%! dbscan (ones (3,2), 1, 1, "Distance", "precomputed")
%!error <dbscan: no other Name-Value arguments are allowed with the 'precomputed' distance.> ...
%! dbscan (ones (3,3), 1, 1, "Distance", "precomputed", "P", 3)
