## Copyright (C) 2016 Nan Zhou <zhnanx@gmail.com>
## Copyright (C) 2021 Stefano Guidoni <ilguido@users.sf.net>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {} silhouette (@var{X}, @var{clust})
## @deftypefnx {statistics} {[@var{si}, @var{h}] =} silhouette (@var{X}, @var{clust})
## @deftypefnx {statistics} {[@var{si}, @var{h}] =} silhouette (@dots{}, @var{Metric}, @var{MetricArg})
##
## Compute the silhouette values of clustered data and show them on a plot.
##
## @var{X} is a n-by-p matrix of n data points in a p-dimensional space.  Each
## datapoint is assigned to a cluster using @var{clust}, a vector of n elements,
## one cluster assignment for each data point.
##
## Each silhouette value of @var{si}, a vector of size n, is a measure of the
## likelihood that a data point is accurately classified to the right cluster.
## Defining "a" as the mean distance between a point and the other points from
## its cluster, and "b" as the mean distance between that point and the points
## from other clusters, the silhouette value of the i-th point is:
##
## @tex
## \def\frac#1#2{{\begingroup#1\endgroup\over#2}}
## $$ S_i = \frac{b_i - a_i}{max(a_1,b_i)} $$
## @end tex
## @ifnottex
## @verbatim
##          bi - ai
## Si =  ------------
##        max(ai,bi)
## @end verbatim
## @end ifnottex
##
## Each element of @var{si} ranges from -1, minimum likelihood of a correct
## classification, to 1, maximum likelihood.
##
## Optional input value @var{Metric} is the metric used to compute the distances
## between data points. Since @code{silhouette} uses @code{pdist} to compute
## these distances, @var{Metric} is similar to the @var{Distance} input argument
## of @code{pdist} and it can be:
## @itemize @bullet
## @item A known distance metric defined as a string: @qcode{euclidean},
## @qcode{squaredeuclidean} (default), @qcode{seuclidean}, @qcode{mahalanobis},
## @qcode{cityblock}, @qcode{minkowski}, @qcode{chebychev}, @qcode{cosine},
## @qcode{correlation}, @qcode{hamming}, @qcode{jaccard}, or @qcode{spearman}.
##
## @item A vector as those created by @code{pdist}. In this case @var{X} does
## nothing.
##
## @item A function handle that is passed to @code{pdist} with @var{MetricArg}
## as optional inputs.
## @end itemize
##
## Optional return value @var{h} is a handle to the silhouette plot.
##
## @strong{Reference}
## Peter J. Rousseeuw, Silhouettes: a Graphical Aid to the Interpretation and
## Validation of Cluster Analysis. 1987. doi:10.1016/0377-0427(87)90125-7
## @end deftypefn
##
## @seealso{dendrogram, evalclusters, kmeans, linkage, pdist}

function [si, h] = silhouette (X, clust, metric = "squaredeuclidean", varargin)
  ## check the input parameters
  if (nargin < 2)
    print_usage ();
  endif

  ## Check for last argument being 'DoNotPlot' to prevent from opening a figure
  ## Undocumented feature to avoid issues with failing tests.  It is only used by
  ## SilhouetteEvaluation class.
  DisplayPlot = true;
  if (numel (varargin) > 0)
    if (ischar (varargin{end}))
      if (strcmp (varargin{end}, "DoNotPlot"))
        DisplayPlot = false;
        varargin{end} = [];
      endif
    endif
  endif

  n = size (clust, 1);

  ## check size
  if (! isempty (X))
    if (size (X, 1) != n)
      error ("First dimension of X <%d> doesn't match that of clust <%d>",...
        size (X, 1), n);
    endif
  endif

  ## check metric
  if (ischar (metric))
    metric = lower (metric);
    switch (metric)
      case "sqeuclidean"
        metric = "squaredeuclidean";
      case {"euclidean", "squaredeuclidean", "seuclidean", "mahalanobis", ...
            "cityblock", "minkowski", "chebychev", "cosine", "correlation", ...
            "hamming", "jaccard", "spearman"}
        ;
      otherwise
        error ("silhouette: invalid metric '%s'", metric);
    endswitch
  elseif (isnumeric (metric) && isvector (metric))
    ## X can be omitted when using this
    distMatrix = squareform (metric);
    if (size (distMatrix, 1) != n)
      error ("First dimension of X <%d> doesn't match that of clust <%d>",...
        size (distMatrix, 1), n);
    endif
  endif

  ## main
  si = zeros(n, 1);
  clusterIDs = unique (clust); # eg [1; 2; 3; 4]
  m = length (clusterIDs);

  ## if only one cluster is defined, the silhouette value is not defined
  if (m == 1)
    si = NaN * ones (n, 1);
    return;
  endif

  ## Precompute cluster membership masks and counts
  clusterMask = false (n, m);
  clusterCount = zeros (m, 1);
  for jjj = 1:m
    clusterMask(:, jjj) = (clust == clusterIDs(jjj));
    clusterCount(jjj) = sum (clusterMask(:, jjj));
  endfor

  ## Map each point to its cluster index (1..m)
  pointClusterIdx = zeros (n, 1);
  for jjj = 1:m
    pointClusterIdx(clusterMask(:, jjj)) = jjj;
  endfor

  ## Use precomputed distance matrix if provided, otherwise compute on-the-fly
  if (exist ('distMatrix', 'var'))
    ## Use full precomputed distance matrix
    for iii = 1:n
      myCluster = pointClusterIdx(iii);
      myMask = clusterMask(:, myCluster);
      myCount = clusterCount(myCluster);

      ## Singleton cluster
      if (myCount == 1)
        si(iii) = 1;
        continue;
      endif

      ## a(i): mean distance to own cluster (excluding self)
      a_i = sum (distMatrix(iii, myMask)) / (myCount - 1);

      ## b(i): minimum mean distance to other clusters
      b_i = Inf;
      for jjj = 1:m
        if (jjj != myCluster)
          b_i = min (b_i, mean (distMatrix(iii, clusterMask(:, jjj))));
        endif
      endfor

      ## s(i)
      si(iii) = (b_i - a_i) / max (a_i, b_i);
    endfor
  else
    ## Compute distances row-by-row to avoid O(n^2) memory
    for iii = 1:n
      ## Compute distances from point iii to all other points
      diffs = X - X(iii, :);
      switch (metric)
        case "squaredeuclidean"
          dists = sum (diffs .^ 2, 2);
        case "euclidean"
          dists = sqrt (sum (diffs .^ 2, 2));
        case "cityblock"
          dists = sum (abs (diffs), 2);
        case "chebychev"
          dists = max (abs (diffs), [], 2);
        otherwise
          ## Fall back to pdist2-style computation for other metrics
          dists = pdist2 (X(iii, :), X, metric, varargin{:})';
      endswitch

      myCluster = pointClusterIdx(iii);
      myMask = clusterMask(:, myCluster);
      myCount = clusterCount(myCluster);

      ## Singleton cluster
      if (myCount == 1)
        si(iii) = 1;
        continue;
      endif

      ## a(i): mean distance to own cluster (excluding self, dist=0)
      a_i = sum (dists(myMask)) / (myCount - 1);

      ## b(i): minimum mean distance to other clusters
      b_i = Inf;
      for jjj = 1:m
        if (jjj != myCluster)
          b_i = min (b_i, mean (dists(clusterMask(:, jjj))));
        endif
      endfor

      ## s(i)
      si(iii) = (b_i - a_i) / max (a_i, b_i);
    endfor
  endif

  ## plot
  ## a poor man silhouette graph
  if (DisplayPlot)
    vBarsc = zeros (m, 1);
    vPadding = [0; 0; 0; 0];
    Bars = vPadding;

    for i = 1 : m
      vBar = si(find (clust == clusterIDs(i)));
      vBarsc(i) = length (Bars) + (length (vBar) / 2);
      Bars = [Bars; (sort (vBar, "descend")); vPadding];
    endfor

    figure();
    h = barh (Bars, "hist", "facecolor", [0 0.4471 0.7412]);

    xlabel ("Silhouette Value");
    ylabel ("Cluster");
    set (gca, "ytick", vBarsc, "yticklabel", clusterIDs);
    ylim ([0 (length (Bars))]);
    axis ("ij");
  endif
endfunction


%!demo
%! load fisheriris;
%! X = meas(:,3:4);
%! cidcs = kmeans (X, 3, "Replicates", 5);
%! silhouette (X, cidcs);
%! y_labels(cidcs([1 51 101])) = unique (species);
%! set (gca, "yticklabel", y_labels);
%! title ("Fisher's iris data");

## Test input validation
%!error silhouette ();
%!error silhouette ([1 2; 1 1]);
%!error <X .* doesn't match .* clust> silhouette ([1 2; 1 1], [1 2 3]');
%!error <invalid metric> silhouette ([1 2; 1 1], [1 2]', "xxx");
