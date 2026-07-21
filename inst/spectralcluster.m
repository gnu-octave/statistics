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
## @deftypefn  {statistics} {@var{idx} =} spectralcluster (@var{X}, @var{k})
## @deftypefnx {statistics} {@var{idx} =} spectralcluster (@var{S}, @var{k}, @qcode{'Distance'}, @qcode{'precomputed'})
## @deftypefnx {statistics} {[@var{idx}, @var{V}] =} spectralcluster (@dots{})
## @deftypefnx {statistics} {[@var{idx}, @var{V}, @var{D}] =} spectralcluster (@dots{})
## @deftypefnx {statistics} {[@dots{}] =} spectralcluster (@dots{}, @var{name}, @var{value})
##
## Partition observations into @var{k} clusters using spectral clustering.
##
## @code{@var{idx} = spectralcluster (@var{X}, @var{k})} partitions the
## @math{N*P} numeric matrix @var{X} into @var{k} clusters and returns the
## @math{N*1} vector @var{idx} of cluster indices.  Rows of @var{X} correspond
## to observations and columns to features.  Spectral clustering builds a
## similarity graph over the observations, embeds them with the eigenvectors of
## the graph Laplacian, and clusters that embedding, which lets it recover
## clusters that are not linearly separable in the original space.
##
## @code{[@var{idx}, @var{V}, @var{D}] = spectralcluster (@dots{})} also returns
## the @math{N*k} matrix @var{V} whose columns are the eigenvectors associated
## with the @var{k} smallest eigenvalues of the Laplacian, and the @math{k*1}
## vector @var{D} of those eigenvalues.  The signs of the eigenvectors, and the
## basis within a repeated eigenvalue, are arbitrary.
##
## Additional parameters can be specified by @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.22 0.76
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Distance'} @tab the distance metric used to build the
## similarity graph, one of @qcode{'euclidean'} (default), @qcode{'seuclidean'},
## @qcode{'mahalanobis'}, @qcode{'cityblock'}, @qcode{'minkowski'},
## @qcode{'chebychev'}, @qcode{'cosine'}, @qcode{'correlation'},
## @qcode{'hamming'}, @qcode{'jaccard'}, @qcode{'spearman'}, or a function
## handle accepted by @code{pdist2}, or the string @qcode{'precomputed'} to
## interpret the first input as an @math{N*N} similarity matrix.
##
## @item @qcode{'SimilarityGraph'} @tab @qcode{'knn'} (default) to connect each
## observation to its nearest neighbors, or @qcode{'epsilon'} to connect
## observations that are within a fixed radius.
##
## @item @qcode{'NumNeighbors'} @tab the number of nearest neighbors for the
## @qcode{'knn'} graph, a positive integer.  The default is
## @code{ceil (log (@var{N}))}.
##
## @item @qcode{'KNNGraphType'} @tab @qcode{'complete'} (default) to connect
## @math{i} and @math{j} when either is a nearest neighbor of the other, or
## @qcode{'mutual'} to connect them only when each is a nearest neighbor of the
## other.
##
## @item @qcode{'Radius'} @tab the radius for the @qcode{'epsilon'} graph, a
## nonnegative scalar.  Required when @qcode{'SimilarityGraph'} is
## @qcode{'epsilon'}.
##
## @item @qcode{'KernelScale'} @tab the positive scale factor @var{sigma} in
## the Gaussian similarity kernel @code{exp (-(dist / sigma)^2)} applied to the
## graph edges.  The default is 1.
##
## @item @qcode{'LaplacianNormalization'} @tab @qcode{'randomwalk'} (default),
## @qcode{'symmetric'}, or @qcode{'none'}, selecting how the graph Laplacian is
## normalized before the eigendecomposition.
##
## @item @qcode{'ClusterMethod'} @tab @qcode{'kmeans'} (default) or
## @qcode{'kmedoids'} to cluster the eigenvector embedding.
##
## @item @qcode{'P'} @tab the Minkowski exponent (default 2), used only with the
## @qcode{'minkowski'} distance.
##
## @item @qcode{'Cov'} @tab the covariance matrix used only with the
## @qcode{'mahalanobis'} distance.
##
## @item @qcode{'Scale'} @tab the scaling vector used only with the
## @qcode{'seuclidean'} distance.
## @end multitable
##
## @seealso{kmeans, kmedoids, dbscan, linkage, pdist2}
## @end deftypefn

function [idx, V, D] = spectralcluster (X, k, varargin)

  ## Check number of input arguments
  if (nargin < 2)
    error ("spectralcluster: too few input arguments.");
  endif

  ## Validate X
  if (! isnumeric (X) || ! isreal (X) || ndims (X) != 2 || isempty (X))
    error ("spectralcluster: X must be a nonempty real numeric matrix.");
  endif
  N = rows (X);

  ## Validate k
  if (! isscalar (k) || ! isnumeric (k) || ! isreal (k) || k < 1
                     || fix (k) != k)
    error ("spectralcluster: K must be a positive integer scalar.");
  endif
  if (k > N)
    error ("spectralcluster: K cannot exceed the number of observations in X.");
  endif

  ## Defaults
  distance     = "euclidean";
  simgraph     = "knn";
  numneighbors = [];
  knntype      = "complete";
  radius       = [];
  kernelscale  = 1;
  laplnorm     = "randomwalk";
  clustermethod = "kmeans";
  Pexp = 2;
  Cov  = [];
  Scl  = [];

  ## Parse Name-Value pairs
  if (mod (numel (varargin), 2) != 0)
    error ("spectralcluster: each NAME must be followed by a VALUE.");
  endif
  while (numel (varargin) > 0)
    name = varargin{1};
    val  = varargin{2};
    if (! ischar (name))
      error ("spectralcluster: optional argument names must be strings.");
    endif
    switch (tolower (name))
      case "distance"
        distance = val;
      case "similaritygraph"
        simgraph = tolower (val);
      case "numneighbors"
        numneighbors = val;
      case "knngraphtype"
        knntype = tolower (val);
      case "radius"
        radius = val;
      case "kernelscale"
        kernelscale = val;
      case "laplaciannormalization"
        laplnorm = tolower (val);
      case "clustermethod"
        clustermethod = tolower (val);
      case "p"
        Pexp = val;
      case "cov"
        Cov = val;
      case "scale"
        Scl = val;
      otherwise
        error ("spectralcluster: unknown parameter name '%s'.", name);
    endswitch
    varargin(1:2) = [];
  endwhile

  ## Validate option choices
  if (! any (strcmp (simgraph, {"knn", "epsilon"})))
    error ("spectralcluster: SIMILARITYGRAPH must be 'knn' or 'epsilon'.");
  endif
  if (! any (strcmp (knntype, {"complete", "mutual"})))
    error ("spectralcluster: KNNGRAPHTYPE must be 'complete' or 'mutual'.");
  endif
  if (! any (strcmp (laplnorm, {"randomwalk", "symmetric", "none"})))
    error (strcat ("spectralcluster: LAPLACIANNORMALIZATION must be", ...
                   " 'randomwalk', 'symmetric', or 'none'."));
  endif
  if (! any (strcmp (clustermethod, {"kmeans", "kmedoids"})))
    error ("spectralcluster: CLUSTERMETHOD must be 'kmeans' or 'kmedoids'.");
  endif
  if (! isscalar (kernelscale) || ! isnumeric (kernelscale)
                               || kernelscale <= 0)
    error ("spectralcluster: KERNELSCALE must be a positive scalar.");
  endif

  ## Build the N-by-N similarity matrix S.
  if (ischar (distance) && strcmpi (distance, "precomputed"))
    if (rows (X) != columns (X))
      error (strcat ("spectralcluster: X must be a square similarity", ...
                     " matrix for the 'precomputed' distance."));
    endif
    S = X;
  else
    ## Pairwise distances in the requested metric.
    metrics = {"euclidean", "seuclidean", "mahalanobis", "cityblock", ...
               "minkowski", "chebychev", "cosine", "correlation", ...
               "hamming", "jaccard", "spearman"};
    if (ischar (distance))
      dname = tolower (distance);
      if (! any (strcmp (dname, metrics)))
        error ("spectralcluster: unsupported distance metric '%s'.", distance);
      endif
      switch (dname)
        case "minkowski"
          Dmat = pdist2 (X, X, dname, Pexp);
        case "mahalanobis"
          if (isempty (Cov))
            Dmat = pdist2 (X, X, dname);
          else
            Dmat = pdist2 (X, X, dname, Cov);
          endif
        case "seuclidean"
          if (isempty (Scl))
            Dmat = pdist2 (X, X, dname);
          else
            Dmat = pdist2 (X, X, dname, Scl);
          endif
        otherwise
          Dmat = pdist2 (X, X, dname);
      endswitch
    elseif (is_function_handle (distance))
      Dmat = pdist2 (X, X, distance);
    else
      error (strcat ("spectralcluster: DISTANCE must be a metric name or", ...
                     " a function handle."));
    endif

    ## Adjacency of the similarity graph.
    A = false (N);
    if (strcmp (simgraph, "knn"))
      if (isempty (numneighbors))
        numneighbors = max (1, ceil (log (N)));
      endif
      if (! isscalar (numneighbors) || ! isnumeric (numneighbors)
          || numneighbors < 1 || fix (numneighbors) != numneighbors)
        error (strcat ("spectralcluster: NUMNEIGHBORS must be a positive", ...
                       " integer scalar."));
      endif
      nn = min (numneighbors, N - 1);
      for i = 1:N
        d = Dmat(i,:);
        d(i) = Inf;
        [~, ord] = sort (d);
        A(i, ord(1:nn)) = true;
      endfor
      if (strcmp (knntype, "complete"))
        A = A | A';
      else
        A = A & A';
      endif
    else
      if (isempty (radius))
        error (strcat ("spectralcluster: RADIUS is required for the", ...
                       " 'epsilon' similarity graph."));
      endif
      if (! isscalar (radius) || ! isnumeric (radius) || radius < 0)
        error ("spectralcluster: RADIUS must be a nonnegative scalar.");
      endif
      A = Dmat <= radius;
      A(1:(N + 1):end) = false;         # drop self-edges on the diagonal
    endif

    ## Gaussian kernel on the graph edges; zero elsewhere.
    S = zeros (N);
    S(A) = exp (-(Dmat(A) / kernelscale) .^ 2);
  endif

  ## Graph Laplacian and its k smallest eigenpairs.
  deg = sum (S, 2);
  Dg  = diag (deg);
  if (! strcmp (laplnorm, "none") && any (deg <= 0))
    error (strcat ("spectralcluster: the similarity graph has isolated", ...
                   " points; increase NumNeighbors or Radius."));
  endif

  L = Dg - S;
  switch (laplnorm)
    case "none"
      [Vec, Val] = eig (L);
    case "randomwalk"
      [Vec, Val] = eig (L, Dg);
    case "symmetric"
      Dh = diag (1 ./ sqrt (deg));
      Ls = Dh * L * Dh;
      Ls = (Ls + Ls') / 2;              # enforce exact symmetry for real eig
      [Vec, Val] = eig (Ls);
  endswitch

  ev = real (diag (Val));
  [ev, order] = sort (ev);
  D = ev(1:k);
  V = real (Vec(:, order(1:k)));

  if (strcmp (laplnorm, "symmetric"))
    rn = sqrt (sum (V .^ 2, 2));         # row-normalize the embedding
    rn(rn == 0) = 1;
    V = V ./ rn;
  endif

  ## Cluster the eigenvector embedding.
  if (strcmp (clustermethod, "kmeans"))
    idx = kmeans (V, k, "Replicates", 5, "EmptyAction", "singleton");
  else
    idx = kmedoids (V, k, "Replicates", 5);
  endif

endfunction

%!demo
%! ## Two concentric rings are not separable by kmeans but are by spectral
%! ## clustering.
%! t = linspace (0, 2*pi, 100)';
%! Xin  = [cos(t), sin(t)] + randn (100, 2) * 0.05;
%! Xout = 4 * [cos(t), sin(t)] + randn (100, 2) * 0.05;
%! X = [Xin; Xout];
%! idx = spectralcluster (X, 2, "NumNeighbors", 10);
%! gscatter (X(:,1), X(:,2), idx);
%! axis equal;
%! title ("spectralcluster: two concentric rings");

## Exact anchor: precomputed similarity, random-walk Laplacian (MATLAB R2024a)
%!test
%! S = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0.1 0 0; 0 0 0.1 0 1 1; ...
%!      0 0 0 1 0 1; 0 0 0 1 1 0];
%! [idx, V, D] = spectralcluster (S, 2, "Distance", "precomputed");
%! assert (abs (D(1)) < 1e-9);
%! assert_equal (D(2), 0.0314065796348158, 1e-12);
%! assert_equal (size (V), [6 2]);
%! assert (idx(1) == idx(2) && idx(2) == idx(3));
%! assert (idx(4) == idx(5) && idx(5) == idx(6));
%! assert (idx(1) != idx(4));

## Symmetric normalization has the same spectrum as random-walk
%!test
%! S = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0.1 0 0; 0 0 0.1 0 1 1; ...
%!      0 0 0 1 0 1; 0 0 0 1 1 0];
%! [~, ~, D] = spectralcluster (S, 2, "Distance", "precomputed", ...
%!                              "LaplacianNormalization", "symmetric");
%! assert_equal (D(2), 0.0314065796348156, 1e-12);

## Unnormalized Laplacian ('none') spectrum
%!test
%! S = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0.1 0 0; 0 0 0.1 0 1 1; ...
%!      0 0 0 1 0 1; 0 0 0 1 1 0];
%! [~, ~, D] = spectralcluster (S, 2, "Distance", "precomputed", ...
%!                              "LaplacianNormalization", "none");
%! assert_equal (D(2), 0.0637708504262773, 1e-12);

## Precomputed path graph splits at the weak middle edge; eigenvalue is 2/7
%!test
%! Sc = [0 0.5 0 0; 0.5 0 0.2 0; 0 0.2 0 0.5; 0 0 0.5 0];
%! [idx, ~, D] = spectralcluster (Sc, 2, "Distance", "precomputed");
%! assert_equal (D(2), 2/7, 1e-12);
%! assert (idx(1) == idx(2));
%! assert (idx(3) == idx(4));
%! assert (idx(1) != idx(3));

## Full knn pipeline (graph + Gaussian kernel + Laplacian) matches MATLAB
%!test
%! X = [0 0; 0.1 0; 0 0.1; 5 0; 5.1 0; 5 0.1];
%! [idx, ~, D] = spectralcluster (X, 2, "NumNeighbors", 3, "KernelScale", 2);
%! assert_equal (D(2), 0.00358001836237202, 1e-11);
%! assert (idx(1) == idx(2) && idx(2) == idx(3));
%! assert (idx(4) == idx(5) && idx(5) == idx(6));
%! assert (idx(1) != idx(4));

## Default NumNeighbors is ceil (log (N)): for N = 8 that is 3
%!test
%! X = (1:8)';
%! [~, ~, Ddef] = spectralcluster (X, 2);
%! [~, ~, D3]   = spectralcluster (X, 2, "NumNeighbors", 3);
%! assert_equal (Ddef(2), D3(2), 1e-12);
%! assert_equal (Ddef(2), 0.113482181378817, 1e-11);

## Default path clusters two well-separated blobs
%!test
%! X = [0 0; 0.2 0; 0 0.2; 0.2 0.2; 10 10; 10.2 10; 10 10.2; 10.2 10.2];
%! idx = spectralcluster (X, 2);
%! assert (numel (unique (idx(1:4))) == 1);
%! assert (numel (unique (idx(5:8))) == 1);
%! assert (idx(1) != idx(5));

## kmedoids embedding-clustering and the epsilon graph both work
%!test
%! X = [0 0; 0.2 0; 0 0.2; 0.2 0.2; 10 10; 10.2 10; 10 10.2; 10.2 10.2];
%! idx = spectralcluster (X, 2, "ClusterMethod", "kmedoids");
%! assert (numel (unique (idx(1:4))) == 1);
%! assert (numel (unique (idx(5:8))) == 1);
%!test
%! X = [0 0; 0.2 0; 0 0.2; 0.2 0.2; 10 10; 10.2 10; 10 10.2; 10.2 10.2];
%! [idx, ~, D] = spectralcluster (X, 2, "SimilarityGraph", "epsilon", ...
%!                                "Radius", 1);
%! assert_equal (numel (idx), 8);
%! assert (idx(1) != idx(5));

## Three clusters: mutual knn graph and output shapes
%!test
%! X = [0 0; 0.2 0; 0 0.2; 10 10; 10.2 10; 10 10.2; 0 10; 0.2 10; 0 10.2];
%! [idx, V, D] = spectralcluster (X, 3, "NumNeighbors", 2, ...
%!                                "KNNGraphType", "mutual");
%! assert_equal (size (V), [9 3]);
%! assert_equal (size (D), [3 1]);
%! assert (numel (unique (idx(1:3))) == 1);
%! assert (numel (unique (idx(4:6))) == 1);
%! assert (numel (unique (idx(7:9))) == 1);

## Test input validation
%!error <spectralcluster: too few input arguments.> spectralcluster (1)
%!error <spectralcluster: X must be a nonempty real numeric matrix.> ...
%! spectralcluster ([], 2)
%!error <spectralcluster: K must be a positive integer scalar.> ...
%! spectralcluster (ones (4,2), 0)
%!error <spectralcluster: K cannot exceed the number of observations in X.> ...
%! spectralcluster (ones (3,2), 4)
%!error <spectralcluster: each NAME must be followed by a VALUE.> ...
%! spectralcluster (ones (4,2), 2, "Distance")
%!error <spectralcluster: unknown parameter name 'foo'.> ...
%! spectralcluster (ones (4,2), 2, "foo", "bar")
%!error <spectralcluster: unsupported distance metric 'taxicab'.> ...
%! spectralcluster (ones (4,2), 2, "Distance", "taxicab")
%!error <spectralcluster: SIMILARITYGRAPH must be 'knn' or 'epsilon'.> ...
%! spectralcluster (ones (4,2), 2, "SimilarityGraph", "tree")
%!error <spectralcluster: KNNGRAPHTYPE must be 'complete' or 'mutual'.> ...
%! spectralcluster (ones (4,2), 2, "KNNGraphType", "partial")
%!error <spectralcluster: CLUSTERMETHOD must be 'kmeans' or 'kmedoids'.> ...
%! spectralcluster (ones (4,2), 2, "ClusterMethod", "dbscan")
%!error <spectralcluster: KERNELSCALE must be a positive scalar.> ...
%! spectralcluster (ones (4,2), 2, "KernelScale", 0)
%!error <spectralcluster: RADIUS is required for the 'epsilon' similarity graph.> ...
%! spectralcluster (ones (4,2), 2, "SimilarityGraph", "epsilon")
%!error <spectralcluster: X must be a square similarity matrix for the 'precomputed' distance.> ...
%! spectralcluster (ones (4,2), 2, "Distance", "precomputed")
%!error <spectralcluster: the similarity graph has isolated points; increase NumNeighbors or Radius.> ...
%! spectralcluster ([0 0; 10 10], 2, "SimilarityGraph", "epsilon", "Radius", 1)
