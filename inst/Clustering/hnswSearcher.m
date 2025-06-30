## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
## Copyright (C) 2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef hnswSearcher
## -*- texinfo -*-
## @deftp {Class} hnswSearcher
##
## Hierarchical Navigable Small World (HNSW) nearest neighbor searcher class.
##
## The @code{hnswSearcher} class implements the HNSW algorithm for efficient
## nearest neighbor queries.  It stores training data and supports various
## distance metrics for performing searches.  The HNSW algorithm builds a
## multilayer graph structure that enables fast approximate nearest neighbor
## searches by navigating through the graph.  It facilitates a nearest neighbor
## search using @code{knnsearch}.
##
## You can either use the @code{hnswSearcher} class constructor or the
## @code{createns} function to create an @qcode{hnswSearcher} object.
##
## @seealso{createns, ExhaustiveSearcher, KDTreeSearcher, knnsearch}
## @end deftp

  properties (SetAccess = private, Hidden)
    ## -*- texinfo -*-
    ## @deftp {Property} HNSWGraph
    ##
    ## The HNSW graph structure built from the training data.  This property is
    ## private and cannot be modified after object creation.
    ##
    ## @end deftp
    HNSWGraph
  endproperties

  properties (SetAccess = private)
    ## -*- texinfo -*-
    ## @deftp {Property} Distance
    ##
    ## Distance metric used for searches, specified as a character vector (e.g.,
    ## @qcode{"euclidean"}, @qcode{"minkowski"}, @qcode{"cityblock"}). Default
    ## is @qcode{"euclidean"}. Supported metrics align with those in
    ## @code{pdist2}.
    ##
    ## @end deftp
    Distance = 'euclidean'

    ## -*- texinfo -*-
    ## @deftp {Property} DistParameter
    ##
    ## Parameter for the distance metric, with type and value depending on
    ## @qcode{Distance}:
    ##
    ## @itemize
    ## @item For @qcode{"minkowski"}, a positive scalar exponent (default 2).
    ## @item For @qcode{"seuclidean"}, a nonnegative vector of scaling factors
    ## matching the number of columns in @qcode{X} (default is standard
    ## deviation of @qcode{X}).
    ## @item For @qcode{"mahalanobis"}, a positive definite covariance matrix
    ## matching the dimensions of @qcode{X} (default is @code{cov (@var{X})}).
    ## @item Empty for other metrics.
    ## @end itemize
    ##
    ## @end deftp
    DistParameter = []

    ## -*- texinfo -*-
    ## @deftp {Property} MaxNumLinksPerNode
    ##
    ## Maximum number of neighbors per node in the HNSW graph. Affects graph
    ## connectivity and search accuracy. Default is 16.
    ##
    ## @end deftp
    MaxNumLinksPerNode = 16

    ## -*- texinfo -*-
    ## @deftp {Property} TrainSetSize
    ##
    ## Size of the dynamic candidate list during graph construction. Higher
    ## values improve accuracy at the cost of construction time. Default is 200.
    ##
    ## @end deftp
    TrainSetSize = 200

    ## -*- texinfo -*-
    ## @deftp {Property} X
    ##
    ## Training data, specified as an @math{NxP} numeric matrix where each row
    ## is an observation and each column is a feature.  This property is private
    ## and cannot be modified after object creation.
    ##
    ## @end deftp
    X = []
  endproperties

  methods (Hidden)

    ## Custom display
    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    ## Custom display
    function disp (this)
      if (isscalar (this))
        fprintf ("\n  hnswSearcher with properties:\n\n");
        fprintf ("%+25s: %d\n", 'MaxNumLinksPerNode', this.MaxNumLinksPerNode);
        fprintf ("%+25s: %d\n", 'TrainSetSize', this.TrainSetSize);
        fprintf ("%+25s: '%s'\n", 'Distance', this.Distance);
        if (! isempty (this.DistParameter))
          if (isscalar (this.DistParameter))
            fprintf ("%+25s: %g\n", 'DistParameter', this.DistParameter);
          elseif (isvector (this.DistParameter))
            fprintf ("%+25s: %s\n", 'DistParameter', ...
                     mat2str (this.DistParameter));
          else
            fprintf ("%+25s: [%dx%d %s]\n", 'DistParameter', ...
                     size (this.DistParameter), class (this.DistParameter));
          endif
        else
          fprintf ("%+25s: []\n", 'DistParameter');
        endif
        fprintf ("%+25s: [%dx%d %s]\n", 'X', size (this.X), class (this.X));
      endif
    endfunction

    ## Class specific subscripted reference
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error ("hnswSearcher.subsref: () indexing not supported.");
        case '{}'
          error ("hnswSearcher.subsref: {} indexing not supported.");
        case '.'
          if (! ischar (s.subs))
            error (strcat ("hnswSearcher.subsref: property", ...
                           " name must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("hnswSearcher.subsref: unrecognized", ...
                           " property: '%s'."), s.subs);
          end_try_catch
      endswitch
      ## Chained references
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

    ## Class specific subscripted assignment
    function this = subsasgn (this, s, val)
      switch s.type
        case '()'
          error ("hnswSearcher.subsasgn: () indexing not supported.");
        case '{}'
          error ("hnswSearcher.subsasgn: {} indexing not supported.");
        case '.'
          if (! ischar (s.subs))
            error (strcat ("hnswSearcher.subsasgn: property", ...
                           " name must be a character vector."));
          endif
          switch (s.subs)
            case 'X'
              error (strcat ("hnswSearcher.subsasgn: 'X' is", ...
                             " read-only and cannot be modified."));
            case 'HNSWGraph'
              error (strcat ("hnswSearcher.subsasgn: 'HNSWGraph' is", ...
                             " read-only and cannot be modified."));
            case 'Distance'
              error (strcat ("hnswSearcher.subsasgn: 'Distance' is", ...
                             " read-only and cannot be modified."));
            case 'DistParameter'
              error (strcat ("hnswSearcher.subsasgn: 'DistParameter'", ...
                             " is read-only and cannot be modified."));
            case 'MaxNumLinksPerNode'
              error (strcat ("hnswSearcher.subsasgn: 'MaxNumLinksPerNode'", ...
                             " is read-only and cannot be modified."));
            case 'TrainSetSize'
              error (strcat ("hnswSearcher.subsasgn: 'TrainSetSize'", ...
                             " is read-only and cannot be modified."));
            otherwise
              error (strcat ("hnswSearcher.subsasgn:", ...
                             " unrecognized property: '%s'."), s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods

    ## -*- texinfo -*-
    ## @deftypefn  {hnswSearcher} {@var{obj} =} hnswSearcher (@var{X})
    ## @deftypefnx {hnswSearcher} {@var{obj} =} hnswSearcher (@var{X}, @var{name}, @var{value})
    ##
    ## Create an @qcode{hnswSearcher} object for approximate nearest neighbor
    ## searches.
    ##
    ## @code{@var{obj} = hnswSearcher (@var{X})} constructs an
    ## @qcode{hnswSearcher} object with training data @var{X} using the
    ## default @qcode{"euclidean"} distance metric. @var{X} must be an
    ## @math{NxP} numeric matrix, where rows represent observations and columns
    ## represent features.
    ##
    ## @code{@var{obj} = hnswSearcher (@var{X}, @var{name}, @var{value})}
    ## allows customization through name-value pairs:
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"Distance"} @tab @tab Distance metric, specified as a
    ## character vector (e.g., @qcode{"euclidean"}, @qcode{"minkowski"},
    ## @qcode{"cityblock"}). Default is @qcode{"euclidean"}. See @code{pdist2}
    ## for supported metrics.
    ##
    ## @item @qcode{"P"} @tab @tab Minkowski distance exponent, a positive
    ## scalar. Valid only when @qcode{"Distance"} is @qcode{"minkowski"}.
    ## Default is 2.
    ##
    ## @item @qcode{"Scale"} @tab @tab Nonnegative vector of scaling factors
    ## matching the number of columns in @var{X}. Valid only when
    ## @qcode{"Distance"} is @qcode{"seuclidean"}. Default is @code{std (X)}.
    ##
    ## @item @qcode{"Cov"} @tab @tab Positive definite covariance matrix
    ## matching the number of columns in @var{X}. Valid only when
    ## @qcode{"Distance"} is @qcode{"mahalanobis"}. Default is @code{cov (X)}.
    ##
    ## @item @qcode{"MaxNumLinksPerNode"} @tab @tab Maximum number of neighbors
    ## per node in the HNSW graph, a positive integer. Default is 16.
    ##
    ## @item @qcode{"TrainSetSize"} @tab @tab Size of the dynamic candidate
    ## list during graph construction, a positive integer. Default is 200.
    ## @end multitable
    ##
    ## @seealso{hnswSearcher, knnsearch, createns, pdist2}
    ## @end deftypefn
    function obj = hnswSearcher (X, varargin)
      ## Initial input validation
      if (nargin < 1)
        error ("hnswSearcher: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error ("hnswSearcher: Name-Value arguments must be in pairs.");
      endif

      ## Validate X
      if (isempty (X))
        error ("hnswSearcher: X cannot be empty.");
      endif
      if (! (isnumeric (X) && ismatrix (X) && all (isfinite (X(:)))))
        error ("hnswSearcher: X must be a finite numeric matrix.");
      endif

      obj.X = X;
      N = size (X, 1);

      ## Default values
      Distance = "euclidean";
      P = [];
      S = [];
      C = [];
      MaxNumLinksPerNode = min (16, N);
      TrainSetSize = min (200, N);
      TrainSetSize = max (TrainSetSize, MaxNumLinksPerNode);

      ## Parse optional parameters
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case "distance"
            Distance = varargin{2};
          case "p"
            P = varargin{2};
          case "scale"
            S = varargin{2};
          case "cov"
            C = varargin{2};
          case "maxnumlinkspernode"
            MaxNumLinksPerNode = varargin{2};
          case "trainsetsize"
            TrainSetSize = varargin{2};
          otherwise
            error ("hnswSearcher: invalid parameter name: '%s'.", varargin{1});
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Validate Distance
      valid_metrics = {'euclidean', 'minkowski', 'seuclidean', ...
                       'mahalanobis', 'cityblock', 'manhattan', ...
                       'chebychev', 'cosine', 'correlation', ...
                       'spearman', 'hamming', 'jaccard'};
      if (ischar (Distance))
        if (! any (strcmpi (valid_metrics, Distance)))
          error ("hnswSearcher: unsupported distance metric '%s'.", Distance);
        endif
        obj.Distance = Distance;
      else
        error ("hnswSearcher: 'Distance' must be a string.");
      endif

      ## Set DistParameter
      if (strcmpi (obj.Distance, "minkowski"))
        if (isempty (P))
          obj.DistParameter = 2;
        else
          if (! (isscalar (P) && isnumeric (P) && P > 0 && isfinite (P)))
            error ("hnswSearcher: 'P' must be a positive finite scalar.");
          endif
          obj.DistParameter = P;
        endif
      elseif (strcmpi (obj.Distance, "seuclidean"))
        if (isempty (S))
          obj.DistParameter = std (X, [], 1);
        else
          if (! (isvector (S) && isnumeric (S) && all (S >= 0)
                                               && all (isfinite (S))
                                               && length (S) == columns (X)))
            error (strcat ("hnswSearcher: 'Scale' must be a", ...
                           " nonnegative vector matching X columns."));
          endif
          obj.DistParameter = S;
        endif
      elseif (strcmpi (obj.Distance, "mahalanobis"))
        if (isempty (C))
          obj.DistParameter = cov (X);
        else
          if (! (ismatrix (C) && isnumeric (C) && all (isfinite (C)(:))
                                               && rows (C) == columns (C)
                                               && rows (C) == columns (X)))
            error (strcat ("hnswSearcher: 'Cov' must be a square", ...
                           " matrix matching X columns."));
          endif
          if (! issymmetric (C))
            error (strcat ("hnswSearcher: 'Cov' must be symmetric", ...
                           " for mahalanobis."));
          endif
          [~, p] = chol (C);
          if (p != 0)
            error (strcat ("hnswSearcher: 'Cov' must be positive", ...
                           " definite for mahalanobis."));
          endif
          obj.DistParameter = C;
        endif
      else
        obj.DistParameter = [];
      endif

      ## Validate MaxNumLinksPerNode and TrainSetSize
      if (! (isscalar (MaxNumLinksPerNode) &&
             isnumeric (MaxNumLinksPerNode) &&
             MaxNumLinksPerNode > 0 &&
             MaxNumLinksPerNode == fix (MaxNumLinksPerNode)))
        error (strcat ("hnswSearcher: 'MaxNumLinksPerNode'", ...
                       " must be a positive integer."));
      endif
      if (! (isscalar (TrainSetSize) &&
             isnumeric (TrainSetSize) &&
             TrainSetSize > 0 &&
             TrainSetSize == fix (TrainSetSize)))
        error ("hnswSearcher: 'TrainSetSize' must be a positive integer.");
      endif
      if (TrainSetSize > N)
        error (strcat ("hnswSearcher: 'TrainSetSize' cannot", ...
                       " exceed the number of rows in X."));
      endif
      if (MaxNumLinksPerNode > TrainSetSize)
        error (strcat ("hnswSearcher: 'MaxNumLinksPerNode'", ...
                       " cannot exceed 'TrainSetSize'."));
      endif
      obj.MaxNumLinksPerNode = MaxNumLinksPerNode;
      obj.TrainSetSize = TrainSetSize;

      ## Build HNSW graph
      obj.HNSWGraph = build_hnsw (X, obj.Distance, obj.DistParameter, ...
                                  MaxNumLinksPerNode, TrainSetSize);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {hnswSearcher} {[@var{idx}, @var{D}] =} knnsearch (@var{obj}, @var{Y})
    ## @deftypefnx {hnswSearcher} {[@var{idx}, @var{D}] =} knnsearch (@var{obj}, @var{Y}, @var{name}, @var{value})
    ##
    ## Find the nearest neighbors in the training data to query points.
    ##
    ## @code{[@var{idx}, @var{D}] = knnsearch (@var{obj}, @var{Y})} returns the
    ## indices @var{idx} and distances @var{D} of the nearest neighbor in
    ## @var{obj.X} to each point in @var{Y}, using the distance metric specified
    ## in @var{obj.Distance}.
    ##
    ## @itemize
    ## @item @var{obj} is an @qcode{hnswSearcher} object.
    ## @item @var{Y} is an @math{MxP} numeric matrix of query points, where
    ## @math{P} must match the number of columns in @var{obj.X}.
    ## @item @var{idx} contains the indices of the nearest neighbors in
    ## @var{obj.X}.
    ## @item @var{D} contains the corresponding distances.
    ## @end itemize
    ##
    ## @code{[@var{idx}, @var{D}] = knnsearch (@var{obj}, @var{Y}, @var{name}, @var{value})}
    ## allows additional options via name-value pairs:
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"K"} @tab @tab A positive integer specifying the number of
    ## nearest neighbors to find. Default is 1.
    ##
    ## @item @qcode{"SearchSetSize"} @tab @tab A positive integer specifying the
    ## size of the candidate list of nearest neighbors for a single query point
    ## during the search process.  Default is @qcode{max (10, @var{C})}, where
    ## @var{C} is the number of columns in @var{obj.X}.  @qcode{"SearchSetSize"}
    ## must be at least @var{C} and no more than the number of rows in training
    ## data @var{obj.X}.
    ## @end multitable
    ##
    ## @seealso{hnswSearcher, pdist2}
    ## @end deftypefn
    function [idx, D] = knnsearch (obj, Y, varargin)
      ## Initial input validation
      if (nargin < 2)
        error ("hnswSearcher.knnsearch: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("hnswSearcher.knnsearch:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      ## Get training data size
      [N, C] = size (obj.X);

      ## Validate Y
      if (isempty (Y))
        error ("hnswSearcher.knnsearch: Y cannot be empty.");
      endif
      if (! (isnumeric (Y) && ismatrix (Y) && all (isfinite (Y)(:))))
        error ("hnswSearcher.knnsearch: Y must be a finite numeric matrix.");
      endif
      if (C != size (Y, 2))
        error (strcat ("hnswSearcher.knnsearch: Y must have the same", ...
                       " number of columns as the training data in OBJ.X."));
      endif

      ## Default values
      K = 1;
      SearchSetSize = max (10, C);

      ## Parse options
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case "k"
            K = varargin{2};
            if (! (isscalar (K) && isnumeric (K) &&
                   K >= 1 && K == fix (K) && isfinite (K)))
              error ("hnswSearcher.knnsearch: 'K' must be a positive integer.");
            endif
          case "searchsetsize"
            SearchSetSize = varargin{2};
            if (! (isscalar (SearchSetSize) &&
                   isnumeric (SearchSetSize) &&
                   SearchSetSize >= 1 &&
                   SearchSetSize == fix (SearchSetSize) &&
                   isfinite (SearchSetSize)))
              error (strcat ("hnswSearcher.knnsearch: 'SearchSetSize' must", ...
                             " be a positive integer."));
            endif
            if (SearchSetSize < C || SearchSetSize > N)
              error (strcat ("hnswSearcher.knnsearch: 'SearchSetSize' must", ...
                             " be at least the number of features but no", ...
                             " more than the sample size of the training data."));
            endif
          otherwise
            error (strcat ("hnswSearcher.knnsearch: invalid", ...
                           " parameter name: '%s'."), varargin{1});
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Search HNSW graph
      idx = cell (rows (Y), 1);
      D = cell (rows (Y), 1);
      for i = 1:rows (Y)
        [temp_idx, temp_D] = search_hnsw (obj.HNSWGraph, Y(i,:), obj.X, ...
                                          obj.Distance, obj.DistParameter, ...
                                          K, SearchSetSize);
        [sorted_D, sort_idx] = sort (temp_D);
        idx{i} = temp_idx(sort_idx);
        D{i} = sorted_D;
      endfor

      idx_mat = NaN (rows (Y), K);
      D_mat = NaN (rows (Y), K);
      for i = 1:rows (Y)
        len = min (K, length (idx{i}));
        idx_mat(i, 1:len) = idx{i}(1:len);
        D_mat(i, 1:len) = D{i}(1:len);
      endfor
      idx = idx_mat;
      D = D_mat;
    endfunction

  endmethods

endclassdef

## Private Function to Build HNSW graph
function graph = build_hnsw (X, dist, param, MaxNumLinksPerNode, TrainSetSize)
  N = size (X, 1);
  if (N < 1)
    error ("build_hnsw: X must have at least one point.");
  endif
  max_layers = floor (log2 (max (N, 2))) + 1;
  graph.layers = cell (max_layers, 1);
  graph.entry_point = 1;
  mL = 1 / log (MaxNumLinksPerNode);

  ## Initialize graph with empty adjacency lists
  for L = 1:max_layers
    graph.layers{L} = cell (N, 1);
  endfor

  ## Handle single-point case
  if (N == 1)
    graph.layers{1}{1} = [];
    return;
  endif

  ## Add points to graph
  for i = 1:N
    layer = min (max_layers - 1, floor (-log (rand ()) * mL));
    for L = 0:layer
      if (i == 1)
        graph.layers{L+1}{i} = [];
        continue;
      endif
      ## Find nearest neighbors in current layer
      [neighbors, dists] = search_hnsw_layer (graph, X(i,:), X, dist, param, ...
                                              MaxNumLinksPerNode, TrainSetSize, L);
      graph.layers{L+1}{i} = neighbors;
      ## Update neighbors' connections
      for j = neighbors
        if (length (graph.layers{L+1}{j}) < MaxNumLinksPerNode)
          graph.layers{L+1}{j} = [graph.layers{L+1}{j}, i];
        else
          ## Select MaxNumLinksPerNode closest neighbors
          all_neighbors = [graph.layers{L+1}{j}, i];
          dists_j = pdist2 (X(all_neighbors,:), X(j,:), dist, param);
          [~, sort_idx] = sort (dists_j);
          graph.layers{L+1}{j} = all_neighbors(sort_idx(1:MaxNumLinksPerNode));
        endif
      endfor
      if (L == layer && i > 1)
        graph.entry_point = i;
      endif
    endfor
  endfor
endfunction

## Private Function to Search HNSW graph for k nearest neighbors
function [indices, distances] = search_hnsw (graph, Y, X, dist, param, ...
                                             K, SearchSetSize)
  max_layers = length (graph.layers);
  current_point = graph.entry_point;
  candidates = current_point;
  distances = pdist2 (X(current_point,:), Y, dist, param);

  ## Navigate to the lowest layer
  for L = max_layers:-1:2
    [new_candidates, new_dists] = search_hnsw_layer (graph, Y, X, ...
                                                     dist, param, ...
                                                     1, SearchSetSize, L-1);
    candidates = new_candidates;
    distances = new_dists;
    current_point = candidates(1);
  endfor

  ## Search in the base layer
  [indices, distances] = search_hnsw_layer (graph, Y, X, dist, param, ...
                                            K, SearchSetSize, 0);
endfunction

## Private Function Search a single HNSW layer
function [indices, distances] = search_hnsw_layer (graph, Y, X, dist, param, ...
                                                   Points, SetSize, L)
  visited = false (size (X, 1), 1);
  candidates = [graph.entry_point];
  dists = pdist2 (X(graph.entry_point,:), Y, dist, param);
  visited(graph.entry_point) = true;
  best_candidates = candidates;
  best_dists = dists;

  while (! isempty (candidates) && ! isempty (dists))
    [~, idx] = min (dists);
    closest = candidates(idx);
    candidates(idx) = [];
    dists(idx) = [];

    if (length (best_dists) > 0 && ! isempty (dists)
                                && best_dists(end) < min (dists))
      break;
    endif

    neighbors = graph.layers{L+1}{closest};
    for n = neighbors
      if (! visited(n))
        visited(n) = true;
        d = pdist2 (X(n,:), Y, dist, param);

        candidates = [candidates, n];
        dists = [dists, d];
        best_candidates = [best_candidates, n];
        best_dists = [best_dists, d];
        if (length (best_dists) > SetSize)
          [best_dists, sort_idx] = sort (best_dists);
          best_candidates = best_candidates(sort_idx);
          best_dists = best_dists(1:SetSize);
          best_candidates = best_candidates(1:SetSize);
        endif

      endif
    endfor
  endwhile

  if (length (best_dists) > Points)
    [best_dists, sort_idx] = sort (best_dists);
    best_candidates = best_candidates(sort_idx);
    indices = best_candidates(1:Points);
    distances = best_dists(1:Points);
  else
    indices = best_candidates;
    distances = best_dists;
  endif
endfunction

%!demo
%! ## Create an hnswSearcher with Euclidean distance
%! X = [1, 2; 3, 4; 5, 6];
%! obj = hnswSearcher (X);
%! ## Find the nearest neighbor to [2, 3]
%! Y = [2, 3];
%! [idx, D] = knnsearch (obj, Y, "K", 1);
%! disp ("Nearest neighbor index:");
%! disp (idx);
%! disp ("Distance:");
%! disp (D);

%!demo
%! ## Create an hnswSearcher with Minkowski distance (P=3)
%! X = [0, 0; 1, 0; 2, 0];
%! obj = hnswSearcher (X, "Distance", "minkowski", "P", 3);
%! ## Find the nearest neighbor to [1, 0]
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, "K", 1);
%! disp ("Nearest neighbor index:");
%! disp (idx);
%! disp ("Distance:");
%! disp (D);

## Test Cases

%!test
%! ## Basic constructor with default Euclidean
%! X = [1, 2; 3, 4; 5, 6];
%! obj = hnswSearcher (X);
%! assert (obj.X, X);
%! assert (obj.Distance, "euclidean");
%! assert (isempty (obj.DistParameter));

%!test
%! ## Minkowski distance with custom P
%! X = [0, 0; 1, 1; 2, 2];
%! obj = hnswSearcher (X, "Distance", "minkowski", "P", 3);
%! assert (obj.Distance, "minkowski");
%! assert (obj.DistParameter, 3);

%!test
%! ## Seuclidean distance with custom Scale
%! X = [1, 2; 3, 4; 5, 6];
%! S = [1, 2];
%! obj = hnswSearcher (X, "Distance", "seuclidean", "Scale", S);
%! assert (obj.Distance, "seuclidean");
%! assert (obj.DistParameter, S);

%!test
%! ## Mahalanobis distance with custom Cov
%! X = [1, 2; 3, 4; 5, 6];
%! C = [1, 0; 0, 1];
%! obj = hnswSearcher (X, "Distance", "mahalanobis", "Cov", C);
%! assert (obj.Distance, "mahalanobis");
%! assert (obj.DistParameter, C);

%!test
%! ## knnsearch with Euclidean distance
%! X = [1, 2; 3, 4; 5, 6];
%! obj = hnswSearcher (X);
%! Y = [2, 3];
%! [idx, D] = knnsearch (obj, Y, "K", 1);
%! assert (ismember (idx, [2]));
%! assert (abs (D - sqrt(2)) < 1e-2);

%!test
%! ## knnsearch with Cityblock distance
%! X = [0, 0; 1, 1; 2, 2];
%! obj = hnswSearcher (X, "Distance", "cityblock");
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, "K", 1);
%! assert (ismember (idx, [1, 2]));
%! assert (abs (D - 1) < 1e-2);

%!test
%! ## knnsearch with Chebychev distance
%! X = [1, 1; 2, 3; 4, 2];
%! obj = hnswSearcher (X, "Distance", "chebychev");
%! Y = [2, 2];
%! [idx, D] = knnsearch (obj, Y, "K", 1);
%! assert (ismember (idx, [1, 2]));
%! assert (abs (D - 1) < 1e-2);

%!test
%! ## knnsearch with Minkowski P=3
%! X = [0, 0; 1, 0; 2, 0];
%! obj = hnswSearcher (X, "Distance", "minkowski", "P", 3);
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, "K", 1);
%! assert (ismember (idx, [2]));
%! assert (abs (D - 0) < 1e-2);

%!test
%! ## Diverse dataset with Euclidean
%! X = [0, 10; 5, 5; 10, 0];
%! obj = hnswSearcher (X);
%! Y = [5, 5];
%! [idx, D] = knnsearch (obj, Y, "K", 1);
%! assert (ismember (idx, [2]));
%! assert (abs (D - 0) < 1e-2);

%!test
%! ## High-dimensional data with Cityblock
%! X = [1, 2, 3; 4, 5, 6; 7, 8, 9];
%! obj = hnswSearcher (X, "Distance", "cityblock");
%! Y = [4, 5, 6];
%! [idx, D] = knnsearch (obj, Y, "K", 1);
%! assert (ismember (idx, [2]));
%! assert (abs (D - 0) < 1e-2);

## Test Input Validation

%!error<hnswSearcher: too few input arguments.> ...
%! hnswSearcher ()
%!error<hnswSearcher: Name-Value arguments must be in pairs.> ...
%! hnswSearcher (ones(3,2), "Distance")
%!error<hnswSearcher: X cannot be empty.> ...
%! hnswSearcher ([])
%!error<hnswSearcher: X must be a finite numeric matrix.> ...
%! hnswSearcher ("abc")
%!error<hnswSearcher: X must be a finite numeric matrix.> ...
%! hnswSearcher ([1; Inf; 3])
%!error<hnswSearcher: invalid parameter name: 'foo'.> ...
%! hnswSearcher (ones(3,2), "foo", "bar")
%!error<hnswSearcher: unsupported distance metric 'invalid'.> ...
%! hnswSearcher (ones(3,2), "Distance", "invalid")
%!error<hnswSearcher: 'Distance' must be a string.> ...
%! hnswSearcher (ones(3,2), "Distance", 1)
%!error<hnswSearcher: 'P' must be a positive finite scalar.> ...
%! hnswSearcher (ones(3,2), "Distance", "minkowski", "P", -1)
%!error<hnswSearcher: 'Scale' must be a nonnegative vector matching X columns.> ...
%! hnswSearcher (ones(3,2), "Distance", "seuclidean", "Scale", [-1, 1])
%!error<hnswSearcher: 'Cov' must be a square matrix matching X columns.> ...
%! hnswSearcher (ones(3,2), "Distance", "mahalanobis", "Cov", ones(3,3))
%!error<hnswSearcher: 'Cov' must be symmetric for mahalanobis.> ...
%! hnswSearcher (ones(3,2), "Distance", "mahalanobis", "Cov", [1, 2; 3, 4])
%!error<hnswSearcher: 'Cov' must be positive definite for mahalanobis.> ...
%! hnswSearcher (ones(3,2), "Distance", "mahalanobis", "Cov", -eye(2))
%!error<hnswSearcher: 'MaxNumLinksPerNode' must be a positive integer.> ...
%! hnswSearcher (ones(3,2), "MaxNumLinksPerNode", 0)
%!error<hnswSearcher: 'TrainSetSize' must be a positive integer.> ...
%! hnswSearcher (ones(3,2), "TrainSetSize", -1)
%!error<hnswSearcher: 'TrainSetSize' cannot exceed the number of rows in X.> ...
%! hnswSearcher (ones(3,2), "TrainSetSize", 4)
%!error<hnswSearcher: 'TrainSetSize' cannot exceed the number of rows in X.> ...
%! hnswSearcher (ones(3,2), "MaxNumLinksPerNode", 200, "TrainSetSize", 100)

%!error<hnswSearcher.knnsearch: too few input arguments.> ...
%! knnsearch (hnswSearcher (ones(3,2)))
%!error<hnswSearcher.knnsearch: Name-Value arguments must be in pairs.> ...
%! knnsearch (hnswSearcher (ones(3,2)), ones(3,2), "K")
%!error<hnswSearcher.knnsearch: Y cannot be empty.> ...
%! knnsearch (hnswSearcher (ones(3,2)), [])
%!error<hnswSearcher.knnsearch: Y must be a finite numeric matrix.> ...
%! knnsearch (hnswSearcher (ones(3,2)), "abc")
%!error<hnswSearcher.knnsearch: Y must have the same number of columns as the training data in OBJ.X.> ...
%! knnsearch (hnswSearcher (ones(3,2)), ones(3,3))
%!error<hnswSearcher.knnsearch: 'K' must be a positive integer.> ...
%! knnsearch (hnswSearcher (ones(3,2)), ones(3,2), "K", 0)
%!error<hnswSearcher.knnsearch: invalid parameter name: 'foo'.> ...
%! knnsearch (hnswSearcher (ones(3,2)), ones(3,2), "foo", "bar")
%!error<hnswSearcher.subsref: \(\) indexing not supported.> ...
%! obj = hnswSearcher (ones(3,2)); obj(1)
%!error<hnswSearcher.subsref: {} indexing not supported.> ...
%! obj = hnswSearcher (ones(3,2)); obj{1}
%!error<hnswSearcher.subsref: unrecognized property: 'invalid'.> ...
%! obj = hnswSearcher (ones(3,2)); obj.invalid

%!error<hnswSearcher.subsasgn: \(\) indexing not supported.> ...
%! obj = hnswSearcher (ones(3,2)); obj(1) = 1
%!error<hnswSearcher.subsasgn: {} indexing not supported.> ...
%! obj = hnswSearcher (ones(3,2)); obj{1} = 1
%!error<hnswSearcher.subsasgn: 'X' is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2)); obj.X = 1


%!error<hnswSearcher.subsasgn: 'HNSWGraph' is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2)); obj.HNSWGraph = 1
%!error<hnswSearcher.subsasgn: 'Distance' is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2)); obj.Distance = "invalid"
%!error<hnswSearcher.subsasgn: 'Distance' is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2)); obj.Distance = 1
%!error<hnswSearcher.subsasgn: 'DistParameter' is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2), "Distance", "minkowski"); obj.DistParameter = -1
%!error<hnswSearcher.subsasgn: 'DistParameter' is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2), "Distance", "seuclidean"); obj.DistParameter = [-1, 1]
%!error<pdist2: covariance matrix for mahalanobis distance must be symmetric and positive definite.> ...
%! obj = hnswSearcher (ones(3,2), "Distance", "mahalanobis"); obj.DistParameter = ones(3,3)
%!error<pdist2: covariance matrix for mahalanobis distance must be symmetric and positive definite.> ...
%! obj = hnswSearcher (ones(3,2), "Distance", "mahalanobis"); obj.DistParameter = -eye(2)
%!error<hnswSearcher.subsasgn: 'DistParameter' is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2)); obj.DistParameter = 1
%!error<hnswSearcher.subsasgn: 'MaxNumLinksPerNode' is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2)); obj.MaxNumLinksPerNode = 0
%!error<hnswSearcher.subsasgn: 'TrainSetSize' is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2)); obj.TrainSetSize = -1
%!error<hnswSearcher.subsasgn: unrecognized property: 'efSearch'.> ...
%! obj = hnswSearcher (ones(3,2)); obj.efSearch = 1.5
%!error<hnswSearcher.subsasgn: unrecognized property: 'invalid'.> ...
%! obj = hnswSearcher (ones(3,2)); obj.invalid = 1
