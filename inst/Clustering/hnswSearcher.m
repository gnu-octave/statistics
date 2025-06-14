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
## search using @code{knnsearch} or a radius search using @code{rangesearch}.
##
## You can either use the @code{hnswSearcher} class constructor or the
## @code{createns} function to create an @qcode{hnswSearcher} object.
##
## @seealso{createns, ExhaustiveSearcher, KDTreeSearcher, knnsearch,
## rangesearch, pdist2}
## @end deftp

	properties (SetAccess = private)
    ## -*- texinfo -*-
    ## @deftp {Property} X
    ##
    ## Training data, specified as an @math{NxP} numeric matrix where each row
    ## is an observation and each column is a feature.  This property is private
    ## and cannot be modified after object creation.
    ##
    ## @end deftp
    X = []

    ## -*- texinfo -*-
    ## @deftp {Property} HNSWGraph
    ##
    ## The HNSW graph structure built from the training data.  This property is
    ## private and cannot be modified after object creation.
    ##
    ## @end deftp
    HNSWGraph
  endproperties

  properties
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
    ## @deftp {Property} M
    ##
    ## Maximum number of neighbors per node in the HNSW graph. Affects graph
    ## connectivity and search accuracy. Default is 16.
    ##
    ## @end deftp
    M = 16

    ## -*- texinfo -*-
    ## @deftp {Property} efConstruction
    ##
    ## Size of the dynamic candidate list during graph construction. Higher
    ## values improve accuracy at the cost of construction time. Default is 100.
    ##
    ## @end deftp
    efConstruction = 100

    ## -*- texinfo -*-
    ## @deftp {Property} efSearch
    ##
    ## Size of the dynamic candidate list during search. Higher values improve
    ## accuracy at the cost of search time. Default is 50.
    ##
    ## @end deftp
    efSearch = 50
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
        fprintf ("%+25s: %d\n", 'M', this.M);
        fprintf ("%+25s: %d\n", 'efConstruction', this.efConstruction);
        fprintf ("%+25s: %d\n", 'efSearch', this.efSearch);
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
              error (strcat ("hnswSearcher.subsasgn: X is", ...
                             " read-only and cannot be modified."));
            case 'HNSWGraph'
              error (strcat ("hnswSearcher.subsasgn: HNSWGraph is", ...
                             " read-only and cannot be modified."));
            case 'Distance'
              vm = {'euclidean', 'minkowski', 'seuclidean', 'mahalanobis', ...
                    'cityblock', 'manhattan', 'chebychev', 'cosine', ...
                    'correlation', 'spearman', 'hamming', 'jaccard'};
              if (ischar (val))
                if (! any (strcmpi (vm, val)))
                  error (strcat ("hnswSearcher.subsasgn:", ...
                                 " unsupported distance metric '%s'."), val);
                endif
                this.Distance = val;
              else
                error (strcat ("hnswSearcher.subsasgn: Distance", ...
                               " must be a string."));
              endif
            case 'DistParameter'
              if (strcmpi (this.Distance, "minkowski"))
                if (! (isscalar (val) && isnumeric (val)
                                      && val > 0 && isfinite (val)))
                  error (strcat ("hnswSearcher.subsasgn:", ...
                                 " DistParameter must be a positive", ...
                                 " finite scalar for minkowski."));
                endif
              elseif (strcmpi (this.Distance, "seuclidean"))
                if (! (isvector (val) && isnumeric (val) && all (val >= 0)
                                      && all (isfinite (val))
                                      && length (val) == columns (this.X)))
                  error (strcat ("hnswSearcher.subsasgn:", ...
                                 " DistParameter must be a nonnegative", ...
                                 " vector matching X columns."));
                endif
              elseif (strcmpi (this.Distance, "mahalanobis"))
                if (! (ismatrix (val) && isnumeric (val)
                                      && all (isfinite (val(:)))
                                      && rows (val) == columns (val)
                                      && rows (val) == columns (this.X)))
                  error (strcat ("hnswSearcher.subsasgn:", ...
                                 " DistParameter must be a square", ...
                                 " matrix matching X columns."));
                endif
                if (! issymmetric (val))
                  error (strcat ("hnswSearcher.subsasgn:", ...
                                 " DistParameter must be symmetric", ...
                                 " for mahalanobis."));
                endif
                [~, p] = chol (val);
                if (p != 0)
                  error (strcat ("hnswSearcher.subsasgn:", ...
                                 " DistParameter must be positive", ...
                                 " definite for mahalanobis."));
                endif
              else
                if (! isempty (val))
                  error (strcat ("hnswSearcher.subsasgn:", ...
                                 " DistParameter must be empty for this", ...
                                 " distance metric."));
                endif
              endif
              this.DistParameter = val;
            case 'M'
              if (! (isscalar (val) && isnumeric (val)
                                    && val > 0 && val == fix (val)))
                error (strcat ("hnswSearcher.subsasgn: M", ...
                               " must be a positive integer."));
              endif
              this.M = val;
            case 'efConstruction'
              if (! (isscalar (val) && isnumeric (val)
                                    && val > 0 && val == fix (val)))
                error (strcat ("hnswSearcher.subsasgn: efConstruction", ...
                               " must be a positive integer."));
              endif
              this.efConstruction = val;
            case 'efSearch'
              if (! (isscalar (val) && isnumeric (val)
                                    && val > 0 && val == fix (val)))
                error (strcat ("hnswSearcher.subsasgn: efSearch", ...
                               " must be a positive integer."));
              endif
              this.efSearch = val;
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
    ## @item @qcode{"M"} @tab @tab Maximum number of neighbors per node in the
    ## HNSW graph, a positive integer. Default is 16.
    ##
    ## @item @qcode{"efConstruction"} @tab @tab Size of the dynamic candidate
    ## list during graph construction, a positive integer. Default is 100.
    ##
    ## @item @qcode{"efSearch"} @tab @tab Size of the dynamic candidate list
    ## during search, a positive integer. Default is 50.
    ## @end multitable
    ##
    ## @seealso{hnswSearcher, knnsearch, rangesearch, createns, pdist2}
    ## @end deftypefn
    function obj = hnswSearcher (X, varargin)
      if (nargin < 1)
        error ("hnswSearcher: too few input arguments.");
      endif

      if (mod (numel (varargin), 2) != 0)
        error ("hnswSearcher: Name-Value arguments must be in pairs.");
      endif

      if (! (isnumeric (X) && ismatrix (X) && all (isfinite (X)(:))))
        error ("hnswSearcher: X must be a finite numeric matrix.");
      endif

      obj.X = X;

      ## Default values
      Distance = "euclidean";
      P = [];
      S = [];
      C = [];
      M = 16;
      efConstruction = 100;
      efSearch = 50;

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
          case "m"
            M = varargin{2};
          case "efconstruction"
            efConstruction = varargin{2};
          case "efsearch"
            efSearch = varargin{2};
          otherwise
            error (strcat ("hnswSearcher: invalid parameter", ...
                           " name: '%s'."), varargin{1});
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
        error ("hnswSearcher: Distance must be a string.");
      endif

      ## Set DistParameter
      if (strcmpi (obj.Distance, "minkowski"))
        if (isempty (P))
          obj.DistParameter = 2;
        else
          if (! (isscalar (P) && isnumeric (P) && P > 0 && isfinite (P)))
            error ("hnswSearcher: P must be a positive finite scalar.");
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
            error (strcat ("hnswSearcher: Scale must be a", ...
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
            error (strcat ("hnswSearcher: Cov must be a square", ...
                           " matrix matching X columns."));
          endif
          if (! issymmetric (C))
            error (strcat ("hnswSearcher: Cov must be symmetric", ...
                           " for mahalanobis."));
          endif
          [~, p] = chol (C);
          if (p != 0)
            error (strcat ("hnswSearcher: Cov must be positive", ...
                           " definite for mahalanobis."));
          endif
          obj.DistParameter = C;
        endif
      else
        obj.DistParameter = [];
      endif

      ## Validate and set HNSW parameters
      if (! (isscalar (M) && isnumeric (M) && M > 0 && M == fix (M)))
        error ("hnswSearcher: M must be a positive integer.");
      endif
      obj.M = M;

      if (! (isscalar (efConstruction) &&
            isnumeric (efConstruction) &&
            efConstruction > 0 && 
            efConstruction == fix (efConstruction)))
        error ("hnswSearcher: efConstruction must be a positive integer.");
      endif
      obj.efConstruction = efConstruction;

      if (! (isscalar (efSearch) && isnumeric (efSearch)
                                 && efSearch > 0
                                 && efSearch == fix (efSearch)))
        error ("hnswSearcher: efSearch must be a positive integer.");
      endif
      obj.efSearch = efSearch;

      ## Build HNSW graph
      obj.HNSWGraph = hnswSearcher.__build_hnsw__ (X, obj.Distance, ...
                                                  obj.DistParameter, M, ...
                                                  efConstruction);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {hnswSearcher} {[@var{idx}, @var{D}] =} knnsearch (@var{obj}, @var{Y}, @var{K})
    ## @deftypefnx {hnswSearcher} {[@var{idx}, @var{D}] =} knnsearch (@var{obj}, @var{Y}, @var{K}, @var{name}, @var{value})
    ##
    ## Find the @math{K} nearest neighbors in the training data to query points.
    ##
    ## @code{[@var{idx}, @var{D}] = knnsearch (@var{obj}, @var{Y}, @var{K})}
    ## returns the indices @var{idx} and distances @var{D} of the @math{K}
    ## nearest neighbors in @var{obj.X} to each point in @var{Y}, using the
    ## distance metric specified in @var{obj.Distance}.
    ##
    ## @itemize
    ## @item @var{obj} is an @qcode{hnswSearcher} object.
    ## @item @var{Y} is an @math{MxP} numeric matrix of query points, where
    ## @math{P} must match the number of columns in @var{obj.X}.
    ## @item @var{K} is a positive integer specifying the number of nearest
    ## neighbors to find.
    ## @end itemize
    ##
    ## @code{[@var{idx}, @var{D}] = knnsearch (@var{obj}, @var{Y}, @var{K}, @var{name}, @var{value})}
    ## allows additional options via name-value pairs:
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"IncludeTies"} @tab @tab Logical flag indicating whether to
    ## include all neighbors tied with the @math{K}th smallest distance. Default
    ## is @qcode{false}. If @qcode{true}, @var{idx} and @var{D} are cell arrays.
    ##
    ## @item @qcode{"SortIndices"} @tab @tab Logical flag indicating whether to
    ## sort the indices by distance. Default is @qcode{true}.
    ## @end multitable
    ##
    ## @var{idx} contains the indices of the nearest neighbors in @var{obj.X}.
    ## @var{D} contains the corresponding distances.
    ##
    ## @seealso{hnswSearcher, rangesearch, pdist2}
    ## @end deftypefn
    function [idx, D] = knnsearch (obj, Y, K, varargin)
      if (nargin < 3)
        error ("hnswSearcher.knnsearch: too few input arguments.");
      endif

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("hnswSearcher.knnsearch:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      if (! (isnumeric (Y) && ismatrix (Y) && all (isfinite (Y)(:))))
        error ("hnswSearcher.knnsearch: Y must be a finite numeric matrix.");
      endif

      if (size (obj.X, 2) != size (Y, 2))
        error (strcat ("hnswSearcher.knnsearch:", ...
                       " number of columns in X and Y must match."));
      endif

      if (! (isscalar (K) && isnumeric (K) && K >= 1
                                           && K == fix (K)
                                           && isfinite (K)))
        error ("hnswSearcher.knnsearch: K must be a positive integer.");
      endif

      ## Validate DistParameter for the distance metric
      if (! strcmpi (obj.Distance, "minkowski") &&
          ! strcmpi (obj.Distance, "seuclidean") &&
          ! strcmpi (obj.Distance, "mahalanobis"))
        if (! isempty (obj.DistParameter))
          error (strcat ("hnswSearcher.knnsearch: DistParameter", ...
                         " must be empty for distance metric '%s'."), ...
                         obj.Distance);
        endif
      endif

      ## Parse options
      IncludeTies = false;
      SortIndices = true;
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case "includeties"
            IncludeTies = varargin{2};
            if (! (islogical (IncludeTies) && isscalar (IncludeTies)))
              error (strcat ("hnswSearcher.knnsearch:", ...
                             " IncludeTies must be a logical scalar."));
            endif
          case "sortindices"
            SortIndices = varargin{2};
            if (! (islogical (SortIndices) && isscalar (SortIndices)))
              error (strcat ("hnswSearcher.knnsearch:", ...
                             " SortIndices must be a logical scalar."));
            endif
          otherwise
            error (strcat ("hnswSearcher.knnsearch: invalid", ...
                           " parameter name: '%s'."), varargin{1});
        endswitch
        varargin (1:2) = [];
      endwhile

      idx = cell (rows (Y), 1);
      D = cell (rows (Y), 1);
      for i = 1:rows (Y)
        [temp_idx, temp_D] = hnswSearcher.__search_hnsw__ (obj.HNSWGraph, ...
                                                           Y(i,:), K, ...
                                                           obj.X, ...
                                                           obj.Distance, ...
                                                           obj.DistParameter, ...
                                                           obj.efSearch, ...
                                                           IncludeTies);
        if (SortIndices)
          [sorted_D, sort_idx] = sort (temp_D);
          idx{i} = temp_idx(sort_idx);
          D{i} = sorted_D;
        else
          idx{i} = temp_idx;
          D{i} = temp_D;
        endif
      endfor

      if (! IncludeTies)
        idx_mat = NaN (rows (Y), K);
        D_mat = NaN (rows (Y), K);
        for i = 1:rows (Y)
          len = min (K, length (idx{i}));
          idx_mat(i, 1:len) = idx{i}(1:len);
          D_mat(i, 1:len) = D{i}(1:len);
        endfor
        idx = idx_mat;
        D = D_mat;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {hnswSearcher} {[@var{idx}, @var{D}] =} rangesearch (@var{obj}, @var{Y}, @var{r})
    ## @deftypefnx {hnswSearcher} {[@var{idx}, @var{D}] =} rangesearch (@var{obj}, @var{Y}, @var{r}, @var{name}, @var{value})
    ##
    ## Find all neighbors within a specified radius of query points.
    ##
    ## @code{[@var{idx}, @var{D}] = rangesearch (@var{obj}, @var{Y}, @var{r})}
    ## returns the indices @var{idx} and distances @var{D} of all points in
    ## @var{obj.X} within radius @var{r} of each point in @var{Y}, using the
    ## distance metric specified in @var{obj.Distance}.
    ##
    ## @itemize
    ## @item @var{obj} is an @qcode{hnswSearcher} object.
    ## @item @var{Y} is an @math{MxP} numeric matrix of query points, where
    ## @math{P} must match the number of columns in @var{obj.X}.
    ## @item @var{r} is a nonnegative scalar specifying the search radius.
    ## @end itemize
    ##
    ## @code{[@var{idx}, @var{D}] = rangesearch (@var{obj}, @var{Y}, @var{r}, @var{name}, @var{value})}
    ## allows additional options via name-value pairs:
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"SortIndices"} @tab @tab Logical flag indicating whether to
    ## sort the indices by distance. Default is @qcode{true}.
    ## @end multitable
    ##
    ## @var{idx} and @var{D} are cell arrays where each cell contains the
    ## indices and distances for one query point in @var{Y}.
    ##
    ## @seealso{hnswSearcher, knnsearch, pdist2}
    ## @end deftypefn
    function [idx, D] = rangesearch (obj, Y, r, varargin)
      if (nargin < 3)
        error ("hnswSearcher.rangesearch: too few input arguments.");
      endif

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("hnswSearcher.rangesearch:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      if (! (isnumeric (Y) && ismatrix (Y) && all (isfinite (Y)(:))))
        error (strcat ("hnswSearcher.rangesearch:", ...
                       " Y must be a finite numeric matrix."));
      endif

      if (size (obj.X, 2) != size (Y, 2))
        error (strcat ("hnswSearcher.rangesearch:", ...
                       " number of columns in X and Y must match."));
      endif

      if (! (isscalar (r) && isnumeric (r) && r >= 0 && isfinite (r)))
        error (strcat ("hnswSearcher.rangesearch:", ...
                       " R must be a nonnegative finite scalar."));
      endif

      ## Validate DistParameter for the distance metric
      if (! strcmpi (obj.Distance, "minkowski") && 
          ! strcmpi (obj.Distance, "seuclidean") &&
          ! strcmpi (obj.Distance, "mahalanobis"))
        if (! isempty (obj.DistParameter))
          error (strcat ("hnswSearcher.rangesearch: DistParameter", ...
                         " must be empty for distance metric '%s'."), ...
                         obj.Distance);
        endif
      endif

      ## Parse options
      SortIndices = true;
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case "sortindices"
            SortIndices = varargin{2};
            if (! (islogical (SortIndices) && isscalar (SortIndices)))
              error (strcat ("hnswSearcher.rangesearch:", ...
                             " SortIndices must be a logical scalar."));
            endif
          otherwise
            error (strcat ("hnswSearcher.rangesearch:", ...
                           " invalid parameter name: '%s'."), varargin{1});
        endswitch
        varargin (1:2) = [];
      endwhile

      idx = cell (rows (Y), 1);
      D = cell (rows (Y), 1);
      for i = 1:rows (Y)
        [idx{i}, D{i}] = hnswSearcher.__search_hnsw_range__ (obj.HNSWGraph, ...
                                                             Y(i,:), r, ...
                                                             obj.X, ...
                                                             obj.Distance, ...
                                                             obj.DistParameter, ...
                                                             obj.efSearch);
        if (SortIndices)
          [sorted_D, sort_idx] = sort (D{i});
          D{i} = sorted_D;
          idx{i} = idx{i}(sort_idx);
        endif
      endfor
    endfunction

  endmethods

  methods (Static, Access = private)

    ## Build HNSW graph
    function graph = __build_hnsw__ (X, dist, distparam, M, efConstruction)
      n = size (X, 1);
      if (n < 1)
        error ("hnswSearcher.__build_hnsw__: X must have at least one point.");
      endif
      max_layers = floor (log2 (max (n, 2))) + 1;
      graph.layers = cell (max_layers, 1);
      graph.entry_point = 1;
      mL = 1 / log (M);

      ## Initialize graph with empty adjacency lists
      for l = 1:max_layers
        graph.layers{l} = cell (n, 1);
      endfor

      ## Handle single-point case
      if (n == 1)
        graph.layers{1}{1} = [];
        return;
      endif

      ## Add points to graph
      for i = 1:n
        layer = min (max_layers - 1, floor (-log (rand ()) * mL));
        for l = 0:layer
          if (i == 1)
            graph.layers{l+1}{i} = [];
            continue;
          endif
          ## Find nearest neighbors in current layer
          [neighbors, dists] = hnswSearcher.__search_hnsw_layer__ (graph, ...
                                                                   X(i,:), ...
                                                                   M, X, ...
                                                                   dist, ...
                                                                   distparam, ...
                                                                   efConstruction, ...
                                                                   l);
          graph.layers{l+1}{i} = neighbors;
          ## Update neighbors' connections
          for j = neighbors
            if (length (graph.layers{l+1}{j}) < M)
              graph.layers{l+1}{j} = [graph.layers{l+1}{j}, i];
            else
              ## Select M closest neighbors
              all_neighbors = [graph.layers{l+1}{j}, i];
              dists_j = pdist2 (X(all_neighbors,:), X(j,:), dist, distparam);
              [~, sort_idx] = sort (dists_j);
              graph.layers{l+1}{j} = all_neighbors(sort_idx(1:M));
            endif
          endfor
          if (l == layer && i > 1)
            graph.entry_point = i;
          endif
        endfor
      endfor
    endfunction

    ## Search HNSW graph for k nearest neighbors
    function [indices, distances] = __search_hnsw__ (graph, query, k, X, ...
                                                     dist, distparam, ...
                                                     efSearch, include_ties)
      max_layers = length (graph.layers);
      current_point = graph.entry_point;
      candidates = current_point;
      distances = pdist2 (X(current_point,:), query, dist, distparam);

      ## Navigate to the lowest layer
      for l = max_layers:-1:2
        [new_candidates, new_dists] = hnswSearcher.__search_hnsw_layer__ (graph, ...
                                                                          query, ...
                                                                          1, X, ...
                                                                          dist, ...
                                                                          distparam, ...
                                                                          efSearch, ...
                                                                          l-1);
        candidates = new_candidates;
        distances = new_dists;
        current_point = candidates(1);
      endfor

      ## Search in the base layer
      [indices, distances] = hnswSearcher.__search_hnsw_layer__ (graph, ...
                                                                 query, k, X,...
                                                                  dist, ...
                                                                 distparam, ...
                                                                 efSearch, 0);

      if (include_ties)
        r = distances(end) + 1e-10;
        [indices, distances] = hnswSearcher.__search_hnsw_layer__ ( ...
                               graph, query, Inf, X, dist, distparam, ...
                               efSearch, 0, r);
      endif
    endfunction

    ## Search HNSW graph for points within radius r
    function [indices, distances] = __search_hnsw_range__ (graph, query, r, ...
                                                           X, dist, ...
                                                           distparam, efSearch)
      max_layers = length (graph.layers);
      current_point = graph.entry_point;
      candidates = current_point;
      distances = pdist2 (X(current_point,:), query, dist, distparam);

      ## Navigate to the lowest layer
      for l = max_layers:-1:2
        [new_candidates, new_dists] = hnswSearcher.__search_hnsw_layer__ (graph, ...
                                                                          query, ...
                                                                          1, X, ...
                                                                          dist, ...
                                                                          distparam, ...
                                                                          efSearch, ...
                                                                          l-1);
        candidates = new_candidates;
        distances = new_dists;
        current_point = candidates(1);
      endfor

      ## Search in the base layer with radius
      [indices, distances] = hnswSearcher.__search_hnsw_layer__ (graph, ...
                                                                 query, ...
                                                                 Inf, X, ...
                                                                 dist, ...
                                                                 distparam, ...
                                                                 efSearch, ...
                                                                 0, r);
    endfunction

    ## Search a single HNSW layer
    function [indices, distances] = __search_hnsw_layer__ (graph, query, k, ...
                                                           X, dist, ...
                                                           distparam, ...
                                                           efSearch, layer, r)
      if (nargin < 9)
        r = Inf;
      endif
      visited = false (size (X, 1), 1);
      candidates = [graph.entry_point];
      dists = pdist2 (X(graph.entry_point,:), query, dist, distparam);
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

        neighbors = graph.layers{layer+1}{closest};
        for n = neighbors
          if (! visited(n))
            visited(n) = true;
            d = pdist2 (X(n,:), query, dist, distparam);
            if (d <= r)
              candidates = [candidates, n];
              dists = [dists, d];
              best_candidates = [best_candidates, n];
              best_dists = [best_dists, d];
              if (length (best_dists) > efSearch && k != Inf)
                [best_dists, sort_idx] = sort (best_dists);
                best_candidates = best_candidates(sort_idx);
                best_dists = best_dists(1:efSearch);
                best_candidates = best_candidates(1:efSearch);
              endif
            endif
          endif
        endfor
      endwhile

      if (k == Inf)
        indices = best_candidates(best_dists <= r);
        distances = best_dists(best_dists <= r);
      else
        if (length (best_dists) > k)
          [best_dists, sort_idx] = sort (best_dists);
          best_candidates = best_candidates(sort_idx);
          indices = best_candidates(1:k);
          distances = best_dists(1:k);
        else
          indices = best_candidates;
          distances = best_dists;
        endif
      endif
    endfunction

  endmethods

endclassdef

%!demo
%! ## Create an hnswSearcher with Euclidean distance
%! X = [1, 2; 3, 4; 5, 6];
%! obj = hnswSearcher (X);
%! ## Find the nearest neighbor to [2, 3]
%! Y = [2, 3];
%! [idx, D] = knnsearch (obj, Y, 1);
%! disp ("Nearest neighbor index:");
%! disp (idx);
%! disp ("Distance:");
%! disp (D);
%! ## Find all points within radius 2
%! [idx, D] = rangesearch (obj, Y, 2);
%! disp ("Indices within radius:");
%! disp (idx);
%! disp ("Distances:");
%! disp (D);

%!demo
%! ## Create an hnswSearcher with Minkowski distance (P=3)
%! X = [0, 0; 1, 0; 2, 0];
%! obj = hnswSearcher (X, "Distance", "minkowski", "P", 3);
%! ## Find the nearest neighbor to [1, 0]
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, 1);
%! disp ("Nearest neighbor index:");
%! disp (idx);
%! disp ("Distance:");
%! disp (D);

%!demo
%! rng(42);
%! disp('Demonstrating hnswSearcher');
%!
%! n = 100;
%! mu1 = [0.3, 0.3];
%! mu2 = [0.7, 0.7];
%! sigma = 0.1;
%! X1 = mu1 + sigma * randn (n / 2, 2);
%! X2 = mu2 + sigma * randn (n / 2, 2);
%! X = [X1; X2];
%!
%! obj = hnswSearcher(X, 'M', 10, 'efConstruction', 50, 'efSearch', 20);
%!
%! Y = [0.3, 0.3; 0.7, 0.7; 0.5, 0.5];
%!
%! K = 5;
%! [idx, D] = knnsearch (obj, Y, K);
%!
%! disp ('For the first query point:');
%! disp (['Query point: ', num2str(Y(1,:))]);
%! disp ('Indices of nearest neighbors:');
%! disp (idx(1,:));
%! disp ('Distances:');
%! disp (D(1,:));
%!
%! figure;
%! scatter (X(:,1), X(:,2), 36, 'b', 'filled'); # Training points
%! hold on;
%! scatter (Y(:,1), Y(:,2), 36, 'r', 'filled'); # Query points
%! for i = 1:size (Y, 1)
%!     query = Y(i,:);
%!     neighbors = X(idx(i,:), :);
%!     for j = 1:K
%!         plot ([query(1), neighbors(j,1)], [query(2), neighbors(j,2)], 'k-');
%!     endfor
%! endfor
%! hold off;
%! title ('K Nearest Neighbors with hnswSearcher');
%! xlabel ('X1');
%! ylabel ('X2');
%!
%! r = 0.15;
%! [idx, D] = rangesearch (obj, Y, r);
%!
%! disp ('For the first query point in rangesearch:');
%! disp (['Query point: ', num2str(Y(1,:))]);
%! disp ('Indices of points within radius:');
%! disp (idx{1});
%! disp ('Distances:');
%! disp (D{1});
%!
%! figure;
%! scatter (X(:,1), X(:,2), 36, 'b', 'filled');
%! hold on;
%! scatter (Y(:,1), Y(:,2), 36, 'r', 'filled');
%! theta = linspace (0, 2 * pi, 100);
%! for i = 1:size (Y, 1)
%!     center = Y(i,:);
%!     x_circle = center(1) + r * cos (theta);
%!     y_circle = center(2) + r * sin (theta);
%!     plot (x_circle, y_circle, 'g-');
%!     ## Highlight points within radius
%!     if (! isempty (idx{i}))
%!       in_radius = X(idx{i}, :);
%!       scatter (in_radius(:,1), in_radius(:,2), 36, 'g', 'filled');
%!     endif
%! endfor
%! hold off;
%! title ('Points within Radius with hnswSearcher');
%! xlabel ('X1');
%! ylabel ('X2');

## Test Cases

%!test
%! ## Basic constructor with default Euclidean
%! X = [1, 2; 3, 4; 5, 6];
%! obj = hnswSearcher (X);
%! assert (obj.X, X);
%! assert (obj.Distance, "euclidean");
%! assert (isempty (obj.DistParameter));
%! assert (obj.M, 16);
%! assert (obj.efConstruction, 100);
%! assert (obj.efSearch, 50);

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
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (ismember (idx, [2]));
%! assert (abs (D - sqrt(2)) < 1e-2);

%!test
%! ## knnsearch with Cityblock distance
%! X = [0, 0; 1, 1; 2, 2];
%! obj = hnswSearcher (X, "Distance", "cityblock");
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (ismember (idx, [1, 2]));
%! assert (abs (D - 1) < 1e-2);

%!test
%! ## knnsearch with Chebychev distance
%! X = [1, 1; 2, 3; 4, 2];
%! obj = hnswSearcher (X, "Distance", "chebychev");
%! Y = [2, 2];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (ismember (idx, [1, 2]));
%! assert (abs (D - 1) < 1e-2);

%!test
%! ## knnsearch with Minkowski P=3
%! X = [0, 0; 1, 0; 2, 0];
%! obj = hnswSearcher (X, "Distance", "minkowski", "P", 3);
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (ismember (idx, [2]));
%! assert (abs (D - 0) < 1e-2);

%!test
%! ## knnsearch with IncludeTies
%! X = [0, 0; 1, 0; 0, 1];
%! obj = hnswSearcher (X);
%! Y = [0.5, 0];
%! [idx, D] = knnsearch (obj, Y, 1, "IncludeTies", true);
%! assert (iscell (idx));
%! assert (all (ismember (idx{1}, [1, 2])));
%! assert (all (abs (D{1} - 0.5) < 1e-2));

%!test
%! ## rangesearch with Euclidean
%! X = [1, 1; 2, 2; 3, 3];
%! obj = hnswSearcher (X);
%! Y = [0, 0];
%! [idx, D] = rangesearch (obj, Y, 2);
%! assert (all (ismember (idx{1}, [1])));
%! assert (all (abs (D{1} - sqrt(2)) < 1e-2));

%!test
%! ## rangesearch with Cityblock
%! X = [0, 0; 1, 1; 2, 2];
%! obj = hnswSearcher (X, "Distance", "cityblock");
%! Y = [0, 0];
%! [idx, D] = rangesearch (obj, Y, 1);
%! assert (all (ismember (idx{1}, [1])));
%! assert (all (abs (D{1} - 0) < 1e-2));

%!test
%! ## rangesearch with Chebychev
%! X = [1, 1; 2, 3; 4, 2];
%! obj = hnswSearcher (X, "Distance", "chebychev");
%! Y = [2, 2];
%! [idx, D] = rangesearch (obj, Y, 1);
%! assert (all (ismember (idx{1}, [1, 2])));
%! assert (all (abs (D{1} - 1) < 1e-2));

%!test
%! ## rangesearch with Minkowski P=3
%! X = [0, 0; 1, 0; 2, 0];
%! obj = hnswSearcher (X, "Distance", "minkowski", "P", 3);
%! Y = [1, 0];
%! [idx, D] = rangesearch (obj, Y, 1);
%! assert (all (ismember (idx{1}, [1, 2, 3])));
%! assert (all (abs (D{1} - [0, 1, 1]) < 1e-2));

%!test
%! ## Diverse dataset with Euclidean
%! X = [0, 10; 5, 5; 10, 0];
%! obj = hnswSearcher (X);
%! Y = [5, 5];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (ismember (idx, [2]));
%! assert (abs (D - 0) < 1e-2);

%!test
%! ## High-dimensional data with Cityblock
%! X = [1, 2, 3; 4, 5, 6; 7, 8, 9];
%! obj = hnswSearcher (X, "Distance", "cityblock");
%! Y = [4, 5, 6];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (ismember (idx, [2]));
%! assert (abs (D - 0) < 1e-2);

## Test Input Validation

%!error<hnswSearcher: too few input arguments.> ...
%! hnswSearcher ()
%!error<hnswSearcher: Name-Value arguments must be in pairs.> ...
%! hnswSearcher (ones(3,2), "Distance")
%!error<hnswSearcher: X must be a finite numeric matrix.> ...
%! hnswSearcher ("abc")
%!error<hnswSearcher: X must be a finite numeric matrix.> ...
%! hnswSearcher ([1; Inf; 3])
%!error<hnswSearcher: invalid parameter name: 'foo'.> ...
%! hnswSearcher (ones(3,2), "foo", "bar")
%!error<hnswSearcher: unsupported distance metric 'invalid'.> ...
%! hnswSearcher (ones(3,2), "Distance", "invalid")
%!error<hnswSearcher: Distance must be a string.> ...
%! hnswSearcher (ones(3,2), "Distance", 1)
%!error<hnswSearcher: P must be a positive finite scalar.> ...
%! hnswSearcher (ones(3,2), "Distance", "minkowski", "P", -1)
%!error<hnswSearcher: Scale must be a nonnegative vector matching X columns.> ...
%! hnswSearcher (ones(3,2), "Distance", "seuclidean", "Scale", [-1, 1])
%!error<hnswSearcher: Cov must be a square matrix matching X columns.> ...
%! hnswSearcher (ones(3,2), "Distance", "mahalanobis", "Cov", ones(3,3))
%!error<hnswSearcher: Cov must be positive definite.> ...
%! hnswSearcher (ones(3,2), "Distance", "mahalanobis", "Cov", [1, 2; 3, 4])
%!error<hnswSearcher: Cov must be positive definite.> ...
%! hnswSearcher (ones(3,2), "Distance", "mahalanobis", "Cov", -eye(2))
%!error<hnswSearcher: M must be a positive integer.> ...
%! hnswSearcher (ones(3,2), "M", 0)
%!error<hnswSearcher: efConstruction must be a positive integer.> ...
%! hnswSearcher (ones(3,2), "efConstruction", -1)
%!error<hnswSearcher: efSearch must be a positive integer.> ...
%! hnswSearcher (ones(3,2), "efSearch", 1.5)

%!error<hnswSearcher.knnsearch: too few input arguments.> ...
%! knnsearch (hnswSearcher (ones(3,2)))
%!error<hnswSearcher.knnsearch: Name-Value arguments must be in pairs.> ...
%! knnsearch (hnswSearcher (ones(3,2)), ones(3,2), 1, "IncludeTies")
%!error<hnswSearcher.knnsearch: Y must be a finite numeric matrix.> ...
%! knnsearch (hnswSearcher (ones(3,2)), "abc", 1)
%!error<hnswSearcher.knnsearch: number of columns in X and Y must match.> ...
%! knnsearch (hnswSearcher (ones(3,2)), ones(3,3), 1)
%!error<hnswSearcher.knnsearch: K must be a positive integer.> ...
%! knnsearch (hnswSearcher (ones(3,2)), ones(3,2), 0)
%!error<hnswSearcher.knnsearch: invalid parameter name: 'foo'.> ...
%! knnsearch (hnswSearcher (ones(3,2)), ones(3,2), 1, "foo", "bar")
%!error<hnswSearcher.knnsearch: IncludeTies must be a logical scalar.> ...
%! knnsearch (hnswSearcher (ones(3,2)), ones(3,2), 1, "IncludeTies", 1)
%!error<hnswSearcher.knnsearch: SortIndices must be a logical scalar.> ...
%! knnsearch (hnswSearcher (ones(3,2)), ones(3,2), 1, "SortIndices", 1)

%!error<hnswSearcher.rangesearch: too few input arguments.> ...
%! rangesearch (hnswSearcher (ones(3,2)))
%!error<hnswSearcher.rangesearch: Name-Value arguments must be in pairs.> ...
%! rangesearch (hnswSearcher (ones(3,2)), ones(3,2), 1, "SortIndices")
%!error<hnswSearcher.rangesearch: Y must be a finite numeric matrix.> ...
%! rangesearch (hnswSearcher (ones(3,2)), "abc", 1)
%!error<hnswSearcher.rangesearch: number of columns in X and Y must match.> ...
%! rangesearch (hnswSearcher (ones(3,2)), ones(3,3), 1)
%!error<hnswSearcher.rangesearch: R must be a nonnegative finite scalar.> ...
%! rangesearch (hnswSearcher (ones(3,2)), ones(3,2), -1)
%!error<hnswSearcher.rangesearch: invalid parameter name: 'foo'.> ...
%! rangesearch (hnswSearcher (ones(3,2)), ones(3,2), 1, "foo", "bar")
%!error<hnswSearcher.rangesearch: SortIndices must be a logical scalar.> ...
%! rangesearch (hnswSearcher (ones(3,2)), ones(3,2), 1, "SortIndices", 1)

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
%!error<hnswSearcher.subsasgn: X is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2)); obj.X = 1
%!error<hnswSearcher.subsasgn: HNSWGraph is read-only and cannot be modified.> ...
%! obj = hnswSearcher (ones(3,2)); obj.HNSWGraph = 1
%!error<hnswSearcher.subsasgn: unsupported distance metric 'invalid'.> ...
%! obj = hnswSearcher (ones(3,2)); obj.Distance = "invalid"
%!error<hnswSearcher.subsasgn: Distance must be a string.> ...
%! obj = hnswSearcher (ones(3,2)); obj.Distance = 1
%!error<hnswSearcher.subsasgn: DistParameter must be a positive finite scalar for minkowski.> ...
%! obj = hnswSearcher (ones(3,2), "Distance", "minkowski"); obj.DistParameter = -1
%!error<hnswSearcher.subsasgn: DistParameter must be a nonnegative vector matching X columns.> ...
%! obj = hnswSearcher (ones(3,2), "Distance", "seuclidean"); obj.DistParameter = [-1, 1]
%!error<pdist2: covariance matrix for mahalanobis distance must be symmetric and positive definite.> ...
%! obj = hnswSearcher (ones(3,2), "Distance", "mahalanobis"); obj.DistParameter = ones(3,3)
%!error<pdist2: covariance matrix for mahalanobis distance must be symmetric and positive definite.> ...
%! obj = hnswSearcher (ones(3,2), "Distance", "mahalanobis"); obj.DistParameter = -eye(2)
%!error<hnswSearcher.subsasgn: DistParameter must be empty for this distance metric.> ...
%! obj = hnswSearcher (ones(3,2)); obj.DistParameter = 1
%!error<hnswSearcher.subsasgn: M must be a positive integer.> ...
%! obj = hnswSearcher (ones(3,2)); obj.M = 0
%!error<hnswSearcher.subsasgn: efConstruction must be a positive integer.> ...
%! obj = hnswSearcher (ones(3,2)); obj.efConstruction = -1
%!error<hnswSearcher.subsasgn: efSearch must be a positive integer.> ...
%! obj = hnswSearcher (ones(3,2)); obj.efSearch = 1.5
%!error<hnswSearcher.subsasgn: unrecognized property: 'invalid'.> ...
%! obj = hnswSearcher (ones(3,2)); obj.invalid = 1