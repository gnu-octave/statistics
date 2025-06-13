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
    ## @qcode{"euclidean"}, @qcode{"minkowski"}). Default is
    ## @qcode{"euclidean"}. Supported metrics align with those in @code{pdist2}.
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
      if (numel (s) > 1)
        error ("hnswSearcher.subsasgn: chained subscripts not allowed.");
      endif
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
    ## character vector (e.g., @qcode{"euclidean"}, @qcode{"minkowski"}).
    ## Default is @qcode{"euclidean"}. See @code{pdist2} for supported metrics.
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
          [~, p] = chol (C);
          if (p != 0)
            error ("hnswSearcher: Cov must be positive definite.");
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

      if (! (isscalar (efConstruction) && isnumeric (efConstruction)
             && efConstruction > 0
             && (efConstruction == fix (efConstruction))))
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