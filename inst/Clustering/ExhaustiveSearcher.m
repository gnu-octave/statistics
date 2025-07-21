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

classdef ExhaustiveSearcher
## -*- texinfo -*-
## @deftp {Class} ExhaustiveSearcher
##
## Exhaustive nearest neighbor searcher class.
##
## The @code{ExhaustiveSearcher} class implements an exhaustive search algorithm
## for nearest neighbor queries.  It stores training data and supports various
## distance metrics along with their parameter values for performing an
## exhaustive search.  The exhaustive search algorithm computes the distance
## from each query point to all the points in the training data and facilitates
## a nearest neighbor search using @code{knnsearch} or a radius search using
## @code{rangesearch}.
##
## You can either use the @code{ExhaustiveSearcher} class constructor or the
## @code{createns} function to create an @qcode{ExhaustiveSearcher} object.
##
## @seealso{createns, KDTreeSearcher, hnswSearcher, knnsearch, rangesearch}
## @end deftp

  properties (SetAccess = private)
    ## -*- texinfo -*-
    ## @deftp {Property} X
    ##
    ## Training data, specified as an @math{NxP} numeric matrix where each row
    ## is an observation and each column is a feature. This property is private
    ## and cannot be modified after object creation.
    ##
    ## @end deftp
    X = []
  endproperties

  properties
    ## -*- texinfo -*-
    ## @deftp {Property} Distance
    ##
    ## Distance metric used for searches, specified as a character vector (e.g.,
    ## @qcode{"euclidean"}, @qcode{"minkowski"}) or a function handle to a
    ## custom distance function. Default is @qcode{"euclidean"}.  Supported
    ## metrics align with those in @code{pdist2}.
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
    ## @item Empty for other metrics or custom functions.
    ## @end itemize
    ##
    ## @end deftp
    DistParameter = []
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
        fprintf ("\n  ExhaustiveSearcher with properties:\n\n");
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
          error ("ExhaustiveSearcher.subsref: () indexing not supported.");
        case '{}'
          error ("ExhaustiveSearcher.subsref: {} indexing not supported.");
        case '.'
          if (! ischar (s.subs))
            error (strcat ("ExhaustiveSearcher.subsref: property", ...
                           " name must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("ExhaustiveSearcher.subsref: unrecognized", ...
                           " property: '%s'"), s.subs);
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
        error ("ExhaustiveSearcher.subsasgn: chained subscripts not allowed.");
      endif
      switch s.type
        case '()'
          error ("ExhaustiveSearcher.subsasgn: () indexing not supported.");
        case '{}'
          error ("ExhaustiveSearcher.subsasgn: {} indexing not supported.");
        case '.'
          if (! ischar (s.subs))
            error (strcat ("ExhaustiveSearcher.subsasgn: property", ...
                           " name must be a character vector."));
          endif
          switch (s.subs)
            case 'X'
              error (strcat ("ExhaustiveSearcher.subsasgn: X is", ...
                             " read-only and cannot be modified."));
            case 'Distance'
              vm = {'euclidean', 'minkowski', 'seuclidean', 'mahalanobis', ...
                    'cityblock', 'manhattan', 'chebychev', 'cosine', ...
                    'correlation', 'spearman', 'hamming', 'jaccard'};
              if (ischar (val))
                if (! any (strcmpi (vm, val)))
                  error (strcat ("ExhaustiveSearcher.subsasgn:", ...
                                 " unsupported distance metric '%s'."), val);
                endif
                this.Distance = val;
              elseif (isa (val, "function_handle"))
                try
                  D = val (this.X(1,:), this.X);
                catch
                  error (strcat ("ExhaustiveSearcher.subsasgn:", ...
                                 " invalid distance function handle."));
                end_try_catch
                if (! isvector (D) || length (D) != rows (this.X))
                  error (strcat ("ExhaustiveSearcher.subsasgn: custom", ...
                                 " distance function output invalid."));
                endif
                this.Distance = val;
              else
                error (strcat ("ExhaustiveSearcher.subsasgn: Distance", ...
                               " must be a string or function handle."));
              endif
            case 'DistParameter'
              if (strcmpi (this.Distance, "minkowski"))
                if (! (isscalar (val) && isnumeric (val)
                                      && val > 0 && isfinite (val)))
                  error (strcat ("ExhaustiveSearcher.subsasgn:", ...
                                 " DistParameter must be a positive", ...
                                 " finite scalar for minkowski."));
                endif
              elseif (strcmpi (this.Distance, "seuclidean"))
                if (! (isvector (val) && isnumeric (val) && all (val >= 0)
                                      && all (isfinite (val))
                                      && length (val) == columns (this.X)))
                  error (strcat ("ExhaustiveSearcher.subsasgn:", ...
                                 " DistParameter must be a nonnegative", ...
                                 " vector matching X columns."));
                endif
              elseif (strcmpi (this.Distance, "mahalanobis"))
                if (! (ismatrix (val) && isnumeric (val)
                                      && all (isfinite (val(:)))
                                      && rows (val) == columns (val)
                                      && rows (val) == columns (this.X)))
                  error (strcat ("ExhaustiveSearcher.subsasgn:", ...
                                 " DistParameter must be a square", ...
                                 " matrix matching X columns."));
                endif
                [~, p] = chol (val);
                if (p != 0)
                  error (strcat ("ExhaustiveSearcher.subsasgn:", ...
                                 " DistParameter must be positive", ...
                                 " definite for mahalanobis."));
                endif
              else
                if (! isempty (val))
                  error (strcat ("ExhaustiveSearcher.subsasgn:", ...
                                 " DistParameter must be empty for this", ...
                                 " distance metric."));
                endif
              endif
              this.DistParameter = val;
            otherwise
              error (strcat ("ExhaustiveSearcher.subsasgn:", ...
                             " unrecognized property: '%s'"), s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods

    ## -*- texinfo -*-
    ## @deftypefn  {ExhaustiveSearcher} {@var{obj} =} ExhaustiveSearcher (@var{X})
    ## @deftypefnx {ExhaustiveSearcher} {@var{obj} =} ExhaustiveSearcher (@var{X}, @var{name}, @var{value})
    ##
    ## Create an @qcode{ExhaustiveSearcher} object for nearest neighbor
    ## searches.
    ##
    ## @code{@var{obj} = ExhaustiveSearcher (@var{X})} constructs an
    ## @qcode{ExhaustiveSearcher} object with training data @var{X} using the
    ## default @qcode{"euclidean"} distance metric. @var{X} must be an
    ## @math{NxP} numeric matrix, where rows represent observations and columns
    ## represent features.
    ##
    ## @code{@var{obj} = ExhaustiveSearcher (@var{X}, @var{name}, @var{value})}
    ## allows customization through name-value pairs:
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"Distance"} @tab @tab Distance metric, specified as a
    ## character vector (e.g., @qcode{"euclidean"}, @qcode{"minkowski"}) or a
    ## function handle.  Default is @qcode{"euclidean"}.  See @code{pdist2} for
    ## supported metrics.
    ##
    ## @item @qcode{"P"} @tab @tab a positive scalar specifying the exponent for
    ## the Minkowski distance.  Valid only when @qcode{"Distance"} is
    ## @qcode{"minkowski"}.  Default is 2.
    ##
    ## @item @qcode{"Scale"} @tab @tab a nonnegative vector with the same number
    ## of elements as the columns in @var{X} specifying the scale parameter for
    ## the standardized Euclidean distance.  Valid only when
    ## @qcode{"Distance"} is @qcode{"seuclidean"}.  Default is @code{std (X)}.
    ##
    ## @item @qcode{"Cov"} @tab @tab a positive definite matrix matching the
    ## number of columns in @var{X} specifying the covariance matrix for the
    ## Mahalanobis distance.  Valid only when @qcode{"Distance"} is
    ## @qcode{"mahalanobis"}.  Default is @code{cov (X)}.
    ## @end multitable
    ##
    ## @seealso{ExhaustiveSearcher, knnsearch, rangesearch, pdist2}
    ## @end deftypefn
    function obj = ExhaustiveSearcher (X, varargin)
      if (nargin < 1)
        error ("ExhaustiveSearcher: too few input arguments.");
      endif

      if (mod (numel (varargin), 2) != 0)
        error ("ExhaustiveSearcher: Name-Value arguments must be in pairs.");
      endif

      if (! (isnumeric (X) && ismatrix (X) && all (isfinite (X)(:))))
        error ("ExhaustiveSearcher: X must be a finite numeric matrix.");
      endif

      obj.X = X;

      ## Default values for optional parameters
      Distance = "euclidean";
      P = [];
      S = [];
      C = [];

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
          otherwise
            error (strcat ("ExhaustiveSearcher: invalid parameter", ...
                           " name: '%s'."), varargin{1});
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Validate and set distance metric
      valid_metrics = {'euclidean', 'minkowski', 'seuclidean', ...
                       'mahalanobis', 'cityblock', 'manhattan', ...
                       'chebychev', 'cosine', 'correlation', ...
                       'spearman', 'hamming', 'jaccard'};
      if (ischar (Distance))
        if (! any (strcmpi (valid_metrics, Distance)))
          error ("ExhaustiveSearcher: unsupported distance metric '%s'.", ...
                 Distance);
        endif
        obj.Distance = Distance;
      elseif (isa (Distance, "function_handle"))
        try
          D = Distance (X(1,:), X);
        catch
          error ("ExhaustiveSearcher: invalid distance function handle.");
        end_try_catch
        if (! isvector (D) || length (D) != rows (X))
          error (strcat ("ExhaustiveSearcher: custom distance", ...
                         " function output invalid."));
        endif
        obj.Distance = Distance;
      else
        error (strcat ("ExhaustiveSearcher: Distance must", ...
                       " be a string or function handle."));
      endif

      ## Set DistParameter based on Distance
      if (strcmpi (obj.Distance, "minkowski"))
        if (isempty (P))
          obj.DistParameter = 2;
        else
          if (! (isscalar (P) && isnumeric (P) && P > 0 && isfinite (P)))
            error ("ExhaustiveSearcher: P must be a positive finite scalar.");
          endif
          obj.DistParameter = P;
        endif
      elseif (strcmpi (obj.Distance, "seuclidean"))
        if (isempty (S))
          obj.DistParameter = std (X, [], 1);
        else
          if (! (isvector (S) && isnumeric (S) && all (S >= 0) && ...
                 all (isfinite (S)) && length (S) == columns (X)))
            error (strcat ("ExhaustiveSearcher: Scale must be a", ...
                           " nonnegative vector matching X columns."));
          endif
          obj.DistParameter = S;
        endif
      elseif (strcmpi (obj.Distance, "mahalanobis"))
        if (isempty (C))
          obj.DistParameter = cov (X);
        else
          if (! (ismatrix (C) && isnumeric (C) && all (isfinite (C)(:)) && ...
                 rows (C) == columns (C) && rows (C) == columns (X)))
            error (strcat ("ExhaustiveSearcher: Cov must be a square", ...
                           " matrix matching X columns."));
          endif
          [~, p] = chol (C);
          if (p != 0)
            error ("ExhaustiveSearcher: Cov must be positive definite.");
          endif
          obj.DistParameter = C;
        endif
      else
        obj.DistParameter = [];
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExhaustiveSearcher} {[@var{idx}, @var{D}] =} knnsearch (@var{obj}, @var{Y})
    ## @deftypefnx {ExhaustiveSearcher} {[@var{idx}, @var{D}] =} knnsearch (@var{obj}, @var{Y}, @var{name}, @var{value})
    ##
    ## Find the @math{K} nearest neighbors in the training data to query points.
    ##
    ## @code{[@var{idx}, @var{D}] = knnsearch (@var{obj}, @var{Y})} returns the
    ## indices @var{idx} and distances @var{D} of the nearest neighbor in
    ## @var{obj.X} to each point in @var{Y}, using the distance metric specified
    ## in @var{obj.Distance}.
    ##
    ## @itemize
    ## @item @var{obj} is an @qcode{ExhaustiveSearcher} object.
    ## @item @var{Y} is an @math{MxP} numeric matrix of query points, where
    ## @math{P} must match the number of columns in @var{obj.X}.
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
    ## @item @qcode{"IncludeTies"} @tab @tab Logical flag indicating whether to
    ## include all neighbors tied with the @math{K}th smallest distance. Default
    ## is @qcode{false}. If @qcode{true}, @var{idx} and @var{D} are cell arrays.
    ## @end multitable
    ##
    ## @var{idx} contains the indices of the nearest neighbors in @var{obj.X}.
    ## @var{D} contains the corresponding distances.
    ##
    ## @seealso{ExhaustiveSearcher, rangesearch, pdist2}
    ## @end deftypefn
    function [idx, D] = knnsearch (obj, Y, varargin)
      if (nargin < 2)
        error ("ExhaustiveSearcher.knnsearch: too few input arguments.");
      endif

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ExhaustiveSearcher.knnsearch:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      if (! (isnumeric (Y) && ismatrix (Y) && all (isfinite (Y)(:))))
        error (strcat ("ExhaustiveSearcher.knnsearch: Y", ...
                       " must be a finite numeric matrix."));
      endif

      if (size (obj.X, 2) != size (Y, 2))
        error (strcat ("ExhaustiveSearcher.knnsearch: number", ...
                       " of columns in X and Y must match."));
      endif

      K = 1;
      IncludeTies = false;

      ## Parse options
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case "k"
            K = varargin{2};
            if (! (isscalar (K) && isnumeric (K) && K >= 1
                                && K == fix (K) && isfinite (K)))
              error (strcat ("ExhaustiveSearcher.knnsearch: K", ...
                             " must be a positive integer."));
            endif
          case "includeties"
            IncludeTies = varargin{2};
            if (! (islogical (IncludeTies) && isscalar (IncludeTies)))
              error (strcat ("ExhaustiveSearcher.knnsearch:", ...
                             " IncludeTies must be a logical scalar."));
            endif
          otherwise
            error (strcat ("ExhaustiveSearcher.knnsearch:", ...
                           " invalid parameter name: '%s'."), varargin{1});
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Compute distance matrix
      D_mat = pdist2 (obj.X, Y, obj.Distance, obj.DistParameter);
      D_mat = reshape (D_mat', size (Y, 1), size (obj.X, 1));

      if (K == 1 && ! IncludeTies)
        [D, idx] = min (D_mat, [], 2);
      else
        [sorted_D, sorted_idx] = sort (D_mat, 2);
        if (IncludeTies)
          idx = cell (rows (Y), 1);
          D = cell (rows (Y), 1);
          for i = 1:rows (Y)
            if (K > columns (sorted_D))
              kth_dist = sorted_D(i, end);
            else
              kth_dist = sorted_D(i, K);
            endif
            tie_idx = find (sorted_D(i, :) <= kth_dist);
            idx{i} = sorted_idx(i, tie_idx);
            D{i} = sorted_D(i, tie_idx);
          endfor
        else
          idx = sorted_idx(:, 1:K);
          D = sorted_D(:, 1:K);
        endif
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ExhaustiveSearcher} {[@var{idx}, @var{D}] =} rangesearch (@var{obj}, @var{Y}, @var{r})
    ## @deftypefnx {ExhaustiveSearcher} {[@var{idx}, @var{D}] =} rangesearch (@var{obj}, @var{Y}, @var{r}, @var{name}, @var{value})
    ##
    ## Find all neighbors within a specified radius of query points.
    ##
    ## @code{[@var{idx}, @var{D}] = rangesearch (@var{obj}, @var{Y}, @var{r})}
    ## returns the indices @var{idx} and distances @var{D} of all points in
    ## @var{obj.X} within radius @var{r} of each point in @var{Y}, using the
    ## distance metric specified in @var{obj.Distance}.
    ##
    ## @itemize
    ## @item @var{obj} is an @qcode{ExhaustiveSearcher} object.
    ## @item @var{Y} is an @math{MxP} numeric matrix of query points, where
    ## @math{P} must match the number of columns in @var{obj.X}.
    ## @item @var{r} is a nonnegative scalar specifying the search radius.
    ## @end itemize
    ##
    ## @code{[@var{idx}, @var{D}] = rangesearch (@dots{}, @var{name},
    ## @var{value})} allows additional options via name-value pairs:
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
    ## @seealso{ExhaustiveSearcher, knnsearch, pdist2}
    ## @end deftypefn
    function [idx, D] = rangesearch (obj, Y, r, varargin)
      if (nargin < 3)
        error ("ExhaustiveSearcher.rangesearch: too few input arguments.");
      endif

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ExhaustiveSearcher.rangesearch:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      if (! (isnumeric (Y) && ismatrix (Y) && all (isfinite (Y)(:))))
        error (strcat ("ExhaustiveSearcher.rangesearch: Y", ...
                       " must be a finite numeric matrix."));
      endif

      if (size (obj.X, 2) != size (Y, 2))
        error (strcat ("ExhaustiveSearcher.rangesearch: number", ...
                       " of columns in X and Y must match."));
      endif

      if (! (isscalar (r) && isnumeric (r) && r >= 0 && isfinite (r)))
        error (strcat ("ExhaustiveSearcher.rangesearch: R", ...
                       " must be a nonnegative finite scalar."));
      endif

      ## Parse options
      SortIndices = true;
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case "sortindices"
            SortIndices = varargin{2};
            if (! (islogical (SortIndices) && isscalar (SortIndices)))
              error (strcat ("ExhaustiveSearcher.rangesearch:", ...
                             " SortIndices must be a logical scalar."));
            endif
          otherwise
            error (strcat ("ExhaustiveSearcher.rangesearch: invalid", ...
                           " parameter name: '%s'."), varargin{1});
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Compute distance matrix
      D_mat = pdist2 (obj.X, Y, obj.Distance, obj.DistParameter);
      D_mat = reshape (D_mat', size (Y, 1), size (obj.X, 1));

      idx = cell (rows (Y), 1);
      D = cell (rows (Y), 1);
      for i = 1:rows (Y)
        within_r = find (D_mat(i, :) <= r);
        if (SortIndices)
          [sorted_D, sort_idx] = sort (D_mat(i, within_r));
          idx{i} = within_r(sort_idx);
          D{i} = sorted_D;
        else
          idx{i} = within_r;
          D{i} = D_mat(i, within_r);
        endif
      endfor
    endfunction

  endmethods

endclassdef

## Demo Examples

%!demo
%! ## Demo to verify implementation using fisheriris dataset
%! load fisheriris
%! rng('default');
%! numSamples = size (meas, 1);
%! queryIndices = [20, 95, 123, 136, 138];
%! dataPoints = meas(~ismember (1:numSamples, queryIndices), :);
%! queryPoints = meas(queryIndices, :);
%! searchModel = ExhaustiveSearcher (dataPoints, 'Distance', 'mahalanobis')
%! mahalanobisParam = searchModel.DistParameter
%! searchRadius = 3;
%! nearestNeighbors = knnsearch (searchModel, queryPoints, "K", 2)
%! neighborsInRange = rangesearch (searchModel, queryPoints, searchRadius)

%!demo
%! ## Create an ExhaustiveSearcher with Euclidean distance
%! X = [1, 2; 3, 4; 5, 6];
%! obj = ExhaustiveSearcher (X);
%! ## Find the nearest neighbor to [2, 3]
%! Y = [2, 3];
%! [idx, D] = knnsearch (obj, Y);
%! disp ("Nearest neighbor index:"); disp (idx);
%! disp ("Distance:"); disp (D);
%! ## Find all points within radius 2
%! [idx, D] = rangesearch (obj, Y, 2);
%! disp ("Indices within radius:"); disp (idx);
%! disp ("Distances:"); disp (D);

%!demo
%! ## Create an ExhaustiveSearcher with Minkowski distance (P=1)
%! X = [0, 0; 1, 0; 0, 1];
%! obj = ExhaustiveSearcher (X, "Distance", "minkowski", "P", 1);
%! ## Find the 2 nearest neighbors to [0.5, 0.5]
%! Y = [0.5, 0.5];
%! [idx, D] = knnsearch (obj, Y, "K", 2);
%! disp ("Nearest neighbor indices:"); disp (idx);
%! disp ("Distances:"); disp (D);

%!demo
%! rng(42);
%! disp('Demonstrating ExhaustiveSearcher');
%!
%! n = 100;
%! mu1 = [0.3, 0.3];
%! mu2 = [0.7, 0.7];
%! sigma = 0.1;
%! X1 = mu1 + sigma * randn(n/2, 2);
%! X2 = mu2 + sigma * randn(n/2, 2);
%! X = [X1; X2];
%!
%! obj = ExhaustiveSearcher(X);
%!
%! Y = [0.3, 0.3; 0.7, 0.7; 0.5, 0.5];
%!
%! K = 5;
%! [idx, D] = knnsearch(obj, Y, "K", K);
%!
%! disp('For the first query point:');
%! disp(['Query point: ', num2str(Y(1,:))]);
%! disp('Indices of nearest neighbors:');
%! disp(idx(1,:));
%! disp('Distances:');
%! disp(D(1,:));
%!
%! figure;
%! scatter(X(:,1), X(:,2), 36, 'b', 'filled'); % Training points
%! hold on;
%! scatter(Y(:,1), Y(:,2), 36, 'r', 'filled'); % Query points
%! for i = 1:size(Y,1)
%!     query = Y(i,:);
%!     neighbors = X(idx(i,:), :);
%!     for j = 1:K
%!         plot([query(1), neighbors(j,1)], [query(2), neighbors(j,2)], 'k-');
%!     end
%! end
%! hold off;
%! title('K Nearest Neighbors with ExhaustiveSearcher');
%! xlabel('X1');
%! ylabel('X2');
%!
%! r = 0.15;
%! [idx, D] = rangesearch(obj, Y, r);
%!
%! disp('For the first query point in rangesearch:');
%! disp(['Query point: ', num2str(Y(1,:))]);
%! disp('Indices of points within radius:');
%! disp(idx{1});
%! disp('Distances:');
%! disp(D{1});
%!
%! figure;
%! scatter(X(:,1), X(:,2), 36, 'b', 'filled');
%! hold on;
%! scatter(Y(:,1), Y(:,2), 36, 'r', 'filled');
%! theta = linspace(0, 2*pi, 100);
%! for i = 1:size(Y,1)
%!     center = Y(i,:);
%!     x_circle = center(1) + r * cos(theta);
%!     y_circle = center(2) + r * sin(theta);
%!     plot(x_circle, y_circle, 'g-');
%!     % Highlight points within radius
%!     if ~isempty(idx{i})
%!         in_radius = X(idx{i}, :);
%!         scatter(in_radius(:,1), in_radius(:,2), 36, 'g', 'filled');
%!     end
%! end
%! hold off;
%! title('Points within Radius with ExhaustiveSearcher');
%! xlabel('X1');
%! ylabel('X2');

## Test Cases

%!test
%! ## Basic constructor with default Euclidean
%! X = [1, 2; 3, 4; 5, 6];
%! obj = ExhaustiveSearcher (X);
%! assert (obj.X, X)
%! assert (obj.Distance, "euclidean")
%! assert (isempty (obj.DistParameter))

%!test
%! ## Minkowski distance with custom P
%! X = [1, 2; 3, 4];
%! obj = ExhaustiveSearcher (X, "Distance", "minkowski", "P", 3);
%! assert (obj.Distance, "minkowski")
%! assert (obj.DistParameter, 3)

%!test
%! ## Seuclidean distance with custom Scale
%! X = [1, 2; 3, 4; 5, 6];
%! S = [1, 2];
%! obj = ExhaustiveSearcher (X, "Distance", "seuclidean", "Scale", S);
%! assert (obj.Distance, "seuclidean")
%! assert (obj.DistParameter, S)

%!test
%! ## Mahalanobis distance with custom Cov
%! X = [1, 2; 3, 4; 5, 6];
%! C = [1, 0; 0, 1];
%! obj = ExhaustiveSearcher (X, "Distance", "mahalanobis", "Cov", C);
%! assert (obj.Distance, "mahalanobis")
%! assert (obj.DistParameter, C)

%!test
%! ## knnsearch with Euclidean distance
%! X = [1, 2; 3, 4; 5, 6];
%! obj = ExhaustiveSearcher (X);
%! Y = [2, 3];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (idx, 1)
%! assert (D, sqrt(2), 1e-10)

%!test
%! ## knnsearch with Cityblock distance
%! X = [0, 0; 1, 1; 2, 2];
%! obj = ExhaustiveSearcher (X, "Distance", "cityblock");
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (idx, 1)
%! assert (D, 1, 1e-10)

%!test
%! ## knnsearch with Chebychev distance
%! X = [1, 1; 2, 3; 4, 2];
%! obj = ExhaustiveSearcher (X, "Distance", "chebychev");
%! Y = [2, 2];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 1)
%! assert (D, 1, 1e-10)

%!test
%! ## knnsearch with Cosine distance
%! X = [1, 0; 0, 1; 1, 1];
%! obj = ExhaustiveSearcher (X, "Distance", "cosine");
%! Y = [1, 0.5];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 3)
%! assert (D < 0.1, true)

%!test
%! ## knnsearch with Minkowski P=1 (Manhattan)
%! X = [0, 0; 1, 0; 0, 1];
%! obj = ExhaustiveSearcher (X, "Distance", "minkowski", "P", 1);
%! Y = [0.5, 0.5];
%! [idx, D] = knnsearch (obj, Y, "K", 2, "IncludeTies", true);
%! assert (iscell (idx))
%! assert (idx{1}, [1, 2, 3])
%! assert (D{1}, [1, 1, 1], 1e-10)

%!test
%! ## rangesearch with Seuclidean
%! X = [1, 1; 2, 2; 3, 3];
%! S = [1, 1];
%! obj = ExhaustiveSearcher (X, "Distance", "seuclidean", "Scale", S);
%! Y = [0, 0];
%! [idx, D] = rangesearch (obj, Y, 2);
%! assert (idx{1}, [1])
%! assert (D{1}, [sqrt(2)], 1e-10)

%!test
%! ## rangesearch with Mahalanobis
%! X = [1, 1; 2, 2; 3, 3];
%! C = [1, 0; 0, 1];
%! obj = ExhaustiveSearcher (X, "Distance", "mahalanobis", "Cov", C);
%! Y = [0, 0];
%! [idx, D] = rangesearch (obj, Y, 3, "SortIndices", false);
%! assert (idx{1}, [1, 2])
%! assert (D{1}, [sqrt(2), sqrt(8)], 1e-10)

%!test
%! ## rangesearch with Hamming distance
%! X = [0, 1; 1, 0; 1, 1];
%! obj = ExhaustiveSearcher (X, "Distance", "hamming");
%! Y = [0, 0];
%! [idx, D] = rangesearch (obj, Y, 0.5);
%! assert (idx{1}, [1, 2])
%! assert (D{1}, [0.5, 0.5], 1e-10)

%!test
%! ## Custom distance function
%! X = [1, 2; 3, 4];
%! custom_dist = @(x, y) sum(abs(x - y));
%! obj = ExhaustiveSearcher (X, "Distance", custom_dist);
%! Y = [2, 3];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 1)
%! assert (D, 2, 1e-10)

%!test
%! ## Euclidean with high-dimensional data
%! X = [1, 2, 3; 4, 5, 6; 7, 8, 9; 10, 11, 12];
%! obj = ExhaustiveSearcher (X);
%! Y = [5, 6, 7];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 2)
%! assert (D, sqrt(3), 1e-10)

%!test
%! ## Minkowski P=3 with scaled data
%! X = [0, 1; 2, 3; 4, 5] * 10;
%! obj = ExhaustiveSearcher (X, "Distance", "minkowski", "P", 3);
%! Y = [20, 30];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 2)
%! assert (D, 0, 1e-10)

%!test
%! ## Seuclidean with custom scales on diverse data
%! X = [1, 10; 2, 20; 3, 30];
%! S = [1, 5];
%! obj = ExhaustiveSearcher (X, "Distance", "seuclidean", "Scale", S);
%! Y = [1.5, 15];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 1)
%! assert (D, sqrt((0.5/1)^2 + (5/5)^2), 1e-10)

%!test
%! ## Mahalanobis with correlated data
%! X = [1, 1; 2, 1.5; 3, 2];
%! C = [1, 0.5; 0.5, 1];
%! obj = ExhaustiveSearcher (X, "Distance", "mahalanobis", "Cov", C);
%! Y = [2, 1.5];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 2)
%! assert (D, 0, 1e-10)

%!test
%! ## Cityblock with sparse data
%! X = [0, 0, 1; 1, 0, 0; 0, 1, 0];
%! obj = ExhaustiveSearcher (X, "Distance", "cityblock");
%! Y = [0, 0, 0];
%! [idx, D] = rangesearch (obj, Y, 1);
%! assert (idx{1}, [1, 2, 3])
%! assert (D{1}, [1, 1, 1], 1e-10)

%!test
%! ## Chebychev with extreme values
%! X = [0, 100; 50, 50; 100, 0];
%! obj = ExhaustiveSearcher (X, "Distance", "chebychev");
%! Y = [60, 60];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 2)
%! assert (D, 10, 1e-10)

%!test
%! ## Cosine with normalized data
%! X = [1, 0; 0, 1; 1/sqrt(2), 1/sqrt(2)];
%! obj = ExhaustiveSearcher (X, "Distance", "cosine");
%! Y = [1, 1];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 3)
%! assert (D < 0.1, true)

%!test
%! ## Correlation with time-series-like data
%! X = [1, 2, 3; 2, 4, 6; 1, 1, 1];
%! obj = ExhaustiveSearcher (X, "Distance", "correlation");
%! Y = [1.5, 3, 4.5];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 1)
%! assert (D < 0.1, true)

%!test
%! ## Spearman with ranked data
%! X = [1, 2, 3; 3, 2, 1; 2, 1, 3];
%! obj = ExhaustiveSearcher (X, "Distance", "spearman");
%! Y = [1, 2, 3];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 1)
%! assert (D, 0, 1e-10)

%!test
%! ## Jaccard with binary sparse data
%! X = [1, 0, 0; 0, 1, 0; 1, 1, 0];
%! obj = ExhaustiveSearcher (X, "Distance", "jaccard");
%! Y = [1, 0, 0];
%! [idx, D] = knnsearch (obj, Y);
%! assert (idx, 1)
%! assert (D, 0, 1e-10)

%!test
%! obj = ExhaustiveSearcher (ones(3,2));
%! assert (obj.X, ones(3,2))
%! assert (obj.Distance, "euclidean")
%! assert (isempty (obj.DistParameter))

%!test
%! obj = ExhaustiveSearcher (ones(3,2));
%! obj.Distance = "minkowski";
%! assert (obj.Distance, "minkowski")

%!test
%! obj = ExhaustiveSearcher (ones(3,2), "Distance", "minkowski");
%! obj.DistParameter = 3;
%! assert (obj.DistParameter, 3)

%!test
%! obj = ExhaustiveSearcher (ones(3,2), "Distance", "seuclidean");
%! obj.DistParameter = [1, 2];
%! assert (obj.DistParameter, [1, 2])

%!test
%! obj = ExhaustiveSearcher (ones(3,2), "Distance", "mahalanobis");
%! obj.DistParameter = eye(2);
%! assert (obj.DistParameter, eye(2))

## Test Input Validation

%!error<ExhaustiveSearcher: too few input arguments.> ...
%! ExhaustiveSearcher ()
%!error<ExhaustiveSearcher: Name-Value arguments must be in pairs.> ...
%! ExhaustiveSearcher (ones(3,2), "Distance")
%!error<ExhaustiveSearcher: X must be a finite numeric matrix.> ...
%! ExhaustiveSearcher ("abc")
%!error<ExhaustiveSearcher: X must be a finite numeric matrix.> ...
%! ExhaustiveSearcher ([1; Inf; 3])
%!error<ExhaustiveSearcher: invalid parameter name: 'foo'.> ...
%! ExhaustiveSearcher (ones(3,2), "foo", "bar")
%!error<ExhaustiveSearcher: unsupported distance metric 'invalid'.> ...
%! ExhaustiveSearcher (ones(3,2), "Distance", "invalid")
%!error<ExhaustiveSearcher: invalid distance function handle.> ...
%! ExhaustiveSearcher (ones(3,2), "Distance", @(x) x)
%!error<ExhaustiveSearcher: Distance must be a string or function handle.> ...
%! ExhaustiveSearcher (ones(3,2), "Distance", 1)
%!error<ExhaustiveSearcher: P must be a positive finite scalar.> ...
%! ExhaustiveSearcher (ones(3,2), "Distance", "minkowski", "P", -1)
%!error<ExhaustiveSearcher: Scale must be a nonnegative vector matching X columns.> ...
%! ExhaustiveSearcher (ones(3,2), "Distance", "seuclidean", "Scale", [-1, 1])
%!error<ExhaustiveSearcher: Cov must be a square matrix matching X columns.> ...
%! ExhaustiveSearcher (ones(3,2), "Distance", "mahalanobis", "Cov", ones(3,3))
%!error<ExhaustiveSearcher: Cov must be positive definite.> ...
%! ExhaustiveSearcher (ones(3,2), "Distance", "mahalanobis", "Cov", -eye(2))

%!error<ExhaustiveSearcher.knnsearch: too few input arguments.> ...
%! knnsearch (ExhaustiveSearcher (ones(3,2)))
%!error<ExhaustiveSearcher.knnsearch: Name-Value arguments must be in pairs.> ...
%! knnsearch (ExhaustiveSearcher (ones(3,2)), ones(3,2), "IncludeTies")
%!error<ExhaustiveSearcher.knnsearch: Y must be a finite numeric matrix.> ...
%! knnsearch (ExhaustiveSearcher (ones(3,2)), "abc")
%!error<ExhaustiveSearcher.knnsearch: number of columns in X and Y must match.> ...
%! knnsearch (ExhaustiveSearcher (ones(3,2)), ones(3,3))
%!error<ExhaustiveSearcher.knnsearch: K must be a positive integer.> ...
%! knnsearch (ExhaustiveSearcher (ones(3,2)), ones(3,2), "K", 0)
%!error<ExhaustiveSearcher.knnsearch: invalid parameter name: 'foo'.> ...
%! knnsearch (ExhaustiveSearcher (ones(3,2)), ones(3,2), "foo", "bar")
%!error<ExhaustiveSearcher.knnsearch: IncludeTies must be a logical scalar.> ...
%! knnsearch (ExhaustiveSearcher (ones(3,2)), ones(3,2), "IncludeTies", 1)

%!error<ExhaustiveSearcher.rangesearch: too few input arguments.> ...
%! rangesearch (ExhaustiveSearcher (ones(3,2)))
%!error<ExhaustiveSearcher.rangesearch: Name-Value arguments must be in pairs.> ...
%! rangesearch (ExhaustiveSearcher (ones(3,2)), ones(3,2), 1, "SortIndices")
%!error<ExhaustiveSearcher.rangesearch: Y must be a finite numeric matrix.> ...
%! rangesearch (ExhaustiveSearcher (ones(3,2)), "abc", 1)
%!error<ExhaustiveSearcher.rangesearch: number of columns in X and Y must match.> ...
%! rangesearch (ExhaustiveSearcher (ones(3,2)), ones(3,3), 1)
%!error<ExhaustiveSearcher.rangesearch: R must be a nonnegative finite scalar.> ...
%! rangesearch (ExhaustiveSearcher (ones(3,2)), ones(3,2), -1)
%!error<ExhaustiveSearcher.rangesearch: invalid parameter name: 'foo'.> ...
%! rangesearch (ExhaustiveSearcher (ones(3,2)), ones(3,2), 1, "foo", "bar")
%!error<ExhaustiveSearcher.rangesearch: SortIndices must be a logical scalar.> ...
%! rangesearch (ExhaustiveSearcher (ones(3,2)), ones(3,2), 1, "SortIndices", 1)

%!error<ExhaustiveSearcher.subsref: \(\) indexing not supported.> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj(1)
%!error<ExhaustiveSearcher.subsref: {} indexing not supported.> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj{1}
%!error<dynamic structure field names must be strings> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj.(1)
%!error<ExhaustiveSearcher.subsref: unrecognized property: 'invalid'> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj.invalid
%!error<ExhaustiveSearcher.subsasgn: \(\) indexing not supported.> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj(1) = 1
%!error<ExhaustiveSearcher.subsasgn: {} indexing not supported.> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj{1} = 1
%!error<ExhaustiveSearcher.subsasgn: chained subscripts not allowed.> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj.X.Y = 1
%!error<dynamic structure field names must be strings> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj.(1) = 1
%!error<ExhaustiveSearcher.subsasgn: X is read-only and cannot be modified.> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj.X = 1
%!error<ExhaustiveSearcher.subsasgn: unsupported distance metric 'invalid'.> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj.Distance = "invalid"
%!error<ExhaustiveSearcher.subsasgn: invalid distance function handle.> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj.Distance = @(x) x
%!error<ExhaustiveSearcher.subsasgn: custom distance function output invalid.> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj.Distance = @(x, y) [1; 1]
%!error<ExhaustiveSearcher.subsasgn: Distance must be a string or function handle.> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj.Distance = 1
%!error<ExhaustiveSearcher.subsasgn: DistParameter must be a positive finite scalar for minkowski.> ...
%! obj = ExhaustiveSearcher (ones(3,2), "Distance", "minkowski"); obj.DistParameter = -1
%!error<ExhaustiveSearcher.subsasgn: DistParameter must be a nonnegative vector matching X columns.> ...
%! obj = ExhaustiveSearcher (ones(3,2), "Distance", "seuclidean"); obj.DistParameter = [-1, 1]
%!error<ExhaustiveSearcher.subsasgn: DistParameter must be a square matrix matching X columns.> ...
%! obj = ExhaustiveSearcher (ones(3,2), "Distance", "mahalanobis"); obj.DistParameter = ones(3,3)
%!error<ExhaustiveSearcher.subsasgn: DistParameter must be positive definite for mahalanobis.> ...
%! obj = ExhaustiveSearcher (ones(3,2), "Distance", "mahalanobis"); obj.DistParameter = -eye(2)
%!error<ExhaustiveSearcher.subsasgn: DistParameter must be empty for this distance metric.> ...
%! obj = ExhaustiveSearcher (ones(3,2), "Distance", "euclidean"); obj.DistParameter = 1
%!error<ExhaustiveSearcher.subsasgn: unrecognized property: 'invalid'> ...
%! obj = ExhaustiveSearcher (ones(3,2)); obj.invalid = 1
