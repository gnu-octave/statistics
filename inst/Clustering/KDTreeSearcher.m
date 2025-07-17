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

classdef KDTreeSearcher
## -*- texinfo -*-
## @deftp {Class} KDTreeSearcher
##
## KD-tree nearest neighbor searcher class.
##
## The @code{KDTreeSearcher} class implements a KD-tree search algorithm for
## nearest neighbor queries.  It stores training data and supports various
## distance metrics along with their parameter values for performing a KD-tree
## search.  The KD-tree algorithm partitions the training data into a
## hierarchical tree structure and performs search operations by traversing the
## tree to reduce the number of distance computations.  It facilitates a nearest
## neighborsearch using @code{knnsearch} or a radius search using
## @code{rangesearch}.
##
## You can either use the @code{KDTreeSearcher} class constructor or the
## @code{createns} function to create an @qcode{KDTreeSearcher} object.
##
## @seealso{createns, ExhaustiveSearcher, hnswSearcher, knnsearch, rangesearch}
## @end deftp

  properties (SetAccess = private, Hidden)

    ## -*- texinfo -*-
    ## @deftp {Property} KDTree
    ##
    ## The KD-tree structure built from the training data.  This property is
    ## private and cannot be modified after object creation.
    ##
    ## @end deftp
    KDTree
  endproperties

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
    ## @deftp {Property} BucketSize
    ##
    ## The maximum number of data points in the leaf node of the KD-tree.
    ## Default is 50.
    ##
    ## @end deftp
    BucketSize = 50
  endproperties

  properties
    ## -*- texinfo -*-
    ## @deftp {Property} Distance
    ##
    ## Distance metric used for searches, specified as a character vector (e.g.,
    ## @qcode{"euclidean"}, @qcode{"minkowski"}).  Default is
    ## @qcode{"euclidean"}.  Supported metrics are @qcode{"euclidean"},
    ## @qcode{"cityblock"}, @qcode{"minkowski"}, and @qcode{"chebychev"}.
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
    ## @item Empty for other metrics (@qcode{"euclidean"}, @qcode{"cityblock"},
    ## @qcode{"chebychev"}). Attempting to set a non-empty value for these
    ## metrics will result in an error.
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
        fprintf ("\n  KDTreeSearcher with properties:\n\n");
        fprintf ("%+25s: %d\n", 'BucketSize', this.BucketSize);
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
          error ("KDTreeSearcher.subsref: () indexing not supported.");
        case '{}'
          error ("KDTreeSearcher.subsref: {} indexing not supported.");
        case '.'
          if (! ischar (s.subs))
            error (strcat ("KDTreeSearcher.subsref: property", ...
                           " name must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("KDTreeSearcher.subsref: unrecognized", ...
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
        error ("KDTreeSearcher.subsasgn: chained subscripts not allowed.");
      endif
      switch s.type
        case '()'
          error ("KDTreeSearcher.subsasgn: () indexing not supported.");
        case '{}'
          error ("KDTreeSearcher.subsasgn: {} indexing not supported.");
        case '.'
          if (! ischar (s.subs))
            error (strcat ("KDTreeSearcher.subsasgn: property", ...
                           " name must be a character vector."));
          endif
          switch (s.subs)
            case 'X'
              error (strcat ("KDTreeSearcher.subsasgn: 'X' is", ...
                             " read-only and cannot be modified."));
            case 'KDTree'
              error (strcat ("KDTreeSearcher.subsasgn: 'KDTree' is", ...
                             " read-only and cannot be modified."));
            case 'BucketSize'
              error (strcat ("KDTreeSearcher.subsasgn: 'BucketSize'", ...
                             " is read-only and cannot be modified."));
            case 'Distance'
              allowed_distances = {'euclidean', 'cityblock', 'minkowski', ...
                                   'chebychev'};
              if (ischar (val))
                if (! any (strcmpi (allowed_distances, val)))
                  error (strcat ("KDTreeSearcher.subsasgn:", ...
                                 " unsupported distance metric '%s'."), val);
                endif
                this.Distance = val;
              else
                error (strcat ("KDTreeSearcher.subsasgn: 'Distance'", ...
                               " must be a string."));
              endif
            case 'DistParameter'
              if (strcmpi (this.Distance, "minkowski"))
                if (! (isscalar (val) && isnumeric (val)
                                      && val > 0 && isfinite (val)))
                  error (strcat ("KDTreeSearcher.subsasgn: 'DistParameter'", ...
                                 " must be a positive finite scalar for", ...
                                 " Minkowski distance."));
                endif
                this.DistParameter = val;
              else
                if (! isempty (val))
                  error (strcat ("KDTreeSearcher.subsasgn: 'DistParameter'", ...
                                 " must be empty for this distance metric."));
                endif
                this.DistParameter = val;
              endif
            otherwise
              error ("KDTreeSearcher.subsasgn: unrecognized property: '%s'.",...
                     s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods

    ## -*- texinfo -*-
    ## @deftypefn  {KDTreeSearcher} {@var{obj} =} KDTreeSearcher (@var{X})
    ## @deftypefnx {KDTreeSearcher} {@var{obj} =} KDTreeSearcher (@var{X}, @var{name}, @var{value})
    ##
    ## Create a @qcode{KDTreeSearcher} object for nearest neighbor searches.
    ##
    ## @code{@var{obj} = KDTreeSearcher (@var{X})} constructs a
    ## @qcode{KDTreeSearcher} object with training data @var{X} using the
    ## default @qcode{"euclidean"} distance metric. @var{X} must be an
    ## @math{NxP} numeric matrix, where rows represent observations and columns
    ## represent features.
    ##
    ## @code{@var{obj} = KDTreeSearcher (@var{X}, @var{name}, @var{value})}
    ## allows customization through name-value pairs:
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"Distance"} @tab @tab Distance metric, specified as a
    ## character vector (@qcode{"euclidean"}, @qcode{"cityblock"},
    ## @qcode{"minkowski"}, @qcode{"chebychev"}).  Default is
    ## @qcode{"euclidean"}.
    ##
    ## @item @qcode{"P"} @tab @tab Minkowski distance exponent, a positive
    ## scalar.  Valid only when @qcode{"Distance"} is @qcode{"minkowski"}.
    ## Default is 2.
    ##
    ## @item @qcode{"BucketSize"} @tab @tab Maximum number of data points in the
    ## leaf node of the KD-tree, a positive integer.  Default is 50.
    ## @end multitable
    ##
    ## You can also create a @qcode{KDTreeSearcher} object using the
    ## @code{createns} function.
    ##
    ## @seealso{KDTreeSearcher, knnsearch, rangesearch, createns}
    ## @end deftypefn
    function obj = KDTreeSearcher (X, varargin)
      if (nargin < 1)
        error ("KDTreeSearcher: too few input arguments.");
      endif

      if (mod (numel (varargin), 2) != 0)
        error ("KDTreeSearcher: Name-Value arguments must be in pairs.");
      endif

      if (! (isnumeric (X) && ismatrix (X) && all (isfinite (X)(:))))
        error ("KDTreeSearcher: X must be a finite numeric matrix.");
      endif

      obj.X = X;

      ## Default values
      Distance = "euclidean";
      P = 2;
      BucketSize = 50;

      ## Parse optional parameters
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case "distance"
            Distance = varargin{2};
          case "p"
            P = varargin{2};
          case "bucketsize"
            BucketSize = varargin{2};
          otherwise
            error (strcat ("KDTreeSearcher: invalid parameter", ...
                           " name: '%s'."), varargin{1});
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Validate Distance
      allowed_distances = {"euclidean", "cityblock", "minkowski", "chebychev"};
      if (ischar (Distance))
        if (! any (strcmpi (allowed_distances, Distance)))
          error ("KDTreeSearcher: unsupported distance metric '%s'.", Distance);
        endif
        obj.Distance = Distance;
      else
        error ("KDTreeSearcher: Distance must be a string.");
      endif

      ## Set DistParameter
      if (strcmpi (obj.Distance, "minkowski"))
        if (! (isscalar (P) && isnumeric (P) && P > 0 && isfinite (P)))
          error ("KDTreeSearcher: P must be a positive finite scalar.");
        endif
        obj.DistParameter = P;
      else
        obj.DistParameter = [];
      endif

      ## Set BucketSize
      if (! (isscalar (BucketSize) && isnumeric (BucketSize)
                                   && BucketSize > 0
                                   && BucketSize == fix (BucketSize)))
        error ("KDTreeSearcher: BucketSize must be a positive integer.");
      endif
      obj.BucketSize = BucketSize;

      ## Build KDTree
      obj.KDTree = build_kdtree (1:size(X,1), 0, X, BucketSize);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {KDTreeSearcher} {[@var{idx}, @var{D}] =} knnsearch (@var{obj}, @var{Y}, @var{K})
    ## @deftypefnx {KDTreeSearcher} {[@var{idx}, @var{D}] =} knnsearch (@var{obj}, @var{Y}, @var{K}, @var{name}, @var{value})
    ##
    ## Find the @math{K} nearest neighbors in the training data to query points.
    ##
    ## @code{[@var{idx}, @var{D}] = knnsearch (@var{obj}, @var{Y}, @var{K})}
    ## returns the indices @var{idx} and distances @var{D} of the @math{K}
    ## nearest neighbors in @var{obj.X} to each point in @var{Y}, using the
    ## distance metric specified in @var{obj.Distance}.
    ##
    ## @itemize
    ## @item @var{obj} is a @qcode{KDTreeSearcher} object.
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
    ## @seealso{KDTreeSearcher, rangesearch}
    ## @end deftypefn
    function [idx, D] = knnsearch (obj, Y, K, varargin)
      if (nargin < 3)
        error ("KDTreeSearcher.knnsearch: too few input arguments.");
      endif

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("KDTreeSearcher.knnsearch:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      if (! (isnumeric (Y) && ismatrix (Y) && all (isfinite (Y)(:))))
        error ("KDTreeSearcher.knnsearch: Y must be a finite numeric matrix.");
      endif

      if (size (obj.X, 2) != size (Y, 2))
        error (strcat ("KDTreeSearcher.knnsearch:", ...
                       " number of columns in X and Y must match."));
      endif

      if (! (isscalar (K) && isnumeric (K) && K >= 1
                          && K == fix (K) && isfinite (K)))
        error ("KDTreeSearcher.knnsearch: K must be a positive integer.");
      endif

      ## Parse options
      IncludeTies = false;
      SortIndices = true;
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case "includeties"
            IncludeTies = varargin{2};
            if (! (islogical (IncludeTies) && isscalar (IncludeTies)))
              error (strcat ("KDTreeSearcher.knnsearch:", ...
                             " IncludeTies must be a logical scalar."));
            endif
          case "sortindices"
            SortIndices = varargin{2};
            if (! (islogical (SortIndices) && isscalar (SortIndices)))
              error (strcat ("KDTreeSearcher.knnsearch:", ...
                             " SortIndices must be a logical scalar."));
            endif
          otherwise
            error (strcat ("KDTreeSearcher.knnsearch: invalid", ...
                           " parameter name: '%s'."), varargin{1});
        endswitch
        varargin (1:2) = [];
      endwhile

      if (IncludeTies)
        idx = cell (rows (Y), 1);
        D = cell (rows (Y), 1);
        for i = 1:rows (Y)
          [temp_idx, temp_D] = search_kdtree (obj.KDTree, Y(i,:), K, obj.X, ...
                                              obj.Distance, obj.DistParameter, ...
                                              false);
          r = temp_D(end) + 1e-10; # Add small epsilon to capture ties
          [idx{i}, D{i}] = search_kdtree (obj.KDTree, Y(i,:), Inf, obj.X, ...
                                          obj.Distance, obj.DistParameter, ...
                                          true, r);
          if (SortIndices)
            [sorted_D, sort_idx] = sortrows ([D{i}(:), idx{i}(:)]);
            D{i} = sorted_D;
            idx{i} = idx{i}(sort_idx);
          endif
        endfor
      else
        idx = zeros (rows (Y), K);
        D = zeros (rows (Y), K);
        for i = 1:rows (Y)
          [temp_idx, temp_D] = search_kdtree (obj.KDTree, Y(i,:), K, obj.X, ...
                                              obj.Distance, obj.DistParameter, ...
                                              false);
          if (SortIndices)
            [sorted_D, sort_idx] = sortrows ([temp_D', temp_idx']);
            idx(i,:) = temp_idx(sort_idx);
            D(i,:) = sorted_D(:,1)';
          else
            idx(i,:) = temp_idx;
            D(i,:) = temp_D;
          endif
        endfor
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {KDTreeSearcher} {[@var{idx}, @var{D}] =} rangesearch (@var{obj}, @var{Y}, @var{r})
    ## @deftypefnx {KDTreeSearcher} {[@var{idx}, @var{D}] =} rangesearch (@var{obj}, @var{Y}, @var{r}, @var{name}, @var{value})
    ##
    ## Find all neighbors within a specified radius of query points.
    ##
    ## @code{[@var{idx}, @var{D}] = rangesearch (@var{obj}, @var{Y}, @var{r})}
    ## returns the indices @var{idx} and distances @var{D} of all points in
    ## @var{obj.X} within radius @var{r} of each point in @var{Y}, using the
    ## distance metric specified in @var{obj.Distance}.
    ##
    ## @itemize
    ## @item @var{obj} is a @qcode{KDTreeSearcher} object.
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
    ## @seealso{KDTreeSearcher, knnsearch}
    ## @end deftypefn
    function [idx, D] = rangesearch (obj, Y, r, varargin)
      if (nargin < 3)
        error ("KDTreeSearcher.rangesearch: too few input arguments.");
      endif

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("KDTreeSearcher.rangesearch:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      if (! (isnumeric (Y) && ismatrix (Y) && all (isfinite (Y)(:))))
        error (strcat ("KDTreeSearcher.rangesearch:", ...
                       " Y must be a finite numeric matrix."));
      endif

      if (size (obj.X, 2) != size (Y, 2))
        error (strcat ("KDTreeSearcher.rangesearch:", ...
                       " number of columns in X and Y must match."));
      endif

      if (! (isscalar (r) && isnumeric (r) && r >= 0 && isfinite (r)))
        error (strcat ("KDTreeSearcher.rangesearch:", ...
                       " R must be a nonnegative finite scalar."));
      endif

      ## Parse options
      SortIndices = true;
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case "sortindices"
            SortIndices = varargin{2};
            if (! (islogical (SortIndices) && isscalar (SortIndices)))
              error (strcat ("KDTreeSearcher.rangesearch:", ...
                             " SortIndices must be a logical scalar."));
            endif
          otherwise
            error (strcat ("KDTreeSearcher.rangesearch:", ...
                           " invalid parameter name: '%s'."), varargin{1});
        endswitch
        varargin (1:2) = [];
      endwhile

      idx = cell (rows (Y), 1);
      D = cell (rows (Y), 1);
      for i = 1:rows (Y)
        [idx{i}, D{i}] = search_kdtree (obj.KDTree, Y(i,:), Inf, obj.X, ...
                                        obj.Distance, obj.DistParameter, ...
                                        true, r);
        if (SortIndices)
          [sorted_D, sort_idx] = sortrows ([D{i}(:), idx{i}(:)]);
          D{i} = sorted_D;
          idx{i} = idx{i}(sort_idx);
        endif
      endfor
    endfunction

  endmethods

endclassdef

## Private functions:

## Function to Build KD-tree
function node = build_kdtree (indices, depth, X, bucket_size)
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
    left_node = build_kdtree (left_indices, depth + 1, X, bucket_size);
    right_node = build_kdtree (right_indices, depth + 1, X, bucket_size);
    node = struct ('axis', axis, 'split_value', split_value, ...
                   'left', left_node, 'right', right_node);
  endif
endfunction

## Function Search KD-tree
function [indices, distances] = search_kdtree (node, query, k, X, dist, ...
                                               distparam, is_range, r)
  if (nargin < 8)
    r = Inf;
  endif
  if (strcmpi (dist, "minkowski"))
    if (! (isscalar (distparam) && isnumeric (distparam) ...
                                && distparam > 0 && isfinite (distparam)))
      error (strcat("search_kdtree: distparam must be a positive finite", ...
                    " scalar for minkowski."));
    endif
  else
    if (! isempty (distparam))
      error (strcat("search_kdtree: distparam must be empty for", ...
                    " non-minkowski metrics."));
    endif
  endif
  indices = zeros(1, 0);
  distances = zeros(1, 0);
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
        indices = horzcat (indices, leaf_indices(mask));
        distances = horzcat (distances, dists(mask)');
      else
        indices = horzcat (indices, leaf_indices);
        distances = horzcat (distances, dists');
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

## Demo Examples

%!demo
%! ## Demo to verify implementation using fisheriris dataset
%! load fisheriris
%! numSamples = size (meas, 1);
%! queryIndices = [1, 23, 46, 63, 109];
%! dataIndices = ~ismember (1:numSamples, queryIndices);
%! queryPoints = meas(queryIndices, :);
%! dataPoints = meas(dataIndices, :);
%! searchRadius = 0.3;
%! kdTree = KDTreeSearcher (dataPoints, 'Distance', 'minkowski')
%! nearestNeighbors = knnsearch (kdTree, queryPoints, 2)
%! neighborsInRange = rangesearch (kdTree, queryPoints, searchRadius)

%!demo
%! ## Create a KDTreeSearcher with Euclidean distance
%! X = [1, 2; 3, 4; 5, 6];
%! obj = KDTreeSearcher (X);
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
%! ## Create a KDTreeSearcher with Minkowski distance (P=3)
%! X = [0, 0; 1, 0; 2, 0];
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 3);
%! ## Find the nearest neighbor to [1, 0]
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, 1);
%! disp ("Nearest neighbor index:");
%! disp (idx);
%! disp ("Distance:");
%! disp (D);

%!demo
%! rng(42);
%! disp('Demonstrating KDTreeSearcher');
%!
%! n = 100;
%! mu1 = [0.3, 0.3];
%! mu2 = [0.7, 0.7];
%! sigma = 0.1;
%! X1 = mu1 + sigma * randn (n / 2, 2);
%! X2 = mu2 + sigma * randn (n / 2, 2);
%! X = [X1; X2];
%!
%! obj = KDTreeSearcher(X);
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
%! title ('K Nearest Neighbors with KDTreeSearcher');
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
%! hold off
%! title ('Points within Radius with KDTreeSearcher');
%! xlabel ('X1');
%! ylabel ('X2');

## Test Cases

%!test
%! load fisheriris
%! X = meas;
%! obj = KDTreeSearcher (X);
%! Y = X(1:5,:);
%! [idx, D] = knnsearch (obj, Y, 3);
%! assert (idx, [[1, 18, 5]; [2, 35, 46]; [3, 48, 4]; [4, 48, 30]; [5, 38, 1]])
%! assert (D, [[0, 0.1000, 0.1414]; [0,  0.1414,  0.1414]; [0, 0.1414, 0.2449];
%!             [0, 0.1414, 0.1732]; [0, 0.1414, 0.1414]], 5e-5)

%!test
%! load fisheriris
%! X = meas;
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 3);
%! Y = X(10:15,:);
%! [idx, D] = knnsearch (obj, Y, 2);
%! assert (idx, [[10, 35]; [11, 49]; [12, 30]; [13, 2]; [14, 39]; [15, 34]])
%! assert (D, [[0, 0.1000]; [0, 0.1000]; [0, 0.2080]; [0, 0.1260]; [0, 0.2154];
%!             [0, 0.3503]], 5e-5)

%!test
%! load fisheriris
%! X = meas;
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! Y = X(20:25,:);
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (idx, [20; 21; 22; 23; 24; 25])
%! assert (D, [0; 0; 0; 0; 0; 0])

%!test
%! load fisheriris
%! X = meas;
%! obj = KDTreeSearcher (X, "Distance", "chebychev");
%! Y = X(30:35,:);
%! [idx, D] = knnsearch (obj, Y, 4);
%! assert (idx, [[30, 31, 4, 12]; [31, 30, 10, 35]; [32, 21, 37, 28];
%!               [33, 20, 34, 47]; [34, 16, 15, 33]; [35, 10, 2, 26]])
%! assert (D, [[0, 0.1000, 0.1000, 0.2000]; [0, 0.1000, 0.1000, 0.1000];
%!             [0, 0.2000, 0.2000, 0.2000]; [0, 0.3000, 0.3000, 0.3000];
%!             [0, 0.2000, 0.3000, 0.3000]; [0, 0.1000, 0.1000, 0.1000]], 5e-15)

%!test
%! load fisheriris
%! X = meas;
%! obj = KDTreeSearcher (X, "BucketSize", 20);
%! Y = X(40:45,:);
%! [idx, D] = knnsearch (obj, Y, 2);
%! assert (idx, [[40, 8]; [41, 18]; [42, 9]; [43, 39]; [44, 27]; [45, 47]])
%! assert (D, [[0, 0.1000]; [0, 0.1414]; [0, 0.6245]; [0, 0.2000]; [0, 0.2236];
%!             [0, 0.3606]], 4.7e-5)

%!test
%! load fisheriris
%! X = meas;
%! obj = KDTreeSearcher (X);
%! Y = X(50:55,:);
%! [idx, D] = knnsearch (obj, Y, 3, "IncludeTies", true);
%! assert (idx, {[50, 8, 40]; [51, 53, 87]; [52, 57, 76]; [53, 51, 87]; [54, ...
%!                90, 81]; [55, 59, 76]})
%! assert (D, {[0, 0.1414, 0.1732]; [0, 0.2646, 0.3317]; [0, 0.2646, 0.3162];
%!             [0, 0.2646, 0.2828]; [0, 0.2000, 0.3000]; [0, 0.2449, 0.3162]}, 5e-5)

%!test
%! load fisheriris
%! X = meas;
%! obj = KDTreeSearcher (X);
%! Y = X(60:65,:);
%! [idx, D] = rangesearch (obj, Y, 0.4);
%! assert (idx, {[60, 90]; [61, 94]; [62, 97, 79, 96, 100, 89, 98, 72]; [63];
%!               [64, 92, 74, 79]; [65]})
%! assert (D, {[0, 0.3873]; [0, 0.3606];
%!             [0, 0.3000, 0.3317, 0.3606, 0.3606, 0.3742, 0.3873, 0.4000]; [0];
%!             [0, 0.1414, 0.2236, 0.2449]; [0]}, 5e-5)

%!test
%! load fisheriris
%! X = meas;
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! Y = X(70:72,:);
%! [idx, D] = rangesearch (obj, Y, 1.0);
%! assert (idx, {[70, 81, 90, 82, 83, 93, 54, 68, 95, 80, 91, 100, 60, 65, 89, 63]; [71, 139, ...
%!                128, 150, 127, 57, 86, 64, 79, 92, 124]; [72, 100, 98, 83, 93, 97, 75, 68, ...
%!                62, 89, 95, 74, 56, 90, 79, 92, 96, 64, 63, 65]})
%! assert (D, {[0, 0.3000, 0.4000, 0.5000, 0.5000, 0.5000, 0.6000, 0.7000, 0.7000, ...
%!              0.7000, 0.8000, 0.8000, 0.9000, 0.9000, 0.9000, 0.9000];
%!             [0, 0.3000, 0.5000, 0.5000, 0.7000, 0.8000, 0.8000, 1.0000, 1.0000, ...
%!              1.0000, 1];
%!             [0, 0.5000, 0.5000, 0.6000, 0.6000, 0.7000, 0.7000, 0.8000, 0.8000, ...
%!              0.8000, 0.8000, 0.8000, 0.9000, 0.9000, 0.9000, 0.9000, 0.9000, ...
%!              0.9000, 1.0000, 1]}, 5e-5)

%!test
%! load fisheriris
%! X = meas;
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 3);
%! Y = X(80:85,:);
%! [idx, D] = rangesearch (obj, Y, 0.8);
%! assert (idx, {[80, 82, 81, 65, 70, 83, 93, 90, 54, 63, 68, 72, 100, 60, 89, 99, 95, 94, 97, ...
%!                96]; [81, 82, 70, 54, 90, 93, 80, 83, 60, 68, 95, 100, 65, 63, 97, 61, 91, ...
%!                94, 89, 96, 72, 58, 62, 56]; [82, 81, 70, 80, 54, 90, 93, 83, 68, 60, 65, 63, ...
%!                100, 95, 94, 61, 58, 97, 89, 72, 96, 91, 99]; [83, 93, 100, 68, 70, 72, 95, ...
%!                90, 97, 65, 89, 96, 81, 82, 80, 62, 54, 98, 63, 91, 56, 60, 79, 67, 75, 88, 85, ...
%!                92, 69]; [84, 134, 102, 143, 150, 124, 128, 73, 127, 139, 147, 64, 112, ...
%!                114, 120, 74, 135, 122, 92, 71, 104, 138, 148, 117, 79, 55, 56, 57, 67, 111, ...
%!                129, 69, 78, 59, 52, 133, 85, 88, 87]; [85, 67, 56, 97, 95, 89, 96, 91, 100, ...
%!                62, 71, 122, 79, 60, 107, 90, 139, 93, 68, 86, 83, 92, 64, 150, 102, 143, 74, ...
%!                114, 70, 128, 84, 54, 72]})
%! assert (D, {[0, 0.2884, 0.3530, 0.3826, 0.4062, 0.4198, 0.5117, 0.5440, 0.5718, ...
%!              0.6000, 0.6018, 0.6073, 0.6308, 0.6333, 0.6753, 0.7000, 0.7192, ...
%!              0.7230, 0.7350, 0.7459]; [0, 0.1260, 0.1442, 0.2571, 0.2571, 0.3530, ...
%!              0.3530, 0.3826, 0.4344, 0.4344, 0.4642, 0.4747, 0.5217, 0.5217, ...
%!              0.5896, 0.6009, 0.6082, 0.6316, 0.6316, 0.6611, 0.6664, 0.6993, ...
%!              0.7417, 0.7507]; [0, 0.1260, 0.2224, 0.2884, 0.3803, 0.3803, 0.4121, ...
%!              0.4121, 0.4905, 0.5013, 0.5360, 0.5429, 0.5463, 0.5646, 0.5749, ...
%!              0.5819, 0.6542, 0.6581, 0.6753, 0.6938, 0.7094, 0.7107, 0.7423]; [0, ...
%!              0.1260, 0.2224, 0.2520, 0.2571, 0.3107, 0.3302, 0.3332, 0.3332, ...
%!              0.3530, 0.3530, 0.3803, 0.3826, 0.4121, 0.4198, 0.4344, 0.4531, ...
%!              0.5155, 0.5217, 0.5348, 0.6028, 0.6073, 0.6374, 0.6527, 0.6611, ...
%!              0.6804, 0.6938, 0.7399, 0.7560]; [0, 0.3072, 0.3271, 0.3271, 0.3302, ...
%!              0.3503, 0.3530, 0.3530, 0.3530, 0.3958, 0.3979, 0.4327, 0.4626, ...
%!              0.4642, 0.5027, 0.5066, 0.5130, 0.5155, 0.5440, 0.5440, 0.5518, ...
%!              0.5848, 0.6009, 0.6073, 0.6082, 0.6316, 0.6471, 0.6746, 0.6753, ...
%!              0.6797, 0.6804, 0.7047, 0.7192, 0.7218, 0.7405, 0.7405, 0.7719, ...
%!              0.7725, 0.7786]; [0, 0.2000, 0.3503, 0.3979, 0.4121, 0.4309, 0.4327, ...
%!              0.4531, 0.4747, 0.5337, 0.5718, 0.5896, 0.6009, 0.6316, 0.6366, ...
%!              0.6374, 0.6463, 0.6542, 0.6542, 0.6550, 0.6938, 0.7014, 0.7067, ...
%!              0.7166, 0.7186, 0.7186, 0.7281, 0.7380, 0.7447, 0.7571, 0.7719, ...
%!              0.7813, 0.7851]}, 5e-5)

%!test
%! ## Constructor with single-point dataset
%! X = [0, 0];
%! obj = KDTreeSearcher (X);
%! assert (obj.X, X);
%! assert (obj.Distance, "euclidean");
%! assert (isempty (obj.DistParameter));
%! assert (obj.BucketSize, 50);

%!test
%! ## Constructor with duplicate points
%! X = [0, 0; 0, 0; 1, 0];
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! assert (obj.X, X);
%! assert (obj.Distance, "cityblock");

%!test
%! ## Constructor with 3D data
%! X = [0, 0, 0; 1, 0, 0; 0, 1, 0];
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 3);
%! assert (obj.X, X);
%! assert (obj.DistParameter, 3);

%!test
%! ## knnsearch with grid, K = 1
%! X = [0, 0; 0, 1; 1, 0; 1, 1];
%! obj = KDTreeSearcher (X, "Distance", "euclidean");
%! Y = [0.5, 0.5];
%! [idx, D] = knnsearch (obj, Y, 1);
%! D_true = pdist2 (X, Y, "euclidean");
%! assert (D, min (D_true), 1e-10);
%! assert (any (idx == find (D_true == min (D_true))));

%!test
%! ## knnsearch with IncludeTies, all points equidistant
%! X = [0, 0; 0, 1; 1, 0; 1, 1];
%! obj = KDTreeSearcher (X);
%! Y = [0.5, 0.5];
%! [idx, D] = knnsearch (obj, Y, 1, "IncludeTies", true);
%! D_true = pdist2 (X, Y, "euclidean");
%! expected_idx = find (D_true == min (D_true));
%! assert (sort (idx{1}(:)), sort (expected_idx));
%! assert (D{1}(:)', repmat (min (D_true), 1, 4), 1e-10);

%!test
%! ## rangesearch with line dataset
%! X = [0, 0; 1, 0; 2, 0; 3, 0];
%! obj = KDTreeSearcher (X);
%! Y = [1.5, 0];
%! r = 1;
%! [idx, D] = rangesearch (obj, Y, r);
%! D_true = pdist2 (X, Y, "euclidean");
%! expected_idx = find (D_true <= r);
%! assert (sort (idx{1}(:)), sort (expected_idx));
%! assert (D{1}', sort (D_true(expected_idx)), 1e-10);

%!test
%! ## knnsearch with duplicates
%! X = [0, 0; 0, 0; 1, 0];
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! Y = [0, 0];
%! [idx, D] = knnsearch (obj, Y, 1, "IncludeTies", true);
%! assert (sort (idx{1}(:))', [1, 2]);
%! assert (D{1}, [0, 0], 1e-10);

%!test
%! ## rangesearch with 3D data
%! X = [0, 0, 0; 1, 0, 0; 0, 1, 0];
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! Y = [0, 0, 0];
%! r = 1;
%! [idx, D] = rangesearch (obj, Y, r);
%! assert (sort (idx{1}(:))', [1, 2, 3]);
%! assert (D{1}, [0, 1, 1], 1e-10);

%!test
%! ## knnsearch with P = 2 (Euclidean equivalent)
%! X = [0, 0; 1, 1];
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 2);
%! Y = [0, 1];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (idx, 1);
%! assert (D, 1, 1e-10);

%!test
%! ## rangesearch with P = 3
%! X = [0, 0; 1, 0; 0, 1];
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 3);
%! Y = [0.5, 0.5];
%! r = 0.8;
%! [idx, D] = rangesearch (obj, Y, r);
%! D_true = pdist2 (X, Y, "minkowski", 3);
%! expected_idx = find (D_true <= r);
%! assert (sort (idx{1}(:)), sort (expected_idx));
%! assert (D{1}', sort (D_true(expected_idx)), 1e-10);

%!test
%! ## knnsearch with P = 4, random data
%! X = rand (5, 2);
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 4);
%! Y = rand (1, 2);
%! [idx, D] = knnsearch (obj, Y, 3);
%! D_true = pdist2 (X, Y, "minkowski", 4);
%! [sorted_D, sort_idx] = sort (D_true);
%! assert (idx', sort_idx(1:3));
%! assert (D', sorted_D(1:3), 1e-10);

%!test
%! ## knnsearch with all same points
%! X = [1, 1; 1, 1; 1, 1];
%! obj = KDTreeSearcher (X, "Distance", "chebychev");
%! Y = [1, 1];
%! [idx, D] = knnsearch (obj, Y, 1, "IncludeTies", true);
%! assert (sort (idx{1}(:))', [1, 2, 3]);
%! assert (D{1}, [0, 0, 0], 1e-10);

%!test
%! ## rangesearch with grid
%! X = [0, 0; 0, 1; 1, 0; 1, 1];
%! obj = KDTreeSearcher (X, "Distance", "chebychev");
%! Y = [0.5, 0.5];
%! r = 0.5;
%! [idx, D] = rangesearch (obj, Y, r);
%! D_true = pdist2 (X, Y, "chebychev");
%! expected_idx = find (D_true <= r);
%! assert (sort (idx{1}(:)), sort (expected_idx));
%! assert (D{1}', D_true(expected_idx), 1e-10);

%!test
%! ## Changing Distance and verifying search
%! X = [0,0; 1,0];
%! obj = KDTreeSearcher(X, "Distance", "euclidean");
%! Y = [0,1];
%! [idx, D] = knnsearch(obj, Y, 1);
%! assert(D, 1, 1e-10);
%! obj.Distance = "chebychev";
%! [idx, D] = knnsearch(obj, Y, 1);
%! assert(D, 1, 1e-10);

%!test
%! ## Changing DistParameter for minkowski
%! X = [0,0; 1,0];
%! obj = KDTreeSearcher(X, "Distance", "minkowski", "P", 1);
%! Y = [0,1];
%! [idx, D] = knnsearch(obj, Y, 1);
%! assert(D, 1, 1e-10);
%! obj.DistParameter = 3;
%! [idx, D] = knnsearch(obj, Y, 1);
%! assert(D, 1, 1e-10);

%!test
%! ## Different BucketSize values
%! X = rand(20,2);
%! obj1 = KDTreeSearcher(X, "BucketSize", 5);
%! obj2 = KDTreeSearcher(X, "BucketSize", 15);
%! Y = rand(1,2);
%! [idx1, D1] = knnsearch(obj1, Y, 3);
%! [idx2, D2] = knnsearch(obj2, Y, 3);
%! assert(idx1, idx2);
%! assert(D1, D2, 1e-10);

%!test
%! ## Basic constructor with default Euclidean
%! X = [1, 2; 3, 4; 5, 6];
%! obj = KDTreeSearcher (X);
%! assert (obj.X, X);
%! assert (obj.Distance, "euclidean");
%! assert (isempty (obj.DistParameter));
%! assert (obj.BucketSize, 50);

%!test
%! ## Minkowski distance with custom P
%! X = [0, 0; 1, 1; 2, 2];
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 3);
%! assert (obj.Distance, "minkowski");
%! assert (obj.DistParameter, 3);

%!test
%! ## Cityblock distance
%! X = [0, 0; 1, 0; 0, 1];
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! assert (obj.Distance, "cityblock");
%! assert (isempty (obj.DistParameter));

%!test
%! ## Chebychev distance
%! X = [1, 1; 2, 3; 4, 2];
%! obj = KDTreeSearcher (X, "Distance", "chebychev");
%! assert (obj.Distance, "chebychev");
%! assert (isempty (obj.DistParameter));

%!test
%! ## knnsearch with Euclidean distance
%! X = [1, 2; 3, 4; 5, 6];
%! obj = KDTreeSearcher (X);
%! Y = [2, 3];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (idx, 1);
%! assert (D, sqrt(2), 1e-10);

%!test
%! ## knnsearch with Cityblock distance
%! X = [0, 0; 1, 1; 2, 2];
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (ismember (idx, [1, 2]));
%! assert (D, 1, 1e-10);

%!test
%! ## knnsearch with Chebychev distance
%! X = [1, 1; 2, 3; 4, 2];
%! obj = KDTreeSearcher (X, "Distance", "chebychev");
%! Y = [2, 2];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (ismember (idx, [1, 2]));
%! assert (D, 1, 1e-10);

%!test
%! ## knnsearch with Minkowski P=3
%! X = [0, 0; 1, 0; 2, 0];
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 3);
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (idx, 2);
%! assert (D, 0, 1e-10);

%!test
%! ## knnsearch with IncludeTies
%! X = [0, 0; 1, 0; 0, 1];
%! obj = KDTreeSearcher (X);
%! Y = [0.5, 0];
%! [idx, D] = knnsearch (obj, Y, 1, "IncludeTies", true);
%! assert (iscell (idx));
%! assert (sort (idx{1}(:))', [1, 2]);
%! assert (sort (D{1}(:)), [0.5; 0.5], 1e-10);

%!test
%! ## rangesearch with Euclidean
%! X = [1, 1; 2, 2; 3, 3];
%! obj = KDTreeSearcher (X);
%! Y = [0, 0];
%! [idx, D] = rangesearch (obj, Y, 2);
%! assert (idx{1}, [1]);
%! assert (D{1}, [sqrt(2)], 1e-10);

%!test
%! ## rangesearch with Cityblock
%! X = [0, 0; 1, 1; 2, 2];
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! Y = [0, 0];
%! [idx, D] = rangesearch (obj, Y, 1);
%! assert (idx{1}, [1]);
%! assert (D{1}, [0], 1e-10);

%!test
%! ## rangesearch with Chebychev
%! X = [1, 1; 2, 3; 4, 2];
%! obj = KDTreeSearcher (X, "Distance", "chebychev");
%! Y = [2, 2];
%! [idx, D] = rangesearch (obj, Y, 1);
%! assert (sort (idx{1}(:))', [1, 2]);
%! assert (sort (D{1}(:))', [1, 1], 1e-10);

%!test
%! ## rangesearch with Minkowski P=3
%! X = [0, 0; 1, 0; 2, 0];
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 3);
%! Y = [1, 0];
%! [idx, D] = rangesearch (obj, Y, 1);
%! assert (sort (idx{1}(:))', [1, 2, 3]);
%! assert (sort (D{1}(:))', [0, 1, 1], 1e-10);

%!test
%! ## Diverse dataset with Euclidean
%! X = [0, 10; 5, 5; 10, 0];
%! obj = KDTreeSearcher (X);
%! Y = [5, 5];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (idx, 2);
%! assert (D, 0, 1e-10);

%!test
%! ## High-dimensional data with Cityblock
%! X = [1, 2, 3; 4, 5, 6; 7, 8, 9];
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! Y = [4, 5, 6];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (idx, 2);
%! assert (D, 0, 1e-10);

## Test Input Validation

%!error<KDTreeSearcher: too few input arguments.> ...
%! KDTreeSearcher ()
%!error<KDTreeSearcher: Name-Value arguments must be in pairs.> ...
%! KDTreeSearcher (ones(3,2), "Distance")
%!error<KDTreeSearcher: X must be a finite numeric matrix.> ...
%! KDTreeSearcher ("abc")
%!error<KDTreeSearcher: X must be a finite numeric matrix.> ...
%! KDTreeSearcher ([1; Inf; 3])
%!error<KDTreeSearcher: invalid parameter name: 'foo'.> ...
%! KDTreeSearcher (ones(3,2), "foo", "bar")
%!error<KDTreeSearcher: unsupported distance metric 'invalid'.> ...
%! KDTreeSearcher (ones(3,2), "Distance", "invalid")
%!error<KDTreeSearcher: Distance must be a string.> ...
%! KDTreeSearcher (ones(3,2), "Distance", 1)
%!error<KDTreeSearcher: P must be a positive finite scalar.> ...
%! KDTreeSearcher (ones(3,2), "Distance", "minkowski", "P", -1)
%!error<KDTreeSearcher: BucketSize must be a positive integer.> ...
%! KDTreeSearcher (ones(3,2), "BucketSize", 0)
%!error<KDTreeSearcher: BucketSize must be a positive integer.> ...
%! KDTreeSearcher(ones(3,2), "BucketSize", -1)

%!error<KDTreeSearcher.knnsearch: too few input arguments.> ...
%! knnsearch (KDTreeSearcher (ones(3,2)))
%!error<KDTreeSearcher.knnsearch: Name-Value arguments must be in pairs.> ...
%! knnsearch (KDTreeSearcher (ones(3,2)), ones(3,2), 1, "IncludeTies")
%!error<KDTreeSearcher.knnsearch: Y must be a finite numeric matrix.> ...
%! knnsearch (KDTreeSearcher (ones(3,2)), "abc", 1)
%!error<KDTreeSearcher.knnsearch: number of columns in X and Y must match.> ...
%! knnsearch (KDTreeSearcher (ones(3,2)), ones(3,3), 1)
%!error<KDTreeSearcher.knnsearch: K must be a positive integer.> ...
%! knnsearch (KDTreeSearcher (ones(3,2)), ones(3,2), 0)
%!error<KDTreeSearcher.knnsearch: K must be a positive integer.> ...
%! obj = KDTreeSearcher(ones(3,2)); knnsearch(obj, ones(1,2), Inf)
%!error<KDTreeSearcher.knnsearch: invalid parameter name: 'foo'.> ...
%! knnsearch (KDTreeSearcher (ones(3,2)), ones(3,2), 1, "foo", "bar")
%!error<KDTreeSearcher.knnsearch: IncludeTies must be a logical scalar.> ...
%! knnsearch (KDTreeSearcher (ones(3,2)), ones(3,2), 1, "IncludeTies", 1)
%!error<KDTreeSearcher.knnsearch: SortIndices must be a logical scalar.> ...
%! knnsearch (KDTreeSearcher (ones(3,2)), ones(3,2), 1, "SortIndices", 1)

%!error<KDTreeSearcher.rangesearch: too few input arguments.> ...
%! rangesearch (KDTreeSearcher (ones(3,2)))
%!error<KDTreeSearcher.rangesearch: Name-Value arguments must be in pairs.> ...
%! rangesearch (KDTreeSearcher (ones(3,2)), ones(3,2), 1, "SortIndices")
%!error<KDTreeSearcher.rangesearch: Y must be a finite numeric matrix.> ...
%! rangesearch (KDTreeSearcher (ones(3,2)), "abc", 1)
%!error<KDTreeSearcher.rangesearch: number of columns in X and Y must match.> ...
%! rangesearch (KDTreeSearcher (ones(3,2)), ones(3,3), 1)
%!error<KDTreeSearcher.rangesearch: R must be a nonnegative finite scalar.> ...
%! rangesearch (KDTreeSearcher (ones(3,2)), ones(3,2), -1)
%!error<KDTreeSearcher.rangesearch: R must be a nonnegative finite scalar.> ...
%! obj = KDTreeSearcher(ones(3,2)); rangesearch(obj, ones(1,2), Inf)
%!error<KDTreeSearcher.rangesearch: invalid parameter name: 'foo'.> ...
%! rangesearch (KDTreeSearcher (ones(3,2)), ones(3,2), 1, "foo", "bar")
%!error<KDTreeSearcher.rangesearch: SortIndices must be a logical scalar.> ...
%! rangesearch (KDTreeSearcher (ones(3,2)), ones(3,2), 1, "SortIndices", 1)

%!error<KDTreeSearcher.subsref: \(\) indexing not supported.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj(1)
%!error<KDTreeSearcher.subsref: {} indexing not supported.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj{1}
%!error<KDTreeSearcher.subsref: unrecognized property: 'invalid'.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.invalid

%!error<KDTreeSearcher.subsasgn: \(\) indexing not supported.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj(1) = 1
%!error<KDTreeSearcher.subsasgn: {} indexing not supported.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj{1} = 1
%!error<KDTreeSearcher.subsasgn: chained subscripts not allowed.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.X.Y = 1
%!error<KDTreeSearcher.subsasgn: 'X' is read-only and cannot be modified.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.X = 1
%!error<KDTreeSearcher.subsasgn: 'KDTree' is read-only and cannot be modified.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.KDTree = 1
%!error<KDTreeSearcher.subsasgn: unsupported distance metric 'invalid'.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.Distance = "invalid"
%!error<KDTreeSearcher.subsasgn: 'Distance' must be a string.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.Distance = 1
%!error<KDTreeSearcher.subsasgn: 'DistParameter' must be a positive finite scalar for Minkowski distance.> ...
%! obj = KDTreeSearcher (ones(3,2), "Distance", "minkowski"); obj.DistParameter = -1
%!error<KDTreeSearcher.subsasgn: 'DistParameter' must be empty for this distance metric.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.DistParameter = 1
%!error<KDTreeSearcher.subsasgn: 'BucketSize' is read-only and cannot be modified.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.BucketSize = 0
%!error<KDTreeSearcher.subsasgn: 'BucketSize' is read-only and cannot be modified.> ...
%! obj = KDTreeSearcher(ones(3,2)); obj.BucketSize = -1
%!error<KDTreeSearcher.subsasgn: 'BucketSize' is read-only and cannot be modified.> ...
%! obj = KDTreeSearcher(ones(3,2)); obj.BucketSize = 1.5
%!error<KDTreeSearcher.subsasgn: unrecognized property: 'invalid'.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.invalid = 1
