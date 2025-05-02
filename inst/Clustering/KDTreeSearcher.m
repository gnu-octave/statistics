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
## @seealso{createns, ExhaustiveSearcher, knnsearch, rangesearch, pdist2}
## @end deftp

classdef KDTreeSearcher

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
    ## @deftp {Property} KDTree
    ##
    ## The KD-tree structure built from the training data.  This property is
    ## private and cannot be modified after object creation.
    ##
    ## @end deftp
    KDTree
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
    ## @qcode{"chebychev"}).
    ## @end itemize
    ##
    ## @end deftp
    DistParameter = []

    ## -*- texinfo -*-
    ## @deftp {Property} BucketSize
    ##
    ## The maximum number of data points in the leaf node of the KD-tree.
    ## Default is 50.
    ##
    ## @end deftp
    BucketSize = 50
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
        fprintf ("\n  KDTreeSearcher\n\n");
        fprintf ("%+25s: [%dx%d %s]\n", 'X', size (this.X), class (this.X));
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
        fprintf ("%+25s: %d\n", 'BucketSize', this.BucketSize);
        fprintf ("%+25s: [KD-tree structure]\n", 'KDTree');
      else
        sz = size (this);
        fprintf ("\n  %s KDTreeSearcher array\n\n", mat2str (sz));
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
              error (strcat ("KDTreeSearcher.subsasgn: X is", ...
                             " read-only and cannot be modified."));
            case 'KDTree'
              error (strcat ("KDTreeSearcher.subsasgn: KDTree is", ...
                             " read-only and cannot be modified."));
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
                error (strcat ("KDTreeSearcher.subsasgn: Distance", ...
                               " must be a string."));
              endif
            case 'DistParameter'
              if (strcmpi (this.Distance, "minkowski"))
                if (! (isscalar (val) && isnumeric (val)
                                      && val > 0 && isfinite (val)))
                  error (strcat ("KDTreeSearcher.subsasgn:", ...
                                 " DistParameter must be a positive", ...
                                 " finite scalar for minkowski."));
                endif
                this.DistParameter = val;
              else
                if (! isempty (val))
                  error (strcat ("KDTreeSearcher.subsasgn: DistParameter", ...
                                 " must be empty for this distance metric."));
                endif
                this.DistParameter = val;
              endif
            case 'BucketSize'
              if (! (isscalar (val) && isnumeric (val)
                                    && val > 0 && val == fix (val)))
                error (strcat ("KDTreeSearcher.subsasgn: BucketSize", ...
                               " must be a positive integer."));
              endif
              this.BucketSize = val;
            otherwise
              error (strcat ("KDTreeSearcher.subsasgn:", ...
                             " unrecognized property: '%s'."), s.subs);
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
      obj.KDTree = KDTreeSearcher.buildkdtree (X, BucketSize);
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

      ## Perform search
      if (! IncludeTies)
        idx = zeros (rows (Y), K);
        D = zeros (rows (Y), K);
        for i = 1:rows (Y)
          NN = KDTreeSearcher.findkdtree (obj.KDTree, Y(i,:), K, ...
                                           obj.Distance, obj.DistParameter);
          D_temp = pdist2 (obj.X(NN,:), Y(i,:), obj.Distance, ...
                           obj.DistParameter);
          [sorted_D, sort_idx] = sort (D_temp);
          NN_sorted = NN(sort_idx);
          if (SortIndices)
            idx(i,:) = NN_sorted(1:K);
            D(i,:) = sorted_D(1:K);
          else
            idx(i,:) = NN(1:K);
            D(i,:) = D_temp(1:K);
          endif
        endfor
      else
        idx = cell (rows (Y), 1);
        D = cell (rows (Y), 1);
        for i = 1:rows (Y)
          NN = KDTreeSearcher.findkdtree (obj.KDTree, Y(i,:), K, ...
                                           obj.Distance, obj.DistParameter);
          D_temp = pdist2 (obj.X(NN,:), Y(i,:), obj.Distance, ...
                           obj.DistParameter);
          [sorted_D, sort_idx] = sort (D_temp);
          if (K <= length (sorted_D))
            kth_dist = sorted_D(K);
          else
            kth_dist = sorted_D(end);
          endif
          [idx_temp, D_temp] = rangesearch (obj, Y(i,:), kth_dist, ...
                                            "SortIndices", SortIndices);
          idx{i} = idx_temp{1};
          D{i} = D_temp{1};
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

      ## Perform search
      k = rows (obj.X);  # Get all points
      idx = cell (rows (Y), 1);
      D = cell (rows (Y), 1);
      for i = 1:rows (Y)
        NN = KDTreeSearcher.findkdtree (obj.KDTree, Y(i,:), k, ...
                                         obj.Distance, obj.DistParameter);
        D_temp = pdist2 (obj.X(NN,:), Y(i,:), obj.Distance, ...
                         obj.DistParameter);
        within_r = find (D_temp <= r);
        if (SortIndices)
          [sorted_D, sort_idx] = sort (D_temp(within_r));
          idx{i} = NN(within_r(sort_idx));
          D{i} = sorted_D;
        else
          idx{i} = NN(within_r);
          D{i} = D_temp(within_r);
        endif
      endfor
    endfunction

  endmethods

  methods (Static, Access = private)

    function ret = buildkdtree (X, BS)
      [val, r] = sort (X(:,1));
      ret = struct ("data", X, "root", ...
                     KDTreeSearcher.buildkdtree_recur (X, r, 1, BS));
    endfunction

    function ret = buildkdtree_recur (X, r, d, BS)
      count = length (r);
      dimen = size (X, 2);
      if (count == 1)
        ret = struct ("point", r(1), "dimen", d);
      else
        mid = ceil (count / 2);
        ret = struct ("point", r(mid), "dimen", d);
        d = mod (d, dimen) + 1;
        ## Build left sub tree
        if (mid > 1)
          left = r(1:mid-1);
          left_points = X(left,d);
          [val, left_idx] = sort (left_points);
          leftr = left(left_idx);
          ret.left = KDTreeSearcher.buildkdtree_recur (X, leftr, d, BS);
        endif
        ## Build right sub tree
        if (count > mid)
          right = r(mid+1:count);
          right_points = X(right,d);
          [val, right_idx] = sort (right_points);
          rightr = right(right_idx);
          ret.right = KDTreeSearcher.buildkdtree_recur (X, rightr, d, BS);
        endif
      endif
    endfunction

    function nn = findkdtree (tree, p, k, dist, distparam)
      X = tree.data;
      root = tree.root;
      nn = KDTreeSearcher.findkdtree_recur (X, root, p, [], k, ...
                                            dist, distparam);
    endfunction

    function nn = findkdtree_recur (X, node, p, nn, k, dist, distparam)
      point = node.point;
      d = node.dimen;
      if (X(point,d) > p(d))
        ## Search in left sub tree
        if (isfield (node, "left"))
          nn = KDTreeSearcher.findkdtree_recur (X, node.left, p, nn, k, ...
                                                 dist, distparam);
        endif
        ## Add current point if necessary
        farthest = KDTreeSearcher.kdtree_cand_farthest (X, p, nn, ...
                                                         dist, distparam);
        if (length (nn) < k || ...
            pdist2 (X(point,:), p, dist, distparam) <= ...
            pdist2 (X(farthest,:), p, dist, distparam))
          nn = KDTreeSearcher.kdtree_cand_insert (X, p, nn, k, point, dist, ...
                                                   distparam);
        endif
        ## Search in right sub tree if necessary
        farthest = KDTreeSearcher.kdtree_cand_farthest (X, p, nn, dist, ...
                                                         distparam);
        radius = pdist2 (X(farthest,:), p, dist, distparam);
        if (isfield (node, "right") && ...
            (length (nn) < k || p(d) + radius > X(point,d)))
          nn = KDTreeSearcher.findkdtree_recur (X, node.right, p, nn, k, ...
                                                 dist, distparam);
        endif
      else
        ## Search in right sub tree
        if (isfield (node, "right"))
          nn = KDTreeSearcher.findkdtree_recur (X, node.right, p, nn, k, ...
                                                 dist, distparam);
        endif
        ## Add current point if necessary
        farthest = KDTreeSearcher.kdtree_cand_farthest (X, p, nn, dist, ...
                                                         distparam);
        if (length (nn) < k || ...
            pdist2 (X(point,:), p, dist, distparam) <= ...
            pdist2 (X(farthest,:), p, dist, distparam))
          nn = KDTreeSearcher.kdtree_cand_insert (X, p, nn, k, point, dist, ...
                                                   distparam);
        endif
        ## Search in left sub tree if necessary
        farthest = KDTreeSearcher.kdtree_cand_farthest (X, p, nn, dist, ...
                                                         distparam);
        radius = pdist2 (X(farthest,:), p, dist, distparam);
        if (isfield (node, "left") && ...
            (length (nn) < k || p(d) - radius <= X(point,d)))
          nn = KDTreeSearcher.findkdtree_recur (X, node.left, p, nn, k, ...
                                                 dist, distparam);
        endif
      endif
    endfunction

    function farthest = kdtree_cand_farthest (X, p, cand, dist, distparam)
      if (isempty (cand))
        farthest = [];
      else
        D = pdist2 (X(cand,:), p, dist, distparam);
        [~, index] = max (D);
        farthest = cand(index);
      endif
    endfunction

    function inserted = kdtree_cand_insert (X, p, cand, k, point, dist, ...
                                             distparam)
      if (length (cand) < k)
        inserted = [cand; point];
      else
        farthest = KDTreeSearcher.kdtree_cand_farthest (X, p, cand, dist, ...
                                                         distparam);
        if (pdist2 (X(point,:), p, dist, distparam) < ...
            pdist2 (X(farthest,:), p, dist, distparam))
          cand(find (cand == farthest)) = point;
          inserted = cand;
        else
          inserted = cand;
        endif
      endif
    endfunction

  endmethods

endclassdef

## Demo Examples

%!demo
%! ## Create a KDTreeSearcher with Euclidean distance
%! X = [1, 2; 3, 4; 5, 6];
%! obj = KDTreeSearcher (X);
%! ## Find the nearest neighbor to [2, 3]
%! Y = [2, 3];
%! [idx, D] = knnsearch (obj, Y, 1);
%! disp ("Nearest neighbor index:"); disp (idx);
%! disp ("Distance:"); disp (D);
%! ## Find all points within radius 2
%! [idx, D] = rangesearch (obj, Y, 2);
%! disp ("Indices within radius:"); disp (idx);
%! disp ("Distances:"); disp (D);

%!demo
%! ## Create a KDTreeSearcher with Minkowski distance (P=3)
%! X = [0, 0; 1, 0; 2, 0];
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 3);
%! ## Find the nearest neighbor to [1, 0]
%! Y = [1, 0];
%! [idx, D] = knnsearch (obj, Y, 1);
%! disp ("Nearest neighbor index:"); disp (idx);
%! disp ("Distance:"); disp (D);

## Test Cases

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
%! assert (D{1}, sort (D_true(expected_idx)), 1e-10);

%!test
%! ## knnsearch with duplicates
%! X = [0, 0; 0, 0; 1, 0];
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! Y = [0, 0];
%! [idx, D] = knnsearch (obj, Y, 1, "IncludeTies", true);
%! assert (sort (idx{1}(:))', [1, 2]);
%! assert (D{1}', [0, 0], 1e-10);

%!test
%! ## rangesearch with 3D data
%! X = [0, 0, 0; 1, 0, 0; 0, 1, 0];
%! obj = KDTreeSearcher (X, "Distance", "cityblock");
%! Y = [0, 0, 0];
%! r = 1;
%! [idx, D] = rangesearch (obj, Y, r);
%! assert (sort (idx{1}(:))', [1, 2, 3]);
%! assert (D{1}', [0, 1, 1], 1e-10);

%!test
%! ## knnsearch with P = 2 (Euclidean equivalent)
%! X = [0, 0; 1, 1];
%! obj = KDTreeSearcher (X, "Distance", "minkowski", "P", 2);
%! Y = [0, 1];
%! [idx, D] = knnsearch (obj, Y, 1);
%! assert (idx, 2);
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
%! assert (D{1}, sort (D_true(expected_idx)), 1e-10);

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
%! assert (D{1}', [0, 0, 0], 1e-10);

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
%! assert (D{1}, D_true(expected_idx), 1e-10);

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
%! obj = KDTreeSearcher (ones(3,2)); obj.(1)
%!error<KDTreeSearcher.subsref: unrecognized property: 'invalid'.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.invalid

%!error<KDTreeSearcher.subsasgn: \(\) indexing not supported.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj(1) = 1
%!error<KDTreeSearcher.subsasgn: {} indexing not supported.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj{1} = 1
%!error<KDTreeSearcher.subsasgn: chained subscripts not allowed.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.X.Y = 1
%! obj = KDTreeSearcher (ones(3,2)); obj.(1) = 1
%!error<KDTreeSearcher.subsasgn: X is read-only and cannot be modified.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.X = 1
%!error<KDTreeSearcher.subsasgn: KDTree is read-only and cannot be modified.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.KDTree = 1
%!error<KDTreeSearcher.subsasgn: unsupported distance metric 'invalid'.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.Distance = "invalid"
%!error<KDTreeSearcher.subsasgn: Distance must be a string.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.Distance = 1
%!error<KDTreeSearcher.subsasgn: DistParameter must be a positive finite scalar for minkowski.> ...
%! obj = KDTreeSearcher (ones(3,2), "Distance", "minkowski"); obj.DistParameter = -1
%!error<KDTreeSearcher.subsasgn: DistParameter must be empty for this distance metric.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.DistParameter = 1
%!error<KDTreeSearcher.subsasgn: BucketSize must be a positive integer.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.BucketSize = 0
%!error<KDTreeSearcher.subsasgn: BucketSize must be a positive integer.> ...
%! obj = KDTreeSearcher(ones(3,2)); obj.BucketSize = -1
%!error<KDTreeSearcher.subsasgn: BucketSize must be a positive integer.> ...
%! obj = KDTreeSearcher(ones(3,2)); obj.BucketSize = 1.5
%!error<KDTreeSearcher.subsasgn: unrecognized property: 'invalid'.> ...
%! obj = KDTreeSearcher (ones(3,2)); obj.invalid = 1