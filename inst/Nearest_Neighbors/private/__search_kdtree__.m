## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the Free
## Software Foundation; either version 3 of the License, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{indices}, @var{distances}] =} __search_kdtree__ (@var{node}, @var{query}, @var{k}, @var{X}, @var{dist}, @var{distparam}, @var{is_range})
## @deftypefnx {statistics} {[@var{indices}, @var{distances}] =} __search_kdtree__ (@dots{}, @var{r})
##
## Search a Kd-tree built by @code{__build_kdtree__} for one query point.
## Shared by @code{knnsearch}, @code{rangesearch} and @code{KDTreeSearcher}.
##
## With @var{is_range} false the @var{k} nearest neighbours are returned; with
## it true, every neighbour within the radius @var{r}.  @var{dist} names the
## metric and @var{distparam} is the exponent for @qcode{'minkowski'} and empty
## otherwise.
##
## The metric is resolved to a function of a leaf's rows once per call rather
## than dispatched per leaf, which is what makes the walk affordable: a leaf
## holds a few dozen points and a general distance call costs far more than the
## arithmetic it performs on them.
##
## @end deftypefn

function [indices, distances] = __search_kdtree__ (node, query, k, X, dist, ...
                                                   distparam, is_range, r)

  if (nargin < 8)
    r = Inf;
  endif

  if (strcmpi (dist, 'minkowski'))
    if (! (isscalar (distparam) && isnumeric (distparam)
                                && distparam > 0 && isfinite (distparam)))
      error (strcat ("__search_kdtree__: DISTPARAM must be a positive", ...
                     " finite scalar for minkowski."));
    endif
  elseif (! isempty (distparam))
    error (strcat ("__search_kdtree__: DISTPARAM must be empty for", ...
                   " non-minkowski metrics."));
  endif

  ## The metric, resolved once.  'manhattan' is the documented alias of
  ## 'cityblock' and the search accepts it wherever the callers let it through.
  switch (lower (dist))
    case 'euclidean'
      compute_dists = @(leaf_X) sqrt (sum ((leaf_X - query) .^ 2, 2));
    case {'cityblock', 'manhattan'}
      compute_dists = @(leaf_X) sum (abs (leaf_X - query), 2);
    case 'chebychev'
      compute_dists = @(leaf_X) max (abs (leaf_X - query), [], 2);
    case 'minkowski'
      p = distparam;
      compute_dists = @(leaf_X) sum (abs (leaf_X - query) .^ p, 2) .^ (1 / p);
    otherwise
      error ("__search_kdtree__: unsupported distance metric '%s'.", dist);
  endswitch

  indices = [];
  distances = [];
  search (node, 0);

  function search (node, depth)

    if (isempty (node))
      return;
    endif

    if (isfield (node, 'indices'))
      leaf_indices = node.indices;
      dists = compute_dists (X(leaf_indices,:));
      if (is_range)
        mask = dists <= r;
        indices = [indices; leaf_indices(mask)'];
        distances = [distances; dists(mask)];
      elseif (length (distances) >= k)
        ## The list is already full, so only a candidate beating its worst
        ## member is worth merging, and most leaves supply none.
        mask = dists < distances(end);
        if (any (mask))
          indices = [indices; leaf_indices(mask)'];
          distances = [distances; dists(mask)];
          [distances, sort_idx] = sort (distances);
          indices = indices(sort_idx);
          distances = distances(1:k);
          indices = indices(1:k);
        endif
      else
        indices = [indices; leaf_indices'];
        distances = [distances; dists];
        if (length (distances) >= k)
          [distances, sort_idx] = sort (distances);
          indices = indices(sort_idx);
          if (length (distances) > k)
            distances = distances(1:k);
            indices = indices(1:k);
          endif
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

      ## The far side can only hold something better if the splitting plane is
      ## itself closer than the worst neighbour held so far.
      plane_dist = abs (query(axis) - split_value);
      if (is_range)
        if (plane_dist <= r)
          search (further, depth + 1);
        endif
      elseif (length (distances) < k || plane_dist < distances(end))
        search (further, depth + 1);
      endif
    endif

  endfunction

endfunction
