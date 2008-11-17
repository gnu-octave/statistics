## Copyright (C) 2006, 2008  Bill Denney  <bill@denney.ws>
## Copyright (C) 2008  Francesco Potortì  <pot@gnu.org>
##
## This software is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3, or (at your option)
## any later version.
##
## This software is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this software; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{y} =} linkage (@var{x})
## @deftypefnx {Function File} {@var{y} =} linkage (@var{x}, @var{method})
##
## Return clusters generated from a distance vector created by the pdist
## function.
##
## @var{x} is the dissimilarity matrix relative to @var{n} observations,
## formatted as a @math{(n-1)*n/2}x1 vector as produced by @code{pdist}.
## @code{linkage} starts by putting each observation into a singleton
## cluster and numbering those from 1 to @var{n}.  Then it merges two
## clusters, chosen according to @var{method}, to create a new cluster
## numbered @var{n+1}, and so on until all observations are grouped into
## a single cluster numbered @var{2*n-1}.  Row @var{m} of the
## @math{m-1}x3 output matrix relates to cluster @math{n+m}: the first
## two columns are the numbers of the two component clusters and column
## 3 contains their distance.
##
## Methods define the way the distance between two clusters is computed:
##
## @table @samp
## @item "single" (default)
## Distance between two clusters is the minimum distance between two
## elements belonging each to one cluster.  Produces a cluster tree
## known as minimum spanning tree.
##
## @item "complete"
## Furthest distance between two elements belonging each to one cluster.
##
## @item "average"
## Unweighted pair group method with averaging (UPGMA)
## The mean distance between all pair of elements each belonging to one
## cluster.
##
## @item "weighted"
## Weighted pair group method with averaging (WPGMA)
##
## @item "centroid"
## Unweighted Pair-Group Method using Centroids (UPGMC)
## Assumes Euclidean metric.  The distance between cluster centroids,
## each centroid being the center of mass of a cluster.
##
## @item "median"
## Weighted pair-group method using centroids (WPGMC).
## Assumes Euclidean metric.  Distance between cluster centroids.  When
## two clusters are joined together, the new centroid is the midpoint
## between the joined centroids.
##
## @item "ward"
## Inner squared distance (minimum variance)
##
## @end table
##
## @seealso{cluster,pdist}
## @end deftypefn

## Author: Bill Denney <denney@...>

function y = linkage (x, method)

  ## check the input
  if (nargin < 1) || (nargin > 2)
    print_usage ();
  elseif (nargin < 2)
    method = "single";
  endif

  if (isempty (x))
    error ("linkage: x cannot be empty");
  elseif (~ isvector (x))
    error ("linkage: x must be a vector");
  endif

  ## Function distance must return a scalar from a vector and a row from
  ## a matrix.

  method = lower (method);
  switch (method)
    case "single"
      dist = @min;
    case "complete"
      dist = @max;
    case "median"
      dist = @mediandist;	# see below
    case { "weighted", "average", "centroid", "ward" }
      error ("linkage: %s is not yet implemented", method);
    otherwise
      error ("linkage: %s: unknown method", method);
  endswitch

  dissim = squareform (x, "tomatrix"); # dissimilarity matrix in square format
  n = rows (dissim);		       # the number of observations
  diagidx = sub2ind ([n,n], 1:n, 1:n); # indices of diagonal elements
  dissim(diagidx) = Inf;	# consider a cluster as far from itself
  cname = 1:n;			# cluster names in dissim
  y = zeros (n-1, 3);		# clusters from n+1 to 2*n-1
  for yidx = 1:n-1
    ## Find the two nearest clusters
    ## For equal-distance nodes, the order in which nodes are added to
    ## clusters is arbitrary.  The following code chooses the nodes so
    ## to mostly get the same ordering as in Matlab.  Note that the
    ## weighted methods can produce different clusterings depending on
    ## the order in which elements are added to clusters.
    [r c] = find (dissim == min (dissim(:)), 1, "last");
    ## Here is the new cluster
    y(yidx, :) = [cname(r) cname(c) dissim(r, c)];
    ## Put it in place of the first one and remove the second
    cname(r) = yidx + n;
    cname(c) = [];
    ## Same for the dissimilarities, take care of the diagonal element
    d = dist (dissim([r c], :));
    d(r) = Inf;
    dissim(r,:) = d;
    dissim(:,r) = d';
    dissim(c,:) = [];
    dissim(:,c) = [];
  endfor

  ## Check that distances are monotonically increasing
  if (any (diff (y(:,3)) < 0))
    warning ("clustering",
	     "linkage: cluster distances do not monotonically increase\n\
	you should maybe use a method different from \"%s\"", method);
  endif

endfunction

## Take two row vectors, which are the distances of clusters I and J
## from the others.  Column J of second row contains Inf, column J of
## first row contains distance between clusters I and J.  The centroid
## of the new cluster is midway between the old ones.  Use the law of
## cosines to find distances of the new ones from all the others.
function y = mediandist (x)
  interdist = x(1, x(2,:) == Inf); # distance between component clusters
  y = sqrt (sumsq (x) / 2 - interdist^2 / 4);
endfunction

%!shared x, y, t
%! x = [3 1.7; 1 1; 2 3; 2 2.5; 1.2 1; 1.1 1.5; 3 1];
%! y = reshape(mod(magic(6),5),[],3);
%! t = 1e-6;
%!assert (cond (linkage (pdist (x))),             56.981591, t);
%!assert (cond (linkage (pdist (y))),             27.771542, t);
%!assert (cond (linkage (pdist (x), "single")),   56.981591, t);
%!assert (cond (linkage (pdist (y), "single")),   27.771542, t);
%!assert (cond (linkage (pdist (x), "complete")), 26.681691, t);
%!assert (cond (linkage (pdist (y), "complete")), 21.726318, t);
%!assert (cond (linkage (pdist (x), "median")),   39.032841, t);
%!assert (cond (linkage (pdist (y), "median")),   26.159331, t);
