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
## Produce a hierarchical clustering dendrogram from a distance vector
## created by the @code{pdist} function.
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
## When two clusters A and B are joined together, the new distance to a
## cluster C is the mean between distances A-C and B-C.
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
## NOT IMPLEMENTED
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

  methods = { "single", "complete", "average", "weighted", "centroid", "median" };
  method = lower (method);
  switch (method)

    case (methods)
      distfunction = ...
      {(@(x) min(x))			       # single
       (@(x) max(x))			       # complete
       (@(x,q) sum(dmult(q,x)) / sum(q))       # average
       (@(x) mean(x))			       # weighted
       (@(x,q) massdist(x,q))		       # centroid
       (@(x) massdist(x))};		       # median
      dist = distfunction {strcmp (method, methods)};
      dissim = squareform (x, "tomatrix"); # dissimilarity NxN matrix
      n = rows (dissim);		   # the number of observations
      diagidx = sub2ind ([n,n], 1:n, 1:n); # indices of diagonal elements
      dissim(diagidx) = Inf;	# consider a cluster as far from itself
      ## For equal-distance nodes, the order in which clusters are
      ## merged is arbitrary, but some methods can produce different
      ## clusterings depending on it.  Rotating the initial matrix
      ## produces an ordering more similar to Matlab's.
      dissim = rot90 (dissim, 2);
      cname = n:-1:1;		# cluster names in dissim
      weight = ones (1, n);	# cluster weights
      y = zeros (n-1, 3);	# clusters from n+1 to 2*n-1
      for yidx = 1:n-1
	## Find the two nearest clusters
	[m midx] = min (dissim(:));
	[r, c] = ind2sub (size (dissim), midx);
	## Here is the new cluster
	y(yidx, :) = [cname(r) cname(c) dissim(r, c)];
	## Put it in place of the first one and remove the second
	cname(r) = n + yidx;
	cname(c) = [];
	## Compute the new distances
	d = dist (dissim([r c], :), weight([r c]));
	d(r) = Inf;		# take care of the diagonal element
	## Put distances in place of the first ones, remove the second ones
	dissim(r,:) = d;
	dissim(:,r) = d';
	dissim(c,:) = [];
	dissim(:,c) = [];
	## The new weight is the sum of the components' weights
	weight(r) += weight(c);
	weight(c) = [];
      endfor
      ## Sort the cluster numbers, as Matlab does
      y(:,1:2) = sort (y(:,1:2), 2);

    case "ward"
      error ("linkage: %s is not yet implemented", method);

    otherwise
      error ("linkage: %s: unknown method", method);
  endswitch

  ## Check that distances are monotonically increasing
  if (any (diff (y(:,3)) < 0))
    warning ("clustering",
	     "linkage: cluster distances do not monotonically increase\n\
	you should probably use a method different from \"%s\"", method);
  endif

endfunction

## Take two row vectors, which are the Euclidean distances of clusters I
## and J from the others.  Column J of second row contains Inf, column J
## of first row contains the distance between clusters I and J.  The
## centroid of the new cluster is on the segment joining the old ones. W
## are the weights of clusters I and J.  Use the law of cosines to find
## the distances of the new cluster from all the others.
function y = massdist (x, w = [1 1])
  c = x(1, x(2,:) == Inf);	# distance between component clusters
  w /= sum (w);			# ratio of distance position
  q2 = w(2);
  y = sqrt (w(1)*x(1,:).^2 + q2*(x(2,:).^2 + (q2-1)*c^2));
endfunction


%!shared x, y, t
%! x = [3 1.7; 1 1; 2 3; 2 2.5; 1.2 1; 1.1 1.5; 3 1];
%! y = reshape(mod(magic(6),5),[],3);
%! t = 1e-6;
%! disp ("linkage: should emit 4 warnings\n\t about the \"centroid\" and \"median\" methods");
%!assert (cond (linkage (pdist (x))),             55.787151, t);
%!assert (cond (linkage (pdist (y))),             34.119045, t);
%!assert (cond (linkage (pdist (x), "complete")), 27.506710, t);
%!assert (cond (linkage (pdist (y), "complete")), 21.793345, t);
%!assert (cond (linkage (pdist (x), "average")),  35.766804, t);
%!assert (cond (linkage (pdist (y), "average")),  27.045012, t);
%!assert (cond (linkage (pdist (x), "weighted")), 36.257913, t);
%!assert (cond (linkage (pdist (y), "weighted")), 27.412889, t);
%!assert (cond (linkage (pdist (x), "centroid")), 39.104461, t);
%!assert (cond (linkage (pdist (y), "centroid")), 27.457477, t);
%!assert (cond (linkage (pdist (x), "median")),   39.671458, t);
%!assert (cond (linkage (pdist (y), "median")),   27.683325, t);
