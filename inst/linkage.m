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
## Methods can be:
##
## @table @samp
## @item "single" (default)
## Shortest distance between two clusters (aka a minimum spanning tree)
##
## @item "complete"
## Furthest distance between two clusters
##
## @item "average"
## Unweighted pair group method with averaging (UPGMA)
##
## @item "weighted"
## Weighted pair group method with averaging (WPGMA)
## Same as Median, its formal definition does not require Euclidean
## metric. (Is this true?)
##
## @item "centroid"
## Centroid distance
##
## @item "median"
## Weighted pair-group method using centroids (WPGMC)
## To be used with Euclidean metric. (Is this true?)
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

  ## Function findfxn must return a scalar from a vector and a row from
  ## a matrix.

  method = lower (method);
  switch (method)
    case "single"
      ## this is just a minimal spanning tree
      findfxn = @min;
    case "complete"
      findfxn = @max;
    case { "median", "weighted", "average", "centroid", "ward" }
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
    [m midx] = min (dissim(:));	# the min distance and its first index
    [r, c] = ind2sub (size (dissim), midx); # row and col number
    ## Here is the new cluster
    y(yidx, :) = [cname(r) cname(c) dissim(r, c)];
    ## Add it as a new cluster index and remove the old ones
    cname(r) = yidx + n;
    cname(c) = [];
    ## Add the new dissimilarities and remove the old ones
    newdissim = findfxn (dissim([r c], :));
    newdissim(r) = Inf;
    dissim(r,:) = newdissim;
    dissim(:,r) = newdissim';
    dissim(c,:) = [];
    dissim(:,c) = [];
  endfor

endfunction

%!shared x, y, t
%! x = [3 1.7; 1 1; 2 3; 2 2.5; 1.2 1; 1.1 1.5; 3 1];
%! y = reshape(mod(magic(6),5),[],3);
%! t = 1e-6;
%!assert (cond (linkage (pdist (x))),             66.534612, t);
%!assert (cond (linkage (pdist (y))),             34.945071, t);
%!assert (cond (linkage (pdist (x), "single")),   66.534612, t);
%!assert (cond (linkage (pdist (y), "single")),   34.945071, t);
%!assert (cond (linkage (pdist (x), "complete")), 27.071750, t);
%!assert (cond (linkage (pdist (y), "complete")), 20.296516, t);
