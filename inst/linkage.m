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
## @var{x} is the dissimilarity matrix relative to @var{n} observations,
## formatted as a @math{(n-1)*n/2}x1 vector as produced by @code{pdist}.
## @code{linkage} starts by putting each observation into a singleton
## cluster and numbering those from 1 to @var{n}.  Then it merges two
## clusters to create a new cluster numbered @var{n+1}, and so on
## until all observations are grouped into a single cluster numbered
## @var{2*n-1}.  Row @var{m} of the @math{m-1}x3 output matrix relates
## to cluster @math{n+m}: the first two columns are the numbers of the
## two component clusters and column 3 contains their distance.
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
      error ("linkage: %s is not yet implemented", method);
      findfxn = @max;
    case { "median", "weighted" }
      findfxn = @mean;
    case { "average", "centroid", "ward" }
      error ("linkage: %s is not yet implemented", method);
    otherwise
      error ("linkage: %s: unknown method", method);
  endswitch

  dissim = squareform (x, "tomatrix");
  startsize = size (dissim, 1);
  y = zeros (startsize - 1, 3);
  cnameidx = 1:startsize;
  for yidx = 1:startsize-1
    ## Find the two nearest clusters.
    available = logical(tril (ones(size(dissim))) - eye(size(dissim)));
    [r, c] = find (min (dissim(available)) == dissim, 1);
    ## Here is the new cluster.
    y(yidx, :) = [cnameidx(r) cnameidx(c) dissim(r, c)];
    ## Add it as a new cluster index and remove the old ones.
    cnameidx(r) = yidx + startsize;
    cnameidx(c) = [];
    ## Update the dissimilarity matrix
    newdissim = findfxn (dissim([r c], :));
    dissim(r,:) = newdissim;
    dissim(:,r) = newdissim';
    dissim(r,r) = 0;
    dissim(c,:) = [];
    dissim(:,c) = [];
  endfor

endfunction

%!shared xy, t
%! xy = [3 1.7; 1 1; 2 3; 2 2.5; 1.2 1; 1.1 1.5; 3 1];
%! t = 1e-6;
%!assert (cond (linkage (pdist (xy))),             66.534612, t);
%!assert (cond (linkage (pdist (xy), "single")),   66.534612, t);
%!assert (cond (linkage (pdist (xy), "complete")), 27.071750, t);
