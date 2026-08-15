## Copyright (C) 2008 Francesco Potortì <pot@gnu.org>
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
## @deftypefn  {statistics} {@var{y} =} linkage (@var{d})
## @deftypefnx {statistics} {@var{y} =} linkage (@var{d}, @var{method})
## @deftypefnx {statistics} {@var{y} =} linkage (@var{x})
## @deftypefnx {statistics} {@var{y} =} linkage (@var{x}, @var{method})
## @deftypefnx {statistics} {@var{y} =} linkage (@var{x}, @var{method}, @var{metric})
## @deftypefnx {statistics} {@var{y} =} linkage (@var{x}, @var{method}, @var{arglist})
##
## Produce a hierarchical clustering dendrogram.
##
## @var{d} is the dissimilarity matrix relative to n observations,
## formatted as a @math{(n-1)*n/2}x1 vector as produced by @code{pdist}.
## Alternatively, @var{x} contains data formatted for input to
## @code{pdist}, @var{metric} is a metric for @code{pdist} and
## @var{arglist} is a cell array containing arguments that are passed to
## @code{pdist}.
##
## @code{linkage} starts by putting each observation into a singleton
## cluster and numbering those from 1 to n.  Then it merges two
## clusters, chosen according to @var{method}, to create a new cluster
## numbered n+1, and so on until all observations are grouped into
## a single cluster numbered 2(n-1).  Row k of the
## (m-1)x3 output matrix relates to cluster n+k: the first
## two columns are the numbers of the two component clusters and column
## 3 contains their distance.
##
## When several pairs of clusters are equally close, which of them is merged
## first is not determined by the data, and the cluster numbers in the first two
## columns are therefore implementation-defined.  Only column 3, the sequence of
## merge distances, is reproducible across implementations, and even that is so
## only for the methods whose recomputation rule does not depend on the merge
## order (@qcode{"weighted"}, @qcode{"centroid"} and @qcode{"median"} do depend
## on it).  Code that must be portable should read column 3, or the cluster
## assignment obtained from @code{cluster}, rather than the raw numbering.
##
## @var{method} defines the way the distance between two clusters is
## computed and how they are recomputed when two clusters are merged:
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
## Unweighted pair group method with averaging (UPGMA).
## The mean distance between all pair of elements each belonging to one
## cluster.
##
## @item "weighted"
## Weighted pair group method with averaging (WPGMA).
## When two clusters A and B are joined together, the new distance to a
## cluster C is the mean between distances A-C and B-C.
##
## @item "centroid"
## Unweighted Pair-Group Method using Centroids (UPGMC).
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
## Ward's sum of squared deviations about the group mean (ESS).
## Also known as minimum variance or inner squared distance.
## Assumes Euclidean metric.  How much the moment of inertia of the
## merged cluster exceeds the sum of those of the individual clusters.
## @end table
##
## @strong{Reference}
## Ward, J. H. Hierarchical Grouping to Optimize an Objective Function
## J. Am. Statist. Assoc. 1963, 58, 236-244,
## @url{http://iv.slis.indiana.edu/sw/data/ward.pdf}.
##
## @seealso{pdist,squareform}
## @end deftypefn

function dgram = linkage (d, method = 'single', distarg, savememory)

  ## check the input
  if (nargin == 4) && (strcmpi (savememory, 'savememory'))
    warning ("Octave:linkage_savemem", ...
             'linkage: option ''savememory'' not implemented');
  elseif (nargin < 1) || (nargin > 3)
    print_usage ();
  endif

  if (isempty (d))
    error ("linkage: d cannot be empty");
  endif

  methods = struct ...
  ('name', { 'single'; 'complete'; 'average'; 'weighted';
            'centroid'; 'median'; 'ward' },
   'distfunc', {(@(x) min(x))                                     # single
                (@(x) max(x))                                     # complete
                (@(x,i,j,w) sum(diag(w([i,j]))*x)/sum(w([i,j])))  # average
                (@(x) mean(x))                                    # weighted
                (@massdist)                                       # centroid
                (@(x,i) massdist(x,i))                            # median
                (@inertialdist)                                   # ward
   });
  mask = strcmp (lower (method), {methods.name});
  if (! any (mask))
    error ("linkage: %s: unknown method", method);
  endif
  dist = {methods.distfunc}{mask};

  if (nargin >= 3 && ! isvector (d))
    if (ischar (distarg))
      d = pdist (d, distarg);
    elseif (iscell (distarg))
      d = pdist (d, distarg{:});
    else
      print_usage ();
    endif
  elseif (nargin < 3)
    if (! isvector (d))
      d = pdist (d);
    endif
  else
    print_usage ();
  endif

  d = squareform (d, 'tomatrix');       # dissimilarity NxN matrix
  n = rows (d);                         # the number of observations
  diagidx = sub2ind ([n,n], 1:n, 1:n);  # indices of diagonal elements
  d(diagidx) = Inf;                     # consider a cluster as far from itself
  ## For equal-distance nodes, the order in which clusters are
  ## merged is arbitrary.  Rotating the initial matrix produces an
  ## ordering similar to Matlab's.
  cname = n:-1:1;                       # cluster names in d
  d = rot90 (d, 2);                     # exchange low and high cluster numbers
  weight = ones (1, n);                 # cluster weights
  dgram = zeros (n-1, 3);               # clusters from n+1 to 2*n-1
  sz = n;                               # current matrix size (avoid size calls)
  mcase = find (mask);                  # pre-compute method case
  for cluster = n+1 : 2*n-1
    ## Find the two nearest clusters
    [~, midx] = min (d(:));
    ## Compute row/column indices directly (faster than ind2sub)
    c = ceil (midx / sz);
    r = midx - (c - 1) * sz;
    ## Here is the new cluster
    dgram(cluster-n, :) = [cname(r) cname(c) d(r, c)];
    ## Put it in place of the first one and remove the second
    cname(r) = cluster;
    cname(c) = [];
    ## Compute the new distances.
    ## (Octave-7+ needs switch stmt to avoid 'called with too many inputs' err.)
    d_rc = d([r, c], :);                # cache row slice
    switch mcase
      case {1, 2, 4}                    # 1 arg
        newd = dist(d_rc);
      case {3, 5, 7}                    # 4 args
        newd = dist(d_rc, r, c, weight);
      case 6                            # 2 args
        newd = dist(d_rc, r);
      otherwise
    endswitch
    newd(r) = Inf;                      # Take care of the diagonal element
    ## Put distances in place of the first ones, remove the second ones
    d(r,:) = newd;
    d(:,r) = newd';
    d(c,:) = [];
    d(:,c) = [];
    sz -= 1;                            # update cached size
    ## The new weight is the sum of the components' weights
    weight(r) += weight(c);
    weight(c) = [];
  endfor
  ## Sort the cluster numbers, as Matlab does
  dgram(:,1:2) = sort (dgram(:,1:2), 2);

  ## Check that distances are monotonically increasing
  if (any (diff (dgram(:,3)) < 0))
    warning ("Octave:clustering",
             "linkage: cluster distances do not monotonically increase\n\
        you should probably use a method different from \"%s\"", method);
  endif

endfunction


## Take two row vectors, which are the Euclidean distances of clusters I
## and J from the others.  Column I of second row contains the distance
## between clusters I and J.  The centre of gravity of the new cluster
## is on the segment joining the old ones. W are the weights of all
## clusters. Use the law of cosines to find the distances of the new
## cluster from all the others.
function y = massdist (x, i, j, w)
  x .^= 2;                              # Squared Euclidean distances
  if (nargin == 2)                      # Median distance
    qi = 0.5;                           # Equal weights ("weighted")
  else                                  # Centroid distance
    qi = 1 / (1 + w(j) / w(i));         # Proportional weights ("unweighted")
  endif
  y = sqrt (qi * x(1, :) + (1 - qi) * (x(2, :) - qi * x(2, i)));
endfunction


## Take two row vectors, which are the inertial distances of clusters I
## and J from the others.  Column I of second row contains the inertial
## distance between clusters I and J. The centre of gravity of the new
## cluster K is on the segment joining I and J.  W are the weights of
## all clusters.  Convert inertial to Euclidean distances, then use the
## law of cosines to find the Euclidean distances of K from all the
## other clusters, convert them back to inertial distances and return
## them.
function y = inertialdist (x, i, j, w)
  wi = w(i);                            # The cluster
  wj = w(j);                            # weights.
  s = [wi + w;                          # Sum of weights for
  wj + w];                              # all cluster pairs.
  p = [wi * w;                          # Product of weights for
  wj * w];                              # all cluster pairs.
  x = x.^2 .* s ./ p;                   # Convert inertial dist. to squared Eucl.
  sij = wi + wj;                        # Sum of weights of I and J
  qi = wi / sij;                        # Normalise the weight of I
  ## Squared Euclidean distances between all clusters and new cluster K
  x = qi * x(1, :) + (1 - qi) * (x(2, :) - qi * x(2, i));
  y = sqrt (x * sij .* w ./ (sij + w)); # convert Eucl. dist. to inertial
endfunction


%!shared x, t
%! x = reshape (mod (magic (6),5), [], 3);
%! t = 1e-6;

## Z(:,3) holds the merge heights and is the meaningful output; Z(:,1:2) are
## cluster labels whose values depend on how ties in the distance matrix are
## broken.  This fixture has only 20 distinct distances among its 66 pairs, so
## the labels are implementation-defined and are deliberately not asserted.
## Do not reintroduce cond (Z) here: it mixes the labels into the metric, is
## blind to the merge order (a row permutation cannot move a matrix's singular
## values) and ranges over 32.6 to 95.5 on this fixture under nothing but a
## different tie-break.  The heights below are exact against MATLAB R2024a.

%!test
%! Z = linkage (pdist (x));
%! assert_equal (Z(:,3), [1; 1; 1; 1.414214; 1.414214; 1.414214; ...
%!                        2.236068; 2.236068; 2.236068; 2.236068; 3], t);

%!test
%! Z = linkage (pdist (x), 'complete');
%! assert_equal (Z(:,3), [1; 1; 1.414214; 1.414214; 1.414214; 2.236068; ...
%!                        2.236068; 3.162278; 3.741657; 4.690416; 6], t);

%!test
%! Z = linkage (pdist (x), 'average');
%! assert_equal (Z(:,3), [1; 1; 1.207107; 1.414214; 1.414214; 1.962117; ...
%!                        2.236068; 2.948887; 3.081139; 3.515667; ...
%!                        4.177650], t);

## 'weighted' (McQuitty) feeds the Lance-Williams recursion with whichever pair
## merged first, so its later heights follow the tie-break.  Three pairs of this
## fixture lie at distance exactly 1 and observation 10 belongs to two of them,
## so at most two can merge as singletons and the choice is forced to be
## arbitrary: we take {9,10} then {1,7}, MATLAB takes {1,7} then {4,10}, and the
## heights part company from row 8.  Both replay correctly through the WPGMA
## recursion.  Only the first seven merges are tie-break invariant, so only they
## are asserted by value; a tie-free fixture would restore the rest.
%!test
%! Z = linkage (pdist (x), 'weighted');
%! assert_equal (Z(1:7,3), [1; 1; 1.207107; 1.414214; 1.414214; ...
%!                          2.030604; 2.236068], t);
%! assert_equal (all (diff (Z(:,3)) >= -eps), true);

%! lastwarn (); # Clear last warning before the test
%!warning <cluster distances> linkage (pdist (x), 'centroid');

## Regression values only -- 'centroid' and 'median' were not measured against
## MATLAB, so these pin our own behaviour rather than parity.  Both build a
## non-monotonic tree, which is what the warning above reports.
%!test
%! warning off Octave:clustering
%! Z = linkage (pdist (x), 'centroid');
%! assert_equal (Z(:,3), [1; 1; 1.118034; 1.414214; 1.224745; 1.885618; ...
%!                        2.236068; 2.708013; 2.980378; 3.041381; ...
%!                        3.529418], t);
%! warning on Octave:clustering

%!warning <cluster distances> linkage (pdist (x), 'median');

%!test
%! warning off Octave:clustering
%! Z = linkage (pdist (x), 'median');
%! assert_equal (Z(:,3), [1; 1; 1.118034; 1.414214; 1.224745; 1.952562; ...
%!                        2.236068; 2.452677; 3.041381; 3.057394; ...
%!                        3.163667], t);
%! warning on Octave:clustering

%!test
%! Z = linkage (pdist (x), 'ward');
%! assert_equal (Z(:,3), [1; 1; 1.290994; 1.414214; 1.414214; 2.236068; ...
%!                        2.309401; 3.511885; 4.690416; 5.228129; ...
%!                        7.713624], t);

## All four ward call forms must return the same tree, labels included.  MATLAB
## does not manage this -- its raw-data and distance-vector paths break the ties
## differently and return different trees for this fixture, at identical merge
## heights.  Ours agree exactly, so assert the whole of Z.
%!test
%! Z = linkage (pdist (x), 'ward');
%! assert_equal (linkage (x, 'ward', 'euclidean'), Z);
%! assert_equal (linkage (x, 'ward', {'euclidean'}), Z);
%! assert_equal (linkage (x, 'ward', {'minkowski', 2}), Z);

## Structural validity of Z, which no cond value could ever check.  The merge
## bookkeeping is shared by every method, so one method exercises it: each of
## the 2*n-2 labels is consumed exactly once and every row is sorted.
%!test
%! Z = linkage (pdist (x));
%! assert_equal (sort ([Z(:,1); Z(:,2)])', 1:22);
%! assert_equal (all (Z(:,1) < Z(:,2)), true);

## Additional tests for method/metric combinations
%!test
%! y = [1 2; 3 5; 4 6; 7 8; 9 11];
%! L = linkage (y, 'single', 'cityblock');
%! assert_equal (size (L), [4, 3]);
%! assert_equal (all (L(:,3) >= 0), true);  # distances non-negative

%!test
%! y = [1 2; 3 5; 4 6; 7 8; 9 11];
%! L = linkage (y, 'complete', 'cityblock');
%! assert_equal (size (L), [4, 3]);
%! assert_equal (all (diff (L(:,3)) >= -eps), true);  # monotonically increasing

%!test
%! y = [1 2; 3 5; 4 6; 7 8; 9 11];
%! L = linkage (y, 'average', 'chebychev');
%! assert_equal (size (L), [4, 3]);
%! assert_equal (all (L(:,3) >= 0), true);

%!test
%! y = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
%! L = linkage (y, 'weighted', {'minkowski', 3});
%! assert_equal (size (L), [3, 3]);
%! assert_equal (all (L(:,3) >= 0), true);

%!test
%! y = [1 0 1; 0 1 1; 1 1 0; 0 0 1];
%! L = linkage (y, 'single', 'cosine');
%! assert_equal (size (L), [3, 3]);
%! assert_equal (all (L(:,3) >= 0), true);

%!test
%! y = [1 2 3; 2 3 4; 5 6 7];
%! L = linkage (y, 'complete', 'correlation');
%! assert_equal (size (L), [2, 3]);
%! assert_equal (all (L(:,3) >= 0), true);

## Test with 2 observations (minimal case)
%!test
%! y = [1 2; 3 4];
%! L = linkage (y, 'single', 'euclidean');
%! assert_equal (size (L), [1, 3]);
%! assert_equal (L(1,1:2), [1, 2]);

## Test output structure: cluster indices are valid
%!test
%! y = rand (6, 3);
%! L = linkage (y, 'average', 'euclidean');
%! assert_equal (all (L(:,1) >= 1 & L(:,1) <= 11), true);  # valid cluster refs
%! assert_equal (all (L(:,2) >= 1 & L(:,2) <= 11), true);
%! assert_equal (all (L(:,1) < L(:,2)), true);  # sorted within rows
