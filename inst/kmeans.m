## Copyright (C) 2011 Soren Hauberg <soren@hauberg.org>
## Copyright (C) 2012 Daniel Ward <dwa012@gmail.com>
## Copyright (C) 2015-2016 Lachlan Andrew <lachlanbis@gmail.com>
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
## @deftypefn {} {[@var{idx}, @var{centers}, @var{sumd}, @var{dist}] =} kmeans (@var{data}, @var{k}, @var{param1}, @var{value1}, @dots{})
## Perform a @var{k}-means clustering of the @var{N}x@var{D} table @var{data}.
## If parameter @qcode{start} is specified, then @var{k} may be empty
## in which case @var{k} is set to the number of rows of @var{start}.
##
## The outputs are:
## @table @code
## @item @var{idx}
## An @var{N}x1 vector whose @var{i}th element is the class to which row @var{i}
## of @var{data} is assigned.
## 
## @item @var{centers}
## A @var{K}x@var{D} array whose @var{i}th row is the centroid of cluster
## @var{i}.
##
## @item @var{sumd}
## A @var{k}x1 vector whose @var{i}th entry is the sum of the distances
## from samples in cluster @var{i} to centroid @var{i}.
##
## @item @var{dist}
## An @var{N}x@var{k} matrix whose @var{i}@var{j}th element is
## the distance from sample @var{i} to centroid @var{j}.
## @end table 
## 
## The following parameters may be placed in any order.  Each parameter
## must be followed by its value.
## @table @code
## @item @var{Start}
## The initialization method for the centroids.
##  @table @code
##  @item @code{plus}
##        (Default) The k-means++ algorithm.
##  @item @code{sample}
#         A subset of @var{k} rows from @var{data},
##        sampled uniformly without replacement.
##  @item @code{cluster}
##        Perform a pilot clustering on 10% of the rows of @var{data}.
##  @item @code{uniform}
##        Each component of each centroid is drawn uniformly
##        from the interval between the maximum and minimum values of that
##        component within @var{data}.
##        This performs poorly and is implemented only for Matlab compatibility.
##  @item A
##        A @var{k}x@var{D}x@var{r} matrix, where @var{r} is the number of
##        replicates.
##  @end table
##
## @item @var{Replicates}
## An positive integer specifying the number of independent clusterings to
## perform.
## The output values are the values for the best clustering, i.e.,
## the one with the smallest value of @var{sumd}.
## If @var{Start} is numeric, then @var{Replicates} defaults to
# (and must equal) the size of the third dimension of @var{Start}.
## Otherwise it defaults to 1.
## 
## @item @var{MaxIter}
## The maximum number of iterations to perform for each replicate.
## If the maximum change of any centroid is less than 0.001, then
## the replicate terminates even if @var{MaxIter} iterations have no occurred.
## The default is 100.
##
## @item @var{Distance}
## The distance measure used for partitioning and calculating centroids.
##  @table @code
##  @item @qcode{sqeuclidean}
##  The squared Euclidean distance, i.e.,
##  the sum of the squares of the differences between corresponding components.
##  In this case, the centroid is the arithmetic mean of all samples in
##  its cluster.
##  This is the only distance for which this algorithm is truly "k-means".
##
##  @item @qcode{cityblock}
##  The sum metric, or L1 distance, i.e.,
##  the sum of the absolute differences between corresponding components.
##  In this case, the centroid is the median of all samples in its cluster.
##  This gives the k-medians algorithm.
##
##  @item @qcode{cosine}
##  (Documentation incomplete.)
##
## @item @qcode{correlation}
##  (Documentation incomplete.)
##
##  @item @qcode{hamming}
##  The number of components in which the sample and the centroid differ.
##  In this case, the centroid is the median of all samples in its cluster.
##  Unlike Matlab, Octave allows non-logical @var{data}.
## 
##  @end table
##
## @item @var{EmptyAction}
## What to do when a centroid is not the closest to any data sample.
##  @table @code
##  @item @qcode{error}
##        (Default) Throw an error.
##  @item @qcode{singleton}
##        Select the row of @var{data} that has the highest error and
##        use that as the new centroid.
##  @item @qcode{drop}
##        Remove the centroid, and continue computation with one fewer centroid.
##        The dimensions of the outputs @var{centroids} and @var{d}
##        are unchanged, with values for omitted centroids replaced by NA.
##        
##  @end table
## @end table
##
## Example:
##
##  [~,c] = kmeans (rand(10, 3), 2, "emptyaction", "singleton");
##
## @seealso{linkage}
## @end deftypefn

function [classes, centers, sumd, D] = kmeans (data, k, varargin)
  [reg, prop] = parseparams (varargin);

  ## defaults for options
  emptyaction = "error";
  start       = "plus";
  replicates  = 1;
  max_iter    = 100;
  distance    = "sqeuclidean";

  replicates_set_explicitly = false;

  ## Remove rows containing NaN / NA, but record which rows are used
  data_idx = ! any (isnan (data), 2);
  original_rows = rows (data);
  data = data(data_idx,:);

  #used for getting the number of samples
  n_rows = rows (data);

  #used for convergence of the centroids
  err = 1;

  ## Input checking, validate the matrix
  if (! isnumeric (data) || ! ismatrix (data) || ! isreal (data))
    error ("kmeans: first input argument must be a DxN real data matrix");
  elseif (! isnumeric (k))
    error ("kmeans: second argument must be numeric");
  endif

  ## Parse options
  while (length (prop) > 0)
    if (length (prop) < 2)
      error ("kmeans: Option '%s' has no argument", prop{1});
    endif
    switch (lower (prop{1}))
      case "emptyaction" emptyaction = prop{2};
      case "start"       start       = prop{2};
      case "maxiter"     max_iter    = prop{2};
      case "distance"    distance    = prop{2};

      case "replicates"  replicates  = prop{2};
                         replicates_set_explicitly = true;

      case {"display", "onlinephase", "options"}
        warning ("kmeans: Ignoring unimplemented option '%s'", prop{1});

      otherwise
        error ("kmeans: Unknown option %s", prop{1});
    endswitch
    prop = {prop{3:end}};
  endwhile

  ## Process options

  ## check for the 'emptyaction' property
  switch (emptyaction)
    case {"singleton", "error", "drop"}
      ;
    otherwise
      d = [", " disp(emptyaction)] (1:end-1);  # strip trailing \n
      if (length (d) > 20)
        d = "";
      endif
      error ("kmeans: unsupported empty cluster action parameter%s", d);
  endswitch

  ## check for the 'replicates' property
  if (! isnumeric (replicates) || ! isscalar (replicates)
     || ! isreal (replicates) || replicates < 1)
    d = [", " disp(replicates)] (1:end-1);     # strip trailing \n
    if (length (d) > 20)
      d = "";
    endif
    error ("kmeans: invalid number of replicates%s", d);
  endif

  ## check for the 'MaxIter' property
  if (! isnumeric (max_iter) || ! isscalar (max_iter)
     || ! isreal (max_iter) || max_iter < 1)
    d = [", " disp(max_iter)] (1:end-1);       # strip trailing \n
    if (length (d) > 20)
      d = "";
    endif
    error ("kmeans: invalid MaxIter%s", d);
  endif

  ## check for the 'start' property
  switch (lower (start))
    case {"sample", "plus", "cluster"}
      start = lower (start);
    case {"uniform"}
      start = "uniform";
      min_data = min (data);
      range = max (data) - min_data;
    otherwise
      if (! isnumeric (start))
        d = [", " disp(start)] (1:end-1);       # strip trailing \n
        if (length (d) > 20)
          d = "";
        endif
        error ("kmeans: invalid start parameter%s", d);
      endif
      if (isempty (k))
        k = rows (start);
      elseif (rows (start) != k)
        error ("kmeans: Number of initializers (%d) should match number of centroids (%d)", rows (start), k);
      endif
      if (replicates_set_explicitly)
        if (replicates != size (start, 3))
           error ("kmeans: The third dimension of the initializer (%d) should match the number of replicates (%d)", size (start, 3), replicates);
        endif
      else
        replicates = size (start, 3);
      endif
  endswitch

  ## check for the 'distance' property
  ## dist  returns the distance btwn each row of matrix x and a row vector c
  switch (lower (distance))
    case "sqeuclidean"
      dist = @(x, c) (sumsq (bsxfun (@minus, x, c), 2));
      centroid  = @(x) (mean (x,1));
    case "cityblock"
      dist = @(x, c) (sum (abs (bsxfun (@minus, x, c)), 2));
      centroid  = @(x) (median (x,1));
    case "cosine"
        ## Pre-normalize all data.
        ## (when Octave implements normr, will use  data = normr (data) )
      for i = 1:rows (data)
        data(i,:) = data(i,:) / sqrt (sumsq (data(i,:)));
      endfor
      dist = @(x, c) (1 - (x * c') ./ sqrt (sumsq (c)));
      centroid = @(x) (mean (x,1));   ## already normalized
    case "correlation"
        ## Pre-normalize all data.
      data = data - mean (data, 2);
        ## (when Octave implements normr, will use  data = normr (data) )
      for i = 1:rows (data)
        data(i,:) = data(i,:) / sqrt (sumsq (data(i,:)));
      endfor

      dist = @(x, c) (1 - (x * (c-mean (c))') ./ sqrt (sumsq (c-mean (c))));
      centroid  = @(x) (mean (x,1));   ## already normalized
    case "hamming"
      dist = @(x, c) (sum (bsxfun (@ne, x, c), 2));
      centroid  = @(x) (median (x,1));
    otherwise
      error ("kmeans: unsupported distance parameter %s", distance);
  endswitch

  ## Done processing options
  ########################################

  ## Now that  k  has been set (possibly by 'replicates' option), check/use it.

  if (! isscalar (k))
    error ("kmeans: second input argument must be a scalar");
  endif

  ## used to hold the distances from each sample to each class
  D = zeros (n_rows, k);

  best = Inf;
  best_centers = [];
  for rep = 1:replicates
    ## check for the 'start' property
    switch (lower (start))
      case "sample"
        idx = randperm (n_rows, k);
        centers = data(idx, :);
      case "plus"                  # k-means++, by Arthur and Vassilios(?)
        centers(1,:) = data(randi (n_rows),:);
	d = inf (n_rows, 1);       # Distance to nearest centroid so far
	for i = 2:k
	  d = min (d, dist (data, centers(i-1, :)));
	  centers(i,:) = data(find (cumsum (d) > rand * sum (d), 1), :);
	endfor
      case "cluster"
        idx = randperm (n_rows, max (k, ceil (n_rows/10)));
        [~, centers] = kmeans (data(idx,:), k, "start", "sample",
                               "distance", distance);
      case "uniform"
	# vectorised 'min_data + range .* rand'
        centers = bsxfun (@plus, min_data,
                          bsxfun (@times, range, rand (k, columns (data))));
      otherwise
        centers = start(:,:,rep);
    endswitch

    ## Run the algorithm
    iter = 1;

        ## Classify once before the loop; to set sumd, and  if  max_iter == 0
    ## Compute distances and classify
    for i = 1:k
      D (:, i) = dist (data, centers(i, :));
    endfor
    [~, classes] = min (D, [], 2);
    sumd = obj_cost (D, classes);

    while (err > 0.001 && iter++ <= max_iter)
      ## Calculate new centroids
      replaced_centroids = [];        ## Used by "emptyaction = singleton"
      for i = 1:k
        ## Get binary vector indicating membership in cluster i
        membership = (classes == i);

        ## Check for empty clusters
        if (! any (membership))
          switch emptyaction
            ## if 'singleton', then find the point that is the
            ## farthest from any centroid (and not replacing an empty cluster
            ## from earlier in this pass) and add it to the empty cluster
            case 'singleton'
             available = setdiff(1:n_rows, replaced_centroids);
             [~, idx] = max (min (D(available,:)'));
             idx = available(idx);
             replaced_centroids = [replaced_centroids, idx];

             classes(idx) = i;
             membership(idx)=1;

           ## if 'drop' then set C and D to NA
           case 'drop'
            centers(i,:) = NA;
            D(i,:) = NA;

           ## if 'error' then throw the error
            otherwise
              error ("kmeans: empty cluster created");
          endswitch
       endif ## end check for empty clusters

        ## update the centroids
        if (any (membership))      ## if we didn't "drop" the cluster
          centers(i, :) = centroid (data(membership, :));
        endif
      endfor

      ## Compute distances
      for i = 1:k
        D (:, i) = dist (data, centers(i, :));
      endfor

      ## Classify
      [~, classes] = min (D, [], 2);

      ## calculate the difference in the sum of distances
      new_sumd = obj_cost (D, classes);
      err  = sum (sumd - new_sumd);
      ## update the current sum of distances
      sumd = new_sumd;
    endwhile
    if (sum (sumd) < sum (best) || isinf (best))
      best = sumd;
      best_centers = centers;
    endif
  endfor
  centers = best_centers;
  sumd = best';

  final_classes = NA (original_rows,1);
  final_classes(data_idx) = classes;        ## other positions already NaN / NA
  classes = final_classes;
endfunction

## calculate the sum of within-class distances
function obj = obj_cost (D, classes)
  obj = zeros (1,columns (D));
  for i = 1:columns (D)
    idx = (classes == i);
    obj(i) = sum (D(idx,i));
  end
endfunction

## Test input parsing
%!error kmeans (rand (3,2), 4);

%!test
%! samples = 4;
%! dims = 3;
%! k = 2;
%! [cls, c, d, z] = kmeans (rand (samples,dims), k, "start", rand (k,dims, 5),
%!                          "emptyAction", "singleton");
%! assert (size (cls), [samples, 1]);
%! assert (size (c), [k, dims]);
%! assert (size (d), [k, 1]);
%! assert (size (z), [samples, k]);

%!test
%! samples = 4;
%! dims = 3;
%! k = 2;
%! [cls, c, d, z] = kmeans (rand (samples,dims), [], "start", rand (k,dims, 5),
%!                          "emptyAction", "singleton");
%! assert (size (cls), [samples, 1]);
%! assert (size (c), [k, dims]);
%! assert (size (d), [k, 1]);
%! assert (size (z), [samples, k]);

%!test
%! kmeans (rand (4,3), 2, "start", rand (2,3, 5), "replicates", 5,
%!         "emptyAction", "singleton");

%!error kmeans (rand (4,3), 2, "start", rand (2,3, 5), "replicates", 1);

%!error kmeans (rand (4,3), 2, "start", rand (2,2));

%!test
%! kmeans (rand (3,4), 2, "start", "sample", "emptyAction", "singleton");
%!test
%! kmeans (rand (3,4), 2, "start", "plus", "emptyAction", "singleton");
%!test
%! kmeans (rand (3,4), 2, "start", "cluster", "emptyAction", "singleton");
%!test
%! kmeans (rand (3,4), 2, "start", "uniform", "emptyAction", "singleton");

%!error kmeans (rand (3,4), 2, "start", "normal");

%!error kmeans (rand (4,3), 2, "replicates", i);
%!error kmeans (rand (4,3), 2, "replicates", -1);
%!error kmeans (rand (4,3), 2, "replicates", []);
%!error kmeans (rand (4,3), 2, "replicates", [1 2]);
%!error kmeans (rand (4,3), 2, "replicates", "one");

%!error kmeans (rand (4,3), 2, "MAXITER", i);
%!error kmeans (rand (4,3), 2, "MaxIter", -1);
%!error kmeans (rand (4,3), 2, "maxiter", []);
%!error kmeans (rand (4,3), 2, "maxiter", [1 2]);
%!error kmeans (rand (4,3), 2, "maxiter", "one");

%!test
%! kmeans (rand (4,3), 2, "distance", "sqeuclidean", "emptyAction", "singleton");

%!test
%! kmeans (rand (4,3), 2, "distance", "cityblock", "emptyAction", "singleton");

%!test
%! kmeans (rand (4,3), 2, "distance", "cosine", "emptyAction", "singleton");

%!test
%! kmeans (rand (4,3), 2, "distance", "correlation", "emptyAction", "singleton");

%!test
%! kmeans (rand (4,3), 2, "distance", "hamming", "emptyAction", "singleton");

%!error kmeans (rand (4,3), 2, "distance", "manhattan");

%!error <empty cluster created> kmeans ([1 0; 1.1 0], 2, "start", eye(2), "emptyaction", "error");

%!test
%! kmeans ([1 0; 1.1 0], 2, "start", eye(2), "emptyaction", "singleton");

%!test
%! [cls, c] = kmeans ([1 0; 2 0], 2, "start", [8,0;0,8], "emptyaction", "drop");
%! assert (cls, [1; 1]);
%! assert (c, [1.5, 0; NA, NA]);

%!error kmeans ([1 0; 1.1 0], 2, "start", eye(2), "emptyaction", "panic");

%!demo
%! ## Generate a two-cluster problem
%! C1 = randn (100, 2) + 1;
%! C2 = randn (100, 2) - 1;
%! data = [C1; C2];
%!
%! ## Perform clustering
%! [idx, centers] = kmeans (data, 2);
%!
%! ## Plot the result
%! figure;
%! plot (data (idx==1, 1), data (idx==1, 2), 'ro');
%! hold on;
%! plot (data (idx==2, 1), data (idx==2, 2), 'bs');
%! plot (centers (:, 1), centers (:, 2), 'kv', 'markersize', 10);
%! hold off;
