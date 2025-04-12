## Copyright (C) 2011 Soren Hauberg <soren@hauberg.org>
## Copyright (C) 2012 Daniel Ward <dwa012@gmail.com>
## Copyright (C) 2015-2016 Lachlan Andrew <lachlanbis@gmail.com>
## Copyright (C) 2016 Michael Bentley <mikebentley15@gmail.com>
## Copyright (C) 2021 Stefano Guidoni <ilguido@users.sf.net>
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{idx} =} kmeans (@var{data}, @var{k})
## @deftypefnx {statistics} {[@var{idx}, @var{centers}] =} kmeans (@var{data}, @var{k})
## @deftypefnx {statistics} {[@var{idx}, @var{centers}, @var{sumd}] =} kmeans (@var{data}, @var{k})
## @deftypefnx {statistics} {[@var{idx}, @var{centers}, @var{sumd}, @var{dist}] =} kmeans (@var{data}, @var{k})
## @deftypefnx {statistics} {[@dots{}] =} kmeans (@var{data}, @var{k}, @var{param1}, @var{value1}, @dots{})
## @deftypefnx {statistics} {[@dots{}] =} kmeans (@var{data}, [], @qcode{"start"}, @var{start}, @dots{})
##
## Perform a @var{k}-means clustering of the @math{NxD} matrix @var{data}.
##
## If parameter @qcode{"start"} is specified, then @var{k} may be empty
## in which case @var{k} is set to the number of rows of @var{start}.
##
## The outputs are:
##
## @multitable @columnfractions 0.15 0.05 0.8
## @item @var{idx} @tab @tab An @math{Nx1} vector whose @math{i}-th element is
## the class to which row @math{i} of @var{data} is assigned.
##
## @item @var{centers} @tab @tab A @math{KxD} array whose @math{i}-th row is the
## centroid of cluster @math{i}.
##
## @item @var{sumd} @tab @tab A @math{kx1} vector whose @math{i}-th entry is the
## sum of the distances from samples in cluster @math{i} to centroid @math{i}.
##
## @item @var{dist} @tab @tab An @math{Nxk} matrix whose @math{i}@math{j}-th
## element is the distance from sample @math{i} to centroid @math{j}.
## @end multitable
##
## The following parameters may be placed in any order.  Each parameter
## must be followed by its value, as in Name-Value pairs.
##
## @multitable @columnfractions 0.15 0.02 0.83
## @headitem Name @tab @tab Description
## @item @qcode{"Start"} @tab @tab The initialization method for the centroids.
## @end multitable
##
## @multitable @columnfractions 0.04 0.19 0.02 0.75
## @headitem @tab Value @tab @tab Description
## @item @tab @qcode{"plus"} @tab @tab The k-means++ algorithm.  (Default)
## @item @tab @qcode{"sample"} @tab @tab A subset of @math{k} rows from
## @var{data}, sampled uniformly without replacement.
## @item @tab @qcode{"cluster"} @tab @tab Perform a pilot clustering on 10% of
## the rows of @var{data}.
## @item @tab @qcode{"uniform"} @tab @tab Each component of each centroid is
## drawn uniformly from the interval between the maximum and minimum values of
## that component within @var{data}.  This performs poorly and is implemented
## only for Matlab compatibility.
## @item @tab @var{numeric matrix} @tab @tab A @math{kxD} matrix of centroid
## starting locations.  The rows correspond to seeds.
## @item @tab @var{numeric array} @tab @tab A @math{kxDxr} array of centroid
## starting locations.  The third dimension invokes replication of the
## clustering routine.  Page @math{r} contains the set of seeds for replicate
## @math{r}.  @qcode{kmeans} infers the number of replicates (specified by the
## @qcode{"Replicates"} Name-Value pair argument) from the size of the third
## dimension.
## @end multitable
##
## @multitable @columnfractions 0.15 0.02 0.838
## @headitem Name @tab @tab Description
## @item @qcode{"Distance"} @tab @tab The distance measure used for partitioning
## and calculating centroids.
## @end multitable
##
## @multitable @columnfractions 0.04 0.19 0.02 0.75
## @headitem @tab Value @tab @tab Description
## @item @tab @qcode{"sqeuclidean"} @tab @tab The squared Euclidean distance.
## i.e. the sum of the squares of the differences between corresponding
## components.  In this case, the centroid is the arithmetic mean of all samples
## in its cluster.  This is the only distance for which this algorithm is truly
## "k-means".
## @item @tab @qcode{"cityblock"} @tab @tab The sum metric, or L1 distance,
## i.e. the sum of the absolute differences between corresponding components.
## In this case, the centroid is the median of all samples in its cluster.
## This gives the k-medians algorithm.
## @item @tab @qcode{"cosine"} @tab @tab One minus the cosine of the included
## angle between points (treated as vectors). Each centroid is the mean of the
## points in that cluster, after normalizing those points to unit Euclidean
## length.
## @item @tab @qcode{"correlation"} @tab @tab One minus the sample correlation
## between points (treated as sequences of values).  Each centroid is the
## component-wise mean of the points in that cluster, after centering and
## normalizing those points to zero mean and unit standard deviation.
## @item @tab @qcode{"hamming"} @tab @tab The number of components in which the
## sample and the centroid differ.  In this case, the centroid is the median of
## all samples in its cluster.  Unlike Matlab, Octave allows non-logical
## @var{data}.
## @end multitable
##
## @multitable @columnfractions 0.15 0.02 0.838
## @headitem Name @tab @tab Description
## @item @qcode{"EmptyAction"} @tab @tab What to do when a centroid is not the
## closest to any data sample.
## @end multitable
##
## @multitable @columnfractions 0.04 0.19 0.02 0.75
## @headitem @tab Value @tab @tab Description
## @item @tab @qcode{"error"} @tab @tab Throw an error.
## @item @tab @qcode{"singleton"} @tab @tab (Default) Select the row of
## @var{data} that has the highest error and use that as the new centroid.
## @item @tab @qcode{"drop"} @tab @tab Remove the centroid, and continue
## computation with one fewer centroid.  The dimensions of the outputs
## @var{centroids} and @var{d} are unchanged, with values for omitted centroids
## replaced by NaN.
## @end multitable
##
## @multitable @columnfractions 0.15 0.02 0.838
## @headitem Name @tab @tab Description
## @item @qcode{"Display"} @tab @tab Display a text summary.
## @end multitable
##
## @multitable @columnfractions 0.04 0.19 0.02 0.75
## @headitem @tab Value @tab @tab Description
## @item @tab @qcode{"off"} @tab @tab (Default) Display no summary.
## @item @tab @qcode{"final"} @tab @tab Display a summary for each clustering
## operation.
## @item @tab @qcode{"iter"} @tab @tab Display a summary for each iteration of a
## clustering operation.
## @end multitable
##
## @multitable @columnfractions 0.15 0.02 0.838
## @headitem Name @tab @tab Value
## @item @qcode{"Replicates"} @tab @tab A positive integer specifying the number
## of independent clusterings to perform.  The output values are the values for
## the best clustering, i.e., the one with the smallest value of @var{sumd}.
## If @var{Start} is numeric, then @var{Replicates} defaults to
## (and must equal) the size of the third dimension of @var{Start}.
## Otherwise it defaults to 1.
## @item @qcode{"MaxIter"} @tab @tab The maximum number of iterations to perform
## for each replicate.  If the maximum change of any centroid is less than
## 0.001, then the replicate terminates even if @var{MaxIter} iterations have no
## occurred.  The default is 100.
## @end multitable
##
## Example:
##
## [~,c] = kmeans (rand(10, 3), 2, "emptyaction", "singleton");
##
## @seealso{linkage}
## @end deftypefn

function [classes, centers, sumd, D] = kmeans (data, k, varargin)
  [reg, prop] = parseparams (varargin);

  ## defaults for options
  emptyaction = "singleton";
  start       = "plus";
  replicates  = 1;
  max_iter    = 100;
  distance    = "sqeuclidean";
  display     = "off";

  replicates_set_explicitly = false;

  ## Remove rows containing NaN / NA, but record which rows are used
  data_idx      = ! any (isnan (data), 2);
  original_rows = rows (data);
  data          = data(data_idx,:);

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
      case "emptyaction"
        emptyaction = prop{2};
      case "start"
        start = prop{2};
      case "maxiter"
        max_iter = prop{2};
      case "distance"
        distance = prop{2};
      case "replicates"
        replicates = prop{2};
        replicates_set_explicitly = true;
      case "display"
        display = prop{2};
      case {"onlinephase", "options"}
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
        error (["kmeans: Number of initializers (%d) " ...
                "should match number of centroids (%d)"], rows (start), k);
      endif
      if (replicates_set_explicitly)
        if (replicates != size (start, 3))
           error (["kmeans: The third dimension of the initializer (%d) " ...
                   "should match the number of replicates (%d)"], ...
                   size (start, 3), replicates);
        endif
      else
        replicates = size (start, 3);
      endif
  endswitch

  ## check for the 'distance' property
  ## dist  returns the distance btwn each row of matrix x and a row vector c
  switch (lower (distance))
    case "sqeuclidean"
      dist     = @(x, c) sumsq (bsxfun (@minus, x, c), 2);
      centroid = @(x) mean (x, 1);
    case "cityblock"
      dist     = @(x, c) sum (abs (bsxfun (@minus, x, c)), 2);
      centroid = @(x) median (x, 1);
    case "cosine"
        ## Pre-normalize all data.
        ## (when Octave implements normr, will use  data = normr (data) )
      for i = 1:rows (data)
        data(i,:) = data(i,:) / sqrt (sumsq (data(i,:)));
      endfor
      dist     = @(x, c) 1 - (x * c') ./ sqrt (sumsq (c));
      centroid = @(x) mean (x, 1);   ## already normalized
    case "correlation"
      ## Pre-normalize all data.
      data = data - mean (data, 2);
      ## (when Octave implements normr, will use  data = normr (data) )
      for i = 1:rows (data)
        data(i,:) = data(i,:) / sqrt (sumsq (data(i,:)));
      endfor
      dist     = @(x, c) 1 - (x * (c - mean (c))') ...
                          ./ sqrt (sumsq (c - mean (c)));
      centroid = @(x) mean (x, 1);   ## already normalized
    case "hamming"
      dist     = @(x, c) sum (bsxfun (@ne, x, c), 2);
      centroid = @(x) median (x, 1);
    otherwise
      error ("kmeans: unsupported distance parameter %s", distance);
  endswitch

  ## check for the 'display' property
  if (! strcmp (display, "off"))
    display = lower (display);
    switch (display)
      case {"off", "final"} ;
      case "iter"
        printf ("%6s\t%6s\t%8s\t%12s\n", "iter", "phase", "num", "sum");
      otherwise
        error ("kmeans: invalid display parameter %s", display);
    endswitch
  endif


  ## Done processing options
  ########################################

  ## Now that  k  has been set (possibly by 'replicates' option), check/use it.
  if (! isscalar (k))
    error ("kmeans: second input argument must be a scalar");
  endif

  ## used to hold the distances from each sample to each class
  D = zeros (n_rows, k);

  best         = Inf;
  best_centers = [];
  for rep = 1:replicates
    ## keep track of the number of data points that change class
    old_classes = zeros (rows (data), 1);
    n_changes = -1;

    ## check for the 'start' property
    switch (lower (start))
      case "sample"
        idx     = randperm (n_rows, k);
        centers = data(idx, :);
      case "plus"                  # k-means++, by Arthur and Vassilios(?)
        centers(1,:) = data(randi (n_rows),:);
        d            = inf (n_rows, 1);    # Distance to nearest centroid so far
        for i = 2:k
          d            = min (d, dist (data, centers(i - 1, :)));
          centers(i,:) = data(find (cumsum (d) > rand * sum (d), 1), :);
        endfor
      case "cluster"
        idx          = randperm (n_rows, max (k, ceil (n_rows / 10)));
        [~, centers] = kmeans (data(idx,:), k, "start", "sample", ...
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
    [D, classes, sumd] = update_dist (data, centers, D, k, dist);

    while (err > 0.001 && iter <= max_iter && n_changes != 0)
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
             available          = setdiff (1:n_rows, replaced_centroids);
             [~, idx]           = max (min (D(available,:)'));
             idx                = available(idx);
             replaced_centroids = [replaced_centroids, idx];

             classes(idx)    = i;
             membership(idx) = 1;

           ## if 'drop' then set C and D to NA
           case 'drop'
            centers(i,:) = NA;
            D(i,:)       = NA;

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

      ## Compute distances, classes and sums
      [D, classes, new_sumd] = update_dist (data, centers, D, k, dist);
      ## calculate the difference in the sum of distances
      err  = sum (sumd - new_sumd);
      ## update the current sum of distances
      sumd = new_sumd;
      ## compute the number of class changes
      n_changes = sum (old_classes != classes);
      old_classes = classes;

      ## display iteration status
      if (strcmp (display, "iter"))
        printf ("%6d\t%6d\t%8d\t%12.3f\n", (iter), 1, ...
          n_changes, sum (sumd));
      endif
      iter++;
    endwhile
    ## throw a warning if the algorithm did not converge
    if (iter > max_iter && err > 0.001 && n_changes != 0)
      warning ("kmeans: failed to converge in %d iterations", max_iter);
    endif

    if (sum (sumd) < sum (best) || isinf (best))
      best = sumd;
      best_centers = centers;
    endif

    ## display final results
    if (strcmp (display, "final"))
      printf ("Replicate %d, %d iterations, total sum of distances = %.3f.\n", ...
        rep, iter, sum (sumd));
    endif
  endfor
  centers = best_centers;
  ## Compute final distances, classes and sums
  [D, classes, sumd] = update_dist (data, centers, D, k, dist);

  ## display final results
  if (strcmp (display, "final") || strcmp (display, "iter"))
    printf ("Best total sum of distances = %.3f\n", sum (sumd));
  endif

  ## Return with equal size as inputs
  if (original_rows != rows (data))
    final           = NA (original_rows,1);
    final(data_idx) = classes;        ## other positions already NaN / NA
    classes         = final;
  endif

endfunction

## Update distances, classes and sums
function [D, classes, sumd] = update_dist (data, centers, D, k, dist)
    for i = 1:k
      D (:, i) = dist (data, centers(i, :));
    endfor
    [~, classes] = min (D, [], 2);
    ## calculate the sum of within-class distances
    sumd = zeros (k, 1);
    for i = 1:k
      sumd(i) = sum (D(classes == i,i));
    endfor
endfunction

%!demo
%! ## Generate a two-cluster problem
%! randn ("seed", 31)  # for reproducibility
%! C1 = randn (100, 2) + 1;
%! randn ("seed", 32)  # for reproducibility
%! C2 = randn (100, 2) - 1;
%! data = [C1; C2];
%!
%! ## Perform clustering
%! rand ("seed", 1)  # for reproducibility
%! [idx, centers] = kmeans (data, 2);
%!
%! ## Plot the result
%! figure;
%! plot (data (idx==1, 1), data (idx==1, 2), "ro");
%! hold on;
%! plot (data (idx==2, 1), data (idx==2, 2), "bs");
%! plot (centers (:, 1), centers (:, 2), "kv", "markersize", 10);
%! title ("A simple two-clusters example");
%! hold off;

%!demo
%! ## Cluster data using k-means clustering, then plot the cluster regions
%! ## Load Fisher's iris data set and use the petal lengths and widths as
%! ## predictors
%!
%! load fisheriris
%! X = meas(:,3:4);
%!
%! plot (X(:,1), X(:,2), "k*", "MarkerSize", 5);
%! title ("Fisher's Iris Data");
%! xlabel ("Petal Lengths (cm)");
%! ylabel ("Petal Widths (cm)");
%!
%! ## Cluster the data. Specify k = 3 clusters
%! rand ("seed", 1)  # for reproducibility
%! [idx, C] = kmeans (X, 3);
%! x1 = min (X(:,1)):0.01:max (X(:,1));
%! x2 = min (X(:,2)):0.01:max (X(:,2));
%! [x1G, x2G] = meshgrid (x1, x2);
%! XGrid = [x1G(:), x2G(:)];
%!
%! idx2Region = kmeans (XGrid, 3, "MaxIter", 10, "Start", C);
%! figure;
%! gscatter (XGrid(:,1), XGrid(:,2), idx2Region, ...
%!           [0, 0.75, 0.75; 0.75, 0, 0.75; 0.75, 0.75, 0], "..");
%! hold on;
%! plot (X(:,1), X(:,2), "k*", "MarkerSize", 5);
%! title ("Fisher's Iris Data");
%! xlabel ("Petal Lengths (cm)");
%! ylabel ("Petal Widths (cm)");
%! legend ("Region 1", "Region 2", "Region 3", "Data", "Location", "SouthEast");
%! hold off

%!demo
%! ## Partition Data into Two Clusters
%!
%! randn ("seed", 1)  # for reproducibility
%! r1 = randn (100, 2) * 0.75 + ones (100, 2);
%! randn ("seed", 2)  # for reproducibility
%! r2 = randn (100, 2) * 0.5 - ones (100, 2);
%! X = [r1; r2];
%!
%! plot (X(:,1), X(:,2), ".");
%! title ("Randomly Generated Data");
%! rand ("seed", 1)  # for reproducibility
%! [idx, C] = kmeans (X, 2, "Distance", "cityblock", ...
%!                          "Replicates", 5, "Display", "final");
%! figure;
%! plot (X(idx==1,1), X(idx==1,2), "r.", "MarkerSize", 12);
%! hold on
%! plot(X(idx==2,1), X(idx==2,2), "b.", "MarkerSize", 12);
%! plot (C(:,1), C(:,2), "kx", "MarkerSize", 15, "LineWidth", 3);
%! legend ("Cluster 1", "Cluster 2", "Centroids", "Location", "NorthWest");
%! title ("Cluster Assignments and Centroids");
%! hold off

%!demo
%! ## Assign New Data to Existing Clusters
%!
%! ## Generate a training data set using three distributions
%! randn ("seed", 5)  # for reproducibility
%! r1 = randn (100, 2) * 0.75 + ones (100, 2);
%! randn ("seed", 7)  # for reproducibility
%! r2 = randn (100, 2) * 0.5 - ones (100, 2);
%! randn ("seed", 9)  # for reproducibility
%! r3 = randn (100, 2) * 0.75;
%! X = [r1; r2; r3];
%!
%! ## Partition the training data into three clusters by using kmeans
%!
%! rand ("seed", 1)  # for reproducibility
%! [idx, C] = kmeans (X, 3);
%!
%! ## Plot the clusters and the cluster centroids
%!
%! gscatter (X(:,1), X(:,2), idx, "bgm", "***");
%! hold on
%! plot (C(:,1), C(:,2), "kx");
%! legend ("Cluster 1", "Cluster 2", "Cluster 3", "Cluster Centroid")
%!
%! ## Generate a test data set
%! randn ("seed", 25)  # for reproducibility
%! r1 = randn (100, 2) * 0.75 + ones (100, 2);
%! randn ("seed", 27)  # for reproducibility
%! r2 = randn (100, 2) * 0.5 - ones (100, 2);
%! randn ("seed", 29)  # for reproducibility
%! r3 = randn (100, 2) * 0.75;
%! Xtest = [r1; r2; r3];
%!
%! ## Classify the test data set using the existing clusters
%! ## Find the nearest centroid from each test data point by using pdist2
%!
%! D = pdist2 (C, Xtest, "euclidean");
%! [group, ~] = find (D == min (D));
%!
%! ## Plot the test data and label the test data using idx_test with gscatter
%!
%! gscatter (Xtest(:,1), Xtest(:,2), group, "bgm", "ooo");
%! box on;
%! legend ("Cluster 1", "Cluster 2", "Cluster 3", "Cluster Centroid", ...
%!         "Data classified to Cluster 1", "Data classified to Cluster 2", ...
%!         "Data classified to Cluster 3", "Location", "NorthWest");
%! title ("Assign New Data to Existing Clusters");

## Test output
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
%! [cls, c] = kmeans ([1 0; 2 0], 2, "start", [8,0;0,8], "emptyaction", "drop");
%! assert (cls, [1; 1]);
%! assert (c, [1.5, 0; NA, NA]);

%!test
%! kmeans (rand (4,3), 2, "start", rand (2,3, 5), "replicates", 5,
%!         "emptyAction", "singleton");
%!test
%! kmeans (rand (3,4), 2, "start", "sample", "emptyAction", "singleton");
%!test
%! kmeans (rand (3,4), 2, "start", "plus", "emptyAction", "singleton");
%!test
%! kmeans (rand (3,4), 2, "start", "cluster", "emptyAction", "singleton");
%!test
%! kmeans (rand (3,4), 2, "start", "uniform", "emptyAction", "singleton");
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
%!test
%! kmeans ([1 0; 1.1 0], 2, "start", eye(2), "emptyaction", "singleton");

## Test input validation
%!error kmeans (rand (3,2), 4);
%!error kmeans ([1 0; 1.1 0], 2, "start", eye(2), "emptyaction", "panic");
%!error kmeans (rand (4,3), 2, "start", rand (2,3, 5), "replicates", 1);
%!error kmeans (rand (4,3), 2, "start", rand (2,2));
%!error kmeans (rand (4,3), 2, "distance", "manhattan");
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
%!error <empty cluster created> kmeans ([1 0; 1.1 0], 2, "start", eye(2), "emptyaction", "error");
