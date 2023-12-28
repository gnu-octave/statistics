## Copyright (C) 2021 Stefano Guidoni <ilguido@users.sf.net>
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

classdef GapEvaluation < ClusterCriterion
  ## -*- texinfo -*-
  ## @deftypefn {Function File} {@var{eva} =} evalclusters (@var{x}, @var{clust}, @qcode{gap})
  ## @deftypefnx {Function File} {@var{eva} =} evalclusters (@dots{}, @qcode{Name}, @qcode{Value})
  ##
  ## A gap object to evaluate clustering solutions.
  ##
  ## A @code{GapEvaluation} object is a @code{ClusterCriterion}
  ## object used to evaluate clustering solutions using the gap criterion,
  ## which is a mathematical formalization of the elbow method.
  ##
  ## List of public properties specific to @code{SilhouetteEvaluation}:
  ## @table @code
  ## @item @qcode{B}
  ## the number of reference datasets to generate.
  ##
  ## @item @qcode{Distance}
  ## a valid distance metric name, or a function handle as accepted by the
  ## @code{pdist} function.
  ##
  ## @item @qcode{ExpectedLogW}
  ## a vector of the expected values for the logarithm of the within clusters
  ## dispersion.
  ##
  ## @item @qcode{LogW}
  ## a vector of the values of the logarithm of the within clusters dispersion.
  ##
  ## @item @qcode{ReferenceDistribution}
  ## a valid name for the reference distribution, namely: @code{PCA} (default)
  ## or @code{uniform}.
  ##
  ## @item @qcode{SE}
  ## a vector of the standard error of the expected values for the logarithm
  ## of the within clusters dispersion.
  ##
  ## @item @qcode{SearchMethod}
  ## a valid name for the search method to use: @code{globalMaxSE} (default) or
  ## @code{firstMaxSE}.
  ##
  ## @item @qcode{StdLogW}
  ## a vector of the standard deviation of the expected values for the logarithm
  ## of the within clusters dispersion.
  ## @end table
  ##
  ## The best solution according to the gap criterion depends on the chosen
  ## search method.  When the search method is @code{globalMaxSE}, the chosen
  ## gap value is the smaller one which is inside a standard error from the
  ## max gap value; when the search method is @code{firstMaxSE}, the chosen
  ## gap value is the first one which is inside a standard error from the next
  ## gap value.
  ## @end deftypefn
  ##
  ## @seealso{ClusterCriterion, evalclusters}

  properties (GetAccess = public, SetAccess = private)
    B = 0; # number of reference datasets
    Distance = ""; # pdist parameter
    ReferenceDistribution = ""; # distribution to use as reference
    SearchMethod = ""; # the method do identify the optimal number of clusters
    ExpectedLogW = []; # expected value for the natural logarithm of W
    LogW = []; # natural logarithm of W
    SE = []; # standard error for the natural logarithm of W
    StdLogW = []; # standard deviation of the natural logarithm of W
  endproperties

  properties (Access = protected)
    DistanceVector = []; # vector of pdist distances
    mExpectedLogW = []; # the result of the Monte-Carlo simulations
  endproperties

  methods (Access = public)
    ## constructor
    function this = GapEvaluation (x, clust, KList, b = 100, ...
                    distanceMetric = "sqeuclidean", ...
                    referenceDistribution = "pca", searchMethod = "globalmaxse")
      this@ClusterCriterion(x, clust, KList);

      ## parsing the distance criterion
      if (ischar (distanceMetric))
        if (any (strcmpi (distanceMetric, {"sqeuclidean", "euclidean", ...
                 "cityblock", "cosine", "correlation", "hamming", "jaccard"})))
          this.Distance = lower (distanceMetric);

          ## kmeans can use only a subset
          if (strcmpi (clust, "kmeans") && any (strcmpi (this.Distance, ...
              {"euclidean", "jaccard"})))
            error (["GapEvaluation: invalid distance criterion '%s' "...
                    "for 'kmeans'"], distanceMetric);
          endif
        else
          error ("GapEvaluation: unknown distance criterion '%s'", ...
                 distanceMetric);
        endif
      elseif (isa (distanceMetric, "function_handle"))
        this.Distance = distanceMetric;

        ## kmeans cannot use a function handle
        if (strcmpi (clust, "kmeans"))
          error ("GapEvaluation: invalid distance criterion for 'kmeans'");
        endif
      elseif (isvector (distanceMetric) && isnumeric (distanceMetric))
        this.Distance = "";
        this.DistanceVector = distanceMetric; # the validity check is delegated

        ## kmeans cannot use a distance vector
        if (strcmpi (clust, "kmeans"))
          error (["GapEvaluation: invalid distance criterion for "...
                  "'kmeans'"]);
        endif
      else
        error ("GapEvaluation: invalid distance metric");
      endif

      ## B: number of Monte-Carlo iterations
      if (! isnumeric (b) || ! isscalar (b) || b != floor (b) || b < 1)
        error ("GapEvaluation: b must a be positive integer number");
      endif
      this.B = b;

      ## reference distribution
      if (! ischar (referenceDistribution) || ! any (strcmpi ...
          (referenceDistribution, {"pca", "uniform"})))
        error (["GapEvaluation: the reference distribution must be either" ...
                "'PCA' or 'uniform'"]);
      elseif (strcmpi (referenceDistribution, "pca"))
        warning (["GapEvaluation: 'PCA' distribution not implemented, " ...
                  "using 'uniform'"]);
      endif
      this.ReferenceDistribution = lower (referenceDistribution);

      if (! ischar (searchMethod) || ! any (strcmpi (searchMethod, ...
          {"globalmaxse", "firstmaxse"})))
        error (["evalclusters: the search method must be either" ...
                "'globalMaxSE' or 'firstMaxSE'"]);
      endif
      this.SearchMethod = lower (searchMethod);

      ## a matrix to store the results from the Monte-Carlo runs
      this.mExpectedLogW = zeros (this.B, length (this.InspectedK));

      this.CriterionName = "gap";
      this.evaluate(this.InspectedK); # evaluate the list of cluster numbers
    endfunction

    ## set functions

    ## addK
    ## add new cluster sizes to evaluate
    function this = addK (this, K)
      addK@ClusterCriterion(this, K);

      ## if we have new data, we need a new evaluation
      if (this.OptimalK == 0)
        mExpectedLogW_tmp = zeros (this.B, length (this.InspectedK));
        pS = 0; # position shift
        for iter = 1 : length (this.InspectedK)
          ## reorganize all the arrays according to the new list
          ## of cluster numbers
          if (any (this.InspectedK(iter) == K))
            pS += 1;
          else
            mExpectedLogW_tmp(:, iter) = this.mExpectedLogW(:, iter - pS);
          endif
        endfor
        this.mExpectedLogW = mExpectedLogW_tmp;

        this.evaluate(K); # evaluate just the new cluster numbers
      endif
    endfunction

    ## compact
    ## ...
    function this = compact (this)
      # FIXME: stub!
      warning ("GapEvaluation: compact is unavailable");
    endfunction

    ## plot
    ## plot the CriterionValues against InspectedK, show the standard deviation
    ## and return a handle to the plot
    function h = plot (this)
      yLabel = sprintf ("%s value", this.CriterionName);
      h = gca ();
      hold on;
      errorbar (this.InspectedK, this.CriterionValues, this.StdLogW);
      plot (this.InspectedK, this.CriterionValues, "bo");
      plot (this.OptimalK, this.CriterionValues(this.OptimalIndex), "b*");
      xlabel ("number of clusters");
      ylabel (yLabel);
      hold off;
    endfunction
  endmethods

  methods (Access = protected)
    ## evaluate
    ## do the evaluation
    function this = evaluate (this, K)
      ## Monte-Carlo runs
      for mcrun = 1 : (this.B + 1)
        ## use complete observations only
        UsableX = this.X(find (this.Missing == false), :);

        ## the last run use tha actual data,
        ## the others are Monte-Carlo runs with reconstructed data
        if (mcrun <= this.B)
          ## uniform distribution
          colMins = min (UsableX);
          colMaxs = max (UsableX);
          for col = 1 : columns (UsableX)
            UsableX(:, col) = colMins(col) + rand (this.NumObservations, 1) *...
                              (colMaxs(col) - colMins(col));
          endfor
        endif

        if (! isempty (this.ClusteringFunction))
          ## build the clusters
          for iter = 1 : length (this.InspectedK)
            ## do it only for the specified K values
            if (any (this.InspectedK(iter) == K))
              if (isa (this.ClusteringFunction, "function_handle"))
                ## custom function
                ClusteringSolution = ...
                  this.ClusteringFunction(UsableX, this.InspectedK(iter));
                if (ismatrix (ClusteringSolution) && ...
                    rows (ClusteringSolution) == this.NumObservations && ...
                    columns (ClusteringSolution) == this.P)
                  ## the custom function returned a matrix:
                  ## we take the index of the maximum value for every row
                  [~, this.ClusteringSolutions(:, iter)] = ...
                    max (ClusteringSolution, [], 2);
                elseif (iscolumn (ClusteringSolution) &&
                        length (ClusteringSolution) == this.NumObservations)
                  this.ClusteringSolutions(:, iter) = ClusteringSolution;
                elseif (isrow (ClusteringSolution) &&
                        length (ClusteringSolution) == this.NumObservations)
                  this.ClusteringSolutions(:, iter) = ClusteringSolution';
                else
                  error (["GapEvaluation: invalid return value from " ...
                          "custom clustering function"]);
                endif
                this.ClusteringSolutions(:, iter) = ...
                  this.ClusteringFunction(UsableX, this.InspectedK(iter));
              else
                switch (this.ClusteringFunction)
                  case "kmeans"
                    this.ClusteringSolutions(:, iter) = kmeans (UsableX, ...
                      this.InspectedK(iter), "Distance", this.Distance, ...
                      "EmptyAction", "singleton", "Replicates", 5);

                  case "linkage"
                    if (! isempty (this.Distance))
                      ## use clusterdata
                      Distance_tmp = this.Distance;
                      LinkageMethod = "average"; # for non euclidean methods
                      if (strcmpi (this.Distance, "sqeuclidean"))
                        ## pdist uses different names for its algorithms
                        Distance_tmp = "squaredeuclidean";
                        LinkageMethod = "ward";
                      elseif (strcmpi (this.Distance, "euclidean"))
                        LinkageMethod = "ward";
                      endif
                      this.ClusteringSolutions(:, iter) = clusterdata ...
                        (UsableX, "MaxClust", this.InspectedK(iter), ...
                        "Distance", Distance_tmp, "Linkage", LinkageMethod);
                    else
                      ## use linkage
                      Z = linkage (this.DistanceVector, "average");
                      this.ClusteringSolutions(:, iter) = ...
                           cluster (Z, "MaxClust", this.InspectedK(iter));
                    endif

                  case "gmdistribution"
                    gmm = fitgmdist (UsableX, this.InspectedK(iter), ...
                          "SharedCov", true, "Replicates", 5);
                    this.ClusteringSolutions(:, iter) = cluster (gmm, UsableX);

                  otherwise
                    ## this should not happen
                    error (["GapEvaluation: unexpected error, " ...
                           "report this bug"]);
                endswitch
              endif
            endif
          endfor
        endif

        ## get the gap values for every clustering
        distance_pdist = this.Distance;
        if (strcmpi (distance_pdist, "sqeuclidean"))
          distance_pdist = "squaredeuclidean";
        endif

        ## compute LogW
        for iter = 1 : length (this.InspectedK)
          ## do it only for the specified K values
          if (any (this.InspectedK(iter) == K))
            wk = 0;
            for r = 1 : this.InspectedK(iter)
              vIndicesR = find (this.ClusteringSolutions(:, iter) == r);
              nr = length (vIndicesR);
              Dr = pdist (UsableX(vIndicesR, :), distance_pdist);
              wk += sum (Dr) / (2 * nr);
            endfor
            if (mcrun <= this.B)
              this.mExpectedLogW(mcrun, iter) = log (wk);
            else
              this.LogW(iter) = log (wk);
            endif
          endif
        endfor
      endfor

      this.ExpectedLogW = mean (this.mExpectedLogW);
      this.SE = sqrt ((1 + 1 / this.B) * sumsq (this.mExpectedLogW - ...
                                               this.ExpectedLogW) / this.B);
      this.StdLogW = std (this.mExpectedLogW);
      this.CriterionValues = this.ExpectedLogW - this.LogW;

      this.OptimalIndex = this.gapSearch ();
      this.OptimalK = this.InspectedK(this.OptimalIndex(1));
      this.OptimalY = this.ClusteringSolutions(:, this.OptimalIndex(1));
    endfunction

    ## gapSearch
    ## find the best solution according to the gap method
    function ind = gapSearch (this)
      if (strcmpi (this.SearchMethod, "globalmaxse"))
        [gapmax, indgp] = max (this.CriterionValues);
        for iter = 1 : length (this.InspectedK)
          ind = iter;
          if (this.CriterionValues(iter) > (gapmax - this.SE(indgp)))
            return
          endif
        endfor
      elseif (strcmpi (this.SearchMethod, "firstmaxse"))
        for iter = 1 : (length (this.InspectedK) - 1)
          ind = iter;
          if (this.CriterionValues(iter) > (this.CriterionValues(iter + 1) - ...
                                            this.SE(iter + 1)))
            return
          endif
        endfor
      else
        ## this should not happen
        error (["GapEvaluation: unexpected error, please report this bug"]);
      endif
    endfunction
  endmethods
endclassdef
