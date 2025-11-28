## Copyright (C) 2021 Stefano Guidoni <ilguido@users.sf.net>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef CalinskiHarabaszEvaluation < ClusterCriterion
  ## -*- texinfo -*-
  ## @deftp {statistics} CalinskiHarabaszEvaluation
  ##
  ## Calinski-Harabasz clustering evaluation.
  ##
  ## A @code{CalinskiHarabaszEvaluation} object contains the results of
  ## evaluating clustering solutions using the Calinski-Harabasz criterion.
  ##
  ## The Calinski-Harabasz index (also known as the Variance Ratio Criterion) is
  ## determined by the ratio of the between-cluster sum of squares (SSB) to the
  ## within-cluster sum of squares (SSW).  A higher Calinski-Harabasz index
  ## value indicates a better clustering solution, implying that clusters are
  ## dense and well-separated.
  ##
  ## Create a @code{CalinskiHarabaszEvaluation} object by using the
  ## @code{evalclusters} function with the @qcode{"CalinskiHarabasz"} criterion.
  ##
  ## @seealso{evalclusters, ClusterCriterion, DaviesBouldinEvaluation,
  ## GapEvaluation, SilhouetteEvaluation}
  ## @end deftp

  properties (Access = protected)
    Centroids = {}; # a list of the centroids for every solution
  endproperties

  methods (Access = public)

    ## constructor
    function this = CalinskiHarabaszEvaluation (x, clust, KList)
      this@ClusterCriterion(x, clust, KList);

      this.CriterionName = "CalinskiHarabasz";
      this.evaluate(this.InspectedK); # evaluate the list of cluster numbers
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CalinskiHarabaszEvaluation} {@var{obj} =} addK (@var{obj}, @var{K})
    ##
    ## Add new cluster sizes for evaluation.
    ##
    ## @code{addK (@var{obj}, @var{K})} evaluates clustering solutions for the
    ## number of clusters specified in the vector @var{K} and adds them to the
    ## @code{CalinskiHarabaszEvaluation} object @var{obj}.
    ##
    ## @seealso{CalinskiHarabaszEvaluation, evalclusters}
    ## @end deftypefn
    function this = addK (this, K)
      addK@ClusterCriterion(this, K);

      ## if we have new data, we need a new evaluation
      if (this.OptimalK == 0)
        Centroids_tmp = {};
        pS = 0; # position shift of the elements of Centroids
        for iter = 1 : length (this.InspectedK)
          ## reorganize Centroids according to the new list of cluster numbers
          if (any (this.InspectedK(iter) == K))
            pS += 1;
          else
            Centroids_tmp{iter} = this.Centroids{iter - pS};
          endif
        endfor
        this.Centroids = Centroids_tmp;
        this.evaluate(K); # evaluate just the new cluster numbers
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CalinskiHarabaszEvaluation} {} plot (@var{obj})
    ## @deftypefnx {CalinskiHarabaszEvaluation} {@var{h} =} plot (@var{obj})
    ##
    ## Plot the clustering evaluation results.
    ##
    ## @code{plot (@var{obj})} plots the Calinski-Harabasz criterion values
    ## against the number of clusters.  The optimal number of clusters is marked
    ## with an asterisk.
    ##
    ## @code{@var{h} = plot (@var{obj})} additionally returns the handle to the
    ## plot axes.
    ##
    ## @seealso{CalinskiHarabaszEvaluation, evalclusters}
    ## @end deftypefn
    function h = plot (this)
      yLabel = sprintf ("%s value", this.CriterionName);
      h = gca ();
      hold on;
      plot (this.InspectedK, this.CriterionValues, "bo-");
      plot (this.OptimalK, this.CriterionValues(this.OptimalIndex), "b*");
      xlabel ("number of clusters");
      ylabel (yLabel);
      hold off;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CalinskiHarabaszEvaluation} {@var{obj} =} compact (@var{obj})
    ##
    ## Create a compact clustering evaluation object.
    ##
    ## This method is not yet implemented for
    ## @code{CalinskiHarabaszEvaluation} objects.
    ##
    ## @seealso{CalinskiHarabaszEvaluation, evalclusters}
    ## @end deftypefn
    function this = compact (this)
      warning (strcat ("CalinskiHarabaszEvaluation.compact:",...
                       " this method is not yet implemented."));
    endfunction

  endmethods

  methods (Access = protected)
    ## evaluate
    ## do the evaluation
    function this = evaluate (this, K)
      ## use complete observations only
      UsableX = this.X(find (this.Missing == false), :);
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
                error (strcat ("CalinskiHarabaszEvaluation: invalid return", ...
                               " value from custom clustering function."));
              endif
              this.ClusteringSolutions(:, iter) = ...
                this.ClusteringFunction(UsableX, this.InspectedK(iter));
            else
              switch (this.ClusteringFunction)
                case "kmeans"
                  [this.ClusteringSolutions(:, iter), this.Centroids{iter}] =...
                    kmeans (UsableX, this.InspectedK(iter),  ...
                    "Distance", "sqeuclidean", "EmptyAction", "singleton", ...
                    "Replicates", 5);

                case "linkage"
                  ## use clusterdata
                  this.ClusteringSolutions(:, iter) = clusterdata (UsableX, ...
                    "MaxClust", this.InspectedK(iter), ...
                    "Distance", "euclidean", "Linkage", "ward");
                  this.Centroids{iter} = this.computeCentroids (UsableX, iter);

                case "gmdistribution"
                  gmm = fitgmdist (UsableX, this.InspectedK(iter), ...
                        "SharedCov", true, "Replicates", 5);
                  this.ClusteringSolutions(:, iter) = cluster (gmm, UsableX);
                  this.Centroids{iter} = gmm.mu;

                otherwise
                  error (strcat ("CalinskiHarabaszEvaluation:" ...
                                 " unexpected error, report this bug."));
              endswitch
            endif
          endif
        endfor
      endif

      ## get the criterion values for every clustering solution
      for iter = 1 : length (this.InspectedK)
        ## do it only for the specified K values
        if (any (this.InspectedK(iter) == K))
          ## not defined for one cluster
          if (this.InspectedK(iter) == 1)
            this.CriterionValues(iter) = NaN;
            continue;
          endif

          ## Caliński-Harabasz index
          ## reference: calinhara function from the fpc package of R,
          ##            by Christian Hennig
          ##            https://CRAN.R-project.org/package=fpc
          W = zeros (columns (UsableX)); # between clusters covariance
          for i = 1 : this.InspectedK(iter)
            vIndicesI = find (this.ClusteringSolutions(:, iter) == i);
            ni = length (vIndicesI); # size of cluster i
            if (ni == 1)
              ## if the cluster has just one member the covariance is zero
              continue;
            endif
            ## weighted update of the covariance matrix
            W += cov (UsableX(vIndicesI, :)) * (ni - 1);
          endfor
          S = (this.NumObservations - 1) * cov (UsableX); # within clusters cov.
          B = S - W; # between clusters means

          ## tr(B) / tr(W) * (N-k) / (k-1)
          this.CriterionValues(iter) = (this.NumObservations - ...
            this.InspectedK(iter)) * trace (B) / ...
            ((this.InspectedK(iter) - 1) * trace (W));
        endif
      endfor

      [~, this.OptimalIndex] = max (this.CriterionValues);
      this.OptimalK = this.InspectedK(this.OptimalIndex(1));
      this.OptimalY = this.ClusteringSolutions(:, this.OptimalIndex(1));
    endfunction
  endmethods

  methods (Access = private)
    ## computeCentroids
    ## compute the centroids if they are not available by other means
    function C = computeCentroids (this, X, index)
      C = zeros (this.InspectedK(index), columns (X));

      for iter = 1 : this.InspectedK(index)
        vIndicesI = find (this.ClusteringSolutions(:, index) == iter);
        C(iter, :) = mean (X(vIndicesI, :));
      endfor
    endfunction
  endmethods
endclassdef

%!test
%! load fisheriris
%! eva = evalclusters (meas, "kmeans", "calinskiharabasz", "KList", [1:6]);
%! assert (class (eva), "CalinskiHarabaszEvaluation");
