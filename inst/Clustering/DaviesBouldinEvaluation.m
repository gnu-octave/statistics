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

classdef DaviesBouldinEvaluation < ClusterCriterion
## -*- texinfo -*-
## @deftp {statistics} DaviesBouldinEvaluation
##
## Davies-Bouldin object to evaluate clustering solutions
##
## A @code{DaviesBouldinEvaluation} object is a @code{ClusterCriterion}
## object used to evaluate clustering solutions using the Davies-Bouldin
## criterion.
##
## The Davies-Bouldin criterion is based on the ratio between the distances
## between clusters and within clusters — distances between centroids and
## distances between each datapoint and its centroid.
##
## The best solution according to the Davies-Bouldin criterion is the one
## that produces the lowest Davies-Bouldin value.
##
## @seealso{evalclusters, ClusterCriterion, CalinskiHarabaszEvaluation,
## GapEvaluation, SilhouetteEvaluation}
## @end deftp

  properties (Access = private, Hidden)
    Centroids = {}; # a list of the centroids for every solution
  endproperties

  methods (Access = public)

    ## constructor
    function this = DaviesBouldinEvaluation (x, clust, KList)
      this@ClusterCriterion(x, clust, KList);

      this.CriterionName = "DaviesBouldin";
      this.evaluate(this.InspectedK); # evaluate the list of cluster numbers
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {DaviesBouldinEvaluation} {@var{obj} =} addK (@var{obj}, @var{K})
    ##
    ## Add new cluster numbers to inspect in the DaviesBouldinEvaluation object.
    ##
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
    ## @deftypefn  {DaviesBouldinEvaluation} {} plot (@var{obj})
    ## @deftypefnx {DaviesBouldinEvaluation} {@var{h} =} plot (@var{obj})
    ##
    ## Plot Davies-Bouldin evaluation results.
    ##
    ## Plot the CriterionValues against InspectedK from the
    ## DaviesBouldinEvaluation ClusterCriterion to the current plot.
    ## Returns an axes handle if requested.
    ##
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
    ## @deftypefn  {DaviesBouldinEvaluation} {@var{obj} =} compact (@var{obj})
    ##
    ## Return a compact DaviesBouldinEvaluation object (not implemented).
    ##
    ## @end deftypefn
    function this = compact (this)
      warning (["DaviesBouldinEvaluation.compact: this"...
                " method is not yet implemented."]);
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
                error (["DaviesBouldinEvaluation: invalid return value "...
                        "from custom clustering function"]);
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
                  error (["DaviesBouldinEvaluation: unexpected error, " ...
                          "report this bug"]);
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

          ## Davies-Bouldin value
          ## an evaluation of the ratio between within-cluster and
          ## between-cluster distances

          ## mean distances between cluster members and their centroid
          vD = zeros (this.InspectedK(iter), 1);
          for i = 1 : this.InspectedK(iter)
            vIndicesI = find (this.ClusteringSolutions(:, iter) == i);
            vD(i) = mean (vecnorm (UsableX(vIndicesI, :) - ...
                                   this.Centroids{iter}(i, :), 2, 2));
          endfor

          ## within-to-between cluster distance ratio
          Dij = zeros (this.InspectedK(iter));
          for i = 1 : (this.InspectedK(iter) - 1)
            for j = (i + 1) : this.InspectedK(iter)
              ## centroid to centroid distance
              dij = vecnorm (this.Centroids{iter}(i, :) - ...
                             this.Centroids{iter}(j, :));
              ## within-to-between cluster distance ratio for clusters i and j
              Dij(i, j) = (vD(i) + vD(j)) / dij;
            endfor
          endfor

          ## ( max_j D1j + max_j D2j + ... + max_j Dkj) / k
          this.CriterionValues(iter) = sum (max (Dij(i, :), [], 2)) / ...
                                            this.InspectedK(iter);
        endif
      endfor

      [~, this.OptimalIndex] = min (this.CriterionValues);
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
%! eva = evalclusters (meas, "kmeans", "DaviesBouldin", "KList", [1:6]);
%! assert (class (eva), "DaviesBouldinEvaluation");
