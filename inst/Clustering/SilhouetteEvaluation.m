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

classdef SilhouetteEvaluation < ClusterCriterion
  ## -*- texinfo -*-
  ## @deftp {statistics} SilhouetteEvaluation
  ##
  ## Silhouette evaluation for clustering
  ##
  ## The @code{SilhouetteEvaluation} class implements an object to evaluate
  ## clustering solutions using the silhouette criterion.  A
  ## @code{SilhouetteEvaluation} object is a @code{ClusterCriterion} object
  ## that computes silhouette values for clustering solutions and selects the
  ## best number of clusters as the one with the highest average silhouette
  ## value.
  ##
  ## Create a @code{SilhouetteEvaluation} object by using the
  ## @code{evalclusters} function or the class constructor.
  ##
  ## List of public properties specific to @code{SilhouetteEvaluation}:
  ## @table @code
  ## @item @qcode{Distance}
  ## A valid distance metric name (string), a function handle, or a numeric
  ## vector as returned by @code{pdist}.  This specifies how pairwise
  ## distances are computed.
  ##
  ## @item @qcode{ClusterPriors}
  ## A character vector specifying how to evaluate silhouette values across
  ## clusters: @qcode{"empirical"} (default) uses empirical cluster priors,
  ## or @qcode{"equal"} treats clusters equally.
  ##
  ## @item @qcode{ClusterSilhouettes}
  ## A cell array containing silhouette values for each observation for each
  ## inspected cluster number.
  ## @end table
  ##
  ## The best clustering solution according to the silhouette criterion is the
  ## one that yields the highest average silhouette value.
  ##
  ## @seealso{evalclusters, ClusterCriterion, CalinskiHarabaszEvaluation,
  ## DaviesBouldinEvaluation, GapEvaluation}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {SilhouetteEvaluation} {property:} Distance
    ##
    ## Distance measure
    ##
    ## A string naming a distance metric, a function handle that computes
    ## distances, or a numeric vector as produced by @code{pdist}.  This
    ## property is read-only.
    ##
    ## @end deftp
    Distance = "";

    ## -*- texinfo -*-
    ## @deftp {SilhouetteEvaluation} {property:} ClusterPriors
    ##
    ## Cluster prior handling
    ##
    ## Specifies how cluster-level silhouette aggregation is computed.  Valid
    ## values are @qcode{"empirical"} (default) and @qcode{"equal"}.  This
    ## property is read-only.
    ##
    ## @end deftp
    ClusterPriors = "";

    ## -*- texinfo -*-
    ## @deftp {SilhouetteEvaluation} {property:} ClusterSilhouettes
    ##
    ## Silhouette values
    ##
    ## A cell array where each element contains the silhouette values for the
    ## observations of a given clustering (corresponding to an inspected K).
    ## This property is read-only.
    ##
    ## @end deftp
    ClusterSilhouettes = {};
  endproperties

  properties (Access = protected)
    ## -*- texinfo -*-
    ## @deftp {SilhouetteEvaluation} {property:} DistanceVector
    ##
    ## Precomputed distance vector
    ##
    ## If a numeric vector is supplied as the distance metric it is stored
    ## here and used instead of computing distances via @code{pdist}.  This
    ## property is read-only.
    ##
    ## @end deftp
    DistanceVector = [];
  endproperties

  methods (Access = public)

    ## constructor
    ## -*- texinfo -*-
    ## @deftypefn  {statistics} {@var{obj} =} SilhouetteEvaluation (@var{x}, @var{clust}, @var{KList})
    ## @deftypefnx {statistics} {@var{obj} =} SilhouetteEvaluation (@dots{}, @var{Name}, @var{Value})
    ##
    ## Create a @code{SilhouetteEvaluation} object to evaluate clustering
    ## solutions for data @var{x} using clustering method @var{clust} over the
    ## list of cluster numbers @var{KList}.
    ##
    ## @itemize
    ## @item
    ## @var{x} is an @math{NxP} numeric matrix of observations (rows) and
    ## predictors (columns).
    ## @item
    ## @var{clust} is a string naming the clustering method (for example
    ## @qcode{"kmeans"}, @qcode{"linkage"}, or a custom function handle).
    ## @item
    ## @var{KList} is a vector of positive integers specifying the cluster
    ## numbers to inspect.
    ## @end itemize
    ##
    ## Optional name-value pairs:
    ##
    ## @multitable @columnfractions 0.20 0.02 0.78
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"Distance"} @tab @tab Distance metric name, function handle,
    ## or numeric pdist vector.  Default: @qcode{"sqeuclidean"}.
    ##
    ## @item @qcode{"ClusterPriors"} @tab @tab Either @qcode{"empirical"}
    ## (default) or @qcode{"equal"}.
    ## @end multitable
    ##
    ## @seealso{silhouette, evalclusters, ClusterCriterion}
    ## @end deftypefn
    function this = SilhouetteEvaluation (x, clust, KList, ...
                    distanceMetric = "sqeuclidean", clusterPriors = "empirical")
      this@ClusterCriterion(x, clust, KList);

      ## parsing the distance criterion
      if (ischar (distanceMetric))
        if (any (strcmpi (distanceMetric, {"sqeuclidean", ...
                  "euclidean", "cityblock", "cosine", "correlation", ...
                  "hamming", "jaccard"})))
          this.Distance = lower (distanceMetric);

          ## kmeans can use only a subset
          if (strcmpi (clust, "kmeans") && any (strcmpi (this.Distance, ...
              {"euclidean", "jaccard"})))
            error (["SilhouetteEvaluation: invalid distance criterion '%s' "...
                    "for 'kmeans'"], distanceMetric);
          endif
        else
          error ("SilhouetteEvaluation: unknown distance criterion '%s'", ...
                 distanceMetric);
        endif
      elseif (isa (distanceMetric, "function_handle"))
        this.Distance = distanceMetric;

        ## kmeans cannot use a function handle
        if (strcmpi (clust, "kmeans"))
          error (["SilhouetteEvaluation: invalid distance criterion for "...
                  "'kmeans'"]);
        endif
      elseif (isvector (distanceMetric) && isnumeric (distanceMetric))
        this.Distance = "";
        this.DistanceVector = distanceMetric; # the validity check is delegated

        ## kmeans cannot use a distance vector
        if (strcmpi (clust, "kmeans"))
          error (["SilhouetteEvaluation: invalid distance criterion for "...
                  "'kmeans'"]);
        endif
      else
        error ("SilhouetteEvaluation: invalid distance metric");
      endif

      ## parsing the prior probabilities of each cluster
      if (ischar (distanceMetric))
        if (any (strcmpi (clusterPriors, {"empirical", "equal"})))
          this.ClusterPriors = lower (clusterPriors);
        else
          error (["SilhouetteEvaluation: unknown prior probability criterion"...
                 " '%s'"], clusterPriors);
        endif
      else
        error ("SilhouetteEvaluation: invalid prior probabilities");
      endif

      this.CriterionName = "silhouette";
      this.evaluate(this.InspectedK); # evaluate the list of cluster numbers
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {SilhouetteEvaluation} {@var{obj} =} addK (@var{obj}, @var{K})
    ##
    ## Add a new cluster number to inspect and re-evaluate silhouette values
    ## for newly added cluster numbers.
    ##
    ## @end deftypefn
    function this = addK (this, K)
      addK@ClusterCriterion(this, K);

      ## if we have new data, we need a new evaluation
      if (this.OptimalK == 0)
        ClusterSilhouettes_tmp = {};
        pS = 0; # position shift of the elements of ClusterSilhouettes
        for iter = 1 : length (this.InspectedK)
          ## reorganize ClusterSilhouettes according to the new list
          ## of cluster numbers
          if (any (this.InspectedK(iter) == K))
            pS += 1;
          else
            ClusterSilhouettes_tmp{iter} = this.ClusterSilhouettes{iter - pS};
          endif
        endfor
        this.ClusterSilhouettes = ClusterSilhouettes_tmp;
        this.evaluate(K); # evaluate just the new cluster numbers
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {SilhouetteEvaluation} {} plot (@var{obj})
    ## @deftypefnx {SilhouetteEvaluation} {@var{h} =} plot (@var{obj})
    ##
    ## Plot the silhouette evaluation results.
    ##
    ## Plot the criterion values (average silhouette) against inspected cluster
    ## numbers (@code{InspectedK}) for the given @var{obj}.  Optionally returns
    ## the axis handle for the plot.
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
    ## @deftypefn {SilhouetteEvaluation} {@var{obj} =} compact (@var{obj})
    ##
    ## Return a compact SilhouetteEvaluation object (not implemented yet).
    ##
    ## @end deftypefn
    function this = compact (this)
      warning (["SilhouetteEvaluation.compact: this"...
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
                error (["SilhouetteEvaluation: invalid return value from " ...
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
                    this.ClusteringSolutions(:, iter) = clusterdata (UsableX,...
                      "MaxClust", this.InspectedK(iter), ...
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
                  error (["SilhouetteEvaluation: unexpected error, " ...
                         "report this bug"]);
              endswitch
            endif
          endif
        endfor
      endif

      ## get the silhouette values for every clustering
      for iter = 1 : length (this.InspectedK)
        ## do it only for the specified K values
        if (any (this.InspectedK(iter) == K))
          ## Custom call to silhouette to avoid plotting any figures
          this.ClusterSilhouettes{iter} = silhouette (UsableX, ...
                                          this.ClusteringSolutions(:, iter), ...
                                          "sqeuclidean", "DoNotPlot");
          if (strcmpi (this.ClusterPriors, "empirical"))
            this.CriterionValues(iter) = mean (this.ClusterSilhouettes{iter});
          else
            ## equal
            this.CriterionValues(iter) = 0;
            si = this.ClusterSilhouettes{iter};
            for k = 1 : this.InspectedK(iter)
              this.CriterionValues(iter) += mean (si(find ...
                                     (this.ClusteringSolutions(:, iter) == k)));
            endfor
            this.CriterionValues(iter) /= this.InspectedK(iter);
          endif
        endif
      endfor

      [~, this.OptimalIndex] = max (this.CriterionValues);
      this.OptimalK = this.InspectedK(this.OptimalIndex(1));
      this.OptimalY = this.ClusteringSolutions(:, this.OptimalIndex(1));
    endfunction
  endmethods
endclassdef

%!test
%! load fisheriris
%! eva = evalclusters (meas, "kmeans", "silhouette", "KList", [1:6]);
%! assert (class (eva), "SilhouetteEvaluation");
