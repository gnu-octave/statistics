## Copyright (C) 2021 Stefano Guidoni <ilguido@users.sf.net>
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
## @deftypefn {Function File} {@var{eva} =} evalclusters (@var{x}, @var{clust}, @var{criterion})
## @deftypefnx {Function File} {@var{eva} =} evalclusters (@dots{}, @qcode{Name}, @qcode{Value})
##
## Create a clustering evaluation object to find the optimal number of clusters.
##
## @code{evalclusters} creates a clustering evaluation object to evaluate the
## optimal number of clusters for data @var{x}, using criterion @var{criterion}.
## The input data @var{x} is a matrix with @code{n} observations of @code{p}
## variables.
## The evaluation criterion @var{criterion} is one of the following:
## @table @code
## @item @qcode{CalinskiHarabasz} 
## to create a @code{CalinskiHarabaszEvaluation} object.
##
## @item @qcode{DaviesBouldin} 
## to create a @code{DaviesBouldinEvaluation} object.
##
## @item @qcode{gap} 
## to create a @code{GapEvaluation} object.
##
## @item @qcode{silhouette} 
## to create a @code{SilhouetteEvaluation} object.
##
## @end table
## The clustering algorithm @var{clust} is one of the following:
## @table @code
## @item @qcode{kmeans} 
## to cluster the data using @code{kmeans} with @code{EmptyAction} set to
## @code{singleton} and @code{Replicates} set to 5.
##
## @item @qcode{linkage} 
## to cluster the data using @code{clusterdata} with @code{linkage} set to
## @code{Ward}.
##
## @item @qcode{gmdistribution} 
## to cluster the data using @code{fitgmdist} with @code{SharedCov} set to
## @code{true} and @code{Replicates} set to 5.
##
## @end table
## If the @var{criterion} is @code{CalinskiHarabasz}, @code{DaviesBouldin}, or
## @code{silhouette}, @var{clust} can also be a function handle to a function
## of the form @code{c = clust(x, k)}, where @var{x} is the input data,
## @var{k} the number of clusters to evaluate and @var{c} the clustering result.
## The clustering result can be either an array of size @code{n} with @code{k}
## different integer values, or a matrix of size @code{n} by @code{k} with a
## likelihood value assigned to each one of the @code{n} observations for each
## one of the @var{k} clusters.  In the latter case, each observation is assigned
## to the cluster with the higher value.
## If the @var{criterion} is @code{CalinskiHarabasz}, @code{DaviesBouldin}, or
## @code{silhouette}, @var{clust} can also be a matrix of size @code{n} by 
## @code{k}, where @code{k} is the number of proposed clustering solutions, so
## that each column of @var{clust} is a clustering solution.
##
## In addition to the obligatory @var{x}, @var{clust} and @var{criterion} inputs
## there is a number of optional arguments, specified as pairs of @code{Name} 
## and @code{Value} options.  The known @code{Name} arguments are:
## @table @code
## @item @qcode{KList} 
## a vector of positive integer numbers, that is the cluster sizes to evaluate.
## This option is necessary, unless @var{clust} is a matrix of proposed 
## clustering solutions.
##
## @item @qcode{Distance} 
## a distance metric as accepted by the chosen @var{clust}.  It can be the
## name of the distance metric as a string or a function handle.  When
## @var{criterion} is @code{silhouette}, it can be a vector as created by 
## function @code{pdist}.  Valid distance metric strings are: @code{sqEuclidean}
## (default), @code{Euclidean}, @code{cityblock}, @code{cosine}, 
## @code{correlation}, @code{Hamming}, @code{Jaccard}.
## Only used by @code{silhouette} and @code{gap} evaluation.
##
## @item @qcode{ClusterPriors} 
## the prior probabilities of each cluster, which can be either @code{empirical}
## (default), or @code{equal}.  When @code{empirical} the silhouette value is
## the average of the silhouette values of all points; when @code{equal} the
## silhouette value is the average of the average silhouette value of each
## cluster.  Only used by @code{silhouette} evaluation.
##
## @item @qcode{B} 
## the number of reference datasets generated from the reference distribution.
## Only used by @code{gap} evaluation.
##
## @item @qcode{ReferenceDistribution} 
## the reference distribution used to create the reference data.  It can be
## @code{PCA} (default) for a distribution based on the principal components of
## @var{X}, or @code{uniform} for a uniform distribution based on the range of
## the observed data.  @code{PCA} is currently not implemented.
## Only used by @code{gap} evaluation.
##
## @item @qcode{SearchMethod} 
## the method for selecting the optimal value with a @code{gap} evaluation.  It 
## can be either @code{globalMaxSE} (default) for selecting the smallest number
## of clusters which is inside the standard error of the maximum gap value, or
## @code{firstMaxSE} for selecting the first number of clusters which is inside
## the standard error of the following cluster number.
## Only used by @code{gap} evaluation.
##
## @end table
##
## Output @var{eva} is a clustering evaluation object.
##
## @end deftypefn
##
## @seealso{CalinskiHarabaszEvaluation, DaviesBouldinEvaluation, GapEvaluation, 
## SilhouetteEvaluation}

function cc = evalclusters (x, clust, criterion, varargin)

  ## input check
  if (nargin < 3)
    print_usage ();
  endif

  ## parsing input data
  if ((! ismatrix (x)) || (! isnumeric (x)))
    error ("evalclusters: 'x' must be a numeric matrix");
  endif

  ## useful values for input check
  n = rows (x);
  p = columns (x);

  ## parsing the clustering algorithm
  if (ischar (clust))
    clust = lower (clust);
    if (! any (strcmpi (clust, {"kmeans", "linkage", "gmdistribution"})))
      error ("evalclusters: unknown clustering algorithm '%s'", clust);
    endif
  elseif (! isscalar (clust))
    if ((! isnumeric (clust))  || (length (size (clust)) != 2) || ...
        (rows (clust) != n))
      error ("evalclusters: invalid matrix of clustering solutions");
    endif
  elseif (! isa (clust, "function_handle"))
    error ("evalclusters: invalid argument for 'clust'");
  endif

  ## parsing the criterion parameter
  ## we check the rest later, as the check depends on the chosen criterion
  if (! ischar (criterion))
    error ("evalclusters: invalid criterion, it must be a string");
  else
    criterion = lower (criterion);
    if (! any (strcmpi (criterion, {"calinskiharabasz", "daviesbouldin", ...
                                   "silhouette", "gap"})))
      error ("evalclusters: unknown criterion '%s'", criterion);
    endif
  endif

  ## some default value
  klist = [];
  distance = "sqeuclidean";
  clusterpriors = "empirical";
  b = 100;
  referencedistribution = "pca";
  searchmethod = "globalmaxse";

  ## parse the name/value pairs
  pair_index = 1;
  while (pair_index < (nargin - 3))
    ## type check
    if (! ischar (varargin{pair_index}))
      error ("evalclusters: invalid property, string expected");
    endif
    
    ## now parse the parameter
    switch (lower (varargin{pair_index}))
      case "klist"
        ## klist must be an array of positive interger numbers;
        ## there is a special case when it can be empty, but that is not the
        ## suggested way to use it (it is better to omit it instead)
        if (isempty (varargin{pair_index + 1}))
          if (ischar (clust) || isa (clust, "function_handle"))
            error (["evalclusters: 'KList' can be empty only when 'clust' "...
                    "is a matrix"]);
          endif
        elseif ((! isnumeric (varargin{pair_index + 1})) || ...
            (! isvector (varargin{pair_index + 1})) || ...
            any (find (varargin{pair_index + 1} < 1)) || ...
            any (floor (varargin{pair_index + 1}) != varargin{pair_index + 1}))
          error ("evalclusters: 'KList' must be an array of positive integers");
        endif
        klist = varargin{pair_index + 1};

      case "distance"
        ## used by silhouette and gap
        if (! (strcmpi (criterion, "silhouette") || strcmpi (criterion, "gap")))
          error (["evalclusters: distance metric cannot be used with '%s'"...
                  " criterion"], criterion);
        endif
        if (ischar (varargin{pair_index + 1}))
          if (! any (strcmpi (varargin{pair_index + 1}, ...
                    {"sqeuclidean", "euclidean", "cityblock", "cosine", ...
                     "correlation", "hamming", "jaccard"})))
            error ("evalclusters: unknown distance criterion '%s'", ...
                    varargin{pair_index + 1});
          endif
        elseif (! isa (varargin{pair_index + 1}, "function_handle") ||
                ! ((isvector (varargin{pair_index + 1}) && ...
                    isnumeric (varargin{pair_index + 1}))))
          error ("evalclusters: invalid distance metric");
        endif
        distance = varargin{pair_index + 1};

      case "clusterpriors"
        ## used by silhouette evaluation
        if (! strcmpi (criterion, "silhouette"))
          error (["evalclusters: cluster prior probabilities cannot be used "...
                  "with '%s' criterion"], criterion);
        endif
        if (any (strcmpi (varargin{pair_index + 1}, {"empirical", "equal"})))
          clusterpriors = lower (varargin{pair_index + 1});
        else
          error ("evalclusters: invalid cluster prior probabilities value");
        endif

      case "b"
        ## used by gap evaluation
        if (! isnumeric (varargin{pair_index + 1}) || ...
            ! isscalar (varargin{pair_index + 1}) || ...
            varargin{pair_index + 1} != floor (varargin{pair_index + 1}) || ...
            varargin{pair_index + 1} < 1)
          error ("evalclusters: b must a be positive integer number");
        endif
        b = varargin{pair_index + 1};

      case "referencedistribution"
        ## used by gap evaluation
        if (! ischar (varargin{pair_index + 1}) || any (strcmpi ...
            (varargin{pair_index + 1}, {"pca", "uniform"})))
          error (["evalclusters: the reference distribution must be either" ...
                  "'PCA' or 'uniform'"]);
        endif
        referencedistribution = lower (varargin{pair_index + 1});

      case "searchmethod"
        ## used by gap evaluation
        if (! ischar (varargin{pair_index + 1}) || any (strcmpi ...
            (varargin{pair_index + 1}, {"globalmaxse", "uniform"})))
          error (["evalclusters: the search method must be either" ...
                  "'globalMaxSE' or 'firstmaxse'"]);
        endif
        searchmethod = lower (varargin{pair_index + 1});

      otherwise
        error ("evalclusters: unknown property %s", varargin{pair_index});

    endswitch

    pair_index += 2;
  endwhile 

  ## check if there are parameters without a value or a name left
  if (nargin - 2 - pair_index)
    if (ischar (varargin{pair_index}))
      error ("evalclusters: invalid parameter '%s'", varargin{pair_index});
    else 
      error ("evalclusters: invalid parameter '%d'", varargin{pair_index});
    endif
  endif 

  ## another check on klist
  if (isempty (klist) && (ischar (clust) || isa (clust, "function_handle")))
    error (["evalclusters: 'KList' can be empty only when 'clust' ", ...
            "is a matrix"]);
  endif

  ## main
  switch (lower (criterion))
    case "calinskiharabasz"
      ## further compatibility checks between the chosen parameters are 
      ## delegated to the class constructor
      if (isempty (klist))
        klist = 1 : columns (clust);
      endif
      cc = CalinskiHarabaszEvaluation (x, clust, klist);

    case "daviesbouldin"
      ## further compatibility checks between the chosen parameters are 
      ## delegated to the class constructor
      if (isempty (klist))
        klist = 1 : columns (clust);
      endif
      cc = DaviesBouldinEvaluation (x, clust, klist);

    case "silhouette"
      ## further compatibility checks between the chosen parameters are 
      ## delegated to the class constructor
      if (isempty (klist))
        klist = 1 : columns (clust);
      endif
      cc = SilhouetteEvaluation (x, clust, klist, distance, clusterpriors);

    case "gap"
      ## gap cannot be used with a pre-computed solution, i.e. a matrix for
      ## 'clust', and klist must be specified
      if (isnumeric (clust))
        error (["evalclusters: 'clust' must be a clustering algorithm when "... 
                "using the gap criterion"]);
      endif
      if (isempty (klist))
        error (["evalclusters: 'klist' cannot be empty when using the gap " ...
                "criterion"]);
      endif
      cc = GapEvaluation (x, clust, klist, b, distance, ...
                          referencedistribution, searchmethod);

    otherwise
      error ("evalclusters: invalid criterion '%s'", criterion);

  endswitch

endfunction


## input tests
%!error evalclusters ()
%!error evalclusters ([1 1;0 1])
%!error evalclusters ([1 1;0 1], "kmeans")
%!error <'x' must be a numeric*> evalclusters ("abc", "kmeans", "gap")
%!error <unknown clustering*> evalclusters ([1 1;0 1], "xxx", "gap")
%!error <invalid matrix*> evalclusters ([1 1;0 1], [1 2], "gap")
%!error <invalid argument*> evalclusters ([1 1;0 1], 1.2, "gap")
%!error <invalid criterion*> evalclusters ([1 1;0 1], [1; 2], 123)
%!error <unknown criterion*> evalclusters ([1 1;0 1], [1; 2], "xxx")
%!error <'KList' can be empty*> evalclusters ([1 1;0 1], "kmeans", "gap")
%!error <invalid parameter*> evalclusters ([1 1;0 1], [1; 2], "gap", 1)
%!error <invalid property*> evalclusters ([1 1;0 1], [1; 2], "gap", 1, 1)
%!error <unknown property*> evalclusters ([1 1;0 1], [1; 2], "gap", "xxx", 1)
%!error <'KList'*> evalclusters ([1 1;0 1], [1; 2], "gap", "KList", [-1 0])
%!error <'KList'*> evalclusters ([1 1;0 1], [1; 2], "gap", "KList", [1 .5])
%!error <'KList'*> evalclusters ([1 1;0 1], [1; 2], "gap", "KList", [1 1; 1 1])
%!error <unknown distance*> evalclusters ([1 1;0 1], [1; 2], "gap", ...
%!                                        "distance", "a")
%!error <distance metric*> evalclusters ([1 1;0 1], [1; 2], "daviesbouldin", ...
%!                                       "distance", "a")
%!error <cluster prior*> evalclusters ([1 1;0 1], [1; 2], "gap", ...
%!                                     "clusterpriors", "equal")
%!error <invalid cluster prior*> evalclusters ([1 1;0 1], [1; 2], ...
%!                                         "silhouette", "clusterpriors", "xxx")
%!error <'clust' must be a clustering*> evalclusters ([1 1;0 1], [1; 2], "gap")

## demonstration
%!demo
%! load fisheriris;
%! eva = evalclusters(meas, "kmeans", "calinskiharabasz", "KList", [1:6])
