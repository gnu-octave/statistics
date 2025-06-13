## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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

classdef hnswSearcher
## -*- texinfo -*-
## @deftp {Class} hnswSearcher
##
## Hierarchical Navigable Small World (HNSW) nearest neighbor searcher class.
##
## The @code{hnswSearcher} class implements the HNSW algorithm for efficient
## nearest neighbor queries.  It stores training data and supports various
## distance metrics for performing searches.  The HNSW algorithm builds a
## multilayer graph structure that enables fast approximate nearest neighbor
## searches by navigating through the graph.  It facilitates a nearest neighbor
## search using @code{knnsearch} or a radius search using @code{rangesearch}.
##
## You can either use the @code{hnswSearcher} class constructor or the
## @code{createns} function to create an @qcode{hnswSearcher} object.
##
## @seealso{createns, ExhaustiveSearcher, KDTreeSearcher, knnsearch,
## rangesearch, pdist2}
## @end deftp

	properties (SetAccess = private)
    ## -*- texinfo -*-
    ## @deftp {Property} X
    ##
    ## Training data, specified as an @math{NxP} numeric matrix where each row
    ## is an observation and each column is a feature.  This property is private
    ## and cannot be modified after object creation.
    ##
    ## @end deftp
    X = []

    ## -*- texinfo -*-
    ## @deftp {Property} HNSWGraph
    ##
    ## The HNSW graph structure built from the training data.  This property is
    ## private and cannot be modified after object creation.
    ##
    ## @end deftp
    HNSWGraph
  endproperties

  properties
    ## -*- texinfo -*-
    ## @deftp {Property} Distance
    ##
    ## Distance metric used for searches, specified as a character vector (e.g.,
    ## @qcode{"euclidean"}, @qcode{"minkowski"}). Default is
    ## @qcode{"euclidean"}. Supported metrics align with those in @code{pdist2}.
    ##
    ## @end deftp
    Distance = 'euclidean'

    ## -*- texinfo -*-
    ## @deftp {Property} DistParameter
    ##
    ## Parameter for the distance metric, with type and value depending on
    ## @qcode{Distance}:
    ##
    ## @itemize
    ## @item For @qcode{"minkowski"}, a positive scalar exponent (default 2).
    ## @item For @qcode{"seuclidean"}, a nonnegative vector of scaling factors
    ## matching the number of columns in @qcode{X} (default is standard
    ## deviation of @qcode{X}).
    ## @item For @qcode{"mahalanobis"}, a positive definite covariance matrix
    ## matching the dimensions of @qcode{X} (default is @code{cov (@var{X})}).
    ## @item Empty for other metrics.
    ## @end itemize
    ##
    ## @end deftp
    DistParameter = []

    ## -*- texinfo -*-
    ## @deftp {Property} M
    ##
    ## Maximum number of neighbors per node in the HNSW graph. Affects graph
    ## connectivity and search accuracy. Default is 16.
    ##
    ## @end deftp
    M = 16

    ## -*- texinfo -*-
    ## @deftp {Property} efConstruction
    ##
    ## Size of the dynamic candidate list during graph construction. Higher
    ## values improve accuracy at the cost of construction time. Default is 100.
    ##
    ## @end deftp
    efConstruction = 100

    ## -*- texinfo -*-
    ## @deftp {Property} efSearch
    ##
    ## Size of the dynamic candidate list during search. Higher values improve
    ## accuracy at the cost of search time. Default is 50.
    ##
    ## @end deftp
    efSearch = 50
  endproperties