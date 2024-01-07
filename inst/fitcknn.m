## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{obj} =} fitcknn (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} fitcknn (@dots{}, @var{name}, @var{value})
##
## Fit a k-Nearest Neighbor classification model.
##
## @code{@var{obj} = fitcknn (@var{X}, @var{Y})} returns a k-Nearest Neighbor
## classification model, @var{obj}, with @var{X} being the predictor data,
## and @var{Y} the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the kNN model.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain any type of
## categorical data. @var{Y} must have same numbers of Rows as @var{X}.
## @item
## @end itemize
##
## @code{@var{obj} = fitcknn (@dots{}, @var{name}, @var{value})} returns a
## k-Nearest Neighbor classification model with additional options specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @qcode{"PredictorNames"} @tab @tab A cell array of character vectors
## specifying the predictor variable names.  The variable names are assumed to
## be in the same order as they appear in the training data @var{X}.
##
## @item @qcode{"ResponseName"} @tab @tab A character vector specifying the name
## of the response variable.
##
## @item @qcode{"ClassNames"} @tab @tab A cell array of character vectors
## specifying the names of the classes in the training data @var{Y}.
##
## @item @qcode{"BreakTies"} @tab @tab Tie-breaking algorithm used by predict
## when multiple classes have the same smallest cost. By default, ties occur
## when multiple classes have the same number of nearest points among the
## @math{k} nearest neighbors. The available options are specified by the
## following character arrays:
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"smallest"} @tab This is the default and it favors the
## class with the smallest index among the tied groups, i.e. the one that
## appears first in the training labelled data.
## @item @tab @qcode{"nearest"} @tab This favors the class with the nearest
## neighbor among the tied groups, i.e. the class with the closest member point
## according to the distance metric used.
## @item @tab @qcode{"nearest"} @tab This randomly picks one class among the
## tied groups.
## @end multitable
##
## @multitable @columnfractions 0.18 0.02 0.8
## @item @qcode{"BucketSize"} @tab @tab The maximum number of data points in the
## leaf node of the Kd-tree and it must be a positive integer.  By default, it
## is 50. This argument is meaningful only when the selected search method is
## @qcode{"kdtree"}.
##
## @item @qcode{"Cost"} @tab @tab A @math{NxR} numeric matrix containing
## misclassification cost for the corresponding instances in @var{X} where
## @math{R} is the number of unique categories in @var{Y}.  If an instance is
## correctly classified into its category the cost is calculated to be 1, If
## not then 0. cost matrix can be altered use @code{@var{obj.cost} = somecost}.
## default value @qcode{@var{cost} = ones(rows(X),numel(unique(Y)))}.
##
## @item @qcode{obj.Prior} @tab @tab A numeric vector specifying the prior
## probabilities for each class.  The order of the elements in @qcode{Prior}
## corresponds to the order of the classes in @qcode{ClassNames}.
##
## @item @qcode{"NumNeighbors"} @tab @tab A positive integer value specifying
## the number of nearest neighbors to be found in the kNN search.  By default,
## it is 1.
##
## @item @qcode{"Exponent"} @tab @tab A positive scalar (usually an integer)
## specifying the Minkowski distance exponent.  This argument is only valid when
## the selected distance metric is @qcode{"minkowski"}.  By default it is 2.
##
## @item @qcode{"Scale"} @tab @tab A nonnegative numeric vector specifying the
## scale parameters for the standardized Euclidean distance.  The vector length
## must be equal to the number of columns in @var{X}.  This argument is only
## valid when the selected distance metric is @qcode{"seuclidean"}, in which
## case each coordinate of @var{X} is scaled by the corresponding element of
## @qcode{"scale"}, as is each query point in @var{Y}.  By default, the scale
## parameter is the standard deviation of each coordinate in @var{X}.  If a
## variable in @var{X} is constant, i.e. zero variance, this value is forced
## to 1 to avoid division by zero.  This is the equivalent of this variable not
## being standardized.
##
## @item @qcode{"Cov"} @tab @tab A square matrix with the same number of columns
## as @var{X} specifying the covariance matrix for computing the mahalanobis
## distance.  This must be a positive definite matrix matching.  This argument
## is only valid when the selected distance metric is @qcode{"mahalanobis"}.
##
## @item @qcode{"Distance"} @tab @tab is the distance metric used by
## @code{knnsearch} as specified below:
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"euclidean"} @tab Euclidean distance.
## @item @tab @qcode{"seuclidean"} @tab standardized Euclidean distance.  Each
## coordinate difference between the rows in @var{X} and the query matrix
## @var{Y} is scaled by dividing by the corresponding element of the standard
## deviation computed from @var{X}.  To specify a different scaling, use the
## @qcode{"Scale"} name-value argument.
## @item @tab @qcode{"cityblock"} @tab City block distance.
## @item @tab @qcode{"chebychev"} @tab Chebychev distance (maximum coordinate
## difference).
## @item @tab @qcode{"minkowski"} @tab Minkowski distance.  The default exponent
## is 2.  To specify a different exponent, use the @qcode{"P"} name-value
## argument.
## @item @tab @qcode{"mahalanobis"} @tab Mahalanobis distance, computed using a
## positive definite covariance matrix.  To change the value of the covariance
## matrix, use the @qcode{"Cov"} name-value argument.
## @item @tab @qcode{"cosine"} @tab Cosine distance.
## @item @tab @qcode{"correlation"} @tab One minus the sample linear correlation
## between observations (treated as sequences of values).
## @item @tab @qcode{"spearman"} @tab One minus the sample Spearman's rank
## correlation between observations (treated as sequences of values).
## @item @tab @qcode{"hamming"} @tab Hamming distance, which is the percentage
## of coordinates that differ.
## @item @tab @qcode{"jaccard"} @tab One minus the Jaccard coefficient, which is
## the percentage of nonzero coordinates that differ.
## @item @tab @var{@@distfun} @tab Custom distance function handle.  A distance
## function of the form @code{function @var{D2} = distfun (@var{XI}, @var{YI})},
## where @var{XI} is a @math{1xP} vector containing a single observation in
## @math{P}-dimensional space, @var{YI} is an @math{NxP} matrix containing an
## arbitrary number of observations in the same @math{P}-dimensional space, and
## @var{D2} is an @math{NxP} vector of distances, where @qcode{(@var{D2}k)} is
## the distance between observations @var{XI} and @qcode{(@var{YI}k,:)}.
## @end multitable
##
## @multitable @columnfractions 0.18 0.02 0.8
## @item @qcode{"DistanceWeight"} @tab @tab A distance weighting function,
## specified either as a function handle, which accepts a matrix of nonnegative
## distances and returns a matrix the same size containing nonnegative distance
## weights, or one of the following values: @qcode{"equal"}, which corresponds
## to no weighting; @qcode{"inverse"}, which corresponds to a weight equal to
## @math{1/distance}; @qcode{"squaredinverse"}, which corresponds to a weight
## equal to @math{1/distance^2}.
##
## @item @qcode{"IncludeTies"} @tab @tab A boolean flag to indicate if the
## returned values should contain the indices that have same distance as the
## @math{K^th} neighbor.  When @qcode{false}, @code{knnsearch} chooses the
## observation with the smallest index among the observations that have the same
## distance from a query point.  When @qcode{true}, @code{knnsearch} includes
## all nearest neighbors whose distances are equal to the @math{K^th} smallest
## distance in the output arguments. To specify @math{K}, use the @qcode{"K"}
## name-value pair argument.
##
## @item @qcode{"NSMethod"} @tab @tab is the nearest neighbor search method used
## by @code{knnsearch} as specified below.
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"kdtree"} @tab Creates and uses a Kd-tree to find nearest
## neighbors.  @qcode{"kdtree"} is the default value when the number of columns
## in @var{X} is less than or equal to 10, @var{X} is not sparse, and the
## distance metric is @qcode{"euclidean"}, @qcode{"cityblock"},
## @qcode{"manhattan"}, @qcode{"chebychev"}, or @qcode{"minkowski"}.  Otherwise,
## the default value is @qcode{"exhaustive"}.  This argument is only valid when
## the distance metric is one of the four aforementioned metrics.
## @item @tab @qcode{"exhaustive"} @tab Uses the exhaustive search algorithm by
## computing the distance values from all the points in @var{X} to each point in
## @var{Y}.
## @end multitable
##
## @seealso{@@ClassificationKNN/predict, @@ClassificationKNN/ClassificationKNN,
## knnsearch, pdist2}
## @end deftypefn

function obj = fitcknn (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitcknn: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitcknn: Name-Value arguments must be in pairs.");
  endif

  ## Check predictor data and labels have equal rows
  if (rows (X) != rows (Y))
    error ("fitcknn: number of rows in X and Y must be equal.");
  endif
  ## Parse arguments to class def function
  obj = ClassificationKNN (X, Y, varargin{:});

endfunction


## Test Output
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcknn (x, y);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 1})
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcknn (x, y, "NSMethod", "exhaustive");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 1})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! k = 10;
%! a = fitcknn (x, y, "NumNeighbors" ,k);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 10})
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = ones (4, 11);
%! y = ["a"; "a"; "b"; "b"];
%! k = 10;
%! a = fitcknn (x, y, "NumNeighbors" ,k);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 10})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! k = 10;
%! a = fitcknn (x, y, "NumNeighbors" ,k, "NSMethod", "exhaustive");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 10})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! k = 10;
%! a = fitcknn (x, y, "NumNeighbors" ,k, "Distance", "hamming");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 10})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "hamming"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! weights = ones (4,1);
%! a = fitcknn (x, y, "Standardize", 1);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 1})
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.Standardize}, {true})
%! assert ({a.Sigma}, {std(x, [], 1)})
%! assert ({a.Mu}, {[3.75, 4.25, 4.75]})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! weights = ones (4,1);
%! a = fitcknn (x, y, "Standardize", false);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 1})
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.Standardize}, {false})
%! assert ({a.Sigma}, {[]})
%! assert ({a.Mu}, {[]})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! s = ones (1, 3);
%! a = fitcknn (x, y, "Scale" , s, "Distance", "seuclidean");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.DistParameter}, {s})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "seuclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! C = cov (x);
%! a = fitcknn (x, y, "Cov" , C, "distance", "mahalanobis");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.DistParameter}, {C})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "mahalanobis"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcknn (x, y, "Exponent" , 5, "Distance", "minkowski");
%! assert (class (a), "ClassificationKNN");
%! assert (a.DistParameter, 5)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "minkowski"})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcknn (x, y, "Exponent" , 5, "Distance", "minkowski", ...
%!                    "NSMethod", "exhaustive");
%! assert (class (a), "ClassificationKNN");
%! assert (a.DistParameter, 5)
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "minkowski"})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcknn (x, y, "BucketSize" , 20, "distance", "mahalanobis");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "mahalanobis"})
%! assert ({a.BucketSize}, {20})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcknn (x, y, "IncludeTies", true);
%! assert (class (a), "ClassificationKNN");
%! assert (a.IncludeTies, true);
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcknn (x, y);
%! assert (class (a), "ClassificationKNN");
%! assert (a.IncludeTies, false);
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcknn (x, y);
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, [0.5; 0.5])
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! prior = [0.5; 0.5];
%! a = fitcknn (x, y, "Prior", "empirical");
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, prior)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "a"; "b"];
%! prior = [0.75; 0.25];
%! a = fitcknn (x, y, "Prior", "empirical");
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, prior)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "a"; "b"];
%! prior = [0.5; 0.5];
%! a = fitcknn (x, y, "Prior", "uniform");
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, prior)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! cost = eye (2);
%! a = fitcknn (x, y, "Cost", cost);
%! assert (class (a), "ClassificationKNN")
%! assert (a.Cost, [1, 0; 0, 1])
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! cost = eye (2);
%! a = fitcknn (x, y, "Cost", cost, "Distance", "hamming" );
%! assert (class (a), "ClassificationKNN")
%! assert (a.Cost, [1, 0; 0, 1])
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "hamming"})
%! assert ({a.BucketSize}, {50})

## Test input validation
%!error<fitcknn: too few arguments.> fitcknn ()
%!error<fitcknn: too few arguments.> fitcknn (ones (4,1))
%!error<fitcknn: Name-Value arguments must be in pairs.>
%! fitcknn (ones (4,2), ones (4, 1), "K")
%!error<fitcknn: number of rows in X and Y must be equal.>
%! fitcknn (ones (4,2), ones (3, 1))
%!error<fitcknn: number of rows in X and Y must be equal.>
%! fitcknn (ones (4,2), ones (3, 1), "K", 2)
