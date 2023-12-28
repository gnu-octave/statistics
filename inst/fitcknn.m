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
## @deftypefn  {statistics} {@var{obj} =} fitcknn
## @deftypefnx {statistics} {@var{obj} =} fitcknn (@var{X}, @var{Y})
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
## @code{@var{obj} = fitcknn ()} returns an empty k-Nearest Neighbor
## classification model and produces a warning.
##
## @code{@var{obj} = fitcknn (@dots{}, @var{name}, @var{value})} returns a
## k-Nearest Neighbor classification model with additional options specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab @var{Name} @tab @var{Value}
##
## @item @tab @qcode{"K"} @tab is the number of nearest neighbors to be found
## in the kNN search.  It must be a positive integer value and by default it is
## 1. @var{"K"} can be changed by using @code{@var{obj.K} = 10}.
##
## @item @tab @qcode{"weights"} @tab is a @math{Nx1} numeric non-negative matrix
## of the observational weights, each row in @var{weights} corresponds to the
## row in @var{Y} and indicates the relative importance or weight to be
## considered in calculating the Nearest-neighbour, negative values are removed
##  before calculations if weights are specified. this propery cannot be
## changed via object. default value @qcode{@var{weight} = ones(rows(Y),1)}.
##
## @item @tab @qcode{"cost"} @tab is a @math{NxR} numeric matrix containing
## misclassification cost for the corresponding instances in @var{X} where
## @var{R} is the number of unique categories in @var{Y}. If an instance is
## correctly classified into its category the cost is calculated to be 1, If
## not then 0. cost matrix can be altered use @code{@var{obj.cost} = somecost}.
## default value @qcode{@var{cost} = ones(rows(X),numel(unique(Y)))}.
##
## @item @tab @qcode{"exponent"} @tab is the Minkowski distance exponent and
## it must be a positive scalar.  This argument is only valid when the selected
## distance metric is @qcode{"minkowski"}.  By default it is 2.
##
## @item @tab @qcode{"scale"} @tab is the scale parameter for the standardized
## Euclidean distance and it must be a nonnegative numeric vector of equal
## length to the number of columns in @var{X}.  This argument is only valid when
## the selected distance metric is @qcode{"seuclidean"}, in which case each
## coordinate of @var{X} is scaled by the corresponding element of
## @qcode{"scale"}, as is each query point in @var{Y}.  By default, the scale
## parameter is the standard deviation of each coordinate in @var{X}.
##
## @item @tab @qcode{"Prior"} @tab is the matrix of prior probablities for each
## unique class in @var{Y}. @var{Prior} can be changed via object using
## @code{@var{obj.Prior} = somePrior}
##
## @item @tab @qcode{"cov"} @tab is the covariance matrix for computing the
## mahalanobis distance and it must be a positive definite matrix matching the
## the number of columns in @var{X}.  This argument is only valid when the
## selected distance metric is @qcode{"mahalanobis"}.
##
## @item @tab @qcode{"BucketSize"} @tab is the maximum number of data points in
## the leaf node of the Kd-tree and it must be a positive integer.  This
## argument is only valid when the selected search method is @qcode{"kdtree"}.
##
## @item @tab @qcode{"Distance"} @tab is the distance metric used by
## @code{knnsearch} as specified below:
## @end multitable
##
## @multitable @columnfractions 0.1 0.25 0.65
## @item @tab @qcode{"euclidean"} @tab Euclidean distance.
## @item @tab @qcode{"seuclidean"} @tab standardized Euclidean distance.  Each
## coordinate difference between the rows in @var{X} and the query matrix
## @var{Y} is scaled by dividing by the corresponding element of the standard
## deviation computed from @var{X}.  To specify a different scaling, use the
## @qcode{"scale"} name-value argument.
## @item @tab @qcode{"cityblock"} @tab City block distance.
## @item @tab @qcode{"chebychev"} @tab Chebychev distance (maximum coordinate
## difference).
## @item @tab @qcode{"minkowski"} @tab Minkowski distance.  The default exponent
## is 2.  To specify a different exponent, use the @qcode{"P"} name-value
## argument.
## @item @tab @qcode{"mahalanobis"} @tab Mahalanobis distance, computed using a
## positive definite covariance matrix.  To change the value of the covariance
## matrix, use the @qcode{"cov"} name-value argument.
## @item @tab @qcode{"cosine"} @tab Cosine distance.
## @item @tab @qcode{"correlation"} @tab One minus the sample linear correlation
## between observations (treated as sequences of values).
## @item @tab @qcode{"spearman"} @tab One minus the sample Spearman's rank
## correlation between observations (treated as sequences of values).
## @item @tab @qcode{"hamming"} @tab Hamming distance, which is the percentage
## of coordinates that differ.
## @item @tab @qcode{"jaccard"} @tab One minus the Jaccard coefficient, which is
## the percentage of nonzero coordinates that differ.
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"NSMethod"} @tab is the nearest neighbor search method used
## by @code{knnsearch} as specified below.
## @end multitable
##
## @multitable @columnfractions 0.1 0.25 0.65
## @item @tab @qcode{"kdtree"} @tab Creates and uses a Kd-tree to find nearest
## neighbors.  @qcode{"kdtree"} is the default value when the number of columns
## in @var{X} is less than or equal to 10, @var{X} is not sparse, and the
## distance metric is @qcode{"euclidean"}, @qcode{"cityblock"},
## @qcode{"chebychev"}, or @qcode{"minkowski"}.  Otherwise, the default value is
## @qcode{"exhaustive"}.  This argument is only valid when the distance metric
## is one of the four aforementioned metrics.
## @item @tab @qcode{"exhaustive"} @tab Uses the exhaustive search algorithm by
## computing the distance values from all the points in @var{X} to each point in
## @var{Y}.
## @end multitable
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab @qcode{"IncludeTies"} @tab is a boolean flag to indicate if the
## returned values should contain the indices that have same distance as the
## @math{K^th} neighbor.  When @qcode{false}, @code{knnsearch} chooses the
## observation with the smallest index among the observations that have the same
## distance from a query point.  When @qcode{true}, @code{knnsearch} includes
## all nearest neighbors whose distances are equal to the @math{K^th} smallest
## distance in the output arguments. To specify @math{K}, use the @qcode{"K"}
## name-value pair argument.
##
## @item @tab @qcode{"breakties"} @tab is the method to break ties in kNN
## algorithm while deciding the label for an observation.
##
## @item @tab @qcode{"classnames"} @tab is a vector of unique classes in
## training data @var{Y}.If not given as name-value pair, the function
## calculates and stores the unique values in @var{Y}. @var{"classnames"} cannot
## be changed via @var{obj}.
##
## @item @tab @qcode{"NosClasses"} @tab is a numeric value that stores the
## number of unique classes in @var{Y}.If not given as name-value pair, the
## function calculates and stores the number of unique values in @var{Y}.
## @var{"NosClasses"} cannot be changed via @var{obj}.
##
## @item @tab @qcode{"Xclass"} @tab must be a @math{MxP} numeric matrix of
## query/new points that are to be classified into the class labels.
## @var{Xclass} must have same numbers of columns as @var{X}. @var{Xclass}
## can be changed via object @var{obj} using @code{@var{obj.Xclass} = newXclass}.
## to predict labels for new set of values while keeping the rest of the model
## same with all its properties.
##
## @end multitable
##
## @seealso{knnpredict, @ClassificationKNN/predict}
## @end deftypefn

function obj = fitcknn (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2 && nargin != 0)
    error ("fitcknn: Too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitcknn: Name-Value arguments must be in pairs.");
  endif

  ## Create object
  if (nargin == 0)
    ## Return empty object and issue a warning
    obj = ClassificationKNN ();
    warning ("fitcknn: No Argument given, Created Object will be empty.");
  else
    ## Check predictor data and labels have equal rows
    if (rows (X) != rows (Y))
      error ("fitcknn: number of rows in X and Y must be equal.");
    endif
    ## Parse arguments to class def function
    obj = ClassificationKNN (X, Y, varargin{:});
  endif

endfunction


%!demo
%! ## Find 10 nearest neighbour of a point using different distance metrics
%! ## and compare the results by plotting
%!
%! load fisheriris
%! x = meas;
%! y = species;
%! xnew = [5, 3, 5, 1.45];
%!
%! ## create an object
%! a = fitcknn (x, y, "Xclass" , xnew, "k", 5)
%!
%! ## predict labels for points in xnew
%! predict (a)
%!
%! ## change properties keeping training data and predict again
%! a.distance = "hamming";
%! a.k = 10;
%! predict (a)

## Test Output
%!test
%! warning("off");
%! a = fitcknn ();
%! assert (class (a), "ClassificationKNN");
%! assert ({a.Breakties, a.Includeties, a.NN, a.NSmethod}, {[],[],[],[]})
%! assert ({a.NosClasses, a.NumObsv, a.Scale, a.Score, a.cost}, {[],[],[],[],[]})
%! assert ({a.X, a.Xclass, a.Y, a.bucketsize, a.classNames}, {[],[],[],[],[]})
%! assert ({a.cov, a.distance, a.P, a.k, a.label}, {[],[],[],[],[]})
%! assert ({a.prior, a.standardize, a.weights}, {[],[],[]})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ["a";"a";"b";"b"];
%! a = fitcknn (x, y);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ["a";"a";"b";"b"], 1})
%! assert ({a.NSmethod, a.distance, a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ['a';'a';'b';'b'];
%! k = 10;
%! a = fitcknn (x, y, "K" ,k);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ["a";"a";"b";"b"], 10})
%! assert ({a.NSmethod, a.distance, a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ["a";"a";"b";"b"];
%! weights = ones (4,1);
%! a = fitcknn (x, y, "weights" , weights);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ["a";"a";"b";"b"], 1})
%! assert (a.weights, ones (4,1))
%! assert ({a.NSmethod, a.distance, a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ["a";"a";"b";"b"];
%! a = fitcknn (x, y, "P" , 10);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ["a";"a";"b";"b"], 1})
%! assert (a.P, 10)
%! assert ({a.NSmethod, a.distance, a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ["a";"a";"b";"b"];
%! cov = rand (4,1);
%! a = fitcknn (x, y, "cov" , cov, "distance", "mahalanobis");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ["a";"a";"b";"b"], 1})
%! assert (a.cov, cov)
%! assert ({a.NSmethod, a.distance, a.bucketsize},{"exhaustive","mahalanobis",50})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ["a";"a";"b";"b"];
%! a = fitcknn (x, y, "bucketsize" , 20, "distance", "mahalanobis");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ["a";"a";"b";"b"], 1})
%! assert ({a.NSmethod,a.distance,a.bucketsize},{"exhaustive","mahalanobis",20})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ["a";"a";"b";"b"];
%! a = fitcknn (x, y, "standardize", true);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ["a";"a";"b";"b"], 1})
%! assert (a.standardize, true);
%! assert ({a.NSmethod,a.distance,a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ["a";"a";"b";"b"'];
%! a = fitcknn (x, y, "includeties", true);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ["a";"a";"b";"b"], 1})
%! assert (a.Includeties, true);
%! assert ({a.NSmethod,a.distance,a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ["a";"a";"b";"b"];
%! cost = ones (4,2);
%! a = fitcknn (x, y, "cost", cost );
%! assert (class (a), "ClassificationKNN")
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ["a";"a";"b";"b"], 1})
%! assert (a.cost, [1,1;1,1;1,1;1,1])
%! assert ({a.NSmethod,a.distance,a.bucketsize},{"exhaustive","euclidean",50})

%!test
%! x = [1,2,3;4,5,6;7,8,9;3,2,1];
%! y = ["a";"a";"b";"b"];
%! scale = [1,2,3,4];
%! a = fitcknn (x, y, "scale", scale );
%! assert (class (a), "ClassificationKNN")
%! assert ({a.X, a.Y, a.k},{[1,2,3;4,5,6;7,8,9;3,2,1], ["a";"a";"b";"b"], 1})
%! assert (a.Scale, [1,2,3,4])
%! assert ({a.NSmethod,a.distance,a.bucketsize},{"exhaustive","euclidean",50})

## Test input validation
%!error<fitcknn: Too few arguments.>fitcknn (ones (4,1));
%!error<fitcknn: Name-Value arguments must be in pairs.>
%! fitcknn (ones (4,2), ones (4, 1), "K");
%!error<fitcknn: number of rows in X and Y must be equal.>
%! fitcknn (ones (4,2), ones (3, 1));
%!error<fitcknn: number of rows in X and Y must be equal.>
%! fitcknn (ones (4,2), ones (3, 1), "K", 2);
%!warning<fitcknn: No Argument given, Created Object will be empty.> fitcknn ();
