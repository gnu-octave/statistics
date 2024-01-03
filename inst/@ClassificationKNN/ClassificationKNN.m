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

classdef ClassificationKNN
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{obj} =} ClassificationKNN (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} ClassificationKNN (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{ClassificationKNN} class object containing a k-Nearest
## Neighbor classification model.
##
## @code{@var{obj} = ClassificationKNN (@var{X}, @var{Y})} returns a
## ClassificationKNN object, with @var{X} as the predictor data and @var{Y}
## containing the class labels of observations in @var{X}.
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
## @code{@var{obj} = ClassificationKNN (@dots{}, @var{name}, @var{value})}
## returns a ClassificationKNN object with parameters specified by
## @qcode{Name-Value} pair arguments.  Type @code{help fitcknn} for more info.
##
## A @qcode{ClassificationKNN} object, @var{obj}, stores the labelled training
## data and various parameters for the k-Nearest Neighbor classification model,
## which can be accessed in the following fields:
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Field} @tab @tab @var{Description}
##
## @item @qcode{obj.X} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Y} @tab @tab Class labels, specified as a logical or
## numeric vector, or cell array of character vectors.  Each value in @var{Y} is
## the observed class label for the corresponding row in @var{X}.
##
## @item @qcode{obj.NumObservations} @tab @tab Number of observations used in
## training the ClassificationKNN model, specified as a positive integer scalar.
## This number can be less than the number of rows in the training data because
## rows containing @qcode{NaN} values are not part of the fit.
##
## @item @qcode{obj.RowsUsed} @tab @tab Rows of the original training data
## used in fitting the ClassificationKNN model, specified as a numerical vector.
## If you want to use this vector for indexing the training data in @var{X}, you
## have to convert it to a logical vector, i.e
## @qcode{X = obj.X(logical (obj.RowsUsed), :);}
##
## @item @qcode{obj.Standardize} @tab @tab A boolean flag indicating whether
## the data in @var{X} have been standardized prior to training.
##
## @item @qcode{obj.Sigma} @tab @tab Predictor standard deviations, specified
## as a numeric vector of the same length as the columns in @var{X}.  If the
## predictor variables have not been standardized, then @qcode{"obj.Sigma"} is
## empty.
##
## @item @qcode{obj.Mu} @tab @tab Predictor means, specified as a numeric
## vector of the same length as the columns in @var{X}.  If the predictor
## variables have not been standardized, then @qcode{"obj.Mu"} is empty.
##
## @item @qcode{obj.NumPredictors} @tab @tab The number of predictors
## (variables) in @var{X}.
##
## @item @qcode{obj.PredictorNames} @tab @tab Predictor variable names,
## specified as a cell array of character vectors.  The variable names are in
## the same order in which they appear in the training data @var{X}.
##
## @item @qcode{obj.ResponseName} @tab @tab Response variable name, specified
## as a character vector.
##
## @item @qcode{obj.ClassNames} @tab @tab Names of the classes in the training
## data @var{Y} with duplicates removed, specified as a cell array of character
## vectors.
##
## @item @qcode{obj.Prior} @tab @tab Prior probabilities for each class,
## specified as a numeric vector.  The order of the elements in @qcode{Prior}
## corresponds to the order of the classes in @qcode{ClassNames}.
##
## @item @qcode{obj.Cost} @tab @tab Cost of the misclassification of a point,
## specified as a square matrix. @qcode{Cost(i,j)} is the cost of classifying a
## point into class @qcode{j} if its true class is @qcode{i} (that is, the rows
## correspond to the true class and the columns correspond to the predicted
## class).  The order of the rows and columns in @qcode{Cost} corresponds to the
## order of the classes in @qcode{ClassNames}.  The number of rows and columns
## in @qcode{Cost} is the number of unique classes in the response.  By default,
## @qcode{Cost(i,j) = 1} if @qcode{i != j}, and @qcode{Cost(i,j) = 0} if
## @qcode{i = j}.  In other words, the cost is 0 for correct classification and
## 1 for incorrect classification.
##
## @item @qcode{obj.NumNeighbors} @tab @tab Number of nearest neighbors in
## @var{X} used to classify each point during prediction, specified as a
## positive integer value.
##
## @item @qcode{obj.Distance} @tab @tab Distance metric, specified as a
## character vector.  The allowable distance metric names depend on the choice
## of the neighbor-searcher method.  See the available distance metrics in
## @code{knnseaarch} for more info.
##
## @item @qcode{obj.DistanceWeight} @tab @tab Distance weighting function,
## specified as a function handle, which accepts a matrix of nonnegative
## distances, and returns a matrix the same size containing nonnegative distance
## weights.
##
## @item @qcode{obj.DistParameter} @tab @tab Parameter for the distance
## metric, specified either as a positive definite covariance matrix (when the
## distance metric is @qcode{"mahalanobis"}, or a positive scalar as the
## Minkowski distance exponent (when the distance metric is @qcode{"minkowski"},
## or a vector of positive scale values with length equal to the number of
## columns of @var{X} (when the distance metric is @qcode{"seuclidean"}.  For
## any other distance metric, the value of @qcode{DistParameter} is empty.
##
## @item @qcode{obj.NSMethod} @tab @tab Nearest neighbor search method,
## specified as either @qcode{"kdtree"}, which creates and uses a Kd-tree to
## find nearest neighbors, or @qcode{"exhaustive"}, which uses the exhaustive
## search algorithm by computing the distance values from all points in @var{X}
## to find nearest neighbors.
##
## @item @qcode{obj.IncludeTies} @tab @tab A boolean flag indicating whether
## prediction includes all the neighbors whose distance values are equal to the
## @math{k-th} smallest distance.  If @qcode{IncludeTies} is @qcode{true},
## prediction includes all of these neighbors.  Otherwise, prediction uses
## exactly @math{k} neighbors.
##
## @item @qcode{obj.BucketSize} @tab @tab Maximum number of data points in the
## leaf node of the Kd-tree, specified as positive integer value. This argument
## is meaningful only when @qcode{NSMethod} is @qcode{"kdtree"}.
##
## @end multitable
##
## @seealso{fitcknn, @@ClassificationKNN/predict}
## @end deftypefn

  properties (Access = public)

    X = [];                   # Predictor data
    Y = [];                   # Class labels

    NumObservations = [];     # Number of observations in training dataset
    RowsUsed        = [];     # Rows used in fitting
    Standardize     = [];     # Flag to standardize predictors
    Sigma           = [];     # Predictor standard deviations
    Mu              = [];     # Predictor means

    NumPredictors   = [];     # Number of predictors
    PredictorNames  = [];     # Predictor variables names
    ResponseName    = [];     # Response variable name
    ClassNames      = [];     # Names of classes in Y
    Prior           = [];     # Prior probability for each class
    Cost            = [];     # Cost of misclassification

    NumNeighbors    = [];     # Number of nearest neighbors
    Distance        = [];     # Distance metric
    DistanceWeight  = [];     # Distance weighting function
    DistParameter   = [];     # Parameter for distance metric
    NSMethod        = [];     # Nearest neighbor search method
    IncludeTies     = [];     # Flag for handling ties
    BucketSize      = [];     # Maximum data points in each node

  endproperties

  properties (Acesss = protected)


  endproperties

  methods (Access = public)
    ## Class object contructor
    function this = ClassificationKNN (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationKNN: too few input arguments.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error ("ClassificationKNN: number of rows in X and Y must be equal.");
      endif

      ## Assign original X and Y data to the ClassificationKNN object
      this.X = X;
      this.Y = Y;

      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);

      ## Set default values before parsing optional parameters
      Standardize     = false;  # Flag to standardize predictors
      PredictorNames  = [];     # Predictor variables names
      ResponseName    = [];     # Response variable name
      ClassNames      = [];     # Names of classes in Y
      Prior           = [];     # Prior probability for each class
      Cost            = [];     # Cost of misclassification
      Scale           = [];     # Distance scale for 'seuclidean'
      Cov             = [];     # Covariance matrix for 'mahalanobis'
      Exponent        = [];     # Exponent for 'minkowski'
      NumNeighbors    = [];     # Number of nearest neighbors
      Distance        = [];     # Distance metric
      DistanceWeight  = [];     # Distance weighting function
      DistParameter   = [];     # Parameter for distance metric
      NSMethod        = [];     # Nearest neighbor search method
      IncludeTies     = false;  # Flag for handling ties
      BucketSize      = 50;     # Maximum data points in each node

      ## Number of parameters for Standardize, Scale, Cov (maximum 1 allowed)
      SSC = 0;

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "standardize"
            if (SSC < 1)
              Standardize = varargin{2};
              if (! (Standardize == true || Standardize == false))
                error (strcat (["ClassificationKNN: Standardize must"], ...
                               [" be either true or false."]));
              endif
            else
              error (strcat (["ClassificationKNN: Standardize cannot"], ...
                             [" simultaneously be specified with either"], ...
                             [" Scale or Cov."]));
            endif

          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat (["ClassificationKNN: PredictorNames must"], ...
                             [" be supplied as a cellstring array."]));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat (["ClassificationKNN: PredictorNames must"], ...
                             [" have the same number of columns as X."]));
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error ("ClassificationKNN: ResponseName must be a char string.");
            endif

          case "classnames"
            ClassNames = varargin{2};
            if (! iscellstr (ClassNames))
              error (strcat (["ClassificationKNN: ClassNames must"], ...
                               [" be a cellstring array."]));
            endif
            ## Check that all class names are available in gnY
            if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                 ClassNames, "UniformOutput", false))))
              error ("ClassificationKNN: not all ClassNames are present in Y");
            endif

          case "prior"
            Prior = varargin{2};
            if (! ((isnumeric (Prior) && isvector (Prior)) ||
                  (strcmpi (Prior, "empirical") || strcmpi (Prior, "uniform"))))
              error (strcat (["ClassificationKNN: Prior must be either a"], ...
                             [" numeric vector or a string."]));
            endif

          case "cost"
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error ("ClassificationKNN: Cost must be a numeric square matrix.");
            endif

          case "numneighbors"
            NumNeighbors = varargin{2};
            if (! (isnumeric (NumNeighbors) && isscalar (NumNeighbors) &&
                                           fix (NumNeighbors) == NumNeighbors))
              error (strcat (["ClassificationKNN: NumNeighbors must be a"], ...
                             [" positive integer."]));
            endif

          case "distance"
            Distance = varargin{2};
            DMs = {"euclidean", "seuclidean", "mahalanobis", "minkowski", ...
                   "cityblock", "manhattan", "chebychev", "cosine", ...
                   "correlation", "spearman", "hamming", "jaccard"};
            if (! any (strcmpi (DMs, Distance)))
              error ("ClassificationKNN: unknown distance metric.");
            endif

          case "distanceweight"
            DistanceWeight = varargin{2};
            DMs = {"equal", "inverse", "squareinverse"};
            if (is_function_handle (DistanceWeight))
              m = eye (5);
              if (! isequal (size (m), size (DistanceWeight(m))))
                error (strcat (["ClassificationKNN: function handle for"], ...
                               [" distance weight must return the same"], ...
                               [" size as its input."]));
              endif
              this.DistanceWeight = DistanceWeight;
            else
              if (! any (strcmpi (DMs, DistanceWeight)))
                error ("ClassificationKNN: invalid distance weight.");
              endif
              if (strcmpi ("equal", DistanceWeight))
                this.DistanceWeight = @(x) x;
              endif
              if (strcmpi ("inverse", DistanceWeight))
                this.DistanceWeight = @(x) x.^(-1);
              endif
              if (strcmpi ("squareinverse", DistanceWeight))
                this.DistanceWeight = @(x) x.^(-2);
              endif
            endif

          case "scale"
            if (SSC < 1)
              Scale = varargin{2};
              if (! (isnumeric (Scale) && isvector (Scale)))
                error ("ClassificationKNN: Scale must be a numeric vector.");
              endif
            else
              error (strcat (["ClassificationKNN: Scale cannot"], ...
                             [" simultaneously be specified with either"], ...
                             [" Standardize or Cov."]));
            endif

          case "cov"
            if (SSC < 1)
              Cov = varargin{2};
              try
                chol (Cov);
              catch ME
                error (strcat (["ClassificationKNN: Cov must be a"], ...
                               [" symmetric positive definite matrix."]));
              end_try_catch
            else
              error (strcat (["ClassificationKNN: Cov cannot"], ...
                             [" simultaneously be specified with either"], ...
                             [" Standardize or Scale."]));
            endif

          case "exponent"
            Exponent = varargin{2};
            if (! (isnumeric (Exponent) && isscalar (Exponent) &&
                           Exponent > 0 && fix (Exponent) == Exponent))
              error ("ClassificationKNN: Exponent must be a positive integer.");
            endif

          case "nsmethod"
            NSMethod = varargin{2};
            NSM = {"kdtree", "exhaustive"};
            if (! any (strcmpi (NSM, NSMethod)))
              error (strcat (["ClassificationKNN: NSMethod must"], ...
                             [" be either kdtree or exhaustive."]));
            endif

          case "includeties"
            IncludeTies = varargin{2};
            if (! (IncludeTies == true || IncludeTies == false))
              error (strcat (["ClassificationKNN: IncludeTies must"], ...
                             [" be either true or false."]));
            endif

          case "bucketsize"
            BucketSize = varargin{2};
            if (! (isnumeric (BucketSize) && isscalar (BucketSize) &&
                           BucketSize > 0 && fix (BucketSize) == BucketSize))
              error (strcat (["ClassificationKNN: BucketSize must be a"], ...
                             [" positive integer."]));
            endif

          otherwise
            error (strcat (["ClassificationKNN: invalid parameter name"],...
                           [" in optional pair arguments."]));

        endswitch
        varargin (1:2) = [];
      endwhile

      ## Get number of variables in training data
      ndims_X = columns (X);

      ## Assign the number of predictors to the ClassificationKNN object
      this.NumPredictors = ndims_X;

      ## Generate default predictors and response variabe names (if necessary)
      if (isempty (PredictorNames))
        for i = 1:ndims_X
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = "Y";
      endif

      ## Assign predictors and response variable names
      this.PredictorNames = PredictorNames;
      this.ResponseName   = ResponseName;

      ## Handle class names
      if (isempty (ClassNames))
        ClassNames = gnY;
      else
        ru = logical (zeros (size (Y)));
        for i = 1:numel (ClassNames)
          ac = find (gnY, ClassNames{i});
          ru = ru | ac;
        endfor
        X = X(ru);
        Y = Y(ru);
        gY = gY(ru);
      endif

      ## Remove missing values from X and Y
      RowsUsed  = ! logical (sum (isnan ([X, gY]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

      ## Renew groups in Y
      [gY, gnY, glY] = grp2idx (Y);
      this.ClassNames = gnY;

      ## Check X contains valid data
      if (! (isnumeric (X) && isfinite (X)))
        error ("ClassificationKNN: invalid values in X.");
      endif

      ## Assign the number of observations and their correspoding indices
      ## on the original data, which will be used for training the model,
      ## to the ClassificationKNN object
      this.NumObservations = rows (X);
      this.RowsUsed = cast (RowsUsed, "double");

      ## Handle Standardize flag
      if (Standardize)
        this.Standardize = true;
        this.Sigma = std (X, [], 1);
        this.Mu = mean (X, 1);
      else
        this.Standardize = false;
        this.Sigma = [];
        this.Mu = [];
      endif

      ## Handle Prior and Cost
      if (isempty (Prior) || strcmpi ("uniform", Prior))
        this.Prior = ones (size (gnY)) ./ numel (gnY);
      elseif (strcmpi ("empirical", Prior))
        pr = [];
        for i = 1:numel (gnY)
          pr = [pr; sum(gY==i)];
        endfor
        this.Prior = pr ./ sum (pr);
      elseif (isnumeric (Prior))
        if (numel (gnY) != numel (Prior))
          error (strcat (["ClassificationKNN: the elements in Prior do not"], ...
                         [" correspond to selected classes in Y."]));
        endif
        this.Prior = Prior ./ sum (Prior);
      endif
      if (isempty (Cost))
        this.Cost = eye (numel (gnY));
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (strcat (["ClassificationKNN: the number of rows and"], ...
                         [" columns in Cost must correspond to selected"], ...
                         [" classes in Y."]));
        endif
        this.Cost = Cost;
      endif

      ## Get number of neighbors
      if (isempty (NumNeighbors))
        this.NumNeighbors = 1;
      else
        this.NumNeighbors = NumNeighbors;
      endif

      ## Get distance metric
      if (isempty (Distance))
        Distance = "euclidean";
      endif
      this.Distance = Distance;

      ## Get distance weight
      if (isempty (DistanceWeight))
        this.DistanceWeight = @(x) x;
      endif

      ## Handle distance metric parameters (Scale, Cov, Exponent)
      if (! isempty (Scale))
        if (! strcmpi (Distance, "seuclidean"))
          error (strcat (["ClassificationKNN: Scale is only valid when"], ...
                         [" distance metric is 'seuclidean'."]));
        endif
        if (numel (Scale) != ndims_X)
          error (strcat (["ClassificationKNN: Scale vector must have"], ...
                         [" equal length to the number of columns in X."]));
        endif
        if (any (Scale < 0))
          error (strcat (["ClassificationKNN: Scale vector must contain"], ...
                         [" nonnegative scalar values."]));
        endif
        this.DistParameter = Scale;
      else
        if (strcmpi (Distance, "seuclidean"))
          this.DistParameter = std (X, [], 1);
        endif
      endif
      if (! isempty (Cov))
        if (! strcmpi (Distance, "mahalanobis"))
          error (strcat (["ClassificationKNN: Cov is only valid when"], ...
                         [" distance metric is 'mahalanobis'."]));
        endif
        if (columns (Cov) != ndims_X)
          error (strcat (["ClassificationKNN: Cov matrix must have"], ...
                         [" equal columns as X."]));
        endif
        this.DistParameter = Cov;
      else
        if (strcmpi (Distance, "mahalanobis"))
          this.DistParameter = cov (X);
        endif
      endif
      if (! isempty (Exponent))
        if (! strcmpi (Distance, "minkowski"))
          error (strcat (["ClassificationKNN: Exponent is only valid when"], ...
                         [" distance metric is 'minkowski'."]));
        endif
        this.DistParameter = Exponent;
      else
        if (strcmpi (Distance, "minkowski"))
          this.DistParameter = 2;
        endif
      endif

      ## Get Nearest neighbor search method
      kdm = {"euclidean", "cityblock", "manhattan", "minkowski", "chebychev"};
      if (! isempty (NSMethod))
        if (strcmpi ("kdtree", NSMethod) && any (strcmpi (kdm, Distance)))
          error (strcat (["ClassificationKNN: kdtree method is only valid"], ...
                         [" for euclidean, cityblock, manhattan,"], ...
                         [" minkowski, and chebychev distance metrics."]));
        endif
        this.NSMethod = NSMethod;
      else
        if (any (strcmpi (kdm, Distance)) && ndims_X <= 10)
          this.NSMethod = "kdtree";
        else
          this.NSMethod = "exhaustive";
        endif
      endif

      ## Assign IncludeTies and BucketSize properties
      this.IncludeTies = IncludeTies;
      this.BucketSize = BucketSize;

    endfunction

  endmethods

endclassdef


%!demo
%! ## Create a ClassificationKNN object using fisheriris dataset
%! load fisheriris
%! x = meas;
%! y = species;
%! a = ClassificationKNN (x, y, "NumNeighbors", 5)

## Test constructor with NSMethod and NumNeighbors parameters
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = ClassificationKNN (x, y);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 1})
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = ClassificationKNN (x, y, "NSMethod", "exhaustive");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 1})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! k = 10;
%! a = ClassificationKNN (x, y, "NumNeighbors" ,k);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 10})
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = ones (4, 11);
%! y = ["a"; "a"; "b"; "b"];
%! k = 10;
%! a = ClassificationKNN (x, y, "NumNeighbors" ,k);
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 10})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! k = 10;
%! a = ClassificationKNN (x, y, "NumNeighbors" ,k, "NSMethod", "exhaustive");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 10})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! k = 10;
%! a = ClassificationKNN (x, y, "NumNeighbors" ,k, "Distance", "hamming");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.X, a.Y, a.NumNeighbors}, {x, y, 10})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "hamming"})
%! assert ({a.BucketSize}, {50})

## Test constructor with Standardize and DistParameter parameters
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! weights = ones (4,1);
%! a = ClassificationKNN (x, y, "standardize", 1);
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
%! a = ClassificationKNN (x, y, "standardize", false);
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
%! a = ClassificationKNN (x, y, "scale" , s, "Distance", "seuclidean");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.DistParameter}, {s})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "seuclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! C = cov (x);
%! a = ClassificationKNN (x, y, "cov" , C, "distance", "mahalanobis");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.DistParameter}, {C})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "mahalanobis"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = ClassificationKNN (x, y, "Exponent" , 5, "Distance", "minkowski");
%! assert (class (a), "ClassificationKNN");
%! assert (a.DistParameter, 5)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "minkowski"})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = ClassificationKNN (x, y, "Exponent" , 5, "Distance", "minkowski", ...
%!                        "NSMethod", "exhaustive");
%! assert (class (a), "ClassificationKNN");
%! assert (a.DistParameter, 5)
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "minkowski"})

## Test constructor with BucketSize and IncludeTies parameters
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = ClassificationKNN (x, y, "bucketsize" , 20, "distance", "mahalanobis");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "mahalanobis"})
%! assert ({a.BucketSize}, {20})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = ClassificationKNN (x, y, "includeties", true);
%! assert (class (a), "ClassificationKNN");
%! assert (a.IncludeTies, true);
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = ClassificationKNN (x, y);
%! assert (class (a), "ClassificationKNN");
%! assert (a.IncludeTies, false);
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})

## Test constructor with Prior and Cost parameters
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = ClassificationKNN (x, y);
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, [0.5; 0.5])
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! prior = [0.5; 0.5];
%! a = ClassificationKNN (x, y, 'prior', "empirical");
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, prior)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "a"; "b"];
%! prior = [0.75; 0.25];
%! a = ClassificationKNN (x, y, 'prior', "empirical");
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, prior)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "a"; "b"];
%! prior = [0.5; 0.5];
%! a = ClassificationKNN (x, y, 'prior', "uniform");
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, prior)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! cost = eye (2);
%! a = ClassificationKNN (x, y, 'cost', cost);
%! assert (class (a), "ClassificationKNN")
%! assert (a.Cost, [1, 0; 0, 1])
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! cost = eye (2);
%! a = ClassificationKNN (x, y, 'cost', cost, "Distance", "hamming" );
%! assert (class (a), "ClassificationKNN")
%! assert (a.Cost, [1, 0; 0, 1])
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "hamming"})
%! assert ({a.BucketSize}, {50})

## Test input validation
%!error<ClassificationKNN: too few input arguments.> ClassificationKNN ()
%!error<ClassificationKNN: too few input arguments.> ...
%! ClassificationKNN (ones(4, 1))
%!error<ClassificationKNN: number of rows in X and Y must be equal.> ...
%! ClassificationKNN(ones (4,2), ones (1,4))
%!error<ClassificationKNN: invalid values in X.> ...
%! ClassificationKNN([1;2;3;'a';4], ones (5,1))
%!error<ClassificationKNN: invalid values in X.> ...
%! ClassificationKNN([1;2;3;Inf;4], ones (5,1))
%!error<ClassificationKNN: Standardize must be either true or false.> ...
%! ClassificationKNN(ones (5,3), ones (5,1), "standardize", "a")
%!error<ClassificationKNN: Standardize cannot simultaneously be specified with> ...
%! ClassificationKNN(ones (5,2), ones (5,1), "standardize", "a", "scale", [1 1])
%!error<ClassificationKNN: ClassificationKNN: PredictorNames must be supplied as a cellstring array.> ...
%! ClassificationKNN(ones (5,2), ones (5,1), "PredictorNames", ["A"])
%!error<ClassificationKNN: ClassificationKNN: PredictorNames must be supplied as a cellstring array.> ...
%! ClassificationKNN(ones (5,2), ones (5,1), "PredictorNames", "A")
%!error<ClassificationKNN: ClassificationKNN: PredictorNames must have same number of columns as X.> ...
%! ClassificationKNN(ones (5,2), ones (5,1), "PredictorNames", {"A", "B", "C"})
%!error<ClassificationKNN: ClassificationKNN: ResponseName must be a char string.> ...
%! ClassificationKNN(ones (5,2), ones (5,1), "ResponseName", {"Y"})
%!error<ClassificationKNN: ClassificationKNN: ResponseName must be a char string.> ...
%! ClassificationKNN(ones (5,2), ones (5,1), "ResponseName", 1)


%!error<ClassificationKNN: invalid NAME in optional pairs of arguments.> ...
%! ClassificationKNN(ones (4,2), ones (4,1), "some","some")
%!error<ClassificationKNN: invalid value of k.> ...
%! ClassificationKNN(ones (4,2),ones (4,1), "K", NaN)
%!error<ClassificationKNN: invalid value of k.> ...
%! ClassificationKNN(ones (4,2),ones (4,1), "K", -5)
%!error<ClassificationKNN: invalid value of k.> ...
%! ClassificationKNN(ones (4,2),ones (4,1), "K", 3.14)
%!error<ClassificationKNN: Weights has invalid observarions.> ...
%! ClassificationKNN(ones(4,2),ones(4,1), "weights", ones(2,2))
%!error<ClassificationKNN: Weights has invalid observarions.> ...
%! ClassificationKNN(ones(4,2),ones(4,1), "weights", [1; 2; 3; "a"; 4])
%!error<ClassificationKNN: Invalid value of Minkowski Exponent.> ...
%! ClassificationKNN(ones(4,2),ones(4,1), "P", -5)
%!error<ClassificationKNN: Invalid value in cost or the size of cost.> ...
%! ClassificationKNN(ones(4,2),ones(4,1), "cost", ones (2,2))
%!error<ClassificationKNN: Invalid value in cost or the size of cost.> ...
%! ClassificationKNN(ones(4,2),ones(4,1), "cost", [1; 2; 3; "a"; 4])
%!error<ClassificationKNN: Invalid value in Scale or the size of scale.> ...
%! ClassificationKNN(ones(4,2),ones(4,1), "Scale", [1; 2; 3; "a"; 4])
%!error<ClassificationKNN: Invalid value in Scale or the size of scale.> ...
%! ClassificationKNN(ones(4,2),ones(4,1), "Scale", zeros (3,3))
%!error<ClassificationKNN: Invalid value in COV, COV can only be given for mahalanobis distance.> ...
%! ClassificationKNN(ones(4,2),ones(4,1), "cov", ones (4, 1))
%!error<ClassificationKNN: Invalid value of bucketsize.> ...
%! ClassificationKNN(ones(4,2),ones(4,1), "bucketsize", -5)
%!error<ClassificationKNN: value of standardize must be boolean.> ...
%! ClassificationKNN(ones(4,2),ones(4,1), "standardize", 5)
