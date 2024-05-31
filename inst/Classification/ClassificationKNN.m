## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2023-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @item @qcode{obj.BreakTies} @tab @tab Tie-breaking algorithm used by predict
## when multiple classes have the same smallest cost, specified as one of the
## following character arrays: @qcode{"smallest"} (default), which favors the
## class with the smallest index among the tied groups, i.e. the one that
## appears first in the training labelled data.  @qcode{"nearest"}, which favors
## the class with the nearest neighbor among the tied groups, i.e. the class
## with the closest member point according to the distance metric used.
## @qcode{"random"}, which randomly picks one class among the tied groups.
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
## @math{k^th} smallest distance.  If @qcode{IncludeTies} is @qcode{true},
## prediction includes all of these neighbors.  Otherwise, prediction uses
## exactly @math{k} neighbors.
##
## @item @qcode{obj.BucketSize} @tab @tab Maximum number of data points in the
## leaf node of the Kd-tree, specified as positive integer value. This argument
## is meaningful only when @qcode{NSMethod} is @qcode{"kdtree"}.
##
## @end multitable
##
## @seealso{fitcknn, knnsearch, rangesearch, pdist2}
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
    BreakTies       = [];     # Tie-breaking algorithm
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

  methods (Access = public)

    ## Class object constructor
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
      BreakTies       = [];     # Tie-breaking algorithm
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
              SSC += 1;
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
              error ("ClassificationKNN: ResponseName must be a character array.");
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
              error ("ClassificationKNN: not all ClassNames are present in Y.");
            endif

          case "breakties"
            BreakTies = varargin{2};
            if (! ischar (BreakTies))
              error ("ClassificationKNN: BreakTies must be a character array.");
            endif
            ## Check that all class names are available in gnY
            BTs = {"smallest", "nearest", "random"};
            if (! any (strcmpi (BTs, BreakTies)))
              error ("ClassificationKNN: invalid value for BreakTies.");
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
                   NumNeighbors > 0 && fix (NumNeighbors) == NumNeighbors))
              error (strcat (["ClassificationKNN: NumNeighbors must be a"], ...
                             [" positive integer."]));
            endif

          case "distance"
            Distance = varargin{2};
            DMs = {"euclidean", "seuclidean", "mahalanobis", "minkowski", ...
                   "cityblock", "manhattan", "chebychev", "cosine", ...
                   "correlation", "spearman", "hamming", "jaccard"};
            if (ischar (Distance))
              if (! any (strcmpi (DMs, Distance)))
                error ("ClassificationKNN: unsupported distance metric.");
              endif
            elseif (is_function_handle (Distance))
              ## Check the input output sizes of the user function
              D2 = [];
              try
                D2 = Distance (X(1,:), Y);
              catch ME
                error (strcat (["ClassificationKNN: invalid function"], ...
                               [" handle for distance metric."]));
              end_try_catch
              Yrows = rows (Y);
              if (! isequal (size (D2), [Yrows, 1]))
                error (strcat (["ClassificationKNN: custom distance"], ...
                               [" function produces wrong output size."]));
              endif
            else
              error ("ClassificationKNN: invalid distance metric.");
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
              SSC += 1;
            else
              error (strcat (["ClassificationKNN: Scale cannot"], ...
                             [" simultaneously be specified with either"], ...
                             [" Standardize or Cov."]));
            endif

          case "cov"
            if (SSC < 1)
              Cov = varargin{2};
              [~, p] = chol (Cov);
              if (p != 0)
                error (strcat (["ClassificationKNN: Cov must be a"], ...
                               [" symmetric positive definite matrix."]));
              endif
              SSC += 1;
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
            if (! ischar (NSMethod))
              error ("ClassificationKNN: NSMethod must be a character array.");
            endif
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
        this.Sigma(this.Sigma == 0) = 1;  # predictor is constant
        this.Mu = mean (X, 1);
      else
        this.Standardize = false;
        this.Sigma = [];
        this.Mu = [];
      endif

      ## Handle BreakTies
      if (isempty (BreakTies))
        this.BreakTies = "smallest";
      else
        this.BreakTies = BreakTies;
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
        this.Cost = cast (! eye (numel (gnY)), "double");
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
                         [" distance metric is seuclidean."]));
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
          if (Standardize)
            this.DistParameter = ones (1, ndims_X);
          else
            this.DistParameter = std (X, [], 1);
          endif
        endif
      endif
      if (! isempty (Cov))
        if (! strcmpi (Distance, "mahalanobis"))
          error (strcat (["ClassificationKNN: Cov is only valid when"], ...
                         [" distance metric is mahalanobis."]));
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
                         [" distance metric is minkowski."]));
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
        if (strcmpi ("kdtree", NSMethod) && (! any (strcmpi (kdm, Distance))))
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

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKNN} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationKNN} {[@var{label}, @var{score}, @var{cost}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the kNN algorithm from a
    ## k-Nearest Neighbor classification model.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the matrix of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## predictor data in @code{obj.X} and corresponding labels, @code{obj.Y},
    ## stored in the k-Nearest Neighbor classification model, @var{obj}.
    ##
    ## @var{XC} must be an @math{MxP} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the kNN model in
    ## @var{obj}.
    ##
    ## @code{[@var{label}, @var{score}, @var{cost}] = predict (@var{obj}, @var{XC})}
    ## also returns @var{score}, which contains the predicted class scores or
    ## posterior probabilities for each instance of the corresponding unique
    ## classes, and @var{cost}, which is a matrix containing the expected cost
    ## of the classifications.
    ##
    ## @seealso{fitcknn, ClassificationKNN}
    ## @end deftypefn
    function [label, score, cost] = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("ClassificationKNN.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("ClassificationKNN.predict: XC is empty.");
      elseif (columns (this.X) != columns (XC))
        error (strcat (["ClassificationKNN.predict: XC must have the same"], ...
                       [" number of features (columns) as in the kNN model."]));
      endif

      ## Get training data and labels
      X = this.X(logical (this.RowsUsed),:);
      Y = this.Y(logical (this.RowsUsed),:);

      ## Standardize (if necessary)
      if (this.Standardize)
        X = (X - this.Mu) ./ this.Sigma;
        XC = (XC - this.Mu) ./ this.Sigma;
      endif

      ## Train kNN
      if (strcmpi (this.Distance, "seuclidean"))
        [idx, dist] = knnsearch (X, XC, "k", this.NumNeighbors, ...
                      "NSMethod", this.NSMethod, "Distance", "seuclidean", ...
                      "Scale", this.DistParameter, "sortindices", true, ...
                      "includeties", this.IncludeTies, ...
                      "bucketsize", this.BucketSize);

      elseif (strcmpi (this.Distance, "mahalanobis"))
        [idx, dist] = knnsearch (X, XC, "k", this.NumNeighbors, ...
                      "NSMethod", this.NSMethod, "Distance", "mahalanobis", ...
                      "cov", this.DistParameter, "sortindices", true, ...
                      "includeties", this.IncludeTies, ...
                      "bucketsize", this.BucketSize);

      elseif (strcmpi (this.Distance, "minkowski"))
        [idx, dist] = knnsearch (X, XC, "k", this.NumNeighbors, ...
                      "NSMethod", this.NSMethod, "Distance", "minkowski", ...
                      "P", this.DistParameter, "sortindices", true, ...
                      "includeties",this.IncludeTies, ...
                      "bucketsize", this.BucketSize);

      else
        [idx, dist] = knnsearch (X, XC, "k", this.NumNeighbors, ...
                      "NSMethod", this.NSMethod, "Distance", this.Distance, ...
                      "sortindices", true, "includeties", this.IncludeTies, ...
                      "bucketsize", this.BucketSize);
      endif

      ## Make prediction
      label = {};
      score = [];
      cost  = [];

      ## Get IDs and labels for each point in training data
      [gY, gnY, glY] = grp2idx (Y);

      ## Evaluate the K nearest neighbours for each new point
      for i = 1:rows (idx)

        ## Get K nearest neighbours
        if (this.IncludeTies)
          NN_idx = idx{i};
          NNdist = dist{i};
        else
          NN_idx = idx(i,:);
          NNdist = dist(i,:);
        endif
        k = numel (NN_idx);
        kNNgY = gY(NN_idx);

        ## Count frequency for each class
        for c = 1:numel (this.ClassNames)
          freq(c) = sum (kNNgY == c) / k;
        endfor

        ## Get label according to BreakTies
        if (strcmpi (this.BreakTies, "smallest"))
          [~, idl] = max (freq);
        else
          idl = find (freq == max (freq));
          tgn = numel (idl);
          if (tgn > 1)
            if (strcmpi (this.BreakTies, "nearest"))
              for t = 1:tgn
                tgs(t) = find (gY(NN_idx) == idl(t));
              endfor
              [~, idm] = min (tgs);
              idl = idl(idm);
            else      # "random"
              idl = idl(randperm (numel (idl))(1));
            endif
          endif
        endif
        label = [label; gnY{idl}];

        ## Calculate score and cost
        score = [score; freq];
        cost = [cost; 1-freq];

      endfor

    endfunction

  endmethods

endclassdef


%!demo
%! ## Create a k-nearest neighbor classifier for Fisher's iris data with k = 5.
%! ## Evaluate some model predictions on new data.
%!
%! load fisheriris
%! x = meas;
%! y = species;
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 5, "Standardize", 1);
%! [label, score, cost] = predict (obj, xc)

%!demo
%! ## Train a k-nearest neighbor classifier for k = 10
%! ## and plot the decision boundaries.
%!
%! load fisheriris
%! idx = ! strcmp (species, "setosa");
%! X = meas(idx,3:4);
%! Y = cast (strcmpi (species(idx), "virginica"), "double");
%! obj = fitcknn (X, Y, "Standardize", 1, "NumNeighbors", 10, "NSMethod", "exhaustive")
%! x1 = [min(X(:,1)):0.03:max(X(:,1))];
%! x2 = [min(X(:,2)):0.02:max(X(:,2))];
%! [x1G, x2G] = meshgrid (x1, x2);
%! XGrid = [x1G(:), x2G(:)];
%! pred = predict (obj, XGrid);
%! gidx = logical (str2num (cell2mat (pred)));
%!
%! figure
%! scatter (XGrid(gidx,1), XGrid(gidx,2), "markerfacecolor", "magenta");
%! hold on
%! scatter (XGrid(!gidx,1), XGrid(!gidx,2), "markerfacecolor", "red");
%! plot (X(Y == 0, 1), X(Y == 0, 2), "ko", X(Y == 1, 1), X(Y == 1, 2), "kx");
%! xlabel ("Petal length (cm)");
%! ylabel ("Petal width (cm)");
%! title ("5-Nearest Neighbor Classifier Decision Boundary");
%! legend ({"Versicolor Region", "Virginica Region", ...
%!         "Sampled Versicolor", "Sampled Virginica"}, ...
%!         "location", "northwest")
%! axis tight
%! hold off


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
%! a = ClassificationKNN (x, y, "Standardize", 1);
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
%! a = ClassificationKNN (x, y, "Standardize", false);
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
%! a = ClassificationKNN (x, y, "Scale" , s, "Distance", "seuclidean");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.DistParameter}, {s})
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "seuclidean"})
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
%! a = ClassificationKNN (x, y, "BucketSize" , 20, "distance", "mahalanobis");
%! assert (class (a), "ClassificationKNN");
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "mahalanobis"})
%! assert ({a.BucketSize}, {20})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = ClassificationKNN (x, y, "IncludeTies", true);
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
%! a = ClassificationKNN (x, y, "Prior", "empirical");
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, prior)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "a"; "b"];
%! prior = [0.75; 0.25];
%! a = ClassificationKNN (x, y, "Prior", "empirical");
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, prior)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "a"; "b"];
%! prior = [0.5; 0.5];
%! a = ClassificationKNN (x, y, "Prior", "uniform");
%! assert (class (a), "ClassificationKNN")
%! assert (a.Prior, prior)
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! cost = eye (2);
%! a = ClassificationKNN (x, y, "Cost", cost);
%! assert (class (a), "ClassificationKNN")
%! assert (a.Cost, [1, 0; 0, 1])
%! assert ({a.NSMethod, a.Distance}, {"kdtree", "euclidean"})
%! assert ({a.BucketSize}, {50})
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! cost = eye (2);
%! a = ClassificationKNN (x, y, "Cost", cost, "Distance", "hamming" );
%! assert (class (a), "ClassificationKNN")
%! assert (a.Cost, [1, 0; 0, 1])
%! assert ({a.NSMethod, a.Distance}, {"exhaustive", "hamming"})
%! assert ({a.BucketSize}, {50})

## Test input validation for constructor
%!error<ClassificationKNN: too few input arguments.> ClassificationKNN ()
%!error<ClassificationKNN: too few input arguments.> ...
%! ClassificationKNN (ones(4, 1))
%!error<ClassificationKNN: number of rows in X and Y must be equal.> ...
%! ClassificationKNN (ones (4,2), ones (1,4))
%!error<ClassificationKNN: Standardize must be either true or false.> ...
%! ClassificationKNN (ones (5,3), ones (5,1), "standardize", "a")
%!error<ClassificationKNN: Standardize cannot simultaneously be specified with> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "scale", [1 1], "standardize", true)
%!error<ClassificationKNN: PredictorNames must be supplied as a cellstring array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "PredictorNames", ["A"])
%!error<ClassificationKNN: PredictorNames must be supplied as a cellstring array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "PredictorNames", "A")
%!error<ClassificationKNN: PredictorNames must have the same number of columns as X.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "PredictorNames", {"A", "B", "C"})
%!error<ClassificationKNN: ResponseName must be a character array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "ResponseName", {"Y"})
%!error<ClassificationKNN: ResponseName must be a character array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "ResponseName", 1)
%!error<ClassificationKNN: ClassNames must be a cellstring array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "ClassNames", 1)
%!error<ClassificationKNN: ClassNames must be a cellstring array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "ClassNames", ["1"])
%!error<ClassificationKNN: not all ClassNames are present in Y.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "ClassNames", {"1", "2"})
%!error<ClassificationKNN: BreakTies must be a character array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BreakTies", 1)
%!error<ClassificationKNN: BreakTies must be a character array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BreakTies", {"1"})
%!error<ClassificationKNN: invalid value for BreakTies.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BreakTies", "some")
%!error<ClassificationKNN: Prior must be either a numeric vector or a string.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Prior", {"1", "2"})
%!error<ClassificationKNN: Cost must be a numeric square matrix.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cost", [1, 2])
%!error<ClassificationKNN: Cost must be a numeric square matrix.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cost", "string")
%!error<ClassificationKNN: Cost must be a numeric square matrix.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cost", {eye(2)})
%!error<ClassificationKNN: NumNeighbors must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NumNeighbors", 0)
%!error<ClassificationKNN: NumNeighbors must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NumNeighbors", 15.2)
%!error<ClassificationKNN: NumNeighbors must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NumNeighbors", "asd")
%!error<ClassificationKNN: unsupported distance metric.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Distance", "somemetric")
%!error<ClassificationKNN: invalid function handle for distance metric.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Distance", ...
%!                    @(v,m)sqrt(repmat(v,rows(m),1)-m,2))
%!error<ClassificationKNN: custom distance function produces wrong output size.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Distance", ...
%!                    @(v,m)sqrt(sum(sumsq(repmat(v,rows(m),1)-m,2))))
%!error<ClassificationKNN: invalid distance metric.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Distance", [1 2 3])
%!error<ClassificationKNN: invalid distance metric.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Distance", {"mahalanobis"})
%!error<ClassificationKNN: invalid distance metric.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Distance", logical (5))
%!error<ClassificationKNN: function handle for distance weight must return the> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "DistanceWeight", @(x)sum(x))
%!error<ClassificationKNN: invalid distance weight.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "DistanceWeight", "text")
%!error<ClassificationKNN: invalid distance weight.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "DistanceWeight", [1 2 3])
%!error<ClassificationKNN: Scale must be a numeric vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Scale", "scale")
%!error<ClassificationKNN: Scale must be a numeric vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Scale", {[1 2 3]})
%!error<ClassificationKNN: Scale cannot simultaneously be specified with> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "standardize", true, "scale", [1 1])
%!error<ClassificationKNN: Cov must be a symmetric positive definite matrix.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cov", ones (2), "Distance", "mahalanobis")
%!error<ClassificationKNN: Cov cannot simultaneously be specified with> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "scale", [1 1], "Cov", ones (2))
%!error<ClassificationKNN: Exponent must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Exponent", 12.5)
%!error<ClassificationKNN: Exponent must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Exponent", -3)
%!error<ClassificationKNN: Exponent must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Exponent", "three")
%!error<ClassificationKNN: Exponent must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Exponent", {3})
%!error<ClassificationKNN: NSMethod must be a character array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NSMethod", {"kdtree"})
%!error<ClassificationKNN: NSMethod must be a character array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NSMethod", 3)
%!error<ClassificationKNN: NSMethod must be either kdtree or exhaustive.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NSMethod", "some")
%!error<ClassificationKNN: IncludeTies must be either true or false.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "IncludeTies", "some")
%!error<ClassificationKNN: BucketSize must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BucketSize", 42.5)
%!error<ClassificationKNN: BucketSize must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BucketSize", -50)
%!error<ClassificationKNN: BucketSize must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BucketSize", "some")
%!error<ClassificationKNN: BucketSize must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BucketSize", {50})
%!error<ClassificationKNN: invalid parameter name in optional pair arguments.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "some", "some")
%!error<ClassificationKNN: invalid values in X.> ...
%! ClassificationKNN ([1;2;3;'a';4], ones (5,1))
%!error<ClassificationKNN: invalid values in X.> ...
%! ClassificationKNN ([1;2;3;Inf;4], ones (5,1))
%!error<ClassificationKNN: the elements in Prior do not correspond to selected> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Prior", [1 2])
%!error<ClassificationKNN: the number of rows and columns in Cost must> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cost", [1 2; 1 3])
%!error<ClassificationKNN: Scale is only valid when distance metric is seuclidean> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Scale", [1 1])
%!error<ClassificationKNN: Scale vector must have equal length to the number of> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Scale", [1 1 1], "Distance", "seuclidean")
%!error<ClassificationKNN: Scale vector must contain nonnegative scalar values.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Scale", [1 -1], "Distance", "seuclidean")
%!error<ClassificationKNN: Cov is only valid when distance metric is mahalanobis> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cov", eye (2))
%!error<ClassificationKNN: Cov matrix must have equal columns as X.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cov", eye (3), "Distance", "mahalanobis")
%!error<ClassificationKNN: Exponent is only valid when distance metric is minkowski.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Exponent", 3)
%!error<ClassificationKNN: kdtree method is only valid for euclidean, cityblock> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Distance", "hamming", "NSMethod", "kdtree")

## Test output for predict method
%!shared x, y
%! load fisheriris
%! x = meas;
%! y = species;
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "versicolor"; "virginica"})
%! assert (s, [1, 0, 0; 0, 1, 0; 0, 0, 1])
%! assert (c, [0, 1, 1; 1, 0, 1; 1, 1, 0])
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 5, "Standardize", 1);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"; "versicolor"; "virginica"})
%! assert (s, [0.4, 0.6, 0; 0, 1, 0; 0, 0, 1])
%! assert (c, [0.6, 0.4, 1; 1, 0, 1; 1, 1, 0])
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "mahalanobis");
%! [l, s, c] = predict (obj, xc);
%! assert (s, [0.3, 0.7, 0; 0, 0.9, 0.1; 0.2, 0.2, 0.6], 1e-4)
%! assert (c, [0.7, 0.3, 1; 1, 0.1, 0.9; 0.8, 0.8, 0.4], 1e-4)
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "cosine");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "versicolor"; "virginica"})
%! assert (s, [1, 0, 0; 0, 1, 0; 0, 0.3, 0.7], 1e-4)
%! assert (c, [0, 1, 1; 1, 0, 1; 1, 0.7, 0.3], 1e-4)
%!test
%! xc = [5.2, 4.1, 1.5,	0.1; 5.1,	3.8, 1.9,	0.4; ...
%!         5.1, 3.8, 1.5, 0.3; 4.9, 3.6, 1.4, 0.1];
%! obj = fitcknn (x, y, "NumNeighbors", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "setosa"; "setosa"; "setosa"})
%! assert (s, [1, 0, 0; 1, 0, 0; 1, 0, 0; 1, 0, 0])
%! assert (c, [0, 1, 1; 0, 1, 1; 0, 1, 1; 0, 1, 1])
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.6, 0.4], 1e-4)
%! assert (c, [1, 0.4, 0.6], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "minkowski", "Exponent", 5);
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.5, 0.5], 1e-4)
%! assert (c, [1, 0.5, 0.5], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "jaccard");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.9, 0.1, 0], 1e-4)
%! assert (c, [0.1, 0.9, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "mahalanobis");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0.1000, 0.5000, 0.4000], 1e-4)
%! assert (c, [0.9000, 0.5000, 0.6000], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 5, "distance", "jaccard");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.8, 0.2, 0], 1e-4)
%! assert (c, [0.2, 0.8, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 5, "distance", "seuclidean");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 1, 0], 1e-4)
%! assert (c, [1, 0, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "chebychev");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.7, 0.3], 1e-4)
%! assert (c, [1, 0.3, 0.7], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "cityblock");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 0.6, 0.4], 1e-4)
%! assert (c, [1, 0.4, 0.6], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "cosine");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"virginica"})
%! assert (s, [0, 0.1, 0.9], 1e-4)
%! assert (c, [1, 0.9, 0.1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "correlation");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"virginica"})
%! assert (s, [0, 0.1, 0.9], 1e-4)
%! assert (c, [1, 0.9, 0.1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 30, "distance", "spearman");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"versicolor"})
%! assert (s, [0, 1, 0], 1e-4)
%! assert (c, [1, 0, 1], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 30, "distance", "hamming");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.4333, 0.3333, 0.2333], 1e-4)
%! assert (c, [0.5667, 0.6667, 0.7667], 1e-4)
%!test
%! xc = [5, 3, 5, 1.45];
%! obj = fitcknn (x, y, "NumNeighbors", 5, "distance", "hamming");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"})
%! assert (s, [0.8, 0.2, 0], 1e-4)
%! assert (c, [0.2, 0.8, 1], 1e-4)
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "correlation");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa"; "versicolor"; "virginica"})
%! assert (s, [1, 0, 0; 0, 1, 0; 0, 0.4, 0.6], 1e-4)
%! assert (c, [0, 1, 1; 1, 0, 1; 1, 0.6, 0.4], 1e-4)
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcknn (x, y, "NumNeighbors", 10, "distance", "hamming");
%! [l, s, c] = predict (obj, xc);
%! assert (l, {"setosa";"setosa";"setosa"})
%! assert (s, [0.9, 0.1, 0; 1, 0, 0; 0.5, 0, 0.5], 1e-4)
%! assert (c, [0.1, 0.9, 1; 0, 1, 1; 0.5, 1, 0.5], 1e-4)

## Test input validation for predict method
%!error<ClassificationKNN.predict: too few input arguments.> ...
%! predict (ClassificationKNN (ones (4,2), ones (4,1)))
%!error<ClassificationKNN.predict: XC is empty.> ...
%! predict (ClassificationKNN (ones (4,2), ones (4,1)), [])
%!error<ClassificationKNN.predict: XC must have the same number of features> ...
%! predict (ClassificationKNN (ones (4,2), ones (4,1)), 1)
