## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
## Copyright (C) 2023-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @multitable @columnfractions 0.23 0.02 0.75
## @headitem @var{Field} @tab @tab @var{Description}
##
## @item @qcode{X} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{Y} @tab @tab Class labels, specified as a logical or
## numeric vector, or cell array of character vectors.  Each value in @var{Y} is
## the observed class label for the corresponding row in @var{X}.
##
## @item @qcode{NumObservations} @tab @tab Number of observations used in
## training the ClassificationKNN model, specified as a positive integer scalar.
## This number can be less than the number of rows in the training data because
## rows containing @qcode{NaN} values are not part of the fit.
##
## @item @qcode{RowsUsed} @tab @tab Rows of the original training data
## used in fitting the ClassificationKNN model, specified as a numerical vector.
## If you want to use this vector for indexing the training data in @var{X}, you
## have to convert it to a logical vector, i.e
## @qcode{X = obj.X(logical (obj.RowsUsed), :);}
##
## @item @qcode{Standardize} @tab @tab A boolean flag indicating whether
## the data in @var{X} have been standardized prior to training.
##
## @item @qcode{Sigma} @tab @tab Predictor standard deviations, specified
## as a numeric vector of the same length as the columns in @var{X}.  If the
## predictor variables have not been standardized, then @qcode{"obj.Sigma"} is
## empty.
##
## @item @qcode{Mu} @tab @tab Predictor means, specified as a numeric
## vector of the same length as the columns in @var{X}.  If the predictor
## variables have not been standardized, then @qcode{"obj.Mu"} is empty.
##
## @item @qcode{NumPredictors} @tab @tab The number of predictors
## (variables) in @var{X}.
##
## @item @qcode{PredictorNames} @tab @tab Predictor variable names,
## specified as a cell array of character vectors.  The variable names are in
## the same order in which they appear in the training data @var{X}.
##
## @item @qcode{ResponseName} @tab @tab Response variable name, specified
## as a character vector.
##
## @item @qcode{ClassNames} @tab @tab Names of the classes in the training
## data @var{Y} with duplicates removed, specified as a cell array of character
## vectors.
##
## @item @qcode{BreakTies} @tab @tab Tie-breaking algorithm used by predict
## when multiple classes have the same smallest cost, specified as one of the
## following character arrays: @qcode{"smallest"} (default), which favors the
## class with the smallest index among the tied groups, i.e. the one that
## appears first in the training labelled data.  @qcode{"nearest"}, which favors
## the class with the nearest neighbor among the tied groups, i.e. the class
## with the closest member point according to the distance metric used.
## @qcode{"random"}, which randomly picks one class among the tied groups.
##
## @item @qcode{Prior} @tab @tab Prior probabilities for each class,
## specified as a numeric vector.  The order of the elements in @qcode{Prior}
## corresponds to the order of the classes in @qcode{ClassNames}.
##
## @item @qcode{Cost} @tab @tab Cost of the misclassification of a point,
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
## @item @qcode{ScoreTransform} @tab @tab A function_handle which is used
## for transforming the kNN prediction score into a posterior probability.  By
## default, it is @qcode{'none'}, in which case the @code{predict} and
## @code{resubPredict} methods return the prediction scores.
##
## @item @qcode{NumNeighbors} @tab @tab Number of nearest neighbors in
## @var{X} used to classify each point during prediction, specified as a
## positive integer value.
##
## @item @qcode{Distance} @tab @tab Distance metric, specified as a
## character vector.  The allowable distance metric names depend on the choice
## of the neighbor-searcher method.  See the available distance metrics in
## @code{knnseaarch} for more info.
##
## @item @qcode{DistanceWeight} @tab @tab Distance weighting function,
## specified as a function handle, which accepts a matrix of nonnegative
## distances, and returns a matrix the same size containing nonnegative distance
## weights.
##
## @item @qcode{DistParameter} @tab @tab Parameter for the distance
## metric, specified either as a positive definite covariance matrix (when the
## distance metric is @qcode{"mahalanobis"}, or a positive scalar as the
## Minkowski distance exponent (when the distance metric is @qcode{"minkowski"},
## or a vector of positive scale values with length equal to the number of
## columns of @var{X} (when the distance metric is @qcode{"seuclidean"}.  For
## any other distance metric, the value of @qcode{DistParameter} is empty.
##
## @item @qcode{NSMethod} @tab @tab Nearest neighbor search method,
## specified as either @qcode{"kdtree"}, which creates and uses a Kd-tree to
## find nearest neighbors, or @qcode{"exhaustive"}, which uses the exhaustive
## search algorithm by computing the distance values from all points in @var{X}
## to find nearest neighbors.
##
## @item @qcode{IncludeTies} @tab @tab A boolean flag indicating whether
## prediction includes all the neighbors whose distance values are equal to the
## @math{k^th} smallest distance.  If @qcode{IncludeTies} is @qcode{true},
## prediction includes all of these neighbors.  Otherwise, prediction uses
## exactly @math{k} neighbors.
##
## @item @qcode{BucketSize} @tab @tab Maximum number of data points in the
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
    NumPredictors   = [];     # Number of predictors
    PredictorNames  = [];     # Predictor variables names
    ResponseName    = [];     # Response variable name
    ClassNames      = [];     # Names of classes in Y
    Prior           = [];     # Prior probability for each class
    Cost            = [];     # Cost of misclassification

    ScoreTransform  = [];     # Transformation for classification scores

    Standardize     = [];     # Flag to standardize predictors
    Sigma           = [];     # Predictor standard deviations
    Mu              = [];     # Predictor means

    BreakTies       = [];     # Tie-breaking algorithm
    NumNeighbors    = [];     # Number of nearest neighbors
    Distance        = [];     # Distance metric
    DistanceWeight  = [];     # Distance weighting function
    DistParameter   = [];     # Parameter for distance metric
    NSMethod        = [];     # Nearest neighbor search method
    IncludeTies     = [];     # Flag for handling ties
    BucketSize      = [];     # Maximum data points in each node

  endproperties

  methods (Hidden)

    ## Custom display
    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    ## Custom display
    function disp (this)
      fprintf ("\n  ClassificationKNN\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      if (iscellstr (this.ClassNames))
        str = repmat ({"'%s'"}, 1, numel (this.ClassNames));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.ClassNames{:});
      elseif (ischar (this.ClassNames))
        str = repmat ({"'%s'"}, 1, rows (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, cellstr (this.ClassNames){:});
      else # single, double, logical
        str = repmat ({"%d"}, 1, numel (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, this.ClassNames);
      endif
      fprintf ("%+25s: %s\n", 'ClassNames', str);
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumPredictors', this.NumPredictors);
      fprintf ("%+25s: '%s'\n", 'Distance', this.Distance);
      fprintf ("%+25s: '%s'\n", 'NSMethod', this.NSMethod);
      fprintf ("%+25s: %d\n", 'NumNeighbors', this.NumNeighbors);
    endfunction

    ## Class specific subscripted reference
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("Invalid () indexing for referencing values", ...
                         " in a ClassificationKNN object."));
        case '{}'
          error (strcat ("Invalid {} indexing for referencing values", ...
                         " in a ClassificationKNN object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("ClassificationKNN.subsref: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("ClassificationKNN.subref:", ...
                           " unrecongized property: '%s'"), s.subs);
          end_try_catch
      endswitch
      ## Chained references
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

    ## Class specific subscripted assignment
    function this = subsasgn (this, s, val)
      if (numel (s) > 1)
        error (strcat ("ClassificationKNN.subsasgn:", ...
                       " chained subscripts not allowed."));
      endif
      switch s.type
        case '()'
          error (strcat ("Invalid () indexing for assigning values", ...
                         " to a ClassificationKNN object."));
        case '{}'
          error (strcat ("Invalid {} indexing for assigning values", ...
                         " to a ClassificationKNN object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("ClassificationKNN.subsasgn: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          switch (s.subs)
            case 'ScoreTransform'
              name = "ClassificationKNN";
              this.ScoreTransform = parseScoreTransform (val, name);
            otherwise
              error (strcat ("ClassificationKNN.subsasgn:", ...
                             " unrecongized or read-only property: '%s'"), ...
                             s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods (Access = public)

    ## Constructor
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
      Standardize     = false;
      PredictorNames  = [];
      ResponseName    = [];
      ClassNames      = [];
      Prior           = [];
      Cost            = [];
      Scale           = [];     # Distance scale for 'seuclidean'
      Cov             = [];     # Covariance matrix for 'mahalanobis'
      Exponent        = [];     # Exponent for 'minkowski'
      BreakTies       = [];
      NumNeighbors    = [];
      Distance        = [];
      DistanceWeight  = [];
      DistParameter   = [];
      NSMethod        = [];
      IncludeTies     = false;
      BucketSize      = 50;
      this.ScoreTransform  = 'none';

      ## Number of parameters for Standardize, Scale, Cov (maximum 1 allowed)
      SSC = 0;

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "standardize"
            if (SSC < 1)
              Standardize = varargin{2};
              if (! (Standardize == true || Standardize == false))
                error (strcat ("ClassificationKNN: 'Standardize' must", ...
                               " be either true or false."));
              endif
              SSC += 1;
            else
              error (strcat ("ClassificationKNN: 'Standardize' cannot", ...
                             " simultaneously be specified with either", ...
                             " Scale or Cov."));
            endif

          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("ClassificationKNN: 'PredictorNames' must", ...
                             " be supplied as a cellstring array."));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat ("ClassificationKNN: 'PredictorNames' must", ...
                             " have the same number of columns as X."));
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("ClassificationKNN: 'ResponseName'", ...
                             " must be a character vector."));
            endif

          case "classnames"
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames) ||
                   islogical (ClassNames) || ischar (ClassNames)))
              error (strcat ("ClassificationKNN: 'ClassNames' must be a", ...
                             " cell array of character vectors, a logical", ...
                             " vector, a numeric vector, or a character array."));
            endif
            ## Check that all class names are available in gnY
            if (iscellstr (ClassNames) || ischar (ClassNames))
              ClassNames = cellstr (ClassNames);
              if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                   ClassNames, "UniformOutput", false))))
                error (strcat ("ClassificationKNN: not all 'ClassNames'", ...
                               " are present in Y."));
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, "UniformOutput", false))))
                error (strcat ("ClassificationKNN: not all 'ClassNames'", ...
                               " are present in Y."));
              endif
            endif

          case "prior"
            Prior = varargin{2};
            if (! ((isnumeric (Prior) && isvector (Prior)) ||
                  (strcmpi (Prior, "empirical") || strcmpi (Prior, "uniform"))))
              error (strcat ("ClassificationKNN: 'Prior' must be either", ...
                             " a numeric vector or a character vector."));
            endif

          case "cost"
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error (strcat ("ClassificationKNN: 'Cost' must be", ...
                             " a numeric square matrix."));
            endif

          case "scoretransform"
            name = "ClassificationKNN";
            this.ScoreTransform = parseScoreTransform (varargin{2}, name);

          case "breakties"
            BreakTies = varargin{2};
            if (! ischar (BreakTies))
              error (strcat ("ClassificationKNN: 'BreakTies'", ...
                             " must be a character vector."));
            endif
            ## Check that all class names are available in gnY
            BTs = {"smallest", "nearest", "random"};
            if (! any (strcmpi (BTs, BreakTies)))
              error ("ClassificationKNN: invalid value for 'BreakTies'.");
            endif

          case "numneighbors"
            NumNeighbors = varargin{2};
            if (! (isnumeric (NumNeighbors) && isscalar (NumNeighbors) &&
                   NumNeighbors > 0 && fix (NumNeighbors) == NumNeighbors))
              error (strcat ("ClassificationKNN: 'NumNeighbors'", ...
                             " must be a positive integer."));
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
                error (strcat ("ClassificationKNN: invalid function", ...
                               " handle for distance metric."));
              end_try_catch
              Yrows = rows (Y);
              if (! isequal (size (D2), [Yrows, 1]))
                error (strcat ("ClassificationKNN: custom distance", ...
                               " function produces wrong output size."));
              endif
            else
              error ("ClassificationKNN: invalid distance metric.");
            endif

          case "distanceweight"
            DistanceWeight = varargin{2};
            DMs = {"equal", "inverse", "squareinverse"};
            if (is_function_handle (DistanceWeight))
              m = eye (5);
              if (! isequal (size (m), size (DistanceWeight (m))))
                error (strcat ("ClassificationKNN: function handle for", ...
                               " distance weight must return the same", ...
                               " size as its input."));
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
                error ("ClassificationKNN: 'Scale' must be a numeric vector.");
              endif
              SSC += 1;
            else
              error (strcat ("ClassificationKNN: 'Scale' cannot", ...
                             " simultaneously be specified with either", ...
                             " 'Standardize' or 'Cov'."));
            endif

          case "cov"
            if (SSC < 1)
              Cov = varargin{2};
              [~, p] = chol (Cov);
              if (p != 0)
                error (strcat ("ClassificationKNN: 'Cov' must be a", ...
                               " symmetric positive definite matrix."));
              endif
              SSC += 1;
            else
              error (strcat ("ClassificationKNN: 'Cov' cannot", ...
                             " simultaneously be specified with either", ...
                             " 'Standardize' or 'Scale'."));
            endif

          case "exponent"
            Exponent = varargin{2};
            if (! (isnumeric (Exponent) && isscalar (Exponent) &&
                           Exponent > 0 && fix (Exponent) == Exponent))
              error (strcat ("ClassificationKNN: 'Exponent'", ...
                             " must be a positive integer."));
            endif

          case "nsmethod"
            NSMethod = varargin{2};
            NSM = {"kdtree", "exhaustive"};
            if (! ischar (NSMethod))
              error (strcat ("ClassificationKNN: 'NSMethod' must", ...
                             " be a character vector."));
            endif
            if (! any (strcmpi (NSM, NSMethod)))
              error (strcat ("ClassificationKNN: 'NSMethod' must", ...
                             " be either 'kdtree' or 'exhaustive'."));
            endif

          case "includeties"
            IncludeTies = varargin{2};
            if (! (IncludeTies == true || IncludeTies == false))
              error (strcat ("ClassificationKNN: 'IncludeTies'", ...
                             " must be either true or false."));
            endif

          case "bucketsize"
            BucketSize = varargin{2};
            if (! (isnumeric (BucketSize) && isscalar (BucketSize) &&
                           BucketSize > 0 && fix (BucketSize) == BucketSize))
              error (strcat ("ClassificationKNN: 'BucketSize'", ...
                             " must be a positive integer."));
            endif

          otherwise
            error (strcat ("ClassificationKNN: invalid parameter",...
                           " name in optional pair arguments."));

        endswitch
        varargin (1:2) = [];
      endwhile

      ## Generate default predictors and response variabe names (if necessary)
      NumPredictors = columns (X);
      if (isempty (PredictorNames))
        for i = 1:NumPredictors
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = "Y";
      endif

      ## Assign predictors and response variable names
      this.NumPredictors  = NumPredictors;
      this.PredictorNames = PredictorNames;
      this.ResponseName   = ResponseName;

      ## Handle class names
      if (! isempty (ClassNames))
        if (iscellstr (ClassNames))
          ru = find (! ismember (gnY, ClassNames));
        else
          ru = find (! ismember (glY, ClassNames));
        endif
        for i = 1:numel (ru)
          gY(gY == ru(i)) = NaN;
        endfor
      endif

      ## Remove missing values from X and Y
      RowsUsed  = ! logical (sum (isnan ([X, gY]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

      ## Renew groups in Y, get classes ordered, keep the same type
      [this.ClassNames, gnY, gY] = unique (Y);

      ## Check X contains valid data
      if (! (isnumeric (X) && isfinite (X)))
        error ("ClassificationKNN: invalid values in X.");
      endif

      ## Assign the number of observations and their corresponding indices
      ## on the original data, which will be used for training the model,
      ## to the ClassificationKNN object
      this.NumObservations = sum (RowsUsed);
      this.RowsUsed = RowsUsed;

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
      if (strcmpi ("uniform", Prior))
        this.Prior = ones (size (gnY)) ./ numel (gnY);
      elseif (isempty (Prior) || strcmpi ("empirical", Prior))
        pr = [];
        for i = 1:numel (gnY)
          pr = [pr; sum(gY==i)];
        endfor
        this.Prior = pr ./ sum (pr);
      elseif (isnumeric (Prior))
        if (numel (gnY) != numel (Prior))
          error (strcat ("ClassificationKNN: the elements in 'Prior'", ...
                         " do not correspond to selected classes in Y."));
        endif
        this.Prior = Prior ./ sum (Prior);
      endif
      if (isempty (Cost))
        this.Cost = cast (! eye (numel (gnY)), "double");
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (strcat ("ClassificationKNN: the number of rows", ...
                         " and columns in 'Cost' must correspond", ...
                         " to selected classes in Y."));
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
          error (strcat ("ClassificationKNN: 'Scale' is only valid", ...
                         " when distance metric is seuclidean."));
        endif
        if (numel (Scale) != NumPredictors)
          error (strcat ("ClassificationKNN: 'Scale' vector must have", ...
                         " equal length to the number of columns in X."));
        endif
        if (any (Scale < 0))
          error (strcat ("ClassificationKNN: 'Scale' vector must", ...
                         " contain nonnegative scalar values."));
        endif
        this.DistParameter = Scale;
      else
        if (strcmpi (Distance, "seuclidean"))
          if (Standardize)
            this.DistParameter = ones (1, NumPredictors);
          else
            this.DistParameter = std (X, [], 1);
          endif
        endif
      endif
      if (! isempty (Cov))
        if (! strcmpi (Distance, "mahalanobis"))
          error (strcat ("ClassificationKNN: 'Cov' is only valid", ...
                         " when distance metric is 'mahalanobis'."));
        endif
        if (columns (Cov) != NumPredictors)
          error (strcat ("ClassificationKNN: 'Cov' matrix", ...
                         " must have equal columns as X."));
        endif
        this.DistParameter = Cov;
      else
        if (strcmpi (Distance, "mahalanobis"))
          this.DistParameter = cov (X);
        endif
      endif
      if (! isempty (Exponent))
        if (! strcmpi (Distance, "minkowski"))
          error (strcat ("ClassificationKNN: 'Exponent' is only", ...
                         " valid when distance metric is 'minkowski'."));
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
          error (strcat ("ClassificationKNN: 'kdtree' method is only va", ...
                         "lid for 'euclidean', 'cityblock', 'manhattan',", ...
                         " 'minkowski', and 'chebychev' distance metrics."));
        endif
        this.NSMethod = NSMethod;
      else
        if (any (strcmpi (kdm, Distance)) && NumPredictors <= 10)
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
    ## @deftypefn  {ClassificationKNN} {@var{labels} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationKNN} {[@var{labels}, @var{scores}, @var{cost}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the kNN algorithm from a
    ## k-Nearest Neighbor classification model.
    ##
    ## @code{@var{labels} = predict (@var{obj}, @var{XC})} returns the matrix of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## predictor data in @code{obj.X} and corresponding labels, @code{obj.Y},
    ## stored in the k-Nearest Neighbor classification model, @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationKNN} class object.
    ## @item
    ## @var{XC} must be an @math{MxP} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the SVM model in
    ## @var{obj}.
    ## @end itemize
    ##
    ## @code{[@var{labels}, @var{scores}, @var{cost}] = predict (@var{obj},
    ## @var{XC})} also returns @var{scores}, which contains the predicted class
    ## scores or posterior probabilities for each instance of the corresponding
    ## unique classes, and @var{cost}, which is a matrix containing the expected
    ## cost of the classifications.  By default, @var{scores} returns the
    ## posterior probabilities for KNN models, unless a specific ScoreTransform
    ## function has been specified.  See @code{fitcknn} for more info.
    ##
    ## Note! @code{predict} is explicitly using @qcode{"exhaustive"} as the
    ## nearest search method due to the very slow implementation of
    ## @qcode{"kdtree"} in the @code{knnsearch} function.
    ##
    ## @seealso{fitcknn, ClassificationKNN, knnsearch}
    ## @end deftypefn
    function [labels, scores, cost] = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("ClassificationKNN.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("ClassificationKNN.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("ClassificationKNN.predict:", ...
                       " XC must have the same number of", ...
                       " predictors as the trained model."));
      endif

      ## Get training data and labels
      X = this.X(this.RowsUsed,:);
      Y = this.Y(this.RowsUsed,:);

      ## Standardize (if necessary)
      if (this.Standardize)
        X = (X - this.Mu) ./ this.Sigma;
        XC = (XC - this.Mu) ./ this.Sigma;
      endif

      ## Train kNN
      if (strcmpi (this.Distance, "seuclidean"))
        [idx, dist] = knnsearch (X, XC, "k", this.NumNeighbors, ...
                      "NSMethod", "exhaustive", "Distance", "seuclidean", ...
                      "Scale", this.DistParameter, "sortindices", true, ...
                      "includeties", this.IncludeTies, ...
                      "bucketsize", this.BucketSize);

      elseif (strcmpi (this.Distance, "mahalanobis"))
        [idx, dist] = knnsearch (X, XC, "k", this.NumNeighbors, ...
                      "NSMethod", "exhaustive", "Distance", "mahalanobis", ...
                      "cov", this.DistParameter, "sortindices", true, ...
                      "includeties", this.IncludeTies, ...
                      "bucketsize", this.BucketSize);

      elseif (strcmpi (this.Distance, "minkowski"))
        [idx, dist] = knnsearch (X, XC, "k", this.NumNeighbors, ...
                      "NSMethod", "exhaustive", "Distance", "minkowski", ...
                      "P", this.DistParameter, "sortindices", true, ...
                      "includeties",this.IncludeTies, ...
                      "bucketsize", this.BucketSize);

      else
        [idx, dist] = knnsearch (X, XC, "k", this.NumNeighbors, ...
                      "NSMethod", "exhaustive", "Distance", this.Distance, ...
                      "sortindices", true, "includeties", this.IncludeTies, ...
                      "bucketsize", this.BucketSize);
      endif

      ## Make prediction
      if (iscellstr (this.ClassNames))
        labels = {};
      elseif (ischar (this.ClassNames))
        labels = '';
      else
        labels = [];
      endif
      scores = [];
      cost  = [];

      ## Get IDs of labels for each point in training data
      [~, ~, gY] = unique (Y);

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
        for c = 1:rows (this.ClassNames)
          freq(c) = sum (kNNgY == c) / k;
        endfor

        ## Get labels according to BreakTies
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
        labels = [labels; this.ClassNames(idl,:)];

        ## Calculate scores and cost
        scores = [scores; freq];
        cost = [cost; 1-freq];

        ## Apply ScoreTransform (if applicable)
        if (! strcmp (this.ScoreTransform, "none"))
          scores = this.ScoreTransform (scores);
        endif

      endfor

      ## Convert double to logical if ClassNames are logical
      if (islogical (this.ClassNames))
        labels = logical (labels);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKNN} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationKNN} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute loss for a trained ClassificationKNN object.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} computes the loss,
    ## @var{L}, using the default loss function @qcode{'mincost'}.
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{ClassificationKNN} object trained on @code{X} and
    ## @code{Y}.
    ## @item
    ## @code{X} must be a @math{NxP} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}.  @var{Y} must have same
    ## numbers of Rows as @var{X}.
    ## @end itemize
    ##
    ## @code{@var{L} = loss (@dots{}, @var{name}, @var{value})} allows
    ## additional options specified by @var{name}-@var{value} pairs:
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"LossFun"} @tab @tab Specifies the loss function to use.
    ## Can be a function handle with four input arguments (C, S, W, Cost)
    ## which returns a scalar value or one of:
    ## 'binodeviance', 'classifcost', 'classiferror', 'exponential',
    ## 'hinge', 'logit','mincost', 'quadratic'.
    ## @itemize
    ## @item
    ## @code{C} is a logical matrix of size @math{NxK}, where @math{N} is the
    ## number of observations and @math{K} is the number of classes.
    ## The element @code{C(i,j)} is true if the class label of the i-th
    ## observation is equal to the j-th class.
    ## @item
    ## @code{S} is a numeric matrix of size @math{NxK}, where each element
    ## represents the classification score for the corresponding class.
    ## @item
    ## @code{W} is a numeric vector of length @math{N}, representing
    ## the observation weights.
    ## @item
    ## @code{Cost} is a @math{KxK} matrix representing the misclassification
    ## costs.
    ## @end itemize
    ##
    ## @item @qcode{"Weights"} @tab @tab Specifies observation weights, must be
    ## a numeric vector of length equal to the number of rows in X.
    ## Default is @code{ones (size (X, 1))}. loss normalizes the weights so that
    ## observation weights in each class sum to the prior probability of that
    ## class. When you supply Weights, loss computes the weighted
    ## classification loss.
    ##
    ## @end multitable
    ##
    ## @seealso{fitcknn, ClassificationKNN}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationKNN.loss: too few input arguments.");
      elseif (mod (nargin - 3, 2) != 0)
        error (strcat ("ClassificationKNN.loss: name-value", ...
                       " arguments must be in pairs."));
      elseif (nargin > 7)
        error ("ClassificationKNN.loss: too many input arguments.");
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationKNN.loss: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat ("ClassificationKNN.loss: X must have the same", ...
                       " number of predictors as the trained model."));
      endif

      ## Default values
      LossFun = 'mincost';
      Weights = [];

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error ("ClassificationKNN.loss: Y must be of a valid type.");
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat ("ClassificationKNN.loss: Y must have", ...
                       " the same number of rows as X."));
      endif

      ## Parse name-value arguments
      while (numel (varargin) > 0)
        Value = varargin{2};
        switch (tolower (varargin{1}))
          case 'lossfun'
            if (isa (Value, 'function_handle'))
              ## Check if the loss function is valid
              if (nargin (Value) != 4)
                error (strcat ("ClassificationKNN.loss: custom loss function", ...
                               " must accept exactly four input arguments."));
              endif
              try
                n = 1;
                K = 2;
                C_test = false (n, K);
                S_test = zeros (n, K);
                W_test = ones (n, 1);
                Cost_test = ones (K) - eye (K);
                test_output = Value (C_test, S_test, W_test, Cost_test);
                if (! isscalar (test_output))
                  error (strcat ("ClassificationKNN.loss: custom loss", ...
                                 " function must return a scalar value."));
                endif
              catch
                error (strcat ("ClassificationKNN.loss: custom loss", ...
                               " function is not valid or does not", ...
                               " produce correct output."));
              end_try_catch
              LossFun = Value;
            elseif (ischar (Value) && any (strcmpi (Value, {"binodeviance", ...
                "classifcost", "classiferror", "exponential", "hinge", ...
                "logit", "mincost", "quadratic"})))
              LossFun = Value;
            else
              error ("ClassificationKNN.loss: invalid loss function.");
            endif
          case 'weights'
            if (isnumeric (Value) && isvector (Value))
              if (numel (Value) != size (X ,1))
                error ("ClassificationKNN.loss: size of Weights must", ...
                       " be equal to the number of rows in X.");
              elseif (numel (Value) == size (X, 1))
                Weights = Value;
              endif
            else
              error ("ClassificationKNN.loss: invalid Weights.");
            endif
          otherwise
            error ("ClassificationKNN.loss: invalid name-value arguments.");
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Check for missing values in X
      if (! isa (LossFun, 'function_handle'))
        lossfun = tolower (LossFun);
        if (! strcmp (lossfun, 'mincost') && ! strcmp (lossfun, ...
         'classiferror') && ! strcmp (lossfun, 'classifcost') ...
          && any (isnan (X(:))))
            L = NaN;
            return;
        endif
      endif

      ## If Y is a char array convert it to a cell array of character vectors
      classes = this.ClassNames;
      if (ischar (Y) && ischar (classes))
        Y = cellstr (Y);
        classes = cellstr (classes);
      endif

      ## Check that Y is of the same type as ClassNames
      if (! strcmp (class (Y), class (classes)))
        error (strcat ("ClassificationKNN.loss: Y must be the", ...
                       " same data type as the model's ClassNames."));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), this.ClassNames)))
        error (strcat ("ClassificationKNN.loss: Y must contain only", ...
                       " the classes in model's ClassNames."));
      endif

      ## Set default weights if not specified
      if (isempty (Weights))
        Weights = ones (size (X, 1), 1);
      endif

      ## Normalize Weights
      K = numel (classes);
      class_prior_probs = this.Prior;
      norm_weights = zeros (size (Weights));
      for i = 1:K
        class_idx = ismember (Y, classes(i));
        if (sum (Weights(class_idx)) > 0)
          norm_weights(class_idx) = ...
          Weights(class_idx) * class_prior_probs(i) / sum (Weights(class_idx));
        endif
      endfor
      Weights = norm_weights / sum (norm_weights);

      ## Number of observations
      n = size (X, 1);

      ## Predict classification scores
      [label, scores] = predict (this, X);
      if (ischar (label))
        label = cellstr (label);
      endif

      ## C is vector of K-1 zeros, with 1 in the
      ## position corresponding to the true class
      C = false (n, K);
      for i = 1:n
        class_idx = find (ismember (classes, Y(i)));
        C(i, class_idx) = true;
      endfor
      Y_new = C';

      ## Compute the loss using custom loss function
      if (isa (LossFun, 'function_handle'))
        L = LossFun (C, scores, Weights, this.Cost);
        return;
      endif

      ## Compute the scalar classification score for each observation
      m_j = zeros (n, 1);
      for i = 1:n
        m_j(i) = scores(i,:) * Y_new(:,i);
      endfor

      ## Compute the loss
      switch (tolower (LossFun))
        case 'binodeviance'
          b = log (1 + exp (-2 * m_j));
          L = (Weights') * b;
        case 'hinge'
          h = max (0, 1 - m_j);
          L = (Weights') * h;
        case 'exponential'
          e = exp (-m_j);
          L = (Weights') * e;
        case 'logit'
          l = log (1 + exp (-m_j));
          L = (Weights') * l;
        case 'quadratic'
          q = (1 - m_j) .^ 2;
          L = (Weights') * q;
        case 'classiferror'
          L = 0;
          for i = 1:n
            L = L + Weights(i) * (! isequal (Y(i), label(i)));
          endfor
        case 'mincost'
          Cost = this.Cost;
          L = 0;
          for i = 1:n
            f_Xj = scores(i, :);
            gamma_jk = f_Xj * Cost;
            [~, min_cost_class] = min (gamma_jk);
            cj = Cost(find (ismember (classes, Y(i))), min_cost_class);
            L = L + Weights(i) * cj;
          endfor
        case 'classifcost'
          Cost = this.Cost;
          L = 0;
          for i = 1:n
            y_idx = find (ismember (classes, Y(i)));
            y_hat_idx = find (ismember (classes, label(i)));
            L = L + Weights(i) * Cost(y_idx, y_hat_idx);
          endfor
        otherwise
          error ("ClassificationKNN.loss: invalid loss function.");
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationKNN} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns
    ## the classification margins for @var{obj} with data @var{X} and
    ## classification @var{Y}. @var{m} is a numeric vector of length size (X,1).
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{ClassificationKNN} object trained on @code{X}
    ## and @code{Y}.
    ## @item
    ## @code{X} must be a @math{NxP} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of Rows as @var{X}.
    ## @end itemize
    ##
    ## The classification margin for each observation is the difference between
    ## the classification score for the true class and the maximal
    ## classification score for the false classes.
    ##
    ## @seealso{fitcknn, ClassificationKNN}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationKNN.margin: too few input arguments.");
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationKNN.margin: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat ("ClassificationKNN.margin: X must have the same", ...
                       " number of predictors as the trained model."));
      endif

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error ("ClassificationKNN.margin: Y must be of a valid type.");
      endif

      ## Validate X
      valid_types = {'single', 'double'};
      if (! (any (strcmp (class (X), valid_types))))
        error ("ClassificationKNN.margin: X must be of a valid type.");
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat ("ClassificationKNN.margin: Y must have", ...
                       " the same number of rows as X."));
      endif

      ## If Y is a char array convert it to a cell array of character vectors
      classes = this.ClassNames;
      if (ischar (Y) && ischar (classes))
        Y = cellstr (Y);
        classes = cellstr (classes);
      endif

      ## Check that Y is of the same type as ClassNames
      if (! strcmp (class (Y), class (classes)))
        error (strcat ("ClassificationKNN.margin: Y must be the", ...
                       " same data type as the model's ClassNames."));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), classes)))
        error (strcat ("ClassificationKNN.margin: Y must contain", ...
                       " only the classes in model's ClassNames."));
      endif

      ## Number of Observations
      n = size (X, 1);

      ## Initialize the margin vector
      m = zeros (n, 1);

      ## Calculate the classification scores
      [~, scores] = predict (this, X);

      ## Loop over each observation to compute the margin
      for i = 1:n
        ## True class index
        true_class_idx = find (ismember (classes, Y(i)));

        ## Score for the true class
        true_class_score = scores(i, true_class_idx);

        ## Get the maximal score for the false classes
        scores(i, true_class_idx) = -Inf;              ## Temporarily
        max_false_class_score = max (scores(i, :));
        if (max_false_class_score == -Inf)
          m = NaN;
          return;
        endif
        scores(i, true_class_idx) = true_class_score;  ## Restore

        ## Calculate the margin
        m(i) = true_class_score - max_false_class_score;
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKNN} {@var{[pd, x, y]} =} partialDependence (@var{obj}, @var{Vars}, @var{Labels})
    ## @deftypefnx {ClassificationKNN} {@var{[pd, x, y]} =} partialDependence (@dots{}, @var{Data})
    ## @deftypefnx {ClassificationKNN} {@var{[pd, x, y]} =} partialDependence (@dots{}, @var{name}, @var{value})
    ##
    ## Compute partial dependence for a trained ClassificationKNN object.
    ##
    ## @code{@var{[pd, x, y]} = partialDependence (@var{obj}, @var{Vars},
    ## @var{Labels})}
    ## computes the partial dependence of the classification scores on the
    ## variables @var{Vars} for the specified class @var{Labels}.
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a trained @var{ClassificationKNN} object.
    ## @item
    ## @code{Vars} is a vector of positive integers, character vector,
    ## string array, or cell array of character
    ## vectors representing predictor variables (it can be indices of
    ## predictor variables in @var{obj.X}).
    ## @item
    ## @code{Labels} is a character vector, logical vector, numeric vector,
    ## or cell array of character vectors representing class
    ## labels. (column vector)
    ## @end itemize
    ##
    ## @code{@var{[pd, x, y]} = partialDependence (@dots{}, @var{Data})}
    ## specifies new predictor data to use for computing the partial dependence.
    ##
    ## @code{@var{[pd, x, y]} = partialDependence (@dots{}, @var{name},
    ## @var{value})} allows additional options specified by name-value pairs:
    ##
    ## @multitable @columnfractions 0.32 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"NumObservationsToSample"} @tab @tab Number of
    ## observations to sample. Must be a positive integer. Defaults to the
    ## number of observations in the training data.
    ## @item @qcode{"QueryPoints"} @tab @tab Points at which to evaluate
    ## the partial dependence.
    ## Must be a numeric column vector, numeric two-column matrix, or
    ## cell array of character column vectors.
    ## @item @qcode{"UseParallel"} @tab @tab Logical value indicating
    ## whether to perform computations in parallel.
    ## Defaults to @code{false}.
    ## @end multitable
    ##
    ## @subheading Return Values
    ## @itemize
    ## @item @code{pd}: Partial dependence values.
    ## @item @code{x}: Query points for the first predictor variable in Vars.
    ## @item @code{y}: Query points for the second predictor variable in
    ## Vars (if applicable).
    ## @end itemize
    ##
    ## @seealso{fitcknn, ClassificationKNN}
    ## @end deftypefn

    function [pd, x, y] = partialDependence (this, Vars, Labels, varargin)
      if (nargin < 3)
        error ("ClassificationKNN.partialDependence: too few input arguments.");
      endif

      ## Validate Vars
      if (isnumeric (Vars))
        if (! all (Vars > 0) || ! (numel (Vars) == 1 || numel (Vars) == 2))
          error ("ClassificationKNN.partialDependence: VARS must be a", ...
                 " positive integer or vector of two positive integers.");
        endif
      elseif (iscellstr (Vars))
        if (! (numel (Vars) == 1 || numel (Vars) == 2))
          error (strcat ("ClassificationKNN.partialDependence: VARS must", ...
                         " be a string array or cell array of one or two", ...
                         " character vectors."));
        endif
        Vars = cellfun (@(v) find (strcmp (this.PredictorNames, v)), Vars);
      elseif (ischar (Vars))
        Vars = find (strcmp (this.PredictorNames, Vars));
        if (isempty (Vars))
          error (strcat ("ClassificationKNN.partialDependence: VARS", ...
                         " must match one of the predictor names."));
        endif
      else
        error (strcat ("ClassificationKNN.partialDependence: VARS", ...
                       " must be a string, or cell array."));
      endif

      ## Validate Labels
      if (! (ischar (Labels) || islogical (Labels) || ...
          isnumeric (Labels) || iscellstr (Labels) || islogical (Labels)))
        error ("ClassificationKNN.partialDependence: invalid type for LABELS.");
      endif

      ## If Labels is a char array convert it to a cell array of character vectors
      classes = this.ClassNames;
      if (ischar (Labels) && ischar (classes))
        Labels = cellstr (Labels);
        classes = cellstr (classes);
      endif

      ## Check that Y is of the same type as ClassNames
      if (! strcmp (class (Labels), class (classes)))
        error (strcat ("ClassificationKNN.margin: LABELS must be the", ...
                       " same data type as the model's ClassNames."));
      endif

      ## Additional validation to match ClassNames
      if (! all (ismember (Labels, classes)))
        error (strcat ("ClassificationKNN.partialDependence: LABELS must", ...
                       " match the class names in the model's ClassNames."));
      endif

      ## Default values
      Data = this.X;
      UseParallel = false;
      NumObservationsToSample = size (Data, 1);
      QueryPoints = [];

      ## Check for Data and other optional arguments
      if (nargin > 3)
        if (size (varargin{1}) == size (this.X))
          Data = varargin{1};
          ## Ensure Data consistency
          if (! all (size (Data, 2) == numel (this.PredictorNames)))
            error (strcat ("ClassificationKNN.partialDependence: DATA must", ...
                           " have the same number and order of columns as", ...
                           " the predictor variables."));
          endif

          ## Ensure Name-Value pairs are even length
          if (mod (nargin - 4, 2) != 0)
            error (strcat ("ClassificationKNN.partialDependence:", ...
                           " name-value arguments must be in pairs."));
          endif

          ## Set the number of observations to sample
          NumObservationsToSample = size (Data, 1);
          idx = 2;
        else
          ## Ensure Name-Value pairs are even length
          if (mod (nargin - 3, 2) != 0)
            error (strcat ("ClassificationKNN.partialDependence:", ...
                           " name-value arguments must be in pairs."));
          endif
          idx = 1;
        endif

        ## Handle name-value pair arguments
        for i = idx:2:length (varargin)
          if (! ischar (varargin{i}))
            error (strcat ("ClassificationKNN.partialDependence: name", ...
                           " arguments must be strings."));
          endif
          Value = varargin{i+1};
          ## Parse name-value pairs
          switch (lower (varargin{i}))
            case 'numobservationstosample'
              if (! isnumeric (Value) || Value <= 0 || Value != round (Value))
                error (strcat ("ClassificationKNN.partialDependence:", ...
                               " NumObservationsToSample must be a", ...
                               " positive integer."));
              endif
              NumObservationsToSample = Value;
              if (Value > size (Data, 1))
                NumObservationsToSample = size (Data, 1);
              endif
            case 'querypoints'
              if (! isnumeric (Value) && ! iscell (Value))
                error (strcat ("ClassificationKNN.partialDependence:", ...
                               " QueryPoints must be a numeric column", ...
                               " vector, numeric two-column matrix, or", ...
                               " cell array of character column vectors."));
              endif
              QueryPoints = Value;
            case 'useparallel'
              if (! islogical (UseParallel))
                error (strcat ("ClassificationKNN.partialDependence:", ...
                               " UseParallel must be a logical value."));
              endif
              UseParallel = Value;
            otherwise
              error (strcat ("ClassificationKNN.partialDependence:", ...
                             " name-value pair argument not recognized."));
          endswitch
        endfor
      endif

      ## Sample observations if needed
      if (NumObservationsToSample < size (Data, 1))
        Data = datasample (Data, NumObservationsToSample, 'Replace', false);
      endif

      ## Generate QueryPoints if not specified
      if (isempty (QueryPoints))
        if (numel (Vars) == 1)
          if (isnumeric (Data(:, Vars)))
            QueryPoints = linspace(min (Data(:, Vars)), ...
                                max (Data(:, Vars)), 100)';
          else
            QueryPoints = unique (Data(:, Vars));
          endif
        else
          QueryPoints = cell (1, numel (Vars));
          for j = 1:numel (Vars)
            if (isnumeric (Data(:, Vars(j))))
              QueryPoints{j} = linspace(min (Data(:, Vars(j))), ...
                                max (Data(:, Vars(j))), 100)';
            else
              QueryPoints{j} = unique (Data(:, Vars(j)));
            endif
          endfor
        endif
      endif

      ## Prepare grid points for predictions
      if (numel (Vars) == 1)
        gridPoints = QueryPoints;
      else
        if (ischar (QueryPoints))
          [X1, X2] = meshgrid (QueryPoints(1), QueryPoints(2));
        else
          [X1, X2] = meshgrid (QueryPoints{1}, QueryPoints{2});
        endif
        gridPoints = [X1(:), X2(:)];
      endif

      ## Predict responses for the grid points
      numClasses = numel (classes);
      numQueryPoints = size (gridPoints, 1);
      predictions = zeros (numQueryPoints, numClasses);

      if (UseParallel)
        parfor i = 1:numQueryPoints
          tempData = Data;
          for j = 1:numel (Vars)
            tempData(:, Vars(j)) = repmat (gridPoints(i, j), ...
                                    NumObservationsToSample, 1);
          endfor
          [~, scores] = predict (this, tempData);
          predictions(i, :) = mean (scores, 1);
        endparfor
      else
        for i = 1:numQueryPoints
          tempData = Data;
          for j = 1:numel (Vars)
            tempData(:, Vars(j)) = repmat (gridPoints(i, j), ...
                                    NumObservationsToSample, 1);
          endfor
          [~, scores] = predict (this, tempData);
          predictions(i, :) = mean (scores, 1);
        endfor
      endif

      ## Compute partial dependence
      if (numel (Vars) == 1)
        if (numel (Labels) == 1)
          classIndex = find (ismember (classes, Labels));
          pd = predictions(:, classIndex)';
        else
          pd = zeros (numel (Labels), numel (QueryPoints));
          for j = 1:numel (Labels)
            classIndex = find (ismember (classes, Labels(j)));
            pd(j, :) = predictions(:, classIndex)';
          endfor
        endif
        x = QueryPoints;
        y = [];
      else
        if (numel (Labels) == 1)
          classIndex = find (ismember (classes, Labels));
          pd = reshape (predictions(:, classIndex), numel (QueryPoints{1}), ...
                        numel (QueryPoints{2}));
        else
          pd = zeros (numel (Labels), numel (QueryPoints{1}), ...
                      numel (QueryPoints{2}));
          for j = 1:numel (Labels)
            classIndex = find (ismember (classes, Labels(j)));
            pd(j, :, :) = reshape (predictions(:, classIndex), ...
                                   numel (QueryPoints{1}), ...
                                   numel (QueryPoints{2}));
          endfor
        endif
        x = QueryPoints{1};
        y = QueryPoints{2};
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKNN} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationKNN} {@var{CVMdl} =} crossval (@dots{}, @var{Name}, @var{Value})
    ##
    ## Cross Validate a ClassificationKNN object.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} returns a cross-validated model
    ## object, @var{CVMdl}, from a trained model, @var{obj}, using 10-fold
    ## cross-validation by default.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj}, @var{name}, @var{value})}
    ## specifies additional name-value pair arguments to customize the
    ## cross-validation process.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"KFold"} @tab @tab Specify the number of folds to use in
    ## k-fold cross-validation.  @code{"KFold", @var{k}}, where @var{k} is an
    ## integer greater than 1.
    ##
    ## @item @qcode{"Holdout"} @tab @tab Specify the fraction of the data to
    ## hold out for testing.  @code{"Holdout", @var{p}}, where @var{p} is a
    ## scalar in the range @math{(0,1)}.
    ##
    ## @item @qcode{"Leaveout"} @tab @tab Specify whether to perform
    ## leave-one-out cross-validation.  @code{"Leaveout", @var{Value}}, where
    ## @var{Value} is 'on' or 'off'.
    ##
    ## @item @qcode{"CVPartition"} @tab @tab Specify a @qcode{cvpartition}
    ## object used for cross-validation.  @code{"CVPartition", @var{cv}}, where
    ## @code{isa (@var{cv}, "cvpartition")} = 1.
    ##
    ## @end multitable
    ##
    ## @seealso{fitcknn, ClassificationKNN, cvpartition,
    ## ClassificationPartitionedModel}
    ## @end deftypefn

    function CVMdl = crossval (this, varargin)
      ## Check input
      if (nargin < 1)
        error ("ClassificationKNN.crossval: too few input arguments.");
      endif

      if (numel (varargin) == 1)
        error (strcat ("ClassificationKNN.crossval: Name-Value", ...
                       " arguments must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("ClassificationKNN.crossval: specify only one", ...
                       " of the optional Name-Value paired arguments."));
      endif

      ## Add default values
      numFolds    = 10;
      Holdout     = [];
      Leaveout    = 'off';
      CVPartition = [];

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'kfold'
            numFolds = varargin{2};
            if (! (isnumeric (numFolds) && isscalar (numFolds)
                   && (numFolds == fix (numFolds)) && numFolds > 1))
              error (strcat ("ClassificationKNN.crossval: 'KFold' must", ...
                             " be an integer value greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("ClassificationKNN.crossval: 'Holdout' must", ...
                             " be a numeric value between 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat ("ClassificationKNN.crossval: 'Leaveout'", ...
                             " must be either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("ClassificationKNN.crossval: 'CVPartition'",...
                             " must be a 'cvpartition' object."));
            endif

          otherwise
            error (strcat ("ClassificationKNN.crossval: invalid",...
                           " parameter name in optional paired arguments."));
          endswitch
        varargin (1:2) = [];
      endwhile

      ## Determine the cross-validation method to use
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (this.Y, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (this.Y, 'LeaveOut');
      else
        partition = cvpartition (this.Y, 'KFold', numFolds);
      endif

      ## Create a cross-validated model object
      CVMdl = ClassificationPartitionedModel (this, partition);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKNN} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationKNN object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## ClassificationKNN object into an Octave binary file, the name of which is
    ## specified in @var{filename}, along with an extra variable, which defines
    ## the type classification object these variables constitute.  Use
    ## @code{loadmodel} in order to load a classification object into Octave's
    ## workspace.
    ##
    ## @seealso{loadmodel, fitcknn, ClassificationKNN}
    ## @end deftypefn

    function savemodel (this, fname)
      ## Generate variable for class name
      classdef_name = "ClassificationKNN";

      ## Create variables from model properties
      X = this.X;
      Y = this.Y;
      NumObservations = this.NumObservations;
      RowsUsed        = this.RowsUsed;
      Standardize     = this.Standardize;
      Sigma           = this.Sigma;
      Mu              = this.Mu;
      NumPredictors   = this.NumPredictors;
      PredictorNames  = this.PredictorNames;
      ResponseName    = this.ResponseName;
      ClassNames      = this.ClassNames;
      Prior           = this.Prior;
      Cost            = this.Cost;
      ScoreTransform  = this.ScoreTransform;
      BreakTies       = this.BreakTies;
      NumNeighbors    = this.NumNeighbors;
      Distance        = this.Distance;
      DistanceWeight  = this.DistanceWeight;
      DistParameter   = this.DistParameter;
      NSMethod        = this.NSMethod;
      IncludeTies     = this.IncludeTies;
      BucketSize      = this.BucketSize;

      ## Save classdef name and all model properties as individual variables
      save ("-binary", fname, "classdef_name", "X", "Y", "NumObservations", ...
            "RowsUsed", "Standardize", "Sigma", "Mu", "NumPredictors", ...
            "PredictorNames", "ResponseName", "ClassNames", "Prior", "Cost", ...
            "ScoreTransform", "BreakTies", "NumNeighbors", "Distance", ...
            "DistanceWeight", "DistParameter", "NSMethod", "IncludeTies", ...
            "BucketSize");
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationKNN object
      mdl = ClassificationKNN (1, 1);

      ## Check that fieldnames in DATA match properties in ClassificationKNN
      names = fieldnames (data);
      props = fieldnames (mdl);
      if (! isequal (sort (names), sort (props)))
        error ("ClassificationKNN.load_model: invalid model in '%s'.", filename)
      endif

      ## Copy data into object
      for i = 1:numel (props)
        mdl.(props{i}) = data.(props{i});
      endfor
    endfunction

  endmethods

endclassdef

## Helper functions for ScoreTransform
function out = ismax (score)
  out = score;
  out(score == max (score)) = 1;
  out(score != max (score)) = 0;
endfunction

function out = symmetricismax (score)
  out = score;
  out(score == max (score)) = 1;
  out(score != max (score)) = -1;
endfunction

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
%! load fisheriris
%! x = meas;
%! y = species;
%! obj = fitcknn (x, y, "NumNeighbors", 5, "Standardize", 1);
%!
%! ## Create a cross-validated model
%! CVMdl = crossval (obj)

%!demo
%! load fisheriris
%! x = meas;
%! y = species;
%! covMatrix = cov (x);
%!
%! ## Fit the k-NN model using the 'mahalanobis' distance
%! ## and the custom covariance matrix
%! obj = fitcknn(x, y, 'NumNeighbors', 5, 'Distance','mahalanobis', ...
%! 'Cov', covMatrix);
%!
%! ## Create a partition model using cvpartition
%! Partition = cvpartition (size (x, 1), 'kfold', 12);
%!
%! ## Create cross-validated model using 'cvPartition' name-value argument
%! CVMdl = crossval (obj, 'cvPartition', Partition)
%!
%! ## Access the trained model from first fold of cross-validation
%! CVMdl.Trained{1}

%!demo
%! X = [1, 2; 3, 4; 5, 6];
%! Y = {'A'; 'B'; 'A'};
%! model = fitcknn (X, Y);
%! customLossFun = @(C, S, W, Cost) sum (W .* sum (abs (C - S), 2));
%! ## Calculate loss using custom loss function
%! L = loss (model, X, Y, 'LossFun', customLossFun)

%!demo
%! X = [1, 2; 3, 4; 5, 6];
%! Y = {'A'; 'B'; 'A'};
%! model = fitcknn (X, Y);
%! ## Calculate loss using 'mincost' loss function
%! L = loss (model, X, Y, 'LossFun', 'mincost')

%!demo
%! X = [1, 2; 3, 4; 5, 6];
%! Y = ['1'; '2'; '3'];
%! model = fitcknn (X, Y);
%! X_test = [3, 3; 5, 7];
%! Y_test = ['1'; '2'];
%! ## Specify custom Weights
%! W = [1; 2];
%! L = loss (model, X_test, Y_test, 'LossFun', 'logit', 'Weights', W);

%!demo
%! load fisheriris
%! mdl = fitcknn (meas, species);
%! X = mean (meas);
%! Y = {'versicolor'};
%! m = margin (mdl, X, Y)

%!demo
%! X = [1, 2; 4, 5; 7, 8; 3, 2];
%! Y = [2; 1; 3; 2];
%! ## Train the model
%! mdl = fitcknn (X, Y);
%! ## Specify Vars and Labels
%! Vars = 1;
%! Labels = 2;
%! ## Calculate partialDependence
%! [pd, x, y] = partialDependence (mdl, Vars, Labels);

%!demo
%! X = [1, 2; 4, 5; 7, 8; 3, 2];
%! Y = [2; 1; 3; 2];
%! ## Train the model
%! mdl = fitcknn (X, Y);
%! ## Specify Vars and Labels
%! Vars = 1;
%! Labels = 1;
%! queryPoints = [linspace(0, 1, 3)', linspace(0, 1, 3)'];
%! ## Calculate partialDependence using queryPoints
%! [pd, x, y] = partialDependence (mdl, Vars, Labels, 'QueryPoints', ...
%! queryPoints)

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
%!test
%! x = [1, 2; 3, 4; 5,6; 5, 8];
%! y = {'9'; '9'; '6'; '7'};
%! a = ClassificationKNN (x, y);
%! assert (a.Prior, [0.25; 0.25; 0.5])

## Test constructor with ClassNames parameter
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! ClassNames = {'setosa', 'versicolor', 'virginica'};
%! a = ClassificationKNN (x, y, 'ClassNames', ClassNames);
%! assert (a.ClassNames, ClassNames')

## Test input validation for constructor
%!error<ClassificationKNN: too few input arguments.> ClassificationKNN ()
%!error<ClassificationKNN: too few input arguments.> ...
%! ClassificationKNN (ones(4, 1))
%!error<ClassificationKNN: number of rows in X and Y must be equal.> ...
%! ClassificationKNN (ones (4,2), ones (1,4))
%!error<ClassificationKNN: 'Standardize' must be either true or false.> ...
%! ClassificationKNN (ones (5,3), ones (5,1), "standardize", "a")
%!error<ClassificationKNN: 'Standardize' cannot simultaneously be specified with> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "scale", [1 1], "standardize", true)
%!error<ClassificationKNN: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "PredictorNames", ["A"])
%!error<ClassificationKNN: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "PredictorNames", "A")
%!error<ClassificationKNN: 'PredictorNames' must have the same number of columns as X.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "PredictorNames", {"A", "B", "C"})
%!error<ClassificationKNN: 'ResponseName' must be a character vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "ResponseName", {"Y"})
%!error<ClassificationKNN: 'ResponseName' must be a character vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "ResponseName", 1)
%!error<ClassificationKNN: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationKNN (ones(10,2), ones (10,1), "ClassNames", @(x)x)
%!error<ClassificationKNN: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationKNN (ones(10,2), ones (10,1), "ClassNames", {1})
%!error<ClassificationKNN: not all 'ClassNames' are present in Y.> ...
%! ClassificationKNN (ones(10,2), ones (10,1), "ClassNames", [1, 2])
%!error<ClassificationKNN: not all 'ClassNames' are present in Y.> ...
%! ClassificationKNN (ones(5,2), ['a';'b';'a';'a';'b'], "ClassNames", ['a';'c'])
%!error<ClassificationKNN: not all 'ClassNames' are present in Y.> ...
%! ClassificationKNN (ones(5,2), {'a';'b';'a';'a';'b'}, "ClassNames", {'a','c'})
%!error<ClassificationKNN: not all 'ClassNames' are present in Y.> ...
%! ClassificationKNN (ones(10,2), logical (ones (10,1)), "ClassNames", [true, false])
%!error<ClassificationKNN: 'BreakTies' must be a character vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BreakTies", 1)
%!error<ClassificationKNN: 'BreakTies' must be a character vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BreakTies", {"1"})
%!error<ClassificationKNN: invalid value for 'BreakTies'.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BreakTies", "some")
%!error<ClassificationKNN: 'Prior' must be either a numeric vector or a character vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Prior", {"1", "2"})
%!error<ClassificationKNN: 'Cost' must be a numeric square matrix.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cost", [1, 2])
%!error<ClassificationKNN: 'Cost' must be a numeric square matrix.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cost", "string")
%!error<ClassificationKNN: 'Cost' must be a numeric square matrix.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cost", {eye(2)})
%!error<ClassificationKNN: 'NumNeighbors' must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NumNeighbors", 0)
%!error<ClassificationKNN: 'NumNeighbors' must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NumNeighbors", 15.2)
%!error<ClassificationKNN: 'NumNeighbors' must be a positive integer.> ...
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
%!error<ClassificationKNN: 'Scale' must be a numeric vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Scale", "scale")
%!error<ClassificationKNN: 'Scale' must be a numeric vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Scale", {[1 2 3]})
%!error<ClassificationKNN: 'Scale' cannot simultaneously be specified with> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "standardize", true, "scale", [1 1])
%!error<ClassificationKNN: 'Cov' must be a symmetric positive definite matrix.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cov", ones (2), "Distance", "mahalanobis")
%!error<ClassificationKNN: 'Cov' cannot simultaneously be specified with> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "scale", [1 1], "Cov", ones (2))
%!error<ClassificationKNN: 'Exponent' must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Exponent", 12.5)
%!error<ClassificationKNN: 'Exponent' must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Exponent", -3)
%!error<ClassificationKNN: 'Exponent' must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Exponent", "three")
%!error<ClassificationKNN: 'Exponent' must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Exponent", {3})
%!error<ClassificationKNN: 'NSMethod' must be a character vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NSMethod", {"kdtree"})
%!error<ClassificationKNN: 'NSMethod' must be a character vector.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NSMethod", 3)
%!error<ClassificationKNN: 'NSMethod' must be either 'kdtree' or 'exhaustive'.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "NSMethod", "some")
%!error<ClassificationKNN: 'IncludeTies' must be either true or false.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "IncludeTies", "some")
%!error<ClassificationKNN: 'BucketSize' must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BucketSize", 42.5)
%!error<ClassificationKNN: 'BucketSize' must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BucketSize", -50)
%!error<ClassificationKNN: 'BucketSize' must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BucketSize", "some")
%!error<ClassificationKNN: 'BucketSize' must be a positive integer.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "BucketSize", {50})
%!error<ClassificationKNN: invalid parameter name in optional pair arguments.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "some", "some")
%!error<ClassificationKNN: invalid values in X.> ...
%! ClassificationKNN ([1;2;3;'a';4], ones (5,1))
%!error<ClassificationKNN: invalid values in X.> ...
%! ClassificationKNN ([1;2;3;Inf;4], ones (5,1))
%!error<ClassificationKNN: the elements in 'Prior' do not correspond to selected classes in Y.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Prior", [1 2])
%!error<ClassificationKNN: the number of rows and columns in 'Cost' must correspond to selected classes in Y.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cost", [1 2; 1 3])
%!error<ClassificationKNN: 'Scale' is only valid when distance metric is seuclidean.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Scale", [1 1])
%!error<ClassificationKNN: 'Scale' vector must have equal length to the number of columns in X.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Scale", [1 1 1], "Distance", "seuclidean")
%!error<ClassificationKNN: 'Scale' vector must contain nonnegative scalar values.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Scale", [1 -1], "Distance", "seuclidean")
%!error<ClassificationKNN: 'Cov' is only valid when distance metric is 'mahalanobis'.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cov", eye (2))
%!error<ClassificationKNN: 'Cov' matrix must have equal columns as X.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Cov", eye (3), "Distance", "mahalanobis")
%!error<ClassificationKNN: 'Exponent' is only valid when distance metric is 'minkowski'.> ...
%! ClassificationKNN (ones (5,2), ones (5,1), "Exponent", 3)
%!error<ClassificationKNN: 'kdtree' method is only valid for 'euclidean', 'cityblock', 'manhattan', 'minkowski', and 'chebychev' distance metrics.> ...
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
%! xc = [5.2, 4.1, 1.5, 0.1; 5.1, 3.8, 1.9, 0.4; ...
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
%!error<ClassificationKNN.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (ClassificationKNN (ones (4,2), ones (4,1)), 1)

## Test output for loss method
%!test
%! load fisheriris
%! model = fitcknn (meas, species, 'NumNeighbors', 5);
%! X = mean (meas);
%! Y = {'versicolor'};
%! L = loss (model, X, Y);
%! assert (L, 0)
%!test
%! load fisheriris
%! model = fitcknn (meas, species, 'NumNeighbors', 5);
%! L = loss (model, meas, species, 'LossFun', 'binodeviance');
%! assert (L, 0.1413, 1e-4)
%!test
%! load fisheriris
%! model = fitcknn (meas, species);
%! L = loss (model, meas, species, 'LossFun', 'binodeviance');
%! assert (L, 0.1269, 1e-4)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = {'A'; 'B'; 'A'};
%! model = fitcknn (X, Y);
%! X_test = [1, 6; 3, 3];
%! Y_test = {'A'; 'B'};
%! L = loss (model, X_test, Y_test);
%! assert (abs (L - 0.6667) > 1e-5)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = {'A'; 'B'; 'A'};
%! model = fitcknn (X, Y);
%! X_with_nan = [1, 2; NaN, 4];
%! Y_test = {'A'; 'B'};
%! L = loss (model, X_with_nan, Y_test);
%! assert (abs (L - 0.3333) < 1e-4)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = {'A'; 'B'; 'A'};
%! model = fitcknn (X, Y);
%! X_with_nan = [1, 2; NaN, 4];
%! Y_test = {'A'; 'B'};
%! L = loss (model, X_with_nan, Y_test, 'LossFun', 'logit');
%! assert (isnan (L))
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = {'A'; 'B'; 'A'};
%! model = fitcknn (X, Y);
%! customLossFun = @(C, S, W, Cost) sum (W .* sum (abs (C - S), 2));
%! L = loss (model, X, Y, 'LossFun', customLossFun);
%! assert (L, 0)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = [1; 2; 1];
%! model = fitcknn (X, Y);
%! L = loss (model, X, Y, 'LossFun', 'classiferror');
%! assert (L, 0)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = [true; false; true];
%! model = fitcknn (X, Y);
%! L = loss (model, X, Y, 'LossFun', 'binodeviance');
%! assert (abs (L - 0.1269) < 1e-4)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = ['1'; '2'; '1'];
%! model = fitcknn (X, Y);
%! L = loss (model, X, Y, 'LossFun', 'classiferror');
%! assert (L, 0)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = ['1'; '2'; '3'];
%! model = fitcknn (X, Y);
%! X_test = [3, 3];
%! Y_test = ['1'];
%! L = loss (model, X_test, Y_test, 'LossFun', 'quadratic');
%! assert (L, 1)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = ['1'; '2'; '3'];
%! model = fitcknn (X, Y);
%! X_test = [3, 3; 5, 7];
%! Y_test = ['1'; '2'];
%! L = loss (model, X_test, Y_test, 'LossFun', 'classifcost');
%! assert (L, 1)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = ['1'; '2'; '3'];
%! model = fitcknn (X, Y);
%! X_test = [3, 3; 5, 7];
%! Y_test = ['1'; '2'];
%! L = loss (model, X_test, Y_test, 'LossFun', 'hinge');
%! assert (L, 1)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = ['1'; '2'; '3'];
%! model = fitcknn (X, Y);
%! X_test = [3, 3; 5, 7];
%! Y_test = ['1'; '2'];
%! W = [1; 2];
%! L = loss (model, X_test, Y_test, 'LossFun', 'logit', 'Weights', W);
%! assert (abs (L - 0.6931) < 1e-4)

## Test input validation for loss method
%!error<ClassificationKNN.loss: too few input arguments.> ...
%! loss (ClassificationKNN (ones (4,2), ones (4,1)))
%!error<ClassificationKNN.loss: too few input arguments.> ...
%! loss (ClassificationKNN (ones (4,2), ones (4,1)), ones (4,2))
%!error<ClassificationKNN.loss: X is empty.> ...
%! loss (ClassificationKNN (ones (40,2), randi ([1, 2], 40, 1)), [], zeros (2))
%!error<ClassificationKNN.loss: X must have the same number of predictors as the trained model.> ...
%! loss (ClassificationKNN (ones (40,2), randi ([1, 2], 40, 1)), 1, zeros (2))
%!error<ClassificationKNN.loss: name-value arguments must be in pairs.> ...
%! loss (ClassificationKNN (ones (4,2), ones (4,1)), ones (4,2), ...
%!        ones (4,1), 'LossFun')
%!error<ClassificationKNN.loss: Y must have the same number of rows as X.> ...
%! loss (ClassificationKNN (ones (4,2), ones (4,1)), ones (4,2), ones (3,1))
%!error<ClassificationKNN.loss: invalid loss function.> ...
%! loss (ClassificationKNN (ones (4,2), ones (4,1)), ones (4,2), ...
%!        ones (4,1), 'LossFun', 'a')
%!error<ClassificationKNN.loss: invalid Weights.> ...
%! loss (ClassificationKNN (ones (4,2), ones (4,1)), ones (4,2), ...
%!        ones (4,1), 'Weights', 'w')

## Test output for margin method
%!test
%! load fisheriris
%! mdl = fitcknn (meas, species, 'NumNeighbors', 5);
%! X = mean (meas);
%! Y = {'versicolor'};
%! m = margin (mdl, X, Y);
%! assert (m, 1)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = [1; 2; 3];
%! mdl = fitcknn (X, Y);
%! m = margin (mdl, X, Y);
%! assert (m, [1; 1; 1])
%!test
%! X = [7, 8; 9, 10];
%! Y = ['1'; '2'];
%! mdl = fitcknn (X, Y);
%! m = margin (mdl, X, Y);
%! assert (m, [1; 1])
%!test
%! X = [11, 12];
%! Y = {'1'};
%! mdl = fitcknn (X, Y);
%! m = margin (mdl, X, Y);
%! assert (isnan (m))
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = [1; 2; 3];
%! mdl = fitcknn (X, Y);
%! X1 = [15, 16];
%! Y1 = [1];
%! m = margin (mdl, X1, Y1);
%! assert (m, -1)

## Test input validation for margin method
%!error<ClassificationKNN.margin: too few input arguments.> ...
%! margin (ClassificationKNN (ones (4,2), ones (4,1)))
%!error<ClassificationKNN.margin: too few input arguments.> ...
%! margin (ClassificationKNN (ones (4,2), ones (4,1)), ones (4,2))
%!error<ClassificationKNN.margin: X is empty.> ...
%! margin (ClassificationKNN (ones (40,2), randi ([1, 2], 40, 1)), [], zeros (2))
%!error<ClassificationKNN.margin: X must have the same number of predictors as the trained model.> ...
%! margin (ClassificationKNN (ones (40,2), randi ([1, 2], 40, 1)), 1, zeros (2))
%!error<ClassificationKNN.margin: Y must have the same number of rows as X.> ...
%! margin (ClassificationKNN (ones (4,2), ones (4,1)), ones (4,2), ones (3,1))

## Test output for partialDependence
%!shared X, Y, mdl
%! X = [1, 2; 4, 5; 7, 8; 3, 2];
%! Y = [2; 1; 3; 2];
%! mdl = fitcknn (X, Y);
%!test
%! Vars = 1;
%! Labels = 2;
%! [pd, x, y] = partialDependence (mdl, Vars, Labels);
%! pdm = [0.7500, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000];
%! assert (pd, pdm)
%!test
%! Vars = 1;
%! Labels = 2;
%! [pd, x, y] = partialDependence (mdl, Vars, Labels, ...
%! 'NumObservationsToSample', 5);
%! pdm = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
%! 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
%! 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
%! 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
%! 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%! assert (abs (pdm - pd) < 1)
%!test
%! Vars = 1;
%! Labels = 2;
%! [pd, x, y] = partialDependence (mdl, Vars, Labels, 'UseParallel', true);
%! pdm = [0.7500, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000];
%! assert (pd, pdm)
%!test
%! Vars = [1, 2];
%! Labels = 1;
%! queryPoints = {linspace(0, 1, 3)', linspace(0, 1, 3)'};
%! [pd, x, y] = partialDependence (mdl, Vars, Labels, 'QueryPoints', ...
%!                            queryPoints, 'UseParallel', true);
%! pdm = [0, 0, 0; 0, 0, 0; 0, 0, 0];
%! assert (pd, pdm)
%!test
%! Vars = 1;
%! Labels = [1; 2];
%! [pd, x, y] = partialDependence (mdl, Vars, Labels);
%! pdm = [0.2500, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.2500, 0.2500, 0.2500, ...
%! 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, ...
%! 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, ...
%! 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, ...
%! 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, ...
%! 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, 0.2500, ...
%! 0.2500, 0.2500; 0.7500, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, 0.5000, ...
%! 0.5000, 0.5000, 0.5000];
%! assert (pd, pdm)
%!test
%! Vars = [1, 2];
%! Labels = [1; 2];
%! queryPoints = {linspace(0, 1, 3)', linspace(0, 1, 3)'};
%! [pd, x, y] = partialDependence (mdl, Vars, Labels, 'QueryPoints', queryPoints);
%! pdm(:,:,1) = [0, 0, 0; 1, 1, 1];
%! pdm(:,:,2) = [0, 0, 0; 1, 1, 1];
%! pdm(:,:,3) = [0, 0, 0; 1, 1, 1];
%! assert (pd, pdm)
%!test
%! X1 = [1; 2; 4; 5; 7; 8; 3; 2];
%! X2 = ['2'; '3'; '1'; '3'; '1'; '3'; '2'; '2'];
%! X = [X1, double(X2)];
%! Y = [1; 2; 3; 3; 2; 1; 2; 1];
%! mdl = fitcknn (X, Y, 'ClassNames', {'1', '2', '3'});
%! Vars = 1;
%! Labels = 1;
%! [pd, x, y] = partialDependence (mdl, Vars, Labels);
%! pdm = [1.0000, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, ...
%! 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, ...
%! 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
%! 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
%! 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3750, ...
%! 0.3750, 0.3750, 0.3750, 0.3750, 0.3750, 0.3750, 0.3750, 0.3750, 0.3750, ...
%! 0.3750, 0.3750, 0.3750, 0.3750, 0.7500, 0.7500, 0.7500, 0.7500, 0.7500, ...
%! 0.7500, 0.7500, 0.7500];
%! assert (pd, pdm)
%!test
%! X1 = [1; 2; 4; 5; 7; 8; 3; 2];
%! X2 = ['2'; '3'; '1'; '3'; '1'; '3'; '2'; '2'];
%! X = [X1, double(X2)];
%! Y = [1; 2; 3; 3; 2; 1; 2; 1];
%! predictorNames = {'Feature1', 'Feature2'};
%! mdl = fitcknn (X, Y, 'PredictorNames', predictorNames);
%! Vars = 'Feature1';
%! Labels = 1;
%! [pd, x, y] = partialDependence (mdl, Vars, Labels);
%! pdm = [1.0000, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, ...
%! 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, ...
%! 0.6250, 0.6250, 0.6250, 0.6250, 0.6250, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
%! 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
%! 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3750, ...
%! 0.3750, 0.3750, 0.3750, 0.3750, 0.3750, 0.3750, 0.3750, 0.3750, 0.3750, ...
%! 0.3750, 0.3750, 0.3750, 0.3750, 0.7500, 0.7500, 0.7500, 0.7500, 0.7500, ...
%! 0.7500, 0.7500, 0.7500];
%! assert (pd, pdm)
%!test
%! X1 = [1; 2; 4; 5; 7; 8; 3; 2];
%! X2 = ['2'; '3'; '1'; '3'; '1'; '3'; '2'; '2'];
%! X = [X1, double(X2)];
%! Y = [1; 2; 3; 3; 2; 1; 2; 1];
%! predictorNames = {'Feature1', 'Feature2'};
%! mdl = fitcknn (X, Y, 'PredictorNames', predictorNames);
%! new_X1 = [10; 5; 6; 8; 9; 20; 35; 6];
%! new_X2 = ['2'; '2'; '1'; '2'; '1'; '3'; '3'; '2'];
%! new_X = [new_X1, double(new_X2)];
%! Vars = 'Feature1';
%! Labels = 1;
%! [pd, x, y] = partialDependence (mdl, Vars, Labels, new_X);
%! pdm = [0, 0, 0, 0, 0, 0.2500, 0.2500, 0.2500, 0.2500, 0.7500, 0.7500, ...
%! 0.7500, 0.7500, 0.7500, 0.7500, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, ...
%! 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, ...
%! 1.0000, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
%! 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
%! 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
%! assert (pd, pdm)

## Test input validation for partialDependence method
%!error<ClassificationKNN.partialDependence: too few input arguments.> ...
%! partialDependence (ClassificationKNN (ones (4,2), ones (4,1)))
%!error<ClassificationKNN.partialDependence: too few input arguments.> ...
%! partialDependence (ClassificationKNN (ones (4,2), ones (4,1)), 1)
%!error<ClassificationKNN.partialDependence: name-value arguments must be in pairs.> ...
%! partialDependence (ClassificationKNN (ones (4,2), ones (4,1)), 1, ...
%! ones (4,1), 'NumObservationsToSample')
%!error<ClassificationKNN.partialDependence: name-value arguments must be in pairs.> ...
%! partialDependence (ClassificationKNN (ones (4,2), ones (4,1)), 1, ...
%! ones (4,1), 2)

## Test output for crossval method
%!shared x, y, obj
%! load fisheriris
%! x = meas;
%! y = species;
%! covMatrix = cov (x);
%! obj = fitcknn (x, y, 'NumNeighbors', 5, 'Distance', ...
%!      'mahalanobis', 'Cov', covMatrix);
%!test
%! CVMdl = crossval (obj);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 10)
%! assert (CVMdl.ModelParameters.NumNeighbors == 5)
%! assert (strcmp (CVMdl.ModelParameters.Distance, "mahalanobis"))
%! assert (class (CVMdl.Trained{1}), "ClassificationKNN")
%! assert (!CVMdl.ModelParameters.Standardize)
%!test
%! CVMdl = crossval (obj, "KFold", 5);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 5)
%! assert (CVMdl.ModelParameters.NumNeighbors == 5)
%! assert (strcmp (CVMdl.ModelParameters.Distance, "mahalanobis"))
%! assert (class (CVMdl.Trained{1}), "ClassificationKNN")
%! assert (CVMdl.ModelParameters.Standardize == obj.Standardize)
%!test
%! obj = fitcknn (x, y, "NumNeighbors", 5, "Distance", "cityblock");
%! CVMdl = crossval (obj, "HoldOut", 0.2);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.ModelParameters.NumNeighbors == 5)
%! assert (strcmp (CVMdl.ModelParameters.Distance, "cityblock"))
%! assert (class (CVMdl.Trained{1}), "ClassificationKNN")
%! assert (CVMdl.ModelParameters.Standardize == obj.Standardize)
%!test
%! obj = fitcknn (x, y, "NumNeighbors", 10, "Distance", "cityblock");
%! CVMdl = crossval (obj, "LeaveOut", 'on');
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.ModelParameters.NumNeighbors == 10)
%! assert (strcmp (CVMdl.ModelParameters.Distance, "cityblock"))
%! assert (class (CVMdl.Trained{1}), "ClassificationKNN")
%! assert (CVMdl.ModelParameters.Standardize == obj.Standardize)
%!test
%! obj = fitcknn (x, y, "NumNeighbors", 10, "Distance", "cityblock");
%! partition = cvpartition (y, 'KFold', 3);
%! CVMdl = crossval (obj, 'cvPartition', partition);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert (CVMdl.KFold == 3)
%! assert (CVMdl.ModelParameters.NumNeighbors == 10)
%! assert (strcmp (CVMdl.ModelParameters.Distance, "cityblock"))
%! assert (class (CVMdl.Trained{1}), "ClassificationKNN")
%! assert (CVMdl.ModelParameters.Standardize == obj.Standardize)

## Test input validation for crossval method
%!error<ClassificationKNN.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (ClassificationKNN (ones (4,2), ones (4,1)), "kfold")
%!error<ClassificationKNN.crossval: specify only one of the optional Name-Value paired arguments.>...
%! crossval (ClassificationKNN (ones (4,2), ones (4,1)), "kfold", 12, "holdout", 0.2)
%!error<ClassificationKNN.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationKNN (ones (4,2), ones (4,1)), "kfold", 'a')
%!error<ClassificationKNN.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationKNN (ones (4,2), ones (4,1)), "holdout", 2)
%!error<ClassificationKNN.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (ClassificationKNN (ones (4,2), ones (4,1)), "leaveout", 1)
%!error<ClassificationKNN.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (ClassificationKNN (ones (4,2), ones (4,1)), "cvpartition", 1)

