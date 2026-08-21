## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
## Copyright (C) 2024 Pallav Purbia <pallavpurbia@gmail.com>
## Copyright (C) 2024-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef ClassificationPartitionedModel
  ## -*- texinfo -*-
  ## @deftp {statistics} ClassificationPartitionedModel
  ##
  ## Cross-validated classification model
  ##
  ## The @code{ClassificationPartitionedModel} class stores cross-validated
  ## classification models trained on different partitions of the data.
  ## It can predict responses for observations not used for training using
  ## the @code{kfoldPredict} method.
  ##
  ## Create a @code{ClassificationPartitionedModel} object by using the
  ## @code{crossval} function.
  ##
  ## @seealso{crossval}
  ## @end deftp

  properties
    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} BinEdges
    ##
    ## Bin edges
    ##
    ## A cell array specifying the bin edges for binned predictors.
    ## This property is read-only.
    ##
    ## @end deftp
    BinEdges                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} CategoricalPredictors
    ##
    ## Indices of categorical predictors
    ##
    ## A vector of positive integers specifying the indices of categorical
    ## predictors.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix containing the unstandardized predictor data.  Each
    ## column of @var{X} represents one predictor (variable), and each row
    ## represents one observation.  This property is read-only.
    ##
    ## @end deftp
    X                            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Y
    ##
    ## Class labels
    ##
    ## Specified as a logical or numeric column vector, or as a character array
    ## or a cell array of character vectors with the same number of rows as the
    ## predictor data.  Each row in @var{Y} is the observed class label for
    ## the corresponding row in @var{X}.  This property is read-only.
    ##
    ## @end deftp
    Y                            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} ClassNames
    ##
    ## Names of classes in the response variable
    ##
    ## An array of unique values of the response variable @var{Y}, which has the
    ## same data types as the data in @var{Y}.  This property is read-only.
    ## @qcode{ClassNames} can have any of the following datatypes:
    ##
    ## @itemize
    ## @item Cell array of character vectors
    ## @item Character array
    ## @item Logical vector
    ## @item Numeric vector
    ## @end itemize
    ##
    ## @end deftp
    ClassNames                   = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Cost
    ##
    ## Cost of Misclassification
    ##
    ## A square matrix specifying the cost of misclassification of a point.
    ## @qcode{Cost(i,j)} is the cost of classifying a point into class @qcode{j}
    ## if its true class is @qcode{i} (that is, the rows correspond to the true
    ## class and the columns correspond to the predicted class).  The order of
    ## the rows and columns in @qcode{Cost} corresponds to the order of the
    ## classes in @qcode{ClassNames}.  The number of rows and columns in
    ## @qcode{Cost} is the number of unique classes in the response.  By
    ## default, @qcode{Cost(i,j) = 1} if @qcode{i != j}, and
    ## @qcode{Cost(i,j) = 0} if @qcode{i = j}.  In other words, the cost is 0
    ## for correct classification and 1 for incorrect classification.
    ##
    ## @end deftp
    Cost                         = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} CrossValidatedModel
    ##
    ## Cross-validated model class
    ##
    ## A character vector specifying the class of the cross-validated model.
    ## This field contains the type of model that was used for the training,
    ## e.g., @qcode{'ClassificationKNN'}.  This property is read-only.
    ##
    ## @end deftp
    CrossValidatedModel          = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} KFold
    ##
    ## Number of cross-validated folds
    ##
    ## A positive integer value specifying the number of cross-validated folds.
    ## This property is read-only.
    ##
    ## @end deftp
    KFold                        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} ModelParameters
    ##
    ## Model parameters
    ##
    ## A structure containing the model parameters used during training.
    ## This includes any model-specific parameters that were configured prior
    ## to training.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters              = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer value specifying the number of observations in the
    ## training dataset used for training the cross-validated model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations              = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Partition
    ##
    ## Partition configuration
    ##
    ## A @code{cvpartition} object specifying the partition configuration used
    ## for cross-validation.  This field stores the cvpartition instance that
    ## describes how the data was split into training and validation sets.
    ## This property is read-only.
    ##
    ## @end deftp
    Partition                    = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which they appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames               = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Prior
    ##
    ## Prior probability for each class
    ##
    ## A numeric vector specifying the prior probabilities for each class.  The
    ## order of the elements in @qcode{Prior} corresponds to the order of the
    ## classes in @qcode{ClassNames}.  This property is read-only.
    ##
    ## @end deftp
    Prior                        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName                 = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} ScoreTransform
    ##
    ## Transformation function for classification scores
    ##
    ## Specified as a function handle for transforming the classification
    ## scores.  This property is read-only.
    ##
    ## @end deftp
    ScoreTransform               = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Standardize
    ##
    ## Standardize predictors flag
    ##
    ## A logical scalar specifying whether to standardize the predictors.
    ## This property is read-only.
    ##
    ## @end deftp
    Standardize                  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Trained
    ##
    ## Models trained on each fold
    ##
    ## A cell array of models trained on each fold.  Each cell contains a model
    ## trained on the minus-one fold of the data (all but one fold used for
    ## training and the remaining fold used for validation).  This property is
    ## read-only.
    ##
    ## @end deftp
    Trained                      = [];
  endproperties

  properties(Access = private, Hidden)
    STname = 'none';
  endproperties

  methods(Access = public)
    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedModel} {@var{this} =} ClassificationPartitionedModel (@var{Mdl}, @var{Partition})
    ##
    ## Create a @code{ClassificationPartitionedModel} class object for
    ## cross-validation of classification models.
    ##
    ## @code{@var{this} = ClassificationPartitionedModel (@var{Mdl},
    ## @var{Partition})} returns a ClassificationPartitionedModel object, with
    ## @var{Mdl} as the trained classification model object and
    ## @var{Partition} as the partitioning object obtained using
    ## @code{cvpartition}
    ## function.
    ##
    ## @seealso{cvpartition}
    ## @end deftypefn
    function this = ClassificationPartitionedModel (Mdl, Partition)

      ## Check input arguments
      if (nargin < 2)
        error ("ClassificationPartitionedModel: too few input arguments.");
      endif

      ## Check for valid Classification object
      validTypes = {'ClassificationDiscriminant', 'ClassificationGAM', ...
                    'ClassificationKNN', 'ClassificationNeuralNetwork', ...
                    'ClassificationSVM'};
      if (! any (strcmp (class (Mdl), validTypes)))
        error ("ClassificationPartitionedModel: unsupported model type.");
      endif

      ## Check for valid cvpartition object
      if (! strcmp (class (Partition), 'cvpartition'))
        error ("ClassificationPartitionedModel: invalid 'cvpartition' object.");
      endif

      ## Set properties.  The rows dropped for missing values are outside the
      ## partition, so they are dropped here too and every index below, the
      ## partition's included, refers to the same set of observations.
      this.X = Mdl.X;
      this.Y = Mdl.Y;
      this.KFold = Partition.NumTestSets;
      this.Trained = cell (this.KFold, 1);
      this.ClassNames = Mdl.ClassNames;
      this.ResponseName = Mdl.ResponseName;
      this.NumObservations = rows (this.X);
      this.PredictorNames = Mdl.PredictorNames;
      this.Partition = Partition;
      this.CrossValidatedModel = class (Mdl);
      this.ScoreTransform = Mdl.ScoreTransform;
      this.STname = Mdl.STname;
      ## Every classifier reports a prior and a cost now, so they are carried
      ## whatever was cross validated; this used to name three of the five.
      this.Prior = Mdl.Prior;
      this.Cost = Mdl.Cost;
      if (ismember (class (Mdl), validTypes(3:5)))
        this.Standardize = Mdl.Standardize;
      endif

      ## Switch Classification object types
      switch (this.CrossValidatedModel)

        case 'ClassificationDiscriminant'
          ## Arguments to pass in fitcdiscr
          args = {};
          ## List of acceptable parameters for fitcdiscr
          DiscrParams = {'PredictorNames', 'ResponseName', 'ClassNames', ...
                         'Cost', 'DiscrimType', 'Gamma'};
          ## Set parameters
          for i = 1:numel (DiscrParams)
            paramName = DiscrParams{i};
            paramValue = Mdl.(paramName);
            if (! isempty (paramValue))
              args = [args, {paramName, paramValue}];
            endif
          endfor
          ## Add 'FillCoeffs' parameter
          if (isempty (Mdl.Coeffs))
            args = [args, {'FillCoeffs', 'off'}];
          endif

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            tmp = fitcdiscr (this.X(idx, :), this.Y(idx,:), args{:});
            this.Trained{k} = compact (tmp);
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct ();
          paramList = {'DiscrimType', 'FillCoeffs', 'Gamma'};
          for i = 1:numel (paramList)
            paramName = paramList{i};
            if (isprop (Mdl, paramName))
              params.(paramName) = Mdl.(paramName);
            endif
          endfor
          this.ModelParameters = params;

        case 'ClassificationGAM'
          ## Arguments to pass in fitcgam
          args = {};
          ## List of acceptable parameters for fitcdiscr
          GAMparams = {'PredictorNames', 'ResponseName', 'ClassNames', ...
                       'Cost', 'Formula', 'Interactions', 'Knots', 'Order', ...
                       'LearningRate', 'NumIterations'};
          ## Set parameters
          for i = 1:numel (GAMparams)
            paramName = GAMparams{i};
            paramValue = Mdl.(paramName);
            if (! isempty (paramValue))
              args = [args, {paramName, paramValue}];
            endif
          endfor

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            tmp = fitcgam (this.X(idx, :), this.Y(idx,:), args{:});
            this.Trained{k} = compact (tmp);
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct ();
          paramList = {'Formula', 'Interactions', 'Knots', 'Order', 'DoF', ...
                       'LearningRate', 'NumIterations'};
          for i = 1:numel (paramList)
            paramName = paramList{i};
            if (isprop (Mdl, paramName))
              params.(paramName) = Mdl.(paramName);
            endif
          endfor
          this.ModelParameters = params;

        case 'ClassificationKNN'
          ## Arguments to pass in fitcknn
          args = {};
          ## List of acceptable parameters for fitcknn
          KNNparams = {'PredictorNames', 'ResponseName', 'ClassNames', ...
                       'Prior', 'Cost', 'ScoreTransform', 'BreakTies', ...
                       'NSMethod', 'BucketSize', 'NumNeighbors', 'Exponent', ...
                       'Scale', 'Cov', 'Distance', 'DistanceWeight', ...
                       'IncludeTies'};
          ## Set parameters
          for i = 1:numel (KNNparams)
            paramName = KNNparams{i};
            if (isprop (Mdl, paramName))
              paramValue = Mdl.(paramName);
              if (! isempty (paramValue))
                args = [args, {paramName, paramValue}];
              endif
            else
              switch (paramName)
                case 'Cov'
                  if (strcmpi (Mdl.Distance, 'mahalanobis') && ...
                      (! isempty (Mdl.DistParameter)))
                    args = [args, {'Cov', Mdl.DistParameter}];
                  endif
                case 'Exponent'
                  if (strcmpi (Mdl.Distance,'minkowski') && ...
                      (! isempty (Mdl.DistParameter)))
                    args = [args, {'Exponent', Mdl.DistParameter}];
                  endif
                case 'Scale'
                  if (strcmpi (Mdl.Distance,'seuclidean') && ...
                      (! isempty (Mdl.DistParameter)))
                    args = [args, {'Scale', Mdl.DistParameter}];
                  endif
              endswitch
            endif
          endfor

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            this.Trained{k} = fitcknn (this.X(idx, :), this.Y(idx,:), args{:});
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct ();
          paramList = {'NumNeighbors', 'Distance', 'DistParameter', ...
                       'NSMethod', 'DistanceWeight', 'Standardize'};
          for i = 1:numel (paramList)
            paramName = paramList{i};
            if (isprop (Mdl, paramName))
              params.(paramName) = Mdl.(paramName);
            endif
          endfor
          this.ModelParameters = params;

        case 'ClassificationNeuralNetwork'
          ## Arguments to pass in fitcnet
          args = {};
          ## List of acceptable parameters for fitcnet
          NNparams = {'PredictorNames', 'ResponseName', 'ClassNames', ...
                      'ScoreTransform', 'Standardize', 'LayerSizes', ...
                      'Activations', 'OutputLayerActivation', ...
                      'LearningRate', 'IterationLimit', 'DisplayInfo'};
          ## Set parameters
          for i = 1:numel (NNparams)
            paramName = NNparams{i};
            paramValue = Mdl.(paramName);
            if (! isempty (paramValue))
              args = [args, {paramName, paramValue}];
            endif
          endfor

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            tmp = fitcnet (this.X(idx, :), this.Y(idx,:), args{:});
            this.Trained{k} = compact (tmp);
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct ();
          paramList = {'LayerSizes', 'Activations', 'OutputLayerActivation', ...
                       'LearningRate', 'IterationLimit', 'Solver'};
          for i = 1:numel (paramList)
            paramName = paramList{i};
            if (isprop (Mdl, paramName))
              params.(paramName) = Mdl.(paramName);
            endif
          endfor
          this.ModelParameters = params;

        case 'ClassificationSVM'
          ## Get ModelParameters structure from ClassificationSVM object
          params = Mdl.ModelParameters;

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            ## Pass all arguments directly to fitcsvm
            tmp = fitcsvm (this.X(idx, :), this.Y(idx,:), ...
                           'Standardize', Mdl.Standardize, ...
                           'PredictorNames', Mdl.PredictorNames, ...
                           'ResponseName', Mdl.ResponseName, ...
                           'ClassNames', Mdl.ClassNames, ...
                           'SVMtype', params.SVMtype, ...
                           'KernelFunction', params.KernelFunction, ...
                           'PolynomialOrder', params.PolynomialOrder, ...
                           'KernelScale', params.KernelScale, ...
                           'KernelOffset', params.KernelOffset, ...
                           'BoxConstraint', params.BoxConstraint, ...
                           'Nu', params.Nu, ...
                           'CacheSize', params.CacheSize, ...
                           'Tolerance', params.Tolerance, ...
                           'Shrinking', params.Shrinking);
            this.Trained{k} = compact (tmp);
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          this.ModelParameters = params;

      endswitch
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedModel} {@var{label} =} kfoldPredict (@var{this})
    ## @deftypefnx {ClassificationPartitionedModel} {[@var{label}, @var{score}, @var{cost}] =} kfoldPredict (@var{this})
    ##
    ## Predict responses for observations not used for training in a
    ## cross-validated classification model.
    ##
    ## @code{@var{[label, Score, Cost]} = kfoldPredict (@var{this})}
    ## returns the predicted class labels, classification scores, and
    ## classification costs for the data used
    ## to train the cross-validated model @var{this}.
    ##
    ## @var{this} is a @code{ClassificationPartitionedModel} object.
    ## The function predicts the response for each observation that was
    ## held out during training in the cross-validation process.
    ##
    ## @multitable @columnfractions 0.28 0.7
    ## @headitem @var{Output} @tab @var{Description}
    ##
    ## @item @qcode{label} @tab Predicted class labels, returned as a
    ## vector or cell array. The type of @var{label} matches the type of
    ## @var{Y} in the original training data. Each element of @var{label}
    ## corresponds to the predicted class
    ## label for the corresponding row in @var{X}.
    ##
    ## @item @qcode{Score} @tab Classification scores, returned as a
    ## numeric matrix. Each row of @var{Score} corresponds to an observation,
    ## and each column corresponds to a class. The value in row @var{i} and
    ## column @var{j} is the
    ## classification score for class @var{j} for observation @var{i}.
    ##
    ## @item @qcode{Cost} @tab Classification costs, returned as a
    ## numeric matrix. Each row of @var{Cost} corresponds to an observation,
    ## and each column corresponds to a class. The value in row @var{i}
    ## and column @var{j} is the classification cost for class @var{j} for
    ## observation @var{i}. This output is optional and only returned if
    ## requested.
    ## @end multitable
    ##
    ## @seealso{ClassificationKNN, ClassificationSVM,
    ## ClassificationPartitionedModel}
    ## @end deftypefn
    function [label, Score, Cost] = kfoldPredict (this)

      ## Input validation
      ## Models whose compact form reports no cost.  CompactClassificationGAM
      ## returns a label and a score only, so asking it for a third output
      ## raised and kfoldPredict could not run on a cross-validated GAM.
      no_cost_models = {'ClassificationGAM', 'ClassificationNeuralNetwork', ...
                        'ClassificationSVM'};
      no_cost = any (strcmp (this.CrossValidatedModel, no_cost_models));
      if (no_cost && nargout > 2)
        error (strcat ("ClassificationPartitionedModel.kfoldPredict:", ...
                       " 'Cost' output is not supported for %s cross", ...
                       " validated models."), this.CrossValidatedModel);
      endif

      ## Initialize the label vector based on the type of Y
      if (iscellstr (this.Y))
        label = cell (this.NumObservations, 1);
      elseif (islogical (this.Y))
        label = false (this.NumObservations, 1);
      elseif (isnumeric (this.Y))
        label = zeros (this.NumObservations, 1);
      elseif (ischar (this.Y))
        label = char (zeros (this.NumObservations, size (this.Y, 2)));
      endif

      ## Initialize the score and cost matrices
      Score = nan (this.NumObservations, numel (this.ClassNames));
      Cost = nan (this.NumObservations, numel (this.ClassNames));

      ## Predict label, score, and cost (if applicable) for each KFold partition
      for k = 1:this.KFold

        ## Get data and trained model for this fold
        testIdx = test (this.Partition, k);
        model = this.Trained{k};

        ## Train
        if (no_cost)
          [predictedLabel, score] = predict (model, this.X(testIdx, :));
        else
          [predictedLabel, score, cost] = predict (model, this.X(testIdx, :));
        endif

        ## Convert cell array of labels to appropriate type (if applicable)
        if (iscell (predictedLabel))
          if (isnumeric (this.Y))
            predictedLabel = cellfun (@str2num, predictedLabel);
          elseif (islogical (this.Y))
            predictedLabel = cellfun (@logical, predictedLabel);
          elseif (iscellstr (this.Y))
            predictedLabel = predictedLabel;
          endif
        endif

        ## Get labels, score, and cost (if applicable)
        label(testIdx) = predictedLabel;
        Score(testIdx, :) = score;
        if (nargout > 2)
          Cost(testIdx, :) = cost;
        endif

      endfor

      ## Handle single fold case (holdout)
      if (this.KFold == 1)
        testIdx = test (this.Partition, 1);
        y = grp2idx (this.Y);
        label(testIdx) = this.Y(find (y == mode (y), 1));
        Score(testIdx, :) = NaN;
        Cost(testIdx, :) = NaN;
        return;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedModel} {@var{L} =} kfoldLoss (@var{obj})
    ## @deftypefnx {ClassificationPartitionedModel} {@var{L} =} kfoldLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the cross-validated classification loss.
    ##
    ## @code{@var{L} = kfoldLoss (@var{obj})} returns the fraction of
    ## observations the folds misclassify, each answered for by the fold's
    ## model that did not see it, which is what @code{kfoldPredict} returns.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationPartitionedModel} object.
    ## @end itemize
    ##
    ## @code{@var{L} = kfoldLoss (@dots{}, @var{name}, @var{value})} accepts
    ## the following @qcode{Name-Value} pairs.
    ##
    ## @multitable @columnfractions 0.24 0.76
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab @qcode{'classiferror'}, the default,
    ## @qcode{'classifcost'}, @qcode{'mincost'}, or a function handle called
    ## as
    ## @code{@var{lossfun} (@var{C}, @var{S}, @var{W}, @var{Cost})}, where
    ## @var{C} is a logical matrix with one true per row marking the true
    ## class, @var{S} the scores, @var{W} the weights and @var{Cost} the
    ## misclassification cost.
    ##
    ## @item @qcode{'Mode'} @tab @qcode{'average'}, the default, which returns
    ## one number over the observations of every fold asked for, or
    ## @qcode{'individual'}, which returns one number per fold.
    ##
    ## @item @qcode{'Folds'} @tab A vector of fold indices to restrict the
    ## loss to.  It defaults to every fold.
    ## @end multitable
    ##
    ## @seealso{ClassificationPartitionedModel, kfoldPredict}
    ## @end deftypefn
    function L = kfoldLoss (this, varargin)

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      ## Defaults, then the optional pairs
      LossFun = 'classiferror';
      Mode    = 'average';
      Folds   = 1:this.KFold;
      ## Only the losses whose value has been measured against R2024a are
      ## offered.  The margin-based ones are not: ours came out 2x MATLAB's
      ## for hinge and 4x for quadratic, the square of the same factor, which
      ## says its margin is not the score difference this would use.  Shipping
      ## a number that close to right and not right is worse than not
      ## shipping it.
      names   = {'classifcost', 'classiferror', 'mincost'};
      while (numel (varargin) > 0)
        if (! (ischar (varargin{1}) && isrow (varargin{1})))
          error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                         " parameter name must be a character vector."));
        endif
        switch (tolower (varargin{1}))

          case 'lossfun'
            LossFun = varargin{2};
            if (! (is_function_handle (LossFun) ||
                   (ischar (LossFun) && isrow (LossFun))))
              error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                             " 'LossFun' must be a character vector or a", ...
                             " function handle."));
            endif
            if (ischar (LossFun))
              LossFun = tolower (LossFun);
              if (! any (strcmp (LossFun, names)))
                error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                               " unsupported 'LossFun' value."));
              endif
            endif

          case 'mode'
            Mode = varargin{2};
            if (! (ischar (Mode) && isrow (Mode) &&
                   any (strcmpi (Mode, {'average', 'individual'}))))
              error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                             " 'Mode' must be either 'average' or", ...
                             " 'individual'."));
            endif

          case 'folds'
            Folds = varargin{2};
            if (! (isnumeric (Folds) && isvector (Folds)
                   && all (Folds == fix (Folds))
                   && all (Folds >= 1) && all (Folds <= this.KFold)))
              error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                             " 'Folds' must be a vector of fold indices", ...
                             " between 1 and KFold."));
            endif

          otherwise
            error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                           " invalid parameter name in optional paired", ...
                           " arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Every observation answered for by the fold that held it out
      [label, Score] = kfoldPredict (this);
      classes = this.ClassNames;
      K = numel (classes);
      n = rows (this.X);

      ## Index of the true and the predicted class, so the cost matrix and
      ## the score matrix can be addressed by the same column.
      true_idx = zeros (n, 1);
      pred_idx = zeros (n, 1);
      for i = 1:n
        t = find (ismember (classes, this.Y(i)));
        p = find (ismember (classes, label(i)));
        if (! isempty (t))
          true_idx(i) = t(1);
        endif
        if (! isempty (p))
          pred_idx(i) = p(1);
        endif
      endfor

      if (strcmpi (Mode, 'individual'))
        L = nan (numel (Folds), 1);
        for i = 1:numel (Folds)
          idx = test (this.Partition, Folds(i));
          L(i) = foldLoss_ (this, idx, Score, true_idx, pred_idx, LossFun, K);
        endfor
      else
        idx = false (n, 1);
        for i = 1:numel (Folds)
          idx = idx | test (this.Partition, Folds(i));
        endfor
        L = foldLoss_ (this, idx, Score, true_idx, pred_idx, LossFun, K);
      endif

    endfunction

    ## Loss over the observations selected by IDX, weighted uniformly and
    ## normalized over that selection, so a subset of folds is an average
    ## rather than a sum.
    function L = foldLoss_ (this, idx, Score, true_idx, pred_idx, LossFun, K)
      if (! any (idx) || ! any (true_idx(idx)))
        L = NaN;
        return;
      endif
      S = Score(idx, :);
      ti = true_idx(idx);
      pi_ = pred_idx(idx);
      m = sum (idx);
      W = ones (m, 1) / m;

      if (is_function_handle (LossFun))
        C = false (m, K);
        C(sub2ind ([m, K], (1:m)', ti)) = true;
        L = LossFun (C, S, W, this.Cost);
        if (! (isnumeric (L) && isscalar (L)))
          error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                         " 'LossFun' must return a numeric scalar."));
        endif
        return;
      endif

      switch (LossFun)
        case 'classiferror'
          L = sum (W .* (pi_ != ti));

        case 'classifcost'
          L = 0;
          for i = 1:m
            L = L + W(i) * this.Cost(ti(i), pi_(i));
          endfor

        case 'mincost'
          L = 0;
          for i = 1:m
            [~, k] = min (S(i,:) * this.Cost);
            L = L + W(i) * this.Cost(ti(i), k);
          endfor

      endswitch
    endfunction

  endmethods

endclassdef


%!demo
%!
%! load fisheriris
%! x = meas;
%! y = species;
%!
%! ## Create a KNN classifier model
%! obj = fitcknn (x, y, 'NumNeighbors', 5, 'Standardize', 1);
%!
%! ## Create a partition for 5-fold cross-validation
%! partition = cvpartition (y, 'KFold', 5);
%!
%! ## Create the ClassificationPartitionedModel object
%! cvModel = crossval (obj, 'cvPartition', partition)

%!demo
%!
%! load fisheriris
%! x = meas;
%! y = species;
%!
%! ## Create a KNN classifier model
%! obj = fitcknn (x, y, 'NumNeighbors', 5, 'Standardize', 1);
%!
%! ## Create the ClassificationPartitionedModel object
%! cvModel = crossval (obj);
%!
%! ## Predict the class labels for the observations not used for training
%! [label, score, cost] = kfoldPredict (cvModel);
%! fprintf ("Cross-validated accuracy = %1.2f%% (%d/%d)\n", ...
%!          sum (strcmp (label, y)) / numel (y) *100, ...
%!          sum (strcmp (label, y)), numel (y))

## Tests
%!test
%! load fisheriris
%! a = fitcdiscr (meas, species, 'gamma', 0.3);
%! cvModel = crossval (a, 'KFold', 5);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (cvModel.NumObservations, 150);
%! assert_equal (numel (cvModel.Trained), 5);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationDiscriminant");
%! assert_equal (cvModel.CrossValidatedModel, "ClassificationDiscriminant");
%! assert_equal (cvModel.KFold, 5);
%!test
%! load fisheriris
%! a = fitcdiscr (meas, species, 'gamma', 0.5, 'fillcoeffs', 'off');
%! cvModel = crossval (a, 'HoldOut', 0.3);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal ({cvModel.X, cvModel.Y}, {meas, species});
%! assert_equal (cvModel.NumObservations, 150);
%! assert_equal (numel (cvModel.Trained), 1);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationDiscriminant");
%! assert_equal (cvModel.CrossValidatedModel, "ClassificationDiscriminant");
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ['a'; 'a'; 'b'; 'b'];
%! a = fitcgam (x, y, 'Interactions', 'all');
%! cvModel = crossval (a, 'KFold', 2);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (numel (cvModel.Trained), 2);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationGAM");
%! assert_equal (cvModel.CrossValidatedModel, "ClassificationGAM");
%! assert_equal (cvModel.KFold, 2);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ['a'; 'a'; 'b'; 'b'];
%! a = fitcgam (x, y);
%! cvModel = crossval (a, 'LeaveOut', 'on');
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal ({cvModel.X, cvModel.Y}, {x, y});
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (numel (cvModel.Trained), 4);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationGAM");
%! assert_equal (cvModel.CrossValidatedModel, "ClassificationGAM");
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ['a'; 'a'; 'b'; 'b'];
%! a = fitcknn (x, y);
%! partition = cvpartition (y, 'KFold', 2);
%! cvModel = ClassificationPartitionedModel (a, partition);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (class (cvModel.Trained{1}), "ClassificationKNN");
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (cvModel.ModelParameters.NumNeighbors, 1);
%! assert_equal (cvModel.ModelParameters.NSMethod, "kdtree");
%! assert_equal (cvModel.ModelParameters.Distance, "euclidean");
%! assert_equal (! cvModel.ModelParameters.Standardize, true);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ['a'; 'a'; 'b'; 'b'];
%! a = fitcknn (x, y, 'NSMethod', 'exhaustive');
%! partition = cvpartition (y, 'HoldOut', 0.2);
%! cvModel = ClassificationPartitionedModel (a, partition);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (class (cvModel.Trained{1}), "ClassificationKNN");
%! assert_equal ({cvModel.X, cvModel.Y}, {x, y});
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (cvModel.ModelParameters.NumNeighbors, 1);
%! assert_equal (cvModel.ModelParameters.NSMethod, "exhaustive");
%! assert_equal (cvModel.ModelParameters.Distance, "euclidean");
%! assert_equal (! cvModel.ModelParameters.Standardize, true);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ['a'; 'a'; 'b'; 'b'];
%! k = 2;
%! a = fitcknn (x, y, 'NumNeighbors' ,k);
%! partition = cvpartition (numel (y), 'LeaveOut');
%! cvModel = ClassificationPartitionedModel (a, partition);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (class (cvModel.Trained{1}), "ClassificationKNN");
%! assert_equal ({cvModel.X, cvModel.Y}, {x, y});
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (cvModel.ModelParameters.NumNeighbors, k);
%! assert_equal (cvModel.ModelParameters.NSMethod, "kdtree");
%! assert_equal (cvModel.ModelParameters.Distance, "euclidean");
%! assert_equal (! cvModel.ModelParameters.Standardize, true);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = {'a'; 'a'; 'b'; 'b'};
%! a = fitcnet (x, y, 'IterationLimit', 50);
%! cvModel = crossval (a, 'KFold', 2);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (numel (cvModel.Trained), 2);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationNeuralNetwork");
%! assert_equal (cvModel.CrossValidatedModel, "ClassificationNeuralNetwork");
%! assert_equal (cvModel.KFold, 2);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = {'a'; 'a'; 'b'; 'b'};
%! a = fitcnet (x, y, 'LayerSizes', [5, 3]);
%! cvModel = crossval (a, 'LeaveOut', 'on');
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal ({cvModel.X, cvModel.Y}, {x, y});
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (numel (cvModel.Trained), 4);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationNeuralNetwork");
%! assert_equal (cvModel.CrossValidatedModel, "ClassificationNeuralNetwork");
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%! SVMModel = fitcsvm (x,y);
%! CVMdl = crossval (SVMModel, 'KFold', 5);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (CVMdl.KFold == 5, true)
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationSVM")
%! assert_equal (CVMdl.CrossValidatedModel, "ClassificationSVM");
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%! obj = fitcsvm (x, y);
%! CVMdl = crossval (obj, 'HoldOut', 0.2);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationSVM")
%! assert_equal (CVMdl.CrossValidatedModel, "ClassificationSVM");
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%! obj = fitcsvm (x, y);
%! CVMdl = crossval (obj, 'LeaveOut', 'on');
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationSVM")
%! assert_equal (CVMdl.CrossValidatedModel, "ClassificationSVM");

## Test input validation for ClassificationPartitionedModel
%!error<ClassificationPartitionedModel: too few input arguments.> ...
%! ClassificationPartitionedModel ()
%!error<ClassificationPartitionedModel: too few input arguments.> ...
%! ClassificationPartitionedModel (ClassificationKNN (ones (4,2), ones (4,1)))
%!error<ClassificationPartitionedModel: unsupported model type.> ...
%! ClassificationPartitionedModel (RegressionGAM (ones (40,2), ...
%! randi ([1, 2], 40, 1)), cvpartition (randi ([1, 2], 40, 1), 'Holdout', 0.3))
%!error<ClassificationPartitionedModel: invalid 'cvpartition' object.> ...
%! ClassificationPartitionedModel (ClassificationKNN (ones (4,2), ...
%! ones (4,1)), 'Holdout')

## Test for kfoldPredict
%!test
%! load fisheriris
%! a = fitcdiscr (meas, species, 'gamma', 0.5, 'fillcoeffs', 'off');
%! cvModel = crossval (a, 'Kfold', 4);
%! [label, score, cost] = kfoldPredict (cvModel);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal ({cvModel.X, cvModel.Y}, {meas, species});
%! assert_equal (cvModel.NumObservations, 150);
%!# assert_equal (label, {"b"; "b"; "a"; "a"});
%!# assert_equal (score, [4.5380e-01, 5.4620e-01; 2.4404e-01, 7.5596e-01; ...
%!#         9.9392e-01, 6.0844e-03; 9.9820e-01, 1.8000e-03], 1e-4);
%!# assert_equal (cost, [5.4620e-01, 4.5380e-01; 7.5596e-01, 2.4404e-01; ...
%!#         6.0844e-03, 9.9392e-01; 1.8000e-03, 9.9820e-01], 1e-4);
%!test
%! x = ones (4, 11);
%! y = {'a'; 'a'; 'b'; 'b'};
%! k = 3;
%! a = fitcknn (x, y, 'NumNeighbors', k);
%! partition = cvpartition (numel (y), 'LeaveOut');
%! cvModel = ClassificationPartitionedModel (a, partition);
%! [label, score, cost] = kfoldPredict (cvModel);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal ({cvModel.X, cvModel.Y}, {x, y});
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (cvModel.ModelParameters.NumNeighbors, k);
%! assert_equal (cvModel.ModelParameters.NSMethod, "exhaustive");
%! assert_equal (cvModel.ModelParameters.Distance, "euclidean");
%! assert_equal (! cvModel.ModelParameters.Standardize, true);
%! assert_equal (label, {'b'; 'b'; 'a'; 'a'});
%! assert_equal (score, [0.3333, 0.6667; 0.3333, 0.6667; 0.6667, 0.3333; ...
%!          0.6667, 0.3333], 1e-4);
%! assert_equal (cost, [0.6667, 0.3333; 0.6667, 0.3333; 0.3333, 0.6667; ...
%!          0.3333, 0.6667], 1e-4);

## Test input validation for kfoldPredict
%!error<ClassificationPartitionedModel.kfoldPredict: 'Cost' output is not supported for ClassificationSVM cross validated models.> ...
%! [label, score, cost] = kfoldPredict (crossval (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1))))
%!error<ClassificationPartitionedModel.kfoldPredict: 'Cost' output is not supported for ClassificationNeuralNetwork cross validated models.> ...
%! [label, score, cost] = kfoldPredict (crossval (ClassificationNeuralNetwork (ones (40,2), randi ([1, 2], 40, 1))))

## The partition, the stored data and NumObservations all describe the same
## set of observations.  A row dropped for a missing value is outside all
## three: before this the partition covered it, NumObservations did not, and
## the folds trained on rows their own fit then discarded.
%!test
%! randn ('seed', 42);
%! lab = double (randn (40, 1) > 0) + 1;
%! X = [randn(40, 2); NaN, 1; 2, NaN];
%! Y = [lab; 1; 2];
%! Mdl = fitcsvm (X, Y);
%! assert_equal (Mdl.NumObservations, 42);
%! assert_equal (Mdl.RowsUsed, []);
%! CVMdl = crossval (Mdl, 'KFold', 4);
%! assert_equal (CVMdl.NumObservations, 42);
%! assert_equal (rows (CVMdl.X), 42);
%! assert_equal (rows (CVMdl.Y), 42);
%! assert_equal (CVMdl.Partition.NumObservations, 42);
%! assert_equal (numel (kfoldPredict (CVMdl)), 42);

## The folds stay stratified, the response being passed to cvpartition
## rather than a bare count.
%!test
%! randn ('seed', 7);
%! X = randn (60, 2);
%! Y = [ones(20, 1); 2 * ones(40, 1)];
%! CVMdl = crossval (fitcsvm (X, Y), 'KFold', 4);
%! for k = 1:4
%!   idx = test (CVMdl.Partition, k);
%!   assert_equal (sum (CVMdl.Y(idx) == 1) >= 4, true);
%!   assert_equal (sum (CVMdl.Y(idx) == 2) >= 8, true);
%! endfor

## kfoldLoss defaults to the misclassification rate of the out-of-fold
## answers, which is the identity R2024a shows, and holds for every model
## type the class accepts.
%!test
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! rand ('seed', 42);
%! cvp = cvpartition (Y, 'KFold', 5);
%! for f = {'fitcsvm', 'fitcnet', 'fitcknn', 'fitcdiscr', 'fitcgam'}
%!   CVMdl = crossval (feval (f{1}, X, Y), 'CVPartition', cvp);
%!   assert_equal (kfoldLoss (CVMdl), ...
%!                 mean (! strcmp (kfoldPredict (CVMdl), CVMdl.Y)), 1e-12);
%! endfor

## Leave-one-out on a two-class fixture agrees with R2024a to the digit for
## every loss offered: it reports 0.02 for all three.
%!test
%! load fisheriris
%! idx = [51:75, 101:125];
%! CVMdl = crossval (fitcdiscr (meas(idx,:), species(idx)), 'Leaveout', 'on');
%! assert_equal (CVMdl.KFold, 50);
%! assert_equal (kfoldLoss (CVMdl), 0.02, 1e-12);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', 'classifcost'), 0.02, 1e-12);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', 'mincost'), 0.02, 1e-12);

## 'individual' is one number per fold, 'Folds' restricts to those named, and
## a handle is called with the class indicator, the scores, the weights and
## the cost, its weights summing to one as MATLAB's do.
%!test
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! rand ('seed', 42);
%! cvp = cvpartition (Y, 'KFold', 5);
%! CVMdl = crossval (fitcsvm (X, Y), 'CVPartition', cvp);
%! L = kfoldLoss (CVMdl, 'Mode', 'individual');
%! assert_equal (size (L), [5, 1]);
%! lab = kfoldPredict (CVMdl);
%! k2 = test (cvp, 2);
%! assert_equal (L(2), mean (! strcmp (lab(k2), CVMdl.Y(k2))), 1e-12);
%! sel = test (cvp, 1) | test (cvp, 3);
%! assert_equal (kfoldLoss (CVMdl, 'Folds', [1, 3]), ...
%!               mean (! strcmp (lab(sel), CVMdl.Y(sel))), 1e-12);
%! assert_equal (numel (kfoldLoss (CVMdl, 'Folds', [2, 4], ...
%!                                 'Mode', 'individual')), 2);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', @(C, S, W, Cost) sum (W)), ...
%!               1, 1e-12);

## Cost reaches the cost-aware losses.
%!test
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! rand ('seed', 42);
%! CVMdl = crossval (fitcsvm (X, Y, 'Cost', [0, 4; 1, 0]), 'KFold', 5);
%! assert_equal (CVMdl.Cost, [0, 4; 1, 0]);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', 'classifcost') >= ...
%!               kfoldLoss (CVMdl, 'LossFun', 'classiferror'), true);

## Test input validation for kfoldLoss
%!shared CVK
%! load fisheriris
%! rand ('seed', 42);
%! CVK = crossval (fitcsvm (meas(51:150,:), species(51:150)), 'KFold', 4);
%!error<ClassificationPartitionedModel.kfoldLoss: Name-Value arguments must be in pairs.> ...
%! kfoldLoss (CVK, 'Mode')
%!error<ClassificationPartitionedModel.kfoldLoss: parameter name must be a character vector.> ...
%! kfoldLoss (CVK, 5, 1)
%!error<ClassificationPartitionedModel.kfoldLoss: 'LossFun' must be a character vector or a function handle.> ...
%! kfoldLoss (CVK, 'LossFun', 5)
%!error<ClassificationPartitionedModel.kfoldLoss: unsupported 'LossFun' value.> ...
%! kfoldLoss (CVK, 'LossFun', 'hinge')
%!error<ClassificationPartitionedModel.kfoldLoss: 'LossFun' must return a numeric scalar.> ...
%! kfoldLoss (CVK, 'LossFun', @(C, S, W, Cost) [1, 2])
%!error<ClassificationPartitionedModel.kfoldLoss: 'Mode' must be either 'average' or 'individual'.> ...
%! kfoldLoss (CVK, 'Mode', 'nope')
%!error<ClassificationPartitionedModel.kfoldLoss: 'Folds' must be a vector of fold indices between 1 and KFold.> ...
%! kfoldLoss (CVK, 'Folds', 0)
%!error<ClassificationPartitionedModel.kfoldLoss: invalid parameter name in optional paired arguments.> ...
%! kfoldLoss (CVK, 'Nope', 1)
