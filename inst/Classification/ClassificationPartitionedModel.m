## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
## Copyright (C) 2024 Pallav Purbia <pallavpurbia@gmail.com>
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

classdef ClassificationPartitionedModel
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{CVMdl} =} ClassificationPartitionedModel (@var{Mdl}, @var{Partition})
##
## Create a @qcode{ClassificationPartitionedModel} class for cross-validation
## of classification models.
##
## @code{@var{CVMdl} = ClassificationPartitionedModel (@var{Mdl},
## @var{Partition})} returns a ClassificationPartitionedModel object, with
## @var{Mdl} as the trained ClassificationKNN or ClassificationSVM object and
## @var{Partition} as the partitioning object obtained using cvpartition
## function.
##
## A @qcode{ClassificationPartitionedModel} object, @var{CVMdl}, stores the
## classification models trained on cross-validated folds
## and various parameters for the cross-validated model,
## which can be accessed in the following fields:
##
## @multitable @columnfractions 0.32 0.02 0.7
## @headitem @var{Field} @tab @tab @var{Description}
##
## @item @qcode{CVMdl.X} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{CVMdl.Y} @tab @tab Class labels, specified as a logical or
## numeric vector, or cell array of character vectors.  Each value in @var{Y} is
## the observed class label for the corresponding row in @var{X}.
##
## @item @qcode{CVMdl.ClassNames} @tab @tab Names of the classes in the training
## data @var{Y} with duplicates removed, specified as a cell array of character
## vectors.
##
## @item @qcode{CVMdl.Cost} @tab @tab Cost of the misclassification of a point,
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
## @item @qcode{CVMdl.CrossValidatedModel} @tab @tab Class of the
## cross-validated model, specified as a character vector. This field
## contains the type of model that was
## used for the training, e.g., @qcode{"ClassificationKNN"}.
##
## @item @qcode{CVMdl.KFold} @tab @tab Number of cross-validated folds,
## specified as a positive interger scalar. Represents how many folds the
## data was divided into for cross-validation purposes.
##
## @item @qcode{CVMdl.ModelParameters} @tab @tab Model parameters used during
## training, specified as a structure. This includes any model-specific
## parameters that were configured prior to training, such as
## @qcode{NumNeighbors} or @qcode{Distance} in the case of a KNN model.
##
## @item @qcode{CVMdl.NumObservations} @tab @tab Number of observations used in
## training the ClassificationKNN model, specified as a positive integer scalar.
## This number can be less than the number of rows in the training data because
## rows containing @qcode{NaN} values are not part of the fit.
##
## @item @qcode{CVMdl.Partition} @tab @tab Partition configuration used for
## cross-validation, specified as a cvpartition object. This field stores the
## cvpartition instance that describes how the data was split into training and
## validation sets.
##
## @item @qcode{CVMdl.PredictorNames} @tab @tab Predictor variable names,
## specified as a cell array of character vectors.  The variable names are in
## the same order in which they appear in the training data @var{X}.
##
## @item @qcode{CVMdl.Prior} @tab @tab Prior probabilities for each class,
## specified as a numeric vector.  The order of the elements in @qcode{Prior}
## corresponds to the order of the classes in @qcode{ClassNames}.
##
## @item @qcode{CVMdl.ResponseName} @tab @tab Response variable name, specified
## as a character vector.
##
## @item @qcode{CVMdl.Trained} @tab @tab Models trained on each fold,
## specified as a cell array of models. Each cell contains a model trained on
## the minus-one fold of the data (all but one fold used for training and the
## remaining fold used for validation).
##
## @end multitable
##
## @seealso{cvpartition, ClassificationDiscriminant, ClassificationGAM,
## ClassificationKNN, ClassificationNeuralNetwork, ClassificationSVM}
## @end deftypefn

  properties
    BinEdges                     = [];
    CategoricalPredictors        = [];
    X                            = [];
    Y                            = [];
    ClassNames                   = [];
    Cost                         = [];
    CrossValidatedModel          = [];
    KFold                        = [];
    ModelParameters              = [];
    NumObservations              = [];
    Partition                    = [];
    PredictorNames               = [];
    Prior                        = [];
    ResponseName                 = [];
    ScoreTransform               = [];
    Standardize                  = [];
    Trained                      = [];
  endproperties

  methods (Access = public)
    ## Constructor to initialize the partitioned model
    function this = ClassificationPartitionedModel (Mdl, Partition)

      ## Check input arguments
      if (nargin < 2)
        error ("ClassificationPartitionedModel: too few input arguments.");
      endif

      ## Check for valid object types
      validTypes = {'ClassificationDiscriminant', 'ClassificationGAM', ...
                    'ClassificationKNN', 'ClassificationNeuralNetwork', ...
                    'ClassificationSVM'};
      if (! any (strcmp (class (Mdl), validTypes)))
        error ("ClassificationPartitionedModel: unsupported model type.");
      endif

      ## Set properties
      this.X = Mdl.X;
      this.Y = Mdl.Y;
      this.KFold = get (Partition, "NumTestSets");
      this.Trained = cell (this.KFold, 1);
      this.ClassNames = Mdl.ClassNames;
      this.ResponseName = Mdl.ResponseName;
      this.NumObservations = Mdl.NumObservations;
      this.PredictorNames = Mdl.PredictorNames;
      this.Partition = Partition;
      this.CrossValidatedModel = class (Mdl);
      this.Prior = Mdl.Prior;
      this.Cost = Mdl.Cost;
      is_valid = {'ClassificationKNN', 'ClassificationNeuralNetwork', ...
                  'ClassificationSVM'};
      if (any (strcmpi (class (Mdl), is_valid)))
        this.Standardize = Mdl.Standardize;
        this.ScoreTransform = Mdl.ScoreTransform;
      endif

      ## Switch Classification object types
      switch (this.CrossValidatedModel)

        case "ClassificationDiscriminant"
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
            this.Trained{k} = fitcdiscr (this.X(idx, :), this.Y(idx), args{:});
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct();
          paramList = {'DiscrimType', 'FillCoeffs', 'Gamma'};
          for i = 1:numel (paramList)
            paramName = paramList{i};
            if (isprop (Mdl, paramName))
              params.(paramName) = Mdl.(paramName);
            endif
          endfor
          this.ModelParameters = params;

        case "ClassificationGAM"
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
            this.Trained{k} = fitcgam (this.X(idx, :), this.Y(idx), args{:});
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct();
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
            this.Trained{k} = fitcknn (this.X(idx, :), this.Y(idx), args{:});
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct();
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
                      'Prior', 'Cost', 'ScoreTransform', 'Standardize', ...
                      'LayerSizes', 'Activations'};
          ## Set parameters
          for i = 1:numel (NNparams)
            paramName = NNparams{i};
            paramValue = Mdl.(paramName);
            if (! isempty (paramValue))
              args = [args, {paramName, paramValue}];
            endif
          endfor
          NNparams = {'LayerBiasesInitializer', 'LayerWeightsInitializer', ...
                      'IterationLimit', 'LossTolerance', 'StepTolerance'};
          for i = 1:numel (NNparams)
            paramName = NNparams{i};
            paramValue = Mdl.ModelParameters.(paramName);
            if (! isempty (paramValue))
              args = [args, {paramName, paramValue}];
            endif
          endfor

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            this.Trained{k} = fitcnet (this.X(idx, :), this.Y(idx), args{:});
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params.LayerSizes = Mdl.ModelParameters.LayerSizes;
          params.Activations = Mdl.ModelParameters.Activations;
          for i = 1:numel (NNparams)
            paramName = NNparams{i};
            paramValue = Mdl.ModelParameters.(paramName);
            params.(paramName) = Mdl.ModelParameters.(paramName);
          endfor
          this.ModelParameters = params;

        case 'ClassificationSVM'
          ## Get ModelParameters structure from ClassificationKNN object
          params = Mdl.ModelParameters;

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            ## Pass all arguments directly to fitcsvm
            this.Trained{k} = fitcsvm (this.X(idx, :), this.Y(idx), ...
                              'Standardize', Mdl.Standardize, ...
                              'PredictorNames', Mdl.PredictorNames, ...
                              'ResponseName', Mdl.ResponseName, ...
                              'ClassNames', Mdl.ClassNames, ...
                              'Prior', Mdl.Prior, ...
                              'Cost', Mdl.Cost, ...
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
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          this.ModelParameters = params;

      endswitch
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedModel} {@var{label} =} kfoldPredict (@var{CVMdl})
    ## @deftypefnx {ClassificationPartitionedModel} {[@var{label}, @var{score}, @var{cost}] =} kfoldPredict (@var{CVMdl})
    ##
    ## Predict responses for observations not used for training in a
    ## cross-validated classification model.
    ##
    ## @code{@var{[label, Score, Cost]} = kfoldPredict (@var{CVMdl})}
    ## returns the predicted class labels, classification scores, and
    ## classification costs for the data used
    ## to train the cross-validated model @var{CVMdl}.
    ##
    ## @var{CVMdl} is a @code{ClassificationPartitionedModel} object.
    ## The function predicts the response for each observation that was
    ## held out during training in the cross-validation process.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Output} @tab @tab @var{Description}
    ##
    ## @item @qcode{label} @tab @tab Predicted class labels, returned as a
    ## vector or cell array. The type of @var{label} matches the type of
    ## @var{Y} in the original training data. Each element of @var{label}
    ## corresponds to the predicted class
    ## label for the corresponding row in @var{X}.
    ##
    ## @item @qcode{Score} @tab @tab Classification scores, returned as a
    ## numeric matrix. Each row of @var{Score} corresponds to an observation,
    ## and each column corresponds to a class. The value in row @var{i} and
    ## column @var{j} is the
    ## classification score for class @var{j} for observation @var{i}.
    ##
    ## @item @qcode{Cost} @tab @tab Classification costs, returned as a
    ## numeric matrix. Each row of @var{Cost} corresponds to an observation,
    ## and each column corresponds to a class. The value in row @var{i}
    ## and column @var{j} is the classification cost for class @var{j} for
    ## observation @var{i}. This output is
    ## optional and only returned if requested.
    ##
    ## @end multitable
    ##
    ## @seealso{ClassificationKNN, ClassificationSVM,
    ## ClassificationPartitionedModel}
    ## @end deftypefn

    function [label, Score, Cost] = kfoldPredict (this)

      ## Input validation
      no_cost_models = {'ClassificationNeuralNetwork', 'ClassificationSVM'};
      no_cost = any (strcmp (this.CrossValidatedModel, no_cost_models));
      if (no_cost && nargout > 2)
        error (strcat (["ClassificationPartitionedModel.kfoldPredict:"], ...
                       [" 'Cost' output is not supported for"], ...
                       [" ClassificationSVM cross validated models."]));
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

      ## Handle single fold case (holdout)
      if (this.KFold == 1)
        testIdx = test (this.Partition, 1);
        label(testIdx) = mode (this.Y);
        Score(testIdx, :) = NaN;
        Cost(testIdx, :) = NaN;
        return;
      endif

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
%! obj = fitcknn (x, y, "NumNeighbors", 5, "Standardize", 1);
%!
%! ## Create a partition for 5-fold cross-validation
%! partition = cvpartition (y, "KFold", 5);
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
%! obj = fitcknn (x, y, "NumNeighbors", 5, "Standardize", 1);
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
%! a = fitcdiscr (meas, species, "gamma", 0.3);
%! cvModel = crossval (a, "KFold", 5);
%! assert (class (cvModel), "ClassificationPartitionedModel");
%! assert (cvModel.NumObservations, 150);
%! assert (numel (cvModel.Trained), 5);
%! assert (cvModel.CrossValidatedModel, "ClassificationDiscriminant");
%! assert (cvModel.KFold, 5);
%!test
%! load fisheriris
%! a = fitcdiscr (meas, species, "gamma", 0.5, "fillcoeffs", "off");
%! cvModel = crossval (a, "HoldOut", 0.3);
%! assert (class (cvModel), "ClassificationPartitionedModel");
%! assert ({cvModel.X, cvModel.Y}, {meas, species});
%! assert (cvModel.NumObservations, 150);
%! assert (numel (cvModel.Trained), 1);
%! assert (cvModel.CrossValidatedModel, "ClassificationDiscriminant");
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcgam (x, y, "Interactions", "all");
%! cvModel = crossval (a, "KFold", 5);
%! assert (class (cvModel), "ClassificationPartitionedModel");
%! assert (cvModel.NumObservations, 4);
%! assert (numel (cvModel.Trained), 5);
%! assert (cvModel.CrossValidatedModel, "ClassificationGAM");
%! assert (cvModel.KFold, 5);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcgam (x, y);
%! cvModel = crossval (a, "LeaveOut", "on");
%! assert (class (cvModel), "ClassificationPartitionedModel");
%! assert ({cvModel.X, cvModel.Y}, {x, y});
%! assert (cvModel.NumObservations, 4);
%! assert (numel (cvModel.Trained), 4);
%! assert (cvModel.CrossValidatedModel, "ClassificationGAM");
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcknn (x, y);
%! partition = cvpartition (y, "KFold", 5);
%! cvModel = ClassificationPartitionedModel (a, partition);
%! assert (class (cvModel), "ClassificationPartitionedModel");
%! assert (cvModel.NumObservations, 4);
%! assert (cvModel.ModelParameters.NumNeighbors, 1);
%! assert (cvModel.ModelParameters.NSMethod, "kdtree");
%! assert (cvModel.ModelParameters.Distance, "euclidean");
%! assert (! cvModel.ModelParameters.Standardize);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcknn (x, y, "NSMethod", "exhaustive");
%! partition = cvpartition (y, "HoldOut", 0.2);
%! cvModel = ClassificationPartitionedModel (a, partition);
%! assert (class (cvModel), "ClassificationPartitionedModel");
%! assert ({cvModel.X, cvModel.Y}, {x, y});
%! assert (cvModel.NumObservations, 4);
%! assert (cvModel.ModelParameters.NumNeighbors, 1);
%! assert (cvModel.ModelParameters.NSMethod, "exhaustive");
%! assert (cvModel.ModelParameters.Distance, "euclidean");
%! assert (! cvModel.ModelParameters.Standardize);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcnet (x, y);
%! cvModel = crossval (a, "KFold", 5);
%! assert (class (cvModel), "ClassificationPartitionedModel");
%! assert (cvModel.NumObservations, 4);
%! assert (numel (cvModel.Trained), 5);
%! assert (cvModel.CrossValidatedModel, "ClassificationNeuralNetwork");
%! assert (cvModel.KFold, 5);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcnet (x, y, "LayerSizes", [5, 3]);
%! cvModel = crossval (a, "LeaveOut", "on");
%! assert (class (cvModel), "ClassificationPartitionedModel");
%! assert ({cvModel.X, cvModel.Y}, {x, y});
%! assert (cvModel.NumObservations, 4);
%! assert (numel (cvModel.Trained), 4);
%! assert (cvModel.CrossValidatedModel, "ClassificationNeuralNetwork");
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! k = 3;
%! a = fitcknn (x, y, "NumNeighbors" ,k);
%! partition = cvpartition (y, "LeaveOut");
%! cvModel = ClassificationPartitionedModel (a, partition);
%! assert (class (cvModel), "ClassificationPartitionedModel");
%! assert ({cvModel.X, cvModel.Y}, {x, y});
%! assert (cvModel.NumObservations, 4);
%! assert (cvModel.ModelParameters.NumNeighbors, k);
%! assert (cvModel.ModelParameters.NSMethod, "kdtree");
%! assert (cvModel.ModelParameters.Distance, "euclidean");
%! assert (! cvModel.ModelParameters.Standardize);
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%! SVMModel = fitcsvm (x,y);
%! CVMdl = crossval (SVMModel, "KFold", 5);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 5)
%! assert (class (CVMdl.Trained{1}), "ClassificationSVM")
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%! obj = fitcsvm (x, y);
%! CVMdl = crossval (obj, "HoldOut", 0.2);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "ClassificationSVM")
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%! obj = fitcsvm (x, y);
%! CVMdl = crossval (obj, "LeaveOut", 'on');
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "ClassificationSVM")

## Test input validation for ClassificationPartitionedModel
%!error<ClassificationPartitionedModel: too few input arguments.> ...
%! ClassificationPartitionedModel ()
%!error<ClassificationPartitionedModel: too few input arguments.> ...
%! ClassificationPartitionedModel (ClassificationKNN (ones (4,2), ones (4,1)))
%!error<ClassificationPartitionedModel: unsupported model type.> ...
%! ClassificationPartitionedModel (RegressionGAM (ones (40,2), ...
%! randi ([1, 2], 40, 1)), cvpartition (randi ([1, 2], 40, 1), 'Holdout', 0.3))

## Test for kfoldPredict
%!test
%! load fisheriris
%! a = fitcdiscr (meas, species, "gamma", 0.5, "fillcoeffs", "off");
%! cvModel = crossval (a, "Kfold", 4);
%! [label, score, cost] = kfoldPredict (cvModel);
%! assert (class(cvModel), "ClassificationPartitionedModel");
%! assert ({cvModel.X, cvModel.Y}, {meas, species});
%! assert (cvModel.NumObservations, 150);
%!# assert (label, {"b"; "b"; "a"; "a"});
%!# assert (score, [4.5380e-01, 5.4620e-01; 2.4404e-01, 7.5596e-01; ...
%!#         9.9392e-01, 6.0844e-03; 9.9820e-01, 1.8000e-03], 1e-4);
%!# assert (cost, [5.4620e-01, 4.5380e-01; 7.5596e-01, 2.4404e-01; ...
%!#         6.0844e-03, 9.9392e-01; 1.8000e-03, 9.9820e-01], 1e-4);
%!test
%! x = ones(4, 11);
%! y = {"a"; "a"; "b"; "b"};
%! k = 3;
%! a = fitcknn (x, y, "NumNeighbors", k);
%! partition = cvpartition (y, "LeaveOut");
%! cvModel = ClassificationPartitionedModel (a, partition);
%! [label, score, cost] = kfoldPredict (cvModel);
%! assert (class(cvModel), "ClassificationPartitionedModel");
%! assert ({cvModel.X, cvModel.Y}, {x, y});
%! assert (cvModel.NumObservations, 4);
%! assert (cvModel.ModelParameters.NumNeighbors, k);
%! assert (cvModel.ModelParameters.NSMethod, "exhaustive");
%! assert (cvModel.ModelParameters.Distance, "euclidean");
%! assert (! cvModel.ModelParameters.Standardize);
%! assert (label, {"b"; "b"; "a"; "a"});
%! assert (score, [0.3333, 0.6667; 0.3333, 0.6667; 0.6667, 0.3333; ...
%!          0.6667, 0.3333], 1e-4);
%! assert (cost, [0.6667, 0.3333; 0.6667, 0.3333; 0.3333, 0.6667; ...
%!          0.3333, 0.6667], 1e-4);

## Test input validation for kfoldPredict
%!error<ClassificationPartitionedModel.kfoldPredict: 'Cost' output is not supported for ClassificationSVM cross validated models.> ...
%! [label, score, cost] = kfoldPredict (crossval (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1))))
