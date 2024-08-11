## Copyright (C) 2024 Pallav Purbia <pallavpurbia@gmail.com>
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

classdef ClassificationNeuralNetwork
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{obj} =} ClassificationNeuralNetwork (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} ClassificationNeuralNetwork (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{ClassificationNeuralNetwork} class object containing a Neural
## Network classification model.
##
## @code{@var{obj} = ClassificationNeuralNetwork (@var{X}, @var{Y})} returns a
## ClassificationNeuralNetwork object, with @var{X} as the predictor data and
## @var{Y} containing the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the model.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain any type of
## categorical data. @var{Y} must have same numbers of rows as @var{X}.
## @end itemize
##
## @code{@var{obj} = ClassificationNeuralNetwork (@dots{}, @var{name},
## @var{value})} returns a ClassificationNeuralNetwork object with parameters
## specified by @qcode{Name-Value} pair arguments. Type @code{help fitcnet}
## for more info.
##
## A @qcode{ClassificationNeuralNetwork} object, @var{obj}, stores the labelled
## training data and various parameters for the Neural Network classification
## model, which can be accessed in the following fields:
##
## @multitable @columnfractions 0.32 0.02 0.7
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
## training the model, specified as a positive integer scalar. This number can
## be less than the number of rows in the training data because rows containing
## @qcode{NaN} values are not part of the fit.
##
## @item @qcode{obj.RowsUsed} @tab @tab Rows of the original training data
## used in fitting the ClassificationNeuralNetwork model, specified as a
## numerical vector. If you want to use this vector for indexing the training
## data in @var{X}, you have to convert it to a logical vector, i.e
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
## @item @qcode{obj.ClassNames} @tab @tab Names of the classes in the class
## labels, @var{Y}, used for fitting the ClassificationNeuralNetwork model.
## @qcode{ClassNames} are of the same type as the class labels in @var{Y}.
##
## @item @qcode{obj.LayerSizes} @tab @tab Sizes of the fully connected layers
## in the neural network model, returned as a positive integer vector. The ith
## element of LayerSizes is the number of outputs in the ith fully connected
## layer of the neural network model. LayerSizes does not include the size of
## the final fully connected layer. This layer always has K outputs, where K
## is the number of classes in Y.
##
## @item @qcode{obj.LayerWeights} @tab @tab A cell array containing the learned
## weights for the fully connected layers. Each element of the cell array
## represents the weights for the corresponding fully connected layer.
##
## @item @qcode{obj.LayerBiases} @tab @tab A cell array containing the learned
## biases for the fully connected layers. Each element of the cell array
## corresponds to the biases for the respective fully connected layer,
## including the final layer.
##
## @item @qcode{obj.Activations} @tab @tab  A character vector specifying the
## activation function used for the fully connected layers of the neural network
## model.
##
## @item @qcode{obj.OutputLayerActivation} @tab @tab Activation function for the
## final fully connected layer, returned as 'softmax'.
##
## @item @qcode{obj.ModelParameters} @tab @tab A structure containing the
## parameters used to train the Neural Network classifier model with the
## following fields: @code{LayerSizes}, @code{Activations},
## @code{LayerWeightsInitializer}, @code{LayerBiasesInitializer},
## @code{IterationLimit}, @code{LossTolerance}, and @code{StepTolerance}. Type
## @code{help fitcnet} for more info on their usage and default values.
##
## @item @qcode{obj.ConvergenceInfo} @tab @tab A structure containing the
## Convergence info of the Neural Network classifier model with the following
## fields:
## @multitable @columnfractions 0.05 0.30 0.75
## @headitem @tab @var{Value} @tab @var{Description}
## @item @tab @qcode{"Iterations"} @tab The count of training iterations
## completed during the neural network model's training process.
## @item @tab @qcode{"TrainingLoss"} @tab The loss value recorded for the model
## after training.
## @item @tab @qcode{"Time"} @tab The cumulative time taken for all iterations,
## measured in seconds.
## @end multitable
##
## @item @qcode{obj.Solver} @tab @tab Solver used to train the neural network
## model, returned as 'Gradient Search'.
##
## @item @qcode{obj.ScoreTransform} @tab @tab A function_handle which is used
## for transforming the Neural Network prediction score into a posterior
## probability.  By default, it is @qcode{'none'}, in which case the
## @code{predict} and @code{resubPredict} methods return the prediction scores.
##
## @end multitable
##
## @seealso{fitcnet}
## @end deftypefn

  properties (Access = public)

    X                     = [];  # Predictor data
    Y                     = [];  # Class labels

    NumObservations       = [];  # Number of observations in training dataset
    RowsUsed              = [];  # Rows used in fitting
    Standardize           = [];  # Flag to standardize predictors
    Sigma                 = [];  # Predictor standard deviations
    Mu                    = [];  # Predictor means

    NumPredictors         = [];  # Number of predictors
    PredictorNames        = [];  # Predictor variables names
    ResponseName          = [];  # Response variable name
    ClassNames            = [];  # Names of classes in Y
    Prior                 = [];  # Prior probability for each class
    Cost                  = [];  # Cost of misclassification

    ScoreTransform        = [];  # Transformation for classification scores

    LayerSizes            = [];  # Size of fully connected layers
    LayerWeights          = [];  # Learned layer weights
    LayerBiases           = [];  # Learned layer biases
    Activations           = [];  # Activation function for fully connected layer
    OutputLayerActivation = [];  # Activation function for final connected layer

    ModelParameters       = [];  # Model parameters
    ConvergenceInfo       = [];  # Training history
    Solver                = [];  # Solver used

  endproperties

  properties (Access = private)

    LayerWeightsInitializer = [];
    LayerBiasesInitializer  = [];
    IterationLimit          = [];
    LossTolerance           = [];
    StepTolerance           = [];

  endproperties

  methods (Access = public)

    ## Class object constructor
    function this = ClassificationNeuralNetwork (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationNeuralNetwork: too few input arguments.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error (strcat (["ClassificationNeuralNetwork: number of rows in X"], ...
                       [" and Y must be equal."]));
      endif

      ## Assign original X and Y data to the ClassificationNeuralNetwork object
      this.X = X;
      this.Y = Y;

      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);

      ## Set default values before parsing optional parameters
      Standardize             = false;
      ResponseName            = [];
      PredictorNames          = [];
      ClassNames              = [];
      Prior                   = [];
      Cost                    = [];
      LayerSizes              = 10;
      Activations             = 'relu';
      LayerWeightsInitializer = 'glorot';
      LayerBiasesInitializer  = 'zeros';
      IterationLimit          = 1e3;
      LossTolerance           = 1e-6;
      StepTolerance           = 1e-6;
      ConvergenceInfo         = struct();
      this.ScoreTransform     = 'none';

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "standardize"
            Standardize = varargin{2};
            if (! (Standardize == true || Standardize == false))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'Standardize' must be either true or false."]));
            endif

          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'PredictorNames' must be supplied as a"], ...
                             [" cellstring array."]));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'PredictorNames' must have the same"], ...
                             [" number of columns as X."]));
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'ResponseName' must be a character vector."]));
            endif

          case "classnames"
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames)
                                          || islogical (ClassNames)))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'ClassNames' must be a cellstring,"], ...
                             [" logical or numeric vector."]));
            endif
            ## Check that all class names are available in gnY
            if (iscellstr (ClassNames))
              if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                   ClassNames, "UniformOutput", false))))
                error (strcat (["ClassificationNeuralNetwork: not all"], ...
                               [" 'ClassNames' are present in Y."]));
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, "UniformOutput", false))))
                error (strcat (["ClassificationNeuralNetwork: not all"], ...
                               [" 'ClassNames' are present in Y."]));
              endif
            endif

          case "prior"
            Prior = varargin{2};
            if (! ((isnumeric (Prior) && isvector (Prior)) ||
                  (strcmpi (Prior, "empirical") || strcmpi (Prior, "uniform"))))
              error (strcat (["ClassificationNeuralNetwork: 'Prior' must"], ...
                             [" be either a numeric vector or a character"], ...
                             [" vector."]));
            endif

          case "cost"
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error (strcat (["ClassificationNeuralNetwork: 'Cost' must"], ...
                             [" be a numeric square matrix."]));
            endif

          case 'layersizes'
            LayerSizes = varargin{2};
            if (! (isnumeric(LayerSizes) && isvector(LayerSizes)
              && all(LayerSizes > 0) && all(mod(LayerSizes, 1) == 0)))
              error (strcat (["ClassificationNeuralNetwork: 'LayerSizes'"], ...
                             [" must be a positive integer vector."]));
            endif

          case 'activations'
            Activations = varargin{2};
            if (!(ischar(Activations)))
              error (strcat (["ClassificationNeuralNetwork: 'Activations'"], ...
                             [" must be a string."]));
            endif
            if (ischar(Activations))
              if (! any (strcmpi (tolower(Activations), {"relu", "tanh", ...
                "sigmoid", "none"})))
              error (strcat (["ClassificationNeuralNetwork: unsupported"], ...
                             [" Activation function."]));
              endif
            endif
            Activations = tolower(Activations);

          case 'layerweightsinitializer'
            LayerWeightsInitializer = varargin{2};
            if (!(ischar(LayerWeightsInitializer)))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'LayerWeightsInitializer' must be a string."]));
            endif
            if (ischar (LayerWeightsInitializer))
              if (! any (strcmpi (tolower(LayerWeightsInitializer), ...
                                 {"glorot", "he"})))
              error (strcat (["ClassificationNeuralNetwork: unsupported"], ...
                             [" 'LayerWeightsInitializer' function."]));
              endif
            endif

          case 'layerbiasesinitializer'
            LayerBiasesInitializer = varargin{2};
            if (!(ischar(LayerBiasesInitializer)))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'LayerBiasesInitializer' must be a string."]));
            endif
            if (ischar(LayerBiasesInitializer))
              if (! any (strcmpi (tolower(LayerBiasesInitializer), {"zeros", ...
                    "ones"})))
              error (strcat (["ClassificationNeuralNetwork: unsupported"], ...
                             [" 'LayerBiasesInitializer' function."]));
              endif
            endif

          case 'iterationlimit'
            IterationLimit = varargin{2};
            if (! (isnumeric(IterationLimit) && isscalar(IterationLimit)
              && (IterationLimit > 0) && mod(IterationLimit, 1) == 0))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'IterationLimit' must be a positive"], ...
                             [" integer."]));
            endif

          case 'losstolerance'
            LossTolerance = varargin{2};
            if (! (isnumeric(LossTolerance) && isscalar(LossTolerance)
              && (LossTolerance >= 0)))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'LossTolerance' must be a non-negative"], ...
                             [" scalar."]));
            endif

          case 'steptolerance'
            StepTolerance = varargin{2};
            if (! (isnumeric(StepTolerance) && isscalar(StepTolerance)
              && (StepTolerance >= 0)))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'StepTolerance' must be a non-negative"], ...
                             [" scalar."]));
            endif

          case "scoretransform"
            ScoreTransform = varargin{2};
            stList = {"doublelogit", "invlogit", "ismax", "logit", "none", ...
                      "identity", "sign", "symmetric", "symmetricismax", ...
                      "symmetriclogit"};
            if (! (ischar (ScoreTransform) ||
                   strcmp (class (ScoreTransform), "function_handle")))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" 'ScoreTransform' must be a character"], ...
                             [" vector or a function handle."]));
            endif
            if (! ismember (ScoreTransform, stList))
              error (strcat (["ClassificationNeuralNetwork: unrecognized"], ...
                             [" 'ScoreTransform' function."]));
            endif
            ## Handle ScoreTransform here
            if (is_function_handle (ScoreTransform))
              m = eye (5);
              if (! isequal (size (m), size (ScoreTransform (m))))
                error (strcat (["ClassificationNeuralNetwork: function"], ...
                               [" handle for 'ScoreTransform' must return"], ...
                               [" the same size as its input."]));
              endif
              this.ScoreTransform = ScoreTransform;
            else
              if (strcmpi ("doublelogit", ScoreTransform))
                this.ScoreTransform = @(x) 1 ./ (1 + exp .^ (-2 * x));
              elseif (strcmpi ("invlogit", ScoreTransform))
                this.ScoreTransform = @(x) log (x ./ (1 - x));
              elseif (strcmpi ("ismax", ScoreTransform))
                this.ScoreTransform = eval (sprintf ("@(x) ismax (x)"));
              elseif (strcmpi ("logit", ScoreTransform))
                this.ScoreTransform = @(x) 1 ./ (1 + exp .^ (-x));
              elseif (strcmpi ("identity", ScoreTransform))
                this.ScoreTransform = 'none';
              elseif (strcmpi ("sign", ScoreTransform))
                this.ScoreTransform = @(x) sign (x);
              elseif (strcmpi ("symmetric", ScoreTransform))
                this.ScoreTransform = @(x) 2 * x - 1;
              elseif (strcmpi ("symmetricismax", ScoreTransform))
                this.ScoreTransform = eval (sprintf ("@(x) symmetricismax (x)"));
              elseif (strcmpi ("symmetriclogit", ScoreTransform))
                this.ScoreTransform = @(x) 2 ./ (1 + exp .^ (-x)) - 1;
              endif
            endif

          otherwise
            error (strcat (["ClassificationNeuralNetwork: invalid"],...
                           [" parameter name in optional pair arguments."]));

        endswitch
        varargin (1:2) = [];
      endwhile

      ## Get number of variables in training data
      ndims_X = columns (X);

      ## Assign the number of predictors to the object
      this.NumPredictors = ndims_X;

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

      ## Renew groups in Y
      [gY, gnY, glY] = grp2idx (Y);
      this.ClassNames = gnY;  # Keep the same type as Y

      ## Force Y into numeric
      if (! isnumeric (Y))
        Y = gY;
      endif

      ## Check X contains valid data
      if (! (isnumeric (X) && isfinite (X)))
        error ("ClassificationNeuralNetwork: invalid values in X.");
      endif

      ## Assign the number of observations and their correspoding indices
      ## on the original data, which will be used for training the model,
      ## to the ClassificationNeuralNetwork object
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
          error (strcat (["ClassificationNeuralNetwork: the elements in"], ...
                         [" 'Prior' do not correspond to selected classes"], ...
                         [" in Y."]));
        endif
        this.Prior = Prior ./ sum (Prior);
      endif
      if (isempty (Cost))
        this.Cost = cast (! eye (numel (gnY)), "double");
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (strcat (["ClassificationNeuralNetwork: the number of"], ...
                         [" rows and columns in 'Cost' must correspond to"], ...
                         [" the selected classes in Y."]));
        endif
        this.Cost = Cost;
      endif

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
      this.ResponseName = ResponseName;
      this.LayerSizes = LayerSizes;
      this.Activations = Activations;
      this.OutputLayerActivation = "softmax";
      this.Solver = "Gradient Search";
      this.LayerWeightsInitializer = LayerWeightsInitializer;
      this.LayerBiasesInitializer = LayerBiasesInitializer;
      this.IterationLimit = IterationLimit;
      this.LossTolerance = LossTolerance;
      this.StepTolerance = StepTolerance;

      this = parameter_initializer(this,LayerWeightsInitializer, ...
                                   LayerBiasesInitializer);
      options = optimset('MaxIter', IterationLimit, 'TolFun', LossTolerance, ...
                         'TolX', StepTolerance);
      initialThetaVec = this.vectorize_parameters();


      ## Start timing the training process
      tic;

      [optThetaVec, cost] = fminunc(@(thetaVec) ...
                 costFunction(thetaVec, this, X, Y), initialThetaVec, options);

      ## Store the total time spent in the training history
      ConvergenceInfo.Time = toc;

      ConvergenceInfo.Iterations = IterationLimit;
      ConvergenceInfo.TrainingLoss = cost;

      [this.LayerWeights, this.LayerBiases] = ...
                                          reshape_parameters(this, optThetaVec);

      ## Populate ModelParameters structure
      params = struct();
      paramList = {'LayerSizes', 'Activations', 'LayerWeightsInitializer', ...
                   'LayerBiasesInitializer', 'IterationLimit', ...
                   'LossTolerance', 'StepTolerance'};
      for i = 1:numel (paramList)
        paramName = paramList{i};
        if (isprop (this, paramName))
          params.(paramName) = this.(paramName);
        endif
      endfor
      this.ModelParameters = params;
      this.ConvergenceInfo = ConvergenceInfo;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{labels} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationNeuralNetwork} {[@var{labels}, @var{scores}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the Neural Network
    ## classification object.
    ##
    ## @code{@var{labels} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## predictor data in @code{obj.X} and corresponding labels, @code{obj.Y},
    ## stored in the Neural Network classification model, @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationNeuralNetwork} class object.
    ## @item
    ## @var{X} must be an @math{MxP} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the SVM model in
    ## @var{obj}.
    ## @end itemize
    ##
    ## @code{[@var{labels}, @var{scores}] = predict (@var{obj}, @var{XC}} also
    ## returns @var{scores}, which represent the probability of each label
    ## belonging to a specific class. For each observation in X, the predicted
    ## class label is the one with the highest score among all classes.
    ## Alternatively, @var{scores} can contain the posterior probabilities if
    ## the ScoreTransform has been previously set.
    ##
    ## @seealso{fitcnet, ClassificationNeuralNetwork}
    ## @end deftypefn

    function [labels, scores] = predict(this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("ClassificationNeuralNetwork.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("ClassificationNeuralNetwork.predict: XC is empty.");
      elseif (columns(this.X) != columns(XC))
        error (strcat (["ClassificationNeuralNetwork.predict: XC must"], ...
                       [" have the same number of features as in the"], ...
                       [" Neural Network model."]));
      endif

      ## Standardize (if necessary)
      if (this.Standardize)
        XC = (XC - this.Mu) ./ this.Sigma;
      endif

      ## Forward propagation
      A = XC;
      for i = 1:length(this.LayerSizes)+1
        Z = A * this.LayerWeights{i}' + this.LayerBiases{i}';
        if (i <= length(this.LayerSizes))
          [A, z] = this.Layer_Activation(Z);
        else
          A = this.softmax(Z);
        endif
      endfor

      scores = A;

      # Get the predicted labels (class with highest probability)
      [~, labels] = max(A, [], 2);
      labels = this.ClassNames(labels);

      if (nargout > 1)
        ## Apply ScoreTransform to return probability estimates
        if (! strcmp (this.ScoreTransform, "none"))
          f = this.ScoreTransform;
          if (! strcmp (class (f), "function_handle"))
            error (strcat (["ClassificationNeuralNetwork.predict:"], ...
                           [" 'ScoreTransform' must be a"], ...
                           [" 'function_handle' object."]));
          endif
          scores = f (scores);
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{labels} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationNeuralNetwork} {[@var{labels}, @var{scores}] =} resubPredict (@var{obj})
    ##
    ## Classify the training data using the trained Neural Network
    ## classification object.
    ##
    ## @code{@var{labels} = resubPredict (@var{obj})} returns the vector of
    ## labels predicted for the corresponding instances in the training data,
    ## using the predictor data in @code{obj.X} and corresponding labels,
    ## @code{obj.Y}, stored in the Neural Network classification model,
    ## @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationNeuralNetwork} class object.
    ## @end itemize
    ##
    ## @code{[@var{labels}, @var{scores}] = resubPredict (@var{obj}, @var{XC})}
    ## also returns @var{scores}, which represent the probability of each label
    ## belonging to a specific class. For each observation in X, the predicted
    ## class label is the one with the highest score among all classes.
    ## Alternatively, @var{scores} can contain the posterior probabilities if
    ## the ScoreTransform has been previously set.
    ##
    ## @seealso{fitcnet, ClassificationNeuralNetwork}
    ## @end deftypefn

    function [labels, scores] = resubPredict (this)

      ## Get used rows (if necessary)
      if (sum (this.RowsUsed) != rows (this.X))
        RowsUsed = logical (this.RowsUsed);
        X = this.X(RowsUsed);
        Y = this.Y(RowsUsed);
      else
        X = this.X;
        Y = this.Y;
      endif

      ## Standardize (if necessary)
      if (this.Standardize)
        X = (X - this.Mu) ./ this.Sigma;
      endif

      ## Predict labels from new data
      ## Forward propagation
      A = X;
      for i = 1:length(this.LayerSizes)+1
        Z = A * this.LayerWeights{i}' + this.LayerBiases{i}';
        if (i <= length(this.LayerSizes))
          [A, z] = this.Layer_Activation(Z);
        else
          A = this.softmax(Z);
        endif
      endfor

      scores = A;

      # Get the predicted labels (class with highest probability)
      [~, labels] = max(A, [], 2);
      labels = this.ClassNames(labels);

      if (nargout > 1)
        ## Apply ScoreTransform to return probability estimates
        if (! strcmp (this.ScoreTransform, "none"))
          f = this.ScoreTransform;
          if (! strcmp (class (f), "function_handle"))
            error (strcat (["ClassificationNeuralNetwork.resubPredict:"], ...
                           [" 'ScoreTransform' must be a"], ...
                           [" 'function_handle' object."]));
          endif
          scores = f (scores);
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationNeuralNetwork} {@var{CVMdl} =} crossval (@dots{}, @var{Name}, @var{Value})
    ##
    ## Cross Validate a Neural Network classification object.
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
    ## @seealso{fitcnet, ClassificationNeuralNetwork, cvpartition,
    ## ClassificationPartitionedModel}
    ## @end deftypefn

    function CVMdl = crossval (this, varargin)

      ## Check for sufficient input arguments
      if (nargin < 1)
        error ("ClassificationNeuralNetwork.crossval: too few input arguments.");
      endif

      if (numel (varargin) == 1)
        error (strcat (["ClassificationNeuralNetwork.crossval: Name-Value"], ...
                       [" arguments must be in pairs."]));
      elseif (numel (varargin) > 2)
        error (strcat (["ClassificationNeuralNetwork.crossval: specify"], ...
                       [" only one of the optional Name-Value paired"], ...
                       [" arguments."]));
      endif

      ## Set default values before parsing optional parameters
      numSamples  = size (this.X, 1);
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
              error (strcat (["ClassificationNeuralNetwork.crossval:"], ...
                             [" 'KFold' must be an integer value greater"], ...
                             [" than 1."]));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat (["ClassificationNeuralNetwork.crossval:"], ...
                             [" 'Holdout' must be a numeric value between"], ...
                             [" 0 and 1."]));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat (["ClassificationNeuralNetwork.crossval:"], ...
                             [" 'Leaveout' must be either 'on' or 'off'."]));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (!(isa (CVPartition, 'cvpartition')))
              error (strcat (["ClassificationNeuralNetwork.crossval:"], ...
                             [" 'CVPartition' must be a 'cvpartition'"], ...
                             [" object."]));
            endif

          otherwise
            error (strcat (["ClassificationNeuralNetwork.crossval: invalid"],...
                           [" parameter name in optional paired arguments."]));
          endswitch
        varargin (1:2) = [];
      endwhile

      ## Determine the cross-validation method to use
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (numSamples, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (numSamples, 'LeaveOut');
      else
        partition = cvpartition (numSamples, 'KFold', numFolds);
      endif

      ## Create a cross-validated model object
      CVMdl = ClassificationPartitionedModel (this, partition);

    endfunction

   endmethods

  ## Helper functions
  methods (Access = private)
    ## Activation function
    function [A, z] = Layer_Activation(this,z)
      switch this.Activations
        case 'relu'
          A = max(0, z);
        case 'tanh'
          A = tanh(z);
        case 'sigmoid'
          A = 1 ./ (1 + exp(-z));
        case 'none'
          A = z;
      endswitch
    endfunction

    ## Activation Gradient
    function dz = Activation_Gradient(this, dA, z)
      switch this.Activations
        case 'relu'
          A = max(0, z);
          dz = dA .* (A > 0);
        case 'tanh'
          A = tanh(z);
          dz = dA .* (1 - A.^2);
        case 'sigmoid'
          A = 1 ./ (1 + exp(-z));
          dz = dA .* A .* (1 - A);
        case 'none'
          A = z;
          dz = dA;
      endswitch
    endfunction

    ## Softmax activation function
    function softmax = softmax(this,z)
      exp_z = exp(z - max(z, [], 2));  # Stability improvement
      softmax = exp_z ./ sum(exp_z, 2);
    endfunction

    ## Initialize weights and biases based on the specified initializers
    function this = parameter_initializer(this, LayerWeightsInitializer, ...
                                          LayerBiasesInitializer)
      numLayers = numel(this.LayerSizes);
      inputSize = this.NumPredictors;

      for i = 1:numLayers
        if (i == 1)
          prevLayerSize = inputSize;
        else
          prevLayerSize = this.LayerSizes(i-1);
        endif

        layerSize = this.LayerSizes(i);

        ## Initialize weights
        if (strcmpi(LayerWeightsInitializer, 'glorot'))
          limit = sqrt(6 / (prevLayerSize + layerSize));
          this.LayerWeights{i} = 2 * limit * (rand(layerSize, ...
                                                   prevLayerSize) - 0.5);
        elseif (strcmpi(LayerWeightsInitializer, 'he'))
          limit = sqrt(2 / prevLayerSize);
          this.LayerWeights{i} = limit * randn(layerSize, prevLayerSize);
        endif

        ## Initialize biases
        if (strcmpi(LayerBiasesInitializer, 'zeros'))
          this.LayerBiases{i} = zeros(layerSize, 1);
        elseif (strcmpi(LayerBiasesInitializer, 'ones'))
          this.LayerBiases{i} = ones(layerSize, 1);
        endif
      endfor

      ## Initialize the weights and biases for the output layer
      outputLayerSize = numel(this.ClassNames);
      prevLayerSize = this.LayerSizes(end);

      ## Initialize output layer weights
      if (strcmpi(LayerWeightsInitializer, 'glorot'))
        limit = sqrt(6 / (prevLayerSize + outputLayerSize));
        this.LayerWeights{end+1} = 2 * limit * (rand(outputLayerSize, ...
                                                prevLayerSize) - 0.5);
      elseif (strcmpi(LayerWeightsInitializer, 'he'))
        limit = sqrt(2 / prevLayerSize);
        this.LayerWeights{end+1} = limit * randn(outputLayerSize, ...
                                                 prevLayerSize);
      endif

      ## Initialize output layer biases
      if (strcmpi(LayerBiasesInitializer, 'zeros'))
        this.LayerBiases{end+1} = zeros(outputLayerSize, 1);
      elseif (strcmpi(LayerBiasesInitializer, 'ones'))
        this.LayerBiases{end+1} = ones(outputLayerSize, 1);
      endif

    endfunction

    ## One Hot Vector Encoder
    function one_hot_vector = one_hot_encoder(this, Y)
       one_hot_vector = bsxfun(@eq, Y(:), 1:max(Y));
    endfunction

    ## Cross Entropy Loss function
    function loss = compute_cross_entropy_loss(this, Y_pred, Y_true)

      ## One-hot encode the true labels
      one_hot_Y_true = this.one_hot_encoder(Y_true);

      ## Number of observations
      m = size(Y_true, 1);

      ## Adding a small value to Y_pred to avoid log(0)
      Y_pred = Y_pred + eps;

      ## Compute the cross-entropy loss
      loss = -sum(sum(one_hot_Y_true .* log(Y_pred))) / m;
    endfunction

    ## Vectorize the parameters so that they can be used for fminunc
    function thetaVec = vectorize_parameters(this)
      thetaVec = [];
      for i = 1:numel(this.LayerWeights)
        thetaVec = [thetaVec; this.LayerWeights{i}(:)];
        thetaVec = [thetaVec; this.LayerBiases{i}(:)];
      endfor
    endfunction

    ## Cost Function
    function [J, gradVec] = costFunction(thetaVec, this, X, Y)
      ## Reshape Parameters
      [LayerWeights, LayerBiases] = reshape_parameters(this, thetaVec);

      ## Initialize storage for forward propagation
      Zs = cell(numel(LayerWeights), 1);
      As = cell(numel(LayerWeights), 1);

      ## Forward propagation
      A = X;
      for i = 1:length(this.LayerSizes)+1
        Zs{i} = A * this.LayerWeights{i}' + this.LayerBiases{i}';

        if (i <= length(this.LayerSizes))
          [A, z] = this.Layer_Activation(Zs{i});
        else
          A = this.softmax(Zs{i});
        endif
        As{i} = A;
      endfor

      ## Compute Loss
      J = this.compute_cross_entropy_loss(A, Y);

      ## Backward Propagation
      m = size(X, 1);
      dA = A - this.one_hot_encoder(Y);
      dWs = cell(numel(LayerWeights), 1);
      dBs = cell(numel(LayerWeights), 1);

      for i = numel(LayerWeights):-1:1
        if i == numel(LayerWeights)
          dZ = dA;
        else
          dA = dZ * LayerWeights{i+1};
          dZ = this.Activation_Gradient(dA, Zs{i});
        endif

        if i == 1
          dWs{i} = (dZ' * X) / m;
        else
          dWs{i} = (dZ' * As{i-1}) / m;
        endif
        dBs{i} = sum(dZ, 1)' / m;
      endfor

      ## Vectorize Gradients
      gradVec = [];
      for i = 1:numel(dWs)
        gradVec = [gradVec; dWs{i}(:)];
        gradVec = [gradVec; dBs{i}(:)];
      endfor

    endfunction

    ## Converts thetaVec to LayerWeights and LayerBiases
    function [LayerWeights, LayerBiases] = reshape_parameters(this, thetaVec)
      LayerWeights = cell(numel(this.LayerSizes) + 1, 1);
      LayerBiases = cell(numel(this.LayerSizes) + 1, 1);

      offset = 0;
      inputSize = this.NumPredictors;

      for i = 1:numel(this.LayerSizes)
        layerSize = this.LayerSizes(i);

        if i == 1
          prevLayerSize = inputSize;
        else
          prevLayerSize = this.LayerSizes(i-1);
        endif

        w_size = layerSize * prevLayerSize;
        b_size = layerSize;

        LayerWeights{i} = reshape(thetaVec(offset + 1:offset + w_size), ...
                                  layerSize, prevLayerSize);
        offset = offset + w_size;
        LayerBiases{i} = reshape(thetaVec(offset + 1:offset + b_size), ...
                                 layerSize, 1);
        offset = offset + b_size;
      endfor

      outputLayerSize = numel(this.ClassNames);
      prevLayerSize = this.LayerSizes(end);
      w_size = outputLayerSize * prevLayerSize;
      b_size = outputLayerSize;

      LayerWeights{end} = reshape(thetaVec(offset + 1:offset + w_size), ...
                                  outputLayerSize, prevLayerSize);
      offset = offset + w_size;
      LayerBiases{end} = reshape(thetaVec(offset + 1:offset + b_size), ...
                                 outputLayerSize, 1);
    endfunction

  endmethods
endclassdef

## Test input validation for constructor
%!error<ClassificationNeuralNetwork: too few input arguments.> ...
%! ClassificationNeuralNetwork ()
%!error<ClassificationNeuralNetwork: too few input arguments.> ...
%! ClassificationNeuralNetwork (ones(10,2))
%!error<ClassificationNeuralNetwork: number of rows in X and Y must be equal.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (5,1))
%!error<ClassificationNeuralNetwork: 'Standardize' must be either true or false.> ...
%! ClassificationNeuralNetwork (ones (5,3), ones (5,1), "standardize", "a")
%!error<ClassificationNeuralNetwork: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), "PredictorNames", ["A"])
%!error<ClassificationNeuralNetwork: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), "PredictorNames", "A")
%!error<ClassificationNeuralNetwork: 'PredictorNames' must have the same number of columns as X.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), "PredictorNames", {"A", "B", "C"})
%!error<ClassificationNeuralNetwork: 'ResponseName' must be a character vector.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), "ResponseName", {"Y"})
%!error<ClassificationNeuralNetwork: 'ResponseName' must be a character vector.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), "ResponseName", 1)
%!error<ClassificationNeuralNetwork: 'ClassNames' must be a cellstring, logical or numeric vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "ClassNames", @(x)x)
%!error<ClassificationNeuralNetwork: 'ClassNames' must be a cellstring, logical or numeric vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "ClassNames", ['a'])
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "ClassNames", [1, 2])
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones(5,2), {'a';'b';'a';'a';'b'}, "ClassNames", {'a','c'})
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones(10,2), logical (ones (10,1)), "ClassNames", [true, false])
%!error<ClassificationNeuralNetwork: 'Prior' must be either a numeric vector or a character vector.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), "Prior", {"1", "2"})
%!error<ClassificationNeuralNetwork: 'Cost' must be a numeric square matrix.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), "Cost", [1, 2])
%!error<ClassificationNeuralNetwork: 'Cost' must be a numeric square matrix.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), "Cost", "string")
%!error<ClassificationNeuralNetwork: 'Cost' must be a numeric square matrix.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), "Cost", {eye(2)})
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", -1)
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", 0.5)
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", [1,-2])
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", [10,20,30.5])
%!error<ClassificationNeuralNetwork: 'Activations' must be a string.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "Activations", 123)
%!error<ClassificationNeuralNetwork: unsupported Activation function.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "Activations", "unsupported_type")
%!error<ClassificationNeuralNetwork: 'LayerWeightsInitializer' must be a string.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerWeightsInitializer", 123)
%!error<ClassificationNeuralNetwork: unsupported 'LayerWeightsInitializer' function.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerWeightsInitializer", "unsupported_type")
%!error<ClassificationNeuralNetwork: 'LayerBiasesInitializer' must be a string.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerBiasesInitializer", 123)
%!error<ClassificationNeuralNetwork: unsupported 'LayerBiasesInitializer' function.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerBiasesInitializer", "unsupported_type")
%!error<ClassificationNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "IterationLimit", -1)
%!error<ClassificationNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "IterationLimit", 0.5)
%!error<ClassificationNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "IterationLimit", [1,2])
%!error<ClassificationNeuralNetwork: 'LossTolerance' must be a non-negative scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LossTolerance", -1)
%!error<ClassificationNeuralNetwork: 'LossTolerance' must be a non-negative scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LossTolerance", [1,2])
%!error<ClassificationNeuralNetwork: 'StepTolerance' must be a non-negative scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "StepTolerance", -1)
%!error<ClassificationNeuralNetwork: 'StepTolerance' must be a non-negative scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "StepTolerance", [1,2])
%!error<ClassificationNeuralNetwork: 'ScoreTransform' must be a character vector or a function handle.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "ScoreTransform", [1,2])
%!error<ClassificationNeuralNetwork: unrecognized 'ScoreTransform' function.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "ScoreTransform", "unsupported_type")
%!error<ClassificationNeuralNetwork: invalid parameter name in optional pair arguments.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "some", "some")
%!error<ClassificationNeuralNetwork: invalid values in X.> ...
%! ClassificationNeuralNetwork ([1;2;3;'a';4], ones (5,1))
%!error<ClassificationNeuralNetwork: invalid values in X.> ...
%! ClassificationNeuralNetwork ([1;2;3;Inf;4], ones (5,1))
%!error<ClassificationNeuralNetwork: the elements in 'Prior' do not correspond to selected classes in Y.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), "Prior", [0,1])
%!error<ClassificationNeuralNetwork: the elements in 'Prior' do not correspond to selected classes in Y.> ...
%! ClassificationNeuralNetwork (ones (5,2), [1;1;2;2;3], "ClassNames", [1,2], "Prior", [0,0.4,0.6])
%!error<ClassificationNeuralNetwork: the number of rows and columns in 'Cost' must correspond to the selected classes in Y.> ...
%! ClassificationNeuralNetwork (ones (5,2), [1;1;2;2;3], "ClassNames", [1,2], "Cost", ones (3))

## Test output for predict method
%!shared x, y, x_train, x_test, y_train, y_test, objST
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));

## Test input validation for predict method
%!error<ClassificationNeuralNetwork.predict: too few input arguments.> ...
%! predict (ClassificationNeuralNetwork (ones (4,2), ones (4,1)))
%!error<ClassificationNeuralNetwork.predict: XC is empty.> ...
%! predict (ClassificationNeuralNetwork (ones (4,2), ones (4,1)), [])
%!error<ClassificationNeuralNetwork.predict: XC must have the same number of features> ...
%! predict (ClassificationNeuralNetwork (ones (4,2), ones (4,1)), 1)
%!test
%! objST = fitcnet (x, y);
%! objST.ScoreTransform = "a";
%!error<ClassificationNeuralNetwork.predict: 'ScoreTransform' must be a 'function_handle' object.> ...
%! [labels, scores] = predict (objST, x);

## Test input validation for resubPredict method
%!error<ClassificationNeuralNetwork.resubPredict: 'ScoreTransform' must be a 'function_handle' object.> ...
%! [labels, scores] = resubPredict (objST);

## Test output for crossval method
%!test
%! Mdl = fitcnet(x,y);
%! CVMdl = crossval (Mdl, "KFold", 5);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 5)
%! assert (class (CVMdl.Trained{1}), "ClassificationNeuralNetwork")
%!test
%! obj = fitcnet (x, y);
%! CVMdl = crossval (obj, "HoldOut", 0.2);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "ClassificationNeuralNetwork")
%!test
%! obj = fitcnet (x, y);
%! CVMdl = crossval (obj, "LeaveOut", 'on');
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "ClassificationNeuralNetwork")

## Test input validation for crossval method
%!error<ClassificationNeuralNetwork.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "KFold")
%!error<ClassificationNeuralNetwork.crossval: specify only one of the optional Name-Value paired arguments.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), ...
%! "KFold", 5, "leaveout", 'on')
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "KFold", 'a')
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "KFold", 1)
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "KFold", -1)
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "KFold", 11.5)
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "KFold", [1,2])
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 'a')
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 11.5)
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "Holdout", -1)
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 0)
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 1)
%!error<ClassificationNeuralNetwork.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "Leaveout", 1)
%!error<ClassificationNeuralNetwork.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "CVPartition", 1)
%!error<ClassificationNeuralNetwork.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "CVPartition", 'a')
%!error<ClassificationNeuralNetwork.crossval: invalid parameter name in optional paired arguments> ...
%! crossval (ClassificationNeuralNetwork (ones (40,2),randi([1, 2], 40, 1)), "some", "some")
