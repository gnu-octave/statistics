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
## @multitable @columnfractions 0.02 0.35 0.7
## @headitem @tab @var{Field} @tab @var{Description}
##
## @item @tab @qcode{"obj.X"} @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @tab @qcode{"obj.Y"} @tab Class labels, specified as a logical or
## numeric vector, or cell array of character vectors.  Each value in @var{Y} is
## the observed class label for the corresponding row in @var{X}.
##
## @item @tab @qcode{"obj.NumClasses"} @tab The number of classes in the
## classification problem.
##
## @item @tab @qcode{"obj.ClassNames"} @tab The labels for each class in the
## classification problem. It provides the actual names or identifiers for the
## classes being predicted by the model.
##
## @item @tab @qcode{"obj.NumObservations"} @tab Number of observations used in
## training the model, specified as a positive integer scalar. This number can
## be less than the number of rows in the training data because rows containing
## @qcode{NaN} values are not part of the fit.
##
## @item @tab @qcode{"obj.LayerSizes"} @tab Response variable name, specified
## as a character vector.
##
## @item @tab @qcode{"obj.LayerWeights"} @tab Response variable name, specified
## as a character vector.
##
## @item @tab @qcode{"obj.LayerBiases"} @tab Response variable name, specified
## as a character vector.
##
## @item @tab @qcode{"obj.Activations"} @tab Response variable name, specified
## as a character vector.
##
## @item @tab @qcode{"obj.OutputLayerActivation"} @tab Response variable name, specified
## as a character vector.
##
## @item @tab @qcode{"obj.ModelParameters"} @tab Response variable name, specified
## as a character vector.
##   Weight                        = [];     # Weights of observations used to train this model
##
##    NumObservations                   = [];     # Number of observations in training dataset
##    PredictorNames                    = [];     # Predictor variables names
##    ResponseName                      = [];     # Response variable name
##    RowsUsed                          = [];     # Rows used in fitting
##    NumIterations                     = [];     # Number of iterations taken by optimization
##
## @item @tab @qcode{"obj.ConvergenceInfo"} @tab Response variable name, specified
## as a character vector.
##
## @item @tab @qcode{"obj.TrainingHistory"} @tab Response variable name, specified
## as a character vector.
##
## @item @tab @qcode{"obj.Solver"} @tab Response variable name, specified
## as a character vector.
##
## @end multitable
##
## @seealso{fitcnet}
## @end deftypefn

  properties (Access = public)

    X                     = [];  # Predictor data
    Y                     = [];  # Class labels

    NumClasses            = [];
    ClassNames            = [];  # Names of classes in Y
    NumObservations       = [];  # Number of observations in training dataset

    LayerSizes            = [];  # Size of fully connected layers
    LayerWeights          = [];  # Learned layer weights
    LayerBiases           = [];  # Learned layer biases
    Activations           = [];  # Activation function for fully connected layer
    OutputLayerActivation = [];  # Activation function for final connected layer

    ModelParameters       = [];  # Model parameters
    ConvergenceInfo       = [];  # Convergence Information
    TrainingHistory       = [];  # Training history
    Solver                = [];  # Solver used

  endproperties

##  properties (Access = private)
##
##    gY                    = [];  # Numeric indices
##    gnY                   = [];  # Unique group names
##    glY                   = [];  # Original group labels
##
##  endproperties

  methods (Access = public)

    ## Class object constructor
    function this = ClassificationNeuralNetwork (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationNeuralNetwork: too few input arguments.");
      endif

      ## Get training sample size and number of variables in training data
      nsample = rows (X);                    # Number of samples in X
      ndims_X = columns (X);                 # Number of dimensions in X

      ## Check correspodence between predictors and response
      if (nsample != rows (Y))
        error (strcat (["ClassificationNeuralNetwork: number of rows in X"], ...
                       [" and Y must be equal."]));
      endif

      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);
##      this.gY = gY;
##      this.gnY = gnY;
##      this.glY = glY;

      ## Remove missing values from X and Y
      RowsUsed  = ! logical (sum (isnan ([X, gY]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

##      ## Validate that Y is numeric
##      if (!isnumeric(Y))
##        error ("ClassificationNeuralNetwork: Y must be a numeric array.");
##      endif

      ## Set default values before parsing optional parameters
      LayerSizes              = 10;
      Activations             = 'relu';
      LayerWeightsInitializer = 'glorot';
      LayerBiasesInitializer  = 'zeros';
      InitialStepSize         = [];
      IterationLimit          = 1e3;
      GradientTolerance       = 1e-6;
      LossTolerance           = 1e-6;
      StepTolerance           = 1e-6;
      Weights                 = ones ( this.NumObservations, 1);

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

           case 'layersizes'
            LayerSizes = varargin{2};
            if (! (isnumeric(LayerSizes) && isvector(LayerSizes)
              && all(LayerSizes > 0) && all(mod(LayerSizes, 1) == 0)))
              error (strcat (["ClassificationNeuralNetwork: LayerSizes"], ...
                             [" must be a positive integer vector."]));
            endif

           case 'activations'
            Activations = varargin{2};
            if (!(ischar(Activations)))
              error (strcat (["ClassificationNeuralNetwork: Activations"], ...
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
                             [" LayerWeightsInitializer must be a string."]));
            endif
            if (ischar(LayerWeightsInitializer))
              if (! any (strcmpi (tolower(LayerWeightsInitializer), ...
                                 {"glorot", "he"})))
              error (strcat (["ClassificationNeuralNetwork: unsupported"], ...
                             [" LayerWeightsInitializer function."]));
              endif
            endif

          case 'layerbiasesinitializer'
            LayerBiasesInitializer = varargin{2};
            if (!(ischar(LayerBiasesInitializer)))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" LayerBiasesInitializer must be a string."]));
            endif
            if (ischar(LayerBiasesInitializer))
              if (! any (strcmpi (tolower(LayerBiasesInitializer), {"zeros", ...
                    "ones"})))
              error (strcat (["ClassificationNeuralNetwork: unsupported"], ...
                             [" LayerBiasesInitializer function."]));
              endif
            endif

          case 'initialstepsize'
            InitialStepSize = varargin{2};
            if ( !(isscalar(InitialStepSize) && (InitialStepSize > 0)))
              error (strcat (["ClassificationNeuralNetwork:"], ...
                             [" InitialStepSize must be a positive scalar."]));
            endif

          case 'iterationlimit'
            IterationLimit = varargin{2};
            if (! (isnumeric(IterationLimit) && isscalar(IterationLimit)
              && (IterationLimit > 0) && mod(IterationLimit, 1) == 0))
              error (strcat (["ClassificationNeuralNetwork: IterationLimit"], ...
                             [" must be a positive integer."]));
            endif

          case 'gradienttolerance'
            GradientTolerance = varargin{2};
            if (! (isnumeric(GradientTolerance) && isscalar(GradientTolerance)
              && (GradientTolerance >= 0)))
              error (strcat (["ClassificationNeuralNetwork: GradientTolerance"], ...
                             [" must be a non-negative scalar."]));
            endif

          case 'losstolerance'
            LossTolerance = varargin{2};
            if (! (isnumeric(LossTolerance) && isscalar(LossTolerance)
              && (LossTolerance >= 0)))
              error (strcat (["ClassificationNeuralNetwork: LossTolerance"], ...
                             [" must be a non-negative scalar."]));
            endif

          case 'steptolerance'
            StepTolerance = varargin{2};
            if (! (isnumeric(StepTolerance) && isscalar(StepTolerance)
              && (StepTolerance >= 0)))
              error (strcat (["ClassificationNeuralNetwork: StepTolerance"], ...
                             [" must be a non-negative scalar."]));
            endif

          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric(Weights) && isvector(Weights)
              && (Weights >= 0)))
              error (strcat (["ClassificationNeuralNetwork: Weights must"], ...
                             [" be a non-negative vector."]));
            endif

          otherwise
            error (strcat (["ClassificationNeuralNetwork: invalid"],...
                           [" parameter name in optional pair arguments."]));

        endswitch
        varargin (1:2) = [];
      endwhile

      this.X = X;
      this.Y = Y;
      this.NumObservations = rows (X);
      this.ClassNames = gnY;
      this.NumClasses = numel (this.ClassNames);
      this.LayerSizes = LayerSizes;
      this.Activations = Activations;
      this.OutputLayerActivation = "softmax";
      this.Solver = "LBFGS";

      inputSize = columns(this.X);
      ## Input layer + Hidden Layer + Output layer
      layerSizes = [inputSize, LayerSizes, this.NumClasses];

      for i = 1:length(layerSizes) - 1
        switch LayerWeightsInitializer
          case 'glorot'
            epsilon = sqrt(2 / (layerSizes(i) + layerSizes(i+1)));
          case 'he'
            epsilon = sqrt(2 / layerSizes(i));
        endswitch

        this.LayerWeights{i} = epsilon * randn(layerSizes(i), layerSizes(i+1));

        switch LayerBiasesInitializer
          case 'zeros'
            this.LayerBiases{i} = zeros(1, layerSizes(i+1));
          case 'ones'
            this.LayerBiases{i} = ones(1, layerSizes(i+1));
        endswitch
      endfor


##      disp(forward(this,X));

##      ## Here we will call the training function
##      this.train();

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{label} =} predict (@var{obj}, @var{X})
    ##
    ## Classify new data points into categories using the Neural Network
    ## classification object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{X})} returns the vector of
    ## labels predicted for the corresponding instances in @var{X}, using the
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
    ## @seealso{fitcnet, ClassificationNeuralNetwork}
    ## @end deftypefn

    function predictions = predict (this, X)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("ClassificationNeuralNetwork.predict: too few input arguments.");
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationNeuralNetwork.predict: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat (["ClassificationNeuralNetwork.predict: X must have the same"], ...
                       [" number of features (columns) as in the Neural Network model."]));
      endif

      ## Forward pass to get network output
      output = this.forward(X);

      ## Get the class with the highest probability for each observation
      [~, predicted_labels] = max(output, [], 2);
      disp(predicted_labels);

      ## Convert numeric labels back to original categorical labels
      predictions = this.reverse_grp2idx(predicted_labels);

    endfunction

   endmethods

  ## Helper functions
  methods (Access = private)

    ## Function to reverse grp2idx
    function original_labels = reverse_grp2idx(this, indices)
      ## Preallocate a cell array for the original labels
      original_labels = cell(size(indices));

      ## Map each index back to its corresponding group name
      for i = 1:length(indices)
        original_labels{i} = this.ClassNames{indices(i)};
      endfor
      ##Convert cell array to a column vector
      original_labels = original_labels(:);
    endfunction

    ## Activation function
    function activation_fn = Activation_fn(this, z)
      switch this.Activations
        case 'relu'
          activation_fn = max(0, z);
        case 'tanh'
          activation_fn = tanh(z);
        case 'sigmoid'
          activation_fn = 1 ./ (1 + exp(-z));
        case 'none'
          activation_fn = z;
      endswitch
    endfunction

    ## Forward pass function
    function output = forward(this, X)
      a = X;  # Input layer activations
      # Total layers including output layer
      numLayers = numel(this.LayerSizes) + 1;

      printf("numlayers: ");
      disp(numLayers);

      for i = 1:numLayers - 1
        # Linear combination
        z = a * this.LayerWeights{i} + this.LayerBiases{i};
        if i < numLayers - 1
          # Activation function for hidden layers
          a = this.Activation_fn(z);
        else
          # Softmax activation for output layer
          a = this.softmax(z);
        endif
      endfor

      output = a;  # Output of the network
    endfunction

    ## Softmax activation function
    function softmax = softmax(this, z)
      exp_z = exp(z - max(z, [], 2));  # Stability improvement
      softmax = exp_z ./ sum(exp_z, 2);
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
%!error<ClassificationNeuralNetwork: Y must be a numeric array.> ...
%! ClassificationNeuralNetwork (ones(5,2), ['A';'B';'A';'C';'B'])
%!error<ClassificationNeuralNetwork: LayerSizes must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", -1)
%!error<ClassificationNeuralNetwork: LayerSizes must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", 0.5)
%!error<ClassificationNeuralNetwork: LayerSizes must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", [1,-2])
%!error<ClassificationNeuralNetwork: LayerSizes must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", [10,20,30.5])
%!error<ClassificationNeuralNetwork: Activations must be a string.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "Activations", 123)
%!error<ClassificationNeuralNetwork: unsupported Activation function.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "Activations", "unsupported_type")
%!error<ClassificationNeuralNetwork: LayerWeightsInitializer must be a string.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerWeightsInitializer", 123)
%!error<ClassificationNeuralNetwork: unsupported LayerWeightsInitializer function.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerWeightsInitializer", "unsupported_type")
%!error<ClassificationNeuralNetwork: LayerBiasesInitializer must be a string.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerBiasesInitializer", 123)
%!error<ClassificationNeuralNetwork: unsupported LayerBiasesInitializer function.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerBiasesInitializer", "unsupported_type")
%!error<ClassificationNeuralNetwork: InitialStepSize must be a positive scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "InitialStepSize", -1)
%!error<ClassificationNeuralNetwork: InitialStepSize must be a positive scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "InitialStepSize", 0)
%!error<ClassificationNeuralNetwork: InitialStepSize must be a positive scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "InitialStepSize", [1, 2])
%!error<ClassificationNeuralNetwork: InitialStepSize must be a positive scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "InitialStepSize", "invalid")
%!error<ClassificationNeuralNetwork: IterationLimit must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "IterationLimit", -1)
%!error<ClassificationNeuralNetwork: IterationLimit must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "IterationLimit", 0.5)
%!error<ClassificationNeuralNetwork: IterationLimit must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "IterationLimit", [1,2])
%!error<ClassificationNeuralNetwork: GradientTolerance must be a non-negative scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "GradientTolerance", -1)
%!error<ClassificationNeuralNetwork: GradientTolerance must be a non-negative scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "GradientTolerance", [1,2])
%!error<ClassificationNeuralNetwork: LossTolerance must be a non-negative scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LossTolerance", -1)
%!error<ClassificationNeuralNetwork: LossTolerance must be a non-negative scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LossTolerance", [1,2])
%!error<ClassificationNeuralNetwork: StepTolerance must be a non-negative scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "StepTolerance", -1)
%!error<ClassificationNeuralNetwork: StepTolerance must be a non-negative scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "StepTolerance", [1,2])
%!error<ClassificationNeuralNetwork: Weights must be a non-negative vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "Weights", -1)
%!error<ClassificationNeuralNetwork: Weights must be a non-negative vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "Weights", [1,-2])
%!error<ClassificationNeuralNetwork: invalid parameter name in optional pair argument> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "some", "some")
