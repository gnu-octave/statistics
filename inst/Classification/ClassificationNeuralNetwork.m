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
## training the model, specified as a positive integer scalar. This number can
## be less than the number of rows in the training data because rows containing
## @qcode{NaN} values are not part of the fit.
##
## @item @qcode{RowsUsed} @tab @tab Rows of the original training data
## used in fitting the ClassificationNeuralNetwork model, specified as a
## numerical vector. If you want to use this vector for indexing the training
## data in @var{X}, you have to convert it to a logical vector, i.e
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
## @item @qcode{ClassNames} @tab @tab Names of the classes in the class
## labels, @var{Y}, used for fitting the ClassificationNeuralNetwork model.
## @qcode{ClassNames} are of the same type as the class labels in @var{Y}.
##
## @item @qcode{LayerSizes} @tab @tab Sizes of the fully connected layers
## in the neural network model, returned as a positive integer vector. The ith
## element of LayerSizes is the number of outputs in the ith fully connected
## layer of the neural network model. LayerSizes does not include the size of
## the final fully connected layer. This layer always has K outputs, where K
## is the number of classes in Y.
##
## @item @qcode{Activations} @tab @tab  A character vector or a cell array of
## character vector specifying the activation functions used in the hidden
## layers of the neural network.
##
## @item @qcode{OutputLayerActivation} @tab @tab  A character vector specifying
## the activation function of the output layer the neural network.
##
## @item @qcode{LearningRate} @tab @tab A positive scalar value defining the
## learning rate used by the gradient descend algorithm during training.
##
## @item @qcode{IterationLimit} @tab @tab A positive scalar value defining
## the number of epochs for training the model.
##
## @item @qcode{DisplayInfo} @tab @tab A boolean flag indicating whether to
## print information during training.
##
## @item @qcode{ModelParameters} @tab @tab A structure containing the
## parameters used to train the Neural Network classifier model containing the
## fields @code{LayerWeights} and @code{Activations} as generated by the
## @code{fcnntrain} function.
##
## @item @qcode{ConvergenceInfo} @tab @tab A structure containing the
## Convergence info of the Neural Network classifier model with the following
## fields:
##
## @multitable @columnfractions 0.05 0.30 0.75
## @headitem @tab @var{Fields} @tab @var{Description}
## @item @tab @qcode{Accuracy} @tab The prediction accuracy at each
## iteration during the neural network model's training process.
## @item @tab @qcode{TrainingLoss} @tab The loss value recorded at each
## iteration during the neural network model's training process.
## @item @tab @qcode{Time} @tab The cumulative time taken for all iterations,
## measured in seconds.
## @end multitable
##
## @item @qcode{Solver} @tab @tab Solver used to train the neural network
## model, returned as 'Gradient Search'.
##
## @item @qcode{ScoreTransform} @tab @tab A function_handle which is used
## for transforming the Neural Network prediction score into a posterior
## probability.  By default, it is @qcode{'none'}, in which case the
## @code{predict} and @code{resubPredict} methods return the prediction scores.
##
## @end multitable
##
## @seealso{fitcnet, fcnntrain, fcnnpredict}
## @end deftypefn

  properties (Access = public)

    X                     = [];  # Predictor data
    Y                     = [];  # Class labels

    NumObservations       = [];  # Number of observations in training dataset
    RowsUsed              = [];  # Rows used in fitting
    NumPredictors         = [];  # Number of predictors
    PredictorNames        = [];  # Predictor variables names
    ResponseName          = [];  # Response variable name
    ClassNames            = [];  # Names of classes in Y

    ScoreTransform        = [];  # Transformation for classification scores

    Standardize           = [];  # Flag to standardize predictors
    Sigma                 = [];  # Predictor standard deviations
    Mu                    = [];  # Predictor means

    LayerSizes            = [];  # Size of fully connected layers
    Activations           = [];  # Activation functions for hidden layers
    OutputLayerActivation = [];  # Activation function for output layer
    LearningRate          = [];  # Learning rate for gradient descend
    IterationLimit        = [];  # Number of training epochs

    ModelParameters       = [];  # Model parameters
    ConvergenceInfo       = [];  # Training history
    DisplayInfo           = [];  # Display information during training
    Solver                = [];  # Solver used

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
      fprintf ("\n  ClassificationNeuralNetwork\n\n");
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
      str = repmat ({"%d"}, 1, numel (this.LayerSizes));
      str = strcat ('[', strjoin (str, ' '), ']');
      str = sprintf (str, this.LayerSizes);
      fprintf ("%+25s: %s\n", 'LayerSizes', str);
      if (iscellstr (this.Activations))
        str = repmat ({"'%s'"}, 1, numel (this.Activations));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.Activations{:});
        fprintf ("%+25s: %s\n", 'Activations', str);
      else # character vector
        fprintf ("%+25s: '%s'\n", 'Activations', this.Activations);
      endif
      fprintf ("%+25s: '%s'\n", 'OutputLayerActivation', ...
               this.OutputLayerActivation);
      fprintf ("%+25s: '%s'\n", 'Solver', this.Solver);
    endfunction

    ## Class specific subscripted reference
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("Invalid () indexing for referencing values", ...
                         " in a ClassificationNeuralNetwork object."));
        case '{}'
          error (strcat ("Invalid {} indexing for referencing values", ...
                         " in a ClassificationNeuralNetwork object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("ClassificationNeuralNetwork.subsref: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("ClassificationNeuralNetwork.subref:", ...
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
        error (strcat ("ClassificationNeuralNetwork.subsasgn:", ...
                       " chained subscripts not allowed."));
      endif
      switch s.type
        case '()'
          error (strcat ("Invalid () indexing for assigning values", ...
                         " to a ClassificationNeuralNetwork object."));
        case '{}'
          error (strcat ("Invalid {} indexing for assigning values", ...
                         " to a ClassificationNeuralNetwork object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("ClassificationNeuralNetwork.subsasgn: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          switch (s.subs)
            case 'ScoreTransform'
              name = "ClassificationNeuralNetwork";
              this.ScoreTransform = parseScoreTransform (val, name);
            otherwise
              error (strcat ("ClassificationNeuralNetwork.subsasgn:", ...
                             " unrecongized or read-only property: '%s'"), ...
                             s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods (Access = public)

    ## Constructor
    function this = ClassificationNeuralNetwork (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationNeuralNetwork: too few input arguments.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error (strcat ("ClassificationNeuralNetwork: number of", ...
                       " rows in X and Y must be equal."));
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
      LayerSizes              = 10;
      Activations             = 'sigmoid';
      OutputLayerActivation   = 'sigmoid';
      LearningRate            = 0.01;
      IterationLimit          = 1000;
      DisplayInfo             = false;
      this.ScoreTransform     = 'none';
      this.Solver = "Gradient Descend";

      ## Supported activation functions
      acList = {"linear", "sigmoid", "relu", "tanh", "softmax", ...
                          "lrelu", "prelu", "elu", "gelu"};
      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "standardize"
            Standardize = varargin{2};
            if (! (Standardize == true || Standardize == false))
              error (strcat ("ClassificationNeuralNetwork:", ...
                             " 'Standardize' must be either true or false."));
            endif

          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("ClassificationNeuralNetwork: 'PredictorNames'", ...
                             " must be supplied as a cellstring array."));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat ("ClassificationNeuralNetwork: 'PredictorNames'", ...
                             " must have the same number of columns as X."));
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("ClassificationNeuralNetwork: 'ResponseName'", ...
                             " must be a character vector."));
            endif

          case "classnames"
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames) ||
                   islogical (ClassNames) || ischar (ClassNames)))
              error (strcat ("ClassificationNeuralNetwork: 'ClassNames'", ...
                             " must be a cell array of character vectors,", ...
                             " a logical vector, a numeric vector,", ...
                             " or a character array."));
            endif
            ## Check that all class names are available in gnY
            if (iscellstr (ClassNames))
              ClassNames = cellstr (ClassNames);
              if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                   ClassNames, "UniformOutput", false))))
                error (strcat ("ClassificationNeuralNetwork: not all", ...
                               " 'ClassNames' are present in Y."));
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, "UniformOutput", false))))
                error (strcat ("ClassificationNeuralNetwork: not all", ...
                               " 'ClassNames' are present in Y."));
              endif
            endif

          case "scoretransform"
            name = "ClassificationNeuralNetwork";
            this.ScoreTransform = parseScoreTransform (varargin{2}, name);

          case 'layersizes'
            LayerSizes = varargin{2};
            if (! (isnumeric(LayerSizes) && isvector(LayerSizes)
              && all(LayerSizes > 0) && all(mod(LayerSizes, 1) == 0)))
              error (strcat ("ClassificationNeuralNetwork: 'LayerSizes'", ...
                             " must be a positive integer vector."));
            endif

          case 'learningrate'
            LearningRate = varargin{2};
            if (! (isnumeric(LearningRate) && isscalar (LearningRate) &&
                   LearningRate > 0))
              error (strcat ("ClassificationNeuralNetwork:", ...
                             " 'LearningRate' must be a positive scalar."));
            endif

          case 'activations'
            Activations = varargin{2};
            if (! (ischar (Activations) || iscellstr (Activations)))
              error (strcat ("ClassificationNeuralNetwork: 'Activations'", ...
                        " must be a character vector or a cellstring vector."));
            endif
            if (ischar (Activations))
              if (! any (strcmpi (Activations, acList)))
                error (strcat ("ClassificationNeuralNetwork: unsupported", ...
                               " 'Activation' function."));
              endif
            else
              if (! all (cell2mat (cellfun (@(x) any (strcmpi (x, acList)),
                                   Activations, "UniformOutput", false))))
                error (strcat ("ClassificationNeuralNetwork: unsupported", ...
                               " 'Activation' functions."));
              endif
            endif
            Activations = tolower (Activations);

          case 'outputlayeractivation'
            OutputLayerActivation = varargin{2};
            if (! (ischar (OutputLayerActivation)))
              error (strcat ("ClassificationNeuralNetwork:", ...
                       " 'OutputLayerActivation' must be a character vector."));
            endif
            if (! any (strcmpi (OutputLayerActivation, acList)))
              error (strcat ("ClassificationNeuralNetwork: unsupported", ...
                             " 'OutputLayerActivation' function."));
            endif
            OutputLayerActivation = tolower (OutputLayerActivation);

          case 'iterationlimit'
            IterationLimit = varargin{2};
            if (! (isnumeric(IterationLimit) && isscalar(IterationLimit)
              && (IterationLimit > 0) && mod(IterationLimit, 1) == 0))
              error (strcat ("ClassificationNeuralNetwork:", ...
                             " 'IterationLimit' must be a positive integer."));
            endif

          case "displayinfo"
            DisplayInfo = varargin{2};
            if (! (DisplayInfo == true || DisplayInfo == false))
              error (strcat ("ClassificationNeuralNetwork: 'DisplayInfo'", ...
                             " must be either true or false."));
            endif

          otherwise
            error (strcat ("ClassificationNeuralNetwork: invalid",...
                           " parameter name in optional pair arguments."));

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
        error ("ClassificationNeuralNetwork: invalid values in X.");
      endif

      ## Assign the number of observations and their correspoding indices
      ## on the original data, which will be used for training the model,
      ## to the ClassificationNeuralNetwork object
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

      ## Store training parameters
      this.LayerSizes = LayerSizes;
      this.Activations = Activations;
      this.OutputLayerActivation = OutputLayerActivation;
      this.LearningRate = LearningRate;
      this.IterationLimit = IterationLimit;
      this.DisplayInfo = DisplayInfo;

      ## Encode activations for fcnntrain (expand if needed)
      nlayers = numel (LayerSizes);
      if (ischar (Activations))
        ActivationCodes = ones (1, nlayers) * activationCode (Activations);
      elseif (nlayers != numel (Activations))
        error (strcat ("ClassificationNeuralNetwork: 'Activations'", ...
                       " vector does not match the number of layers."));
      else
        ActivationCodes = [];
        for i = 1:nlayers
          code = activationCode (Activations{i});
          ActivationCodes = [ActivationCodes, code];
        endfor
      endif
      code = activationCode (OutputLayerActivation);
      ActivationCodes = [ActivationCodes, code];

      ## Start the training process
      NumThreads = nproc ();
      Alpha = 0.01;  # used for ReLU and ELU activation layers
      cnn_timer_ = tic;
      Mdl = fcnntrain (X, gY, LayerSizes, ActivationCodes, NumThreads, ...
                       Alpha, LearningRate, IterationLimit, DisplayInfo);

      ## Store training time, Iterations, and Loss
      ConvergenceInfo.Time = toc (cnn_timer_);
      ConvergenceInfo.Accuracy = Mdl.Accuracy;
      ConvergenceInfo.TrainingLoss = Mdl.Loss;

      ## Remove redundant fields
      Mdl = rmfield (Mdl, "Accuracy");
      Mdl = rmfield (Mdl, "Loss");

      ## Save ModelParameters and ConvergenceInfo
      this.ModelParameters = Mdl;
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
    ## trained neural network classification model in @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationNeuralNetwork} class object.
    ## @item
    ## @var{X} must be an @math{MxP} numeric matrix with the same number of
    ## predictors @math{P} as the corresponding predictors of the trained neural
    ## network model in @var{obj}.
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

    function [labels, scores] = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("ClassificationNeuralNetwork.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("ClassificationNeuralNetwork.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("ClassificationNeuralNetwork.predict: XC must have", ...
                       " the same number of predictors as the trained model."));
      endif

      ## Standardize (if necessary)
      if (this.Standardize)
        XC = (XC - this.Mu) ./ this.Sigma;
      endif

      ## Predict labels from new data
      NumThreads = nproc ();
      [labels, scores] = fcnnpredict (this.ModelParameters, XC, NumThreads);

      # Get class labels
      labels = this.ClassNames(labels);

      if (nargout > 1)
        ## Apply ScoreTransform to return probability estimates
        if (! strcmp (this.ScoreTransform, "none"))
          scores = this.ScoreTransform (scores);
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

      ## Get used rows
      X = this.X(RowsUsed);

      ## Standardize (if necessary)
      if (this.Standardize)
        X = (X - this.Mu) ./ this.Sigma;
      endif

      ## Predict labels from existing data
      NumThreads = nproc ();
      [labels, scores] = fcnnpredict (this.ModelParameters, X, NumThreads);

      # Get class labels
      labels = this.ClassNames(labels);

      if (nargout > 1)
        ## Apply ScoreTransform to return probability estimates
        if (! strcmp (this.ScoreTransform, "none"))
          scores = this.ScoreTransform (scores);
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
        error (strcat ("ClassificationNeuralNetwork.crossval: Name-Value", ...
                       " arguments must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("ClassificationNeuralNetwork.crossval:", ...
                       " specify only one of the optional", ...
                       " Name-Value paired arguments."));
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
              error (strcat ("ClassificationNeuralNetwork.crossval:", ...
                             " 'KFold' must be an integer value", ...
                             " greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("ClassificationNeuralNetwork.crossval:", ...
                             " 'Holdout' must be a numeric value", ...
                             " between 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat ("ClassificationNeuralNetwork.crossval:", ...
                             " 'Leaveout' must be either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (!(isa (CVPartition, 'cvpartition')))
              error (strcat ("ClassificationNeuralNetwork.crossval:", ...
                             " 'CVPartition' must be a 'cvpartition' object."));
            endif

          otherwise
            error (strcat ("ClassificationNeuralNetwork.crossval: invalid",...
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
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{CVMdl} =} compact (@var{obj})
    ##
    ## Create a CompactClassificationNeuralNetwork object.
    ##
    ## @code{@var{CVMdl} = compact (@var{obj})} creates a compact version of the
    ## ClassificationNeuralNetwork object, @var{obj}.
    ##
    ## @seealso{fitcnet, ClassificationNeuralNetwork,
    ## CompactClassificationNeuralNetwork}
    ## @end deftypefn

    function CVMdl = compact (this)
      ## Greate a compact model
      CVMdl = CompactClassificationNeuralNetwork (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationNeuralNetwork object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## ClassificationNeuralNetwork object into an Octave binary file, the name
    ## of which is specified in @var{filename}, along with an extra variable,
    ## which defines the type classification object these variables constitute.
    ## Use @code{loadmodel} in order to load a classification object into
    ## Octave's workspace.
    ##
    ## @seealso{loadmodel, fitcnet, ClassificationNeuralNetwork, cvpartition,
    ## ClassificationPartitionedModel}
    ## @end deftypefn

    function savemodel (this, fname)
      ## Generate variable for class name
      classdef_name = "ClassificationNeuralNetwork";

      ## Create variables from model properties
      X = this.X;
      Y = this.Y;
      NumObservations         = this.NumObservations;
      RowsUsed                = this.RowsUsed;
      NumPredictors           = this.NumPredictors;
      PredictorNames          = this.PredictorNames;
      ResponseName            = this.ResponseName;
      ClassNames              = this.ClassNames;
      ScoreTransform          = this.ScoreTransform;
      Standardize             = this.Standardize;
      Sigma                   = this.Sigma;
      Mu                      = this.Mu;
      LayerSizes              = this.LayerSizes;
      Activations             = this.Activations;
      OutputLayerActivation   = this.OutputLayerActivation;
      LearningRate            = this.LearningRate;
      IterationLimit          = this.IterationLimit;
      ModelParameters         = this.ModelParameters;
      ConvergenceInfo         = this.ConvergenceInfo;
      DislayInfo              = this.DislayInfo;
      Solver                  = this.Solver;

      ## Save classdef name and all model properties as individual variables
      save ("-binary", fname, "classdef_name", "X", "Y", "NumObservations", ...
            "RowsUsed", "NumPredictors", "PredictorNames", "ResponseName", ...
            "ClassNames", "ScoreTransform", "Standardize", "Sigma", "Mu", ...
            "LayerSizes", "Activations", "OutputLayerActivation", ...
            "LearningRate", "IterationLimit", "Solver", "ModelParameters", ...
            "ConvergenceInfo", "DislayInfo");
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationNeuralNetwork object
      mdl = ClassificationNeuralNetwork (1, 1);

      ## Get fieldnames from DATA (including private properties)
      names = fieldnames (data);

      ## Copy data into object
      for i = 1:numel (names)
        ## Check fieldnames in DATA match properties in ClassificationNeuralNetwork
        try
          mdl.(names{i}) = data.(names{i});
        catch
          error (strcat ("ClassificationNeuralNetwork.load_model:", ...
                         " invalid model in '%s'."), filename)
        end_try_catch
      endfor
    endfunction

  endmethods

endclassdef

function numCode = activationCode (strCode)
  switch (strCode)
    case "linear"
      numCode = 0;
    case "sigmoid"
      numCode = 1;
    case "relu"
      numCode = 2;
    case "tanh"
      numCode = 3;
    case "softmax"
      numCode = 4;
    case {"lrelu", "prelu"}
      numCode = 5;
    case "elu"
      numCode = 6;
    case "gelu"
      numCode = 7;
  endswitch
endfunction

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
%!error<ClassificationNeuralNetwork: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "ClassNames", @(x)x)
%!error<ClassificationNeuralNetwork: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "ClassNames", {1})
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones (10,1), "ClassNames", [1, 2])
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones(5,2), ['a';'b';'a';'a';'b'], "ClassNames", ['a';'c'])
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones(5,2), {'a';'b';'a';'a';'b'}, "ClassNames", {'a','c'})
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones(10,2), logical (ones (10,1)), "ClassNames", [true, false])
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", -1)
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", 0.5)
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", [1,-2])
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", [10,20,30.5])
%!error<ClassificationNeuralNetwork: 'LearningRate' must be a positive scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LearningRate", -0.1)
%!error<ClassificationNeuralNetwork: 'LearningRate' must be a positive scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LearningRate", [0.1, 0.01])
%!error<ClassificationNeuralNetwork: 'LearningRate' must be a positive scalar.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LearningRate", "a")
%!error<ClassificationNeuralNetwork: 'Activations' must be a character vector or a cellstring vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "Activations", 123)
%!error<ClassificationNeuralNetwork: unsupported 'Activation' function.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "Activations", "unsupported_type")
%!error<ClassificationNeuralNetwork: unsupported 'Activation' functions.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "LayerSizes", [10, 5], ...
%! "Activations", {"sigmoid", "unsupported_type"})
%!error<ClassificationNeuralNetwork: 'Activations' vector does not match the number of layers.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "Activations", {"sigmoid", "relu", "softmax"})
%!error<ClassificationNeuralNetwork: 'OutputLayerActivation' must be a character vector.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "OutputLayerActivation", 123)
%!error<ClassificationNeuralNetwork: unsupported 'OutputLayerActivation' function.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "OutputLayerActivation", "unsupported_type")
%!error<ClassificationNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "IterationLimit", -1)
%!error<ClassificationNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "IterationLimit", 0.5)
%!error<ClassificationNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones(10,2), ones(10,1), "IterationLimit", [1,2])
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

## Test input validation for subsasgn method
%!shared x, y, objST, Mdl
%! load fisheriris
%! x = meas;
%! y = grp2idx (species);
%! Mdl = fitcnet (x, y, "IterationLimit", 100);
%!error<ClassificationNeuralNetwork: unrecognized 'ScoreTransform' function.> ...
%! Mdl.ScoreTransform = "a";

## Test input validation for predict method
%!error<ClassificationNeuralNetwork.predict: too few input arguments.> ...
%! predict (Mdl)
%!error<ClassificationNeuralNetwork.predict: XC is empty.> ...
%! predict (Mdl, [])
%!error<ClassificationNeuralNetwork.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (Mdl, 1)

## Test output for crossval method
%!test
%! CVMdl = crossval (Mdl, "KFold", 5);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 5)
%! assert (class (CVMdl.Trained{1}), "CompactClassificationNeuralNetwork")
%! assert (CVMdl.CrossValidatedModel, "ClassificationNeuralNetwork")
%!test
%! CVMdl = crossval (Mdl, "HoldOut", 0.2);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "CompactClassificationNeuralNetwork")
%! assert (CVMdl.CrossValidatedModel, "ClassificationNeuralNetwork")

## Test input validation for crossval method
%!error<ClassificationNeuralNetwork.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (Mdl, "KFold")
%!error<ClassificationNeuralNetwork.crossval: specify only one of the optional Name-Value paired arguments.> ...
%! crossval (Mdl, "KFold", 5, "leaveout", 'on')
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (Mdl, "KFold", 'a')
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (Mdl, "KFold", 1)
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (Mdl, "KFold", -1)
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (Mdl, "KFold", 11.5)
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (Mdl, "KFold", [1,2])
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (Mdl, "Holdout", 'a')
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (Mdl, "Holdout", 11.5)
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (Mdl, "Holdout", -1)
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (Mdl, "Holdout", 0)
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (Mdl, "Holdout", 1)
%!error<ClassificationNeuralNetwork.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (Mdl, "Leaveout", 1)
%!error<ClassificationNeuralNetwork.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (Mdl, "CVPartition", 1)
%!error<ClassificationNeuralNetwork.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (Mdl, "CVPartition", 'a')
%!error<ClassificationNeuralNetwork.crossval: invalid parameter name in optional paired arguments> ...
%! crossval (Mdl, "some", "some")
