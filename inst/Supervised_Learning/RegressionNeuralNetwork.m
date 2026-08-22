## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{obj} =} RegressionNeuralNetwork (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} RegressionNeuralNetwork (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{RegressionNeuralNetwork} object containing a neural network
## regression model.
##
## @code{@var{obj} = RegressionNeuralNetwork (@var{X}, @var{Y})} returns a
## neural network regression model, @var{obj}, with @var{X} being the predictor
## data and @var{Y} the continuous response of the observations in @var{X}.
##
## @itemize
## @item
## @var{X} must be an @math{NxP} numeric matrix of predictor data, where rows
## correspond to observations and columns to features.
## @item
## @var{Y} must be an @math{Nx1} numeric vector holding the response of the
## corresponding predictor data in @var{X}.  @var{Y} must have the same number
## of rows as @var{X}.
## @end itemize
##
## The network is trained against the mean squared error, and its output layer
## applies the identity, so a prediction is an unrestricted real number rather
## than a score over classes.  This is the only difference in the engine
## between this class and @code{ClassificationNeuralNetwork}; everything else,
## the layer sizes, the activations, the learning rate and the initialisation,
## behaves identically.
##
## @code{@var{obj} = RegressionNeuralNetwork (@dots{}, @var{name},
## @var{value})} returns a model with additional options specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.32 0.68
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Standardize'} @tab A logical scalar specifying whether the
## predictor data should be centred and scaled before training.  The same
## transformation is applied by @code{predict}.  The default is @qcode{false}.
##
## @item @qcode{'PredictorNames'} @tab A cell array of character vectors
## naming the predictors, in the order they appear in @var{X}.
##
## @item @qcode{'ResponseName'} @tab A character vector naming the response.
## The default is @qcode{'Y'}.
##
## @item @qcode{'ResponseTransform'} @tab A character vector naming one of the
## supported transformations, or a function handle, applied to the predicted
## response by @code{predict} and @code{resubPredict}.  The default is
## @qcode{'none'}.
##
## @item @qcode{'LayerSizes'} @tab A positive integer vector specifying the
## number of units in each fully connected hidden layer.  The default is 10,
## one hidden layer of ten units.
##
## @item @qcode{'Activations'} @tab A character vector or cell array of
## character vectors specifying the activation of the hidden layers.  The
## supported functions are @qcode{'linear'}, @qcode{'sigmoid'},
## @qcode{'relu'}, @qcode{'tanh'}, @qcode{'lrelu'}, @qcode{'prelu'},
## @qcode{'elu'}, @qcode{'gelu'} and @qcode{'none'}.  The default is
## @qcode{'relu'}.
##
## @item @qcode{'OutputLayerActivation'} @tab A character vector specifying
## the activation of the output layer.  The default is @qcode{'none'}, the
## identity, which is what a regression output calls for.  The supported
## values are the same as for @qcode{'Activations'}.
##
## @item @qcode{'LearningRate'} @tab A positive scalar specifying the learning
## rate for gradient descent.  The default is 0.003.  A larger rate can drive
## every unit of a hidden layer negative, after which a rectifier passes no
## gradient and the network stops training.
##
## @item @qcode{'IterationLimit'} @tab A positive integer specifying the
## maximum number of training iterations.  The default is 1000.
##
## @item @qcode{'DisplayInfo'} @tab A logical scalar specifying whether to
## print information during training.  The default is @qcode{false}.
## @end multitable
##
## The supported values for @qcode{'ResponseTransform'} are:
##
## @multitable @columnfractions 0.3 0.7
## @headitem @var{Value} @tab @var{Description}
## @item @qcode{'none'} @tab @math{x} (no transformation)
## @item @qcode{'identity'} @tab @math{x} (no transformation)
## @item @qcode{'exp'} @tab @math{exp (x)}
## @item @qcode{'log'} @tab @math{log (x)}
## @end multitable
##
## @seealso{fitrnet, ClassificationNeuralNetwork, fcnntrain, fcnnpredict}
## @end deftypefn

classdef RegressionNeuralNetwork

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} X
    ##
    ## Predictor data
    ##
    ## An @math{NxP} numeric matrix, as it was supplied to the constructor,
    ## before any rows carrying missing values were dropped.  This property is
    ## read-only.
    ##
    ## @end deftp
    X                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} Y
    ##
    ## Response data
    ##
    ## An @math{Nx1} numeric vector, as it was supplied to the constructor.
    ## This property is read-only.
    ##
    ## @end deftp
    Y                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} NumObservations
    ##
    ## Number of observations used to train the model
    ##
    ## A positive integer scalar, counting only the rows that survived the
    ## removal of missing values.  This property is read-only.
    ##
    ## @end deftp
    NumObservations       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} RowsUsed
    ##
    ## Rows used for fitting
    ##
    ## A logical column vector with the same length as the observations in the
    ## original predictor data @var{X}, true for each row that was used for
    ## fitting the RegressionNeuralNetwork model.  It is empty, @qcode{[]},
    ## when every observation was used, so a non-empty value means that rows
    ## holding missing values were dropped.  This property is read-only.
    ##
    ## @end deftp
    RowsUsed              = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer scalar, the number of columns of @code{X}.  This
    ## property is read-only.
    ##
    ## @end deftp
    NumPredictors         = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors, one per column of @code{X}.  This
    ## property is read-only.
    ##
    ## @end deftp
    PredictorNames        = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    ResponseName          = [];


    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} Sigma
    ##
    ## Standard deviation of the predictors
    ##
    ## A row vector with one entry per predictor, used for standardization.
    ## Empty when the predictor data were not standardized.  This property is
    ## read-only.
    ##
    ## @end deftp
    Sigma                 = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} Mu
    ##
    ## Mean of the predictors
    ##
    ## A row vector with one entry per predictor, used for standardization.
    ## Empty when the predictor data were not standardized.  This property is
    ## read-only.
    ##
    ## @end deftp
    Mu                    = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} LayerSizes
    ##
    ## Sizes of the fully connected hidden layers
    ##
    ## A row vector of positive integers, one per hidden layer.  It does not
    ## include the output layer, whose width is the number of responses.  This
    ## property is read-only.
    ##
    ## @end deftp
    LayerSizes            = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} Activations
    ##
    ## Activation functions of the hidden layers
    ##
    ## A character vector, applying to every hidden layer, or a cell array of
    ## character vectors with one entry per hidden layer.  This property is
    ## read-only.
    ##
    ## @end deftp
    Activations           = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} OutputLayerActivation
    ##
    ## Activation function of the output layer
    ##
    ## A character vector.  The default, @qcode{'none'}, applies the identity,
    ## so a prediction is an unrestricted real number.  This property is
    ## read-only.
    ##
    ## @end deftp
    OutputLayerActivation = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} LearningRate
    ##
    ## Learning rate for gradient descent
    ##
    ## A positive scalar value defining the learning rate used by the gradient
    ## descent algorithm during training.  This property is read-only.
    ##
    ## @end deftp
    LearningRate          = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} IterationLimit
    ##
    ## Maximum number of training iterations
    ##
    ## A positive integer scalar.  This property is read-only.
    ##
    ## @end deftp
    IterationLimit        = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} ModelParameters
    ##
    ## Parameters of the trained network
    ##
    ## A structure as returned by @code{fcnntrain}, holding the layer weights
    ## and the activation codes, and consumed by @code{fcnnpredict}.  This
    ## property is read-only.
    ##
    ## @end deftp
    ModelParameters       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} ConvergenceInfo
    ##
    ## Information recorded during training
    ##
    ## A structure with the fields @code{Time}, the seconds training took, and
    ## @code{TrainingLoss}, the mean squared error of the network at the end of
    ## each iteration.  This property is read-only.
    ##
    ## @end deftp
    ConvergenceInfo       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} DisplayInfo
    ##
    ## Whether training printed its progress
    ##
    ## A logical scalar.  This property is read-only.
    ##
    ## @end deftp
    DisplayInfo           = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} Solver
    ##
    ## Solver used to train the network
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    Solver                = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} LayerWeights
    ##
    ## Weights the network learned
    ##
    ## A cell array with one entry per layer, the output layer included.
    ## @code{LayerWeights@{i@}} has one row per unit of layer @math{i} and one
    ## column per input to that layer.  This property is read-only.
    ##
    ## @end deftp
    LayerWeights          = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} LayerBiases
    ##
    ## Biases the network learned
    ##
    ## A cell array with one entry per layer, the output layer included.
    ## @code{LayerBiases@{i@}} is a column with one entry per unit of layer
    ## @math{i}.  This property is read-only.
    ##
    ## @end deftp
    LayerBiases           = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} TrainingHistory
    ##
    ## Iteration by iteration record of the fit
    ##
    ## A @code{table} with the variables @code{Iteration} and
    ## @code{TrainingLoss}, one row per training iteration.  This property is
    ## read-only.
    ##
    ## @end deftp
    TrainingHistory       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector with one entry per training observation.  It
    ## defaults to a uniform weight for every observation.  This property is
    ## read-only.
    ##
    ## @end deftp
    W                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector of column indices into @code{X} naming the predictors
    ## treated as categorical, and empty when none is.  This property is
    ## read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the model expanded them
    ##
    ## A cell array of character vectors.  It matches @code{PredictorNames}
    ## unless a categorical predictor was expanded into indicator variables.
    ## This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};
    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## A cell array with one entry per predictor, holding that predictor's bin
    ## edges where the learner discretized it before fitting.  It is empty here
    ## and stays empty: this learner fits the predictors as they are, and
    ## MATLAB's reports an empty cell for it as well.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    BinEdges        = {};

  endproperties

  ## Properties a user may set after the model is built.  Each one is
  ## validated by its set method below.
  properties (GetAccess = public, SetAccess = public)
    ## -*- texinfo -*-
    ## @deftp {RegressionNeuralNetwork} {property} ResponseTransform
    ##
    ## Transformation applied to the predicted response
    ##
    ## A function handle, applied by @code{predict} and @code{resubPredict} to
    ## the network's output.  It defaults to the identity and may be set after
    ## construction, either to a handle or to the name of a supported
    ## transformation.
    ##
    ## @end deftp
    ResponseTransform     = 'none';
  endproperties

  ## Readable by the counterpart class, which copies it, and kept out of
  ## the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    RTfun = @(y) y;
  endproperties

  ## Set methods for the properties a user may assign.
  methods (Hidden)

    function this = set.ResponseTransform (this, val)
        name = 'RegressionNeuralNetwork';
        [this.RTfun, this.ResponseTransform] = ...
              parseResponseTransform (val, name);
    endfunction

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
      fprintf ("\n  RegressionNeuralNetwork\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumPredictors', this.NumPredictors);
      str = repmat ({'%d'}, 1, numel (this.LayerSizes));
      str = strcat ('[', strjoin (str, ' '), ']');
      str = sprintf (str, this.LayerSizes);
      fprintf ("%+25s: %s\n", 'LayerSizes', str);
      if (iscellstr (this.Activations))
        str = repmat ({'''%s'''}, 1, numel (this.Activations));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.Activations{:});
        fprintf ("%+25s: %s\n", 'Activations', str);
      else # character vector
        fprintf ("%+25s: '%s'\n", 'Activations', this.Activations);
      endif
      fprintf ("%+25s: '%s'\n", 'OutputLayerActivation', ...
               this.OutputLayerActivation);
      fprintf ("%+25s: '%s'\n", 'ResponseTransform', this.ResponseTransform);
      fprintf ("%+25s: '%s'\n", 'Solver', this.Solver);
    endfunction

  endmethods

  methods(Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionNeuralNetwork} {@var{obj} =} RegressionNeuralNetwork (@var{X}, @var{Y})
    ## @deftypefnx {RegressionNeuralNetwork} {@var{obj} =} RegressionNeuralNetwork (@dots{}, @var{name}, @var{value})
    ##
    ## Create a @qcode{RegressionNeuralNetwork} object containing a neural
    ## network regression model.
    ##
    ## See the class documentation for the accepted @qcode{Name-Value} pairs.
    ##
    ## @seealso{fitrnet, RegressionNeuralNetwork}
    ## @end deftypefn
    function this = RegressionNeuralNetwork (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("RegressionNeuralNetwork: too few input arguments.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error (strcat ("RegressionNeuralNetwork: number of", ...
                       " rows in X and Y must be equal."));
      endif

      ## The response is continuous, so it must be numeric: a cellstr or a
      ## categorical response belongs to a classifier, and accepting one here
      ## would report a number for a question the model cannot answer.
      if (! (isnumeric (Y) && isreal (Y)))
        error (strcat ("RegressionNeuralNetwork: Y must be a", ...
                       " real numeric vector."));
      endif
      if (! (isvector (Y) || isempty (Y)))
        error ("RegressionNeuralNetwork: Y must be a vector.");
      endif

      ## Assign original X and Y data to the RegressionNeuralNetwork object
      this.X = X;
      this.Y = Y;

      ## Set default values before parsing optional parameters
      Standardize             = false;
      ResponseName            = [];
      PredictorNames          = [];
      LayerSizes              = 10;
      Activations             = 'relu';
      OutputLayerActivation   = 'none';
      LearningRate            = 0.003;
      IterationLimit          = 1000;
      DisplayInfo             = false;
      this.Solver = 'Gradient Descent';

      ## Supported activation functions.  'none' is MATLAB's name for the
      ## identity and is what a regression output layer wants.
      acList = {'linear', 'none', 'sigmoid', 'relu', 'tanh', ...
                'lrelu', 'prelu', 'elu', 'gelu'};

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'standardize'
            Standardize = varargin{2};
            if (! (Standardize == true || Standardize == false))
              error (strcat ("RegressionNeuralNetwork:", ...
                             " 'Standardize' must be either true or false."));
            endif

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("RegressionNeuralNetwork: 'PredictorNames'", ...
                             " must be supplied as a cellstring array."));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat ("RegressionNeuralNetwork: 'PredictorNames'", ...
                             " must have the same number of columns as X."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("RegressionNeuralNetwork: 'ResponseName'", ...
                             " must be a character vector."));
            endif

          case 'responsetransform'
            name = 'RegressionNeuralNetwork';
            [this.RTfun, this.ResponseTransform] = ...
                  parseResponseTransform (varargin{2}, name);

          case 'layersizes'
            LayerSizes = varargin{2};
            if (! (isnumeric (LayerSizes) && isvector (LayerSizes)
              && all (LayerSizes > 0) && all (mod (LayerSizes, 1) == 0)))
              error (strcat ("RegressionNeuralNetwork: 'LayerSizes'", ...
                             " must be a positive integer vector."));
            endif

          case 'learningrate'
            LearningRate = varargin{2};
            if (! (isnumeric (LearningRate) && isscalar (LearningRate) &&
                   LearningRate > 0))
              error (strcat ("RegressionNeuralNetwork:", ...
                             " 'LearningRate' must be a positive scalar."));
            endif

          case 'activations'
            Activations = varargin{2};
            if (! (ischar (Activations) || iscellstr (Activations)))
              error (strcat ("RegressionNeuralNetwork: 'Activations'", ...
                        " must be a character vector or a cellstring vector."));
            endif
            if (ischar (Activations))
              if (! any (strcmpi (Activations, acList)))
                error (strcat ("RegressionNeuralNetwork: unsupported", ...
                               " 'Activation' function."));
              endif
            else
              if (! all (cell2mat (cellfun (@(x) any (strcmpi (x, acList)),
                                   Activations, 'UniformOutput', false))))
                error (strcat ("RegressionNeuralNetwork: unsupported", ...
                               " 'Activation' functions."));
              endif
            endif
            Activations = tolower (Activations);

          case 'outputlayeractivation'
            OutputLayerActivation = varargin{2};
            if (! (ischar (OutputLayerActivation)))
              error (strcat ("RegressionNeuralNetwork:", ...
                       " 'OutputLayerActivation' must be a character vector."));
            endif
            if (! any (strcmpi (OutputLayerActivation, acList)))
              error (strcat ("RegressionNeuralNetwork: unsupported", ...
                             " 'OutputLayerActivation' function."));
            endif
            OutputLayerActivation = tolower (OutputLayerActivation);

          case 'iterationlimit'
            IterationLimit = varargin{2};
            if (! (isnumeric (IterationLimit) && isscalar (IterationLimit)
              && (IterationLimit > 0) && mod (IterationLimit, 1) == 0))
              error (strcat ("RegressionNeuralNetwork:", ...
                             " 'IterationLimit' must be a positive integer."));
            endif

          case 'displayinfo'
            DisplayInfo = varargin{2};
            if (! (DisplayInfo == true || DisplayInfo == false))
              error (strcat ("RegressionNeuralNetwork: 'DisplayInfo'", ...
                             " must be either true or false."));
            endif

          otherwise
            error (strcat ("RegressionNeuralNetwork: invalid",...
                           " parameter name in optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Generate default predictors and response variable names (if necessary)
      NumPredictors = columns (X);
      if (isempty (PredictorNames))
        for i = 1:NumPredictors
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = 'Y';
      endif

      ## Assign predictors and response variable names
      this.NumPredictors  = NumPredictors;
      this.PredictorNames = PredictorNames;
      this.ExpandedPredictorNames = PredictorNames;
      this.ResponseName   = ResponseName;

      ## An observation is dropped only when its response is missing.  A row
      ## whose predictors hold missing values is kept and reported as used,
      ## while the fit below draws on the complete observations alone.
      RowsUsed  = ! isnan (Y(:));
      Yret      = Y(RowsUsed);
      Xret      = X(RowsUsed, :);
      this.X    = Xret;
      this.Y    = Yret;
      cobs      = ! any (isnan (Xret), 2);
      Y         = Yret(cobs);
      X         = Xret(cobs, :);

      ## Check X and Y contain valid data
      if (! (isnumeric (X) && isfinite (X)))
        error ("RegressionNeuralNetwork: invalid values in X.");
      endif
      if (isempty (Y))
        error ("RegressionNeuralNetwork: Y cannot be empty.");
      endif
      if (! all (isfinite (Y)))
        error ("RegressionNeuralNetwork: invalid values in Y.");
      endif

      ## Assign the number of observations and their corresponding indices
      ## on the original data, which will be used for training the model
      this.NumObservations = rows (this.X);
      ## RowsUsed is left empty when every observation was used, as in MATLAB
      if (all (RowsUsed))
        this.RowsUsed = [];
      else
        this.RowsUsed = RowsUsed;
      endif

      ## Every observation carries the same weight
      this.W = ones (this.NumObservations, 1) / this.NumObservations;

      ## No predictor is treated as categorical, so the expanded names are
      ## the predictor names themselves.
      this.CategoricalPredictors = [];

      ## Handle the Standardize option.  The network must be trained on the
      ## scale it predicts on, so X is transformed here as well as in
      ## predict.
      if (Standardize)
        this.Sigma = std (X, [], 1);
        this.Sigma(this.Sigma == 0) = 1;  # predictor is constant
        this.Mu = mean (X, 1);
        X = (X - this.Mu) ./ this.Sigma;
      else
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
        ActivationCodes = ones (1, nlayers) * activationCode_ (Activations);
      elseif (nlayers != numel (Activations))
        error (strcat ("RegressionNeuralNetwork: 'Activations'", ...
                       " vector does not match the number of layers."));
      else
        ActivationCodes = [];
        for i = 1:nlayers
          code = activationCode_ (Activations{i});
          ActivationCodes = [ActivationCodes, code];
        endfor
      endif
      code = activationCode_ (OutputLayerActivation);
      ActivationCodes = [ActivationCodes, code];

      ## Start the training process.  LossFunction 2 is the mean squared error
      ## over a continuous response, which is what regression trains against.
      NumThreads = nproc ();
      Alpha = 0.01;  # used for ReLU and ELU activation layers
      rnn_timer_ = tic;
      Mdl = fcnntrain (X, Y(:), LayerSizes, ActivationCodes, NumThreads, ...
                       Alpha, LearningRate, IterationLimit, DisplayInfo, 2);

      ## Store training time and Loss
      ConvergenceInfo.Time = toc (rnn_timer_);
      ConvergenceInfo.TrainingLoss = Mdl.Loss;

      ## Remove redundant field
      Mdl = rmfield (Mdl, 'Loss');

      ## Save ModelParameters and ConvergenceInfo
      this.ModelParameters = Mdl;
      this.ConvergenceInfo = ConvergenceInfo;

      ## fcnntrain packs each neuron as [weights, bias] in one row, so the
      ## last column of every layer's matrix is its bias.
      nlay = numel (Mdl.LayerWeights);
      this.LayerWeights = cell (1, nlay);
      this.LayerBiases = cell (1, nlay);
      for i = 1:nlay
        Wb = Mdl.LayerWeights{i};
        this.LayerWeights{i} = Wb(:, 1:end-1);
        this.LayerBiases{i} = Wb(:, end);
      endfor

      ## Iteration by iteration record of the fit
      iter = (1:numel (ConvergenceInfo.TrainingLoss))';
      this.TrainingHistory = table (iter, ConvergenceInfo.TrainingLoss(:), ...
                                    'VariableNames', ...
                                    {'Iteration', 'TrainingLoss'});

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionNeuralNetwork} {@var{yFit} =} predict (@var{obj}, @var{XC})
    ##
    ## Predict the response for new data with a neural network regression
    ## model.
    ##
    ## @code{@var{yFit} = predict (@var{obj}, @var{XC})} returns a column
    ## vector holding the predicted response for each row of @var{XC}, using
    ## the network stored in @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionNeuralNetwork} class object.
    ## @item
    ## @var{XC} must be a numeric matrix with the same number of predictors as
    ## the data the model was trained on.
    ## @end itemize
    ##
    ## The transformation named by @code{ResponseTransform} is applied to the
    ## network's output before it is returned.
    ##
    ## @seealso{RegressionNeuralNetwork, fitrnet}
    ## @end deftypefn
    function yFit = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("RegressionNeuralNetwork.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("RegressionNeuralNetwork.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("RegressionNeuralNetwork.predict: XC must have", ...
                       " the same number of predictors as the trained model."));
      endif

      ## Standardize (if necessary)
      if (! isempty (this.Mu))
        XC = (XC - this.Mu) ./ this.Sigma;
      endif

      ## The network's output is its second return value: the first is an
      ## index of the largest output, which a single regression unit makes
      ## constant and meaningless.
      NumThreads = nproc ();
      [~, yFit] = fcnnpredict (this.ModelParameters, XC, NumThreads);

      ## Apply ResponseTransform
      yFit = this.RTfun (yFit);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionNeuralNetwork} {@var{yFit} =} resubPredict (@var{obj})
    ##
    ## Predict the response of the training data with a neural network
    ## regression model.
    ##
    ## @code{@var{yFit} = resubPredict (@var{obj})} returns a column vector
    ## holding the predicted response for every observation the model was
    ## trained on, that is the rows of @code{obj.X} selected by
    ## @code{obj.RowsUsed}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionNeuralNetwork} class object.
    ## @end itemize
    ##
    ## @seealso{RegressionNeuralNetwork, fitrnet}
    ## @end deftypefn
    function yFit = resubPredict (this)

      ## Get used rows
      XC = this.X;

      ## Standardize (if necessary)
      if (! isempty (this.Mu))
        XC = (XC - this.Mu) ./ this.Sigma;
      endif

      NumThreads = nproc ();
      [~, yFit] = fcnnpredict (this.ModelParameters, XC, NumThreads);

      ## Apply ResponseTransform
      yFit = this.RTfun (yFit);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionNeuralNetwork} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionNeuralNetwork} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the regression loss of a neural network model.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## weighted mean squared error between the response @var{Y} and the
    ## response the model predicts for @var{X}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionNeuralNetwork} class object.
    ## @item
    ## @var{X} must be a numeric matrix with the same number of predictors as
    ## the data the model was trained on.
    ## @item
    ## @var{Y} must be a numeric vector with as many rows as @var{X}.
    ## @end itemize
    ##
    ## @code{@var{L} = loss (@dots{}, @var{name}, @var{value})} accepts the
    ## following @qcode{Name-Value} pairs.
    ##
    ## @multitable @columnfractions 0.28 0.72
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab @qcode{'mse'}, the default, or a function
    ## handle called as @code{@var{lossfun} (@var{Y}, @var{yFit}, @var{W})}
    ## and returning a scalar.
    ##
    ## @item @qcode{'Weights'} @tab A numeric vector of observation weights
    ## with one entry per row of @var{X}.  It defaults to a uniform weight.
    ## The weights are normalized to sum to one before the loss is formed, so
    ## scaling them all by the same factor leaves the loss unchanged.
    ## @end multitable
    ##
    ## @seealso{RegressionNeuralNetwork, fitrnet}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("RegressionNeuralNetwork.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionNeuralNetwork.loss: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, 'loss');

      ## Defaults, then the optional pairs
      LossFun = 'mse';
      W = [];
      args = varargin;
      keep = true (1, numel (args));
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("RegressionNeuralNetwork.loss: parameter name", ...
                         " must be a character vector."));
        endif
        if (strcmpi (args{i}, 'lossfun'))
          LossFun = args{i+1};
          if (! (is_function_handle (LossFun) ||
                 (ischar (LossFun) && isrow (LossFun))))
            error (strcat ("RegressionNeuralNetwork.loss: 'LossFun' must", ...
                           " be a character vector or a function handle."));
          endif
          if (ischar (LossFun) && ! strcmpi (LossFun, 'mse'))
            error (strcat ("RegressionNeuralNetwork.loss: unsupported", ...
                           " 'LossFun' value."));
          endif
          keep(i:i+1) = false;
        endif
      endfor
      W = getWeights_ (this, args(keep), rows (X), 'loss');

      ## Weights are normalized to sum to one, as MATLAB does, so a loss is
      ## a weighted average rather than a weighted sum.
      W = W(:) / sum (W);
      yFit = predict (this, X);
      Y = Y(:);

      if (is_function_handle (LossFun))
        L = LossFun (Y, yFit, W);
        if (! (isnumeric (L) && isscalar (L)))
          error (strcat ("RegressionNeuralNetwork.loss: 'LossFun' must", ...
                         " return a numeric scalar."));
        endif
      else
        L = sum (W .* (Y - yFit) .^ 2);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionNeuralNetwork} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {RegressionNeuralNetwork} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the resubstitution regression loss of a neural network model.
    ##
    ## @code{@var{L} = resubLoss (@var{obj})} returns the weighted mean
    ## squared error of the model on the data it was trained on.  It accepts
    ## the same @qcode{Name-Value} pairs as @code{loss}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionNeuralNetwork} class object.
    ## @end itemize
    ##
    ## @seealso{RegressionNeuralNetwork, fitrnet}
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      used = true (rows (this.X), 1);
      X = this.X(used, :);
      Y = this.Y(used);
      L = loss (this, X, Y, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionNeuralNetwork} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {RegressionNeuralNetwork} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross validate a neural network regression model.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} returns a
    ## @qcode{RegressionPartitionedModel} holding one refit of @var{obj} per
    ## fold of a ten-fold partition, or of an @math{n}-fold one where the
    ## model has fewer than ten observations.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionNeuralNetwork} class object.
    ## @end itemize
    ##
    ## @code{@var{CVMdl} = crossval (@dots{}, @var{name}, @var{value})}
    ## accepts one, and only one, of the following @qcode{Name-Value} pairs.
    ##
    ## @multitable @columnfractions 0.28 0.72
    ## @headitem @var{Name} @tab @var{Value}
    ## @item @qcode{'KFold'} @tab An integer greater than 1, the number of
    ## folds.
    ## @item @qcode{'Holdout'} @tab A scalar in @math{(0, 1)}, the fraction
    ## of observations held out for testing.
    ## @item @qcode{'Leaveout'} @tab @qcode{'on'} or @qcode{'off'}, whether
    ## to hold out one observation at a time.
    ## @item @qcode{'CVPartition'} @tab A @code{cvpartition} object over as
    ## many observations as the model was trained on.
    ## @end multitable
    ##
    ## @seealso{RegressionNeuralNetwork, RegressionPartitionedModel,
    ## cvpartition}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)

      if (numel (varargin) == 1)
        error (strcat ("RegressionNeuralNetwork.crossval: Name-Value", ...
                       " arguments must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("RegressionNeuralNetwork.crossval: specify only", ...
                       " one of the optional Name-Value paired arguments."));
      endif

      if (this.NumObservations < 10)
        numFolds  = this.NumObservations;
      else
        numFolds  = 10;
      endif
      Holdout     = [];
      Leaveout    = 'off';
      CVPartition = [];

      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'kfold'
            numFolds = varargin{2};
            if (! (isnumeric (numFolds) && isscalar (numFolds)
                   && (numFolds == fix (numFolds)) && numFolds > 1))
              error (strcat ("RegressionNeuralNetwork.crossval: 'KFold'", ...
                             " must be an integer value greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("RegressionNeuralNetwork.crossval: 'Holdout'", ...
                             " must be a numeric value between 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat ("RegressionNeuralNetwork.crossval: 'Leaveout'", ...
                             " must be either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("RegressionNeuralNetwork.crossval:", ...
                             " 'CVPartition' must be a 'cvpartition'", ...
                             " object."));
            endif

          otherwise
            error (strcat ("RegressionNeuralNetwork.crossval: invalid", ...
                           " parameter name in optional paired arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Determine the cross-validation method to use.  The partition is
      ## built over the observations actually trained on, so its indices and
      ## the partitioned model's rows are the same set.
      n = this.NumObservations;
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (n, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (n, 'LeaveOut');
      else
        partition = cvpartition (n, 'KFold', numFolds);
      endif

      ## Create a cross-validated model object
      CVMdl = RegressionPartitionedModel (this, partition);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionNeuralNetwork} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Create a @qcode{CompactRegressionNeuralNetwork} object.
    ##
    ## @code{@var{CMdl} = compact (@var{obj})} returns a compact version of
    ## the @qcode{RegressionNeuralNetwork} object @var{obj}, which keeps the
    ## trained network but drops the training data, so it predicts identically
    ## while carrying no observations.
    ##
    ## @seealso{fitrnet, RegressionNeuralNetwork,
    ## CompactRegressionNeuralNetwork}
    ## @end deftypefn
    function CMdl = compact (this)
      ## Create a compact model
      CMdl = CompactRegressionNeuralNetwork (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionNeuralNetwork} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a neural network regression model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves every property of
    ## the @qcode{RegressionNeuralNetwork} object @var{obj} into
    ## @var{filename} in binary format, so that it can be read back with
    ## @code{loadmodel}.
    ##
    ## @seealso{loadmodel, RegressionNeuralNetwork, fitrnet}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error ("RegressionNeuralNetwork.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("RegressionNeuralNetwork.savemodel: FNAME must be", ...
                       " a character vector."));
      endif
      ## Generate variable for class name
      classdef_name = 'RegressionNeuralNetwork';

      ## Create variables from model properties
      X = this.X;
      Y = this.Y;
      NumObservations         = this.NumObservations;
      RowsUsed                = this.RowsUsed;
      BinEdges            = this.BinEdges;
      NumPredictors           = this.NumPredictors;
      PredictorNames          = this.PredictorNames;
      ResponseName            = this.ResponseName;
      ResponseTransform       = this.ResponseTransform;
      Sigma                   = this.Sigma;
      Mu                      = this.Mu;
      LayerSizes              = this.LayerSizes;
      Activations             = this.Activations;
      OutputLayerActivation   = this.OutputLayerActivation;
      LearningRate            = this.LearningRate;
      IterationLimit          = this.IterationLimit;
      ModelParameters         = this.ModelParameters;
      ConvergenceInfo         = this.ConvergenceInfo;
      DisplayInfo             = this.DisplayInfo;
      Solver                  = this.Solver;
      LayerWeights            = this.LayerWeights;
      LayerBiases             = this.LayerBiases;
      W                       = this.W;
      CategoricalPredictors   = this.CategoricalPredictors;
      ExpandedPredictorNames  = this.ExpandedPredictorNames;
      RTfun                  = this.RTfun;

      ## TrainingHistory is a table, and Octave cannot save a classdef object
      ## to a binary file, so it is left out here and rebuilt on loading from
      ## ConvergenceInfo, which holds the same numbers as a plain vector.

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'NumObservations', ...
            'RowsUsed', 'BinEdges', 'NumPredictors', 'PredictorNames', ...
            'ResponseName', ...
            'ResponseTransform', 'Sigma', 'Mu', ...
            'LayerSizes', 'Activations', 'OutputLayerActivation', ...
            'LearningRate', 'IterationLimit', 'Solver', 'ModelParameters', ...
            'ConvergenceInfo', 'DisplayInfo', 'LayerWeights', 'LayerBiases', ...
            'W', 'CategoricalPredictors', 'ExpandedPredictorNames', 'RTfun');
    endfunction

  endmethods

  methods (Access = private)

    ## Shared validation for the assessment methods, so each reports under
    ## its own name.
    function [X, Y] = checkXY_ (this, X, Y, caller)
      if (isempty (X))
        error ("RegressionNeuralNetwork.%s: X is empty.", caller);
      elseif (this.NumPredictors != columns (X))
        error (strcat ("RegressionNeuralNetwork.%s: X must have the", ...
                       " same number of predictors as the trained model."), ...
               caller);
      endif
      if (isempty (Y))
        error ("RegressionNeuralNetwork.%s: Y is empty.", caller);
      elseif (! (isnumeric (Y) && isreal (Y)))
        error (strcat ("RegressionNeuralNetwork.%s: Y must be a real", ...
                       " numeric vector."), caller);
      elseif (rows (X) != numel (Y))
        error (strcat ("RegressionNeuralNetwork.%s: Y must have the", ...
                       " same number of rows as X."), caller);
      endif
    endfunction

    ## Pull a "Weights" pair out of the optional arguments, defaulting to a
    ## uniform weight, and reject any other name.
    function W = getWeights_ (this, args, n, caller)
      W = ones (n, 1);
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("RegressionNeuralNetwork.%s: parameter name", ...
                         " must be a character vector."), caller);
        endif
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
          if (! (isnumeric (W) && isvector (W)))
            error (strcat ("RegressionNeuralNetwork.%s: 'Weights'", ...
                           " must be a numeric vector."), caller);
          endif
          if (numel (W) != n)
            error (strcat ("RegressionNeuralNetwork.%s: size of", ...
                           " 'Weights' must equal the number of", ...
                           " rows in X."), caller);
          endif
        else
          error (strcat ("RegressionNeuralNetwork.%s: invalid", ...
                         " parameter name in optional paired", ...
                         " arguments."), caller);
        endif
      endfor
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a RegressionNeuralNetwork object
      mdl = RegressionNeuralNetwork (1, 1);

      ## Get fieldnames from DATA (including private properties)
      names = fieldnames (data);
      ## The set methods for these read other properties, and one of them
      ## rebuilds Coeffs, so they are assigned once everything else is in
      ## place rather than in the order the file happens to list them.
      late = ismember (names, {'Cost', 'Prior', 'ScoreTransform', ...
                               'ResponseTransform'});
      names = [names(! late); names(late)];

      ## Copy data into object
      for i = 1:numel (names)
        ## Check fieldnames in DATA match properties in RegressionNeuralNetwork
        try
          mdl.(names{i}) = data.(names{i});
        catch
          error (strcat ("RegressionNeuralNetwork.load_model:", ...
                         " invalid model in '%s'."), filename)
        end_try_catch
      endfor

      ## Rebuild the TrainingHistory table, which savemodel cannot write out
      iter = (1:numel (mdl.ConvergenceInfo.TrainingLoss))';
      mdl.TrainingHistory = table (iter, ...
                                   mdl.ConvergenceInfo.TrainingLoss(:), ...
                                   'VariableNames', ...
                                   {'Iteration', 'TrainingLoss'});
    endfunction

  endmethods

endclassdef

## Map an activation name to the code fcnntrain expects.  'none' and 'linear'
## are the same identity map; MATLAB spells the output layer's 'none'.
function numCode = activationCode_ (strCode)
  switch (strCode)
    case {'linear', 'none'}
      numCode = 0;
    case 'sigmoid'
      numCode = 1;
    case 'relu'
      numCode = 2;
    case 'tanh'
      numCode = 3;
    case {'lrelu', 'prelu'}
      numCode = 5;
    case 'elu'
      numCode = 6;
    case 'gelu'
      numCode = 7;
    otherwise
      error (strcat ("RegressionNeuralNetwork: misspelling or unsupported", ...
                     " 'Activation' function: '%s'."), strCode);
  endswitch
endfunction


## A fitted model carries the defaults MATLAB documents for fitrnet.
%!test
%! rand ('seed', 42);
%! X = linspace (-1, 1, 40)';
%! Y = 2 * X + 0.5;
%! Mdl = RegressionNeuralNetwork (X, Y, 'IterationLimit', 50);
%! assert_equal (class (Mdl), 'RegressionNeuralNetwork');
%! assert_equal (Mdl.LayerSizes, 10);
%! assert_equal (Mdl.Activations, 'relu');
%! assert_equal (Mdl.OutputLayerActivation, 'none');
%! assert_equal (Mdl.LearningRate, 0.003);
%! assert_equal (Mdl.IterationLimit, 50);
%! assert_equal (isempty (Mdl.Mu), true);
%! assert_equal (Mdl.Solver, 'Gradient Descent');
%! assert_equal (Mdl.NumObservations, 40);
%! assert_equal (Mdl.NumPredictors, 1);
%! assert_equal (Mdl.ResponseName, 'Y');
%! assert_equal (Mdl.PredictorNames, {'x1'});

## The output layer has one unit, and the identity leaves it unbounded, so a
## prediction is a real number rather than a score.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = linspace (-2, 2, 60)';
%! Y = 30 * X + 100;
%! Mdl = RegressionNeuralNetwork (X, Y, 'IterationLimit', 400);
%! assert_equal (rows (Mdl.LayerWeights{end}), 1);
%! yFit = predict (Mdl, X);
%! assert_equal (size (yFit), [60, 1]);
%! assert_equal (max (yFit) > 50, true);

## A network recovers a smooth function to within the noise on it.
%!test
%! rand ('seed', 7); randn ('seed', 7);
%! X = linspace (-2, 2, 80)';
%! Y = 3 * X.^2 - 1 + randn (80, 1) * 0.05;
%! Mdl = RegressionNeuralNetwork (X, Y, 'LayerSizes', [12, 12], ...
%!                                'IterationLimit', 600);
%! assert_equal (sqrt (resubLoss (Mdl)) < 0.2, true);

## The recorded history is the network's own loss, not a running average.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = linspace (0, 1, 50)';
%! Y = 4 * X - 2;
%! Mdl = RegressionNeuralNetwork (X, Y, 'IterationLimit', 200);
%! h = Mdl.TrainingHistory;
%! assert_equal (class (h), 'table');
%! assert_equal (h.Properties.VariableNames, {'Iteration', 'TrainingLoss'});
%! assert_equal (rows (h), 200);
%! assert_equal (h.Iteration', 1:200);
%! assert_equal (h.TrainingLoss(end) < h.TrainingLoss(1), true);
%! assert_equal (h.TrainingLoss(end), resubLoss (Mdl), 1e-12);

## ConvergenceInfo carries the same loss and the time the fit took.
%!test
%! rand ('seed', 42);
%! X = linspace (0, 1, 30)';
%! Mdl = RegressionNeuralNetwork (X, 2 * X, 'IterationLimit', 25);
%! assert_equal (fieldnames (Mdl.ConvergenceInfo), {'Time'; 'TrainingLoss'});
%! assert_equal (numel (Mdl.ConvergenceInfo.TrainingLoss), 25);
%! assert_equal (Mdl.ConvergenceInfo.Time > 0, true);

## predict on the training rows is resubPredict.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = [randn(30, 2); randn(30, 2) + 3];
%! Y = X(:,1) - 2 * X(:,2);
%! Mdl = RegressionNeuralNetwork (X, Y, 'IterationLimit', 100);
%! assert_equal (predict (Mdl, X), resubPredict (Mdl));

## Standardize trains on the scale it predicts on.  A model trained on raw
## data and asked about standardized data is not merely worse, it is wrong.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = [randn(60, 1), randn(60, 1) * 1000];
%! Y = X(:,1) + X(:,2) / 1000;
%! Mdl = RegressionNeuralNetwork (X, Y, 'Standardize', true, ...
%!                                'IterationLimit', 300);
%! assert_equal (size (Mdl.Mu), [1, 2]);
%! assert_equal (size (Mdl.Sigma), [1, 2]);
%! assert_equal (predict (Mdl, X), resubPredict (Mdl));
%! assert_equal (sqrt (resubLoss (Mdl)) < std (Y), true);

## A constant predictor gets a unit scale rather than a division by zero.
%!test
%! rand ('seed', 42);
%! X = [linspace(0, 1, 20)', ones(20, 1)];
%! Mdl = RegressionNeuralNetwork (X, X(:,1), 'Standardize', true, ...
%!                                'IterationLimit', 50);
%! assert_equal (Mdl.Sigma(2), 1);
%! assert_equal (all (isfinite (resubPredict (Mdl))), true);

## Layer sizes and per-layer activations reach the network.
%!test
%! rand ('seed', 42);
%! X = linspace (0, 1, 30)';
%! Mdl = RegressionNeuralNetwork (X, 2 * X, 'LayerSizes', [4, 6], ...
%!                                'Activations', {'tanh', 'sigmoid'}, ...
%!                                'IterationLimit', 20);
%! assert_equal (Mdl.LayerSizes, [4, 6]);
%! assert_equal (Mdl.Activations, {'tanh', 'sigmoid'});
%! assert_equal (numel (Mdl.LayerWeights), 3);
%! assert_equal (size (Mdl.LayerWeights{1}), [4, 1]);
%! assert_equal (size (Mdl.LayerWeights{2}), [6, 4]);
%! assert_equal (size (Mdl.LayerWeights{3}), [1, 6]);
%! assert_equal (Mdl.ModelParameters.Activations, [3, 1, 0]);

## 'none' and 'linear' name the same identity output.
%!test
%! rand ('seed', 42);
%! X = linspace (0, 1, 20)';
%! M1 = RegressionNeuralNetwork (X, 2 * X, 'OutputLayerActivation', 'none', ...
%!                               'IterationLimit', 10);
%! M2 = RegressionNeuralNetwork (X, 2 * X, 'OutputLayerActivation', ...
%!                               'linear', 'IterationLimit', 10);
%! assert_equal (M1.ModelParameters.Activations(end), 0);
%! assert_equal (M2.ModelParameters.Activations(end), 0);

## Rows carrying a missing value are dropped from both X and Y.
%!test
%! rand ('seed', 42);
%! X = [linspace(0, 1, 12)'; NaN; 0.5];
%! Y = [2 * linspace(0, 1, 12)'; 1; NaN];
%! Mdl = RegressionNeuralNetwork (X, Y, 'IterationLimit', 20);
%! assert_equal (Mdl.NumObservations, 13);
%! assert_equal (sum (Mdl.RowsUsed), 13);
%! assert_equal (Mdl.RowsUsed(13:14), [true; false]);
%! assert_equal (numel (resubPredict (Mdl)), 13);

## Observation weights default to a uniform weight summing to one.
%!test
%! rand ('seed', 42);
%! X = linspace (0, 1, 25)';
%! Mdl = RegressionNeuralNetwork (X, 2 * X, 'IterationLimit', 10);
%! assert_equal (size (Mdl.W), [25, 1]);
%! assert_equal (sum (Mdl.W), 1, 1e-12);
%! assert_equal (Mdl.W, ones (25, 1) / 25, 1e-12);

## loss defaults to the weighted mean squared error, and weights are
## normalized, so scaling every weight leaves the loss alone.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = linspace (0, 1, 30)';
%! Y = 3 * X + 1;
%! Mdl = RegressionNeuralNetwork (X, Y, 'IterationLimit', 100);
%! yFit = predict (Mdl, X);
%! assert_equal (loss (Mdl, X, Y), mean ((Y - yFit) .^ 2), 1e-12);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'mse'), loss (Mdl, X, Y), 1e-12);
%! w = rand (30, 1) + 0.1;
%! assert_equal (loss (Mdl, X, Y, 'Weights', w), ...
%!               loss (Mdl, X, Y, 'Weights', 7 * w), 1e-12);
%! assert_equal (loss (Mdl, X, Y, 'Weights', w), ...
%!               sum ((w / sum (w)) .* (Y - yFit) .^ 2), 1e-12);

## loss takes a function handle of the response, the fit and the weights.
%!test
%! rand ('seed', 42);
%! X = linspace (0, 1, 20)';
%! Y = 2 * X;
%! Mdl = RegressionNeuralNetwork (X, Y, 'IterationLimit', 50);
%! f = @(y, yf, w) sum (w .* abs (y - yf));
%! yFit = predict (Mdl, X);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', f), ...
%!               mean (abs (Y - yFit)), 1e-12);

## resubLoss is loss on the training data.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = linspace (0, 1, 21)';
%! Y = [3 * linspace(0, 1, 20)'; NaN];
%! Mdl = RegressionNeuralNetwork (X, Y, 'IterationLimit', 80);
%! Xu = X(Mdl.RowsUsed, :);
%! Yu = Y(Mdl.RowsUsed);
%! assert_equal (resubLoss (Mdl), loss (Mdl, Xu, Yu), 1e-12);
%! assert_equal (resubLoss (Mdl, 'Weights', ones (20, 1)), ...
%!               loss (Mdl, Xu, Yu), 1e-12);

## ResponseTransform is applied to the prediction, by name or by handle.
%!test
%! rand ('seed', 42);
%! X = linspace (0, 1, 20)';
%! Y = 2 * X + 1;
%! Mdl = RegressionNeuralNetwork (X, Y, 'IterationLimit', 50);
%! raw = predict (Mdl, X);
%! Mdl.ResponseTransform = 'exp';
%! assert_equal (predict (Mdl, X), exp (raw), 1e-12);
%! Mdl.ResponseTransform = @(y) 2 * y;
%! assert_equal (predict (Mdl, X), 2 * raw, 1e-12);
%! Mdl.ResponseTransform = 'none';
%! assert_equal (predict (Mdl, X), raw, 1e-12);

## The transform set at construction is the one predict uses.
%!test
%! rand ('seed', 42);
%! X = linspace (0, 1, 20)';
%! Mdl = RegressionNeuralNetwork (X, 2 * X, 'ResponseTransform', 'identity', ...
%!                                'IterationLimit', 20);
%! assert_equal (class (Mdl.ResponseTransform), 'char');
%! assert_equal (Mdl.RTfun (5), 5);

## A saved model comes back carrying its own numbers.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = linspace (0, 1, 30)';
%! Y = 4 * X - 1;
%! Mdl = RegressionNeuralNetwork (X, Y, 'LayerSizes', [6, 4], ...
%!                                'IterationLimit', 60);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'RegressionNeuralNetwork');
%! assert_equal (M2.LayerWeights, Mdl.LayerWeights);
%! assert_equal (M2.LayerBiases, Mdl.LayerBiases);
%! assert_equal (M2.NumObservations, Mdl.NumObservations);
%! assert_equal (M2.LayerSizes, Mdl.LayerSizes);
%! assert_equal (M2.ResponseName, Mdl.ResponseName);
%! assert_equal (M2.W, Mdl.W);
%! assert_equal (predict (M2, X), predict (Mdl, X));
%! assert_equal (rows (M2.TrainingHistory), 60);
%! assert_equal (M2.TrainingHistory.TrainingLoss, ...
%!               Mdl.TrainingHistory.TrainingLoss);

## compact drops the training data but predicts identically.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = [randn(40, 2); randn(40, 2) + 2];
%! Y = X(:,1) - X(:,2);
%! Mdl = RegressionNeuralNetwork (X, Y, 'IterationLimit', 100);
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), 'CompactRegressionNeuralNetwork');
%! assert_equal (predict (CMdl, X), predict (Mdl, X));
%! assert_equal (loss (CMdl, X, Y), loss (Mdl, X, Y));

## crossval returns a partitioned model holding one fit per fold.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (30, 2);
%! Y = X(:,1) - X(:,2);
%! Mdl = fitrnet (X, Y, 'IterationLimit', 20);
%! CVMdl = crossval (Mdl, 'KFold', 3);
%! assert_equal (class (CVMdl), 'RegressionPartitionedModel');
%! assert_equal (CVMdl.KFold, 3);
%! assert_equal (CVMdl.CrossValidatedModel, 'RegressionNeuralNetwork');
%! assert_equal (numel (kfoldPredict (CVMdl)), 30);
%! assert_equal (isfinite (kfoldLoss (CVMdl)), true);

## Test input validation for crossval
%!error<RegressionNeuralNetwork.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (fitrnet (randn (12, 2), randn (12, 1), ...
%!                    'IterationLimit', 20), 'KFold')
%!error<RegressionNeuralNetwork.crossval: specify only one of the optional Name-Value paired arguments.> ...
%! crossval (fitrnet (randn (12, 2), randn (12, 1), ...
%!                    'IterationLimit', 20), ...
%!           'KFold', 3, 'Leaveout', 'on')
%!error<RegressionNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (fitrnet (randn (12, 2), randn (12, 1), ...
%!                    'IterationLimit', 20), 'KFold', 1)
%!error<RegressionNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (fitrnet (randn (12, 2), randn (12, 1), ...
%!                    'IterationLimit', 20), 'Holdout', 1)
%!error<RegressionNeuralNetwork.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (fitrnet (randn (12, 2), randn (12, 1), ...
%!                    'IterationLimit', 20), 'Leaveout', 1)
%!error<RegressionNeuralNetwork.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (fitrnet (randn (12, 2), randn (12, 1), ...
%!                    'IterationLimit', 20), 'CVPartition', 1)
%!error<RegressionNeuralNetwork.crossval: invalid parameter name in optional paired arguments.> ...
%! crossval (fitrnet (randn (12, 2), randn (12, 1), ...
%!                    'IterationLimit', 20), 'Nope', 1)

## Test input validation for the constructor
%!error<RegressionNeuralNetwork: too few input arguments.> ...
%! RegressionNeuralNetwork ()
%!error<RegressionNeuralNetwork: too few input arguments.> ...
%! RegressionNeuralNetwork (ones (10, 2))
%!error<RegressionNeuralNetwork: number of rows in X and Y must be equal.> ...
%! RegressionNeuralNetwork (ones (10, 2), ones (5, 1))
%!error<RegressionNeuralNetwork: Y must be a real numeric vector.> ...
%! RegressionNeuralNetwork (ones (5, 2), {'a'; 'b'; 'c'; 'd'; 'e'})
%!error<RegressionNeuralNetwork: Y must be a real numeric vector.> ...
%! RegressionNeuralNetwork (ones (5, 2), complex (ones (5, 1)))
%!error<RegressionNeuralNetwork: Y must be a vector.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 3))
%!error<RegressionNeuralNetwork: invalid values in X.> ...
%! RegressionNeuralNetwork ([1; Inf; 3], [1; 2; 3])
%!error<RegressionNeuralNetwork: invalid values in Y.> ...
%! RegressionNeuralNetwork ([1; 2; 3], [1; Inf; 3])
%!error<RegressionNeuralNetwork: 'Standardize' must be either true or false.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'Standardize', 'yes')
%!error<RegressionNeuralNetwork: 'PredictorNames' must be supplied as a cellstring array.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'PredictorNames', 'a')
%!error<RegressionNeuralNetwork: 'PredictorNames' must have the same number of columns as X.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'PredictorNames', {'a'})
%!error<RegressionNeuralNetwork: 'ResponseName' must be a character vector.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'ResponseName', 5)
%!error<RegressionNeuralNetwork: 'ResponseTransform' must be a character vector or a function handle.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'ResponseTransform', 5)
%!error<RegressionNeuralNetwork: unrecognized 'ResponseTransform' function.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), ...
%!                          'ResponseTransform', 'nope')
%!error<RegressionNeuralNetwork: function handle for 'ResponseTransform' must return the same size as its input.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), ...
%!                          'ResponseTransform', @(y) [y; y])
%!error<RegressionNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'LayerSizes', -1)
%!error<RegressionNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'LayerSizes', 2.5)
%!error<RegressionNeuralNetwork: 'LearningRate' must be a positive scalar.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'LearningRate', 0)
%!error<RegressionNeuralNetwork: 'LearningRate' must be a positive scalar.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), ...
%!                          'LearningRate', [0.1, 0.2])
%!error<RegressionNeuralNetwork: 'Activations' must be a character vector or a cellstring vector.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'Activations', 5)
%!error<RegressionNeuralNetwork: unsupported 'Activation' function.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'Activations', 'softmax')
%!error<RegressionNeuralNetwork: unsupported 'Activation' functions.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), ...
%!                          'Activations', {'relu', 'nope'})
%!error<RegressionNeuralNetwork: 'OutputLayerActivation' must be a character vector.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), ...
%!                          'OutputLayerActivation', 5)
%!error<RegressionNeuralNetwork: unsupported 'OutputLayerActivation' function.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), ...
%!                          'OutputLayerActivation', 'softmax')
%!error<RegressionNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'IterationLimit', 0)
%!error<RegressionNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'IterationLimit', 2.5)
%!error<RegressionNeuralNetwork: 'DisplayInfo' must be either true or false.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'DisplayInfo', 'yes')
%!error<RegressionNeuralNetwork: invalid parameter name in optional pair arguments.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'Prior', 1)
%!error<RegressionNeuralNetwork: 'Activations' vector does not match the number of layers.> ...
%! RegressionNeuralNetwork (ones (5, 2), ones (5, 1), 'LayerSizes', [4, 4], ...
%!                          'Activations', {'relu', 'relu', 'relu'})

## Test input validation for predict
%!error<RegressionNeuralNetwork.predict: too few input arguments.> ...
%! predict (RegressionNeuralNetwork ([1; 2; 3], [1; 2; 3], ...
%!                                   'IterationLimit', 5))
%!error<RegressionNeuralNetwork.predict: XC is empty.> ...
%! predict (RegressionNeuralNetwork ([1; 2; 3], [1; 2; 3], ...
%!                                   'IterationLimit', 5), [])
%!error<RegressionNeuralNetwork.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (RegressionNeuralNetwork ([1; 2; 3], [1; 2; 3], ...
%!                                   'IterationLimit', 5), ones (2, 3))

## Test input validation for loss
%!shared RNNMdl
%! rand ('seed', 42);
%! RNNMdl = RegressionNeuralNetwork ([1; 2; 3; 4], [2; 4; 6; 8], ...
%!                                   'IterationLimit', 10);
%!error<RegressionNeuralNetwork.loss: too few input arguments.> ...
%! loss (RNNMdl)
%!error<RegressionNeuralNetwork.loss: too few input arguments.> ...
%! loss (RNNMdl, [1; 2])
%!error<RegressionNeuralNetwork.loss: Name-Value arguments must be in pairs.> ...
%! loss (RNNMdl, [1; 2], [2; 4], 'Weights')
%!error<RegressionNeuralNetwork.loss: X is empty.> ...
%! loss (RNNMdl, [], [2; 4])
%!error<RegressionNeuralNetwork.loss: X must have the same number of predictors as the trained model.> ...
%! loss (RNNMdl, ones (2, 3), [2; 4])
%!error<RegressionNeuralNetwork.loss: Y is empty.> ...
%! loss (RNNMdl, [1; 2], [])
%!error<RegressionNeuralNetwork.loss: Y must be a real numeric vector.> ...
%! loss (RNNMdl, [1; 2], {'a'; 'b'})
%!error<RegressionNeuralNetwork.loss: Y must have the same number of rows as X.> ...
%! loss (RNNMdl, [1; 2], [2; 4; 6])
%!error<RegressionNeuralNetwork.loss: 'LossFun' must be a character vector or a function handle.> ...
%! loss (RNNMdl, [1; 2], [2; 4], 'LossFun', 5)
%!error<RegressionNeuralNetwork.loss: unsupported 'LossFun' value.> ...
%! loss (RNNMdl, [1; 2], [2; 4], 'LossFun', 'mae')
%!error<RegressionNeuralNetwork.loss: 'LossFun' must return a numeric scalar.> ...
%! loss (RNNMdl, [1; 2], [2; 4], 'LossFun', @(y, yf, w) [1, 2])
%!error<RegressionNeuralNetwork.loss: 'Weights' must be a numeric vector.> ...
%! loss (RNNMdl, [1; 2], [2; 4], 'Weights', {'a'})
%!error<RegressionNeuralNetwork.loss: size of 'Weights' must equal the number of rows in X.> ...
%! loss (RNNMdl, [1; 2], [2; 4], 'Weights', [1; 2; 3])
%!error<RegressionNeuralNetwork.loss: invalid parameter name in optional paired arguments.> ...
%! loss (RNNMdl, [1; 2], [2; 4], 'Nope', 1)
%!error<RegressionNeuralNetwork.loss: parameter name must be a character vector.> ...
%! loss (RNNMdl, [1; 2], [2; 4], 5, 1)

## Test input validation for savemodel
%!error<RegressionNeuralNetwork.savemodel: too few input arguments.> ...
%! savemodel (RNNMdl)
%!error<RegressionNeuralNetwork.savemodel: FNAME must be a character vector.> ...
%! savemodel (RNNMdl, 5)

%!error<RegressionNeuralNetwork: unrecognized 'ResponseTransform' function.> ...
%! RNNMdl.ResponseTransform = 'nope';

## RowsUsed is empty when every observation was used.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! Y = meas(:,1);
%! Mdl = fitrnet (X, Y, 'IterationLimit', 20);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (class (Mdl.RowsUsed), 'double');
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (rows (Mdl.X), 150);
%! assert_equal (rows (Mdl.W), 150);

## A missing response drops its observation and RowsUsed marks it.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! Y = meas(:,1);
%! Y(5) = NaN;
%! Mdl = fitrnet (X, Y, 'IterationLimit', 20);
%! assert_equal (class (Mdl.RowsUsed), 'logical');
%! assert_equal (size (Mdl.RowsUsed), [150, 1]);
%! assert_equal (sum (Mdl.RowsUsed), 149);
%! assert_equal (Mdl.RowsUsed(5), false);
%! assert_equal (Mdl.NumObservations, 149);
%! assert_equal (rows (Mdl.X), 149);
%! assert_equal (rows (Mdl.W), 149);

## A missing predictor keeps its observation, so RowsUsed stays empty.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! X(3,2) = NaN;
%! Y = meas(:,1);
%! Mdl = fitrnet (X, Y, 'IterationLimit', 20);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (rows (Mdl.X), 150);
%! assert_equal (sum (isnan (Mdl.X(:))), 1);

## Standardizing summarizes the complete observations.  With no classes the
## weights are uniform over them.  Values from MATLAB R2024a.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! X(7,2) = NaN; X(120,3) = NaN;
%! Mdl = fitrnet (X, meas(:,1), 'Standardize', true, 'IterationLimit', 20);
%! assert_equal (Mdl.Mu, [3.0608108108108096, 3.7655405405405395, ...
%!                        1.203378378378378], 1e-13);
%! assert_equal (Mdl.Sigma, [0.43214937187299296, 1.7636045663278643, ...
%!                           0.76339873235638711], 1e-13);

## A fitted model survives savemodel and loadmodel: the properties come
## back as they were and it predicts the same.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! Y = meas(:,1);
%! Mdl = fitrnet (X, Y, 'IterationLimit', 20);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'RegressionNeuralNetwork');
%! assert_equal (M2.NumObservations, Mdl.NumObservations);
%! assert_equal (M2.PredictorNames, Mdl.PredictorNames);
%! assert_equal (class (M2.ResponseTransform), class (Mdl.ResponseTransform));
%! assert_equal (predict (M2, X(1:5,:)), predict (Mdl, X(1:5,:)), 1e-12);

## BinEdges is an empty cell, which is what MATLAB reports for this
## learner as well: it fits the predictors as they are.
%!test
%! load fisheriris
%! Mdl = fitrnet (meas(:,1:3), meas(:,4), 'IterationLimit', 10);
%! assert_equal (class (Mdl.BinEdges), 'cell');
%! assert_equal (Mdl.BinEdges, {});
