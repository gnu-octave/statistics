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
  ## @deftp {statistics} ClassificationNeuralNetwork
  ##
  ## Neural network classification
  ##
  ## The @code{ClassificationNeuralNetwork} class implements a neural network
  ## classifier object, which can predict responses for new data using the
  ## @code{predict} method.
  ##
  ## Neural network classification is a machine learning method that uses
  ## interconnected nodes in multiple layers to learn complex patterns in data.
  ## It processes inputs through hidden layers with activation functions to
  ## produce classification outputs.
  ##
  ## Create a @code{ClassificationNeuralNetwork} object by using the
  ## @code{fitcnet} function or the class constructor.
  ##
  ## @seealso{fitcnet}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix containing the unstandardized predictor data.  Each
    ## column of @var{X} represents one predictor (variable), and each row
    ## represents one observation.  This property is read-only.
    ##
    ## @end deftp
    X                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} Y
    ##
    ## Class labels
    ##
    ## Specified as a logical or numeric column vector, or as a character array
    ## or a cell array of character vectors with the same number of rows as the
    ## predictor data.  Each row in @var{Y} is the observed class label for
    ## the corresponding row in @var{X}.  This property is read-only.
    ##
    ## @end deftp
    Y                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer value specifying the number of observations in the
    ## training dataset used for training the ClassificationNeuralNetwork model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations       = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} RowsUsed
    ##
    ## Rows used for fitting
    ##
    ## A logical column vector with the same length as the observations in the
    ## original predictor data @var{X}, true for each row that was used for
    ## fitting the ClassificationNeuralNetwork model.  It is empty, @qcode{[]},
    ## when every observation was used, so a non-empty value means that rows
    ## holding missing values were dropped.  This property is read-only.
    ##
    ## @end deftp
    RowsUsed              = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer value specifying the number of predictors in the
    ## training dataset used for training the ClassificationNeuralNetwork model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumPredictors         = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which they appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName          = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} ClassNames
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
    ClassNames            = [];


    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} Standardize
    ##
    ## Flag to standardize predictors
    ##
    ## A boolean flag indicating whether the predictor data has been
    ## standardized prior to training. When @qcode{true}, the predictors are
    ## centered and scaled to have zero mean and unit variance. This property is
    ## read-only.
    ##
    ## @end deftp
    Standardize           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} Sigma
    ##
    ## Predictor standard deviations
    ##
    ## A numeric vector containing the standard deviations of the predictors
    ## used for standardization.  Empty if @qcode{Standardize} is @qcode{false}.
    ## This property is read-only.
    ##
    ## Only observations with no missing predictor enter the estimate, and
    ## they are weighted so that each class keeps the share of the
    ## observation weight it carried before any row was set aside.
    ##
    ## @end deftp
    Sigma                 = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} Mu
    ##
    ## Predictor means
    ##
    ## A numeric vector containing the means of the predictors used for
    ## standardization.  Empty if @qcode{Standardize} is @qcode{false}.
    ## This property is read-only.
    ##
    ## Only observations with no missing predictor enter the estimate, and
    ## they are weighted so that each class keeps the share of the
    ## observation weight it carried before any row was set aside.
    ##
    ## @end deftp
    Mu                    = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} LayerSizes
    ##
    ## Sizes of fully connected layers
    ##
    ## A positive integer vector specifying the sizes of the fully connected
    ## layers in the neural network model. The i-th element of
    ## @qcode{LayerSizes} is the number of outputs in the i-th fully connected
    ## layer of the neural network model. @qcode{LayerSizes} does not include
    ## the size of the final fully connected layer. This layer always has K
    ## outputs, where K is the number of classes in Y. This property is
    ## read-only.
    ##
    ## @end deftp
    LayerSizes            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} Activations
    ##
    ## Activation functions for hidden layers
    ##
    ## A character vector or cell array of character vectors specifying the
    ## activation functions used in the hidden layers of the neural network.
    ## Supported activation functions include: @qcode{'linear'},
    ## @qcode{'sigmoid'}, @qcode{'relu'}, @qcode{'tanh'}, @qcode{'softmax'},
    ## @qcode{'lrelu'}, @qcode{'prelu'}, @qcode{'elu'}, and @qcode{'gelu'}.
    ## This property is read-only.
    ##
    ## @end deftp
    Activations           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} OutputLayerActivation
    ##
    ## Activation function for output layer
    ##
    ## A character vector specifying the activation function of the output layer
    ## of the neural network.  Supported activation functions are the same as
    ## for the @qcode{Activations} property.  The default, @qcode{softmax},
    ## reports a probability over the classes; the network is then trained
    ## against cross entropy rather than the mean squared error.  This
    ## property is read-only.
    ##
    ## @end deftp
    OutputLayerActivation = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} LearningRate
    ##
    ## Learning rate for gradient descent
    ##
    ## A positive scalar value defining the learning rate used by the gradient
    ## descent algorithm during training.  This property is read-only.
    ##
    ## @end deftp
    LearningRate          = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} IterationLimit
    ##
    ## Maximum number of training iterations
    ##
    ## A positive integer value defining the maximum number of epochs for
    ## training the model.  This property is read-only.
    ##
    ## @end deftp
    IterationLimit        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} ModelParameters
    ##
    ## Neural network model parameters
    ##
    ## A structure containing the parameters used to train the neural network
    ## classifier model, including layer weights and activations as generated by
    ## the @code{fcnntrain} function.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters       = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} ConvergenceInfo
    ##
    ## Training convergence information
    ##
    ## A structure containing convergence information of the neural network
    ## classifier model with the following fields:
    ##
    ## @itemize
    ## @item @qcode{Accuracy} - The prediction accuracy at each iteration
    ## during training
    ## @item @qcode{TrainingLoss} - The loss value recorded at each iteration
    ## during training
    ## @item @qcode{Time} - The cumulative time taken for all iterations in
    ## seconds
    ## @end itemize
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ConvergenceInfo       = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} DisplayInfo
    ##
    ## Display training information flag
    ##
    ## A boolean flag indicating whether to print information during training.
    ## This property is read-only.
    ##
    ## @end deftp
    DisplayInfo           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} Solver
    ##
    ## Solver used for training
    ##
    ## A character vector specifying the solver algorithm used to train the
    ## neural network model.  Currently only @qcode{'Gradient Descent'} is
    ## supported.  This property is read-only.
    ##
    ## @end deftp
    Solver                = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} LayerWeights
    ##
    ## Learned weights of each fully connected layer
    ##
    ## A cell array holding one weight matrix per layer, the output layer
    ## included.  @code{LayerWeights@{i@}} has one row per neuron of layer
    ## @math{i} and one column per input it receives.  This property is
    ## read-only.
    ##
    ## @end deftp
    LayerWeights          = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} LayerBiases
    ##
    ## Learned bias of each fully connected layer
    ##
    ## A cell array holding one column vector per layer, the output layer
    ## included, with one entry per neuron of that layer.  This property is
    ## read-only.
    ##
    ## @end deftp
    LayerBiases           = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} TrainingHistory
    ##
    ## Iteration by iteration record of training
    ##
    ## A table with one row per iteration, holding the iteration number, the
    ## training loss and the training accuracy recorded at it.  This property
    ## is read-only.
    ##
    ## @end deftp
    TrainingHistory       = [];


    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} Prior
    ##
    ## Prior probability of each class
    ##
    ## A numeric vector with one entry per class, in the order of
    ## @code{ClassNames}, summing to one.  It defaults to the relative
    ## frequency of each class in the training data.  Change it with
    ## @code{setPrior}.
    ##
    ## Specified as a row vector with one entry per class, in the order of
    ## @qcode{ClassNames}, and rescaled to sum to one.  It may be given as
    ## @qcode{'empirical'}, @qcode{'uniform'}, a numeric vector, or a
    ## structure with @qcode{ClassNames} and @qcode{ClassProbs} fields, which
    ## assigns each probability by class name rather than by position.
    ##
    ## @end deftp
    Prior                 = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector with one entry per training observation.  It
    ## defaults to a uniform weight for every observation.  This property is
    ## read-only.
    ##
    ## Each class carries its prior spread evenly over its own observations,
    ## so an observation of a class weighs @qcode{Prior} for that class
    ## divided by the number of observations it holds.
    ##
    ## @end deftp
    W                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} CategoricalPredictors
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
    ## @deftp {ClassificationNeuralNetwork} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the model expanded them
    ##
    ## A cell array of character vectors.  It matches @code{PredictorNames}
    ## unless a categorical predictor was expanded into indicator variables.
    ## This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};
  endproperties

  ## Properties a user may set after the model is built.  Each one is
  ## validated by its set method below.
  properties (GetAccess = public, SetAccess = public)
    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} Cost
    ##
    ## Cost of misclassification
    ##
    ## A numeric matrix with one row and one column per class, where
    ## @code{Cost(i,j)} is the cost of classifying an observation of class
    ## @math{i} as class @math{j}.  The default has zeros on the diagonal and
    ## ones elsewhere.  Change it with @code{setCost}.
    ##
    ## @end deftp
    Cost                  = [];
    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} ScoreTransform
    ##
    ## Transformation function for classification scores
    ##
    ## Specified as a function handle for transforming the classification
    ## scores.  Add or change the @qcode{ScoreTransform} property using dot
    ## notation as in:
    ##
    ## @itemize
    ## @item @qcode{@var{obj}.ScoreTransform = 'function_name'}
    ## @item @qcode{@var{obj}.ScoreTransform = @@function_handle}
    ## @end itemize
    ##
    ## When specified as a character vector, it can be any of the following
    ## built-in functions.  Nevertheless, the @qcode{ScoreTransform} property
    ## always stores their function handle equivalent.
    ##
    ## @multitable @columnfractions 0.2 0.75
    ## @headitem @var{Value} @tab @var{Description}
    ## @item @qcode{'doublelogit'} @tab @math{1 ./ (1 + exp (-2 * x))}
    ## @item @qcode{'invlogit'} @tab @math{log (x ./ (1 - x))}
    ## @item @qcode{'ismax'} @tab Sets the score for the class with the
    ## largest score to 1, and for all other classes to 0
    ## @item @qcode{'logit'} @tab @math{1 ./ (1 + exp (-x))}
    ## @item @qcode{'none'} @tab @math{x} (no transformation)
    ## @item @qcode{'identity'} @tab @math{x} (no transformation)
    ## @item @qcode{'sign'} @tab
    ## @math{-1 for x < 0, 0 for x = 0, 1 for x >
    ## 0}
    ## @item @qcode{'symmetric'} @tab @math{2 * x - 1}
    ## @item @qcode{'symmetricismax'} @tab Sets the score for the class
    ## with the largest score to 1, and for all other classes to -1
    ## @item @qcode{'symmetriclogit'} @tab @math{2 ./ (1 + exp (-x)) - 1}
    ## @end multitable
    ##
    ## @end deftp
    ScoreTransform        = @(x) x;
  endproperties

  ## Readable by the counterpart class, which copies it, and kept out of
  ## the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    STname = 'none';
  endproperties

  ## Set methods for the properties a user may assign.
  methods

    function this = set.ScoreTransform (this, val)
        name = 'ClassificationNeuralNetwork';
        [this.ScoreTransform, this.STname] = parseScoreTransform (val, ...
                                                                  name);
    endfunction

    function this = set.Cost (this, val)
      gnY = this.ClassNames;
      if (isempty (val))
        this.Cost = cast (! eye (numel (gnY)), 'double');
      else
        if (numel (gnY) != sqrt (numel (val)))
          error (strcat ("ClassificationNeuralNetwork: the number of rows and", ...
                         " columns in 'Cost' must correspond to the", ...
                         " selected classes in Y."));
        endif
        this.Cost = val;
      endif
    endfunction

  endmethods

  methods(Hidden)

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
        str = repmat ({'''%s'''}, 1, numel (this.ClassNames));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.ClassNames{:});
      elseif (ischar (this.ClassNames))
        str = repmat ({'''%s'''}, 1, rows (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, cellstr (this.ClassNames){:});
      else # single, double, logical
        str = repmat ({'%d'}, 1, numel (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, this.ClassNames);
      endif
      fprintf ("%+25s: %s\n", 'ClassNames', str);
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.STname);
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
      fprintf ("%+25s: '%s'\n", 'Solver', this.Solver);
    endfunction



  endmethods

  methods(Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {statistics} {@var{obj} =} ClassificationNeuralNetwork (@var{X}, @var{Y})
    ## @deftypefnx {statistics} {@var{obj} =} ClassificationNeuralNetwork (@dots{}, @var{name}, @var{value})
    ##
    ## Create a @qcode{ClassificationNeuralNetwork} class object containing a
    ## neural network classification model.
    ##
    ## @code{@var{obj} = ClassificationNeuralNetwork (@var{X}, @var{Y})} returns
    ## a ClassificationNeuralNetwork object, with @var{X} as the predictor data
    ## and @var{Y} containing the class labels of observations in @var{X}.
    ##
    ## @itemize
    ## @item
    ## @code{X} must be a @math{N*P} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.  @var{X} will be used to train the neural network model.
    ## @item
    ## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}.  @var{Y} can contain any type
    ## of categorical data. @var{Y} must have the same number of rows as
    ## @var{X}.
    ## @end itemize
    ##
    ## @code{@var{obj} = ClassificationNeuralNetwork (@dots{}, @var{name},
    ## @var{value})} returns a ClassificationNeuralNetwork object with
    ## parameters specified by the following @qcode{@var{name}, @var{value}}
    ## paired input arguments:
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'PredictorNames'} @tab A cell array of character
    ## vectors specifying the names of the predictors. The length of this array
    ## must match the number of columns in @var{X}.
    ##
    ## @item @qcode{'ResponseName'} @tab A character vector specifying the
    ## name of the response variable.
    ##
    ## @item @qcode{'ClassNames'} @tab Names of the classes in the class
    ## labels, @var{Y}, used for fitting the neural network model.
    ## @qcode{ClassNames} are of the same type as the class labels in @var{Y}.
    ##
    ## @item @qcode{'ScoreTransform'} @tab A user-defined function handle
    ## or a character vector specifying one of the following builtin functions
    ## specifying the transformation applied to predicted classification scores.
    ## Supported values include @qcode{'doublelogit'}, @qcode{'invlogit'},
    ## @qcode{'ismax'}, @qcode{'logit'}, @qcode{'none'}, @qcode{'identity'},
    ## @qcode{'sign'}, @qcode{'symmetric'}, @qcode{'symmetricismax'}, and
    ## @qcode{'symmetriclogit'}.
    ##
    ## @item @qcode{'Standardize'} @tab A logical scalar specifying whether
    ## to standardize the predictor data.  When @qcode{true}, the predictors are
    ## centered and scaled to have zero mean and unit variance.
    ##
    ## @item @qcode{'LayerSizes'} @tab A positive integer vector specifying
    ## the sizes of the fully connected layers in the neural network.  The
    ## default is 10.
    ##
    ## @item @qcode{'Activations'} @tab A character vector or cell array of
    ## character vectors specifying the activation functions for the hidden
    ## layers.  Supported values include @qcode{'linear'}, @qcode{'sigmoid'},
    ## @qcode{'relu'}, @qcode{'tanh'}, @qcode{'softmax'}, @qcode{'lrelu'},
    ## @qcode{'prelu'}, @qcode{'elu'}, and @qcode{'gelu'}.  The default is
    ## @qcode{'relu'}, whose gradient is one wherever a unit is active and so
    ## does not shrink as it passes back through the layers, where a sigmoid
    ## multiplies it by at most a quarter at every one.
    ##
    ## @item @qcode{'OutputLayerActivation'} @tab A character vector
    ## specifying the activation function for the output layer.  Supported
    ## values are the same as for @qcode{'Activations'}.  The default is
    ## @qcode{'softmax'}, which makes the scores a probability over the
    ## classes and trains the network against cross entropy; any other value
    ## trains it against the mean squared error.
    ##
    ## @item @qcode{'LearningRate'} @tab A positive scalar specifying the
    ## learning rate for gradient descent.  The default is 0.003.  A larger
    ## rate can drive every unit of a hidden layer negative, after which a
    ## rectifier passes no gradient and the network stops training.
    ##
    ## @item @qcode{'IterationLimit'} @tab A positive integer specifying
    ## the maximum number of training iterations.  The default is 1000.
    ##
    ## @item @qcode{'DisplayInfo'} @tab A logical scalar specifying whether
    ## to display training information.  The default is @qcode{false}.
    ## @end multitable
    ##
    ## @seealso{fitcnet}
    ## @end deftypefn
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
      Activations             = 'relu';
      OutputLayerActivation   = 'softmax';
      LearningRate            = 0.003;
      IterationLimit          = 1000;
      DisplayInfo             = false;
      this.Solver = 'Gradient Descent';

      ## Supported activation functions
      acList = {'linear', 'none', 'sigmoid', 'relu', 'tanh', 'softmax', ...
                          'lrelu', 'prelu', 'elu', 'gelu'};
      ## Parse extra parameters
      Prior = [];
      Cost  = [];
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'standardize'
            Standardize = varargin{2};
            if (! (Standardize == true || Standardize == false))
              error (strcat ("ClassificationNeuralNetwork:", ...
                             " 'Standardize' must be either true or false."));
            endif

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("ClassificationNeuralNetwork: 'PredictorNames'", ...
                             " must be supplied as a cellstring array."));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat ("ClassificationNeuralNetwork: 'PredictorNames'", ...
                             " must have the same number of columns as X."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("ClassificationNeuralNetwork: 'ResponseName'", ...
                             " must be a character vector."));
            endif

          case 'classnames'
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
                                   ClassNames, 'UniformOutput', false))))
                error (strcat ("ClassificationNeuralNetwork: not all", ...
                               " 'ClassNames' are present in Y."));
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, 'UniformOutput', false))))
                error (strcat ("ClassificationNeuralNetwork: not all", ...
                               " 'ClassNames' are present in Y."));
              endif
            endif

          case 'scoretransform'
            name = 'ClassificationNeuralNetwork';
            [this.ScoreTransform, this.STname] = parseScoreTransform ...
                                                 (varargin{2}, name);

          case 'prior'
            Prior = varargin{2};

          case 'cost'
            Cost = varargin{2};

          case 'layersizes'
            LayerSizes = varargin{2};
            if (! (isnumeric (LayerSizes) && isvector (LayerSizes)
              && all (LayerSizes > 0) && all (mod (LayerSizes, 1) == 0)))
              error (strcat ("ClassificationNeuralNetwork: 'LayerSizes'", ...
                             " must be a positive integer vector."));
            endif

          case 'learningrate'
            LearningRate = varargin{2};
            if (! (isnumeric (LearningRate) && isscalar (LearningRate) &&
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
                                   Activations, 'UniformOutput', false))))
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
            if (! (isnumeric (IterationLimit) && isscalar (IterationLimit)
              && (IterationLimit > 0) && mod (IterationLimit, 1) == 0))
              error (strcat ("ClassificationNeuralNetwork:", ...
                             " 'IterationLimit' must be a positive integer."));
            endif

          case 'displayinfo'
            DisplayInfo = varargin{2};
            if (! (DisplayInfo == true || DisplayInfo == false))
              error (strcat ("ClassificationNeuralNetwork: 'DisplayInfo'", ...
                             " must be either true or false."));
            endif

          otherwise
            error (strcat ("ClassificationNeuralNetwork: invalid",...
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

      ## An observation is dropped only when its response is missing.  A row
      ## whose predictors hold missing values is kept and reported as used,
      ## while the fit below draws on the complete observations alone.
      RowsUsed  = ! isnan (gY);
      Yret      = Y(RowsUsed);
      Xret      = X(RowsUsed, :);
      this.X    = Xret;
      this.Y    = Yret;
      cobs      = ! any (isnan (Xret), 2);
      Y         = Yret(cobs);
      X         = Xret(cobs, :);

      ## Renew groups in Y over the retained observations, so a class held
      ## only by a row with missing predictors is still a class of the model
      [this.ClassNames, gnY, gret] = unique (Yret);
      gY = gret(cobs);

      ## Check X contains valid data
      if (! (isnumeric (X) && isfinite (X)))
        error ("ClassificationNeuralNetwork: invalid values in X.");
      endif

      ## Assign the number of observations and their corresponding indices
      ## on the original data, which will be used for training the model,
      ## to the ClassificationNeuralNetwork object
      this.NumObservations = rows (this.X);
      ## RowsUsed is left empty when every observation was used, as in MATLAB
      if (all (RowsUsed))
        this.RowsUsed = [];
      else
        this.RowsUsed = RowsUsed;
      endif

      ## Cost, Prior and the observation weights.  Cost defaults to zero on
      ## the diagonal and one elsewhere, Prior to the class frequencies of the
      ## training data, and each class spreads its prior over its own
      ## observations.  Prior cannot be changed once the model is built, so
      ## the option here is the only way to weigh the classes differently.
      nclasses = numel (this.ClassNames);
      if (isempty (Cost))
        this.Cost = ones (nclasses) - eye (nclasses);
      else
        if (! (isnumeric (Cost) && issquare (Cost)
               && rows (Cost) == nclasses))
          error (strcat ("ClassificationNeuralNetwork: 'Cost' must be a", ...
                         " numeric square matrix with one row and column", ...
                         " per class."));
        endif
        this.Cost = Cost;
      endif
      if (isstruct (Prior))
        Prior = priorFromStruct (Prior, this.ClassNames, ...
                                 'ClassificationNeuralNetwork');
      endif
      if (isempty (Prior) || (ischar (Prior) && strcmpi (Prior, 'empirical')))
        this.Prior = accumarray (gY(:), 1, [nclasses, 1])' / numel (gY);
      elseif (ischar (Prior) && strcmpi (Prior, 'uniform'))
        this.Prior = ones (1, nclasses) / nclasses;
      elseif (isnumeric (Prior) && isreal (Prior) && isvector (Prior)
              && numel (Prior) == nclasses && all (Prior >= 0)
              && sum (Prior) > 0)
        this.Prior = Prior(:)' / sum (Prior);
      else
        error (strcat ("ClassificationNeuralNetwork: 'Prior' must be", ...
                       " 'empirical', 'uniform', or a non-negative numeric", ...
                       " vector with one entry per class."));
      endif
      this.W = priorWeights (this.Prior, gY, this.NumObservations);

      ## No predictor is treated as categorical, so the expanded names are
      ## the predictor names themselves.
      this.CategoricalPredictors = [];

      ## Handle Standardize flag
      if (Standardize)
        this.Standardize = true;
        ## Mu and Sigma weight the complete observations so that each class
        ## keeps the share of the observation weight it carried before any row
        ## was set aside, which is what MATLAB reports.
        sw = zeros (rows (X), 1);
        for k = 1:numel (gnY)
          ck = (gY == k);
          if (any (ck))
            sw(ck) = (sum (gret == k) / numel (gret)) / sum (ck);
          endif
        endfor
        sw = sw / sum (sw);
        this.Mu = sum (sw .* X, 1);
        Zs = X - this.Mu;
        this.Sigma = sqrt (sum (sw .* Zs .^ 2, 1) / (1 - sum (sw .^ 2)));
        this.Sigma(this.Sigma == 0) = 1;  # predictor is constant
        ## Train on the scale the model predicts on: predict, resubPredict
        ## and loss all standardize their input from Mu and Sigma, so the
        ## training data must be standardized here as well.
        X = (X - this.Mu) ./ this.Sigma;
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
      ## A softmax output reports a probability over the classes, and the
      ## loss that belongs with it is cross entropy: paired that way the two
      ## gradients compose to y - t.  Any other output layer keeps the mean
      ## squared error.
      LossFunction = double (strcmp (OutputLayerActivation, 'softmax'));
      Mdl = fcnntrain (X, gY, LayerSizes, ActivationCodes, NumThreads, ...
                       Alpha, LearningRate, IterationLimit, DisplayInfo, ...
                       LossFunction);

      ## Store training time, Iterations, and Loss
      ConvergenceInfo.Time = toc (cnn_timer_);
      ConvergenceInfo.Accuracy = Mdl.Accuracy;
      ConvergenceInfo.TrainingLoss = Mdl.Loss;

      ## Remove redundant fields
      Mdl = rmfield (Mdl, 'Accuracy');
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
      this.TrainingHistory = table (iter, ...
                                    ConvergenceInfo.TrainingLoss(:), ...
                                    ConvergenceInfo.Accuracy(:), ...
                                    'VariableNames', ...
                                    {'Iteration', 'TrainingLoss', ...
                                     'TrainingAccuracy'});

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationNeuralNetwork} {[@var{label}, @var{score}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the neural network
    ## classification model from a ClassificationNeuralNetwork object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## predictor data in @code{obj.X} and corresponding labels, @code{obj.Y},
    ## stored in the ClassificationNeuralNetwork model, @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationNeuralNetwork} class object.
    ## @item
    ## @var{XC} must be an @math{M*P} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the neural network
    ## model in @var{obj}.
    ## @end itemize
    ##
    ## @code{[@var{label}, @var{score}] = predict (@var{obj}, @var{XC})} also
    ## returns @var{score}, which contains the predicted class scores or
    ## posterior probabilities for each instance of the corresponding unique
    ## classes.
    ##
    ## The @var{score} matrix contains the classification scores for each class.
    ## For each observation in @var{XC}, the predicted class label is the one
    ## with the highest score among all classes.  If the @qcode{ScoreTransform}
    ## property is set to a transformation function, the scores are transformed
    ## accordingly before being returned.
    ##
    ## @seealso{ClassificationNeuralNetwork, fitcnet}
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

      ## Apply ScoreTransform
      scores = this.ScoreTransform (scores);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{label} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationNeuralNetwork} {[@var{label}, @var{score}] =} resubPredict (@var{obj})
    ##
    ## Classify the training data using the trained neural network
    ## classification object.
    ##
    ## @code{@var{label} = resubPredict (@var{obj})} returns the vector of
    ## labels predicted for the corresponding instances in the training data,
    ## using the predictor data in @code{obj.X} and corresponding labels,
    ## @code{obj.Y}, stored in the neural network classification model,
    ## @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationNeuralNetwork} class object.
    ## @end itemize
    ##
    ## @code{[@var{label}, @var{score}] = resubPredict (@var{obj})} also
    ## returns @var{score}, which contains the predicted class scores or
    ## posterior probabilities for each instance of the corresponding unique
    ## classes.
    ##
    ## @seealso{ClassificationNeuralNetwork, fitcnet}
    ## @end deftypefn
    function [labels, scores] = resubPredict (this)

      ## Get used rows
      X = this.X;

      ## Standardize (if necessary)
      if (this.Standardize)
        X = (X - this.Mu) ./ this.Sigma;
      endif

      ## Predict labels from existing data
      NumThreads = nproc ();
      [labels, scores] = fcnnpredict (this.ModelParameters, X, NumThreads);

      # Get class labels
      labels = this.ClassNames(labels);

      ## Apply ScoreTransform
      scores = this.ScoreTransform (scores);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin of a neural network classifier.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns a column
    ## vector holding, for each row of @var{X}, the score the model gives its
    ## true class in @var{Y} less the largest score it gives any other class.
    ## A positive margin means the observation is classified correctly, and
    ## the larger it is the more confidently so.
    ##
    ## @seealso{ClassificationNeuralNetwork, edge, loss, predict}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationNeuralNetwork.margin: too few input arguments.");
      endif

      [X, Y] = checkXY_ (this, X, Y, "margin");

      [~, scores] = predict (this, X);
      classes = this.ClassNames;
      m = zeros (rows (X), 1);
      for i = 1:rows (X)
        idx = find (ismember (classes, Y(i)));
        if (isempty (idx))
          m(i) = NaN;
          continue;
        endif
        true_score = scores(i, idx);
        scores(i, idx) = -Inf;
        m(i) = true_score - max (scores(i,:));
        scores(i, idx) = true_score;
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationNeuralNetwork} {@var{e} =} edge (@dots{}, @qcode{"Weights"}, @var{w})
    ##
    ## Classification edge of a neural network classifier.
    ##
    ## @code{@var{e} = edge (@var{obj}, @var{X}, @var{Y})} returns the mean of
    ## the classification margins over the rows of @var{X}.
    ##
    ## @code{@var{e} = edge (@dots{}, @qcode{"Weights"}, @var{w})} takes the
    ## weighted mean instead, @var{w} holding one weight per row of @var{X}.
    ## The weights are normalised to sum to one before they are applied.
    ##
    ## @seealso{ClassificationNeuralNetwork, margin, loss, predict}
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationNeuralNetwork.edge: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationNeuralNetwork.edge: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, "edge");
      W = getWeights_ (this, varargin, rows (X), "edge");

      m = margin (this, X, Y);
      e = sum (W(:) .* m) / sum (W);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationNeuralNetwork} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of a neural network classifier.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## proportion of the rows of @var{X} the model misclassifies against the
    ## true labels @var{Y}.
    ##
    ## @code{@var{L} = loss (@dots{}, @var{name}, @var{value})} accepts the
    ## following name-value pairs:
    ##
    ## @itemize
    ## @item
    ## @qcode{"LossFun"} selects the loss.  Supported values are
    ## @qcode{"mincost"}, the default, @qcode{"binodeviance"},
    ## @qcode{"classifcost"}, @qcode{"classiferror"}, @qcode{"crossentropy"},
    ## @qcode{"exponential"}, @qcode{"hinge"}, @qcode{"logit"} and
    ## @qcode{"quadratic"}.  @qcode{"mincost"} assigns each observation to
    ## the class of least expected cost and charges what that assignment
    ## costs, so it reads the scores as a posterior; @qcode{"classifcost"}
    ## charges what the model's own prediction costs.  @qcode{"crossentropy"}
    ## is defined for a network only.  Note that the default differs from the
    ## other classifiers in this package, which default to
    ## @qcode{"classiferror"}, and follows MATLAB's for this class.
    ##
    ## @item
    ## @qcode{"Weights"} holds one weight per row of @var{X}, normalised to
    ## sum to one before it is applied.
    ## @end itemize
    ##
    ## @seealso{ClassificationNeuralNetwork, margin, edge, predict}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationNeuralNetwork.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationNeuralNetwork.loss: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, "loss");

      ## Parse optional arguments
      LossFun = 'mincost';
      lossnames = {'binodeviance', 'classifcost', 'classiferror', ...
                   'crossentropy', 'exponential', 'hinge', 'logit', ...
                   'mincost', 'quadratic'};
      args = varargin;
      keep = true (1, numel (args));
      for i = 1:2:numel (args)
        if (strcmpi (args{i}, 'lossfun'))
          LossFun = args{i+1};
          if (! (ischar (LossFun) && isrow (LossFun)))
            error (strcat ("ClassificationNeuralNetwork.loss: 'LossFun'", ...
                           " must be a character vector."));
          endif
          LossFun = tolower (LossFun);
          if (! any (strcmpi (LossFun, lossnames)))
            error (strcat ("ClassificationNeuralNetwork.loss: unsupported", ...
                           " Loss function."));
          endif
          keep(i:i+1) = false;
        endif
      endfor
      W = getWeights_ (this, args(keep), rows (X), "loss");
      W = W(:) / sum (W);

      [label, scores] = predict (this, X);
      classes = this.ClassNames;
      K = numel (classes);

      ## Membership of the true class, as a +1/-1 indicator per class
      Yind = zeros (rows (X), K);
      for i = 1:rows (X)
        idx = find (ismember (classes, Y(i)));
        if (isempty (idx))
          L = NaN;
          return;
        endif
        Yind(i, idx) = 1;
      endfor

      ## The scalar score of the true class of each observation
      mj = sum (scores .* Yind, 2);

      switch (LossFun)
        case 'classiferror'
          wrong = zeros (rows (X), 1);
          for i = 1:rows (X)
            wrong(i) = ! isequal (label(i), Y(i));
          endfor
          L = sum (W .* wrong);
        case 'binodeviance'
          L = sum (W .* log (1 + exp (-2 * mj)));
        case 'hinge'
          L = sum (W .* max (0, 1 - mj));
        case 'exponential'
          L = sum (W .* exp (-mj));
        case 'logit'
          L = sum (W .* log (1 + exp (-mj)));
        case 'quadratic'
          L = sum (W .* (1 - mj) .^ 2);
        case 'mincost'
          ## Each observation is assigned to the class of least expected
          ## cost, and charged what that assignment actually costs given its
          ## true class.
          L = 0;
          for i = 1:rows (X)
            [~, k] = min (scores(i,:) * this.Cost);
            true_idx = find (ismember (classes, Y(i)));
            L = L + W(i) * this.Cost(true_idx, k);
          endfor
        case 'classifcost'
          ## What the model's own prediction costs, given the true class
          L = 0;
          for i = 1:rows (X)
            true_idx = find (ismember (classes, Y(i)));
            pred_idx = find (ismember (classes, label(i)));
            L = L + W(i) * this.Cost(true_idx, pred_idx);
          endfor
        case 'crossentropy'
          ## Defined for a network only, whose scores are a posterior.  The
          ## weights are rescaled to sum to the number of observations, as
          ## MATLAB documents, and the sum is taken over classes as well.
          Wn = W * rows (X);
          L = -sum (Wn .* log (max (mj, realmin))) / (K * rows (X));
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{m} =} resubMargin (@var{obj})
    ##
    ## Classification margin of a neural network classifier on its training
    ## data.
    ##
    ## @code{@var{m} = resubMargin (@var{obj})} is @code{margin} applied to
    ## the observations the model was fitted on.
    ##
    ## @seealso{ClassificationNeuralNetwork, margin}
    ## @end deftypefn
    function m = resubMargin (this)
      if (nargin != 1)
        error (strcat ("ClassificationNeuralNetwork.resubMargin:", ...
                       " too many input arguments."));
      endif
      m = margin (this, this.X, this.Y);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{e} =} resubEdge (@var{obj})
    ##
    ## Classification edge of a neural network classifier on its training
    ## data.
    ##
    ## @code{@var{e} = resubEdge (@var{obj})} is @code{edge} applied to the
    ## observations the model was fitted on, weighted by @code{obj.W}.
    ##
    ## @seealso{ClassificationNeuralNetwork, edge}
    ## @end deftypefn
    function e = resubEdge (this)
      if (nargin != 1)
        error (strcat ("ClassificationNeuralNetwork.resubEdge:", ...
                       " too many input arguments."));
      endif
      e = edge (this, this.X, this.Y, 'Weights', this.W);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {ClassificationNeuralNetwork} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of a neural network classifier on its training
    ## data.
    ##
    ## @code{@var{L} = resubLoss (@var{obj})} is @code{loss} applied to the
    ## observations the model was fitted on, weighted by @code{obj.W}.  It
    ## takes the same @qcode{"LossFun"} name-value pair.
    ##
    ## @seealso{ClassificationNeuralNetwork, loss}
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationNeuralNetwork.resubLoss: Name-Value", ...
                       " arguments must be in pairs."));
      endif
      L = loss (this, this.X, this.Y, varargin{:}, 'Weights', this.W);
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
    ## @multitable @columnfractions 0.28 0.7
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'KFold'} @tab Specify the number of folds to use in
    ## k-fold cross-validation.  @code{"KFold", @var{k}}, where @var{k} is an
    ## integer greater than 1.
    ##
    ## @item @qcode{'Holdout'} @tab Specify the fraction of the data to
    ## hold out for testing.  @code{"Holdout", @var{p}}, where @var{p} is a
    ## scalar in the range @math{(0,1)}.
    ##
    ## @item @qcode{'Leaveout'} @tab Specify whether to perform
    ## leave-one-out cross-validation.  @code{"Leaveout", @var{Value}}, where
    ## @var{Value} is 'on' or 'off'.
    ##
    ## @item @qcode{'CVPartition'} @tab Specify a @qcode{cvpartition}
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
      if (this.NumObservations < 10)
        numFolds  = this.NumObservations;
      else
        numFolds  = 10;
      endif
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
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("ClassificationNeuralNetwork.crossval:", ...
                             " 'CVPartition' must be a 'cvpartition' object."));
            endif

          otherwise
            error (strcat ("ClassificationNeuralNetwork.crossval: invalid",...
                           " parameter name in optional paired arguments."));
          endswitch
        varargin(1:2) = [];
      endwhile

      ## Determine the cross-validation method to use.  The partition covers
      ## the observations actually trained on: a row dropped for a missing
      ## value is not one the folds can use, and including it would leave the
      ## partition, the stored data and NumObservations disagreeing.  The
      ## response is passed rather than a count so the folds stay stratified.
      Yused = this.Y;
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (Yused, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (this.NumObservations, 'LeaveOut');
      else
        partition = cvpartition (Yused, 'KFold', numFolds);
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
      ## Create a compact model
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
    ## Use
    ## @code{loadmodel} in order to load a classification object into Octave's
    ## workspace.
    ##
    ## @seealso{loadmodel, fitcnet, ClassificationNeuralNetwork}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error ("ClassificationNeuralNetwork.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error ("ClassificationNeuralNetwork.savemodel: FNAME must be a character vector.");
      endif
      ## Generate variable for class name
      classdef_name = 'ClassificationNeuralNetwork';

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
      DisplayInfo             = this.DisplayInfo;
      Solver                  = this.Solver;
      LayerWeights            = this.LayerWeights;
      LayerBiases             = this.LayerBiases;
      Cost                    = this.Cost;
      Prior                   = this.Prior;
      W                       = this.W;
      CategoricalPredictors   = this.CategoricalPredictors;
      ExpandedPredictorNames  = this.ExpandedPredictorNames;
      STname                  = this.STname;

      ## TrainingHistory is a table, and Octave cannot save a classdef object
      ## to a binary file, so it is left out here and rebuilt on loading from
      ## ConvergenceInfo, which holds the same numbers as plain vectors.

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'NumObservations', ...
            'RowsUsed', 'NumPredictors', 'PredictorNames', 'ResponseName', ...
            'ClassNames', 'ScoreTransform', 'Standardize', 'Sigma', 'Mu', ...
            'LayerSizes', 'Activations', 'OutputLayerActivation', ...
            'LearningRate', 'IterationLimit', 'Solver', 'ModelParameters', ...
            'ConvergenceInfo', 'DisplayInfo', 'LayerWeights', 'LayerBiases', ...
            'Cost', 'Prior', 'W', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'STname');
    endfunction

  endmethods

  methods (Access = private)

    ## Shared validation for the assessment methods, so each reports under
    ## its own name.
    function [X, Y] = checkXY_ (this, X, Y, caller)
      if (isempty (X))
        error ("ClassificationNeuralNetwork.%s: X is empty.", caller);
      elseif (columns (this.X) != columns (X))
        error (strcat ("ClassificationNeuralNetwork.%s: X must have the", ...
                       " same number of predictors as the trained model."), ...
               caller);
      endif
      if (isempty (Y))
        error ("ClassificationNeuralNetwork.%s: Y is empty.", caller);
      elseif (rows (X) != rows (Y))
        error (strcat ("ClassificationNeuralNetwork.%s: Y must have the", ...
                       " same number of rows as X."), caller);
      endif
    endfunction

    ## Pull a "Weights" pair out of the optional arguments, defaulting to a
    ## uniform weight, and reject any other name.
    function W = getWeights_ (this, args, n, caller)
      W = ones (n, 1);
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("ClassificationNeuralNetwork.%s: parameter name", ...
                         " must be a character vector."), caller);
        endif
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
          if (! (isnumeric (W) && isvector (W)))
            error (strcat ("ClassificationNeuralNetwork.%s: 'Weights'", ...
                           " must be a numeric vector."), caller);
          endif
          if (numel (W) != n)
            error (strcat ("ClassificationNeuralNetwork.%s: size of", ...
                           " 'Weights' must equal the number of", ...
                           " rows in X."), caller);
          endif
        else
          error (strcat ("ClassificationNeuralNetwork.%s: invalid", ...
                         " parameter name in optional paired", ...
                         " arguments."), caller);
        endif
      endfor
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{obj} =} setCost (@var{obj}, @var{cost})
    ##
    ## Set the misclassification cost of a neural network classifier.
    ##
    ## @code{@var{obj} = setCost (@var{obj}, @var{cost})} replaces the
    ## @code{Cost} property.  @var{cost} must be square with one row and one
    ## column per class, in the order of @code{ClassNames}.  An empty
    ## @var{cost} restores the default, zero on the diagonal and one
    ## elsewhere.
    ##
    ## @seealso{ClassificationNeuralNetwork, setPrior, loss}
    ## @end deftypefn
    function this = setCost (this, Cost)
      if (nargin != 2)
        error ("ClassificationNeuralNetwork.setCost: too few input arguments.");
      endif
      K = numel (this.ClassNames);
      if (isempty (Cost))
        this.Cost = ones (K) - eye (K);
        return;
      endif
      if (! (isnumeric (Cost) && isreal (Cost)
             && isequal (size (Cost), [K, K])))
        error (strcat ("ClassificationNeuralNetwork.setCost: 'Cost' must", ...
                       " be a square numeric matrix with one row and one", ...
                       " column per class."));
      endif
      this.Cost = Cost;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {@var{obj} =} setPrior (@var{obj}, @var{prior})
    ##
    ## Set the class prior probabilities of a neural network classifier.
    ##
    ## @code{@var{obj} = setPrior (@var{obj}, @var{prior})} replaces the
    ## @code{Prior} property.  @var{prior} may be @qcode{"empirical"}, the
    ## class frequencies of the training data and the default,
    ## @qcode{"uniform"}, or a numeric vector with one entry per class in the
    ## order of @code{ClassNames}.  A numeric vector is normalised to sum to
    ## one.
    ##
    ## @seealso{ClassificationNeuralNetwork, setCost, loss}
    ## @end deftypefn
    function this = setPrior (this, Prior)
      if (nargin != 2)
        error (strcat ("ClassificationNeuralNetwork.setPrior:", ...
                       " too few input arguments."));
      endif
      K = numel (this.ClassNames);
      if (ischar (Prior) && strcmpi (Prior, 'uniform'))
        this.Prior = ones (1, K) / K;
      elseif (isempty (Prior)
              || (ischar (Prior) && strcmpi (Prior, 'empirical')))
        [~, ~, gY] = unique (this.Y);
        this.Prior = accumarray (gY(:), 1, [K, 1])' / numel (gY);
      elseif (isnumeric (Prior) && isreal (Prior) && isvector (Prior)
              && numel (Prior) == K && all (Prior >= 0) && sum (Prior) > 0)
        this.Prior = Prior(:)' / sum (Prior);
      else
        error (strcat ("ClassificationNeuralNetwork.setPrior: 'Prior' must", ...
                       " be 'empirical', 'uniform', or a non-negative", ...
                       " numeric vector with one entry per class."));
      endif
    endfunction

  endmethods

  methods(Static, Hidden)

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

      ## Rebuild the TrainingHistory table, which savemodel cannot write out
      iter = (1:numel (mdl.ConvergenceInfo.TrainingLoss))';
      mdl.TrainingHistory = table (iter, ...
                                   mdl.ConvergenceInfo.TrainingLoss(:), ...
                                   mdl.ConvergenceInfo.Accuracy(:), ...
                                   'VariableNames', ...
                                   {'Iteration', 'TrainingLoss', ...
                                    'TrainingAccuracy'});
    endfunction

  endmethods

endclassdef

function numCode = activationCode (strCode)
  switch (strCode)
    ## 'none' is MATLAB's name for the identity, which fitcnet has always
    ## documented as available; 'linear' is this package's older spelling of
    ## the same map.
    case {'linear', 'none'}
      numCode = 0;
    case 'sigmoid'
      numCode = 1;
    case 'relu'
      numCode = 2;
    case 'tanh'
      numCode = 3;
    case 'softmax'
      numCode = 4;
    case {'lrelu', 'prelu'}
      numCode = 5;
    case 'elu'
      numCode = 6;
    case 'gelu'
      numCode = 7;
    otherwise
      error (strcat ("ClassificationNeuralNetwork: misspelling or unsupported", ...
                     " 'Activation' function: '%s'."), strCode);
  endswitch
endfunction

## Test input validation for constructor
%!error<ClassificationNeuralNetwork: too few input arguments.> ...
%! ClassificationNeuralNetwork ()
%!error<ClassificationNeuralNetwork: too few input arguments.> ...
%! ClassificationNeuralNetwork (ones (10,2))
%!error<ClassificationNeuralNetwork: number of rows in X and Y must be equal.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (5,1))
%!error<ClassificationNeuralNetwork: 'Standardize' must be either true or false.> ...
%! ClassificationNeuralNetwork (ones (5,3), ones (5,1), 'standardize', 'a')
%!error<ClassificationNeuralNetwork: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), 'PredictorNames', ['A'])
%!error<ClassificationNeuralNetwork: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), 'PredictorNames', 'A')
%!error<ClassificationNeuralNetwork: 'PredictorNames' must have the same number of columns as X.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), 'PredictorNames', {'A', 'B', 'C'})
%!error<ClassificationNeuralNetwork: 'ResponseName' must be a character vector.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), 'ResponseName', {'Y'})
%!error<ClassificationNeuralNetwork: 'ResponseName' must be a character vector.> ...
%! ClassificationNeuralNetwork (ones (5,2), ones (5,1), 'ResponseName', 1)
%!error<ClassificationNeuralNetwork: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'ClassNames', @(x)x)
%!error<ClassificationNeuralNetwork: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'ClassNames', {1})
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'ClassNames', [1, 2])
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones (5,2), ['a';'b';'a';'a';'b'], 'ClassNames', ['a';'c'])
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones (5,2), {'a';'b';'a';'a';'b'}, 'ClassNames', {'a','c'})
%!error<ClassificationNeuralNetwork: not all 'ClassNames' are present in Y.> ...
%! ClassificationNeuralNetwork (ones (10,2), logical (ones (10,1)), 'ClassNames', [true, false])
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'LayerSizes', -1)
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'LayerSizes', 0.5)
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'LayerSizes', [1,-2])
%!error<ClassificationNeuralNetwork: 'LayerSizes' must be a positive integer vector.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'LayerSizes', [10,20,30.5])
%!error<ClassificationNeuralNetwork: 'LearningRate' must be a positive scalar.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'LearningRate', -0.1)
%!error<ClassificationNeuralNetwork: 'LearningRate' must be a positive scalar.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'LearningRate', [0.1, 0.01])
%!error<ClassificationNeuralNetwork: 'LearningRate' must be a positive scalar.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'LearningRate', 'a')
%!error<ClassificationNeuralNetwork: 'Activations' must be a character vector or a cellstring vector.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'Activations', 123)
%!error<ClassificationNeuralNetwork: unsupported 'Activation' function.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'Activations', 'unsupported_type')
%!error<ClassificationNeuralNetwork: unsupported 'Activation' functions.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'LayerSizes', [10, 5], ...
%! 'Activations', {'sigmoid', 'unsupported_type'})
%!error<ClassificationNeuralNetwork: 'Activations' vector does not match the number of layers.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'Activations', {'sigmoid', 'relu', 'softmax'})
%!error<ClassificationNeuralNetwork: 'OutputLayerActivation' must be a character vector.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'OutputLayerActivation', 123)
%!error<ClassificationNeuralNetwork: unsupported 'OutputLayerActivation' function.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'OutputLayerActivation', 'unsupported_type')
%!error<ClassificationNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'IterationLimit', -1)
%!error<ClassificationNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'IterationLimit', 0.5)
%!error<ClassificationNeuralNetwork: 'IterationLimit' must be a positive integer.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'IterationLimit', [1,2])
%!error<ClassificationNeuralNetwork: 'ScoreTransform' must be a character vector or a function handle.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'ScoreTransform', [1,2])
%!error<ClassificationNeuralNetwork: unrecognized 'ScoreTransform' function.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'ScoreTransform', 'unsupported_type')
%!error<ClassificationNeuralNetwork: invalid parameter name in optional pair arguments.> ...
%! ClassificationNeuralNetwork (ones (10,2), ones (10,1), 'some', 'some')
%!error<ClassificationNeuralNetwork: invalid values in X.> ...
%! ClassificationNeuralNetwork ([1;2;3;'a';4], ones (5,1))
%!error<ClassificationNeuralNetwork: invalid values in X.> ...
%! ClassificationNeuralNetwork ([1;2;3;Inf;4], ones (5,1))

## Test input validation for subsasgn method
%!shared x, y, objST, Mdl
%! load fisheriris
%! x = meas;
%! y = grp2idx (species);
%! Mdl = fitcnet (x, y, 'IterationLimit', 100);
%!error<ClassificationNeuralNetwork: unrecognized 'ScoreTransform' function.> ...
%! Mdl.ScoreTransform = 'a';

## Test input validation for predict method
%!error<ClassificationNeuralNetwork.predict: too few input arguments.> ...
%! predict (Mdl)
%!error<ClassificationNeuralNetwork.predict: XC is empty.> ...
%! predict (Mdl, [])
%!error<ClassificationNeuralNetwork.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (Mdl, 1)

## Test output for crossval method
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (Mdl, 'KFold', 5);
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (CVMdl.KFold == 5, true)
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationNeuralNetwork")
%! assert_equal (CVMdl.CrossValidatedModel, "ClassificationNeuralNetwork")
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (Mdl, 'HoldOut', 0.2);
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationNeuralNetwork")
%! assert_equal (CVMdl.CrossValidatedModel, "ClassificationNeuralNetwork")

## Test input validation for crossval method
%!error<ClassificationNeuralNetwork.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (Mdl, 'KFold')
%!error<ClassificationNeuralNetwork.crossval: specify only one of the optional Name-Value paired arguments.> ...
%! crossval (Mdl, 'KFold', 5, 'leaveout', 'on')
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (Mdl, 'KFold', 'a')
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (Mdl, 'KFold', 1)
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (Mdl, 'KFold', -1)
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (Mdl, 'KFold', 11.5)
%!error<ClassificationNeuralNetwork.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (Mdl, 'KFold', [1,2])
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (Mdl, 'Holdout', 'a')
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (Mdl, 'Holdout', 11.5)
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (Mdl, 'Holdout', -1)
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (Mdl, 'Holdout', 0)
%!error<ClassificationNeuralNetwork.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (Mdl, 'Holdout', 1)
%!error<ClassificationNeuralNetwork.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (Mdl, 'Leaveout', 1)
%!error<ClassificationNeuralNetwork.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (Mdl, 'CVPartition', 1)
%!error<ClassificationNeuralNetwork.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (Mdl, 'CVPartition', 'a')
%!error<ClassificationNeuralNetwork.crossval: invalid parameter name in optional paired arguments> ...
%! crossval (Mdl, 'some', 'some')
## A saved and reloaded model carries the trained parameters, not those of
## the placeholder object load_model starts from.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], ...
%!                'LayerSizes', [3, 2], 'IterationLimit', 20);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (Mdl2.LayerWeights, Mdl.LayerWeights);
%! assert_equal (Mdl2.LayerBiases, Mdl.LayerBiases);

## Cost, Prior, W and the predictor bookkeeping survive a save and load.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'IterationLimit', 20);
%! Mdl = setCost (Mdl, [0, 3; 5, 0]);
%! Mdl = setPrior (Mdl, [0.25, 0.75]);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (Mdl2.Cost, [0, 3; 5, 0]);
%! assert_equal (Mdl2.Prior, [0.25, 0.75]);
%! assert_equal (Mdl2.W, Mdl.W);
%! assert_equal (Mdl2.ExpandedPredictorNames, Mdl.ExpandedPredictorNames);
%! assert_equal (Mdl2.CategoricalPredictors, Mdl.CategoricalPredictors);

## TrainingHistory comes back a table, rebuilt from ConvergenceInfo because
## Octave cannot write a classdef object to a binary file.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'IterationLimit', 25);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (istable (Mdl2.TrainingHistory), true);
%! assert_equal (height (Mdl2.TrainingHistory), 25);
%! assert_equal (table2cell (Mdl2.TrainingHistory), ...
%!               table2cell (Mdl.TrainingHistory));

## A reloaded model predicts exactly what the original did.
%!test
%! load fisheriris
%! Mdl = fitcnet (meas, species, 'IterationLimit', 20);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mdl2 = loadmodel (fname);
%! delete (fname);
%! [label, score] = predict (Mdl, meas);
%! [label2, score2] = predict (Mdl2, meas);
%! assert_equal (label2, label);
%! assert_equal (score2, score);

## A non-default ScoreTransform survives a save and load.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'IterationLimit', 20);
%! Mdl.ScoreTransform = 'symmetric';
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (Mdl2.ScoreTransform (0.25), -0.5);

%!error <ClassificationNeuralNetwork.savemodel: too few input arguments.> ...
%! savemodel (ClassificationNeuralNetwork ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]))
%!error <ClassificationNeuralNetwork.savemodel: FNAME must be a character vector.> ...
%! savemodel (ClassificationNeuralNetwork ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), 1)
%!error <ClassificationNeuralNetwork.savemodel: FNAME must be a character vector.> ...
%! savemodel (ClassificationNeuralNetwork ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), ['ab'; 'cd'])

## A trained network returns a posterior, not a fraction of one.  A backward
## pass that settles at half the target leaves every label right and every
## score halved, so the scores and their row sums are what catch it.
%!test
%! rand ('seed', 42);
%! randn ('seed', 42);
%! X = [randn(40, 2) * 0.3 + 3; randn(40, 2) * 0.3 - 3];
%! Y = [ones(40, 1); 2 * ones(40, 1)];
%! Mdl = fitcnet (X, Y, 'LayerSizes', [8, 8], 'IterationLimit', 400);
%! [label, score] = predict (Mdl, [3, 3; -3, -3]);
%! assert_equal (label, [1; 2]);
%! assert_equal (all (abs (sum (score, 2) - 1) < 0.1), true);
%! assert_equal (max (score(1,:)) > 0.8, true);
%! assert_equal (max (score(2,:)) > 0.8, true);

## The defaults follow MATLAB: rectified hidden layers and a softmax output.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]);
%! assert_equal (Mdl.Activations, 'relu');
%! assert_equal (Mdl.OutputLayerActivation, 'softmax');

## Training reports a history that was recorded, one entry per iteration.
%!test
%! rand ('seed', 42);
%! randn ('seed', 42);
%! X = [randn(30, 2) * 0.4 + 2; randn(30, 2) * 0.4 - 2];
%! Y = [ones(30, 1); 2 * ones(30, 1)];
%! Mdl = fitcnet (X, Y, 'IterationLimit', 100);
%! loss = Mdl.ConvergenceInfo.TrainingLoss;
%! acc = Mdl.ConvergenceInfo.Accuracy;
%! assert_equal (numel (loss), 100);
%! assert_equal (numel (acc), 100);
%! assert_equal (any (loss != 0), true);
%! assert_equal (loss(end) < loss(1), true);
%! assert_equal (mean (predict (Mdl, X) == Y) > 0.95, true);

## Every activation trains to a usable posterior on separable data.
%!test
%! rand ('seed', 7);
%! randn ('seed', 7);
%! X = [randn(30, 2) * 0.3 + 3; randn(30, 2) * 0.3 - 3];
%! Y = [ones(30, 1); 2 * ones(30, 1)];
%! names = {'linear', 'sigmoid', 'relu', 'tanh', 'lrelu', 'elu', 'gelu'};
%! for k = 1:numel (names)
%!   Mdl = fitcnet (X, Y, 'LayerSizes', 8, 'Activations', names{k}, ...
%!                  'IterationLimit', 300);
%!   [label, score] = predict (Mdl, [3, 3; -3, -3]);
%!   assert_equal (label, [1; 2]);
%!   assert_equal (all (isfinite (score(:))), true);
%! endfor

## The trained parameters are reachable, one matrix and one bias per layer.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], ...
%!                'LayerSizes', [3, 2], 'IterationLimit', 20);
%! assert_equal (numel (Mdl.LayerWeights), 3);
%! assert_equal (numel (Mdl.LayerBiases), 3);
%! assert_equal (size (Mdl.LayerWeights{1}), [3, 2]);
%! assert_equal (size (Mdl.LayerWeights{2}), [2, 3]);
%! assert_equal (size (Mdl.LayerWeights{3}), [2, 2]);
%! assert_equal (size (Mdl.LayerBiases{1}), [3, 1]);
%! assert_equal (size (Mdl.LayerBiases{3}), [2, 1]);

## Cost, Prior, W and the predictor bookkeeping carry their defaults.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'IterationLimit', 20);
%! assert_equal (Mdl.Cost, [0, 1; 1, 0]);
%! assert_equal (Mdl.Prior, [0.5, 0.5]);
%! assert_equal (size (Mdl.W), [4, 1]);
%! assert_equal (sum (Mdl.W), 1, 1e-12);
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (Mdl.ExpandedPredictorNames, Mdl.PredictorNames);

## TrainingHistory holds one row per iteration.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'IterationLimit', 25);
%! assert_equal (istable (Mdl.TrainingHistory), true);
%! assert_equal (height (Mdl.TrainingHistory), 25);
%! assert_equal (Mdl.TrainingHistory.Properties.VariableNames, ...
%!               {'Iteration', 'TrainingLoss', 'TrainingAccuracy'});
%! assert_equal (Mdl.TrainingHistory.Iteration', 1:25);

## A margin is positive wherever the model is right, and edge averages it.
%!test
%! rand ('seed', 42);
%! randn ('seed', 42);
%! X = [randn(30, 2) * 0.4 + 2; randn(30, 2) * 0.4 - 2];
%! Y = [ones(30, 1); 2 * ones(30, 1)];
%! Mdl = fitcnet (X, Y, 'IterationLimit', 200);
%! m = margin (Mdl, X, Y);
%! assert_equal (size (m), [60, 1]);
%! assert_equal (all (m > 0), true);
%! assert_equal (edge (Mdl, X, Y), mean (m), 1e-12);
%! assert_equal (resubMargin (Mdl), m, 1e-12);
%! assert_equal (resubEdge (Mdl), edge (Mdl, X, Y, 'Weights', Mdl.W), 1e-12);

## Every documented loss function is accepted and finite.
%!test
%! rand ('seed', 42);
%! randn ('seed', 42);
%! X = [randn(30, 2) * 0.4 + 2; randn(30, 2) * 0.4 - 2];
%! Y = [ones(30, 1); 2 * ones(30, 1)];
%! Mdl = fitcnet (X, Y, 'IterationLimit', 200);
%! names = {'binodeviance', 'classifcost', 'classiferror', 'crossentropy', ...
%!          'exponential', 'hinge', 'logit', 'mincost', 'quadratic'};
%! for k = 1:numel (names)
%!   L = loss (Mdl, X, Y, 'LossFun', names{k});
%!   assert_equal (isscalar (L) && isfinite (L) && L >= 0, true);
%! endfor
%! assert_equal (loss (Mdl, X, Y), loss (Mdl, X, Y, 'LossFun', 'mincost'));
%! assert_equal (resubLoss (Mdl), loss (Mdl, X, Y, 'Weights', Mdl.W), 1e-12);

## setPrior and setCost replace the properties they name.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 1; 2], 'IterationLimit', 20);
%! assert_equal (Mdl.Prior, [0.75, 0.25], 1e-12);
%! assert_equal (setPrior (Mdl, 'uniform').Prior, [0.5, 0.5]);
%! assert_equal (setPrior (Mdl, [3, 1]).Prior, [0.75, 0.25], 1e-12);
%! assert_equal (setPrior (Mdl, 'empirical').Prior, Mdl.Prior, 1e-12);
%! assert_equal (setCost (Mdl, [0, 2; 1, 0]).Cost, [0, 2; 1, 0]);
%! assert_equal (setCost (Mdl, []).Cost, [0, 1; 1, 0]);

%!error<ClassificationNeuralNetwork.margin: too few input arguments.> ...
%! margin (fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [1, 2])
%!error<ClassificationNeuralNetwork.loss: X is empty.> ...
%! loss (fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [], [1])
%!error<ClassificationNeuralNetwork.loss: X must have the same number of predictors> ...
%! loss (fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [1, 2, 3], [1])
%!error<ClassificationNeuralNetwork.loss: Y must have the same number of rows as X.> ...
%! loss (fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [1, 2], [1; 2])
%!error<ClassificationNeuralNetwork.loss: unsupported Loss function.> ...
%! loss (fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [1, 2], 1, 'LossFun', 'bogus')
%!error<ClassificationNeuralNetwork.edge: invalid parameter name in optional paired arguments.> ...
%! edge (fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [1, 2], 1, 'bogus', 1)
%!error<ClassificationNeuralNetwork.setCost: 'Cost' must be a square numeric matrix> ...
%! setCost (fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [0, 1, 2])
%!error<ClassificationNeuralNetwork.setPrior: 'Prior' must be 'empirical', 'uniform',> ...
%! setPrior (fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), 'bogus')

## RowsUsed is empty when every observation was used.
%!test
%! load fisheriris
%! X = meas;
%! Y = grp2idx (species);
%! Mdl = fitcnet (X, Y, 'IterationLimit', 20);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (class (Mdl.RowsUsed), 'double');
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (rows (Mdl.X), 150);
%! assert_equal (rows (Mdl.W), 150);

## A missing response drops its observation and RowsUsed marks it.
%!test
%! load fisheriris
%! X = meas;
%! Y = grp2idx (species);
%! Y(5) = NaN;
%! Mdl = fitcnet (X, Y, 'IterationLimit', 20);
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
%! X = meas;
%! X(3,2) = NaN;
%! Y = grp2idx (species);
%! Mdl = fitcnet (X, Y, 'IterationLimit', 20);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (rows (Mdl.X), 150);
%! assert_equal (sum (isnan (Mdl.X(:))), 1);

## Mu and Sigma are empty unless the predictors are standardized.
%!test
%! load fisheriris
%! Mdl = fitcnet (meas, species, 'IterationLimit', 20);
%! assert_equal (Mdl.Mu, []);
%! assert_equal (Mdl.Sigma, []);

## Standardizing weights the complete observations by each class's original
## share of the observation weight.  Values from MATLAB R2024a.
%!test
%! load fisheriris
%! X = meas;
%! X(3,2) = NaN; X(17,4) = NaN; X(140,1) = NaN;
%! Mdl = fitcnet (X, species, 'Standardize', true, 'IterationLimit', 20);
%! assert_equal (Mdl.Mu, [5.8405997732426291, 3.0547817460317455, ...
%!                        3.7612840136054406, 1.1980799319727886], 1e-13);
%! assert_equal (Mdl.Sigma, [0.82803317153591371, 0.43533400398915184, ...
%!                           1.7640762592568813, 0.76297286301694878], 1e-13);

## fitcnet takes Prior and Cost, and the prior reweights the observations.
## Values from R2024a.
%!test
%! load fisheriris
%! i3 = [1:50, 51:80, 101:120];
%! Mdl = fitcnet (meas(i3,:), species(i3), 'IterationLimit', 20, ...
%!                'Prior', [0.2, 0.3, 0.5]);
%! assert_equal (Mdl.Prior, [0.2, 0.3, 0.5], 1e-14);
%! assert_equal (Mdl.W(1), 0.004, 1e-14);
%! assert_equal (Mdl.W(51), 0.01, 1e-14);
%! assert_equal (Mdl.W(81), 0.025, 1e-14);

## Cost may be assigned after the fit.
%!test
%! load fisheriris
%! i3 = [1:50, 51:80, 101:120];
%! Mdl = fitcnet (meas(i3,:), species(i3), 'IterationLimit', 20);
%! Mdl.Cost = [0, 2, 3; 1, 0, 1; 1, 1, 0];
%! assert_equal (Mdl.Cost, [0, 2, 3; 1, 0, 1; 1, 1, 0]);

## Prior is fixed at construction and cannot be assigned afterwards.  The
## refusal comes from the property attributes, so the message is core
## Octave's and is not pinned here.
%!error ...
%! load fisheriris; ...
%! Mdl = fitcnet (meas, species, 'IterationLimit', 20); ...
%! Mdl.Prior = [0.2, 0.3, 0.5];
