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

  properties (Access = public)
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
    ## original predictor data @var{X} specifying which rows have been used for
    ## fitting the ClassificationNeuralNetwork model.  This property is
    ## read-only.
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
    ## variables.  The names are in the order in which the appear in the
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
    ## built-in functions.  Nevertherless, the @qcode{ScoreTransform} property
    ## always stores their function handle equivalent.
    ##
    ## @multitable @columnfractions 0.2 0.05 0.75
    ## @headitem @var{Value} @tab @tab @var{Description}
    ## @item @qcode{"doublelogit"} @tab @tab @math{1 ./ (1 + exp .^ (-2 * x))}
    ## @item @qcode{"invlogit"} @tab @tab @math{log (x ./ (1 - x))}
    ## @item @qcode{"ismax"} @tab @tab Sets the score for the class with the
    ## largest score to 1, and for all other classes to 0
    ## @item @qcode{"logit"} @tab @tab @math{1 ./ (1 + exp .^ (-x))}
    ## @item @qcode{"none"} @tab @tab @math{x} (no transformation)
    ## @item @qcode{"identity"} @tab @tab @math{x} (no transformation)
    ## @item @qcode{"sign"} @tab @tab @math{-1 for x < 0, 0 for x = 0, 1 for x > 0}
    ## @item @qcode{"symmetric"} @tab @tab @math{2 * x + 1}
    ## @item @qcode{"symmetricismax"} @tab @tab Sets the score for the class
    ## with the largest score to 1, and for all other classes to -1
    ## @item @qcode{"symmetriclogit"} @tab @tab @math{2 ./ (1 + exp .^ (-x)) - 1}
    ## @end multitable
    ##
    ## @end deftp
    ScoreTransform        = @(x) x;

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} Standardize
    ##
    ## Flag to standardize predictors
    ##
    ## A boolean flag indicating whether the predictor data has been standardized
    ## prior to training.  When @qcode{true}, the predictors are centered and
    ## scaled to have zero mean and unit variance.  This property is read-only.
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
    ## @end deftp
    Mu                    = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNeuralNetwork} {property} LayerSizes
    ##
    ## Sizes of fully connected layers
    ##
    ## A positive integer vector specifying the sizes of the fully connected
    ## layers in the neural network model.  The i-th element of @qcode{LayerSizes}
    ## is the number of outputs in the i-th fully connected layer of the neural
    ## network model.  @qcode{LayerSizes} does not include the size of the final
    ## fully connected layer.  This layer always has K outputs, where K is the
    ## number of classes in Y.  This property is read-only.
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
    ## Supported activation functions include: @qcode{"linear"},
    ## @qcode{"sigmoid"}, @qcode{"relu"}, @qcode{"tanh"}, @qcode{"softmax"},
    ## @qcode{"lrelu"}, @qcode{"prelu"}, @qcode{"elu"}, and @qcode{"gelu"}.
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
    ## for the @qcode{Activations} property.  This property is read-only.
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
    ## neural network model.  Currently only @qcode{"Gradient Descend"} is
    ## supported.  This property is read-only.
    ##
    ## @end deftp
    Solver                = [];
  endproperties

  properties (Access = private, Hidden)
    STname = 'none';
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
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.STname);
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
            error (strcat ("ClassificationNeuralNetwork.subsref:", ...
                           " unrecognized property: '%s'"), s.subs);
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
              [this.ScoreTransform, this.STname] = parseScoreTransform (val, ...
                                                                        name);
            otherwise
              error (strcat ("ClassificationNeuralNetwork.subsasgn:", ...
                             " unrecognized or read-only property: '%s'"), ...
                             s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods (Access = public)

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
    ## @code{X} must be a @math{NxP} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.  @var{X} will be used to train the neural network model.
    ## @item
    ## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels
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
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'PredictorNames'} @tab @tab A cell array of character
    ## vectors specifying the names of the predictors. The length of this array
    ## must match the number of columns in @var{X}.
    ##
    ## @item @qcode{'ResponseName'} @tab @tab A character vector specifying the
    ## name of the response variable.
    ##
    ## @item @qcode{'ClassNames'} @tab @tab Names of the classes in the class
    ## labels, @var{Y}, used for fitting the neural network model.
    ## @qcode{ClassNames} are of the same type as the class labels in @var{Y}.
    ##
    ## @item @qcode{'ScoreTransform'} @tab @tab A user-defined function handle
    ## or a character vector specifying one of the following builtin functions
    ## specifying the transformation applied to predicted classification scores.
    ## Supported values include @qcode{'doublelogit'}, @qcode{'invlogit'},
    ## @qcode{'ismax'}, @qcode{'logit'}, @qcode{'none'}, @qcode{'identity'},
    ## @qcode{'sign'}, @qcode{'symmetric'}, @qcode{'symmetricismax'}, and
    ## @qcode{'symmetriclogit'}.
    ##
    ## @item @qcode{'Standardize'} @tab @tab A logical scalar specifying whether
    ## to standardize the predictor data.  When @qcode{true}, the predictors are
    ## centered and scaled to have zero mean and unit variance.
    ##
    ## @item @qcode{'LayerSizes'} @tab @tab A positive integer vector specifying
    ## the sizes of the fully connected layers in the neural network.  The
    ## default is 10.
    ##
    ## @item @qcode{'Activations'} @tab @tab A character vector or cell array of
    ## character vectors specifying the activation functions for the hidden
    ## layers.  Supported values include @qcode{'linear'}, @qcode{'sigmoid'},
    ## @qcode{'relu'}, @qcode{'tanh'}, @qcode{'softmax'}, @qcode{'lrelu'},
    ## @qcode{'prelu'}, @qcode{'elu'}, and @qcode{'gelu'}.  The default is
    ## @qcode{'sigmoid'}.
    ##
    ## @item @qcode{'OutputLayerActivation'} @tab @tab A character vector
    ## specifying the activation function for the output layer.  Supported
    ## values are the same as for @qcode{'Activations'}.  The default is
    ## @qcode{'sigmoid'}.
    ##
    ## @item @qcode{'LearningRate'} @tab @tab A positive scalar specifying the
    ## learning rate for gradient descent.  The default is 0.01.
    ##
    ## @item @qcode{'IterationLimit'} @tab @tab A positive integer specifying
    ## the maximum number of training iterations.  The default is 1000.
    ##
    ## @item @qcode{'DisplayInfo'} @tab @tab A logical scalar specifying whether
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
            [this.ScoreTransform, this.STname] = parseScoreTransform ...
                                                 (varargin{2}, name);

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

      ## Generate default predictors and response variable names (if necessary)
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

      ## Assign the number of observations and their corresponding indices
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
    ## @var{XC} must be an @math{MxP} numeric matrix with the same number of
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
      X = this.X(this.RowsUsed, :);

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
        partition = cvpartition (this.NumObservations, 'LeaveOut');
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
      ## Create a compact model
      CVMdl = CompactClassificationNeuralNetwork (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationNeuralNetwork object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## ClassificationNeuralNetwork object into an Octave binary file, the name of
    ## which is specified in @var{filename}, along with an extra variable, which
    ## defines the type classification object these variables constitute.  Use
    ## @code{loadmodel} in order to load a classification object into Octave's
    ## workspace.
    ##
    ## @seealso{loadmodel, fitcnet, ClassificationNeuralNetwork}
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
      DisplayInfo             = this.DisplayInfo;
      Solver                  = this.Solver;
      STname                  = this.STname;

      ## Save classdef name and all model properties as individual variables
      save ("-binary", fname, "classdef_name", "X", "Y", "NumObservations", ...
            "RowsUsed", "NumPredictors", "PredictorNames", "ResponseName", ...
            "ClassNames", "ScoreTransform", "Standardize", "Sigma", "Mu", ...
            "LayerSizes", "Activations", "OutputLayerActivation", ...
            "LearningRate", "IterationLimit", "Solver", "ModelParameters", ...
            "ConvergenceInfo", "DisplayInfo", "STname");
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
    otherwise
      error (strcat ("ClassificationNeuralNetwork: misspelling or unsupported", ...
                     " 'Activation' function: '%s'."), strCode);
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
%! status = warning;
%! warning ('off');
%! rand ("seed", 23);
%! CVMdl = crossval (Mdl, "KFold", 5);
%! warning (status);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 5)
%! assert (class (CVMdl.Trained{1}), "CompactClassificationNeuralNetwork")
%! assert (CVMdl.CrossValidatedModel, "ClassificationNeuralNetwork")
%!test
%! status = warning;
%! warning ('off');
%! rand ("seed", 23);
%! CVMdl = crossval (Mdl, "HoldOut", 0.2);
%! warning (status);
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
