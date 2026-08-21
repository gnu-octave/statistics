## Copyright (C) 2024-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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

classdef CompactClassificationNeuralNetwork
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactClassificationNeuralNetwork
  ##
  ## Compact neural network classification
  ##
  ## The @code{CompactClassificationNeuralNetwork} class implements a compact
  ## version of the neural network classifier object, which can predict
  ## responses for new data using the @code{predict} method, but does not store
  ## the training data.
  ##
  ## A compact neural network classification model is a smaller version of the
  ## full @code{ClassificationNeuralNetwork} model that does not include the
  ## training data.  It consumes less memory than the full model, but cannot
  ## perform tasks that require the training data, such as cross-validation.
  ##
  ## Create a @code{CompactClassificationNeuralNetwork} object by using the
  ## @code{compact} method on a @code{ClassificationNeuralNetwork} object.
  ##
  ## @seealso{ClassificationNeuralNetwork, fitcnet}
  ## @end deftp

  properties(Access = public)
    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer value specifying the number of predictors in the
    ## training dataset used for training the neural network model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumPredictors         = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} PredictorNames
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
    ## @deftp {CompactClassificationNeuralNetwork} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName          = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} ClassNames
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
    ## @deftp {CompactClassificationNeuralNetwork} {property} ScoreTransform
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
    ## @item @qcode{'invlogit'} @tab @math{1 ./ (1 + exp (-x))}
    ## @item @qcode{'ismax'} @tab Sets the score for the class with the
    ## largest score to 1, and for all other classes to 0
    ## @item @qcode{'logit'} @tab @math{log (x ./ (1 - x))}
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

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} Standardize
    ##
    ## Flag to standardize predictors
    ##
    ## A boolean flag indicating whether the predictor data has been
    ## standardized prior to training.  When @qcode{true}, the predictors are
    ## centered and scaled to have zero mean and unit variance.  This property
    ## is read-only.
    ##
    ## @end deftp
    Standardize           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} Sigma
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
    ## @deftp {CompactClassificationNeuralNetwork} {property} Mu
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
    ## @deftp {CompactClassificationNeuralNetwork} {property} LayerSizes
    ##
    ## Sizes of fully connected layers
    ##
    ## A positive integer vector specifying the sizes of the fully connected
    ## layers in the neural network model.  The i-th element of
    ## @qcode{LayerSizes} is the number of outputs in the i-th fully connected
    ## layer of the neural network model.  @qcode{LayerSizes} does not include
    ## the size of the final fully connected layer.  This layer always has K
    ## outputs, where K is the number of classes in Y.  This property is
    ## read-only.
    ##
    ## @end deftp
    LayerSizes            = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} Activations
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
    ## @deftp {CompactClassificationNeuralNetwork} {property} OutputLayerActivation
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
    ## @deftp {CompactClassificationNeuralNetwork} {property} LearningRate
    ##
    ## Learning rate for gradient descent
    ##
    ## A positive scalar value defining the learning rate used by the gradient
    ## descent algorithm during training.  This property is read-only.
    ##
    ## @end deftp
    LearningRate          = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} IterationLimit
    ##
    ## Maximum number of training iterations
    ##
    ## A positive integer value defining the maximum number of epochs for
    ## training the model.  This property is read-only.
    ##
    ## @end deftp
    IterationLimit        = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} ModelParameters
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
    ## @deftp {CompactClassificationNeuralNetwork} {property} ConvergenceInfo
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
    ## @deftp {CompactClassificationNeuralNetwork} {property} DisplayInfo
    ##
    ## Display training information flag
    ##
    ## A boolean flag indicating whether to print information during training.
    ## This property is read-only.
    ##
    ## @end deftp
    DisplayInfo           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} Solver
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
    ## @deftp {CompactClassificationNeuralNetwork} {property} LayerWeights
    ##
    ## Learned weights of each fully connected layer
    ##
    ## A cell array holding one matrix per layer, the output layer included,
    ## with one row per neuron of that layer and one column per input it
    ## takes.  This property is read-only.
    ##
    ## @end deftp
    LayerWeights          = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} LayerBiases
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
    ## @deftp {CompactClassificationNeuralNetwork} {property} Cost
    ##
    ## Cost of misclassification
    ##
    ## A numeric matrix with one row and one column per class, where
    ## @code{Cost(i,j)} is the cost of classifying an observation of class
    ## @math{i} as class @math{j}.  It is taken from the model this object was
    ## compacted from.  This property is read-only.
    ##
    ## @end deftp
    Cost                  = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} Prior
    ##
    ## Prior probability of each class
    ##
    ## A numeric vector with one entry per class, in the order of
    ## @code{ClassNames}, summing to one.  It is taken from the model this
    ## object was compacted from.  This property is read-only.
    ##
    ## @end deftp
    Prior                 = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector holding the column of each predictor treated as
    ## categorical, and empty when none is.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNeuralNetwork} {property} ExpandedPredictorNames
    ##
    ## Names of the expanded predictor variables
    ##
    ## A cell array of character vectors naming the predictors as the model
    ## sees them.  It matches @code{PredictorNames} unless a categorical
    ## predictor was expanded into dummy variables.  This property is
    ## read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};
  endproperties

  properties(Access = private, Hidden)
    STname = 'none';
  endproperties

  methods(Hidden)

    ## constructor
    function this = CompactClassificationNeuralNetwork (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), 'ClassificationNeuralNetwork'))
        error (strcat ("CompactClassificationNeuralNetwork: invalid", ...
                       " classification object."));
      endif

      ## Save properties to compact model
      this.NumPredictors         = Mdl.NumPredictors;
      this.PredictorNames        = Mdl.PredictorNames;
      this.ResponseName          = Mdl.ResponseName;
      this.ClassNames            = Mdl.ClassNames;

      this.ScoreTransform        = Mdl.ScoreTransform;
      this.STname                = Mdl.STname;

      this.Standardize           = Mdl.Standardize;
      this.Sigma                 = Mdl.Sigma;
      this.Mu                    = Mdl.Mu;

      this.LayerSizes            = Mdl.LayerSizes;
      this.Activations           = Mdl.Activations;
      this.OutputLayerActivation = Mdl.OutputLayerActivation;
      this.LearningRate          = Mdl.LearningRate;
      this.IterationLimit        = Mdl.IterationLimit;

      this.ModelParameters       = Mdl.ModelParameters;
      this.ConvergenceInfo       = Mdl.ConvergenceInfo;
      this.DisplayInfo           = Mdl.DisplayInfo;
      this.Solver                = Mdl.Solver;

      this.LayerWeights          = Mdl.LayerWeights;
      this.LayerBiases           = Mdl.LayerBiases;
      this.Cost                  = Mdl.Cost;
      this.Prior                 = Mdl.Prior;
      this.CategoricalPredictors = Mdl.CategoricalPredictors;
      this.ExpandedPredictorNames = Mdl.ExpandedPredictorNames;

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
      fprintf ("\n  CompactClassificationNeuralNetwork\n\n");
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

    ## Class specific subscripted reference
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("Invalid () indexing for referencing values", ...
                         " in a CompactClassificationNeuralNetwork object."));
        case '{}'
          error (strcat ("Invalid {} indexing for referencing values", ...
                         " in a CompactClassificationNeuralNetwork object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("CompactClassificationNeuralNetwork.subsref: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("CompactClassificationNeuralNetwork.subsref:", ...
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
        error (strcat ("CompactClassificationNeuralNetwork.subsasgn:", ...
                       " chained subscripts not allowed."));
      endif
      switch s.type
        case '()'
          error (strcat ("Invalid () indexing for assigning values", ...
                         " to a CompactClassificationNeuralNetwork object."));
        case '{}'
          error (strcat ("Invalid {} indexing for assigning values", ...
                         " to a CompactClassificationNeuralNetwork object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("CompactClassificationNeuralNetwork.subsasgn: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          switch (s.subs)
            case 'ScoreTransform'
              name = 'CompactClassificationNeuralNetwork';
              [this.ScoreTransform, this.STname] = parseScoreTransform (val, ...
                                                                        name);
            otherwise
              error (strcat ("CompactClassificationNeuralNetwork.subsasgn:", ...
                             " unrecognized or read-only property: '%s'"), ...
                             s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods(Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationNeuralNetwork} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationNeuralNetwork} {[@var{label}, @var{score}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the neural network
    ## classification model from a CompactClassificationNeuralNetwork object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## neural network model stored in the CompactClassificationNeuralNetwork
    ## model, @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactClassificationNeuralNetwork} class
    ## object.
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
    ## @seealso{CompactClassificationNeuralNetwork,
    ## ClassificationNeuralNetwork, fitcnet}
    ## @end deftypefn

    function [labels, scores] = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error (strcat ("CompactClassificationNeuralNetwork.predict:", ...
                       " too few input arguments."));
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("CompactClassificationNeuralNetwork.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("CompactClassificationNeuralNetwork.predict:", ...
                       " XC must have the same number of predictors", ...
                       " as the trained neural network model."));
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
    ## @deftypefn  {CompactClassificationNeuralNetwork} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin of a compact neural network classifier.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns a column
    ## vector holding, for each row of @var{X}, the score the model gives its
    ## true class in @var{Y} less the largest score it gives any other class.
    ## A positive margin means the observation is classified correctly, and
    ## the larger it is the more confidently so.
    ##
    ## @seealso{CompactClassificationNeuralNetwork,
    ## ClassificationNeuralNetwork, edge, loss, predict}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error (strcat ("CompactClassificationNeuralNetwork.margin:", ...
                       " too few input arguments."));
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
    ## @deftypefn  {CompactClassificationNeuralNetwork} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationNeuralNetwork} {@var{e} =} edge (@dots{}, @qcode{"Weights"}, @var{w})
    ##
    ## Classification edge of a compact neural network classifier.
    ##
    ## @code{@var{e} = edge (@var{obj}, @var{X}, @var{Y})} returns the mean of
    ## the classification margins over the rows of @var{X}.
    ##
    ## @code{@var{e} = edge (@dots{}, @qcode{"Weights"}, @var{w})} takes the
    ## weighted mean instead, with one weight per row of @var{X}.
    ##
    ## @seealso{CompactClassificationNeuralNetwork,
    ## ClassificationNeuralNetwork, margin, loss, predict}
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error (strcat ("CompactClassificationNeuralNetwork.edge:", ...
                       " too few input arguments."));
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactClassificationNeuralNetwork.edge:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, "edge");
      W = getWeights_ (this, varargin, rows (X), "edge");

      m = margin (this, X, Y);
      e = sum (W(:) .* m) / sum (W);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationNeuralNetwork} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationNeuralNetwork} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of a compact neural network classifier.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the loss of
    ## the model on the rows of @var{X} against the true labels @var{Y}.
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
    ## @seealso{CompactClassificationNeuralNetwork,
    ## ClassificationNeuralNetwork, margin, edge, predict}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error (strcat ("CompactClassificationNeuralNetwork.loss:", ...
                       " too few input arguments."));
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactClassificationNeuralNetwork.loss:", ...
                       " Name-Value arguments must be in pairs."));
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
            error (strcat ("CompactClassificationNeuralNetwork.loss:", ...
                           " 'LossFun' must be a character vector."));
          endif
          LossFun = tolower (LossFun);
          if (! any (strcmpi (LossFun, lossnames)))
            error (strcat ("CompactClassificationNeuralNetwork.loss:", ...
                           " unsupported Loss function."));
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
    ## @deftypefn  {CompactClassificationNeuralNetwork} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a CompactClassificationNeuralNetwork object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## CompactClassificationNeuralNetwork object into an Octave binary file, the
    ## name of which is specified in @var{filename}, along with an extra
    ## variable, which defines the type classification object these variables
    ## constitute.  Use @code{loadmodel} in order to load a classification
    ## object into Octave's workspace.
    ##
    ## @seealso{loadmodel, fitcnet, ClassificationNeuralNetwork}
    ## @end deftypefn

    function savemodel (this, fname)
      if (nargin < 2)
        error ("CompactClassificationNeuralNetwork.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error ("CompactClassificationNeuralNetwork.savemodel: FNAME must be a character vector.");
      endif
      ## Generate variable for class name
      classdef_name = 'CompactClassificationNeuralNetwork';

      ## Create variables from model properties
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
      CategoricalPredictors   = this.CategoricalPredictors;
      ExpandedPredictorNames  = this.ExpandedPredictorNames;
      STname                  = this.STname;

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'NumPredictors', ...
            'PredictorNames', 'ResponseName', 'ClassNames', ...
            'ScoreTransform', 'Standardize', 'Sigma', 'Mu', 'LayerSizes', ...
            'Activations', 'OutputLayerActivation', 'LearningRate', ...
            'IterationLimit', 'ModelParameters', 'ConvergenceInfo', ...
            'DisplayInfo', 'Solver', 'LayerWeights', 'LayerBiases', ...
            'Cost', 'Prior', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'STname');
    endfunction

  endmethods

  methods (Access = private)

    ## Shared validation for the assessment methods, so each reports under
    ## its own name.
    function [X, Y] = checkXY_ (this, X, Y, caller)
      if (isempty (X))
        error ("CompactClassificationNeuralNetwork.%s: X is empty.", caller);
      elseif (this.NumPredictors != columns (X))
        error (strcat ("CompactClassificationNeuralNetwork.%s: X must have", ...
                       " the same number of predictors as the trained", ...
                       " neural network model."), caller);
      endif
      if (isempty (Y))
        error ("CompactClassificationNeuralNetwork.%s: Y is empty.", caller);
      elseif (rows (X) != rows (Y))
        error (strcat ("CompactClassificationNeuralNetwork.%s: Y must have", ...
                       " the same number of rows as X."), caller);
      endif
    endfunction

    ## Pull a "Weights" pair out of the optional arguments, defaulting to a
    ## uniform weight, and reject any other name.
    function W = getWeights_ (this, args, n, caller)
      W = ones (n, 1);
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("CompactClassificationNeuralNetwork.%s: parameter", ...
                         " name must be a character vector."), caller);
        endif
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
          if (! (isnumeric (W) && isvector (W)))
            error (strcat ("CompactClassificationNeuralNetwork.%s:", ...
                           " 'Weights' must be a numeric vector."), caller);
          endif
          if (numel (W) != n)
            error (strcat ("CompactClassificationNeuralNetwork.%s: size of", ...
                           " 'Weights' must equal the number of", ...
                           " rows in X."), caller);
          endif
        else
          error (strcat ("CompactClassificationNeuralNetwork.%s: invalid", ...
                         " parameter name in optional paired", ...
                         " arguments."), caller);
        endif
      endfor
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationNeuralNetwork object
      mdl = CompactClassificationNeuralNetwork ();

      ## Get fieldnames from DATA (including private properties)
      names = fieldnames (data);

      ## Copy data into object
      for i = 1:numel (names)
        ## Check fieldnames in DATA match properties in
        ## CompactClassificationNeuralNetwork
        try
          mdl.(names{i}) = data.(names{i});
        catch
          error (strcat ("CompactClassificationNeuralNetwork.load_model:", ...
                         " invalid model in '%s'."), filename)
        end_try_catch
      endfor
    endfunction

  endmethods

endclassdef

%!demo
%! ## Train a neural network classifier and take its compact version, which
%! ## drops the training data but predicts identically.
%!
%! load fisheriris
%! X = meas;
%! Y = species;
%!
%! Mdl = fitcnet (X, Y, 'IterationLimit', 100)
%! CMdl = compact (Mdl)
%!
%! ## The compact model keeps no training data
%! isprop (Mdl, 'X')
%! isprop (CMdl, 'X')
%!
%! ## and predicts the same labels
%! isequal (predict (Mdl, X), predict (CMdl, X))

## Test input validation for constructor
%!error<CompactClassificationNeuralNetwork: invalid classification object.> ...
%! CompactClassificationNeuralNetwork (1)

## Test output for predict method
%!shared x, y, CMdl
%! load fisheriris
%! x = meas;
%! y = grp2idx (species);
%! Mdl = fitcnet (x, y, 'IterationLimit', 100);
%! CMdl = compact (Mdl);

## Test input validation for predict method
%!error<CompactClassificationNeuralNetwork.predict: too few input arguments.> ...
%! predict (CMdl)
%!error<CompactClassificationNeuralNetwork.predict: XC is empty.> ...
%! predict (CMdl, [])
%!error<CompactClassificationNeuralNetwork.predict: XC must have the same number of predictors as the trained neural network model.> ...
%! predict (CMdl, 1)

## Test input validation for assigning a new ScoreTransform
%!error<CompactClassificationNeuralNetwork: unrecognized 'ScoreTransform' function.> ...
%! CMdl.ScoreTransform = 'a';
## The compact model carries the trained parameters and the class bookkeeping.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], ...
%!                'LayerSizes', [3, 2], 'IterationLimit', 20);
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.LayerWeights, Mdl.LayerWeights);
%! assert_equal (CMdl.LayerBiases, Mdl.LayerBiases);
%! assert_equal (CMdl.Cost, Mdl.Cost);
%! assert_equal (CMdl.Prior, Mdl.Prior);
%! assert_equal (CMdl.CategoricalPredictors, Mdl.CategoricalPredictors);
%! assert_equal (CMdl.ExpandedPredictorNames, Mdl.ExpandedPredictorNames);

## margin, edge and loss agree with the model the object was compacted from.
%!test
%! load fisheriris
%! Mdl = fitcnet (meas, species, 'IterationLimit', 20);
%! CMdl = compact (Mdl);
%! assert_equal (margin (CMdl, meas, species), margin (Mdl, meas, species));
%! assert_equal (edge (CMdl, meas, species), edge (Mdl, meas, species));
%! assert_equal (loss (CMdl, meas, species), loss (Mdl, meas, species));

## Every loss function agrees with the full model's.
%!test
%! load fisheriris
%! Mdl = fitcnet (meas, species, 'IterationLimit', 20);
%! CMdl = compact (Mdl);
%! names = {'binodeviance', 'classifcost', 'classiferror', 'crossentropy', ...
%!          'exponential', 'hinge', 'logit', 'mincost', 'quadratic'};
%! for k = 1:numel (names)
%!   assert_equal (loss (CMdl, meas, species, 'LossFun', names{k}), ...
%!                 loss (Mdl, meas, species, 'LossFun', names{k}));
%! endfor

## A weighted edge weights the margins it averages.
%!test
%! load fisheriris
%! CMdl = compact (fitcnet (meas, species, 'IterationLimit', 20));
%! w = [ones(50, 1); 2 * ones(50, 1); 3 * ones(50, 1)];
%! m = margin (CMdl, meas, species);
%! assert_equal (edge (CMdl, meas, species, 'Weights', w), ...
%!               sum (w .* m) / sum (w), 1e-12);

## Test input validation for margin method
%!error<CompactClassificationNeuralNetwork.margin: too few input arguments.> ...
%! margin (CMdl)
%!error<CompactClassificationNeuralNetwork.margin: too few input arguments.> ...
%! margin (CMdl, x)
%!error<CompactClassificationNeuralNetwork.margin: X is empty.> ...
%! margin (CMdl, [], y)
%!error<CompactClassificationNeuralNetwork.margin: X must have the same number of predictors as the trained neural network model.> ...
%! margin (CMdl, 1, y)
%!error<CompactClassificationNeuralNetwork.margin: Y is empty.> ...
%! margin (CMdl, x, [])
%!error<CompactClassificationNeuralNetwork.margin: Y must have the same number of rows as X.> ...
%! margin (CMdl, x, y(1:10))

## Test input validation for edge method
%!error<CompactClassificationNeuralNetwork.edge: too few input arguments.> ...
%! edge (CMdl, x)
%!error<CompactClassificationNeuralNetwork.edge: Name-Value arguments must be in pairs.> ...
%! edge (CMdl, x, y, 'Weights')
%!error<CompactClassificationNeuralNetwork.edge: invalid parameter name in optional paired arguments.> ...
%! edge (CMdl, x, y, 'LossFun', 'hinge')
%!error<CompactClassificationNeuralNetwork.edge: 'Weights' must be a numeric vector.> ...
%! edge (CMdl, x, y, 'Weights', 'a')
%!error<CompactClassificationNeuralNetwork.edge: size of 'Weights' must equal the number of rows in X.> ...
%! edge (CMdl, x, y, 'Weights', [1, 2, 3])

## Test input validation for loss method
%!error<CompactClassificationNeuralNetwork.loss: too few input arguments.> ...
%! loss (CMdl, x)
%!error<CompactClassificationNeuralNetwork.loss: Name-Value arguments must be in pairs.> ...
%! loss (CMdl, x, y, 'LossFun')
%!error<CompactClassificationNeuralNetwork.loss: 'LossFun' must be a character vector.> ...
%! loss (CMdl, x, y, 'LossFun', 1)
%!error<CompactClassificationNeuralNetwork.loss: unsupported Loss function.> ...
%! loss (CMdl, x, y, 'LossFun', 'nonsense')

## A saved and reloaded compact model carries every property it holds.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], ...
%!                'LayerSizes', [3, 2], 'IterationLimit', 20);
%! CMdl = compact (Mdl);
%! fname = tempname ();
%! savemodel (CMdl, fname);
%! CMdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (CMdl2.ModelParameters, CMdl.ModelParameters);
%! assert_equal (CMdl2.LayerSizes, CMdl.LayerSizes);
%! assert_equal (CMdl2.ClassNames, CMdl.ClassNames);
%! assert_equal (CMdl2.ConvergenceInfo, CMdl.ConvergenceInfo);

## A reloaded compact model predicts exactly what the full model did.
%!test
%! load fisheriris
%! Mdl = fitcnet (meas, species, 'IterationLimit', 20);
%! fname = tempname ();
%! savemodel (compact (Mdl), fname);
%! CMdl2 = loadmodel (fname);
%! delete (fname);
%! [label, score] = predict (Mdl, meas);
%! [label2, score2] = predict (CMdl2, meas);
%! assert_equal (label2, label);
%! assert_equal (score2, score);

## A non-default ScoreTransform survives compacting and a save and load.
%!test
%! Mdl = fitcnet ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'IterationLimit', 20);
%! Mdl.ScoreTransform = 'symmetric';
%! fname = tempname ();
%! savemodel (compact (Mdl), fname);
%! CMdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (CMdl2.ScoreTransform (0.25), -0.5);

%!error <CompactClassificationNeuralNetwork.savemodel: too few input arguments.> ...
%! savemodel (CompactClassificationNeuralNetwork ())
%!error <CompactClassificationNeuralNetwork.savemodel: FNAME must be a character vector.> ...
%! savemodel (CompactClassificationNeuralNetwork (), 1)
%!error <CompactClassificationNeuralNetwork.savemodel: FNAME must be a character vector.> ...
%! savemodel (CompactClassificationNeuralNetwork (), ['ab'; 'cd'])
