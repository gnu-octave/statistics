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

  properties (Access = public)
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
    ## @multitable @columnfractions 0.2 0.05 0.75
    ## @headitem @var{Value} @tab @tab @var{Description}
    ## @item @qcode{"doublelogit"} @tab @tab @math{1 ./ (1 + exp (-2 * x))}
    ## @item @qcode{"invlogit"} @tab @tab @math{1 ./ (1 + exp (-x))}
    ## @item @qcode{"ismax"} @tab @tab Sets the score for the class with the
    ## largest score to 1, and for all other classes to 0
    ## @item @qcode{"logit"} @tab @tab @math{log (x ./ (1 - x))}
    ## @item @qcode{"none"} @tab @tab @math{x} (no transformation)
    ## @item @qcode{"identity"} @tab @tab @math{x} (no transformation)
    ## @item @qcode{"sign"} @tab @tab @math{-1 for x < 0, 0 for x = 0, 1 for x > 0}
    ## @item @qcode{"symmetric"} @tab @tab @math{2 * x - 1}
    ## @item @qcode{"symmetricismax"} @tab @tab Sets the score for the class
    ## with the largest score to 1, and for all other classes to -1
    ## @item @qcode{"symmetriclogit"} @tab @tab @math{2 ./ (1 + exp (-x)) - 1}
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
    ## Supported activation functions include: @qcode{"linear"},
    ## @qcode{"sigmoid"}, @qcode{"relu"}, @qcode{"tanh"}, @qcode{"softmax"},
    ## @qcode{"lrelu"}, @qcode{"prelu"}, @qcode{"elu"}, and @qcode{"gelu"}.
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
    ## neural network model.  Currently only @qcode{"Gradient Descent"} is
    ## supported.  This property is read-only.
    ##
    ## @end deftp
    Solver                = [];
  endproperties

  methods (Hidden)

    ## constructor
    function this = CompactClassificationNeuralNetwork (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), "ClassificationNeuralNetwork"))
        error (strcat ("CompactClassificationNeuralNetwork: invalid", ...
                       " classification object."));
      endif

      ## Save properties to compact model
      this.NumPredictors         = Mdl.NumPredictors;
      this.PredictorNames        = Mdl.PredictorNames;
      this.ResponseName          = Mdl.ResponseName;
      this.ClassNames            = Mdl.ClassNames;

      this.ScoreTransform        = Mdl.ScoreTransform;

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
        str = repmat ({"'%s'"}, 1, numel (this.ClassNames));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.ClassNames{:});
      else # numeric
        str = repmat ({"%d"}, 1, numel (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, this.ClassNames);
      endif
      fprintf ("%+25s: '%s'\n", 'ClassNames', str);
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("%+25s: '%d'\n", 'NumPredictors', this.NumPredictors);
      str = repmat ({"%d"}, 1, numel (this.LayerSizes));
      str = strcat ('[', strjoin (str, ' '), ']');
      str = sprintf (str, this.LayerSizes);
      fprintf ("%+25s: '%s'\n", 'LayerSizes', str);
      if (iscellstr (this.Activations))
        str = repmat ({"'%s'"}, 1, numel (this.Activations));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.Activations{:});
        fprintf ("%+25s: '%s'\n", 'Activations', str);
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
              name = "CompactClassificationNeuralNetwork";
              this.ScoreTransform = parseScoreTransform (val, name);
            otherwise
              error (strcat ("CompactClassificationNeuralNetwork.subsasgn:", ...
                             " unrecognized or read-only property: '%s'"), ...
                             s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods (Access = public)

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
      [labels, scores] = fcnnpredict (this.ModelParameters, XC);

      # Get class labels
      labels = this.ClassNames(labels);

      if (nargout > 1)
        ## Apply ScoreTransform to return probability estimates
        if (! strcmp (this.ScoreTransform, "none"))
          f = this.ScoreTransform;
          if (! strcmp (class (f), "function_handle"))
            error (strcat ("CompactClassificationNeuralNetwork.predict:", ...
                           " 'ScoreTransform' must be a", ...
                           " 'function_handle' object."));
          endif
          scores = f (scores);
        endif
      endif

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
      ## Generate variable for class name
      classdef_name = "CompactClassificationNeuralNetwork";

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

      ## Save classdef name and all model properties as individual variables
      save ("-binary", fname, "classdef_name", "NumPredictors", ...
            "PredictorNames", "ResponseName", "ClassNames", ...
            "ScoreTransform", "Standardize", "Sigma", "Mu", "LayerSizes", ...
            "Activations", "OutputLayerActivation", "LearningRate", ...
            "IterationLimit", "ModelParameters", "ConvergenceInfo", ...
            "DisplayInfo", "Solver");
    endfunction

  endmethods

  methods (Static, Hidden)

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
%! ## Create a neural network classifier and its compact version
%! # and compare their size
%!
%! load fisheriris
%! X = meas;
%! Y = species;
%!
%! Mdl = fitcnet (X, Y, 'ClassNames', unique (species))
%! CMdl = crossval (Mdl)

## Test input validation for constructor
%!error<CompactClassificationDiscriminant: invalid classification object.> ...
%! CompactClassificationDiscriminant (1)

## Test output for predict method
%!shared x, y, CMdl
%! load fisheriris
%! x = meas;
%! y = grp2idx (species);
%! Mdl = fitcnet (x, y, "IterationLimit", 100);
%! CMdl = compact (Mdl);

## Test input validation for predict method
%!error<CompactClassificationNeuralNetwork.predict: too few input arguments.> ...
%! predict (CMdl)
%!error<CompactClassificationNeuralNetwork.predict: XC is empty.> ...
%! predict (CMdl, [])
%!error<CompactClassificationNeuralNetwork.predict: XC must have the same number of predictors as the trained neural network.> ...
%! predict (CMdl, 1)

## Test input validation for assigning a new ScoreTransform
%!error<CompactClassificationNeuralNetwork: unrecognized 'ScoreTransform' function.> ...
%! CMdl.ScoreTransform = "a";
