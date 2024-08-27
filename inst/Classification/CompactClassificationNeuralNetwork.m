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

classdef CompactClassificationNeuralNetwork
## -*- texinfo -*-
## @deftypefn {statistics} CompactClassificationNeuralNetwork
##
## A @qcode{CompactClassificationNeuralNetwork} object is a compact version of a
## discriminant analysis model, @qcode{CompactClassificationNeuralNetwork}.
##
## The @qcode{CompactClassificationDiscriminant} does not include the training
## data resulting to a smaller classifier size, which can be used for making
## predictions from new data, but not for tasks such as cross validation.  It
## can only be created from a @qcode{ClassificationNeuralNetwork} model by using
## the @code{compact} object method.
##
## The available methods for a @qcode{CompactClassificationNeuralNetwork} object
## are:
## @itemize
## @item
## @code{predict}
## @item
## @code{savemodel}
## @end itemize
##
## @seealso{fitcdiscr, ClassificationDiscriminant}
## @end deftypefn

  properties (Access = public)

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

    ## constructor
    function this = CompactClassificationNeuralNetwork (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), "ClassificationNeuralNetwork"))
        error (strcat (["CompactClassificationNeuralNetwork: invalid"], ...
                       [" classification object."]));
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

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationNeuralNetwork} {@var{labels} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationNeuralNetwork} {[@var{labels}, @var{scores}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the Neural Network
    ## classification object.
    ##
    ## @code{@var{labels} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## trained neural network classification compact model in @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactClassificationNeuralNetwork} class
    ## object.
    ## @item
    ## @var{X} must be an @math{MxP} numeric matrix with the same number of
    ## predictors @math{P} as the corresponding predictors of the trained neural
    ## network compact model in @var{obj}.
    ## @end itemize
    ##
    ## @code{[@var{labels}, @var{scores}] = predict (@var{obj}, @var{XC}} also
    ## returns @var{scores}, which represent the probability of each label
    ## belonging to a specific class. For each observation in X, the predicted
    ## class label is the one with the highest score among all classes.
    ## Alternatively, @var{scores} can contain the posterior probabilities if
    ## the ScoreTransform has been previously set.
    ##
    ## @seealso{fitcnet, ClassificationNeuralNetwork,
    ## CompactClassificationNeuralNetwork}
    ## @end deftypefn

    function [labels, scores] = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error (strcat (["CompactClassificationNeuralNetwork.predict:"], ...
                       [" too few input arguments."]));
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("CompactClassificationNeuralNetwork.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat (["CompactClassificationNeuralNetwork.predict:"], ...
                       [" XC must have the same number of predictors"], ...
                       [" as the trained neural network model."]));
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
            error (strcat (["CompactClassificationNeuralNetwork.predict:"], ...
                           [" 'ScoreTransform' must be a"], ...
                           [" 'function_handle' object."]));
          endif
          scores = f (scores);
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNeuralNetwork} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationNeuralNetwork object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves a
    ## ClassificationNeuralNetwork object into a file defined by @var{filename}.
    ##
    ## @seealso{loadmodel, fitcnet, ClassificationNeuralNetwork, cvpartition,
    ## ClassificationPartitionedModel}
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
      DislayInfo              = this.DislayInfo;
      Solver                  = this.Solver;

      ## Save classdef name and all model properties as individual variables
      save (fname, "classdef_name", "NumPredictors", "PredictorNames", ...
            "ResponseName", "ClassNames", "ScoreTransform", "Standardize", ...
            "Sigma", "Mu", "LayerSizes", "Activations", ...
            "OutputLayerActivation", "LearningRate", "IterationLimit", ...
            "Solver", "ModelParameters", "ConvergenceInfo", "DislayInfo");
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
          error (strcat (["CompactClassificationNeuralNetwork.load_model:"], ...
                         [" invalid model in '%s'."]), filename)
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
%! CMdl = crossval (Mdl);
%!
%! whos ('Mdl', 'CMdl')

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
%!test
%! CMdl.ScoreTransform = "a";
%!error<CompactClassificationNeuralNetwork.predict: 'ScoreTransform' must be a 'function_handle' object.> ...
%! [labels, scores] = predict (CMdl, x);
