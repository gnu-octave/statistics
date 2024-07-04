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

    X                                 = [];     # Predictor data
    Y                                 = [];     # Class labels

    NumClasses                        = [];
    ClassNames                        = [];
    NumObservations                   = [];

    LayerSizes                        = [];
    LayerWeights                      = [];
    LayerBiases                       = [];
    Activations                       = [];
    OutputLayerActivation             = [];

    ModelParameters                   = [];
    ConvergenceInfo                   = [];
    TrainingHistory                   = [];
    Solver                            = [];

  endproperties

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

      ## Remove missing values from X and Y
      RowsUsed  = ! logical (sum (isnan ([X, gY]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

      this.NumObservations = rows (X);

      ## Validate that Y is numeric
      if (!isnumeric(Y))
        error ("ClassificationNeuralNetwork: Y must be a numeric array.");
      endif

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
