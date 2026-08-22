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
## @deftypefn {statistics} {@var{obj} =} CompactRegressionNeuralNetwork (@var{Mdl})
##
## Create a @qcode{CompactRegressionNeuralNetwork} object, a neural network
## regression model that has dropped its training data.
##
## @code{@var{obj} = CompactRegressionNeuralNetwork (@var{Mdl})} returns the
## compact form of the @qcode{RegressionNeuralNetwork} object @var{Mdl}.  It is
## normally reached through @code{compact (@var{Mdl})} rather than called
## directly.
##
## The compact model keeps what is needed to answer about new data, the layer
## weights and biases, the activations, the standardization and the response
## transform, and drops what only describes the fit: the predictor and response
## data, the observation weights, the rows used, the number of observations and
## the iteration by iteration training history.  @code{predict} and
## @code{loss} therefore agree with the full model to the last digit, while
## @code{resubPredict} and @code{resubLoss} do not exist here, there being no
## training data left to resubstitute.
##
## @seealso{RegressionNeuralNetwork, fitrnet}
## @end deftypefn

classdef CompactRegressionNeuralNetwork

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer scalar.  This property is read-only.
    ##
    ## @end deftp
    NumPredictors         = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    ResponseName          = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} Sigma
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
    ## @deftp {CompactRegressionNeuralNetwork} {property} Mu
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
    ## @deftp {CompactRegressionNeuralNetwork} {property} LayerSizes
    ##
    ## Sizes of the fully connected hidden layers
    ##
    ## A row vector of positive integers, one per hidden layer.  This property
    ## is read-only.
    ##
    ## @end deftp
    LayerSizes            = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} Activations
    ##
    ## Activation functions of the hidden layers
    ##
    ## A character vector, or a cell array of character vectors with one entry
    ## per hidden layer.  This property is read-only.
    ##
    ## @end deftp
    Activations           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} OutputLayerActivation
    ##
    ## Activation function of the output layer
    ##
    ## A character vector.  @qcode{'none'} applies the identity, so a
    ## prediction is an unrestricted real number.  This property is read-only.
    ##
    ## @end deftp
    OutputLayerActivation = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} LearningRate
    ##
    ## Learning rate used to train the network
    ##
    ## A positive scalar value.  This property is read-only.
    ##
    ## @end deftp
    LearningRate          = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} IterationLimit
    ##
    ## Maximum number of training iterations
    ##
    ## A positive integer scalar.  This property is read-only.
    ##
    ## @end deftp
    IterationLimit        = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} ModelParameters
    ##
    ## Parameters of the trained network
    ##
    ## A structure as returned by @code{fcnntrain} and consumed by
    ## @code{fcnnpredict}.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters       = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} ConvergenceInfo
    ##
    ## Information recorded during training
    ##
    ## A structure with the fields @code{Time} and @code{TrainingLoss}, carried
    ## over from the model this one was compacted from.  This property is
    ## read-only.
    ##
    ## @end deftp
    ConvergenceInfo       = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} DisplayInfo
    ##
    ## Whether training printed its progress
    ##
    ## A logical scalar.  This property is read-only.
    ##
    ## @end deftp
    DisplayInfo           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} Solver
    ##
    ## Solver used to train the network
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    Solver                = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} LayerWeights
    ##
    ## Weights the network learned
    ##
    ## A cell array with one entry per layer, the output layer included.  This
    ## property is read-only.
    ##
    ## @end deftp
    LayerWeights          = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} LayerBiases
    ##
    ## Biases the network learned
    ##
    ## A cell array with one entry per layer, the output layer included.  This
    ## property is read-only.
    ##
    ## @end deftp
    LayerBiases           = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector of column indices, and empty when none is.  This
    ## property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the model expanded them
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};
  endproperties

  ## Properties a user may set after the model is built.  Each one is
  ## validated by its set method below.
  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionNeuralNetwork} {property} ResponseTransform
    ##
    ## Transformation applied to the predicted response
    ##
    ## A function handle, applied by @code{predict} to the network's output.
    ## It may be set after construction, either to a handle or to the name of
    ## a supported transformation.
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
      name = 'CompactRegressionNeuralNetwork';
      [this.RTfun, this.ResponseTransform] = parseResponseTransform (val, name);
    endfunction

    ## constructor
    function this = CompactRegressionNeuralNetwork (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), 'RegressionNeuralNetwork'))
        error (strcat ("CompactRegressionNeuralNetwork: invalid", ...
                       " regression object."));
      endif

      ## Save properties to compact model.  The training data, the observation
      ## weights, the rows used, the observation count and the training
      ## history are deliberately left behind: they describe the fit, not the
      ## model, and keeping them is what "compact" exists to avoid.
      this.NumPredictors         = Mdl.NumPredictors;
      this.PredictorNames        = Mdl.PredictorNames;
      this.ResponseName          = Mdl.ResponseName;

      this.ResponseTransform     = Mdl.ResponseTransform;
      this.RTfun                = Mdl.RTfun;

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
      fprintf ("\n  CompactRegressionNeuralNetwork\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
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

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {CompactRegressionNeuralNetwork} {@var{yFit} =} predict (@var{obj}, @var{XC})
    ##
    ## Predict the response for new data with a compact neural network
    ## regression model.
    ##
    ## @code{@var{yFit} = predict (@var{obj}, @var{XC})} returns a column
    ## vector holding the predicted response for each row of @var{XC}.  It
    ## agrees with the full model this object was compacted from.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactRegressionNeuralNetwork} class object.
    ## @item
    ## @var{XC} must be a numeric matrix with the same number of predictors as
    ## the data the model was trained on.
    ## @end itemize
    ##
    ## @seealso{CompactRegressionNeuralNetwork, RegressionNeuralNetwork}
    ## @end deftypefn
    function yFit = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error (strcat ("CompactRegressionNeuralNetwork.predict:", ...
                       " too few input arguments."));
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("CompactRegressionNeuralNetwork.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("CompactRegressionNeuralNetwork.predict:", ...
                       " XC must have the same number of predictors", ...
                       " as the trained neural network model."));
      endif

      ## Standardize (if necessary)
      if (! isempty (this.Mu))
        XC = (XC - this.Mu) ./ this.Sigma;
      endif

      ## The network's output is the second return value; the first is an
      ## index of the largest output, constant for a single regression unit.
      NumThreads = nproc ();
      [~, yFit] = fcnnpredict (this.ModelParameters, XC, NumThreads);

      ## Apply ResponseTransform
      yFit = this.RTfun (yFit);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionNeuralNetwork} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactRegressionNeuralNetwork} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the regression loss of a compact neural network model.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## weighted mean squared error between the response @var{Y} and the
    ## response the model predicts for @var{X}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactRegressionNeuralNetwork} class object.
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
    ## The weights are normalized to sum to one before the loss is formed.
    ## @end multitable
    ##
    ## @seealso{CompactRegressionNeuralNetwork, RegressionNeuralNetwork}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error (strcat ("CompactRegressionNeuralNetwork.loss:", ...
                       " too few input arguments."));
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactRegressionNeuralNetwork.loss: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, 'loss');

      ## Defaults, then the optional pairs
      LossFun = 'mse';
      args = varargin;
      keep = true (1, numel (args));
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("CompactRegressionNeuralNetwork.loss: parameter", ...
                         " name must be a character vector."));
        endif
        if (strcmpi (args{i}, 'lossfun'))
          LossFun = args{i+1};
          if (! (is_function_handle (LossFun) ||
                 (ischar (LossFun) && isrow (LossFun))))
            error (strcat ("CompactRegressionNeuralNetwork.loss: 'LossFun'", ...
                           " must be a character vector or a function", ...
                           " handle."));
          endif
          if (ischar (LossFun) && ! strcmpi (LossFun, 'mse'))
            error (strcat ("CompactRegressionNeuralNetwork.loss:", ...
                           " unsupported 'LossFun' value."));
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
          error (strcat ("CompactRegressionNeuralNetwork.loss: 'LossFun'", ...
                         " must return a numeric scalar."));
        endif
      else
        L = sum (W .* (Y - yFit) .^ 2);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactRegressionNeuralNetwork} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a compact neural network regression model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves every property of
    ## the @qcode{CompactRegressionNeuralNetwork} object @var{obj} into
    ## @var{filename} in binary format, so that it can be read back with
    ## @code{loadmodel}.
    ##
    ## @seealso{loadmodel, CompactRegressionNeuralNetwork}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error (strcat ("CompactRegressionNeuralNetwork.savemodel:", ...
                       " too few input arguments."));
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("CompactRegressionNeuralNetwork.savemodel: FNAME", ...
                       " must be a character vector."));
      endif
      ## Generate variable for class name
      classdef_name = 'CompactRegressionNeuralNetwork';

      ## Create variables from model properties
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
      CategoricalPredictors   = this.CategoricalPredictors;
      ExpandedPredictorNames  = this.ExpandedPredictorNames;
      RTfun                  = this.RTfun;

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'NumPredictors', ...
            'PredictorNames', 'ResponseName', 'ResponseTransform', ...
            'Sigma', 'Mu', 'LayerSizes', ...
            'Activations', 'OutputLayerActivation', 'LearningRate', ...
            'IterationLimit', 'ModelParameters', 'ConvergenceInfo', ...
            'DisplayInfo', 'Solver', 'LayerWeights', 'LayerBiases', ...
            'CategoricalPredictors', 'ExpandedPredictorNames', 'RTfun');
    endfunction

  endmethods

  methods (Access = private)

    ## Shared validation for the assessment methods, so each reports under
    ## its own name.
    function [X, Y] = checkXY_ (this, X, Y, caller)
      if (isempty (X))
        error ("CompactRegressionNeuralNetwork.%s: X is empty.", caller);
      elseif (this.NumPredictors != columns (X))
        error (strcat ("CompactRegressionNeuralNetwork.%s: X must have the", ...
                       " same number of predictors as the trained model."), ...
               caller);
      endif
      if (isempty (Y))
        error ("CompactRegressionNeuralNetwork.%s: Y is empty.", caller);
      elseif (! (isnumeric (Y) && isreal (Y)))
        error (strcat ("CompactRegressionNeuralNetwork.%s: Y must be a", ...
                       " real numeric vector."), caller);
      elseif (rows (X) != numel (Y))
        error (strcat ("CompactRegressionNeuralNetwork.%s: Y must have the", ...
                       " same number of rows as X."), caller);
      endif
    endfunction

    ## Pull a "Weights" pair out of the optional arguments, defaulting to a
    ## uniform weight, and reject any other name.
    function W = getWeights_ (this, args, n, caller)
      W = ones (n, 1);
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("CompactRegressionNeuralNetwork.%s: parameter", ...
                         " name must be a character vector."), caller);
        endif
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
          if (! (isnumeric (W) && isvector (W)))
            error (strcat ("CompactRegressionNeuralNetwork.%s: 'Weights'", ...
                           " must be a numeric vector."), caller);
          endif
          if (numel (W) != n)
            error (strcat ("CompactRegressionNeuralNetwork.%s: size of", ...
                           " 'Weights' must equal the number of", ...
                           " rows in X."), caller);
          endif
        else
          error (strcat ("CompactRegressionNeuralNetwork.%s: invalid", ...
                         " parameter name in optional paired", ...
                         " arguments."), caller);
        endif
      endfor
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a CompactRegressionNeuralNetwork object
      mdl = CompactRegressionNeuralNetwork ();

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
        ## Check fieldnames in DATA match the class properties
        try
          mdl.(names{i}) = data.(names{i});
        catch
          error (strcat ("CompactRegressionNeuralNetwork.load_model:", ...
                         " invalid model in '%s'."), filename)
        end_try_catch
      endfor
    endfunction

  endmethods

endclassdef


## The compact model keeps what answers about new data and drops the fit.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = [randn(40, 2); randn(40, 2) + 3];
%! Y = X(:,1) - 2 * X(:,2);
%! Mdl = fitrnet (X, Y, 'IterationLimit', 100);
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), 'CompactRegressionNeuralNetwork');
%! kept = {'NumPredictors', 'PredictorNames', 'ResponseName', ...
%!         'ResponseTransform', 'Sigma', 'Mu', 'LayerSizes', ...
%!         'Activations', 'OutputLayerActivation', 'LayerWeights', ...
%!         'LayerBiases', 'CategoricalPredictors', 'ExpandedPredictorNames'};
%! assert_equal (all (ismember (kept, properties (CMdl))), true);
%! dropped = {'X', 'Y', 'W', 'NumObservations', 'RowsUsed', 'TrainingHistory'};
%! assert_equal (any (ismember (dropped, properties (CMdl))), false);

## A compact model predicts exactly what the model it came from predicts.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = [randn(50, 2); randn(50, 2) + 2];
%! Y = 3 * X(:,1) + X(:,2);
%! Mdl = fitrnet (X, Y, 'LayerSizes', [8, 6], 'IterationLimit', 200);
%! CMdl = compact (Mdl);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));
%! assert_equal (loss (CMdl, X, Y), loss (Mdl, X, Y));
%! assert_equal (loss (CMdl, X, Y), resubLoss (Mdl), 1e-12);
%! assert_equal (CMdl.LayerWeights, Mdl.LayerWeights);
%! assert_equal (CMdl.LayerBiases, Mdl.LayerBiases);

## Standardization travels with the compact model, so it predicts on the
## scale the network was trained on.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = [randn(60, 1), randn(60, 1) * 1000];
%! Y = X(:,1) + X(:,2) / 1000;
%! Mdl = fitrnet (X, Y, 'Standardize', true, 'IterationLimit', 200);
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.Mu, Mdl.Mu);
%! assert_equal (CMdl.Sigma, Mdl.Sigma);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));

## loss takes the same options as the full model's.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = linspace (0, 1, 30)';
%! Y = 3 * X + 1;
%! CMdl = compact (fitrnet (X, Y, 'IterationLimit', 100));
%! yFit = predict (CMdl, X);
%! assert_equal (loss (CMdl, X, Y), mean ((Y - yFit) .^ 2), 1e-12);
%! assert_equal (loss (CMdl, X, Y, 'LossFun', 'mse'), loss (CMdl, X, Y), 1e-12);
%! w = rand (30, 1) + 0.1;
%! assert_equal (loss (CMdl, X, Y, 'Weights', w), ...
%!               loss (CMdl, X, Y, 'Weights', 5 * w), 1e-12);
%! f = @(y, yf, ww) sum (ww .* abs (y - yf));
%! assert_equal (loss (CMdl, X, Y, 'LossFun', f), mean (abs (Y - yFit)), 1e-12);

## ResponseTransform travels, and can be replaced on the compact model.
%!test
%! rand ('seed', 42);
%! X = linspace (0, 1, 20)';
%! Mdl = fitrnet (X, 2 * X + 1, 'ResponseTransform', 'exp', ...
%!                'IterationLimit', 50);
%! CMdl = compact (Mdl);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));
%! CMdl.ResponseTransform = 'none';
%! assert_equal (predict (CMdl, X), log (predict (Mdl, X)), 1e-12);

## A saved compact model comes back carrying its own numbers.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = linspace (0, 1, 30)';
%! Y = 4 * X - 1;
%! CMdl = compact (fitrnet (X, Y, 'LayerSizes', [6, 4], ...
%!                          'IterationLimit', 60));
%! fname = tempname ();
%! savemodel (CMdl, fname);
%! C2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (C2), 'CompactRegressionNeuralNetwork');
%! assert_equal (C2.LayerWeights, CMdl.LayerWeights);
%! assert_equal (C2.LayerBiases, CMdl.LayerBiases);
%! assert_equal (C2.NumPredictors, CMdl.NumPredictors);
%! assert_equal (C2.LayerSizes, CMdl.LayerSizes);
%! assert_equal (predict (C2, X), predict (CMdl, X));

## Test input validation for the constructor
%!error<CompactRegressionNeuralNetwork: invalid regression object.> ...
%! CompactRegressionNeuralNetwork (1)
%!error<CompactRegressionNeuralNetwork: invalid regression object.> ...
%! CompactRegressionNeuralNetwork (fitcnet (ones (4, 2), [1; 1; 2; 2]))

## Test input validation for predict and loss
%!shared CRNN
%! rand ('seed', 42);
%! CRNN = compact (fitrnet ([1; 2; 3; 4], [2; 4; 6; 8], 'IterationLimit', 10));
%!error<CompactRegressionNeuralNetwork.predict: too few input arguments.> ...
%! predict (CRNN)
%!error<CompactRegressionNeuralNetwork.predict: XC is empty.> ...
%! predict (CRNN, [])
%!error<CompactRegressionNeuralNetwork.predict: XC must have the same number of predictors as the trained neural network model.> ...
%! predict (CRNN, ones (2, 3))
%!error<CompactRegressionNeuralNetwork.loss: too few input arguments.> ...
%! loss (CRNN)
%!error<CompactRegressionNeuralNetwork.loss: Name-Value arguments must be in pairs.> ...
%! loss (CRNN, [1; 2], [2; 4], 'Weights')
%!error<CompactRegressionNeuralNetwork.loss: X is empty.> ...
%! loss (CRNN, [], [2; 4])
%!error<CompactRegressionNeuralNetwork.loss: X must have the same number of predictors as the trained model.> ...
%! loss (CRNN, ones (2, 3), [2; 4])
%!error<CompactRegressionNeuralNetwork.loss: Y is empty.> ...
%! loss (CRNN, [1; 2], [])
%!error<CompactRegressionNeuralNetwork.loss: Y must be a real numeric vector.> ...
%! loss (CRNN, [1; 2], {'a'; 'b'})
%!error<CompactRegressionNeuralNetwork.loss: Y must have the same number of rows as X.> ...
%! loss (CRNN, [1; 2], [2; 4; 6])
%!error<CompactRegressionNeuralNetwork.loss: 'LossFun' must be a character vector or a function handle.> ...
%! loss (CRNN, [1; 2], [2; 4], 'LossFun', 5)
%!error<CompactRegressionNeuralNetwork.loss: unsupported 'LossFun' value.> ...
%! loss (CRNN, [1; 2], [2; 4], 'LossFun', 'mae')
%!error<CompactRegressionNeuralNetwork.loss: 'LossFun' must return a numeric scalar.> ...
%! loss (CRNN, [1; 2], [2; 4], 'LossFun', @(y, yf, w) [1, 2])
%!error<CompactRegressionNeuralNetwork.loss: 'Weights' must be a numeric vector.> ...
%! loss (CRNN, [1; 2], [2; 4], 'Weights', {'a'})
%!error<CompactRegressionNeuralNetwork.loss: size of 'Weights' must equal the number of rows in X.> ...
%! loss (CRNN, [1; 2], [2; 4], 'Weights', [1; 2; 3])
%!error<CompactRegressionNeuralNetwork.loss: invalid parameter name in optional paired arguments.> ...
%! loss (CRNN, [1; 2], [2; 4], 'Nope', 1)

## Test input validation for savemodel
%!error<CompactRegressionNeuralNetwork.savemodel: too few input arguments.> ...
%! savemodel (CRNN)
%!error<CompactRegressionNeuralNetwork.savemodel: FNAME must be a character vector.> ...
%! savemodel (CRNN, 5)

%!error<CompactRegressionNeuralNetwork: unrecognized 'ResponseTransform' function.> ...
%! CRNN.ResponseTransform = 'nope';
%!error<CompactRegressionNeuralNetwork: function handle for 'ResponseTransform' must return the same size as its input.> ...
%! CRNN.ResponseTransform = @(y) [y; y];

## A fitted model survives savemodel and loadmodel: the properties come
## back as they were and it predicts the same.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! Y = meas(:,1);
%! Mdl = compact (fitrnet (X, Y, 'IterationLimit', 20));
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'CompactRegressionNeuralNetwork');
%! assert_equal (M2.PredictorNames, Mdl.PredictorNames);
%! assert_equal (class (M2.ResponseTransform), class (Mdl.ResponseTransform));
%! assert_equal (predict (M2, X(1:5,:)), predict (Mdl, X(1:5,:)), 1e-12);
