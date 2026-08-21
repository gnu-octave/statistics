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
## @deftypefn {statistics} {@var{obj} =} CompactRegressionSVM (@var{Mdl})
##
## Create a @qcode{CompactRegressionSVM} object, a support vector regression
## model that has dropped its training data.
##
## @code{@var{obj} = CompactRegressionSVM (@var{Mdl})} returns the compact
## form of the @qcode{RegressionSVM} object @var{Mdl}.  It is normally reached
## through @code{compact (@var{Mdl})} rather than called directly.
##
## The compact model keeps what is needed to answer about new data, the
## support vectors and their coefficients, the intercept, the kernel, the
## standardization and the response transform, and drops what only describes
## the fit: the predictor and response data, the observation weights, the rows
## used, the observation count, and which training rows became support
## vectors.  @code{predict} and @code{loss} therefore agree with the full
## model to the last digit, while @code{resubPredict} and @code{resubLoss} do
## not exist here, there being no training data left to resubstitute.
##
## @seealso{RegressionSVM, fitrsvm}
## @end deftypefn

classdef CompactRegressionSVM

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer scalar.  This property is read-only.
    ##
    ## @end deftp
    NumPredictors         = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    ResponseName          = [];


    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} Epsilon
    ##
    ## Half-width of the insensitive tube
    ##
    ## A non-negative scalar, carried over from the model this one was
    ## compacted from.  It is what the @qcode{'epsiloninsensitive'} loss
    ## charges against.  This property is read-only.
    ##
    ## @end deftp
    Epsilon               = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} Standardize
    ##
    ## Whether the predictor data was standardized
    ##
    ## A logical scalar.  When true @code{predict} centres and scales its
    ## input by @code{Mu} and @code{Sigma}, as training did.  This property is
    ## read-only.
    ##
    ## @end deftp
    Standardize           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} Sigma
    ##
    ## Standard deviation of the predictors
    ##
    ## A row vector with one entry per predictor, used for standardization.
    ## Empty if @qcode{Standardize} is @qcode{false}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Sigma                 = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} Mu
    ##
    ## Mean of the predictors
    ##
    ## A row vector with one entry per predictor, used for standardization.
    ## Empty if @qcode{Standardize} is @qcode{false}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Mu                    = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} ModelParameters
    ##
    ## Parameters the model was fitted with
    ##
    ## A structure holding the SVM formulation, the kernel and its parameters,
    ## the box constraint, @code{Epsilon} and the solver settings.  This
    ## property is read-only.
    ##
    ## @end deftp
    ModelParameters       = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} Model
    ##
    ## The trained LIBSVM model
    ##
    ## A structure as returned by @code{svmtrain} and consumed by
    ## @code{svmpredict}.  This property is read-only.
    ##
    ## @end deftp
    Model                 = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} Alpha
    ##
    ## Dual coefficients of the support vectors
    ##
    ## A numeric column vector with one entry per support vector, signed, as
    ## in the model this one was compacted from.  This property is read-only.
    ##
    ## @end deftp
    Alpha                 = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} Beta
    ##
    ## Primal coefficients, one per predictor
    ##
    ## A numeric column vector, equal to
    ## @code{obj.SupportVectors' * obj.Alpha}, and empty for any kernel other
    ## than linear.  This property is read-only.
    ##
    ## @end deftp
    Beta                  = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} Bias
    ##
    ## Intercept of the fitted function
    ##
    ## A numeric scalar.  With a linear kernel the prediction is
    ## @code{X * obj.Beta + obj.Bias}.  This property is read-only.
    ##
    ## @end deftp
    Bias                  = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} SupportVectors
    ##
    ## The support vectors themselves
    ##
    ## A numeric matrix with one row per support vector, on the scale the
    ## model was trained on.  This property is read-only.
    ##
    ## @end deftp
    SupportVectors        = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector of column indices, and empty when none is.  This
    ## property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionSVM} {property} ExpandedPredictorNames
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
    ## @deftp {CompactRegressionSVM} {property} ResponseTransform
    ##
    ## Transformation applied to the predicted response
    ##
    ## A function handle, applied by @code{predict} to the model's output.  It
    ## may be set after construction, either to a handle or to the name of a
    ## supported transformation.
    ##
    ## @end deftp
    ResponseTransform     = @(y) y;
  endproperties

  ## Readable by the counterpart class, which copies it, and kept out of
  ## the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    RTname = 'none';
  endproperties

  ## Set methods for the properties a user may assign.
  methods

    function this = set.ResponseTransform (this, val)
      name = 'CompactRegressionSVM';
      [this.ResponseTransform, this.RTname] = ...
            parseResponseTransform (val, name);
    endfunction

  endmethods

  methods(Hidden)

    ## constructor
    function this = CompactRegressionSVM (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), 'RegressionSVM'))
        error ("CompactRegressionSVM: invalid regression object.");
      endif

      ## Save properties to compact model.  The training data, the observation
      ## weights, the rows used, the observation count and IsSupportVector are
      ## deliberately left behind: each is sized to the training set, which is
      ## what "compact" exists not to carry.
      this.NumPredictors         = Mdl.NumPredictors;
      this.PredictorNames        = Mdl.PredictorNames;
      this.ResponseName          = Mdl.ResponseName;

      this.ResponseTransform     = Mdl.ResponseTransform;
      this.RTname                = Mdl.RTname;

      this.Epsilon               = Mdl.Epsilon;
      this.Standardize           = Mdl.Standardize;
      this.Sigma                 = Mdl.Sigma;
      this.Mu                    = Mdl.Mu;

      this.ModelParameters       = Mdl.ModelParameters;
      this.Model                 = Mdl.Model;

      this.Alpha                 = Mdl.Alpha;
      this.Beta                  = Mdl.Beta;
      this.Bias                  = Mdl.Bias;
      this.SupportVectors        = Mdl.SupportVectors;
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
      fprintf ("\n  CompactRegressionSVM\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %d\n", 'NumPredictors', this.NumPredictors);
      fprintf ("%+25s: '%s'\n", 'ResponseTransform', this.RTname);
      fprintf ("%+25s: %g\n", 'Epsilon', this.Epsilon);
      fprintf ("%+25s: [%dx1 double]\n", 'Alpha', numel (this.Alpha));
      if (! isempty (this.Beta))
        fprintf ("%+25s: [%dx1 double]\n", 'Beta', numel (this.Beta));
      endif
      fprintf ("%+25s: %f\n", 'Bias', this.Bias);
      fprintf ("%+25s: '%s'\n", 'KernelFunction', ...
               this.ModelParameters.KernelFunction);
      fprintf ("%+25s: [%dx%d double]\n", 'SupportVectors', ...
               rows (this.SupportVectors), columns (this.SupportVectors));
    endfunction



  endmethods

  methods(Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {CompactRegressionSVM} {@var{yFit} =} predict (@var{obj}, @var{XC})
    ##
    ## Predict the response for new data with a compact support vector
    ## regression model.
    ##
    ## @code{@var{yFit} = predict (@var{obj}, @var{XC})} returns a column
    ## vector holding the predicted response for each row of @var{XC}.  It
    ## agrees with the full model this object was compacted from.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactRegressionSVM} class object.
    ## @item
    ## @var{XC} must be a numeric matrix with the same number of predictors as
    ## the data the model was trained on.
    ## @end itemize
    ##
    ## @seealso{CompactRegressionSVM, RegressionSVM}
    ## @end deftypefn
    function yFit = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("CompactRegressionSVM.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("CompactRegressionSVM.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("CompactRegressionSVM.predict: XC must have the", ...
                       " same number of predictors as the trained model."));
      endif

      ## Standardize (if necessary)
      if (this.Standardize)
        XC = (XC - this.Mu) ./ this.Sigma;
      endif

      ## LIBSVM returns the fitted response as its first output for a
      ## regression model, there being no label to decide.
      yFit = svmpredict (zeros (rows (XC), 1), XC, this.Model, '-q');

      ## Apply ResponseTransform
      yFit = this.ResponseTransform (yFit);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionSVM} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactRegressionSVM} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the regression loss of a compact support vector machine model.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## weighted mean squared error between the response @var{Y} and the
    ## response the model predicts for @var{X}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactRegressionSVM} class object.
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
    ## @item @qcode{'LossFun'} @tab @qcode{'mse'}, the default,
    ## @qcode{'epsiloninsensitive'}, or a function handle called as
    ## @code{@var{lossfun} (@var{Y}, @var{yFit}, @var{W})} returning a scalar.
    ##
    ## @item @qcode{'Weights'} @tab A numeric vector of observation weights
    ## with one entry per row of @var{X}.  It defaults to a uniform weight.
    ## The weights are normalized to sum to one before the loss is formed.
    ## @end multitable
    ##
    ## @seealso{CompactRegressionSVM, RegressionSVM}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("CompactRegressionSVM.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactRegressionSVM.loss: Name-Value arguments", ...
                       " must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, 'loss');

      ## Defaults, then the optional pairs
      LossFun = 'mse';
      args = varargin;
      keep = true (1, numel (args));
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("CompactRegressionSVM.loss: parameter name must", ...
                         " be a character vector."));
        endif
        if (strcmpi (args{i}, 'lossfun'))
          LossFun = args{i+1};
          if (! (is_function_handle (LossFun) ||
                 (ischar (LossFun) && isrow (LossFun))))
            error (strcat ("CompactRegressionSVM.loss: 'LossFun' must be", ...
                           " a character vector or a function handle."));
          endif
          if (ischar (LossFun) && ! any (strcmpi (LossFun, ...
                                         {'mse', 'epsiloninsensitive'})))
            error (strcat ("CompactRegressionSVM.loss: unsupported", ...
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
          error (strcat ("CompactRegressionSVM.loss: 'LossFun' must", ...
                         " return a numeric scalar."));
        endif
      elseif (strcmpi (LossFun, 'epsiloninsensitive'))
        L = sum (W .* max (0, abs (Y - yFit) - this.Epsilon));
      else
        L = sum (W .* (Y - yFit) .^ 2);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactRegressionSVM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a compact support vector regression model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves every property of
    ## the @qcode{CompactRegressionSVM} object @var{obj} into @var{filename}
    ## in binary format, so that it can be read back with @code{loadmodel}.
    ##
    ## @seealso{loadmodel, CompactRegressionSVM}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error ("CompactRegressionSVM.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("CompactRegressionSVM.savemodel: FNAME must be", ...
                       " a character vector."));
      endif
      ## Generate variable for class name
      classdef_name = 'CompactRegressionSVM';

      ## Create variables from model properties
      NumPredictors           = this.NumPredictors;
      PredictorNames          = this.PredictorNames;
      ResponseName            = this.ResponseName;
      ResponseTransform       = this.ResponseTransform;
      Epsilon                 = this.Epsilon;
      Standardize             = this.Standardize;
      Sigma                   = this.Sigma;
      Mu                      = this.Mu;
      ModelParameters         = this.ModelParameters;
      Model                   = this.Model;
      Alpha                   = this.Alpha;
      Beta                    = this.Beta;
      Bias                    = this.Bias;
      SupportVectors          = this.SupportVectors;
      CategoricalPredictors   = this.CategoricalPredictors;
      ExpandedPredictorNames  = this.ExpandedPredictorNames;
      RTname                  = this.RTname;

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'NumPredictors', ...
            'PredictorNames', 'ResponseName', 'ResponseTransform', ...
            'Epsilon', 'Standardize', 'Sigma', 'Mu', 'ModelParameters', ...
            'Model', 'Alpha', 'Beta', 'Bias', 'SupportVectors', ...
            'CategoricalPredictors', 'ExpandedPredictorNames', 'RTname');
    endfunction

  endmethods

  methods (Access = private)

    ## Shared validation for the assessment methods, so each reports under
    ## its own name.
    function [X, Y] = checkXY_ (this, X, Y, caller)
      if (isempty (X))
        error ("CompactRegressionSVM.%s: X is empty.", caller);
      elseif (this.NumPredictors != columns (X))
        error (strcat ("CompactRegressionSVM.%s: X must have the same", ...
                       " number of predictors as the trained model."), caller);
      endif
      if (isempty (Y))
        error ("CompactRegressionSVM.%s: Y is empty.", caller);
      elseif (! (isnumeric (Y) && isreal (Y)))
        error (strcat ("CompactRegressionSVM.%s: Y must be a real", ...
                       " numeric vector."), caller);
      elseif (rows (X) != numel (Y))
        error (strcat ("CompactRegressionSVM.%s: Y must have the same", ...
                       " number of rows as X."), caller);
      endif
    endfunction

    ## Pull a "Weights" pair out of the optional arguments, defaulting to a
    ## uniform weight, and reject any other name.
    function W = getWeights_ (this, args, n, caller)
      W = ones (n, 1);
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("CompactRegressionSVM.%s: parameter name must", ...
                         " be a character vector."), caller);
        endif
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
          if (! (isnumeric (W) && isvector (W)))
            error (strcat ("CompactRegressionSVM.%s: 'Weights' must be", ...
                           " a numeric vector."), caller);
          endif
          if (numel (W) != n)
            error (strcat ("CompactRegressionSVM.%s: size of 'Weights'", ...
                           " must equal the number of rows in X."), caller);
          endif
        else
          error (strcat ("CompactRegressionSVM.%s: invalid parameter name", ...
                         " in optional paired arguments."), caller);
        endif
      endfor
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a CompactRegressionSVM object
      mdl = CompactRegressionSVM ();

      ## Get fieldnames from DATA (including private properties)
      names = fieldnames (data);

      ## Copy data into object
      for i = 1:numel (names)
        ## Check fieldnames in DATA match the class properties
        try
          mdl.(names{i}) = data.(names{i});
        catch
          error (strcat ("CompactRegressionSVM.load_model: invalid", ...
                         " model in '%s'."), filename)
        end_try_catch
      endfor
    endfunction

  endmethods

endclassdef


## The compact model keeps what answers about new data and drops the fit.
%!test
%! randn ('seed', 42);
%! X = randn (50, 2);
%! Y = X(:,1) - 2 * X(:,2);
%! CMdl = compact (RegressionSVM (X, Y));
%! assert_equal (class (CMdl), 'CompactRegressionSVM');
%! kept = {'NumPredictors', 'PredictorNames', 'ResponseName', ...
%!         'ResponseTransform', 'Epsilon', 'Standardize', 'Sigma', 'Mu', ...
%!         'Alpha', 'Beta', 'Bias', 'SupportVectors', ...
%!         'CategoricalPredictors', 'ExpandedPredictorNames'};
%! assert_equal (all (ismember (kept, properties (CMdl))), true);
%! dropped = {'X', 'Y', 'W', 'NumObservations', 'RowsUsed', ...
%!            'IsSupportVector'};
%! assert_equal (any (ismember (dropped, properties (CMdl))), false);

## A compact model predicts exactly what the model it came from predicts.
%!test
%! randn ('seed', 42);
%! X = randn (60, 3);
%! Y = X * [1; -2; 0.5] + 3;
%! Mdl = RegressionSVM (X, Y);
%! CMdl = compact (Mdl);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));
%! assert_equal (loss (CMdl, X, Y), loss (Mdl, X, Y));
%! assert_equal (loss (CMdl, X, Y), resubLoss (Mdl), 1e-12);
%! assert_equal (CMdl.Alpha, Mdl.Alpha);
%! assert_equal (CMdl.Beta, Mdl.Beta);
%! assert_equal (CMdl.Bias, Mdl.Bias);
%! assert_equal (CMdl.SupportVectors, Mdl.SupportVectors);
%! assert_equal (CMdl.Epsilon, Mdl.Epsilon);

## With a linear kernel the compact model is a plain linear function.
%!test
%! randn ('seed', 42);
%! X = randn (40, 2);
%! Y = X * [3; -1] + 2;
%! CMdl = compact (RegressionSVM (X, Y));
%! assert_equal (X * CMdl.Beta + CMdl.Bias, predict (CMdl, X), 1e-8);
%! assert_equal (CMdl.Beta, CMdl.SupportVectors' * CMdl.Alpha, 1e-12);

## Standardization travels with the compact model.
%!test
%! randn ('seed', 42);
%! X = [randn(60, 1), randn(60, 1) * 1000];
%! Y = X(:,1) + X(:,2) / 1000;
%! Mdl = RegressionSVM (X, Y, 'Standardize', true, 'Epsilon', 0.01);
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.Standardize, true);
%! assert_equal (CMdl.Mu, Mdl.Mu);
%! assert_equal (CMdl.Sigma, Mdl.Sigma);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));

## A non-linear kernel compacts too, and keeps Beta empty.
%!test
%! randn ('seed', 42);
%! X = randn (40, 2);
%! Y = sum (X .^ 2, 2);
%! Mdl = RegressionSVM (X, Y, 'KernelFunction', 'rbf');
%! CMdl = compact (Mdl);
%! assert_equal (isempty (CMdl.Beta), true);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));

## loss takes the same options as the full model's.
%!test
%! randn ('seed', 42);
%! X = randn (30, 2);
%! Y = X(:,1) * 3 + 1;
%! Mdl = RegressionSVM (X, Y, 'Epsilon', 0.5);
%! CMdl = compact (Mdl);
%! yFit = predict (CMdl, X);
%! assert_equal (loss (CMdl, X, Y), mean ((Y - yFit) .^ 2), 1e-12);
%! assert_equal (loss (CMdl, X, Y, 'LossFun', 'epsiloninsensitive'), ...
%!               mean (max (0, abs (Y - yFit) - 0.5)), 1e-12);
%! w = rand (30, 1) + 0.1;
%! assert_equal (loss (CMdl, X, Y, 'Weights', w), ...
%!               loss (CMdl, X, Y, 'Weights', 4 * w), 1e-12);
%! f = @(y, yf, ww) sum (ww .* abs (y - yf));
%! assert_equal (loss (CMdl, X, Y, 'LossFun', f), ...
%!               mean (abs (Y - yFit)), 1e-12);

## ResponseTransform travels, and can be replaced on the compact model.
%!test
%! X = [linspace(0, 1, 20)', linspace(1, 2, 20)'];
%! Mdl = RegressionSVM (X, 2 * X(:,1) + 1, 'ResponseTransform', 'exp');
%! CMdl = compact (Mdl);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));
%! CMdl.ResponseTransform = 'none';
%! assert_equal (predict (CMdl, X), log (predict (Mdl, X)), 1e-12);

## A saved compact model comes back carrying its own numbers.
%!test
%! randn ('seed', 42);
%! X = randn (30, 2);
%! Y = X(:,1) * 4 - 1;
%! CMdl = compact (RegressionSVM (X, Y, 'Standardize', true));
%! fname = tempname ();
%! savemodel (CMdl, fname);
%! C2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (C2), 'CompactRegressionSVM');
%! assert_equal (C2.Alpha, CMdl.Alpha);
%! assert_equal (C2.Bias, CMdl.Bias);
%! assert_equal (C2.Epsilon, CMdl.Epsilon);
%! assert_equal (C2.SupportVectors, CMdl.SupportVectors);
%! assert_equal (predict (C2, X), predict (CMdl, X));

## Test input validation for the constructor
%!error<CompactRegressionSVM: invalid regression object.> ...
%! CompactRegressionSVM (1)
%!error<CompactRegressionSVM: invalid regression object.> ...
%! CompactRegressionSVM (fitcsvm (ones (4, 2), [1; 1; 2; 2]))

## Test input validation for predict and loss
%!shared CRSVM
%! CRSVM = compact (RegressionSVM ([1, 1; 2, 1; 3, 2; 4, 2], [2; 4; 6; 8]));
%!error<CompactRegressionSVM.predict: too few input arguments.> ...
%! predict (CRSVM)
%!error<CompactRegressionSVM.predict: XC is empty.> ...
%! predict (CRSVM, [])
%!error<CompactRegressionSVM.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (CRSVM, ones (2, 3))
%!error<CompactRegressionSVM.loss: too few input arguments.> ...
%! loss (CRSVM)
%!error<CompactRegressionSVM.loss: Name-Value arguments must be in pairs.> ...
%! loss (CRSVM, [1, 1; 2, 1], [2; 4], 'Weights')
%!error<CompactRegressionSVM.loss: X is empty.> ...
%! loss (CRSVM, [], [2; 4])
%!error<CompactRegressionSVM.loss: X must have the same number of predictors as the trained model.> ...
%! loss (CRSVM, ones (2, 3), [2; 4])
%!error<CompactRegressionSVM.loss: Y is empty.> ...
%! loss (CRSVM, [1, 1; 2, 1], [])
%!error<CompactRegressionSVM.loss: Y must be a real numeric vector.> ...
%! loss (CRSVM, [1, 1; 2, 1], {'a'; 'b'})
%!error<CompactRegressionSVM.loss: Y must have the same number of rows as X.> ...
%! loss (CRSVM, [1, 1; 2, 1], [2; 4; 6])
%!error<CompactRegressionSVM.loss: 'LossFun' must be a character vector or a function handle.> ...
%! loss (CRSVM, [1, 1; 2, 1], [2; 4], 'LossFun', 5)
%!error<CompactRegressionSVM.loss: unsupported 'LossFun' value.> ...
%! loss (CRSVM, [1, 1; 2, 1], [2; 4], 'LossFun', 'mae')
%!error<CompactRegressionSVM.loss: 'LossFun' must return a numeric scalar.> ...
%! loss (CRSVM, [1, 1; 2, 1], [2; 4], 'LossFun', @(y, yf, w) [1, 2])
%!error<CompactRegressionSVM.loss: 'Weights' must be a numeric vector.> ...
%! loss (CRSVM, [1, 1; 2, 1], [2; 4], 'Weights', {'a'})
%!error<CompactRegressionSVM.loss: size of 'Weights' must equal the number of rows in X.> ...
%! loss (CRSVM, [1, 1; 2, 1], [2; 4], 'Weights', [1; 2; 3])
%!error<CompactRegressionSVM.loss: invalid parameter name in optional paired arguments.> ...
%! loss (CRSVM, [1, 1; 2, 1], [2; 4], 'Nope', 1)

## Test input validation for savemodel
%!error<CompactRegressionSVM.savemodel: too few input arguments.> ...
%! savemodel (CRSVM)
%!error<CompactRegressionSVM.savemodel: FNAME must be a character vector.> ...
%! savemodel (CRSVM, 5)

%!error<CompactRegressionSVM: unrecognized 'ResponseTransform' function.> ...
%! CRSVM.ResponseTransform = 'nope';
