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
## @deftypefn {statistics} {@var{obj} =} RegressionPartitionedModel (@var{Mdl}, @var{Partition})
##
## Create a @qcode{RegressionPartitionedModel} object, a regression model
## cross validated over a partition of its training data.
##
## @code{@var{obj} = RegressionPartitionedModel (@var{Mdl}, @var{Partition})}
## refits @var{Mdl} once per fold of @var{Partition}, each time on the
## observations that fold holds out of its test set, and stores the compact
## form of every fit in @code{Trained}.  It is normally reached through
## @code{crossval (@var{Mdl})} rather than called directly.
##
## @itemize
## @item
## @var{Mdl} must be a @qcode{RegressionGAM}, a
## @qcode{RegressionNeuralNetwork}, or a @qcode{RegressionSVM} object.
## @item
## @var{Partition} must be a @qcode{cvpartition} object over as many
## observations as @var{Mdl} was trained on.
## @end itemize
##
## Every observation is held out by exactly one fold under @math{k}-fold or
## leave-one-out partitioning, so @code{kfoldPredict} can answer for it with a
## model that never saw it.  Under a holdout partition only the test set is
## answered for, and the rest come back @code{NaN}.
##
## @seealso{crossval, cvpartition, RegressionGAM, RegressionNeuralNetwork,
## RegressionSVM}
## @end deftypefn

classdef RegressionPartitionedModel

  properties (GetAccess = public, SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix holding the observations the model was trained on,
    ## the rows carrying missing values already removed.  This property is
    ## read-only.
    ##
    ## @end deftp
    X                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} Y
    ##
    ## Response data
    ##
    ## A numeric column vector with one entry per row of @code{X}.  This
    ## property is read-only.
    ##
    ## @end deftp
    Y                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector with one entry per observation.  This property
    ## is read-only.
    ##
    ## @end deftp
    W                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## A cell array with one entry per predictor, holding that predictor's bin
    ## edges where the learner discretized it before fitting.  It is carried
    ## over from the model that was cross validated, and is empty whenever that
    ## model did no binning, which is every learner this package implements:
    ## MATLAB fills it only for its generalized additive model, which bins
    ## because it is built from boosted trees where ours is built from splines.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    BinEdges              = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} CrossValidatedModel
    ##
    ## Name of the cross-validated model
    ##
    ## A character vector holding the short name of the learner that was
    ## cross validated, as MATLAB reports it: @qcode{'GAM'}, @qcode{'GP'},
    ## @qcode{'NeuralNetwork'} or @qcode{'SVM'}.  It is not the class name of
    ## that learner, and the classification side uses the same names.  This
    ## property is read-only.
    ##
    ## @end deftp
    CrossValidatedModel   = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector of column indices, and empty when none is.  This
    ## property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    ResponseName          = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer scalar.  This property is read-only.
    ##
    ## @end deftp
    NumObservations       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} ModelParameters
    ##
    ## Parameters the folds were fitted with
    ##
    ## A structure, carried over from the model that was cross validated.
    ## This property is read-only.
    ##
    ## @end deftp
    ModelParameters       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} Trained
    ##
    ## The models fitted to each fold
    ##
    ## A cell array with one compact model per fold, each fitted on the
    ## observations its fold holds out of the test set.  This property is
    ## read-only.
    ##
    ## @end deftp
    Trained               = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} KFold
    ##
    ## Number of folds
    ##
    ## A positive integer scalar.  This property is read-only.
    ##
    ## @end deftp
    KFold                 = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} Partition
    ##
    ## The partition the folds came from
    ##
    ## A @code{cvpartition} object.  This property is read-only.
    ##
    ## @end deftp
    Partition             = [];

  endproperties

  ## Properties a user may set.  Each one is validated by its set method.
  properties (GetAccess = public, SetAccess = public)
    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} ResponseTransform
    ##
    ## Transformation applied to the predicted response
    ##
    ## A function handle, carried over from the model that was cross
    ## validated.  This property is read-only.
    ##
    ## @end deftp
    ResponseTransform     = 'none';
  endproperties

  ## Copied from the parent model and kept out of the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    RTfun = @(y) y;
  endproperties

  ## Set methods for the properties a user may assign.
  methods (Hidden)

    function this = set.ResponseTransform (this, val)
      [f, nm] = parseResponseTransform (val, 'RegressionPartitionedModel');
      this.ResponseTransform = nm;
      this.RTfun = f;
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
      fprintf ("\n  RegressionPartitionedModel\n\n");
      fprintf ("%+25s: '%s'\n", 'CrossValidatedModel', ...
               this.CrossValidatedModel);
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'KFold', this.KFold);
      fprintf ("%+25s: '%s'\n", 'ResponseTransform', this.ResponseTransform);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {RegressionPartitionedModel} {@var{obj} =} RegressionPartitionedModel (@var{Mdl}, @var{Partition})
    ##
    ## Create a @qcode{RegressionPartitionedModel} object.
    ##
    ## See the class documentation for what it holds and how it is reached.
    ##
    ## @seealso{crossval, RegressionPartitionedModel}
    ## @end deftypefn
    function this = RegressionPartitionedModel (Mdl, Partition)

      ## Check input arguments
      if (nargin < 2)
        error ("RegressionPartitionedModel: too few input arguments.");
      endif

      ## Check for valid Regression object
      validTypes = {'RegressionGAM', 'RegressionGP', ...
                    'RegressionNeuralNetwork', 'RegressionSVM'};
      if (! any (strcmp (class (Mdl), validTypes)))
        error ("RegressionPartitionedModel: unsupported model type.");
      endif

      ## Check for valid cvpartition object
      if (! strcmp (class (Partition), 'cvpartition'))
        error (strcat ("RegressionPartitionedModel: invalid", ...
                       " 'cvpartition' object."));
      endif

      ## The partition indexes the observations actually used for training,
      ## so the rows dropped for missing values are removed here as well and
      ## every index below refers to the same set.
      X = Mdl.X;
      Y = Mdl.Y;
      Y = Y(:);

      if (Partition.NumObservations != rows (X))
        error (strcat ("RegressionPartitionedModel: 'cvpartition' object", ...
                       " must be defined over the %d observations the", ...
                       " model was trained on."), rows (X));
      endif

      ## Set properties
      this.X = X;
      this.Y = Y;
      this.W = Mdl.W;
      this.BinEdges = Mdl.BinEdges;
      this.KFold = Partition.NumTestSets;
      this.Trained = cell (this.KFold, 1);
      this.ResponseName = Mdl.ResponseName;
      this.NumObservations = rows (X);
      this.PredictorNames = Mdl.PredictorNames;
      this.CategoricalPredictors = Mdl.CategoricalPredictors;
      this.Partition = Partition;
      ## MATLAB stores a short name here, shared with the classification
      ## side: 'GAM', 'GP', 'NeuralNetwork', 'SVM'.  All measured.
      this.CrossValidatedModel = strrep (class (Mdl), 'Regression', '');
      this.ResponseTransform = Mdl.ResponseTransform;
      this.RTfun = Mdl.RTfun;
      ## A GAM carries no fitting parameter struct, so the property keeps its
      ## default rather than being read off a model that does not declare it.
      if (! strcmp (this.CrossValidatedModel, 'GAM'))
        this.ModelParameters = Mdl.ModelParameters;
      endif

      ## Switch Regression object types
      switch (this.CrossValidatedModel)

        case 'GAM'
          ## Knots, Order and DoF are three views of one parameterisation and
          ## the constructor accepts any two, recomputing the third, so only
          ## Knots and Order are passed on.
          args = {};
          GAMparams = {'PredictorNames', 'ResponseName', 'Formula', ...
                       'Knots', 'Order', 'Tol'};
          for i = 1:numel (GAMparams)
            paramName = GAMparams{i};
            paramValue = Mdl.(paramName);
            if (! isempty (paramValue))
              args = [args, {paramName, paramValue}];
            endif
          endfor

          ## Interactions now holds the fitted pairs, which the constructor
          ## does not take as a specification.  The term matrix does, and it
          ## reproduces the parent's terms exactly rather than re-selecting
          ## them.  A formula names its own terms and is passed instead, so
          ## this must not be passed alongside one.
          if (isempty (Mdl.Formula) && ! isempty (Mdl.IntMatrix))
            args = [args, {'Interactions', Mdl.IntMatrix}];
          endif
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            tmp = fitrgam (X(idx, :), Y(idx), args{:});
            this.Trained{k} = compact (tmp);
          endfor

        case 'NeuralNetwork'
          ## Computed before the cell literal: inside braces a space before
          ## the paren would split the call from its argument.
          stdz = ! isempty (Mdl.Mu);
          args = {'LayerSizes', Mdl.LayerSizes, ...
                  'Activations', Mdl.Activations, ...
                  'OutputLayerActivation', Mdl.OutputLayerActivation, ...
                  'LearningRate', Mdl.LearningRate, ...
                  'IterationLimit', Mdl.IterationLimit, ...
                  'Standardize', stdz, ...
                  'ResponseName', Mdl.ResponseName, ...
                  'PredictorNames', Mdl.PredictorNames};
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            tmp = fitrnet (X(idx, :), Y(idx), args{:});
            this.Trained{k} = compact (tmp);
          endfor

        case 'GP'
          p = Mdl.ModelParameters;
          args = {'KernelFunction', Mdl.KernelFunction, ...
                  'BasisFunction', Mdl.BasisFunction, ...
                  'FitMethod', p.FitMethod, ...
                  'PredictMethod', p.PredictMethod, ...
                  'Optimizer', p.Optimizer, ...
                  'ConstantSigma', p.ConstantSigma, ...
                  'SigmaLowerBound', p.SigmaLowerBound, ...
                  'Standardize', p.Standardize, ...
                  'ResponseName', Mdl.ResponseName, ...
                  'PredictorNames', Mdl.PredictorNames};
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            tmp = fitrgp (X(idx, :), Y(idx), args{:});
            this.Trained{k} = compact (tmp);
          endfor

        case 'SVM'
          p = Mdl.ModelParameters;
          stdz = ! isempty (Mdl.Mu);
          args = {'SVMtype', p.SVMtype, ...
                  'KernelFunction', p.KernelFunction, ...
                  'PolynomialOrder', p.PolynomialOrder, ...
                  'KernelScale', p.KernelScale, ...
                  'KernelOffset', p.KernelOffset, ...
                  'BoxConstraint', p.BoxConstraint, ...
                  'Epsilon', p.Epsilon, 'Nu', p.Nu, ...
                  'CacheSize', p.CacheSize, ...
                  'Tolerance', p.Tolerance, ...
                  'Shrinking', p.Shrinking, ...
                  'Standardize', stdz, ...
                  'ResponseName', Mdl.ResponseName, ...
                  'PredictorNames', Mdl.PredictorNames};
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            tmp = fitrsvm (X(idx, :), Y(idx), args{:});
            this.Trained{k} = compact (tmp);
          endfor

      endswitch

      ## No fold carries the transform: the parent applies it once to the
      ## assembled prediction.  Cleared here rather than at each fit call so
      ## that a backing added later cannot reintroduce a double application
      ## by inheriting a class default of its own.
      ##
      ## Safe because the transform here is a caller's preference and not
      ## fitted content.
      for k = 1:this.KFold
        T = this.Trained{k};
        T.ResponseTransform = 'none';
        this.Trained{k} = T;
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionPartitionedModel} {@var{yFit} =} kfoldPredict (@var{obj})
    ##
    ## Predict the response of every observation from the fold that held it
    ## out.
    ##
    ## @code{@var{yFit} = kfoldPredict (@var{obj})} returns a column vector
    ## with one entry per observation, each predicted by the fold's model that
    ## did not see it during training.  An observation no fold tests, which a
    ## holdout partition leaves outside its test set, comes back @code{NaN}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionPartitionedModel} class object.
    ## @end itemize
    ##
    ## @seealso{RegressionPartitionedModel, kfoldLoss}
    ## @end deftypefn
    function yFit = kfoldPredict (this)

      yFit = nan (this.NumObservations, 1);

      for k = 1:this.KFold
        testIdx = test (this.Partition, k);
        if (! any (testIdx))
          continue;
        endif
        yFit(testIdx) = predict (this.Trained{k}, this.X(testIdx, :));
      endfor

      ## As on the classification side, the folds never carry the transform and
      ## it is applied once to the assembled predictions.  MathWorks documents
      ## ResponseTransform as the function for transforming the predicted
      ## response values and documents assigning it by dot notation, so this is
      ## what the property is for; it used to be stored and never read.
      yFit = this.RTfun (yFit);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionPartitionedModel} {@var{L} =} kfoldLoss (@var{obj})
    ## @deftypefnx {RegressionPartitionedModel} {@var{L} =} kfoldLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the cross-validated regression loss.
    ##
    ## @code{@var{L} = kfoldLoss (@var{obj})} returns the weighted mean
    ## squared error between the response and the out-of-fold predictions of
    ## @code{kfoldPredict}, over every observation some fold tests.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionPartitionedModel} class object.
    ## @end itemize
    ##
    ## @code{@var{L} = kfoldLoss (@dots{}, @var{name}, @var{value})} accepts
    ## the following @qcode{Name-Value} pairs.
    ##
    ## @multitable @columnfractions 0.24 0.76
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab @qcode{'mse'}, the default,
    ## @qcode{'epsiloninsensitive'}, or a function handle called as
    ## @code{@var{lossfun} (@var{Y}, @var{yFit}, @var{W})} returning a scalar.
    ## The @math{epsilon}-insensitive loss belongs to a support vector model
    ## and is refused for any other, there being no tube to measure against.
    ##
    ## @item @qcode{'Mode'} @tab @qcode{'average'}, the default, which returns
    ## one number over the observations of every fold asked for, or
    ## @qcode{'individual'}, which returns one number per fold.
    ##
    ## @item @qcode{'Folds'} @tab A vector of fold indices to restrict the
    ## loss to.  It defaults to every fold.
    ## @end multitable
    ##
    ## @seealso{RegressionPartitionedModel, kfoldPredict}
    ## @end deftypefn
    function L = kfoldLoss (this, varargin)

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionPartitionedModel.kfoldLoss: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      ## Defaults, then the optional pairs
      LossFun = 'mse';
      Mode    = 'average';
      Folds   = 1:this.KFold;
      while (numel (varargin) > 0)
        if (! (ischar (varargin{1}) && isrow (varargin{1})))
          error (strcat ("RegressionPartitionedModel.kfoldLoss: parameter", ...
                         " name must be a character vector."));
        endif
        switch (tolower (varargin{1}))

          case 'lossfun'
            LossFun = varargin{2};
            if (! (is_function_handle (LossFun) ||
                   (ischar (LossFun) && isrow (LossFun))))
              error (strcat ("RegressionPartitionedModel.kfoldLoss:", ...
                             " 'LossFun' must be a character vector or", ...
                             " a function handle."));
            endif
            if (ischar (LossFun) && ! any (strcmpi (LossFun, ...
                                           {'mse', 'epsiloninsensitive'})))
              error (strcat ("RegressionPartitionedModel.kfoldLoss:", ...
                             " unsupported 'LossFun' value."));
            endif

          case 'mode'
            Mode = varargin{2};
            if (! (ischar (Mode) && isrow (Mode) &&
                   any (strcmpi (Mode, {'average', 'individual'}))))
              error (strcat ("RegressionPartitionedModel.kfoldLoss:", ...
                             " 'Mode' must be either 'average' or", ...
                             " 'individual'."));
            endif

          case 'folds'
            Folds = varargin{2};
            if (! (isnumeric (Folds) && isvector (Folds)
                   && all (Folds == fix (Folds))
                   && all (Folds >= 1) && all (Folds <= this.KFold)))
              error (strcat ("RegressionPartitionedModel.kfoldLoss:", ...
                             " 'Folds' must be a vector of fold indices", ...
                             " between 1 and KFold."));
            endif

          otherwise
            error (strcat ("RegressionPartitionedModel.kfoldLoss: invalid", ...
                           " parameter name in optional paired arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## The insensitive tube is a property of a support vector model, so any
      ## other cross-validated model has nothing to measure against.  Answering
      ## anyway would report a number for a quantity the model does not have.
      if (ischar (LossFun) && strcmpi (LossFun, 'epsiloninsensitive')
          && ! strcmp (this.CrossValidatedModel, 'SVM'))
        error (strcat ("RegressionPartitionedModel.kfoldLoss: the", ...
                       " 'epsiloninsensitive' loss applies to a", ...
                       " RegressionSVM model only."));
      endif

      yFit = kfoldPredict (this);

      if (strcmpi (Mode, 'individual'))
        L = nan (numel (Folds), 1);
        for i = 1:numel (Folds)
          idx = test (this.Partition, Folds(i));
          L(i) = foldLoss_ (this, idx, yFit, LossFun, Folds(i));
        endfor
      else
        idx = false (this.NumObservations, 1);
        for i = 1:numel (Folds)
          idx = idx | test (this.Partition, Folds(i));
        endfor
        L = foldLoss_ (this, idx, yFit, LossFun, Folds(1));
      endif

    endfunction

  endmethods

  methods (Access = private)

    ## Loss over the observations selected by IDX, weighted by W normalized
    ## over that selection so a subset is an average rather than a sum.
    function L = foldLoss_ (this, idx, yFit, LossFun, fold)
      if (! any (idx))
        L = NaN;
        return;
      endif
      y = this.Y(idx);
      f = yFit(idx);
      w = this.W(idx);
      w = w(:) / sum (w);
      if (is_function_handle (LossFun))
        L = LossFun (y(:), f(:), w);
        if (! (isnumeric (L) && isscalar (L)))
          error (strcat ("RegressionPartitionedModel.kfoldLoss: 'LossFun'", ...
                         " must return a numeric scalar."));
        endif
      elseif (strcmpi (LossFun, 'epsiloninsensitive'))
        ## Every fold was fitted with the same Epsilon, so any trained model
        ## reports it; the one belonging to this fold is used for clarity.
        eps_ = this.Trained{fold}.Epsilon;
        L = sum (w .* max (0, abs (y(:) - f(:)) - eps_));
      else
        L = sum (w .* (y(:) - f(:)) .^ 2);
      endif
    endfunction

  endmethods

endclassdef

## crossval builds one compact model per fold, over the observations used.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (40, 2);
%! Y = X(:,1) - 2 * X(:,2);
%! Mdl = fitrnet (X, Y, 'IterationLimit', 30);
%! CVMdl = crossval (Mdl, 'KFold', 4);
%! assert_equal (class (CVMdl), 'RegressionPartitionedModel');
%! assert_equal (CVMdl.KFold, 4);
%! assert_equal (numel (CVMdl.Trained), 4);
%! assert_equal (class (CVMdl.Trained{1}), 'CompactRegressionNeuralNetwork');
%! assert_equal (CVMdl.CrossValidatedModel, 'NeuralNetwork');
%! assert_equal (CVMdl.NumObservations, 40);
%! assert_equal (CVMdl.ResponseName, 'Y');
%! assert_equal (class (CVMdl.Partition), 'cvpartition');

## The same for a support vector model.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (40, 2);
%! Y = X(:,1) + X(:,2);
%! CVMdl = crossval (fitrsvm (X, Y), 'KFold', 4);
%! assert_equal (class (CVMdl.Trained{1}), 'CompactRegressionSVM');
%! assert_equal (CVMdl.CrossValidatedModel, 'SVM');

## The same for a generalized additive model, which carries no parameter
## struct, so that property keeps its default.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrgam (meas(1:20,1:3), meas(1:20,4)), 'KFold', 4);
%! assert_equal (class (CVMdl.Trained{1}), 'CompactRegressionGAM');
%! assert_equal (CVMdl.CrossValidatedModel, 'GAM');
%! assert_equal (CVMdl.NumObservations, 20);
%! assert_equal (CVMdl.ModelParameters, []);

%!error<RegressionPartitionedModel: unsupported model type.> ...
%! RegressionPartitionedModel (1, cvpartition (10, 'KFold', 2))

## Every observation is answered for by the fold that held it out.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (30, 2);
%! Y = X(:,1) * 2;
%! CVMdl = crossval (fitrsvm (X, Y), 'KFold', 5);
%! yFit = kfoldPredict (CVMdl);
%! assert_equal (size (yFit), [30, 1]);
%! assert_equal (any (isnan (yFit)), false);

## A holdout partition tests only its test set, and leaves the rest NaN.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (40, 2);
%! Y = X(:,1);
%! CVMdl = crossval (fitrsvm (X, Y), 'Holdout', 0.25);
%! yFit = kfoldPredict (CVMdl);
%! assert_equal (CVMdl.KFold, 1);
%! assert_equal (sum (! isnan (yFit)), sum (test (CVMdl.Partition, 1)));
%! assert_equal (sum (isnan (yFit)), sum (training (CVMdl.Partition, 1)));

## Leaveout holds out one observation at a time.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (12, 2);
%! Y = X(:,1) + 1;
%! CVMdl = crossval (fitrsvm (X, Y), 'Leaveout', 'on');
%! assert_equal (CVMdl.KFold, 12);
%! assert_equal (any (isnan (kfoldPredict (CVMdl))), false);

## A cvpartition object is taken as given.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (30, 2);
%! Y = X(:,2);
%! cvp = cvpartition (30, 'KFold', 3);
%! CVMdl = crossval (fitrsvm (X, Y), 'CVPartition', cvp);
%! assert_equal (CVMdl.KFold, 3);
%! assert_equal (CVMdl.Partition.NumObservations, 30);

## kfoldLoss is the weighted mean squared error of the out-of-fold answers.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (40, 2);
%! Y = X(:,1) - X(:,2);
%! CVMdl = crossval (fitrsvm (X, Y), 'KFold', 4);
%! yFit = kfoldPredict (CVMdl);
%! assert_equal (kfoldLoss (CVMdl), mean ((Y - yFit) .^ 2), 1e-12);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', 'mse'), ...
%!               kfoldLoss (CVMdl), 1e-12);

## 'individual' returns one loss per fold, and each is that fold's own.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (40, 2);
%! Y = X(:,1) * 3;
%! CVMdl = crossval (fitrsvm (X, Y), 'KFold', 4);
%! yFit = kfoldPredict (CVMdl);
%! L = kfoldLoss (CVMdl, 'Mode', 'individual');
%! assert_equal (size (L), [4, 1]);
%! idx = test (CVMdl.Partition, 2);
%! assert_equal (L(2), mean ((Y(idx) - yFit(idx)) .^ 2), 1e-12);

## 'Folds' restricts the loss to the folds named.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (40, 2);
%! Y = X(:,2) - 1;
%! CVMdl = crossval (fitrsvm (X, Y), 'KFold', 4);
%! yFit = kfoldPredict (CVMdl);
%! idx = test (CVMdl.Partition, 1) | test (CVMdl.Partition, 3);
%! assert_equal (kfoldLoss (CVMdl, 'Folds', [1, 3]), ...
%!               mean ((Y(idx) - yFit(idx)) .^ 2), 1e-12);
%! assert_equal (numel (kfoldLoss (CVMdl, 'Folds', [2, 4], ...
%!                                 'Mode', 'individual')), 2);

## The epsilon-insensitive loss uses the tube the folds were fitted with.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (40, 2);
%! Y = X(:,1) * 2;
%! CVMdl = crossval (fitrsvm (X, Y, 'Epsilon', 0.4), 'KFold', 4);
%! yFit = kfoldPredict (CVMdl);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', 'epsiloninsensitive'), ...
%!               mean (max (0, abs (Y - yFit) - 0.4)), 1e-12);

## kfoldLoss takes a function handle of the response, the fit and the weights.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = randn (30, 2);
%! Y = X(:,1);
%! CVMdl = crossval (fitrsvm (X, Y), 'KFold', 3);
%! yFit = kfoldPredict (CVMdl);
%! f = @(y, yf, w) sum (w .* abs (y - yf));
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', f), ...
%!               mean (abs (Y - yFit)), 1e-12);

## Rows dropped for missing values are outside the partition entirely.
%!test
%! rand ('seed', 42); randn ('seed', 42);
%! X = [randn(20, 2); NaN, 1];
%! Y = [randn(20, 1); 3];
%! Mdl = fitrsvm (X, Y);
%! CVMdl = crossval (Mdl, 'KFold', 4);
%! assert_equal (CVMdl.NumObservations, 21);
%! assert_equal (rows (CVMdl.X), 21);
%! assert_equal (numel (kfoldPredict (CVMdl)), 21);

## Test input validation for the constructor
%!error<RegressionPartitionedModel: too few input arguments.> ...
%! RegressionPartitionedModel ()
%!error<RegressionPartitionedModel: too few input arguments.> ...
%! RegressionPartitionedModel (fitrsvm (ones (5, 2), [1; 2; 3; 4; 5]))
%!error<RegressionPartitionedModel: unsupported model type.> ...
%! RegressionPartitionedModel (5, cvpartition (5, 'KFold', 2))
%!error<RegressionPartitionedModel: unsupported model type.> ...
%! RegressionPartitionedModel (fitcnet (ones (4, 2), [1; 1; 2; 2]), ...
%!                             cvpartition (4, 'KFold', 2))
%!error<RegressionPartitionedModel: invalid 'cvpartition' object.> ...
%! RegressionPartitionedModel (fitrsvm (ones (5, 2), [1; 2; 3; 4; 5]), 5)
%!error<RegressionPartitionedModel: 'cvpartition' object must be defined over the 5 observations the model was trained on.> ...
%! RegressionPartitionedModel (fitrsvm (ones (5, 2), [1; 2; 3; 4; 5]), ...
%!                             cvpartition (9, 'KFold', 3))

## Test input validation for kfoldLoss
%!shared CVR
%! rand ('seed', 42); randn ('seed', 42);
%! CVR = crossval (fitrsvm (randn (20, 2), randn (20, 1)), 'KFold', 4);
%!error<RegressionPartitionedModel.kfoldLoss: Name-Value arguments must be in pairs.> ...
%! kfoldLoss (CVR, 'Mode')
%!error<RegressionPartitionedModel.kfoldLoss: parameter name must be a character vector.> ...
%! kfoldLoss (CVR, 5, 1)
%!error<RegressionPartitionedModel.kfoldLoss: 'LossFun' must be a character vector or a function handle.> ...
%! kfoldLoss (CVR, 'LossFun', 5)
%!error<RegressionPartitionedModel.kfoldLoss: unsupported 'LossFun' value.> ...
%! kfoldLoss (CVR, 'LossFun', 'mae')
%!error<RegressionPartitionedModel.kfoldLoss: 'LossFun' must return a numeric scalar.> ...
%! kfoldLoss (CVR, 'LossFun', @(y, yf, w) [1, 2])
%!error<RegressionPartitionedModel.kfoldLoss: 'Mode' must be either 'average' or 'individual'.> ...
%! kfoldLoss (CVR, 'Mode', 'nope')
%!error<RegressionPartitionedModel.kfoldLoss: 'Folds' must be a vector of fold indices between 1 and KFold.> ...
%! kfoldLoss (CVR, 'Folds', 0)
%!error<RegressionPartitionedModel.kfoldLoss: 'Folds' must be a vector of fold indices between 1 and KFold.> ...
%! kfoldLoss (CVR, 'Folds', 9)
%!error<RegressionPartitionedModel.kfoldLoss: invalid parameter name in optional paired arguments.> ...
%! kfoldLoss (CVR, 'Nope', 1)

## The insensitive tube belongs to a support vector model only.
%!error<RegressionPartitionedModel.kfoldLoss: the 'epsiloninsensitive' loss applies to a RegressionSVM model only.> ...
%! kfoldLoss (crossval (fitrnet (randn (20, 2), randn (20, 1), ...
%!                      'IterationLimit', 5), 'KFold', 4), ...
%!            'LossFun', 'epsiloninsensitive')

## BinEdges is an empty cell, and a cell rather than an empty matrix, matching
## the classification counterpart and MATLAB.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrsvm (meas(:,1:3), meas(:,4)), 'KFold', 3);
%! assert_equal (class (CVMdl.BinEdges), 'cell');
%! assert_equal (CVMdl.BinEdges, {});

## No fold carries the transform, whichever model was cross validated.  The
## parent applies it once to the assembled prediction, so a fold that kept one
## of its own would apply it twice.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrgam (meas(:,1:3), meas(:,4), ...
%!                            'ResponseTransform', 'exp'), 'KFold', 3);
%! assert_equal (CVMdl.ResponseTransform, 'exp');
%! assert_equal (CVMdl.Trained{1}.ResponseTransform, 'none');

%!test
%! load fisheriris
%! CVMdl = crossval (fitrnet (meas(:,1:3), meas(:,4), ...
%!                            'ResponseTransform', 'exp'), 'KFold', 3);
%! assert_equal (CVMdl.ResponseTransform, 'exp');
%! assert_equal (CVMdl.Trained{1}.ResponseTransform, 'none');

%!test
%! load fisheriris
%! CVMdl = crossval (fitrgp (meas(:,1:3), meas(:,4), ...
%!                           'ResponseTransform', 'exp'), 'KFold', 3);
%! assert_equal (CVMdl.ResponseTransform, 'exp');
%! assert_equal (CVMdl.Trained{1}.ResponseTransform, 'none');

## ResponseTransform is applied to the assembled predictions, not carried into
## the folds.  MathWorks documents it as the function for transforming the
## predicted response values; it used to be stored here and never read.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrsvm (meas(:,1:3), meas(:,4)), 'KFold', 3);
%! y0 = kfoldPredict (CVMdl);
%! CVMdl.ResponseTransform = @(x) x + 100;
%! y1 = kfoldPredict (CVMdl);
%! assert_equal (CVMdl.Trained{1}.ResponseTransform, 'none');
%! assert_equal (y1, y0 + 100, 1e-12);

## And it reaches kfoldLoss, which is computed from those predictions.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrsvm (meas(:,1:3), meas(:,4)), 'KFold', 3);
%! before = kfoldLoss (CVMdl);
%! CVMdl.ResponseTransform = @(x) x + 100;
%! assert (kfoldLoss (CVMdl) > before);

## 'none' is the identity, so the default transforms nothing.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrnet (meas(:,1:3), meas(:,4)), 'KFold', 3);
%! y0 = kfoldPredict (CVMdl);
%! CVMdl.ResponseTransform = 'none';
%! assert_equal (kfoldPredict (CVMdl), y0);

## foldLoss_ is a private helper and stays out of the method list.
%!test
%! assert_equal (any (strcmp (methods ("RegressionPartitionedModel"), ...
%!                           "foldLoss_")), false);
