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

  properties (GetAccess = public, SetAccess = protected)

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
    ## @deftp {RegressionPartitionedModel} {property} ModelParameters
    ##
    ## Parameters the folds were fitted with
    ##
    ## A structure holding the parameters the folds were fitted with, carried
    ## through from the learner that was cross validated, beside
    ## @qcode{NLearn}, the number of folds, and the @qcode{Version},
    ## @qcode{Method} and @qcode{Type} tags of this class.  The
    ## learner's own tags are replaced rather than kept, so a cross-validated
    ## SVM reports @qcode{Method} as @qcode{'PartitionedModel'} and not
    ## @qcode{'SVM'}.
    ##
    ## @strong{Deviation from MATLAB.}  MATLAB reports the parameter record of
    ## the cross-validation @emph{ensemble} here rather than of the learner,
    ## so it says nothing at all about how the folds were fitted: of its
    ## eighteen fields only the fold count, its partitioner and a fit template
    ## carry anything, and the rest are boosting settings left inert.  Nor can
    ## the parameters be reached through the folds, a compact model carrying
    ## none in MATLAB.  This class reports the fit instead, which is strictly
    ## more than MATLAB offers, and everything MATLAB's record does carry is
    ## published here as the @qcode{KFold}, @qcode{Partition}, @qcode{X},
    ## @qcode{Y}, @qcode{W} and @qcode{CrossValidatedModel} properties.
    ##
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
    ## @deftp {RegressionPartitionedModel} {property} IsStandardDeviationFit
    ##
    ## Whether the folds fitted a standard deviation model
    ##
    ## A logical scalar for a generalized additive model backing, taken from
    ## the model that was cross validated, and empty for every other backing.
    ##
    ## MATLAB carries this on @code{RegressionPartitionedGAM}, one of five
    ## per-learner partitioned classes this package deliberately does not have
    ## (see @code{crossval}).  With one class serving every backing the
    ## property has to be declared for all of them, so it is empty where it
    ## does not apply.  It is placed last rather than first, where MATLAB's
    ## subclass shows it, because that subclass also moves
    ## @qcode{ResponseTransform} to the end and no single order can match both
    ## of MATLAB's classes; matching the general one and appending is the only
    ## coherent choice.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    IsStandardDeviationFit = [];


    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedModel} {property} NumTrainedPerFold
    ##
    ## How many trees each fold fitted
    ##
    ## A scalar structure with fields @qcode{PredictorTrees} and
    ## @qcode{InteractionTrees}, each a row with one entry per fold, for a
    ## generalized additive model backing, and empty for every other.
    ##
    ## It reports what each fold actually fitted, which the budget in
    ## @qcode{ModelParameters} does not: a phase stops early when it can no
    ## longer improve the fit, and the folds need not stop at the same place.
    ##
    ## MATLAB carries this on its per-learner partitioned GAM classes, which
    ## this package deliberately does not have (see @code{crossval}), so like
    ## @qcode{IsStandardDeviationFit} it is declared here for every backing
    ## and left empty where it does not apply.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumTrainedPerFold = [];

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
      ## The learner's parameters, under this class's own tags.  MATLAB
      ## reports an EnsembleParams here instead and so says nothing about the
      ## fit; see the ModelParameters property for the deviation.
      this.ModelParameters = partitionedModelParams (Mdl, this.KFold, ...
                               'PartitionedModel', 'regression');
      if (any (strcmp (properties (Mdl), 'IsStandardDeviationFit')))
        this.IsStandardDeviationFit = Mdl.IsStandardDeviationFit;
      endif


      ## Switch Regression object types
      switch (this.CrossValidatedModel)

        case 'GAM'
          ## Knots, Order and DoF are three views of one parameterisation and
          ## the constructor accepts any two, recomputing the third, so only
          ## Knots and Order are passed on.
          args = {};
          ## Which parameters a fold takes depends on which engine fitted
          ## the parent: the two have disjoint argument surfaces and each
          ## refuses the other's, so the fold is refitted with its own.
          if (strcmp (Mdl.FitMethod, 'boostedtrees'))
            GAMparams = {'PredictorNames', 'ResponseName'};
          else
            GAMparams = {'PredictorNames', 'ResponseName', 'Formula', ...
                         'Knots', 'Order', 'Tol'};
          endif
          args = [args, {'FitMethod', Mdl.FitMethod}];
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
                  'IterationLimit', Mdl.IterationLimit, ...
                  'Standardize', stdz, ...
                  'ResponseName', Mdl.ResponseName, ...
                  'PredictorNames', Mdl.PredictorNames};
          ## As in ClassificationPartitionedModel: the folds are trained by
          ## the parent's solver, and a learning rate goes only with 'sgd'.
          if (strcmp (Mdl.Solver, 'LBFGS'))
            args = [args, {'Solver', 'lbfgs'}];
          else
            args = [args, {'Solver', 'sgd', 'LearningRate', Mdl.LearningRate}];
          endif
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
          ## The polynomial order is recorded for the polynomial kernel
          ## alone, so it is passed only when the parent carries one.
          korder = {};
          if (! isempty (p.KernelPolynomialOrder))
            korder = {'PolynomialOrder', p.KernelPolynomialOrder};
          endif
          args = {'SVMtype', p.SVMtype, ...
                  'KernelFunction', p.KernelFunction, ...
                  korder{:}, ...
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

      ## Gathered across the folds, which is the shape MATLAB reports: one
      ## structure carrying a row per phase rather than a structure per fold.
      if (strcmp (this.CrossValidatedModel, 'GAM'))
        pt = zeros (1, this.KFold);
        it = zeros (1, this.KFold);
        boosted = true;
        for k = 1:this.KFold
          nt = this.Trained{k}.NumTrainedTrees;
          ## A spline fit counts no trees, so there is nothing to report and
          ## the property stays empty rather than claiming zero of them.
          if (isempty (nt))
            boosted = false;
            break;
          endif
          pt(k) = nt.PredictorTrees;
          it(k) = nt.InteractionTrees;
        endfor
        if (boosted)
          this.NumTrainedPerFold = struct ('PredictorTrees', pt, ...
                                           'InteractionTrees', it);
        endif
      endif

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
    ## @deftypefn  {RegressionPartitionedModel} {@var{yFit} =} kfoldPredict (@var{obj})
    ## @deftypefnx {RegressionPartitionedModel} {[@var{yFit}, @var{ySD}, @var{yInt}] =} kfoldPredict (@var{obj})
    ## @deftypefnx {RegressionPartitionedModel} {[@dots{}] =} kfoldPredict (@dots{}, @qcode{'Alpha'}, @var{alpha})
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
    ## @code{[@var{yFit}, @var{ySD}, @var{yInt}] = kfoldPredict (@var{obj})}
    ## also returns the standard deviation @var{ySD} of each predicted
    ## response and the two-column matrix @var{yInt} of prediction intervals,
    ## each answered for by the fold that held the observation out.  A
    ## @qcode{RegressionGP} backing is the only one that fits the uncertainty
    ## its predictions carry, so any other raises here.  An untested
    ## observation is @code{NaN} in all three.
    ##
    ## @code{[@dots{}] = kfoldPredict (@dots{}, @qcode{'Alpha'}, @var{alpha})}
    ## sets the significance level of the prediction intervals, which default
    ## to 95 per cent at an @var{alpha} of 0.05.
    ##
    ## @var{ySD} does not follow @qcode{ResponseTransform} and the other two
    ## outputs do, the same rule @code{RegressionGP.predict} applies: a
    ## predicted response and an interval endpoint are on the response scale
    ## and a standard deviation is not.
    ##
    ## @seealso{RegressionPartitionedModel, kfoldLoss}
    ## @end deftypefn
    function [yFit, ySD, yInt] = kfoldPredict (this, varargin)

      ## Only the GP fits the uncertainty around its predictions, so it is the
      ## only backing whose folds can be asked for one.  MATLAB refuses the
      ## network and the SVM outright and refuses the GAM unless it was fitted
      ## with 'FitStandardDeviation', an option this package does not offer, so
      ## every backing but the GP is refused here under one message.
      is_gp = strcmp (this.CrossValidatedModel, 'GP');
      if (nargout > 1 && ! is_gp)
        error (strcat ("RegressionPartitionedModel.kfoldPredict: a", ...
                       " standard deviation and a prediction interval are", ...
                       " only available for a cross-validated",  ...
                       " RegressionGP."));
      endif

      ## Validated here rather than in the fold's predict, so the message
      ## names the method the user called.
      if (numel (varargin) > 0 && ! is_gp)
        error (strcat ("RegressionPartitionedModel.kfoldPredict:", ...
                       " optional arguments are only accepted for a", ...
                       " cross-validated RegressionGP."));
      endif
      CIAlpha = 0.05;
      while (numel (varargin) > 0)
        if (numel (varargin) < 2)
          error (strcat ("RegressionPartitionedModel.kfoldPredict:", ...
                         " optional arguments must be given in Name-Value", ...
                         " pairs."));
        endif
        switch (lower (varargin{1}))
          case 'alpha'
            CIAlpha = varargin{2};
            if (! (isnumeric (CIAlpha) && isscalar (CIAlpha) && ...
                   CIAlpha >= 0 && CIAlpha <= 1))
              error (strcat ("RegressionPartitionedModel.kfoldPredict:", ...
                             " 'Alpha' must be a scalar between 0 and 1."));
            endif
          otherwise
            error (strcat ("RegressionPartitionedModel.kfoldPredict:", ...
                           " invalid NAME in optional pairs of arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      yFit = nan (this.NumObservations, 1);
      if (nargout > 1)
        ySD = nan (this.NumObservations, 1);
        yInt = nan (this.NumObservations, 2);
      endif

      for k = 1:this.KFold
        testIdx = test (this.Partition, k);
        if (! any (testIdx))
          continue;
        endif
        if (nargout > 2)
          [yFit(testIdx), ySD(testIdx), yInt(testIdx,:)] = ...
                predict (this.Trained{k}, this.X(testIdx, :), 'Alpha', CIAlpha);
        elseif (nargout > 1)
          [yFit(testIdx), ySD(testIdx)] = ...
                predict (this.Trained{k}, this.X(testIdx, :));
        else
          yFit(testIdx) = predict (this.Trained{k}, this.X(testIdx, :));
        endif
      endfor

      ## As on the classification side, the folds never carry the transform and
      ## it is applied once to the assembled predictions.  MathWorks documents
      ## ResponseTransform as the function for transforming the predicted
      ## response values and documents assigning it by dot notation, so this is
      ## what the property is for; it used to be stored and never read.
      yFit = this.RTfun (yFit);
      if (nargout > 2)
        yInt = this.RTfun (yInt);
      endif

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

    ## -*- texinfo -*-
    ## @deftypefn {RegressionPartitionedModel} {@var{vals} =} kfoldfun (@var{obj}, @var{fun})
    ##
    ## Apply a function to each fold of a cross-validated model.
    ##
    ## @code{@var{vals} = kfoldfun (@var{obj}, @var{fun})} calls @var{fun} once
    ## per fold and returns a @math{K*M} numeric matrix whose row @math{k} is
    ## what @var{fun} returned for fold @math{k}.
    ##
    ## @var{fun} is a function handle taking seven inputs and returning a
    ## numeric vector of the same length every time it is called:
    ##
    ## @example
    ## @var{testvals} = @var{fun} (@var{M}, @var{Xtrain}, @var{Ytrain}, @var{Wtrain}, @dots{}
    ##                @var{Xtest}, @var{Ytest}, @var{Wtest})
    ## @end example
    ##
    ## @var{M} is the model the fold was fitted with, taken from
    ## @code{@var{obj}.Trained@{k@}}; @var{Xtrain}, @var{Ytrain} and
    ## @var{Wtrain} are the predictors, response and weights of the
    ## observations that fold was trained on, and @var{Xtest}, @var{Ytest} and
    ## @var{Wtest} those of the observations it held out.
    ##
    ## @seealso{RegressionPartitionedModel, kfoldPredict, kfoldLoss}
    ## @end deftypefn

    function vals = kfoldfun (this, fun)

      if (nargin < 2)
        error ("RegressionPartitionedModel.kfoldfun: too few input arguments.");
      endif
      if (! is_function_handle (fun))
        error ("RegressionPartitionedModel.kfoldfun: FUN must be a function handle.");
      endif

      vals = [];
      for k = 1:this.KFold

        trIdx = training (this.Partition, k);
        teIdx = test (this.Partition, k);

        tv = fun (this.Trained{k}, this.X(trIdx,:), this.Y(trIdx,:), ...
                  this.W(trIdx), this.X(teIdx,:), this.Y(teIdx,:), ...
                  this.W(teIdx));

        ## The returned values become one row, so they have to be numeric and
        ## of one length: a fold answering with a different width could not be
        ## stacked with the others, and finding that out at the concatenation
        ## would name neither the fold nor the reason.
        if (! ((isnumeric (tv) || islogical (tv)) && isvector (tv)))
          error (strcat ("RegressionPartitionedModel.kfoldfun: FUN must", ...
                         " return a numeric vector; fold %d returned a %s."), ...
                 k, class (tv));
        endif
        tv = tv(:)';
        if (! isempty (vals) && numel (tv) != columns (vals))
          error (strcat ("RegressionPartitionedModel.kfoldfun: FUN must", ...
                         " return the same number of values for every fold;", ...
                         " fold %d returned %d where the first returned", ...
                         " %d."), k, numel (tv), columns (vals));
        endif
        vals = [vals; tv];

      endfor

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
%! ## The learner's thirteen fields come through, with NLearn added and the
%! ## three tags reissued for this class, so fourteen in all.
%! assert_equal (isstruct (CVMdl.ModelParameters), true);
%! assert_equal (numfields (CVMdl.ModelParameters), 14);
%! assert_equal (isfield (CVMdl.ModelParameters, 'NumTreesPerPredictor'), true);
%! assert_equal (CVMdl.ModelParameters.Method, 'PartitionedModel');

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

## A GP backing answers for the uncertainty of its out-of-fold predictions.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrgp (meas(:,1:3), meas(:,4)), 'KFold', 5);
%! [yFit, ySD, yInt] = kfoldPredict (CVMdl);
%! assert_equal (size (yFit), [150, 1]);
%! assert_equal (size (ySD), [150, 1]);
%! assert_equal (size (yInt), [150, 2]);
%! assert_equal (all (ySD > 0), true);
%! assert_equal (yFit, kfoldPredict (CVMdl), 1e-12);

## The interval is the prediction plus and minus a normal quantile of the
## standard deviation, and 'Alpha' sets which quantile.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrgp (meas(:,1:3), meas(:,4)), 'KFold', 5);
%! [yFit, ySD, yInt] = kfoldPredict (CVMdl);
%! z = norminv (0.975);
%! assert_equal (yInt, [yFit - z * ySD, yFit + z * ySD], 1e-12);
%! [~, ~, yInt90] = kfoldPredict (CVMdl, 'Alpha', 0.10);
%! z90 = norminv (0.95);
%! assert_equal (yInt90, [yFit - z90 * ySD, yFit + z90 * ySD], 1e-12);
%! assert_equal (all (diff (yInt90, 1, 2) < diff (yInt, 1, 2)), true);

## Each observation is answered for by the fold that held it out, in all three
## outputs, and one no fold tests is NaN in all three.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrgp (meas(:,1:3), meas(:,4)), 'KFold', 5);
%! [yFit, ySD, yInt] = kfoldPredict (CVMdl);
%! idx = test (CVMdl.Partition, 2);
%! [p, s, i] = predict (CVMdl.Trained{2}, meas(idx,1:3));
%! assert_equal (p, yFit(idx), 1e-12);
%! assert_equal (s, ySD(idx), 1e-12);
%! assert_equal (i, yInt(idx,:), 1e-12);

%!test
%! load fisheriris
%! CVMdl = crossval (fitrgp (meas(:,1:3), meas(:,4)), 'Holdout', 0.3);
%! [yFit, ySD, yInt] = kfoldPredict (CVMdl);
%! untested = isnan (yFit);
%! assert_equal (any (untested), true);
%! assert_equal (isnan (ySD), untested);
%! assert_equal (isnan (yInt), [untested, untested]);

## ResponseTransform reaches the two outputs on the response scale and not the
## standard deviation, the rule RegressionGP.predict applies.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrgp (meas(:,1:3), meas(:,4)), 'KFold', 5);
%! [yFit, ySD, yInt] = kfoldPredict (CVMdl);
%! CVMdl.ResponseTransform = @(x) 2 * x;
%! [yFit2, ySD2, yInt2] = kfoldPredict (CVMdl);
%! assert_equal (yFit2, 2 * yFit, 1e-12);
%! assert_equal (ySD2, ySD, 1e-12);
%! assert_equal (yInt2, 2 * yInt, 1e-12);

## Measured against R2024a over the same folds, an explicit CustomPartition
## naming them in both engines so the assembled values are comparable at all.
%!test
%! load fisheriris
%! cvp = cvpartition ('CustomPartition', repmat ((1:5)', 30, 1));
%! CVMdl = crossval (fitrgp (meas(:,1:3), meas(:,4)), 'CVPartition', cvp);
%! [yFit, ySD, yInt] = kfoldPredict (CVMdl);
%! assert_equal (yFit(1), 0.21819014717001, 1e-9);
%! assert_equal (ySD(1), 0.17986522125163, 1e-9);
%! assert_equal (yInt(1,:), [-0.13433920855451, 0.57071950289454], 1e-9);
%! assert_equal (sum (abs (yFit)), 179.61497470708, 1e-6);
%! assert_equal (sum (abs (ySD)), 27.747043040728, 1e-6);
%! [~, ~, yInt90] = kfoldPredict (CVMdl, 'Alpha', 0.10);
%! assert_equal (yInt90(1,:), [-0.077661814368161, 0.51404210870819], 1e-9);

## No other backing fits the uncertainty around its predictions, so none of
## them may be asked for it.
%!error<RegressionPartitionedModel.kfoldPredict: a standard deviation and a prediction interval are only available for a cross-validated RegressionGP.> ...
%! load fisheriris; ...
%! CVMdl = crossval (fitrsvm (meas(:,1:3), meas(:,4)), 'KFold', 3); ...
%! [yFit, ySD] = kfoldPredict (CVMdl);
%!error<RegressionPartitionedModel.kfoldPredict: a standard deviation and a prediction interval are only available for a cross-validated RegressionGP.> ...
%! load fisheriris; ...
%! CVMdl = crossval (fitrgam (meas(:,1:3), meas(:,4)), 'KFold', 3); ...
%! [yFit, ySD, yInt] = kfoldPredict (CVMdl);
%!error<RegressionPartitionedModel.kfoldPredict: optional arguments are only accepted for a cross-validated RegressionGP.> ...
%! load fisheriris; ...
%! CVMdl = crossval (fitrsvm (meas(:,1:3), meas(:,4)), 'KFold', 3); ...
%! kfoldPredict (CVMdl, 'Alpha', 0.1);
%!error<RegressionPartitionedModel.kfoldPredict: optional arguments must be given in Name-Value pairs.> ...
%! load fisheriris; ...
%! CVMdl = crossval (fitrgp (meas(:,1:3), meas(:,4)), 'KFold', 3); ...
%! kfoldPredict (CVMdl, 'Alpha');
%!error<RegressionPartitionedModel.kfoldPredict: 'Alpha' must be a scalar between 0 and 1.> ...
%! load fisheriris; ...
%! CVMdl = crossval (fitrgp (meas(:,1:3), meas(:,4)), 'KFold', 3); ...
%! kfoldPredict (CVMdl, 'Alpha', 2);
%!error<RegressionPartitionedModel.kfoldPredict: invalid NAME in optional pairs of arguments.> ...
%! load fisheriris; ...
%! CVMdl = crossval (fitrgp (meas(:,1:3), meas(:,4)), 'KFold', 3); ...
%! kfoldPredict (CVMdl, 'Bogus', 1);

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

## A fold is refitted with the argument set of the engine that fitted the
## parent, so a tree-fitted parent must not hand its folds Knots and Order.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,2:4), meas(:,1), 'FitMethod', 'boostedtrees');
%! CVMdl = crossval (Mdl, 'KFold', 3);
%! assert_equal (class (CVMdl), 'RegressionPartitionedModel');
%! assert_equal (numel (CVMdl.Trained), 3);
%! assert_equal (CVMdl.Trained{1}.FitMethod, 'boostedtrees');
%! assert_equal (numel (kfoldPredict (CVMdl)), rows (meas));

## A spline-fitted parent is unaffected: its folds take the spline parameters
## and are fitted by the same engine.  They do not report the parameters
## themselves, a compact model keeping no record of its fitting.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(1:60,1:3), meas(1:60,4), 'FitMethod', 'splines', ...
%!                'Knots', 4);
%! CVMdl = crossval (Mdl, 'KFold', 3);
%! assert_equal (CVMdl.Trained{1}.FitMethod, 'splines');
%! assert_equal (isfield (CVMdl.Trained{1}.BaseModel, 'Intercept'), true);

%!test
%! ## kfoldfun hands over seven arguments, the fold's model first.
%! load fisheriris
%! CV = crossval (fitrgam (meas(:,2:4), meas(:,1)), "KFold", 3);
%! seen = kfoldfun (CV, @(M, Xtr, Ytr, Wtr, Xte, Yte, Wte) ...
%!                      [rows(Xtr), rows(Yte), columns(Xtr)]);
%! assert_equal (size (seen), [3, 3]);
%! assert_equal (seen(:,1) + seen(:,2), repmat (150, 3, 1));
%! assert_equal (seen(:,3), repmat (3, 3, 1));

%!test
%! ## The use it exists for: a held-out mean squared error per fold.
%! load fisheriris
%! CV = crossval (fitrgam (meas(:,2:4), meas(:,1)), "KFold", 3);
%! f = @(M, Xtr, Ytr, Wtr, Xte, Yte, Wte) mean ((predict (M, Xte) - Yte) .^ 2);
%! mse = kfoldfun (CV, f);
%! assert_equal (size (mse), [3, 1]);
%! assert_equal (all (mse > 0), true);

%!error<RegressionPartitionedModel.kfoldfun: too few input arguments.> ...
%! kfoldfun (crossval (fitrgam (ones (8, 2), (1:8)'), "KFold", 2))
%!error<RegressionPartitionedModel.kfoldfun: FUN must be a function handle.> ...
%! kfoldfun (crossval (fitrgam (ones (8, 2), (1:8)'), "KFold", 2), "nope")

%!test
%! ## The property order is MATLAB's, measured on R2024a.
%! load fisheriris
%! CVMdl = crossval (fitrsvm (meas(:,2:4), meas(:,1)), "KFold", 3);
%! assert_equal (sort (properties (CVMdl)), ...
%!               sort ({'ResponseTransform'; 'CrossValidatedModel'; ...
%!                      'PredictorNames'; 'CategoricalPredictors'; ...
%!                      'ResponseName'; 'NumObservations'; 'X'; 'Y'; 'W'; ...
%!                      'ModelParameters'; 'Trained'; 'KFold'; 'Partition'; ...
%!                      'BinEdges'; 'IsStandardDeviationFit'; ...
%!                      'NumTrainedPerFold'}));

%!test
%! ## IsStandardDeviationFit comes from the model for a GAM backing and is
%! ## empty for every other, one class serving all of them.
%! load fisheriris
%! CVg = crossval (fitrgam (meas(:,2:4), meas(:,1)), "KFold", 3);
%! assert_equal (islogical (CVg.IsStandardDeviationFit), true);
%! assert_equal (CVg.IsStandardDeviationFit, false);
%! CVs = crossval (fitrsvm (meas(:,2:4), meas(:,1)), "KFold", 3);
%! assert_equal (isempty (CVs.IsStandardDeviationFit), true);

%!test
%! ## NumTrainedPerFold reports what each fold actually fitted.  On this
%! ## fixture every fold uses its whole budget, which is what R2024a reports
%! ## for it as well.
%! load fisheriris
%! CVMdl = crossval (fitrgam (meas(:,2:4), meas(:,1)), "KFold", 3);
%! n = CVMdl.NumTrainedPerFold;
%! assert_equal (sort (fieldnames (n)), {"InteractionTrees"; "PredictorTrees"});
%! assert_equal (n.PredictorTrees, [300, 300, 300]);
%! assert_equal (n.InteractionTrees, [0, 0, 0]);

%!test
%! ## Empty for a backing that fits no trees.
%! load fisheriris
%! CVMdl = crossval (fitrsvm (meas(:,2:4), meas(:,1)), "KFold", 3);
%! assert_equal (isempty (CVMdl.NumTrainedPerFold), true);

## A polynomial-kernel SVM refits its folds through the recorded order, which
## the parent carries only under that kernel.
%!test
%! load fisheriris
%! Mdl = fitrsvm (meas(:,2:4), meas(:,1), 'KernelFunction', 'polynomial', ...
%!                'PolynomialOrder', 2);
%! CVMdl = crossval (Mdl, 'KFold', 3);
%! assert_equal (CVMdl.ModelParameters.KernelPolynomialOrder, 2);
%! assert_equal (numel (CVMdl.Trained), 3);

## ModelParameters carries the learner's parameters under this class's own
## tags.
%!test
%! load fisheriris
%! CVMdl = crossval (fitrsvm (meas(:,2:4), meas(:,1)), 'KFold', 3);
%! MP = CVMdl.ModelParameters;
%! assert_equal (MP.Method, 'PartitionedModel');
%! assert_equal (MP.Type, 'regression');
%! assert_equal (MP.Version, 1);
%! assert_equal (MP.NLearn, 3);
%! assert_equal (MP.SVMtype, 'eps_svr');

%!test
%! load fisheriris
%! MP = crossval (fitrgp (meas(:,2:4), meas(:,1)), 'KFold', 3).ModelParameters;
%! assert_equal (MP.FitMethod, 'Exact');
%! assert_equal (MP.Method, 'PartitionedModel');

## Every documented response transform reaches the response that is reported.
%!test
%! load fisheriris
%! Mdl = crossval (fitrsvm (meas(:,2:4), meas(:,1)), 'KFold', 3);
%! Mdl.ResponseTransform = 'none';
%! raw = kfoldPredict (Mdl);
%! T = {'identity', @(x) x; 'exp', @(x) exp (x); 'log', @(x) log (x)};
%! for i = 1:rows (T)
%!   Mdl.ResponseTransform = T{i,1};
%!   yhat = kfoldPredict (Mdl);
%!   assert_equal (yhat, T{i,2}(raw), 1e-12);
%! endfor

## A function handle is taken as given and applied to the response.
%!test
%! load fisheriris
%! Mdl = crossval (fitrsvm (meas(:,2:4), meas(:,1)), 'KFold', 3);
%! Mdl.ResponseTransform = 'none';
%! raw = kfoldPredict (Mdl);
%! Mdl.ResponseTransform = @(x) x .^ 2;
%! yhat = kfoldPredict (Mdl);
%! assert_equal (yhat, raw .^ 2, 1e-12);
