## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftp {statistics} RegressionKernel
##
## Gaussian kernel regression model for large data.
##
## A @qcode{RegressionKernel} object maps the predictors into a randomized
## feature space whose inner product approximates a Gaussian kernel, and then
## fits a linear model there.  A kernel regression is therefore as nonlinear
## as a support vector machine with a Gaussian kernel, while costing what a
## linear fit costs: nothing of size @math{NxN} is ever formed.
##
## The expansion is the random Fourier basis of Rahimi and Recht, drawn once
## when the model is fitted and kept with it, so @code{predict} maps new data
## through the same basis.  MATLAB approximates the same kernel by the
## Fastfood construction, which reaches the same distribution more cheaply;
## the two are interchangeable in distribution but not draw by draw, and the
## draws come from different generators in any case, so the predictions of a
## model fitted here and one fitted in MATLAB differ even from the same seed.
## What does not differ is what they estimate.
##
## Like @qcode{RegressionLinear} the object holds no copy of the training
## data.  It does hold the basis and the coefficients, so it is bounded by
## the number of expansion dimensions rather than by the number of
## observations.
##
## Create a @qcode{RegressionKernel} object with @code{fitrkernel}.
##
## @seealso{fitrkernel, RegressionLinear, RegressionSVM}
## @end deftp

classdef RegressionKernel

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} Epsilon
    ##
    ## Half the width of the epsilon-insensitive band
    ##
    ## A nonnegative scalar for a support vector machine, and empty for a
    ## least squares fit, which has no such band.  It defaults to the
    ## interquartile range of the response over 13.49, an estimate of its
    ## standard deviation, or to @qcode{0.1} when that range is zero.  This
    ## property is read-only.
    ##
    ## @end deftp
    Epsilon                = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} BoxConstraint
    ##
    ## Box constraint of the support vector machine
    ##
    ## A positive scalar.  It is the reciprocal of the product of
    ## @qcode{Lambda} and the number of observations, so setting either of
    ## the two in the constructor fixes the other, and giving both is an
    ## error.  This property is read-only.
    ##
    ## @end deftp
    BoxConstraint          = 1;

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} ResponseTransform
    ##
    ## Transformation applied to the predicted response
    ##
    ## A character vector, or the text of the function handle that was
    ## supplied.  Assigning to it accepts either.
    ##
    ## @end deftp
    ResponseTransform      = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors with one name per column of the
    ## training data, defaulting to @qcode{'x1'}, @qcode{'x2'} and so on.
    ## This property is read-only.
    ##
    ## @end deftp
    PredictorNames         = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A row vector of column indices, empty when every predictor is
    ## numeric.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors  = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} ResponseName
    ##
    ## Name of the response
    ##
    ## A character vector, defaulting to @qcode{'Y'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    ResponseName           = 'Y';

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the fit saw them
    ##
    ## A cell array of character vectors.  These name the original
    ## predictors, not the expansion dimensions, which have no names.  This
    ## property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} NumExpansionDimensions
    ##
    ## Number of dimensions of the expanded space
    ##
    ## A positive integer scalar.  It defaults to
    ## @code{2 .^ ceil (min (log2 (@var{p}) + 5, 15))} for @var{p}
    ## predictors, so four predictors give 128 dimensions.  More dimensions
    ## approximate the kernel more closely and cost proportionally more.
    ## This property is read-only.
    ##
    ## @end deftp
    NumExpansionDimensions = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} FittedLoss
    ##
    ## Loss function the fit minimized
    ##
    ## @qcode{'epsiloninsensitive'} for a support vector machine and
    ## @qcode{'mse'} for a least squares fit.  This property is read-only.
    ##
    ## @end deftp
    FittedLoss             = 'epsiloninsensitive';

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} Lambda
    ##
    ## Regularization strength
    ##
    ## A nonnegative scalar, the reciprocal of the product of
    ## @qcode{BoxConstraint} and the number of observations.  This property
    ## is read-only.
    ##
    ## @end deftp
    Lambda                 = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} ModelParameters
    ##
    ## Fitting options, as they were given
    ##
    ## A structure holding every parameter of the fit, with the
    ## @qcode{'auto'} values as they were given rather than as they were
    ## resolved.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters        = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} Regularization
    ##
    ## Penalty on the coefficients
    ##
    ## Always @qcode{'ridge (L2)'}: a kernel model fits in the expanded
    ## space, where a lasso penalty has nothing to select.  This property is
    ## read-only.
    ##
    ## @end deftp
    Regularization         = 'ridge (L2)';

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} KernelScale
    ##
    ## Scale of the Gaussian kernel
    ##
    ## A positive scalar dividing every predictor before the expansion, so a
    ## larger scale makes the kernel wider and the fit smoother.  This
    ## property is read-only.
    ##
    ## @end deftp
    KernelScale            = 1;

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} Learner
    ##
    ## Linear model fitted in the expanded space
    ##
    ## Either @qcode{'svm'} or @qcode{'leastsquares'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Learner                = 'svm';

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} Mu
    ##
    ## Predictor means used to standardize
    ##
    ## A row vector with one element per predictor, or empty when the model
    ## was fitted without standardizing.  This property is read-only.
    ##
    ## @end deftp
    Mu                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionKernel} {property} Sigma
    ##
    ## Predictor standard deviations used to standardize
    ##
    ## A row vector with one element per predictor, or empty when the model
    ## was fitted without standardizing.  This property is read-only.
    ##
    ## @end deftp
    Sigma                  = [];

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## The callable behind ResponseTransform.
    RTfun                  = @(y) y;

    ## The random basis and the coefficients fitted in the space it spans.
    ## MATLAB keeps all three out of sight, so they are hidden here too, but
    ## they are the model: without the basis the coefficients index nothing.
    Basis_                 = [];
    Beta_                  = [];
    Bias_                  = [];

    ## Number of predictors and of observations the fit saw, neither of them
    ## a property of MATLAB's class.
    NumPredictors_         = [];
    NumObservations_       = [];

    ## What the fit reported, so that fitrkernel can hand it back as its
    ## second output.
    FitInfo_               = [];

  endproperties

  methods (Access = public)

    ## Custom setter, so that assigning a name or a handle updates both the
    ## text the property reports and the callable predict uses.
    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = ...
                     parseResponseTransform (val, 'RegressionKernel');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionKernel} {@var{obj} =} RegressionKernel (@var{X}, @var{Y})
    ## @deftypefnx {RegressionKernel} {@var{obj} =} RegressionKernel (@dots{}, @var{name}, @var{value})
    ##
    ## Fit a Gaussian kernel regression model.
    ##
    ## @code{@var{obj} = RegressionKernel (@var{X}, @var{Y})} fits a support
    ## vector machine in a randomized Gaussian kernel space to the @math{NxP}
    ## predictor matrix @var{X} and the @math{Nx1} continuous response
    ## @var{Y}.
    ##
    ## @code{@var{obj} = RegressionKernel (@dots{}, @var{name}, @var{value})}
    ## takes the following @qcode{Name-Value} pairs.
    ##
    ## @multitable @columnfractions 0.28 0.72
    ## @headitem Name @tab Value
    ##
    ## @item @qcode{'Learner'} @tab @qcode{'svm'}, the default, or
    ## @qcode{'leastsquares'}.
    ##
    ## @item @qcode{'Epsilon'} @tab Half the width of the insensitive band,
    ## a nonnegative scalar or @qcode{'auto'}, which is the interquartile
    ## range of @var{Y} over 13.49.  It applies to a support vector machine
    ## alone.
    ##
    ## @item @qcode{'NumExpansionDimensions'} @tab @qcode{'auto'}, the
    ## default, or a positive integer.
    ##
    ## @item @qcode{'KernelScale'} @tab @qcode{1} by default, a positive
    ## scalar, or @qcode{'auto'}, which takes the median distance between the
    ## observations.
    ##
    ## @item @qcode{'Lambda'} @tab @qcode{'auto'}, the default, which is the
    ## reciprocal of the number of observations, or a nonnegative scalar.  It
    ## cannot be given beside @qcode{'BoxConstraint'}.
    ##
    ## @item @qcode{'BoxConstraint'} @tab A positive scalar, @qcode{1} by
    ## default.  It applies to a support vector machine alone.
    ##
    ## @item @qcode{'Standardize'} @tab Whether to centre and scale the
    ## predictors, false by default.
    ##
    ## @item @qcode{'BetaTolerance'} @tab Relative tolerance on the
    ## coefficients, @qcode{1e-4} by default.
    ##
    ## @item @qcode{'GradientTolerance'} @tab Absolute tolerance on the
    ## gradient's infinity norm, @qcode{1e-6} by default.
    ##
    ## @item @qcode{'IterationLimit'} @tab Largest number of iterations,
    ## @qcode{1000} by default.
    ##
    ## @item @qcode{'HessianHistorySize'} @tab Number of curvature pairs the
    ## solver keeps, @qcode{15} by default.
    ##
    ## @item @qcode{'BlockSize'} @tab Memory the expansion may occupy, in
    ## megabytes, @qcode{4e3} by default.
    ##
    ## @item @qcode{'ResponseTransform'} @tab A transformation applied to the
    ## predicted response, named or given as a function handle.
    ##
    ## @item @qcode{'Weights'} @tab One nonnegative weight per observation.
    ##
    ## @item @qcode{'PredictorNames'} @tab One name per predictor.
    ##
    ## @item @qcode{'ResponseName'} @tab A name for the response.
    ##
    ## @item @qcode{'CategoricalPredictors'} @tab Indices of the categorical
    ## predictors.
    ## @end multitable
    ##
    ## The fit is always by limited-memory BFGS, the only solver MATLAB
    ## offers a kernel model, and always under a ridge penalty.
    ##
    ## @seealso{fitrkernel, RegressionLinear}
    ## @end deftypefn
    function this = RegressionKernel (X, Y, varargin)

      if (nargin < 2)
        error ("RegressionKernel: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionKernel: optional arguments must be", ...
                       " given in Name-Value pairs."));
      endif

      ## Defaults
      Learner                = 'svm';
      EpsilonIn              = 'auto';
      EpsilonGiven           = false;
      NumDimsIn              = 'auto';
      KernelScaleIn          = 1;
      LambdaIn               = 'auto';
      BoxConstraint          = 1;
      BoxGiven               = false;
      LambdaGiven            = false;
      Standardize            = false;
      BetaTolerance          = 1e-4;
      GradientTolerance      = 1e-6;
      IterationLimit         = 1000;
      HessianHistorySize     = 15;
      BlockSize              = 4e3;
      Verbose                = 0;
      ResponseTransform      = 'none';
      Weights                = [];
      PredictorNames         = {};
      ResponseName           = 'Y';
      CategoricalPredictors  = [];

      while (numel (varargin) > 0)
        switch (lower (varargin{1}))

          case 'learner'
            Learner = varargin{2};
            if (! (ischar (Learner)
                   && any (strcmpi (Learner, {'svm', 'leastsquares'}))))
              error (strcat ("RegressionKernel: 'Learner' must be either", ...
                             " 'svm' or 'leastsquares'."));
            endif
            Learner = lower (Learner);

          case 'epsilon'
            EpsilonIn = varargin{2};
            EpsilonGiven = true;
            if (! ((ischar (EpsilonIn) && strcmpi (EpsilonIn, 'auto'))
                   || (isnumeric (EpsilonIn) && isscalar (EpsilonIn)
                       && isreal (EpsilonIn) && EpsilonIn >= 0)))
              error (strcat ("RegressionKernel: 'Epsilon' must be 'auto'", ...
                             " or a nonnegative scalar."));
            endif

          case 'numexpansiondimensions'
            NumDimsIn = varargin{2};
            if (! ((ischar (NumDimsIn) && strcmpi (NumDimsIn, 'auto'))
                   || (isnumeric (NumDimsIn) && isscalar (NumDimsIn)
                       && isreal (NumDimsIn) && NumDimsIn > 0
                       && fix (NumDimsIn) == NumDimsIn)))
              error (strcat ("RegressionKernel:", ...
                             " 'NumExpansionDimensions' must be 'auto'", ...
                             " or a positive integer scalar."));
            endif

          case 'kernelscale'
            KernelScaleIn = varargin{2};
            if (! ((ischar (KernelScaleIn)
                    && strcmpi (KernelScaleIn, 'auto'))
                   || (isnumeric (KernelScaleIn)
                       && isscalar (KernelScaleIn)
                       && isreal (KernelScaleIn) && KernelScaleIn > 0)))
              error (strcat ("RegressionKernel: 'KernelScale' must be", ...
                             " 'auto' or a positive scalar."));
            endif

          case 'lambda'
            LambdaIn = varargin{2};
            LambdaGiven = true;
            if (! ((ischar (LambdaIn) && strcmpi (LambdaIn, 'auto'))
                   || (isnumeric (LambdaIn) && isscalar (LambdaIn)
                       && isreal (LambdaIn) && LambdaIn >= 0
                       && isfinite (LambdaIn))))
              error (strcat ("RegressionKernel: 'Lambda' must be 'auto'", ...
                             " or a nonnegative finite scalar."));
            endif

          case 'boxconstraint'
            BoxConstraint = varargin{2};
            BoxGiven = true;
            if (! (isnumeric (BoxConstraint) && isscalar (BoxConstraint)
                   && isreal (BoxConstraint) && BoxConstraint > 0
                   && isfinite (BoxConstraint)))
              error (strcat ("RegressionKernel: 'BoxConstraint' must be", ...
                             " a positive finite scalar."));
            endif

          case 'standardize'
            Standardize = varargin{2};
            if (! (islogical (Standardize) || (isnumeric (Standardize)
                   && isscalar (Standardize)
                   && any (Standardize == [0, 1]))))
              error (strcat ("RegressionKernel: 'Standardize' must be", ...
                             " either true or false."));
            endif
            Standardize = logical (Standardize);

          case 'betatolerance'
            BetaTolerance = varargin{2};
            if (! (isnumeric (BetaTolerance) && isscalar (BetaTolerance)
                   && isreal (BetaTolerance) && BetaTolerance >= 0))
              error (strcat ("RegressionKernel: 'BetaTolerance' must be", ...
                             " a nonnegative scalar."));
            endif

          case 'gradienttolerance'
            GradientTolerance = varargin{2};
            if (! (isnumeric (GradientTolerance)
                   && isscalar (GradientTolerance)
                   && isreal (GradientTolerance) && GradientTolerance >= 0))
              error (strcat ("RegressionKernel: 'GradientTolerance' must", ...
                             " be a nonnegative scalar."));
            endif

          case 'iterationlimit'
            IterationLimit = varargin{2};
            if (! (isnumeric (IterationLimit) && isscalar (IterationLimit)
                   && isreal (IterationLimit) && IterationLimit > 0
                   && fix (IterationLimit) == IterationLimit))
              error (strcat ("RegressionKernel: 'IterationLimit' must be", ...
                             " a positive integer scalar."));
            endif

          case 'hessianhistorysize'
            HessianHistorySize = varargin{2};
            if (! (isnumeric (HessianHistorySize)
                   && isscalar (HessianHistorySize)
                   && isreal (HessianHistorySize) && HessianHistorySize > 0
                   && fix (HessianHistorySize) == HessianHistorySize))
              error (strcat ("RegressionKernel: 'HessianHistorySize'", ...
                             " must be a positive integer scalar."));
            endif

          case 'blocksize'
            BlockSize = varargin{2};
            if (! (isnumeric (BlockSize) && isscalar (BlockSize)
                   && isreal (BlockSize) && BlockSize > 0))
              error (strcat ("RegressionKernel: 'BlockSize' must be a", ...
                             " positive scalar."));
            endif

          case 'verbose'
            Verbose = varargin{2};
            if (! (isnumeric (Verbose) && isscalar (Verbose)
                   && isreal (Verbose) && any (Verbose == [0, 1])))
              error ("RegressionKernel: 'Verbose' must be 0 or 1.");
            endif

          case 'responsetransform'
            ResponseTransform = varargin{2};

          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isreal (Weights)
                   && isvector (Weights) && all (Weights >= 0)))
              error (strcat ("RegressionKernel: 'Weights' must be a", ...
                             " vector of nonnegative values."));
            endif

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! (iscellstr (PredictorNames) && isvector (PredictorNames)))
              error (strcat ("RegressionKernel: 'PredictorNames' must be", ...
                             " a cell array of character vectors."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! (ischar (ResponseName) && isrow (ResponseName)))
              error (strcat ("RegressionKernel: 'ResponseName' must be a", ...
                             " character vector."));
            endif

          case 'categoricalpredictors'
            CategoricalPredictors = varargin{2};
            if (! ((isnumeric (CategoricalPredictors)
                    && isvector (CategoricalPredictors)
                    && all (fix (CategoricalPredictors)
                            == CategoricalPredictors)
                    && all (CategoricalPredictors > 0))
                   || islogical (CategoricalPredictors)
                   || isempty (CategoricalPredictors)))
              error (strcat ("RegressionKernel:", ...
                             " 'CategoricalPredictors' must be a vector", ...
                             " of positive integers or a logical vector."));
            endif

          otherwise
            error (strcat ("RegressionKernel: invalid parameter name in", ...
                           " optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      if (LambdaGiven && BoxGiven)
        error (strcat ("RegressionKernel: 'Lambda' and 'BoxConstraint'", ...
                       " cannot be given together, one being the", ...
                       " reciprocal of the other times the number of", ...
                       " observations."));
      endif
      if (BoxGiven && strcmp (Learner, 'leastsquares'))
        error (strcat ("RegressionKernel: 'BoxConstraint' applies to a", ...
                       " support vector machine only."));
      endif

      ## Validate the data and resolve the weights.  The four linear and
      ## kernel regression models share that opening, so it lives in one
      ## place.
      F = regFrame (X, Y, Weights, 'RegressionKernel');
      X = F.X;
      Y = F.Y;
      W = F.W;
      n = F.n;
      p = F.p;

      ## Epsilon belongs to the insensitive band, so it means nothing to a
      ## least squares fit and is refused there rather than ignored.
      if (strcmp (Learner, 'leastsquares'))
        if (EpsilonGiven)
          error (strcat ("RegressionKernel: 'Epsilon' applies to a", ...
                         " support vector machine only."));
        endif
        Epsilon = [];
      elseif (ischar (EpsilonIn))
        r = iqr (Y);
        if (r == 0)
          Epsilon = 0.1;
        else
          Epsilon = r / 13.49;
        endif
      else
        Epsilon = EpsilonIn;
      endif

      ## Standardize before anything is measured off the predictors
      if (Standardize)
        this.Mu = mean (X, 1);
        this.Sigma = std (X, 0, 1);
        this.Sigma(this.Sigma == 0) = 1;
        X = (X - this.Mu) ./ this.Sigma;
      endif

      ## Resolve the expansion
      if (ischar (NumDimsIn))
        m = 2 .^ ceil (min (log2 (p) + 5, 15));
      else
        m = NumDimsIn;
      endif
      if (ischar (KernelScaleIn))
        sigma = autoKernelScale (X);
      else
        sigma = KernelScaleIn;
      endif

      ## Resolve Lambda and the box constraint from whichever was given
      if (BoxGiven)
        Lambda = 1 / (n * BoxConstraint);
      elseif (LambdaGiven && ! ischar (LambdaIn))
        Lambda = LambdaIn;
        BoxConstraint = 1 / (n * Lambda);
      else
        Lambda = 1 / n;
        BoxConstraint = 1 / (n * Lambda);
      endif

      ## Draw the basis, map the data through it, and fit a linear model
      ## there.
      basis = kernelBasis (p, m, sigma);
      T = kernelExpand (X, basis);

      P = struct ();
      P.Learner = Learner;
      P.LossFunction = 'epsiloninsensitive';
      if (strcmp (Learner, 'leastsquares'))
        P.LossFunction = 'mse';
      endif
      P.Epsilon = Epsilon;
      P.Regularization = 'ridge';
      P.Lambda = Lambda;
      P.Solver = 'lbfgs';
      P.FitBias = true;
      P.PostFitBias = false;
      P.BetaTolerance = BetaTolerance;
      P.GradientTolerance = GradientTolerance;
      P.DeltaGradientTolerance = [];
      P.IterationLimit = IterationLimit;
      P.PassLimit = 1;
      P.BatchSize = 10;
      P.BatchLimit = [];
      P.LearnRate = 1;
      P.OptimizeLearnRate = true;
      P.TruncationPeriod = 10;
      P.NumCheckConvergence = 5;
      P.HessianHistorySize = HessianHistorySize;
      P.InitialBeta = zeros (m, 1);
      P.InitialBias = sum (W .* Y);

      [Beta, Bias, S] = linearSolve (T, Y, W, P);

      ## Fill in the model
      this.Epsilon = Epsilon;
      this.BoxConstraint = BoxConstraint;
      this.PredictorNames = PredictorNames;
      if (isempty (this.PredictorNames))
        this.PredictorNames = ...
                     arrayfun (@(k) sprintf ("x%d", k), 1:p, ...
                               'UniformOutput', false);
      elseif (numel (this.PredictorNames) != p)
        error (strcat ("RegressionKernel: 'PredictorNames' must have one", ...
                       " name per predictor."));
      endif
      this.ExpandedPredictorNames = this.PredictorNames;
      this.CategoricalPredictors = CategoricalPredictors;
      this.ResponseName = ResponseName;
      this.NumExpansionDimensions = m;
      this.FittedLoss = P.LossFunction;
      this.Lambda = Lambda;
      this.KernelScale = sigma;
      this.Learner = Learner;
      this.ResponseTransform = ResponseTransform;
      this.Basis_ = basis;
      this.Beta_ = Beta;
      this.Bias_ = Bias;
      this.NumPredictors_ = p;
      this.NumObservations_ = n;

      this.ModelParameters = struct ('BetaTolerance', BetaTolerance, ...
        'BlockSize', BlockSize, 'BoxConstraint', BoxConstraint, ...
        'Epsilon', EpsilonIn, 'NumExpansionDimensions', NumDimsIn, ...
        'GradientTolerance', GradientTolerance, ...
        'HessianHistorySize', HessianHistorySize, ...
        'IterationLimit', IterationLimit, 'KernelScale', KernelScaleIn, ...
        'Lambda', LambdaIn, 'Learner', Learner, ...
        'LossFunction', P.LossFunction, 'Stream', [], ...
        'VerbosityLevel', Verbose, 'StandardizeData', Standardize, ...
        'Version', 1, 'Method', 'Kernel', 'Type', 'regression');

      this.FitInfo_ = kernelFitInfo (S, P.LossFunction, Lambda, ...
                                     BetaTolerance, GradientTolerance);

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionKernel} {@var{yFit} =} predict (@var{obj}, @var{XC})
    ##
    ## Predict the response of new observations.
    ##
    ## @code{@var{yFit} = predict (@var{obj}, @var{XC})} maps each row of
    ## @var{XC} through the model's own random basis and returns the
    ## predicted response, with @qcode{ResponseTransform} applied.
    ##
    ## @end deftypefn
    function yFit = predict (this, XC)

      if (nargin < 2)
        error ("RegressionKernel.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("RegressionKernel.predict: XC is empty.");
      endif
      if (! (isnumeric (XC) && isreal (XC) && ismatrix (XC)))
        error ("RegressionKernel.predict: invalid values in XC.");
      endif
      if (columns (XC) != this.NumPredictors_)
        error (strcat ("RegressionKernel.predict: XC must have the same", ...
                       " number of predictors as the trained model."));
      endif

      if (! isempty (this.Mu))
        XC = (XC - this.Mu) ./ this.Sigma;
      endif
      yFit = this.RTfun (kernelExpand (XC, this.Basis_) * this.Beta_ ...
                         + this.Bias_);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionKernel} {@var{l} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionKernel} {@var{l} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss on new data.
    ##
    ## @code{@var{l} = loss (@var{obj}, @var{X}, @var{Y})} returns the mean
    ## squared error.
    ##
    ## @code{@var{l} = loss (@dots{}, @var{name}, @var{value})} takes
    ## @qcode{'LossFun'}, either @qcode{'mse'} or
    ## @qcode{'epsiloninsensitive'}, and @qcode{'Weights'}.  The
    ## epsilon-insensitive loss needs a band to be insensitive within, so it
    ## is offered by a support vector machine alone.
    ##
    ## @end deftypefn
    function l = loss (this, X, Y, varargin)

      if (nargin < 3)
        error ("RegressionKernel.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionKernel.loss: optional arguments must", ...
                       " be given in Name-Value pairs."));
      endif

      LossFun = 'mse';
      Weights = [];
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case 'lossfun'
            LossFun = varargin{2};
            if (! (ischar (LossFun) && any (strcmpi (LossFun, ...
                                     {'mse', 'epsiloninsensitive'}))))
              error (strcat ("RegressionKernel.loss: 'LossFun' must be", ...
                             " either 'mse' or 'epsiloninsensitive'."));
            endif
            LossFun = lower (LossFun);
          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isreal (Weights)
                   && isvector (Weights) && all (Weights >= 0)))
              error (strcat ("RegressionKernel.loss: 'Weights' must be a", ...
                             " vector of nonnegative values."));
            endif
          otherwise
            error (strcat ("RegressionKernel.loss: invalid parameter", ...
                           " name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      if (strcmp (LossFun, 'epsiloninsensitive') && isempty (this.Epsilon))
        error (strcat ("RegressionKernel.loss: the", ...
                       " 'epsiloninsensitive' loss applies to a support", ...
                       " vector machine only."));
      endif

      Y = Y(:);
      if (! (isnumeric (Y) && isreal (Y)))
        error ("RegressionKernel.loss: invalid values in Y.");
      endif
      if (rows (X) != numel (Y))
        error (strcat ("RegressionKernel.loss: number of rows in X and Y", ...
                       " must be equal."));
      endif
      if (isempty (Weights))
        w = ones (numel (Y), 1);
      else
        w = Weights(:);
        if (numel (w) != numel (Y))
          error (strcat ("RegressionKernel.loss: 'Weights' must have one", ...
                         " element per observation."));
        endif
      endif
      w = w / sum (w);

      r = Y - predict (this, X);
      if (strcmp (LossFun, 'mse'))
        l = sum (w .* (r .^ 2));
      else
        l = sum (w .* max (0, abs (r) - this.Epsilon));
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionKernel} {@var{obj} =} resume (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionKernel} {@var{obj} =} resume (@dots{}, @var{name}, @var{value})
    ##
    ## Continue fitting a kernel regression model.
    ##
    ## @code{@var{obj} = resume (@var{obj}, @var{X}, @var{Y})} restarts the
    ## optimization from the coefficients the model already carries, through
    ## the basis it already holds.  It takes @qcode{'BetaTolerance'},
    ## @qcode{'GradientTolerance'} and @qcode{'IterationLimit'}, each
    ## defaulting to what the model was fitted with, and @qcode{'Weights'}.
    ##
    ## @var{X} and @var{Y} must be the data the model was fitted to; the
    ## object keeps no copy of them, which is what makes it small.  Neither
    ## does it keep the observation weights, so a model fitted with
    ## @qcode{'Weights'} must be given them again here or it will resume
    ## against uniform ones.  MATLAB behaves the same way: measured on
    ## R2024a, resuming a weighted fit without passing the weights back
    ## reaches the objective of the @emph{unweighted} fit.
    ##
    ## @end deftypefn
    function this = resume (this, X, Y, varargin)

      if (nargin < 3)
        error ("RegressionKernel.resume: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionKernel.resume: optional arguments must", ...
                       " be given in Name-Value pairs."));
      endif

      BetaTolerance = this.ModelParameters.BetaTolerance;
      GradientTolerance = this.ModelParameters.GradientTolerance;
      IterationLimit = this.ModelParameters.IterationLimit;
      Weights = [];
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isreal (Weights)
                   && isvector (Weights) && all (Weights >= 0)))
              error (strcat ("RegressionKernel.resume: 'Weights' must", ...
                             " be a vector of nonnegative values."));
            endif
          case 'betatolerance'
            BetaTolerance = varargin{2};
            if (! (isnumeric (BetaTolerance) && isscalar (BetaTolerance)
                   && isreal (BetaTolerance) && BetaTolerance >= 0))
              error (strcat ("RegressionKernel.resume: 'BetaTolerance'", ...
                             " must be a nonnegative scalar."));
            endif
          case 'gradienttolerance'
            GradientTolerance = varargin{2};
            if (! (isnumeric (GradientTolerance)
                   && isscalar (GradientTolerance)
                   && isreal (GradientTolerance) && GradientTolerance >= 0))
              error (strcat ("RegressionKernel.resume:", ...
                             " 'GradientTolerance' must be a", ...
                             " nonnegative scalar."));
            endif
          case 'iterationlimit'
            IterationLimit = varargin{2};
            if (! (isnumeric (IterationLimit) && isscalar (IterationLimit)
                   && isreal (IterationLimit) && IterationLimit > 0
                   && fix (IterationLimit) == IterationLimit))
              error (strcat ("RegressionKernel.resume:", ...
                             " 'IterationLimit' must be a positive", ...
                             " integer scalar."));
            endif
          otherwise
            error (strcat ("RegressionKernel.resume: invalid parameter", ...
                           " name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
        error ("RegressionKernel.resume: invalid values in X.");
      endif
      if (columns (X) != this.NumPredictors_)
        error (strcat ("RegressionKernel.resume: X must have the same", ...
                       " number of predictors as the trained model."));
      endif
      Y = Y(:);
      if (! (isnumeric (Y) && isreal (Y)))
        error ("RegressionKernel.resume: invalid values in Y.");
      endif
      if (rows (X) != numel (Y))
        error (strcat ("RegressionKernel.resume: number of rows in X and", ...
                       " Y must be equal."));
      endif

      if (! isempty (this.Mu))
        X = (X - this.Mu) ./ this.Sigma;
      endif
      T = kernelExpand (X, this.Basis_);
      if (isempty (Weights))
        W = ones (numel (Y), 1);
      else
        W = Weights(:);
        if (numel (W) != numel (Y))
          error (strcat ("RegressionKernel.resume: 'Weights' must have", ...
                         " one element per observation."));
        endif
      endif
      W = W / sum (W);

      P = struct ();
      P.Learner = this.Learner;
      P.LossFunction = this.FittedLoss;
      P.Epsilon = this.Epsilon;
      P.Regularization = 'ridge';
      P.Lambda = this.Lambda;
      P.Solver = 'lbfgs';
      P.FitBias = true;
      P.PostFitBias = false;
      P.BetaTolerance = BetaTolerance;
      P.GradientTolerance = GradientTolerance;
      P.DeltaGradientTolerance = [];
      P.IterationLimit = IterationLimit;
      P.PassLimit = 1;
      P.BatchSize = 10;
      P.BatchLimit = [];
      P.LearnRate = 1;
      P.OptimizeLearnRate = true;
      P.TruncationPeriod = 10;
      P.NumCheckConvergence = 5;
      P.HessianHistorySize = this.ModelParameters.HessianHistorySize;
      P.InitialBeta = this.Beta_;
      P.InitialBias = this.Bias_;

      [Beta, Bias, S] = linearSolve (T, Y, W, P);

      this.Beta_ = Beta;
      this.Bias_ = Bias;
      this.ModelParameters.BetaTolerance = BetaTolerance;
      this.ModelParameters.GradientTolerance = GradientTolerance;
      this.ModelParameters.IterationLimit = IterationLimit;
      this.FitInfo_ = kernelFitInfo (S, this.FittedLoss, this.Lambda, ...
                                     BetaTolerance, GradientTolerance);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionKernel} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a kernel regression model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves the model
    ## @var{obj} into @var{filename} in a form @code{loadmodel} can read
    ## back, the random basis included.
    ##
    ## @end deftypefn
    function savemodel (obj, fname)

      classdef_name = 'RegressionKernel';
      Epsilon = obj.Epsilon;
      BoxConstraint = obj.BoxConstraint;
      ResponseTransform = obj.ResponseTransform;
      PredictorNames = obj.PredictorNames;
      CategoricalPredictors = obj.CategoricalPredictors;
      ResponseName = obj.ResponseName;
      ExpandedPredictorNames = obj.ExpandedPredictorNames;
      NumExpansionDimensions = obj.NumExpansionDimensions;
      FittedLoss = obj.FittedLoss;
      Lambda = obj.Lambda;
      ModelParameters = obj.ModelParameters;
      Regularization = obj.Regularization;
      KernelScale = obj.KernelScale;
      Learner = obj.Learner;
      Mu = obj.Mu;
      Sigma = obj.Sigma;
      Basis_ = obj.Basis_;
      Beta_ = obj.Beta_;
      Bias_ = obj.Bias_;
      NumPredictors_ = obj.NumPredictors_;
      NumObservations_ = obj.NumObservations_;

      save ('-binary', fname, 'classdef_name', 'Epsilon', ...
            'BoxConstraint', 'ResponseTransform', 'PredictorNames', ...
            'CategoricalPredictors', 'ResponseName', ...
            'ExpandedPredictorNames', 'NumExpansionDimensions', ...
            'FittedLoss', 'Lambda', 'ModelParameters', 'Regularization', ...
            'KernelScale', 'Learner', 'Mu', 'Sigma', 'Basis_', 'Beta_', ...
            'Bias_', 'NumPredictors_', 'NumObservations_');

    endfunction

  endmethods

  methods (Access = public, Hidden)

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        printf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      printf ("\n  RegressionKernel\n\n");
      printf ("%+26s: '%s'\n", 'ResponseName', this.ResponseName);
      printf ("%+26s: '%s'\n", 'Learner', this.Learner);
      printf ("%+26s: %d\n", 'NumExpansionDimensions', ...
              this.NumExpansionDimensions);
      printf ("%+26s: %g\n", 'KernelScale', this.KernelScale);
      printf ("%+26s: %g\n", 'Lambda', this.Lambda);
      printf ("%+26s: %g\n", 'BoxConstraint', this.BoxConstraint);
      if (! isempty (this.Epsilon))
        printf ("%+26s: %g\n", 'Epsilon', this.Epsilon);
      endif
      printf ("\n");
    endfunction

    ## What the fit reported, which fitrkernel returns as its second output.
    function S = fitInfo_ (this)
      S = this.FitInfo_;
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      mdl = RegressionKernel ([0; 1], [0; 1]);
      fields = fieldnames (data);
      for k = 1:numel (fields)
        mdl.(fields{k}) = data.(fields{k});
      endfor

    endfunction

  endmethods

endclassdef

%!demo
%! ## Fit fuel consumption through a randomized Gaussian kernel, which
%! ## bends where a linear model cannot.
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! Mdl = RegressionKernel (X(ok,:), MPG(ok))
%! yFit = predict (Mdl, X(find (ok, 3),:))

%!demo
%! ## Standardizing matters more here than for a linear fit, since the
%! ## kernel measures one distance across predictors of every scale.
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! plain = RegressionKernel (X(ok,:), MPG(ok));
%! scaled = RegressionKernel (X(ok,:), MPG(ok), 'Standardize', true);
%! plainLoss = loss (plain, X(ok,:), MPG(ok))
%! scaledLoss = loss (scaled, X(ok,:), MPG(ok))

%!shared X, Y
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! X = X(ok,:);
%! Y = MPG(ok);

%!test
%! ## The model reports the surface MATLAB reports
%! Mdl = RegressionKernel (X, Y);
%! assert_equal (class (Mdl), 'RegressionKernel');
%! assert_equal (Mdl.Learner, 'svm');
%! assert_equal (Mdl.FittedLoss, 'epsiloninsensitive');
%! assert_equal (Mdl.Regularization, 'ridge (L2)');
%! assert_equal (Mdl.ResponseTransform, 'none');
%! assert_equal (Mdl.KernelScale, 1);
%! assert_equal (Mdl.BoxConstraint, 1);
%! assert_equal (Mdl.NumExpansionDimensions, 128);
%! assert_equal (Mdl.Mu, []);
%! assert_equal (Mdl.Sigma, []);
%! assert_equal (Mdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'});

%!test
%! ## The properties are the ones MATLAB lists, in its order
%! Mdl = RegressionKernel (X, Y);
%! assert_equal (properties (Mdl), {'Epsilon'; 'BoxConstraint'; ...
%!                                  'ResponseTransform'; ...
%!                                  'PredictorNames'; ...
%!                                  'CategoricalPredictors'; ...
%!                                  'ResponseName'; ...
%!                                  'ExpandedPredictorNames'; ...
%!                                  'NumExpansionDimensions'; ...
%!                                  'FittedLoss'; 'Lambda'; ...
%!                                  'ModelParameters'; 'Regularization'; ...
%!                                  'KernelScale'; 'Learner'; 'Mu'; 'Sigma'});

%!test
%! ## Epsilon defaults to the interquartile range over 13.49, as it does for
%! ## the linear model
%! Mdl = RegressionKernel (X, Y);
%! assert_equal (Mdl.Epsilon, 0.926612305411416, 1e-12);
%! assert_equal (Mdl.Epsilon, iqr (Y) / 13.49, 1e-15);

%!test
%! ## Least squares has no insensitive band at all
%! Mdl = RegressionKernel (X, Y, 'Learner', 'leastsquares');
%! assert_equal (Mdl.Epsilon, []);
%! assert_equal (Mdl.FittedLoss, 'mse');

%!test
%! ## Lambda and the box constraint are reciprocal through the number of
%! ## observations, and either one may be the one that is given
%! Mb = RegressionKernel (X, Y, 'BoxConstraint', 4);
%! assert_equal (Mb.Lambda, 1 / (93 * 4), 1e-15);
%! Ml = RegressionKernel (X, Y, 'Lambda', 0.02);
%! assert_equal (Ml.BoxConstraint, 1 / (93 * 0.02), 1e-12);

%!test
%! ## Lambda defaults to the reciprocal of the observations that were used
%! Mdl = RegressionKernel (X, Y);
%! assert_equal (Mdl.Lambda, 1 / 93, 1e-15);

%!test
%! ## The default expansion is MATLAB's, two to the power of five more than
%! ## the base two logarithm of the predictors
%! assert_equal (RegressionKernel (X(:,1:2), Y).NumExpansionDimensions, 64);
%! assert_equal (RegressionKernel (X, Y, ...
%!               'NumExpansionDimensions', 50).NumExpansionDimensions, 50);

%!test
%! ## An odd expansion cannot be paired throughout, and the dimension left
%! ## over still comes back as a dimension
%! Mdl = RegressionKernel (X, Y, 'NumExpansionDimensions', 51);
%! assert_equal (Mdl.NumExpansionDimensions, 51);
%! assert_equal (numel (predict (Mdl, X)), 93);

%!test
%! ## Standardizing records the means and deviations of the predictors
%! Mdl = RegressionKernel (X, Y, 'Standardize', true);
%! assert_equal (Mdl.Mu, mean (X), 1e-12);
%! assert_equal (Mdl.Sigma, std (X), 1e-12);

%!test
%! ## A kernel fit follows the response it was given
%! Mdl = RegressionKernel (X, Y, 'Learner', 'leastsquares', ...
%!                         'Standardize', true, 'Lambda', 1e-4);
%! assert_equal (loss (Mdl, X, Y) < var (Y), true);
%! assert_equal (corr (predict (Mdl, X), Y) > 0.5, true);

%!test
%! ## Predicting through the model's own basis gives the same answer every
%! ## time it is asked
%! Mdl = RegressionKernel (X, Y);
%! assert_equal (predict (Mdl, X(1:10,:)), predict (Mdl, X(1:10,:)));

%!test
%! ## resume continues from the coefficients the model already holds, so it
%! ## cannot leave the objective higher than it found it
%! Mdl = RegressionKernel (X, Y, 'IterationLimit', 3);
%! before = Mdl.FitInfo_.ObjectiveValue;
%! Mdl = resume (Mdl, X, Y, 'IterationLimit', 500);
%! assert_equal (class (Mdl), 'RegressionKernel');
%! assert_equal (Mdl.FitInfo_.ObjectiveValue <= before, true);
%! assert_equal (Mdl.ModelParameters.IterationLimit, 500);

%!test
%! ## The fit information is MATLAB's kernel structure
%! Mdl = RegressionKernel (X, Y);
%! F = Mdl.FitInfo_;
%! assert_equal (fieldnames (F), {'Solver'; 'LossFunction'; 'Lambda'; ...
%!                                'BetaTolerance'; 'GradientTolerance'; ...
%!                                'ObjectiveValue'; 'GradientMagnitude'; ...
%!                                'RelativeChangeInBeta'; 'FitTime'; ...
%!                                'History'});
%! assert_equal (F.Solver, 'LBFGS-fast');
%! assert_equal (F.LossFunction, 'epsiloninsensitive');

%!test
%! ## A response transform reaches predict
%! Mdl = RegressionKernel (X, Y, 'Learner', 'leastsquares');
%! plain = predict (Mdl, X(1:5,:));
%! Mdl.ResponseTransform = 'exp';
%! assert_equal (predict (Mdl, X(1:5,:)), exp (plain), 1e-12);

%!test
%! ## A row with a missing predictor or a missing response is dropped
%! Xn = X;
%! Xn(3,2) = NaN;
%! Mdl = RegressionKernel (Xn, Y);
%! assert_equal (Mdl.Lambda, 1 / 92, 1e-15);

%!test
%! ## A saved model reads back as the same model, the random basis included
%! Mdl = RegressionKernel (X, Y);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mnew = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (Mnew), 'RegressionKernel');
%! assert_equal (predict (Mnew, X(1:5,:)), predict (Mdl, X(1:5,:)));

%!test
%! ## resume keeps no observation weights, because the model keeps no data:
%! ## passing them back restores the weighted fit, and omitting them
%! ## continues against uniform ones
%! Mdl = RegressionKernel (X, Y, 'IterationLimit', 5, 'Weights', (1:93)');
%! kept = resume (Mdl, X, Y, 'IterationLimit', 400, 'Weights', (1:93)');
%! lost = resume (Mdl, X, Y, 'IterationLimit', 400);
%! assert_equal (kept.FitInfo_.ObjectiveValue ...
%!               != lost.FitInfo_.ObjectiveValue, true);

## Test input validation
%!error<RegressionKernel: too few input arguments.> RegressionKernel (ones (5, 2))
%!error<RegressionKernel: optional arguments must be given in Name-Value pairs.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Learner')
%!error<RegressionKernel: 'Learner' must be either 'svm' or 'leastsquares'.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Learner', 'logistic')
%!error<RegressionKernel: 'Epsilon' must be 'auto' or a nonnegative scalar.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Epsilon', -1)
%!error<RegressionKernel: 'Epsilon' applies to a support vector machine only.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Learner', 'leastsquares', ...
%!                     'Epsilon', 1)
%!error<RegressionKernel: 'NumExpansionDimensions' must be 'auto' or a positive integer scalar.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'NumExpansionDimensions', 2.5)
%!error<RegressionKernel: 'KernelScale' must be 'auto' or a positive scalar.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'KernelScale', 0)
%!error<RegressionKernel: 'Lambda' must be 'auto' or a nonnegative finite scalar.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Lambda', Inf)
%!error<RegressionKernel: 'BoxConstraint' must be a positive finite scalar.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'BoxConstraint', -2)
%!error<RegressionKernel: 'Lambda' and 'BoxConstraint' cannot be given together, one being the reciprocal of the other times the number of observations.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Lambda', 0.1, ...
%!                     'BoxConstraint', 2)
%!error<RegressionKernel: 'BoxConstraint' applies to a support vector machine only.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Learner', 'leastsquares', ...
%!                     'BoxConstraint', 2)
%!error<RegressionKernel: 'Standardize' must be either true or false.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Standardize', 'yes')
%!error<RegressionKernel: 'Verbose' must be 0 or 1.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Verbose', 2)
%!error<RegressionKernel: invalid parameter name in optional pair arguments.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Nonsense', 1)
%!error<RegressionKernel: invalid values in X.> RegressionKernel ({1, 2; 3, 4}, [1; 2])
%!error<RegressionKernel: X is empty.> RegressionKernel ([], [])
%!error<RegressionKernel: invalid values in Y.> RegressionKernel (ones (10, 2), {1, 2})
%!error<RegressionKernel: number of rows in X and Y must be equal.> ...
%! RegressionKernel (ones (10, 2), ones (3, 1))
%!error<RegressionKernel: 'Weights' must have one element per observation.> ...
%! RegressionKernel (ones (10, 2), ones (10, 1), 'Weights', ones (3, 1))
%!error<RegressionKernel.predict: XC is empty.> ...
%! predict (RegressionKernel (ones (10, 2), ones (10, 1)), [])
%!error<RegressionKernel.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (RegressionKernel (ones (10, 2), ones (10, 1)), ones (3, 5))
%!error<RegressionKernel.loss: 'LossFun' must be either 'mse' or 'epsiloninsensitive'.> ...
%! loss (RegressionKernel (ones (10, 2), ones (10, 1)), ones (10, 2), ...
%!                     ones (10, 1), 'LossFun', 'hinge')
%!error<RegressionKernel.loss: the 'epsiloninsensitive' loss applies to a support vector machine only.> ...
%! loss (RegressionKernel (ones (10, 2), ones (10, 1), 'Learner', ...
%!                     'leastsquares'), ones (10, 2), ones (10, 1), ...
%!                     'LossFun', 'epsiloninsensitive')
%!error<RegressionKernel.resume: too few input arguments.> ...
%! resume (RegressionKernel (ones (10, 2), ones (10, 1)), ones (10, 2))
%!error<RegressionKernel.resume: X must have the same number of predictors as the trained model.> ...
%! resume (RegressionKernel (ones (10, 2), ones (10, 1)), ones (10, 5), ...
%!                     ones (10, 1))
%!error<RegressionKernel.resume: invalid parameter name in optional pair arguments.> ...
%! resume (RegressionKernel (ones (10, 2), ones (10, 1)), ones (10, 2), ...
%!                     ones (10, 1), 'Nonsense', 1)
%!error<RegressionKernel.resume: 'Weights' must be a vector of nonnegative values.> ...
%! resume (RegressionKernel (ones (10, 2), ones (10, 1)), ones (10, 2), ...
%!                     ones (10, 1), 'Weights', -1)
