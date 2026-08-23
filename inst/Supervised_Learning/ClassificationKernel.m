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
## @deftp {statistics} ClassificationKernel
##
## Gaussian kernel binary classifier for large data.
##
## A @qcode{ClassificationKernel} object maps the predictors into a
## randomized feature space whose inner product approximates a Gaussian
## kernel, and then fits a linear model there.  A kernel classifier is
## therefore as nonlinear as a support vector machine with a Gaussian kernel,
## while costing what a linear fit costs: nothing of size @math{NxN} is ever
## formed.
##
## The expansion is the random Fourier basis of Rahimi and Recht, drawn once
## when the model is fitted and kept with it, so @code{predict} maps new data
## through the same basis.  MATLAB approximates the same kernel by the
## Fastfood construction, which reaches the same distribution more cheaply;
## the two are interchangeable in distribution but not draw by draw, and the
## draws come from different generators in any case, so the scores of a model
## fitted here and one fitted in MATLAB differ even from the same seed.
## What does not differ is what they estimate.
##
## Like @qcode{ClassificationLinear} the object holds no copy of the training
## data.  It does hold the basis and the coefficients, so it is bounded by
## the number of expansion dimensions rather than by the number of
## observations.
##
## Create a @qcode{ClassificationKernel} object with @code{fitckernel}.
##
## @seealso{fitckernel, ClassificationLinear, ClassificationSVM}
## @end deftp

classdef ClassificationKernel

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} BoxConstraint
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

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} ClassNames
    ##
    ## Names of the two classes
    ##
    ## A column of the same type as the response supplied to the constructor.
    ## The second of the two is the positive class, the one a positive score
    ## belongs to.  This property is read-only.
    ##
    ## @end deftp
    ClassNames             = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} Prior
    ##
    ## Prior probability of each class
    ##
    ## A numeric row vector with one element per class, in the order of
    ## @qcode{ClassNames} and summing to one.  It defaults to the class
    ## frequencies of the training data.  This property is read-only.
    ##
    ## @end deftp
    Prior                  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} Cost
    ##
    ## Cost of misclassifying an observation
    ##
    ## A square numeric matrix with one row and one column per class, whose
    ## @math{(i,j)} element is the cost of classifying an observation of
    ## class @math{i} into class @math{j}.  It defaults to one everywhere
    ## except the diagonal, which is zero.  This property is read-only, as it
    ## is in MATLAB; a cost matrix is given to the constructor instead.
    ##
    ## The costs are folded into the prior before the observations are
    ## weighted, so a class that is costlier to misclassify weighs more in
    ## the fit.  They are read again by the @qcode{'mincost'} and
    ## @qcode{'classifcost'} losses.
    ##
    ## @end deftp
    Cost                   = [];

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} ScoreTransform
    ##
    ## Transformation applied to the predicted scores
    ##
    ## A character vector naming a transformation, or the text of the
    ## function handle that was supplied.  Assigning to it accepts either.
    ## It defaults to @qcode{'logit'} for a logistic learner and to
    ## @qcode{'none'} for a support vector machine.
    ##
    ## @end deftp
    ScoreTransform         = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} PredictorNames
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
    ## @deftp {ClassificationKernel} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A row vector of column indices, empty when every predictor is
    ## numeric.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} ResponseName
    ##
    ## Name of the response
    ##
    ## A character vector, defaulting to @qcode{'Y'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    ResponseName           = 'Y';

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} ExpandedPredictorNames
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
    ## @deftp {ClassificationKernel} {property} NumExpansionDimensions
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
    ## @deftp {ClassificationKernel} {property} FittedLoss
    ##
    ## Loss function the fit minimized
    ##
    ## @qcode{'hinge'} for a support vector machine and @qcode{'logit'} for a
    ## logistic regression.  This property is read-only.
    ##
    ## @end deftp
    FittedLoss             = 'hinge';

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} Lambda
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
    ## @deftp {ClassificationKernel} {property} ModelParameters
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
    ## @deftp {ClassificationKernel} {property} Regularization
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
    ## @deftp {ClassificationKernel} {property} KernelScale
    ##
    ## Scale of the Gaussian kernel
    ##
    ## A positive scalar dividing every predictor before the expansion, so a
    ## larger scale makes the kernel wider and the classifier smoother.  This
    ## property is read-only.
    ##
    ## @end deftp
    KernelScale            = 1;

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} Learner
    ##
    ## Linear model fitted in the expanded space
    ##
    ## Either @qcode{'svm'} or @qcode{'logistic'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Learner                = 'svm';

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} Mu
    ##
    ## Predictor means used to standardize
    ##
    ## A row vector with one element per predictor, or empty when the model
    ## was fitted without standardizing.  This property is read-only.
    ##
    ## @end deftp
    Mu                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationKernel} {property} Sigma
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

    ## The callable behind ScoreTransform.
    STfun                  = @(s) s;

    ## The random basis and the coefficients fitted in the space it spans.
    ## MATLAB keeps all three out of sight, so they are hidden here too, but
    ## they are the model: without the basis the coefficients index nothing.
    Basis_                 = [];
    Beta_                  = [];
    Bias_                  = [];

    ## Number of predictors and of observations the fit saw.  Neither is a
    ## property of MATLAB's class, but predict needs the first to validate
    ## its input and resume needs the second to keep Lambda and BoxConstraint
    ## reciprocal.
    NumPredictors_         = [];
    NumObservations_       = [];

    ## What the fit reported, so that fitckernel can hand it back as its
    ## second output.
    FitInfo_               = [];

  endproperties

  methods (Access = public)

    ## Custom setter, so that assigning a name or a handle updates both the
    ## text the property reports and the callable predict uses.
    function this = set.ScoreTransform (this, val)
      [this.STfun, this.ScoreTransform] = ...
                     parseScoreTransform (val, 'ClassificationKernel');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKernel} {@var{obj} =} ClassificationKernel (@var{X}, @var{Y})
    ## @deftypefnx {ClassificationKernel} {@var{obj} =} ClassificationKernel (@dots{}, @var{name}, @var{value})
    ##
    ## Fit a Gaussian kernel binary classifier.
    ##
    ## @code{@var{obj} = ClassificationKernel (@var{X}, @var{Y})} fits a
    ## support vector machine in a randomized Gaussian kernel space to the
    ## @math{NxP} predictor matrix @var{X} and the @math{Nx1} response
    ## @var{Y}, which must name exactly two classes.
    ##
    ## @code{@var{obj} = ClassificationKernel (@dots{}, @var{name},
    ## @var{value})} takes the following @qcode{Name-Value} pairs.
    ##
    ## @multitable @columnfractions 0.28 0.72
    ## @headitem Name @tab Value
    ##
    ## @item @qcode{'Learner'} @tab @qcode{'svm'}, the default, or
    ## @qcode{'logistic'}.
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
    ## @item @qcode{'ClassNames'} @tab The classes to keep, given in the type
    ## of @var{Y}.
    ##
    ## @item @qcode{'Cost'} @tab A square misclassification cost matrix.
    ##
    ## @item @qcode{'Prior'} @tab @qcode{'empirical'}, the default,
    ## @qcode{'uniform'}, a vector of probabilities, or a structure with
    ## @qcode{ClassNames} and @qcode{ClassProbs} fields.
    ##
    ## @item @qcode{'ScoreTransform'} @tab A transformation applied to the
    ## scores, named or given as a function handle.
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
    ## @seealso{fitckernel, ClassificationLinear}
    ## @end deftypefn
    function this = ClassificationKernel (X, Y, varargin)

      if (nargin < 2)
        error ("ClassificationKernel: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationKernel: optional arguments must", ...
                       " be given in Name-Value pairs."));
      endif

      ## Defaults
      Learner                = 'svm';
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
      ClassNames             = [];
      CostIn                 = [];
      Prior                  = [];
      ScoreTransform         = [];
      Weights                = [];
      PredictorNames         = {};
      ResponseName           = 'Y';
      CategoricalPredictors  = [];

      while (numel (varargin) > 0)
        switch (lower (varargin{1}))

          case 'learner'
            Learner = varargin{2};
            if (! (ischar (Learner)
                   && any (strcmpi (Learner, {'svm', 'logistic'}))))
              error (strcat ("ClassificationKernel: 'Learner' must be", ...
                             " either 'svm' or 'logistic'."));
            endif
            Learner = lower (Learner);

          case 'numexpansiondimensions'
            NumDimsIn = varargin{2};
            if (! ((ischar (NumDimsIn) && strcmpi (NumDimsIn, 'auto'))
                   || (isnumeric (NumDimsIn) && isscalar (NumDimsIn)
                       && isreal (NumDimsIn) && NumDimsIn > 0
                       && fix (NumDimsIn) == NumDimsIn)))
              error (strcat ("ClassificationKernel:", ...
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
              error (strcat ("ClassificationKernel: 'KernelScale' must", ...
                             " be 'auto' or a positive scalar."));
            endif

          case 'lambda'
            LambdaIn = varargin{2};
            LambdaGiven = true;
            if (! ((ischar (LambdaIn) && strcmpi (LambdaIn, 'auto'))
                   || (isnumeric (LambdaIn) && isscalar (LambdaIn)
                       && isreal (LambdaIn) && LambdaIn >= 0
                       && isfinite (LambdaIn))))
              error (strcat ("ClassificationKernel: 'Lambda' must be", ...
                             " 'auto' or a nonnegative finite scalar."));
            endif

          case 'boxconstraint'
            BoxConstraint = varargin{2};
            BoxGiven = true;
            if (! (isnumeric (BoxConstraint) && isscalar (BoxConstraint)
                   && isreal (BoxConstraint) && BoxConstraint > 0
                   && isfinite (BoxConstraint)))
              error (strcat ("ClassificationKernel: 'BoxConstraint' must", ...
                             " be a positive finite scalar."));
            endif

          case 'standardize'
            Standardize = varargin{2};
            if (! (islogical (Standardize) || (isnumeric (Standardize)
                   && isscalar (Standardize)
                   && any (Standardize == [0, 1]))))
              error (strcat ("ClassificationKernel: 'Standardize' must", ...
                             " be either true or false."));
            endif
            Standardize = logical (Standardize);

          case 'betatolerance'
            BetaTolerance = varargin{2};
            if (! (isnumeric (BetaTolerance) && isscalar (BetaTolerance)
                   && isreal (BetaTolerance) && BetaTolerance >= 0))
              error (strcat ("ClassificationKernel: 'BetaTolerance' must", ...
                             " be a nonnegative scalar."));
            endif

          case 'gradienttolerance'
            GradientTolerance = varargin{2};
            if (! (isnumeric (GradientTolerance)
                   && isscalar (GradientTolerance)
                   && isreal (GradientTolerance) && GradientTolerance >= 0))
              error (strcat ("ClassificationKernel: 'GradientTolerance'", ...
                             " must be a nonnegative scalar."));
            endif

          case 'iterationlimit'
            IterationLimit = varargin{2};
            if (! (isnumeric (IterationLimit) && isscalar (IterationLimit)
                   && isreal (IterationLimit) && IterationLimit > 0
                   && fix (IterationLimit) == IterationLimit))
              error (strcat ("ClassificationKernel: 'IterationLimit'", ...
                             " must be a positive integer scalar."));
            endif

          case 'hessianhistorysize'
            HessianHistorySize = varargin{2};
            if (! (isnumeric (HessianHistorySize)
                   && isscalar (HessianHistorySize)
                   && isreal (HessianHistorySize) && HessianHistorySize > 0
                   && fix (HessianHistorySize) == HessianHistorySize))
              error (strcat ("ClassificationKernel:", ...
                             " 'HessianHistorySize' must be a positive", ...
                             " integer scalar."));
            endif

          case 'blocksize'
            BlockSize = varargin{2};
            if (! (isnumeric (BlockSize) && isscalar (BlockSize)
                   && isreal (BlockSize) && BlockSize > 0))
              error (strcat ("ClassificationKernel: 'BlockSize' must be", ...
                             " a positive scalar."));
            endif

          case 'verbose'
            Verbose = varargin{2};
            if (! (isnumeric (Verbose) && isscalar (Verbose)
                   && isreal (Verbose) && any (Verbose == [0, 1])))
              error (strcat ("ClassificationKernel: 'Verbose' must be 0", ...
                             " or 1."));
            endif

          case 'classnames'
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames)
                   || islogical (ClassNames) || ischar (ClassNames)))
              error (strcat ("ClassificationKernel: 'ClassNames' must be", ...
                             " a cell array of character vectors, a", ...
                             " logical vector, a numeric vector, or a", ...
                             " character array."));
            endif

          case 'cost'
            CostIn = varargin{2};
            if (! (isnumeric (CostIn) && isreal (CostIn)
                   && ismatrix (CostIn) && ndims (CostIn) == 2
                   && rows (CostIn) == columns (CostIn)))
              error (strcat ("ClassificationKernel: 'Cost' must be a", ...
                             " square numeric matrix."));
            endif

          case 'prior'
            Prior = varargin{2};
            if (! ((ischar (Prior) && any (strcmpi (Prior, {'empirical', ...
                                                            'uniform'})))
                   || (isnumeric (Prior) && isreal (Prior)
                       && isvector (Prior) && all (Prior >= 0))
                   || (isstruct (Prior) && isscalar (Prior))))
              error (strcat ("ClassificationKernel: 'Prior' must be", ...
                             " 'empirical', 'uniform', a vector of", ...
                             " nonnegative values, or a structure."));
            endif

          case 'scoretransform'
            ScoreTransform = varargin{2};

          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isreal (Weights)
                   && isvector (Weights) && all (Weights >= 0)))
              error (strcat ("ClassificationKernel: 'Weights' must be a", ...
                             " vector of nonnegative values."));
            endif

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! (iscellstr (PredictorNames) && isvector (PredictorNames)))
              error (strcat ("ClassificationKernel: 'PredictorNames'", ...
                             " must be a cell array of character", ...
                             " vectors."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! (ischar (ResponseName) && isrow (ResponseName)))
              error (strcat ("ClassificationKernel: 'ResponseName' must", ...
                             " be a character vector."));
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
              error (strcat ("ClassificationKernel:", ...
                             " 'CategoricalPredictors' must be a vector", ...
                             " of positive integers or a logical", ...
                             " vector."));
            endif

          otherwise
            error (strcat ("ClassificationKernel: invalid parameter name", ...
                           " in optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Lambda and the box constraint are reciprocal, so naming both
      ## overdetermines the fit rather than describing it twice.
      if (LambdaGiven && BoxGiven)
        error (strcat ("ClassificationKernel: 'Lambda' and", ...
                       " 'BoxConstraint' cannot be given together, one", ...
                       " being the reciprocal of the other times the", ...
                       " number of observations."));
      endif
      if (BoxGiven && strcmp (Learner, 'logistic'))
        error (strcat ("ClassificationKernel: 'BoxConstraint' applies to", ...
                       " a support vector machine only."));
      endif

      ## Validate the data, resolve the classes, the prior, the cost and
      ## the observation weights.  The four linear and kernel classifiers
      ## share that opening, so it lives in one place.
      F = classFrame (X, Y, ClassNames, Prior, CostIn, Weights, ...
                      'ClassificationKernel');
      X = F.X;
      y = F.y;
      W = F.W;
      n = F.n;
      p = F.p;
      this.ClassNames = F.ClassNames;
      this.Prior = F.Prior;
      this.Cost = F.Cost;

      ## Standardize before anything is measured off the predictors, so the
      ## kernel scale and the expansion both see the same data predict will.
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
      ## there.  The basis is what makes the fit nonlinear; everything after
      ## it is the linear machinery.
      basis = kernelBasis (p, m, sigma);
      T = kernelExpand (X, basis);

      P = struct ();
      P.Learner = Learner;
      P.LossFunction = 'hinge';
      if (strcmp (Learner, 'logistic'))
        P.LossFunction = 'logit';
      endif
      P.Epsilon = [];
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
      P.InitialBias = 0;

      [Beta, Bias, S] = linearSolve (T, y, W, P);

      ## Fill in the model
      this.BoxConstraint = BoxConstraint;
      this.PredictorNames = PredictorNames;
      if (isempty (this.PredictorNames))
        this.PredictorNames = ...
                     arrayfun (@(k) sprintf ("x%d", k), 1:p, ...
                               'UniformOutput', false);
      elseif (numel (this.PredictorNames) != p)
        error (strcat ("ClassificationKernel: 'PredictorNames' must have", ...
                       " one name per predictor."));
      endif
      this.ExpandedPredictorNames = this.PredictorNames;
      this.CategoricalPredictors = CategoricalPredictors;
      this.ResponseName = ResponseName;
      this.NumExpansionDimensions = m;
      this.FittedLoss = P.LossFunction;
      this.Lambda = Lambda;
      this.KernelScale = sigma;
      this.Learner = Learner;
      this.Basis_ = basis;
      this.Beta_ = Beta;
      this.Bias_ = Bias;
      this.NumPredictors_ = p;
      this.NumObservations_ = n;

      if (isempty (ScoreTransform))
        if (strcmp (Learner, 'logistic'))
          ScoreTransform = 'logit';
        else
          ScoreTransform = 'none';
        endif
      endif
      this.ScoreTransform = ScoreTransform;

      this.ModelParameters = struct ('BetaTolerance', BetaTolerance, ...
        'BlockSize', BlockSize, 'BoxConstraint', BoxConstraint, ...
        'Epsilon', 'auto', 'NumExpansionDimensions', NumDimsIn, ...
        'GradientTolerance', GradientTolerance, ...
        'HessianHistorySize', HessianHistorySize, ...
        'IterationLimit', IterationLimit, 'KernelScale', KernelScaleIn, ...
        'Lambda', LambdaIn, 'Learner', Learner, ...
        'LossFunction', P.LossFunction, 'Stream', [], ...
        'VerbosityLevel', Verbose, 'StandardizeData', Standardize, ...
        'Version', 1, 'Method', 'Kernel', 'Type', 'classification');

      this.FitInfo_ = kernelFitInfo (S, P.LossFunction, Lambda, ...
                                     BetaTolerance, GradientTolerance);

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKernel} {@var{labels} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationKernel} {[@var{labels}, @var{scores}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new observations.
    ##
    ## @code{@var{labels} = predict (@var{obj}, @var{XC})} maps each row of
    ## @var{XC} through the model's own random basis and returns the class of
    ## largest score.
    ##
    ## @code{[@var{labels}, @var{scores}] = predict (@var{obj}, @var{XC})}
    ## also returns the @math{Nx2} scores, whose columns follow
    ## @qcode{ClassNames}, after @qcode{ScoreTransform} has been applied.
    ##
    ## @end deftypefn
    function [labels, scores] = predict (this, XC)

      if (nargin < 2)
        error ("ClassificationKernel.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("ClassificationKernel.predict: XC is empty.");
      endif
      if (! (isnumeric (XC) && isreal (XC) && ismatrix (XC)))
        error ("ClassificationKernel.predict: invalid values in XC.");
      endif
      if (columns (XC) != this.NumPredictors_)
        error (strcat ("ClassificationKernel.predict: XC must have the", ...
                       " same number of predictors as the trained model."));
      endif

      f = rawScore (this, XC);
      labels = labelsFromIndex (this.ClassNames, 1 + (f > 0));

      if (nargout > 1)
        scores = this.STfun ([-f, f]);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKernel} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin of each observation.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns the
    ## score of the true class less the score of the other one.  A positive
    ## margin is a correct classification.
    ##
    ## @end deftypefn
    function m = margin (this, X, Y)

      if (nargin < 3)
        error ("ClassificationKernel.margin: too few input arguments.");
      endif
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("ClassificationKernel.margin: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("ClassificationKernel.margin: number of rows in X", ...
                       " and Y must be equal."));
      endif

      [~, s] = predict (this, X);
      n = rows (s);
      strue = s(sub2ind (size (s), (1:n)', gY));
      sother = s(sub2ind (size (s), (1:n)', 3 - gY));
      m = strue - sother;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKernel} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationKernel} {@var{e} =} edge (@dots{}, @qcode{'Weights'}, @var{W})
    ##
    ## Weighted mean of the classification margins.
    ##
    ## The weights are normalized within each class to that class's prior
    ## before they are applied.
    ##
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      if (nargin < 3)
        error ("ClassificationKernel.edge: too few input arguments.");
      endif
      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       'ClassificationKernel', 'edge');
      e = sum (W .* margin (this, X, Y));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKernel} {@var{l} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationKernel} {@var{l} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss on new data.
    ##
    ## @code{@var{l} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## misclassification rate.
    ##
    ## @code{@var{l} = loss (@dots{}, @var{name}, @var{value})} takes
    ## @qcode{'LossFun'}, one of @qcode{'binodeviance'},
    ## @qcode{'classifcost'}, @qcode{'classiferror'}, @qcode{'exponential'},
    ## @qcode{'hinge'}, @qcode{'logit'}, @qcode{'mincost'} and
    ## @qcode{'quadratic'}, and @qcode{'Weights'}.
    ##
    ## @end deftypefn
    function l = loss (this, X, Y, varargin)

      if (nargin < 3)
        error ("ClassificationKernel.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationKernel.loss: optional arguments", ...
                       " must be given in Name-Value pairs."));
      endif

      LossFun = 'classiferror';
      Weights = [];
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))
          case 'lossfun'
            LossFun = varargin{2};
            valid = {'binodeviance', 'classifcost', 'classiferror', ...
                     'exponential', 'hinge', 'logit', 'mincost', ...
                     'quadratic'};
            if (! (ischar (LossFun) && any (strcmpi (LossFun, valid))))
              error (strcat ("ClassificationKernel.loss: 'LossFun' must", ...
                             " be 'binodeviance', 'classifcost',", ...
                             " 'classiferror', 'exponential', 'hinge',", ...
                             " 'logit', 'mincost', or 'quadratic'."));
            endif
            LossFun = lower (LossFun);
          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isreal (Weights)
                   && isvector (Weights) && all (Weights >= 0)))
              error (strcat ("ClassificationKernel.loss: 'Weights' must", ...
                             " be a vector of nonnegative values."));
            endif
          otherwise
            error (strcat ("ClassificationKernel.loss: invalid parameter", ...
                           " name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("ClassificationKernel.loss: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("ClassificationKernel.loss: number of rows in X", ...
                       " and Y must be equal."));
      endif
      if (isempty (Weights))
        w = ones (numel (gY), 1);
      else
        w = Weights(:);
        if (numel (w) != numel (gY))
          error (strcat ("ClassificationKernel.loss: 'Weights' must have", ...
                         " one element per observation."));
        endif
      endif
      w = w / sum (w);

      [~, s] = predict (this, X);
      l = classificationLoss (LossFun, s, gY, w, this.Cost);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKernel} {@var{obj} =} resume (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationKernel} {@var{obj} =} resume (@dots{}, @var{name}, @var{value})
    ##
    ## Continue fitting a kernel classifier.
    ##
    ## @code{@var{obj} = resume (@var{obj}, @var{X}, @var{Y})} restarts the
    ## optimization from the coefficients the model already carries, through
    ## the basis it already holds, and returns the model it reaches.  It
    ## takes @qcode{'BetaTolerance'}, @qcode{'GradientTolerance'} and
    ## @qcode{'IterationLimit'}, each defaulting to what the model was
    ## fitted with, and @qcode{'Weights'}.
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
        error ("ClassificationKernel.resume: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationKernel.resume: optional arguments", ...
                       " must be given in Name-Value pairs."));
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
              error (strcat ("ClassificationKernel.resume: 'Weights'", ...
                             " must be a vector of nonnegative values."));
            endif
          case 'betatolerance'
            BetaTolerance = varargin{2};
            if (! (isnumeric (BetaTolerance) && isscalar (BetaTolerance)
                   && isreal (BetaTolerance) && BetaTolerance >= 0))
              error (strcat ("ClassificationKernel.resume:", ...
                             " 'BetaTolerance' must be a nonnegative", ...
                             " scalar."));
            endif
          case 'gradienttolerance'
            GradientTolerance = varargin{2};
            if (! (isnumeric (GradientTolerance)
                   && isscalar (GradientTolerance)
                   && isreal (GradientTolerance) && GradientTolerance >= 0))
              error (strcat ("ClassificationKernel.resume:", ...
                             " 'GradientTolerance' must be a", ...
                             " nonnegative scalar."));
            endif
          case 'iterationlimit'
            IterationLimit = varargin{2};
            if (! (isnumeric (IterationLimit) && isscalar (IterationLimit)
                   && isreal (IterationLimit) && IterationLimit > 0
                   && fix (IterationLimit) == IterationLimit))
              error (strcat ("ClassificationKernel.resume:", ...
                             " 'IterationLimit' must be a positive", ...
                             " integer scalar."));
            endif
          otherwise
            error (strcat ("ClassificationKernel.resume: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      [T, y, W] = resumeData (this, X, Y, Weights, 'resume');

      P = resumeOptions (this, BetaTolerance, GradientTolerance, ...
                         IterationLimit);
      P.InitialBeta = this.Beta_;
      P.InitialBias = this.Bias_;
      [Beta, Bias, S] = linearSolve (T, y, W, P);

      this.Beta_ = Beta;
      this.Bias_ = Bias;
      this.ModelParameters.BetaTolerance = BetaTolerance;
      this.ModelParameters.GradientTolerance = GradientTolerance;
      this.ModelParameters.IterationLimit = IterationLimit;
      this.FitInfo_ = kernelFitInfo (S, this.FittedLoss, this.Lambda, ...
                                     BetaTolerance, GradientTolerance);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationKernel} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a kernel classifier to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves the model
    ## @var{obj} into @var{filename} in a form @code{loadmodel} can read
    ## back, the random basis included.
    ##
    ## @end deftypefn
    function savemodel (obj, fname)

      classdef_name = 'ClassificationKernel';
      BoxConstraint = obj.BoxConstraint;
      ClassNames = obj.ClassNames;
      Prior = obj.Prior;
      Cost = obj.Cost;
      ScoreTransform = obj.ScoreTransform;
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

      save ('-binary', fname, 'classdef_name', 'BoxConstraint', ...
            'ClassNames', 'Prior', 'Cost', 'ScoreTransform', ...
            'PredictorNames', 'CategoricalPredictors', 'ResponseName', ...
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
      printf ("\n  ClassificationKernel\n\n");
      printf ("%+26s: '%s'\n", 'ResponseName', this.ResponseName);
      printf ("%+26s: %s\n", 'ClassNames', ...
              classNameListing (this.ClassNames));
      printf ("%+26s: '%s'\n", 'Learner', this.Learner);
      printf ("%+26s: %d\n", 'NumExpansionDimensions', ...
              this.NumExpansionDimensions);
      printf ("%+26s: %g\n", 'KernelScale', this.KernelScale);
      printf ("%+26s: %g\n", 'Lambda', this.Lambda);
      printf ("%+26s: %g\n", 'BoxConstraint', this.BoxConstraint);
      printf ("\n");
    endfunction

    ## What the fit reported, which fitckernel returns as its second output.
    function S = fitInfo_ (this)
      S = this.FitInfo_;
    endfunction

  endmethods

  methods (Access = private)

    ## The raw model value of each row of XC, standardized and expanded
    ## through the model's own basis first.
    function f = rawScore (this, XC)

      if (! isempty (this.Mu))
        XC = (XC - this.Mu) ./ this.Sigma;
      endif
      f = kernelExpand (XC, this.Basis_) * this.Beta_ + this.Bias_;

    endfunction

    ## The expanded predictors, the signed response and the weights of a
    ## resume call, validated the way the constructor validated them.
    function [T, y, W] = resumeData (this, X, Y, Weights, caller)

      if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
        error ("ClassificationKernel.%s: invalid values in X.", caller);
      endif
      if (columns (X) != this.NumPredictors_)
        error (strcat ("ClassificationKernel.%s: X must have the same", ...
                       " number of predictors as the trained model."), ...
               caller);
      endif
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("ClassificationKernel.%s: %s", caller, errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("ClassificationKernel.%s: number of rows in X and", ...
                       " Y must be equal."), caller);
      endif

      y = -ones (numel (gY), 1);
      y(gY == 2) = 1;

      if (isempty (Weights))
        Weights = ones (numel (gY), 1);
      else
        Weights = Weights(:);
        if (numel (Weights) != numel (gY))
          error (strcat ("ClassificationKernel.%s: 'Weights' must have", ...
                         " one element per observation."), caller);
        endif
      endif

      adjPrior = this.Prior .* sum (this.Cost, 2)';
      adjPrior = adjPrior / sum (adjPrior);
      W = zeros (numel (gY), 1);
      for k = 1:numel (this.Prior)
        idx = (gY == k);
        tot = sum (Weights(idx));
        if (tot > 0)
          W(idx) = Weights(idx) / tot * adjPrior(k);
        endif
      endfor
      W = W / sum (W);

      if (! isempty (this.Mu))
        X = (X - this.Mu) ./ this.Sigma;
      endif
      T = kernelExpand (X, this.Basis_);

    endfunction

    ## The option structure linearSolve takes, filled from the model.
    function P = resumeOptions (this, BetaTol, GradTol, IterLimit)

      P = struct ();
      P.Learner = this.Learner;
      P.LossFunction = this.FittedLoss;
      P.Epsilon = [];
      P.Regularization = 'ridge';
      P.Lambda = this.Lambda;
      P.Solver = 'lbfgs';
      P.FitBias = true;
      P.PostFitBias = false;
      P.BetaTolerance = BetaTol;
      P.GradientTolerance = GradTol;
      P.DeltaGradientTolerance = [];
      P.IterationLimit = IterLimit;
      P.PassLimit = 1;
      P.BatchSize = 10;
      P.BatchLimit = [];
      P.LearnRate = 1;
      P.OptimizeLearnRate = true;
      P.TruncationPeriod = 10;
      P.NumCheckConvergence = 5;
      P.HessianHistorySize = this.ModelParameters.HessianHistorySize;

    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      mdl = ClassificationKernel (zeros (2, 1), [0; 1]);
      fields = fieldnames (data);
      for k = 1:numel (fields)
        mdl.(fields{k}) = data.(fields{k});
      endfor

    endfunction

  endmethods

endclassdef

%!demo
%! ## Separate the two overlapping iris species through a randomized
%! ## Gaussian kernel, and read the model the fit produced.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationKernel (X, Y)
%! predict (Mdl, X([1, 51],:))

%!demo
%! ## The box constraint and the regularization strength are two names for
%! ## the same quantity: setting either fixes the other.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationKernel (X, Y, 'BoxConstraint', 4);
%! Mdl.BoxConstraint
%! Mdl.Lambda

%!test
%! ## The model reports the surface MATLAB reports
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationKernel (X, Y);
%! assert_equal (class (Mdl), 'ClassificationKernel');
%! assert_equal (Mdl.ClassNames, {'versicolor'; 'virginica'});
%! assert_equal (Mdl.Learner, 'svm');
%! assert_equal (Mdl.FittedLoss, 'hinge');
%! assert_equal (Mdl.Regularization, 'ridge (L2)');
%! assert_equal (Mdl.ScoreTransform, 'none');
%! assert_equal (Mdl.KernelScale, 1);
%! assert_equal (Mdl.BoxConstraint, 1);
%! assert_equal (Mdl.Lambda, 0.01);
%! assert_equal (Mdl.NumExpansionDimensions, 128);
%! assert_equal (Mdl.Mu, []);
%! assert_equal (Mdl.Sigma, []);
%! assert_equal (Mdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'});

%!test
%! ## The properties are the ones MATLAB lists, in its order
%! load fisheriris
%! Mdl = ClassificationKernel (meas(51:end,:), species(51:end));
%! assert_equal (properties (Mdl), {'BoxConstraint'; 'ClassNames'; ...
%!                                  'Prior'; 'Cost'; 'ScoreTransform'; ...
%!                                  'PredictorNames'; ...
%!                                  'CategoricalPredictors'; ...
%!                                  'ResponseName'; ...
%!                                  'ExpandedPredictorNames'; ...
%!                                  'NumExpansionDimensions'; ...
%!                                  'FittedLoss'; 'Lambda'; ...
%!                                  'ModelParameters'; 'Regularization'; ...
%!                                  'KernelScale'; 'Learner'; 'Mu'; 'Sigma'});

%!test
%! ## The default expansion is MATLAB's, two to the power of five more than
%! ## the base two logarithm of the predictors, capped at fifteen
%! X = randn (40, 2);
%! Y = [ones(20, 1); 2 * ones(20, 1)];
%! assert_equal (ClassificationKernel (X, Y).NumExpansionDimensions, 64);
%! assert_equal (ClassificationKernel (randn (40, 32), Y) ...
%!               .NumExpansionDimensions, 1024);

%!test
%! ## Lambda and the box constraint are reciprocal through the number of
%! ## observations, and either one may be the one that is given
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mb = ClassificationKernel (X, Y, 'BoxConstraint', 3);
%! assert_equal (Mb.Lambda, 1 / 300, 1e-15);
%! assert_equal (Mb.BoxConstraint, 3);
%! Ml = ClassificationKernel (X, Y, 'Lambda', 0.05);
%! assert_equal (Ml.Lambda, 0.05);
%! assert_equal (Ml.BoxConstraint, 0.2, 1e-15);

%!test
%! ## Standardizing records the means and deviations MATLAB records
%! load fisheriris
%! Mdl = ClassificationKernel (meas(51:end,:), species(51:end), ...
%!                             'Standardize', true);
%! assert_equal (Mdl.Mu, [6.262, 2.872, 4.906, 1.676], 1e-12);
%! assert_equal (Mdl.Sigma, [0.662834440074967, 0.332751006494695, ...
%!                           0.82557846264289, 0.424768504986284], 1e-12);

%!test
%! ## A logistic learner fits the deviance and reports posteriors
%! load fisheriris
%! X = meas(51:end,:);
%! Mdl = ClassificationKernel (X, species(51:end), 'Learner', 'logistic');
%! assert_equal (Mdl.FittedLoss, 'logit');
%! assert_equal (Mdl.ScoreTransform, 'logit');
%! [~, score] = predict (Mdl, X(1:5,:));
%! assert_equal (sum (score, 2), ones (5, 1), 1e-12);

%!test
%! ## A support vector machine leaves the scores untransformed, so they are
%! ## a value and its negative
%! load fisheriris
%! X = meas(51:end,:);
%! Mdl = ClassificationKernel (X, species(51:end));
%! [~, score] = predict (Mdl, X(1:5,:));
%! assert_equal (score(:,1), -score(:,2), 1e-12);

%!test
%! ## The fit separates the two species it was given, whatever basis it drew
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationKernel (X, Y);
%! assert_equal (loss (Mdl, X, Y) < 0.15, true);
%! assert_equal (edge (Mdl, X, Y) > 0, true);

%!test
%! ## margin is the true class score less the other, and the labels follow
%! ## the sign of the raw score
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationKernel (X, Y);
%! m = margin (Mdl, X, Y);
%! assert_equal (size (m), [100, 1]);
%! assert_equal (mean (m > 0), 1 - loss (Mdl, X, Y), 1e-12);

%!test
%! ## Predicting through the model's own basis is what makes it a model at
%! ## all: the same rows give the same scores every time it is asked
%! load fisheriris
%! X = meas(51:end,:);
%! Mdl = ClassificationKernel (X, species(51:end));
%! [~, s1] = predict (Mdl, X(1:10,:));
%! [~, s2] = predict (Mdl, X(1:10,:));
%! assert_equal (s1, s2);

%!test
%! ## A wider kernel gives a smoother rule, so it fits the training data
%! ## less closely than a narrow one does
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mnarrow = ClassificationKernel (X, Y, 'KernelScale', 0.5, ...
%!                                 'Lambda', 1e-4);
%! Mwide = ClassificationKernel (X, Y, 'KernelScale', 20, 'Lambda', 1e-4);
%! assert_equal (loss (Mnarrow, X, Y) <= loss (Mwide, X, Y), true);

%!test
%! ## resume continues from the coefficients the model already holds, so it
%! ## cannot leave the objective higher than it found it
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationKernel (X, Y, 'IterationLimit', 3);
%! before = Mdl.FitInfo_.ObjectiveValue;
%! Mdl = resume (Mdl, X, Y, 'IterationLimit', 500);
%! assert_equal (class (Mdl), 'ClassificationKernel');
%! assert_equal (Mdl.FitInfo_.ObjectiveValue <= before, true);
%! assert_equal (Mdl.ModelParameters.IterationLimit, 500);

%!test
%! ## The fit information is MATLAB's kernel structure, not its linear one
%! load fisheriris
%! Mdl = ClassificationKernel (meas(51:end,:), species(51:end));
%! F = Mdl.FitInfo_;
%! assert_equal (fieldnames (F), {'Solver'; 'LossFunction'; 'Lambda'; ...
%!                                'BetaTolerance'; 'GradientTolerance'; ...
%!                                'ObjectiveValue'; 'GradientMagnitude'; ...
%!                                'RelativeChangeInBeta'; 'FitTime'; ...
%!                                'History'});
%! assert_equal (F.Solver, 'LBFGS-fast');
%! assert_equal (F.LossFunction, 'hinge');
%! assert_equal (F.Lambda, 0.01);

%!test
%! ## A cost matrix reaches the fit through the prior, as it does for the
%! ## linear classifier
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationKernel (X, Y, 'Cost', [0, 4; 1, 0]);
%! assert_equal (Mdl.Cost, [0, 4; 1, 0]);
%! assert_equal (Mdl.Prior, [0.5, 0.5]);

%!test
%! ## A saved model reads back as the same model, the random basis included
%! load fisheriris
%! X = meas(51:end,:);
%! Mdl = ClassificationKernel (X, species(51:end));
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mnew = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (Mnew), 'ClassificationKernel');
%! [~, s1] = predict (Mdl, X(1:5,:));
%! [~, s2] = predict (Mnew, X(1:5,:));
%! assert_equal (s1, s2);

%!test
%! ## The labels come back in the type the response was given in, a
%! ## character matrix included
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mchar = ClassificationKernel (X, char (Y));
%! assert_equal (size (Mchar.ClassNames), [2, 10]);
%! assert_equal (size (predict (Mchar, X(1:3,:))), [3, 10]);
%! Mnum = ClassificationKernel (X, double (strcmp (Y, 'virginica')));
%! assert_equal (Mnum.ClassNames, [0; 1]);
%! assert_equal (class (predict (Mnum, X(1:3,:))), 'double');

%!test
%! ## A character matrix response names one class per row, and every method
%! ## answers the same through it as through the equivalent cell array
%! load fisheriris
%! X = meas(51:end,:);
%! Yc = species(51:end);
%! Ym = char (Yc);
%! rand ('seed', 11); randn ('seed', 11);
%! Mc = ClassificationKernel (X, Yc);
%! rand ('seed', 11); randn ('seed', 11);
%! Mm = ClassificationKernel (X, Ym);
%! assert_equal (cellstr (Mm.ClassNames), Mc.ClassNames);
%! assert_equal (cellstr (predict (Mm, X)), predict (Mc, X));
%! assert_equal (margin (Mm, X, Ym), margin (Mc, X, Yc));
%! assert_equal (edge (Mm, X, Ym), edge (Mc, X, Yc));
%! assert_equal (loss (Mm, X, Ym), loss (Mc, X, Yc));

%!test
%! ## 'ClassNames' selects a subset when it is given as a character matrix
%! load fisheriris
%! Mdl = ClassificationKernel (meas, char (species), 'ClassNames', ...
%!                             char ({'versicolor', 'virginica'}));
%! assert_equal (size (Mdl.ClassNames), [2, 10]);
%! assert_equal (cellstr (Mdl.ClassNames), {'versicolor'; 'virginica'});

%!test
%! ## resume keeps no observation weights, because the model keeps no data:
%! ## a weighted fit resumed without them continues against uniform ones,
%! ## and passing them back restores the weighted fit.  MATLAB behaves the
%! ## same way, measured on R2024a.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! w = (1:100)';
%! rand ('seed', 4); randn ('seed', 4);
%! Mdl = ClassificationKernel (X, Y, 'IterationLimit', 5, 'Weights', w);
%! kept = resume (Mdl, X, Y, 'IterationLimit', 400, 'Weights', w);
%! lost = resume (Mdl, X, Y, 'IterationLimit', 400);
%! assert_equal (kept.FitInfo_.ObjectiveValue ...
%!               != lost.FitInfo_.ObjectiveValue, true);
%! rand ('seed', 4); randn ('seed', 4);
%! direct = ClassificationKernel (X, Y, 'IterationLimit', 405, 'Weights', w);
%! assert_equal (kept.FitInfo_.ObjectiveValue, ...
%!               direct.FitInfo_.ObjectiveValue, 1e-3);

%!test
%! ## A character matrix response survives a round trip through savemodel
%! ## and loadmodel, the random basis with it
%! load fisheriris
%! X = meas(51:end,:);
%! Ym = char (species(51:end));
%! Mdl = ClassificationKernel (X, Ym);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mnew = loadmodel (fname);
%! delete (fname);
%! assert_equal (Mnew.ClassNames, Mdl.ClassNames);
%! assert_equal (predict (Mnew, X(1:5,:)), predict (Mdl, X(1:5,:)));

## Test input validation
%!error<ClassificationKernel: too few input arguments.> ...
%! ClassificationKernel (ones (5, 2))
%!error<ClassificationKernel: optional arguments must be given in Name-Value pairs.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Learner')
%!error<ClassificationKernel: 'Learner' must be either 'svm' or 'logistic'.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Learner', ...
%!                     'tree')
%!error<ClassificationKernel: 'NumExpansionDimensions' must be 'auto' or a positive integer scalar.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'NumExpansionDimensions', 0)
%!error<ClassificationKernel: 'KernelScale' must be 'auto' or a positive scalar.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'KernelScale', -1)
%!error<ClassificationKernel: 'Lambda' must be 'auto' or a nonnegative finite scalar.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Lambda', -1)
%!error<ClassificationKernel: 'BoxConstraint' must be a positive finite scalar.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'BoxConstraint', 0)
%!error<ClassificationKernel: 'Lambda' and 'BoxConstraint' cannot be given together, one being the reciprocal of the other times the number of observations.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Lambda', ...
%!                     0.1, 'BoxConstraint', 2)
%!error<ClassificationKernel: 'BoxConstraint' applies to a support vector machine only.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Learner', ...
%!                     'logistic', 'BoxConstraint', 2)
%!error<ClassificationKernel: 'Standardize' must be either true or false.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'Standardize', 'yes')
%!error<ClassificationKernel: 'IterationLimit' must be a positive integer scalar.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'IterationLimit', -5)
%!error<ClassificationKernel: 'Verbose' must be 0 or 1.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Verbose', 2)
%!error<ClassificationKernel: invalid parameter name in optional pair arguments.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Nonsense', 1)
%!error<ClassificationKernel: invalid values in X.> ...
%! ClassificationKernel ({1, 2; 3, 4}, [1; 2])
%!error<ClassificationKernel: X is empty.> ClassificationKernel ([], [])
%!error<ClassificationKernel: number of rows in X and Y must be equal.> ...
%! ClassificationKernel (ones (10, 2), [1; 2])
%!error<ClassificationKernel: Y must name exactly two classes, this being a binary model; use fitcecoc for more than two.> ...
%! ClassificationKernel (ones (9, 2), [1; 1; 1; 2; 2; 2; 3; 3; 3])
%!error<ClassificationKernel: the number of rows and columns in 'Cost' must correspond to the classes in Y.> ...
%! ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Cost', ...
%!                     ones (3))
%!error<ClassificationKernel.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (3, 5))
%!error<ClassificationKernel.predict: XC is empty.> ...
%! predict (ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)]), [])
%!error<ClassificationKernel.margin: too few input arguments.> ...
%! margin (ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (3, 2))
%!error<ClassificationKernel.loss: 'LossFun' must be 'binodeviance', 'classifcost', 'classiferror', 'exponential', 'hinge', 'logit', 'mincost', or 'quadratic'.> ...
%! loss (ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (10, 2), [ones(5,1); 2*ones(5,1)], 'LossFun', 'mse')
%!error<ClassificationKernel.resume: too few input arguments.> ...
%! resume (ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (10, 2))
%!error<ClassificationKernel.resume: 'IterationLimit' must be a positive integer scalar.> ...
%! resume (ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'IterationLimit', 0)
%!error<ClassificationKernel.resume: invalid parameter name in optional pair arguments.> ...
%! resume (ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Nonsense', 1)
%!error<ClassificationKernel.resume: 'Weights' must be a vector of nonnegative values.> ...
%! resume (ClassificationKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Weights', -1)
