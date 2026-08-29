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
## @deftp {statistics} RegressionLinear
##
## Linear regression model for high dimensional data.
##
## A @qcode{RegressionLinear} object fits a linear model,
## @code{@var{X} * Beta + Bias}, to a continuous response by minimizing a
## regularized average loss.  The loss is the epsilon-insensitive loss for a
## support vector machine and the squared error for a least squares fit, and
## the penalty is either a ridge or a lasso one.
##
## Unlike the other regression models of this package the object holds no
## copy of the training data: the coefficients, the intercept and the fitting
## options are the whole model.  That is what makes it suited to data with
## more predictors than a kernel matrix could carry, and it is why the class
## has no @code{compact} method and no resubstitution methods.
##
## A vector of regularization strengths fits one model per value in a single
## object.  @qcode{Beta} is then a @math{PxL} matrix and @qcode{Bias} a
## @math{1xL} row, every method returns one column per strength, and
## @code{selectModels} narrows the object down to the strengths worth
## keeping.
##
## Create a @qcode{RegressionLinear} object with @code{fitrlinear}.
##
## @seealso{fitrlinear, RegressionKernel, RegressionSVM}
## @end deftp

classdef RegressionLinear

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} Epsilon
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

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} ResponseTransform
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
    ## @deftp {RegressionLinear} {property} PredictorNames
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
    ## @deftp {RegressionLinear} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A row vector of column indices, empty when every predictor is
    ## numeric.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors  = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} ResponseName
    ##
    ## Name of the response
    ##
    ## A character vector, defaulting to @qcode{'Y'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    ResponseName           = 'Y';

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the fit saw them
    ##
    ## A cell array of character vectors.  It equals @qcode{PredictorNames}
    ## unless categorical predictors were expanded into indicator variables.
    ## This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} Learner
    ##
    ## Linear regression model that was fitted
    ##
    ## Either @qcode{'svm'} or @qcode{'leastsquares'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Learner                = 'svm';

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} Beta
    ##
    ## Fitted linear coefficients
    ##
    ## A @math{Px1} column, or a @math{PxL} matrix with one column per
    ## regularization strength when @qcode{Lambda} holds more than one.  This
    ## property is read-only.
    ##
    ## @end deftp
    Beta                   = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} Bias
    ##
    ## Fitted intercept
    ##
    ## A scalar, or a @math{1xL} row with one element per regularization
    ## strength.  It is zero throughout when the model was fitted with
    ## @qcode{'FitBias'} set to false.  This property is read-only.
    ##
    ## @end deftp
    Bias                   = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} FittedLoss
    ##
    ## Loss function the fit minimized
    ##
    ## @qcode{'epsiloninsensitive'} for a support vector machine and
    ## @qcode{'mse'} for a least squares fit.  This is the loss of the
    ## objective, which is not the loss @code{loss} reports unless it is
    ## asked for.  This property is read-only.
    ##
    ## @end deftp
    FittedLoss             = 'epsiloninsensitive';

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} Lambda
    ##
    ## Regularization strength
    ##
    ## A nonnegative scalar, or a @math{1xL} row of them in ascending order.
    ## It defaults to the reciprocal of the number of observations used to
    ## train the model.  This property is read-only.
    ##
    ## @end deftp
    Lambda                 = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} ModelParameters
    ##
    ## Fitting options, as they were given
    ##
    ## A structure holding every parameter of the fit, including the ones
    ## that a different solver would have used and the @qcode{'auto'} values
    ## before they were resolved.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters        = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionLinear} {property} Regularization
    ##
    ## Penalty on the coefficients
    ##
    ## @qcode{'ridge (L2)'} or @qcode{'lasso (L1)'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Regularization         = 'ridge (L2)';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## The callable behind ResponseTransform.  The public property is the
    ## text MATLAB stores; this is what predict actually applies.
    RTfun                  = @(y) y;

    ## Number of predictors the fit saw, which the class does not report
    ## because MATLAB does not, but predict needs to validate its input.
    NumPredictors_         = [];

    ## What the fit reported, so that fitrlinear can hand it back as its
    ## second output.  MATLAB returns it from the fitting function rather
    ## than storing it on the model, and a constructor cannot have a second
    ## output argument.
    FitInfo_               = [];

  endproperties

  methods (Access = public)

    ## Custom setter, so that assigning a name or a handle updates both the
    ## text the property reports and the callable predict uses.
    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = ...
                     parseResponseTransform (val, 'RegressionLinear');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionLinear} {@var{obj} =} RegressionLinear (@var{X}, @var{Y})
    ## @deftypefnx {RegressionLinear} {@var{obj} =} RegressionLinear (@dots{}, @var{name}, @var{value})
    ##
    ## Fit a linear regression model.
    ##
    ## @code{@var{obj} = RegressionLinear (@var{X}, @var{Y})} fits a linear
    ## support vector machine to the @math{NxP} predictor matrix @var{X} and
    ## the @math{Nx1} continuous response @var{Y}.
    ##
    ## @code{@var{obj} = RegressionLinear (@dots{}, @var{name},
    ## @var{value})} takes the following @qcode{Name-Value} pairs.
    ##
    ## @multitable @columnfractions 0.28 0.72
    ## @headitem Name @tab Value
    ##
    ## @item @qcode{'Learner'} @tab @qcode{'svm'}, the default, or
    ## @qcode{'leastsquares'}.  The first minimizes the epsilon-insensitive
    ## loss and the second the squared error.
    ##
    ## @item @qcode{'Epsilon'} @tab Half the width of the insensitive band,
    ## a nonnegative scalar or @qcode{'auto'}, which is the interquartile
    ## range of @var{Y} over 13.49.  It applies to a support vector machine
    ## alone.
    ##
    ## @item @qcode{'Regularization'} @tab @qcode{'ridge'} or
    ## @qcode{'lasso'}.  It defaults to @qcode{'lasso'} when the solver is
    ## @qcode{'sparsa'} and to @qcode{'ridge'} otherwise.
    ##
    ## @item @qcode{'Lambda'} @tab @qcode{'auto'}, the default, which is the
    ## reciprocal of the number of observations, or a nonnegative scalar, or
    ## a vector of them.  A vector fits one model per value.
    ##
    ## @item @qcode{'Solver'} @tab One of @qcode{'sgd'}, @qcode{'asgd'},
    ## @qcode{'dual'}, @qcode{'bfgs'}, @qcode{'lbfgs'} and @qcode{'sparsa'},
    ## or a cell array of them applied in turn, each warm starting the next.
    ##
    ## @item @qcode{'Beta'} @tab Initial coefficients, a @math{Px1} column or
    ## a @math{PxL} matrix.  It defaults to zeros.
    ##
    ## @item @qcode{'Bias'} @tab Initial intercept, a scalar or a @math{1xL}
    ## row.  It defaults to the weighted mean of @var{Y} for a least squares
    ## fit and to its weighted median for a support vector machine.
    ##
    ## @item @qcode{'FitBias'} @tab Whether to fit an intercept at all, true
    ## by default.
    ##
    ## @item @qcode{'PostFitBias'} @tab Whether to refit the intercept once
    ## the coefficients are settled, false by default.
    ##
    ## @item @qcode{'ObservationsIn'} @tab @qcode{'rows'}, the default, or
    ## @qcode{'columns'}, which transposes @var{X} before fitting.
    ##
    ## @item @qcode{'BetaTolerance'} @tab Relative tolerance on the
    ## coefficients, @qcode{1e-4} by default.
    ##
    ## @item @qcode{'GradientTolerance'} @tab Absolute tolerance on the
    ## gradient's infinity norm, @qcode{1e-6} by default.
    ##
    ## @item @qcode{'DeltaGradientTolerance'} @tab Tolerance on the
    ## complementarity gap of the @qcode{'dual'} solver, @qcode{0.1} by
    ## default.
    ##
    ## @item @qcode{'IterationLimit'} @tab Largest number of iterations,
    ## @qcode{1000} by default.
    ##
    ## @item @qcode{'PassLimit'} @tab Largest number of passes over the data
    ## for the stochastic solvers, @qcode{1} by default, and @qcode{10} for
    ## @qcode{'dual'}.
    ##
    ## @item @qcode{'BatchSize'} @tab Mini-batch size of the stochastic
    ## solvers, @qcode{10} by default.
    ##
    ## @item @qcode{'BatchLimit'} @tab Largest number of mini-batches.
    ##
    ## @item @qcode{'LearnRate'} @tab Step size of the stochastic solvers.
    ##
    ## @item @qcode{'OptimizeLearnRate'} @tab Whether to halve the step size
    ## when the objective rises, true by default.
    ##
    ## @item @qcode{'TruncationPeriod'} @tab Number of mini-batches between
    ## soft thresholdings under a lasso penalty, @qcode{10} by default.
    ##
    ## @item @qcode{'NumCheckConvergence'} @tab Number of passes between
    ## convergence checks of the @qcode{'dual'} solver, @qcode{2} by
    ## default.  MathWorks documents @qcode{5}; R2024a and R2026a both
    ## report @qcode{2}.
    ##
    ## @item @qcode{'HessianHistorySize'} @tab Number of curvature pairs the
    ## quasi-Newton solvers keep, @qcode{15} by default.
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
    ## The default solver is @qcode{'sparsa'} under a lasso penalty.  Under a
    ## ridge penalty it is @qcode{'bfgs'} when there are no more than 100
    ## predictors, and beyond that @qcode{'dual'} for a support vector
    ## machine and @qcode{'sgd'} for a least squares fit.
    ##
    ## @seealso{fitrlinear, RegressionKernel}
    ## @end deftypefn
    function this = RegressionLinear (X, Y, varargin)

      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("RegressionLinear: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionLinear: optional arguments must be", ...
                       " given in Name-Value pairs."));
      endif

      ## Defaults, before the optional arguments are parsed
      Learner                = 'svm';
      EpsilonIn              = 'auto';
      Regularization         = [];
      Lambda                 = 'auto';
      Solver                 = [];
      BetaIn                 = [];
      BiasIn                 = [];
      FitBias                = true;
      PostFitBias            = false;
      ObservationsIn         = 'rows';
      BetaTolerance          = 1e-4;
      GradientTolerance      = 1e-6;
      DeltaGradientTolerance = 0.1;
      IterationLimit         = 1000;
      PassLimit              = [];
      BatchSize              = 10;
      BatchLimit             = [];
      LearnRate              = [];
      OptimizeLearnRate      = true;
      TruncationPeriod       = 10;
      NumCheckConvergence    = 2;
      HessianHistorySize     = 15;
      Verbose                = 0;
      ResponseTransform      = 'none';
      Weights                = [];
      PredictorNames         = {};
      ResponseName           = 'Y';
      CategoricalPredictors  = [];
      EpsilonGiven           = false;

      ## Parse optional parameters
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))

          case 'learner'
            Learner = varargin{2};
            if (! (ischar (Learner)
                   && any (strcmpi (Learner, {'svm', 'leastsquares'}))))
              error (strcat ("RegressionLinear: 'Learner' must be either", ...
                             " 'svm' or 'leastsquares'."));
            endif
            Learner = lower (Learner);

          case 'epsilon'
            EpsilonIn = varargin{2};
            EpsilonGiven = true;
            if (! ((ischar (EpsilonIn) && strcmpi (EpsilonIn, 'auto'))
                   || (isnumeric (EpsilonIn) && isscalar (EpsilonIn)
                       && isreal (EpsilonIn) && EpsilonIn >= 0)))
              error (strcat ("RegressionLinear: 'Epsilon' must be 'auto'", ...
                             " or a nonnegative scalar."));
            endif

          case 'regularization'
            Regularization = varargin{2};
            if (! (ischar (Regularization)
                   && any (strcmpi (Regularization, {'ridge', 'lasso'}))))
              error (strcat ("RegressionLinear: 'Regularization' must be", ...
                             " either 'ridge' or 'lasso'."));
            endif
            Regularization = lower (Regularization);

          case 'lambda'
            Lambda = varargin{2};
            if (! ((ischar (Lambda) && strcmpi (Lambda, 'auto'))
                   || (isnumeric (Lambda) && isreal (Lambda)
                       && isvector (Lambda) && ! isempty (Lambda)
                       && all (Lambda >= 0) && all (isfinite (Lambda)))))
              error (strcat ("RegressionLinear: 'Lambda' must be 'auto'", ...
                             " or a vector of nonnegative finite values."));
            endif

          case 'solver'
            Solver = varargin{2};
            if (ischar (Solver))
              Solver = {Solver};
            endif
            valid = {'sgd', 'asgd', 'dual', 'bfgs', 'lbfgs', 'sparsa'};
            if (! (iscellstr (Solver) && ! isempty (Solver)
                   && all (cellfun (@(s) any (strcmpi (s, valid)), Solver))))
              error (strcat ("RegressionLinear: 'Solver' must be one of", ...
                             " 'sgd', 'asgd', 'dual', 'bfgs', 'lbfgs'", ...
                             " and 'sparsa', or a cell array of them."));
            endif
            Solver = lower (Solver);

          case 'beta'
            BetaIn = varargin{2};
            if (! (isnumeric (BetaIn) && isreal (BetaIn)
                   && ismatrix (BetaIn) && ! isempty (BetaIn)))
              error (strcat ("RegressionLinear: 'Beta' must be a real", ...
                             " numeric matrix."));
            endif

          case 'bias'
            BiasIn = varargin{2};
            if (! (isnumeric (BiasIn) && isreal (BiasIn)
                   && isvector (BiasIn) && ! isempty (BiasIn)))
              error (strcat ("RegressionLinear: 'Bias' must be a real", ...
                             " numeric vector."));
            endif

          case 'fitbias'
            FitBias = varargin{2};
            if (! (islogical (FitBias) || (isnumeric (FitBias)
                   && isscalar (FitBias) && any (FitBias == [0, 1]))))
              error (strcat ("RegressionLinear: 'FitBias' must be either", ...
                             " true or false."));
            endif
            FitBias = logical (FitBias);

          case 'postfitbias'
            PostFitBias = varargin{2};
            if (! (islogical (PostFitBias) || (isnumeric (PostFitBias)
                   && isscalar (PostFitBias) && any (PostFitBias == [0, 1]))))
              error (strcat ("RegressionLinear: 'PostFitBias' must be", ...
                             " either true or false."));
            endif
            PostFitBias = logical (PostFitBias);

          case 'observationsin'
            ObservationsIn = varargin{2};
            if (! (ischar (ObservationsIn)
                   && any (strcmpi (ObservationsIn, {'rows', 'columns'}))))
              error (strcat ("RegressionLinear: 'ObservationsIn' must be", ...
                             " either 'rows' or 'columns'."));
            endif
            ObservationsIn = lower (ObservationsIn);

          case 'betatolerance'
            BetaTolerance = varargin{2};
            if (! (isnumeric (BetaTolerance) && isscalar (BetaTolerance)
                   && isreal (BetaTolerance) && BetaTolerance >= 0))
              error (strcat ("RegressionLinear: 'BetaTolerance' must be", ...
                             " a nonnegative scalar."));
            endif

          case 'gradienttolerance'
            GradientTolerance = varargin{2};
            if (! (isnumeric (GradientTolerance)
                   && isscalar (GradientTolerance)
                   && isreal (GradientTolerance) && GradientTolerance >= 0))
              error (strcat ("RegressionLinear: 'GradientTolerance'", ...
                             " must be a nonnegative scalar."));
            endif

          case 'deltagradienttolerance'
            DeltaGradientTolerance = varargin{2};
            if (! (isnumeric (DeltaGradientTolerance)
                   && isscalar (DeltaGradientTolerance)
                   && isreal (DeltaGradientTolerance)
                   && DeltaGradientTolerance >= 0))
              error (strcat ("RegressionLinear:", ...
                             " 'DeltaGradientTolerance' must be a", ...
                             " nonnegative scalar."));
            endif

          case 'iterationlimit'
            IterationLimit = varargin{2};
            if (! (isnumeric (IterationLimit) && isscalar (IterationLimit)
                   && isreal (IterationLimit) && IterationLimit > 0
                   && fix (IterationLimit) == IterationLimit))
              error (strcat ("RegressionLinear: 'IterationLimit' must", ...
                             " be a positive integer scalar."));
            endif

          case 'passlimit'
            PassLimit = varargin{2};
            if (! (isnumeric (PassLimit) && isscalar (PassLimit)
                   && isreal (PassLimit) && PassLimit > 0
                   && fix (PassLimit) == PassLimit))
              error (strcat ("RegressionLinear: 'PassLimit' must be a", ...
                             " positive integer scalar."));
            endif

          case 'batchsize'
            BatchSize = varargin{2};
            if (! (isnumeric (BatchSize) && isscalar (BatchSize)
                   && isreal (BatchSize) && BatchSize > 0
                   && fix (BatchSize) == BatchSize))
              error (strcat ("RegressionLinear: 'BatchSize' must be a", ...
                             " positive integer scalar."));
            endif

          case 'batchlimit'
            BatchLimit = varargin{2};
            if (! (isnumeric (BatchLimit) && isscalar (BatchLimit)
                   && isreal (BatchLimit) && BatchLimit > 0
                   && fix (BatchLimit) == BatchLimit))
              error (strcat ("RegressionLinear: 'BatchLimit' must be a", ...
                             " positive integer scalar."));
            endif

          case 'learnrate'
            LearnRate = varargin{2};
            if (! (isnumeric (LearnRate) && isscalar (LearnRate)
                   && isreal (LearnRate) && LearnRate > 0))
              error (strcat ("RegressionLinear: 'LearnRate' must be a", ...
                             " positive scalar."));
            endif

          case 'optimizelearnrate'
            OptimizeLearnRate = varargin{2};
            if (! (islogical (OptimizeLearnRate)
                   || (isnumeric (OptimizeLearnRate)
                       && isscalar (OptimizeLearnRate)
                       && any (OptimizeLearnRate == [0, 1]))))
              error (strcat ("RegressionLinear: 'OptimizeLearnRate'", ...
                             " must be either true or false."));
            endif
            OptimizeLearnRate = logical (OptimizeLearnRate);

          case 'truncationperiod'
            TruncationPeriod = varargin{2};
            if (! (isnumeric (TruncationPeriod)
                   && isscalar (TruncationPeriod)
                   && isreal (TruncationPeriod) && TruncationPeriod > 0
                   && fix (TruncationPeriod) == TruncationPeriod))
              error (strcat ("RegressionLinear: 'TruncationPeriod' must", ...
                             " be a positive integer scalar."));
            endif

          case 'numcheckconvergence'
            NumCheckConvergence = varargin{2};
            if (! (isnumeric (NumCheckConvergence)
                   && isscalar (NumCheckConvergence)
                   && isreal (NumCheckConvergence)
                   && NumCheckConvergence > 0
                   && fix (NumCheckConvergence) == NumCheckConvergence))
              error (strcat ("RegressionLinear: 'NumCheckConvergence'", ...
                             " must be a positive integer scalar."));
            endif

          case 'hessianhistorysize'
            HessianHistorySize = varargin{2};
            if (! (isnumeric (HessianHistorySize)
                   && isscalar (HessianHistorySize)
                   && isreal (HessianHistorySize) && HessianHistorySize > 0
                   && fix (HessianHistorySize) == HessianHistorySize))
              error (strcat ("RegressionLinear: 'HessianHistorySize'", ...
                             " must be a positive integer scalar."));
            endif

          case 'verbose'
            Verbose = varargin{2};
            if (! (isnumeric (Verbose) && isscalar (Verbose)
                   && isreal (Verbose) && any (Verbose == [0, 1, 2])))
              error (strcat ("RegressionLinear: 'Verbose' must be 0, 1,", ...
                             " or 2."));
            endif

          case 'responsetransform'
            ResponseTransform = varargin{2};

          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isreal (Weights)
                   && isvector (Weights) && all (Weights >= 0)))
              error (strcat ("RegressionLinear: 'Weights' must be a", ...
                             " vector of nonnegative values."));
            endif

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! (iscellstr (PredictorNames) && isvector (PredictorNames)))
              error (strcat ("RegressionLinear: 'PredictorNames' must", ...
                             " be a cell array of character vectors."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! (ischar (ResponseName) && isrow (ResponseName)))
              error (strcat ("RegressionLinear: 'ResponseName' must be a", ...
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
              error (strcat ("RegressionLinear:", ...
                             " 'CategoricalPredictors' must be a vector", ...
                             " of positive integers or a logical vector."));
            endif

          otherwise
            error (strcat ("RegressionLinear: invalid parameter name in", ...
                           " optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Observations may be given down the columns, which only means the
      ## predictor matrix arrives transposed.
      if (strcmp (ObservationsIn, 'columns'))
        X = X';
      endif

      ## Validate the data and resolve the weights.  The four linear and
      ## kernel regression models share that opening, so it lives in one
      ## place.
      F = regFrame (X, Y, Weights, 'RegressionLinear');
      X = F.X;
      Y = F.Y;
      W = F.W;
      n = F.n;
      p = F.p;

      ## Epsilon belongs to the insensitive band, so it means nothing to a
      ## least squares fit and MATLAB refuses it there rather than ignoring
      ## it.
      if (strcmp (Learner, 'leastsquares'))
        if (EpsilonGiven)
          error (strcat ("RegressionLinear: 'Epsilon' applies to a", ...
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

      ## Resolve the penalty and the solver against one another, since each
      ## has a default that depends on the other.
      if (isempty (Regularization))
        if (! isempty (Solver) && numel (Solver) == 1
            && strcmp (Solver{1}, 'sparsa'))
          Regularization = 'lasso';
        else
          Regularization = 'ridge';
        endif
      endif
      if (isempty (Solver))
        if (strcmp (Regularization, 'lasso'))
          Solver = {'sparsa'};
        elseif (p <= 100)
          Solver = {'bfgs'};
        elseif (strcmp (Learner, 'svm'))
          Solver = {'dual'};
        else
          Solver = {'sgd'};
        endif
      endif
      for k = 1:numel (Solver)
        if (strcmp (Regularization, 'lasso')
            && ! any (strcmp (Solver{k}, {'sgd', 'asgd', 'sparsa'})))
          error (strcat ("RegressionLinear: the '%s' solver fits a ridge", ...
                         " penalty only."), Solver{k});
        endif
        if (strcmp (Regularization, 'ridge') && strcmp (Solver{k}, 'sparsa'))
          error (strcat ("RegressionLinear: the 'sparsa' solver fits a", ...
                         " lasso penalty only."));
        endif
        if (strcmp (Solver{k}, 'dual')
            && strcmp (Learner, 'leastsquares'))
          error (strcat ("RegressionLinear: the 'dual' solver fits an", ...
                         " epsilon-insensitive loss only, so it needs", ...
                         " 'Learner' set to 'svm'."));
        endif
      endfor

      ## Resolve Lambda, which the class reports as a number even when it
      ## was given as 'auto', and which ModelParameters keeps as given.
      LambdaIn = Lambda;
      if (ischar (Lambda))
        Lambda = 1 / n;
      else
        Lambda = sort (Lambda(:)');
      endif
      L = numel (Lambda);

      ## Resolve the starting point
      if (isempty (BetaIn))
        Beta0 = zeros (p, L);
      else
        if (rows (BetaIn) != p)
          error (strcat ("RegressionLinear: 'Beta' must have one row per", ...
                         " predictor."));
        endif
        if (columns (BetaIn) == 1)
          Beta0 = repmat (BetaIn, 1, L);
        elseif (columns (BetaIn) == L)
          Beta0 = BetaIn;
        else
          error (strcat ("RegressionLinear: 'Beta' must have one column,", ...
                         " or one per value of 'Lambda'."));
        endif
      endif
      if (isempty (BiasIn))
        if (strcmp (Learner, 'leastsquares'))
          Bias0 = repmat (sum (W .* Y), 1, L);
        else
          Bias0 = repmat (weightedMedian (Y, W), 1, L);
        endif
      else
        BiasIn = BiasIn(:)';
        if (numel (BiasIn) == 1)
          Bias0 = repmat (BiasIn, 1, L);
        elseif (numel (BiasIn) == L)
          Bias0 = BiasIn;
        else
          error (strcat ("RegressionLinear: 'Bias' must be a scalar, or", ...
                         " hold one value per value of 'Lambda'."));
        endif
      endif

      if (isempty (PassLimit))
        if (any (strcmp (Solver, 'dual')))
          PassLimit = 10;
        else
          PassLimit = 1;
        endif
      endif
      if (isempty (LearnRate))
        ## MATLAB's default: the reciprocal root of one plus the largest
        ## squared length of an observation, so a step never overshoots the
        ## widest row of the data.
        LearnRate = 1 / sqrt (1 + max (sum (X .^ 2, 2)));
      endif

      ## Fit one model per regularization strength, each warm starting the
      ## next, which is what makes an ascending Lambda cheaper than the same
      ## values fitted apart.
      Beta = zeros (p, L);
      Bias = zeros (1, L);
      info = struct ([]);
      P = struct ();
      P.Learner = Learner;
      P.LossFunction = 'epsiloninsensitive';
      if (strcmp (Learner, 'leastsquares'))
        P.LossFunction = 'mse';
      endif
      P.Epsilon = Epsilon;
      P.Regularization = Regularization;
      P.FitBias = FitBias;
      P.PostFitBias = PostFitBias;
      P.BetaTolerance = BetaTolerance;
      P.GradientTolerance = GradientTolerance;
      P.DeltaGradientTolerance = DeltaGradientTolerance;
      P.IterationLimit = IterationLimit;
      P.PassLimit = PassLimit;
      P.BatchSize = BatchSize;
      P.BatchLimit = BatchLimit;
      P.LearnRate = LearnRate;
      P.OptimizeLearnRate = OptimizeLearnRate;
      P.TruncationPeriod = TruncationPeriod;
      P.NumCheckConvergence = NumCheckConvergence;
      P.HessianHistorySize = HessianHistorySize;

      for l = 1:L
        P.Lambda = Lambda(l);
        b = Beta0(:,l);
        b0 = Bias0(l);
        if (l > 1 && isempty (BetaIn))
          b = Beta(:,l-1);
          b0 = Bias(l-1);
        endif
        for k = 1:numel (Solver)
          P.Solver = Solver{k};
          P.InitialBeta = b;
          P.InitialBias = b0;
          [b, b0, S] = linearSolve (X, Y, W, P);
        endfor
        Beta(:,l) = b;
        Bias(l) = b0;
        S.Lambda = Lambda(l);
        if (l == 1)
          info = S;
        else
          info(l) = S;
        endif
      endfor

      ## Fill in the model
      this.Epsilon = Epsilon;
      this.PredictorNames = PredictorNames;
      if (isempty (this.PredictorNames))
        this.PredictorNames = ...
                     arrayfun (@(k) sprintf ("x%d", k), 1:p, ...
                               'UniformOutput', false);
      elseif (numel (this.PredictorNames) != p)
        error (strcat ("RegressionLinear: 'PredictorNames' must have one", ...
                       " name per predictor."));
      endif
      this.ExpandedPredictorNames = this.PredictorNames;
      this.CategoricalPredictors = CategoricalPredictors;
      this.ResponseName = ResponseName;
      this.Learner = Learner;
      this.Beta = Beta;
      this.Bias = Bias;
      this.FittedLoss = P.LossFunction;
      this.Lambda = Lambda;
      this.NumPredictors_ = p;
      this.ResponseTransform = ResponseTransform;
      if (strcmp (Regularization, 'ridge'))
        this.Regularization = 'ridge (L2)';
      else
        this.Regularization = 'lasso (L1)';
      endif

      this.ModelParameters = linearModelParams (P, Solver, LambdaIn, ...
                                                EpsilonIn, Beta0, Bias0, ...
                                                Verbose, 'regression');

      this.FitInfo_ = linearFitInfo (info, BetaTolerance, ...
                                     GradientTolerance, ...
                                     DeltaGradientTolerance, ...
                                     IterationLimit, Solver, PassLimit, ...
                                     BatchLimit);

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionLinear} {@var{yFit} =} predict (@var{obj}, @var{XC})
    ##
    ## Predict the response of new observations.
    ##
    ## @code{@var{yFit} = predict (@var{obj}, @var{XC})} returns one
    ## predicted value per row of @var{XC}, and one column per regularization
    ## strength.  @qcode{ResponseTransform} is applied to the result.
    ##
    ## @end deftypefn
    function yFit = predict (this, XC)

      if (nargin < 2)
        error ("RegressionLinear.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("RegressionLinear.predict: XC is empty.");
      endif
      if (! (isnumeric (XC) && isreal (XC) && ismatrix (XC)))
        error ("RegressionLinear.predict: invalid values in XC.");
      endif
      if (columns (XC) != this.NumPredictors_)
        error (strcat ("RegressionLinear.predict: XC must have the same", ...
                       " number of predictors as the trained model."));
      endif

      yFit = this.RTfun (XC * this.Beta + this.Bias);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionLinear} {@var{l} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionLinear} {@var{l} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss on new data.
    ##
    ## @code{@var{l} = loss (@var{obj}, @var{X}, @var{Y})} returns the mean
    ## squared error, one value per regularization strength.
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
        error ("RegressionLinear.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionLinear.loss: optional arguments must", ...
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
              error (strcat ("RegressionLinear.loss: 'LossFun' must be", ...
                             " either 'mse' or 'epsiloninsensitive'."));
            endif
            LossFun = lower (LossFun);
          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isreal (Weights)
                   && isvector (Weights) && all (Weights >= 0)))
              error (strcat ("RegressionLinear.loss: 'Weights' must be", ...
                             " a vector of nonnegative values."));
            endif
          otherwise
            error (strcat ("RegressionLinear.loss: invalid parameter", ...
                           " name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      if (strcmp (LossFun, 'epsiloninsensitive') && isempty (this.Epsilon))
        error (strcat ("RegressionLinear.loss: the", ...
                       " 'epsiloninsensitive' loss applies to a support", ...
                       " vector machine only."));
      endif

      Y = Y(:);
      if (! (isnumeric (Y) && isreal (Y)))
        error ("RegressionLinear.loss: invalid values in Y.");
      endif
      if (rows (X) != numel (Y))
        error (strcat ("RegressionLinear.loss: number of rows in X and Y", ...
                       " must be equal."));
      endif
      if (isempty (Weights))
        w = ones (numel (Y), 1);
      else
        w = Weights(:);
        if (numel (w) != numel (Y))
          error (strcat ("RegressionLinear.loss: 'Weights' must have one", ...
                         " element per observation."));
        endif
      endif
      w = w / sum (w);

      yFit = predict (this, X);
      L = numel (this.Lambda);
      l = zeros (1, L);
      for k = 1:L
        r = Y - yFit(:,k);
        if (strcmp (LossFun, 'mse'))
          l(k) = sum (w .* (r .^ 2));
        else
          l(k) = sum (w .* max (0, abs (r) - this.Epsilon));
        endif
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionLinear} {@var{sub} =} selectModels (@var{obj}, @var{idx})
    ##
    ## Keep a subset of the fitted regularization strengths.
    ##
    ## @code{@var{sub} = selectModels (@var{obj}, @var{idx})} returns a model
    ## holding only the strengths @var{idx} names, which may be indices into
    ## @qcode{Lambda} or a logical vector over it.
    ##
    ## @end deftypefn
    function sub = selectModels (this, idx)

      if (nargin < 2)
        error (strcat ("RegressionLinear.selectModels: too few input", ...
                       " arguments."));
      endif
      L = numel (this.Lambda);
      if (islogical (idx))
        if (numel (idx) != L)
          error (strcat ("RegressionLinear.selectModels: a logical IDX", ...
                         " must have one element per value of 'Lambda'."));
        endif
        idx = find (idx);
      endif
      if (! (isnumeric (idx) && isreal (idx) && isvector (idx)
             && ! isempty (idx) && all (fix (idx) == idx)
             && all (idx >= 1) && all (idx <= L)))
        error (strcat ("RegressionLinear.selectModels: IDX must hold", ...
                       " integers between 1 and %d."), L);
      endif

      sub = this;
      sub.Lambda = this.Lambda(idx);
      sub.Beta = this.Beta(:,idx);
      sub.Bias = this.Bias(idx);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionLinear} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a linear regression model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves the model
    ## @var{obj} into @var{filename} in a form @code{loadmodel} can read
    ## back.
    ##
    ## @end deftypefn
    function savemodel (obj, fname)

      classdef_name = 'RegressionLinear';
      Epsilon = obj.Epsilon;
      ResponseTransform = obj.ResponseTransform;
      PredictorNames = obj.PredictorNames;
      CategoricalPredictors = obj.CategoricalPredictors;
      ResponseName = obj.ResponseName;
      ExpandedPredictorNames = obj.ExpandedPredictorNames;
      Learner = obj.Learner;
      Beta = obj.Beta;
      Bias = obj.Bias;
      FittedLoss = obj.FittedLoss;
      Lambda = obj.Lambda;
      ModelParameters = obj.ModelParameters;
      Regularization = obj.Regularization;
      NumPredictors_ = obj.NumPredictors_;

      save ('-binary', fname, 'classdef_name', 'Epsilon', ...
            'ResponseTransform', 'PredictorNames', ...
            'CategoricalPredictors', 'ResponseName', ...
            'ExpandedPredictorNames', 'Learner', 'Beta', 'Bias', ...
            'FittedLoss', 'Lambda', 'ModelParameters', 'Regularization', ...
            'NumPredictors_');

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
      printf ("\n  RegressionLinear\n\n");
      printf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      printf ("%+25s: '%s'\n", 'ResponseTransform', this.ResponseTransform);
      printf ("%+25s: [%dx%d double]\n", 'Beta', rows (this.Beta), ...
              columns (this.Beta));
      if (numel (this.Bias) == 1)
        printf ("%+25s: %g\n", 'Bias', this.Bias);
        printf ("%+25s: %g\n", 'Lambda', this.Lambda);
      else
        printf ("%+25s: [1x%d double]\n", 'Bias', numel (this.Bias));
        printf ("%+25s: [1x%d double]\n", 'Lambda', numel (this.Lambda));
      endif
      printf ("%+25s: '%s'\n", 'Learner', this.Learner);
      printf ("\n");
    endfunction

    ## What the fit reported, which fitrlinear returns as its second output.
    function S = fitInfo_ (this)
      S = this.FitInfo_;
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      mdl = RegressionLinear ([0; 1], [0; 1]);
      fields = fieldnames (data);
      for k = 1:numel (fields)
        mdl.(fields{k}) = data.(fields{k});
      endfor

    endfunction

  endmethods

endclassdef

%!shared X, Y
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! X = X(ok,:);
%! Y = MPG(ok);

%!demo
%! ## Fit fuel consumption on four engine measurements and read the
%! ## coefficients and the insensitive band the fit chose for itself.
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! Mdl = RegressionLinear (X(ok,:), MPG(ok))
%! yFit = predict (Mdl, X(find (ok, 3),:))

%!demo
%! ## Least squares has no insensitive band, so Epsilon is empty, and its
%! ## loss is the mean squared error the fit minimizes.
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! Mdl = RegressionLinear (X(ok,:), MPG(ok), 'Learner', 'leastsquares');
%! band = Mdl.Epsilon
%! mse = loss (Mdl, X(ok,:), MPG(ok))

%!test
%! ## The model reports the surface MATLAB reports
%! Mdl = RegressionLinear (X, Y);
%! assert_equal (class (Mdl), 'RegressionLinear');
%! assert_equal (Mdl.Learner, 'svm');
%! assert_equal (Mdl.FittedLoss, 'epsiloninsensitive');
%! assert_equal (Mdl.Regularization, 'ridge (L2)');
%! assert_equal (Mdl.ResponseTransform, 'none');
%! assert_equal (Mdl.ResponseName, 'Y');
%! assert_equal (Mdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (Mdl.ExpandedPredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (size (Mdl.Beta), [4, 1]);

%!test
%! ## The properties are the ones MATLAB lists, in its order
%! Mdl = RegressionLinear (X, Y);
%! assert_equal (sort (properties (Mdl)), ...
%!               sort ({'Epsilon'; 'ResponseTransform'; 'PredictorNames'; ...
%!                      'CategoricalPredictors'; 'ResponseName'; ...
%!                      'ExpandedPredictorNames'; 'Learner'; 'Beta'; ...
%!                      'Bias'; 'FittedLoss'; 'Lambda'; 'ModelParameters'; ...
%!                      'Regularization'}));

%!test
%! ## A least squares ridge fit reproduces R2024a's coefficients
%! Mdl = RegressionLinear (X, Y, 'Learner', 'leastsquares', ...
%!                         'Solver', 'lbfgs', 'BetaTolerance', 0, ...
%!                         'GradientTolerance', 1e-12, ...
%!                         'IterationLimit', 20000);
%! assert_equal (Mdl.Beta, [-0.060147206634135; -0.00667928082241265; ...
%!                          -0.0375377897819825; -0.00608459830980712], 1e-7);
%! assert_equal (Mdl.Bias, 48.1149099265466, 1e-6);
%! assert_equal (Mdl.FittedLoss, 'mse');
%! assert_equal (Mdl.Epsilon, []);
%! assert_equal (loss (Mdl, X, Y), 15.9484134801834, 1e-6);

%!test
%! ## Epsilon defaults to the interquartile range over 13.49, R2024a's own
%! ## estimate of the standard deviation of the response
%! Mdl = RegressionLinear (X, Y);
%! assert_equal (Mdl.Epsilon, 0.926612305411416, 1e-12);
%! assert_equal (Mdl.Epsilon, iqr (Y) / 13.49, 1e-15);

%!test
%! ## A response of no spread has no interquartile range to scale, and the
%! ## band falls back on a tenth
%! Mdl = RegressionLinear (ones (10, 2), 5 * ones (10, 1));
%! assert_equal (Mdl.Epsilon, 0.1);

%!test
%! ## The epsilon-insensitive loss is not differentiable at the edge of the
%! ## band, so the line search gives up short of a minimum and where it stops
%! ## follows the last bits: 2.8273 here, 3.4883 under clang and on macOS,
%! ## against the 2.82717018548797 R2024a reaches.  Assert the identity the
%! ## reported objective holds at whichever point it stops.
%! Mdl = RegressionLinear (X, Y, 'Learner', 'svm', 'Solver', 'lbfgs', ...
%!                         'BetaTolerance', 0, 'GradientTolerance', 1e-12, ...
%!                         'IterationLimit', 20000);
%! assert_equal (Mdl.FitInfo_.Objective, ...
%!               loss (Mdl, X, Y, 'LossFun', 'epsiloninsensitive') ...
%!               + 0.5 * Mdl.Lambda * sum (Mdl.Beta .^ 2), 1e-12);

%!test
%! ## Lambda defaults to the reciprocal of the observations that were used
%! Mdl = RegressionLinear (X, Y);
%! assert_equal (Mdl.Lambda, 1 / 93, 1e-15);

%!test
%! ## A vector of strengths fits one model per value, sorted ascending, and
%! ## every method reports one column per value
%! Mdl = RegressionLinear (X, Y, 'Lambda', [0.1, 0.001, 0.01]);
%! assert_equal (Mdl.Lambda, [0.001, 0.01, 0.1]);
%! assert_equal (size (Mdl.Beta), [4, 3]);
%! assert_equal (size (Mdl.Bias), [1, 3]);
%! assert_equal (size (predict (Mdl, X(1:4,:))), [4, 3]);
%! assert_equal (size (loss (Mdl, X, Y)), [1, 3]);

%!test
%! ## selectModels keeps the strengths it is given and drops the rest
%! Mdl = RegressionLinear (X, Y, 'Lambda', [0.001, 0.01, 0.1]);
%! sub = selectModels (Mdl, [1, 3]);
%! assert_equal (sub.Lambda, [0.001, 0.1]);
%! assert_equal (sub.Beta, Mdl.Beta(:,[1, 3]));

%!test
%! ## A lasso penalty drives coefficients to exactly zero
%! Mdl = RegressionLinear (X, Y, 'Learner', 'leastsquares', ...
%!                         'Regularization', 'lasso', 'Lambda', 100);
%! assert_equal (Mdl.Regularization, 'lasso (L1)');
%! assert_equal (sum (Mdl.Beta == 0) > 0, true);

%!test
%! ## With the penalty all but switched off, a least squares fit is ordinary
%! ## least squares, whichever solver ran.  This is what pins the solvers to
%! ## a value that is known independently of either of them.
%! b = [ones(rows (X), 1), X] \ Y;
%! Ml = RegressionLinear (X, Y, 'Learner', 'leastsquares', ...
%!                        'Regularization', 'lasso', 'Lambda', 1e-8, ...
%!                        'BetaTolerance', 0, 'GradientTolerance', 1e-14, ...
%!                        'IterationLimit', 200000);
%! Mr = RegressionLinear (X, Y, 'Learner', 'leastsquares', ...
%!                        'Solver', 'lbfgs', 'Lambda', 1e-8, ...
%!                        'BetaTolerance', 0, 'GradientTolerance', 1e-14, ...
%!                        'IterationLimit', 200000);
%! ## Relative, so every coefficient is held to the same number of figures:
%! ## an absolute tolerance is far stricter on the largest of them, and how
%! ## close a solver stops to the closed form answer varies with the BLAS.
%! assert_equal (Ml.Beta, b(2:end), -1e-3);
%! assert_equal (Ml.Bias, b(1), 1e-4);
%! assert_equal (Mr.Beta, b(2:end), -1e-3);
%! assert_equal (Mr.Bias, b(1), 1e-4);

%!test
%! ## A response transform reaches predict
%! Mdl = RegressionLinear (X, Y, 'Learner', 'leastsquares');
%! plain = predict (Mdl, X(1:5,:));
%! Mdl.ResponseTransform = 'exp';
%! assert_equal (Mdl.ResponseTransform, 'exp');
%! assert_equal (predict (Mdl, X(1:5,:)), exp (plain), 1e-12);

%!test
%! ## FitBias false leaves the intercept at zero
%! Mdl = RegressionLinear (X, Y, 'FitBias', false);
%! assert_equal (Mdl.Bias, 0);

%!test
%! ## Observations may be given down the columns instead
%! Mr = RegressionLinear (X, Y, 'Learner', 'leastsquares');
%! Mc = RegressionLinear (X', Y, 'Learner', 'leastsquares', ...
%!                        'ObservationsIn', 'columns');
%! assert_equal (Mr.Beta, Mc.Beta, 1e-12);

%!test
%! ## A row with a missing predictor or a missing response is dropped, and
%! ## Lambda follows the count that survived
%! Xn = X;
%! Xn(3,2) = NaN;
%! Mdl = RegressionLinear (Xn, Y);
%! assert_equal (Mdl.Lambda, 1 / 92, 1e-15);

%!test
%! ## A saved model reads back as the same model
%! Mdl = RegressionLinear (X, Y, 'Learner', 'leastsquares');
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mnew = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (Mnew), 'RegressionLinear');
%! assert_equal (Mnew.Beta, Mdl.Beta);
%! assert_equal (predict (Mnew, X(1:5,:)), predict (Mdl, X(1:5,:)));

%!test
%! ## LossTolerance: the engine's loss test is on the objective VALUE, not on
%! ## its change, so linearSolve must switch it off with -Inf.  This fit is
%! ## the one that catches it: an exact linear relation at a vanishing
%! ## penalty drives the objective through 1e-6 long before the coefficients
%! ## are right, and with the engine's default the fit stops five iterations
%! ## early with them out by 1.6e-4.  If this test starts failing, look at
%! ## opt.LossTolerance in linearSolve before anything else.
%! randn ('seed', 7);
%! Xe = randn (60, 3);
%! btrue = [2; -3; 0.5];
%! Ye = Xe * btrue + 4;
%! Mdl = RegressionLinear (Xe, Ye, 'Learner', 'leastsquares', ...
%!                         'Solver', 'lbfgs', 'Lambda', 1e-10, ...
%!                         'BetaTolerance', 0, 'GradientTolerance', 1e-14, ...
%!                         'IterationLimit', 20000);
%! assert_equal (Mdl.Beta, btrue, 1e-8);
%! assert_equal (Mdl.Bias, 4, 1e-8);

%!test
%! ## HessianHistorySize reaches the solver rather than only being recorded:
%! ## a one-pair history takes several times the iterations a fifteen-pair
%! ## one does on the same problem.  Fifteen is MATLAB's documented default
%! ## for fitrlinear, where the engine's own default is ten.
%! bf = RegressionLinear (X, Y);
%! assert_equal (bf.ModelParameters.Solver, {'bfgs'});
%! assert_equal (bf.ModelParameters.HessianHistorySize, []);
%! lb = RegressionLinear (X, Y, 'Solver', 'lbfgs');
%! assert_equal (lb.ModelParameters.HessianHistorySize, 15);
%! opts = {'Learner', 'leastsquares', 'Solver', 'lbfgs', 'BetaTolerance', 0, ...
%!         'GradientTolerance', 1e-10, 'IterationLimit', 5000};
%! short = RegressionLinear (X, Y, opts{:}, 'HessianHistorySize', 1);
%! long = RegressionLinear (X, Y, opts{:}, 'HessianHistorySize', 15);
%! assert_equal (short.FitInfo_.NumIterations ...
%!               > long.FitInfo_.NumIterations, true);

%!test
%! ## The regression counterpart keeps the documented 0.1 where the
%! ## classifier takes 1, confirmed on R2026a; the convergence check count
%! ## is 2 on both, where the documentation says 5
%! Mdl = RegressionLinear (X, Y, 'Solver', 'dual');
%! assert_equal (Mdl.ModelParameters.DeltaGradientTolerance, 0.1);
%! assert_equal (Mdl.ModelParameters.NumCheckConvergence, 2);
%! assert_equal (Mdl.ModelParameters.PassLimit, 10);

## Test input validation
%!error<RegressionLinear: too few input arguments.> RegressionLinear (ones (5, 2))
%!error<RegressionLinear: optional arguments must be given in Name-Value pairs.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Learner')
%!error<RegressionLinear: 'Learner' must be either 'svm' or 'leastsquares'.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Learner', 'logistic')
%!error<RegressionLinear: 'Epsilon' must be 'auto' or a nonnegative scalar.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Epsilon', -1)
%!error<RegressionLinear: 'Epsilon' applies to a support vector machine only.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Learner', 'leastsquares', ...
%!                     'Epsilon', 1)
%!error<RegressionLinear: 'Regularization' must be either 'ridge' or 'lasso'.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Regularization', 'elastic')
%!error<RegressionLinear: 'Lambda' must be 'auto' or a vector of nonnegative finite values.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Lambda', Inf)
%!error<RegressionLinear: 'Solver' must be one of 'sgd', 'asgd', 'dual', 'bfgs', 'lbfgs' and 'sparsa', or a cell array of them.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Solver', 'newton')
%!error<RegressionLinear: 'IterationLimit' must be a positive integer scalar.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'IterationLimit', 0)
%!error<RegressionLinear: invalid parameter name in optional pair arguments.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Nonsense', 1)
%!error<RegressionLinear: invalid values in X.> RegressionLinear ({1, 2; 3, 4}, [1; 2])
%!error<RegressionLinear: X is empty.> RegressionLinear ([], [])
%!error<RegressionLinear: invalid values in Y.> ...
%! RegressionLinear (ones (10, 2), {1, 2})
%!error<RegressionLinear: number of rows in X and Y must be equal.> ...
%! RegressionLinear (ones (10, 2), ones (3, 1))
%!error<RegressionLinear: 'Weights' must have one element per observation.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Weights', ones (3, 1))
%!error<RegressionLinear: the 'sparsa' solver fits a lasso penalty only.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Regularization', 'ridge', ...
%!                     'Solver', 'sparsa')
%!error<RegressionLinear: the 'lbfgs' solver fits a ridge penalty only.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Regularization', 'lasso', ...
%!                     'Solver', 'lbfgs')
%!error<RegressionLinear: the 'dual' solver fits an epsilon-insensitive loss only, so it needs 'Learner' set to 'svm'.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'Learner', 'leastsquares', ...
%!                     'Solver', 'dual')
%!error<RegressionLinear: 'PredictorNames' must have one name per predictor.> ...
%! RegressionLinear (ones (10, 2), ones (10, 1), 'PredictorNames', {'a'})
%!error<RegressionLinear.predict: too few input arguments.> ...
%! predict (RegressionLinear (ones (10, 2), ones (10, 1)))
%!error<RegressionLinear.predict: XC is empty.> ...
%! predict (RegressionLinear (ones (10, 2), ones (10, 1)), [])
%!error<RegressionLinear.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (RegressionLinear (ones (10, 2), ones (10, 1)), ones (3, 5))
%!error<RegressionLinear.loss: too few input arguments.> ...
%! loss (RegressionLinear (ones (10, 2), ones (10, 1)), ones (10, 2))
%!error<RegressionLinear.loss: 'LossFun' must be either 'mse' or 'epsiloninsensitive'.> ...
%! loss (RegressionLinear (ones (10, 2), ones (10, 1)), ones (10, 2), ...
%!                     ones (10, 1), 'LossFun', 'hinge')
%!error<RegressionLinear.loss: the 'epsiloninsensitive' loss applies to a support vector machine only.> ...
%! loss (RegressionLinear (ones (10, 2), ones (10, 1), 'Learner', ...
%!                     'leastsquares'), ones (10, 2), ones (10, 1), ...
%!                     'LossFun', 'epsiloninsensitive')
%!error<RegressionLinear.loss: number of rows in X and Y must be equal.> ...
%! loss (RegressionLinear (ones (10, 2), ones (10, 1)), ones (10, 2), ...
%!                     ones (3, 1))
%!error<RegressionLinear.selectModels: IDX must hold integers between 1 and 2.> ...
%! selectModels (RegressionLinear (ones (10, 2), ones (10, 1), 'Lambda', ...
%!                     [0.1, 0.2]), 5)

## Every documented response transform reaches the response that is reported.
%!test
%! load fisheriris
%! Mdl = fitrlinear (meas(:,2:4), meas(:,1));
%! Mdl.ResponseTransform = 'none';
%! raw = predict (Mdl, meas([1, 60, 120],2:4));
%! T = {'identity', @(x) x; 'exp', @(x) exp (x); 'log', @(x) log (x)};
%! for i = 1:rows (T)
%!   Mdl.ResponseTransform = T{i,1};
%!   yhat = predict (Mdl, meas([1, 60, 120],2:4));
%!   assert_equal (yhat, T{i,2}(raw), 1e-12);
%! endfor

## A function handle is taken as given and applied to the response.
%!test
%! load fisheriris
%! Mdl = fitrlinear (meas(:,2:4), meas(:,1));
%! Mdl.ResponseTransform = 'none';
%! raw = predict (Mdl, meas([1, 60, 120],2:4));
%! Mdl.ResponseTransform = @(x) x .^ 2;
%! yhat = predict (Mdl, meas([1, 60, 120],2:4));
%! assert_equal (yhat, raw .^ 2, 1e-12);
