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
## @deftp {statistics} ClassificationLinear
##
## Linear binary classifier for high dimensional data.
##
## A @qcode{ClassificationLinear} object fits a linear model,
## @code{@var{X} * Beta + Bias}, to a two class problem by minimizing a
## regularized average loss.  The loss is the hinge loss for a support vector
## machine and the deviance for a logistic regression, and the penalty is
## either a ridge or a lasso one.
##
## Unlike the other classifiers of this package the object holds no copy of
## the training data: the coefficients, the intercept and the fitting options
## are the whole model.  That is what makes it suited to data with more
## predictors than an in memory kernel matrix could carry, and it is why the
## class has no @code{compact} method and no resubstitution methods.
##
## A vector of regularization strengths fits one model per value in a single
## object.  @qcode{Beta} is then a @math{PxL} matrix and @qcode{Bias} a
## @math{1xL} row, every method returns one column per strength, and
## @code{selectModels} narrows the object down to the strengths worth
## keeping.
##
## Create a @qcode{ClassificationLinear} object with @code{fitclinear}.
##
## @seealso{fitclinear, ClassificationKernel, ClassificationSVM}
## @end deftp

classdef ClassificationLinear

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationLinear} {property} ClassNames
    ##
    ## Names of the two classes
    ##
    ## A column of the same type as the response supplied to the constructor:
    ## a cell array of character vectors, a numeric vector, a logical vector
    ## or a character matrix.  The second of the two is the positive class,
    ## the one a positive score belongs to.  This property is read-only.
    ##
    ## @end deftp
    ClassNames             = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationLinear} {property} Prior
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
    ## @deftp {ClassificationLinear} {property} Cost
    ##
    ## Cost of misclassifying an observation
    ##
    ## A square numeric matrix with one row and one column per class, whose
    ## @math{(i,j)} element is the cost of classifying an observation of
    ## class @math{i} into class @math{j}.  It defaults to one everywhere
    ## except the diagonal, which is zero.  This property is read-only:
    ## MATLAB refuses an assignment into it on this class, as it does on the
    ## support vector machine, so a cost matrix is given to the constructor
    ## instead.
    ##
    ## The cost matrix takes no part in the fit and none in @code{predict},
    ## which returns the class of largest score.  It is read by the
    ## @qcode{'mincost'} and @qcode{'classifcost'} losses alone.
    ##
    ## @end deftp
    Cost                   = [];

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationLinear} {property} ScoreTransform
    ##
    ## Transformation applied to the predicted scores
    ##
    ## A character vector naming a transformation, or the text of the
    ## function handle that was supplied.  Assigning to it accepts either.
    ## It defaults to @qcode{'logit'} for a logistic learner, which turns the
    ## scores into posterior probabilities, and to @qcode{'none'} for a
    ## support vector machine.
    ##
    ## @end deftp
    ScoreTransform         = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationLinear} {property} PredictorNames
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
    ## @deftp {ClassificationLinear} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A row vector of column indices, empty when every predictor is
    ## numeric.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationLinear} {property} ResponseName
    ##
    ## Name of the response
    ##
    ## A character vector, defaulting to @qcode{'Y'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    ResponseName           = 'Y';

    ## -*- texinfo -*-
    ## @deftp {ClassificationLinear} {property} ExpandedPredictorNames
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
    ## @deftp {ClassificationLinear} {property} Learner
    ##
    ## Linear classification model that was fitted
    ##
    ## Either @qcode{'svm'} or @qcode{'logistic'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Learner                = 'svm';

    ## -*- texinfo -*-
    ## @deftp {ClassificationLinear} {property} Beta
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
    ## @deftp {ClassificationLinear} {property} Bias
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
    ## @deftp {ClassificationLinear} {property} FittedLoss
    ##
    ## Loss function the fit minimized
    ##
    ## @qcode{'hinge'} for a support vector machine and @qcode{'logit'} for a
    ## logistic regression.  This is the loss of the objective, which is not
    ## the loss @code{loss} reports unless it is asked for.  This property is
    ## read-only.
    ##
    ## @end deftp
    FittedLoss             = 'hinge';

    ## -*- texinfo -*-
    ## @deftp {ClassificationLinear} {property} Lambda
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
    ## @deftp {ClassificationLinear} {property} ModelParameters
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
    ## @deftp {ClassificationLinear} {property} Regularization
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

    ## The callable behind ScoreTransform.  The public property is the text
    ## MATLAB stores; this is what predict actually applies.
    STfun                  = @(s) s;

    ## Number of predictors the fit saw, which the class does not report
    ## because MATLAB does not, but predict needs to validate its input.
    NumPredictors_         = [];

    ## What the fit reported, so that fitclinear can hand it back as its
    ## second output.  MATLAB returns it from the fitting function rather
    ## than storing it on the model, and a constructor cannot have a second
    ## output argument.
    FitInfo_               = [];

  endproperties

  methods (Access = public)

    ## Custom setter, so that assigning a name or a handle updates both the
    ## text the property reports and the callable predict uses.
    function this = set.ScoreTransform (this, val)
      [this.STfun, this.ScoreTransform] = ...
                     parseScoreTransform (val, 'ClassificationLinear');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationLinear} {@var{obj} =} ClassificationLinear (@var{X}, @var{Y})
    ## @deftypefnx {ClassificationLinear} {@var{obj} =} ClassificationLinear (@dots{}, @var{name}, @var{value})
    ##
    ## Fit a linear binary classifier.
    ##
    ## @code{@var{obj} = ClassificationLinear (@var{X}, @var{Y})} fits a
    ## linear support vector machine to the @math{NxP} predictor matrix
    ## @var{X} and the @math{Nx1} response @var{Y}, which must name exactly
    ## two classes.
    ##
    ## @code{@var{obj} = ClassificationLinear (@dots{}, @var{name},
    ## @var{value})} takes the following @qcode{Name-Value} pairs.
    ##
    ## @multitable @columnfractions 0.28 0.72
    ## @headitem Name @tab Value
    ##
    ## @item @qcode{'Learner'} @tab @qcode{'svm'}, the default, or
    ## @qcode{'logistic'}.  The first minimizes the hinge loss and the second
    ## the deviance.
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
    ## The default depends on the data and the penalty, as described below.
    ##
    ## @item @qcode{'Beta'} @tab Initial coefficients, a @math{Px1} column or
    ## a @math{PxL} matrix.  It defaults to zeros.
    ##
    ## @item @qcode{'Bias'} @tab Initial intercept, a scalar or a @math{1xL}
    ## row.  It defaults to the weighted average of the class labels for a
    ## logistic learner and to zero for a support vector machine.
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
    ## complementarity gap of the @qcode{'dual'} solver, @qcode{1} by
    ## default for a hinge loss.  MathWorks documents @qcode{0.1}, which is
    ## the default of the @emph{regression} counterpart; R2024a and R2026a
    ## both report @qcode{1} here.
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
    ## report @qcode{2}, so the documentation is stale rather than the
    ## releases being inconsistent.
    ##
    ## @item @qcode{'HessianHistorySize'} @tab Number of curvature pairs the
    ## quasi-Newton solvers keep, @qcode{15} by default.
    ##
    ## @item @qcode{'ClassNames'} @tab The classes to keep, given in the type
    ## of @var{Y}.  Observations of any other class are dropped.
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
    ## The default solver is @qcode{'sparsa'} under a lasso penalty.  Under a
    ## ridge penalty it is @qcode{'bfgs'} when there are no more than 100
    ## predictors, and beyond that @qcode{'dual'} for a support vector
    ## machine and @qcode{'sgd'} for a logistic regression.
    ##
    ## @seealso{fitclinear, ClassificationKernel}
    ## @end deftypefn
    function this = ClassificationLinear (X, Y, varargin)

      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationLinear: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationLinear: optional arguments must", ...
                       " be given in Name-Value pairs."));
      endif

      ## Defaults, before the optional arguments are parsed.  Anything left
      ## empty here is resolved once the data is known.
      Learner                = 'svm';
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
      DeltaGradientTolerance = 1;
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
      ClassNames             = [];
      CostIn                 = [];
      Prior                  = [];
      ScoreTransform         = [];
      Weights                = [];
      PredictorNames         = {};
      ResponseName           = 'Y';
      CategoricalPredictors  = [];

      ## Parse optional parameters
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))

          case 'learner'
            Learner = varargin{2};
            if (! (ischar (Learner)
                   && any (strcmpi (Learner, {'svm', 'logistic'}))))
              error (strcat ("ClassificationLinear: 'Learner' must be", ...
                             " either 'svm' or 'logistic'."));
            endif
            Learner = lower (Learner);

          case 'regularization'
            Regularization = varargin{2};
            if (! (ischar (Regularization)
                   && any (strcmpi (Regularization, {'ridge', 'lasso'}))))
              error (strcat ("ClassificationLinear: 'Regularization'", ...
                             " must be either 'ridge' or 'lasso'."));
            endif
            Regularization = lower (Regularization);

          case 'lambda'
            Lambda = varargin{2};
            if (! ((ischar (Lambda) && strcmpi (Lambda, 'auto'))
                   || (isnumeric (Lambda) && isreal (Lambda)
                       && isvector (Lambda) && ! isempty (Lambda)
                       && all (Lambda >= 0) && all (isfinite (Lambda)))))
              error (strcat ("ClassificationLinear: 'Lambda' must be", ...
                             " 'auto' or a vector of nonnegative", ...
                             " finite values."));
            endif

          case 'solver'
            Solver = varargin{2};
            if (ischar (Solver))
              Solver = {Solver};
            endif
            valid = {'sgd', 'asgd', 'dual', 'bfgs', 'lbfgs', 'sparsa'};
            if (! (iscellstr (Solver) && ! isempty (Solver)
                   && all (cellfun (@(s) any (strcmpi (s, valid)), Solver))))
              error (strcat ("ClassificationLinear: 'Solver' must be one", ...
                             " of 'sgd', 'asgd', 'dual', 'bfgs',", ...
                             " 'lbfgs' and 'sparsa', or a cell array of", ...
                             " them."));
            endif
            Solver = lower (Solver);

          case 'beta'
            BetaIn = varargin{2};
            if (! (isnumeric (BetaIn) && isreal (BetaIn)
                   && ismatrix (BetaIn) && ! isempty (BetaIn)))
              error (strcat ("ClassificationLinear: 'Beta' must be a", ...
                             " real numeric matrix."));
            endif

          case 'bias'
            BiasIn = varargin{2};
            if (! (isnumeric (BiasIn) && isreal (BiasIn)
                   && isvector (BiasIn) && ! isempty (BiasIn)))
              error (strcat ("ClassificationLinear: 'Bias' must be a", ...
                             " real numeric vector."));
            endif

          case 'fitbias'
            FitBias = varargin{2};
            if (! (islogical (FitBias) || (isnumeric (FitBias)
                   && isscalar (FitBias) && any (FitBias == [0, 1]))))
              error (strcat ("ClassificationLinear: 'FitBias' must be", ...
                             " either true or false."));
            endif
            FitBias = logical (FitBias);

          case 'postfitbias'
            PostFitBias = varargin{2};
            if (! (islogical (PostFitBias) || (isnumeric (PostFitBias)
                   && isscalar (PostFitBias) && any (PostFitBias == [0, 1]))))
              error (strcat ("ClassificationLinear: 'PostFitBias' must", ...
                             " be either true or false."));
            endif
            PostFitBias = logical (PostFitBias);

          case 'observationsin'
            ObservationsIn = varargin{2};
            if (! (ischar (ObservationsIn)
                   && any (strcmpi (ObservationsIn, {'rows', 'columns'}))))
              error (strcat ("ClassificationLinear: 'ObservationsIn'", ...
                             " must be either 'rows' or 'columns'."));
            endif
            ObservationsIn = lower (ObservationsIn);

          case 'betatolerance'
            BetaTolerance = varargin{2};
            if (! (isnumeric (BetaTolerance) && isscalar (BetaTolerance)
                   && isreal (BetaTolerance) && BetaTolerance >= 0))
              error (strcat ("ClassificationLinear: 'BetaTolerance'", ...
                             " must be a nonnegative scalar."));
            endif

          case 'gradienttolerance'
            GradientTolerance = varargin{2};
            if (! (isnumeric (GradientTolerance)
                   && isscalar (GradientTolerance)
                   && isreal (GradientTolerance) && GradientTolerance >= 0))
              error (strcat ("ClassificationLinear:", ...
                             " 'GradientTolerance' must be a", ...
                             " nonnegative scalar."));
            endif

          case 'deltagradienttolerance'
            DeltaGradientTolerance = varargin{2};
            if (! (isnumeric (DeltaGradientTolerance)
                   && isscalar (DeltaGradientTolerance)
                   && isreal (DeltaGradientTolerance)
                   && DeltaGradientTolerance >= 0))
              error (strcat ("ClassificationLinear:", ...
                             " 'DeltaGradientTolerance' must be a", ...
                             " nonnegative scalar."));
            endif

          case 'iterationlimit'
            IterationLimit = varargin{2};
            if (! (isnumeric (IterationLimit) && isscalar (IterationLimit)
                   && isreal (IterationLimit) && IterationLimit > 0
                   && fix (IterationLimit) == IterationLimit))
              error (strcat ("ClassificationLinear: 'IterationLimit'", ...
                             " must be a positive integer scalar."));
            endif

          case 'passlimit'
            PassLimit = varargin{2};
            if (! (isnumeric (PassLimit) && isscalar (PassLimit)
                   && isreal (PassLimit) && PassLimit > 0
                   && fix (PassLimit) == PassLimit))
              error (strcat ("ClassificationLinear: 'PassLimit' must be", ...
                             " a positive integer scalar."));
            endif

          case 'batchsize'
            BatchSize = varargin{2};
            if (! (isnumeric (BatchSize) && isscalar (BatchSize)
                   && isreal (BatchSize) && BatchSize > 0
                   && fix (BatchSize) == BatchSize))
              error (strcat ("ClassificationLinear: 'BatchSize' must be", ...
                             " a positive integer scalar."));
            endif

          case 'batchlimit'
            BatchLimit = varargin{2};
            if (! (isnumeric (BatchLimit) && isscalar (BatchLimit)
                   && isreal (BatchLimit) && BatchLimit > 0
                   && fix (BatchLimit) == BatchLimit))
              error (strcat ("ClassificationLinear: 'BatchLimit' must", ...
                             " be a positive integer scalar."));
            endif

          case 'learnrate'
            LearnRate = varargin{2};
            if (! (isnumeric (LearnRate) && isscalar (LearnRate)
                   && isreal (LearnRate) && LearnRate > 0))
              error (strcat ("ClassificationLinear: 'LearnRate' must be", ...
                             " a positive scalar."));
            endif

          case 'optimizelearnrate'
            OptimizeLearnRate = varargin{2};
            if (! (islogical (OptimizeLearnRate)
                   || (isnumeric (OptimizeLearnRate)
                       && isscalar (OptimizeLearnRate)
                       && any (OptimizeLearnRate == [0, 1]))))
              error (strcat ("ClassificationLinear:", ...
                             " 'OptimizeLearnRate' must be either true", ...
                             " or false."));
            endif
            OptimizeLearnRate = logical (OptimizeLearnRate);

          case 'truncationperiod'
            TruncationPeriod = varargin{2};
            if (! (isnumeric (TruncationPeriod)
                   && isscalar (TruncationPeriod)
                   && isreal (TruncationPeriod) && TruncationPeriod > 0
                   && fix (TruncationPeriod) == TruncationPeriod))
              error (strcat ("ClassificationLinear:", ...
                             " 'TruncationPeriod' must be a positive", ...
                             " integer scalar."));
            endif

          case 'numcheckconvergence'
            NumCheckConvergence = varargin{2};
            if (! (isnumeric (NumCheckConvergence)
                   && isscalar (NumCheckConvergence)
                   && isreal (NumCheckConvergence)
                   && NumCheckConvergence > 0
                   && fix (NumCheckConvergence) == NumCheckConvergence))
              error (strcat ("ClassificationLinear:", ...
                             " 'NumCheckConvergence' must be a positive", ...
                             " integer scalar."));
            endif

          case 'hessianhistorysize'
            HessianHistorySize = varargin{2};
            if (! (isnumeric (HessianHistorySize)
                   && isscalar (HessianHistorySize)
                   && isreal (HessianHistorySize) && HessianHistorySize > 0
                   && fix (HessianHistorySize) == HessianHistorySize))
              error (strcat ("ClassificationLinear:", ...
                             " 'HessianHistorySize' must be a positive", ...
                             " integer scalar."));
            endif

          case 'verbose'
            Verbose = varargin{2};
            if (! (isnumeric (Verbose) && isscalar (Verbose)
                   && isreal (Verbose) && any (Verbose == [0, 1, 2])))
              error (strcat ("ClassificationLinear: 'Verbose' must be", ...
                             " 0, 1, or 2."));
            endif

          case 'classnames'
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames)
                   || islogical (ClassNames) || ischar (ClassNames)))
              error (strcat ("ClassificationLinear: 'ClassNames' must", ...
                             " be a cell array of character vectors, a", ...
                             " logical vector, a numeric vector, or a", ...
                             " character array."));
            endif

          case 'cost'
            CostIn = varargin{2};
            if (! (isnumeric (CostIn) && isreal (CostIn)
                   && ismatrix (CostIn) && ndims (CostIn) == 2
                   && rows (CostIn) == columns (CostIn)))
              error (strcat ("ClassificationLinear: 'Cost' must be a", ...
                             " square numeric matrix."));
            endif

          case 'prior'
            Prior = varargin{2};
            if (! ((ischar (Prior) && any (strcmpi (Prior, {'empirical', ...
                                                            'uniform'})))
                   || (isnumeric (Prior) && isreal (Prior)
                       && isvector (Prior) && all (Prior >= 0))
                   || (isstruct (Prior) && isscalar (Prior))))
              error (strcat ("ClassificationLinear: 'Prior' must be", ...
                             " 'empirical', 'uniform', a vector of", ...
                             " nonnegative values, or a structure."));
            endif

          case 'scoretransform'
            ScoreTransform = varargin{2};

          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isreal (Weights)
                   && isvector (Weights) && all (Weights >= 0)))
              error (strcat ("ClassificationLinear: 'Weights' must be a", ...
                             " vector of nonnegative values."));
            endif

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! (iscellstr (PredictorNames)
                   && isvector (PredictorNames)))
              error (strcat ("ClassificationLinear: 'PredictorNames'", ...
                             " must be a cell array of character", ...
                             " vectors."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! (ischar (ResponseName) && isrow (ResponseName)))
              error (strcat ("ClassificationLinear: 'ResponseName' must", ...
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
              error (strcat ("ClassificationLinear:", ...
                             " 'CategoricalPredictors' must be a vector", ...
                             " of positive integers or a logical", ...
                             " vector."));
            endif

          otherwise
            error (strcat ("ClassificationLinear: invalid parameter", ...
                           " name in optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Observations may be given down the columns, which only means the
      ## predictor matrix arrives transposed.
      if (strcmp (ObservationsIn, 'columns'))
        X = X';
      endif

      ## Validate the data, resolve the classes, the prior, the cost and
      ## the observation weights.  The four linear and kernel classifiers
      ## share that opening, so it lives in one place.
      F = classFrame (X, Y, ClassNames, Prior, CostIn, Weights, ...
                      'ClassificationLinear');
      X = F.X;
      y = F.y;
      W = F.W;
      n = F.n;
      p = F.p;
      this.ClassNames = F.ClassNames;
      this.Prior = F.Prior;
      this.Cost = F.Cost;

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
          error (strcat ("ClassificationLinear: the '%s' solver fits a", ...
                         " ridge penalty only."), Solver{k});
        endif
        if (strcmp (Regularization, 'ridge') && strcmp (Solver{k}, 'sparsa'))
          error (strcat ("ClassificationLinear: the 'sparsa' solver", ...
                         " fits a lasso penalty only."));
        endif
        if (strcmp (Solver{k}, 'dual') && strcmp (Learner, 'logistic'))
          error (strcat ("ClassificationLinear: the 'dual' solver fits", ...
                         " a hinge loss only, so it needs 'Learner' set", ...
                         " to 'svm'."));
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
          error (strcat ("ClassificationLinear: 'Beta' must have one", ...
                         " row per predictor."));
        endif
        if (columns (BetaIn) == 1)
          Beta0 = repmat (BetaIn, 1, L);
        elseif (columns (BetaIn) == L)
          Beta0 = BetaIn;
        else
          error (strcat ("ClassificationLinear: 'Beta' must have one", ...
                         " column, or one per value of 'Lambda'."));
        endif
      endif
      if (isempty (BiasIn))
        if (strcmp (Learner, 'logistic'))
          Bias0 = repmat (sum (W .* y), 1, L);
        else
          Bias0 = zeros (1, L);
        endif
      else
        BiasIn = BiasIn(:)';
        if (numel (BiasIn) == 1)
          Bias0 = repmat (BiasIn, 1, L);
        elseif (numel (BiasIn) == L)
          Bias0 = BiasIn;
        else
          error (strcat ("ClassificationLinear: 'Bias' must be a", ...
                         " scalar, or hold one value per value of", ...
                         " 'Lambda'."));
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
      P.LossFunction = 'hinge';
      if (strcmp (Learner, 'logistic'))
        P.LossFunction = 'logit';
      endif
      P.Epsilon = [];
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
          [b, b0, S] = linearSolve (X, y, W, P);
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
      this.PredictorNames = PredictorNames;
      if (isempty (this.PredictorNames))
        this.PredictorNames = ...
                     arrayfun (@(k) sprintf ("x%d", k), 1:p, ...
                               'UniformOutput', false);
      elseif (numel (this.PredictorNames) != p)
        error (strcat ("ClassificationLinear: 'PredictorNames' must", ...
                       " have one name per predictor."));
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
      if (strcmp (Regularization, 'ridge'))
        this.Regularization = 'ridge (L2)';
      else
        this.Regularization = 'lasso (L1)';
      endif

      if (isempty (ScoreTransform))
        if (strcmp (Learner, 'logistic'))
          ScoreTransform = 'logit';
        else
          ScoreTransform = 'none';
        endif
      endif
      this.ScoreTransform = ScoreTransform;

      this.ModelParameters = linearModelParams (P, Solver, LambdaIn, ...
                                                [], Beta0, Bias0, ...
                                                Verbose, 'classification');

      this.FitInfo_ = linearFitInfo (info, BetaTolerance, ...
                                     GradientTolerance, ...
                                     DeltaGradientTolerance, ...
                                     IterationLimit, Solver, PassLimit, ...
                                     BatchLimit);

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationLinear} {@var{labels} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationLinear} {[@var{labels}, @var{scores}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new observations.
    ##
    ## @code{@var{labels} = predict (@var{obj}, @var{XC})} returns the class
    ## of largest score for each row of @var{XC}, in the type of the
    ## response the model was fitted to.  With @math{L} regularization
    ## strengths @var{labels} has one column per strength.
    ##
    ## @code{[@var{labels}, @var{scores}] = predict (@var{obj}, @var{XC})}
    ## also returns the scores, an @math{Nx2} matrix whose columns follow
    ## @qcode{ClassNames}, or an @math{Nx2xL} array with more than one
    ## strength.  The scores are @math{-f} and @math{+f} for the raw model
    ## value @math{f}, after @qcode{ScoreTransform} has been applied.
    ##
    ## @end deftypefn
    function [labels, scores] = predict (this, XC)

      if (nargin < 2)
        error ("ClassificationLinear.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("ClassificationLinear.predict: XC is empty.");
      endif
      if (! (isnumeric (XC) && isreal (XC) && ismatrix (XC)))
        error ("ClassificationLinear.predict: invalid values in XC.");
      endif
      if (columns (XC) != this.NumPredictors_)
        error (strcat ("ClassificationLinear.predict: XC must have the", ...
                       " same number of predictors as the trained model."));
      endif

      f = XC * this.Beta + this.Bias;
      ## A character matrix holds one name per row and cannot hold a grid of
      ## them, so a response given that way has nowhere to put the labels of
      ## more than one regularization strength.  Refused rather than
      ## returned in a shape nothing can read.
      if (ischar (this.ClassNames) && numel (this.Lambda) > 1)
        error (strcat ("ClassificationLinear.predict: a character matrix", ...
                       " response cannot carry labels for more than one", ...
                       " value of 'Lambda'; give the response as a cell", ...
                       " array of character vectors."));
      endif
      labels = labelsFromIndex (this.ClassNames, 1 + (f > 0));

      if (nargout > 1)
        L = numel (this.Lambda);
        scores = zeros (rows (XC), 2, L);
        for l = 1:L
          scores(:,:,l) = this.STfun ([-f(:,l), f(:,l)]);
        endfor
        if (L == 1)
          scores = scores(:,:,1);
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationLinear} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin of each observation.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns the
    ## score of the true class less the score of the other one, one row per
    ## observation and one column per regularization strength.  A positive
    ## margin is a correct classification.
    ##
    ## @end deftypefn
    function m = margin (this, X, Y)

      if (nargin < 3)
        error ("ClassificationLinear.margin: too few input arguments.");
      endif
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("ClassificationLinear.margin: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("ClassificationLinear.margin: number of rows in", ...
                       " X and Y must be equal."));
      endif

      [~, scores] = predict (this, X);
      L = numel (this.Lambda);
      m = zeros (rows (X), L);
      for l = 1:L
        s = scores(:,:,l);
        strue = s(sub2ind (size (s), (1:rows (s))', gY));
        sother = s(sub2ind (size (s), (1:rows (s))', 3 - gY));
        m(:,l) = strue - sother;
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationLinear} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationLinear} {@var{e} =} edge (@dots{}, @qcode{'Weights'}, @var{W})
    ##
    ## Weighted mean of the classification margins.
    ##
    ## @code{@var{e} = edge (@var{obj}, @var{X}, @var{Y})} returns one value
    ## per regularization strength.  The weights are normalized within each
    ## class to that class's prior before they are applied.
    ##
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      if (nargin < 3)
        error ("ClassificationLinear.edge: too few input arguments.");
      endif
      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       'ClassificationLinear', 'edge');
      m = margin (this, X, Y);
      e = sum (W .* m, 1);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationLinear} {@var{l} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationLinear} {@var{l} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss on new data.
    ##
    ## @code{@var{l} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## misclassification rate, one value per regularization strength.
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
        error ("ClassificationLinear.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationLinear.loss: optional arguments", ...
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
              error (strcat ("ClassificationLinear.loss: 'LossFun' must", ...
                             " be 'binodeviance', 'classifcost',", ...
                             " 'classiferror', 'exponential', 'hinge',", ...
                             " 'logit', 'mincost', or 'quadratic'."));
            endif
            LossFun = lower (LossFun);
          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isreal (Weights)
                   && isvector (Weights) && all (Weights >= 0)))
              error (strcat ("ClassificationLinear.loss: 'Weights' must", ...
                             " be a vector of nonnegative values."));
            endif
          otherwise
            error (strcat ("ClassificationLinear.loss: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("ClassificationLinear.loss: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("ClassificationLinear.loss: number of rows in X", ...
                       " and Y must be equal."));
      endif
      if (isempty (Weights))
        w = ones (numel (gY), 1);
      else
        w = Weights(:);
        if (numel (w) != numel (gY))
          error (strcat ("ClassificationLinear.loss: 'Weights' must", ...
                         " have one element per observation."));
        endif
      endif
      w = w / sum (w);

      [~, scores] = predict (this, X);
      L = numel (this.Lambda);
      l = zeros (1, L);
      for k = 1:L
        s = scores(:,:,k);
        l(k) = classificationLoss (LossFun, s, gY, w, this.Cost);
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationLinear} {@var{sub} =} selectModels (@var{obj}, @var{idx})
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
        error (strcat ("ClassificationLinear.selectModels: too few", ...
                       " input arguments."));
      endif
      L = numel (this.Lambda);
      if (islogical (idx))
        if (numel (idx) != L)
          error (strcat ("ClassificationLinear.selectModels: a logical", ...
                         " IDX must have one element per value of", ...
                         " 'Lambda'."));
        endif
        idx = find (idx);
      endif
      if (! (isnumeric (idx) && isreal (idx) && isvector (idx)
             && ! isempty (idx) && all (fix (idx) == idx)
             && all (idx >= 1) && all (idx <= L)))
        error (strcat ("ClassificationLinear.selectModels: IDX must", ...
                       " hold integers between 1 and %d."), L);
      endif

      sub = this;
      sub.Lambda = this.Lambda(idx);
      sub.Beta = this.Beta(:,idx);
      sub.Bias = this.Bias(idx);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationLinear} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a linear classifier to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves the model
    ## @var{obj} into @var{filename} in a form @code{loadmodel} can read
    ## back.
    ##
    ## @end deftypefn
    function savemodel (obj, fname)

      classdef_name = 'ClassificationLinear';
      ClassNames = obj.ClassNames;
      Prior = obj.Prior;
      Cost = obj.Cost;
      ScoreTransform = obj.ScoreTransform;
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

      save ('-binary', fname, 'classdef_name', 'ClassNames', 'Prior', ...
            'Cost', 'ScoreTransform', 'PredictorNames', ...
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
      printf ("\n  ClassificationLinear\n\n");
      printf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      printf ("%+25s: %s\n", 'ClassNames', classNameListing (this.ClassNames));
      printf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
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

    ## What the fit reported, which fitclinear returns as its second output.
    function S = fitInfo_ (this)
      S = this.FitInfo_;
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      mdl = ClassificationLinear (zeros (2, 1), [0; 1]);
      fields = fieldnames (data);
      for k = 1:numel (fields)
        mdl.(fields{k}) = data.(fields{k});
      endfor

    endfunction

  endmethods

endclassdef

%!demo
%! ## Separate the two overlapping iris species with a linear classifier and
%! ## read the posterior probability it gives each observation.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationLinear (X, Y, 'Learner', 'logistic')
%! [label, score] = predict (Mdl, X([1, 51],:))

%!demo
%! ## One object can hold a whole regularization path.  A stronger penalty
%! ## shrinks the coefficients, and every method reports one column per
%! ## strength.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationLinear (X, Y, 'Lambda', [0.001, 0.01, 0.1]);
%! Mdl.Beta
%! loss (Mdl, X, Y)

%!test
%! ## The model reports the surface MATLAB reports
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationLinear (X, Y);
%! assert_equal (class (Mdl), 'ClassificationLinear');
%! assert_equal (Mdl.ClassNames, {'versicolor'; 'virginica'});
%! assert_equal (Mdl.Learner, 'svm');
%! assert_equal (Mdl.FittedLoss, 'hinge');
%! assert_equal (Mdl.Regularization, 'ridge (L2)');
%! assert_equal (Mdl.ScoreTransform, 'none');
%! assert_equal (Mdl.ResponseName, 'Y');
%! assert_equal (Mdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (Mdl.ExpandedPredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (Mdl.Prior, [0.5, 0.5]);
%! assert_equal (Mdl.Cost, [0, 1; 1, 0]);
%! assert_equal (Mdl.Lambda, 0.01);
%! assert_equal (size (Mdl.Beta), [4, 1]);

%!test
%! ## The properties are the ones MATLAB lists, in its order
%! load fisheriris
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end));
%! assert_equal (sort (properties (Mdl)), ...
%!               sort ({'ClassNames'; 'Prior'; 'Cost'; 'ScoreTransform'; ...
%!                      'PredictorNames'; 'CategoricalPredictors'; ...
%!                      'ResponseName'; 'ExpandedPredictorNames'; ...
%!                      'Learner'; 'Beta'; 'Bias'; 'FittedLoss'; 'Lambda'; ...
%!                      'ModelParameters'; 'Regularization'}));

%!test
%! ## A logistic ridge fit reproduces R2024a's coefficients
%! load fisheriris
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end), ...
%!                             'Learner', 'logistic', 'Solver', 'lbfgs', ...
%!                             'BetaTolerance', 0, 'GradientTolerance', ...
%!                                                 1e-12, ...
%!                             'IterationLimit', 20000);
%! assert_equal (Mdl.Beta, [-0.394433478724131; -0.513277404049627; ...
%!                          2.9307513837558; 2.41703218835114], 1e-7);
%! assert_equal (Mdl.Bias, -14.4307581801051, 1e-7);
%! assert_equal (Mdl.FittedLoss, 'logit');
%! assert_equal (Mdl.ScoreTransform, 'logit');

%!test
%! ## The default score transform follows the learner, and only the logistic
%! ## one turns the scores into posterior probabilities
%! load fisheriris
%! X = meas(51:end,:);
%! Ml = ClassificationLinear (X, species(51:end), 'Learner', 'logistic');
%! Ms = ClassificationLinear (X, species(51:end), 'Learner', 'svm');
%! assert_equal (Ml.ScoreTransform, 'logit');
%! assert_equal (Ms.ScoreTransform, 'none');
%! [~, sl] = predict (Ml, X(1:5,:));
%! [~, ss] = predict (Ms, X(1:5,:));
%! assert_equal (sum (sl, 2), ones (5, 1), 1e-12);
%! assert_equal (ss(:,1), -ss(:,2), 1e-12);

%!test
%! ## The scores of a logistic fit are R2024a's posteriors
%! load fisheriris
%! X = meas(51:end,:);
%! Mdl = ClassificationLinear (X, species(51:end), 'Learner', 'logistic', ...
%!                             'Solver', 'lbfgs', 'BetaTolerance', 0, ...
%!                             'GradientTolerance', 1e-12, ...
%!                             'IterationLimit', 20000);
%! [label, score] = predict (Mdl, X(1:2,:));
%! assert_equal (label, {'versicolor'; 'versicolor'});
%! assert_equal (score, [0.842361345526662, 0.157638654473338; ...
%!                       0.856151985655854, 0.143848014344146], 1e-7);

%!test
%! ## margin is the true class score less the other, and edge their weighted
%! ## mean.  Both halves are asserted: the values against R2024a, and the
%! ## relationship against the scores this model actually reports.  The
%! ## values alone would not catch a margin computed consistently wrongly,
%! ## since a wrong-but-stable margin agrees with its own past output
%! ## forever; the relationship alone would not catch our agreeing with
%! ## ourselves about the wrong quantity.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationLinear (X, Y, 'Learner', 'logistic', ...
%!                             'Solver', 'lbfgs', 'BetaTolerance', 0, ...
%!                             'GradientTolerance', 1e-12, ...
%!                             'IterationLimit', 20000);
%! assert_equal (margin (Mdl, X(1,:), Y(1)), 0.684722691053324, 1e-7);
%! assert_equal (edge (Mdl, X, Y), 0.726067240395032, 1e-7);
%! ## the invariant, independent of any solver output
%! [~, score] = predict (Mdl, X);
%! virg = strcmp (Y, 'virginica');
%! strue = score(:,1);
%! strue(virg) = score(virg,2);
%! sother = score(:,2);
%! sother(virg) = score(virg,1);
%! assert_equal (margin (Mdl, X, Y), strue - sother, 1e-15);
%! assert_equal (edge (Mdl, X, Y), mean (margin (Mdl, X, Y)), 1e-15);

%!test
%! ## Every loss reproduces R2024a, and each is a function of the score the
%! ## model gives the true class rather than of the margin
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationLinear (X, Y, 'Learner', 'logistic', ...
%!                             'Solver', 'lbfgs', 'BetaTolerance', 0, ...
%!                             'GradientTolerance', 1e-12, ...
%!                             'IterationLimit', 20000);
%! assert_equal (loss (Mdl, X, Y), 0.04, 1e-12);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'hinge'), ...
%!               0.136966379802484, 1e-7);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'logit'), ...
%!               0.354333320227268, 1e-7);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'binodeviance'), ...
%!               0.170050306584128, 1e-7);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'exponential'), ...
%!               0.426899191073676, 1e-7);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'quadratic'), ...
%!               0.040701761220091, 1e-7);

%!test
%! ## Under the default cost both cost based losses reduce to the error
%! ## rate: the class of least expected cost is then the class of largest
%! ## score, which is what predict returns
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationLinear (X, Y, 'Learner', 'logistic');
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'classifcost'), ...
%!               loss (Mdl, X, Y), 1e-12);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'mincost'), ...
%!               loss (Mdl, X, Y), 1e-12);

%!test
%! ## A cost that is not symmetric moves the least cost assignment away from
%! ## the largest score, so the two cost based losses part company
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationLinear (X, Y, 'Learner', 'logistic', ...
%!                             'Cost', [0, 2; 3, 0]);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'classifcost') > 0, true);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'mincost') > 0, true);

%!test
%! ## A cost matrix reaches the fit through the prior, which is what MATLAB
%! ## does: a fit costing four times as much to miss the first class equals
%! ## one given the prior that scaling implies
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mc = ClassificationLinear (X, Y, 'Cost', [0, 4; 1, 0]);
%! Mp = ClassificationLinear (X, Y, 'Prior', [0.8, 0.2]);
%! assert_equal (Mc.Beta, Mp.Beta, 1e-12);
%! assert_equal (Mc.Bias, Mp.Bias, 1e-12);
%! assert_equal (Mc.Prior, [0.5, 0.5]);

%!test
%! ## Lambda defaults to the reciprocal of the observations that were used
%! load fisheriris
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end));
%! assert_equal (Mdl.Lambda, 1 / 100, 1e-15);

%!test
%! ## A vector of strengths fits one model per value, sorted ascending, and
%! ## every method reports one column per value
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationLinear (X, Y, 'Lambda', [0.1, 0.001, 0.01]);
%! assert_equal (Mdl.Lambda, [0.001, 0.01, 0.1]);
%! assert_equal (size (Mdl.Beta), [4, 3]);
%! assert_equal (size (Mdl.Bias), [1, 3]);
%! [label, score] = predict (Mdl, X(1:4,:));
%! assert_equal (size (label), [4, 3]);
%! assert_equal (size (score), [4, 2, 3]);
%! assert_equal (size (margin (Mdl, X, Y)), [100, 3]);
%! assert_equal (size (edge (Mdl, X, Y)), [1, 3]);
%! assert_equal (size (loss (Mdl, X, Y)), [1, 3]);

%!test
%! ## A stronger penalty shrinks the coefficients
%! load fisheriris
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end), ...
%!                             'Learner', 'logistic', 'Lambda', [0.001, 1]);
%! assert_equal (norm (Mdl.Beta(:,1)) > norm (Mdl.Beta(:,2)), true);

%!test
%! ## selectModels keeps the strengths it is given and drops the rest
%! load fisheriris
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end), ...
%!                             'Lambda', [0.001, 0.01, 0.1]);
%! sub = selectModels (Mdl, [1, 3]);
%! assert_equal (sub.Lambda, [0.001, 0.1]);
%! assert_equal (size (sub.Beta), [4, 2]);
%! assert_equal (sub.Beta, Mdl.Beta(:,[1, 3]));
%! sub = selectModels (Mdl, logical ([0, 1, 0]));
%! assert_equal (sub.Lambda, 0.01);

%!test
%! ## A lasso penalty drives coefficients to exactly zero, which a ridge
%! ## penalty never does
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Ml = ClassificationLinear (X, Y, 'Learner', 'logistic', ...
%!                            'Regularization', 'lasso', 'Lambda', 0.05);
%! Mr = ClassificationLinear (X, Y, 'Learner', 'logistic', ...
%!                            'Regularization', 'ridge', 'Lambda', 0.05);
%! assert_equal (Ml.Regularization, 'lasso (L1)');
%! assert_equal (sum (Ml.Beta == 0) > 0, true);
%! assert_equal (any (Mr.Beta == 0), false);

%!test
%! ## The default solver follows the penalty and the width of the data
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! [~, F] = fitclinear (X, Y);
%! assert_equal (F.Solver, {'bfgs'});
%! [~, F] = fitclinear (X, Y, 'Regularization', 'lasso');
%! assert_equal (F.Solver, {'sparsa'});

%!test
%! ## FitBias false leaves the intercept at zero
%! load fisheriris
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end), ...
%!                             'FitBias', false);
%! assert_equal (Mdl.Bias, 0);

%!test
%! ## Observations may be given down the columns instead
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mr = ClassificationLinear (X, Y, 'Learner', 'logistic');
%! Mc = ClassificationLinear (X', Y, 'Learner', 'logistic', ...
%!                            'ObservationsIn', 'columns');
%! assert_equal (Mr.Beta, Mc.Beta, 1e-12);

%!test
%! ## A row with a missing predictor or a missing response is dropped, and
%! ## Lambda follows the count that survived
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! X(3,2) = NaN;
%! Mdl = ClassificationLinear (X, Y);
%! assert_equal (Mdl.Lambda, 1 / 99, 1e-15);

%!test
%! ## ScoreTransform is settable and reaches predict
%! load fisheriris
%! X = meas(51:end,:);
%! Mdl = ClassificationLinear (X, species(51:end));
%! Mdl.ScoreTransform = 'logit';
%! assert_equal (Mdl.ScoreTransform, 'logit');
%! [~, score] = predict (Mdl, X(1:3,:));
%! assert_equal (sum (score, 2), ones (3, 1), 1e-12);

%!test
%! ## A saved model reads back as the same model
%! load fisheriris
%! X = meas(51:end,:);
%! Mdl = ClassificationLinear (X, species(51:end), 'Learner', 'logistic');
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mnew = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (Mnew), 'ClassificationLinear');
%! assert_equal (Mnew.Beta, Mdl.Beta);
%! assert_equal (Mnew.Bias, Mdl.Bias);
%! assert_equal (predict (Mnew, X(1:5,:)), predict (Mdl, X(1:5,:)));

%!test
%! ## The labels come back in the type the response was given in, a
%! ## character matrix included: its rows name the classes, and indexing it
%! ## by element rather than by row would flatten it into single characters
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mcell = ClassificationLinear (X, Y);
%! assert_equal (class (predict (Mcell, X(1:3,:))), 'cell');
%! Mchar = ClassificationLinear (X, char (Y));
%! assert_equal (size (Mchar.ClassNames), [2, 10]);
%! assert_equal (size (predict (Mchar, X(1:3,:))), [3, 10]);
%! assert_equal (loss (Mchar, X, char (Y)), loss (Mcell, X, Y), 1e-12);
%! Ynum = double (strcmp (Y, 'virginica'));
%! Mnum = ClassificationLinear (X, Ynum);
%! assert_equal (class (predict (Mnum, X(1:3,:))), 'double');
%! assert_equal (Mnum.ClassNames, [0; 1]);
%! Mlog = ClassificationLinear (X, logical (Ynum));
%! assert_equal (class (predict (Mlog, X(1:3,:))), 'logical');

%!test
%! ## Observation weights reach the fit, and the prior follows them.  Both
%! ## numbers are R2024a's for the weights 1 to 100.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! Mdl = ClassificationLinear (X, Y, 'Learner', 'logistic', ...
%!                             'Solver', 'lbfgs', 'Weights', (1:100)', ...
%!                             'BetaTolerance', 0, ...
%!                             'GradientTolerance', 1e-12, ...
%!                             'IterationLimit', 20000);
%! assert_equal (Mdl.Prior, [0.252475247524752, 0.747524752475248], 1e-12);
%! assert_equal (Mdl.Beta, [-0.0907861095243143; -0.344147850956068; ...
%!                          2.73449884653132; 2.20876095012392], 1e-7);
%! assert_equal (Mdl.Bias, -14.6291356484861, 1e-6);
%! assert_equal (Mdl.FitInfo_.Objective, 0.208158687811647, 1e-10);

%!test
%! ## A character matrix response names one class per row, and a model
%! ## fitted from one is the model fitted from the equivalent cell array.
%! ## The check is worth its length: a linear index into a character matrix
%! ## flattens the names into single characters, and so does an ismember
%! ## between two character matrices, and neither shows up at the property
%! ## that is usually looked at.
%! load fisheriris
%! X = meas(51:end,:);
%! Yc = species(51:end);
%! Ym = char (Yc);
%! Mc = ClassificationLinear (X, Yc, 'Learner', 'logistic');
%! Mm = ClassificationLinear (X, Ym, 'Learner', 'logistic');
%! assert_equal (cellstr (Mm.ClassNames), Mc.ClassNames);
%! assert_equal (Mm.Prior, Mc.Prior);
%! assert_equal (Mm.Beta, Mc.Beta);
%! assert_equal (Mm.Bias, Mc.Bias);

%!test
%! ## Every method answers the same through a character matrix response as
%! ## through the cell array it came from
%! load fisheriris
%! X = meas(51:end,:);
%! Yc = species(51:end);
%! Ym = char (Yc);
%! Mc = ClassificationLinear (X, Yc, 'Learner', 'logistic');
%! Mm = ClassificationLinear (X, Ym, 'Learner', 'logistic');
%! assert_equal (cellstr (predict (Mm, X)), predict (Mc, X));
%! assert_equal (margin (Mm, X, Ym), margin (Mc, X, Yc));
%! assert_equal (edge (Mm, X, Ym), edge (Mc, X, Yc));
%! assert_equal (edge (Mm, X, Ym, 'Weights', (1:100)'), ...
%!               edge (Mc, X, Yc, 'Weights', (1:100)'));
%! assert_equal (loss (Mm, X, Ym), loss (Mc, X, Yc));
%! assert_equal (loss (Mm, X, Ym, 'LossFun', 'hinge'), ...
%!               loss (Mc, X, Yc, 'LossFun', 'hinge'));
%! assert_equal (loss (Mm, X, Ym, 'LossFun', 'mincost'), ...
%!               loss (Mc, X, Yc, 'LossFun', 'mincost'));

%!test
%! ## 'ClassNames' selects a subset when it is given as a character matrix
%! ## too, which needs both sides of the comparison turned into names
%! load fisheriris
%! Mcell = ClassificationLinear (meas, species, ...
%!                               'ClassNames', {'versicolor', 'virginica'});
%! Mchar = ClassificationLinear (meas, char (species), 'ClassNames', ...
%!                               char ({'versicolor', 'virginica'}));
%! assert_equal (cellstr (Mchar.ClassNames), Mcell.ClassNames);
%! assert_equal (Mchar.Beta, Mcell.Beta);
%! assert_equal (Mchar.Lambda, Mcell.Lambda);

%!test
%! ## Names of unequal length are padded by the character matrix and the
%! ## padding is not part of the name
%! load fisheriris
%! X = meas(51:end,:);
%! Y = [repmat({'ab'}, 50, 1); repmat({'abcd'}, 50, 1)];
%! Mcell = ClassificationLinear (X, Y);
%! Mchar = ClassificationLinear (X, char (Y));
%! assert_equal (cellstr (Mchar.ClassNames), Mcell.ClassNames);
%! assert_equal (cellstr (predict (Mchar, X)), predict (Mcell, X));

%!test
%! ## The default fit stops on the coefficients, as MATLAB's does, which is
%! ## only true once the engine offers BetaTolerance.  Before it did, this
%! ## fit ran on to the gradient tolerance and reported NaN here.
%! load fisheriris
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end), ...
%!                             'Learner', 'logistic');
%! assert_equal (Mdl.FitInfo_.BetaTolerance, 1e-4);
%! assert_equal (Mdl.FitInfo_.TerminationCode, 1);
%! assert_equal (Mdl.FitInfo_.TerminationStatus, ...
%!               {'Tolerance on coefficients satisfied.'});
%! assert_equal (isfinite (Mdl.FitInfo_.RelativeChangeInBeta), true);

%!test
%! ## A line search that cannot improve the objective says so, rather than
%! ## reporting the iteration limit it never reached.  Code and wording
%! ## measured on R2024a.  Where the search gives up varies by platform, so
%! ## the count is only asserted to be short of the limit.
%! load fisheriris
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end), ...
%!                             'Learner', 'logistic', 'Lambda', 1e-10, ...
%!                             'BetaTolerance', 0, 'GradientTolerance', 0);
%! assert_equal (Mdl.FitInfo_.TerminationCode, -11);
%! assert_equal (Mdl.FitInfo_.TerminationStatus, ...
%!               {'Unable to find a step decreasing the objective.'});
%! S = Mdl.fitInfo_ ();
%! assert_equal (S.NumIterations < S.IterationLimit, true);

%!test
%! ## A tolerance of exactly zero switches its test off and the quantity it
%! ## governs comes back NaN, which is what MATLAB reports.  Not "never
%! ## satisfied": not computed.
%! load fisheriris
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end), ...
%!                             'BetaTolerance', 0, 'GradientTolerance', 1e-8);
%! assert_equal (isnan (Mdl.FitInfo_.RelativeChangeInBeta), true);
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end), ...
%!                             'GradientTolerance', 0);
%! assert_equal (isnan (Mdl.FitInfo_.GradientNorm), true);
%! assert_equal (Mdl.FitInfo_.GradientTolerance, 0);

%!test
%! ## A character matrix response survives a round trip through savemodel
%! ## and loadmodel, names and all
%! load fisheriris
%! X = meas(51:end,:);
%! Ym = char (species(51:end));
%! Mdl = ClassificationLinear (X, Ym, 'Learner', 'logistic');
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mnew = loadmodel (fname);
%! delete (fname);
%! assert_equal (Mnew.ClassNames, Mdl.ClassNames);
%! assert_equal (size (Mnew.ClassNames), [2, 10]);
%! assert_equal (predict (Mnew, X(1:5,:)), predict (Mdl, X(1:5,:)));
%! assert_equal (loss (Mnew, X, Ym), loss (Mdl, X, Ym), 1e-12);

%!test
%! ## ModelParameters reports every option but gives a value only to the
%! ## ones the chosen solver can act on, which is MATLAB's behaviour and was
%! ## measured one fit per solver on R2024a.  Two of these are easy to get
%! ## wrong: HessianHistorySize is empty for 'bfgs' and 15 for 'lbfgs', a
%! ## full quasi-Newton method having no limited memory to size, and the two
%! ## solvers that run no gradient test report a tolerance of zero.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! b = ClassificationLinear (X, Y, 'Solver', 'bfgs').ModelParameters;
%! assert_equal (numel (fieldnames (b)), 28);
%! assert_equal (b.HessianHistorySize, []);
%! assert_equal (b.BatchSize, []);
%! assert_equal (b.PassLimit, []);
%! assert_equal (b.IterationLimit, 1000);
%! assert_equal (b.LineSearch, 'strongwolfe');
%! l = ClassificationLinear (X, Y, 'Solver', 'lbfgs').ModelParameters;
%! assert_equal (l.HessianHistorySize, 15);
%! g = ClassificationLinear (X, Y, 'Solver', 'sgd').ModelParameters;
%! assert_equal (g.BatchSize, 10);
%! assert_equal (g.IterationLimit, []);
%! assert_equal (g.GradientTolerance, 0);
%! assert_equal (g.LineSearch, []);
%! d = ClassificationLinear (X, Y, 'Solver', 'dual').ModelParameters;
%! assert_equal (isempty (d.DeltaGradientTolerance), false);
%! assert_equal (d.PassLimit, 10);
%! assert_equal (d.HessianHistorySize, []);
%! assert_equal (d.GradientTolerance, 0);

%!test
%! ## The default learning rate is MATLAB's own formula, the reciprocal root
%! ## of one plus the largest squared observation norm.  R2024a reports
%! ## 0.0896365435911655 on this fixture.
%! load fisheriris
%! X = meas(51:end,:);
%! Mdl = ClassificationLinear (X, species(51:end), 'Solver', 'sgd');
%! assert_equal (Mdl.ModelParameters.LearnRate, ...
%!               1 / sqrt (1 + max (sum (X .^ 2, 2))), 1e-15);
%! assert_equal (Mdl.ModelParameters.LearnRate, 0.0896365435911655, 1e-12);

%!test
%! ## The dual solver's two defaults, measured on R2024a and confirmed on
%! ## R2026a through MATLAB Online.  Both differ from what MathWorks
%! ## documents, and the first differs between the two learners: the
%! ## complementarity tolerance is 1 for a hinge loss and 0.1 for an
%! ## epsilon-insensitive one, which is the 0.1 the documentation quotes.
%! load fisheriris
%! Mdl = ClassificationLinear (meas(51:end,:), species(51:end), ...
%!                             'Solver', 'dual');
%! assert_equal (Mdl.ModelParameters.DeltaGradientTolerance, 1);
%! assert_equal (Mdl.ModelParameters.NumCheckConvergence, 2);
%! assert_equal (Mdl.ModelParameters.PassLimit, 10);

## Test input validation
%!error<ClassificationLinear: too few input arguments.> ...
%! ClassificationLinear (ones (5, 2))
%!error<ClassificationLinear: optional arguments must be given in Name-Value pairs.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Learner')
%!error<ClassificationLinear: 'Learner' must be either 'svm' or 'logistic'.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Learner', ...
%!                     'tree')
%!error<ClassificationLinear: 'Regularization' must be either 'ridge' or 'lasso'.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'Regularization', 'elastic')
%!error<ClassificationLinear: 'Lambda' must be 'auto' or a vector of nonnegative finite values.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Lambda', -1)
%!error<ClassificationLinear: 'Solver' must be one of 'sgd', 'asgd', 'dual', 'bfgs', 'lbfgs' and 'sparsa', or a cell array of them.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Solver', ...
%!                     'newton')
%!error<ClassificationLinear: 'FitBias' must be either true or false.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'FitBias', ...
%!                     'yes')
%!error<ClassificationLinear: 'ObservationsIn' must be either 'rows' or 'columns'.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'ObservationsIn', 'pages')
%!error<ClassificationLinear: 'BetaTolerance' must be a nonnegative scalar.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'BetaTolerance', -1)
%!error<ClassificationLinear: 'IterationLimit' must be a positive integer scalar.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'IterationLimit', 2.5)
%!error<ClassificationLinear: 'Verbose' must be 0, 1, or 2.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Verbose', 3)
%!error<ClassificationLinear: invalid parameter name in optional pair arguments.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Nonsense', 1)
%!error<ClassificationLinear: invalid values in X.> ...
%! ClassificationLinear ({1, 2; 3, 4}, [1; 2])
%!error<ClassificationLinear: X is empty.> ClassificationLinear ([], [])
%!error<ClassificationLinear: number of rows in X and Y must be equal.> ...
%! ClassificationLinear (ones (10, 2), [1; 2])
%!error<ClassificationLinear: Y must name exactly two classes, this being a binary model; use fitcecoc for more than two.> ...
%! ClassificationLinear (ones (9, 2), [1; 1; 1; 2; 2; 2; 3; 3; 3])
%!error<ClassificationLinear: 'Prior' must have one entry per class.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Prior', ...
%!                     [0.2, 0.3, 0.5])
%!error<ClassificationLinear: the number of rows and columns in 'Cost' must correspond to the classes in Y.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Cost', ...
%!                     ones (3))
%!error<ClassificationLinear: the 'sparsa' solver fits a lasso penalty only.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'Regularization', 'ridge', 'Solver', 'sparsa')
%!error<ClassificationLinear: the 'lbfgs' solver fits a ridge penalty only.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'Regularization', 'lasso', 'Solver', 'lbfgs')
%!error<ClassificationLinear: the 'dual' solver fits a hinge loss only, so it needs 'Learner' set to 'svm'.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Learner', ...
%!                     'logistic', 'Solver', 'dual')
%!error<ClassificationLinear: 'PredictorNames' must have one name per predictor.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'PredictorNames', {'a', 'b', 'c'})
%!error<ClassificationLinear: 'Beta' must have one row per predictor.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Beta', ...
%!                     ones (3, 1))
%!error<ClassificationLinear: 'Bias' must be a scalar, or hold one value per value of 'Lambda'.> ...
%! ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Bias', [1, 2])
%!error<ClassificationLinear.predict: too few input arguments.> ...
%! predict (ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)]))
%!error<ClassificationLinear.predict: XC is empty.> ...
%! predict (ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)]), [])
%!error<ClassificationLinear.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (3, 5))
%!error<ClassificationLinear.margin: too few input arguments.> ...
%! margin (ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (3, 2))
%!error<ClassificationLinear.loss: 'LossFun' must be 'binodeviance', 'classifcost', 'classiferror', 'exponential', 'hinge', 'logit', 'mincost', or 'quadratic'.> ...
%! loss (ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (10, 2), [ones(5,1); 2*ones(5,1)], 'LossFun', 'mse')
%!error<ClassificationLinear.loss: invalid parameter name in optional pair arguments.> ...
%! loss (ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (10, 2), [ones(5,1); 2*ones(5,1)], 'Nonsense', 1)
%!error<ClassificationLinear.loss: Y must hold only classes the model was trained on.> ...
%! loss (ClassificationLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)]), ...
%!                     ones (10, 2), 3*ones (10, 1))
%!error<ClassificationLinear.selectModels: IDX must hold integers between 1 and 3.> ...
%! selectModels (ClassificationLinear (ones (10, 2), ...
%!                     [ones(5,1); 2*ones(5,1)], 'Lambda', [0.1, 0.2, 0.3]), 5)
%!error<ClassificationLinear.selectModels: a logical IDX must have one element per value of 'Lambda'.> ...
%! selectModels (ClassificationLinear (ones (10, 2), ...
%!                     [ones(5,1); 2*ones(5,1)], 'Lambda', [0.1, 0.2, ...
%!                     0.3]), logical ([1, 0]))
%!error<ClassificationLinear.predict: a character matrix response cannot carry labels for more than one value of 'Lambda'; give the response as a cell array of character vectors.> ...
%! predict (ClassificationLinear (ones (10, 2), ...
%!                     ['a';'a';'a';'a';'a';'b';'b';'b';'b';'b'], 'Lambda', ...
%!                     [0.1, 0.2]), ones (3, 2))

## Every documented score transform reaches the scores that are reported, and
## none of them moves the label: a transform reshapes what is reported, not
## what is decided.
%!test
%! load fisheriris
%! Mdl = fitclinear (meas, strcmp (species, 'setosa'));
%! Mdl.ScoreTransform = 'none';
%! [label, raw] = predict (Mdl, meas([1, 60, 120],:));
%! T = {'identity', @(x) x; 'doublelogit', @(x) 1 ./ (1 + exp (-2 * x)); ...
%!      'invlogit', @(x) log (x ./ (1 - x)); ...
%!      'logit', @(x) 1 ./ (1 + exp (-x)); ...
%!      'sign', @(x) sign (x); 'symmetric', @(x) 2 * x - 1; ...
%!      'symmetriclogit', @(x) 2 ./ (1 + exp (-x)) - 1};
%! for i = 1:rows (T)
%!   Mdl.ScoreTransform = T{i,1};
%!   [l, s] = predict (Mdl, meas([1, 60, 120],:));
%!   assert_equal (s, T{i,2}(raw), 1e-12);
%!   assert_equal (l, label);
%! endfor
%! ## ismax marks the largest score of each observation, ties to the first.
%! [~, k] = max (raw, [], 2);
%! e = zeros (size (raw));
%! e(sub2ind (size (raw), (1:rows (raw))', k)) = 1;
%! Mdl.ScoreTransform = 'ismax';
%! [~, s] = predict (Mdl, meas([1, 60, 120],:));
%! assert_equal (s, e);
%! Mdl.ScoreTransform = 'symmetricismax';
%! [~, s] = predict (Mdl, meas([1, 60, 120],:));
%! assert_equal (s, 2 * e - 1);

## A function handle is taken as given and applied to the scores.
%!test
%! load fisheriris
%! Mdl = fitclinear (meas, strcmp (species, 'setosa'));
%! Mdl.ScoreTransform = 'none';
%! [label, raw] = predict (Mdl, meas([1, 60, 120],:));
%! Mdl.ScoreTransform = @(x) x .^ 2;
%! [l, s] = predict (Mdl, meas([1, 60, 120],:));
%! assert_equal (s, raw .^ 2, 1e-12);
%! assert_equal (l, label);
