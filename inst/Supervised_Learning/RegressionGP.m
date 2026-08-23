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
## @deftypefn  {statistics} {@var{obj} =} RegressionGP (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} RegressionGP (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{RegressionGP} object containing a Gaussian process
## regression model.
##
## @code{@var{obj} = RegressionGP (@var{X}, @var{Y})} returns a Gaussian
## process regression model, @var{obj}, with @var{X} being the predictor data
## and @var{Y} the continuous response of the observations in @var{X}.
##
## @itemize
## @item
## @var{X} must be an @math{NxP} numeric matrix of predictor data, where rows
## correspond to observations and columns to features.
## @item
## @var{Y} must be an @math{Nx1} numeric vector holding the response of the
## corresponding predictor data in @var{X}.  @var{Y} must have the same number
## of rows as @var{X}.
## @end itemize
##
## A Gaussian process places a prior over functions, given by the covariance
## function, and conditions it on the observations.  The response is modelled
## as @math{H*Beta} plus a draw from that process plus independent noise of
## standard deviation @qcode{Sigma}, where @math{H} is the explicit basis.  The
## covariance parameters and @qcode{Sigma} are estimated by maximizing the log
## marginal likelihood, and @qcode{Beta} follows from them in closed form as
## the generalized least squares estimate.
##
## @code{@var{obj} = RegressionGP (@dots{}, @var{name}, @var{value})} returns a
## model with additional options specified by @qcode{Name-Value} pair
## arguments listed below.
##
## @multitable @columnfractions 0.32 0.68
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'KernelFunction'} @tab A character vector naming the covariance
## function, or a function handle taking two matrices of points and a parameter
## vector.  The default is @qcode{'squaredexponential'}.  The supported names
## are listed below.
##
## @item @qcode{'KernelParameters'} @tab A numeric vector of initial values for
## the covariance parameters.  Its length depends on the covariance function.
## These are starting values for the optimization, not fixed values.
##
## @item @qcode{'BasisFunction'} @tab A character vector naming the explicit
## basis, one of @qcode{'none'}, @qcode{'constant'}, @qcode{'linear'} or
## @qcode{'pureQuadratic'}, or a function handle taking @var{X} and returning
## the basis matrix.  The default is @qcode{'constant'}.
##
## @item @qcode{'Beta'} @tab A numeric vector of basis coefficients.  These are
## used as known values only when @qcode{'FitMethod'} is @qcode{'none'}.
##
## @item @qcode{'Sigma'} @tab A positive scalar, the initial value of the noise
## standard deviation.  The default is @code{std (@var{Y}) / sqrt (2)}.
##
## @item @qcode{'ConstantSigma'} @tab A logical scalar.  When @qcode{true} the
## noise standard deviation is held at its initial value instead of being
## estimated.  The default is @qcode{false}.
##
## @item @qcode{'SigmaLowerBound'} @tab A positive scalar bounding the noise
## standard deviation from below.  The default is
## @code{1e-2 * std (@var{Y})}.
##
## @item @qcode{'FitMethod'} @tab A character vector, either @qcode{'exact'} to
## estimate the parameters or @qcode{'none'} to keep them at their initial
## values.  The default is @qcode{'exact'}.
##
## @item @qcode{'PredictMethod'} @tab A character vector.  Only @qcode{'exact'}
## is implemented, which is also the only method under which a standard
## deviation and a prediction interval are available.
##
## @item @qcode{'Optimizer'} @tab A character vector naming the optimizer used
## to maximize the log marginal likelihood, either @qcode{'quasinewton'} or
## @qcode{'fminunc'}, which are the same thing, or @qcode{'fminsearch'}.  The
## default is @qcode{'quasinewton'}.
##
## @item @qcode{'Standardize'} @tab A logical scalar specifying whether the
## predictor data should be centred and scaled before training.  The same
## transformation is applied by @code{predict}.  The default is @qcode{false}.
##
## @item @qcode{'Weights'} @tab An @math{Nx1} numeric vector of non-negative
## observation weights.  The default is a vector of ones.
##
## @item @qcode{'PredictorNames'} @tab A cell array of character vectors
## naming the predictors, in the order they appear in @var{X}.
##
## @item @qcode{'ResponseName'} @tab A character vector naming the response.
## The default is @qcode{'Y'}.
##
## @item @qcode{'ResponseTransform'} @tab A character vector or a function
## handle applied to the response the model predicts.  The default is
## @qcode{'none'}.
## @end multitable
##
## The supported values for @qcode{'KernelFunction'} are:
##
## @multitable @columnfractions 0.4 0.6
## @headitem @var{Value} @tab @var{Parameters}
## @item @qcode{'exponential'} @tab @qcode{[SigmaL; SigmaF]}
## @item @qcode{'squaredexponential'} @tab @qcode{[SigmaL; SigmaF]}
## @item @qcode{'matern32'} @tab @qcode{[SigmaL; SigmaF]}
## @item @qcode{'matern52'} @tab @qcode{[SigmaL; SigmaF]}
## @item @qcode{'rationalquadratic'} @tab @qcode{[SigmaL; AlphaRQ; SigmaF]}
## @item @qcode{'ardexponential'} @tab @qcode{[LengthScale1; @dots{}; SigmaF]}
## @item @qcode{'ardsquaredexponential'} @tab
## @qcode{[LengthScale1; @dots{}; SigmaF]}
## @item @qcode{'ardmatern32'} @tab @qcode{[LengthScale1; @dots{}; SigmaF]}
## @item @qcode{'ardmatern52'} @tab @qcode{[LengthScale1; @dots{}; SigmaF]}
## @item @qcode{'ardrationalquadratic'} @tab
## @qcode{[LengthScale1; @dots{}; AlphaRQ; SigmaF]}
## @end multitable
##
## The automatic relevance determination kernels carry one length scale per
## predictor, so a predictor the response does not depend on is given a large
## length scale and stops contributing.
##
## The supported values for @qcode{'ResponseTransform'} are:
##
## @multitable @columnfractions 0.3 0.7
## @headitem @var{Value} @tab @var{Description}
## @item @qcode{'none'} @tab @math{x} (no transformation)
## @item @qcode{'identity'} @tab @math{x} (no transformation)
## @item @qcode{'exp'} @tab @math{exp (x)}
## @item @qcode{'log'} @tab @math{log (x)}
## @end multitable
##
## Two deviations from MATLAB are deliberate and documented.  The distance
## between points is accumulated one predictor at a time instead of by the
## expanded form MATLAB uses by default, because the expanded form does not
## return exactly zero for a point against itself and the rough kernels
## amplify that residue through their square root.  The approximate fitting
## and prediction methods, @qcode{'sd'}, @qcode{'sr'}, @qcode{'fic'} and
## @qcode{'bcd'}, together with the active set options that serve them, are
## not implemented and are refused rather than silently ignored.
##
## @seealso{fitrgp, CompactRegressionGP, RegressionSVM, RegressionGAM}
## @end deftypefn

classdef RegressionGP

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} X
    ##
    ## Predictor data
    ##
    ## An @math{NxP} numeric matrix, as it was supplied to the constructor.
    ## This property is read-only.
    ##
    ## @end deftp
    X                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} Y
    ##
    ## Response data
    ##
    ## An @math{Nx1} numeric vector, as it was supplied to the constructor.
    ## This property is read-only.
    ##
    ## @end deftp
    Y                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} NumObservations
    ##
    ## Number of observations used to train the model
    ##
    ## A positive integer scalar, counting only the rows that survived the
    ## removal of missing values.  This property is read-only.
    ##
    ## @end deftp
    NumObservations       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} RowsUsed
    ##
    ## Rows of the original data used to train the model
    ##
    ## A logical vector with one element per row of the data as supplied, true
    ## where the row was used.  It is empty when no row was dropped.  This
    ## property is read-only.
    ##
    ## @end deftp
    RowsUsed              = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} W
    ##
    ## Observation weights
    ##
    ## An @math{Nx1} numeric vector, one weight per observation used to train
    ## the model.  This property is read-only.
    ##
    ## @end deftp
    W                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} PredictorNames
    ##
    ## Predictor variable names
    ##
    ## A cell array of character vectors, one per column of @qcode{X}.  This
    ## property is read-only.
    ##
    ## @end deftp
    PredictorNames        = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} ExpandedPredictorNames
    ##
    ## Expanded predictor variable names
    ##
    ## A cell array of character vectors.  It differs from
    ## @qcode{PredictorNames} only where a categorical predictor has been
    ## expanded into indicator variables.  This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    ResponseName          = 'Y';

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A vector of positive integers indexing the columns of @qcode{X} that
    ## hold categorical predictors, or empty when none does.  This property is
    ## read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## A cell array with one entry per predictor, holding that predictor's bin
    ## edges where the model discretized it before fitting.  It is empty here
    ## and stays empty: a Gaussian process takes its predictors as they are.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    BinEdges              = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} FitMethod
    ##
    ## Method used to estimate the parameters
    ##
    ## @qcode{'Exact'} when the covariance parameters and the noise were
    ## estimated by maximizing the log marginal likelihood, and @qcode{'None'}
    ## when they were kept at their initial values.  This property is
    ## read-only.
    ##
    ## @end deftp
    FitMethod             = 'Exact';

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} BasisFunction
    ##
    ## Explicit basis of the model
    ##
    ## @qcode{'None'}, @qcode{'Constant'}, @qcode{'Linear'},
    ## @qcode{'PureQuadratic'}, or the function handle that was supplied.  This
    ## property is read-only.
    ##
    ## @end deftp
    BasisFunction         = 'Constant';

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} Beta
    ##
    ## Estimated coefficients of the explicit basis
    ##
    ## A numeric vector with one element per basis term, empty when the basis
    ## is @qcode{'None'}.  This property is read-only.
    ##
    ## @end deftp
    Beta                  = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} Sigma
    ##
    ## Estimated noise standard deviation
    ##
    ## A positive scalar.  This property is read-only.
    ##
    ## @end deftp
    Sigma                 = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} LogLikelihood
    ##
    ## Maximized log marginal likelihood
    ##
    ## A scalar, or empty when @qcode{FitMethod} is @qcode{'None'} and nothing
    ## was maximized.  This property is read-only.
    ##
    ## @end deftp
    LogLikelihood         = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} ModelParameters
    ##
    ## Parameters used to train the model
    ##
    ## A structure holding the options the fit was performed under.  MATLAB
    ## returns an object of its own class here; a structure carries the same
    ## information and is what every other learner in this package returns.
    ## This property is read-only.
    ##
    ## @end deftp
    ModelParameters       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} KernelFunction
    ##
    ## Form of the covariance function
    ##
    ## A character vector naming the covariance function, or the function
    ## handle that was supplied.  This property is read-only.
    ##
    ## @end deftp
    KernelFunction        = 'SquaredExponential';

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} KernelInformation
    ##
    ## Covariance function and its parameters
    ##
    ## A structure with fields @qcode{Name}, @qcode{KernelParameters} and
    ## @qcode{KernelParameterNames}, the last naming each parameter in the
    ## order they are stored.  This property is read-only.
    ##
    ## @end deftp
    KernelInformation     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} PredictMethod
    ##
    ## Method used to make predictions
    ##
    ## @qcode{'Exact'}.  This property is read-only.
    ##
    ## @end deftp
    PredictMethod         = 'Exact';

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} Alpha
    ##
    ## Weights the predictions are made from
    ##
    ## An @math{Nx1} numeric vector.  A prediction is the basis term plus the
    ## covariance between the new point and the active set, weighted by these.
    ## This property is read-only.
    ##
    ## @end deftp
    Alpha                 = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} ActiveSetVectors
    ##
    ## Subset of the training data used for predictions
    ##
    ## An @math{MxP} numeric matrix, standardized where the model standardized
    ## its predictors.  It is the whole of the training data, since only the
    ## exact method is implemented.  This property is read-only.
    ##
    ## @end deftp
    ActiveSetVectors      = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} ActiveSetMethod
    ##
    ## Method used to select the active set
    ##
    ## @qcode{'Random'}.  This property is read-only.
    ##
    ## @end deftp
    ActiveSetMethod       = 'Random';

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} ActiveSetSize
    ##
    ## Size of the active set
    ##
    ## A positive integer scalar.  This property is read-only.
    ##
    ## @end deftp
    ActiveSetSize         = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} IsActiveSetVector
    ##
    ## Which observations are in the active set
    ##
    ## A logical vector with one element per training observation.  This
    ## property is read-only.
    ##
    ## @end deftp
    IsActiveSetVector     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} PredictorLocation
    ##
    ## Means the predictors were centred by
    ##
    ## A @math{1xP} numeric vector when the model standardized its predictors,
    ## and empty when it did not.  This property is read-only.
    ##
    ## @end deftp
    PredictorLocation     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} PredictorScale
    ##
    ## Standard deviations the predictors were scaled by
    ##
    ## A @math{1xP} numeric vector when the model standardized its predictors,
    ## and empty when it did not.  This property is read-only.
    ##
    ## @end deftp
    PredictorScale        = [];

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {RegressionGP} {property} ResponseTransform
    ##
    ## Transformation applied to the predicted response
    ##
    ## A character vector, or the text of the function handle that was
    ## supplied.  Assigning to it accepts either.
    ##
    ## @end deftp
    ResponseTransform     = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## The callable behind ResponseTransform.  The public property is the text
    ## MATLAB stores; this is what predict actually applies.
    RTfun                 = @(y) y;

  endproperties

  methods (Access = public)

    ## Custom setter, so that assigning a name or a handle updates both the
    ## text the property reports and the callable predict uses.
    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = ...
                     parseResponseTransform (val, 'RegressionGP');
    endfunction

    ## Class constructor
    function this = RegressionGP (X, Y, varargin)

      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("RegressionGP: too few input arguments.");
      endif

      ## Validate X and Y and drop rows with missing values
      [X, Y, RowsUsed] = this.checkXY_ (X, Y, 'RegressionGP');
      [n, p] = size (X);

      ## Defaults
      KernelFunction   = 'squaredexponential';
      KernelParameters = [];
      BasisFunction    = 'constant';
      BetaIn           = [];
      SigmaIn          = [];
      ConstantSigma    = false;
      SigmaLowerBound  = [];
      FitMethod        = 'exact';
      PredictMethod    = 'exact';
      Optimizer        = 'quasinewton';
      Standardize      = false;
      Weights          = [];
      PredictorNames   = {};
      ResponseName     = 'Y';
      ResponseTransform = 'none';
      CategoricalPredictors = [];

      ## Parse optional parameters
      while (numel (varargin) > 0)
        if (numel (varargin) < 2)
          error (strcat ("RegressionGP: optional arguments must be given", ...
                         " in Name-Value pairs."));
        endif
        switch (lower (varargin{1}))

          case 'kernelfunction'
            KernelFunction = varargin{2};
            if (! (ischar (KernelFunction) || ...
                   is_function_handle (KernelFunction)))
              error (strcat ("RegressionGP: 'KernelFunction' must be a", ...
                             " character vector or a function handle."));
            endif
            if (ischar (KernelFunction) && ...
                ! any (strcmpi (KernelFunction, kernelNames ())))
              error ("RegressionGP: unsupported 'KernelFunction' value.");
            endif

          case 'kernelparameters'
            KernelParameters = varargin{2};
            if (! (isnumeric (KernelParameters) && ...
                   isvector (KernelParameters) && ...
                   all (KernelParameters > 0)))
              error (strcat ("RegressionGP: 'KernelParameters' must be a", ...
                             " vector of positive values."));
            endif
            KernelParameters = KernelParameters(:);

          case 'basisfunction'
            BasisFunction = varargin{2};
            if (! (ischar (BasisFunction) || ...
                   is_function_handle (BasisFunction)))
              error (strcat ("RegressionGP: 'BasisFunction' must be a", ...
                             " character vector or a function handle."));
            endif
            if (ischar (BasisFunction) && ...
                ! any (strcmpi (BasisFunction, ...
                                {'none', 'constant', 'linear', ...
                                 'purequadratic'})))
              error ("RegressionGP: unsupported 'BasisFunction' value.");
            endif

          case 'beta'
            BetaIn = varargin{2};
            if (! (isnumeric (BetaIn) && isvector (BetaIn)))
              error ("RegressionGP: 'Beta' must be a numeric vector.");
            endif
            BetaIn = BetaIn(:);

          case 'sigma'
            SigmaIn = varargin{2};
            if (! (isnumeric (SigmaIn) && isscalar (SigmaIn) && SigmaIn > 0))
              error ("RegressionGP: 'Sigma' must be a positive scalar.");
            endif

          case 'constantsigma'
            ConstantSigma = varargin{2};
            if (! (islogical (ConstantSigma) && isscalar (ConstantSigma)))
              error (strcat ("RegressionGP: 'ConstantSigma' must be a", ...
                             " logical scalar."));
            endif

          case 'sigmalowerbound'
            SigmaLowerBound = varargin{2};
            if (! (isnumeric (SigmaLowerBound) && ...
                   isscalar (SigmaLowerBound) && SigmaLowerBound > 0))
              error (strcat ("RegressionGP: 'SigmaLowerBound' must be a", ...
                             " positive scalar."));
            endif

          case 'fitmethod'
            FitMethod = varargin{2};
            if (! ischar (FitMethod))
              error (strcat ("RegressionGP: 'FitMethod' must be a", ...
                             " character vector."));
            endif
            if (any (strcmpi (FitMethod, {'sd', 'sr', 'fic'})))
              error (strcat ("RegressionGP: the approximate fitting", ...
                             " methods are not implemented; 'FitMethod'", ...
                             " must be either 'exact' or 'none'."));
            endif
            if (! any (strcmpi (FitMethod, {'exact', 'none'})))
              error ("RegressionGP: unsupported 'FitMethod' value.");
            endif

          case 'predictmethod'
            PredictMethod = varargin{2};
            if (! ischar (PredictMethod))
              error (strcat ("RegressionGP: 'PredictMethod' must be a", ...
                             " character vector."));
            endif
            if (any (strcmpi (PredictMethod, {'bcd', 'sd', 'sr', 'fic'})))
              error (strcat ("RegressionGP: the approximate prediction", ...
                             " methods are not implemented;", ...
                             " 'PredictMethod' must be 'exact'."));
            endif
            if (! strcmpi (PredictMethod, 'exact'))
              error ("RegressionGP: unsupported 'PredictMethod' value.");
            endif

          case 'optimizer'
            Optimizer = varargin{2};
            if (! ischar (Optimizer))
              error (strcat ("RegressionGP: 'Optimizer' must be a", ...
                             " character vector."));
            endif
            if (strcmpi (Optimizer, 'fmincon'))
              error (strcat ("RegressionGP: 'fmincon' is not available in", ...
                             " core Octave; use 'quasinewton' or", ...
                             " 'fminsearch'."));
            endif
            if (! any (strcmpi (Optimizer, {'quasinewton', 'fminunc', ...
                                            'lbfgs', 'fminsearch'})))
              error ("RegressionGP: unsupported 'Optimizer' value.");
            endif

          case 'standardize'
            Standardize = varargin{2};
            if (! (islogical (Standardize) && isscalar (Standardize)))
              error (strcat ("RegressionGP: 'Standardize' must be a", ...
                             " logical scalar."));
            endif

          case 'weights'
            Weights = varargin{2};

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! (iscellstr (PredictorNames) && ...
                   numel (PredictorNames) == p))
              error (strcat ("RegressionGP: 'PredictorNames' must be a", ...
                             " cell array of character vectors with one", ...
                             " name per column of X."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("RegressionGP: 'ResponseName' must be a", ...
                             " character vector."));
            endif

          case 'responsetransform'
            ResponseTransform = varargin{2};

          case 'categoricalpredictors'
            CategoricalPredictors = varargin{2};

          otherwise
            error (strcat ("RegressionGP: invalid parameter name in", ...
                           " optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Store the data and its description
      this.X = X;
      this.Y = Y;
      this.NumObservations = n;
      this.RowsUsed = RowsUsed;
      this.W = this.getWeights_ (Weights, n, 'RegressionGP');
      if (isempty (PredictorNames))
        PredictorNames = arrayfun (@(k) sprintf ('x%d', k), 1:p, ...
                                   'UniformOutput', false);
      endif
      this.PredictorNames = PredictorNames;
      this.ExpandedPredictorNames = PredictorNames;
      this.ResponseName = ResponseName;
      this.CategoricalPredictors = CategoricalPredictors;
      this.BinEdges = {};
      this.ResponseTransform = ResponseTransform;

      ## Standardize the predictors, if asked.  The location and scale are
      ## kept so that predict can apply the same transformation.
      XS = X;
      if (Standardize)
        this.PredictorLocation = mean (X, 1);
        s = std (X, 0, 1);
        s(s == 0) = 1;
        this.PredictorScale = s;
        XS = (X - this.PredictorLocation) ./ this.PredictorScale;
      endif

      ## Initial values, which are the documented defaults where the caller
      ## gave none
      ## Every documented default is a multiple of std (Y), which a constant
      ## response makes zero, and the logarithm of zero is not a starting
      ## point.  A unit spread is used instead: the fit is degenerate either
      ## way, but it returns a model rather than a NaN.
      sy = std (Y);
      if (sy == 0)
        sy = 1;
      endif
      if (isempty (SigmaIn))
        SigmaIn = sy / sqrt (2);
      endif
      if (isempty (SigmaLowerBound))
        SigmaLowerBound = 1e-2 * sy;
      endif
      if (isempty (KernelParameters))
        KernelParameters = defaultKernelParameters (XS, Y, KernelFunction);
      endif

      ## Fit
      [theta, sigmaN, beta, LL] = ...
          gpFit (XS, Y, KernelFunction, KernelParameters, BasisFunction, ...
                 BetaIn, SigmaIn, ConstantSigma, SigmaLowerBound, ...
                 FitMethod, Optimizer);

      ## Store the fitted model, in the form MATLAB stores it
      this.FitMethod = properName (FitMethod);
      this.PredictMethod = properName (PredictMethod);
      this.BasisFunction = properName (BasisFunction);
      this.KernelFunction = properName (KernelFunction);
      kpnames = kernelParameterNames (KernelFunction, columns (X));
      KI = struct ();
      KI.Name = properName (KernelFunction);
      KI.KernelParameters = theta;
      KI.KernelParameterNames = kpnames;
      this.KernelInformation = KI;
      this.Beta = beta;
      this.Sigma = sigmaN;
      this.LogLikelihood = LL;
      this.ActiveSetVectors = XS;
      this.ActiveSetSize = n;
      this.IsActiveSetVector = true (n, 1);
      this.ModelParameters = struct ('FitMethod', properName (FitMethod), ...
                                     'PredictMethod', ...
                                     properName (PredictMethod), ...
                                     'BasisFunction', ...
                                     properName (BasisFunction), ...
                                     'KernelFunction', ...
                                     properName (KernelFunction), ...
                                     'Standardize', Standardize, ...
                                     'ConstantSigma', ConstantSigma, ...
                                     'SigmaLowerBound', SigmaLowerBound, ...
                                     'Optimizer', Optimizer);

      ## The prediction weights
      K = gpCovariance (XS, XS, KernelFunction, theta);
      H = gpBasis (XS, this.BasisFunction);
      A = K + sigmaN^2 * eye (n);
      if (isempty (beta))
        this.Alpha = A \ Y;
      else
        this.Alpha = A \ (Y - H * beta);
      endif

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGP} {@var{yFit} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {RegressionGP} {[@var{yFit}, @var{ySD}, @var{yInt}] =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {RegressionGP} {[@dots{}] =} predict (@dots{}, @qcode{'Alpha'}, @var{alpha})
    ##
    ## Predict the response for new data with a Gaussian process model.
    ##
    ## @code{@var{yFit} = predict (@var{obj}, @var{XC})} returns the predicted
    ## response of the @qcode{RegressionGP} model @var{obj} at the points in
    ## @var{XC}, which must have as many columns as the model has predictors.
    ##
    ## @code{[@var{yFit}, @var{ySD}, @var{yInt}] = predict (@dots{})} also
    ## returns the standard deviation of each predicted response and the
    ## prediction intervals.  The standard deviation is that of a new
    ## @emph{response}, so it carries the noise as well as the uncertainty of
    ## the latent function, and the interval is the normal quantile of the
    ## level times it.
    ##
    ## @code{[@dots{}] = predict (@dots{}, @qcode{'Alpha'}, @var{alpha})} sets
    ## the significance level of the intervals, so that they are
    ## @math{100 * (1 - @var{alpha})} per cent intervals.  @var{alpha} must be
    ## a scalar in the range @math{[0, 1]} and defaults to @math{0.05}.
    ##
    ## @end deftypefn
    function [yFit, ySD, yInt] = predict (this, XC, varargin)

      if (nargin < 2)
        error ("RegressionGP.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("RegressionGP.predict: XC is empty.");
      endif
      if (columns (XC) != columns (this.X))
        error (strcat ("RegressionGP.predict: XC must have the same", ...
                       " number of predictors as the trained model."));
      endif

      CIAlpha = 0.05;
      while (numel (varargin) > 0)
        if (numel (varargin) < 2)
          error (strcat ("RegressionGP.predict: optional arguments must", ...
                         " be given in Name-Value pairs."));
        endif
        switch (lower (varargin{1}))
          case 'alpha'
            CIAlpha = varargin{2};
            if (! (isnumeric (CIAlpha) && isscalar (CIAlpha) && ...
                   CIAlpha >= 0 && CIAlpha <= 1))
              error (strcat ("RegressionGP.predict: 'Alpha' must be a", ...
                             " scalar between 0 and 1."));
            endif
          otherwise
            error (strcat ("RegressionGP.predict: invalid NAME in optional", ...
                           " pairs of arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      M = this.predictModel_ (CIAlpha);
      if (nargout < 2)
        yFit = this.RTfun (gpPredict (XC, M));
      elseif (nargout < 3)
        [yFit, ySD] = gpPredict (XC, M);
        yFit = this.RTfun (yFit);
      else
        [yFit, ySD, yInt] = gpPredict (XC, M);
        yFit = this.RTfun (yFit);
        yInt = this.RTfun (yInt);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGP} {@var{yFit} =} resubPredict (@var{obj})
    ## @deftypefnx {RegressionGP} {[@var{yFit}, @var{ySD}, @var{yInt}] =} resubPredict (@var{obj})
    ##
    ## Predict the response of the training data with a Gaussian process model.
    ##
    ## @code{@var{yFit} = resubPredict (@var{obj})} returns the response the
    ## @qcode{RegressionGP} model @var{obj} predicts at its own training data,
    ## and the further outputs are those of @code{predict}.
    ##
    ## @end deftypefn
    function [yFit, ySD, yInt] = resubPredict (this)

      if (nargout < 2)
        yFit = this.predict (this.X);
      elseif (nargout < 3)
        [yFit, ySD] = this.predict (this.X);
      else
        [yFit, ySD, yInt] = this.predict (this.X);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGP} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionGP} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the regression loss of a Gaussian process model.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the mean
    ## squared error of the model @var{obj} on the data @var{X} and @var{Y}.
    ##
    ## @code{@var{L} = loss (@dots{}, @var{name}, @var{value})} accepts
    ## @qcode{'LossFun'}, either @qcode{'mse'}, @qcode{'mae'},
    ## @qcode{'epsiloninsensitive'} or a function handle taking the observed
    ## and the predicted response, and @qcode{'Weights'}, a vector of
    ## non-negative observation weights.
    ##
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      if (nargin < 3)
        error ("RegressionGP.loss: too few input arguments.");
      endif
      [X, Y] = this.checkXY_ (X, Y, 'RegressionGP.loss');

      LossFun = 'mse';
      Weights = ones (rows (X), 1);
      Epsilon = 0;
      while (numel (varargin) > 0)
        if (numel (varargin) < 2)
          error (strcat ("RegressionGP.loss: optional arguments must be", ...
                         " given in Name-Value pairs."));
        endif
        switch (lower (varargin{1}))
          case 'lossfun'
            LossFun = varargin{2};
            if (! (ischar (LossFun) || is_function_handle (LossFun)))
              error (strcat ("RegressionGP.loss: 'LossFun' must be a", ...
                             " character vector or a function handle."));
            endif
            if (ischar (LossFun) && ...
                ! any (strcmpi (LossFun, {'mse', 'mae', ...
                                          'epsiloninsensitive'})))
              error ("RegressionGP.loss: unsupported 'LossFun' value.");
            endif
          case 'weights'
            Weights = varargin{2};
            if (! (isnumeric (Weights) && isvector (Weights) && ...
                   numel (Weights) == rows (X) && all (Weights >= 0)))
              error (strcat ("RegressionGP.loss: 'Weights' must be a", ...
                             " vector of non-negative values with one", ...
                             " element per observation."));
            endif
            Weights = Weights(:);
          case 'epsilon'
            Epsilon = varargin{2};
          otherwise
            error (strcat ("RegressionGP.loss: invalid NAME in optional", ...
                           " pairs of arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      yFit = this.predict (X);
      if (is_function_handle (LossFun))
        L = LossFun (Y, yFit);
        return;
      endif

      switch (lower (LossFun))
        case 'mse'
          L = sum (Weights .* (Y - yFit) .^ 2) / sum (Weights);
        case 'mae'
          L = sum (Weights .* abs (Y - yFit)) / sum (Weights);
        case 'epsiloninsensitive'
          e = max (0, abs (Y - yFit) - Epsilon);
          L = sum (Weights .* e) / sum (Weights);
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGP} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {RegressionGP} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the resubstitution loss of a Gaussian process model.
    ##
    ## @code{@var{L} = resubLoss (@var{obj})} returns the loss of the model
    ## @var{obj} on the data it was trained on, and accepts the same
    ## Name-Value pairs as @code{loss}.
    ##
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      L = this.loss (this.X, this.Y, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGP} {[@var{loores}, @var{neff}] =} postFitStatistics (@var{obj})
    ##
    ## Compute the leave-one-out residuals of a Gaussian process model.
    ##
    ## @code{[@var{loores}, @var{neff}] = postFitStatistics (@var{obj})}
    ## returns the @math{Nx1} vector of leave-one-out residuals of the model
    ## @var{obj}, and the number of effective parameters the fit uses.  Neither
    ## requires refitting the model: both follow from the factorization the fit
    ## already produced.
    ##
    ## The coefficients of the explicit basis are treated as estimated, which
    ## is what @qcode{FitMethod} @qcode{'Exact'} makes them, while the
    ## covariance parameters and the noise are treated as known.
    ##
    ## @end deftypefn
    function [loores, neff] = postFitStatistics (this)

      if (! strcmpi (this.PredictMethod, 'Exact'))
        error (strcat ("RegressionGP.postFitStatistics: post-fit", ...
                       " statistics are available only when PredictMethod", ...
                       " is 'Exact'."));
      endif

      n = this.NumObservations;
      K = gpCovariance (this.ActiveSetVectors, this.ActiveSetVectors, ...
                        this.KernelFunction, ...
                        this.KernelInformation.KernelParameters);
      A = K + this.Sigma^2 * eye (n);
      H = gpBasis (this.ActiveSetVectors, this.BasisFunction);
      [L, p] = chol (A, 'lower');
      if (p != 0)
        error (strcat ("RegressionGP.postFitStatistics: the covariance", ...
                       " matrix is not positive definite."));
      endif
      Li = L \ eye (n);
      Ai = Li' * Li;

      if (isempty (H))
        P = Ai;
      else
        P = Ai - Ai * H * ((H' * Ai * H) \ (H' * Ai));
      endif

      loores = (P * this.Y) ./ diag (P);
      neff = columns (H) + trace (K * P);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGP} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {RegressionGP} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross validate a Gaussian process model.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} returns a
    ## @qcode{RegressionPartitionedModel} built from the model @var{obj} by
    ## ten-fold cross validation.
    ##
    ## @code{@var{CVMdl} = crossval (@dots{}, @var{name}, @var{value})} accepts
    ## @qcode{'KFold'}, @qcode{'Holdout'}, @qcode{'Leaveout'} and
    ## @qcode{'CVPartition'}, of which at most one may be given.
    ##
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)

      numFolds = 10;
      Holdout = [];
      Leaveout = 'off';
      CVPartition = [];
      given = 0;

      while (numel (varargin) > 0)
        if (numel (varargin) < 2)
          error (strcat ("RegressionGP.crossval: optional arguments must", ...
                         " be given in Name-Value pairs."));
        endif
        switch (lower (varargin{1}))
          case 'kfold'
            numFolds = varargin{2};
            if (! (isnumeric (numFolds) && isscalar (numFolds) && ...
                   numFolds == fix (numFolds) && numFolds > 1))
              error (strcat ("RegressionGP.crossval: 'KFold' must be an", ...
                             " integer value greater than 1."));
            endif
            given++;
          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && ...
                   Holdout > 0 && Holdout < 1))
              error (strcat ("RegressionGP.crossval: 'Holdout' must be a", ...
                             " numeric value between 0 and 1."));
            endif
            given++;
          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout) && ...
                   any (strcmpi (Leaveout, {'on', 'off'}))))
              error (strcat ("RegressionGP.crossval: 'Leaveout' must be", ...
                             " either 'on' or 'off'."));
            endif
            given++;
          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("RegressionGP.crossval: 'CVPartition' must", ...
                             " be a cvpartition object."));
            endif
            given++;
          otherwise
            error (strcat ("RegressionGP.crossval: invalid parameter name", ...
                           " in optional paired arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      if (given > 1)
        error (strcat ("RegressionGP.crossval: you can use only one of", ...
                       " 'KFold', 'Holdout', 'Leaveout', or", ...
                       " 'CVPartition' options."));
      endif

      n = this.NumObservations;
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (n, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (n, 'LeaveOut');
      else
        partition = cvpartition (n, 'KFold', numFolds);
      endif

      CVMdl = RegressionPartitionedModel (this, partition);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGP} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Return a compact Gaussian process regression model.
    ##
    ## @code{@var{CMdl} = compact (@var{obj})} returns a
    ## @qcode{CompactRegressionGP} object holding what is needed to predict
    ## and nothing else: the training data, the response and everything that
    ## describes them are dropped.
    ##
    ## @end deftypefn
    function CMdl = compact (this)
      CMdl = CompactRegressionGP (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGP} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a Gaussian process model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves the model @var{obj}
    ## into @var{filename} in a form @code{loadmodel} can read back.
    ##
    ## @end deftypefn
    function savemodel (obj, fname)

      classdef_name = 'RegressionGP';
      X = obj.X;
      Y = obj.Y;
      NumObservations = obj.NumObservations;
      RowsUsed = obj.RowsUsed;
      W = obj.W;
      PredictorNames = obj.PredictorNames;
      ExpandedPredictorNames = obj.ExpandedPredictorNames;
      ResponseName = obj.ResponseName;
      CategoricalPredictors = obj.CategoricalPredictors;
      BinEdges = obj.BinEdges;
      FitMethod = obj.FitMethod;
      BasisFunction = obj.BasisFunction;
      Beta = obj.Beta;
      Sigma = obj.Sigma;
      LogLikelihood = obj.LogLikelihood;
      ModelParameters = obj.ModelParameters;
      KernelFunction = obj.KernelFunction;
      KernelInformation = obj.KernelInformation;
      PredictMethod = obj.PredictMethod;
      Alpha = obj.Alpha;
      ActiveSetVectors = obj.ActiveSetVectors;
      ActiveSetMethod = obj.ActiveSetMethod;
      ActiveSetSize = obj.ActiveSetSize;
      IsActiveSetVector = obj.IsActiveSetVector;
      PredictorLocation = obj.PredictorLocation;
      PredictorScale = obj.PredictorScale;
      ResponseTransform = obj.ResponseTransform;

      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'NumObservations', ...
            'RowsUsed', 'W', 'PredictorNames', 'ExpandedPredictorNames', ...
            'ResponseName', 'CategoricalPredictors', 'BinEdges', ...
            'FitMethod', 'BasisFunction', 'Beta', 'Sigma', 'LogLikelihood', ...
            'ModelParameters', 'KernelFunction', 'KernelInformation', ...
            'PredictMethod', 'Alpha', 'ActiveSetVectors', ...
            'ActiveSetMethod', 'ActiveSetSize', 'IsActiveSetVector', ...
            'PredictorLocation', 'PredictorScale', 'ResponseTransform');

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
      printf ("\n  RegressionGP\n\n");
      printf ("%25s: '%s'\n", 'ResponseName', this.ResponseName);
      printf ("%25s: %d\n", 'NumObservations', this.NumObservations);
      printf ("%25s: %d\n", 'NumPredictors', columns (this.X));
      printf ("%25s: '%s'\n", 'KernelFunction', ...
              nameOf (this.KernelFunction));
      printf ("%25s: '%s'\n", 'BasisFunction', nameOf (this.BasisFunction));
      printf ("%25s: '%s'\n", 'FitMethod', this.FitMethod);
      printf ("%25s: '%s'\n", 'PredictMethod', this.PredictMethod);
      printf ("%25s: %g\n", 'Sigma', this.Sigma);
      if (! isempty (this.LogLikelihood))
        printf ("%25s: %g\n", 'LogLikelihood', this.LogLikelihood);
      endif
      printf ("\n");
    endfunction

    ## The structure gpPredict consumes.  It is assembled here so that the
    ## compact class and this one predict through the same code.
    function M = predictModel_ (this, CIAlpha)
      M = struct ('X', this.ActiveSetVectors, 'Alpha', this.Alpha, ...
                  'KernelFunction', this.KernelFunction, ...
                  'Theta', this.KernelInformation.KernelParameters, ...
                  'BasisFunction', this.BasisFunction, 'Beta', this.Beta, ...
                  'Sigma', this.Sigma, 'Location', this.PredictorLocation, ...
                  'Scale', this.PredictorScale, 'CIAlpha', CIAlpha);
    endfunction

  endmethods

  methods (Access = private)

    ## Validate the predictor and response data, and drop any row that is not
    ## complete in both.
    function [X, Y, used] = checkXY_ (this, X, Y, caller)

      if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
        error ("%s: invalid values in X.", caller);
      endif
      if (! (isnumeric (Y) && isreal (Y) && isvector (Y)))
        error ("%s: invalid values in Y.", caller);
      endif
      Y = Y(:);
      if (rows (X) != rows (Y))
        error ("%s: number of rows in X and Y must be equal.", caller);
      endif
      if (isempty (X))
        error ("%s: X is empty.", caller);
      endif

      used = ! (any (isnan (X), 2) | isnan (Y));
      if (all (used))
        used = [];
      else
        X = X(used, :);
        Y = Y(used);
        if (isempty (Y))
          error ("%s: no complete observations in the data.", caller);
        endif
      endif

    endfunction

    ## Validate observation weights, defaulting to a vector of ones.
    function W = getWeights_ (this, Weights, n, caller)
      if (isempty (Weights))
        W = ones (n, 1);
        return;
      endif
      if (! (isnumeric (Weights) && isvector (Weights) && ...
             numel (Weights) == n && all (Weights >= 0)))
        error (strcat ("%s: 'Weights' must be a vector of non-negative", ...
                       " values with one element per observation."), caller);
      endif
      W = Weights(:);
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      mdl = RegressionGP (data.X, data.Y);
      fields = fieldnames (data);
      for k = 1:numel (fields)
        mdl.(fields{k}) = data.(fields{k});
      endfor

    endfunction

  endmethods

endclassdef

## The proper name of a value, as the model reports it: MATLAB accepts the
## lowercase spelling and stores the capitalized one.
function s = properName (v)
  if (is_function_handle (v))
    s = v;
    return;
  endif
  known = {'none', 'None'; 'exact', 'Exact'; 'constant', 'Constant'; ...
           'linear', 'Linear'; 'purequadratic', 'PureQuadratic'; ...
           'random', 'Random'; ...
           'exponential', 'Exponential'; ...
           'squaredexponential', 'SquaredExponential'; ...
           'matern32', 'Matern32'; 'matern52', 'Matern52'; ...
           'rationalquadratic', 'RationalQuadratic'; ...
           'ardexponential', 'ARDExponential'; ...
           'ardsquaredexponential', 'ARDSquaredExponential'; ...
           'ardmatern32', 'ARDMatern32'; 'ardmatern52', 'ARDMatern52'; ...
           'ardrationalquadratic', 'ARDRationalQuadratic'};
  idx = find (strcmpi (v, known(:,1)), 1);
  if (isempty (idx))
    s = v;
  else
    s = known{idx, 2};
  endif
endfunction

function s = nameOf (v)
  if (is_function_handle (v))
    s = func2str (v);
  else
    s = v;
  endif
endfunction

function names = kernelNames ()
  names = {'exponential', 'squaredexponential', 'matern32', 'matern52', ...
           'rationalquadratic', 'ardexponential', 'ardsquaredexponential', ...
           'ardmatern32', 'ardmatern52', 'ardrationalquadratic'};
endfunction

## The names of a kernel's parameters, in the order they are stored.
function names = kernelParameterNames (kern, d)
  if (is_function_handle (kern))
    names = {};
    return;
  endif
  if (strncmpi (kern, 'ard', 3))
    names = arrayfun (@(k) sprintf ('LengthScale%d', k), 1:d, ...
                      'UniformOutput', false)';
    if (strcmpi (kern, 'ardrationalquadratic'))
      names = [names; {'AlphaRQ'}];
    endif
    names = [names; {'SigmaF'}];
  else
    if (strcmpi (kern, 'rationalquadratic'))
      names = {'SigmaL'; 'AlphaRQ'; 'SigmaF'};
    else
      names = {'SigmaL'; 'SigmaF'};
    endif
  endif
endfunction

## The documented initial values of the covariance parameters.
function theta0 = defaultKernelParameters (X, Y, kern)
  if (is_function_handle (kern))
    error (strcat ("RegressionGP: 'KernelParameters' must be given when", ...
                   " 'KernelFunction' is a function handle."));
  endif
  sf = std (Y) / sqrt (2);
  if (sf == 0)
    sf = 1 / sqrt (2);
  endif
  if (strncmpi (kern, 'ard', 3))
    theta0 = std (X, 0, 1)';
    if (strcmpi (kern, 'ardrationalquadratic'))
      theta0 = [theta0; 1];
    endif
    theta0 = [theta0; sf];
  else
    theta0 = mean (std (X, 0, 1));
    if (strcmpi (kern, 'rationalquadratic'))
      theta0 = [theta0; 1];
    endif
    theta0 = [theta0; sf];
  endif
  theta0(theta0 <= 0) = 1;
endfunction

## Covariance between two sets of points, built in or supplied.
function K = gpCovariance (A, B, kern, theta)
  if (is_function_handle (kern))
    K = kern (A, B, theta);
  else
    K = gpKernel (A, B, kern, theta);
  endif
endfunction

## Fit the model.  With FitMethod 'none' the initial values are kept and only
## Beta is computed; with 'exact' the covariance parameters and the noise are
## estimated by maximizing the log marginal likelihood, with Beta following in
## closed form at every step, so the optimizer never carries it.
function [theta, sigmaN, beta, LL] = gpFit (X, Y, kern, theta0, basis, ...
                                            betaIn, sigma0, constSigma, ...
                                            sigmaLB, fitMethod, optimizer)

  H = gpBasis (X, basis);

  if (strcmpi (fitMethod, 'none'))
    theta = theta0;
    sigmaN = sigma0;
    if (isempty (betaIn))
      beta = zeros (columns (H), 1);
    else
      if (numel (betaIn) != columns (H))
        error (strcat ("RegressionGP: 'Beta' must have one element per", ...
                       " basis term."));
      endif
      beta = betaIn;
    endif
    LL = [];
    return;
  endif

  if (sigma0 <= sigmaLB)
    sigma0 = sigmaLB * (1 + 1e-6);
  endif

  ## The optimizer works on an unconstrained scale: the covariance parameters
  ## through their logarithm, the noise as the lower bound plus an exponential,
  ## which is the parameterization MATLAB documents.
  if (constSigma)
    u0 = log (theta0);
  else
    u0 = [log(theta0); log(sigma0 - sigmaLB)];
  endif

  f = @(u) gpObjective (u, X, Y, kern, H, sigmaLB, constSigma, sigma0);

  if (strcmpi (optimizer, 'fminsearch'))
    opts = optimset ('MaxFunEvals', 5000, 'MaxIter', 2000, 'TolX', 1e-10, ...
                     'TolFun', 1e-10);
    u = fminsearch (@(u) f(u), u0, opts);
  else
    opts = optimset ('GradObj', 'on', 'MaxFunEvals', 5000, ...
                     'MaxIter', 1000, 'TolX', 1e-12, 'TolFun', 1e-12);
    u = fminunc (f, u0, opts);
  endif

  if (constSigma)
    theta = exp (u);
    sigmaN = sigma0;
  else
    theta = exp (u(1:end-1));
    sigmaN = sigmaLB + exp (u(end));
  endif

  [nll, beta] = gpNegLogLik (u, X, Y, kern, H, sigmaLB, constSigma, sigma0);
  LL = -nll;

endfunction

## The objective as an optimizer expects it: the value first and the gradient
## second, where gpNegLogLik returns Beta in between for the caller that wants
## it.
function [nll, g] = gpObjective (u, X, Y, kern, H, sigmaLB, constSigma, sigma0)
  if (nargout > 1)
    [nll, ~, g] = gpNegLogLik (u, X, Y, kern, H, sigmaLB, constSigma, sigma0);
  else
    nll = gpNegLogLik (u, X, Y, kern, H, sigmaLB, constSigma, sigma0);
  endif
endfunction

## The negative profile log marginal likelihood and its gradient.  Beta sits at
## its generalized least squares optimum throughout, so by the envelope theorem
## it contributes nothing to the gradient.
function [nll, beta, g] = gpNegLogLik (u, X, Y, kern, H, sigmaLB, ...
                                       constSigma, sigma0)

  n = rows (X);
  if (constSigma)
    theta = exp (u);
    sigmaN = sigma0;
  else
    theta = exp (u(1:end-1));
    sigmaN = sigmaLB + exp (u(end));
  endif

  K = gpCovariance (X, X, kern, theta);
  A = K + sigmaN^2 * eye (n);

  [L, p] = chol (A, 'lower');
  if (p != 0)
    nll = 1e10;
    beta = zeros (columns (H), 1);
    g = zeros (size (u));
    return;
  endif

  Li_Y = L \ Y;
  if (isempty (H))
    beta = zeros (0, 1);
    v = Li_Y;
  else
    ## Beta is the least squares solution of the whitened system, taken
    ## directly rather than through the normal equations: forming Li_H' * Li_H
    ## squares the condition number, and a basis whose columns are nearly
    ## dependent then warns from inside the objective on every iteration.
    Li_H = L \ H;
    beta = Li_H \ Li_Y;
    v = Li_Y - Li_H * beta;
  endif
  nll = 0.5 * (v' * v) + sum (log (diag (L))) + 0.5 * n * log (2*pi);

  if (nargout < 3)
    return;
  endif

  r = Y;
  if (! isempty (H))
    r = Y - H * beta;
  endif
  ## The inverse is taken through the factor that is already to hand: inv ()
  ## on a covariance whose noise has been driven small warns about a
  ## conditioning the caller can do nothing about, and the triangular solves
  ## are both quieter and cheaper.
  Li = L \ eye (n);
  Ai = Li' * Li;
  a = Ai * r;
  W = Ai - a * a';

  g = zeros (size (u));
  if (! is_function_handle (kern))
    dK = gpKernelDeriv (X, kern, theta, K);
    for j = 1:numel (theta)
      g(j) = 0.5 * sum (sum (W .* dK{j})) * theta(j);
    endfor
  else
    ## A supplied kernel has no derivative to hand, so its parameters are
    ## differenced instead.
    for j = 1:numel (theta)
      h = 1e-6 * max (1, abs (u(j)));
      up = u; up(j) += h;
      um = u; um(j) -= h;
      g(j) = (gpNegLogLik (up, X, Y, kern, H, sigmaLB, constSigma, sigma0) - ...
              gpNegLogLik (um, X, Y, kern, H, sigmaLB, constSigma, sigma0)) ...
             / (2*h);
    endfor
  endif
  if (! constSigma)
    g(end) = 0.5 * trace (W) * 2 * sigmaN * (sigmaN - sigmaLB);
  endif

endfunction

## Derivatives of the covariance matrix with respect to each of its parameters.
function dK = gpKernelDeriv (X, kern, theta, K)

  d = columns (X);
  ard = (numel (kern) > 3 && strncmpi (kern, 'ARD', 3));
  base = lower (strrep (lower (kern), 'ard', ''));

  if (ard)
    L = theta(1:d)(:)';
    if (strcmp (base, 'rationalquadratic'))
      alphaRQ = theta(d+1);
      sigmaF = theta(d+2);
    else
      sigmaF = theta(d+1);
    endif
  else
    L = repmat (theta(1), 1, d);
    if (strcmp (base, 'rationalquadratic'))
      alphaRQ = theta(2);
      sigmaF = theta(3);
    else
      sigmaF = theta(2);
    endif
  endif

  n = rows (X);
  S = cell (1, d);
  r2 = zeros (n, n);
  for j = 1:d
    S{j} = ((X(:,j) - X(:,j)') / L(j)) .^ 2;
    r2 += S{j};
  endfor
  r = sqrt (r2);
  zero = (r2 == 0);

  switch (base)
    case 'exponential'
      dKdr2 = -K ./ (2 * r);
      dKdr2(zero) = 0;
    case 'squaredexponential'
      dKdr2 = -K / 2;
    case 'matern32'
      s = sqrt (3) * r;
      dKdr2 = -sigmaF^2 * (s .^ 2) .* exp (-s) ./ (2 * r2);
      dKdr2(zero) = -sigmaF^2 * 3 / 2;
    case 'matern52'
      s = sqrt (5) * r;
      dKdr2 = -sigmaF^2 * (s .^ 2) .* (1 + s) .* exp (-s) ./ (6 * r2);
      dKdr2(zero) = -sigmaF^2 * 5 / 6;
    case 'rationalquadratic'
      Q = 1 + r2 / (2 * alphaRQ);
      dKdr2 = -sigmaF^2 * Q .^ (-alphaRQ - 1) / 2;
  endswitch

  dK = {};
  if (ard)
    for j = 1:d
      dK{end+1} = dKdr2 .* (-2 * S{j} / L(j));
    endfor
  else
    dK{end+1} = dKdr2 .* (-2 * r2 / L(1));
  endif
  if (strcmp (base, 'rationalquadratic'))
    Q = 1 + r2 / (2 * alphaRQ);
    dK{end+1} = K .* (-log (Q) + r2 ./ (2 * alphaRQ * Q));
  endif
  dK{end+1} = 2 * K / sigmaF;

endfunction

## Every expected value below was measured on MATLAB R2024a.  The tolerances
## are set by what this implementation actually achieves: the covariance
## parameters and the noise agree to about nine digits, because the optimizer
## here converges to a log marginal likelihood a little above MATLAB's, and
## every prediction inherits that.

%!demo
%! ## Fit a Gaussian process to noisy observations of a smooth function and
%! ## show the prediction interval widening away from the data.
%! x = [0.1; 0.3; 0.4; 0.7; 0.9; 1.4; 1.6; 1.9]';
%! x = x(:);
%! y = sin (3 * x) + 0.05 * cos (11 * x);
%! Mdl = fitrgp (x, y);
%! xq = linspace (-0.2, 2.2, 200)';
%! [yq, ~, yint] = predict (Mdl, xq);
%! figure ('visible', 'off');
%! hold on;
%! plot (xq, yint(:,1), 'r:');
%! plot (xq, yint(:,2), 'r:');
%! plot (xq, yq, 'b-');
%! plot (x, y, 'ko');
%! hold off;
%! xlabel ('x');
%! ylabel ('y');
%! title ('Gaussian process regression with 95% prediction interval');

%!test
%! ## The model reports the surface MATLAB reports
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.1 * cos (7*x);
%! Mdl = RegressionGP (x, y);
%! assert_equal (class (Mdl), 'RegressionGP');
%! assert_equal (Mdl.NumObservations, 20);
%! assert_equal (Mdl.FitMethod, 'Exact');
%! assert_equal (Mdl.PredictMethod, 'Exact');
%! assert_equal (Mdl.BasisFunction, 'Constant');
%! assert_equal (Mdl.KernelFunction, 'SquaredExponential');
%! assert_equal (Mdl.ResponseName, 'Y');
%! assert_equal (Mdl.PredictorNames, {'x1'});
%! assert_equal (Mdl.ResponseTransform, 'none');

%!test
%! ## BinEdges is an empty cell, as it is for every learner that does no
%! ## binning: a Gaussian process takes its predictors as they are
%! x = linspace (0, 1, 20)';
%! Mdl = RegressionGP (x, sin (2*pi*x));
%! assert_equal (class (Mdl.BinEdges), 'cell');
%! assert_equal (Mdl.BinEdges, {});

%!test
%! ## The fitted covariance parameters and noise match R2024a
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.1 * cos (7*x);
%! Mdl = RegressionGP (x, y);
%! assert_equal (Mdl.KernelInformation.KernelParameters, ...
%!               [0.386370514454926; 1.505132511329997], 1e-8);
%! assert_equal (Mdl.Beta, -0.072781321631854, 1e-7);
%! assert_equal (Mdl.Sigma, 0.006877618228324, 1e-7);
%! assert_equal (Mdl.LogLikelihood, 47.574057509534278, 1e-4);

%!test
%! ## predict reproduces the R2024a values, and the standard deviation is
%! ## that of a response, so it carries the noise
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.1 * cos (7*x);
%! Mdl = RegressionGP (x, y);
%! xq = [0.05; 0.33; 0.5; 0.77; 0.95];
%! [yp, ysd, yint] = predict (Mdl, xq);
%! assert_equal (yp, [0.404401855406521; 0.809355242891053; ...
%!                    -0.093708552144171; -0.928420803239529; ...
%!                    -0.217104024154071], 1e-7);
%! assert_equal (ysd, [0.008117035802727; 0.007747238553279; ...
%!                     0.007727093294465; 0.007796728197045; ...
%!                     0.008117035802754], 1e-6);
%! assert_equal (size (yint), [5, 2]);
%! assert (all (yint(:,1) < yp));
%! assert (all (yint(:,2) > yp));

%!test
%! ## The interval is the normal quantile of the level times the standard
%! ## deviation, so a looser level gives a narrower interval
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.1 * cos (7*x);
%! Mdl = RegressionGP (x, y);
%! xq = [0.05; 0.33; 0.5];
%! [yp, ysd, yint] = predict (Mdl, xq);
%! assert_equal (yint(:,2) - yp, norminv (0.975) * ysd, 1e-12);
%! [~, ~, yint90] = predict (Mdl, xq, 'Alpha', 0.10);
%! assert_equal (yint90(:,2) - yp, norminv (0.95) * ysd, 1e-12);
%! assert (all (yint90(:,2) - yint90(:,1) < yint(:,2) - yint(:,1)));

%!test
%! ## FitMethod 'none' keeps the documented initial values, estimates nothing
%! ## and has no likelihood to report
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.1 * cos (7*x);
%! Mdl = RegressionGP (x, y, 'FitMethod', 'none');
%! assert_equal (Mdl.FitMethod, 'None');
%! assert_equal (Mdl.Sigma, std (y) / sqrt (2), 1e-14);
%! assert_equal (Mdl.KernelInformation.KernelParameters, ...
%!               [mean(std (x)); std(y) / sqrt(2)], 1e-14);
%! assert_equal (Mdl.Beta, 0);
%! assert_equal (Mdl.LogLikelihood, []);

%!test
%! ## The parameters of each covariance function are named as MATLAB names
%! ## them, and there is one length scale per predictor for the ARD kernels
%! X = [linspace(0, 1, 15)', cos(linspace(0, 3, 15))'];
%! y = X(:,1) .^ 2 + 0.3 * X(:,2);
%! M1 = RegressionGP (X, y, 'KernelFunction', 'squaredexponential');
%! assert_equal (M1.KernelInformation.KernelParameterNames, ...
%!               {'SigmaL'; 'SigmaF'});
%! M2 = RegressionGP (X, y, 'KernelFunction', 'rationalquadratic');
%! assert_equal (M2.KernelInformation.KernelParameterNames, ...
%!               {'SigmaL'; 'AlphaRQ'; 'SigmaF'});
%! M3 = RegressionGP (X, y, 'KernelFunction', 'ardmatern32');
%! assert_equal (M3.KernelInformation.KernelParameterNames, ...
%!               {'LengthScale1'; 'LengthScale2'; 'SigmaF'});
%! assert_equal (numel (M3.KernelInformation.KernelParameters), 3);

%!test
%! ## The name of the covariance function is stored capitalized, as MATLAB
%! ## stores it, whatever spelling the caller used
%! x = linspace (0, 1, 12)';
%! y = cos (3*x);
%! Mdl = RegressionGP (x, y, 'KernelFunction', 'ardmatern52');
%! assert_equal (Mdl.KernelFunction, 'ARDMatern52');
%! assert_equal (Mdl.KernelInformation.Name, 'ARDMatern52');

%!test
%! ## The explicit basis has one coefficient per term, and none has none
%! X = [linspace(0, 1, 15)', cos(linspace(0, 3, 15))'];
%! y = X(:,1) .^ 2 + 0.3 * X(:,2);
%! assert_equal (numel (RegressionGP (X, y, 'BasisFunction', 'none').Beta), 0);
%! assert_equal (numel (RegressionGP (X, y, ...
%!                                    'BasisFunction', 'constant').Beta), 1);
%! assert_equal (numel (RegressionGP (X, y, ...
%!                                    'BasisFunction', 'linear').Beta), 3);
%! assert_equal (numel (RegressionGP (X, y, ...
%!                        'BasisFunction', 'pureQuadratic').Beta), 5);

%!test
%! ## Standardizing records the location and scale it used, and predict
%! ## applies the same transformation
%! X = [linspace(0, 10, 20)', linspace(-5, 5, 20)'];
%! y = 0.3 * X(:,1) - 0.2 * X(:,2);
%! Mdl = RegressionGP (X, y, 'Standardize', true);
%! assert_equal (Mdl.PredictorLocation, mean (X, 1), 1e-14);
%! assert_equal (Mdl.PredictorScale, std (X, 0, 1), 1e-14);
%! assert_equal (predict (Mdl, X(1:3,:)), y(1:3), 1e-3);

%!test
%! ## Without standardizing there is no location or scale to report
%! x = linspace (0, 1, 12)';
%! Mdl = RegressionGP (x, cos (3*x));
%! assert_equal (Mdl.PredictorLocation, []);
%! assert_equal (Mdl.PredictorScale, []);

%!test
%! ## A held noise is not estimated
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! Mdl = RegressionGP (x, y, 'ConstantSigma', true, 'Sigma', 0.3);
%! assert_equal (Mdl.Sigma, 0.3, 1e-14);

%!test
%! ## The noise cannot go below its lower bound
%! x = linspace (0, 1, 15)';
%! y = cos (3*x);
%! Mdl = RegressionGP (x, y, 'SigmaLowerBound', 0.05);
%! assert (Mdl.Sigma >= 0.05);

%!test
%! ## The active set of an exactly fitted model is the whole training data
%! x = linspace (0, 1, 15)';
%! Mdl = RegressionGP (x, cos (3*x));
%! assert_equal (Mdl.ActiveSetSize, 15);
%! assert_equal (Mdl.ActiveSetVectors, x);
%! assert_equal (Mdl.IsActiveSetVector, true (15, 1));
%! assert_equal (Mdl.ActiveSetMethod, 'Random');

%!test
%! ## resubPredict is predict on the training data
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! Mdl = RegressionGP (x, y);
%! assert_equal (resubPredict (Mdl), predict (Mdl, x), 1e-14);

%!test
%! ## The default loss is the mean squared error, and it is small for a model
%! ## that interpolates its own training data closely
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.1 * cos (7*x);
%! Mdl = RegressionGP (x, y);
%! assert_equal (loss (Mdl, x, y), mean ((y - predict (Mdl, x)) .^ 2), 1e-14);
%! assert_equal (resubLoss (Mdl), loss (Mdl, x, y), 1e-14);
%! assert (resubLoss (Mdl) < 1e-4);

%!test
%! ## The other loss functions, and weights
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! Mdl = RegressionGP (x, y);
%! r = y - predict (Mdl, x);
%! assert_equal (loss (Mdl, x, y, 'LossFun', 'mae'), mean (abs (r)), 1e-14);
%! w = linspace (1, 2, 15)';
%! assert_equal (loss (Mdl, x, y, 'Weights', w), ...
%!               sum (w .* r .^ 2) / sum (w), 1e-14);
%! mae = @(a, b) mean (abs (a - b));
%! assert_equal (loss (Mdl, x, y, 'LossFun', mae), mean (abs (r)), 1e-14);

%!test
%! ## The leave-one-out residuals and the effective parameter count match
%! ## R2024a, and neither refits the model
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.1 * cos (7*x);
%! Mdl = RegressionGP (x, y);
%! [loores, neff] = postFitStatistics (Mdl);
%! assert_equal (size (loores), [20, 1]);
%! assert_equal (loores(1:4), [0.006483807954955; -0.002405391043015; ...
%!                             -0.000899188010473; 0.000975297024791], 1e-6);
%! assert_equal (neff, 7.186830203018070, 1e-5);

%!test
%! ## The leave-one-out residuals are larger than the fitted residuals, since
%! ## each leaves out the observation it is about
%! x = linspace (0, 1, 20)';
%! y = sin (2*pi*x) + 0.1 * cos (7*x);
%! Mdl = RegressionGP (x, y);
%! loores = postFitStatistics (Mdl);
%! assert (all (abs (loores) >= abs (y - predict (Mdl, x)) - 1e-12));

%!test
%! ## compact keeps what predicts and drops the rest, and predicts the same
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! Mdl = RegressionGP (x, y);
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), 'CompactRegressionGP');
%! assert_equal (predict (CMdl, x), predict (Mdl, x), 1e-14);
%! [y1, s1] = predict (Mdl, x);
%! [y2, s2] = predict (CMdl, x);
%! assert_equal (s1, s2, 1e-14);

%!test
%! ## crossval returns a partitioned model over the observations trained on
%! x = linspace (0, 1, 20)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! Mdl = RegressionGP (x, y);
%! CVMdl = crossval (Mdl, 'KFold', 4);
%! assert_equal (class (CVMdl), 'RegressionPartitionedModel');
%! assert_equal (CVMdl.KFold, 4);
%! assert_equal (CVMdl.CrossValidatedModel, 'RegressionGP');
%! assert_equal (numel (CVMdl.Trained), 4);
%! assert_equal (class (CVMdl.Trained{1}), 'CompactRegressionGP');
%! assert_equal (size (kfoldPredict (CVMdl)), [20, 1]);

%!test
%! ## A model saved and loaded predicts what it predicted before
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! Mdl = RegressionGP (x, y);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (Mdl2), 'RegressionGP');
%! assert_equal (Mdl2.NumObservations, Mdl.NumObservations);
%! assert_equal (Mdl2.Beta, Mdl.Beta);
%! assert_equal (Mdl2.Sigma, Mdl.Sigma);
%! assert_equal (Mdl2.KernelInformation.KernelParameters, ...
%!               Mdl.KernelInformation.KernelParameters);
%! assert_equal (Mdl2.LogLikelihood, Mdl.LogLikelihood);
%! assert_equal (predict (Mdl2, x), predict (Mdl, x), 1e-14);

%!test
%! ## A response transform is applied to the prediction and reported as text
%! x = linspace (0, 1, 12)';
%! y = cos (3*x);
%! Mdl = RegressionGP (x, y, 'ResponseTransform', 'exp');
%! assert_equal (Mdl.ResponseTransform, 'exp');
%! Mdl2 = RegressionGP (x, y);
%! assert_equal (predict (Mdl, x), exp (predict (Mdl2, x)), 1e-12);

%!test
%! ## A row with a missing value is dropped, and RowsUsed says which
%! x = linspace (0, 1, 12)';
%! y = cos (3*x);
%! x(4) = NaN;
%! Mdl = RegressionGP (x, y);
%! assert_equal (Mdl.NumObservations, 11);
%! assert_equal (Mdl.RowsUsed, [true(3,1); false; true(8,1)]);

%!test
%! ## With no missing value RowsUsed is empty, as it is for every learner
%! x = linspace (0, 1, 12)';
%! Mdl = RegressionGP (x, cos (3*x));
%! assert_equal (Mdl.RowsUsed, []);

%!test
%! ## Observation weights default to one apiece
%! x = linspace (0, 1, 12)';
%! Mdl = RegressionGP (x, cos (3*x));
%! assert_equal (Mdl.W, ones (12, 1));

%!test
%! ## A supplied covariance function is used, and reproduces the built-in one
%! ## it imitates
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! kfcn = @(a, b, th) th(2)^2 * exp (-0.5 * ((a - b') / th(1)) .^ 2);
%! M1 = RegressionGP (x, y, 'KernelFunction', kfcn, ...
%!                    'KernelParameters', [0.4; 1.2], 'FitMethod', 'none');
%! M2 = RegressionGP (x, y, 'KernelFunction', 'squaredexponential', ...
%!                    'KernelParameters', [0.4; 1.2], 'FitMethod', 'none');
%! assert_equal (predict (M1, x), predict (M2, x), 1e-12);

## Test input validation for the constructor
%!error<RegressionGP: too few input arguments.> RegressionGP (ones (5, 2))
%!error<RegressionGP: number of rows in X and Y must be equal.> ...
%! RegressionGP (ones (5, 2), ones (4, 1))
%!error<RegressionGP: invalid values in X.> ...
%! RegressionGP ('a', ones (5, 1))
%!error<RegressionGP: invalid values in Y.> ...
%! RegressionGP (ones (5, 2), 'a')
%!error<RegressionGP: optional arguments must be given in Name-Value pairs.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'Standardize')
%!error<RegressionGP: invalid parameter name in optional pair arguments.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'bogus', 1)
%!error<RegressionGP: 'KernelFunction' must be a character vector or a function handle.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'KernelFunction', 5)
%!error<RegressionGP: unsupported 'KernelFunction' value.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'KernelFunction', 'bogus')
%!error<RegressionGP: 'KernelParameters' must be a vector of positive values.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'KernelParameters', [-1, 2])
%!error<RegressionGP: 'BasisFunction' must be a character vector or a function handle.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'BasisFunction', 5)
%!error<RegressionGP: unsupported 'BasisFunction' value.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'BasisFunction', 'bogus')
%!error<RegressionGP: 'Sigma' must be a positive scalar.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'Sigma', -1)
%!error<RegressionGP: 'ConstantSigma' must be a logical scalar.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'ConstantSigma', 5)
%!error<RegressionGP: 'SigmaLowerBound' must be a positive scalar.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'SigmaLowerBound', 0)
%!error<RegressionGP: the approximate fitting methods are not implemented; 'FitMethod' must be either 'exact' or 'none'.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'FitMethod', 'fic')
%!error<RegressionGP: unsupported 'FitMethod' value.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'FitMethod', 'bogus')
%!error<RegressionGP: the approximate prediction methods are not implemented; 'PredictMethod' must be 'exact'.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'PredictMethod', 'bcd')
%!error<RegressionGP: 'fmincon' is not available in core Octave; use 'quasinewton' or 'fminsearch'.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'Optimizer', 'fmincon')
%!error<RegressionGP: unsupported 'Optimizer' value.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'Optimizer', 'bogus')
%!error<RegressionGP: 'Standardize' must be a logical scalar.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'Standardize', 5)
%!error<RegressionGP: 'PredictorNames' must be a cell array of character vectors with one name per column of X.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'PredictorNames', {'a'})
%!error<RegressionGP: 'ResponseName' must be a character vector.> ...
%! RegressionGP (ones (5, 2), ones (5, 1), 'ResponseName', 5)

## Test input validation for the predict method
%!error<RegressionGP.predict: too few input arguments.> ...
%! predict (RegressionGP (ones (5, 2), ones (5, 1)))
%!error<RegressionGP.predict: XC is empty.> ...
%! predict (RegressionGP (ones (5, 2), ones (5, 1)), [])
%!error<RegressionGP.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (RegressionGP (ones (5, 2), ones (5, 1)), ones (3, 3))
%!error<RegressionGP.predict: 'Alpha' must be a scalar between 0 and 1.> ...
%! predict (RegressionGP (ones (5, 2), ones (5, 1)), ones (3, 2), 'Alpha', 2)
%!error<RegressionGP.predict: invalid NAME in optional pairs of arguments.> ...
%! predict (RegressionGP (ones (5, 2), ones (5, 1)), ones (3, 2), 'bogus', 1)

## Test input validation for the loss method
%!error<RegressionGP.loss: too few input arguments.> ...
%! loss (RegressionGP (ones (5, 2), ones (5, 1)), ones (3, 2))
%!error<RegressionGP.loss: unsupported 'LossFun' value.> ...
%! loss (RegressionGP (ones (5, 2), ones (5, 1)), ones (3, 2), ...
%!       ones (3, 1), 'LossFun', 'bogus')
%!error<RegressionGP.loss: invalid NAME in optional pairs of arguments.> ...
%! loss (RegressionGP (ones (5, 2), ones (5, 1)), ones (3, 2), ...
%!       ones (3, 1), 'bogus', 1)

## Test input validation for the crossval method
%!error<RegressionGP.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (RegressionGP (ones (10, 2), ones (10, 1)), 'KFold', 1)
%!error<RegressionGP.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (RegressionGP (ones (10, 2), ones (10, 1)), 'Holdout', 2)
%!error<RegressionGP.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (RegressionGP (ones (10, 2), ones (10, 1)), 'Leaveout', 1)
%!error<RegressionGP.crossval: 'CVPartition' must be a cvpartition object.> ...
%! crossval (RegressionGP (ones (10, 2), ones (10, 1)), 'CVPartition', 1)
%!error<RegressionGP.crossval: you can use only one of 'KFold', 'Holdout', 'Leaveout', or 'CVPartition' options.> ...
%! crossval (RegressionGP (ones (10, 2), ones (10, 1)), 'KFold', 3, ...
%!           'Holdout', 0.2)
%!error<RegressionGP.crossval: invalid parameter name in optional paired arguments.> ...
%! crossval (RegressionGP (ones (10, 2), ones (10, 1)), 'bogus', 1)
