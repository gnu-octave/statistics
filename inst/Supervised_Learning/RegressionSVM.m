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
## @deftypefn  {statistics} {@var{obj} =} RegressionSVM (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} RegressionSVM (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{RegressionSVM} object containing a support vector machine
## regression model.
##
## @code{@var{obj} = RegressionSVM (@var{X}, @var{Y})} returns a support vector
## regression model, @var{obj}, with @var{X} being the predictor data and
## @var{Y} the continuous response of the observations in @var{X}.
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
## The model is fitted by @math{epsilon}-insensitive regression: errors smaller
## than @qcode{Epsilon} cost nothing, so only the observations outside that
## tube become support vectors.  @qcode{Epsilon} defaults to
## @code{iqr (@var{Y}) / 13.49}, a robust estimate of a tenth of the response's
## standard deviation, which is what MATLAB uses.
##
## @code{@var{obj} = RegressionSVM (@dots{}, @var{name}, @var{value})} returns
## a model with additional options specified by @qcode{Name-Value} pair
## arguments listed below.
##
## @multitable @columnfractions 0.32 0.68
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Standardize'} @tab A logical scalar specifying whether the
## predictor data should be centred and scaled before training.  The same
## transformation is applied by @code{predict}.  The default is @qcode{false}.
##
## @item @qcode{'PredictorNames'} @tab A cell array of character vectors
## naming the predictors, in the order they appear in @var{X}.
##
## @item @qcode{'ResponseName'} @tab A character vector naming the response.
## The default is @qcode{'Y'}.
##
## @item @qcode{'ResponseTransform'} @tab A character vector naming one of the
## supported transformations, or a function handle, applied to the predicted
## response by @code{predict} and @code{resubPredict}.  The default is
## @qcode{'none'}.
##
## @item @qcode{'Epsilon'} @tab A non-negative scalar, the half-width of the
## insensitive tube.  The default is @code{iqr (@var{Y}) / 13.49}, or
## @math{0.1} where that is zero.
##
## @item @qcode{'BoxConstraint'} @tab A positive scalar bounding the dual
## coefficients, the cost of an error outside the tube.  The default is 1.
##
## @item @qcode{'KernelFunction'} @tab A character vector naming the kernel,
## one of @qcode{'linear'}, the default, @qcode{'rbf'}, @qcode{'gaussian'},
## @qcode{'polynomial'} or @qcode{'sigmoid'}.
##
## @item @qcode{'PolynomialOrder'} @tab A positive integer, the order of the
## polynomial kernel.  The default is 3.  It is ignored by every other kernel.
##
## @item @qcode{'KernelScale'} @tab A positive scalar dividing the predictors
## before the kernel is applied.  The default is 1.
##
## @item @qcode{'KernelOffset'} @tab A non-negative scalar added to the kernel
## value.  The default is 0.
##
## @item @qcode{'SVMtype'} @tab A character vector selecting the formulation,
## either @qcode{'eps_svr'}, the default, or @qcode{'nu_svr'}.  MATLAB fits
## only the @math{epsilon} form; @qcode{'nu_svr'} is an Octave extension, in
## which @qcode{Nu} bounds the fraction of support vectors and @qcode{Epsilon}
## is determined by the fit rather than given.
##
## @item @qcode{'Nu'} @tab A scalar in @math{(0, 1]} used by
## @qcode{'nu_svr'}.  The default is 0.5.
##
## @item @qcode{'CacheSize'} @tab A positive scalar, the kernel cache in
## megabytes.  The default is 1000.
##
## @item @qcode{'Tolerance'} @tab A non-negative scalar, the tolerance of the
## termination criterion.  The default is @math{1e-6}.
##
## @item @qcode{'Shrinking'} @tab Either 0 or 1, whether to use the shrinking
## heuristic.  The default is 1.
## @end multitable
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
## @seealso{fitrsvm, ClassificationSVM, RegressionNeuralNetwork}
## @end deftypefn

classdef RegressionSVM

  properties(Access = public)

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} X
    ##
    ## Predictor data
    ##
    ## An @math{NxP} numeric matrix, as it was supplied to the constructor.
    ## This property is read-only.
    ##
    ## @end deftp
    X                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} Y
    ##
    ## Response data
    ##
    ## An @math{Nx1} numeric vector, as it was supplied to the constructor.
    ## This property is read-only.
    ##
    ## @end deftp
    Y                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} NumObservations
    ##
    ## Number of observations used to train the model
    ##
    ## A positive integer scalar, counting only the rows that survived the
    ## removal of missing values.  This property is read-only.
    ##
    ## @end deftp
    NumObservations       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} RowsUsed
    ##
    ## Rows of the original data used for training
    ##
    ## A logical column vector with one entry per row of @code{X}, true where
    ## the row was used.  This property is read-only.
    ##
    ## @end deftp
    RowsUsed              = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer scalar.  This property is read-only.
    ##
    ## @end deftp
    NumPredictors         = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    ResponseName          = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} ResponseTransform
    ##
    ## Transformation applied to the predicted response
    ##
    ## A function handle, applied by @code{predict} and @code{resubPredict} to
    ## the model's output.  It defaults to the identity and may be set after
    ## construction, either to a handle or to the name of a supported
    ## transformation.
    ##
    ## @end deftp
    ResponseTransform     = @(y) y;

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} Epsilon
    ##
    ## Half-width of the insensitive tube
    ##
    ## A non-negative scalar.  An error smaller than @code{Epsilon} costs
    ## nothing, so only observations outside the tube become support vectors.
    ## This property is read-only.
    ##
    ## @end deftp
    Epsilon               = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} Standardize
    ##
    ## Whether the predictor data was standardized
    ##
    ## A logical scalar.  When true the training data was centred and scaled
    ## by @code{Mu} and @code{Sigma}, and @code{predict} applies the same
    ## transformation.  This property is read-only.
    ##
    ## @end deftp
    Standardize           = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} Sigma
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
    ## @deftp {RegressionSVM} {property} Mu
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
    ## @deftp {RegressionSVM} {property} ModelParameters
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
    ## @deftp {RegressionSVM} {property} Model
    ##
    ## The trained LIBSVM model
    ##
    ## A structure as returned by @code{svmtrain} and consumed by
    ## @code{svmpredict}.  This property is read-only.
    ##
    ## @end deftp
    Model                 = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} Alpha
    ##
    ## Dual coefficients of the support vectors
    ##
    ## A numeric column vector with one entry per support vector, holding the
    ## difference of the two multipliers each observation carries.  Unlike a
    ## classifier's, these are signed: there are no labels to take the sign
    ## into, so an observation above the tube and one below it are told apart
    ## by the sign of its coefficient.  This property is read-only.
    ##
    ## @end deftp
    Alpha                 = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} Beta
    ##
    ## Primal coefficients, one per predictor
    ##
    ## A numeric column vector, equal to
    ## @code{obj.SupportVectors' * obj.Alpha}.  It exists only for a linear
    ## kernel; for any other kernel there is no primal representation and this
    ## is empty.  This property is read-only.
    ##
    ## @end deftp
    Beta                  = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} Bias
    ##
    ## Intercept of the fitted function
    ##
    ## A numeric scalar.  With a linear kernel the prediction is
    ## @code{X * obj.Beta + obj.Bias}.  This property is read-only.
    ##
    ## @end deftp
    Bias                  = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} IsSupportVector
    ##
    ## Which training observations are support vectors
    ##
    ## A logical column vector with one entry per training observation.  This
    ## property is read-only.
    ##
    ## @end deftp
    IsSupportVector       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} SupportVectors
    ##
    ## The support vectors themselves
    ##
    ## A numeric matrix with one row per support vector, on the scale the
    ## model was trained on, standardized where @code{Standardize} is true.
    ## This property is read-only.
    ##
    ## @end deftp
    SupportVectors        = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector of column indices, and empty when none is.  This
    ## property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the model expanded them
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionSVM} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector with one entry per training observation.  This
    ## property is read-only.
    ##
    ## @end deftp
    W                     = [];
  endproperties

  properties(Access = private, Hidden)
    RTname = 'none';
  endproperties

  methods(Hidden)

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
      fprintf ("\n  RegressionSVM\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
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
      if (this.Standardize)
        fprintf ("%+25s: [1x%d double]\n", 'Mu', numel (this.Mu));
        fprintf ("%+25s: [1x%d double]\n", 'Sigma', numel (this.Sigma));
      endif
    endfunction

    ## Class specific subscripted reference
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("Invalid () indexing for referencing values", ...
                         " in a RegressionSVM object."));
        case '{}'
          error (strcat ("Invalid {} indexing for referencing values", ...
                         " in a RegressionSVM object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("RegressionSVM.subsref: '.' indexing", ...
                           " argument must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("RegressionSVM.subsref: unrecognized", ...
                           " property: '%s'"), s.subs);
          end_try_catch
      endswitch
      ## Chained references
      if (! isempty (chain_s))
        out = subsref (out, chain_s);
      endif
      varargout{1} = out;
    endfunction

    ## Class specific subscripted assignment
    function this = subsasgn (this, s, val)
      if (numel (s) > 1)
        error ("RegressionSVM.subsasgn: chained subscripts not allowed.");
      endif
      switch s.type
        case '()'
          error (strcat ("Invalid () indexing for assigning values", ...
                         " to a RegressionSVM object."));
        case '{}'
          error (strcat ("Invalid {} indexing for assigning values", ...
                         " to a RegressionSVM object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("RegressionSVM.subsasgn: '.' indexing", ...
                           " argument must be a character vector."));
          endif
          switch (s.subs)
            case 'ResponseTransform'
              [this.ResponseTransform, this.RTname] = ...
                                        parseResponseTransform___ (val);
            otherwise
              error (strcat ("RegressionSVM.subsasgn: unrecognized or", ...
                             " read-only property: '%s'"), s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods(Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionSVM} {@var{obj} =} RegressionSVM (@var{X}, @var{Y})
    ## @deftypefnx {RegressionSVM} {@var{obj} =} RegressionSVM (@dots{}, @var{name}, @var{value})
    ##
    ## Create a @qcode{RegressionSVM} object containing a support vector
    ## machine regression model.
    ##
    ## See the class documentation for the accepted @qcode{Name-Value} pairs.
    ##
    ## @seealso{fitrsvm, RegressionSVM}
    ## @end deftypefn
    function this = RegressionSVM (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("RegressionSVM: too few input arguments.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error ("RegressionSVM: number of rows in X and Y must be equal.");
      endif

      ## The response is continuous, so it must be numeric
      if (! (isnumeric (Y) && isreal (Y)))
        error ("RegressionSVM: Y must be a real numeric vector.");
      endif
      if (! (isvector (Y) || isempty (Y)))
        error ("RegressionSVM: Y must be a vector.");
      endif

      ## Assign original X and Y data to the RegressionSVM object
      this.X = X;
      this.Y = Y;

      ## Set default values before parsing optional parameters
      SVMtype                 = 'eps_svr';
      KernelFunction          = 'linear';
      KernelScale             = 1;
      KernelOffset            = 0;
      PolynomialOrder         = 3;
      BoxConstraint           = 1;
      Epsilon                 = [];
      Nu                      = 0.5;
      CacheSize               = 1000;
      Tolerance               = 1e-6;
      Shrinking               = 1;
      Standardize             = false;
      ResponseName            = [];
      PredictorNames          = [];

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'standardize'
            Standardize = varargin{2};
            if (! (Standardize == true || Standardize == false))
              error (strcat ("RegressionSVM: 'Standardize' must", ...
                             " be either true or false."));
            endif

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("RegressionSVM: 'PredictorNames' must", ...
                             " be supplied as a cellstring array."));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat ("RegressionSVM: 'PredictorNames' must", ...
                             " have the same number of columns as X."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("RegressionSVM: 'ResponseName' must", ...
                             " be a character vector."));
            endif

          case 'responsetransform'
            [this.ResponseTransform, this.RTname] = ...
                                      parseResponseTransform___ (varargin{2});

          case 'svmtype'
            SVMtype = varargin{2};
            if (! (ischar (SVMtype) && isrow (SVMtype)))
              error ("RegressionSVM: 'SVMtype' must be a character vector.");
            endif
            SVMtype = tolower (SVMtype);
            if (! any (strcmp (SVMtype, {'eps_svr', 'nu_svr'})))
              error ("RegressionSVM: unsupported 'SVMtype'.");
            endif

          case 'epsilon'
            Epsilon = varargin{2};
            if (! (isnumeric (Epsilon) && isscalar (Epsilon) && Epsilon >= 0))
              error (strcat ("RegressionSVM: 'Epsilon' must be a", ...
                             " non-negative scalar."));
            endif

          case 'kernelfunction'
            KernelFunction = varargin{2};
            if (! ischar (KernelFunction))
              error (strcat ("RegressionSVM: 'KernelFunction' must", ...
                             " be a character vector."));
            endif
            KernelFunction = tolower (KernelFunction);
            if (! any (strcmpi (KernelFunction, ...
                       {'linear', 'rbf', 'gaussian', 'polynomial', 'sigmoid'})))
              error ("RegressionSVM: unsupported Kernel function.");
            endif

          case 'polynomialorder'
            PolynomialOrder = varargin{2};
            if (! (isnumeric (PolynomialOrder) && isscalar (PolynomialOrder)
                   && PolynomialOrder > 0 && mod (PolynomialOrder, 1) == 0))
              error (strcat ("RegressionSVM: 'PolynomialOrder' must", ...
                             " be a positive integer."));
            endif

          case 'kernelscale'
            KernelScale = varargin{2};
            if (! (isscalar (KernelScale) && KernelScale > 0))
              error (strcat ("RegressionSVM: 'KernelScale' must", ...
                             " be a positive scalar."));
            endif

          case 'kerneloffset'
            KernelOffset = varargin{2};
            if (! (isnumeric (KernelOffset) && isscalar (KernelOffset)
                                            && KernelOffset >= 0))
              error (strcat ("RegressionSVM: 'KernelOffset' must", ...
                             " be a non-negative scalar."));
            endif

          case 'boxconstraint'
            BoxConstraint = varargin{2};
            if (! (isscalar (BoxConstraint) && BoxConstraint > 0))
              error (strcat ("RegressionSVM: 'BoxConstraint' must", ...
                             " be a positive scalar."));
            endif

          case 'nu'
            Nu = varargin{2};
            if (! (isscalar (Nu) && Nu > 0 && Nu <= 1))
              error (strcat ("RegressionSVM: 'Nu' must be a positive", ...
                             " scalar in the range 0 < Nu <= 1."));
            endif

          case 'cachesize'
            CacheSize = varargin{2};
            if (! (isscalar (CacheSize) && CacheSize > 0))
              error ("RegressionSVM: 'CacheSize' must be a positive scalar.");
            endif

          case 'tolerance'
            Tolerance = varargin{2};
            if (! (isscalar (Tolerance) && Tolerance >= 0))
              error ("RegressionSVM: 'Tolerance' must be a positive scalar.");
            endif

          case 'shrinking'
            Shrinking = varargin{2};
            if (! (ismember (Shrinking, [0, 1]) && isscalar (Shrinking)))
              error ("RegressionSVM: 'Shrinking' must be either 0 or 1.");
            endif

          otherwise
            error (strcat ("RegressionSVM: invalid parameter name", ...
                           " in optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Get number of variables in training data
      ndims_X = columns (X);
      this.NumPredictors = ndims_X;

      ## Generate default predictors and response variable names
      if (isempty (PredictorNames))
        for i = 1:ndims_X
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = 'Y';
      endif
      this.PredictorNames = PredictorNames;
      this.ExpandedPredictorNames = PredictorNames;
      this.ResponseName   = ResponseName;
      this.CategoricalPredictors = [];

      ## Remove missing values from X and Y
      RowsUsed  = ! logical (sum (isnan ([X, Y(:)]), 2));
      Y         = Y(RowsUsed);
      X         = X(RowsUsed, :);

      ## Check X and Y contain valid data
      if (! (isnumeric (X) && isfinite (X)))
        error ("RegressionSVM: invalid values in X.");
      endif
      if (isempty (Y))
        error ("RegressionSVM: Y cannot be empty.");
      endif
      if (! all (isfinite (Y)))
        error ("RegressionSVM: invalid values in Y.");
      endif

      this.NumObservations = sum (RowsUsed);
      this.RowsUsed = RowsUsed;
      this.W = ones (this.NumObservations, 1);

      ## Handle Standardize flag.  The model must be fitted on the scale it
      ## predicts on: predict and resubPredict standardize their input from
      ## Mu and Sigma, so the training data is standardized here as well.
      if (Standardize)
        this.Standardize = true;
        this.Sigma = std (X, [], 1);
        this.Sigma(this.Sigma == 0) = 1;  # predictor is constant
        this.Mu = mean (X, 1);
        X = (X - this.Mu) ./ this.Sigma;
      else
        this.Standardize = false;
        this.Sigma = [];
        this.Mu = [];
      endif

      ## Epsilon defaults to a robust tenth of the response's spread, which
      ## is what MATLAB uses.  A constant response gives a zero interquartile
      ## range, and a zero-width tube would make every observation a support
      ## vector, so it falls back to LIBSVM's own default.
      if (isempty (Epsilon))
        Epsilon = iqr (Y) / 13.49;
        if (Epsilon == 0)
          Epsilon = 0.1;
        endif
      endif
      this.Epsilon = Epsilon;

      ## Set svmtrain parameters for SVMtype and KernelFunction
      switch (SVMtype)
        case 'eps_svr'
          s = 3;
        case 'nu_svr'
          s = 4;
      endswitch
      switch (KernelFunction)
        case 'linear'
          t = 0;
        case 'polynomial'
          t = 1;
        case {'rbf', 'gaussian'}
          t = 2;
        case 'sigmoid'
          t = 3;
      endswitch

      ## Set svmtrain parameters for gamma
      g = KernelScale / ndims_X;

      ## Build options string for svmtrain function
      str_options = strcat ("-s %d -t %d -g %.16g -d %d -r %.16g", ...
                            " -c %.16g -n %.16g -p %.16g -m %.16g", ...
                            " -e %e -h %d -q");
      svm_options = sprintf (str_options, s, t, g, PolynomialOrder, ...
                             KernelOffset, BoxConstraint, Nu, Epsilon, ...
                             CacheSize, Tolerance, Shrinking);

      ## Train the SVM model using svmtrain from libsvm
      Model = svmtrain (Y, X, svm_options);
      this.Model = Model;

      ## Populate the model properties.  For regression LIBSVM's sv_coef is
      ## already the difference of the two multipliers, so it is signed and
      ## kept as it stands: unlike a classifier there are no labels to carry
      ## the sign into.
      this.Alpha = Model.sv_coef;

      ## LIBSVM evaluates sum_i coef_i * K(sv_i, x) - rho, so the intercept
      ## MATLAB reports is the negated rho.  Measured against svmpredict.
      this.Bias = -Model.rho;

      ## BETA holds the primal coefficients, one per predictor, and exists
      ## only for a linear kernel; for any other there is no primal
      ## representation and MATLAB leaves it empty.
      if (t == 0)
        this.Beta = Model.SVs' * this.Alpha;
      else
        this.Beta = [];
      endif

      this.IsSupportVector = false (this.NumObservations, 1);
      this.IsSupportVector(Model.sv_indices) = true;
      this.SupportVectors = Model.SVs;

      ## Populate ModelParameters structure
      this.ModelParameters = struct ('SVMtype', SVMtype, 'BoxConstraint', ...
                                     BoxConstraint, 'CacheSize', CacheSize, ...
                                     'KernelScale', KernelScale, ...
                                     'KernelOffset', KernelOffset, ...
                                     'KernelFunction', KernelFunction, ...
                                     'PolynomialOrder', PolynomialOrder, ...
                                     'Epsilon', Epsilon, 'Nu', Nu, ...
                                     'Tolerance', Tolerance, ...
                                     'Shrinking', Shrinking);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionSVM} {@var{yFit} =} predict (@var{obj}, @var{XC})
    ##
    ## Predict the response for new data with a support vector regression
    ## model.
    ##
    ## @code{@var{yFit} = predict (@var{obj}, @var{XC})} returns a column
    ## vector holding the predicted response for each row of @var{XC}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionSVM} class object.
    ## @item
    ## @var{XC} must be a numeric matrix with the same number of predictors as
    ## the data the model was trained on.
    ## @end itemize
    ##
    ## The transformation named by @code{ResponseTransform} is applied to the
    ## model's output before it is returned.
    ##
    ## @seealso{RegressionSVM, fitrsvm}
    ## @end deftypefn
    function yFit = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("RegressionSVM.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("RegressionSVM.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("RegressionSVM.predict: XC must have the same", ...
                       " number of predictors as the trained model."));
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
    ## @deftypefn {RegressionSVM} {@var{yFit} =} resubPredict (@var{obj})
    ##
    ## Predict the response of the training data with a support vector
    ## regression model.
    ##
    ## @code{@var{yFit} = resubPredict (@var{obj})} returns a column vector
    ## holding the predicted response for every observation the model was
    ## trained on.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionSVM} class object.
    ## @end itemize
    ##
    ## @seealso{RegressionSVM, fitrsvm}
    ## @end deftypefn
    function yFit = resubPredict (this)
      yFit = predict (this, this.X(this.RowsUsed, :));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionSVM} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionSVM} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the regression loss of a support vector machine model.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## weighted mean squared error between the response @var{Y} and the
    ## response the model predicts for @var{X}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionSVM} class object.
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
    ## The @math{epsilon}-insensitive loss charges nothing for an error inside
    ## the tube, @code{max (0, abs (@var{Y} - @var{yFit}) - Epsilon)}, which is
    ## the quantity the fit itself minimizes.
    ##
    ## @item @qcode{'Weights'} @tab A numeric vector of observation weights
    ## with one entry per row of @var{X}.  It defaults to a uniform weight.
    ## The weights are normalized to sum to one before the loss is formed.
    ## @end multitable
    ##
    ## @seealso{RegressionSVM, fitrsvm}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("RegressionSVM.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionSVM.loss: Name-Value arguments", ...
                       " must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, 'loss');

      ## Defaults, then the optional pairs
      LossFun = 'mse';
      args = varargin;
      keep = true (1, numel (args));
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("RegressionSVM.loss: parameter name must be", ...
                         " a character vector."));
        endif
        if (strcmpi (args{i}, 'lossfun'))
          LossFun = args{i+1};
          if (! (is_function_handle (LossFun) ||
                 (ischar (LossFun) && isrow (LossFun))))
            error (strcat ("RegressionSVM.loss: 'LossFun' must be a", ...
                           " character vector or a function handle."));
          endif
          if (ischar (LossFun) && ! any (strcmpi (LossFun, ...
                                         {'mse', 'epsiloninsensitive'})))
            error ("RegressionSVM.loss: unsupported 'LossFun' value.");
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
          error (strcat ("RegressionSVM.loss: 'LossFun' must return", ...
                         " a numeric scalar."));
        endif
      elseif (strcmpi (LossFun, 'epsiloninsensitive'))
        L = sum (W .* max (0, abs (Y - yFit) - this.Epsilon));
      else
        L = sum (W .* (Y - yFit) .^ 2);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionSVM} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {RegressionSVM} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the resubstitution regression loss of a support vector machine
    ## model.
    ##
    ## @code{@var{L} = resubLoss (@var{obj})} returns the weighted mean
    ## squared error of the model on the data it was trained on.  It accepts
    ## the same @qcode{Name-Value} pairs as @code{loss}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionSVM} class object.
    ## @end itemize
    ##
    ## @seealso{RegressionSVM, fitrsvm}
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      X = this.X(this.RowsUsed, :);
      Y = this.Y(this.RowsUsed);
      L = loss (this, X, Y, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionSVM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a support vector regression model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves every property of
    ## the @qcode{RegressionSVM} object @var{obj} into @var{filename} in
    ## binary format, so that it can be read back with @code{loadmodel}.
    ##
    ## @seealso{loadmodel, RegressionSVM, fitrsvm}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error ("RegressionSVM.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("RegressionSVM.savemodel: FNAME must be", ...
                       " a character vector."));
      endif
      ## Generate variable for class name
      classdef_name = 'RegressionSVM';

      ## Create variables from model properties
      X                       = this.X;
      Y                       = this.Y;
      NumObservations         = this.NumObservations;
      RowsUsed                = this.RowsUsed;
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
      IsSupportVector         = this.IsSupportVector;
      SupportVectors          = this.SupportVectors;
      CategoricalPredictors   = this.CategoricalPredictors;
      ExpandedPredictorNames  = this.ExpandedPredictorNames;
      W                       = this.W;
      RTname                  = this.RTname;

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'X', 'Y', ...
            'NumObservations', 'RowsUsed', 'NumPredictors', ...
            'PredictorNames', 'ResponseName', 'ResponseTransform', ...
            'Epsilon', 'Standardize', 'Sigma', 'Mu', 'ModelParameters', ...
            'Model', 'Alpha', 'Beta', 'Bias', 'IsSupportVector', ...
            'SupportVectors', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'W', 'RTname');
    endfunction

  endmethods

  methods (Access = private)

    ## Shared validation for the assessment methods, so each reports under
    ## its own name.
    function [X, Y] = checkXY_ (this, X, Y, caller)
      if (isempty (X))
        error ("RegressionSVM.%s: X is empty.", caller);
      elseif (this.NumPredictors != columns (X))
        error (strcat ("RegressionSVM.%s: X must have the same number", ...
                       " of predictors as the trained model."), caller);
      endif
      if (isempty (Y))
        error ("RegressionSVM.%s: Y is empty.", caller);
      elseif (! (isnumeric (Y) && isreal (Y)))
        error (strcat ("RegressionSVM.%s: Y must be a real numeric", ...
                       " vector."), caller);
      elseif (rows (X) != numel (Y))
        error (strcat ("RegressionSVM.%s: Y must have the same number", ...
                       " of rows as X."), caller);
      endif
    endfunction

    ## Pull a "Weights" pair out of the optional arguments, defaulting to a
    ## uniform weight, and reject any other name.
    function W = getWeights_ (this, args, n, caller)
      W = ones (n, 1);
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("RegressionSVM.%s: parameter name must be", ...
                         " a character vector."), caller);
        endif
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
          if (! (isnumeric (W) && isvector (W)))
            error (strcat ("RegressionSVM.%s: 'Weights' must be a", ...
                           " numeric vector."), caller);
          endif
          if (numel (W) != n)
            error (strcat ("RegressionSVM.%s: size of 'Weights' must", ...
                           " equal the number of rows in X."), caller);
          endif
        else
          error (strcat ("RegressionSVM.%s: invalid parameter name in", ...
                         " optional paired arguments."), caller);
        endif
      endfor
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a RegressionSVM object
      mdl = RegressionSVM ([1; 2; 3], [1; 2; 3]);

      ## Get fieldnames from DATA (including private properties)
      names = fieldnames (data);

      ## Copy data into object
      for i = 1:numel (names)
        ## Check fieldnames in DATA match the class properties
        try
          mdl.(names{i}) = data.(names{i});
        catch
          error ("RegressionSVM.load_model: invalid model in '%s'.", filename)
        end_try_catch
      endfor
    endfunction

  endmethods

endclassdef

## Resolve a ResponseTransform, given either as a name or as a handle, into
## the handle that is applied and the name that is displayed and saved.
function [rt, name] = parseResponseTransform___ (val)
  if (is_function_handle (val))
    v = (1:5)';
    if (! isequal (size (v), size (val (v))))
      error (strcat ("RegressionSVM: function handle for", ...
                     " 'ResponseTransform' must return the same size", ...
                     " as its input."));
    endif
    rt = val;
    name = 'custom function handle';
  elseif (ischar (val) && isrow (val))
    name = tolower (val);
    switch (name)
      case {'none', 'identity'}
        rt = @(y) y;
      case 'exp'
        rt = @(y) exp (y);
      case 'log'
        rt = @(y) log (y);
      otherwise
        error (strcat ("RegressionSVM: unrecognized", ...
                       " 'ResponseTransform' function."));
    endswitch
  else
    error (strcat ("RegressionSVM: 'ResponseTransform' must be a", ...
                   " character vector or a function handle."));
  endif
endfunction

## A fitted model carries the defaults MATLAB documents for fitrsvm.
%!test
%! X = [linspace(0, 1, 40)', linspace(2, 3, 40)'];
%! Y = 2 * X(:,1) + 0.5;
%! Mdl = RegressionSVM (X, Y);
%! assert_equal (class (Mdl), 'RegressionSVM');
%! assert_equal (Mdl.ModelParameters.KernelFunction, 'linear');
%! assert_equal (Mdl.ModelParameters.SVMtype, 'eps_svr');
%! assert_equal (Mdl.ModelParameters.BoxConstraint, 1);
%! assert_equal (Mdl.ModelParameters.KernelScale, 1);
%! assert_equal (Mdl.ModelParameters.PolynomialOrder, 3);
%! assert_equal (Mdl.Standardize, false);
%! assert_equal (Mdl.NumObservations, 40);
%! assert_equal (Mdl.NumPredictors, 2);
%! assert_equal (Mdl.ResponseName, 'Y');
%! assert_equal (Mdl.PredictorNames, {'x1', 'x2'});

## Epsilon defaults to a robust tenth of the response's spread.
%!test
%! randn ('seed', 42);
%! X = randn (60, 2);
%! Y = X(:,1) * 3 + randn (60, 1);
%! Mdl = RegressionSVM (X, Y);
%! assert_equal (Mdl.Epsilon, iqr (Y) / 13.49, 1e-12);
%! M2 = RegressionSVM (X, Y, 'Epsilon', 0.25);
%! assert_equal (M2.Epsilon, 0.25);

## A constant response has no spread, so Epsilon falls back rather than
## collapsing the tube to zero width.
%!test
%! X = [linspace(0, 1, 20)', linspace(1, 2, 20)'];
%! Mdl = RegressionSVM (X, ones (20, 1));
%! assert_equal (Mdl.Epsilon, 0.1);

## With a linear kernel the model is a plain linear function of its input.
%!test
%! randn ('seed', 42);
%! X = randn (50, 3);
%! Y = X * [2; -1; 0.5] + 1;
%! Mdl = RegressionSVM (X, Y);
%! assert_equal (size (Mdl.Beta), [3, 1]);
%! assert_equal (isscalar (Mdl.Bias), true);
%! assert_equal (X * Mdl.Beta + Mdl.Bias, resubPredict (Mdl), 1e-8);
%! assert_equal (Mdl.Beta, Mdl.SupportVectors' * Mdl.Alpha, 1e-12);

## A non-linear kernel has no primal representation, so Beta is empty.
%!test
%! randn ('seed', 42);
%! X = randn (40, 2);
%! Y = sum (X .^ 2, 2);
%! Mdl = RegressionSVM (X, Y, 'KernelFunction', 'rbf');
%! assert_equal (isempty (Mdl.Beta), true);
%! assert_equal (numel (resubPredict (Mdl)), 40);

## Alpha is signed, one entry per support vector, and IsSupportVector marks
## the training rows they came from.
%!test
%! randn ('seed', 42);
%! X = randn (60, 2);
%! Y = X(:,1) - X(:,2) + randn (60, 1) * 0.5;
%! Mdl = RegressionSVM (X, Y);
%! nsv = rows (Mdl.SupportVectors);
%! assert_equal (size (Mdl.Alpha), [nsv, 1]);
%! assert_equal (any (Mdl.Alpha < 0), true);
%! assert_equal (class (Mdl.IsSupportVector), 'logical');
%! assert_equal (numel (Mdl.IsSupportVector), 60);
%! assert_equal (sum (Mdl.IsSupportVector), nsv);

## A wider tube is fitted by fewer support vectors.
%!test
%! randn ('seed', 42);
%! X = randn (60, 2);
%! Y = X(:,1) * 2 + randn (60, 1);
%! narrow = RegressionSVM (X, Y, 'Epsilon', 0.05);
%! wide   = RegressionSVM (X, Y, 'Epsilon', 3);
%! assert_equal (sum (wide.IsSupportVector) < ...
%!               sum (narrow.IsSupportVector), true);

## Standardize fits on the scale it predicts on.  A model fitted on raw data
## and asked about standardized data is not merely worse, it is wrong.
%!test
%! randn ('seed', 42);
%! X = [randn(60, 1), randn(60, 1) * 1000];
%! Y = X(:,1) + X(:,2) / 1000;
%! Mdl = RegressionSVM (X, Y, 'Standardize', true, 'Epsilon', 0.01);
%! assert_equal (Mdl.Standardize, true);
%! assert_equal (size (Mdl.Mu), [1, 2]);
%! assert_equal (size (Mdl.Sigma), [1, 2]);
%! assert_equal (predict (Mdl, X), resubPredict (Mdl));
%! assert_equal (sqrt (resubLoss (Mdl)) < std (Y), true);

## A constant predictor gets a unit scale rather than a division by zero.
%!test
%! X = [linspace(0, 1, 20)', ones(20, 1)];
%! Mdl = RegressionSVM (X, X(:,1), 'Standardize', true);
%! assert_equal (Mdl.Sigma(2), 1);
%! assert_equal (all (isfinite (resubPredict (Mdl))), true);

## Rows carrying a missing value are dropped from both X and Y.
%!test
%! X = [linspace(0, 1, 12)', linspace(1, 2, 12)'; NaN, 1; 0.5, 1];
%! Y = [2 * linspace(0, 1, 12)'; 1; NaN];
%! Mdl = RegressionSVM (X, Y);
%! assert_equal (Mdl.NumObservations, 12);
%! assert_equal (sum (Mdl.RowsUsed), 12);
%! assert_equal (Mdl.RowsUsed(13:14), [false; false]);
%! assert_equal (numel (resubPredict (Mdl)), 12);

## predict on the training rows is resubPredict.
%!test
%! randn ('seed', 42);
%! X = randn (40, 2);
%! Y = X(:,1) + 2 * X(:,2);
%! Mdl = RegressionSVM (X, Y);
%! assert_equal (predict (Mdl, X), resubPredict (Mdl));

## loss defaults to the weighted mean squared error, and weights are
## normalized, so scaling every weight leaves the loss alone.
%!test
%! randn ('seed', 42);
%! X = randn (30, 2);
%! Y = X(:,1) * 3 + 1;
%! Mdl = RegressionSVM (X, Y);
%! yFit = predict (Mdl, X);
%! assert_equal (loss (Mdl, X, Y), mean ((Y - yFit) .^ 2), 1e-12);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'mse'), loss (Mdl, X, Y), 1e-12);
%! w = rand (30, 1) + 0.1;
%! assert_equal (loss (Mdl, X, Y, 'Weights', w), ...
%!               loss (Mdl, X, Y, 'Weights', 7 * w), 1e-12);
%! assert_equal (loss (Mdl, X, Y, 'Weights', w), ...
%!               sum ((w / sum (w)) .* (Y - yFit) .^ 2), 1e-12);

## The epsilon-insensitive loss charges nothing inside the tube.
%!test
%! randn ('seed', 42);
%! X = randn (30, 2);
%! Y = X(:,1) * 3 + 1;
%! Mdl = RegressionSVM (X, Y, 'Epsilon', 0.5);
%! yFit = predict (Mdl, X);
%! L = loss (Mdl, X, Y, 'LossFun', 'epsiloninsensitive');
%! assert_equal (L, mean (max (0, abs (Y - yFit) - 0.5)), 1e-12);
%! assert_equal (L <= loss (Mdl, X, Y, 'LossFun', 'mse') + 1, true);
%! assert_equal (L >= 0, true);

## loss takes a function handle of the response, the fit and the weights.
%!test
%! randn ('seed', 42);
%! X = randn (20, 2);
%! Y = X(:,1) + 1;
%! Mdl = RegressionSVM (X, Y);
%! f = @(y, yf, w) sum (w .* abs (y - yf));
%! assert_equal (loss (Mdl, X, Y, 'LossFun', f), ...
%!               mean (abs (Y - predict (Mdl, X))), 1e-12);

## resubLoss is loss on the training data.
%!test
%! randn ('seed', 42);
%! X = [randn(20, 2); NaN, 1];
%! Y = [randn(20, 1); 2];
%! Mdl = RegressionSVM (X, Y);
%! Xu = X(Mdl.RowsUsed, :);
%! Yu = Y(Mdl.RowsUsed);
%! assert_equal (resubLoss (Mdl), loss (Mdl, Xu, Yu), 1e-12);
%! assert_equal (resubLoss (Mdl, 'LossFun', 'epsiloninsensitive'), ...
%!               loss (Mdl, Xu, Yu, 'LossFun', 'epsiloninsensitive'), 1e-12);

## ResponseTransform is applied to the prediction, by name or by handle.
%!test
%! X = [linspace(0, 1, 20)', linspace(1, 2, 20)'];
%! Y = 2 * X(:,1) + 1;
%! Mdl = RegressionSVM (X, Y);
%! raw = predict (Mdl, X);
%! Mdl.ResponseTransform = 'exp';
%! assert_equal (predict (Mdl, X), exp (raw), 1e-12);
%! Mdl.ResponseTransform = @(y) 2 * y;
%! assert_equal (predict (Mdl, X), 2 * raw, 1e-12);
%! Mdl.ResponseTransform = 'none';
%! assert_equal (predict (Mdl, X), raw, 1e-12);

## nu-SVR is an Octave extension the LIBSVM engine already provides.
%!test
%! randn ('seed', 42);
%! X = randn (50, 2);
%! Y = X(:,1) - X(:,2);
%! Mdl = RegressionSVM (X, Y, 'SVMtype', 'nu_svr', 'Nu', 0.3);
%! assert_equal (Mdl.ModelParameters.SVMtype, 'nu_svr');
%! assert_equal (Mdl.ModelParameters.Nu, 0.3);
%! assert_equal (Mdl.Model.Parameters(1), 4);
%! assert_equal (all (isfinite (resubPredict (Mdl))), true);

## Every kernel trains and predicts finitely.
%!test
%! randn ('seed', 42);
%! X = randn (40, 2);
%! Y = X(:,1) + X(:,2);
%! names = {'linear', 'rbf', 'gaussian', 'polynomial', 'sigmoid'};
%! for k = 1:numel (names)
%!   Mdl = RegressionSVM (X, Y, 'KernelFunction', names{k});
%!   assert_equal (all (isfinite (resubPredict (Mdl))), true);
%! endfor

## A saved model comes back carrying its own numbers.
%!test
%! randn ('seed', 42);
%! X = randn (30, 2);
%! Y = X(:,1) * 4 - 1;
%! Mdl = RegressionSVM (X, Y, 'Standardize', true);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'RegressionSVM');
%! assert_equal (M2.Alpha, Mdl.Alpha);
%! assert_equal (M2.Beta, Mdl.Beta);
%! assert_equal (M2.Bias, Mdl.Bias);
%! assert_equal (M2.Epsilon, Mdl.Epsilon);
%! assert_equal (M2.SupportVectors, Mdl.SupportVectors);
%! assert_equal (M2.NumObservations, Mdl.NumObservations);
%! assert_equal (predict (M2, X), predict (Mdl, X));

## Test input validation for the constructor
%!error<RegressionSVM: too few input arguments.> ...
%! RegressionSVM ()
%!error<RegressionSVM: too few input arguments.> ...
%! RegressionSVM (ones (10, 2))
%!error<RegressionSVM: number of rows in X and Y must be equal.> ...
%! RegressionSVM (ones (10, 2), ones (5, 1))
%!error<RegressionSVM: Y must be a real numeric vector.> ...
%! RegressionSVM (ones (5, 2), {'a'; 'b'; 'c'; 'd'; 'e'})
%!error<RegressionSVM: Y must be a vector.> ...
%! RegressionSVM (ones (5, 2), ones (5, 3))
%!error<RegressionSVM: invalid values in X.> ...
%! RegressionSVM ([1, 1; Inf, 1; 3, 1], [1; 2; 3])
%!error<RegressionSVM: invalid values in Y.> ...
%! RegressionSVM ([1, 1; 2, 1; 3, 1], [1; Inf; 3])
%!error<RegressionSVM: 'Standardize' must be either true or false.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'Standardize', 'yes')
%!error<RegressionSVM: 'PredictorNames' must be supplied as a cellstring array.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'PredictorNames', 'a')
%!error<RegressionSVM: 'PredictorNames' must have the same number of columns as X.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'PredictorNames', {'a'})
%!error<RegressionSVM: 'ResponseName' must be a character vector.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'ResponseName', 5)
%!error<RegressionSVM: 'ResponseTransform' must be a character vector or a function handle.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'ResponseTransform', 5)
%!error<RegressionSVM: unrecognized 'ResponseTransform' function.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'ResponseTransform', 'nope')
%!error<RegressionSVM: 'SVMtype' must be a character vector.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'SVMtype', 5)
%!error<RegressionSVM: unsupported 'SVMtype'.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'SVMtype', 'c_svc')
%!error<RegressionSVM: 'Epsilon' must be a non-negative scalar.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'Epsilon', -1)
%!error<RegressionSVM: 'KernelFunction' must be a character vector.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'KernelFunction', 5)
%!error<RegressionSVM: unsupported Kernel function.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'KernelFunction', 'nope')
%!error<RegressionSVM: 'PolynomialOrder' must be a positive integer.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'PolynomialOrder', 2.5)
%!error<RegressionSVM: 'KernelScale' must be a positive scalar.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'KernelScale', 0)
%!error<RegressionSVM: 'KernelOffset' must be a non-negative scalar.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'KernelOffset', -1)
%!error<RegressionSVM: 'BoxConstraint' must be a positive scalar.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'BoxConstraint', 0)
%!error<RegressionSVM: 'Nu' must be a positive scalar in the range 0 < Nu <= 1.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'Nu', 0)
%!error<RegressionSVM: 'CacheSize' must be a positive scalar.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'CacheSize', 0)
%!error<RegressionSVM: 'Tolerance' must be a positive scalar.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'Tolerance', -1)
%!error<RegressionSVM: 'Shrinking' must be either 0 or 1.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'Shrinking', 2)
%!error<RegressionSVM: invalid parameter name in optional pair arguments.> ...
%! RegressionSVM (ones (5, 2), ones (5, 1), 'Prior', 1)

## Test input validation for predict and loss
%!shared RSVM
%! RSVM = RegressionSVM ([1, 1; 2, 1; 3, 2; 4, 2], [2; 4; 6; 8]);
%!error<RegressionSVM.predict: too few input arguments.> ...
%! predict (RSVM)
%!error<RegressionSVM.predict: XC is empty.> ...
%! predict (RSVM, [])
%!error<RegressionSVM.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (RSVM, ones (2, 3))
%!error<RegressionSVM.loss: too few input arguments.> ...
%! loss (RSVM)
%!error<RegressionSVM.loss: Name-Value arguments must be in pairs.> ...
%! loss (RSVM, [1, 1; 2, 1], [2; 4], 'Weights')
%!error<RegressionSVM.loss: X is empty.> ...
%! loss (RSVM, [], [2; 4])
%!error<RegressionSVM.loss: X must have the same number of predictors as the trained model.> ...
%! loss (RSVM, ones (2, 3), [2; 4])
%!error<RegressionSVM.loss: Y is empty.> ...
%! loss (RSVM, [1, 1; 2, 1], [])
%!error<RegressionSVM.loss: Y must be a real numeric vector.> ...
%! loss (RSVM, [1, 1; 2, 1], {'a'; 'b'})
%!error<RegressionSVM.loss: Y must have the same number of rows as X.> ...
%! loss (RSVM, [1, 1; 2, 1], [2; 4; 6])
%!error<RegressionSVM.loss: 'LossFun' must be a character vector or a function handle.> ...
%! loss (RSVM, [1, 1; 2, 1], [2; 4], 'LossFun', 5)
%!error<RegressionSVM.loss: unsupported 'LossFun' value.> ...
%! loss (RSVM, [1, 1; 2, 1], [2; 4], 'LossFun', 'mae')
%!error<RegressionSVM.loss: 'LossFun' must return a numeric scalar.> ...
%! loss (RSVM, [1, 1; 2, 1], [2; 4], 'LossFun', @(y, yf, w) [1, 2])
%!error<RegressionSVM.loss: 'Weights' must be a numeric vector.> ...
%! loss (RSVM, [1, 1; 2, 1], [2; 4], 'Weights', {'a'})
%!error<RegressionSVM.loss: size of 'Weights' must equal the number of rows in X.> ...
%! loss (RSVM, [1, 1; 2, 1], [2; 4], 'Weights', [1; 2; 3])
%!error<RegressionSVM.loss: invalid parameter name in optional paired arguments.> ...
%! loss (RSVM, [1, 1; 2, 1], [2; 4], 'Nope', 1)

## Test input validation for savemodel
%!error<RegressionSVM.savemodel: too few input arguments.> ...
%! savemodel (RSVM)
%!error<RegressionSVM.savemodel: FNAME must be a character vector.> ...
%! savemodel (RSVM, 5)

## Test subscripted reference and assignment
%!error<Invalid \(\) indexing for referencing values in a RegressionSVM object.> ...
%! RSVM(1)
%!error<Invalid \{\} indexing for referencing values in a RegressionSVM object.> ...
%! RSVM{1}
%!error<RegressionSVM.subsref: unrecognized property: 'Nope'> ...
%! RSVM.Nope
%!error<Invalid \(\) indexing for assigning values to a RegressionSVM object.> ...
%! RSVM(1) = 1;
%!error<RegressionSVM.subsasgn: unrecognized or read-only property: 'Epsilon'> ...
%! RSVM.Epsilon = 5;
%!error<RegressionSVM.subsasgn: chained subscripts not allowed.> ...
%! RSVM.ResponseTransform.foo = 1;
%!error<RegressionSVM: unrecognized 'ResponseTransform' function.> ...
%! RSVM.ResponseTransform = 'nope';
