## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2023-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef RegressionGAM
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{obj} =} RegressionGAM (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} RegressionGAM (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{RegressionGAM} class object containing a Generalized Additive
## Model (GAM) for regression.
##
## A @qcode{RegressionGAM} class object can store the predictors and response
## data along with various parameters for the GAM model.  It is recommended to
## use the @code{fitrgam} function to create a @qcode{RegressionGAM} object.
##
## @code{@var{obj} = RegressionGAM (@var{X}, @var{Y})} returns an object of
## class RegressionGAM, with matrix @var{X} containing the predictor data and
## vector @var{Y} containing the continuous response data.
##
## @itemize
## @item
## @var{X} must be a @math{N*P} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the GAM model.
## @item
## @var{Y} must be @math{N*1} numeric vector containing the response data
## corresponding to the predictor data in @var{X}. @var{Y} must have same
## number of rows as @var{X}.
## @end itemize
##
## @code{@var{obj} = RegressionGAM (@dots{}, @var{name}, @var{value})} returns
## an object of class RegressionGAM with additional properties specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.2 0.75
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'predictors'} @tab Predictor Variable names, specified as
## a row vector cell of strings with the same length as the columns in @var{X}.
## If omitted, the program will generate default variable names
## @qcode{(x1, x2, ..., xn)} for each column in @var{X}.
##
## @item @qcode{'responsename'} @tab Response Variable Name, specified as
## a string.  If omitted, the default value is @qcode{'Y'}.
##
## @item @qcode{'formula'} @tab a model specification given as a string in
## the form @qcode{'Y ~ terms'} where @qcode{Y} represents the response variable
## and @qcode{terms} the predictor variables.  The formula can be used to
## specify a subset of variables for training model.  For example:
## @qcode{'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3'} specifies four linear terms
## for the first four columns of for predictor data, and @qcode{x1:x2} and
## @qcode{x2:x3} specify the two interaction terms for 1st-2nd and 3rd-4th
## columns respectively.  Only these terms will be used for training the model,
## but @var{X} must have at least as many columns as referenced in the formula.
## If Predictor Variable names have been defined, then the terms in the formula
## must reference to those.  When @qcode{'formula'} is specified, all terms used
## for training the model are referenced in the @qcode{IntMatrix} field of the
## @var{obj} class object as a matrix containing the column indexes for each
## term including both the predictors and the interactions used.
##
## @item @qcode{'interactions'} @tab a logical matrix, a positive integer
## scalar, or the string @qcode{'all'} for defining the interactions between
## predictor variables.  When given a logical matrix, it must have the same
## number of columns as @var{X} and each row corresponds to a different
## interaction term combining the predictors indexed as @qcode{true}.  Each
## interaction term is appended as a column vector after the available predictor
## column in @var{X}.  When @qcode{'all'} is defined, then all possible
## combinations of interactions are appended in @var{X} before training.  At the
## moment, parsing a positive integer has the same effect as the @qcode{'all'}
## option.  When @qcode{'interactions'} is specified, only the interaction terms
## appended to @var{X} are referenced in the @qcode{IntMatrix} field of the
## @var{obj} class object.
##
## @item @qcode{'knots'} @tab a scalar or a row vector with the same
## columns as @var{X}.  It defines the knots for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 5, hence expanded to @qcode{ones (1, columns (X)) * 5}.  You can
## parse a row vector with different number of knots for each predictor
## variable to be fitted with, although not recommended.
##
## @item @qcode{'order'} @tab a scalar or a row vector with the same
## columns as @var{X}.  It defines the order of the polynomial when training the
## GAM.  As a scalar, it is expanded to a row vector.  The default values is 3,
## hence expanded to @qcode{ones (1, columns (X)) * 3}.  You can parse a row
## vector with different number of polynomial order for each predictor variable
## to be fitted with, although not recommended.
##
## @item @qcode{'dof'} @tab a scalar or a row vector with the same columns
## as @var{X}.  It defines the degrees of freedom for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 8, hence expanded to @qcode{ones (1, columns (X)) * 8}.  You can
## parse a row vector with different degrees of freedom for each predictor
## variable to be fitted with, although not recommended.
##
## @item @qcode{'tol'} @tab a positive scalar to set the tolerance for
## convergence during training. By default, it is set to @qcode{1e-3}.
## @end multitable
##
## You can parse either a @qcode{'formula'} or an @qcode{'interactions'}
## optional parameter.  Parsing both parameters will result an error.
## Accordingly, you can only pass up to two parameters among @qcode{'knots'},
## @qcode{'order'}, and @qcode{'dof'} to define the required polynomial for
## training the GAM model.
##
## @seealso{fitrgam, regress, regress_gp}
## @end deftypefn

  properties(Access = public)

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix with one row per observation and one column per
    ## predictor of the training data.  This property is read-only.
    ##
    ## @end deftp
    X                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} Y
    ##
    ## Response data
    ##
    ## A numeric column vector with one entry per observation of the
    ## training data.  This property is read-only.
    ##
    ## @end deftp
    Y                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer, the number of observations of the training data
    ## the model was fitted on, rows with missing values excluded.  This
    ## property is read-only.
    ##
    ## @end deftp
    NumObservations       = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} RowsUsed
    ##
    ## Rows used for fitting
    ##
    ## A logical column vector with the same length as the observations in the
    ## original predictor data @var{X}, true for each row that was used for
    ## fitting the RegressionGAM model.  It is empty, @qcode{[]},
    ## when every observation was used, so a non-empty value means that rows
    ## holding missing values were dropped.  This property is read-only.
    ##
    ## @end deftp
    RowsUsed              = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer, the number of predictors of the training data.
    ## This property is read-only.
    ##
    ## @end deftp
    NumPredictors         = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} PredictorNames
    ##
    ## Names of the predictor variables
    ##
    ## A cell array of character vectors naming the predictors, in the order
    ## they appear in the training data.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector naming the response variable @var{Y}.  This
    ## property is read-only.
    ##
    ## @end deftp
    ResponseName          = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector holding the column of each predictor treated as
    ## categorical, and empty when none is.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} ExpandedPredictorNames
    ##
    ## Names of the expanded predictor variables
    ##
    ## A cell array of character vectors naming the predictors as the model
    ## sees them.  It matches @code{PredictorNames} unless a categorical
    ## predictor was expanded into dummy variables.  This property is
    ## read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector with one entry per observation used for
    ## training, normalised to sum to one.  This property is read-only.
    ##
    ## @end deftp
    W                     = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} ResponseTransform
    ##
    ## Transformation applied to the predicted response
    ##
    ## A function handle applied to the response the model predicts.  Add or
    ## change it using dot notation, as in
    ## @qcode{@var{obj}.ResponseTransform = 'log'} or
    ## @qcode{@var{obj}.ResponseTransform = @@function_handle}.  It defaults
    ## to @qcode{'none'}, the identity.
    ##
    ## @end deftp
    ResponseTransform     = @(x) x;

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} Intercept
    ##
    ## Intercept of the fitted model
    ##
    ## A numeric scalar, the mean of the response, which every additive term
    ## is measured against.  This property is read-only.
    ##
    ## @end deftp
    Intercept             = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} Formula
    ##
    ## Formula of the model
    ##
    ## A character vector naming the response and the terms of the model, as
    ## in @qcode{'Y ~ x1 + x2 + x1:x2'}, or empty when the model was not
    ## given one.  This property is read-only.
    ##
    ## @end deftp
    Formula               = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} Interactions
    ##
    ## Interaction terms of the model
    ##
    ## A logical matrix, a matrix of column indices, a count of terms, or the
    ## character vector @qcode{'all'}, describing the interaction terms asked
    ## for.  This property is read-only.
    ##
    ## @end deftp
    Interactions          = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} Knots
    ##
    ## Knots of the spline fitting
    ##
    ## A numeric vector with one entry per predictor, the number of breaks
    ## the spline of that predictor is fitted over.  This property is
    ## read-only.
    ##
    ## @end deftp
    Knots                 = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} Order
    ##
    ## Order of the spline fitting
    ##
    ## A numeric vector with one entry per predictor, the polynomial order of
    ## the spline of that predictor.  This property is read-only.
    ##
    ## @end deftp
    Order                 = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} DoF
    ##
    ## Degrees of freedom of the spline fitting
    ##
    ## A numeric vector with one entry per predictor, the sum of its number
    ## of knots and its order.  This property is read-only.
    ##
    ## @end deftp
    DoF                   = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} Tol
    ##
    ## Tolerance for convergence
    ##
    ## A positive scalar, the largest change in the residual sum of squares
    ## of a backfitting cycle that counts as converged.  This property is
    ## read-only.
    ##
    ## @end deftp
    Tol                   = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} IsStandardDeviationFit
    ##
    ## Flag for a fitted standard deviation model
    ##
    ## A boolean flag, always @qcode{false}, as this class estimates the
    ## standard deviation of a prediction from the residuals of the fit
    ## rather than fitting a model for it.  This property is read-only.
    ##
    ## @end deftp
    IsStandardDeviationFit = false;

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} BaseModel
    ##
    ## Model without interaction terms
    ##
    ## A structure holding the intercept, the piecewise polynomial of each
    ## predictor, the number of backfitting cycles, the residuals and the
    ## residual sum of squares of the model fitted without interaction
    ## terms.  This property is read-only.
    ##
    ## @end deftp
    BaseModel             = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} ModelwInt
    ##
    ## Model with interaction terms
    ##
    ## A structure of the same fields as @code{BaseModel}, for the model
    ## fitted with the interaction terms, and empty when none was asked for.
    ## This property is read-only.
    ##
    ## @end deftp
    ModelwInt             = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} IntMatrix
    ##
    ## Interaction terms applied to the predictor data
    ##
    ## A logical matrix or a matrix of column indices, one row per
    ## interaction term, describing the terms appended to the predictor data.
    ## This property is read-only.
    ##
    ## @end deftp
    IntMatrix             = [];

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
      fprintf ("\n  RegressionGAM\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumPredictors', this.NumPredictors);
      fprintf ("%+25s: '%s'\n", 'ResponseTransform', this.RTname);
      fprintf ("%+25s: %g\n", 'Intercept', this.Intercept);
      str = repmat ({'%d'}, 1, numel (this.Knots));
      str = strcat ('[', strjoin (str, ' '), ']');
      fprintf ("%+25s: %s\n", 'Knots', sprintf (str, this.Knots));
      str = repmat ({'%d'}, 1, numel (this.Order));
      str = strcat ('[', strjoin (str, ' '), ']');
      fprintf ("%+25s: %s\n", 'Order', sprintf (str, this.Order));
      fprintf ("%+25s: %g\n", 'Tol', this.Tol);
    endfunction

    ## Class specific subscripted reference
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("Invalid () indexing for referencing values", ...
                         " in a RegressionGAM object."));
        case '{}'
          error (strcat ("Invalid {} indexing for referencing values", ...
                         " in a RegressionGAM object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("RegressionGAM.subsref: '.' indexing argument", ...
                           " must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("RegressionGAM.subsref: unrecognized", ...
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
        error (strcat ("RegressionGAM.subsasgn: chained subscripts", ...
                       " not allowed."));
      endif
      switch s.type
        case '()'
          error (strcat ("Invalid () indexing for assigning values", ...
                         " to a RegressionGAM object."));
        case '{}'
          error (strcat ("Invalid {} indexing for assigning values", ...
                         " to a RegressionGAM object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("RegressionGAM.subsasgn: '.' indexing", ...
                           " argument must be a character vector."));
          endif
          switch (s.subs)
            case 'ResponseTransform'
              [this.ResponseTransform, this.RTname] = ...
                        parseResponseTransform (val, 'RegressionGAM');
            otherwise
              error (strcat ("RegressionGAM.subsasgn: unrecognized or", ...
                             " read-only property: '%s'"), s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods(Access = public)

    ## Class object constructor
    function this = RegressionGAM (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("RegressionGAM: too few input arguments.");
      endif

      ## Get training sample size and number of variables in training data
      nsample = rows (X);
      ndims_X = columns (X);

      ## Check correspondence between predictors and response
      if (nsample != rows (Y))
        error ("RegressionGAM: number of rows in X and Y must be equal.");
      endif

      ## Set default values before parsing optional parameters
      PredictorNames = {};                    # Predictor variable names
      ResponseName   = [];                    # Response variable name
      Formula        = [];                    # Formula for GAM model
      Interactions   = [];                    # Interaction terms
      DoF            = ones (1, ndims_X) * 8; # Degrees of freedom
      Order          = ones (1, ndims_X) * 3; # Order of spline
      Knots          = ones (1, ndims_X) * 5; # Knots
      Tol            = 1e-3;                  # Tolerance for convergence
      ResponseTransform = @(y) y;             # Transform of the prediction
      RTname            = 'none';             # and its name

      ## Number of parameters for Knots, DoF, Order (maximum 2 allowed)
      KOD = 0;
      ## Number of parameters for Formula, Interactions (maximum 1 allowed)
      F_I = 0;

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case {'predictors', 'predictornames'}
            PredictorNames = varargin{2};
            if (! isempty (PredictorNames))
              if (! iscellstr (PredictorNames))
                error (strcat ("RegressionGAM: PredictorNames must", ...
                               " be a cellstring array."));
              elseif (columns (PredictorNames) != columns (X))
                error (strcat ("RegressionGAM: PredictorNames must", ...
                               " have same number of columns as X."));
              endif
            endif

          case 'responsetransform'
            [ResponseTransform, RTname] = ...
                      parseResponseTransform (varargin{2}, 'RegressionGAM');

          case 'responsename'
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error ("RegressionGAM: ResponseName must be a char string.");
            endif

          case 'formula'
            if (F_I < 1)
              Formula = varargin{2};
              if (! ischar (Formula) && ! islogical (Formula))
                error ("RegressionGAM: Formula must be a string.");
              endif
              F_I += 1;
            else
              error ("RegressionGAM: Interactions have been already defined.");
            endif

          case 'interactions'
            if (F_I < 1)
              tmp = varargin{2};
              if (isnumeric (tmp) && isscalar (tmp)
                                  && tmp == fix (tmp) && tmp >= 0)
                Interactions = tmp;
              elseif (islogical (tmp))
                Interactions = tmp;
              elseif (ischar (tmp) && strcmpi (tmp, 'all'))
                Interactions = tmp;
              else
                error ("RegressionGAM: invalid Interactions parameter.");
              endif
              F_I += 1;
            else
              error ("RegressionGAM: Formula has been already defined.");
            endif

          case 'knots'
            if (KOD < 2)
              Knots = varargin{2};
              if (! isnumeric (Knots) || ! (isscalar (Knots) ||
                  isequal (size (Knots), [1, ndims_X])))
                error ("RegressionGAM: invalid value for Knots.");
              endif
              DoF = Knots + Order;
              Order = DoF - Knots;
              KOD += 1;
            else
              error ("RegressionGAM: DoF and Order have been set already.");
            endif

          case 'order'
            if (KOD < 2)
              Order = varargin{2};
              if (! isnumeric (Order) || ! (isscalar (Order) ||
                  isequal (size (Order), [1, ndims_X])))
                error ("RegressionGAM: invalid value for Order.");
              endif
              DoF = Knots + Order;
              Knots = DoF - Order;
              KOD += 1;
            else
              error ("RegressionGAM: DoF and Knots have been set already.");
            endif

          case 'dof'
            if (KOD < 2)
              DoF = varargin{2};
              if (! isnumeric (DoF) ||
                  ! (isscalar (DoF) || isequal (size (DoF), [1, ndims_X])))
                error ("RegressionGAM: invalid value for DoF.");
              endif
              Knots = DoF - Order;
              Order = DoF - Knots;
              KOD += 1;
            else
              error ("RegressionGAM: Knots and Order have been set already.");
            endif

          case 'tol'
            Tol = varargin{2};
            if (! (isnumeric (Tol) && isscalar (Tol) && (Tol > 0)))
              error ("RegressionGAM: Tolerance must be a Positive scalar.");
            endif

          otherwise
            error (strcat ("RegressionGAM: invalid parameter name", ...
                           " in optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Assign original X and Y data to the RegressionGAM object
      this.X = X;
      this.Y = Y;

      ## An observation is dropped only when its response is missing.  A row
      ## whose predictors hold missing values is kept and reported as used,
      ## while the fit below draws on the complete observations alone.
      RowsUsed  = ! isnan (Y(:));
      Yret      = Y(RowsUsed);
      Xret      = X(RowsUsed, :);
      this.X    = Xret;
      this.Y    = Yret;
      cobs      = ! any (isnan (Xret), 2);
      Y         = Yret(cobs);
      X         = Xret(cobs, :);

      ## Check X and Y contain valid data
      if (! isnumeric (X) || ! all (isfinite (X(:))))
        error ("RegressionGAM: invalid values in X.");
      endif
      if (! isnumeric (Y) || ! all (isfinite (Y(:))))
        error ("RegressionGAM: invalid values in Y.");
      endif

      ## Assign the number of observations and their corresponding indices
      ## on the original data, which will be used for training the model,
      ## to the RegressionGAM object
      this.NumObservations = rows (this.X);
      ## RowsUsed is left empty when every observation was used, as in MATLAB
      if (all (RowsUsed))
        this.RowsUsed = [];
      else
        this.RowsUsed = RowsUsed;
      endif

      ## Assign the number of original predictors to the RegressionGAM object
      this.NumPredictors = ndims_X;

      ## Generate default predictors and response variable names (if necessary)
      if (isempty (PredictorNames))
        for i = 1:ndims_X
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = 'Y';
      endif

      ## Assign predictors and response variable names
      this.PredictorNames = PredictorNames;
      this.ResponseName   = ResponseName;

      ## A scalar 'Knots', 'Order' or 'DoF' applies to every predictor, which
      ## the validation above accepts but nothing expanded, so the fit indexed
      ## past the end of a scalar for the second predictor onwards.
      if (isscalar (Knots))
        Knots = repmat (Knots, 1, ndims_X);
      endif
      if (isscalar (Order))
        Order = repmat (Order, 1, ndims_X);
      endif
      DoF = Knots + Order;

      ## Assign remaining optional parameters
      this.Formula      = Formula;
      this.Interactions = Interactions;
      this.Knots        = Knots;
      this.Order        = Order;
      this.DoF          = DoF;
      this.Tol          = Tol;
      this.ResponseTransform = ResponseTransform;
      this.RTname            = RTname;

      ## Bookkeeping MATLAB reports alongside the fit
      this.CategoricalPredictors  = [];
      this.ExpandedPredictorNames = PredictorNames;
      this.W = ones (this.NumObservations, 1) / this.NumObservations;
      this.IsStandardDeviationFit = false;

      ## Fit the basic model
      Inter = mean (Y);
      [iter, param, res, RSS] = this.fitGAM (X, Y, Inter, Knots, Order);
      this.BaseModel.Intercept  = Inter;
      this.BaseModel.Parameters = param;
      this.BaseModel.Iterations = iter;
      this.BaseModel.Residuals  = res;
      this.BaseModel.RSS        = RSS;
      this.Intercept            = Inter;

      ## Handle interaction terms (if given)
      if (F_I > 0)
        if (isempty (this.Formula))
          ## Analyze Interactions optional parameter
          this.IntMatrix = this.parseInteractions ();
          ## Append interaction terms to the predictor matrix
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = X(:,tindex);
            Xinter = ones (this.NumObservations, 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append interaction terms
            X = [X, Xinter];
          endfor

        else
          ## Analyze Formula optional parameter
          this.IntMatrix = this.parseFormula ();
          ## Add selected predictors and interaction terms
          XN = [];
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = X(:,tindex);
            Xinter = ones (this.NumObservations, 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append selected predictors and interaction terms
            XN = [XN, Xinter];
          endfor
          X = XN;
        endif

        ## Update length of Knots, Order, and DoF vectors to match
        ## the columns of X with the interaction terms
        Knots = ones (1, columns (X)) * Knots(1); # Knots
        Order = ones (1, columns (X)) * Order(1); # Order of spline
        DoF   = ones (1, columns (X)) * DoF(1);   # Degrees of freedom

        ## Fit the model with interactions
        [iter, param, res, RSS] = this.fitGAM (X, Y, Inter, Knots, Order);
        this.ModelwInt.Intercept  = Inter;
        this.ModelwInt.Parameters = param;
        this.ModelwInt.Iterations = iter;
        this.ModelwInt.Residuals  = res;
        this.ModelwInt.RSS        = RSS;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGAM} {@var{yFit} =} predict (@var{obj}, @var{Xfit})
    ## @deftypefnx {RegressionGAM} {@var{yFit} =} predict (@dots{}, @var{Name}, @var{Value})
    ## @deftypefnx {RegressionGAM} {[@var{yFit}, @var{ySD}, @var{yInt}] =} predict (@dots{})
    ##
    ## Predict new data points using generalized additive model regression
    ## object.
    ##
    ## @code{@var{yFit} = predict (@var{obj}, @var{Xfit}} returns a vector of
    ## predicted responses, @var{yFit}, for the predictor data in matrix
    ## @var{Xfit} based on the Generalized Additive Model in @var{obj}.
    ## @var{Xfit} must have the same number of features/variables as the
    ## training data in @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{RegressionGAM} class object.
    ## @end itemize
    ##
    ## @code{[@var{yFit}, @var{ySD}, @var{yInt}] = predict (@var{obj},
    ## @var{Xfit}}
    ## also returns the standard deviations, @var{ySD}, and prediction
    ## intervals,
    ## @var{yInt}, of the response variable @var{yFit}, evaluated at each
    ## observation in the predictor data @var{Xfit}.
    ##
    ## @code{@var{yFit} = predict (@dots{}, @var{Name}, @var{Value})} returns
    ## the
    ## aforementioned results with additional properties specified by
    ## @qcode{Name-Value} pair arguments listed below.
    ##
    ## @multitable @columnfractions 0.28 0.7
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'alpha'} @tab significance level of the prediction
    ## intervals @var{yInt}, specified as scalar in range @qcode{[0,1]}. The
    ## default value is 0.05, which corresponds to 95% prediction intervals.
    ##
    ## @item @qcode{'includeinteractions'} @tab a boolean flag to include
    ## interactions to predict new values based on @var{Xfit}.  By default,
    ## @qcode{'includeinteractions'} is @qcode{true} when the GAM model in
    ## @var{obj}
    ## contains a @qcode{obj.Formula} or @qcode{obj.Interactions} fields.
    ## Otherwise, is set to @qcode{false}. If set to @qcode{true} when no
    ## interactions are present in the trained model, it will result to an
    ## error. If set to
    ## @qcode{false} when using a model that includes interactions, the
    ## predictions
    ## will be made on the basic model without any interaction terms. This way
    ## you can make predictions from the same GAM model without having to
    ## retrain it.
    ## @end multitable
    ##
    ## @seealso{fitrgam, RegressionGAM}
    ## @end deftypefn
    function [yFit, ySD, yInt] = predict (this, Xfit, varargin)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("RegressionGAM.predict: too few arguments.");
      endif

      ## Check for valid XC
      if (isempty (Xfit))
        error ("RegressionGAM.predict: Xfit is empty.");
      elseif (columns (this.X) != columns (Xfit))
        error (strcat ("@RegressionGAM/predict: Xfit must have the same", ...
                       " number of features (columns) as in the GAM model."));
      endif

      ## Clean Xfit data
      notnansf  = ! logical (sum (isnan (Xfit), 2));
      Xfit      = Xfit(notnansf, :);

      ## Default values for Name-Value Pairs
      alpha = 0.05;
      if (isempty (this.IntMatrix))
        incInt = false;
      else
        incInt = true;
      endif

      ## Parse optional arguments
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'includeinteractions'
            tmpInt = varargin{2};
            if (! islogical (tmpInt) || (tmpInt != 0 && tmpInt != 1))
              error (strcat ("RegressionGAM.predict: includeinteractions", ...
                             " must be a logical value."));
            endif
            ## Check model for interactions
            if (tmpInt && isempty (this.IntMatrix))
              error (strcat ("RegressionGAM.predict: trained model", ...
                             " does not include any interactions."));
            endif
            incInt = tmpInt;

          case 'alpha'
            alpha = varargin{2};
            if (! (isnumeric (alpha) && isscalar (alpha)
                                      && alpha > 0 && alpha < 1))
              error (strcat ("RegressionGAM.predict: alpha must be a", ...
                             " scalar value between 0 and 1."));
            endif

          otherwise
            error (strcat ("RegressionGAM.predict: invalid NAME in", ...
                          " optional pairs of arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      ## Choose whether interactions must be included
      if (incInt)
        if (! isempty (this.Interactions))
          ## Append interaction terms to the predictor matrix
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = Xfit(:,tindex);
            Xinter = ones (rows (Xfit), 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append interaction terms
            Xfit = [Xfit, Xinter];
          endfor
        else
          ## Add selected predictors and interaction terms
          XN = [];
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = Xfit(:,tindex);
            Xinter = ones (rows (Xfit), 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append selected predictors and interaction terms
            XN = [XN, Xinter];
          endfor
          Xfit = XN;
        endif
        ## Get parameters and intercept vectors from model with interactions
        params = this.ModelwInt.Parameters;
        Interc = this.ModelwInt.Intercept;
        ## Update length of DoF vector
        DoF = ones (1, columns (Xfit)) * this.DoF(1);
      else
        ## Get parameters and intercept vectors from base model
        params = this.BaseModel.Parameters;
        Interc = this.BaseModel.Intercept;
        ## Get DoF from model
        DoF = this.DoF;
      endif


      ## Predict values from testing data
      yFit = predict_val (params, Xfit, Interc);
      yFit = this.ResponseTransform (yFit);

      ## Predict Standard Deviation and Intervals of estimated data if requested
      if (nargout > 1)
        ## Ensure that RowsUsed in the model are selected
        used = true (rows (this.X), 1);
        Y = this.Y(used);
        X = this.X(used, :);
        ## Predict response from training predictor data with the trained model
        yrs = predict_val (params, X , Interc);
        yrs_fit = predict_val (params, Xfit, Interc);

        ## Get the residuals between predicted and actual response data
        rs     = Y - yrs;
        var_rs = var (rs);

        t_mul  = tinv (1 - alpha / 2, this.DoF);

        ySD    = sqrt (var_rs) * ones (rows (yFit), 1);

        if (nargout > 2)
          moe    = t_mul(1) * ySD;
          lower  = this.ResponseTransform (yrs_fit - moe);
          upper  = this.ResponseTransform (yrs_fit + moe);
          yInt   = [lower, upper];
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGAM} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionGAM} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss of a generalized additive model.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the weighted
    ## mean squared error of the model on the rows of @var{X} against the true
    ## response @var{Y}.
    ##
    ## @code{@var{L} = loss (@dots{}, @var{name}, @var{value})} accepts the
    ## following name-value pairs:
    ##
    ## @itemize
    ## @item
    ## @qcode{"LossFun"} selects the loss, either @qcode{"mse"}, the default,
    ## or a function handle taking the true response, the predicted response
    ## and the weights, and returning a numeric scalar.
    ##
    ## @item
    ## @qcode{"Weights"} holds one weight per row of @var{X}, normalised to
    ## sum to one before it is applied.
    ## @end itemize
    ##
    ## @seealso{RegressionGAM, fitrgam, predict}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("RegressionGAM.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionGAM.loss: Name-Value arguments must", ...
                       " be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, 'loss');

      ## Defaults, then the optional pairs
      LossFun = 'mse';
      args = varargin;
      keep = true (1, numel (args));
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("RegressionGAM.loss: parameter name must be a", ...
                         " character vector."));
        endif
        if (strcmpi (args{i}, 'lossfun'))
          LossFun = args{i+1};
          if (! (is_function_handle (LossFun) ||
                 (ischar (LossFun) && isrow (LossFun))))
            error (strcat ("RegressionGAM.loss: 'LossFun' must be a", ...
                           " character vector or a function handle."));
          endif
          if (ischar (LossFun) && ! strcmpi (LossFun, 'mse'))
            error ("RegressionGAM.loss: unsupported 'LossFun' value.");
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
          error (strcat ("RegressionGAM.loss: 'LossFun' must return a", ...
                         " numeric scalar."));
        endif
      else
        L = sum (W .* (Y - yFit) .^ 2);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGAM} {@var{yFit} =} resubPredict (@var{obj})
    ##
    ## Predict the training response with the model it was fitted on.
    ##
    ## @code{@var{yFit} = resubPredict (@var{obj})} is @code{predict} applied
    ## to the observations the model was fitted on.
    ##
    ## @seealso{RegressionGAM, fitrgam, predict}
    ## @end deftypefn
    function yFit = resubPredict (this)
      used = true (rows (this.X), 1);
      yFit = predict (this, this.X(used, :));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGAM} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {RegressionGAM} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss of a generalized additive model on its training data.
    ##
    ## @code{@var{L} = resubLoss (@var{obj})} returns the weighted mean
    ## squared error of the model on the data it was fitted on.  It accepts
    ## the same @qcode{Name-Value} pairs as @code{loss}.
    ##
    ## @seealso{RegressionGAM, fitrgam, loss}
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      used = true (rows (this.X), 1);
      L = loss (this, this.X(used, :), this.Y(used), varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGAM} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Create a @qcode{CompactRegressionGAM} object.
    ##
    ## @code{@var{CMdl} = compact (@var{obj})} returns a compact version of
    ## the model, which predicts as it does but keeps no training data.
    ##
    ## @seealso{RegressionGAM, CompactRegressionGAM, fitrgam}
    ## @end deftypefn
    function CMdl = compact (this)
      CMdl = CompactRegressionGAM (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGAM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a RegressionGAM object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves a RegressionGAM
    ## object into a file defined by @var{filename}.
    ##
    ## @seealso{loadmodel, fitrgam, RegressionGAM}
    ## @end deftypefn

    function savemodel (obj, fname)
      if (nargin < 2)
        error ("RegressionGAM.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error ("RegressionGAM.savemodel: FNAME must be a character vector.");
      endif
      ## Generate variable for class name
      classdef_name = 'RegressionGAM';

      ## Create variables from model properties
      X = obj.X;
      Y = obj.Y;
      NumObservations     = obj.NumObservations;
      RowsUsed            = obj.RowsUsed;
      NumPredictors       = obj.NumPredictors;
      PredictorNames      = obj.PredictorNames;
      ResponseName        = obj.ResponseName;

      Formula             = obj.Formula;
      Interactions        = obj.Interactions;
      Knots               = obj.Knots;
      Order               = obj.Order;
      DoF                 = obj.DoF;
      Tol                 = obj.Tol;

      BaseModel           = obj.BaseModel;
      ModelwInt           = obj.ModelwInt;
      IntMatrix           = obj.IntMatrix;

      CategoricalPredictors  = obj.CategoricalPredictors;
      ExpandedPredictorNames = obj.ExpandedPredictorNames;
      W                      = obj.W;
      ResponseTransform      = obj.ResponseTransform;
      Intercept              = obj.Intercept;
      IsStandardDeviationFit = obj.IsStandardDeviationFit;
      RTname                 = obj.RTname;

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'NumObservations', ...
            'RowsUsed', 'NumPredictors', 'PredictorNames', 'ResponseName', ...
            'Formula', 'Interactions', 'Knots', 'Order', 'DoF', 'Tol', ...
            'BaseModel', 'ModelwInt', 'IntMatrix', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'W', 'ResponseTransform', ...
            'Intercept', 'IsStandardDeviationFit', 'RTname');
    endfunction

  endmethods

  ## Helper functions
  methods(Access = private)

    ## Shared validation for the assessment methods, so each reports under
    ## its own name.
    function [X, Y] = checkXY_ (this, X, Y, caller)
      if (isempty (X))
        error ("RegressionGAM.%s: X is empty.", caller);
      elseif (this.NumPredictors != columns (X))
        error (strcat ("RegressionGAM.%s: X must have the same number of", ...
                       " predictors as the trained model."), caller);
      endif
      if (isempty (Y))
        error ("RegressionGAM.%s: Y is empty.", caller);
      elseif (rows (X) != rows (Y))
        error (strcat ("RegressionGAM.%s: Y must have the same number of", ...
                       " rows as X."), caller);
      endif
    endfunction

    ## Pull a "Weights" pair out of the optional arguments, defaulting to a
    ## uniform weight, and reject any other name.
    function W = getWeights_ (this, args, n, caller)
      W = ones (n, 1);
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("RegressionGAM.%s: parameter name must be a", ...
                         " character vector."), caller);
        endif
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
          if (! (isnumeric (W) && isvector (W)))
            error (strcat ("RegressionGAM.%s: 'Weights' must be a numeric", ...
                           " vector."), caller);
          endif
          if (numel (W) != n)
            error (strcat ("RegressionGAM.%s: size of 'Weights' must equal", ...
                           " the number of rows in X."), caller);
          endif
        else
          error (strcat ("RegressionGAM.%s: invalid parameter name in", ...
                         " optional paired arguments."), caller);
        endif
      endfor

    endfunction

    ## Determine interactions from Interactions optional parameter
    function intMat = parseInteractions (this)
      if (islogical (this.Interactions))
        ## Check that interaction matrix corresponds to predictors
        if (numel (this.PredictorNames) != columns (this.Interactions))
          error (strcat ("RegressionGAM: columns in 'Interactions'", ...
                         " matrix must equal to the number of predictors."));
        endif
        intMat = this.Interactions;
      elseif (isnumeric (this.Interactions))
        ## Need to measure the effect of all interactions to keep the best
        ## performing. Just check that the given number is not higher than
        ## p*(p-1)/2, where p is the number of predictors.
        p = this.NumPredictors;
        if (this.Interactions > p * (p - 1) / 2)
          error (strcat ("RegressionGAM: number of interaction terms", ...
                         " requested is larger than all possible", ...
                         " combinations of predictors in X."));
        endif
        ## The pairs are not ranked by how much each contributes, so the
        ## first ones asked for are taken in the order nchoosek lists them.
        intMat = pairTerms (p)(1:this.Interactions, :);
      elseif (strcmpi (this.Interactions, 'all'))
        ## Calculate all p*(p-1)/2 interaction terms
        intMat = pairTerms (this.NumPredictors);
      endif
    endfunction

    ## Determine interactions from formula
    function intMat = parseFormula (this)
      intMat = [];
      ## Check formula for syntax
      if (isempty (strfind (this.Formula, '~')))
        error ("RegressionGAM: invalid syntax in Formula.");
      endif
      ## Split formula and keep predictor terms
      formulaParts = strsplit (this.Formula, '~');
      ## Check there is some string after '~'
      if (numel (formulaParts) < 2)
        error ("RegressionGAM: no predictor terms in Formula.");
      endif
      predictorString = strtrim (formulaParts{2});
      if (isempty (predictorString))
        error ("RegressionGAM: no predictor terms in Formula.");
      endif
      ## Split additive terms (between + sign)
      aterms = strtrim (strsplit (predictorString, '+'));
      ## Process all terms
      for i = 1:numel (aterms)
        ## Find individual terms (string missing ':')
        if (isempty (strfind (aterms(i), ':'){:}))
          ## Search PredictorNames to associate with column in X
          sterms = strcmp (this.PredictorNames, aterms(i));
          ## Append to interactions matrix
          intMat = [intMat; sterms];
        else
          ## Split interaction terms (string contains ':')
          mterms = strsplit (aterms{i}, ':');
          ## Add each individual predictor to interaction term vector
          iterms = logical (zeros (1, this.NumPredictors));
          for t = 1:numel (mterms)
            iterms = iterms | strcmp (this.PredictorNames, mterms(t));
          endfor
          ## Check that all predictors have been identified
          if (sum (iterms) != t)
            error ("RegressionGAM: some predictors have not been identified.");
          endif
          ## Append to interactions matrix
          intMat = [intMat; iterms];
        endif
      endfor
      ## Check that all terms have been identified
      if (! all (sum (intMat, 2) > 0))
        error ("RegressionGAM: some terms have not been identified.");
      endif
    endfunction

    ## Fit the model
    function [iter, param, res, RSS] = fitGAM (this, X, Y, Inter, Knots, Order)
      ## Initialize variables
      converged = false;
      iter      = 0;
      RSS       = zeros (1, columns (X));
      res       = Y - Inter;
      ns        = rows (X);
      Tol       = this.Tol;

      ## Every backfitting cycle refits the same predictor at the same knots
      ## and order and varies only the partial residual, so each predictor's
      ## design is fixed and its basis and factorisation are built once.
      tmpl  = cell (1, columns (X));
      cache = cell (1, columns (X));
      for j = 1:columns (X)
        [tmpl{j}, cache{j}] = splineCache (X(:,j), Knots(j), Order(j));
      endfor

      ## What each predictor last contributed, so a cycle can take it back out
      ## of the partial residual without evaluating the spline again.
      contrib = cell (1, columns (X));

      ## Start training
      while (! (converged || iter > 1000))
        iter += 1;

        ## Calculate residuals to fit spline
        for j = 1:columns (X)

          ## Take this predictor's own contribution back out
          if (iter > 1)
            res = res + contrib{j};
          endif

          ## Fit an spline to the data
          if (isempty (cache{j}))
            gk = splinefit (X(:,j), res, Knots(j), 'order', Order(j));
            gk_pred = ppval (gk, X(:,j));
            gk_pred = gk_pred(:);
          else
            b = cache{j};
            c = b.R \ (b.Q' * res);
            gk = tmpl{j};
            gk.coefs = reshape (sum (b.C .* c, 1), ...
                                [tmpl{j}.pieces, tmpl{j}.order]);
            gk_pred = b.F * c;
          endif

          ## This might be wrong! We need to check this out
          RSSk(j) = abs (sum (abs (Y - gk_pred - Inter)) .^ 2) / ns;
          param(j) = gk;
          contrib{j} = gk_pred;
          res = res - gk_pred;
        endfor

        ## Check if RSS is less than the tolerance
        if (all (abs (RSS - RSSk) <= Tol))
          converged = true;
        endif

        ## Update RSS
        RSS = RSSk;
      endwhile

    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a RegressionGAM object
      mdl = RegressionGAM (1, 1);

      ## Get fieldnames from DATA (including private properties)
      names = fieldnames (data);

      ## Copy data into object
      for i = 1:numel (names)
        ## Check that fieldnames in DATA match properties in RegressionGAM
        try
          mdl.(names{i}) = data.(names{i});
        catch
          error ("RegressionGAM.load_model: invalid model in '%s'.", ...
                 filename)
        end_try_catch
      endfor

      ## A model written before RowsUsed became a mask stored it as a
      ## double, which is a valid subscript for nothing.  An empty RowsUsed
      ## means every observation was used and stays an empty double.
      if (! isempty (mdl.RowsUsed))
        mdl.RowsUsed = logical (mdl.RowsUsed);
      endif
    endfunction

  endmethods

endclassdef

## Helper function for making prediction of new data based on GAM model
function ypred = predict_val (params, X, intercept)
  [nsample, ndims_X] = size (X);
  ypred = ones (nsample, 1) * intercept;
  ## Add the remaining terms
  for j = 1:ndims_X
    ypred = ypred + ppval (params(j), X(:,j));
  endfor
endfunction

%!demo
%! ## Train a RegressionGAM Model for synthetic values
%! f1 = @(x) cos (3 * x);
%! f2 = @(x) x .^ 3;
%! x1 = 2 * rand (50, 1) - 1;
%! x2 = 2 * rand (50, 1) - 1;
%! y = f1(x1) + f2(x2);
%! y = y + y .* 0.2 .* rand (50,1);
%! X = [x1, x2];
%! a = fitrgam (X, y, 'tol', 1e-3)

%!demo
%! ## Declare two different functions
%! f1 = @(x) cos (3 * x);
%! f2 = @(x) x .^ 3;
%!
%! ## Generate 80 samples for f1 and f2
%! x = [-4*pi:0.1*pi:4*pi-0.1*pi]';
%! X1 = f1(x);
%! X2 = f2(x);
%!
%! ## Create a synthetic response by adding noise
%! rand ('seed', 3);
%! Ytrue = X1 + X2;
%! Y = Ytrue + Ytrue .* 0.2 .* rand (80,1);
%!
%! ## Assemble predictor data
%! X = [X1, X2];
%!
%! ## Train the GAM and test on the same data
%! a = fitrgam (X, Y, 'order', [5, 5]);
%! [ypred, ySDsd, yInt] = predict (a, X);
%!
%! ## Plot the results
%! figure
%! [sortedY, indY] = sort (Ytrue);
%! plot (sortedY, 'r-');
%! xlim ([0, 80]);
%! hold on
%! plot (ypred(indY), 'g+')
%! plot (yInt(indY,1), 'k:')
%! plot (yInt(indY,2), 'k:')
%! xlabel ('Predictor samples');
%! ylabel ('Response');
%! title ('actual vs predicted values for function f1(x) = cos (3x) ');
%! legend ({'Theoretical Response', 'Predicted Response', 'Prediction Intervals'});
%!
%! ## Use 30% Holdout partitioning for training and testing data
%! C = cvpartition (80, 'HoldOut', 0.3);
%! [ypred, ySDsd, yInt] = predict (a, X(test (C),:));
%!
%! ## Plot the results
%! figure
%! [sortedY, indY] = sort (Ytrue(test (C)));
%! plot (sortedY, 'r-');
%! xlim ([0, sum(test(C))]);
%! hold on
%! plot (ypred(indY), 'g+')
%! plot (yInt(indY,1),'k:')
%! plot (yInt(indY,2),'k:')
%! xlabel ('Predictor samples');
%! ylabel ('Response');
%! title ('actual vs predicted values for function f1(x) = cos (3x) ');
%! legend ({'Theoretical Response', 'Predicted Response', 'Prediction Intervals'});

## Test constructor
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [1; 2; 3; 4];
%! a = RegressionGAM (x, y);
%! assert_equal ({a.X, a.Y}, {x, y})
%! assert_equal ({a.BaseModel.Intercept}, {2.5000})
%! assert_equal ({a.Knots, a.Order, a.DoF}, {[5, 5, 5], [3, 3, 3], [8, 8, 8]})
%! assert_equal ({a.NumObservations, a.NumPredictors}, {4, 3})
%! assert_equal ({a.ResponseName, a.PredictorNames}, {'Y', {'x1', 'x2', 'x3'}})
%! assert_equal ({a.Formula}, {[]})
%!test
%! x = [1, 2, 3, 4; 4, 5, 6, 7; 7, 8, 9, 1; 3, 2, 1, 2];
%! y = [1; 2; 3; 4];
%! pnames = {'A', 'B', 'C', 'D'};
%! formula = 'Y ~ A + B + C + D + A:C';
%! intMat = logical ([1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1;1,0,1,0]);
%! a = RegressionGAM (x, y, 'predictors', pnames, 'formula', formula);
%! assert_equal (a.IntMatrix, double (intMat))
%! assert_equal ({a.ResponseName, a.PredictorNames}, {'Y', pnames})
%! assert_equal (a.Formula, formula)

%!test
%! ## Test that predict() executes correctly when interactions are present
%! X = [1, 2; 3, 4; 5, 6; 7, 8];
%! Y = [10; 20; 30; 40];
%! mdl = RegressionGAM (X, Y, 'formula', 'Y ~ x1 + x2 + x1:x2');
%! ypred = predict (mdl, X);
%! assert_equal (isnumeric (ypred), true);
%! assert_equal (size (ypred), [4, 1]);
%! [ypred2, ySD, yInt] = predict (mdl, X, 'includeinteractions', true);
%! assert_equal (size (ypred2), [4, 1]);
%! assert_equal (size (ySD), [4, 1]);
%! assert_equal (size (yInt), [4, 2]);

%!test
%! ## Verify ySD is based on training residual variance
%! X = (1:10)';
%! Y = [2; 1; 4; 3; 6; 5; 8; 7; 10; 9];
%! mdl = RegressionGAM (X, Y);
%! y_train = predict (mdl, X);
%! rs = Y - y_train;
%! expected_ySD = sqrt (var (rs));
%! [~, ySD] = predict (mdl, X(1:4,:));
%! assert_equal (ySD, expected_ySD * ones (4, 1), 1e-10);

%!test
%! ## Verify ySD remains the same for one or more prediction points
%! X = (1:10)';
%! Y = [2; 1; 4; 3; 6; 5; 8; 7; 10; 9];
%! mdl = RegressionGAM (X, Y);
%! y_train = predict (mdl, X);
%! expected_ySD = sqrt (var (Y - y_train));
%! [~, ySD_1] = predict (mdl, X(1,:));
%! [~, ySD_3] = predict (mdl, X(1:3,:));
%! assert_equal (ySD_1, expected_ySD, 1e-10);
%! assert_equal (ySD_3, expected_ySD * ones (3, 1), 1e-10);

## Test input validation for constructor
%!error<RegressionGAM: too few input arguments.> RegressionGAM ()
%!error<RegressionGAM: too few input arguments.> RegressionGAM (ones (10,2))
%!error<RegressionGAM: number of rows in X and Y must be equal.> ...
%! RegressionGAM (ones (10,2), ones (5,1))
%!error<RegressionGAM: invalid values in X.> ...
%! RegressionGAM ([1;2;3;'a';4], ones (5,1))
%!error<RegressionGAM: invalid parameter name in optional pair arguments.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'some', 'some')
%!error<RegressionGAM: Formula must be a string.>
%! RegressionGAM (ones (10,2), ones (10,1), 'formula', {'y~x1+x2'})
%!error<RegressionGAM: Formula must be a string.>
%! RegressionGAM (ones (10,2), ones (10,1), 'formula', [0, 1, 0])
%!error<RegressionGAM: invalid syntax in Formula.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'formula', 'something')
%!error<RegressionGAM: no predictor terms in Formula.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'formula', 'something~')
%!error<RegressionGAM: no predictor terms in Formula.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'formula', 'something~')
%!error<RegressionGAM: some predictors have not been identified> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'formula', 'something~x1:')
%!error<RegressionGAM: invalid Interactions parameter.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'interactions', 'some')
%!error<RegressionGAM: invalid Interactions parameter.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'interactions', -1)
%!error<RegressionGAM: invalid Interactions parameter.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'interactions', [1 2 3 4])
%!error<RegressionGAM: number of interaction terms requested is larger than> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'interactions', 3)
%!error<RegressionGAM: Formula has been already defined.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'formula', 'y ~ x1 + x2', 'interactions', 1)
%!error<RegressionGAM: Interactions have been already defined.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'interactions', 1, 'formula', 'y ~ x1 + x2')
%!error<RegressionGAM: invalid value for Knots.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'knots', 'a')
%!error<RegressionGAM: DoF and Order have been set already.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'order', 3, 'dof', 2, 'knots', 5)
%!error<RegressionGAM: invalid value for DoF.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'dof', 'a')
%!error<RegressionGAM: Knots and Order have been set already.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'knots', 5, 'order', 3, 'dof', 2)
%!error<RegressionGAM: invalid value for Order.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'order', 'a')
%!error<RegressionGAM: DoF and Knots have been set already.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'knots', 5, 'dof', 2, 'order', 2)
%!error<RegressionGAM: Tolerance must be a Positive scalar.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'tol', -1)
%!error<RegressionGAM: ResponseName must be a char string.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'responsename', -1)
%!error<RegressionGAM: PredictorNames must be a cellstring array.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'predictors', -1)
%!error<RegressionGAM: PredictorNames must be a cellstring array.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'predictors', ['a','b','c'])
%!error<RegressionGAM: PredictorNames must have same number of columns as X.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'predictors', {'a','b','c'})

## Test input validation for predict method
%!error<RegressionGAM.predict: too few arguments.> ...
%! predict (RegressionGAM (ones (10,1), ones (10,1)))
%!error<RegressionGAM.predict: Xfit is empty.> ...
%! predict (RegressionGAM (ones (10,1), ones (10,1)), [])
%!error<RegressionGAM.predict: Xfit must have the same number of features> ...
%! predict (RegressionGAM (ones (10,2), ones (10,1)), 2)
%!error<RegressionGAM.predict: invalid NAME in optional pairs of arguments.> ...
%! predict (RegressionGAM (ones (10,2), ones (10,1)), ones (10,2), 'some', 'some')
%!error<RegressionGAM.predict: includeinteractions must be a logical value.> ...
%! predict (RegressionGAM (ones (10,2), ones (10,1)), ones (10,2), 'includeinteractions', 'some')
%!error<RegressionGAM.predict: includeinteractions must be a logical value.> ...
%! predict (RegressionGAM (ones (10,2), ones (10,1)), ones (10,2), 'includeinteractions', 5)
%!error<RegressionGAM.predict: alpha must be a scalar value between 0 and 1.> ...
%! predict (RegressionGAM (ones (10,2), ones (10,1)), ones (10,2), 'alpha', 5)
%!error<RegressionGAM.predict: alpha must be a scalar value between 0 and 1.> ...
%! predict (RegressionGAM (ones (10,2), ones (10,1)), ones (10,2), 'alpha', -1)
%!error<RegressionGAM.predict: alpha must be a scalar value between 0 and 1.> ...
%! predict (RegressionGAM (ones (10,2), ones (10,1)), ones (10,2), 'alpha', 'a')
%!error <RegressionGAM.savemodel: too few input arguments.> ...
%! savemodel (RegressionGAM ([1, 2; 2, 3; 3, 4; 4, 5], [1; 2; 3; 4]))
%!error <RegressionGAM.savemodel: FNAME must be a character vector.> ...
%! savemodel (RegressionGAM ([1, 2; 2, 3; 3, 4; 4, 5], [1; 2; 3; 4]), 1)
%!error <RegressionGAM.savemodel: FNAME must be a character vector.> ...
%! savemodel (RegressionGAM ([1, 2; 2, 3; 3, 4; 4, 5], [1; 2; 3; 4]), ['ab'; 'cd'])

## The bookkeeping MATLAB reports alongside the fit is present.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,1:3), meas(:,4));
%! assert_equal (Mdl.Intercept, Mdl.BaseModel.Intercept);
%! assert_equal (size (Mdl.W), [Mdl.NumObservations, 1]);
%! assert_equal (sum (Mdl.W), 1, 1e-12);
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (Mdl.ExpandedPredictorNames, Mdl.PredictorNames);
%! assert_equal (Mdl.IsStandardDeviationFit, false);

## A scalar Knots, Order or DoF applies to every predictor.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,1:3), meas(:,4), 'Knots', 4);
%! assert_equal (Mdl.Knots, [4, 4, 4]);
%! assert_equal (Mdl.Order, [3, 3, 3]);
%! assert_equal (Mdl.DoF, [7, 7, 7]);

## A single non-finite value in X or Y is refused.
%!error<RegressionGAM: invalid values in X.> ...
%! RegressionGAM ([1, 2; Inf, 4; 5, 6; 7, 8], [1; 2; 3; 4])
%!error<RegressionGAM: invalid values in Y.> ...
%! RegressionGAM ([1, 2; 3, 4; 5, 6; 7, 8], [1; 2; Inf; 4])

## 'Interactions' names every pairwise term, one row per term.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,1:3), meas(:,4), 'Interactions', 'all');
%! assert_equal (Mdl.IntMatrix, logical ([1, 1, 0; 1, 0, 1; 0, 1, 1]));
%! assert_equal (sum (Mdl.IntMatrix(:)), 6);

## An assigned ResponseTransform reaches the predicted response.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! Mdl = fitrgam (X, meas(:,4));
%! y0 = predict (Mdl, X);
%! Mdl.ResponseTransform = 'exp';
%! assert_equal (predict (Mdl, X), exp (y0), 1e-12);

## resubPredict and resubLoss are predict and loss on the training data.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! Y = meas(:,4);
%! Mdl = fitrgam (X, Y);
%! assert_equal (resubPredict (Mdl), predict (Mdl, X));
%! assert_equal (resubLoss (Mdl), loss (Mdl, X, Y));

## loss is the weighted mean squared error, and takes a function of its own.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! Y = meas(:,4);
%! Mdl = fitrgam (X, Y);
%! r = Y - predict (Mdl, X);
%! assert_equal (loss (Mdl, X, Y), mean (r .^ 2), 1e-12);
%! w = [ones(75, 1); 3 * ones(75, 1)];
%! assert_equal (loss (Mdl, X, Y, 'Weights', w), ...
%!               sum (w .* r .^ 2) / sum (w), 1e-12);
%! mae = @(y, yfit, wt) sum (wt .* abs (y - yfit));
%! assert_equal (loss (Mdl, X, Y, 'LossFun', mae), mean (abs (r)), 1e-12);

## A saved and reloaded model carries every property and predicts alike.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! Mdl = fitrgam (X, meas(:,4), 'Interactions', 'all');
%! Mdl.ResponseTransform = 'exp';
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (Mdl2), 'RegressionGAM');
%! assert_equal (Mdl2.Intercept, Mdl.Intercept);
%! assert_equal (Mdl2.W, Mdl.W);
%! assert_equal (Mdl2.IntMatrix, Mdl.IntMatrix);
%! assert_equal (Mdl2.ResponseTransform (2), exp (2), 1e-12);
%! assert_equal (predict (Mdl2, X), predict (Mdl, X));

## Test input validation for loss method
%!shared xr, yr, Mr
%! load fisheriris
%! xr = meas(:,1:3);
%! yr = meas(:,4);
%! Mr = fitrgam (xr, yr);
%!error<RegressionGAM.loss: too few input arguments.> ...
%! loss (Mr, xr)
%!error<RegressionGAM.loss: Name-Value arguments must be in pairs.> ...
%! loss (Mr, xr, yr, 'Weights')
%!error<RegressionGAM.loss: X is empty.> ...
%! loss (Mr, [], yr)
%!error<RegressionGAM.loss: X must have the same number of predictors as the trained model.> ...
%! loss (Mr, 1, yr)
%!error<RegressionGAM.loss: Y must have the same number of rows as X.> ...
%! loss (Mr, xr, yr(1:10))
%!error<RegressionGAM.loss: unsupported 'LossFun' value.> ...
%! loss (Mr, xr, yr, 'LossFun', 'mad')
%!error<RegressionGAM.subsasgn: unrecognized or read-only property: 'Knots'> ...
%! Mr.Knots = 3;
%!error<RegressionGAM: unrecognized 'ResponseTransform' function.> ...
%! Mr.ResponseTransform = 'nonsense';

## RowsUsed is empty when every observation was used.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! Y = meas(:,1);
%! Mdl = fitrgam (X, Y);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (class (Mdl.RowsUsed), 'double');
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (rows (Mdl.X), 150);
%! assert_equal (rows (Mdl.W), 150);

## A missing response drops its observation and RowsUsed marks it.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! Y = meas(:,1);
%! Y(5) = NaN;
%! Mdl = fitrgam (X, Y);
%! assert_equal (class (Mdl.RowsUsed), 'logical');
%! assert_equal (size (Mdl.RowsUsed), [150, 1]);
%! assert_equal (sum (Mdl.RowsUsed), 149);
%! assert_equal (Mdl.RowsUsed(5), false);
%! assert_equal (Mdl.NumObservations, 149);
%! assert_equal (rows (Mdl.X), 149);
%! assert_equal (rows (Mdl.W), 149);

## A missing predictor keeps its observation, so RowsUsed stays empty.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! X(3,2) = NaN;
%! Y = meas(:,1);
%! Mdl = fitrgam (X, Y);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (rows (Mdl.X), 150);
%! assert_equal (sum (isnan (Mdl.X(:))), 1);
