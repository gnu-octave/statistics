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

classdef CompactRegressionGAM
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactRegressionGAM
  ##
  ## Compact generalized additive model regression
  ##
  ## The @code{CompactRegressionGAM} class implements a compact version of the
  ## generalized additive model regression object, which predicts responses for
  ## new data with the @code{predict} method but does not store the training
  ## data.
  ##
  ## A compact model consumes less memory than the full @code{RegressionGAM}
  ## model, but cannot perform tasks that need the training data, such as
  ## computing a resubstitution loss or the standard deviation of a prediction.
  ##
  ## Create a @code{CompactRegressionGAM} object by using the @code{compact}
  ## method on a @code{RegressionGAM} object.
  ##
  ## The engine that fitted the model is carried over in @code{FitMethod},
  ## and the compact model predicts by the same scheme the full one did.
  ## Under @qcode{'boostedtrees'}, the default, the fit is described by
  ## @code{TreeModel}, @code{BinEdges} and @code{PairDetectionBinEdges}.
  ## Under @qcode{'splines'} it is described by @code{Formula},
  ## @code{BaseModel}, @code{ModelwInt} and @code{IntMatrix}, which MATLAB's
  ## compact model does not carry.
  ## Whichever fitted the model, the other set is empty.  A standard
  ## deviation is available from the spline engine alone.
  ##
  ## @seealso{RegressionGAM, fitrgam}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer, the number of predictors of the training data.
    ## This property is read-only.
    ##
    ## @end deftp
    NumPredictors         = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} PredictorNames
    ##
    ## Names of the predictor variables
    ##
    ## A cell array of character vectors naming the predictors, in the order
    ## they appear in the training data.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector naming the response variable @var{Y}.  This
    ## property is read-only.
    ##
    ## @end deftp
    ResponseName          = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector holding the column of each predictor treated as
    ## categorical, and empty when none is.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} ExpandedPredictorNames
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
    ## @deftp {CompactRegressionGAM} {property} Intercept
    ##
    ## Intercept of the fitted model
    ##
    ## A numeric scalar, the mean of the response, which every additive term
    ## is measured against.  This property is read-only.
    ##
    ## @end deftp
    Intercept             = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} Formula
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
    ## @deftp {CompactRegressionGAM} {property} Interactions
    ##
    ## Two-way interaction terms of the fitted model
    ##
    ## A @math{Kx2} matrix of predictor index pairs, one row per two-way term
    ## the model carries, and @code{zeros (0, 2)} when it carries none.  It
    ## reports what was fitted rather than what was asked for, so a count of
    ## terms, @qcode{'all'}, a logical matrix and a formula all leave the same
    ## kind of value behind.  This property is read-only.
    ##
    ## A main effect names one predictor and a higher-order term names three
    ## or more, and neither has a two-column form, so neither appears here.
    ## @code{IntMatrix} remains the complete record of every term fitted.
    ##
    ## @end deftp
    Interactions          = zeros (0, 2);

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} IsStandardDeviationFit
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
    ## @deftp {CompactRegressionGAM} {property} BaseModel
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
    ## @deftp {CompactRegressionGAM} {property} ModelwInt
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
    ## @deftp {CompactRegressionGAM} {property} IntMatrix
    ##
    ## Every term the model fits
    ##
    ## A logical matrix with one row per term and one column per predictor,
    ## true wherever the term multiplies that predictor.  A row naming one
    ## predictor is a main effect, two an interaction, and three or more a
    ## higher-order term.  This property is read-only.
    ##
    ## It is the complete record, where @code{Interactions} reports only the
    ## two-way terms, in the form MATLAB reports them.  It is also the form
    ## the @qcode{'Interactions'} option takes back, so passing it to the
    ## constructor rebuilds a model over the same terms.
    ##
    ## @end deftp
    IntMatrix             = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} BinEdges
    ##
    ## Bin edges of the fitted shape functions, empty under the spline
    ## engine.  This property is read-only.
    ##
    ## @end deftp
    BinEdges = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} PairDetectionBinEdges
    ##
    ## Coarse bin edges the interaction terms are held on, empty when the
    ## model carries none.  This property is read-only.
    ##
    ## @end deftp
    PairDetectionBinEdges = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} FitMethod
    ##
    ## Which engine fitted the model, @qcode{'boostedtrees'} or
    ## @qcode{'splines'}.  This property is read-only.
    ##
    ## @end deftp
    FitMethod = 'boostedtrees';

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} TreeModel
    ##
    ## The fitted shape functions and interaction surfaces, empty under the
    ## spline engine.  This property is read-only.
    ##
    ## @end deftp
    TreeModel = [];

  endproperties

  ## Properties a user may set after the model is built.  Each one is
  ## validated by its set method below.
  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionGAM} {property} ResponseTransform
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

  endproperties

  ## Readable by the counterpart class, which copies it, and kept out of
  ## the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## Carried from the fitted model so a fold, which is stored
    ## compact, can still say how many trees it fitted.
    NumTrainedTrees = [];
    RTfun = @(y) y;
  endproperties

  ## Set methods for the properties a user may assign.
  methods (Hidden)

    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = parseResponseTransform ...
                                             (val, 'CompactRegressionGAM');
    endfunction

    ## constructor
    function this = CompactRegressionGAM (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), 'RegressionGAM'))
        error ("CompactRegressionGAM: invalid regression object.");
      endif

      ## Save properties to compact model
      this.NumPredictors          = Mdl.NumPredictors;
      this.PredictorNames         = Mdl.PredictorNames;
      this.ResponseName           = Mdl.ResponseName;
      this.CategoricalPredictors  = Mdl.CategoricalPredictors;
      this.ExpandedPredictorNames = Mdl.ExpandedPredictorNames;
      this.ResponseTransform      = Mdl.ResponseTransform;
      this.Intercept              = Mdl.Intercept;
      this.Formula                = Mdl.Formula;
      this.Interactions           = Mdl.Interactions;
      this.IsStandardDeviationFit = Mdl.IsStandardDeviationFit;
      this.BaseModel              = Mdl.BaseModel;
      this.ModelwInt              = Mdl.ModelwInt;
      this.IntMatrix              = Mdl.IntMatrix;
      this.RTfun                 = Mdl.RTfun;
      this.FitMethod             = Mdl.FitMethod;
      this.TreeModel             = Mdl.TreeModel;
      this.BinEdges              = Mdl.BinEdges;
      this.PairDetectionBinEdges = Mdl.PairDetectionBinEdges;
      this.NumTrainedTrees = Mdl.NumTrainedTrees;

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
      fprintf ("\n  CompactRegressionGAM\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %d\n", 'NumPredictors', this.NumPredictors);
      fprintf ("%+25s: '%s'\n", 'ResponseTransform', this.ResponseTransform);
      fprintf ("%+25s: %g\n", 'Intercept', this.Intercept);
    endfunction

  endmethods

  methods (Access = public)
    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionGAM} {@var{yFit} =} predict (@var{obj}, @var{Xfit})
    ## @deftypefnx {CompactRegressionGAM} {@var{yFit} =} predict (@dots{}, @var{Name}, @var{Value})
    ## @deftypefnx {CompactRegressionGAM} {[@var{yFit}, @var{ySD}, @var{yInt}] =} predict (@dots{})
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
    ## @var{obj} must be a @qcode{CompactRegressionGAM} class object.
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
    function yFit = predict (this, Xfit, varargin)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("CompactRegressionGAM.predict: too few arguments.");
      endif

      ## Check for valid XC
      if (isempty (Xfit))
        error ("CompactRegressionGAM.predict: Xfit is empty.");
      elseif (this.NumPredictors != columns (Xfit))
        error (strcat ("CompactRegressionGAM.predict: Xfit must have the", ...
                       " same number of features (columns) as in the GAM model."));
      endif

      ## Clean Xfit data
      notnansf  = ! logical (sum (isnan (Xfit), 2));
      Xfit      = Xfit(notnansf, :);

      ## Default values for Name-Value Pairs
      alpha = 0.05;
      hasInt = ! isempty (this.IntMatrix);
      if (strcmp (this.FitMethod, 'boostedtrees') && ! isempty (this.TreeModel))
        hasInt = ! isempty (this.TreeModel.Pairs);
      endif
      if (! hasInt)
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
              error (strcat ("CompactRegressionGAM.predict:", ...
                             " includeinteractions must be a logical value."));
            endif
            ## Check model for interactions
            if (tmpInt && ! hasInt)
              error (strcat ("CompactRegressionGAM.predict: trained model", ...
                             " does not include any interactions."));
            endif
            incInt = tmpInt;

          case 'alpha'
            alpha = varargin{2};
            if (! (isnumeric (alpha) && isscalar (alpha)
                                      && alpha > 0 && alpha < 1))
              error (strcat ("CompactRegressionGAM.predict: alpha must be", ...
                             " a scalar value between 0 and 1."));
            endif

          otherwise
            error (strcat ("CompactRegressionGAM.predict: invalid NAME in", ...
                          " optional pairs of arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      ## Choose whether interactions must be included
      ## The boosted-tree engine keeps its fit as step functions over bins, so
      ## a term is a lookup rather than a spline evaluation.
      if (strcmp (this.FitMethod, 'boostedtrees') && ! isempty (this.TreeModel))
        ## Excluding the interactions means excluding the constant they
        ## handed the intercept as well.
        interc = this.Intercept;
        if (! incInt && isfield (this.TreeModel, 'PairIntercept'))
          interc = interc - this.TreeModel.PairIntercept;
        endif
        if (! incInt || isempty (this.TreeModel.Pairs))
          yFit = gamboostpredict (this.BinEdges, ...
                                  this.TreeModel.ShapeValues, Xfit, ...
                                  interc);
        else
          yFit = gamboostpredict (this.BinEdges, ...
                                  this.TreeModel.ShapeValues, Xfit, ...
                                  interc, 0, ...
                                  this.PairDetectionBinEdges, ...
                                  this.TreeModel.PairValues, ...
                                  this.TreeModel.Pairs);
        endif
        yFit = this.RTfun (yFit);
        return;
      endif

      if (incInt)
        ## Which construction path the model took: an interaction
        ## list appends its terms to the predictors, a formula
        ## names every term the model has and replaces them.
        if (isempty (this.Formula))
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
      else
        ## Get parameters and intercept vectors from base model
        params = this.BaseModel.Parameters;
        Interc = this.BaseModel.Intercept;
      endif

      ## Predict values from testing data
      yFit = predict_val (params, Xfit, Interc);
      yFit = this.RTfun (yFit);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionGAM} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactRegressionGAM} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
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
    ## @seealso{CompactRegressionGAM, RegressionGAM, fitrgam, predict}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("CompactRegressionGAM.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactRegressionGAM.loss: Name-Value arguments", ...
                       " must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, 'loss');

      ## Defaults, then the optional pairs
      LossFun = 'mse';
      args = varargin;
      keep = true (1, numel (args));
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("CompactRegressionGAM.loss: parameter name must", ...
                         " be a character vector."));
        endif
        if (strcmpi (args{i}, 'lossfun'))
          LossFun = args{i+1};
          if (! (is_function_handle (LossFun) ||
                 (ischar (LossFun) && isrow (LossFun))))
            error (strcat ("CompactRegressionGAM.loss: 'LossFun' must be a", ...
                           " character vector or a function handle."));
          endif
          if (ischar (LossFun) && ! strcmpi (LossFun, 'mse'))
            error ("CompactRegressionGAM.loss: unsupported 'LossFun' value.");
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
          error (strcat ("CompactRegressionGAM.loss: 'LossFun' must", ...
                         " return a numeric scalar."));
        endif
      else
        L = sum (W .* (Y - yFit) .^ 2);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionGAM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a CompactRegressionGAM object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## CompactRegressionGAM object into an Octave binary file, the name of
    ## which is specified in @var{filename}, along with an extra variable which
    ## defines the type of object these variables constitute.  Use
    ## @code{loadmodel} in order to load the object back into Octave.
    ##
    ## @seealso{loadmodel, fitrgam, RegressionGAM}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error ("CompactRegressionGAM.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("CompactRegressionGAM.savemodel: FNAME must be a", ...
                       " character vector."));
      endif
      ## Generate variable for class name
      classdef_name = 'CompactRegressionGAM';

      ## Create variables from model properties
      NumPredictors          = this.NumPredictors;
      PredictorNames         = this.PredictorNames;
      ResponseName           = this.ResponseName;
      CategoricalPredictors  = this.CategoricalPredictors;
      ExpandedPredictorNames = this.ExpandedPredictorNames;
      ResponseTransform      = this.ResponseTransform;
      Intercept              = this.Intercept;
      Formula                = this.Formula;
      Interactions           = this.Interactions;
      IsStandardDeviationFit = this.IsStandardDeviationFit;
      BaseModel              = this.BaseModel;
      ModelwInt              = this.ModelwInt;
      IntMatrix              = this.IntMatrix;
      RTfun                 = this.RTfun;
      FitMethod             = this.FitMethod;
      TreeModel             = this.TreeModel;
      BinEdges              = this.BinEdges;
      PairDetectionBinEdges = this.PairDetectionBinEdges;

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'NumPredictors', ...
            'PredictorNames', 'ResponseName', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'ResponseTransform', 'Intercept', ...
            'Formula', 'Interactions', 'IsStandardDeviationFit', ...
            'BaseModel', 'ModelwInt', 'IntMatrix', 'RTfun', 'FitMethod', ...
            'TreeModel', 'BinEdges', 'PairDetectionBinEdges');
    endfunction

  endmethods

  methods(Access = private)

    ## Shared validation for the assessment methods, so each reports under
    ## its own name.
    function [X, Y] = checkXY_ (this, X, Y, caller)
      if (isempty (X))
        error ("CompactRegressionGAM.%s: X is empty.", caller);
      elseif (this.NumPredictors != columns (X))
        error (strcat ("CompactRegressionGAM.%s: X must have the same", ...
                       " number of predictors as the trained model."), caller);
      endif
      if (isempty (Y))
        error ("CompactRegressionGAM.%s: Y is empty.", caller);
      elseif (rows (X) != rows (Y))
        error (strcat ("CompactRegressionGAM.%s: Y must have the same", ...
                       " number of rows as X."), caller);
      endif
    endfunction

    ## Pull a "Weights" pair out of the optional arguments, defaulting to a
    ## uniform weight, and reject any other name.
    function W = getWeights_ (this, args, n, caller)
      W = ones (n, 1);
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("CompactRegressionGAM.%s: parameter name must be", ...
                         " a character vector."), caller);
        endif
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
          if (! (isnumeric (W) && isvector (W)))
            error (strcat ("CompactRegressionGAM.%s: 'Weights' must be a", ...
                           " numeric vector."), caller);
          endif
          if (numel (W) != n)
            error (strcat ("CompactRegressionGAM.%s: size of 'Weights'", ...
                           " must equal the number of rows in X."), caller);
          endif
        else
          error (strcat ("CompactRegressionGAM.%s: invalid parameter name", ...
                         " in optional paired arguments."), caller);
        endif
      endfor
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a CompactRegressionGAM object
      mdl = CompactRegressionGAM ();

      ## Get fieldnames from DATA (including private properties)
      names = fieldnames (data);
      ## The set methods for these read other properties, and one of them
      ## rebuilds Coeffs, so they are assigned once everything else is in
      ## place rather than in the order the file happens to list them.
      late = ismember (names, {'Cost', 'Prior', 'ScoreTransform', ...
                               'ResponseTransform'});
      names = [names(! late); names(late)];

      ## Copy data into object
      for i = 1:numel (names)
        try
          mdl.(names{i}) = data.(names{i});
        catch
          error ("CompactRegressionGAM.load_model: invalid model in '%s'.", ...
                 filename);
        end_try_catch
      endfor
    endfunction

  endmethods

endclassdef

## Helper function
function ypred = predict_val (params, X, intercept)
  ## The shared prediction engine evaluates every additive term and adds the
  ## intercept.
  ypred = gampredict (params, X, intercept);
endfunction

%!demo
%! ## Take the compact version of a fitted model and predict with it
%!
%! load fisheriris
%! X = meas(:,1:3);
%! Y = meas(:,4);
%!
%! mdl = fitrgam (X, Y)
%! cmdl = compact (mdl)

## Test input validation for constructor
%!error<CompactRegressionGAM: invalid regression object.> ...
%! CompactRegressionGAM (1)

## The compact model carries what MATLAB's compact model reports.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,1:3), meas(:,4));
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), 'CompactRegressionGAM');
%! assert_equal (CMdl.Intercept, Mdl.Intercept);
%! assert_equal (CMdl.CategoricalPredictors, Mdl.CategoricalPredictors);
%! assert_equal (CMdl.ExpandedPredictorNames, Mdl.ExpandedPredictorNames);
%! assert_equal (CMdl.IsStandardDeviationFit, false);
%! assert_equal (isprop (CMdl, 'X'), false);

## predict and loss agree with the model it was compacted from.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! Y = meas(:,4);
%! Mdl = fitrgam (X, Y, 'Interactions', 'all');
%! CMdl = compact (Mdl);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));
%! assert_equal (loss (CMdl, X, Y), loss (Mdl, X, Y));

## A saved and reloaded compact model carries every property it holds.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! CMdl = compact (fitrgam (X, meas(:,4)));
%! fname = tempname ();
%! savemodel (CMdl, fname);
%! CMdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (CMdl2), 'CompactRegressionGAM');
%! assert_equal (CMdl2.Intercept, CMdl.Intercept);
%! assert_equal (predict (CMdl2, X), predict (CMdl, X));

## An assigned ResponseTransform reaches the predicted response.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! CMdl = compact (fitrgam (X, meas(:,4)));
%! y0 = predict (CMdl, X);
%! CMdl.ResponseTransform = 'log';
%! assert_equal (predict (CMdl, X), log (y0), 1e-12);

## Test input validation
%!shared xc, yc, CMr
%! load fisheriris
%! xc = meas(:,1:3);
%! yc = meas(:,4);
%! CMr = compact (fitrgam (xc, yc));
%!error<CompactRegressionGAM.predict: too few arguments.> ...
%! predict (CMr)
%!error<CompactRegressionGAM.predict: Xfit is empty.> ...
%! predict (CMr, [])
%!error<CompactRegressionGAM.loss: too few input arguments.> ...
%! loss (CMr, xc)
%!error<CompactRegressionGAM.loss: unsupported 'LossFun' value.> ...
%! loss (CMr, xc, yc, 'LossFun', 'mad')
%!error<CompactRegressionGAM.savemodel: too few input arguments.> ...
%! savemodel (CompactRegressionGAM ())
%!error<CompactRegressionGAM.savemodel: FNAME must be a character vector.> ...
%! savemodel (CompactRegressionGAM (), 1)

## A fitted model survives savemodel and loadmodel: the properties come
## back as they were and it predicts the same.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! Y = meas(:,1);
%! Mdl = compact (fitrgam (X, Y, 'FitMethod', 'splines'));
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'CompactRegressionGAM');
%! assert_equal (M2.PredictorNames, Mdl.PredictorNames);
%! assert_equal (class (M2.ResponseTransform), class (Mdl.ResponseTransform));
%! assert_equal (M2.BaseModel.Parameters(1).coefs, ...
%!               Mdl.BaseModel.Parameters(1).coefs);
%! assert_equal (predict (M2, X(1:5,:)), predict (Mdl, X(1:5,:)), 1e-12);

## The same round trip under the boosted-tree engine.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! Y = meas(:,1);
%! Mdl = compact (fitrgam (X, Y, 'FitMethod', 'boostedtrees'));
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (M2.FitMethod, 'boostedtrees');
%! assert_equal (M2.TreeModel.ShapeValues, Mdl.TreeModel.ShapeValues);
%! assert_equal (predict (M2, X(1:5,:)), predict (Mdl, X(1:5,:)), 1e-12);

## A compacted tree-fitted model predicts as the full model does.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! CMdl = compact (fitrgam (X, meas(:,1), 'FitMethod', 'boostedtrees'));
%! assert_equal (CMdl.FitMethod, 'boostedtrees');
%! assert_equal (numel (CMdl.BinEdges), 3);
%! assert_equal (numel (predict (CMdl, X)), rows (X));
