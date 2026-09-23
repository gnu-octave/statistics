## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2023-2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

classdef RegressionGAM < PredictiveModel
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
## @item @qcode{'Predictors'} @tab Predictor Variable names, specified as
## a row vector cell of strings with the same length as the columns in @var{X}.
## If omitted, the program will generate default variable names
## @qcode{(x1, x2, ..., xn)} for each column in @var{X}.
##
## @item @qcode{'responsename'} @tab Response Variable Name, specified as
## a string.  If omitted, the default value is @qcode{'Y'}.
##
## @item @qcode{'formula'} @tab (spline option) a model specification given as a
## string in
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
## @item @qcode{'knots'} @tab (spline option) a scalar or a row vector with the
## same
## columns as @var{X}.  It defines the knots for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 5, hence expanded to @qcode{ones (1, columns (X)) * 5}.  You can
## parse a row vector with different number of knots for each predictor
## variable to be fitted with, although not recommended.
##
## @item @qcode{'order'} @tab (spline option) a scalar or a row vector with the
## same
## columns as @var{X}.  It defines the order of the polynomial when training the
## GAM.  As a scalar, it is expanded to a row vector.  The default values is 3,
## hence expanded to @qcode{ones (1, columns (X)) * 3}.  You can parse a row
## vector with different number of polynomial order for each predictor variable
## to be fitted with, although not recommended.
##
## @item @qcode{'dof'} @tab (spline option) a scalar or a row vector with the
## same columns
## as @var{X}.  It defines the degrees of freedom for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 8, hence expanded to @qcode{ones (1, columns (X)) * 8}.  You can
## parse a row vector with different degrees of freedom for each predictor
## variable to be fitted with, although not recommended.
##
## @item @qcode{'tol'} @tab (spline option) a positive scalar to set the
## tolerance for
## convergence during training. By default, it is set to @qcode{1e-3}.
##
## @end multitable
##
## A row marked @qcode{(spline option)} belongs to the spline
## engine and requires @qcode{'FitMethod', 'splines'}; passing one
## under the default boosted-tree engine is an error rather than
## being ignored.  The boosted-tree engine's own options are
## documented under @code{fitrgam}.
##
## You can parse either a @qcode{'formula'} or an @qcode{'interactions'}
## optional parameter.  Parsing both parameters will result an error.
## Accordingly, you can only pass up to two parameters among @qcode{'knots'},
## @qcode{'order'}, and @qcode{'dof'} to define the required polynomial for
## training the GAM model.
##
## Two weak learners are available, selected by @code{FitMethod}.
##
## @qcode{'boostedtrees'}, the default, boosts one shallow decision tree per
## predictor in each round, which is the scheme MATLAB's generalized additive
## model uses.  A second phase then boosts trees over pairs of predictors,
## where interactions are asked for.
##
## @qcode{'splines'} boosts a smoothing spline per predictor until the
## residual sum of squares changes by less than @qcode{'Tol'}.  It has no
## MATLAB counterpart and is an Octave extension, kept because a smooth
## additive fit is a genuinely different and often better answer than a
## staircase of stumps.  A standard deviation and a prediction interval are
## available from it alone.
##
## The two take different arguments, and an argument meant for one is refused
## by the other rather than ignored.
##
## The choice is visible in the properties.  @code{Knots}, @code{Order},
## @code{DoF}, @code{Formula}, @code{Tol}, @code{BaseModel},
## @code{ModelwInt} and @code{IntMatrix} describe a spline fit and are empty
## under the boosted-tree engine, while @code{ModelParameters},
## @code{ReasonForTermination}, @code{BinEdges},
## @code{PairDetectionBinEdges} and @code{TreeModel} describe a tree fit and
## are empty under the spline engine.
##
## Fitted values are not expected to equal MATLAB's even under
## @qcode{'boostedtrees'}.  The stopping rule and the step-reduction limit are
## not recoverable from anything MATLAB reports, so this engine documents its
## own; what the two share is the estimator and the reported surface, not the
## arithmetic.
##
## @seealso{fitrgam, regress, regress_gp}
## @end deftypefn

  properties (GetAccess = public, SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix with one row per observation and one column per
    ## predictor of the training data.  This property is read-only.
    ##
    ## Where the model was fitted from a table, the predictors are the coded
    ## matrix and not the table: a variable holding levels is stored as its
    ## level codes, and the coding is kept with the model.
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
    ## the model was fitted on, rows with a missing response excluded.  This
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
    ## training, the @qcode{'Weights'} normalised to sum to one, and equal
    ## when none were given.  This property is read-only.
    ##
    ## @end deftp
    W                     = [];

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
    ## @deftp {RegressionGAM} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## A cell array with one entry per predictor, holding that predictor's bin
    ## edges where the model discretized it before fitting.  It is empty here
    ## and stays empty: this generalized additive model is built from splines,
    ## which take the predictors as they are, where MATLAB's is built from
    ## boosted trees and bins them.  That difference is described in the class
    ## documentation.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    BinEdges        = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} PairDetectionBinEdges
    ##
    ## Bin edges used to detect interactions
    ##
    ## A cell array with one row vector per predictor, holding the coarse cut
    ## points the residuals of the predictor phase were laid on while pairs
    ## were being tested.  The grid is eight equal-frequency bins whatever the
    ## sample size, as MATLAB's is.  It is empty when the model carries no
    ## interaction terms, and empty throughout under the spline engine, which
    ## does not bin.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    PairDetectionBinEdges = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} ModelParameters
    ##
    ## Parameters the model was fitted with
    ##
    ## A structure holding the fitting parameters.  Under the boosted-tree
    ## engine it carries MATLAB's own fields, with @qcode{Type} reading
    ## @qcode{'regression'}; under the spline engine it describes that scheme
    ## instead, since none of the tree vocabulary applies to it.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ModelParameters = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} ReasonForTermination
    ##
    ## Why each fitting phase stopped
    ##
    ## A structure with the fields @qcode{PredictorTrees} and
    ## @qcode{InteractionTrees}, each saying why that phase ended.  A phase
    ## that never ran reports an empty character vector.  It is empty under
    ## the spline engine, which has no tree budget to exhaust.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ReasonForTermination = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} FitMethod
    ##
    ## Which engine fitted the model
    ##
    ## A character vector, either @qcode{'boostedtrees'} or
    ## @qcode{'splines'}.  The default is @qcode{'boostedtrees'}, the scheme
    ## MATLAB's generalized additive model uses.  @qcode{'splines'} selects
    ## the penalised-spline engine, an Octave extension with no MATLAB
    ## counterpart and the scheme this class fitted before version 1.9.0.  The
    ## two engines take different arguments and an argument meant for one is
    ## refused by the other rather than ignored.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    FitMethod = 'boostedtrees';

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} TreeModel
    ##
    ## The fitted shape functions and interaction surfaces
    ##
    ## A structure with fields @qcode{ShapeValues}, @qcode{PairValues} and
    ## @qcode{Pairs}, holding what the boosted-tree engine fitted.  MATLAB
    ## exposes no equivalent, reporting its bin edges but never the values on
    ## them, so this is an Octave extension.  It is empty under the spline
    ## engine, whose fit lives in @code{BaseModel} and @code{ModelwInt}.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    TreeModel = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionGAM} {property} HyperparameterOptimizationResults
    ##
    ## Results of the hyperparameter optimization
    ##
    ## @strong{Always empty.}  It is declared for MATLAB compatibility, where
    ## it holds what an automatic search over the hyperparameters found.  This
    ## class fits the parameters it is given and runs no such search, so there
    ## is nothing to report.  This property is read-only.
    ##
    ## @end deftp
    HyperparameterOptimizationResults = [];

  endproperties

  ## Properties a user may set after the model is built.  Each one is
  ## validated by its set method below.
  properties (GetAccess = public, SetAccess = public)
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
  endproperties

  ## Readable by the counterpart class, which copies it, and kept out of
  ## the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    RTfun = @(y) y;

    ## How many trees each boosting phase actually fitted, which the budget
    ## in ModelParameters does not say: a phase may stop early.  Hidden
    ## because MATLAB reports it on the partitioned classes and not on the
    ## model, and that is where this is read from.
    NumTrainedTrees = [];
  endproperties

  ## Set methods for the properties a user may assign.
  methods (Hidden)

    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = parseResponseTransform ...
                                             (val, 'RegressionGAM');
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
      fprintf ("\n  RegressionGAM\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumPredictors', this.NumPredictors);
      fprintf ("%+25s: '%s'\n", 'ResponseTransform', this.ResponseTransform);
      fprintf ("%+25s: %g\n", 'Intercept', this.Intercept);
      str = repmat ({'%d'}, 1, numel (this.Knots));
      str = strcat ('[', strjoin (str, ' '), ']');
      fprintf ("%+25s: %s\n", 'Knots', sprintf (str, this.Knots));
      str = repmat ({'%d'}, 1, numel (this.Order));
      str = strcat ('[', strjoin (str, ' '), ']');
      fprintf ("%+25s: %s\n", 'Order', sprintf (str, this.Order));
      fprintf ("%+25s: %g\n", 'Tol', this.Tol);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGAM} {@var{obj} =} RegressionGAM (@var{X}, @var{Y})
    ## @deftypefnx {RegressionGAM} {@var{obj} =} RegressionGAM (@var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {RegressionGAM} {@var{obj} =} RegressionGAM (@var{Tbl}, @var{formula})
    ## @deftypefnx {RegressionGAM} {@var{obj} =} RegressionGAM (@var{Tbl}, @var{Y})
    ## @deftypefnx {RegressionGAM} {@var{obj} =} RegressionGAM (@dots{}, @var{name}, @var{value})
    ##
    ## Fit a generalized additive model for regression.
    ##
    ## @var{X} is an @math{N*P} numeric matrix of predictor data, one
    ## observation per row, and @var{Y} is the continuous response of those
    ## @math{N} observations.  The fit runs at construction, so @var{obj}
    ## arrives fitted.
    ##
    ## The @var{name}/@var{value} pairs the fit accepts, and the validation
    ## each one is held to, are listed in @code{help RegressionGAM}.
    ## @code{fitrgam} is the documented way to reach this constructor and
    ## takes the same pairs.
    ##
    ## @end deftypefn
    function this = RegressionGAM (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("RegressionGAM: too few input arguments.");
      endif

      ## A table names its own predictors and says which hold levels
      [this, X, Y, varargin] = resolveTable (this, 'RegressionGAM', ...
                                             X, Y, varargin);

      ## Get training sample size and number of variables in training data
      nsample = rows (X);
      ndims_X = columns (X);

      ## Check correspondence between predictors and response
      if (nsample != rows (Y))
        error ("RegressionGAM: number of rows in X and Y must be equal.");
      endif

      ## The response transform and the callable it names, unless one is
      ## given below
      ResponseTransform = 'none';
      RTfun             = @(y) y;

      ## Parse optional paired arguments
      optNames = {'PredictorNames', 'predictors', 'ResponseTransform', ...
                  'ResponseName', 'Formula', 'Interactions', 'Knots', ...
                  'Order', 'DoF', 'Tol', 'FitMethod', ...
                  'NumTreesPerPredictor', 'NumTreesPerInteraction', ...
                  'MaxNumSplitsPerPredictor', 'MaxNumSplitsPerInteraction', ...
                  'InitialLearnRateForPredictors', ...
                  'InitialLearnRateForInteractions', ...
                  'CategoricalPredictors', 'Weights', 'Verbose', 'NumPrint', ...
                  'MaxPValue'};
      ## An empty default stands for one resolved once the options are
      ## known, so that giving an option can be told apart from leaving it
      ## out: 'Knots', 'Order' and 'DoF' are 5, 3 and 8 per predictor, any
      ## two determining the third, and 'Tol' is 1e-3.  The boosted-tree
      ## options take MATLAB's own defaults, which ModelParameters reports and
      ## so are part of the surface being matched and not ours to improve:
      ## 300 and 100 trees, 1 and 4 splits, learning rates of 1, a
      ## 'MaxPValue' of 1, 'Verbose' 0 and 'NumPrint' 10.
      dfValues = {{}, [], [], [], [], [], [], [], [], [], 'boostedtrees', ...
                  [], [], [], [], [], [], [], [], [], [], []};
      [PredictorNames, PredictorsIn, RTin, ResponseName, Formula, ...
       Interactions, Knots, Order, DoF, Tol, FitMethod, ...
       NumTreesPerPredictor, NumTreesPerInteraction, ...
       MaxNumSplitsPerPredictor, MaxNumSplitsPerInteraction, ...
       InitialLearnRateForPredictors, InitialLearnRateForInteractions, ...
       CategoricalPredictors, Weights, Verbose, NumPrint, MaxPValue, args] = ...
                 parsePairedArguments (optNames, dfValues, varargin(:));

      ## Validate optional paired arguments
      if (! isempty (ResponseName) && ! ischar (ResponseName))
        error ("RegressionGAM: ResponseName must be a char string.");
      endif
      if (! isempty (Formula) && ! ischar (Formula) && ! islogical (Formula))
        error ("RegressionGAM: Formula must be a string.");
      endif
      if (! isempty (Interactions) &&
          ! ((isnumeric (Interactions) && isscalar (Interactions)
              && Interactions == fix (Interactions) && Interactions >= 0)
             || islogical (Interactions)
             || (ischar (Interactions) && strcmpi (Interactions, 'all'))))
        error ("RegressionGAM: invalid Interactions parameter.");
      endif
      if (! isempty (Knots)
          && (! isnumeric (Knots) || ! (isscalar (Knots) ||
              isequal (size (Knots), [1, ndims_X]))))
        error ("RegressionGAM: invalid value for Knots.");
      endif
      if (! isempty (Order)
          && (! isnumeric (Order) || ! (isscalar (Order) ||
              isequal (size (Order), [1, ndims_X]))))
        error ("RegressionGAM: invalid value for Order.");
      endif
      if (! isempty (DoF)
          && (! isnumeric (DoF) ||
              ! (isscalar (DoF) || isequal (size (DoF), [1, ndims_X]))))
        error ("RegressionGAM: invalid value for DoF.");
      endif
      if (! isempty (Tol) &&
          ! (isnumeric (Tol) && isscalar (Tol) && (Tol > 0)))
        error ("RegressionGAM: Tolerance must be a Positive scalar.");
      endif
      if (! (ischar (FitMethod) && isrow (FitMethod)) ||
          ! any (strcmpi (FitMethod, {'boostedtrees', 'splines'})))
        error (strcat ("RegressionGAM: 'FitMethod' must be", ...
                       " 'boostedtrees' or 'splines'."));
      endif
      FitMethod = tolower (FitMethod);
      if (! isempty (NumTreesPerPredictor)
          && (! isnumeric (NumTreesPerPredictor) ||
              ! isscalar (NumTreesPerPredictor) ||
              NumTreesPerPredictor < 1 ||
              fix (NumTreesPerPredictor) != NumTreesPerPredictor))
        error (strcat ("RegressionGAM: 'NumTreesPerPredictor'", ...
                       " must be a positive integer value."));
      endif
      if (! isempty (NumTreesPerInteraction)
          && (! isnumeric (NumTreesPerInteraction) ||
              ! isscalar (NumTreesPerInteraction) ||
              NumTreesPerInteraction < 1 ||
              fix (NumTreesPerInteraction) != NumTreesPerInteraction))
        error (strcat ("RegressionGAM: 'NumTreesPerInteraction'", ...
                       " must be a positive integer value."));
      endif
      if (! isempty (MaxNumSplitsPerPredictor)
          && (! isnumeric (MaxNumSplitsPerPredictor) ||
              ! isscalar (MaxNumSplitsPerPredictor) ||
              MaxNumSplitsPerPredictor < 1 ||
              fix (MaxNumSplitsPerPredictor) != MaxNumSplitsPerPredictor))
        error (strcat ("RegressionGAM:", ...
                       " 'MaxNumSplitsPerPredictor' must be a", ...
                       " positive integer value."));
      endif
      if (! isempty (MaxNumSplitsPerInteraction)
          && (! isnumeric (MaxNumSplitsPerInteraction) ||
              ! isscalar (MaxNumSplitsPerInteraction) ||
              MaxNumSplitsPerInteraction < 1 ||
              fix (MaxNumSplitsPerInteraction) != MaxNumSplitsPerInteraction))
        error (strcat ("RegressionGAM:", ...
                       " 'MaxNumSplitsPerInteraction' must be a", ...
                       " positive integer value."));
      endif
      if (! isempty (InitialLearnRateForPredictors)
          && (! isnumeric (InitialLearnRateForPredictors) ||
              ! isscalar (InitialLearnRateForPredictors) ||
              InitialLearnRateForPredictors <= 0 ||
              InitialLearnRateForPredictors > 1))
        error (strcat ("RegressionGAM:", ...
                       " 'InitialLearnRateForPredictors' must be", ...
                       " greater than 0 and at most 1."));
      endif
      if (! isempty (InitialLearnRateForInteractions)
          && (! isnumeric (InitialLearnRateForInteractions) ||
              ! isscalar (InitialLearnRateForInteractions) ||
              InitialLearnRateForInteractions <= 0 ||
              InitialLearnRateForInteractions > 1))
        error (strcat ("RegressionGAM:", ...
                       " 'InitialLearnRateForInteractions' must be", ...
                       " greater than 0 and at most 1."));
      endif
      if (! isempty (Weights) &&
          ! (isnumeric (Weights) && isreal (Weights)
             && isvector (Weights)))
        error ("RegressionGAM: 'Weights' must be a numeric vector.");
      endif
      if (! isempty (Weights) && numel (Weights) != rows (X))
        error (strcat ("RegressionGAM: 'Weights' must have one", ...
                       " element per row of X."));
      endif
      if (! isempty (Weights)
          && (any (Weights(:) < 0) || ! all (isfinite (Weights(:)))))
        error (strcat ("RegressionGAM: 'Weights' must hold", ...
                       " finite non-negative values."));
      endif
      if (! isempty (Verbose)
          && (! isnumeric (Verbose) || ! isscalar (Verbose) || Verbose < 0
              || fix (Verbose) != Verbose))
        error (strcat ("RegressionGAM: 'Verbose' must be a", ...
                       " non-negative integer value."));
      endif
      if (! isempty (NumPrint)
          && (! isnumeric (NumPrint) || ! isscalar (NumPrint)
              || NumPrint < 1 || fix (NumPrint) != NumPrint))
        error (strcat ("RegressionGAM: 'NumPrint' must be a positive", ...
                       " integer value."));
      endif
      if (! isempty (MaxPValue)
          && (! isnumeric (MaxPValue) || ! isscalar (MaxPValue) ||
              MaxPValue < 0 || MaxPValue > 1))
        error (strcat ("RegressionGAM: 'MaxPValue' must be", ...
                       " between 0 and 1."));
      endif

      ## 'Predictors' is another name for 'PredictorNames'
      if (isempty (PredictorNames))
        PredictorNames = PredictorsIn;
      endif
      if (! isempty (PredictorNames))
        if (! iscellstr (PredictorNames))
          error (strcat ("RegressionGAM: PredictorNames must be a", ...
                         " cellstring array."));
        elseif (columns (PredictorNames) != columns (X))
          error (strcat ("RegressionGAM: PredictorNames must have same", ...
                         " number of columns as X."));
        endif
      endif
      if (! isempty (RTin))
        [RTfun, ResponseTransform] = parseResponseTransform (RTin, ...
                                                             'RegressionGAM');
      endif

      if (! isempty (args))
        error ("RegressionGAM: invalid optional paired argument.");
      endif

      ## An argument belongs to one engine or the other, and asking for one
      ## the chosen engine cannot honour is refused rather than ignored.
      ## An option was given when it is not empty, its default being filled
      ## in below.
      splineOnly = {'knots', Knots; 'order', Order; 'dof', DoF; ...
                    'formula', Formula; 'tol', Tol};
      treeOnly = {'numtreesperpredictor', NumTreesPerPredictor; ...
                  'numtreesperinteraction', NumTreesPerInteraction; ...
                  'maxnumsplitsperpredictor', MaxNumSplitsPerPredictor; ...
                  'maxnumsplitsperinteraction', MaxNumSplitsPerInteraction; ...
                  'initiallearnrateforpredictors', ...
                  InitialLearnRateForPredictors; ...
                  'initiallearnrateforinteractions', ...
                  InitialLearnRateForInteractions; ...
                  'maxpvalue', MaxPValue; 'verbose', Verbose; ...
                  'numprint', NumPrint; 'weights', Weights; ...
                  'categoricalpredictors', CategoricalPredictors};
      if (strcmp (FitMethod, 'boostedtrees'))
        clash = splineOnly(! cellfun (@isempty, splineOnly(:,2)), 1);
        if (! isempty (clash))
          error (strcat ("RegressionGAM: '", clash{1}, "' is a parameter", ...
                         " of the spline engine and cannot be used with", ...
                         " 'FitMethod' 'boostedtrees'."));
        endif
      else
        clash = treeOnly(! cellfun (@isempty, treeOnly(:,2)), 1);
        if (! isempty (clash))
          error (strcat ("RegressionGAM: '", clash{1}, "' is a parameter", ...
                         " of the boosted-tree engine and cannot be used", ...
                         " with 'FitMethod' 'splines'."));
        endif
      endif

      ## A model is described either by a formula or by its interactions
      if (! isempty (Formula) && ! isempty (Interactions))
        error (strcat ("RegressionGAM: 'Formula' and 'Interactions'", ...
                       " cannot be given together."));
      endif

      ## Any two of 'Knots', 'Order' and 'DoF' determine the third, DoF being
      ## Knots plus Order; one given alone keeps the default of another.
      if (! isempty (Knots) && ! isempty (Order) && ! isempty (DoF))
        error (strcat ("RegressionGAM: at most two of 'Knots', 'Order'", ...
                       " and 'DoF' may be given."));
      endif
      if (isempty (DoF))
        if (isempty (Knots))
          Knots = ones (1, ndims_X) * 5;
        endif
        if (isempty (Order))
          Order = ones (1, ndims_X) * 3;
        endif
        DoF = Knots + Order;
      elseif (isempty (Knots))
        if (isempty (Order))
          Order = ones (1, ndims_X) * 3;
        endif
        Knots = DoF - Order;
      else
        Order = DoF - Knots;
      endif

      ## The defaults of the remaining engine options
      if (isempty (Tol))
        Tol = 1e-3;
      endif
      if (isempty (NumTreesPerPredictor))
        NumTreesPerPredictor = 300;
      endif
      if (isempty (NumTreesPerInteraction))
        NumTreesPerInteraction = 100;
      endif
      if (isempty (MaxNumSplitsPerPredictor))
        MaxNumSplitsPerPredictor = 1;
      endif
      if (isempty (MaxNumSplitsPerInteraction))
        MaxNumSplitsPerInteraction = 4;
      endif
      if (isempty (InitialLearnRateForPredictors))
        InitialLearnRateForPredictors = 1;
      endif
      if (isempty (InitialLearnRateForInteractions))
        InitialLearnRateForInteractions = 1;
      endif
      if (isempty (MaxPValue))
        MaxPValue = 1;
      endif
      if (isempty (Verbose))
        Verbose = 0;
      endif
      if (isempty (NumPrint))
        NumPrint = 10;
      endif

      ## Assign original X and Y data to the RegressionGAM object
      this.X = X;
      this.Y = Y;

      ## An observation is dropped only when its response is missing.  A row
      ## whose predictors hold missing values is kept and reported as used,
      ## and the boosted-tree fit uses it, as R2024a does, the row sitting in
      ## the root of each tree it cannot be split by; the spline fit below
      ## draws on the complete observations alone.
      RowsUsed  = ! isnan (Y(:));
      Yret      = Y(RowsUsed);
      Xret      = X(RowsUsed, :);
      if (isempty (Weights))
        Wret    = ones (rows (Xret), 1);
      else
        Wret    = Weights(RowsUsed);
        Wret    = Wret(:);
      endif
      this.X    = Xret;
      this.Y    = Yret;
      if (strcmp (FitMethod, 'boostedtrees'))
        cobs    = true (rows (Xret), 1);
      else
        cobs    = ! any (isnan (Xret), 2);
      endif
      Y         = Yret(cobs);
      X         = Xret(cobs, :);

      ## Check X and Y contain valid data
      if (! isnumeric (X) || ! all (isfinite (X(:)) | isnan (X(:))))
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
      this.RTfun            = RTfun;

      ## Bookkeeping MATLAB reports alongside the fit
      ## A categorical predictor is fitted on the indices of its levels,
      ## which the engine is handed in place of the values.
      this.CategoricalPredictors = [];
      if (! isempty (CategoricalPredictors))
        [CP, errmsg] = dummyCoding (X, CategoricalPredictors, ...
                                    this.PredictorNames);
        if (! isempty (errmsg))
          error ("RegressionGAM: %s", errmsg);
        endif
        this.CategoricalPredictors = CP.Index;
        Levels = cell (1, columns (X));
        Levels(CP.Index) = CP.Levels(CP.Index);
        this.TreeModel = struct ('CategoricalLevels', {Levels});
      endif
      this.ExpandedPredictorNames = PredictorNames;
      ## The weights over their sum, as MATLAB reports them; the boosted fit
      ## sees them on the rows it is fitted to.
      this.W = Wret / sum (Wret);
      this.IsStandardDeviationFit = false;

      this.FitMethod = FitMethod;

      if (strcmp (FitMethod, 'boostedtrees'))
        ## The spline parameters describe a scheme that did not run, so they
        ## are left empty rather than reporting numbers nothing used.
        this.Knots = [];
        this.Order = [];
        this.DoF   = [];
        this.Tol   = [];
        this = this.fitBoosted (X, Y, Interactions, ...
                                NumTreesPerPredictor, ...
                                NumTreesPerInteraction, ...
                                MaxNumSplitsPerPredictor, ...
                                MaxNumSplitsPerInteraction, ...
                                InitialLearnRateForPredictors, ...
                                InitialLearnRateForInteractions, MaxPValue, ...
                                Verbose, NumPrint);
        return;
      endif

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
      if (! isempty (Formula) || ! isempty (Interactions))
        this = this.fitModelwInt (X, Y, Inter, Knots, Order, DoF);
      endif

      ## The property MATLAB reports is the two-way terms the fitted model
      ## carries, as predictor index pairs, whatever form they were asked for
      ## in.  The term matrix stays the complete record: it also holds the
      ## main effects a formula names and any term above two predictors,
      ## neither of which has a two-column form.
      this.Interactions = interactionPairs (this.IntMatrix);

      ## The spline scheme has no tree vocabulary to report, so its parameter
      ## struct describes itself instead.
      this.ModelParameters = struct ('Knots', this.Knots, ...
                                     'Order', this.Order, ...
                                     'DoF', this.DoF, ...
                                     'Formula', this.Formula, ...
                                     'Interactions', this.Interactions, ...
                                     'Tol', this.Tol);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGAM} {@var{obj} =} addInteractions (@var{obj}, @var{interactions})
    ##
    ## Add interaction terms to a fitted model.
    ##
    ## @code{@var{obj} = addInteractions (@var{obj}, @var{interactions})} fits
    ## the interaction terms named by @var{interactions} on top of the terms
    ## the model already carries and returns the updated model.  The univariate
    ## fit is left alone, so @code{predict} with
    ## @qcode{'IncludeInteractions'} set @qcode{false} answers exactly as it
    ## answered before.
    ##
    ## @var{interactions} takes the forms the constructor's
    ## @qcode{'Interactions'} option takes: a nonnegative integer count of
    ## terms, a logical matrix with a column per predictor, or @qcode{'all'}.
    ##
    ## A model already carrying interaction terms is not extended, which is
    ## what MATLAB refuses too.  A model fitted from a @qcode{'Formula'} names
    ## every term it has, interactions among them, and is refused for the same
    ## reason.
    ##
    ## Which terms a count selects is this implementation's own: they are
    ## taken in the order @code{nchoosek} lists the pairs, where MATLAB ranks
    ## them by how much each contributes.  The constructor's option chooses
    ## the same way, so the two agree with each other.
    ##
    ## @seealso{fitrgam, RegressionGAM}
    ## @end deftypefn
    function this = addInteractions (this, interactions)

      if (nargin != 2)
        print_usage ();
      endif

      ## Which store already holds interaction terms depends on the engine.
      hasInt = ! isempty (this.IntMatrix);
      if (strcmp (this.FitMethod, 'boostedtrees') && ! isempty (this.TreeModel))
        hasInt = ! isempty (this.TreeModel.Pairs);
      endif
      if (hasInt)
        error (strcat ("RegressionGAM.addInteractions: adding interaction", ...
                       " terms to a model that already includes them is", ...
                       " not supported."));
      endif

      if (! ((isnumeric (interactions) && isscalar (interactions)
              && interactions == fix (interactions) && interactions >= 0)
             || islogical (interactions)
             || (ischar (interactions) && strcmpi (interactions, 'all'))))
        error (strcat ("RegressionGAM.addInteractions: invalid", ...
                       " 'Interactions' parameter."));
      endif

      ## Under the boosted-tree engine the interaction phase simply runs now,
      ## starting from the predictor phase the model already carries, so a
      ## model with interactions added is the model the constructor would
      ## have built had it been asked for them.
      if (strcmp (this.FitMethod, 'boostedtrees'))
        cobs = true (rows (this.X), 1);
        [Xfit, E, cat] = gamCatCode (this.TreeModel, this.X(cobs, :), ...
                                     this.BinEdges);
        Yfit = this.Y(cobs);
        MP = this.ModelParameters;
        f = gamboostpredict (E, this.TreeModel.ShapeValues, ...
                             Xfit, this.Intercept);
        res = Yfit - f;

        wanted = -1;
        pairs = zeros (0, 2);
        if (ischar (interactions))
          wanted = Inf;
        elseif (isscalar (interactions) && ! isempty (interactions))
          wanted = interactions;
        elseif (! isempty (interactions))
          pairs = interactionPairs (logical (interactions));
        endif

        if (wanted > 0 && columns (Xfit) > 1)
          S = gamboostpairs (Xfit, res, cat);
          pval = 1 - fcdf (S.F, S.DF1, S.DF2);
          pval(S.DF1 <= 0) = 1;
          [pval, ord] = sort (pval);
          ranked = S.Pairs(ord, :);
          ranked = ranked(pval <= MP.MaxPValue, :);
          if (isfinite (wanted) && rows (ranked) > wanted)
            ranked = ranked(1:wanted, :);
          endif
          pairs = ranked;
          ## As in fitBoosted: a single split is no interaction.
          if (MP.MaxNumSplitsPerInteraction < 2)
            pairs = zeros (0, 2);
          endif
          this.PairDetectionBinEdges = S.BinEdges(:);
          this.PairDetectionBinEdges(cat) = {[]};
        endif

        reason = this.ReasonForTermination;
        ntrees = this.NumTrainedTrees;
        if (! isempty (pairs))
          I = gamboostinter (Xfit, Yfit, f, 2, pairs, ...
                             MP.NumTreesPerInteraction, ...
                             MP.InitialLearnRateForInteractions, ...
                             MP.MaxNumSplitsPerInteraction, this.W(cobs), cat);
          ## A pair tree splitting on one predictor alone is a main effect and
          ## adds nothing, so a fit made only of such trees keeps no pair, as
          ## R2024a keeps none.
          if (I.NumTrees == 0)
            pairs = zeros (0, 2);
            warning (strcat ("RegressionGAM: model does not include", ...
                             " interaction terms because all interaction", ...
                             " terms have p-values greater than the", ...
                             " 'MaxPValue' value, or the software was", ...
                             " unable to improve the model fit."));
          endif
        endif
        if (! isempty (pairs))
          this.Intercept = this.Intercept + I.Intercept;
          if (wanted <= 0)
            this.PairDetectionBinEdges = I.PairBinEdges(:);
            this.PairDetectionBinEdges(cat) = {[]};
          endif
          this.TreeModel.PairValues = I.PairValues;
          this.TreeModel.PairEdges = I.PairEdges;
          this.TreeModel.PairMissing = I.PairMissing;
          this.TreeModel.PairIntercept = I.Intercept;
          reason.InteractionTrees = I.ReasonForTermination;
          ntrees.InteractionTrees = I.NumTrees;
        endif
        this.TreeModel.Pairs = pairs;
      if (isempty (pairs))
        this.TreeModel.PairIntercept = 0;
      endif
        this.Interactions = pairs;
        this.ReasonForTermination = reason;
        this.NumTrainedTrees = ntrees;
        MP.Interactions = interactions;
        this.ModelParameters = MP;
        return;
      endif

      ## parseInteractions reads the specification from the property, which
      ## afterwards holds the pairs the fit settled on, exactly as the
      ## constructor leaves it.
      this.Interactions = interactions;
      this.IntMatrix = this.parseInteractions ();

      ## The fit sees the complete observations, prepared as the constructor
      ## prepares them.  Knots, Order and DoF are held unexpanded on the
      ## object and are widened to the interaction columns by fitModelwInt.
      cobs = ! any (isnan (this.X), 2);
      Xfit = this.X(cobs, :);
      Yfit = this.Y(cobs);
      this = this.fitModelwInt (Xfit, Yfit, mean (Yfit), this.Knots, ...
                                this.Order, this.DoF);
      this.Interactions = interactionPairs (this.IntMatrix);

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
    ## training data in @var{obj}.  Every row is predicted.  Under boosted
    ## trees a missing value adds nothing from a main effect, and an
    ## interaction term takes the value its trees give a row missing that
    ## predictor, so a row missing every predictor predicts the intercept;
    ## under splines a row holding a missing value is predicted as
    ## @code{NaN}.
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
    ##
    ## The new data may be a table, whose variables are matched to the
    ## predictors the model was fitted on by name and not by position: one
    ## the model was not fitted on is passed over, one it needs and cannot
    ## find is named, and a value holding a level is coded as that level
    ## was coded at fitting.
    ## @seealso{fitrgam, RegressionGAM}
    ## @end deftypefn
    function [yFit, ySD, yInt] = predict (this, Xfit, varargin)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("RegressionGAM.predict: too few arguments.");
      endif

      ## A table is read by the names the model was fitted on
      Xfit = tableColumns (this, 'RegressionGAM.predict', Xfit);

      ## Check for valid XC
      if (isempty (Xfit))
        error ("RegressionGAM.predict: Xfit is empty.");
      elseif (columns (this.X) != columns (Xfit))
        error (strcat ("@RegressionGAM/predict: Xfit must have the same", ...
                       " number of features (columns) as in the GAM model."));
      endif

      ## Default values for Name-Value Pairs
      alpha = 0.05;
      ## Which store holds the interaction terms depends on the engine: the
      ## spline scheme keeps them as extra columns described by IntMatrix,
      ## the boosted-tree scheme as surfaces over predictor pairs.
      hasInt = ! isempty (this.IntMatrix);
      if (strcmp (this.FitMethod, 'boostedtrees') && ! isempty (this.TreeModel))
        hasInt = ! isempty (this.TreeModel.Pairs);
      endif
      if (! hasInt)
        incInt = false;
      else
        incInt = true;
      endif

      ## Parse optional paired arguments; interactions are included when the
      ## model has them
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionGAM.predict: optional arguments must be", ...
                       " given in Name-Value pairs."));
      endif
      optNames = {'IncludeInteractions', 'Alpha'};
      dfValues = {incInt, alpha};
      [incInt, alpha, args] = ...
                 parsePairedArguments (optNames, dfValues, varargin(:));

      ## Validate optional paired arguments
      if (! islogical (incInt) || (incInt != 0 && incInt != 1))
        error (strcat ("RegressionGAM.predict: includeinteractions must be", ...
                       " a logical value."));
      endif
      if (incInt && ! hasInt)
        error (strcat ("RegressionGAM.predict: trained model does not", ...
                       " include any interactions."));
      endif
      if (! (isnumeric (alpha) && isscalar (alpha) && alpha > 0 && alpha < 1))
        error (strcat ("RegressionGAM.predict: alpha must be a scalar", ...
                       " value between 0 and 1."));
      endif

      if (! isempty (args))
        error ("RegressionGAM.predict: invalid optional paired argument.");
      endif

      ## Choose whether interactions must be included.  The reshaping is done
      ## by gamTerms rather than inline, because the training data has to be
      ## reshaped exactly the same way further down: a model built with
      ## interactions or with a formula has a term for every column of the
      ## matrix it was fitted on, and that is no longer the stored X.
      ## The boosted-tree engine keeps its fit as step functions over bins, so
      ## a term is a lookup rather than a spline evaluation and the whole
      ## prediction is one call.  Everything after it is shared.
      if (strcmp (this.FitMethod, 'boostedtrees') && ! isempty (this.TreeModel))
        ## Excluding the interactions means excluding the constant they
        ## handed the intercept as well.
        [Xb, E] = gamCatCode (this.TreeModel, Xfit, this.BinEdges);
        interc = this.Intercept;
        if (! incInt && isfield (this.TreeModel, 'PairIntercept'))
          interc = interc - this.TreeModel.PairIntercept;
        endif
        if (! incInt || isempty (this.TreeModel.Pairs))
          yFit = gamboostpredict (E, ...
                                  this.TreeModel.ShapeValues, Xb, ...
                                  interc);
        else
          [PE, PM] = gamPairEdges (this.TreeModel, this.PairDetectionBinEdges);
          yFit = gamboostpredict (E, ...
                                  this.TreeModel.ShapeValues, Xb, ...
                                  interc, 0, PE, ...
                                  this.TreeModel.PairValues, ...
                                  this.TreeModel.Pairs, PM);
        endif
        yFit = this.RTfun (yFit);
        if (nargout > 1)
          error (strcat ("RegressionGAM.predict: a standard deviation is", ...
                         " only available from a model fitted with", ...
                         " 'FitMethod' 'splines'."));
        endif
        return;
      endif

      if (incInt)
        Xfit = gamTerms (Xfit, this.IntMatrix, isempty (this.Formula));
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
      yFit = this.RTfun (yFit);

      ## Predict Standard Deviation and Intervals of estimated data if requested
      if (nargout > 1)
        ## Ensure that RowsUsed in the model are selected
        used = true (rows (this.X), 1);
        Y = this.Y(used);
        X = this.X(used, :);
        ## Reshape the training data as the prediction data was reshaped.  It
        ## used to be passed as stored, so the terms and the columns disagreed
        ## whenever the model was built with interactions or with a formula:
        ## with more terms than columns the residuals came from a truncated
        ## model and ySD was several times too large, and with fewer terms than
        ## columns the call raised out of bound.
        if (incInt)
          X = gamTerms (X, this.IntMatrix, isempty (this.Formula));
        endif
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
          lower  = this.RTfun (yrs_fit - moe);
          upper  = this.RTfun (yrs_fit + moe);
          yInt   = [lower, upper];
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGAM} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionGAM} {@var{L} =} loss (@var{obj}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {RegressionGAM} {@var{L} =} loss (@var{obj}, @var{Tbl})
    ## @deftypefnx {RegressionGAM} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss of a generalized additive model.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the weighted
    ## mean squared error of the model on the rows of @var{X} against the true
    ## response @var{Y}.
    ##
    ## @var{X} may also be a table @var{Tbl}, whose variables are matched to
    ## the predictors the model was fitted on by name and not by position.
    ## @code{loss (@var{obj}, @var{Tbl}, @var{ResponseVarName})} takes the
    ## response from the variable @var{ResponseVarName} names, and
    ## @code{loss (@var{obj}, @var{Tbl})} from the variable the model was
    ## fitted on.  The response may also be given beside the table as
    ## @var{Y}.
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
    ## @seealso{RegressionGAM, fitrgam, RegressionGAM.predict}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3 && ! (nargin > 1 && istable (X)))
        error ("RegressionGAM.loss: too few input arguments.");
      endif

      ## A table carries the response: named in the call, given beside
      ## the table, or the variable the model was fitted on
      if (nargin < 3)
        Y = [];
      endif
      [X, Y, varargin] = tableResponse (this, 'loss', X, Y, varargin, ...
                                        nargin > 2);
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionGAM.loss: Name-Value arguments must", ...
                       " be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, 'loss');

      ## Parse optional paired arguments; an empty 'Weights' stands for
      ## uniform weights
      optNames = {'LossFun', 'Weights'};
      dfValues = {'mse', []};
      [LossFun, W, args] = ...
                 parsePairedArguments (optNames, dfValues, varargin(:));

      ## Validate optional paired arguments
      if (! (is_function_handle (LossFun) ||
             (ischar (LossFun) && isrow (LossFun))))
        error (strcat ("RegressionGAM.loss: 'LossFun' must be a character", ...
                       " vector or a function handle."));
      endif
      if (ischar (LossFun) && ! strcmpi (LossFun, 'mse'))
        error ("RegressionGAM.loss: unsupported 'LossFun' value.");
      endif
      if (! isempty (W) && ! (isnumeric (W) && isvector (W)))
        error ("RegressionGAM.loss: 'Weights' must be a numeric vector.");
      endif
      if (! isempty (W) && numel (W) != rows (X))
        error (strcat ("RegressionGAM.loss: size of 'Weights' must equal", ...
                       " the number of rows in X."));
      endif

      if (! isempty (args))
        error ("RegressionGAM.loss: invalid optional paired argument.");
      endif
      if (isempty (W))
        W = ones (rows (X), 1);
      endif

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
    ## @seealso{RegressionGAM, fitrgam, RegressionGAM.predict}
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
    ## @seealso{RegressionGAM, fitrgam, RegressionGAM.loss}
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      used = true (rows (this.X), 1);
      L = loss (this, this.X(used, :), this.Y(used), varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionGAM} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {RegressionGAM} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross validate a Generalized Additive Model regression object.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} returns a cross-validated
    ## model object, @var{CVMdl}, from a trained model, @var{obj}, using
    ## 10-fold cross-validation by default.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj}, @var{name}, @var{value})}
    ## specifies additional name-value pair arguments to customize the
    ## cross-validation process.
    ##
    ## @multitable @columnfractions 0.28 0.7
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'KFold'} @tab Specify the number of folds to use in
    ## k-fold cross-validation.  @code{"KFold", @var{k}}, where @var{k} is an
    ## integer greater than 1.
    ##
    ## @item @qcode{'Holdout'} @tab Specify the fraction of the data to
    ## hold out for testing.  @code{"Holdout", @var{p}}, where @var{p} is a
    ## scalar in the range @math{(0,1)}.
    ##
    ## @item @qcode{'Leaveout'} @tab Specify whether to perform
    ## leave-one-out cross-validation.  @code{"Leaveout", @var{Value}}, where
    ## @var{Value} is 'on' or 'off'.
    ##
    ## @item @qcode{'CVPartition'} @tab Specify a @qcode{cvpartition}
    ## object used for cross-validation.  @code{"CVPartition", @var{cv}},
    ## where @code{isa (@var{cv}, "cvpartition")} = 1.
    ##
    ## @end multitable
    ##
    ## @seealso{fitrgam, RegressionGAM, cvpartition,
    ## RegressionPartitionedModel}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)

      if (numel (varargin) == 1)
        error (strcat ("RegressionGAM.crossval: Name-Value", ...
                       " arguments must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("RegressionGAM.crossval: specify only", ...
                       " one of the optional Name-Value paired arguments."));
      endif

      if (this.NumObservations < 10)
        numFolds  = this.NumObservations;
      else
        numFolds  = 10;
      endif
      Holdout     = [];
      Leaveout    = 'off';
      CVPartition = [];

      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'kfold'
            numFolds = varargin{2};
            if (! (isnumeric (numFolds) && isscalar (numFolds)
                   && (numFolds == fix (numFolds)) && numFolds > 1))
              error (strcat ("RegressionGAM.crossval: 'KFold'", ...
                             " must be an integer value greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("RegressionGAM.crossval: 'Holdout'", ...
                             " must be a numeric value between 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat ("RegressionGAM.crossval: 'Leaveout'", ...
                             " must be either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("RegressionGAM.crossval: 'CVPartition'", ...
                             " must be a 'cvpartition' object."));
            endif

          otherwise
            error (strcat ("RegressionGAM.crossval: invalid", ...
                           " parameter name in optional paired arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Determine the cross-validation method to use.  The partition is
      ## built over the observations actually trained on, so its indices and
      ## the partitioned model's rows are the same set.
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

      ## Create a cross-validated model object
      CVMdl = RegressionPartitionedModel (this, partition);

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
      BinEdges            = obj.BinEdges;
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
      RTfun                 = obj.RTfun;
      FitMethod             = obj.FitMethod;
      TreeModel             = obj.TreeModel;
      PairDetectionBinEdges = obj.PairDetectionBinEdges;
      ModelParameters       = obj.ModelParameters;
      ReasonForTermination  = obj.ReasonForTermination;

      ## Save classdef name and all model properties as individual variables
      HyperparameterOptimizationResults = obj.HyperparameterOptimizationResults;
      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'NumObservations', ...
            'RowsUsed', 'BinEdges', 'NumPredictors', 'PredictorNames', ...
            'ResponseName', ...
            'Formula', 'Interactions', 'Knots', 'Order', 'DoF', 'Tol', ...
            'BaseModel', 'ModelwInt', 'IntMatrix', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'W', 'ResponseTransform', ...
            'Intercept', 'IsStandardDeviationFit', 'RTfun', ...
            'FitMethod', 'TreeModel', 'PairDetectionBinEdges', ...
            'ModelParameters', 'ReasonForTermination', ...
            'HyperparameterOptimizationResults');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionGAM} {@var{Mdl} =} resume (@var{obj}, @var{numTrees})
    ##
    ## Resume training a generalized additive model.
    ##
    ## @code{@var{Mdl} = resume (@var{obj}, @var{numTrees})} adds
    ## @var{numTrees} more trees to @var{obj} and returns the result.  The
    ## original model is not modified.
    ##
    ## Training continues in the phase that ran last, which is what MATLAB
    ## does: a model carrying interaction terms gains interaction trees and
    ## its predictor shape functions are left alone, while a model without
    ## them gains predictor trees.  A round starts at its initial learning
    ## rate whatever its number, so the model this returns is the model a
    ## single fit of the combined budget would have produced.
    ##
    ## @var{numTrees} must be a positive integer scalar.  Resuming raises
    ## where there is nothing left to gain, rather than returning the model
    ## unchanged, and it is not available under
    ## @qcode{'FitMethod', 'splines'}: a backfit that has converged to its
    ## tolerance has no budget to extend.
    ##
    ## @seealso{RegressionGAM, fitrgam, RegressionGAM.addInteractions}
    ## @end deftypefn
    function Mdl = resume (this, numTrees)

      if (nargin < 2)
        error ("RegressionGAM.resume: Not enough input arguments.");
      endif
      if (! strcmp (this.FitMethod, 'boostedtrees'))
        error (strcat ("RegressionGAM.resume: resuming is available", ...
                       " only under 'FitMethod', 'boostedtrees'; a spline", ...
                       " backfit stops at its tolerance and has no budget", ...
                       " to extend."));
      endif
      if (! (isnumeric (numTrees) && isscalar (numTrees) && isreal (numTrees)
             && numTrees > 0 && fix (numTrees) == numTrees))
        error (strcat ("RegressionGAM.resume: NUMTREES must be a", ...
                       " positive integer scalar."));
      endif

      ## The rows the fit saw, which are all of them.
      cobs = true (rows (this.X), 1);
      [X, E, cat] = gamCatCode (this.TreeModel, this.X(cobs, :), ...
                               this.BinEdges);
      Y = this.Y(cobs);

      Mdl = this;
      MP = this.ModelParameters;
      reason = this.ReasonForTermination;
      ntrees = this.NumTrainedTrees;

      if (isempty (this.TreeModel.Pairs))
        ## No interaction phase ever ran, so the predictor phase is the one
        ## still open.  The engine is handed the prediction reached so far and
        ## returns the increment to add to it.
        f = gamboostpredict (E, this.TreeModel.ShapeValues, X, ...
                             this.Intercept);
        M = gamboosttrain (X, Y, 2, numTrees, ...
                           MP.InitialLearnRateForPredictors, ...
                           MP.MaxNumSplitsPerPredictor, 0, MP.NumPrint, ...
                           f(:), this.W(cobs), cat);
        if (M.NumTrees == 0)
          error (strcat ("RegressionGAM.resume: unable to resume", ...
                         " training because the software was unable to", ...
                         " improve the model fit."));
        endif
        sv = this.TreeModel.ShapeValues;
        for j = 1:numel (sv)
          sv{j} = sv{j} + M.ShapeValues{j};
        endfor
        Mdl.TreeModel.ShapeValues = sv;
        Mdl.Intercept = this.Intercept + M.Intercept;
        MP.NumTreesPerPredictor = MP.NumTreesPerPredictor + M.NumTrees;
        reason.PredictorTrees = M.ReasonForTermination;
        ntrees.PredictorTrees = ntrees.PredictorTrees + M.NumTrees;
      else
        ## The interaction phase ran last, so it is the one extended.  The
        ## running prediction includes the surfaces already fitted, and the
        ## new ones are added to them.
        [PE, PM] = gamPairEdges (this.TreeModel, this.PairDetectionBinEdges);
        f = gamboostpredict (E, this.TreeModel.ShapeValues, X, ...
                             this.Intercept, 0, PE, ...
                             this.TreeModel.PairValues, ...
                             this.TreeModel.Pairs, PM);
        I = gamboostinter (X, Y, f(:), 2, this.TreeModel.Pairs, numTrees, ...
                           MP.InitialLearnRateForInteractions, ...
                           MP.MaxNumSplitsPerInteraction, this.W(cobs), cat);
        if (I.NumTrees == 0)
          error (strcat ("RegressionGAM.resume: unable to resume", ...
                         " training because the software was unable to", ...
                         " improve the model fit."));
        endif
        ## The new trees cut where they cut, so the surfaces are added on the
        ## union of both grids.
        [pe, pv, pm] = gamPairMerge (PE, this.TreeModel.PairValues, PM, ...
                                     I.PairEdges, I.PairValues, ...
                                     I.PairMissing);
        Mdl.TreeModel.PairEdges = pe;
        Mdl.TreeModel.PairValues = pv;
        Mdl.TreeModel.PairMissing = pm;
        Mdl.TreeModel.PairIntercept = this.TreeModel.PairIntercept ...
                                      + I.Intercept;
        Mdl.Intercept = this.Intercept + I.Intercept;
        MP.NumTreesPerInteraction = MP.NumTreesPerInteraction + I.NumTrees;
        reason.InteractionTrees = I.ReasonForTermination;
        ntrees.InteractionTrees = ntrees.InteractionTrees + I.NumTrees;
      endif

      Mdl.ModelParameters = MP;
      Mdl.ReasonForTermination = reason;
      Mdl.NumTrainedTrees = ntrees;

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

    ## Drive the boosted-tree engine: the predictor phase, then a search for
    ## interactions worth adding, then the interaction phase over whichever
    ## pairs survived.  The two phases share a running fit, so the second
    ## continues from the prediction the first left rather than starting over.
    function this = fitBoosted (this, X, Y, Interactions, NTP, NTI, MSP, ...
                                MSI, LRP, LRI, MaxPValue, Verb, NPrint)

      ## Method 2 boosts the squared error, which is what a regression fits,
      ## weighted by W over the rows the fit sees.
      Wfit = this.W;
      [X, ~, cat, Levels] = gamCatCode (this.TreeModel, X);
      M = gamboosttrain (X, Y, 2, NTP, LRP, MSP, Verb, NPrint, [], Wfit, ...
                         cat);
      f = gamboostpredict (M.BinEdges, M.ShapeValues, X, M.Intercept);

      this.BinEdges  = M.BinEdges(:);   ## a column cell, as MATLAB reports it
      this.BinEdges(cat) = {[]};
      this.Intercept = M.Intercept;
      reason = struct ('PredictorTrees', M.ReasonForTermination, ...
                       'InteractionTrees', '');
      ntrees = struct ('PredictorTrees', M.NumTrees, ...
                       'InteractionTrees', 0);
      pairs = zeros (0, 2);
      pairValues = {};
      pairEdges = {};
      pairMissing = {};
      pairShift = 0;

      wanted = -1;
      if (ischar (Interactions))
        wanted = Inf;
      elseif (isscalar (Interactions) && ! isempty (Interactions))
        wanted = Interactions;
      elseif (! isempty (Interactions))
        pairs = interactionPairs (logical (Interactions));
      endif

      if (wanted > 0 && columns (X) > 1)
        S = gamboostpairs (X, M.Residuals, cat);
        ## The F ratio becomes a probability through the package's own fcdf,
        ## which is verified against MATLAB; the engine deliberately does not
        ## carry a second incomplete beta of its own.
        pval = 1 - fcdf (S.F, S.DF1, S.DF2);
        pval(S.DF1 <= 0) = 1;
        [pval, ord] = sort (pval);
        ranked = S.Pairs(ord, :);
        ranked = ranked(pval <= MaxPValue, :);
        if (isfinite (wanted) && rows (ranked) > wanted)
          ranked = ranked(1:wanted, :);
        endif
        pairs = ranked;
        ## A pair tree allowed a single split can only cut along one
        ## predictor, which is a main effect and no interaction, so no pair is
        ## kept.  R2024a selects none and warns as below.
        if (MSI < 2)
          pairs = zeros (0, 2);
        endif
        this.PairDetectionBinEdges = S.BinEdges(:);
        this.PairDetectionBinEdges(cat) = {[]};
        if (isempty (pairs))
          warning (strcat ("RegressionGAM: model does not include", ...
                           " interaction terms because all interaction", ...
                           " terms have p-values greater than the", ...
                           " 'MaxPValue' value, or the software was unable", ...
                           " to improve the model fit."));
        endif
      endif

      if (! isempty (pairs))
        I = gamboostinter (X, Y, f, 2, pairs, NTI, LRI, MSI, Wfit, cat);
        ## A pair tree splitting on one predictor alone is a main effect and
        ## adds nothing, so a fit made only of such trees keeps no pair, as
        ## R2024a keeps none.
        if (I.NumTrees == 0)
          pairs = zeros (0, 2);
          warning (strcat ("RegressionGAM: model does not include", ...
                           " interaction terms because all interaction", ...
                           " terms have p-values greater than the", ...
                           " 'MaxPValue' value, or the software was unable", ...
                           " to improve the model fit."));
        endif
      endif
      if (! isempty (pairs))
        this.Intercept = this.Intercept + I.Intercept;
        pairShift = I.Intercept;
        ## A requested pair list skips detection, so the grid it would have
        ## used comes from the interaction phase instead.
        if (wanted <= 0)
          this.PairDetectionBinEdges = I.PairBinEdges(:);
          this.PairDetectionBinEdges(cat) = {[]};
        endif
        pairValues = I.PairValues;
        pairEdges = I.PairEdges;
        pairMissing = I.PairMissing;
        reason.InteractionTrees = I.ReasonForTermination;
        ntrees.InteractionTrees = I.NumTrees;
      endif

      this.Interactions = pairs;
      this.ReasonForTermination = reason;
      this.NumTrainedTrees = ntrees;
      ## The constant the interaction surfaces gave up when they were
      ## recentred is kept apart from the predictor phase's intercept.  The
      ## Intercept property still reports their sum, as MATLAB's does, but
      ## predicting without the interactions has to take this part back out
      ## or it would answer with a constant the main effects never earned.
      this.TreeModel = struct ('ShapeValues', {M.ShapeValues}, ...
                               'PairValues', {pairValues}, ...
                               'PairEdges', {pairEdges}, ...
                               'PairMissing', {pairMissing}, ...
                               'Pairs', pairs, ...
                               'PairIntercept', pairShift, ...
                               'CategoricalLevels', {Levels});

      if (ischar (Interactions))
        request = Interactions;
      elseif (isempty (Interactions))
        request = 0;
      else
        request = Interactions;
      endif
      this.ModelParameters = struct ( ...
        'NumPrint', NPrint, ...
        'MaxPValue', MaxPValue, ...
        'InitialLearnRateForPredictors', LRP, ...
        'InitialLearnRateForInteractions', LRI, ...
        'NumTreesPerPredictor', NTP, ...
        'NumTreesPerInteraction', NTI, ...
        'MaxNumSplitsPerPredictor', MSP, ...
        'MaxNumSplitsPerInteraction', MSI, ...
        'VerbosityLevel', Verb, ...
        'Interactions', request, ...
        'Version', 1, ...
        'Method', 'GAM', ...
        'Type', 'regression');

    endfunction

    ## Determine interactions from Interactions optional parameter
    ## Fit the model that carries the interaction terms.  The constructor and
    ## addInteractions both arrive here with IntMatrix already decided and the
    ## predictors and response prepared as the fit wants them, so the two
    ## cannot drift: a model given its interactions after the fact is the
    ## model it would have been had they been asked for at the outset.
    function this = fitModelwInt (this, X, Y, Inter, Knots, Order, DoF)

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
    endfunction

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
      ## The fit is performed by the shared spline engine, which builds and
      ## factorises each predictor's design once and reduces a backfitting
      ## cycle to two products against the factors per term.
      Mdl = gamtrain (X, Y, Knots, Order, 2, Inter, this.Tol, 1000);
      iter  = Mdl.Iterations;
      param = Mdl.Parameters;
      res   = Mdl.Residuals;
      RSS   = Mdl.RSS;
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a RegressionGAM object
      mdl = RegressionGAM (1, 1);

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

## Reshape a predictor matrix into the terms a model was fitted on.  With
## APPEND true the interaction columns are added to the predictors, which is
## what an 'Interactions' model was fitted on; with it false they replace them,
## which is what a 'Formula' model was fitted on.  Every column of the result
## carries one term of the model, which is what predict_val walks.
function XA = gamTerms (X, IntMatrix, append)

  if (append)
    XA = X;
  else
    XA = [];
  endif
  for i = 1:rows (IntMatrix)
    tindex = logical (IntMatrix(i,:));
    Xterms = X(:,tindex);
    Xinter = ones (rows (X), 1);
    for c = 1:sum (tindex)
      Xinter = Xinter .* Xterms(:,c);
    endfor
    XA = [XA, Xinter];
  endfor

endfunction

## Helper function for making prediction of new data based on GAM model
function ypred = predict_val (params, X, intercept)
  ## The shared prediction engine evaluates every additive term and adds the
  ## intercept.
  ypred = gampredict (params, X, intercept);
endfunction

%!demo
%! ## Train a RegressionGAM Model for synthetic values
%! rng (42);
%! f1 = @(x) cos (3 * x);
%! f2 = @(x) x .^ 3;
%! x1 = 2 * rand (50, 1) - 1;
%! x2 = 2 * rand (50, 1) - 1;
%! y = f1(x1) + f2(x2);
%! y = y + y .* 0.2 .* rand (50,1);
%! X = [x1, x2];
%! a = fitrgam (X, y)

%!demo
%! ## Declare two different functions
%! rng (42);
%! f1 = @(x) cos (3 * x);
%! f2 = @(x) x .^ 3;
%!
%! ## Generate 80 samples for f1 and f2
%! x = [-4*pi:0.1*pi:4*pi-0.1*pi]';
%! X1 = f1(x);
%! X2 = f2(x);
%!
%! ## Create a synthetic response by adding noise
%! Ytrue = X1 + X2;
%! Y = Ytrue + Ytrue .* 0.2 .* rand (80,1);
%!
%! ## Assemble predictor data
%! X = [X1, X2];
%!
%! ## Train the GAM and test on the same data
%! ## A standard deviation and a prediction interval come from the spline
%! ## engine, which fits one; the boosted-tree engine reports none.
%! a = fitrgam (X, Y, 'FitMethod', 'splines', 'order', [5, 5]);
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
%! a = RegressionGAM (x, y, 'FitMethod', 'splines');
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
%! a = RegressionGAM (x, y, 'FitMethod', 'splines', ...
%!                    'predictors', pnames, 'formula', formula);
%! assert_equal (a.IntMatrix, double (intMat))
%! assert_equal ({a.ResponseName, a.PredictorNames}, {'Y', pnames})
%! assert_equal (a.Formula, formula)

%!test
%! ## Test that predict() executes correctly when interactions are present
%! X = [1, 2; 3, 4; 5, 6; 7, 8];
%! Y = [10; 20; 30; 40];
%! mdl = RegressionGAM (X, Y, 'FitMethod', 'splines', ...
%!                      'formula', 'Y ~ x1 + x2 + x1:x2');
%! ypred = predict (mdl, X);
%! assert_equal (isnumeric (ypred), true);
%! assert_equal (size (ypred), [4, 1]);
%! [ypred2, ySD, yInt] = predict (mdl, X, 'includeinteractions', true);
%! assert_equal (size (ypred2), [4, 1]);
%! assert_equal (size (ySD), [4, 1]);
%! assert_equal (size (yInt), [4, 2]);
%! ## The three terms fit these four points exactly, so the residuals are zero
%! ## and so is ySD.  This block asserted the three sizes and nothing else,
%! ## which is why it stayed green while ySD was computed from two of the
%! ## three terms: a size is true whatever the number is.
%! assert_equal (ypred2, [10; 20; 30; 40], 1e-10);
%! assert_equal (ySD, zeros (4, 1), 1e-10);

%!test
%! ## Verify ySD is based on training residual variance
%! X = (1:10)';
%! Y = [2; 1; 4; 3; 6; 5; 8; 7; 10; 9];
%! mdl = RegressionGAM (X, Y, 'FitMethod', 'splines');
%! y_train = predict (mdl, X);
%! rs = Y - y_train;
%! expected_ySD = sqrt (var (rs));
%! [~, ySD] = predict (mdl, X(1:4,:));
%! assert_equal (ySD, expected_ySD * ones (4, 1), 1e-10);

%!test
%! ## Verify ySD remains the same for one or more prediction points
%! X = (1:10)';
%! Y = [2; 1; 4; 3; 6; 5; 8; 7; 10; 9];
%! mdl = RegressionGAM (X, Y, 'FitMethod', 'splines');
%! y_train = predict (mdl, X);
%! expected_ySD = sqrt (var (Y - y_train));
%! [~, ySD_1] = predict (mdl, X(1,:));
%! [~, ySD_3] = predict (mdl, X(1:3,:));
%! assert_equal (ySD_1, expected_ySD, 1e-10);
%! assert_equal (ySD_3, expected_ySD * ones (3, 1), 1e-10);

## Test input validation for constructor
## Interactions reports the two-way terms the fitted model carries, as
## predictor index pairs, matching what R2024a's fitrgam returns.
%!test
%! k = (1:60)';
%! X = [mod(k*7,11)-5, mod(k*3,11)-5, mod(k*5,11)-5];
%! y = X(:,1) .* X(:,2) + 0.5 * X(:,3);
%! Mdl = fitrgam (X, y, "Interactions", "all");
%! assert_equal (Mdl.Interactions, [1, 2; 1, 3; 2, 3]);

## No interactions is an empty list of pairs, keeping its two columns.
%!test
%! k = (1:60)';
%! X = [mod(k*7,11)-5, mod(k*3,11)-5, mod(k*5,11)-5];
%! Mdl = fitrgam (X, X(:,1) + X(:,3));
%! assert_equal (size (Mdl.Interactions), [0, 2]);

## A formula's main effects are terms of the model but not interactions.
%!test
%! k = (1:60)';
%! X = [mod(k*7,11)-5, mod(k*3,11)-5, mod(k*5,11)-5];
%! y = X(:,1) .* X(:,2) + 0.5 * X(:,3);
%! Mdl = fitrgam (X, y, "FitMethod", "splines", ...
%!                "Formula", "Y ~ x1 + x2 + x1:x2");
%! assert_equal (Mdl.Interactions, [1, 2]);
%! assert_equal (compact (Mdl).Interactions, [1, 2]);

## addInteractions fits the interaction terms onto a model that already has
## its univariate ones.  The result is the model that would have been fitted
## had the terms been asked for at the outset: both go through one private
## method, so the two cannot drift apart.
%!test
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! Xai = meas(bai,2:4); Yai = meas(bai,1);
%! Aai = addInteractions (fitrgam (Xai, Yai), "all");
%! Bai = fitrgam (Xai, Yai, "Interactions", "all");
%! assert_equal (Aai.Interactions, Bai.Interactions);
%! assert_equal (Aai.ModelwInt, Bai.ModelwInt);
%! assert_equal (predict (Aai, Xai), predict (Bai, Xai));

## The univariate fit is left alone, which is what MATLAB leaves alone too.
%!test
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! Xai = meas(bai,2:4); Yai = meas(bai,1);
%! Cai = fitrgam (Xai, Yai);
%! Aai = addInteractions (Cai, "all");
%! assert_equal (predict (Aai, Xai, "IncludeInteractions", false), ...
%!               predict (Cai, Xai));

## A count and a logical matrix name terms as the constructor's option does.
%!test
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! Xai = meas(bai,2:4); Yai = meas(bai,1);
%! Aai = addInteractions (fitrgam (Xai, Yai), 2);
%! assert_equal (size (Aai.Interactions), [2, 2]);
%! assert_equal (sort (Aai.Interactions, 2), Aai.Interactions);
%! Lai = addInteractions (fitrgam (Xai, Yai), logical ([1 1 0; 0 1 1]));
%! assert_equal (Lai.Interactions, [1, 2; 2, 3]);

## A model that already carries interaction terms is not extended, and a
## model fitted from a formula names every term it has, interactions among
## them, so it is refused for the same reason.  R2024a refuses both.
## resume continues the phase that ran last, and a model with no interactions
## has only the predictor phase open.  Resuming reproduces the model a single
## fit of the combined budget would have produced, which is the oracle every
## test here uses.
%!test
%! load fisheriris
%! X = meas;
%! Y = meas(:,1) + 0.3 * meas(:,3);
%! A = fitrgam (X, Y, 'NumTreesPerPredictor', 5);
%! B = resume (A, 10);
%! C = fitrgam (X, Y, 'NumTreesPerPredictor', 15);
%! assert_equal (B.ModelParameters.NumTreesPerPredictor, 15);
%! assert_equal (predict (B, X), predict (C, X), 1e-10);

%!test
%! ## A model carrying interactions gains interaction trees, and its predictor
%! ## shape functions are left where they were.
%! load fisheriris
%! X = meas;
%! Y = meas(:,1) + 0.3 * meas(:,3);
%! A = fitrgam (X, Y, 'NumTreesPerPredictor', 5, 'Interactions', 3, ...
%!              'NumTreesPerInteraction', 4);
%! B = resume (A, 10);
%! assert_equal (B.ModelParameters.NumTreesPerPredictor, 5);
%! assert_equal (B.ModelParameters.NumTreesPerInteraction, 14);
%! assert_equal (B.TreeModel.ShapeValues, A.TreeModel.ShapeValues);
%! C = fitrgam (X, Y, 'NumTreesPerPredictor', 5, 'Interactions', 3, ...
%!              'NumTreesPerInteraction', 14);
%! assert_equal (predict (B, X), predict (C, X), 1e-10);

%!test
%! ## The selected pairs survive, and resuming twice accumulates.
%! load fisheriris
%! X = meas;
%! Y = meas(:,1) + 0.3 * meas(:,3);
%! A = fitrgam (X, Y, 'NumTreesPerPredictor', 5, 'Interactions', 3, ...
%!              'NumTreesPerInteraction', 4);
%! B = resume (resume (A, 10), 6);
%! assert_equal (B.Interactions, A.Interactions);
%! assert_equal (B.ModelParameters.NumTreesPerInteraction, 20);
%! assert_equal (B.ModelParameters.NumTreesPerPredictor, 5);

%!test
%! ## The model handed in is not modified.
%! load fisheriris
%! X = meas;
%! Y = meas(:,1) + 0.3 * meas(:,3);
%! A = fitrgam (X, Y, 'NumTreesPerPredictor', 5);
%! B = resume (A, 10);
%! assert_equal (A.ModelParameters.NumTreesPerPredictor, 5);
%! assert_equal (B.ModelParameters.NumTreesPerPredictor, 15);

%!test  # a row missing a predictor is kept, and its spline term has no value
%! x = linspace (0, 1, 30)';
%! X = [x, cos(4 * x)];
%! y = sin (3 * x) + X(:,2);
%! Mdl = RegressionGAM (X, y, 'FitMethod', 'splines');
%! yFit = predict (Mdl, [0.5, 0.2; NaN, 0.2; 0.5, 0.2]);
%! assert_equal (size (yFit), [3, 1]);
%! assert_equal (isnan (yFit)', [false, true, false]);

%!test  # MATLAB parity: a missing value adds nothing from its term
%! k = (1:200)';
%! X = [sin(k), cos(2 * k), mod(k, 5)];
%! y = 2 * sin (k) + X(:,2) .^ 2 + 0.3 * X(:,3);
%! Q = [0.5, 0.2, 1; NaN, 0.2, 1; 0.5, NaN, 1; NaN, NaN, NaN; 0.1, 0.2, 1; ...
%!      NaN, 0.7, 3; NaN, NaN, 1];
%! Mdl = RegressionGAM (X, y);
%! yFit = predict (Mdl, Q);
%! assert_equal (yFit, [1.4399946345; 0.675842807654; 1.60393628463; ...
%!                      1.09800542249; 1.33125338416; 1.62383199199; ...
%!                      0.839784457781], 1e-10);
%! assert_equal (yFit(4), Mdl.Intercept);

%!test  # MATLAB parity: categorical predictors are split by sets of levels
%! k = (0:119)';
%! c1 = mod (k, 3) + 1;
%! x2 = sin (k);
%! c3 = 10 * (mod (floor (k / 2), 2) + 1);
%! X = [c1, x2, c3];
%! y = 5 * (c1 == 2) + 0.5 * x2 - 3 * (c3 == 20) + 0.1 * cos (k);
%! Q = [1, 0, 10; 2, 0, 10; 3, 0, 10; 1, 0, 20; 1, 0.5, 10; 4, 0, 10; ...
%!      NaN, 0, 10; 2.5, 0, 20];
%! Mdl = RegressionGAM (X, y, 'CategoricalPredictors', [1, 3], ...
%!                      'NumTreesPerPredictor', 5, ...
%!                      'MaxNumSplitsPerPredictor', 2);
%! assert_equal (Mdl.CategoricalPredictors, [1, 3]);
%! assert_equal (predict (Mdl, Q), [-0.1559786927; 4.811738827; ...
%!                                  -0.1849962119; -3.091005709; ...
%!                                  0.2642714677; 1.490254641; ...
%!                                  1.490254641; -1.444772376], 1e-9);

%!test  # MATLAB parity: a categorical predictor reports no bin edges
%! k = (0:119)';
%! c1 = mod (k, 3) + 1;
%! x2 = sin (k);
%! c3 = 10 * (mod (floor (k / 2), 2) + 1);
%! X = [c1, x2, c3];
%! y = 5 * (c1 == 2) + 0.5 * x2 - 3 * (c3 == 20) + 0.1 * cos (k);
%! Q = [1, 0, 10; 2, 0, 10; 3, 0, 10; 1, 0, 20; 1, 0.5, 10; 4, 0, 10; ...
%!      NaN, 0, 10; 2.5, 0, 20];
%! Mdl = RegressionGAM (X, y, 'CategoricalPredictors', [1, 3], ...
%!                      'NumTreesPerPredictor', 1);
%! assert_equal (Mdl.Intercept, 0.1666859425, 1e-10);
%! assert_equal (predict (Mdl, Q(6:8,:)), [1.312029301; 1.312029301; ...
%!                                         -1.640617244], 1e-9);
%! assert_equal (Mdl.BinEdges{1}, []);
%! assert_equal (Mdl.BinEdges{3}, []);
%! assert_equal (numel (Mdl.BinEdges{2}), 119);

%!test  # MATLAB parity: pair trees split categorical predictors by level sets
%! k = (0:149)';
%! c1 = 1 + mod (floor (k * 7 / 11), 4);
%! x2 = sin (k);
%! c3 = 10 * (1 + (mod (k, 7) > 3)) + 10 * (mod (k, 13) == 0);
%! X = [c1, x2, c3];
%! y = 2 * (c1 == 2) - 1.5 * (c1 == 4) + 0.5 * x2 - 3 * (c3 == 20) ...
%!     + (c3 == 30) + (c1 == 3) .* (c3 == 20) + 0.1 * cos (k) ...
%!     + 0.05 * sin (3 * k);
%! [A, B] = ndgrid ([1, 2, 3, 4], [10, 20, 30]);
%! Q = [A(:), zeros(12, 1), B(:)];
%! Mdl = RegressionGAM (X, y, 'CategoricalPredictors', [1, 3], ...
%!                      'Interactions', logical ([1, 0, 1]), ...
%!                      'NumTreesPerInteraction', 5);
%! d = predict (Mdl, Q) - predict (Mdl, Q, 'IncludeInteractions', false);
%! assert_equal (d, [0.03701741101; 0.06202240762; -0.1556975657; ...
%!                   0.08736164669; -0.07100596369; -0.08113441688; ...
%!                   0.3496023158; -0.06710888279; 0.06027597832; ...
%!                   0.05025393807; -0.06227465812; 0.01256864956], 1e-9);

%!test  # MATLAB parity: a pair tree using one predictor alone adds nothing
%! k = (0:149)';
%! c1 = 1 + mod (floor (k * 7 / 11), 4);
%! x2 = sin (k);
%! c3 = 10 * (1 + (mod (k, 7) > 3)) + 10 * (mod (k, 13) == 0);
%! X = [c1, x2, c3];
%! y = 2 * (c1 == 2) - 1.5 * (c1 == 4) + 0.5 * x2 - 3 * (c3 == 20) ...
%!     + (c3 == 30) + (c1 == 3) .* (c3 == 20) + 0.1 * cos (k) ...
%!     + 0.05 * sin (3 * k);
%! g = linspace (-1, 1, 5)';
%! [A, B, C] = ndgrid ([1, 2, 3, 4, NaN], g, [10, 20, 30, NaN]);
%! Q = [A(:), B(:), C(:)];
%! Q = Q([1, 7, 13, 19, 28, 34, 45, 52, 59, 67],:);
%! Mdl = RegressionGAM (X, y, 'CategoricalPredictors', [1, 3], ...
%!                      'Interactions', 'all', 'NumTreesPerInteraction', 5);
%! assert_equal (Mdl.Interactions, [1, 3; 2, 3; 1, 2]);
%! d = predict (Mdl, Q) - predict (Mdl, Q, 'IncludeInteractions', false);
%! ## Five rounds, not twenty: a fit long enough to compound its rounding
%! ## parts company with a toolchain that contracts multiply-add pairs, and
%! ## then no tolerance holds, one element changing sign.  Measured on R2024a,
%! ## 2026-09-16, and agreeing to thirteen digits or better.
%! assert_equal (d, [0.14496260366867; 0.00289829354386018; ...
%!                   -0.125368787144964; 0.0691339453924824; ...
%!                   0.5053322741541; -0.0581889291447562; ...
%!                   0.00447142646759335; 0.15228433993982; ...
%!                   0.0385020050553703; 0.0838177432700107], 1e-9);

%!test  # MATLAB parity: a pair whose trees split one predictor alone is dropped
%! k = (0:119)';
%! c1 = mod (k, 3) + 1;
%! x2 = sin (k);
%! c3 = 10 * (mod (floor (k / 2), 2) + 1);
%! X = [c1, x2, c3];
%! y = 5 * (c1 == 2) + 0.5 * x2 - 3 * (c3 == 20) + 0.1 * cos (k);
%! warning ('off', 'all', 'local');
%! Mdl = RegressionGAM (X, y, 'CategoricalPredictors', [1, 3], ...
%!                      'Interactions', logical ([0, 1, 1]), ...
%!                      'NumTreesPerInteraction', 5);
%! assert_equal (size (Mdl.Interactions), [0, 2]);
%!warning<RegressionGAM: model does not include interaction terms because all interaction terms have p-values greater than the 'MaxPValue' value, or the software was unable to improve the model fit.> ...
%! k = (0:119)';
%! RegressionGAM ([mod(k, 3) + 1, sin(k), 10 * (mod (floor (k / 2), 2) + 1)], ...
%!                5 * (mod (k, 3) == 1) + 0.5 * sin (k) ...
%!                - 3 * (mod (floor (k / 2), 2) == 1) + 0.1 * cos (k), ...
%!                'CategoricalPredictors', [1, 3], ...
%!                'Interactions', logical ([0, 1, 1]), ...
%!                'NumTreesPerInteraction', 5);

%!test  # MATLAB parity: 'CategoricalPredictors', 'all'
%! k = (0:59)';
%! X = [mod(k, 3) + 1, mod(k, 4)];
%! y = X(:,1) .^ 2 - X(:,2) + 0.1 * cos (k);
%! Mdl = RegressionGAM (X, y, 'CategoricalPredictors', 'all', ...
%!                      'Interactions', 'all', 'NumTreesPerInteraction', 2);
%! assert_equal (Mdl.CategoricalPredictors, [1, 2]);
%! assert_equal (Mdl.BinEdges, {[]; []});
%! assert_equal (Mdl.PairDetectionBinEdges, {[]; []});

%!test  # MATLAB parity: observation weights enter the fit
%! k = (1:200)';
%! X = [sin(k), cos(2 * k), mod(k, 5)];
%! y = 2 * sin (k) + X(:,2) .^ 2 + 0.3 * X(:,3);
%! w = 1 + mod (k, 3);
%! Q = [0.5, 0.2, 1; 0.1, 0.7, 3; -0.4, -0.2, 0];
%! Mdl = RegressionGAM (X, y, 'Weights', w);
%! assert_equal (Mdl.W, w / sum (w), 1e-15);
%! assert_equal (Mdl.Intercept, 1.09607817728, 1e-10);
%! assert_equal (predict (Mdl, Q), [1.32799066689; 2.22459947488; ...
%!                                  -0.705872562832], 1e-10);
%! M5 = RegressionGAM (X, y, 'Weights', 5 * w);
%! assert_equal (predict (M5, Q), predict (Mdl, Q), 1e-10);

%!test  # MATLAB parity: a row missing a predictor is fitted, not dropped
%! k = (1:200)';
%! X = [sin(k), cos(2 * k), mod(k, 5)];
%! y = 2 * sin (k) + X(:,2) .^ 2 + 0.3 * X(:,3);
%! X(1:5,1) = NaN;
%! Q = [0.5, 0.2, 1; NaN, 0.2, 1; 0.5, NaN, 1; NaN, NaN, NaN; 0.1, 0.2, 1; ...
%!      NaN, 0.7, 3; NaN, NaN, 1];
%! Mdl = RegressionGAM (X, y);
%! assert_equal (Mdl.NumObservations, 200);
%! assert_equal (Mdl.Intercept, 1.28924934314, 1e-10);
%! assert_equal (predict (Mdl, Q), [1.39295651846; 0.938063487647; ...
%!                                  1.51847439653; 1.28924934314; ...
%!                                  1.06837329164; 1.96571223429; ...
%!                                  1.06358136571], 1e-10);

%!test  # MATLAB parity: a pair tree of one split selects no interaction
%! k = (1:200)';
%! X = [sin(k), cos(2 * k), mod(k, 5)];
%! y = 2 * sin (k) + X(:,2) .^ 2 + 0.3 * X(:,3);
%! S = warning ('off', 'all');
%! unwind_protect
%!   Mdl = RegressionGAM (X, y, 'Interactions', 1, ...
%!                        'MaxNumSplitsPerInteraction', 1);
%! unwind_protect_cleanup
%!   warning (S);
%! end_unwind_protect
%! assert_equal (Mdl.Interactions, zeros (0, 2));
%! assert_equal (Mdl.Intercept, 1.09800542249, 1e-10);

%!test  # MATLAB parity: pair trees are fitted to the rows, on their own grid
%! k = (1:200)';
%! X = [sin(k), cos(2 * k), mod(k, 5)];
%! y = 2 * sin (k) + X(:,2) .^ 2 + 0.3 * X(:,3);
%! P = [0, -0.99, 0; 0, -0.42, 2; 0, 0.33, 1; 0, 0.78, 4];
%! Mdl = RegressionGAM (X, y, 'Interactions', 1, ...
%!                      'NumTreesPerInteraction', 1, ...
%!                      'MaxNumSplitsPerInteraction', 4);
%! assert_equal (Mdl.Interactions, [2, 3]);
%! pc = predict (Mdl, P) - predict (Mdl, P, 'IncludeInteractions', false);
%! assert_equal (pc, [-0.014545869; 0.0020626874; 0.013767507; ...
%!                    -0.0043992362], 1e-8);

%!test  # resuming the interaction phase adds trees a longer run would add
%! k = (1:200)';
%! X = [sin(k), cos(2 * k), mod(k, 5)];
%! y = 2 * sin (k) + X(:,2) .^ 2 + 0.3 * X(:,3);
%! P = [0, -0.99, 0; 0, -0.42, 2; 0, 0.33, 1; 0, 0.78, 4];
%! M1 = RegressionGAM (X, y, 'Interactions', 1, ...
%!                     'NumTreesPerInteraction', 1, ...
%!                     'MaxNumSplitsPerInteraction', 4);
%! M2 = RegressionGAM (X, y, 'Interactions', 1, ...
%!                     'NumTreesPerInteraction', 2, ...
%!                     'MaxNumSplitsPerInteraction', 4);
%! assert_equal (predict (resume (M1, 1), P), predict (M2, P), 1e-12);

%!test  # MATLAB parity: a missing value stops at the node that splits on it
%! k = (1:200)';
%! X = [sin(k), cos(2 * k), mod(k, 5)];
%! y = 2 * sin (k) + X(:,2) .^ 2 + 0.3 * X(:,3);
%! X(1:5,1) = NaN;
%! X(6:10,3) = NaN;
%! P = [-0.99, 0, 0; -0.33, 0, 2; 0.33, 0, 4; 0.99, 0, 1; ...
%!      NaN, 0, 0; -0.9, 0, NaN; 0.9, 0, NaN; NaN, 0, NaN];
%! Mdl = RegressionGAM (X, y, 'Interactions', logical ([1, 0, 1]), ...
%!                      'NumTreesPerInteraction', 1, ...
%!                      'MaxNumSplitsPerInteraction', 4);
%! assert_equal (Mdl.Interactions, [1, 3]);
%! assert_equal (Mdl.Intercept, 1.27402212251, 1e-10);
%! pc = predict (Mdl, P) - predict (Mdl, P, 'IncludeInteractions', false);
%! assert_equal (pc, [0.004044014047; -0.01436832448; 0.01780996877; ...
%!                    -0.01436832448; -0.0001224643025; -0.002123518423; ...
%!                    -0.002123518423; -0.0001224643025], 1e-9);

%!test  # resuming keeps what a row missing a pair predictor takes
%! k = (1:200)';
%! X = [sin(k), cos(2 * k), mod(k, 5)];
%! y = 2 * sin (k) + X(:,2) .^ 2 + 0.3 * X(:,3);
%! X(1:5,1) = NaN;
%! X(6:10,3) = NaN;
%! P = [-0.99, 0, 0; 0.33, 0, 4; NaN, 0, 0; -0.9, 0, NaN; NaN, 0, NaN];
%! M1 = RegressionGAM (X, y, 'Interactions', logical ([1, 0, 1]), ...
%!                     'NumTreesPerInteraction', 1, ...
%!                     'MaxNumSplitsPerInteraction', 4);
%! M2 = RegressionGAM (X, y, 'Interactions', logical ([1, 0, 1]), ...
%!                     'NumTreesPerInteraction', 2, ...
%!                     'MaxNumSplitsPerInteraction', 4);
%! assert_equal (predict (resume (M1, 1), P), predict (M2, P), 1e-12);

%!test  # MATLAB parity: the detection grid, with values missing
%! k = (1:300)';
%! X = [sin(k), cos(2 * k), mod(k, 5), sin(0.7 * k), cos(0.3 * k)];
%! y = X(:,1) + X(:,2) .^ 2 + 0.3 * X(:,3) + 0.8 * X(:,1) .* X(:,4) ...
%!     + 0.5 * X(:,2) .* X(:,3) + 0.3 * X(:,4) .* X(:,5) + 0.05 * sin (13 * k);
%! X(1:15,1) = NaN;
%! X(16:30,3) = NaN;
%! X(31:45,4) = NaN;
%! Mdl = RegressionGAM (X, y, 'Interactions', 'all', ...
%!                      'NumTreesPerInteraction', 1);
%! assert_equal (Mdl.PairDetectionBinEdges{4}(:)', ...
%!               [-0.92412955668638763, -0.67185165069256547, ...
%!                -0.37673775931105463, 0.064103966449078689, ...
%!                0.43097453606213998, 0.73708146543750352, ...
%!                0.93347437171529057], 1e-15);

%!error<RegressionGAM: 'CategoricalPredictors' indices must not exceed the number of predictors.> ...
%! RegressionGAM (ones (10, 2), (1:10)', 'CategoricalPredictors', 3)
%!error<RegressionGAM: a logical 'CategoricalPredictors' must have one element per predictor.> ...
%! RegressionGAM (ones (10, 2), (1:10)', 'CategoricalPredictors', true)
%!error<RegressionGAM: 'CategoricalPredictors' must be a vector of positive integers, a logical vector, a character matrix, a string array, a cell array of character vectors or 'all'.> ...
%! RegressionGAM (ones (10, 2), (1:10)', 'CategoricalPredictors', 0.5)
%!error<RegressionGAM: 'categoricalpredictors' is a parameter of the boosted-tree engine and cannot be used with 'FitMethod' 'splines'.> ...
%! RegressionGAM (ones (10, 2), (1:10)', 'FitMethod', 'splines', ...
%!                'CategoricalPredictors', 1)
%!error<RegressionGAM: 'Weights' must be a numeric vector.> ...
%! RegressionGAM (ones (10, 2), (1:10)', 'Weights', 'a')
%!error<RegressionGAM: 'Weights' must have one element per row of X.> ...
%! RegressionGAM (ones (10, 2), (1:10)', 'Weights', [1, 2])
%!error<RegressionGAM: 'Weights' must hold finite non-negative values.> ...
%! RegressionGAM (ones (10, 2), (1:10)', 'Weights', -ones (10, 1))
%!error<RegressionGAM: 'weights' is a parameter of the boosted-tree engine and cannot be used with 'FitMethod' 'splines'.> ...
%! RegressionGAM ([(1:10)', mod((1:10)', 3)], (1:10)', ...
%!                'Weights', ones (10, 1), 'FitMethod', 'splines')
%!error<RegressionGAM.resume: Not enough input arguments.> ...
%! load fisheriris; ...
%! resume (fitrgam (meas, meas(:,1), 'NumTreesPerPredictor', 5))
%!error<RegressionGAM.resume: resuming is available only under 'FitMethod', 'boostedtrees'; a spline backfit stops at its tolerance and has no budget to extend.> ...
%! load fisheriris; ...
%! resume (fitrgam (meas, meas(:,1), 'FitMethod', 'splines'), 5)
%!error<RegressionGAM.resume: NUMTREES must be a positive integer scalar.> ...
%! load fisheriris; ...
%! resume (fitrgam (meas, meas(:,1), 'NumTreesPerPredictor', 5), 0)
%!error<RegressionGAM.resume: NUMTREES must be a positive integer scalar.> ...
%! load fisheriris; ...
%! resume (fitrgam (meas, meas(:,1), 'NumTreesPerPredictor', 5), 2.5)
%!error<RegressionGAM.resume: NUMTREES must be a positive integer scalar.> ...
%! load fisheriris; ...
%! resume (fitrgam (meas, meas(:,1), 'NumTreesPerPredictor', 5), [1, 2])

%!error<RegressionGAM.resume: unable to resume training because the software was unable to improve the model fit.> ...
%! X = [ones(20,1)*[1, 2]; ones(20,1)*[3, 4]]; ...
%! resume (fitrgam (X, ones (40, 1), 'NumTreesPerPredictor', 5), 5)

%!error<RegressionGAM.addInteractions: adding interaction terms to a model that already includes them is not supported.> ...
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! Mai = fitrgam (meas(bai,2:4), meas(bai,1), "Interactions", 2);
%! addInteractions (Mai, "all")
%!error<RegressionGAM.addInteractions: adding interaction terms to a model that already includes them is not supported.> ...
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! addInteractions (fitrgam (meas(bai,2:4), meas(bai,1), ...
%!                  "FitMethod", "splines", ...
%!                  "Formula", "Y ~ x1 + x2 + x1:x2"), "all")
%!error<RegressionGAM.addInteractions: invalid 'Interactions' parameter.> ...
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! addInteractions (fitrgam (meas(bai,2:4), meas(bai,1)), {1})

%!error<RegressionGAM: too few input arguments.> RegressionGAM ()
%!error<RegressionGAM: too few input arguments.> RegressionGAM (ones (10,2))
%!error<RegressionGAM: number of rows in X and Y must be equal.> ...
%! RegressionGAM (ones (10,2), ones (5,1))
%!error<RegressionGAM: invalid values in X.> ...
%! RegressionGAM ([1;2;3;'a';4], ones (5,1))
%!error<RegressionGAM: invalid optional paired argument.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'some', 'some')
%!error<RegressionGAM: Formula must be a string.>
%! RegressionGAM (ones (10,2), ones (10,1), 'formula', {'y~x1+x2'})
%!error<RegressionGAM: Formula must be a string.>
%! RegressionGAM (ones (10,2), ones (10,1), 'formula', [0, 1, 0])
%!error<RegressionGAM: invalid syntax in Formula.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'FitMethod', 'splines', ...
%!                'formula', 'something')
%!error<RegressionGAM: no predictor terms in Formula.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'FitMethod', 'splines', ...
%!                'formula', 'something~')
%!error<RegressionGAM: no predictor terms in Formula.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'FitMethod', 'splines', ...
%!                'formula', 'something~')
%!error<RegressionGAM: some predictors have not been identified> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'FitMethod', 'splines', ...
%!                'formula', 'something~x1:')
%!error<RegressionGAM: invalid Interactions parameter.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'interactions', 'some')
%!error<RegressionGAM: invalid Interactions parameter.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'interactions', -1)
%!error<RegressionGAM: invalid Interactions parameter.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'interactions', [1 2 3 4])
%!error<RegressionGAM: number of interaction terms requested is larger than> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'FitMethod', 'splines', ...
%!                'interactions', 3)
%!error<RegressionGAM: 'Formula' and 'Interactions' cannot be given together.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'FitMethod', 'splines', ...
%!                'formula', 'y ~ x1 + x2', 'interactions', 1)
%!error<RegressionGAM: 'Formula' and 'Interactions' cannot be given together.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'FitMethod', 'splines', ...
%!                'interactions', 1, 'formula', 'y ~ x1 + x2')
%!error<RegressionGAM: invalid value for Knots.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'knots', 'a')
%!error<RegressionGAM: at most two of 'Knots', 'Order' and 'DoF' may be given.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'FitMethod', 'splines', ...
%!                'order', 3, 'dof', 2, 'knots', 5)
%!error<RegressionGAM: invalid value for DoF.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'dof', 'a')
%!error<RegressionGAM: at most two of 'Knots', 'Order' and 'DoF' may be given.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'FitMethod', 'splines', ...
%!                'knots', 5, 'order', 3, 'dof', 2)
%!error<RegressionGAM: invalid value for Order.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'order', 'a')
%!error<RegressionGAM: at most two of 'Knots', 'Order' and 'DoF' may be given.> ...
%! RegressionGAM (ones (10,2), ones (10,1), 'FitMethod', 'splines', ...
%!                'knots', 5, 'dof', 2, 'order', 2)
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
%!error<RegressionGAM.predict: invalid optional paired argument.> ...
%! predict (RegressionGAM (ones (10,2), ones (10,1)), ones (10,2), 'some', 'some')
%!error<RegressionGAM.predict: optional arguments must be given in Name-Value pairs.> ...
%! predict (RegressionGAM (ones (10,2), ones (10,1)), ones (10,2), 'Alpha')
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
%! Mdl = fitrgam (meas(:,1:3), meas(:,4), 'FitMethod', 'splines');
%! assert_equal (Mdl.Intercept, Mdl.BaseModel.Intercept);
%! assert_equal (size (Mdl.W), [Mdl.NumObservations, 1]);
%! assert_equal (sum (Mdl.W), 1, 1e-12);
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (Mdl.ExpandedPredictorNames, Mdl.PredictorNames);
%! assert_equal (Mdl.IsStandardDeviationFit, false);

## A scalar Knots, Order or DoF applies to every predictor.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,1:3), meas(:,4), 'FitMethod', 'splines', ...
%!                'Knots', 4);
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
%! Mdl = fitrgam (meas(:,1:3), meas(:,4), 'FitMethod', 'splines', ...
%!                'Interactions', 'all');
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
%!error<RegressionGAM.loss: invalid optional paired argument.> ...
%! loss (Mr, xr, yr, 'Bogus', 1)
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

## A fitted model survives savemodel and loadmodel: the properties come
## back as they were and it predicts the same.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! Y = meas(:,1);
%! Mdl = fitrgam (X, Y, 'FitMethod', 'splines');
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'RegressionGAM');
%! assert_equal (M2.NumObservations, Mdl.NumObservations);
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
%! Mdl = fitrgam (X, Y, 'FitMethod', 'boostedtrees');
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (M2.FitMethod, 'boostedtrees');
%! assert_equal (M2.TreeModel.ShapeValues, Mdl.TreeModel.ShapeValues);
%! assert_equal (M2.BinEdges, Mdl.BinEdges);
%! assert_equal (predict (M2, X(1:5,:)), predict (Mdl, X(1:5,:)), 1e-12);

## crossval refits one compact model per fold over the observations used.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,1:3), meas(:,4));
%! CVMdl = crossval (Mdl, 'KFold', 3);
%! assert_equal (class (CVMdl), 'RegressionPartitionedModel');
%! assert_equal (CVMdl.CrossValidatedModel, 'GAM');
%! assert_equal (class (CVMdl.Trained{1}), 'CompactRegressionGAM');
%! assert_equal (CVMdl.KFold, 3);
%! assert_equal (numel (CVMdl.Trained), 3);
%! assert_equal (CVMdl.NumObservations, 150);

%!test
%! load fisheriris
%! CVMdl = crossval (fitrgam (meas(1:20,1:3), meas(1:20,4)));
%! assert_equal (CVMdl.KFold, 10);

%!test
%! load fisheriris
%! CVMdl = crossval (fitrgam (meas(1:20,1:3), meas(1:20,4)), 'Holdout', 0.25);
%! assert_equal (CVMdl.KFold, 1);

%!test
%! load fisheriris
%! CVMdl = crossval (fitrgam (meas(1:12,1:3), meas(1:12,4)), 'Leaveout', 'on');
%! assert_equal (CVMdl.KFold, 12);

%!test
%! load fisheriris
%! cvp = cvpartition (20, 'KFold', 4);
%! Mdl = fitrgam (meas(1:20,1:3), meas(1:20,4));
%! assert_equal (crossval (Mdl, 'CVPartition', cvp).KFold, 4);

## The fold models are refitted with the spline parameterisation of the model
## they came from.  The compact fold does not report the parameters themselves.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(1:20,1:3), meas(1:20,4), 'FitMethod', 'splines', ...
%!                'Knots', 6, 'Order', 3);
%! CVMdl = crossval (Mdl, 'KFold', 3);
%! assert_equal (CVMdl.Trained{1}.FitMethod, 'splines');
%! assert_equal (isfield (CVMdl.Trained{1}.BaseModel, 'Intercept'), true);

## Held-out error exceeds the resubstitution error of the same fit.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,1:3), meas(:,4));
%! assert (kfoldLoss (crossval (Mdl, 'KFold', 5)) > resubLoss (Mdl));

%!shared cvobj
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1; 4, 5, 6];
%! y = [1; 2; 3; 4; 5];
%! cvobj = fitrgam (x, y);
%!error<RegressionGAM.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (cvobj, 'kfold')
%!error<RegressionGAM.crossval: specify only one of the optional Name-Value paired arguments.> ...
%! crossval (cvobj, 'kfold', 3, 'holdout', 0.2)
%!error<RegressionGAM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (cvobj, 'kfold', 'a')
%!error<RegressionGAM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (cvobj, 'holdout', 2)
%!error<RegressionGAM.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (cvobj, 'leaveout', 1)
%!error<RegressionGAM.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (cvobj, 'cvpartition', 1)
%!error<RegressionGAM.crossval: invalid parameter name in optional paired arguments.> ...
%! crossval (cvobj, 'bogus', 1)

## BinEdges is an empty cell, as MATLAB reports for every learner that
## does no binning.  MATLAB's own generalized additive model fills it,
## being boosted trees where this one is splines.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,1:3), meas(:,4));
%! assert_equal (class (Mdl.BinEdges), 'cell');
%! assert_equal (numel (Mdl.BinEdges), 3);
%! Msp = fitrgam (meas(:,1:3), meas(:,4), 'FitMethod', 'splines');
%! assert_equal (Msp.BinEdges, {});

## ySD is computed from the model's own terms, not from the stored predictors.
## The two disagree whenever the model was reshaped, and both directions were
## wrong: an 'Interactions' model has more terms than X has columns and the
## residuals came from a truncated model, and a 'Formula' model can have fewer
## and the call raised.

%!test
%! load fisheriris
%! X = meas(:,1:3);
%! Y = meas(:,4);
%! Mdl = fitrgam (X, Y, 'FitMethod', 'splines', 'Interactions', 'all');
%! [~, ySD] = predict (Mdl, X(1:4,:));
%! ## Six terms against three stored columns.  Truncating to three gave
%! ## 0.9379137529, close to seven times the answer.
%! assert_equal (ySD, 0.136050646147 * ones (4, 1), 1e-9);

## The same number, derived rather than pinned: ySD is the residual standard
## deviation of the full model over the training data.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! Y = meas(:,4);
%! Mdl = fitrgam (X, Y, 'FitMethod', 'splines', 'Interactions', 'all');
%! Xa = X;
%! for i = 1:rows (Mdl.IntMatrix)
%!   t = logical (Mdl.IntMatrix(i,:));
%!   Xt = X(:,t);
%!   Xi = ones (rows (X), 1);
%!   for c = 1:sum (t)
%!     Xi = Xi .* Xt(:,c);
%!   endfor
%!   Xa = [Xa, Xi];
%! endfor
%! yr = ones (rows (Xa), 1) * Mdl.ModelwInt.Intercept;
%! for j = 1:columns (Xa)
%!   yr = yr + ppval (Mdl.ModelwInt.Parameters(j), Xa(:,j));
%! endfor
%! [~, ySD] = predict (Mdl, X(1:4,:));
%! assert_equal (ySD, sqrt (var (Y - yr)) * ones (4, 1), 1e-12);

## A formula that selects fewer terms than there are predictors.  Asking this
## model for a standard deviation used to raise out of bound, so the second
## and third outputs of predict were unreachable on any formula-built model of
## this shape.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! Y = meas(:,4);
%! Mdl = fitrgam (X, Y, 'FitMethod', 'splines', 'Formula', 'Y ~ x1 + x2');
%! assert_equal (numel (Mdl.ModelwInt.Parameters), 2);
%! assert_equal (columns (Mdl.X), 3);
%! [~, ySD, yInt] = predict (Mdl, X(1:4,:));
%! assert_equal (ySD, 0.348625969903 * ones (4, 1), 1e-9);
%! assert_equal (size (yInt), [4, 2]);

## yFit never had this defect and must not acquire one: the prediction matrix
## was always reshaped, and still is.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! Y = meas(:,4);
%! Mdl = fitrgam (X, Y, 'FitMethod', 'splines', 'Interactions', 'all');
%! yA = predict (Mdl, X(1:3,:));
%! [yB, ~] = predict (Mdl, X(1:3,:));
%! assert_equal (yA, yB);
%! assert_equal (yA, [0.255203457138; 0.210404303824; 0.162599243283], 1e-9);

## The boosted-tree engine is reachable by name while the default is still the
## spline engine, and it reports the surface MATLAB reports.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,2:4), meas(:,1), 'FitMethod', 'boostedtrees');
%! assert_equal (Mdl.FitMethod, 'boostedtrees');
%! assert_equal (numel (Mdl.BinEdges), 3);
%! assert_equal (numel (fieldnames (Mdl.ModelParameters)), 13);
%! assert_equal (Mdl.ModelParameters.Type, 'regression');

## A regression intercept is the response mean and boosting leaves it there,
## where a classifier's is fitted and moves.  The asymmetry is real and no
## MATLAB-facing property reports it, so it is pinned here.
%!test
%! load fisheriris
%! y = meas(:,1);
%! Mdl = fitrgam (meas(:,2:4), y, 'FitMethod', 'boostedtrees');
%! assert_equal (Mdl.Intercept, mean (y), 1e-12);

## Interactions are detected on their own coarse grid.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,2:4), meas(:,1), 'FitMethod', 'boostedtrees', ...
%!                'Interactions', 'all');
%! assert_equal (rows (Mdl.Interactions), 3);
%! assert_equal (numel (Mdl.PairDetectionBinEdges{1}), 7);

## A tree-fitted model predicts, and compact and loadmodel predict alike.
%!test
%! load fisheriris
%! X = meas(:,2:4);
%! Mdl = fitrgam (X, meas(:,1), 'FitMethod', 'boostedtrees');
%! CMdl = compact (Mdl);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (predict (M2, X), predict (Mdl, X));

## The spline engine is unchanged and still reachable by name.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,2:4), meas(:,1), 'FitMethod', 'splines');
%! assert_equal (Mdl.FitMethod, 'splines');
%! assert_equal (Mdl.BinEdges, {});
%! assert_equal (Mdl.Knots, [5, 5, 5]);

## A standard deviation is a spline-engine capability: the boosted-tree engine
## fits no second model for it and says so rather than returning a wrong one.
%!error<RegressionGAM.predict: a standard deviation is only available from a model fitted with 'FitMethod' 'splines'.> ...
%! load fisheriris
%! Mdl = fitrgam (meas(:,2:4), meas(:,1), 'FitMethod', 'boostedtrees');
%! [y, ySD] = predict (Mdl, meas(1:4,2:4));

## An argument belonging to the other engine is refused, not ignored.
%!error<RegressionGAM: 'FitMethod' must be 'boostedtrees' or 'splines'.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'nonsense')
%!error<RegressionGAM: 'knots' is a parameter of the spline engine and cannot be used with 'FitMethod' 'boostedtrees'.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'boostedtrees', 'Knots', 4)
%!error<RegressionGAM: 'maxpvalue' is a parameter of the boosted-tree engine and cannot be used with 'FitMethod' 'splines'.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'splines', 'MaxPValue', 0.5)
%!error<RegressionGAM: 'NumTreesPerPredictor' must be a positive integer value.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'boostedtrees', ...
%!         'NumTreesPerPredictor', 0)
%!error<RegressionGAM: 'NumTreesPerInteraction' must be a positive integer value.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'boostedtrees', ...
%!         'NumTreesPerInteraction', 1.5)
%!error<RegressionGAM: 'MaxNumSplitsPerPredictor' must be a positive integer value.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'boostedtrees', ...
%!         'MaxNumSplitsPerPredictor', -1)
%!error<RegressionGAM: 'MaxNumSplitsPerInteraction' must be a positive integer value.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'boostedtrees', ...
%!         'MaxNumSplitsPerInteraction', 'a')
%!error<RegressionGAM: 'InitialLearnRateForPredictors' must be greater than 0 and at most 1.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'boostedtrees', ...
%!         'InitialLearnRateForPredictors', 0)
%!error<RegressionGAM: 'InitialLearnRateForInteractions' must be greater than 0 and at most 1.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'boostedtrees', ...
%!         'InitialLearnRateForInteractions', 2)
%!error<RegressionGAM: 'Verbose' must be a non-negative integer value.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'boostedtrees', 'Verbose', -1)
%!error<RegressionGAM: 'NumPrint' must be a positive integer value.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'boostedtrees', 'NumPrint', 0)
%!error<RegressionGAM: 'MaxPValue' must be between 0 and 1.> ...
%! fitrgam ([1;2;3;4], [1;2;3;4], 'FitMethod', 'boostedtrees', 'MaxPValue', 2)

## HyperparameterOptimizationResults is declared for MATLAB compatibility and
## stays empty, this class running no search over its hyperparameters.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,1:3), meas(:,4));
%! assert_equal (isempty (Mdl.HyperparameterOptimizationResults), true);

## Every documented response transform reaches the response that is reported.
%!test
%! load fisheriris
%! Mdl = fitrgam (meas(:,2:4), meas(:,1));
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
%! Mdl = fitrgam (meas(:,2:4), meas(:,1));
%! Mdl.ResponseTransform = 'none';
%! raw = predict (Mdl, meas([1, 60, 120],2:4));
%! Mdl.ResponseTransform = @(x) x .^ 2;
%! yhat = predict (Mdl, meas([1, 60, 120],2:4));
%! assert_equal (yhat, raw .^ 2, 1e-12);

## A table at loss
%!test  # the response is named, left out, or given beside the table
%! load fisheriris
%! X = meas(:,2:3);
%! y = meas(:,1);
%! T = table (X(:,1), X(:,2), 'VariableNames', {'SW', 'PL'});
%! T.SL = y;
%! Mdl = fitrgam (T, 'SL');
%! a = loss (Mdl, X, y);
%! assert_equal (loss (Mdl, T(:,1:2), y), a);
%! assert_equal (loss (Mdl, T, 'SL'), a);
%! assert_equal (loss (Mdl, T), a);

## Any two of 'Knots', 'Order' and 'DoF' determine the third, whichever
## order they are given in
%!test
%! X = linspace (0, 1, 50)';
%! Y = sin (6 * X);
%! M = fitrgam (X, Y, 'FitMethod', 'splines', 'Knots', 5, 'DoF', 9);
%! assert_equal ([M.Knots, M.Order, M.DoF], [5, 4, 9]);
%! M = fitrgam (X, Y, 'FitMethod', 'splines', 'DoF', 9, 'Knots', 5);
%! assert_equal ([M.Knots, M.Order, M.DoF], [5, 4, 9]);
%! M = fitrgam (X, Y, 'FitMethod', 'splines', 'Order', 2, 'DoF', 9);
%! assert_equal ([M.Knots, M.Order, M.DoF], [7, 2, 9]);
