## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
## Copyright (C) 2024-2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef ClassificationGAM
  ## -*- texinfo -*-
  ## @deftp {statistics} ClassificationGAM
  ##
  ## Generalized additive model classification
  ##
  ## The @code{ClassificationGAM} class implements a gradient boosting
  ## algorithm for classification.  This approach allows the model to capture
  ## non-linear relationships between predictors and the binary response
  ## variable.
  ##
  ## Generalized additive model classification is a statistical method that
  ## extends linear models by allowing non-linear relationships between each
  ## predictor and the response variable through smooth functions.  It combines
  ## the interpretability of linear models with the flexibility of
  ## non-parametric methods.
  ##
  ## Create a @code{ClassificationGAM} object by using the @code{fitcgam}
  ## function or the class constructor.
  ##
  ## Two weak learners are available, selected by @code{FitMethod}.
  ##
  ## @qcode{'boostedtrees'}, the default, boosts one shallow decision tree per
  ## predictor in each round, which is the scheme MATLAB's generalized
  ## additive model uses.  A second phase then boosts trees over pairs of
  ## predictors, where interactions are asked for.
  ##
  ## @qcode{'splines'} boosts a smoothing spline per predictor over
  ## @code{NumIterations} passes.  It has no MATLAB counterpart and is an
  ## Octave extension, kept because a smooth additive fit is a genuinely
  ## different and often better answer than a staircase of stumps.
  ##
  ## The two take different arguments, and an argument meant for one is
  ## refused by the other rather than ignored.
  ##
  ## The choice is visible in the properties.  @code{Knots}, @code{Order},
  ## @code{DoF}, @code{Formula}, @code{LearningRate}, @code{NumIterations},
  ## @code{BaseModel}, @code{ModelwInt} and @code{IntMatrix} describe a spline
  ## fit and are empty under the boosted-tree engine, while
  ## @code{ModelParameters}, @code{ReasonForTermination}, @code{BinEdges},
  ## @code{PairDetectionBinEdges} and @code{TreeModel} describe a tree fit and
  ## are empty under the spline engine.
  ##
  ## Fitted values are not expected to equal MATLAB's even under
  ## @qcode{'boostedtrees'}.  The stopping rule and the step-reduction limit
  ## are not recoverable from anything MATLAB reports, so this engine
  ## documents its own; what the two share is the estimator and the reported
  ## surface, not the arithmetic.
  ##
  ## @seealso{fitcgam}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix containing the unstandardized predictor data.  Each
    ## column of @var{X} represents one predictor (variable), and each row
    ## represents one observation.  This property is read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} Y
    ##
    ## Class labels
    ##
    ## Specified as a logical or numeric column vector, or as a character array
    ## or a cell array of character vectors with the same number of rows as the
    ## predictor data.  Each row in @var{Y} is the observed class label for
    ## the corresponding row in @var{X}.  This property is read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer value specifying the number of observations in the
    ## training dataset used for training the ClassificationGAM model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} RowsUsed
    ##
    ## Rows used for fitting
    ##
    ## A logical column vector with the same length as the observations in the
    ## original predictor data @var{X}, true for each row that was used for
    ## fitting the ClassificationGAM model.  It is empty, @qcode{[]},
    ## when every observation was used, so a non-empty value means that rows
    ## holding missing values were dropped.  This property is read-only.
    ##
    ## @end deftp
    RowsUsed        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer value specifying the number of predictors in the
    ## training dataset used for training the ClassificationGAM model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumPredictors   = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which they appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames  = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName    = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} ClassNames
    ##
    ## Names of classes in the response variable
    ##
    ## An array of unique values of the response variable @var{Y}, which has the
    ## same data types as the data in @var{Y}.  This property is read-only.
    ## @qcode{ClassNames} can have any of the following datatypes:
    ##
    ## @itemize
    ## @item Cell array of character vectors
    ## @item Character array
    ## @item Logical vector
    ## @item Numeric vector
    ## @end itemize
    ##
    ## @end deftp
    ClassNames      = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} Prior
    ##
    ## Prior probability for each class
    ##
    ## A 2-element numeric vector specifying the prior probabilities for each
    ## class.  The order of the elements in @qcode{Prior} corresponds to the
    ## order of the classes in @qcode{ClassNames}.  This property is read-only.
    ##
    ## Specified as a row vector with one entry per class, in the order of
    ## @qcode{ClassNames}, and rescaled to sum to one.  It may be given as
    ## @qcode{'empirical'}, @qcode{'uniform'}, a numeric vector, or a
    ## structure with @qcode{ClassNames} and @qcode{ClassProbs} fields, which
    ## assigns each probability by class name rather than by position.
    ##
    ## @end deftp
    Prior           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} Formula
    ##
    ## Model specification formula
    ##
    ## A character vector specifying the model formula in the form
    ## @qcode{'Y ~ terms'} where @qcode{Y} represents the response variable and
    ## @qcode{terms} specifies the predictor variables and interaction terms.
    ## This property is read-only.
    ##
    ## @end deftp
    Formula         = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} Interactions
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
    Interactions    = zeros (0, 2);

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} Knots
    ##
    ## Knots for spline fitting
    ##
    ## A scalar or row vector specifying the number of knots for each predictor
    ## variable in the spline fitting.  This property is read-only.
    ##
    ## @end deftp
    Knots           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} Order
    ##
    ## Order of spline fitting
    ##
    ## A scalar or row vector specifying the order of the spline for each
    ## predictor variable.  This property is read-only.
    ##
    ## @end deftp
    Order           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} DoF
    ##
    ## Degrees of freedom for spline fitting
    ##
    ## A scalar or row vector specifying the degrees of freedom for each
    ## predictor variable in the spline fitting.  This property is read-only.
    ##
    ## @end deftp
    DoF             = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} LearningRate
    ##
    ## Learning rate for gradient boosting
    ##
    ## A scalar value between 0 and 1 specifying the learning rate used in the
    ## gradient boosting algorithm.  This property is read-only.
    ##
    ## @end deftp
    LearningRate    = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} NumIterations
    ##
    ## Maximum number of iterations
    ##
    ## A positive integer specifying the maximum number of iterations for the
    ## gradient boosting algorithm.  This property is read-only.
    ##
    ## @end deftp
    NumIterations   = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} Intercept
    ##
    ## Intercept of the fitted model
    ##
    ## A numeric scalar, the log-odds of the response mean, which every
    ## additive term is measured against.  This property is read-only.
    ##
    ## @end deftp
    Intercept       = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector with one entry per observation used for
    ## training, normalised to sum to one.  This property is read-only.
    ##
    ## Each class carries its prior spread evenly over its own observations,
    ## so an observation of a class weighs @qcode{Prior} for that class
    ## divided by the number of observations it holds.
    ##
    ## @end deftp
    W               = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector holding the column of each predictor treated as
    ## categorical, and empty when none is.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} ExpandedPredictorNames
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
    ## @deftp {ClassificationGAM} {property} BaseModel
    ##
    ## Base model parameters
    ##
    ## A structure containing the parameters of the base model without any
    ## interaction terms.  The base model represents the generalized additive
    ## model with only the main effects (predictor terms) included.
    ## This property is read-only.
    ##
    ## @end deftp
    BaseModel = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} ModelwInt
    ##
    ## Model parameters with interactions
    ##
    ## A structure containing the parameters of the model that includes
    ## interaction terms.  This model extends the base model by adding
    ## interaction terms between predictors.  This property is read-only.
    ##
    ## @end deftp
    ModelwInt = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} IntMatrix
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
    IntMatrix = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} BinEdges
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
    ## @deftp {ClassificationGAM} {property} PairDetectionBinEdges
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
    ## @deftp {ClassificationGAM} {property} ModelParameters
    ##
    ## Parameters the model was fitted with
    ##
    ## A structure holding the fitting parameters.  Under the boosted-tree
    ## engine it carries MATLAB's own fields: @qcode{NumPrint},
    ## @qcode{MaxPValue}, @qcode{InitialLearnRateForPredictors},
    ## @qcode{InitialLearnRateForInteractions},
    ## @qcode{NumTreesPerPredictor}, @qcode{NumTreesPerInteraction},
    ## @qcode{MaxNumSplitsPerPredictor}, @qcode{MaxNumSplitsPerInteraction},
    ## @qcode{VerbosityLevel}, @qcode{Interactions}, @qcode{Version},
    ## @qcode{Method} and @qcode{Type}.  @qcode{Interactions} here is the
    ## request as it was made, a count or @qcode{'all'}, where the
    ## @qcode{Interactions} property of the model is the pairs actually
    ## selected.
    ##
    ## Under the spline engine it describes that scheme instead, carrying
    ## @qcode{Knots}, @qcode{Order}, @qcode{DoF}, @qcode{Formula},
    ## @qcode{Interactions}, @qcode{LearningRate} and @qcode{NumIterations},
    ## since none of the tree vocabulary applies to it.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ModelParameters = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} ReasonForTermination
    ##
    ## Why each fitting phase stopped
    ##
    ## A structure with the fields @qcode{PredictorTrees} and
    ## @qcode{InteractionTrees}, each a character vector saying why that phase
    ## of the fit ended: that it trained the trees it was asked for, or that
    ## it could no longer improve the model.  A phase that never ran reports
    ## an empty character vector, which is what a model with no interaction
    ## terms shows for the second field.
    ##
    ## It is empty under the spline engine, which has no tree budget to
    ## exhaust.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ReasonForTermination = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} FitMethod
    ##
    ## Which engine fitted the model
    ##
    ## A character vector, either @qcode{'boostedtrees'} or
    ## @qcode{'splines'}.  The default is @qcode{'boostedtrees'}, which is the
    ## scheme MATLAB's generalized additive model uses and the one the
    ## tree-shaped properties above describe.
    ##
    ## @qcode{'splines'} selects the penalised-spline engine instead, which is
    ## an Octave extension with no MATLAB counterpart.  It is the scheme this
    ## class fitted before version 1.9.0, and it is kept because a smooth
    ## additive fit is a genuinely different and often better answer than a
    ## staircase of stumps.  The two engines take different arguments and an
    ## argument meant for one is refused by the other rather than ignored.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    FitMethod = 'boostedtrees';

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} TreeModel
    ##
    ## The fitted shape functions and interaction surfaces
    ##
    ## A structure holding what the boosted-tree engine fitted, with fields
    ## @qcode{ShapeValues}, one column vector per predictor giving that
    ## predictor's contribution in each of its bins, @qcode{PairValues}, one
    ## matrix per selected pair, and @qcode{Pairs}, the predictor indices those
    ## matrices belong to.  A shape function is a step function, so these are
    ## the whole of the fit however many trees produced them.
    ##
    ## MATLAB exposes no equivalent: it reports the bin edges but never the
    ## values on them, so its shape functions can only be reached through
    ## @code{predict}.  This property is an Octave extension, and it is empty
    ## under the spline engine, whose fit lives in @code{BaseModel} and
    ## @code{ModelwInt}.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    TreeModel = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} HyperparameterOptimizationResults
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
    ## @deftp {ClassificationGAM} {property} Cost
    ##
    ## Cost of Misclassification
    ##
    ## A square matrix specifying the cost of misclassification of a point.
    ## @qcode{Cost(i,j)} is the cost of classifying a point into class @qcode{j}
    ## if its true class is @qcode{i} (that is, the rows correspond to the true
    ## class and the columns correspond to the predicted class).  The order of
    ## the rows and columns in @qcode{Cost} corresponds to the order of the
    ## classes in @qcode{ClassNames}.  The number of rows and columns in
    ## @qcode{Cost} is the number of unique classes in the response.  By
    ## default, @qcode{Cost(i,j) = 1} if @qcode{i != j}, and
    ## @qcode{Cost(i,j) = 0} if @qcode{i = j}.  In other words, the cost is 0
    ## for correct classification and 1 for incorrect classification.
    ##
    ## Add or change the @qcode{Cost} property using dot notation as in:
    ## @itemize
    ## @item @qcode{@var{obj}.Cost = @var{costMatrix}}
    ## @end itemize
    ##
    ##
    ## A cost may also be given as a struct with the fields
    ## @qcode{ClassNames} and @qcode{ClassificationCosts}, which names the
    ## order its own matrix is written in.  That matrix is permuted into the
    ## order of @qcode{ClassNames} above, so a caller need not know which
    ## order the classes were sorted into.  It must name every class.
    ##
    ## A cost must be floating point, not sparse, not complex, non-negative
    ## and zero down its diagonal, and must hold no @qcode{NaN} or
    ## @qcode{Inf}.  A @code{single} is widened to @code{double}.
    ## @end deftp
    Cost            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationGAM} {property} ScoreTransform
    ##
    ## Transformation function for classification scores
    ##
    ## Specified as a function handle for transforming the classification
    ## scores.  Add or change the @qcode{ScoreTransform} property using dot
    ## notation as in:
    ##
    ## @itemize
    ## @item @qcode{@var{obj}.ScoreTransform = 'function_name'}
    ## @item @qcode{@var{obj}.ScoreTransform = @@function_handle}
    ## @end itemize
    ##
    ## When specified as a character vector, it can be any of the following
    ## built-in functions.  Nevertheless, the @qcode{ScoreTransform} property
    ## always stores their function handle equivalent.
    ##
    ## @multitable @columnfractions 0.2 0.75
    ## @headitem @var{Value} @tab @var{Description}
    ## @item @qcode{'doublelogit'} @tab @math{1 ./ (1 + exp (-2 * x))}
    ## @item @qcode{'invlogit'} @tab @math{log (x ./ (1 - x))}
    ## @item @qcode{'ismax'} @tab Sets the score for the class with the
    ## largest score to 1, and for all other classes to 0
    ## @item @qcode{'logit'} @tab @math{1 ./ (1 + exp (-x))}
    ## @item @qcode{'none'} @tab @math{x} (no transformation)
    ## @item @qcode{'identity'} @tab @math{x} (no transformation)
    ## @item @qcode{'sign'} @tab
    ## @math{-1 for x < 0, 0 for x = 0, 1 for x >
    ## 0}
    ## @item @qcode{'symmetric'} @tab @math{2 * x - 1}
    ## @item @qcode{'symmetricismax'} @tab Sets the score for the class
    ## with the largest score to 1, and for all other classes to -1
    ## @item @qcode{'symmetriclogit'} @tab @math{2 ./ (1 + exp (-x)) - 1}
    ## @end multitable
    ##
    ## The default is @qcode{'logit'}, as in MATLAB.  This model's raw
    ## score is a log-odds, reported as the pair @math{[-f, f]} whose two
    ## columns sum to zero, and the transform is what turns it into the
    ## posterior probabilities that sum to one.  Every transform therefore
    ## composes on the log-odds and not on the probabilities, so
    ## @qcode{'none'} returns the log-odds themselves.
    ##
    ## @end deftp
    ScoreTransform  = 'logit';

  endproperties

  ## Readable by the counterpart class, which copies it, and kept out of
  ## the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    STfun = @(x) 1 ./ (1 + exp (-x));

    ## How many trees each boosting phase actually fitted, which the
    ## budget in ModelParameters does not say: a phase may stop early.
    ## Hidden because MATLAB reports it on the partitioned classes and
    ## not on the model, and that is where this is read from.
    NumTrainedTrees = [];
  endproperties

  ## Set methods for the properties a user may assign.
  methods (Hidden)

    function this = set.Cost (this, val)
      gnY = this.ClassNames;
      if (isempty (val))
        this.Cost = cast (! eye (classCount (gnY)), 'double');
      else
        ## Everything a cost must be, and the struct form, which
        ## is permuted into this model's class order.
        [val, errmsg] = costMatrix (val, gnY);
        if (! isempty (errmsg))
          error ("ClassificationGAM: %s", errmsg);
        endif
        this.Cost = val;
      endif
    endfunction

    function this = set.ScoreTransform (this, val)
      [f, nm] = parseScoreTransform (val, 'ClassificationGAM');
      this.ScoreTransform = nm;
      this.STfun = f;
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
      fprintf ("\n  ClassificationGAM\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      if (iscellstr (this.ClassNames))
        str = repmat ({'''%s'''}, 1, numel (this.ClassNames));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.ClassNames{:});
      elseif (ischar (this.ClassNames))
        str = repmat ({'''%s'''}, 1, rows (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, cellstr (this.ClassNames){:});
      else # single, double, logical
        str = repmat ({'%d'}, 1, numel (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, this.ClassNames);
      endif
      fprintf ("%+25s: %s\n", 'ClassNames', str);
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumPredictors', this.NumPredictors);
      if (! isempty (this.Formula))
        fprintf ("%+25s: '%s'\n", 'Formula', this.Formula);
      endif
      if (! isempty (this.Interactions))
        fprintf ("%+25s: [%dx%d %s]\n", 'Interactions', ...
                 size (this.Interactions, 1), size (this.Interactions, 2), ...
                 class (this.Interactions));
      endif
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {statistics} {@var{obj} =} ClassificationGAM (@var{X}, @var{Y})
    ## @deftypefnx {statistics} {@var{obj} =} ClassificationGAM (@dots{}, @var{name}, @var{value})
    ##
    ## Create a @qcode{ClassificationGAM} class object containing a generalized
    ## additive classification model.
    ##
    ## @code{@var{obj} = ClassificationGAM (@var{X}, @var{Y})} returns
    ## a ClassificationGAM object, with @var{X} as the predictor data
    ## and @var{Y} containing the class labels of observations in @var{X}.
    ##
    ## @itemize
    ## @item
    ## @code{X} must be a @math{N*P} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.  @var{X} will be used to train the GAM model.
    ## @item
    ## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}.  @var{Y} can contain any type
    ## of categorical data. @var{Y} must have the same number of rows as
    ## @var{X}.
    ## @end itemize
    ##
    ## @code{@var{obj} = ClassificationGAM (@dots{}, @var{name},
    ## @var{value})} returns a ClassificationGAM object with parameters
    ## specified by the following @qcode{@var{name}, @var{value}} paired input
    ## arguments:
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'PredictorNames'} @tab A cell array of character
    ## vectors specifying the names of the predictors. The length of this array
    ## must match the number of columns in @var{X}.
    ##
    ## @item @qcode{'ResponseName'} @tab A character vector specifying the
    ## name of the response variable.
    ##
    ## @item @qcode{'ClassNames'} @tab Names of the classes in the class
    ## labels, @var{Y}, used for fitting the GAM model.
    ## @qcode{ClassNames} are of the same type as the class labels in @var{Y}.
    ##
    ## @item @qcode{'Cost'} @tab An @math{N*R} numeric matrix containing
    ## misclassification cost for the corresponding instances in @var{X}, where
    ## @math{R} is the number of unique categories in @var{Y}.  If an instance
    ## is correctly classified into its category the cost is calculated to be 1,
    ## otherwise 0. The cost matrix can be altered by using
    ## @code{@var{Mdl}.cost = somecost}.  By default, its value is
    ## @qcode{@var{cost} = ones (rows (X), numel (unique (Y)))}.
    ##
    ## @item @qcode{'Prior'} @tab A numeric vector specifying the prior
    ## probabilities for each class.  The order of the elements in @qcode{Prior}
    ## corresponds to the order of the classes in @qcode{ClassNames}.
    ## Alternatively, you can specify @qcode{'empirical'} to use the empirical
    ## class probabilities or @qcode{'uniform'} to assume equal class
    ## probabilities.
    ##
    ## @item @qcode{'ScoreTransform'} @tab A user-defined function handle
    ## or a character vector specifying one of the following builtin functions
    ## specifying the transformation applied to predicted classification scores.
    ## Supported values include @qcode{'doublelogit'}, @qcode{'invlogit'},
    ## @qcode{'ismax'}, @qcode{'logit'}, @qcode{'none'}, @qcode{'identity'},
    ## @qcode{'sign'}, @qcode{'symmetric'}, @qcode{'symmetricismax'}, and
    ## @qcode{'symmetriclogit'}.
    ##
    ## @item @qcode{'Formula'} @tab (spline option) A character vector
    ## specifying the model
    ## formula in the form @qcode{'Y ~ terms'} where @qcode{Y} represents the
    ## response variable and @qcode{terms} specifies the predictor variables and
    ## interaction terms.
    ##
    ## @item @qcode{'Interactions'} @tab A logical matrix, a positive
    ## integer scalar, or the string @qcode{'all'} for defining the interactions
    ## between predictor variables.
    ##
    ## @item @qcode{'Knots'} @tab (spline option) A scalar or row vector
    ## specifying the
    ## number of knots for each predictor variable in the spline fitting.
    ##
    ## @item @qcode{'Order'} @tab (spline option) A scalar or row vector
    ## specifying the
    ## order of the spline for each predictor variable.
    ##
    ## @item @qcode{'DoF'} @tab (spline option) A scalar or row vector
    ## specifying the
    ## degrees of freedom for each predictor variable in the spline fitting.
    ##
    ## @item @qcode{'LearningRate'} @tab (spline option) A scalar value between
    ## 0 and 1
    ## specifying the learning rate used in the gradient boosting algorithm.
    ##
    ## @item @qcode{'NumIterations'} @tab (spline option) A positive integer
    ## specifying
    ## the maximum number of iterations for the gradient boosting algorithm.
    ## @end multitable
    ##
    ## A row marked @qcode{(spline option)} belongs to the spline
    ## engine and requires @qcode{'FitMethod', 'splines'}; passing one
    ## under the default boosted-tree engine is an error rather than
    ## being ignored.  The boosted-tree engine's own options are
    ## documented under @code{fitcgam}.
    ##
    ## @seealso{fitcgam}
    ## @end deftypefn
    function this = ClassificationGAM (X, Y, varargin)

      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationGAM: too few input arguments.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error ("ClassificationGAM: number of rows in X and Y must be equal.");
      endif

      nsample = rows (X);
      ndims_X = columns (X);

      ## Assign original X and Y data
      this.X = X;
      this.Y = Y;

      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);

      ## Set default values before parsing optional parameters
      PredictorNames = {};
      ResponseName   = [];
      Formula        = [];
      Interactions   = [];
      ClassNames     = [];
      Prior          = 'empirical';
      DoF            = ones (1, ndims_X) * 8;
      Order          = ones (1, ndims_X) * 3;
      Knots          = ones (1, ndims_X) * 5;
      LearningRate   = 0.1;
      NumIterations  = 100;
      Cost           = [];

      ## Boosted-tree defaults, MATLAB's own.  They are reported through
      ## ModelParameters, so they are part of the surface being matched and
      ## are not ours to improve: every other boosted additive model shrinks
      ## far harder than a step of 1, scikit-learn, gbm and mboost defaulting
      ## to 0.1 and the Explainable Boosting Machine to 0.015.  The docstring
      ## says so rather than the default being quietly changed.
      FitMethod                       = 'boostedtrees';
      NumTreesPerPredictor            = 300;
      NumTreesPerInteraction          = 100;
      MaxNumSplitsPerPredictor        = 1;
      MaxNumSplitsPerInteraction      = 4;
      InitialLearnRateForPredictors   = 1;
      InitialLearnRateForInteractions = 1;
      MaxPValue                       = 1;
      Verbose                         = 0;
      NumPrint                        = 10;

      ## Every name the caller asked for, so an argument meant for the other
      ## engine is refused instead of quietly doing nothing.
      namesGiven = {};

      ## Number of parameters for Knots, DoF, Order (maximum 2 allowed)
      KOD = 0;
      ## Number of parameters for Formula, Interactions (maximum 1 allowed)
      F_I = 0;

      ## Parse extra parameters
      while (numel (varargin) > 0)
        namesGiven{end+1} = tolower (varargin{1});
        switch (tolower (varargin {1}))

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("ClassificationGAM: 'PredictorNames'", ...
                             " must be supplied as a cellstring array."));
            elseif (numel (PredictorNames) != columns (X))
              error (strcat ("ClassificationGAM: 'PredictorNames'", ...
                             " must equal the number of columns in X."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("ClassificationGAM: 'ResponseName'", ...
                             " must be a character vector."));
            endif

          case 'classnames'
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames) ||
                   islogical (ClassNames) || ischar (ClassNames)))
              error (strcat ("ClassificationGAM: 'ClassNames' must be a", ...
                             " cell array of character vectors, a logical", ...
                             " vector, a numeric vector, or a character array."));
            endif
            ## Check that all class names are available in gnY
            if (iscellstr (ClassNames) || ischar (ClassNames))
              ClassNames = cellstr (ClassNames);
              if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                   ClassNames, 'UniformOutput', false))))
                error (strcat ("ClassificationGAM: not all 'ClassNames'", ...
                               " are present in Y."));
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, 'UniformOutput', false))))
                error (strcat ("ClassificationGAM: not all 'ClassNames'", ...
                               " are present in Y."));
              endif
            endif

          case 'prior'
            Prior = varargin{2};
            if (! (isstruct (Prior) || isnumeric (Prior) || ischar (Prior)))
              error (strcat ("ClassificationGAM: 'Prior' must be", ...
                             " a numeric vector or a string."));
            endif
            if (ischar (Prior) && ! any (strcmpi (Prior, {'empirical', 'uniform'})))
              error (strcat ("ClassificationGAM: 'Prior' must be", ...
                             " 'empirical', 'uniform', or a numeric vector."));
            endif
            if (isnumeric (Prior) && numel (Prior) != 2 && ! isstruct (Prior))
              error ("ClassificationGAM: 'Prior' must be a 2-element vector.");
            endif

          case 'cost'
            Cost = varargin{2};
            ## A struct carrying its own class order is a cost too,
            ## and is resolved by the property's own set method.
            if (! (isstruct (Cost)
                   || (isnumeric (Cost) && issquare (Cost))))
              error (strcat ("ClassificationGAM: 'Cost' must be", ...
                             " a numeric square matrix."));
            endif

          case 'scoretransform'
            name = 'ClassificationGAM';
            this.ScoreTransform = varargin{2};

          case 'formula'
            if (F_I < 1)
              Formula = varargin{2};
              if (! ischar (Formula) && ! islogical (Formula))
                error ("ClassificationGAM: 'Formula' must be a string.");
              endif
              F_I += 1;
            else
              error (strcat ("ClassificationGAM: 'Interactions'", ...
                             " have already been defined."));
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
                error ("ClassificationGAM: invalid 'Interactions' parameter.");
              endif
              F_I += 1;
            else
              error ("ClassificationGAM: 'Formula' has already been defined.");
            endif

          case 'knots'
            if (KOD < 2)
              Knots = varargin{2};
              if (! isnumeric (Knots) || ! (isscalar (Knots) ||
                  isequal (size (Knots), [1, ndims_X])))
                error ("ClassificationGAM: invalid value for 'Knots'.");
              endif
              DoF = Knots + Order;
              Order = DoF - Knots;
              KOD += 1;
            else
              error (strcat ("ClassificationGAM: 'DoF' and 'Order'", ...
                             " have been set already."));
            endif

          case 'order'
            if (KOD < 2)
              Order = varargin{2};
              if (! isnumeric (Order) || ! (isscalar (Order) ||
                  isequal (size (Order), [1, ndims_X])))
                error ("ClassificationGAM: invalid value for 'Order'.");
              endif
              DoF = Knots + Order;
              Knots = DoF - Order;
              KOD += 1;
            else
              error (strcat ("ClassificationGAM: 'DoF' and 'Knots'", ...
                             " have been set already."));
            endif

          case 'dof'
            if (KOD < 2)
              DoF = varargin{2};
              if (! isnumeric (DoF) ||
                  ! (isscalar (DoF) || isequal (size (DoF), [1, ndims_X])))
                error ("ClassificationGAM: invalid value for 'DoF'.");
              endif
              Knots = DoF - Order;
              Order = DoF - Knots;
              KOD += 1;
            else
              error (strcat ("ClassificationGAM: 'Knots' and 'Order'", ...
                             " have been set already."));
            endif

          case 'learningrate'
            LearningRate = varargin{2};
            if (LearningRate > 1 || LearningRate <= 0)
              error (strcat ("ClassificationGAM: 'LearningRate'", ...
                             " must be between 0 and 1."));
            endif

          case 'numiterations'
            NumIterations = varargin{2};
            if (! isnumeric (NumIterations) || NumIterations <= 0)
              error (strcat ("ClassificationGAM: 'NumIterations'", ...
                             " must be a positive integer value."));
            endif

          case 'fitmethod'
            FitMethod = varargin{2};
            if (! (ischar (FitMethod) && isrow (FitMethod)) ||
                ! any (strcmpi (FitMethod, {'boostedtrees', 'splines'})))
              error (strcat ("ClassificationGAM: 'FitMethod' must be", ...
                             " 'boostedtrees' or 'splines'."));
            endif
            FitMethod = tolower (FitMethod);

          case 'numtreesperpredictor'
            NumTreesPerPredictor = varargin{2};
            if (! isnumeric (NumTreesPerPredictor) ||
                ! isscalar (NumTreesPerPredictor) ||
                NumTreesPerPredictor < 1 ||
                fix (NumTreesPerPredictor) != NumTreesPerPredictor)
              error (strcat ("ClassificationGAM:", ...
                             " 'NumTreesPerPredictor' must be a positive", ...
                             " integer value."));
            endif

          case 'numtreesperinteraction'
            NumTreesPerInteraction = varargin{2};
            if (! isnumeric (NumTreesPerInteraction) ||
                ! isscalar (NumTreesPerInteraction) ||
                NumTreesPerInteraction < 1 ||
                fix (NumTreesPerInteraction) != NumTreesPerInteraction)
              error (strcat ("ClassificationGAM:", ...
                             " 'NumTreesPerInteraction' must be a", ...
                             " positive integer value."));
            endif

          case 'maxnumsplitsperpredictor'
            MaxNumSplitsPerPredictor = varargin{2};
            if (! isnumeric (MaxNumSplitsPerPredictor) ||
                ! isscalar (MaxNumSplitsPerPredictor) ||
                MaxNumSplitsPerPredictor < 1 ||
                fix (MaxNumSplitsPerPredictor) != MaxNumSplitsPerPredictor)
              error (strcat ("ClassificationGAM:", ...
                             " 'MaxNumSplitsPerPredictor' must be a", ...
                             " positive integer value."));
            endif

          case 'maxnumsplitsperinteraction'
            MaxNumSplitsPerInteraction = varargin{2};
            if (! isnumeric (MaxNumSplitsPerInteraction) ||
                ! isscalar (MaxNumSplitsPerInteraction) ||
                MaxNumSplitsPerInteraction < 1 ||
                fix (MaxNumSplitsPerInteraction) != MaxNumSplitsPerInteraction)
              error (strcat ("ClassificationGAM:", ...
                             " 'MaxNumSplitsPerInteraction' must be a", ...
                             " positive integer value."));
            endif

          case 'initiallearnrateforpredictors'
            InitialLearnRateForPredictors = varargin{2};
            if (! isnumeric (InitialLearnRateForPredictors) ||
                ! isscalar (InitialLearnRateForPredictors) ||
                InitialLearnRateForPredictors <= 0 ||
                InitialLearnRateForPredictors > 1)
              error (strcat ("ClassificationGAM:", ...
                             " 'InitialLearnRateForPredictors' must be", ...
                             " greater than 0 and at most 1."));
            endif

          case 'initiallearnrateforinteractions'
            InitialLearnRateForInteractions = varargin{2};
            if (! isnumeric (InitialLearnRateForInteractions) ||
                ! isscalar (InitialLearnRateForInteractions) ||
                InitialLearnRateForInteractions <= 0 ||
                InitialLearnRateForInteractions > 1)
              error (strcat ("ClassificationGAM:", ...
                             " 'InitialLearnRateForInteractions' must be", ...
                             " greater than 0 and at most 1."));
            endif

          case 'verbose'
            Verbose = varargin{2};
            if (! isnumeric (Verbose) || ! isscalar (Verbose) || Verbose < 0
                || fix (Verbose) != Verbose)
              error (strcat ("ClassificationGAM: 'Verbose' must be a", ...
                             " non-negative integer value."));
            endif

          case 'numprint'
            NumPrint = varargin{2};
            if (! isnumeric (NumPrint) || ! isscalar (NumPrint)
                || NumPrint < 1 || fix (NumPrint) != NumPrint)
              error (strcat ("ClassificationGAM: 'NumPrint' must be a", ...
                             " positive integer value."));
            endif

          case 'maxpvalue'
            MaxPValue = varargin{2};
            if (! isnumeric (MaxPValue) || ! isscalar (MaxPValue) ||
                MaxPValue < 0 || MaxPValue > 1)
              error (strcat ("ClassificationGAM: 'MaxPValue' must be", ...
                             " between 0 and 1."));
            endif

          otherwise
            error (strcat ("ClassificationGAM: invalid parameter", ...
                           " name in optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## An argument belongs to one engine or the other, and asking for one
      ## the chosen engine cannot honour is refused rather than ignored.  The
      ## alternative is the trap MATLAB sets with a name it accepts and never
      ## reads: the caller gets a fit that quietly disregarded what it asked
      ## for.
      splineOnly = {'knots', 'order', 'dof', 'formula', ...
                    'learningrate', 'numiterations'};
      treeOnly = {'numtreesperpredictor', 'numtreesperinteraction', ...
                  'maxnumsplitsperpredictor', 'maxnumsplitsperinteraction', ...
                  'initiallearnrateforpredictors', ...
                  'initiallearnrateforinteractions', 'maxpvalue', ...
                  'verbose', 'numprint'};
      if (strcmp (FitMethod, 'boostedtrees'))
        clash = intersect (namesGiven, splineOnly);
        if (! isempty (clash))
          error (strcat ("ClassificationGAM: '", clash{1}, "' is a", ...
                         " parameter of the spline engine and cannot be", ...
                         " used with 'FitMethod' 'boostedtrees'."));
        endif
      else
        clash = intersect (namesGiven, treeOnly);
        if (! isempty (clash))
          error (strcat ("ClassificationGAM: '", clash{1}, "' is a", ...
                         " parameter of the boosted-tree engine and cannot", ...
                         " be used with 'FitMethod' 'splines'."));
        endif
      endif

      ## Generate default predictors and response variable names (if necessary)
      if (isempty (PredictorNames))
        for i = 1:columns (X)
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = 'Y';
      endif

      ## Assign predictors and response variable names
      this.PredictorNames = PredictorNames;
      this.ResponseName   = ResponseName;

      ## Handle class names
      if (! isempty (ClassNames))
        ## Anything textual is matched as whole names, gnY being grp2idx's
        ## own cellstr of them.  A character matrix is not a cellstr, and
        ## ismember between two of them compares character by character, so
        ## it would answer a question nobody asked.
        if (iscellstr (ClassNames) || ischar (ClassNames))
          ru = find (! ismember (gnY, cellstr (ClassNames)));
        else
          ru = find (! ismember (glY, ClassNames));
        endif
        for i = 1:numel (ru)
          gY(gY == ru(i)) = NaN;
        endfor
      endif

      ## An observation is dropped only when its response is missing.  A row
      ## whose predictors hold missing values is kept and reported as used,
      ## while the fit below draws on the complete observations alone.
      RowsUsed  = ! isnan (gY);
      ## Index the rows and not the elements: a response naming its
      ## classes in the rows of a character matrix has one column per
      ## character, and a linear index flattens the names into single
      ## letters.  Every other accepted type is a column, for which
      ## the two forms agree.
      Yret      = Y(RowsUsed, :);
      Xret      = X(RowsUsed, :);
      this.X    = Xret;
      this.Y    = Yret;
      cobs      = ! any (isnan (Xret), 2);
      Y         = Yret(cobs, :);
      X         = Xret(cobs, :);

      ## Renew groups in Y.  The third output of grp2idx holds the levels in
      ## the type of Y, where the second is always a cell array of character
      ## vectors, so a numeric or logical response keeps its own type.  grp2idx
      ## orders a character or cell response by first appearance where MATLAB
      ## sorts the classes, so the levels are sorted and the indices remapped.
      [gY, gnY, glY] = grp2idx (Y);
      if (ischar (glY))
        [glY, sidx] = sortrows (glY);
      else
        [glY, sidx] = sort (glY);
      endif
      remap(sidx) = 1:numel (sidx);
      gY = remap(gY)(:);
      gnY = gnY(sidx);
      this.ClassNames = glY;

      ## Check that we are dealing only with binary classification
      if (numel (gnY) > 2)
        error ("ClassificationGAM: can only be used for binary classification.");
      endif

      ## Force Y into the 0/1 coding the fitter needs.  This used to be done
      ## only for a response that was not already numeric, so a numeric one
      ## reached the fitter as it was given: a response of 1 and 2 seeded the
      ## model with log (mean (Y) / (1 - mean (Y))) at a mean of 1.5, whose
      ## argument is negative, and the intercept came back NaN with every
      ## score along with it.  gY indexes the sorted class names, so gY - 1 is
      ## the coding for any binary response and reproduces a 0/1 one exactly.
      Y = gY - 1;

      this.NumObservations = rows (this.X);
      ## RowsUsed is left empty when every observation was used, as in MATLAB
      if (all (RowsUsed))
        this.RowsUsed = [];
      else
        this.RowsUsed = RowsUsed;
      endif

      ## Assign the number of original predictors to the ClassificationGAM object
      this.NumPredictors = ndims_X;

      ## Assign Cost and compute Prior
      this.Cost = Cost;
      if (isstruct (Prior))
        Prior = priorFromStruct (Prior, this.ClassNames, ...
                                 'ClassificationGAM');
      endif
      if (ischar (Prior))
        if (strcmpi (Prior, 'uniform'))
          this.Prior = [0.5, 0.5];
        elseif (strcmpi (Prior, 'empirical'))
          counts = histc (gY, 1:2);
          this.Prior = counts(:)' / sum (counts);
        endif
      else
        ## Numeric prior - normalize to sum to 1
        this.Prior = Prior(:)' / sum (Prior);
      endif

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
      this.Formula       = Formula;
      this.Interactions  = Interactions;
      this.Knots         = Knots;
      this.Order         = Order;
      this.DoF           = DoF;
      this.LearningRate  = LearningRate;
      this.NumIterations = NumIterations;

      this.FitMethod = FitMethod;

      ## Bookkeeping MATLAB reports alongside the fit, whichever engine ran
      this.W                     = priorWeights (this.Prior, gY, ...
                                                 this.NumObservations);
      this.CategoricalPredictors = [];
      this.ExpandedPredictorNames = this.PredictorNames;

      if (strcmp (FitMethod, 'boostedtrees'))
        ## The spline parameters describe a scheme that did not run, so they
        ## are left empty rather than reporting numbers nothing used.
        this.Knots         = [];
        this.Order         = [];
        this.DoF           = [];
        this.LearningRate  = [];
        this.NumIterations = [];
        ## The engine wants the response as zeros and ones.  Y itself is
        ## only converted when it was not numeric to begin with, so a numeric
        ## label of 1 and 2 would arrive unconverted; the group index is the
        ## encoding that is always right.
        this = this.fitBoosted (X, gY(:) - 1, Interactions, ...
                                NumTreesPerPredictor, ...
                                NumTreesPerInteraction, ...
                                MaxNumSplitsPerPredictor, ...
                                MaxNumSplitsPerInteraction, ...
                                InitialLearnRateForPredictors, ...
                                InitialLearnRateForInteractions, MaxPValue, ...
                                Verbose, NumPrint);
      else
        ## Fit the basic model
        Inter = mean (Y);
        [iter, param, res, RSS, intercept] = this.fitGAM (X, Y, Inter, ...
                                                          Knots, Order, ...
                                                          LearningRate, ...
                                                          NumIterations);
        this.BaseModel.Intercept  = intercept;
        this.BaseModel.Parameters = param;
        this.BaseModel.Iterations = iter;
        this.BaseModel.Residuals  = res;
        this.BaseModel.RSS        = RSS;
        this.Intercept            = intercept;

        ## Handle interaction terms (if given)
        if (F_I > 0)
          this = this.fitModelwInt (X, Y, Inter, Knots, Order, DoF, ...
                                    LearningRate, NumIterations);
        endif

        ## The property MATLAB reports is the two-way terms the fitted model
        ## carries, as predictor index pairs, whatever form they were asked
        ## for in.  The term matrix stays the complete record: it also holds
        ## the main effects a formula names and any term above two
        ## predictors, neither of which has a two-column form.
        this.Interactions = interactionPairs (this.IntMatrix);

        ## The spline scheme has no tree vocabulary to report, so its
        ## parameter struct describes itself instead.
        this.ModelParameters = struct ('Knots', this.Knots, ...
                                       'Order', this.Order, ...
                                       'DoF', this.DoF, ...
                                       'Formula', this.Formula, ...
                                       'Interactions', this.Interactions, ...
                                       'LearningRate', this.LearningRate, ...
                                       'NumIterations', this.NumIterations);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{obj} =} addInteractions (@var{obj}, @var{interactions})
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
    ## @seealso{fitcgam, ClassificationGAM}
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
        error (strcat ("ClassificationGAM.addInteractions: adding", ...
                       " interaction terms to a model that already", ...
                       " includes them is not supported."));
      endif

      if (! ((isnumeric (interactions) && isscalar (interactions)
              && interactions == fix (interactions) && interactions >= 0)
             || islogical (interactions)
             || (ischar (interactions) && strcmpi (interactions, 'all'))))
        error (strcat ("ClassificationGAM.addInteractions: invalid", ...
                       " 'Interactions' parameter."));
      endif

      ## Under the boosted-tree engine the interaction phase simply runs now,
      ## starting from the predictor phase the model already carries, which is
      ## the same order the constructor would have run them in.  The predictor
      ## phase is not refitted, so the main effects are untouched and a model
      ## with interactions added is the model the constructor would have built
      ## had it been asked for them.
      if (strcmp (this.FitMethod, 'boostedtrees'))
        cobs = ! any (isnan (this.X), 2);
        Xfit = this.X(cobs, :);
        [~, ~, gY] = uniqueLabels (this.Y(cobs, :));
        Yfit = gY(:) - 1;
        MP = this.ModelParameters;
        lrInter = MP.InitialLearnRateForInteractions;
        this = this.fitBoostedInteractions (Xfit, Yfit, interactions, ...
                                            MP.NumTreesPerInteraction, ...
                                            MP.MaxNumSplitsPerInteraction, ...
                                            lrInter, MP.MaxPValue);
        return;
      endif

      ## parseInteractions reads the specification from the property, which
      ## afterwards holds the pairs the fit settled on, exactly as the
      ## constructor leaves it.
      this.Interactions = interactions;
      this.IntMatrix = this.parseInteractions ();

      ## The fit sees the complete observations and the response in the
      ## coding the boosting works in, prepared as the constructor prepares
      ## them.  Knots, Order and DoF are held unexpanded on the object and
      ## are widened to the interaction columns by fitModelwInt.
      cobs = ! any (isnan (this.X), 2);
      Xfit = this.X(cobs, :);
      [~, ~, gY] = uniqueLabels (this.Y(cobs, :));
      Yfit = gY(:) - 1;
      this = this.fitModelwInt (Xfit, Yfit, mean (Yfit), this.Knots, ...
                                this.Order, this.DoF, this.LearningRate, ...
                                this.NumIterations);
      this.Interactions = interactionPairs (this.IntMatrix);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationGAM} {[@var{label}, @var{score}] =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationGAM} {[@var{label}, @var{score}] =} predict (@dots{}, @qcode{'IncludeInteractions'}, @var{includeInteractions})
    ##
    ## Predict labels for new data using the Generalized Additive Model (GAM)
    ## stored in a ClassificationGAM object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the predicted
    ## labels for the data in @var{XC} based on the model stored in the
    ## ClassificationGAM object, @var{obj}.
    ##
    ## @code{[@var{label}, @var{score}] = predict (@var{obj}, @var{XC})} also
    ## returns @var{score}, which contains the predicted class scores or
    ## posterior probabilities for each observation.
    ##
    ## @code{[@var{label}, @var{score}] = predict (@var{obj}, @var{XC},
    ## 'IncludeInteractions', @var{includeInteractions})} allows you to specify
    ## whether interaction terms should be included when making predictions.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationGAM} class object.
    ## @item
    ## @var{XC} must be an @math{M*P} numeric matrix where each row is an
    ## observation and each column corresponds to a predictor variable.
    ## @item
    ## @var{includeInteractions} is a logical scalar indicating whether to
    ## include interaction terms in the predictions.
    ## @end itemize
    ##
    ## @seealso{ClassificationGAM, fitcgam}
    ## @end deftypefn

    function [labels, scores] = predict (this, XC, varargin)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("ClassificationGAM.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("ClassificationGAM.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("ClassificationGAM.predict:", ...
                       " XC must have the same number of", ...
                       " predictors as the trained model."));
      endif

      ## Clean XC data
      notnansf  = ! logical (sum (isnan (XC), 2));
      XC        = XC(notnansf, :);

      ## Default values for Name-Value Pairs
      ## Which store holds the interaction terms depends on the engine: the
      ## spline scheme keeps them as extra columns described by IntMatrix,
      ## the boosted-tree scheme as surfaces over predictor pairs.
      hasInt = ! isempty (this.IntMatrix);
      if (strcmp (this.FitMethod, 'boostedtrees') && ! isempty (this.TreeModel))
        hasInt = ! isempty (this.TreeModel.Pairs);
      endif
      incInt = hasInt;
      Cost = this.Cost;

      ## Parse optional arguments
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'includeinteractions'
            tmpInt = varargin{2};
            if (! islogical (tmpInt) || (tmpInt != 0 && tmpInt != 1))
              error (strcat ("ClassificationGAM.predict:", ...
                             " includeinteractions must be a logical value."));
            endif
            ## Check model for interactions
            if (tmpInt && ! hasInt)
              error (strcat ("ClassificationGAM.predict: trained model", ...
                             " does not include any interactions."));
            endif
            incInt = tmpInt;

          otherwise
            error (strcat ("ClassificationGAM.predict: invalid NAME in", ...
                           " optional pairs of arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      ## The boosted-tree engine keeps its fit as step functions over bins,
      ## so a term is a lookup rather than a spline evaluation and the whole
      ## prediction is one call.  It shares everything after it: the cost
      ## matrix, the label, and the transform.
      ## An empty TreeModel means no tree fit is present, whatever FitMethod
      ## says: a default-constructed object has one, and so would a model
      ## saved before the engine existed.  Such an object falls through to
      ## the spline path rather than indexing into nothing.
      if (strcmp (this.FitMethod, 'boostedtrees') && ! isempty (this.TreeModel))
        ## Excluding the interactions means excluding the constant they
        ## handed the intercept as well.
        interc = this.Intercept;
        if (! incInt && isfield (this.TreeModel, 'PairIntercept'))
          interc = interc - this.TreeModel.PairIntercept;
        endif
        if (! incInt || isempty (this.TreeModel.Pairs))
          scores = gamboostpredict (this.BinEdges, ...
                                    this.TreeModel.ShapeValues, XC, ...
                                    interc);
        else
          scores = gamboostpredict (this.BinEdges, ...
                                    this.TreeModel.ShapeValues, XC, ...
                                    interc, 0, ...
                                    this.PairDetectionBinEdges, ...
                                    this.TreeModel.PairValues, ...
                                    this.TreeModel.Pairs);
        endif
        scores = [-scores, scores];

        post = 1 ./ (1 + exp (-scores));
        numObservations = size (XC, 1);
        CE = zeros (numObservations, 2);
        for k = 1:2
          for i = 1:2
            CE(:, k) = CE(:, k) + post(:, i) * Cost(k, i);
          endfor
        endfor
        [~, minIdx] = min (CE, [], 2);
        labels = labelsFromIndex (this.ClassNames, minIdx);
        scores = this.STfun (scores);
        return;
      endif

      ## Choose whether interactions must be included
      if (incInt)
        ## Which construction path the model took: an interaction
        ## list appends its terms to the predictors, a formula
        ## names every term the model has and replaces them.
        if (isempty (this.Formula))
          ## Append interaction terms to the predictor matrix
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = XC(:,tindex);
            Xinter = ones (rows (XC), 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append interaction terms
            XC = [XC, Xinter];
          endfor
        else
          ## Add selected predictors and interaction terms
          XN = [];
          for i = 1:rows (this.IntMatrix)
            tindex = logical (this.IntMatrix(i,:));
            Xterms = XC(:,tindex);
            Xinter = ones (rows (XC), 1);
            for c = 1:sum (tindex)
              Xinter = Xinter .* Xterms(:,c);
            endfor
            ## Append selected predictors and interaction terms
            XN = [XN, Xinter];
          endfor
          XC = XN;
        endif
        ## Get parameters and intercept vectors from model with interactions
        params = this.ModelwInt.Parameters;
        Interc = this.ModelwInt.Intercept;
      else
        ## Get parameters and intercept vectors from base model
        params = this.BaseModel.Parameters;
        Interc = this.BaseModel.Intercept;
      endif

      ## Predict the raw score from testing data
      scores = predict_val (params, XC, Interc);

      ## Expected misclassification cost is defined on the posteriors, which
      ## the score becomes only under the logistic link, so the label is
      ## taken from those and not from the score the caller is handed.
      post = 1 ./ (1 + exp (-scores));

      ## Compute the expected misclassification cost matrix
      numObservations = size (XC, 1);
      CE = zeros (numObservations, 2);

      for k = 1:2
        for i = 1:2
          CE(:, k) = CE(:, k) + post(:, i) * Cost(k, i);
        endfor
      endfor

      ## Select the class with the minimum expected misclassification cost
      [~, minIdx] = min (CE, [], 2);
      labels = labelsFromIndex (this.ClassNames, minIdx);

      ## Apply ScoreTransform
      scores = this.STfun (scores);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationGAM} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross Validate a Generalized Additive Model classification object.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} returns a cross-validated model
    ## object, @var{CVMdl}, from a trained model, @var{obj}, using 10-fold
    ## cross-validation by default.
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
    ## object used for cross-validation.  @code{"CVPartition", @var{cv}}, where
    ## @code{isa (@var{cv}, "cvpartition")} = 1.
    ##
    ## @end multitable
    ##
    ## @seealso{fitcgam, ClassificationGAM, cvpartition,
    ## ClassificationPartitionedModel}
    ## @end deftypefn

    function CVMdl = crossval (this, varargin)
      ## Check input
      if (nargin < 1)
        error ("ClassificationGAM.crossval: too few input arguments.");
      endif

      if (numel (varargin) == 1)
        error (strcat ("ClassificationGAM.crossval: Name-Value", ...
                       " arguments must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("ClassificationGAM.crossval: specify only", ...
                       " one of the optional Name-Value paired arguments."));
      endif

      ## Add default values
      if (this.NumObservations < 10)
        numFolds  = this.NumObservations;
      else
        numFolds  = 10;
      endif
      Holdout     = [];
      Leaveout    = 'off';
      CVPartition = [];

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'kfold'
            numFolds = varargin{2};
            if (! (isnumeric (numFolds) && isscalar (numFolds)
                   && (numFolds == fix (numFolds)) && numFolds > 1))
              error (strcat ("ClassificationGAM.crossval: 'KFold'", ...
                             " must be an integer value greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("ClassificationGAM.crossval: 'Holdout'", ...
                             " must be a numeric value between 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat ("ClassificationGAM.crossval: 'Leaveout'", ...
                             " must be either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("ClassificationGAM.crossval: 'CVPartition'",...
                             " must be a 'cvpartition' object."));
            endif

          otherwise
            error (strcat ("ClassificationGAM.crossval: invalid",...
                           " parameter name in optional paired arguments."));
          endswitch
        varargin(1:2) = [];
      endwhile

      ## Determine the cross-validation method to use.  The partition covers
      ## the observations actually trained on: a row dropped for a missing
      ## value is not one the folds can use, and including it would leave the
      ## partition, the stored data and NumObservations disagreeing.  The
      ## response is passed rather than a count so the folds stay stratified.
      Yused = this.Y;
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (Yused, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (this.NumObservations, 'LeaveOut');
      else
        partition = cvpartition (Yused, 'KFold', numFolds);
      endif

      ## Create a cross-validated model object
      CVMdl = ClassificationPartitionedModel (this, partition);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{CVMdl} =} compact (@var{obj})
    ##
    ## Create a CompactClassificationGAM object.
    ##
    ## @code{@var{CVMdl} = compact (@var{obj})} creates a compact version of the
    ## ClassificationGAM object, @var{obj}.
    ##
    ## @seealso{fitcgam, ClassificationGAM, CompactClassificationGAM}
    ## @end deftypefn
    function CVMdl = compact (this)
      ## Create a compact model
      CVMdl = CompactClassificationGAM (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin of a generalized additive model.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns a column
    ## vector holding, for each row of @var{X}, the score the model gives its
    ## true class in @var{Y} less the score it gives the other class.  A
    ## positive margin means the observation is classified correctly, and the
    ## larger it is the more confidently so.
    ##
    ## @seealso{ClassificationGAM, edge, loss, predict}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationGAM.margin: too few input arguments.");
      endif

      [X, Y] = checkXY_ (this, X, Y, "margin");

      [~, scores] = predict (this, X);
      classes = this.ClassNames;
      m = zeros (rows (X), 1);
      ## Resolve every observation's class once, rather than once per
      ## iteration: the lookup does not depend on i.
      [gYidx, ~] = labelIndices (classes, Y);
      for i = 1:rows (X)
        idx = gYidx(i);
        if (isempty (idx))
          m(i) = NaN;
          continue;
        endif
        true_score = scores(i, idx);
        scores(i, idx) = -Inf;
        m(i) = true_score - max (scores(i,:));
        scores(i, idx) = true_score;
      endfor


    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationGAM} {@var{e} =} edge (@dots{}, @qcode{"Weights"}, @var{w})
    ##
    ## Classification edge of a generalized additive model.
    ##
    ## @code{@var{e} = edge (@var{obj}, @var{X}, @var{Y})} returns the mean of
    ## the classification margins over the rows of @var{X}.
    ##
    ## @code{@var{e} = edge (@dots{}, @qcode{"Weights"}, @var{w})} takes the
    ## weighted mean instead, with one weight per row of @var{X}.
    ##
    ## @seealso{ClassificationGAM, margin, loss, predict}
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationGAM.edge: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationGAM.edge: Name-Value arguments", ...
                       " must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, "edge");

      ## The weights are normalized within each class to that class's prior,
      ## which is what the oracle does and is not the same as dividing by
      ## their total.  This used to divide by the total.
      ## The weights are parsed before anything is computed, so a bad
      ## Name-Value pair is reported as such rather than after a margin.
      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       "ClassificationGAM", "edge");
      m = margin (this, X, Y);
      e = sum (W .* m(:)) / sum (W);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationGAM} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of a generalized additive model.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the loss of
    ## the model on the rows of @var{X} against the true labels @var{Y}.
    ##
    ## @code{@var{L} = loss (@dots{}, @var{name}, @var{value})} accepts the
    ## following name-value pairs:
    ##
    ## @itemize
    ## @item
    ## @qcode{"LossFun"} selects the loss.  Supported values are
    ## @qcode{"mincost"}, the default, @qcode{"binodeviance"},
    ## @qcode{"classifcost"}, @qcode{"classiferror"}, @qcode{"exponential"},
    ## @qcode{"hinge"}, @qcode{"logit"} and @qcode{"quadratic"}.
    ## @qcode{"mincost"} assigns each observation to the class of least
    ## expected cost and charges what that assignment costs, so it reads the
    ## scores as a posterior, which is what this model returns;
    ## @qcode{"classifcost"} charges what the model's own prediction costs.
    ##
    ## @item
    ## @qcode{"Weights"} holds one weight per row of @var{X}, normalised to
    ## sum to one before it is applied.
    ## @end itemize
    ##
    ## @seealso{ClassificationGAM, margin, edge, predict}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationGAM.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationGAM.loss: Name-Value arguments", ...
                       " must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, "loss");

      ## Parse optional arguments
      LossFun = 'mincost';
      lossnames = {'binodeviance', 'classifcost', 'classiferror', ...
                   'exponential', 'hinge', 'logit', 'mincost', 'quadratic'};
      args = varargin;
      keep = true (1, numel (args));
      for i = 1:2:numel (args)
        if (strcmpi (args{i}, 'lossfun'))
          LossFun = args{i+1};
          if (! (ischar (LossFun) && isrow (LossFun)))
            error (strcat ("ClassificationGAM.loss: 'LossFun' must be", ...
                           " a character vector."));
          endif
          LossFun = tolower (LossFun);
          if (! any (strcmpi (LossFun, lossnames)))
            error ("ClassificationGAM.loss: unsupported Loss function.");
          endif
          keep(i:i+1) = false;
        endif
      endfor
      W = getWeights_ (this, args(keep), rows (X), "loss");
      W = W(:) / sum (W);

      [label, scores] = predict (this, X);
      classes = this.ClassNames;

      ## Membership of the true class, as an indicator per class
      Yind = zeros (rows (X), classCount (classes));
      ## Resolve every observation's class once, rather than once per
      ## iteration: the lookup does not depend on i.
      [gYidx, ~] = labelIndices (classes, Y);
      for i = 1:rows (X)
        idx = gYidx(i);
        if (isempty (idx))
          L = NaN;
          return;
        endif
        Yind(i, idx) = 1;
      endfor

      ## The scalar score of the true class of each observation
      mj = sum (scores .* Yind, 2);

      switch (LossFun)
        case 'classiferror'
          wrong = zeros (rows (X), 1);
          for i = 1:rows (X)
            wrong(i) = ! isequal (label(i), Y(i));
          endfor
          L = sum (W .* wrong);
        case 'binodeviance'
          L = sum (W .* log (1 + exp (-2 * mj)));
        case 'hinge'
          L = sum (W .* max (0, 1 - mj));
        case 'exponential'
          L = sum (W .* exp (-mj));
        case 'logit'
          L = sum (W .* log (1 + exp (-mj)));
        case 'quadratic'
          L = sum (W .* (1 - mj) .^ 2);
        case 'mincost'
          ## Each observation is assigned to the class of least expected
          ## cost, and charged what that assignment actually costs given its
          ## true class.
          L = 0;
          ## Resolve every observation's class once, rather than once per
          ## iteration: the lookup does not depend on i.
          [gYidx, ~] = labelIndices (classes, Y);
          for i = 1:rows (X)
            [~, k] = min (scores(i,:) * this.Cost);
            true_idx = gYidx(i);
            L = L + W(i) * this.Cost(true_idx, k);
          endfor
        case 'classifcost'
          ## What the model's own prediction costs, given the true class
          L = 0;
          ## Resolve every observation's class once, rather than once per
          ## iteration: the lookup does not depend on i.
          [gYidx, ~] = labelIndices (classes, Y);
          for i = 1:rows (X)
            true_idx = gYidx(i);
            pred_idx = find (ismember (classes, label(i)));
            L = L + W(i) * this.Cost(true_idx, pred_idx);
          endfor
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{label} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationGAM} {[@var{label}, @var{score}] =} resubPredict (@var{obj})
    ##
    ## Classify the training data with the generalized additive model it was
    ## fitted on.
    ##
    ## @code{@var{label} = resubPredict (@var{obj})} is @code{predict} applied
    ## to the observations the model was fitted on.
    ##
    ## @seealso{ClassificationGAM, predict}
    ## @end deftypefn
    function [labels, scores] = resubPredict (this)
      used = true (rows (this.X), 1);
      [labels, scores] = predict (this, this.X(used, :));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{m} =} resubMargin (@var{obj})
    ##
    ## Classification margin of a generalized additive model on its training
    ## data.
    ##
    ## @seealso{ClassificationGAM, margin}
    ## @end deftypefn
    function m = resubMargin (this)
      used = true (rows (this.X), 1);
      m = margin (this, this.X(used, :), this.Y(used));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{e} =} resubEdge (@var{obj})
    ##
    ## Classification edge of a generalized additive model on its training
    ## data.
    ##
    ## @seealso{ClassificationGAM, edge}
    ## @end deftypefn
    function e = resubEdge (this)
      used = true (rows (this.X), 1);
      e = edge (this, this.X(used, :), this.Y(used));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {ClassificationGAM} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of a generalized additive model on its training
    ## data.
    ##
    ## @seealso{ClassificationGAM, loss}
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      used = true (rows (this.X), 1);
      L = loss (this, this.X(used, :), this.Y(used), varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationGAM object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## ClassificationGAM object into an Octave binary file, the name of which is
    ## specified in @var{filename}, along with an extra variable, which defines
    ## the type classification object these variables constitute.  Use
    ## @code{loadmodel} in order to load a classification object into Octave's
    ## workspace.
    ##
    ## @seealso{loadmodel, fitcgam, ClassificationGAM}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error ("ClassificationGAM.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error ("ClassificationGAM.savemodel: FNAME must be a character vector.");
      endif
      ## Generate variable for class name
      classdef_name = 'ClassificationGAM';

      ## Create variables from model properties
      X = this.X;
      Y = this.Y;
      NumObservations = this.NumObservations;
      RowsUsed        = this.RowsUsed;
      NumPredictors   = this.NumPredictors;
      PredictorNames  = this.PredictorNames;
      BinEdges        = this.BinEdges;
      ResponseName    = this.ResponseName;
      ClassNames      = this.ClassNames;
      Prior           = this.Prior;
      Cost            = this.Cost;
      ScoreTransform  = this.ScoreTransform;
      Formula         = this.Formula;
      Interactions    = this.Interactions;
      Knots           = this.Knots;
      Order           = this.Order;
      DoF             = this.DoF;
      BaseModel       = this.BaseModel;
      ModelwInt       = this.ModelwInt;
      IntMatrix       = this.IntMatrix;
      Intercept       = this.Intercept;
      W               = this.W;
      CategoricalPredictors  = this.CategoricalPredictors;
      ExpandedPredictorNames = this.ExpandedPredictorNames;
      STfun          = this.STfun;
      FitMethod             = this.FitMethod;
      TreeModel             = this.TreeModel;
      ModelParameters       = this.ModelParameters;
      ReasonForTermination  = this.ReasonForTermination;
      PairDetectionBinEdges = this.PairDetectionBinEdges;

      ## Save classdef name and all model properties as individual variables
      LearningRate    = this.LearningRate;
      NumIterations   = this.NumIterations;

      HyperparameterOptimizationResults = this.HyperparameterOptimizationResults;
      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'NumObservations', ...
            'RowsUsed', 'NumPredictors', 'PredictorNames', 'BinEdges', ...
            'ResponseName', ...
            'ClassNames', 'Prior', 'Cost', 'ScoreTransform', 'Formula', ...
            'Interactions', 'Knots', 'Order', 'DoF', 'BaseModel', ...
            'ModelwInt', 'IntMatrix', 'LearningRate', 'NumIterations', ...
            'Intercept', 'W', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'STfun', 'FitMethod', ...
            'TreeModel', 'ModelParameters', 'ReasonForTermination', ...
            'PairDetectionBinEdges', ...
            'HyperparameterOptimizationResults');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationGAM} {@var{Mdl} =} resume (@var{obj}, @var{numTrees})
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
    ## @seealso{ClassificationGAM, fitcgam, addInteractions}
    ## @end deftypefn
    function Mdl = resume (this, numTrees)

      if (nargin < 2)
        error ("ClassificationGAM.resume: Not enough input arguments.");
      endif
      if (! strcmp (this.FitMethod, 'boostedtrees'))
        error (strcat ("ClassificationGAM.resume: resuming is available", ...
                       " only under 'FitMethod', 'boostedtrees'; a spline", ...
                       " backfit stops at its tolerance and has no budget", ...
                       " to extend."));
      endif
      if (! (isnumeric (numTrees) && isscalar (numTrees) && isreal (numTrees)
             && numTrees > 0 && fix (numTrees) == numTrees))
        error (strcat ("ClassificationGAM.resume: NUMTREES must be a", ...
                       " positive integer scalar."));
      endif

      ## The rows the fit saw, coded as the engine takes them.
      cobs = ! any (isnan (this.X), 2);
      X = this.X(cobs, :);
      [~, ~, gY] = uniqueLabels (this.Y(cobs, :));
      Y = gY(:) - 1;

      Mdl = this;
      MP = this.ModelParameters;
      reason = this.ReasonForTermination;
      ntrees = this.NumTrainedTrees;

      if (isempty (this.TreeModel.Pairs))
        ## No interaction phase ever ran, so the predictor phase is the one
        ## still open.  The engine is handed the prediction reached so far and
        ## returns the increment to add to it.
        f = gamboostpredict (this.BinEdges, this.TreeModel.ShapeValues, X, ...
                             this.Intercept);
        M = gamboosttrain (X, Y, 1, numTrees, ...
                           MP.InitialLearnRateForPredictors, ...
                           MP.MaxNumSplitsPerPredictor, 0, MP.NumPrint, f(:));
        if (M.NumTrees == 0)
          error (strcat ("ClassificationGAM.resume: unable to resume", ...
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
        f = gamboostpredict (this.BinEdges, this.TreeModel.ShapeValues, X, ...
                             this.Intercept, 0, ...
                             this.PairDetectionBinEdges, ...
                             this.TreeModel.PairValues, ...
                             this.TreeModel.Pairs);
        I = gamboostinter (X, Y, f(:), 1, this.TreeModel.Pairs, numTrees, ...
                           MP.InitialLearnRateForInteractions, ...
                           MP.MaxNumSplitsPerInteraction);
        if (I.NumTrees == 0)
          error (strcat ("ClassificationGAM.resume: unable to resume", ...
                         " training because the software was unable to", ...
                         " improve the model fit."));
        endif
        pv = this.TreeModel.PairValues;
        for k = 1:numel (pv)
          pv{k} = pv{k} + I.PairValues{k};
        endfor
        Mdl.TreeModel.PairValues = pv;
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

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationGAM object
      mdl = ClassificationGAM (1, 1);

      ## Copy the saved data into the object.  Iterate over what was
      ## saved rather than over fieldnames (mdl): a private property such
      ## as STfun is written out by savemodel but is not reported by
      ## fieldnames, so comparing the two sets could never match and every
      ## load failed.  Assignment is legal here because this is a method of
      ## the class itself.
      names = fieldnames (data);
      ## The set methods for these read other properties, and one of them
      ## rebuilds Coeffs, so they are assigned once everything else is in
      ## place rather than in the order the file happens to list them.
      late = ismember (names, {'Cost', 'Prior', 'ScoreTransform', ...
                               'ResponseTransform'});
      names = [names(! late); names(late)];
      for i = 1:numel (names)
        try
          mdl.(names{i}) = data.(names{i});
        catch
          error ("ClassificationGAM.load_model: invalid model in '%s'.", filename)
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

  ## Helper functions
  methods(Access = private)

    ## Shared validation for the assessment methods, so each reports under
    ## its own name.
    function [X, Y] = checkXY_ (this, X, Y, caller)
      if (isempty (X))
        error ("ClassificationGAM.%s: X is empty.", caller);
      elseif (this.NumPredictors != columns (X))
        error (strcat ("ClassificationGAM.%s: X must have the same number", ...
                       " of predictors as the trained model."), caller);
      endif
      if (isempty (Y))
        error ("ClassificationGAM.%s: Y is empty.", caller);
      elseif (rows (X) != rows (Y))
        error (strcat ("ClassificationGAM.%s: Y must have the same number", ...
                       " of rows as X."), caller);
      endif
    endfunction

    ## Pull a "Weights" pair out of the optional arguments, defaulting to a
    ## uniform weight, and reject any other name.
    function W = getWeights_ (this, args, n, caller)
      W = ones (n, 1);
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("ClassificationGAM.%s: parameter name must be", ...
                         " a character vector."), caller);
        endif
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
          if (! (isnumeric (W) && isvector (W)))
            error (strcat ("ClassificationGAM.%s: 'Weights' must be a", ...
                           " numeric vector."), caller);
          endif
          if (numel (W) != n)
            error (strcat ("ClassificationGAM.%s: size of 'Weights' must", ...
                           " equal the number of rows in X."), caller);
          endif
        else
          error (strcat ("ClassificationGAM.%s: invalid parameter name in", ...
                         " optional paired arguments."), caller);
        endif
      endfor
    endfunction

    ## Determine interactions from Interactions optional parameter
    ## Fit the model that carries the interaction terms.  The constructor and
    ## Run the interaction phase alone, over a model whose predictor phase is
    ## already fitted.  Shared by addInteractions and by the constructor, so
    ## that asking for interactions after the fact gives the same model as
    ## asking for them up front.
    function this = fitBoostedInteractions (this, X, Y, Interactions, NTI, ...
                                            MSI, LRI, MaxPValue)

      f = gamboostpredict (this.BinEdges, this.TreeModel.ShapeValues, X, ...
                           this.Intercept);

      ## Residuals of the predictor phase, which is what pairs are tested on.
      res = Y - 1 ./ (1 + exp (-f));

      wanted = -1;
      pairs = zeros (0, 2);
      if (ischar (Interactions))
        wanted = Inf;
      elseif (isscalar (Interactions) && ! isempty (Interactions))
        wanted = Interactions;
      elseif (! isempty (Interactions))
        pairs = interactionPairs (logical (Interactions));
      endif

      if (wanted > 0 && columns (X) > 1)
        S = gamboostpairs (X, res);
        pval = 1 - fcdf (S.F, S.DF1, S.DF2);
        pval(S.DF1 <= 0) = 1;
        [pval, ord] = sort (pval);
        ranked = S.Pairs(ord, :);
        ranked = ranked(pval <= MaxPValue, :);
        if (isfinite (wanted) && rows (ranked) > wanted)
          ranked = ranked(1:wanted, :);
        endif
        pairs = ranked;
        this.PairDetectionBinEdges = S.BinEdges(:);
        if (isempty (pairs))
          warning (strcat ("ClassificationGAM: model does not include", ...
                           " interaction terms because all interaction", ...
                           " terms have p-values greater than the", ...
                           " 'MaxPValue' value, or the software was unable", ...
                           " to improve the model fit."));
        endif
      endif

      reason = this.ReasonForTermination;
      ntrees = this.NumTrainedTrees;
      if (! isempty (pairs))
        I = gamboostinter (X, Y, f, 1, pairs, NTI, LRI, MSI);
        this.Intercept = this.Intercept + I.Intercept;
        this.PairDetectionBinEdges = I.PairBinEdges(:);
        this.TreeModel.PairValues = I.PairValues;
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
      MP = this.ModelParameters;
      if (ischar (Interactions))
        MP.Interactions = Interactions;
      elseif (isempty (Interactions))
        MP.Interactions = 0;
      else
        MP.Interactions = Interactions;
      endif
      this.ModelParameters = MP;

    endfunction

    ## Drive the boosted-tree engine: the predictor phase, then a search for
    ## interactions worth adding, then the interaction phase over whichever
    ## pairs survived.  The two phases share a running fit, so the second
    ## continues from the prediction the first left rather than starting over,
    ## which is how MATLAB's own trace behaves.
    function this = fitBoosted (this, X, Y, Interactions, NTP, NTI, MSP, ...
                                MSI, LRP, LRI, MaxPValue, Verb, NPrint)

      ## The predictor phase.
      M = gamboosttrain (X, Y, 1, NTP, LRP, MSP, Verb, NPrint);
      f = gamboostpredict (M.BinEdges, M.ShapeValues, X, M.Intercept);

      this.BinEdges  = M.BinEdges(:);   ## a column cell, as MATLAB reports it
      this.Intercept = M.Intercept;
      reason = struct ('PredictorTrees', M.ReasonForTermination, ...
                       'InteractionTrees', '');
      ntrees = struct ('PredictorTrees', M.NumTrees, ...
                       'InteractionTrees', 0);
      pairs = zeros (0, 2);
      pairValues = {};
      pairShift = 0;

      ## How many pairs were asked for.  A count or 'all' means search; an
      ## explicit Nx2 matrix names the pairs outright and skips the test, as
      ## MATLAB takes a matrix at its word.
      wanted = -1;
      if (ischar (Interactions))
        wanted = Inf;
      elseif (isscalar (Interactions) && ! isempty (Interactions))
        wanted = Interactions;
      elseif (! isempty (Interactions))
        ## A matrix names the terms outright.  This class has always taken
        ## that as a term matrix, one row per term and one column per
        ## predictor, so it is converted to the pairs the engine works in
        ## rather than being mistaken for a two-column list of indices.
        pairs = interactionPairs (logical (Interactions));
      endif

      if (wanted > 0 && columns (X) > 1)
        S = gamboostpairs (X, M.Residuals);
        ## The F ratio becomes a probability through the package's own fcdf,
        ## which is verified against MATLAB; the engine deliberately does not
        ## carry a second incomplete beta of its own.
        pval = 1 - fcdf (S.F, S.DF1, S.DF2);
        pval(S.DF1 <= 0) = 1;
        [pval, ord] = sort (pval);
        ranked = S.Pairs(ord, :);
        keep = pval <= MaxPValue;
        ranked = ranked(keep, :);
        if (isfinite (wanted) && rows (ranked) > wanted)
          ranked = ranked(1:wanted, :);
        endif
        pairs = ranked;
        this.PairDetectionBinEdges = S.BinEdges(:);
        if (isempty (pairs))
          warning (strcat ("ClassificationGAM: model does not include", ...
                           " interaction terms because all interaction", ...
                           " terms have p-values greater than the", ...
                           " 'MaxPValue' value, or the software was unable", ...
                           " to improve the model fit."));
        endif
      endif

      if (! isempty (pairs))
        I = gamboostinter (X, Y, f, 1, pairs, NTI, LRI, MSI);
        this.Intercept = this.Intercept + I.Intercept;
        pairShift = I.Intercept;
        this.PairDetectionBinEdges = I.PairBinEdges(:);
        pairValues = I.PairValues;
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
                               'Pairs', pairs, ...
                               'PairIntercept', pairShift);

      ## MATLAB reports the request rather than the selection here, and the
      ## selection through the Interactions property.  NumPrint and
      ## VerbosityLevel are reported at their defaults: this engine keeps no
      ## printed trace, so they are stated rather than accepted and ignored.
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
        'Type', 'classification');

    endfunction

    ## addInteractions both arrive here with IntMatrix already decided and the
    ## predictors and response prepared as the fit wants them, so the two
    ## cannot drift: a model given its interactions after the fact is the
    ## model it would have been had they been asked for at the outset.
    function this = fitModelwInt (this, X, Y, Inter, Knots, Order, DoF, ...
                                  LearningRate, NumIterations)

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
      [iter, param, res, RSS, intercept] = this.fitGAM (X, Y, Inter, Knots, ...
                                                        Order, LearningRate, ...
                                                        NumIterations);
      this.ModelwInt.Intercept  = intercept;
      this.ModelwInt.Parameters = param;
      this.ModelwInt.Iterations = iter;
      this.ModelwInt.Residuals  = res;
      this.ModelwInt.RSS        = RSS;
    endfunction

    function intMat = parseInteractions (this)
      if (islogical (this.Interactions))
        ## Check that interaction matrix corresponds to predictors
        if (numel (this.PredictorNames) != columns (this.Interactions))
          error (strcat ("ClassificationGAM: columns in 'Interactions'", ...
                         " matrix must equal to the number of predictors."));
        endif
        intMat = this.Interactions;
      elseif (isnumeric (this.Interactions))
        ## Need to measure the effect of all interactions to keep the best
        ## performing. Just check that the given number is not higher than
        ## p*(p-1)/2, where p is the number of predictors.
        p = this.NumPredictors;
        if (this.Interactions > p * (p - 1) / 2)
          error (strcat ("ClassificationGAM: number of interaction terms", ...
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
        error ("ClassificationGAM: invalid syntax in 'Formula'.");
      endif
      ## Split formula and keep predictor terms
      formulaParts = strsplit (this.Formula, '~');
      ## Check there is some string after '~'
      if (numel (formulaParts) < 2)
        error ("ClassificationGAM: no predictor terms in 'Formula'.");
      endif
      predictorString = strtrim (formulaParts{2});
      if (isempty (predictorString))
        error ("ClassificationGAM: no predictor terms in 'Formula'.");
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
            error (strcat ("ClassificationGAM: some predictors", ...
                           " have not been identified."));
          endif
          ## Append to interactions matrix
          intMat = [intMat; iterms];
        endif
      endfor
      ## Check that all terms have been identified
      if (! all (sum (intMat, 2) > 0))
        error ("ClassificationGAM: some terms have not been identified.");
      endif
    endfunction

    ## Fit the model
    function [iter, param, res, RSS, intercept] = fitGAM (this, X, Y, Inter, ...
                                    Knots, Order, learning_rate, num_iterations)
      ## The fit is performed by the shared spline engine, which builds and
      ## factorises each predictor's design once and reduces a boosting round
      ## to two products against the factors.
      Mdl = gamtrain (X, Y, Knots, Order, 1, Inter, learning_rate, ...
                      num_iterations);
      iter      = Mdl.Iterations;
      param     = Mdl.Parameters;
      res       = Mdl.Residuals;
      RSS       = Mdl.RSS;
      intercept = Mdl.Intercept;
    endfunction

    ## Set cost

  endmethods

endclassdef

## Helper function
function scores = predict_val (params, XC, intercept)
  ## The shared prediction engine evaluates every additive term and sums
  ## them.  That sum is the log-odds of the second class, so the raw score of
  ## the first is its negative and the pair sums to zero, as MATLAB's does.
  f = gampredict (params, XC, intercept, 0);
  scores = [-f, f];
endfunction

%!demo
%! ## Train a GAM classifier for binary classification
%! ## using specific data and plot the decision boundaries.
%!
%! ## Define specific data
%! X = [1, 2; 2, 3; 3, 3; 4, 5; 5, 5; ...
%!     6, 7; 7, 8; 8, 8; 9, 9; 10, 10];
%! Y = [0; 0; 0; 0; 0; ...
%!     1; 1; 1; 1; 1];
%!
%! ## Train the GAM model
%! obj = fitcgam (X, Y, 'Interactions', 'all')
%!
%! ## Create a grid of values for prediction
%! x1 = [min(X(:,1)):0.1:max(X(:,1))];
%! x2 = [min(X(:,2)):0.1:max(X(:,2))];
%! [x1G, x2G] = meshgrid (x1, x2);
%! XGrid = [x1G(:), x2G(:)];
%! [labels, score] = predict (obj, XGrid);

## Test constructor
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [0; 0; 1; 1];
%! PredictorNames = {'Feature1', 'Feature2', 'Feature3'};
%! a = ClassificationGAM (x, y, 'FitMethod', 'splines', ...
%!                       'PredictorNames', PredictorNames);
%! assert_equal (class (a), "ClassificationGAM");
%! assert_equal ({a.X, a.Y, a.NumObservations}, {x, y, 4})
%! assert_equal ({a.NumPredictors, a.ResponseName}, {3, 'Y'})
%! assert_equal (a.ClassNames, [0; 1])
%! assert_equal (a.PredictorNames, PredictorNames)
%! assert_equal (a.BaseModel.Intercept, 0)
%!test
%! load fisheriris
%! inds = strcmp (species,'versicolor') | strcmp (species,'virginica');
%! X = meas(inds, :);
%! Y = species(inds, :)';
%! Y = strcmp (Y, 'virginica')';
%! a = ClassificationGAM (X, Y, 'FitMethod', 'splines', ...
%!                       'Formula', 'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3');
%! assert_equal (class (a), "ClassificationGAM");
%! assert_equal ({a.X, a.Y, a.NumObservations}, {X, Y, 100})
%! assert_equal ({a.NumPredictors, a.ResponseName}, {4, 'Y'})
%! assert_equal (a.ClassNames, logical ([0; 1]))
%! assert_equal (a.Formula, 'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3')
%! assert_equal (a.PredictorNames, {'x1', 'x2', 'x3', 'x4'})
%! assert_equal (a.ModelwInt.Intercept, 0)
%!test
%! X = [2, 3, 5; 4, 6, 8; 1, 2, 3; 7, 8, 9; 5, 4, 3];
%! Y = [0; 1; 0; 1; 1];
%! a = ClassificationGAM (X, Y, 'FitMethod', 'splines', ...
%!                       'Knots', [4, 4, 4], 'Order', [3, 3, 3]);
%! assert_equal (class (a), "ClassificationGAM");
%! assert_equal ({a.X, a.Y, a.NumObservations}, {X, Y, 5})
%! assert_equal ({a.NumPredictors, a.ResponseName}, {3, 'Y'})
%! assert_equal (a.ClassNames, [0; 1])
%! assert_equal (a.PredictorNames, {'x1', 'x2', 'x3'})
%! assert_equal (a.Knots, [4, 4, 4])
%! assert_equal (a.Order, [3, 3, 3])
%! assert_equal (a.DoF, [7, 7, 7])
%! assert_equal (a.BaseModel.Intercept, 0.4055, 1e-1)

## Test Prior calculation
%!test
%! ## Test uniform prior
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = [0; 0; 1; 1];
%! a = ClassificationGAM (x, y, 'Prior', 'uniform');
%! assert_equal (a.Prior, [0.5, 0.5], 1e-6);
%!test
%! ## Test empirical prior
%! x = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10];
%! y = [0; 0; 0; 1; 1];
%! a = ClassificationGAM (x, y, 'Prior', 'empirical');
%! assert_equal (a.Prior, [0.6, 0.4], 1e-6);
%!test
%! ## Test numeric prior
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = [0; 0; 1; 1];
%! a = ClassificationGAM (x, y, 'Prior', [0.7, 0.3]);
%! assert_equal (a.Prior, [0.7, 0.3], 1e-6);
%!test
%! ## Test default prior (empirical)
%! x = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10; 11, 12];
%! y = [0; 0; 0; 1; 1; 1];
%! a = ClassificationGAM (x, y);
%! assert_equal (a.Prior, [0.5, 0.5], 1e-6);
%!test
%! ## Test prior normalization
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = [0; 0; 1; 1];
%! a = ClassificationGAM (x, y, 'Prior', [2, 1]);
%! assert_equal (a.Prior, [2/3, 1/3], 1e-6);

## Test input validation for Prior
## Interactions reports the two-way terms the fitted model carries, as
## predictor index pairs.  R2024a's GAM with 'Interactions', 'all' over three
## predictors returns [1 2; 1 3; 2 3], which is what this matches.
%!test
%! k = (1:60)';
%! X = [mod(k*7,11)-5, mod(k*3,11)-5, mod(k*5,11)-5];
%! y = double (X(:,1).*X(:,2) > 0) + 1;
%! Mdl = fitcgam (X, y, "Interactions", "all");
%! assert_equal (Mdl.Interactions, [1, 2; 1, 3; 2, 3]);

## No interactions is an empty list of pairs, keeping its two columns, and
## not an empty matrix of no width.
%!test
%! k = (1:60)';
%! X = [mod(k*7,11)-5, mod(k*3,11)-5, mod(k*5,11)-5];
%! y = double (X(:,1).*X(:,2) > 0) + 1;
%! Mdl = fitcgam (X, y);
%! assert_equal (size (Mdl.Interactions), [0, 2]);
%! assert_equal (class (Mdl.Interactions), "double");

## The same kind of value follows whatever form the request took.
%!test
%! k = (1:60)';
%! X = [mod(k*7,11)-5, mod(k*3,11)-5, mod(k*5,11)-5];
%! y = double (X(:,1).*X(:,2) > 0) + 1;
%! Mc = fitcgam (X, y, "Interactions", 2);
%! Ml = fitcgam (X, y, "Interactions", logical ([1, 1, 0; 0, 1, 1]));
%! assert_equal (Mc.Interactions, [1, 2; 1, 3]);
%! assert_equal (Ml.Interactions, [1, 2; 2, 3]);

## A formula names its main effects as terms of the model, and a main effect
## is not an interaction: the term matrix holds all three, Interactions the
## one two-way term among them.
%!test
%! k = (1:60)';
%! X = [mod(k*7,11)-5, mod(k*3,11)-5, mod(k*5,11)-5];
%! y = double (X(:,1).*X(:,2) > 0) + 1;
%! Mdl = fitcgam (X, y, "FitMethod", "splines", ...
%!                "Formula", "Y ~ x1 + x2 + x1:x2");
%! assert_equal (Mdl.Interactions, [1, 2]);
%! assert_equal (rows (Mdl.IntMatrix), 3);

## The compact model carries the same pairs.
%!test
%! k = (1:60)';
%! X = [mod(k*7,11)-5, mod(k*3,11)-5, mod(k*5,11)-5];
%! y = double (X(:,1).*X(:,2) > 0) + 1;
%! Mdl = fitcgam (X, y, "Interactions", "all");
%! assert_equal (compact (Mdl).Interactions, Mdl.Interactions);

## A response naming its classes in the rows of a character matrix is one
## of the documented types and MATLAB accepts it on every classifier.  The
## whole surface below was broken and untested, which is why it stayed so.
%!test
%! load fisheriris
%! bch = ! strcmp (species, "setosa");
%! Xch = meas(bch,:); Ycell = species(bch); Ych = char (Ycell);
%! rand ("state", 1); randn ("state", 1); Mc = fitcgam (Xch, Ych);
%! rand ("state", 1); randn ("state", 1); Ms = fitcgam (Xch, Ycell);
%! assert_equal (size (Mc.ClassNames), [2, 10]);
%! assert_equal (cellstr (Mc.ClassNames), Ms.ClassNames);

## predict returns whole names, not their first letters.
%!test
%! load fisheriris
%! bch = ! strcmp (species, "setosa");
%! Xch = meas(bch,:); Ycell = species(bch); Ych = char (Ycell);
%! rand ("state", 1); randn ("state", 1); Mc = fitcgam (Xch, Ych);
%! rand ("state", 1); randn ("state", 1); Ms = fitcgam (Xch, Ycell);
%! pch = predict (Mc, Xch);
%! assert_equal (columns (pch), 10);
%! assert_equal (cellstr (pch), predict (Ms, Xch));

## loss, margin and edge read a character response as the same response.
%!test
%! load fisheriris
%! bch = ! strcmp (species, "setosa");
%! Xch = meas(bch,:); Ycell = species(bch); Ych = char (Ycell);
%! rand ("state", 1); randn ("state", 1); Mc = fitcgam (Xch, Ych);
%! rand ("state", 1); randn ("state", 1); Ms = fitcgam (Xch, Ycell);
%! assert_equal (loss (Mc, Xch, Ych), loss (Ms, Xch, Ycell), 1e-12);
%! assert_equal (margin (Mc, Xch, Ych), margin (Ms, Xch, Ycell), 1e-12);
%! assert_equal (edge (Mc, Xch, Ych), edge (Ms, Xch, Ycell), 1e-12);

## A character matrix pads its rows out to the longest name, and the padding
## is part of the name: R2024a reports ClassNames of ['ab  '; 'abcd'].
%!test
%! Xpad = [1 2; 3 4; 1.1 2.1; 3.1 4.1; 1.2 2.2; 3.2 4.2];
%! Ypad = char ({"ab", "abcd", "ab", "abcd", "ab", "abcd"});
%! rand ("state", 1); randn ("state", 1);
%! Mp = fitcgam (Xpad, Ypad);
%! assert_equal (size (Mp.ClassNames), [2, 4]);
%! assert_equal (Mp.ClassNames(1,:), "ab  ");

## A row dropped for a missing predictor is the only case that exercises
## indexing the response by row rather than by element.
%!test
%! load fisheriris
%! bch = ! strcmp (species, "setosa");
%! Xch = meas(bch,:); Ycell = species(bch); Ych = char (Ycell);
%! Xmiss = Xch; Xmiss(3,2) = NaN;
%! rand ("state", 1); randn ("state", 1); Md = fitcgam (Xmiss, Ych);
%! rand ("state", 1); randn ("state", 1); Ms = fitcgam (Xmiss, Ycell);
%! assert_equal (size (Md.ClassNames), [2, 10]);
%! assert_equal (cellstr (Md.ClassNames), Ms.ClassNames);

## ClassNames may itself be given as a character matrix, which selects the
## classes by whole name: ismember between two character matrices compares
## them character by character and would select by letter.
%!test
%! load fisheriris
%! rand ("state", 1); randn ("state", 1);
%! Mf = fitcgam (meas, char (species), ...
%!               "ClassNames", char ({"versicolor", "virginica"}));
%! assert_equal (rows (Mf.ClassNames), 2);
%! assert_equal (cellstr (Mf.ClassNames), {"versicolor"; "virginica"});

## A model fitted from a character response comes back off disk unchanged.
%!test
%! load fisheriris
%! bch = ! strcmp (species, "setosa");
%! Xch = meas(bch,:); Ycell = species(bch); Ych = char (Ycell);
%! rand ("state", 1); randn ("state", 1); Mc = fitcgam (Xch, Ych);
%! fname = tempname ();
%! savemodel (Mc, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (M2.ClassNames, Mc.ClassNames);
%! assert_equal (predict (M2, Xch), predict (Mc, Xch));

## crossval carries a character response through cvpartition and back.
%!test
%! load fisheriris
%! bch = ! strcmp (species, "setosa");
%! Xch = meas(bch,:); Ycell = species(bch); Ych = char (Ycell);
%! rand ("state", 1); randn ("state", 1); Mc = fitcgam (Xch, Ych);
%! rand ("state", 1); randn ("state", 1); Ms = fitcgam (Xch, Ycell);
%! rand ("state", 2); cvc = crossval (Mc, "KFold", 3);
%! rand ("state", 2); cvs = crossval (Ms, "KFold", 3);
%! assert_equal (cellstr (kfoldPredict (cvc)), kfoldPredict (cvs));

## addInteractions fits the interaction terms onto a model that already has
## its univariate ones.  The result is the model that would have been fitted
## had the terms been asked for at the outset: both go through one private
## method, so the two cannot drift apart.
%!test
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! Xai = meas(bai,2:4); Yai = species(bai);
%! Aai = addInteractions (fitcgam (Xai, Yai), "all");
%! Bai = fitcgam (Xai, Yai, "Interactions", "all");
%! assert_equal (Aai.Interactions, Bai.Interactions);
%! assert_equal (Aai.ModelwInt, Bai.ModelwInt);
%! assert_equal (predict (Aai, Xai), predict (Bai, Xai));

## The univariate fit is left alone, which is what MATLAB leaves alone too.
%!test
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! Xai = meas(bai,2:4); Yai = species(bai);
%! Cai = fitcgam (Xai, Yai);
%! Aai = addInteractions (Cai, "all");
%! assert_equal (predict (Aai, Xai, "IncludeInteractions", false), ...
%!               predict (Cai, Xai));

## A count and a logical matrix name terms as the constructor's option does.
%!test
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! Xai = meas(bai,2:4); Yai = species(bai);
%! Aai = addInteractions (fitcgam (Xai, Yai, 'FitMethod', 'splines'), 2);
%! assert_equal (Aai.Interactions, [1, 2; 1, 3]);
%! Lai = addInteractions (fitcgam (Xai, Yai, 'FitMethod', 'splines'), ...
%!                        logical ([1 1 0; 0 1 1]));
%! assert_equal (Lai.Interactions, [1, 2; 2, 3]);

## The boosted-tree engine names the same terms from the same specifications,
## but a count takes the pairs the interaction search ranked highest rather
## than the first in index order, so the two engines select differently.
##
## Only the top of the ranking is pinned, and deliberately.  On this fixture
## one pair carries an interaction and the other two do not: (2,3) scores
## F 3.25 at p 5.5e-06, while (1,2) and (1,3) score F 0.5599 and F 0.4331,
## both below one, meaning their cells explain less than the noise within
## them.  Their p values sit at 0.974 and 0.998, so which of the two the
## search ranks second is decided by 2% of the scale at the very top of it
## and turns over on any change that moves the predictor-phase residuals at
## all.  Pinning that order pins nothing about the search, and it broke when
## the tree learner gained its minimum leaf size.  What is worth holding is
## that the pair with signal comes first, that a count does not simply take
## the first pair in index order, and that 'all' names every pair.
%!test
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! Xai = meas(bai,2:4); Yai = species(bai);
%! Aai = addInteractions (fitcgam (Xai, Yai), 2);
%! assert_equal (rows (Aai.Interactions), 2);
%! assert_equal (Aai.Interactions(1,:), [2, 3]);
%! Lai = addInteractions (fitcgam (Xai, Yai), logical ([1 1 0; 0 1 1]));
%! assert_equal (Lai.Interactions, [1, 2; 2, 3]);
%! All = addInteractions (fitcgam (Xai, Yai), "all");
%! assert_equal (All.Interactions(1,:), [2, 3]);
%! assert_equal (sortrows (All.Interactions), [1, 2; 1, 3; 2, 3]);

## A model that already carries interaction terms is not extended, and a
## model fitted from a formula names every term it has, interactions among
## them, so it is refused for the same reason.  R2024a refuses both.
## resume continues the phase that ran last, and a model with no interactions
## has only the predictor phase open.  Resuming reproduces the model a single
## fit of the combined budget would have produced, which is the oracle every
## test here uses.
%!test
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! A = fitcgam (X, Y, 'NumTreesPerPredictor', 5);
%! B = resume (A, 10);
%! C = fitcgam (X, Y, 'NumTreesPerPredictor', 15);
%! assert_equal (B.ModelParameters.NumTreesPerPredictor, 15);
%! [~, sB] = predict (B, X);
%! [~, sC] = predict (C, X);
%! assert_equal (sB, sC, 1e-12);

%!test
%! ## A model carrying interactions gains interaction trees, and its predictor
%! ## shape functions are left where they were.
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! A = fitcgam (X, Y, 'NumTreesPerPredictor', 5, 'Interactions', 3, ...
%!              'NumTreesPerInteraction', 4);
%! B = resume (A, 10);
%! assert_equal (B.ModelParameters.NumTreesPerPredictor, 5);
%! assert_equal (B.ModelParameters.NumTreesPerInteraction, 14);
%! assert_equal (B.TreeModel.ShapeValues, A.TreeModel.ShapeValues);
%! C = fitcgam (X, Y, 'NumTreesPerPredictor', 5, 'Interactions', 3, ...
%!              'NumTreesPerInteraction', 14);
%! [~, sB] = predict (B, X);
%! [~, sC] = predict (C, X);
%! assert_equal (sB, sC, 1e-12);

%!test
%! ## The selected pairs survive, and resuming twice accumulates.
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! A = fitcgam (X, Y, 'NumTreesPerPredictor', 5, 'Interactions', 3, ...
%!              'NumTreesPerInteraction', 4);
%! B = resume (resume (A, 10), 6);
%! assert_equal (B.Interactions, A.Interactions);
%! assert_equal (B.ModelParameters.NumTreesPerInteraction, 20);
%! assert_equal (B.ModelParameters.NumTreesPerPredictor, 5);

%!test
%! ## The model handed in is not modified.
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! A = fitcgam (X, Y, 'NumTreesPerPredictor', 5);
%! B = resume (A, 10);
%! assert_equal (A.ModelParameters.NumTreesPerPredictor, 5);
%! assert_equal (B.ModelParameters.NumTreesPerPredictor, 15);

%!error<ClassificationGAM.resume: Not enough input arguments.> ...
%! load fisheriris; ...
%! resume (fitcgam (meas(51:150,:), species(51:150), 'NumTreesPerPredictor', 5))
%!error<ClassificationGAM.resume: resuming is available only under 'FitMethod', 'boostedtrees'; a spline backfit stops at its tolerance and has no budget to extend.> ...
%! load fisheriris; ...
%! resume (fitcgam (meas(51:150,:), species(51:150), 'FitMethod', 'splines'), 5)
%!error<ClassificationGAM.resume: NUMTREES must be a positive integer scalar.> ...
%! load fisheriris; ...
%! resume (fitcgam (meas(51:150,:), species(51:150), 'NumTreesPerPredictor', 5), 0)
%!error<ClassificationGAM.resume: NUMTREES must be a positive integer scalar.> ...
%! load fisheriris; ...
%! resume (fitcgam (meas(51:150,:), species(51:150), 'NumTreesPerPredictor', 5), 2.5)
%!error<ClassificationGAM.resume: NUMTREES must be a positive integer scalar.> ...
%! load fisheriris; ...
%! resume (fitcgam (meas(51:150,:), species(51:150), 'NumTreesPerPredictor', 5), [1, 2])

%!error<ClassificationGAM.addInteractions: adding interaction terms to a model that already includes them is not supported.> ...
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! Mai = fitcgam (meas(bai,2:4), species(bai), "Interactions", 2);
%! addInteractions (Mai, "all")
%!error<ClassificationGAM.addInteractions: adding interaction terms to a model that already includes them is not supported.> ...
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! addInteractions (fitcgam (meas(bai,2:4), species(bai), ...
%!                  "FitMethod", "splines", ...
%!                  "Formula", "Y ~ x1 + x2 + x1:x2"), "all")
%!error<ClassificationGAM.addInteractions: invalid 'Interactions' parameter.> ...
%! load fisheriris
%! bai = ! strcmp (species, "setosa");
%! addInteractions (fitcgam (meas(bai,2:4), species(bai)), {1})

%!error<ClassificationGAM: 'Prior' must be a 2-element vector.> ...
%! ClassificationGAM (ones (4,2), ones (4,1), 'Prior', [1])
%!error<ClassificationGAM: 'Prior' must be a 2-element vector.> ...
%! ClassificationGAM (ones (4,2), ones (4,1), 'Prior', [1, 2, 3])
%!error<ClassificationGAM: 'Prior' must be a numeric vector or a string.> ...
%! ClassificationGAM (ones (4,2), ones (4,1), 'Prior', {1, 2})
%!error<ClassificationGAM: 'Prior' must be> ...
%! ClassificationGAM (ones (4,2), ones (4,1), 'Prior', 'invalid')

## Test input validation for constructor
%!error<ClassificationGAM: too few input arguments.> ClassificationGAM ()
%!error<ClassificationGAM: too few input arguments.> ...
%! ClassificationGAM (ones (4, 1))
%!error<ClassificationGAM: number of rows in X and Y must be equal.> ...
%! ClassificationGAM (ones (4,2), ones (1,4))
%!error<ClassificationGAM: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), 'PredictorNames', ['A'])
%!error<ClassificationGAM: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), 'PredictorNames', 'A')
%!error<ClassificationGAM: 'PredictorNames' must equal the number of columns in X.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), 'PredictorNames', {'A', 'B', 'C'})
%!error<ClassificationGAM: 'ResponseName' must be a character vector.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), 'ResponseName', {'Y'})
%!error<ClassificationGAM: 'ResponseName' must be a character vector.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), 'ResponseName', 1)
%!error<ClassificationGAM: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationGAM (ones (10,2), ones (10,1), 'ClassNames', @(x)x)
%!error<ClassificationGAM: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationGAM (ones (10,2), ones (10,1), 'ClassNames', {1})
%!error<ClassificationGAM: not all 'ClassNames' are present in Y.> ...
%! ClassificationGAM (ones (10,2), ones (10,1), 'ClassNames', [1, 2])
%!error<ClassificationGAM: not all 'ClassNames' are present in Y.> ...
%! ClassificationGAM (ones (5,2), ['a';'b';'a';'a';'b'], 'ClassNames', ['a';'c'])
%!error<ClassificationGAM: not all 'ClassNames' are present in Y.> ...
%! ClassificationGAM (ones (5,2), {'a';'b';'a';'a';'b'}, 'ClassNames', {'a','c'})
%!error<ClassificationGAM: not all 'ClassNames' are present in Y.> ...
%! ClassificationGAM (ones (10,2), logical (ones (10,1)), 'ClassNames', [true, false])
%!error<ClassificationGAM: 'Cost' must be a numeric square matrix.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), 'Cost', [1, 2])
%!error<ClassificationGAM: the number of rows and columns in 'Cost' must correspond to selected classes in Y.> ...
%! Mdl = fitcgam ([1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1], [0; 0; 1; 1]);
%! Mdl.Cost = 1:4;
%!error<ClassificationGAM: 'Cost' must be a numeric square matrix.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), 'Cost', 'string')
%!error<ClassificationGAM: 'Cost' must be a numeric square matrix.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), 'Cost', {eye(2)})

## Test predict method
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10];
%! y = [1; 0; 1; 0; 1];
%! a = ClassificationGAM (x, y, 'FitMethod', 'splines', ...
%!                       'interactions', 'all');
%! l = [1; 0; 1; 0; 1];
%! s = [0.0334, 0.9666; 0.9648, 0.0352; 0.0334, 0.9666; ...
%!      0.9648, 0.0352; 0.0334, 0.9666];
%! [labels, scores] = predict (a, x);
%! assert_equal (class (a), "ClassificationGAM");
%! assert_equal ({a.X, a.Y, a.NumObservations}, {x, y, 5})
%! assert_equal ({a.NumPredictors, a.ResponseName}, {2, 'Y'})
%! assert_equal (a.ClassNames, [0; 1])
%! assert_equal (a.PredictorNames, {'x1', 'x2'})
%! assert_equal (a.ModelwInt.Intercept, 0.4055, 1e-1)
%! assert_equal (labels, l)
%! assert_equal (scores, s, 1e-1)
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [0; 0; 1; 1];
%! interactions = [false, true, false; true, false, true; false, true, false];
%! a = fitcgam (x, y, 'FitMethod', 'splines', ...
%!               'learningrate', 0.2, 'interactions', interactions);
%! [label, score] = predict (a, x, 'includeinteractions', true);
%! l = [0; 0; 1; 1];
%! s = [0.9725, 0.0275; 0.9895, 0.0105; 0.0070, 0.9930; 0.0238, 0.9762];
%! assert_equal (class (a), "ClassificationGAM");
%! assert_equal ({a.X, a.Y, a.NumObservations}, {x, y, 4})
%! assert_equal ({a.NumPredictors, a.ResponseName}, {3, 'Y'})
%! assert_equal (a.ClassNames, [0; 1])
%! assert_equal (a.PredictorNames, {'x1', 'x2', 'x3'})
%! assert_equal (a.ModelwInt.Intercept, 0)
%! assert_equal (label, l)
%! assert_equal (score, s, 1e-1)

## Test input validation for predict method
%!error<ClassificationGAM.predict: too few input arguments.> ...
%! predict (ClassificationGAM (ones (4,2), ones (4,1)))
%!error<ClassificationGAM.predict: XC is empty.> ...
%! predict (ClassificationGAM (ones (4,2), ones (4,1)), [])
%!error<ClassificationGAM.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (ClassificationGAM (ones (4,2), ones (4,1)), 1)

## Test crossval method
## A numeric response is coded 0/1 for the fitter whatever its own labels
## are.  A response of 1 and 2 used to reach the fitter as it was given: the
## seed is log (mean (Y) / (1 - mean (Y))), whose argument is negative at a
## mean of 1.5, so the intercept came back NaN and took every score with it,
## and predict answered class 1 for every row.
%!test
%! Xn = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! Mn = fitcgam (Xn, [1; 1; 2; 2], 'FitMethod', 'splines');
%! assert_equal (Mn.BaseModel.Intercept, 0, 1e-12);
%! [label, score] = predict (Mn, Xn);
%! assert_equal (label, [1; 1; 2; 2]);
%! assert_equal (all (isfinite (score(:))), true);

## Labels nowhere near 0 and 1 fit the same way, the coding being the class
## index rather than the label.
%!test
%! Xn = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! Mn = fitcgam (Xn, [5; 5; 9; 9], 'FitMethod', 'splines');
%! assert_equal (Mn.BaseModel.Intercept, 0, 1e-12);
%! assert_equal (predict (Mn, Xn), [5; 5; 9; 9]);

## A response already coded 0 and 1 is unchanged by that coding.
%!test
%! Xn = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! M0 = fitcgam (Xn, [0; 0; 1; 1], 'FitMethod', 'splines');
%! assert_equal (M0.BaseModel.Intercept, 0, 1e-12);
%! assert_equal (predict (M0, Xn), [0; 0; 1; 1]);

## The boosted engine fits its predictors by local scoring: the working
## response and the weights are formed once per cycle and held fixed across
## it, so the first predictor sees the residual it would see fitted alone and
## each one after sees what the earlier ones took out.  Measured against
## R2024a over two predictors at one tree each, on the two overlapping species
## with every third label flipped.
%!test
%! load fisheriris
%! X = meas(51:150,1:2);
%! Y = species(51:150);
%! f = 1:3:100;
%! t = Y(f);
%! t(strcmp (t, 'versicolor')) = {'ZZ'};
%! t(strcmp (t, 'virginica')) = {'versicolor'};
%! t(strcmp (t, 'ZZ')) = {'virginica'};
%! Y(f) = t;
%! o = {'NumTreesPerPredictor', 1, 'MaxNumSplitsPerPredictor', 1, ...
%!      'Interactions', 0};
%! g1 = (4.8:0.1:8.0)';
%! g2 = (1.9:0.1:3.9)';
%! M = fitcgam (X, Y, o{:});
%! M.ScoreTransform = 'none';
%! [~, s0] = predict (M, [4.8, 1.9]);
%! [~, sA] = predict (M, [g1, repmat(1.9, numel (g1), 1)]);
%! [~, sB] = predict (M, [repmat(4.8, numel (g2), 1), g2]);
%! ## the first predictor is what it is fitted alone, the second is not
%! assert_equal (max (sA(:,2)) - s0(1,2), 1.19047619047619, 1e-12);
%! assert_equal (max (sB(:,2)) - s0(1,2), 1.04232804232804, 1e-12);
%! A1 = fitcgam (X(:,1), Y, o{:});
%! A1.ScoreTransform = 'none';
%! [~, a0] = predict (A1, 4.8);
%! [~, aA] = predict (A1, g1);
%! assert_equal (max (aA(:,2)) - a0(1,2), 1.19047619047619, 1e-12);
%! A2 = fitcgam (X(:,2), Y, o{:});
%! A2.ScoreTransform = 'none';
%! [~, b0] = predict (A2, 1.9);
%! [~, bB] = predict (A2, g2);
%! assert_equal (max (bB(:,2)) - b0(1,2), 1.33333333333333, 1e-12);

## Swapping the predictors swaps which of them is fitted first, so both shape
## functions change.  R2024a again.
%!test
%! load fisheriris
%! X = meas(51:150,1:2);
%! Y = species(51:150);
%! f = 1:3:100;
%! t = Y(f);
%! t(strcmp (t, 'versicolor')) = {'ZZ'};
%! t(strcmp (t, 'virginica')) = {'versicolor'};
%! t(strcmp (t, 'ZZ')) = {'virginica'};
%! Y(f) = t;
%! o = {'NumTreesPerPredictor', 1, 'MaxNumSplitsPerPredictor', 1, ...
%!      'Interactions', 0};
%! S = fitcgam (X(:,[2, 1]), Y, o{:});
%! S.ScoreTransform = 'none';
%! g1 = (4.8:0.1:8.0)';
%! g2 = (1.9:0.1:3.9)';
%! [~, t0] = predict (S, [1.9, 4.8]);
%! [~, tA] = predict (S, [repmat(1.9, numel (g1), 1), g1]);
%! [~, tB] = predict (S, [g2, repmat(4.8, numel (g2), 1)]);
%! assert_equal (max (tB(:,2)) - t0(1,2), 1.33333333333333, 1e-12);
%! assert_equal (max (tA(:,2)) - t0(1,2), 1.04497354497354, 1e-12);

## The intercept is the working response's own weighted mean and not the base
## log-odds of the response.  A balanced fixture cannot tell them apart, both
## being zero, so this one is unbalanced: 42 against 28, whose log-odds is
## -0.40546510810816 and whose first-cycle intercept R2024a reports as exactly
## -0.4.  Two and three cycles are pinned too, the weights being refreshed at
## the start of each.
%!test
%! load fisheriris
%! ii = [51:100, 101:120]';
%! X = meas(ii,1:2);
%! Y = species(ii);
%! f = 1:4:numel (ii);
%! t = Y(f);
%! t(strcmp (t, 'versicolor')) = {'ZZ'};
%! t(strcmp (t, 'virginica')) = {'versicolor'};
%! t(strcmp (t, 'ZZ')) = {'virginica'};
%! Y(f) = t;
%! assert_equal (sum (strcmp (Y, 'virginica')), 28);
%! M1 = fitcgam (X, Y, 'NumTreesPerPredictor', 1, ...
%!               'MaxNumSplitsPerPredictor', 1, 'Interactions', 0);
%! assert_equal (M1.Intercept, -0.4, 1e-12);
%! M2 = fitcgam (X, Y, 'NumTreesPerPredictor', 2, ...
%!               'MaxNumSplitsPerPredictor', 1, 'Interactions', 0);
%! assert_equal (M2.Intercept, -0.43634347574616, 1e-11);
%! M3 = fitcgam (X, Y, 'NumTreesPerPredictor', 3, ...
%!               'MaxNumSplitsPerPredictor', 1, 'Interactions', 0);
%! assert_equal (M3.Intercept, -0.4460203997322, 1e-11);

## The coding of a numeric response is a property of the class, not of either
## engine: the boosted-tree engine reads labels of 1 and 2, or 5 and 9, the
## same way and returns them in their own values.
##
## Four observations cannot be split at all under the tree learner's minimum
## leaf size, so the boosted fit here is a constant and every row comes back
## as the first class.  R2024a answers exactly the same, [1 1 1 1], [5 5 5 5]
## and [0 0 0 0], and takes its first split at ten observations as this does.
## What the test pins is therefore the coding and not the separation: a
## response of 5 and 9 comes back as 5 rather than as a 0/1 index.  The
## spline engine above, which does not fit trees, still separates these four.
%!test
%! Xn = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! assert_equal (predict (fitcgam (Xn, [1; 1; 2; 2]), Xn), [1; 1; 1; 1]);
%! assert_equal (predict (fitcgam (Xn, [5; 5; 9; 9]), Xn), [5; 5; 5; 5]);
%! assert_equal (predict (fitcgam (Xn, [0; 0; 1; 1]), Xn), [0; 0; 0; 0]);

## Ten observations is where a boosted fit first splits, two leaves of the
## minimum five.  Measured against R2024a, which does the same.
%!test
%! Y9 = repmat ({'a'}, 9, 1); Y9(6:9) = {'b'};
%! M9 = fitcgam ((1:9)', Y9, 'NumTreesPerPredictor', 1, ...
%!               'MaxNumSplitsPerPredictor', 1, 'Interactions', 0);
%! M9.ScoreTransform = 'none';
%! [~, s9] = predict (M9, (1:9)');
%! assert_equal (max (s9(:,2)) - min (s9(:,2)), 0, 1e-12);
%! Y10 = repmat ({'a'}, 10, 1); Y10(6:10) = {'b'};
%! M10 = fitcgam ((1:10)', Y10, 'NumTreesPerPredictor', 1, ...
%!                'MaxNumSplitsPerPredictor', 1, 'Interactions', 0);
%! M10.ScoreTransform = 'none';
%! [~, s10] = predict (M10, (1:10)');
%! assert_equal (max (s10(:,2)) - min (s10(:,2)) > 1, true);

%!shared x, y, obj
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1; 4, 5, 6];
%! y = [0; 0; 1; 1; 0];
%! obj = fitcgam (x, y);
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (obj);
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (CVMdl.KFold == 5, true)
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationGAM")
%! assert_equal (CVMdl.CrossValidatedModel, "GAM")
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (obj, 'KFold', 2);
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (CVMdl.KFold == 2, true)
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationGAM")
%! assert_equal (CVMdl.CrossValidatedModel, "GAM")
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (obj, 'HoldOut', 0.2);
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationGAM")
%! assert_equal (CVMdl.CrossValidatedModel, "GAM")
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! partition = cvpartition (y, 'KFold', 3);
%! warning (status);
%! CVMdl = crossval (obj, 'cvPartition', partition);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal (CVMdl.KFold == 3, true)
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationGAM")
%! assert_equal (CVMdl.CrossValidatedModel, "GAM")

## Test input validation for crossval method
%!error<ClassificationGAM.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (obj, 'kfold')
%!error<ClassificationGAM.crossval: specify only one of the optional Name-Value paired arguments.>...
%! crossval (obj, 'kfold', 12, 'holdout', 0.2)
%!error<ClassificationGAM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (obj, 'kfold', 'a')
%!error<ClassificationGAM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (obj, 'holdout', 2)
%!error<ClassificationGAM.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (obj, 'leaveout', 1)
%!error<ClassificationGAM.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (obj, 'cvpartition', 1)
%!error <ClassificationGAM.savemodel: too few input arguments.> ...
%! savemodel (ClassificationGAM ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]))
%!error <ClassificationGAM.savemodel: FNAME must be a character vector.> ...
%! savemodel (ClassificationGAM ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), 1)
%!error <ClassificationGAM.savemodel: FNAME must be a character vector.> ...
%! savemodel (ClassificationGAM ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), ['ab'; 'cd'])

## A ScoreTransform can be assigned, and is stored as a function handle.
%!test
%! Mdl = fitcgam ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]);
%! Mdl.ScoreTransform = 'symmetric';
%! assert_equal (class (Mdl.ScoreTransform), 'char');
%! assert_equal (Mdl.ScoreTransform, 'symmetric');
%! assert_equal (Mdl.STfun (0.25), -0.5);

## The new properties MATLAB reports are carried and saved.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds), 'FitMethod', 'splines');
%! assert_equal (Mdl.Intercept, Mdl.BaseModel.Intercept);
%! assert_equal (size (Mdl.W), [Mdl.NumObservations, 1]);
%! assert_equal (sum (Mdl.W), 1, 1e-12);
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (Mdl.ExpandedPredictorNames, Mdl.PredictorNames);

## ClassNames keeps the type of the response, as MATLAB has it.
%!test
%! assert_equal (fitcgam ([1;2;3;4], [7;3;7;3]).ClassNames, [3; 7]);
%! assert_equal (fitcgam ([1;2;3;4], logical ([1;0;1;0])).ClassNames, ...
%!               logical ([0; 1]));
%! assert_equal (fitcgam ([1;2;3;4], {'b';'a';'b';'a'}).ClassNames, ...
%!               {'a'; 'b'});
%! assert_equal (fitcgam ([1;2;3;4], ['b';'a';'b';'a']).ClassNames, ['a';'b']);

## An assigned ScoreTransform reaches the scores predict returns, and it
## composes on the raw log-odds rather than on the posteriors, as measured
## on MATLAB R2024a.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds));
%! Mdl.ScoreTransform = 'none';
%! [~, raw] = predict (Mdl, meas(inds,:));
%! Mdl.ScoreTransform = 'symmetric';
%! [~, s1] = predict (Mdl, meas(inds,:));
%! assert_equal (s1, 2 * raw - 1, 1e-12);

## The edge is the mean of the margins, and weights reweight that mean.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Y = species(inds);
%! Mdl = fitcgam (X, Y);
%! m = margin (Mdl, X, Y);
%! assert_equal (edge (Mdl, X, Y), mean (m), 1e-12);
%! ## The weights are normalized within each class to that class's prior,
%! ## not divided by their total, so a class keeps the influence its prior
%! ## gives it however the weights inside it are spread.
%! w = [ones(50, 1); 3 * ones(50, 1)];
%! wn = w;
%! wn(1:50) = w(1:50) / sum (w(1:50)) * Mdl.Prior(1);
%! wn(51:100) = w(51:100) / sum (w(51:100)) * Mdl.Prior(2);
%! assert_equal (edge (Mdl, X, Y, 'Weights', w), ...
%!               sum (wn .* m) / sum (wn), 1e-12);
%! ## Scaling a whole class leaves the edge where it was.
%! assert_equal (edge (Mdl, X, Y, 'Weights', w), ...
%!               edge (Mdl, X, Y, 'Weights', [ones(50,1); ones(50,1)]), 1e-12);

## The resubstitution methods are the assessment methods on the training data.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Y = species(inds);
%! Mdl = fitcgam (X, Y);
%! assert_equal (resubPredict (Mdl), predict (Mdl, X));
%! assert_equal (resubMargin (Mdl), margin (Mdl, X, Y));
%! assert_equal (resubEdge (Mdl), edge (Mdl, X, Y));
%! assert_equal (resubLoss (Mdl), loss (Mdl, X, Y));

## Every loss function returns a finite scalar.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Y = species(inds);
%! Mdl = fitcgam (X, Y);
%! names = {'binodeviance', 'classifcost', 'classiferror', 'exponential', ...
%!          'hinge', 'logit', 'mincost', 'quadratic'};
%! for k = 1:numel (names)
%!   L = loss (Mdl, X, Y, 'LossFun', names{k});
%!   assert_equal (isscalar (L) && isfinite (L), true);
%! endfor

## A 'Cost' matching the number of classes is accepted.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds));
%! Mdl.Cost = [0, 2; 5, 0];
%! assert_equal (Mdl.Cost, [0, 2; 5, 0]);

## A saved and reloaded model carries every property, and predicts alike.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Mdl = fitcgam (X, species(inds), 'FitMethod', 'splines', ...
%!                'Interactions', 'all', 'NumIterations', 20);
%! Mdl.ScoreTransform = 'symmetric';
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! Mdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (Mdl2.Intercept, Mdl.Intercept);
%! assert_equal (Mdl2.W, Mdl.W);
%! assert_equal (Mdl2.ExpandedPredictorNames, Mdl.ExpandedPredictorNames);
%! assert_equal (Mdl2.ModelwInt.Parameters(1).coefs, ...
%!               Mdl.ModelwInt.Parameters(1).coefs);
%! assert_equal (Mdl2.ScoreTransform, 'symmetric');
%! assert_equal (Mdl2.STfun (0.25), -0.5);
%! [label, score] = predict (Mdl, X);
%! [label2, score2] = predict (Mdl2, X);
%! assert_equal (label2, label);
%! assert_equal (score2, score);

%!shared x, y, Mdl
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! x = meas(inds,:);
%! y = species(inds);
%! Mdl = fitcgam (x, y);

## Test input validation for margin method
%!error<ClassificationGAM.margin: too few input arguments.> ...
%! margin (Mdl, x)
%!error<ClassificationGAM.margin: X is empty.> ...
%! margin (Mdl, [], y)
%!error<ClassificationGAM.margin: X must have the same number of predictors as the trained model.> ...
%! margin (Mdl, 1, y)
%!error<ClassificationGAM.margin: Y is empty.> ...
%! margin (Mdl, x, [])
%!error<ClassificationGAM.margin: Y must have the same number of rows as X.> ...
%! margin (Mdl, x, y(1:10))

## Test input validation for edge method
%!error<ClassificationGAM.edge: too few input arguments.> ...
%! edge (Mdl, x)
%!error<ClassificationGAM.edge: Name-Value arguments must be in pairs.> ...
%! edge (Mdl, x, y, 'Weights')
%!error<ClassificationGAM.edge: invalid parameter name in optional paired arguments.> ...
%! edge (Mdl, x, y, 'LossFun', 'hinge')
%!error<ClassificationGAM.edge: 'Weights' must be a numeric vector.> ...
%! edge (Mdl, x, y, 'Weights', 'a')
%!error<ClassificationGAM.edge: size of 'Weights' must equal the number of rows in X.> ...
%! edge (Mdl, x, y, 'Weights', [1, 2, 3])

## Test input validation for loss method
%!error<ClassificationGAM.loss: too few input arguments.> ...
%! loss (Mdl, x)
%!error<ClassificationGAM.loss: Name-Value arguments must be in pairs.> ...
%! loss (Mdl, x, y, 'LossFun')
%!error<ClassificationGAM.loss: 'LossFun' must be a character vector.> ...
%! loss (Mdl, x, y, 'LossFun', 1)
%!error<ClassificationGAM.loss: unsupported Loss function.> ...
%! loss (Mdl, x, y, 'LossFun', 'nonsense')

## RowsUsed is empty when every observation was used.
%!test
%! load fisheriris
%! X = meas(1:100,:);
%! Y = grp2idx (species(1:100));
%! Mdl = fitcgam (X, Y);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (class (Mdl.RowsUsed), 'double');
%! assert_equal (Mdl.NumObservations, 100);
%! assert_equal (rows (Mdl.X), 100);
%! assert_equal (rows (Mdl.W), 100);

## A missing response drops its observation and RowsUsed marks it.
%!test
%! load fisheriris
%! X = meas(1:100,:);
%! Y = grp2idx (species(1:100));
%! Y(5) = NaN;
%! Mdl = fitcgam (X, Y);
%! assert_equal (class (Mdl.RowsUsed), 'logical');
%! assert_equal (size (Mdl.RowsUsed), [100, 1]);
%! assert_equal (sum (Mdl.RowsUsed), 99);
%! assert_equal (Mdl.RowsUsed(5), false);
%! assert_equal (Mdl.NumObservations, 99);
%! assert_equal (rows (Mdl.X), 99);
%! assert_equal (rows (Mdl.W), 99);

## A missing predictor keeps its observation, so RowsUsed stays empty.
%!test
%! load fisheriris
%! X = meas(1:100,:);
%! X(3,2) = NaN;
%! Y = grp2idx (species(1:100));
%! Mdl = fitcgam (X, Y);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (Mdl.NumObservations, 100);
%! assert_equal (rows (Mdl.X), 100);
%! assert_equal (sum (isnan (Mdl.X(:))), 1);

## The prior reweights the observations of each class.  Values from R2024a.
%!test
%! load fisheriris
%! i2 = [1:50, 51:80];
%! Mdl = fitcgam (meas(i2,:), species(i2));
%! assert_equal (Mdl.Prior, [0.625, 0.375], 1e-14);
%! assert_equal (Mdl.W(1), 0.0125, 1e-14);
%! Mdl = fitcgam (meas(i2,:), species(i2), 'Prior', 'uniform');
%! assert_equal (Mdl.W(1), 0.01, 1e-14);
%! assert_equal (Mdl.W(51), 1/60, 1e-14);

## A fitted model survives savemodel and loadmodel: the properties come
## back as they were and it predicts the same.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds));
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'ClassificationGAM');
%! assert_equal (M2.NumObservations, Mdl.NumObservations);
%! assert_equal (M2.PredictorNames, Mdl.PredictorNames);
%! assert_equal (class (M2.ScoreTransform), class (Mdl.ScoreTransform));
%! assert_equal (predict (M2, meas(1:5,:)), predict (Mdl, meas(1:5,:)));

## The default transform is 'logit', as MATLAB's is, so the scores predict
## reports are posterior probabilities that sum to one.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds));
%! assert_equal (Mdl.ScoreTransform, 'logit');
%! [~, scores] = predict (Mdl, meas(1:6,:));
%! assert_equal (sum (scores, 2), ones (6, 1), 1e-12);
%! assert_equal (all (scores(:) >= 0 & scores(:) <= 1), true);

## The untransformed score is the log-odds pair whose columns sum to zero,
## and the default transform is what maps it to the posteriors.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds));
%! [~, post] = predict (Mdl, meas(1:6,:));
%! Mdl.ScoreTransform = 'none';
%! [~, raw] = predict (Mdl, meas(1:6,:));
%! assert_equal (sum (raw, 2), zeros (6, 1), 1e-12);
%! assert_equal (1 ./ (1 + exp (-raw)), post, 1e-12);

## BinEdges reports the cut points the boosted-tree engine binned each
## predictor at, as MATLAB's generalized additive model does.  Under the
## spline engine, which does no binning, it stays the empty cell MATLAB
## reports for every learner that does none.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds));
%! assert_equal (class (Mdl.BinEdges), 'cell');
%! assert_equal (numel (Mdl.BinEdges), 4);
%! assert_equal (numel (Mdl.BinEdges{1}), 27);
%! Msp = fitcgam (meas(inds,:), species(inds), 'FitMethod', 'splines');
%! assert_equal (Msp.BinEdges, {});

## The shared cost guard is in force here too, and the struct form is
## permuted into this model's class order.  The battery is on
## ClassificationDiscriminant.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds));
%! S = struct ('ClassNames', {{'versicolor'; 'setosa'}}, ...
%!             'ClassificationCosts', [0, 1; 2, 0]);
%! Mdl.Cost = S;
%! assert_equal (Mdl.Cost, [0, 2; 1, 0]);

%!error<ClassificationGAM: 'Cost' must have zeros on its diagonal.> ...
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds));
%! Mdl.Cost = ones (2);

## The boosted-tree engine is reachable by name while the default is still
## the spline engine, and it reports the surface MATLAB reports.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds), 'FitMethod', 'boostedtrees');
%! assert_equal (Mdl.FitMethod, 'boostedtrees');
%! assert_equal (numel (Mdl.BinEdges), 4);
%! assert_equal (numel (Mdl.BinEdges{1}), 27);
%! assert_equal (numel (fieldnames (Mdl.ModelParameters)), 13);
%! assert_equal (Mdl.ModelParameters.Type, 'classification');
%! assert_equal (Mdl.ModelParameters.Method, 'GAM');

## Its defaults are MATLAB's own, reported through ModelParameters.
%!test
%! Mdl = fitcgam ([1, 2; 2, 3; 3, 4; 4, 5; 5, 6; 6, 7], [0;0;0;1;1;1], ...
%!                'FitMethod', 'boostedtrees');
%! MP = Mdl.ModelParameters;
%! assert_equal (MP.NumTreesPerPredictor, 300);
%! assert_equal (MP.NumTreesPerInteraction, 100);
%! assert_equal (MP.MaxNumSplitsPerPredictor, 1);
%! assert_equal (MP.MaxNumSplitsPerInteraction, 4);
%! assert_equal (MP.InitialLearnRateForPredictors, 1);
%! assert_equal (MP.InitialLearnRateForInteractions, 1);
%! assert_equal (MP.MaxPValue, 1);
%! assert_equal (MP.NumPrint, 10);
%! assert_equal (MP.VerbosityLevel, 0);

## Each phase reports why it stopped; a phase that never ran reports nothing.
%!test
%! Mdl = fitcgam ([1, 2; 2, 3; 3, 4; 4, 5; 5, 6; 6, 7], [0;0;0;1;1;1], ...
%!                'FitMethod', 'boostedtrees');
%! assert_equal (isfield (Mdl.ReasonForTermination, 'PredictorTrees'), true);
%! assert_equal (isfield (Mdl.ReasonForTermination, 'InteractionTrees'), true);
%! assert_equal (Mdl.ReasonForTermination.InteractionTrees, '');

## Interactions are detected, held on their own coarse grid, and the second
## phase reports its own termination.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds), 'FitMethod', 'boostedtrees', ...
%!                'Interactions', 'all');
%! assert_equal (rows (Mdl.Interactions), 6);
%! assert_equal (columns (Mdl.Interactions), 2);
%! assert_equal (numel (Mdl.PairDetectionBinEdges{1}), 7);
%! assert_equal (! isempty (Mdl.ReasonForTermination.InteractionTrees), true);

## A tree-fitted model predicts, and its scores are posteriors summing to one
## under the default transform.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Mdl = fitcgam (X, species(inds), 'FitMethod', 'boostedtrees');
%! [label, score] = predict (Mdl, X);
%! assert_equal (numel (label), rows (X));
%! assert_equal (sum (score, 2), ones (rows (X), 1), 1e-12);

## compact carries the engine and every property it needs to predict alike.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Mdl = fitcgam (X, species(inds), 'FitMethod', 'boostedtrees');
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.FitMethod, 'boostedtrees');
%! assert_equal (CMdl.BinEdges, Mdl.BinEdges);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));

## savemodel and loadmodel carry the tree fit, so a reloaded model predicts
## identically rather than coming back with an empty engine.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Mdl = fitcgam (X, species(inds), 'FitMethod', 'boostedtrees');
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (M2.FitMethod, 'boostedtrees');
%! assert_equal (M2.TreeModel.ShapeValues, Mdl.TreeModel.ShapeValues);
%! assert_equal (predict (M2, X), predict (Mdl, X));

## The spline engine is unchanged and still reachable by name.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds), 'FitMethod', 'splines');
%! assert_equal (Mdl.FitMethod, 'splines');
%! assert_equal (Mdl.BinEdges, {});
%! assert_equal (isempty (Mdl.TreeModel), true);
%! assert_equal (Mdl.Knots, [5, 5, 5, 5]);

## An argument belonging to the other engine is refused, not ignored.
%!error<ClassificationGAM: 'FitMethod' must be 'boostedtrees' or 'splines'.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'nonsense')
%!error<ClassificationGAM: 'FitMethod' must be 'boostedtrees' or 'splines'.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 5)
%!error<ClassificationGAM: 'knots' is a parameter of the spline engine and cannot be used with 'FitMethod' 'boostedtrees'.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'boostedtrees', 'Knots', 4)
%!error<ClassificationGAM: 'maxpvalue' is a parameter of the boosted-tree engine and cannot be used with 'FitMethod' 'splines'.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'splines', 'MaxPValue', 0.5)
%!error<ClassificationGAM: 'NumTreesPerPredictor' must be a positive integer value.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'boostedtrees', ...
%!          'NumTreesPerPredictor', 0)
%!error<ClassificationGAM: 'NumTreesPerInteraction' must be a positive integer value.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'boostedtrees', ...
%!          'NumTreesPerInteraction', 1.5)
%!error<ClassificationGAM: 'MaxNumSplitsPerPredictor' must be a positive integer value.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'boostedtrees', ...
%!          'MaxNumSplitsPerPredictor', -1)
%!error<ClassificationGAM: 'MaxNumSplitsPerInteraction' must be a positive integer value.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'boostedtrees', ...
%!          'MaxNumSplitsPerInteraction', 'a')
%!error<ClassificationGAM: 'InitialLearnRateForPredictors' must be greater than 0 and at most 1.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'boostedtrees', ...
%!          'InitialLearnRateForPredictors', 0)
%!error<ClassificationGAM: 'InitialLearnRateForInteractions' must be greater than 0 and at most 1.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'boostedtrees', ...
%!          'InitialLearnRateForInteractions', 2)
%!error<ClassificationGAM: 'Verbose' must be a non-negative integer value.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'boostedtrees', 'Verbose', -1)
%!error<ClassificationGAM: 'NumPrint' must be a positive integer value.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'boostedtrees', 'NumPrint', 0)
%!error<ClassificationGAM: 'MaxPValue' must be between 0 and 1.> ...
%! fitcgam ([1;2;3;4], [0;0;1;1], 'FitMethod', 'boostedtrees', 'MaxPValue', 2)

## HyperparameterOptimizationResults is declared for MATLAB compatibility and
## stays empty, this class running no search over its hyperparameters.
%!test
%! load fisheriris
%! b = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(b,:), species(b));
%! assert_equal (isempty (Mdl.HyperparameterOptimizationResults), true);
