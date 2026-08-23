## Copyright (C) 2024-2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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

classdef CompactClassificationGAM
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactClassificationGAM
  ##
  ## Compact generalized additive model classification
  ##
  ## The @code{CompactClassificationGAM} class is a compact version of a
  ## Generalized Additive Model classifier, @code{ClassificationGAM}.  It does
  ## not include the training data, resulting in a smaller classifier size that
  ## can be used for making predictions from new data, but not for tasks such as
  ## cross validation.
  ##
  ## A @code{CompactClassificationGAM} object can only be created from a
  ## @code{ClassificationGAM} model by using the @code{compact} method.
  ##
  ## The weak learner here is a smoothing spline, one per predictor, where
  ## MATLAB's generalized additive model boosts shallow decision trees, so
  ## predictions and scores will not agree with MATLAB's on the same data.
  ## The two are different estimators of the same additive structure, and the
  ## difference is deliberate.  The spline fit is described by @code{Knots},
  ## @code{Order}, @code{DoF}, @code{Formula}, @code{LearningRate},
  ## @code{NumIterations}, @code{BaseModel}, @code{ModelwInt} and
  ## @code{IntMatrix}, which MATLAB's compact model does not carry.
  ##
  ## @seealso{ClassificationGAM, fitcgam}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} NumPredictors
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
    ## @deftp {CompactClassificationGAM} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which they appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames  = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName    = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} ClassNames
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
    ## @deftp {CompactClassificationGAM} {property} Prior
    ##
    ## Prior probability for each class
    ##
    ## A 2-element numeric vector specifying the prior probabilities for each
    ## class.  The order of the elements in @qcode{Prior} corresponds to the
    ## order of the classes in @qcode{ClassNames}.  This property is read-only.
    ##
    ## @end deftp
    Prior           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} Formula
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
    ## @deftp {CompactClassificationGAM} {property} Interactions
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
    ## @deftp {CompactClassificationGAM} {property} Knots
    ##
    ## Knots for spline fitting
    ##
    ## A scalar or row vector specifying the number of knots for each predictor
    ## variable in the spline fitting.  This property is read-only.
    ##
    ## @end deftp
    Knots           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} Order
    ##
    ## Order of spline fitting
    ##
    ## A scalar or row vector specifying the order of the spline for each
    ## predictor variable.  This property is read-only.
    ##
    ## @end deftp
    Order           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} DoF
    ##
    ## Degrees of freedom for spline fitting
    ##
    ## A scalar or row vector specifying the degrees of freedom for each
    ## predictor variable in the spline fitting.  This property is read-only.
    ##
    ## @end deftp
    DoF             = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} LearningRate
    ##
    ## Learning rate for gradient boosting
    ##
    ## A scalar value between 0 and 1 specifying the learning rate used in the
    ## gradient boosting algorithm.  This property is read-only.
    ##
    ## @end deftp
    LearningRate    = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} NumIterations
    ##
    ## Maximum number of iterations
    ##
    ## A positive integer specifying the maximum number of iterations for the
    ## gradient boosting algorithm.  This property is read-only.
    ##
    ## @end deftp
    NumIterations   = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} BaseModel
    ##
    ## Base model parameters
    ##
    ## A structure containing the parameters of the base model without any
    ## interaction terms.  The base model represents the generalized additive
    ## model with only the main effects (predictor terms) included.
    ## This property is read-only.
    ##
    ## @end deftp
    BaseModel       = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} ModelwInt
    ##
    ## Model parameters with interactions
    ##
    ## A structure containing the parameters of the model that includes
    ## interaction terms.  This model extends the base model by adding
    ## interaction terms between predictors.  This property is read-only.
    ##
    ## @end deftp
    ModelwInt       = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} IntMatrix
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
    IntMatrix       = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} Intercept
    ##
    ## Intercept of the fitted model
    ##
    ## A numeric scalar, the log-odds of the response mean, which every
    ## additive term is measured against.  This property is read-only.
    ##
    ## @end deftp
    Intercept       = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector holding the column of each predictor treated as
    ## categorical, and empty when none is.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} ExpandedPredictorNames
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
  endproperties

  ## Properties a user may set after the model is built.  Each one is
  ## validated by its set method below.
  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} Cost
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
    ## @end deftp
    Cost            = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationGAM} {property} ScoreTransform
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
    ## The default is @qcode{'none'} because this model's @code{predict}
    ## already returns posterior probabilities, which sum to one across the
    ## classes.  MATLAB's generalized additive model is built from boosted
    ## trees whose raw scores are log-odds, so it defaults to
    ## @qcode{'logit'} and reports probabilities only after applying it.
    ## Both return the same kind of quantity while naming a different
    ## transform, so do not set @qcode{'logit'} here expecting to match
    ## MATLAB: applying it to numbers that are already probabilities
    ## leaves them summing to about 1.23 rather than one.
    ##
    ## @end deftp
    ScoreTransform  = 'none';

  endproperties

  ## Readable by the counterpart class, which copies it, and kept out of
  ## the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    STfun = @(x) x;
  endproperties

  ## Set methods for the properties a user may assign.
  methods (Hidden)

    function this = set.Cost (this, val)
      gnY = this.ClassNames;
      if (isempty (val))
        this.Cost = cast (! eye (classCount (gnY)), 'double');
      else
        K = classCount (gnY);
        if (! isequal (size (val), [K, K]))
          error (strcat ("CompactClassificationGAM: the number", ...
                         " of rows and columns in 'Cost' must", ...
                         " correspond to selected classes in Y."));
        endif
        this.Cost = val;
      endif
    endfunction

    function this = set.ScoreTransform (this, val)
      [f, nm] = parseScoreTransform (val, 'CompactClassificationGAM');
      this.ScoreTransform = nm;
      this.STfun = f;
    endfunction

    ## Class object constructor
    function this = CompactClassificationGAM (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), 'ClassificationGAM'))
        error ("CompactClassificationGAM: invalid classification object.");
      endif

      ## Save properties to compact model
      this.NumPredictors   = Mdl.NumPredictors;
      this.PredictorNames  = Mdl.PredictorNames;
      this.ResponseName    = Mdl.ResponseName;
      this.ClassNames      = Mdl.ClassNames;

      this.Cost            = Mdl.Cost;
      this.Prior           = Mdl.Prior;
      this.ScoreTransform  = Mdl.ScoreTransform;
      this.STfun          = Mdl.STfun;

      this.Formula         = Mdl.Formula;
      this.Interactions    = Mdl.Interactions;
      this.Knots           = Mdl.Knots;
      this.Order           = Mdl.Order;
      this.DoF             = Mdl.DoF;
      this.LearningRate    = Mdl.LearningRate;
      this.NumIterations   = Mdl.NumIterations;
      this.BaseModel       = Mdl.BaseModel;
      this.ModelwInt       = Mdl.ModelwInt;
      this.IntMatrix       = Mdl.IntMatrix;
      this.Intercept       = Mdl.Intercept;
      this.CategoricalPredictors  = Mdl.CategoricalPredictors;
      this.ExpandedPredictorNames = Mdl.ExpandedPredictorNames;

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
      fprintf ("\n  CompactClassificationGAM\n\n");
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
    ## @deftypefn  {CompactClassificationGAM} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationGAM} {[@var{label}, @var{score}] =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationGAM} {[@var{label}, @var{score}] =} predict (@dots{}, @qcode{'IncludeInteractions'}, @var{includeInteractions})
    ##
    ## Predict labels for new data using the Generalized Additive Model (GAM)
    ## stored in a CompactClassificationGAM object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the predicted
    ## labels for the data in @var{XC} based on the model stored in the
    ## CompactClassificationGAM object, @var{obj}.
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
    ## @var{obj} must be a @qcode{CompactClassificationGAM} class object.
    ## @item
    ## @var{XC} must be an @math{M*P} numeric matrix where each row is an
    ## observation and each column corresponds to a predictor variable.
    ## @item
    ## @var{includeInteractions} is a logical scalar indicating whether to
    ## include interaction terms in the predictions.
    ## @end itemize
    ##
    ## @seealso{CompactClassificationGAM, ClassificationGAM, fitcgam}
    ## @end deftypefn
    function [labels, scores] = predict (this, XC, varargin)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("CompactClassificationGAM.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("CompactClassificationGAM.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("CompactClassificationGAM.predict: XC must have", ...
                       " the same number of features as the trained model."));
      endif

      ## Clean XC data
      notnansf  = ! logical (sum (isnan (XC), 2));
      XC        = XC(notnansf, :);

      ## Default values for Name-Value Pairs
      incInt = ! isempty (this.IntMatrix);
      Cost = this.Cost;

      ## Parse optional arguments
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'includeinteractions'
            tmpInt = varargin{2};
            if (! islogical (tmpInt) || (tmpInt != 0 && tmpInt != 1))
              error (strcat ("CompactClassificationGAM.predict:", ...
                     " includeinteractions must be a logical value."));
            endif
            ## Check model for interactions
            if (tmpInt && isempty (this.IntMatrix))
              error (strcat ("CompactClassificationGAM.predict: trained", ...
                             " model does not include any interactions."));
            endif
            incInt = tmpInt;

          otherwise
            error (strcat ("CompactClassificationGAM.predict: invalid", ...
                           " NAME in optional pairs of arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

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

      ## Predict probabilities from testing data
      scores = predict_val (params, XC, Interc);

      ## Compute the expected misclassification cost matrix
      numObservations = size (XC, 1);
      CE = zeros (numObservations, 2);

      for k = 1:2
        for i = 1:2
          CE(:, k) = CE(:, k) + scores(:, i) * Cost(k, i);
        endfor
      endfor

      ## Select the class with the minimum expected misclassification cost
      [~, minIdx] = min (CE, [], 2);
      labels = labelsFromIndex (this.ClassNames, minIdx);

      ## Apply ScoreTransform
      scores = this.STfun (scores);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationGAM} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin of a compact generalized additive model.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns a column
    ## vector holding, for each row of @var{X}, the score the model gives its
    ## true class in @var{Y} less the score it gives the other class.  A
    ## positive margin means the observation is classified correctly, and the
    ## larger it is the more confidently so.
    ##
    ## @seealso{CompactClassificationGAM, ClassificationGAM, edge, loss,
    ## predict}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("CompactClassificationGAM.margin: too few input arguments.");
      endif

      [X, Y] = checkXY_ (this, X, Y, "margin");

      [~, scores] = predict (this, X);
      classes = this.ClassNames;
      m = zeros (rows (X), 1);
      for i = 1:rows (X)
        idx = labelIndex (classes, Y, i);
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
    ## @deftypefn  {CompactClassificationGAM} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationGAM} {@var{e} =} edge (@dots{}, @qcode{"Weights"}, @var{w})
    ##
    ## Classification edge of a compact generalized additive model.
    ##
    ## @code{@var{e} = edge (@var{obj}, @var{X}, @var{Y})} returns the mean of
    ## the classification margins over the rows of @var{X}.
    ##
    ## @code{@var{e} = edge (@dots{}, @qcode{"Weights"}, @var{w})} takes the
    ## weighted mean instead, with one weight per row of @var{X}.
    ##
    ## @seealso{CompactClassificationGAM, ClassificationGAM, margin, loss,
    ## predict}
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("CompactClassificationGAM.edge: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactClassificationGAM.edge: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      [X, Y] = checkXY_ (this, X, Y, "edge");

      ## The weights are normalized within each class to that class's prior,
      ## which is what the oracle does and is not the same as dividing by
      ## their total.  This used to divide by the total.
      ## The weights are parsed before anything is computed, so a bad
      ## Name-Value pair is reported as such rather than after a margin.
      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       "CompactClassificationGAM", "edge");
      m = margin (this, X, Y);
      e = sum (W .* m(:)) / sum (W);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationGAM} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationGAM} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of a compact generalized additive model.
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
    ## @seealso{CompactClassificationGAM, ClassificationGAM, margin, edge,
    ## predict}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("CompactClassificationGAM.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactClassificationGAM.loss: Name-Value", ...
                       " arguments must be in pairs."));
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
            error (strcat ("CompactClassificationGAM.loss: 'LossFun'", ...
                           " must be a character vector."));
          endif
          LossFun = tolower (LossFun);
          if (! any (strcmpi (LossFun, lossnames)))
            error ("CompactClassificationGAM.loss: unsupported Loss function.");
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
      for i = 1:rows (X)
        idx = labelIndex (classes, Y, i);
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
          for i = 1:rows (X)
            [~, k] = min (scores(i,:) * this.Cost);
            true_idx = labelIndex (classes, Y, i);
            L = L + W(i) * this.Cost(true_idx, k);
          endfor
        case 'classifcost'
          ## What the model's own prediction costs, given the true class
          L = 0;
          for i = 1:rows (X)
            true_idx = labelIndex (classes, Y, i);
            pred_idx = find (ismember (classes, label(i)));
            L = L + W(i) * this.Cost(true_idx, pred_idx);
          endfor
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationGAM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a CompactClassificationGAM object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## CompactClassificationGAM object into an Octave binary file, the name of
    ## which is specified in @var{filename}, along with an extra variable,
    ## which defines the type classification object these variables constitute.
    ## Use @code{loadmodel} in order to load a classification object into
    ## Octave's workspace.
    ##
    ## @seealso{loadmodel, fitcgam, ClassificationGAM, CompactClassificationGAM}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error ("CompactClassificationGAM.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error ("CompactClassificationGAM.savemodel: FNAME must be a character vector.");
      endif
      ## Generate variable for class name
      classdef_name = 'CompactClassificationGAM';

      ## Create variables from model properties
      NumPredictors   = this.NumPredictors;
      PredictorNames  = this.PredictorNames;
      ResponseName    = this.ResponseName;
      ClassNames      = this.ClassNames;
      Prior           = this.Prior;
      Cost            = this.Cost;
      ScoreTransform  = this.ScoreTransform;
      Intercept       = this.Intercept;
      CategoricalPredictors  = this.CategoricalPredictors;
      ExpandedPredictorNames = this.ExpandedPredictorNames;
      STfun          = this.STfun;
      Formula         = this.Formula;
      Interactions    = this.Interactions;
      Knots           = this.Knots;
      Order           = this.Order;
      DoF             = this.DoF;
      LearningRate    = this.LearningRate;
      NumIterations   = this.NumIterations;
      BaseModel       = this.BaseModel;
      ModelwInt       = this.ModelwInt;
      IntMatrix       = this.IntMatrix;

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'NumPredictors', ...
            'PredictorNames', 'ResponseName', 'ClassNames', 'Prior', 'Cost', ...
            'ScoreTransform', 'STfun', 'Intercept', ...
            'CategoricalPredictors', 'ExpandedPredictorNames', ...
            'Formula', 'Interactions', 'Knots', ...
            'Order', 'DoF', 'BaseModel', 'ModelwInt', 'IntMatrix', ...
            'LearningRate', 'NumIterations');
    endfunction

  endmethods

  methods (Access = private)

    ## Shared validation for the assessment methods, so each reports under
    ## its own name.
    function [X, Y] = checkXY_ (this, X, Y, caller)
      if (isempty (X))
        error ("CompactClassificationGAM.%s: X is empty.", caller);
      elseif (this.NumPredictors != columns (X))
        error (strcat ("CompactClassificationGAM.%s: X must have the", ...
                       " same number of predictors as the trained", ...
                       " model."), caller);
      endif
      if (isempty (Y))
        error ("CompactClassificationGAM.%s: Y is empty.", caller);
      elseif (rows (X) != rows (Y))
        error (strcat ("CompactClassificationGAM.%s: Y must have the", ...
                       " same number of rows as X."), caller);
      endif
    endfunction

    ## Pull a "Weights" pair out of the optional arguments, defaulting to a
    ## uniform weight, and reject any other name.
    function W = getWeights_ (this, args, n, caller)
      W = ones (n, 1);
      for i = 1:2:numel (args)
        if (! (ischar (args{i}) && isrow (args{i})))
          error (strcat ("CompactClassificationGAM.%s: parameter name", ...
                         " must be a character vector."), caller);
        endif
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
          if (! (isnumeric (W) && isvector (W)))
            error (strcat ("CompactClassificationGAM.%s: 'Weights'", ...
                           " must be a numeric vector."), caller);
          endif
          if (numel (W) != n)
            error (strcat ("CompactClassificationGAM.%s: size of", ...
                           " 'Weights' must equal the number of rows", ...
                           " in X."), caller);
          endif
        else
          error (strcat ("CompactClassificationGAM.%s: invalid parameter", ...
                         " name in optional paired arguments."), caller);
        endif
      endfor
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a CompactClassificationGAM object
      mdl = CompactClassificationGAM ();

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
          error ("CompactClassificationGAM.load_model: invalid model in '%s'.", ...
                 filename)
        end_try_catch
      endfor
    endfunction

  endmethods

  methods(Access = private)

    ## Set cost

  endmethods

endclassdef

## Helper function
function scores = predict_val (params, XC, intercept)
  ## The shared prediction engine evaluates every additive term and maps the
  ## sum through the logistic link, returning both class probabilities.
  scores = gampredict (params, XC, intercept, 1);
endfunction

%!demo
%! ## Create a generalized additive model classifier and its compact version
%! # and compare their size
%!
%! load fisheriris
%! X = meas;
%! Y = species;
%!
%! Mdl = fitcdiscr (X, Y, 'ClassNames', unique (species))
%! CMdl = crossval (Mdl)

## Test constructor
%!test
%! Mdl = CompactClassificationGAM ();
%! assert_equal (class (Mdl), "CompactClassificationGAM")
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [0; 0; 1; 1];
%! PredictorNames = {'Feature1', 'Feature2', 'Feature3'};
%! Mdl = fitcgam (x, y, 'PredictorNames', PredictorNames);
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), "CompactClassificationGAM");
%! assert_equal ({CMdl.NumPredictors, CMdl.ResponseName}, {3, 'Y'})
%! assert_equal (CMdl.ClassNames, [0; 1])
%! assert_equal (CMdl.PredictorNames, PredictorNames)
%! assert_equal (CMdl.BaseModel.Intercept, 0)
%!test
%! load fisheriris
%! inds = strcmp (species,'versicolor') | strcmp (species,'virginica');
%! X = meas(inds, :);
%! Y = species(inds, :)';
%! Y = strcmp (Y, 'virginica')';
%! Mdl = fitcgam (X, Y, 'Formula', 'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3');
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), "CompactClassificationGAM");
%! assert_equal ({CMdl.NumPredictors, CMdl.ResponseName}, {4, 'Y'})
%! assert_equal (CMdl.ClassNames, logical ([0; 1]))
%! assert_equal (CMdl.Formula, 'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3')
%! assert_equal (CMdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'})
%! assert_equal (CMdl.ModelwInt.Intercept, 0)
%!test
%! X = [2, 3, 5; 4, 6, 8; 1, 2, 3; 7, 8, 9; 5, 4, 3];
%! Y = [0; 1; 0; 1; 1];
%! Mdl = fitcgam (X, Y, 'Knots', [4, 4, 4], 'Order', [3, 3, 3]);
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), "CompactClassificationGAM");
%! assert_equal ({CMdl.NumPredictors, CMdl.ResponseName}, {3, 'Y'})
%! assert_equal (CMdl.ClassNames, [0; 1])
%! assert_equal (CMdl.PredictorNames, {'x1', 'x2', 'x3'})
%! assert_equal (CMdl.Knots, [4, 4, 4])
%! assert_equal (CMdl.Order, [3, 3, 3])
%! assert_equal (CMdl.DoF, [7, 7, 7])
%! assert_equal (CMdl.BaseModel.Intercept, 0.4055, 1e-1)

## Test input validation for constructor
## A compact model keeps the character class names and predicts whole names.
%!test
%! load fisheriris
%! bch = ! strcmp (species, "setosa");
%! Xch = meas(bch,:); Ycell = species(bch); Ych = char (Ycell);
%! rand ("state", 1); randn ("state", 1); Cc = compact (fitcgam (Xch, Ych));
%! rand ("state", 1); randn ("state", 1); Cs = compact (fitcgam (Xch, Ycell));
%! assert_equal (cellstr (Cc.ClassNames), Cs.ClassNames);
%! assert_equal (cellstr (predict (Cc, Xch)), predict (Cs, Xch));

## and reads one back in its assessment methods.
%!test
%! load fisheriris
%! bch = ! strcmp (species, "setosa");
%! Xch = meas(bch,:); Ycell = species(bch); Ych = char (Ycell);
%! rand ("state", 1); randn ("state", 1); Cc = compact (fitcgam (Xch, Ych));
%! rand ("state", 1); randn ("state", 1); Cs = compact (fitcgam (Xch, Ycell));
%! assert_equal (loss (Cc, Xch, Ych), loss (Cs, Xch, Ycell), 1e-12);

%!error<CompactClassificationGAM: invalid classification object.> ...
%! CompactClassificationGAM (1)

## Test predict method
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10];
%! y = [1; 0; 1; 0; 1];
%! Mdl = fitcgam (x, y, 'interactions', 'all');
%! CMdl = compact (Mdl);
%! l = [1; 0; 1; 0; 1];
%! s = [0.0334, 0.9666; 0.9648, 0.0352; 0.0334, 0.9666; ...
%!      0.9648, 0.0352; 0.0334, 0.9666];
%! [labels, scores] = predict (CMdl, x);
%! assert_equal (class (CMdl), "CompactClassificationGAM");
%! assert_equal ({CMdl.NumPredictors, CMdl.ResponseName}, {2, 'Y'})
%! assert_equal (CMdl.ClassNames, [0; 1])
%! assert_equal (CMdl.PredictorNames, {'x1', 'x2'})
%! assert_equal (CMdl.ModelwInt.Intercept, 0.4055, 1e-1)
%! assert_equal (labels, l)
%! assert_equal (scores, s, 1e-1)
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [0; 0; 1; 1];
%! interactions = [false, true, false; true, false, true; false, true, false];
%! Mdl = fitcgam (x, y, 'learningrate', 0.2, 'interactions', interactions);
%! CMdl = compact (Mdl);
%! [label, score] = predict (CMdl, x, 'includeinteractions', true);
%! l = [0; 0; 1; 1];
%! s = [0.9725, 0.0275; 0.9895, 0.0105; 0.0070, 0.9930; 0.0238, 0.9762];
%! assert_equal (class (CMdl), "CompactClassificationGAM");
%! assert_equal ({CMdl.NumPredictors, CMdl.ResponseName}, {3, 'Y'})
%! assert_equal (CMdl.ClassNames, [0; 1])
%! assert_equal (CMdl.PredictorNames, {'x1', 'x2', 'x3'})
%! assert_equal (CMdl.ModelwInt.Intercept, 0)
%! assert_equal (label, l)
%! assert_equal (score, s, 1e-1)

## Test input validation for predict method
%!shared CMdl
%! Mdl = fitcgam (ones (4,2), ones (4,1));
%! CMdl = compact (Mdl);
%!error<CompactClassificationGAM.predict: too few input arguments.> ...
%! predict (CMdl)
%!error<CompactClassificationGAM.predict: XC is empty.> ...
%! predict (CMdl, [])
%!error<CompactClassificationGAM.predict: XC must have the same number of features as the trained model.> ...
%! predict (CMdl, 1)
%!error <CompactClassificationGAM.savemodel: too few input arguments.> ...
%! savemodel (CompactClassificationGAM ())
%!error <CompactClassificationGAM.savemodel: FNAME must be a character vector.> ...
%! savemodel (CompactClassificationGAM (), 1)
%!error <CompactClassificationGAM.savemodel: FNAME must be a character vector.> ...
%! savemodel (CompactClassificationGAM (), ['ab'; 'cd'])

## A ScoreTransform can be assigned, and is stored as a function handle.
%!test
%! CMdl = compact (fitcgam ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]));
%! CMdl.ScoreTransform = 'symmetric';
%! assert_equal (class (CMdl.ScoreTransform), 'char');
%! assert_equal (CMdl.ScoreTransform, 'symmetric');
%! assert_equal (CMdl.STfun (0.25), -0.5);

## The compact model carries what MATLAB's compact model reports.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcgam (meas(inds,:), species(inds));
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.Intercept, Mdl.Intercept);
%! assert_equal (CMdl.CategoricalPredictors, Mdl.CategoricalPredictors);
%! assert_equal (CMdl.ExpandedPredictorNames, Mdl.ExpandedPredictorNames);

## margin, edge and loss agree with the model it was compacted from.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Y = species(inds);
%! Mdl = fitcgam (X, Y);
%! CMdl = compact (Mdl);
%! assert_equal (margin (CMdl, X, Y), margin (Mdl, X, Y));
%! assert_equal (edge (CMdl, X, Y), edge (Mdl, X, Y));
%! assert_equal (loss (CMdl, X, Y), loss (Mdl, X, Y));

## A saved and reloaded compact model carries every property it holds.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! CMdl = compact (fitcgam (meas(inds,:), species(inds)));
%! fname = tempname ();
%! savemodel (CMdl, fname);
%! CMdl2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (CMdl2.Intercept, CMdl.Intercept);
%! assert_equal (CMdl2.ClassNames, CMdl.ClassNames);
%! assert_equal (CMdl2.BaseModel.Parameters(1).coefs, ...
%!               CMdl.BaseModel.Parameters(1).coefs);
%! assert_equal (predict (CMdl2, meas(inds,:)), predict (CMdl, meas(inds,:)));

%!shared x2, y2, CM
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! x2 = meas(inds,:);
%! y2 = species(inds);
%! CM = compact (fitcgam (x2, y2));

## Test input validation for the assessment methods
%!error<CompactClassificationGAM.margin: too few input arguments.> ...
%! margin (CM, x2)
%!error<CompactClassificationGAM.margin: X is empty.> ...
%! margin (CM, [], y2)
%!error<CompactClassificationGAM.edge: Name-Value arguments must be in pairs.> ...
%! edge (CM, x2, y2, 'Weights')
%!error<CompactClassificationGAM.loss: unsupported Loss function.> ...
%! loss (CM, x2, y2, 'LossFun', 'nonsense')

## A fitted model survives savemodel and loadmodel: the properties come
## back as they were and it predicts the same.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = compact (fitcgam (meas(inds,:), species(inds)));
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'CompactClassificationGAM');
%! assert_equal (M2.PredictorNames, Mdl.PredictorNames);
%! assert_equal (class (M2.ScoreTransform), class (Mdl.ScoreTransform));
%! assert_equal (predict (M2, meas(1:5,:)), predict (Mdl, meas(1:5,:)));

%!error<CompactClassificationGAM: the number of rows and columns in 'Cost' must correspond to selected classes in Y.> ...
%! Mdl = fitcgam ([1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1], [0; 0; 1; 1]);
%! CMdl = compact (Mdl);
%! CMdl.Cost = 1:4;
