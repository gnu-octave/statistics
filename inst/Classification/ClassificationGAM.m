## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
## Copyright (C) 2024-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{obj} =} ClassificationGAM (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} ClassificationGAM (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{ClassificationGAM} class object containing a generalized
## additive classification model.
##
## The @qcode{ClassificationGAM} class implements a gradient boosting algorithm
## for classification, using spline fitting as the weak learner. This approach
## allows the model to capture non-linear relationships between predictors and
## the binary response variable.
##
## @code{@var{obj} = ClassificationGAM (@var{X}, @var{Y})} returns a
## ClassificationGAM object, with @var{X} as the predictor data and @var{Y}
## containing the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of predictor data where rows
## correspond to observations and columns correspond to features or variables.
## @item
## @code{Y} is @math{Nx1} numeric vector containing binary class labels,
## typically 0 or 1.
## @end itemize
##
## @code{@var{obj} = ClassificationGAM (@dots{}, @var{name}, @var{value})}
## returns a ClassificationGAM object with parameters specified by
## @qcode{Name-Value} pair arguments. Type @code{help fitcgam} for more info.
##
## A @qcode{ClassificationGAM} object, @var{obj}, stores the labeled training
## data and various parameters for the Generalized Additive Model (GAM) for
## classification, which can be accessed in the following fields:
##
## @multitable @columnfractions 0.23 0.02 0.75
## @headitem @var{Field} @tab @tab @var{Description}
##
## @item @qcode{X} @tab @tab Predictor data, specified as a
## numeric matrix. Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{Y} @tab @tab Class labels, specified as numeric vector
## of 0's and 1's. Each value in @var{Y} is the observed class
## label for the corresponding row in @var{X}.
##
## @item @qcode{BaseModel} @tab @tab A structure containing the parameters
## of the base model without any interaction terms. The base model represents
## the generalized additive model (GAM) with only the main effects (predictor
## terms) included.
##
## @item @qcode{ModelwInt} @tab @tab A structure containing the parameters
## of the model that includes interaction terms. This model extends the base
## model by adding interaction terms between predictors, as specified by the
## @qcode{Interactions} property.
##
## @item @qcode{IntMatrix} @tab @tab A logical matrix or a matrix of
## column indices that describes the interaction terms applied to the predictor
## data.
##
## @item @qcode{NumObservations} @tab @tab Number of observations used in
## training the ClassificationGAM model, specified as a positive integer scalar.
## This number can be less than the number of rows in the training data because
## rows containing @qcode{NaN} values are not part of the fit.
##
## @item @qcode{RowsUsed} @tab @tab Rows of the original training data
## used in fitting the ClassificationGAM model, specified as a numerical vector.
## If you want to use this vector for indexing the training data in @var{X}, you
## have to convert it to a logical vector, i.e
## @qcode{X = obj.X(logical (obj.RowsUsed), :);}
##
## @item @qcode{NumPredictors} @tab @tab The number of predictors
## (variables) in @var{X}.
##
## @item @qcode{PredictorNames} @tab @tab Predictor variable names,
## specified as a cell array of character vectors. The variable names are in
## the same order in which they appear in the training data @var{X}.
##
## @item @qcode{ResponseName} @tab @tab Response variable name, specified
## as a character vector.
##
## @item @qcode{ClassNames} @tab @tab Names of the classes in the training
## data @var{Y} with duplicates removed, specified as a cell array of character
## vectors.
##
## @item @qcode{Cost} @tab @tab Cost of the misclassification of a point,
## specified as a square matrix. @qcode{Cost(i,j)} is the cost of classifying a
## point into class @qcode{j} if its true class is @qcode{i} (that is, the rows
## correspond to the true class and the columns correspond to the predicted
## class).  The order of the rows and columns in @qcode{Cost} corresponds to the
## order of the classes in @qcode{ClassNames}.  The number of rows and columns
## in @qcode{Cost} is the number of unique classes in the response.  By default,
## @qcode{Cost(i,j) = 1} if @qcode{i != j}, and @qcode{Cost(i,j) = 0} if
## @qcode{i = j}.  In other words, the cost is 0 for correct classification and
## 1 for incorrect classification.
##
## @item @qcode{Formula} @tab @tab A model specification given as a string
## in the form @qcode{"Y ~ terms"} where @qcode{Y} represents the response
## variable and @qcode{terms} the predictor variables. The formula can be used
## to specify a subset of variables for training model. For example:
## @qcode{"Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3"} specifies four linear terms
## for the first four columns of for predictor data, and @qcode{x1:x2} and
## @qcode{x2:x3} specify the two interaction terms for 1st-2nd and 3rd-4th
## columns respectively.  Only these terms will be used for training the model,
## but @var{X} must have at least as many columns as referenced in the formula.
## If Predictor Variable names have been defined, then the terms in the formula
## must reference to those.  When @qcode{"formula"} is specified, all terms used
## for training the model are referenced in the @qcode{IntMatrix} field of the
## @var{obj} class object as a matrix containing the column indexes for each
## term including both the predictors and the interactions used.
##
## @item @qcode{Interactions} @tab @tab A logical matrix, a positive integer
## scalar, or the string @qcode{"all"} for defining the interactions between
## predictor variables.  When given a logical matrix, it must have the same
## number of columns as @var{X} and each row corresponds to a different
## interaction term combining the predictors indexed as @qcode{true}.  Each
## interaction term is appended as a column vector after the available predictor
## column in @var{X}.  When @qcode{"all"} is defined, then all possible
## combinations of interactions are appended in @var{X} before training.  At the
## moment, parsing a positive integer has the same effect as the @qcode{"all"}
## option.  When @qcode{"interactions"} is specified, only the interaction terms
## appended to @var{X} are referenced in the @qcode{IntMatrix} field of the
## @var{obj} class object.
##
## @item @qcode{Knots} @tab @tab A scalar or a row vector with the same
## columns as @var{X}.  It defines the knots for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 5, hence expanded to @qcode{ones (1, columns (X)) * 5}.  You can
## parse a row vector with different number of knots for each predictor
## variable to be fitted with, although not recommended.
##
## @item @qcode{Order} @tab @tab A scalar or a row vector with the same
## columns as @var{X}.  It defines the order of the polynomial when training the
## GAM.  As a scalar, it is expanded to a row vector.  The default values is 3,
## hence expanded to @qcode{ones (1, columns (X)) * 3}.  You can parse a row
## vector with different number of polynomial order for each predictor variable
## to be fitted with, although not recommended.
##
## @item @qcode{DoF} @tab @tab A scalar or a row vector with the same
## columns as @var{X}.  It defines the degrees of freedom for fitting a
## polynomial when training the GAM.  As a scalar, it is expanded to a row
## vector.  The default value is 8, hence expanded to
## @qcode{ones (1, columns (X)) * 8}. You can parse a row vector with different
## degrees of freedom for each predictor variable to be fitted with,
## although not recommended.
##
## @end multitable
##
## You can parse either a @qcode{"Formula"} or an @qcode{"Interactions"}
## optional parameter.  Parsing both parameters will result an error.
## Accordingly, you can only pass up to two parameters among @qcode{"Knots"},
## @qcode{"Order"}, and @qcode{"DoF"} to define the required polynomial for
## training the GAM model.
##
## @seealso{fitcgam, CompactClassificationGAM, ClassificationPartitionedModel}
## @end deftypefn

  properties (Access = public)

    X = [];                 # Predictor data
    Y = [];                 # Class labels

    NumObservations = [];   # Number of observations in training dataset
    RowsUsed        = [];   # Rows used in fitting
    NumPredictors   = [];   # Number of predictors
    PredictorNames  = [];   # Predictor variable names
    ResponseName    = [];   # Response variable name
    ClassNames      = [];   # Names of classes in Y
    Prior           = [];   # Prior probability for each class
    Cost            = [];   # Cost of Misclassification

    ScoreTransform  = [];   # Transformation for classification scores

    Formula         = [];   # Formula for GAM model
    Interactions    = [];   # Number or matrix of interaction terms

    Knots           = [];   # Knots of spline fitting
    Order           = [];   # Order of spline fitting
    DoF             = [];   # Degrees of freedom for fitting spline

    LearningRate    = [];   # Learning rate for training
    NumIterations   = [];   # Max number of iterations for training

    BaseModel = [];         # Base model parameters (no interactions)
    ModelwInt = [];         # Model parameters with interactions
    IntMatrix = [];         # Interactions matrix applied to predictor data

  endproperties

  methods (Access = public)

    ## Class object constructor
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
      DoF            = ones (1, ndims_X) * 8;
      Order          = ones (1, ndims_X) * 3;
      Knots          = ones (1, ndims_X) * 5;
      LearningRate   = 0.1;
      NumIterations  = 100;
      Cost           = [];
      this.ScoreTransform  = 'none';

      ## Number of parameters for Knots, DoF, Order (maximum 2 allowed)
      KOD = 0;
      ## Number of parameters for Formula, Interactions (maximum 1 allowed)
      F_I = 0;

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("ClassificationGAM: 'PredictorNames'", ...
                             " must be supplied as a cellstring array."));
            elseif (numel (PredictorNames) != columns (X))
              error (strcat ("ClassificationGAM: 'PredictorNames'", ...
                             " must equal the number of columns in X."));
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("ClassificationGAM: 'ResponseName'", ...
                             " must be a character vector."));
            endif

          case "classnames"
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
                                   ClassNames, "UniformOutput", false))))
                error (strcat ("ClassificationGAM: not all 'ClassNames'", ...
                               " are present in Y."));
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, "UniformOutput", false))))
                error (strcat ("ClassificationGAM: not all 'ClassNames'", ...
                               " are present in Y."));
              endif
            endif

          case "cost"
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error (strcat ("ClassificationGAM: 'Cost' must be", ...
                             " a numeric square matrix."));
            endif

          case "scoretransform"
            name = "ClassificationGAM";
            this.ScoreTransform = parseScoreTransform (varargin{2}, name);

          case "formula"
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

          case "interactions"
            if (F_I < 1)
              tmp = varargin{2};
              if (isnumeric (tmp) && isscalar (tmp)
                                  && tmp == fix (tmp) && tmp >= 0)
                Interactions = tmp;
              elseif (islogical (tmp))
                Interactions = tmp;
              elseif (ischar (tmp) && strcmpi (tmp, "all"))
                Interactions = tmp;
              else
                error ("ClassificationGAM: invalid 'Interactions' parameter.");
              endif
              F_I += 1;
            else
              error ("ClassificationGAM: 'Formula' has already been defined.");
            endif

          case "knots"
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

          case "order"
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

          case "dof"
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

          case "learningrate"
            LearningRate = varargin{2};
            if (LearningRate > 1 || LearningRate <= 0)
              error (strcat ("ClassificationGAM: 'LearningRate'", ...
                             " must be between 0 and 1."));
            endif

          case "numiterations"
            NumIterations = varargin{2};
            if (! isnumeric (NumIterations) || NumIterations <= 0)
              error (strcat ("ClassificationGAM: 'NumIterations'", ...
                             " must be a positive integer value."));
            endif

          otherwise
            error (strcat ("ClassificationGAM: invalid parameter", ...
                           " name in optional pair arguments."));

        endswitch
        varargin (1:2) = [];
      endwhile

      ## Generate default predictors and response variable names (if necessary)
      if (isempty (PredictorNames))
        for i = 1:columns (X)
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = "Y";
      endif

      ## Assign predictors and response variable names
      this.PredictorNames = PredictorNames;
      this.ResponseName   = ResponseName;

      ## Handle class names
      if (! isempty (ClassNames))
        if (iscellstr (ClassNames))
          ru = find (! ismember (gnY, ClassNames));
        else
          ru = find (! ismember (glY, ClassNames));
        endif
        for i = 1:numel (ru)
          gY(gY == ru(i)) = NaN;
        endfor
      endif

      ## Remove nans from X and Y
      RowsUsed  = ! logical (sum (isnan ([X, gY]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

      ## Renew groups in Y
      [gY, gnY, glY] = grp2idx (Y);
      this.ClassNames = gnY;

      ## Check that we are dealing only with binary classification
      if (numel (gnY) > 2)
        error ("ClassificationGAM: can only be used for binary classification.");
      endif

      ## Force Y into numeric
      if (! isnumeric (Y))
        Y = gY - 1;
      endif

      this.NumObservations = rows (X);
      this.RowsUsed = cast (RowsUsed, "double");

      ## Assign the number of original predictors to the ClassificationGAM object
      this.NumPredictors = ndims_X;

      if (isempty (Cost))
        this.Cost = cast (! eye (numel (gnY)), "double");
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (strcat ("ClassificationGAM: the number of rows", ...
                         " and columns in 'Cost' must correspond", ...
                         " to selected classes in Y."));
        endif
        this.Cost = Cost;
      endif

      ## Assign remaining optional parameters
      this.Formula       = Formula;
      this.Interactions  = Interactions;
      this.Knots         = Knots;
      this.Order         = Order;
      this.DoF           = DoF;
      this.LearningRate  = LearningRate;
      this.NumIterations = NumIterations;

      ## Fit the basic model
      Inter = mean (Y);
      [iter, param, res, RSS, intercept] = this.fitGAM (X, Y, Inter, Knots, ...
                                                        Order, LearningRate, ...
                                                        NumIterations);
      this.BaseModel.Intercept  = intercept;
      this.BaseModel.Parameters = param;
      this.BaseModel.Iterations = iter;
      this.BaseModel.Residuals  = res;
      this.BaseModel.RSS        = RSS;

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
        [iter, param, res, RSS, intercept] = this.fitGAM (X, Y, Inter, Knots, ...
                                                          Order, LearningRate, ...
                                                          NumIterations);
        this.ModelwInt.Intercept  = intercept;
        this.ModelwInt.Parameters = param;
        this.ModelwInt.Iterations = iter;
        this.ModelwInt.Residuals  = res;
        this.ModelwInt.RSS        = RSS;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationGAM} {@var{label} =} predict (@dots{}, @qcode{'IncludeInteractions'}, @var{includeInteractions})
    ## @deftypefnx {ClassificationGAM} {[@var{label}, @var{score}] =} predict (@dots{})
    ##
    ## Predict labels for new data using the Generalized Additive Model (GAM)
    ## stored in a ClassificationGAM object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the predicted
    ## labels for the data in @var{X} based on the model stored in the
    ## ClassificationGAM object, @var{obj}.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC}, 'IncludeInteractions',
    ## @var{includeInteractions})} allows you to specify whether interaction
    ## terms should be included when making predictions.
    ##
    ## @code{[@var{label}, @var{score}] = predict (@dots{})} also returns
    ## @var{score}, which contains the predicted class scores or posterior
    ## probabilities for each observation.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationGAM} class object.
    ## @item
    ## @var{XC} must be an @math{MxP} numeric matrix where each row is an
    ## observation and each column corresponds to a predictor variable.
    ## @item
    ## @var{includeInteractions} is a 'true' or 'false' indicating whether to
    ## include interaction terms in the predictions.
    ## @end itemize
    ##
    ## @seealso{fitcgam, ClassificationGAM}
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
      XC        = XC (notnansf, :);

      ## Default values for Name-Value Pairs
      incInt = ! isempty (this.IntMatrix);
      Cost = this.Cost;

      ## Parse optional arguments
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "includeinteractions"
            tmpInt = varargin{2};
            if (! islogical (tmpInt) || (tmpInt != 0 && tmpInt != 1))
              error (strcat ("ClassificationGAM.predict:", ...
                             " includeinteractions must be a logical value."));
            endif
            ## Check model for interactions
            if (tmpInt && isempty (this.IntMatrix))
              error (strcat ("ClassificationGAM.predict: trained model", ...
                             " does not include any interactions."));
            endif
            incInt = tmpInt;

          otherwise
            error (strcat ("ClassificationGAM.predict: invalid NAME in", ...
                           " optional pairs of arguments."));
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Choose whether interactions must be included
      if (incInt)
        if (! isempty (this.Interactions))
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
      labels = this.ClassNames (minIdx);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationGAM} {@var{CVMdl} =} crossval (@dots{}, @var{Name}, @var{Value})
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
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"KFold"} @tab @tab Specify the number of folds to use in
    ## k-fold cross-validation.  @code{"KFold", @var{k}}, where @var{k} is an
    ## integer greater than 1.
    ##
    ## @item @qcode{"Holdout"} @tab @tab Specify the fraction of the data to
    ## hold out for testing.  @code{"Holdout", @var{p}}, where @var{p} is a
    ## scalar in the range @math{(0,1)}.
    ##
    ## @item @qcode{"Leaveout"} @tab @tab Specify whether to perform
    ## leave-one-out cross-validation.  @code{"Leaveout", @var{Value}}, where
    ## @var{Value} is 'on' or 'off'.
    ##
    ## @item @qcode{"CVPartition"} @tab @tab Specify a @qcode{cvpartition}
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
        varargin (1:2) = [];
      endwhile

      ## Determine the cross-validation method to use
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (this.Y, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (this.NumObservations, 'LeaveOut');
      else
        partition = cvpartition (this.Y, 'KFold', numFolds);
      endif

      ## Create a cross-validated model object
      CVMdl = ClassificationPartitionedModel (this, partition);
    endfunction

  endmethods

  ## Helper functions
  methods (Access = private)

    ## Determine interactions from Interactions optional parameter
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
        ## Get all combinations except all zeros
        allMat = flip (fullfact(p)([2:end],:), 2);
        ## Only keep interaction terms
        iterms = find (sum (allMat, 2) != 1);
        intMat = allMat(iterms);
      elseif (strcmpi (this.Interactions, "all"))
        p = this.NumPredictors;
        ## Calculate all p*(p-1)/2 interaction terms
        allMat = flip (fullfact(p)([2:end],:), 2);
        ## Only keep interaction terms
        iterms = find (sum (allMat, 2) != 1);
        intMat = allMat(iterms);
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
            error (strcat (["ClassificationGAM: some predictors"], ...
                           [" have not been identified."]));
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
      ## Initialize variables
      [n_samples, n_features] = size (X);
      RSS = zeros (1, n_features);

      ## Initialize model predictions with the intercept (log-odds)
      p = Inter;
      intercept = log (p / (1 - p));
      f = intercept * ones (n_samples, 1);

      ## Start boosting iterations
      for iter = 1:num_iterations
        ## Compute the gradient
        y_pred = 1 ./ (1 + exp (-f));  ## Sigmoid function
        gradient = Y - y_pred;         ## Negative gradient of log-loss

        ## Initialize a variable to store predictions for this iteration
        f_new = zeros (n_samples, 1);

        for j = 1:n_features
          ## Fit a spline to the gradient for feature X_j
          spline_model = splinefit (X(:, j), gradient, Knots(j), ...
                                             "order", Order(j));

          ## Predict using the fitted spline
          spline_pred = ppval (spline_model, X(:, j));

          ## Store the spline model parameters
          param(j) = spline_model;

          ## Update the model predictions
          f_new = f_new + learning_rate * spline_pred;
        endfor

        ## Update the overall model predictions
        f = f + f_new ;
      endfor

      ## Final residuals and RSS calculation
      res = Y  - 1 ./ (1 + exp (-f));
      RSS = sum (res .^ 2);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationGAM} {@var{CVMdl} =} compact (@var{obj})
    ##
    ## Create a CompactClassificationGAM object.
    ##
    ## @code{@var{CVMdl} = compact (@var{obj})} creates a compact version of the
    ## ClassificationGAM object, @var{obj}.
    ##
    ## @seealso{fitcdiscr, ClassificationGAM, CompactClassificationGAM}
    ## @end deftypefn

    function CVMdl = compact (this)
      ## Create a compact model
      CVMdl = CompactClassificationGAM (this);
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
    ## @seealso{loadmodel, fitcgam, ClassificationGAM, cvpartition,
    ## ClassificationPartitionedModel}
    ## @end deftypefn

    function savemodel (this, fname)
      ## Generate variable for class name
      classdef_name = "ClassificationGAM";

      ## Create variables from model properties
      X = this.X;
      Y = this.Y;
      NumObservations = this.NumObservations;
      RowsUsed        = this.RowsUsed;
      NumPredictors   = this.NumPredictors;
      PredictorNames  = this.PredictorNames;
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

      ## Save classdef name and all model properties as individual variables
      save ("-binary", fname, "classdef_name", "X", "Y", "NumObservations", ...
            "RowsUsed", "NumPredictors", "PredictorNames", "ResponseName", ...
            "ClassNames", "Prior", "Cost", "ScoreTransform", "Formula", ...
            "Interactions", "Knots", "Order", "DoF", "BaseModel", ...
            "ModelwInt", "IntMatrix");
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationGAM object
      mdl = ClassificationGAM (1, 1);

      ## Check that fieldnames in DATA match properties in ClassificationGAM
      names = fieldnames (data);
      props = fieldnames (mdl);
      if (! isequal (sort (names), sort (props)))
        error ("ClassificationGAM.load_model: invalid model in '%s'.", filename)
      endif

      ## Copy data into object
      for i = 1:numel (props)
        mdl.(props{i}) = data.(props{i});
      endfor
    endfunction

  endmethods

endclassdef

## Helper function
function scores = predict_val (params, XC, intercept)
  [nsample, ndims_X] = size (XC);
  ypred = ones (nsample, 1) * intercept;

  ## Add the remaining terms
  for j = 1:ndims_X
    ypred = ypred + ppval (params(j), XC (:,j));
  endfor

  ## Apply the sigmoid function to get probabilities
  pos_prob = 1 ./ (1 + exp (-ypred));
  neg_prob = 1 - pos_prob;

  scores = [neg_prob, pos_prob];
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
%! obj = fitcgam (X, Y, "Interactions", "all")
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
%! a = ClassificationGAM (x, y, "PredictorNames", PredictorNames);
%! assert (class (a), "ClassificationGAM");
%! assert ({a.X, a.Y, a.NumObservations}, {x, y, 4})
%! assert ({a.NumPredictors, a.ResponseName}, {3, "Y"})
%! assert (a.ClassNames, {'0'; '1'})
%! assert (a.PredictorNames, PredictorNames)
%! assert (a.BaseModel.Intercept, 0)
%!test
%! load fisheriris
%! inds = strcmp (species,'versicolor') | strcmp (species,'virginica');
%! X = meas(inds, :);
%! Y = species(inds, :)';
%! Y = strcmp (Y, 'virginica')';
%! a = ClassificationGAM (X, Y, 'Formula', 'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3');
%! assert (class (a), "ClassificationGAM");
%! assert ({a.X, a.Y, a.NumObservations}, {X, Y, 100})
%! assert ({a.NumPredictors, a.ResponseName}, {4, "Y"})
%! assert (a.ClassNames, {'0'; '1'})
%! assert (a.Formula, 'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3')
%! assert (a.PredictorNames, {'x1', 'x2', 'x3', 'x4'})
%! assert (a.ModelwInt.Intercept, 0)
%!test
%! X = [2, 3, 5; 4, 6, 8; 1, 2, 3; 7, 8, 9; 5, 4, 3];
%! Y = [0; 1; 0; 1; 1];
%! a = ClassificationGAM (X, Y, 'Knots', [4, 4, 4], 'Order', [3, 3, 3]);
%! assert (class (a), "ClassificationGAM");
%! assert ({a.X, a.Y, a.NumObservations}, {X, Y, 5})
%! assert ({a.NumPredictors, a.ResponseName}, {3, "Y"})
%! assert (a.ClassNames, {'0'; '1'})
%! assert (a.PredictorNames, {'x1', 'x2', 'x3'})
%! assert (a.Knots, [4, 4, 4])
%! assert (a.Order, [3, 3, 3])
%! assert (a.DoF, [7, 7, 7])
%! assert (a.BaseModel.Intercept, 0.4055, 1e-1)

## Test input validation for constructor
%!error<ClassificationGAM: too few input arguments.> ClassificationGAM ()
%!error<ClassificationGAM: too few input arguments.> ...
%! ClassificationGAM (ones(4, 1))
%!error<ClassificationGAM: number of rows in X and Y must be equal.> ...
%! ClassificationGAM (ones (4,2), ones (1,4))
%!error<ClassificationGAM: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), "PredictorNames", ["A"])
%!error<ClassificationGAM: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), "PredictorNames", "A")
%!error<ClassificationGAM: 'PredictorNames' must equal the number of columns in X.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), "PredictorNames", {"A", "B", "C"})
%!error<ClassificationGAM: 'ResponseName' must be a character vector.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), "ResponseName", {"Y"})
%!error<ClassificationGAM: 'ResponseName' must be a character vector.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), "ResponseName", 1)
%!error<ClassificationGAM: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationGAM (ones(10,2), ones (10,1), "ClassNames", @(x)x)
%!error<ClassificationGAM: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationGAM (ones(10,2), ones (10,1), "ClassNames", {1})
%!error<ClassificationGAM: not all 'ClassNames' are present in Y.> ...
%! ClassificationGAM (ones(10,2), ones (10,1), "ClassNames", [1, 2])
%!error<ClassificationGAM: not all 'ClassNames' are present in Y.> ...
%! ClassificationGAM (ones(5,2), ['a';'b';'a';'a';'b'], "ClassNames", ['a';'c'])
%!error<ClassificationGAM: not all 'ClassNames' are present in Y.> ...
%! ClassificationGAM (ones(5,2), {'a';'b';'a';'a';'b'}, "ClassNames", {'a','c'})
%!error<ClassificationGAM: not all 'ClassNames' are present in Y.> ...
%! ClassificationGAM (ones(10,2), logical (ones (10,1)), "ClassNames", [true, false])
%!error<ClassificationGAM: 'Cost' must be a numeric square matrix.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), "Cost", [1, 2])
%!error<ClassificationGAM: 'Cost' must be a numeric square matrix.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), "Cost", "string")
%!error<ClassificationGAM: 'Cost' must be a numeric square matrix.> ...
%! ClassificationGAM (ones (5,2), ones (5,1), "Cost", {eye(2)})

## Test predict method
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10];
%! y = [1; 0; 1; 0; 1];
%! a = ClassificationGAM (x, y, "interactions", "all");
%! l = {'0'; '0'; '0'; '0'; '0'};
%! s = [0.3760, 0.6240; 0.4259, 0.5741; 0.3760, 0.6240; ...
%!      0.4259, 0.5741; 0.3760, 0.6240];
%! [labels, scores] = predict (a, x);
%! assert (class (a), "ClassificationGAM");
%! assert ({a.X, a.Y, a.NumObservations}, {x, y, 5})
%! assert ({a.NumPredictors, a.ResponseName}, {2, "Y"})
%! assert (a.ClassNames, {'1'; '0'})
%! assert (a.PredictorNames, {'x1', 'x2'})
%! assert (a.ModelwInt.Intercept, 0.4055, 1e-1)
%! assert (labels, l)
%! assert (scores, s, 1e-1)
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [0; 0; 1; 1];
%! interactions = [false, true, false; true, false, true; false, true, false];
%! a = fitcgam (x, y, "learningrate", 0.2, "interactions", interactions);
%! [label, score] = predict (a, x, "includeinteractions", true);
%! l = {'0'; '0'; '1'; '1'};
%! s = [0.5106, 0.4894; 0.5135, 0.4865; 0.4864, 0.5136; 0.4847, 0.5153];
%! assert (class (a), "ClassificationGAM");
%! assert ({a.X, a.Y, a.NumObservations}, {x, y, 4})
%! assert ({a.NumPredictors, a.ResponseName}, {3, "Y"})
%! assert (a.ClassNames, {'0'; '1'})
%! assert (a.PredictorNames, {'x1', 'x2', 'x3'})
%! assert (a.ModelwInt.Intercept, 0)
%! assert (label, l)
%! assert (score, s, 1e-1)

## Test input validation for predict method
%!error<ClassificationGAM.predict: too few input arguments.> ...
%! predict (ClassificationGAM (ones (4,2), ones (4,1)))
%!error<ClassificationGAM.predict: XC is empty.> ...
%! predict (ClassificationGAM (ones (4,2), ones (4,1)), [])
%!error<ClassificationGAM.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (ClassificationGAM (ones (4,2), ones (4,1)), 1)

## Test crossval method
%!shared x, y, obj
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1; 4, 5, 6];
%! y = [0; 0; 1; 1; 0];
%! obj = fitcgam (x, y);
%!test
%! CVMdl = crossval (obj);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 4)
%! assert (class (CVMdl.Trained{1}), "CompactClassificationGAM")
%! assert (CVMdl.CrossValidatedModel, "ClassificationGAM")
%!test
%! CVMdl = crossval (obj, "KFold", 2);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 2)
%! assert (class (CVMdl.Trained{1}), "CompactClassificationGAM")
%! assert (CVMdl.CrossValidatedModel, "ClassificationGAM")
%!test
%! CVMdl = crossval (obj, "HoldOut", 0.2);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "CompactClassificationGAM")
%! assert (CVMdl.CrossValidatedModel, "ClassificationGAM")
%!test
%! partition = cvpartition (y, 'KFold', 3);
%! CVMdl = crossval (obj, 'cvPartition', partition);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert (CVMdl.KFold == 3)
%! assert (class (CVMdl.Trained{1}), "CompactClassificationGAM")
%! assert (CVMdl.CrossValidatedModel, "ClassificationGAM")

## Test input validation for crossval method
%!error<ClassificationGAM.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (obj, "kfold")
%!error<ClassificationGAM.crossval: specify only one of the optional Name-Value paired arguments.>...
%! crossval (obj, "kfold", 12, "holdout", 0.2)
%!error<ClassificationGAM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (obj, "kfold", 'a')
%!error<ClassificationGAM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (obj, "holdout", 2)
%!error<ClassificationGAM.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (obj, "leaveout", 1)
%!error<ClassificationGAM.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (obj, "cvpartition", 1)
