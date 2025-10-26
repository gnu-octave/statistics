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

classdef ClassificationDiscriminant
## -*- texinfo -*-
## @deftp {statistics} ClassificationDiscriminant
##
## Discriminant analysis classification
##
## The @code{ClassificationDiscriminant} class implements a linear discriminant
## analysis classifier object, which can predict responses for new data using
## the @code{predict} method.
##
## Discriminant analysis classification is a statistical method used to classify
## observations into predefined groups based on their characteristics.  It
## estimates the parameters of different distributions for each class and
## predicts the class of new observations by finding the one with the smallest
## misclassification cost.
##
## Create a @code{ClassificationDiscriminant} object by using the
## @code{fitcdiscr} function or the class constructor.
##
## @seealso{fitcdiscr}
## @end deftp

  properties (Access = public)
    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} X
    ##
    ## Predictor data
    ##
    ## A numerix matrix containing the unstandardized predictor data.  Each
    ## column of @var{X} represents one predictor (variable), and each row
    ## represents one observation.  This property is read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Y
    ##
    ## Class labels
    ##
    ## Specified as a logical or numeric vector, or cell array of character
    ## vectors.  Each value in @var{Y} is the observed class label for the
    ## corresponding row in @var{X}.  This property is read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer value specifying the number of observations in the
    ## training dataset used for training the ClassificationDiscriminant model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} RowsUsed
    ##
    ## Rows used for fitting
    ##
    ## A logical column vector with the same length as the observations in the
    ## original predictor data @var{X} specifying which rows have been used for
    ## fitting the ClassificationDiscriminant model.  This property is
    ## read-only.
    ##
    ## @end deftp
    RowsUsed        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer value specifying the number of predictors in the
    ## training dataset used for training the ClassificationDiscriminant model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumPredictors   = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which the appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames  = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName    = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} ClassNames
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
    ## @deftp {ClassificationDiscriminant} {property} Cost
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
    ## Add or change the @qcode{Cost} using dot notation as in:
    ## @itemize
    ## @item @qcode{@var{obj}.Cost = @var{costMatrix}}
    ## @end itemize
    ##
    ## @end deftp
    Cost            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Prior
    ##
    ## Prior probability for each class
    ##
    ## A numeric vector specifying the prior probabilities for each class.  The
    ## order of the elements in @qcode{Prior} corresponds to the order of the
    ## classes in @qcode{ClassNames}.
    ##
    ## Add or change the @qcode{Prior} using dot notation as in:
    ## @itemize
    ## @item @qcode{@var{obj}.Prior = @var{priorVector}}
    ## @end itemize
    ##
    ## @end deftp
    Prior           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} ScoreTransform
    ##
    ## Transformation function for classification scores
    ##
    ## Specified as a character vector representing a built-in function or as a
    ## function handle for transforming the classification scores.  The
    ## following built-in functions are supported:
    ##
    ## @itemize
    ## @item @qcode{'doublelogit'}
    ## @item @qcode{'invlogit'}
    ## @item @qcode{'ismax'}
    ## @item @qcode{'logit'}
    ## @item @qcode{'none'}
    ## @item @qcode{'identity'}
    ## @item @qcode{'sign'}
    ## @item @qcode{'symmetric'}
    ## @item @qcode{'symmetricismax'}
    ## @item @qcode{'symmetriclogit'}
    ## @end itemize
    ##
    ## Add or change the @qcode{ScoreTransform} using dot notation as in:
    ## @itemize
    ## @item @qcode{@var{obj}.ScoreTransform = 'function_name'}
    ## @item @qcode{@var{obj}.ScoreTransform = @function_handle}
    ## @end itemize
    ##
    ## @end deftp
    ScoreTransform  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Sigma
    ##
    ## Within-class covariance
    ##
    ## A numeric array specifying the within-class covariance.  For linear
    ## discriminant type (currently supported) this is a @math{PxP} matrix,
    ## where @math{P} is the number of predictors in @var{X}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Sigma           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Mu
    ##
    ## Class means
    ##
    ## A @math {KxP} numeric matrix specifying the mean of the multivariate
    ## normal distribution of each corresponding class, where @math{K} is the
    ## number of classes and @math{P} is the number of predictors in @var{X}.
    ## This property is read-only.
    ##
    ## @end deftp
    Mu              = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Coeffs
    ##
    ## Coefficient matrices
    ##
    ## A @math {KxK} structure containing the coeeficient matrices, where
    ## @math{K} is the number of classes.  If the @qcode{'FillCoeffs'} parameter
    ## was set to @qcode{'off'} in either the @code{fitcdiscr} function or the
    ## @code{ClassificationDiscriminant} constructor, then @qcode{Coeffs} is
    ## empty @qcode{([])}.  This property is read-only.
    ##
    ## @qcode{Coeffs(i,j)} contains the coefficients of the linear (currently
    ## supported) boundaries between the classes @code{i} and @code{j} in the
    ## following fields:
    ##
    ## @itemize
    ## @item @qcode{DiscrimType} - A character vector
    ## @item @qcode{Class1} - @qcode{@var{ClassNames}(i)}
    ## @item @qcode{Class2} - @qcode{@var{ClassNames}(j)}
    ## @item @qcode{Const} - A scalar
    ## @item @qcode{Linear} - A vector with length as the number of predictors.
    ## @end itemize
    ##
    ## @end deftp
    Coeffs          = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Delta
    ##
    ## Delta threshold
    ##
    ## A nonnegative scalar specifying the threshold for linear discriminant
    ## model.  Currently unimplemented and fixed to 0.  This property is
    ## read-only.
    ##
    ## @end deftp
    Delta           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} DiscrimType
    ##
    ## Discriminant type
    ##
    ## A character vector specifying the type discriminant model.  Currently
    ## only linear discriminant models are supported.  This property is
    ## read-only.
    ##
    ## @end deftp
    DiscrimType     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Gamma
    ##
    ## Gamma regularization parameter
    ##
    ## A scalar value ranging from 0 to 1, specifying the Gamma regularization
    ## parameter.  This property is read-only.
    ##
    ## @end deftp
    Gamma           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} MinGamma
    ##
    ## Minimum value for Gamma regularization parameter
    ##
    ## A scalar value ranging from 0 to 1, specifying the minimum value that the
    ## Gamma regularization parameter can have.  This property is read-only.
    ##
    ## @end deftp
    MinGamma        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} LogDetSigma
    ##
    ## A scalar value specifying the logarithm of the determinant of the
    ## within-class covariance matrix.  This property is read-only.
    ##
    ## @end deftp
    LogDetSigma     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} XCentered
    ##
    ## A matrix of the same size as @var{X} and the values in @var{X} with the
    ## corresponding class means subtracted.  This property is read-only.
    ##
    ## @end deftp
    XCentered       = [];

  endproperties

  methods (Hidden)

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
      fprintf ("\n  ClassificationDiscriminant\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      if (iscellstr (this.ClassNames))
        str = repmat ({"'%s'"}, 1, numel (this.ClassNames));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.ClassNames{:});
      elseif (ischar (this.ClassNames))
        str = repmat ({"'%s'"}, 1, rows (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, cellstr (this.ClassNames){:});
      else # single, double, logical
        str = repmat ({"%d"}, 1, numel (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, this.ClassNames);
      endif
      fprintf ("%+25s: %s\n", 'ClassNames', str);
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumPredictors', this.NumPredictors);
      fprintf ("%+25s: '%s'\n", 'DiscrimType', this.DiscrimType);
      fprintf ("%+25s: [%dx%d double]\n", 'Mu', size (this.Mu));
      fprintf ("%+25s: [%dx%d struct]\n\n", 'Coeffs', size (this.Sigma));
    endfunction

    ## Class specific subscripted reference
    function varargout = subsref (this, s)
      chain_s = s(2:end);
      s = s(1);
      switch (s.type)
        case '()'
          error (strcat ("Invalid () indexing for referencing values", ...
                         " in a ClassificationDiscriminant object."));
        case '{}'
          error (strcat ("Invalid {} indexing for referencing values", ...
                         " in a ClassificationDiscriminant object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("ClassificationDiscriminant.subsref: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("ClassificationDiscriminant.subsref:", ...
                           " unrecognized property: '%s'"), s.subs);
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
        error (strcat ("ClassificationDiscriminant.subsasgn:", ...
                       " chained subscripts not allowed."));
      endif
      switch s.type
        case '()'
          error (strcat ("Invalid () indexing for assigning values", ...
                         " to a ClassificationDiscriminant object."));
        case '{}'
          error (strcat ("Invalid {} indexing for assigning values", ...
                         " to a ClassificationDiscriminant object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("ClassificationDiscriminant.subsasgn: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          switch (s.subs)
            case 'Cost'
              this.Cost = setCost (this, val);
            case 'Prior'
              this.Prior = setPrior (this, val);
            case 'ScoreTransform'
              name = "ClassificationDiscriminant";
              this.ScoreTransform = parseScoreTransform (val, name);
            otherwise
              error (strcat ("ClassificationDiscriminant.subsasgn:", ...
                             " unrecognized or read-only property: '%s'"), ...
                             s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods (Access = public)
    ## -*- texinfo -*-
    ## @deftypefn  {statistics} {@var{obj} =} ClassificationDiscriminant (@var{X}, @var{Y})
    ## @deftypefnx {statistics} {@var{obj} =} ClassificationDiscriminant (@dots{}, @var{name}, @var{value})
    ##
    ## Create a @qcode{ClassificationDiscriminant} class object containing a
    ## discriminant analysis model.
    ##
    ## @code{@var{obj} = ClassificationDiscriminant (@var{X}, @var{Y})} returns
    ## a ClassificationDiscriminant object, with @var{X} as the predictor data
    ## and @var{Y} containing the class labels of observations in @var{X}.
    ##
    ## @itemize
    ## @item
    ## @code{X} must be a @math{NxP} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @var{X} will be used to train the discriminant model.
    ## @item
    ## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} can contain any type
    ## of categorical data. @var{Y} must have the same number of rows as
    ## @var{X}.
    ## @end itemize
    ##
    ## @code{@var{obj} = ClassificationDiscriminant (@dots{}, @var{name},
    ## @var{value})} returns a ClassificationDiscriminant object with parameters
    ## specified by the following@qcode{@var{name}, @var{value}} pair arguments:
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'PredictorNames'} @tab @tab A cell array of character
    ## vectors specifying the names of the predictors. The length of this array
    ## must match the number of columns in @var{X}.
    ##
    ## @item @qcode{'ResponseName'} @tab @tab A character vector specifying the
    ## name of the response variable.
    ##
    ## @item @qcode{'ClassNames'} @tab @tab Names of the classes in the class
    ## labels, @var{Y}, used for fitting the Discriminant model.
    ## @qcode{ClassNames} are of the same type as the class labels in @var{Y}.
    ##
    ## @item @qcode{'Prior'} @tab @tab A numeric vector specifying the prior
    ## probabilities for each class.  The order of the elements in @qcode{Prior}
    ## corresponds to the order of the classes in @qcode{ClassNames}.
    ## Alternatively, you can specify @qcode{"empirical"} to use the empirical
    ## class probabilities or @qcode{"uniform"} to assume equal class
    ## probabilities.
    ##
    ## @item @qcode{'Cost'} @tab @tab An @math{NxR} numeric matrix containing
    ## misclassification cost for the corresponding instances in @var{X}, where
    ## @math{R} is the number of unique categories in @var{Y}.  If an instance
    ## is correctly classified into its category the cost is calculated to be 1,
    ## otherwise 0. The cost matrix can be altered by using
    ## @code{@var{Mdl}.cost = somecost}.  By default, its value is
    ## @qcode{@var{cost} = ones (rows (X), numel (unique (Y)))}.
    ##
    ## @item @qcode{'DiscrimType'} @tab @tab A character vector or string scalar
    ## specifying the type of discriminant analysis to perform. The only
    ## supported value is @qcode{'linear'}.
    ##
    ## @item @qcode{'FillCoeffs'} @tab @tab A character vector or string scalar
    ## with values @qcode{'on'} or @qcode{'off'} specifying whether to fill the
    ## coefficients after fitting. If set to @qcode{"on"}, the coefficients are
    ## computed during model fitting, which can be useful for prediction.
    ##
    ## @item @qcode{'Gamma'} @tab @tab A numeric scalar specifying the
    ## regularization parameter for the covariance matrix. It adjusts the linear
    ## discriminant analysis to make the model more stable in the presence of
    ## multicollinearity or small sample sizes. A value of 0 corresponds to no
    ## regularization, while a value of 1 corresponds to
    ## a completely regularized model.
    ## @end multitable
    ##
    ## @seealso{fitcdiscr}
    ## @end deftypefn
    function this = ClassificationDiscriminant (X, Y, varargin)

      ## Check for appropriate number of input arguments
      if (nargin < 2)
        error ("ClassificationDiscriminant: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationDiscriminant: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      ## Validate X
      if (! isnumeric (X))
        error ("ClassificationDiscriminant: X must be a numeric matrix.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error (strcat ("ClassificationDiscriminant: number", ...
                       " of rows in X and Y must be equal."));
      endif

      ## Assign original X and Y data
      this.X = X;
      this.Y = Y;

      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);

      ## Set default values before parsing optional parameters
      ClassNames           = [];
      Cost                 = [];
      DiscrimType          = "linear";
      Gamma                = 0;
      Delta                = 0;
      NumPredictors        = [];
      PredictorNames       = {};
      ResponseName         = 'Y';
      Prior                = "empirical";
      FillCoeffs           = "on";
      this.ScoreTransform  = 'none';

      ## Parse optional parameters
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))

          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("ClassificationDiscriminant: 'PredictorNames'", ...
                             " must be supplied as a cellstring array."));
            elseif (numel (PredictorNames) != columns (X))
              error (strcat ("ClassificationDiscriminant: 'PredictorNames'", ...
                             " must equal the number of columns in X."));
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("ClassificationDiscriminant: 'ResponseName'", ...
                             " must be a character vector."));
            endif

          case "classnames"
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames) ||
                   islogical (ClassNames) || ischar (ClassNames)))
              error (strcat ("ClassificationDiscriminant: 'ClassNames'", ...
                             " must be a cell array of character vectors,", ...
                             " a logical vector, a numeric vector,", ...
                             " or a character array."));
            endif
            ## Check that all class names are available in gnY
            if (iscellstr (ClassNames) || ischar (ClassNames))
              ClassNames = cellstr (ClassNames);
              if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                   ClassNames, "UniformOutput", false))))
                error (strcat ("ClassificationDiscriminant: not all", ...
                               " 'ClassNames' are present in Y."));
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, "UniformOutput", false))))
                error (strcat ("ClassificationDiscriminant: not all", ...
                               " 'ClassNames' are present in Y."));
              endif
            endif

          case "prior"
            Prior = varargin{2};
            if (! ((isnumeric (Prior) && isvector (Prior)) ||
                  (strcmpi (Prior, "empirical") || strcmpi (Prior, "uniform"))))
              error (strcat ("ClassificationDiscriminant: 'Prior' must", ...
                             " be either a numeric or a character vector."));
            endif

          case "cost"
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error (strcat ("ClassificationDiscriminant: 'Cost'", ...
                             " must be a numeric square matrix."));
            endif

          case "scoretransform"
            name = "ClassificationDiscriminant";
            this.ScoreTransform = parseScoreTransform (varargin{2}, name);

          case "discrimtype"
            DiscrimType = tolower (varargin{2});
            if (! (strcmpi (DiscrimType, "linear")))
                error (strcat ("ClassificationDiscriminant: unsupported", ...
                               " discriminant type."));
            endif

          case "fillcoeffs"
            FillCoeffs = tolower (varargin{2});
            if (! any (strcmpi (FillCoeffs, {"on", "off"})))
              error (strcat ("ClassificationDiscriminant: 'FillCoeffs'", ...
                             " must be 'on' or 'off'."));
            endif

          case "gamma"
            Gamma = varargin{2};
            if (Gamma >= 1 || Gamma < 0)
              error (strcat ("ClassificationDiscriminant: 'Gamma'", ...
                             " must be between 0 and 1."));
            endif

          otherwise
            error (strcat ("ClassificationDiscriminant: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Generate default predictors and response variable names (if necessary)
      NumPredictors = columns (X);
      if (isempty (PredictorNames))
        for i = 1:NumPredictors
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = "Y";
      endif

      ## Assign predictors and response variable names
      this.NumPredictors  = NumPredictors;
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

      ## Remove missing values from X and Y
      RowsUsed  = ! logical (sum (isnan ([X, gY]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

      ## Renew groups in Y, get classes ordered, keep the same type
      [this.ClassNames, gnY, gY] = unique (Y);

      ## Check X contains valid data
      if (! (isnumeric (X) && isfinite (X)))
        error ("ClassificationDiscriminant: invalid values in X.");
      endif

      ## Assign the number of observations and their corresponding indices
      ## on the original data, which will be used for training the model,
      ## to the ClassificationNeuralNetwork object
      this.NumObservations = sum (RowsUsed);
      this.RowsUsed = RowsUsed;

      ## Handle Cost and Prior
      this = setCost (this, Cost, gnY);
      this = setPrior (this, Prior, gnY, gY);

      ## Assign DiscrimType
      this.DiscrimType = DiscrimType;
      this.Delta = Delta;
      this.Gamma = Gamma;

      num_classes = rows (this.ClassNames);
      num_features = columns (X);
      this.Mu = zeros (num_classes, num_features);
      for i = 1:num_classes
        this.Mu(i, :) = mean (X(gY == i, :), 1);
      endfor

      ## Center the predictors (XCentered)
      this.XCentered = zeros (size (X));
      for i = 1:rows (X)
        class_idx = gY(i);
        this.XCentered(i, :) = X(i, :) - this.Mu(class_idx, :);
      endfor

      ## Calculate Within-class covariance (Sigma)
      if (strcmp (this.DiscrimType, "linear"))
        this.Sigma = zeros (num_features);
        for i = 1:num_classes
          Xi = X(gY == i, :) - this.Mu(i, :);
          this.Sigma = this.Sigma + (Xi' * Xi);
        endfor
        this.Sigma = this.Sigma / (this.NumObservations - num_classes);

        ## Check for predictors wih zero within-class variance
        zwcv = find (diag (this.Sigma) == 0);
        if (! isempty (zwcv))
          msg = strcat ("ClassificationDiscriminant: Predictor", ...
                        " '%s' has zero within-class variance.");
          error (msg, PredictorNames{zwcv(1)});
        endif

        D = diag (diag (this.Sigma));

        ## FIX ME: MinGamma calculation is not same as Matlab.
        ## Instead of using (det (sigma) > 0) as a criterion (see code below),
        ## maybe Matlab is using some threshold.
        ## Also linear search might not be best here.

        ## Regularize Sigma
        this.Sigma = (this.Sigma * (1 - this.Gamma)) + (D * this.Gamma);
        this.MinGamma = 0;
        ## Calculate the MinGamma
        if (det (this.Sigma) <= 0)
          gamma = 0;
          step = 1e-15;
          sigma = this.Sigma;
          while (true)
            sigma = (sigma * (1 - gamma)) + (D * gamma);
            if (det (sigma) > 0)
              minGamma = gamma;
              break;
            endif
            gamma = gamma + step;
            if (gamma > 1)
              error (strcat ("ClassificationDiscriminant: failed to", ...
                             " find 'MinGamma' within reasonable range."));
            endif
          endwhile

          this.MinGamma = minGamma;
          if (this.Gamma < minGamma)
            this.Gamma = minGamma;
          endif
          this.Sigma = (this.Sigma * (1 - this.Gamma)) + (D * this.Gamma);
        endif
      endif

      ## Calculate log determinant of Sigma
      if (strcmp (this.DiscrimType, "linear"))
        this.LogDetSigma = log (det (this.Sigma));
      endif

      if (strcmpi (FillCoeffs, "on"))
        ## Calculate coefficients
        switch (this.DiscrimType)
          case "linear"
            this.Coeffs = struct();
            for i = 1:num_classes
              for j = 1:num_classes
                this.Coeffs(i, j).DiscrimType = "";
                this.Coeffs(i, j).Const = [];
                this.Coeffs(i, j).Linear = [];
                this.Coeffs(i, j).Class1 = this.ClassNames(i,:);
                this.Coeffs(i, j).Class2 = this.ClassNames(j,:);
                if (i != j)
                  A = (this.Mu(i, :) - this.Mu(j, :)) / this.Sigma;
                  K = log (this.Prior(i) ...
                      / this.Prior(j)) - 0.5 * (this.Mu(i, :) ...
                      / this.Sigma * this.Mu(i, :)') + 0.5 * (this.Mu(j, :) ...
                      / this.Sigma * this.Mu(j, :)');
                  this.Coeffs(i, j).DiscrimType = this.DiscrimType;
                  this.Coeffs(i, j).Linear = A';
                  this.Coeffs(i, j).Const = K;
                endif
              endfor
            endfor
        endswitch
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationDiscriminant} {[@var{label}, @var{score}, @var{cost}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the discriminant
    ## analysis model from a ClassificationDiscriminant object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## predictor data in @code{obj.X} and corresponding labels, @code{obj.Y},
    ## stored in the ClassificationDiscriminant model, @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationDiscriminant} class object.
    ## @item
    ## @var{XC} must be an @math{MxP} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the discriminant
    ## model in @var{obj}.
    ## @end itemize
    ##
    ## @code{[@var{label}, @var{score}, @var{cost}] = predict (@var{obj},
    ## @var{XC})} also returns @var{score}, which contains the predicted class
    ## scores or posterior probabilities for each instance of the corresponding
    ## unique classes, and @var{cost}, which is a matrix containing the expected
    ## cost of the classifications.
    ##
    ## The @var{score} matrix contains the posterior probabilities for each
    ## class, calculated using the multivariate normal probability density
    ## function and the prior probabilities of each class.  These scores are
    ## normalized to ensure they sum to 1 for each observation.
    ##
    ## The @var{cost} matrix contains the expected classification cost for each
    ## class, computed based on the posterior probabilities and the specified
    ## misclassification costs.
    ##
    ## @seealso{ClassificationDiscriminant, fitcdiscr}
    ## @end deftypefn
    function [label, score, cost] = predict (this, XC)
      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("ClassificationDiscriminant.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("ClassificationDiscriminant.predict: XC is empty.");
      elseif (columns (this.X) != columns (XC))
        error (strcat ("ClassificationDiscriminant.predict: XC must have ", ...
                       " the same number of predictors as the trained model."));
      endif

      ## Get training data and labels
      X = this.X(this.RowsUsed,:);
      Y = this.Y(this.RowsUsed,:);

      numObservations = rows (XC);
      numClasses = rows (this.ClassNames);
      score = zeros (numObservations, numClasses);
      cost = zeros (numObservations, numClasses);

      ## Calculate discriminant score (posterior probabilities)
      for i = 1:numClasses
        for j = 1:numObservations
          P_x_given_k = mvnpdf (XC(j, :), this.Mu(i, :), this.Sigma);
          score(j, i) = P_x_given_k * this.Prior(i);
        endfor
      endfor

      ## Normalize score to get posterior probabilities
      scoreSum = sum (score, 2);
      score = bsxfun (@rdivide, score, scoreSum);

      ## Handle numerical issues
      score(isnan (score)) = 0;

      ## Calculate expected classification cost
      for i = 1:numClasses
        cost(:, i) = sum (bsxfun (@times, score, this.Cost(:, i)'), 2);
      endfor

      ## Predict the class labels based on the minimum cost
      [~, minIdx] = min (cost, [], 2);
      label = this.ClassNames(minIdx,:);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationDiscriminant} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute loss for a trained ClassificationDiscriminant object.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} computes the loss,
    ## @var{L}, using the default loss function @qcode{'mincost'}.
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{ClassificationDiscriminant} object trained on
    ## @code{X} and @code{Y}.
    ## @item
    ## @code{X} must be a @math{NxP} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of Rows as @var{X}.
    ## @end itemize
    ##
    ## @code{@var{L} = loss (@dots{}, @var{name}, @var{value})} allows
    ## additional options specified by @var{name}-@var{value} pairs:
    ##
    ## @multitable @columnfractions 0.18 0.02 0.8
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"LossFun"} @tab @tab Specifies the loss function to use.
    ## Can be a function handle with four input arguments (C, S, W, Cost)
    ## which returns a scalar value or one of:
    ## 'binodeviance', 'classifcost', 'classiferror', 'exponential',
    ## 'hinge', 'logit','mincost', 'quadratic'.
    ## @itemize
    ## @item
    ## @code{C} is a logical matrix of size @math{NxK}, where @math{N} is the
    ## number of observations and @math{K} is the number of classes.
    ## The element @code{C(i,j)} is true if the class label of the i-th
    ## observation is equal to the j-th class.
    ## @item
    ## @code{S} is a numeric matrix of size @math{NxK}, where each element
    ## represents the classification score for the corresponding class.
    ## @item
    ## @code{W} is a numeric vector of length @math{N}, representing
    ## the observation weights.
    ## @item
    ## @code{Cost} is a @math{KxK} matrix representing the misclassification
    ## costs.
    ## @end itemize
    ##
    ## @item @qcode{"Weights"} @tab @tab Specifies observation weights, must be
    ## a numeric vector of length equal to the number of rows in X.
    ## Default is @code{ones (size (X, 1))}. loss normalizes the weights so that
    ## observation weights in each class sum to the prior probability of that
    ## class. When you supply Weights, loss computes the weighted
    ## classification loss.
    ##
    ## @end multitable
    ##
    ## @seealso{ClassificationDiscriminant}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationDiscriminant.loss: too few input arguments.");
      elseif (mod (nargin - 3, 2) != 0)
        error (strcat ("ClassificationDiscriminant.loss: name-value", ...
                       " arguments must be in pairs."));
      elseif (nargin > 7)
        error ("ClassificationDiscriminant.loss: too many input arguments.");
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationDiscriminant.loss: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat ("ClassificationDiscriminant.loss: X must have the", ...
                       " same number of predictors as the trained model."));
      endif

      ## Default values
      LossFun = 'mincost';
      Weights = [];

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error ("ClassificationDiscriminant.loss: Y must be of a valid type.");
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat ("ClassificationDiscriminant.loss: Y must", ...
                       " have the same number of rows as X."));
      endif

      ## Parse name-value arguments
      while (numel (varargin) > 0)
        Value = varargin{2};
        switch (tolower (varargin{1}))
          case 'lossfun'
            lf_opt = {"binodeviance", "classifcost", "classiferror", ...
                      "exponential", "hinge","logit", "mincost", "quadratic"};
            if (isa (Value, 'function_handle'))
              ## Check if the loss function is valid
              if (nargin (Value) != 4)
                error (strcat ("ClassificationDiscriminant.loss: custom", ...
                               " loss function must accept exactly four", ...
                               " input arguments."));
              endif
              try
                n = 1;
                K = 2;
                C_test = false (n, K);
                S_test = zeros (n, K);
                W_test = ones (n, 1);
                Cost_test = ones (K) - eye (K);
                test_output = Value (C_test, S_test, W_test, Cost_test);
                if (! isscalar (test_output))
                  error (strcat ("ClassificationDiscriminant.loss:", ...
                                 " custom loss function must return", ...
                                 " a scalar value."));
                endif
              catch
                error (strcat ("ClassificationDiscriminant.loss: custom", ...
                               " loss function is not valid or does not", ...
                               " produce correct output."));
              end_try_catch
              LossFun = Value;
            elseif (ischar (Value) && any (strcmpi (Value, lf_opt)))
              LossFun = Value;
            else
              error ("ClassificationDiscriminant.loss: invalid loss function.");
            endif

          case 'weights'
            if (isnumeric (Value) && isvector (Value))
              if (numel (Value) != size (X ,1))
                error (strcat ("ClassificationDiscriminant.loss: number", ...
                               " of 'Weights' must be equal to the", ...
                               " number of rows in X."));
              elseif (numel (Value) == size (X, 1))
                Weights = Value;
              endif
            else
              error ("ClassificationDiscriminant.loss: invalid 'Weights'.");
            endif

          otherwise
            error (strcat ("ClassificationDiscriminant.loss: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Check for missing values in X
      if (! isa (LossFun, 'function_handle'))
        lossfun = tolower (LossFun);
        if (! strcmp (lossfun, 'mincost') && ! strcmp (lossfun, 'classiferror')
            && ! strcmp (lossfun, 'classifcost') && any (isnan (X(:))))
          L = NaN;
          return;
        endif
      endif

      ## If Y is a char array convert it to a cell array of character vectors
      classes = this.ClassNames;
      if (ischar (Y) && ischar (classes))
        Y = cellstr (Y);
        classes = cellstr (classes);
      endif

      ## Check that Y is of the same type as ClassNames
      if (! strcmp (class (Y), class (classes)))
        error (strcat ("ClassificationDiscriminant.loss: Y must be", ...
                       " the same data type as the model's ClassNames."));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), classes)))
        error (strcat ("ClassificationDiscriminant.loss: Y must contain", ...
                       " only the classes in model's ClassNames."));
      endif

      ## Set default weights if not specified
      if (isempty (Weights))
        Weights = ones (size (X, 1), 1);
      endif

      ## Normalize Weights
      K = numel (classes);
      class_prior_probs = this.Prior;
      norm_weights = zeros (size (Weights));
      for i = 1:K
        class_idx = ismember (Y, classes(i));
        if (sum (Weights(class_idx)) > 0)
          norm_weights(class_idx) = ...
          Weights(class_idx) * class_prior_probs(i) / sum (Weights(class_idx));
        endif
      endfor
      Weights = norm_weights / sum (norm_weights);

      ## Number of observations
      n = size (X, 1);

      ## Predict classification scores
      [label, scores] = predict (this, X);

      ## C is vector of K-1 zeros, with 1 in the
      ## position corresponding to the true class
      C = false (n, K);
      for i = 1:n
        class_idx = find (ismember (classes, Y(i)));
        C(i, class_idx) = true;
      endfor
      Y_new = C';

      ## Compute the loss using custom loss function
      if (isa (LossFun, 'function_handle'))
        L = LossFun (C, scores, Weights, this.Cost);
        return;
      endif

      ## Compute the scalar classification score for each observation
      m_j = zeros (n, 1);
      for i = 1:n
        m_j(i) = scores(i,:) * Y_new(:,i);
      endfor

      ## Compute the loss
      switch (tolower (LossFun))
        case 'binodeviance'
          b = log (1 + exp (-2 * m_j));
          L = (Weights') * b;
        case 'hinge'
          h = max (0, 1 - m_j);
          L = (Weights') * h;
        case 'exponential'
          e = exp (-m_j);
          L = (Weights') * e;
        case 'logit'
          l = log (1 + exp (-m_j));
          L = (Weights') * l;
        case 'quadratic'
          q = (1 - m_j) .^ 2;
          L = (Weights') * q;
        case 'classiferror'
          L = 0;
          for i = 1:n
            L = L + Weights(i) * (! isequal (Y(i), label(i)));
          endfor
        case 'mincost'
          Cost = this.Cost;
          L = 0;
          for i = 1:n
            f_Xj = scores(i, :);
            gamma_jk = f_Xj * Cost;
            [~, min_cost_class] = min (gamma_jk);
            cj = Cost(find (ismember (classes, Y(i))), min_cost_class);
            L = L + Weights(i) * cj;
          endfor
        case 'classifcost'
          Cost = this.Cost;
          L = 0;
          for i = 1:n
            y_idx = find (ismember (classes, Y(i)));
            y_hat_idx = find (ismember (classes, label(i)));
            L = L + Weights(i) * Cost(y_idx, y_hat_idx);
          endfor
        otherwise
          error ("ClassificationDiscriminant.loss: invalid loss function.");
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationDiscriminant} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margins for discriminant analysis classifier.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns
    ## the classification margins for @var{obj} with data @var{X} and
    ## classification @var{Y}. @var{m} is a numeric vector of length size (X,1).
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{ClassificationDiscriminant} object trained on @code{X}
    ## and @code{Y}.
    ## @item
    ## @code{X} must be a @math{NxP} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of Rows as @var{X}.
    ## @end itemize
    ##
    ## The classification margin for each observation is the difference between
    ## the classification score for the true class and the maximal
    ## classification score for the false classes.
    ##
    ## @seealso{fitcdiscr, ClassificationDiscriminant}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationDiscriminant.margin: too few input arguments.");
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationDiscriminant.margin: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat ("ClassificationDiscriminant.margin: X must have the", ...
                       " same number of predictors as the trained model."));
      endif

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error ("ClassificationDiscriminant.margin: Y must be of a valid type.");
      endif

      ## Validate X
      valid_types = {'single', 'double'};
      if (! (any (strcmp (class (X), valid_types))))
        error ("ClassificationDiscriminant.margin: X must be of a valid type.");
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat ("ClassificationDiscriminant.margin: Y must", ...
                       " have the same number of rows as X."));
      endif

      ## If Y is a char array convert it to a cell array of character vectors
      classes = this.ClassNames;
      if (ischar (Y) && ischar (classes))
        Y = cellstr (Y);
        classes = cellstr (classes);
      endif

      ## Check that Y is of the same type as ClassNames
      if (! strcmp (class (Y), class (classes)))
        error (strcat ("ClassificationDiscriminant.margin: Y must be", ...
                       " the same data type as the model's ClassNames."));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), classes)))
        error (strcat ("ClassificationDiscriminant.margin: Y must", ...
                       " contain only the classes in model's ClassNames."));
      endif

      ## Number of Observations
      n = size (X, 1);

      ## Initialize the margin vector
      m = zeros (n, 1);

      ## Calculate the classification scores
      [~, scores] = predict (this, X);

      ## Loop over each observation to compute the margin
      for i = 1:n
        ## True class index
        true_class_idx = find (ismember (classes, Y(i)));

        ## Score for the true class
        true_class_score = scores(i, true_class_idx);

        ## Get the maximal score for the false classes
        scores(i, true_class_idx) = -Inf;              # Temporarily
        max_false_class_score = max (scores(i, :));
        if (max_false_class_score == -Inf)
          m = NaN;
          return;
        endif
        scores(i, true_class_idx) = true_class_score;  # Restore

        ## Calculate the margin
        m(i) = true_class_score - max_false_class_score;
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationDiscriminant} {@var{CVMdl} =} crossval (@dots{}, @var{Name}, @var{Value})
    ##
    ## Cross Validate a Discriminant classification object.
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
    ## @seealso{fitcdiscr, ClassificationDiscriminant, cvpartition,
    ## ClassificationPartitionedModel}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)
      ## Check input
      if (nargin < 1)
        error ("ClassificationDiscriminant.crossval: too few input arguments.");
      endif

      if (numel (varargin) == 1)
        error (strcat ("ClassificationDiscriminant.crossval: Name-Value", ...
                       " arguments must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("ClassificationDiscriminant.crossval: specify only", ...
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
              error (strcat ("ClassificationDiscriminant.crossval: 'KFold'", ...
                             " must be an integer value greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("ClassificationDiscriminant.crossval: 'Holdout'", ...
                             " must be a numeric value between 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat ("ClassificationDiscriminant.crossval:", ...
                             " 'Leaveout' must be either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("ClassificationDiscriminant.crossval:",...
                             " 'CVPartition' must be a 'cvpartition' object."));
            endif

          otherwise
            error (strcat ("ClassificationDiscriminant.crossval: invalid",...
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
        partition = cvpartition (numel (this.Y), 'LeaveOut');
      else
        partition = cvpartition (this.Y, 'KFold', numFolds);
      endif

      ## Create a cross-validated model object
      CVMdl = ClassificationPartitionedModel (this, partition);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{CVMdl} =} compact (@var{obj})
    ##
    ## Create a CompactClassificationDiscriminant object.
    ##
    ## @code{@var{CVMdl} = compact (@var{obj})} creates a compact version of the
    ## ClassificationDiscriminant object, @var{obj}.
    ##
    ## @seealso{fitcdiscr, ClassificationDiscriminant,
    ## CompactClassificationDiscriminant}
    ## @end deftypefn
    function CVMdl = compact (this)
      ## Create a compact model
      CVMdl = CompactClassificationDiscriminant (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationDiscriminant object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## ClassificationDiscriminant object into an Octave binary file, the name of
    ## which is specified in @var{filename}, along with an extra variable, which
    ## defines the type classification object these variables constitute.  Use
    ## @code{loadmodel} in order to load a classification object into Octave's
    ## workspace.
    ##
    ## @seealso{loadmodel, fitcdiscr, ClassificationDiscriminant}
    ## @end deftypefn
    function savemodel (this, fname)
      ## Generate variable for class name
      classdef_name = "ClassificationDiscriminant";

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
      Sigma           = this.Sigma;
      Mu              = this.Mu;
      Coeffs          = this.Coeffs;
      Delta           = this.Delta;
      DiscrimType     = this.DiscrimType;
      Gamma           = this.Gamma;
      MinGamma        = this.MinGamma;
      LogDetSigma     = this.LogDetSigma;
      XCentered       = this.XCentered;

      ## Save classdef name and all model properties as individual variables
      save ("-binary", fname, "classdef_name", "X", "Y", "NumObservations", ...
            "RowsUsed", "NumPredictors", "PredictorNames", "ResponseName", ...
            "ClassNames", "ScoreTransform", "Prior", "Cost", "Sigma", "Mu", ...
            "Coeffs", "Delta", "DiscrimType", "Gamma", "MinGamma", ...
            "LogDetSigma", "XCentered");
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationDiscriminant object
      mdl = ClassificationDiscriminant (1, 1);

      ## Check that fieldnames in DATA match properties in ClassificationDiscriminant
      names = fieldnames (data);
      props = fieldnames (mdl);
      if (! isequal (sort (names), sort (props)))
        msg = "ClassificationDiscriminant.load_model: invalid model in '%s'.";
        error (msg, filename);
      endif

      ## Copy data into object
      for i = 1:numel (props)
        mdl.(props{i}) = data.(props{i});
      endfor
    endfunction

  endmethods

  methods (Access = private)

    function this = setCost (this, Cost, gnY = [])
      if (isempty (gnY))
        [~, gnY, gY] = unique (this.Y(this.RowsUsed));
      endif
      if (isempty (Cost))
        this.Cost = cast (! eye (numel (gnY)), "double");
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (strcat ("ClassificationDiscriminant: the number", ...
                         " of rows and columns in 'Cost' must", ...
                         " correspond to selected classes in Y."));
        endif
        this.Cost = Cost;
      endif

    endfunction

    function this = setPrior (this, Prior, gnY = [], gY = [])
      if (isempty (gnY) || isempty (gY))
        [~, gnY, gY] = unique (this.Y(this.RowsUsed));
      endif
      ## Set prior
      if (strcmpi ("uniform", Prior))
        this.Prior = ones (size (gnY)) ./ numel (gnY);
      elseif (isempty (Prior) || strcmpi ("empirical", Prior))
        pr = [];
        for i = 1:numel (gnY)
          pr = [pr; sum(gY==i)];
        endfor
        this.Prior = pr ./ sum (pr);
      elseif (isnumeric (Prior))
        if (numel (gnY) != numel (Prior))
          error (strcat ("ClassificationDiscriminant: the elements", ...
                         " in 'Prior' do not correspond to the", ...
                         " selected classes in Y."));
        endif
        this.Prior = Prior ./ sum (Prior);
      endif

      ## Recalculate the Const field in the Coeffs structure
      if (! isempty (this.Coeffs))
        num_classes = rows (this.ClassNames);
        ## Calculate coefficients
        switch (this.DiscrimType)
          case "linear"
            for i = 1:num_classes
              for j = 1:num_classes
                if (i != j)
                  K = log (this.Prior(i) ...
                      / this.Prior(j)) - 0.5 * (this.Mu(i, :) ...
                      / this.Sigma * this.Mu(i, :)') + 0.5 * (this.Mu(j, :) ...
                      / this.Sigma * this.Mu(j, :)');
                  this.Coeffs(i, j).Const = K;
                endif
              endfor
            endfor
        endswitch
      endif
    endfunction

  endmethods

endclassdef

%!demo
%! ## Create discriminant classifier
%! ## Evaluate some model predictions on new data.
%!
%! load fisheriris
%! x = meas;
%! y = species;
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcdiscr (x, y);
%! [label, score, cost] = predict (obj, xc);

%!demo
%! load fisheriris
%! model = fitcdiscr (meas, species);
%! X = mean (meas);
%! Y = {'versicolor'};
%! ## Compute loss for discriminant model
%! L = loss (model, X, Y)

%!demo
%! load fisheriris
%! mdl = fitcdiscr (meas, species);
%! X = mean (meas);
%! Y = {'versicolor'};
%! ## Margin for discriminant model
%! m = margin (mdl, X, Y)

%!demo
%! load fisheriris
%! x = meas;
%! y = species;
%! obj = fitcdiscr (x, y, "gamma", 0.4);
%! ## Cross-validation for discriminant model
%! CVMdl = crossval (obj)

## Test constructor
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! PredictorNames = {'Sepal Length', 'Sepal Width', 'Petal Length', 'Petal Width'};
%! Mdl = ClassificationDiscriminant (x, y, "PredictorNames", PredictorNames);
%! sigma = [0.265008, 0.092721, 0.167514, 0.038401; ...
%!          0.092721, 0.115388, 0.055244, 0.032710; ...
%!          0.167514, 0.055244, 0.185188, 0.042665; ...
%!          0.038401, 0.032710, 0.042665, 0.041882];
%! mu = [5.0060, 3.4280, 1.4620, 0.2460; ...
%!       5.9360, 2.7700, 4.2600, 1.3260; ...
%!       6.5880, 2.9740, 5.5520, 2.0260];
%! xCentered = [ 9.4000e-02,  7.2000e-02, -6.2000e-02, -4.6000e-02; ...
%!              -1.0600e-01, -4.2800e-01, -6.2000e-02, -4.6000e-02; ...
%!              -3.0600e-01, -2.2800e-01, -1.6200e-01, -4.6000e-02];
%! assert (class (Mdl), "ClassificationDiscriminant");
%! assert ({Mdl.X, Mdl.Y, Mdl.NumObservations}, {x, y, 150})
%! assert ({Mdl.DiscrimType, Mdl.ResponseName}, {"linear", "Y"})
%! assert ({Mdl.Gamma, Mdl.MinGamma}, {0, 0}, 1e-15)
%! assert (Mdl.ClassNames, unique (species))
%! assert (Mdl.Sigma, sigma, 1e-6)
%! assert (Mdl.Mu, mu, 1e-14)
%! assert (Mdl.XCentered([1:3],:), xCentered, 1e-14)
%! assert (Mdl.LogDetSigma, -9.9585, 1e-4)
%! assert (Mdl.PredictorNames, PredictorNames)
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! Mdl = ClassificationDiscriminant (x, y, "Gamma", 0.5);
%! sigma = [0.265008, 0.046361, 0.083757, 0.019201; ...
%!          0.046361, 0.115388, 0.027622, 0.016355; ...
%!          0.083757, 0.027622, 0.185188, 0.021333; ...
%!          0.019201, 0.016355, 0.021333, 0.041882];
%! mu = [5.0060, 3.4280, 1.4620, 0.2460; ...
%!       5.9360, 2.7700, 4.2600, 1.3260; ...
%!       6.5880, 2.9740, 5.5520, 2.0260];
%! xCentered = [ 9.4000e-02,  7.2000e-02, -6.2000e-02, -4.6000e-02; ...
%!              -1.0600e-01, -4.2800e-01, -6.2000e-02, -4.6000e-02; ...
%!              -3.0600e-01, -2.2800e-01, -1.6200e-01, -4.6000e-02];
%! assert (class (Mdl), "ClassificationDiscriminant");
%! assert ({Mdl.X, Mdl.Y, Mdl.NumObservations}, {x, y, 150})
%! assert ({Mdl.DiscrimType, Mdl.ResponseName}, {"linear", "Y"})
%! assert ({Mdl.Gamma, Mdl.MinGamma}, {0.5, 0})
%! assert (Mdl.ClassNames, unique (species))
%! assert (Mdl.Sigma, sigma, 1e-6)
%! assert (Mdl.Mu, mu, 1e-14)
%! assert (Mdl.XCentered([1:3],:), xCentered, 1e-14)
%! assert (Mdl.LogDetSigma, -8.6884, 1e-4)

## Test input validation for constructor
%!shared X, Y, MODEL
%! X = rand (10,2);
%! Y = [ones(5,1);2*ones(5,1)];
%! MODEL = ClassificationDiscriminant (X, Y);
%!error<ClassificationDiscriminant: too few input arguments.> ClassificationDiscriminant ()
%!error<ClassificationDiscriminant: too few input arguments.> ...
%! ClassificationDiscriminant (ones(4, 1))
%!error<ClassificationDiscriminant: Name-Value arguments must be in pairs.> ...
%! ClassificationDiscriminant (X, Y, "prior")
%!error<ClassificationDiscriminant: number of rows in X and Y must be equal.> ...
%! ClassificationDiscriminant (ones (4,2), ones (1,4))
%!error<ClassificationDiscriminant: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationDiscriminant (X, Y, "PredictorNames", ["A"])
%!error<ClassificationDiscriminant: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationDiscriminant (X, Y, "PredictorNames", "A")
%!error<ClassificationDiscriminant: 'PredictorNames' must equal the number of columns in X.> ...
%! ClassificationDiscriminant (X, Y, "PredictorNames", {"A", "B", "C"})
%!error<ClassificationDiscriminant: 'ResponseName' must be a character vector.> ...
%! ClassificationDiscriminant (X, Y, "ResponseName", {"Y"})
%!error<ClassificationDiscriminant: 'ResponseName' must be a character vector.> ...
%! ClassificationDiscriminant (X, Y, "ResponseName", 1)
%!error<ClassificationDiscriminant: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationDiscriminant (X, Y, "ClassNames", @(x)x)
%!error<ClassificationDiscriminant: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationDiscriminant (X, Y, "ClassNames", {1})
%!error<ClassificationDiscriminant: not all 'ClassNames' are present in Y.> ...
%! ClassificationDiscriminant (X, ones (10,1), "ClassNames", [1, 2])
%!error<ClassificationDiscriminant: not all 'ClassNames' are present in Y.> ...
%! ClassificationDiscriminant ([1;2;3;4;5], ['a';'b';'a';'a';'b'], "ClassNames", ['a';'c'])
%!error<ClassificationDiscriminant: not all 'ClassNames' are present in Y.> ...
%! ClassificationDiscriminant ([1;2;3;4;5], {'a';'b';'a';'a';'b'}, "ClassNames", {'a','c'})
%!error<ClassificationDiscriminant: not all 'ClassNames' are present in Y.> ...
%! ClassificationDiscriminant (X, logical (ones (10,1)), "ClassNames", [true, false])
%!error<ClassificationDiscriminant: 'Prior' must be either a numeric or a character vector.> ...
%! ClassificationDiscriminant (X, Y, "Prior", {"1", "2"})
%!error<ClassificationDiscriminant: the elements in 'Prior' do not correspond to the selected classes in Y.> ...
%! ClassificationDiscriminant (X, ones (10,1), "Prior", [1 2])
%!error<ClassificationDiscriminant: 'Cost' must be a numeric square matrix.> ...
%! ClassificationDiscriminant (X, Y, "Cost", [1, 2])
%!error<ClassificationDiscriminant: 'Cost' must be a numeric square matrix.> ...
%! ClassificationDiscriminant (X, Y, "Cost", "string")
%!error<ClassificationDiscriminant: 'Cost' must be a numeric square matrix.> ...
%! ClassificationDiscriminant (X, Y, "Cost", {eye(2)})
%!error<ClassificationDiscriminant: the number of rows and columns in 'Cost' must correspond to selected classes in Y.> ...
%! ClassificationDiscriminant (X, Y, "Cost", ones (3))
%!error<ClassificationDiscriminant: Predictor 'x1' has zero within-class variance.> ...
%! ClassificationDiscriminant (ones (5,2), [1; 1; 2; 2; 2])
%!error<ClassificationDiscriminant: Predictor 'A' has zero within-class variance.> ...
%! ClassificationDiscriminant (ones (5,2), [1; 1; 2; 2; 2], "PredictorNames", {"A", "B"})
%!error<ClassificationDiscriminant: Predictor 'x2' has zero within-class variance.> ...
%! ClassificationDiscriminant ([1,2;2,2;3,2;4,2;5,2], ones (5, 1))
%!error<ClassificationDiscriminant: Predictor 'B' has zero within-class variance.> ...
%! ClassificationDiscriminant ([1,2;2,2;3,2;4,2;5,2], ones (5, 1), "PredictorNames", {"A", "B"})

## Test predict method
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! Mdl = fitcdiscr (meas, species, "Gamma", 0.5);
%! [label, score, cost] = predict (Mdl, [2, 2, 2, 2]);
%! assert (label, {'versicolor'})
%! assert (score, [0, 0.9999, 0.0001], 1e-4)
%! assert (cost, [1, 0.0001, 0.9999], 1e-4)
%! [label, score, cost] = predict (Mdl, [2.5, 2.5, 2.5, 2.5]);
%! assert (label, {'versicolor'})
%! assert (score, [0, 0.6368, 0.3632], 1e-4)
%! assert (cost, [1, 0.3632, 0.6368], 1e-4)
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! xc = [min(x); mean(x); max(x)];
%! Mdl = fitcdiscr (x, y);
%! [label, score, cost] = predict (Mdl, xc);
%! l = {'setosa'; 'versicolor'; 'virginica'};
%! s = [1, 0, 0; 0, 1, 0; 0, 0, 1];
%! c = [0, 1, 1; 1, 0, 1; 1, 1, 0];
%! assert (label, l)
%! assert (score, s, 1e-4)
%! assert (cost, c, 1e-4)

## Test input validation for predict method
%!error<ClassificationDiscriminant.predict: too few input arguments.> ...
%! predict (MODEL)
%!error<ClassificationDiscriminant.predict: XC is empty.> ...
%! predict (MODEL, [])
%!error<ClassificationDiscriminant.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (MODEL, 1)

## Test loss method
%!test
%! load fisheriris
%! model = fitcdiscr (meas, species);
%! x = mean (meas);
%! y = {'versicolor'};
%! L = loss (model, x, y);
%! assert (L, 0)
%!test
%! x = [1, 2; 3, 4; 5, 6];
%! y = {'A'; 'B'; 'A'};
%! model = fitcdiscr (x, y, "Gamma", 0.4);
%! x_test = [1, 6; 3, 3];
%! y_test = {'A'; 'B'};
%! L = loss (model, x_test, y_test);
%! assert (L, 0.3333, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = ['1'; '2'; '3'; '1'];
%! model = fitcdiscr (x, y, "gamma" , 0.5);
%! x_test = [3, 3];
%! y_test = ['1'];
%! L = loss (model, x_test, y_test, 'LossFun', 'quadratic');
%! assert (L, 0.2423, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = ['1'; '2'; '3'; '1'];
%! model = fitcdiscr (x, y, "gamma" , 0.5);
%! x_test = [3, 3; 5, 7];
%! y_test = ['1'; '2'];
%! L = loss (model, x_test, y_test, 'LossFun', 'classifcost');
%! assert (L, 0.3333, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = ['1'; '2'; '3'; '1'];
%! model = fitcdiscr (x, y, "gamma" , 0.5);
%! x_test = [3, 3; 5, 7];
%! y_test = ['1'; '2'];
%! L = loss (model, x_test, y_test, 'LossFun', 'hinge');
%! assert (L, 0.5886, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = ['1'; '2'; '3'; '1'];
%! model = fitcdiscr (x, y, "gamma" , 0.5);
%! x_test = [3, 3; 5, 7];
%! y_test = ['1'; '2'];
%! W = [1; 2];
%! L = loss (model, x_test, y_test, 'LossFun', 'logit', 'Weights', W);
%! assert (L, 0.5107, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6];
%! y = {'A'; 'B'; 'A'};
%! model = fitcdiscr (x, y, "gamma" , 0.5);
%! x_with_nan = [1, 2; NaN, 4];
%! y_test = {'A'; 'B'};
%! L = loss (model, x_with_nan, y_test);
%! assert (L, 0.3333, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6];
%! y = {'A'; 'B'; 'A'};
%! model = fitcdiscr (x, y);
%! x_with_nan = [1, 2; NaN, 4];
%! y_test = {'A'; 'B'};
%! L = loss (model, x_with_nan, y_test, 'LossFun', 'logit');
%! assert (isnan (L))
%!test
%! x = [1, 2; 3, 4; 5, 6];
%! y = {'A'; 'B'; 'A'};
%! model = fitcdiscr (x, y);
%! customLossFun = @(C, S, W, Cost) sum (W .* sum (abs (C - S), 2));
%! L = loss (model, x, y, 'LossFun', customLossFun);
%! assert (L, 0.8889, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6];
%! y = [1; 2; 1];
%! model = fitcdiscr (x, y);
%! L = loss (model, x, y, 'LossFun', 'classiferror');
%! assert (L, 0.3333, 1e-4)

## Test input validation for loss method
%!error<ClassificationDiscriminant.loss: too few input arguments.> ...
%! loss (MODEL)
%!error<ClassificationDiscriminant.loss: too few input arguments.> ...
%! loss (MODEL, ones (4,2))
%!error<ClassificationDiscriminant.loss: X is empty.> ...
%! loss (MODEL, [], zeros (2))
%!error<ClassificationDiscriminant.loss: X must have the same number of predictors as the trained model.> ...
%! loss (MODEL, 1, zeros (2))
%!error<ClassificationDiscriminant.loss: name-value arguments must be in pairs.> ...
%! loss (MODEL, ones (4,2), ones (4,1), 'LossFun')
%!error<ClassificationDiscriminant.loss: Y must have the same number of rows as X.> ...
%! loss (MODEL, ones (4,2), ones (3,1))
%!error<ClassificationDiscriminant.loss: invalid loss function.> ...
%! loss (MODEL, ones (4,2), ones (4,1), 'LossFun', 'a')
%!error<ClassificationDiscriminant.loss: invalid 'Weights'.> ...
%! loss (MODEL, ones (4,2), ones (4,1), 'Weights', 'w')

## Test margin method
%! load fisheriris
%! mdl = fitcdiscr (meas, species);
%! X = mean (meas);
%! Y = {'versicolor'};
%! m = margin (mdl, X, Y);
%! assert (m, 1, 1e-6)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = [1; 2; 1];
%! mdl = fitcdiscr (X, Y, "gamma", 0.5);
%! m = margin (mdl, X, Y);
%! assert (m, [0.3333; -0.3333; 0.3333], 1e-4)

## Test input validation for margin method
%!error<ClassificationDiscriminant.margin: too few input arguments.> ...
%! margin (MODEL)
%!error<ClassificationDiscriminant.margin: too few input arguments.> ...
%! margin (MODEL, ones (4,2))
%!error<ClassificationDiscriminant.margin: X is empty.> ...
%! margin (MODEL, [], zeros (2))
%!error<ClassificationDiscriminant.margin: X must have the same number of predictors as the trained model.> ...
%! margin (MODEL, 1, zeros (2))
%!error<ClassificationDiscriminant.margin: Y must have the same number of rows as X.> ...
%! margin (MODEL, ones (4,2), ones (3,1))

## Test crossval method
%!shared x, y, obj
%! load fisheriris
%! x = meas;
%! y = species;
%! obj = fitcdiscr (x, y, "gamma", 0.4);
%!test
%! CVMdl = crossval (obj);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 10)
%! assert (class (CVMdl.Trained{1}), "CompactClassificationDiscriminant")
%! assert (CVMdl.CrossValidatedModel, "ClassificationDiscriminant")
%!test
%! CVMdl = crossval (obj, "KFold", 3);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 3)
%! assert (class (CVMdl.Trained{1}), "CompactClassificationDiscriminant")
%! assert (CVMdl.CrossValidatedModel, "ClassificationDiscriminant")
%!test
%! CVMdl = crossval (obj, "HoldOut", 0.2);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "CompactClassificationDiscriminant")
%! assert (CVMdl.CrossValidatedModel, "ClassificationDiscriminant")
%!test
%! CVMdl = crossval (obj, "LeaveOut", 'on');
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "CompactClassificationDiscriminant")
%! assert (CVMdl.CrossValidatedModel, "ClassificationDiscriminant")
%!test
%! partition = cvpartition (y, 'KFold', 3);
%! CVMdl = crossval (obj, 'cvPartition', partition);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert (CVMdl.KFold == 3)
%! assert (class (CVMdl.Trained{1}), "CompactClassificationDiscriminant")
%! assert (CVMdl.CrossValidatedModel, "ClassificationDiscriminant")

## Test input validation for crossval method
%!error<ClassificationDiscriminant.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (obj, "kfold")
%!error<ClassificationDiscriminant.crossval: specify only one of the optional Name-Value paired arguments.>...
%! crossval (obj, "kfold", 12, "holdout", 0.2)
%!error<ClassificationDiscriminant.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (obj, "kfold", 'a')
%!error<ClassificationDiscriminant.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (obj, "holdout", 2)
%!error<ClassificationDiscriminant.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (obj, "leaveout", 1)
%!error<ClassificationDiscriminant.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (obj, "cvpartition", 1)
