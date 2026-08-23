## Copyright (C) 2024 Pallav Purbia <pallavpurbia@gmail.com>
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

classdef ClassificationSVM
  ## -*- texinfo -*-
  ## @deftp {statistics} ClassificationSVM
  ##
  ## Support Vector Machine classification
  ##
  ## The @code{ClassificationSVM} class implements a Support Vector Machine
  ## classifier object for one-class or two-class problems, which can predict
  ## responses for new data using the @code{predict} method.
  ##
  ## Support Vector Machine classification is a supervised learning method used
  ## for classification tasks.  It works by finding the optimal hyperplane that
  ## separates classes in the feature space with the maximum margin.  For
  ## non-linearly separable data, it uses kernel functions to map data to a
  ## higher-dimensional space where separation is possible.
  ##
  ## Create a @code{ClassificationSVM} object by using the @code{fitcsvm}
  ## function or the class constructor.
  ##
  ## @seealso{fitcsvm}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix containing the unstandardized predictor data.  Each
    ## column of @var{X} represents one predictor (variable), and each row
    ## represents one observation.  This property is read-only.
    ##
    ## @end deftp
    X                   = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} Y
    ##
    ## Class labels
    ##
    ## Specified as a logical or numeric column vector, or as a character array
    ## or a cell array of character vectors with the same number of rows as the
    ## predictor data.  Each row in @var{Y} is the observed class label for
    ## the corresponding row in @var{X}.  This property is read-only.
    ##
    ## @end deftp
    Y                   = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer value specifying the number of observations in the
    ## training dataset used for training the ClassificationSVM model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} RowsUsed
    ##
    ## Rows used for fitting
    ##
    ## A logical column vector with the same length as the observations in the
    ## original predictor data @var{X}, true for each row that was used for
    ## fitting the ClassificationSVM model.  It is empty, @qcode{[]},
    ## when every observation was used, so a non-empty value means that rows
    ## holding missing values were dropped.  This property is read-only.
    ##
    ## @end deftp
    RowsUsed            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer value specifying the number of predictors in the
    ## training dataset used for training the ClassificationSVM model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumPredictors       = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which they appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames      = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} ClassNames
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
    ClassNames          = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} Sigma
    ##
    ## Predictor standard deviations
    ##
    ## A numeric vector of the same length as the columns in @var{X} containing
    ## the standard deviations of predictor variables.  If the predictor
    ## variables have not been standardized, then @qcode{Sigma} is empty.
    ## This property is read-only.
    ##
    ## Only observations with no missing predictor enter the estimate, and
    ## they are weighted so that each class keeps the share of the
    ## observation weight it carried before any row was set aside.
    ##
    ## @end deftp
    Sigma               = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} Mu
    ##
    ## Predictor means
    ##
    ## A numeric vector of the same length as the columns in @var{X} containing
    ## the means of predictor variables.  If the predictor variables have not
    ## been standardized, then @qcode{Mu} is empty.  This property is read-only.
    ##
    ## Only observations with no missing predictor enter the estimate, and
    ## they are weighted so that each class keeps the share of the
    ## observation weight it carried before any row was set aside.
    ##
    ## @end deftp
    Mu                  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} ModelParameters
    ##
    ## SVM training parameters
    ##
    ## A structure containing the parameters used to train the SVM model with
    ## the following fields: @code{SVMtype}, @code{BoxConstraint},
    ## @code{CacheSize}, @code{KernelScale}, @code{KernelOffset},
    ## @code{KernelFunction}, @code{PolynomialOrder}, @code{Nu},
    ## @code{Tolerance}, and @code{Shrinking}.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} Alpha
    ##
    ## Trained classifier coefficients
    ##
    ## The coefficients of the trained SVM classifier specified as an @math{s*1}
    ## numeric vector, where @math{s} is the number of support vectors equal to
    ## @qcode{sum (obj.IsSupportVector)}.  They are the magnitudes of the dual
    ## coefficients and are never negative; the class each belongs to is given
    ## by the corresponding entry of @qcode{SupportVectorLabels}.
    ## @qcode{Alpha} is populated for every kernel function.  This property is
    ## read-only.
    ##
    ## @end deftp
    Alpha               = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} Beta
    ##
    ## Linear predictor coefficients
    ##
    ## The linear predictor coefficients specified as a @math{p*1} numeric
    ## vector, where @math{p} is the number of predictors.  @qcode{Beta} is
    ## the primal representation of the fitted hyperplane and exists only when
    ## the SVM classifier was trained with a @qcode{'linear'} kernel function;
    ## for any other kernel there is no such representation and @qcode{Beta} is
    ## empty.  It equals
    ## @qcode{obj.SupportVectors' * (obj.Alpha .* obj.SupportVectorLabels)}.
    ## This property is read-only.
    ##
    ## @end deftp
    Beta                = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} Bias
    ##
    ## Bias term
    ##
    ## The bias term specified as a scalar.  This property is read-only.
    ##
    ## @end deftp
    Bias                = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} IsSupportVector
    ##
    ## Support vector indicator
    ##
    ## An @math{N*1} logical vector that flags whether a corresponding
    ## observation in the predictor data matrix is a Support Vector.  @math{N}
    ## is the number of observations in the training data.  This property is
    ## read-only.
    ##
    ## @end deftp
    IsSupportVector     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} SupportVectorLabels
    ##
    ## Support vector class labels
    ##
    ## The support vector class labels specified as an @math{s*1} numeric
    ## vector, where @math{s} is the number of support vectors equal to
    ## @qcode{sum (obj.IsSupportVector)}.  A value of +1 in
    ## @code{SupportVectorLabels} indicates that the corresponding support
    ## vector
    ## belongs to the positive class @qcode{(ClassNames@{2@})}.  A value of -1
    ## indicates that the corresponding support vector belongs to the negative
    ## class @qcode{(ClassNames@{1@})}.  This property is read-only.
    ##
    ## @end deftp
    SupportVectorLabels = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} SupportVectors
    ##
    ## Support vectors
    ##
    ## The support vectors of the trained SVM classifier specified an @math{s*p}
    ## numeric matrix, where @math{s} is the number of support vectors equal to
    ## @qcode{sum (obj.IsSupportVector)}, and @math{p} is the number of
    ## predictor
    ## variables in the predictor data.  This property is read-only.
    ##
    ## @end deftp
    SupportVectors      = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} Prior
    ##
    ## Prior probabilities of the classes
    ##
    ## A numeric row vector with one entry per class, in the order of
    ## @code{ClassNames}, summing to one.  It defaults to the class
    ## frequencies of the training data.  This property is read-only.
    ##
    ## Specified as a row vector with one entry per class, in the order of
    ## @qcode{ClassNames}, and rescaled to sum to one.  It may be given as
    ## @qcode{'empirical'}, @qcode{'uniform'}, a numeric vector, or a
    ## structure with @qcode{ClassNames} and @qcode{ClassProbs} fields, which
    ## assigns each probability by class name rather than by position.
    ##
    ## @end deftp
    Prior               = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} Cost
    ##
    ## Cost of misclassification
    ##
    ## A numeric square matrix, where @code{Cost(i,j)} is the cost of
    ## classifying an observation of class @math{i} as class @math{j}.  It
    ## defaults to zero on the diagonal and one elsewhere.  This property is
    ## read-only.
    ##
    ## @end deftp
    Cost                = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector with one entry per training observation,
    ## normalized to sum to one, as MATLAB reports it.  This property is
    ## read-only.
    ##
    ## Each class carries its prior spread evenly over its own observations,
    ## so an observation of a class weighs @qcode{Prior} for that class
    ## divided by the number of observations it holds.
    ##
    ## @end deftp
    W                   = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} KernelParameters
    ##
    ## Parameters of the kernel function
    ##
    ## A structure with fields @qcode{Function} and @qcode{Scale}, and
    ## @qcode{Order} for a polynomial kernel.  @qcode{Function} names the
    ## kernel as MATLAB names it, so a radial basis kernel reports
    ## @qcode{'gaussian'} whichever spelling was given; the kernel the fit was
    ## handed is unchanged in @qcode{ModelParameters}.  This property is
    ## read-only.
    ##
    ## @end deftp
    KernelParameters = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} BoxConstraints
    ##
    ## Box constraints
    ##
    ## A numeric column vector with one entry per observation, holding the
    ## box constraint the fit applied to it.  It is @qcode{BoxConstraint} for
    ## every observation unless @qcode{Prior} or @qcode{Cost} reweighted the
    ## classes, in which case each class is scaled by the weight it carried
    ## into the fit, normalized so the weights average to one.  This property
    ## is read-only.
    ##
    ## @end deftp
    BoxConstraints = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} OutlierFraction
    ##
    ## Expected fraction of outliers in the training data
    ##
    ## A scalar in @code{[0, 1)}, zero unless one was asked for.
    ##
    ## @strong{Deviation from MATLAB.}  The value is reported as it was given,
    ## but it reaches the fit by a different route: MATLAB removes outliers
    ## iteratively and reports @qcode{Solver} as @qcode{'ISDA'}, where a
    ## nonzero fraction here selects LIBSVM's @math{\nu}-SVC, in which
    ## @math{\nu} bounds the fraction of margin errors.  The two agree on what
    ## the number means and not on how the fit reaches it.  This property is
    ## read-only.
    ##
    ## @end deftp
    OutlierFraction = 0;

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} Nu
    ##
    ## Nu parameter for one-class learning
    ##
    ## A positive scalar, and empty unless the model is a one-class learner,
    ## which is what MATLAB reports.  This property is read-only.
    ##
    ## @end deftp
    Nu = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector of column indices into @code{X} naming the predictors
    ## treated as categorical, and empty when none is.  This property is
    ## read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the model expanded them
    ##
    ## A cell array of character vectors.  It matches @code{PredictorNames}
    ## unless a categorical predictor was expanded into indicator variables.
    ## This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};
    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## A cell array with one entry per predictor, holding that predictor's bin
    ## edges where the learner discretized it before fitting.  It is empty here
    ## and stays empty: this learner fits the predictors as they are, and
    ## MATLAB's reports an empty cell for it as well.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    BinEdges        = {};

  endproperties

  ## The LIBSVM structure the engine works in.  It is ours alone, with no
  ## MATLAB counterpart, so it is kept out of the property listing while
  ## staying readable for anyone who needs the raw model.
  properties (Hidden, GetAccess = public, SetAccess = protected)
    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} Model
    ##
    ## Trained SVM model
    ##
    ## A structure containing the trained model in @qcode{'libsvm'} format.
    ## This property is read-only.
    ##
    ## It is the engine's own structure and has no MATLAB counterpart,
    ## so it is kept out of @code{properties} and out of the online
    ## documentation.  Reading it works exactly as it always did.
    ##
    ## @end deftp
    Model               = [];
  endproperties

  ## Properties a user may set after the model is built.  Each one is
  ## validated by its set method below.
  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationSVM} {property} ScoreTransform
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
    ## @end deftp
    ScoreTransform      = 'none';

  endproperties

  ## Readable by the counterpart class, which copies it, and kept out of
  ## the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    STfun = @(x) x;
  endproperties

  ## Set methods for the properties a user may assign.
  methods (Hidden)

    function this = set.ScoreTransform (this, val)
        name = 'ClassificationSVM';
        try
          [this.STfun, this.ScoreTransform] = parseScoreTransform ...
                                               (val, name);
        catch
          error (strcat ("ClassificationSVM.subsasgn: 'ScoreTransform'", ...
                         " must be a 'function_handle' object."));
        end_try_catch
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
      fprintf ("\n  ClassificationSVM\n\n");
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
      fprintf ("%+25s: [%dx1 double]\n", 'Alpha', numel (this.Alpha));
      if (! isempty (this.Beta))
        fprintf ("%+25s: [%dx1 double]\n", 'Beta', numel (this.Beta));
      endif
      fprintf ("%+25s: %f\n", 'Bias', this.Bias);
      fprintf ("%+25s: [1x1 struct]\n", 'KernelParameters');
      if (! isempty (this.Mu))
        if (numel (this.Mu) < 6)
          str = repmat ({'''%0.4f'''}, 1, numel (this.Mu));
          str = strcat ('[', strjoin (str, ' '), ']');
          out = sprintf (str, this.Mu);
          fprintf ("%+25s: %s\n", 'Mu', out);
          out = sprintf (str, this.Sigma);
          fprintf ("%+25s: %s\n", 'Sigma', out);
        else
          fprintf ("%+25s: [1x%d double]\n", 'Mu', numel (this.Mu));
          fprintf ("%+25s: [1x%d double]\n", 'Sigma', numel (this.Sigma));
        endif
      endif
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {statistics} {@var{obj} =} ClassificationSVM (@var{X}, @var{Y})
    ## @deftypefnx {statistics} {@var{obj} =} ClassificationSVM (@dots{}, @var{name}, @var{value})
    ##
    ## Create a @qcode{ClassificationSVM} class object containing a Support
    ## Vector Machine classification model for one-class or two-class problems.
    ##
    ## @code{@var{obj} = ClassificationSVM (@var{X}, @var{Y})} returns a
    ## ClassificationSVM object, with @var{X} as the predictor data and @var{Y}
    ## containing the class labels of observations in @var{X}.
    ##
    ## @itemize
    ## @item
    ## @code{X} must be a @math{N*P} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.  @var{X} will be used to train the SVM model.
    ## @item
    ## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}.  @var{Y} can be either
    ## numeric, logical, or cell array of character vectors.  It must have same
    ## numbers of rows as @var{X}.
    ## @end itemize
    ##
    ## @code{@var{obj} = ClassificationSVM (@dots{}, @var{name}, @var{value})}
    ## returns a ClassificationSVM object with parameters specified by the
    ## following @qcode{@var{name}, @var{value}} paired input arguments:
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
    ## labels, @var{Y}, used for fitting the SVM model. @qcode{ClassNames} are
    ## of the same type as the class labels in @var{Y}.
    ##
    ## @item @qcode{'ScoreTransform'} @tab A user-defined function handle
    ## or a character vector specifying one of the following builtin functions
    ## specifying the transformation applied to predicted classification scores.
    ## Supported values include @qcode{'doublelogit'}, @qcode{'invlogit'},
    ## @qcode{'ismax'}, @qcode{'logit'}, @qcode{'none'}, @qcode{'identity'},
    ## @qcode{'sign'}, @qcode{'symmetric'}, @qcode{'symmetricismax'}, and
    ## @qcode{'symmetriclogit'}.
    ##
    ## @item @qcode{'Standardize'} @tab A logical scalar specifying whether
    ## to standardize the predictor variables.  Default is @qcode{false}.
    ##
    ## @item @qcode{'SVMtype'} @tab A character vector specifying the type
    ## of SVM to use.  Supported values are @qcode{'c_svc'} (C-support vector
    ## classification), @qcode{'nu_svc'} (nu-support vector classification), and
    ## @qcode{'one_class_svm'} (one-class SVM).
    ##
    ## @item @qcode{'KernelFunction'} @tab A character vector specifying
    ## the kernel function to use.  Supported values are @qcode{'linear'},
    ## @qcode{'rbf'} or @qcode{'gaussian'}, @qcode{'polynomial'}, and
    ## @qcode{'sigmoid'}.
    ##
    ## @item @qcode{'PolynomialOrder'} @tab A positive integer specifying
    ## the order of the polynomial kernel function.  Default is 3.
    ##
    ## @item @qcode{'KernelScale'} @tab A positive scalar specifying the
    ## kernel scale parameter.  Default is 1.
    ##
    ## @item @qcode{'KernelOffset'} @tab A non-negative scalar specifying
    ## the kernel offset parameter.  Default is 0.
    ##
    ## @item @qcode{'BoxConstraint'} @tab A positive scalar specifying the
    ## box constraint parameter.  Default is 1.
    ##
    ## @item @qcode{'Nu'} @tab A positive scalar in the range (0,1]
    ## specifying the nu parameter for nu-SVM and one-class SVM. Default is 0.5.
    ##
    ## @item @qcode{'CacheSize'} @tab A positive scalar specifying the
    ## cache size in MB.  Default is 1000.
    ##
    ## @item @qcode{'Tolerance'} @tab A positive scalar specifying the
    ## tolerance of termination criterion.  Default is 1e-6.
    ##
    ## @item @qcode{'Shrinking'} @tab Either 0 or 1 specifying whether to
    ## use the shrinking heuristics.  Default is 1.
    ##
    ## @item @qcode{'OutlierFraction'} @tab A positive scalar in the range
    ## [0,1) specifying the fraction of outliers for one-class SVM.
    ## @end multitable
    ##
    ## @seealso{fitcsvm}
    ## @end deftypefn
    function this = ClassificationSVM (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationSVM: too few input arguments.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error ("ClassificationSVM: number of rows in X and Y must be equal.");
      endif

      ## Assign original X and Y data to the ClassificationSVM object
      this.X = X;
      this.Y = Y;

      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);

      ## Set default values before parsing optional parameters
      SVMtype                 = 'c_svc';
      KernelFunction          = [];
      KernelScale             = 1;
      KernelOffset            = 0;
      PolynomialOrder         = 3;
      BoxConstraint           = 1;
      Nu                      = 0.5;
      OutlierFraction         = 0;
      CacheSize               = 1000;
      Tolerance               = 1e-6;
      Shrinking               = 1;
      Standardize             = false;
      ResponseName            = [];
      PredictorNames          = [];
      ClassNames              = [];
      Prior                   = [];
      Cost                    = [];

      ## Parse extra parameters
      SVMtype_override = true;
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'standardize'
            Standardize = varargin{2};
            if (! (Standardize == true || Standardize == false))
              error (strcat ("ClassificationSVM: 'Standardize' must", ...
                             " be either true or false."));
            endif

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("ClassificationSVM: 'PredictorNames' must", ...
                             " be supplied as a cellstring array."));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat ("ClassificationSVM: 'PredictorNames' must", ...
                             " have the same number of columns as X."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("ClassificationSVM: 'ResponseName' must", ...
                             " be a character vector."));
            endif

          case 'classnames'
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames) ||
                   islogical (ClassNames) || ischar (ClassNames)))
              error (strcat ("ClassificationSVM: 'ClassNames' must be a", ...
                             " cell array of character vectors, a logical", ...
                             " vector, a numeric vector, or a character array."));
            endif
            ## Check that all class names are available in gnY
            if (iscellstr (ClassNames) || ischar (ClassNames))
              ClassNames = cellstr (ClassNames);
              if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                   ClassNames, 'UniformOutput', false))))
                error (strcat ("ClassificationSVM: not all 'ClassNames'", ...
                               " are present in Y."));
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, 'UniformOutput', false))))
                error (strcat ("ClassificationSVM: not all 'ClassNames'", ...
                               " are present in Y."));
              endif
            endif

          case 'prior'
            Prior = varargin{2};
            if (! (isstruct (Prior)
                   || (isnumeric (Prior) && isvector (Prior) && all (Prior >= 0)
                       && any (Prior > 0))
                   || (ischar (Prior)
                       && any (strcmpi (Prior, {'empirical', 'uniform'})))))
              error (strcat ("ClassificationSVM: 'Prior' must be a", ...
                             " non-negative numeric vector, 'empirical'", ...
                             " or 'uniform'."));
            endif

          case 'cost'
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost) && all (Cost(:) >= 0)))
              error (strcat ("ClassificationSVM: 'Cost' must be a", ...
                             " non-negative square matrix."));
            endif

          case 'scoretransform'
            name = 'ClassificationSVM';
            [this.STfun, this.ScoreTransform] = parseScoreTransform ...
                                                 (varargin{2}, name);

          case 'svmtype'
            SVMtype = varargin{2};
            SVMtype_override = false;
            if (! any (strcmp (SVMtype, {'c_svc', 'nu_svc', 'one_class_svm'})))
              error (strcat ("ClassificationSVM: 'SVMtype' must be", ...
                             " 'c_svc', 'nu_svc', or 'one_class_svm'."));
            endif

          case 'outlierfraction'
            Nu = varargin{2};
            OutlierFraction = Nu;
            if (! (isscalar (Nu) && Nu >= 0 && Nu < 1))
              error (strcat ("ClassificationSVM: 'OutlierFraction' must", ...
                             " be a positive scalar in the range 0 =<", ...
                             " OutlierFraction < 1."));
            endif
            if (Nu > 0)
              SVMtype = 'nu_svc';
            endif

          case 'kernelfunction'
            KernelFunction = varargin{2};
            if (! ischar (KernelFunction))
              error (strcat ("ClassificationSVM: 'KernelFunction' must", ...
                             " be a character vector."));
            endif
            KernelFunction = tolower (KernelFunction);
            if (! any (strcmpi (KernelFunction, ...
                       {'linear', 'rbf', 'gaussian', 'polynomial', 'sigmoid'})))
              error ("ClassificationSVM: unsupported Kernel function.");
            endif

          case 'polynomialorder'
            PolynomialOrder = varargin{2};
            if (! (isnumeric (PolynomialOrder) && isscalar (PolynomialOrder)
                   && PolynomialOrder > 0 && mod (PolynomialOrder, 1) == 0))
              error (strcat ("ClassificationSVM: 'PolynomialOrder' must", ...
                             " be a positive integer."));
            endif

          case 'kernelscale'
            KernelScale = varargin{2};
            if (! (isscalar (KernelScale) && KernelScale > 0))
              error (strcat ("ClassificationSVM: 'KernelScale'", ...
                             " must be a positive scalar."));
            endif

          case 'kerneloffset'
            KernelOffset = varargin{2};
            if (! (isnumeric (KernelOffset) && isscalar (KernelOffset)
                                            && KernelOffset >= 0))
              error (strcat ("ClassificationSVM: 'KernelOffset' must", ...
                             " be a non-negative scalar."));
            endif

          case 'boxconstraint'
            BoxConstraint = varargin{2};
            if (! (isscalar (BoxConstraint) && BoxConstraint > 0))
              error (strcat ("ClassificationSVM: 'BoxConstraint' must", ...
                             " be a positive scalar."));
            endif

          case 'nu'
            Nu = varargin{2};
            if (SVMtype_override)
              SVMtype = 'one_class_svm';
            endif
            if (! (isscalar (Nu) && Nu > 0 && Nu <= 1))
              error (strcat ("ClassificationSVM: 'Nu' must be a positive", ...
                             " scalar in the range 0 < Nu <= 1."));
            endif

          case 'cachesize'
            CacheSize = varargin{2};
            if (! (isscalar (CacheSize) && CacheSize > 0))
              error (strcat ("ClassificationSVM: 'CacheSize' must", ...
                             " be a positive scalar."));
            endif

          case 'tolerance'
            Tolerance = varargin{2};
            if (! (isscalar (Tolerance) && Tolerance >= 0))
              error (strcat ("ClassificationSVM: 'Tolerance' must", ...
                             " be a positive scalar."));
            endif

          case 'shrinking'
            Shrinking = varargin{2};
            if (! (ismember (Shrinking, [0, 1]) && isscalar (Shrinking)))
              error ("ClassificationSVM: 'Shrinking' must be either 0 or 1.");
            endif

          otherwise
            error (strcat ("ClassificationSVM: invalid parameter name", ...
                           " in optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Get number of variables in training data
      ndims_X = columns (X);

      ## Assign the number of predictors to the ClassificationSVM object
      this.NumPredictors = ndims_X;

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

      ## An observation is dropped only when its response is missing.  A row
      ## whose predictors hold missing values is kept and reported as used,
      ## while the fit below draws on the complete observations alone.
      RowsUsed  = ! isnan (gY);
      Yret      = Y(RowsUsed);
      Xret      = X(RowsUsed, :);
      this.X    = Xret;
      this.Y    = Yret;
      cobs      = ! any (isnan (Xret), 2);
      Y         = Yret(cobs);
      X         = Xret(cobs, :);

      ## Renew groups in Y over the retained observations, so a class held
      ## only by a row with missing predictors is still a class of the model
      [gret, gnY, glY] = grp2idx (Yret);
      gY = gret(cobs);
      nclasses = numel (gnY);
      this.ClassNames = glY;  # Keep the same type as Y

      ## Resolve Prior and Cost against the classes that survived.  Prior
      ## defaults to the frequencies of the training data and Cost to zero on
      ## the diagonal and one elsewhere, which is what MATLAB reports.
      if (isstruct (Prior))
        Prior = priorFromStruct (Prior, this.ClassNames, ...
                                 'ClassificationSVM');
      endif
      freq = accumarray (gY(:), 1, [nclasses, 1])' / numel (gY);
      if (isempty (Prior) || (ischar (Prior) && strcmpi (Prior, 'empirical')))
        Prior = freq;
      elseif (ischar (Prior) && strcmpi (Prior, 'uniform'))
        Prior = ones (1, nclasses) / nclasses;
      else
        if (numel (Prior) != nclasses)
          error (strcat ("ClassificationSVM: 'Prior' must have one entry", ...
                         " per class."));
        endif
        Prior = Prior(:)' / sum (Prior);
      endif
      if (isempty (Cost))
        Cost = ones (nclasses) - eye (nclasses);
      elseif (rows (Cost) != nclasses)
        error (strcat ("ClassificationSVM: the number of rows and columns", ...
                       " in 'Cost' must correspond to the classes in Y."));
      endif
      this.Prior = Prior;
      this.Cost = Cost;

      ## If only one class available, force 'SVMtype' to 'one_class_svm'
      if (nclasses == 1)
        if (! SVMtype_override && ! strcmp (SVMtype, 'one_class_svm'))
          error (strcat ("ClassificationSVM: cannot train a binary", ...
                         " problem with only one class available."));
        endif
        SVMtype = 'one_class_svm';
        if (isempty (KernelFunction))
          KernelFunction = 'rbf';
        endif
      else
        if (isempty (KernelFunction))
          KernelFunction = 'linear';
        endif
      endif

      ## Check that we are dealing only with one-class or binary classification
      if (nclasses > 2)
        error (strcat ("ClassificationSVM: can only be used for", ...
                       " one-class or two-class learning."));
      endif

      ## Force Y into numeric
      if (! isnumeric (Y))
        Y = gY;
      endif

      ## Force Y labels to -1 and +1 to avoid numeric issues with different
      ## compiling options; see https://github.com/cjlin1/libsvm/issues/220
      if (nclasses == 2)
        Y(Y == 2) = -1;
      endif

      ## Check X contains valid data
      if (! (isnumeric (X) && isfinite (X)))
        error ("ClassificationSVM: invalid values in X.");
      endif

      ## Assign the number of observations and their corresponding indices
      ## on the original data, which will be used for training the model,
      ## to the ClassificationSVM object
      this.NumObservations = rows (this.X);
      ## RowsUsed is left empty when every observation was used, as in MATLAB
      if (all (RowsUsed))
        this.RowsUsed = [];
      else
        this.RowsUsed = RowsUsed;
      endif

      ## Handle the Standardize option.  The model must be fitted on the
      ## scale it predicts on: predict and resubPredict standardize their
      ## input from Mu and Sigma, so the training data is standardized here
      ## as well.
      ## The support vectors are therefore stored standardized too, which is
      ## the scale svmpredict receives.
      if (Standardize)
        ## Mu and Sigma weight the complete observations so that each class
        ## keeps the share of the observation weight it carried before any row
        ## was set aside, which is what MATLAB reports.
        sw = zeros (rows (X), 1);
        for k = 1:numel (gnY)
          ck = (gY == k);
          if (any (ck))
            sw(ck) = (sum (gret == k) / numel (gret)) / sum (ck);
          endif
        endfor
        sw = sw / sum (sw);
        this.Mu = sum (sw .* X, 1);
        Zs = X - this.Mu;
        this.Sigma = sqrt (sum (sw .* Zs .^ 2, 1) / (1 - sum (sw .^ 2)));
        this.Sigma(this.Sigma == 0) = 1;  # predictor is constant
        X = (X - this.Mu) ./ this.Sigma;
      else
        this.Sigma = [];
        this.Mu = [];
      endif

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

      ## No predictor is treated as categorical, so the expanded names are the
      ## predictor names themselves, and every observation carries the same
      ## weight, normalized to sum to one as MATLAB reports it.
      this.CategoricalPredictors = [];
      this.ExpandedPredictorNames = PredictorNames;
      this.W = priorWeights (this.Prior, gY, this.NumObservations);

      ## Set svmtrain parameters for SVMtype and KernelFunction
      switch (SVMtype)
        case 'c_svc'
          s = 0;
        case 'nu_svc'
          s = 1;
        case 'one_class_svm'
          s = 2;
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

      ## svmpredict:
      ##    '-s':  SVMtype
      ##    '-t':  KernelFunction
      ##    '-g':  Gamma
      ##    '-d':  PolynomialOrder
      ##    '-r':  KernelOffset
      ##    '-c':  BoxConstraint
      ##    '-n':  Nu
      ##    '-m':  CacheSize
      ##    '-e':  Tolerance
      ##    '-h':  Shrinking

      ## Build options string for svmtrain function
      str_options = strcat ("-s %d -t %d -g %f -d %d -r %f", ...
                            " -c %f -n %f -m %f -e %e -h %d -q");
      svm_options = sprintf (str_options, s, t, g, PolynomialOrder, ...
                             KernelOffset, BoxConstraint, Nu, ...
                             CacheSize, Tolerance, Shrinking);

      ## Prior and Cost enter the fit through LIBSVM's per-class weights,
      ## which scale the box constraint of each class: the weight of class i
      ## is the cost it carries, Prior(i) times the cost of getting it wrong,
      ## divided by how often it actually occurs.  The default prior is the
      ## observed frequency and the default cost is one, so the default weight
      ## is one for every class and the fit is the same as with no weights.
      cw = ones (1, nclasses);
      if (nclasses == 2 && ! strcmp (SVMtype, 'one_class_svm'))
        cw = (Prior .* sum (Cost, 2)') ./ max (freq, eps);
        cw = cw / mean (cw);
        if (any (abs (cw - 1) > 1e-12))
          ## LIBSVM names the classes by the labels it was given, which are
          ## +1 and -1 here, in the order of ClassNames.
          svm_options = sprintf ("%s -w1 %.16g -w-1 %.16g", svm_options, ...
                                 cw(1), cw(2));
        endif
      endif

      ## Train the SVM model using svmtrain from libsvm
      Model = svmtrain (Y, X, svm_options);
      this.Model = Model;

      ## Populate ClassificationSVM object properties.  LIBSVM returns the
      ## dual coefficients already multiplied by the class sign, whereas
      ## MATLAB keeps their magnitudes in ALPHA and the sign in
      ## SUPPORTVECTORLABELS, so the two are separated here.
      this.Alpha = abs (Model.sv_coef);
      this.Bias = Model.rho;

      ## One label per support vector, in the order of SupportVectors, taking
      ## the sign from the coefficients themselves.  LIBSVM's sign is opposite
      ## to the labelling MATLAB reports.
      this.SupportVectorLabels = -sign (Model.sv_coef);

      ## BETA holds the primal coefficients, one per predictor, and exists
      ## only for a linear kernel; for any other kernel there is no primal
      ## representation and MATLAB leaves it empty.
      if (t == 0)
        this.Beta = Model.SVs' * (this.Alpha .* this.SupportVectorLabels);
      else
        this.Beta = [];
      endif

      this.IsSupportVector = false (this.NumObservations, 1);
      this.IsSupportVector(Model.sv_indices) = true;
      this.SupportVectors = Model.SVs;

      ## The kernel, the per-observation box constraints and the two
      ## one-class parameters, in the shapes MATLAB reports them.  The box
      ## constraints are the scalar scaled by the very weights the fit was
      ## given above, so they describe what was solved rather than what was
      ## asked for.
      this.KernelParameters = svmKernelParams (KernelFunction, KernelScale, ...
                                               PolynomialOrder);
      this.BoxConstraints = BoxConstraint * cw(gY)(:);
      this.OutlierFraction = OutlierFraction;
      if (strcmp (SVMtype, 'one_class_svm'))
        this.Nu = Nu;
      endif

      ## Populate ModelParameters structure
      params = struct ('SVMtype', SVMtype, 'BoxConstraint', BoxConstraint, ...
                       'CacheSize', CacheSize, 'KernelScale', KernelScale, ...
                       'KernelOffset', KernelOffset, 'KernelFunction', ...
                        KernelFunction, 'PolynomialOrder', PolynomialOrder, ...
                       'Nu', Nu, 'Tolerance',Tolerance, ...
                       'Shrinking', Shrinking);
      this.ModelParameters = params;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{obj} =} discardSupportVectors (@var{obj})
    ##
    ## Discard the support vectors of a linear SVM model.
    ##
    ## @code{@var{obj} = discardSupportVectors (@var{obj})} empties
    ## @code{Alpha}, @code{SupportVectors} and
    ## @code{SupportVectorLabels}, leaving @code{Beta} and @code{Bias} to
    ## decide every prediction.  A linear kernel needs nothing else, so the
    ## returned model predicts what it predicted before while carrying one
    ## vector in place of many.
    ##
    ## The kernel must be linear.  Under any other the support vectors are
    ## part of the decision function and cannot be dropped.  Discarding twice
    ## is not an error and changes nothing.
    ##
    ## @seealso{fitcsvm, ClassificationSVM, CompactClassificationSVM}
    ## @end deftypefn
    function this = discardSupportVectors (this)

      if (nargin != 1)
        print_usage ();
      endif

      if (! strcmpi (this.ModelParameters.KernelFunction, 'linear'))
        error (strcat ("ClassificationSVM.discardSupportVectors: you", ...
                       " cannot discard support vectors for a non-linear", ...
                       " kernel."));
      endif

      ## The engine keeps its own copy of the support vectors, so emptying
      ## the properties alone would free nothing.  Collapsing the model onto
      ## the single vector that decides it leaves every scoring path as it
      ## was, svmpredict going on being the engine over one vector.
      this.Model = discardSVs (this.Model);
      this.Alpha = [];
      this.SupportVectors = [];
      this.SupportVectorLabels = [];
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationSVM} {[@var{label}, @var{score}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the Support Vector Machine
    ## classification model from a ClassificationSVM object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## predictor data in @code{obj.X} and corresponding labels, @code{obj.Y},
    ## stored in the ClassificationSVM model, @var{obj}.  For one-class SVM
    ## model, +1 or -1 is returned.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationSVM} class object.
    ## @item
    ## @var{XC} must be an @math{M*P} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the SVM model in
    ## @var{obj}.
    ## @end itemize
    ##
    ## @code{[@var{label}, @var{score}] = predict (@var{obj}, @var{XC})} also
    ## returns @var{score}, which contains the decision values for each each
    ## prediction.  Alternatively, @var{score} can contain the posterior
    ## probabilities if the ScoreTransform has been previously set using the
    ## @code{fitPosterior} method.
    ##
    ## @seealso{ClassificationSVM, fitcsvm, ClassificationSVM.fitPosterior}
    ## @end deftypefn
    function [labels, scores] = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("ClassificationSVM.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("ClassificationSVM.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("ClassificationSVM.predict:", ...
                       " XC must have the same number of", ...
                       " predictors as the trained model."));
      endif

      ## Standardize (if necessary)
      if (! isempty (this.Mu))
        XC = (XC - this.Mu) ./ this.Sigma;
      endif

      ## Predict labels and scores from new data
      [out, ~, scores] = svmpredict (ones (rows (XC), 1), XC, this.Model, '-q');

      ## Expand scores for two classes
      if (numel (this.ClassNames) == 2)
        scores = [scores, -scores];
      endif

      ## Translate labels to classnames
      if (iscellstr (this.Y))
        labels = cell (rows (XC), 1);
        labels(out==1) = this.ClassNames{1};
        labels(out!=1) = this.ClassNames{2};
      elseif (islogical (this.Y))
        labels = false (rows (XC), 1);
      elseif (isnumeric (this.Y))
        labels = zeros (rows (XC), 1);
      elseif (ischar (this.Y))
        labels = char (zeros (rows (XC), size (this.Y, 2)));
      endif
      if (! iscellstr (this.Y))
        labels(out==1) = this.ClassNames(1);
        labels(out!=1) = this.ClassNames(2);
      endif

      ## Apply ScoreTransform
      scores = this.STfun (scores);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{label} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationSVM} {[@var{label}, @var{score}] =} resubPredict (@var{obj})
    ##
    ## Classify the training data using the trained Support Vector Machine
    ## classification object.
    ##
    ## @code{@var{label} = resubPredict (@var{obj})} returns the vector of
    ## labels predicted for the corresponding instances in the training data,
    ## using the predictor data in @code{obj.X} and corresponding labels,
    ## @code{obj.Y}, stored in the Support Vector Machine classification model,
    ## @var{obj}. For one-class model, +1 or -1 is returned.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationSVM} class object.
    ## @end itemize
    ##
    ## @code{[@var{label}, @var{scores}] = resubPredict (@var{obj}} also
    ## returns @var{scores}, which contains the decision values for each each
    ## prediction.   Alternatively, @var{scores} can contain the posterior
    ## probabilities if the ScoreTransform has been previously set using the
    ## @code{fitPosterior} method.
    ##
    ## @seealso{fitcsvm, ClassificationSVM.fitPosterior}
    ## @end deftypefn
    function [labels, scores] = resubPredict (this)

      ## X and Y hold exactly the observations the model retained
      X = this.X;
      Y = this.Y;

      ## Standardize (if necessary)
      if (! isempty (this.Mu))
        X = (X - this.Mu) ./ this.Sigma;
      endif

      ## Predict labels and scores from new data
      [out, ~, scores] = svmpredict (ones (rows (X), 1), X, this.Model, '-q');

      ## Expand scores for two classes
      if (numel (this.ClassNames) == 2)
        scores = [scores, -scores];
      endif

      ## Translate labels to classnames
      if (iscellstr (this.Y))
        labels = cell (rows (X), 1);
        labels(out==1) = this.ClassNames{1};
        labels(out!=1) = this.ClassNames{2};
      elseif (islogical (this.Y))
        labels = false (rows (X), 1);
      elseif (isnumeric (this.Y))
        labels = zeros (rows (X), 1);
      elseif (ischar (this.Y))
        labels = char (zeros (rows (X), size (this.Y, 2)));
      endif
      if (! iscellstr (this.Y))
        labels(out==1) = this.ClassNames(1);
        labels(out!=1) = this.ClassNames(2);
      endif

      ## Apply ScoreTransform
      scores = this.STfun (scores);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margins for Support Vector Machine classifier.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns
    ## the classification margins for @var{obj} with data @var{X} and
    ## classification @var{Y}. @var{m} is a numeric vector of length size (X,1).
    ##
    ## @itemize
    ## @item
    ## @var{obj} is a @var{ClassificationSVM} object trained on @code{X}
    ## and @code{Y}.
    ## @item
    ## @var{X} must be a @math{N*P} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @var{Y} is @math{N*1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of Rows as @var{X}.
    ## @end itemize
    ##
    ## The classification margin for each observation is the difference between
    ## the classification score for the true class and the maximal
    ## classification score for the false classes.
    ##
    ## @seealso{fitcsvm, ClassificationSVM}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationSVM.margin: too few input arguments.");
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationSVM.margin: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat ("ClassificationSVM.margin: X must have the same", ...
                       " number of predictors as the trained model."));
      endif

      ## Check for valid Y
      if (isempty (Y))
        error ("ClassificationSVM.margin: Y is empty.");
      elseif (rows (X) != rows (Y))
        error (strcat ("ClassificationSVM.margin: Y must have", ...
                       " the same number of rows as X."));
      endif

      ## Y may be the class labels, which is what this method documents and
      ## what MATLAB accepts, or already the +1/-1 coding the solver works in.
      ## It used to be the latter only, so passing the labels the docstring
      ## promises reached LIBSVM as a cell array and raised its own message.
      Ypm = svmPlusMinus (Y, this.ClassNames);
      [~, ~, dec_values_L] = svmpredict (Ypm, X, this.Model, '-q');
      m = 2 * Ypm .* dec_values_L;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationSVM} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute loss for a trained ClassificationSVM object.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} computes the loss,
    ## @var{L}, using the default loss function @qcode{'classiferror'}.
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{ClassificationSVM} object trained on
    ## @code{X} and @code{Y}.
    ## @item
    ## @code{X} must be a @math{N*P} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of Rows as @var{X}.
    ## @end itemize
    ##
    ## @code{@var{L} = loss (@dots{}, @var{name}, @var{value})} allows
    ## additional options specified by @var{name}-@var{value} pairs:
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab Specifies the loss function to use.
    ## Can be a function handle with four input arguments (C, S, W, Cost)
    ## which returns a scalar value or one of:
    ## 'binodeviance', 'classifcost', 'classiferror', 'exponential',
    ## 'hinge', 'logit','mincost', 'quadratic'.
    ## @itemize
    ## @item
    ## @code{C} is a logical matrix of size @math{N*K}, where @math{N} is the
    ## number of observations and @math{K} is the number of classes.
    ## The element @code{C(i,j)} is true if the class label of the i-th
    ## observation is equal to the j-th class.
    ## @item
    ## @code{S} is a numeric matrix of size @math{N*K}, where each element
    ## represents the classification score for the corresponding class.
    ## @item
    ## @code{W} is a numeric vector of length @math{N}, representing
    ## the observation weights.
    ## @item
    ## @code{Cost} is a @math{K*K} matrix representing the misclassification
    ## costs.
    ## @end itemize
    ##
    ## @item @qcode{'Weights'} @tab Specifies observation weights, must be
    ## a numeric vector of length equal to the number of rows in X.
    ## Default is @code{ones (size (X, 1))}. loss normalizes the weights so that
    ## observation weights in each class sum to the prior probability of that
    ## class. When you supply Weights, loss computes the weighted
    ## classification loss.
    ##
    ## @end multitable
    ##
    ## @seealso{ClassificationSVM}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationSVM.loss: too few input arguments.");
      endif

      if (mod (nargin, 2) == 0)
        error ("ClassificationSVM.loss: Name-Value arguments must be in pairs.");
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationSVM.loss: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat ("ClassificationSVM.loss: X must have the same", ...
                       " number of predictors as the trained model."));
      endif

      ## Check for valid Y
      if (isempty (Y))
        error ("ClassificationSVM.loss: Y is empty.");
      elseif (rows (X)!= rows (Y))
        error (strcat ("ClassificationSVM.loss: Y must have the same", ...
                       " number of rows as X."));
      endif

      ## Set default values before parsing optional parameters
      LossFun = 'classiferror';
      Weights = ones (size (X, 1), 1);

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'lossfun'
            LossFun = varargin{2};
            if (! (ischar (LossFun)))
              error (strcat ("ClassificationSVM.loss: 'LossFun'", ...
                             " must be a character vector."));
            endif
            LossFun = tolower (LossFun);
            if (! any (strcmpi (LossFun, {'binodeviance', 'classiferror', ...
                                          'classifcost', 'exponential', ...
                                          'hinge', 'logit', 'mincost', ...
                                          'quadratic'})))
              error ("ClassificationSVM.loss: unsupported Loss function.");
            endif

          case 'weights'
            Weights = varargin{2};
            ## Validate if weights is a numeric vector
            if (! (isnumeric (Weights) && isvector (Weights)))
              error (strcat ("ClassificationSVM.loss: 'Weights'", ...
                             " must be a numeric vector."));
            endif

            ## Check if the size of weights matches the number of rows in X
            if (numel (Weights) != size (X, 1))
              error (strcat ("ClassificationSVM.loss: size of 'Weights'", ...
                             " must be equal to the number of rows in X."));
            endif

          otherwise
            error (strcat ("ClassificationSVM.loss: invalid parameter", ...
                           " name in optional pair arguments."));
          endswitch
        varargin(1:2) = [];
      endwhile

      ## Compute the classification score
      [~, ~, dec_values_L] = svmpredict (Y, X, this.Model, '-q');

        ## Compute the margin
        margin = Y .* dec_values_L;

        ## Compute the loss based on the specified loss function
        switch (LossFun)
          case 'classiferror'
            L = mean ((margin <= 0) .* Weights);

          case 'hinge'
            L = mean (max (0, 1 - margin) .* Weights);

          case 'logit'
            L = mean (log (1 + exp (-margin)) .* Weights);

          case 'exponential'
            L = mean (exp (-margin) .* Weights);

          case 'quadratic'
            L = mean (((1 - margin) .^2) .* Weights);

          case 'binodeviance'
            L = mean (log (1 + exp (-2 * margin)) .* Weights);

          case 'mincost'
            ## Each observation is assigned to the class of least expected
            ## cost, and charged what that assignment actually costs given
            ## its true class.  Y is the +1/-1 coding the margin above uses,
            ## in which +1 is the first of ClassNames and -1 the second.
            [~, scores] = predict (this, X);
            true_idx = ones (rows (X), 1);
            true_idx(Y == -1) = 2;
            L = 0;
            for i = 1:rows (X)
              [~, k] = min (scores(i,:) * this.Cost);
              L = L + Weights(i) * this.Cost(true_idx(i), k);
            endfor
            L = L / rows (X);

          case 'classifcost'
            ## What the model's own prediction costs, given the true class
            pred_idx = ones (rows (X), 1);
            pred_idx(dec_values_L <= 0) = 2;
            true_idx = ones (rows (X), 1);
            true_idx(Y == -1) = 2;
            L = 0;
            for i = 1:rows (X)
              L = L + Weights(i) * this.Cost(true_idx(i), pred_idx(i));
            endfor
            L = L / rows (X);

          otherwise
            error ("ClassificationSVM.loss: unsupported Loss function.");
        endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {ClassificationSVM} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute resubstitution loss for a trained ClassificationSVM object.
    ##
    ## @code{@var{L} = resubLoss (@var{obj})} computes the resubstitution loss,
    ## @var{L}, using the default loss function @qcode{'classiferror'}.
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{ClassificationSVM} object trained on
    ## @code{X} and @code{Y}.
    ## @end itemize
    ##
    ## @code{@var{L} = resubLoss (@dots{}, @var{name}, @var{value})} allows
    ## additional options specified by @var{name}-@var{value} pairs:
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab Specifies the loss function to use.
    ## Can be a function handle with four input arguments (C, S, W, Cost)
    ## which returns a scalar value or one of:
    ## 'binodeviance', 'classifcost', 'classiferror', 'exponential',
    ## 'hinge', 'logit','mincost', 'quadratic'.
    ## @itemize
    ## @item
    ## @code{C} is a logical matrix of size @math{N*K}, where @math{N} is the
    ## number of observations and @math{K} is the number of classes.
    ## The element @code{C(i,j)} is true if the class label of the i-th
    ## observation is equal to the j-th class.
    ## @item
    ## @code{S} is a numeric matrix of size @math{N*K}, where each element
    ## represents the classification score for the corresponding class.
    ## @item
    ## @code{W} is a numeric vector of length @math{N}, representing
    ## the observation weights.
    ## @item
    ## @code{Cost} is a @math{K*K} matrix representing the misclassification
    ## costs.
    ## @end itemize
    ##
    ## @item @qcode{'Weights'} @tab Specifies observation weights, must be
    ## a numeric vector of length equal to the number of rows in X.
    ## Default is @code{ones (size (X, 1))}. loss normalizes the weights so that
    ## observation weights in each class sum to the prior probability of that
    ## class. When you supply Weights, loss computes the weighted
    ## classification loss.
    ##
    ## @end multitable
    ##
    ## @seealso{ClassificationSVM}
    ## @end deftypefn
    function L = resubLoss (this, varargin)

      if (mod (nargin, 2) != 1)
        error (strcat ("ClassificationSVM.resubLoss: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      ## Set default values before parsing optional parameters
      LossFun = 'classiferror';
      Weights = ones (size (this.X, 1), 1);

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))

          case 'lossfun'
            LossFun = varargin{2};
            if (! ischar (LossFun))
              error (strcat ("ClassificationSVM.resubLoss: 'LossFun'", ...
                             " must be a character vector."));
            endif
            LossFun = tolower (LossFun);
            if (! any (strcmpi (LossFun, {'binodeviance', 'classiferror', ...
                                          'classifcost', 'exponential', ...
                                          'hinge', 'logit', 'mincost', ...
                                          'quadratic'})))
              error (strcat ("ClassificationSVM.resubLoss: unsupported", ...
                             " Loss function."));
            endif

          case 'weights'
            Weights = varargin{2};
            ## Validate if weights is a numeric vector
            if (! (isnumeric (Weights) && isvector (Weights)))
              error (strcat ("ClassificationSVM.resubLoss: 'Weights'", ...
                             " must be a numeric vector."));
            endif

            ## Check if the size of weights matches the number of rows in X
            if (numel (Weights) != size (this.X, 1))
              error (strcat ("ClassificationSVM.resubLoss: size", ...
                             " of 'Weights' must be equal to the", ...
                             " number of rows in X."));
            endif

          otherwise
            error (strcat ("ClassificationSVM.resubLoss: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      ## The loss of the model on its own training data.  This used to
      ## recompute every loss here from this.Y, which holds the labels as
      ## they were given: the margin wants the +1/-1 coding, so a 1/2 coded
      ## response scored 0.49 where the answer was 0.01, and a logical one
      ## did not reach LIBSVM at all.  The rows dropped for missing values
      ## were included as well.
      used = true (rows (this.X), 1);
      Xu = this.X(used, :);
      gY = grp2idx (this.Y(used, :));
      Ypm = ones (rows (Xu), 1);
      Ypm(gY == 2) = -1;
      L = loss (this, Xu, Ypm, 'LossFun', LossFun, 'Weights', Weights);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationSVM} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross Validate a Support Vector Machine classification object.
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
    ## @seealso{fitcsvm, ClassificationSVM, cvpartition,
    ## ClassificationPartitionedModel}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)

      ## Check for sufficient input arguments
      if (nargin < 1)
        error ("ClassificationSVM.crossval: too few input arguments.");
      endif

      if (numel (varargin) == 1)
        error (strcat ("ClassificationSVM.crossval: Name-Value arguments", ...
                       " must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("ClassificationSVM.crossval: specify only one of", ...
                       " the optional Name-Value paired arguments."));
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
              error (strcat ("ClassificationSVM.crossval: 'KFold' must", ...
                             " be an integer value greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("ClassificationSVM.crossval: 'Holdout' must", ...
                             " be a numeric value between 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat ("ClassificationSVM.crossval: 'Leaveout'", ...
                             " must be either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("ClassificationSVM.crossval: 'CVPartition'",...
                             " must be a 'cvpartition' object."));
            endif

          otherwise
            error (strcat ("ClassificationSVM.crossval: invalid",...
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
    ## @deftypefn  {ClassificationSVM} {@var{CVMdl} =} compact (@var{obj})
    ##
    ## Create a CompactClassificationSVM object.
    ##
    ## @code{@var{CVMdl} = compact (@var{obj})} creates a compact version of the
    ## ClassificationSVM object, @var{obj}.
    ##
    ## @seealso{fitcsvm, ClassificationSVM, CompactClassificationSVM}
    ## @end deftypefn
    function CVMdl = compact (this)
      ## Create a compact model
      CVMdl = CompactClassificationSVM (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationSVM} {@var{e} =} edge (@dots{}, @qcode{"Weights"}, @var{w})
    ##
    ## Classification edge, the mean of the classification margins.
    ##
    ## @code{@var{e} = edge (@var{obj}, @var{X}, @var{Y})} reduces the vector
    ## that @code{margin} returns to a single number, the mean margin over the
    ## rows of @var{X}.  It says how far the model puts the true class ahead of
    ## its nearest rival on average, so a larger edge is a better model, and
    ## unlike a loss it is not bounded above and rewards confidence rather than
    ## bare correctness.
    ##
    ## @code{@var{e} = edge (@dots{}, @qcode{"Weights"}, @var{w})} takes the
    ## weighted mean instead, with one weight per row of @var{X}.
    ##
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      if (nargin < 3)
        error ("ClassificationSVM.edge: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationSVM.edge: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      ## The weights are parsed before anything is computed, so a bad
      ## Name-Value pair is reported as such rather than after a margin.
      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       "ClassificationSVM", "edge");
      m = margin (this, X, Y);
      e = sum (W .* m(:)) / sum (W);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationSVM} {@var{m} =} resubMargin (@var{obj})
    ##
    ## Classification margins of the model on its own training data.
    ##
    ## @code{@var{m} = resubMargin (@var{obj})} is @code{margin} applied to the
    ## observations the model was fitted on, one number per observation.  Being
    ## a resubstitution quantity it is optimistic by construction.
    ##
    ## @end deftypefn
    function m = resubMargin (this)
      m = margin (this, this.X, this.Y);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationSVM} {@var{e} =} resubEdge (@var{obj})
    ##
    ## Classification edge of the model on its own training data.
    ##
    ## @code{@var{e} = resubEdge (@var{obj})} is @code{edge} applied to the
    ## observations the model was fitted on, the mean of @code{resubMargin}.
    ##
    ## @end deftypefn
    function e = resubEdge (this)
      e = edge (this, this.X, this.Y);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationSVM object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## ClassificationSVM object into an Octave binary file, the name of which is
    ## specified in @var{filename}, along with an extra variable, which defines
    ## the type classification object these variables constitute.  Use
    ## @code{loadmodel} in order to load a classification object into Octave's
    ## workspace.
    ##
    ## @seealso{loadmodel, fitcsvm, ClassificationSVM}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error ("ClassificationSVM.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error ("ClassificationSVM.savemodel: FNAME must be a character vector.");
      endif
      ## Generate variable for class name
      classdef_name = 'ClassificationSVM';

      ## Create variables from model properties
      X = this.X;
      Y = this.Y;
      NumObservations     = this.NumObservations;
      RowsUsed            = this.RowsUsed;
      NumPredictors       = this.NumPredictors;
      PredictorNames      = this.PredictorNames;
      ResponseName        = this.ResponseName;
      ClassNames          = this.ClassNames;
      ScoreTransform      = this.ScoreTransform;
      Sigma               = this.Sigma;
      BinEdges        = this.BinEdges;
      Mu                  = this.Mu;
      ModelParameters     = this.ModelParameters;
      KernelParameters    = this.KernelParameters;
      BoxConstraints      = this.BoxConstraints;
      OutlierFraction     = this.OutlierFraction;
      Nu                  = this.Nu;
      Model               = this.Model;
      Alpha               = this.Alpha;
      Beta                = this.Beta;
      Bias                = this.Bias;
      IsSupportVector     = this.IsSupportVector;
      SupportVectorLabels = this.SupportVectorLabels;
      SupportVectors      = this.SupportVectors;
      W                   = this.W;
      Prior               = this.Prior;
      Cost                = this.Cost;
      CategoricalPredictors  = this.CategoricalPredictors;
      ExpandedPredictorNames = this.ExpandedPredictorNames;

      ## Save classdef name and all model properties as individual variables
      STfun          = this.STfun;

      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'NumObservations', ...
            'RowsUsed', 'NumPredictors', 'PredictorNames', 'BinEdges', ...
            'ResponseName', ...
            'ClassNames', 'ScoreTransform', 'Sigma', 'Mu',  ...
            'ModelParameters', 'Model', 'Alpha', 'Beta', 'Bias', ...
            'IsSupportVector', 'SupportVectorLabels', 'SupportVectors', ...
            'W', 'Prior', 'Cost', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'KernelParameters', ...
            'BoxConstraints', 'OutlierFraction', 'Nu', 'STfun');
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationSVM object
      mdl = ClassificationSVM (1, 1);

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
          error ("ClassificationSVM.load_model: invalid model in '%s'.", filename)
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

%!demo
%! ## Create a Support Vector Machine classifier and determine margin for test
%! ## data.
%! load fisheriris
%! rng (1);  ## For reproducibility
%!
%! ## Select indices of the non-setosa species
%! inds = ! strcmp (species, 'setosa');
%!
%!  ## Select features and labels for non-setosa species
%! X = meas(inds, 3:4);
%! Y = grp2idx (species(inds));
%!
%! ##  Convert labels to +1 and -1
%! unique_classes = unique (Y);
%! Y(Y == unique_classes(1)) = -1;
%! Y(Y == unique_classes(2)) = 1;
%!
%! ## Partition data for training and testing
%! cv = cvpartition (Y, 'HoldOut', 0.15);
%! X_train = X(training(cv), :);
%! Y_train = Y(training(cv));
%! X_test = X(test (cv), :);
%! Y_test = Y(test (cv));
%!
%! ## Train the SVM model
%! CVSVMModel = fitcsvm (X_train, Y_train);
%!
%! ## Calculate margins
%! m = margin (CVSVMModel, X_test, Y_test);
%! disp (m);

%!demo
%! ## Create a Support Vector Machine classifier and determine loss for test
%! ## data.
%! load fisheriris
%! rng (1);  ## For reproducibility
%!
%!  ## Select indices of the non-setosa species
%! inds = ! strcmp (species, 'setosa');
%!
%!  ## Select features and labels for non-setosa species
%! X = meas(inds, 3:4);
%! Y = grp2idx (species(inds));
%!
%! ##  Convert labels to +1 and -1
%! unique_classes = unique (Y);
%! Y(Y == unique_classes(1)) = -1;
%! Y(Y == unique_classes(2)) = 1;
%!
%! ## Randomly partition the data into training and testing sets
%! cv = cvpartition (Y, 'HoldOut', 0.3); # 30% data for testing, 60% for training
%!
%! X_train = X(training(cv), :);
%! Y_train = Y(training(cv));
%!
%! X_test = X(test (cv), :);
%! Y_test = Y(test (cv));
%!
%! ## Train the SVM model
%! SVMModel = fitcsvm (X_train, Y_train);
%!
%! ## Calculate loss
%!
%! L = loss (SVMModel,X_test,Y_test,'LossFun','binodeviance')
%! L = loss (SVMModel,X_test,Y_test,'LossFun','classiferror')
%! L = loss (SVMModel,X_test,Y_test,'LossFun','exponential')
%! L = loss (SVMModel,X_test,Y_test,'LossFun','hinge')
%! L = loss (SVMModel,X_test,Y_test,'LossFun','logit')
%! L = loss (SVMModel,X_test,Y_test,'LossFun','quadratic')

## Test output of constructor
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1; 4, 5, 6; 7, 8, 9; ...
%! 3, 2, 1; 4, 5, 6; 7, 8, 9; 3, 2, 1; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [1; 2; 3; 4; 2; 3; 4; 2; 3; 4; 2; 3; 4];
%! a = ClassificationSVM (x, y, 'ClassNames', [1, 2]);
%! assert_equal (class (a), "ClassificationSVM");
%! m = logical ([1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0]');
%! assert_equal (a.RowsUsed, m);
%! assert_equal ({a.X, a.Y}, {x(m,:), y(m)})
%! assert_equal (a.NumObservations, 5)
%! assert_equal ({a.ResponseName, a.PredictorNames}, {'Y', {'x1', 'x2', 'x3'}})
%! assert_equal ({a.ClassNames, a.ModelParameters.SVMtype}, {[1; 2], 'c_svc'})
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! a = ClassificationSVM (x, y);
%! assert_equal (class (a), "ClassificationSVM");
%! assert_equal ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, 'linear'})
%! assert_equal (a.ModelParameters.BoxConstraint, 1)
%! assert_equal (a.ClassNames, [-1; 1])
%! assert_equal (a.ModelParameters.KernelOffset, 0)
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! a = ClassificationSVM (x, y, 'KernelFunction', 'rbf', 'BoxConstraint', 2, ...
%! 'KernelOffset', 2);
%! assert_equal (class (a), "ClassificationSVM");
%! assert_equal ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, 'rbf'})
%! assert_equal (a.ModelParameters.BoxConstraint, 2)
%! assert_equal (a.ModelParameters.KernelOffset, 2)
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! a = ClassificationSVM (x, y, 'KernelFunction', 'polynomial', ...
%! 'PolynomialOrder', 3);
%! assert_equal (class (a), "ClassificationSVM");
%! assert_equal ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, 'polynomial'})
%! assert_equal (a.ModelParameters.PolynomialOrder, 3)

## RowsUsed is the logical mask its documentation describes, so the documented
## use of it works.  It was a double, for which a 0 is not a subscript at all.
%!test
%! randn ('seed', 42);
%! lab = double (randn (40, 1) > 0) + 1;
%! X = [randn(40, 3); NaN, 1, 1];
%! Y = [lab; 1];
%! Mdl = fitcsvm (X, Y);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (rows (Mdl.X), Mdl.NumObservations);

## resubPredict answers about the rows it was trained on, all their columns.
## Selecting them with a bare mask took the first column only, and the model
## answered about data it was never given, without complaint.
%!test
%! randn ('seed', 42);
%! lab = double (randn (40, 1) > 0) + 1;
%! X = [randn(40, 3); NaN, 1, 1];
%! Y = [lab; 1];
%! Mdl = fitcsvm (X, Y);
%! assert_equal (Mdl.NumObservations, 41);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (resubPredict (Mdl), predict (Mdl, Mdl.X));

## A model saved before RowsUsed became a mask still reads back as one.
%!test
%! randn ('seed', 42);
%! X = randn (31, 2);
%! Y = [double(randn (30, 1) > 0) + 1; NaN];
%! Mdl = fitcsvm (X, Y);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! d = load (fname);
%! d.RowsUsed = double (d.RowsUsed);
%! save ('-binary', fname, '-struct', 'd');
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2.RowsUsed), 'logical');
%! assert_equal (rows (M2.X(M2.RowsUsed, :)), M2.NumObservations);

## Prior defaults to the class frequencies and Cost to zero on the diagonal
## and one elsewhere, both measured against R2024a.
%!test
%! load fisheriris
%! Yb = strcmp (species, 'setosa');
%! Mdl = fitcsvm (meas, Yb);
%! assert_equal (Mdl.Prior, [2/3, 1/3], 1e-12);
%! assert_equal (Mdl.Cost, [0, 1; 1, 0]);
%! Yu = [repmat({'a'}, 120, 1); repmat({'b'}, 30, 1)];
%! M2 = fitcsvm (meas, Yu);
%! assert_equal (M2.Prior, [0.8, 0.2], 1e-12);
%! assert_equal (fitcsvm (meas, Yu, 'Prior', 'uniform').Prior, [0.5, 0.5]);
%! assert_equal (fitcsvm (meas, Yu, 'Prior', [3, 1]).Prior, [0.75, 0.25]);

## Prior and Cost reach the fit, as they do in MATLAB, through the per-class
## box constraint.  Asking for the defaults explicitly changes nothing.
%!test
%! load fisheriris
%! Yu = [repmat({'a'}, 120, 1); repmat({'b'}, 30, 1)];
%! base = fitcsvm (meas, Yu);
%! assert_equal (isequal (fitcsvm (meas, Yu, 'Prior', [0.5, 0.5]).Alpha, ...
%!                        base.Alpha), false);
%! assert_equal (isequal (fitcsvm (meas, Yu, 'Cost', [0, 5; 1, 0]).Alpha, ...
%!                        base.Alpha), false);
%! same = fitcsvm (meas, Yu, 'Prior', 'empirical', 'Cost', [0, 1; 1, 0]);
%! assert_equal (same.Alpha, base.Alpha);
%! assert_equal (same.Bias, base.Bias);

## The cost-aware losses agree with R2024a, which returns 0.01 for all three
## on this fixture and 0.04 once an error is charged four times over.
%!test
%! load fisheriris
%! X = meas(51:150,:);
%! Yb = strcmp (species(51:150), 'versicolor');
%! Ypm = ones (100, 1);
%! Ypm(Yb) = -1;
%! Mdl = fitcsvm (X, Yb, 'KernelFunction', 'linear');
%! assert_equal (loss (Mdl, X, Ypm, 'LossFun', 'classiferror'), 0.01, 1e-12);
%! assert_equal (loss (Mdl, X, Ypm, 'LossFun', 'classifcost'), 0.01, 1e-12);
%! assert_equal (loss (Mdl, X, Ypm, 'LossFun', 'mincost'), 0.01, 1e-12);
%! Mc = fitcsvm (X, Yb, 'KernelFunction', 'linear', 'Cost', [0, 4; 1, 0]);
%! assert_equal (loss (Mc, X, Ypm, 'LossFun', 'classifcost'), 0.04, 1e-12);
%! assert_equal (loss (Mc, X, Ypm, 'LossFun', 'mincost'), 0.04, 1e-12);

## resubLoss is the loss of the model on its own training rows, whatever type
## the labels are.  It formed the margin from the labels as given, so a 1/2
## coded response scored 0.49 where the answer is 0.01, and a logical one did
## not reach LIBSVM at all.
%!test
%! load fisheriris
%! X = meas(51:150,:);
%! Yb = strcmp (species(51:150), 'versicolor');
%! assert_equal (resubLoss (fitcsvm (X, Yb, 'KernelFunction', 'linear')), ...
%!               0.01, 1e-12);
%! assert_equal (resubLoss (fitcsvm (X, double (Yb) + 1, ...
%!                                   'KernelFunction', 'linear')), 0.01, 1e-12);
%! assert_equal (resubLoss (fitcsvm (X, species(51:150), ...
%!                                   'KernelFunction', 'linear')), 0.01, 1e-12);

## The compact model carries both, and its cost-aware losses agree with the
## model it was compacted from.
%!test
%! load fisheriris
%! X = meas(51:150,:);
%! Yb = strcmp (species(51:150), 'versicolor');
%! Ypm = ones (100, 1);
%! Ypm(Yb) = -1;
%! Mdl = fitcsvm (X, Yb, 'KernelFunction', 'linear', 'Cost', [0, 4; 1, 0]);
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.Prior, Mdl.Prior);
%! assert_equal (CMdl.Cost, Mdl.Cost);
%! for f = {'classiferror', 'classifcost', 'mincost'}
%!   assert_equal (loss (CMdl, X, Ypm, 'LossFun', f{1}), ...
%!                 loss (Mdl, X, Ypm, 'LossFun', f{1}), 1e-12);
%! endfor

## discardSupportVectors empties what R2024a empties and keeps what it
## keeps: Alpha and the support vectors go, Beta, Bias and IsSupportVector
## stay, and the class is unchanged.
%!test
%! load fisheriris
%! keep = ! strcmp (species, "setosa");
%! X = meas(keep,:); y = species(keep);
%! Mdl = fitcsvm (X, y, "KernelFunction", "linear");
%! D = discardSupportVectors (Mdl);
%! assert_equal (class (D), "ClassificationSVM");
%! assert_equal (isempty (D.Alpha), true);
%! assert_equal (isempty (D.SupportVectors), true);
%! assert_equal (isempty (D.SupportVectorLabels), true);
%! assert_equal (D.Beta, Mdl.Beta);
%! assert_equal (D.Bias, Mdl.Bias);
%! assert_equal (D.IsSupportVector, Mdl.IsSupportVector);

## A linear decision needs only Beta and Bias, so the model predicts what it
## predicted before.
%!test
%! load fisheriris
%! keep = ! strcmp (species, "setosa");
%! X = meas(keep,:); y = species(keep);
%! Mdl = fitcsvm (X, y, "KernelFunction", "linear");
%! D = discardSupportVectors (Mdl);
%! assert_equal (predict (D, X), predict (Mdl, X), 1e-10);

## The saving is real rather than cosmetic: the engine keeps its own copy of
## the support vectors, and it collapses to the one vector that decides a
## linear model.  Emptying the properties alone would free nothing.
%!test
%! load fisheriris
%! keep = ! strcmp (species, "setosa");
%! X = meas(keep,:); y = species(keep);
%! Mdl = fitcsvm (X, y, "KernelFunction", "linear");
%! D = discardSupportVectors (Mdl);
%! assert_equal (rows (Mdl.Model.SVs) > 1, true);
%! assert_equal (rows (D.Model.SVs), 1);
%! assert_equal (predict (discardSupportVectors (D), X), predict (D, X));

%!error<ClassificationSVM.discardSupportVectors: you cannot discard support vectors for a non-linear kernel.> ...
%! load fisheriris
%! keep = ! strcmp (species, "setosa");
%! X = meas(keep,:); y = species(keep);
%! discardSupportVectors (fitcsvm (X, y, "KernelFunction", "rbf"))

%!error<ClassificationSVM: 'Prior' must be a non-negative numeric vector, 'empirical' or 'uniform'.> ...
%! fitcsvm (ones (10,2), [1;1;1;1;1;2;2;2;2;2], 'Prior', 'nope')
%!error<ClassificationSVM: 'Cost' must be a non-negative square matrix.> ...
%! fitcsvm (ones (10,2), [1;1;1;1;1;2;2;2;2;2], 'Cost', [0, 1])
%!error<ClassificationSVM: 'Prior' must have one entry per class.> ...
%! fitcsvm (ones (10,2), [1;1;1;1;1;2;2;2;2;2], 'Prior', [1, 1, 1])
%!error<ClassificationSVM: the number of rows and columns in 'Cost' must correspond to the classes in Y.> ...
%! fitcsvm (ones (10,2), [1;1;1;1;1;2;2;2;2;2], 'Cost', eye (3))

## The model reports the observation weights and the expanded predictor
## names MATLAB reports, W normalized to sum to one.
%!test
%! load fisheriris
%! Yb = strcmp (species, 'setosa');
%! Mdl = fitcsvm (meas, Yb);
%! assert_equal (size (Mdl.W), [150, 1]);
%! assert_equal (sum (Mdl.W), 1, 1e-12);
%! assert_equal (Mdl.W, ones (150, 1) / 150, 1e-12);
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (Mdl.ExpandedPredictorNames, Mdl.PredictorNames);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (M2.W, Mdl.W);
%! assert_equal (M2.ExpandedPredictorNames, Mdl.ExpandedPredictorNames);

## Standardize fits on the scale it predicts on.  A model fitted on raw data
## and asked about standardized data is not merely worse, it is wrong: on a
## problem whose second predictor is a thousand times the first it answered at
## chance, where the same data standardized by hand reaches 0.9.
%!test
%! rand ('seed', 42);
%! randn ('seed', 42);
%! n = 80;
%! A = [randn(n,1) + 1, (randn(n,1) + 1) * 1000];
%! B = [randn(n,1) - 1, (randn(n,1) - 1) * 1000];
%! X = [A; B];
%! Y = [ones(n,1); 2 * ones(n,1)];
%! Mdl = fitcsvm (X, Y, 'Standardize', true, 'KernelFunction', 'rbf');
%! assert_equal (mean (resubPredict (Mdl) == Y) > 0.85, true);
%! assert_equal (predict (Mdl, X), resubPredict (Mdl));
%! Mu = mean (X, 1);
%! Sg = std (X, [], 1);
%! byhand = fitcsvm ((X - Mu) ./ Sg, Y, 'Standardize', false, ...
%!                   'KernelFunction', 'rbf');
%! assert_equal (mean (resubPredict (Mdl) == Y), ...
%!               mean (resubPredict (byhand) == Y));

## Test input validation for constructor
%!error<ClassificationSVM: too few input arguments.> ClassificationSVM ()
%!error<ClassificationSVM: too few input arguments.> ...
%! ClassificationSVM (ones (10,2))
%!error<ClassificationSVM: number of rows in X and Y must be equal.> ...
%! ClassificationSVM (ones (10,2), ones (5,1))
%!error<ClassificationSVM: 'Standardize' must be either true or false.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'Standardize', 'a')
%!error<ClassificationSVM: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'PredictorNames', ['x1';'x2'])
%!error<ClassificationSVM: 'PredictorNames' must have the same number of columns as X.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'PredictorNames', {'x1','x2','x3'})
%!error<ClassificationSVM: 'ResponseName' must be a character vector.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'ResponseName', {'Y'})
%!error<ClassificationSVM: 'ResponseName' must be a character vector.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'ResponseName', 21)
%!error<ClassificationSVM: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'ClassNames', @(x)x)
%!error<ClassificationSVM: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'ClassNames', {1})
%!error<ClassificationSVM: not all 'ClassNames' are present in Y.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'ClassNames', [1, 2])
%!error<ClassificationSVM: not all 'ClassNames' are present in Y.> ...
%! ClassificationSVM (ones (5,2), ['a';'b';'a';'a';'b'], 'ClassNames', ['a';'c'])
%!error<ClassificationSVM: not all 'ClassNames' are present in Y.> ...
%! ClassificationSVM (ones (5,2), {'a';'b';'a';'a';'b'}, 'ClassNames', {'a','c'})
%!error<ClassificationSVM: not all 'ClassNames' are present in Y.> ...
%! ClassificationSVM (ones (10,2), logical (ones (10,1)), 'ClassNames', [true, false])
%!error<ClassificationSVM: 'SVMtype' must be 'c_svc', 'nu_svc', or 'one_class_svm'.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'svmtype', 123)
%!error<ClassificationSVM: 'SVMtype' must be 'c_svc', 'nu_svc', or 'one_class_svm'.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'svmtype', 'some_type')
%!error<ClassificationSVM: 'OutlierFraction' must be a positive scalar in the range 0 =< OutlierFraction < 1.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'OutlierFraction', -1)
%!error<ClassificationSVM: 'KernelFunction' must be a character vector.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'KernelFunction', 123)
%!error<ClassificationSVM: unsupported Kernel function.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'KernelFunction', 'fcn')
%!error<ClassificationSVM: 'PolynomialOrder' must be a positive integer.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'PolynomialOrder', -1)
%!error<ClassificationSVM: 'PolynomialOrder' must be a positive integer.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'PolynomialOrder', 0.5)
%!error<ClassificationSVM: 'PolynomialOrder' must be a positive integer.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'PolynomialOrder', [1,2])
%!error<ClassificationSVM: 'KernelScale' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'KernelScale', -1)
%!error<ClassificationSVM: 'KernelScale' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'KernelScale', 0)
%!error<ClassificationSVM: 'KernelScale' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'KernelScale', [1, 2])
%!error<ClassificationSVM: 'KernelScale' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'KernelScale', 'invalid')
%!error<ClassificationSVM: 'KernelOffset' must be a non-negative scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'KernelOffset', -1)
%!error<ClassificationSVM: 'KernelOffset' must be a non-negative scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'KernelOffset', [1,2])
%!error<ClassificationSVM: 'BoxConstraint' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'BoxConstraint', -1)
%!error<ClassificationSVM: 'BoxConstraint' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'BoxConstraint', 0)
%!error<ClassificationSVM: 'BoxConstraint' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'BoxConstraint', [1, 2])
%!error<ClassificationSVM: 'BoxConstraint' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'BoxConstraint', 'invalid')
%!error<ClassificationSVM: 'Nu' must be a positive scalar in the range 0 < Nu <= 1.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'nu', -0.5)
%!error<ClassificationSVM: 'Nu' must be a positive scalar in the range 0 < Nu <= 1.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'nu', 0)
%!error<ClassificationSVM: 'Nu' must be a positive scalar in the range 0 < Nu <= 1.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'nu', 1.5)
%!error<ClassificationSVM: 'CacheSize' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'CacheSize', -1)
%!error<ClassificationSVM: 'CacheSize' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'CacheSize', [1,2])
%!error<ClassificationSVM: 'Tolerance' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'Tolerance', -0.1)
%!error<ClassificationSVM: 'Tolerance' must be a positive scalar.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'Tolerance', [0.1,0.2])
%!error<ClassificationSVM: 'Shrinking' must be either 0 or 1.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'shrinking', 2)
%!error<ClassificationSVM: 'Shrinking' must be either 0 or 1.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'shrinking', -1)
%!error<ClassificationSVM: 'Shrinking' must be either 0 or 1.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'shrinking', [1 0])
%!error<ClassificationSVM: invalid parameter name in optional pair arguments.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'invalid_name', 'c_svc')
%!error<ClassificationSVM: cannot train a binary problem with only one class available.> ...
%! ClassificationSVM (ones (10,2), ones (10,1), 'SVMtype', 'c_svc')
%!error<ClassificationSVM: can only be used for one-class or two-class learning.> ...
%! ClassificationSVM (ones (10,2), [1;1;1;1;2;2;2;2;3;3])
%!error<ClassificationSVM: invalid values in X.> ...
%! ClassificationSVM ([ones(9,2);2,Inf], ones (10,1))

## Test output for predict method
%!shared x, y, x_train, x_test, y_train, y_test, objST
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcsvm (x, y, 'KernelFunction', 'rbf', 'Tolerance', 1e-7);
%! assert_equal (isempty (obj.Beta), true)
%! assert_equal (sum (obj.IsSupportVector), numel (obj.Alpha))
%! [label, score] = predict (obj, xc);
%! assert_equal (label, [1; 2; 2]);
%! assert_equal (score(:,1), [0.99285; -0.080296; -0.93694], 2e-5);
%! assert_equal (score(:,1), -score(:,2), eps)
%!test
%! obj = fitcsvm (x, y);
%! assert_equal (obj.Beta, [2.182926829268275; 2.253658536585344], 1e-5)
%! assert_equal (sum (obj.IsSupportVector), numel (obj.Alpha))
%! assert_equal (numel (obj.Alpha), 24)
%! assert_equal (obj.Bias, -14.415, 1e-3)
%! xc = [min(x); mean(x); max(x)];
%! label = predict (obj, xc);
%! assert_equal (label, [1; 2; 2]);

## Values below are R2024a's, measured 2026-08-17.
%!test
%! ## a linear kernel has a primal representation, one coefficient per predictor
%! obj = fitcsvm (x, y);
%! assert_equal (size (obj.Beta), [2, 1]);
%! assert_equal (obj.Beta, [2.182926829268275; 2.253658536585344], 1e-5);
%! assert_equal (obj.Beta, ...
%!               obj.SupportVectors' * (obj.Alpha .* obj.SupportVectorLabels));
%!test
%! ## a nonlinear kernel has none, but keeps the dual coefficients
%! obj = fitcsvm (x, y, 'KernelFunction', 'rbf');
%! assert_equal (isempty (obj.Beta), true);
%! assert_equal (numel (obj.Alpha), sum (obj.IsSupportVector));
%!test
%! ## the dual coefficients are magnitudes; their class is in the labels
%! obj = fitcsvm (x, y);
%! assert_equal (any (obj.Alpha < 0), false);
%! assert_equal (max (obj.Alpha) <= 1, true);
%! assert_equal (size (obj.SupportVectorLabels), [24, 1]);
%! assert_equal (unique (obj.SupportVectorLabels)', [-1, 1]);
%!test
%! ## the support vector indicator is logical, one entry per observation
%! obj = fitcsvm (x, y);
%! assert_equal (class (obj.IsSupportVector), 'logical');
%! assert_equal (size (obj.IsSupportVector), [100, 1]);
%! assert_equal (sum (obj.IsSupportVector), 24);

## A single-observation query used to corrupt the heap and abort the
## interpreter; each row on its own must agree with the batch answer.
%!test
%! obj = fitcsvm (x, y);
%! xc = [min(x); mean(x); max(x)];
%! batch = predict (obj, xc);
%! for i = 1:rows (xc)
%!   assert_equal (predict (obj, xc(i,:)), batch(i));
%! endfor
%!test
%! obj = fitcsvm (x, y, 'KernelFunction', 'rbf', 'Tolerance', 1e-7);
%! xc = [min(x); mean(x); max(x)];
%! [bl, bs] = predict (obj, xc);
%! [l1, s1] = predict (obj, xc(1,:));
%! assert_equal (l1, bl(1));
%! assert_equal (size (l1), [1, 1]);
%! assert_equal (size (s1), [1, 2]);
%! assert_equal (s1, bs(1,:), 2e-5);
%!test
%! obj = compact (fitcsvm (x, y));
%! xc = [min(x); mean(x); max(x)];
%! batch = predict (obj, xc);
%! assert_equal (predict (obj, xc(1,:)), batch(1));

## Test input validation for predict method
%!error<ClassificationSVM.predict: too few input arguments.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)))
%!error<ClassificationSVM.predict: XC is empty.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), [])
%!error<ClassificationSVM.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), 1)
%!test
%! objST = fitcsvm (x, y);
%!error<ClassificationSVM.subsasgn: 'ScoreTransform' must be a 'function_handle' object.> ...
%! objST.ScoreTransform = 'a';
%! [labels, scores] = predict (objST, x);

## Test input validation for resubPredict method
%! [labels, scores] = resubPredict (objST);

## Test output for margin method
%!test
%! rand ('seed', 1);
%! CVSVMModel = fitcsvm (x, y, 'KernelFunction', 'rbf', 'HoldOut', 0.15, ...
%!                       'Tolerance', 1e-7);
%! obj = CVSVMModel.Trained{1};
%! testInds = test (CVSVMModel.Partition);
%! ## Every one of these fifteen is classified correctly, so every margin is
%! ## positive.  They used to read -4.0000 downwards for the second class:
%! ## the margin was formed from the response as given, so a 1/2 coding
%! ## scaled that class by four instead of negating it, and the model looked
%! ## as though it misclassified every observation of it.
%! expected_margin = [2.0000;  0.8579;  1.6690;  3.4141;  3.4552; ...
%!                    2.6605;  3.5251;  2.0000;  3.1705;  3.2256; ...
%!                    1.5266;  3.7527;  0.8350;  2.8113;  3.6820];
%! computed_margin = margin (obj, x(testInds,:), y(testInds,:));
%! assert_equal (computed_margin, expected_margin, 1e-4);
%! assert (all (computed_margin > 0));

## Test input validation for margin method
%!error<ClassificationSVM.margin: too few input arguments.> ...
%! margin (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)))
%!error<ClassificationSVM.margin: too few input arguments.> ...
%! margin (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2))
%!error<ClassificationSVM.margin: X is empty.> ...
%! margin (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), [], zeros (2))
%!error<ClassificationSVM.margin: X must have the same number of predictors as the trained model.> ...
%! margin (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), 1, zeros (2))
%!error<ClassificationSVM.margin: Y is empty.> ...
%! margin (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), [])
%!error<ClassificationSVM.margin: Y must have the same number of rows as X.> ...
%! margin (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), 1)

## Test output for loss method
%!test
%! rand ('seed', 1);
%! CVSVMModel = fitcsvm (x, y, 'KernelFunction', 'rbf', 'HoldOut', 0.15);
%! obj = CVSVMModel.Trained{1};
%! testInds = test (CVSVMModel.Partition);
%! L1 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'binodeviance');
%! L2 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'classiferror');
%! L3 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'exponential');
%! L4 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'hinge');
%! L5 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'logit');
%! L6 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'quadratic');
%! assert_equal (L1, 2.8711, 1e-4);
%! assert_equal (L2, 0.5333, 1e-4);
%! assert_equal (L3, 10.9685, 1e-4);
%! assert_equal (L4, 1.9827, 1e-4);
%! assert_equal (L5, 1.5849, 1e-4);
%! assert_equal (L6, 7.6739, 1e-4);

## Test input validation for loss method
%!error<ClassificationSVM.loss: too few input arguments.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)))
%!error<ClassificationSVM.loss: too few input arguments.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2))
%!error<ClassificationSVM.loss: Name-Value arguments must be in pairs.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), 'LossFun')
%!error<ClassificationSVM.loss: X is empty.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), [], zeros (2))
%!error<ClassificationSVM.loss: X must have the same number of predictors as the trained model.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), 1, zeros (2))
%!error<ClassificationSVM.loss: Y is empty.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), [])
%!error<ClassificationSVM.loss: Y must have the same number of rows as X.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), 1)
%!error<ClassificationSVM.loss: 'LossFun' must be a character vector.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), 'LossFun', 1)
%!error<ClassificationSVM.loss: unsupported Loss function.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), 'LossFun', 'some')
%!error<ClassificationSVM.loss: 'Weights' must be a numeric vector.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), 'Weights', ['a','b'])
%!error<ClassificationSVM.loss: 'Weights' must be a numeric vector.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), 'Weights', 'a')
%!error<ClassificationSVM.loss: size of 'Weights' must be equal to the number> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), 'Weights', [1,2,3])
%!error<ClassificationSVM.loss: size of 'Weights' must be equal to the number> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), 'Weights', 3)
%!error<ClassificationSVM.loss: invalid parameter name in optional pair arg> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), 'some', 'some')

## Test input validation for resubLoss method
%!error<ClassificationSVM.resubLoss: Name-Value arguments must be in pairs.> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), 'LossFun')
%!error<ClassificationSVM.resubLoss: 'LossFun' must be a character vector.> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), 'LossFun', 1)
%!error<ClassificationSVM.resubLoss: unsupported Loss function.> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), 'LossFun', 'some')
%!error<ClassificationSVM.resubLoss: 'Weights' must be a numeric vector.> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), 'Weights', ['a','b'])
%!error<ClassificationSVM.resubLoss: 'Weights' must be a numeric vector.> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), 'Weights', 'a')
%!error<ClassificationSVM.resubLoss: size of 'Weights' must be equal to the n> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), 'Weights', [1,2,3])
%!error<ClassificationSVM.resubLoss: size of 'Weights' must be equal to the n> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), 'Weights', 3)
%!error<ClassificationSVM.resubLoss: invalid parameter name in optional pai> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), 'some', 'some')

## Test output for crossval method
%!test
%! SVMModel = fitcsvm (x, y);
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (SVMModel, 'KFold', 5);
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (CVMdl.KFold == 5, true)
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationSVM")
%! assert_equal (CVMdl.CrossValidatedModel, "SVM")
%!test
%! obj = fitcsvm (x, y);
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (obj, 'HoldOut', 0.2);
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationSVM")
%! assert_equal (CVMdl.CrossValidatedModel, "SVM")
%!test
%! obj = fitcsvm (x, y);
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (obj, 'LeaveOut', 'on');
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationSVM")
%! assert_equal (CVMdl.CrossValidatedModel, "SVM")

## Test input validation for crossval method
%!error<ClassificationSVM.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'KFold')
%!error<ClassificationSVM.crossval: specify only one of the optional Name-Value paired arguments.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), ...
%! 'KFold', 5, 'leaveout', 'on')
%!error<ClassificationSVM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'KFold', 'a')
%!error<ClassificationSVM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'KFold', 1)
%!error<ClassificationSVM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'KFold', -1)
%!error<ClassificationSVM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'KFold', 11.5)
%!error<ClassificationSVM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'KFold', [1,2])
%!error<ClassificationSVM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'Holdout', 'a')
%!error<ClassificationSVM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'Holdout', 11.5)
%!error<ClassificationSVM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'Holdout', -1)
%!error<ClassificationSVM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'Holdout', 0)
%!error<ClassificationSVM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'Holdout', 1)
%!error<ClassificationSVM.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'Leaveout', 1)
%!error<ClassificationSVM.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'CVPartition', 1)
%!error<ClassificationSVM.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'CVPartition', 'a')
%!error<ClassificationSVM.crossval: invalid parameter name in optional paired arguments> ...
%! crossval (ClassificationSVM (ones (40,2),randi ([1, 2], 40, 1)), 'some', 'some')
%!error <ClassificationSVM.savemodel: too few input arguments.> ...
%! savemodel (ClassificationSVM ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]))
%!error <ClassificationSVM.savemodel: FNAME must be a character vector.> ...
%! savemodel (ClassificationSVM ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), 1)
%!error <ClassificationSVM.savemodel: FNAME must be a character vector.> ...
%! savemodel (ClassificationSVM ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), ['ab'; 'cd'])

## RowsUsed is empty when every observation was used.
%!test
%! load fisheriris
%! X = meas(1:100,:);
%! Y = grp2idx (species(1:100));
%! Mdl = fitcsvm (X, Y);
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
%! Mdl = fitcsvm (X, Y);
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
%! Mdl = fitcsvm (X, Y);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (Mdl.NumObservations, 100);
%! assert_equal (rows (Mdl.X), 100);
%! assert_equal (sum (isnan (Mdl.X(:))), 1);

## Mu and Sigma are empty unless the predictors are standardized.
%!test
%! load fisheriris
%! Mdl = fitcsvm (meas(1:100,:), species(1:100));
%! assert_equal (Mdl.Mu, []);
%! assert_equal (Mdl.Sigma, []);

## Standardizing weights the complete observations by each class's original
## share of the observation weight.  Values from MATLAB R2024a.
%!test
%! load fisheriris
%! X = meas(1:100,:);
%! X(3,2) = NaN; X(17,4) = NaN;
%! Mdl = fitcsvm (X, species(1:100), 'Standardize', true);
%! assert_equal (Mdl.Mu, [5.4700833333333314, 3.0964583333333331, ...
%!                        2.864374999999999, 0.78487499999999988], 1e-13);
%! assert_equal (Mdl.Sigma, [0.64239212471819018, 0.47709034264660688, ...
%!                           1.4464268842305728, 0.56625850081877471], 1e-13);

## The prior reweights the observations of each class.  Values from R2024a.
%!test
%! load fisheriris
%! i2 = [1:50, 51:80];
%! Mdl = fitcsvm (meas(i2,:), species(i2));
%! assert_equal (Mdl.Prior, [0.625, 0.375], 1e-14);
%! assert_equal (Mdl.W(1), 0.0125, 1e-14);
%! Mdl = fitcsvm (meas(i2,:), species(i2), 'Prior', 'uniform');
%! assert_equal (Mdl.W(1), 0.01, 1e-14);
%! assert_equal (Mdl.W(51), 1/60, 1e-14);

## Neither Prior nor Cost may be assigned after the fit.  The refusal comes
## from the property attributes, so the message is core Octave's and is not
## pinned here.
%!error ...
%! load fisheriris; ...
%! Mdl = fitcsvm (meas(1:80,:), species(1:80)); ...
%! Mdl.Prior = [0.5, 0.5];
%!error ...
%! load fisheriris; ...
%! Mdl = fitcsvm (meas(1:80,:), species(1:80)); ...
%! Mdl.Cost = [0, 2; 1, 0];

## A fitted model survives savemodel and loadmodel: the properties come
## back as they were and it predicts the same.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcsvm (meas(inds,:), species(inds));
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'ClassificationSVM');
%! assert_equal (M2.NumObservations, Mdl.NumObservations);
%! assert_equal (M2.PredictorNames, Mdl.PredictorNames);
%! assert_equal (class (M2.ScoreTransform), class (Mdl.ScoreTransform));
%! assert_equal (predict (M2, meas(1:5,:)), predict (Mdl, meas(1:5,:)));

## edge, resubMargin and resubEdge.  The absolute value is not pinned to the
## oracle here: this class fits through LIBSVM where MATLAB uses SMO, and the
## two disagree on the support vector set, so the edge reads 3.6942705869
## against MATLAB's 3.6935623843 on the iris pair.  What is pinned is what the
## methods mean.

%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Y = species(inds);
%! Mdl = fitcsvm (X, Y);
%! assert_equal (edge (Mdl, X, Y), mean (margin (Mdl, X, Y)), 1e-12);
%! assert_equal (resubEdge (Mdl), edge (Mdl, X, Y), 1e-12);
%! assert_equal (resubMargin (Mdl), margin (Mdl, X, Y), 1e-12);

## margin takes the class labels, as it documents, and the +1/-1 coding the
## solver works in gives the same answer.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Y = species(inds);
%! Mdl = fitcsvm (X, Y);
%! Ypm = ones (100, 1);
%! Ypm(strcmp (Y, Mdl.ClassNames{2})) = -1;
%! assert_equal (margin (Mdl, X, Y), margin (Mdl, X, Ypm), 1e-12);

## The compact model answers the same edge as the full one.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! X = meas(inds,:);
%! Y = species(inds);
%! Mdl = fitcsvm (X, Y);
%! assert_equal (edge (compact (Mdl), X, Y), edge (Mdl, X, Y), 1e-12);

%!error<ClassificationSVM.edge: too few input arguments.> ...
%! load fisheriris; ...
%! inds = ! strcmp (species, 'virginica'); ...
%! edge (fitcsvm (meas(inds,:), species(inds)), meas(inds,:))

## BinEdges is an empty cell, which is what MATLAB reports for this
## learner as well: it fits the predictors as they are.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! Mdl = fitcsvm (meas(inds,:), species(inds));
%! assert_equal (class (Mdl.BinEdges), 'cell');
%! assert_equal (Mdl.BinEdges, {});

## KernelParameters, BoxConstraints, OutlierFraction and Nu, all measured on
## MATLAB R2024a.  A radial basis kernel reports MATLAB's name for it.
%!test
%! load fisheriris
%! b = ismember (species, {'setosa', 'versicolor'});
%! Mdl = fitcsvm (meas(b,:), species(b));
%! assert_equal (Mdl.KernelParameters, struct ('Function', 'linear', ...
%!                                             'Scale', 1));
%! assert_equal (Mdl.BoxConstraints, ones (100, 1));
%! assert_equal (Mdl.OutlierFraction, 0);
%! assert_equal (Mdl.Nu, []);

%!test
%! load fisheriris
%! b = ismember (species, {'setosa', 'versicolor'});
%! Mdl = fitcsvm (meas(b,:), species(b), 'KernelFunction', 'rbf', ...
%!                'KernelScale', 2.5, 'BoxConstraint', 3);
%! assert_equal (Mdl.KernelParameters.Function, 'gaussian');
%! assert_equal (Mdl.KernelParameters.Scale, 2.5);
%! assert_equal (unique (Mdl.BoxConstraints), 3);

%!test
%! load fisheriris
%! b = ismember (species, {'setosa', 'versicolor'});
%! Mdl = fitcsvm (meas(b,:), species(b), 'KernelFunction', 'polynomial', ...
%!                'PolynomialOrder', 3);
%! assert_equal (fieldnames (Mdl.KernelParameters), ...
%!               {'Function'; 'Scale'; 'Order'});
%! assert_equal (Mdl.KernelParameters.Order, 3);

## A cost matrix reweights the classes, and the box constraints report what
## the fit actually solved: 1.6 and 0.4 against MATLAB's 1.6 and 0.4.
%!test
%! load fisheriris
%! b = ismember (species, {'setosa', 'versicolor'});
%! Mdl = fitcsvm (meas(b,:), species(b), 'Cost', [0, 4; 1, 0]);
%! assert_equal (Mdl.BoxConstraints(1), 1.6, 1e-12);
%! assert_equal (Mdl.BoxConstraints(51), 0.4, 1e-12);
%! assert_equal (mean (Mdl.BoxConstraints), 1, 1e-12);

%!test
%! load fisheriris
%! b = ismember (species, {'setosa', 'versicolor'});
%! Mdl = fitcsvm (meas(b,:), species(b), 'OutlierFraction', 0.05);
%! assert_equal (Mdl.OutlierFraction, 0.05);

%!test
%! load fisheriris
%! b = ismember (species, {'setosa', 'versicolor'});
%! Mdl = fitcsvm (meas(b,:), species(b));
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (M2.KernelParameters, Mdl.KernelParameters);
%! assert_equal (M2.BoxConstraints, Mdl.BoxConstraints);
