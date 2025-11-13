## Copyright (C) 2024-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef CompactClassificationSVM
## -*- texinfo -*-
## @deftp {statistics} CompactClassificationSVM
##
## Compact Support Vector Machine classification
##
## The @code{CompactClassificationSVM} class implements a compact version of a
## Support Vector Machine classifier object for one-class or two-class problems,
## which can predict responses for new data using the @code{predict} method.
##
## A @code{CompactClassificationSVM} object is a compact version of a support
## vector machine model, @code{ClassificationSVM}.  It does not include the
## training data resulting in a smaller classifier size, which can be used for
## making predictions from new data, but not for tasks such as cross validation.
## It can only be created from a @code{ClassificationSVM} model by using the
## @code{compact} object method.
##
## @seealso{ClassificationSVM}
## @end deftp

  properties (Access = public)
    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer value specifying the number of predictors in the
    ## training dataset used for training the SVM model.  This property is
    ## read-only.
    ##
    ## @end deftp
    NumPredictors       = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which they appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames      = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName        = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} ClassNames
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
    ## @deftp {CompactClassificationSVM} {property} ScoreTransform
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
    ## @multitable @columnfractions 0.2 0.05 0.75
    ## @headitem @var{Value} @tab @tab @var{Description}
    ## @item @qcode{"doublelogit"} @tab @tab @math{1 ./ (1 + exp .^ (-2 * x))}
    ## @item @qcode{"invlogit"} @tab @tab @math{1 ./ (1 + exp .^ (-x))}
    ## @item @qcode{"ismax"} @tab @tab Sets the score for the class with the
    ## largest score to 1, and for all other classes to 0
    ## @item @qcode{"logit"} @tab @tab @math{log (x ./ (1 - x))}
    ## @item @qcode{"none"} @tab @tab @math{x} (no transformation)
    ## @item @qcode{"identity"} @tab @tab @math{x} (no transformation)
    ## @item @qcode{"sign"} @tab @tab @math{-1 for x < 0, 0 for x = 0, 1 for x > 0}
    ## @item @qcode{"symmetric"} @tab @tab @math{2 * x + 1}
    ## @item @qcode{"symmetricismax"} @tab @tab Sets the score for the class
    ## with the largest score to 1, and for all other classes to -1
    ## @item @qcode{"symmetriclogit"} @tab @tab @math{2 ./ (1 + exp .^ (-x)) - 1}
    ## @end multitable
    ##
    ## @end deftp
    ScoreTransform      = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} Standardize
    ##
    ## Flag to standardize predictors
    ##
    ## A boolean flag indicating whether the data in @var{X} have been
    ## standardized prior to training.  This property is read-only.
    ##
    ## @end deftp
    Standardize         = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} Sigma
    ##
    ## Predictor standard deviations
    ##
    ## A numeric vector of the same length as the columns in @var{X} containing
    ## the standard deviations of predictor variables.  If the predictor
    ## variables have not been standardized, then @qcode{Sigma} is empty.
    ## This property is read-only.
    ##
    ## @end deftp
    Sigma               = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} Mu
    ##
    ## Predictor means
    ##
    ## A numeric vector of the same length as the columns in @var{X} containing
    ## the means of predictor variables.  If the predictor variables have not
    ## been standardized, then @qcode{Mu} is empty.  This property is read-only.
    ##
    ## @end deftp
    Mu                  = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} ModelParameters
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
    ## @deftp {CompactClassificationSVM} {property} Model
    ##
    ## Trained SVM model
    ##
    ## A structure containing the trained model in @qcode{'libsvm'} format.
    ## This property is read-only.
    ##
    ## @end deftp
    Model               = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} Alpha
    ##
    ## Trained classifier coefficients
    ##
    ## The coefficients of the trained SVM classifier specified as an @math{sx1}
    ## numeric vector, where @math{s} is the number of support vectors equal to
    ## @qcode{sum (obj.IsSupportVector)}.  If the SVM classifier was trained
    ## with a kernel function other than @qcode{'linear'}, then @qcode{Alpha} is
    ## empty.  This property is read-only.
    ##
    ## @end deftp
    Alpha               = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} Beta
    ##
    ## Linear predictor coefficients
    ##
    ## The linear predictor coefficients specified as an @math{sx1} numeric
    ## vector, where @math{s} is the number of support vectors equal to
    ## @qcode{sum (obj.IsSupportVector)}.  If the SVM classifier was trained
    ## with a @qcode{'linear'} kernel function, then @qcode{Beta} is empty.
    ## This property is read-only.
    ##
    ## @end deftp
    Beta                = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} Bias
    ##
    ## Bias term
    ##
    ## The bias term specified as a scalar.  This property is read-only.
    ##
    ## @end deftp
    Bias                = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} IsSupportVector
    ##
    ## Support vector indicator
    ##
    ## An @math{Nx1} logical vector that flags whether a corresponding
    ## observation in the predictor data matrix is a Support Vector.  @math{N}
    ## is the number of observations in the training data.  This property is
    ## read-only.
    ##
    ## @end deftp
    IsSupportVector     = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} SupportVectorLabels
    ##
    ## Support vector class labels
    ##
    ## The support vector class labels specified as an @math{sx1} numeric
    ## vector, where @math{s} is the number of support vectors equal to
    ## @qcode{sum (obj.IsSupportVector)}.  A value of +1 in
    ## @code{SupportVectorLabels} indicates that the corresponding support
    ## vector belongs to the positive class @qcode{(ClassNames@{2@})}.  A value
    ## of -1 indicates that the corresponding support vector belongs to the
    ## negative class @qcode{(ClassNames@{1@})}.  This property is read-only.
    ##
    ## @end deftp
    SupportVectorLabels = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationSVM} {property} SupportVectors
    ##
    ## Support vectors
    ##
    ## The support vectors of the trained SVM classifier specified an @math{sxp}
    ## numeric matrix, where @math{s} is the number of support vectors equal to
    ## @qcode{sum (obj.IsSupportVector)}, and @math{p} is the number of
    ## predictor variables in the predictor data.  This property is read-only.
    ##
    ## @end deftp
    SupportVectors      = [];
  endproperties

  methods (Hidden)

    ## constructor
    function this = CompactClassificationSVM (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), "ClassificationSVM"))
        error (strcat (["CompactClassificationSVM: invalid"], ...
                       [" classification object."]));
      endif

      ## Save properties to compact model
      this.NumPredictors         = Mdl.NumPredictors;
      this.PredictorNames        = Mdl.PredictorNames;
      this.ResponseName          = Mdl.ResponseName;
      this.ClassNames            = Mdl.ClassNames;

      this.ScoreTransform        = Mdl.ScoreTransform;

      this.Standardize           = Mdl.Standardize;
      this.Sigma                 = Mdl.Sigma;
      this.Mu                    = Mdl.Mu;

      this.ModelParameters       = Mdl.ModelParameters;
      this.Model                 = Mdl.Model;

      this.Alpha                 = Mdl.Alpha;
      this.Beta                  = Mdl.Beta;
      this.Bias                  = Mdl.Bias;
      this.IsSupportVector       = Mdl.IsSupportVector;
      this.SupportVectorLabels   = Mdl.SupportVectorLabels;
      this.SupportVectors        = Mdl.SupportVectors;

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationSVM} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationSVM} {[@var{label}, @var{score}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the Support Vector Machine
    ## classification model from a CompactClassificationSVM object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## predictor data in the CompactClassificationSVM model, @var{obj}.  For
    ## one-class SVM model, +1 or -1 is returned.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactClassificationSVM} class object.
    ## @item
    ## @var{XC} must be an @math{MxP} numeric matrix with the same number of
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
    ## @seealso{CompactClassificationSVM, ClassificationSVM}
    ## @end deftypefn

    function [labels, scores] = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("CompactClassificationSVM.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("CompactClassificationSVM.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat (["CompactClassificationSVM.predict: XC must have"], ...
                       [" the same number of predictors as the trained"], ...
                       [" SVM model."]));
      endif

      ## Standardize (if necessary)
      if (this.Standardize)
        XC = (XC - this.Mu) ./ this.Sigma;
      endif

      ## Predict labels and scores from new data
      [out, ~, scores] = svmpredict (ones (rows (XC), 1), XC, this.Model, '-q');

      ## Expand scores for two classes
      if (numel (this.ClassNames) == 2)
        scores = [scores, -scores];
      endif

      ## Translate labels to classnames
      if (iscellstr (this.ClassNames))
        labels = cell (rows (XC), 1);
        labels(out==1) = this.ClassNames{1};
        labels(out!=1) = this.ClassNames{2};
      elseif (islogical (this.ClassNames))
        labels = false (rows (XC), 1);
      elseif (isnumeric (this.ClassNames))
        labels = zeros (rows (XC), 1);
      elseif (ischar (this.ClassNames))
        labels = char (zeros (rows (XC), size (this.ClassNames, 2)));
      endif
      if (! iscellstr (this.ClassNames))
        labels(out==1) = this.ClassNames(1);
        labels(out!=1) = this.ClassNames(2);
      endif

      if (nargout > 1)
        ## Apply ScoreTransform to return probability estimates
        if (! strcmp (this.ScoreTransform, "none"))
          f = this.ScoreTransform;
          if (! strcmp (class (f), "function_handle"))
            error (strcat (["CompactClassificationSVM.predict: 'Score"], ...
                           ["Transform' must be a 'function_handle' object."]));
          endif
          scores = f (scores);
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationSVM} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margins for Support Vector Machine classifier.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns
    ## the classification margins for @var{obj} with data @var{X} and
    ## classification @var{Y}.  @var{m} is a numeric vector of length size (X,1).
    ##
    ## @itemize
    ## @item
    ## @var{obj} is a @var{CompactClassificationSVM} object.
    ## @item
    ## @var{X} must be a @math{NxP} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @var{Y} is @math{Nx1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}.  @var{Y} must have same
    ## numbers of Rows as @var{X}.
    ## @end itemize
    ##
    ## The classification margin for each observation is the difference between
    ## the classification score for the true class and the maximal
    ## classification score for the false classes.
    ##
    ## @seealso{CompactClassificationSVM}
    ## @end deftypefn

    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("CompactClassificationSVM.margin: too few input arguments.");
      endif

      ## Check for valid X
      if (isempty (X))
        error ("CompactClassificationSVM.margin: X is empty.");
      elseif (this.NumPredictors != columns (X))
        error (strcat (["CompactClassificationSVM.margin: X must"], ...
                       [" have the same number of predictors as"], ...
                       [" the trained SVM model."]));
      endif

      ## Check for valid Y
      if (isempty (Y))
        error ("CompactClassificationSVM.margin: Y is empty.");
      elseif (rows (X) != rows (Y))
        error (strcat (["CompactClassificationSVM.margin: Y must have"], ...
                       [" the same number of rows as X."]));
      endif

      [~, ~, dec_values_L] = svmpredict (Y, X, this.Model, '-q');
      m = 2 * Y .* dec_values_L;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationSVM} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationSVM} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute loss for a trained CompactClassificationSVM object.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} computes the loss,
    ## @var{L}, using the default loss function @qcode{'classiferror'}.
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{CompactClassificationSVM} object.
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
    ## @seealso{CompactClassificationSVM}
    ## @end deftypefn

    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("CompactClassificationSVM.loss: too few input arguments.");
      endif

      if (mod (nargin, 2) == 0)
        error (strcat (["CompactClassificationSVM.loss: Name-Value"], ...
                       [" arguments must be in pairs."]));
      endif

      ## Check for valid X
      if (isempty (X))
        error ("CompactClassificationSVM.loss: X is empty.");
      elseif (this.NumPredictors != columns (X))
        error (strcat (["CompactClassificationSVM.loss: X must"], ...
                       [" have the same number of predictors as"], ...
                       [" the trained SVM model."]));
      endif

      ## Check for valid Y
      if (isempty (Y))
        error ("CompactClassificationSVM.loss: Y is empty.");
      elseif (rows (X)!= rows (Y))
        error (strcat (["CompactClassificationSVM.loss: Y must have"], ...
                       [" the same number of rows as X."]));
      endif

      ## Set default values before parsing optional parameters
      LossFun = 'classiferror';
      Weights = ones (size (X, 1), 1);

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "lossfun"
            LossFun = varargin{2};
            if (! (ischar (LossFun)))
              error (strcat (["CompactClassificationSVM.loss: 'LossFun'"], ...
                             [" must be a character vector."]));
            endif
            LossFun = tolower (LossFun);
            if (! any (strcmpi (LossFun, {"binodeviance", "classiferror", ...
                                          "exponential", "hinge", "logit", ...
                                          "quadratic"})))
              error (strcat (["CompactClassificationSVM.loss:"], ...
                             [" unsupported Loss function."]));
            endif

          case "weights"
            Weights = varargin{2};
            ## Validate if weights is a numeric vector
            if(! (isnumeric (Weights) && isvector (Weights)))
              error (strcat (["CompactClassificationSVM.loss: 'Weights'"], ...
                             [" must be a numeric vector."]));
            endif

            ## Check if the size of weights matches the number of rows in X
            if (numel (Weights) != size (X, 1))
              error (strcat (["CompactClassificationSVM.loss: size of"], ...
                             [" 'Weights' must be equal to the number"], ...
                             [" of rows in X."]));
            endif

          otherwise
            error (strcat (["CompactClassificationSVM.loss: invalid"], ...
                           [" parameter name in optional pair arguments."]));
          endswitch
        varargin (1:2) = [];
      endwhile

      ## Compute the classification score
      [~, ~, dec_values_L] = svmpredict (Y, X, this.Model, '-q');

        ## Compute the margin
        margin = Y .* dec_values_L;

        ## Compute the loss based on the specified loss function
        switch (LossFun)
          case "classiferror"
            L = mean ((margin <= 0) .* Weights);

          case "hinge"
            L = mean (max (0, 1 - margin) .* Weights);

          case "logit"
            L = mean (log (1 + exp (-margin)) .* Weights);

          case "exponential"
            L = mean (exp (-margin) .* Weights);

          case "quadratic"
            L = mean (((1 - margin) .^2) .* Weights);

          case "binodeviance"
            L = mean (log (1 + exp (-2 * margin)) .* Weights);

          otherwise
            error ("CompactClassificationSVM.loss: unsupported Loss function.");
        endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationSVM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a CompactClassificationSVM object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## CompactClassificationSVM object into an Octave binary file, the name of
    ## which is specified in @var{filename}, along with an extra variable,
    ## which defines the type classification object these variables constitute.
    ## Use @code{loadmodel} in order to load a classification object into
    ## Octave's workspace.
    ##
    ## @seealso{loadmodel, ClassificationSVM, CompactClassificationSVM}
    ## @end deftypefn

    function savemodel (this, fname)
      ## Generate variable for class name
      classdef_name = "CompactClassificationSVM";

      ## Create variables from model properties
      NumPredictors       = this.NumPredictors;
      PredictorNames      = this.PredictorNames;
      ResponseName        = this.ResponseName;
      ClassNames          = this.ClassNames;
      ScoreTransform      = this.ScoreTransform;
      Standardize         = this.Standardize;
      Sigma               = this.Sigma;
      Mu                  = this.Mu;
      ModelParameters     = this.ModelParameters;
      Model               = this.Model;
      Alpha               = this.Alpha;
      Beta                = this.Beta;
      Bias                = this.Bias;
      IsSupportVector     = this.IsSupportVector;
      SupportVectorLabels = this.SupportVectorLabels;
      SupportVectors      = this.SupportVectors;

      ## Save classdef name and all model properties as individual variables
      save ("-binary", fname, "classdef_name", "NumPredictors", ...
            "PredictorNames", "ResponseName", "ClassNames", ...
            "ScoreTransform", "Standardize", "Sigma", "Mu", ...
            "ModelParameters", "Model", "Alpha", "Beta", "Bias", ...
            "IsSupportVector", "SupportVectorLabels", "SupportVectors");
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationSVM object
      mdl = CompactClassificationSVM ();

      ## Check that fieldnames in DATA match properties in
      ## CompactClassificationSVM
      names = fieldnames (data);
      props = fieldnames (mdl);
      if (! isequal (sort (names), sort (props)))
        msg = "CompactClassificationSVM.load_model: invalid model in '%s'.";
        error (msg, filename);
      endif

      ## Copy data into object
      for i = 1:numel (props)
        mdl.(props{i}) = data.(props{i});
      endfor
    endfunction

  endmethods

endclassdef

%!demo
%! ## Create a support vectors machine classifier and its compact version
%! # and compare their size
%!
%! load fisheriris
%! X = meas;
%! Y = species;
%!
%! selected_classes = unique (Y)(randperm (3, 2));
%! selected_indices = ismember (Y, selected_classes);
%! X_selected = X(selected_indices, :);
%! Y_selected = Y(selected_indices);
%! Mdl = fitcsvm (X_selected, Y_selected, 'ClassNames', selected_classes);
%! CMdl = crossval (Mdl)

## Test input validation for constructor
%!error<CompactClassificationSVM: invalid classification object.> ...
%! CompactClassificationSVM (1)

## Test output for predict method
%!shared x, y, CMdl
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%!test
%! xc = [min(x); mean(x); max(x)];
%! Mdl = fitcsvm (x, y, 'KernelFunction', 'rbf', 'Tolerance', 1e-7);
%! CMdl = compact (Mdl);
%! assert (isempty (CMdl.Alpha), true)
%! assert (sum (CMdl.IsSupportVector), numel (CMdl.Beta))
%! [label, score] = predict (CMdl, xc);
%! assert (label, [1; 2; 2]);
%! assert (score(:,1), [0.99285; -0.080296; -0.93694], 1e-5);
%! assert (score(:,1), -score(:,2), eps)
%!test
%! Mdl = fitcsvm (x, y);
%! CMdl = compact (Mdl);
%! assert (isempty (CMdl.Beta), true)
%! assert (sum (CMdl.IsSupportVector), numel (CMdl.Alpha))
%! assert (numel (CMdl.Alpha), 24)
%! assert (CMdl.Bias, -14.415, 1e-3)
%! xc = [min(x); mean(x); max(x)];
%! label = predict (CMdl, xc);
%! assert (label, [1; 2; 2]);

## Test input validation for predict method
%!error<CompactClassificationSVM.predict: too few input arguments.> ...
%! predict (CMdl)
%!error<CompactClassificationSVM.predict: XC is empty.> ...
%! predict (CMdl, [])
%!error<CompactClassificationSVM.predict: XC must have the same number of predictors as the trained SVM model.> ...
%! predict (CMdl, 1)
%!test
%! CMdl.ScoreTransform = "a";
%!error<CompactClassificationSVM.predict: 'ScoreTransform' must be a 'function_handle' object.> ...
%! [labels, scores] = predict (CMdl, x);

## Test output for margin method
%!test
%! rand ("seed", 1);
%! C = cvpartition (y, 'HoldOut', 0.15);
%! Mdl = fitcsvm (x(training (C),:), y(training (C)), ...
%!                'KernelFunction', 'rbf', 'Tolerance', 1e-7);
%! CMdl = compact (Mdl);
%! testInds = test (C);
%! expected_margin = [2.0000;  0.8579;  1.6690;  3.4141;  3.4552; ...
%!                    2.6605;  3.5251; -4.0000; -6.3411; -6.4511; ...
%!                   -3.0532; -7.5054; -1.6700; -5.6227; -7.3640];
%! computed_margin = margin (CMdl, x(testInds,:), y(testInds,:));
%! assert (computed_margin, expected_margin, 1e-4);

## Test input validation for margin method
%!error<CompactClassificationSVM.margin: too few input arguments.> ...
%! margin (CMdl)
%!error<CompactClassificationSVM.margin: too few input arguments.> ...
%! margin (CMdl, zeros (2))
%!error<CompactClassificationSVM.margin: X is empty.> ...
%! margin (CMdl, [], 1)
%!error<CompactClassificationSVM.margin: X must have the same number of predictors as the trained SVM model.> ...
%! margin (CMdl, 1, 1)
%!error<CompactClassificationSVM.margin: Y is empty.> ...
%! margin (CMdl, [1, 2], [])
%!error<CompactClassificationSVM.margin: Y must have the same number of rows as X.> ...
%! margin (CMdl, [1, 2], [1; 2])

## Test output for loss method
%!test
%! rand ("seed", 1);
%! C = cvpartition (y, 'HoldOut', 0.15);
%! Mdl = fitcsvm (x(training (C),:), y(training (C)), ...
%!                'KernelFunction', 'rbf', 'Tolerance', 1e-7);
%! CMdl = compact (Mdl);
%! testInds = test (C);
%! L1 = loss (CMdl, x(testInds,:), y(testInds,:), 'LossFun', 'binodeviance');
%! L2 = loss (CMdl, x(testInds,:), y(testInds,:), 'LossFun', 'classiferror');
%! L3 = loss (CMdl, x(testInds,:), y(testInds,:), 'LossFun', 'exponential');
%! L4 = loss (CMdl, x(testInds,:), y(testInds,:), 'LossFun', 'hinge');
%! L5 = loss (CMdl, x(testInds,:), y(testInds,:), 'LossFun', 'logit');
%! L6 = loss (CMdl, x(testInds,:), y(testInds,:), 'LossFun', 'quadratic');
%! assert (L1, 2.8711, 1e-4);
%! assert (L2, 0.5333, 1e-4);
%! assert (L3, 10.9685, 1e-4);
%! assert (L4, 1.9827, 1e-4);
%! assert (L5, 1.5849, 1e-4);
%! assert (L6, 7.6739, 1e-4);

## Test input validation for loss method
%!error<CompactClassificationSVM.loss: too few input arguments.> ...
%! loss (CMdl)
%!error<CompactClassificationSVM.loss: too few input arguments.> ...
%! loss (CMdl, zeros (2))
%!error<CompactClassificationSVM.loss: Name-Value arguments must be in pairs.> ...
%! loss (CMdl, [1, 2], 1, "LossFun")
%!error<CompactClassificationSVM.loss: X is empty.> ...
%! loss (CMdl, [], zeros (2))
%!error<CompactClassificationSVM.loss: X must have the same number of predictors as the trained SVM model.> ...
%! loss (CMdl, 1, zeros (2))
%!error<CompactClassificationSVM.loss: Y is empty.> ...
%! loss (CMdl, [1, 2], [])
%!error<CompactClassificationSVM.loss: Y must have the same number of rows as X.> ...
%! loss (CMdl, [1, 2], [1; 2])
%!error<CompactClassificationSVM.loss: 'LossFun' must be a character vector.> ...
%! loss (CMdl, [1, 2], 1, "LossFun", 1)
%!error<CompactClassificationSVM.loss: unsupported Loss function.> ...
%! loss (CMdl, [1, 2], 1, "LossFun", "some")
%!error<CompactClassificationSVM.loss: 'Weights' must be a numeric vector.> ...
%! loss (CMdl, [1, 2], 1, "Weights", ['a', 'b'])
%!error<CompactClassificationSVM.loss: 'Weights' must be a numeric vector.> ...
%! loss (CMdl, [1, 2], 1, "Weights", 'a')
%!error<CompactClassificationSVM.loss: size of 'Weights' must be equal to the number of rows in X.> ...
%! loss (CMdl, [1, 2], 1, "Weights", [1, 2])
%!error<CompactClassificationSVM.loss: invalid parameter name in optional pair arguments.> ...
%! loss (CMdl, [1, 2], 1, "some", "some")
