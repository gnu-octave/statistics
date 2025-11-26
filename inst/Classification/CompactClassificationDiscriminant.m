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

classdef CompactClassificationDiscriminant
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactClassificationDiscriminant
  ##
  ## Compact discriminant analysis classification
  ##
  ## The @code{CompactClassificationDiscriminant} class implements a compact
  ## version of a linear discriminant analysis classifier object, which can
  ## predict responses for new data using the @code{predict} method but does not
  ## store the training data.
  ##
  ## A @code{CompactClassificationDiscriminant} object is a compact version of a
  ## discriminant analysis model, @code{ClassificationDiscriminant}.  It does
  ## not include the training data resulting in a smaller classifier size, which
  ## can be used for making predictions from new data, but not for tasks such as
  ## cross validation.  It can only be created from a
  ## @code{ClassificationDiscriminant} model by using the @code{compact} object
  ## method.
  ##
  ## Create a @code{CompactClassificationDiscriminant} object by using the
  ## @code{compact} method of a @code{ClassificationDiscriminant} object.
  ##
  ## @seealso{fitcdiscr, ClassificationDiscriminant}
  ## @end deftp

  properties (Access = public)
    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer value specifying the number of predictors in the
    ## training dataset used for training the CompactClassificationDiscriminant
    ## model.  This property is read-only.
    ##
    ## @end deftp
    NumPredictors   = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} PredictorNames
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
    ## @deftp {CompactClassificationDiscriminant} {property:} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName    = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} ClassNames
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
    ## @deftp {CompactClassificationDiscriminant} {property:} Prior
    ##
    ## Prior probability for each class
    ##
    ## A numeric vector specifying the prior probabilities for each class.  The
    ## order of the elements in @qcode{Prior} corresponds to the order of the
    ## classes in @qcode{ClassNames}.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Prior           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} Cost
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
    ## This property is read-only.
    ##
    ## @end deftp
    Cost            = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} ScoreTransform
    ##
    ## Transformation function for classification scores
    ##
    ## Specified as a function handle for transforming the classification
    ## scores.  This property is read-only.
    ##
    ## When specified as a character vector, it can be any of the following
    ## built-in functions.  Nevertheless, the @qcode{ScoreTransform} property
    ## always stores their function handle equivalent.
    ##
    ## @multitable @columnfractions 0.2 0.05 0.75
    ## @headitem @var{Value} @tab @tab @var{Description}
    ## @item @qcode{"doublelogit"} @tab @tab @math{1 ./ (1 + exp (-2 * x))}
    ## @item @qcode{"invlogit"} @tab @tab @math{log (x ./ (1 - x))}
    ## @item @qcode{"ismax"} @tab @tab Sets the score for the class with the
    ## largest score to 1, and for all other classes to 0
    ## @item @qcode{"logit"} @tab @tab @math{1 ./ (1 + exp (-x))}
    ## @item @qcode{"none"} @tab @tab @math{x} (no transformation)
    ## @item @qcode{"identity"} @tab @tab @math{x} (no transformation)
    ## @item @qcode{"sign"} @tab @tab @math{-1 for x < 0, 0 for x = 0, 1 for x > 0}
    ## @item @qcode{"symmetric"} @tab @tab @math{2 * x - 1}
    ## @item @qcode{"symmetricismax"} @tab @tab Sets the score for the class
    ## with the largest score to 1, and for all other classes to -1
    ## @item @qcode{"symmetriclogit"} @tab @tab @math{2 ./ (1 + exp (-x)) - 1}
    ## @end multitable
    ##
    ## @end deftp
    ScoreTransform  = @(x) x;

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} Sigma
    ##
    ## Within-class covariance
    ##
    ## A numeric array specifying the within-class covariance. For linear
    ## discriminant type (currently supported) this is a @math{PxP} matrix,
    ## where @math{P} is the number of predictors.  This property is read-only.
    ##
    ## @end deftp
    Sigma           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} Mu
    ##
    ## Class means
    ##
    ## A @math{KxP} numeric matrix specifying the mean of the multivariate
    ## normal distribution of each corresponding class, where @math{K} is the
    ## number of classes and @math{P} is the number of predictors.  This property
    ## is read-only.
    ##
    ## @end deftp
    Mu              = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} Coeffs
    ##
    ## Coefficient matrices
    ##
    ## A @math{KxK} structure containing the coefficient matrices, where
    ## @math{K} is the number of classes.  If the @qcode{'FillCoeffs'} parameter
    ## was set to @qcode{'off'} in the original @code{ClassificationDiscriminant}
    ## model, then @qcode{Coeffs} is empty @qcode{([])}.  This property is
    ## read-only.
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
    ## @deftp {CompactClassificationDiscriminant} {property:} Delta
    ##
    ## Delta threshold
    ##
    ## A nonnegative scalar specifying the threshold for linear discriminant
    ## model. Currently unimplemented and fixed to 0.  This property is
    ## read-only.
    ##
    ## @end deftp
    Delta           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} DiscrimType
    ##
    ## Discriminant type
    ##
    ## A character vector specifying the type discriminant model. Currently
    ## only linear discriminant models are supported.  This property is
    ## read-only.
    ##
    ## @end deftp
    DiscrimType     = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} Gamma
    ##
    ## Gamma regularization parameter
    ##
    ## A scalar value ranging from 0 to 1, specifying the Gamma regularization
    ## parameter.  This property is read-only.
    ##
    ## @end deftp
    Gamma           = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} MinGamma
    ##
    ## Minimum value for Gamma regularization parameter
    ##
    ## A scalar value ranging from 0 to 1, specifying the minimum value that the
    ## Gamma regularization parameter can have.  This property is read-only.
    ##
    ## @end deftp
    MinGamma        = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property:} LogDetSigma
    ##
    ## Logarithm of the determinant of the within-class covariance matrix
    ##
    ## A scalar value specifying the logarithm of the determinant of the
    ## within-class covariance matrix.  This property is read-only.
    ##
    ## @end deftp
    LogDetSigma     = [];
  endproperties

  properties (Access = private, Hidden)
    STname = 'none';
  endproperties

  methods (Hidden)

    ## Constructor
    function this = CompactClassificationDiscriminant (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), "ClassificationDiscriminant"))
        error (strcat ("CompactClassificationDiscriminant: invalid", ...
                       " classification object."));
      endif

      ## Save properties to compact model
      this.NumPredictors   = Mdl.NumPredictors;
      this.PredictorNames  = Mdl.PredictorNames;
      this.ResponseName    = Mdl.ResponseName;
      this.ClassNames      = Mdl.ClassNames;

      this.Cost            = Mdl.Cost;
      this.Prior           = Mdl.Prior;
      this.ScoreTransform  = Mdl.ScoreTransform;
      this.STname          = Mdl.STname;

      this.Sigma           = Mdl.Sigma;
      this.Mu              = Mdl.Mu;
      this.Coeffs          = Mdl.Coeffs;
      this.Delta           = Mdl.Delta;
      this.DiscrimType     = Mdl.DiscrimType;
      this.Gamma           = Mdl.Gamma;
      this.MinGamma        = Mdl.MinGamma;
      this.LogDetSigma     = Mdl.LogDetSigma;

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
      fprintf ("\n  CompactClassificationDiscriminant\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      if (iscellstr (this.ClassNames))
        str = repmat ({"'%s'"}, 1, numel (this.ClassNames));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.ClassNames{:});
      else # numeric
        str = repmat ({"%d"}, 1, numel (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, this.ClassNames);
      endif
      fprintf ("%+25s: %s\n", 'ClassNames', str);
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.STname);
      fprintf ("%+25s: '%d'\n", 'NumPredictors', this.NumPredictors);
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
                         " in a CompactClassificationDiscriminant object."));
        case '{}'
          error (strcat ("Invalid {} indexing for referencing values", ...
                         " in a CompactClassificationDiscriminant object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("CompactClassificationDiscriminant.subsref: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          try
            out = this.(s.subs);
          catch
            error (strcat ("CompactClassificationDiscriminant.subsref:", ...
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
        error (strcat ("CompactClassificationDiscriminant.subsasgn:", ...
                       " chained subscripts not allowed."));
      endif
      switch s.type
        case '()'
          error (strcat ("Invalid () indexing for assigning values", ...
                         " to a CompactClassificationDiscriminant object."));
        case '{}'
          error (strcat ("Invalid {} indexing for assigning values", ...
                         " to a CompactClassificationDiscriminant object."));
        case '.'
          if (! ischar (s.subs))
            error (strcat ("CompactClassificationDiscriminant.subsasgn: '.'", ...
                           " indexing argument must be a character vector."));
          endif
          switch (s.subs)
            case 'ScoreTransform'
              name = "CompactClassificationDiscriminant";
              [this.ScoreTransform, this.STname] = parseScoreTransform ...
                                                   (varargin{2}, name);
            otherwise
              error (strcat ("CompactClassificationDiscriminant.subsasgn:", ...
                             " unrecognized or read-only property: '%s'"), ...
                             s.subs);
          endswitch
      endswitch
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationDiscriminant} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationDiscriminant} {[@var{label}, @var{score}, @var{cost}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the discriminant
    ## analysis model from a CompactClassificationDiscriminant object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## corresponding labels from the trained @qcode{ClassificationDiscriminant},
    ## model, @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactClassificationDiscriminant} class object.
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
    ## @seealso{CompactClassificationDiscriminant, fitcdiscr}
    ## @end deftypefn
    function [label, score, cost] = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error (strcat ("CompactClassificationDiscriminant.predict:", ...
                       " too few input arguments."));
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("CompactClassificationDiscriminant.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("CompactClassificationDiscriminant.predict: XC", ...
                       " must have the same number of features as the", ...
                       " trained model."));
      endif

      ## Initialize matrices
      numObservations = rows (XC);
      numClasses = numel (this.ClassNames);
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
      label = this.ClassNames(minIdx);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationDiscriminant} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationDiscriminant} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute loss for a trained CompactClassificationDiscriminant object.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} computes the loss,
    ## @var{L}, using the default loss function @qcode{'mincost'}.
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{CompactClassificationDiscriminant} object.
    ## @item
    ## @code{X} must be a @math{NxP} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of rows as @var{X}.
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
    ## @seealso{CompactClassificationDiscriminant}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error (strcat ("CompactClassificationDiscriminant.loss:", ...
                       " too few input arguments."));
      elseif (mod (nargin - 3, 2) != 0)
        error (strcat ("CompactClassificationDiscriminant.loss:", ...
                       " name-value arguments must be in pairs."));
      elseif (nargin > 7)
        error (strcat ("CompactClassificationDiscriminant.loss:", ...
                       " too many input arguments."));
      endif

      ## Default values
      LossFun = 'mincost';
      Weights = [];

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error (strcat ("CompactClassificationDiscriminant.loss:", ...
                       " Y must be of a valid type."));
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat ("CompactClassificationDiscriminant.loss: Y must", ...
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
                error (strcat ("CompactClassificationDiscriminant.loss:", ...
                               " custom loss function must accept", ...
                               " exactly four input arguments."));
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
                  error (strcat ("CompactClassificationDiscriminant.loss:", ...
                                 " custom loss function must return", ...
                                 " a scalar value."));
                endif
              catch
                error (strcat ("CompactClassificationDiscriminant.loss:", ...
                               " custom loss function is not valid or", ...
                               " does not produce correct output."));
              end_try_catch
              LossFun = Value;
            elseif (ischar (Value) && any (strcmpi (Value, lf_opt)))
              LossFun = Value;
            else
              error (strcat ("CompactClassificationDiscriminant.loss:", ...
                             " invalid loss function."));
            endif

          case 'weights'
            if (isnumeric (Value) && isvector (Value))
              if (numel (Value) != size (X ,1))
                error (strcat ("CompactClassificationDiscriminant.loss:", ...
                               " number of 'Weights' must be equal to", ...
                               " the number of rows in X."));
              elseif (numel (Value) == size (X, 1))
                Weights = Value;
              endif
            else
              error (strcat ("CompactClassificationDiscriminant.loss:", ...
                             " invalid 'Weights'."));
            endif

          otherwise
            error (strcat ("CompactClassificationDiscriminant.loss:", ...
                           " invalid parameter name in optional pair", ...
                           " arguments."));
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

      ## Convert Y to a cell array of strings
      if (ischar (Y))
        Y = cellstr (Y);
      elseif (isnumeric (Y))
        Y = cellstr (num2str (Y));
      elseif (islogical (Y))
        Y = cellstr (num2str (double (Y)));
      elseif (iscell (Y))
        Y = cellfun (@num2str, Y, 'UniformOutput', false);
      else
        error (strcat ("CompactClassificationDiscriminant.loss: Y must be", ...
                       " a numeric, logical, char, string, or cell array."));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), this.ClassNames)))
        error (strcat ("CompactClassificationDiscriminant.loss: Y must", ...
                       " contain only the classes in ClassNames."));
      endif

      ## Set default weights if not specified
      if (isempty (Weights))
        Weights = ones (size (X, 1), 1);
      endif

      ## Normalize Weights
      unique_classes = this.ClassNames;
      class_prior_probs = this.Prior;
      norm_weights = zeros (size (Weights));
      for i = 1:numel (unique_classes)
        class_idx = ismember (Y, unique_classes{i});
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
      K = numel (this.ClassNames);
      C = false (n, K);
      for i = 1:n
        class_idx = find (ismember (this.ClassNames, Y{i}));
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
            cj = Cost(find (ismember (this.ClassNames, Y(i))), min_cost_class);
            L = L + Weights(i) * cj;
          endfor
        case 'classifcost'
          Cost = this.Cost;
          L = 0;
          for i = 1:n
            y_idx = find (ismember (this.ClassNames, Y(i)));
            y_hat_idx = find (ismember (this.ClassNames, label(i)));
            L = L + Weights(i) * Cost(y_idx, y_hat_idx);
          endfor
        otherwise
          error (strcat ("CompactClassificationDiscriminant.loss:", ...
                         " invalid loss function."));
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationDiscriminant} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margins for discriminant analysis classifier.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns
    ## the classification margins for @var{obj} with data @var{X} and
    ## classification @var{Y}. @var{m} is a numeric vector of length size (X,1).
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{CompactClassificationDiscriminant} object.
    ## @item
    ## @code{X} must be a @math{NxP} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of rows as @var{X}.
    ## @end itemize
    ##
    ## The classification margin for each observation is the difference between
    ## the classification score for the true class and the maximal
    ## classification score for the false classes.
    ##
    ## @seealso{fitcdiscr, CompactClassificationDiscriminant}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error (strcat ("CompactClassificationDiscriminant.margin:", ...
                       " too few input arguments."));
      endif

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error (strcat ("CompactClassificationDiscriminant.margin:", ...
                       " Y must be of a valid type."));
      endif

      ## Validate X
      valid_types = {'single', 'double'};
      if (! (any (strcmp (class (X), valid_types))))
        error (strcat ("CompactClassificationDiscriminant.margin:", ...
                       " X must be of a valid type."));
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat ("CompactClassificationDiscriminant.margin: Y must", ...
                       " have the same number of rows as X."));
      endif

      ## Convert Y to a cell array of strings
      if (ischar (Y))
        Y = cellstr (Y);
      elseif (isnumeric (Y))
        Y = cellstr (num2str (Y));
      elseif (islogical (Y))
        Y = cellstr (num2str (double (Y)));
      elseif (iscell (Y))
        Y = cellfun (@num2str, Y, 'UniformOutput', false);
      else
        error (strcat ("CompactClassificationDiscriminant.margin: Y must", ...
                       " be a numeric, logical, char, string, or cell array."));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), this.ClassNames)))
        error (strcat ("CompactClassificationDiscriminant.margin: Y must", ...
                       " contain only the classes in ClassNames."));
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
        true_class_idx = find (ismember (this.ClassNames, Y{i}));

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
    ## @deftypefn  {CompactClassificationDiscriminant} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a CompactClassificationDiscriminant object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## CompactClassificationDiscriminant object into an Octave binary file, the
    ## name of which is specified in @var{filename}, along with an extra
    ## variable, which defines the type classification object these variables
    ## constitute. Use @code{loadmodel} in order to load a classification object
    ## into Octave's workspace.
    ##
    ## @seealso{loadmodel, fitcdiscr, ClassificationDiscriminant}
    ## @end deftypefn
    function savemodel (this, fname)
      ## Generate variable for class name
      classdef_name = "CompactClassificationDiscriminant";

      ## Create variables from model properties
      NumPredictors   = this.NumPredictors;
      PredictorNames  = this.PredictorNames;
      ResponseName    = this.ResponseName;
      ClassNames      = this.ClassNames;
      Prior           = this.Prior;
      Cost            = this.Cost;
      ScoreTransform  = this.ScoreTransform;
      STname          = this.STname;
      Sigma           = this.Sigma;
      Mu              = this.Mu;
      Coeffs          = this.Coeffs;
      Delta           = this.Delta;
      DiscrimType     = this.DiscrimType;
      Gamma           = this.Gamma;
      MinGamma        = this.MinGamma;
      LogDetSigma     = this.LogDetSigma;

      ## Save classdef name and all model properties as individual variables
      save ("-binary", fname, "classdef_name", "NumPredictors", ...
            "PredictorNames", "ResponseName", "ClassNames", "Prior", ...
            "Cost", "ScoreTransform", "STname", "Sigma", "Mu", "Coeffs", ...
            "Delta", "DiscrimType", "Gamma", "MinGamma", "LogDetSigma");
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a CompactClassificationDiscriminant object
      mdl = CompactClassificationDiscriminant ();

      ## Check that fieldnames in DATA match properties in
      ## CompactClassificationDiscriminant
      names = fieldnames (data);
      props = fieldnames (mdl);
      if (! isequal (sort (names), sort (props)))
        msg = strcat ("CompactClassificationDiscriminant.load_model:", ...
                      " invalid model in '%s'.");
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
%! ## Create a discriminant analysis classifier and its compact version
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
%! load fisheriris
%! x = meas;
%! y = species;
%! PredictorNames = {'Sepal Length', 'Sepal Width', 'Petal Length', 'Petal Width'};
%! Mdl = fitcdiscr (x, y, "PredictorNames", PredictorNames);
%! CMdl = compact (Mdl);
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
%! assert (class (CMdl), "CompactClassificationDiscriminant");
%! assert ({CMdl.DiscrimType, CMdl.ResponseName}, {"linear", "Y"})
%! assert ({CMdl.Gamma, CMdl.MinGamma}, {0, 0}, 1e-15)
%! assert (CMdl.ClassNames, unique (species))
%! assert (CMdl.Sigma, sigma, 1e-6)
%! assert (CMdl.Mu, mu, 1e-14)
%! assert (CMdl.LogDetSigma, -9.9585, 1e-4)
%! assert (CMdl.PredictorNames, PredictorNames)
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! Mdl = fitcdiscr (x, y, "Gamma", 0.5);
%! CMdl = compact (Mdl);
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
%! assert (class (CMdl), "CompactClassificationDiscriminant");
%! assert ({CMdl.DiscrimType, CMdl.ResponseName}, {"linear", "Y"})
%! assert ({CMdl.Gamma, CMdl.MinGamma}, {0.5, 0})
%! assert (CMdl.ClassNames, unique (species))
%! assert (CMdl.Sigma, sigma, 1e-6)
%! assert (CMdl.Mu, mu, 1e-14)
%! assert (CMdl.LogDetSigma, -8.6884, 1e-4)

## Test input validation for constructor
%!error<CompactClassificationDiscriminant: invalid classification object.> ...
%! CompactClassificationDiscriminant (1)

## Test predict method
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! Mdl = fitcdiscr (meas, species, "Gamma", 0.5);
%! CMdl = compact (Mdl);
%! [label, score, cost] = predict (CMdl, [2, 2, 2, 2]);
%! assert (label, {'versicolor'})
%! assert (score, [0, 0.9999, 0.0001], 1e-4)
%! assert (cost, [1, 0.0001, 0.9999], 1e-4)
%! [label, score, cost] = predict (CMdl, [2.5, 2.5, 2.5, 2.5]);
%! assert (label, {'versicolor'})
%! assert (score, [0, 0.6368, 0.3632], 1e-4)
%! assert (cost, [1, 0.3632, 0.6368], 1e-4)
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! xc = [min(x); mean(x); max(x)];
%! Mdl = fitcdiscr (x, y);
%! CMdl = compact (Mdl);
%! [label, score, cost] = predict (CMdl, xc);
%! l = {'setosa'; 'versicolor'; 'virginica'};
%! s = [1, 0, 0; 0, 1, 0; 0, 0, 1];
%! c = [0, 1, 1; 1, 0, 1; 1, 1, 0];
%! assert (label, l)
%! assert (score, s, 1e-4)
%! assert (cost, c, 1e-4)

%!shared MODEL
%! X = rand (10,2);
%! Y = [ones(5,1);2*ones(5,1)];
%! MODEL = compact (ClassificationDiscriminant (X, Y));

## Test input validation for predict method
%!error<CompactClassificationDiscriminant.predict: too few input arguments.> ...
%! predict (MODEL)
%!error<CompactClassificationDiscriminant.predict: XC is empty.> ...
%! predict (MODEL, [])
%!error<CompactClassificationDiscriminant.predict: XC must have the same number of features as the trained model.> ...
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
%!error<CompactClassificationDiscriminant.loss: too few input arguments.> ...
%! loss (MODEL)
%!error<CompactClassificationDiscriminant.loss: too few input arguments.> ...
%! loss (MODEL, ones (4,2))
%!error<CompactClassificationDiscriminant.loss: name-value arguments must be in pairs.> ...
%! loss (MODEL, ones (4,2), ones (4,1), 'LossFun')
%!error<CompactClassificationDiscriminant.loss: Y must have the same number of rows as X.> ...
%! loss (MODEL, ones (4,2), ones (3,1))
%!error<CompactClassificationDiscriminant.loss: invalid loss function.> ...
%! loss (MODEL, ones (4,2), ones (4,1), 'LossFun', 'a')
%!error<CompactClassificationDiscriminant.loss: invalid 'Weights'.> ...
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
%!error<CompactClassificationDiscriminant.margin: too few input arguments.> ...
%! margin (MODEL)
%!error<CompactClassificationDiscriminant.margin: too few input arguments.> ...
%! margin (MODEL, ones (4,2))
%!error<CompactClassificationDiscriminant.margin: Y must have the same number of rows as X.> ...
%! margin (MODEL, ones (4,2), ones (3,1))
