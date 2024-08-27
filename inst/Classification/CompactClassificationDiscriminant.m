## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {statistics} CompactClassificationDiscriminant
##
## A @qcode{CompactClassificationDiscriminant} object is a compact version of a
## discriminant analysis model, @qcode{ClassificationDiscriminant}.
##
## The @qcode{CompactClassificationDiscriminant} does not include the training
## data resulting to a smaller classifier size, which can be used for making
## predictions from new data, but not for tasks such as cross validation.  It
## can only be created from a @qcode{ClassificationDiscriminant} model by using
## the @code{compact} object method.
##
## The available methods for a @qcode{CompactClassificationDiscriminant} object
## are:
## @itemize
## @item
## @code{predict}
## @item
## @code{loss}
## @item
## @code{margin}
## @item
## @code{savemodel}
## @end itemize
##
## @seealso{fitcdiscr, compact, ClassificationDiscriminant}
## @end deftypefn

  properties (Access = public)

    NumPredictors   = [];     # Number of predictors
    PredictorNames  = [];     # Predictor variables names
    ResponseName    = [];     # Response variable name
    ClassNames      = [];     # Names of classes in Y
    Prior           = [];     # Prior probability for each class
    Cost            = [];     # Cost of Misclassification

    ScoreTransform  = [];     # Transformation for classification scores

    Sigma           = [];     # Within-class covariance
    Mu              = [];     # Class means
    Coeffs          = [];     # Coefficient matrices
    Delta           = [];     # Threshold for linear discriminant model
    DiscrimType     = [];     # Discriminant type
    Gamma           = [];     # Gamma regularization parameter
    MinGamma        = [];     # Minmum value of Gamma
    LogDetSigma     = [];     # Log of det of within-class covariance matrix
    XCentered       = [];     # X data with class means subtracted

  endproperties

  methods (Hidden)

    ## constructor
    function this = CompactClassificationDiscriminant (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), "ClassificationDiscriminant"))
        error (strcat (["CompactClassificationDiscriminant: invalid"], ...
                       [" classification object."]));
      endif

      ## Save properties to compact model
      this.NumPredictors   = Mdl.NumPredictors;
      this.PredictorNames  = Mdl.PredictorNames;
      this.ResponseName    = Mdl.ResponseName;
      this.ClassNames      = Mdl.ClassNames;
      this.Prior           = Mdl.Prior;
      this.Cost            = Mdl.Cost;

      this.ScoreTransform  = Mdl.ScoreTransform;

      this.Sigma           = Mdl.Sigma;
      this.Mu              = Mdl.Mu;
      this.Coeffs          = Mdl.Coeffs;
      this.Delta           = Mdl.Delta;
      this.DiscrimType     = Mdl.DiscrimType;
      this.Gamma           = Mdl.Gamma;
      this.MinGamma        = Mdl.MinGamma;
      this.LogDetSigma     = Mdl.LogDetSigma;
      this.XCentered       = Mdl.XCentered;

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
        error ("CompactClassificationDiscriminant.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("CompactClassificationDiscriminant.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat (["CompactClassificationDiscriminant.predict: XC"], ...
                       [" must have the same number of features as the"], ...
                       [" trained model."]));
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
    ## @code{obj} is a @var{CompactClassificationDiscriminant} object trained on
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
    ## @seealso{CompactClassificationDiscriminant}
    ## @end deftypefn

    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                       [" too few input arguments."]));
      elseif (mod (nargin - 3, 2) != 0)
        error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                       [" name-value arguments must be in pairs."]));
      elseif (nargin > 7)
        error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                       [" too many input arguments."]));
      endif

      ## Default values
      LossFun = 'mincost';
      Weights = [];

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                       [" Y must be of a valid type."]));
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat (["CompactClassificationDiscriminant.loss: Y must"], ...
                       [" have the same number of rows as X."]));
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
                error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                               [" custom loss function must accept"], ...
                               [" exactly four input arguments."]));
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
                  error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                                 [" custom loss function must return"], ...
                                 [" a scalar value."]));
                endif
              catch
                error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                               [" custom loss function is not valid or"], ...
                               [" does not produce correct output."]));
              end_try_catch
              LossFun = Value;
            elseif (ischar (Value) && any (strcmpi (Value, lf_opt)))
              LossFun = Value;
            else
              error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                             [" invalid loss function."]));
            endif

          case 'weights'
            if (isnumeric (Value) && isvector (Value))
              if (numel (Value) != size (X ,1))
                error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                               [" number of 'Weights' must be equal to"], ...
                               [" the number of rows in X."]));
              elseif (numel (Value) == size (X, 1))
                Weights = Value;
              endif
            else
              error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                             [" invalid 'Weights'."]));
            endif

          otherwise
            error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                           [" invalid parameter name in optional pair"], ...
                           [" arguments."]));
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
        error (strcat (["CompactClassificationDiscriminant.loss: Y must be"], ...
                       [" a numeric, logical, char, string, or cell array."]));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), this.ClassNames)))
        error (strcat (["CompactClassificationDiscriminant.loss: Y must"], ...
                       [" contain only the classes in ClassNames."]));
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
          error (strcat (["CompactClassificationDiscriminant.loss:"], ...
                         [" invalid loss function."]));
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationDiscriminant} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns
    ## the classification margins for @var{obj} with data @var{X} and
    ## classification @var{Y}. @var{m} is a numeric vector of length size (X,1).
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{CompactClassificationDiscriminant} object trained on @code{X}
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
    ## @seealso{fitcdiscr, CompactClassificationDiscriminant}
    ## @end deftypefn

    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error (strcat (["CompactClassificationDiscriminant.margin:"], ...
                       [" too few input arguments."]));
      endif

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error (strcat (["CompactClassificationDiscriminant.margin:"], ...
                       [" Y must be of a valid type."]));
      endif

      ## Validate X
      valid_types = {'single', 'double'};
      if (! (any (strcmp (class (X), valid_types))))
        error (strcat (["CompactClassificationDiscriminant.margin:"], ...
                       [" X must be of a valid type."]));
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat (["CompactClassificationDiscriminant.margin: Y must"], ...
                       [" have the same number of rows as X."]));
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
        error (strcat (["CompactClassificationDiscriminant.margin: Y must be"], ...
                       [" a numeric, logical, char, string, or cell array."]));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), this.ClassNames)))
        error (strcat (["CompactClassificationDiscriminant.margin: Y must"], ...
                       [" contain only the classes in ClassNames."]));
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
    ## @code{savemodel (@var{obj}, @var{filename})} saves a
    ## CompactClassificationDiscriminant object into a file defined by @var{filename}.
    ##
    ## @seealso{loadmodel, fitcdiscr, CompactClassificationDiscriminant}
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
      save (fname, "classdef_name", "NumPredictors", "PredictorNames", ...
            "ResponseName", "ClassNames", "Prior", "Cost", "ScoreTransform", ...
            "Sigma", "Mu", "Coeffs", "Delta", "DiscrimType", "Gamma", ...
            "MinGamma", "LogDetSigma", "XCentered");
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
        msg = strcat (["CompactClassificationDiscriminant.load_model:"], ...
                      [" invalid model in '%s'."]);
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
%! CMdl = crossval (Mdl);
%!
%! whos ('Mdl', 'CMdl')

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
%! assert (CMdl.XCentered([1:3],:), xCentered, 1e-14)
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
%! assert (CMdl.XCentered([1:3],:), xCentered, 1e-14)
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
