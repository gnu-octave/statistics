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

classdef CompactClassificationSVM
## -*- texinfo -*-
## @deftypefn {statistics} CompactClassificationSVM
##
## A @qcode{CompactClassificationSVM} object is a compact version of a support
## vectors machine model, @qcode{CompactClassificationSVM}.
##
## The @qcode{CompactClassificationSVM} does not include the training data
## resulting to a smaller classifier size, which can be used for making
## predictions from new data, but not for tasks such as cross validation.  It
## can only be created from a @qcode{ClassificationSVM} model by using the
## @code{compact} object method.
##
## The available methods for a @qcode{CompactClassificationSVM} object
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
## @seealso{fitcsvm, ClassificationSVM}
## @end deftypefn

  properties (Access = public)

    NumPredictors       = [];    # Number of predictors
    PredictorNames      = [];    # Predictor variables names
    ResponseName        = [];    # Response variable name
    ClassNames          = [];    # Names of classes in Y
    Prior               = [];    # Prior probability for each class
    Cost                = [];    # Cost of misclassification

    ScoreTransform      = [];    # Transformation for classification scores

    Standardize         = [];    # Flag to standardize predictors
    Sigma               = [];    # Predictor standard deviations
    Mu                  = [];    # Predictor means

    ModelParameters     = [];    # SVM parameters
    Model               = [];    # Stores the 'libsvm' trained model

    Alpha               = [];    # Trained classifier coefficients
    Beta                = [];    # Linear predictor coefficients
    Bias                = [];    # Bias term

    IsSupportVector     = [];    # Indices of Support vectors
    SupportVectorLabels = [];    # Support vector class labels
    SupportVectors      = [];    # Support vectors

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
      this.Prior                 = Mdl.Prior;
      this.Cost                  = Mdl.Cost;

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
    ## @deftypefn  {CompactClassificationSVM} {@var{labels} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationSVM} {[@var{labels}, @var{scores}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the Support Vector Machine
    ## classification object.
    ##
    ## @code{@var{labels} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## trained Support Vector Machine classification compact model, @var{obj}.
    ## For one-class SVM model, +1 or -1 is returned.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactClassificationSVM} class object.
    ## @item
    ## @var{XC} must be an @math{MxP} numeric matrix with the same number of
    ## predictors @math{P} as the corresponding predictors of the SVM model in
    ## @var{obj}.
    ## @end itemize
    ##
    ## @code{[@var{labels}, @var{scores}] = predict (@var{obj}, @var{XC}} also
<<<<<<< Updated upstream
    ## returns @var{scores}, which contains the desicion values for each each
=======
    ## returns @var{scores}, which contains the decision values for each
>>>>>>> Stashed changes
    ## prediction.   Alternatively, @var{scores} can contain the posterior
    ## probabilities if the ScoreTransform has been previously set using the
    ## @code{fitPosterior} method.
    ##
    ## @seealso{fitcsvm, ClassificationSVM.fitPosterior}
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
    ## Determine the classification margins for a Support Vector Machine
    ## classification object.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns the
    ## classification margins for the trained support vector machine (SVM)
    ## classifier @var{obj} using the sample data in @var{X} and the class
    ## labels in @var{Y}. It supports only binary classifier models. The
    ## classification margin is commonly defined as @var{m} = @var{y}f(@var{x}),
    ## where @var{f(x)} is the classification score and @var{y} is the true
    ## class label corresponding to @var{x}. A greater margin indicates a better
    ## model.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a binary class @qcode{CompactClassificationSVM} object.
    ## @item
    ## @var{X} must be an @math{MxP} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the SVM model in
    ## @var{obj}.
    ## @item
    ## @var{Y} must be @math{Mx1} numeric vector containing the class labels
    ## corresponding to the predictor data in @var{X}. @var{Y} must have same
    ## number of rows as @var{X}.
    ## @end itemize
    ##
    ## @seealso{fitcsvm, CompactClassificationSVM}
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
    ## Determine the classification error for a Support Vector Machine
    ## classifier.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## predictive accuracy of support vector machine (SVM) classification models.
    ## Comparing the same type of loss across multiple models allows you to
    ## identify which model is more accurate, with a lower loss indicating
    ## superior predictive performance. It supports only binary classifier
    ## models.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a binary class @qcode{CompactClassificationSVM} object.
    ## @item
    ## @var{X} must be an @math{MxP} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the SVM model in
    ## @var{obj}.
    ## @item
    ## @var{Y} must be @math{Mx1} numeric vector containing the class labels
    ## corresponding to the predictor data in @var{X}. @var{Y} must have same
    ## number of rows as @var{X}.
    ## @end itemize
    ##
    ## @code{@var{L} = loss (@dots{}, @var{Name}, @var{Value})} returns the
    ## aforementioned results with additional properties specified by
    ## @qcode{Name-Value} pair arguments listed below.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"LossFun"} @tab @tab Loss function, specified as a built-in
    ## loss function name. It accepts the following options: (Default is
    ## 'classiferror')
    ##
    ## @itemize
    ##
    ## @item 'binodeviance': Binomial deviance:
    ## The binomial deviance loss function is used to evaluate the performance
    ## of a binary classifier. It is calculated as:
    ## @math{L = \sum_{j=1}^{n} w_j \log \{1 + \exp [-2m_j]\}}
    ##
    ## @item 'classiferror': Misclassification rate in decimal
    ## The classification error measures the fraction of misclassified instances
    ## out of the total instances. It is calculated as:
    ## @math{L = \frac{1}{n} \sum_{j=1}^{n} \mathbb{I}(m_j \leq 0)}
    ##
    ## @item 'exponential': Exponential loss:
    ## The exponential loss function is used to penalize misclassified instances
    ## exponentially. It is calculated as:
    ## @math{L = \sum_{j=1}^{n} w_j \exp [-m_j]}
    ##
    ## @item 'hinge': Hinge loss:
    ## The hinge loss function is often used for maximum-margin classification,
    ## particularly for support vector machines. It is calculated as:
    ## @math{L = \sum_{j=1}^{n} w_j \max (0, 1 - m_j)}
    ##
    ## @item 'logit': Logistic loss:
    ## The logistic loss function, also known as log loss, measures the
    ## performance of a classification model where the prediction is a
    ## probability value. It is calculated as:
    ## @math{L = \sum_{j=1}^{n} w_j \log \{1 + \exp [-m_j]\}}
    ##
    ## @item 'quadratic': Quadratic loss:
    ## The quadratic loss function penalizes the square of the margin.
    ## It is calculated as:
    ## @math{L = \sum_{j=1}^{n} w_j (1 - m_j)^2}
    ##
    ## @end itemize
    ##
    ## @item @qcode{"Weights"} @tab @tab Specified as a numeric vector which
    ## weighs each observation (row) in X. The size of Weights must be equal
    ## to the number of rows in X. The default value is: ones(size(X,1),1)
    ##
    ## @end multitable
    ##
    ## @seealso{fitcsvm, ClassificationSVM}
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
    ## CompactClassificationSVM object into an Octave binary file, the name
    ## of which is specified in @var{filename}, along with an extra variable,
    ## which defines the type classification object these variables constitute.
    ## Use @code{loadmodel} in order to load a classification object into
    ## Octave's workspace.
    ##
    ## @seealso{loadmodel, fitcsvm, ClassificationSVM, CompactClassificationSVM}
    ## @end deftypefn

    function savemodel (this, fname)
      ## Generate variable for class name
      classdef_name = "CompactClassificationSVM";

      ## Create variables from model properties
      NumPredictors       = this.NumPredictors;
      PredictorNames      = this.PredictorNames;
      ResponseName        = this.ResponseName;
      ClassNames          = this.ClassNames;
      Prior               = this.Prior;
      Cost                = this.Cost;
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
            "PredictorNames", "ResponseName", "ClassNames", "Prior", ...
            "Cost", "ScoreTransform", "Standardize", "Sigma", "Mu", ...
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
%! Mdl = fitcsvm (X, Y, 'ClassNames', unique (species))
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
