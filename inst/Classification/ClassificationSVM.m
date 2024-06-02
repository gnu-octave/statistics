## Copyright (C) 2024 Pallav Purbia <pallavpurbia@gmail.com>
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

classdef ClassificationSVM
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{obj} =} ClassificationSVM (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} ClassificationSVM (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{ClassificationSVM} class object containing a Support Vector
## Machine classification model.
##
## @code{@var{obj} = ClassificationSVM (@var{X}, @var{Y})} returns a
## ClassificationSVM object, with @var{X} as the predictor data and @var{Y}
## containing the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the SVM model.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain any type of
## categorical data. @var{Y} must have same numbers of Rows as @var{X}.
## @end itemize
##
## @code{@var{obj} = ClassificationSVM (@dots{}, @var{name}, @var{value})}
## returns a ClassificationSVM object with parameters specified by
## @qcode{Name-Value} pair arguments.  Type @code{help fitcsvm} for more info.
##
## A @qcode{ClassificationSVM} object, @var{obj}, stores the labelled training
## data and various parameters for the Support Vector machine classification
## model, which can be accessed in the following fields:
##
## @multitable @columnfractions 0.02 0.35 0.7
## @headitem @tab @var{Field} @tab @var{Description}
##
## @item @tab @qcode{"obj.X"} @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @tab @qcode{"obj.Y"} @tab Class labels, specified as a logical or
## numeric vector, or cell array of character vectors.  Each value in @var{Y} is
## the observed class label for the corresponding row in @var{X}.
##
## @item @tab @qcode{"obj.CrossValidationAccuracy"} @tab This property is
## computed during the cross-validation process and represents the proportion of
## correctly classified instances in the validation sets. It is calculated by
## performing a k-fold cross-validation on the training dataset. The dataset is
## divided into k subsets (folds), and the SVM model is trained on k-1 subsets
## while the remaining subset is used for validation. This process is repeated
## k times, each time using a different subset for validation. The accuracy is
## the average of the accuracies from each fold.
##
## @item @tab @qcode{"obj.ModelParameters"} @tab  This field contains the
## parameters used to train the SVM model, such as C, gamma, kernel type, etc.
## These parameters define the behavior and performance of the SVM. For example,
## 'C' controls the trade-off between achieving a low training error and a low
## testing error, 'gamma' defines the influence of a single training example,
## and 'kernel type' specifies the type of transformation applied to the input.
##
## @item @tab @qcode{"obj.NumClasses"} @tab The number of classes in the
## classification problem. For a binary classification, NumClasses is 2. In the
## case of a one-class SVM, NumClasses is also considered as 2 because the
## one-class SVM tries to separate data from one class against all other
## possible instances.
##
## @item @tab @qcode{"obj.SupportVectorCount"} @tab The total number of support
## vectors in the model. Support vectors are the data points that lie closest to
## the decision surface (or hyperplane) and are most difficult to classify. They
## are critical elements of the training dataset as they directly influence the
## position and orientation of the decision surface.
##
## @item @tab @qcode{"obj.Rho"} @tab Rho is the bias term in the decision
## function @math{sgn(w^Tx - rho)}. It represents the offset of the hyperplane
## from the origin. In other words, it is the value that helps to determine the
## decision boundary in the feature space, allowing the SVM to make
## classifications.
##
## @item @tab @qcode{"obj.ClassNames"} @tab The labels for each class in the
## classification problem. It provides the actual names or identifiers for the
## classes being predicted by the model. This field is empty for one-class SVM
## because it only involves a single class during training and testing.
##
## @item @tab @qcode{"obj.SupportVectorIndices"} @tab Indices of the support
## vectors in the training dataset. This field indicates the positions of the
## support vectors within the original training data. It helps in identifying
## which data points are the most influential in constructing the decision
## boundary.
##
## @item @tab @qcode{"obj.ProbA"} @tab Pairwise probability estimates for binary
## classification problem. This field is empty if the Probability_estimates is
## set to 0 or in one-class SVM. It is part of the pairwise coupling method used
## to estimate the probability that a data point belongs to a particular class.
##
## @item @tab @qcode{"obj.ProbB"} @tab Pairwise probability estimates for binary
## classification problem. This field is empty if the Probability_estimates is
## set to 0 or in one-class SVM. Similar to ProbA, this field is used in
## conjunction with ProbA to provide probability estimates of class memberships.
##
## @item @tab @qcode{"obj.SupportVectorPerClass"} @tab The number of support
## vectors for each class. This field provides a count of how many support
## vectors belong to each class. This field is empty for one-class SVM because
## it does not categorize support vectors by class.
##
## @item @tab @qcode{"obj.SupportVectorCoef"} @tab Coefficients for the support
## vectors in the decision functions. It contains all the @math{alpha_i * y_i},
## where alpha_i are the Lagrange multipliers and y_i are the class labels.
## These coefficients are used to scale the influence of each support vector on
## the decision boundary.
##
## @item @tab @qcode{"obj.SupportVectors"} @tab It contains all the support
## vectors. Support vectors are the critical elements of the training data that
## are used to define the position of the decision boundary in the feature
## space. They are the data points that are most informative for the
## classification task.
##
## @item @tab @qcode{"obj.Solver"} @tab This field specifies the algorithm used
## for training the SVM model, such as SMO or ISDA. It determines how the
## optimization problem is solved to find the support vectors and model
## parameters.
##
## @end multitable
##
## @seealso{fitcsvm, svmtrain, svmpredict}
## @end deftypefn

  properties (Access = public)

    X                       = [];     # Predictor data
    Y                       = [];     # Class labels

    CrossValidationAccuracy = [];     # Cross Validation Accuracy of the model
    ModelParameters         = [];     # SVM parameters.
    NumClasses              = [];     # Number of classes in Y
    ClassNames              = [];     # Names of classes in Y
    Rho                     = [];     # Bias term

    SupportVectors          = [];     # Support vectors
    SupportVectorIndices    = [];     # Indices of Support vectors
    SupportVectorCount      = [];     # Total number of Support vectors
    SupportVectorPerClass   = [];     # Number of Support vectors for each class
    SupportVectorCoef       = [];     # Coefficients of support vectors in the decision functions

    ProbA                   = [];     # Pairwise probability estimates
    ProbB                   = [];     # Pairwise probability estimates
    Solver                  = "SMO";  # Solver used

  endproperties

  properties (Access = private)

    Model                   = [];

  endproperties

  methods (Access = public)

    ## Class object constructor
    function this = ClassificationSVM (X, Y, varargin)
      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationSVM: too few input arguments.");
      endif

      ## Get training sample size and number of variables in training data
      nsample = rows (X);                    #Number of samples in X
      ndims_X = columns (X);                 #Number of dimensions in X

##      For debugging
##      disp(nsample);
##      disp(ndims_X);
##      disp(rows(Y));

      ## Check correspodence between predictors and response
      if (nsample != rows (Y))
        error ("ClassificationSVM: number of rows in X and Y must be equal.");
      endif

      ## Validate that Y is numeric
      if (!isnumeric(Y))
        error ("ClassificationSVM: Y must be a numeric array.");
      endif

      SVMtype                 = 'C_SVC';
      KernelFunction          = 'rbf';
      PolynomialOrder         = 3;
      Gamma                   = 1 / (ndims_X);
      KernelOffset            = 0;
      BoxConstraint           = 1;
      Nu                      = 0.5;
      CacheSize               = 100;
      Tolerance               = 1e-3;
      Shrinking               = 1;
      ProbabilityEstimates    = 0;
      Weight                  = [];
      KFold                   = 10;

      weight_given            = "no";

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "svmtype"
            SVMtype = varargin{2};
            if (!(ischar(SVMtype)))
              error("ClassificationSVM: SVMtype must be a string.");
            endif
            if (ischar(SVMtype))
              if (! any (strcmpi (tolower(SVMtype), {"c_svc", "nu_svc",  ...
                "one_class_svm"})))
              error ("ClassificationSVM: unsupported SVMtype.");
              endif
            endif

          case "kernelfunction"
            KernelFunction = varargin{2};
            if (!(ischar(KernelFunction)))
              error("ClassificationSVM: KernelFunction must be a string.");
            endif
            if (ischar(KernelFunction))
              if (! any (strcmpi (tolower(KernelFunction), {"linear", "rbf", ...
                "polynomial", "sigmoid"})))
              error ("ClassificationSVM: unsupported Kernel function.");
              endif
            endif

          case "polynomialorder"
            PolynomialOrder = varargin{2};
            if (! (isnumeric(PolynomialOrder) && isscalar(PolynomialOrder)
              && (PolynomialOrder > 0) && mod(PolynomialOrder, 1) == 0))
              error (strcat(["ClassificationSVM: PolynomialOrder must be a"], ...
              [" positive integer."]));
            endif

          case "gamma"
            Gamma = varargin{2};
            if ( !(isscalar(Gamma) && (Gamma > 0)))
              error ("ClassificationSVM: Gamma must be a positive scalar.");
            endif

          case "kerneloffset"
            KernelOffset = varargin{2};
            if (! (isnumeric(KernelOffset) && isscalar(KernelOffset)
              && (KernelOffset >= 0)))
              error (strcat(["ClassificationSVM: KernelOffset must be a non"], ...
              ["-negative scalar."]));
            endif

          case "boxconstraint"
            BoxConstraint = varargin{2};
            if ( !(isscalar(BoxConstraint) && (BoxConstraint > 0)))
              error ("ClassificationSVM: BoxConstraint must be a positive scalar.");
            endif

          case "nu"
            Nu = varargin{2};
            if ( !(isscalar(Nu) && (Nu > 0) && (Nu <= 1)))
              error (strcat(["ClassificationSVM: Nu must be a positive scalar"], ...
              [" in the range 0 < Nu <= 1."]));
            endif

          case "cachesize"
            CacheSize = varargin{2};
            if ( !(isscalar(CacheSize) && CacheSize > 0))
              error ("ClassificationSVM: CacheSize must be a positive scalar.");
            endif

          case "tolerance"
            Tolerance = varargin{2};
            if ( !(isscalar(Tolerance) && (Tolerance >= 0)))
              error ("ClassificationSVM: Tolerance must be a positive scalar.");
            endif

          case "shrinking"
            Shrinking = varargin{2};
            if ( !(ismember(Shrinking, [0, 1]) && isscalar(Shrinking)))
              error ("ClassificationSVM: Shrinking must be either 0 or 1.");
            endif

          case "probabilityestimates"
            ProbabilityEstimates = varargin{2};
            if ( !(ismember(ProbabilityEstimates, [0, 1]) && isscalar(ProbabilityEstimates)))
              error ( strcat(["ClassificationSVM: ProbabilityEstimates must be"], ...
             [" either 0 or 1."]));
            endif

          case "weight"

            Weight = varargin{2};
            weight_given = "yes";

            ## Check if weights is a structure
            if (!isstruct(Weight))
                error("ClassificationSVM: Weights must be provided as a structure.");
            endif

            ## Get the field names of the structure
            weightFields = fieldnames(Weight);

            ## Iterate over each field and validate
            for i = 1:length(weightFields)
                ## The field name must be a numeric string
                classLabel = str2double(weightFields{i});
                if (isnan(classLabel))
                    error("ClassificationSVM: Class labels in the weight structure must be numeric.");
                endif

                ## The field value must be a numeric scalar
                if (!(isnumeric(Weight.(weightFields{i})) && isscalar(Weight.(weightFields{i}))))
                    error("ClassificationSVM: Weights in the weight structure must be numeric scalars.");
                endif
            endfor

          case "kfold"
            KFold = varargin{2};
            if ( !(isnumeric(KFold) && isscalar(KFold) && (KFold > 1)
              && mod(KFold, 1) == 0 ))
              error ("ClassificationSVM: KFold must be a positive integer greater than 1.");
            endif

          otherwise
            error (strcat (["ClassificationSVM: invalid parameter name"],...
                           [" in optional pair arguments."]));

        endswitch
        varargin (1:2) = [];
      endwhile

      SVMtype = tolower(SVMtype);
        switch (SVMtype)
          case "c_svc"
            s = 0;
          case "nu_svc"
            s = 1;
          case "one_class_svm"
            s = 2;
        endswitch
      KernelFunction = tolower(KernelFunction);
        switch (KernelFunction)
          case "linear"
            t = 0;
          case "polynomial"
            t = 1;
          case "rbf"
            t = 2;
          case "sigmoid"
            t = 3;
        endswitch
      d = PolynomialOrder;
      g = Gamma;
      r = KernelOffset;
      c = BoxConstraint;
      n = Nu;
      m = CacheSize;
      e = Tolerance;
      h = Shrinking;
      b = ProbabilityEstimates;
      v = KFold;

    ## Assign properties
    this.X = X;
    this.Y = Y;

##    This is kept if the output of Model.parameter is not sufficient
##    this.ModelParameters = struct('SVMtype', SVMtype, ...
##                                 'KernelFunction', KernelFunction, ...
##                                 'PolynomialOrder', PolynomialOrder, ...
##                                 'Gamma', Gamma, ...
##                                 'KernelOffset', KernelOffset, ...
##                                 'BoxConstraint', BoxConstraint, ...
##                                 'Nu', Nu, ...
##                                 'CacheSize', CacheSize, ...
##                                 'Tolerance', Tolerance, ...
##                                 'Shrinking', Shrinking, ...
##                                 'ProbabilityEstimates', ProbabilityEstimates, ...
##                                 'Weight', Weight, ...
##                                 'KFold', KFold);

    ## svmpredict:
    ##    '-s':  SVMtype
    ##    '-t':  KernelFunction
    ##    '-d':  PolynomialOrder
    ##    '-g':  Gamma
    ##    '-r':  KernelOffset
    ##    '-c':  BoxConstraint
    ##    '-n':  Nu
    ##    '-m':  CacheSize
    ##    '-e':  Tolerance
    ##    '-h':  Shrinking
    ##    '-b':  ProbabilityEstimates
    ##    '-wi':  Weight
    ##    '-v':  KFold

    if (strcmp(weight_given, "yes"))
      ## Initialize an empty string
      weight_options = '';

      ## Get the field names of the structure
      weightFields = fieldnames(Weight);

      ## Iterate over each field and format it into the LIBSVM weight string
      for i = 1:length(weightFields)
          ## Convert the field name to a numeric value (class label)
          classLabel = str2double(weightFields{i});

          ## Get the corresponding weight
          wi = Weight.(weightFields{i});

          ## Append to the weight_options string
          weight_options = strcat(weight_options, sprintf(' -w%d %.2f ', classLabel, wi));
      endfor

      ## Remove the trailing space
      weight_options = strtrim(weight_options);

      ## Train the SVM model using svmtrain
      svm_options = sprintf(strcat(["-s %d -t %d -d %d -g %f -r %f -c %f -n %f"], ...
                                   [" -m %f -e %f -h %d -b %d %s -q"]), ...
                            s, t, d, g, r, c, n, m, e, h, b, weight_options);

    elseif (strcmp(weight_given, "no"))

      ## Train the SVM model using svmtrain
      svm_options = sprintf(strcat(["-s %d -t %d -d %d -g %f -r %f -c %f -n %f"], ...
                                  [" -m %f -e %f -h %d -b %d -q"]), ...
                            s, t, d, g, r, c, n, m, e, h, b);
    endif

##    disp(svm_options); ## For debugging

    svm_options_with_kfold = strcat(svm_options, sprintf(" -v %d ", v));

##    disp(svm_options_with_kfold); ## For debugging

    Model = svmtrain(Y, X, svm_options);

    this.Model =  Model;
    this.CrossValidationAccuracy = svmtrain(Y, X, svm_options_with_kfold);
    this.ModelParameters = Model.Parameters;
    this.NumClasses = Model.nr_class;
    this.SupportVectorCount = Model.totalSV;
    this.Rho = Model.rho;
    this.ClassNames = Model.Label;
    this.SupportVectorIndices = Model.sv_indices;
    this.ProbA = Model.ProbA;
    this.ProbB = Model.ProbB;
    this.SupportVectorPerClass = Model.nSV;
    this.SupportVectorCoef = Model.sv_coef;
    this.SupportVectors = Model.SVs;
    this.Solver = "SMO";

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationSVM} {[@var{label}, @var{decision_values}] =} predict (@dots{}, @var{name}, @var{value})
    ## @deftypefnx {ClassificationSVM} {[@var{label}, @var{prob_estimates}] =} predict (@dots{})
    ##
    ## Classify new data points into categories using the Support Vector Machine
    ## classification object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## predictor data in @code{obj.X} and corresponding labels, @code{obj.Y},
    ## stored in the Support Vector Machine classification model, @var{obj}.
    ## For one-class model, +1 or -1 is returned.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationSVM} class object.
    ##
    ## @item
    ## @var{XC} must be an @math{MxP} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the SVM model in
    ## @var{obj}.
    ##
    ## @end itemize
    ##
    ## @code{[@var{label}, @var{prob_estimates}] = predict (@var{obj}, @var{XC}, "ProbabilityEstimates", 1)}
    ## also returns @var{prob_estimates}, If k is the number of classes in
    ## training data, each row contains k values indicating the probability that
    ## the testing instance is in each class.
    ##
    ## @code{[@var{label}, @var{decision_values}] = predict (@var{obj}, @var{XC}, "ProbabilityEstimates", 0)}
    ## also returns @var{decision_values},  If k is the number of classes in
    ## training data, each row includes results of predicting k(k-1)/2
    ## binary-class SVMs.  For classification, k = 1 is a special case.
    ## Decision value +1 is returned for each testing instance, instead of
    ## an empty vector.
    ##
    ##
    ## @code{@var{label} = predict (@dots{}, @var{Name}, @var{Value})} returns the
    ## aforementioned results with additional properties specified by
    ## @qcode{Name-Value} pair arguments listed below.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"ProbabilityEstimates"} @tab @tab Specifies whether to
    ## output Probability Estimates or Decision Values. It accepts either
    ## 0 or 1. The default value is 0.
    ##
    ## @itemize
    ##
    ## @item
    ## 0 return decision values.
    ##
    ## @item
    ## 1 return probability estimates.
    ##
    ## @end itemize
    ##
    ## @end multitable
    ##
    ## @seealso{fitcsvm, ClassificationSVM}
    ## @end deftypefn

    function [label, value] = predict (this, XC, varargin)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("ClassificationSVM.predict: too few input arguments.");
      endif

      if (mod (nargin, 2) != 0)
        error ("ClassificationSVM.predict: Name-Value arguments must be in pairs.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("ClassificationSVM.predict: XC is empty.");
      elseif (columns (this.X) != columns (XC))
        error (strcat (["ClassificationSVM.predict: XC must have the same"], ...
                       [" number of features (columns) as in the SVM model."]));
      endif

      b = 0;  ## Default: return decision values.

      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "probabilityestimates"
            b = varargin{2};
            if ( !(ismember(b, [0, 1]) && isscalar(b)))
              error (strcat (["ClassificationSVM.predict: ProbabilityEstimates"], ...
                             [" must be either 1 or 0."]));
            endif

          otherwise
            error (strcat (["ClassificationSVM.predict: invalid parameter name"],...
                           [" in optional pair arguments."]));
          endswitch
        varargin (1:2) = [];
      endwhile

      predict_options = sprintf(strcat(["-b %d -q"]), b);

      [predict_label_L, accuracy_L, dec_values_L] = svmpredict(ones(rows(XC),1), XC, this.Model, predict_options);

      if (nargout > 0)
        label = predict_label_L;
        if (nargout > 1)
          value = dec_values_L;
        endif
      endif

    endfunction

   endmethods

endclassdef


%!demo
%! ## Create a Support Vector Machine classifier for Fisher's iris data.
%!
%! load fisheriris
%! X = meas;                   # Feature matrix
%! Y = species;                # Class labels
%! ## Convert species to numerical labels
%! ## 'setosa' -> 1, 'versicolor' -> 2, 'virginica' -> 3
%! Y = grp2idx(Y);
%!
%! ## Randomly partition the data into training and testing sets
%! cv = cvpartition(Y, 'HoldOut', 0.33);  # 33% data for testing, 67% for training
%!
%! X_train = X(training(cv), :);
%! Y_train = Y(training(cv));
%!
%! X_test = X(test(cv), :);
%! Y_test = Y(test(cv));
%!
%! svm_obj = fitcsvm(X_train, Y_train,"svmtype",'c_svc',"kernelfunction",'rbf');
%!
%! ## Predict the labels for the test set
%! predicted_labels = svm_obj.predict(X_test);
%!
%! ## Calculate the accuracy
%! accuracy = sum(predicted_labels == Y_test) / length(Y_test) * 100;
%! printf('Prediction Accuracy = %d%%\n', accuracy);

## Test input validation for constructor
%!error<ClassificationSVM: too few input arguments.> ClassificationSVM ()
%!error<ClassificationSVM: too few input arguments.> ClassificationSVM (ones(10,2))
%!error<ClassificationSVM: number of rows in X and Y must be equal.> ...
%! ClassificationSVM (ones(10,2), ones (5,1))
%!error<ClassificationSVM: Y must be a numeric array.> ...
%! ClassificationSVM (ones(5,2), ['A';'B';'A';'C';'B'])
%!error<ClassificationSVM: SVMtype must be a string.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "svmtype", 123)
%!error<ClassificationSVM: unsupported SVMtype.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "svmtype", "unsupported_type")
%!error<ClassificationSVM: KernelFunction must be a string.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "kernelfunction", 123)
%!error<ClassificationSVM: unsupported Kernel function.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "kernelfunction", "unsupported_function")
%!error<ClassificationSVM: PolynomialOrder must be a positive integer.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "polynomialorder", -1)
%!error<ClassificationSVM: PolynomialOrder must be a positive integer.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "polynomialorder", 0.5)
%!error<ClassificationSVM: PolynomialOrder must be a positive integer.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "polynomialorder", [1,2])
%!error<ClassificationSVM: Gamma must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "Gamma", -1)
%!error<ClassificationSVM: Gamma must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "Gamma", 0)
%!error<ClassificationSVM: Gamma must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "Gamma", [1, 2])
%!error<ClassificationSVM: Gamma must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "Gamma", "invalid")
%!error<ClassificationSVM: KernelOffset must be a non-negative scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "kerneloffset", -1)
%!error<ClassificationSVM: KernelOffset must be a non-negative scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "kerneloffset", [1,2])
%!error<ClassificationSVM: Nu must be a positive scalar in the range 0 < Nu <= 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "nu", -0.5)
%!error<ClassificationSVM: Nu must be a positive scalar in the range 0 < Nu <= 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "nu", 1.5)
%!error<ClassificationSVM: CacheSize must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "cachesize", -1)
%!error<ClassificationSVM: CacheSize must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "cachesize", [1,2])
%!error<ClassificationSVM: Tolerance must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "tolerance", -0.1)
%!error<ClassificationSVM: Tolerance must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "tolerance", [0.1,0.2])
%!error<ClassificationSVM: Shrinking must be either 0 or 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "shrinking", 2)
%!error<ClassificationSVM: Shrinking must be either 0 or 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "shrinking", -1)
%!error<ClassificationSVM: Shrinking must be either 0 or 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "shrinking", [1 0])
%!error<ClassificationSVM: ProbabilityEstimates must be either 0 or 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "probabilityestimates", 2)
%!error<ClassificationSVM: ProbabilityEstimates must be either 0 or 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "probabilityestimates", -1)
%!error<ClassificationSVM: ProbabilityEstimates must be either 0 or 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "probabilityestimates", [0 1])
%!error<ClassificationSVM: Weights must be provided as a structure.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "weight", 1)
%!error<ClassificationSVM: Class labels in the weight structure must be numeric.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "weight", struct('a',15, '2',7))
%!error<ClassificationSVM: Weights in the weight structure must be numeric scalars.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "weight", struct('1', 'a', '2',7))
%!error<ClassificationSVM: BoxConstraint must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "boxconstraint", -1)
%!error<ClassificationSVM: BoxConstraint must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "boxconstraint", [1,2])
%!error<ClassificationSVM: KFold must be a positive integer greater than 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "kfold", 1)
%!error<ClassificationSVM: KFold must be a positive integer greater than 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "kfold", -1)
%!error<ClassificationSVM: KFold must be a positive integer greater than 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "kfold", 0.5)
%!error<ClassificationSVM: KFold must be a positive integer greater than 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "kfold", [1,2])
%!error<ClassificationSVM: invalid parameter name in optional pair arguments.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "some", "some")

## Test input validation for predict method
%!error<ClassificationSVM.predict: too few input arguments.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)))
%!error<ClassificationSVM.predict: Name-Value arguments must be in pairs.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), zeros(2,2), "ProbabilityEstimates")
%!error<ClassificationSVM.predict: XC is empty.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), [])
%!error<ClassificationSVM.predict: XC must have the same number of features> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), 1)
%!error<ClassificationSVM.predict: ProbabilityEstimates must be either 1 or 0.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), zeros(2,2),"ProbabilityEstimates", "some")
%!error<ClassificationSVM.predict: ProbabilityEstimates must be either 1 or 0.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), zeros(2,2),"ProbabilityEstimates", 3)
%!error<ClassificationSVM.predict: ProbabilityEstimates must be either 1 or 0.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), zeros(2,2),"ProbabilityEstimates", [1 0])
%!error<ClassificationSVM.predict: invalid parameter name in optional pair arguments.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), zeros(2,2), "some", "some")
