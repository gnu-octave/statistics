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
## @item @tab @qcode{"obj.KernelOffset"} @tab It is the offset to be added
## to the kernel function's result. This parameter is used to shift the decision
## boundary of the SVM. It is particularly useful for certain kernel functions
## like the polynomial kernel, where it can help control the flexibility and
## complexity of the model. The value of `KernelOffset` is a numeric scalar.
##
## @item @tab @qcode{"obj.BoxConstraint"} @tab Stores the box constraint
## parameter, also known as C, which controls the trade-off between achieving
## a low training error and a low testing error that is, it regulates the
## balance between the margin size and classification error. Higher values of
## `BoxConstraint` aim to fit the training data more accurately, potentially at
## the risk of overfitting. The value of `BoxConstraint` is a positive numeric
## scalar.
##
## @item @tab @qcode{"obj.CacheSize"} @tab It is the size of the kernel
## cache in megabytes. The cache is used to store the kernel matrix, which can
## speed up the computation of the SVM by avoiding recalculations. A larger
## cache size can improve training speed, especially for large datasets, but
## it also increases memory usage. The value of `CacheSize` is a positive
## numeric scalar.
##
## @item @tab @qcode{"obj.SVMtype"} @tab Stores the type of Support Vector
## Machine (SVM) to be used for classification. The `SVMtype` property
## determines the formulation of the SVM, such as C-SVM or ν-SVM. Each type has
## its own characteristics and usage scenarios. The value of `SVMtype` is a
## string indicating the type of SVM.
##
## @item @tab @qcode{"obj.ClassNames"} @tab The labels for each class in the
## classification problem. It provides the actual names or identifiers for the
## classes being predicted by the model. This field is empty for one-class SVM
## because it only involves a single class during training and testing.
##
## @item @tab @qcode{"obj.NumObservations"} @tab Number of observations used in
## training the ClassificationSVM model, specified as a positive integer scalar.
## This number can be less than the number of rows in the training data because
## rows containing @qcode{NaN} values are not part of the fit.
##
## @item @tab @qcode{"obj.ResponseName"} @tab Response variable name, specified
## as a character vector.
##
## @item @tab @qcode{"obj.PredictorNames"} @tab Predictor variable names,
## specified as a cell array of character vectors.  The variable names are in
## the same order in which they appear in the training data @var{X}.
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
    NumObservations         = [];     # Number of Observations
    ResponseName            = [];     # Response variable name
    PredictorNames          = [];     # Predictor variable name

    Rho                     = [];     # Bias term
    KernelOffset            = [];     # Kernel Offset
    BoxConstraint           = [];     # Box Constraint
    CacheSize               = [];     # Cache Size
    SVMtype                 = [];     # Type of SVM

    SupportVectors          = [];     # Support vectors
    SupportVectorIndices    = [];     # Indices of Support vectors
    SupportVectorCount      = [];     # Total number of Support vectors
    SupportVectorPerClass   = [];     # Number of Support vectors for each class
    SupportVectorCoef       = [];     # Coefficients of support vectors in the
                                      # decision functions

    ProbA                   = [];     # Pairwise probability estimates
    ProbB                   = [];     # Pairwise probability estimates
    Solver                  = "SMO";  # Type of solver used

  endproperties

  properties (Access = private)

    Model                   = [];     # Stores the trained model
    KernelFunction          = [];     # Kernel Function
    PolynomialOrder         = [];     # Order of Polynomial in kernel
    Gamma                   = [];     # Gamma
    Nu                      = [];     # Nu
    Tolerance               = [];     # Tolerance for convergence
    Shrinking               = [];     # Shrinking
    ProbabilityEstimates    = [];     # Probability Estimates
    Weight                  = [];     # Weight of each predictor

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

      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);

      ## Set default values before parsing optional parameters
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
      NumObservations         = [];
      ResponseName            = [];
      PredictorNames          = [];

      weight_given            = "no";

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "svmtype"
            SVMtype = varargin{2};
            if (!(ischar(SVMtype)))
              error ("ClassificationSVM: SVMtype must be a string.");
            endif
            if (ischar(SVMtype))
              if (! any (strcmpi (tolower(SVMtype), {"c_svc", "nu_svc", ...
                "one_class_svm"})))
              error ("ClassificationSVM: unsupported SVMtype.");
              endif
            endif

          case "kernelfunction"
            KernelFunction = varargin{2};
            if (!(ischar(KernelFunction)))
              error ("ClassificationSVM: KernelFunction must be a string.");
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
              error (strcat (["ClassificationSVM: PolynomialOrder must be"], ...
                            [" a positive integer."]));
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
              error (strcat (["ClassificationSVM: KernelOffset must be a"], ...
                             [" non-negative scalar."]));
            endif

          case "boxconstraint"
            BoxConstraint = varargin{2};
            if ( !(isscalar(BoxConstraint) && (BoxConstraint > 0)))
              error (strcat (["ClassificationSVM: BoxConstraint must be a"], ...
                             [" positive scalar."]));
            endif

          case "nu"
            Nu = varargin{2};
            if ( !(isscalar(Nu) && (Nu > 0) && (Nu <= 1)))
              error (strcat (["ClassificationSVM: Nu must be a positive"], ...
                             [" scalar in the range 0 < Nu <= 1."]));
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
            if ( !(ismember(ProbabilityEstimates, [0, 1]) && ...
                   isscalar(ProbabilityEstimates)))
              error (strcat (["ClassificationSVM: ProbabilityEstimates"], ...
                             [" must be either 0 or 1."]));
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error ("ClassificationSVM: ResponseName must be a character array.");
            endif

          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat (["ClassificationSVM: PredictorNames must"], ...
                             [" be supplied as a cellstring array."]));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat (["ClassificationSVM: PredictorNames must"], ...
                             [" have the same number of columns as X."]));
            endif

          case "weight"

            Weight = varargin{2};
            weight_given = "yes";

            ## Check if weights is a structure
            if (!isstruct(Weight))
                error (strcat (["ClassificationSVM: Weights must be"], ...
                               [" provided as a structure."]));
            endif

            ## Get the field names of the structure
            weightFields = fieldnames(Weight);

            ## Iterate over each field and validate
            for i = 1:length(weightFields)
                ## The field name must be a numeric string
                classLabel = str2double(weightFields{i});
                if (isnan(classLabel))
                    error (strcat (["ClassificationSVM: Class labels in"], ...
                                   [" the weight structure must be numeric."]));
                endif

                ## The field value must be a numeric scalar
                if (!(isnumeric(Weight.(weightFields{i})) && ...
                      isscalar(Weight.(weightFields{i}))))
                    error (strcat (["ClassificationSVM: Weights in the"], ...
                           [" weight structure must be numeric scalars."]));
                endif
            endfor

          case "kfold"
            KFold = varargin{2};
            if ( !(isnumeric(KFold) && isscalar(KFold) && (KFold > 1) ...
              && mod(KFold, 1) == 0 ))
              error (strcat (["ClassificationSVM: KFold must be a"], ...
                             [" positive integer greater than 1."]));
            endif

          otherwise
            error (strcat (["ClassificationSVM: invalid parameter name"], ...
                           [" in optional pair arguments."]));

        endswitch
        varargin (1:2) = [];
      endwhile

      ## Remove missing values from X and Y
      RowsUsed  = ! logical (sum (isnan ([X, gY]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

      this.NumObservations = rows (X);

      ## Generate default predictors and response variabe names (if necessary)
      if (isempty (PredictorNames))
        for i = 1:ndims_X
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = "Y";
      endif

      ## Assign predictors and response variable names
      this.PredictorNames = PredictorNames;
      this.ResponseName   = ResponseName;

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
    this.KernelOffset            = KernelOffset;
    this.BoxConstraint           = BoxConstraint;
    this.CacheSize               = CacheSize;
    this.SVMtype                 = SVMtype;
    this.KernelFunction          = KernelFunction;
    this.PolynomialOrder         = PolynomialOrder;
    this.Gamma                   = Gamma;
    this.Nu                      = Nu;
    this.Tolerance               = Tolerance;
    this.Shrinking               = Shrinking;
    this.ProbabilityEstimates    = ProbabilityEstimates;
    this.Weight                  = Weight;

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
          weight_options = strcat(weight_options, sprintf(' -w%d %.2f ', ...
                                  classLabel, wi));
      endfor

      ## Remove the trailing space
      weight_options = strtrim(weight_options);

      ## Train the SVM model using svmtrain
      svm_options = sprintf(strcat(["-s %d -t %d -d %d -g %f -r %f -c %f"], ...
                                   [" -n %f -m %f -e %f -h %d -b %d %s -q"]), ...
                               s, t, d, g, r, c, n, m, e, h, b, weight_options);

    elseif (strcmp(weight_given, "no"))

      ## Train the SVM model using svmtrain
      svm_options = sprintf(strcat(["-s %d -t %d -d %d -g %f -r %f -c %f"], ...
                                  [" -n %f -m %f -e %f -h %d -b %d -q"]), ...
                               s, t, d, g, r, c, n, m, e, h, b);
    endif

    ##  disp(svm_options); ## For debugging

    svm_options_with_kfold = strcat(svm_options, sprintf(" -v %d ", v));

    ##  disp(svm_options_with_kfold); ## For debugging

    Model = svmtrain(Y, X, svm_options);

    this.Model =  Model;
    this.CrossValidationAccuracy = svmtrain(Y, X, svm_options_with_kfold);
    this.ModelParameters = struct(...
                                  'KernelFunction', this.KernelFunction,
                                  'PolynomialOrder', this.PolynomialOrder,
                                  'Gamma', this.Gamma,
                                  'Nu', this.Nu,
                                  'Tolerance', this.Tolerance,
                                  'Shrinking', this.Shrinking,
                                  'ProbabilityEstimates', this.ProbabilityEstimates,
                                  'Weight', this.Weight);
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
    ## @item
    ## @var{XC} must be an @math{MxP} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the SVM model in
    ## @var{obj}.
    ## @end itemize
    ##
    ## @code{[@var{label}, @var{prob_estimates}] = predict (@var{obj}, @var{XC},
    ## "ProbabilityEstimates", 1)}
    ## also returns @var{prob_estimates}, If k is the number of classes in
    ## training data, each row contains k values indicating the probability that
    ## the testing instance is in each class.
    ##
    ## @code{[@var{label}, @var{decision_values}] = predict (@var{obj}, @var{XC},
    ## "ProbabilityEstimates", 0)}
    ## also returns @var{decision_values},  If k is the number of classes in
    ## training data, each row includes results of predicting k(k-1)/2
    ## binary-class SVMs.  For classification, k = 1 is a special case.
    ## Decision value +1 is returned for each testing instance, instead of
    ## an empty vector.
    ##
    ##
    ## @code{@var{label} = predict (@dots{}, @var{Name}, @var{Value})} returns
    ## the aforementioned results with additional properties specified by
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
    ## @item
    ## 0 return decision values.
    ## @item
    ## 1 return probability estimates.
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
        error (strcat (["ClassificationSVM.predict: Name-Value arguments"], ...
                       [" must be in pairs."]));
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

      predict_options = sprintf("-b %d -q", b);

      [predict_label_L, ~, dec_values_L] = svmpredict(ones(rows(XC),1), XC, ...
                                                   this.Model, predict_options);

      if (nargout > 0)
        label = predict_label_L;
        if (nargout > 1)
          value = dec_values_L;
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
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
    ## @var{obj} must be a binary class @qcode{ClassificationSVM} object.
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
    ## @seealso{fitcsvm, ClassificationSVM}
    ## @end deftypefn

    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationSVM.margin: too few input arguments.");
      endif

      ## Check if binary classifier model or not.
      if (this.NumClasses != 2)
        error (strcat(["ClassificationSVM.margin: only binary classifier SVM"], ...
                      [" model is supported."]));
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationSVM.margin: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat (["ClassificationSVM.margin: X must have the same"], ...
                       [" number of features (columns) as in the SVM model."]));
      endif

      ## Check for valid Y
      if (isempty (Y))
        error ("ClassificationSVM.margin: Y is empty.");
      elseif (rows (X)!= rows (Y))
        error (strcat (["ClassificationSVM.margin: Y must have the same"], ...
                       [" number of rows as X."]));
      endif

      [~, ~, dec_values_L] = svmpredict(Y, X, this.Model, '-q');
      m = 2 * Y .* dec_values_L;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationSVM} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
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
    ## @var{obj} must be a binary class @qcode{ClassificationSVM} object.
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
        error ("ClassificationSVM.loss: too few input arguments.");
      endif

      if (mod (nargin, 2) == 0)
        error ("ClassificationSVM.loss: Name-Value arguments must be in pairs.");
      endif

      ## Check if binary classifier model or not.
      if (this.NumClasses != 2)
        error (strcat (["ClassificationSVM.loss: only binary classifier SVM"], ...
                       [" model is supported."]));
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationSVM.loss: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat (["ClassificationSVM.loss: X must have the same"], ...
                       [" number of features (columns) as in the SVM model."]));
      endif

      ## Check for valid Y
      if (isempty (Y))
        error ("ClassificationSVM.loss: Y is empty.");
      elseif (rows (X)!= rows (Y))
        error (strcat (["ClassificationSVM.loss: Y must have the same"], ...
                       [" number of rows as X."]));
      endif

      ## Set default values before parsing optional parameters
      LossFun                 = 'classiferror';
      Weights                 = ones(size(X,1),1);

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "lossfun"
            LossFun = varargin{2};
            if (!(ischar(LossFun)))
              error ("ClassificationSVM.loss: LossFun must be a string.");
            endif
            if (ischar(LossFun))
              if (! any (strcmpi (tolower(LossFun), {"binodeviance",  ...
                "classiferror", "exponential", "hinge", "logit", "quadratic"})))
              error ("ClassificationSVM.loss: unsupported Loss function.");
              endif
            endif

          case "weights"
            Weights = varargin{2};
            ## Validate if weights is a numeric vector
            if(!(isnumeric(Weights) && isvector(Weights)))
              error (strcat (["ClassificationSVM.loss: Weights must be a"], ...
                             [" numeric vector."]));
            endif

            ## Check if the size of weights matches the number of rows in X
            if (numel(Weights) != size(X, 1))
              error(strcat (["ClassificationSVM.loss: size of Weights"], ...
                            [" must be equal to the number of rows in X."]));
            endif

           otherwise
            error (strcat (["ClassificationSVM.loss: invalid parameter"], ...
                           [" name in optional pair arguments."]));
          endswitch
        varargin (1:2) = [];
      endwhile

      ## Compute the classification score
      [~, ~, dec_values_L] = svmpredict(Y, X, this.Model, '-q');

        ## Compute the margin
        margin = Y .* dec_values_L;

        ## Compute the loss based on the specified loss function
        switch tolower(LossFun)
          case "classiferror"
            L = mean((margin <= 0) .* Weights);

          case "hinge"
            L = mean(max(0, 1 - margin) .* Weights);

          case "logit"
            L = mean(log(1 + exp(-margin)) .* Weights);

          case "exponential"
            L = mean(exp(-margin) .* Weights);

          case "quadratic"
            L = mean((1 - margin).^2 .* Weights);

          case "binodeviance"
            L = mean(log(1 + exp(-2 * margin)) .* Weights);

          otherwise
            error("ClassificationSVM.loss: unsupported Loss function.");
        endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{label} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationSVM} {[@var{label}, @var{decision_values}] =} resubPredict (@dots{}, @var{name}, @var{value})
    ## @deftypefnx {ClassificationSVM} {[@var{label}, @var{prob_estimates}] =} resubPredict (@dots{})
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
    ## @code{[@var{label}, @var{prob_estimates}] = resubPredict (@var{obj},
    ## "ProbabilityEstimates", 1)}
    ## also returns @var{prob_estimates}. If k is the number of classes in the
    ## training data, each row contains k values indicating the probability that
    ## the training instance is in each class.
    ##
    ## @code{[@var{label}, @var{decision_values}] = resubPredict (@var{obj},
    ## "ProbabilityEstimates", 0)}
    ## also returns @var{decision_values}.  If k is the number of classes in the
    ## training data, each row includes results of predicting k(k-1)/2
    ## binary-class SVMs.  For classification, k = 1 is a special case.
    ## Decision value +1 is returned for each training instance, instead of
    ## an empty vector.
    ##
    ##
    ## @code{@var{label} = resubPredict (@dots{}, @var{Name}, @var{Value})}
    ## returns the aforementioned results with additional properties specified
    ## by @qcode{Name-Value} pair arguments listed below.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"ProbabilityEstimates"} @tab @tab Specifies whether to
    ## output Probability Estimates or Decision Values. It accepts either
    ## 0 or 1. The default value is 0.
    ##
    ## @itemize
    ## @item
    ## 0 return decision values.
    ## @item
    ## 1 return probability estimates.
    ## @end itemize
    ##
    ## @end multitable
    ##
    ## @seealso{fitcsvm, ClassificationSVM}
    ## @end deftypefn

    function [label, value] = resubPredict (this, varargin)

      if (mod (nargin, 2) != 1)
        error (strcat (["ClassificationSVM.resubPredict: Name-Value"], ...
                       [" arguments must be in pairs."]));
      endif

      b = 0;  ## Default: return decision values.

      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "probabilityestimates"
            b = varargin{2};
            if ( !(ismember(b, [0, 1]) && isscalar(b)))
              error (strcat (["ClassificationSVM.resubPredict:"], ...
                             [" ProbabilityEstimates must be either 1 or 0."]));
            endif

          otherwise
            error (strcat (["ClassificationSVM.resubPredict: invalid"],...
                           [" parameter name in optional pair arguments."]));
          endswitch
        varargin (1:2) = [];
      endwhile

      predict_options = sprintf("-b %d -q", b);

      [predict_label_L, ~, dec_values_L] = svmpredict(this.Y, this.X, ...
                                                   this.Model, predict_options);

      if (nargout > 0)
        label = predict_label_L;
        if (nargout > 1)
          value = dec_values_L;
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {ClassificationSVM} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the resubstitution classification loss for the trained Support
    ## Vector Machine classification object.
    ##
    ## @code{@var{L} = resubLoss (@var{obj})} returns the classification loss by
    ## resubstitution (L), or the in-sample classification loss, for the trained
    ## classification model @var{obj} using the training data stored in
    ## @code{obj.X} and the corresponding class labels stored in @code{obj.Y}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a binary class @qcode{ClassificationSVM} object.
    ## @end itemize
    ##
    ## @code{@var{l} = resubLoss (@dots{}, @var{Name}, @var{Value})} returns the
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

    function L = resubLoss(this, varargin)

      if (mod(nargin, 2) != 1)
        error (strcat (["ClassificationSVM.resubLoss: Name-Value arguments"], ...
                       [" must be in pairs."]));
      endif

      ## Check if binary classifier model or not.
      if (this.NumClasses != 2)
        error (strcat(["ClassificationSVM.resubloss: only binary"], ...
                      [" classifier SVM model is supported."]));
      endif

      ## Set default values before parsing optional parameters
      LossFun                        = 'classiferror';
      Weights                        = ones(size(this.X, 1), 1);

      ## Parse extra parameters
      while (numel(varargin) > 0)
        switch (tolower(varargin{1}))

          case "lossfun"
            LossFun = varargin{2};
            if (!ischar(LossFun))
              error("ClassificationSVM.resubLoss: LossFun must be a string.");
            endif
            if (ischar(LossFun))
              if (!any(strcmpi(tolower(LossFun), {"binodeviance", ...
                "classiferror", "exponential", "hinge", "logit", "quadratic"})))
                error (strcat (["ClassificationSVM.resubLoss: unsupported"], ...
                               [" Loss function."]));
              endif
            endif

          case "weights"
            Weights = varargin{2};
            ## Validate if weights is a numeric vector
            if (!(isnumeric(Weights) && isvector(Weights)))
              error (strcat (["ClassificationSVM.resubLoss: Weights must"], ...
                             [" be a numeric vector."]));
            endif

            ## Check if the size of weights matches the number of rows in X
            if (numel(Weights) != size(this.X, 1))
              error (strcat (["ClassificationSVM.resubLoss: size of Weights"], ...
                             [" must be equal to the number of rows in X."]));
            endif

          otherwise
            error (strcat (["ClassificationSVM.resubLoss: invalid parameter"], ...
                           [" name in optional pair arguments."]));
        endswitch
        varargin(1:2) = [];
      endwhile

      ## Compute the classification score
      [~, ~, dec_values_L] = svmpredict(this.Y, this.X, this.Model, '-q');

      ## Compute the margin
      margin = this.Y .* dec_values_L;

        ## Compute the loss based on the specified loss function
        switch tolower(LossFun)
          case "classiferror"
            L = mean((margin <= 0) .* Weights);

          case "hinge"
            L = mean(max(0, 1 - margin) .* Weights);

          case "logit"
            L = mean(log(1 + exp(-margin)) .* Weights);

          case "exponential"
            L = mean(exp(-margin) .* Weights);

          case "quadratic"
            L = mean((1 - margin).^2 .* Weights);

          case "binodeviance"
            L = mean(log(1 + exp(-2 * margin)) .* Weights);

          otherwise
            error ("ClassificationSVM.resubloss: unsupported Loss function.");
        endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{m} =} resubMargin (@var{obj})
    ##
    ## Determine the resubstitution margins for a Support Vector Machine
    ## classification object.
    ##
    ## @code{@var{m} = resubMargin (@var{obj})} returns the resubstitution
    ## classification margins for the trained support vector machine (SVM)
    ## classifier @var{obj} using the training data stored in @code{obj.X} and
    ## the corresponding class labels stored in @code{obj.Y}. It supports only
    ## binary classifier models. The classification margin is commonly defined
    ## as @var{m} = @var{y}f(@var{x}), where @var{f(x)} is the classification
    ## score and @var{y} is the true class label corresponding to @var{x}. A
    ## greater margin indicates a better model.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a binary class @qcode{ClassificationSVM} object.
    ## @end itemize
    ##
    ## @seealso{fitcsvm, ClassificationSVM}
    ## @end deftypefn

    function m = resubMargin (this)

      ## Check if binary classifier model or not.
      if (this.NumClasses != 2)
        error (strcat(["ClassificationSVM.resubMargin: only binary"], ...
                      [" classifier SVM model is supported."]));
      endif

      ## Get the decision values for the training data
      [~, ~, dec_values_L] = svmpredict(this.Y, this.X, this.Model, '-q');
      m = this.Y .* dec_values_L;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationSVM} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross Validates Support Vector Machine classification object.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} returns a cross-validated model
    ## object, @var{CVMdl}, from a trained model, @var{obj}, using 10-fold
    ## cross-validation by default.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj}, @var{name}, @var{value})}
    ## specifies additional name-value pair arguments to customise the
    ## cross-validation process.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{"KFold"} @tab @tab Specify the number of folds to use in
    ## k-fold cross-validation. @code{'kfold', @var{k}} where @var{k} is an
    ## integer greater than 1.
    ##
    ## @item @qcode{"Holdout"} @tab @tab Specify the fraction of the data to
    ## hold out for testing. @code{'holdout', @var{p}} where @var{p} is a scalar
    ## in the range (0,1).
    ##
    ## @item @qcode{"CVPartition"} @tab @tab Specify the fraction of the data to
    ## hold out for testing. @code{'holdout', @var{p}} where @var{p} is a scalar
    ## in the range (0,1).
    ##
    ## @item @qcode{"Leaveout"} @tab @tab Specify whether to perform
    ## leave-one-out cross-validation. @code{'leaveout', @var{Value}} where
    ## @var{Value} is 'on' or 'off'.
    ##
    ## @end multitable
    ##
    ## @seealso{cvpartition, fitcsvm, ClassificationSVM}
    ## @end deftypefn

    function CVMdl = crossval (this, varargin)

      ## Check for sufficient input arguments
      if (nargin < 1)
        error ("ClassificationSVM.crossval: too few input arguments.");
      endif

      if (numel (varargin) == 1)
        error (strcat (["ClassificationSVM.crossval: Name-Value arguments"], ...
                       [" must be in pairs."]));
      elseif (numel (varargin) > 2)
        error (strcat (["ClassificationSVM.crossval: specify only one of"], ...
                       [" the Name-Value arguments."]));
      endif

      ## Set default values before parsing optional parameters
      numSamples  = size (this.X, 1);
      numFolds    = 10;
      Holdout     = [];
      Leaveout    = 'off';
      CVPartition = [];

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'kfold'
            numFolds = varargin{2};
            if (!((isnumeric (numFolds) && isscalar (numFolds)
              && (numFolds == fix (numFolds)) && numFolds > 1 )))
              error (strcat (["ClassificationSVM.crossval: KFold"],...
             [" should be an integer value greater than 1."]));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (!(isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
              && Holdout < 1))
              error (strcat (["ClassificationSVM.crossval: Holdout should"], ...
                             [" be a numeric value in the range 0 to 1."]));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (!(ischar (Leaveout)
              && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat (["ClassificationSVM.crossval: Leaveout should"], ...
                             [" be either 'on' or 'off'."]));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (!(isa (CVPartition, 'cvpartition')))
              error (strcat (["ClassificationSVM.crossval: CVPartition"],...
                             [" should be a cvPartition object."]));
            endif

           otherwise
            error (strcat (["ClassificationSVM.crossval: invalid"],...
                           [" parameter name in optional pair arguments."]));
          endswitch
        varargin (1:2) = [];
      endwhile

      ## Determine the cross-validation method to use
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (numSamples, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (numSamples, 'LeaveOut');
      else
        partition = cvpartition (numSamples, 'KFold', numFolds);
      endif

      ## Create a cross-validated model object
      CVMdl = ClassificationPartitionedModel (this, partition);

      endfunction

   endmethods

endclassdef


%!demo
%! ## Create a Support Vector Machine classifier for Fisher's iris data and
%! ## predict labels for test data.
%!
%! load fisheriris
%! X = meas;                   # Feature matrix
%! Y = species;                # Class labels
%! ## Convert species to numerical labels
%! ## 'setosa' -> 1, 'versicolor' -> 2, 'virginica' -> 3
%! Y = grp2idx(Y);
%!
%! rng(1); ## For reproducibility
%!
%! ## Randomly partition the data into training and testing sets
%! cv = cvpartition(Y, 'HoldOut', 0.3); # 30% data for testing, 60% for training
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

%!demo
%! ## Create a Support Vector Machine classifier for Fisher's iris data and plot
%! ## the support vectors.
%!
%! load fisheriris
%! inds = !strcmp(species,'setosa');
%! X = meas(inds,3:4);              # Feature matrix
%! Y = grp2idx(species(inds));      # Class Labels
%!
%! SVMModel = fitcsvm(X,Y)
%!
%! sv = SVMModel.SupportVectors;
%! figure
%! gscatter(X(:,1),X(:,2),Y)
%! hold on
%! plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
%! legend('versicolor','virginica','Support Vector')
%! hold off

%!demo
%! ## Create a Support Vector Machine classifier and determine margin for test
%! ## data.
%! pkg load statistics
%! load fisheriris
%! rng(1);  ## For reproducibility
%!
%! ## Select indices of the non-setosa species
%! inds = !strcmp(species, 'setosa');
%!
%!  ## Select features and labels for non-setosa species
%! X = meas(inds, 3:4);
%! Y = grp2idx(species(inds));
%!
%! ##  Convert labels to +1 and -1
%! unique_classes = unique(Y);
%! Y(Y == unique_classes(1)) = -1;
%! Y(Y == unique_classes(2)) = 1;
%!
%! ## Partition data for training and testing
%! cv = cvpartition(Y, 'HoldOut', 0.15);
%! X_train = X(training(cv), :);
%! Y_train = Y(training(cv));
%! X_test = X(test(cv), :);
%! Y_test = Y(test(cv));
%!
%! ## Train the SVM model
%! CVSVMModel = fitcsvm(X_train, Y_train);
%!
%! ## Calculate margins
%! m = margin(CVSVMModel, X_test, Y_test);
%! disp(m);

%!demo
%! ## Create a Support Vector Machine classifier and determine loss for test
%! ## data.
%! pkg load statistics
%! load fisheriris
%! rng(1);  ## For reproducibility
%!
%!  ## Select indices of the non-setosa species
%! inds = !strcmp(species, 'setosa');
%!
%!  ## Select features and labels for non-setosa species
%! X = meas(inds, 3:4);
%! Y = grp2idx(species(inds));
%!
%! ##  Convert labels to +1 and -1
%! unique_classes = unique(Y);
%! Y(Y == unique_classes(1)) = -1;
%! Y(Y == unique_classes(2)) = 1;
%!
%! ## Randomly partition the data into training and testing sets
%! cv = cvpartition(Y, 'HoldOut', 0.3); # 30% data for testing, 60% for training
%!
%! X_train = X(training(cv), :);
%! Y_train = Y(training(cv));
%!
%! X_test = X(test(cv), :);
%! Y_test = Y(test(cv));
%!
%! ## Train the SVM model
%! SVMModel = fitcsvm(X_train, Y_train);
%!
%! ## Calculate loss
%!
%! L = loss(SVMModel,X_test,Y_test,'LossFun','binodeviance')
%! L = loss(SVMModel,X_test,Y_test,'LossFun','classiferror')
%! L = loss(SVMModel,X_test,Y_test,'LossFun','exponential')
%! L = loss(SVMModel,X_test,Y_test,'LossFun','hinge')
%! L = loss(SVMModel,X_test,Y_test,'LossFun','logit')
%! L = loss(SVMModel,X_test,Y_test,'LossFun','quadratic')

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
%!error<ClassificationSVM: PredictorNames must be supplied as a cellstring array.> ...
%! ClassificationSVM (ones (5,2), ones (5,1), "PredictorNames", ["A"])
%!error<ClassificationSVM: PredictorNames must be supplied as a cellstring array.> ...
%! ClassificationSVM (ones (5,2), ones (5,1), "PredictorNames", "A")
%!error<ClassificationSVM: PredictorNames must have the same number of columns as X.> ...
%! ClassificationSVM (ones (5,2), ones (5,1), "PredictorNames", {"A", "B", "C"})
%!error<ClassificationSVM: ResponseName must be a character array.> ...
%! ClassificationSVM (ones (5,2), ones (5,1), "ResponseName", {"Y"})
%!error<ClassificationSVM: ResponseName must be a character array.> ...
%! ClassificationSVM (ones (5,2), ones (5,1), "ResponseName", 1)
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

## Test input validation for margin method
%!error<ClassificationSVM.margin: too few input arguments.> ...
%! margin (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)))
%!error<ClassificationSVM.margin: too few input arguments.> ...
%! margin (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2))
%!error<ClassificationSVM.margin: only binary classifier SVM model is supported.> ...
%! margin (ClassificationSVM (ones (40,2),randi([1, 3], 40, 1)), zeros(2,2), ones(2,1))
%!error<ClassificationSVM.margin: X is empty.> ...
%! margin (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), [], zeros(2,2))
%!error<ClassificationSVM.margin: X must have the same number of features> ...
%! margin (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), 1, zeros(2,2))
%!error<ClassificationSVM.margin: Y is empty.> ...
%! margin (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), [])
%!error<ClassificationSVM.margin: Y must have the same number of rows as X.> ...
%! margin (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), 1)

## Test input validation for loss method
%!error<ClassificationSVM.loss: too few input arguments.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)))
%!error<ClassificationSVM.loss: too few input arguments.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2))
%!error<ClassificationSVM.loss: Name-Value arguments must be in pairs.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), ones(2,1), "LossFun")
%!error<ClassificationSVM.loss: only binary classifier SVM model is supported.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 3], 40, 1)), zeros(2,2), ones(2,1))
%!error<ClassificationSVM.loss: X is empty.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), [], zeros(2,2))
%!error<ClassificationSVM.loss: X must have the same number of features> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), 1, zeros(2,2))
%!error<ClassificationSVM.loss: Y is empty.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), [])
%!error<ClassificationSVM.loss: Y must have the same number of rows as X.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), 1)
%!error<ClassificationSVM.loss: LossFun must be a string.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), ones(2,1), "LossFun", 1)
%!error<ClassificationSVM.loss: unsupported Loss function.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), ones(2,1), "LossFun", "some")
%!error<ClassificationSVM.loss: Weights must be a numeric vector.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), ones(2,1), "Weights", ['a','b'])
%!error<ClassificationSVM.loss: Weights must be a numeric vector.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), ones(2,1), "Weights", 'a')
%!error<ClassificationSVM.loss: size of Weights must be equal to the number of rows in X.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), ones(2,1), "Weights", [1,2,3])
%!error<ClassificationSVM.loss: size of Weights must be equal to the number of rows in X.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), ones(2,1), "Weights", 3)
%!error<ClassificationSVM.loss: invalid parameter name in optional pair arguments.> ...
%! loss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), zeros(2,2), ones(2,1), "some", "some")

## Test input validation for resubPredict method
%!error<ClassificationSVM.resubPredict: Name-Value arguments must be in pairs.> ...
%! resubPredict (ClassificationSVM (ones (40,2), ones (40,1)), "ProbabilityEstimates")
%!error<ClassificationSVM.resubPredict: ProbabilityEstimates must be either 1 or 0.> ...
%! resubPredict (ClassificationSVM (ones (40,2), ones (40,1)), "ProbabilityEstimates", "some")
%!error<ClassificationSVM.resubPredict: ProbabilityEstimates must be either 1 or 0.> ...
%! resubPredict (ClassificationSVM (ones (40,2), ones (40,1)), "ProbabilityEstimates", 3)
%!error<ClassificationSVM.resubPredict: ProbabilityEstimates must be either 1 or 0.> ...
%! resubPredict (ClassificationSVM (ones (40,2), ones (40,1)), "ProbabilityEstimates", [1 0])
%!error<ClassificationSVM.resubPredict: invalid parameter name in optional pair arguments.> ...
%! resubPredict (ClassificationSVM (ones (40,2), ones (40,1)), "some", "some")

## Test input validation for resubLoss method
%!error<ClassificationSVM.resubLoss: Name-Value arguments must be in pairs.> ...
%! resubLoss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "LossFun")
%!error<ClassificationSVM.resubloss: only binary classifier SVM model is supported.> ...
%! resubLoss (ClassificationSVM (ones (40,2),randi([1, 3], 40, 1)))
%!error<ClassificationSVM.resubLoss: LossFun must be a string.> ...
%! resubLoss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "LossFun", 1)
%!error<ClassificationSVM.resubLoss: unsupported Loss function.> ...
%! resubLoss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "LossFun", "some")
%!error<ClassificationSVM.resubLoss: Weights must be a numeric vector.> ...
%! resubLoss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Weights", ['a','b'])
%!error<ClassificationSVM.resubLoss: Weights must be a numeric vector.> ...
%! resubLoss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Weights", 'a')
%!error<ClassificationSVM.resubLoss: size of Weights must be equal to the number of rows in X.> ...
%! resubLoss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Weights", [1,2,3])
%!error<ClassificationSVM.resubLoss: size of Weights must be equal to the number of rows in X.> ...
%! resubLoss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Weights", 3)
%!error<ClassificationSVM.resubLoss: invalid parameter name in optional pair arguments.> ...
%! resubLoss (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "some", "some")

## Test input validation for resubMargin method
%!error<ClassificationSVM.resubMargin: only binary classifier SVM model is supported.> ...
%! resubMargin (ClassificationSVM (ones (40,2),randi([1, 3], 40, 1)))

## Test input validation for crossval method
%!error<ClassificationSVM.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold")
%!error<ClassificationSVM.crossval: specify only one of the Name-Value arguments.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", 5, "leaveout", 'on')
%!error<ClassificationSVM.crossval: KFold should be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", 'a')
%!error<ClassificationSVM.crossval: KFold should be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", 1)
%!error<ClassificationSVM.crossval: KFold should be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", -1)
%!error<ClassificationSVM.crossval: KFold should be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", 11.5)
%!error<ClassificationSVM.crossval: KFold should be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", [1,2])
%!error<ClassificationSVM.crossval: Holdout should be a numeric value in the range 0 to 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 'a')
%!error<ClassificationSVM.crossval: Holdout should be a numeric value in the range 0 to 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 11.5)
%!error<ClassificationSVM.crossval: Holdout should be a numeric value in the range 0 to 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Holdout", -1)
%!error<ClassificationSVM.crossval: Holdout should be a numeric value in the range 0 to 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 0)
%!error<ClassificationSVM.crossval: Holdout should be a numeric value in the range 0 to 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 1)
%!error<ClassificationSVM.crossval: Leaveout should be either 'on' or 'off'.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Leaveout", 1)
%!error<ClassificationSVM.crossval: CVPartition should be a cvPartition object.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "CVPartition", 1)
%!error<ClassificationSVM.crossval: CVPartition should be a cvPartition object.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "CVPartition", 'a')
%!error<ClassificationSVM.crossval: invalid parameter name in optional pair arguments.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "some", "some")
