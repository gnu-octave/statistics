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
## Machine classification model for one-class or two-class problems.
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
## corresponding predictor data in @var{X}.  @var{Y} can be either numeric,
## logical, or cell array of character vectors.  It must have same numbers of
## rows as @var{X}.
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
## @multitable @columnfractions 0.23 0.02 0.75
## @headitem @var{Field} @tab @tab @var{Description}
##
## @item @qcode{X} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{Y} @tab @tab Class labels, specified as a logical or
## numeric vector, or cell array of character vectors.  Each value in @var{Y} is
## the observed class label for the corresponding row in @var{X}.
##
## @item  @qcode{NumObservations} @tab@tab Number of observations used in
## training the ClassificationSVM model, specified as a positive integer scalar.
## This number can be less than the number of rows in the training data because
## rows containing @qcode{NaN} values are not part of the fit.
##
## @item @qcode{RowsUsed} @tab @tab Rows of the original training data
## used in fitting the ClassificationSVM model, specified as a numerical vector.
## If you want to use this vector for indexing the training data in @var{X}, you
## have to convert it to a logical vector, i.e
## @qcode{X = obj.X(logical (obj.RowsUsed), :);}
##
## @item @qcode{Standardize} @tab @tab A boolean flag indicating whether
## the data in @var{X} have been standardized prior to training.
##
## @item @qcode{Sigma} @tab @tab Predictor standard deviations, specified
## as a numeric vector of the same length as the columns in @var{X}.  If the
## predictor variables have not been standardized, then @qcode{"obj.Sigma"} is
## empty.
##
## @item @qcode{Mu} @tab @tab Predictor means, specified as a numeric
## vector of the same length as the columns in @var{X}.  If the predictor
## variables have not been standardized, then @qcode{"obj.Mu"} is empty.
##
## @item @qcode{NumPredictors} @tab @tab The number of predictors
## (variables) in @var{X}.
##
## @item @qcode{PredictorNames} @tab @tab Predictor variable names,
## specified as a cell array of character vectors.  The variable names are in
## the same order in which they appear in the training data @var{X}.
##
## @item @qcode{ResponseName} @tab @tab Response variable name, specified
## as a character vector.
##
## @item @qcode{ClassNames} @tab @tab Names of the classes in the class
## labels, @var{Y}, used for fitting the SVM model.  @qcode{ClassNames} are of
## the same type as the class labels in @var{Y}.
##
## @item @qcode{Prior} @tab @tab Prior probabilities for each class,
## specified as a numeric vector.  The order of the elements in @qcode{Prior}
## corresponds to the order of the classes in @qcode{ClassNames}.
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
## @item @qcode{ScoreTransform} @tab @tab A function_handle which is used
## for transforming the SVM prediction score into a posterior probability.  By
## default, it is @qcode{'none'}, in which case the @code{predict} and
## @code{resubPredict} methods return the prediction scores.  Use the
## @code{fitPosterior} method to compute the appropriate @code{ScoreTransform},
## in which case the @code{predict} and @code{resubPredict} methods return the
## posterior probabilities.
##
## @item @qcode{ModelParameters} @tab @tab  A structure containing the
## parameters used to train the SVM model with the following fields:
## @code{SVMtype}, @code{BoxConstraint}, @code{CacheSize}, @code{KernelScale},
## @code{KernelOffset}, @code{KernelFunction}, @code{PolynomialOrder},
## @code{Nu}, @code{Tolerance}, and @code{Shrinking}.  Type @code{help fitcsvm}
## for more info on their usage and default values.
##
## @item @qcode{Model} @tab @tab  A structure containing the trained model in
## @qcode{'libsvm'} format.
##
## @item @qcode{Alpha} @tab @tab The coefficients of the trained SVM
## classifier specified as an @math{sx1} numeric vector, where @math{s} is the
## number of support vectors equal to @qcode{sum (obj.IsSupportVector)}.  If the
## SVM classifier was trained with a @qcode{'linear'} kernel function, then
## @qcode{obj.Alpha} is left empty.
##
## @item @qcode{Beta} @tab @tab The linear predictor coefficients specified
## as an @math{sx1} numeric vector, where @math{s} is the number of support
## vectors equal to @qcode{sum (obj.IsSupportVector)}.  If the SVM classifier
## was trained with a kernel function other than @qcode{'linear'}, then
## @qcode{obj.Beta} is left empty.
##
## @item @qcode{Bias} @tab @tab The bias term specified as a scalar.
##
## @item @qcode{IsSupportVector} @tab @tab Support vector indicator,
## specified as an @math{Nx1} logical vector that flags whether a corresponding
## observation in the predictor data matrix is a Support Vector.  @math{N} is
## the number of observations in the training data (see @code{NumObservations}).
##
## @item @qcode{SupportVectorLabels} @tab @tab The support vector class
## labels specified as an @math{sx1} numeric vector, where @math{s} is the
## number of support vectors equal to @qcode{sum (obj.IsSupportVector)}.  A
## value of +1 in @code{SupportVectorLabels} indicates that the corresponding
## support vector belongs to the positive class @qcode{(ClassNames@{2@})}.  A
## value of -1 indicates that the corresponding support vector belongs to the
## negative class @qcode{(ClassNames@{1@})}.
##
## @item @qcode{SupportVectors} @tab @tab The support vectors of the
## trained SVM classifier specified an @math{sxp} numeric matrix, where @math{s}
## is the number of support vectors equal to @qcode{sum (obj.IsSupportVector)},
## and @math{p} is the number of predictor variables in the predictor data.
##
## @end multitable
##
## @seealso{fitcsvm, svmtrain, svmpredict}
## @end deftypefn

  properties (Access = public)

    X                   = [];    # Predictor data
    Y                   = [];    # Class labels

    NumObservations     = [];    # Number of observations in training dataset
    RowsUsed            = [];    # Rows used in fitting
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

  methods (Access = public)

    ## Class object constructor
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
      CacheSize               = 1000;
      Tolerance               = 1e-6;
      Shrinking               = 1;
      Standardize             = false;
      ResponseName            = [];
      PredictorNames          = [];
      ClassNames              = [];
      Prior                   = [];
      Cost                    = [];
      this.ScoreTransform     = 'none';

      ## Parse extra parameters
      SVMtype_override = true;
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "standardize"
            Standardize = varargin{2};
            if (! (Standardize == true || Standardize == false))
              error (strcat (["ClassificationSVM: 'Standardize' must"], ...
                             [" be either true or false."]));
            endif

          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat (["ClassificationSVM: 'PredictorNames' must"], ...
                             [" be supplied as a cellstring array."]));
            elseif (columns (PredictorNames) != columns (X))
              error (strcat (["ClassificationSVM: 'PredictorNames' must"], ...
                             [" have the same number of columns as X."]));
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat (["ClassificationSVM: 'ResponseName' must"], ...
                             [" be a character vector."]));
            endif

          case "classnames"
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames)
                                          || islogical (ClassNames)))
              error (strcat (["ClassificationSVM: 'ClassNames' must be a"], ...
                             [" cellstring, logical or numeric vector."]));
            endif
            ## Check that all class names are available in gnY
            if (iscellstr (ClassNames))
              if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                   ClassNames, "UniformOutput", false))))
                error (strcat (["ClassificationSVM: not all 'ClassNames'"], ...
                               [" are present in Y."]));
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, "UniformOutput", false))))
                error (strcat (["ClassificationSVM: not all 'ClassNames'"], ...
                               [" are present in Y."]));
              endif
            endif

          case "prior"
            Prior = varargin{2};
            if (! ((isnumeric (Prior) && isvector (Prior)) ||
                  (strcmpi (Prior, "empirical") || strcmpi (Prior, "uniform"))))
              error (strcat (["ClassificationSVM: 'Prior' must be either"], ...
                             [" a numeric vector or a character vector."]));
            endif

          case "cost"
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error (strcat (["ClassificationSVM: 'Cost' must be"], ...
                             [" a numeric square matrix."]));
            endif

          case "scoretransform"
            name = "ClassificationSVM";
            this.ScoreTransform = parseScoreTransform (varargin{2}, name);

          case "svmtype"
            SVMtype = varargin{2};
            SVMtype_override = false;
            if (! any (strcmp (SVMtype, {"c_svc", "nu_svc", "one_class_svm"})))
              error (strcat (["ClassificationSVM: 'SVMtype' must be"], ...
                             [" 'c_svc', 'nu_svc', or 'one_class_svm'."]));
            endif

          case "outlierfraction"
            Nu = varargin{2};
            if (! (isscalar (Nu) && Nu >= 0 && Nu < 1))
              error (strcat (["ClassificationSVM: 'OutlierFraction' must"], ...
                             [" be a positive scalar in the range 0 =<"], ...
                             [" OutlierFraction < 1."]));
            endif
            if (Nu > 0)
              SVMtype = 'nu_svc';
            endif

          case "kernelfunction"
            KernelFunction = varargin{2};
            if (! ischar (KernelFunction))
              error (strcat (["ClassificationSVM: 'KernelFunction' must"], ...
                             [" be a character vector."]));
            endif
            KernelFunction = tolower (KernelFunction);
            if (! any (strcmpi (KernelFunction, ...
                       {"linear", "rbf", "gaussian", "polynomial", "sigmoid"})))
              error ("ClassificationSVM: unsupported Kernel function.");
            endif

          case "polynomialorder"
            PolynomialOrder = varargin{2};
            if (! (isnumeric (PolynomialOrder) && isscalar (PolynomialOrder)
                   && PolynomialOrder > 0 && mod (PolynomialOrder, 1) == 0))
              error (strcat (["ClassificationSVM: 'PolynomialOrder' must"], ...
                             [" be a positive integer."]));
            endif

          case "kernelscale"
            KernelScale = varargin{2};
            if (! (isscalar (KernelScale) && KernelScale > 0))
              error (strcat (["ClassificationSVM: 'KernelScale'"], ...
                             [" must be a positive scalar."]));
            endif

          case "kerneloffset"
            KernelOffset = varargin{2};
            if (! (isnumeric (KernelOffset) && isscalar (KernelOffset)
                                            && KernelOffset >= 0))
              error (strcat (["ClassificationSVM: 'KernelOffset' must"], ...
                             [" be a non-negative scalar."]));
            endif

          case "boxconstraint"
            BoxConstraint = varargin{2};
            if (! (isscalar (BoxConstraint) && BoxConstraint > 0))
              error (strcat (["ClassificationSVM: 'BoxConstraint' must"], ...
                             [" be a positive scalar."]));
            endif

          case "nu"
            Nu = varargin{2};
            if (SVMtype_override)
              SVMtype = 'one_class_svm';
            endif
            if (! (isscalar (Nu) && Nu > 0 && Nu <= 1))
              error (strcat (["ClassificationSVM: 'Nu' must be a positive"], ...
                             [" scalar in the range 0 < Nu <= 1."]));
            endif

          case "cachesize"
            CacheSize = varargin{2};
            if (! (isscalar (CacheSize) && CacheSize > 0))
              error (strcat (["ClassificationSVM: 'CacheSize' must"], ...
                             [" be a positive scalar."]));
            endif

          case "tolerance"
            Tolerance = varargin{2};
            if (! (isscalar (Tolerance) && Tolerance >= 0))
              error (strcat (["ClassificationSVM: 'Tolerance' must"], ...
                             [" be a positive scalar."]));
            endif

          case "shrinking"
            Shrinking = varargin{2};
            if (! (ismember (Shrinking, [0, 1]) && isscalar (Shrinking)))
              error ("ClassificationSVM: 'Shrinking' must be either 0 or 1.");
            endif

          otherwise
            error (strcat (["ClassificationSVM: invalid parameter name"], ...
                           [" in optional pair arguments."]));

        endswitch
        varargin (1:2) = [];
      endwhile

      ## Get number of variables in training data
      ndims_X = columns (X);

      ## Assign the number of predictors to the ClassificationKNN object
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

      ## Remove missing values from X and Y
      RowsUsed  = ! logical (sum (isnan ([X, gY]), 2));
      Y         = Y (RowsUsed);
      X         = X (RowsUsed, :);

      ## Renew groups in Y
      [gY, gnY, glY] = grp2idx (Y);
      nclasses = numel (gnY);
      this.ClassNames = glY;  # Keep the same type as Y

      ## If only one class available, force 'SVMtype' to 'one_class_svm'
      if (nclasses == 1)
        if (! SVMtype_override && ! strcmp (SVMtype, 'one_class_svm'))
          error (strcat (["ClassificationSVM: cannot train a binary"], ...
                         [" problem with only one class available."]));
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
        error (strcat (["ClassificationSVM: can only be used for"], ...
                       [" one-class or two-class learning."]));
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

      ## Assign the number of observations and their correspoding indices
      ## on the original data, which will be used for training the model,
      ## to the ClassificationSVM object
      this.NumObservations = rows (X);
      this.RowsUsed = cast (RowsUsed, "double");

      ## Handle Standardize flag
      if (Standardize)
        this.Standardize = true;
        this.Sigma = std (X, [], 1);
        this.Sigma(this.Sigma == 0) = 1;  # predictor is constant
        this.Mu = mean (X, 1);
      else
        this.Standardize = false;
        this.Sigma = [];
        this.Mu = [];
      endif

      ## Handle Prior and Cost
      if (strcmpi ("uniform", Prior))
        this.Prior = ones (size (gnY)) ./ nclasses;
      elseif (isempty (Prior) || strcmpi ("empirical", Prior))
        pr = [];
        for i = 1:nclasses
          pr = [pr; sum(gY==i)];
        endfor
        this.Prior = pr ./ sum (pr);
      elseif (isnumeric (Prior))
        if (nclasses != numel (Prior))
          error (strcat (["ClassificationSVM: the elements in 'Prior' do"], ...
                         [" not correspond to selected classes in Y."]));
        endif
        this.Prior = Prior ./ sum (Prior);
      endif
      if (isempty (Cost))
        this.Cost = cast (! eye (nclasses), "double");
      else
        if (nclasses != sqrt (numel (Cost)))
          error (strcat (["ClassificationSVM: the number of rows and"], ...
                         [" columns in 'Cost' must correspond to"], ...
                         [" the selected classes in Y."]));
        endif
        this.Cost = Cost;
      endif

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

      ## Set svmtrain parameters for SVMtype and KernelFunction
      switch (SVMtype)
        case "c_svc"
          s = 0;
        case "nu_svc"
          s = 1;
        case "one_class_svm"
          s = 2;
      endswitch
      switch (KernelFunction)
        case "linear"
          t = 0;
        case "polynomial"
          t = 1;
        case {"rbf", "gaussian"}
          t = 2;
        case "sigmoid"
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
      str_options = strcat (["-s %d -t %d -g %f -d %d -r %f"], ...
                            [" -c %f -n %f -m %f -e %e -h %d -q"]);
      svm_options = sprintf (str_options, s, t, g, PolynomialOrder, ...
                             KernelOffset, BoxConstraint, Nu, ...
                             CacheSize, Tolerance, Shrinking);

      ## Train the SVM model using svmtrain from libsvm
      Model = svmtrain (Y, X, svm_options);
      this.Model = Model;

      ## Populate ClassificationSVM object properties
      if (t == 0)   # linear kernel
        this.Alpha = Model.sv_coef;
      else          # other kernels
        this.Beta = Model.sv_coef;
      endif
      this.Bias = Model.rho;

      this.IsSupportVector = zeros (this.NumObservations, 1);
      this.IsSupportVector(Model.sv_indices) = 1;
      this.SupportVectorLabels = zeros (size (Model.sv_indices));
      ## Handle one class
      if (isempty (Model.nSV))
        this.SupportVectorLabels(Model.sv_indices) = -1;
      else
        idx = Model.nSV(1);
        this.SupportVectorLabels(Model.sv_indices([1:idx])) = -1;
        this.SupportVectorLabels(Model.sv_indices([idx+1:end])) = 1;
      endif
      this.SupportVectors = Model.SVs;

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
    ## @deftypefn  {ClassificationSVM} {@var{labels} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationSVM} {[@var{labels}, @var{scores}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the Support Vector Machine
    ## classification object.
    ##
    ## @code{@var{labels} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## trained Support Vector Machine classification model, @var{obj}.
    ## For one-class SVM model, +1 or -1 is returned.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationSVM} class object.
    ## @item
    ## @var{XC} must be an @math{MxP} numeric matrix with the same number of
    ## predictors @math{P} as the corresponding predictors of the SVM model in
    ## @var{obj}.
    ## @end itemize
    ##
    ## @code{[@var{labels}, @var{scores}] = predict (@var{obj}, @var{XC}} also
    ## returns @var{scores}, which contains the desicion values for each each
    ## prediction.   Alternatively, @var{scores} can contain the posterior
    ## probabilities if the ScoreTransform has been previously set using the
    ## @code{fitPosterior} method.
    ##
    ## @seealso{fitcsvm, ClassificationSVM.fitPosterior}
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
        error (strcat (["ClassificationSVM.predict:"], ...
                       [" XC must have the same number of"], ...
                       [" predictors as the trained model."]));
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

      if (nargout > 1)
        ## Apply ScoreTransform to return probability estimates
        if (! strcmp (this.ScoreTransform, "none"))
          f = this.ScoreTransform;
          if (! strcmp (class (f), "function_handle"))
            error (strcat (["ClassificationSVM.predict: 'ScoreTransform'"], ...
                           [" must be a 'function_handle' object."]));
          endif
          scores = f (scores);
        endif
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{labels} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationSVM} {[@var{labels}, @var{score}] =} resubPredict (@var{obj})
    ##
    ## Classify the training data using the trained Support Vector Machine
    ## classification object.
    ##
    ## @code{@var{labels} = resubPredict (@var{obj})} returns the vector of
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
    ## @code{[@var{labels}, @var{scores}] = resubPredict (@var{obj}} also
    ## returns @var{scores}, which contains the desicion values for each each
    ## prediction.   Alternatively, @var{scores} can contain the posterior
    ## probabilities if the ScoreTransform has been previously set using the
    ## @code{fitPosterior} method.
    ##
    ## @seealso{fitcsvm, ClassificationSVM.fitPosterior}
    ## @end deftypefn

    function [labels, scores] = resubPredict (this)

      ## Get used rows (if necessary)
      if (sum (this.RowsUsed) != rows (this.X))
        RowsUsed = logical (this.RowsUsed);
        X = this.X(RowsUsed);
        Y = this.Y(RowsUsed);
      else
        X = this.X;
        Y = this.Y;
      endif

      ## Standardize (if necessary)
      if (this.Standardize)
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

      if (nargout > 1)
        ## Apply ScoreTransform to return probability estimates
        if (! strcmp (this.ScoreTransform, "none"))
          f = this.ScoreTransform;
          if (! strcmp (class (f), "function_handle"))
            error (strcat (["ClassificationSVM.resubPredict: 'Score"], ...
                           ["Transform' must be a 'function_handle' object."]));
          endif
          scores = f (scores);
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

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationSVM.margin: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat (["ClassificationSVM.margin: X must have the same"], ...
                       [" number of predictors as the trained model."]));
      endif

      ## Check for valid Y
      if (isempty (Y))
        error ("ClassificationSVM.margin: Y is empty.");
      elseif (rows (X) != rows (Y))
        error (strcat (["ClassificationSVM.margin: Y must have"], ...
                       [" the same number of rows as X."]));
      endif

      [~, ~, dec_values_L] = svmpredict (Y, X, this.Model, '-q');
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

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationSVM.loss: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat (["ClassificationSVM.loss: X must have the same"], ...
                       [" number of predictors as the trained model."]));
      endif

      ## Check for valid Y
      if (isempty (Y))
        error ("ClassificationSVM.loss: Y is empty.");
      elseif (rows (X)!= rows (Y))
        error (strcat (["ClassificationSVM.loss: Y must have the same"], ...
                       [" number of rows as X."]));
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
              error (strcat (["ClassificationSVM.loss: 'LossFun'"], ...
                             [" must be a character vector."]));
            endif
            LossFun = tolower (LossFun);
            if (! any (strcmpi (LossFun, {"binodeviance", "classiferror", ...
                                          "exponential", "hinge", "logit", ...
                                          "quadratic"})))
              error ("ClassificationSVM.loss: unsupported Loss function.");
            endif

          case "weights"
            Weights = varargin{2};
            ## Validate if weights is a numeric vector
            if(! (isnumeric (Weights) && isvector (Weights)))
              error (strcat (["ClassificationSVM.loss: 'Weights'"], ...
                             [" must be a numeric vector."]));
            endif

            ## Check if the size of weights matches the number of rows in X
            if (numel (Weights) != size (X, 1))
              error (strcat (["ClassificationSVM.loss: size of 'Weights'"], ...
                            [" must be equal to the number of rows in X."]));
            endif

          otherwise
            error (strcat (["ClassificationSVM.loss: invalid parameter"], ...
                           [" name in optional pair arguments."]));
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
            error ("ClassificationSVM.loss: unsupported Loss function.");
        endswitch

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

    function L = resubLoss (this, varargin)

      if (mod(nargin, 2) != 1)
        error (strcat (["ClassificationSVM.resubLoss: Name-Value"], ...
                       [" arguments must be in pairs."]));
      endif

      ## Set default values before parsing optional parameters
      LossFun = 'classiferror';
      Weights = ones (size (this.X, 1), 1);

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))

          case "lossfun"
            LossFun = varargin{2};
            if (! ischar (LossFun))
              error (strcat (["ClassificationSVM.resubLoss: 'LossFun'"], ...
                             [" must be a character vector."]));
            endif
            LossFun = tolower (LossFun);
            if (! any (strcmpi (LossFun, {"binodeviance", "classiferror", ...
                                          "exponential", "hinge", "logit", ...
                                          "quadratic"})))
              error (strcat (["ClassificationSVM.resubLoss: unsupported"], ...
                             [" Loss function."]));
            endif

          case "weights"
            Weights = varargin{2};
            ## Validate if weights is a numeric vector
            if (! (isnumeric (Weights) && isvector (Weights)))
              error (strcat (["ClassificationSVM.resubLoss: 'Weights'"], ...
                             [" must be a numeric vector."]));
            endif

            ## Check if the size of weights matches the number of rows in X
            if (numel (Weights) != size (this.X, 1))
              error (strcat (["ClassificationSVM.resubLoss: size"], ...
                             [" of 'Weights' must be equal to the"], ...
                             [" number of rows in X."]));
            endif

          otherwise
            error (strcat (["ClassificationSVM.resubLoss: invalid"], ...
                           [" parameter name in optional pair arguments."]));
        endswitch
        varargin(1:2) = [];
      endwhile

      ## Compute the classification score
      [~, ~, dec_values_L] = svmpredict (this.Y, this.X, this.Model, '-q');

      ## Compute the margin
      margin = this.Y .* dec_values_L;

        ## Compute the loss based on the specified loss function
        switch tolower(LossFun)
          case "classiferror"
            L = mean ((margin <= 0) .* Weights);

          case "hinge"
            L = mean (max (0, 1 - margin) .* Weights);

          case "logit"
            L = mean (log (1 + exp (-margin)) .* Weights);

          case "exponential"
            L = mean (exp (-margin) .* Weights);

          case "quadratic"
            L = mean ((1 - margin).^2 .* Weights);

          case "binodeviance"
            L = mean (log (1 + exp (-2 * margin)) .* Weights);

          otherwise
            error ("ClassificationSVM.resubloss: unsupported Loss function.");
        endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationSVM} {@var{CVMdl} =} crossval (@dots{}, @var{Name}, @var{Value})
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
    ## @seealso{fitcsvm, ClassificationSVM, cvpartition,
    ## ClassificationPartitionedModel}
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
                       [" the optional Name-Value paired arguments."]));
      endif

      ## Add default values
      numFolds    = 10;
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
              error (strcat (["ClassificationSVM.crossval: 'KFold' must"], ...
                             [" be an integer value greater than 1."]));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat (["ClassificationSVM.crossval: 'Holdout' must"], ...
                             [" be a numeric value between 0 and 1."]));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat (["ClassificationSVM.crossval: 'Leaveout'"], ...
                             [" must be either 'on' or 'off'."]));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (!(isa (CVPartition, 'cvpartition')))
              error (strcat (["ClassificationSVM.crossval: 'CVPartition'"],...
                             [" must be a 'cvpartition' object."]));
            endif

          otherwise
            error (strcat (["ClassificationSVM.crossval: invalid"],...
                           [" parameter name in optional paired arguments."]));
          endswitch
        varargin (1:2) = [];
      endwhile

      ## Determine the cross-validation method to use
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (this.Y, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (this.Y, 'LeaveOut');
      else
        partition = cvpartition (this.Y, 'KFold', numFolds);
      endif

      ## Create a cross-validated model object
      CVMdl = ClassificationPartitionedModel (this, partition);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{Mdl} =} fitPosterior (@var{obj})
    ## @deftypefnx {ClassificationSVM} {@var{CVMdl} =} fitPosterior (@var{obj}, @var{name}, @var{value})
    ##
    ## Fit posterior probabilities to a Support Vector Machine model.
    ##
    ## @code{@var{Mdl} = fitPosterior (@var{obj})} returns the ClassificationSVM
    ## object, @var{Mdl}, from an already trained SVM model, @var{obj}, after
    ## fitting a posterior probabilities ScoreTransform.
    ##
    ## @code{@var{CVMdl} = fitPosterior (@var{obj}, @var{name}, @var{value})}
    ## returns the ClassificationPartitionedModel, @var{CVMdl}, from an already
    ## trained SVM model, @var{obj}, after fitting a posterior probabilities
    ## ScoreTransform. Use the additional name-value pair arguments to customize
    ## the cross-validation process.
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

    function CVMdl = fitPosterior (this, varargin)

      ## Check for ScoreTransform and emit reset warning


      ## Cross-validate SVM model and get labels and scores
      CVMdl = crossval (this, varargin{:});
      [~, score] = kfoldPredict (CVMdl);

      ## Get class labels at 0 and 1
      Y = grp2idx (CVMdl.Y(logical (this.RowsUsed))) - 1;

      ## Get prior probability for second class
      prior = this.Prior(2);

      ## Determine perfect separation or overlapping
      ub = max (score(Y==0));
      lb = min (score(Y==1));

      if (ub <= lb)
        warning ("ClassificationSVM.fitPosterior: PerfectSeparation.");
        f = eval (sprintf ('@(S) ClassificationSVM.step (S, %e, %e, %e)', ...
                           ub, lb, prior));
      else
        coeff = glmfit (score(:,2), Y, 'binomial', 'link', 'logit');
        f = eval (sprintf ('@(S) ClassificationSVM.sigmoid (S, %e, %e)', ...
                           -coeff(2), -coeff(1)));
      endif

      ## Decide returning model type
      if (isempty (varargin))
        this.ScoreTransform = f;
        CVMdl = this;
      else
        CVMdl.ScoreTransform = f;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {@var{CVMdl} =} compact (@var{obj})
    ##
    ## Create a CompactClassificationSVM object.
    ##
    ## @code{@var{CVMdl} = compact (@var{obj})} creates a compact version of the
    ## ClassificationSVM object, @var{obj}.
    ##
    ## @seealso{fitcnet, ClassificationSVM, CompactClassificationSVM}
    ## @end deftypefn

    function CVMdl = compact (this)
      ## Greate a compact model
      CVMdl = CompactClassificationSVM (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationSVM} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationSVM object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves a ClassificationSVM
    ## object into a file defined by @var{filename}.
    ##
    ## @seealso{loadmodel, fitcsvm, ClassificationSVM}
    ## @end deftypefn

    function savemodel (this, fname)
      ## Generate variable for class name
      classdef_name = "ClassificationSVM";

      ## Create variables from model properties
      X = this.X;
      Y = this.Y;
      NumObservations     = this.NumObservations;
      RowsUsed            = this.RowsUsed;
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
      save (fname, "classdef_name", "X", "Y", "NumObservations", "RowsUsed", ...
            "NumPredictors", "PredictorNames", "ResponseName", "ClassNames", ...
            "Prior", "Cost", "ScoreTransform", "Standardize", "Sigma", "Mu", ...
            "ModelParameters", "Model", "Alpha", "Beta", "Bias", ...
            "IsSupportVector", "SupportVectorLabels", "SupportVectors");
    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationSVM object
      mdl = ClassificationSVM (1, 1);

      ## Check that fieldnames in DATA match properties in ClassificationSVM
      names = fieldnames (data);
      props = fieldnames (mdl);
      if (! isequal (sort (names), sort (props)))
        error ("ClassificationSVM.load_model: invalid model in '%s'.", filename)
      endif

      ## Copy data into object
      for i = 1:numel (props)
        mdl.(props{i}) = data.(props{i});
      endfor
    endfunction

    ## Helper functions for fitPosterior
    function prob = step (score, ub, lb, prior)
      prob = zeros (size (score));
      prob(score > lb) = 1;
      prob(score >= ub & score <= lb) = prior;
    endfunction

    function prob = sigmoid (score, a, b)
      prob = zeros (size (score));
      prob = 1 ./ (1 + exp (-a * score + b));
    endfunction

  endmethods

endclassdef

%!demo
%! ## Create a Support Vector Machine classifier and determine margin for test
%! ## data.
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

## Test output of constructor
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1; 4, 5, 6; 7, 8, 9; ...
%! 3, 2, 1; 4, 5, 6; 7, 8, 9; 3, 2, 1; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [1; 2; 3; 4; 2; 3; 4; 2; 3; 4; 2; 3; 4];
%! a = ClassificationSVM (x, y, "ClassNames", [1, 2]);
%! assert (class (a), "ClassificationSVM");
%! assert (a.RowsUsed, [1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0]');
%! assert ({a.X, a.Y}, {x, y})
%! assert (a.NumObservations, 5)
%! assert ({a.ResponseName, a.PredictorNames}, {"Y", {"x1", "x2", "x3"}})
%! assert ({a.ClassNames, a.ModelParameters.SVMtype}, {[1; 2], "c_svc"})
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! a = ClassificationSVM (x, y);
%! assert (class (a), "ClassificationSVM");
%! assert ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, "linear"})
%! assert (a.ModelParameters.BoxConstraint, 1)
%! assert (a.ClassNames, [1; -1])
%! assert (a.ModelParameters.KernelOffset, 0)
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! a = ClassificationSVM (x, y, "KernelFunction", "rbf", "BoxConstraint", 2, ...
%! "KernelOffset", 2);
%! assert (class (a), "ClassificationSVM");
%! assert ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, "rbf"})
%! assert (a.ModelParameters.BoxConstraint, 2)
%! assert (a.ModelParameters.KernelOffset, 2)
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! a = ClassificationSVM (x, y, "KernelFunction", "polynomial", ...
%! "PolynomialOrder", 3);
%! assert (class (a), "ClassificationSVM");
%! assert ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, "polynomial"})
%! assert (a.ModelParameters.PolynomialOrder, 3)

## Test input validation for constructor
%!error<ClassificationSVM: too few input arguments.> ClassificationSVM ()
%!error<ClassificationSVM: too few input arguments.> ...
%! ClassificationSVM (ones(10,2))
%!error<ClassificationSVM: number of rows in X and Y must be equal.> ...
%! ClassificationSVM (ones(10,2), ones (5,1))
%!error<ClassificationSVM: 'Standardize' must be either true or false.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "Standardize", 'a')
%!error<ClassificationSVM: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "PredictorNames", ['x1';'x2'])
%!error<ClassificationSVM: 'PredictorNames' must have the same number of columns as X.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "PredictorNames", {'x1','x2','x3'})
%!error<ClassificationSVM: 'ResponseName' must be a character vector.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "ResponseName", {'Y'})
%!error<ClassificationSVM: 'ResponseName' must be a character vector.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "ResponseName", 21)
%!error<ClassificationSVM: 'ClassNames' must be a cellstring, logical or numeric vector.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "ClassNames", @(x)x)
%!error<ClassificationSVM: 'ClassNames' must be a cellstring, logical or numeric vector.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "ClassNames", ['a'])
%!error<ClassificationSVM: not all 'ClassNames' are present in Y.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "ClassNames", [1, 2])
%!error<ClassificationSVM: not all 'ClassNames' are present in Y.> ...
%! ClassificationSVM (ones(5,2), {'a';'b';'a';'a';'b'}, "ClassNames", {'a','c'})
%!error<ClassificationSVM: not all 'ClassNames' are present in Y.> ...
%! ClassificationSVM (ones(10,2), logical (ones (10,1)), "ClassNames", [true, false])
%!error<ClassificationSVM: 'Prior' must be either a numeric vector or a character vector.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "Prior", {"asd"})
%!error<ClassificationSVM: 'Prior' must be either a numeric vector or a character vector.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "Prior", ones (2))
%!error<ClassificationSVM: 'Cost' must be a numeric square matrix.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "Cost", [1:4])
%!error<ClassificationSVM: 'Cost' must be a numeric square matrix.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "Cost", {0,1;1,0})
%!error<ClassificationSVM: 'Cost' must be a numeric square matrix.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "Cost", 'a')
%!error<ClassificationSVM: 'SVMtype' must be 'c_svc', 'nu_svc', or 'one_class_svm'.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "svmtype", 123)
%!error<ClassificationSVM: 'SVMtype' must be 'c_svc', 'nu_svc', or 'one_class_svm'.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "svmtype", 'some_type')
%!error<ClassificationSVM: 'OutlierFraction' must be a positive scalar in the range 0 =< OutlierFraction < 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "OutlierFraction", -1)
%!error<ClassificationSVM: 'KernelFunction' must be a character vector.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "KernelFunction", 123)
%!error<ClassificationSVM: unsupported Kernel function.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "KernelFunction", "fcn")
%!error<ClassificationSVM: 'PolynomialOrder' must be a positive integer.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "PolynomialOrder", -1)
%!error<ClassificationSVM: 'PolynomialOrder' must be a positive integer.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "PolynomialOrder", 0.5)
%!error<ClassificationSVM: 'PolynomialOrder' must be a positive integer.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "PolynomialOrder", [1,2])
%!error<ClassificationSVM: 'KernelScale' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "KernelScale", -1)
%!error<ClassificationSVM: 'KernelScale' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "KernelScale", 0)
%!error<ClassificationSVM: 'KernelScale' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "KernelScale", [1, 2])
%!error<ClassificationSVM: 'KernelScale' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "KernelScale", "invalid")
%!error<ClassificationSVM: 'KernelOffset' must be a non-negative scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "KernelOffset", -1)
%!error<ClassificationSVM: 'KernelOffset' must be a non-negative scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "KernelOffset", [1,2])
%!error<ClassificationSVM: 'BoxConstraint' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "BoxConstraint", -1)
%!error<ClassificationSVM: 'BoxConstraint' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "BoxConstraint", 0)
%!error<ClassificationSVM: 'BoxConstraint' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "BoxConstraint", [1, 2])
%!error<ClassificationSVM: 'BoxConstraint' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "BoxConstraint", "invalid")
%!error<ClassificationSVM: 'Nu' must be a positive scalar in the range 0 < Nu <= 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "nu", -0.5)
%!error<ClassificationSVM: 'Nu' must be a positive scalar in the range 0 < Nu <= 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "nu", 0)
%!error<ClassificationSVM: 'Nu' must be a positive scalar in the range 0 < Nu <= 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "nu", 1.5)
%!error<ClassificationSVM: 'CacheSize' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "CacheSize", -1)
%!error<ClassificationSVM: 'CacheSize' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "CacheSize", [1,2])
%!error<ClassificationSVM: 'Tolerance' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "Tolerance", -0.1)
%!error<ClassificationSVM: 'Tolerance' must be a positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "Tolerance", [0.1,0.2])
%!error<ClassificationSVM: 'Shrinking' must be either 0 or 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "shrinking", 2)
%!error<ClassificationSVM: 'Shrinking' must be either 0 or 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "shrinking", -1)
%!error<ClassificationSVM: 'Shrinking' must be either 0 or 1.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "shrinking", [1 0])
%!error<ClassificationSVM: invalid parameter name in optional pair arguments.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "invalid_name", 'c_svc')
%!error<ClassificationSVM: cannot train a binary problem with only one class available.> ...
%! ClassificationSVM (ones(10,2), ones(10,1), "SVMtype", 'c_svc')
%!error<ClassificationSVM: can only be used for one-class or two-class learning.> ...
%! ClassificationSVM (ones(10,2), [1;1;1;1;2;2;2;2;3;3])
%!error<ClassificationSVM: invalid values in X.> ...
%! ClassificationSVM ([ones(9,2);2,Inf], ones(10,1))
%!error<ClassificationSVM: the elements in 'Prior' do not correspond to selected classes in Y.> ...
%! ClassificationSVM (ones (5,2), ones (5,1), "Prior", [0,1])
%!error<ClassificationSVM: the elements in 'Prior' do not correspond to selected classes in Y.> ...
%! ClassificationSVM (ones (5,2), [1;1;2;2;3], "ClassNames", [1,2], "Prior", [0,0.4,0.6])
%!error<ClassificationSVM: the number of rows and columns in 'Cost' must correspond to the selected classes in Y.> ...
%! ClassificationSVM (ones (5,2), [1;1;2;2;3], "ClassNames", [1,2], "Cost", ones (3))

## Test output for predict method
%!shared x, y, x_train, x_test, y_train, y_test, objST
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%!test
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcsvm (x, y, 'KernelFunction', 'rbf', 'Tolerance', 1e-7);
%! assert (isempty (obj.Alpha), true)
%! assert (sum (obj.IsSupportVector), numel (obj.Beta))
%! [label, score] = predict (obj, xc);
%! assert (label, [1; 2; 2]);
%! assert (score(:,1), [0.99285; -0.080296; -0.93694], 1e-5);
%! assert (score(:,1), -score(:,2), eps)
%! obj = fitPosterior (obj);
%! [label, probs] = predict (obj, xc);
%! assert (probs(:,2), [0.97555; 0.428164; 0.030385], 1e-5);
%! assert (probs(:,1) + probs(:,2), [1; 1; 1], 0.05)
%!test
%! obj = fitcsvm (x, y);
%! assert (isempty (obj.Beta), true)
%! assert (sum (obj.IsSupportVector), numel (obj.Alpha))
%! assert (numel (obj.Alpha), 24)
%! assert (obj.Bias, -14.415, 1e-3)
%! xc = [min(x); mean(x); max(x)];
%! label = predict (obj, xc);
%! assert (label, [1; 2; 2]);

## Test input validation for predict method
%!error<ClassificationSVM.predict: too few input arguments.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)))
%!error<ClassificationSVM.predict: XC is empty.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), [])
%!error<ClassificationSVM.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (ClassificationSVM (ones (40,2), ones (40,1)), 1)
%!test
%! objST = fitcsvm (x, y);
%! objST.ScoreTransform = "a";
%!error<ClassificationSVM.predict: 'ScoreTransform' must be a 'function_handle' object.> ...
%! [labels, scores] = predict (objST, x);

## Test input validation for resubPredict method
%!error<ClassificationSVM.resubPredict: 'ScoreTransform' must be a 'function_handle' object.> ...
%! [labels, scores] = resubPredict (objST);

## Test output for margin method
%!test
%! rand ("seed", 1);
%! CVSVMModel = fitcsvm (x, y, 'KernelFunction', 'rbf', 'HoldOut', 0.15, ...
%!                       'Tolerance', 1e-7);
%! obj = CVSVMModel.Trained{1};
%! testInds = test (CVSVMModel.Partition);
%! expected_margin = [2.0000;  0.8579;  1.6690;  3.4141;  3.4552; ...
%!                    2.6605;  3.5251; -4.0000; -6.3411; -6.4511; ...
%!                   -3.0532; -7.5054; -1.6700; -5.6227; -7.3640];
%! computed_margin = margin (obj, x(testInds,:), y(testInds,:));
%! assert (computed_margin, expected_margin, 1e-4);

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
%! rand ("seed", 1);
%! CVSVMModel = fitcsvm (x, y, 'KernelFunction', 'rbf', 'HoldOut', 0.15);
%! obj = CVSVMModel.Trained{1};
%! testInds = test (CVSVMModel.Partition);
%! L1 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'binodeviance');
%! L2 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'classiferror');
%! L3 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'exponential');
%! L4 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'hinge');
%! L5 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'logit');
%! L6 = loss (obj, x(testInds,:), y(testInds,:), 'LossFun', 'quadratic');
%! assert (L1, 2.8711, 1e-4);
%! assert (L2, 0.5333, 1e-4);
%! assert (L3, 10.9685, 1e-4);
%! assert (L4, 1.9827, 1e-4);
%! assert (L5, 1.5849, 1e-4);
%! assert (L6, 7.6739, 1e-4);

## Test input validation for loss method
%!error<ClassificationSVM.loss: too few input arguments.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)))
%!error<ClassificationSVM.loss: too few input arguments.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2))
%!error<ClassificationSVM.loss: Name-Value arguments must be in pairs.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones(2,1), "LossFun")
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
%! ones (2,1), "LossFun", 1)
%!error<ClassificationSVM.loss: unsupported Loss function.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), "LossFun", "some")
%!error<ClassificationSVM.loss: 'Weights' must be a numeric vector.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), "Weights", ['a','b'])
%!error<ClassificationSVM.loss: 'Weights' must be a numeric vector.> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), "Weights", 'a')
%!error<ClassificationSVM.loss: size of 'Weights' must be equal to the number> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), "Weights", [1,2,3])
%!error<ClassificationSVM.loss: size of 'Weights' must be equal to the number> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), "Weights", 3)
%!error<ClassificationSVM.loss: invalid parameter name in optional pair arg> ...
%! loss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), zeros (2), ...
%! ones (2,1), "some", "some")

## Test input validation for resubLoss method
%!error<ClassificationSVM.resubLoss: Name-Value arguments must be in pairs.> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), "LossFun")
%!error<ClassificationSVM.resubLoss: 'LossFun' must be a character vector.> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), "LossFun", 1)
%!error<ClassificationSVM.resubLoss: unsupported Loss function.> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), "LossFun", "some")
%!error<ClassificationSVM.resubLoss: 'Weights' must be a numeric vector.> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), "Weights", ['a','b'])
%!error<ClassificationSVM.resubLoss: 'Weights' must be a numeric vector.> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), "Weights", 'a')
%!error<ClassificationSVM.resubLoss: size of 'Weights' must be equal to the n> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), "Weights", [1,2,3])
%!error<ClassificationSVM.resubLoss: size of 'Weights' must be equal to the n> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), "Weights", 3)
%!error<ClassificationSVM.resubLoss: invalid parameter name in optional pai> ...
%! resubLoss (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1)), "some", "some")

## Test output for crossval method
%!test
%! SVMModel = fitcsvm (x,y);
%! CVMdl = crossval (SVMModel, "KFold", 5);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 5)
%! assert (class (CVMdl.Trained{1}), "ClassificationSVM")
%!test
%! obj = fitcsvm (x, y);
%! CVMdl = crossval (obj, "HoldOut", 0.2);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "ClassificationSVM")
%!test
%! obj = fitcsvm (x, y);
%! CVMdl = crossval (obj, "LeaveOut", 'on');
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "ClassificationSVM")

## Test input validation for crossval method
%!error<ClassificationSVM.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold")
%!error<ClassificationSVM.crossval: specify only one of the optional Name-Value paired arguments.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), ...
%! "KFold", 5, "leaveout", 'on')
%!error<ClassificationSVM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", 'a')
%!error<ClassificationSVM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", 1)
%!error<ClassificationSVM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", -1)
%!error<ClassificationSVM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", 11.5)
%!error<ClassificationSVM.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "KFold", [1,2])
%!error<ClassificationSVM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 'a')
%!error<ClassificationSVM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 11.5)
%!error<ClassificationSVM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Holdout", -1)
%!error<ClassificationSVM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 0)
%!error<ClassificationSVM.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Holdout", 1)
%!error<ClassificationSVM.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "Leaveout", 1)
%!error<ClassificationSVM.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "CVPartition", 1)
%!error<ClassificationSVM.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "CVPartition", 'a')
%!error<ClassificationSVM.crossval: invalid parameter name in optional paired arguments> ...
%! crossval (ClassificationSVM (ones (40,2),randi([1, 2], 40, 1)), "some", "some")
