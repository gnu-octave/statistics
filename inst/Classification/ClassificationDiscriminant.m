## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
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

classdef ClassificationDiscriminant
## -*- texinfo -*-
## @deftypefn  {statistics} {@var{obj} =} ClassificationDiscriminant (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} ClassificationDiscriminant (@dots{}, @var{name}, @var{value})
##
## Create a @qcode{ClassificationDiscriminant} class object containing a
## discriminant analysis model.
##
## @code{@var{obj} = ClassificationDiscriminant (@var{X}, @var{Y})} returns a
## ClassificationDiscriminant object, with @var{X} as the predictor data
## and @var{Y} containing the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{NxP} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the discriminant model.
## @item
## @code{Y} is @math{Nx1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}. @var{Y} can contain any type of
## categorical data. @var{Y} must have the same number of rows as @var{X}.
## @end itemize
##
## @code{@var{obj} = ClassificationDiscriminant (@dots{}, @var{name},
## @var{value})} returns a ClassificationDiscriminant object with parameters
## specified by @qcode{Name-Value} pair arguments.
## Type @code{help ClassificationDiscriminant} for more info.
##
## A @qcode{ClassificationDiscriminant} object, @var{obj}, stores the labeled
## training data and various parameters for the discriminant analysis model,
## which can be accessed in the following fields:
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Field} @tab @tab @var{Description}
##
## @item @qcode{obj.X} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Y} @tab @tab Class labels, specified as a logical,
## numeric vector, or cell array of character vectors. Each value in @var{Y}
## is the observed class label for the corresponding row in @var{X}.
##
## @item @qcode{obj.NumObservations} @tab @tab Number of observations used in
## training the ClassificationDiscriminant model, specified as a positive
## integer scalar.
##
## @item @qcode{obj.RowsUsed} @tab @tab Rows of the original training data
## used in fitting the ClassificationDiscriminant model, specified as a
## numerical vector.
##
## @item @qcode{obj.PredictorNames} @tab @tab Predictor variable names,
## specified as a cell array of character vectors. The variable names are in
## the same order in which they appear in the training data @var{X}.
##
## @item @qcode{obj.ResponseName} @tab @tab Response variable name, specified
## as a character vector.
##
## @item @qcode{obj.ClassNames} @tab @tab Names of the classes in the training
## data @var{Y} with duplicates removed, specified as a cell array of character
## vectors.
##
## @item @qcode{obj.Prior} @tab @tab Prior probabilities for each class,
## specified as a numeric vector.  The order of the elements in @qcode{Prior}
## corresponds to the order of the classes in @qcode{ClassNames}.
##
## @item @qcode{obj.Cost} @tab @tab Cost of the misclassification of a point,
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
## @item @qcode{obj.Sigma} @tab @tab Within-class covariance matrix, specified
## as a numeric matrix. For 'linear' discriminant type matrix is of size
## @math{pxp}, where p is the number of predictors.
##
## @item @qcode{obj.Mu} @tab @tab Class means, specified as a @math{Kxp}
## real matrix. K is the number of classes, and p is the number of
## predictors.
##
## @item @qcode{obj.Coeffs} @tab @tab Coefficient matrices, specified as a
## struct array.
##
## @item @qcode{obj.Delta} @tab @tab Threshold for linear discriminant model,
## specified as a numeric scalar.
##
## @item @qcode{obj.DiscrimType} @tab @tab Discriminant type, specified as a
## character vector.
##
## @item @qcode{obj.Gamma} @tab @tab Gamma regularization parameter, specified
## as a numeric scalar.
##
## @item @qcode{obj.MinGamma} @tab @tab Minimum value of Gamma so that the
## correlation matrix is invertible, specified as nonnegative scalar.
##
## @item @qcode{obj.LogDetSigma} @tab @tab Logarithm of the determinant of the
## within-class covariance matrix. For linear discriminant analysis it is
## specified as a numeric scalar.
##
## @item @qcode{obj.XCentered} @tab @tab X data with class means
## subtracted, returned as a real matrix.
##
## @end multitable
##
## @seealso{fitcdiscr}
## @end deftypefn

  properties (Access = public)
    X = [];                   # Predictor data
    Y = [];                   # Class labels

    NumObservations = [];     # Number of observations in training dataset
    RowsUsed        = [];     # Rows used in fitting
    NumPredictors   = [];     # Number of predictors
    PredictorNames  = [];     # Predictor variables names
    ResponseName    = [];     # Response variable name
    ClassNames      = [];     # Names of classes in Y
    Prior           = [];     # Prior probability for each class
    Cost            = [];     # Cost of Misclassification

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

  methods (Access = public)

    ## constructor
    function this = ClassificationDiscriminant (X, Y, varargin)

      ## Check for sufficient number of input arguments
      if (nargin < 2)
        error ("ClassificationDiscriminant: too few input arguments.");
      endif

      ## Validate X
      if (! isnumeric (X))
        error ("ClassificationDiscriminant: X must be a numeric matrix.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error (["ClassificationDiscriminant: number of rows ", ...
                "in X and Y must be equal."]);
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
      PredictorNames       = [];
      ResponseName         = 'Y';
      Prior                = "empirical";
      FillCoeffs           = "on";

      ## Parse optional parameters
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))

          case "predictornames"
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat (["ClassificationDiscriminant: 'PredictorNames'"], ...
                             [" must be supplied as a cellstring array."]));
            elseif (numel (PredictorNames) != columns (X))
              error (strcat (["ClassificationDiscriminant: 'PredictorNames'"], ...
                             [" must equal the number of columns in X."]));
            endif

          case "responsename"
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat (["ClassificationDiscriminant: 'ResponseName'"], ...
                             [" must be a character vector."]));
            endif

          case "classnames"
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames)
                                          || islogical (ClassNames)))
              error (strcat (["ClassificationDiscriminant: 'ClassNames'"], ...
                             [" must be a cellstring, logical or numeric"], ...
                             [" vector."]));
            endif
            ## Check that all class names are available in gnY
            msg = strcat (["ClassificationDiscriminant: not all"], ...
                          [" 'ClassNames' are present in Y."]);
            if (iscellstr (ClassNames))
              if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                   ClassNames, "UniformOutput", false))))
                error (msg);
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, "UniformOutput", false))))
                error (msg);
              endif
            endif

          case "prior"
            Prior = varargin{2};
            if (! ((isnumeric (Prior) && isvector (Prior)) ||
                  (strcmpi (Prior, "empirical") || strcmpi (Prior, "uniform"))))
              error (strcat (["ClassificationDiscriminant: 'Prior' must"], ...
                             [" be either a numeric or a character vector."]));
            endif

          case "cost"
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error (strcat (["ClassificationDiscriminant: 'Cost'"], ...
                             [" must be a numeric square matrix."]));
            endif

          case "discrimtype"
            DiscrimType = tolower (varargin{2});
            if (! (strcmpi (DiscrimType, "linear")))
                error (strcat (["ClassificationDiscriminant: unsupported "], ...
                               [" discriminant type."]));
            endif

          case "fillcoeffs"
            FillCoeffs = tolower (varargin{2});
            if (! any (strcmpi (FillCoeffs, {"on", "off"})))
              error (strcat (["ClassificationDiscriminant: 'FillCoeffs'"], ...
                             [" must be 'on' or 'off'."]));
            endif

          case "gamma"
            Gamma = varargin{2};
            if (Gamma >= 1 || Gamma < 0)
              error (strcat (["ClassificationDiscriminant: 'Gamma'"], ...
                             [" must be between 0 and 1."]));
            endif

          otherwise
            error (strcat (["ClassificationDiscriminant: invalid"], ...
                           [" parameter name in optional pair arguments."]));
        endswitch
        varargin (1:2) = [];
      endwhile

      ## Generate default predictors and response variabe names (if necessary)
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

      ## Renew groups in Y
      [gY, gnY, glY] = grp2idx (Y);
      this.ClassNames = gnY;

      ## Check X contains valid data
      if (! (isnumeric (X) && isfinite (X)))
        error ("ClassificationDiscriminant: invalid values in X.");
      endif

      this.NumObservations = rows (X);
      this.RowsUsed = cast (RowsUsed, "double");

      ## Handle Prior and Cost
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
          error (strcat (["ClassificationDiscriminant: the elements"], ...
                         [" in 'Prior' do not correspond to the"], ...
                         [" selected classes in Y."]));
        endif
        this.Prior = Prior ./ sum (Prior);
      endif
      if (isempty (Cost))
        this.Cost = cast (! eye (numel (gnY)), "double");
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (strcat (["ClassificationDiscriminant: the number"], ...
                         [" of rows and columns in 'Cost' must"], ...
                         [" correspond to selected classes in Y."]));
        endif
        this.Cost = Cost;
      endif

      ## Assign DiscrimType
      this.DiscrimType = DiscrimType;
      this.Delta = Delta;
      this.Gamma = Gamma;

      num_classes = numel (this.ClassNames);
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

        ## Check for predictors zero within-class variance
        zwcv = find (diag (this.Sigma) == 0);
        if (! isempty (zwcv))
          msg = strcat (["ClassificationDiscriminant: Predictor"], ...
                        [" '%s' has zero within-class variance."]);
          error (msg, PredictorNames{zwcv(1)});
        endif

        D = diag (diag (this.Sigma));

        ## fix me: MinGamma calculation is not same as Matlab.
        ## Instead of using (det (sigma) > 0) this criteria (line no: 391)
        ## may be Matlab is using some threshold.
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
              error (["ClassificationDiscriminant: failed to ", ...
                      "find 'MinGamma' within reasonable range."]);
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
                this.Coeffs(i, j).Class1 = this.ClassNames{i};
                this.Coeffs(i, j).Class2 = this.ClassNames{j};
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
        error (strcat (["ClassificationDiscriminant.predict: XC must have"], ...
                       [" the same number of features as the trained model."]));
      endif

      ## Get training data and labels
      X = this.X(logical (this.RowsUsed),:);
      Y = this.Y(logical (this.RowsUsed),:);

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
        error (strcat (["ClassificationDiscriminant.loss: name-value"], ...
                       [" arguments must be in pairs."]));
      elseif (nargin > 7)
        error ("ClassificationDiscriminant.loss: too many input arguments.");
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
        error (strcat (["ClassificationDiscriminant.loss: Y must"], ...
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
                error (strcat (["ClassificationDiscriminant.loss: custom"], ...
                               [" loss function must accept exactly four"], ...
                               [" input arguments."]));
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
                  error (strcat (["ClassificationDiscriminant.loss:"], ...
                                 [" custom loss function must return"], ...
                                 [" a scalar value."]));
                endif
              catch
                error (strcat (["ClassificationDiscriminant.loss: custom"], ...
                               [" loss function is not valid or does not"], ...
                               [" produce correct output."]));
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
                error (strcat (["ClassificationDiscriminant.loss: number"], ...
                               [" of 'Weights' must be equal to the"], ...
                               [" number of rows in X."]));
              elseif (numel (Value) == size (X, 1))
                Weights = Value;
              endif
            else
              error ("ClassificationDiscriminant.loss: invalid 'Weights'.");
            endif

          otherwise
            error (strcat (["ClassificationDiscriminant.loss: invalid"], ...
                           [" parameter name in optional pair arguments."]));
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
        error (strcat (["ClassificationDiscriminant.loss: Y must be a"], ...
                       [" numeric, logical, char, string, or cell array."]));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), this.ClassNames)))
        error (strcat (["ClassificationDiscriminant.loss: Y must contain"], ...
                       [" only the classes in ClassNames."]));
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
          error ("ClassificationDiscriminant.loss: invalid loss function.");
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationDiscriminant} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
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
        error (strcat (["ClassificationDiscriminant.margin: Y must"], ...
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
        error (strcat (["ClassificationDiscriminant.margin: Y must be"], ...
                       [" a numeric, logical, char, string, or cell array."]));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), this.ClassNames)))
        error (strcat (["ClassificationDiscriminant.margin: Y must"], ...
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
        error (strcat (["ClassificationDiscriminant.crossval: Name-Value"], ...
                       [" arguments must be in pairs."]));
      elseif (numel (varargin) > 2)
        error (strcat (["ClassificationDiscriminant.crossval: specify only"], ...
                       [" one of the optional Name-Value paired arguments."]));
      endif

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
            if (! (isnumeric (numFolds) && isscalar (numFolds)
                   && (numFolds == fix (numFolds)) && numFolds > 1))
              error (strcat (["ClassificationDiscriminant.crossval: 'KFold'"], ...
                             [" must be an integer value greater than 1."]));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat (["ClassificationDiscriminant.crossval: 'Holdout'"], ...
                             [" must be a numeric value between 0 and 1."]));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat (["ClassificationDiscriminant.crossval:"], ...
                             [" 'Leaveout' must be either 'on' or 'off'."]));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (!(isa (CVPartition, 'cvpartition')))
              error (strcat (["ClassificationDiscriminant.crossval:"],...
                             [" 'CVPartition' must be a 'cvpartition' object."]));
            endif

          otherwise
            error (strcat (["ClassificationDiscriminant.crossval: invalid"],...
                           [" parameter name in optional paired arguments."]));
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

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationDiscriminant object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves a
    ## ClassificationDiscriminant object into a file defined by @var{filename}.
    ##
    ## @seealso{loadmodel, fitcdiscr, ClassificationDiscriminant}
    ## @end deftypefn

    function savemodel (obj, fname)
      ## Generate variable for class name
      classdef_name = "ClassificationDiscriminant";

      ## Create variables from model properties
      X = obj.X;
      Y = obj.Y;
      NumObservations = obj.NumObservations;
      RowsUsed        = obj.RowsUsed;
      NumPredictors   = obj.NumPredictors;
      PredictorNames  = obj.PredictorNames;
      ResponseName    = obj.ResponseName;
      ClassNames      = obj.ClassNames;
      Prior           = obj.Prior;
      Cost            = obj.Cost;
      Sigma           = obj.Sigma;
      Mu              = obj.Mu;
      Coeffs          = obj.Coeffs;
      Delta           = obj.Delta;
      DiscrimType     = obj.DiscrimType;
      Gamma           = obj.Gamma;
      MinGamma        = obj.MinGamma;
      LogDetSigma     = obj.LogDetSigma;
      XCentered       = obj.XCentered;

      ## Save classdef name and all model properties as individual variables
      save (fname, "classdef_name", "X", "Y", "NumObservations", "RowsUsed", ...
            "NumPredictors", "PredictorNames", "ResponseName", "ClassNames", ...
            "Prior", "Cost", "Sigma", "Mu", "Coeffs", "Delta", ...
            "DiscrimType", "Gamma", "MinGamma", "LogDetSigma", "XCentered");
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

## Test Constructor
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! PredictorNames = {'Feature1', 'Feature2', 'Feature3'};
%! a = ClassificationDiscriminant (x, y, "PredictorNames", PredictorNames);
%! sigma = [6.2500, 8.2500, 10.2500; ...
%!          8.2500, 11.2500, 14.2500; ...
%!          10.2500, 14.2500, 18.2500];
%! mu = [2.5000, 3.5000, 4.5000; ...
%!       5.0000, 5.0000, 5.0000];
%! xCentered = [-1.5000, -1.5000, -1.5000; ...
%!              1.5000, 1.5000, 1.5000; ...
%!              2.0000, 3.0000, 4.0000; ...
%!             -2.0000, -3.0000, -4.0000];
%! assert (class (a), "ClassificationDiscriminant");
%! assert ({a.X, a.Y, a.NumObservations}, {x, y, 4})
%! assert ({a.DiscrimType, a.ResponseName}, {"linear", "Y"})
%! assert ({a.Gamma, a.MinGamma}, {1e-15, 1e-15})
%! assert (a.ClassNames, {'a'; 'b'})
%! assert (a.Sigma, sigma, 1e-11)
%! assert (a.Mu, mu)
%! assert (a.XCentered, xCentered)
%! assert (a.LogDetSigma, -29.369, 1e-4)
%! assert (a.PredictorNames, PredictorNames)
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = ClassificationDiscriminant (x, y, "Gamma", 0.5);
%! sigma = [6.2500, 4.1250, 5.1250; ...
%!          4.1250, 11.2500, 7.1250; ...
%!          5.1250, 7.1250, 18.2500];
%! mu = [2.5000, 3.5000, 4.5000; ...
%!       5.0000, 5.0000, 5.0000];
%! xCentered = [-1.5000, -1.5000, -1.5000; ...
%!              1.5000, 1.5000, 1.5000; ...
%!              2.0000, 3.0000, 4.0000; ...
%!             -2.0000, -3.0000, -4.0000];
%! assert (class (a), "ClassificationDiscriminant");
%! assert ({a.X, a.Y, a.NumObservations}, {x, y, 4})
%! assert ({a.DiscrimType, a.ResponseName}, {"linear", "Y"})
%! assert ({a.Gamma, a.MinGamma}, {0.5, 0})
%! assert (a.ClassNames, {'a'; 'b'})
%! assert (a.Sigma, sigma)
%! assert (a.Mu, mu)
%! assert (a.XCentered, xCentered)
%! assert (a.LogDetSigma, 6.4940, 1e-4)

## Test input validation for constructor
%!shared X, Y, MODEL
%! X = rand (10,2);
%! Y = [ones(5,1);2*ones(5,1)];
%! MODEL = ClassificationDiscriminant (X, Y);
%!error<ClassificationDiscriminant: too few input arguments.> ClassificationDiscriminant ()
%!error<ClassificationDiscriminant: too few input arguments.> ...
%! ClassificationDiscriminant (ones(4, 1))
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
%!error<ClassificationDiscriminant: 'ClassNames' must be a cellstring, logical or numeric vector.> ...
%! ClassificationDiscriminant (X, Y, "ClassNames", @(x)x)
%!error<ClassificationDiscriminant: 'ClassNames' must be a cellstring, logical or numeric vector.> ...
%! ClassificationDiscriminant (X, Y, "ClassNames", ['a'])
%!error<ClassificationDiscriminant: not all 'ClassNames' are present in Y.> ...
%! ClassificationDiscriminant (X, ones (10,1), "ClassNames", [1, 2])
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
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ["a"; "a"; "b"; "b"];
%! a = fitcdiscr (x, y, "Gamma", 0.5);
%! [label, score, cost] = predict (a, x);
%! l = {'a'; 'a'; 'b'; 'b'};
%! s = [0.7642, 0.2358; 0.5011, 0.4989; ...
%!      0.2375, 0.7625; 0.4966, 0.5034];
%! c = [0.2358, 0.7642; 0.4989, 0.5011; ...
%!      0.7625, 0.2375; 0.5034, 0.4966];
%! assert (label, l)
%! assert (score, s, 1e-4)
%! assert (cost, c, 1e-4)
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcdiscr (x, y);
%! [label, score, cost] = predict (obj, xc);
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
%!error<ClassificationDiscriminant.predict: XC must have the same number of features as the trained model.> ...
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
%! assert (class (CVMdl.Trained{1}), "ClassificationDiscriminant")
%!test
%! CVMdl = crossval (obj, "KFold", 5);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (CVMdl.KFold == 5)
%! assert (class (CVMdl.Trained{1}), "ClassificationDiscriminant")
%!test
%! CVMdl = crossval (obj, "HoldOut", 0.2);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "ClassificationDiscriminant")
%!test
%! CVMdl = crossval (obj, "LeaveOut", 'on');
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert (class (CVMdl.Trained{1}), "ClassificationDiscriminant")
%!test
%! partition = cvpartition (size (x, 1), 'KFold', 3);
%! CVMdl = crossval (obj, 'cvPartition', partition);
%! assert (class (CVMdl), "ClassificationPartitionedModel")
%! assert (CVMdl.KFold == 3)
%! assert (class (CVMdl.Trained{1}), "ClassificationDiscriminant")

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
