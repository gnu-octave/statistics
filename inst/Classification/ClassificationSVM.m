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
## classification model.
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
## @item
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
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Field} @tab @tab @var{Description}
##
## @item @qcode{obj.X} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Y} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.NumObservations} @tab @tab Class labels, specified as a logical or
## numeric vector, or cell array of character vectors.  Each value in @var{Y} is
## the observed class label for the corresponding row in @var{X}.
##
## @item @qcode{obj.W} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.ModelParameters} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.PredictorNames} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.ExpandedPredictorNames} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.ResponseName} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.ClassNames} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Cost} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Prior} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.ScoreTransform} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Alpha} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Beta} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Bias} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.KernelParameters} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Mu} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Sigma} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.SupportVectors} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.SupportVectorLabels} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.BoxConstraint} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.CacheInfo} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.ConvergenceInfo} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Gradient} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.IsSupportVector} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Nu} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.NumIterations} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.OutlierFraction} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.ShrinkagePeriod} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.Solver} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.RowsUsed} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @item @qcode{obj.HyperparameterOptimizationResults} @tab @tab Unstandardized predictor data, specified as a
## numeric matrix.  Each column of @var{X} represents one predictor (variable),
## and each row represents one observation.
##
## @end multitable
##
## @seealso{fitcsvm, svmtrain, svmpredict}
## @end deftypefn

  properties (Access = public)

    X                                 = [];     # Predictor data
    Y                                 = [];     # Class labels
    W                                 = [];     # Weights of observations used to train this model

    NumObservations                   = [];     # Number of observations in training dataset
    PredictorNames                    = [];     # Predictor variables names
    ResponseName                      = [];     # Response variable name
    RowsUsed                          = [];     # Rows used in fitting
    Mu                                = [];     # Predictor means
    Sigma                             = [];     # Predictor standard deviations

    ModelParameters                   = [];     # SVM parameters.
    ExpandedPredictorNames            = [];     # Expanded predictor names
    ClassNames                        = [];     # Names of classes in Y
    Cost                              = [];     # Cost of misclassification
    Prior                             = [];     # Prior probability for each class
    ScoreTransform                    = [];     # Transformation applied to predicted classification scores

    Alpha                             = [];     # Coefficients obtained by solving the dual problem
    Beta                              = [];     # Coefficients for the primal linear problem
    Bias                              = [];     # Bias term
    KernelParameters                  = [];     # Kernel parameters

    SupportVectors                    = [];     # Support vectors
    SupportVectorLabels               = [];     # Support vector labels (+1 and -1)
    IsSupportVector                   = [];     # Indices of support vectors in the training data
    BoxConstraints                    = [];     # Box constraints
    CacheInfo                         = [];     # Cache information
    ConvergenceInfo                   = [];     # Convergence information

    Gradient                          = [];     # Gradient values in the training data
    Nu                                = [];     # Nu parameter for one-class learning
    NumIterations                     = [];     # Number of iterations taken by optimization
    OutlierFraction                   = [];     # Expected fraction of outliers in the training data
    ShrinkagePeriod                   = [];     # Number of iterations between reductions of the active set
    Solver                            = [];     # Solver used
    HyperparameterOptimizationResults = [];     # An object or table describing the results of hyperparameter optimization

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

      ## Check if X is a table
      if (istable(X))
        if (ischar(Y) && ismember(Y, X.Properties.VariableNames))
          ## Use the variable Y from the table as the response
          Y = X.(Y);
          X(:, Y) = [];
        elseif (isstring(Y))          ## if formula is given as input
          parts = strsplit(Y, '~');
          if (numel(parts) != 2)
              error("ClassificationSVM: Formula must be of the form 'y ~ x1 + x2 + ...'");
          endif
          responseVar = strtrim(parts{1});
          predictorStr = strtrim(parts{2});
          predictorVars = strsplit(predictorStr, '+');

          if (!ismember(responseVar, X.Properties.VariableNames))
             error("ClassificationSVM: Response variable not found in table.");
          endif
          for i = 1:numel(predictorVars)
              if (!ismember(predictorVars{i}, X.Properties.VariableNames))
                  error("ClassificationSVM: Predictor variable not found in table.");
              end
          endif
          ## Extract response variable
          Y = X.(responseVar);

          ## Extract predictor variables
          X = X(:, predictorVars);

        else
          error('ClassificationSVM: Invalid Y.');
      endif

      ## Check correspodence between predictors and response
      if (nsample != rows (Y))
        error ("ClassificationSVM: number of rows in X and Y must be equal.");
      endif

      ## Check if it's one-class or two-class learning
      if (numel(unique(Y)) == 1)
       learning_class = 1;
      elseif(numel(unique(Y)) == 2)
       learning_class = 2;
      else
       error ("ClassificationSVM: SVM only supports one class or two class learning.");
      endif

      ## Set default values before parsing optional parameters
      if (learning_class == 1)                #Default values for one class learning
       Alpha                  = 0.5 * ones(size(X,1),1);
       KernelFunction         = 'gaussian';
      elseif(learning_class == 2)             #Default values for two class learning
       Alpha                  = zeros(size(X,1),1);
       KernelFunction         = 'linear';
      endif

      BoxConstraint           = 1;
      CacheSize               = 1000;
      CategoricalPredictors   = ;
      ClassNames              = ;
      ClipAlphas              = true;
      Cost                    = ;
      CrossVal                = 'off';
      CVPartition             = ;
      Holdout                 = [];
      KFold                   = 10;
      Leaveout                = ;
      GapTolerance            = 0;
      DeltaGradientTolerance  = [];
      KKTTolerance            = [];
      IterationLimit          = 1e6;
      KernelScale             = 1;
      KernelOffset            = ;
      OptimizeHyperparameters = 'none';
      PolynomialOrder         = 3;
      Nu                      = 0.5;
      NumPrint                = 1000;
      OutlierFraction         = 0;
      PredictorNames          = {};           #Predictor variable names
      Prior                   = ;
      RemoveDuplicates        = false;
      ResponseName            = 'Y';           #Response variable name
      ScoreTransform          = 'none';
      Solver                  = ;
      ShrinkagePeriod         = 0;
      Standardize             = false;
      Verbose                 = 0;
      Weights                 = ones(size(X,1),1);


      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "alpha"
          Alpha = varargin{2};
            if (!isvector(alpha))
              error ("ClassificationSVM: Alpha must be a vector.");
            elseif (size(alpha, 1) != rows(X))
              error ("ClassificationSVM: Alpha must have one element per row of X.");
            elseif (any(alpha < 0))
              error ("ClassificationSVM: Alpha must be non-negative.");
            endif

          case "boxconstraint"
            BoxConstraint = varargin{2};
            if ( !(isscalar(BoxConstraint) && BoxConstraint > 0))
              error ("ClassificationSVM: BoxConstraint must be a positive scalar.");
            endif

          case "cachesize"
            CacheSize = varargin{2};
            if ( isscalar(CacheSize))
              if (CacheSize <= 0)
              error ("ClassificationSVM: CacheSize must be a positive scalar.");
            elseif (isstring(CacheSize) && tolower(CacheSize) != "maximal")
              error ("ClassificationSVM: unidentified CacheSize.");
            else
              error (strcat(["ClassificationSVM: CacheSize must be either"], ...
              [" a positive scalar or a string 'maximal'."]));
            endif

          case "categoricalpredictors"
            if (! ((isnumeric (CategoricalPredictors) && isvector (CategoricalPredictors)) ||
                  (strcmpi (CategoricalPredictors, "all") || )))
              error (strcat (["ClassificationSVM: CategoricalPredictors must be either a"], ...
                             [" numeric vector or a string."]));
            endif

          case "classnames"

          case "clipalphas"
            ClipAlphas = tolower(varargin{2});
            if (! (islogical (ClipAlphas) && isscalar (ClipAlphas)))
              error ("ClassificationSVM: ClipAlphas must be a logical scalar.");
            endif

          case "cost"
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error ("ClassificationSVM: Cost must be a numeric square matrix.");
            endif

          case "crossval"

          case "cvpartition"

          case "holdout"
            Holdout = varargin{2};
            if (! isnumeric(Holdout) && isscalar(Holdout))
              error ("ClassificationSVM: Holdout must be a numeric scalar.");
            endif
            if (Holdout < 0 || Holdout >1)
              error ("ClassificationSVM: Holdout must be between 0 and 1.");
            endif


          case "kfold"
            KFold = varargin{2};
            if (! isnumeric(KFold))
              error ("ClassificationSVM: KFold must be a numeric value.");
            endif

          case "leaveout"

          case "gaptolerance"
            GapTolerance = varargin{2};
            if (! isnumeric(GapTolerance) && isscalar(GapTolerance))
              error ("ClassificationSVM: GapTolerance must be a numeric scalar.");
            endif
            if (GapTolerance < 0)
              error ("ClassificationSVM: GapTolerance must be non-negative scalar.");
            endif

          case "deltagradienttolerance"
            DeltaGradientTolerance = varargin{2};
            if (! isnumeric(DeltaGradientTolerance))
              error (strcat(["ClassificationSVM: DeltaGradientTolerance must ], ...
              ["be a numeric value."]));
            endif
            if (GapTolerance < 0)
              error (strcat(["ClassificationSVM: DeltaGradientTolerance must ], ...
              ["be non-negative scalar."]));
            endif


          case "kkttolerance"

          case "iterationlimit"
            IterationLimit = varargin{2};
            if (! isnumeric(IterationLimit) && isscalar(IterationLimit)...
              && IterationLimit >= 0)
              error ("ClassificationSVM: IterationLimit must be a positive number.");
            endif

          case "kernelfunction"
            KernelFunction = varargin{2};
            if (!(ischar(kernelfunction) || isa(KernelFunction, 'function_handle')))
              error("ClassificationSVM: KernelFunction must be a string or a function handle.");
            endif
            if (ischar(kernelfunction))
              if (! any (strcmpi (KernelFunction, {"linear", "gaussian", "rbf", ...
                "polynomial"})))
              error ("ClassificationSVM: unsupported Kernel function.");
            endif

          case "kernelscale"

          case "kerneloffset"
            KernelOffset = varargin{2};
            if (! isnumeric(KernelOffset) && isscalar(KernelOffset)...
              && KernelOffset >= 0)
              error ("ClassificationSVM: KernelOffset must be a non-negative scalar.");
            endif

          case "optimizehyperparameters"

          case "polynomialorder"
            PolynomialOrder = varargin{2};
            if (! isnumeric(PolynomialOrder))
              error ("ClassificationSVM: PolynomialOrder must be a positive integer.");
            endif

          case "nu"
            Nu = varargin{2};
            if ( !((isscalar(Nu) && Nu > 0 ))
              error ("ClassificationSVM: Nu must be positive scalar.");
            endif

          case "numprint"
            NumPrint = varargin{2};
            if ( !((isscalar(NumPrint) && NumPrint >= 0 ))
              error ("ClassificationSVM: NumPrint must be non-negative scalar.");
            endif

          case "outlierfraction"
            OutlierFraction = varargin{2};
            if (! (isscalar(OutlierFraction) && OutlierFraction >= 0 && OutlierFraction <= 1)
              error (strcat(["ClassificationSVM: OutlierFraction must be a scalar"], ...
              [" between 0 and 1."]));
            endif
            if (OutlierFraction > 0 && learning_class == 2 && isempty(Solver))
              Solver = 'isda';
            elseif( isempty(Solver))
              Solver = 'SMO';
            endif

          case "predictornames"
            PredictorNames = varargin{2};
            if (! isempty (PredictorNames))
              if (! iscellstr (PredictorNames))
                error (strcat (["ClassificationSVM: PredictorNames must"], ...
                               [" be a cellstring array."]));
              elseif (columns (PredictorNames) != columns (X))
                error (strcat (["ClassificationSVM: PredictorNames must"], ...
                               [" have same number of columns as X."]));
              endif
            endif

          case "prior"
            Prior = varargin{2};
            if ( isstring(Prior))
              Prior = tolower(Prior);
              if (! any (strcmpi (Prior, {"empirical", "uniform"})))
                error ("ClassificationSVM: Unsupported Prior.");
              endif
            elseif(! isstruct (Prior) || ! isfield (Prior, "ClassProbs") ...
                     || ! isfield (Prior, "ClassNames"))
              error (strcat (["ClassificationSVM: Prior must be a structure"], ...
                     [" with 'ClassProbs', and 'ClassNames' fields present."]));
            endif

          case "removeduplicates"
            RemoveDuplicates = tolower(varargin{2});
            if (! (islogical (RemoveDuplicates) && isscalar (RemoveDuplicates)))
              error ("ClassificationSVM: RemoveDuplicates must be a logical scalar.");
            endif

          case "responsename"
            if (! ischar (varargin{2}))
              error ("ClassificationSVM: ResponseName must be a char string.");
            endif
            ResponseName = tolower(varargin{2});

          case "scoretransform"
            ScoreTransform = varargin{2};
            if (! any (strcmpi (ScoreTransform, {"symmetric", "invlogit", "ismax", ...
              "symmetricismax", "none", "logit", "doublelogit", "symmetriclogit", ...
              "sign"})))
            error ("ClassificationSVM: unsupported ScoreTransform function handle.");
            endif

          case "solver"
            Solver = tolower(varargin{2});
            if (! any (strcmpi (Solver, {"smo", "isda", "l1qp"})))
              error ("ClassificationSVM: Unsupported Solver.");
            endif
            if(Solver=='smo')
              if(isempty(KernelOffset))
                KernelOffset = 0;
              elseif(isempty(DeltaGradientTolerance))
                DeltaGradientTolerance = 1e-3;
              elseif(isempty(KKTTolerance))
                KKTTolerance = 0;
            elseif(Solver=='isda')
              if(isempty(KernelOffset))
                KernelOffset = 0.1;
              elseif(isempty(DeltaGradientTolerance))
                DeltaGradientTolerance = 0;
              elseif(isempty(KKTTolerance))
                KKTTolerance = 1e-3;
            endif

          case "shrinkageperiod"

          case "standardize"
            Standardize = tolower(varargin{2});
            if (! (islogical (Standardize) && isscalar (Standardize)))
              error ("ClassificationSVM: Standardize must be a logical scalar.");
            endif

          case "verbose"
            Verbose = varargin{2};
            if (! isscalar (Verbose) || !any (isequal (Verbose, {0, 1, 2})))
              error ("ClassificationSVM: Verbose must be either 0, 1 or 2.");
            endif


          case "weights"

          otherwise
            error (strcat (["ClassificationSVM: invalid parameter name"],...
                           [" in optional pair arguments."]));

        endswitch
        varargin (1:2) = [];
      endwhile
    endfunction
   endmethods
endclassdef


## Test input validation for constructor
%!error<ClassificationSVM: too few input arguments.> ClassificationSVM ()
%!error<ClassificationSVM: too few input arguments.> ClassificationSVM (ones(10,2))
%!error<ClassificationSVM: Y must be of the form 'y ~ x1 + x2 + ...'>
%! ClassificationSVM (table([1,2],[9,8],[5,6], 'VariableNames', {'y', 'x1', 'x2'}), 'y x1 + x2');
%!error<ClassificationSVM: Y must be of the form 'y ~ x1 + x2 + ...'>
%! ClassificationSVM (table([1,2],[9,8],[5,6], 'VariableNames', {'y', 'x1', 'x2'}), 'x1 + x2');
%!error<ClassificationSVM: Y must be of the form 'y ~ x1 + x2 + ...'>
%! ClassificationSVM (table([1,2],[9,8],[5,6], 'VariableNames', {'y', 'x1', 'x2'}), 'y ~ ');
%!error<ClassificationSVM: Response variable not found in table.>
%! ClassificationSVM (table([1,2],[9,8],[5,6], 'VariableNames', {'y', 'x1', 'x2'}), 'y1 ~ x1 + x2');
%!error<ClassificationSVM: Predictor variable not found in table.>
%! ClassificationSVM (table([1,2],[9,8],[5,6], 'VariableNames', {'y', 'x1', 'x2'}), 'y ~ x1 + x3');
%!error<ClassificationSVM: Invalid Y.>
%! ClassificationSVM (table([1,2],[9,8],[5,6], 'VariableNames', {'y', 'x1', 'x2'}), 1);
%!error<ClassificationSVM: number of rows in X and Y must be equal.> ...
%! ClassificationSVM (ones(10,2), ones (5,1))
%!error<ClassificationSVM: SVM only supports one class or two class learning.>
%! ClassificationSVM (ones(10,2), ones (10,3))
%!error<ClassificationSVM: Alpha must be a vector.>
%! ClassificationSVM (ones(10,2), ones (10,1), "Alpha", 1)
%!error<ClassificationSVM: Alpha must have one element per row of X.>
%! ClassificationSVM (ones(10,2), ones (10,1), "Alpha", ones(5,1))
%!error<ClassificationSVM: Alpha must be non-negative.>
%! ClassificationSVM (ones(10,2), ones (10,1), "Alpha", -1)
%!error<ClassificationSVM: BoxConstraint must be a positive scalar.>
%! ClassificationSVM (ones(10,2), ones (10,1), "BoxConstraint", -1)
%!error<ClassificationSVM: CacheSize must be a positive scalar.>
%! ClassificationSVM (ones(10,2), ones (10,1), "CacheSize", -100)
%!error<ClassificationSVM: unidentified CacheSize.>
%! ClassificationSVM (ones(10,2), ones (10,1), "CacheSize", 'some')
%!error<ClassificationSVM: CacheSize must be either a positive scalar or a string 'maximal'>
%! ClassificationSVM (ones(10,2), ones (10,1), "CacheSize", [1,2])

%!error<ClassificationSVM: KernelFunction must be a string or a function handle.>
%! ClassificationSVM (ones(10,2), ones (10,1), "KernelFunction",[1,2])
%!error<ClassificationSVM: unsupported Kernel function.>
%! ClassificationSVM (ones(10,2), ones (10,1), "KernelFunction","some")



%!error<ClassificationSVM: PredictorNames must be a cellstring array.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "PredictorNames", -1)
%!error<ClassificationSVM: PredictorNames must be a cellstring array.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "PredictorNames", ['a','b','c'])
%!error<ClassificationSVM: PredictorNames must have same number of columns as X.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "PredictorNames", {'a','b','c'})
%!error<ClassificationSVM: invalid parameter name in optional pair arguments.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "some", "some")

%!error<ClassificationSVM: invalid values in X.> ...
%! ClassificationSVM ([1;2;3;"a";4], ones (5,1))

%!error<ClassificationSVM: Formula must be a string.>
%! ClassificationSVM (ones(10,2), ones (10,1), "formula", {"y~x1+x2"})
%!error<ClassificationSVM: Formula must be a string.>
%! ClassificationSVM (ones(10,2), ones (10,1), "formula", [0, 1, 0])
%!error<ClassificationSVM: invalid syntax in Formula.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "formula", "something")
%!error<ClassificationSVM: no predictor terms in Formula.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "formula", "something~")
%!error<ClassificationSVM: no predictor terms in Formula.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "formula", "something~")
%!error<ClassificationSVM: some predictors have not been identified> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "formula", "something~x1:")
%!error<ClassificationSVM: invalid Interactions parameter.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "interactions", "some")
%!error<ClassificationSVM: invalid Interactions parameter.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "interactions", -1)
%!error<ClassificationSVM: invalid Interactions parameter.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "interactions", [1 2 3 4])
%!error<ClassificationSVM: number of interaction terms requested is larger than> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "interactions", 3)
%!error<ClassificationSVM: Formula has been already defined.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "formula", "y ~ x1 + x2", "interactions", 1)
%!error<ClassificationSVM: Interactions have been already defined.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "interactions", 1, "formula", "y ~ x1 + x2")
%!error<ClassificationSVM: invalid value for Knots.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "knots", "a")
%!error<ClassificationSVM: DoF and Order have been set already.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "order", 3, "dof", 2, "knots", 5)
%!error<ClassificationSVM: invalid value for DoF.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "dof", 'a')
%!error<ClassificationSVM: Knots and Order have been set already.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "knots", 5, "order", 3, "dof", 2)
%!error<ClassificationSVM: invalid value for Order.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "order", 'a')
%!error<ClassificationSVM: DoF and Knots have been set already.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "knots", 5, "dof", 2, "order", 2)
%!error<ClassificationSVM: Tolerance must be a Positive scalar.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "tol", -1)
%!error<ClassificationSVM: ResponseName must be a char string.> ...
%! ClassificationSVM (ones(10,2), ones (10,1), "responsename", -1)


## Test input validation for predict method
%!error<ClassificationSVM.predict: too few arguments.> ...
%! predict (ClassificationSVM (ones(10,1), ones(10,1)))
%!error<ClassificationSVM.predict: Xfit is empty.> ...
%! predict (ClassificationSVM (ones(10,1), ones(10,1)), [])
%!error<ClassificationSVM.predict: Xfit must have the same number of features> ...
%! predict (ClassificationSVM(ones(10,2), ones(10,1)), 2)
%!error<ClassificationSVM.predict: invalid NAME in optional pairs of arguments.> ...
%! predict (ClassificationSVM(ones(10,2), ones(10,1)), ones (10,2), "some", "some")
%!error<ClassificationSVM.predict: includeinteractions must be a logical value.> ...
%! predict (ClassificationSVM(ones(10,2), ones(10,1)), ones (10,2), "includeinteractions", "some")
%!error<ClassificationSVM.predict: includeinteractions must be a logical value.> ...
%! predict (ClassificationSVM(ones(10,2), ones(10,1)), ones (10,2), "includeinteractions", 5)





