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
    Solver                  = 'SMO';  # Solver used

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

      ## Check correspodence between predictors and response
      if (nsample != rows (Y))
        error ("ClassificationSVM: number of rows in X and Y must be equal.");
      endif


      SVMtype                 = 'C_SVC';
      KernelFunction          = 'rbf';
      PolynomialOrder         = 3;
      Gamma                   = 1 / (ndims_X);
      KernelOffset            = 0;
      Nu                      = 0.5;
      CacheSize               = 100;
      Epsilon                 = 1e-3;
      Shrinking               = 1;
      ProbabilityEstimates    = 0;
      Weight                  = 1;
      BoxConstraint           = 1;
      KFold                   = 10;



      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case "svmtype"
            SVMtype = tolower(varargin{2});
            if (!(ischar(SVMtype)))
              error("ClassificationSVM: SVMtype must be a string.");
            endif
            if (ischar(svmtype))
              if (! any (strcmpi (svmtype, {"c_svc", "nu_svc",  ...
                "one_class_svc"})))
              error ("ClassificationSVM: unsupported SVMtype.");
              endif
            endif

          case "kernelfunction"
            KernelFunction = tolower(varargin{2});
            if (!(ischar(kernelfunction)))
              error("ClassificationSVM: KernelFunction must be a string.");
            endif
            if (ischar(kernelfunction))
              if (! any (strcmpi (KernelFunction, {"linear", "gaussian", "rbf", ...
                "polynomial", "sigmoid", "precomputed"})))
              error ("ClassificationSVM: unsupported Kernel function.");
              endif
            endif

          case "polynomialorder"
            PolynomialOrder = varargin{2};
            if (! (isnumeric(PolynomialOrder) && isscalar(PolynomialOrder)
              && (PolynomialOrder > 0) && mod(PolynomialOrder, 1) == 0))
              error (strcat(["ClassificationSVM: PolynomialOrder must be a "], ...
              ["positive integer"]));
            endif

          case "kerneloffset"
            KernelOffset = varargin{2};
            if (! (isnumeric(KernelOffset) && isscalar(KernelOffset)
              && KernelOffset >= 0))
              error (strcat(["ClassificationSVM: KernelOffset must be a non"], ...
              ["-negative scalar."]));
            endif

          case "nu"
            Nu = varargin{2};
            if ( !((isscalar(Nu) && (Nu > 0) && (Nu <= 1))))
              error (strcat(["ClassificationSVM: Nu must be positive scalar "], ...
              ["in the range (0,1]."]));
            endif

          case "cachesize"
            CacheSize = varargin{2};
            if ( !(isscalar(CacheSize) && CacheSize > 0))
              error ("ClassificationSVM: CacheSize must be a positive scalar.");
            endif

          case "epsilon"
            epsilon = varargin{2};
            if ( !(isscalar(epsilon) && epsilon >= 0))
              error ("ClassificationSVM: Epsilon must be a positive scalar.");
            endif

          case "shrinking"
            Shrinking = varargin{2};
            if ( !ismember(Shrinking, [0, 1]) )
              error ("ClassificationSVM: Shrinking must be either 0 or 1.");
            endif

          case "probabilityestimates"
            ProbabilityEstimates = varargin{2};
            if ( !ismember(ProbabilityEstimates, [0, 1]) )
              error ( strcat(["ClassificationSVM: ProbabilityEstimates must be"], ...
             [" either 0 or 1."]));
            endif

          case "weight"
            Weight = varargin{2};
            if ( !(isscalar(Weight) && Weight > 0))
              error ("ClassificationSVM: Weight must be a positive scalar.");
            endif

          case "boxconstraint"
            BoxConstraint = varargin{2};
            if ( !(isscalar(BoxConstraint) && BoxConstraint > 0))
              error ("ClassificationSVM: BoxConstraint must be a positive scalar.");
            endif

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





