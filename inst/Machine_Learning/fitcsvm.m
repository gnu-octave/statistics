## Copyright (C) 2024 Pallav Purbia <pallavpurbia@gmail.com>
## Copyright (C) 2024-2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{Mdl} =} fitcsvm (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcsvm (@var{Tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{Mdl} =} fitcsvm (@var{Tbl}, @var{formula})
## @deftypefnx {statistics} {@var{Mdl} =} fitcsvm (@var{Tbl}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitcsvm (@dots{}, @var{name}, @var{value})
##
## Fit a Support Vector Machine classification model.
##
## @code{@var{Mdl} = fitcsvm (@var{X}, @var{Y})} returns a Support Vector
## Machine classification model, @var{Mdl}, with @var{X} being the predictor
## data, and @var{Y} the class labels of observations in @var{X}.
##
## @itemize
## @item
## @code{X} must be a @math{N*P} numeric matrix of predictor data where rows
## correspond to observations and columns correspond to features or variables.
## @item
## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels of
## corresponding predictor data in @var{X}.  @var{Y} can be numerical, logical,
## char array or cell array of character vectors.  @var{Y} must have same number
## of rows as @var{X}.
## @end itemize
##
## @code{@var{Mdl} = fitcsvm (@dots{}, @var{name}, @var{value})} returns a
## Support Vector Machine model with additional options specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @subheading Model Parameters
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Standardize'} @tab A boolean flag indicating whether
## the data in @var{X} should be standardized prior to training.
##
## @item @qcode{'CategoricalPredictors'} @tab The predictors whose values are
## levels, as indices, as a logical vector with one element per predictor, or as
## @qcode{'all'}.  Each is dummy coded in its place, one column of zeros and
## ones per level seen in training, named as in @qcode{'x1 == 2'} in
## @code{ExpandedPredictorNames}, and the coded columns are not standardized.
## An observation holding a level the training data did not has no score.
## A predictor may be named rather than indexed, as a character matrix of one
## padded name per row, a string array or a cellstr; a name must match an entry
## of @qcode{'PredictorNames'} exactly, its case included.
##
## @item @qcode{'PredictorNames'} @tab A cell array of character vectors
## specifying the predictor variable names.  The variable names are assumed to
## be in the same order as they appear in the training data @var{X}.
##
## @item @qcode{'ResponseName'} @tab A character vector specifying the name
## of the response variable.
##
## @item @qcode{'ClassNames'} @tab Names of the classes in the class
## labels, @var{Y}, used for fitting the kNN model.  @qcode{ClassNames} are of
## the same type as the class labels in @var{Y}.  The model keeps the classes
## in this order; by default they are sorted.
##
## @item @qcode{'SVMtype'} @tab Specifies the type of SVM used for training
## the @code{ClassificationSVM} model.  By default, the type of SVM is defined
## by setting other parameters and/or by the data itself.  Setting the
## @qcode{'SVMtype'} parameter overrides the default behavior and it accepts the
## following options:
## @end multitable
##
## @multitable @columnfractions 0.2 0.75
## @headitem @var{Value} @tab @var{Description}
## @item @qcode{'C_SVC'} @tab It is the standard SVM formulation for
## classification tasks.  It aims to find the optimal hyperplane that separates
## different classes by maximizing the margin between them while allowing some
## misclassifications.  The parameter @qcode{'C'} controls the trade-off between
## maximizing the margin and minimizing the classification error.  It is the
## default type, unless otherwise specified.
## @item @qcode{'nu_SVC'} @tab It is a variation of the standard SVM that
## introduces a parameter @math{ν} (nu) as an upper bound on the fraction of
## margin errors and a lower bound on the fraction of support vectors.  This
## formulation provides more control over the number of support vectors and the
## margin errors, making it useful for specific classification scenarios.  It is
## the default type, when the @qcode{'OutlierFraction'} parameter is set.
## @item @qcode{'one_class_SVM'} @tab It is used for anomaly detection and
## novelty detection tasks. It aims to separate the data points of a single
## class from the origin in a high-dimensional feature space. This method is
## particularly useful for identifying outliers or unusual patterns in the data.
## It is the default type, when the @qcode{'Nu'} parameter is set or when there
## is a single class in @var{Y}.  When @qcode{'one_class_SVM'} is set by the
## @qcode{'SVMtype'} pair argument, @var{Y} has no effect and any classes are
## ignored.
## @end multitable
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'OutlierFraction'} @tab The expected proportion of outliers
## in the training data, specified as a scalar value in the range @math{[0,1]}.
## When specified, the type of SVM model is switched to @qcode{'nu_SVC'} and
## @qcode{'OutlierFraction'} defines the @math{ν} (nu) parameter.
##
## @item @qcode{'KernelFunction'} @tab A character vector specifying the
## method for computing elements of the Gram matrix.  The available kernel
## functions are @qcode{'gaussian'} or @qcode{'rbf'}, @qcode{'linear'},
## @qcode{'polynomial'}, and @qcode{'sigmoid'}.  For one-class learning, the
## default Kernel function is @qcode{'rbf'}. For two-class learning the default
## is @qcode{'linear'}.
##
## @item @qcode{'PolynomialOrder'} @tab A positive integer that specifies
## the order of polynomial in kernel function.  The default value is 3.  Unless
## the @qcode{'KernelFunction'} is set to @qcode{'polynomial'}, this parameter
## is ignored.
##
## @item @qcode{'KernelScale'} @tab A positive scalar that specifies a
## scaling factor for the @math{γ} (gamma) parameter, which can be seen as the
## inverse of the radius of influence of samples selected by the model as
## support vectors.  The @math{γ} (gamma) parameter is computed as
## @math{gamma = @qcode{KernelScale} / (number of features)}.  The default value
## for @qcode{'KernelScale'} is 1.
##
## @item @qcode{'KernelOffset'} @tab A nonnegative scalar that specifies
## the @math{coef0} in kernel function. For the polynomial kernel, it influences
## the polynomial's shift, and for the sigmoid kernel, it affects the hyperbolic
## tangent's shift. The default value for @qcode{'KernelOffset'} is 0.
##
## @item @qcode{'Weights'} @tab A single or double vector of nonnegative
## observation weights, one per row of @var{X}.  Each observation's box
## constraint is @math{n} times @qcode{BoxConstraint} times its weight, the
## weights scaled so that each class carries its prior times the cost of
## misclassifying it.  An empirical prior sums the weights per class,
## standardization uses weighted means and standard deviations, and a row of
## zero weight is left out.  The model's @code{W} keeps the class of the
## weights, while every computation runs in double.
##
## @item @qcode{'BoxConstraint'} @tab A positive scalar that specifies the
## upper bound of the Lagrange multipliers, i.e. the parameter C, which is used
## for training @qcode{'C_SVC'} and @qcode{'one_class_SVM'} type of models.  It
## determines the trade-off between maximizing the margin and minimizing the
## classification error. The default value for @qcode{'BoxConstraint'} is 1.
##
## @item @qcode{'Nu'} @tab A positive scalar, in the range @math{(0,1]}
## that specifies the parameter @math{ν} (nu) for training @qcode{'nu_SVC'} and
## @qcode{'one_class_SVM'} type of models.  Unless overridden by setting the
## @qcode{'SVMtype'} parameter, setting the @qcode{'Nu'} parameter always forces
## the training model type to @qcode{'one_class_SVM'}, in which case, the number
## of classes in @var{Y} is ignored.  The default value for @qcode{'Nu'} is 1.
##
## @item @qcode{'CacheSize'} @tab A positive scalar that specifies the
## memory requirements (in MB) for storing the Gram matrix. The default is 1000.
##
## @item @qcode{'Tolerance'} @tab A nonnegative scalar that specifies
## the tolerance of termination criterion. The default value is 1e-6.
##
## @item @qcode{'Shrinking'} @tab Specifies whether to use shrinking
## heuristics. It accepts either 0 or 1. The default value is 1.
## @end multitable
##
## @subheading Cross Validation Options
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Crossval'} @tab Cross-validation flag specified as
## @qcode{'on'} or @qcode{'off'}.  If @qcode{'on'} is specified, a 10-fold
## cross validation is performed and a @code{ClassificationPartitionedModel} is
## returned in @var{Mdl}.  To override this cross-validation setting, use only
## one of the following Name-Value pair arguments.
##
## @item @qcode{'CVPartition'} @tab A @code{cvpartition} object that
## specifies the type of cross-validation and the indexing for the training and
## validation sets.  A @code{ClassificationPartitionedModel} is returned in
## @var{Mdl} and the trained model is stored in the @code{Trained} property.
##
## @item @qcode{'Holdout'} @tab Fraction of the data used for holdout
## validation, specified as a scalar value in the range @math{[0,1]}.  When
## specified, a randomly selected percentage is reserved as validation data and
## the remaining set is used for training.  The trained model is stored in the
## @code{Trained} property of the @code{ClassificationPartitionedModel} returned
## in @var{Mdl}.  @qcode{'Holdout'} partitioning attempts to ensure that each
## partition represents the classes proportionately.
##
## @item @qcode{'KFold'} @tab Number of folds to use in the cross-validated
## model, specified as a positive integer value greater than 1.  When specified,
## then the data is randomly partitioned in @math{k} sets and for each set, the
## set is reserved as validation data while the remaining @math{k-1} sets are
## used for training.  The trained models are stored in the @code{Trained}
## property of the @code{ClassificationPartitionedModel} returned in @var{Mdl}.
## @qcode{'KFold'} partitioning attempts to ensure that each partition
## represents the classes proportionately.
##
## @item @qcode{'Leaveout'} @tab Leave-one-out cross-validation flag
## specified as @qcode{'on'} or @qcode{'off'}.  If @qcode{'on'} is specified,
## then for each of the @math{n} observations (where @math{n} is the number of
## observations, excluding missing observations, specified in the
## @code{NumObservations} property of the model), one observation is reserved as
## validation data while the remaining observations are used for training.  The
## trained models are stored in the @code{Trained} property of the
## @code{ClassificationPartitionedModel} returned in @var{Mdl}.
## @end multitable
##
## @seealso{ClassificationSVM, ClassificationPartitionedModel, svmtrain,
## svmpredict}
## @end deftypefn

function Mdl = fitcsvm (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitcsvm: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitcsvm: Name-Value arguments must be in pairs.");
  endif

  ## Check predictor data and labels have equal rows
  if (! istable (X) && rows (X) != rows (Y))
    error ("fitcsvm: number of rows in X and Y must be equal.");
  endif

  ## Check optional input parameters for cross-validation options
  cv_opt = false;
  cv_arg = 0;
  args = {};
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))

      case 'crossval'
        CrossVal = varargin{2};
        if (! any (strcmp (CrossVal, {'off', 'on'})))
          error ("fitcsvm: 'CrossVal' must be either 'off' or 'on'.");
        endif
        if (strcmp (CrossVal, 'on'))
          cv_opt = true;
        endif

      case 'kfold'
        Name = 'KFold';
        Value = varargin{2};
        cv_arg += 1;
        cv_opt = true;

      case 'holdout'
        Name = 'Holdout';
        Value = varargin{2};
        cv_arg += 1;
        cv_opt = true;

      case 'leaveout'
        Name = 'Leaveout';
        Value = varargin{2};
        cv_arg += 1;
        cv_opt = true;

      case 'cvpartition'
        Name = 'CVPartition';
        Value = varargin{2};
        cv_arg += 1;
        cv_opt = true;

      otherwise
        args = [args, {varargin{1}, varargin{2}}];
      endswitch
    varargin(1:2) = [];
  endwhile

  ## Check for multiple cross-validation paired arguments
  if (cv_arg > 1)
    error (strcat ("fitcsvm: You can use only one cross-validation", ...
                   " name-value pair argument at a time to create a", ...
                   " cross-validated model."));
  endif

  ## Parse arguments to classdef constructor
  Mdl = ClassificationSVM (X, Y, args{:});

  ## If cross validation has been requested,
  ## return a ClassificationPartitionedModel
  if (cv_opt)
    if (cv_arg)
      Mdl = crossval (Mdl, Name, Value);
    else
      Mdl = crossval (Mdl);
    endif
  endif

endfunction

%!demo
%! ## Use a subset of Fisher's iris data set
%!
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! X = meas(inds, [3,4]);
%! Y = species(inds);
%!
%! ## Train a linear SVM classifier
%! SVMModel = fitcsvm (X, Y)
%!
%! ## Plot a scatter diagram of the data and circle the support vectors.
%! sv = SVMModel.SupportVectors;
%! figure
%! gscatter (X(:,1), X(:,2), Y)
%! hold on
%! plot (sv(:,1), sv(:,2), 'ko', 'MarkerSize', 10)
%! legend ('versicolor', 'virginica', 'Support Vector')
%! hold off

## Test constructor
%!demo
%! ## Fit from a table, and predict on one
%!
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! X = meas(inds,:);
%! T = table (X(:,1), X(:,2), X(:,3), X(:,4), ...
%!            'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! T.Species = categorical (species(inds));
%!
%! ## A column holding levels is a categorical predictor without being named
%! ## one
%! T.Wide = categorical (X(:,2) > 2.9, [false true], {'narrow', 'wide'});
%!
%! ## The response is named by its column, and everything else is a predictor
%! Mdl = fitcsvm (T, 'Species');
%! Mdl.PredictorNames
%! Mdl.CategoricalPredictors
%!
%! ## predict reads a table by the names the model was fitted on, so the
%! ## columns may come in any order
%! label = predict (Mdl, T(1:5, [6, 5, 4, 3, 2, 1]));
%! label'

%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = {'a'; 'a'; 'b'; 'b'};
%! a = fitcsvm (x, y);
%! assert_equal (class (a), "ClassificationSVM");
%! assert_equal ({a.X, a.Y}, {x, y})
%! assert_equal (a.NumObservations, 4)
%! assert_equal ({a.ResponseName, a.PredictorNames}, {'Y', {'x1', 'x2', 'x3'}})
%! assert_equal (a.ModelParameters.SVMtype, "c_svc")
%! assert_equal (a.ClassNames, {'a'; 'b'})

## Test Output
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! a = fitcsvm (x, y);
%! assert_equal (class (a), "ClassificationSVM");
%! assert_equal ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, 'linear'})
%! assert_equal (a.ModelParameters.BoxConstraint, 1)
%! assert_equal (a.ModelParameters.KernelOffset, 0)
%! assert_equal (a.ClassNames, [-1; 1])
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! a = fitcsvm (x, y, 'KernelFunction', 'rbf', 'BoxConstraint', 2, ...
%! 'KernelOffset', 2);
%! assert_equal (class (a), "ClassificationSVM");
%! assert_equal ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, 'rbf'})
%! assert_equal (a.ModelParameters.BoxConstraint, 2)
%! assert_equal (a.ModelParameters.KernelOffset, 2)
%! assert_equal (isempty (a.Alpha), false)
%! assert_equal (isempty (a.Beta), true)
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! a = fitcsvm (x, y, 'KernelFunction', 'polynomial', 'PolynomialOrder', 3);
%! assert_equal (class (a), "ClassificationSVM");
%! assert_equal ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, 'polynomial'})
%! assert_equal (a.ModelParameters.KernelPolynomialOrder, 3)
%! assert_equal (isempty (a.Alpha), false)
%! assert_equal (isempty (a.Beta), true)
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! a = fitcsvm (x, y, 'KernelFunction', 'linear', 'PolynomialOrder', 3);
%! assert_equal (class (a), "ClassificationSVM");
%! assert_equal ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, 'linear'})
%! assert_equal (isempty (a.ModelParameters.KernelPolynomialOrder), true)
%! assert_equal (isempty (a.Alpha), false)
%! assert_equal (isempty (a.Beta), false)
%! assert_equal (size (a.Beta), [2, 1])
%!test
%! x = [1, 2; 2, 3; 3, 4; 4, 5; 2, 3; 3, 4; 2, 3; 3, 4; 2, 3; 3, 4];
%! y = [1; 1; -1; -1; 1; -1; -1; -1; -1; -1];
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! a = fitcsvm (x, y, 'KernelFunction', 'linear', 'CrossVal', 'on');
%! warning (status);
%! assert_equal (class (a), "ClassificationPartitionedModel");
%! assert_equal ({a.X, a.Y, a.ModelParameters.KernelFunction}, {x, y, 'linear'})
%! assert_equal (isempty (a.ModelParameters.KernelPolynomialOrder), true)
%! assert_equal (isempty (a.Trained{1}.Alpha), false)
%! assert_equal (isempty (a.Trained{1}.Beta), false)

## Test input validation
%!error<fitcsvm: too few arguments.> fitcsvm ()
%!error<fitcsvm: too few arguments.> fitcsvm (ones (4,1))
%!error<fitcsvm: Name-Value arguments must be in pairs.>
%! fitcsvm (ones (4,2), ones (4, 1), 'KFold')
%!error<fitcsvm: number of rows in X and Y must be equal.>
%! fitcsvm (ones (4,2), ones (3, 1))
%!error<fitcsvm: number of rows in X and Y must be equal.>
%! fitcsvm (ones (4,2), ones (3, 1), 'KFold', 2)
%!error <fitcsvm: 'CrossVal' must be either 'off' or 'on'.>
%! fitcsvm (ones (4,2), ones (4, 1), 'CrossVal', 2)
%!error <fitcsvm: 'CrossVal' must be either 'off' or 'on'.>
%! fitcsvm (ones (4,2), ones (4, 1), 'CrossVal', 'a')
%!error <fitcsvm: You can use only one cross-validation name-value pair argument> ...
%! fitcsvm (ones (4,2), ones (4, 1), 'KFold', 10, 'Holdout', 0.3)

%!test  # 'Leaveout' leaves one observation out of each fold
%! load fisheriris
%! k = [51:70, 101:120];
%! CVMdl = fitcsvm (meas(k,:), species(k), 'Leaveout', 'on');
%! assert_equal (CVMdl.KFold, 40);

%!test  # MATLAB parity: classes given as text are sorted
%! load fisheriris
%! k = [101:150, 51:100];
%! Mdl = fitcsvm (meas(k,:), species(k));
%! assert_equal (Mdl.ClassNames, {'versicolor'; 'virginica'});

%!test  # MATLAB parity: a given ClassNames order is kept
%! load fisheriris
%! Mdl = fitcsvm (meas(51:150,:), species(51:150), ...
%!             'ClassNames', {'virginica'; 'versicolor'});
%! assert_equal (Mdl.ClassNames, {'virginica'; 'versicolor'});

%!test  # the score columns follow a given ClassNames order
%! load fisheriris
%! Mdl = fitcsvm (meas(51:150,:), species(51:150), ...
%!             'ClassNames', {'virginica'; 'versicolor'});
%! [label, s] = predict (Mdl, meas(51,:));
%! assert_equal (label, {'versicolor'});
%! assert_equal (s(2) > s(1), true);

%!error<ClassificationSVM: not all 'ClassNames' are present in Y.> ...
%! load fisheriris
%! fitcsvm (meas(51:150,:), species(51:150), 'ClassNames', [3, 2])

## Table input
%!shared fcsT
%! load fisheriris
%! fcsI = ! strcmp (species, 'setosa');
%! fcsX = meas(fcsI,:);
%! fcsT = table (fcsX(:,1), fcsX(:,2), fcsX(:,3), fcsX(:,4), ...
%!               'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! fcsT.Species = categorical (species(fcsI));
%! fcsT.Wide = categorical (fcsX(:,2) > 2.9, [false true], ...
%!                          {'narrow', 'wide'});

%!test  # the response is named by a column and the rest are predictors
%! Mdl = fitcsvm (fcsT, 'Species');
%! assert_equal (Mdl.PredictorNames, {'SL', 'SW', 'PL', 'PW', 'Wide'});
%! assert_equal (Mdl.ResponseName, 'Species');
%! assert_equal (Mdl.CategoricalPredictors, 5);

%!test  # a model formula names the response and the predictors together
%! Mdl = fitcsvm (fcsT, 'Species ~ PW + SL');
%! assert_equal (Mdl.PredictorNames, {'PW', 'SL'});
%! assert_equal (isempty (Mdl.CategoricalPredictors), true);

%!test  # the response may be given beside a table of predictors
%! Mdl = fitcsvm (fcsT(:,1:4), fcsT.Species);
%! assert_equal (Mdl.PredictorNames, {'SL', 'SW', 'PL', 'PW'});
%! assert_equal (Mdl.ResponseName, 'Y');

%!test  # predict takes a table, matched by name and not by position
%! Mdl = fitcsvm (fcsT, 'Species');
%! a = predict (Mdl, fcsT);
%! assert_equal (class (a), 'categorical');
%! assert_equal (predict (Mdl, fcsT(:, [6, 5, 4, 3, 2, 1])), a);

%!test  # the levels travel with the model when it is made compact
%! Mdl = fitcsvm (fcsT, 'Species');
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.PredictorLevels, Mdl.PredictorLevels);
%! assert_equal (predict (CMdl, fcsT), predict (Mdl, fcsT));

%!error<ClassificationSVM: the table holds no variable 'NoSuch'.> ...
%! fitcsvm (fcsT, 'NoSuch')

%!error<ClassificationSVM: a model formula holds main effects only, so no products, powers or wildcards.> ...
%! fitcsvm (fcsT, 'Species ~ SL*PW')

%!error<ClassificationSVM.predict: the table holds no predictor 'SW'.> ...
%! predict (fitcsvm (fcsT, 'Species'), fcsT(:, [1, 3, 4, 5, 6]))
