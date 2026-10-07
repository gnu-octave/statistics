## Copyright (C) 2023 Mohammed Azmat Khan <azmat.dev0@gmail.com>
## Copyright (C) 2023-2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{obj} =} fitrgam (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} fitrgam (@var{Tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{obj} =} fitrgam (@var{Tbl}, @var{formula})
## @deftypefnx {statistics} {@var{obj} =} fitrgam (@var{Tbl}, @var{Y})
## @deftypefnx {statistics} {@var{obj} =} fitrgam (@var{X}, @var{Y}, @var{name}, @var{value})
##
## Fit a Generalized Additive Model (GAM) for regression.
##
## @code{@var{obj} = fitrgam (@var{X}, @var{Y})} returns an object of
## class RegressionGAM, with matrix @var{X} containing the predictor data and
## vector @var{Y} containing the continuous response data.
##
## @itemize
## @item
## @var{X} must be a @math{N*P} numeric matrix of input data where rows
## correspond to observations and columns correspond to features or variables.
## @var{X} will be used to train the GAM model.
## @item
## @var{Y} must be @math{N*1} numeric vector containing the response data
## corresponding to the predictor data in @var{X}. @var{Y} must have same
## number of rows as @var{X}.
## @end itemize
##
## @var{Tbl} may stand in place of @var{X}, the response named by one of its
## variables, written into a model formula holding main effects,
## @qcode{'Y ~ x1 + x2'}, or given beside it as @var{Y}.  A variable holding
## levels rather than numbers is a categorical predictor without being named
## one, and @qcode{'CategoricalPredictors'} adds to that set rather than
## replacing it.  @code{PredictorNames} and @code{ResponseName} come from the
## table.
##
## @code{@var{obj} = fitrgam (@dots{}, @var{name}, @var{value})} returns
## an object of class RegressionGAM with additional properties specified by
## @qcode{Name-Value} pair arguments listed below.
##
## @multitable @columnfractions 0.2 0.75
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'FitMethod'} @tab A character vector selecting the weak
## learner, either @qcode{'boostedtrees'} or @qcode{'splines'}.  The default
## is @qcode{'boostedtrees'}, which boosts one shallow decision tree per
## predictor and is the scheme MATLAB uses.  @qcode{'splines'} boosts a
## smoothing spline per predictor instead and is an Octave extension.  The
## two take different options and an option meant for one is refused by the
## other rather than ignored, so the rows below say which engine each
## belongs to.
##
## @item @qcode{'predictors'} @tab Predictor Variable names, specified as
## a row vector cell of strings with the same length as the columns in @var{X}.
## If omitted, the program will generate default variable names
## @qcode{(x1, x2, ..., xn)} for each column in @var{X}.
##
## @item @qcode{'responsename'} @tab Response Variable Name, specified as
## a string.  If omitted, the default value is @qcode{'Y'}.
##
## @item @qcode{'formula'} @tab (spline option) a model specification given as a
## string in
## the form @qcode{'Y ~ terms'} where @qcode{Y} represents the response variable
## and @qcode{terms} the predictor variables.  The formula can be used to
## specify a subset of variables for training model.  For example:
## @qcode{'Y ~ x1 + x2 + x3 + x4 + x1:x2 + x2:x3'} specifies four linear terms
## for the first four columns of for predictor data, and @qcode{x1:x2} and
## @qcode{x2:x3} specify the two interaction terms for 1st-2nd and 3rd-4th
## columns respectively.  Only these terms will be used for training the model,
## but @var{X} must have at least as many columns as referenced in the formula.
## If Predictor Variable names have been defined, then the terms in the formula
## must reference to those.  When @qcode{'formula'} is specified, all terms used
## for training the model are referenced in the @qcode{IntMatrix} field of the
## @var{obj} class object as a matrix containing the column indexes for each
## term including both the predictors and the interactions used.
##
## @item @qcode{'interactions'} @tab a logical matrix, a positive integer
## scalar, or the string @qcode{'all'} for defining the interactions between
## predictor variables.  When given a logical matrix, it must have the same
## number of columns as @var{X} and each row corresponds to a different
## interaction term combining the predictors indexed as @qcode{true}.  Each
## interaction term is appended as a column vector after the available predictor
## column in @var{X}.  When @qcode{'all'} is defined, then all possible
## combinations of interactions are appended in @var{X} before training.  At the
## moment, parsing a positive integer has the same effect as the @qcode{'all'}
## option.  When @qcode{'interactions'} is specified, only the interaction terms
## appended to @var{X} are referenced in the @qcode{IntMatrix} field of the
## @var{obj} class object.
##
## @item @qcode{'knots'} @tab (spline option) a scalar or a row vector with the
## same
## columns as @var{X}.  It defines the knots for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 5, hence expanded to @qcode{ones (1, columns (X)) * 5}.  You can
## parse a row vector with different number of knots for each predictor
## variable to be fitted with, although not recommended.
##
## @item @qcode{'order'} @tab (spline option) a scalar or a row vector with the
## same
## columns as @var{X}.  It defines the order of the polynomial when training the
## GAM.  As a scalar, it is expanded to a row vector.  The default values is 3,
## hence expanded to @qcode{ones (1, columns (X)) * 3}.  You can parse a row
## vector with different number of polynomial order for each predictor variable
## to be fitted with, although not recommended.
##
## @item @qcode{'dof'} @tab (spline option) a scalar or a row vector with the
## same columns
## as @var{X}.  It defines the degrees of freedom for fitting a polynomial when
## training the GAM.  As a scalar, it is expanded to a row vector.  The default
## value is 8, hence expanded to @qcode{ones (1, columns (X)) * 8}.  You can
## parse a row vector with different degrees of freedom for each predictor
## variable to be fitted with, although not recommended.
##
## @item @qcode{'tol'} @tab (spline option) a positive scalar to set the
## tolerance for
## convergence during training. By default, it is set to @qcode{1e-3}.
## @end multitable
##
## The rows above marked as spline options require
## @qcode{'FitMethod', 'splines'}.  The remaining options belong to the
## boosted-tree engine and require @qcode{'FitMethod', 'boostedtrees'}, which
## is the default.
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'NumTreesPerPredictor'} @tab A positive integer, the number of
## boosting rounds of the predictor phase.  It is a budget rather than a
## count: a fit that stops improving ends earlier and reports so.  The default
## is 300.
##
## @item @qcode{'NumTreesPerInteraction'} @tab A positive integer, the same
## budget for the interaction phase.  The default is 100.
##
## @item @qcode{'MaxNumSplitsPerPredictor'} @tab A positive integer, the
## largest number of splits any one predictor tree may make.  The default is
## 1, which makes each tree a stump.
##
## @item @qcode{'MaxNumSplitsPerInteraction'} @tab The same limit for a tree
## over a pair of predictors.  The default is 4.
##
## @item @qcode{'InitialLearnRateForPredictors'} @tab A value greater than 0
## and at most 1, the step a round of the predictor phase starts at.  A round
## that fails to improve the fit is retried at half the step, so this is an
## initial value rather than a fixed one.  The default is 1.
##
## @item @qcode{'InitialLearnRateForInteractions'} @tab The same for the
## interaction phase.  The default is 1.
##
## @item @qcode{'MaxPValue'} @tab A value between 0 and 1.  A candidate pair
## of predictors is kept only if its interaction test gives a @math{p}-value
## no larger than this.  The default is 1, which keeps every pair asked for.
##
## @item @qcode{'Verbose'} @tab A non-negative integer.  Greater than zero
## prints a trace of the fit.  The default is 0.
##
## @item @qcode{'NumPrint'} @tab A positive integer, how often the trace
## reports: the first round and then every @var{NumPrint} rounds.  The
## default is 10.
##
## @item @qcode{'CategoricalPredictors'} @tab The predictors to treat as
## categorical: a vector of column indices, a logical vector with one element
## per predictor, or @qcode{'all'}.  A tree splits a categorical predictor
## into two sets of levels, and a level not seen in training predicts as a
## missing value.  Its @qcode{BinEdges} and @qcode{PairDetectionBinEdges} are
## empty.  The default is none.
## A predictor may be named rather than indexed, as a character matrix of one
## padded name per row, a string array or a cellstr; a name must match an entry
## of @qcode{'PredictorNames'} exactly, its case included.
##
## @item @qcode{'Weights'} @tab A single or double vector of non-negative
## observation weights, one per row of @var{X}.  Only their proportions matter.
## The model's @code{W} keeps the class of the weights, while every computation
## runs in double, so the predictions are double where MATLAB returns single.
##
## @end multitable
##
## You can parse either a @qcode{'formula'} or an @qcode{'interactions'}
## optional parameter.  Parsing both parameters will result an error.
## Accordingly, you can only pass up to two parameters among @qcode{'knots'},
## @qcode{'order'}, and @qcode{'dof'} to define the required polynomial for
## training the GAM model.
##
## @seealso{RegressionGAM, regress, regress_gp}
## @end deftypefn

function obj = fitrgam (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitrgam: too few arguments.");
  endif
  if (mod (nargin, 2) != 0)
    error ("fitrgam: Name-Value arguments must be in pairs.");
  endif

  ## Check predictor data and labels have equal rows
  if (! istable (X) && rows (X) != rows (Y))
    error ("fitrgam: number of rows in X and Y must be equal.");
  endif
  ## Parse arguments to class def function
  obj = RegressionGAM (X, Y, varargin{:});

endfunction

%!demo
%! rng (42);
%! # Train a RegressionGAM Model for synthetic values
%!
%! f1 = @(x) cos (3 *x);
%! f2 = @(x) x .^ 3;
%!
%! # generate x1 and x2 for f1 and f2
%! x1 = 2 * rand (50, 1) - 1;
%! x2 = 2 * rand (50, 1) - 1;
%!
%! # calculate y
%! y = f1(x1) + f2(x2);
%!
%! # add noise
%! y = y + y .* 0.2 .* rand (50,1);
%! X = [x1, x2];
%!
%! # create an object
%! a = fitrgam (X, y, 'FitMethod', 'splines', 'tol', 1e-3)


## Test constructor
%!demo
%! ## Fit from a table, and predict on one
%!
%! load fisheriris
%! T = table (meas(:,2), meas(:,3), meas(:,4), meas(:,1), ...
%!            'VariableNames', {'SW', 'PL', 'PW', 'SL'});
%!
%! ## A column holding levels is a categorical predictor without being named
%! ## one
%! T.Wide = categorical (meas(:,2) > 3, [false true], {'narrow', 'wide'});
%!
%! ## A model formula names the response and the predictors together, and
%! ## holds main effects only
%! Mdl = fitrgam (T, 'SL ~ PL + Wide');
%! Mdl.PredictorNames
%! Mdl.ResponseName
%!
%! ## predict matches the table's variables by name, so a column the model
%! ## was not fitted on is passed over
%! yFit = predict (Mdl, T(1:5,:));
%! yFit'

%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = [1; 2; 3; 4];
%! a = fitrgam (x, y, 'FitMethod', 'splines');
%! assert_equal ({a.X, a.Y}, {x, y})
%! assert_equal ({a.BaseModel.Intercept}, {2.5000})
%! assert_equal ({a.Knots, a.Order, a.DoF}, {[5, 5, 5], [3, 3, 3], [8, 8, 8]})
%! assert_equal ({a.NumObservations, a.NumPredictors}, {4, 3})
%! assert_equal ({a.ResponseName, a.PredictorNames}, {'Y', {'x1', 'x2', 'x3'}})
%! assert_equal ({a.Formula}, {[]})
%!test
%! x = [1, 2, 3, 4; 4, 5, 6, 7; 7, 8, 9, 1; 3, 2, 1, 2];
%! y = [1; 2; 3; 4];
%! pnames = {'A', 'B', 'C', 'D'};
%! formula = 'Y ~ A + B + C + D + A:C';
%! intMat = logical ([1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1;1,0,1,0]);
%! a = fitrgam (x, y, 'FitMethod', 'splines', ...
%!              'predictors', pnames, 'formula', formula);
%! assert_equal (a.IntMatrix, double (intMat))
%! assert_equal ({a.ResponseName, a.PredictorNames}, {'Y', pnames})
%! assert_equal (a.Formula, formula)

## Test input validation
%!error<fitrgam: too few arguments.> fitrgam ()
%!error<fitrgam: too few arguments.> fitrgam (ones (10,2))
%!error<fitrgam: Name-Value arguments must be in pairs.>
%! fitrgam (ones (4,2), ones (4, 1), 'K')
%!error<fitrgam: number of rows in X and Y must be equal.>
%! fitrgam (ones (4,2), ones (3, 1))
%!error<fitrgam: number of rows in X and Y must be equal.>
%! fitrgam (ones (4,2), ones (3, 1), 'K', 2)

## Table input
%!shared frgT
%! load fisheriris
%! frgT = table (meas(:,2), meas(:,3), meas(:,4), meas(:,1), ...
%!               'VariableNames', {'SW', 'PL', 'PW', 'SL'});
%! frgT.Wide = categorical (meas(:,2) > 3, [false true], {'narrow', 'wide'});

%!test  # the response is named by a column and the rest are predictors
%! Mdl = fitrgam (frgT, 'SL');
%! assert_equal (Mdl.PredictorNames, {'SW', 'PL', 'PW', 'Wide'});
%! assert_equal (Mdl.ResponseName, 'SL');
%! assert_equal (Mdl.CategoricalPredictors, 4);

%!test  # a model formula names the response and the predictors together
%! Mdl = fitrgam (frgT, 'SL ~ PL + Wide');
%! assert_equal (Mdl.PredictorNames, {'PL', 'Wide'});
%! assert_equal (Mdl.CategoricalPredictors, 2);

%!test  # the response may be given beside a table of predictors
%! Mdl = fitrgam (frgT(:,1:3), frgT.SL);
%! assert_equal (Mdl.PredictorNames, {'SW', 'PL', 'PW'});
%! assert_equal (Mdl.ResponseName, 'Y');

%!test  # predict takes a table, matched by name and not by position
%! Mdl = fitrgam (frgT, 'SL');
%! a = predict (Mdl, frgT);
%! assert_equal (numel (a), 150);
%! assert_equal (predict (Mdl, frgT(:, [5, 4, 3, 2, 1])), a);

%!test  # the levels travel with the model when it is made compact
%! Mdl = fitrgam (frgT, 'SL');
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.PredictorLevels, Mdl.PredictorLevels);
%! assert_equal (predict (CMdl, frgT), predict (Mdl, frgT));

%!error<RegressionGAM: the table holds no variable 'NoSuch'.> ...
%! fitrgam (frgT, 'NoSuch')

%!error<RegressionGAM: a model formula holds main effects only, so no products, powers or wildcards.> ...
%! fitrgam (frgT, 'SL ~ PL*PW')

%!error<RegressionGAM.predict: the table holds no predictor 'PL'.> ...
%! predict (fitrgam (frgT, 'SL'), frgT(:, [1, 3, 4, 5]))
