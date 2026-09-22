## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{Mdl} =} fitrgp (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitrgp (@var{Tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{Mdl} =} fitrgp (@var{Tbl}, @var{formula})
## @deftypefnx {statistics} {@var{Mdl} =} fitrgp (@var{Tbl}, @var{Y})
## @deftypefnx {statistics} {@var{Mdl} =} fitrgp (@dots{}, @var{name}, @var{value})
##
## Fit a Gaussian process regression model.
##
## @code{@var{Mdl} = fitrgp (@var{X}, @var{Y})} returns a @qcode{RegressionGP}
## object fitted to the predictor data @var{X} and the continuous response
## @var{Y}, where @var{X} is an @math{NxP} numeric matrix and @var{Y} an
## @math{Nx1} numeric vector with as many rows as @var{X}.
##
## @code{@var{Mdl} = fitrgp (@dots{}, @var{name}, @var{value})} passes the
## given @qcode{Name-Value} pairs to the model.  They are documented under
## @code{RegressionGP}, and the ones most often wanted are
## @qcode{'KernelFunction'}, @qcode{'BasisFunction'}, @qcode{'Standardize'},
## @qcode{'Sigma'} and @qcode{'FitMethod'}.
##
## When any of @qcode{'CrossVal'}, @qcode{'KFold'}, @qcode{'Holdout'},
## @qcode{'Leaveout'} or @qcode{'CVPartition'} is given, a cross validated
## model is returned instead, as a @qcode{RegressionPartitionedModel}.  Only
## one of them may be given at a time.
##
## @seealso{RegressionGP, CompactRegressionGP, RegressionPartitionedModel}
## @end deftypefn

function Mdl = fitrgp (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitrgp: too few input arguments.");
  endif
  if (mod (numel (varargin), 2) != 0)
    error ("fitrgp: name-value arguments must be in pairs.");
  endif

  ## Pull out the cross validation options, which the model does not take:
  ## they say what to do with the model once it exists, not how to fit it.
  CrossVal = false;
  cvArgs = {};
  given = 0;
  args = {};
  while (numel (varargin) > 0)
    switch (lower (varargin{1}))
      case 'crossval'
        val = varargin{2};
        if (! (ischar (val) && any (strcmpi (val, {'on', 'off'}))))
          error ("fitrgp: 'CrossVal' must be either 'on' or 'off'.");
        endif
        CrossVal = strcmpi (val, 'on');
      case {'kfold', 'holdout', 'leaveout', 'cvpartition'}
        cvArgs = [cvArgs, varargin(1:2)];
        CrossVal = true;
        given++;
      otherwise
        args = [args, varargin(1:2)];
    endswitch
    varargin(1:2) = [];
  endwhile

  if (given > 1)
    error (strcat ("fitrgp: you can use only one of 'KFold', 'Holdout',", ...
                   " 'Leaveout', or 'CVPartition' options."));
  endif

  Mdl = RegressionGP (X, Y, args{:});

  if (CrossVal)
    Mdl = crossval (Mdl, cvArgs{:});
  endif

endfunction

%!demo
%! ## Fit a Gaussian process to a noisy sine and predict on a fine grid.
%! x = linspace (0, 2*pi, 30)';
%! y = sin (x) + 0.1 * cos (7*x);
%! Mdl = fitrgp (x, y)
%! xq = linspace (0, 2*pi, 5)';
%! [yq, ysd] = predict (Mdl, xq)

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
%! ## The response is named by its column, and everything else is a predictor
%! Mdl = fitrgp (T, 'SL');
%! Mdl.PredictorNames
%! Mdl.ResponseName
%! Mdl.CategoricalPredictors
%!
%! ## A model formula names them instead, holding main effects only
%! Mdl2 = fitrgp (T, 'SL ~ PL + Wide');
%! Mdl2.PredictorNames
%!
%! ## predict reads a table by the names the model was fitted on, so the
%! ## columns may come in any order
%! yFit = predict (Mdl, T(1:5, [5, 4, 3, 2, 1]));
%! yFit'

%!test
%! ## fitrgp returns the model the class constructor returns
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! M1 = fitrgp (x, y);
%! M2 = RegressionGP (x, y);
%! assert_equal (class (M1), 'RegressionGP');
%! assert_equal (M1.Beta, M2.Beta);
%! assert_equal (M1.Sigma, M2.Sigma);
%! assert_equal (predict (M1, x), predict (M2, x), 1e-14);

%!test
%! ## The options reach the model
%! x = linspace (0, 1, 15)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! Mdl = fitrgp (x, y, 'KernelFunction', 'matern52', ...
%!               'BasisFunction', 'linear', 'ResponseName', 'temp');
%! assert_equal (Mdl.KernelFunction, 'Matern52');
%! assert_equal (Mdl.BasisFunction, 'Linear');
%! assert_equal (Mdl.ResponseName, 'temp');
%! assert_equal (numel (Mdl.Beta), 2);

%!test
%! ## A cross validation option returns a partitioned model instead
%! x = linspace (0, 1, 20)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! CVMdl = fitrgp (x, y, 'KFold', 4);
%! assert_equal (class (CVMdl), 'RegressionPartitionedModel');
%! assert_equal (CVMdl.KFold, 4);
%! assert_equal (CVMdl.CrossValidatedModel, 'GP');

%!test
%! ## 'CrossVal' on gives the ten folds it defaults to
%! x = linspace (0, 1, 20)';
%! y = cos (3*x) + 0.1 * sin (11*x);
%! CVMdl = fitrgp (x, y, 'CrossVal', 'on');
%! assert_equal (class (CVMdl), 'RegressionPartitionedModel');
%! assert_equal (CVMdl.KFold, 10);

%!test
%! ## 'CrossVal' off is the model itself
%! x = linspace (0, 1, 15)';
%! y = cos (3*x);
%! assert_equal (class (fitrgp (x, y, 'CrossVal', 'off')), 'RegressionGP');

## Test input validation
%!error<fitrgp: too few input arguments.> fitrgp (ones (5, 2))
%!error<fitrgp: name-value arguments must be in pairs.> ...
%! fitrgp (ones (5, 2), ones (5, 1), 'Standardize')
%!error<fitrgp: 'CrossVal' must be either 'on' or 'off'.> ...
%! fitrgp (ones (5, 2), ones (5, 1), 'CrossVal', 5)
%!error<fitrgp: you can use only one of 'KFold', 'Holdout', 'Leaveout', or 'CVPartition' options.> ...
%! fitrgp (ones (20, 2), ones (20, 1), 'KFold', 3, 'Holdout', 0.2)

## Table input
%!shared fgpT
%! load fisheriris
%! fgpT = table (meas(:,2), meas(:,3), meas(:,4), meas(:,1), ...
%!               'VariableNames', {'SW', 'PL', 'PW', 'SL'});
%! fgpT.Wide = categorical (meas(:,2) > 3, [false true], {'narrow', 'wide'});

%!test  # the response is named by a column and the rest are predictors
%! Mdl = fitrgp (fgpT, 'SL');
%! assert_equal (Mdl.PredictorNames, {'SW', 'PL', 'PW', 'Wide'});
%! assert_equal (Mdl.ResponseName, 'SL');
%! assert_equal (Mdl.CategoricalPredictors, 4);

%!test  # a model formula names the response and the predictors together
%! Mdl = fitrgp (fgpT, 'SL ~ PL + Wide');
%! assert_equal (Mdl.PredictorNames, {'PL', 'Wide'});
%! assert_equal (Mdl.CategoricalPredictors, 2);

%!test  # the response may be given beside a table of predictors
%! Mdl = fitrgp (fgpT(:,1:3), fgpT.SL);
%! assert_equal (Mdl.PredictorNames, {'SW', 'PL', 'PW'});
%! assert_equal (Mdl.ResponseName, 'Y');

%!test  # predict takes a table, matched by name and not by position
%! Mdl = fitrgp (fgpT, 'SL');
%! a = predict (Mdl, fgpT);
%! assert_equal (numel (a), 150);
%! assert_equal (predict (Mdl, fgpT(:, [5, 4, 3, 2, 1])), a);

%!test  # the levels travel with the model when it is made compact
%! Mdl = fitrgp (fgpT, 'SL');
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.PredictorLevels, Mdl.PredictorLevels);
%! assert_equal (predict (CMdl, fgpT), predict (Mdl, fgpT));

%!error<RegressionGP: the table holds no variable 'NoSuch'.> ...
%! fitrgp (fgpT, 'NoSuch')

%!error<RegressionGP: a model formula holds main effects only, so no products, powers or wildcards.> ...
%! fitrgp (fgpT, 'SL ~ PL*PW')

%!error<RegressionGP.predict: the table holds no predictor 'PL'.> ...
%! predict (fitrgp (fgpT, 'SL'), fgpT(:, [1, 3, 4, 5]))
