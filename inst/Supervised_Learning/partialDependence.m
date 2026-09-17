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
## @deftypefn  {statistics} {@var{pd} =} partialDependence (@var{Mdl}, @var{Vars})
## @deftypefnx {statistics} {@var{pd} =} partialDependence (@var{Mdl}, @var{Vars}, @var{Labels})
## @deftypefnx {statistics} {@var{pd} =} partialDependence (@dots{}, @var{Data})
## @deftypefnx {statistics} {@var{pd} =} partialDependence (@var{fun}, @var{Vars}, @var{Data})
## @deftypefnx {statistics} {@var{pd} =} partialDependence (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{pd}, @var{x}, @var{y}] =} partialDependence (@dots{})
##
## Compute partial dependence.
##
## @code{@var{pd} = partialDependence (@var{Mdl}, @var{Vars})} returns the
## partial dependence of the response of the regression model @var{Mdl} on the
## predictors named by @var{Vars}, averaged over the observations @var{Mdl}
## was fitted on.  @var{Vars} names one predictor or two, by column index or
## by name, and a model that does not keep its observations must be given
## them as @var{Data}.
##
## @code{@var{pd} = partialDependence (@var{Mdl}, @var{Vars}, @var{Labels})}
## does the same for a classification model, averaging the score of each class
## named by @var{Labels} rather than a response.  @var{Labels} is required for
## such a model and refused for any other.
##
## @code{@var{pd} = partialDependence (@var{fun}, @var{Vars}, @var{Data})}
## takes a function handle in place of a model.  @var{fun} is called with a
## matrix of observations and answers with one row for each, and @var{Data} is
## then required.
##
## @var{pd} is a @math{1xnumX} vector for a regression model varying one
## predictor and a @math{numYxnumX} matrix for two, where @var{numX} and
## @var{numY} count the query points of the first and second.  For a
## classification model it gains a leading dimension of one row per class,
## giving @math{numxnumX} and @math{numxnumYxnumX}.
##
## @var{x} and @var{y} hold the query points of the first and the second
## predictor, @var{y} empty where only one was named.  Where a predictor is
## categorical they are its levels.
##
## @multitable @columnfractions 0.28 0.02 0.7
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{'QueryPoints'} @tab @tab The points to answer at, as a vector
## for one predictor and, for two, either a matrix of one column each or a
## cell holding a vector for each, which is how they may differ in length.
## The default is 100 points evenly spaced between the smallest and the
## largest value the predictor takes among the observations averaged over,
## and the levels themselves where it is categorical.
##
## @item @qcode{'NumObservationsToSample'} @tab @tab How many observations to
## draw, without replacement, from those averaged over.  The default is all of
## them, and so is a number reaching or exceeding how many there are.  The
## default query points span what was drawn.
##
## @item @qcode{'CategoricalPredictors'} @tab @tab The predictors whose values
## are levels, taken as by every learner of this package.  It applies only to
## a function handle, a model being asked for its own.
##
## @item @qcode{'IncludeInteractions'} @tab @tab Whether the interaction terms
## of a generalized additive model are included.  It applies only to such a
## model, and the default is the model's own.
##
## @item @qcode{'IncludeIntercept'} @tab @tab Whether the intercept of a
## generalized additive model is included, @code{true} by default.  Excluding
## it takes the intercept off the result.
##
## @item @qcode{'OutputColumns'} @tab @tab Which of the columns a function
## handle answers with are wanted, as indices or as @qcode{'all'}, which is
## the default.  It applies only to a function handle.
## @end multitable
##
## @code{'UseParallel'} and @code{'PredictionForMissingValue'} are not
## implemented and are refused rather than ignored.
##
## @seealso{plotPartialDependence, PredictiveModel}
## @end deftypefn

function [pd, x, y] = partialDependence (Mdl, Vars, varargin)

  if (nargin < 2)
    error ("partialDependence: too few input arguments.");
  endif
  [pd, x, y] = pdCompute ('partialDependence', Mdl, Vars, varargin{:});

endfunction

%!shared pdX, pdYr, pdYc, pdY3, pdQ, pdLo, pdHi
%! x1 = repmat ([1;2;3], 8, 1);
%! x2 = reshape (repmat (1:8, 3, 1), 24, 1);
%! pdX = [x1, x2];
%! pdYr = 10 * x2 + x1;
%! pdYc = (x2 > 4) + 1;
%! pdY3 = repmat ({'a'; 'b'; 'c'}, 8, 1);
%! pdQ = [1; 2; 3];
%! lo1 = repmat ([1;2;3], 4, 1);
%! one = ones (12, 1);
%! pdLo = [lo1, one];
%! pdHi = [lo1, 8 * one];

## MATLAB parity: the values of a linear model, whose fit we reproduce exactly
%!test
%! Mdl = fitlm (pdX, pdYr);
%! assert_equal (partialDependence (Mdl, 1, 'QueryPoints', pdQ), ...
%!               [46, 47, 48], 1e-10);

%!test  # and over the observations it is handed rather than its own
%! Mdl = fitlm (pdX, pdYr);
%! assert_equal (partialDependence (Mdl, 1, pdLo, 'QueryPoints', pdQ), ...
%!               [11, 12, 13], 1e-10);

%!test  # MATLAB parity: two variables give one row per point of the second
%! Mdl = fitlm (pdX, pdYr);
%! [pd, x, y] = partialDependence (Mdl, [1, 2], 'QueryPoints', {pdQ, [2;5]});
%! assert_equal (pd, [21, 22, 23; 51, 52, 53], 1e-10);
%! assert_equal (size (pd), [2, 3]);
%! assert_equal (x', [1, 2, 3]);
%! assert_equal (y', [2, 5]);

## MATLAB parity: a tree answers over the distribution it was fitted on
%!test
%! Mdl = fitrtree (pdX, pdYr);
%! assert_equal (partialDependence (Mdl, 1, 'QueryPoints', pdQ), ...
%!               [47, 47, 47], 1e-10);

%!test  # and does not move when the observations it is handed are replaced
%! Mdl = fitrtree (pdX, pdYr);
%! lo = partialDependence (Mdl, 1, pdLo, 'QueryPoints', pdQ);
%! hi = partialDependence (Mdl, 1, pdHi, 'QueryPoints', pdQ);
%! assert_equal (lo, hi);
%! assert_equal (lo, [47, 47, 47], 1e-10);

%!test  # a compact tree keeps no observations and answers just the same
%! Mdl = compact (fitrtree (pdX, pdYr));
%! assert_equal (partialDependence (Mdl, 1, pdLo, 'QueryPoints', pdQ), ...
%!               [47, 47, 47], 1e-10);

%!test  # the traversal equals averaging the tree over its training data
%! Mdl = fitrtree (pdX, pdYr);
%! pd = partialDependence (Mdl, 1, 'QueryPoints', pdQ);
%! avg = zeros (1, 3);
%! for k = 1:3
%!   Z = pdX;
%!   Z(:,1) = pdQ(k);
%!   avg(k) = mean (predict (Mdl, Z));
%! endfor
%! assert_equal (pd, avg, 1e-10);

## MATLAB parity: a categorical predictor is answered at its levels
%!test
%! xc = repmat ([1;2;3;4], 6, 1);
%! x2 = reshape (repmat (1:8, 3, 1), 24, 1);
%! y = 10 * (xc == 1) + 20 * (xc == 2) + 30 * (xc == 3) + 40 * (xc == 4) + x2;
%! Mdl = fitrtree ([xc, x2], y, 'CategoricalPredictors', 1);
%! [pd, x] = partialDependence (Mdl, 1);
%! assert_equal (x', [1, 2, 3, 4]);
%! assert_equal (pd, [14, 73/3, 104/3, 45], 1e-10);

%!test  # and takes its levels whatever query points are asked for
%! xc = repmat ([1;2;3;4], 6, 1);
%! x2 = reshape (repmat (1:8, 3, 1), 24, 1);
%! y = 10 * (xc == 1) + 20 * (xc == 2) + 30 * (xc == 3) + 40 * (xc == 4) + x2;
%! Mdl = fitrtree ([xc, x2], y, 'CategoricalPredictors', 1);
%! [a, ax] = partialDependence (Mdl, 1);
%! [b, bx] = partialDependence (Mdl, 1, 'QueryPoints', [1; 2]);
%! assert_equal (a, b);
%! assert_equal (ax, bx);

## MATLAB parity: a generalized additive model answers over its own
%!test
%! Mdl = fitrgam (pdX, pdYr);
%! lo = partialDependence (Mdl, 1, pdLo, 'QueryPoints', pdQ);
%! hi = partialDependence (Mdl, 1, pdHi, 'QueryPoints', pdQ);
%! assert_equal (lo, hi);
%! assert_equal (lo, [46, 47, 48], 1e-6);

%!test  # excluding the intercept takes it off the result
%! Mdl = fitrgam (pdX, pdYr);
%! with = partialDependence (Mdl, 1, 'QueryPoints', pdQ);
%! without = partialDependence (Mdl, 1, 'QueryPoints', pdQ, ...
%!                              'IncludeIntercept', false);
%! assert_equal (with - Mdl.Intercept, without, 1e-10);

## MATLAB parity: a bagged ensemble is answered over its trees, a boosted
## one over the observations
%!test
%! Mdl = fitrensemble (pdX, pdYr, 'Method', 'Bag', 'NumLearningCycles', 5);
%! lo = partialDependence (Mdl, 1, pdLo, 'QueryPoints', pdQ);
%! hi = partialDependence (Mdl, 1, pdHi, 'QueryPoints', pdQ);
%! assert_equal (lo, hi);

%!test  # a compacted bagged ensemble carries no Method and is still bagged
%! Mdl = compact (fitrensemble (pdX, pdYr, 'Method', 'Bag', ...
%!                              'NumLearningCycles', 5));
%! lo = partialDependence (Mdl, 1, pdLo, 'QueryPoints', pdQ);
%! hi = partialDependence (Mdl, 1, pdHi, 'QueryPoints', pdQ);
%! assert_equal (lo, hi);

%!test  # a boosted ensemble moves with the observations it is handed
%! Mdl = fitrensemble (pdX, pdYr, 'Method', 'LSBoost', 'NumLearningCycles', 5);
%! lo = partialDependence (Mdl, 1, pdLo, 'QueryPoints', pdQ);
%! hi = partialDependence (Mdl, 1, pdHi, 'QueryPoints', pdQ);
%! assert_equal (isequal (lo, hi), false);

## Shapes
%!test  # a regression model varying one predictor gives a row
%! Mdl = fitrsvm (pdX, pdYr);
%! assert_equal (size (partialDependence (Mdl, 1, 'QueryPoints', pdQ)), [1, 3]);

%!test  # a classifier gives one row per class named
%! Mdl = fitctree (pdX, pdY3);
%! pd3 = partialDependence (Mdl, 1, {'a', 'b', 'c'}, 'QueryPoints', pdQ);
%! pd1 = partialDependence (Mdl, 1, 'b', 'QueryPoints', pdQ);
%! assert_equal (size (pd3), [3, 3]);
%! assert_equal (size (pd1), [1, 3]);

%!test  # and two predictors add a dimension between them
%! Mdl = fitctree (pdX, pdY3);
%! pd = partialDependence (Mdl, [1, 2], {'a', 'b'}, ...
%!                         'QueryPoints', {pdQ, [2; 5]});
%! assert_equal (size (pd), [2, 2, 3]);

%!test  # MATLAB parity: the rows follow the order the classes were named in
%! Mdl = fitctree (pdX, pdY3);
%! abc = partialDependence (Mdl, 1, {'a', 'b', 'c'}, 'QueryPoints', pdQ);
%! ca = partialDependence (Mdl, 1, {'c', 'a'}, 'QueryPoints', pdQ);
%! assert_equal (ca(1,:), abc(3,:));
%! assert_equal (ca(2,:), abc(1,:));

## Query points and sampling
%!test  # MATLAB parity: a hundred points spanning what the observations hold
%! Mdl = fitrsvm (pdX, pdYr);
%! [~, x] = partialDependence (Mdl, 1);
%! assert_equal (numel (x), 100);
%! assert_equal ([min(x), max(x)], [1, 3]);

%!test  # the query points come from the observations handed over
%! Mdl = fitrtree (pdX, pdYr);
%! wide = [repmat([10;15;20], 4, 1), 4 * ones(12, 1)];
%! [~, x] = partialDependence (Mdl, 1, wide);
%! assert_equal ([min(x), max(x)], [10, 20]);

%!test  # a second output is empty where only one predictor was named
%! Mdl = fitrsvm (pdX, pdYr);
%! [~, ~, y] = partialDependence (Mdl, 1, 'QueryPoints', pdQ);
%! assert_equal (isempty (y), true);

%!test  # sampling draws from the observations and the answer stays finite
%! Mdl = fitrsvm (pdX, pdYr);
%! pd = partialDependence (Mdl, 1, 'QueryPoints', pdQ, ...
%!                         'NumObservationsToSample', 6);
%! assert_equal (size (pd), [1, 3]);
%! assert_equal (all (isfinite (pd)), true);

%!test  # asking for more observations than there are takes all of them
%! Mdl = fitrsvm (pdX, pdYr);
%! a = partialDependence (Mdl, 1, 'QueryPoints', pdQ);
%! b = partialDependence (Mdl, 1, 'QueryPoints', pdQ, ...
%!                        'NumObservationsToSample', 500);
%! assert_equal (a, b);

## Naming the predictors
%!test  # a predictor may be named rather than indexed
%! Mdl = fitrsvm (pdX, pdYr, 'PredictorNames', {'a', 'b'});
%! assert_equal (partialDependence (Mdl, 'a', 'QueryPoints', pdQ), ...
%!               partialDependence (Mdl, 1, 'QueryPoints', pdQ));

%!test  # and two may be, in a cellstr
%! Mdl = fitrsvm (pdX, pdYr, 'PredictorNames', {'a', 'b'});
%! byname = partialDependence (Mdl, {'a', 'b'}, 'QueryPoints', {pdQ, [2;5]});
%! byidx = partialDependence (Mdl, [1, 2], 'QueryPoints', {pdQ, [2;5]});
%! assert_equal (byname, byidx);

## A function handle in place of a model
%!test
%! f = @(Z) 2 * Z(:,1) + Z(:,2);
%! assert_equal (partialDependence (f, 1, pdX, 'QueryPoints', pdQ), ...
%!               [2 + 4.5, 4 + 4.5, 6 + 4.5], 1e-10);

%!test  # and its columns may be picked
%! f = @(Z) [Z(:,1), 2 * Z(:,1)];
%! pd = partialDependence (f, 1, pdX, 'QueryPoints', pdQ, ...
%!                         'OutputColumns', 2);
%! assert_equal (pd, [2, 4, 6], 1e-10);

## The two ways of writing the call answer alike
%!test
%! Mdl = fitrsvm (pdX, pdYr);
%! [a, ax, ay] = partialDependence (Mdl, 1, 'QueryPoints', pdQ);
%! [b, bx, by] = Mdl.partialDependence (1, 'QueryPoints', pdQ);
%! assert_equal (a, b);
%! assert_equal (ax, bx);
%! assert_equal (ay, by);

## Input validation
%!error<partialDependence: too few input arguments.> ...
%! partialDependence (1)

%!error<partialDependence: MDL must be a fitted model that predicts, or a function handle.> ...
%! partialDependence (42, 1)

%!error<partialDependence: LABELS is required for a classification model.> ...
%! partialDependence (fitctree (pdX, pdYc), 1)

%!error<partialDependence: LABELS applies only to a classification model.> ...
%! partialDependence (fitrsvm (pdX, pdYr), 1, {'a'})

%!error<partialDependence: LABELS names a class the model was not fitted with.> ...
%! partialDependence (fitctree (pdX, pdY3), 1, {'z'}, 'QueryPoints', pdQ)

%!error<partialDependence: DATA is required for a model that does not keep the observations it was fitted on.> ...
%! partialDependence (compact (fitrsvm (pdX, pdYr)), 1)

%!error<partialDependence: partial dependence of a generalized additive model is taken over the observations it was fitted on, which a compact one does not keep; use the model it was compacted from.> ...
%! partialDependence (compact (fitrgam (pdX, pdYr)), 1)

%!error<partialDependence: DATA must be a real numeric matrix.> ...
%! partialDependence (fitctree (pdX, pdYc), 1, 1, {1, 2})

%!error<partialDependence: VARS must index the predictors of the model.> ...
%! partialDependence (fitrsvm (pdX, pdYr), 5)

%!error<partialDependence: VARS does not name a predictor: 'z'> ...
%! partialDependence (fitrsvm (pdX, pdYr), 'z')

%!error<partialDependence: VARS must name two different predictors.> ...
%! partialDependence (fitrsvm (pdX, pdYr), [1, 1])

%!error<partialDependence: VARS must name one or two predictors.> ...
%! partialDependence (fitrsvm (pdX, pdYr), [1, 2, 1])

%!error<partialDependence: 'UseParallel' is not implemented.> ...
%! partialDependence (fitrsvm (pdX, pdYr), 1, 'UseParallel', true)

%!error<partialDependence: 'PredictionForMissingValue' is not implemented.> ...
%! partialDependence (fitrsvm (pdX, pdYr), 1, 'PredictionForMissingValue', 0)

%!error<partialDependence: 'IncludeInteractions' applies only to a generalized additive model.> ...
%! partialDependence (fitrsvm (pdX, pdYr), 1, 'IncludeInteractions', true)

%!error<partialDependence: 'IncludeIntercept' applies only to a generalized additive model.> ...
%! partialDependence (fitrsvm (pdX, pdYr), 1, 'IncludeIntercept', false)

%!error<partialDependence: 'CategoricalPredictors' applies only to a function handle, a model carrying its own.> ...
%! partialDependence (fitrsvm (pdX, pdYr), 1, 'CategoricalPredictors', 1)

%!error<partialDependence: 'OutputColumns' applies only to a function handle.> ...
%! partialDependence (fitrsvm (pdX, pdYr), 1, 'OutputColumns', 1)

%!error<partialDependence: 'QueryPoints' must be a vector for a single variable.> ...
%! partialDependence (fitrsvm (pdX, pdYr), 1, 'QueryPoints', ones (3, 3))

%!error<partialDependence: 'QueryPoints' must have one column per variable, or be a cell holding a vector for each.> ...
%! partialDependence (fitrsvm (pdX, pdYr), [1, 2], 'QueryPoints', ones (3, 3))

%!error<partialDependence: 'NumObservationsToSample' must be a positive integer.> ...
%! partialDependence (fitrsvm (pdX, pdYr), 1, 'NumObservationsToSample', 2.5)

%!error<partialDependence: DATA is required when the model is a function handle.> ...
%! partialDependence (@(Z) Z(:,1), 1)
