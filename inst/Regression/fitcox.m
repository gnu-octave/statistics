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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{mdl} =} fitcox (@var{X}, @var{T})
## @deftypefnx {statistics} {@var{mdl} =} fitcox (@var{tbl}, @var{respvar})
## @deftypefnx {statistics} {@var{mdl} =} fitcox (@dots{}, @var{Name}, @var{Value})
##
## Fit a Cox proportional hazards regression model.
##
## @code{@var{mdl} = fitcox (@var{X}, @var{T})} fits the Cox proportional
## hazards model
##
## @tex
## $$ h(x_i, t) = h_0(t)\exp\left(\sum_{j=1}^{p} x_{ij} b_j\right) $$
## @end tex
## @ifnottex
## @math{h(x_i, t) = h_0(t) exp (x_i' b)}
## @end ifnottex
##
## to the @math{n}-by-@math{p} numeric matrix of predictors @var{X} and the
## @math{n}-by-1 vector of event times @var{T}, and returns a @code{CoxModel}
## object.  @var{T} may instead be an @math{n}-by-2 matrix whose rows give a
## @math{(start, stop]} interval of exposure, the counting process form, in
## which an observation joins the risk set only after its start time.
##
## @math{h_0(t)} is the baseline hazard, which is left unspecified: the
## coefficients are estimated by maximizing the Cox partial likelihood, which
## does not involve it.  @strong{@var{X} must not contain a column of ones},
## the model having no constant term, since any constant is absorbed into that
## baseline.
##
## @code{@var{mdl} = fitcox (@var{tbl}, @var{respvar})} takes the data from the
## table @var{tbl}, using the variable named @var{respvar} as the response and
## every other variable as a predictor.  A @code{categorical} variable is
## encoded as indicator columns, one per level bar the first, which the baseline
## hazard carries.
##
## The following @var{Name}/@var{Value} pairs are accepted:
##
## @multitable @columnfractions 0.25 0.75
## @headitem Name @tab Value
## @item @qcode{"Baseline"} @tab The @var{X} values at which the baseline hazard
## is computed, either a scalar or a 1-by-@math{p} vector.  The default is the
## mean of each numeric predictor and zero for each indicator column of a
## categorical predictor, taken within each stratum.  Pass @qcode{0} for a
## hazard relative to the origin.  The coefficients do not depend on this
## choice.
## @item @qcode{"Beta"} @tab The starting value of the iteration, a vector of
## length @math{p}.  The default is @code{0.01 ./ std (@var{X})}.
## @item @qcode{"CategoricalPredictors"} @tab The predictors to treat as
## categorical, given as column indices, a logical vector, or a cell array of
## predictor names.  Table variables of class @code{categorical} are detected
## without this argument.
## @item @qcode{"Censoring"} @tab A logical or 0/1 vector of length @math{n},
## where 1 marks an observation right-censored at its recorded time.  The
## default is a vector of zeros, so every observation is a recorded event.
## @item @qcode{"Frequency"} @tab A vector of length @math{n} of non-negative
## values giving the number of observations each row represents, or a weight.
## The default is a vector of ones.
## @item @qcode{"OptimizationOptions"} @tab A structure of iteration settings,
## as built by @code{statset ("fitcox")}.  The fields used are
## @qcode{"MaxIter"}, @qcode{"TolX"} and @qcode{"Display"}.
## @item @qcode{"PredictorNames"} @tab A cell array of @math{p} predictor names.
## The default is @qcode{"X1"}, @qcode{"X2"}, and so on, or the table variable
## names.
## @item @qcode{"Stratification"} @tab A vector of length @math{n} of stratum
## labels.  Each stratum carries its own baseline hazard and its own risk sets,
## while the coefficients are shared across all of them.
## @item @qcode{"TieBreakMethod"} @tab The method of handling tied event times,
## either @qcode{"breslow"} (default) or @qcode{"efron"}.
## @end multitable
##
## @code{fitcox} is the object interface to @code{coxphfit}, which fits the same
## model and returns the estimates as plain arrays.  The two agree exactly; the
## object additionally reports the proportional hazards assumption tests and
## carries the @code{survival}, @code{hazardratio}, @code{coefci},
## @code{linhyptest} and @code{plotSurvival} methods.
##
## @seealso{CoxModel, coxphfit, ecdf, statset}
## @end deftypefn

function mdl = fitcox (X, T, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  mdl = CoxModel (X, T, varargin{:});

endfunction

%!demo
%! ## Fit a Cox proportional hazards model to right-censored survival times
%! X = [2 0; 5 1; 3 0; 8 1; 4 0; 7 1; 6 0; 9 1; 5 0; 10 1];
%! T = [4; 6; 8; 11; 13; 16; 18; 21; 25; 30];
%! C = [0; 0; 1; 0; 0; 1; 0; 0; 1; 0];
%! mdl = fitcox (X, T, 'Censoring', C)

%!demo
%! ## The hazard of each observation relative to the average one
%! X = [2 0; 5 1; 3 0; 8 1; 4 0; 7 1; 6 0; 9 1; 5 0; 10 1];
%! T = [4; 6; 8; 11; 13; 16; 18; 21; 25; 30];
%! mdl = fitcox (X, T);
%! hazardratio (mdl, X)

%!demo
%! ## Each stratum carries its own baseline hazard
%! X = [2 0; 5 1; 3 0; 8 1; 4 0; 7 1; 6 0; 9 1; 5 0; 10 1];
%! T = [4; 6; 8; 11; 13; 16; 18; 21; 25; 30];
%! S = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];
%! mdl = fitcox (X, T, 'Stratification', S);
%! mdl.Baseline

%!shared X, T, C, S
%! X = [2 0; 5 1; 3 0; 8 1; 4 0; 7 1; 6 0; 9 1; 5 0; 10 1];
%! T = [4; 6; 8; 11; 13; 16; 18; 21; 25; 30];
%! C = [0; 0; 1; 0; 0; 1; 0; 0; 1; 0];
%! S = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];

## fitcox returns a CoxModel object
%!test
%! mdl = fitcox (X, T);
%! assert_equal (class (mdl), 'CoxModel');
%! assert_equal (mdl.Coefficients.Beta, ...
%!               [-1.3886093196382836; 4.3814437183613322], 1e-8);

## It fits the same model as coxphfit
%!test
%! mdl = fitcox (X, T, 'Censoring', C);
%! [b, logl] = coxphfit (X, T, 'Censoring', C);
%! assert_equal (mdl.Coefficients.Beta, b, 1e-12);
%! assert_equal (mdl.LogLikelihood, logl, 1e-12);

## The name-value pairs are fitcox's, not coxphfit's
%!test
%! mdl = fitcox (X, T, 'Censoring', C, 'TieBreakMethod', 'efron');
%! b = coxphfit (X, T, 'Censoring', C, 'Ties', 'efron');
%! assert_equal (mdl.Coefficients.Beta, b, 1e-12);

%!test
%! mdl = fitcox (X, T, 'Censoring', C, 'Stratification', S);
%! b = coxphfit (X, T, 'Censoring', C, 'Strata', S);
%! assert_equal (mdl.Coefficients.Beta, b, 1e-12);

%!test
%! mdl = fitcox (X, T, 'Beta', [0.1; -0.1]);
%! b = coxphfit (X, T, 'B0', [0.1; -0.1]);
%! assert_equal (mdl.Coefficients.Beta, b, 1e-12);

## A table names its own variables
%!test
%! tbl = table (X(:,1), X(:,2), T, 'VariableNames', {'age', 'trt', 'time'});
%! mdl = fitcox (tbl, 'time');
%! assert_equal (mdl.PredictorNames, {'age', 'trt'});
%! assert_equal (mdl.ResponseName, 'time');
%! assert_equal (char (mdl.Formula), 'time ~ age + trt');
%! assert_equal (mdl.Coefficients.Beta, ...
%!               [-1.3886093196382836; 4.3814437183613322], 1e-8);

## A categorical predictor is reference coded against its first level
%!test
%! g = categorical ({'a'; 'b'; 'a'; 'b'; 'a'; 'b'; 'a'; 'b'; 'a'; 'b'});
%! tbl = table (X(:,1), g, T, 'VariableNames', {'age', 'grp', 'time'});
%! mdl = fitcox (tbl, 'time');
%! assert_equal (mdl.PredictorNames, {'age', 'grp'});
%! assert_equal (mdl.Coefficients.Properties.RowNames, {'age'; 'grp_b'});
%! assert_equal (mdl.Coefficients.Beta, ...
%!               [-1.3886093196382836; 4.3814437183613322], 1e-8);
%! assert_equal (mdl.Baseline, [5.9, 0], 1e-12);
%! assert_equal (mdl.VariableInfo.IsCategorical, [false; true; false]);

## Errors
%!error<Invalid call to fitcox> fitcox (1)
%!error<CoxModel: 'time' is not a variable of the data table.> ...
%! fitcox (table ([1; 2], [3; 4]), 'time')
