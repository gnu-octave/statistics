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
## @deftypefn  {statistics} {@var{mdl} =} fitmnr (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{mdl} =} fitmnr (@var{Tbl}, @var{Y})
## @deftypefnx {statistics} {@var{mdl} =} fitmnr (@var{Tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{mdl} =} fitmnr (@var{Tbl}, @var{Formula})
## @deftypefnx {statistics} {@var{mdl} =} fitmnr (@dots{}, @var{name}, @var{value})
##
## Fit a multinomial regression model.
##
## @code{@var{mdl} = fitmnr (@var{X}, @var{Y})} fits a multinomial regression
## of the response @var{Y} on the columns of the @math{N}-by-@math{P} numeric
## predictor matrix @var{X}, and returns a @code{MultinomialRegression} object.
## An intercept is always fitted, so @var{X} holds no column of ones.
##
## @var{Y} is an @math{N}-by-1 numeric, logical, string or categorical vector,
## a cell array of character vectors, or a character matrix, each element or
## row the category of one observation; or an @math{N}-by-@math{K} matrix of
## counts, @code{@var{Y}(i,j)} observations of category @math{j} at the
## predictor values @code{@var{X}(i,:)}, as @code{mnrfit} documents.
##
## @code{@var{mdl} = fitmnr (@var{Tbl}, @var{Y})} takes the predictors from
## the variables of the table @var{Tbl}, and the response from @var{Y}.
## @code{@var{mdl} = fitmnr (@var{Tbl}, @var{ResponseVarName})} takes the
## response from the variable of @var{Tbl} that @var{ResponseVarName} names, a
## character vector or a string scalar, and the predictors from the others.
## @code{@var{mdl} = fitmnr (@var{Tbl}, @var{Formula})} takes both from the
## Wilkinson formula @var{Formula}, such as @qcode{'y ~ x1 + x2*x3'}.  A
## numeric variable of a table is a continuous predictor and any other is
## categorical, coded as one indicator per level bar the first.
##
## The following @var{name}, @var{value} pairs are accepted:
##
## @multitable @columnfractions 0.3 0.65
## @headitem Name @tab Value
## @item @qcode{'ModelType'} @tab @qcode{'nominal'} (default), a model
## contrasting each category with the last; @qcode{'ordinal'}, a model of the
## cumulative probabilities of ordered categories; or @qcode{'hierarchical'},
## a model of each category against the ones after it.
## @item @qcode{'Link'} @tab The link function of an ordinal or hierarchical
## model: @qcode{'logit'} (default), @qcode{'probit'}, @qcode{'comploglog'}
## or @qcode{'loglog'}.  A nominal model always uses the logit link.
## @item @qcode{'IncludeClassInteractions'} @tab @qcode{true} to give each
## category its own slopes, @qcode{false} for slopes shared by every
## category.  The default is @qcode{true} for nominal and hierarchical models
## and @qcode{false} for ordinal ones.  A nominal model without interactions
## fits its intercepts alone.
## @item @qcode{'EstimateDispersion'} @tab @qcode{true} to estimate the
## dispersion and scale the standard errors by it, @qcode{false} (default) for
## the theoretical dispersion of 1.
## @item @qcode{'Weights'} @tab An @math{N}-by-1 vector of nonnegative
## observation weights, each multiplying the counts of its observation; ones
## by default.
## @item @qcode{'CategoricalPredictors'} @tab The predictors to treat as
## categorical: a vector of indices, a logical vector with one element per
## predictor, a cell array or string array of names, or @qcode{'all'}.
## @item @qcode{'PredictorNames'} @tab The names of the predictors: for a
## matrix @var{X}, one per column, @qcode{x1}, @qcode{x2}, @dots{} by default;
## for a table, the variables to use as predictors.  A formula names its own.
## @item @qcode{'ResponseName'} @tab The name of the response, @qcode{'y'} by
## default.  A response variable name or a formula names its own.
## @item @qcode{'IterationLimit'} @tab The largest number of iterations of the
## fit, 100 by default.
## @item @qcode{'Tolerance'} @tab The convergence tolerance of the fit, the
## largest change of a coefficient relative to its value, @code{1e-6} by
## default.
## @end multitable
##
## A row holding a missing value of a predictor it uses, of the response, or
## of the weights is left out of the fit.
##
## @seealso{MultinomialRegression, mnrfit, mnrval, fitglm}
## @end deftypefn

function mdl = fitmnr (varargin)

  if (nargin < 2)
    error ("fitmnr: too few input arguments.");
  endif
  mdl = MultinomialRegression (varargin{:});

endfunction

%!demo
%! ## A nominal model of the iris species on the length of the sepal, and
%! ## the probabilities it gives three new flowers
%! load fisheriris
%! mdl = fitmnr (meas(:,1), species, 'PredictorNames', {'SepalLength'})
%! [label, probs] = predict (mdl, [5; 6; 7])

%!demo
%! ## An ordinal model, the species ordered by size, and the cumulative
%! ## probabilities it gives three new flowers
%! load fisheriris
%! mdl = fitmnr (meas(:,1), species, 'ModelType', 'ordinal')
%! [~, cumprob] = predict (mdl, [5; 6; 7], 'ProbabilityType', 'cumulative')

## Test output
%!test
%! load fisheriris
%! mdl = fitmnr (meas(:,1), species, 'ModelType', 'ordinal');
%! assert_equal (class (mdl), 'MultinomialRegression');
%!test
%! load fisheriris
%! mdl = fitmnr (meas(:,1), species, 'ModelType', 'ordinal');
%! B = mnrfit (meas(:,1), species, 'model', 'ordinal');
%! assert_equal (mdl.Coefficients.Value, B, 1e-12);

## Test input validation
%!error<fitmnr: too few input arguments.> fitmnr (1)
