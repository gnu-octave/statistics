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
## @deftypefn  {statistics} {@var{mdl} =} fitglm (@var{X}, @var{y})
## @deftypefnx {statistics} {@var{mdl} =} fitglm (@var{X}, @var{y}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} fitglm (@var{tbl})
## @deftypefnx {statistics} {@var{mdl} =} fitglm (@var{tbl}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} fitglm (@dots{}, @var{Name}, @var{Value})
##
## Fit a generalized linear regression model.
##
## @code{@var{mdl} = fitglm (@var{X}, @var{y})} fits a generalized linear model
## of the response vector @var{y} on the columns of the @math{n}-by-@math{p}
## numeric predictor matrix @var{X}, and returns a @code{GeneralizedLinearModel}
## object.  @code{@var{mdl} = fitglm (@var{tbl})} instead takes the predictors
## and response from the table @var{tbl} (the last column is the response unless
## overridden).  By default the response is @qcode{'normal'} with an identity
## link, an intercept is included, and the model is additive in the predictors.
##
## @var{modelspec} selects the model terms.  It is either a Wilkinson formula
## string (e.g.@: @qcode{'y ~ x1 + x2*x3'}), a keyword (@qcode{'constant'},
## @qcode{'linear'}, @qcode{'interactions'}, @qcode{'purequadratic'},
## @qcode{'quadratic'}, or @qcode{'full'}), or a terms matrix.
##
## The following @var{Name}/@var{Value} pairs are accepted:
##
## @multitable @columnfractions 0.2 0.75
## @headitem Name @tab Value
## @item @qcode{'Distribution'} @tab the response distribution:
## @qcode{'normal'} (default), @qcode{'binomial'}, @qcode{'poisson'},
## @qcode{'gamma'}, or @qcode{'inverse gaussian'}.
## @item @qcode{'Link'} @tab the link function.  Defaults to the canonical link
## of the distribution; accepts any link name understood by @code{glmfit} or a
## numeric exponent for a power link.
## @item @qcode{'Weights'} @tab a vector of nonnegative observation weights.
## @item @qcode{'Offset'} @tab a vector added as a fixed term to the linear
## predictor.
## @item @qcode{'BinomialSize'} @tab for the @qcode{'binomial'} distribution,
## the number of trials (a scalar or a per-observation vector); @var{y} holds
## the proportion of successes.
## @item @qcode{'Intercept'} @tab a logical value (default @qcode{true}) whether
## to include an intercept term.
## @item @qcode{'DispersionFlag'} @tab a logical value forcing the dispersion
## parameter to be estimated (@qcode{true}) or held at 1 (@qcode{false}).
## @item @qcode{'CategoricalVars'} @tab predictors to treat as categorical (a
## logical vector, numeric indices, or a cell array of names).
## @item @qcode{'Exclude'} @tab observations to exclude from the fit (a logical
## vector or numeric indices).
## @item @qcode{'VarNames'} @tab a cell array of @math{p + 1} variable names
## (predictors followed by the response) for numeric @var{X}.
## @item @qcode{'PredictorVars'}, @qcode{'ResponseVar'} @tab for table input,
## the predictor and response variable names.
## @end multitable
##
## @seealso{GeneralizedLinearModel, fitlm, glmfit, glmval, lassoglm}
## @end deftypefn

function mdl = fitglm (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  arg1 = varargin{1};
  if (istable (arg1))
    [modelspec, nv] = split_modelspec (varargin(2:end));
    mdl = GeneralizedLinearModel (arg1, [], modelspec, nv{:});
  else
    if (nargin < 2)
      print_usage ();
    endif
    [modelspec, nv] = split_modelspec (varargin(3:end));
    mdl = GeneralizedLinearModel (arg1, varargin{2}, modelspec, nv{:});
  endif

endfunction

## Split an optional leading model specification from the Name/Value pairs.
function [modelspec, nv] = split_modelspec (rest)
  modelspec = 'linear';
  nv        = rest;
  if (! isempty (rest))
    a = rest{1};
    if ((ischar (a) && ! is_param_name (a)) || isnumeric (a))
      modelspec = a;
      nv        = rest(2:end);
    endif
  endif
endfunction

## True if S names one of fitglm's Name/Value parameters.
function tf = is_param_name (s)
  tf = ischar (s) && any (strcmpi (s, {'Distribution', 'Link', 'Weights', ...
       'Offset', 'BinomialSize', 'Intercept', 'DispersionFlag', ...
       'CategoricalVars', 'Exclude', 'VarNames', 'PredictorVars', ...
       'ResponseVar'}));
endfunction

%!demo
%! ## Poisson regression of counts on two predictors.
%! X = [0.1, 1.2; 0.4, 0.7; 1.1, 0.2; 1.5, 1.9; 0.3, 0.5; 1.8, 1.1; 0.9, 0.3];
%! y = [1; 0; 2; 3; 1; 4; 2];
%! mdl = fitglm (X, y, 'Distribution', 'poisson')

%!demo
%! ## Logistic regression with an interaction, specified by a formula.
%! X = [0.1, 1.2; 0.4, 0.7; 1.1, 0.2; 1.5, 1.9; 0.3, 0.5; 1.8, 1.1; 0.9, 0.3];
%! y = [0; 0; 1; 1; 0; 1; 1];
%! tbl = array2table ([X, y], 'VariableNames', {'x1', 'x2', 'y'});
%! mdl = fitglm (tbl, 'y ~ x1 + x2 + x1:x2', 'Distribution', 'binomial')

## Test input validation
%!error<Invalid call> fitglm ()
%!error<GeneralizedLinearModel: X must be a real matrix.> ...
%! fitglm ("a", [1;2])
%!error<GeneralizedLinearModel: unknown distribution 'wibble'.> ...
%! fitglm ([1, 2; 3, 4], [1; 0], 'Distribution', 'wibble')
%!error<GeneralizedLinearModel: unknown parameter name 'foo'.> ...
%! fitglm ([1, 2; 3, 4], [1; 0], 'linear', 'foo', 1)
