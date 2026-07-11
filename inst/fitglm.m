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
## @deftypefnx {statistics} {@var{mdl} =} fitglm (@var{X}, @var{y}, @var{Name}, @var{Value})
##
## Fit a generalized linear regression model.
##
## @code{@var{mdl} = fitglm (@var{X}, @var{y})} fits a generalized linear model
## of the response vector @var{y} on the columns of the @math{n}-by-@math{p}
## numeric predictor matrix @var{X}, and returns a @code{GeneralizedLinearModel}
## object.  By default the response is @qcode{'normal'} with an identity link and
## an intercept is included.
##
## @code{@var{mdl} = fitglm (@var{X}, @var{y}, @var{Name}, @var{Value})} accepts
## the following @var{Name}/@var{Value} pairs:
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
## the number of trials (a scalar or a per-observation vector); @var{y} holds the
## proportion of successes.
## @item @qcode{'Intercept'} @tab a logical value (default @qcode{true}) whether
## to include an intercept term.
## @item @qcode{'DispersionFlag'} @tab a logical value forcing the dispersion
## parameter to be estimated (@qcode{true}) or held at 1 (@qcode{false}).  The
## default depends on the distribution.
## @item @qcode{'VarNames'} @tab a cell array of @math{p + 1} variable names
## (predictors followed by the response).
## @end multitable
##
## @seealso{GeneralizedLinearModel, fitlm, glmfit, glmval, lassoglm}
## @end deftypefn

function mdl = fitglm (X, y, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  ## An optional positional model specification precedes the Name/Value pairs.
  ## Full Wilkinson formulae and table inputs are handled in a later phase; for
  ## now only the additive 'linear' specification is supported for numeric X.
  args = varargin;
  if (! isempty (args) && (ischar (args{1}) || isnumeric (args{1})) ...
      && ! is_param_name (args{1}))
    modelspec = args{1};
    args(1) = [];
    if (! (ischar (modelspec) && strcmpi (modelspec, 'linear')))
      error (strcat ("fitglm: only the 'linear' model specification is", ...
                     " currently supported for numeric X."));
    endif
  endif

  mdl = GeneralizedLinearModel (X, y, args{:});

endfunction

## True if S names one of fitglm's Name/Value parameters.
function tf = is_param_name (s)
  tf = ischar (s) && any (strcmpi (s, {'Distribution', 'Link', 'Weights', ...
       'Offset', 'BinomialSize', 'Intercept', 'DispersionFlag', 'VarNames'}));
endfunction

%!demo
%! ## Poisson regression of counts on two predictors.
%! X = [0.1, 1.2; 0.4, 0.7; 1.1, 0.2; 1.5, 1.9; 0.3, 0.5; 1.8, 1.1; 0.9, 0.3];
%! y = [1; 0; 2; 3; 1; 4; 2];
%! mdl = fitglm (X, y, 'Distribution', 'poisson')

## Test input validation
%!error<Invalid call> fitglm (1)
%!error<GeneralizedLinearModel: X must be a real matrix.> fitglm ("a", [1;2])
%!error<GeneralizedLinearModel: Y must be a real vector.> ...
%! fitglm ([1, 2; 3, 4], "a")
%!error<GeneralizedLinearModel: unknown distribution 'wibble'.> ...
%! fitglm ([1, 2; 3, 4], [1; 0], 'Distribution', 'wibble')
%!error<fitglm: only the 'linear' model specification> ...
%! fitglm ([1, 2; 3, 4], [1; 0], 'quadratic')
