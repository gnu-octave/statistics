## Copyright (C) 2022 Andrew Penn <A.C.Penn@sussex.ac.uk>
## Copyright (C) 2026 Avanish Salunke <avanishsalunke16@gmail.com>
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
## @deftypefn  {statistics} {@var{mdl} =} fitlm (@var{X}, @var{y})
## @deftypefnx {statistics} {@var{mdl} =} fitlm (@var{tbl})
## @deftypefnx {statistics} {@var{mdl} =} fitlm (@var{tbl}, @var{ResponseVarName})
## @deftypefnx {statistics} {@var{mdl} =} fitlm (@var{tbl}, @var{y})
## @deftypefnx {statistics} {@var{mdl} =} fitlm (@dots{}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} fitlm (@dots{}, @var{Name}, @var{Value}, @dots{})
##
## Fit a linear regression model to data and return a @code{LinearModel}
## object.
##
## The returned object stores the fitted coefficients, their standard errors,
## t-statistics, and p-values, summary statistics of the fit (@math{R^2},
## RMSE, F-statistic, etc.), and the residuals and diagnostics of the fit, and
## exposes methods such as @code{predict}, @code{plotResiduals},
## @code{coefTest}, @code{addTerms}, and @code{removeTerms} for further
## analysis of the fitted model.
##
## @subheading Basic Syntax
##
## @code{@var{mdl} = fitlm (@var{X}, @var{y})} fits a linear regression model
## of the response @var{y} to the predictor data @var{X}.  Unless removed via
## the @qcode{'Intercept'} option, the fitted model contains a constant
## (intercept) term and one linear term for every column of @var{X}.
##
## @itemize
## @item
## @var{X} is an @math{NxP} numeric or logical matrix of predictor data, where
## rows correspond to observations and columns correspond to variables.  By
## default, the predictors are named @qcode{'x1'}, @qcode{'x2'}, @dots{},
## @qcode{'xP'}.
## @item
## @var{X} can also be a categorical vector of length @math{N}, representing a
## single categorical predictor.  In this case @var{y} must be supplied as the
## next argument, and the predictor is named @qcode{'x1'} by default.
## @item
## @var{y} is an @math{Nx1} numeric or logical vector of response values, and
## must have the same number of observations (rows) as @var{X}.  By default,
## the response is named @qcode{'y'}.
## @end itemize
##
## @code{@var{mdl} = fitlm (@var{tbl})} fits a linear regression model using
## the variables contained in the table (or dataset) @var{tbl}.  By default,
## the last variable in @var{tbl} is used as the response and all other
## variables are used as predictors.  Variables that are @code{categorical}
## arrays, cell arrays of character vectors, or logical arrays are
## automatically treated as categorical predictors.
##
## @code{@var{mdl} = fitlm (@var{tbl}, @var{ResponseVarName})} fits a model
## using the variable named @var{ResponseVarName} in @var{tbl} as the
## response, and all remaining variables in @var{tbl} as predictors.
##
## @code{@var{mdl} = fitlm (@var{tbl}, @var{y})} fits a model using the
## variables in @var{tbl} as predictors and the external numeric vector
## @var{y} as the response.  @var{y} must have @code{height (@var{tbl})}
## elements.
##
## @subheading Model Specification
##
## @code{@var{mdl} = fitlm (@dots{}, @var{modelspec})} additionally specifies
## the terms of the model to fit, using any of the input combinations shown
## above.  @var{modelspec} can be any of the following.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Value} @tab @tab @var{Description}
##
## @item @qcode{'constant'} @tab @tab Model contains only an intercept term.
##
## @item @qcode{'linear'} @tab @tab Model contains an intercept and one term
## for each predictor variable.  This is the default when @var{modelspec} is
## not specified.
##
## @item @qcode{'interactions'} @tab @tab Model contains an intercept, all
## linear terms, and all pairwise products of distinct predictor variables
## (no squared terms).
##
## @item @qcode{'purequadratic'} @tab @tab Model contains an intercept, all
## linear terms, and all squared terms.
##
## @item @qcode{'quadratic'} @tab @tab Model contains an intercept, all linear
## terms, all pairwise products of distinct predictor variables, and all
## squared terms.
##
## @item @qcode{'full'} @tab @tab Model contains an intercept and all terms up
## to and including the full @math{P}-way interaction of the predictor
## variables, i.e. every combination of one or more distinct predictors.
##
## @item terms matrix @tab @tab A @math{TxP} or @math{Tx(P+1)} numeric matrix,
## where @math{T} is the number of terms and @math{P} is the number of
## predictor variables.  Each row represents one term, and the value in
## column @math{j} is the exponent to which predictor @math{j} is raised in
## that term; a row of all zeros represents the intercept.  If a
## @math{Tx(P+1)} matrix is supplied, its last column (representing the
## response variable) must be all zeros.
##
## @item Wilkinson formula @tab @tab A character vector of the form
## @qcode{'y ~ terms'} describing the response and predictor terms using
## Wilkinson notation.  The variable name to the left of @qcode{'~'} is used
## as the response, overriding any response implied elsewhere in the call.
## @end multitable
##
## When @var{modelspec} is given as a Wilkinson formula, the following
## operators may be used on its right-hand side to build up @code{terms}:
##
## @multitable @columnfractions 0.12 0.38 0.5
## @headitem Operator @tab Meaning @tab Example
## @item @code{+} @tab add a term @tab @qcode{'x1 + x2'} adds @code{x1} and
## @code{x2} as separate terms
## @item @code{-} @tab remove a term @tab @qcode{'x1*x2 - x1:x2'} removes the
## interaction, leaving only @code{x1} and @code{x2}
## @item @code{*} @tab cross two terms @tab @qcode{'x1*x2'} expands to
## @code{x1}, @code{x2}, @code{x1:x2}
## @item @code{:} @tab interaction only @tab @qcode{'x1:x2'} adds only the
## interaction term between @code{x1} and @code{x2}
## @item @code{^} @tab power / crossing limit @tab @qcode{'x^2'} adds
## @code{x} and @code{x^2}; @qcode{'(x1+x2)^2'} expands to @code{x1},
## @code{x2}, @code{x1:x2}
## @item @code{-1} @tab remove intercept @tab @qcode{'x1 + x2 - 1'} fits the
## model without a constant term
## @end multitable
##
## A formula includes an intercept term by default; append @qcode{'- 1'} to
## the formula to omit it.  For a categorical predictor, @code{fitlm}
## generates the necessary indicator (dummy) variables automatically from the
## formula, so a formula does not need to be changed when the underlying
## design matrix changes.
##
## @subheading Options
##
## @code{@var{mdl} = fitlm (@dots{}, @var{Name}, @var{Value}, @dots{})}
## specifies additional options using one or more @qcode{Name-Value} pair
## arguments, which may be combined with @var{modelspec} or used on their own.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{'Intercept'} @tab @tab A logical scalar indicating whether to
## include a constant (intercept) term in the model.  Default is @qcode{true}.
## This option only applies when @var{modelspec} is a character vector model
## name (or omitted); it is ignored when @var{modelspec} is a terms matrix or
## a Wilkinson formula, where the intercept is instead controlled by the
## matrix/formula itself.
##
## @item @qcode{'Weights'} @tab @tab A numeric vector of nonnegative
## observation weights, with one element per observation, used to fit a
## weighted least squares model.  Default is a vector of ones, i.e. an
## unweighted ordinary least squares fit.
##
## @item @qcode{'Exclude'} @tab @tab A numeric or logical vector specifying
## observations to exclude from the fit, given as row indices into the
## original data or as a logical mask the same length as the number of
## observations.  Excluded observations, together with any observation that
## contains a missing (@qcode{NaN}) value in a predictor or the response, are
## recorded in the @code{ObservationInfo} property of the fitted model but do
## not contribute to the fitted coefficients or summary statistics.
##
## @item @qcode{'CategoricalVars'} @tab @tab Specifies which predictor
## variables are treated as categorical, given as a vector of column indices,
## a logical vector, or a cell array of variable names (only valid for table
## input).  Each categorical predictor with @math{L} distinct categories is
## expanded into @math{L-1} indicator (dummy) variables, using the first
## category (in sorted or original order) as the reference level that is
## omitted from the design matrix.  Variables that are already
## @code{categorical} arrays or cell arrays of character vectors are always
## treated as categorical, regardless of this option.
##
## @item @qcode{'VarNames'} @tab @tab A cell array of character vectors
## naming the predictor and response variables, listed in order with the
## response variable name last, e.g. @code{@{"x1", "x2", "y"@}} for two
## predictors.  Only applies when @var{X} and @var{y} (or a categorical
## vector and @var{y}) are supplied directly, since table variables already
## carry their own names.  By default, predictors are named @qcode{'x1'},
## @qcode{'x2'}, etc. and the response is named @qcode{'y'}.
##
## @item @qcode{'ResponseVar'} @tab @tab A character vector naming the
## response variable, used to override the response variable name that would
## otherwise be inferred (the last table variable, or @qcode{'y'} for matrix
## input).
##
## @item @qcode{'PredictorVars'} @tab @tab A cell array of character vectors
## naming which variables in @var{tbl} to use as predictors.  By default, all
## variables in @var{tbl} other than the response variable are used as
## predictors.
##
## @item @qcode{'RobustOpts'} @tab @tab Selects ordinary least squares or
## robust regression fitting.  This value can be @qcode{'off'} (default,
## ordinary least squares), @qcode{'on'} (robust fitting using the
## @qcode{'bisquare'} weighting function), the name of one of the weighting
## functions below, a function handle for a custom weighting function, or a
## scalar structure with fields @qcode{RobustWgtFun} and @qcode{Tune}
## specifying the weighting function and its tuning constant.  Robust fitting
## uses Iteratively Reweighted Least Squares (IRLS), refitting the model with
## updated observation weights until the coefficients converge.  Supported
## weighting function names: @qcode{'andrews'}, @qcode{'bisquare'},
## @qcode{'cauchy'}, @qcode{'fair'}, @qcode{'huber'}, @qcode{'logistic'},
## @qcode{'ols'}, @qcode{'talwar'}, @qcode{'welsch'}, each with its own default
## tuning constant.
## @end multitable
##
## @subheading Algorithm
##
## @code{fitlm} solves the (weighted) least squares problem by applying a
## pivoted QR decomposition to the design matrix, which remains numerically
## stable even when predictors are collinear; coefficients corresponding to
## columns beyond the numerically detected rank of the design matrix are set
## to zero.  Robust fits refine this ordinary least squares solution using
## IRLS as described above.  Observations with missing values in any variable
## used by the model, or explicitly excluded via @qcode{'Exclude'}, are
## omitted from the fit entirely and flagged in @code{ObservationInfo}, but
## are otherwise not counted as errors.
##
## @var{mdl} is returned as a @code{LinearModel} object.  If
## @qcode{'RobustOpts'} is anything other than @qcode{'off'}, the returned
## model is a robust fit rather than an ordinary least squares fit, and its
## @code{Robust} property is populated accordingly.
##
## @seealso{LinearModel}
## @end deftypefn

function mdl = fitlm (varargin)

  if (nargin < 1)
    error ("fitlm: Not enough input arguments.");
  endif

  ## List of Name-Value keys used to check if the response variable y is missing.
  nv_keys = {'varnames', 'intercept', 'responsevar', 'predictorvars', ...
             'categoricalvars', 'exclude', 'weights', 'robustopts'};
  is_nv = @(s) (ischar (s) || isstring (s)) && ...
               any (strcmpi (char (s), nv_keys));

  arg1 = varargin{1};
  rest = varargin(2:end);

  if (isa (arg1, 'categorical'))
    if (! isvector (arg1))
      error (strcat ("fitlm: Predictor variables must be numeric vectors,", ...
                     " numeric matrices, or categorical vectors."));
    endif
    if (isempty (rest) || is_nv(rest{1}))
      error ("fitlm: Y argument is required unless X is a dataset or table.");
    endif
    y_arg = rest{1};
    if (numel (y_arg) != size (arg1, 1))
      error ("fitlm: Predictor and response variables must have the same length.");
    endif
    if (! isvector (y_arg) || (! isnumeric (y_arg) && ! islogical (y_arg)))
      error ("fitlm: Response variable must be a numeric vector.");
    endif

    pred_name = 'x1';
    resp_name = 'y';
    tail      = rest(2:end);
    keep      = true (1, numel (tail));
    for k = 1:2:numel (tail)-1
      if (ischar (tail{k}) && strcmpi (tail{k}, 'VarNames') && iscell (tail{k+1}))
        vn = tail{k+1};
        if (numel (vn) >= 1); pred_name = vn{1}; endif
        if (numel (vn) >= 2); resp_name = vn{2}; endif
        keep(k:k+1) = false;
      endif
    endfor
    tail = tail(keep);

    tbl = table (arg1(:), double (y_arg(:)), 'VariableNames', {pred_name, resp_name});
    mdl = fitlm (tbl, tail{:});
    return;
  endif

  if (istable (arg1))

    response  = [];
    modelspec = [];
    nv_args   = {};

    if (! isempty (rest))

      ## Even length starting with Name-Value key means all are Name-Value pairs
      if (mod (numel (rest), 2) == 0 && is_nv(rest{1}))
        nv_args = rest;

      else
        arg2       = rest{1};
        after_arg2 = rest(2:end);
        n_rows     = height (arg1);
        n_cols     = width  (arg1);
        col_names  = arg1.Properties.VariableNames;

        if (ischar (arg2) || isstring (arg2))
          s = char (arg2);

          if (any (s == '~'))
            ## Wilkinson formula string
            modelspec = s;
            if (mod (numel (after_arg2), 2) != 0)
              error ("fitlm: Name-Value arguments must be in pairs.");
            endif
            nv_args = after_arg2;

          elseif (any (strcmp (s, col_names)))
            ## Response variable name; ODD/EVEN for remainder
            response = s;
            [modelspec, nv_args] = lm_split_args (after_arg2);

          else
            ## Modelspec keyword or invalid string; LinearModel will validate
            modelspec = s;
            if (mod (numel (after_arg2), 2) != 0)
              error ("fitlm: Name-Value arguments must be in pairs.");
            endif
            nv_args = after_arg2;
          endif

        elseif (isnumeric (arg2) || islogical (arg2))
          [nr2, nc2] = size (arg2);

          if (isempty (arg2))
            error (strcat ("fitlm: The terms matrix must have one column", ...
                           " for each variable in the dataset or table."));

          elseif (nc2 == n_cols)
            ## Column count matches table, so it is a terms matrix
            if (! any (all (double (arg2) == 0, 1)))
              error (strcat ("fitlm: Cannot determine the response", ...
                             " variable from the terms matrix."));
            endif
            modelspec = double (arg2);
            if (mod (numel (after_arg2), 2) != 0)
              error ("fitlm: Name-Value arguments must be in pairs.");
            endif
            nv_args = after_arg2;

          elseif (nc2 == 1 && nr2 == n_rows)
            ## Single column matching table height is an external y vector
            response = double (arg2(:));
            [modelspec, nv_args] = lm_split_args (after_arg2);

          else
            error ("fitlm: Predictor and response variables must have the same length.");
          endif

        else
          error ("fitlm: invalid second argument for table input.");
        endif
      endif
    endif

    mdl = LinearModel (arg1, response, modelspec, nv_args{:});

  elseif ((isnumeric (arg1) || islogical (arg1)) && ismatrix (arg1))

    n = size (arg1, 1);

    ## Ensure the response variable y is provided.
    if (isempty (rest) || is_nv(rest{1}))
      error ("fitlm: Y argument is required unless X is a dataset or table.");
    endif

    arg2 = rest{1};

    ## Ensure the response variable has the correct length.
    if (max ([size(arg2), 0]) != n)
      error ("fitlm: Predictor and response variables must have the same length.");
    endif

    if (! isvector (arg2))
      error ("fitlm: Response variable must be a numeric vector.");
    endif

    if (! isnumeric (arg2) && ! islogical (arg2))
      error ("fitlm: Response variable must be a numeric vector.");
    endif

    y = double (arg2(:));
    [modelspec, nv_args] = lm_split_args (rest(2:end));

    mdl = LinearModel (arg1, y, modelspec, nv_args{:});

  else
    error (strcat ("fitlm: Predictor variables must be numeric vectors,", ...
                   " numeric matrices, or categorical vectors."));
  endif

endfunction


## If count is odd, the first element is the modelspec and the rest are Name-Value pairs.
## If count is even, there is no modelspec and all elements are Name-Value pairs.
function [modelspec, nv_args] = lm_split_args (remaining)
  if (isempty (remaining))
    modelspec = [];
    nv_args   = {};
  elseif (mod (numel (remaining), 2) == 1)
    modelspec = remaining{1};
    nv_args   = remaining(2:end);
  else
    modelspec = [];
    nv_args   = remaining;
  endif
endfunction

%!demo
%! y =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! X = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%!
%! mdl = fitlm (X, y, 'linear', 'CategoricalVars', 1)

%!demo
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! brands = {'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'};
%! popper = {'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; ...
%!           'air', 'air', 'air'; 'air', 'air', 'air'; 'air', 'air', 'air'};
%!
%! T = table (brands(:), popper(:), 'VariableNames', {'brands', 'popper'});
%! mdl = fitlm (T, popcorn(:), 'interactions')

%!test
%! y =  [ 8.706 10.362 11.552  6.941 10.983 10.092  6.421 14.943 15.931 ...
%!        22.968 18.590 16.567 15.944 21.637 14.492 17.965 18.851 22.891 ...
%!        22.028 16.884 17.252 18.325 25.435 19.141 21.238 22.196 18.038 ...
%!        22.628 31.163 26.053 24.419 32.145 28.966 30.207 29.142 33.212 ...
%!        25.694 ]';
%! X = [1 1 1 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5]';
%! fitlm (X, y, 'CategoricalVars', 1);
%! fitlm (X, y, 'constant', 'CategoricalVars', 1);
%! fitlm (X, y, 'linear', 'CategoricalVars', 1);
%! mdl = fitlm (X, y, 'linear', 'CategoricalVars', 1);
%! assert (mdl.Coefficients.Estimate(1), 10, 1e-04);
%! assert (mdl.Coefficients.Estimate(2), 7.99999999999999, 1e-09);
%! assert (mdl.Coefficients.Estimate(3), 8.99999999999999, 1e-09);
%! assert (mdl.Coefficients.Estimate(4), 11.0001428571429, 1e-09);
%! assert (mdl.Coefficients.Estimate(5), 19.0001111111111, 1e-09);
%! assert (mdl.Coefficients.SE(1), 1.01775379540949, 1e-09);
%! assert (mdl.Coefficients.SE(2), 1.64107868458008, 1e-09);
%! assert (mdl.Coefficients.SE(3), 1.43932122062479, 1e-09);
%! assert (mdl.Coefficients.SE(4), 1.48983900477565, 1e-09);
%! assert (mdl.Coefficients.SE(5), 1.3987687997822, 1e-09);
%! assert (mdl.Coefficients.tStat(1), 9.82555903510687, 1e-09);
%! assert (mdl.Coefficients.tStat(2), 4.87484242844031, 1e-09);
%! assert (mdl.Coefficients.tStat(3), 6.25294748040552, 1e-09);
%! assert (mdl.Coefficients.tStat(4), 7.38344399756088, 1e-09);
%! assert (mdl.Coefficients.tStat(5), 13.5834536158296, 1e-09);
%! assert (mdl.Coefficients.pValue(2), 2.85812420217862e-05, 1e-12);
%! assert (mdl.Coefficients.pValue(3), 5.22936741204002e-07, 1e-06);
%! assert (mdl.Coefficients.pValue(4), 2.12794763209106e-08, 1e-07);
%! assert (mdl.Coefficients.pValue(5), 7.82091664406755e-15, 1e-08);

%!test
%! popcorn = [5.5, 4.5, 3.5; 5.5, 4.5, 4.0; 6.0, 4.0, 3.0; ...
%!            6.5, 5.0, 4.0; 7.0, 5.5, 5.0; 7.0, 5.0, 4.5];
%! brands = bsxfun (@times, ones (6, 1), [1, 2, 3]);
%! popper = bsxfun (@times, [1; 1; 1; 2; 2; 2], ones (1, 3));
%! X = [brands(:), popper(:)];
%! mdl = fitlm (X, popcorn(:), 'interactions', 'CategoricalVars', [1, 2]);
%! assert (mdl.Coefficients.Estimate(1),  5.66666666666667, 1e-09);
%! assert (mdl.Coefficients.Estimate(2), -1.33333333333333, 1e-09);
%! assert (mdl.Coefficients.Estimate(3), -2.16666666666667, 1e-09);
%! assert (mdl.Coefficients.Estimate(4),  1.16666666666667, 1e-09);
%! assert (mdl.Coefficients.Estimate(6), -0.333333333333334, 1e-09);
%! assert (mdl.Coefficients.Estimate(7), -0.166666666666667, 1e-09);
%! assert (mdl.Coefficients.SE(1), 0.215165741455965, 1e-09);
%! assert (mdl.Coefficients.SE(2), 0.304290309725089, 1e-09);
%! assert (mdl.Coefficients.SE(3), 0.304290309725089, 1e-09);
%! assert (mdl.Coefficients.SE(4), 0.304290309725089, 1e-09);
%! assert (mdl.Coefficients.SE(6), 0.43033148291193, 1e-09);
%! assert (mdl.Coefficients.SE(7), 0.43033148291193, 1e-09);
%! assert (mdl.Coefficients.tStat(1),  26.3362867542108,   1e-09);
%! assert (mdl.Coefficients.tStat(2),  -4.38178046004138,  1e-09);
%! assert (mdl.Coefficients.tStat(3),  -7.12039324756724,  1e-09);
%! assert (mdl.Coefficients.tStat(4),   3.83405790253621,  1e-09);
%! assert (mdl.Coefficients.tStat(6),  -0.774596669241495, 1e-09);
%! assert (mdl.Coefficients.tStat(7),  -0.387298334620748, 1e-09);
%! assert (mdl.Coefficients.pValue(1), 5.49841502258254e-12, 1e-09);
%! assert (mdl.Coefficients.pValue(2), 0.000893505495903642, 1e-09);
%! assert (mdl.Coefficients.pValue(3), 1.21291454302428e-05, 1e-09);
%! assert (mdl.Coefficients.pValue(4), 0.00237798044119407,  1e-09);
%! assert (mdl.Coefficients.pValue(6), 0.453570536021938,    1e-09);
%! assert (mdl.Coefficients.pValue(7), 0.705316781644046,    1e-09);
%! brands = {'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'; ...
%!           'Gourmet', 'National', 'Generic'};
%! popper = {'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; 'oil', 'oil', 'oil'; ...
%!           'air', 'air', 'air'; 'air', 'air', 'air'; 'air', 'air', 'air'};
%! T = table (brands(:), popper(:), 'VariableNames', {'brands', 'popper'});
%! mdl = fitlm (T, popcorn(:), 'interactions');

%!test
%! load carsmall
%! X = [Weight, Horsepower, Acceleration];
%! fitlm (X, MPG, 'constant');
%! mdl = fitlm (X, MPG, 'linear');
%! assert (mdl.Coefficients.Estimate(1),  47.9767628118615,     1e-09);
%! assert (mdl.Coefficients.Estimate(2),  -0.00654155878851796, 1e-09);
%! assert (mdl.Coefficients.Estimate(3),  -0.0429433065881864,  1e-09);
%! assert (mdl.Coefficients.Estimate(4),  -0.0115826516894871,  1e-09);
%! assert (mdl.Coefficients.SE(1), 3.87851641748551,            1e-09);
%! assert (mdl.Coefficients.SE(2), 0.00112741016370336,         1e-09);
%! assert (mdl.Coefficients.SE(3), 0.0243130608813806,          1e-09);
%! assert (mdl.Coefficients.SE(4), 0.193325043113178,           1e-09);
%! assert (mdl.Coefficients.tStat(1),  12.369874881944,         1e-09);
%! assert (mdl.Coefficients.tStat(2),  -5.80228828790225,       1e-09);
%! assert (mdl.Coefficients.tStat(3),  -1.76626492228599,       1e-09);
%! assert (mdl.Coefficients.tStat(4),  -0.0599128364487485,     1e-09);
%! assert (mdl.Coefficients.pValue(1), 4.89570341688996e-21,    1e-09);
%! assert (mdl.Coefficients.pValue(2), 9.87424814144e-08,       1e-09);
%! assert (mdl.Coefficients.pValue(3), 0.0807803098213114,      1e-09);
%! assert (mdl.Coefficients.pValue(4), 0.952359384151778,       1e-09);

%!shared X, y, yl, T1, T2, T3, C
%! X  = [1 2; 3 4; 5 6];
%! y  = [2; 4; 5];
%! yl = logical ([1; 0; 1]);
%! T1 = table ([1;2;3], [4;5;6], 'VariableNames', {'x1','x2'});
%! T2 = table ([1;2;3], [4;5;6], 'VariableNames', {'x1','y'});
%! T3 = table ([1;2;3], [4;5;6], [2;4;5], 'VariableNames', {'x1','x2','y'});
%! C  = categorical ({'a';'b';'a'});

%!test
%! assert (class (fitlm (X, y)), 'LinearModel');
%!test
%! assert (class (fitlm (X, yl)), 'LinearModel');
%!test
%! assert (class (fitlm (X, y, 'linear')), 'LinearModel');
%!test
%! assert (class (fitlm (X, y, [1 0; 0 1])), 'LinearModel');
%!test
%! assert (class (fitlm (X, y, 'Intercept', false)), 'LinearModel');
%!test
%! assert (class (fitlm (X, y, 'linear', 'Weights', [1;2;1])), 'LinearModel');

%!test
%! mdl = fitlm (C, y);
%! assert (class (mdl), 'LinearModel');
%! assert (mdl.VariableNames, {'x1', 'y'});
%!test
%! mdl = fitlm (C, y, 'VarNames', {'grp', 'score'});
%! assert (mdl.VariableNames, {'grp', 'score'});
%!test
%! assert (class (fitlm (C, y, 'Intercept', false)), 'LinearModel');

%!test
%! assert (class (fitlm (T2)), 'LinearModel');
%!test
%! assert (class (fitlm (T3)), 'LinearModel');
%!test
%! assert (class (fitlm (T3, 'Exclude', [2])), 'LinearModel');
%!test
%! assert (class (fitlm (T2, 'y')), 'LinearModel');
%!test
%! assert (class (fitlm (T3, 'x1')), 'LinearModel');
%!test
%! assert (class (fitlm (T3, 'y ~ x1 + x2')), 'LinearModel');
%!test
%! assert (class (fitlm (T1, 'linear')), 'LinearModel');
%!test
%! assert (class (fitlm (T1, [2;4;5])), 'LinearModel');
%!test
%! assert (class (fitlm (T2, [0 0; 1 0])), 'LinearModel');
%!test
%! assert (class (fitlm (T2, 'y', 'linear', 'Intercept', false)), 'LinearModel');

%!error <Not enough input arguments> fitlm ()
%!error <Predictor variables must be numeric vectors, numeric matrices, or categorical vectors> ...
%! fitlm ('hello', y)
%!error <Predictor variables must be numeric vectors, numeric matrices, or categorical vectors> ...
%! fitlm (struct ('a', 1), [1;2])
%!error <Predictor variables must be numeric vectors, numeric matrices, or categorical vectors> ...
%! fitlm (categorical ([1 2; 3 4]), y)
%!error <Y argument is required unless X is a dataset or table> ...
%! fitlm (C)
%!error <Y argument is required unless X is a dataset or table> ... 
%! fitlm (C, 'Intercept', false)
%!error <Predictor and response variables must have the same length> ...
%! fitlm (C, [1;2])
%!error <Response variable must be a numeric vector> fitlm (C, {'a';'b';'a'})
%!error <Y argument is required unless X is a dataset or table> fitlm (X)
%!error <Y argument is required unless X is a dataset or table> ... 
%! fitlm (X, 'Weights', [1;1;1])
%!error <Predictor and response variables must have the same length> ...
%! fitlm (X, [])
%!error <Predictor and response variables must have the same length> ...
%! fitlm (X, [1;2])
%!error <Response variable must be a numeric vector> fitlm (X, ones (3, 2))
%!error <Response variable must be a numeric vector> fitlm (X, {'1';'2';'3'})
%!error <The terms matrix must have one column for each variable in the dataset or table> ...
%! fitlm (T1, [])
%!error <Cannot determine the response variable from the terms matrix> ...
%! fitlm (T1, ones (1, 2))
%!error <Predictor and response variables must have the same length> ...
%! fitlm (T1, ones (4, 1))
%!error <Predictor and response variables must have the same length> ...
%! fitlm (T1, ones (2, 3))
%!error <invalid second argument for table input> fitlm (T1, {1, 2})
%!error <Name-Value arguments must be in pairs> fitlm (T2, 'y ~ x1', 'linear')
%!error <Name-Value arguments must be in pairs> fitlm (T1, 'linear', 'Weights')
%!error <Name-Value arguments must be in pairs> fitlm (T2, [0 0; 1 0], 'Weights')