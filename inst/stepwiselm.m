## Copyright (C) 2026 Jayant Chauhan <0001jayant@gmail.com>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and / or modify it under
##the terms of the GNU General Public License as published by the Free Software
##Foundation; either version 3 of the License, or (at your option) any later
##version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see < http: // www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{mdl} =} stepwiselm (@var{X}, @var{y})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@var{X}, @var{y}, @var{modelspec})
## @deftypefnx {statistics} {@var{mdl} =} stepwiselm (@dots{}, @var{Name}, @var{Value})
##
## Stepwise linear regression with automatic predictor selection.
##
## @code{stepwiselm} fits a linear regression model by iteratively adding
## and removing predictor terms.  Starting from a lower-bound model, the
## function adds the term that most improves the fit (while its addition
## p-value stays below @code{PEnter}) and removes the term whose p-value
## exceeds @code{PRemove}, repeating until the model stabilises or
## @code{NSteps} iterations are exhausted.
##
## @var{X} is an @var{n}-by-@var{p} numeric matrix of predictor variables.
## @var{y} is an @var{n}-by-1 numeric response vector.
##
## @var{modelspec} specifies the starting model and must be one of:
## @itemize
## @item @qcode{"constant"} -- intercept only (no predictors)
## @item @qcode{"linear"} (default) -- intercept plus all linear terms
## @item @qcode{"interactions"} -- linear terms plus all pairwise products
## @item @qcode{"full"} -- all main effects and higher-order interactions
## @end itemize
##
## @subheading Name-Value Options
##
## @table @asis
## @item @qcode{"Upper"}
## Upper-bound model: maximum complexity allowed.  Accepts the same
## keywords as @var{modelspec}.  Default is @qcode{"linear"}.
##
## @item @qcode{"Lower"}
## Lower-bound model: terms that are always kept in the model.
## Default is @qcode{"constant"} (intercept only).
##
## @item @qcode{"Criterion"}
## Selection criterion.  One of @qcode{"sse"} (p-value, default),
## @qcode{"aic"}, @qcode{"bic"}, @qcode{"rsquared"},
## @qcode{"adjrsquared"}.
##
## @item @qcode{"PEnter"}
## Threshold for adding a term.  For @qcode{"sse"}: maximum p-value for
## entry; default @code{0.05}, must be in @code{(0,1)}.
## For @qcode{"aic"}/@qcode{"bic"}: maximum AIC/BIC change allowed for
## entry; default @code{0} (add term only if criterion strictly decreases).
## For @qcode{"rsquared"}: minimum R@sup{2} increase required; default
## @code{0.1}.  For @qcode{"adjrsquared"}: default @code{0}.
## Must satisfy @code{PEnter < PRemove} for all criteria.
## @strong{Note:} for non-@qcode{"sse"} criteria, PEnter and PRemove are
## accepted but currently ignored internally; the greedy (best-move) rule
## is used instead.  This matches MATLAB's default-threshold behaviour.
##
## @item @qcode{"PRemove"}
## Threshold for removing a term.  For @qcode{"sse"}: minimum p-value to
## trigger removal; default @code{0.10}, must be in @code{(0,1)}, and
## must satisfy @code{PRemove > PEnter} strictly.
## For @qcode{"aic"}/@qcode{"bic"}: default @code{0.01}.
## For @qcode{"rsquared"}: default @code{0.05}.
## For @qcode{"adjrsquared"}: default @code{-0.05}.
## Must satisfy @code{PEnter < PRemove} for all criteria.
##
## @item @qcode{"NSteps"}
## Maximum number of stepwise iterations.  Default @code{Inf}.
##
## @item @qcode{"Verbose"}
## Display level: @code{0} (silent), @code{1} (final summary, default),
## @code{2} (per-step detail).
##
## @item @qcode{"VarNames"}
## Cell array of length @var{p}+1 giving names for each predictor
## (columns 1 to @var{p}) and the response (column @var{p}+1).
## Defaults to @{@qcode{"x1"},@dots{},@qcode{"xp"},@qcode{"y"}@}.
##
## @item @qcode{"Intercept"}
## Logical scalar;
if @code{true} (default) include an intercept term.
## @end table
##
## @subheading Return Value
##
## @var{mdl} is a structure with fields mirroring MATLAB's
## @code{LinearModel} object:
##
## @itemize
## @item @var{mdl}.Coefficients -- struct with Estimate, SE, tStat,
## pValue for each coefficient in the final model.
## @item @var{mdl}.CoefficientNames -- cell array of coefficient names.
## @item @var{mdl}.Fitted -- fitted values (@var{n}-by-1).
## @item @var{mdl}.Residuals -- struct with Raw, Pearson, Standardized.
## @item @var{mdl}.Rsquared -- struct with Ordinary and Adjusted.
## @item @var{mdl}.ModelCriterion -- struct with AIC, BIC, AICc.
## @item @var{mdl}.RMSE, @var{mdl}.MSE -- error magnitude measures.
## @item @var{mdl}.SSE, @var{mdl}.SST, @var{mdl}.SSR -- sums of squares.
## @item @var{mdl}.DFE -- residual degrees of freedom.
## @item @var{mdl}.NumObservations -- number of (non-NaN) observations.
## @item @var{mdl}.NumCoefficients -- total estimated coefficients.
## @item @var{mdl}.NumPredictors -- predictors in the final model.
## @item @var{mdl}.PredictorNames -- all input predictor names.
## @item @var{mdl}.SelectedPredictors -- names of selected predictors.
## @item @var{mdl}.ResponseName -- name of the response variable.
## @item @var{mdl}.Formula -- formula string of the final model.
## @item @var{mdl}.History -- struct tracking the stepwise steps.
## @item @var{mdl}.NumRemovedNaN -- observations dropped due to NaN.
## @end itemize
##
## @seealso{
  stepwisefit, fitlm, regress
}
##@end deftypefn

function mdl = stepwiselm (X, y, varargin)

  ## Input validation
  if (nargin < 2)
    error ("stepwiselm: at least two input arguments required");
  endif

  if (! isnumeric (X) || ! ismatrix (X))
    error ("stepwiselm: X must be a numeric matrix");
  endif

  if (! isnumeric (y) || ! isvector (y))
    error ("stepwiselm: y must be a numeric vector");
  endif

  y = y(:);

  if (rows (X) != rows (y))
    error ("stepwiselm: X and y must have the same number of rows");
  endif

  [~, p] = size (X);

  ## ── Optional positional
      modelspec ──────────────────────────────── valid_ms = {
          "constant", "linear", "interactions", "full"};
  opt_keys = {"Upper",      "Lower",   "Criterion", "PEnter",   "PRemove",
              ... "NSteps", "Verbose", "VarNames",  "Intercept"};

  optStart = 1;
  modelspec = "constant";

  if (!isempty(varargin))
    v1 = varargin{1};
  if (ischar(v1) && any(strcmpi(v1, valid_ms)))
    ##Valid modelspec keyword provided as 3rd positional argument modelspec =
        lower(v1);
  optStart = 2;
  elseif(ischar(v1) && !any(strcmpi(v1, opt_keys)))##Unknown string — not a
      valid modelspec and not a known option key.##Let it fall into the Name
      - Value parser so pairedArgs can flag it##as an
        unrecognized option with the standard error message.optStart = 1;
  elseif(!ischar(v1) && !isnumeric(v1))
      error("stepwiselm: third argument must be a modelspec string");
  endif endif

  ## Parse Name-Value pairs
  rest = varargin(optStart:end);

  dfVals = {"interactions", "constant", "sse", 0.05, 0.10, Inf, 1, {}, true};
  [Upper, Lower, Criterion, PEnter, PRemove, NSteps, Verbose, ...
   VarNames, Intercept, extras] = pairedArgs (opt_keys, dfVals, rest(:));

  if (! isempty (extras))
    if (ischar (extras{1}))
      error ("stepwiselm: unrecognized option '%s'", extras{1});
    else
      error ("stepwiselm: unrecognized input arguments");
    endif
  endif