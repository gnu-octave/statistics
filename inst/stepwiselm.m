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

  ## Optional positional modelspec
  valid_ms = {"constant", "linear", "interactions", "full"};
  opt_keys = {"Upper", "Lower", "Criterion", "PEnter", "PRemove",
              "NSteps", "Verbose", "VarNames", "Intercept"};

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

  ## Semantic validation
  Upper     = lower (Upper);
  Lower     = lower (Lower);
  modelspec = lower (modelspec);
  Criterion = lower (Criterion);

  if (! any (strcmp (Upper, valid_ms)))
    error ("stepwiselm: Upper must be constant, linear, interactions, full");
  endif
  if (! any (strcmp (Lower, valid_ms)))
    error ("stepwiselm: Lower must be constant, linear, interactions, full");
  endif
  if (! any (strcmp (modelspec, valid_ms)))
    error ("stepwiselm: modelspec must be constant, linear, interactions, full");
  endif

  mrank = struct ("constant", 0, "linear", 1, "interactions", 2, "full", 3);
  if (mrank.(Lower) > mrank.(Upper))
    error ("stepwiselm: Lower model cannot be more complex than Upper");
  endif

  ## Clamp starting model between Lower and Upper
  if (mrank.(modelspec) < mrank.(Lower))
    modelspec = Lower;
  endif
  if (mrank.(modelspec) > mrank.(Upper))
    modelspec = Upper;
  endif

  valid_crit = {"sse", "aic", "bic", "rsquared", "adjrsquared"};
  if (! any (strcmp (Criterion, valid_crit)))
    error ("stepwiselm: Criterion must be sse, aic, bic, rsquared, adjrsquared");
  endif

  ## PEnter / PRemove validation is criterion-aware
  ## For SSE: both must be probabilities in (0,1).
  ## For AIC/BIC/R2/AdjR2: can be any finite scalar (change thresholds).
  ## For ALL criteria:PEnter < PRemove strictly enforced.
  if (! (isscalar (PEnter) && isnumeric (PEnter) && isfinite (PEnter)))
    error ("stepwiselm: PEnter must be a finite scalar");
  endif
  if (! (isscalar (PRemove) && isnumeric (PRemove) && isfinite (PRemove)))
    error ("stepwiselm: PRemove must be a finite scalar");
  endif
  if (strcmp (Criterion, "sse"))
    if (PEnter <= 0 || PEnter >= 1)
      error ("stepwiselm: PEnter must be a scalar in (0,1) for SSE criterion");
    endif
    if (PRemove <= 0 || PRemove >= 1)
      error ("stepwiselm: PRemove must be a scalar in (0,1) for SSE criterion");
    endif
  endif
  ## Universal rule: PRemove > PEnter (confirmed for AIC in sm3.m output)
  if (PRemove <= PEnter)
    error ("stepwiselm: PRemove must be strictly greater than PEnter");
  endif
  if (! (isscalar (NSteps) && isnumeric (NSteps) && NSteps > 0))
    error ("stepwiselm: NSteps must be a positive scalar or Inf");
  endif
  if (! (isscalar (Verbose) && isnumeric (Verbose) && any (Verbose == [0 1 2])))
    error ("stepwiselm: Verbose must be 0, 1, or 2");
  endif
  if (! (isscalar (Intercept) && islogical (Intercept)))
    error ("stepwiselm: Intercept must be a logical scalar");
  endif

  ## Default variable names
  if (isempty (VarNames))
    VarNames = cell (1, p + 1);
    for k = 1:p
      VarNames{k} = sprintf ("x%d", k);
    endfor
    VarNames{p + 1} = "y";
  else
    if (! iscell (VarNames) || numel (VarNames) != p + 1)
      error ("stepwiselm: VarNames must be a cell array of length p+1");
    endif
  endif

  PredNames    = VarNames(1:p);
  ResponseName = VarNames{p + 1};

  ## Drop NaN rows
  wasnan = any (isnan ([X, y]), 2);
  Xc = X(! wasnan, :);
  yc = y(! wasnan);
  n  = rows (Xc);

  if (n < 2)
    error ("stepwiselm: fewer than 2 complete observations");
  endif

  ## Build augmented design matrix for Upper bound
  [Xaug, termNames] = swlm_build_terms (Xc, PredNames, Upper);
  p_aug = columns (Xaug);

  ## Lower bound: forced-in terms
  [~, lowerNames] = swlm_build_terms (Xc, PredNames, Lower);
  keep_mask = ismember (termNames, lowerNames);

  ## Starting model
  [~, startNames] = swlm_build_terms (Xc, PredNames, modelspec);
  inmodel_mask = ismember (termNames, startNames) | keep_mask;

  ## Stepwise selection
  if (strcmp (Criterion, "sse"))
    [finalmodel, history] = swlm_pval_step ( ...
        Xaug, yc, logical (inmodel_mask), logical (keep_mask), ...
        PEnter, PRemove, NSteps, Verbose, termNames, Intercept);
  else
    [finalmodel, history] = swlm_crit_step ( ...
        Xaug, yc, logical (inmodel_mask), logical (keep_mask), ...
        Criterion, NSteps, Verbose, termNames, Intercept);
  endif
