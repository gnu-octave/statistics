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
## @deftypefn  {statistics} {@var{beta} =} nlinfit (@var{X}, @var{y}, @var{modelfun}, @var{beta0})
## @deftypefnx {statistics} {@var{beta} =} nlinfit (@dots{}, @var{options})
## @deftypefnx {statistics} {@var{beta} =} nlinfit (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{beta}, @var{R}, @var{J}, @var{CovB}, @var{MSE}, @var{ErrorModelInfo}] =} nlinfit (@dots{})
##
## Fit a nonlinear regression model.
##
## @code{@var{beta} = nlinfit (@var{X}, @var{y}, @var{modelfun}, @var{beta0})}
## estimates the coefficients of the nonlinear regression model
## @code{@var{y} = @var{modelfun} (@var{beta}, @var{X})} by iteratively
## minimizing the (possibly weighted) sum of squared residuals, starting from
## the initial coefficient vector @var{beta0}.  The fit uses the
## Levenberg-Marquardt algorithm with a numerically computed Jacobian.
##
## @itemize
## @item @var{X} is a matrix of predictor values.  @code{nlinfit} does not
## interpret the columns of @var{X}; the array is passed unchanged as the second
## argument of @var{modelfun}, so its shape is whatever @var{modelfun} expects.
## @item @var{y} is a numeric vector of responses, one element per observation.
## @item @var{modelfun} is a function handle @code{@@(@var{b}, @var{X})}
## returning a vector of fitted responses the same size as @var{y}.
## @item @var{beta0} is a numeric vector of initial values for the coefficients.
## @end itemize
##
## Additional options are given either as a statset-style @var{options}
## structure or as @qcode{Name}/@qcode{Value} pairs (or both).  The supported
## options are:
##
## @multitable @columnfractions 0.2 0.78
## @headitem @var{Name} @tab @var{Value}
## @item @qcode{'Weights'} @tab A vector of nonnegative observation weights, or
## a function handle @code{@@(@var{yhat})} returning such a vector.  Weighted
## least squares is used.
## @item @qcode{'ErrorModel'} @tab The form of the error variance:
## @qcode{'constant'} (default), @qcode{'proportional'}, or @qcode{'combined'}.
## @item @qcode{'ErrorParameters'} @tab Initial values for the error-model
## parameters.
## @item @qcode{'RobustWgtFun'} @tab The name of a robust weight function
## (@qcode{'andrews'}, @qcode{'bisquare'}, @qcode{'cauchy'}, @qcode{'fair'},
## @qcode{'huber'}, @qcode{'logistic'}, @qcode{'talwar'}, or @qcode{'welsch'}),
## enabling robust iteratively reweighted least squares.
## @item @qcode{'Tune'} @tab The tuning constant for the robust weight function.
## @item @qcode{'Options'} @tab A statset-style structure whose @qcode{MaxIter},
## @qcode{TolFun}, @qcode{TolX}, and @qcode{DerivStep} fields override the
## corresponding defaults.
## @end multitable
##
## The remaining outputs describe the converged fit: @var{R} is the vector of
## raw residuals @code{@var{y} - @var{modelfun} (@var{beta}, @var{X})}, @var{J}
## is the Jacobian of @var{modelfun} with respect to @var{beta} at the solution,
## @var{CovB} is the estimated covariance matrix of the coefficients, @var{MSE}
## is the mean squared error, and @var{ErrorModelInfo} is a structure describing
## the fitted error model.
##
## @subheading Algorithm
##
## The coefficients are estimated by the Levenberg-Marquardt algorithm using a
## numerically computed (forward-difference) Jacobian.  For an ordinary or
## weighted fit the coefficient covariance is @code{@var{CovB} = @var{MSE} *
## inv (@var{J}' * @var{W} * @var{J})}, where @var{W} is the diagonal matrix of
## observation weights and @var{MSE} is the weighted residual sum of squares
## divided by the error degrees of freedom @math{n - p} (with @math{p}
## coefficients).  A non-constant @qcode{'ErrorModel'} is fitted by generalized
## least squares, re-deriving the observation weights from the fitted values
## each iteration; the @qcode{'proportional'} model weights each observation by
## the inverse squared fitted value, and @var{MSE} then estimates the
## proportionality constant of the variance.
##
## For a robust fit (@qcode{'RobustWgtFun'}) the coefficients are found by
## iteratively reweighted least squares applied to leverage-adjusted residuals
## (the leverage is taken from the ordinary fit and held fixed).  The robust
## coefficient covariance follows the Street-Carroll-Ruppert convention, the
## same one used by @code{robustfit}: @code{@var{CovB} = s^2 * inv (@var{J}' *
## @var{J})} and @math{@var{MSE} = s^2}, where the scale @code{s} blends the
## ordinary-fit scale @code{ols_s} with the robust scale @code{robust_s} at the
## solution as @math{s^2 = (p^2 * ols_s^2 + n * robust_s^2) / (n + p^2)}, taken
## to be at least @code{robust_s}.
##
## @seealso{fitnlm, nlparci, nlpredci, NonLinearModel, robustfit}
## @end deftypefn

function [beta, R, J, CovB, MSE, ErrorModelInfo] = nlinfit (X, y, modelfun, ...
                                                            beta0, varargin)

  if (nargin < 4)
    print_usage ();
  endif

  ## Validate the core inputs.
  if (! is_function_handle (modelfun))
    error ("nlinfit: MODELFUN must be a function handle.");
  endif
  if (! (isnumeric (beta0) && isvector (beta0) && isreal (beta0)))
    error ("nlinfit: BETA0 must be a real numeric vector.");
  endif
  if (! (isnumeric (y) && isreal (y)))
    error ("nlinfit: Y must be a real numeric vector.");
  endif

  beta0 = beta0(:);
  y     = y(:);
  n     = numel (y);
  p     = numel (beta0);

  ## Parse the options structure and the Name/Value pairs.
  opts = parse_options (varargin);

  ## Resolve observation weights (constant unless a vector/handle is given).
  if (isempty (opts.Weights))
    w = ones (n, 1);
  elseif (is_function_handle (opts.Weights))
    w = [];                    # deferred: depends on the fitted values
  else
    w = check_weights (opts.Weights, n);
  endif

  ## Fit.  A robust weight function triggers iteratively reweighted least
  ## squares; a non-constant error model triggers generalized least squares
  ## with variance-dependent weights; both wrap the weighted LM solver.
  dfe = n - p;
  if (! isempty (opts.RobustWgtFun))
    ## Robust fit: MSE and CovB use the Street-Carroll-Ruppert robust scale
    ## (as in robustfit), CovB = MSE * inv (J' * J) at the robust solution.
    [wgtfun, tune] = robust_weight_function (opts.RobustWgtFun, opts.Tune);
    [beta, R, J, ols_s, adj] = robust_fit (X, y, modelfun, beta0, w, ...
                                           wgtfun, tune, opts);
    [MSE, CovB] = robust_covariance (R, J, adj, ols_s, wgtfun, tune, n, p);
  else
    if (! strcmp (opts.ErrorModel, "constant"))
      [beta, R, J, extraw] = errormodel_fit (X, y, modelfun, beta0, w, opts);
    else
      [beta, R, J] = lm_fit (X, y, modelfun, beta0, w, opts);
      extraw = ones (n, 1);
    endif

    ## Effective per-observation weights combine the user weights (which may be
    ## a function of the fitted values) with any error-model weights.
    if (is_function_handle (opts.Weights))
      w = check_weights (opts.Weights (modelfun (beta, X)), n);
    endif
    weff = w .* extraw;

    ## Mean squared error and CovB = MSE * inv (J' * Weff * J).
    if (dfe > 0)
      MSE = sum (weff .* R .^ 2) / dfe;
    else
      MSE = NaN;
    endif
    Jw   = sqrt (weff) .* J;
    CovB = MSE * pinv (Jw' * Jw);
    CovB = (CovB + CovB') / 2;              # symmetrise against round-off
  endif

  ErrorModelInfo = error_model_info (opts, MSE, p);

endfunction

## ---------------------------------------------------------------------------
## Generalized least squares for a non-constant error model.  The error weights
## depend on the fitted values and are refreshed each outer iteration.
##   proportional:  var_i = phi * mu_i^2          -> weight 1 / mu_i^2
##   combined:      var_i = phi * (a + |mu_i|)^2   -> weight 1 / (a + |mu_i|)^2
function [beta, R, J, ew] = errormodel_fit (X, y, modelfun, beta0, wuser, opts)

  if (isempty (wuser))
    wuser = ones (numel (y), 1);
  endif
  beta = beta0;
  ew   = ones (numel (y), 1);

  for outer = 1:opts.MaxIter
    ewprev = ew;
    [beta, R, J] = lm_fit (X, y, modelfun, beta, wuser .* ew, opts);
    mu = y - R;
    ew = errormodel_weights (opts.ErrorModel, mu);
    if (max (abs (ew - ewprev) ./ (ewprev + eps)) < opts.TolX)
      break;
    endif
  endfor

endfunction

## ---------------------------------------------------------------------------
## Per-observation error-model weights from the current fitted values MU.
function ew = errormodel_weights (errmodel, mu)

  am = abs (mu);
  am(am < eps) = eps;                      # guard against division by zero
  switch (errmodel)
    case "proportional"
      ew = 1 ./ am .^ 2;
    case "combined"
      ## a is a small offset relative to the scale of the fitted values.
      a  = 0.5 * mean (am);
      ew = 1 ./ (a + am) .^ 2;
    otherwise
      ew = ones (size (mu));
  endswitch

endfunction

## ---------------------------------------------------------------------------
## Plain (optionally weighted) Levenberg-Marquardt nonlinear least squares.
## Returns the coefficients, raw residuals, and Jacobian at the solution.
function [beta, R, J] = lm_fit (X, y, modelfun, beta0, w, opts)

  if (isempty (w))
    w = ones (numel (y), 1);
  endif
  sw = sqrt (w);

  beta   = beta0;
  yhat   = modelfun (beta, X);
  R      = y - yhat;
  sse    = sum (w .* R .^ 2);
  lambda = 1e-2;                           # Marquardt damping
  J      = nlfun_jacobian (modelfun, beta, X, opts.DerivStep);

  for iter = 1:opts.MaxIter
    Jw   = sw .* J;
    Rw   = sw .* R;
    JtJ  = Jw' * Jw;
    Jtr  = Jw' * Rw;
    diagJtJ = diag (diag (JtJ));

    ## Try damped steps, increasing lambda until the objective decreases.
    stepok = false;
    for inner = 1:20
      A     = JtJ + lambda * diagJtJ;
      delta = pinv (A) * Jtr;
      bnew  = beta + delta;
      rnew  = y - modelfun (bnew, X);
      ssenew = sum (w .* rnew .^ 2);
      if (isfinite (ssenew) && ssenew < sse)
        stepok = true;
        break;
      endif
      lambda = lambda * 10;
    endfor

    if (! stepok)
      break;                               # cannot improve: converged/stalled
    endif

    ## Accept the step and relax the damping.
    relchg = abs (ssenew - sse) / (sse + eps);
    stepsz = max (abs (delta) ./ (abs (beta) + eps));
    beta   = bnew;
    R      = rnew;
    sse    = ssenew;
    lambda = max (lambda / 10, 1e-12);
    J      = nlfun_jacobian (modelfun, beta, X, opts.DerivStep);

    if (relchg < opts.TolFun || stepsz < opts.TolX)
      break;
    endif
  endfor

endfunction

## ---------------------------------------------------------------------------
## Robust iteratively reweighted least squares around LM_FIT, following the
## robustfit scheme: residuals are leverage-adjusted before scaling, and the
## leverage (from the ordinary fit) is held fixed.  Returns the ordinary-fit
## scale OLS_S and the leverage adjustment ADJ for the covariance step.
function [beta, R, J, ols_s, adj] = robust_fit (X, y, modelfun, beta0, w, ...
                                                wgtfun, tune, opts)

  if (isempty (w))
    w = ones (numel (y), 1);
  endif
  n = numel (y);
  p = numel (beta0);

  ## Ordinary weighted fit: OLS scale and the (fixed) leverage adjustment.
  [beta, R, J] = lm_fit (X, y, modelfun, beta0, w, opts);
  ols_s = norm (R) / sqrt (max (n - p, 1));
  Jw    = sqrt (w) .* J;
  h     = min (0.9999, diag (Jw * pinv (Jw' * Jw) * Jw'));
  adj   = 1 ./ sqrt (1 - h);

  for iter = 1:opts.MaxIter
    betaprev = beta;
    s = madsigma (R .* adj, p);
    if (s == 0)
      break;
    endif
    rw = max (wgtfun (R .* adj ./ (s * tune)), 0);
    [beta, R, J] = lm_fit (X, y, modelfun, beta, w .* rw, opts);
    if (all (abs (beta - betaprev) <= sqrt (eps) * max (abs (beta), ...
                                                        abs (betaprev))))
      break;
    endif
  endfor

endfunction

## ---------------------------------------------------------------------------
## Street-Carroll-Ruppert robust coefficient covariance (matching robustfit).
## The scale blends the ordinary-fit scale with the robust scale, and
## CovB = s^2 * inv (J' * J) at the robust solution.
function [MSE, CovB] = robust_covariance (R, J, adj, ols_s, wgtfun, tune, n, p)

  radj  = R .* adj;
  mad_s = madsigma (radj, p);
  if (mad_s == 0)
    mad_s = 1;
  endif
  z = radj ./ (mad_s * tune);
  w = wgtfun (z);
  if (all (w == 1))
    robust_s = ols_s;
  else
    ## psi = z .* w and its derivative (central difference, as robustfit does
    ## for user weight functions).
    psi  = z .* w;
    d    = 1e-6;
    psip = ((z + d) .* wgtfun (z + d) - (z - d) .* wgtfun (z - d)) / (2 * d);
    K    = 1 + (p / n) * var (psip) / mean (psip) ^ 2;
    robust_s = tune * mad_s * sqrt (mean (psi .^ 2)) / mean (psip) * K;
  endif
  s = sqrt ((p ^ 2 * ols_s ^ 2 + n * robust_s ^ 2) / (n + p ^ 2));
  s = max (s, robust_s);

  MSE  = s ^ 2;
  CovB = MSE * pinv (J' * J);
  CovB = (CovB + CovB') / 2;

endfunction

## ---------------------------------------------------------------------------
## Robust scale from the median absolute deviation of the residuals, dropping
## the P-1 smallest to account for the fitted parameters (matches robustfit).
function s = madsigma (r, p)

  rs = sort (abs (r));
  s = median (rs(max (1, p):end)) / 0.6745;

endfunction

## ---------------------------------------------------------------------------
## Return the robust weight function handle and its default tuning constant.
function [wgtfun, tune] = robust_weight_function (name, tune)

  switch (lower (name))
    case 'andrews'
      deftune = 1.339;
      wgtfun  = @(u) (abs (u) < pi) .* sin (u) ./ (u + (u == 0));
    case 'bisquare'
      deftune = 4.685;
      wgtfun  = @(u) (abs (u) < 1) .* (1 - u .^ 2) .^ 2;
    case 'cauchy'
      deftune = 2.385;
      wgtfun  = @(u) 1 ./ (1 + u .^ 2);
    case 'fair'
      deftune = 1.400;
      wgtfun  = @(u) 1 ./ (1 + abs (u));
    case 'huber'
      deftune = 1.345;
      wgtfun  = @(u) 1 ./ max (1, abs (u));
    case 'logistic'
      deftune = 1.205;
      wgtfun  = @(u) tanh (u) ./ (u + (u == 0));
    case 'talwar'
      deftune = 2.795;
      wgtfun  = @(u) 1 .* (abs (u) < 1);
    case 'welsch'
      deftune = 2.985;
      wgtfun  = @(u) exp (- (u .^ 2));
    otherwise
      error ("nlinfit: unknown RobustWgtFun '%s'.", name);
  endswitch
  if (isempty (tune))
    tune = deftune;
  endif

endfunction

## ---------------------------------------------------------------------------
## Assemble the ErrorModelInfo structure returned as the sixth output.  The
## Scheffe dimension for simultaneous prediction is p for a non-constant error
## model and p + 1 for the constant model (the extra observation error term).
function info = error_model_info (opts, MSE, p)

  info = struct ();
  info.ErrorModel       = opts.ErrorModel;
  info.ErrorParameters  = sqrt (MSE);
  if (strcmp (opts.ErrorModel, "proportional"))
    info.ErrorVariance  = @(x) MSE * abs (x) .^ 2;
  else
    info.ErrorVariance  = @(x) MSE * ones (size (x, 1), 1);
  endif
  info.MSE              = MSE;
  info.ScheffeSimPred   = p + strcmp (opts.ErrorModel, "constant");
  info.WeightFunction   = is_function_handle (opts.Weights);
  info.FixedWeights     = (! isempty (opts.Weights) ...
                           && ! is_function_handle (opts.Weights));
  info.RobustWeightFunction = ! isempty (opts.RobustWgtFun);

endfunction

## ---------------------------------------------------------------------------
## Validate a numeric weight vector.
function w = check_weights (w, n)

  if (! (isnumeric (w) && isvector (w) && isreal (w)))
    error ("nlinfit: WEIGHTS must be a real numeric vector.");
  endif
  w = w(:);
  if (numel (w) != n)
    error ("nlinfit: WEIGHTS must have one element per observation.");
  endif
  if (any (w < 0) || any (isnan (w)))
    error ("nlinfit: WEIGHTS must be nonnegative.");
  endif

endfunction

## ---------------------------------------------------------------------------
## Parse an optional leading statset structure followed by Name/Value pairs.
function opts = parse_options (args)

  ## Defaults.
  opts = struct ("Weights", [], "ErrorModel", "constant", ...
                 "ErrorParameters", [], "RobustWgtFun", [], "Tune", [], ...
                 "MaxIter", 200, "TolFun", 1e-8, "TolX", 1e-8, ...
                 "DerivStep", eps ^ (1/3));

  ## An options structure may lead the Name/Value pairs.
  if (! isempty (args) && isstruct (args{1}))
    opts = merge_statset (opts, args{1});
    args = args(2:end);
  endif

  if (mod (numel (args), 2) != 0)
    error ("nlinfit: Name/Value arguments must come in pairs.");
  endif

  for k = 1:2:numel (args)
    name = args{k};
    val  = args{k+1};
    if (! ischar (name))
      error ("nlinfit: parameter names must be character vectors.");
    endif
    switch (lower (name))
      case 'weights'
        opts.Weights = val;
      case 'errormodel'
        opts.ErrorModel = lower (val);
      case 'errorparameters'
        opts.ErrorParameters = val;
      case 'robustwgtfun'
        opts.RobustWgtFun = val;
      case 'tune'
        opts.Tune = val;
      case 'options'
        opts = merge_statset (opts, val);
      otherwise
        error ("nlinfit: unknown parameter name '%s'.", name);
    endswitch
  endfor

  if (! any (strcmp (opts.ErrorModel, ...
                     {'constant', 'proportional', 'combined'})))
    error ("nlinfit: unknown ErrorModel '%s'.", opts.ErrorModel);
  endif

endfunction

## ---------------------------------------------------------------------------
## Copy the recognised statset fields from S into OPTS.
function opts = merge_statset (opts, s)

  if (! isstruct (s))
    error ("nlinfit: OPTIONS must be a structure.");
  endif
  fn = fieldnames (s);
  for k = 1:numel (fn)
    val = s.(fn{k});
    if (isempty (val))
      continue;
    endif
    switch (lower (fn{k}))
      case 'maxiter'
        opts.MaxIter = val;
      case 'tolfun'
        opts.TolFun = val;
      case 'tolx'
        opts.TolX = val;
      case 'derivstep'
        opts.DerivStep = val;
      case {'robustwgtfun', 'wgtfun'}
        opts.RobustWgtFun = val;
      case 'tune'
        opts.Tune = val;
    endswitch
  endfor

endfunction

%!demo
%! ## Fit an exponential growth model y = b1 * exp (b2 * x).
%! x = [1:10]';
%! y = [2.1;2.9;4.2;5.3;7.1;9.4;12.8;16.5;22.1;29.8];
%! modelfun = @(b, x) b(1) .* exp (b(2) .* x);
%! beta = nlinfit (x, y, modelfun, [1; 0.3])

%!demo
%! ## Robust fitting downweights a gross outlier (5th point corrupted).
%! x = [1:10]';
%! y = [2.1;2.9;4.2;5.3;7.1;9.4;12.8;16.5;22.1;29.8];
%! y(5) = 30;
%! modelfun = @(b, x) b(1) .* exp (b(2) .* x);
%! beta_ols = nlinfit (x, y, modelfun, [1; 0.3]);
%! beta_rob = nlinfit (x, y, modelfun, [1; 0.3], 'RobustWgtFun', 'bisquare');
%! [beta_ols, beta_rob]

%!shared x, y, modelfun, beta0
%! x = [1;2;3;4;5;6;7;8;9;10];
%! y = [2.1;2.9;4.2;5.3;7.1;9.4;12.8;16.5;22.1;29.8];
%! modelfun = @(b, x) b(1) .* exp (b(2) .* x);
%! beta0 = [1; 0.3];

## Values verified against MATLAB's nlinfit.
%!test
%! [beta, R, J, CovB, MSE] = nlinfit (x, y, modelfun, beta0);
%! assert_equal (beta, [1.683746951; 0.286911092], 1e-6);
%! assert_equal (MSE, 0.029221494, 1e-7);
%! assert_equal (sqrt (diag (CovB)), [0.035194484; 0.002350838], 1e-7);
%!test
%! ## The returned Jacobian equals the analytic model Jacobian.
%! [beta, R, J] = nlinfit (x, y, modelfun, beta0);
%! Ja = [exp(beta(2).*x), beta(1).*x.*exp(beta(2).*x)];
%! assert_equal (J, Ja, -1e-4);
%!test
%! ## Proportional error model reweights by 1 / mu^2 and changes the fit.
%! [beta, R, J, CovB, MSE, EMI] = nlinfit (x, y, modelfun, beta0, ...
%!                                         "ErrorModel", "proportional");
%! assert_equal (beta, [1.650639591; 0.289917266], 1e-6);
%! assert_equal (MSE, 9.988359e-4, 1e-9);
%! assert_equal (EMI.ErrorModel, "proportional");
%! assert_equal (EMI.ScheffeSimPred, 2);
%!test
%! ## Constant error model: ErrorModelInfo fields.
%! [beta, R, J, CovB, MSE, EMI] = nlinfit (x, y, modelfun, beta0);
%! assert_equal (EMI.ErrorModel, "constant");
%! assert_equal (EMI.ScheffeSimPred, 3);
%! assert_equal (EMI.ErrorParameters, 0.170942956, 1e-7);
%!test
%! ## Observation weights.
%! w = (1:10)';
%! [beta, R, J, CovB, MSE] = nlinfit (x, y, modelfun, beta0, "Weights", w);
%! assert_equal (beta, [1.677485713; 0.287323059], 1e-6);
%! assert_equal (MSE, 0.178838470, 1e-6);
%!test
%! ## Robust fitting resists a gross outlier: the estimate stays close to the
%! ## clean-data fit, far closer than ordinary least squares.  Coefficients,
%! ## MSE, and covariance match MATLAB's nlinfit (Street-Carroll-Ruppert scale).
%! yo = y; yo(5) = 30;
%! [br, Rr, Jr, Cr, Mr] = nlinfit (x, yo, modelfun, beta0, ...
%!                                 "RobustWgtFun", "bisquare");
%! bo = nlinfit (x, yo, modelfun, beta0);
%! bc = nlinfit (x, y, modelfun, beta0);
%! assert (max (abs (br - bc)) < max (abs (bo - bc)));
%! assert_equal (br, [1.679154218; 0.287198922], 1e-5);
%! assert_equal (Mr, 15.892255, 1e-2);
%! assert_equal (Cr, [0.671544, -0.044220; -0.044220, 0.003012], 1e-3);
%!test
%! ## The huber robust fit also matches MATLAB's coefficients.
%! yo = y; yo(5) = 30;
%! bh = nlinfit (x, yo, modelfun, beta0, "RobustWgtFun", "huber");
%! assert_equal (bh, [1.716586997; 0.284865506], 1e-5);

## Test input validation
%!error<Invalid call> nlinfit (1, 2, @(b, x) x)
%!error<nlinfit: MODELFUN must be a function handle.> nlinfit (1, 2, 3, 4)
%!error<nlinfit: BETA0 must be a real numeric vector.> ...
%! nlinfit ([1;2], [1;2], @(b, x) b(1) * ones (2, 1), "a")
%!error<nlinfit: unknown parameter name 'foo'.> ...
%! nlinfit ([1;2], [1;2], @(b, x) b(1) * ones (2, 1), 1, "foo", 1)
%!error<nlinfit: unknown ErrorModel 'bad'.> ...
%! nlinfit ([1;2], [1;2], @(b, x) b(1) * ones (2, 1), 1, "ErrorModel", "bad")
%!error<nlinfit: unknown RobustWgtFun 'bad'.> ...
%! nlinfit ([1;2], [1;2], @(b, x) b(1) * ones (2, 1), 1, "RobustWgtFun", "bad")
