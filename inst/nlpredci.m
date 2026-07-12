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
## @deftypefn  {statistics} {[@var{ypred}, @var{delta}] =} nlpredci (@var{modelfun}, @var{X}, @var{beta}, @var{resid}, @qcode{'Jacobian'}, @var{J})
## @deftypefnx {statistics} {[@var{ypred}, @var{delta}] =} nlpredci (@var{modelfun}, @var{X}, @var{beta}, @var{resid}, @qcode{'Covar'}, @var{CovB})
## @deftypefnx {statistics} {[@var{ypred}, @var{delta}] =} nlpredci (@dots{}, @var{Name}, @var{Value})
##
## Confidence intervals for predictions of a nonlinear regression.
##
## @code{[@var{ypred}, @var{delta}] = nlpredci (@var{modelfun}, @var{X},
## @var{beta}, @var{resid}, @qcode{'Jacobian'}, @var{J})} returns the predicted
## responses @var{ypred} of the model @code{@var{modelfun} (@var{beta},
## @var{X})} at the new predictor values @var{X}, together with the half-widths
## @var{delta} of the @math{100 (1 - @var{alpha})%} confidence intervals, so
## that @code{@var{ypred} - @var{delta}} and @code{@var{ypred} + @var{delta}}
## bound the response.  @var{beta}, @var{resid} (the residuals) and @var{J} (the
## Jacobian) come from @code{nlinfit}.
##
## Instead of the Jacobian, an estimated coefficient covariance may be supplied
## with @code{@qcode{'Covar'}, @var{CovB}}.  The following @qcode{Name}/
## @qcode{Value} pairs are also accepted:
##
## @multitable @columnfractions 0.2 0.78
## @headitem @var{Name} @tab @var{Value}
## @item @qcode{'MSE'} @tab The mean squared error from @code{nlinfit}, required
## with @qcode{'Covar'} for observation (prediction) intervals.
## @item @qcode{'PredOpt'} @tab @qcode{'curve'} (default) for confidence
## intervals on the fitted curve, or @qcode{'observation'} for prediction
## intervals on a new observation.
## @item @qcode{'SimOpt'} @tab @qcode{'off'} (default) for pointwise intervals,
## or @qcode{'on'} for simultaneous (Scheffe) intervals.
## @item @qcode{'Alpha'} @tab The significance level; the interval has
## confidence @math{100 (1 - @var{alpha})%} (default @var{alpha} = 0.05).
## @end multitable
##
## @subheading Algorithm
##
## Each half-width is @code{@var{delta} = c * sqrt (v)}.  The variance @code{v}
## of the fitted curve is @code{diag (@var{Jnew} * @var{V} * @var{Jnew}')},
## where @var{V} is the coefficient covariance (either @var{CovB}, or
## @code{@var{MSE} * inv (@var{J}' * @var{J})} when a Jacobian is supplied) and
## @var{Jnew} is the Jacobian of @var{modelfun} at @var{X}; an
## @qcode{'observation'} interval adds the error variance @var{MSE} to @code{v}.
## The critical value @code{c} is the Student's @math{t} quantile at
## @math{1 - @var{alpha}/2} with the error degrees of freedom for a pointwise
## interval, or the Scheffe value @code{sqrt (k * finv (1 - @var{alpha}, k,
## dfe))} for a simultaneous interval, where @math{k} is the number of
## coefficients (plus one for an observation interval).
##
## @seealso{nlinfit, nlparci, fitnlm, NonLinearModel}
## @end deftypefn

function [ypred, delta] = nlpredci (modelfun, X, beta, resid, varargin)

  if (nargin < 5)
    print_usage ();
  endif

  if (! is_function_handle (modelfun))
    error ("nlpredci: MODELFUN must be a function handle.");
  endif
  if (! (isnumeric (beta) && isvector (beta) && isreal (beta)))
    error ("nlpredci: BETA must be a real numeric vector.");
  endif
  if (! (isnumeric (resid) && isreal (resid)))
    error ("nlpredci: RESID must be a real numeric vector.");
  endif

  beta  = beta(:);
  resid = resid(:);
  p     = numel (beta);
  dfe   = numel (resid) - p;
  if (dfe <= 0)
    error ("nlpredci: not enough residuals to estimate the error variance.");
  endif

  ## Parse the covariance source and options.  A bare Jacobian may follow the
  ## residuals (legacy), otherwise 'Jacobian'/'Covar' keywords select it.
  J        = [];
  CovB     = [];
  MSE      = [];
  predopt  = 'curve';
  simopt   = 'off';
  alpha    = 0.05;
  args     = varargin;
  if (! ischar (args{1}))
    J    = args{1};
    args = args(2:end);
  endif
  if (mod (numel (args), 2) != 0)
    error ("nlpredci: Name/Value arguments must come in pairs.");
  endif
  for k = 1:2:numel (args)
    switch (lower (args{k}))
      case 'jacobian'
        J = args{k+1};
      case 'covar'
        CovB = args{k+1};
      case 'mse'
        MSE = args{k+1};
      case 'predopt'
        predopt = lower (args{k+1});
      case 'simopt'
        simopt = lower (args{k+1});
      case 'alpha'
        alpha = args{k+1};
      case {'weights', 'errormodelinfo'}
        ## accepted for compatibility; not needed for the covariance path
      otherwise
        error ("nlpredci: unknown parameter name '%s'.", args{k});
    endswitch
  endfor

  if (! (isscalar (alpha) && isreal (alpha) && alpha > 0 && alpha < 1))
    error ("nlpredci: ALPHA must be a scalar in the range (0, 1).");
  endif
  if (! any (strcmp (predopt, {'curve', 'observation'})))
    error ("nlpredci: PREDOPT must be 'curve' or 'observation'.");
  endif
  if (! any (strcmp (simopt, {'on', 'off'})))
    error ("nlpredci: SIMOPT must be 'on' or 'off'.");
  endif

  ## Fitted responses and the Jacobian of the model at the new points.
  ypred = modelfun (beta, X);
  ypred = ypred(:);
  Jnew  = nlfun_jacobian (modelfun, beta, X, eps ^ (1/3));

  ## Coefficient covariance and the error variance for observation intervals.
  if (! isempty (CovB))
    Vbeta = CovB;
    if (isempty (MSE))
      errvar = sum (resid .^ 2) / dfe;
    else
      errvar = MSE;
    endif
  elseif (! isempty (J))
    if (isempty (MSE))
      errvar = sum (resid .^ 2) / dfe;
    else
      errvar = MSE;
    endif
    Vbeta = errvar * pinv (J' * J);
  else
    error ("nlpredci: a Jacobian or a covariance matrix is required.");
  endif

  ## Variance of the fitted curve, plus the error variance for observations.
  varpred = sum ((Jnew * Vbeta) .* Jnew, 2);
  varpred = max (varpred, 0);
  if (strcmp (predopt, 'observation'))
    varpred = varpred + errvar;
  endif

  ## Critical value: Student's t (pointwise) or Scheffe (simultaneous).  The
  ## Scheffe dimension is the number of coefficients, plus one for a new
  ## observation's own error term.
  if (strcmp (simopt, 'on'))
    nsim = p + strcmp (predopt, 'observation');
    crit = sqrt (nsim * finv (1 - alpha, nsim, dfe));
  else
    crit = tinv (1 - alpha / 2, dfe);
  endif

  delta = crit * sqrt (varpred);

endfunction

%!demo
%! ## Prediction intervals for an exponential fit at three new x-values.
%! x = [1:10]';
%! y = [2.1;2.9;4.2;5.3;7.1;9.4;12.8;16.5;22.1;29.8];
%! modelfun = @(b, x) b(1) .* exp (b(2) .* x);
%! [beta, R, J] = nlinfit (x, y, modelfun, [1; 0.3]);
%! [ypred, delta] = nlpredci (modelfun, [2.5; 5.5; 8.5], beta, R, 'Jacobian', J)

%!shared modelfun, beta, R, J, CovB, MSE, xp
%! x = [1;2;3;4;5;6;7;8;9;10];
%! y = [2.1;2.9;4.2;5.3;7.1;9.4;12.8;16.5;22.1;29.8];
%! modelfun = @(b, xx) b(1) .* exp (b(2) .* xx);
%! [beta, R, J, CovB, MSE] = nlinfit (x, y, modelfun, [1; 0.3]);
%! xp = [2.5; 5.5; 8.5];

## Values verified against MATLAB's nlpredci.
%!test
%! [yp, dc] = nlpredci (modelfun, xp, beta, R, "Jacobian", J);
%! assert_equal (yp, [3.449741734; 8.158274148; 19.293455050], 1e-6);
%! assert_equal (dc, [0.120614865; 0.160334882; 0.171532624], 1e-6);
%!test
%! ## Observation (prediction) intervals exceed the curve intervals by MSE.
%! [yp, dc] = nlpredci (modelfun, xp, beta, R, "Jacobian", J);
%! [yp, dp] = nlpredci (modelfun, xp, beta, R, "Jacobian", J, ...
%!                      "PredOpt", "observation");
%! assert_equal (dp, [0.412235094; 0.425555051; 0.429899138], 1e-6);
%! assert (all (dp > dc));
%!test
%! ## Simultaneous (Scheffe) intervals are wider than the pointwise ones.
%! [yp, dc] = nlpredci (modelfun, xp, beta, R, "Jacobian", J);
%! [yp, ds] = nlpredci (modelfun, xp, beta, R, "Jacobian", J, "SimOpt", "on");
%! assert_equal (ds, [0.156197123; 0.207634833; 0.222135990], 1e-6);
%! assert (all (ds > dc));
%!test
%! ## The covariance form matches the Jacobian form.
%! [yp, dj] = nlpredci (modelfun, xp, beta, R, "Jacobian", J);
%! [yp, dv] = nlpredci (modelfun, xp, beta, R, "Covar", CovB, "MSE", MSE);
%! assert_equal (dj, dv, 1e-8);

## Test input validation
%!error<Invalid call> nlpredci (@(b, x) x, [1;2], [1;2], [1;2])
%!error<nlpredci: MODELFUN must be a function handle.> ...
%! nlpredci (1, [1;2], [1;2], [1;2;3], "Jacobian", eye (2))
%!error<nlpredci: PREDOPT must be> ...
%! nlpredci (@(b, x) x, [1;2], [1;2], [1;2;3], "Jacobian", eye (2), ...
%!           "PredOpt", "bad")
%!error<nlpredci: unknown parameter name 'foo'.> ...
%! nlpredci (@(b, x) x, [1;2], [1;2], [1;2;3], "foo", 1)
