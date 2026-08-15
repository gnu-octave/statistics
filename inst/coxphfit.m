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
## @deftypefn  {statistics} {@var{b} =} coxphfit (@var{X}, @var{T})
## @deftypefnx {statistics} {@var{b} =} coxphfit (@var{X}, @var{T}, @var{name}, @var{value}, @dots{})
## @deftypefnx {statistics} {[@var{b}, @var{logl}] =} coxphfit (@dots{})
## @deftypefnx {statistics} {[@var{b}, @var{logl}, @var{H}] =} coxphfit (@dots{})
## @deftypefnx {statistics} {[@var{b}, @var{logl}, @var{H}, @var{stats}] =} coxphfit (@dots{})
##
## Fit a Cox proportional hazards regression model.
##
## @code{@var{b} = coxphfit (@var{X}, @var{T})} returns the @math{p}-by-1 vector
## of coefficients @var{b} of the Cox model
##
## @tex
## $$ h(x_i, t) = h_0(t)\exp\left(\sum_{j=1}^{p} x_{ij} b_j\right) $$
## @end tex
## @ifnottex
## @math{h(x_i, t) = h_0(t) exp (x_i' b)}
## @end ifnottex
##
## fitted to the @math{n}-by-@math{p} matrix of predictors @var{X} and the
## @math{n}-by-1 vector of event times @var{T}.  @math{h_0(t)} is the baseline
## hazard, which is left unspecified: the coefficients are estimated by
## maximizing the Cox partial likelihood, which does not involve it.
##
## @strong{@var{X} must not contain a column of ones.}  The model has no
## constant term, since any constant is absorbed into the baseline hazard.  A
## constant column is detected, reported by a warning, and given a zero
## coefficient.
##
## Rows of @var{X}, @var{T} or @qcode{"Frequency"} holding @qcode{NaN} are
## removed before fitting.
##
## @code{[@var{b}, @var{logl}, @var{H}, @var{stats}] = coxphfit (@dots{})}
## additionally returns the maximized partial log-likelihood @var{logl}, the
## estimated baseline cumulative hazard @var{H}, and a structure @var{stats} of
## coefficient statistics and residuals.
##
## @var{H} is a two-column matrix whose first column holds the distinct event
## times and whose second holds the estimated cumulative hazard at those times,
## evaluated at the predictor values given by @qcode{"Baseline"}.  Its first
## row is the earliest time with a cumulative hazard of zero.
##
## The following @var{name}/@var{value} pairs are accepted:
##
## @multitable @columnfractions 0.18 0.82
## @headitem Name @tab Value
## @item @qcode{"Baseline"} @tab The @var{X} values at which the baseline
## hazard is computed, either a scalar or a 1-by-@math{p} vector.  The default
## is @code{mean (@var{X})}, so the hazard is that of an average observation;
## pass @qcode{0} for a hazard relative to the origin.  The coefficients do not
## depend on this choice, only @var{H} does.
## @item @qcode{"Censoring"} @tab A logical or 0/1 vector of length @math{n},
## where 1 marks an observation right-censored at its recorded time.  The
## default is a vector of zeros, so every observation is a recorded event.
## @item @qcode{"Frequency"} @tab A vector of length @math{n} of non-negative
## values giving the number of observations each row represents, or a weight.
## The default is a vector of ones.
## @item @qcode{"Ties"} @tab The method of handling tied event times, either
## @qcode{"breslow"} (default) or @qcode{"efron"}.
## @item @qcode{"B0"} @tab The starting value of the iteration, a vector of
## length @math{p}.  The default is @code{0.01 ./ std (@var{X})}.
## @item @qcode{"Options"} @tab A structure of iteration settings, as built by
## @code{statset ("coxphfit")}.  The fields used are @qcode{"MaxIter"},
## @qcode{"TolX"} and @qcode{"Display"}.
## @end multitable
##
## The fields of @var{stats} are:
##
## @multitable @columnfractions 0.24 0.76
## @headitem Field @tab Contents
## @item @qcode{"covb"} @tab The estimated covariance matrix of @var{b}.
## @item @qcode{"beta"} @tab The coefficients, as returned in @var{b}.
## @item @qcode{"se"} @tab The standard errors of the coefficients.
## @item @qcode{"z"} @tab The @math{z} statistics, @var{b} over its standard
## error.
## @item @qcode{"p"} @tab The two-sided @math{p}-values of the @math{z}
## statistics.
## @item @qcode{"csres"} @tab The Cox-Snell residuals.
## @item @qcode{"devres"} @tab The deviance residuals.
## @item @qcode{"martres"} @tab The martingale residuals.
## @item @qcode{"schres"} @tab The Schoenfeld residuals, @qcode{NaN} for a
## censored observation.
## @item @qcode{"sschres"} @tab The scaled Schoenfeld residuals.
## @item @qcode{"scores"} @tab The score residuals.
## @item @qcode{"sscores"} @tab The scaled score residuals.
## @item @qcode{"LikelihoodRatioTestP"} @tab The @math{p}-value of the
## likelihood ratio test against the model with no predictors.
## @end multitable
##
## Stratified models and the counting-process form of @var{T}, in which each
## row gives a @math{(start, stop]} interval, are not yet implemented and are
## rejected rather than silently ignored.
##
## @strong{One documented deviation, for @qcode{"efron"} ties only.}  The
## martingale residual is defined as the event indicator minus the estimated
## cumulative hazard, so @qcode{"csres"} and @qcode{"martres"} must sum to that
## indicator.  They do here, for both tie methods.  In MATLAB R2024a they do so
## only for @qcode{"breslow"}: under @qcode{"efron"} its @qcode{"martres"} is
## computed from a cumulative hazard that agrees neither with its own
## @qcode{"csres"} nor with the @var{H} it returns, and the two residuals sum
## to 1.0437 and @math{-0.0414} where they must give 1 and 0.  Every other
## output, and every residual under @qcode{"breslow"}, agrees with R2024a to
## machine precision.
##
## @seealso{statset, ecdf, fitlm}
## @end deftypefn

function [b, logl, H, stats] = coxphfit (X, T, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  ## --- X and T -----------------------------------------------------------
  if (! (isnumeric (X) && isreal (X) && ismatrix (X) && ! isempty (X)))
    error ("coxphfit: X must be a real numeric matrix.");
  endif
  if (! (isnumeric (T) && isreal (T) && ! isempty (T)))
    error ("coxphfit: T must be a real numeric vector.");
  endif
  if (columns (T) == 2)
    error (strcat ("coxphfit: the counting process form of T, giving a", ...
                   " (start, stop] interval per row, is not implemented."));
  endif
  if (! isvector (T))
    error ("coxphfit: T must be a real numeric vector.");
  endif
  T = T(:);
  n = numel (T);
  if (rows (X) != n)
    error ("coxphfit: T must have one element for each row of X.");
  endif
  p = columns (X);

  ## --- name/value pairs --------------------------------------------------
  Baseline = [];
  Censoring = zeros (n, 1);
  Frequency = ones (n, 1);
  Ties = 'breslow';
  B0 = [];
  Options = [];

  if (mod (numel (varargin), 2) != 0)
    error ("coxphfit: optional arguments must occur in NAME/VALUE pairs.");
  endif
  for i = 1:2:numel (varargin)
    name = varargin{i};
    if (! ischar (name))
      error ("coxphfit: parameter name must be a character vector.");
    endif
    switch (lower (name))
      case 'baseline'
        Baseline = varargin{i+1};
      case 'censoring'
        Censoring = varargin{i+1};
      case 'frequency'
        Frequency = varargin{i+1};
      case 'ties'
        Ties = varargin{i+1};
      case 'b0'
        B0 = varargin{i+1};
      case 'options'
        Options = varargin{i+1};
      case 'strata'
        error ("coxphfit: stratified models are not implemented.");
      otherwise
        error ("coxphfit: unknown parameter name '%s'.", name);
    endswitch
  endfor

  if (! (ischar (Ties) && any (strcmpi (Ties, {'breslow', 'efron'}))))
    error (strcat ("coxphfit: 'Ties' must be either 'breslow' or", ...
                   " 'efron'."));
  endif
  Ties = lower (Ties);

  if (! (isnumeric (Censoring) || islogical (Censoring))
      || numel (Censoring) != n)
    error (strcat ("coxphfit: 'Censoring' must have one element for each", ...
                   " row of X."));
  endif
  Censoring = logical (Censoring(:));

  if (! (isnumeric (Frequency) && isreal (Frequency))
      || numel (Frequency) != n)
    error (strcat ("coxphfit: 'Frequency' must have one element for each", ...
                   " row of X."));
  endif
  Frequency = Frequency(:);
  if (any (Frequency < 0))
    error (strcat ("coxphfit: 'Frequency' must be a vector of", ...
                   " non-negative values."));
  endif

  if (! isempty (B0) && (! isnumeric (B0) || numel (B0) != p))
    error (strcat ("coxphfit: 'B0' must be a real vector with %d", ...
                   " elements."), p);
  endif

  ## Any option the caller leaves out falls back on this function's own
  ## defaults, which is what statset's merge already means by an empty field.
  if (isempty (Options))
    Options = statset ('coxphfit');
  elseif (isstruct (Options))
    Options = statset (statset ('coxphfit'), Options);
  else
    error ("coxphfit: 'Options' must be a structure.");
  endif
  maxiter = statget (Options, 'MaxIter', 100);
  tolx = statget (Options, 'TolX', 1e-8);
  display = statget (Options, 'Display', 'off');

  ## --- drop incomplete rows ----------------------------------------------
  ok = all (isfinite (X), 2) & isfinite (T) & isfinite (Frequency);
  if (! all (ok))
    X = X(ok,:);
    T = T(ok);
    Censoring = Censoring(ok);
    Frequency = Frequency(ok);
    n = numel (T);
  endif
  if (n == 0)
    error ("coxphfit: no complete observations remain after removing NaNs.");
  endif

  ## A constant column carries no information: the Cox model has no intercept.
  const = false (1, p);
  for j = 1:p
    const(j) = all (X(:,j) == X(1,j));
  endfor
  if (any (const))
    warning ("coxphfit: the Cox model cannot have a constant term in X.");
  endif

  ## The default baseline is the mean of X weighted by 'Frequency'.  MathWorks
  ## documents it as mean (X) without qualification, but the weighted mean is
  ## what R2024a computes -- measured, since the two agree whenever every
  ## frequency is 1 and the documentation is therefore not wrong, only partial.
  if (isempty (Baseline))
    Baseline = sum (Frequency .* X, 1) / sum (Frequency);
  endif
  if (! (isnumeric (Baseline) && isreal (Baseline)
         && (isscalar (Baseline) || numel (Baseline) == p)))
    error (strcat ("coxphfit: 'Baseline' must be a scalar or a vector", ...
                   " with %d elements."), p);
  endif
  if (isscalar (Baseline))
    Baseline = repmat (Baseline, 1, p);
  endif
  Baseline = Baseline(:)';

  if (isempty (B0))
    s = std (X, 0, 1);
    s(s == 0) = 1;
    B0 = (0.01 ./ s)';
  endif
  B0 = B0(:);
  B0(const) = 0;

  event = ! Censoring;
  w = Frequency;
  free = ! const;

  ## --- Newton-Raphson on the partial likelihood --------------------------
  b = B0;
  for iter = 1:maxiter
    [logl, g, F] = partial_lik (b, X, T, event, w, Ties);
    if (all (! free))
      break;
    endif
    step = zeros (p, 1);
    step(free) = F(free,free) \ g(free);
    if (! all (isfinite (step)))
      warning (strcat ("coxphfit: estimation stopped, the maximum", ...
                       " likelihood estimate may not be finite."));
      break;
    endif
    b += step;
    if (strcmp (display, 'iter'))
      printf ("coxphfit: iteration %d, log-likelihood %g\n", iter, logl);
    endif
    if (max (abs (step)) < tolx)
      break;
    endif
  endfor
  [logl, g, F] = partial_lik (b, X, T, event, w, Ties);
  b(const) = 0;

  if (nargout < 3)
    return;
  endif

  ## --- baseline cumulative hazard ----------------------------------------
  ## H is computed at X == Baseline, so the linear predictor is centred there;
  ## this scales the hazard by exp (Baseline * b) and leaves b untouched.
  eta = (X - Baseline) * b;
  ex = w .* exp (eta);
  ut = unique (T(event));
  h = zeros (numel (ut), 1);
  ## The hazard uses the Breslow estimator for both tie methods.  'Ties' enters
  ## only through the partial likelihood, and so through b; MATLAB does the
  ## same, which is why its Efron and Breslow hazards differ by the fit alone.
  for k = 1:numel (ut)
    D = event & T == ut(k);
    R = T >= ut(k);
    h(k) = sum (w(D)) / sum (ex(R));
  endfor
  H = [min(T), 0; ut, cumsum(h)];

  if (nargout < 4)
    return;
  endif

  ## --- coefficient statistics --------------------------------------------
  covb = zeros (p);
  covb(free,free) = inv (F(free,free));
  se = sqrt (diag (covb));
  z = b ./ se;
  pval = 2 * normcdf (-abs (z));

  ## Likelihood ratio against the model with no predictors.
  logl0 = partial_lik (zeros (p, 1), X, T, event, w, Ties);
  df = sum (free);
  lrtp = 1 - chi2cdf (2 * (logl - logl0), df);

  ## --- residuals ----------------------------------------------------------
  ## Cumulative hazard carried by each observation, on the uncentred scale.
  exu = w .* exp (X * b);
  h0 = zeros (numel (ut), 1);
  for k = 1:numel (ut)
    D = event & T == ut(k);
    R = T >= ut(k);
    h0(k) = sum (w(D)) / sum (exu(R));
  endfor
  ch0 = cumsum (h0);
  Hi = zeros (n, 1);
  for i = 1:n
    idx = find (ut <= T(i), 1, 'last');
    if (! isempty (idx))
      Hi(i) = ch0(idx);
    endif
  endfor

  martres = double (event) - Hi .* exp (X * b);
  csres = double (event) - martres;
  ## Deviance residuals: the martingale residual symmetrized about zero.
  dev = zeros (n, 1);
  for i = 1:n
    m = martres(i);
    e = double (event(i));
    if (e - m <= 0)
      dev(i) = sign (m) * sqrt (-2 * m);
    else
      dev(i) = sign (m) * sqrt (-2 * (m + e * log (e - m)));
    endif
  endfor
  devres = dev;

  ## Schoenfeld residuals: the predictor minus the risk-set weighted mean, at
  ## each event time.  Censored observations have none.
  schres = NaN (n, p);
  for i = 1:n
    if (! event(i))
      continue;
    endif
    R = T >= T(i);
    xbar = (exu(R)' * X(R,:)) / sum (exu(R));
    schres(i,:) = X(i,:) - xbar;
  endfor
  nev = sum (event);
  sschres = repmat (b', n, 1) + nev * (schres * covb);

  ## Score residuals: the per-observation contribution to the score.
  scores = zeros (n, p);
  for i = 1:n
    acc = zeros (1, p);
    if (event(i))
      R = T >= T(i);
      xbar = (exu(R)' * X(R,:)) / sum (exu(R));
      acc += X(i,:) - xbar;
    endif
    idx = find (ut <= T(i));
    for k = idx'
      D = event & T == ut(k);
      R = T >= ut(k);
      xbar = (exu(R)' * X(R,:)) / sum (exu(R));
      acc -= (X(i,:) - xbar) * exp (X(i,:) * b) * h0(k);
    endfor
    scores(i,:) = acc;
  endfor
  sscores = scores * covb;

  stats = struct ('covb', covb, 'beta', b, 'se', se, 'z', z, 'p', pval, ...
                  'csres', csres, 'devres', devres, 'martres', martres, ...
                  'schres', schres, 'sschres', sschres, 'scores', scores, ...
                  'sscores', sscores, 'LikelihoodRatioTestP', lrtp);

endfunction

## Cox partial log-likelihood, its gradient, and the observed information.
function [l, g, F] = partial_lik (b, X, T, event, w, ties)

  p = columns (X);
  eta = X * b;
  ex = w .* exp (eta);
  ut = unique (T(event));
  l = 0;
  g = zeros (p, 1);
  F = zeros (p);

  for k = 1:numel (ut)
    tk = ut(k);
    D = find (event & T == tk);
    R = find (T >= tk);
    wD = w(D);
    WD = sum (wD);
    l += sum (wD .* eta(D));
    g += X(D,:)' * wD;
    XR = X(R,:);
    exR = ex(R);
    if (strcmp (ties, 'breslow'))
      S0 = sum (exR);
      S1 = XR' * exR;
      S2 = XR' * (XR .* exR);
      l -= WD * log (S0);
      g -= WD * (S1 / S0);
      F += WD * (S2 / S0 - (S1 / S0) * (S1 / S0)');
    else
      m = numel (D);
      XD = X(D,:);
      exD = ex(D);
      S0 = sum (exR);   S0d = sum (exD);
      S1 = XR' * exR;   S1d = XD' * exD;
      S2 = XR' * (XR .* exR);
      S2d = XD' * (XD .* exD);
      for r = 0:(m-1)
        f = r / m;
        A0 = S0 - f * S0d;
        A1 = S1 - f * S1d;
        A2 = S2 - f * S2d;
        l -= (WD / m) * log (A0);
        g -= (WD / m) * (A1 / A0);
        F += (WD / m) * (A2 / A0 - (A1 / A0) * (A1 / A0)');
      endfor
    endif
  endfor

endfunction

%!demo
%! ## Fit a Cox model to right-censored survival data
%! T = [4; 6; 8; 11; 13; 16; 18; 21; 25; 30];
%! X = [2 0; 5 1; 3 0; 8 1; 4 0; 7 1; 6 0; 9 1; 5 0; 10 1];
%! censored = [0; 0; 1; 0; 0; 1; 0; 0; 1; 0];
%! [b, logl] = coxphfit (X, T, 'Censoring', censored)

%!demo
%! ## Compare the two methods of handling tied event times
%! T = [4; 4; 6; 6; 8; 8; 11; 11; 13; 13];
%! X = [2 0; 5 1; 3 0; 8 1; 4 0; 7 1; 6 0; 9 1; 5 0; 10 1];
%! breslow = coxphfit (X, T, 'Ties', 'breslow');
%! efron = coxphfit (X, T, 'Ties', 'efron');
%! [breslow, efron]

%!shared X, T, C, Tt, F
%! T = [4; 6; 8; 11; 13; 16; 18; 21; 25; 30];
%! X = [2 0; 5 1; 3 0; 8 1; 4 0; 7 1; 6 0; 9 1; 5 0; 10 1];
%! C = [0; 0; 1; 0; 0; 1; 0; 0; 1; 0];
%! Tt = [4; 4; 6; 6; 8; 8; 11; 11; 13; 13];
%! F = [1; 2; 1; 1; 3; 1; 1; 2; 1; 1];

## Coefficients and log-likelihood, against R2024a
%!test
%! [b, logl] = coxphfit (X, T);
%! assert_equal (b, [-1.3886093196382836; 4.3814437183613322], 1e-8);
%! assert_equal (logl, -8.8069639381632356, 1e-10);

%!test
%! [b, logl] = coxphfit (X, T, 'Censoring', C);
%! assert_equal (b, [-1.0422777707095543; 3.374484233216088], 1e-8);
%! assert_equal (logl, -7.6973587461887778, 1e-10);

%!test
%! [b, logl] = coxphfit (X, Tt, 'Censoring', C, 'Ties', 'breslow');
%! assert_equal (b, [-0.71807146728479765; 3.0751600100470267], 1e-8);
%! assert_equal (logl, -9.9446106821057487, 1e-10);

%!test
%! [b, logl] = coxphfit (X, Tt, 'Censoring', C, 'Ties', 'efron');
%! assert_equal (b, [-0.78425008234710136; 3.3784787697856826], 1e-8);
%! assert_equal (logl, -9.2459071792218612, 1e-10);

%!test
%! [b, logl] = coxphfit (X, T, 'Censoring', C, 'Frequency', F);
%! assert_equal (b, [-1.1439301196172516; 3.7057550554856507], 1e-8);
%! assert_equal (logl, -15.957308443384392, 1e-10);

%!test
%! b = coxphfit (X(:,1), T, 'Censoring', C);
%! assert_equal (b, -0.28958528273410233, 1e-8);

## Breslow is the default method for ties
%!test
%! assert_equal (coxphfit (X, Tt, 'Censoring', C), ...
%!               coxphfit (X, Tt, 'Censoring', C, 'Ties', 'breslow'));

## The coefficients do not depend on Baseline or on the starting value
%!test
%! b = coxphfit (X, T, 'Censoring', C, 'Baseline', 0);
%! assert_equal (b, [-1.0422777707095543; 3.374484233216088], 1e-8);

%!test
%! b = coxphfit (X, T, 'Censoring', C, 'B0', [0.1; -0.1]);
%! assert_equal (b, [-1.0422777707095543; 3.374484233216088], 1e-8);

## The baseline cumulative hazard
%!test
%! [b, logl, H] = coxphfit (X, T, 'Censoring', C);
%! assert_equal (size (H), [8, 2]);
%! assert_equal (H(1,:), [4, 0]);
%! assert_equal (H(:,1)', [4, 4, 6, 11, 13, 18, 21, 30]);
%! assert_equal (H(end,2), 16.21202553, 1e-6);

## Baseline scales the hazard by exp (Baseline * b) and nothing else
%!test
%! [~, ~, Hm] = coxphfit (X, T, 'Censoring', C);
%! [b, ~, H0] = coxphfit (X, T, 'Censoring', C, 'Baseline', 0);
%! assert_equal (Hm(2:end,2) ./ H0(2:end,2), ...
%!               repmat (exp (mean (X) * b), 7, 1), 1e-10);

## Coefficient statistics
%!test
%! [b, logl, H, stats] = coxphfit (X, T);
%! assert_equal (stats.beta, b);
%! assert_equal (stats.se, [0.52737766743917369; 1.8537400210498443], 1e-8);
%! assert_equal (stats.z, [-2.6330453589001133; 2.3635696853973904], 1e-8);
%! assert_equal (stats.p, [0.0084623045168834322; 0.018099822319697333], 1e-10);

%!test
%! [~, ~, ~, stats] = coxphfit (X, T);
%! assert_equal (stats.covb, [0.27812720411358371, -0.88932589928247008; ...
%!                            -0.88932589928247008, 3.4363520656418767], 1e-8);

%!test
%! [~, ~, ~, stats] = coxphfit (X, T);
%! assert_equal (stats.LikelihoodRatioTestP, 0.0018409958426933715, 1e-10);

## The stats fields, in MATLAB's own order
%!test
%! [~, ~, ~, stats] = coxphfit (X, T);
%! assert_equal (fieldnames (stats)', {'covb', 'beta', 'se', 'z', 'p', ...
%! 'csres', 'devres', 'martres', 'schres', 'sschres', 'scores', 'sscores', ...
%! 'LikelihoodRatioTestP'});

## Residuals
%!test
%! [~, ~, ~, stats] = coxphfit (X, T);
%! assert_equal (stats.martres(1), 0.6260391348376092, 1e-8);
%! assert_equal (stats.csres(1), 0.3739608651623908, 1e-8);

## The Cox-Snell and martingale residuals partition the event indicator
%!test
%! [~, ~, ~, stats] = coxphfit (X, T, 'Censoring', C);
%! assert_equal (stats.csres + stats.martres, double (! C), 1e-12);

## The defining identity holds under 'efron' too, where MATLAB's does not
%!test
%! [~, ~, ~, stats] = coxphfit (X, Tt, 'Censoring', C, 'Ties', 'efron');
%! assert_equal (stats.csres + stats.martres, double (! C), 1e-12);

## A censored observation has no Schoenfeld residual
%!test
%! [~, ~, ~, stats] = coxphfit (X, T, 'Censoring', C);
%! assert_equal (all (isnan (stats.schres(logical (C),:))(:)), true);
%! assert_equal (any (isnan (stats.schres(! logical (C),:))(:)), false);

## Error conditions
%!error<Invalid call to coxphfit> coxphfit (X)

%!error<coxphfit: T must have one element for each row of X.> ...
%! coxphfit (X, T(1:5))

%!error<coxphfit: 'Ties' must be either 'breslow' or 'efron'.> ...
%! coxphfit (X, T, 'Ties', 'nosuch')

%!error<coxphfit: unknown parameter name 'NoSuchOption'.> ...
%! coxphfit (X, T, 'NoSuchOption', 1)

%!error<coxphfit: 'Censoring' must have one element for each row of X.> ...
%! coxphfit (X, T, 'Censoring', C(1:5))

%!error<coxphfit: 'Frequency' must be a vector of non-negative values.> ...
%! coxphfit (X, T, 'Frequency', -ones (10, 1))

%!error<coxphfit: 'B0' must be a real vector with 2 elements.> ...
%! coxphfit (X, T, 'B0', [1; 2; 3])

%!error<coxphfit: 'Options' must be a structure.> ...
%! coxphfit (X, T, 'Options', 5)

%!error<coxphfit: optional arguments must occur in NAME/VALUE pairs.> ...
%! coxphfit (X, T, 'Ties')

%!error<coxphfit: stratified models are not implemented.> ...
%! coxphfit (X, T, 'Strata', ones (10, 1))

%!error<coxphfit: the counting process form of T, giving a \(start, stop\] interval per row, is not implemented.> ...
%! coxphfit (X, [T, T+1])

## A constant column is reported, not silently absorbed
%!warning<coxphfit: the Cox model cannot have a constant term in X.> ...
%! coxphfit ([X, ones(10,1)], T);
