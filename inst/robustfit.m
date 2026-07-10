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
## @deftypefn  {statistics} {@var{b} =} robustfit (@var{X}, @var{y})
## @deftypefnx {statistics} {@var{b} =} robustfit (@var{X}, @var{y}, @var{wfun})
## @deftypefnx {statistics} {@var{b} =} robustfit (@var{X}, @var{y}, @var{wfun}, @var{tune})
## @deftypefnx {statistics} {@var{b} =} robustfit (@var{X}, @var{y}, @var{wfun}, @var{tune}, @var{const})
## @deftypefnx {statistics} {[@var{b}, @var{stats}] =} robustfit (@dots{})
##
## Robust linear regression.
##
## @code{@var{b} = robustfit (@var{X}, @var{y})} returns the coefficient vector
## @var{b} of a linear regression of the response @var{y} on the predictors
## @var{X}, fitted by robust M-estimation (iteratively reweighted least squares)
## so that outlying observations are downweighted.  A column of ones is added to
## @var{X} by default, so @code{@var{b}(1)} is the intercept.
##
## @code{@var{b} = robustfit (@var{X}, @var{y}, @var{wfun}, @var{tune},
## @var{const})} selects the weight function @var{wfun}, its tuning constant
## @var{tune}, and whether a constant term is included.  @var{wfun} is one of
## @qcode{'bisquare'} (default), @qcode{'andrews'}, @qcode{'cauchy'},
## @qcode{'fair'}, @qcode{'huber'}, @qcode{'logistic'}, @qcode{'ols'},
## @qcode{'talwar'}, @qcode{'welsch'}, or a function handle @code{@@(r)} giving
## the weights as a function of the scaled residual.  @var{tune} defaults to the
## value that gives 95% efficiency for each weight function.  @var{const} is
## @qcode{'on'} (default) to include a constant term or @qcode{'off'} to omit.
##
## @code{[@var{b}, @var{stats}] = robustfit (@dots{})} also returns a structure
## @var{stats} with fields @code{ols_s}, @code{robust_s}, @code{mad_s},
## @code{s}, @code{se}, @code{covb}, @code{coeffcorr}, @code{t}, @code{p},
## @code{w}, @code{R}, @code{dfe}, @code{h}, and @code{resid}.  The coefficients
## and the fields @code{ols_s}, @code{mad_s}, @code{dfe}, @code{h}, @code{w},
## and @code{resid} match MATLAB.  The standard errors and quantities derived
## from them (@code{se}, @code{t}, @code{p}, @code{covb}) agree with MATLAB to
## within a small fraction of a percent; @code{robust_s} is the
## Street-Carroll-Ruppert robust scale estimate, which may differ from MATLAB's
## by a few percent (it has a negligible effect on the standard errors).
##
## @seealso{regress, fitlm}
## @end deftypefn

function [b, stats] = robustfit (X, y, wfun, tune, const)

  if (nargin < 2)
    print_usage ();
  endif
  if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
    error ("robustfit: X must be a real matrix.");
  endif
  if (! (isnumeric (y) && isreal (y) && isvector (y)))
    error ("robustfit: Y must be a real vector.");
  endif
  y = y(:);
  if (isvector (X) && rows (X) == 1 && columns (X) == numel (y))
    X = X(:);
  endif
  if (rows (X) != numel (y))
    error ("robustfit: X and Y must have the same number of observations.");
  endif
  if (nargin < 3 || isempty (wfun))
    wfun = "bisquare";
  endif
  if (ischar (wfun))
    wfun = lower (wfun);
    if (! any (strcmp (wfun, {"bisquare", "andrews", "cauchy", "fair", ...
                              "huber", "logistic", "ols", "talwar", "welsch"})))
      error ("robustfit: unrecognised weight function '%s'.", wfun);
    endif
  elseif (! is_function_handle (wfun))
    error ("robustfit: WFUN must be a name or a function handle.");
  endif
  if (nargin < 4 || isempty (tune))
    tune = default_tune (wfun);
  elseif (! (isnumeric (tune) && isscalar (tune) && isreal (tune) && tune > 0))
    error ("robustfit: TUNE must be a positive scalar.");
  endif
  if (nargin < 5 || isempty (const))
    const = "on";
  endif
  if (! (ischar (const) && any (strcmpi (const, {"on", "off"}))))
    error ("robustfit: CONST must be 'on' or 'off'.");
  endif
  addconst = strcmpi (const, "on");

  ## Drop observations with missing values.
  ok = ! (any (isnan (X), 2) | isnan (y));
  n0 = numel (y);
  X = X(ok,:);  y = y(ok);
  if (addconst)
    X = [ones(rows (X), 1), X];
  endif
  [n, p] = size (X);
  if (n <= p)
    error ("robustfit: not enough observations for the number of parameters.");
  endif

  ## Ordinary least squares start, leverages, and the OLS scale.
  [Q, R, perm] = qr (X, 0);
  b = zeros (p, 1);
  b(perm) = R \ (Q' * y);
  hlev = min (0.9999, sum (Q .^ 2, 2));
  adj = 1 ./ sqrt (1 - hlev);
  dfe = n - p;
  ols_s = norm (y - X * b) / sqrt (dfe);

  ## Iteratively reweighted least squares.
  bprev = b;
  for iter = 1:50
    r = y - X * b;
    radj = r .* adj;
    s = madsigma (radj, p);
    if (s == 0)
      s = 1;
    endif
    w = weightfun (radj / (s * tune), wfun);
    sw = sqrt (w);
    b(perm) = (X(:,perm) .* sw) \ (y .* sw);
    if (all (abs (b - bprev) <= sqrt (eps) * max (abs (b), abs (bprev))))
      break;
    endif
    bprev = b;
  endfor

  if (nargout < 2)
    return;
  endif

  ## Robust scale and coefficient covariance.
  r = y - X * b;
  radj = r .* adj;
  mad_s = madsigma (radj, p);
  [w, psi, psip] = weightfun (radj / (mad_s * tune), wfun);
  if (all (w == 1))
    robust_s = ols_s;
  else
    K = 1 + (p / n) * var (psip) / mean (psip) ^ 2;
    robust_s = tune * mad_s * sqrt (mean (psi .^ 2)) / mean (psip) * K;
  endif
  ## Combine the OLS and robust scales (larger of the two, blended by n and p).
  s = sqrt ((p ^ 2 * ols_s ^ 2 + n * robust_s ^ 2) / (n + p ^ 2));
  s = max (s, robust_s);

  covb = s ^ 2 * inv (X' * X);
  se = sqrt (diag (covb));
  se_outer = se * se';
  coeffcorr = covb ./ se_outer;
  t = b ./ se;
  pval = 2 * tcdf (-abs (t), dfe);

  stats.ols_s = ols_s;
  stats.robust_s = robust_s;
  stats.mad_s = mad_s;
  stats.s = s;
  stats.se = se;
  stats.covb = covb;
  stats.coeffcorr = coeffcorr;
  stats.t = t;
  stats.p = pval;
  stats.w = unfilter (w, ok, n0);
  stats.R = R;
  stats.dfe = dfe;
  stats.h = unfilter (hlev, ok, n0);
  stats.resid = unfilter (r, ok, n0);

endfunction

## Default tuning constant (95% efficiency) for each named weight function.
function t = default_tune (wfun)
  if (is_function_handle (wfun))
    t = 1;
    return;
  endif
  switch (wfun)
    case "andrews";      t = 1.339;
    case "bisquare";     t = 4.685;
    case "cauchy";       t = 2.385;
    case "fair";         t = 1.400;
    case "huber";        t = 1.345;
    case "logistic";     t = 1.205;
    case "ols";          t = 1;
    case "talwar";       t = 2.795;
    case "welsch";       t = 2.985;
  endswitch
endfunction

## Robust scale from the median absolute deviation of the residuals, dropping
## the P-1 smallest to account for the fitted parameters.
function s = madsigma (r, p)
  rs = sort (abs (r));
  s = median (rs(max (1, p):end)) / 0.6745;
endfunction

## Weight W (and, when requested, PSI = z.*W and its derivative PSIP) of the
## scaled residual Z for the given weight function.
function [w, psi, psip] = weightfun (z, wfun)
  if (is_function_handle (wfun))
    w = wfun (z);
    if (nargout > 1)
      d = 1e-6;
      psi = z .* w;
      psip = ((z + d) .* wfun (z + d) - (z - d) .* wfun (z - d)) / (2 * d);
    endif
    return;
  endif
  a = abs (z);
  switch (wfun)
    case "andrews"
      in = a < pi;
      w = (sin (z) ./ (z + (z == 0))) .* in;
      w(z == 0) = 1;
      if (nargout > 1); psi = sin (z) .* in;  psip = cos (z) .* in;  endif
    case "bisquare"
      in = a < 1;
      w = (1 - z .^ 2) .^ 2 .* in;
      if (nargout > 1)
        psi = z .* w;  psip = (1 - z .^ 2) .* (1 - 5 * z .^ 2) .* in;
      endif
    case "cauchy"
      w = 1 ./ (1 + z .^ 2);
      if (nargout > 1)
        psi = z .* w;  psip = (1 - z .^ 2) ./ (1 + z .^ 2) .^ 2;
      endif
    case "fair"
      w = 1 ./ (1 + a);
      if (nargout > 1); psi = z .* w;  psip = 1 ./ (1 + a) .^ 2;  endif
    case "huber"
      w = 1 ./ max (1, a);
      if (nargout > 1); psi = z .* w;  psip = double (a < 1);  endif
    case "logistic"
      th = tanh (z);
      w = th ./ (z + (z == 0));
      w(z == 0) = 1;
      if (nargout > 1); psi = th;  psip = 1 - th .^ 2;  endif
    case "ols"
      w = ones (size (z));
      if (nargout > 1); psi = z;  psip = ones (size (z));  endif
    case "talwar"
      in = a < 1;
      w = double (in);
      if (nargout > 1); psi = z .* in;  psip = double (in);  endif
    case "welsch"
      e = exp (- z .^ 2);
      w = e;
      if (nargout > 1); psi = z .* e;  psip = (1 - 2 * z .^ 2) .* e;  endif
  endswitch
endfunction

## Scatter a per-observation vector back to the original length, with NaN for
## observations that were dropped for missing values.
function out = unfilter (v, ok, n0)
  if (all (ok))
    out = v;
  else
    out = NaN (n0, 1);
    out(ok) = v;
  endif
endfunction

%!demo
%! ## Robust fit is resistant to an outlier that pulls the OLS line
%! x = (1:10)';
%! y = 2 * x + 1;
%! y(10) = 0;                       # an outlier
%! b_ols = regress (y, [ones(10,1), x]);
%! b_rob = robustfit (x, y);
%! plot (x, y, "o", x, [ones(10,1) x]*b_ols, "r-", ...
%!       x, [ones(10,1) x]*b_rob, "b-");
%! legend ("data", "OLS", "robust", "location", "northwest");

%!shared X, y
%! X = [1;2;3;4;5;6;7;8;9;10];
%! y = [3.1;5.2;6.9;9.1;11.0;12.9;15.2;17.1;19.0;5.0];

%!test  # MATLAB parity: bisquare coefficients and exact stats fields
%! [b, st] = robustfit (X, y);
%! assert (b, [1.08223788958791; 1.9947584179332], 1e-8);
%! assert (st.ols_s, 4.58673317428878, 1e-8);
%! assert (st.mad_s, 0.219326216641263, 1e-8);
%! assert (st.dfe, 8);
%! assert (st.h(1), 0.345454545454545, 1e-9);
%! assert (st.w(10), 0, 1e-10);
%! assert (st.resid(10), -16.0298220689199, 1e-6);

%!test  # MATLAB parity: standard errors, t, p within a fraction of a percent
%! [b, st] = robustfit (X, y);
%! assert (st.s, 2.45430591997485, 5e-3);
%! assert (st.se, [1.67661012843903; 0.270210188642742], 2e-3);
%! assert (st.t, [0.64549168064224; 7.38224723483898], 5e-3);
%! assert (st.p, [0.536678473587217; 7.75017584465372e-05], 5e-3);
%! assert (st.covb, [2.81102152278434, -0.401574503254905; ...
%!                   -0.401574503254905, 0.0730135460463465], 1e-2);

%!test  # weight-function coefficients match MATLAB
%! assert (robustfit (X, y, "huber"), ...
%!         [1.13942072329047; 1.97894586334502], 1e-6);
%! assert (robustfit (X, y, "andrews"), ...
%!         [1.08226666865865; 1.99475428605778], 1e-6);
%! assert (robustfit (X, y, "cauchy"), ...
%!         [1.08721258646419; 1.99357974946387], 1e-6);
%! assert (robustfit (X, y, "fair"), ...
%!         [1.25156505474148; 1.94736889966414], 1e-6);
%! assert (robustfit (X, y, "logistic"), ...
%!         [1.14596301338789; 1.9773188568744], 1e-6);
%! assert (robustfit (X, y, "talwar"), [1.08055555555555; 1.995], 1e-6);
%! assert (robustfit (X, y, "welsch"), ...
%!         [1.08260376495793; 1.99470613677464], 1e-6);

%!test  # 'ols' weight reproduces ordinary least squares
%! assert (robustfit (X, y, "ols"), regress (y, [ones(10,1), X]), 1e-10);

%!test  # custom tuning constant and const='off'
%! assert (robustfit (X, y, "bisquare", 3), ...
%!         [1.08512259456964; 1.99435302345156], 1e-6);
%! [b, st] = robustfit (X, y, "bisquare", 4.685, "off");
%! assert (b, 2.16460032959659, 1e-6);
%! assert (st.dfe, 9);

%!test  # a planted outlier is downweighted relative to OLS
%! x = (1:20)';
%! yy = 3 * x - 5;  yy(7) = yy(7) + 100;
%! b = robustfit (x, yy);
%! assert (b, [-5; 3], 0.1);

%!test  # missing observations are dropped, weights padded with NaN
%! x = (1:10)';  yy = 2*x + 1;  yy(4) = NaN;
%! [b, st] = robustfit (x, yy);
%! assert (b, [1; 2], 1e-8);
%! assert (isnan (st.w(4)));
%! assert (st.dfe, 7);

## Test input validation
%!error <Invalid call to robustfit> robustfit (1)
%!error <robustfit: X and Y must have the same number of observations.> ...
%! robustfit ([1;2;3], [1;2])
%!error <robustfit: unrecognised weight function 'foo'.> ...
%! robustfit (X, y, "foo")
%!error <robustfit: TUNE must be a positive scalar.> ...
%! robustfit (X, y, "huber", -1)
%!error <robustfit: CONST must be 'on' or 'off'.> ...
%! robustfit (X, y, "huber", 1.345, "maybe")
