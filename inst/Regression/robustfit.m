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
## Street-Carroll-Ruppert robust scale estimate and differs from MATLAB's by
## about 1.5%, measured, with a negligible effect on the standard errors.
##
## That difference is left in place deliberately.  The squared influence is
## averaged here over @math{n}, where the estimator as it is usually published
## averages over @math{n-p}; taking that published form moves the result
## further from MATLAB rather than closer, so MATLAB implements neither, and
## matching it would mean reproducing an undocumented variant.  The same scale
## serves @code{nlinfit}, which is why its robust @code{MSE} and @code{CovB}
## carry the same difference.
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
    tune = robusttune (wfun);
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
    s = madsigma (radj, p, y);
    w = robustwfun (radj / (s * tune), wfun);
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

  ## Robust scale and coefficient covariance.  The reported weights stay the
  ## ones the last iteration actually used, which came from the floored
  ## scale; the scale below is computed from the unfloored one, which is what
  ## MATLAB reports as mad_s and what its robust_s is built from.  Recomputing
  ## the weights here and reporting those instead graded them out of rounding
  ## noise on a fit that is already exact.
  r = y - X * b;
  radj = r .* adj;
  mad_s = madsigma (radj, p);
  [wc, psi, psip] = robustwfun (radj / (mad_s * tune), wfun);
  if (all (wc == 1))
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
%! assert_equal (b, [1.08223788958791; 1.9947584179332], 1e-8);
%! assert_equal (st.ols_s, 4.58673317428878, 1e-8);
%! assert_equal (st.mad_s, 0.219326216641263, 1e-8);
%! assert_equal (st.dfe, 8);
%! assert_equal (st.h(1), 0.345454545454545, 1e-9);
%! assert_equal (st.w(10), 0, 1e-10);
%! assert_equal (st.resid(10), -16.0298220689199, 1e-6);

%!test  # MATLAB parity: standard errors, t, p within a fraction of a percent
%! [b, st] = robustfit (X, y);
%! assert_equal (st.s, 2.45430591997485, 5e-3);
%! assert_equal (st.se, [1.67661012843903; 0.270210188642742], 2e-3);
%! assert_equal (st.t, [0.64549168064224; 7.38224723483898], 5e-3);
%! assert_equal (st.p, [0.536678473587217; 7.75017584465372e-05], 5e-3);
%! assert_equal (st.covb, [2.81102152278434, -0.401574503254905; ...
%!                   -0.401574503254905, 0.0730135460463465], 1e-2);

%!test  # weight-function coefficients match MATLAB
%! assert_equal (robustfit (X, y, "huber"), ...
%!         [1.13942072329047; 1.97894586334502], 1e-6);
%! assert_equal (robustfit (X, y, "andrews"), ...
%!         [1.08226666865865; 1.99475428605778], 1e-6);
%! assert_equal (robustfit (X, y, "cauchy"), ...
%!         [1.08721258646419; 1.99357974946387], 1e-6);
%! assert_equal (robustfit (X, y, "fair"), ...
%!         [1.25156505474148; 1.94736889966414], 1e-6);
%! assert_equal (robustfit (X, y, "logistic"), ...
%!         [1.14596301338789; 1.9773188568744], 1e-6);
%! assert_equal (robustfit (X, y, "talwar"), [1.08055555555555; 1.995], 1e-6);
%! assert_equal (robustfit (X, y, "welsch"), ...
%!         [1.08260376495793; 1.99470613677464], 1e-6);

%!test  # 'ols' weight reproduces ordinary least squares
%! assert_equal (robustfit (X, y, "ols"), regress (y, [ones(10,1), X]), 1e-10);

%!test  # custom tuning constant and const='off'
%! assert_equal (robustfit (X, y, "bisquare", 3), ...
%!         [1.08512259456964; 1.99435302345156], 1e-6);
%! [b, st] = robustfit (X, y, "bisquare", 4.685, "off");
%! assert_equal (b, 2.16460032959659, 1e-6);
%! assert_equal (st.dfe, 9);

%!test  # MATLAB parity: an exactly zero residual weighs 1, not 0
%! ## andrews and logistic are w = sin (z) / z and tanh (z) / z, both 0/0 at
%! ## the origin.  With no constant term the middle residual is exactly zero
%! ## whatever the coefficient, so this reaches the limit every time.
%! xz = [-1; 0; 1];
%! yz = [-2; 0; 3];
%! [bz, sz] = robustfit (xz, yz, "andrews", [], "off");
%! assert_equal (sz.resid(2), 0);
%! assert_equal (sz.w, [0.958241991423609; 1; 0.958241991423609], 1e-12);
%! [bz, sz] = robustfit (xz, yz, "logistic", [], "off");
%! assert_equal (sz.resid(2), 0);
%! assert_equal (sz.w, [0.907175968590982; 1; 0.907175968590982], 1e-12);

%!test  # a planted outlier is downweighted relative to OLS
%! x = (1:20)';
%! yy = 3 * x - 5;  yy(7) = yy(7) + 100;
%! b = robustfit (x, yy);
%! assert_equal (b, [-5; 3], 0.1);

%!test  # missing observations are dropped, weights padded with NaN
%! x = (1:10)';  yy = 2*x + 1;  yy(4) = NaN;
%! [b, st] = robustfit (x, yy);
%! assert_equal (b, [1; 2], 1e-8);
%! assert_equal (isnan (st.w(4)), true);
%! assert_equal (st.dfe, 7);

## A fit that is already exact keeps every weight at 1.  The scale that
## weights an iteration is floored at 1e-6 * std (Y), so residuals that are
## nothing but rounding noise cannot grade the weights out of it; R2024a
## returns ones here and used to differ from us by 0.24 on the outermost.
%!test
%! xf = [-2; -1; 0; 1; 2];
%! yf = 3 * xf;
%! [bf, stats] = robustfit (xf, yf, 'andrews', [], 'off');
%! assert_equal (bf, 3, 1e-12);
%! assert_equal (stats.w, ones (5, 1), 1e-8);

## Just below the floor the weights are graded, but only slightly, and these
## are R2024a's own numbers.
%!test
%! xf = [-2; -1; 0; 1; 2];
%! yf = 3 * xf + 1e-6 * [1; -1; 0; 1; -1];
%! [bf, stats] = robustfit (xf, yf, 'andrews', [], 'off');
%! assert_equal (stats.w', [0.997523, 0.993403, 1, 0.993403, 0.997523], 1e-6);

## Well above the floor nothing is floored and the weights are the ordinary
## ones, again R2024a's.
%!test
%! xf = [-2; -1; 0; 1; 2];
%! yf = 3 * xf + 1e-3 * [1; -1; 0; 1; -1];
%! [bf, stats] = robustfit (xf, yf, 'andrews', [], 'off');
%! assert_equal (stats.w', [0.958242, 0.868203, 1, 0.868203, 0.958242], 1e-6);

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
