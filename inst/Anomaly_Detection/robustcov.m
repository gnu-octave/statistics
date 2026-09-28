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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{sig} =} robustcov (@var{X})
## @deftypefnx {statistics} {[@var{sig}, @var{mu}] =} robustcov (@var{X})
## @deftypefnx {statistics} {[@var{sig}, @var{mu}, @var{mah}] =} robustcov (@var{X})
## @deftypefnx {statistics} {[@var{sig}, @var{mu}, @var{mah}, @var{outliers}] =} robustcov (@var{X})
## @deftypefnx {statistics} {[@var{sig}, @var{mu}, @var{mah}, @var{outliers}, @var{s}] =} robustcov (@var{X})
## @deftypefnx {statistics} {[@dots{}] =} robustcov (@dots{}, @var{name}, @var{value})
##
## Robust multivariate covariance and mean estimate.
##
## @code{@var{sig} = robustcov (@var{X})} returns a robust estimate @var{sig} of
## the covariance matrix of the @math{N*P} data matrix @var{X}, computed so that
## it is not distorted by outlying observations.  Rows of @var{X} are
## observations and columns are variables.  Rows containing @code{NaN} values
## are removed.
##
## @code{[@var{sig}, @var{mu}, @var{mah}, @var{outliers}, @var{s}] = robustcov
## (@dots{})} also returns the robust mean @var{mu} (@math{1*P}), the robust
## Mahalanobis distances @var{mah} (@math{N*1}) of each observation from the
## estimated distribution, a logical vector @var{outliers} (@math{N*1}) flagging
## observations whose distance exceeds @code{sqrt (chi2inv (0.975, @var{P}))},
## and a structure @var{s} holding the estimate metadata.
##
## Additional parameters can be specified by @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.22 0.76
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Method'} @tab the estimator, either @qcode{'fmcd'} (default,
## the Fast Minimum Covariance Determinant algorithm) or @qcode{'ogk'} (the
## Orthogonalized Gnanadesikan-Kettenring estimator).  @qcode{'olivehawkins'} is
## not implemented.
##
## @item @qcode{'OutlierFraction'} @tab the maximum fraction of outliers, a
## scalar in @math{[0, 0.5]} (default 0.5), used to set the size of the
## elemental subsets in @qcode{'fmcd'}.
##
## @item @qcode{'NumTrials'} @tab the number of random elemental subsets drawn
## by @qcode{'fmcd'}, a positive integer (default 500).
##
## @item @qcode{'BiasCorrection'} @tab a logical scalar (default @code{true})
## that applies the small-sample bias correction to the @qcode{'fmcd'} estimate.
##
## @item @qcode{'NumOGKIterations'} @tab the number of orthogonalization
## iterations for @qcode{'ogk'}, a positive integer (default 2).
##
## @item @qcode{'UnivariateEstimator'} @tab the robust univariate
## location/scale estimator used by @qcode{'ogk'}, either @qcode{'tauscale'}
## (default) or @qcode{'qn'}.
## @end multitable
##
## @strong{Note on reproducibility.} @qcode{'fmcd'} draws random subsets, so its
## exact estimate depends on the random number generator and is not identical to
## MATLAB's on data where the optimal subset is ambiguous; on well-separated
## data both converge to the same estimate.  For @qcode{'fmcd'} with
## @qcode{'BiasCorrection'} enabled, the small-sample factor uses the published
## Pison-Van Aelst-Willems asymptotic formula, which differs from MATLAB's
## tabulated simulation values by up to about 1.6% for very small samples.
##
## @seealso{mahal, cov, mad, dbscan}
## @end deftypefn

function [sig, mu, mah, outliers, s] = robustcov (X, varargin)

  ## Check number of input arguments
  if (nargin < 1)
    error ("robustcov: too few input arguments.");
  endif

  ## Validate X and drop rows with NaN values
  if (! isnumeric (X) || ! isreal (X) || ndims (X) != 2 || isempty (X))
    error ("robustcov: X must be a nonempty real numeric matrix.");
  endif
  X = X(! any (isnan (X), 2), :);
  [n, p] = size (X);

  ## Defaults
  method         = "fmcd";
  outlierfrac    = 0.5;
  numtrials      = 500;
  biascorrection = true;
  numogkiter     = 2;
  univestimator  = "tauscale";

  ## Parse Name-Value pairs
  if (mod (numel (varargin), 2) != 0)
    error ("robustcov: each NAME must be followed by a VALUE.");
  endif
  while (numel (varargin) > 0)
    name = varargin{1};
    val  = varargin{2};
    if (! ischar (name))
      error ("robustcov: optional argument names must be strings.");
    endif
    switch (tolower (name))
      case "method"
        method = tolower (val);
      case "outlierfraction"
        outlierfrac = val;
      case "numtrials"
        numtrials = val;
      case "biascorrection"
        biascorrection = logical (val);
      case "numogkiterations"
        numogkiter = val;
      case "univariateestimator"
        univestimator = tolower (val);
      otherwise
        error ("robustcov: unknown parameter name '%s'.", name);
    endswitch
    varargin(1:2) = [];
  endwhile

  ## Validate options
  if (strcmp (method, "olivehawkins"))
    error ("robustcov: the 'olivehawkins' method is not implemented.");
  endif
  if (! any (strcmp (method, {"fmcd", "ogk"})))
    error ("robustcov: METHOD must be 'fmcd' or 'ogk'.");
  endif
  if (! isscalar (outlierfrac) || ! isnumeric (outlierfrac)
      || outlierfrac < 0 || outlierfrac > 0.5)
    error ("robustcov: OUTLIERFRACTION must be a scalar in [0, 0.5].");
  endif
  if (! any (strcmp (univestimator, {"tauscale", "qn"})))
    error ("robustcov: UNIVARIATEESTIMATOR must be 'tauscale' or 'qn'.");
  endif
  if (n <= p)
    error ("robustcov: X must have more rows than columns.");
  endif

  ## Robust cutoff for reweighting and for flagging outliers.
  cutoff = chi2inv (0.975, p);

  if (strcmp (method, "fmcd"))
    ## Raw Fast-MCD subset and its robust distances.
    h = floor (2 * floor ((n + p + 1) / 2) - n ...
               + 2 * (n - floor ((n + p + 1) / 2)) * (1 - outlierfrac));
    h = max (min (h, n), p + 1);
    [Traw, Sraw] = fastmcd (X, h, numtrials);
    Sraw = mcdcons (p, h / n) * mcdcnp2raw (p, n) * Sraw;
    ## Reweighting: keep observations within the 0.975 cutoff.
    d2 = mahal2 (X, Traw, Sraw);
    keep = d2 <= cutoff;
    mu = mean (X(keep, :));
    C  = cov (X(keep, :));
    hk = sum (keep);
    fac = mcdcons (p, hk / n);
    if (biascorrection)
      fac *= mcdcnp2rew (p, n);
    endif
    sig = fac * C;
  else
    ## Orthogonalized Gnanadesikan-Kettenring, then MCD-style reweighting.
    [Togk, Sogk] = ogk (X, numogkiter, univestimator);
    d2 = mahal2 (X, Togk, Sogk);
    ## Scale so that the median squared distance matches the chi-square median,
    ## making the reweighting cutoff meaningful.  Maronna and Zamar reweight at
    ## the 0.9 quantile, not the 0.975 used to flag outliers.
    scale = median (d2) / chi2inv (0.5, p);
    d2 = d2 / scale;
    keep = d2 <= chi2inv (0.9, p);
    mu = mean (X(keep, :));
    sig = cov (X(keep, :), 1);          # 1/N normalization, as MATLAB's OGK
  endif

  ## Final robust Mahalanobis distances and outlier flags.
  md2 = mahal2 (X, mu, sig);
  mah = sqrt (md2);
  outliers = md2 > cutoff;

  if (nargout > 4)
    if (strcmp (method, "fmcd"))
      s = struct ("BiasCorrection", biascorrection, "NumTrials", numtrials, ...
                  "OutlierFraction", outlierfrac, "Mu", mu, "Sigma", sig, ...
                  "Method", "fmcd", "Distances", mah, "Outliers", outliers);
    else
      s = struct ("NumOGKIterations", numogkiter, ...
                  "UnivariateScale", univestimator, "Mu", mu, "Sigma", sig, ...
                  "Method", "ogk", "Distances", mah, "Outliers", outliers);
    endif
  endif

endfunction

## Squared Mahalanobis distances of the rows of X to (T, S).
function d2 = mahal2 (X, T, S)
  R = chol (S);
  Z = (X - T) / R;
  d2 = sum (Z .^ 2, 2);
endfunction

## One MCD concentration step (C-step): keep the h closest points, recompute.
function [T, S, H, dt] = cstep (X, T, S, h)
  d2 = mahal2 (X, T, S);
  [~, ord] = sort (d2);
  H = ord(1:h);
  T = mean (X(H, :));
  S = cov (X(H, :));
  dt = det (S);
endfunction

## Fast-MCD raw location and scatter of the best h-subset.
function [Tbest, Sbest] = fastmcd (X, h, nsamp)
  [n, p] = size (X);
  bestdet = Inf;
  Tbest = mean (X);
  Sbest = cov (X);
  ## Candidate starts: random elemental subsets plus the classical estimate.
  cand = cell (nsamp + 1, 1);
  for t = 1:nsamp
    idx = randperm (n, p + 1);
    T = mean (X(idx, :));
    S = cov (X(idx, :));
    k = p + 1;
    while (rcond (S) < 1e-12 && k < h)
      k += 1;
      idx = randperm (n, k);
      T = mean (X(idx, :));
      S = cov (X(idx, :));
    endwhile
    cand{t} = {T, S};
  endfor
  cand{nsamp + 1} = {mean(X), cov(X)};
  ## Two C-steps from each start, then keep the best few and iterate them.
  dets = inf (nsamp + 1, 1);
  sols = cell (nsamp + 1, 1);
  for t = 1:(nsamp + 1)
    T = cand{t}{1};
    S = cand{t}{2};
    if (rcond (S) < 1e-12)
      continue;
    endif
    [T, S] = cstep (X, T, S, h);
    [T, S, H, dt] = cstep (X, T, S, h);
    dets(t) = dt;
    sols{t} = {T, S};
  endfor
  [~, ord] = sort (dets);
  nbest = min (10, sum (isfinite (dets)));
  for j = 1:nbest
    t = ord(j);
    T = sols{t}{1};
    S = sols{t}{2};
    dt_prev = Inf;
    for it = 1:50
      [T, S, H, dt] = cstep (X, T, S, h);
      if (dt >= dt_prev * (1 - 1e-12))
        break;
      endif
      dt_prev = dt;
    endfor
    if (dt < bestdet)
      bestdet = dt;
      Tbest = T;
      Sbest = S;
    endif
  endfor
endfunction

## MCD consistency factor  c = alpha / F_{chi2_{p+2}}(chi2inv(alpha, p)).
function c = mcdcons (p, alpha)
  if (alpha >= 1)
    c = 1;
  else
    c = alpha / chi2cdf (chi2inv (alpha, p), p + 2);
  endif
endfunction

## Raw MCD small-sample correction factor (Pison et al.), alpha = 0.5.
function f = mcdcnp2raw (p, n)
  if (p == 1)
    fp = 1 - exp (0.262024211897096) / n ^ 0.604756680630497;
  elseif (p == 2)
    fp = 1 - exp (0.673292623522027) / n ^ 0.691365864961895;
  else
    coeff = [-1.42764571687802, 1.26263336932151, 2; ...
             -1.06141115981725, 1.28907991440387, 3];
    fp = mcdcnp2_solve (coeff, p, n);
  endif
  f = 1 / fp;
endfunction

## Reweighted MCD small-sample correction factor (Pison et al.), alpha = 0.5.
function f = mcdcnp2rew (p, n)
  if (p == 1)
    fp = 1 - exp (1.11098143415027) / n ^ 1.5182890270453;
  elseif (p == 2)
    fp = 1 - exp (3.11101712909049) / n ^ 1.91401056721863;
  else
    coeff = [-1.02842572724793, 1.67659883081926, 2; ...
             -0.26800273450853, 1.35968562893582, 3];
    fp = mcdcnp2_solve (coeff, p, n);
  endif
  f = 1 / fp;
endfunction

## Shared p > 2 interpolation for the Pison correction coefficients.  Each row
## of COEFF is [coefficient, exponent, reference-dimension] for the two anchor
## dimensions (p = 2 and p = 3), following the robustbase covMcd formulas.
function fp = mcdcnp2_solve (coeff, p, n)
  a   = coeff(:, 1);
  b   = coeff(:, 2);
  dim = coeff(:, 3);
  y = log (- a ./ p .^ b);
  A = [1, - log(dim(1) * p ^ 2); 1, - log(dim(2) * p ^ 2)];
  c = A \ y;
  fp = 1 - exp (c(1)) / n ^ c(2);
endfunction

## Orthogonalized Gnanadesikan-Kettenring robust location and scatter.
## Maintains the affine invariant  x = A * z + centre  where z is the current
## (rotated) coordinate of an observation, so the final diagonal scatter and
## location in z-space map straight back to the original coordinates.
function [T, C] = ogk (X, niter, estimator)
  [n, p] = size (X);
  Z = X;
  A = eye (p);
  centre = zeros (p, 1);
  for iter = 1:niter
    ## Robust marginal location/scale of the current coordinates.
    m = zeros (p, 1);
    d = zeros (1, p);
    for j = 1:p
      [d(j), m(j)] = uniscale (Z(:, j), estimator);
    endfor
    d(d == 0) = 1;
    Y = (Z - m') ./ d;                 # standardized coordinates
    ## Pairwise robust covariance (Gnanadesikan-Kettenring).
    U = eye (p);
    for j = 1:p
      for k = (j + 1):p
        sp = uniscale (Y(:, j) + Y(:, k), estimator);
        sm = uniscale (Y(:, j) - Y(:, k), estimator);
        U(j, k) = 0.25 * (sp ^ 2 - sm ^ 2);
        U(k, j) = U(j, k);
      endfor
    endfor
    ## Orthogonalizing rotation E; update the affine map and rotate the data.
    [E, ~] = eig ((U + U') / 2);
    centre = A * m + centre;           # uses the pre-update A
    A = A * diag (d) * E;
    Z = Y * E;
  endfor
  ## Robust location/scale in the final rotated coordinates map straight back.
  gm = zeros (p, 1);
  gs = zeros (1, p);
  for j = 1:p
    [gs(j), gm(j)] = uniscale (Z(:, j), estimator);
  endfor
  T = (A * gm + centre)';
  C = A * diag (gs .^ 2) * A';
  C = (C + C') / 2;
endfunction

## Robust univariate scale (and location) via tau-scale or Qn.  The tau-scale
## thresholds apply to residuals over the raw MAD.  A constant factor on the
## scale cancels in the reweighted estimate, so the tau-scale carries no
## consistency factor.
function [s, m] = uniscale (x, estimator)
  x = x(:);
  med = median (x);
  s0 = median (abs (x - med));
  if (s0 == 0)
    s = 0;
    m = med;
    return;
  endif
  if (strcmp (estimator, "qn"))
    m = med;
    s = qn (x);
  else
    ## Tau-scale (Maronna & Zamar): bisquare-weighted mean, bounded rho scale.
    c1 = 4.5;
    r = (x - med) / s0;
    w = (1 - (r / c1) .^ 2) .^ 2;
    w(abs (r) > c1) = 0;
    m = sum (w .* x) / sum (w);
    c2 = 3.0;
    r2 = (x - m) / s0;
    rho = min (r2 .^ 2, c2 ^ 2);
    s = s0 * sqrt (mean (rho));
  endif
endfunction

## Qn robust scale estimator (Croux & Rousseeuw).
function s = qn (x)
  x = x(:);
  n = numel (x);
  d = abs (x - x.');
  diffs = sort (d(tril (true (n), -1)));
  h = floor (n / 2) + 1;
  k = h * (h - 1) / 2;
  if (k < 1)
    k = 1;
  endif
  s = 2.2219 * diffs(k);
endfunction

%!demo
%! ## Robust covariance is unaffected by a cluster of outliers.
%! rng (42);
%! X = [randn(80,2); 8 + randn(10,2)];
%! [sig, mu, mah, outliers] = robustcov (X);
%! gscatter (X(:,1), X(:,2), outliers);
%! title ("robustcov: inliers vs. flagged outliers");

## FMCD without bias correction is bit-exact: consistency factor c(p, |J|/n)
%!test
%! X = [0.1 0.2; -0.3 0.5; 0.4 -0.1; 0.2 0.3; -0.2 -0.4; 0.5 0.1; -0.1 0.2; ...
%!      0.3 -0.3; -0.4 0.1; 0.2 -0.2; 0.1 0.4; -0.3 -0.2; 5 5; -6 4];
%! [sig, mu, mah, ol] = robustcov (X, "BiasCorrection", 0);
%! assert_equal (mu, [0.041666666666667 0.05], 1e-12);
%! assert_equal (sig, [0.130395818314974 -0.012781705320642; ...
%!                     -0.012781705320642 0.122435282545100], 1e-9);
%! assert_equal (ol, logical ([0;0;0;0;0;0;0;0;0;0;0;0;1;1]));

## OGK is bit-exact on separated data: covariance of the retained set (1/N)
%!test
%! X = [0.1 0.2; -0.3 0.5; 0.4 -0.1; 0.2 0.3; -0.2 -0.4; 0.5 0.1; -0.1 0.2; ...
%!      0.3 -0.3; -0.4 0.1; 0.2 -0.2; 0.1 0.4; -0.3 -0.2; 5 5; -6 4];
%! [sig, mu, mah, ol] = robustcov (X, "Method", "ogk");
%! assert_equal (mu, [0.041666666666667 0.05], 1e-12);
%! assert_equal (sig, [0.080763888888889 -0.007916666666667; ...
%!                     -0.007916666666667 0.075833333333333], 1e-9);
%! assert_equal (ol, logical ([0;0;0;0;0;0;0;0;0;0;0;0;1;1]));

## FMCD with default bias correction inflates the covariance (formula factor)
%!test
%! X = [0.1 0.2; -0.3 0.5; 0.4 -0.1; 0.2 0.3; -0.2 -0.4; 0.5 0.1; -0.1 0.2; ...
%!      0.3 -0.3; -0.4 0.1; 0.2 -0.2; 0.1 0.4; -0.3 -0.2; 5 5; -6 4];
%! sig0 = robustcov (X, "BiasCorrection", 0);
%! sig1 = robustcov (X);
%! r = sig1 ./ sig0;
%! # scalar inflation factor
%! assert_equal (all (abs (r - r(1)) < 1e-12, 'all'), true);
%! assert_equal (r(1) > 1, true);

## mah is the robust Mahalanobis distance, outliers use the 0.975 cutoff
%!test
%! X = [0.1 0.2; -0.3 0.5; 0.4 -0.1; 0.2 0.3; -0.2 -0.4; 0.5 0.1; -0.1 0.2; ...
%!      0.3 -0.3; -0.4 0.1; 0.2 -0.2; 0.1 0.4; -0.3 -0.2; 5 5; -6 4];
%! [sig, mu, mah, ol] = robustcov (X);
%! d2 = sum (((X - mu) / chol (sig)) .^ 2, 2);
%! assert_equal (mah, sqrt (d2), 1e-12);
%! assert_equal (ol, mah > sqrt (chi2inv (0.975, 2)));

## The metadata structure carries the expected fields
%!test
%! X = [randn(40,2); 10 10; -10 8];
%! [sig, mu, mah, ol, s] = robustcov (X);
%! assert_equal (isfield (s, "Method") && strcmp (s.Method, "fmcd"), true);
%! assert_equal (isfield (s, "Sigma") && isfield (s, "Mu"), true);
%! assert_equal (isfield (s, "Distances") && isfield (s, "Outliers"), true);
%! [~, ~, ~, ~, s2] = robustcov (X, "Method", "ogk");
%! assert_equal (strcmp (s2.Method, "ogk"), true);
%! assert_equal (isfield (s2, "NumOGKIterations"), true);

## OGK is bit-exact in 3-D: cube vertices have identity covariance (1/N)
%!test
%! X = [1 1 1; 1 1 -1; 1 -1 1; 1 -1 -1; -1 1 1; -1 1 -1; -1 -1 1; -1 -1 -1; ...
%!      10 10 10; -10 8 -9];
%! [sig, mu, mah, ol] = robustcov (X, "Method", "ogk");
%! assert_equal (sig, eye (3), 1e-12);
%! assert_equal (mu, [0 0 0], 1e-12);
%! assert_equal (ol, logical ([0;0;0;0;0;0;0;0;1;1]));

## OGK on moderate contamination: expected values from MATLAB R2024a
%!shared A, B, C
%! A = [0.34 0.6; 0.32 0.96; -0.49 0.23; -0.84 -0.8; 2.77 -0.6; -0.21 -1.31; ...
%!      -1.57 -2.15; -0.76 1.76; -0.57 0.68; -0.88 -0.18; -1.39 -0.6; ...
%!      1.75 0.58; -0.54 0.15; -0.52 -0.7; -1.7 0.16; 1.26 -0.11; ...
%!      -1.05 0.84; -0.25 1.72; 1.09 -1.17; 2.23 2.41; 3.22 3.23; ...
%!      1.76 4.47; 1.9 2.17; 1.23 1.13];
%! B = [1.39 -0.25 -0.72; -0.31 -0.24 1.55; -0.02 -0.07 -1.59; ...
%!      1.5 -0.84 0.91; -0.04 0.02 -0.01; -1.04 -1.34 -0.22; ...
%!      -0.3 -0.25 0.62; -0.42 -2.01 -0.04; -0.81 -0.14 1.13; ...
%!      0.91 0.5 1.91; -0.06 -2.11 -2; 0.08 -1.65 -0.98; -0.85 0.25 -0.14; ...
%!      -0.12 0.12 0.66; -0.84 -1.26 -0.92; 0.12 -0.3 0.54; ...
%!      -0.3 -0.19 0.87; -0.22 0.14 0.76; 0 -1.45 1.23; 2.96 0.65 3.32; ...
%!      3.93 0.25 -0.76; 1.93 3.02 0.77; 3.35 1.08 2.56; 2.19 1.23 3];
%! C = [0.29 -0.61 -1.36; -0.17 -0.49 1.17; 0.69 0.11 -0.14; ...
%!      0.42 0.66 -1.19; 0.58 -0.41 -0.62; -0.77 -1.08 1.15; ...
%!      -1.01 -0.89 -0.12; -0.17 0.4 -0.69; -1.13 2.84 0.86; ...
%!      1.51 -0.72 -0.72; 0.64 2.95 -0.08; -0.12 -1.92 -0.55; ...
%!      1.57 -1.36 -1.85; -0.02 0.51 -0.32; -1.18 -0.25 -1.38; ...
%!      -2.82 -0.44 -0.05; 0.82 -1.14 -0.01; -0.25 0.81 -1.19; ...
%!      1.42 -0.67 -0.98; -0.14 2.88 1.61; 2.25 1.46 1.94; ...
%!      -0.18 3.07 1.89; 2.67 1.94 2.91; 2.39 -0.03 1.95];
%!test
%! [sig, mu, ~, ol] = robustcov (A, "Method", "ogk");
%! assert_equal (sig, [1.722767148760330 1.020670247933884; ...
%!                     1.020670247933884 1.646486776859504], 1e-12);
%! assert_equal (mu, [0.116818181818182 0.436363636363636], 1e-12);
%! assert_equal (find (ol)', [5, 22]);
%!test
%! [sig, mu, ~, ol] = robustcov (A, "Method", "ogk", ...
%!                               "UnivariateEstimator", "qn");
%! assert_equal (sig, [1.324408616780046 0.636798412698413; ...
%!                     0.636798412698413 1.335555555555556], 1e-12);
%! assert_equal (mu, [-0.030952380952381 0.303333333333333], 1e-12);
%! assert_equal (find (ol)', [5, 21, 22]);
%!test
%! [sig, mu, ~, ol] = robustcov (B, "Method", "ogk");
%! assert_equal (sig, [0.4791385802469 0.0781234567901 0.0867438271605; ...
%!                     0.0781234567901 0.5096617283951 0.2086302469136; ...
%!                     0.0867438271605 0.2086302469136 0.8511209876543], ...
%!               1e-12);
%! assert_equal (mu, [-0.0705555555556 -0.4977777777778 0.3088888888889], ...
%!               1e-12);
%! assert_equal (find (ol)', 20:24);
%!test
%! [sig, mu, ~, ol] = robustcov (B, "Method", "ogk", ...
%!                               "UnivariateEstimator", "qn");
%! assert_equal (sig, [0.673841 0.264091 0.378851; ...
%!                     0.264091 0.737886 0.606261; ...
%!                     0.378851 0.606261 1.394296], 1e-12);
%! assert_equal (mu, [0.043 -0.492 0.328], 1e-12);
%! assert_equal (find (ol)', 20:23);
%!test
%! [sig, mu, ~, ol] = robustcov (C, "Method", "ogk", "NumOGKIterations", 1);
%! assert_equal (sig, [0.64969 -0.302265 -0.327265; ...
%!                     -0.302265 2.27578475 0.799012; ...
%!                     -0.327265 0.799012 1.068354], 1e-12);
%! assert_equal (mu, [0.14 0.2345 -0.226], 1e-12);
%! assert_equal (find (ol)', [16, 21, 23, 24]);
%!test
%! [sig, mu, ~, ol] = robustcov (C, "Method", "ogk", "NumOGKIterations", 1, ...
%!                               "UnivariateEstimator", "qn");
%! assert_equal (sig, [1.159454320988 -0.2871074074074 -0.3650679012346; ...
%!                     -0.2871074074074 1.065891666667 0.1486175925926; ...
%!                     -0.3650679012346 0.1486175925926 0.7041015432099], ...
%!               1e-12);
%! assert_equal (mu, [-0.0188888888889 -0.2583333333333 -0.4438888888889], ...
%!               1e-12);
%! assert_equal (find (ol)', [9, 11, 20:24]);

## FMCD in 3-D returns a symmetric positive-definite estimate (p > 2 path)
%!test
%! X = [1.2 0.4 2.1; 0.8 1.1 1.9; 1.5 0.9 2.3; 0.9 0.7 1.8; 1.1 1.3 2.0; ...
%!      1.4 0.6 2.2; 0.7 1.0 1.7; 1.3 0.8 2.4; 1.0 1.2 1.6; 0.7 0.9 2.0; ...
%!      15 16 17; -9 12 -10];
%! [sig, mu, mah, ol] = robustcov (X);
%! assert_equal (size (sig), [3 3]);
%! assert_equal (sig, sig', 1e-12);
%! assert_equal (all (eig (sig) > 0, 'all'), true);
%! assert_equal (ol(11:12), logical ([1;1]));

## Rows with NaN values are dropped before estimation
%!test
%! X = [randn(30,2); NaN 1; 2 NaN];
%! [sig, mu, mah] = robustcov (X);
%! assert_equal (numel (mah), 30);

## Test input validation
%!error <robustcov: too few input arguments.> robustcov ()
%!error <robustcov: X must be a nonempty real numeric matrix.> robustcov ([])
%!error <robustcov: X must be a nonempty real numeric matrix.> robustcov ("a")
%!error <robustcov: each NAME must be followed by a VALUE.> ...
%! robustcov (ones (5,2), "Method")
%!error <robustcov: unknown parameter name 'foo'.> ...
%! robustcov (ones (5,2), "foo", "bar")
%!error <robustcov: the 'olivehawkins' method is not implemented.> ...
%! robustcov (ones (5,2), "Method", "olivehawkins")
%!error <robustcov: METHOD must be 'fmcd' or 'ogk'.> ...
%! robustcov (ones (5,2), "Method", "bogus")
%!error <robustcov: OUTLIERFRACTION must be a scalar in .0, 0.5..> ...
%! robustcov (ones (5,2), "OutlierFraction", 0.8)
%!error <robustcov: UNIVARIATEESTIMATOR must be 'tauscale' or 'qn'.> ...
%! robustcov (ones (5,2), "Method", "ogk", "UnivariateEstimator", "mad")
%!error <robustcov: X must have more rows than columns.> robustcov (ones (2,5))
