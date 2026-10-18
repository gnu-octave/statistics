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
## are removed with a warning.
##
## @code{[@var{sig}, @var{mu}, @var{mah}, @var{outliers}, @var{s}] = robustcov
## (@dots{})} also returns the robust mean @var{mu} (@math{1*P}), the robust
## Mahalanobis distances @var{mah} (@math{N*1}) of each observation from the
## estimated distribution, a logical vector @var{outliers} (@math{N*1}) flagging
## observations whose distance exceeds @code{sqrt (chi2inv (0.975, @var{P}))},
## and a structure @var{s} holding the estimate metadata.  @var{mah} and
## @var{outliers} have one row per row of @var{X}; for removed rows, @var{mah}
## is @code{NaN} and @var{outliers} is @code{false}.
##
## Additional parameters can be specified by @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.22 0.76
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Method'} @tab the estimator, one of @qcode{'fmcd'} (default,
## the Fast Minimum Covariance Determinant algorithm), @qcode{'ogk'} (the
## Orthogonalized Gnanadesikan-Kettenring estimator) or
## @qcode{'olivehawkins'} (the Olive-Hawkins concentration algorithm).  If any
## of @qcode{'ReweightingMethod'}, @qcode{'StartMethod'} or
## @qcode{'NumConcentrationSteps'} is given and @qcode{'Method'} is not, the
## method is @qcode{'olivehawkins'}.
##
## @item @qcode{'OutlierFraction'} @tab the maximum fraction of outliers, a
## scalar in @math{[0, 0.5]} (default 0.5), used to set the size of the
## subsets in @qcode{'fmcd'} and @qcode{'olivehawkins'}.  It is not valid for
## @qcode{'ogk'}.
##
## @item @qcode{'NumTrials'} @tab a positive integer.  For @qcode{'fmcd'} it is
## the number of random elemental subsets drawn (default 500).  For
## @qcode{'olivehawkins'} it is the number of starts and must equal the number
## of starts in @qcode{'StartMethod'} (default 2); with @qcode{'elemental'}
## alone it is the number of random starts.  It is not valid for
## @qcode{'ogk'}.
##
## @item @qcode{'BiasCorrection'} @tab a logical scalar (default @code{true})
## that applies the small-sample bias correction to the @qcode{'fmcd'} estimate.
## Only valid for @qcode{'fmcd'}.
##
## @item @qcode{'NumOGKIterations'} @tab the number of orthogonalization
## iterations for @qcode{'ogk'}, a positive integer (default 2).  Only valid
## for @qcode{'ogk'}.
##
## @item @qcode{'UnivariateEstimator'} @tab the robust univariate
## location/scale estimator used by @qcode{'ogk'}, either @qcode{'tauscale'}
## (default) or @qcode{'qn'}.  Only valid for @qcode{'ogk'}.
##
## @item @qcode{'ReweightingMethod'} @tab the reweighting used by
## @qcode{'olivehawkins'}, either @qcode{'rfch'} (default) or @qcode{'rmvn'}.
##
## @item @qcode{'NumConcentrationSteps'} @tab the number of concentration
## steps taken from each start by @qcode{'olivehawkins'}, a positive integer
## (default 10).
##
## @item @qcode{'StartMethod'} @tab the starts of @qcode{'olivehawkins'}: one of
## @qcode{'classical'}, @qcode{'medianball'} or @qcode{'elemental'}, a function
## handle @code{[T, C] = f (X)} returning a @math{1*P} location and a
## @math{P*P} scatter, or a cell array of these.  The default is
## @code{@{'classical', 'medianball'@}}.
## @end multitable
##
## @strong{Note on reproducibility.} @qcode{'fmcd'} and the @qcode{'elemental'}
## start of @qcode{'olivehawkins'} draw random subsets, so their exact estimate
## depends on the random number generator and is not identical to
## MATLAB's on data where the optimal subset is ambiguous; on well-separated
## data both converge to the same estimate.  For @qcode{'fmcd'} with
## @qcode{'BiasCorrection'} enabled, the small-sample factor uses the published
## Pison-Van Aelst-Willems asymptotic formula, which differs from MATLAB's
## tabulated simulation values by up to about 1.6% for very small samples.
##
## @strong{Singular subsets.} If @qcode{'olivehawkins'} finds a singular
## subset, a warning is issued and the selected estimate is returned without
## reweighting.  A singular estimate gives zero robust distances.  When
## several starts end in singular subsets, the @qcode{'medianball'} start is
## selected if present, otherwise the first one listed in @qcode{'StartMethod'}.
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
  nanrows = any (isnan (X), 2);
  if (any (nanrows))
    warning ("robustcov: rows with missing data were removed.");
  endif
  X = X(! nanrows, :);
  [n, p] = size (X);

  ## Defaults
  method         = "fmcd";
  outlierfrac    = 0.5;
  numtrials      = [];
  biascorrection = true;
  numogkiter     = 2;
  univestimator  = "tauscale";
  rewmethod      = "rfch";
  numconc        = 10;
  startmethod    = [];
  given          = {};

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
      case "reweightingmethod"
        rewmethod = tolower (val);
      case "numconcentrationsteps"
        numconc = val;
      case "startmethod"
        startmethod = val;
        if (ischar (val))
          startmethod = tolower (val);
        endif
      otherwise
        error ("robustcov: unknown parameter name '%s'.", name);
    endswitch
    given{end+1} = tolower (name);
    varargin(1:2) = [];
  endwhile

  ## Validate options
  if (! any (strcmp (method, {"fmcd", "ogk", "olivehawkins"})))
    error ("robustcov: METHOD must be 'fmcd', 'ogk' or 'olivehawkins'.");
  endif
  isoh = ismember (given, {"reweightingmethod", "startmethod", ...
                           "numconcentrationsteps"});
  if (any (isoh))
    if (! any (strcmp (given, "method")))
      method = "olivehawkins";
    elseif (! strcmp (method, "olivehawkins"))
      error ("robustcov: %s is only valid for the 'olivehawkins' method.", ...
             upper (given{find (isoh, 1)}));
    endif
  endif
  if (any (strcmp (given, "biascorrection")) && ! strcmp (method, "fmcd"))
    error ("robustcov: BIASCORRECTION is only valid for the 'fmcd' method.");
  endif
  if (any (ismember (given, {"numogkiterations", "univariateestimator"})) ...
      && ! strcmp (method, "ogk"))
    error (["robustcov: NUMOGKITERATIONS and UNIVARIATEESTIMATOR are only ", ...
            "valid for the 'ogk' method."]);
  endif
  if (strcmp (method, "ogk") ...
      && any (ismember (given, {"numtrials", "outlierfraction"})))
    error (["robustcov: NUMTRIALS and OUTLIERFRACTION are not valid for the ", ...
            "'ogk' method."]);
  endif
  if (! isscalar (outlierfrac) || ! isnumeric (outlierfrac)
      || outlierfrac < 0 || outlierfrac > 0.5)
    error ("robustcov: OUTLIERFRACTION must be a scalar in [0, 0.5].");
  endif
  if (! any (strcmp (univestimator, {"tauscale", "qn"})))
    error ("robustcov: UNIVARIATEESTIMATOR must be 'tauscale' or 'qn'.");
  endif
  if (! any (strcmp (rewmethod, {"rfch", "rmvn"})))
    error ("robustcov: REWEIGHTINGMETHOD must be 'rfch' or 'rmvn'.");
  endif
  if (! isposint (numconc))
    error ("robustcov: NUMCONCENTRATIONSTEPS must be a positive integer.");
  endif
  if (! isempty (numtrials) && ! isposint (numtrials))
    error ("robustcov: NUMTRIALS must be a positive integer.");
  endif
  if (n < 2 * p)
    error ("robustcov: X must have at least twice as many rows as columns.");
  endif

  if (strcmp (method, "olivehawkins"))
    if (ischar (startmethod) && strcmp (startmethod, "elemental"))
      if (isempty (numtrials))
        numtrials = 2;
      endif
      starts = repmat ({"elemental"}, 1, numtrials);
    else
      starts = oh_start_list (startmethod);
      if (isempty (numtrials))
        numtrials = numel (starts);
      elseif (numtrials != numel (starts))
        error (["robustcov: NUMTRIALS must equal the number of starts ", ...
                "specified in STARTMETHOD."]);
      endif
    endif
  elseif (isempty (numtrials))
    numtrials = 500;
  endif

  ## Robust cutoff for reweighting and for flagging outliers.
  cutoff = chi2inv (0.975, p);

  if (strcmp (method, "fmcd"))
    ## Raw Fast-MCD subset and its robust distances.
    h = subsetsize (n, p, outlierfrac);
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
  elseif (strcmp (method, "ogk"))
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
  else
    h = subsetsize (n, p, outlierfrac);
    hs = [h, h];
    if (outlierfrac == 0)
      hs(2) = max (n - ceil (n / 5), p + 1);
    endif
    [mu, sig] = olivehawkins (X, hs, numconc, starts, rewmethod);
  endif

  ## Final robust Mahalanobis distances and outlier flags.
  md2 = mahal2 (X, mu, sig);
  mah = double (sqrt (md2));
  outliers = md2 > cutoff;

  ## Restore one row per row of the input for rows removed as missing.
  if (any (nanrows))
    dall = zeros (numel (nanrows), 1);
    oall = false (numel (nanrows), 1);
    if (ispdmat (sig))
      dall(nanrows) = NaN;
      dall(! nanrows) = mah;
      oall(! nanrows) = outliers;
    endif
    mah = dall;
    outliers = oall;
  endif

  if (nargout > 4)
    if (strcmp (method, "fmcd"))
      s = struct ("BiasCorrection", biascorrection, "NumTrials", numtrials, ...
                  "OutlierFraction", outlierfrac, "Mu", mu, "Sigma", sig, ...
                  "Method", "fmcd", "Distances", mah, "Outliers", outliers);
    elseif (strcmp (method, "ogk"))
      s = struct ("NumOGKIterations", numogkiter, ...
                  "UnivariateScale", univestimator, "Mu", mu, "Sigma", sig, ...
                  "Method", "ogk", "Distances", mah, "Outliers", outliers);
    else
      s = struct ("OutlierFraction", outlierfrac, ...
                  "reweightingMethod", rewmethod, "NumTrials", numtrials, ...
                  "Mu", mu, "Sigma", sig, "Method", "olivehawkins", ...
                  "Distances", mah, "Outliers", outliers);
    endif
  endif

endfunction

## Squared Mahalanobis distances of the rows of X to (T, S).
function d2 = mahal2 (X, T, S)
  if (! ispdmat (S))
    d2 = zeros (rows (X), 1);
    return;
  endif
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

## Subset size shared by the FMCD and Olive-Hawkins estimators.
function h = subsetsize (n, p, outlierfrac)
  h = floor (2 * floor ((n + p + 1) / 2) - n ...
             + 2 * (n - floor ((n + p + 1) / 2)) * (1 - outlierfrac));
  h = max (min (h, n), p + 1);
endfunction

function tf = isposint (x)
  tf = isnumeric (x) && isscalar (x) && isreal (x) && isfinite (x) ...
       && x >= 1 && x == fix (x);
endfunction

## True for a finite, numerically positive definite matrix.
function tf = ispdmat (S)
  tf = all (isfinite (S(:)));
  if (tf)
    [~, flag] = chol (S);
    tf = (flag == 0 && rank (S) == rows (S));
  endif
endfunction

## Normalize STARTMETHOD into a cell array of start names and function handles.
function starts = oh_start_list (startmethod)
  names = {"classical", "medianball", "elemental"};
  msg = ["robustcov: STARTMETHOD must be a name, a function handle, or a ", ...
         "cell array of them."];
  if (isnumeric (startmethod) && isempty (startmethod))
    starts = {"classical", "medianball"};
  elseif (ischar (startmethod) || is_function_handle (startmethod))
    starts = {startmethod};
  elseif (iscell (startmethod) && ! isempty (startmethod))
    starts = startmethod(:).';
  else
    error (msg);
  endif
  for i = 1:numel (starts)
    if (ischar (starts{i}))
      starts{i} = tolower (starts{i});
      if (! any (strcmp (starts{i}, names)))
        error ("robustcov: invalid start name '%s' in STARTMETHOD.", starts{i});
      endif
    elseif (! is_function_handle (starts{i}))
      error (msg);
    endif
  endfor
endfunction

## Initial location and scatter of one Olive-Hawkins start.
function [T, S] = oh_start (X, spec, h)
  [n, p] = size (X);
  if (is_function_handle (spec))
    [T, S] = spec (X);
    if (! isequal (size (T), [1, p]) || ! isequal (size (S), [p, p]))
      error (["robustcov: STARTMETHOD function must return a 1xP location ", ...
              "and a PxP scatter."]);
    endif
  elseif (strcmp (spec, "classical"))
    T = mean (X);
    S = cov (X);
  elseif (strcmp (spec, "medianball"))
    m = median (X, 1);
    d = sqrt (sum ((X - m) .^ 2, 2));
    idx = d <= median (d);
    T = mean (X(idx, :));
    S = cov (X(idx, :));
  else
    k = p + 1;
    idx = randperm (n, k);
    S = cov (X(idx, :));
    while (! ispdmat (S) && k < h)
      k += 1;
      idx = randperm (n, k);
      S = cov (X(idx, :));
    endwhile
    T = mean (X(idx, :));
  endif
endfunction

## K concentration steps from (T, S); a start with a finite singular scatter
## is kept as it is.
function [T, S, ok, sing] = oh_conc (X, T, S, h, k)
  ok = false;
  sing = false;
  if (! ispdmat (S))
    sing = true;
    ok = all (isfinite ([T(:); S(:)]));
    return;
  endif
  for i = 1:k
    [~, ord] = sort (mahal2 (X, T, S));
    idx = ord(1:h);
    Tn = mean (X(idx, :));
    Sn = cov (X(idx, :));
    if (! ispdmat (Sn))
      sing = true;
      return;
    endif
    T = Tn;
    S = Sn;
    ok = true;
  endfor
endfunction

## Olive-Hawkins estimate: attractors from each start, selection, reweighting.
function [T, S] = olivehawkins (X, hs, k, starts, rewmethod)
  nst = numel (starts);
  Ts = cell (1, nst);
  Ss = cell (1, nst);
  ok = false (1, nst);
  sing = false;
  for i = 1:nst
    hi = hs(1 + ! (ischar (starts{i}) && strcmp (starts{i}, "classical")));
    [T0, S0] = oh_start (X, starts{i}, hi);
    [Ts{i}, Ss{i}, ok(i), si] = oh_conc (X, T0, S0, hi, k);
    sing = sing || si;
  endfor
  if (! any (ok))
    error ("robustcov: unable to find a positive definite covariance estimate.");
  endif
  if (sing)
    warning ("robustcov:singularSubset", ...
             ["robustcov: a singular subset was found while running the ", ...
              "Olive-Hawkins estimate."]);
  endif
  idx = find (ok);
  if (any (cellfun (@(t) ischar (t) && strcmp (t, "medianball"), starts)))
    m = median (X, 1);
    r = median (sqrt (sum ((X - m) .^ 2, 2)));
    inball = cellfun (@(t) norm (t - m), Ts(idx)) <= r;
    if (any (inball))
      idx = idx(inball);
    endif
  endif
  dets = zeros (1, numel (idx));
  for i = 1:numel (idx)
    if (ispdmat (Ss{idx(i)}))
      dets(i) = det (Ss{idx(i)});
    endif
  endfor
  ismb = cellfun (@(t) ischar (t) && strcmp (t, "medianball"), starts(idx));
  ord = [find(ismb), find(! ismb)];
  [~, j] = min (dets(ord));
  j = idx(ord(j));
  T = Ts{j};
  S = Ss{j};
  if (! sing)
    [T, S] = oh_reweight (X, T, S, rewmethod);
  endif
endfunction

## Two reweighting steps; RFCH rescales by the median distance, RMVN by the
## truncated normal consistency factor.
function [T, S] = oh_reweight (X, T, S, rewmethod)
  [n, p] = size (X);
  cp = 0.975;
  cut = chi2inv (cp, p);
  c50 = chi2inv (0.5, p);
  S = oh_rescale (X, T, S, c50);
  for stage = 1:2
    keep = mahal2 (X, T, S) <= cut;
    if (sum (keep) <= p)
      break;
    endif
    Tn = mean (X(keep, :));
    Sn = cov (X(keep, :));
    if (! ispdmat (Sn))
      break;
    endif
    ref = c50;
    if (strcmp (rewmethod, "rmvn"))
      ref = chi2inv (0.5 * cp * n / sum (keep), p);
    endif
    T = Tn;
    S = oh_rescale (X, Tn, Sn, ref);
  endfor
endfunction

## Scale S so that the median squared distance of X to (T, S) equals REF.
function S = oh_rescale (X, T, S, ref)
  scale = median (mahal2 (X, T, S)) / ref;
  if (scale > 0 && isfinite (scale))
    S *= scale;
  endif
endfunction

%!demo
%! ## Robust covariance is unaffected by a cluster of outliers.
%! rng (42);
%! X = [randn(80,2); 8 + randn(10,2)];
%! [sig, mu, mah, outliers] = robustcov (X);
%! gscatter (X(:,1), X(:,2), outliers);
%! title ("robustcov: inliers vs. flagged outliers");

%!demo
%! ## Olive-Hawkins robust distances expose outliers that the classical
%! ## Mahalanobis distances mask.
%! rng (42);
%! X = [randn(80,2); 6 + randn(10,2)];
%! [~, ~, mah, outliers] = robustcov (X, "Method", "olivehawkins");
%! dcl = sqrt (mahal (X, X));
%! plot (dcl(! outliers), mah(! outliers), "bo", ...
%!       dcl(outliers), mah(outliers), "r+");
%! xlabel ("Classical Mahalanobis distance");
%! ylabel ("Robust distance");
%! legend ("inliers", "flagged outliers");

%!demo
%! ## Tolerance ellipses of the two Olive-Hawkins reweighting methods.
%! rng (42);
%! X = [randn(80,2); 6 + randn(10,2)];
%! [s1, m1] = robustcov (X, "ReweightingMethod", "rfch");
%! [s2, m2] = robustcov (X, "ReweightingMethod", "rmvn");
%! t = linspace (0, 2 * pi, 100)';
%! c = sqrt (chi2inv (0.975, 2));
%! e1 = m1 + c * [cos(t), sin(t)] * chol (s1);
%! e2 = m2 + c * [cos(t), sin(t)] * chol (s2);
%! plot (X(:,1), X(:,2), "k.", e1(:,1), e1(:,2), "b-", e2(:,1), e2(:,2), "r-");
%! axis equal;
%! legend ("data", "rfch", "rmvn");

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
%!shared A, B, C, D, E, F, G, H
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
%! D = [A(1:12,:); A(13:24,:) + 8];
%! E = [A(1:14,:); A(15:24,:) + 6];
%! F = A(1:23,:);
%! G = A(1:17,:);
%! H = [B, C(:,1)];
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

## Olive-Hawkins values
%!test
%! [sig, mu, ~, ol] = robustcov (A, "Method", "olivehawkins");
%! assert_equal (sig, [2.18619861641093 1.05116184773748; ...
%!                     1.05116184773748 2.20459884562983], 1e-12);
%! assert_equal (mu, [-0.0309523809523809 0.303333333333333], 1e-12);
%! assert_equal (find (ol)', 22);
%!test
%! [sig, mu, ~, ol] = robustcov (B, "Method", "olivehawkins");
%! assert_equal (sig, [0.539535973554938 0.0742685780314062 0.248554707369894; ...
%!                     0.0742685780314062 0.97500208295236 0.643292674523965; ...
%!                     0.248554707369894 0.643292674523965 1.65058254051043], ...
%!               1e-12);
%! assert_equal (mu, [-0.151111111111111 -0.601111111111111 0.237777777777778], ...
%!               1e-12);
%! assert_equal (find (ol)', 20:24);
%!test
%! [sig, mu, ~, ol] = robustcov (D, "Method", "olivehawkins");
%! assert_equal (sig, [32.4372591926917 31.9759974697833; ...
%!                     31.9759974697833 33.1062365099099], 1e-11);
%! assert_equal (mu, [3.88 4.32238095238095], 1e-12);
%! assert_equal (find (ol)', 5);
%!test
%! [sig, mu, ~, ol] = robustcov (E, "Method", "olivehawkins");
%! assert_equal (sig, [4.44755827554569 0.827837857580459; ...
%!                     0.827837857580459 3.27753304453576], 1e-12);
%! assert_equal (mu, [-0.185 -0.0985714285714285], 1e-12);
%! assert_equal (find (ol)', 15:24);
%!test
%! [sig, mu, ~, ol] = robustcov (G, "Method", "olivehawkins");
%! assert_equal (sig, [0.671430029617765 0.214180155408852; ...
%!                     0.214180155408852 1.06777070102528], 1e-12);
%! assert_equal (mu, [-0.573333333333333 -0.0313333333333333], 1e-12);
%! assert_equal (find (ol)', [5 12]);
%!test
%! [sig, mu, ~, ol] = robustcov (H, "Method", "olivehawkins");
%! assert_equal (sig, [1.4773032249017 0.633509046093723 0.903971467669543 ...
%!                     1.01676000987605; ...
%!                     0.633509046093723 1.05742783733659 0.955810531417285 ...
%!                     0.610833205323175; ...
%!                     0.903971467669543 0.955810531417285 2.01745154899905 ...
%!                     0.768192889806327; ...
%!                     1.01676000987605 0.610833205323175 0.768192889806327 ...
%!                     1.41197078963033], 1e-12);
%! assert_equal (mu, [0.2045 -0.423 0.429 0.409], 1e-12);
%! assert_equal (find (ol)', [16 20 21 22]);

## Olive-Hawkins with ReweightingMethod rmvn
%!test
%! [sig, mu, ~, ol] = robustcov (A, "Method", "olivehawkins", ...
%!                               "ReweightingMethod", "rmvn");
%! assert_equal (sig, [1.86045727877771 0.894539817250082; ...
%!                     0.894539817250082 1.8761158928325], 1e-12);
%! assert_equal (mu, [-0.0309523809523809 0.303333333333333], 1e-12);
%! assert_equal (find (ol)', [5 22]);
%!test
%! [sig, mu, ~, ol] = robustcov (B, "Method", "olivehawkins", ...
%!                               "ReweightingMethod", "rmvn");
%! assert_equal (sig, [0.244711990865239 0.107713539485174 0.194534248297532; ...
%!                     0.107713539485174 0.63579655102699 0.366819588781046; ...
%!                     0.194534248297532 0.366819588781046 0.717288266600757], ...
%!               1e-12);
%! assert_equal (mu, [-0.276 -0.52 0.464], 1e-12);
%! assert_equal (find (ol)', [1 3 4 11 20:24]);
%!test
%! [sig, mu, ~, ol] = robustcov (D, "Method", "olivehawkins", ...
%!                               "ReweightingMethod", "rmvn");
%! assert_equal (sig, [23.5704669292852 22.5233205014285; ...
%!                     22.5233205014285 21.9411394241276], 1e-11);
%! assert_equal (mu, [3.840625 3.975], 1e-12);
%! assert_equal (find (ol)', [5 8 15 17 18 19 22]);
%!test
%! [sig, mu, ~, ol] = robustcov (G, "Method", "olivehawkins", ...
%!                               "ReweightingMethod", "rmvn");
%! assert_equal (sig, [0.608331988611545 0.356119261505905; ...
%!                     0.356119261505905 1.67877003559776], 1e-12);
%! assert_equal (mu, [-0.704285714285714 -0.0257142857142857], 1e-12);
%! assert_equal (find (ol)', [5 12 16]);
%!test
%! [sig, mu, ~, ol] = robustcov (H, "Method", "olivehawkins", ...
%!                               "ReweightingMethod", "rmvn");
%! assert_equal (sig, [0.800488564971132 0.378310709774622 0.451742172174145 ...
%!                     0.643235580770011; ...
%!                     0.378310709774622 0.894583657008321 0.78836418951942 ...
%!                     0.334654019100639; ...
%!                     0.451742172174145 0.78836418951942 1.74955325492995 ...
%!                     0.504616991965037; ...
%!                     0.643235580770011 0.334654019100639 0.504616991965037 ...
%!                     0.995733910134607], 1e-12);
%! assert_equal (mu, [0.0883333333333333 -0.543888888888889 0.342222222222222 ...
%!                    0.218888888888889], 1e-12);
%! assert_equal (find (ol)', [13 16 20:23]);
%!test
%! [sig, mu] = robustcov (A(:,1), "ReweightingMethod", "rmvn");
%! assert_equal (size (sig), [1 1]);
%! assert_equal (size (mu), [1 1]);
%! assert_equal (sig > 0, true);
%!test
%! i = (1:80)';
%! X = [norminv(0.001 + 0.998 * mod (i * 0.6180339887, 1)), ...
%!      norminv(0.001 + 0.998 * mod (i * 0.7548776662, 1)), ...
%!      norminv(0.001 + 0.998 * mod (i * 0.5698402909, 1))];
%! r = robustcov (X, "Method", "olivehawkins") ./ cov (X);
%! assert_equal (r(1,1), 1.08696964749719, 1e-10);
%! r = robustcov (X, "Method", "olivehawkins", "ReweightingMethod", "rmvn") ./ cov (X);
%! assert_equal (r(1,1), 1.1180955548406, 1e-10);
%!test
%! i = (1:40)';
%! X = [norminv(0.001 + 0.998 * mod (i * 0.6180339887, 1)), ...
%!      norminv(0.001 + 0.998 * mod (i * 0.7548776662, 1))];
%! X(1:6,:) += 7;
%! [sig, ~, ~, ol] = robustcov (X, "Method", "olivehawkins");
%! r = sig ./ cov (X(! ol,:));
%! assert_equal (sum (ol), 6);
%! assert_equal (r(1,1), 1.15911540848535, 1e-10);
%! [sig, ~, ~, ol] = robustcov (X, "Method", "olivehawkins", ...
%!                              "ReweightingMethod", "rmvn");
%! r = sig ./ cov (X(! ol,:));
%! assert_equal (sum (ol), 6);
%! assert_equal (r(1,1), 0.942767404122863, 1e-10);

## NumConcentrationSteps is an exact step count
%!test
%! [sig, mu] = robustcov (B, "Method", "olivehawkins", ...
%!                        "NumConcentrationSteps", 1);
%! assert_equal (sig, [1.25376010044175 0.491174146513676 0.923105291911069; ...
%!                     0.491174146513676 0.929536611388943 0.893493059181142; ...
%!                     0.923105291911069 0.893493059181142 2.11540419879994], ...
%!               1e-12);
%! assert_equal (mu, [0.181904761904762 -0.437619047619048 0.47047619047619], ...
%!               1e-12);
%! sig2 = robustcov (B, "Method", "olivehawkins", "NumConcentrationSteps", 2);
%! assert_equal (sig2, robustcov (B, "Method", "olivehawkins"), 1e-12);
%!test
%! [sig, mu] = robustcov (D, "Method", "olivehawkins", ...
%!                        "NumConcentrationSteps", 2);
%! assert_equal (sig, [25.2931432445733 24.9466427944171; ...
%!                     24.9466427944171 26.3210109498855], 1e-11);
%! assert_equal (mu, [4.36217391304348 4.78565217391304], 1e-12);
%! sig3 = robustcov (D, "Method", "olivehawkins", "NumConcentrationSteps", 3);
%! assert_equal (sig3, robustcov (D, "Method", "olivehawkins"), 1e-12);
%!test
%! sig1 = robustcov (A, "Method", "olivehawkins", "NumConcentrationSteps", 1);
%! sig30 = robustcov (A, "Method", "olivehawkins", "NumConcentrationSteps", 30);
%! assert_equal (sig1, sig30, 1e-12);

## OutlierFraction sets the concentration subset size
%!test
%! [sig, mu] = robustcov (A, "Method", "olivehawkins", "OutlierFraction", 0);
%! assert_equal (sig, [1.91585537598172 0.594684519490792; ...
%!                     0.594684519490792 1.56440812034418], 1e-12);
%! assert_equal (mu, [0.0963636363636364 0.262272727272727], 1e-12);
%!test
%! [sig, mu] = robustcov (A, "Method", "olivehawkins", "OutlierFraction", 0.1);
%! assert_equal (sig, [2.17025568191781 0.963936502972127; ...
%!                     0.963936502972127 1.81121588386622], 1e-12);
%! assert_equal (mu, [0.232173913043478 0.391304347826087], 1e-12);
%! sig4 = robustcov (A, "Method", "olivehawkins", "OutlierFraction", 0.4);
%! assert_equal (sig4, sig, 1e-12);
%!test
%! sig = robustcov (A, "Method", "olivehawkins", "OutlierFraction", 0.05);
%! assert_equal (sig(1,1), 2.07915771918762, 1e-12);
%!test
%! sig = robustcov (C, "Method", "olivehawkins", "OutlierFraction", 0.2);
%! assert_equal (sig(1,1), 0.963904697885623, 1e-12);
%!test
%! sig = robustcov (F, "Method", "olivehawkins", "OutlierFraction", 0.3);
%! assert_equal (sig(1,1), 2.263922006869, 1e-12);
%!test
%! sig = robustcov (H, "Method", "olivehawkins", "OutlierFraction", 0);
%! assert_equal (sig(1,1), 0.90133965942982, 1e-12);
%!test
%! [sig, mu] = robustcov (A(1:10,:), "Method", "olivehawkins", ...
%!                        "OutlierFraction", 0);
%! assert_equal (sig, [0.531998699152933 0.65445152173685; ...
%!                     0.65445152173685 1.54235378327137], 1e-12);
%! assert_equal (mu, [-0.4875 -0.24625], 1e-12);
%!test
%! sig = robustcov (A(1:17,:), "Method", "olivehawkins", "OutlierFraction", 0);
%! assert_equal (sig, robustcov (G, "Method", "olivehawkins"), 1e-12);
%!test
%! [sig, mu] = robustcov (C(1:12,:), "Method", "olivehawkins", ...
%!                        "OutlierFraction", 0);
%! assert_equal (sig, [0.682305684670065 0.202111073398596 -0.396039314208948; ...
%!                     0.202111073398596 0.712912692984757 -0.230146533982739; ...
%!                     -0.396039314208948 -0.230146533982739 0.940572931634535], ...
%!               1e-12);
%! assert_equal (mu, [0.125 -0.495 -0.307], 1e-12);
%!test
%! [sig, mu] = robustcov (C(1:15,:), "Method", "olivehawkins", ...
%!                        "OutlierFraction", 0.02);
%! assert_equal (sig, [0.714900001270072 -0.158312562788954 -0.345298199782001; ...
%!                     -0.158312562788954 1.80384391946054 0.31843983688034; ...
%!                     -0.345298199782001 0.31843983688034 0.788227595223814], ...
%!               1e-12);
%! assert_equal (mu, [0.0753333333333333 -0.0173333333333333 -0.389333333333333], ...
%!               1e-12);
%!test
%! [sig, mu] = robustcov (A, "Method", "olivehawkins", "OutlierFraction", 0.02);
%! assert_equal (sig, [2.23760836399026 1.59153534794663; ...
%!                     1.59153534794663 2.86229325333207], 1e-12);
%! assert_equal (mu, [0.188260869565217 0.611739130434783], 1e-12);

%!test
%! [sig, mu] = robustcov (A, "Method", "olivehawkins", ...
%!                        "StartMethod", "classical", "OutlierFraction", 0);
%! assert_equal (sig, [2.07915771918762 1.14430744580376; ...
%!                     1.14430744580376 2.35957206487605], 1e-12);
%! assert_equal (mu, [0.295833333333333 0.56125], 1e-12);
%! [sig, mu] = robustcov (A(1:10,:), "Method", "olivehawkins", ...
%!                        "StartMethod", "classical", "OutlierFraction", 0);
%! assert_equal (sig, [0.528700494767184 0.548289625422536; ...
%!                     0.548289625422536 2.14245039772679], 1e-12);
%! assert_equal (mu, [-0.517777777777778 -0.0233333333333333], 1e-12);
%!test
%! sig1 = robustcov (A, "Method", "olivehawkins", "StartMethod", "medianball", ...
%!                   "OutlierFraction", 0, "NumConcentrationSteps", 1);
%! sig10 = robustcov (A, "Method", "olivehawkins", "StartMethod", "medianball", ...
%!                    "OutlierFraction", 0, "NumConcentrationSteps", 10);
%! assert_equal (sig1, sig10, 1e-12);
%! assert_equal (sig1, robustcov (A, "Method", "olivehawkins", ...
%!                                "OutlierFraction", 0), 1e-12);

## Starts and attractors
%!test
%! sig = robustcov (A, "Method", "olivehawkins", "StartMethod", "classical");
%! assert_equal (sig, [1.91585537598172 0.594684519490792; ...
%!                     0.594684519490792 1.56440812034418], 1e-12);
%! sigmb = robustcov (A, "Method", "olivehawkins", "StartMethod", "medianball");
%! assert_equal (sigmb, robustcov (A, "Method", "olivehawkins"), 1e-12);
%!test
%! sig = robustcov (A, "Method", "olivehawkins");
%! sig1 = robustcov (A, "Method", "olivehawkins", ...
%!                   "StartMethod", {"classical", "medianball"});
%! sig2 = robustcov (A, "Method", "olivehawkins", ...
%!                   "StartMethod", {"medianball", "classical"});
%! assert_equal (sig1, sig, 1e-12);
%! assert_equal (sig2, sig, 1e-12);
%!test
%! f = @(x) deal (mean (x), cov (x));
%! sig = robustcov (A, "Method", "olivehawkins", "StartMethod", f);
%! assert_equal (sig, robustcov (A, "Method", "olivehawkins", ...
%!                               "StartMethod", "classical"), 1e-12);
%!test
%! g = @(x) deal (median (x), eye (columns (x)));
%! s1 = robustcov (B, "Method", "olivehawkins", "StartMethod", "medianball");
%! s2 = robustcov (B, "Method", "olivehawkins", "StartMethod", g);
%! assert_equal (max (abs (s1(:) - s2(:))), 0.354744, 1e-6);
%! s1 = robustcov (H, "Method", "olivehawkins", "StartMethod", "medianball");
%! s2 = robustcov (H, "Method", "olivehawkins", "StartMethod", g);
%! assert_equal (max (abs (s1(:) - s2(:))), 0.749506, 1e-6);
%! s1 = robustcov (G, "Method", "olivehawkins", "StartMethod", "medianball");
%! s2 = robustcov (G, "Method", "olivehawkins", "StartMethod", g);
%! assert_equal (s1, s2, 1e-12);
%!test
%! rand ("state", 1);
%! sig1 = robustcov (A, "StartMethod", "elemental", "NumTrials", 5);
%! rand ("state", 1);
%! sig2 = robustcov (A, "StartMethod", "elemental", "NumTrials", 5);
%! assert_equal (sig1, sig2);
%! assert_equal (sig1, sig1', 1e-12);
%! assert_equal (all (eig (sig1) > 0), true);
%!test
%! [~, ~, ~, ~, s] = robustcov (A, "StartMethod", "elemental");
%! assert_equal (s.NumTrials, 2);
%! [sig, ~, ~, ~, s] = robustcov (A, "StartMethod", "elemental", ...
%!                                "NumTrials", 1);
%! assert_equal (s.NumTrials, 1);
%! assert_equal (size (sig), [2 2]);

## Method resolution and the metadata structure
%!test
%! [~, ~, ~, ~, s] = robustcov (A, "StartMethod", "medianball");
%! assert_equal (s.Method, "olivehawkins");
%! [~, ~, ~, ~, s] = robustcov (A, "NumTrials", 2);
%! assert_equal (s.Method, "fmcd");
%! sig = robustcov (A, "ReweightingMethod", "rmvn");
%! assert_equal (sig, robustcov (A, "Method", "olivehawkins", ...
%!                               "ReweightingMethod", "rmvn"), 1e-12);
%!test
%! [sig, mu, mah, ol, s] = robustcov (A, "Method", "olivehawkins");
%! assert_equal (s.Sigma, sig);
%! assert_equal (s.Mu, mu);
%! assert_equal (s.Distances, mah);
%! assert_equal (s.Outliers, ol);
%! assert_equal (fieldnames (s)', {"OutlierFraction", "reweightingMethod", ...
%!                                 "NumTrials", "Mu", "Sigma", "Method", ...
%!                                 "Distances", "Outliers"});
%! assert_equal (s.NumTrials, 2);
%! assert_equal (s.reweightingMethod, "rfch");
%!test
%! [~, mu, ~, ~, s] = robustcov (A, "Method", "olivehawkins", ...
%!                               "ReweightingMethod", "rmvn", ...
%!                               "NumConcentrationSteps", 5, ...
%!                               "StartMethod", "medianball", ...
%!                               "NumTrials", 1, "OutlierFraction", 0.3);
%! assert_equal (mu, [0.116818181818182 0.436363636363636], 1e-12);
%! assert_equal (s.NumTrials, 1);
%! assert_equal (s.OutlierFraction, 0.3);
%! assert_equal (s.reweightingMethod, "rmvn");

## Olive-Hawkins on small, univariate, degenerate and non-finite data
%!test
%! x1 = A(:,1);
%! [sig, mu, mah, ol] = robustcov (x1, "Method", "olivehawkins");
%! assert_equal (sig, 2.64311810032498, 1e-12);
%! assert_equal (mu, 0.295833333333333, 1e-12);
%! assert_equal (mah(1:3), [0.0271666657367651; 0.0148647793653997; ...
%!                          0.483361618674896], 1e-12);
%! assert_equal (ol, false (24, 1));
%!test
%! [sig, mu, mah, ol] = robustcov (A(1:6,:), "Method", "olivehawkins");
%! assert_equal (sig, [0.471176873257531 0.554010079012926; ...
%!                     0.554010079012926 0.774995718278673], 1e-12);
%! assert_equal (mu, [-0.1675 0.2475], 1e-12);
%! assert_equal (mah, [1.01449908743213; 0.813258155753472; 1.13104603226515; ...
%!                     1.22201620084884; 12.9622207511509; 4.28863012812016], ...
%!               1e-11);
%! assert_equal (find (ol)', [5 6]);
%!warning <singular subset> ...
%! [sig, mu, mah, ol] = robustcov (A(1:4,:), "Method", "olivehawkins");
%! assert_equal (sig, [0.34445 0.15355; 0.15355 0.06845], 1e-12);
%! assert_equal (mu, [-0.075 0.415], 1e-12);
%! assert_equal (mah, zeros (4, 1));
%! assert_equal (ol, false (4, 1));
%!warning <singular subset> ...
%! [sig, mu, mah] = robustcov (B(1:6,:), "Method", "olivehawkins");
%! assert_equal (sig, [0.34013333 0.44226667 -0.20316667; ...
%!                     0.44226667 0.57843333 -0.21868333; ...
%!                     -0.20316667 -0.21868333 0.73623333], 1e-7);
%! assert_equal (mu, [-0.36666667 -0.46333333 -0.60666667], 1e-7);
%! assert_equal (mah, zeros (6, 1));
%!test
%! [sig, mu] = robustcov (A(1:4,:), "Method", "olivehawkins", ...
%!                        "StartMethod", "classical");
%! assert_equal (sig, [0.21566688 0.14285887; 0.14285887 0.12814338], 1e-7);
%! assert_equal (mu, [0.056666667 0.59666667], 1e-7);
%!warning <singular subset> ...
%! robustcov (ones (10,2), "Method", "olivehawkins");
%!warning <singular subset> ...
%! robustcov ([A; Inf 1], "Method", "olivehawkins");
%!test
%! warning ("off", "robustcov:singularSubset", "local");
%! [sig, mu, mah, ol] = robustcov (ones (10,2), "Method", "olivehawkins");
%! assert_equal (sig, zeros (2));
%! assert_equal (mu, [1 1]);
%! assert_equal (mah, zeros (10, 1));
%! assert_equal (ol, false (10, 1));
%!test
%! warning ("off", "robustcov:singularSubset", "local");
%! [sig, mu, mah, ol] = robustcov ([A; Inf 1], "Method", "olivehawkins");
%! assert_equal (sig, [0.561683516483517 0.264930769230769; ...
%!                     0.264930769230769 0.682857692307692], 1e-12);
%! assert_equal (mu, [-0.507142857142857 0.425], 1e-12);
%! assert_equal (mah(25), Inf);
%! assert_equal (find (ol)', [5 7 12 16 19 20 21 22 23 25]);
%!warning <robustcov: a singular subset was found while running the Olive-Hawkins estimate.> ...
%! Xc = [A, A(:,1)];
%! [sig, mu, mah, ol] = robustcov (Xc, "Method", "olivehawkins");
%! assert_equal (sig, [0.19422652 0.049432576 0.19422652; ...
%!                     0.049432576 0.92357197 0.049432576; ...
%!                     0.19422652 0.049432576 0.19422652], 1e-7);
%! assert_equal (mu, [-0.45416667 0.32916667 -0.45416667], 1e-7);
%! assert_equal (mah, zeros (24, 1));
%! assert_equal (ol, false (24, 1));
%!warning <robustcov: a singular subset was found while running the Olive-Hawkins estimate.> ...
%! [sig, mu] = robustcov ([A, 2 * A(:,1)], "Method", "olivehawkins");
%! assert_equal (sig, [0.19422652 0.049432576 0.38845303; ...
%!                     0.049432576 0.92357197 0.098865152; ...
%!                     0.38845303 0.098865152 0.77690606], 1e-7);
%! assert_equal (mu, [-0.45416667 0.32916667 -0.90833333], 1e-7);
%!warning <robustcov: a singular subset was found while running the Olive-Hawkins estimate.> ...
%! [sig, mu] = robustcov ([A, A(:,1) + A(:,2)], "Method", "olivehawkins");
%! assert_equal (sig, [0.57159015 -0.25405606 0.31753409; ...
%!                     -0.25405606 0.78003333 0.52597727; ...
%!                     0.31753409 0.52597727 0.84351136], 1e-7);
%! assert_equal (mu, [-0.17083333 0.39833333 0.2275], 1e-7);
%!warning <robustcov: a singular subset was found while running the Olive-Hawkins estimate.> ...
%! Xc = [A, A(:,1)];
%! sig1 = robustcov (Xc, "Method", "olivehawkins", ...
%!                   "StartMethod", {"classical", "medianball"});
%! sig2 = robustcov (Xc, "Method", "olivehawkins", ...
%!                   "StartMethod", {"medianball", "classical"});
%! assert_equal (sig1, sig2);
%! assert_equal (sig1(1,1), 0.19422652, 1e-7);
%!warning <robustcov: a singular subset was found while running the Olive-Hawkins estimate.> ...
%! Xc = [A, A(:,1)];
%! hX = @(x) deal (median (x), ones (size (x, 2)));
%! [sig, mu] = robustcov (Xc, "Method", "olivehawkins", ...
%!                        "StartMethod", {hX, "medianball"});
%! assert_equal (sig(1,1), 0.19422652, 1e-7);
%! assert_equal (mu, [-0.45416667 0.32916667 -0.45416667], 1e-7);
%!warning <robustcov: a singular subset was found while running the Olive-Hawkins estimate.> ...
%! Xc = [A, A(:,1)];
%! hX = @(x) deal (median (x), ones (size (x, 2)));
%! [sig, mu] = robustcov (Xc, "Method", "olivehawkins", ...
%!                        "StartMethod", {"classical", hX});
%! assert_equal (sig(1,1), 2.0378775, 1e-7);
%! assert_equal (mu, [0.29583333 0.56125 0.29583333], 1e-7);
%!warning <robustcov: a singular subset was found while running the Olive-Hawkins estimate.> ...
%! Xc = [A, A(:,1)];
%! hX = @(x) deal (median (x), ones (size (x, 2)));
%! [sig, mu] = robustcov (Xc, "Method", "olivehawkins", ...
%!                        "StartMethod", {hX, "classical"});
%! assert_equal (sig, ones (3));
%! assert_equal (mu, [-0.23 0.405 -0.23], 1e-12);
%!test
%! [sig, mu, ~, ol] = robustcov (single (A), "Method", "olivehawkins");
%! assert_equal (class (sig), "single");
%! assert_equal (class (mu), "single");
%! assert_equal (class (ol), "logical");

## mah is double for single input and keeps the values of the double run
%!test
%! X = [0.1 0.2; -0.3 0.5; 0.4 -0.1; 0.2 0.3; -0.2 -0.4; 0.5 0.1; -0.1 0.2; ...
%!      0.3 -0.3; -0.4 0.1; 0.2 -0.2; 0.1 0.4; -0.3 -0.2; 5 5; -6 4];
%! for m = {"fmcd", "ogk", "olivehawkins"}
%!   [~, ~, mah0, ol0] = robustcov (X, "Method", m{1});
%!   [~, ~, mah, ol] = robustcov (single (X), "Method", m{1});
%!   assert_equal (class (mah), "double");
%!   assert_equal (mah, mah0, 1e-3);
%!   assert_equal (ol, ol0);
%! endfor

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

## Rows with NaN values are removed with a warning and padded in the outputs
%!warning <robustcov: rows with missing data were removed.> ...
%! Xn = A;
%! Xn(3,1) = NaN;
%! keep = [1:2, 4:24];
%! for m = {"ogk", "olivehawkins"}
%!   [sig, mu, mah, ol, s] = robustcov (Xn, "Method", m{1});
%!   [sig0, mu0, mah0, ol0] = robustcov (A(keep,:), "Method", m{1});
%!   assert_equal (find (isnan (mah)), 3);
%!   assert_equal (size (mah), [24 1]);
%!   assert_equal (ol(3), false);
%!   assert_equal (sig, sig0);
%!   assert_equal (mu, mu0);
%!   assert_equal (mah(keep), mah0);
%!   assert_equal (ol(keep), ol0);
%!   assert_equal (size (s.Distances), [24 1]);
%!   assert_equal (size (s.Outliers), [24 1]);
%! endfor
%!warning <robustcov: rows with missing data were removed.> ...
%! Xn = A;
%! Xn(3,1) = NaN;
%! [~, ~, mah, ol, s] = robustcov (Xn);
%! assert_equal (find (isnan (mah)), 3);
%! assert_equal (size (mah), [24 1]);
%! assert_equal (ol(3), false);
%! assert_equal (size (s.Outliers), [24 1]);
%!warning <robustcov: rows with missing data were removed.> ...
%! Xn = A(1:12,:);
%! Xn(3,:) = NaN;
%! Xn(7,2) = NaN;
%! [sig, ~, mah, ol] = robustcov (Xn, "Method", "ogk");
%! assert_equal (find (isnan (mah))', [3 7]);
%! assert_equal (sum (ol), 2);
%! assert_equal (robustcov (Xn, "Method", "ogk"), sig);
%!error <robustcov: X must have at least twice as many rows as columns.> ...
%! robustcov ([NaN(9,1), A(1:9,2); A(10:12,:)])
%!warning <robustcov: rows with missing data were removed.> ...
%! Xn = A(1:12,:);
%! Xn(1:8,1) = NaN;
%! [~, ~, mah, ol] = robustcov (Xn, "Method", "ogk");
%! assert_equal (find (isnan (mah))', 1:8);
%! assert_equal (mah(9:12), [1.41421; 1.41421; 1.41421; 20.939], 1e-4);
%! assert_equal (find (ol)', 12);
%!warning <robustcov: a singular subset was found while running the Olive-Hawkins estimate.> ...
%! Xn = [A(1:12,:), A(1:12,1)];
%! Xn(3,1) = NaN;
%! [~, ~, mah, ol] = robustcov (Xn, "Method", "olivehawkins");
%! assert_equal (mah, zeros (12, 1));
%! assert_equal (ol, false (12, 1));
%! Xn = A(1:12,:);
%! Xn(1:8,1) = NaN;
%! [~, ~, mah, ol] = robustcov (Xn, "Method", "olivehawkins");
%! assert_equal (mah, zeros (12, 1));
%! assert_equal (ol, false (12, 1));

## Test input validation
%!error <robustcov: too few input arguments.> robustcov ()
%!error <robustcov: X must be a nonempty real numeric matrix.> robustcov ([])
%!error <robustcov: X must be a nonempty real numeric matrix.> robustcov ("a")
%!error <robustcov: each NAME must be followed by a VALUE.> ...
%! robustcov (ones (5,2), "Method")
%!error <robustcov: unknown parameter name 'foo'.> ...
%! robustcov (ones (5,2), "foo", "bar")
%!error <robustcov: METHOD must be 'fmcd', 'ogk' or 'olivehawkins'.> ...
%! robustcov (ones (5,2), "Method", "bogus")
%!error <robustcov: OUTLIERFRACTION must be a scalar in .0, 0.5..> ...
%! robustcov (ones (5,2), "OutlierFraction", 0.8)
%!error <robustcov: UNIVARIATEESTIMATOR must be 'tauscale' or 'qn'.> ...
%! robustcov (ones (5,2), "Method", "ogk", "UnivariateEstimator", "mad")
%!error <robustcov: X must have at least twice as many rows as columns.> ...
%! robustcov (ones (2,5))
%!error <robustcov: REWEIGHTINGMETHOD is only valid for the 'olivehawkins' method.> ...
%! robustcov (A, "Method", "ogk", "ReweightingMethod", "rmvn")
%!error <robustcov: STARTMETHOD is only valid for the 'olivehawkins' method.> ...
%! robustcov (A, "Method", "fmcd", "StartMethod", "medianball")
%!error <robustcov: NUMCONCENTRATIONSTEPS is only valid for the 'olivehawkins' method.> ...
%! robustcov (A, "Method", "fmcd", "NumConcentrationSteps", 3)
%!error <robustcov: NUMTRIALS and OUTLIERFRACTION are not valid for the 'ogk' method.> ...
%! robustcov (A, "Method", "ogk", "OutlierFraction", 0.25)
%!error <robustcov: BIASCORRECTION is only valid for the 'fmcd' method.> ...
%! robustcov (A, "Method", "olivehawkins", "BiasCorrection", 0)
%!error <robustcov: NUMOGKITERATIONS and UNIVARIATEESTIMATOR are only valid for the 'ogk' method.> ...
%! robustcov (A, "Method", "olivehawkins", "NumOGKIterations", 1)
%!error <robustcov: REWEIGHTINGMETHOD must be 'rfch' or 'rmvn'.> ...
%! robustcov (A, "Method", "olivehawkins", "ReweightingMethod", "bogus")
%!error <robustcov: invalid start name 'bogus' in STARTMETHOD.> ...
%! robustcov (A, "Method", "olivehawkins", "StartMethod", "bogus")
%!error <robustcov: STARTMETHOD must be a name, a function handle, or a cell array of them.> ...
%! robustcov (A, "Method", "olivehawkins", "StartMethod", 5)
%!error <robustcov: STARTMETHOD function must return a 1xP location and a PxP scatter.> ...
%! robustcov (A, "StartMethod", @(x) deal (mean (x)', cov (x)))
%!error <robustcov: NUMCONCENTRATIONSTEPS must be a positive integer.> ...
%! robustcov (A, "Method", "olivehawkins", "NumConcentrationSteps", 0)
%!error <robustcov: NUMTRIALS must be a positive integer.> ...
%! robustcov (A, "Method", "olivehawkins", "NumTrials", 2.5)
%!error <robustcov: NUMTRIALS must equal the number of starts specified in STARTMETHOD.> ...
%! robustcov (A, "StartMethod", {"classical", "medianball"}, "NumTrials", 3)
%!error <robustcov: OUTLIERFRACTION must be a scalar in .0, 0.5..> ...
%! robustcov (A, "Method", "olivehawkins", "OutlierFraction", [0.1 0.2])
%!error <robustcov: STARTMETHOD must be a name, a function handle, or a cell array of them.> ...
%! robustcov (A, "StartMethod", {})
%!error <robustcov: STARTMETHOD must be a name, a function handle, or a cell array of them.> ...
%! robustcov (A, "StartMethod", {"classical", 5})
