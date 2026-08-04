## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{lambda} =} factoran (@var{X}, @var{m})
## @deftypefnx {statistics} {[@var{lambda}, @var{psi}] =} factoran (@var{X}, @var{m})
## @deftypefnx {statistics} {[@var{lambda}, @var{psi}, @var{T}] =} factoran (@var{X}, @var{m})
## @deftypefnx {statistics} {[@var{lambda}, @var{psi}, @var{T}, @var{stats}] =} factoran (@var{X}, @var{m})
## @deftypefnx {statistics} {[@var{lambda}, @var{psi}, @var{T}, @var{stats}, @var{F}] =} factoran (@var{X}, @var{m})
## @deftypefnx {statistics} {[@dots{}] =} factoran (@dots{}, @var{Name}, @var{Value})
##
## Common factor analysis.
##
## @code{@var{lambda} = factoran (@var{X}, @var{m})} fits the common factor
## model with @var{m} common factors to the @math{N}-by-@math{P} data matrix
## @var{X}, whose rows are observations and columns are variables, and returns
## the @math{P}-by-@var{m} matrix @var{lambda} of factor loadings.  The model is
##
## @example
## @var{x} = @var{mu} + @var{lambda} * @var{f} + @var{e}
## @end example
##
## @noindent
## where @var{f} are the common factors and @var{e} the variable-specific
## errors, uncorrelated with each other and with @var{f}.  The analysis is
## carried out on the correlation matrix, so the loadings are in standardized
## units and a variable's communality @code{sum (@var{lambda}(i,:) .^ 2)} and
## its specific variance @code{@var{psi}(i)} sum to one.
##
## @code{[@var{lambda}, @var{psi}] = factoran (@dots{})} also returns the
## @math{P}-by-1 vector @var{psi} of specific variances.
##
## @code{[@var{lambda}, @var{psi}, @var{T}] = factoran (@dots{})} also returns
## the @var{m}-by-@var{m} rotation matrix @var{T} that was applied to the
## loadings.  It is the identity when @qcode{'Rotate'} is @qcode{'none'}.
##
## @code{[@var{lambda}, @var{psi}, @var{T}, @var{stats}] = factoran (@dots{})}
## also returns a structure @var{stats} with the fields
##
## @multitable @columnfractions 0.15 0.8
## @headitem Field @tab Description
## @item @qcode{loglike} @tab the maximized log-likelihood, up to a constant.
## @item @qcode{dfe} @tab the error degrees of freedom,
## @code{((@var{p} - @var{m})^2 - @var{p} - @var{m}) / 2}.
## @item @qcode{chisq} @tab the likelihood ratio statistic testing @var{m}
## common factors against an unrestricted covariance.
## @item @qcode{p} @tab the significance of @var{chisq}.
## @end multitable
##
## @noindent
## The last two are present only when they can be computed: they are omitted
## when the degrees of freedom are not positive, when a specific variance has
## reached its lower bound (a Heywood case, where the likelihood is on the
## boundary), and for the @qcode{'paf'} extraction, which does not maximize a
## likelihood.  Test for them with @code{isfield} before using them.
##
## @code{[@var{lambda}, @var{psi}, @var{T}, @var{stats}, @var{F}] = factoran
## (@dots{})} also returns the @math{N}-by-@var{m} matrix @var{F} of predicted
## factor scores.  Scores are not available from a covariance matrix, only from
## data.
##
## @subheading Name-Value pairs
##
## @multitable @columnfractions 0.18 0.77
## @headitem Name @tab Value
## @item @qcode{'Extraction'} @tab How the loadings are estimated, either
## @qcode{'ml'} (default) for maximum likelihood or @qcode{'paf'} for principal
## axis factoring.  @strong{An Octave extension}; MATLAB fits by maximum
## likelihood only.  See the note below.
##
## @item @qcode{'Xtype'} @tab Whether @var{X} holds @qcode{'data'} (default) or
## a @qcode{'covariance'} (or correlation) matrix.
##
## @item @qcode{'Nobs'} @tab The number of observations behind a covariance
## matrix.  Required for @var{stats} when @qcode{'Xtype'} is
## @qcode{'covariance'}.
##
## @item @qcode{'Delta'} @tab The lower bound on the specific variances, a
## scalar in @code{[0, 1)} (default @code{0.005}).  Bounding them away from zero
## keeps the likelihood finite; a solution that reaches the bound is a Heywood
## case and is reported by a warning.
##
## @item @qcode{'Rotate'} @tab The rotation applied to the loadings, passed to
## @code{rotatefactors}: @qcode{'varimax'} (default), @qcode{'none'},
## @qcode{'quartimax'}, @qcode{'equamax'}, @qcode{'parsimax'},
## @qcode{'orthomax'}, or @qcode{'promax'}.
##
## @item @qcode{'Normalize'} @tab Whether the rotation normalizes the rows of
## the loadings (Kaiser normalization), @qcode{'on'} (default) or @qcode{'off'}.
##
## @item @qcode{'Power'} @tab The exponent of the @qcode{'promax'} target, a
## scalar not less than 1 (default 4).
##
## @item @qcode{'Scores'} @tab How @var{F} is predicted, either @qcode{'wls'}
## (default, also named @qcode{'Bartlett'}) or @qcode{'regression'} (also named
## @qcode{'Thomson'}).
##
## @item @qcode{'Maxit'} @tab The iteration limit of the extraction (default
## @code{500}).
##
## @item @qcode{'Tolerance'} @tab The convergence tolerance of the extraction
## (default @code{1e-8}).
## @end multitable
##
## @subheading Choosing the extraction
##
## The two extractions fit the same model but estimate it differently, and they
## answer to different circumstances.
##
## @qcode{'ml'} maximizes the likelihood of a multivariate normal, and is what
## MATLAB's @code{factoran} does.  Use it when you want the likelihood ratio
## test in @var{stats} to decide how many factors the data support, and when
## the data are plausibly normal.  It can fail to converge, or push a specific
## variance to zero, when the model asks for more factors than the data hold.
##
## @qcode{'paf'} iterates communalities on the reduced correlation matrix.  It
## makes no distributional assumption and is stable where maximum likelihood
## struggles, which is why it remains available here, but it provides no
## likelihood and therefore no test: @var{stats} carries only @qcode{dfe}.
##
## The two agree closely when the model fits the data well and diverge when it
## does not, so a large difference between them is itself informative.
##
## @subheading Number of factors
##
## @var{m} must leave the model identified, that is
## @code{(@var{p} - @var{m})^2 >= @var{p} + @var{m}}.  With six variables at
## most three factors can be fitted, and only the smaller counts leave degrees
## of freedom to test.
##
## @subheading References
## @enumerate
## @item
## Lawley, D. N., and Maxwell, A. E., @cite{Factor Analysis as a Statistical
## Method}, 2nd Edition, Butterworths, 1971.
##
## @item
## Joreskog, K. G., "Some contributions to maximum likelihood factor analysis",
## @cite{Psychometrika} 32(4), 443-482, 1967.
##
## @item
## Harman, H. H., @cite{Modern Factor Analysis}, 3rd Edition, University of
## Chicago Press, 1976.
## @end enumerate
##
## @seealso{rotatefactors, pca, pcacov, princomp, barttest}
## @end deftypefn

function [lambda, psi, T, stats, F] = factoran (X, m, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  ## ------------------------------------------------------------------ ##
  ## Options
  ## ------------------------------------------------------------------ ##
  extraction = "ml";
  xtype      = "data";
  nobs       = [];
  delta      = 0.005;
  rotate     = "varimax";
  normalize  = "on";
  power      = 4;
  scores     = "wls";
  maxit      = 500;
  tolerance  = 1e-8;

  if (mod (numel (varargin), 2) != 0)
    error ("factoran: Name-Value arguments must come in pairs.");
  endif
  for k = 1:2:numel (varargin)
    name = varargin{k};
    if (! ischar (name))
      error ("factoran: parameter name must be a character vector.");
    endif
    val = varargin{k+1};
    switch (lower (name))
      case "extraction"
        extraction = check_choice ("Extraction", val, {"ml", "paf"});
      case "xtype"
        xtype = check_choice ("Xtype", val, {"data", "covariance"});
      case "nobs"
        nobs = val;
      case "delta"
        delta = val;
      case "rotate"
        rotate = check_choice ("Rotate", val, {"none", "varimax", ...
                  "quartimax", "equamax", "parsimax", "orthomax", "promax"});
      case "normalize"
        normalize = check_choice ("Normalize", val, {"on", "off"});
      case "power"
        power = val;
      case "scores"
        scores = check_choice ("Scores", val, ...
                  {"wls", "bartlett", "regression", "thomson"});
      case "maxit"
        maxit = val;
      case "tolerance"
        tolerance = val;
      otherwise
        error ("factoran: unknown parameter name '%s'.", name);
    endswitch
  endfor

  if (! (isscalar (delta) && isnumeric (delta) && isreal (delta) ...
         && delta >= 0 && delta < 1))
    error ("factoran: DELTA must be a scalar in the range [0, 1).");
  endif
  if (! (isscalar (maxit) && isnumeric (maxit) && maxit >= 1))
    error ("factoran: MAXIT must be a positive integer.");
  endif
  if (! (isscalar (tolerance) && isnumeric (tolerance) && tolerance > 0))
    error ("factoran: TOLERANCE must be a positive scalar.");
  endif
  ## The two spellings of each score predictor are the same thing
  if (strcmp (scores, "bartlett"))
    scores = "wls";
  elseif (strcmp (scores, "thomson"))
    scores = "regression";
  endif

  ## ------------------------------------------------------------------ ##
  ## Data or covariance
  ## ------------------------------------------------------------------ ##
  if (! (isnumeric (X) && isreal (X) && ismatrix (X) && ndims (X) == 2))
    error ("factoran: X must be a real numeric matrix.");
  endif
  isdata = strcmp (xtype, "data");
  if (isdata)
    if (any (isnan (X(:))))
      error ("factoran: X must not contain missing values.");
    endif
    [n, p] = size (X);
    if (n < 2)
      error ("factoran: X must have at least two observations.");
    endif
    S = corr_from_data (X);
  else
    [p, pc] = size (X);
    if (p != pc)
      error ("factoran: a covariance matrix must be square.");
    endif
    if (any (any (abs (X - X') > sqrt (eps) * max (1, max (abs (X(:)))))))
      error ("factoran: a covariance matrix must be symmetric.");
    endif
    d = sqrt (diag (X));
    if (any (d <= 0))
      error ("factoran: a covariance matrix must have positive diagonal.");
    endif
    S = X ./ (d * d');
    S = (S + S') / 2;
    n = nobs;
  endif

  if (! (isscalar (m) && isnumeric (m) && isreal (m) && m >= 0 ...
         && m == fix (m)))
    error ("factoran: M must be a non-negative integer.");
  endif
  dfe = ((p - m) ^ 2 - p - m) / 2;
  if (dfe < 0)
    error (strcat ("factoran: M is too large for %d variables; the model", ...
                   " is not identified beyond %d factors."), ...
           p, max_factors (p));
  endif

  ## ------------------------------------------------------------------ ##
  ## Extraction
  ## ------------------------------------------------------------------ ##
  heywood = false;
  if (m == 0)
    ## No common factors: the model is a diagonal covariance, which on a
    ## correlation matrix is the identity, so the fit tests independence.
    L = zeros (p, 0);
    psi = ones (p, 1);
    fmin = -log (det (S));
    converged = true;
  elseif (strcmp (extraction, "ml"))
    [L, psi, fmin, converged] = ml_extract (S, m, delta, maxit, tolerance);
    if (! converged)
      warning (strcat ("factoran: maximum likelihood did not converge in", ...
                       " %d iterations."), maxit);
    endif
    heywood = any (psi <= delta * (1 + 1e-8));
  else
    [L, psi, converged] = paf_extract (S, m, maxit, tolerance);
    if (! converged)
      warning (strcat ("factoran: principal axis factoring did not", ...
                       " converge in %d iterations."), maxit);
    endif
    fmin = NaN;
    heywood = any (psi <= 0);
  endif
  if (heywood)
    warning (strcat ("factoran: some specific variances are at their lower", ...
                     " bound; the fit is a Heywood case."));
  endif

  ## ------------------------------------------------------------------ ##
  ## Statistics.  chisq and p exist only for a likelihood fit that is inside
  ## the parameter space with degrees of freedom to spare.
  ## ------------------------------------------------------------------ ##
  if (nargout > 3)
    stats = struct ("loglike", -fmin, "dfe", dfe);
    if ((strcmp (extraction, "ml") || m == 0) && dfe > 0 && ! heywood)
      if (isempty (n))
        warning (strcat ("factoran: NOBS is needed to test the fit when X", ...
                         " is a covariance matrix."));
      else
        stats.chisq = (n - 1 - (2 * p + 5) / 6 - 2 * m / 3) * fmin;
        stats.p = 1 - chi2cdf (stats.chisq, dfe);
      endif
    endif
  endif

  ## ------------------------------------------------------------------ ##
  ## Rotation
  ## ------------------------------------------------------------------ ##
  if (strcmp (rotate, "none") || m <= 1)
    lambda = L;
    T = eye (m);
  else
    rargs = {"Method", rotate, "Normalize", normalize};
    if (strcmp (rotate, "promax"))
      rargs = [rargs, {"Power", power}];
    endif
    [lambda, T] = rotatefactors (L, rargs{:});
  endif

  ## ------------------------------------------------------------------ ##
  ## Scores
  ## ------------------------------------------------------------------ ##
  if (nargout > 4)
    if (! isdata)
      error ("factoran: factor scores need data, not a covariance matrix.");
    endif
    Z = zscore_cols (X);
    if (strcmp (scores, "wls"))
      ## Bartlett: weighted least squares in the specific variances
      W = lambda ./ psi;
      F = Z * W / (lambda' * W);
    else
      ## Thomson: regression of the factors on the variables
      Sigma = lambda * lambda' + diag (psi);
      F = Z * (Sigma \ lambda);
    endif
  endif

endfunction

## Largest number of factors leaving the model identified
function mx = max_factors (p)
  mx = 0;
  for k = 1:p
    if (((p - k) ^ 2 - p - k) / 2 >= 0)
      mx = k;
    endif
  endfor
endfunction

## Validate a string-valued option against its accepted values
function out = check_choice (name, val, choices)
  if (! ischar (val))
    error ("factoran: %s must be a character vector.", name);
  endif
  idx = find (strcmpi (val, choices), 1);
  if (isempty (idx))
    error ("factoran: '%s' is not a valid value for %s.", val, name);
  endif
  out = choices{idx};
endfunction

## Correlation matrix of the columns of X
function S = corr_from_data (X)
  Z = zscore_cols (X);
  S = (Z' * Z) / (rows (X) - 1);
  S = (S + S') / 2;
  S(1:(columns (S) + 1):end) = 1;
endfunction

## Centre and scale the columns to unit variance
function Z = zscore_cols (X)
  s = std (X, 0, 1);
  if (any (s == 0))
    error ("factoran: X must not have a constant column.");
  endif
  Z = (X - mean (X, 1)) ./ s;
endfunction

## The profile discrepancy of the maximum likelihood factor model at a given
## vector of specific variances, and the loadings it implies.  Concentrating
## the likelihood on PSI this way is Joreskog's formulation: for a fixed PSI
## the loadings follow from the eigenvectors of the scaled correlation matrix,
## so only PSI has to be searched over.
function [f, L] = ml_discrepancy (psi, S, m)
  sq = sqrt (psi(:));
  Sstar = S ./ (sq * sq');
  Sstar = (Sstar + Sstar') / 2;
  [V, D] = eig (Sstar);
  [theta, idx] = sort (diag (D), "descend");
  V = V(:, idx);
  theta = max (theta, eps);
  f = sum (theta(m+1:end) - log (theta(m+1:end)) - 1);
  if (nargout > 1)
    L = (sq * ones (1, m)) .* V(:, 1:m) .* sqrt (max (theta(1:m)' - 1, 0));
  endif
endfunction

## Maximum likelihood extraction: minimize the profile discrepancy over the
## specific variances, held inside [DELTA, 1] by a logistic reparametrization
## so the search itself is unconstrained.
function [L, psi, fmin, converged] = ml_extract (S, m, delta, maxit, tol)
  p = rows (S);
  ## Joreskog's starting values
  Sinv = pinv (S);
  psi0 = (1 - 0.5 * m / p) ./ diag (Sinv);
  psi0 = min (max (psi0, delta + 1e-6), 1 - 1e-6);
  u0 = logit ((psi0 - delta) / (1 - delta));
  obj = @(u) ml_discrepancy (delta + (1 - delta) * sigmoid (u), S, m);
  opts = optimset ("MaxIter", maxit * 20, "MaxFunEvals", maxit * 200, ...
                   "TolX", tol, "TolFun", tol);
  [uhat, fmin, flag] = fminsearch (obj, u0, opts);
  converged = (flag == 1);
  psi = delta + (1 - delta) * sigmoid (uhat);
  psi = min (max (psi(:), delta), 1);

  ## Polish on the stationarity condition of the model, that the fitted
  ## covariance reproduces the unit diagonal.  The simplex search lands close
  ## to the optimum but leaves the condition satisfied only to about 1e-9,
  ## which shows up as the communalities and the specific variances not quite
  ## summing to one; iterating the condition drives it to machine precision.
  for it = 1:maxit
    [~, L] = ml_discrepancy (psi, S, m);
    psinew = min (max (1 - sum (L .^ 2, 2), delta), 1);
    if (max (abs (psinew - psi)) < eps * 8)
      psi = psinew;
      break;
    endif
    psi = psinew;
  endfor

  [fmin, L] = ml_discrepancy (psi, S, m);
  L = sign_convention (L);
endfunction

function y = sigmoid (u)
  y = 1 ./ (1 + exp (-u));
endfunction

function u = logit (y)
  y = min (max (y, 1e-12), 1 - 1e-12);
  u = log (y ./ (1 - y));
endfunction

## Principal axis factoring: iterate the communalities on the reduced
## correlation matrix until they stop moving.
function [L, psi, converged] = paf_extract (S, m, maxit, tol)
  p = rows (S);
  h2 = ones (p, 1);
  converged = false;
  L = zeros (p, m);
  for it = 1:maxit
    Rstar = S - diag (1 - h2);
    [V, D] = eig ((Rstar + Rstar') / 2);
    [ev, idx] = sort (diag (D), "descend");
    V = V(:, idx);
    L = V(:, 1:m) .* sqrt (max (ev(1:m)', 0));
    h2new = sum (L .^ 2, 2);
    if (max (abs (h2new - h2)) < tol)
      converged = true;
      h2 = h2new;
      break;
    endif
    h2 = h2new;
  endfor
  psi = max (1 - h2, 0);
  L = sign_convention (L);
endfunction

## Make each column's largest-magnitude loading positive, so that a solution
## does not change sign between runs for no reason.
function L = sign_convention (L)
  [~, idx] = max (abs (L), [], 1);
  for j = 1:columns (L)
    if (L(idx(j), j) < 0)
      L(:, j) = -L(:, j);
    endif
  endfor
endfunction


%!demo
%! ## Six measured variables built from two underlying factors, plus noise.
%! ## Factor analysis recovers the structure without being told it: the first
%! ## three variables load on one factor and the last three on the other.
%! randn ("seed", 42);
%! F = randn (300, 2);
%! X = F * [0.8 0.1; 0.7 0.2; 0.75 0.15; 0.15 0.8; 0.2 0.7; 0.1 0.75]' ...
%!     + 0.6 * randn (300, 6);
%! lambda = factoran (X, 2);
%! printf ("loadings on the two rotated factors:\n");
%! disp (round (lambda * 1000) / 1000);

%!demo
%! ## How many factors do the data support?  The likelihood ratio test in
%! ## stats answers it.  These data were built from two factors, and the test
%! ## rejects one factor while accepting two.  Three factors leave no degrees
%! ## of freedom, so there is nothing left to test with.
%! randn ("seed", 42);
%! F = randn (300, 2);
%! X = F * [0.8 0.1; 0.7 0.2; 0.75 0.15; 0.15 0.8; 0.2 0.7; 0.1 0.75]' ...
%!     + 0.6 * randn (300, 6);
%! for m = 1:3
%!   [~, ~, ~, stats] = factoran (X, m);
%!   if (isfield (stats, "p"))
%!     printf ("%d factor(s): chisq = %8.3f, dfe = %d, p = %.4f\n", ...
%!             m, stats.chisq, stats.dfe, stats.p);
%!   else
%!     printf ("%d factor(s): nothing to test against (dfe = %d)\n", ...
%!             m, stats.dfe);
%!   endif
%! endfor

%!demo
%! ## Rotation decides how a fit is presented, not how good it is.  The
%! ## unrotated solution puts most of the variance on a general first factor;
%! ## varimax turns it so each variable loads mainly on one factor, which is
%! ## easier to read.  The specific variances are untouched either way.
%! randn ("seed", 42);
%! F = randn (300, 2);
%! X = F * [0.8 0.1; 0.7 0.2; 0.75 0.15; 0.15 0.8; 0.2 0.7; 0.1 0.75]' ...
%!     + 0.6 * randn (300, 6);
%! [Lnone, psi_none] = factoran (X, 2, "Rotate", "none");
%! [Lvari, psi_vari, T] = factoran (X, 2, "Rotate", "varimax");
%! printf ("unrotated:\n");  disp (round (Lnone * 100) / 100);
%! printf ("varimax:\n");    disp (round (Lvari * 100) / 100);
%! printf ("the rotation matrix takes one to the other: %d\n", ...
%!         max (max (abs (Lnone * T - Lvari))) < 1e-10);
%! printf ("specific variances unchanged: %d\n", ...
%!         max (abs (psi_none - psi_vari)) < 1e-10);

%!demo
%! ## Factor scores place each observation on the factors, so they can be
%! ## plotted or used as inputs downstream.  The two predictors optimise
%! ## different things and are not equal, but they agree closely on the
%! ## ordering of observations.
%! randn ("seed", 42);
%! F = randn (300, 2);
%! X = F * [0.8 0.1; 0.7 0.2; 0.75 0.15; 0.15 0.8; 0.2 0.7; 0.1 0.75]' ...
%!     + 0.6 * randn (300, 6);
%! [~, ~, ~, ~, Fwls] = factoran (X, 2, "Scores", "wls");
%! [~, ~, ~, ~, Freg] = factoran (X, 2, "Scores", "regression");
%! printf ("first three observations, weighted least squares:\n");
%! disp (round (Fwls(1:3,:) * 1000) / 1000);
%! printf ("first three observations, regression:\n");
%! disp (round (Freg(1:3,:) * 1000) / 1000);
%! printf ("the two agree on factor 1 to a correlation of %.4f\n", ...
%!         corr (Fwls(:,1), Freg(:,1)));

%!demo
%! ## Two ways to estimate the same model.  Maximum likelihood is the default
%! ## and is what MATLAB does; principal axis factoring is an Octave
%! ## extension that assumes no distribution.  They agree closely when the
%! ## model fits, so a large gap between them is a warning about the fit.
%! ## Only maximum likelihood carries a test.
%! randn ("seed", 42);
%! F = randn (300, 2);
%! X = F * [0.8 0.1; 0.7 0.2; 0.75 0.15; 0.15 0.8; 0.2 0.7; 0.1 0.75]' ...
%!     + 0.6 * randn (300, 6);
%! Lml = factoran (X, 2, "Extraction", "ml");
%! Lpaf = factoran (X, 2, "Extraction", "paf");
%! printf ("largest loading difference between the extractions: %.4f\n", ...
%!         max (abs (abs (Lml(:)) - abs (Lpaf(:)))));
%! [~, ~, ~, sml] = factoran (X, 2, "Extraction", "ml");
%! [~, ~, ~, spaf] = factoran (X, 2, "Extraction", "paf");
%! printf ("ml  reports: %s\n", strjoin (fieldnames (sml)', ", "));
%! printf ("paf reports: %s\n", strjoin (fieldnames (spaf)', ", "));


## Reference values below are MATLAB R2024a's, on data generated by a
## deterministic stream so that both engines see identical bytes.
%!shared X, Lref, Pref, Tref
%! s = 7;  u = zeros (3200, 1);
%! for k = 1:3200
%!   s = mod (16807 * s, 2147483647);
%!   u(k) = s / 2147483647;
%! endfor
%! u = min (max (u, 1e-12), 1 - 1e-12);
%! a = u(1:2:end);  b = u(2:2:end);  r = sqrt (-2 * log (a));
%! z = zeros (3200, 1);
%! z(1:2:end) = r .* cos (2 * pi * b);
%! z(2:2:end) = r .* sin (2 * pi * b);
%! Z = reshape (z(1:1600), 200, 8);
%! X = Z(:,1:2) * [0.7 0.1; 0.6 0.2; 0.65 0.15; ...
%!                 0.15 0.7; 0.2 0.65; 0.1 0.6]' + 0.75 * Z(:,3:8);
%! Lref = [0.23744646899123, 0.62572191026907; ...
%!         0.1948126299518,  0.57314527421474; ...
%!         0.038102881086025, 0.63721347868091; ...
%!         0.76535682403077, 0.15782050760078; ...
%!         0.54869037618166, 0.1041708573065; ...
%!         0.52204369118599, 0.15800651296731];
%! Pref = [0.55209126537284; 0.63355253385654; 0.59250715304029; ...
%!         0.38932161929025; 0.68808730357355; 0.70250432635276];
%! Tref = [0.78348287808137, 0.62141337268628; ...
%!         -0.62141337268628, 0.78348287808137];

%!test  # the maximum likelihood fit, its rotation and its statistics
%! [L, psi, T, stats] = factoran (X, 2);
%! assert_equal (L, Lref, 1e-6);
%! assert_equal (psi, Pref, 1e-6);
%! assert_equal (T, Tref, 1e-6);
%! assert_equal (stats.loglike, -0.0043677316938693, 1e-9);
%! assert_equal (stats.dfe, 4);
%! assert_equal (stats.chisq, 0.8509797250222, 1e-6);
%! assert_equal (stats.p, 0.93148578937045, 1e-8);

%!test  # a single factor, where the model does not fit and the test says so
%! [L, psi, T, stats] = factoran (X, 1);
%! assert_equal (L, [0.58734985138995; 0.52528039923063; 0.42927173838555; ...
%!                   0.5943390974845; 0.47167200636347; ...
%!                   0.49418848792238], 1e-6);
%! assert_equal (stats.dfe, 9);
%! assert_equal (stats.chisq, 51.130175056732, 1e-5);
%! assert_equal (stats.p < 1e-6, true);
%! assert_equal (T, 1);

%!test  # rotation changes the loadings but not the fit
%! [Ln, psin, Tn] = factoran (X, 2, "Rotate", "none");
%! assert_equal (Ln, [0.57486720553951, 0.34268999200789; ...
%!                    0.50879249789022, 0.32799033558027; ...
%!                    0.42582593184473, 0.47556821038442; ...
%!                    0.69771574115811, -0.3519532998141; ...
%!                    0.49462267888081, -0.25934745412885; ...
%!                    0.50719965378403, -0.20060953329425], 1e-6);
%! assert_equal (Tn, eye (2));
%! assert_equal (psin, Pref, 1e-6);
%! [~, psiv] = factoran (X, 2, "Rotate", "varimax");
%! assert_equal (psiv, psin, 1e-10);

%!test  # the rotation matrix is the one that was applied
%! [Ln, ~, ~] = factoran (X, 2, "Rotate", "none");
%! [Lv, ~, T] = factoran (X, 2, "Rotate", "varimax");
%! assert_equal (Ln * T, Lv, 1e-8);

%!test  # quartimax, another orthogonal rotation
%! L = factoran (X, 2, "Rotate", "quartimax");
%! assert_equal (L, [0.24417499900219, 0.62312703719984; ...
%!                   0.20097710825146, 0.57101284407823; ...
%!                   0.044966808556795, 0.63676591702752; ...
%!                   0.76701294807206, 0.14956442825637; ...
%!                   0.54978098993064, 0.098252529419005; ...
%!                   0.52371594497243, 0.15237218456401], 1e-5);

%!test  # communality and specific variance partition the unit variance
%! [L, psi] = factoran (X, 2);
%! assert_equal (sum (L .^ 2, 2) + psi, ones (6, 1), 1e-10);

%!test  # both score predictors, and the default is the weighted least squares
%! [~, ~, ~, ~, Fw] = factoran (X, 2, "Scores", "wls");
%! assert_equal (Fw(1:2,:), [0.8685910477393, 3.3629807100209; ...
%!                           -0.03713794537533, -1.9311992150794], 1e-5);
%! [~, ~, ~, ~, Fr] = factoran (X, 2, "Scores", "regression");
%! assert_equal (Fr(1:2,:), [0.94397703443356, 2.2276228684448; ...
%!                           -0.22623256144569, -1.2312105486235], 1e-5);
%! [~, ~, ~, ~, Fd] = factoran (X, 2);
%! assert_equal (Fd, Fw, 1e-12);
%! [~, ~, ~, ~, Fb] = factoran (X, 2, "Scores", "Bartlett");
%! assert_equal (Fb, Fw, 1e-12);
%! [~, ~, ~, ~, Ft] = factoran (X, 2, "Scores", "Thomson");
%! assert_equal (Ft, Fr, 1e-12);

%!test  # a correlation matrix gives the same fit as the data behind it
%! [L, psi, T, stats] = factoran (corr (X), 2, "Xtype", "covariance", ...
%!                                "Nobs", 200);
%! assert_equal (L, Lref, 1e-6);
%! assert_equal (psi, Pref, 1e-6);
%! assert_equal (stats.chisq, 0.85097972502237, 1e-6);

%!test  # principal axis factoring is available and agrees where the fit is good
%! Lml = factoran (X, 2, "Extraction", "ml");
%! Lpaf = factoran (X, 2, "Extraction", "paf");
%! assert_equal (size (Lpaf), [6, 2]);
%! assert_equal (max (abs (abs (Lml(:)) - abs (Lpaf(:)))) < 0.15, true);

%!test  # but it reports no likelihood, so no test comes with it
%! [~, ~, ~, s] = factoran (X, 2, "Extraction", "paf");
%! assert_equal (isnan (s.loglike), true);
%! assert_equal (s.dfe, 4);
%! assert_equal (isfield (s, "chisq"), false);
%! assert_equal (isfield (s, "p"), false);

%!test  # the degrees of freedom follow the variable and factor counts
%! for m = 1:3
%!   [~, ~, ~, s] = factoran (X, m);
%!   assert_equal (s.dfe, ((6 - m) ^ 2 - 6 - m) / 2);
%! endfor

%!test  # with no degrees of freedom there is nothing to test
%! [~, ~, ~, s] = factoran (X, 3);
%! assert_equal (s.dfe, 0);
%! assert_equal (isfield (s, "chisq"), false);

%!test  # no common factors is a model too, and it tests independence
%! [L, psi, T, stats] = factoran (X, 0);
%! assert_equal (size (L), [6, 0]);
%! assert_equal (psi, ones (6, 1));
%! assert_equal (stats.dfe, 15);
%! assert_equal (stats.loglike, log (det (corr (X))), 1e-10);
%! assert_equal (stats.p < 1e-6, true);

%!test  # the loadings do not change sign between runs
%! L1 = factoran (X, 2);
%! L2 = factoran (X, 2);
%! assert_equal (L1, L2);

## Test input validation
%!error <Invalid call to factoran> factoran (1)
%!error <factoran: M must be a non-negative integer.> ...
%! factoran (rand (20, 5), 1.5)
%!error <factoran: M must be a non-negative integer.> ...
%! factoran (rand (20, 5), -1)
%!error <factoran: M is too large for 6 variables> factoran (rand (30, 6), 4)
%!error <factoran: X must be a real numeric matrix.> factoran ("abc", 1)
%!error <factoran: X must not contain missing values.> ...
%! factoran ([rand(20, 4); NaN(1, 4)], 1)
%!error <factoran: X must not have a constant column.> ...
%! factoran ([rand(20, 3), ones(20, 1)], 1)
%!error <factoran: X must have at least two observations.> ...
%! factoran (rand (1, 5), 1)
%!error <factoran: Name-Value arguments must come in pairs.> ...
%! factoran (rand (20, 5), 1, "Rotate")
%!error <factoran: unknown parameter name 'bogus'.> ...
%! factoran (rand (20, 5), 1, "bogus", 1)
%!error <'nosuch' is not a valid value for Extraction.> ...
%! factoran (rand (20, 5), 1, "Extraction", "nosuch")
%!error <'nosuch' is not a valid value for Rotate.> ...
%! factoran (rand (20, 5), 1, "Rotate", "nosuch")
%!error <'nosuch' is not a valid value for Scores.> ...
%! factoran (rand (20, 5), 1, "Scores", "nosuch")
%!error <factoran: DELTA must be a scalar in the range \[0, 1\).> ...
%! factoran (rand (20, 5), 1, "Delta", 1)
%!error <factoran: a covariance matrix must be square.> ...
%! factoran (rand (5, 4), 1, "Xtype", "covariance")
%!error <factoran: factor scores need data, not a covariance matrix.> ...
%! [a, b, c, d, e] = factoran (corr (rand (20, 5)), 1, "Xtype", "covariance");
