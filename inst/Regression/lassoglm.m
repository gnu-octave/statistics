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
## @deftypefn  {statistics} {@var{B} =} lassoglm (@var{X}, @var{y})
## @deftypefnx {statistics} {@var{B} =} lassoglm (@var{X}, @var{y}, @var{distr})
## @deftypefnx {statistics} {@var{B} =} lassoglm (@var{X}, @var{y}, @var{distr}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{B}, @var{FitInfo}] =} lassoglm (@dots{})
##
## Lasso and elastic-net regularized generalized linear model regression.
##
## @code{@var{B} = lassoglm (@var{X}, @var{y}, @var{distr})} returns fitted
## least-squares regression coefficients for a generalized linear model of the
## response @var{y} on the predictor data @var{X}, penalized by the lasso (L1)
## or elastic-net penalty.  @var{X} is an @math{n}-by-@math{p} numeric matrix of
## @math{p} predictors at each of @math{n} observations, and @var{y} is a
## numeric vector of @math{n} responses.  @var{distr} names the distribution of
## the response: @qcode{'normal'} (default), @qcode{'binomial'},
## @qcode{'poisson'}, @qcode{'gamma'}, or @qcode{'inverse gaussian'}.  The
## canonical link function of the chosen distribution is used unless overridden
## by the @qcode{'Link'} option.
##
## @var{B} is a @math{p}-by-@math{L} matrix, where @math{L} is the number of
## regularization (@qcode{'Lambda'}) values used; column @math{k} holds the
## coefficients for the @math{k}-th value of @math{lambda}, in order of
## ascending @math{lambda}.
##
## @code{[@var{B}, @var{FitInfo}] = lassoglm (@dots{})} also returns a structure
## @var{FitInfo} with information about the fitted models:
##
## @multitable @columnfractions 0.2 0.75
## @item @qcode{Intercept} @tab a @math{1}-by-@math{L} vector of intercept terms
## @item @qcode{Lambda} @tab the @math{1}-by-@math{L} vector of @math{lambda}
## values, in ascending order
## @item @qcode{Alpha} @tab the elastic-net mixing value used
## @item @qcode{DF} @tab the number of nonzero coefficients in each column of
## @var{B}
## @item @qcode{Deviance} @tab the deviance of the fitted model at each
## @math{lambda}; when cross-validation is requested this is instead the
## cross-validated mean deviance
## @end multitable
##
## When cross-validation is requested (see @qcode{'CV'} below), @var{FitInfo}
## additionally contains @qcode{SE} (standard error of the cross-validated
## deviance), @qcode{LambdaMinDeviance} and @qcode{IndexMinDeviance} (the
## @math{lambda} with minimum cross-validated deviance and its index), and
## @qcode{Lambda1SE} and @qcode{Index1SE} (the largest @math{lambda} within one
## standard error of that minimum).
##
## @code{lassoglm} accepts the following @var{Name}/@var{Value} pairs:
##
## @multitable @columnfractions 0.2 0.75
## @headitem Name @tab Value
## @item @qcode{'Alpha'} @tab the elastic-net mixing parameter, a scalar in
## @math{(0, 1]}.  @qcode{'Alpha'} = 1 is the lasso penalty (default);
## values towards 0 approach ridge regression.
## @item @qcode{'Lambda'} @tab a vector of non-negative @math{lambda} values.
## @item @qcode{'Standardize'} @tab a logical value (default @qcode{true})
## specifying whether the predictors are standardized before fitting.
## @item @qcode{'Weights'} @tab a vector of non-negative observation weights.
## @item @qcode{'Size'} @tab for the @qcode{'binomial'} distribution, the number
## of trials (a scalar or a per-observation vector); @var{y} holds the number of
## successes.  Default is 1 (Bernoulli responses).
## @item @qcode{'Link'} @tab the link function to use instead of the family's
## canonical link.  Accepts any link name understood by @code{glmfit} (e.g.
## @qcode{'log'}, @qcode{'probit'}) or a numeric exponent for a power link.
## @item @qcode{'Offset'} @tab a numeric vector, one value per observation,
## added as a fixed term to the linear predictor (not penalized or fitted).
## @item @qcode{'RelTol'} @tab convergence tolerance for the coordinate descent.
## @item @qcode{'MaxIter'} @tab maximum number of iterations.
## @item @qcode{'DFmax'} @tab maximum number of nonzero coefficients.
## @item @qcode{'Intercept'} @tab a logical value (default @qcode{true}) whether
## to fit an intercept term.
## @item @qcode{'PredictorNames'} @tab a cell array of predictor names.
## @item @qcode{'CV'} @tab the number of folds @math{K} for @math{K}-fold
## cross-validation, or a @code{cvpartition} object.  The fold assignment is
## random, so the selected @math{lambda} values are not reproducible without a
## fixed random seed.
## @item @qcode{'MCReps'} @tab the number of Monte-Carlo repetitions of the
## cross-validation (default 1).
## @end multitable
##
## @seealso{lasso, glmfit, glmval, cvpartition}
## @end deftypefn

function [B, FitInfo] = lassoglm (X, y, distr, varargin)

  if (nargin < 2)
    print_usage ();
  endif
  if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
    error ("lassoglm: X must be a real matrix.");
  endif
  if (! (isnumeric (y) && isreal (y) && isvector (y)))
    error ("lassoglm: Y must be a real vector.");
  endif
  y = y(:);
  if (rows (X) != numel (y))
    error ("lassoglm: X and Y must have the same number of observations.");
  endif
  if (nargin < 3 || isempty (distr))
    distr = 'normal';
  endif
  if (! (ischar (distr) && isrow (distr)))
    error ("lassoglm: DISTR must be a character vector.");
  endif
  distr = lower (distr);
  if (! any (strcmp (distr, {'normal', 'binomial', 'poisson', 'gamma', ...
                             'inverse gaussian'})))
    error ("lassoglm: unknown distribution '%s'.", distr);
  endif

  ## Defaults and Name-Value parsing.
  alpha = 1;
  lambda = [];
  numlambda = 100;
  lambdaratio = [];
  standardize = true;
  weights = [];
  reltol = 1e-4;
  maxiter = 1e5;
  dfmax = columns (X);
  intercept = true;
  prednames = {};
  binomsize = [];
  linkarg = [];
  offset = [];
  cvarg = [];
  mcreps = 1;
  usecov = false;
  if (mod (numel (varargin), 2) != 0)
    error ("lassoglm: optional arguments must be Name-Value pairs.");
  endif
  for k = 1:2:numel (varargin)
    if (! ischar (varargin{k}))
      error ("lassoglm: parameter names must be character vectors.");
    endif
    switch (lower (varargin{k}))
      case "alpha";          alpha = varargin{k+1};
      case "lambda";         lambda = varargin{k+1};
      case "numlambda";      numlambda = varargin{k+1};
      case "lambdaratio";    lambdaratio = varargin{k+1};
      case "standardize";    standardize = logical (varargin{k+1});
      case "weights";        weights = varargin{k+1};
      case "size";           binomsize = varargin{k+1};
      case "link";           linkarg = varargin{k+1};
      case "offset";         offset = varargin{k+1};
      case "reltol";         reltol = varargin{k+1};
      case "maxiter";        maxiter = varargin{k+1};
      case "dfmax";          dfmax = varargin{k+1};
      case "intercept";      intercept = logical (varargin{k+1});
      case "predictornames"; prednames = varargin{k+1};
      case "cv";             cvarg = varargin{k+1};
      case "mcreps";         mcreps = varargin{k+1};
      case "usecovariance";  usecov = logical (varargin{k+1});
      otherwise
        error ("lassoglm: unknown parameter name '%s'.", varargin{k});
    endswitch
  endfor
  if (! (isnumeric (alpha) && isscalar (alpha) && alpha > 0 && alpha <= 1))
    error ("lassoglm: 'Alpha' must be a scalar in (0, 1].");
  endif
  if (! isempty (lambda) && (! isnumeric (lambda) || any (lambda < 0)))
    error ("lassoglm: 'Lambda' must be a vector of non-negative values.");
  endif

  ## Drop observations with missing values.
  ok = ! (any (isnan (X), 2) | isnan (y));
  X = X(ok,:);  y = y(ok);
  [n, p] = size (X);

  ## Number of binomial trials per observation (Size), default one.
  if (strcmp (distr, 'binomial'))
    if (isempty (binomsize))
      N = ones (n, 1);
    elseif (isscalar (binomsize))
      N = binomsize * ones (n, 1);
    else
      N = binomsize(:);
      N = N(ok);
    endif
    if (! (isnumeric (N) && all (N > 0)))
      error ("lassoglm: 'Size' must be a positive scalar or vector.");
    endif
  else
    N = ones (n, 1);
  endif

  ## Prior observation weights.  PW keeps the raw weights for the deviance sum;
  ## O is normalised to sum to one and sets the scale of the penalty (as in
  ## LASSO for the Gaussian case).
  if (isempty (weights))
    pw = ones (n, 1);
  else
    if (! (isnumeric (weights) && isvector (weights) && numel (weights) == n
           && all (weights >= 0) && any (weights > 0)))
      error ("lassoglm: 'Weights' must be a non-negative vector, one per observation.");
    endif
    pw = weights(:);
  endif
  o = pw / sum (pw);

  ## Offset: a fixed per-observation term added to the linear predictor.
  if (isempty (offset))
    off = zeros (n, 1);
  else
    if (! (isnumeric (offset) && isvector (offset)
           && numel (offset) == numel (ok)))
      error ("lassoglm: 'Offset' must be a numeric vector, one per observation.");
    endif
    off = offset(:);
    off = off(ok);
  endif

  ## Response-domain checks for the chosen family.
  switch (distr)
    case 'binomial'
      if (any (y < 0) || any (y > N))
        error (strcat ("lassoglm: for the 'binomial' distribution Y must be", ...
                       " between 0 and 'Size'."));
      endif
    case 'poisson'
      if (any (y < 0))
        error ("lassoglm: for the 'poisson' distribution Y must be non-negative.");
      endif
    case {'gamma', 'inverse gaussian'}
      if (any (y <= 0))
        error ("lassoglm: for the '%s' distribution Y must be positive.", distr);
      endif
  endswitch

  ## Link, variance, and (unit) deviance functions for the chosen family.  A
  ## user-supplied 'Link' overrides the family's canonical link.
  [flink, dlink, ilink, varfun, devfun, mulims] = glmfamily (distr);
  if (! isempty (linkarg))
    [flink, dlink, ilink, errmsg] = getlinkfunctions (linkarg);
    if (! isempty (errmsg))
      error ("lassoglm: %s", errmsg);
    endif
  endif

  ## Standardise the predictors using the prior weights (fixed for all fits).
  mux = o' * X;
  Xc = X - mux;
  if (standardize)
    sig = sqrt (o' * (Xc .^ 2));
    sig(sig == 0) = 1;
  else
    sig = ones (1, p);
  endif
  Xs = Xc ./ sig;

  ## Intercept-only (null) fit -- offset-aware.  Supplies the starting
  ## intercept B0, the null deviance (the path-trim scale), and the null-model
  ## fitted mean MU0 used for the LAMBDA_MAX gradient.
  [b0, mu0] = glm_null_fit (off, y, N, o, flink, dlink, ilink, varfun, ...
                            mulims, intercept, maxiter, reltol);
  nulldev = sum (pw .* devfun (mu0, y ./ N, N));

  ## Regularization path.
  soft = @(z, g) sign (z) .* max (abs (z) - g, 0);
  userlambda = ! isempty (lambda);
  if (userlambda)
    lampath = sort (lambda(:)', "descend");
  else
    deta0 = dlink (mu0);
    grad = Xs' * (o .* (y ./ N - mu0) ./ (deta0 .* varfun (mu0, N)));
    lambda_max = max (abs (grad)) / alpha;
    if (isempty (lambdaratio))
      lambdaratio = ifelse (n < p, 1e-2, 1e-4);
    endif
    lampath = lambda_max * lambdaratio .^ ((0:numlambda-1) / (numlambda - 1));
  endif

  ## Fit along the path, warm-started from the previous (larger-lambda) fit.
  L = numel (lampath);
  Bstd = zeros (p, L);
  icept_eta = zeros (1, L);
  bs = zeros (p, 1);
  nkeep = L;
  for kk = 1:L
    lam = lampath(kk);
    for outer = 1:maxiter
      bouter = bs;
      ## IRLS: quadratic (working-response) approximation at the current fit.
      eta = off + b0 + Xs * bs;
      mu = ilink (eta);
      mu = max (min (mu, mulims(2)), mulims(1));
      deta = dlink (mu);
      ## Working response with the offset removed (fitted by intercept + Xs).
      zc = eta - off + (y ./ N - mu) .* deta;
      W = o ./ (deta .^ 2 .* varfun (mu, N));
      sumW = sum (W);
      ## Inner coordinate descent on the penalized weighted least squares.
      for it = 1:maxiter
        bprev = bs;
        if (intercept)
          b0 = sum (W .* (zc - Xs * bs)) / sumW;
        endif
        for j = 1:p
          rj = zc - b0 - Xs * bs + Xs(:,j) * bs(j);
          rho = sum (W .* Xs(:,j) .* rj);
          dj = sum (W .* Xs(:,j) .^ 2);
          bs(j) = soft (rho, lam * alpha) / (dj + lam * (1 - alpha));
        endfor
        if (max (abs (bs - bprev)) <= reltol * max (max (abs (bs)), 1))
          break;
        endif
      endfor
      if (max (abs (bs - bouter)) <= reltol * max (max (abs (bs)), 1))
        break;
      endif
    endfor
    ## Snap coordinate-descent residue at the KKT boundary to exact zero, so a
    ## coefficient that is numerically negligible (e.g. ~1e-16 at LAMBDA_MAX)
    ## does not count towards the degrees of freedom.
    bs(abs (bs) < 1e-9) = 0;
    Bstd(:,kk) = bs;
    icept_eta(kk) = b0;
    ## Trim the default path once the model saturates (DFmax) or the fit is
    ## essentially exact (deviance a negligible fraction of the null deviance).
    if (! userlambda)
      if (sum (bs != 0) > dfmax)
        nkeep = kk - 1;
        break;
      endif
      eta = off + b0 + Xs * bs;
      mu = max (min (ilink (eta), mulims(2)), mulims(1));
      devk = sum (pw .* devfun (mu, y ./ N, N));
      if (devk < 1e-3 * nulldev)
        nkeep = kk;
        break;
      endif
    endif
  endfor
  Bstd = Bstd(:,1:nkeep);
  icept_eta = icept_eta(1:nkeep);
  lampath = lampath(1:nkeep);

  ## Unwind to the original predictor scale; assemble ascending-lambda outputs.
  B = flip (Bstd ./ sig', 2);
  icept_eta = flip (icept_eta);
  lampath = flip (lampath);
  Lk = columns (B);
  icept = zeros (1, Lk);
  dev = zeros (1, Lk);
  for kk = 1:Lk
    icept(kk) = icept_eta(kk) - mux * B(:,kk);
    eta = off + icept(kk) + X * B(:,kk);
    mu = ilink (eta);
    mu = max (min (mu, mulims(2)), mulims(1));
    dev(kk) = sum (pw .* devfun (mu, y ./ N, N));
  endfor

  FitInfo.Intercept = icept;
  FitInfo.Lambda = lampath;
  FitInfo.Alpha = alpha;
  FitInfo.DF = sum (B != 0, 1);
  FitInfo.Deviance = dev;

  ## Cross-validation of the mean deviance along the fitted path.  The fold
  ## assignment is random, so the selected lambdas are self-consistent but not
  ## reproducible across runs (or identical to MATLAB) without a fixed seed.
  if (! isempty (cvarg))
    lam = FitInfo.Lambda;
    folddev = [];
    for rep = 1:mcreps
      if (isnumeric (cvarg))
        cvp = cvpartition (n, "KFold", cvarg);
      elseif (rep == 1)
        cvp = cvarg;
      else
        cvp = repartition (cvp);
      endif
      for f = 1:cvp.NumTestSets
        tr = training (cvp, f);
        te = test (cvp, f);
        folddev(end+1,:) = lassoglm_folddev_ ( ...
            X(tr,:), y(tr), N(tr), off(tr), pw(tr), ...
            X(te,:), y(te), N(te), off(te), lam, distr, alpha, ...
            standardize, intercept, reltol, maxiter, linkarg);
      endfor
    endfor
    cvdev = mean (folddev, 1);
    cvse = std (folddev, 0, 1) / sqrt (rows (folddev));
    FitInfo.Deviance = cvdev;
    FitInfo.SE = cvse;
    [~, imin] = min (cvdev);
    FitInfo.LambdaMinDeviance = lam(imin);
    FitInfo.IndexMinDeviance = imin;
    i1se = max (find (cvdev <= cvdev(imin) + cvse(imin)));
    FitInfo.Lambda1SE = lam(i1se);
    FitInfo.Index1SE = i1se;
  endif

  FitInfo.PredictorNames = prednames;
  FitInfo.UseCovariance = usecov;

endfunction

## Family link, variance, deviance, mean-initialiser, and mean-support limits.
## The link triplet is taken from the shared private helper GETLINKFUNCTIONS;
## the variance and deviance formulas mirror those in GLMFIT.
function [flink, dlink, ilink, varfun, devfun, mulims] = glmfamily (distr)
  switch (distr)
    case 'normal'
      [flink, dlink, ilink] = getlinkfunctions ('identity');
      varfun = @(mu, N) ones (size (mu));
      devfun = @(mu, y, N) (y - mu) .^ 2;
      mulims = [-Inf, Inf];
    case 'binomial'
      [flink, dlink, ilink] = getlinkfunctions ('logit');
      varfun = @(mu, N) mu .* (1 - mu) ./ N;
      devfun = @(mu, y, N) 2 * N .* (y .* log ((y + (y == 0)) ./ mu) + ...
               (1 - y) .* log ((1 - y + (y == 1)) ./ (1 - mu)));
      seps = sqrt (eps);
      mulims = [seps, 1 - seps];
    case 'poisson'
      [flink, dlink, ilink] = getlinkfunctions ('log');
      varfun = @(mu, N) mu;
      devfun = @(mu, y, N) 2 * (y .* log ((y + (y == 0)) ./ mu) - (y - mu));
      mulims = [realmin, Inf];
    case 'gamma'
      [flink, dlink, ilink] = getlinkfunctions ('reciprocal');
      varfun = @(mu, N) mu .^ 2;
      devfun = @(mu, y, N) 2 * (- log (y ./ mu) + (y - mu) ./ mu);
      mulims = [realmin, Inf];
    case 'inverse gaussian'
      [flink, dlink, ilink] = getlinkfunctions (-2);
      varfun = @(mu, N) mu .^ 3;
      devfun = @(mu, y, N) ((y - mu) ./ mu) .^ 2 ./ y;
      mulims = [realmin, Inf];
  endswitch
endfunction

## Intercept-only (null) fit by IRLS, honouring the offset.  Returns the fitted
## intercept B0 (linear-predictor scale) and the fitted mean vector MU0.  When
## there is no intercept the model is just the offset.
function [b0, mu0] = glm_null_fit (off, y, N, o, flink, dlink, ilink, ...
                                   varfun, mulims, intercept, maxiter, reltol)
  if (! intercept)
    b0 = 0;
    mu0 = max (min (ilink (off), mulims(2)), mulims(1));
    return;
  endif
  ## Initialise from the marginal mean of the response.
  b0 = flink (o' * (y ./ N));
  for it = 1:maxiter
    eta = off + b0;
    mu = max (min (ilink (eta), mulims(2)), mulims(1));
    deta = dlink (mu);
    zc = eta - off + (y ./ N - mu) .* deta;
    W = o ./ (deta .^ 2 .* varfun (mu, N));
    b0new = sum (W .* zc) / sum (W);
    if (abs (b0new - b0) <= reltol * max (abs (b0), 1))
      b0 = b0new;
      break;
    endif
    b0 = b0new;
  endfor
  mu0 = max (min (ilink (off + b0), mulims(2)), mulims(1));
endfunction

## Mean held-out deviance (per test observation, one value per LAMBDA) for a
## single cross-validation fold: fit on the training rows at the fixed LAMBDA
## path, then score the deviance on the test rows.
function dev = lassoglm_folddev_ (Xtr, ytr, Ntr, offtr, wtr, Xte, yte, Nte, ...
                                  offte, lam, distr, alpha, standardize, ...
                                  intercept, reltol, maxiter, linkarg)
  args = {"Lambda", lam, "Alpha", alpha, "Standardize", standardize, ...
          "Intercept", intercept, "RelTol", reltol, "MaxIter", maxiter, ...
          "Weights", wtr, "Offset", offtr};
  if (strcmp (distr, "binomial"))
    args = [args, {"Size", Ntr}];
  endif
  if (! isempty (linkarg))
    args = [args, {"Link", linkarg}];
  endif
  [Btr, Ftr] = lassoglm (Xtr, ytr, distr, args{:});
  [flink, dlink, ilink, varfun, devfun, mulims] = glmfamily (distr);
  if (! isempty (linkarg))
    [flink, dlink, ilink] = getlinkfunctions (linkarg);
  endif
  dev = zeros (1, columns (Btr));
  for k = 1:columns (Btr)
    eta = offte + Ftr.Intercept(k) + Xte * Btr(:,k);
    mu = max (min (ilink (eta), mulims(2)), mulims(1));
    dev(k) = mean (devfun (mu, yte ./ Nte, Nte));
  endfor
endfunction

function out = ifelse (cond, a, b)
  if (cond)
    out = a;
  else
    out = b;
  endif
endfunction

%!demo
%! ## Logistic regression with a lasso penalty on a small binary dataset.
%! X = [0.1, 1.2; 0.4, 0.7; 1.1, 0.2; 1.5, 1.9; 0.3, 0.5; 1.8, 1.1];
%! y = [0; 0; 1; 1; 0; 1];
%! [B, FitInfo] = lassoglm (X, y, 'binomial', 'Lambda', [0.01, 0.1])

%!demo
%! ## Poisson (count) regression with an elastic-net penalty.
%! X = [0.1, 1.2; 0.4, 0.7; 1.1, 0.2; 1.5, 1.9; 0.3, 0.5; 1.8, 1.1; 0.9, 0.3];
%! y = [1; 0; 2; 3; 1; 4; 2];
%! [B, FitInfo] = lassoglm (X, y, 'poisson', 'Lambda', [0.05, 0.2], ...
%!                          'Alpha', 0.6)

%!demo
%! ## Choosing the penalty by cross-validation.  Passing a 'CV' fold count makes
%! ## lassoglm cross-validate the deviance along the lambda path.  It then
%! ## reports the lambda with the lowest mean deviance (LambdaMinDeviance) and
%! ## the largest lambda within one standard error of it (Lambda1SE) -- the
%! ## sparser "one-standard-error" model that is often preferred.
%! ##
%! ## Here only the first two of eight predictors truly drive the response, so a
%! ## good fit should keep few nonzero coefficients.  The fold assignment is
%! ## random; the seed below just makes the printed numbers reproducible.
%! rand ('seed', 42);
%! randn ('seed', 42);
%! X = randn (60, 8);
%! beta = [1.5; -2; zeros(6, 1)];
%! y = double (rand (60, 1) < 1 ./ (1 + exp (- X * beta)));
%! [B, FitInfo] = lassoglm (X, y, 'binomial', 'CV', 5);
%! printf ('LambdaMinDeviance = %.4f (%d nonzero)\n', ...
%!         FitInfo.LambdaMinDeviance, FitInfo.DF(FitInfo.IndexMinDeviance));
%! printf ('Lambda1SE         = %.4f (%d nonzero)\n', ...
%!         FitInfo.Lambda1SE, FitInfo.DF(FitInfo.Index1SE));
%! ## Coefficients of the one-standard-error model.
%! coef_1SE = B(:, FitInfo.Index1SE)

%!shared X, yb, yp, yn
%! X = [ 0.12, -0.16,  1.35,  0.48;  0.39,  0.72, -0.91, -0.12; ...
%!       0.55, -0.11, -0.28, -0.28; -1.32,  0.80, -1.25,  2.14; ...
%!      -0.24,  0.77, -2.62,  0.66;  0.05,  0.70, -1.55,  0.06; ...
%!       1.05,  0.84, -0.16,  0.54;  0.55, -0.50, -0.25,  0.07; ...
%!       0.31, -0.58, -0.77, -0.50;  1.09,  0.44, -0.47, -0.26; ...
%!       1.08, -1.67,  2.08,  0.29; -1.52,  1.51, -2.46, -1.41; ...
%!      -0.55,  1.40, -0.59, -0.17;  1.35, -0.16, -0.38, -0.26];
%! yb = [1 1 1 0 0 1 1 1 1 1 1 0 0 1]';
%! yp = [1 0 1 1 0 1 1 2 1 1 1 0 0 0]';
%! yn = [0.7 -0.86 0.21 -1.38 -0.81 -1.16 0.12 1.33 0.87 0.65 4.18 -4.57 ...
%!       -2.76 1.51]';

## Test results (values verified against MATLAB's lassoglm)
%!test
%! ## Gaussian lassoglm coincides with lasso and matches MATLAB.
%! [Bg, Fg] = lassoglm (X, yn, "normal", "Lambda", [0.3, 0.1, 0.02]);
%! [Bl, Fl] = lasso (X, yn, "Lambda", [0.3, 0.1, 0.02]);
%! assert_equal (Bg, Bl, 1e-4);
%! assert_equal (Fg.Intercept, Fl.Intercept, 1e-4);
%!test
%! ## Poisson regression: intercept and deviance are exact to MATLAB.
%! [B, F] = lassoglm (X, yp, "poisson", "Lambda", [0.2, 0.05, 0.01]);
%! assert_equal (F.Intercept, [-0.4084524, -0.3186332, -0.2972851], 1e-6);
%! assert_equal (F.Deviance, [7.068890, 7.267237, 8.458039], 1e-5);
%! assert_equal (F.DF, [3, 3, 1]);
%!test
%! ## Binomial (logistic) regression: intercept, deviance, and DF.
%! [B, F] = lassoglm (X, yb, "binomial", "Lambda", [0.2, 0.1, 0.05, 0.01]);
%! assert_equal (F.Intercept, ...
%!               [2.282662, 1.388003, 0.9649387, 0.8340267], 1e-4);
%! assert_equal (F.Deviance, [1.147154, 3.739176, 6.057663, 9.585824], 1e-4);
%! assert_equal (F.DF, [3, 3, 3, 1]);
%!test
%! ## Default path: DF is zero at lambda_max and lambda is ascending.
%! [B, F] = lassoglm (X, yb, "binomial");
%! assert_equal (F.DF(end), 0);
%! assert_equal (issorted (F.Lambda), true);
%! assert_equal (all (F.Deviance(1:end-1) <= F.Deviance(2:end)), true);
%!test
%! ## Cross-validation adds the expected FitInfo fields and the 1-SE rule
%! ## selects a lambda no smaller than the minimum-deviance lambda.
%! cvp = cvpartition (14, "KFold", 5);
%! [B, F] = lassoglm (X, yb, "binomial", "CV", cvp);
%! assert_equal (isfield (F, "SE") && isfield (F, "LambdaMinDeviance"), true);
%! assert_equal (numel (F.SE), numel (F.Lambda));
%! assert_equal (F.Lambda1SE >= F.LambdaMinDeviance, true);
%! assert_equal (F.Index1SE >= F.IndexMinDeviance, true);

## Test input validation
%!error<Invalid call> lassoglm (1)
%!error<lassoglm: X must be a real matrix.> lassoglm ("a", [1;2])
%!error<lassoglm: Y must be a real vector.> lassoglm ([1, 2; 3, 4], "a")
%!error<lassoglm: X and Y must have the same number of observations.> ...
%! lassoglm ([1, 2; 3, 4], [1; 2; 3])
%!error<lassoglm: unknown distribution 'wibble'.> ...
%! lassoglm ([1, 2; 3, 4], [1; 0], 'wibble')
%!error<lassoglm: for the 'poisson' distribution Y must be non-negative.> ...
%! lassoglm ([1, 2; 3, 4], [1; -1], 'poisson')
%!error<lassoglm: for the 'gamma' distribution Y must be positive.> ...
%! lassoglm ([1, 2; 3, 4], [1; 0], 'gamma')
%!error<lassoglm: for the 'binomial' distribution Y must be between 0 and 'Size'.> ...
%! lassoglm ([1, 2; 3, 4], [1; 3], 'binomial')
%!error<lassoglm: 'Offset' must be a numeric vector, one per observation.> ...
%! lassoglm ([1, 2; 3, 4], [1; 0], 'binomial', 'Offset', [1, 2, 3])
%!error<lassoglm: 'Alpha' must be a scalar in .0, 1..> ...
%! lassoglm ([1, 2; 3, 4], [1; 0], 'binomial', 'Alpha', 0)
%!error<lassoglm: 'Lambda' must be a vector of non-negative values.> ...
%! lassoglm ([1, 2; 3, 4], [1; 0], 'binomial', 'Lambda', -1)
%!error<lassoglm: optional arguments must be Name-Value pairs.> ...
%! lassoglm ([1, 2; 3, 4], [1; 0], 'binomial', 'Lambda')
%!error<lassoglm: unknown parameter name 'foo'.> ...
%! lassoglm ([1, 2; 3, 4], [1; 0], 'binomial', 'foo', 1)
