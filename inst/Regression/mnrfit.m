## Copyright (C) 2024 Andrew C Penn
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{B} =} mnrfit (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{B} =} mnrfit (@var{X}, @var{Y}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{B}, @var{dev}] =} mnrfit (@dots{})
## @deftypefnx {statistics} {[@var{B}, @var{dev}, @var{stats}] =} mnrfit (@dots{})
##
## Fit a multinomial logistic regression model.
##
## Nominal models are fitted with a baseline-category multinomial logit, using
## the last category of @var{Y} as the reference.  Ordinal models are fitted
## with a cumulative link model and hierarchical models with a sequential
## (continuation-ratio) link model, both honouring the @qcode{'link'} option
## below.  Nominal models always use the logit link.
##
## @code{@var{B} = mnrfit (@var{X}, @var{Y})}  returns a matrix, @var{B}, of
## coefficient estimates for a multinomial logistic regression of the nominal
## responses in @var{Y} on the predictors in @var{X}.  @var{X} is an @math{N*P}
## numeric matrix the observations on predictor variables, where @math{N}
## corresponds to the number of observations and @math{P} corresponds to
## predictor variables.  @var{Y} contains the response category labels and it
## can be an @math{N*K} numeric or logical matrix of 1s and 0s, with a column
## for each of the @math{K} categories, or an @math{N*1} vector of positive
## integers, logical values, character vectors in a cell array, strings, or
## categorical values.  A cell array, string, or categorical @var{Y} can also be
## a row vector.  @var{Y} can also be a character matrix with each row holding
## the label of an observation of @var{X}.
##
## The categories of a numeric or logical @var{Y} are taken in sorted order,
## those of a cell array, string, or character matrix in order of first
## appearance, and those of a categorical @var{Y} in the order of its
## categories.  A category without observations is ignored wherever it lies
## among the categories of a categorical @var{Y}; MATLAB ignores only a last one
## and fits coefficients for any other.  The last category is the reference
## category of a nominal model.  An observation whose label is missing (an empty
## character vector, @qcode{<missing>}, or @qcode{<undefined>}), or whose
## predictors hold a @qcode{NaN}, is left out of the fit, and its rows of the
## residuals in @var{stats} are @qcode{NaN}.
##
## A @var{Y} with a single category leaves nothing to fit: a nominal or
## hierarchical model returns a @var{B} with no columns, a deviance of 0 and
## residuals of 0, and an ordinal model is an error.
##
## @code{@var{B} = mnrfit (@var{X}, @var{Y}, @var{name}, @var{value})} returns a
## matrix, @var{B}, of coefficient estimates for a multinomial model fit with
## additional parameters specified @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'model'} @tab The type of model to fit: @qcode{'nominal'}
## (default) for a baseline-category model, @qcode{'ordinal'} for a cumulative
## model, or @qcode{'hierarchical'} for a sequential (continuation-ratio) model.
##
## @item @qcode{'link'} @tab The link function for ordinal and hierarchical
## models: @qcode{'logit'} (default), @qcode{'probit'}, @qcode{'comploglog'}, or
## @qcode{'loglog'}.  Nominal models always use the logit link.
##
## @item @qcode{'estdisp'} @tab @qcode{'on'} to estimate a dispersion
## parameter, scaling the coefficient standard errors by it and testing the
## coefficients against the @math{t} distribution, or @qcode{'off'} (default)
## for the theoretical dispersion of @math{1}.
##
## @item @qcode{'display'} @tab A flag to enable/disable displaying
## information about the fitted model.  Default is @qcode{'off'}.
## @end multitable
##
## @code{[@var{B}, @var{dev}, @var{stats}] = mnrfit (@dots{})} also returns the
## deviance of the fit, @var{dev}, and a structure @var{stats} with the fitted
## coefficients @qcode{'beta'} (same as @var{B}), their standard errors
## @qcode{'se'}, covariance matrix @qcode{'covb'}, correlation matrix
## @qcode{'coeffcorr'}, error degrees of freedom @qcode{'dfe'}, the coefficient
## @math{t} statistics @qcode{'t'} and @math{p}-values @qcode{'p'}, the
## dispersion parameters @qcode{'s'}, @qcode{'sfit'}, and @qcode{'estdisp'}, and
## the raw, Pearson, and deviance residuals @qcode{'resid'}, @qcode{'residp'},
## and @qcode{'residd'}.
##
## @seealso{mnrval, logistic_regression}
## @end deftypefn

function [B, DEV, STATS] = mnrfit (X, Y, varargin)

  ## Check input arguments X and Y
  if (nargin < 2)
    error ("mnrfit: too few input arguments.");
  endif
  if (! isnumeric (X))
    error ("mnrfit: Predictors must be numeric.")
  endif
  if (isscalar (X) || (ndims (X) > 2))
    error ("mnrfit: Predictors must be a vector or a 2D matrix.")
  endif
  if (isscalar (Y) || (ndims (Y) > 2))
    error ("mnrfit: Response must be a vector or a 2D matrix.")
  endif
  ## A row vector of text or categories holds one label per observation
  if ((iscellstr (Y) || isstring (Y) || iscategorical (Y)) && isrow (Y))
    Y = Y(:);
  endif
  [N, P] = size (X);
  [n, K] = size (Y);
  if (N == 1)   ## if X is a row vector, make it a column vector
    X = X(:);
    N = P;
    P = 1;
  endif
  if (n != N)
    error ("mnrfit: Y must have the same number of rows as X.")
  endif
  if (! (isnumeric (Y) || islogical (Y) || ischar (Y) || iscellstr (Y) ||
         isstring (Y) || iscategorical (Y)))
    error (strcat ("mnrfit: Response labels must be a character array,", ...
                   " a cell vector of strings, a string or categorical", ...
                   " vector, or a vector or matrix of doubles, singles or", ...
                   " logical values."));
  endif

  ## Check supplied parameters
  if (mod (numel (varargin), 2) != 0)
    error ("mnrfit: optional arguments must be in pairs.")
  endif
  MODELTYPE = 'nominal';
  DISPLAY = 'off';
  LINK = 'logit';
  ESTDISP = 'off';
  while (numel (varargin) > 0)
    name = varargin{1};
    value = varargin{2};
    switch (lower (name))
      case 'model'
        MODELTYPE = value;
      case 'display'
        DISPLAY = value;
      case 'link'
        LINK = value;
      case 'estdisp'
        ESTDISP = value;
      otherwise
        warning (sprintf ("mnrfit: parameter %s will be ignored", name));
    endswitch
    varargin(1:2) = [];
  endwhile
  LINK = lower (LINK);
  if (! any (strcmp (LINK, {'logit', 'probit', 'comploglog', 'loglog'})))
    error ("mnrfit: unrecognised 'link' value.");
  endif
  if (! (ischar (ESTDISP) && any (strcmpi (ESTDISP, {'on', 'off'}))))
    error ("mnrfit: 'estdisp' must be 'on' or 'off'.");
  endif
  estdisp = strcmpi (ESTDISP, 'on');

  ## Evaluate display input argument
  switch (lower (DISPLAY))
    case 'on'
      dispopt = true;
    case 'off'
      dispopt = false;
  endswitch

  ## Categorize Y if it holds text or categories.  Text is numbered by first
  ## appearance and a categorical by its categories, as in MATLAB.  Each row of
  ## a character array is one label, so its columns are not categories.
  istext = ischar (Y) || iscellstr (Y) || isstring (Y) || iscategorical (Y);
  if (istext)
    if (ischar (Y))
      Y = cellstr (Y);
    elseif (K > 1)
      if (iscellstr (Y))
        ytype = 'cellstr';
      else
        ytype = class (Y);
      endif
      error ("mnrfit: Y must be a column vector when given as %s.", ytype);
    endif
    ## Get groups in Y: "", <missing> and <undefined> are NaN
    YN = grp2idx (Y);
    ## Remove missing values from X and Y
    RowsUsed  = ! logical (sum (isnan ([X, YN]), 2));
    Y         = Y(RowsUsed);
    X         = X(RowsUsed, :);
    ## Renew groups in Y, as a category may have lost all its observations,
    ## and number them without gaps, as a categorical may hold unused ones
    [UY, ~, YN] = unique (grp2idx (Y));
    n = numel (UY);
  endif

  if (K > 1 && ! istext)
    ## So far, if K > 1, Y must be a matrix of logical, singles or doubles
    if (! all (all (Y == 0 | Y == 1)))
      error ("mnrfit: Y must contain only 1 and 0 when given as a 2D matrix.");
    endif
    ## Convert Y to a vector of positive integer categories
    Y = sum (bsxfun (@times, (1:K), Y), 2);
  endif

  ## Categorize Y in all other cases
  if (! istext)
    RowsUsed  = ! logical (sum (isnan ([X, Y]), 2));
    Y         = Y(RowsUsed);
    X         = X(RowsUsed, :);
    [UY, ~, YN] = unique (Y);  ## find unique categories in the response
    n = numel (UY);            ## number of unique response categories
  endif

  if (isempty (Y))
    error ("mnrfit: X and Y must contain at least one observation.");
  endif

  if (isnumeric (Y))
    if (! (all (Y > 0) && all (fix (Y) == Y)))
      error ("mnrfit: Y must contain positive integer category numbers.")
    endif
  endif

  ## Fit the requested model type
  switch (lower (MODELTYPE))
    case 'nominal'
      if (! strcmp (LINK, 'logit'))
        error ("mnrfit: nominal models support only the logit link.");
      endif
      if (n == 1)
        [B, DEV, STATS] = mnrfit_empty_ (X, estdisp);
      else
        [B, DEV, STATS] = mnrfit_nominal_ (X, YN, n, dispopt);
      endif
    case 'ordinal'
      if (n == 1)
        error ("mnrfit: an ordinal model needs at least two categories.");
      elseif (strcmp (LINK, 'logit'))
        [INTERCEPT, SLOPE, DEV, ~, ~, ~, S] = logistic_regression (YN - 1, X, ...
                                                                   dispopt);
        B = cat (1, INTERCEPT, SLOPE);
        se = S.se;
        STATS = struct ('beta', B, ...
                        'dfe', rows (X) * (n - 1) - numel (B), ...
                        's', 1, ...
                        'sfit', 1, ...
                        'estdisp', false, ...
                        'coeffcorr', S.coeffcorr, ...
                        'covb', S.cov, ...
                        'se', se, ...
                        't', B ./ se, ...
                        'p', 2 * normcdf (- abs (B ./ se)), ...
                        'resid', [], ...
                        'residp', [], ...
                        'residd', []);
      else
        [B, DEV, STATS] = mnrfit_ordinal_ (X, YN, n, LINK);
      endif
    case 'hierarchical'
      if (n == 1)
        [B, DEV, STATS] = mnrfit_empty_ (X, estdisp);
      else
        [B, DEV, STATS] = mnrfit_hierarchical_ (X, YN, n, LINK);
      endif
    otherwise
      error ("mnrfit: model type not recognised.");
  endswitch

  ## Residuals and dispersion, computed from the fitted category probabilities
  ## (common to every model type; individual responses have sample size 1)
  if (n > 1)
    pihat = mnrval (B, X, 'model', lower (MODELTYPE), 'link', LINK);
    [STATS.resid, STATS.residp, STATS.residd] = ...
                                    mnrfit_residuals_ (pihat, YN, n);
    STATS.estdisp = estdisp;
  endif
  if (! isempty (STATS.dfe) && STATS.dfe > 0)
    ## Standard Pearson dispersion estimate
    STATS.sfit = sqrt (sum (STATS.resid(:) .^ 2 ./ pihat(:)) / STATS.dfe);
    if (estdisp)
      ## Scale the coefficient covariance and standard errors by the estimated
      ## dispersion and test the coefficients against the t distribution.
      STATS.s = STATS.sfit;
      STATS.covb = STATS.covb * STATS.s ^ 2;
      STATS.se = STATS.se * STATS.s;
      STATS.t = STATS.beta ./ STATS.se;
      STATS.p = 2 * tcdf (- abs (STATS.t), STATS.dfe);
    endif
  endif

  ## Observations left out keep their rows in the residuals, as NaN
  if (! all (RowsUsed))
    for name = {'resid', 'residp', 'residd'}
      R = NaN (numel (RowsUsed), columns (STATS.(name{1})));
      R(RowsUsed,:) = STATS.(name{1});
      STATS.(name{1}) = R;
    endfor
  endif

endfunction

## Fit a baseline-category multinomial logit model by Newton-Raphson.  The last
## category is the reference; B is a (P+1)-by-(K-1) matrix whose column j holds
## the intercept and slopes contrasting category j against the reference.
function [B, dev, stats] = mnrfit_nominal_ (X, YN, k, dispopt)

  [nobs, p] = size (X);
  Z = [ones(nobs, 1), X];               ## design with intercept, nobs-by-q
  q = p + 1;
  ncat = k - 1;                         ## non-reference categories (1..k-1)

  ## Response indicator matrix for categories 1..k-1 (reference = category k)
  Yind = double (YN(:) == (1:ncat));    ## nobs-by-ncat

  ## Newton-Raphson on the multinomial logit log-likelihood
  beta = zeros (q, ncat);
  maxiter = 100;
  tol = 1e-8;
  converged = false;
  for iter = 1:maxiter
    [P, ~] = mnrfit_softmax_ (Z, beta, nobs);
    grad = Z' * (Yind - P);             ## q-by-ncat
    ## Hessian of the log-likelihood (block form, negative definite)
    H = zeros (q * ncat);
    for a = 1:ncat
      for b = 1:ncat
        if (a == b)
          w = P(:,a) .* (1 - P(:,a));
        else
          w = - P(:,a) .* P(:,b);
        endif
        H((a-1)*q+(1:q), (b-1)*q+(1:q)) = - (Z' * (Z .* w));
      endfor
    endfor
    ## A vanishing Hessian condition signals (quasi-)separable data with no
    ## finite maximum likelihood estimate; stop and keep the current estimate.
    if (rcond (H) < eps)
      break;
    endif
    step = - (H \ grad(:));             ## Newton step (maximise the likelihood)
    if (! all (isfinite (step)))
      break;
    endif
    beta(:) = beta(:) + step;
    if (max (abs (step)) < tol)
      converged = true;
      break;
    endif
  endfor
  if (! converged)
    warning ("mnrfit: iteration limit reached; results may be unreliable.");
  endif
  B = beta;

  ## Deviance = -2 * log-likelihood (saturated log-likelihood is 0 for
  ## individual responses)
  [~, Pall] = mnrfit_softmax_ (Z, beta, nobs);
  idx = sub2ind ([nobs, k], (1:nobs)', YN(:));
  dev = -2 * sum (log (Pall(idx)));

  ## Coefficient covariance from the inverse negative Hessian at the MLE
  covb = inv (-H);
  se_v = sqrt (diag (covb));
  se = reshape (se_v, q, ncat);
  coeffcorr = covb ./ (se_v * se_v');
  tstat = B ./ se;
  pval = 2 * normcdf (- abs (tstat));

  stats = struct ('beta', B, ...
                  'dfe', nobs * ncat - numel (B), ...
                  's', 1, ...
                  'sfit', 1, ...
                  'estdisp', false, ...
                  'coeffcorr', coeffcorr, ...
                  'covb', covb, ...
                  'se', se, ...
                  't', tstat, ...
                  'p', pval, ...
                  'resid', [], ...
                  'residp', [], ...
                  'residd', []);

endfunction

## Fit a hierarchical (sequential / continuation-ratio) model.  Category j is
## contrasted against the later categories using only the observations with
## Y >= j, giving K-1 independent binary logit fits.  B is a (P+1)-by-(K-1)
## matrix whose column j holds the intercept and slopes of the conditional
## logit for category j.
function [B, dev, stats] = mnrfit_hierarchical_ (X, YN, k, link)

  YN = YN(:);
  p = columns (X);
  q = p + 1;
  ncat = k - 1;
  B = zeros (q, ncat);
  se = zeros (q, ncat);
  covb = zeros (q * ncat);
  dev = 0;
  usedobs = 0;
  for j = 1:ncat
    mask = (YN >= j);
    ## Each stage is a binary GLM of category j against the later categories
    ysucc = double (YN(mask) == j);
    [Bj, devj, sj] = glmfit (X(mask, :), ysucc, 'binomial', 'link', link);
    B(:, j) = Bj;
    se(:, j) = sj.se;
    covb((j-1)*q+(1:q), (j-1)*q+(1:q)) = sj.covb;
    dev += devj;
    usedobs += sum (mask);
  endfor

  se_v = se(:);
  coeffcorr = covb ./ (se_v * se_v');
  tstat = B ./ se;
  pval = 2 * normcdf (- abs (tstat));

  stats = struct ('beta', B, ...
                  'dfe', usedobs - numel (B), ...
                  's', 1, ...
                  'sfit', 1, ...
                  'estdisp', false, ...
                  'coeffcorr', coeffcorr, ...
                  'covb', covb, ...
                  'se', se, ...
                  't', tstat, ...
                  'p', pval, ...
                  'resid', [], ...
                  'residp', [], ...
                  'residd', []);

endfunction

## Numerically stable baseline-category softmax.  Returns the nobs-by-(K-1)
## matrix P of non-reference category probabilities and, optionally, the full
## nobs-by-K matrix Pall including the reference category in the last column.
function [P, Pall] = mnrfit_softmax_ (Z, beta, nobs)
  eta = Z * beta;                       ## nobs-by-(k-1)
  m = max ([eta, zeros(nobs, 1)], [], 2);   ## row max including baseline 0
  ee = exp (eta - m);
  base = exp (-m);                      ## reference category, unnormalised
  den = base + sum (ee, 2);
  P = ee ./ den;
  if (nargout > 1)
    Pall = [P, base ./ den];
  endif
endfunction

## Raw, Pearson, and deviance residuals for individual multinomial responses.
## PIHAT is the nobs-by-K matrix of fitted category probabilities and YN holds
## the observed category indices.  resid and residp are nobs-by-K; residd is
## nobs-by-1 (the per-observation deviance contribution).
function [resid, residp, residd] = mnrfit_residuals_ (pihat, YN, k)
  nobs = rows (pihat);
  yind = double (YN(:) == (1:k));          ## observed indicators, nobs-by-K
  resid = yind - pihat;
  residp = resid ./ sqrt (pihat .* (1 - pihat));
  residp(isnan (residp)) = 0;              ## 0/0 at a saturated probability -> 0
  pobs = pihat (sub2ind ([nobs, k], (1:nobs)', YN(:)));
  residd = -2 * log (max (pobs, realmin));
endfunction

## The fit of a response with a single category, which leaves nothing to fit:
## no coefficients, a deviance of 0, and residuals of 0, as MATLAB returns it.
function [B, dev, stats] = mnrfit_empty_ (X, estdisp)
  nobs = rows (X);
  B = zeros (columns (X) + 1, 0);
  dev = 0;
  if (estdisp)
    s = NaN;
  else
    s = 1;
  endif
  stats = struct ('beta', B, ...
                  'dfe', 0, ...
                  's', s, ...
                  'sfit', NaN, ...
                  'estdisp', estdisp, ...
                  'coeffcorr', [], ...
                  'covb', [], ...
                  'se', B, ...
                  't', B, ...
                  'p', B, ...
                  'resid', zeros (nobs, 1), ...
                  'residp', zeros (nobs, 1), ...
                  'residd', zeros (nobs, 1));
endfunction

## Fit a cumulative-link (proportional-odds) ordinal model for a non-logit link
## by Newton-Raphson with numerical derivatives.  B is a (K-1+P)-by-1 vector of
## K-1 thresholds followed by the P shared slopes, matching the logit path.
function [B, dev, stats] = mnrfit_ordinal_ (X, YN, k, link)

  nobs = rows (X);
  p = columns (X);
  ncat = k - 1;
  npar = ncat + p;
  ilink = mnrfit_ilink_ (link);

  ## Initial thresholds from cumulative category frequencies, slopes at zero
  counts = accumarray (YN(:), 1, [k, 1]);
  cumfreq = cumsum (counts(1:ncat)) / nobs;
  cumfreq = min (max (cumfreq, 1e-3), 1 - 1e-3);
  params = [mnrfit_flink_(link, cumfreq)(:); zeros(p, 1)];

  ll_fn = @(pp) mnrfit_ord_ll_ (pp, X, YN, k, ilink);
  maxiter = 200;
  tol = 1e-9;
  converged = false;
  for iter = 1:maxiter
    [g, H] = mnrfit_numderiv_ (ll_fn, params);
    if (rcond (H) < eps || ! all (isfinite (g(:))))
      break;
    endif
    step = - (H \ g);
    if (! all (isfinite (step)))
      break;
    endif
    ## Damp large steps for the nonlinear links
    if (max (abs (step)) > 5)
      step *= 5 / max (abs (step));
    endif
    params += step;
    if (max (abs (step)) < tol)
      converged = true;
      break;
    endif
  endfor
  if (! converged)
    warning ("mnrfit: iteration limit reached; results may be unreliable.");
  endif

  B = params;
  dev = -2 * ll_fn (params);
  covb = inv (-H);
  se = sqrt (diag (covb));
  stats = struct ('beta', B, ...
                  'dfe', nobs * ncat - npar, ...
                  's', 1, ...
                  'sfit', 1, ...
                  'estdisp', false, ...
                  'coeffcorr', covb ./ (se * se'), ...
                  'covb', covb, ...
                  'se', se, ...
                  't', B ./ se, ...
                  'p', 2 * normcdf (- abs (B ./ se)), ...
                  'resid', [], ...
                  'residp', [], ...
                  'residd', []);

endfunction

## Log-likelihood of the cumulative-link ordinal model at parameter vector PP
## ([thresholds; slopes]); ILINK is the inverse link.
function ll = mnrfit_ord_ll_ (pp, X, YN, k, ilink)
  ncat = k - 1;
  theta = pp(1:ncat)(:).';
  beta = pp(ncat+1:end);
  gamma = ilink (X * beta + theta);          ## cumulative probabilities
  Pall = [gamma(:,1), diff(gamma, 1, 2), 1 - gamma(:,end)];
  nobs = rows (X);
  pobs = max (Pall(sub2ind ([nobs, k], (1:nobs)', YN(:))), realmin);
  ll = sum (log (pobs));
endfunction

## Central-difference gradient and Hessian of a scalar function FN at X.
function [g, H] = mnrfit_numderiv_ (fn, x)
  n = numel (x);
  h = 1e-5 * max (abs (x), 1);
  g = zeros (n, 1);
  for i = 1:n
    xp = x; xp(i) += h(i);
    xm = x; xm(i) -= h(i);
    g(i) = (fn (xp) - fn (xm)) / (2 * h(i));
  endfor
  H = zeros (n);
  for i = 1:n
    for j = i:n
      xpp = x; xpp(i) += h(i); xpp(j) += h(j);
      xpm = x; xpm(i) += h(i); xpm(j) -= h(j);
      xmp = x; xmp(i) -= h(i); xmp(j) += h(j);
      xmm = x; xmm(i) -= h(i); xmm(j) -= h(j);
      H(i,j) = (fn (xpp) - fn (xpm) - fn (xmp) + fn (xmm)) / (4 * h(i) * h(j));
      H(j,i) = H(i,j);
    endfor
  endfor
endfunction

## Inverse link (mean function) handle; matches the links used by mnrval.
function f = mnrfit_ilink_ (link)
  switch (link)
    case 'logit'
      f = @(e) 1 ./ (1 + exp (-e));
    case 'probit'
      f = @(e) normcdf (e);
    case 'comploglog'
      f = @(e) 1 - exp (-exp (e));
    case 'loglog'
      f = @(e) exp (-exp (e));
  endswitch
endfunction

## Forward link, used only to seed the initial thresholds.
function e = mnrfit_flink_ (link, g)
  switch (link)
    case 'logit'
      e = log (g ./ (1 - g));
    case 'probit'
      e = norminv (g);
    case 'comploglog'
      e = log (- log (1 - g));
    case 'loglog'
      e = log (- log (g));
  endswitch
endfunction

## Test nominal multinomial logit fitting
%!test  # nominal MLE reproduces the observed category totals
%! X = [-2; -1; 0; 1; 2; -2; -1; 0; 1; 2; -1.5; 1.5];
%! Y = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
%! [B, dev, stats] = mnrfit (X, Y, 'model', 'nominal');
%! assert_equal (size (B), [2, 2]);
%! P = mnrval (B, X);
%! assert_equal (sum (P, 2), ones (12, 1), 1e-10);
%! assert_equal (sum (P, 1), [sum(Y == 1), sum(Y == 2), sum(Y == 3)], 1e-6);

%!test  # binary nominal agrees with the ordinal cumulative-logit fit
%! X = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%! Y = [1; 1; 1; 1; 2; 1; 2; 2; 2; 2];
%! assert_equal (mnrfit (X, Y, 'model', 'nominal'), ...
%!         mnrfit (X, Y, 'model', 'ordinal'), 1e-5);

## Test hierarchical (sequential) fitting
%!test  # hierarchical fit yields valid probabilities through mnrval
%! X = [-2; -1; 0; 1; 2; -2; -1; 0; 1; 2; -1.5; 1.5];
%! Y = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
%! [B, dev, stats] = mnrfit (X, Y, 'model', 'hierarchical');
%! assert_equal (size (B), [2, 2]);
%! P = mnrval (B, X, 'model', 'hierarchical');
%! assert_equal (sum (P, 2), ones (12, 1), 1e-10);
%! assert_equal (all (P(:) >= 0 & P(:) <= 1), true);

%!test  # first hierarchical stage is the binary logit of category 1 vs. rest
%! X = [-2; -1; 0; 1; 2; -2; -1; 0; 1; 2; -1.5; 1.5];
%! Y = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
%! Bh = mnrfit (X, Y, 'model', 'hierarchical');
%! Bb = mnrfit (X, 2 - double (Y == 1), 'model', 'nominal');
%! assert_equal (Bh(:,1), Bb, 1e-8);

## Test link functions for ordinal and hierarchical models
%!test  # non-logit ordinal fits round-trip through mnrval with valid probs
%! X = [-2; -1; 0; 1; 2; -1.5; 0.5; 1.2; -0.7; 0.3; -0.4; 0.8];
%! Y = [1; 1; 2; 2; 3; 1; 2; 3; 3; 2; 1; 3];
%! for lk = {'probit', 'comploglog', 'loglog'}
%!   B = mnrfit (X, Y, 'model', 'ordinal', 'link', lk{1});
%!   P = mnrval (B, X, 'model', 'ordinal', 'link', lk{1});
%!   assert_equal (sum (P, 2), ones (12, 1), 1e-9);
%!   assert_equal (all (P(:) >= -1e-12 & P(:) <= 1 + 1e-12), true);
%! endfor

%!test  # non-logit hierarchical fit runs and round-trips through mnrval
%! X = [-2; -1; 0; 1; 2; -1.5; 0.5; 1.2; -0.7; 0.3; -0.4; 0.8];
%! Y = [1; 1; 2; 2; 3; 1; 2; 3; 3; 2; 1; 3];
%! B = mnrfit (X, Y, 'model', 'hierarchical', 'link', 'probit');
%! P = mnrval (B, X, 'model', 'hierarchical', 'link', 'probit');
%! assert_equal (sum (P, 2), ones (12, 1), 1e-9);

## Test residual outputs
%!test  # residuals have the right shape and are internally consistent
%! X = [-2; -1; 0; 1; 2; -2; -1; 0; 1; 2; -1.5; 1.5];
%! Y = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
%! [B, dev, stats] = mnrfit (X, Y, 'model', 'nominal');
%! assert_equal (size (stats.resid), [12, 3]);
%! assert_equal (size (stats.residp), [12, 3]);
%! assert_equal (size (stats.residd), [12, 1]);
%! assert_equal (sum (stats.resid, 2), zeros (12, 1), 1e-10);
%! assert_equal (sum (stats.residd), dev, 1e-8);

%!test  # dispersion estimate sfit = sqrt (Pearson X2 / dfe) matches MATLAB
%! X = [-2; -1; 0; 1; 2; -2; -1; 0; 1; 2; -1.5; 1.5];
%! Y = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
%! [~, ~, sn] = mnrfit (X, Y, 'model', 'nominal');
%! assert_equal (sn.sfit, 1.0954, 1e-3);
%! [~, ~, so] = mnrfit (X, Y, 'model', 'ordinal');
%! assert_equal (so.sfit, 1.0691, 1e-3);

## Test the EstDisp option
%!test  # EstDisp on scales se/covb by the dispersion and uses t-based p-values
%! X = [-2; -1; 0; 1; 2; -2; -1; 0; 1; 2; -1.5; 1.5];
%! Y = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
%! [~, ~, s0] = mnrfit (X, Y, 'model', 'nominal');
%! [~, ~, s1] = mnrfit (X, Y, 'model', 'nominal', 'estdisp', 'on');
%! assert_equal (s0.s, 1);
%! assert_equal (s1.s, s1.sfit);
%! assert_equal (s1.se, s0.se * s1.sfit, 1e-12);
%! assert_equal (s1.covb, s0.covb * s1.sfit ^ 2, 1e-10);
%! assert_equal (s1.p, 2 * tcdf (- abs (s1.t), s1.dfe), 1e-12);

## Test text and categorical responses against MATLAB R2026b.  Text responses
## number their categories by first appearance and categorical ones by their
## categories, the last category being the reference; a missing label drops
## its observation, and a category without observations is ignored.
%!shared X, sp, B1, Brev, Bmiss
%! load fisheriris
%! X = meas(:,1);
%! sp = species;
%! B1 = [38.7590, 12.6771; -6.8464, -2.0307];
%! Brev = [-38.7590, -26.0819; 6.8464, 4.8157];
%! Bmiss = [38.4986, 12.6788; -6.8013, -2.0310];
%!test
%! assert_equal (mnrfit (X, sp), B1, 1e-4);
%! assert_equal (mnrfit (flipud (X), flipud (sp)), Brev, 1e-4);
%!test
%! assert_equal (mnrfit (X, string (sp)), B1, 1e-4);
%! assert_equal (mnrfit (flipud (X), flipud (string (sp))), Brev, 1e-4);
%!test
%! assert_equal (mnrfit (X, char (sp)), B1, 1e-4);
%! ## labels of different lengths, padded with blanks in the char matrix
%! y = {'a'; 'bcd'; 'ef'}([2 * ones(50,1); ones(50,1); 3 * ones(50,1)]);
%! assert_equal (mnrfit (X, char (y)), mnrfit (X, y), 1e-12);
%!test
%! assert_equal (mnrfit (X, categorical (sp)), B1, 1e-4);
%! c = categorical (sp, {'virginica', 'setosa', 'versicolor'});
%! assert_equal (mnrfit (X, c), [-12.6771, 26.0819; 2.0307, -4.8157], 1e-4);
%! assert_equal (mnrfit (X, addcats (categorical (sp), 'extra')), B1, 1e-4);
%!test
%! c = categorical (sp);
%! c(1) = missing;
%! [B, dev] = mnrfit (X, c);
%! assert_equal (B, Bmiss, 1e-4);
%! assert_equal (dev, 181.6322, 1e-4);
%!test
%! s = string (sp);
%! s(1) = missing;
%! [B, dev] = mnrfit (X, s);
%! assert_equal (B, Bmiss, 1e-4);
%! assert_equal (dev, 181.6322, 1e-4);
%!test
%! c = categorical (sp, {'setosa', 'versicolor', 'virginica'}, 'Ordinal', true);
%! B = mnrfit (X, c, 'model', 'ordinal');
%! assert_equal (B, [18.5691; 21.6898; -3.4524], 1e-4);
%!test
%! c = categorical (sp, {'extra', 'setosa', 'versicolor', 'virginica'});
%! assert_equal (mnrfit (X, c), B1, 1e-4);
%!test
%! assert_equal (mnrfit (X, sp'), mnrfit (X, sp));
%! assert_equal (mnrfit (X, string (sp)'), mnrfit (X, sp));
%! assert_equal (mnrfit (X, categorical (sp)'), mnrfit (X, sp));

## Test the residuals of observations left out against MATLAB R2024a
%!test
%! c = categorical (sp);
%! c(1) = missing;
%! [~, ~, s] = mnrfit (X, c);
%! assert_equal (size (s.resid), [150, 3]);
%! assert_equal (s.resid(1:2,:), [NaN, NaN, NaN; 0.084568478433509675, ...
%!               -0.079375981712630134, -0.0051924967208795622], -1e-10);
%! assert_equal (s.residp(1:2,:), [NaN, NaN, NaN; 0.3039424270334104, ...
%!               -0.29363201424385138, -0.072246795596986726], -1e-10);
%! assert_equal (s.residd(1:2), [NaN; 0.17671943321891109], -1e-10);
%!test
%! Xn = X;
%! Xn(5) = NaN;
%! [~, ~, s] = mnrfit (Xn, grp2idx (sp));
%! assert_equal (size (s.resid), [150, 3]);
%! assert_equal (s.resid(4:6,:), [0.020809847448850061, ...
%!               -0.020095179210191657, -0.00071466823865832662; ...
%!               NaN, NaN, NaN; 0.52496547239191615, ...
%!               -0.44466304411132351, -0.080302428280592628], -1e-10);

## Test a response with a single category against MATLAB R2024a
%!test
%! [B, dev, s] = mnrfit (X(1:50), sp(1:50));
%! assert_equal (B, zeros (2, 0));
%! assert_equal (dev, 0);
%! assert_equal ([s.dfe, s.s, s.sfit], [0, 1, NaN]);
%! assert_equal (s.resid, zeros (50, 1));
%! assert_equal (s.residd, zeros (50, 1));
%!test
%! [B, dev] = mnrfit (X(1:50), sp(1:50), 'model', 'hierarchical');
%! assert_equal (B, zeros (2, 0));
%! assert_equal (dev, 0);
%!test
%! [~, ~, s] = mnrfit (X(1:50), sp(1:50), 'estdisp', 'on');
%! assert_equal ([s.s, s.sfit], [NaN, NaN]);
%!test
%! y = sp(1:50);
%! y{2} = '';
%! [~, ~, s] = mnrfit (X(1:50), y);
%! assert_equal (s.resid(1:3), [0; NaN; 0]);

## Test input validation
%!error<mnrfit: too few input arguments.> mnrfit (ones (50,1))
%!error<mnrfit: Predictors must be numeric.> ...
%! mnrfit ({1 ;2 ;3 ;4 ;5}, ones (5,1))
%!error<mnrfit: Predictors must be a vector or a 2D matrix.> ...
%! mnrfit (ones (50, 4, 2), ones (50, 1))
%!error<mnrfit: Response must be a vector or a 2D matrix.> ...
%! mnrfit (ones (50, 4), ones (50, 1, 3))
%!error<mnrfit: Y must have the same number of rows as X.> ...
%! mnrfit (ones (50, 4), ones (45,1))
%!error<mnrfit: Response labels must be a character array, a cell vector> ...
%! mnrfit (ones (5, 4), {1 ;2 ;3 ;4 ;5})
%!error<mnrfit: optional arguments must be in pairs.> ...
%! mnrfit (ones (5, 4), ones (5, 1), 'model')
%!error<mnrfit: Y must be a column vector when given as cellstr.> ...
%! mnrfit (ones (5, 4), {'q','q';'w','w';'q','q';'w','w';'q','q'})
%!error<mnrfit: Y must be a column vector when given as categorical.> ...
%! mnrfit (ones (5, 4), categorical ([1, 2; 1, 2; 1, 2; 1, 2; 1, 2]))
%!error<mnrfit: Y must be a column vector when given as string.> ...
%! mnrfit (ones (5, 4), string ({'q','q';'w','w';'q','q';'w','w';'q','q'}))
%!error<mnrfit: Y must contain only 1 and 0 when given as a 2D matrix.> ...
%! mnrfit (ones (5, 4), [1, 2; 1, 2; 1, 2; 1, 2; 1, 2])
%!error<mnrfit: Y must contain positive integer category numbers.> ...
%! mnrfit (ones (5, 4), [1; -1; 1; 2; 1])
%!error<mnrfit: model type not recognised.> ...
%! mnrfit (ones (5, 4), [1; 2; 3; 2; 1], 'model', 'whatever')
%!error<mnrfit: unrecognised 'link' value.> ...
%! mnrfit (ones (5, 4), [1; 2; 1; 2; 1], 'link', 'cauchit')
%!error<mnrfit: nominal models support only the logit link.> ...
%! mnrfit (ones (5, 4), [1; 2; 1; 2; 1], 'model', 'nominal', 'link', 'probit')
%!error<mnrfit: 'estdisp' must be 'on' or 'off'.> ...
%! mnrfit (ones (5, 4), [1; 2; 1; 2; 1], 'estdisp', 'maybe')
%!error<mnrfit: an ordinal model needs at least two categories.> ...
%! mnrfit ((1:5)', {'a'; 'a'; 'a'; 'a'; 'a'}, 'model', 'ordinal')
%!error<mnrfit: Y must have the same number of rows as X.> ...
%! mnrfit ((1:5)', [1, 2, 1, 2, 1])

## Edge cases with empty arrays
%!error <mnrfit: X and Y must contain at least one observation.> ...
%! mnrfit ([], [])
%!error <mnrfit: X and Y must contain at least one observation.> ...
%! mnrfit (zeros (0, 3), zeros (0, 1))
