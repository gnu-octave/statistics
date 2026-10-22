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
## can be an @math{N*1} vector of positive integers, logical values, character
## vectors in a cell array, strings, or categorical values, each element the
## category of one observation.  A cell array, string, or categorical @var{Y}
## can also be a row vector.  @var{Y} can also be a character matrix with each
## row holding the label of an observation of @var{X}.
##
## @var{Y} can also be an @math{N*K} numeric or logical matrix of counts, with
## a column for each of the @math{K} categories: @code{@var{Y}(i,j)} is the
## number of observations in category @math{j} at the predictor values in
## @code{@var{X}(i,:)}, so that a row of @var{X} stands for the
## @code{sum (@var{Y}(i,:))} observations it was seen with.  The counts must be
## nonnegative integers.  A row of zeros stands for no observations and takes
## no part in the estimates, though it counts in @qcode{'dfe'}, as in MATLAB.
## A matrix of 1s and 0s with a single 1 in each row is the same as a vector
## of one category per observation.
##
## The categories of a numeric or logical @var{Y} are taken in sorted order,
## those of a cell array, string, or character matrix in order of first
## appearance, and those of a categorical @var{Y} in the order of its
## categories.  A category without observations is ignored wherever it lies
## among the categories of a categorical @var{Y}; MATLAB ignores only a last one
## and fits coefficients for any other.  Likewise the categories of a numeric
## @var{Y} are the values it holds, where MATLAB numbers them from 1 to the
## largest and fits coefficients for any value missing below it, and a column
## of zeros in a matrix @var{Y} is no category, where MATLAB fits coefficients
## for it.  The last
## category is the reference category of a nominal model.  An observation whose
## label is missing (an empty character vector, @qcode{<missing>}, or
## @qcode{<undefined>}), whose row of a matrix @var{Y} holds a @qcode{NaN}, or
## whose predictors hold a @qcode{NaN}, is left out of the fit, and its rows of
## the residuals in @var{stats} are @qcode{NaN}.
##
## A @var{Y} with a single category leaves nothing to fit: a nominal or
## hierarchical model returns a @var{B} with no columns, a deviance of 0 and
## residuals of 0, the raw and Pearson residuals of a hierarchical model having
## no columns, and an ordinal model is an error.
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
## @item @qcode{'interactions'} @tab @qcode{'on'} to give each category its
## own slopes, or @qcode{'off'} for slopes shared by every category, the
## intercepts always differing.  The default is @qcode{'on'} for nominal and
## hierarchical models and @qcode{'off'} for ordinal ones.  With
## @qcode{'on'}, @var{B} is a @math{(P+1)*(K-1)} matrix whose column @math{j}
## holds the intercept and slopes of category @math{j}; with @qcode{'off'},
## @var{B} is a @math{(K-1+P)*1} vector of the @math{K-1} intercepts followed by
## the @math{P} shared slopes.  A nominal model without interactions has no
## slopes to share, so it fits the @math{K-1} intercepts alone, with a warning.
##
## @item @qcode{'weights'} @tab An @math{N*1} vector of nonnegative
## observation weights, each multiplying the counts of its row of @var{Y}.
## The default is a vector of ones.  A row whose weight is @qcode{NaN} is left
## out of the fit.
##
## @item @qcode{'iterationlimit'} @tab The largest number of iterations of the
## fit, a positive integer.  The default is 100.
##
## @item @qcode{'tolerance'} @tab The convergence tolerance of the fit, a
## positive scalar: the fit stops once no coefficient changes by more than this
## fraction of its value.  The default is @code{1e-6}.
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
## The raw and Pearson residuals are taken on what the model fits: an
## @math{N*K} matrix over the categories of a nominal model, and an
## @math{N*(K-1)} matrix over the cumulative probabilities of an ordinal model
## or over the stages of a hierarchical one, where a row is @qcode{NaN} at the
## stages that none of its observations reach.  A raw residual is the count
## minus its fitted value, the count times the fitted probability, and a
## Pearson residual divides it by the binomial standard deviation of the count;
## a weight multiplies both counts.  A row of zeros has a Pearson residual of
## @qcode{NaN}.  The deviance residuals are an @math{N*1} vector holding the
## contribution of each row to the deviance, so they sum to @var{dev}; for an
## observation of one category it is @math{-2} times the log of the fitted
## probability of that category.  For a hierarchical model MATLAB sums the
## signed deviance residuals of the stages instead, which is @qcode{NaN} for
## every observation in the first category.  @qcode{'sfit'} is the square root
## of the Pearson statistic over the categories divided by @qcode{'dfe'}, and it
## differs from MATLAB's, which for ordinal and hierarchical models is taken
## over the residuals in @qcode{'residp'}.
##
## @qcode{'dfe'} is the number of rows of @var{X} times @math{K-1}, less the
## number of coefficients; for a hierarchical model each stage counts only the
## rows with observations that reach it.
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
  INTERACTIONS = [];            ## resolved to a default per model below
  WEIGHTS = [];
  MAXITER = 100;
  TOL = 1e-6;
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
      case 'interactions'
        INTERACTIONS = value;
      case 'weights'
        WEIGHTS = value;
      case 'iterationlimit'
        MAXITER = value;
      case 'tolerance'
        TOL = value;
      otherwise
        warning (sprintf ("mnrfit: parameter %s will be ignored", name));
    endswitch
    varargin(1:2) = [];
  endwhile
  MODELTYPE = lower (MODELTYPE);
  if (! any (strcmp (MODELTYPE, {'nominal', 'ordinal', 'hierarchical'})))
    error ("mnrfit: model type not recognised.");
  endif
  LINK = lower (LINK);
  if (! any (strcmp (LINK, {'logit', 'probit', 'comploglog', 'loglog'})))
    error ("mnrfit: unrecognised 'link' value.");
  endif
  if (! (ischar (ESTDISP) && any (strcmpi (ESTDISP, {'on', 'off'}))))
    error ("mnrfit: 'estdisp' must be 'on' or 'off'.");
  endif
  estdisp = strcmpi (ESTDISP, 'on');
  if (isempty (INTERACTIONS))
    if (strcmp (MODELTYPE, 'ordinal'))
      INTERACTIONS = 'off';
    else
      INTERACTIONS = 'on';
    endif
  endif
  if (! (ischar (INTERACTIONS) && any (strcmpi (INTERACTIONS, {'on', 'off'}))))
    error ("mnrfit: 'interactions' must be 'on' or 'off'.");
  endif
  interact = strcmpi (INTERACTIONS, 'on');
  if (! (isnumeric (MAXITER) && isscalar (MAXITER) && MAXITER >= 1
                             && fix (MAXITER) == MAXITER))
    error ("mnrfit: 'iterationlimit' must be a positive integer.");
  endif
  if (! (isnumeric (TOL) && isscalar (TOL) && TOL > 0))
    error ("mnrfit: 'tolerance' must be a positive scalar.");
  endif
  if (isempty (WEIGHTS))
    W = ones (N, 1);
  else
    if (! (isnumeric (WEIGHTS) && isreal (WEIGHTS) && isvector (WEIGHTS)
                               && numel (WEIGHTS) == N
                               && all (WEIGHTS(! isnan (WEIGHTS)) >= 0)))
      error (strcat ("mnrfit: 'weights' must be a vector of nonnegative", ...
                     " values with one element per observation."));
    endif
    W = double (WEIGHTS(:));
  endif

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
  iscount = K > 1 && ! istext;
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
    Ymiss = isnan (YN);
  elseif (iscount)
    ## A matrix Y holds the number of observations in each category for each
    ## row of X; a row holding a NaN is missing
    YN = double (Y);
    Yok = isnan (YN) | (YN >= 0 & fix (YN) == YN);
    if (! all (Yok(:)))
      error (strcat ("mnrfit: Y must contain nonnegative integer counts", ...
                     " when given as a 2D matrix."));
    endif
    Ymiss = any (isnan (YN), 2);
  else
    YN = double (Y);
    Ymiss = isnan (YN);
  endif

  ## Remove missing values from X, Y and the weights
  RowsUsed = ! (any (isnan (X), 2) | Ymiss | isnan (W));
  X = X(RowsUsed,:);
  W = W(RowsUsed);
  YN = YN(RowsUsed,:);
  if (isempty (YN))
    error ("mnrfit: X and Y must contain at least one observation.");
  endif
  if (isnumeric (Y) && ! iscount)
    if (! (all (YN > 0) && all (fix (YN) == YN)))
      error ("mnrfit: Y must contain positive integer category numbers.")
    endif
  endif

  ## Hold the response as counts, one column per category.  A category without
  ## observations is ignored, as a categorical may hold unused ones and a
  ## category may have lost all its observations.
  if (iscount)
    C = YN(:, any (YN > 0, 1));
  else
    [~, ~, YN] = unique (YN);
    C = double (YN(:) == (1:max (YN)));
  endif
  n = columns (C);
  if (n == 0)
    error ("mnrfit: X and Y must contain at least one observation.");
  endif

  ## Fit the requested model type
  switch (MODELTYPE)
    case 'nominal'
      if (! strcmp (LINK, 'logit'))
        error ("mnrfit: nominal models support only the logit link.");
      endif
      if (n == 1)
        STATS = mnrfit_empty_ (X, estdisp, 1);
      else
        if (! interact)
          warning (strcat ("mnrfit: a nominal model without interactions", ...
                           " fits the intercepts alone."));
        endif
        STATS = mnrfit_nominal_ (X, C, W, interact, MAXITER, TOL);
      endif
    case 'ordinal'
      if (n == 1)
        error ("mnrfit: an ordinal model needs at least two categories.");
      endif
      STATS = mnrfit_ordinal_ (X, C, W, LINK, interact, MAXITER, TOL);
    case 'hierarchical'
      if (n == 1)
        STATS = mnrfit_empty_ (X, estdisp, 0);
      else
        STATS = mnrfit_hierarchical_ (X, C, W, LINK, interact, MAXITER, TOL);
      endif
  endswitch
  B = STATS.beta;

  ## Residuals, deviance and dispersion, computed from the fitted category
  ## probabilities (common to every model type)
  DEV = 0;
  if (n > 1)
    pihat = mnrfit_pihat_ (B, X, MODELTYPE, interact, LINK);
    [STATS.resid, STATS.residp, STATS.residd] = ...
                    mnrfit_residuals_ (pihat, C, W, MODELTYPE);
    DEV = sum (STATS.residd);
    STATS.estdisp = estdisp;
  endif
  if (! isempty (STATS.dfe) && STATS.dfe > 0)
    ## Standard Pearson dispersion estimate, over the categories
    e = sum (C, 2) .* pihat;
    X2 = W .* (C - e) .^ 2 ./ e;
    STATS.sfit = sqrt (sum (X2(e > 0)) / STATS.dfe);
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
## the intercept and slopes contrasting category j against the reference, or
## without interactions a (K-1)-by-1 vector of intercepts alone.  C holds the
## counts of each category and W multiplies them.
function stats = mnrfit_nominal_ (X, C, W, interact, maxiter, tol)

  nobs = rows (X);
  if (interact)
    Z = [ones(nobs, 1), X];             ## design with intercept, nobs-by-q
  else
    Z = ones (nobs, 1);
  endif
  q = columns (Z);
  ncat = columns (C) - 1;               ## non-reference categories (1..k-1)
  wC = W .* C(:,1:ncat);                ## weighted counts, nobs-by-ncat
  wn = W .* sum (C, 2);                 ## weighted trials, nobs-by-1

  ## Newton-Raphson on the multinomial logit log-likelihood
  beta = zeros (q, ncat);
  converged = false;
  for iter = 1:maxiter
    [grad, H] = mnrfit_nominal_derivs_ (Z, beta, wC, wn);
    ## A vanishing Hessian condition signals (quasi-)separable data with no
    ## finite maximum likelihood estimate; stop and keep the current estimate.
    if (rcond (H) < eps)
      break;
    endif
    step = - (H \ grad(:));             ## Newton step (maximise the likelihood)
    if (! all (isfinite (step)))
      break;
    endif
    old = beta(:);
    beta(:) = old + step;
    if (all (abs (step) <= tol * max (sqrt (eps), abs (old))))
      converged = true;
      break;
    endif
  endfor
  if (! converged)
    warning ("mnrfit: iteration limit reached; results may be unreliable.");
  endif
  [~, H] = mnrfit_nominal_derivs_ (Z, beta, wC, wn);
  if (interact)
    B = beta;
  else
    B = beta(:);
  endif

  ## Coefficient covariance from the inverse negative Hessian at the MLE
  covb = inv (-H);
  se_v = sqrt (diag (covb));
  se = reshape (se_v, size (B));
  tstat = B ./ se;
  stats = mnrfit_stats_ (B, nobs * ncat - numel (B), covb, se, tstat);

endfunction

## Gradient and Hessian of the multinomial logit log-likelihood at BETA
function [grad, H] = mnrfit_nominal_derivs_ (Z, beta, wC, wn)
  [nobs, q] = size (Z);
  ncat = columns (beta);
  P = mnrfit_softmax_ (Z, beta, nobs);
  grad = Z' * (wC - wn .* P);           ## q-by-ncat
  ## Hessian of the log-likelihood (block form, negative definite)
  H = zeros (q * ncat);
  for a = 1:ncat
    for b = 1:ncat
      if (a == b)
        w = P(:,a) .* (1 - P(:,a));
      else
        w = - P(:,a) .* P(:,b);
      endif
      H((a-1)*q+(1:q), (b-1)*q+(1:q)) = - (Z' * (Z .* (wn .* w)));
    endfor
  endfor
endfunction

## Fit a sequential (continuation-ratio) model.  Stage j is a binary model of
## category j against the later categories, over the trials that reach it.
## With interactions each stage has its own coefficients, and B is a
## (P+1)-by-(K-1) matrix; without, the stages share the slopes and B is a
## (K-1+P)-by-1 vector of K-1 intercepts followed by the P slopes.
function stats = mnrfit_hierarchical_ (X, C, W, link, interact, maxiter, tol)

  [nobs, p] = size (X);
  ncat = columns (C) - 1;
  ## The trials reaching each stage, those not in an earlier category
  m = sum (C, 2) - [zeros(nobs, 1), cumsum(C(:,1:ncat-1), 2)];
  opts = struct ('MaxIter', maxiter, 'TolX', tol);
  if (interact)
    q = p + 1;
    B = zeros (q, ncat);
    se = zeros (q, ncat);
    covb = zeros (q * ncat);
    for j = 1:ncat
      r = m(:,j) > 0;
      [Bj, ~, sj] = glmfit (X(r,:), [C(r,j), m(r,j)], 'binomial', ...
                            'link', link, 'weights', W(r), 'options', opts);
      B(:,j) = Bj;
      se(:,j) = sj.se;
      covb((j-1)*q+(1:q), (j-1)*q+(1:q)) = sj.covb;
    endfor
  else
    ## One binary model over every stage, with an intercept for each stage
    [ii, jj] = find (m > 0);
    idx = sub2ind (size (m), ii, jj);
    Z = [double(jj == (1:ncat)), X(ii,:)];
    [B, ~, sj] = glmfit (Z, [C(idx), m(idx)], 'binomial', 'link', link, ...
                         'constant', 'off', 'weights', W(ii), 'options', opts);
    se = sj.se;
    covb = sj.covb;
  endif
  tstat = B ./ se;
  stats = mnrfit_stats_ (B, sum (m(:) > 0) - numel (B), covb, se, tstat);

endfunction

## Fit a cumulative link model by Fisher scoring.  Without interactions the
## categories share the slopes and B is a (K-1+P)-by-1 vector of K-1
## thresholds followed by the P slopes; with interactions each cumulative
## probability has its own slopes and B is a (P+1)-by-(K-1) matrix.
function stats = mnrfit_ordinal_ (X, C, W, link, interact, maxiter, tol)

  [nobs, p] = size (X);
  ncat = columns (C) - 1;
  [F, f] = mnrfit_links_ (link);

  ## The derivative of each cumulative linear predictor by the parameters
  if (interact)
    npar = ncat * (p + 1);
  else
    npar = ncat + p;
  endif
  D = cell (1, ncat);
  for j = 1:ncat
    D{j} = zeros (nobs, npar);
    if (interact)
      D{j}(:,(j-1)*(p+1)+(1:p+1)) = [ones(nobs, 1), X];
    else
      D{j}(:,j) = 1;
      D{j}(:,ncat+1:end) = X;
    endif
  endfor

  ## Initial thresholds from the cumulative category frequencies, slopes at 0
  tot = sum (W .* C, 1);
  cumfreq = cumsum (tot(1:ncat)) / sum (tot);
  cumfreq = min (max (cumfreq, 1e-3), 1 - 1e-3);
  theta = mnrfit_flink_ (link, cumfreq)(:);
  if (interact)
    psi = reshape ([theta'; zeros(p, ncat)], [], 1);
  else
    psi = [theta; zeros(p, 1)];
  endif

  [ll, g, I] = mnrfit_ordinal_score_ (psi, D, C, W, F, f);
  converged = false;
  for iter = 1:maxiter
    if (rcond (I) < eps || ! all (isfinite (g)))
      break;
    endif
    step = I \ g;
    if (! all (isfinite (step)))
      break;
    endif
    ## Halve a step that leaves a probability outside (0, 1) or lowers the
    ## likelihood
    for halving = 1:30
      [ll_new, g_new, I_new] = mnrfit_ordinal_score_ (psi + step, D, C, ...
                                                      W, F, f);
      if (ll_new >= ll - 1e-10 * abs (ll))
        break;
      endif
      step /= 2;
    endfor
    if (ll_new < ll - 1e-10 * abs (ll))
      break;
    endif
    old = psi;
    psi = old + step;
    ll = ll_new;
    g = g_new;
    I = I_new;
    if (all (abs (step) <= tol * max (sqrt (eps), abs (old))))
      converged = true;
      break;
    endif
  endfor
  if (! converged)
    warning ("mnrfit: iteration limit reached; results may be unreliable.");
  endif

  if (interact)
    B = reshape (psi, p + 1, ncat);
  else
    B = psi;
  endif
  covb = inv (I);
  se = reshape (sqrt (diag (covb)), size (B));
  tstat = B ./ se;
  stats = mnrfit_stats_ (B, nobs * ncat - npar, covb, se, tstat);

endfunction

## Log-likelihood, score and expected information of a cumulative link model
## at the parameters PSI.  D{j} is the derivative of the j-th cumulative linear
## predictor by the parameters.  A probability outside (0, 1) gives a
## log-likelihood of -Inf.
function [ll, g, I] = mnrfit_ordinal_score_ (psi, D, C, W, F, f)
  [nobs, k] = size (C);
  ncat = k - 1;
  npar = numel (psi);
  eta = zeros (nobs, ncat);
  for j = 1:ncat
    eta(:,j) = D{j} * psi;
  endfor
  gamma = F (eta);                       ## cumulative probabilities
  P = diff ([zeros(nobs, 1), gamma, ones(nobs, 1)], 1, 2);
  g = zeros (npar, 1);
  I = zeros (npar);
  if (! all (P(:) > 0))
    ll = -Inf;
    return;
  endif
  wC = W .* C;
  ll = sum (wC(wC > 0) .* log (P(wC > 0)));
  wn = W .* sum (C, 2);
  dens = f (eta);
  Aprev = zeros (nobs, npar);
  for c = 1:k
    if (c < k)
      A = dens(:,c) .* D{c};
    else
      A = zeros (nobs, npar);
    endif
    dP = A - Aprev;                      ## derivative of category c
    g += dP' * (wC(:,c) ./ P(:,c));
    I += dP' * (dP .* (wn ./ P(:,c)));
    Aprev = A;
  endfor
endfunction

## The STATS structure of a fit, before its residuals and dispersion
function stats = mnrfit_stats_ (B, dfe, covb, se, tstat)
  se_v = se(:);
  stats = struct ('beta', B, ...
                  'dfe', dfe, ...
                  's', 1, ...
                  'sfit', 1, ...
                  'estdisp', false, ...
                  'coeffcorr', covb ./ (se_v * se_v'), ...
                  'covb', covb, ...
                  'se', se, ...
                  't', tstat, ...
                  'p', 2 * normcdf (- abs (tstat)), ...
                  'resid', [], ...
                  'residp', [], ...
                  'residd', []);
endfunction

## Softmax over the K-1 linear predictors and the reference category (fixed at
## 0).  P is nobs-by-(K-1) for categories 1..K-1; Pall appends the reference.
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

## Fitted category probabilities of the model B on the predictors X
function pihat = mnrfit_pihat_ (B, X, model, interact, link)
  if (strcmp (model, 'nominal') && ! interact)
    ## The intercepts alone: every observation has the same probabilities
    e = exp ([B(:)', 0] - max ([B(:)', 0]));
    pihat = repmat (e / sum (e), rows (X), 1);
  else
    if (interact)
      onoff = 'on';
    else
      onoff = 'off';
    endif
    pihat = mnrval (B, X, 'model', model, 'link', link, 'interactions', onoff);
  endif
endfunction

## Raw, Pearson, and deviance residuals.  PIHAT is the nobs-by-K matrix of
## fitted category probabilities, C the counts of each category, and W
## multiplies them.  resid and residp are taken on what MODEL fits:
## nobs-by-K over the categories of a nominal model, nobs-by-(K-1) over the
## cumulative probabilities of an ordinal one and over the stages of a
## hierarchical one, NaN at a stage no trials reach.  residd is nobs-by-1, the
## contribution of each observation to the deviance.
function [resid, residp, residd] = mnrfit_residuals_ (pihat, C, W, model)
  [nobs, k] = size (C);
  wC = W .* C;
  wn = W .* sum (C, 2);
  switch (model)
    case 'nominal'
      obs = wC;
      fit = pihat;
      m = repmat (wn, 1, k);
    case 'ordinal'
      ## The last cumulative probability is always 1, so it has no residual
      obs = cumsum (wC(:,1:k-1), 2);
      fit = cumsum (pihat(:,1:k-1), 2);
      m = repmat (wn, 1, k - 1);
    case 'hierarchical'
      ## Stage j is category j against the later ones, over the trials that
      ## reach it
      obs = wC(:,1:k-1);
      fit = pihat(:,1:k-1) ./ ...
            (1 - [zeros(nobs, 1), cumsum(pihat(:,1:k-2), 2)]);
      m = wn - [zeros(nobs, 1), cumsum(wC(:,1:k-2), 2)];
  endswitch
  resid = obs - m .* fit;
  residp = resid ./ sqrt (m .* fit .* (1 - fit));
  residp(isnan (residp)) = 0;              ## 0/0 at a saturated probability -> 0
  residp(m == 0) = NaN;                    ## no trials, no residual
  if (strcmp (model, 'hierarchical'))
    resid(m == 0) = NaN;
  endif
  d = wC .* log (wC ./ max (wn .* pihat, realmin));
  d(wC == 0) = 0;
  residd = 2 * sum (d, 2);
endfunction

## The fit of a response with a single category, which leaves nothing to fit:
## no coefficients, a deviance of 0, and residuals of 0, as MATLAB returns it.
## NCOL is the number of columns of the raw and Pearson residuals: one for the
## single category, or none for the stages of a hierarchical model.
function stats = mnrfit_empty_ (X, estdisp, ncol)
  nobs = rows (X);
  B = zeros (columns (X) + 1, 0);
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
                  'resid', zeros (nobs, ncol), ...
                  'residp', zeros (nobs, ncol), ...
                  'residd', zeros (nobs, 1));
endfunction

## Inverse link (cumulative probability) of the linear predictor and its
## derivative, matching the links used by mnrval.
function [F, f] = mnrfit_links_ (link)
  switch (link)
    case 'logit'
      F = @(e) 1 ./ (1 + exp (-e));
      f = @(e) 1 ./ (1 + exp (-e)) .* (1 - 1 ./ (1 + exp (-e)));
    case 'probit'
      F = @(e) normcdf (e);
      f = @(e) normpdf (e);
    case 'comploglog'
      F = @(e) 1 - exp (-exp (e));
      f = @(e) exp (e - exp (e));
    case 'loglog'
      F = @(e) exp (-exp (e));
      f = @(e) - exp (e - exp (e));
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

## Test ordinal and hierarchical residuals against MATLAB R2024a
%!test
%! [~, ~, s] = mnrfit ((1:9)', [1 1 2 1 2 3 2 3 3]', 'model', 'ordinal');
%! assert_equal (s.resid(3:6,:), [-0.640017262401357, 0.0119172052798397; ...
%!               0.662144725149099, 0.0403306249191799; -0.12773091630328, ...
%!               0.12773091630328; -0.04033062491918, -0.662144725149098], ...
%!               1e-7);
%!test
%! [~, ~, s] = mnrfit ((1:9)', [1 1 2 1 2 3 2 3 3]', 'model', 'ordinal');
%! assert_equal (s.residp(3:6,:), [-1.33338328387092, 0.109822301711886; ...
%!               1.39994543884798, 0.205001320144029; -0.382668552795833, ...
%!               0.382668552795833; -0.205001320144029, -1.39994543884798], ...
%!               1e-7);
%!test
%! [~, ~, s] = mnrfit ((1:9)', [1 1 2 1 2 3 2 3 3]', 'model', 'hierarchical');
%! assert_equal (s.resid(1:5,:), [0.041296779926102, NaN; ...
%!               0.132552307459812, NaN; -0.648478668920583, ...
%!               0.0180292890394411; 0.657879876493689, NaN; ...
%!               -0.127852690165849, 0.156132139470677], 1e-10);
%!test
%! [~, ~, s] = mnrfit ((1:9)', [1 1 2 1 2 3 2 3 3]', 'model', 'hierarchical');
%! assert_equal (s.residp(1:5,:), [0.207546780834513, NaN; ...
%!               0.390905695445191, NaN; -1.35822589269996, ...
%!               0.135500230327833; 1.38670471123006, NaN; ...
%!               -0.382877647171603, 0.430139109910011], 1e-10);
%!test
%! [~, dev, s] = mnrfit ((1:9)', [1 1 2 1 2 3 2 3 3]', 'model', 'hierarchical');
%! assert_equal (sum (s.residd), dev, 1e-10);
%!test
%! [~, ~, s] = mnrfit ((1:9)', [1 1 2 1 2 2 1 2 2]', 'model', 'ordinal');
%! assert_equal (size (s.residp), [9, 1]);
%!test
%! [~, ~, s] = mnrfit ((1:9)', [1 1 2 1 2 2 1 2 2]', 'model', 'hierarchical');
%! assert_equal (size (s.residp), [9, 1]);
%!test
%! [~, ~, s] = mnrfit ([1 NaN 3:9]', [1 1 2 1 2 3 2 3 3]', ...
%!                     'model', 'hierarchical');
%! assert_equal (s.resid(1:3,:), [0.0823310009353325, NaN; NaN, NaN; ...
%!               -0.566038550416822, 0.0180292890394411], 1e-10);

## Test an indicator matrix response with a NaN row against MATLAB R2024a
%!test
%! Yi = double ([1 1 2 1 2 3 2 3 3]' == 1:3);
%! Bk = mnrfit ([1:2, 4:9]', Yi([1:2, 4:9],:));
%! Yi(3,:) = NaN;
%! assert_equal (mnrfit ((1:9)', Yi), Bk);
%!test
%! Yi = double ([1 1 2 1 2 3 2 3 3]' == 1:3);
%! Bk = mnrfit ([1:2, 4:9]', Yi([1:2, 4:9],:));
%! Yi(3,1) = NaN;
%! assert_equal (mnrfit ((1:9)', Yi), Bk);
%!test
%! Yi = double ([1 1 2 1 2 3 2 3 3]' == 1:3);
%! Yi(3,:) = NaN;
%! [~, ~, s] = mnrfit ((1:9)', Yi);
%! assert_equal (s.dfe, 12);
%!test
%! Yi = double ([1 1 2 1 2 3 2 3 3]' == 1:3);
%! Yi(3,:) = NaN;
%! [~, ~, s] = mnrfit ((1:9)', Yi);
%! assert_equal (s.resid(3,:), [NaN, NaN, NaN]);

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
%! [~, ~, s] = mnrfit (X(1:50), sp(1:50), 'model', 'hierarchical');
%! assert_equal (s.resid, zeros (50, 0));
%!test
%! [~, ~, s] = mnrfit (X(1:50), sp(1:50), 'model', 'hierarchical');
%! assert_equal (s.residd, zeros (50, 1));
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
%!error<mnrfit: Y must contain nonnegative integer counts when given as a 2D matrix.> ...
%! mnrfit ((1:3)', [1, 2.5; 1, 2; 2, 1])
%!error<mnrfit: Y must contain nonnegative integer counts when given as a 2D matrix.> ...
%! mnrfit ((1:3)', [1, -2; 1, 2; 2, 1])
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

## Test count matrices, weights and interactions against MATLAB R2024a
%!shared X6, Yc, x, y, wy
%! X6 = (1:6)';
%! Yc = [3, 1, 0; 2, 2, 1; 1, 2, 2; 0, 2, 3; 1, 0, 4; 0, 1, 1];
%! x = [0.2, 1.1, 1.9, 3.2, 3.8, 5.1, 5.9, 7.2, 8.1, 8.8, 1.5, 4.4, 6.6, ...
%!      2.7, 9.3]';
%! y = [1, 1, 2, 1, 2, 2, 3, 2, 3, 3, 2, 3, 2, 1, 3]';
%! wy = [1, 2, 1, 1, 3, 1, 2, 1, 1, 2, 1, 1, 3, 1, 2]';
%!test
%! B = mnrfit (X6, Yc);
%! assert_equal (B, [3.11521381121846, 1.72595553252332; ...
%!                   -1.14291919345809, -0.557616323768001], -1e-10);
%!test
%! [~, dev, s] = mnrfit (X6, Yc);
%! assert_equal ([dev, s.dfe], [7.29816119718767, 8], 1e-10);
%!test
%! [~, ~, s] = mnrfit (X6, Yc);
%! assert_equal (s.resid(1:2,:), [0.479032510699473, -0.12827209523121, ...
%!               -0.350760415468263; -0.232219981795573, ...
%!               0.206186037368912, 0.0260339444266604], 1e-10);
%!test
%! [~, ~, s] = mnrfit (X6, Yc);
%! assert_equal (s.residp(1:2,:), [0.496160860213629, ...
%!               -0.142522560672393, -0.620060541159726; ...
%!               -0.208905708865958, 0.192247718709891, ...
%!               0.0293977673685418], 1e-10);
%!test
%! [~, ~, s] = mnrfit (X6, Yc);
%! assert_equal (s.residd, [0.802442532667939; 0.0485724721972971; ...
%!               0.109354412004761; 1.36850192992456; 3.70158681524588; ...
%!               1.26770303514723], 1e-10);
%!test
%! B = mnrfit (X6, Yc, 'model', 'ordinal');
%! assert_equal (B, [1.29939112375683; 3.04351150827255; ...
%!                   -0.797767453558435], -1e-6);
%!test
%! [~, dev, s] = mnrfit (X6, Yc, 'model', 'ordinal');
%! assert_equal ([dev, s.dfe], [7.34011610110504, 9], 1e-10);
%!test
%! [~, ~, s] = mnrfit (X6, Yc, 'model', 'ordinal');
%! assert_equal (s.resid(1:2,:), [0.508636704845893, 0.382868827000934; ...
%!               -0.132502184597763, -0.0484341779702326], 1e-6);
%!test
%! B = mnrfit (X6, Yc, 'model', 'hierarchical');
%! assert_equal (B, [1.4485769644502, 1.98445711156869; ...
%!                   -0.847495943829511, -0.629407372230075], -1e-10);
%!test
%! [~, dev, s] = mnrfit (X6, Yc, 'model', 'hierarchical');
%! assert_equal ([dev, s.dfe], [7.01155523089269, 8], 1e-10);
%!test
%! [~, ~, s] = mnrfit (X6, Yc, 'model', 'hierarchical');
%! assert_equal (s.resid(1:2,:), [0.416385649011944, 0.205046027351086; ...
%!               -0.193530522479006, -0.0215448777590356], 1e-10);
%!test
%! Xe = repelem (X6, sum (Yc, 2));
%! ye = repelem (repmat ((1:3)', 6, 1), reshape (Yc', [], 1));
%! assert_equal (mnrfit (X6, Yc), mnrfit (Xe, ye), 1e-10);
%!test
%! Yz = Yc;
%! Yz(2,:) = 0;
%! [~, dev, s] = mnrfit (X6, Yz);
%! assert_equal ([dev, s.dfe], [7.22221543051311, 8], 1e-10);
%!test
%! Yz = Yc;
%! Yz(2,:) = 0;
%! [~, ~, s] = mnrfit (X6, Yz);
%! assert_equal (s.residp(2,:), [NaN, NaN, NaN]);
%!test
%! Y2 = [1, 0, 0; 0, 1, 0; 1, 1, 0; 0, 0, 1; 0, 1, 0; 0, 0, 1];
%! assert_equal (mnrfit (X6, Y2), [7.55936380401886, 5.12548052211949; ...
%!               -2.20811162886661, -1.12377716339686], -1e-10);
%!test
%! Yu = [Yc(:,1), zeros(6,1), Yc(:,2:3)];
%! assert_equal (mnrfit (X6, Yu), mnrfit (X6, Yc));
%!test
%! Yr = Yc;
%! Yr(1,:) = [2, 0, 0];
%! [~, ~, s] = mnrfit (X6, Yr, 'model', 'hierarchical');
%! assert_equal (s.dfe, 7);
%!test
%! Yr = Yc;
%! Yr(1,:) = [2, 0, 0];
%! [~, ~, s] = mnrfit (X6, Yr, 'model', 'hierarchical');
%! assert_equal (s.residp(1,:), [0.914729298681304, NaN], 1e-10);
%!test
%! B = mnrfit (X6, Yc, 'weights', [1, 2, 1, 3, 1, 2]');
%! assert_equal (B, [3.26764731943353, 1.37782734242591; ...
%!                   -1.31700484949595, -0.429541892229879], -1e-10);
%!test
%! [~, ~, s] = mnrfit (X6, Yc, 'weights', [1, 2, 1, 3, 1, 2]');
%! assert_equal (s.residp(1,:), [0.369677091482785, 0.0317726223417427, ...
%!               -0.645010798142362], 1e-10);
%!test
%! [~, dev] = mnrfit (X6, Yc, 'weights', [1, 2, 1, 3, 1, 2]');
%! assert_equal (dev, 10.2534960048411, 1e-10);
%!test
%! B = mnrfit (X6, Yc, 'model', 'ordinal', 'interactions', 'on');
%! assert_equal (B, [1.29967633179764, 3.0431268128545; ...
%!                   -0.797892338817538, -0.797669756913824], -1e-6);
%!test
%! B = mnrfit (X6, Yc, 'model', 'hierarchical', 'interactions', 'off');
%! assert_equal (B, [1.18480730583952; 2.39950317498638; ...
%!                   -0.745002595688337], -1e-10);
%!test
%! [~, ~, s] = mnrfit (X6, Yc, 'model', 'hierarchical', 'interactions', 'off');
%! assert_equal (s.dfe, 9);
%!test
%! S = warning ('off', 'all');
%! B = mnrfit (x, y, 'interactions', 'off');
%! warning (S);
%! assert_equal (B, [log(4/5); log(6/5)], 1e-10);
%!warning<mnrfit: a nominal model without interactions fits the intercepts alone.> ...
%! mnrfit (x, y, 'interactions', 'off');
%!test
%! B = mnrfit (x, y, 'model', 'ordinal', 'link', 'probit');
%! assert_equal (B, [1.07863009080549; 3.04642129916916; ...
%!                   -0.485126029810669], -1e-6);
%!test
%! B = mnrfit (x, y, 'model', 'hierarchical', 'link', 'comploglog');
%! assert_equal (B, [1.11606258630227, 2.15097058822628; ...
%!                   -0.720736258358282, -0.42477426118437], -1e-10);
%!test
%! B = mnrfit (x, y, 'model', 'ordinal', 'interactions', 'on', ...
%!             'link', 'loglog');
%! assert_equal (B, [-1.64348124603758, -4.37698558082175; ...
%!                   0.536188854191686, 0.608581711015173], -1e-6);
%!test
%! B = mnrfit (x, y, 'weights', wy);
%! assert_equal (B, [7.59751957764148, 5.03769507149626; ...
%!                   -1.88801982284467, -0.77411323264876], -1e-10);
%!test
%! B = mnrfit (x, y, 'model', 'ordinal', 'weights', wy);
%! assert_equal (B, [2.05249214200963; 6.21644116058936; ...
%!                   -0.945758358953641], -1e-6);
%!test
%! B = mnrfit (x, y, 'model', 'hierarchical', 'weights', wy);
%! assert_equal (B, [2.64758153523882, 5.0288032858106; ...
%!                   -1.17838215436985, -0.771428066114368], -1e-10);
%!test
%! B = mnrfit (X6, Yc, 'tolerance', 1e-2);
%! assert_equal (B, mnrfit (X6, Yc), 1e-4);
%!warning<mnrfit: iteration limit reached; results may be unreliable.> ...
%! mnrfit (X6, Yc, 'iterationlimit', 2);
%!error<mnrfit: 'interactions' must be 'on' or 'off'.> ...
%! mnrfit (X6, Yc, 'interactions', true)
%!error<mnrfit: 'weights' must be a vector of nonnegative values with one element per observation.> ...
%! mnrfit (X6, Yc, 'weights', ones (5, 1))
%!error<mnrfit: 'weights' must be a vector of nonnegative values with one element per observation.> ...
%! mnrfit (X6, Yc, 'weights', [1; 1; -1; 1; 1; 1])
%!error<mnrfit: 'iterationlimit' must be a positive integer.> ...
%! mnrfit (X6, Yc, 'iterationlimit', 2.5)
%!error<mnrfit: 'tolerance' must be a positive scalar.> ...
%! mnrfit (X6, Yc, 'tolerance', 0)
