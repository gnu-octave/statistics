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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
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
## with a cumulative logit model, and hierarchical models with a sequential
## (continuation-ratio) logit.  Only the logit link is currently available.
##
## @code{@var{B} = mnrfit (@var{X}, @var{Y})}  returns a matrix, @var{B}, of
## coefficient estimates for a multinomial logistic regression of the nominal
## responses in @var{Y} on the predictors in @var{X}.  @var{X} is an @math{NxP}
## numeric matrix the observations on predictor variables, where @math{N}
## corresponds to the number of observations and @math{P} corresponds to
## predictor variables.  @var{Y} contains the response category labels and it
## either be an @math{NxP} categorical or numerical matrix (containing only 1s
## and 0s) or an @math{Nx1} numeric vector with positive integer values, a cell
## array of character vectors and a logical vector.  @var{Y} can also be defined
## as a character matrix with each row corresponding to an observation of
## @var{X}.
##
## @code{@var{B} = mnrfit (@var{X}, @var{Y}, @var{name}, @var{value})} returns a
## matrix, @var{B}, of coefficient estimates for a multinomial model fit with
## additional parameters specified @qcode{Name-Value} pair arguments.
##
## @multitable @columnfractions 0.18 0.02 0.8
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{'model'} @tab @tab The type of model to fit: @qcode{'nominal'}
## (default) for a baseline-category model, @qcode{'ordinal'} for a cumulative
## logit model, or @qcode{'hierarchical'} for a sequential logit model.
##
## @item @qcode{'display'} @tab @tab A flag to enable/disable displaying
## information about the fitted model.  Default is @qcode{'off'}.
## @end multitable
##
## @code{[@var{B}, @var{dev}, @var{stats}] = mnrfit (@dots{})} also returns the
## deviance of the fit, @var{dev}, and a structure @var{stats} with the fitted
## coefficients @qcode{'beta'} (same as @var{B}), their standard errors
## @qcode{'se'}, covariance matrix @qcode{'covb'}, and correlation matrix
## @qcode{'coeffcorr'}.  For nominal models @var{stats} also includes the error
## degrees of freedom @qcode{'dfe'} and the coefficient @math{t} statistics
## @qcode{'t'} and @math{p}-values @qcode{'p'}.
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
  if (! (isnumeric (Y) || islogical (Y) || ischar (Y) || iscellstr (Y)))
    error (strcat ("mnrfit: Response labels must be a character array,", ...
                   " a cell vector of strings, \nor a vector or", ...
                   " matrix of doubles, singles or logical values."));
  endif

  ## Check supplied parameters
  if (mod (numel (varargin), 2) != 0)
    error ("mnrfit: optional arguments must be in pairs.")
  endif
  MODELTYPE = 'nominal';
  DISPLAY = 'off';
  while (numel (varargin) > 0)
    name = varargin{1};
    value = varargin{2};
    switch (lower (name))
      case 'model'
        MODELTYPE = value;
      case 'display'
        DISPLAY = value;
      otherwise
        warning (sprintf ("mnrfit: parameter %s will be ignored", name));
    endswitch
    varargin(1:2) = [];
  endwhile

  ## Evaluate display input argument
  switch (lower (DISPLAY))
    case 'on'
      dispopt = true;
    case 'off'
      dispopt = false;
  endswitch

  ## Categorize Y if it is a cellstring array
  if (iscellstr (Y))
    if (K > 1)
      error ("mnrfit: Y must be a column vector when given as cellstr.");
    endif
    ## Get groups in Y
    [YN, ~, UY] = grp2idx (Y);  # this will also catch "" as missing values
    ## Remove missing values from X and Y
    RowsUsed  = ! logical (sum (isnan ([X, YN]), 2));
    Y         = Y(RowsUsed);
    X         = X(RowsUsed, :);
    ## Renew groups in Y
    [YN, ~, UY] = grp2idx (Y);  # in case a category is removed due to NaNs in X
    n = numel (UY);
  endif

  ## Categorize Y if it is a character array
  if (ischar (Y))
    ## Get groups in Y
    [YN, ~, UY] = grp2idx (Y);  # this will also catch "" as missing values
    ## Remove missing values from X and Y
    RowsUsed  = ! logical (sum (isnan ([X, YN]), 2));
    Y         = Y(RowsUsed);
    X         = X(RowsUsed, :);
    ## Renew groups in Y
    [YN, ~, UY] = grp2idx (Y);  # in case a category is removed due to NaNs in X
    n = numel (UY);
  endif

  if (K > 1)
    ## So far, if K > 1, Y must be a matrix of logical, singles or doubles
    if (! all (all (Y == 0 | Y == 1)))
      error ("mnrfit: Y must contain only 1 and 0 when given as a 2D matrix.");
    endif
    ## Convert Y to a vector of positive integer categories
    Y = sum (bsxfun (@times, (1:K), Y), 2);
  endif

  ## Categorize Y in all other cases
  if (! iscellstr (Y))
    RowsUsed  = ! logical (sum (isnan ([X, Y]), 2));
    Y         = Y(RowsUsed);
    X         = X(RowsUsed, :);
    [UY, ~, YN] = unique (Y);  ## find unique categories in the response
    n = numel (UY);            ## number of unique response categories
  endif

  if (isnumeric (Y))
    if (! (all (Y > 0) && all (fix (Y) == Y)))
      error ("mnrfit: Y must contain positive integer category numbers.")
    endif
  endif

  ## Fit the requested model type
  switch (lower (MODELTYPE))
    case 'nominal'
      [B, DEV, STATS] = mnrfit_nominal_ (X, YN, n, dispopt);
    case 'ordinal'
      [INTERCEPT, SLOPE, DEV, ~, ~, ~, S] = logistic_regression (YN - 1, X, ...
                                                                 dispopt);
      B = cat (1, INTERCEPT, SLOPE);
      STATS = struct ('beta', B, ...
                      'dfe', [], ...      ## Not used
                      's', [], ...        ## Not used
                      'sfit', [], ...     ## Not used
                      'estdisp', [], ...  ## Not used
                      'coeffcorr', S.coeffcorr, ...
                      'covb', S.cov, ...
                      'se', S.se, ...
                      't', [], ...        ## Not used
                      'p', [], ...        ## Not used
                      'resid', [], ...    ## Not used
                      'residp', [], ...   ## Not used
                      'residd', []);      ## Not used
    case 'hierarchical'
      [B, DEV, STATS] = mnrfit_hierarchical_ (X, YN, n, dispopt);
    otherwise
      error ("mnrfit: model type not recognised.");
  endswitch

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
function [B, dev, stats] = mnrfit_hierarchical_ (X, YN, k, dispopt)

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
    ## Category j is class 1 (the reference class 2 is the later categories)
    YNsub = 2 - double (YN(mask) == j);
    [Bj, devj, sj] = mnrfit_nominal_ (X(mask, :), YNsub, 2, dispopt);
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

## Test nominal multinomial logit fitting
%!test  # nominal MLE reproduces the observed category totals
%! X = [-2; -1; 0; 1; 2; -2; -1; 0; 1; 2; -1.5; 1.5];
%! Y = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
%! [B, dev, stats] = mnrfit (X, Y, 'model', 'nominal');
%! assert (size (B), [2, 2]);
%! P = mnrval (B, X);
%! assert (sum (P, 2), ones (12, 1), 1e-10);
%! assert (sum (P, 1), [sum(Y == 1), sum(Y == 2), sum(Y == 3)], 1e-6);

%!test  # binary nominal agrees with the ordinal cumulative-logit fit
%! X = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%! Y = [1; 1; 1; 1; 2; 1; 2; 2; 2; 2];
%! assert (mnrfit (X, Y, 'model', 'nominal'), ...
%!         mnrfit (X, Y, 'model', 'ordinal'), 1e-5);

## Test hierarchical (sequential) fitting
%!test  # hierarchical fit yields valid probabilities through mnrval
%! X = [-2; -1; 0; 1; 2; -2; -1; 0; 1; 2; -1.5; 1.5];
%! Y = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
%! [B, dev, stats] = mnrfit (X, Y, 'model', 'hierarchical');
%! assert (size (B), [2, 2]);
%! P = mnrval (B, X, 'model', 'hierarchical');
%! assert (sum (P, 2), ones (12, 1), 1e-10);
%! assert (all (P(:) >= 0 & P(:) <= 1));

%!test  # first hierarchical stage is the binary logit of category 1 vs. rest
%! X = [-2; -1; 0; 1; 2; -2; -1; 0; 1; 2; -1.5; 1.5];
%! Y = [1; 2; 3; 1; 2; 3; 1; 2; 3; 1; 2; 3];
%! Bh = mnrfit (X, Y, 'model', 'hierarchical');
%! Bb = mnrfit (X, 2 - double (Y == 1), 'model', 'nominal');
%! assert (Bh(:,1), Bb, 1e-8);

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
%!error<mnrfit: Y must contain only 1 and 0 when given as a 2D matrix.> ...
%! mnrfit (ones (5, 4), [1, 2; 1, 2; 1, 2; 1, 2; 1, 2])
%!error<mnrfit: Y must contain positive integer category numbers.> ...
%! mnrfit (ones (5, 4), [1; -1; 1; 2; 1])
%!error<mnrfit: model type not recognised.> ...
%! mnrfit (ones (5, 4), [1; 2; 3; 2; 1], 'model', 'whatever')

