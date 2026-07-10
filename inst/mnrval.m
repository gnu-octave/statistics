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
## @deftypefn  {statistics} {@var{pihat} =} mnrval (@var{B}, @var{X})
## @deftypefnx {statistics} {[@var{pihat}, @var{dlo}, @var{dhi}] =} mnrval (@var{B}, @var{X}, @var{stats})
## @deftypefnx {statistics} {@var{yhat} =} mnrval (@var{B}, @var{X}, @var{ssize})
## @deftypefnx {statistics} {[@var{yhat}, @var{dlo}, @var{dhi}] =} mnrval (@var{B}, @var{X}, @var{ssize}, @var{stats})
## @deftypefnx {statistics} {[@dots{}] =} mnrval (@dots{}, @var{name}, @var{value})
##
## Predict values for a multinomial logistic regression model.
##
## @code{@var{pihat} = mnrval (@var{B}, @var{X})} returns the predicted
## category probabilities @var{pihat} of a multinomial logistic regression with
## coefficients @var{B}, evaluated at the predictor values in @var{X}.  @var{X}
## is an @math{NxP} numeric matrix of @math{N} observations on @math{P}
## predictors.  @var{pihat} is an @math{NxK} matrix, where @math{K} is the number
## of response categories and each row sums to one.  @var{B} is the coefficient
## matrix returned by @code{mnrfit} (see below for its shape under each model).
##
## @code{mnrval} is the prediction companion of @code{mnrfit}.  Unlike the
## current @code{mnrfit}, which only fits ordinal and two-category nominal
## models, @code{mnrval} evaluates all three model types, so a coefficient
## matrix @var{B} obtained elsewhere (e.g.@: MATLAB) can be used for prediction.
##
## @code{@var{yhat} = mnrval (@var{B}, @var{X}, @var{ssize})} returns predicted
## category counts instead of probabilities, for the sample sizes in @var{ssize}
## (a scalar or an @math{Nx1} vector).
##
## @code{[@var{pihat}, @var{dlo}, @var{dhi}] = mnrval (@var{B}, @var{X},
## @var{stats})} also returns @math{95%} confidence bounds on the predictions.
## @var{stats} is the structure returned by @code{mnrfit}; its @qcode{'covb'}
## field (the coefficient covariance matrix) is required.  The confidence
## interval for each prediction is @code{[@var{pihat} - @var{dlo}, @var{pihat} +
## @var{dhi}]}.  The bounds are nonsimultaneous and apply to the fitted values,
## not to new observations.
##
## The following @qcode{Name-Value} pairs control the model:
##
## @multitable @columnfractions 0.18 0.8
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'model'} @tab The model type: @qcode{'nominal'} (default),
## @qcode{'ordinal'}, or @qcode{'hierarchical'}.
##
## @item @qcode{'interactions'} @tab @qcode{'on'} to include category-specific
## coefficients, or @qcode{'off'} for a common set of coefficients with
## category-specific intercepts only.  Default is @qcode{'on'} for nominal and
## hierarchical models and @qcode{'off'} for ordinal models.  With
## @qcode{'interactions','on'}, @var{B} is a @math{(P+1)x(K-1)} matrix.  With
## @qcode{'interactions','off'}, @var{B} is a @math{(K-1+P)x1} vector holding the
## @math{K-1} intercepts followed by the @math{P} common slopes.
##
## @item @qcode{'link'} @tab The link function for ordinal and hierarchical
## models: @qcode{'logit'} (default), @qcode{'probit'}, @qcode{'comploglog'}, or
## @qcode{'loglog'}.  Nominal models always use the multinomial logit link.
##
## @item @qcode{'type'} @tab The kind of probability returned:
## @qcode{'category'} (default, @math{NxK} category probabilities),
## @qcode{'cumulative'} (@math{Nx(K-1)} cumulative probabilities of the first
## @math{K-1} categories), or @qcode{'conditional'} (@math{Nx(K-1)} conditional
## probabilities of each category given membership in that or a later category).
##
## @item @qcode{'confidence'} @tab The confidence level for @var{dlo} and
## @var{dhi}, a scalar in the range @math{(0,1)}.  Default is @math{0.95}.
## @end multitable
##
## @seealso{mnrfit, glmval, logistic_regression}
## @end deftypefn

function [pihat, dlo, dhi] = mnrval (B, X, varargin)

  if (nargin < 2)
    error ("mnrval: too few input arguments.");
  endif
  if (! isnumeric (B) || ! isreal (B))
    error ("mnrval: B must be a real numeric matrix.");
  endif
  if (! isnumeric (X) || ! isreal (X) || ndims (X) > 2)
    error ("mnrval: X must be a real numeric 2D matrix.");
  endif

  ## Peel off the optional positional arguments (SSIZE and/or STATS), which
  ## precede any Name-Value pairs.  SSIZE is numeric, STATS is a struct.
  ssize = [];
  stats = [];
  args = varargin;
  while (numel (args) > 0 && ! ischar (args{1}))
    a = args{1};
    if (isstruct (a))
      stats = a;
    elseif (isnumeric (a))
      ssize = a;
    else
      error ("mnrval: invalid optional positional argument.");
    endif
    args(1) = [];
  endwhile

  ## Parse Name-Value pairs
  if (mod (numel (args), 2) != 0)
    error ("mnrval: optional arguments must be in Name-Value pairs.");
  endif
  model = 'nominal';
  interactions = [];              ## resolved to a default per model below
  link = 'logit';
  type = 'category';
  conf = 0.95;
  while (numel (args) > 0)
    name = args{1};
    value = args{2};
    switch (lower (name))
      case 'model'
        model = value;
      case 'interactions'
        interactions = value;
      case 'link'
        link = value;
      case 'type'
        type = value;
      case 'confidence'
        conf = value;
      otherwise
        error ("mnrval: unknown parameter name '%s'.", name);
    endswitch
    args(1:2) = [];
  endwhile

  ## Validate and normalise the options
  model = lower (model);
  if (! any (strcmp (model, {'nominal', 'ordinal', 'hierarchical'})))
    error ("mnrval: unrecognised 'model' value.");
  endif
  if (isempty (interactions))
    if (strcmp (model, 'ordinal'))
      interactions = 'off';
    else
      interactions = 'on';
    endif
  endif
  interactions = lower (interactions);
  if (! any (strcmp (interactions, {'on', 'off'})))
    error ("mnrval: 'interactions' must be 'on' or 'off'.");
  endif
  link = lower (link);
  if (! any (strcmp (link, {'logit', 'probit', 'comploglog', 'loglog'})))
    error ("mnrval: unrecognised 'link' value.");
  endif
  if (strcmp (model, 'nominal') && ! strcmp (link, 'logit'))
    error ("mnrval: nominal models use the multinomial logit link only.");
  endif
  type = lower (type);
  if (! any (strcmp (type, {'category', 'cumulative', 'conditional'})))
    error ("mnrval: unrecognised 'type' value.");
  endif
  if (! (isscalar (conf) && isreal (conf) && conf > 0 && conf < 1))
    error ("mnrval: 'confidence' must be a scalar in the range (0,1).");
  endif

  [n, P] = size (X);

  ## Determine K-1 and validate the shape of B against X and 'interactions'
  if (strcmp (interactions, 'on'))
    if (rows (B) != P + 1)
      error (strcat ("mnrval: with 'interactions','on', B must have one", ...
                     " more row than the number of columns of X."));
    endif
    Km1 = columns (B);
  else
    B = B(:);
    if (numel (B) <= P)
      error (strcat ("mnrval: with 'interactions','off', B must have more", ...
                     " elements than the number of columns of X."));
    endif
    Km1 = numel (B) - P;
  endif
  if (Km1 < 1)
    error ("mnrval: B implies fewer than two response categories.");
  endif

  ## Validate SSIZE and STATS if supplied
  if (! isempty (ssize))
    if (! (isscalar (ssize) || (isvector (ssize) && numel (ssize) == n)))
      error ("mnrval: SSIZE must be a scalar or an N-element vector.");
    endif
    ssize = ssize(:);
  endif
  if (nargout > 1 && isempty (stats))
    error ("mnrval: STATS is required to compute confidence bounds.");
  endif

  ## Point predictions
  pihat = mnr_predict (B, X, model, interactions, link, type, ssize, P);

  ## Confidence bounds by the delta method on the coefficient covariance
  if (nargout > 1)
    if (! isfield (stats, 'covb'))
      error ("mnrval: STATS must contain a 'covb' field.");
    endif
    covb = stats.covb;
    b0 = B(:);
    np = numel (b0);
    if (! isequal (size (covb), [np, np]))
      error ("mnrval: size of STATS.covb does not match the number of coefficients.");
    endif
    ncol = columns (pihat);
    ## Central-difference Jacobian of the returned quantity with respect to B
    J = zeros (n * ncol, np);
    h = 1e-6 * max (abs (b0), 1);
    for j = 1:np
      bp = b0;  bp(j) += h(j);
      bm = b0;  bm(j) -= h(j);
      fp = mnr_predict (reshape (bp, size (B)), X, model, interactions, ...
                        link, type, ssize, P);
      fm = mnr_predict (reshape (bm, size (B)), X, model, interactions, ...
                        link, type, ssize, P);
      J(:,j) = (fp(:) - fm(:)) / (2 * h(j));
    endfor
    v = sum ((J * covb) .* J, 2);         ## variance of each output element
    se = reshape (sqrt (max (v, 0)), n, ncol);
    z = norminv (1 - (1 - conf) / 2);
    dlo = z * se;
    dhi = z * se;
  endif

endfunction

## Core prediction: returns the requested quantity (category/cumulative/
## conditional probabilities, optionally scaled to counts by SSIZE).
function out = mnr_predict (B, X, model, interactions, link, type, ssize, P)

  n = rows (X);

  ## Linear predictors ETA (n-by-(K-1))
  if (strcmp (interactions, 'on'))
    Km1 = columns (B);
    eta = [ones(n, 1), X] * B;
  else
    Km1 = numel (B) - P;
    icept = B(1:Km1)(:).';
    slope = B(Km1+1:end);
    if (P > 0)
      eta = X * slope + repmat (icept, n, 1);
    else
      eta = repmat (icept, n, 1);
    endif
  endif

  ## Inverse link (mean function) for ordinal and hierarchical models
  switch (link)
    case 'logit'
      G = 1 ./ (1 + exp (-eta));
    case 'probit'
      G = normcdf (eta);
    case 'comploglog'
      G = 1 - exp (-exp (eta));
    case 'loglog'
      G = exp (-exp (eta));
  endswitch

  ## Category probabilities PCAT (n-by-K)
  switch (model)
    case 'nominal'
      ee = exp (eta);
      den = 1 + sum (ee, 2);
      pcat = [ee ./ den, 1 ./ den];
    case 'ordinal'
      pcat = [G(:,1), diff(G, 1, 2), 1 - G(:,end)];
    case 'hierarchical'
      ## G holds the conditional probabilities P(Y = j | Y >= j)
      surv = cumprod ([ones(n, 1), 1 - G], 2);    ## surv(:,j) = P(Y >= j)
      pcat = [G .* surv(:,1:Km1), surv(:,end)];
  endswitch

  ## Requested probability type
  switch (type)
    case 'category'
      out = pcat;
    case 'cumulative'
      c = cumsum (pcat, 2);
      out = c(:,1:end-1);
    case 'conditional'
      surv = fliplr (cumsum (fliplr (pcat), 2));   ## surv(:,j) = P(Y >= j)
      out = pcat(:,1:end-1) ./ surv(:,1:end-1);
  endswitch

  ## Scale to counts when sample sizes are supplied
  if (! isempty (ssize))
    out = out .* ssize;
  endif

endfunction

%!demo
%! ## Fit an ordinal model and predict the category probabilities
%! X = [1; 2; 3; 4; 5; 6; 7; 8];
%! Y = [1; 1; 1; 2; 2; 2; 3; 3];
%! B = mnrfit (X, Y, 'model', 'ordinal');
%! pihat = mnrval (B, X, 'model', 'ordinal')

## Round-trip against logistic_regression: mnrval must reproduce the fitted
## probabilities of the ordinal model that mnrfit wraps.
%!test
%! X = [1.489381332449196, 1.1534152241851305; ...
%!      1.8110085304863965, 0.9449666896938425; ...
%!      -0.04453299665130296, 0.34278203449678646; ...
%!      -0.36616019468850347, 1.130254275908322; ...
%!       0.15339143291005095, -0.7921044310668951; ...
%!      -1.6031878794469698, -1.8343471035233376; ...
%!      -0.14349521143198166, -0.6762996896828459; ...
%!      -0.4403818557740143, -0.7921044310668951; ...
%!      -0.7372685001160434, -0.027793137932169563; ...
%!      -0.11875465773681024, 0.5512305689880763];
%! Y = [1;1;1;1;1;0;0;0;0;0];
%! B = mnrfit (X, Y + 1, 'model', 'ordinal');
%! [~, ~, ~, ~, ~, Pref] = logistic_regression (Y, X, false);
%! assert (mnrval (B, X, 'model', 'ordinal'), Pref, 1e-10);

## For a binary response the default nominal call equals the ordinal call
%!test
%! X = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%! Y = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];
%! B = mnrfit (X, Y, 'model', 'ordinal');
%! p_nom = mnrval (B, X);
%! p_ord = mnrval (B, X, 'model', 'ordinal');
%! assert (p_nom, p_ord, 1e-12);

## Every row of the category probabilities sums to one, for all three models
%!test
%! X = randn (20, 2);
%! Bnom = [0.5, -0.2; 1.0, 0.3; -0.4, 0.1];      # (P+1)-by-(K-1), K = 3
%! assert (sum (mnrval (Bnom, X), 2), ones (20, 1), 1e-12);
%! assert (sum (mnrval (Bnom, X, 'model', 'hierarchical'), 2), ...
%!         ones (20, 1), 1e-12);
%! Bord = [-1; 1; 0.5; -0.3];                     # (K-1+P)-by-1, K = 3, P = 2
%! assert (sum (mnrval (Bord, X, 'model', 'ordinal'), 2), ...
%!         ones (20, 1), 1e-12);

## Probabilities lie in [0, 1]
%!test
%! X = randn (15, 2);
%! Bord = [-1; 0.8; 0.5; -0.3];
%! p = mnrval (Bord, X, 'model', 'ordinal');
%! assert (all (p(:) >= 0 & p(:) <= 1));

## Cumulative type equals the running sum of the category probabilities
%!test
%! X = randn (12, 2);
%! Bord = [-0.5; 1.2; 0.4; -0.2];
%! pc = mnrval (Bord, X, 'model', 'ordinal', 'type', 'category');
%! cu = mnrval (Bord, X, 'model', 'ordinal', 'type', 'cumulative');
%! assert (cu, cumsum (pc(:,1:end-1), 2), 1e-12);

## Predicted counts equal probabilities times the sample size
%!test
%! X = randn (10, 2);
%! Bord = [-0.5; 1.0; 0.4; -0.2];
%! p = mnrval (Bord, X, 'model', 'ordinal');
%! y = mnrval (Bord, X, 100, 'model', 'ordinal');
%! assert (y, p * 100, 1e-10);

## Confidence bounds have the same size as the prediction and are non-negative
%!test
%! X = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%! Y = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];
%! [B, ~, stats] = mnrfit (X, Y, 'model', 'ordinal');
%! [p, dlo, dhi] = mnrval (B, X, stats, 'model', 'ordinal');
%! assert (size (dlo), size (p));
%! assert (size (dhi), size (p));
%! assert (all (dlo(:) >= 0) && all (dhi(:) >= 0));

## Every link produces category probabilities that sum to one
%!test
%! X = randn (10, 1);
%! Bord = [-0.5; 0.7; 0.4];
%! for lnk = {'logit', 'probit', 'comploglog', 'loglog'}
%!   p = mnrval (Bord, X, 'model', 'ordinal', 'link', lnk{1});
%!   assert (sum (p, 2), ones (10, 1), 1e-10);
%! endfor

## Increasing links with ordered intercepts keep probabilities in [0, 1]
%!test
%! X = randn (10, 1);
%! Bord = [-0.5; 0.7; 0.4];
%! for lnk = {'logit', 'probit', 'comploglog'}
%!   p = mnrval (Bord, X, 'model', 'ordinal', 'link', lnk{1});
%!   assert (all (p(:) >= 0 & p(:) <= 1));
%! endfor

## The loglog link is decreasing, so a valid ordinal model needs its
## intercepts in decreasing order
%!test
%! X = randn (10, 1);
%! Bll = [1.0; -0.5; 0.4];
%! p = mnrval (Bll, X, 'model', 'ordinal', 'link', 'loglog');
%! assert (sum (p, 2), ones (10, 1), 1e-10);
%! assert (all (p(:) >= 0 & p(:) <= 1));

## Test input validation
%!error<mnrval: too few input arguments.> mnrval (1)
%!error<mnrval: B must be a real numeric matrix.> mnrval ({1}, ones (3, 1))
%!error<mnrval: X must be a real numeric 2D matrix.> ...
%! mnrval ([1; 1], ones (3, 3, 3))
%!error<mnrval: optional arguments must be in Name-Value pairs.> ...
%! mnrval ([1; 1], ones (3, 1), 'model')
%!error<mnrval: unknown parameter name 'foo'.> ...
%! mnrval ([1; 1], ones (3, 1), 'foo', 'bar')
%!error<mnrval: unrecognised 'model' value.> ...
%! mnrval ([1; 1], ones (3, 1), 'model', 'whatever')
%!error<mnrval: 'interactions' must be 'on' or 'off'.> ...
%! mnrval ([1; 1], ones (3, 1), 'interactions', 'maybe')
%!error<mnrval: unrecognised 'link' value.> ...
%! mnrval ([1; 1], ones (3, 1), 'model', 'ordinal', 'link', 'foo')
%!error<mnrval: nominal models use the multinomial logit link only.> ...
%! mnrval ([1; 1], ones (3, 1), 'link', 'probit')
%!error<mnrval: unrecognised 'type' value.> ...
%! mnrval ([1; 1], ones (3, 1), 'type', 'foo')
%!error<mnrval: 'confidence' must be a scalar in the range .0,1..> ...
%! mnrval ([1; 1], ones (3, 1), 'confidence', 2)
%!error<mnrval: with 'interactions','on', B must have one more row> ...
%! mnrval ([1; 1; 1], ones (3, 1))
%!error<mnrval: STATS is required to compute confidence bounds.> ...
%! [p, dlo] = mnrval ([1; 1], ones (3, 1));
%!error<mnrval: STATS must contain a 'covb' field.> ...
%! [p, dlo] = mnrval ([1; 1], ones (3, 1), struct ('x', 1));
