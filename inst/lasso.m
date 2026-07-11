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
## @deftypefn  {statistics} {@var{B} =} lasso (@var{X}, @var{y})
## @deftypefnx {statistics} {[@var{B}, @var{FitInfo}] =} lasso (@var{X}, @var{y})
## @deftypefnx {statistics} {[@dots{}] =} lasso (@dots{}, @var{Name}, @var{Value})
##
## Lasso and elastic-net regularized least-squares regression.
##
## @code{@var{B} = lasso (@var{X}, @var{y})} fits a series of regularized linear
## models of the response @var{y} on the predictor matrix @var{X} by lasso, over
## a sequence of values of the regularization parameter @var{Lambda}.  @var{B}
## is a @math{P*L} matrix whose column @math{k} holds the coefficient estimates
## for the @math{k}-th @var{Lambda}, in ascending order of @var{Lambda}.
##
## @code{[@var{B}, @var{FitInfo}] = lasso (@dots{})} additionally returns a
## structure @var{FitInfo} with fields @code{Intercept}, @code{Lambda},
## @code{Alpha}, @code{DF} (number of non-zero coefficients), and @code{MSE}
## (mean squared error), one entry per value of @var{Lambda}.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.18 0.82
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'Alpha'} @tab The elastic-net mixing parameter in
## @math{(0, 1]}.  @math{1} (default) is the lasso penalty; smaller values add a
## ridge penalty.
##
## @item @qcode{'Lambda'} @tab A vector of non-negative regularization
## parameters.  By default a geometric sequence of @qcode{'NumLambda'} values is
## used, from the smallest value that drives all coefficients to zero down to
## @qcode{'LambdaRatio'} times that value.
##
## @item @qcode{'NumLambda'} @tab The number of @var{Lambda} values in the
## default sequence (default @math{100}).
##
## @item @qcode{'LambdaRatio'} @tab The ratio of the smallest to the largest
## @var{Lambda} in the default sequence (default @math{1e-4}, or @math{1e-2}
## when the number of observations is below the number of predictors).
##
## @item @qcode{'Standardize'} @tab Whether to standardize @var{X} to zero mean
## and unit variance before fitting (default @qcode{true}).  Coefficients are
## always returned on the original scale.
##
## @item @qcode{'Weights'} @tab A vector of non-negative observation weights.
##
## @item @qcode{'RelTol'} @tab Convergence tolerance for the coordinate descent
## (default @math{1e-4}).
##
## @item @qcode{'MaxIter'} @tab Maximum number of coordinate-descent iterations
## (default @math{1e5}).
##
## @item @qcode{'DFmax'} @tab The maximum number of non-zero coefficients; the
## default sequence stops once this is exceeded.
##
## @item @qcode{'Intercept'} @tab Whether to fit a constant term (default
## @qcode{true}).
##
## @item @qcode{'PredictorNames'} @tab A cell array of predictor names, kept in
## @var{FitInfo}.
##
## @item @qcode{'CV'} @tab The number of folds @math{K} for @math{K}-fold
## cross-validation of the mean squared error, or a @code{cvpartition} object.
##
## @item @qcode{'MCReps'} @tab The number of Monte-Carlo repetitions of the
## cross-validation (default @math{1}).
## @end multitable
##
## When @qcode{'CV'} is used, @var{FitInfo}@code{.MSE} is the cross-validated
## error, plus @code{SE}, @code{LambdaMinMSE}, @code{IndexMinMSE},
## @code{Lambda1SE}, and @code{Index1SE}, which report the @var{Lambda} with the
## lowest error and the largest @var{Lambda} within one standard error of it.
## The fold assignment is random, so these selections are not reproducible
## without fixing the random seed.
##
## @seealso{ridge, regress, lassoglm}
## @end deftypefn

function [B, FitInfo] = lasso (X, y, varargin)

  if (nargin < 2)
    print_usage ();
  endif
  if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
    error ("lasso: X must be a real matrix.");
  endif
  if (! (isnumeric (y) && isreal (y) && isvector (y)))
    error ("lasso: Y must be a real vector.");
  endif
  y = y(:);
  if (rows (X) != numel (y))
    error ("lasso: X and Y must have the same number of observations.");
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
  cvarg = [];
  mcreps = 1;
  if (mod (numel (varargin), 2) != 0)
    error ("lasso: optional arguments must be Name-Value pairs.");
  endif
  for k = 1:2:numel (varargin)
    if (! ischar (varargin{k}))
      error ("lasso: parameter names must be character vectors.");
    endif
    switch (lower (varargin{k}))
      case "alpha";          alpha = varargin{k+1};
      case "lambda";         lambda = varargin{k+1};
      case "numlambda";      numlambda = varargin{k+1};
      case "lambdaratio";    lambdaratio = varargin{k+1};
      case "standardize";    standardize = logical (varargin{k+1});
      case "weights";        weights = varargin{k+1};
      case "reltol";         reltol = varargin{k+1};
      case "maxiter";        maxiter = varargin{k+1};
      case "dfmax";          dfmax = varargin{k+1};
      case "intercept";      intercept = logical (varargin{k+1});
      case "predictornames"; prednames = varargin{k+1};
      case "cv";             cvarg = varargin{k+1};
      case "mcreps";         mcreps = varargin{k+1};
      otherwise
        error ("lasso: unknown parameter name '%s'.", varargin{k});
    endswitch
  endfor
  if (! (isnumeric (alpha) && isscalar (alpha) && alpha > 0 && alpha <= 1))
    error ("lasso: 'Alpha' must be a scalar in (0, 1].");
  endif
  if (! isempty (lambda) && (! isnumeric (lambda) || any (lambda < 0)))
    error ("lasso: 'Lambda' must be a vector of non-negative values.");
  endif

  ## Drop observations with missing values.
  ok = ! (any (isnan (X), 2) | isnan (y));
  X = X(ok,:);  y = y(ok);
  [n, p] = size (X);

  ## Observation weights, normalised to sum to one.
  if (isempty (weights))
    w = ones (n, 1) / n;
  else
    if (! (isnumeric (weights) && isvector (weights) && numel (weights) == n
           && all (weights >= 0) && any (weights > 0)))
      error ("lasso: 'Weights' must be a non-negative vector, one per observation.");
    endif
    w = weights(:) / sum (weights);
  endif

  ## Centre and (optionally) standardise the predictors; centre the response.
  if (intercept)
    mux = w' * X;
    muy = w' * y;
  else
    mux = zeros (1, p);
    muy = 0;
  endif
  Xc = X - mux;
  if (standardize)
    sig = sqrt (w' * (Xc .^ 2));
    sig(sig == 0) = 1;
  else
    sig = ones (1, p);
  endif
  Xs = Xc ./ sig;
  yc = y - muy;

  ## Column "norms" d_j and the response correlations for the path bounds.
  d = (w' * (Xs .^ 2))';
  xy = Xs' * (w .* yc);
  lambda_max = max (abs (xy)) / alpha;

  ## Regularization path.
  userlambda = ! isempty (lambda);
  if (userlambda)
    lampath = sort (lambda(:)', "descend");   # fit large -> small (warm start)
  else
    if (isempty (lambdaratio))
      lambdaratio = ifelse (n < p, 1e-2, 1e-4);
    endif
    lampath = lambda_max * lambdaratio .^ ((0:numlambda-1) / (numlambda - 1));
  endif

  nullmse = w' * (yc .^ 2);
  soft = @(z, g) sign (z) .* max (abs (z) - g, 0);

  ## Fit along the path, warm-started from the previous solution.
  L = numel (lampath);
  Bstd = zeros (p, L);
  bs = zeros (p, 1);
  nkeep = L;
  for kk = 1:L
    lam = lampath(kk);
    for it = 1:maxiter
      bprev = bs;
      for j = 1:p
        rj = (w .* (yc - Xs * bs + Xs(:,j) * bs(j)));
        rho = Xs(:,j)' * rj;
        bs(j) = soft (rho, lam * alpha) / (d(j) + lam * (1 - alpha));
      endfor
      if (max (abs (bs - bprev)) <= reltol * max (max (abs (bs)), 1))
        break;
      endif
    endfor
    Bstd(:,kk) = bs;
    ## Trim the default path once the fit is essentially exact or too dense.
    if (! userlambda)
      bo = bs ./ sig';
      mse = w' * ((yc - Xs * bs) .^ 2);
      if (sum (bs != 0) > dfmax)
        nkeep = kk - 1;
        break;
      elseif (mse < 1e-3 * nullmse)
        nkeep = kk;
        break;
      endif
    endif
  endfor
  Bstd = Bstd(:,1:nkeep);
  lampath = lampath(1:nkeep);

  ## Unwind to the original scale; assemble outputs in ascending Lambda order.
  B = flip (Bstd ./ sig', 2);
  lampath = flip (lampath);
  icept = zeros (1, columns (B));
  mse = zeros (1, columns (B));
  for kk = 1:columns (B)
    icept(kk) = muy - mux * B(:,kk);
    mse(kk) = w' * ((y - icept(kk) - X * B(:,kk)) .^ 2);
  endfor

  FitInfo.Intercept = icept;
  FitInfo.Lambda = lampath;
  FitInfo.Alpha = alpha;
  FitInfo.DF = sum (B != 0, 1);
  FitInfo.MSE = mse;
  if (! isempty (prednames))
    FitInfo.PredictorNames = prednames;
  endif

  ## Cross-validation of the mean squared error along the fitted path.  The
  ## fold assignment is random, so the selected lambdas are not reproducible
  ## across runs (or identical to MATLAB) without a fixed random seed.
  if (! isempty (cvarg))
    lam = FitInfo.Lambda;
    foldmse = [];
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
        wtr = w(tr) / sum (w(tr));
        yhat = lasso_foldpredict_ (X(tr,:), y(tr), wtr, X(te,:), lam, alpha, ...
                                   standardize, intercept, reltol, maxiter);
        foldmse(end+1,:) = mean ((y(te) - yhat) .^ 2, 1);
      endfor
    endfor
    cvmse = mean (foldmse, 1);
    cvse = std (foldmse, 0, 1) / sqrt (rows (foldmse));
    FitInfo.MSE = cvmse;
    FitInfo.SE = cvse;
    [~, imin] = min (cvmse);
    FitInfo.LambdaMinMSE = lam(imin);
    FitInfo.IndexMinMSE = imin;
    i1se = max (find (cvmse <= cvmse(imin) + cvse(imin)));
    FitInfo.Lambda1SE = lam(i1se);
    FitInfo.Index1SE = i1se;
  endif

endfunction

## Fit the lasso path on training data and predict at the test rows, for the
## fixed sequence LAM (returns an ntest-by-numel(LAM) matrix of predictions).
function yhat = lasso_foldpredict_ (Xtr, ytr, w, Xte, lam, alpha, ...
                                    standardize, intercept, reltol, maxiter)
  [n, p] = size (Xtr);
  if (intercept)
    mux = w' * Xtr;
    muy = w' * ytr;
  else
    mux = zeros (1, p);
    muy = 0;
  endif
  Xc = Xtr - mux;
  if (standardize)
    sig = sqrt (w' * (Xc .^ 2));
    sig(sig == 0) = 1;
  else
    sig = ones (1, p);
  endif
  Xs = Xc ./ sig;
  yc = ytr - muy;
  d = (w' * (Xs .^ 2))';
  soft = @(z, g) sign (z) .* max (abs (z) - g, 0);
  [lamd, ord] = sort (lam(:)', "descend");     # warm-start from large lambda
  L = numel (lamd);
  Bstd = zeros (p, L);
  bs = zeros (p, 1);
  for kk = 1:L
    for it = 1:maxiter
      bprev = bs;
      for j = 1:p
        rj = w .* (yc - Xs * bs + Xs(:,j) * bs(j));
        bs(j) = soft (Xs(:,j)' * rj, lamd(kk) * alpha) ...
                / (d(j) + lamd(kk) * (1 - alpha));
      endfor
      if (max (abs (bs - bprev)) <= reltol * max (max (abs (bs)), 1))
        break;
      endif
    endfor
    Bstd(:,kk) = bs;
  endfor
  Bo = Bstd ./ sig';
  yhat = zeros (rows (Xte), L);
  yhat(:,ord) = (muy - mux * Bo) + Xte * Bo;
endfunction

function out = ifelse (cond, a, b)
  if (cond)
    out = a;
  else
    out = b;
  endif
endfunction

%!demo
%! ## Lasso path drives coefficients to zero as the penalty grows
%! rand ("seed", 1);
%! X = rand (50, 6);
%! b = [3; 0; -2; 0; 1.5; 0];
%! y = X * b + 0.1 * randn (50, 1);
%! [B, FitInfo] = lasso (X, y);
%! plot (log (FitInfo.Lambda), B');
%! xlabel ("log (Lambda)");  ylabel ("coefficient");

%!shared X, y
%! n = 20;
%! X = zeros (n, 6);
%! for j = 1:6
%!   X(:,j) = mod ((1:n)' * j, 7) + cos ((1:n)' * j);
%! endfor
%! y = X * [3;0;-2;0;1.5;0] + 0.1 * sin ((1:n)' * 3);

%!test  # MATLAB parity: coefficients and intercept at explicit Lambda (lasso)
%! lam = [2 1 0.5 0.2 0.1 0.05 0.01];
%! [B, I] = lasso (X, y, "Lambda", lam);
%! assert_equal (B(:,1), [2.99553372017526; 0; -1.99263231698882; ...
%!                  0.00211715208971491; 1.49324410153201; 0], 1e-3);
%! assert_equal (B(:,7), [2.23897682145602; 0; -0.482753193552206; ...
%!                  0.0267252048528861; 0.118128986514695; 0], 1e-3);
%! assert_equal (I.Intercept(1), 0.00566677333597632, 1e-3);
%! assert_equal (I.Intercept(7), 2.00069783707236, 1e-3);
%! assert_equal (I.DF, [4 4 4 4 4 4 4]);
%! assert_equal (I.MSE(1), 0.00590963740808093, 1e-3);

%!test  # MATLAB parity: Lambda is ascending and columns follow it
%! [B, I] = lasso (X, y, "Lambda", [2 1 0.5 0.2 0.1 0.05 0.01]);
%! assert_equal (I.Lambda, [0.01 0.05 0.1 0.2 0.5 1 2], 1e-12);
%! assert_equal (issorted (I.Lambda), true);

%!test  # MATLAB parity: default path endpoints and length
%! [B, I] = lasso (X, y);
%! assert_equal (I.Lambda(end), 7.18535290566157, 1e-4);
%! assert_equal (I.Lambda(1), 0.119858910419061, 1e-4);
%! assert_equal (numel (I.Lambda), 45);

%!test  # MATLAB parity: elastic net (Alpha = 0.5)
%! B = lasso (X, y, "Lambda", [2 1 0.5 0.2 0.1 0.05 0.01], "Alpha", 0.5);
%! assert_equal (B(:,1), [2.945477; -0.009332; -1.940313; ...
%!                  0.053970; 1.485977; -0.049272], 2e-3);

%!test  # MATLAB parity: Standardize = false
%! [B, I] = lasso (X, y, "Lambda", [2 1 0.5 0.2 0.1 0.05 0.01], ...
%!                 "Standardize", false);
%! assert_equal (B(:,1), [2.99774227822918; 0; -1.99668154187853; ...
%!                  0.00196231656035102; 1.49691209089871; 0], 2e-3);

%!test  # at Lambda -> 0 the lasso solution approaches least squares
%! B = lasso (X, y, "Lambda", [1 0.1 0.001]);   # columns ascend in Lambda
%! bols = regress (y - mean (y), (X - mean (X)));
%! assert_equal (B(:,1), bols, 1e-2);                 # smallest Lambda ~ OLS

%!test  # a large Lambda drives all coefficients to exactly zero
%! B = lasso (X, y, "Lambda", 100);
%! assert_equal (B, zeros (6, 1));

%!test  # DFmax limits the number of non-zero coefficients on the path
%! [B, I] = lasso (X, y, "DFmax", 2);
%! assert_equal (all (I.DF <= 2), true);

%!test  # cross-validation adds the selection fields and they are consistent
%! rand ("seed", 7);
%! [B, I] = lasso (X, y, "CV", 5);
%! assert_equal (numel (I.MSE), numel (I.Lambda));
%! assert_equal (numel (I.SE), numel (I.Lambda));
%! assert_equal (all (I.SE >= 0), true);
%! assert_equal (I.LambdaMinMSE, I.Lambda(I.IndexMinMSE), 1e-12);
%! assert_equal (I.Lambda1SE, I.Lambda(I.Index1SE), 1e-12);
%! assert_equal (I.Index1SE >= I.IndexMinMSE, true);           # 1-SE picks a larger Lambda
%! assert_equal (I.MSE(I.IndexMinMSE), min (I.MSE), 1e-12);

%!test  # cross-validated MSE at the min is within one SE at the 1-SE lambda
%! rand ("seed", 3);
%! [~, I] = lasso (X, y, "CV", 4, "MCReps", 2);
%! thr = I.MSE(I.IndexMinMSE) + I.SE(I.IndexMinMSE);
%! assert_equal (I.MSE(I.Index1SE) <= thr + 1e-9, true);

## Test input validation
%!error <Invalid call to lasso> lasso (1)
%!error <lasso: X and Y must have the same number of observations.> ...
%! lasso ([1 2; 3 4], [1 2 3]')
%!error <lasso: 'Alpha' must be a scalar in .0, 1..> lasso (X, y, "Alpha", 0)
%!error <lasso: 'Lambda' must be a vector of non-negative values.> ...
%! lasso (X, y, "Lambda", [-1 2])
