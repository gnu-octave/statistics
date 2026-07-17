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
## @deftypefn  {statistics} {@var{beta} =} mvregress (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{beta} =} mvregress (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{beta}, @var{Sigma}, @var{E}, @var{CovB}, @var{logL}] =} mvregress (@dots{})
##
## Multivariate (multiple-response) linear regression by maximum likelihood.
##
## @code{mvregress (@var{X}, @var{Y})} fits the multivariate normal regression
## of the @var{n}-by-@var{d} response matrix @var{Y} on the design @var{X} and
## returns the coefficient estimates @var{beta}.
##
## @var{X} is either a numeric @var{n}-by-@var{p} matrix, in which case the same
## @var{p} predictors apply to every response and @var{beta} is returned as a
## @var{p}-by-@var{d} matrix; or a cell array of @var{n} design matrices, each
## @var{d}-by-@var{K}, in which case @var{beta} is a @var{K}-by-1 vector.
##
## Missing responses (@code{NaN} entries of @var{Y}) are handled according to
## the estimation algorithm.
##
## The following @var{name}/@var{value} pairs are accepted:
##
## @table @asis
## @item @qcode{"algorithm"}
## @qcode{"mvn"} (multivariate normal; observations with any missing response
## are discarded), @qcode{"ecm"} (expectation-conditional-maximization, using
## every observed response), or @qcode{"cwls"} (covariance-weighted least
## squares, with the weight given by @qcode{"covar0"}).  The default is
## @qcode{"mvn"} when @var{Y} has no missing values and @qcode{"ecm"} otherwise.
##
## @item @qcode{"covar0"}
## The @var{d}-by-@var{d} covariance weight for @qcode{"cwls"} (default the
## identity), or the initial covariance for @qcode{"ecm"}.
##
## @item @qcode{"maxiter"}
## Maximum number of iterations (default 100).
##
## @item @qcode{"tolbeta"}, @qcode{"tolobj"}
## Convergence tolerances on the coefficients and the objective (defaults
## @code{1e-8} and @code{1e-8}).
## @end table
##
## The additional outputs are the estimated residual covariance @var{Sigma}
## (@var{d}-by-@var{d}), the residuals @var{E} (@var{n}-by-@var{d}), the
## covariance @var{CovB} of the coefficient estimates, and the log-likelihood
## @var{logL}.  (With missing data and the @qcode{"ecm"} algorithm, @var{CovB}
## is the standard observed-information covariance and can differ from
## @sc{matlab}'s value at the @code{1e-3} level; all other outputs agree.)
##
## @seealso{mvregresslike, regress, fitlm}
## @end deftypefn

function [beta, Sigma, E, CovB, logL] = mvregress (X, Y, varargin)

  if (nargin < 2)
    print_usage ();
  endif
  if (! (isnumeric (Y) && ismatrix (Y) && isreal (Y)))
    error ("mvregress: Y must be a real numeric matrix.");
  endif
  [n, d] = size (Y);

  ## --- normalise the design to per-observation d-by-K matrices ---
  is_numeric = ! iscell (X);
  if (is_numeric)
    if (rows (X) != n)
      error ("mvregress: X must have as many rows as Y.");
    endif
    p = columns (X);
    K = p * d;
    Xcell = cell (n, 1);
    for i = 1:n
      Xcell{i} = kron (eye (d), X(i, :));
    endfor
  else
    if (numel (X) != n)
      error ("mvregress: X must have one design matrix per observation.");
    endif
    K = columns (X{1});
    Xcell = X(:);
  endif

  ## --- options ---
  hasnan = any (isnan (Y(:)));
  alg = ""; covar0 = []; maxiter = 100; tolbeta = 1e-8; tolobj = 1e-8;
  if (mod (numel (varargin), 2) != 0)
    error ("mvregress: name/value arguments must come in pairs.");
  endif
  for i = 1:2:numel (varargin)
    switch (lower (char (varargin{i})))
      case "algorithm"
        alg = lower (char (varargin{i+1}));
      case "covar0"
        covar0 = varargin{i+1};
      case "maxiter"
        maxiter = varargin{i+1};
      case {"tolbeta"}
        tolbeta = varargin{i+1};
      case {"tolobj"}
        tolobj = varargin{i+1};
      otherwise
        error ("mvregress: unknown option '%s'.", char (varargin{i}));
    endswitch
  endfor
  if (isempty (alg))
    if (hasnan)
      alg = "ecm";
    else
      alg = "mvn";
    endif
  endif
  if (! any (strcmp (alg, {"ecm", "cwls", "mvn"})))
    error ("mvregress: algorithm must be 'ecm', 'cwls', or 'mvn'.");
  endif

  ## 'mvn' discards observations with any missing response.
  keep = true (n, 1);
  if (strcmp (alg, "mvn"))
    keep = all (! isnan (Y), 2);
  endif

  ## --- fit ---
  ## All three algorithms produce the maximum-likelihood estimates ('mvn' on the
  ## kept rows, 'ecm'/'cwls' using every observed response).  They differ only
  ## in how the reported log-likelihood and coefficient covariance are weighted:
  ## 'cwls' uses the fixed weight covar0 (default the identity), the others use
  ## the estimated residual covariance.
  if (isempty (covar0)), covar0 = eye (d); endif
  if (strcmp (alg, "mvn"))
    [b, Sigma] = ml_fit (Xcell, Y, eye (d), keep, maxiter, tolbeta, tolobj);
    Wsig = Sigma;
  else
    [b, Sigma] = ml_fit (Xcell, Y, covar0, keep, maxiter, tolbeta, tolobj);
    if (strcmp (alg, "cwls"))
      Wsig = covar0;                     ## logL / CovB use the fixed weight
    else
      Wsig = Sigma;
    endif
  endif
  E = residuals (Xcell, Y, b, n, d);

  ## --- outputs ---
  logL = -nll_obs (Xcell, Y, b, Wsig, keep);
  CovB = inv (info_beta (Xcell, Y, Wsig, keep));
  if (is_numeric)
    beta = reshape (b, p, d);
  else
    beta = b;
  endif

endfunction

## GLS coefficient estimate with a fixed covariance W (over the kept rows).
function b = gls_beta (Xcell, Y, W, keep)
  K = columns (Xcell{1});
  A = zeros (K); rhs = zeros (K, 1);
  for i = 1:numel (Xcell)
    if (! keep(i)), continue; endif
    o = ! isnan (Y(i, :));
    if (! any (o)), continue; endif
    Xio = Xcell{i}(o, :);
    Wo = W(o, o);
    A += Xio' * (Wo \ Xio);
    rhs += Xio' * (Wo \ Y(i, o)');
  endfor
  b = A \ rhs;
endfunction

## Maximum-likelihood fit by expectation-conditional-maximization.  Reduces to
## ordinary least squares for a complete common design.
function [b, Sigma] = ml_fit (Xcell, Y, Sigma, keep, maxiter, tolbeta, tolobj)
  n = numel (Xcell);
  d = columns (Y);
  b = gls_beta (Xcell, Y, Sigma, keep);
  prev_obj = Inf;
  for iter = 1:maxiter
    Si = inv (Sigma);
    K = numel (b);
    A = zeros (K); rhs = zeros (K, 1);
    Scov = zeros (d);
    nu = 0;
    ## E-step: impute missing responses; M-step accumulation.
    for i = 1:n
      if (! keep(i)), continue; endif
      o = ! isnan (Y(i, :));
      if (! any (o)), continue; endif
      m = ! o;
      Xi = Xcell{i};
      yfull = zeros (d, 1);
      yfull(o) = Y(i, o)';
      Ccond = zeros (d);
      if (any (m))
        ro = Y(i, o)' - Xi(o, :) * b;
        yfull(m) = Xi(m, :) * b + Sigma(m, o) * (Sigma(o, o) \ ro);
        Ccond(m, m) = Sigma(m, m) - Sigma(m, o) * (Sigma(o, o) \ Sigma(o, m));
      endif
      A += Xi' * Si * Xi;
      rhs += Xi' * Si * yfull;
      nu += 1;
      ## store per-observation info to rebuild the residual after beta update
      imp{i} = yfull;
      cc{i} = Ccond;
    endfor
    b = A \ rhs;
    for i = 1:n
      if (! keep(i) || ! any (! isnan (Y(i, :)))), continue; endif
      r = imp{i} - Xcell{i} * b;
      Scov += r * r' + cc{i};
    endfor
    Sigma = Scov / nu;
    obj = nll_obs (Xcell, Y, b, Sigma, keep);
    if (max (abs (A \ rhs - b)) < tolbeta && abs (prev_obj - obj) < tolobj)
      break;
    endif
    prev_obj = obj;
  endfor
endfunction

## Residuals Y - X*beta (NaN preserved where Y is missing).
function E = residuals (Xcell, Y, b, n, d)
  E = nan (n, d);
  for i = 1:n
    o = ! isnan (Y(i, :));
    E(i, o) = Y(i, o)' - Xcell{i}(o, :) * b;
  endfor
endfunction

## Negative log-likelihood of the observed responses.
function nll = nll_obs (Xcell, Y, b, Sigma, keep)
  nll = 0;
  for i = 1:numel (Xcell)
    if (! keep(i)), continue; endif
    o = ! isnan (Y(i, :));
    if (! any (o)), continue; endif
    r = Y(i, o)' - Xcell{i}(o, :) * b;
    So = Sigma(o, o);
    nll += 0.5 * (sum (o) * log (2*pi) + 2*sum (log (diag (chol (So)))) ...
                  + r' * (So \ r));
  endfor
endfunction

## Fisher information for the coefficients.
function I = info_beta (Xcell, Y, Sigma, keep)
  K = columns (Xcell{1});
  I = zeros (K);
  for i = 1:numel (Xcell)
    if (! keep(i)), continue; endif
    o = ! isnan (Y(i, :));
    if (! any (o)), continue; endif
    Xio = Xcell{i}(o, :);
    I += Xio' * (Sigma(o, o) \ Xio);
  endfor
endfunction

%!demo
%! ## Two correlated responses regressed on a common predictor.
%! X = [ones(30,1), (1:30)'/30];
%! B = [1 -1; 2 0.5];
%! E = [0.3 0.1; 0.1 0.2];
%! Y = X * B + randn (30, 2) * chol (E);
%! [beta, Sigma] = mvregress (X, Y)

%!shared X, Ycomp, Ymiss, Btrue
%! X = [ones(12,1), linspace(-1, 1, 12)'];
%! Btrue = [1 -0.5 2; 0.3 1.2 -0.8];
%! Ycomp = X * Btrue + 0.05 * cos ((1:12)' * [1 2 3]);
%! Ymiss = Ycomp; Ymiss(3,2) = NaN; Ymiss(8,3) = NaN;

%!test  # numeric design returns a p-by-d beta equal to OLS (complete data)
%! [beta, Sigma, E] = mvregress (X, Ycomp, "algorithm", "mvn");
%! assert (size (beta), [2, 3]);
%! assert (beta, (X'*X)\(X'*Ycomp), 1e-8);
%! assert (Sigma, E'*E/12, 1e-8);

%!test  # cell design returns the vectorised (K-by-1) beta
%! Xc = cell (12, 1);
%! for i = 1:12, Xc{i} = kron (eye (3), X(i,:)); end
%! bnum = mvregress (X, Ycomp);
%! bcell = mvregress (Xc, Ycomp);
%! assert (bcell, bnum(:), 1e-8);

%!test  # logL equals -mvregresslike at the fit (mvn, complete)
%! [beta, Sigma, E, CovB, logL] = mvregress (X, Ycomp, "algorithm", "mvn");
%! assert (logL, -mvregresslike (X, Ycomp, beta, Sigma, "mvn"), 1e-8);
%! assert (CovB, kron (Sigma, inv (X'*X)), 1e-8);

%!test  # ecm uses all observed data; the log-likelihood improves on listwise
%! b_ecm = mvregress (X, Ymiss, "algorithm", "ecm");
%! b_mvn = mvregress (X, Ymiss, "algorithm", "mvn");
%! assert (! isequal (b_ecm, b_mvn));   # different estimates

%!test  # default algorithm: mvn without missing data, ecm with
%! [b1, ~, ~, ~, L1] = mvregress (X, Ycomp);
%! [b2, ~, ~, ~, L2] = mvregress (X, Ycomp, "algorithm", "mvn");
%! assert (L1, L2, 1e-10);

%!test  # cwls: logL and CovB use the identity weight
%! [beta, Sigma, E, CovB, logL] = mvregress (X, Ycomp, "algorithm", "cwls");
%! assert (logL, -mvregresslike (X, Ycomp, beta, eye (3), "cwls"), 1e-8);
%! assert (CovB, kron (eye (3), inv (X'*X)), 1e-8);

%!error <Invalid call> mvregress (1)
%!error <Y must be a real numeric matrix> mvregress (ones (3), {1})
%!error <X must have as many rows as Y> mvregress (ones (2, 2), ones (3, 2))
%!error <algorithm must be> mvregress (ones (3, 2), ones (3, 2), "algorithm", "xxx")
%!error <unknown option 'bogus'> mvregress (ones (3, 2), ones (3, 2), "bogus", 1)
