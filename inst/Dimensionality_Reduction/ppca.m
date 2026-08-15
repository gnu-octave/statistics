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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{coeff} =} ppca (@var{Y}, @var{K})
## @deftypefnx {statistics} {[@var{coeff}, @var{score}] =} ppca (@var{Y}, @var{K})
## @deftypefnx {statistics} {[@var{coeff}, @var{score}, @var{pcvar}] =} ppca (@var{Y}, @var{K})
## @deftypefnx {statistics} {[@var{coeff}, @var{score}, @var{pcvar}, @var{mu}] =} ppca (@var{Y}, @var{K})
## @deftypefnx {statistics} {[@var{coeff}, @var{score}, @var{pcvar}, @var{mu}, @var{v}] =} ppca (@var{Y}, @var{K})
## @deftypefnx {statistics} {[@var{coeff}, @var{score}, @var{pcvar}, @var{mu}, @var{v}, @var{S}] =} ppca (@var{Y}, @var{K})
## @deftypefnx {statistics} {[@dots{}] =} ppca (@dots{}, @var{Name}, @var{Value})
##
## Probabilistic principal component analysis.
##
## @code{@var{coeff} = ppca (@var{Y}, @var{K})} fits a probabilistic principal
## component analysis (PPCA) model with @var{K} components to the @math{N * P}
## data matrix @var{Y} (rows are observations, columns are variables) and returns
## the @math{P * @var{K}} matrix @var{coeff} of orthonormal principal component
## coefficients, ordered by decreasing component variance.  @var{Y} may contain
## @code{NaN} values marking missing observations; the model is fitted by an
## expectation-maximization algorithm that accounts for them.  @var{K} must be a
## positive integer smaller than @math{P}.
##
## @code{[@var{coeff}, @var{score}, @var{pcvar}, @var{mu}, @var{v}, @var{S}] =
## ppca (@dots{})} returns further outputs:
##
## @table @var
## @item score
## The @math{N * @var{K}} principal component scores (the data projected onto the
## components; missing entries are reconstructed from the model before
## projection).
##
## @item pcvar
## A @math{@var{K} * 1} vector of the principal component variances (the variance
## explained by each component).
##
## @item mu
## A @math{1 * P} vector of the estimated mean of @var{Y}.
##
## @item v
## The estimated residual (isotropic noise) variance.
##
## @item S
## A structure with the fitted model details: the loadings @qcode{W}, the
## expected scores @qcode{Xexp}, the reconstruction @qcode{Recon}, the number of
## iterations @qcode{NumIter}, and the root-mean-square residual @qcode{RMSResid}.
## @end table
##
## Name/Value pairs control the fit:
##
## @table @asis
## @item @qcode{'W0'}
## A @math{P * @var{K}} initial value for the loadings used by the
## expectation-maximization algorithm (missing-data case).
##
## @item @qcode{'Options'}
## A structure of algorithm options, as returned by @code{statset}, whose
## @qcode{MaxIter}, @qcode{TolFun}, and @qcode{TolX} fields set the maximum number
## of iterations and the convergence tolerances of the
## expectation-maximization algorithm.
## @end table
##
## When @var{Y} has no missing values the model is fitted directly from the
## eigendecomposition of its covariance matrix; @var{coeff}, @var{pcvar}, and
## @var{v} are then the principal component directions, the leading variances, and
## the mean of the trailing variances, respectively.
##
## @seealso{pca, pcacov, pcares, factoran, barttest}
## @end deftypefn

function [coeff, score, pcvar, mu, v, S] = ppca (Y, K, varargin)

  ## Input validation
  if (nargin < 2)
    print_usage ();
  endif
  if (! (isnumeric (Y) && ismatrix (Y) && ndims (Y) == 2))
    error ("ppca: Y must be a numeric matrix.");
  endif
  if (! isreal (Y))
    error ("ppca: Y must be real.");
  endif

  [n, p] = size (Y);

  if (! (isnumeric (K) && isscalar (K) && isreal (K) && K >= 1 ...
         && K == fix (K)))
    error ("ppca: K must be a positive integer.");
  endif
  if (K >= p)
    error ("ppca: K must be smaller than the number of variables.");
  endif

  ## Defaults and Name/Value / Options parsing
  opts = struct ("MaxIter", 1000, "TolFun", 1e-6, "TolX", 1e-6);
  W0 = [];
  args = varargin;
  if (mod (numel (args), 2) != 0)
    error ("ppca: Name/Value arguments must come in pairs.");
  endif
  for k = 1:2:numel (args)
    name = args{k};
    val = args{k+1};
    if (! (ischar (name) && isrow (name)))
      error ("ppca: parameter names must be character vectors.");
    endif
    switch (lower (name))
      case 'w0'
        W0 = val;
      case 'options'
        opts = merge_statset (opts, val);
      otherwise
        error ("ppca: unknown parameter name '%s'.", name);
    endswitch
  endfor

  if (! isempty (W0) && ! isequal (size (W0), [p, K]))
    error ("ppca: 'W0' must be a %d-by-%d matrix.", p, K);
  endif

  ## Fit the model
  if (any (isnan (Y(:))))
    [coeff, W, ex, mu, v, numiter] = ppca_em (Y, K, W0, opts);
  else
    [coeff, W, ex, mu, v, numiter] = ppca_full (Y, K);
  endif

  ## Sign convention: the element of largest magnitude in each column is positive
  for j = 1:K
    [~, im] = max (abs (coeff(:,j)));
    if (coeff(im,j) < 0)
      coeff(:,j) = -coeff(:,j);
    endif
  endfor

  ## Reconstruct missing entries from the model and project onto the components
  obs = ! isnan (Y);
  Yhat = Y;
  Recon = ex * W' + mu;
  Yhat(! obs) = Recon(! obs);
  score = (Yhat - mu) * coeff;
  pcvar = var (score, 0)'(:);

  if (nargout > 5)
    resid = Yhat - (score * coeff' + mu);
    S = struct ("W", coeff .* sqrt (max (pcvar - v, 0))', ...
                "Xexp", score, "Recon", score * coeff' + mu, ...
                "NumIter", numiter, ...
                "RMSResid", sqrt (mean (resid(obs) .^ 2)));
  endif

endfunction

## Direct fit for fully observed data via the covariance eigendecomposition.
function [coeff, W, ex, mu, v, numiter] = ppca_full (Y, K)

  [n, p] = size (Y);
  mu = mean (Y);
  C = cov (Y);
  [U, L] = eig ((C + C') / 2);
  [lambda, idx] = sort (diag (L), "descend");
  U = U(:,idx);
  coeff = U(:,1:K);
  v = mean (lambda(K+1:end));
  ## ML loadings and expected scores (posterior mean) for the S struct
  W = coeff .* sqrt (max (lambda(1:K) - v, 0))';
  M = W' * W + v * eye (K);
  ex = (Y - mu) * W / M;
  numiter = 0;

endfunction

## Expectation-maximization fit for data with missing (NaN) entries.
function [coeff, W, ex, mu, v, iter] = ppca_em (Y, K, W0, opts)

  [n, p] = size (Y);
  obs = ! isnan (Y);

  ## Initialize the mean, a filled data matrix, and the loadings
  mu = zeros (1, p);
  for j = 1:p
    mu(j) = mean (Y(obs(:,j), j));
  endfor
  Yf = Y;
  for j = 1:p
    Yf(! obs(:,j), j) = mu(j);
  endfor
  if (isempty (W0))
    [~, Sv, V] = svd (Yf - mu, 0);
    W = V(:,1:K) * Sv(1:K,1:K) / sqrt (n);
  else
    W = W0;
  endif
  v = 1;
  ex = zeros (n, K);
  Exx = zeros (K, K, n);

  for iter = 1:opts.MaxIter
    Wold = W;
    muold = mu;
    vold = v;

    ## E-step: expected scores and second moments per observation
    for i = 1:n
      o = obs(i,:);
      Wo = W(o,:);
      yo = Y(i,o) - mu(o);
      Minv = inv (v * eye (K) + Wo' * Wo);
      ex(i,:) = (Minv * Wo' * yo')';
      Exx(:,:,i) = v * Minv + ex(i,:)' * ex(i,:);
    endfor

    ## M-step: update the mean and loadings per variable, then the noise
    for j = 1:p
      r = obs(:,j);
      exj = ex(r,:);
      yj = Y(r,j);
      ZtZ = [sum(Exx(:,:,r), 3), sum(exj, 1)'; sum(exj, 1), sum(r)];
      sol = ZtZ \ [exj' * yj; sum(yj)];
      W(j,:) = sol(1:K)';
      mu(j) = sol(K+1);
    endfor

    num = 0;
    cnt = 0;
    for i = 1:n
      o = obs(i,:);
      Wo = W(o,:);
      yo = Y(i,o) - mu(o);
      num += yo * yo' - 2 * ex(i,:) * Wo' * yo' ...
             + trace (Wo' * Wo * Exx(:,:,i));
      cnt += sum (o);
    endfor
    v = num / cnt;

    if (max (abs (W(:) - Wold(:))) < opts.TolX ...
        && max (abs (mu - muold)) < opts.TolX && abs (v - vold) < opts.TolX)
      break;
    endif
  endfor

  ## Orthonormalize the loadings into principal component directions
  [Uw, ~, ~] = svd (W, 0);
  coeff = Uw(:,1:K);

endfunction

## Copy the recognised statset fields from S into OPTS.
function opts = merge_statset (opts, s)

  if (! isstruct (s))
    error ("ppca: 'Options' must be a structure.");
  endif
  fn = fieldnames (s);
  for k = 1:numel (fn)
    val = s.(fn{k});
    if (isempty (val))
      continue;
    endif
    switch (lower (fn{k}))
      case 'maxiter'
        opts.MaxIter = val;
      case 'tolfun'
        opts.TolFun = val;
      case 'tolx'
        opts.TolX = val;
    endswitch
  endfor

endfunction

%!demo
%! ## Fit a two-component PPCA model and reconstruct the data.
%! Y = [ 1.0,  2.0,  0.5;  2.1,  3.9,  1.2; ...
%!      -1.0, -2.2, -0.4; -2.0, -3.8, -1.1; ...
%!       0.5,  1.1,  0.9;  1.6,  2.8, -0.2];
%! [coeff, score, pcvar, mu, v] = ppca (Y, 2);
%! coeff
%! pcvar
%! ## The scores reconstruct the data through the coefficients.
%! max (abs (vec (score * coeff' + mu - Y)))

## Reference values below are from MATLAB R2023b for
## Y = reshape (mod ((1:40)*7, 17), 10, 4) - 8.

%!test
%! Y = reshape (mod ((1:40)*7, 17), 10, 4) - 8;
%! opt = struct ("MaxIter", 2000, "TolFun", 1e-12, "TolX", 1e-12);
%! [coeff, score, pcvar, mu, v] = ppca (Y, 2, "Options", opt);
%! coeff_ref = [ 0.543716564044483, -0.391810693220729; ...
%!               0.599330273587830,  0.314784124203745; ...
%!               0.302509942024958,  0.761619894252268; ...
%!              -0.503649934101907,  0.409060475365608];
%! assert_equal (coeff, coeff_ref, 1e-6);
%! assert_equal (pcvar, [49.558237292069599; 30.960532427249241], 1e-4);
%! assert_equal (mu, [-0.1, 0.2, 0.5, -0.9], 1e-10);
%! assert_equal (v, 10.623960883864141, 1e-4);

%!test
%! Y = reshape (mod ((1:40)*7, 17), 10, 4) - 8;
%! [coeff, score] = ppca (Y, 2);
%! score_ref = [ -2.225140840148114,   4.921963749319239; ...
%!                7.787588722827117,  -7.324026640834169; ...
%!               -5.050861878994857,   1.641002156881092; ...
%!               10.104537612540948,   2.342550721115434; ...
%!               -7.876582917841600,  -1.639959435557055; ...
%!               -1.283233827199200,   6.015617613465288; ...
%!               -1.459120725735530, -11.581703199267448; ...
%!               -4.108954866045941,   2.734656021027141; ...
%!               11.046444625489862,   3.436204585261484; ...
%!               -6.934675904892686,  -0.546305571411006];
%! assert_equal (score, score_ref, 1e-4);

%!test
%! ## coeff is orthonormal and score is the projection of the centred data.
%! Y = reshape (mod ((1:40)*7, 17), 10, 4) - 8;
%! [coeff, score, pcvar, mu] = ppca (Y, 2);
%! assert_equal (coeff' * coeff, eye (2), 1e-10);
%! assert_equal (score, (Y - mu) * coeff, 1e-10);

%!test
%! ## Missing-data (NaN) case fitted by expectation-maximization.
%! Y = reshape (mod ((1:40)*7, 17), 10, 4) - 8;
%! Y(3,2) = NaN;
%! Y(7,4) = NaN;
%! opt = struct ("MaxIter", 5000, "TolFun", 1e-12, "TolX", 1e-12);
%! [coeff, score, pcvar, mu, v] = ppca (Y, 2, "Options", opt);
%! coeff_ref = [ 0.489878122099491, -0.455706110813571; ...
%!               0.631344835584810,  0.254329922477782; ...
%!               0.331227809703440,  0.763974559487245; ...
%!              -0.501708343709497,  0.379461596944783];
%! assert_equal (coeff, coeff_ref, 1e-4);
%! assert_equal (pcvar, [49.932943652722741; 28.917799526089862], 1e-3);
%! assert_equal (mu, [-0.099999999999924, 0.225205294055629, ...
%!                     0.500000000000008, -0.683793705409387], 1e-4);
%! assert_equal (v, 10.684669221332555, 1e-3);

%!test
%! ## Output sizes.
%! Y = reshape (mod ((1:40)*7, 17), 10, 4) - 8;
%! [coeff, score, pcvar, mu, v, S] = ppca (Y, 3);
%! assert_equal (size (coeff), [4, 3]);
%! assert_equal (size (score), [10, 3]);
%! assert_equal (size (pcvar), [3, 1]);
%! assert_equal (size (mu), [1, 4]);
%! assert_equal (isscalar (v), true);
%! assert_equal (isfield (S, "W") && isfield (S, "Recon"), true);

## Test input validation
%!error<Invalid call to ppca> ppca (ones (5, 3))
%!error<ppca: Y must be a numeric matrix.> ppca ({1, 2}, 1)
%!error<ppca: Y must be real.> ppca ([1+2i, 3; 4, 5], 1)
%!error<ppca: K must be a positive integer.> ppca (ones (5, 3), 0)
%!error<ppca: K must be a positive integer.> ppca (ones (5, 3), 1.5)
%!error<ppca: K must be smaller than the number of variables.> ppca (ones (5, 3), 3)
%!error<ppca: Name/Value arguments must come in pairs.> ppca (ones (5, 3), 1, "W0")
%!error<ppca: unknown parameter name 'foo'.> ppca (ones (5, 3), 1, "foo", 1)
%!error<ppca: 'W0' must be a 3-by-1 matrix.> ppca (ones (5, 3), 1, "W0", ones (2, 2))
%!error<ppca: 'Options' must be a structure.> ppca (ones (5, 3), 1, "Options", 5)
