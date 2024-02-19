## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{nlogL} =} nakalike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} nakalike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} nakalike (@var{params}, @var{x}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} nakalike (@var{params}, @var{x}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the Nakagami distribution.
##
## @code{@var{nlogL} = nakalike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the Nakagami
## distribution with (1) scale parameter @var{mu} and (2) shape parameter
## @var{omega} given in the two-element vector @var{params}.
##
## @code{[@var{nlogL}, @var{acov}] = nakalike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{params} are their asymptotic variances.
##
## @code{[@dots{}] = nakalike (@var{params}, @var{x}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = nakalike (@var{params}, @var{x}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} must contain non-negative integer frequencies for the
## corresponding elements in @var{x}.  By default, or if left empty,
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the Nakagami distribution can be found at
## @url{https://en.wikipedia.org/wiki/Nakagami_distribution}
##
## @seealso{nakacdf, nakainv, nakapdf, nakarnd, nakafit}
## @end deftypefn

function [nlogL, acov] = nakalike (params, x, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("nakalike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("nakalike: X must be a vector.");
  endif

  if (length (params) != 2)
    error ("nakalike: PARAMS must be a two-element vector.");
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("nakalike: X and CENSOR vector mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("nakalike: X and FREQ vector mismatch.");
  elseif (any (freq < 0))
    error ("nakalike: FREQ must not contain negative values.");
  elseif (any (fix (freq) != freq))
    error ("nakafit: FREQ must contain integer values.");
  endif

  ## Expand frequency and censor vectors (if necessary)
  if (! all (freq == 1))
    xf = [];
    cf = [];
    for i = 1:numel (freq)
      xf = [xf, repmat(x(i), 1, freq(i))];
      cf = [cf, repmat(censor(i), 1, freq(i))];
    endfor
    x = xf;
    freq = ones (size (x));
    censor = cf;
  endif

  ## Get parameters
  mu = params(1);
  omega = params(2);
  log_a = gammaln (mu);
  log_b = log (omega / mu);

  z = x .^ 2 ./ (omega / mu);
  log_z = log (z);
  L = (mu - 1) .* log_z - z - log_a - log_b + log (2 .* x);

  ## Handle censored data
  n_censored = sum (freq .* censor);
  if (n_censored > 0)
    censored = (censor == 1);
    z_censored = z(censored);
    [S, dS] = dgammainc (z_censored, mu);
    L(censored) = log (S);
  endif
  ## Sum up the neg log likelihood
  nlogL = -sum (freq .* L);


  ## Compute asymptotic covariance
  if (nargout > 1)
    ## Compute first order central differences of the log-likelihood gradient
    dp = 0.0001 .* max (abs (params), 1);

    ngrad_p1 = logl_grad (params + [dp(1), 0], x, censor, freq);
    ngrad_m1 = logl_grad (params - [dp(1), 0], x, censor, freq);
    ngrad_p2 = logl_grad (params + [0, dp(2)], x, censor, freq);
    ngrad_m2 = logl_grad (params - [0, dp(2)], x, censor, freq);

    ## Compute negative Hessian by normalizing the differences by the increment
    nH = [(ngrad_p1(:) - ngrad_m1(:))./(2 * dp(1)), ...
          (ngrad_p2(:) - ngrad_m2(:))./(2 * dp(2))];

    ## Force neg Hessian being symmetric
    nH = 0.5 .* (nH + nH');
    ## Check neg Hessian is positive definite
    [R, p] = chol (nH);
    if (p > 0)
      warning ("nakalike: non positive definite Hessian matrix.");
      acov = NaN (2);
      return
    endif
    ## ACOV estimate is the negative inverse of the Hessian.
    Rinv = inv (R);
    acov = Rinv * Rinv;
  endif

endfunction

## Helper function for computing negative gradient
function ngrad = logl_grad (params, x, censor, freq)
  mu = params(1);
  omega = params(2);
  ## Transform to Gamma parameters
  log_a = gammaln (mu);
  log_b = log (omega / mu);

  z = x .^ 2 ./ (omega / mu);
  log_z = log (z);
  dL1 = log_z - psi (mu);
  dL2 = (z - mu) ./ (omega / mu);

  ## Handle censored data
  n_censored = sum (freq .* censor);
  if (n_censored > 0)
    censored = (censor == 1);
    z_censored = z(censored);
    [S, dS] = dgammainc (z_censored, a);
    dL1(censored) = dS ./ S;
    tmp = mu .* log_z(censored) - log_b - z_censored - log_a;
    dL2(censored) = exp (tmp) ./ S;
  endif
  ngrad = -[sum(freq .* dL1), sum(freq .* dL2)];

  ## Transform back to Nakagami parameters
  ngrad = ngrad * [1, 0; (-omega ./ (mu .^ 2)), (1 ./ mu)];
endfunction

## Compute the incomplete Gamma function with its 1st and 2nd derivatives
function [y, dy, d2y] = dgammainc (x, k)

  ## Initialize return variables
  y = nan (size (x));
  dy = y;
  d2y = y;

  ## Use approximation for K > 2^20
  ulim = 2^20;
  is_lim = find (k > ulim);
  if (! isempty (is_lim))
    x(is_lim) = max (ulim - 1/3 + sqrt (ulim ./ k(is_lim)) .* ...
                     (x(is_lim) - (k(is_lim) - 1/3)), 0);
    k(is_lim) = ulim;
  endif

  ## For x < k+1
  is_lo = find(x < k + 1 & x != 0);
  if (! isempty (is_lo))
    x_lo = x(is_lo);
    k_lo = k(is_lo);
    k_1 = k_lo;
    step = 1;
    d1st = 0;
    d2st = 0;
    stsum = step;
    d1sum = d1st;
    d2sum = d2st;
    while norm (step, "inf") >= 100 * eps (norm (stsum, "inf"))
      k_1 += 1;
      step = step .* x_lo ./ k_1;
      d1st = (d1st .* x_lo - step) ./ k_1;
      d2st = (d2st .* x_lo - 2 .* d1st) ./ k_1;
      stsum = stsum + step;
      d1sum = d1sum + d1st;
      d2sum = d2sum + d2st;
    endwhile
    fklo = exp (-x_lo + k_lo .* log (x_lo) - gammaln (k_lo + 1));
    y_lo = fklo .* stsum;
    ## Fix very small k
    y_lo(x_lo > 0 & y_lo > 1) = 1;
    ## Compute 1st derivative
    dlogfklo = (log (x_lo) - psi (k_lo + 1));
    d1fklo = fklo .* dlogfklo;
    d1y_lo = d1fklo .* stsum + fklo .* d1sum;
    ## Compute 2nd derivative
    d2fklo = d1fklo .* dlogfklo - fklo .* psi (1, k_lo + 1);
    d2y_lo = d2fklo .* stsum + 2 .* d1fklo .* d1sum + fklo .* d2sum;
    ## Considering the upper tail
    y(is_lo) = 1 - y_lo;
    dy(is_lo) = -d1y_lo;
    d2y(is_lo) = -d2y_lo;
  endif

  ## For x >= k+1
  is_hi = find(x >= k+1);
  if (! isempty (is_hi))
    x_hi = x(is_hi);
    k_hi = k(is_hi);
    zc = 0;
    k0 = 0;
    k1 = k_hi;
    x0 = 1;
    x1 = x_hi;
    d1k0 = 0;
    d1k1 = 1;
    d1x0 = 0;
    d1x1 = 0;
    d2k0 = 0;
    d2k1 = 0;
    d2x0 = 0;
    d2x2 = 0;
    kx = k_hi ./ x_hi;
    d1kx = 1 ./ x_hi;
    d2kx = 0;
    start = 1;
    while norm (d2kx - start, "Inf") > 100 * eps (norm (d2kx, "Inf"))
      rescale = 1 ./ x1;
      zc += 1;
      n_k = zc - k_hi;
      d2k0 = (d2k1 + d2k0 .* n_k - 2 .* d1k0) .* rescale;
      d2x0 = (d2x2 + d2x0 .* n_k - 2 .* d1x0) .* rescale;
      d1k0 = (d1k1 + d1k0 .* n_k - k0) .* rescale;
      d1x0 = (d1x1 + d1x0 .* n_k - x0) .* rescale;
      k0 = (k1 + k0 .* n_k) .* rescale;
      x0 = 1 + (x0 .* n_k) .* rescale;
      nrescale = zc .* rescale;
      d2k1 = d2k0 .* x_hi + d2k1 .* nrescale;
      d2x2 = d2x0 .* x_hi + d2x2 .* nrescale;
      d1k1 = d1k0 .* x_hi + d1k1 .* nrescale;
      d1x1 = d1x0 .* x_hi + d1x1 .* nrescale;
      k1 = k0 .* x_hi + k1 .* nrescale;
      x1 = x0 .* x_hi + zc;
      start = d2kx;
      kx = k1 ./ x1;
      d1kx = (d1k1 - kx .* d1x1) ./ x1;
      d2kx = (d2k1 - d1kx .* d1x1 - kx .* d2x2 - d1kx .* d1x1) ./ x1;
    endwhile
    fkhi = exp (-x_hi + k_hi .* log (x_hi) - gammaln (k_hi + 1));
    y_hi = fkhi .* kx;
    ## Compute 1st derivative
    dlogfkhi = (log (x_hi) - psi (k_hi + 1));
    d1fkhi = fkhi .* dlogfkhi;
    d1y_hi = d1fkhi .* kx + fkhi .* d1kx;
    ## Compute 2nd derivative
    d2fkhi = d1fkhi .* dlogfkhi - fkhi .* psi (1, k_hi + 1);
    d2y_hi = d2fkhi .* kx + 2 .* d1fkhi .* d1kx + fkhi .* d2kx;
    ## Considering the upper tail
    y(is_hi) = y_hi;
    dy(is_hi) = d1y_hi;
    d2y(is_hi) = d2y_hi;
  endif

  ## Handle x == 0
  is_x0 = find (x == 0);
  if (! isempty (is_x0))
    ## Considering the upper tail
    y(is_x0) = 1;
    dy(is_x0) = 0;
    d2y(is_x0) = 0;
  endif

  ## Handle k == 0
  is_k0 = find (k == 0);
  if (! isempty (is_k0))
    is_k0x0 = find (k == 0 & x == 0);
    ## Considering the upper tail
    y(is_k0) = 0;
    dy(is_k0x0) = Inf;
    d2y(is_k0x0) = -Inf;
  endif
endfunction

## Test output
%!test
%! nlogL = nakalike ([0.735504, 858.5], [1:50]);
%! assert (nlogL, 202.8689, 1e-4);
%!test
%! nlogL = nakalike ([1.17404, 11], [1:5]);
%! assert (nlogL, 8.6976, 1e-4);
%!test
%! nlogL = nakalike ([1.17404, 11], [1:5], [], [1, 1, 1, 1, 1]);
%! assert (nlogL, 8.6976, 1e-4);
%!test
%! nlogL = nakalike ([1.17404, 11], [1:6], [], [1, 1, 1, 1, 1, 0]);
%! assert (nlogL, 8.6976, 1e-4);

## Test input validation
%!error<nakalike: function called with too few input arguments.> nakalike (3.25)
%!error<nakalike: X must be a vector.> nakalike ([5, 0.2], ones (2))
%!error<nakalike: PARAMS must be a two-element vector.> ...
%! nakalike ([1, 0.2, 3], [1, 3, 5, 7])
%!error<nakalike: X and CENSOR vector mismatch.> ...
%! nakalike ([1.5, 0.2], [1:5], [0, 0, 0])
%!error<nakalike: X and FREQ vector mismatch.> ...
%! nakalike ([1.5, 0.2], [1:5], [0, 0, 0, 0, 0], [1, 1, 1])
%!error<nakalike: X and FREQ vector mismatch.> ...
%! nakalike ([1.5, 0.2], [1:5], [], [1, 1, 1])
%!error<nakalike: FREQ must not contain negative values.> ...
%! nakalike ([1.5, 0.2], [1:5], [], [1, 1, 1, 1, -1])
