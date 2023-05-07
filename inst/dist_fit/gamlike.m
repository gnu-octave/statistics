## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Based on previous work by Martijn van Oosterhout <kleptog@svana.org>
## originally granted to the public domain.
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
## @deftypefn  {statistics} {@var{nlogL} =} gamlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} gamlike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} gamlike (@var{params}, @var{x}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} gamlike (@var{params}, @var{x}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the Gamma distribution.
##
## @code{@var{nlogL} = gamlike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the Gamma distribution
## with (1) shape parameter @var{k} and (2) scale parameter @var{theta} given in
## the two-element vector @var{params}.
##
## @code{[@var{nlogL}, @var{acov}] = gamlike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{acov} are their asymptotic variances.
##
## @code{[@dots{}] = gamlike (@var{params}, @var{x}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = gamlike (@var{params}, @var{x}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## There are two equivalent parameterizations in common use:
## @enumerate
## @item With a shape parameter @math{k} and a scale parameter @math{θ}, which
## is used by @code{gamcdf}.
## @item With a shape parameter @math{α = k} and an inverse scale parameter
## @math{β = 1 / θ}, called a rate parameter.
## @end enumerate
##
## Further information about the Gamma distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gamma_distribution}
##
## @seealso{gamcdf, gampdf, gaminv, gamrnd, gamfit}
## @end deftypefn

function [nlogL, acov] = gamlike (params, x, censor, freq)

  ## Check input arguments and add defaults
  if (nargin < 2)
    error ("gamlike: function called with too few input arguments.");
  endif
  if (numel (params) != 2)
    error ("gamlike: wrong parameters length.");
  endif
  if (! isvector (x))
    error ("gamlike: X must be a vector.");
  endif
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("gamlike: X and CENSOR vectors mismatch.");
  endif
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (isequal (size (x), size (freq)))
    nulls = find (freq == 0);
    if (numel (nulls) > 0)
      x(nulls) = [];
      censor(nulls) = [];
      freq(nulls) = [];
    endif
  else
    error ("gamlike: X and FREQ vectors mismatch.");
  endif

  ## Get K and THETA values
  k = params(1);
  t = params(2);

  ## Parameters K and THETA must be positive, otherwise make them NaN
  k(k <= 0) = NaN;
  t(t <= 0) = NaN;

  ## Data in X must be positive, otherwise make it NaN
  x(x <= 0) = NaN;

  ## Compute the individual log-likelihood terms
  z = x ./ t;
  L = (k - 1) .* log (z) - z - gammaln (k) - log (t);
  n_censored = sum (freq .* censor);
  if (n_censored > 0)
    z_censored = z(logical (censor));
    Scen = gammainc (z_censored, k, "upper");
    L(logical (censor)) = log (Scen);
  endif

  ## Force a log(0)==-Inf for X from extreme right tail
  L(z == Inf) = -Inf;

  ## Neg-log-likelihood is the sum of the individual contributions
  nlogL = -sum (freq .* L);

  ## Compute the negative hessian at the parameter values.
  ## Invert to get the observed information matrix.
  if (nargout == 2)
    ## Calculate all data
    dL11 = -psi (1, k) * ones (size (z), "like", z);
    dL12 = -(1 ./ t) * ones (size (z), "like", z);
    dL22 = -(2 .* z - k) ./ (t .^ 2);
    ## Calculate censored data
    if (n_censored > 0)
      ## Compute derivatives
      [y, dy, d2y] = dgammainc (z_censored, k);
      dlnS = dy ./ y;
      d2lnS = d2y ./ y - dlnS.*dlnS;

      #[dlnS,d2lnS] = dlngamsf(z_censored,k);
      logzcen = log(z_censored);
      tmp = exp(k .* logzcen - z_censored - gammaln(k) - log(t)) ./ Scen;
      dL11(logical (censor)) = d2lnS;
      dL12(logical (censor)) = tmp .* (logzcen - dlnS - psi(0,k));
      dL22(logical (censor)) = tmp .* ((z_censored-1-k)./t - tmp);
    endif
    nH11 = -sum(freq .* dL11);
    nH12 = -sum(freq .* dL12);
    nH22 = -sum(freq .* dL22);
    nH = [nH11 nH12; nH12 nH22];
    if (any (isnan (nH(:))))
      acov = nan (2, "like", nH);
    else
      acov = inv (nH);
    endif
  endif

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
%! [nlogL, acov] = gamlike([2, 3], [2, 3, 4, 5, 6, 7, 8, 9]);
%! assert (nlogL, 19.4426, 1e-4);
%! assert (acov, [2.7819, -5.0073; -5.0073, 9.6882], 1e-4);
%!test
%! [nlogL, acov] = gamlike([2, 3], [5:45]);
%! assert (nlogL, 305.8070, 1e-4);
%! assert (acov, [0.0423, -0.0087; -0.0087, 0.0167], 1e-4);
%!test
%! [nlogL, acov] = gamlike([2, 13], [5:45]);
%! assert (nlogL, 163.2261, 1e-4);
%! assert (acov, [0.2362, -1.6631; -1.6631, 13.9440], 1e-4);

## Test input validation
%!error<gamlike: function called with too few input arguments.> ...
%! gamlike ([12, 15])
%!error<gamlike: wrong parameters length.> gamlike ([12, 15, 3], [1:50])
%!error<gamlike: X must be a vector.> gamlike ([12, 3], ones (10, 2))
%!error<gamlike: X and CENSOR vectors mismatch.> ...
%! gamlike ([12, 15], [1:50], [1, 2, 3])
%!error<gamlike: X and FREQ vectors mismatch.> ...
%! gamlike ([12, 15], [1:50], [], [1, 2, 3])
