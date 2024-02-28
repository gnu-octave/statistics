## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{nlogL} =} ricelike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} ricelike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} ricelike (@var{params}, @var{x}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} ricelike (@var{params}, @var{x}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the Rician distribution.
##
## @code{@var{nlogL} = ricelike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the Rician
## distribution with (1) non-centrality (distance) parameter @math{s} and (2)
## scale parameter @math{sigma} given in the two-element vector @var{params}.
##
## @code{[@var{nlogL}, @var{acov}] = ricelike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{params} are their asymptotic variances.
##
## @code{[@dots{}] = ricelike (@var{params}, @var{x}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = ricelike (@var{params}, @var{x}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the Rician distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rice_distribution}
##
## @seealso{ricecdf, riceinv, ricepdf, ricernd, ricefit, ricestat}
## @end deftypefn

function [nlogL, acov] = ricelike (params, x, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("ricelike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("ricelike: X must be a vector.");
  endif

  if (length (params) != 2)
    error ("ricelike: PARAMS must be a two-element vector.");
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("ricelike: X and CENSOR vector mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("ricelike: X and FREQ vector mismatch.");
  elseif (any (freq < 0))
    error ("ricelike: FREQ must not contain negative values.");
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
  nu = params(1);
  sigma = params(2);

  theta = nu ./ sigma;
  xsigma = x ./ sigma;
  xstheta = xsigma.*theta;
  I_0 = besseli (0, xstheta, 1);
  XNS = (xsigma .^ 2 + theta .^ 2) ./ 2;
  ## Compute log likelihood
  L = -XNS + log (I_0) + log (xsigma ./ sigma) + xstheta;

  ## Handle censored data
  n_censored = sum (freq .* censor);
  if (n_censored > 0)
    censored = (censor == 1);
    xsigma_censored = xsigma(censored);
    Q = marcumQ1 (theta, xsigma_censored);
    L(censored) = log (Q);
  endif
  ## Sum up the neg log likelihood
  nlogL = -sum (freq .* L);

  ## Compute asymptotic covariance
  if (nargout > 1)
    ## Compute first order central differences of the log-likelihood gradient
    dp = 0.0001 .* max (abs (params), 1);

    ngrad_p1 = rice_grad (params + [dp(1), 0], x, censor, freq);
    ngrad_m1 = rice_grad (params - [dp(1), 0], x, censor, freq);
    ngrad_p2 = rice_grad (params + [0, dp(2)], x, censor, freq);
    ngrad_m2 = rice_grad (params - [0, dp(2)], x, censor, freq);

    ## Compute negative Hessian by normalizing the differences by the increment
    nH = [(ngrad_p1(:) - ngrad_m1(:))./(2 * dp(1)), ...
          (ngrad_p2(:) - ngrad_m2(:))./(2 * dp(2))];

    ## Force neg Hessian being symmetric
    nH = 0.5 .* (nH + nH');
    ## Check neg Hessian is positive definite
    [R, p] = chol (nH);
    if (p > 0)
      warning ("ricelike: non positive definite Hessian matrix.");
      acov = NaN (2);
      return
    endif
    ## ACOV estimate is the negative inverse of the Hessian.
    Rinv = inv (R);
    acov = Rinv * Rinv;
  endif

endfunction

## Helper function for computing negative gradient
function ngrad = rice_grad (params, x, censor, freq)

  ## Get parameters
  nu = params(1);
  sigma = params(2);
  theta = nu ./ sigma;
  xsigma = x ./ sigma;
  xstheta = xsigma.*theta;
  I_0 = besseli (0, xstheta, 1);
  XNS = (xsigma .^ 2 + theta .^ 2) ./ 2;

  ## Compute derivatives
  I_1 = besseli(1, xstheta, 1);
  dII = I_1 ./ I_0;
  dL1 = (-theta + dII .* xsigma) ./ sigma;
  dL2 = -2 * (1 - XNS + dII .* xstheta) ./ sigma;

  ## Handle censored data
  n_censored = sum (freq .* censor);
  if (n_censored > 0)
    censored = (censor == 1);
    xsigma_censored = xsigma(censored);
    Q = marcumQ1 (theta, xsigma_censored);
    expt = exp (-XNS(censored) + xstheta(censored));
    dQdtheta = xsigma_censored .* I_1(censored) .* expt;
    dQdz = -xsigma_censored .* I_0(censored) .* expt;
    dtheta1 = 1 ./ sigma;
    dtheta2 = -theta ./ sigma;
    dz2 = -xsigma_censored ./ sigma;
    dL1(censored) = dQdtheta .* dtheta1 ./ Q;
    dL2(censored) = (dQdtheta .* dtheta2 + dQdz .* dz2) ./ Q;
  endif

  ## Compute gradient
  ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
endfunction

## Marcum's "Q" function of order 1
function Q = marcumQ1 (a, b)

  ## Prepare output matrix
  if (isa (a, "single") || isa (b, "single"))
   Q = NaN (size (b), "single");
  else
   Q = NaN (size (b));
  endif

  ## Force marginal cases
  Q(a != Inf & b == 0) = 1;
  Q(a != Inf & b == Inf) = 0;
  Q(a == Inf & b != Inf) = 1;
  z = isnan (Q) & a == 0 & b != Inf;
  if (any(z))
    Q(z) = exp ((-b(z) .^ 2) ./ 2);
  endif

  ## Compute the remaining cases
  z = isnan (Q) & ! isnan (a) & ! isnan (b);
  if (any(z(:)))
    aa = (a(z) .^ 2) ./ 2;
    bb = (b(z) .^ 2) ./ 2;
    eA = exp (-aa);
    eB = bb .* exp (-bb);
    h = eA;
    d = eB .* h;
    s = d;
    j = (d > s.*eps(class(d)));
    k = 1;
    while (any (j))
      eA = aa .* eA ./ k;
      h = h + eA;
      eB = bb .* eB ./ (k + 1);
      d = eB .* h;
      s(j) = s (j) + d(j);
      j = (d > s .* eps (class (d)));
      k = k + 1;
    endwhile
    Q(z) = 1 - s;
  endif
endfunction

## Test output
%!test
%! nlogL = ricelike ([15.3057344, 17.6668458], [1:50]);
%! assert (nlogL, 204.5230311010569, 1e-12);
%!test
%! nlogL = ricelike ([2.312346885, 1.681228265], [1:5]);
%! assert (nlogL, 8.65562164930058, 1e-12);

## Test input validation
%!error<ricelike: function called with too few input arguments.> ricelike (3.25)
%!error<ricelike: X must be a vector.> ricelike ([5, 0.2], ones (2))
%!error<ricelike: PARAMS must be a two-element vector.> ...
%! ricelike ([1, 0.2, 3], [1, 3, 5, 7])
%!error<ricelike: X and CENSOR vector mismatch.> ...
%! ricelike ([1.5, 0.2], [1:5], [0, 0, 0])
%!error<ricelike: X and FREQ vector mismatch.> ...
%! ricelike ([1.5, 0.2], [1:5], [0, 0, 0, 0, 0], [1, 1, 1])
%!error<ricelike: X and FREQ vector mismatch.> ...
%! ricelike ([1.5, 0.2], [1:5], [], [1, 1, 1])
%!error<ricelike: FREQ must not contain negative values.> ...
%! ricelike ([1.5, 0.2], [1:5], [], [1, 1, 1, 0, -1])
