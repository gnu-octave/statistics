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
## @deftypefn  {statistics} {@var{nlogL} =} invglike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} invglike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} invglike (@var{params}, @var{x}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} invglike (@var{params}, @var{x}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the inverse Gaussian distribution.
##
## @code{@var{nlogL} = invglike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the inverse Gaussian
## distribution with (1) scale parameter @var{mu} and (2) shape parameter
## @var{lambda} given in the two-element vector @var{params}.
##
## @code{[@var{nlogL}, @var{acov}] = invglike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{params} are their asymptotic variances.
##
## @code{[@dots{}] = invglike (@var{params}, @var{x}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = invglike (@var{params}, @var{x}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the inverse Gaussian distribution can be found at
## @url{https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}
##
## @seealso{invgcdf, invginv, invgpdf, invgrnd, invgfit}
## @end deftypefn

function [nlogL, acov] = invglike (params, x, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("invglike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("invglike: X must be a vector.");
  endif

  if (any (x < 0))
    error ("invglike: X must have positive values.");
  endif

  if (length (params) != 2)
    error ("invglike: PARAMS must be a two-element vector.");
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("invglike: X and CENSOR vector mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("invglike: X and FREQ vector mismatch.");
  endif

  ## Get parameters
  mu = params(1);
  lambda = params(2);

  L = 0.5 .* log (lambda ./ (2 * pi)) - 1.5 .* log (x) ...
      -lambda .* (x ./ mu - 1) .^ 2 ./ (2 .* x);
  n_censored = sum (freq .* censor);

  ## Handle censored data
  if (n_censored > 0)
    censored = (censor == 1);
    x_censored = x(censored);
    sqrt_lx = sqrt (lambda ./ x_censored);
    z_censored = -(x_censored ./ mu - 1) .* sqrt_lx;
    w_censored = -(x_censored ./ mu + 1) .* sqrt_lx;
    Fz = 0.5 .* erfc (-z_censored ./ sqrt (2));
    Fw = 0.5 .* erfc (-w_censored ./ sqrt (2));
    S_censored = Fz - exp (2 .* lambda ./ mu) .* Fw;
    L(censored) = log (S_censored);
  endif

  ## Sum up the neg log likelihood
  nlogL = -sum (freq .* L);

  ## Compute asymptotic covariance
  if (nargout > 1)
    ## Compute first order central differences of the log-likelihood gradient
    dp = 0.0001 .* max (abs (params), 1);

    ngrad_p1 = invg_grad (params + [dp(1), 0], x, censor, freq);
    ngrad_m1 = invg_grad (params - [dp(1), 0], x, censor, freq);
    ngrad_p2 = invg_grad (params + [0, dp(2)], x, censor, freq);
    ngrad_m2 = invg_grad (params - [0, dp(2)], x, censor, freq);

    ## Compute negative Hessian by normalizing the differences by the increment
    nH = [(ngrad_p1(:) - ngrad_m1(:))./(2 * dp(1)), ...
          (ngrad_p2(:) - ngrad_m2(:))./(2 * dp(2))];

    ## Force neg Hessian being symmetric
    nH = 0.5 .* (nH + nH');
    ## Check neg Hessian is positive definite
    [R, p] = chol (nH);
    if (p > 0)
      warning ("invglike: non positive definite Hessian matrix.");
      acov = NaN (2);
      return
    endif
    ## ACOV estimate is the negative inverse of the Hessian.
    Rinv = inv (R);
    acov = Rinv * Rinv;
  endif

endfunction

## Helper function for computing negative gradient
function ngrad = invg_grad (params, x, censor, freq)
  mu = params(1);
  lambda = params(2);
  dL1 = lambda .* (x - mu) ./ mu .^ 3;
  dL2 = 1 ./ (2 .* lambda) - (x ./ mu - 1) .^ 2 ./ (2 .* x);
  n_censored = sum (freq .* censor);
  if (n_censored > 0)
    censored = (censor == 1);
    x_censored = x(censored);
    sqrt_lx = sqrt (lambda ./ x_censored);
    exp_lmu = exp (2 .* lambda ./ mu);
    z_censored = -(x_censored ./ mu - 1) .* sqrt_lx;
    w_censored = -(x_censored ./ mu + 1) .* sqrt_lx;
    Fw = 0.5 .* erfc (-w_censored ./ sqrt (2));
    fz = exp (-0.5 .* z_censored .^ 2) ./ sqrt (2 .* pi);
    fw = exp (-0.5 .* w_censored .^ 2) ./ sqrt (2 .* pi);
    dS1cen = (fz - exp_lmu .* fw) .* (x_censored ./ mu .^ 2) .* sqrt_lx ...
                              + 2 .* Fw .* exp_lmu .* lambda ./ mu .^ 2;
    dS2cen = 0.5 .* (fz .* z_censored - exp_lmu .* fw .* w_censored) ...
                 ./ lambda - 2 .* Fw .* exp_lmu ./ mu;
    dL1(cen) = dS1cen ./ Scen;
    dL2(cen) = dS2cen ./ Scen;
  endif
  ngrad = -[sum(freq .* dL1), sum(freq .* dL2)];
endfunction


## Test results
%!test
%! nlogL = invglike ([25.5, 19.6973], [1:50]);
%! assert (nlogL, 219.1516, 1e-4);
%!test
%! nlogL = invglike ([3, 8.1081], [1:5]);
%! assert (nlogL, 9.0438, 1e-4);

## Test input validation
%!error<invglike: function called with too few input arguments.> invglike (3.25)
%!error<invglike: X must be a vector.> invglike ([5, 0.2], ones (2))
%!error<invglike: X must have positive values.> invglike ([5, 0.2], [-1, 3])
%!error<invglike: PARAMS must be a two-element vector.> ...
%! invglike ([1, 0.2, 3], [1, 3, 5, 7])
%!error<invglike: X and CENSOR vector mismatch.> ...
%! invglike ([1.5, 0.2], [1:5], [0, 0, 0])
%!error<invglike: X and FREQ vector mismatch.> ...
%! invglike ([1.5, 0.2], [1:5], [0, 0, 0, 0, 0], [1, 1, 1])
%!error<invglike: X and FREQ vector mismatch.> ...
%! invglike ([1.5, 0.2], [1:5], [], [1, 1, 1])
