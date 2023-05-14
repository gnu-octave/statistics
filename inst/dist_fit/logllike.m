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
## @deftypefn  {statistics} {@var{nlogL} =} logllike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} logllike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} logllike (@var{params}, @var{x}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} logllike (@var{params}, @var{x}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the log-logistic distribution.
##
## @code{@var{nlogL} = logllike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the log-logistic
## distribution with (1) scale parameter @var{a} and (2) shape parameter @var{b}
## given in the two-element vector @var{params}.
##
## @code{[@var{nlogL}, @var{acov}] = logllike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{params} are their asymptotic variances.
##
## @code{[@dots{}] = logllike (@var{params}, @var{x}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = logllike (@var{params}, @var{x}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the log-logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-logistic_distribution}
##
## MATLAB compatibility: MATLAB uses an alternative parameterization given by
## the pair @math{μ, s}, i.e. @var{mu} and @var{s}, in analogy with the logistic
## distribution.  Their relation to the @var{a} and @var{b} parameters is given
## below:
##
## @itemize
## @item @qcode{@var{a} = exp (@var{mu})}
## @item @qcode{@var{b} = 1 / @var{s}}
## @end itemize
##
## @seealso{loglcdf, loglinv, loglpdf, loglrnd, loglfit}
## @end deftypefn

function [nlogL, acov] = logllike (params, x, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("logllike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("logllike: X must be a vector.");
  endif

  if (length (params) != 2)
    error ("logllike: PARAMS must be a two-element vector.");
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("logllike: X and CENSOR vector mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("logllike: X and FREQ vector mismatch.");
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
  a = params(1);
  b = params(2);

  z = (log (x) - log (a)) .* b;
  logclogitz = log (1 ./ (1 + exp (z)));
  k = (z > 700);
  if (any (k))
    logclogitz(k) = z(k);
  endif

  L = z + 2 .* logclogitz - log (1 / b) - log (x);
  n_censored = sum (freq .* censor);
  ## Handle censored data
  if (n_censored > 0)
    censored = (censor == 1);
    L(censored) = logclogitz(censored);
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
      warning ("logllike: non positive definite Hessian matrix.");
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
  a = params(1);
  b = params(2);
  z = (log (x) - log (a)) .* b;
  logitz = 1 ./ (1 + exp (-z));
  dL1 = (2 .* logitz - 1) .* b;
  dL2 = z .* dL1 - b;
  n_censored = sum (freq .* censor);
  if (n_censored > 0)
    censored = (censor == 1);
    dL1(censored) = logitz(censored) .* b;
    dL2(censored) = z(censored) .* dL1(censored);
  endif
  ngrad = -[sum(freq .* dL1), sum(freq .* dL2)];
endfunction


## Test results
%!test
%! nlogL = logllike ([exp(3.09717), 1/0.468525], [1:50]);
%! assert (nlogL, 211.2965, 1e-4);
%!test
%! nlogL = logllike ([exp(1.01124), 1/0.336449], [1:5]);
%! assert (nlogL, 9.2206, 1e-4);

## Test input validation
%!error<logllike: function called with too few input arguments.> logllike (3.25)
%!error<logllike: X must be a vector.> logllike ([5, 0.2], ones (2))
%!error<logllike: PARAMS must be a two-element vector.> ...
%! logllike ([1, 0.2, 3], [1, 3, 5, 7])
%!error<logllike: X and CENSOR vector mismatch.> ...
%! logllike ([1.5, 0.2], [1:5], [0, 0, 0])
%!error<logllike: X and FREQ vector mismatch.> ...
%! logllike ([1.5, 0.2], [1:5], [0, 0, 0, 0, 0], [1, 1, 1])
%!error<logllike: X and FREQ vector mismatch.> ...
%! logllike ([1.5, 0.2], [1:5], [], [1, 1, 1])
