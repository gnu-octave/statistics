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
## @deftypefn  {statistics} {@var{nlogL} =} bisalike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} bisalike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} bisalike (@var{params}, @var{x}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} bisalike (@var{params}, @var{x}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the Birnbaum-Saunders distribution.
##
## @code{@var{nlogL} = bisalike (@var{params}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the Birnbaum-Saunders
## distribution with (1) scale parameter @var{beta} and (2) shape parameter
## @var{gamma} given in the two-element vector @var{paramhat}.
##
## @code{[@var{nlogL}, @var{acov}] = bisalike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{params} are their asymptotic variances.
##
## @code{[@dots{}] = bisalike (@var{params}, @var{x}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = bisalike (@var{params}, @var{x}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the Birnbaum-Saunders distribution can be found at
## @url{https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution}
##
## @seealso{bisacdf, bisainv, bisapdf, bisarnd, bisafit, bisastat}
## @end deftypefn

function [nlogL, acov] = bisalike (params, x, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("bisalike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("bisalike: X must be a vector.");
  endif

  if (any (x < 0))
    error ("bisalike: X cannot have negative values.");
  endif

  if (length (params) != 2)
    error ("bisalike: PARAMS must be a two-element vector.");
  endif

  ## Check censor vector
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("bisalike: X and CENSOR vector mismatch.");
  endif

  ## Check frequency vector
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (! isequal (size (x), size (freq)))
    error ("bisalike: X and FREQ vector mismatch.");
  endif

  beta = params(1);
  gamma = params(2);
  z = (sqrt (x ./ beta) - sqrt (beta ./ x)) ./ gamma;
  w = (sqrt (x ./ beta) + sqrt (beta ./ x)) ./ gamma;

  L = -0.5 .* (z .^ 2 + log (2 .* pi)) + log (w) - log (2 .* x);
  n_censored = sum (freq .* censor);

  if (n_censored > 0)
    censored = (censor == 1);
    z_censored = z(censored);
    Scen = 0.5 * erfc (z_censored ./ sqrt(2));
    L(censored) = log (Scen);
  endif
  nlogL = -sum (freq .* L);

  ## Compute asymptotic covariance
  if (nargout > 1)
    ## Compute first order central differences of the log-likelihood gradient
    dp = 0.0001 .* max (abs (params), 1);

    ngrad_p1 = bisa_ngrad (params + [dp(1), 0], x, censor, freq);
    ngrad_m1 = bisa_ngrad (params - [dp(1), 0], x, censor, freq);
    ngrad_p2 = bisa_ngrad (params + [0, dp(2)], x, censor, freq);
    ngrad_m2 = bisa_ngrad (params - [0, dp(2)], x, censor, freq);

    ## Compute negative Hessian by normalizing the differences by the increment
    nH = [(ngrad_p1(:) - ngrad_m1(:))./(2 * dp(1)), ...
          (ngrad_p2(:) - ngrad_m2(:))./(2 * dp(2))];

    ## Force neg Hessian being symmetric
    nH = 0.5 .* (nH + nH');
    ## Check neg Hessian is positive definite
    [R, p] = chol (nH);
    if (p > 0)
      warning ("bisalike: non positive definite Hessian matrix.");
      acov = NaN (2);
      return
    endif
    ## ACOV estimate is the negative inverse of the Hessian.
    Rinv = inv (R);
    acov = Rinv * Rinv;
  endif

endfunction

## Helper function for computing negative gradient
function ngrad = bisa_ngrad (params, x, censor, freq)
  beta = params(1);
  gamma = params(2);
  z = (sqrt (x ./ beta) - sqrt (beta ./ x)) ./ gamma;
  w = (sqrt (x ./ beta) + sqrt (beta ./ x)) ./ gamma;
  logphi = -0.5 .* (z .^ 2 + log (2 .* pi));
  n_censored = sum (freq .* censor);
  if (n_censored > 0)
    censored = (censor == 1);
    z_censored = z(censored);
    Scen = 0.5 * erfc (z_censored ./ sqrt(2));
  endif
  dL1 = (w .^ 2 - 1) .* 0.5 .* z ./ (w .* beta);
  dL2 = (z .^ 2 - 1) ./ gamma;
  if (n_censored > 0)
    phi_censored = exp (logphi(censored));
    wcen = w(censored);
    d1Scen = phi_censored .* 0.5 .* wcen ./ beta;
    d2Scen = phi_censored .* z_censored ./ gamma;
    dL1(censored) = d1Scen ./ Scen;
    dL2(censored) = d2Scen ./ Scen;
  endif
  ngrad = -[sum(freq .* dL1), sum(freq .* dL2)];
endfunction


## Test results
%!test
%! nlogL = bisalike ([16.2649, 1.0156], [1:50]);
%! assert (nlogL, 215.5905, 1e-4);
%!test
%! nlogL = bisalike ([2.5585, 0.5839], [1:5]);
%! assert (nlogL, 8.9950, 1e-4);

## Test input validation
%!error<bisalike: function called with too few input arguments.> bisalike (3.25)
%!error<bisalike: X must be a vector.> bisalike ([5, 0.2], ones (2))
%!error<bisalike: X cannot have negative values.> bisalike ([5, 0.2], [-1, 3])
%!error<bisalike: PARAMS must be a two-element vector.> ...
%! bisalike ([1, 0.2, 3], [1, 3, 5, 7])
%!error<bisalike: X and CENSOR vector mismatch.> ...
%! bisalike ([1.5, 0.2], [1:5], [0, 0, 0])
%!error<bisalike: X and FREQ vector mismatch.> ...
%! bisalike ([1.5, 0.2], [1:5], [0, 0, 0, 0, 0], [1, 1, 1])
%!error<bisalike: X and FREQ vector mismatch.> ...
%! bisalike ([1.5, 0.2], [1:5], [], [1, 1, 1])
