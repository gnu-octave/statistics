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
## Further information about the loglogistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-logistic_distribution}
##
## OCTAVE/MATLAB use an alternative parameterization given by the pair
## @math{μ, σ}, i.e. @var{mu} and @var{sigma}, in analogy with the logistic
## distribution.  Their relation to the @math{α} and @math{b} parameters used
## in Wikipedia are given below:
##
## @itemize
## @item @qcode{@var{mu} = log (@var{a})}
## @item @qcode{@var{sigma} = 1 / @var{a}}
## @end itemize
##
## @seealso{loglcdf, loglinv, loglpdf, loglrnd, loglfit, loglstat}
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

  ## Compute the negative loglikelihood
  nlogL = loglnll (params, x, censor, freq);

  ## Compute the negative hessian and invert to get the information matrix
  if (nargout > 1)
    ei = zeros (1, 2);
    ej = zeros (1, 2);
    nH = zeros (2, 2);
    dp = (eps ^ (1/4)) .* max (abs (params), 1);
    for i = 1:2
      ei(i) = dp(i);
      for j = 1:(i-1)
        ej(j) = dp(j);
        ## Four-point central difference for mixed second partials
        nH(i,j) = loglnll (params+ei+ej, x, censor, freq) ...
                - loglnll (params+ei-ej, x, censor, freq) ...
                - loglnll (params-ei+ej, x, censor, freq) ...
                + loglnll (params-ei-ej, x, censor, freq);
        ej(j) = 0;
      endfor
      ## Five-point central difference for pure second partial
      nH(i,i) = -  loglnll (params+2*ei, x, censor, freq) ...
                + 16 * loglnll (params+ei, x, censor, freq) - 30 * nlogL ...
                + 16 * loglnll (params-ei, x, censor, freq) ...
                - loglnll (params-2*ei, x, censor, freq);
      ei(i) = 0;
    endfor

    ## Fill in the upper triangle
    nH = nH + triu (nH', 1);

    ## Normalize the second differences to get derivative estimates
    nH = nH ./ (4 .* dp(:) * dp(:)' + diag (8 * dp(:) .^ 2));

    ## Check neg Hessian is positive definite
    [R, p] = chol (nH);
    if (p > 0)
      warning ("logllike: non positive definite Hessian matrix.");
      acov = NaN (2);
      return
    endif

    ## ACOV estimate is the negative inverse of the Hessian
    Rinv = inv (R);
    acov = Rinv * Rinv';
  endif

endfunction

## Helper function for computing negative loglikelihood
function nlogL = loglnll (params, x, censor, freq)
  ## Get parameters
  mu = params(1);
  sigma = params(2);
  ## Compute intermediate values
  z = (log (x) - mu) ./ sigma;
  logclogitz = log (1 ./ (1 + exp (z)));
  k = (z > 708);
  if (any (k))
    logclogitz(k) = z(k);
  endif
  L = z + 2 .* logclogitz - log (sigma) - log (x);
  n_censored = sum (freq .* censor);
  ## Handle censored data
  if (n_censored > 0)
    censored = (censor == 1);
    L(censored) = logclogitz(censored);
  endif
  ## Sum up the neg log likelihood
  nlogL = -sum (freq .* L);
endfunction

## Test output
%!test
%! [nlogL, acov] = logllike ([3.09717, 0.468525], [1:50]);
%! assert (nlogL, 211.2965, 1e-4);
%! assert (acov, [0.0131, -0.0007; -0.0007, 0.0031], 1e-4);
%!test
%! [nlogL, acov] = logllike ([1.01124, 0.336449], [1:5]);
%! assert (nlogL, 9.2206, 1e-4);
%! assert (acov, [0.0712, -0.0032; -0.0032, 0.0153], 1e-4);

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
