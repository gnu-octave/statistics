## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR l PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{nlogL} =} tlslike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{acov}] =} tlslike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} tlslike (@var{params}, @var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} tlslike (@var{params}, @var{x}, @var{alpha}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the location-scale Student's T distribution.
##
## @code{@var{nlogL} = tlslike (@var{params}, @var{x})} returns the negative
## log-likelihood of the x in @var{x} corresponding to the location-scale T
## distribution with (1) location parameter @math{mu}, (2) scale parameter
## @math{sigma} and (3) degrees of freedom @math{nu} given in the three-element
## vector @var{params}.
##
## @code{[@var{nlogL}, @var{acov}] = tlslike (@var{params}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{acov}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{acov} are their asymptotic variances.  @var{acov}
## is based on the observed Fisher's information, not the expected information.
##
## @code{[@dots{}] = tlslike (@var{params}, @var{x}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = tlslike (@var{params}, @var{x}, @var{censor},
## @var{freq})} accepts a frequency vector, @var{freq}, of the same size as
## @var{x}.  @var{freq} typically contains integer frequencies for the
## corresponding elements in @var{x}, but may contain any non-integer
## non-negative values.  By default, or if left empty,
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## Further information about the location-scale Student's T distribution can be
## found at @url{https://en.wikipedia.org/wiki/Student%27s_t-distribution#Location-scale_t_distribution}
##
## @seealso{tlscdf, tlsinv, tlspdf, tlsrnd, tlsfit, tlsstat}
## @end deftypefn

function [nlogL, acov] = tlslike (params, x, censor, freq)

  ## Check input arguments and add defaults
  if (nargin < 2)
    error ("tlslike: too few input arguments.");
  endif
  if (numel (params) != 3)
    error ("tlslike: wrong parameters length.");
  endif
  if (! isvector (x))
    error ("tlslike: X must be a vector.");
  endif
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("tlslike: X and CENSOR vectors mismatch.");
  endif
  if (nargin < 4 || isempty (freq))
    freq = ones (size (x));
  elseif (any (freq < 0))
    error ("tlslike: FREQ cannot have negative values.");
  elseif (isequal (size (x), size (freq)))
    nulls = find (freq == 0);
    if (numel (nulls) > 0)
      x(nulls) = [];
      censor(nulls) = [];
      freq(nulls) = [];
    endif
  else
    error ("tlslike: X and FREQ vectors mismatch.");
  endif

  ## Compute the negative log-likelihood
  nlogL = tlsnll (x, params, censor, freq);

  ## Compute the negative hessian and invert to get the information matrix
  if (nargout > 1)
    ei = zeros (1, 3);
    ej = zeros (1, 3);
    nH = zeros (3, 3);
    dp = (eps ^ (1/4)) .* max (abs (params), 1);
    for i = 1:3
      ei(i) = dp(i);
      for j = 1:(i-1)
        ej(j) = dp(j);
        ## Four-point central difference for mixed second partials
        nH(i,j) = tlsnll (x, params+ei+ej, censor, freq) ...
                - tlsnll (x, params+ei-ej, censor, freq) ...
                - tlsnll (x, params-ei+ej, censor, freq) ...
                + tlsnll (x, params-ei-ej, censor, freq);
        ej(j) = 0;
      endfor
      ## Five-point central difference for pure second partial
      nH(i,i) = -  tlsnll (x, params+2*ei, censor, freq) ...
                + 16 * tlsnll (x, params+ei, censor, freq) - 30 * nlogL ...
                + 16 * tlsnll (x, params-ei, censor, freq) ...
                - tlsnll (x, params-2*ei, censor, freq);
      ei(i) = 0;
    endfor

    ## Fill in the upper triangle
    nH = nH + triu (nH', 1);

    ## Normalize the second differences to get derivative estimates
    nH = nH ./ (4 .* dp(:) * dp(:)' + diag (8 * dp(:) .^ 2));

    ## Check neg Hessian is positive definite
    [R, p] = chol (nH);
    if (p > 0)
      warning ("tlslike: non positive definite Hessian matrix.");
      acov = NaN (3);
      return
    endif

    ## ACOV estimate is the negative inverse of the Hessian
    Rinv = inv (R);
    acov = Rinv * Rinv';
  endif

endfunction

## Internal function to calculate negative log likelihood for tlslike
function nlogL = tlsnll (x, params, censor, freq)

  mu = params(1);
  sigma = params(2);
  sigma(sigma <= 0) = NaN;
  nu = params(3);
  nu(nu <= 0) = NaN;

  z = (x - mu) ./ sigma;
  w = nu + (z .^ 2);
  logw = log (w);

  L = - 0.5 .* (nu + 1) .* logw + gammaln (0.5 .* (nu + 1)) ...
      - gammaln (0.5 .* nu) + 0.5 .* nu .* log (nu) ...
      - log (sigma) - 0.5 .* log (pi);

  n_censored = sum (freq .* censor);
  if (n_censored > 0)
    censored = (censor == 1);
    if (nu < 1e7)   # Use incomplete beta function
      S_censored = betainc (nu ./ w(censored), 0.5 .* nu, 0.5) ./ 2;
      S_censored(z(censored) < 0) = 1 - S_censored(z(censored) < 0);
    else            # Use a normal approximation
      S_censored = log(0.5 * erfc(z(censored) ./ sqrt(2)));
    endif
    L(censored) = log(S_censored);
  endif
  nlogL = - sum (freq .* L);

endfunction

## Test output
%!test
%! x = [-1.2352, -0.2741, 0.1726, 7.4356, 1.0392, 16.4165];
%! [nlogL, acov] = tlslike ([0.035893, 0.862711, 0.649261], x);
%! acov_out = [0.2525, 0.0670, 0.0288; ...
%!             0.0670, 0.5724, 0.1786; ...
%!             0.0288, 0.1786, 0.1789];
%! assert (nlogL, 17.9979636579, 1e-10);
%! assert (acov, acov_out, 1e-4);

## Test input validation
%!error<tlslike: too few input arguments.> tlslike ([12, 15, 1]);
%!error<tlslike: wrong parameters length.> tlslike ([12, 15], [1:50]);
%!error<tlslike: X must be a vector.> tlslike ([12, 3, 1], ones (10, 2));
%!error<tlslike: X and CENSOR> tlslike ([12, 15, 1], [1:50], [1, 2, 3]);
%!error<tlslike: X and FREQ> tlslike ([12, 15, 1], [1:50], [], [1, 2, 3]);
%!error<tlslike: FREQ cannot> tlslike ([12, 15, 1], [1:3], [], [1, 2, -3]);
