## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## the two-element vector @var{paramhat}.
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
    L(cens) = log (Scen);
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
      [y, dy, d2y] = dgammainc (z_censored, k, "upper");
      dlnS = dy ./ y;
      d2lnS = d2y ./ y - dlnS.*dlnS;

      #[dlnS,d2lnS] = dlngamsf(z_censored,k);
      logzcen = log(z_censored);
      tmp = exp(k .* logzcen - z_censored - gammaln(k) - log(t)) ./ Scen;
      dL11(cens) = d2lnS;
      dL12(cens) = tmp .* (logzcen - dlnS - psi(0,k));
      dL22(cens) = tmp .* ((z_censored-1-k)./t - tmp);
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
