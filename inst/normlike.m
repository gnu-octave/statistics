## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {statistics} @var{nlogL} = normlike (@var{params}, @var{data})
## @deftypefnx {statistics} [@var{nlogL}, @var{avar}] = normlike (@var{params}, @var{data})
## @deftypefnx {statistics} [@dots{}] = normlike (@var{params}, @var{data}, @var{censor})
## @deftypefnx {statistics} [@dots{}] = normlike (@var{params}, @var{data}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the normal distribution.
##
## @code{@var{nlogL} = normlike (@var{params}, @var{data})} returns the
## negative of the log-likelihood for the normal distribution, evaluated at
## parameters @var{params(1)} = mean and @var{params(2)} = standard deviation,
## given @var{data}.  @var{nlogL} is a scalar.
##
## @code{[@var{nlogL}, @var{avar}] = normlike (@var{params}, @var{data})}
## returns the inverse of Fisher's information matrix, @var{avar}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{avar} are their asymptotic variances.  @var{avar}
## is based on the observed Fisher's information, not the expected information.
##
## @code{[@dots{}] = normlike (@var{params}, @var{data}, @var{censor})} accepts
## a boolean vector of the same size as @var{data} that is 1 for observations
## that are right-censored and 0 for observations that are observed exactly.
##
## @code{[@dots{}] = normlike (@var{params}, @var{data}, @var{censor},
## @var{freq})} accepts a frequency vector of the same size as @var{data}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{data}, but it may contain any non-integer non-negative
## values.  Pass in [] for @var{censor} to use its default value.
##
## @seealso{normcdf, norminv, normpdf, normrnd, normfit, normstat}
## @end deftypefn

function [nlogL, avar] = normlike (params ,data, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("normlike: too few input arguments.");
  endif
  if (! isvector (data))
    error ("normlike: DATA must be a vector.");
  endif
  if (numel (params) != 2)
    error ("normlike: PARAMS must be a two-element vector.");
  endif
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (data));
  elseif (! isequal (size (data), size (censor)))
    error ("normlike: DATA and CENSOR vectors mismatch.");
  endif
  if nargin < 4 || isempty(freq)
    freq = ones (size (data));
  elseif (isequal (size (data), size (freq)))
    nulls = find (freq == 0);
    if (numel (nulls) > 0)
      data(nulls) = [];
      censor(nulls) = [];
      freq(nulls) = [];
    endif
  else
    error ("normlike: DATA and FREQ vectors mismatch.");
  endif

  ## Get mu and sigma values
  mu = params(1);
  sigma = params(2);
  ## sigma must be positive, otherwise make it NaN
  if (sigma <= 0)
    sigma = NaN;
  endif

  ## Compute the individual log-likelihood terms.  Force a log(0)==-Inf for
  ## data from extreme right tail, instead of getting exp(Inf-Inf)==NaN.
  z = (data - mu) ./ sigma;
  L = -0.5 .* z .^ 2 - log (sqrt (2 .* pi) .* sigma);
  if (any (censor))
    censored = censor == 1;
    z_censor = z(censored);
    S_censor = 0.5 * erfc (z_censor / sqrt (2));
    L(censored) = log (S_censor);
  endif
  ## Neg-log-like is the sum of the individual contributions
  nlogL = -sum (freq .* L);

  ## Compute the negative hessian at the parameter values.
  ## Invert to get the observed information matrix.
  if (nargout == 2)
    dL11 = -ones (size (z), class (z));
    dL12 = -2 .* z;
    dL22 = 1 - 3 .* z .^ 2;
    if (any (censor))
      dlogScen = exp (-0.5 .* z_censor .^ 2) ./ (sqrt (2 * pi) .* S_censor);
      d2logScen = dlogScen .* (dlogScen - z_censor);
      dL11(censored) = -d2logScen;
      dL12(censored) = -dlogScen - z_censor .* d2logScen;
      dL22(censored) = -z_censor .* (2 .* dlogScen + z_censor .* d2logScen);
    endif
    nH11 = -sum (freq .* dL11);
    nH12 = -sum (freq .* dL12);
    nH22 = -sum (freq .* dL22);
    avar =  (sigma .^ 2) * [nH22, -nH12; -nH12, nH11] / ...
                           (nH11 * nH22 - nH12 * nH12);
  endif

endfunction

## Test input validation
%!error<normlike: too few input arguments.> normlike ([12, 15]);
%!error<normlike: DATA must be a vector.> normlike ([12, 15], ones (2));
%!error<normlike: PARAMS must be a two-element vector.> ...
%! normlike ([12, 15, 3], [1:50]);
%!error<normlike: DATA and CENSOR vectors mismatch.> ...
%! normlike ([12, 15], [1:50], [1, 2, 3]);
%!error<normlike: DATA and FREQ vectors mismatch.> ...
%! normlike ([12, 15], [1:50], [], [1, 2, 3]);

## Results compared with Matlab
%!test
%! data = 1:50;
%! [nlogL, avar] = normlike ([2.3, 1.2], data);
%! avar_out = [7.5767e-01, -1.8850e-02; -1.8850e-02, 4.8750e-04];
%! assert (nlogL, 13014.95883783327, 1e-10);
%! assert (avar, avar_out, 1e-4);
%!test
%! data = 1:50;
%! [nlogL, avar] = normlike ([2.3, 1.2], data * 0.5);
%! avar_out = [3.0501e-01, -1.5859e-02; -1.5859e-02, 9.1057e-04];
%! assert (nlogL, 2854.802587833265, 1e-10);
%! assert (avar, avar_out, 1e-4);
%!test
%! data = 1:50;
%! [nlogL, avar] = normlike ([21, 15], data);
%! avar_out = [5.460474308300396, -1.600790513833993; ...
%!             -1.600790513833993, 2.667984189723321];
%! assert (nlogL, 206.738325604233, 1e-12);
%! assert (avar, avar_out, 1e-14);
%!test
%! data = 1:50;
%! censor = ones (1, 50);
%! censor([2, 4, 6, 8, 12, 14]) = 0;
%! [nlogL, avar] = normlike ([2.3, 1.2], data, censor);
%! avar_out = [3.0501e-01, -1.5859e-02; -1.5859e-02, 9.1057e-04];
%! assert (nlogL, Inf);
%! assert (avar, [NaN, NaN; NaN, NaN]);
%!test
%! data = 1:50;
%! censor = ones (1, 50);
%! censor([2, 4, 6, 8, 12, 14]) = 0;
%! [nlogL, avar] = normlike ([21, 15], data, censor);
%! avar_out = [24.4824488866131, -10.6649544179636; ...
%!             -10.6649544179636, 6.22827849965737];
%! assert (nlogL, 86.9254371829733, 1e-12);
%! assert (avar, avar_out, 1e-14);
