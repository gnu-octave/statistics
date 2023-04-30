## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{nlogL} =} wbllike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{avar}] =} wbllike (@var{params}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} wbllike (@var{params}, @var{x}, @var{alpha}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} wbllike (@var{params}, @var{x}, @var{alpha}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the Weibull distribution.
##
## @code{@var{nlogL} = wbllike (@var{params}, @var{data})} returns the
## negative of the log-likelihood for the Weibull distribution given the data in
## @var{x}, evaluated at parameters @math{λ} and @math{k} given as a two-element
## vector @qcode{@var{paramhat}([1, 2])}, respectively.
##
## @code{[@var{nlogL}, @var{avar}] = wbllike (@var{params}, @var{data})} also
## returns the inverse of Fisher's information matrix, @var{avar}.  If the input
## parameter values in @var{params} are the maximum likelihood estimates, the
## diagonal elements of @var{avar} are their asymptotic variances.  @var{avar}
## is based on the observed Fisher's information, not the expected information.
##
## @code{[@dots{}] = wbllike (@var{params}, @var{data}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = wbllike (@var{params}, @var{data}, @var{censor},
## @var{freq})} accepts a frequency vector, @var{freq}, of the same size as
## @var{x}.  @var{freq} typically contains integer frequencies for the
## corresponding elements in @var{x}, but may contain any non-integer
## non-negative values.  By default, or if left empty,
## @qcode{@var{freq} = ones (size (@var{x}))}.
##
## @seealso{wblcdf, wblinv, wblpdf, wblrnd, wblfit, wblstat}
## @end deftypefn

function [nlogL, avar] = wbllike (params, x, censor, freq)

  ## Check input arguments and add defaults
  if (nargin < 2)
    error ("wbllike: too few input arguments.");
  endif
  if (numel (params) != 2)
    error ("wbllike: wrong parameters length.");
  endif
  if (! isvector (x))
    error ("wbllike: X must be a vector.");
  endif
  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("wbllike: X and CENSOR vectors mismatch.");
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
    error ("wbllike: X and FREQ vectors mismatch.");
  endif

  ## Get lambda and k parameter values
  l = params(1);
  k = params(2);

  ## Force NaNs for out of range parameters or x.
  l(l <= 0) = NaN;
  k(k <= 0) = NaN;
  x(x < 0) = NaN;

  ## Compute the individual log-likelihood terms
  z = x ./ l;
  logz = log (z);
  expz = exp (k .* logz);
  ilogL = ((k - 1) .* logz + log (k ./ l)) .* (1 - censor) - expz;
  ilogL(z == Inf) = -Inf;

  ## Sum up the individual log-likelihood contributions
  nlogL = -sum (freq .* ilogL);

  ## Compute the negative hessian and invert to get the information matrix.
  if (nargout > 1)
    ucen = (1 - censor);
    nH11 = sum (freq .* (k .* ((1 + k) .* expz - ucen))) ./ l .^ 2;
    nH12 = -sum(freq .* (((1 + k .* logz) .* expz - ucen))) ./ l;
    nH22 = sum(freq .* ((logz .^ 2) .* expz + ucen ./ k .^ 2));
    avar = [nH22, -nH12; -nH12, nH11] / (nH11 * nH22 - nH12 * nH12);
  endif

endfunction

## Results compared with Matlab
%!test
%! x = 1:50;
%! [nlogL, avar] = wbllike ([2.3, 1.2], x);
%! avar_out = [0.0250, 0.0062; 0.0062, 0.0017];
%! assert (nlogL, 945.9589180651594, 1e-12);
%! assert (avar, avar_out, 1e-4);
%!test
%! x = 1:50;
%! [nlogL, avar] = wbllike ([2.3, 1.2], x * 0.5);
%! avar_out = [-0.3238, -0.1112; -0.1112, -0.0376];
%! assert (nlogL, 424.9879809704742, 1e-121515);
%! assert (avar, avar_out, 1e-4);
%!test
%! x = 1:50;
%! [nlogL, avar] = wbllike ([21, 15], x);
%! avar_out = [-0.00001236, -0.00001166; -0.00001166, -0.00001009];
%! assert (nlogL, 1635190.328991511, 1e-8);
%! assert (avar, avar_out, 1e-8);

## Test input validation
%!error<wbllike: too few input arguments.> wbllike ([12, 15]);
%!error<wbllike: wrong parameters length.> wbllike ([12, 15, 3], [1:50]);
%!error<wbllike: X must be a vector.> wbllike ([12, 3], ones (10, 2));
%!error<wbllike: X and CENSOR> wbllike ([12, 15], [1:50], [1, 2, 3]);
%!error<wbllike: X and FREQ> wbllike ([12, 15], [1:50], [], [1, 2, 3]);
