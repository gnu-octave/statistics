## Copyright (C) 2021 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{nlogL} =} explike (@var{mu}, @var{x})
## @deftypefnx {statistics} {[@var{nlogL}, @var{avar}] =} explike (@var{mu}, @var{x})
## @deftypefnx {statistics} {[@dots{}] =} explike (@var{mu}, @var{x}, @var{censor})
## @deftypefnx {statistics} {[@dots{}] =} explike (@var{mu}, @var{x}, @var{censor}, @var{freq})
##
## Negative log-likelihood for the exponential distribution.
##
## @code{@var{nlogL} = explike (@var{mu}, @var{x})} returns the negative
## log likelihood of the data in @var{x} corresponding to the exponential
## distribution with mean parameter @var{mu}.  @var{x} must be a vector of
## non-negative values, otherwise @qcode{NaN} is returned.
##
## @code{[@var{nlogL}, @var{avar}] = explike (@var{mu}, @var{x})} also
## returns the inverse of Fisher's information matrix, @var{avar}.  If the input
## mean parameter, @var{mu}, is the maximum likelihood estimate, @var{avar} is
## its asymptotic variance.
##
## @code{[@dots{}] = explike (@var{mu}, @var{x}, @var{censor})} accepts a
## boolean vector, @var{censor}, of the same size as @var{x} with @qcode{1}s for
## observations that are right-censored and @qcode{0}s for observations that are
## observed exactly.  By default, or if left empty,
## @qcode{@var{censor} = zeros (size (@var{x}))}.
##
## @code{[@dots{}] = explike (@var{mu}, @var{x}, @var{censor}, @var{freq})}
## accepts a frequency vector, @var{freq}, of the same size as @var{x}.
## @var{freq} typically contains integer frequencies for the corresponding
## elements in @var{x}, but it can contain any non-integer non-negative values.
## By default, or if left empty, @qcode{@var{freq} = ones (size (@var{x}))}.
##
## A common alternative parameterization of the exponential distribution is to
## use the parameter @math{λ} defined as the mean number of events in an
## interval as opposed to the parameter @math{μ}, which is the mean wait time
## for an event to occur. @math{λ} and @math{μ} are reciprocals,
## i.e. @math{μ = 1 / λ}.
##
## Further information about the exponential distribution can be found at
## @url{https://en.wikipedia.org/wiki/Exponential_distribution}
##
## @seealso{expcdf, expinv, exppdf, exprnd, expfit, expstat}
## @end deftypefn

function [nlogL, avar] = explike (mu, x, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("explike: function called with too few input arguments.");
  endif

  if (! isvector (x))
    error ("explike: X must be a vector.");
  endif

  if (numel (mu) != 1)
    error ("explike: MU must be a scalar.");
  endif

  ## Return NaNs for non-positive MU or negative values in X
  if (mu <= 0 || any (x(:) < 0))
    nlogL = NaN;
    if (nargout > 1)
      avar = NaN;
    endif
    return
  endif

  if (nargin < 3 || isempty (censor))
    censor = zeros (size (x));
  elseif (! isequal (size (x), size (censor)))
    error ("explike: X and CENSOR vectors mismatch.");
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
    error ("explike: X and FREQ vectors mismatch.");
  endif

  ## Start processing
  numx = numel (x);
  sumz = sum (x .* freq) / mu;
  numc = numx - sum (freq .* censor);

  ## Calculate negative log likelihood
  nlogL = sumz + numc * log (mu);

  ## Optionally calculate the inverse (reciprocal) of the second derivative
  ## of the negative log likelihood with respect to parameter
  if (nargout > 1)
    avar = (mu ^ 2) ./ (2 * sumz - numc);
  endif

endfunction

%!test
%! x = 12;
%! beta = 5;
%! [L, V] = explike (beta, x);
%! expected_L = 4.0094;
%! expected_V = 6.5789;
%! assert (L, expected_L, 0.001);
%! assert (V, expected_V, 0.001);

%!test
%! x = 1:5;
%! beta = 2;
%! [L, V] = explike (beta, x);
%! expected_L = 10.9657;
%! expected_V = 0.4;
%! assert (L, expected_L, 0.001);
%! assert (V, expected_V, 0.001);

## Test input validation
%!error<explike: function called with too few input arguments.> explike ()
%!error<explike: function called with too few input arguments.> explike (2)
%!error<explike: MU must be a scalar.> explike ([12, 3], [1:50])
%!error<explike: X must be a vector.> explike (3, ones (10, 2))
%!error<explike: X and CENSOR vectors mismatch.> ...
%! explike (3, [1:50], [1, 2, 3])
%!error<explike: X and FREQ vectors mismatch.> ...
%! explike (3, [1:50], [], [1, 2, 3])
