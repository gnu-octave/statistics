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
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{mu} is a scalar containing the mean of the exponential distribution
## (@math{@var{mu} = 1 / λ}, where λ is the rate, or inverse scale parameter).
## @item
## @var{x} is the vector of given values.
## @item
## @var{censor} is a boolean vector of the same size as @var{x} with 1 for
## observations that are right-censored and 0 for observations that are observed
## exactly.
## @item
## @var{freq} is a vector of the same size as @var{x} that contains integer
## frequencies for the corresponding elements in @var{x}, but may contain any
## non-integer non-negative values.  Pass in [] for @var{censor} to use its
## default value.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{nlogL} is the negative log-likelihood.
## @item
## @var{avar} is the inverse of the Fisher information matrix.
## (The Fisher information matrix is the second derivative of the negative log
## likelihood with respect to the parameter value.)
## @end itemize
##
## @seealso{expcdf, expinv, exppdf, exprnd, expfit, expstat}
## @end deftypefn

function [nlogL, avar] = explike (mu, x, censor, freq)

  ## Check input arguments
  if (nargin < 2)
    error ("explike: too few input arguments.");
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
%!error explike ();
