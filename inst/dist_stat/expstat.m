## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} expstat (@var{mu})
##
## Compute statistics of the exponential distribution.
##
## @code{[@var{m}, @var{v}] = expstat (@var{mu})} returns the mean and
## variance of the exponential distribution with mean parameter @var{mu}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the same size of the
## input argument.
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
## @seealso{expcdf, expinv, exppdf, exprnd, expfit, explike}
## @end deftypefn

function [m, v] = expstat (mu)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("expstat: function called with too few input arguments.");
  endif

  ## Check for MU being numeric
  if (! isnumeric (mu))
    error ("expstat: MU must be numeric.");
  endif

  ## Check for MU being real
  if (iscomplex (mu))
    error ("expstat: MU must not be complex.");
  endif

  ## Calculate moments
  m = mu;
  v = m .^ 2;

  ## Continue argument check
  k = find (! (mu > 0) | ! (mu < Inf));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<expstat: function called with too few input arguments.> expstat ()
%!error<expstat: MU must be numeric.> expstat ({})
%!error<expstat: MU must be numeric.> expstat ("")
%!error<expstat: MU must not be complex.> expstat (i)

## Output validation tests
%!test
%! mu = 1:6;
%! [m, v] = expstat (mu);
%! assert (m, [1, 2, 3, 4, 5, 6], 0.001);
%! assert (v, [1, 4, 9, 16, 25, 36], 0.001);
