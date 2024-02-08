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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} poisstat (@var{lambda})
##
## Compute statistics of the Poisson distribution.
##
## @code{[@var{m}, @var{v}] = poisstat (@var{lambda})} returns the mean and
## variance of the Poisson distribution with rate parameter @var{lambda}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the same size of the
## input argument.
##
## Further information about the Poisson distribution can be found at
## @url{https://en.wikipedia.org/wiki/Poisson_distribution}
##
## @seealso{poisscdf, poissinv, poisspdf, poissrnd, poissfit, poisslike}
## @end deftypefn

function [m, v] = poisstat (lambda)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("poisstat: function called with too few input arguments.");
  endif

  ## Check for SIGMA being numeric
  if (! isnumeric (lambda))
    error ("poisstat: SIGMA must be numeric.");
  endif

  ## Check for SIGMA being real
  if (iscomplex (lambda))
    error ("poisstat: SIGMA must not be complex.");
  endif

  ## Set moments
  m = lambda;
  v = lambda;

  ## Continue argument check
  k = find (! (lambda > 0) | ! (lambda < Inf));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<poisstat: function called with too few input arguments.> poisstat ()
%!error<poisstat: SIGMA must be numeric.> poisstat ({})
%!error<poisstat: SIGMA must be numeric.> poisstat ("")
%!error<poisstat: SIGMA must not be complex.> poisstat (i)

## Output validation tests
%!test
%! lambda = 1 ./ (1:6);
%! [m, v] = poisstat (lambda);
%! assert (m, lambda);
%! assert (v, lambda);
