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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} gamstat (@var{a}, @var{b})
##
## Compute statistics of the Gamma distribution.
##
## @code{[@var{m}, @var{v}] = gamstat (@var{a}, @var{b})} returns the mean
## and variance of the Gamma distribution with with shape parameter @var{a} and
## scale parameter @var{b}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## OCTAVE/MATLAB use the alternative parameterization given by the pair
## @math{α, β}, i.e. shape @var{a} and scale @var{b}.  In Wikipedia, the two
## common parameterizations use the pairs @math{k, θ}, as shape and scale, and
## @math{α, β}, as shape and rate, respectively.  The parameter names @var{a}
## and @var{b} used here (for MATLAB compatibility) correspond to the parameter
## notation @math{k, θ} instead of the @math{α, β} as reported in Wikipedia.
##
## Further information about the Gamma distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gamma_distribution}
##
## @seealso{gamcdf, gaminv, gampdf, gamrnd, gamfit, gamlike}
## @end deftypefn

function [m, v] = gamstat (a, b)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("gamstat: function called with too few input arguments.");
  endif

  ## Check for A and B being numeric
  if (! (isnumeric (a) && isnumeric (b)))
    error ("gamstat: A and B must be numeric.");
  endif

  ## Check for A and B being real
  if (iscomplex (a) || iscomplex (b))
    error ("gamstat: A and B must not be complex.");
  endif

  ## Check for common size of A and B
  if (! isscalar (a) || ! isscalar (b))
    [retval, a, b] = common_size (a, b);
    if (retval > 0)
      error ("gamstat: A and B must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = a .* b;
  v = a .* (b .^ 2);

  ## Continue argument check
  a = find (! (a > 0) | ! (a < Inf) | ! (b > 0) | ! (b < Inf));
  if (any (a))
    m(a) = NaN;
    v(a) = NaN;
  endif

endfunction

## Input validation tests
%!error<gamstat: function called with too few input arguments.> gamstat ()
%!error<gamstat: function called with too few input arguments.> gamstat (1)
%!error<gamstat: A and B must be numeric.> gamstat ({}, 2)
%!error<gamstat: A and B must be numeric.> gamstat (1, "")
%!error<gamstat: A and B must not be complex.> gamstat (i, 2)
%!error<gamstat: A and B must not be complex.> gamstat (1, i)
%!error<gamstat: A and B must be of common size or scalars.> ...
%! gamstat (ones (3), ones (2))
%!error<gamstat: A and B must be of common size or scalars.> ...
%! gamstat (ones (2), ones (3))

## Output validation tests
%!test
%! a = 1:6;
%! b = 1:0.2:2;
%! [m, v] = gamstat (a, b);
%! expected_m = [1.00, 2.40, 4.20,  6.40,  9.00, 12.00];
%! expected_v = [1.00, 2.88, 5.88, 10.24, 16.20, 24.00];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
%!test
%! a = 1:6;
%! [m, v] = gamstat (a, 1.5);
%! expected_m = [1.50, 3.00, 4.50, 6.00,  7.50,  9.00];
%! expected_v = [2.25, 4.50, 6.75, 9.00, 11.25, 13.50];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
