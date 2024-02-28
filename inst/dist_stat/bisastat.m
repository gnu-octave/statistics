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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} bisastat (@var{beta}, @var{gamma})
##
## Compute statistics of the Birnbaum-Saunders distribution.
##
## @code{[@var{m}, @var{v}] = bisastat (@var{beta}, @var{gamma})} returns the
## mean and variance of the Birnbaum-Saunders distribution with scale parameter
## @var{beta} and shape parameter @var{gamma}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the Birnbaum-Saunders distribution can be found at
## @url{https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution}
##
## @seealso{bisacdf, bisainv, bisapdf, bisarnd, bisafit, bisalike}
## @end deftypefn

function [m, v] = bisastat (beta, gamma)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("bisastat: function called with too few input arguments.");
  endif

  ## Check for BETA and GAMMA being numeric
  if (! (isnumeric (beta) && isnumeric (gamma)))
    error ("bisastat: BETA and GAMMA must be numeric.");
  endif

  ## Check for BETA and GAMMA being real
  if (iscomplex (beta) || iscomplex (gamma))
    error ("bisastat: BETA and GAMMA must not be complex.");
  endif

  ## Check for common size of BETA and GAMMA
  if (! isscalar (beta) || ! isscalar (gamma))
    [retval, beta, gamma] = common_size (beta, gamma);
    if (retval > 0)
      error ("bisastat: BETA and GAMMA must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = beta .* (1 + ((gamma .^ 2) ./ 2));
  v = ((beta .* gamma) .^ 2) .* (1 + ((5 .* (gamma .^ 2)) ./ 4));

  ## Continue argument check
  beta = find (! (beta > 0) | ! (beta < Inf) | ! (gamma > 0) | ! (gamma < Inf));
  if (any (beta))
    m(beta) = NaN;
    v(beta) = NaN;
  endif

endfunction

## Input validation tests
%!error<bisastat: function called with too few input arguments.> bisastat ()
%!error<bisastat: function called with too few input arguments.> bisastat (1)
%!error<bisastat: BETA and GAMMA must be numeric.> bisastat ({}, 2)
%!error<bisastat: BETA and GAMMA must be numeric.> bisastat (1, "")
%!error<bisastat: BETA and GAMMA must not be complex.> bisastat (i, 2)
%!error<bisastat: BETA and GAMMA must not be complex.> bisastat (1, i)
%!error<bisastat: BETA and GAMMA must be of common size or scalars.> ...
%! bisastat (ones (3), ones (2))
%!error<bisastat: BETA and GAMMA must be of common size or scalars.> ...
%! bisastat (ones (2), ones (3))

## Output validation tests
%!test
%! beta = 1:6;
%! gamma = 1:0.2:2;
%! [m, v] = bisastat (beta, gamma);
%! expected_m = [1.50, 3.44, 5.94,  9.12,  13.10, 18];
%! expected_v = [2.25, 16.128, 60.858, 172.032, 409.050, 864];
%! assert (m, expected_m, 1e-2);
%! assert (v, expected_v, 1e-3);
%!test
%! beta = 1:6;
%! [m, v] = bisastat (beta, 1.5);
%! expected_m = [2.125, 4.25, 6.375, 8.5, 10.625, 12.75];
%! expected_v = [8.5781, 34.3125, 77.2031, 137.2500, 214.4531, 308.8125];
%! assert (m, expected_m, 1e-3);
%! assert (v, expected_v, 1e-4);
