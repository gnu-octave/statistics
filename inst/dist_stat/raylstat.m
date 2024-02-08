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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} raylstat (@var{sigma})
##
## Compute statistics of the Rayleigh distribution.
##
## @code{[@var{m}, @var{v}] = raylstat (@var{sigma})} returns the mean and
## variance of the Rayleigh distribution with scale parameter @var{sigma}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the same size of the
## input argument.
##
## Further information about the Rayleigh distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rayleigh_distribution}
##
## @seealso{raylcdf, raylinv, raylpdf, raylrnd, raylfit, rayllike}
## @end deftypefn

function [m, v] = raylstat (sigma)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("raylstat: function called with too few input arguments.");
  endif

  ## Check for SIGMA being numeric
  if (! isnumeric (sigma))
    error ("raylstat: SIGMA must be numeric.");
  endif

  ## Check for SIGMA being real
  if (iscomplex (sigma))
    error ("raylstat: SIGMA must not be complex.");
  endif

  ## Calculate moments
  m = sigma .* sqrt (pi ./ 2);
  v = (2 - pi ./ 2) .* sigma .^ 2;

  ## Continue argument check
  k = find (! (sigma > 0));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<raylstat: function called with too few input arguments.> raylstat ()
%!error<raylstat: SIGMA must be numeric.> raylstat ({})
%!error<raylstat: SIGMA must be numeric.> raylstat ("")
%!error<raylstat: SIGMA must not be complex.> raylstat (i)

## Output validation tests
%!test
%! sigma = 1:6;
%! [m, v] = raylstat (sigma);
%! expected_m = [1.2533, 2.5066, 3.7599, 5.0133, 6.2666, 7.5199];
%! expected_v = [0.4292, 1.7168, 3.8628, 6.8673, 10.7301, 15.4513];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
