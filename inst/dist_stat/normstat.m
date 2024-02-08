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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} normstat (@var{mu}, @var{sigma})
##
## Compute statistics of the normal distribution.
##
## @code{[@var{m}, @var{v}] = normstat (@var{mu}, @var{sigma})} returns the mean
## and variance of the normal distribution with non-centrality (distance)
## parameter @var{mu} and scale parameter @var{sigma}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Normal_distribution}
##
## @seealso{norminv, norminv, normpdf, normrnd, normfit, normlike}
## @end deftypefn

function [m, v] = normstat (mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("normstat: function called with too few input arguments.");
  endif

  ## Check for MU and SIGMA being numeric
  if (! (isnumeric (mu) && isnumeric (sigma)))
    error ("normstat: MU and SIGMA must be numeric.");
  endif

  ## Check for MU and SIGMA being real
  if (iscomplex (mu) || iscomplex (sigma))
    error ("normstat: MU and SIGMA must not be complex.");
  endif

  ## Check for common size of MU and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, mu, sigma] = common_size (mu, sigma);
    if (retval > 0)
      error ("normstat: MU and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = mu;
  v = sigma .* sigma;

  ## Continue argument check
  k = find (! (sigma > 0) | ! (sigma < Inf));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<normstat: function called with too few input arguments.> normstat ()
%!error<normstat: function called with too few input arguments.> normstat (1)
%!error<normstat: MU and SIGMA must be numeric.> normstat ({}, 2)
%!error<normstat: MU and SIGMA must be numeric.> normstat (1, "")
%!error<normstat: MU and SIGMA must not be complex.> normstat (i, 2)
%!error<normstat: MU and SIGMA must not be complex.> normstat (1, i)
%!error<normstat: MU and SIGMA must be of common size or scalars.> ...
%! normstat (ones (3), ones (2))
%!error<normstat: MU and SIGMA must be of common size or scalars.> ...
%! normstat (ones (2), ones (3))

## Output validation tests
%!test
%! mu = 1:6;
%! sigma = 0.2:0.2:1.2;
%! [m, v] = normstat (mu, sigma);
%! expected_v = [0.0400, 0.1600, 0.3600, 0.6400, 1.0000, 1.4400];
%! assert (m, mu);
%! assert (v, expected_v, 0.001);
%!test
%! sigma = 0.2:0.2:1.2;
%! [m, v] = normstat (0, sigma);
%! expected_mn = [0, 0, 0, 0, 0, 0];
%! expected_v = [0.0400, 0.1600, 0.3600, 0.6400, 1.0000, 1.4400];
%! assert (m, expected_mn, 0.001);
%! assert (v, expected_v, 0.001);
