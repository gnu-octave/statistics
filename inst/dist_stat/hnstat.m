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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} hnstat (@var{mu}, @var{sigma})
##
## Compute statistics of the half-normal distribution.
##
## @code{[@var{m}, @var{v}] = hnstat (@var{mu}, @var{sigma})} returns the mean
## and variance of the half-normal distribution with non-centrality (distance)
## parameter @var{mu} and scale parameter @var{sigma}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the half-normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Half-normal_distribution}
##
## @seealso{norminv, norminv, normpdf, normrnd, normfit, normlike}
## @end deftypefn

function [m, v] = hnstat (mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("hnstat: function called with too few input arguments.");
  endif

  ## Check for MU and SIGMA being numeric
  if (! (isnumeric (mu) && isnumeric (sigma)))
    error ("hnstat: MU and SIGMA must be numeric.");
  endif

  ## Check for MU and SIGMA being real
  if (iscomplex (mu) || iscomplex (sigma))
    error ("hnstat: MU and SIGMA must not be complex.");
  endif

  ## Check for common size of MU and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, mu, sigma] = common_size (mu, sigma);
    if (retval > 0)
      error ("hnstat: MU and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = mu + (sigma .* sqrt (2)) ./ sqrt (pi);
  v = sigma .^ 2 .* (1 - 2 / pi);

  ## Continue argument check
  k = find (! (sigma > 0) | ! (sigma < Inf));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<hnstat: function called with too few input arguments.> hnstat ()
%!error<hnstat: function called with too few input arguments.> hnstat (1)
%!error<hnstat: MU and SIGMA must be numeric.> hnstat ({}, 2)
%!error<hnstat: MU and SIGMA must be numeric.> hnstat (1, "")
%!error<hnstat: MU and SIGMA must not be complex.> hnstat (i, 2)
%!error<hnstat: MU and SIGMA must not be complex.> hnstat (1, i)
%!error<hnstat: MU and SIGMA must be of common size or scalars.> ...
%! hnstat (ones (3), ones (2))
%!error<hnstat: MU and SIGMA must be of common size or scalars.> ...
%! hnstat (ones (2), ones (3))

## Output validation tests
%!test
%! [m, v] = hnstat (0, 1);
%! assert (m, 0.7979, 1e-4);
%! assert (v, 0.3634, 1e-4);
%!test
%! [m, v] = hnstat (2, 1);
%! assert (m, 2.7979, 1e-4);
%! assert (v, 0.3634, 1e-4);
%!test
%! [m, v] = hnstat (2, 2);
%! assert (m, 3.5958, 1e-4);
%! assert (v, 1.4535, 1e-4);
%!test
%! [m, v] = hnstat (2, 2.5);
%! assert (m, 3.9947, 1e-4);
%! assert (v, 2.2711, 1e-4);
%!test
%! [m, v] = hnstat (1.5, 0.5);
%! assert (m, 1.8989, 1e-4);
%! assert (v, 0.0908, 1e-4);
%!test
%! [m, v] = hnstat (-1.5, 0.5);
%! assert (m, -1.1011, 1e-4);
%! assert (v, 0.0908, 1e-4);
