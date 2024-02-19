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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} nakastat (@var{mu}, @var{omega})
##
## Compute statistics of the Nakagami distribution.
##
## @code{[@var{m}, @var{v}] = nakastat (@var{mu}, @var{omega})} returns the mean
## and variance of the Nakagami distribution with shape parameter @var{mu} and
## spread parameter @var{omega}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the Nakagami distribution can be found at
## @url{https://en.wikipedia.org/wiki/Normal_distribution}
##
## @seealso{nakacdf, nakainv, nakapdf, nakarnd, nakafit, nakalike}
## @end deftypefn

function [m, v] = nakastat (mu, omega)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("nakastat: function called with too few input arguments.");
  endif

  ## Check for MU and OMEGA being numeric
  if (! (isnumeric (mu) && isnumeric (omega)))
    error ("nakastat: MU and OMEGA must be numeric.");
  endif

  ## Check for MU and OMEGA being real
  if (iscomplex (mu) || iscomplex (omega))
    error ("nakastat: MU and OMEGA must not be complex.");
  endif

  ## Check for common size of MU and OMEGA
  if (! isscalar (mu) || ! isscalar (omega))
    [retval, mu, omega] = common_size (mu, omega);
    if (retval > 0)
      error ("nakastat: MU and OMEGA must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  g = gamma (mu + 0.5) ./ gamma (mu);
  m = g .* sqrt (omega ./ mu);
  v = omega .* (1 - ((1 ./ mu) .* (g .^ 2)));

  ## Continue argument check
  knan = mu < 0.5 | omega <= 0;
  m(knan) = NaN;
  v(knan) = NaN;

endfunction

## Input validation tests
%!error<nakastat: function called with too few input arguments.> nakastat ()
%!error<nakastat: function called with too few input arguments.> nakastat (1)
%!error<nakastat: MU and OMEGA must be numeric.> nakastat ({}, 2)
%!error<nakastat: MU and OMEGA must be numeric.> nakastat (1, "")
%!error<nakastat: MU and OMEGA must not be complex.> nakastat (i, 2)
%!error<nakastat: MU and OMEGA must not be complex.> nakastat (1, i)
%!error<nakastat: MU and OMEGA must be of common size or scalars.> ...
%! nakastat (ones (3), ones (2))
%!error<nakastat: MU and OMEGA must be of common size or scalars.> ...
%! nakastat (ones (2), ones (3))

## Output validation tests
%!test
%! [m, v] = nakastat (1, 1);
%! assert (m, 0.8862269254, 1e-10);
%! assert (v, 0.2146018366, 1e-10);
%!test
%! [m, v] = nakastat (1, 2);
%! assert (m, 1.25331413731, 1e-10);
%! assert (v, 0.42920367321, 1e-10);
%!test
%! [m, v] = nakastat (2, 1);
%! assert (m, 0.93998560299, 1e-10);
%! assert (v, 0.11642706618, 1e-10);
