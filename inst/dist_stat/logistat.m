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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} logistat (@var{mu}, @var{sigma})
##
## Compute statistics of the logistic distribution.
##
## @code{[@var{m}, @var{v}] = logistat (@var{mu}, @var{sigma})} returns the mean
## and variance of the logistic distribution with mean parameter @var{mu} and
## scale parameter @var{sigma}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Logistic_distribution}
##
## @seealso{logicdf, logiinv, logipdf, logirnd, logifit, logilike}
## @end deftypefn

function [m, v] = logistat (mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("logistat: function called with too few input arguments.");
  endif

  ## Check for MU and SIGMA being numeric
  if (! (isnumeric (mu) && isnumeric (sigma)))
    error ("logistat: MU and SIGMA must be numeric.");
  endif

  ## Check for MU and SIGMA being real
  if (iscomplex (mu) || iscomplex (sigma))
    error ("logistat: MU and SIGMA must not be complex.");
  endif

  ## Check for common size of MU and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, mu, sigma] = common_size (mu, sigma);
    if (retval > 0)
      error ("logistat: MU and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = mu;
  v = (sigma .^ 2 .* pi .^ 2) ./ 3;

  ## Continue argument check
  m(sigma <= 0) = NaN;
  v(sigma <= 0) = NaN;

endfunction

## Input validation tests
%!error<logistat: function called with too few input arguments.> logistat ()
%!error<logistat: function called with too few input arguments.> logistat (1)
%!error<logistat: MU and SIGMA must be numeric.> logistat ({}, 2)
%!error<logistat: MU and SIGMA must be numeric.> logistat (1, "")
%!error<logistat: MU and SIGMA must not be complex.> logistat (i, 2)
%!error<logistat: MU and SIGMA must not be complex.> logistat (1, i)
%!error<logistat: MU and SIGMA must be of common size or scalars.> ...
%! logistat (ones (3), ones (2))
%!error<logistat: MU and SIGMA must be of common size or scalars.> ...
%! logistat (ones (2), ones (3))

## Output validation tests
%!test
%! [m, v] = logistat (0, 1);
%! assert (m, 0);
%! assert (v, 3.2899, 0.001);
%!test
%! [m, v] = logistat (0, 0.8);
%! assert (m, 0);
%! assert (v, 2.1055, 0.001);
%!test
%! [m, v] = logistat (1, 0.6);
%! assert (m, 1);
%! assert (v, 1.1844, 0.001);
%!test
%! [m, v] = logistat (0, 0.4);
%! assert (m, 0);
%! assert (v, 0.5264, 0.001);
%!test
%! [m, v] = logistat (-1, 0.2);
%! assert (m, -1);
%! assert (v, 0.1316, 0.001);
