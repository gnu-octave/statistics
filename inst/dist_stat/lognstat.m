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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} lognstat (@var{mu}, @var{sigma})
##
## Compute statistics of the log-normal distribution.
##
## @code{[@var{m}, @var{v}] = lognstat (@var{mu}, @var{sigma})} returns the mean
## and variance of the log-normal distribution with mean parameter @var{mu} and
## standard deviation parameter @var{sigma}, each corresponding to the
## associated normal distribution.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the log-normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-normal_distribution}
##
## @seealso{logncdf, logninv, lognpdf, lognrnd, lognfit, lognlike}
## @end deftypefn

function [m, v] = lognstat (mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("lognstat: function called with too few input arguments.");
  endif

  ## Check for MU and SIGMA being numeric
  if (! (isnumeric (mu) && isnumeric (sigma)))
    error ("lognstat: MU and SIGMA must be numeric.");
  endif

  ## Check for MU and SIGMA being real
  if (iscomplex (mu) || iscomplex (sigma))
    error ("lognstat: MU and SIGMA must not be complex.");
  endif

  ## Check for common size of MU and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, mu, sigma] = common_size (mu, sigma);
    if (retval > 0)
      error ("lognstat: MU and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = exp (mu + (sigma .^ 2) ./ 2);
  v = (exp (sigma .^ 2) - 1) .* exp (2 .* mu + sigma .^ 2);

  ## Continue argument check
  k = find (! (sigma >= 0) | ! (sigma < Inf));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<lognstat: function called with too few input arguments.> lognstat ()
%!error<lognstat: function called with too few input arguments.> lognstat (1)
%!error<lognstat: MU and SIGMA must be numeric.> lognstat ({}, 2)
%!error<lognstat: MU and SIGMA must be numeric.> lognstat (1, "")
%!error<lognstat: MU and SIGMA must not be complex.> lognstat (i, 2)
%!error<lognstat: MU and SIGMA must not be complex.> lognstat (1, i)
%!error<lognstat: MU and SIGMA must be of common size or scalars.> ...
%! lognstat (ones (3), ones (2))
%!error<lognstat: MU and SIGMA must be of common size or scalars.> ...
%! lognstat (ones (2), ones (3))

## Output validation tests
%!test
%! mu = 0:0.2:1;
%! sigma = 0.2:0.2:1.2;
%! [m, v] = lognstat (mu, sigma);
%! expected_m = [1.0202, 1.3231, 1.7860, 2.5093,  3.6693,   5.5845];
%! expected_v = [0.0425, 0.3038, 1.3823, 5.6447, 23.1345, 100.4437];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
%!test
%! sigma = 0.2:0.2:1.2;
%! [m, v] = lognstat (0, sigma);
%! expected_m = [1.0202, 1.0833, 1.1972, 1.3771, 1.6487,  2.0544];
%! expected_v = [0.0425, 0.2036, 0.6211, 1.7002, 4.6708, 13.5936];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
