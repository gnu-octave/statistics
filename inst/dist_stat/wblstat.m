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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} wblstat (@var{lambda}, @var{k})
##
## Compute statistics of the Weibull distribution.
##
## @code{[@var{m}, @var{v}] = wblstat (@var{lambda}, @var{k})} returns the mean
## and variance of the Weibull distribution with scale parameter @var{lambda}
## and shape parameter @var{k}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the Weibull distribution can be found at
## @url{https://en.wikipedia.org/wiki/Weibull_distribution}
##
## @seealso{wblcdf, wblinv, wblpdf, wblrnd, wblfit, wbllike, wblplot}
## @end deftypefn

function [m, v] = wblstat (lambda, k)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("wblstat: function called with too few input arguments.");
  endif

  ## Check for LAMBDA and K being numeric
  if (! (isnumeric (lambda) && isnumeric (k)))
    error ("wblstat: LAMBDA and K must be numeric.");
  endif

  ## Check for LAMBDA and K being real
  if (iscomplex (lambda) || iscomplex (k))
    error ("wblstat: LAMBDA and K must not be complex.");
  endif

  ## Check for common size of LAMBDA and K
  if (! isscalar (lambda) || ! isscalar (k))
    [retval, lambda, k] = common_size (lambda, k);
    if (retval > 0)
      error ("wblstat: LAMBDA and K must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = lambda .* gamma (1 + 1 ./ k);
  v = (lambda .^ 2) .* gamma (1 + 2 ./ k) - m .^ 2;

  ## Continue argument check
  is_nan = find (! (lambda > 0) | ! (lambda < Inf) | ! (k > 0) | ! (k < Inf));
  if (any (is_nan))
    m(is_nan) = NaN;
    v(is_nan) = NaN;
  endif

endfunction

## Input validation tests
%!error<wblstat: function called with too few input arguments.> wblstat ()
%!error<wblstat: function called with too few input arguments.> wblstat (1)
%!error<wblstat: LAMBDA and K must be numeric.> wblstat ({}, 2)
%!error<wblstat: LAMBDA and K must be numeric.> wblstat (1, "")
%!error<wblstat: LAMBDA and K must not be complex.> wblstat (i, 2)
%!error<wblstat: LAMBDA and K must not be complex.> wblstat (1, i)
%!error<wblstat: LAMBDA and K must be of common size or scalars.> ...
%! wblstat (ones (3), ones (2))
%!error<wblstat: LAMBDA and K must be of common size or scalars.> ...
%! wblstat (ones (2), ones (3))

## Output validation tests
%!test
%! lambda = 3:8;
%! k = 1:6;
%! [m, v] = wblstat (lambda, k);
%! expected_m = [3.0000, 3.5449, 4.4649, 5.4384, 6.4272, 7.4218];
%! expected_v = [9.0000, 3.4336, 2.6333, 2.3278, 2.1673, 2.0682];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
%!test
%! k = 1:6;
%! [m, v] = wblstat (6, k);
%! expected_m = [ 6.0000, 5.3174, 5.3579, 5.4384, 5.5090, 5.5663];
%! expected_v = [36.0000, 7.7257, 3.7920, 2.3278, 1.5923, 1.1634];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
