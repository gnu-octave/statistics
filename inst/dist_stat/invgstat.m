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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} invgstat (@var{mu}, @var{lambda})
##
## Compute statistics of the inverse Gaussian distribution.
##
## @code{[@var{m}, @var{v}] = invgstat (@var{mu}, @var{lambda})} returns the
## mean and variance of the inverse Gaussian distribution with mean parameter
## @var{mu} and shape parameter @var{lambda}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the inverse Gaussian distribution can be found at
## @url{https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}
##
## @seealso{invgcdf, invginv, invgpdf, invgrnd, invgfit, invglike}
## @end deftypefn

function [m, v] = invgstat (mu, lambda)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("invgstat: function called with too few input arguments.");
  endif

  ## Check for MU and LAMBDA being numeric
  if (! (isnumeric (mu) && isnumeric (lambda)))
    error ("invgstat: MU and LAMBDA must be numeric.");
  endif

  ## Check for MU and LAMBDA being real
  if (iscomplex (mu) || iscomplex (lambda))
    error ("invgstat: MU and LAMBDA must not be complex.");
  endif

  ## Check for common size of MU and LAMBDA
  if (! isscalar (mu) || ! isscalar (lambda))
    [retval, mu, lambda] = common_size (mu, lambda);
    if (retval > 0)
      error ("invgstat: MU and LAMBDA must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = mu;
  v = (mu .^ 3) ./ lambda;

  ## Continue argument check
  m(lambda <= 0 | mu <= 0) = NaN;
  v(lambda <= 0 | mu <= 0) = NaN;

endfunction

## Input validation tests
%!error<invgstat: function called with too few input arguments.> invgstat ()
%!error<invgstat: function called with too few input arguments.> invgstat (1)
%!error<invgstat: MU and LAMBDA must be numeric.> invgstat ({}, 2)
%!error<invgstat: MU and LAMBDA must be numeric.> invgstat (1, "")
%!error<invgstat: MU and LAMBDA must not be complex.> invgstat (i, 2)
%!error<invgstat: MU and LAMBDA must not be complex.> invgstat (1, i)
%!error<invgstat: MU and LAMBDA must be of common size or scalars.> ...
%! invgstat (ones (3), ones (2))
%!error<invgstat: MU and LAMBDA must be of common size or scalars.> ...
%! invgstat (ones (2), ones (3))

## Output validation tests
%!test
%! [m, v] = invgstat (1, 1);
%! assert (m, 1);
%! assert (v, 1);
%!test
%! [m, v] = invgstat (2, 1);
%! assert (m, 2);
%! assert (v, 8);
%!test
%! [m, v] = invgstat (2, 2);
%! assert (m, 2);
%! assert (v, 4);
%!test
%! [m, v] = invgstat (2, 2.5);
%! assert (m, 2);
%! assert (v, 3.2);
%!test
%! [m, v] = invgstat (1.5, 0.5);
%! assert (m, 1.5);
%! assert (v, 6.75);
