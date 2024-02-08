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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} tlsstat (@var{mu}, @var{sigma}, @var{df})
##
## Compute statistics of the location-scale Student's T distribution.
##
## @code{[@var{m}, @var{v}] = tlsstat (@var{mu}, @var{sigma}, @var{df})} returns
## the mean and variance of the location-scale Student's T distribution with
## location parameter @var{mu}, scale parameter @var{sigma}, and @var{df}
## degrees of freedom.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the location-scale Student's T distribution can be
## found at @url{https://en.wikipedia.org/wiki/Student%27s_t-distribution#Location-scale_t_distribution}
##
## @seealso{tlscdf, tlsinv, tlspdf, tlsrnd, tlsfit, tlslike}
## @end deftypefn

function [m, v] = tlsstat (mu, sigma, df)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("tlsstat: function called with too few input arguments.");
  endif

  ## Check for MU, SIGMA, and DF being numeric
  if (! (isnumeric (mu) && isnumeric (sigma) && isnumeric (df)))
    error ("tlsstat: MU, SIGMA, and DF must be numeric.");
  endif

  ## Check for MU, SIGMA, and DF being real
  if (iscomplex (mu) || iscomplex (sigma) || iscomplex (df))
    error ("tlsstat: MU, SIGMA, and DF must not be complex.");
  endif

  ## Check for common size of MU, SIGMA, and DF
  if (! isscalar (mu) || ! isscalar (sigma) || ! isscalar (df))
    [retval, mu, sigma, df] = common_size (mu, sigma, df);
    if (retval > 0)
      error ("tlsstat: MU, SIGMA, and DF must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = zeros (size (df)) + mu;
  v = sigma .* (df ./ (df - 2));

  ## Continue argument check
  k = find (! (df > 1) | ! (df < Inf));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif
  k = find (! (df > 2) & (df < Inf));
  if (any (k))
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<tlsstat: function called with too few input arguments.> tlsstat ()
%!error<tlsstat: function called with too few input arguments.> tlsstat (1)
%!error<tlsstat: function called with too few input arguments.> tlsstat (1, 2)
%!error<tlsstat: MU, SIGMA, and DF must be numeric.> tlsstat ({}, 2, 3)
%!error<tlsstat: MU, SIGMA, and DF must be numeric.> tlsstat (1, "", 3)
%!error<tlsstat: MU, SIGMA, and DF must be numeric.> tlsstat (1, 2, ["d"])
%!error<tlsstat: MU, SIGMA, and DF must not be complex.> tlsstat (i, 2, 3)
%!error<tlsstat: MU, SIGMA, and DF must not be complex.> tlsstat (1, i, 3)
%!error<tlsstat: MU, SIGMA, and DF must not be complex.> tlsstat (1, 2, i)
%!error<tlsstat: MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlsstat (ones (3), ones (2), 1)
%!error<tlsstat: MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlsstat (ones (2), 1, ones (3))
%!error<tlsstat: MU, SIGMA, and DF must be of common size or scalars.> ...
%! tlsstat (1, ones (2), ones (3))

## Output validation tests
%!test
%! [m, v] = tlsstat (0, 1, 0);
%! assert (m, NaN);
%! assert (v, NaN);
%!test
%! [m, v] = tlsstat (0, 1, 1);
%! assert (m, NaN);
%! assert (v, NaN);
%!test
%! [m, v] = tlsstat (2, 1, 1);
%! assert (m, NaN);
%! assert (v, NaN);
%!test
%! [m, v] = tlsstat (-2, 1, 1);
%! assert (m, NaN);
%! assert (v, NaN);
%!test
%! [m, v] = tlsstat (0, 1, 2);
%! assert (m, 0);
%! assert (v, NaN);
%!test
%! [m, v] = tlsstat (2, 1, 2);
%! assert (m, 2);
%! assert (v, NaN);
%!test
%! [m, v] = tlsstat (-2, 1, 2);
%! assert (m, -2);
%! assert (v, NaN);
%!test
%! [m, v] = tlsstat (0, 2, 2);
%! assert (m, 0);
%! assert (v, NaN);
%!test
%! [m, v] = tlsstat (2, 2, 2);
%! assert (m, 2);
%! assert (v, NaN);
%!test
%! [m, v] = tlsstat (-2, 2, 2);
%! assert (m, -2);
%! assert (v, NaN);
%!test
%! [m, v] = tlsstat (0, 1, 3);
%! assert (m, 0);
%! assert (v, 3);
%!test
%! [m, v] = tlsstat (0, 2, 3);
%! assert (m, 0);
%! assert (v, 6);
%!test
%! [m, v] = tlsstat (2, 1, 3);
%! assert (m, 2);
%! assert (v, 3);
%!test
%! [m, v] = tlsstat (2, 2, 3);
%! assert (m, 2);
%! assert (v, 6);
%!test
%! [m, v] = tlsstat (-2, 1, 3);
%! assert (m, -2);
%! assert (v, 3);
%!test
%! [m, v] = tlsstat (-2, 2, 3);
%! assert (m, -2);
%! assert (v, 6);
