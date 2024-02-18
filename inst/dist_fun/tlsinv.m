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
## @deftypefn {statistics} {@var{x} =} tlsinv (@var{p}, @var{mu}, @var{sigma}, @var{nu})
##
## Inverse of the location-scale Student's T cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the location-scale Student's T distribution with location parameter @var{mu},
## scale parameter @var{sigma}, and @var{nu} degrees of freedom.  The size of
## @var{x} is the common size of @var{p}, @var{mu}, @var{sigma}, and @var{nu}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## Further information about the location-scale Student's T distribution can be
## found at @url{https://en.wikipedia.org/wiki/Student%27s_t-distribution#Location-scale_t_distribution}
##
## @seealso{tlscdf, tlspdf, tlsrnd, tlsfit, tlslike, tlsstat}
## @end deftypefn

function x = tlsinv (p, mu, sigma, nu)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("tlsinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, MU, SIGMA, and NU
  if (! isscalar (p) || ! isscalar (mu) || ! isscalar (sigma) || ! isscalar (nu))
    [retval, p, mu, sigma, nu] = common_size (p, mu, sigma, nu);
    if (retval > 0)
      error ("tlsinv: P, MU, SIGMA, and NU must be of common size or scalars.");
    endif
  endif

  ## Check for P, MU, SIGMA, and NU being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (sigma) || iscomplex (nu))
    error ("tlsinv: P, MU, SIGMA, and NU must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (mu, "single") ||
      isa (sigma, "single") || isa (nu, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Force invalid SIGMA parameter to NaN
  sigma(sigma <= 0) = NaN;

  ## Call tinv to do the work
  x = tinv (p, nu) .* sigma + mu;

  ## Force class type
  x = cast (x, cls);

endfunction

%!demo
%! ## Plot various iCDFs from the location-scale Student's T distribution
%! p = 0.001:0.001:0.999;
%! x1 = tlsinv (p, 0, 1, 1);
%! x2 = tlsinv (p, 0, 2, 2);
%! x3 = tlsinv (p, 3, 2, 5);
%! x4 = tlsinv (p, -1, 3, Inf);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-m")
%! grid on
%! xlim ([0, 1])
%! ylim ([-8, 8])
%! legend ({"mu = 0, sigma = 1, nu = 1", "mu = 0, sigma = 2, nu = 2", ...
%!          "mu = 3, sigma = 2, nu = 5", 'mu = -1, sigma = 3, nu = \infty'}, ...
%!         "location", "southeast")
%! title ("Location-scale Student's T iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (tlsinv (p, 0, 1, ones (1,5)), [NaN -Inf 0 Inf NaN])
%!assert (tlsinv (p, 0, 1, 1), [NaN -Inf 0 Inf NaN], eps)
%!assert (tlsinv (p, 0, 1, [1 0 NaN 1 1]), [NaN NaN NaN Inf NaN], eps)
%!assert (tlsinv ([p(1:2) NaN p(4:5)], 0, 1, 1), [NaN -Inf NaN Inf NaN])

## Test class of input preserved
%!assert (class (tlsinv ([p, NaN], 0, 1, 1)), "double")
%!assert (class (tlsinv (single ([p, NaN]), 0, 1, 1)), "single")
%!assert (class (tlsinv ([p, NaN], single (0), 1, 1)), "single")
%!assert (class (tlsinv ([p, NaN], 0, single (1), 1)), "single")
%!assert (class (tlsinv ([p, NaN], 0, 1, single (1))), "single")

## Test input validation
%!error<tlsinv: function called with too few input arguments.> tlsinv ()
%!error<tlsinv: function called with too few input arguments.> tlsinv (1)
%!error<tlsinv: function called with too few input arguments.> tlsinv (1, 2)
%!error<tlsinv: function called with too few input arguments.> tlsinv (1, 2, 3)
%!error<tlsinv: P, MU, SIGMA, and NU must be of common size or scalars.> ...
%! tlsinv (ones (3), ones (2), 1, 1)
%!error<tlsinv: P, MU, SIGMA, and NU must be of common size or scalars.> ...
%! tlsinv (ones (2), 1, ones (3), 1)
%!error<tlsinv: P, MU, SIGMA, and NU must be of common size or scalars.> ...
%! tlsinv (ones (2), 1, 1, ones (3))
%!error<tlsinv: P, MU, SIGMA, and NU must not be complex.> tlsinv (i, 2, 3, 4)
%!error<tlsinv: P, MU, SIGMA, and NU must not be complex.> tlsinv (2, i, 3, 4)
%!error<tlsinv: P, MU, SIGMA, and NU must not be complex.> tlsinv (2, 2, i, 4)
%!error<tlsinv: P, MU, SIGMA, and NU must not be complex.> tlsinv (2, 2, 3, i)
