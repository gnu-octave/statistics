## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{x} =} nakacdf (@var{x}, @var{mu}, @var{omega})
##
## Inverse of the Nakagami cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the Nakagami distribution with shape parameter @var{mu} and spread parameter
## @var{omega}.  The size of @var{p} is the common size of @var{x}, @var{mu},
## and @var{omega}.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## Both parameters must be positive reals and @qcode{@var{mu} >= 0.5}.  For
## @qcode{@var{mu} < 0.5} or @qcode{@var{omega} <= 0}, @qcode{NaN} is returned.
##
## Further information about the Nakagami distribution can be found at
## @url{https://en.wikipedia.org/wiki/Nakagami_distribution}
##
## @seealso{nakacdf, nakapdf, nakarnd, nakafit, nakalike}
## @end deftypefn

function x = nakainv (p, mu, omega)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("nakainv: function called with too few input arguments.");
  endif

  ## Check for common size of P, MU, and OMEGA
  if (! isscalar (p) || ! isscalar (mu) || ! isscalar (omega))
    [retval, p, mu, omega] = common_size (p, mu, omega);
    if (retval > 0)
      error ("nakainv: P, MU, and OMEGA must be of common size or scalars.");
    endif
  endif

  ## Check for P, MU, and OMEGA being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (omega))
    error ("nakainv: P, MU, and OMEGA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (mu, "single") || isa (omega, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  ## Force invalid parameters and missing data to NaN
  k = isnan (p) | ! (p >= 0) | ! (p <= 1) | ! (mu >= 0.5) | ! (omega > 0);
  x(k) = NaN;

  ## Handle edge cases
  k = (p == 1) & (mu >= 0.5) & (mu < Inf) & (omega > 0) & (omega < Inf);
  x(k) = Inf;

  ## Find normal cases
  k = (0 < p) & (p < 1) & (0.5 <= mu) & (mu < Inf) ...
                        & (0 < omega) & (omega < Inf);

  ## Compute Nakagami iCDF
  if (isscalar (mu) && isscalar(omega))
    m_gamma = mu;
    w_gamma = omega / mu;
    x(k) = gaminv (p(k), m_gamma, w_gamma);
    x(k) = sqrt (x(k));
  else
    m_gamma = mu;
    w_gamma = omega ./ mu;
    x(k) = gaminv (p(k), m_gamma(k), w_gamma(k));
    x(k) = sqrt (x(k));
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the Nakagami distribution
%! p = 0.001:0.001:0.999;
%! x1 = nakainv (p, 0.5, 1);
%! x2 = nakainv (p, 1, 1);
%! x3 = nakainv (p, 1, 2);
%! x4 = nakainv (p, 1, 3);
%! x5 = nakainv (p, 2, 1);
%! x6 = nakainv (p, 2, 2);
%! x7 = nakainv (p, 5, 1);
%! plot (p, x1, "-r", p, x2, "-g", p, x3, "-y", p, x4, "-m", ...
%!       p, x5, "-k", p, x6, "-b", p, x7, "-c")
%! grid on
%! ylim ([0, 3])
%! legend ({"μ = 0.5, ω = 1", "μ = 1, ω = 1", "μ = 1, ω = 2", ...
%!          "μ = 1, ω = 3", "μ = 2, ω = 1", "μ = 2, ω = 2", ...
%!          "μ = 5, ω = 1"}, "location", "northwest")
%! title ("Nakagami iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p, y
%! p = [-Inf, -1, 0, 1/2, 1, 2, Inf];
%! y = [NaN, NaN, 0, 0.83255461115769769, Inf, NaN, NaN];
%!assert (nakainv (p, ones (1,7), ones (1,7)), y, eps)
%!assert (nakainv (p, 1, 1), y, eps)
%!assert (nakainv (p, [1, 1, 1, NaN, 1, 1, 1], 1), [y(1:3), NaN, y(5:7)], eps)
%!assert (nakainv (p, 1, [1, 1, 1, NaN, 1, 1, 1]), [y(1:3), NaN, y(5:7)], eps)
%!assert (nakainv ([p, NaN], 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (nakainv (single ([p, NaN]), 1, 1), single ([y, NaN]))
%!assert (nakainv ([p, NaN], single (1), 1), single ([y, NaN]))
%!assert (nakainv ([p, NaN], 1, single (1)), single ([y, NaN]))

## Test input validation
%!error<nakainv: function called with too few input arguments.> nakainv ()
%!error<nakainv: function called with too few input arguments.> nakainv (1)
%!error<nakainv: function called with too few input arguments.> nakainv (1, 2)
%!error<nakainv: P, MU, and OMEGA must be of common size or scalars.> ...
%! nakainv (ones (3), ones (2), ones(2))
%!error<nakainv: P, MU, and OMEGA must be of common size or scalars.> ...
%! nakainv (ones (2), ones (3), ones(2))
%!error<nakainv: P, MU, and OMEGA must be of common size or scalars.> ...
%! nakainv (ones (2), ones (2), ones(3))
%!error<nakainv: P, MU, and OMEGA must not be complex.> nakainv (i, 4, 3)
%!error<nakainv: P, MU, and OMEGA must not be complex.> nakainv (1, i, 3)
%!error<nakainv: P, MU, and OMEGA must not be complex.> nakainv (1, 4, i)
