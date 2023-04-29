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
## @deftypefn  {statistics} {@var{y} =} nakapdf (@var{x}, @var{mu}, @var{omega})
##
## Nakagami probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Nakagami distribution with shape parameter @var{mu} and
## spread parameter @var{omega}.  The size of @var{p} is the common size of
## @var{x}, @var{mu}, and @var{omega}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Further information about the Nakagami distribution can be found at
## @url{https://en.wikipedia.org/wiki/Nakagami_distribution}
##
## @seealso{nakacdf, nakapdf, nakarnd}
## @end deftypefn

function y = nakapdf (x, mu, omega)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("nakapdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, MU, and OMEGA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (omega))
    [retval, x, mu, omega] = common_size (x, mu, omega);
    if (retval > 0)
      error ("nakapdf: X, MU, and OMEGA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and OMEGA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (omega))
    error ("nakapdf: X, MU, and OMEGA must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (mu, "single") || isa (omega, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Compute Nakagami PDF
  k = isnan (x) | ! (mu > 0.5) | ! (omega > 0);
  y(k) = NaN;

  k = (0 < x) & (x < Inf) & (0 < mu) & (mu < Inf) & (0 < omega) & (omega < Inf);
  if (isscalar (mu) && isscalar(omega))
    y(k) = exp (log (2) + mu*log (mu) - log (gamma (mu)) - ...
               mu*log (omega) + (2*mu-1) * ...
               log (x(k)) - (mu/omega) * x(k).^2);
  else
    y(k) = exp(log(2) + mu(k).*log (mu(k)) - log (gamma (mu(k))) - ...
               mu(k).*log (omega(k)) + (2*mu(k)-1) ...
               .* log (x(k)) - (mu(k)./omega(k)) .* x(k).^2);
  endif

endfunction

%!demo
%! ## Plot various PDFs from the Nakagami distribution
%! x = 0:0.01:3;
%! y1 = nakapdf (x, 0.5, 1);
%! y2 = nakapdf (x, 1, 1);
%! y3 = nakapdf (x, 1, 2);
%! y4 = nakapdf (x, 1, 3);
%! y5 = nakapdf (x, 2, 1);
%! y6 = nakapdf (x, 2, 2);
%! y7 = nakapdf (x, 5, 1);
%! plot (x, y1, "-r", x, y2, "-g", x, y3, "-y", x, y4, "-m", ...
%!       x, y5, "-k", x, y6, "-b", x, y7, "-c")
%! grid on
%! xlim ([0, 3])
%! ylim ([0, 2])
%! legend ({"μ = 0.5, ω = 1", "μ = 1, ω = 1", "μ = 1, ω = 2", ...
%!          "μ = 1, ω = 3", "μ = 2, ω = 1", "μ = 2, ω = 2", ...
%!          "μ = 5, ω = 1"}, "location", "northeast")
%! title ("Nakagami PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test results
%!shared x, y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 0.73575888234288467, 0.073262555554936715, 0];
%!assert (nakapdf (x, ones (1,5), ones (1,5)), y, eps)
%!assert (nakapdf (x, 1, 1), y, eps)
%!assert (nakapdf (x, [1, 1, NaN, 1, 1], 1), [y(1:2), NaN, y(4:5)], eps)
%!assert (nakapdf (x, 1, [1, 1, NaN, 1, 1]), [y(1:2), NaN, y(4:5)], eps)
%!assert (nakapdf ([x, NaN], 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (nakapdf (single ([x, NaN]), 1, 1), single ([y, NaN]))
%!assert (nakapdf ([x, NaN], single (1), 1), single ([y, NaN]))
%!assert (nakapdf ([x, NaN], 1, single (1)), single ([y, NaN]))

## Test input validation
%!error<nakapdf: function called with too few input arguments.> nakapdf ()
%!error<nakapdf: function called with too few input arguments.> nakapdf (1)
%!error<nakapdf: function called with too few input arguments.> ...
%! nakapdf (1, 2)
%!error<nakapdf: function called with too many inputs> nakapdf (1, 2, 3, 4)
%!error<nakapdf: X, MU, and OMEGA must be of common size or scalars.> ...
%! nakapdf (ones (3), ones (2), ones(2))
%!error<nakapdf: X, MU, and OMEGA must be of common size or scalars.> ...
%! nakapdf (ones (2), ones (3), ones(2))
%!error<nakapdf: X, MU, and OMEGA must be of common size or scalars.> ...
%! nakapdf (ones (2), ones (2), ones(3))
%!error<nakapdf: X, MU, and OMEGA must not be complex.> nakapdf (i, 4, 3)
%!error<nakapdf: X, MU, and OMEGA must not be complex.> nakapdf (1, i, 3)
%!error<nakapdf: X, MU, and OMEGA must not be complex.> nakapdf (1, 4, i)
