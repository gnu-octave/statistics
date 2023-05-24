## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{y} =} evpdf (@var{x})
## @deftypefnx {statistics} {@var{y} =} evpdf (@var{x}, @var{mu})
## @deftypefnx {statistics} {@var{y} =} evpdf (@var{x}, @var{mu}, @var{sigma})
##
## Extreme value probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the extreme value distribution (also known as the Gumbel or the type I
## generalized extreme value distribution) with location parameter @var{mu} and
## scale parameter @var{sigma}.  The size of @var{y} is the common size of
## @var{x}, @var{mu} and @var{sigma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Default values are @var{mu} = 0 and @var{sigma} = 1.
##
## The Gumbel distribution is used to model the distribution of the maximum (or
## the minimum) of a number of samples of various distributions.  This version
## is suitable for modeling minima.  For modeling maxima, use the alternative
## Gumbel iCDF, @code{gumbelinv}.
##
## Further information about the Gumbel distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gumbel_distribution}
##
## @seealso{evcdf, evinv, evrnd, evfit, evlike, evstat, gumbelpdf}
## @end deftypefn

function y = evpdf (x, mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("evpdf: function called with too few input arguments.");
  endif

  ## Add defaults (if missing input arguments)
  if (nargin < 2)
    mu = 0;
  endif
  if (nargin < 3)
    sigma = 1;
  endif

  ## Check for common size of X, MU, and SIGMA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (sigma))
    [err, x, mu, sigma] = common_size (x, mu, sigma);
    if (err > 0)
      error ("evpdf: X, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("evpdf: X, MU, and SIGMA must not be complex.");
  endif

  ## Return NaNs for out of range parameters
  sigma(sigma <= 0) = NaN;

  ## Compute pdf of type 1 extreme value distribution
  z = (x - mu) ./ sigma;
  y = exp (z - exp (z)) ./ sigma;

  ## Force 0 for extreme right tail, instead of getting exp (Inf - Inf) = NaN
  y(z == Inf) = 0;

endfunction

%!demo
%! ## Plot various PDFs from the Extreme value distribution
%! x = -10:0.001:10;
%! y1 = evpdf (x, 0.5, 2);
%! y2 = evpdf (x, 1.0, 2);
%! y3 = evpdf (x, 1.5, 3);
%! y4 = evpdf (x, 3.0, 4);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c")
%! grid on
%! ylim ([0, 0.2])
%! legend ({"μ = 0.5, σ = 2", "μ = 1.0, σ = 2", ...
%!          "μ = 1.5, σ = 3", "μ = 3.0, σ = 4"}, "location", "northeast")
%! title ("Extreme value PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y0, y1
%! x = [-5, 0, 1, 2, 3];
%! y0 = [0.0067, 0.3679, 0.1794, 0.0046, 0];
%! y1 = [0.0025, 0.2546, 0.3679, 0.1794, 0.0046];
%!assert (evpdf (x), y0, 1e-4)
%!assert (evpdf (x, zeros (1,5), ones (1,5)), y0, 1e-4)
%!assert (evpdf (x, ones (1,5), ones (1,5)), y1, 1e-4)

## Test input validation
%!error<evpdf: function called with too few input arguments.> evpdf ()
%!error<evpdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! evpdf (ones (3), ones (2), ones (2))
%!error<evpdf: X, MU, and SIGMA must not be complex.> evpdf (i, 2, 2)
%!error<evpdf: X, MU, and SIGMA must not be complex.> evpdf (2, i, 2)
%!error<evpdf: X, MU, and SIGMA must not be complex.> evpdf (2, 2, i)
