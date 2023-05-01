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
## @deftypefn  {statistics} {@var{y} =} gumbelpdf (@var{x})
## @deftypefnx {statistics} {@var{y} =} gumbelpdf (@var{x}, @var{mu})
## @deftypefnx {statistics} {@var{y} =} gumbelpdf (@var{x}, @var{mu}, @var{beta})
##
## Gumbel probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Gumbel distribution (also known as the extreme value or the
## type I generalized extreme value distribution) with location parameter
## @var{mu} and scale parameter @var{beta}.  The size of @var{y} is the common
## size of @var{x}, @var{mu} and @var{beta}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## Default values are @var{mu} = 0, @var{beta} = 1.
##
## The Gumbel distribution is used to model the distribution of the maximum (or
## the minimum) of a number of samples of various distributions.  This version
## is suitable for modeling maxima.  For modeling minima, use the alternative
## extreme value iCDF, @code{evpdf}.
##
## Further information about the Gumbel distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gumbel_distribution}
##
## @seealso{gumbelcdf, gumbelinv, gumbelrnd, gumbelfit, gumbellike, gumbelstat,
## evpdf}
## @end deftypefn

function y = gumbelpdf (x, mu, beta)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("gumbelpdf: too few input arguments.");
  endif

  ## Add defaults (if missing input arguments)
  if (nargin < 2)
    mu = 0;
  endif
  if (nargin < 3)
    beta = 1;
  endif

  ## Check for common size of X, MU, and BETA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (beta))
    [err, x, mu, beta] = common_size (x, mu, beta);
    if (err > 0)
      error ("gumbelpdf: X, MU, and BETA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and BETA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (beta))
    error ("gumbelpdf: X, MU, and BETA must not be complex.");
  endif

  ## Return NaNs for out of range parameters
  beta(beta <= 0) = NaN;

  ## Compute pdf of type 1 extreme value distribution
  z = -(x - mu) ./ beta;
  y = exp (z - exp (z)) ./ beta;

  ## Force 0 for extreme right tail, instead of getting exp (Inf - Inf) = NaN
  y(z == Inf) = 0;

endfunction

%!demo
%! ## Plot various PDFs from the Extreme value distribution
%! x = -5:0.001:20;
%! y1 = gumbelpdf (x, 0.5, 2);
%! y2 = gumbelpdf (x, 1.0, 2);
%! y3 = gumbelpdf (x, 1.5, 3);
%! y4 = gumbelpdf (x, 3.0, 4);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c")
%! grid on
%! ylim ([0, 0.2])
%! legend ({"μ = 0.5, β = 2", "μ = 1.0, β = 2", ...
%!          "μ = 1.5, β = 3", "μ = 3.0, β = 4"}, "location", "northeast")
%! title ("Extreme value PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test results
%!shared x, y0, y1
%! x = [-5, 0, 1, 2, 3];
%! y0 = [0, 0.3679, 0.2547, 0.1182, 0.0474];
%! y1 = [0, 0.1794, 0.3679, 0.2547, 0.1182];
%!assert (gumbelpdf (x), y0, 1e-4)
%!assert (gumbelpdf (x, zeros (1,5), ones (1,5)), y0, 1e-4)
%!assert (gumbelpdf (x, ones (1,5), ones (1,5)), y1, 1e-4)

## Test input validation
%!error<gumbelpdf: too few input arguments.> gumbelpdf ()
%!error<gumbelpdf: X, MU, and BETA must be of common size or scalars.> ...
%! gumbelpdf (ones (3), ones (2), ones (2))
%!error<gumbelpdf: X, MU, and BETA must not be complex.> gumbelpdf (i, 2, 2)
%!error<gumbelpdf: X, MU, and BETA must not be complex.> gumbelpdf (2, i, 2)
%!error<gumbelpdf: X, MU, and BETA must not be complex.> gumbelpdf (2, 2, i)

