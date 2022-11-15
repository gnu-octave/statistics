## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{y} = evpdf (@var{x})
## @deftypefnx {Function File} @var{y} = evpdf (@var{x}, @var{mu})
## @deftypefnx {Function File} @var{y} = evpdf (@var{x}, @var{mu}, @var{sigma})
##
## Extreme value probability density function (pdf).
##
## @code{@var{y} = evpdf (@var{x}, @var{mu}, @var{sigma})} returns the pdf of
## the type 1 extreme value distribution with location parameter @var{mu} and
## scale parameter @var{sigma}.  The size of @var{x} is the common size of
## @var{p}, @var{mu} and @var{sigma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Default values are @var{mu} = 0, @var{sigma} = 1.
##
## The type 1 extreme value distribution is also known as the Gumbel
## distribution.  The version used here is suitable for modeling minima; the
## mirror image of this distribution can be used to model maxima by negating
## @var{x}.  If @var{y} has a Weibull distribution, then
## @code{@var{x} = log (@var{y})} has the type 1 extreme value distribution.
##
## @seealso{evcdf, evinv, evrnd, evfit, evlike, evstat}
## @end deftypefn

function y = evpdf (x, mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("evpdf: too few input arguments.");
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

## Test input validation
%!error<evpdf: too few input arguments.> evpdf ()
%!error<evpdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! evpdf (ones (3), ones (2), ones (2))
%!error<evpdf: X, MU, and SIGMA must not be complex.> evpdf (i, 2, 2)
%!error<evpdf: X, MU, and SIGMA must not be complex.> evpdf (2, i, 2)
%!error<evpdf: X, MU, and SIGMA must not be complex.> evpdf (2, 2, i)

## Test results
%!shared x, y0, y1
%! x = [-5, 0, 1, 2, 3];
%! y0 = [0.0067, 0.3679, 0.1794, 0.0046, 0];
%! y1 = [0.0025, 0.2546, 0.3679, 0.1794, 0.0046];
%!assert (evpdf (x), y0, 1e-4)
%!assert (evpdf (x, zeros (1,5), ones (1,5)), y0, 1e-4)
%!assert (evpdf (x, ones (1,5), ones (1,5)), y1, 1e-4)

