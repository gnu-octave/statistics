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
## @deftypefn  {statistics} [@var{m}, @var{v}] = evstat (@var{mu}, @var{sigma})
##
## Mean and variance of the extreme value distribution.
##
## @code{[@var{m}, @var{v}] = evstat (@var{mu}, @var{sigma})} returns the mean
## and variance of the type 1 extreme value distribution with location parameter
## @var{mu} and scale parameter @var{sigma}.  The sizes of @var{m} and @var{v}
## are the common size of @var{mu} and @var{sigma}.  A scalar input functions as
## a constant matrix of the same size as the other inputs.
##
## The type 1 extreme value distribution is also known as the Gumbel
## distribution.  The version used here is suitable for modeling minima; the
## mirror image of this distribution can be used to model maxima by negating
## @var{x}.  If @var{y} has a Weibull distribution, then
## @code{@var{x} = log (@var{y})} has the type 1 extreme value distribution.
##
## @seealso{evcdf, evinv, evpdf, evrnd, evfit, evlike}
## @end deftypefn

function [m, v] = evstat (mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("evstat: too few input arguments.");
  endif
  ## Check for MU and SIGMA being reals
  if (iscomplex (mu) || iscomplex (sigma))
    error ("evstat: MU and SIGMA must not be complex.");
  endif
  ## Check for common size of MU and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [err, mu, sigma] = common_size (mu, sigma);
    if (err > 0)
      error ("evstat: MU and SIGMA must be of common size or scalars.");
    endif
  endif
  ## Return NaNs for out of range values of SIGMA
  sigma(sigma <= 0) = NaN;
  ## Calculate mean and variance
  m = mu + psi(1) .* sigma;
  v = (pi .* sigma) .^ 2 ./ 6;

endfunction

## Test input validation
%!error<evstat: too few input arguments.> evstat ()
%!error<evstat: too few input arguments.> evstat (1)
%!error<evstat: MU and SIGMA must not be complex.> evstat (i, 2)
%!error<evstat: MU and SIGMA must not be complex.> evstat (2, i)
%!error<evstat: MU and SIGMA must be of common size or scalars.> ...
%! evstat (ones (3), ones (2))

## Test results
%!shared x, y0, y1
%! x = [-5, 0, 1, 2, 3];
%! y0 = [NaN, NaN, 0.4228, 0.8456, 1.2684];
%! y1 = [-5.5772, -3.4633, -3.0405, -2.6177, -2.1949];
%!assert (evstat (x, x), y0, 1e-4)
%!assert (evstat (x, x+6), y1, 1e-4)
%!assert (evstat (x, x-6), NaN (1,5))
