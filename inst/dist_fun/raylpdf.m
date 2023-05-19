## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{y} =} raylpdf (@var{x}, @var{sigma})
##
## Rayleigh probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Rayleigh distribution with scale parameter @var{sigma}.  The size of
## @var{p} is the common size of @var{x} and @var{sigma}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## Further information about the Rayleigh distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rayleigh_distribution}
##
## @seealso{raylcdf, raylinv, raylrnd, raylstat}
## @end deftypefn

function y = raylpdf (x, sigma)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("raylpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X and SIGMA
  if (! isscalar (x) || ! isscalar (sigma))
    [retval, x, sigma] = common_size (x, sigma);
    if (retval > 0)
      error ("raylpdf: X and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X and SIGMA being reals
  if (iscomplex (x) || iscomplex (sigma))
    error ("raylpdf: X and SIGMA must not be complex.");
  endif

  ## Calculate Rayleigh PDF
  y = x .* exp ((-x .^ 2) ./ (2 .* sigma .^ 2)) ./ (sigma .^ 2);

  ## Continue argument check
  k = find (! (x >= 0) | ! (x < Inf) | ! (sigma > 0));
  if (any (k))
    y(k) = NaN;
  endif

endfunction

%!demo
%! ## Plot various PDFs from the Rayleigh distribution
%! x = 0:0.01:10;
%! y1 = raylpdf (x, 0.5);
%! y2 = raylpdf (x, 1);
%! y3 = raylpdf (x, 2);
%! y4 = raylpdf (x, 3);
%! y5 = raylpdf (x, 4);
%! plot (x, y1, "-b", x, y2, "g", x, y3, "-r", x, y4, "-m", x, y5, "-k")
%! grid on
%! ylim ([0, 1.25])
%! legend ({"σ = 0,5", "σ = 1", "σ = 2", ...
%!          "σ = 3", "σ = 4"}, "location", "northeast")
%! title ("Rayleigh PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!test
%! x = 0:0.5:2.5;
%! sigma = 1:6;
%! y = raylpdf (x, sigma);
%! expected_y = [0.0000, 0.1212, 0.1051, 0.0874, 0.0738, 0.0637];
%! assert (y, expected_y, 0.001);
%!test
%! x = 0:0.5:2.5;
%! y = raylpdf (x, 0.5);
%! expected_y = [0.0000, 1.2131, 0.5413, 0.0667, 0.0027, 0.0000];
%! assert (y, expected_y, 0.001);

## Test input validation
%!error<raylpdf: function called with too few input arguments.> raylpdf ()
%!error<raylpdf: function called with too few input arguments.> raylpdf (1)
%!error<raylpdf: X and SIGMA must be of common size or scalars.> ...
%! raylpdf (ones (3), ones (2))
%!error<raylpdf: X and SIGMA must be of common size or scalars.> ...
%! raylpdf (ones (2), ones (3))
%!error<raylpdf: X and SIGMA must not be complex.> raylpdf (i, 2)
%!error<raylpdf: X and SIGMA must not be complex.> raylpdf (2, i)
