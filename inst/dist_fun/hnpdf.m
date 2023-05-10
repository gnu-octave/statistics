## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{y} =} hnpdf (@var{x}, @var{mu}, @var{sigma})
##
## Half-normal probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the half-normal distribution with location parameter @var{mu} and scale
## parameter @var{sigma}.  The size of @var{y} is the common size of @var{x},
## @var{mu}, and @var{sigma}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## The half-normal CDF is only defined for @qcode{@var{x} >= @var{mu}}.
##
## Further information about the half-normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Half-normal_distribution}
##
## @seealso{hncdf, hninv, hnrnd, hnfit, hnlike}
## @end deftypefn

function y = hnpdf (x, mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("hnpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, MU, and SIGMA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(sigma))
    [retval, x, mu, sigma] = common_size (x, mu, sigma);
    if (retval > 0)
      error ("hnpdf: X, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("hnpdf: X, MU, and SIGMA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"));
    y = NaN (size (x), "single");
  else
    y = NaN (size (x));
  endif

  ## Return NaNs for out of range values of SIGMA parameter
  sigma(sigma <= 0) = NaN;

  ## Compute half-normal PDF
  z = (x - mu) ./ sigma;
  y = sqrt (2 / pi) ./ sigma .* exp (-0.5 * z .^ 2);

  ## Force zero for unsupported X
  y(z < 0) = 0;

endfunction

%!demo
%! ## Plot various PDFs from the half-normal distribution
%! x = 0:0.001:10;
%! y1 = hnpdf (x, 0, 1);
%! y2 = hnpdf (x, 0, 2);
%! y3 = hnpdf (x, 0, 3);
%! y4 = hnpdf (x, 0, 5);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c")
%! grid on
%! xlim ([0, 10])
%! ylim ([0, 0.9])
%! legend ({"μ = 0, σ = 1", "μ = 0, σ = 2", ...
%!          "μ = 0, σ = 3", "μ = 0, σ = 5"}, "location", "northeast")
%! title ("Half-normal PDF")
%! xlabel ("values in x")
%! ylabel ("density")

%!demo
%! ## Plot half-normal against normal probability density function
%! x = -5:0.001:5;
%! y1 = hnpdf (x, 0, 1);
%! y2 = normpdf (x);
%! plot (x, y1, "-b", x, y2, "-g")
%! grid on
%! xlim ([-5, 5])
%! ylim ([0, 0.9])
%! legend ({"half-normal with μ = 0, σ = 1", ...
%!          "standart normal (μ = 0, σ = 1)"}, "location", "northeast")
%! title ("Half-normal against standard normal PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-Inf, -1, 0, 1/2, 1, Inf];
%! y = [0, 0, 0.7979, 0.7041, 0.4839, 0];
%!assert (hnpdf ([x, NaN], 0, 1), [y, NaN], 1e-4)
%!assert (hnpdf (x, 0, [-2, -1, 0, 1, 1, 1]), [nan(1,3), y([4:6])], 1e-4)

## Test class of input preserved
%!assert (class (hncdf (single ([x, NaN]), 0, 1)), "single")
%!assert (class (hncdf ([x, NaN], 0, single (1))), "single")
%!assert (class (hncdf ([x, NaN], single (0), 1)), "single")

## Test input validation
%!error<hnpdf: function called with too few input arguments.> hnpdf ()
%!error<hnpdf: function called with too few input arguments.> hnpdf (1)
%!error<hnpdf: function called with too few input arguments.> hnpdf (1, 2)
%!error<hnpdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! hnpdf (1, ones (2), ones (3))
%!error<hnpdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! hnpdf (ones (2), 1, ones (3))
%!error<hnpdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! hnpdf (ones (2), ones (3), 1)
%!error<hnpdf: X, MU, and SIGMA must not be complex.> hnpdf (i, 2, 3)
%!error<hnpdf: X, MU, and SIGMA must not be complex.> hnpdf (1, i, 3)
%!error<hnpdf: X, MU, and SIGMA must not be complex.> hnpdf (1, 2, i)
