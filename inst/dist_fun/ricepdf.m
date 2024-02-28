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
## @deftypefn  {statistics} {@var{y} =} ricepdf (@var{x}, @var{nu}, @var{sigma})
##
## Rician probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Rician distribution with non-centrality (distance) parameter @var{nu}
## and scale parameter @var{sigma}.  The size of @var{y} is the common size of
## @var{x}, @var{nu}, and @var{sigma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Further information about the Rician distribution can be found at
## @url{https://en.wikipedia.org/wiki/Rice_distribution}
##
## @seealso{ricecdf, riceinv, ricernd, ricefit, ricelike, ricestat}
## @end deftypefn

function y = ricepdf (x, nu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("ricepdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, NU, and SIGMA
  if (! isscalar (x) || ! isscalar (nu) || ! isscalar (sigma))
    [retval, x, nu, sigma] = common_size (x, nu, sigma);
    if (retval > 0)
      error ("ricepdf: X, NU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, NU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (nu) || iscomplex (sigma))
    error ("ricepdf: X, NU, and SIGMA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (nu, "single") || isa (sigma, "single"));
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  k = nu < 0 | sigma <= 0 | x < 0 | isnan (x) | isnan (nu) | isnan (sigma);
  y(k) = NaN;

  k = ! k;
  ## Do the math
  x_k = x(k);
  n_k = nu(k);
  s_sq = sigma(k) .^ 2;
  x_s2 = x_k ./ s_sq;
  xnsq = (x_k .^ 2 + n_k .^ 2) ./ (2 .* s_sq);
  epxt = xnsq - x_s2 .* n_k;
  term = exp (-epxt);
  y(k) = x_s2 .* term .* besseli (0, x_s2 .* n_k, 1);

  ## Fix arithmetic overflow due to exponent
  y(epxt > (log(realmax(class(y))))) = 0;

  ## Fix x < 0 -> 0
  y(x < 0) = 0;

endfunction

%!demo
%! ## Plot various PDFs from the Rician distribution
%! x = 0:0.01:8;
%! y1 = ricepdf (x, 0, 1);
%! y2 = ricepdf (x, 0.5, 1);
%! y3 = ricepdf (x, 1, 1);
%! y4 = ricepdf (x, 2, 1);
%! y5 = ricepdf (x, 4, 1);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-m", x, y5, "-k")
%! grid on
%! ylim ([0, 0.65])
%! xlim ([0, 8])
%! legend ({"ν = 0, σ = 1", "ν = 0.5, σ = 1", "ν = 1, σ = 1", ...
%!          "ν = 2, σ = 1", "ν = 4, σ = 1"}, "location", "northeast")
%! title ("Rician PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 0.5 1 2];
%! y = [0 0 0.1073 0.1978 0.2846];
%!assert (ricepdf (x, ones (1, 5), 2 * ones (1, 5)), y, 1e-4)
%!assert (ricepdf (x, 1, 2 * ones (1, 5)), y, 1e-4)
%!assert (ricepdf (x, ones (1, 5), 2), y, 1e-4)
%!assert (ricepdf (x, [0 NaN 1 1 1], 2), [0 NaN y(3:5)], 1e-4)
%!assert (ricepdf (x, 1, 2 * [0 NaN 1 1 1]), [0 NaN y(3:5)], 1e-4)
%!assert (ricepdf ([x, NaN], 1, 2), [y, NaN], 1e-4)

## Test class of input preserved
%!assert (ricepdf (single ([x, NaN]), 1, 2), single ([y, NaN]), 1e-4)
%!assert (ricepdf ([x, NaN], single (1), 2), single ([y, NaN]), 1e-4)
%!assert (ricepdf ([x, NaN], 1, single (2)), single ([y, NaN]), 1e-4)

## Test input validation
%!error<ricepdf: function called with too few input arguments.> ricepdf ()
%!error<ricepdf: function called with too few input arguments.> ricepdf (1)
%!error<ricepdf: function called with too few input arguments.> ricepdf (1,2)
%!error<ricepdf: function called with too many inputs> ricepdf (1,2,3,4)
%!error<ricepdf: X, NU, and SIGMA must be of common size or scalars.> ...
%! ricepdf (ones (3), ones (2), ones (2))
%!error<ricepdf: X, NU, and SIGMA must be of common size or scalars.> ...
%! ricepdf (ones (2), ones (3), ones (2))
%!error<ricepdf: X, NU, and SIGMA must be of common size or scalars.> ...
%! ricepdf (ones (2), ones (2), ones (3))
%!error<ricepdf: X, NU, and SIGMA must not be complex.> ricepdf (i, 2, 2)
%!error<ricepdf: X, NU, and SIGMA must not be complex.> ricepdf (2, i, 2)
%!error<ricepdf: X, NU, and SIGMA must not be complex.> ricepdf (2, 2, i)
