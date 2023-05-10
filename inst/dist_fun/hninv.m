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
## @deftypefn  {statistics} {@var{x} =} hninv (@var{p}, @var{mu}, @var{sigma})
##
## Inverse of the half-normal cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the half-normal distribution with location parameter @var{mu} and scale
## parameter @var{sigma}.  The size of @var{x} is the common size of @var{p},
## @var{mu}, and @var{sigma}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Further information about the half-normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Half-normal_distribution}
##
## @seealso{hncdf, hnpdf, hnrnd, hnfit, hnlike}
## @end deftypefn

function x = hninv (p, mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("hninv: function called with too few input arguments.");
  endif

  ## Check for common size of P, MU, and SIGMA
  if (! isscalar (p) || ! isscalar (mu) || ! isscalar(sigma))
    [retval, p, mu, sigma] = common_size (p, mu, sigma);
    if (retval > 0)
      error ("hninv: P, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (sigma))
    error ("hninv: P, MU, and SIGMA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (mu, "single") || isa (sigma, "single"));
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Return NaNs for out of range values of SIGMA parameter
  sigma(sigma <= 0) = NaN;

  ## Return NaNs for out of range P values
  p(p < 0 | 1 < p) = NaN;

  ## Calculate the quantile of half-normal distribution
  x = sqrt(2) * sigma .* erfinv (p) + mu;

endfunction

%!demo
%! ## Plot various iCDFs from the half-normal distribution
%! p = 0.001:0.001:0.999;
%! x1 = hninv (p, 0, 1);
%! x2 = hninv (p, 0, 2);
%! x3 = hninv (p, 0, 3);
%! x4 = hninv (p, 0, 5);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c")
%! grid on
%! ylim ([0, 10])
%! legend ({"μ = 0, σ = 1", "μ = 0, σ = 2", ...
%!          "μ = 0, σ = 3", "μ = 0, σ = 5"}, "location", "northwest")
%! title ("Half-normal iCDF")
%! xlabel ("probability")
%! ylabel ("x")

## Test output
%!shared p, x
%! p = [0, 0.3829, 0.6827, 1];
%! x = [0, 1/2, 1, Inf];
%!assert (hninv (p, 0, 1), x, 1e-4);
%!assert (hninv (p, 5, 1), x + 5, 1e-4);
%!assert (hninv (p, 0, ones (1,4)), x, 1e-4);
%!assert (hninv (p, 0, [-1, 0, 1, 1]), [NaN, NaN, x(3:4)], 1e-4)

## Test class of input preserved
%!assert (class (hninv (single ([p, NaN]), 0, 1)), "single")
%!assert (class (hninv ([p, NaN], single (0), 1)), "single")
%!assert (class (hninv ([p, NaN], 0, single (1))), "single")

## Test input validation
%!error<hninv: function called with too few input arguments.> hninv (1)
%!error<hninv: function called with too few input arguments.> hninv (1, 2)
%!error<hninv: P, MU, and SIGMA must be of common size or scalars.> ...
%! hninv (1, ones (2), ones (3))
%!error<hninv: P, MU, and SIGMA must be of common size or scalars.> ...
%! hninv (ones (2), 1, ones (3))
%!error<hninv: P, MU, and SIGMA must be of common size or scalars.> ...
%! hninv (ones (2), ones (3), 1)
%!error<hninv: P, MU, and SIGMA must not be complex.> hninv (i, 2, 3)
%!error<hninv: P, MU, and SIGMA must not be complex.> hninv (1, i, 3)
%!error<hninv: P, MU, and SIGMA must not be complex.> hninv (1, 2, i)
