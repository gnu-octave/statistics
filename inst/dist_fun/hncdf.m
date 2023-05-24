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
## @deftypefn  {statistics} {@var{p} =} hncdf (@var{x}, @var{mu}, @var{sigma})
## @deftypefnx {statistics} {@var{p} =} hncdf (@var{x}, @var{mu}, @var{sigma}, @qcode{"upper"})
##
## Half-normal cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the half-normal distribution with location parameter @var{mu} and
## scale parameter @var{sigma}.  The size of @var{p} is the common size of
## @var{x}, @var{mu} and @var{sigma}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## @code{[@dots{}] = hncdf (@var{x}, @var{mu}, @var{sigma}, "upper")} computes
## the upper tail probability of the half-normal distribution with parameters
## @var{mu} and @var{sigma}, at the values in @var{x}.
##
## The half-normal CDF is only defined for @qcode{@var{x} >= @var{mu}}.
##
## Further information about the half-normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Half-normal_distribution}
##
## @seealso{hninv, hnpdf, hnrnd, hnfit, hnlike}
## @end deftypefn

function p = hncdf (x, mu, sigma, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("hncdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("hncdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, MU, and SIGMA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (sigma))
    [err, x, mu, sigma] = common_size (x, mu, sigma);
    if (err > 0)
      error ("hncdf: X, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("hncdf: X, MU, and SIGMA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"))
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Prepare output
  p = zeros (size (x), is_class);

  ## Return NaNs for out of range values of SIGMA parameter
  sigma(sigma <= 0) = NaN;

  ## Calculate (x-mu)/sigma => 0 and force zero below that
  z = (x - mu) ./ sigma;
  z(z < 0) = 0;

  if (uflag)
    p = erfc(z./sqrt(2));
  else
    p = erf(z./sqrt(2));
  endif

endfunction

%!demo
%! ## Plot various CDFs from the half-normal distribution
%! x = 0:0.001:10;
%! p1 = hncdf (x, 0, 1);
%! p2 = hncdf (x, 0, 2);
%! p3 = hncdf (x, 0, 3);
%! p4 = hncdf (x, 0, 5);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c")
%! grid on
%! xlim ([0, 10])
%! legend ({"μ = 0, σ = 1", "μ = 0, σ = 2", ...
%!          "μ = 0, σ = 3", "μ = 0, σ = 5"}, "location", "southeast")
%! title ("Half-normal CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

%!demo
%! ## Plot half-normal against normal cumulative distribution function
%! x = -5:0.001:5;
%! p1 = hncdf (x, 0, 1);
%! p2 = normcdf (x);
%! plot (x, p1, "-b", x, p2, "-g")
%! grid on
%! xlim ([-5, 5])
%! legend ({"half-normal with μ = 0, σ = 1", ...
%!          "standart normal (μ = 0, σ = 1)"}, "location", "southeast")
%! title ("Half-normal against standard normal CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, p1, p1u, y2, y2u, y3, y3u
%! x = [-Inf, -1, 0, 1/2, 1, Inf];
%! p1 = [0, 0, 0, 0.3829, 0.6827, 1];
%! p1u = [1, 1, 1, 0.6171, 0.3173, 0];
%!assert (hncdf (x, zeros (1,6), ones (1,6)), p1, 1e-4)
%!assert (hncdf (x, 0, 1), p1, 1e-4)
%!assert (hncdf (x, 0, ones (1,6)), p1, 1e-4)
%!assert (hncdf (x, zeros (1,6), 1), p1, 1e-4)
%!assert (hncdf (x, 0, [1, 1, 1, NaN, 1, 1]), [p1(1:3), NaN, p1(5:6)], 1e-4)
%!assert (hncdf (x, [0, 0, 0, NaN, 0, 0], 1), [p1(1:3), NaN, p1(5:6)], 1e-4)
%!assert (hncdf ([x(1:3), NaN, x(5:6)], 0, 1), [p1(1:3), NaN, p1(5:6)], 1e-4)
%!assert (hncdf (x, zeros (1,6), ones (1,6), "upper"), p1u, 1e-4)
%!assert (hncdf (x, 0, 1, "upper"), p1u, 1e-4)
%!assert (hncdf (x, 0, ones (1,6), "upper"), p1u, 1e-4)
%!assert (hncdf (x, zeros (1,6), 1, "upper"), p1u, 1e-4)

## Test class of input preserved
%!assert (class (hncdf (single ([x, NaN]), 0, 1)), "single")
%!assert (class (hncdf ([x, NaN], 0, single (1))), "single")
%!assert (class (hncdf ([x, NaN], single (0), 1)), "single")

## Test input validation
%!error<hncdf: function called with too few input arguments.> hncdf ()
%!error<hncdf: function called with too few input arguments.> hncdf (1)
%!error<hncdf: function called with too few input arguments.> hncdf (1, 2)
%!error<hncdf: invalid argument for upper tail.> hncdf (1, 2, 3, "tail")
%!error<hncdf: invalid argument for upper tail.> hncdf (1, 2, 3, 5)
%!error<hncdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! hncdf (ones (3), ones (2), ones(2))
%!error<hncdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! hncdf (ones (2), ones (3), ones(2))
%!error<hncdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! hncdf (ones (2), ones (2), ones(3))
%!error<hncdf: X, MU, and SIGMA must not be complex.> hncdf (i, 2, 3)
%!error<hncdf: X, MU, and SIGMA must not be complex.> hncdf (1, i, 3)
%!error<hncdf: X, MU, and SIGMA must not be complex.> hncdf (1, 2, i)
