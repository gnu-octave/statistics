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
## @deftypefn  {statistics} {@var{p} =} invgcdf (@var{x}, @var{mu}, @var{lambda})
## @deftypefnx {statistics} {@var{p} =} invgcdf (@var{x}, @var{mu}, @var{lambda}, @qcode{"upper"})
##
## Inverse Gaussian cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the inverse Gaussian distribution with scale parameter @var{mu} and
## shape parameter @var{lambda}.  The size of @var{p} is the common size of
## @var{x}, @var{mu} and @var{lambda}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## @code{[@dots{}] = invgcdf(@dots{}, "upper")} computes the upper tail
## probability of the inverse Gaussian distribution.
##
## The inverse Gaussian CDF is only defined for @qcode{@var{mu} > 0} and
## @qcode{@var{lambda} > 0}.
##
## Further information about the inverse Gaussian distribution can be found at
## @url{https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}
##
## @seealso{invginv, invgpdf, invgrnd, invgfit, invglike}
## @end deftypefn

function p = invgcdf (x, mu, lambda, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("invgcdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("invgcdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, MU, and LAMBDA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (lambda))
    [err, x, mu, lambda] = common_size (x, mu, lambda);
    if (err > 0)
      error ("invgcdf: X, MU, and LAMBDA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and LAMBDA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (lambda))
    error ("invgcdf: X, MU, and LAMBDA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (lambda, "single"))
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Prepare output
  p = zeros (size (x), is_class);

  ## Return NaNs for out of range values of MU and LAMBDA parameters
  mu(mu <= 0) = NaN;
  lambda(lambda <= 0) = NaN;

  ## Check for valid support of X
  is_zero = (x <= 0);
  x(is_zero) = realmin;
  is_inf = (x == Inf);

  ## Calculate z1, z2
  z1 = sqrt (lambda ./ x) .* (x ./ mu - 1);
  z2 = -sqrt (lambda ./ x) .* (x ./ mu + 1);

  ## Compute the CDF if the inverse Gaussian
  if (uflag)
    p = 0.5 .* erfc (z1 ./ sqrt (2)) - ...
               exp (2 .* lambda ./ mu) .* 0.5 .* erfc (-z2 ./ sqrt (2));
    p(is_zero) = 1;
    p(is_inf) = 0;
  else
    p = 0.5 .* erfc (-z1 ./ sqrt (2)) + ...
               exp (2 .* lambda ./ mu) .* 0.5 .* erfc (-z2 ./ sqrt (2));
    p(is_zero) = 0;
    p(is_inf) = 1;
  endif

endfunction

%!demo
%! ## Plot various CDFs from the inverse Gaussian distribution
%! x = 0:0.001:3;
%! p1 = invgcdf (x, 1, 0.2);
%! p2 = invgcdf (x, 1, 1);
%! p3 = invgcdf (x, 1, 3);
%! p4 = invgcdf (x, 3, 0.2);
%! p5 = invgcdf (x, 3, 1);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c", x, p5, "-y")
%! grid on
%! xlim ([0, 3])
%! legend ({"μ = 1, σ = 0.2", "μ = 1, σ = 1", "μ = 1, σ = 3", ...
%!          "μ = 3, σ = 0.2", "μ = 3, σ = 1"}, "location", "southeast")
%! title ("Inverse Gaussian CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, p1, p1u, y2, y2u, y3, y3u
%! x = [-Inf, -1, 0, 1/2, 1, Inf];
%! p1 = [0, 0, 0, 0.3650, 0.6681, 1];
%! p1u = [1, 1, 1, 0.6350, 0.3319, 0];
%!assert (invgcdf (x, ones (1,6), ones (1,6)), p1, 1e-4)
%!assert (invgcdf (x, 1, 1), p1, 1e-4)
%!assert (invgcdf (x, 1, ones (1,6)), p1, 1e-4)
%!assert (invgcdf (x, ones (1,6), 1), p1, 1e-4)
%!assert (invgcdf (x, 1, [1, 1, 1, NaN, 1, 1]), [p1(1:3), NaN, p1(5:6)], 1e-4)
%!assert (invgcdf (x, [1, 1, 1, NaN, 1, 1], 1), [p1(1:3), NaN, p1(5:6)], 1e-4)
%!assert (invgcdf ([x(1:3), NaN, x(5:6)], 1, 1), [p1(1:3), NaN, p1(5:6)], 1e-4)
%!assert (invgcdf (x, ones (1,6), ones (1,6), "upper"), p1u, 1e-4)
%!assert (invgcdf (x, 1, 1, "upper"), p1u, 1e-4)
%!assert (invgcdf (x, 1, ones (1,6), "upper"), p1u, 1e-4)
%!assert (invgcdf (x, ones (1,6), 1, "upper"), p1u, 1e-4)

## Test class of input preserved
%!assert (class (invgcdf (single ([x, NaN]), 1, 1)), "single")
%!assert (class (invgcdf ([x, NaN], 1, single (1))), "single")
%!assert (class (invgcdf ([x, NaN], single (1), 1)), "single")

## Test input validation
%!error<invgcdf: function called with too few input arguments.> invgcdf ()
%!error<invgcdf: function called with too few input arguments.> invgcdf (1)
%!error<invgcdf: function called with too few input arguments.> invgcdf (1, 2)
%!error<invgcdf: invalid argument for upper tail.> invgcdf (1, 2, 3, "tail")
%!error<invgcdf: invalid argument for upper tail.> invgcdf (1, 2, 3, 5)
%!error<invgcdf: X, MU, and LAMBDA must be of common size or scalars.> ...
%! invgcdf (ones (3), ones (2), ones(2))
%!error<invgcdf: X, MU, and LAMBDA must be of common size or scalars.> ...
%! invgcdf (ones (2), ones (3), ones(2))
%!error<invgcdf: X, MU, and LAMBDA must be of common size or scalars.> ...
%! invgcdf (ones (2), ones (2), ones(3))
%!error<invgcdf: X, MU, and LAMBDA must not be complex.> invgcdf (i, 2, 3)
%!error<invgcdf: X, MU, and LAMBDA must not be complex.> invgcdf (1, i, 3)
%!error<invgcdf: X, MU, and LAMBDA must not be complex.> invgcdf (1, 2, i)
