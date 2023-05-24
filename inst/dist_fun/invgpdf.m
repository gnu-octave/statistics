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
## @deftypefn  {statistics} {@var{y} =} invgpdf (@var{x}, @var{mu}, @var{lambda})
##
## Inverse Gaussian probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the inverse Gaussian distribution with scale parameter @var{mu} and shape
## parameter @var{lambda}.  The size of @var{y} is the common size of @var{x},
## @var{mu}, and @var{lambda}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## The inverse Gaussian CDF is only defined for @qcode{@var{mu} > 0} and
## @qcode{@var{lambda} > 0}.
##
## Further information about the inverse Gaussian distribution can be found at
## @url{https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}
##
## @seealso{invgcdf, invginv, invgrnd, invgfit, invglike}
## @end deftypefn

function y = invgpdf (x, mu, lambda)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("invgpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, MU, and LAMBDA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(lambda))
    [retval, x, mu, lambda] = common_size (x, mu, lambda);
    if (retval > 0)
      error ("invgpdf: X, MU, and LAMBDA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and LAMBDA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (lambda))
    error ("invgpdf: X, MU, and LAMBDA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (lambda, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Return NaNs for out of range values of MU and LAMBDA parameters
  mu(mu <= 0) = NaN;
  lambda(lambda <= 0) = NaN;

  ## Check for valid support of X
  is_zero = (x <= 0);
  x(is_zero) = realmin;

  ## Compute inverse Gaussian PDF
  y = sqrt (lambda ./ (2 .* pi .* x .^ 3)) .* ...
      exp (-0.5 .* lambda .* (x ./ mu - 2 + mu ./ x) ./ mu);

  ## Force zero for unsupported X but valid parameters
  k0 = is_zero & mu > 0 & lambda > 0;
  y(k0) = 0;

  ## Cast to appropriate class
  y = cast (y, is_class);

endfunction

%!demo
%! ## Plot various PDFs from the inverse Gaussian distribution
%! x = 0:0.001:3;
%! y1 = invgpdf (x, 1, 0.2);
%! y2 = invgpdf (x, 1, 1);
%! y3 = invgpdf (x, 1, 3);
%! y4 = invgpdf (x, 3, 0.2);
%! y5 = invgpdf (x, 3, 1);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c", x, y5, "-y")
%! grid on
%! xlim ([0, 3])
%! ylim ([0, 3])
%! legend ({"μ = 1, σ = 0.2", "μ = 1, σ = 1", "μ = 1, σ = 3", ...
%!          "μ = 3, σ = 0.2", "μ = 3, σ = 1"}, "location", "northeast")
%! title ("Inverse Gaussian PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-Inf, -1, 0, 1/2, 1, Inf];
%! y = [0, 0, 0, 0.8788, 0.3989, 0];
%!assert (invgpdf ([x, NaN], 1, 1), [y, NaN], 1e-4)
%!assert (invgpdf (x, 1, [-2, -1, 0, 1, 1, 1]), [nan(1,3), y([4:6])], 1e-4)

## Test class of input preserved
%!assert (class (hncdf (single ([x, NaN]), 1, 1)), "single")
%!assert (class (hncdf ([x, NaN], 1, single (1))), "single")
%!assert (class (hncdf ([x, NaN], single (1), 1)), "single")

## Test input validation
%!error<invgpdf: function called with too few input arguments.> invgpdf ()
%!error<invgpdf: function called with too few input arguments.> invgpdf (1)
%!error<invgpdf: function called with too few input arguments.> invgpdf (1, 2)
%!error<invgpdf: X, MU, and LAMBDA must be of common size or scalars.> ...
%! invgpdf (1, ones (2), ones (3))
%!error<invgpdf: X, MU, and LAMBDA must be of common size or scalars.> ...
%! invgpdf (ones (2), 1, ones (3))
%!error<invgpdf: X, MU, and LAMBDA must be of common size or scalars.> ...
%! invgpdf (ones (2), ones (3), 1)
%!error<invgpdf: X, MU, and LAMBDA must not be complex.> invgpdf (i, 2, 3)
%!error<invgpdf: X, MU, and LAMBDA must not be complex.> invgpdf (1, i, 3)
%!error<invgpdf: X, MU, and LAMBDA must not be complex.> invgpdf (1, 2, i)
