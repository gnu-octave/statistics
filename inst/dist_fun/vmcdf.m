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
## @deftypefn  {statistics} {@var{p} =} vmcdf (@var{x}, @var{mu}, @var{k})
## @deftypefnx {statistics} {@var{p} =} vmcdf (@var{x}, @var{mu}, @var{k}, @qcode{"upper"})
##
## Von Mises probability density function (PDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the von Mises distribution with location parameter @var{mu} and
## concentration parameter @var{k} on the interval @math{[-pi,pi]}.  The size of
## @var{p} is the common size of @var{x}, @var{mu}, and @var{k}.  A scalar input
## functions as a constant matrix of the same same size as the other inputs.
##
## @code{@var{p} = vmcdf (@var{x}, @var{mu}, @var{k}, "upper")} computes the
## upper tail probability of the von Mises distribution with parameters @var{mu}
## and @var{k}, at the values in @var{x}.
##
## Note: the CDF of the von Mises distribution is not analytic.  Hence, it is
## calculated by integrating its probability density which is expressed as a
## series of Bessel functions.  Balancing between performance and accuracy, the
## integration uses a step of @qcode{1e-5} on the interval @math{[-pi,pi]},
## which results to an accuracy of about 10 significant digits.
##
## Further information about the von Mises distribution can be found at
## @url{https://en.wikipedia.org/wiki/Von_Mises_distribution}
##
## @seealso{vminv, vmpdf, vmrnd}
## @end deftypefn

function p = vmcdf (x, mu, k, uflag)

  ## Check for valid number of input arguments
  if (nargin <  3)
    error ("vmcdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("vmcdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, MU, and K
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (k))
    [retval, x, mu, k] = common_size (x, mu, k);
    if (retval > 0)
      error ("vmcdf: X, MU, and K must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and K being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (k))
    error ("vmcdf: X, MU, and K must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (k, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Evaluate Von Mises CDF by integrating from -PI to PI
  interval = linspace (-pi, pi, 1e5)'; # accurate to >10 significant digits
  f = exp (k .* cos (interval)) ./ (2 .* pi .* besseli (0, k));
  c = cumtrapz (interval, f);
  p = diag (interp1 (interval, c, x - mu, "spline"))';

  ## Force Nan for negative K
  p(k < 0) = NaN;

  ## Apply upper flag (if required)
  if (uflag)
    p = 1 - p;
  endif

endfunction

%!demo
%! ## Plot various CDFs from the von Mises distribution
%! x1 = [-pi:0.1:pi];
%! p1 = vmcdf (x1, 0, 0.5);
%! p2 = vmcdf (x1, 0, 1);
%! p3 = vmcdf (x1, 0, 2);
%! p4 = vmcdf (x1, 0, 4);
%! plot (x1, p1, "-r", x1, p2, "-g", x1, p3, "-b", x1, p4, "-c")
%! grid on
%! xlim ([-pi, pi])
%! legend ({"μ = 0, k = 0.5", "μ = 0, k = 1", ...
%!          "μ = 0, k = 2", "μ = 0, k = 4"}, "location", "northwest")
%! title ("Von Mises CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, p0, p1
%! x = [-pi:pi/2:pi];
%! p0 = [0, 0.10975, 0.5, 0.89025, 1];
%! p1 = [0, 0.03752, 0.5, 0.99622, 1];
%!assert (vmcdf (x, 0, 1), p0, 1e-5)
%!assert (vmcdf (x, 0, 1, "upper"), 1 - p0, 1e-5)
%!assert (vmcdf (x, zeros (1,5), ones (1,5)), p0, 1e-5)
%!assert (vmcdf (x, zeros (1,5), ones (1,5), "upper"), 1 - p0, 1e-5)
%!assert (vmcdf (x, 0, [1 2 3 4 5]), p1, 1e-5)
%!assert (vmcdf (x, 0, [1 2 3 4 5], "upper"), 1 - p1, 1e-5)

## Test class of input preserved
%!assert (isa (vmcdf (single (pi), 0, 1), "single"), true)
%!assert (isa (vmcdf (pi, single (0), 1), "single"), true)
%!assert (isa (vmcdf (pi, 0, single (1)), "single"), true)

## Test input validation
%!error<vmcdf: function called with too few input arguments.> vmcdf ()
%!error<vmcdf: function called with too few input arguments.> vmcdf (1)
%!error<vmcdf: function called with too few input arguments.> vmcdf (1, 2)
%!error<vmcdf: invalid argument for upper tail.> vmcdf (1, 2, 3, "tail")
%!error<vmcdf: invalid argument for upper tail.> vmcdf (1, 2, 3, 4)
%!error<vmcdf: X, MU, and K must be of common size or scalars.> ...
%! vmcdf (ones (3), ones (2), ones (2))
%!error<vmcdf: X, MU, and K must be of common size or scalars.> ...
%! vmcdf (ones (2), ones (3), ones (2))
%!error<vmcdf: X, MU, and K must be of common size or scalars.> ...
%! vmcdf (ones (2), ones (2), ones (3))
%!error<vmcdf: X, MU, and K must not be complex.> vmcdf (i, 2, 2)
%!error<vmcdf: X, MU, and K must not be complex.> vmcdf (2, i, 2)
%!error<vmcdf: X, MU, and K must not be complex.> vmcdf (2, 2, i)
