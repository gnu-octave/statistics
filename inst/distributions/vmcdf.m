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
## @deftypefn  {statistics} {@var{p} =} vmcdf (@var{theta})
## @deftypefnx {statistics} {@var{p} =} vmcdf (@var{theta}, @var{mu})
## @deftypefnx {statistics} {@var{p} =} vmcdf (@var{theta}, @var{mu}, @var{k})
##
## Von Mises probability density function (PDF).
##
## For each element of @var{theta}, compute the probability density function
## (PDF) at @var{theta} of the von Mises distribution with mean direction
## parameter @var{mu} and concentration parameter @var{k} on the interval
## [-pi, pi].  The size of @var{p} is the common size of the input arguments.
## A scalar input functions as a constant matrix of the same same size as the
## other inputs.
##
## Default values are @var{mu} = 0, @var{k} = 1.  The function returns NaN for
## negative @var{k}.
##
## @seealso{vmpdf, vmrnd}
## @end deftypefn

function p = vmcdf (theta, mu = 0, k = 1)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  if (! isscalar (theta) || ! isscalar (mu) || ! isscalar (k))
    [retval, theta, mu, k] = common_size (theta, mu, k);
    if (retval > 0)
      error ("vmpdf: THETA, MU, and K must be of common size or scalars.");
    endif
  endif

  if (iscomplex (theta) || iscomplex (mu) || iscomplex (k))
    error ("vmpdf: THETA, MU, and K must not be complex.");
  endif

  if (isa (theta, "single") || isa (mu, "single") || isa (k, "single"))
    p = zeros (size (theta), "single");
  else
    p = zeros (size (theta));
  endif

  ## Evaluate Von Mises CDF by integrating from -PI to PI
  interval = linspace (-pi, pi, 1e6)';
  f = exp (k .* cos (interval - mu)) ./ (2 .* pi .* besseli (0, k));
  c = cumtrapz (interval, f);
  p = diag (interp1 (interval, c, theta, "spline"))';

  ## Force Nan for negative K
  p(k < 0) = NaN;

endfunction

%!shared theta, p0, p1
%! theta = [-pi:pi/2:pi];
%! p0 = [0, 0.10975, 0.5, 0.89025, 1];
%! p1 = [0, 0.03752, 0.5, 0.99622, 1];
%!assert (vmcdf (theta), p0, 1e-5)
%!assert (vmcdf (theta, zeros (1,5), ones (1,5)), p0, 1e-5)
%!assert (vmcdf (theta, 0, [1 2 3 4 5]), p1, 1e-5)

## Test class of input preserved
%!assert (isa (vmcdf (single (pi), 0, 1), "single"), true)
%!assert (isa (vmcdf (pi, single (0), 1), "single"), true)
%!assert (isa (vmcdf (pi, 0, single (1)), "single"), true)

## Test input validation
%!error vmcdf ()
%!error vmcdf (1, 2, 3, 4)
%!error vmcdf (ones (3), ones (2), ones (2))
%!error vmcdf (ones (2), ones (3), ones (2))
%!error vmcdf (ones (2), ones (2), ones (3))
%!error vmcdf (i, 2, 2)
%!error vmcdf (2, i, 2)
%!error vmcdf (2, 2, i)

