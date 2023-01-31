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
## @deftypefn  {statistics} @var{x} = vmcdf (@var{theta})
## @deftypefnx {statistics} @var{x} = vmcdf (@var{theta}, @var{mu})
## @deftypefnx {statistics} @var{x} = vmcdf (@var{theta}, @var{mu}, @var{k})
##
## Von Mises probability density function (PDF).
##
## For each element of @var{theta}, compute the probability density function
## (PDF) at @var{theta} of the von Mises distribution with parameters @var{mu}
## and @var{k} on the interval [-pi, pi}].  The size of @var{x} is the common
## size of the input arguments.  A scalar input functions as a constant matrix
## of the same same size as the other inputs.  @var{k} must be a real positive
## scalar.
##
## Default values are @var{mu} = 0, @var{k} = 1.
##
## @seealso{vmpdf, vmrnd}
## @end deftypefn

function p = vmcdf (theta, mu = 0, k = 1)

  if (nargin < 1 && nargin > 3)
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
    y = zeros (size (theta), "single");
  else
    y = zeros (size (theta));
  endif

  ## Evaluate Von Mises CDF by integrating from -PI to PI
  interval = linspace (-pi, pi, 1e6);
  f = exp (k .* cos (interval - mu)) ./ (2 .* pi .* besseli (0, k));
  c = cumtrapz (interval, f);
  p = interp1 (interval, c, theta, "spline");

  ## Force Nan for negative K
  p(k <= 0) = NaN;

endfunction
