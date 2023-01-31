## Copyright (C) 2009 Soren Hauberg <soren@hauberg.org>
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
## @deftypefn  {statistics} @var{y} = vmpdf (@var{theta})
## @deftypefnx {statistics} @var{y} = vmpdf (@var{theta}, @var{mu})
## @deftypefnx {statistics} @var{y} = vmpdf (@var{theta}, @var{mu}, @var{k})
##
## Von Mises probability density function (PDF).
##
## For each element of @var{theta}, compute the probability density function
## (PDF) at @var{theta} of the von Mises distribution with mean direction
## parameter @var{mu} and concentration parameter @var{k} on the interval
## [-pi, pi].  The size of @var{y} is the common size of the input arguments.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## Default values are @var{mu} = 0, @var{k} = 1.  The function returns NaN for
## negative @var{k}.
##
## @seealso{vmcdf, vmrnd}
## @end deftypefn

function y = vmpdf (theta, mu = 0, k = 1)

  if (nargin < 1 && nargin > 3)
    print_usage ();
  endif

  if (! isscalar (theta) || ! isscalar (mu) || ! isscalar (k))
    [retval, theta, mu, k] = common_size (theta, mu, k);
    if (retval > 0)
      error ("vmpdf: X, MU, and K must be of common size or scalars.");
    endif
  endif

  if (iscomplex (theta) || iscomplex (mu) || iscomplex (k))
    error ("vmpdf: X, MU, and K must not be complex.");
  endif

  ## Evaluate Von Mises PDF
  Z = 2 .* pi .* besseli (0, k);
  y = exp (k .* cos (theta - mu)) ./ Z;

  ## Force Nan for negative K
  y(k < 0) = NaN;

  ## Preserve class
  if (isa (theta, "single") || isa (mu, "single") || isa (k, "single"))
    y = cast (y, "single");
  else
    y = cast (y, "double");
  endif

endfunction

%!shared theta, p0, p1
%! theta = [-pi:pi/2:pi];
%! p0 = [0.046245, 0.125708, 0.341710, 0.125708, 0.046245];
%! p1 = [0.046245, 0.069817, 0.654958, 0.014082, 0.000039];
%!assert (vmpdf (theta), p0, 1e-5)
%!assert (vmpdf (theta, zeros (1,5), ones (1,5)), p0, 1e-6)
%!assert (vmpdf (theta, 0, [1 2 3 4 5]), p1, 1e-6)

## Test class of input preserved
%!assert (isa (vmpdf (single (pi), 0, 1), "single"), true)
%!assert (isa (vmpdf (pi, single (0), 1), "single"), true)
%!assert (isa (vmpdf (pi, 0, single (1)), "single"), true)

## Test input validation
%!error vmcdf ()
%!error vmcdf (1, 2, 3, 4)
%!error vmcdf (ones (3), ones (2), ones (2))
%!error vmcdf (ones (2), ones (3), ones (2))
%!error vmcdf (ones (2), ones (2), ones (3))
%!error vmcdf (i, 2, 2)
%!error vmcdf (2, i, 2)
%!error vmcdf (2, 2, i)
