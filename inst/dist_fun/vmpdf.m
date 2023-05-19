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
## @deftypefn  {statistics} {@var{y} =} vmpdf (@var{x}, @var{mu}, @var{k})
##
## Von Mises probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the von Mises distribution with location parameter @var{mu} and
## concentration parameter @var{k} on the interval [-pi, pi].  The size of
## @var{y} is the common size of @var{x}, @var{mu}, and @var{k}. A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## Further information about the von Mises distribution can be found at
## @url{https://en.wikipedia.org/wiki/Von_Mises_distribution}
##
## @seealso{vmcdf, vminv, vmrnd}
## @end deftypefn

function y = vmpdf (x, mu, k)

  ## Check for valid number of input arguments
  if (nargin <  3)
    error ("vmpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, MU, and K
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (k))
    [retval, x, mu, k] = common_size (x, mu, k);
    if (retval > 0)
      error ("vmpdf: X, MU, and K must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and K being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (k))
    error ("vmpdf: X, MU, and K must not be complex.");
  endif

  ## Evaluate Von Mises PDF
  Z = 2 .* pi .* besseli (0, k);
  y = exp (k .* cos (x - mu)) ./ Z;

  ## Force Nan for negative K
  y(k < 0) = NaN;

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (k, "single"))
    y = cast (y, "single");
  else
    y = cast (y, "double");
  endif

endfunction

%!demo
%! ## Plot various PDFs from the von Mises distribution
%! x1 = [-pi:0.1:pi];
%! y1 = vmpdf (x1, 0, 0.5);
%! y2 = vmpdf (x1, 0, 1);
%! y3 = vmpdf (x1, 0, 2);
%! y4 = vmpdf (x1, 0, 4);
%! plot (x1, y1, "-r", x1, y2, "-g", x1, y3, "-b", x1, y4, "-c")
%! grid on
%! xlim ([-pi, pi])
%! ylim ([0, 0.8])
%! legend ({"μ = 0, k = 0.5", "μ = 0, k = 1", ...
%!          "μ = 0, k = 2", "μ = 0, k = 4"}, "location", "northwest")
%! title ("Von Mises PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!shared x, y0, y1
%! x = [-pi:pi/2:pi];
%! y0 = [0.046245, 0.125708, 0.341710, 0.125708, 0.046245];
%! y1 = [0.046245, 0.069817, 0.654958, 0.014082, 0.000039];
%!assert (vmpdf (x, 0, 1), y0, 1e-5)
%!assert (vmpdf (x, zeros (1,5), ones (1,5)), y0, 1e-6)
%!assert (vmpdf (x, 0, [1 2 3 4 5]), y1, 1e-6)

## Test class of input preserved
%!assert (isa (vmpdf (single (pi), 0, 1), "single"), true)
%!assert (isa (vmpdf (pi, single (0), 1), "single"), true)
%!assert (isa (vmpdf (pi, 0, single (1)), "single"), true)

## Test input validation
%!error<vmpdf: function called with too few input arguments.> vmpdf ()
%!error<vmpdf: function called with too few input arguments.> vmpdf (1)
%!error<vmpdf: function called with too few input arguments.> vmpdf (1, 2)
%!error<vmpdf: X, MU, and K must be of common size or scalars.> ...
%! vmpdf (ones (3), ones (2), ones (2))
%!error<vmpdf: X, MU, and K must be of common size or scalars.> ...
%! vmpdf (ones (2), ones (3), ones (2))
%!error<vmpdf: X, MU, and K must be of common size or scalars.> ...
%! vmpdf (ones (2), ones (2), ones (3))
%!error<vmpdf: X, MU, and K must not be complex.> vmpdf (i, 2, 2)
%!error<vmpdf: X, MU, and K must not be complex.> vmpdf (2, i, 2)
%!error<vmpdf: X, MU, and K must not be complex.> vmpdf (2, 2, i)
