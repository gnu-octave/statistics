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
## @deftypefn  {statistics} @var{theta} = vmpdf (@var{x})
## @deftypefnx {statistics} @var{theta} = vmpdf (@var{x}, @var{mu})
## @deftypefnx {statistics} @var{theta} = vmpdf (@var{x}, @var{mu}, @var{k})
##
## Von Mises probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the von Mises distribution with parameters @var{mu} and @var{k}
## on the interval [-pi, pi}].  The size of @var{theta} is the common size of
## the input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.  @var{k} must be a real positive scalar.
##
## Default values are @var{mu} = 0, @var{k} = 1.
##
## @seealso{vmcdf, vmrnd}
## @end deftypefn

function theta = vmpdf (x, mu = 0, k = 1)

  if (nargin < 1 && nargin > 3)
    print_usage ();
  endif

  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ("vmpdf: X, MU, and K must be of common size or scalars.");
    endif
  endif

  if (iscomplex (x) || iscomplex (a) || iscomplex (b))
    error ("vmpdf: X, MU, and K must not be complex.");
  endif

  ## Evaluate Von Mises PDF
  Z = 2 .* pi .* besseli (0, k);
  theta = exp (k .* cos (x - mu)) ./ Z;

  ## Force Nan for negative K
  theta(k <= 0) = NaN;

  ## Preserve class
  if (isa (x, "single") || isa (a, "single") || isa (b, "single"))
    theta = cast (theta, "single");
  else
    theta = cast (theta, "double");
  endif

endfunction


