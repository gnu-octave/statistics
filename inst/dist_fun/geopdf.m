## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{y} =} geopdf (@var{x}, @var{ps})
##
## Geometric probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the geometric distribution with probability of success parameter @var{ps}.
## The size of @var{y} is the common size of @var{x} and @var{ps}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## The geometric distribution models the number of failures (@var{x}) of a
## Bernoulli trial with probability @var{ps} before the first success.
##
## Further information about the geometric distribution can be found at
## @url{https://en.wikipedia.org/wiki/Geometric_distribution}
##
## @seealso{geocdf, geoinv, geornd, geofit, geostat}
## @end deftypefn

function y = geopdf (x, ps)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("geopdf: function called with too few input arguments.");
  endif

  ## Check for common size of X and PS
  if (! isscalar (x) || ! isscalar (ps))
    [retval, x, ps] = common_size (x, ps);
    if (retval > 0)
      error ("geopdf: X and PS must be of common size or scalars.");
    endif
  endif

  ## Check for X and PS being reals
  if (iscomplex (x) || iscomplex (ps))
    error ("geopdf: X and PS must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (ps, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  ## Return NaN for out of range parameters
  k = isnan (x) | (x == Inf) | !(ps >= 0) | !(ps <= 1);
  y(k) = NaN;

  ## Get valid instances
  k = (x >= 0) & (x < Inf) & (x == fix (x)) & (ps > 0) & (ps <= 1);

  ## Compute CDF
  if (isscalar (ps))
    y(k) = ps * ((1 - ps) .^ x(k));
  else
    y(k) = ps(k) .* ((1 - ps(k)) .^ x(k));
  endif

endfunction

%!demo
%! ## Plot various PDFs from the geometric distribution
%! x = 0:10;
%! y1 = geopdf (x, 0.2);
%! y2 = geopdf (x, 0.5);
%! y3 = geopdf (x, 0.7);
%! plot (x, y1, "*b", x, y2, "*g", x, y3, "*r")
%! grid on
%! ylim ([0, 0.8])
%! legend ({"ps = 0.2", "ps = 0.5", "ps = 0.7"}, "location", "northeast")
%! title ("Geometric PDF")
%! xlabel ("values in x (number of failures)")
%! ylabel ("density")

## Test output
%!shared x, y
%! x = [-1 0 1 Inf];
%! y = [0, 1/2, 1/4, NaN];
%!assert (geopdf (x, 0.5*ones (1,4)), y)
%!assert (geopdf (x, 0.5), y)
%!assert (geopdf (x, 0.5*[-1 NaN 4 1]), [NaN NaN NaN y(4)])
%!assert (geopdf ([x, NaN], 0.5), [y, NaN])

## Test class of input preserved
%!assert (geopdf (single ([x, NaN]), 0.5), single ([y, NaN]), 5*eps ("single"))
%!assert (geopdf ([x, NaN], single (0.5)), single ([y, NaN]), 5*eps ("single"))

## Test input validation
%!error geopdf ()
%!error geopdf (1)
%!error geopdf (1,2,3)
%!error geopdf (ones (3), ones (2))
%!error geopdf (ones (2), ones (3))
%!error geopdf (i, 2)
%!error geopdf (2, i)
