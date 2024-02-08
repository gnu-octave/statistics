## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} geostat (@var{ps})
##
## Compute statistics of the geometric distribution.
##
## @code{[@var{m}, @var{v}] = geostat (@var{ps})} returns the mean and
## variance of the geometric distribution with probability of success parameter
## @var{ps}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the same size of the
## input argument.
##
## Further information about the geometric distribution can be found at
## @url{https://en.wikipedia.org/wiki/Geometric_distribution}
##
## @seealso{geocdf, geoinv, geopdf, geornd, geofit}
## @end deftypefn

function [m, v] = geostat (ps)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("geostat: function called with too few input arguments.");
  endif

  ## Check for PS being numeric
  if (! isnumeric (ps))
    error ("geostat: PS must be numeric.");
  endif

  ## Check for PS being real
  if (iscomplex (ps))
    error ("geostat: PS must not be complex.");
  endif

  ## Calculate moments
  q = 1 - ps;
  m = q ./ ps;
  v = q ./ (ps .^ 2);

  ## Continue argument check
  k = find (! (ps >= 0) | ! (ps <= 1));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<geostat: function called with too few input arguments.> geostat ()
%!error<geostat: PS must be numeric.> geostat ({})
%!error<geostat: PS must be numeric.> geostat ("")
%!error<geostat: PS must not be complex.> geostat (i)

## Output validation tests
%!test
%! ps = 1 ./ (1:6);
%! [m, v] = geostat (ps);
%! assert (m, [0, 1, 2, 3, 4, 5], 0.001);
%! assert (v, [0, 2, 6, 12, 20, 30], 0.001);
