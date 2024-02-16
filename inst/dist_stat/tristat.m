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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} tristat (@var{a}, @var{b}, @var{c})
##
## Compute statistics of the Triangular distribution.
##
## @code{[@var{m}, @var{v}] = tristat (@var{a}, @var{b}, @var{c})} returns the
## mean and variance of the Triangular distribution with lower limit parameter
## @var{a}, peak location (mode) parameter @var{b}, and upper limit parameter
## @var{c}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Note that the order of the parameter input arguments has been changed after
## statistics version 1.6.3 in order to be MATLAB compatible with the parameters
## used in the TriangularDistribution probability distribution object.  More
## specifically, the positions of the parameters @var{b} and @var{c} have been
## swapped.  As a result, the naming conventions no longer coinside with those
## used in Wikipedia, in which @math{b} denotes the upper limit and @math{c}
## denotes the mode or peak parameter.
##
## Further information about the triangular distribution can be found at
## @url{https://en.wikipedia.org/wiki/Triangular_distribution}
##
## @seealso{tcdf, tinv, tpdf, trnd}
## @end deftypefn

function [m, v] = tristat (a, b, c)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("tristat: function called with too few input arguments.");
  endif

  ## Check for A, B, and C being numeric
  if (! (isnumeric (a) && isnumeric (b) && isnumeric (c)))
    error ("tristat: A, B, and C must be numeric.");
  endif

  ## Check for A, B, and C being real
  if (! (isreal (a) && isreal (b) && isreal (c)))
    error ("tristat: A, B, and C must be real.");
  endif

  ## Calculate moments
  m = (a + b + c) ./ 3;
  v = (a .^ 2 + b .^ 2 + c .^ 2 - a .* b - a .* c - b .* c) ./ 18;

  ## Continue argument check
  k = find (! (a < c) | ! (a <= b & b <= c));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<tristat: function called with too few input arguments.> tristat ()
%!error<tristat: function called with too few input arguments.> tristat (1)
%!error<tristat: function called with too few input arguments.> tristat (1, 2)
%!error<tristat: A, B, and C must be numeric.> tristat ("i", 2, 1)
%!error<tristat: A, B, and C must be numeric.> tristat (0, "d", 1)
%!error<tristat: A, B, and C must be numeric.> tristat (0, 3, {})
%!error<tristat: A, B, and C must be real.> tristat (i, 2, 1)
%!error<tristat: A, B, and C must be real.> tristat (0, i, 1)
%!error<tristat: A, B, and C must be real.> tristat (0, 3, i)

## Output validation tests
%!test
%! a = 1:5;
%! b = 3:7;
%! c = 5:9;
%! [m, v] = tristat (a, b, c);
%! expected_m = [3, 4, 5, 6, 7];
%! assert (m, expected_m);
%! assert (v, ones (1, 5) * (2/3));
