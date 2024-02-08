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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} unifstat (@var{df})
##
## Compute statistics of the continuous uniform cumulative distribution.
##
## @code{[@var{m}, @var{v}] = unifstat (@var{df})} returns the mean and variance
## of the continuous uniform cumulative distribution with parameters @var{a} and
## @var{b}, which define the lower and upper bounds of the interval
## @qcode{[@var{a}, @var{b}]}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the continuous uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Continuous_uniform_distribution}
##
## @seealso{unifcdf, unifinv, unifpdf, unifrnd, unifit}
## @end deftypefn

function [m, v] = unifstat (a, b)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("unifstat: function called with too few input arguments.");
  endif

  ## Check for A and B being numeric
  if (! (isnumeric (a) && isnumeric (b)))
    error ("unifstat: A and B must be numeric.");
  endif

  ## Check for A and B being real
  if (iscomplex (a) || iscomplex (b))
    error ("unifstat: A and B must not be complex.");
  endif

  ## Check for common size of A and B
  if (! isscalar (a) || ! isscalar (b))
    [retval, a, b] = common_size (a, b);
    if (retval > 0)
      error ("unifstat: A and B must be of common size or scalars.");
    endif
  endif

  ## Calculate moments
  m = (a + b) ./ 2;
  v = ((b - a) .^ 2) ./ 12;

  ## Continue argument check
  k = find (! (-Inf < a) | ! (a < b) | ! (b < Inf));
  if (any (k))
    m(k) = NaN;
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<unifstat: function called with too few input arguments.> unifstat ()
%!error<unifstat: function called with too few input arguments.> unifstat (1)
%!error<unifstat: A and B must be numeric.> unifstat ({}, 2)
%!error<unifstat: A and B must be numeric.> unifstat (1, "")
%!error<unifstat: A and B must not be complex.> unifstat (i, 2)
%!error<unifstat: A and B must not be complex.> unifstat (1, i)
%!error<unifstat: A and B must be of common size or scalars.> ...
%! unifstat (ones (3), ones (2))
%!error<unifstat: A and B must be of common size or scalars.> ...
%! unifstat (ones (2), ones (3))

## Output validation tests
%!test
%! a = 1:6;
%! b = 2:2:12;
%! [m, v] = unifstat (a, b);
%! expected_m = [1.5000, 3.0000, 4.5000, 6.0000, 7.5000, 9.0000];
%! expected_v = [0.0833, 0.3333, 0.7500, 1.3333, 2.0833, 3.0000];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
%!test
%! a = 1:6;
%! [m, v] = unifstat (a, 10);
%! expected_m = [5.5000, 6.0000, 6.5000, 7.0000, 7.5000, 8.0000];
%! expected_v = [6.7500, 5.3333, 4.0833, 3.0000, 2.0833, 1.3333];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
