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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} betastat (@var{a}, @var{b})
##
## Compute statistics of the Beta distribution.
##
## @code{[@var{m}, @var{v}] = betastat (@var{a}, @var{b})} returns the mean
## and variance of the Beta distribution with shape parameters @var{a} and
## @var{b}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the Beta distribution can be found at
## @url{https://en.wikipedia.org/wiki/Beta_distribution}
##
## @seealso{betacdf, betainv, betapdf, betarnd, betafit, betalike}
## @end deftypefn

function [m, v] = betastat (a, b)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("betastat: function called with too few input arguments.");
  endif

  ## Check for A and B being numeric
  if (! (isnumeric (a) && isnumeric (b)))
    error ("betastat: A and B must be numeric.");
  endif

  ## Check for A and B being real
  if (iscomplex (a) || iscomplex (b))
    error ("betastat: A and B must not be complex.");
  endif

  ## Check for common size of A and B
  if (! isscalar (a) || ! isscalar (b))
    [retval, a, b] = common_size (a, b);
    if (retval > 0)
      error ("betastat: A and B must be of common size or scalars.");
    endif
  endif

  ## Catch invalid parameters
  k = find (! (a > 0 & b > 0));

  ## Calculate moments
  a_b = a + b;
  m = a ./ (a_b);
  m(k) = NaN;

  if (nargout > 1)
    v = (a .* b) ./ ((a_b .^ 2) .* (a_b + 1));
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<betastat: function called with too few input arguments.> betastat ()
%!error<betastat: function called with too few input arguments.> betastat (1)
%!error<betastat: A and B must be numeric.> betastat ({}, 2)
%!error<betastat: A and B must be numeric.> betastat (1, "")
%!error<betastat: A and B must not be complex.> betastat (i, 2)
%!error<betastat: A and B must not be complex.> betastat (1, i)
%!error<betastat: A and B must be of common size or scalars.> ...
%! betastat (ones (3), ones (2))
%!error<betastat: A and B must be of common size or scalars.> ...
%! betastat (ones (2), ones (3))

## Output validation tests
%!test
%! a = -2:6;
%! b = 0.4:0.2:2;
%! [m, v] = betastat (a, b);
%! expected_m = [NaN NaN NaN 1/2 2/3.2 3/4.4 4/5.6 5/6.8 6/8];
%! expected_v = [NaN NaN NaN 0.0833, 0.0558, 0.0402, 0.0309, 0.0250, 0.0208];
%! assert (m, expected_m, eps*100);
%! assert (v, expected_v, 0.001);
%!test
%! a = -2:1:6;
%! [m, v] = betastat (a, 1.5);
%! expected_m = [NaN NaN NaN 1/2.5 2/3.5 3/4.5 4/5.5 5/6.5 6/7.5];
%! expected_v = [NaN NaN NaN 0.0686, 0.0544, 0.0404, 0.0305, 0.0237, 0.0188];
%! assert (m, expected_m);
%! assert (v, expected_v, 0.001);
%!test
%! a = [14  Inf   10  NaN  10];
%! b = [12    9  NaN  Inf  12];
%! [m, v] = betastat (a, b);
%! expected_m = [14/26 NaN NaN NaN 10/22];
%! expected_v = [168/18252 NaN NaN NaN 120/11132];
%! assert (m, expected_m);
%! assert (v, expected_v);
%!assert (nthargout (1:2, @betastat, 5, []), {[], []})
%!assert (nthargout (1:2, @betastat, [], 5), {[], []})
%!assert (size (betastat (rand (10, 5, 4), rand (10, 5, 4))), [10 5 4])
%!assert (size (betastat (rand (10, 5, 4), 7)), [10 5 4])

