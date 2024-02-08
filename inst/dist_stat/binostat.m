## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
## Copyright (C) 2015 Carnë Draug <carandraug@octave.org>
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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} binostat (@var{n}, @var{ps})
##
## Compute statistics of the binomial distribution.
##
## @code{[@var{m}, @var{v}] = binostat (@var{n}, @var{ps})} returns the mean and
## variance of the binomial distribution with parameters @var{n} and @var{ps},
## where @var{n} is the number of trials and @var{ps} is the probability of
## success.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the
## same size as the other inputs.
##
## Further information about the binomial distribution can be found at
## @url{https://en.wikipedia.org/wiki/Binomial_distribution}
##
## @seealso{binocdf, binoinv, binopdf, binornd, binofit, binolike, binotest}
## @end deftypefn

function [m, v] = binostat (n, ps)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("binostat: function called with too few input arguments.");
  endif

  ## Check for N and PS being numeric
  if (! (isnumeric (n) && isnumeric (ps)))
    error ("binostat: N and PS must be numeric.");
  endif

  ## Check for N and PS being real
  if (iscomplex (n) || iscomplex (ps))
    error ("binostat: N and PS must not be complex.");
  endif

  ## Check for common size of N and PS
  if (! isscalar (n) || ! isscalar (ps))
    [retval, n, ps] = common_size (n, ps);
    if (retval > 0)
      error ("binostat: N and PS must be of common size or scalars.");
    endif
  endif

  ## Catch invalid parameters
  k = find (! (n > 0 & fix (n) == n & ps >= 0 & ps <= 1));

  ## Calculate moments
  m = n .* ps;
  m(k) = NaN;

  if (nargout > 1)
    v = m .* (1 - ps);
    v(k) = NaN;
  endif

endfunction

## Input validation tests
%!error<binostat: function called with too few input arguments.> binostat ()
%!error<binostat: function called with too few input arguments.> binostat (1)
%!error<binostat: N and PS must be numeric.> binostat ({}, 2)
%!error<binostat: N and PS must be numeric.> binostat (1, "")
%!error<binostat: N and PS must not be complex.> binostat (i, 2)
%!error<binostat: N and PS must not be complex.> binostat (1, i)
%!error<binostat: N and PS must be of common size or scalars.> ...
%! binostat (ones (3), ones (2))
%!error<binostat: N and PS must be of common size or scalars.> ...
%! binostat (ones (2), ones (3))

## Output validation tests
%!test
%! n = 1:6;
%! ps = 0:0.2:1;
%! [m, v] = binostat (n, ps);
%! expected_m = [0.00, 0.40, 1.20, 2.40, 4.00, 6.00];
%! expected_v = [0.00, 0.32, 0.72, 0.96, 0.80, 0.00];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
%!test
%! n = 1:6;
%! [m, v] = binostat (n, 0.5);
%! expected_m = [0.50, 1.00, 1.50, 2.00, 2.50, 3.00];
%! expected_v = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50];
%! assert (m, expected_m, 0.001);
%! assert (v, expected_v, 0.001);
%!test
%! n = [-Inf -3 5 0.5 3 NaN 100, Inf];
%! [m, v] = binostat (n, 0.5);
%! assert (isnan (m), [true true false true false true false false])
%! assert (isnan (v), [true true false true false true false false])
%! assert (m(end), Inf);
%! assert (v(end), Inf);
%!assert (nthargout (1:2, @binostat, 5, []), {[], []})
%!assert (nthargout (1:2, @binostat, [], 5), {[], []})
%!assert (size (binostat (randi (100, 10, 5, 4), rand (10, 5, 4))), [10 5 4])
%!assert (size (binostat (randi (100, 10, 5, 4), 7)), [10 5 4])
