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
## @deftypefn  {statistics} {[@var{m}, @var{v}] =} burrstat (@var{lambda}, @var{c}, @var{k})
##
## Compute statistics of the Burr type XII distribution.
##
## @code{[@var{m}, @var{v}] = burrstat (@var{lambda}, @var{c}, @var{k})} returns
## the mean and variance of the Burr type XII distribution with scale parameter
## @var{lambda}, first shape parameter @var{c}, and second shape parameter
## @var{k}.
##
## The size of @var{m} (mean) and @var{v} (variance) is the common size of the
## input arguments.  A scalar input functions as a constant matrix of the same
## size as the other inputs.
##
## Further information about the Burr distribution can be found at
## @url{https://en.wikipedia.org/wiki/Burr_distribution}
##
## @seealso{gevcdf, gevinv, gevpdf, gevrnd, gevfit, gevlike}
## @end deftypefn

function [m, v] = burrstat (lambda, c, k)

  ## Check for is_val number of input arguments
  if (nargin < 3)
    error ("burrstat: function called with too few input arguments.");
  endif

  ## Check for LAMBDA, C, and K being numeric
  if (! (isnumeric (lambda) && isnumeric (c) && isnumeric (k)))
    error ("burrstat: LAMBDA, C, and K must be numeric.");
  endif

  ## Check for LAMBDA, C, and K being real
  if (iscomplex (lambda) || iscomplex (c) || iscomplex (k))
    error ("burrstat: LAMBDA, C, and K must not be complex.");
  endif

  ## Check for common size of LAMBDA, C, and K
  if (! isscalar (lambda) || ! isscalar (c) || ! isscalar (k))
    [retval, lambda, c, k] = common_size (lambda, c, k);
    if (retval > 0)
      error ("burrstat: LAMBDA, C, and K must be of common size or scalars.");
    endif
  endif

  ## Preallocate mean annd variance
  m = v = nan (size (lambda));
  ## Precalculate some values
  c_i = 1 ./ c;
  l_c = lambda .* c_i;
  kci = k - c_i;

  ## Find valid vases
  c_k = c .* k;
  is_val = lambda > 0 & c > 0 & k > 0;

  ## Calculate 1st moment
  is_inf = is_val & c_k <= 1;
  m(is_inf) = Inf;
  no_inf = ! is_inf;
  m(no_inf) = l_c(no_inf) .* beta (c_i(no_inf), kci(no_inf));

  ## Calculate 2nd moment
  is_inf = is_val & c_k <= 2;
  v(is_inf) = Inf;
  no_inf = ! is_inf;
  v(no_inf) = 2 * lambda(no_inf) .* l_c(no_inf) .* ...
              beta (c_i(no_inf) .* 2, kci(no_inf) - c_i(no_inf)) ...
              - m(no_inf) .^ 2;

endfunction

## Input validation tests
%!error<burrstat: function called with too few input arguments.> burrstat ()
%!error<burrstat: function called with too few input arguments.> burrstat (1)
%!error<burrstat: function called with too few input arguments.> burrstat (1, 2)
%!error<burrstat: LAMBDA, C, and K must be numeric.> burrstat ({}, 2, 3)
%!error<burrstat: LAMBDA, C, and K must be numeric.> burrstat (1, "", 3)
%!error<burrstat: LAMBDA, C, and K must be numeric.> burrstat (1, 2, "")
%!error<burrstat: LAMBDA, C, and K must not be complex.> burrstat (i, 2, 3)
%!error<burrstat: LAMBDA, C, and K must not be complex.> burrstat (1, i, 3)
%!error<burrstat: LAMBDA, C, and K must not be complex.> burrstat (1, 2, i)
%!error<burrstat: LAMBDA, C, and K must be of common size or scalars.> ...
%! burrstat (ones (3), ones (2), 3)
%!error<burrstat: LAMBDA, C, and K must be of common size or scalars.> ...
%! burrstat (ones (2), 2, ones (3))
%!error<burrstat: LAMBDA, C, and K must be of common size or scalars.> ...
%! burrstat (1, ones (2), ones (3))

## Output validation tests
%!test
%! [m, v] = burrstat (1, 2, 5);
%! assert (m, 0.4295, 1e-4);
%! assert (v, 0.0655, 1e-4);
%!test
%! [m, v] = burrstat (1, 1, 1);
%! assert (m, Inf);
%! assert (v, Inf);
%!test
%! [m, v] = burrstat (2, 4, 1);
%! assert (m, 2.2214, 1e-4);
%! assert (v, 1.3484, 1e-4);
