## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{x} =} unifinv (@var{p})
## @deftypefnx {statistics} {@var{x} =} unifcdf (@var{p}, @var{a})
## @deftypefnx {statistics} {@var{x} =} unifcdf (@var{p}, @var{a}, @var{b})
##
## Inverse of the continuous uniform cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the continuous uniform distribution on the interval @qcode{[@var{a},
## @var{b}]}.  The size of @var{x} is the common size of @var{p}, @var{a}, and
## @var{b}.  A scalar input functions as a constant matrix of the same size as
## the other inputs.
##
## Further information about the continuous uniform distribution can be found at
## @url{https://en.wikipedia.org/wiki/Continuous_uniform_distribution}
##
## @seealso{unifcdf, unifpdf, unifrnd, unifit, unifstat}
## @end deftypefn

function x = unifinv (p, a, b)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("unifinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, A, and B
  if (! isscalar (p) || ! isscalar (a) || ! isscalar (b))
    [retval, p, a, b] = common_size (p, a, b);
    if (retval > 0)
      error ("unifinv: P, A, and B must be of common size or scalars.");
    endif
  endif

  ## Check for P, A, and B being reals
  if (iscomplex (p) || iscomplex (a) || iscomplex (b))
    error ("unifinv: P, A, and B must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (a, "single") || isa (b, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Calculate continuous uniform iCDF for valid parameter and data range
  k = (p >= 0) & (p <= 1) & (a < b);
  x(k) = a(k) + p(k) .* (b(k) - a(k));

endfunction

%!demo
%! ## Plot various iCDFs from the continuous uniform distribution
%! p = 0.001:0.001:0.999;
%! x1 = unifinv (p, 2, 5);
%! x2 = unifinv (p, 3, 9);
%! plot (p, x1, "-b", p, x2, "-g")
%! grid on
%! xlim ([0, 1])
%! ylim ([0, 10])
%! legend ({"a = 2, b = 5", "a = 3, b = 9"}, "location", "northwest")
%! title ("Continuous uniform iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (unifinv (p, ones (1,5), 2*ones (1,5)), [NaN 1 1.5 2 NaN])
%!assert (unifinv (p, 0, 1), [NaN 1 1.5 2 NaN] - 1)
%!assert (unifinv (p, 1, 2*ones (1,5)), [NaN 1 1.5 2 NaN])
%!assert (unifinv (p, ones (1,5), 2), [NaN 1 1.5 2 NaN])
%!assert (unifinv (p, [1 2 NaN 1 1], 2), [NaN NaN NaN 2 NaN])
%!assert (unifinv (p, 1, 2*[1 0 NaN 1 1]), [NaN NaN NaN 2 NaN])
%!assert (unifinv ([p(1:2) NaN p(4:5)], 1, 2), [NaN 1 NaN 2 NaN])

## Test class of input preserved
%!assert (unifinv ([p, NaN], 1, 2), [NaN 1 1.5 2 NaN NaN])
%!assert (unifinv (single ([p, NaN]), 1, 2), single ([NaN 1 1.5 2 NaN NaN]))
%!assert (unifinv ([p, NaN], single (1), 2), single ([NaN 1 1.5 2 NaN NaN]))
%!assert (unifinv ([p, NaN], 1, single (2)), single ([NaN 1 1.5 2 NaN NaN]))

## Test input validation
%!error<unifinv: function called with too few input arguments.> unifinv ()
%!error<unifinv: function called with too few input arguments.> unifinv (1, 2)
%!error<unifinv: P, A, and B must be of common size or scalars.> ...
%! unifinv (ones (3), ones (2), ones (2))
%!error<unifinv: P, A, and B must be of common size or scalars.> ...
%! unifinv (ones (2), ones (3), ones (2))
%!error<unifinv: P, A, and B must be of common size or scalars.> ...
%! unifinv (ones (2), ones (2), ones (3))
%!error<unifinv: P, A, and B must not be complex.> unifinv (i, 2, 2)
%!error<unifinv: P, A, and B must not be complex.> unifinv (2, i, 2)
%!error<unifinv: P, A, and B must not be complex.> unifinv (2, 2, i)
