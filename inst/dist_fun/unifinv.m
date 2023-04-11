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
## Inverse of the uniform cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the uniform distribution on the interval [@var{a}, @var{b}].
## The size of @var{x} is the common size of the input arguments.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## Default values are @var{a} = 0, @var{b} = 1.
##
## @seealso{unifcdf, unifpdf, unifrnd, unifstat}
## @end deftypefn

function x = unifinv (p, a = 0, b = 1)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  if (! isscalar (p) || ! isscalar (a) || ! isscalar (b))
    [retval, p, a, b] = common_size (p, a, b);
    if (retval > 0)
      error ("unifinv: P, A, and B must be of common size or scalars.");
    endif
  endif

  if (iscomplex (p) || iscomplex (a) || iscomplex (b))
    error ("unifinv: P, A, and B must not be complex.");
  endif

  if (isa (p, "single") || isa (a, "single") || isa (b, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  k = (p >= 0) & (p <= 1) & (a < b);
  x(k) = a(k) + p(k) .* (b(k) - a(k));

endfunction


%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (unifinv (p, ones (1,5), 2*ones (1,5)), [NaN 1 1.5 2 NaN])
%!assert (unifinv (p), [NaN 1 1.5 2 NaN] - 1)
%!assert (unifinv (p, 0), [NaN 1 1.5 2 NaN] - 1)
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
%!error unifinv ()
%!error unifinv (1,2,3,4)
%!error unifinv (ones (3), ones (2), ones (2))
%!error unifinv (ones (2), ones (3), ones (2))
%!error unifinv (ones (2), ones (2), ones (3))
%!error unifinv (i, 2, 2)
%!error unifinv (2, i, 2)
%!error unifinv (2, 2, i)
