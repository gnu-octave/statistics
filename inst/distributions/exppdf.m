## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} @var{y} = exppdf (@var{x}, @var{mu})
##
## Exponential probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the exponential distribution with mean @var{mu}.  The size of
## @var{y} is the common size of @var{x} and @var{mu}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## @seealso{expcdf, expinv, exprnd, expfit, explike, expstat}
## @end deftypefn

function y = exppdf (x, mu)

  if (nargin != 2)
    print_usage ();
  endif

  if (! isscalar (mu) || ! isscalar (mu))
    [retval, x, mu] = common_size (x, mu);
    if (retval > 0)
      error ("exppdf: X and MU must be of common size or scalars.");
    endif
  endif

  if (iscomplex (x) || iscomplex (mu))
    error ("exppdf: X and MU must not be complex.");
  endif

  if (isa (x, "single") || isa (mu, "single"))
    y = zeros (size (x), "single");
  else
    y = zeros (size (x));
  endif

  k = isnan (x) | !(mu > 0);
  y(k) = NaN;

  k = (x >= 0) & (x < Inf) & (mu > 0);
  if (isscalar (mu))
    y(k) = exp (-x(k) / mu) / mu;
  else
    y(k) = exp (-x(k) ./ mu(k)) ./ mu(k);
  endif

endfunction


%!shared x,y
%! x = [-1 0 0.5 1 Inf];
%! y = gampdf (x, 1, 2);
%!assert (exppdf (x, 2*ones (1,5)), y)
%!assert (exppdf (x, 2*[1 0 NaN 1 1]), [y(1) NaN NaN y(4:5)])
%!assert (exppdf ([x, NaN], 2), [y, NaN])

## Test class of input preserved
%!assert (exppdf (single ([x, NaN]), 2), single ([y, NaN]))
%!assert (exppdf ([x, NaN], single (2)), single ([y, NaN]))

## Test input validation
%!error exppdf ()
%!error exppdf (1)
%!error exppdf (1,2,3)
%!error exppdf (ones (3), ones (2))
%!error exppdf (ones (2), ones (3))
%!error exppdf (i, 2)
%!error exppdf (2, i)
