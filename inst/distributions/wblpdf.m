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
## @deftypefn  {statistics} @var{y} = wblinv (@var{x})
## @deftypefnx {statistics} @var{y} = wblinv (@var{x}, @var{lambda})
## @deftypefnx {statistics} @var{y} = wblinv (@var{x}, @var{lambda}, @var{xk})
##
## Weibull probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## at @var{x} of the Weibull distribution with parameters @var{lambda} and
## @var{xk}.  The size of @var{y} is the common size of @var{x}, @var{lambda},
## and @var{xk}.  A scalar input functions as a constant matrix of the same size
## as the other inputs.
##
## Default values are @var{lambda} = 1, @var{xk} = 1.
##
## @seealso{wblcdf, wblinv, wblrnd, wblstat, wblplot}
## @end deftypefn

function y = wblpdf (x, lambda = 1, k = 1)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  if (! isscalar (lambda) || ! isscalar (k))
    [retval, x, lambda, k] = common_size (x, lambda, k);
    if (retval > 0)
      error ("wblpdf: X, LAMBDA, and K must be of common size or scalars.");
    endif
  endif

  if (iscomplex (x) || iscomplex (lambda) || iscomplex (k))
    error ("wblpdf: X, LAMBDA, and K must not be complex.");
  endif

  if (isa (x, "single") || isa (lambda, "single") || isa (k, "single"))
    y = NaN (size (x), "single");
  else
    y = NaN (size (x));
  endif

  ok = ((lambda > 0) & (lambda < Inf) & (k > 0) & (k < Inf));

  xk = (x < 0) & ok;
  y(xk) = 0;

  xk = (x >= 0) & (x < Inf) & ok;
  if (isscalar (lambda) && isscalar (k))
    y(xk) = (k * (lambda .^ -k) ...
              .* (x(xk) .^ (k - 1)) ...
              .* exp (- (x(xk) / lambda) .^ k));
  else
    y(xk) = (k(xk) .* (lambda(xk) .^ -k(xk)) ...
              .* (x(xk) .^ (k(xk) - 1)) ...
              .* exp (- (x(xk) ./ lambda(xk)) .^ k(xk)));
  endif

endfunction


%!shared x,y
%! x = [-1 0 0.5 1 Inf];
%! y = [0, exp(-x(2:4)), NaN];
%!assert (wblpdf (x, ones (1,5), ones (1,5)), y)
%!assert (wblpdf (x, 1, ones (1,5)), y)
%!assert (wblpdf (x, ones (1,5), 1), y)
%!assert (wblpdf (x, [0 NaN Inf 1 1], 1), [NaN NaN NaN y(4:5)])
%!assert (wblpdf (x, 1, [0 NaN Inf 1 1]), [NaN NaN NaN y(4:5)])
%!assert (wblpdf ([x, NaN], 1, 1), [y, NaN])

## Test class of input preserved
%!assert (wblpdf (single ([x, NaN]), 1, 1), single ([y, NaN]))
%!assert (wblpdf ([x, NaN], single (1), 1), single ([y, NaN]))
%!assert (wblpdf ([x, NaN], 1, single (1)), single ([y, NaN]))

## Test input validation
%!error wblpdf ()
%!error wblpdf (1,2,3,4)
%!error wblpdf (ones (3), ones (2), ones (2))
%!error wblpdf (ones (2), ones (3), ones (2))
%!error wblpdf (ones (2), ones (2), ones (3))
%!error wblpdf (i, 2, 2)
%!error wblpdf (2, i, 2)
%!error wblpdf (2, 2, i)
