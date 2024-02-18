## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{y} =} wblinv (@var{x})
## @deftypefnx {statistics} {@var{y} =} wblinv (@var{x}, @var{lambda})
## @deftypefnx {statistics} {@var{y} =} wblinv (@var{x}, @var{lambda}, @var{k})
##
## Weibull probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the Weibull distribution with scale parameter @var{lambda} and shpe
## parameter @var{k}.  The size of @var{y} is the common size of @var{x},
## @var{lambda}, and @var{k}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Default values are @var{lambda} = 1, @var{k} = 1.
##
## Further information about the Weibull distribution can be found at
## @url{https://en.wikipedia.org/wiki/Weibull_distribution}
##
## @seealso{wblcdf, wblinv, wblrnd, wblfit, wbllike, wblstat, wblplot}
## @end deftypefn

function y = wblpdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    error ("wblpdf: invalid number of input arguments.");
  endif

  ## Get extra arguments (if they exist) or add defaults
  if (numel (varargin) > 0)
    lambda = varargin{1};
  else
    lambda = 1;
  endif
  if (numel (varargin) > 1)
    k = varargin{2};
  else
    k = 1;
  endif

  ## Check for common size of X, LAMBDA, and K
  if (! isscalar (lambda) || ! isscalar (k))
    [retval, x, lambda, k] = common_size (x, lambda, k);
    if (retval > 0)
      error ("wblpdf: X, LAMBDA, and K must be of common size or scalars.");
    endif
  endif

  ## Check for X, LAMBDA, and K being reals
  if (iscomplex (x) || iscomplex (lambda) || iscomplex (k))
    error ("wblpdf: X, LAMBDA, and K must not be complex.");
  endif

  ## Check for class type
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

%!demo
%! ## Plot various PDFs from the Weibul distribution
%! x = 0:0.001:2.5;
%! y1 = wblpdf (x, 1, 0.5);
%! y2 = wblpdf (x, 1, 1);
%! y3 = wblpdf (x, 1, 1.5);
%! y4 = wblpdf (x, 1, 5);
%! plot (x, y1, "-b", x, y2, "-r", x, y3, "-m", x, y4, "-g")
%! grid on
%! ylim ([0, 2.5])
%! legend ({"λ = 5, k = 0.5", "λ = 9, k = 1", ...
%!          "λ = 6, k = 1.5", "λ = 2, k = 5"}, "location", "northeast")
%! title ("Weibul PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
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
