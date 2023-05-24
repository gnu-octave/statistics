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
## @deftypefn  {statistics} {@var{y} =} exppdf (@var{x})
## @deftypefnx {statistics} {@var{y} =} exppdf (@var{x}, @var{mu})
##
## Exponential probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
##  of the exponential distribution with mean parameter @var{mu}.  The size of
## @var{y} is the common size of @var{x} and @var{mu}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## Default value for @var{mu} = 1.
##
## A common alternative parameterization of the exponential distribution is to
## use the parameter @math{λ} defined as the mean number of events in an
## interval as opposed to the parameter @math{μ}, which is the mean wait time
## for an event to occur. @math{λ} and @math{μ} are reciprocals,
## i.e. @math{μ = 1 / λ}.
##
## Further information about the exponential distribution can be found at
## @url{https://en.wikipedia.org/wiki/Exponential_distribution}
##
## @seealso{expcdf, expinv, exprnd, expfit, explike, expstat}
## @end deftypefn

function y = exppdf (x, mu)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("exppdf: function called with too few input arguments.");
  endif

  ## Add defaults (if missing input arguments)
  if (nargin < 2)
    mu = 0;
  endif

  ## Check for common size of X and MU
  if (! isscalar (x) || ! isscalar (mu))
    [retval, x, mu] = common_size (x, mu);
    if (retval > 0)
      error ("exppdf: X and MU must be of common size or scalars.");
    endif
  endif

  ## Check for X and MU being reals
  if (iscomplex (x) || iscomplex (mu))
    error ("exppdf: X and MU must not be complex.");
  endif

  ## Check for appropriate class
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

%!demo
%! ## Plot various PDFs from the exponential distribution
%! x = 0:0.01:5;
%! y1 = exppdf (x, 2/3);
%! y2 = exppdf (x, 1.0);
%! y3 = exppdf (x, 2.0);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r")
%! grid on
%! ylim ([0, 1.5])
%! legend ({"μ = 2/3", "μ = 1", "μ = 2"}, "location", "northeast")
%! title ("Exponential PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
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
%!error<exppdf: function called with too few input arguments.> exppdf ()
%!error<exppdf: function called with too many inputs> exppdf (1,2,3)
%!error<exppdf: X and MU must be of common size or scalars.> ...
%! exppdf (ones (3), ones (2))
%!error<exppdf: X and MU must be of common size or scalars.> ...
%! exppdf (ones (2), ones (3))
%!error<exppdf: X and MU must not be complex.> exppdf (i, 2)
%!error<exppdf: X and MU must not be complex.> exppdf (2, i)
