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
## @deftypefn  {statistics} {@var{x} =} geoinv (@var{p}, @var{ps})
##
## Inverse of the geometric cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the geometric distribution with probability of success parameter @var{ps}.
## The size of @var{x} is the common size of @var{p} and @var{ps}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## The geometric distribution models the number of failures (@var{p}) of a
## Bernoulli trial with probability @var{ps} before the first success.
##
## Further information about the geometric distribution can be found at
## @url{https://en.wikipedia.org/wiki/Geometric_distribution}
##
## @seealso{geocdf, geopdf, geornd, geostat}
## @end deftypefn

function x = geoinv (p, ps)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("geoinv: function called with too few input arguments.");
  endif

  ## Check for common size of P and PS
  if (! isscalar (ps) || ! isscalar (ps))
    [retval, p, ps] = common_size (p, ps);
    if (retval > 0)
      error ("geoinv: P and PS must be of common size or scalars.");
    endif
  endif

  ## Check for P and PS being reals
  if (iscomplex (p) || iscomplex (ps))
    error ("geoinv: P and PS must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (ps, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Handle edge cases
  k = (p == 1) & (ps >= 0) & (ps <= 1);
  x(k) = Inf;

  ## Get valid instances
  k = (p >= 0) & (p < 1) & (ps > 0) & (ps <= 1);

  ## Compute iCDF
  if (isscalar (ps))
    x(k) = max (ceil (log (1 - p(k)) / log (1 - ps)) - 1, 0);
  else
    x(k) = max (ceil (log (1 - p(k)) ./ log (1 - ps(k))) - 1, 0);
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the geometric distribution
%! p = 0.001:0.001:0.999;
%! x1 = geoinv (p, 0.2);
%! x2 = geoinv (p, 0.5);
%! x3 = geoinv (p, 0.7);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r")
%! grid on
%! ylim ([0, 10])
%! legend ({"ps = 0.2", "ps = 0.5", "ps = 0.7"}, "location", "northwest")
%! title ("Geometric iCDF")
%! xlabel ("probability")
%! ylabel ("values in x (number of failures)")

## Test output
%!shared p
%! p = [-1 0 0.75 1 2];
%!assert (geoinv (p, 0.5*ones (1,5)), [NaN 0 1 Inf NaN])
%!assert (geoinv (p, 0.5), [NaN 0 1 Inf NaN])
%!assert (geoinv (p, 0.5*[1 -1 NaN 4 1]), [NaN NaN NaN NaN NaN])
%!assert (geoinv ([p(1:2) NaN p(4:5)], 0.5), [NaN 0 NaN Inf NaN])

## Test class of input preserved
%!assert (geoinv ([p, NaN], 0.5), [NaN 0 1 Inf NaN NaN])
%!assert (geoinv (single ([p, NaN]), 0.5), single ([NaN 0 1 Inf NaN NaN]))
%!assert (geoinv ([p, NaN], single (0.5)), single ([NaN 0 1 Inf NaN NaN]))

## Test input validation
%!error<geoinv: function called with too few input arguments.> geoinv ()
%!error<geoinv: function called with too few input arguments.> geoinv (1)
%!error<geoinv: P and PS must be of common size or scalars.> ...
%! geoinv (ones (3), ones (2))
%!error<geoinv: P and PS must be of common size or scalars.> ...
%! geoinv (ones (2), ones (3))
%!error<geoinv: P and PS must not be complex.> ...
%! geoinv (i, 2)
%!error<geoinv: P and PS must not be complex.> ...
%! geoinv (2, i)
