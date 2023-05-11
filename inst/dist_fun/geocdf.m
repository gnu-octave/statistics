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
## @deftypefn  {statistics} {@var{p} =} geocdf (@var{x}, @var{ps})
## @deftypefnx {statistics} {@var{p} =} geocdf (@var{x}, @var{ps}, @qcode{"upper"})
##
## Geometric cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the geometric distribution with probability of success parameter
## @var{ps}.  The size of @var{p} is the common size of @var{x} and @var{ps}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## @code{@var{p} = geocdf (@var{x}, @var{ps}, "upper")} computes the upper tail
## probability of the geometric distribution with parameter @var{ps}, at the
## values in @var{x}.
##
## The geometric distribution models the number of failures (@var{x}) of a
## Bernoulli trial with probability @var{ps} before the first success.
##
## Further information about the geometric distribution can be found at
## @url{https://en.wikipedia.org/wiki/Geometric_distribution}
##
## @seealso{geoinv, geopdf, geornd, geostat}
## @end deftypefn

function p = geocdf (x, ps, uflag)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("geocdf: function called with too few input arguments.");
  endif

  ## Check for common size of X and PS
  if (! isscalar (x) || ! isscalar (ps))
    [retval, x, ps] = common_size (x, ps);
    if (retval > 0)
      error ("geocdf: X and PS must be of common size or scalars.");
    endif
  endif

  ## Check for X and PS being reals
  if (iscomplex (x) || iscomplex (ps))
    error ("geocdf: X and PS must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (ps, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Return NaN for out of range parameters
  k = isnan (x) | ! (ps >= 0) | ! (ps <= 1);
  p(k) = NaN;

  ## Return 1 for valid range parameters when X = Inf
  k = (x == Inf) & (ps >= 0) & (ps <= 1);
  p(k) = 1;

  ## Return 0 for X < 0
  x(x < 0) = -1;

  ## Check for "upper" flag
  if (nargin > 2 && strcmpi (uflag, "upper"))
    uflag = true;
  elseif (nargin > 2  && ! strcmpi (uflag, "upper"))
    error ("geocdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Get valid instances
  k = (x >= 0) & (x < Inf) & (x == fix (x)) & (ps > 0) & (ps <= 1);

  ## Compute CDF
  if (uflag)
    if (any (k))
      p(k) = betainc (ps(k), 1, (x(k)) + 1, "upper");
    endif
  else
    if (isscalar (ps))
      p(k) = 1 - ((1 - ps) .^ (x(k) + 1));
    else
      p(k) = 1 - ((1 - ps(k)) .^ (x(k) + 1));
    endif
  endif

endfunction

%!demo
%! ## Plot various CDFs from the geometric distribution
%! x = 0:10;
%! p1 = geocdf (x, 0.2);
%! p2 = geocdf (x, 0.5);
%! p3 = geocdf (x, 0.7);
%! plot (x, p1, "*b", x, p2, "*g", x, p3, "*r")
%! grid on
%! xlim ([0, 10])
%! legend ({"ps = 0.2", "ps = 0.5", "ps = 0.7"}, "location", "southeast")
%! title ("Geometric CDF")
%! xlabel ("values in x (number of failures)")
%! ylabel ("probability")

## Test output
%!test
%! p = geocdf ([1, 2, 3, 4], 0.25);
%! assert (p(1), 0.4375000000, 1e-14);
%! assert (p(2), 0.5781250000, 1e-14);
%! assert (p(3), 0.6835937500, 1e-14);
%! assert (p(4), 0.7626953125, 1e-14);
%!test
%! p = geocdf ([1, 2, 3, 4], 0.25, "upper");
%! assert (p(1), 0.5625000000, 1e-14);
%! assert (p(2), 0.4218750000, 1e-14);
%! assert (p(3), 0.3164062500, 1e-14);
%! assert (p(4), 0.2373046875, 1e-14);
%!shared x, p
%! x = [-1 0 1 Inf];
%! p = [0 0.5 0.75 1];
%!assert (geocdf (x, 0.5*ones (1,4)), p)
%!assert (geocdf (x, 0.5), p)
%!assert (geocdf (x, 0.5*[-1 NaN 4 1]), [NaN NaN NaN p(4)])
%!assert (geocdf ([x(1:2) NaN x(4)], 0.5), [p(1:2) NaN p(4)])

## Test class of input preserved
%!assert (geocdf ([x, NaN], 0.5), [p, NaN])
%!assert (geocdf (single ([x, NaN]), 0.5), single ([p, NaN]))
%!assert (geocdf ([x, NaN], single (0.5)), single ([p, NaN]))

## Test input validation
%!error<geocdf: function called with too few input arguments.> geocdf ()
%!error<geocdf: function called with too few input arguments.> geocdf (1)
%!error<geocdf: X and PS must be of common size or scalars.> ...
%! geocdf (ones (3), ones (2))
%!error<geocdf: X and PS must be of common size or scalars.> ...
%! geocdf (ones (2), ones (3))
%!error<geocdf: X and PS must not be complex.> geocdf (i, 2)
%!error<geocdf: X and PS must not be complex.> geocdf (2, i)
%!error<geocdf: invalid argument for upper tail.> geocdf (2, 3, "tail")
%!error<geocdf: invalid argument for upper tail.> geocdf (2, 3, 5)
