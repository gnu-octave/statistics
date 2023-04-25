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
## @deftypefn  {statistics} {@var{x} =} betainv (@var{p}, @var{a}, @var{b})
##
## Inverse of the Beta distribution (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Beta distribution with parameters @var{a} and @var{b}.  The
## size of @var{x} is the common size of @var{x}, @var{a} and @var{b}.  A scalar
## input functions as a constant matrix of the same size as the other inputs.
##
## Further information about the Beta distribution can be found at
## @url{https://en.wikipedia.org/wiki/Beta_distribution}
##
## @seealso{betacdf, betapdf, betarnd, betafit, betalike, betastat}
## @end deftypefn

function x = betainv (p, a, b)

  if (nargin != 3)
    print_usage ();
  endif

  if (! isscalar (a) || ! isscalar (b))
    [retval, p, a, b] = common_size (p, a, b);
    if (retval > 0)
      error ("betainv: P, A, and B must be of common size or scalars.");
    endif
  endif

  if (iscomplex (p) || iscomplex (a) || iscomplex (b))
    error ("betainv: P, A, and B must not be complex.");
  endif

  if (isa (p, "single") || isa (a, "single") || isa (b, "single"))
    x = zeros (size (p), "single");
  else
    x = zeros (size (p));
  endif

  k = (p < 0) | (p > 1) | !(a > 0) | !(b > 0) | isnan (p);
  x(k) = NaN;

  k = (p == 1) & (a > 0) & (b > 0);
  x(k) = 1;

  k = find ((p > 0) & (p < 1) & (a > 0) & (b > 0));
  if (! isempty (k))
    if (! isscalar (a) || ! isscalar (b))
      a = a(k);
      b = b(k);
      y = a ./ (a + b);
    else
      y = a / (a + b) * ones (size (k));
    endif
    p = p(k);

    if (isa (y, "single"))
      myeps = eps ("single");
    else
      myeps = eps;
    endif

    l = find (y < myeps);
    if (any (l))
      y(l) = sqrt (myeps) * ones (length (l), 1);
    endif
    l = find (y > 1 - myeps);
    if (any (l))
      y(l) = 1 - sqrt (myeps) * ones (length (l), 1);
    endif

    y_new = y;
    loopcnt = 0;
    do
      y_old = y_new;
      h     = (betacdf (y_old, a, b) - p) ./ betapdf (y_old, a, b);
      y_new = y_old - h;
      ind   = find (y_new <= myeps);
      if (any (ind))
        y_new(ind) = y_old(ind) / 10;
      endif
      ind = find (y_new >= 1 - myeps);
      if (any (ind))
        y_new(ind) = 1 - (1 - y_old(ind)) / 10;
      endif
      h = y_old - y_new;
    until (max (abs (h)) < sqrt (myeps) || ++loopcnt == 40)

    if (loopcnt == 40)
      warning ("betainv: calculation failed to converge for some values.");
    endif

    x(k) = y_new;
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the Beta distribution
%! p = 0.001:0.001:0.999;
%! x1 = betainv (p, 0.5, 0.5);
%! x2 = betainv (p, 5, 1);
%! x3 = betainv (p, 1, 3);
%! x4 = betainv (p, 2, 2);
%! x5 = betainv (p, 2, 5);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c", p, x5, "-m")
%! grid on
%! legend ({"α = β = 0.5", "α = 5, β = 1", "α = 1, β = 3", ...
%!          "α = 2, β = 2", "α = 2, β = 5"}, "location", "southeast")
%! title ("Beta iCDF")
%! xlabel ("probability")
%! ylabel ("x")


%!shared p
%! p = [-1 0 0.75 1 2];
%!assert (betainv (p, ones (1,5), 2*ones (1,5)), [NaN 0 0.5 1 NaN], eps)
%!assert (betainv (p, 1, 2*ones (1,5)), [NaN 0 0.5 1 NaN], eps)
%!assert (betainv (p, ones (1,5), 2), [NaN 0 0.5 1 NaN], eps)
%!assert (betainv (p, [1 0 NaN 1 1], 2), [NaN NaN NaN 1 NaN])
%!assert (betainv (p, 1, 2*[1 0 NaN 1 1]), [NaN NaN NaN 1 NaN])
%!assert (betainv ([p(1:2) NaN p(4:5)], 1, 2), [NaN 0 NaN 1 NaN])

## Test class of input preserved
%!assert (betainv ([p, NaN], 1, 2), [NaN 0 0.5 1 NaN NaN], eps)
%!assert (betainv (single ([p, NaN]), 1, 2), single ([NaN 0 0.5 1 NaN NaN]))
%!assert (betainv ([p, NaN], single (1), 2), single ([NaN 0 0.5 1 NaN NaN]), eps("single"))
%!assert (betainv ([p, NaN], 1, single (2)), single ([NaN 0 0.5 1 NaN NaN]), eps("single"))

## Test input validation
%!error betainv ()
%!error betainv (1)
%!error betainv (1,2)
%!error betainv (1,2,3,4)
%!error betainv (ones (3), ones (2), ones (2))
%!error betainv (ones (2), ones (3), ones (2))
%!error betainv (ones (2), ones (2), ones (3))
%!error betainv (i, 2, 2)
%!error betainv (2, i, 2)
%!error betainv (2, 2, i)
