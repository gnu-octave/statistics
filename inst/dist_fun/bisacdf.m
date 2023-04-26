## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 2018 John Donoghue
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{p} =} bisacdf (@var{x}, @var{a}, @var{b}, @var{mu})
## @deftypefnx {statistics} {@var{p} =} bisacdf (@var{x}, @var{a}, @var{b}, @var{mu}, @qcode{"upper"})
##
## Birnbaum-Saunders cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the Birnbaum-Saunders distribution with shape parameter
## @var{a}, scale parameter @var{b}, and location parameter @var{mu}.  The size
## of @var{p} is the common size of @var{x}, @var{a}, @var{b}, and @var{mu}.  A
## scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## @code{@var{p} = bisacdf (@var{x}, @var{a}, @var{b}, @var{mu}, "upper")}
## computes the upper tail probability of the Birnbaum-Saunders distribution
## with parameters @var{a}, @var{b}, and @var{mu} at the values in @var{x}.
##
## Further information about the Birnbaum-Saunders distribution can be found at
## @url{https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution}
##
## @seealso{bisainv, bisapdf, bisarnd, bisafit, bisalike, bisastat}
## @end deftypefn

function p = bisacdf (x, a, b, mu, uflag)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("bisacdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 4)
    if (! strcmpi (uflag, "upper"))
      error ("bisacdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, A, B, and MU
  if (! isscalar (x) || ! isscalar (a) || ! isscalar (b) ...
                     || ! isscalar (mu))
    [retval, x, a, b, mu] = common_size (x, a, b, mu);
    if (retval > 0)
      error (strcat (["bisacdf: X, A, B, and MU must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, A, B, and MU being reals
  if (iscomplex (x) || iscomplex (a) || iscomplex (b) ...
                    || iscomplex(mu))
    error ("bisacdf: X, A, B, and MU must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (a, "single") || isa (b, "single") ...
                        || isa (mu, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Force NaNs for out of range parameters.
  k = isnan (x) | ! (-Inf < mu) | ! (mu < Inf) ...
                | ! (b > 0) | ! (b < Inf) | ! (a > 0) | ! (a < Inf);
  p(k) = NaN;

  ## Find valid values in parameters and data
  k = (x > mu) & (x <= Inf) & (-Inf < mu) & (mu < Inf) ...
               & (0 < b) & (b < Inf) & (0 < a) & (a < Inf);

  if (isscalar (mu) && isscalar(b) && isscalar(a))
    x_m = x(k) - mu;
    if (uflag)
      z = (-sqrt (x_m ./ b) + sqrt (b ./ x_m)) ./ a;
    else
      z = (sqrt (x_m ./ b) - sqrt (b ./ x_m)) ./ a;
    endif
    p(k) = 0.5 * erfc (-z ./ sqrt (2));
  else
    x_m = x(k) - mu(k);
    if (uflag)
      z = (-sqrt (x_m ./ b(k)) + sqrt (b(k) ./ x_m)) ./ a(k);
    else
      z = (sqrt (x_m ./ b(k)) - sqrt (b(k) ./ x_m)) ./ a(k);
    endif
    p(k) = 0.5 * erfc (-z ./ sqrt (2));
  endif

endfunction

%!demo
%! ## Plot various CDFs from the Birnbaum-Saunders distribution
%! x = 0.01:0.01:10;
%! p1 = bisacdf (x, 0.5, 1, 0);
%! p2 = bisacdf (x, 1, 1, 0);
%! p3 = bisacdf (x, 2, 1, 0);
%! p4 = bisacdf (x, 5, 1, 0);
%! p5 = bisacdf (x, 10, 1, 0);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c", x, p5, "-m")
%! grid on
%! legend ({"α = 0.5, β = 1, μ = 0", "α = 1,    β = 1, μ = 0", ...
%!          "α = 2,    β = 1, μ = 0", "α = 5,    β = 1, μ = 0", ...
%!          "α = 10,  β = 1, μ = 0"}, "location", "southeast")
%! title ("Birnbaum-Saunders CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

%!demo
%! ## Plot various CDFs from the Birnbaum-Saunders distribution
%! x = 0.01:0.01:10;
%! p1 = bisacdf (x, 0.3, 1, 0);
%! p2 = bisacdf (x, 0.3, 2, 0);
%! p3 = bisacdf (x, 0.5, 1, 0);
%! p4 = bisacdf (x, 0.5, 3, 0);
%! p5 = bisacdf (x, 0.5, 5, 0);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c", x, p5, "-m")
%! grid on
%! legend ({"α = 0.3, β = 1, μ = 0", "α = 0.3, β = 2, μ = 0", ...
%!          "α = 0.5, β = 1, μ = 0", "α = 0.5, β = 3, μ = 0", ...
%!          "α = 0.5, β = 5, μ = 0"}, "location", "southeast")
%! title ("Birnbaum-Saunders CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test results
%!shared x,y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 1/2, 0.76024993890652337, 1];
%!assert (bisacdf (x, ones (1,5), ones (1,5), zeros (1,5)), y, eps)
%!assert (bisacdf (x, 1, 1, zeros (1,5)), y, eps)
%!assert (bisacdf (x, 1, ones (1,5), 0), y, eps)
%!assert (bisacdf (x, ones (1,5), 1, 0), y, eps)
%!assert (bisacdf (x, 1, 1, 0), y, eps)
%!assert (bisacdf (x, 1, 1, [0, 0, NaN, 0, 0]), [y(1:2), NaN, y(4:5)], eps)
%!assert (bisacdf (x, 1, [1, 1, NaN, 1, 1], 0), [y(1:2), NaN, y(4:5)], eps)
%!assert (bisacdf (x, [1, 1, NaN, 1, 1], 1, 0), [y(1:2), NaN, y(4:5)], eps)
%!assert (bisacdf ([x, NaN], 1, 1, 0), [y, NaN], eps)

## Test class of input preserved
%!assert (bisacdf (single ([x, NaN]), 1, 1, 0), single ([y, NaN]), eps('single'))
%!assert (bisacdf ([x, NaN], 1, 1, single (0)), single ([y, NaN]), eps('single'))
%!assert (bisacdf ([x, NaN], 1, single (1), 0), single ([y, NaN]), eps('single'))
%!assert (bisacdf ([x, NaN], single (1), 1, 0), single ([y, NaN]), eps('single'))

## Test input validation
%!error<bisacdf: function called with too few input arguments.> bisacdf ()
%!error<bisacdf: function called with too few input arguments.> bisacdf (1)
%!error<bisacdf: function called with too few input arguments.> bisacdf (1, 2)
%!error<bisacdf: function called with too few input arguments.> ...
%! bisacdf (1, 2, 3)
%!error<bisacdf: function called with too many inputs> ...
%! bisacdf (1, 2, 3, 4, 5, 6)
%!error<bisacdf: invalid argument for upper tail.> bisacdf (1, 2, 3, 4, "tail")
%!error bisacdf (ones (3), ones (2), ones(2), ones(2))
%!error bisacdf (ones (2), ones (3), ones(2), ones(2))
%!error bisacdf (ones (2), ones (2), ones(3), ones(2))
%!error bisacdf (ones (2), ones (2), ones(2), ones(3))
%!error bisacdf (i, 4, 3, 2)
%!error bisacdf (1, i, 3, 2)
%!error bisacdf (1, 4, i, 2)
%!error bisacdf (1, 4, 3, i)

