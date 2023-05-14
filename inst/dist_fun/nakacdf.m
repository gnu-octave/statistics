## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1995-2015 Kurt Hornik
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
## @deftypefn  {statistics} {@var{p} =} nakacdf (@var{x}, @var{mu}, @var{omega})
## @deftypefnx {statistics} {@var{p} =} nakacdf (@var{x}, @var{mu}, @var{omega}, @qcode{"upper"})
##
## Nakagami cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the Nakagami distribution with shape parameter @var{mu} and spread
## parameter @var{omega}.  The size of @var{p} is the common size of @var{x},
## @var{mu}, and @var{omega}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Both parameters must be positive reals and @qcode{@var{mu} >= 0.5}.  For
## @qcode{@var{mu} < 0.5} or @qcode{@var{omega} <= 0}, @qcode{NaN} is returned.
##
## @code{@var{p} = nakacdf (@var{x}, @var{mu}, @var{omega}, "upper")} computes
## the upper tail probability of the Nakagami distribution with parameters
## @var{mu} and @var{beta}, at the values in @var{x}.
##
## Further information about the Nakagami distribution can be found at
## @url{https://en.wikipedia.org/wiki/Nakagami_distribution}
##
## @seealso{nakainv, nakapdf, nakarnd, nakafitm, nakalike}
## @end deftypefn

function p = nakacdf (x, mu, omega, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("nakacdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("nakacdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, MU, and OMEGA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (omega))
    [retval, x, mu, omega] = common_size (x, mu, omega);
    if (retval > 0)
      error ("nakacdf: X, MU, and OMEGA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and OMEGA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (omega))
    error ("nakacdf: X, MU, and OMEGA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (omega, "single"))
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Force invalid parameters and missing data to NaN
  k1 = isnan (x) | ! (mu >= 0.5) | ! (omega > 0);
  p(k1) = NaN;

  ## Find normal and edge cases
  k2 = (x == Inf) & (mu >= 0.5) & (mu < Inf) & (omega > 0) & (omega < Inf);
  k = (x > 0) & (x < Inf) & (mu >= 0.5) & (mu < Inf) ...
                          & (omega > 0) & (omega < Inf);

  ## Compute Nakagami CDF
  if (uflag)
    p(k2) = 0;
    left = mu .* ones (size (x));
    right = (mu ./ omega) .* x .^ 2;
    p(k) = gammainc (right(k), left(k), "upper");
  else
    p(k2) = 1;
    left = mu .* ones (size (x));
    right = (mu ./ omega) .* x .^ 2;
    p(k) = gammainc (right(k), left(k));
  endif

endfunction

%!demo
%! ## Plot various CDFs from the Nakagami distribution
%! x = 0:0.01:3;
%! p1 = nakacdf (x, 0.5, 1);
%! p2 = nakacdf (x, 1, 1);
%! p3 = nakacdf (x, 1, 2);
%! p4 = nakacdf (x, 1, 3);
%! p5 = nakacdf (x, 2, 1);
%! p6 = nakacdf (x, 2, 2);
%! p7 = nakacdf (x, 5, 1);
%! plot (x, p1, "-r", x, p2, "-g", x, p3, "-y", x, p4, "-m", ...
%!       x, p5, "-k", x, p6, "-b", x, p7, "-c")
%! grid on
%! xlim ([0, 3])
%! legend ({"μ = 0.5, ω = 1", "μ = 1, ω = 1", "μ = 1, ω = 2", ...
%!          "μ = 1, ω = 3", "μ = 2, ω = 1", "μ = 2, ω = 2", ...
%!          "μ = 5, ω = 1"}, "location", "southeast")
%! title ("Nakagami CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-1, 0, 1, 2, Inf];
%! y = [0, 0, 0.63212055882855778, 0.98168436111126578, 1];
%!assert (nakacdf (x, ones (1,5), ones (1,5)), y, eps)
%!assert (nakacdf (x, 1, 1), y, eps)
%!assert (nakacdf (x, [1, 1, NaN, 1, 1], 1), [y(1:2), NaN, y(4:5)])
%!assert (nakacdf (x, 1, [1, 1, NaN, 1, 1]), [y(1:2), NaN, y(4:5)])
%!assert (nakacdf ([x, NaN], 1, 1), [y, NaN], eps)

## Test class of input preserved
%!assert (nakacdf (single ([x, NaN]), 1, 1), single ([y, NaN]), eps("single"))
%!assert (nakacdf ([x, NaN], single (1), 1), single ([y, NaN]), eps("single"))
%!assert (nakacdf ([x, NaN], 1, single (1)), single ([y, NaN]), eps("single"))

## Test input validation
%!error<nakacdf: function called with too few input arguments.> nakacdf ()
%!error<nakacdf: function called with too few input arguments.> nakacdf (1)
%!error<nakacdf: function called with too few input arguments.> nakacdf (1, 2)
%!error<nakacdf: invalid argument for upper tail.> nakacdf (1, 2, 3, "tail")
%!error<nakacdf: invalid argument for upper tail.> nakacdf (1, 2, 3, 4)
%!error<nakacdf: X, MU, and OMEGA must be of common size or scalars.> ...
%! nakacdf (ones (3), ones (2), ones (2))
%!error<nakacdf: X, MU, and OMEGA must be of common size or scalars.> ...
%! nakacdf (ones (2), ones (3), ones (2))
%!error<nakacdf: X, MU, and OMEGA must be of common size or scalars.> ...
%! nakacdf (ones (2), ones (2), ones (3))
%!error<nakacdf: X, MU, and OMEGA must not be complex.> nakacdf (i, 2, 2)
%!error<nakacdf: X, MU, and OMEGA must not be complex.> nakacdf (2, i, 2)
%!error<nakacdf: X, MU, and OMEGA must not be complex.> nakacdf (2, 2, i)
