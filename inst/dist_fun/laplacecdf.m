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
## @deftypefn  {statistics} {@var{p} =} laplacecdf (@var{x}, @var{mu}, @var{beta})
## @deftypefnx {statistics} {@var{p} =} laplacecdf (@var{x}, @var{mu}, @var{beta}, @qcode{"upper"})
##
## Laplace cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the Laplace distribution with location parameter @var{mu} and scale
## parameter (i.e. "diversity") @var{beta}.  The size of @var{p} is the common
## size of @var{x}, @var{mu}, and @var{beta}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## Both parameters must be reals and @qcode{@var{beta} > 0}.
## For @qcode{@var{beta} <= 0}, @qcode{NaN} is returned.
##
## @code{@var{p} = laplacecdf (@var{x}, @var{mu}, @var{beta}, "upper")} computes
## the upper tail probability of the Laplace distribution with parameters
## @var{mu} and @var{beta}, at the values in @var{x}.
##
## Further information about the Laplace distribution can be found at
## @url{https://en.wikipedia.org/wiki/Laplace_distribution}
##
## @seealso{laplaceinv, laplacepdf, laplacernd}
## @end deftypefn

function p = laplacecdf (x, mu, beta, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("laplacecdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("laplacecdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, MU, and BETA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(beta))
    [retval, x, mu, beta] = common_size (x, mu, beta);
    if (retval > 0)
      error (strcat ("laplacecdf: X, MU, and BETA must be of", ...
                     " common size or scalars."));
    endif
  endif

  ## Check for X, MU, and BETA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (beta))
    error ("laplacecdf: X, MU, and BETA must not be complex.");
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (beta, "single"));
    p = NaN (size (x), "single");
  else
    p = NaN (size (x));
  endif

  ## Find normal and edge cases
  k1 = (x == -Inf) & (beta > 0);
  k2 = (x == Inf) & (beta > 0);
  k = ! k1 & ! k2 & (beta > 0);

  ## Compute Laplace CDF
  if (uflag)
    p(k1) = 1;
    p(k2) = 0;
    p(k) = (1 + sign (-x(k) + mu(k)) .* ...
          (1 - exp (- abs (-x(k) + mu(k)) ./ beta(k)))) ./ 2;
  else
    p(k1) = 0;
    p(k2) = 1;
    p(k) = (1 + sign (x(k) - mu(k)) .* ...
          (1 - exp (- abs (x(k) - mu(k)) ./ beta(k)))) ./ 2;
  endif

endfunction

%!demo
%! ## Plot various CDFs from the Laplace distribution
%! x = -10:0.01:10;
%! p1 = laplacecdf (x, 0, 1);
%! p2 = laplacecdf (x, 0, 2);
%! p3 = laplacecdf (x, 0, 4);
%! p4 = laplacecdf (x, -5, 4);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c")
%! grid on
%! xlim ([-10, 10])
%! legend ({"μ = 0, β = 1", "μ = 0, β = 2", ...
%!          "μ = 0, β = 4", "μ = -5, β = 4"}, "location", "southeast")
%! title ("Laplace CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-Inf, -log(2), 0, log(2), Inf];
%! y = [0, 1/4, 1/2, 3/4, 1];
%!assert (laplacecdf ([x, NaN], 0, 1), [y, NaN])
%!assert (laplacecdf (x, 0, [-2, -1, 0, 1, 2]), [nan(1, 3), 0.75, 1])

## Test class of input preserved
%!assert (laplacecdf (single ([x, NaN]), 0, 1), single ([y, NaN]), eps ("single"))
%!assert (laplacecdf ([x, NaN], single (0), 1), single ([y, NaN]), eps ("single"))
%!assert (laplacecdf ([x, NaN], 0, single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<laplacecdf: function called with too few input arguments.> laplacecdf ()
%!error<laplacecdf: function called with too few input arguments.> laplacecdf (1)
%!error<laplacecdf: function called with too few input arguments.> ...
%! laplacecdf (1, 2)
%!error<laplacecdf: function called with too many inputs> ...
%! laplacecdf (1, 2, 3, 4, 5)
%!error<laplacecdf: invalid argument for upper tail.> laplacecdf (1, 2, 3, "tail")
%!error<laplacecdf: invalid argument for upper tail.> laplacecdf (1, 2, 3, 4)
%!error<laplacecdf: X, MU, and BETA must be of common size or scalars.> ...
%! laplacecdf (ones (3), ones (2), ones (2))
%!error<laplacecdf: X, MU, and BETA must be of common size or scalars.> ...
%! laplacecdf (ones (2), ones (3), ones (2))
%!error<laplacecdf: X, MU, and BETA must be of common size or scalars.> ...
%! laplacecdf (ones (2), ones (2), ones (3))
%!error<laplacecdf: X, MU, and BETA must not be complex.> laplacecdf (i, 2, 2)
%!error<laplacecdf: X, MU, and BETA must not be complex.> laplacecdf (2, i, 2)
%!error<laplacecdf: X, MU, and BETA must not be complex.> laplacecdf (2, 2, i)
