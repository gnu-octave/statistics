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
## @deftypefn  {statistics} {@var{x} =} laplaceinv (@var{p}, @var{mu}, @var{beta})
##
## Inverse of the Laplace cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the Laplace distribution with location parameter @var{mu} and
## scale parameter (i.e. "diversity") @var{beta}.  The size of @var{x} is the
## common size of @var{p}, @var{mu}, and @var{beta}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## Both parameters must be reals and @qcode{@var{beta} > 0}.
## For @qcode{@var{beta} <= 0}, @qcode{NaN} is returned.
##
## Further information about the Laplace distribution can be found at
## @url{https://en.wikipedia.org/wiki/Laplace_distribution}
##
## @seealso{laplaceinv, laplacepdf, laplacernd}
## @end deftypefn

function x = laplaceinv (p, mu, beta)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("laplaceinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, MU, and BETA
  if (! isscalar (p) || ! isscalar (mu) || ! isscalar(beta))
    [retval, p, mu, beta] = common_size (p, mu, beta);
    if (retval > 0)
      error (strcat (["laplaceinv: P, MU, and BETA must be of"], ...
                     [" common size or scalars."]));
    endif
  endif

  ## Check for X, MU, and BETA being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (beta))
    error ("laplaceinv: P, MU, and BETA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (mu, "single") || isa (beta, "single"));
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Compute Laplace iCDF
  k = (p >= 0) & (p <= 1) & (beta > 0);
  x(k) = mu(k) + beta(k) .* ((p(k) < 1/2) .* log (2 .* p(k)) - ...
                             (p(k) > 1/2) .* log (2 .* (1 - p(k))));

endfunction

%!demo
%! ## Plot various iCDFs from the Laplace distribution
%! p = 0.001:0.001:0.999;
%! x1 = cauchyinv (p, 0, 1);
%! x2 = cauchyinv (p, 0, 2);
%! x3 = cauchyinv (p, 0, 4);
%! x4 = cauchyinv (p, -5, 4);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c")
%! grid on
%! ylim ([-10, 10])
%! legend ({"μ = 0, β = 1", "μ = 0, β = 2", ...
%!          "μ = 0, β = 4", "μ = -5, β = 4"}, "location", "northwest")
%! title ("Laplace iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test results
%!shared p, x
%! p = [-1 0 0.5 1 2];
%! x = [NaN, -Inf, 0, Inf, NaN];
%!assert (laplaceinv (p, 0, 1), x)
%!assert (laplaceinv (p, 0, [-2, -1, 0, 1, 2]), [nan(1, 3), Inf, NaN])
%!assert (laplaceinv ([p, NaN], 0, 1), [x, NaN])

## Test class of input preserved
%!assert (laplaceinv (single ([p, NaN]), 0, 1), single ([x, NaN]))
%!assert (laplaceinv ([p, NaN], single (0), 1), single ([x, NaN]))
%!assert (laplaceinv ([p, NaN], 0, single (1)), single ([x, NaN]))

## Test input validation
%!error<laplaceinv: function called with too few input arguments.> laplaceinv ()
%!error<laplaceinv: function called with too few input arguments.> laplaceinv (1)
%!error<laplaceinv: function called with too few input arguments.> ...
%! laplaceinv (1, 2)
%!error<laplaceinv: function called with too many inputs> laplaceinv (1, 2, 3, 4)
%!error<laplaceinv: P, MU, and BETA must be of common size or scalars.> ...
%! laplaceinv (1, ones (2), ones (3))
%!error<laplaceinv: P, MU, and BETA must be of common size or scalars.> ...
%! laplaceinv (ones (2), 1, ones (3))
%!error<laplaceinv: P, MU, and BETA must be of common size or scalars.> ...
%! laplaceinv (ones (2), ones (3), 1)
%!error<laplaceinv: P, MU, and BETA must not be complex.> laplaceinv (i, 2, 3)
%!error<laplaceinv: P, MU, and BETA must not be complex.> laplaceinv (1, i, 3)
%!error<laplaceinv: P, MU, and BETA must not be complex.> laplaceinv (1, 2, i)
