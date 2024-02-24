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
## @deftypefn  {statistics} {@var{x} =} logiinv (@var{p}, @var{mu}, @var{sigma})
##
## Inverse of the logistic cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the logistic distribution with location parameter @var{mu} and scale
## parameter @var{sigma}.  The size of @var{p} is the common size of @var{x},
## @var{mu}, and @var{sigma}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Both parameters must be reals and @qcode{@var{sigma} > 0}.
## For @qcode{@var{sigma} <= 0}, @qcode{NaN} is returned.
##
## Further information about the logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Logistic_distribution}
##
## @seealso{logicdf, logipdf, logirnd, logifit, logilike, logistat}
## @end deftypefn

function x = logiinv (p, mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("logiinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, MU, and SIGMA
  if (! isscalar (p) || ! isscalar (mu) || ! isscalar(sigma))
    [retval, p, mu, sigma] = common_size (p, mu, sigma);
    if (retval > 0)
      error ("logiinv: P, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (sigma))
    error ("logiinv: P, MU, and SIGMA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (mu, "single") || isa (sigma, "single"));
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  k = (p == 0) & (sigma > 0);
  x(k) = -Inf;

  k = (p == 1) & (sigma > 0);
  x(k) = Inf;

  k = (p > 0) & (p < 1) & (sigma > 0);
  x(k) = mu(k) + sigma(k) .* log (p(k) ./ (1 - p(k)));

endfunction

%!demo
%! ## Plot various iCDFs from the logistic distribution
%! p = 0.001:0.001:0.999;
%! x1 = logiinv (p, 5, 2);
%! x2 = logiinv (p, 9, 3);
%! x3 = logiinv (p, 9, 4);
%! x4 = logiinv (p, 6, 2);
%! x5 = logiinv (p, 2, 1);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c", p, x5, "-m")
%! grid on
%! legend ({"μ = 5, σ = 2", "μ = 9, σ = 3", "μ = 9, σ = 4", ...
%!          "μ = 6, σ = 2", "μ = 2, σ = 1"}, "location", "southeast")
%! title ("Logistic iCDF")
%! xlabel ("probability")
%! ylabel ("x")

## Test output
%!test
%! p = [0.01:0.01:0.99];
%! assert (logiinv (p, 0, 1), log (p ./ (1-p)), 25*eps);
%!shared p
%! p = [-1 0 0.5 1 2];
%!assert (logiinv (p, 0, 1), [NaN -Inf 0 Inf NaN])
%!assert (logiinv (p, 0, [-1, 0, 1, 2, 3]), [NaN NaN 0 Inf NaN])

## Test class of input preserved
%!assert (logiinv ([p, NaN], 0, 1), [NaN -Inf 0 Inf NaN NaN])
%!assert (logiinv (single ([p, NaN]), 0, 1), single ([NaN -Inf 0 Inf NaN NaN]))
%!assert (logiinv ([p, NaN], single (0), 1), single ([NaN -Inf 0 Inf NaN NaN]))
%!assert (logiinv ([p, NaN], 0, single (1)), single ([NaN -Inf 0 Inf NaN NaN]))

## Test input validation
%!error<logiinv: function called with too few input arguments.> logiinv ()
%!error<logiinv: function called with too few input arguments.> logiinv (1)
%!error<logiinv: function called with too few input arguments.> ...
%! logiinv (1, 2)
%!error<logiinv: P, MU, and SIGMA must be of common size or scalars.> ...
%! logiinv (1, ones (2), ones (3))
%!error<logiinv: P, MU, and SIGMA must be of common size or scalars.> ...
%! logiinv (ones (2), 1, ones (3))
%!error<logiinv: P, MU, and SIGMA must be of common size or scalars.> ...
%! logiinv (ones (2), ones (3), 1)
%!error<logiinv: P, MU, and SIGMA must not be complex.> logiinv (i, 2, 3)
%!error<logiinv: P, MU, and SIGMA must not be complex.> logiinv (1, i, 3)
%!error<logiinv: P, MU, and SIGMA must not be complex.> logiinv (1, 2, i)
