## Copyright (C) 2023-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{x} =} loglinv (@var{p}, @var{mu}, @var{sigma})
##
## Inverse of the log-logistic cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the log-logistic distribution with mean parameter @var{mu} and scale
## parameter @var{sigma}.  The size of @var{x} is the common size of @var{p},
## @var{mu}, and @var{sigma}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Both parameters, @var{mu} and @var{sigma}, must be positive reals and @var{p}
## is supported in the range @math{[0,1]}, otherwise @qcode{NaN} is returned.
##
## Further information about the loglogistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-logistic_distribution}
##
## OCTAVE/MATLAB use an alternative parameterization given by the pair
## @math{μ, σ}, i.e. @var{mu} and @var{sigma}, in analogy with the logistic
## distribution.  Their relation to the @math{α} and @math{b} parameters used
## in Wikipedia are given below:
##
## @itemize
## @item @qcode{@var{mu} = log (@var{a})}
## @item @qcode{@var{sigma} = 1 / @var{a}}
## @end itemize
##
## @seealso{loglcdf, loglpdf, loglrnd, loglfit, logllike, loglstat}
## @end deftypefn

function x = loglinv (p, mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("loglinv: function called with too few input arguments.");
  endif

  ## Check for common size of P, MU, and SIGMA
  if (! isscalar (p) || ! isscalar (mu) || ! isscalar(sigma))
    [retval, p, mu, sigma] = common_size (p, mu, sigma);
    if (retval > 0)
      error ("loglinv: P, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (p) || iscomplex (mu) || iscomplex (sigma))
    error ("loglinv: P, MU, and SIGMA must not be complex.");
  endif

  ## Check for class type
  if (isa (p, "single") || isa (mu, "single") || isa (sigma, "single"));
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Check for valid points
  k = (p >= 0) & (p <= 1) & (mu >= 0) & (sigma > 0);

  ## Compute the log-logistic iCDF
  x(k) = exp (logit (p(k)) .* sigma(k) + mu(k));

endfunction

%!demo
%! ## Plot various iCDFs from the log-logistic distribution
%! p = 0.001:0.001:0.999;
%! x1 = loglinv (p, log (1), 1/0.5);
%! x2 = loglinv (p, log (1), 1);
%! x3 = loglinv (p, log (1), 1/2);
%! x4 = loglinv (p, log (1), 1/4);
%! x5 = loglinv (p, log (1), 1/8);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r", p, x4, "-c", p, x5, "-m")
%! ylim ([0, 20])
%! grid on
%! legend ({"σ = 2 (β = 0.5)", "σ = 1 (β = 1)", "σ = 0.5 (β = 2)", ...
%!          "σ = 0.25 (β = 4)", "σ = 0.125 (β = 8)"}, "location", "northwest")
%! title ("Log-logistic iCDF")
%! xlabel ("probability")
%! ylabel ("x")
%! text (0.03, 12.5, "μ = 0 (α = 1), values of σ (β) as shown in legend")

## Test output
%!shared p, out1, out2
%! p = [-1, 0, 0.2, 0.5, 0.8, 0.95, 1, 2];
%! out1 = [NaN, 0, 0.25, 1, 4, 19, Inf, NaN];
%! out2 = [NaN, 0, 0.0424732, 2.718282, 173.970037, 18644.695061, Inf, NaN];
%!assert (loglinv (p, 0, 1), out1, 1e-8)
%!assert (loglinv (p, 0, 1), out1, 1e-8)
%!assert (loglinv (p, 1, 3), out2, 1e-6)

## Test class of input preserved
%!assert (class (loglinv (single (1), 2, 3)), "single")
%!assert (class (loglinv (1, single (2), 3)), "single")
%!assert (class (loglinv (1, 2, single (3))), "single")

## Test input validation
%!error<loglinv: function called with too few input arguments.> loglinv (1)
%!error<loglinv: function called with too few input arguments.> loglinv (1, 2)
%!error<loglinv: P, MU, and SIGMA must be of common size or scalars.> ...
%! loglinv (1, ones (2), ones (3))
%!error<loglinv: P, MU, and SIGMA must be of common size or scalars.> ...
%! loglinv (ones (2), 1, ones (3))
%!error<loglinv: P, MU, and SIGMA must be of common size or scalars.> ...
%! loglinv (ones (2), ones (3), 1)
%!error<loglinv: P, MU, and SIGMA must not be complex.> loglinv (i, 2, 3)
%!error<loglinv: P, MU, and SIGMA must not be complex.> loglinv (1, i, 3)
%!error<loglinv: P, MU, and SIGMA must not be complex.> loglinv (1, 2, i)
