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
## @deftypefn  {statistics} {@var{p} =} loglcdf (@var{x}, @var{mu}, @var{sigma})
## @deftypefnx {statistics} {@var{p} =} loglcdf (@var{x}, @var{mu}, @var{sigma}, @qcode{"upper"})
##
## Loglogistic cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the loglogistic distribution with mean parameter @var{mu} and scale
## parameter @var{sigma}.  The size of @var{p} is the common size of @var{x},
## @var{mu}, and @var{sigma}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Both parameters, @var{mu} and @var{sigma}, must be positive reals and @var{x}
## is supported in the range @math{[0,inf)}, otherwise @qcode{NaN} is returned.
##
## @code{@var{p} = loglcdf (@var{x}, @var{mu}, @var{sigma}, "upper")} computes
## the upper tail probability of the log-logistic distribution with parameters
## @var{mu} and @var{sigma}, at the values in @var{x}.
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
## @seealso{loglinv, loglpdf, loglrnd, loglfit, logllike, loglstat}
## @end deftypefn

function p = loglcdf (x, mu, sigma, uflag)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("loglcdf: function called with too few input arguments.");
  endif

  ## Check for valid "upper" flag
  if (nargin > 3)
    if (! strcmpi (uflag, "upper"))
      error ("loglcdf: invalid argument for upper tail.");
    else
      uflag = true;
    endif
  else
    uflag = false;
  endif

  ## Check for common size of X, MU, and SIGMA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(sigma))
    [retval, x, mu, sigma] = common_size (x, mu, sigma);
    if (retval > 0)
      error ("loglcdf: X, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("loglcdf: X, MU, and SIGMA must not be complex.");
  endif

  ## Check for invalid points
  mu(mu < 0) = NaN;
  sigma(sigma <= 0) = NaN;
  x(x <= 0) = realmin;

  ## Compute log-logistic CDF
  z = (log (x) - mu) ./ sigma;
  if (uflag)
    p = 1 ./ (1 + exp (z));
  else
    p = 1 ./ (1 + exp (-z));
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"));
    p = cast (p, "single");
  endif

endfunction

%!demo
%! ## Plot various CDFs from the log-logistic distribution
%! x = 0:0.001:2;
%! p1 = loglcdf (x, log (1), 1/0.5);
%! p2 = loglcdf (x, log (1), 1);
%! p3 = loglcdf (x, log (1), 1/2);
%! p4 = loglcdf (x, log (1), 1/4);
%! p5 = loglcdf (x, log (1), 1/8);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c", x, p5, "-m")
%! legend ({"σ = 2 (β = 0.5)", "σ = 1 (β = 1)", "σ = 0.5 (β = 2)", ...
%!          "σ = 0.25 (β = 4)", "σ = 0.125 (β = 8)"}, "location", "northwest")
%! grid on
%! title ("Log-logistic CDF")
%! xlabel ("values in x")
%! ylabel ("probability")
%! text (0.05, 0.64, "μ = 0 (α = 1), values of σ (β) as shown in legend")

## Test output
%!shared out1, out2
%! out1 = [0, 0.5, 0.66666667, 0.75, 0.8, 0.83333333];
%! out2 = [0, 0.4174, 0.4745, 0.5082, 0.5321, 0.5506];
%!assert (loglcdf ([0:5], 0, 1), out1, 1e-8)
%!assert (loglcdf ([0:5], 0, 1, "upper"), 1 - out1, 1e-8)
%!assert (loglcdf ([0:5], 0, 1), out1, 1e-8)
%!assert (loglcdf ([0:5], 0, 1, "upper"), 1 - out1, 1e-8)
%!assert (loglcdf ([0:5], 1, 3), out2, 1e-4)
%!assert (loglcdf ([0:5], 1, 3, "upper"), 1 - out2, 1e-4)

## Test class of input preserved
%!assert (class (loglcdf (single (1), 2, 3)), "single")
%!assert (class (loglcdf (1, single (2), 3)), "single")
%!assert (class (loglcdf (1, 2, single (3))), "single")

## Test input validation
%!error<loglcdf: function called with too few input arguments.> loglcdf (1)
%!error<loglcdf: function called with too few input arguments.> loglcdf (1, 2)
%!error<loglcdf: invalid argument for upper tail.> ...
%! loglcdf (1, 2, 3, 4)
%!error<loglcdf: invalid argument for upper tail.> ...
%! loglcdf (1, 2, 3, "uper")
%!error<loglcdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! loglcdf (1, ones (2), ones (3))
%!error<loglcdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! loglcdf (1, ones (2), ones (3), "upper")
%!error<loglcdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! loglcdf (ones (2), 1, ones (3))
%!error<loglcdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! loglcdf (ones (2), 1, ones (3), "upper")
%!error<loglcdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! loglcdf (ones (2), ones (3), 1)
%!error<loglcdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! loglcdf (ones (2), ones (3), 1, "upper")
%!error<loglcdf: X, MU, and SIGMA must not be complex.> loglcdf (i, 2, 3)
%!error<loglcdf: X, MU, and SIGMA must not be complex.> loglcdf (i, 2, 3, "upper")
%!error<loglcdf: X, MU, and SIGMA must not be complex.> loglcdf (1, i, 3)
%!error<loglcdf: X, MU, and SIGMA must not be complex.> loglcdf (1, i, 3, "upper")
%!error<loglcdf: X, MU, and SIGMA must not be complex.> loglcdf (1, 2, i)
%!error<loglcdf: X, MU, and SIGMA must not be complex.> loglcdf (1, 2, i, "upper")
