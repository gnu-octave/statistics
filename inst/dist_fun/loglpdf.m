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
## @deftypefn  {statistics} {@var{y} =} loglpdf (@var{x}, @var{mu}, @var{sigma})
##
## Loglogistic probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the loglogistic distribution with mean parameter @var{mu} and scale
## parameter @var{sigma}.  The size of @var{y} is the common size of @var{x},
## @var{mu}, and @var{sigma}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Mean of logarithmic values @var{mu} must be a non-negative real value, scale
## parameter of logarithmic values @var{sigma} must be a positive real value and
## @var{x} is supported in the range @math{[0,Inf)}, otherwise 0 is returned.
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
## @seealso{loglcdf, loglinv, loglrnd, loglfit, logllike, loglstat}
## @end deftypefn

function y = loglpdf (x, mu, sigma)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("loglpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, MU, and SIGMA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar(sigma))
    [retval, x, mu, sigma] = common_size (x, mu, sigma);
    if (retval > 0)
      error ("loglpdf: X, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("loglpdf: X, MU, and SIGMA must not be complex.");
  endif

  ## Check for invalid points
  mu(mu < 0) = NaN;
  sigma(sigma <= 0) = NaN;

  ## Compute log-logistic PDF
  a = exp (mu);
  b = 1./ sigma;
  y = ((b ./ a) .* (x ./ a) .^ (b - 1)) ./ ((1 + (x ./ a) .^ b) .^ 2);
  y(x <= 0) = 0;

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"));
    y = cast (y, "single");
  endif

endfunction

%!demo
%! ## Plot various PDFs from the log-logistic distribution
%! x = 0.001:0.001:2;
%! y1 = loglpdf (x, log (1), 1/0.5);
%! y2 = loglpdf (x, log (1), 1);
%! y3 = loglpdf (x, log (1), 1/2);
%! y4 = loglpdf (x, log (1), 1/4);
%! y5 = loglpdf (x, log (1), 1/8);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", x, y4, "-c", x, y5, "-m")
%! grid on
%! ylim ([0,3])
%! legend ({"σ = 2 (β = 0.5)", "σ = 1 (β = 1)", "σ = 0.5 (β = 2)", ...
%!          "σ = 0.25 (β = 4)", "σ = 0.125 (β = 8)"}, "location", "northeast")
%! title ("Log-logistic PDF")
%! xlabel ("values in x")
%! ylabel ("density")
%! text (0.1, 2.8, "μ = 0 (α = 1), values of σ (β) as shown in legend")

## Test output
%!shared out1, out2
%! out1 = [0, 0, 1, 0.2500, 0.1111, 0.0625, 0.0400, 0.0278, 0];
%! out2 = [0, 0, 0.0811, 0.0416, 0.0278, 0.0207, 0.0165, 0];
%!assert (loglpdf ([-1,0,realmin,1:5,Inf], 0, 1), out1, 1e-4)
%!assert (loglpdf ([-1,0,realmin,1:5,Inf], 0, 1), out1, 1e-4)
%!assert (loglpdf ([-1:5,Inf], 1, 3), out2, 1e-4)

## Test class of input preserved
%!assert (class (loglpdf (single (1), 2, 3)), "single")
%!assert (class (loglpdf (1, single (2), 3)), "single")
%!assert (class (loglpdf (1, 2, single (3))), "single")

## Test input validation
%!error<loglpdf: function called with too few input arguments.> loglpdf (1)
%!error<loglpdf: function called with too few input arguments.> loglpdf (1, 2)
%!error<loglpdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! loglpdf (1, ones (2), ones (3))
%!error<loglpdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! loglpdf (ones (2), 1, ones (3))
%!error<loglpdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! loglpdf (ones (2), ones (3), 1)
%!error<loglpdf: X, MU, and SIGMA must not be complex.> loglpdf (i, 2, 3)
%!error<loglpdf: X, MU, and SIGMA must not be complex.> loglpdf (1, i, 3)
%!error<loglpdf: X, MU, and SIGMA must not be complex.> loglpdf (1, 2, i)
