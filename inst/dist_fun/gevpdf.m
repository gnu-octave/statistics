## Copyright (C) 2012 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## @deftypefn  {statistics} {@var{y} =} gevpdf (@var{x}, @var{k}, @var{sigma}, @var{mu})
##
## Generalized extreme value (GEV) probability density function (PDF).
##
## For each element of @var{x}, compute the probability density function (PDF)
## of the GEV distribution with shape parameter @var{k}, scale parameter
## @var{sigma}, and location parameter @var{mu}.  The size of @var{y} is the
## common size of @var{x}, @var{k}, @var{sigma}, and @var{mu}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## When @qcode{@var{k} < 0}, the GEV is the type III extreme value distribution.
## When @qcode{@var{k} > 0}, the GEV distribution is the type II, or Frechet,
## extreme value distribution.  If @var{W} has a Weibull distribution as
## computed by the @code{wblcdf} function, then @qcode{-@var{W}} has a type III
## extreme value distribution and @qcode{1/@var{W}} has a type II extreme value
## distribution.  In the limit as @var{k} approaches @qcode{0}, the GEV is the
## mirror image of the type I extreme value distribution as computed by the
## @code{evcdf} function.
##
## The mean of the GEV distribution is not finite when @qcode{@var{k} >= 1}, and
## the variance is not finite when @qcode{@var{k} >= 1/2}.  The GEV distribution
## has positive density only for values of @var{x} such that
## @qcode{@var{k} * (@var{x} - @var{mu}) / @var{sigma} > -1}.
##
## Further information about the generalized extreme value distribution can be
## found at
## @url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}
##
## @subheading References
## @enumerate
## @item
## Rolf-Dieter Reiss and Michael Thomas. @cite{Statistical Analysis of Extreme
## Values with Applications to Insurance, Finance, Hydrology and Other Fields}.
## Chapter 1, pages 16-17, Springer, 2007.
## @end enumerate
##
## @seealso{gevcdf, gevinv, gevrnd, gevfit, gevlike, gevstat}
## @end deftypefn

function y = gevpdf (x, k, sigma, mu)

  ## Check for valid number of input arguments
  if (nargin < 4)
    error ("gevpdf: function called with too few input arguments.");
  endif

  ## Check for common size of X, K, SIGMA, and MU
  [retval, x, k, sigma, mu] = common_size (x, k, sigma, mu);
  if (retval > 0)
    error ("gevpdf: X, K, SIGMA, and MU must be of common size or scalars.");
  endif

  ## Check for X, K, SIGMA, and MU being reals
  if (iscomplex (x) || iscomplex (k) || iscomplex (sigma) || iscomplex (mu))
    error ("gevpdf: X, K, SIGMA, and MU must not be complex.");
  endif

  z = 1 + k .* (x - mu) ./ sigma;

  ## Calculate generalized extreme value PDF
  y = exp(-(z .^ (-1 ./ k))) .* (z .^ (-1 - 1 ./ k)) ./ sigma;

  y(z <= 0) = 0;

  ## Use a different formula if k is very close to zero
  inds = (abs (k) < (eps^0.7));
  if (any (inds))
    z = (mu(inds) - x(inds)) ./ sigma(inds);
    y(inds) = exp (z - exp (z)) ./ sigma(inds);
  endif

endfunction

%!demo
%! ## Plot various PDFs from the generalized extreme value distribution
%! x = -1:0.001:10;
%! y1 = gevpdf (x, 1, 1, 1);
%! y2 = gevpdf (x, 0.5, 1, 1);
%! y3 = gevpdf (x, 1, 1, 5);
%! y4 = gevpdf (x, 1, 2, 5);
%! y5 = gevpdf (x, 1, 5, 5);
%! y6 = gevpdf (x, 1, 0.5, 5);
%! plot (x, y1, "-b", x, y2, "-g", x, y3, "-r", ...
%!       x, y4, "-c", x, y5, "-m", x, y6, "-k")
%! grid on
%! xlim ([-1, 10])
%! ylim ([0, 1.1])
%! legend ({"ξ = 1, σ = 1, μ = 1", "ξ = 0.5, σ = 1, μ = 1", ...
%!          "ξ = 1, σ = 1, μ = 5", "ξ = 1, σ = 2, μ = 5", ...
%!          "ξ = 1, σ = 5, μ = 5", "ξ = 1, σ = 0.5, μ = 5"}, ...
%!         "location", "northeast")
%! title ("Generalized extreme value PDF")
%! xlabel ("values in x")
%! ylabel ("density")

## Test output
%!test
%! x = 0:0.5:2.5;
%! sigma = 1:6;
%! k = 1;
%! mu = 0;
%! y = gevpdf (x, k, sigma, mu);
%! expected_y = [0.367879   0.143785   0.088569   0.063898   0.049953   0.040997];
%! assert (y, expected_y, 0.001);
%!test
%! x = -0.5:0.5:2.5;
%! sigma = 0.5;
%! k = 1;
%! mu = 0;
%! y = gevpdf (x, k, sigma, mu);
%! expected_y = [0 0.735759   0.303265   0.159229   0.097350   0.065498   0.047027];
%! assert (y, expected_y, 0.001);
%!test # check for continuity for k near 0
%! x = 1;
%! sigma = 0.5;
%! k = -0.03:0.01:0.03;
%! mu = 0;
%! y = gevpdf (x, k, sigma, mu);
%! expected_y = [0.23820   0.23764   0.23704   0.23641   0.23576   0.23508   0.23438];
%! assert (y, expected_y, 0.001);

## Test input validation
%!error<gevpdf: function called with too few input arguments.> gevpdf ()
%!error<gevpdf: function called with too few input arguments.> gevpdf (1)
%!error<gevpdf: function called with too few input arguments.> gevpdf (1, 2)
%!error<gevpdf: function called with too few input arguments.> gevpdf (1, 2, 3)
%!error<gevpdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevpdf (ones (3), ones (2), ones(2), ones(2))
%!error<gevpdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevpdf (ones (2), ones (3), ones(2), ones(2))
%!error<gevpdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevpdf (ones (2), ones (2), ones(3), ones(2))
%!error<gevpdf: X, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevpdf (ones (2), ones (2), ones(2), ones(3))
%!error<gevpdf: X, K, SIGMA, and MU must not be complex.> gevpdf (i, 2, 3, 4)
%!error<gevpdf: X, K, SIGMA, and MU must not be complex.> gevpdf (1, i, 3, 4)
%!error<gevpdf: X, K, SIGMA, and MU must not be complex.> gevpdf (1, 2, i, 4)
%!error<gevpdf: X, K, SIGMA, and MU must not be complex.> gevpdf (1, 2, 3, i)
