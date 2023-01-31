## Copyright (C) 2012 Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} @var{x} = gevinv (@var{p}, @var{k}, @var{sigma}, @var{mu})
##
## Inverse of the generalized extreme value (GEV) cumulative distribution
## function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF)
## at @var{p} of the GEV distribution with shape parameter @var{k}, scale
## parameter @var{sigma}, and location parameter @var{mu}.  The size of @var{x}
## is the common size of the input arguments.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## Default values for K, SIGMA, and MU are 0, 1, and 0, respectively.
##
## @subheading References
##
## @enumerate
## @item
## Rolf-Dieter Reiss and Michael Thomas. @cite{Statistical Analysis of Extreme
## Values with Applications to Insurance, Finance, Hydrology and Other Fields}.
## Chapter 1, pages 16-17, Springer, 2007.
## @end enumerate
##
## @seealso{gevcdf, gevpdf, gevrnd, gevfit, gevlike, gevstat}
## @end deftypefn

function x = gevinv (p, k = 0, sigma = 1, mu = 0)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("gevinv: too few input arguments.");
  endif

  ## Check for common size of P, K, SIGMA, and MU
  [retval, p, k, sigma, mu] = common_size (p, k, sigma, mu);
  if (retval > 0)
    error ("gevinv: P, K, SIGMA, and MU must be of common size or scalars.");
  endif

  x = p;

  llP = log (-log (p));
  kllP = k .* llP;

  ## Use the Taylor series expansion of the exponential to
  ## avoid roundoff error or dividing by zero when k is small
  ii = (abs(kllP) < 1E-4);
  x(ii) = mu(ii) - sigma(ii) .* llP(ii) .* (1 - kllP(ii) .* (1 - kllP(ii)));
  x(~ii) = mu(~ii) + (sigma(~ii) ./ k(~ii)) .* (exp(-kllP(~ii)) - 1);

endfunction

%!test
%! p = 0.1:0.1:0.9;
%! k = 0;
%! sigma = 1;
%! mu = 0;
%! x = gevinv (p, k, sigma, mu);
%! c = gevcdf(x, k, sigma, mu);
%! assert (c, p, 0.001);

%!test
%! p = 0.1:0.1:0.9;
%! k = 1;
%! sigma = 1;
%! mu = 0;
%! x = gevinv (p, k, sigma, mu);
%! c = gevcdf(x, k, sigma, mu);
%! assert (c, p, 0.001);

%!test
%! p = 0.1:0.1:0.9;
%! k = 0.3;
%! sigma = 1;
%! mu = 0;
%! x = gevinv (p, k, sigma, mu);
%! c = gevcdf(x, k, sigma, mu);
%! assert (c, p, 0.001);

## Test input validation
%!error<gevinv: too few input arguments.> gevinv ()
%!error<gevinv: P, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevinv (ones (3), ones (2))
%!error<gevinv: P, K, SIGMA, and MU must be of common size or scalars.> ...
%! gevinv (ones (3), ones (2), 1)

