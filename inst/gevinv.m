## Copyright (C) 2012 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## @deftypefn {Function File} {@var{X} =} gevinv (@var{P}, @var{k}, @var{sigma}, @var{mu})
## Compute a desired quantile (inverse CDF) of the generalized extreme value (GEV) distribution.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{P} is the desired quantile of the GEV distribution. (Between 0 and 1.)
## @item
## @var{k} is the shape parameter of the GEV distribution. (Also denoted gamma or xi.)
## @item
## @var{sigma} is the scale parameter of the GEV distribution. The elements
## of @var{sigma} must be positive.
## @item
## @var{mu} is the location parameter of the GEV distribution.
## @end itemize
## The inputs must be of common size, or some of them must be scalar.
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{X} is the value corresponding to each quantile of the GEV distribution
## @end itemize
## @subheading References
##
## @enumerate
## @item
## Rolf-Dieter Reiss and Michael Thomas. @cite{Statistical Analysis of Extreme Values with Applications to Insurance, Finance, Hydrology and Other Fields}. Chapter 1, pages 16-17, Springer, 2007.
## @item
## J. R. M. Hosking (2012). @cite{L-moments}. R package, version 1.6. URL: http://CRAN.R-project.org/package=lmom.
##
## @end enumerate
## @seealso{gevcdf, gevfit, gevlike, gevpdf, gevrnd, gevstat}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: Inverse CDF of the generalized extreme value distribution

function [X] = gevinv (P, k = 0, sigma = 1, mu = 0)

  [retval, P, k, sigma, mu] = common_size (P, k, sigma, mu);
  if (retval > 0)
    error ("gevinv: inputs must be of common size or scalars");
  endif

  X = P;
  
  llP = log(-log(P));
  kllP = k .* llP;
  
  ii = (abs(kllP) < 1E-4); #use the Taylor series expansion of the exponential to avoid roundoff error or dividing by zero when k is small
  X(ii) = mu(ii) - sigma(ii) .* llP(ii) .* (1 - kllP(ii) .* (1 - kllP(ii)));
  X(~ii) = mu(~ii) + (sigma(~ii) ./ k(~ii)) .* (exp(-kllP(~ii)) - 1);

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

