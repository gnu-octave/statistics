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
## @deftypefn {Function File} {@var{y} =} gevpdf (@var{x}, @var{k}, @var{sigma}, @var{mu})
## Compute the probability density function of the generalized extreme value (GEV) distribution.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{x} is the support.
##
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
## @var{y} is the probability density of the GEV distribution at each
## element of @var{x} and corresponding parameter values.
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## x = 0:0.5:2.5;
## sigma = 1:6;
## k = 1;
## mu = 0;
## y = gevpdf (x, k, sigma, mu)
## @end group
##
## @group
## y = gevpdf (x, k, 0.5, mu)
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Rolf-Dieter Reiss and Michael Thomas. @cite{Statistical Analysis of Extreme Values with Applications to Insurance, Finance, Hydrology and Other Fields}. Chapter 1, pages 16-17, Springer, 2007.
##
## @end enumerate
## @seealso{gevcdf, gevfit, gevinv, gevlike, gevrnd, gevstat}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: PDF of the generalized extreme value distribution

function y = gevpdf (x, k, sigma, mu)

  # Check arguments
  if (nargin != 4)
    print_usage ();
  endif

  if (isempty (x) || isempty (k) || isempty (sigma) || isempty (mu) || ~ismatrix (x) || ~ismatrix (k) || ~ismatrix (sigma) || ~ismatrix (mu))
    error ("gevpdf: inputs must be numeric matrices");
  endif

  [retval, x, k, sigma, mu] = common_size (x, k, sigma, mu);
  if (retval > 0)
    error ("gevpdf: inputs must be of common size or scalars");
  endif

  z = 1 + k .* (x - mu) ./ sigma;

  # Calculate pdf
  y = exp(-(z .^ (-1 ./ k))) .* (z .^ (-1 - 1 ./ k)) ./ sigma;

  y(z <= 0) = 0;
  
  inds = (abs (k) < (eps^0.7)); %use a different formula if k is very close to zero
  if any(inds)
    z = (mu(inds) - x(inds)) ./ sigma(inds);
    y(inds) = exp(z-exp(z)) ./ sigma(inds);
  endif
  

endfunction

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

%!test #check for continuity for k near 0
%! x = 1;
%! sigma = 0.5;
%! k = -0.03:0.01:0.03;
%! mu = 0;
%! y = gevpdf (x, k, sigma, mu);
%! expected_y = [0.23820   0.23764   0.23704   0.23641   0.23576   0.23508   0.23438];
%! assert (y, expected_y, 0.001);
