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
## @deftypefn {Function File} {[@var{m}, @var{v}] =} gevstat (@var{k}, @var{sigma}, @var{mu})
## Compute the mean and variance of the generalized extreme value (GEV) distribution.
##
## @subheading Arguments
##
## @itemize @bullet
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
## @var{m} is the mean of the GEV distribution
##
## @item
## @var{v} is the variance of the GEV distribution
## @end itemize
## @seealso{gevcdf, gevfit, gevinv, gevlike, gevpdf, gevrnd}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: Moments of the generalized extreme value distribution

function [m, v] = gevstat (k, sigma, mu)

  # Check arguments
  if (nargin < 3)
    print_usage ();
  endif

  if (isempty (k) || isempty (sigma) || isempty (mu) || ~ismatrix (k) || ~ismatrix (sigma) || ~ismatrix (mu))
    error ("gevstat: inputs must be numeric matrices");
  endif

  [retval, k, sigma, mu] = common_size (k, sigma, mu);
  if (retval > 0)
    error ("gevstat: inputs must be of common size or scalars");
  endif

  eg = 0.57721566490153286; %Euler-Mascheroni constant

  m = v = k;
  
  #find the mean  
  m(k >= 1) = Inf;
  m(k == 0) = mu(k == 0) + eg*sigma(k == 0);
  m(k < 1 & k ~= 0) = mu(k < 1 & k ~= 0) + sigma(k < 1 & k ~= 0) .* (gamma(1-k(k < 1 & k ~= 0)) - 1) ./ k(k < 1 & k ~= 0);
  
  #find the variance
  v(k >= 0.5) = Inf;
  v(k == 0) = (pi^2 / 6) * sigma(k == 0) .^ 2;
  v(k < 0.5 & k ~= 0) = (gamma(1-2*k(k < 0.5 & k ~= 0)) - gamma(1-k(k < 0.5 & k ~= 0)).^2) .* (sigma(k < 0.5 & k ~= 0) ./ k(k < 0.5 & k ~= 0)) .^ 2;  
  
  

endfunction

%!test
%! k = [-1 -0.5 0 0.2 0.4 0.5 1];
%! sigma = 2;
%! mu = 1;
%! [m, v] = gevstat (k, sigma, mu);
%! expected_m = [1 1.4551 2.1544 2.6423   3.4460   4.0898      Inf];
%! expected_v = [4 3.4336 6.5797 13.3761   59.3288       Inf       Inf];
%! assert (m, expected_m, -0.001);
%! assert (v, expected_v, -0.001);
