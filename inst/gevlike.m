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
## @deftypefn {Function File} {@var{nlogL}, @var{ACOV}, @var{Grad} =} gevlike (@var{params}, @var{data})
## Compute the negative log-likelihood of data under the generalized extreme value (GEV) distribution with given parameter values.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{params} is the 3-parameter vector [@var{k}, @var{sigma}, @var{mu}], where @var{k} is the shape parameter of the GEV distribution, @var{sigma} is the scale parameter of the GEV distribution, and @var{mu} is the location parameter of the GEV distribution.
## @item
## @var{data} is the vector of given values.
##
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{nlogL} is the negative log-likelihood.
## @item
## @var{ACOV} is the 3 by 3 Fisher information matrix (negative second derivative of the log likelihood with respect to the parameter values)
## @var{Grad} is the 3 by 1 gradient vector (first derivative of the log likelihood with respect to the parameter values)
## 
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## x = -5:-1;
## k = -0.2;
## sigma = 0.3;
## mu = 0.5;
## [L, C] = gevlike ([k sigma mu], x);
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
## @seealso{gevcdf, gevfit, gevinv, gevpdf, gevrnd, gevstat}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: Negative log-likelihood for the generalized extreme value distribution

function [nlogL, ACOV, Grad] = gevlike (params, data)

  # Check arguments
  if (nargin != 2)
    print_usage;
  endif
  
  k = params(1);
  sigma = params(2);
  mu = params(3);
  
  #calculate negative log likelihood
  nlogL = sum(gevnll (data, k, sigma, mu) (:));

  if nargout > 1
  	ACOV = -gevfim (data, k, sigma, mu);
    if nargout > 2
  	  Grad = gevgrad (data, k, sigma, mu);
  	endif
  endif

endfunction


function nlogL = gevnll (x, k, sigma, mu)
#internal function to calculate negative log likelihood for gevlike
#no input checking done

  if k == 0
    z = (mu - x) ./ sigma;
    nlogL = exp(z) - z + log(sigma);	
  else
    z = 1 + k .* (x - mu) ./ sigma;
    nlogL = z .^ (-1 ./ k) + (1 + 1 ./ k) .* log(z) + log(sigma);
    nlogL(z <= 0) = Inf;
  endif

endfunction

function G = gevgrad (x, k, sigma, mu)
#calculate the gradient of the negative log likelihood of data x with respect to the parameters of the generalized extreme value distribution for gevlike
#no input checking done

G = ones(3, 1);

if k == 0 ##use the expressions for first derivatives that are the limits as k --> 0
  a = (x - mu) ./ sigma;
  f = exp(a);
  #k
  g = -(2 * x .* (mu .* (1 - f) - sigma .* f) + 2 .* sigma .* mu .* f + (x.^2 + mu.^2).*(f - 1)) ./ (2 * f .* sigma .^ 2);
  G(1) = sum(g(:));
  
  #sigma
  g = ((a ./ f) - a + 1) ./ sigma;
  G(2) = sum(g(:));
  
  #mu
  g = (1 ./ f - 1) ./ sigma;
  G(3) = sum(g(:));
  
  return
endif

z = 1 + k .* (x - mu) ./ sigma;
if any (z <= 0)
  G(:) = 0; #negative log likelihood is locally infinite
  return
end

#k
a = (x - mu) ./ sigma;
b = k .* a + 1;
c = log(b);
d = 1 ./ k + 1;
g = (c ./ k - a ./ b) ./ (k .* b .^ (1/k)) - c ./ (k .^ 2) + a .* d ./ b;
G(1) = sum(g(:));

#sigma
g = (a .* b .^ (-d) - d .* k .* a ./ b + 1) ./ sigma;
G(2) = sum(g(:));

#mu
g = (b .^ (-d) - d .* k ./ b) ./ sigma;
G(3) = sum(g(:));


endfunction

function ACOV = gevfim (x, k, sigma, mu)
#internal function to calculate the (negative) Fisher information matrix for gevlike
#no input checking done

#calculate the various second derivatives (used Maxima to find the expressions)

ACOV = ones(3);

if k == 0 ##use the expressions for second derivatives that are the limits as k --> 0
  #k, k
  a = (x - mu) ./ sigma;
  f = exp(a);
  der = (x .* (24 * mu .^ 2 .* sigma .* (f - 1) + 24 * mu .* sigma .^ 2 .* f - 12 * mu .^ 3) + x .^ 3 .* (8 * sigma .* (f - 1) - 12*mu) + x .^ 2 .* (-12 * sigma .^ 2 .* f + 24 * mu .* sigma .* (1 - f) + 18 * mu .^ 2) - 12 * mu .^ 2 .* sigma .^ 2 .* f + 8 * mu .^ 3 .* sigma .* (1 - f) + 3 * (x .^ 4 + mu .^ 4)) ./ (12 .* f .* sigma .^ 4);
  ACOV(1, 1) = sum(der(:));

  #sigma, sigma
  der = (sigma .^ -2) .* ...
    ((-2*a ./ f) + ...
    a .^ 2 ./ f + ...
    2*a - 1);
  ACOV(2, 2) = sum(der(:));

  #mu, mu
  der = (sigma .^ -2) ./ f;
  ACOV(3, 3) = sum(der(:));

  #k, sigma
  der =  (x .^2 .* (2*sigma .* (f - 1) - 3*mu) + ...
  x .* (-2 * sigma .^ 2 .* f + 4 * mu .* sigma .* (1 - f) + 3 .* mu .^ 2) + ...
  2 * mu .^ 2 .* sigma .* (f - 1) + 2 * mu * sigma .^ 2 * f + x .^ 3 - mu .^ 3) ...
  ./ (2 .* f .* sigma .^ 4);
  ACOV(1, 2) = ACOV(2, 1) = sum(der(:));

  #k, mu
  der = (x .* (2*sigma .* (f - 1) - 2*mu) - 2 * f .* sigma .^ 2 + ...
  2 .* mu .* sigma .* (1 - f) + x .^ 2 + mu .^ 2)...
  ./ (2 .* f .* sigma .^ 3);
  ACOV(1, 3) = ACOV(3, 1) = sum(der(:));

  #sigma, mu
  der = (sigma + (x - sigma - mu) ./ f) ./ (sigma .^ 3);
  ACOV(2, 3) = ACOV(3, 2) = sum(der(:));

  return
endif


#general case

z = 1 + k .* (x - mu) ./ sigma;

if any (z <= 0)
  ACOV(:) = 0; #negative log likelihood is locally infinite
  return
end

#k, k
a = (x - mu) ./ sigma;
b = k .* a + 1;
c = log(b);
d = 1 ./ k + 1;
der = ((((c ./ k.^2) - (a ./ (k .* b))) .^ 2) ./ (b .^ (1 ./ k))) + ...
  ((-2*c ./ k.^3) + (2*a ./ (k.^2 .* b)) + ((a ./ b) .^ 2 ./ k)) ./ (b .^ (1 ./ k)) + ...
  2*c ./ k.^3 - ...
  (2*a ./ (k.^2 .* b)) - (d .* (a ./ b) .^ 2);
ACOV(1, 1) = sum(der(:));

#sigma, sigma
der = (sigma .^ -2) .* (...
  -2*a .* b .^ (-d) + ...
  d .* k .* a .^ 2 .* (b .^ (-d-1)) + ...
  2 .* d .* k .* a ./ b - ...
  d .* (k .* a ./ b) .^ 2 - 1);
ACOV(2, 2) = sum(der(:));

#mu, mu
der = (d .* (sigma .^ -2)) .*  (...
  k .* (b .^ (-d-1)) - ...
  (k ./ b) .^ 2);
ACOV(3, 3) = sum(der(:));

#k, sigma
der = (a ./ sigma) .* (...
(b .^ (-d)) .* (c ./ k  - a ./ b) ./ k + ...
-a .* (b .^ (-d-1)) + ...
((1 ./ k) - d) ./ b + ...
a .* k .* d ./ (b .^ 2));
ACOV(1, 2) = ACOV(2, 1) = sum(der(:));

#k, mu
der = ( (b .^ (-d)) .* (c ./ k  - a ./ b) ./ k - ...
a .* (b .^ (-d-1)) + ...
((1 ./ k) - d) ./ b +
a .* k .* d ./ (b .^ 2)) ./ sigma;
ACOV(1, 3) = ACOV(3, 1) = sum(der(:));

#sigma, mu
der = ( -(b .^ (-d)) + ...
a .* k .* d .* (b .^ (-d-1)) + ...
(d .* k ./ b) - a .* (k./b).^2 .* d) ./ (sigma .^ 2);
ACOV(2, 3) = ACOV(3, 2) = sum(der(:));

endfunction




%!test
%! x = 1;
%! k = 0.2;
%! sigma = 0.3;
%! mu = 0.5;
%! [L, C] = gevlike ([k sigma mu], x);
%! expected_L = 0.75942;
%! expected_C = -[-0.12547 1.77884 1.06731; 1.77884 16.40761 8.48877; 1.06731 8.48877 0.27979];
%! assert (L, expected_L, 0.001);
%! assert (C, expected_C, 0.001);

%!test
%! x = 1;
%! k = 0;
%! sigma = 0.3;
%! mu = 0.5;
%! [L, C] = gevlike ([k sigma mu], x);
%! expected_L = 0.65157;
%! expected_C = -[0.090036 3.41229 2.047337; 3.412229 24.760027 12.510190; 2.047337 12.510190 2.098618];
%! assert (L, expected_L, 0.001);
%! assert (C, expected_C, 0.001);

%!test
%! x = -5:-1;
%! k = -0.2;
%! sigma = 0.3;
%! mu = 0.5;
%! [L, C] = gevlike ([k sigma mu], x);
%! expected_L = 3786.4;
%! expected_C = [-1.4937e+06 1.0083e+06 -6.1837e+04; 1.0083e+06 -8.1138e+05 4.0917e+04; -6.1837e+04 4.0917e+04 -2.0422e+03];
%! assert (L, expected_L, -0.001);
%! assert (C, expected_C, -0.001);
