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
## @deftypefn {Function File} {@var{paramhat}, @var{paramci} =} gevfit_lmom (@var{data})
## Find an estimator (@var{paramhat}) of the generalized extreme value (GEV) distribution fitting @var{data} using the method of L-moments.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{data} is the vector of given values.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{parmhat} is the 3-parameter maximum-likelihood parameter vector [@var{k}; @var{sigma}; @var{mu}], where @var{k} is the shape parameter of the GEV distribution, @var{sigma} is the scale parameter of the GEV distribution, and @var{mu} is the location parameter of the GEV distribution.
## @item
## @var{paramci} has the approximate 95% confidence intervals of the parameter values (currently not implemented).
## 
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## data = gevrnd (0.1, 1, 0, 100, 1);
## [pfit, pci] = gevfit_lmom (data);
## p1 = gevcdf (data,pfit(1),pfit(2),pfit(3));
## [f, x] = ecdf (data);
## plot(data, p1, 's', x, f)
## @end group
## @end example
## @seealso{gevfit}
## @subheading References
##
## @enumerate
## @item
## Ailliot, P.; Thompson, C. & Thomson, P. Mixed methods for fitting the GEV distribution, Water Resources Research, 2011, 47, W05551
##
## @end enumerate
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: L-moments parameter estimation for the generalized extreme value distribution

function [paramhat, paramci] = gevfit_lmom (data)

  # Check arguments
  if (nargin < 1)
    print_usage;
  endif
  
  # find the L-moments
  data = sort (data(:))';
  n = numel(data);
  L1 = mean(data);
  L2 = sum(data .* (2*(1:n) - n - 1)) / (2*nchoosek(n, 2)); # or mean(triu(data' - data, 1, 'pack')) / 2;
  b = bincoeff((1:n) - 1, 2);
  L3 = sum(data .* (b - 2 * ((1:n) - 1) .* (n - (1:n)) + fliplr(b))) / (3*nchoosek(n, 3));

  #match the moments to the GEV distribution
  #first find k based on L3/L2
  f = @(k) (L3/L2 + 3)/2 - limdiv((1 - 3^(k)), (1 - 2^(k)));
  k = fzero(f, 0);
  
  #next find sigma and mu given k
  if abs(k) < 1E-8
    sigma = L2 / log(2);
    eg = 0.57721566490153286; %Euler-Mascheroni constant
    mu = L1 - sigma * eg;
  else
    sigma = -k*L2 / (gamma(1 - k) * (1 - 2^(k)));
    mu = L1 - sigma * ((gamma(1 - k) - 1) / k);
  endif
  
  paramhat = [k; sigma; mu];
  
  if nargout > 1
    paramci = NaN;
  endif
endfunction

#internal function to accurately evaluate (1 - 3^k)/(1 - 2^k) in the limit as k --> 0
function c = limdiv(a, b)
  # c = ifelse (abs(b) < 1E-8, log(3)/log(2), a ./ b);
  if abs(b) < 1E-8
    c = log(3)/log(2);
  else
    c = a / b;
  endif
endfunction


%!xtest <31070>
%! data = 1:50;
%! [pfit, pci] = gevfit_lmom (data);
%! expected_p = [-0.28 15.01 20.22]';
%! assert (pfit, expected_p, 0.1);
