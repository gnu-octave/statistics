## Copyright (C) 2012-2016 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## @deftypefn {Function File} {@var{paramhat}, @var{paramci} =} gevfit (@var{data}, @var{parmguess})
## Find the maximum likelihood estimator (@var{paramhat}) of the generalized extreme value (GEV) distribution to fit @var{data}.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{data} is the vector of given values.
## @item
## @var{parmguess} is an initial guess for the maximum likelihood parameter vector. If not given, this defaults to @var{k}=0 and @var{sigma}, @var{mu} determined by matching the data mean and standard deviation to their expected values.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{parmhat} is the 3-parameter maximum-likelihood parameter vector [@var{k}; @var{sigma}; @var{mu}], where @var{k} is the shape parameter of the GEV distribution, @var{sigma} is the scale parameter of the GEV distribution, and @var{mu} is the location parameter of the GEV distribution.
## @item
## @var{paramci} has the approximate 95% confidence intervals of the parameter values based on the Fisher information matrix at the maximum-likelihood position.
## 
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## data = 1:50;
## [pfit, pci] = gevfit (data);
## p1 = gevcdf(data,pfit(1),pfit(2),pfit(3));
## plot(data, p1)
## @end group
## @end example
## @seealso{gevcdf, gevinv, gevlike, gevpdf, gevrnd, gevstat}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: Maximum likelihood parameter estimation for the generalized extreme value distribution

function [paramhat, paramci] = gevfit (data, paramguess)

  # Check arguments
  if (nargin < 1)
    print_usage;
  endif

  if (nargin < 2) || isempty(paramguess)
    paramguess = zeros (3, 1);
    paramguess(2) = (sqrt(6)/pi) * std (data);
    paramguess(3) = mean(data) - 0.5772156649*paramguess(2); #expectation involves Euler–Mascheroni constant
  endif

  #cost function to minimize
  f = @(p) gevlike (p, data);
  
  paramhat = fminunc(f, paramguess, optimset("GradObj", "on"));
    
  if nargout > 1
  	[nlogL, ~, ACOV] = gevlike (paramhat, data);
  	param_se = sqrt(diag(inv(ACOV)));
    if any(iscomplex(param_se))
      warning ('gevfit: Fisher information matrix not positive definite; parameter optimization likely did not converge')
      paramci = nan (3, 2);
    else
      paramci(:, 1) = paramhat - 1.96*param_se;
      paramci(:, 2) = paramhat + 1.96*param_se;
    endif
  endif

endfunction

%!test
%! data = 1:50;
%! [pfit, pci] = gevfit (data);
%! expected_p = [-0.44 15.19 21.53]';
%! expected_pu = [-0.13 19.31 26.49]';
%! assert (pfit, expected_p, 0.1);
%! assert (pci(:, 2), expected_pu, 0.1);
