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
## @deftypefn {Function File} {@var{p} =} gevcdf (@var{x}, @var{k}, @var{sigma}, @var{mu})
## Compute the cumulative distribution function of the generalized extreme value (GEV) distribution.
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
## @var{p} is the cumulative distribution of the GEV distribution at each
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
## y = gevcdf (x, k, sigma, mu)
## @end group
##
## @group
## y = gevcdf (x, k, 0.5, mu)
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
## @seealso{gevfit, gevinv, gevlike, gevpdf, gevrnd, gevstat}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: CDF of the generalized extreme value distribution

function p = gevcdf (x, k, sigma, mu)

  # Check arguments
  if (nargin != 4)
    print_usage ();
  endif

  if (isempty (x) || isempty (k) || isempty (sigma) || isempty (mu) || ~ismatrix (x) || ~ismatrix (k) || ~ismatrix (sigma) || ~ismatrix (mu))
    error ("gevcdf: inputs must be a numeric matrices");
  endif

  [retval, x, k, sigma, mu] = common_size (x, k, sigma, mu);
  if (retval > 0)
    error ("gevcdf: inputs must be of common size or scalars");
  endif

  z = 1 + k .* (x - mu) ./ sigma;

  # Calculate pdf
  p = exp(-(z .^ (-1 ./ k)));

  p(z <= 0 & x < mu) = 0;
  p(z <= 0 & x > mu) = 1;  
  
  inds = (k == 0); %use a different formula
  if any(inds)
    z = (mu(inds) - x(inds)) ./ sigma(inds);
    p(inds) = exp(-exp(z));
  endif
  

endfunction

%!test
%! x = 0:0.5:2.5;
%! sigma = 1:6;
%! k = 1;
%! mu = 0;
%! p = gevcdf (x, k, sigma, mu);
%! expected_p = [0.36788   0.44933   0.47237   0.48323   0.48954   0.49367];
%! assert (p, expected_p, 0.001);

%!test
%! x = -0.5:0.5:2.5;
%! sigma = 0.5;
%! k = 1;
%! mu = 0;
%! p = gevcdf (x, k, sigma, mu);
%! expected_p = [0   0.36788   0.60653   0.71653   0.77880   0.81873   0.84648];
%! assert (p, expected_p, 0.001);

%!test #check for continuity for k near 0
%! x = 1;
%! sigma = 0.5;
%! k = -0.03:0.01:0.03;
%! mu = 0;
%! p = gevcdf (x, k, sigma, mu);
%! expected_p = [0.88062   0.87820   0.87580   0.87342   0.87107   0.86874   0.86643];
%! assert (p, expected_p, 0.001);
