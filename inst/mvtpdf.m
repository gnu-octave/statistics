## Copyright (C) 2015 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## @deftypefn {Function File} {@var{p} =} mvtpdf (@var{x}, @var{sigma}, @var{nu})
## Compute the probability density function of the multivariate Student's t distribution.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{x} are the points at which to find the probability, where each row corresponds
## to an observation. (@var{n} by @var{d} matrix)
##
## @item
## @var{sigma} is the scale matrix. (@var{d} by @var{d} symmetric positive definite matrix)
##
## @item
## @var{nu} is the degrees of freedom. (scalar or @var{n} vector)
##
## @end itemize
##
## The distribution is assumed to be centered (zero mean).
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{p} is the probability density for each row of @var{x}. (@var{n} by 1 vector)
## @end itemize
##
## @subheading Examples
##
## @example
## @group
## x = [1 2];
## sigma = [1.0 0.5; 0.5 1.0];
## nu = 4;
## p = mvtpdf (x, sigma, nu)
## @end group
## @end example
##
## @subheading References
##
## @enumerate
## @item
## Michael Roth, On the Multivariate t Distribution, Technical report from Automatic Control at Linkoepings universitet, @url{http://users.isy.liu.se/en/rt/roth/student.pdf}
## @end enumerate
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: PDF of the multivariate Student's t distribution

function [p] = mvtpdf (x, sigma, nu)

  if (nargin != 3)
    print_usage ();
  endif

  # Dimensions
  d = size (sigma, 1);
  n = size (x, 1);

  # Check parameters
  if (size (x, 2) != d)
    error ("mvtpdf: x must have the same number of columns as sigma");
  endif
  if (! isscalar (nu) && (! isvector (nu) || numel (nu) != n))
    error ("mvtpdf: nu must be a scalar or a vector with the same number of rows as x");
  endif
  if (d < 1 || size (sigma, 2) != d || ! issymmetric (sigma))
    error ("mvtpdf: sigma must be nonempty and symmetric");
  endif

  try
    U = chol (sigma);
  catch
    error ("mvtpdf: sigma must be positive definite");
  end_try_catch
  
  nu = nu(:);
  sqrt_det_sigma = prod(diag(U)); #square root of determinant of sigma

  c = (gamma((nu+d)/2) ./ gamma(nu/2)) ./ (sqrt_det_sigma * (nu*pi).^(d/2)); #scale factor for PDF
  p = c ./ ((1 + sumsq(U' \ x') ./ nu') .^ ((nu' + d)/2))'; #note: sumsq(U' \ x') is equivalent to the quadratic form x*inv(sigma)*x'


endfunction

#test results verified with R mvtnorm package dmvt function
%!assert (mvtpdf ([0 0], eye(2), 1), 0.1591549, 1E-7) #dmvt(x = c(0,0), sigma = diag(2), log = FALSE)
%!assert (mvtpdf ([1 0], [1 0.5; 0.5 1], 2), 0.06615947, 1E-7) #dmvt(x = c(1,0), sigma = matrix(c(1, 0.5, 0.5, 1), nrow=2, ncol=2), df = 2, log = FALSE)
%!assert (mvtpdf ([1 0.4 0; 1.2 0.5 0.5; 1.4 0.6 1], [1 0.5 0.3; 0.5 1 0.6; 0.3 0.6 1], [5 6 7]), [0.04713313 0.03722421 0.02069011]', 1E-7) #dmvt(x = c(1,0.4,0), sigma = matrix(c(1, 0.5, 0.3, 0.5, 1, 0.6, 0.3, 0.6, 1), nrow=3, ncol=3), df = 5, log = FALSE); dmvt(x = c(1.2,0.5,0.5), sigma = matrix(c(1, 0.5, 0.3, 0.5, 1, 0.6, 0.3, 0.6, 1), nrow=3, ncol=3), df = 6, log = FALSE); dmvt(x = c(1.4,0.6,1), sigma = matrix(c(1, 0.5, 0.3, 0.5, 1, 0.6, 0.3, 0.6, 1), nrow=3, ncol=3), df = 7, log = FALSE)
