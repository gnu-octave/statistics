## Copyright (C) 2013 Nir Krakauer <nkrakauer@ccny.cuny.edu>
##
## This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with Octave; see the file COPYING.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {} @var{y} = wishpdf (@var{W}, @var{Sigma}, @var{df}, @var{log_y}=false)
## Compute the probability density function of the Wishart distribution
##
## Inputs: A @var{p} x @var{p} matrix @var{W} where to find the PDF. The @var{p} x @var{p} positive definite matrix @var{Sigma} and scalar degrees of freedom parameter @var{df} characterizing the Wishart distribution. (For the density to be finite, need @var{df} > (@var{p} - 1).)
## If the flag @var{log_y} is set, return the log probability density -- this helps avoid underflow when the numerical value of the density is very small
##
## Output: @var{y} is the probability density of Wishart(@var{Sigma}, @var{df}) at @var{W}.
## 
## @seealso{wishrnd, iwishpdf}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: Compute the probability density function of the Wishart distribution

function [y] = wishpdf(W, Sigma, df, log_y=false)

if (nargin < 3)
  print_usage ();
endif

p = size(Sigma, 1);

if (df <= (p - 1))
  error('df too small, no finite densities exist')
endif 

#calculate the logarithm of G_d(df/2), the multivariate gamma function
g = (p * (p-1) / 4) * log(pi);
for i = 1:p
  g = g + log(gamma((df + (1 - i))/2)); #using lngamma_gsl(.) from the gsl package instead of log(gamma(.)) might help avoid underflow/overflow 
endfor

C = chol(Sigma);

#use formulas for determinant of positive definite matrix for better efficiency and numerical accuracy
logdet_W = 2*sum(log(diag(chol(W))));
logdet_Sigma = 2*sum(log(diag(C)));

y = -(df*p)/2 * log(2) - (df/2)*logdet_Sigma - g + ((df - p - 1)/2)*logdet_W - trace(chol2inv(C)*W)/2;

if ~log_y
  y = exp(y);
endif


endfunction

##test results cross-checked against dwish function in R MCMCpack library 
%!assert(wishpdf(4, 3, 3.1), 0.07702496, 1E-7);
%!assert(wishpdf([2 -0.3;-0.3 4], [1 0.3;0.3 1], 4), 0.004529741, 1E-7);
%!assert(wishpdf([6 2 5; 2 10 -5; 5 -5 25], [9 5 5; 5 10 -8; 5 -8 22], 5.1), 4.474865e-10, 1E-15);

%% Test input validation
%!error wishpdf ()
%!error wishpdf (1, 2)
%!error wishpdf (1, 2, 0)

%!error wishpdf (1, 2)
