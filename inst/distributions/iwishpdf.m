## Copyright (C) 2013 Nir Krakauer <nkrakauer@ccny.cuny.edu>
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {} @var{y} = iwishpdf (@var{W}, @var{Tau}, @var{df}, @var{log_y}=false)
## Compute the probability density function of the inverse Wishart distribution.
##
## Inputs: A @var{p} x @var{p} matrix @var{W} where to find the PDF and the
## @var{p} x @var{p} positive definite scale matrix @var{Tau} and scalar degrees
## of freedom parameter @var{df} characterizing the inverse Wishart distribution.
## (For the density to be finite, need @var{df} > (@var{p} - 1).)
## If the flag @var{log_y} is set, return the log probability density -- this
## helps avoid underflow when the numerical value of the density is very small.
##
## Output: @var{y} is the probability density of Wishart(@var{Sigma}, @var{df})
## at @var{W}.
## 
## @seealso{iwishrnd, wishpdf}
## @end deftypefn

function [y] = iwishpdf(W, Tau, df, log_y=false)

  if (nargin < 3)
    print_usage ();
  endif

  p = size(Tau, 1);

  if (df <= (p - 1))
    error('df too small, no finite densities exist')
  endif 

  ## Calculate the logarithm of G_d(df/2), the multivariate gamma function
  g = (p * (p-1) / 4) * log(pi);
  for i = 1:p
    g = g + log(gamma((df + (1 - i))/2));
  endfor

  C = chol(W);

  ## Use formulas for determinant of positive definite matrix for better
  ## efficiency and numerical accuracy
  logdet_W = 2*sum(log(diag(C)));
  logdet_Tau = 2*sum(log(diag(chol(Tau))));

  y = -(df*p)/2 * log(2) + (df/2)*logdet_Tau - g ...
      -((df + p + 1)/2)*logdet_W - trace(Tau*chol2inv(C))/2;

  if ~log_y
    y = exp(y);
  endif
endfunction

## Test results cross-checked against diwish function in R MCMCpack library 
%!assert(iwishpdf(4, 3, 3.1), 0.04226595, 1E-7);
%!assert(iwishpdf([2 -0.3;-0.3 4], [1 0.3;0.3 1], 4), 1.60166e-05, 1E-10);
%!assert(iwishpdf([6 2 5; 2 10 -5; 5 -5 25], [9 5 5; 5 10 -8; 5 -8 22], 5.1), 4.946831e-12, 1E-17);

## Test input validation
%!error iwishpdf ()
%!error iwishpdf (1, 2)
%!error iwishpdf (1, 2, 0)

%!error wishpdf (1, 2)