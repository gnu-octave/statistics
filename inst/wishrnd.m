## Copyright (C) 2013 Nir Krakauer <nkrakauer@ccny.cuny.edu>
##
## This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with Octave; see the file COPYING.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {} [@var{W}[, @var{D}]] = wishrnd (@var{Sigma}, @var{df}[, @var{D}][, @var{n}=1])
## Return a random matrix sampled from the Wishart distribution with given parameters
##
## Inputs: the @var{p} x @var{p} positive definite matrix @var{Sigma} and scalar degrees of freedom parameter @var{df} (and optionally the Cholesky factor @var{D} of @var{Sigma}).
## @var{df} can be non-integer as long as @var{df} > @var{p}
##
## Output: a random @var{p} x @var{p}  matrix @var{W} from the Wishart(@var{Sigma}, @var{df}) distribution. If @var{n} > 1, then @var{W} is @var{p} x @var{p} x @var{n} and holds @var{n} such random matrices. (Optionally, the Cholesky factor @var{D} of @var{Sigma} is also returned.)
##
## Averaged across many samples, the mean of @var{W} should approach @var{df}*@var{Sigma}, and the variance of each element @var{W}_ij should approach @var{df}*(@var{Sigma}_ij^2 + @var{Sigma}_ii*@var{Sigma}_jj)
##
## Reference: Yu-Cheng Ku and Peter Bloomfield (2010), Generating Random Wishart Matrices with Fractional Degrees of Freedom in OX, http://www.gwu.edu/~forcpgm/YuChengKu-030510final-WishartYu-ChengKu.pdf
## 
## @seealso{iwishrnd, wishpdf}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: Compute the probability density function of the Wishart distribution

function [W, D] = wishrnd(Sigma, df, D, n=1)

if (nargin < 3)
  print_usage ();
endif

if nargin < 3 || isempty(D)
  try
    D = chol(Sigma);
  catch
    error('Cholesky decomposition failed; Sigma probably not positive definite')
  end_try_catch
endif

p = size(D, 1);

if df < p
  df = floor(df); #distribution not defined for small noninteger df
  df_isint = 1;
else 
#check for integer degrees of freedom
 df_isint = (df == floor(df));
endif

if ~df_isint
  [ii, jj] = ind2sub([p, p], 1:(p*p));
endif

if n > 1
  W = nan(p, p, n);
endif

for i = 1:n
  if df_isint
    Z = randn(df, p) * D;
    W(:, :, i) = Z'*Z;
  else
    Z = diag(sqrt(chi2rnd(df - (0:(p-1))))); #fill diagonal
    #note: chi2rnd(x) is equivalent to 2*randg(x/2), but the latter seems to offer no performance advantage
    Z(ii > jj) = randn(p*(p-1)/2, 1); #fill lower triangle with normally distributed variates
    Z = D * Z;
    W(:, :, i) = Z*Z';
  endif

  
endfor

endfunction



%!assert(size (wishrnd (1,2,1)), [1, 1]);
%!assert(size (wishrnd ([],2,1)), [1, 1]);
%!assert(size (wishrnd ([3 1; 1 3], 2.00001, [], 1)), [2, 2]);
%!assert(size (wishrnd (eye(2), 2, [], 3)), [2, 2, 3]);

%% Test input validation
%!error wishrnd ()
%!error wishrnd (1)
%!error wishrnd ([-3 1; 1 3],1)
%!error wishrnd ([1; 1],1)
