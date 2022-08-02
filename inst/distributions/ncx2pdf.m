## Copyright (C) 2018 gold holk
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} ncx2pdf (@var{X}, @var{N}, @var{LAMBDA})
## @deftypefnx {Function File} ncx2pdf (@dots{}, @var{TERM})
##   compute the non-central chi square probalitity density function
##   at @var{X} , degree of freedom @var{N} ,
##   and non-centrality parameter @var{LAMBDA} .
##
##   @var{TERM} is the term number of series, default is 32.
##   
## @end deftypefn

## Author: gold holk <gholk@dt13>
## Created: 2018-10-25

function f = ncx2pdf(x, n, lambda, term = 32)
  f = exp(-lambda/2) * arrayfun(@(x) sum_expression([0:term],x,n,lambda), x);
end

function t = sum_expression(j,v,n,l)
  # j is vector, v is scalar.
  numerator = (l/2).^j .* v.^(n/2+j-1) * exp(-v/2);
  denominator = factorial(j) .* 2.^(n/2+j) .* gamma(n/2+j);
  t = sum(numerator ./ denominator);
end


%!assert (ncx2pdf (3, 4, 0), chi2pdf(3, 4), eps)
%!assert (ncx2pdf (5, 3, 1), 0.091858459565020, 1E-15) #compared with Matlab's values
%!assert (ncx2pdf (4, 5, 2), 0.109411958414115, 1E-15)

