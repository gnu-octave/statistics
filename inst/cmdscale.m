## Copyright (C) 2014 Nir Krakauer
##
## This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} @var{Y} = cmdscale (@var{D})
## @deftypefnx{Function File} [@var{Y}, @var{e}] = cmdscale (@var{D})
## Classical multidimensional scaling of a matrix. Also known as principal coordinates analysis.
##
## Given an @var{n} by @var{n} Euclidean distance matrix @var{D}, find @var{n} points in @var{p} dimensional space which have this distance matrix. The coordinates of the points @var{Y} are returned.
## 
## @var{D} should be a full distance matrix (hollow, symmetric, entries obeying the triangle inequality), or can be a vector of length @code{n(n-1)/2} containing the upper triangular elements of the distance matrix (such as that returned by the pdist function). If @var{D} is not a valid distance matrix, points @var{Y} will be returned whose distance matrix approximates @var{D}.
##
## The returned @var{Y} is an @var{n} by @var{p} matrix showing possible coordinates of the points in @var{p} dimensional space (@code{p < n}). Of course, any translation, rotation, or reflection of these would also have the same distance matrix.
## 
## Can also return the eigenvalues @var{e} of @code{(D(1, :) .^ 2 +  D(:, 1) .^ 2 - D .^ 2) / 2}, where the number of positive eigenvalues determines @var{p}.
##
## Reference: Rudolf Mathar (1985), The best Euclidian fit to a given distance matrix in prescribed dimensions, Linear Algebra and its Applications, 67: 1-6, doi: 10.1016/0024-3795(85)90181-8
## 
## @seealso{pdist, squareform}
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>
## Description: Classical multidimensional scaling

function [Y, e] = cmdscale (D)
  
  if isvector (D)
    D = squareform (D);
  endif
  
  warning ("off", "Octave:broadcast","local");
  M = (D(1, :) .^ 2 +  D(:, 1) .^ 2 - D .^ 2) / 2;

  [v e] = eig(M);
  
  e = diag(e);
  pe = (e > 0); #positive eigenvalues
  
  Y = v(:, pe) * diag(sqrt(e(pe)));
 
  
endfunction


%!shared m, n, X, D
%! m = 4; n = 3; X = rand(m, n); D = pdist(X);
%!assert(pdist(cmdscale(D)), D, m*n*eps)
%!assert(pdist(cmdscale(squareform(D))), D, m*n*eps)

