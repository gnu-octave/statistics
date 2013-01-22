## Copyright (C) 2013 Fernando Damian Nieuwveldt <fdnieuwveldt@gmail.com>
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{residuals},@var{reconstructed}]}=pcares(@var{X}, @var{NDIM})
## @itemize @bullet
## @item
## @var{X} : N x P Matrix with N observations and P variables, the variables will be mean centered
## @item
## @var{ndim} : Is a scalar indicating the number of principal components to use and should be <= P
## @end itemize
##
## @subheading References
##
## @enumerate
## @item
## Jolliffe, I. T., Principal Component Analysis, 2nd Edition, Springer, 2002
## 
## @end enumerate
## @end deftypefn

## Author: Fernando Damian Nieuwveldt <fdnieuwveldt@gmail.com>
## Description: Residuals from Principal Components Analysis

function [residuals,reconstructed] = pcares(X,NDIM)

  if (nargin ~= 2)
    error('pcares takes two inputs: The data Matrix X and number of principal components NDIM')
  endif    

  # Mean center data
  Xcentered = bsxfun(@minus,X,mean(X));

  # Apply svd to get the principal component coefficients
  [U,S,V] = svd(Xcentered);

  # Use only the first ndim PCA components
  v = V(:,1:NDIM);

  if (nargout == 2)
    # Calculate the residuals
    residuals = Xcentered - Xcentered * (v*v'); 
				  
    # Reconstructed data using ndim PCA components
    reconstructed = X - residuals;
  else
     # Calculate the residuals
     residuals = Xcentered - Xcentered * (v*v'); 
  endif    
endfunction
%!demo
%! X = [ 7    26     6    60;
%!       1    29    15    52;
%!      11    56     8    20;
%!      11    31     8    47;
%!       7    52     6    33;
%!      11    55     9    22;
%!       3    71    17     6;
%!       1    31    22    44;
%!       2    54    18    22;
%!      21    47     4    26;
%!       1    40    23    34;
%!      11    66     9    12;
%!      10    68     8    12 
%!     ];
%! # As we increase the number of principal components, the norm 
%! # of the residuals matrix will decrease
%! r1 = pcares(X,1);
%! n1 = norm(r1)
%! r2 = pcares(X,2);
%! n2 = norm(r2)
%! r3 = pcares(X,3);
%! n3 = norm(r3)
%! r4 = pcares(X,4);
%! n4 = norm(r4)
    
