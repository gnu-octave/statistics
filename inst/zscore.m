## Copyright (C) 2001 Paul Kienzle
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

## A = zscore (X)
## compute the z-score of each element of X relative to the data in the
## columns of X.  The z-score for a single data point x_i is:
##    (x_i - mean(x))/std(x)
function A = zscore(X,varargin)
  if (nargin != 1 && nargin != 2)
    usage("zscore(X,dim)");
  endif
  if (nargin == 2)
    dim = varargin{1}
  else
    dim = min(find(size(X)>1));
    if isempty(dim), dim=1; endif;
  endif
  
  sz = ones(1,length(size(X)));
  sz(dim) = size(X,dim);
  A = (X - repmat(mean(X,varargin{:}),sz)) ./ repmat(std(X,varargin{:}),sz);
endfunction
