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
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{A} =} zscore (@var{X})
## Compute the @math{z}-score.
##
## Compute the @math{z}-score of each element of @var{X} relative to the data in
## the columns of @var{X}.  The @math{z}-score for a single data point @math{x_i}
## is:
##
## @example
## (x_i - mean(x))/std(x)
## @end example
## @end deftypefn

function [A, mu, sigma] = zscore(X,varargin)
  if (nargin != 1 && nargin != 2)
    print_usage;
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
  
  if (nargout > 1)
    mu = mean (X, dim);
  endif
  if (nargout > 2)
    sigma = std (X, 0, dim);
  endif
endfunction
