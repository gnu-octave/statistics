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
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## a = mad(X)
##    mean absolute deviation of X
function a = mad(X,dim)
  if nargin == 1
    dim = min(find(size(X)>1));
    if isempty(dim), dim=1; endif;
  endif
  if (nargin != 1 || nargin != 2)
    usage("a = mad (X,dim)");
  elseif (prod(size(X)) != size(X,dim))
    sz = ones(1,length(size(X));
    sz(dim) = size(X,dim);
    a = nanmean (abs (X - repmat (nanmean (X, dim), sz)), dim);
  elseif all (size (X) > 1)
    a = nanmean (abs (X - ones(size(X,1),1) * nanmean(X, dim)), dim);
  else
    a = nanmean (abs (X - nanmean(X, dim)), dim);
  endif
endfunction
