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
function a = mad(X)
  if nargin != 1
    usage("a = mad (X)");
  elseif all (size (X) > 1)
    a = mean (abs (X - ones(size(X,1),1) * mean(X)));
  else
    a = mean (abs (X - mean(X)));
  endif
endfunction
