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

## a = trimmean(x,p)
##    mean of x excluding highest and lowest p% of the data
##
## E.g.,
##    mean([-inf 1:9 inf]) is NaN
##    trimmean([-inf 1:9 inf], 10) is 5
function a = trimmean(x, p)
  if nargin != 2
    usage("a = trimmean(x,p)");
  endif
  y = sort(x);
  if size (y,1) == 1, y = y.'; endif
  trim = round(size(y,1)*p*0.01);
  rng = 1+trim : size(y,1)-trim;
  a = mean ( y (rng, :) );
endfunction
