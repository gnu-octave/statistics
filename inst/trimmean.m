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

## a = trimmean(x,p)
##    mean of x excluding highest and lowest p% of the data
##
## E.g.,
##    mean([-inf 1:9 inf]) is NaN
##    trimmean([-inf 1:9 inf], 10) is 5
function a = trimmean(x, p, varargin)
  if (nargin != 2 && nargin != 3)
    usage("a = trimmean(x,p, dim)");
  endif
  y = sort(x, varargin{:});
  sz = size(x);
  if nargin < 3
    dim = min(find(sz>1));
    if isempty(dim), dim=1; endif;
  else
    dim = varargin{1};
  endif
  idx = cell (0);
  for i=1:length(sz), idx{i} = 1:sz(i); end;
  trim = round(sz(dim)*p*0.01);
  idx{dim} = 1+trim : sz(dim)-trim;
  a = mean (y (idx{:}), varargin{:});
endfunction
