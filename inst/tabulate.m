## Copyright (C) 2003 Alberto Terruzzi
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

## usage:  table = tabulate (data,edges)
##
## Frequency table.
##
## table = tabulate (data,edges), for vector data, counts the number of
## values in data that fall between the elements in the edges vector
## (which must contain monotonically non-decreasing values). table is a
## matrix.
## The first column of table is the number of bin, the second
## is the number of instances in each class (absolute frequency). The
## third column contains the percentage of each value (relative
## frequency) and the fourth column contains the cumulative frequency.
## If edges is missed the width of each class is unitary, if edges is a
## scalar then represent the number of classes, or you can define the
## width of each bin.
## table(k,2) will count the value data(i) if edges(k) <= data(i) <
## edges(k+1).  The  last bin will count the value of data(i) if
## edges(k) <= data(i) <=  edges(k+1).  
## Values outside the values in edges are not counted.  Use -inf and inf
## in edges to include all values. 
## Tabulate with no output arguments returns a formatted table in the
## command window. 
##
## Example
##
##   sphere_radius = [1:0.05:2.5];
##   tabulate (sphere_radius)
##
## Tabulate returns 2 bins, the first contains the sphere with radius
## between 1 and 2 mm excluded, and the second one contains the sphere with
## radius between 2 and 3 mm.
##
##   tabulate (sphere_radius,10)
##
## Tabulate returns ten bins.
##
##   tabulate (sphere_radius,[1 1.5 2 2.5])
##
## Tabulate returns three bins, the first contains the sphere with radius
## between 1 and 1.5 mm excluded, the second one contains the sphere with
## radius between 1.5 and 2 mm excluded, and the third contains the sphere with
## radius between 2 and 2.5 mm. 
##
##   bar (table(:,1),table(:,2))
##
## draw histogram.
##
## See also bar and pareto.

## Author: Alberto Terruzzi <t-albert@libero.it>
## Version: 1.0
## Created: 13 February 2003

function table = tabulate (varargin)

if nargin < 1 || nargin > 2
   usage("table = tabulate (data,edges)")
endif

data = varargin{1};
if isvector (data) != 1
  error ("data must be a vector.");
endif
n = length(data);
m = min(data);
M = max(data);

if nargin == 1 edges = 1:1:max(data)+1;
else edges = varargin{2};
end 

if isscalar(edges)
  h=(M-m)/edges;
  edges = [m:h:M];
end

# number of classes
bins=length(edges)-1;
# initialize freqency table
freqtable = zeros(bins,4);

for k=1:1:bins;
  if k != bins
    freqtable(k,2)=length(find (data >= edges(k) & data < edges(k+1)));
  else 
    freqtable(k,2)=length(find (data >= edges(k) & data <= edges(k+1)));
  end
  if k == 1 freqtable (k,4) = freqtable(k,2);
  else freqtable(k,4) = freqtable(k-1,4) + freqtable(k,2); 
  end
end

freqtable(:,1) = edges(1:end-1)(:);
freqtable(:,3) = 100*freqtable(:,2)/n;

if nargout == 0
  disp("     bin     Fa       Fr%        Fc");
  printf("%8g  %5d    %6.2f%%    %5d\n",freqtable');
else table = freqtable;
end

endfunction

