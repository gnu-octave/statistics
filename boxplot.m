## Copyright (C) 2002 Alberto Terruzzi
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

## usage: s = boxplot (data,notched,symbol,vertical,maxwhisker);
##
## The box plot is a graphical display that simultaneously describes several 
## important features of a data set, such as center, spread, departure from 
## symmetry, and identification of observations that lie unusually far from
## the bulk of the data.
##
## data is a matrix with one column for each dataset, or data is a cell
## vector with one cell for each dataset.
## s is a matrix with one column for each of: 
##   minimum, first quartile, median, third quartile, maximum
##
## Example
##
##   title("Grade 3 heights");
##   tics("x",1:2,["girls";"boys"]);
##   axis([0,3]);
##   boxplot({randn(10,1)*5+140, randn(13,1)*8+135});
##

## Author: Alberto Terruzzi <t-albert@libero.it>
## Version: 1.3
## Created: 6 January 2002

function [s] = boxplot (data,notched,symbol,vertical,maxwhisker)

if nargin < 1 || nargin > 6
   usage("s = boxplot (data,notch,symbol,vertical,maxwhisker)")
endif
if nargin < 5, maxwhisker = 1.5; end
if nargin < 4, vertical = 1; end
if nargin < 3, symbol = '+'; end
if nargin < 2, notch = 0; end

if iscell(data), 
  nc = length(data);
else
  if isvector(data), data = data(:); end
  nc = columns(data);
end

whisker_x = ones(2,1)*[1:nc,1:nc];
whisker_y = zeros(2,2*nc);
outliers_x = [];
outliers_y = [];
s = zeros(5,nc);
box = zeros(1,nc);

for i=1:nc
  if iscell(data)
    col = data{i};
  else
    col = data(:,i);
  endif
  col(is_nan_or_na (col)) = [];
  nd = length(col);
  box(i) = nd;
  if (nd > 1)
    s(:,i) = statistics(col)(1:5);
    IQR = maxwhisker*(s(4,i)-s(2,i));
    whisker_y(:,i) = [min(col(col >= s(2,i)-IQR)); s(2,i)];
    whisker_y(:,nc+i) = [max(col(col <= s(4,i)+IQR)); s(4,i)];
    outliers = col(col < s(2,i)-IQR | col > s(4,i)+IQR);
    outliers_x = [outliers_x; i*ones(size(outliers))];
    outliers_y = [outliers_y; outliers];
  else
    s(:,i) = NaN;
  end
end

chop = find(box <= 1);

# Draw a box around the quartiles, with width proportional to the number of
# items in the box.
box *= 0.4/max(box);
quartile_x = ones(5,1)*[1:nc];
quartile_x(1,:) -= box;
quartile_x(2,:) += box;
quartile_x(3,:) += box;
quartile_x(4,:) -= box;
quartile_x(5,:) -= box;
quartile_y = zeros(5,nc);
quartile_y(1,:) = s(2,:);
quartile_y(2,:) = s(2,:);
quartile_y(3,:) = s(4,:);
quartile_y(4,:) = s(4,:);
quartile_y(5,:) = s(2,:);

# Draw a line through the median
median_x = ones(2,1)*[1:nc];
median_x(1,:) -= box;
median_x(2,:) += box;
median_y = zeros(2,nc);
median_y(1,:) = s(3,:);
median_y(2,:) = s(3,:);

# Chop all missing columns
quartile_x(:,chop) = [];
quartile_y(:,chop) = [];
whisker_x(:,chop+nc) = [];
whisker_y(:,chop+nc) = [];
whisker_x(:,chop) = [];
whisker_y(:,chop) = [];
median_x(:,chop) = [];
median_y(:,chop) = [];

# Add caps to the remaining whiskers
cap_x = whisker_x;
cap_x(1,:) -= 0.05;
cap_x(2,:) += 0.05;
cap_y = [whisker_y(1,:); whisker_y(1,:)];

#quartile_x,quartile_y
#whisker_x,whisker_y
#median_x,median_y
#cap_x,cap_y

if vertical
  plot (quartile_x, quartile_y, "b;;",
	whisker_x, whisker_y, "b;;",
	cap_x, cap_y, "b;;",
	median_x, median_y, "r;;",
	outliers_x, outliers_y, [symbol,"r;;"])
else
  plot (quartile_y, quartile_x, "b;;",
	whisker_y, whisker_x, "b;;",
	cap_x, cap_y, "b;;",
	median_y, median_x, "r;;",
	outliers_y, outliers_x, [symbol,"r;;"]);
endif

endfunction
