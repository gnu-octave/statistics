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

## usage:  pareto (data, labels, vertical)
##
## Draw Pareto chart, also called ABC chart.
##
## Pareto chart is a bar graph used to arrange information in such a way that
## priorities for process improvement can be established. It organizes and 
## displays information to show the relative importance of data.
## Pareto diagrams are named after Vilfredo Pareto, an Italian sociologist and
## economist, who invented this method of information presentation in 1897.
## The chart is similar to the histogram or bar chart, except that the bars are
## arranged in decreasing order from left to right along the abscissa.
## The fundamental idea (Pareto principle) behind the use of Pareto diagrams is
## that 80% of an effect is due to 20% of the causes, so for quality improvement
## the first few (as presented on the diagram) contributing causes to a problem
## usually account for the majority of the result. Thus, targeting these "major
## causes" for elimination results in the most cost-effective improvement scheme.
##
## pareto (data) produces a Pareto chart where the values in the vector data are
## drawn as bars in descending order and labels each bar with its index.
## pareto (data, vertical) does the same but if vertical is 0 the bar will be drawn
## orizzontaly.
##
## pareto (data, labels, vertical) or pareto (data, vertical, labels) 
## labels each element of data with the values from labels.
##
## Example
## Suppose that, we want establish which  products makes 80 % of turnover.
##
##   CODES = ["AB4";"BD7";"CF8";"CC5";"AD11";"BB5";"BB3";"AD8";"DF3";"DE7"];
##   Value = [2.35 7.9 2.45 1.1 0.15 13.45 5.4 2.05 0.85  1.65]';
##   SoldUnits = [54723 41114 16939 1576091 168000 687197 120222 168195
##   1084118 55576]';
##
##   pareto (Value.*SoldUnits, CODES, 0);
## 
## Pareto returns three classes: the first contains those codes that
## makes 80 % of turnover.
##
## See also tabulate.

## Author: Alberto Terruzzi <t-albert@libero.it>
## Version: 2.0 (Dicember 2003)
## Created: 26 October 2003

function pareto (varargin)

if nargin < 1 || nargin > 3
   usage("pareto (data,labels,vertical)");
endif

data = varargin{1};
if (min(size(data))~=1)
   error("data must be a vector.");
endif
vertical = 1;

if nargin == 2
   if isscalar(varargin {2}) 
     vertical = varargin {2};
   else
     labels = varargin {2};
   endif
endif

if nargin == 3
   if isscalar(varargin {2}) 
     vertical = varargin {2};
     labels = varargin {3};
   else 
     vertical = varargin {3};
     labels = varargin {2};
   endif
endif

data = data(:);
if exist ("labels") && (length(labels) ~= length(data))
  error("size of labels and data must be the same.");
endif

## sort data in decreasing order
[y,j]=sort(data);
y = flipud(y);
j = flipud(j);

## sort labels
if exist ("labels")
  oldlabels = labels;
  ## use numbers as labels
  if isstr(oldlabels) ~= 1
    labels = [];
    for i=1:1:length(oldlabels)
      labels = [labels; num2str(oldlabels(i))];
    endfor
  endif

  oldlabels = labels;
  for i=1:1:length(oldlabels)
    labels(i,:) = oldlabels(j(i),:);
  endfor
else ## use positions as labels
  labels = [];
  for i=1:1:length(j)
    labels = [labels; num2str(j(i))];
  endfor
endif

## compute the cumulative frequency
cum = cumsum (y);
cum = cum/cum(length(cum))*100;

## seek 80-95 boundary
for i=1:1:length(cum)
  if cum(i) == 80 A = i;
  elseif cum(i) < 80 && cum(i+1) >80
    if abs(cum(i) - 80) < abs(cum(i+1) - 80)
      A = i;
    else A = i+1;
    endif
  endif

  if cum(i) == 95 B = i;
  elseif cum(i) < 95 && cum(i+1) >95
    if abs(cum(i) - 95) < abs(cum(i+1) - 95)
      B = i;
    else B = i+1;
    endif
  endif
endfor

[xb, yb] =  bar (y);
dataA = [xb, yb](1:1:3*A+1,:);
dataAB = [xb, yb](3*A+1:1:3*B+1,:);
dataB = [xb, yb](3*B+1:1:length(yb),:);
datacum = [[1:1:length(cum)]', cum];

## draw pareto chart (bar graph plus Lorentz function)
if vertical ~= 0
  if exist("labels")
    tics ("x", 1:1:length(y), labels);
    tics ("y");
  endif
###
  gplot dataA axes x1y1 with lines 1, dataAB axes x1y1 with lines 7, dataB axes x1y1 with lines 3, datacum axes x1y2 with linespoints 13
  gset y2tics; gset nox2tics; replot;
  legend ("hide");
### plot (xb , yb, "b-;;", [1:1:length(cum)], cum*max(y)/100, "c-*;;");
else
  if exist("labels")
    tics ("y", 1:1:length(y), labels);
    tics ("x");
  endif
###
  gplot fliplr(dataA) axes x1y1 with lines 1, fliplr(dataAB) axes x1y1 with lines 7, fliplr(dataB) axes x1y1 with lines 3, fliplr(datacum) axes x2y1 with linespoints 13
  gset x2tics; gset noy2tics; replot;
  legend ("hide");
###  plot (yb , xb, "b-;;", cum*max(y)/100, [1:1:length(cum)], "r-*;;");
endif

endfunction

