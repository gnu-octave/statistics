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
## pareto (data, vertical) does the same but if vertical is 0 the bar
## will be drawn horizontally.
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
  n = length(data);
  
  ## by default use position number in original vector as a label
  if ~exist ("labels"), labels = [1:n]; end
  
  ## check if using numbers as labels
  if ~isstr(labels), labels = num2str(labels(:)); endif
  
  ## sort labels along with data
  labels = labels(j,:);
  
  ## compute the cumulative frequency
  cum = cumsum (y);

  ## seek 80% and 95% boundaries
  p80 = .8*cum(end);
  p95 = .95*cum(end);
  [m,Aidx] = min(abs(cum-p80));
  [m,Bidx] = min(abs(cum-p95));

  [xb, yb] =  bar (y);
  A = 1:3*Aidx+1;
  AB = 3*Aidx+1:3*Bidx+1;
  B = 3*Bidx+1:length(yb);

  ## draw pareto chart (bar graph plus Lorentz function)
  text;
  if vertical ~= 0
    tics ("x", 1:n, labels);
    tics ("y");
    plot(xb(A),yb(A),"b-;;", xb(AB), yb(AB),"c-;;", xb(B), yb(B),"m-;;",
	 [1:Aidx], cum(1:Aidx), "b-*;;", 
	 [Aidx:Bidx], cum(Aidx:Bidx), "c-*;;",
	 [Bidx:n], cum(Bidx:end), "m-*;;",
	 [1,n+.5],[p80,p80],"b-;;",
	 [1,n+.5],[p95,p95],"m-;;");
    text(.5,p80," 80%",'HorizontalAlignment','left');
    text(.5,p95," 95%",'HorizontalAlignment','left');
  else
    tics ("y", [1:n], labels);
    tics ("x");
    plot(yb(A),xb(A),"b-;;", yb(AB), xb(AB),"c-;;", yb(B), xb(B),"m-;;",
	 cum(1:Aidx), [1:Aidx], "b-*;;", 
	 cum(Aidx:Bidx), [Aidx:Bidx], "c-*;;",
	 cum(Bidx:end), [Bidx:n], "m-*;;",
	 [p80,p80],[.5,n+.5],"b-;;",
	 [p95,p95],[.5,n+.5],"m-;;");
    text(p80,n,"80%",'HorizontalAlignment','right');
    text(p95,n,"95%",'HorizontalAlignment','right');
  endif
  
endfunction

%!demo
%! % Suppose that we want establish which products makes 80 % of turnover.
%! CODES = ["AB4";"BD7";"CF8";"CC5";"AD11";"BB5";"BB3";"AD8";"DF3";"DE7"];
%! Value = [2.35 7.9 2.45 1.1 0.15 13.45 5.4 2.05 0.85  1.65]';
%! SoldUnits = [54723 41114 16939 1576091 168000 687197 120222 168195, ...
%!              1084118 55576]';
%! pareto (Value.*SoldUnits, CODES, 0);
