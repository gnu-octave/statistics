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

## SCATTER Scatter plot.
##   SCATTER(X,Y,S,C) displays colored circles at the locations specified
##   by the vectors X and Y (which must be the same size).  The area of
##   each marker is determined by the values in the vector S (surface area)
##   and the colors of each marker are based on the values in C. S can be a
##   scalar, in which case all the markers are drawn the same size, or a
##   vector the same length as X and Y.
##   When C is a vector the same length as X and Y, the values in C
##   are used as the color of each circle.  
##   RGB color not supported
##
##   usage:  scatter (...)
##   scatter(X,Y) draws the markers in the default color.
##   scatter(X,Y,S) draws the markers with size S.
##   scatter(X,Y,S,C) draws the markers with size S and color C.
##   scatter(...,M) uses the marker M instead of default one.
##   scatter(...,"filled") fills the markers, this feature isn't supported.
##
##   Use PLOT for single color, single marker size scatter plots.
##
##   Example: A=[0:0.5:2*pi];
##            B=[0:0.2:2.4]; B(1,4)=5; B(1,8)=-5;
##            C=[0:1:6, 0:1:5];
##            scatter (A,B,abs(B*0.2),C,"@","filled");

## Author: Alberto Terruzzi <t-albert@libero.it>
## Version: 1.0
## Created: 19 November 2002
##
## ToDo: filled markers, RGB support.

## 2003-03-15 Paul Kienzle
## * cleanup, optimization

function scatter (x,y,varargin)

error(nargchk(2,6,nargin))
filled = 0; 
marker = "+";
c = 1; ## default color
s = 0; ## circle surface

seen=0;
for i=1:length(varargin)
  v = varargin{i};
  if !ischar(v)
    if seen==0, s = v; 
    elseif seen==1, c = v;
    else error("too many numeric arguments");
    endif
    seen++;
  elseif strcmp(v,"filled")
    filled = 1;
  else
    marker = v;
  endif
endfor

## check X and Y vectors
x=x(:); y=y(:); c=c(:); s = s(:);
if length(x) ~= length(y)
  error("X and Y must be vectors of the same length.");
endif

## Map colors in color vector if necessary.
if isempty(c)
  c = ones(size(x));
elseif length(c)==1
  c = c*ones(size(x));
elseif length(c) ~= length(x)
  error("C must be a single color or a vector the same length as X.");
endif

## Scalar expand the marker size if necessary.
if isempty(s)
  s = zeros(size(x));
elseif length(s)==1
  s = s*ones(size(x));
elseif length(s) ~= length(x)
  error("S must be a scalar or a vector the same length as X.");
endif

## Now draw the plot, one patch per point.
clearplot;
hold on;

## set axis limits
rmax=sqrt(max(s)/pi);
ex=0.1*(max(x)-min(x))+rmax; ey=0.1*(max(y)-min(y))+rmax;
axis ([min(x)-ex max(x)+ex, min(y)-ey, max(y)+ey]);

## draw the markers

for color = unique(c)'
  i =  (color == c);
  ## patch(x(i),y(i),marker);
				# patch is too slow!! 
				# delete the next line if you use patch
  plot(x(i), y(i), sprintf(";;%d%s",color,marker));
endfor

## draw the circles
markers = find(s);
if ~isempty(markers)
  axis ("equal");
  for i=markers'
    ## t SIN COS is used to draw circles around point
    t=[0:0.01:2*pi]; SIN=sin(t); COS=cos(t);
    ## plot circles 
    plot(x(i)+sqrt(s(i)/pi)*SIN,y(i)+sqrt(s(i)/pi)*COS,sprintf(";;%d",c(i)));
  endfor
endif

## fill the circles
if filled==1
  disp("'fill' isn't supported.");
endif
hold off;
endfunction
