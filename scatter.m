function scatter (varargin)

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
  ##
  ## Author: Alberto Terruzzi <t-albert@libero.it>
  ## Version: 1.0
  ## Created: 19 November 2002
  ##
  ## ToDo: filled markers, RGB support.

error(nargchk(2,6,nargin))
filled = 0; 
marker = "";
c = 1; ## default color
s = []; ## circle surface

Nin=nargin;
while Nin > 0 & isstr(varargin{Nin})
  if strcmp(varargin{Nin},"filled")
    filled = 1;
  else
    marker=varargin{Nin};
    if filled~=1
      filled = 0; 
    endif
  endif
  Nin = Nin-1;
endwhile

if isempty(marker)
  marker = "+";
endif

switch Nin
case 2  ## scatter(x,y)
  x = varargin{1};
  y = varargin{2};
  if nargin < 3 
    marker = "o";
  endif
case 3  ## scatter(x,y,s)
  x = varargin{1};
  y = varargin{2};
  s = varargin{3};
case 4  ## scatter(x,y,s,c)
  x = varargin{1};
  y = varargin{2};
  s = varargin{3};
  c = varargin{4};
otherwise
   error("Wrong number of input arguments.");
endswitch

## check X and Y vectors
if length(x) ~= length(y) | ...
   length(x) ~= prod(size(x)) | length(y) ~= prod(size(y))
  error("X and Y must be vectors of the same length.");
endif

## Map colors in color vector if necessary.
if prod(size(c))==1
  color = repmat(c,length(x),1); 
elseif length(c)==prod(size(c)) & length(c)==length(x), ## C is a vector 
  color = c;
else
  error("C must be a single color or a vector the same length as X.");
endif

## Scalar expand the marker size if necessary.
if Nin > 2
  if length(s)==1,
    s = repmat(s,length(x),1); 
  elseif length(s)~=prod(size(s)) | length(s)~=length(x)
    error("S must be a scalar or a vector the same length as X.");
  endif
endif

## Now draw the plot, one patch per point.
clearplot;
hold on;
				# delete the following line if you use
				# patch
rmax=sqrt(max(s)/pi);
if size(rmax)==[0 0]
  rmax=0;
endif
ex=0.1*(max(x)-min(x))+rmax; ey=0.1*(max(y)-min(y))+rmax;
axis ([min(x)-ex max(x)+ex, min(y)-ey, max(y)+ey]);
for i=1:length(x)
  ## patch(x(i),y(i),marker);
				# patch is too slow!! 
				# delete the next line if you use patch
  plot(x(i), y(i), sprintf(";;%d%s",color(i),marker));
endfor
if Nin > 2
  axis ("equal");
  for i=1:length(x)
    ## t SIN COS is used to draw circles around point
    t=[0:0.01:2*pi]; SIN=sin(t); COS=cos(t);
    ## plot circles 
    plot(x(i)+sqrt(s(i)/pi)*SIN,y(i)+sqrt(s(i)/pi)*COS,sprintf(";;%d",color(i)));
  endfor
endif
## fill the circles
if filled==1
  disp("'fill' isn't supported.");
endif
hold off;
endfunction
