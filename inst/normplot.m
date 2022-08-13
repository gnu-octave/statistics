## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Based on previous work by Paul Kienzle <pkienzle@users.sf.net>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.
## Author: Paul Kienzle <pkienzle@users.sf.net>
## This program is granted to the public domain.

## -*- texinfo -*-
## @deftypefn {Function File} normplot (@var{x})
## @deftypefnx {Function File} normplot (@var{ax}, @var{x})
## @deftypefnx {Function File} @var{h} = normplot (@dots{})
##
## Produce normal probability plot for each column of @var{x}.
##
## The line joing the 1st and 3rd quantile is drawn on the
## graph.  If the underlying distribution is normal, the
## points will cluster around this line.
##
## Note that this function sets the title, xlabel, ylabel,
## axis, grid, tics and hold properties of the graph.  These 
## need to be cleared before subsequent graphs using 'clf'.
## @end deftypefn

function h = normplot (varargin)
  
  ## Check for valid input arguments
  narginchk (1, 2);
  ## Parse input arguments
  if nargin == 1
    ax = [];
    x = varargin{1};
  else
    ax = varargin{1};
    ## Check that ax is a valid axis handle
    try
      isstruct (get (ax));  
    catch
      error ("normplot: invalid handle %f.", ax);
    end_try_catch
    x = varargin{2};    
  endif
  ## Check that x is a vector or a 2-D matrix
  if isscalar (x) || ndims (x) > 2
    error ("normplot: x must be a vecctor or a 2-D matrix handle.");
  endif
  ## If x is a vector, make it a column vector
  if rows (x) == 1
    x = x(:);
  endif
  ## If ax is empty, create a new axes
  if isempty(ax)
    ax = newplot();
  end
  ## Get number of column vectors in x
  col = size (x, 2);
  ## Process each column and plot data and fit lines
  color = {"blue", "black", "cyan", "green", "magenta", "red", "white", "yellow"};
  hold on;
  for i = 1:col
    xc = x(:,i);
    ## Remove NaNs, get min, max, and range
    xc(isnan(xc)) = [];
    if isempty (xc)
      break;
    endif
    ## Transform data
    row_xc = rows (xc);
    yc = norminv (([1:row_xc]' - 0.5) / row_xc);
    xc = sort(xc);
    ## Find quartiles
    q1x = prctile(xc,25);
    q3x = prctile(xc,75);
    q1y = prctile(yc,25);
    q3y = prctile(yc,75);
    qx = [q1x; q3x];
    qy = [q1y; q3y];
    ## Calculate coordinates and limits for fitting lines
    dx = q3x - q1x;
    dy = q3y - q1y;
    slope = dy ./ dx;
    centerx = (q1x + q3x)/2;
    centery = (q1y + q3y)/2;
    maxx = max (xc);
    minx = min (xc);
    maxy = centery + slope.*(maxx - centerx);
    miny = centery - slope.*(centerx - minx);
    yinter = centery - slope.*(centerx);
    mx = [minx; maxx];
    my = [miny; maxy];
    ## Plot data and corresponding reference lines in the same color,
    ## following the default color order.  Plot reference line first, 
    ## followed by the data, so that data will be on top of reference line.
    h_end(i) = line (ax, mx, my, "LineStyle", "-.", "Marker", "none", ...
                                 "color", color{mod(i,8)});
    h_mid(i) = line (ax, qx, qy, "LineStyle", "-", "Marker", "none", ...
                                 "color", color{mod(i,8)});
    h_dat(i) = line (ax, xc, yc, "LineStyle", "none", "Marker", "+", ...
                                 "color", color{mod(i,8)});
  endfor
  hold off;
  ## Change colors for single column vector
  if i == 1
    set (h_dat, "Color", "b");
    set (h_mid, "Color", "r");
    set (h_end, "Color", "r");
  endif
  ## Bundle handles together if output requested
  if nargout > 0
    h = [h_dat, h_mid, h_end]';
  endif
  ## Plot labels
  title "Normal Probability Plot"
  ylabel "Probability"
  xlabel "Data"
  ## Plot grid
  p = [0.001, 0.003, 0.01, 0.02, 0.05, 0.10, 0.25, 0.5, ...
       0.75, 0.90, 0.95, 0.98, 0.99, 0.997, 0.999];
  label = {"0.001", "0.003", "0.01", "0.02", "0.05", "0.10", "0.25", "0.50", ...
           "0.75", "0.90", "0.95", "0.98", "0.99", "0.997", "0.999"};
  tick  = norminv(p,0,1);
  set (ax, "ytick", tick, "yticklabel", label);
  ## Set view range with a bit of space around data
  range = nanmax (x(:)) - nanmin (x(:));
  if range > 0
    minxaxis = nanmin (x(:)) - 0.025 * range;
    maxxaxis = nanmax (x(:)) + 0.025 * range;
  else
    minxaxis = nanmin (x(:)) - 1;
    maxxaxis = nanmax (x(:)) + 1;
  end
  minyaxis = norminv (0.25 ./ row_xc, 0, 1);
  maxyaxis = norminv ((row_xc - 0.25) ./ row_xc, 0, 1);
  set (ax, "ylim", [minyaxis, maxyaxis], "xlim", [minxaxis, maxxaxis]);
  grid (ax, "on");
  box (ax, "off");
endfunction

%!demo
%! h = normplot([1:20]);
%!demo
%! h = normplot([1:20;5:2:44]');
%!demo
%! ax = newplot();
%! h = normplot(ax, [1:20]);
%! ax = gca;
%! h = normplot(ax, [-10:10]);
%! set (ax, "xlim", [-11, 21]);

%!test
%!error normplot ();
%!error normplot (23);
%!error normplot (23, [1:20]);
%!error normplot (ones(3,4,5));
