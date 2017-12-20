## Copyright (C) 2016 - Juan Pablo Carbajal
##
## This progrm is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## Author: Juan Pablo Carbajal <ajuanpi+dev@gmail.com>

## -*- texinfo -*-
## @defun {@var{h} =} violin (@var{x})
## @defunx {@var{h} =} violin (@dots{}, @var{property}, @var{value}, @dots{})
## @defunx {@var{h} =} violin (@var{hax}, @dots{})
## @defunx {@var{h} =} violin (@dots{}, @asis{"horizontal"})
## Produce a Violin plot of the data @var{x}.
##
## The input data @var{x} can be a N-by-m array containg N observations of m variables.
## It can also be a cell with m elements, for the case in which the varibales
## are not uniformly sampled.
##
## The following @var{property} can be set using @var{property}/@var{value} pairs
## (default values in parenthesis).
## The value of the property can be a scalar indicating that it applies
## to all the variables in the data.
## It can also be a cell/array, indicating the property for each variable.
## In this case it should have m columns (as many as variables).
##
## @table @asis
##
## @item Color
## (@asis{"y"}) Indicates the filling color of the violins.
##
## @item Nbins
## (50) Internally, the function calls @command{hist} to compute the histogram of the data.
## This property indicates how many bins to use. See @command{help hist}
## for more details.
##
## @item SmoothFactor
## (4) The fuction performs simple kernel density estimation and automatically
## finds the bandwith of the kernel function that best approximates the histogram
## using optimization (@command{sqp}).
## The result is in general very noisy. To smooth the result the bandwidth is
## multiplied by the value of this property. The higher the value the smoother
## the violings, but values too high might remove features from the data distribution.
##
## @item Bandwidth
## (NA) If this property is given a value other than NA, it sets the bandwith of the
## kernel function. No optimization is peformed and the property @asis{SmoothFactor}
## is ignored.
##
## @item Width
## (0.5) Sets the maximum width of the violins. Violins are centered at integer axis
## values. The distance between two violin middle axis is 1. Setting a value
## higher thna 1 in this property will cause the violins to overlap.
## @end table
##
## If the string @asis{"Horizontal"} is among the input arguments, the violin
## plot is rendered along the x axis with the variables in the y axis.
##
## The returned structure @var{h} has handles to the plot elements, allowing
## customization of the visualization using set/get functions.
##
## Example:
##
## @example
## title ("Grade 3 heights");
## axis ([0,3]);
## set (gca, "xtick", 1:2, "xticklabel", @{"girls"; "boys"@});
## h = violin (@{randn(100,1)*5+140, randn(130,1)*8+135@}, "Nbins", 10);
## set (h.violin, "linewidth", 2)
## @end example
##
## @seealso{boxplot, hist}
## @end defun

function h = violin (ax, varargin)

  old_hold = ishold ();
  # First argument is not an axis
  if (~ishandle (ax) || ~isscalar (ax))
    x  = ax;
    ax = gca ();
  else
    x = varargin{1};
    varargin(1) = [];
  endif

  ######################
  ## Parse parameters ##
  parser = inputParser ();
  parser.CaseSensitive = false;
  parser.FunctionName = 'violin';

  parser.addParamValue ('Nbins', 50);
  parser.addParamValue ('SmoothFactor', 4);
  parser.addParamValue ('Bandwidth', NA);
  parser.addParamValue ('Width', 0.5);
  parser.addParamValue ('Color', "y");
  parser.addSwitch ('Horizontal');

  parser.parse (varargin{:});
  res = parser.Results;

  c        = res.Color;        # Color of violins
  if (ischar (c)) c = c(:); endif
  nb       = res.Nbins;        # Number of bins in histogram
  sf       = res.SmoothFactor; # Smoothing factor for kernel estimation
  r0       = res.Bandwidth;    # User value for KDE bandwth to prevent optimization
  is_horiz = res.Horizontal;   # Whether the plot must be rotated
  width    = res.Width;        # Width of the violins
  clear parser res
  ######################

  ## Make everything a cell for code simplicity
  if (~iscell (x))
    [N Nc] = size (x);
    x      = mat2cell (x, N, ones (1, Nc));
  else
    Nc = numel (x);
  endif

  try
    [nb, c, sf, r0, width] = to_cell (nb, c, sf, r0, width, Nc);
  catch err
    if strcmp (err.identifier, "to_cell:element_idx")
      n = str2num (err.message);
      txt = {"Nbins", "Color", "SmoothFactor", "Bandwidth", "Width"};
      error ("Octave:invaid-input-arg", ...
             ["options should be scalars or call/array with as many values as" ...
              " numbers of variables in the data (wrong size of %s)."], txt{n});
    else
      rethrow (lasterror())
    endif
  end

  ## Build violins
  [px py mx] = cellfun (@(y,n,s,r)build_polygon(y, n, s, r), ...
                          x, nb, sf, r0, "unif", 0);

  Nc    = 1:numel (px);
  Ncc   = mat2cell (Nc, 1, ones (1, Nc(end)));

  # get hold state
  old_hold = ishold ();

  # Draw plain violins
  tmp      = cellfun (@(x,y,n,u, w)patch(ax, (w * x + n)(:), y(:) ,u), ...
                        px, py, Ncc, c, width);
  h.violin = tmp;

  hold on
  # Overlay mean value
  tmp    = cellfun (@(z,y)plot(ax, z, y,'.k', "markersize", 6), Ncc, mx);
  h.mean = tmp;

  # Overlay median
  Mx       = cellfun (@median, x, "unif", 0);
  tmp      = cellfun (@(z,y)plot(ax, z, y, 'ok'), Ncc, Mx);
  h.median = tmp;

  # Overlay 1nd and 3th quartiles
  LUBU = cellfun (@(x,y)abs(quantile(x,[0.25 0.75])-y), x, Mx, "unif", 0);
  tmp  = cellfun (@(x,y,z)errorbar(ax, x, y, z(1),z(2)), Ncc, Mx, LUBU)(:);
  # Flatten errorbar output handles
  tmp2       = allchild (tmp);
  if (~iscell (tmp2))
    tmp2 = mat2cell (tmp2, ones(length (tmp2), 1), 1);
  endif
  tmp        = mat2cell (tmp, ones (length (tmp), 1), 1);
  tmp        = cellfun (@vertcat, tmp, tmp2, "unif", 0);
  h.quartile = cell2mat (tmp);

  hold off

  # Rotate the plot if it is horizontal
  if (is_horiz)
    structfun (@swap_axes, h);
    set (ax, "ytick", Nc);
  else
    set (ax, "xtick", Nc);
  endif

  if (nargout < 1);
    clear h;
  endif

  # restore hold state
  if (old_hold)
    hold on
  endif
endfunction

function k = kde(x,r)
  k  = mean (stdnormal_pdf (x / r)) / r;
  k /= max (k);
endfunction

function [px py mx] = build_polygon (x, nb, sf, r)
  N  = size (x, 1);
  mx = mean (x);
  sx = std (x);
  X  = (x - mx ) / sx;

  [count bin] = hist (X, nb);
  count      /= max (count);

  Y  = X - bin;
  if isna (r)
    r0 = 1.06 * N^(1/5);
    r  = sqp (r0, @(r)sumsq (kde(Y,r) - count), [], [], 1e-3, 1e2);
  else
    sf = 1;
  endif
  sig = sf * r;

  ## Create violin polygon
  # smooth tails: extend to 1.83 sigmas, i.e. ~99% of data.
  xx  = linspace (0, 1.83 * sig, 5);
  bin = [bin(1)-fliplr(xx) bin bin(end)+xx];
  py  = [bin; fliplr(bin)].'  * sx + mx;

  v  = kde (X-bin, sig).';
  px = [v -flipud(v)];

endfunction

function tf = swap_axes (h)
    tmp  = mat2cell (h(:), ones (length (h),1), 1);
%    tmp  = cellfun (@(x)[x; allchild(x)], tmp, "unif", 0);
    tmpy = cellfun(@(x)get(x, "ydata"), tmp, "unif", 0);
    tmpx = cellfun(@(x)get(x, "xdata"), tmp, "unif", 0);
    cellfun (@(h,x,y)set (h, "xdata", y, "ydata", x), tmp, tmpx, tmpy);
    tf = true;
endfunction

function varargout = to_cell (varargin)

    m = varargin{end};
    varargin(end) = [];

    for i = 1:numel(varargin)
      x  = varargin{i};
      if (isscalar (x)) x = repmat (x, m, 1); endif

      if (iscell (x))
        if (numel(x) ~= m) # no dimension equals m
          error ("to_cell:element_idx", "%d\n",i);
        endif
        varargout{i} = x;
        continue
      endif

      sz = size (x);
      d  = find (sz == m);
      if (isempty (d)) # no dimension equals m
        error ("to_cell:element_idx", "%d\n",i);
      elseif (length (d) == 2)
        #both dims are m, choose 1st
      elseif (d == 1) # 2nd dimension is m --> transpose
        x  = x.';
        sz = fliplr (sz);
      endif

      varargout{i} = mat2cell (x, sz(1), ones (m,1));

    endfor

endfunction

%!demo
%! clf
%! x = zeros (9e2, 10);
%! for i=1:10
%!   x(:,i) = (0.1 * randn (3e2, 3) * (randn (3,1) + 1) + ...
%!          2 * randn (1,3))(:);
%! endfor
%! h = violin (x, "color", "c");
%! axis tight
%! set (h.violin, "linewidth", 2);
%! set (gca, "xgrid", "on");
%! xlabel ("Variables")
%! ylabel ("Values")

%!demo
%! clf
%! data = {randn(100,1)*5+140, randn(130,1)*8+135};
%! subplot (1,2,1)
%! title ("Grade 3 heights - vertical");
%! set (gca, "xtick", 1:2, "xticklabel", {"girls"; "boys"});
%! violin (data, "Nbins", 10);
%! axis tight
%!
%! subplot(1,2,2)
%! title ("Grade 3 heights - horizontal");
%! set (gca, "ytick", 1:2, "yticklabel", {"girls"; "boys"});
%! violin (data, "horizontal", "Nbins", 10);
%! axis tight

%!demo
%! clf
%! data = exprnd (0.1, 500,4);
%! violin (data, "nbins", {5,10,50,100});
%! axis ([0 5 0 max(data(:))])

%!demo
%! clf
%! data = exprnd (0.1, 500,4);
%! violin (data, "color", jet(4));
%! axis ([0 5 0 max(data(:))])

%!demo
%! clf
%! data = repmat(exprnd (0.1, 500,1), 1, 4);
%! violin (data, "width", linspace (0.1,0.5,4));
%! axis ([0 5 0 max(data(:))])

%!demo
%! clf
%! data = repmat(exprnd (0.1, 500,1), 1, 4);
%! violin (data, "nbins", [5,10,50,100], "smoothfactor", [4 4 8 10]);
%! axis ([0 5 0 max(data(:))])
