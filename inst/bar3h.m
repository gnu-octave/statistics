## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
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

## -*- texinfo -*-
## @deftypefn  {statistics} {} bar3h (@var{y})
## @deftypefnx {statistics} {} bar3h (@var{z}, @var{y})
## @deftypefnx {statistics} {} bar3h (@dots{}, @var{width})
## @deftypefnx {statistics} {} bar3h (@dots{}, @var{style})
## @deftypefnx {statistics} {} bar3h (@dots{}, @var{color})
## @deftypefnx {statistics} {} bar3h (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {} bar3h (@var{ax}, @dots{})
## @deftypefnx {statistics} {@var{p} =} bar3h (@dots{})
##
## Plot a horizontal 3D bar graph.
##
## @code{bar3h (@var{y})} plots 3D bar graph for the elements of @var{y}.  Each
## bar corresponds to an element in @var{y}, which can be a scalar, vector, or
## 2D matrix.  By default, each column in @var{y} is considered as a series and
## it is handled as a distinct series of bars.  When @var{y} is a vector, unlike
## MATLAB, which plots it as a single series of bars, Octave distinguishes
## between a row and column vector of @var{y}.  Hence, when @var{y} is column
## vector, it is plotted as a single series of bars (same color), whereas when
## @var{y} is row vector, each bar is plotted as a different group (different
## colors).  For an @math{MxN} matrix, the function plots the bars corresponding
## to each row on the @qcode{z-axis} ranging from @math{1} to @math{M} and each
## column on the @qcode{x-axis} ranging from @math{1} to @math{N}.
##
## @code{bar3h (@var{z}, @var{y})} plots a 3D bar graph of the elements in
## @var{y} at the @qcode{z-values} specified in @var{z}.  It should be noted
## that @var{z} only affects the tick names along the @qcode{z-axis} rather the
## actual values.  If you want to specify non-numerical values for @var{z}, you
## can specify it with the paired @var{name}/@var{value} syntax shown below.
##
## @code{bar3h (@dots{}, @var{width})} sets the width of the bars along the
## @qcode{x-} and @qcode{z-axes} and controls the separation of bars among each
## other. @var{width} can take any value in the range @math{(0,1]}.  By default,
## @var{width} is 0.8 and the bars have a small separation. If width is 1, the
## bars touch one another.  Alternatively, you can define @var{width} as a two-
## element vector using the paired @var{name}/@var{value} syntax shown below, in
## which case you can control the bar separation along each axis independently.
##
## @code{bar3h (@dots{}, @var{style})} specifies the style of the bars, where
## @var{style} can be @qcode{'detached'}, @qcode{'grouped'}, or
## @qcode{'stacked'}. The default style is @qcode {'detached'}.
##
## @code{bar3h (@dots{}, @var{color})} displays all bars using the color
## specified by color.  For example, use @qcode{'red'} or @qcode{'red'} to
## specify all red bars.  When you want to specify colors for several groups,
## @var{color} can be a cellstr vector with each element specifying the color of
## each group.  @var{color} can also be specified as a numerical @math{Mx3}
## matrix, where each row corresponds to a RGB value with its elements in the
## range @math{[0,1]}.  If only one color is specified, then it applies to all
## bars.  If the number of colors equals the number of groups, then each color
## is applied to each group.  If the number of colors equals the number of
## elements in @var{y}, then each individual bar is assigned the particular
## color.  You can also define @var{color} using the paired
## @var{name}/@var{value} syntax shown below.
##
## @code{bar3h (@dots{}, @var{name}, @var{value})} specifies one or more of the
## following name/value pairs:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @headitem @tab Name @tab Value
## @item @tab @qcode{"width"} @tab A two-element vector specifying the width if
## of the bars along the @qcode{x-} and @qcode{z-axes}, respectively.  Each
## element must be in the range @math{(0,1]}.
##
## @item @tab @qcode{"color"} @tab A character or a cellstr vector, or a
## numerical @math{Mx3} matrix following the same conventions as the @var{color}
## input argument.
##
## @item @tab @qcode{"xlabel"} @tab A cellstr vector specifying the group names
## along the @qcode{x-axis}.
##
## @item @tab @qcode{"zlabel"} @tab A cellstr vector specifying the names of the
## bars in the same series along the @qcode{z-axis}.
## @end multitable
##
## @code{bar3h (@var{ax}, @dots{})} can also take an axes handle @var{ax} as a
## first argument in which case it plots into the axes specified by @var{ax}
## instead of into the current axes specified by @code{gca ()}.  The optional
## argument @var{ax} can precede any of the input argument combinations in the
## previous syntaxes.
##
## @code{@var{p} = bar3h (@dots{})} returns a patch handle @var{p}, which can be
## used to set properties of the bars after displaying the 3D bar graph.
##
## @seealso{boxplot, hist3}
## @end deftypefn

function [varargout] = bar3h (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ## Check if first input is an axes handle
  if (isaxes (varargin{1}))
    ax = varargin{1};
    varargin(1) = [];
    new_axes = false;
  else
    new_axes = true;
  endif

  ## Parse input argument Z
  if (numel (varargin) < 1)
    print_usage ();
  endif
  y = varargin{1};
  varargin(1) = [];
  if (! isnumeric (y))
    error ("bar3h: Z must be numeric.");
  endif

  ## Add defaults
  z = [];
  width = 0.8;
  depth = 0.8;
  style = "detached";
  color = [];
  xlabel = [];
  zlabel = [];

  ## Valid Colors for input validation
  vc = {'red', 'r', 'green', 'g', 'blue', 'b', 'cyan', 'c', ...
        'magenta', 'm', 'yellow', 'y', 'black', 'k', 'white', 'w'};

  ## Parse extra input arguments
  while (numel (varargin) > 0)
    tmp = varargin{1};
    if (isnumeric (tmp) && isempty (z))
      if (isvector (tmp) && isvector (y) && numel (tmp) == numel (y))
        z = y;
        y = tmp;
      elseif (all (size (tmp) > 1) && size (tmp, 1) == numel (y) && isvector (y))
        z = y;
        y = tmp;
      elseif (isscalar (tmp))
        z = NaN;
        if (tmp > 0 && tmp <= 1)
          width = tmp;
          depth = tmp;
        else
          error ("bar3h: WIDTH must be a scalar in the range (0,1].");
        endif
      elseif (size (tmp, 2) == 3 && all (tmp(:) >= 0) && all (tmp(:) <= 1))
        z = NaN;
        color = tmp;
      else
        error ("bar3h: inconsistent size in Y and Z input arguments.");
      endif
      varargin(1) = [];
    elseif (isnumeric (tmp) && isscalar (tmp))
      if (tmp > 0 && tmp <= 1)
        width = tmp;
        depth = tmp;
      else
        error ("bar3h: WIDTH must be a scalar in the range (0,1].");
      endif
      varargin(1) = [];
    elseif (isnumeric (tmp))
      if (size (tmp, 2) == 3 && all (tmp(:) >= 0) && all (tmp(:) <= 1))
        color = tmp;
      else
        error (["bar3h: numeric COLOR must be a 1x3 vector of an Nx3 matrix", ...
                " where each value is between 0 and 1 inclusive."]);
      endif
      color = tmp;
      varargin(1) = [];
    elseif (ischar (tmp))
      if (any (strcmpi (tmp, {"detached", "grouped", "stacked"})))
        style = tmp;
        varargin(1) = [];
      elseif (any (strcmpi (tmp, vc)))
        color = tmp;
        varargin(1) = [];
      elseif (strcmpi (tmp, "width"))
        if (numel (varargin) < 2)
          error ("bar3h: missing value for optional argument 'width'.");
        endif
        w = varargin{2};
        if (isscalar (w) && isnumeric (w) && isfinite (w) && w > 0 && w <= 1)
          width = w;
          depth = w;
        elseif (numel (w) == 2 && isnumeric (w) && isfinite (w)
                               && all (w > 0) && all (w <= 1))
          width = w(1);
          depth = w(2);
        else
          error ("bar3h: invalid value for optional argument 'width'.");
        endif
        varargin([1:2]) = [];
      elseif (strcmpi (tmp, "color"))
        if (numel (varargin) < 2)
          error ("bar3h: missing value for optional argument 'color'.");
        endif
        c = varargin{2};
        if (iscellstr (c))
          is_vc = all (cell2mat (cellfun (@(x) any (strcmpi (vc, x)), ...
                                          c, "UniformOutput", false)));
          if (is_vc)
            color = c;
          else
            error ("bar3h: invalid value for optional argument 'color'.");
          endif
        elseif (ischar (c) && isvector (c))
          if (any (strcmpi (c, vc)))
            color = c;
          else
            error ("bar3h: invalid value for optional argument 'color'.");
          endif
        elseif (isnumeric (c))
          if (size (c, 2) == 3 && all (c(:) >= 0) && all (c(:) <= 1))
            color = c;
          else
            error (["bar3h: numeric COLOR must be a 1x3 vector of an Nx3", ...
                    " matrix where each value is between 0 and 1 inclusive."]);
          endif
        else
          error ("bar3h: invalid value for optional argument 'color'.");
        endif
        varargin([1:2]) = [];
      elseif (strcmpi (tmp, "xlabel"))
        if (numel (varargin) < 2)
          error ("bar3h: missing value for optional argument 'xlabel'.");
        endif
        xlabel = varargin{2};
        if (! iscellstr (xlabel))
          error ("bar3h: invalid value for optional argument 'xlabel'.");
        endif
        varargin([1:2]) = [];
      elseif (strcmpi (tmp, "zlabel"))
        if (numel (varargin) < 2)
          error ("bar3h: missing value for optional argument 'zlabel'.");
        endif
        zlabel = varargin{2};
        if (! iscellstr (zlabel))
          error ("bar3h: invalid value for optional argument 'zlabel'.");
        endif
        varargin([1:2]) = [];
      else
        error ("bar3h: invalid optional argument.");
      endif
    elseif (iscellstr (tmp))
      is_vc = all (cell2mat (cellfun (@(x) any (strcmpi (vc, x)), ...
                                      tmp, "UniformOutput", false)));
      if (is_vc)
        color = tmp;
      else
        error ("bar3h: invalid value for optional COLOR argument.");
      endif
      varargin(1) = [];
    else
      error ("bar3h: invalid optional argument.");
    endif
  endwhile

  ## Get number of column bars from y input
  [nz, nx] = size (y);

  ## Check xlabel and zlabel
  if (! isempty (xlabel) && numel (xlabel) != nx)
    error ("bar3h: the elements in 'xlabel' must equal the columns in Z.");
  endif
  if (! isempty (zlabel) && numel (zlabel) != nz)
    error ("bar3h: the elements in 'zlabel' must equal the rows in Z.");
  endif

  ## Check COLOR for valid dimensions
  if (isempty (color))
    if (nx == 1)
      cargs = {'FaceColor', 'b'};
    elseif (nx == 0)
      defc = [0,0,1;1,0,1;0,1,0;1,1,0];
      fvcd = kron (defc([1:nx],:), ones (6 * nz, 1));
      cargs = {'FaceVertexCData', fvcd, 'FaceColor', 'flat'};
    else
      fvcd = kron ((1:nx)', ones (6 * nz, 1));
      cargs = {'FaceVertexCData', fvcd, 'FaceColor', 'flat', ...
               'CDataMapping', 'scaled'};
    endif
  elseif (isnumeric (color))
    if (size (color, 1) == numel (y))
      fvcd = kron (color, ones (6, 1));
      cargs = {'FaceVertexCData', fvcd, 'FaceColor', 'flat'};
    elseif (size (color, 1) == nx)
      fvcd = kron (color, ones (6 * nz, 1));
      cargs = {'FaceVertexCData', fvcd, 'FaceColor', 'flat'};
    elseif (size (color, 1) == nz)
      fvcd = repmat (kron (color, ones (6, 1)), nx, 1);
      cargs = {'FaceVertexCData', fvcd, 'FaceColor', 'flat'};
    elseif (size (color, 1) == 1)
      cargs = {'FaceVertexCData', color, 'FaceColor', 'flat'};
    endif
  elseif (ischar (color))
    cargs = {'FaceColor', color};
  elseif (iscellstr (color))

  endif

  ## Construct a "template" column-bar (8 vertices and 6 faces) centered at
  ## origin, with height = 1, and width along x-axis and depth along z-axis
  hw = width / 2;
  hd = depth / 2;
  ## Scale the bar's base when grouping together
  if (strcmpi (style, "grouped"))
    sc = nx + 1;
    hw = hw / sc;
    hd = hd / sc;
  endif
  [X, Y, Z] = ndgrid ([-hw, hw], [0, 1], [-hd, hd]);
  V = [X(:), Y(:), Z(:)];
  F = [1, 2, 4, 3; 5, 6, 8, 7; 1, 2, 6, 5; 3, 4, 8, 7; 1, 5, 7, 3; 2, 6, 8, 4];

  ## Replicate faces to the required number of bars
  increments = 0:8:8 * (nx * nz - 1);
  F = bsxfun (@plus, F, permute (increments, [1, 3, 2]));
  F = reshape (permute (F, [2, 1, 3]), 4, []).';

  ## Replicate vertices to the required number of bars
  [offsetX, offsetZ] = meshgrid (1:nx, 1:nz);
  offset = [offsetX(:), zeros(numel (y), 1), offsetZ(:)];
  #offset(:,3) = 0;
  V = bsxfun (@plus, V, permute (offset, [3, 2, 1]));
  V = reshape (permute (V, [2, 1, 3]), 3, []).';
  A = NaN;
  if (strcmpi (style, "detached"))
    ## Adjust bar heights according to values in y input
    V(:,2) = V(:,2) .* kron (y(:), ones (8,1));
  elseif (strcmpi (style, "grouped"))
    ## Adjust bar heights according to values in y input
    V(:,2) = V(:,2) .* kron (y(:), ones (8,1));
    ## Move groups along x axis
    V(:,1) = V(:,1) - kron (kron ([0:nx-1], ones (nz, 1))(:), ones (8,1));
    ## Move groups along z axis
    offset = [-nx+1:2:nx-1] * (hd / width);
    V(:,3) = V(:,3) + kron (kron (ones (1, nz), offset)(:), ones (8,1));
    nx = 1;
  elseif (strcmpi (style, "stacked"))
    ## Move groups along x axis
    V(:,1) = V(:,1) - kron (kron ([0:nx-1], ones (nz, 1))(:), ones (8,1));
    ## Adjust bar heights according to values in y input
    ZC = cumsum (y,2);
    Q1 = kron (ZC(:), ones (8,1));
    ZC(:,end) = [];
    ZC = [zeros(nz,1), ZC];
    Q2 = kron (ZC(:), ones (8,1));
    V(:,2) = V(:,2) .* Q1;
    idx = V(:,2) == 0;
    V(idx,2) = Q2(idx);
    ## Collapse x axis
    nx = 1;
  endif

  ## Draw column bars as patches specified by faces/vertices
  if (new_axes)
    ax = gca ();
  endif
  p = patch ('Faces', F, 'Vertices', V, 'EdgeColor', 'k', 'Parent', ax, cargs{:});

  ## Set view port and axes
  view (ax, 3);
  grid (ax, 'on');
  axis tight;
  xlim ([0.5, nx+0.5]);
  zlim ([0.5, nz+0.5]);
  set (ax, 'XTick', 1:nx, 'ZTick', 1:nz, 'Box', 'off', 'YDir', 'reverse');

  ## Fix aspect ratio so that bars appear square when rotating the bar plot
  if (nx > nz)
    set (ax, 'PlotBoxAspectRatio', [1, (sqrt(5)-1)/2, nz/nx]);
  elseif (nx < nz)
    set (ax, 'PlotBoxAspectRatio', [nx/nz, (sqrt(5)-1)/2, 1]);
  else
    set (ax, 'PlotBoxAspectRatio', [1, (sqrt(5)-1)/2, 1]);
  endif

  ## Add tick labels in axes (if requested)
  if (! isempty (xlabel))
    set (ax, 'XTickLabel', xlabel);
  endif
  if (! isempty (zlabel))
    set (ax, 'ZTickLabel', zlabel);
  elseif (! isempty (z) && ! isnan (z))
    zlabel = arrayfun (@(x) sprintf ('%d', x), z, 'UniformOutput', false);
    set (ax, 'ZTickLabel', zlabel);
  endif

  ## Return handle to patch object if requested
  if (nargout > 0)
    varargout{1} = p;
  endif

endfunction

%!demo
%! ## Plotting 5 bars in the same series.
%!
%! y = [50; 40; 30; 20; 10];
%! bar3h (y);

%!demo
%! ## Plotting 5 bars in different groups.
%!
%! y = [50, 40, 30, 20, 10];
%! bar3h (y);

%!demo
%! ## A 3D bar graph with each series corresponding to a column in y.
%!
%! y = [1, 4, 7; 2, 5, 8; 3, 6, 9; 4, 7, 10];
%! bar3h (y);

%!demo
%! ## Specify z-axis locations as tick names. z must be a column vector!
%!
%! z = [1950, 1960, 1970, 1980, 1990]';
%! y = [16, 8, 4, 2, 1]';
%! bar3h (z, y);

%!demo
%! ## Plot 3 series as a grouped plot without any space between the grouped bars
%!
%! y = [70 50 33 10; 75 55 35 15; 80 60 40 20];
%! bar3h (y, 1, 'grouped');

%!demo
%! ## Plot a stacked style 3D bar graph
%!
%! y = [19, 30, 21, 30; 40, 16, 32, 12];
%! b = bar3h (y, 0.5, 'stacked');

## Test input validation
%!error <bar3h: Z must be numeric.> bar3h ("A")
%!error <bar3h: Z must be numeric.> bar3h ({2,3,4,5})
%!error <bar3h: inconsistent size in Y and Z input arguments.> ...
%! bar3h ([1,2,3]', ones (2))
%!error <bar3h: WIDTH must be a scalar in the range> ...
%! bar3h ([1:5], 1.2)
%!error <bar3h: WIDTH must be a scalar in the range> ...
%! bar3h ([1:5]', ones (5), 1.2)
%!error <bar3h: numeric COLOR must be a 1x3 vector of an Nx3 matrix> ...
%! bar3h ([1:5]', ones (5), [0.8, 0.7])
%!error <bar3h: missing value for optional argument 'width'.> ...
%! bar3h (ones (5), 'width')
%!error <bar3h: invalid value for optional argument 'width'.> ...
%! bar3h (ones (5), 'width', 1.2)
%!error <bar3h: invalid value for optional argument 'width'.> ...
%! bar3h (ones (5), 'width', [0.8, 0.8, 0.8])
%!error <bar3h: missing value for optional argument 'color'.> ...
%! bar3h (ones (5), 'color')
%!error <bar3h: numeric COLOR must be a 1x3 vector of an Nx3 matrix> ...
%! bar3h (ones (5), 'color', [0.8, 0.8])
%!error <bar3h: invalid value for optional argument 'color'.> ...
%! bar3h (ones (5), 'color', "brown")
%!error <bar3h: invalid value for optional argument 'color'.> ...
%! bar3h (ones (5), 'color', {"r", "k", "c", "m", "brown"})
%!error <bar3h: missing value for optional argument 'xlabel'.> ...
%! bar3h (ones (5), 'xlabel')
%!error <bar3h: invalid value for optional argument 'xlabel'.> ...
%! bar3h (ones (5), 'xlabel', 4)
%!error <bar3h: missing value for optional argument 'zlabel'.> ...
%! bar3h (ones (5), 'zlabel')
%!error <bar3h: invalid value for optional argument 'zlabel'.> ...
%! bar3h (ones (5), 'zlabel', 4)
%!error <bar3h: invalid optional argument.> bar3h (ones (5), 'this', 4)
%!error <bar3h: the elements in 'xlabel' must equal the columns in Z.> ...
%! bar3h (ones (5), 'xlabel', {"A", "B", "C"})
%!error <bar3h: the elements in 'zlabel' must equal the rows in Z.> ...
%! bar3h (ones (5), 'zlabel', {"A", "B", "C"})
