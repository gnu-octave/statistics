## Copyright (C) 2021 Stefano Guidoni <ilguido@users.sf.net>
## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {} gscatter (@var{x}, @var{y}, @var{g})
## @deftypefnx {statistics} {} gscatter (@var{x}, @var{y}, @var{g}, @var{clr}, @var{sym}, @var{siz})
## @deftypefnx {statistics} {} gscatter (@dots{}, @var{doleg}, @var{xnam}, @var{ynam})
## @deftypefnx {statistics} @var{h} = gscatter (@dots{})
##
## Draw a scatter plot with grouped data.
##
## @code{gscatter} is a utility function to draw a scatter plot of @var{x} and
## @var{y}, according to the groups defined by @var{g}.  Input @var{x} and
## @var{y} are numeric vectors of the same size, while @var{g} is either a
## vector of the same size as @var{x} or a character matrix with the same number
## of rows as the size of @var{x}.  As a vector @var{g} can be numeric, logical,
## a character array, a string array (not implemented), a cell string or cell
## array.
##
## A number of optional inputs change the appearance of the plot:
## @itemize @bullet
## @item @var{"clr"}
## defines the color for each group; if not enough colors are defined by
## @var{"clr"}, @code{gscatter} cycles through the specified colors.  Colors can
## be defined as named colors, as rgb triplets or as indices for the current
## @code{colormap}.  The default value is a different color for each group,
## according to the current @code{colormap}.
##
## @item @var{"sym"}
## is a char array of symbols for each group; if not enough symbols are defined
## by @var{"sym"}, @code{gscatter} cycles through the specified symbols.
##
## @item @var{"siz"}
## is a numeric array of sizes for each group; if not enough sizes are defined
## by @var{"siz"}, @code{gscatter} cycles through the specified sizes.
##
## @item @var{"doleg"}
## is a boolean value to show the legend; it can be either @qcode{on} (default)
## or @qcode{off}.
##
## @item @var{"xnam"}
## is a character array, the name for the x axis.
##
## @item @var{"ynam"}
## is a character array, the name for the y axis.
## @end itemize
##
## Output @var{h} is an array of graphics handles to the @code{line} object of
## each group.
##
## @end deftypefn
##
## @seealso{scatter}

function h = gscatter (varargin)
  ## optional axes handle
  if (isaxes (varargin{1}))
    ## parameter is an axes handle
    hax = varargin{1};
    varargin = varargin(2:end);
    nargin--;
  endif

  ## check the input parameters
  if (nargin < 3)
    print_usage ();
  endif

  ##
  ## necessary parameters
  ##

  ## x coordinates
  if (isvector (varargin{1}) &&
      isnumeric (varargin{1}))
    x = varargin{1};
    n = numel (x);
  else
    error ("gscatter: x must be a numeric vector");
  endif

  ## y coordinates
  if (isvector (varargin{2}) &&
      isnumeric (varargin{2}))
    if (numel (varargin{2}) == n)
      y = varargin{2};
    else
      error ("gscatter: x and y must have the same size");
    endif
  else
    error ("gscatter: y must be a numeric vector");
  endif

  ## groups
  if (isrow (varargin{3}))
    varargin{3} = transpose (varargin{3});
  endif
  if (ismatrix (varargin{3}) && ischar (varargin{3}))
    varargin{3} = cellstr (varargin{3}); # char matrix to cellstr
  elseif (iscell (varargin{3}) && ! iscellstr (varargin{3}))
    varargin{3} = cell2mat (varargin{3}); # numeric cell to vector
  endif
  if (isvector (varargin{3})) # only numeric vectors or cellstr
    if (rows (varargin{3}) == n)
      gv = varargin{3};
      g_names = unique (gv, "rows");
      g_len = numel (g_names);
      if (iscellstr (g_names))
        for i = 1 : g_len
            g(find (strcmp(gv, g_names{i}))) = i;
        endfor
      else
        for i = 1 : g_len
            g(find (gv == g_names(i))) = i;
        endfor
      endif
    else
      error ("gscatter: g must have the same size as x and y");
    endif
  else
    error (["gscatter: g must be a numeric or logical or char vector, "...
      "or a cell or cellstr array, or a char matrix"]);
  endif

  ##
  ## optional parameters
  ##

  ## Note: this parameters are passed as they are to 'line',
  ##       the validity check is delegated to 'line'
  g_col = lines (g_len);
  g_size = 6 * ones (g_len, 1);
  g_sym = repmat ('o', 1, g_len);

  ## optional parameters for legend and axes labels
  do_legend = 1; # legend shown by default
  ## MATLAB compatibility: by default MATLAB uses the variable name as
  ## label for either axis
  mygetname = @(x) inputname(1); # to retrieve the name of a variable
  x_nam = mygetname (varargin{1}); # this should retrieve the name of the var,
  y_nam = mygetname (varargin{2}); # but it does not work

  ## parameters are all in fixed positions
  for i = 4 : nargin
    switch (i)
      case 4
        ## colours
        c_list = varargin{4};
        if (isrow (c_list))
          c_list = transpose (c_list);
        endif
        c_list_len = rows (c_list);

        g_col = repmat (c_list, ceil (g_len / c_list_len));
      case {5, 6}
        ## size and symbols
        s_list = varargin{i};
        s_list_len = length (s_list);

        g_tmp = repmat (s_list, ceil (g_len / s_list_len));
        if (i == 6)
          g_size = g_tmp;
        else
          g_sym = g_tmp;
        endif
      case 7
        ## legend
        switch (lower (varargin{7}))
          case "on"
            do_legend = 1;
          case "off"
            do_legend = 0;
          otherwise
            error ("gscatter: invalid dolegend parameter '%s'", varargin{7});
        endswitch
      case {8, 9}
        ## x and y label
        if (! ischar (varargin{i}) && ! isvector (varargin{i}))
          error ("gscatter: xnam and ynam must be strings");
        endif
        if (i == 8)
          x_nam = varargin{8};
        else
          y_nam = varargin{9};
        endif
    endswitch
  endfor


  ## scatter plot with grouping
  if (! exist ("hax", "var"))
    hax = gca ();
  endif

  ## return value
  h = [];

  hold on;
  for i = 1 : g_len
    idcs = find (g == i);
    h(i) = line (hax, x(idcs), y(idcs), "linestyle", "none", ...
            "markersize", g_size(i), "color", g_col(i,:), "marker", g_sym(i));
  endfor
  if (do_legend)
    if (isnumeric (g_names))
      g_names = num2str (g_names);
    endif
    warning ("off", "Octave:legend:unimplemented-location", "local");
    legend (hax, g_names, "location", "best");
  endif
  xlabel (hax, x_nam);
  ylabel (hax, y_nam);
  hold off;
endfunction

## demonstration
%!demo
%! load fisheriris;
%! X = meas(:,3:4);
%! cidcs = kmeans (X, 3, "Replicates", 5);
%! gscatter (X(:,1), X(:,2), cidcs, [.75 .75 0; 0 .75 .75; .75 0 .75], "os^");
%! title ("Fisher's iris data");

## input tests
%!error gscatter();
%!error gscatter([1]);
%!error gscatter([1], [2]);
%!error <x must be a numeric vector> gscatter('abc', [1 2 3], [1]);
%!error <x and y must have the same size> gscatter([1 2 3], [1 2], [1]);
%!error <y must be a numeric vector> gscatter([1 2 3], 'abc', [1]);
%!error <g must have the same size as x and y> gscatter([1 2], [1 2], [1]);
%!error <invalid dolegend> gscatter([1 2], [1 2], [1 2], 'rb', 'so', 12, 'xxx');
## testing demonstration
%!shared visibility_setting
%! visibility_setting = get (0, "DefaultFigureVisible");
%!test
%! set (0, "DefaultFigureVisible", "off");
%! load fisheriris;
%! X = meas(:,3:4);
%! cidcs = kmeans (X, 3, "Replicates", 5);
%! gscatter (X(:,1), X(:,2), cidcs, [.75 .75 0; 0 .75 .75; .75 0 .75], "os^");
%! title ("Fisher's iris data");
%! set (0, "DefaultFigureVisible", visibility_setting);
