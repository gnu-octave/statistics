## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {} gplotmatrix (@var{x}, @var{y}, @var{group})
## @deftypefnx {statistics} {} gplotmatrix (@var{x}, [], @var{group})
## @deftypefnx {statistics} {} gplotmatrix (@var{x}, @var{y}, @var{group}, @var{clr}, @var{sym}, @var{siz})
## @deftypefnx {statistics} {} gplotmatrix (@dots{}, @var{doleg}, @var{dispopt})
## @deftypefnx {statistics} {} gplotmatrix (@dots{}, @var{doleg}, @var{dispopt}, @var{xnam}, @var{ynam})
## @deftypefnx {statistics} {} gplotmatrix (@var{parent}, @dots{})
## @deftypefnx {statistics} {[@var{h}, @var{ax}, @var{bigax}] =} gplotmatrix (@dots{})
##
## Create a matrix of scatter plots grouped by a categorical variable.
##
## @code{gplotmatrix (@var{x}, @var{y}, @var{group})} creates a matrix of
## scatter plots.  Each subplot in the resulting figure is a scatter plot of a
## column of @var{x} against a column of @var{y}.  If @var{x} is @code{n}-by-@var{p}
## and @var{y} is @code{n}-by-@var{q}, the resulting figure holds a
## @var{q}-by-@var{p} grid of subplots; the subplot in row @var{i} and column
## @var{j} plots @code{@var{x}(:,j)} on the horizontal axis against
## @code{@var{y}(:,i)} on the vertical axis.  Points are grouped and colored
## according to @var{group}, which is a grouping variable (numeric, logical,
## character, string, or cell array of strings) with one entry per row of
## @var{x}.
##
## @code{gplotmatrix (@var{x}, [], @var{group})} is equivalent to
## @code{gplotmatrix (@var{x}, @var{x}, @var{group})} except that the diagonal of
## the @var{p}-by-@var{p} grid is replaced by grouped histograms of the columns
## of @var{x}.
##
## The appearance of the plot is controlled by further positional arguments:
##
## @table @asis
## @item @var{clr}
## Marker colors, given as a character vector of color specifiers (e.g.
## @qcode{"rgb"}) or as a matrix of RGB triplets, one row per group.  Colors
## cycle if fewer are supplied than there are groups.
##
## @item @var{sym}
## Marker symbols, given as a character vector (e.g. @qcode{"o+x"}); defaults to
## @qcode{"."}.  Symbols cycle if fewer are supplied than there are groups.
##
## @item @var{siz}
## Marker sizes, given as a numeric vector.  Sizes cycle if fewer are supplied
## than there are groups.
##
## @item @var{doleg}
## Either @qcode{"on"} (default) to display a legend of the groups or
## @qcode{"off"} to suppress it.
##
## @item @var{dispopt}
## Controls the diagonal of the grid when @var{y} is empty: @qcode{"stairs"}
## (default) for grouped stairstep histograms, @qcode{"hist"} or
## @qcode{"grpbars"} for grouped bar histograms, @qcode{"none"} to leave the
## diagonal empty, or @qcode{"variable"} to write the variable names on the
## diagonal.
##
## @item @var{xnam}, @var{ynam}
## Character vectors or cell arrays of strings giving the names of the columns of
## @var{x} and @var{y}, used to label the outer axes.
## @end table
##
## An optional leading @var{parent} argument (a figure or uipanel handle) selects
## the container for the plot.
##
## The optional outputs are @var{h}, an array of handles to the plotted objects
## with size @var{ny}-by-@var{p}-by-@var{k} (where @var{ny} is the number of rows
## of the grid and @var{k} the number of groups); @var{ax}, the matrix of handles
## to the subplot axes (with an extra row of hidden axes for the diagonal
## histograms); and @var{bigax}, the handle to the invisible enclosing axes used
## for titles and labels.
##
## @seealso{gscatter, plotmatrix, grpstats}
## @end deftypefn

function [h, ax, bigax] = gplotmatrix (varargin)

  ## Optional leading parent (figure or uipanel) handle
  parent = [];
  if (numel (varargin) > 0 && isscalar (varargin{1}) && ishghandle (varargin{1})
      && any (strcmp (get (varargin{1}, "type"), {"figure", "uipanel"})))
    parent = varargin{1};
    varargin(1) = [];
  endif

  if (numel (varargin) < 1)
    print_usage ();
  endif

  ## Positional arguments
  X = varargin{1};
  Y = []; group = []; clr = []; sym = []; siz = [];
  doleg = "on"; dispopt = "stairs"; xnam = []; ynam = [];
  nv = numel (varargin);
  if (nv >= 2), Y = varargin{2}; endif
  if (nv >= 3), group = varargin{3}; endif
  if (nv >= 4), clr = varargin{4}; endif
  if (nv >= 5), sym = varargin{5}; endif
  if (nv >= 6), siz = varargin{6}; endif
  if (nv >= 7 && ! isempty (varargin{7})), doleg = varargin{7}; endif
  if (nv >= 8 && ! isempty (varargin{8})), dispopt = varargin{8}; endif
  if (nv >= 9), xnam = varargin{9}; endif
  if (nv >= 10), ynam = varargin{10}; endif
  if (nv > 10)
    error ("gplotmatrix: too many input arguments.");
  endif

  if (! isnumeric (X) || ! isreal (X) || ! ismatrix (X) || ndims (X) > 2)
    error ("gplotmatrix: X must be a real numeric matrix.");
  endif
  n = rows (X);
  p = columns (X);

  ## Y empty selects the self plot with histograms on the diagonal
  do_hist = isempty (Y);
  if (do_hist)
    Y = X;
  elseif (! isnumeric (Y) || ! isreal (Y) || ndims (Y) > 2)
    error ("gplotmatrix: Y must be a real numeric matrix.");
  elseif (rows (Y) != n)
    error ("gplotmatrix: X and Y must have the same number of rows.");
  endif
  ny = columns (Y);

  ## Grouping variable
  if (isempty (group))
    gidx = ones (n, 1);
    gnames = {"1"};
  else
    [gidx, gnames] = grp2idx (group);
    if (numel (gidx) != n)
      error ("gplotmatrix: GROUP must have one entry per row of X.");
    endif
  endif
  k = numel (gnames);

  ## Per-group color, symbol, and size
  gcol = expand_color (clr, k);
  gsym = expand_sym (sym, k);
  gsiz = expand_siz (siz, k);

  ## Names for the outer axes
  xnam = name_list (xnam, p);
  ynam = name_list (ynam, ny);

  ## Container to draw into
  if (isempty (parent))
    parent = gcf ();
  endif
  clf (parent);

  ## Grid geometry (row 1 at the top, column 1 at the left)
  Lm = 0.10; Bm = 0.10; Wt = 0.86; Ht = 0.86; gap = 0.015;
  cw = Wt / p;
  chh = Ht / ny;

  ax = zeros (ny, p);
  histax = zeros (1, p);
  h = [];
  glines = [];
  legax = [];

  for r = 1:ny
    for c = 1:p
      x0 = Lm + (c - 1) * cw;
      y0 = Bm + (ny - r) * chh;
      pos = [x0 + gap/2, y0 + gap/2, cw - gap, chh - gap];
      a = axes ("parent", parent, "position", pos, "box", "on", ...
                "nextplot", "add");
      ax(r,c) = a;
      if (do_hist && r == c)
        ## Diagonal: grouped histogram drawn in an overlay axes
        ah = axes ("parent", parent, "position", pos, "color", "none", ...
                   "nextplot", "add", "xtick", [], "ytick", []);
        histax(c) = ah;
        for l = 1:k
          h(r,c,l) = draw_diag (ah, X(gidx == l, c), dispopt, ...
                                gcol(l,:), xnam{c});
        endfor
        set (a, "xtick", [], "ytick", []);
      else
        cell_lines = zeros (1, k);
        for l = 1:k
          idx = (gidx == l);
          hl = line (a, X(idx, c), Y(idx, r), "linestyle", "none", ...
                     "marker", gsym(l), "markersize", gsiz(l), ...
                     "color", gcol(l,:));
          h(r,c,l) = hl;
          cell_lines(l) = hl;
        endfor
        if (isempty (glines))
          glines = cell_lines;
          legax = a;
        endif
      endif
      ## Only the outer edges carry tick labels
      if (r != ny)
        set (a, "xticklabel", []);
      endif
      if (c != 1)
        set (a, "yticklabel", []);
      endif
      if (r == ny && ! isempty (xnam{c}))
        xlabel (a, xnam{c});
      endif
      if (c == 1 && ! isempty (ynam{r}))
        ylabel (a, ynam{r});
      endif
    endfor
  endfor

  ## Enclosing invisible axes for titles and overall labels
  bigax = axes ("parent", parent, "position", [Lm, Bm, Wt, Ht], ...
                "visible", "off", "xtick", [], "ytick", []);

  ## Legend of the groups
  if (strcmpi (doleg, "on") && k > 1 && ! isempty (glines))
    warning ("off", "Octave:legend:unimplemented-location", "local");
    legend (legax, glines, gnames, "location", "best");
  endif

  if (nargout == 0)
    clear h ax bigax;
  elseif (do_hist)
    ax = [ax; histax];
  endif

endfunction

## Expand a color specification to a k-by-3 matrix of RGB triplets.
function gcol = expand_color (clr, k)
  if (isempty (clr))
    base = lines (k);
  elseif (ischar (clr))
    base = zeros (numel (clr), 3);
    for i = 1:numel (clr)
      base(i,:) = char2rgb (clr(i));
    endfor
  elseif (isnumeric (clr) && columns (clr) == 3)
    base = clr;
  else
    error ("gplotmatrix: CLR must be a color string or an n-by-3 RGB matrix.");
  endif
  idx = mod (0:k-1, rows (base)) + 1;
  gcol = base(idx, :);
endfunction

function rgb = char2rgb (ch)
  switch (ch)
    case "r", rgb = [1 0 0];
    case "g", rgb = [0 1 0];
    case "b", rgb = [0 0 1];
    case "c", rgb = [0 1 1];
    case "m", rgb = [1 0 1];
    case "y", rgb = [1 1 0];
    case "k", rgb = [0 0 0];
    case "w", rgb = [1 1 1];
    otherwise
      error ("gplotmatrix: unknown color '%s'.", ch);
  endswitch
endfunction

## Expand a marker specification to a k-element character vector.
function gsym = expand_sym (sym, k)
  if (isempty (sym))
    sym = ".";
  elseif (! ischar (sym))
    error ("gplotmatrix: SYM must be a character vector of markers.");
  endif
  idx = mod (0:k-1, numel (sym)) + 1;
  gsym = sym(idx);
endfunction

## Expand a size specification to a k-element numeric vector.
function gsiz = expand_siz (siz, k)
  if (isempty (siz))
    siz = 6;
  elseif (! isnumeric (siz) || ! isreal (siz))
    error ("gplotmatrix: SIZ must be a numeric vector of marker sizes.");
  endif
  idx = mod (0:k-1, numel (siz)) + 1;
  gsiz = siz(idx);
endfunction

## Normalize a name specification to a cell array of p strings.
function nm = name_list (names, p)
  if (isempty (names))
    nm = repmat ({""}, 1, p);
  elseif (ischar (names))
    nm = cellstr (names).';
  elseif (iscellstr (names))
    nm = names(:).';
  else
    error ("gplotmatrix: variable names must be strings.");
  endif
  if (numel (nm) < p)
    nm(end+1:p) = {""};
  endif
endfunction

## Draw one group's contribution to a diagonal histogram cell.
function ho = draw_diag (ah, data, dispopt, col, vname)
  switch (lower (dispopt))
    case "none"
      ho = line (ah, NaN, NaN, "linestyle", "none");
    case "variable"
      ho = text (ah, 0.5, 0.5, vname, "parent", ah, ...
                 "units", "normalized", "horizontalalignment", "center");
    otherwise
      if (isempty (data) || all (isnan (data)))
        ho = line (ah, NaN, NaN, "linestyle", "none");
        return;
      endif
      [nn, xx] = hist (data, 10);
      switch (lower (dispopt))
        case {"hist", "grpbars"}
          ho = bar (ah, xx, nn, 1.0, "facecolor", col, "edgecolor", col);
        otherwise  # "stairs"
          e = xx(1:end-1) + diff (xx) / 2;
          xe = [xx(1)-(xx(2)-xx(1))/2, e, xx(end)+(xx(end)-xx(end-1))/2];
          ho = stairs (ah, xe, [nn, nn(end)], "color", col);
      endswitch
  endswitch
endfunction

%!demo
%! ## Grouped scatter-plot matrix of Fisher's iris measurements.
%!
%! load fisheriris;
%! gplotmatrix (meas, [], species);

%!demo
%! ## Two sets of variables plotted against each other by group.
%!
%! load fisheriris;
%! gplotmatrix (meas(:,1:2), meas(:,3:4), species);

## Test output shapes and scatter orientation
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [10 20; 11 25; 12 21; 13 28; 14 23; 15 29];
%!   g = [1 1 1 2 2 2]';
%!   [h, ax, bigax] = gplotmatrix (X, [], g);
%!   assert (size (h), [2, 2, 2]);
%!   assert (size (ax), [3, 2]);
%!   assert (isscalar (bigax) && isaxes (bigax));
%!   assert (get (h(1,2,1), "xdata"), [20 25 21]);
%!   assert (get (h(1,2,1), "ydata"), [10 11 12]);
%!   assert (get (h(2,1,1), "xdata"), [10 11 12]);
%!   assert (get (h(2,1,1), "ydata"), [20 25 21]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [10 20; 11 25; 12 21; 13 28; 14 23; 15 29];
%!   Y = [100 200 300; 110 250 280; 120 210 260; ...
%!        130 280 240; 140 230 220; 150 290 210];
%!   g = [1 1 1 2 2 2]';
%!   [h, ax] = gplotmatrix (X, Y, g);
%!   assert (size (h), [3, 2, 2]);
%!   assert (size (ax), [3, 2]);
%!   assert (get (h(1,2,1), "xdata"), [20 25 21]);
%!   assert (get (h(1,2,1), "ydata"), [100 110 120]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # runs without a grouping variable
%! hf = figure ("visible", "off");
%! unwind_protect
%!   [h, ax] = gplotmatrix (randn (20, 3), [], []);
%!   assert (size (h), [3, 3]);
%!   assert (size (ax), [4, 3]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error <Invalid call to gplotmatrix> gplotmatrix ()
%!error <gplotmatrix: X must be a real numeric matrix.> gplotmatrix ({1})
%!error <gplotmatrix: X and Y must have the same number of rows.> ...
%! gplotmatrix (ones (5, 2), ones (4, 2), ones (5, 1))
%!error <gplotmatrix: GROUP must have one entry per row of X.> ...
%! gplotmatrix (ones (5, 2), [], ones (4, 1))
%!error <gplotmatrix: too many input arguments.> ...
%! gplotmatrix (ones (5, 2), [], ones (5, 1), "r", ".", 6, "on", "hist", ...
%!              "a", "b", "c")
