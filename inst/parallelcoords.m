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
## @deftypefn  {statistics} {} parallelcoords (@var{x})
## @deftypefnx {statistics} {} parallelcoords (@var{x}, @var{name}, @var{value}, @dots{})
## @deftypefnx {statistics} {} parallelcoords (@var{ax}, @dots{})
## @deftypefnx {statistics} {@var{h} =} parallelcoords (@dots{})
##
## Create a parallel coordinates plot of the multivariate data in @var{x}.
##
## @code{parallelcoords (@var{x})} plots each observation (row) of the
## @code{n}-by-@code{p} matrix @var{x} as a line connecting the values of its
## @var{p} coordinates, which are placed at the equally spaced horizontal
## positions @code{1, 2, @dots{}, p}.
##
## The following name/value pairs are accepted:
##
## @table @asis
## @item @qcode{"Group"}
## A grouping variable (numeric, logical, character, string, or cell array of
## strings) with one entry per row of @var{x}.  Lines are colored by group.
##
## @item @qcode{"Standardize"}
## Controls how the columns of @var{x} are transformed before plotting:
## @qcode{"off"} (default) uses the raw data, @qcode{"on"} centers and scales
## each column to zero mean and unit standard deviation, @qcode{"PCA"} uses the
## principal component scores, and @qcode{"PCAStd"} uses the principal component
## scores of the standardized data.
##
## @item @qcode{"Quantile"}
## A scalar @var{alpha} in the interval @math{(0,1)}.  Instead of one line per
## observation, only three lines per group are drawn: the coordinate-wise median
## and the @var{alpha} and @math{1-}@var{alpha} quantiles of the group.
##
## @item @qcode{"Labels"}
## A character array or cell array of strings giving the tick labels for the
## coordinate axis.
## @end table
##
## @code{parallelcoords (@var{ax}, @dots{})} plots into the axes @var{ax}.
##
## The optional output @var{h} is a vector of handles to the plotted lines: one
## per observation, or three per group when @qcode{"Quantile"} is used.
##
## @seealso{andrewsplot, glyphplot, pca}
## @end deftypefn

function h = parallelcoords (varargin)

  ## Optional leading axes handle
  hax = [];
  if (numel (varargin) > 0 && isaxes (varargin{1}))
    hax = varargin{1};
    varargin(1) = [];
  endif

  if (numel (varargin) < 1)
    print_usage ();
  endif

  X = varargin{1};
  varargin(1) = [];
  if (! isnumeric (X) || ! isreal (X) || ndims (X) > 2)
    error ("parallelcoords: X must be a real numeric matrix.");
  endif
  n = rows (X);

  ## Parse name/value options
  group = [];
  standardize = "off";
  alpha = [];
  labels = {};
  if (mod (numel (varargin), 2) != 0)
    error ("parallelcoords: name/value arguments must come in pairs.");
  endif
  for i = 1:2:numel (varargin)
    name = varargin{i};
    value = varargin{i+1};
    if (! ischar (name))
      error ("parallelcoords: property names must be strings.");
    endif
    switch (lower (name))
      case "group"
        group = value;
      case "standardize"
        standardize = value;
      case "quantile"
        alpha = value;
      case "labels"
        labels = value;
      otherwise
        error ("parallelcoords: unknown property '%s'.", name);
    endswitch
  endfor

  ## Standardize the data
  switch (lower (standardize))
    case "off"
      Z = X;
    case "on"
      Z = zscore (X);
    case "pca"
      [~, Z] = pca (X);
    case "pcastd"
      [~, Z] = pca (zscore (X));
    otherwise
      error ("parallelcoords: invalid Standardize option '%s'.", standardize);
  endswitch
  p = columns (Z);
  cx = 1:p;

  ## Grouping
  if (isempty (group))
    gidx = ones (n, 1);
    gnames = {"1"};
  else
    [gidx, gnames] = grp2idx (group);
    if (numel (gidx) != n)
      error ("parallelcoords: GROUP must have one entry per row of X.");
    endif
  endif
  k = numel (gnames);
  gcol = lines (k);

  if (isempty (hax))
    hax = newplot ();
  else
    newplot (hax);
  endif
  old_hold = ishold (hax);
  hold (hax, "on");

  h = [];
  if (isempty (alpha))
    ## One line per observation, colored by group
    for i = 1:n
      h(end+1) = line (hax, cx, Z(i,:), "color", gcol(gidx(i),:));
    endfor
  else
    ## Coordinate-wise median and alpha / 1-alpha quantile lines per group
    if (! (isscalar (alpha) && isreal (alpha) && alpha > 0 && alpha < 1))
      error ("parallelcoords: Quantile ALPHA must be a scalar in (0,1).");
    endif
    for g = 1:k
      Zg = Z(gidx == g, :);
      med = median (Zg, 1);
      lo = quantile (Zg, alpha, 1);
      hi = quantile (Zg, 1 - alpha, 1);
      h(end+1) = line (hax, cx, med, "color", gcol(g,:), "linewidth", 2);
      h(end+1) = line (hax, cx, lo, "color", gcol(g,:), "linestyle", "--");
      h(end+1) = line (hax, cx, hi, "color", gcol(g,:), "linestyle", "--");
    endfor
  endif

  ## Coordinate axis ticks and labels
  set (hax, "xtick", cx);
  set (hax, "xlim", [1, max(p, 2)]);
  if (! isempty (labels))
    if (ischar (labels))
      labels = cellstr (labels);
    endif
    set (hax, "xticklabel", labels);
  endif
  xlabel (hax, "Coordinate");
  ylabel (hax, "Coordinate Value");

  if (k > 1 && ! isempty (group))
    warning ("off", "Octave:legend:unimplemented-location", "local");
    if (isempty (alpha))
      [~, first] = unique (gidx, "first");
      legend (hax, h(sort (first)), gnames, "location", "best");
    else
      legend (hax, h(1:3:end), gnames, "location", "best");
    endif
  endif

  if (! old_hold)
    hold (hax, "off");
  endif

  if (nargout == 0)
    clear h;
  endif

endfunction

%!demo
%! ## Parallel coordinates plot of Fisher's iris data, grouped by species.
%!
%! load fisheriris;
%! parallelcoords (meas, "Group", species, "Labels", ...
%!                 {"SL", "SW", "PL", "PW"});

%!demo
%! ## The same data with median and quartile lines per species.
%!
%! load fisheriris;
%! parallelcoords (meas, "Group", species, "Quantile", 0.25, ...
%!                 "Standardize", "on");

## Test output
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [1 2 3; 4 5 6; 7 8 9; 2 1 5];
%!   h = parallelcoords (X);
%!   assert (numel (h), 4);
%!   assert (get (h(1), "xdata"), [1 2 3]);
%!   assert (get (h(1), "ydata"), [1 2 3]);
%!   assert (get (h(2), "ydata"), [4 5 6]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # Standardize "on" (z-score each column)
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [1 2 3; 4 5 6; 7 8 9; 2 1 5];
%!   h = parallelcoords (X, "Standardize", "on");
%!   assert (get (h(1), "ydata"), [-0.94491, -0.63246, -1.1], 1e-4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # Standardize "pca"
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [1 2 3; 4 5 6; 7 8 9; 2 1 5];
%!   h = parallelcoords (X, "Standardize", "pca");
%!   assert (get (h(1), "ydata")(1:2), [-4.1082, 0.9670], 1e-4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # Quantile mode: median, alpha, 1-alpha lines
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [1 2 3; 4 5 6; 7 8 9; 2 1 5];
%!   h = parallelcoords (X, "Quantile", 0.25);
%!   assert (numel (h), 3);
%!   assert (get (h(1), "ydata"), [3, 3.5, 5.5], 1e-12);
%!   assert (get (h(2), "ydata"), [1.5, 1.5, 4], 1e-12);
%!   assert (get (h(3), "ydata"), [5.5, 6.5, 7.5], 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # coordinate tick labels
%! hf = figure ("visible", "off");
%! unwind_protect
%!   h = parallelcoords ([1 2 3; 4 5 6], "Labels", {"a", "b", "c"});
%!   assert (get (gca, "xtick"), [1 2 3]);
%!   assert (get (gca, "xticklabel"), {"a"; "b"; "c"});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error <Invalid call to parallelcoords> parallelcoords ()
%!error <parallelcoords: X must be a real numeric matrix.> parallelcoords ({1})
%!error <parallelcoords: name/value arguments must come in pairs.> ...
%! parallelcoords (ones (3, 2), "Group")
%!error <parallelcoords: unknown property 'bogus'.> ...
%! parallelcoords (ones (3, 2), "bogus", 1)
%!error <parallelcoords: invalid Standardize option 'xxx'.> ...
%! parallelcoords (ones (3, 2), "Standardize", "xxx")
%!error <parallelcoords: GROUP must have one entry per row of X.> ...
%! parallelcoords (ones (3, 2), "Group", [1 2])
%!error <parallelcoords: Quantile ALPHA must be a scalar in .0,1..> ...
%! parallelcoords (ones (3, 2), "Quantile", 0)
