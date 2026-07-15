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
## @deftypefn  {statistics} {} scatterhist (@var{x}, @var{y})
## @deftypefnx {statistics} {} scatterhist (@var{x}, @var{y}, @var{name}, @var{value}, @dots{})
## @deftypefnx {statistics} {@var{h} =} scatterhist (@dots{})
##
## Create a scatter plot of @var{x} and @var{y} with marginal histograms.
##
## @code{scatterhist (@var{x}, @var{y})} draws a scatter plot of the vectors
## @var{x} and @var{y} in a central set of axes, with a histogram of @var{x}
## above it and a histogram of @var{y} to its right.  @var{x} and @var{y} must
## be vectors of the same length; @qcode{NaN} values are removed pairwise from
## the scatter plot and individually from each marginal histogram.
##
## The following name/value pairs are accepted:
##
## @table @asis
## @item @qcode{"Group"}
## A grouping variable (numeric, logical, character, string, or cell array of
## strings) with one entry per point.  The scatter points and the marginal
## histograms are separated and colored by group.
##
## @item @qcode{"NBins"}
## The number of bins for the marginal histograms, either a scalar applied to
## both or a two-element vector @code{[nx ny]}.  The default is chosen by Scott's
## rule.
##
## @item @qcode{"Kernel"}
## @qcode{"off"} (default) draws histograms for the marginals; @qcode{"on"} or
## @qcode{"overlay"} draws kernel density estimates instead.
##
## @item @qcode{"Location"}
## Position of the marginal plots: @qcode{"SouthWest"} (default),
## @qcode{"SouthEast"}, @qcode{"NorthEast"}, or @qcode{"NorthWest"}.
##
## @item @qcode{"Legend"}
## @qcode{"on"} or @qcode{"off"} to show or hide the group legend.  The default
## is @qcode{"on"} when a grouping variable is supplied.
##
## @item @qcode{"Marker"}, @qcode{"MarkerSize"}
## The marker symbol(s) and size(s) for the scatter points, cycled over the
## groups.
## @end table
##
## The optional output @var{h} is a three-element vector of axes handles: the
## central scatter axes, the axes of the @var{x} (horizontal) histogram, and the
## axes of the @var{y} (vertical) histogram.
##
## @seealso{gscatter, scatter, hist, ksdensity}
## @end deftypefn

function h = scatterhist (varargin)

  if (numel (varargin) < 2)
    print_usage ();
  endif

  x = varargin{1}(:);
  y = varargin{2}(:);
  varargin(1:2) = [];
  if (! isnumeric (x) || ! isreal (x) || ! isnumeric (y) || ! isreal (y))
    error ("scatterhist: X and Y must be real numeric vectors.");
  endif
  if (numel (x) != numel (y))
    error ("scatterhist: X and Y must have the same length.");
  endif

  ## Parse name/value options
  group = [];
  nbins = [];
  kernel = "off";
  location = "southwest";
  legend_opt = "";
  marker = "o";
  markersize = 6;
  if (mod (numel (varargin), 2) != 0)
    error ("scatterhist: name/value arguments must come in pairs.");
  endif
  for i = 1:2:numel (varargin)
    name = varargin{i};
    value = varargin{i+1};
    if (! ischar (name))
      error ("scatterhist: property names must be strings.");
    endif
    switch (lower (name))
      case "group"
        group = value;
      case "nbins"
        nbins = value;
      case "kernel"
        kernel = value;
      case "location"
        location = value;
      case "legend"
        legend_opt = value;
      case "marker"
        marker = value;
      case "markersize"
        markersize = value;
      otherwise
        error ("scatterhist: unknown property '%s'.", name);
    endswitch
  endfor

  n = numel (x);
  if (isempty (group))
    gidx = ones (n, 1);
    gnames = {"1"};
  else
    [gidx, gnames] = grp2idx (group);
    if (numel (gidx) != n)
      error ("scatterhist: GROUP must have one entry per point.");
    endif
  endif
  k = numel (gnames);
  gcol = lines (k);
  do_kernel = any (strcmpi (kernel, {"on", "overlay"}));

  ## Number of bins for the marginals (Scott's rule by default)
  if (isempty (nbins))
    nbx = scott_nbins (x);
    nby = scott_nbins (y);
  elseif (isscalar (nbins))
    nbx = nby = nbins;
  else
    nbx = nbins(1);
    nby = nbins(2);
  endif

  ## Axes layout: scatter lower-left, x-hist on top, y-hist on the right
  cf = gcf ();
  clf (cf);
  pos_s = [0.10, 0.10, 0.60, 0.60];
  pos_x = [0.10, 0.72, 0.60, 0.20];
  pos_y = [0.72, 0.10, 0.20, 0.60];
  ax_s = axes ("parent", cf, "position", pos_s, "box", "on", "nextplot", "add");
  ax_x = axes ("parent", cf, "position", pos_x, "nextplot", "add");
  ax_y = axes ("parent", cf, "position", pos_y, "nextplot", "add");

  ## Central scatter plot, grouped
  scat = zeros (1, k);
  for g = 1:k
    idx = (gidx == g) & ! isnan (x) & ! isnan (y);
    scat(g) = line (ax_s, x(idx), y(idx), "linestyle", "none", ...
                    "marker", marker(mod (g - 1, numel (marker)) + 1), ...
                    "markersize", markersize(mod (g - 1, numel (markersize)) + 1), ...
                    "color", gcol(g,:));
  endfor
  xlabel (ax_s, inputname (1));
  ylabel (ax_s, inputname (2));

  ## Marginal for x along the top, for y along the right
  multi = (k > 1);
  for g = 1:k
    xg = x(gidx == g & ! isnan (x));
    yg = y(gidx == g & ! isnan (y));
    marginal (ax_x, xg, nbx, gcol(g,:), do_kernel, multi, false);
    marginal (ax_y, yg, nby, gcol(g,:), do_kernel, multi, true);
  endfor

  ## Share the data axes with the scatter and tidy the marginal axes
  set (ax_x, "xlim", get (ax_s, "xlim"), "xtick", [], "ytick", []);
  set (ax_y, "ylim", get (ax_s, "ylim"), "xtick", [], "ytick", []);
  axis (ax_x, "off");
  axis (ax_y, "off");

  ## Legend
  show_legend = (! isempty (group)) && ! strcmpi (legend_opt, "off");
  if (show_legend && k > 1)
    warning ("off", "Octave:legend:unimplemented-location", "local");
    legend (ax_s, scat, gnames, "location", "best");
  endif

  if (nargout > 0)
    h = [ax_s, ax_x, ax_y];
  endif

endfunction

## Scott's rule for the number of histogram bins.
function nb = scott_nbins (v)
  v = v(! isnan (v));
  n = numel (v);
  if (n < 2)
    nb = 1;
    return;
  endif
  bw = 3.5 * std (v) * n ^ (-1/3);
  if (bw <= 0)
    nb = 1;
  else
    nb = max (1, ceil ((max (v) - min (v)) / bw));
  endif
endfunction

## Draw one group's marginal, either as a histogram or a kernel density.
## When horiz is true the marginal runs along the vertical (y) axis.
function marginal (ax, v, nb, col, do_kernel, multi, horiz)
  if (isempty (v))
    return;
  endif
  if (do_kernel)
    [f, u] = ksdensity (v);
    if (horiz)
      line (ax, f, u, "color", col);
    else
      line (ax, u, f, "color", col);
    endif
  else
    [nn, cc] = hist (v, nb);
    if (multi)
      ## stairstep outline
      e = cc(1:end-1) + diff (cc) / 2;
      ce = [cc(1)-(cc(2)-cc(1))/2, e, cc(end)+(cc(end)-cc(end-1))/2];
      ne = [nn, nn(end)];
      if (horiz)
        stairs (ax, ne, ce, "color", col);
      else
        stairs (ax, ce, ne, "color", col);
      endif
    else
      if (horiz)
        barh (ax, cc, nn, 1.0, "facecolor", col, "edgecolor", col);
      else
        bar (ax, cc, nn, 1.0, "facecolor", col, "edgecolor", col);
      endif
    endif
  endif
endfunction

%!demo
%! ## Scatter plot of two iris measurements with marginal histograms by species.
%!
%! load fisheriris;
%! scatterhist (meas(:,1), meas(:,2), "Group", species);

%!demo
%! ## Marginal kernel density estimates instead of histograms.
%!
%! load fisheriris;
%! scatterhist (meas(:,3), meas(:,4), "Group", species, "Kernel", "on");

## Test output
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [2.1 3.4 1.9 5.6 4.2 3.3 2.8 6.1 4.9 3.7 2.2 5.1 4.4 3.9 2.6]';
%!   y = [1.2 2.4 3.1 2.6 4.5 3.3 5.1 2.8 4.0 3.6 1.9 4.4 3.2 2.1 5.0]';
%!   h = scatterhist (x, y);
%!   assert (numel (h), 3);
%!   assert (all (isaxes (h)));
%!   ## the scatter axes hold the data
%!   sc = get (h(1), "children");
%!   assert (get (sc(1), "xdata")(:), x, 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # grouped scatterhist runs and returns three axes
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6]';
%!   y = [2 1 4 3 6 5]';
%!   g = [1 1 1 2 2 2]';
%!   h = scatterhist (x, y, "Group", g);
%!   assert (numel (h), 3);
%!   assert (all (isaxes (h)));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # kernel option runs
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = randn (50, 1);
%!   y = randn (50, 1);
%!   h = scatterhist (x, y, "Kernel", "on");
%!   assert (numel (h), 3);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # NBins accepts a two-element specification
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = randn (40, 1);
%!   y = randn (40, 1);
%!   h = scatterhist (x, y, "NBins", [5 8]);
%!   assert (numel (h), 3);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error <Invalid call to scatterhist> scatterhist (1)
%!error <scatterhist: X and Y must be real numeric vectors.> ...
%! scatterhist ({1}, {2})
%!error <scatterhist: X and Y must have the same length.> ...
%! scatterhist ([1 2 3], [1 2])
%!error <scatterhist: name/value arguments must come in pairs.> ...
%! scatterhist ([1 2 3], [1 2 3], "Group")
%!error <scatterhist: unknown property 'bogus'.> ...
%! scatterhist ([1 2 3], [1 2 3], "bogus", 1)
%!error <scatterhist: GROUP must have one entry per point.> ...
%! scatterhist ([1 2 3], [1 2 3], "Group", [1 2])
