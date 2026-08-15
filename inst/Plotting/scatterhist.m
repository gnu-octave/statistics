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
##
## @qcode{'Parent'} draws into a supplied figure or uipanel instead of the
## current figure.  The container is used as it stands and is never
## cleared, so a scatterhist can be placed alongside other axes.
##
## @qcode{'Location'} names the corner the scatter occupies,
## @qcode{'SouthWest'} by default, or @qcode{'SouthEast'}, @qcode{'NorthEast'}
## or @qcode{'NorthWest'}; the marginals take the two opposite sides.
##
## @qcode{'Direction'} points the marginal bars toward the scatter plot,
## @qcode{'in'} by default, or away from it with @qcode{'out'}, whichever side
## @qcode{'Location'} has placed them on.
##
## @qcode{'PlotGroup'} draws one marginal per group with @qcode{'on'}, or a
## single pooled marginal with @qcode{'off'}.  It defaults to @qcode{'on'} when
## @qcode{'Group'} is given and to @qcode{'off'} otherwise.
##
## @qcode{'Style'} outlines the marginals as @qcode{'bar'} or @qcode{'stairs'},
## defaulting to stairs when the data are grouped and bars when they are not.
##
## @qcode{'Color'} sets the group colours, either as a character vector of
## colour names or as an @math{N}-by-3 matrix of RGB values, cycled over the
## groups.
##
## @qcode{'LineStyle'} and @qcode{'LineWidth'} style the marginal outlines, not
## the scatter markers, and are cycled over the groups.
##
## @qcode{'Bandwidth'} sets the kernel bandwidth used when @qcode{'Kernel'} is
## on: a scalar for all marginals, a pair for @var{x} and @var{y}, or one row
## per group.
##
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

  ## Parse name/value options through the datatypes parser.  It consumes the
  ## names it knows and returns whatever is left, so anything remaining is an
  ## unrecognised property.
  optNames = {"Group", "NBins", "Kernel", "Location", "Legend", "Marker", ...
              "MarkerSize", "Direction", "PlotGroup", "Style", "Color", ...
              "Bandwidth", "LineStyle", "LineWidth", "Parent"};
  dfValues = {[], [], "off", "southwest", "", "o", 6, "in", "", "", [], ...
              [], "-", 0.5, []};
  if (mod (numel (varargin), 2) != 0)
    error ("scatterhist: name/value arguments must come in pairs.");
  endif
  for i = 1:2:numel (varargin)
    if (! ischar (varargin{i}))
      error ("scatterhist: property names must be strings.");
    endif
  endfor
  [group, nbins, kernel, location, legend_opt, marker, markersize, ...
   direction, plotgroup, style, colour, bandwidth, linestyle, linewidth, ...
   parent, rest] = parsePairedArguments (optNames, dfValues, varargin(:));
  if (! isempty (rest))
    error ("scatterhist: unknown property '%s'.", rest{1});
  endif

  ## The parser only splits the pairs; every value is validated here.
  if (! isempty (nbins))
    if (! (isnumeric (nbins) && all (nbins(:) > 0)
           && all (nbins(:) == fix (nbins(:))) && numel (nbins) <= 2))
      error (strcat ("scatterhist: NBINS must be one or two positive", ...
                     " integers."));
    endif
  endif
  if (! (ischar (kernel) && any (strcmpi (kernel, {"off", "on", "overlay"}))))
    error ("scatterhist: KERNEL must be 'off', 'on' or 'overlay'.");
  endif
  if (! (ischar (location) && any (strcmpi (location, {"southwest", ...
         "southeast", "northeast", "northwest"}))))
    error (strcat ("scatterhist: LOCATION must be 'SouthWest', 'SouthEast',", ...
                   " 'NorthEast' or 'NorthWest'."));
  endif
  if (! (ischar (legend_opt) && any (strcmpi (legend_opt, {"", "on", "off"}))))
    error ("scatterhist: LEGEND must be 'on' or 'off'.");
  endif
  if (! ischar (marker))
    error ("scatterhist: MARKER must be a character vector of markers.");
  endif
  if (! (isnumeric (markersize) && all (markersize(:) > 0)))
    error ("scatterhist: MARKERSIZE must be positive numeric values.");
  endif
  if (! (ischar (direction) && any (strcmpi (direction, {"in", "out"}))))
    error ("scatterhist: DIRECTION must be 'in' or 'out'.");
  endif
  if (! isempty (style)
      && ! (ischar (style) && any (strcmpi (style, {"bar", "stairs"}))))
    error ("scatterhist: STYLE must be 'bar' or 'stairs'.");
  endif
  if (! isempty (plotgroup)
      && ! (ischar (plotgroup) && any (strcmpi (plotgroup, {"on", "off"}))))
    error ("scatterhist: PLOTGROUP must be 'on' or 'off'.");
  endif
  if (! isempty (parent))
    if (! (isscalar (parent) && ishghandle (parent)
           && any (strcmp (get (parent, "type"), {"figure", "uipanel"}))))
      error ("scatterhist: PARENT must be a figure or uipanel handle.");
    endif
  endif

  ## Marginal line style.  'Bandwidth' is a scalar, a pair for x and y, or one
  ## row per group; 'LineStyle' and 'LineWidth' cycle over the groups and apply
  ## to the marginals, not to the scatter markers.
  if (! isempty (bandwidth))
    if (! (isnumeric (bandwidth) && all (bandwidth(:) > 0)))
      error ("scatterhist: BANDWIDTH must be positive numeric values.");
    endif
    if (isvector (bandwidth))
      bandwidth = bandwidth(:)';
    endif
  endif
  if (! (ischar (linestyle) || iscellstr (linestyle)))
    error ("scatterhist: LINESTYLE must be a string or a cell array of them.");
  endif
  if (ischar (linestyle))
    linestyle = {linestyle};
  endif
  if (! (isnumeric (linewidth) && all (linewidth(:) > 0)))
    error ("scatterhist: LINEWIDTH must be positive numeric values.");
  endif

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
  if (isempty (colour))
    gcol = lines (k);
  else
    if (ischar (colour))
      colour = colour(:);
      gcol = zeros (numel (colour), 3);
      for i = 1:numel (colour)
        gcol(i,:) = colour_rgb (colour(i));
      endfor
    elseif (isnumeric (colour) && columns (colour) == 3)
      gcol = colour;
    else
      error (strcat ("scatterhist: COLOR must be a character vector of", ...
                     " colour names or an N-by-3 matrix of RGB values."));
    endif
    ## cycle the supplied colours over the groups
    gcol = gcol(mod ((0:k-1), rows (gcol)) + 1, :);
  endif
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

  ## Axes layout.  'Location' names the corner the scatter sits in, so the
  ## marginals take the two opposite sides: with the default "southwest" the
  ## scatter is lower left, the x marginal above it and the y marginal to its
  ## right.
  west = any (strcmpi (location, {"southwest", "northwest"}));
  south = any (strcmpi (location, {"southwest", "southeast"}));
  ## A container supplied through 'Parent' is drawn into as it stands: it
  ## belongs to the caller, so it is never cleared.  Only a figure of our own
  ## choosing is.
  if (isempty (parent))
    cf = gcf ();
    clf (cf);
  else
    cf = parent;
  endif
  if (west)
    xs = 0.10;  ys_marg = 0.72;
  else
    xs = 0.30;  ys_marg = 0.08;
  endif
  if (south)
    ys = 0.10;  xs_marg = 0.72;
  else
    ys = 0.30;  xs_marg = 0.08;
  endif
  pos_s = [xs, ys, 0.60, 0.60];
  pos_x = [xs, xs_marg, 0.60, 0.20];
  pos_y = [ys_marg, ys, 0.20, 0.60];
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

  ## Marginal for x along the top, for y along the right.  'Style' picks the
  ## outline, defaulting to stairs when grouped and bars when not, and
  ## 'PlotGroup' says whether the groups are drawn apart or pooled.
  if (isempty (style))
    multi = (k > 1);
  else
    multi = strcmpi (style, "stairs");
  endif
  if (isempty (plotgroup))
    bygroup = ! isempty (group);
  else
    bygroup = strcmpi (plotgroup, "on");
  endif
  if (bygroup)
    for g = 1:k
      xg = x(gidx == g & ! isnan (x));
      yg = y(gidx == g & ! isnan (y));
      marginal (ax_x, xg, nbx, gcol(g,:), do_kernel, multi, false, ...
                sty (bandwidth, linestyle, linewidth, g, 1));
      marginal (ax_y, yg, nby, gcol(g,:), do_kernel, multi, true, ...
                sty (bandwidth, linestyle, linewidth, g, 2));
    endfor
  else
    marginal (ax_x, x(! isnan (x)), nbx, gcol(1,:), do_kernel, multi, false, ...
              sty (bandwidth, linestyle, linewidth, 1, 1));
    marginal (ax_y, y(! isnan (y)), nby, gcol(1,:), do_kernel, multi, true, ...
              sty (bandwidth, linestyle, linewidth, 1, 2));
  endif

  ## 'Direction' points the marginal bars toward the scatter ("in", as MATLAB
  ## defaults) or away from it ("out").  Which way that is depends on where
  ## 'Location' put the marginals: bars grow along the positive axis, so the x
  ## marginal must be reversed only when it sits above the scatter, and the y
  ## marginal only when it sits to the scatter's right.
  toward = strcmpi (direction, "in");
  set (ax_x, "ydir", ternary (xor (toward, ! south), "reverse", "normal"));
  set (ax_y, "xdir", ternary (xor (toward, ! west), "reverse", "normal"));

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

## Choose between two values without an if block.
function v = ternary (cond, a, b)
  if (cond)
    v = a;
  else
    v = b;
  endif
endfunction

## Pick the bandwidth, line style and line width for group G and variable V
## (1 for x, 2 for y), cycling whatever the user supplied.
function sy = sty (bandwidth, linestyle, linewidth, g, v)
  if (isempty (bandwidth))
    sy.bw = [];
  elseif (isscalar (bandwidth))
    sy.bw = bandwidth;
  elseif (isrow (bandwidth))
    sy.bw = bandwidth(mod (v - 1, numel (bandwidth)) + 1);
  else
    row = bandwidth(mod (g - 1, rows (bandwidth)) + 1, :);
    sy.bw = row(mod (v - 1, numel (row)) + 1);
  endif
  sy.ls = linestyle{mod (g - 1, numel (linestyle)) + 1};
  sy.lw = linewidth(mod (g - 1, numel (linewidth)) + 1);
endfunction

## Map a colour name to its RGB triplet.
function rgb = colour_rgb (c)
  switch (lower (c))
    case "r", rgb = [1, 0, 0];
    case "g", rgb = [0, 1, 0];
    case "b", rgb = [0, 0, 1];
    case "c", rgb = [0, 1, 1];
    case "m", rgb = [1, 0, 1];
    case "y", rgb = [1, 1, 0];
    case "k", rgb = [0, 0, 0];
    case "w", rgb = [1, 1, 1];
    otherwise
      error ("scatterhist: unknown colour '%s'.", c);
  endswitch
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
function marginal (ax, v, nb, col, do_kernel, multi, horiz, sy)
  if (isempty (v))
    return;
  endif
  if (do_kernel)
    if (isempty (sy.bw))
      [f, u] = ksdensity (v);
    else
      [f, u] = ksdensity (v, "Bandwidth", sy.bw);
    endif
    if (horiz)
      line (ax, f, u, "color", col, "linestyle", sy.ls, "linewidth", sy.lw);
    else
      line (ax, u, f, "color", col, "linestyle", sy.ls, "linewidth", sy.lw);
    endif
  else
    [nn, cc] = hist (v, nb);
    if (multi)
      ## stairstep outline
      e = cc(1:end-1) + diff (cc) / 2;
      ce = [cc(1)-(cc(2)-cc(1))/2, e, cc(end)+(cc(end)-cc(end-1))/2];
      ne = [nn, nn(end)];
      if (horiz)
        stairs (ax, ne, ce, "color", col, "linestyle", sy.ls, ...
                "linewidth", sy.lw);
      else
        stairs (ax, ce, ne, "color", col, "linestyle", sy.ls, ...
                "linewidth", sy.lw);
      endif
    else
      if (horiz)
        barh (ax, cc, nn, 1.0, "facecolor", col, "edgecolor", col, ...
              "linestyle", sy.ls, "linewidth", sy.lw);
      else
        bar (ax, cc, nn, 1.0, "facecolor", col, "edgecolor", col, ...
             "linestyle", sy.ls, "linewidth", sy.lw);
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
%!test  # LineWidth and LineStyle apply to the marginals and cycle over groups
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6 7 8]';  y = [2 1 4 3 6 5 8 7]';
%!   g = [1 1 1 1 2 2 2 2]';
%!   h = scatterhist (x, y, "Group", g, "LineWidth", 3);
%!   ch = get (h(2), "children");
%!   assert_equal (get (ch(1), "linewidth"), 3);
%!   assert_equal (get (ch(2), "linewidth"), 3);
%!   h = scatterhist (x, y, "Group", g, "LineStyle", {"--", ":"});
%!   ch = get (h(2), "children");
%!   assert_equal (get (ch(2), "linestyle"), "--");
%!   assert_equal (get (ch(1), "linestyle"), ":");
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # Bandwidth is passed to the kernel density
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6 7 8]';  y = [2 1 4 3 6 5 8 7]';
%!   h = scatterhist (x, y, "Kernel", "on", "Bandwidth", 2);
%!   wide = get (get (h(2), "children")(1), "ydata");
%!   h = scatterhist (x, y, "Kernel", "on");
%!   auto = get (get (h(2), "children")(1), "ydata");
%!   assert_equal (! isequal (wide, auto), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!error <scatterhist: LINEWIDTH must be positive numeric values.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "LineWidth", -1)
%!error <scatterhist: LINESTYLE must be a string or a cell array of them.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "LineStyle", 5)
%!error <scatterhist: BANDWIDTH must be positive numeric values.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Bandwidth", 0)

%!test  # Parent draws into a supplied container without clearing it
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6 7 8]';  y = [2 1 4 3 6 5 8 7]';
%!   pre = axes ("parent", hf, "position", [0.01, 0.01, 0.05, 0.05]);
%!   h = scatterhist (x, y, "Parent", hf);
%!   assert_equal (all (arrayfun (@(a) get (a, "parent"), h) == hf, 'all'), ...
%!                 true);
%!   ## the caller's own axes are left alone
%!   assert_equal (ishghandle (pre), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a uipanel is a container too
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6 7 8]';  y = [2 1 4 3 6 5 8 7]';
%!   hp = uipanel ("parent", hf, "position", [0, 0, 1, 1]);
%!   h = scatterhist (x, y, "Parent", hp);
%!   assert_equal (all (arrayfun (@(a) get (a, "parent"), h) == hp, 'all'), ...
%!                 true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # without Parent the current figure is still claimed and cleared
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6 7 8]';  y = [2 1 4 3 6 5 8 7]';
%!   old = axes ("parent", hf);
%!   h = scatterhist (x, y);
%!   assert_equal (ishghandle (old), false);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!error <scatterhist: PARENT must be a figure or uipanel handle.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Parent", 0)
%!error <scatterhist: PARENT must be a figure or uipanel handle.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Parent", -5)

%!test  # Location moves the scatter to the named corner and the marginals with it
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6 7 8]';  y = [2 1 4 3 6 5 8 7]';
%!   h = scatterhist (x, y, "Location", "SouthWest");
%!   ps = get (h(1), "position");  px = get (h(2), "position");
%!   assert_equal (ps(1:2), [0.10, 0.10], 1e-12);
%!   assert_equal (px(2) > ps(2), true);          # x marginal above the scatter
%!   h = scatterhist (x, y, "Location", "NorthEast");
%!   ps = get (h(1), "position");  px = get (h(2), "position");
%!   assert_equal (ps(1:2), [0.30, 0.30], 1e-12);
%!   assert_equal (px(2) < ps(2), true);          # and below it here
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # Direction is relative to where Location put each marginal
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6 7 8]';  y = [2 1 4 3 6 5 8 7]';
%!   ## "in" points the bars at the scatter from whichever side they sit on
%!   h = scatterhist (x, y, "Location", "SouthWest");
%!   assert_equal (get (h(2), "ydir"), "reverse");
%!   assert_equal (get (h(3), "xdir"), "reverse");
%!   h = scatterhist (x, y, "Location", "NorthEast");
%!   assert_equal (get (h(2), "ydir"), "normal");
%!   assert_equal (get (h(3), "xdir"), "normal");
%!   h = scatterhist (x, y, "Location", "SouthEast");
%!   assert_equal (get (h(2), "ydir"), "reverse");
%!   assert_equal (get (h(3), "xdir"), "normal");
%!   h = scatterhist (x, y, "Location", "NorthWest");
%!   assert_equal (get (h(2), "ydir"), "normal");
%!   assert_equal (get (h(3), "xdir"), "reverse");
%!   ## and "out" inverts each of them
%!   h = scatterhist (x, y, "Location", "NorthEast", "Direction", "out");
%!   assert_equal (get (h(2), "ydir"), "reverse");
%!   assert_equal (get (h(3), "xdir"), "reverse");
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # Direction points the marginals toward the scatter by default
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6 7 8]';  y = [2 1 4 3 6 5 8 7]';
%!   h = scatterhist (x, y);
%!   assert_equal (get (h(2), "ydir"), "reverse");
%!   assert_equal (get (h(3), "xdir"), "reverse");
%!   h = scatterhist (x, y, "Direction", "out");
%!   assert_equal (get (h(2), "ydir"), "normal");
%!   assert_equal (get (h(3), "xdir"), "normal");
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # Color sets the group colours, cycling if fewer than the groups
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6 7 8]';  y = [2 1 4 3 6 5 8 7]';
%!   g = [1 1 1 1 2 2 2 2]';
%!   h = scatterhist (x, y, "Group", g, "Color", "rb");
%!   c = get (get (h(1), "children"), "color");
%!   assert_equal (c{2}, [1, 0, 0]);
%!   assert_equal (c{1}, [0, 0, 1]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # PlotGroup pools the marginals when off
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [1 2 3 4 5 6 7 8]';  y = [2 1 4 3 6 5 8 7]';
%!   g = [1 1 1 1 2 2 2 2]';
%!   h = scatterhist (x, y, "Group", g, "PlotGroup", "off");
%!   assert_equal (numel (get (h(2), "children")), 1);
%!   h = scatterhist (x, y, "Group", g, "PlotGroup", "on");
%!   assert_equal (numel (get (h(2), "children")), 2);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!error <scatterhist: NBINS must be one or two positive integers.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "NBins", 0)
%!error <scatterhist: NBINS must be one or two positive integers.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "NBins", 2.5)
%!error <scatterhist: NBINS must be one or two positive integers.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "NBins", [2, 3, 4])
%!error <scatterhist: KERNEL must be 'off', 'on' or 'overlay'.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Kernel", "sometimes")
%!error <scatterhist: KERNEL must be 'off', 'on' or 'overlay'.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Kernel", 1)
%!error <scatterhist: LOCATION must be 'SouthWest', 'SouthEast', 'NorthEast' or 'NorthWest'.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Location", "middle")
%!error <scatterhist: LEGEND must be 'on' or 'off'.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Legend", "maybe")
%!error <scatterhist: MARKER must be a character vector of markers.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Marker", 7)
%!error <scatterhist: MARKERSIZE must be positive numeric values.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "MarkerSize", 0)
%!error <scatterhist: COLOR must be a character vector of colour names or an N-by-3 matrix of RGB values.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Color", [1, 2])
%!error <scatterhist: unknown colour 'q'.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Color", "q")

%!error <scatterhist: DIRECTION must be 'in' or 'out'.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Direction", "sideways")
%!error <scatterhist: STYLE must be 'bar' or 'stairs'.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "Style", "curvy")
%!error <scatterhist: PLOTGROUP must be 'on' or 'off'.> ...
%! scatterhist ([1 2 3]', [1 2 3]', "PlotGroup", "maybe")

%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [2.1 3.4 1.9 5.6 4.2 3.3 2.8 6.1 4.9 3.7 2.2 5.1 4.4 3.9 2.6]';
%!   y = [1.2 2.4 3.1 2.6 4.5 3.3 5.1 2.8 4.0 3.6 1.9 4.4 3.2 2.1 5.0]';
%!   h = scatterhist (x, y);
%!   assert_equal (numel (h), 3);
%!   assert_equal (all (isaxes (h), 'all'), true);
%!   ## the scatter axes hold the data
%!   sc = get (h(1), "children");
%!   assert_equal (get (sc(1), "xdata")(:), x, 1e-12);
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
%!   assert_equal (numel (h), 3);
%!   assert_equal (all (isaxes (h), 'all'), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # kernel option runs
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = randn (50, 1);
%!   y = randn (50, 1);
%!   h = scatterhist (x, y, "Kernel", "on");
%!   assert_equal (numel (h), 3);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # NBins accepts a two-element specification
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = randn (40, 1);
%!   y = randn (40, 1);
%!   h = scatterhist (x, y, "NBins", [5 8]);
%!   assert_equal (numel (h), 3);
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
