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
## @deftypefn  {statistics} {} andrewsplot (@var{x})
## @deftypefnx {statistics} {} andrewsplot (@var{x}, @var{name}, @var{value}, @dots{})
## @deftypefnx {statistics} {} andrewsplot (@var{ax}, @dots{})
## @deftypefnx {statistics} {@var{h} =} andrewsplot (@dots{})
##
## Create an Andrews plot of the multivariate data in @var{x}.
##
## @code{andrewsplot (@var{x})} plots each observation (row) of the
## @code{n}-by-@code{p} matrix @var{x} as a smooth curve defined by the finite
## Fourier series
##
## @example
## f_i(t) = x_i1/sqrt(2) + x_i2 sin(2*pi*t) + x_i3 cos(2*pi*t)
##          + x_i4 sin(4*pi*t) + x_i5 cos(4*pi*t) + @dots{}
## @end example
##
## @noindent
## evaluated over @math{t} in the interval @math{[0,1]}, where @var{x_ij} is the
## @math{j}-th variable of the @math{i}-th observation.
##
## The following name/value pairs are accepted:
##
## @table @asis
## @item @qcode{"Group"}
## A grouping variable (numeric, logical, character, string, or cell array of
## strings) with one entry per row of @var{x}.  Curves are colored by group.
##
## @item @qcode{"Standardize"}
## Controls how the columns of @var{x} are transformed before the curves are
## computed: @qcode{"off"} (default) uses the raw data, @qcode{"on"} centers and
## scales each column to zero mean and unit standard deviation, @qcode{"PCA"}
## uses the principal component scores, and @qcode{"PCAStd"} uses the principal
## component scores of the standardized data.
##
## @item @qcode{"Quantile"}
## A scalar @var{alpha} in the interval @math{(0,1)}.  Instead of one curve per
## observation, only three curves per group are drawn: the pointwise median and
## the @var{alpha} and @math{1-}@var{alpha} quantiles of the group's curves.
## @end table
##
## @code{andrewsplot (@var{ax}, @dots{})} plots into the axes @var{ax}.
##
## The optional output @var{h} is a vector of handles to the plotted lines: one
## per observation, or three per group when @qcode{"Quantile"} is used.
##
## @seealso{parallelcoords, glyphplot, pca}
## @end deftypefn

function h = andrewsplot (varargin)

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
    error ("andrewsplot: X must be a real numeric matrix.");
  endif
  n = rows (X);

  ## Parse name/value options
  group = [];
  standardize = "off";
  alpha = [];
  if (mod (numel (varargin), 2) != 0)
    error ("andrewsplot: name/value arguments must come in pairs.");
  endif
  for i = 1:2:numel (varargin)
    name = varargin{i};
    value = varargin{i+1};
    if (! ischar (name))
      error ("andrewsplot: property names must be strings.");
    endif
    switch (lower (name))
      case "group"
        group = value;
      case "standardize"
        standardize = value;
      case "quantile"
        alpha = value;
      otherwise
        error ("andrewsplot: unknown property '%s'.", name);
    endswitch
  endfor

  ## Validate before plotting, so a bad option never leaves a stray figure
  if (! isempty (alpha))
    if (! (isscalar (alpha) && isreal (alpha) && alpha > 0 && alpha < 1))
      error ("andrewsplot: Quantile ALPHA must be a scalar in (0,1).");
    endif
  endif

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
      error ("andrewsplot: invalid Standardize option '%s'.", standardize);
  endswitch

  ## Grouping
  if (isempty (group))
    gidx = ones (n, 1);
    gnames = {"1"};
  else
    [gidx, gnames] = grp2idx (group);
    if (numel (gidx) != n)
      error ("andrewsplot: GROUP must have one entry per row of X.");
    endif
  endif
  k = numel (gnames);
  gcol = lines (k);

  ## Evaluate the Andrews curves f_i(t) over t in [0,1]
  t = linspace (0, 1, 1001);
  p = columns (Z);
  F = (Z(:,1) / sqrt (2)) * ones (1, numel (t));
  for j = 2:p
    kf = floor (j / 2);
    if (mod (j, 2) == 0)
      F += Z(:,j) * sin (2 * pi * kf * t);
    else
      F += Z(:,j) * cos (2 * pi * kf * t);
    endif
  endfor

  if (isempty (hax))
    hax = newplot ();
  else
    newplot (hax);
  endif
  old_hold = ishold (hax);
  hold (hax, "on");

  h = [];
  if (isempty (alpha))
    ## One curve per observation, colored by group
    for i = 1:n
      h(end+1) = line (hax, t, F(i,:), "color", gcol(gidx(i),:));
    endfor
  else
    ## Median and alpha / 1-alpha quantile curves per group
    for g = 1:k
      Fg = F(gidx == g, :);
      med = median (Fg, 1);
      lo = quantile (Fg, alpha, 1);
      hi = quantile (Fg, 1 - alpha, 1);
      h(end+1) = line (hax, t, med, "color", gcol(g,:), "linewidth", 2);
      h(end+1) = line (hax, t, lo, "color", gcol(g,:), "linestyle", "--");
      h(end+1) = line (hax, t, hi, "color", gcol(g,:), "linestyle", "--");
    endfor
  endif

  xlabel (hax, "t");
  ylabel (hax, "f(t)");

  if (k > 1 && ! isempty (group))
    warning ("off", "Octave:legend:unimplemented-location", "local");
    ## One representative line handle per group for the legend
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
%! ## Andrews plot of Fisher's iris data, grouped by species.
%!
%! load fisheriris;
%! andrewsplot (meas, "Group", species);

%!demo
%! ## The same data with median and quartile curves per species.
%!
%! load fisheriris;
%! andrewsplot (meas, "Group", species, "Quantile", 0.25);

## Test output
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [1 2 3 4; 5 6 7 8; 2 3 1 4];
%!   h = andrewsplot (X);
%!   assert (numel (h), 3);
%!   t = get (h(1), "xdata");
%!   assert (numel (t), 1001);
%!   assert (t([1, 501, 1001]), [0, 0.5, 1], 1e-12);
%!   y1 = get (h(1), "ydata");
%!   assert (y1([1, 501, 1001]), [3.70711, -2.29289, 3.70711], 1e-4);
%!   y3 = get (h(3), "ydata");
%!   assert (y3([1, 501, 1001]), [2.41421, 0.41421, 2.41421], 1e-4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # Standardize "on" (z-score each column)
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [1 2 3 4; 5 6 7 8; 2 3 1 4];
%!   h = andrewsplot (X, "Standardize", "on");
%!   y1 = get (h(1), "ydata");
%!   assert (y1([1, 501, 1001]), [-0.78436, -0.34792, -0.78436], 1e-4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # Standardize "pca"
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [1 2 3 4; 5 6 7 8; 2 3 1 4];
%!   h = andrewsplot (X, "Standardize", "pca");
%!   y1 = get (h(1), "ydata");
%!   assert (y1([1, 501, 1001]), [-1.76701, -1.76701, -1.76701], 1e-4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # grouping and quantile mode give three curves per group
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [1 2 3 4; 5 6 7 8; 2 3 1 4; 3 1 4 1; 5 9 2 6; 4 2 1 3];
%!   g = [1 1 1 2 2 2]';
%!   h = andrewsplot (X, "Group", g, "Quantile", 0.25);
%!   assert (numel (h), 6);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error <Invalid call to andrewsplot> andrewsplot ()
%!error <andrewsplot: X must be a real numeric matrix.> andrewsplot ({1})
%!error <andrewsplot: name/value arguments must come in pairs.> ...
%! andrewsplot (ones (3, 2), "Group")
%!error <andrewsplot: unknown property 'bogus'.> ...
%! andrewsplot (ones (3, 2), "bogus", 1)
%!error <andrewsplot: invalid Standardize option 'xxx'.> ...
%! andrewsplot (ones (3, 2), "Standardize", "xxx")
%!error <andrewsplot: GROUP must have one entry per row of X.> ...
%! andrewsplot (ones (3, 2), "Group", [1 2])
%!error <andrewsplot: Quantile ALPHA must be a scalar in .0,1..> ...
%! andrewsplot (ones (3, 2), "Quantile", 1.5)
