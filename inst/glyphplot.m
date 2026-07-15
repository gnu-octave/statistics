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
## @deftypefn  {statistics} {} glyphplot (@var{x})
## @deftypefnx {statistics} {} glyphplot (@var{x}, @var{name}, @var{value}, @dots{})
## @deftypefnx {statistics} {@var{g} =} glyphplot (@dots{})
##
## Create a star (glyph) plot of the multivariate data in @var{x}.
##
## @code{glyphplot (@var{x})} draws each observation (row) of the
## @code{n}-by-@code{p} matrix @var{x} as a star glyph, arranged on a grid.  The
## @var{p} spokes of each star radiate from its center at equally spaced angles,
## with lengths proportional to the values of the @var{p} variables; the tips of
## the spokes are joined to form the star perimeter.
##
## The following name/value pairs are accepted:
##
## @table @asis
## @item @qcode{"Glyph"}
## @qcode{"star"} (default) draws star glyphs.  @qcode{"face"} (Chernoff faces)
## is not currently supported.
##
## @item @qcode{"Standardize"}
## How the columns of @var{x} are scaled to spoke lengths: @qcode{"column"}
## (default) scales each column to the range @math{[0,1]}, @qcode{"matrix"}
## scales the whole matrix to @math{[0,1]}, @qcode{"PCA"} uses principal
## component scores scaled to @math{[0,1]}, and @qcode{"off"} uses the values as
## given.  A spoke of relative length 0 is still drawn at 10% of the maximum
## radius so that it remains visible.
##
## @item @qcode{"Grid"}
## A two-element vector @code{[rows cols]} specifying the layout of the glyphs.
## The default is chosen automatically.
##
## @item @qcode{"Centers"}
## An @code{n}-by-2 matrix giving the center coordinates of the glyphs
## explicitly, overriding @qcode{"Grid"}.
##
## @item @qcode{"Radius"}
## The maximum glyph radius (default 0.4).
##
## @item @qcode{"ObsLabels"}
## A character array or cell array of strings labeling the observations.  The
## default is the observation numbers.
## @end table
##
## The optional output @var{g} is an @code{n}-by-3 matrix of handles whose
## columns hold, respectively, the star perimeters, the star spokes, and the
## text labels.
##
## @seealso{andrewsplot, parallelcoords}
## @end deftypefn

function g = glyphplot (varargin)

  ## Optional leading figure handle
  parent = [];
  if (numel (varargin) > 0 && isscalar (varargin{1}) && ishghandle (varargin{1})
      && strcmp (get (varargin{1}, "type"), "figure"))
    parent = varargin{1};
    varargin(1) = [];
  endif

  if (numel (varargin) < 1)
    print_usage ();
  endif

  X = varargin{1};
  varargin(1) = [];
  if (! isnumeric (X) || ! isreal (X) || ndims (X) > 2)
    error ("glyphplot: X must be a real numeric matrix.");
  endif
  n = rows (X);
  p = columns (X);

  ## Parse name/value options
  glyph = "star";
  standardize = "column";
  grid = [];
  centers = [];
  maxr = 0.4;
  obslabels = {};
  if (mod (numel (varargin), 2) != 0)
    error ("glyphplot: name/value arguments must come in pairs.");
  endif
  for i = 1:2:numel (varargin)
    name = varargin{i};
    value = varargin{i+1};
    if (! ischar (name))
      error ("glyphplot: property names must be strings.");
    endif
    switch (lower (name))
      case "glyph"
        glyph = value;
      case "standardize"
        standardize = value;
      case "grid"
        grid = value;
      case "centers"
        centers = value;
      case "radius"
        maxr = value;
      case "obslabels"
        obslabels = value;
      otherwise
        error ("glyphplot: unknown property '%s'.", name);
    endswitch
  endfor

  if (strcmpi (glyph, "face"))
    error ("glyphplot: 'face' (Chernoff) glyphs are not yet supported.");
  elseif (! strcmpi (glyph, "star"))
    error ("glyphplot: Glyph must be 'star' or 'face'.");
  endif

  ## Scale the columns to spoke lengths in [0,1]
  switch (lower (standardize))
    case "column"
      S = normalize01 (X, 1);
    case "matrix"
      S = normalize01 (X(:), 1);
      S = reshape (S, size (X));
    case "pca"
      [~, sc] = pca (X);
      S = normalize01 (sc, 1);
      p = columns (S);
    case "off"
      S = X;
    otherwise
      error ("glyphplot: invalid Standardize option '%s'.", standardize);
  endswitch

  if (ischar (obslabels))
    obslabels = cellstr (obslabels);
  endif
  if (isempty (obslabels))
    obslabels = arrayfun (@num2str, (1:n)', "UniformOutput", false);
  endif

  ## Glyph centers
  if (! isempty (centers))
    if (! isnumeric (centers) || rows (centers) != n || columns (centers) != 2)
      error ("glyphplot: Centers must be an n-by-2 numeric matrix.");
    endif
    cxy = centers;
    nrows = max (centers(:,2));
  else
    if (isempty (grid))
      ncols = ceil (sqrt (n));
      nrows = ceil (n / ncols);
    else
      nrows = grid(1);
      ncols = grid(2);
    endif
    cxy = zeros (n, 2);
    for i = 1:n
      cxy(i,1) = mod (i - 1, ncols) + 1;
      cxy(i,2) = nrows - floor ((i - 1) / ncols);
    endfor
  endif

  hax = newplot ();
  old_hold = ishold (hax);
  hold (hax, "on");

  ang = 2 * pi * (0:p-1) / p;
  g = zeros (n, 3);
  for i = 1:n
    cx = cxy(i,1);
    cy = cxy(i,2);
    r = maxr * (0.1 + 0.9 * S(i,:));
    tx = cx + r .* cos (ang);
    ty = cy + r .* sin (ang);
    ## Perimeter (closed)
    g(i,1) = line (hax, [tx, tx(1)], [ty, ty(1)], "color", [0 0 1]);
    set (g(i,1), "userdata", [cx, cy, 1]);
    ## Spokes from the center to each tip
    sx = reshape ([cx * ones(1, p); tx; nan(1, p)], 1, []);
    sy = reshape ([cy * ones(1, p); ty; nan(1, p)], 1, []);
    g(i,2) = line (hax, sx, sy, "color", [0 0 1]);
    ## Text label below the glyph
    g(i,3) = text (hax, cx, cy - maxr * 1.2, obslabels{i}, ...
                   "horizontalalignment", "center", "verticalalignment", "top");
  endfor

  axis (hax, "equal");
  axis (hax, "off");

  if (! old_hold)
    hold (hax, "off");
  endif

  if (nargout == 0)
    clear g;
  endif

endfunction

## Scale to [0,1] along dimension dim (columns).  Constant columns map to 0.
function S = normalize01 (X, dim)
  lo = min (X, [], dim);
  hi = max (X, [], dim);
  rng = hi - lo;
  rng(rng == 0) = 1;
  S = (X - lo) ./ rng;
endfunction

%!demo
%! ## Star plot of the first few cars in the carsmall data set.
%!
%! load carsmall;
%! X = [Acceleration, Cylinders, Displacement, Horsepower, Weight];
%! glyphplot (X(1:9,:), "ObsLabels", cellstr (num2str ((1:9)')));

## Test output
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [1 4 2; 3 2 5; 5 5 1; 2 1 4];
%!   g = glyphplot (X);
%!   assert (size (g), [4, 3]);
%!   ## star perimeter of observation 1 (column standardize, radius 0.4)
%!   assert (get (g(1,1), "xdata"), [1.04 0.845 0.935 1.04], 1e-4);
%!   assert (get (g(1,1), "ydata"), [2 2.2685 1.8874 2], 1e-4);
%!   assert (get (g(1,1), "userdata"), [1 2 1], 1e-12);
%!   ## observation 2
%!   assert (get (g(2,1), "xdata"), [2.22 1.935 1.8 2.22], 1e-4);
%!   assert (get (g(2,1), "ydata"), [2 2.1126 1.6536 2], 1e-4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test  # Radius option scales the glyph
%! hf = figure ("visible", "off");
%! unwind_protect
%!   X = [1 4 2; 3 2 5; 5 5 1; 2 1 4];
%!   g = glyphplot (X, "Radius", 0.8);
%!   ## perimeter x of obs 1 spoke 1: cx + 0.8*(0.1+0.9*0) = 1 + 0.08
%!   assert (get (g(1,1), "xdata")(1), 1.08, 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error <Invalid call to glyphplot> glyphplot ()
%!error <glyphplot: X must be a real numeric matrix.> glyphplot ({1})
%!error <glyphplot: 'face' .Chernoff. glyphs are not yet supported.> ...
%! glyphplot (ones (3, 3), "Glyph", "face")
%!error <glyphplot: name/value arguments must come in pairs.> ...
%! glyphplot (ones (3, 3), "Radius")
%!error <glyphplot: unknown property 'bogus'.> ...
%! glyphplot (ones (3, 3), "bogus", 1)
%!error <glyphplot: invalid Standardize option 'xxx'.> ...
%! glyphplot (ones (3, 3), "Standardize", "xxx")
%!error <glyphplot: Centers must be an n-by-2 numeric matrix.> ...
%! glyphplot (ones (3, 3), "Centers", [1 2])
