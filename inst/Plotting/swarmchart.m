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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {} swarmchart (@var{x}, @var{y})
## @deftypefnx {statistics} {} swarmchart (@var{x}, @var{y}, @var{sz})
## @deftypefnx {statistics} {} swarmchart (@var{x}, @var{y}, @var{sz}, @var{c})
## @deftypefnx {statistics} {} swarmchart (@dots{}, @qcode{'filled'})
## @deftypefnx {statistics} {} swarmchart (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {} swarmchart (@var{ax}, @dots{})
## @deftypefnx {statistics} {@var{s} =} swarmchart (@dots{})
##
## Draw a swarm chart.
##
## @code{swarmchart (@var{x}, @var{y})} draws the points of @var{y} against
## @var{x}, spread sideways where they would otherwise sit on top of one
## another, so that the shape of each column shows where its values gather.
## @var{sz} sets the marker area and @var{c} its colour, as they do for
## @code{scatter}, and @qcode{'filled'} fills the markers.
##
## @var{s} is the scatter object drawn, and carries @code{XJitter},
## @code{YJitter}, @code{XJitterWidth} and @code{YJitterWidth} besides
## everything a scatter object carries.  Setting any of them spreads the
## points again and redraws.
##
## @multitable @columnfractions 0.24 0.02 0.74
## @headitem @var{Name} @tab @tab @var{Value}
##
## @item @qcode{'XJitter'} @tab @tab How the points are spread sideways:
## @qcode{'density'}, the default here, spreads each point by a uniform draw
## weighted by how many of its neighbours share its value, so a crowded part
## of a column spreads wide and a lone point hardly moves;
## @qcode{'rand'} spreads every point alike; @qcode{'randn'} draws the
## offsets from a normal distribution; @qcode{'none'} leaves them in a line.
##
## @item @qcode{'YJitter'} @tab @tab The same, upwards.  The default is
## @qcode{'none'}.
##
## @item @qcode{'XJitterWidth'} @tab @tab How far a point may be moved
## sideways, in the units of the axis.  The default is nine tenths of the
## smallest gap between distinct values of @var{x}.
##
## @item @qcode{'YJitterWidth'} @tab @tab The same, upwards.
## @end multitable
##
## MATLAB leaves @code{XData} as it was given and spreads the points as it
## draws them.  Octave's scatter object has no such step, so this moves the
## values themselves: @code{XData} holds where each point was drawn rather
## than where it came from.  The difference shows in @code{XData} and
## nowhere else, the chart being the same.
##
## @seealso{scatter, boxchart, gscatter}
## @end deftypefn

function s = swarmchart (varargin)

  if (nargin < 2)
    print_usage ();
  endif

  ax = [];
  if (isscalar (varargin{1}) && ishghandle (varargin{1})
      && strcmp (get (varargin{1}, 'type'), 'axes'))
    ax = varargin{1};
    varargin(1) = [];
  endif
  if (numel (varargin) < 2)
    print_usage ();
  endif

  x = varargin{1};
  y = varargin{2};
  rest = varargin(3:end);
  if (! (isnumeric (x) && isreal (x) && isnumeric (y) && isreal (y)))
    error ("swarmchart: X and Y must be real numeric.");
  endif
  x = double (x(:));
  y = double (y(:));
  if (numel (x) != numel (y))
    error ("swarmchart: X and Y must hold the same number of points.");
  endif

  ## Whatever scatter itself takes comes first, then the pairs of our own
  [opts, rest] = scJitterArgs (rest, x, y);

  if (isempty (ax))
    ax = gca ();
  endif

  xd = scSpread (x, y, opts.XJitter, opts.XJitterWidth);
  yd = scSpread (y, x, opts.YJitter, opts.YJitterWidth);
  s = scatter (ax, xd, yd, rest{:});

  ## The jitter is ours to carry: Octave's scatter has no such properties,
  ## and a listener on each lets the points be spread again after the fact
  addproperty ('XJitter', s, 'string', opts.XJitter);
  addproperty ('YJitter', s, 'string', opts.YJitter);
  addproperty ('XJitterWidth', s, 'double', opts.XJitterWidth);
  addproperty ('YJitterWidth', s, 'double', opts.YJitterWidth);
  for name = {'XJitter', 'YJitter', 'XJitterWidth', 'YJitterWidth'}
    addlistener (s, name{1}, @(h, ~) scRespread (h, x, y));
  endfor

endfunction

## Spread the points again from where they started.
function scRespread (h, x, y)

  set (h, 'XData', scSpread (x, y, get (h, 'XJitter'), ...
                             get (h, 'XJitterWidth')));
  set (h, 'YData', scSpread (y, x, get (h, 'YJitter'), ...
                             get (h, 'YJitterWidth')));

endfunction

## Move each value off its column by as much as the kind of jitter allows.
function out = scSpread (v, other, kind, width)

  out = v;
  if (strcmp (kind, 'none') || width <= 0)
    return;
  endif

  ## Every point sharing a column is spread against the others in it
  col = unique (v(! isnan (v)));
  for k = 1:numel (col)
    take = (v == col(k));
    n = sum (take);
    if (n < 2)
      continue;
    endif
    switch (kind)
      case 'rand'
        off = (rand (n, 1) - 0.5) * width;
      case 'randn'
        off = randn (n, 1) * width / 6;
      otherwise
        w = scDensity (other(take));
        off = (rand (n, 1) - 0.5) * width .* w;
    endswitch
    out(take) = v(take) + off;
  endfor

endfunction

## How crowded each value is among its neighbours, from nothing to one.
## A Gaussian kernel with Silverman's bandwidth, which is what makes a dense
## part of a column spread wide and a lone point stay where it is.
function w = scDensity (v)

  v = v(:);
  n = numel (v);
  s = std (v);
  if (! isfinite (s) || s == 0)
    w = ones (n, 1);
    return;
  endif
  h = 1.06 * s * n^(-1/5);
  if (h <= 0)
    w = ones (n, 1);
    return;
  endif
  d = (v - v') / h;
  w = sum (exp (-0.5 * d.^2), 2);
  w = w - min (w);
  if (max (w) > 0)
    w = w / max (w);
  else
    w = ones (n, 1);
  endif

endfunction

## The pairs that belong to the jitter, taken out of what scatter is given.
function [opts, rest] = scJitterArgs (in, x, y)

  opts = struct ('XJitter', 'density', 'YJitter', 'none', ...
                 'XJitterWidth', scDefaultWidth (x), ...
                 'YJitterWidth', scDefaultWidth (y));
  keep = true (1, numel (in));
  k = 1;
  while (k <= numel (in))
    name = in{k};
    if (ischar (name) && any (strcmpi (name, {'XJitter', 'YJitter', ...
                                              'XJitterWidth', ...
                                              'YJitterWidth'})))
      if (k == numel (in))
        error ("swarmchart: '%s' needs a value.", name);
      endif
      opts = scSetOne (opts, name, in{k+1});
      keep(k:k+1) = false;
      k += 2;
    else
      k++;
    endif
  endwhile
  rest = in(keep);

endfunction

## One jitter option, checked.
function opts = scSetOne (opts, name, value)

  if (! isempty (strfind (lower (name), 'width')))
    if (! (isnumeric (value) && isscalar (value) && isreal (value)
           && value >= 0))
      error ("swarmchart: '%s' must be a nonnegative scalar.", name);
    endif
    value = double (value);
  else
    allowed = {'none', 'density', 'rand', 'randn'};
    if (! (ischar (value) && isrow (value) && any (strcmpi (allowed, value))))
      error (strcat ("swarmchart: '%s' must be 'none', 'density',", ...
                     " 'rand' or 'randn'."), name);
    endif
    value = lower (value);
  endif
  fields = {'XJitter', 'YJitter', 'XJitterWidth', 'YJitterWidth'};
  j = find (strcmpi (fields, name), 1);
  opts.(fields{j}) = value;

endfunction

## Nine tenths of the smallest gap between distinct values, as MATLAB sets it.
function w = scDefaultWidth (v)

  u = unique (v(! isnan (v)));
  if (numel (u) < 2)
    w = 0.9;
  else
    w = 0.9 * min (diff (u));
  endif

endfunction

%!demo
%! ## A swarm chart spreads the points of each column sideways, so the shape
%! ## of a column shows where its values gather.  A crowded part spreads
%! ## wide and a lone point hardly moves.
%! load fisheriris
%! g = grp2idx (species);
%! swarmchart (g, meas(:,1));
%! xlabel ('species');
%! ylabel ('sepal length');

%!demo
%! ## A swarm sits over a box chart readily, the one showing the summary and
%! ## the other every observation behind it.
%! load fisheriris
%! g = grp2idx (species);
%! boxchart (g, meas(:,1));
%! hold on
%! swarmchart (g, meas(:,1), 10, [0.3, 0.3, 0.3]);
%! hold off
%! xlabel ('species');
%! ylabel ('sepal length');

%!demo
%! ## How far the points are spread, and by what rule, may be set afterwards.
%! load fisheriris
%! g = grp2idx (species);
%! s = swarmchart (g, meas(:,1));
%! set (s, 'XJitterWidth', 0.3);

%!shared scX, scY
%! scX = [ones(8, 1); 2 * ones(8, 1)];
%! scY = [1; 1; 1; 2; 2; 3; 4; 9; 5; 5; 5; 6; 6; 7; 8; 20];

## MATLAB parity: a scatter object carrying the jitter properties
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (scX, scY);
%!   assert_equal (size (get (s, 'cdata')), [1, 3]);
%!   ## Octave's scatter builds an hggroup under the gnuplot toolkit
%!   if (! strcmp (graphics_toolkit (), 'gnuplot'))
%!     assert_equal (get (s, 'type'), 'scatter');
%!   endif
%!   assert_equal (get (s, 'XJitter'), 'density');
%!   assert_equal (get (s, 'YJitter'), 'none');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # MATLAB parity: the default width is nine tenths of the smallest gap
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (scX, scY);
%!   assert_equal (get (s, 'XJitterWidth'), 0.9, 1e-12);
%!   s2 = swarmchart ([1; 1; 3; 3], [1; 2; 3; 4]);
%!   assert_equal (get (s2, 'XJitterWidth'), 1.8, 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a column with one value alone falls back to the fixed width
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (ones (5, 1), (1:5)');
%!   assert_equal (get (s, 'XJitterWidth'), 0.9, 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Every point stays in its own band, whatever the draw
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (scX, scY);
%!   xd = get (s, 'XData')(:);
%!   w = get (s, 'XJitterWidth');
%!   assert_equal (numel (xd), numel (scX));
%!   assert_equal (all (abs (xd - scX) <= w / 2 + 1e-12), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the other axis is left alone where it is not jittered
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (scX, scY);
%!   assert_equal (get (s, 'YData')(:), scY);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the points really are spread rather than left in two lines
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (scX, scY);
%!   xd = get (s, 'XData')(:);
%!   assert_equal (numel (unique (xd)) > 2, true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # asking for no jitter leaves the values exactly as they came
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (scX, scY, 'XJitter', 'none');
%!   assert_equal (get (s, 'XData')(:), scX);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # setting the width spreads them again, into the narrower band
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (scX, scY);
%!   set (s, 'XJitterWidth', 0.2);
%!   xd = get (s, 'XData')(:);
%!   assert_equal (all (abs (xd - scX) <= 0.1 + 1e-12), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # and setting the kind to none puts them back where they started
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (scX, scY);
%!   set (s, 'XJitter', 'none');
%!   assert_equal (get (s, 'XData')(:), scX);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the y values may be spread instead, and then they move
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (scX, scY, 'XJitter', 'none', 'YJitter', 'rand', ...
%!                   'YJitterWidth', 0.5);
%!   assert_equal (get (s, 'XData')(:), scX);
%!   assert_equal (isequal (get (s, 'YData')(:), scY), false);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the marker size and colour reach scatter as they would on their own
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = swarmchart (scX, scY, 25, [1, 0, 0]);
%!   assert_equal (get (s, 'SizeData'), 25);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the axes to draw into may be given first, as for any plot
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   a1 = subplot (1, 2, 1);
%!   a2 = subplot (1, 2, 2);
%!   s = swarmchart (a2, scX, scY);
%!   assert_equal (get (s, 'parent'), a2);
%!   assert_equal (isempty (get (a1, 'children')), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Input validation
%!error<Invalid call> swarmchart (1)

%!error<swarmchart: X and Y must hold the same number of points.> ...
%! swarmchart ([1; 2], [1; 2; 3])

%!error<swarmchart: X and Y must be real numeric.> ...
%! swarmchart ({1, 2}, [1; 2])

%!error<swarmchart: 'XJitter' must be 'none', 'density', 'rand' or 'randn'.> ...
%! swarmchart ([1; 2], [1; 2], 'XJitter', 'sideways')

%!error<swarmchart: 'XJitterWidth' must be a nonnegative scalar.> ...
%! swarmchart ([1; 2], [1; 2], 'XJitterWidth', -1)

%!error<swarmchart: 'XJitter' needs a value.> ...
%! swarmchart ([1; 2], [1; 2], 'XJitter')
