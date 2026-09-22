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
## @deftypefn  {statistics} {} boxchart (@var{ydata})
## @deftypefnx {statistics} {} boxchart (@var{xgroupdata}, @var{ydata})
## @deftypefnx {statistics} {} boxchart (@dots{}, @qcode{'GroupByColor'}, @var{cgroupdata})
## @deftypefnx {statistics} {} boxchart (@var{tbl}, @var{yvar})
## @deftypefnx {statistics} {} boxchart (@var{tbl}, @var{xvar}, @var{yvar})
## @deftypefnx {statistics} {} boxchart (@dots{}, @var{name}, @var{value})
## @deftypefnx {statistics} {} boxchart (@var{ax}, @dots{})
## @deftypefnx {statistics} {@var{b} =} boxchart (@dots{})
##
## Draw a box chart.
##
## @code{boxchart (@var{ydata})} draws one box over the observations in
## @var{ydata}.  The box spans the lower and upper quartiles with the median
## across it, the whiskers reach the furthest observations within one and a
## half interquartile ranges of the box, and whatever lies beyond is drawn as
## a marker.
##
## @code{boxchart (@var{xgroupdata}, @var{ydata})} draws one box per distinct
## value of @var{xgroupdata}, which holds one value per observation.
##
## @code{boxchart (@var{tbl}, @var{yvar})} and
## @code{boxchart (@var{tbl}, @var{xvar}, @var{yvar})} take the observations,
## and where asked the grouping, from the named columns of a table.  The
## chart keeps the table, so assigning @code{SourceTable} again re-reads it.
##
## @code{@var{b} = boxchart (@dots{})} returns the
## @code{stats.chart.BoxChart} object drawn, whose properties may be set
## afterwards to redraw it; see @code{stats.chart.BoxChart}.  Where
## @qcode{'GroupByColor'} names a grouping, one object is returned per
## colour group, in the order of its distinct values.
##
## @qcode{'GroupByColor'} cannot be given with a table, which MATLAB R2026a
## refuses as well.
##
## Every other name-value pair sets the property of that name on the object;
## see @code{stats.chart.BoxChart} for what each takes.
##
## @seealso{stats.chart.BoxChart, boxplot, swarmchart}
## @end deftypefn

function varargout = boxchart (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ## An axes may come first, as it may for every plotting function
  ax = [];
  if (isscalar (varargin{1}) && ishghandle (varargin{1})
      && strcmp (get (varargin{1}, 'type'), 'axes'))
    ax = varargin{1};
    varargin(1) = [];
    if (isempty (varargin))
      print_usage ();
    endif
  endif

  [spec, cgroup, args] = bcParse (varargin);

  if (isempty (cgroup))
    b = stats.chart.BoxChart (ax, spec, args);
  else
    ## One chart per colour group, over its own share of the observations
    [lev, ~, idx] = unique (cgroup(:));
    n = numel (lev);
    co = bcColorOrder (ax);
    for k = 1:n
      take = (idx == k);
      sub = spec;
      sub.YData = spec.YData(take);
      if (! isempty (spec.XData))
        sub.XData = spec.XData(take);
      endif
      col = co(mod (k - 1, rows (co)) + 1, :);
      kargs = [{'BoxFaceColor', col, 'BoxEdgeColor', col, ...
                'BoxMedianLineColor', col, 'WhiskerLineColor', col, ...
                'MarkerColor', col}, args];
      ## Octave's classdef has no 'empty', so the array is grown by assigning
      ## its first element rather than preallocated
      if (k == 1)
        b = stats.chart.BoxChart (ax, sub, kargs);
        ## Every group draws into the axes the first of them settled on
        ax = b.Parent;
      else
        b(k) = stats.chart.BoxChart (ax, sub, kargs);
      endif
    endfor
  endif

  ## The chart is handed back only where it was asked for, so a call made
  ## for the drawing alone prints nothing
  if (nargout > 0)
    varargout{1} = b;
  endif

endfunction

## Work out which syntax was used and turn it into what the chart takes.
function [spec, cgroup, args] = bcParse (in)

  spec = struct ('XData', [], 'YData', [], 'SourceTable', [], ...
                 'XVariable', [], 'YVariable', [], ...
                 'XDataMode', 'auto', 'YDataMode', 'manual');
  cgroup = [];

  if (istable (in{1}))
    [spec, args] = bcParseTable (in, spec);
  else
    [spec, args] = bcParseVectors (in, spec);
  endif

  ## 'GroupByColor' is a grouping rather than a property of the chart
  keep = true (1, numel (args));
  for k = 1:2:numel (args)
    if (ischar (args{k}) && strcmpi (args{k}, 'GroupByColor'))
      if (! isempty (spec.SourceTable))
        error (strcat ("boxchart: 'GroupByColor' is not supported with", ...
                       " table input."));
      endif
      cgroup = args{k+1};
      keep(k:k+1) = false;
    endif
  endfor
  args = args(keep);

  if (! isempty (cgroup) && numel (cgroup) != numel (spec.YData))
    error (strcat ("boxchart: 'GroupByColor' must hold one value per", ...
                   " observation."));
  endif

endfunction

## boxchart (tbl, yvar) and boxchart (tbl, xvar, yvar).
function [spec, args] = bcParseTable (in, spec)

  t = in{1};
  rest = in(2:end);
  if (isempty (rest) || ! bcIsName (rest{1}))
    error ("boxchart: a table must be followed by the variable to draw.");
  endif
  if (numel (rest) > 1 && bcIsName (rest{2}))
    xvar = rest{1};
    yvar = rest{2};
    rest = rest(3:end);
  else
    xvar = [];
    yvar = rest{1};
    rest = rest(2:end);
  endif

  spec.SourceTable = t;
  spec.XVariable = xvar;
  spec.YVariable = yvar;
  spec.YDataMode = 'auto';
  spec.XDataMode = 'auto';
  names = t.Properties.VariableNames;
  if (! any (strcmp (names, yvar)))
    error ("boxchart: the table holds no variable '%s'.", yvar);
  endif
  spec.YData = double (t.(yvar)(:));
  if (! isempty (xvar))
    if (! any (strcmp (names, xvar)))
      error ("boxchart: the table holds no variable '%s'.", xvar);
    endif
    spec.XData = t.(xvar);
  endif
  args = rest;

endfunction

## boxchart (ydata) and boxchart (xgroupdata, ydata).
function [spec, args] = bcParseVectors (in, spec)

  if (numel (in) > 1 && isnumeric (in{2}) && ! bcIsName (in{1}))
    spec.XData = in{1};
    spec.YData = double (in{2}(:));
    spec.XDataMode = 'manual';
    args = in(3:end);
  else
    spec.YData = double (in{1}(:));
    args = in(2:end);
  endif
  if (! (isnumeric (spec.YData) && isreal (spec.YData)))
    error ("boxchart: YDATA must be a real numeric vector.");
  endif
  if (! isempty (spec.XData) && numel (spec.XData) != numel (spec.YData))
    error (strcat ("boxchart: XGROUPDATA must hold one value per", ...
                   " observation."));
  endif
  if (mod (numel (args), 2) != 0)
    error ("boxchart: name-value arguments must be in pairs.");
  endif

endfunction

## The colour order the chart will draw with, without making an axes for it.
function co = bcColorOrder (ax)

  if (isempty (ax))
    hf = get (0, 'currentfigure');
    if (! isempty (hf))
      ax = get (hf, 'currentaxes');
    endif
  endif
  if (isempty (ax))
    co = get (0, 'defaultaxescolororder');
  else
    co = get (ax, 'ColorOrder');
  endif

endfunction

## Whether an argument is the name of a table variable.
function tf = bcIsName (x)

  tf = ((ischar (x) && isrow (x)) || (isa (x, 'string') && isscalar (x)));

endfunction

%!demo
%! ## One box over a sample: the box spans the quartiles with the median
%! ## across it, the whiskers reach the furthest values within one and a half
%! ## interquartile ranges, and the rest are drawn as markers.
%! load fisheriris
%! boxchart (meas(:,1));
%! ylabel ('sepal length');

%!demo
%! ## One box per group, and a notch about each median.  Two notches that do
%! ## not overlap mark medians that differ at roughly the five per cent level.
%! load fisheriris
%! g = grp2idx (species);
%! boxchart (g, meas(:,1), 'Notch', 'on');
%! xlabel ('species');
%! ylabel ('sepal length');

%!demo
%! ## The chart keeps what it was drawn from, so its properties may be set
%! ## afterwards and it redraws itself.
%! load fisheriris
%! b = boxchart (meas(:,1));
%! b.Orientation = 'horizontal';
%! b.BoxFaceColor = [0.85, 0.325, 0.098];
%! b.JitterOutliers = 'on';

%!shared bcY, bcG, bcC, bcT
%! bcY = [1; 2; 3; 4; 5; 6; 7; 8; 9; 100];
%! bcG = [1; 1; 1; 1; 1; 2; 2; 2; 2; 2];
%! bcC = {'a'; 'b'; 'a'; 'b'; 'a'; 'b'; 'a'; 'b'; 'a'; 'b'};
%! bcT = table (bcG, bcY, 'VariableNames', {'grp', 'val'});

## MATLAB parity: one box is drawn from the quartiles and the median
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (bcY);
%!   assert_equal (class (b), 'stats.chart.BoxChart');
%!   assert_equal (quantile (bcY, 0.25), 3);
%!   assert_equal (quantile (bcY, 0.75), 8);
%!   assert_equal (median (bcY), 5.5);
%!   assert_equal (sum (isoutlier (bcY, 'quartiles')), 1);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # MATLAB parity: the defaults are the ones R2026a reports
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (bcY);
%!   assert_equal (b.BoxWidth, 0.5);
%!   assert_equal (b.BoxFaceAlpha, 0.2);
%!   assert_equal (b.CapWidth, 0.25);
%!   assert_equal (b.LineWidth, 1);
%!   assert_equal (b.MarkerSize, 6);
%!   assert_equal (b.MarkerStyle, 'o');
%!   assert_equal (b.Notch, 'off');
%!   assert_equal (b.JitterOutliers, 'off');
%!   assert_equal (b.Orientation, 'vertical');
%!   assert_equal (b.WhiskerLineStyle, '-');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a group per distinct value, each with its own box
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   one = boxchart (bcY);
%!   n1 = numel (get (gca (), 'children'));
%!   clf (hf);
%!   two = boxchart (bcG, bcY);
%!   n2 = numel (get (gca (), 'children'));
%!   assert_equal (n2 > n1, true);
%!   assert_equal (two.XDataMode, 'manual');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## MATLAB parity: a table sets both modes to auto and records the columns
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (bcT, 'grp', 'val');
%!   assert_equal (b.YVariable, 'val');
%!   assert_equal (b.XVariable, 'grp');
%!   assert_equal (b.XDataMode, 'auto');
%!   assert_equal (b.YDataMode, 'auto');
%!   assert_equal (b.YData, bcY);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a table may name the observations alone
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (bcT, 'val');
%!   assert_equal (b.YVariable, 'val');
%!   assert_equal (isempty (b.XVariable), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## MATLAB parity: assigning the table again re-reads it
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (bcT, 'grp', 'val');
%!   t = bcT;
%!   t.val(1) = 999;
%!   b.SourceTable = t;
%!   assert_equal (b.YData(1), 999);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## MATLAB parity: one chart per colour group
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (bcG, bcY, 'GroupByColor', bcC);
%!   assert_equal (numel (b), 2);
%!   assert_equal (class (b), 'stats.chart.BoxChart');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # setting a property redraws rather than draws again over the old
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (bcY);
%!   before = numel (get (gca (), 'children'));
%!   b.Notch = 'on';
%!   assert_equal (numel (get (gca (), 'children')), before);
%!   b.BoxFaceColor = 'r';
%!   assert_equal (numel (get (gca (), 'children')), before);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # setting a colour records that it was chosen
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (bcY);
%!   assert_equal (b.BoxFaceColorMode, 'auto');
%!   b.BoxFaceColor = [1, 0, 0];
%!   assert_equal (b.BoxFaceColorMode, 'manual');
%!   assert_equal (b.BoxFaceColor, [1, 0, 0]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a horizontal chart draws the same pieces the other way about
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (bcY);
%!   before = numel (get (gca (), 'children'));
%!   b.Orientation = 'horizontal';
%!   assert_equal (b.Orientation, 'horizontal');
%!   assert_equal (numel (get (gca (), 'children')), before);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the axes to draw into may be given first, as for any plot
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   a1 = subplot (1, 2, 1);
%!   a2 = subplot (1, 2, 2);
%!   b = boxchart (a2, bcY);
%!   assert_equal (b.Parent, a2);
%!   assert_equal (isempty (get (a1, 'children')), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Input validation
%!error<Invalid call> boxchart ()

%!error<boxchart: 'GroupByColor' is not supported with table input.> ...
%! boxchart (bcT, 'grp', 'val', 'GroupByColor', bcC)

%!error<boxchart: the table holds no variable 'nope'.> ...
%! boxchart (bcT, 'nope')

%!error<boxchart: XGROUPDATA must hold one value per observation.> ...
%! boxchart ([1; 2], bcY)

%!error<boxchart: 'GroupByColor' must hold one value per observation.> ...
%! boxchart (bcG, bcY, 'GroupByColor', {'a'; 'b'})

%!error<boxchart: name-value arguments must be in pairs.> ...
%! boxchart (bcY, 'Notch')

%!error<stats.chart.BoxChart: 'BoxWidth' must be a positive scalar.> ...
%! boxchart (bcY, 'BoxWidth', -1)

%!error<stats.chart.BoxChart: 'Notch' must be one of 'on', 'off'.> ...
%! boxchart (bcY, 'Notch', 'maybe')

%!error<stats.chart.BoxChart: 'BoxFaceAlpha' must be a scalar between 0 and 1.> ...
%! boxchart (bcY, 'BoxFaceAlpha', 2)
