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
## @deftypefn  {statistics} {} scatterhistogram (@var{x}, @var{y})
## @deftypefnx {statistics} {} scatterhistogram (@var{tbl}, @var{xvar}, @var{yvar})
## @deftypefnx {statistics} {} scatterhistogram (@var{parent}, @dots{})
## @deftypefnx {statistics} {} scatterhistogram (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {@var{h} =} scatterhistogram (@dots{})
##
## Scatter plot of two samples with the histogram of each along its sides.
##
## @code{scatterhistogram (@var{x}, @var{y})} draws the points
## @code{(@var{x}(i), @var{y}(i))} with a histogram of @var{x} above the
## scatter plot and a histogram of @var{y} to its right, each normalized as a
## probability density, with the bins @code{histcounts} chooses.  @var{x} and
## @var{y} are numeric vectors of the same length; a missing value is left out
## of the drawing.
##
## @code{scatterhistogram (@var{tbl}, @var{xvar}, @var{yvar})} takes the
## samples from the variables @var{xvar} and @var{yvar} of the table
## @var{tbl}, and labels the axes with their names.
##
## With @qcode{'GroupData'}, a vector of categories, text, numbers or logical
## values with one element for each observation, or @qcode{'GroupVariable'},
## a variable of @var{tbl} holding them, each group is drawn in a colour of
## its own, with histograms of its own and a legend naming the groups in the
## order in which each first appears.
##
## @code{scatterhistogram (@var{parent}, @dots{})} places the chart in the
## figure or panel @var{parent} rather than the current figure.  The chart
## takes the place of the current axes of its parent, and is refused where
## @code{hold} is on for them.
##
## @code{scatterhistogram (@dots{}, @var{Name}, @var{Value})} sets
## properties of the chart; see @code{stats.chart.ScatterHistogramChart} for
## all of them, among which @qcode{'HistogramDisplayStyle'},
## @qcode{'stairs'}, @qcode{'bar'} or @qcode{'smooth'} for a kernel density
## estimate, @qcode{'NumBins'}, @qcode{'BinWidths'} and
## @qcode{'ScatterPlotLocation'}.
##
## @code{@var{h} = scatterhistogram (@dots{})} returns the chart, a
## @code{stats.chart.ScatterHistogramChart} object.
##
## @seealso{stats.chart.ScatterHistogramChart, scatterhist, histcounts,
## ksdensity}
## @end deftypefn

function varargout = scatterhistogram (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ## Input validation
  parent = [];
  args = varargin;
  if (isscalar (args{1}) && ! istable (args{1}) && isnumeric (args{1})
      && ishghandle (args{1}) && numel (args) > 1 && ! isName (args{2})
      && ! isscalar (args{2}))
    type = get (args{1}, 'type');
    if (strcmp (type, 'axes'))
      error ("scatterhistogram: a scatter histogram cannot be placed in axes.");
    elseif (! any (strcmp (type, {'figure', 'uipanel'})))
      error ("scatterhistogram: PARENT must be a figure or a panel.");
    endif
    parent = args{1};
    args(1) = [];
  endif

  spec = struct ('XData', [], 'YData', [], 'GroupData', [], ...
                 'SourceTable', [], 'XVariable', '', 'YVariable', '', ...
                 'GroupVariable', '');
  if (numel (args) > 0 && istable (args{1}))
    if (numel (args) < 3 || ! isName (args{2}) || ! isName (args{3}))
      error (strcat ("scatterhistogram: a table must be followed by two", ...
                     " variable names."));
    endif
    spec.SourceTable = args{1};
    spec.XVariable = char (args{2});
    spec.YVariable = char (args{3});
    rest = args(4:end);
  else
    if (numel (args) < 2 || ! isSample (args{1}) || ! isSample (args{2}))
      error (strcat ("scatterhistogram: the data must be two numeric", ...
                     " vectors or a table and two variable names."));
    endif
    if (numel (args{1}) != numel (args{2}))
      error ("scatterhistogram: X and Y must be vectors of the same length.");
    endif
    spec.XData = args{1};
    spec.YData = args{2};
    rest = args(3:end);
  endif

  ## Name-value pairs, the names resolved to the properties they set
  if (mod (numel (rest), 2) != 0)
    error (strcat ("scatterhistogram: optional arguments must be in Name,", ...
                   " Value pairs."));
  endif
  props = {'XData', 'YData', 'GroupData', 'SourceTable', 'XVariable', ...
           'YVariable', 'GroupVariable', 'HistogramDisplayStyle', ...
           'NumBins', 'BinWidths', 'ScatterPlotLocation', ...
           'ScatterPlotProportion', 'XHistogramDirection', ...
           'YHistogramDirection', 'XLimits', 'YLimits', 'Color', ...
           'LineStyle', 'LineWidth', 'MarkerStyle', 'MarkerSize', ...
           'MarkerFilled', 'MarkerAlpha', 'LegendVisible', 'LegendTitle', ...
           'Title', 'XLabel', 'YLabel', 'FontName', 'FontSize', 'Position', ...
           'InnerPosition', 'OuterPosition', 'PositionConstraint', 'Units', ...
           'Visible'};
  for k = 1:2:numel (rest)
    i = [];
    if (isName (rest{k}))
      i = find (strcmpi (char (rest{k}), props), 1);
    endif
    if (isempty (i))
      error ("scatterhistogram: invalid optional paired argument.");
    endif
    rest{k} = props{i};
  endfor
  ## The groups belong with the data, so that the chart reads them at once
  for f = {'GroupData', 'GroupVariable'}
    k = find (strcmp (rest(1:2:end), f{1}), 1, 'last');
    if (! isempty (k))
      if (strcmp (f{1}, 'GroupVariable') && isempty (spec.SourceTable))
        error ("scatterhistogram: 'GroupVariable' needs a table.");
      endif
      spec.(f{1}) = rest{2*k};
      rest(2*k-1:2*k) = [];
    endif
  endfor

  h = stats.chart.ScatterHistogramChart (parent, spec, rest);

  ## The chart is handed back only where it was asked for
  if (nargout > 0)
    varargout{1} = h;
  endif

endfunction

## A name: a character vector or a string scalar
function tf = isName (v)
  tf = (ischar (v) && isrow (v)) || (isa (v, 'string') && isscalar (v));
endfunction

## A sample: a real numeric or logical vector
function tf = isSample (v)
  tf = (isnumeric (v) || islogical (v)) && isreal (v) ...
       && (isvector (v) || isempty (v));
endfunction

%!demo
%! ## Two measurements of the same flowers, and how each is spread.
%! load fisheriris
%! scatterhistogram (meas(:,1), meas(:,2));
%! xlabel ('');

%!demo
%! ## Grouped by species, each with its own histograms and colour.
%! load fisheriris
%! scatterhistogram (meas(:,3), meas(:,4), 'GroupData', species, ...
%!                   'XLabel', 'Petal length', 'YLabel', 'Petal width');

%!demo
%! ## The same as kernel density estimates, the histograms below and to the
%! ## left of the scatter plot.
%! load fisheriris
%! scatterhistogram (meas(:,3), meas(:,4), 'GroupData', species, ...
%!                   'HistogramDisplayStyle', 'smooth', ...
%!                   'ScatterPlotLocation', 'NorthEast');

%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   scatterhistogram ([1, 2, 3], [4, 5, 7]);
%!   tag = 'stats.chart.ScatterHistogramChart';
%!   ax = findall (hf, 'type', 'axes', 'tag', tag);
%!   assert_equal (numel (ax), 3);
%!   assert_equal (ismember (gca (), ax), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   t = table ([1; 2; 3; 4], [2; 4; 1; 3], {'a'; 'b'; 'a'; 'b'}, ...
%!              'VariableNames', {'X', 'Y', 'G'});
%!   h = scatterhistogram (t, 'X', 'Y', 'GroupVariable', 'G');
%!   assert_equal (h.XLabel, 'X');
%!   assert_equal (h.YLabel, 'Y');
%!   assert_equal (h.LegendTitle, 'G');
%!   assert_equal (h.GroupData, {'a'; 'b'; 'a'; 'b'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   plot (1:3);
%!   pax = gca ();
%!   h = scatterhistogram ([1, 2, 3], [4, 5, 7]);
%!   assert_equal (ishghandle (pax), false);
%!   plot (1:3);
%!   assert_equal (isempty (h.Parent), true);
%!   tag = 'stats.chart.ScatterHistogramChart';
%!   assert_equal (isempty (findall (hf, 'tag', tag)), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!error<Invalid call to scatterhistogram> scatterhistogram ()
%!error<scatterhistogram: the data must be two numeric vectors or a table and two variable names.> ...
%! scatterhistogram ([1, 2, 3])
%!error<scatterhistogram: the data must be two numeric vectors or a table and two variable names.> ...
%! scatterhistogram ('a', 'b')
%!error<scatterhistogram: X and Y must be vectors of the same length.> ...
%! scatterhistogram ([1, 2, 3], [1, 2])
%!error<scatterhistogram: a table must be followed by two variable names.> ...
%! scatterhistogram (table ([1; 2]), 'Var1')
%!error<scatterhistogram: optional arguments must be in Name, Value pairs.> ...
%! scatterhistogram ([1, 2, 3], [4, 5, 6], 'Title')
%!error<scatterhistogram: invalid optional paired argument.> ...
%! scatterhistogram ([1, 2, 3], [4, 5, 6], 'Foo', 1)
%!error<scatterhistogram: 'GroupVariable' needs a table.> ...
%! scatterhistogram ([1, 2, 3], [4, 5, 6], 'GroupVariable', 'G')
%!error<scatterhistogram: a scatter histogram cannot be placed in axes.> ...
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   scatterhistogram (axes (hf), [1, 2, 3], [4, 5, 6]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!error<stats.chart.ScatterHistogramChart: a scatter histogram cannot be added to axes on which hold is on.> ...
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   plot (1:3);
%!   hold on;
%!   scatterhistogram ([1, 2, 3], [4, 5, 6]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
