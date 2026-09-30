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
## @deftypefn  {statistics} {} heatmap (@var{cdata})
## @deftypefnx {statistics} {} heatmap (@var{xvalues}, @var{yvalues}, @var{cdata})
## @deftypefnx {statistics} {} heatmap (@var{tbl}, @var{xvar}, @var{yvar})
## @deftypefnx {statistics} {} heatmap (@var{parent}, @dots{})
## @deftypefnx {statistics} {} heatmap (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {@var{h} =} heatmap (@dots{})
##
## Heatmap of a matrix, or of the rows of a table counted or aggregated over
## two of its variables.
##
## @code{heatmap (@var{cdata})} draws the matrix @var{cdata} as a grid of
## coloured cells, each carrying its value, with a colour bar beside them.
## The columns are named @qcode{'1'}, @qcode{'2'}, @dots{}, and so are the
## rows.
##
## @code{heatmap (@var{xvalues}, @var{yvalues}, @var{cdata})} names the columns
## with @var{xvalues} and the rows with @var{yvalues}, vectors of distinct
## numbers, text or categories with one element for each column and each row
## of @var{cdata}.
##
## @code{heatmap (@var{tbl}, @var{xvar}, @var{yvar})} makes a column of each
## distinct value of the table variable @var{xvar} and a row of each distinct
## value of @var{yvar}, and colours each cell by the number of rows of
## @var{tbl} falling in it.  With @qcode{'ColorVariable'} naming a numeric
## variable, each cell holds the mean of that variable over its rows instead,
## or what @qcode{'ColorMethod'} asks for: @qcode{'count'}, @qcode{'mean'},
## @qcode{'median'}, @qcode{'sum'}, @qcode{'min'}, @qcode{'max'}, or
## @qcode{'none'} where each cell holds at most one row.  The categories of a
## @code{categorical} variable are taken in their order, those no row uses
## included; numbers, text and logical values are sorted.  A row whose
## @var{xvar} or @var{yvar} is missing is left out.
##
## @code{heatmap (@var{parent}, @dots{})} places the heatmap in the figure or
## panel @var{parent} rather than the current figure.  The heatmap takes the
## place of the current axes of its parent, and is refused where
## @code{hold} is on for them.
##
## @code{heatmap (@dots{}, @var{Name}, @var{Value})} sets properties of the
## heatmap; see @code{stats.chart.HeatmapChart} for all of them.
##
## @code{@var{h} = heatmap (@dots{})} returns the heatmap, a
## @code{stats.chart.HeatmapChart} object, whose methods @code{sortx},
## @code{sorty}, @code{xlim} and @code{ylim} reorder and restrict its columns
## and rows.
##
## @seealso{stats.chart.HeatmapChart, imagesc, crosstab}
## @end deftypefn

function varargout = heatmap (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ## Input validation
  parent = [];
  args = varargin;
  if (isscalar (args{1}) && ! istable (args{1}) && isnumeric (args{1}) ...
      && ishghandle (args{1}))
    type = get (args{1}, 'type');
    if (strcmp (type, 'axes'))
      error ("heatmap: a heatmap cannot be placed in axes.");
    elseif (! any (strcmp (type, {'figure', 'uipanel'})))
      error ("heatmap: PARENT must be a figure or a panel.");
    endif
    parent = args{1};
    args(1) = [];
  endif
  if (isempty (args))
    print_usage ();
  endif

  spec = struct ('ColorData', [], 'XData', {{}}, 'YData', {{}}, ...
                 'SourceTable', [], 'XVariable', '', 'YVariable', '', ...
                 'ColorVariable', '');
  if (istable (args{1}))
    if (numel (args) < 3 || ! isName (args{2}) || ! isName (args{3}))
      error ("heatmap: a table must be followed by two variable names.");
    endif
    spec.SourceTable = args{1};
    spec.XVariable = char (args{2});
    spec.YVariable = char (args{3});
    rest = args(4:end);
  else
    npos = find (cellfun (@isName, args), 1) - 1;
    if (isempty (npos))
      npos = numel (args);
    endif
    if (npos == 1)
      cdata = args{1};
    elseif (npos == 3)
      cdata = args{3};
    else
      error (strcat ("heatmap: the data must be a matrix, two vectors and", ...
                     " a matrix, or a table and two variable names."));
    endif
    if (! ((isnumeric (cdata) || islogical (cdata)) && isreal (cdata)
           && ndims (cdata) == 2))
      error ("heatmap: CDATA must be a 2-dimensional real numeric matrix.");
    endif
    spec.ColorData = cdata;
    if (npos == 3)
      if (! isvector (args{1}) || numel (args{1}) != columns (cdata))
        error (strcat ("heatmap: XVALUES must hold one value for each", ...
                       " column of CDATA."));
      endif
      if (! isvector (args{2}) || numel (args{2}) != rows (cdata))
        error (strcat ("heatmap: YVALUES must hold one value for each", ...
                       " row of CDATA."));
      endif
      spec.XData = args{1};
      spec.YData = args{2};
    endif
    rest = args(npos+1:end);
  endif

  ## Name-value pairs, the names resolved to the properties they set
  if (mod (numel (rest), 2) != 0)
    error ("heatmap: optional arguments must be in Name, Value pairs.");
  endif
  props = {'ColorData', 'XData', 'YData', 'XDisplayData', 'YDisplayData', ...
           'XDisplayLabels', 'YDisplayLabels', 'XLimits', 'YLimits', ...
           'SourceTable', 'XVariable', 'YVariable', 'ColorVariable', ...
           'ColorMethod', 'ColorScaling', 'ColorLimits', 'Colormap', ...
           'ColorbarVisible', 'MissingDataColor', 'MissingDataLabel', ...
           'CellLabelFormat', 'CellLabelColor', 'GridVisible', 'Title', ...
           'XLabel', 'YLabel', 'FontName', 'FontSize', 'FontColor', ...
           'Interpreter', 'Position', 'InnerPosition', 'OuterPosition', ...
           'PositionConstraint', ...
           'Units', 'Visible'};
  for k = 1:2:numel (rest)
    i = [];
    if (isName (rest{k}))
      i = find (strcmpi (char (rest{k}), props), 1);
    endif
    if (isempty (i))
      error ("heatmap: invalid optional paired argument.");
    endif
    rest{k} = props{i};
  endfor
  ## The colour variable belongs with the table, so that the chart reads
  ## the table once with everything it needs
  k = find (strcmp (rest(1:2:end), 'ColorVariable'), 1, 'last');
  if (! isempty (k))
    if (isempty (spec.SourceTable))
      error ("heatmap: 'ColorVariable' needs a table.");
    endif
    spec.ColorVariable = rest{2*k};
    rest(2*k-1:2*k) = [];
  endif

  h = stats.chart.HeatmapChart (parent, spec, rest);

  ## The chart is handed back only where it was asked for
  if (nargout > 0)
    varargout{1} = h;
  endif

endfunction

## A name: a character vector or a string scalar
function tf = isName (v)
  tf = (ischar (v) && isrow (v)) || (isa (v, 'string') && isscalar (v));
endfunction

%!demo
%! ## A matrix as a heatmap: each cell carries its value and its colour.
%! heatmap (magic (5));

%!demo
%! ## Columns and rows may be named, and the colours changed.
%! h = heatmap ({'Mon', 'Tue', 'Wed', 'Thu', 'Fri'}, {'am', 'pm'}, ...
%!              [12, 15, 9, 20, 18; 7, 11, 14, 6, 9]);
%! h.Title = 'Calls answered';
%! h.Colormap = gray (64);

%!demo
%! ## From a table: the cars of each origin and number of cylinders, and
%! ## then their mean mileage.
%! load carsmall
%! t = table (categorical (cellstr (Origin)), Cylinders, MPG, ...
%!            'VariableNames', {'Origin', 'Cylinders', 'MPG'});
%! subplot (1, 2, 1);
%! heatmap (t, 'Cylinders', 'Origin');
%! subplot (1, 2, 2);
%! heatmap (t, 'Cylinders', 'Origin', 'ColorVariable', 'MPG');

%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   heatmap (magic (3));
%!   ax = findall (hf, 'tag', 'stats.chart.HeatmapChart');
%!   assert_equal (numel (ax), 1);
%!   assert_equal (gca (), ax);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ({'a', 'b'}, [1, 2, 3], [1, 2; 3, 4; 5, 6], 'Title', 'T', ...
%!                'colormap', gray (4));
%!   assert_equal (h.XData, {'a'; 'b'});
%!   assert_equal (h.YData, {'1'; '2'; '3'});
%!   assert_equal (h.Title, 'T');
%!   assert_equal (h.Colormap, gray (4));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   t = table (categorical ({'a';'b';'a'}), categorical ({'u';'u';'v'}), ...
%!              [1;2;3], 'VariableNames', {'X', 'Y', 'V'});
%!   h = heatmap (t, "X", "Y", 'ColorVariable', 'V', 'ColorMethod', 'sum');
%!   assert_equal (h.ColorData, [1, 2; 3, 0]);
%!   assert_equal (h.Title, 'Sum of V');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   plot (1:3);
%!   pax = gca ();
%!   h = heatmap (magic (3));
%!   assert_equal (ishghandle (pax), false);
%!   h2 = heatmap (magic (4));
%!   assert_equal (numel (findall (hf, 'tag', 'stats.chart.HeatmapChart')), 1);
%!   assert_equal (isempty (h.Parent), true);
%!   plot (1:3);
%!   assert_equal (isempty (h2.Parent), true);
%!   tag = 'stats.chart.HeatmapChart';
%!   assert_equal (isempty (findall (hf, 'tag', tag)), true);
%!   assert_equal (numel (findall (gca (), 'type', 'line')), 1);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   left = subplot (1, 2, 1);
%!   plot (1:3);
%!   subplot (1, 2, 2);
%!   outer = get (gca (), 'outerposition');
%!   h = heatmap (magic (3));
%!   assert_equal (h.OuterPosition, outer, -1e-12);
%!   assert_equal (ishghandle (left), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!error<Invalid call to heatmap> heatmap ()
%!error<heatmap: a heatmap cannot be placed in axes.> ...
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   heatmap (axes (hf), magic (3));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!error<heatmap: a table must be followed by two variable names.> ...
%! heatmap (table ([1; 2]), 'Var1')
%!error<heatmap: the data must be a matrix, two vectors and a matrix, or a table and two variable names.> ...
%! heatmap ([1, 2], [1, 2])
%!error<heatmap: the data must be a matrix, two vectors and a matrix, or a table and two variable names.> ...
%! heatmap ('abc')
%!error<heatmap: CDATA must be a 2-dimensional real numeric matrix.> ...
%! heatmap ({1, 2})
%!error<heatmap: CDATA must be a 2-dimensional real numeric matrix.> ...
%! heatmap (ones (2, 2, 2))
%!error<heatmap: XVALUES must hold one value for each column of CDATA.> ...
%! heatmap ({'a', 'b'}, {'r1'}, [1, 2, 3])
%!error<heatmap: YVALUES must hold one value for each row of CDATA.> ...
%! heatmap ({'a', 'b', 'c'}, {'r1', 'r2'}, [1, 2, 3])
%!error<heatmap: optional arguments must be in Name, Value pairs.> ...
%! heatmap (magic (3), 'Title')
%!error<heatmap: invalid optional paired argument.> ...
%! heatmap (magic (3), 'Foo', 1)
%!error<heatmap: 'ColorVariable' needs a table.> ...
%! heatmap (magic (3), 'ColorVariable', 'V')
%!error<stats.chart.HeatmapChart: a heatmap cannot be added to axes on which hold is on.> ...
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   plot (1:3);
%!   hold on;
%!   heatmap (magic (3));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
