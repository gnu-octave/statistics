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
## @deftypefn  {statistics} {} parallelplot (@var{data})
## @deftypefnx {statistics} {} parallelplot (@var{tbl})
## @deftypefnx {statistics} {} parallelplot (@var{parent}, @dots{})
## @deftypefnx {statistics} {} parallelplot (@dots{}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {@var{h} =} parallelplot (@dots{})
##
## Parallel coordinates plot of observations of several variables.
##
## @code{parallelplot (@var{data})} draws each column of the numeric matrix
## @var{data} as a vertical ruler, a coordinate, and each row as a line
## joining its values across the rulers.  By default each ruler has limits of
## its own, from the smallest to the largest value of its column.
##
## @code{parallelplot (@var{tbl})} draws the variables of the table
## @var{tbl} as coordinates, named after them.  A variable of categories,
## text or logical values places its distinct values evenly along its ruler.
##
## @qcode{'CoordinateData'} chooses the columns of @var{data} to draw, and
## @qcode{'CoordinateVariables'} the variables of @var{tbl}, by names,
## indices or a logical vector.  With @qcode{'GroupData'}, a vector with one
## element for each observation, or @qcode{'GroupVariable'}, a variable of
## @var{tbl}, each group is drawn in a colour of its own, with a legend
## naming the groups in the order in which each first appears.
## @qcode{'DataNormalization'} places the values as they are on rulers of
## their own, @qcode{'range'}, or on one shared scale: as they are,
## @qcode{'none'}; as z-scores, @qcode{'zscore'}; divided by the standard
## deviation, @qcode{'scale'}; less the mean, @qcode{'center'}; divided by
## the 2-norm, @qcode{'norm'}.
##
## @code{parallelplot (@var{parent}, @dots{})} places the chart in the
## figure or panel @var{parent} rather than the current figure.  The chart
## takes the place of the current axes of its parent, and is refused where
## @code{hold} is on for them.
##
## @code{parallelplot (@dots{}, @var{Name}, @var{Value})} sets properties of
## the chart; see @code{stats.chart.ParallelCoordinatesPlot} for all of them.
##
## @code{@var{h} = parallelplot (@dots{})} returns the chart, a
## @code{stats.chart.ParallelCoordinatesPlot} object.
##
## @seealso{stats.chart.ParallelCoordinatesPlot, parallelcoords}
## @end deftypefn

function varargout = parallelplot (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  ## Input validation
  parent = [];
  args = varargin;
  if (isscalar (args{1}) && ! istable (args{1}) && isnumeric (args{1})
      && ishghandle (args{1}) && numel (args) > 1
      && (istable (args{2}) || (isnumeric (args{2}) && ! isscalar (args{2}))))
    type = get (args{1}, 'type');
    if (strcmp (type, 'axes'))
      error (strcat ("parallelplot: a parallel coordinates plot cannot be", ...
                     " placed in axes."));
    elseif (! any (strcmp (type, {'figure', 'uipanel'})))
      error ("parallelplot: PARENT must be a figure or a panel.");
    endif
    parent = args{1};
    args(1) = [];
  endif

  spec = struct ('Data', [], 'SourceTable', [], 'CoordinateData', [], ...
                 'CoordinateVariables', {{}}, 'GroupData', [], ...
                 'GroupVariable', '');
  if (istable (args{1}))
    spec.SourceTable = args{1};
  elseif ((isnumeric (args{1}) || islogical (args{1})) && isreal (args{1})
          && ndims (args{1}) == 2)
    spec.Data = args{1};
  elseif (isnumeric (args{1}))
    error ("parallelplot: DATA must be a 2-dimensional real numeric matrix.");
  else
    error ("parallelplot: the data must be a table or a numeric matrix.");
  endif
  rest = args(2:end);

  ## Name-value pairs, the names resolved to the properties they set
  if (mod (numel (rest), 2) != 0)
    error ("parallelplot: optional arguments must be in Name, Value pairs.");
  endif
  props = {'Data', 'CoordinateData', 'SourceTable', ...
           'CoordinateVariables', 'GroupData', 'GroupVariable', ...
           'DataNormalization', 'CoordinateTickLabels', 'Jitter', 'Color', ...
           'LineStyle', 'LineWidth', 'LineAlpha', 'MarkerStyle', ...
           'MarkerSize', 'LegendVisible', 'LegendTitle', ...
           'CoordinateLabel', 'DataLabel', 'Title', 'FontName', 'FontSize', ...
           'Position', 'InnerPosition', 'OuterPosition', ...
           'PositionConstraint', 'Units', 'Visible'};
  for k = 1:2:numel (rest)
    i = [];
    if (isName (rest{k}))
      i = find (strcmpi (char (rest{k}), props), 1);
    endif
    if (isempty (i))
      error ("parallelplot: invalid optional paired argument.");
    endif
    rest{k} = props{i};
  endfor
  ## What chooses and groups the observations belongs with the data
  for f = {'CoordinateData', 'CoordinateVariables', 'GroupData', ...
           'GroupVariable'}
    k = find (strcmp (rest(1:2:end), f{1}), 1, 'last');
    if (! isempty (k))
      fromTable = any (strcmp (f{1}, {'CoordinateVariables', ...
                                      'GroupVariable'}));
      if (fromTable && isempty (spec.SourceTable))
        error ("parallelplot: '%s' needs a table.", f{1});
      elseif (! fromTable && ! isempty (spec.SourceTable)
              && strcmp (f{1}, 'CoordinateData'))
        error ("parallelplot: 'CoordinateData' needs a matrix.");
      endif
      spec.(f{1}) = rest{2*k};
      rest(2*k-1:2*k) = [];
    endif
  endfor
  if (! isempty (spec.SourceTable) && ! isempty (spec.GroupData))
    ## GroupData with a table groups its rows as with a matrix
    rest(end+1:end+2) = {'GroupData', spec.GroupData};
    spec.GroupData = [];
  endif

  h = stats.chart.ParallelCoordinatesPlot (parent, spec, rest);

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
%! ## The four measurements of each iris flower, joined across four rulers.
%! load fisheriris
%! parallelplot (meas);

%!demo
%! ## Grouped by species, with the columns named, as z-scores on one scale.
%! load fisheriris
%! parallelplot (meas, 'GroupData', species, 'CoordinateTickLabels', ...
%!               {'Sepal length', 'Sepal width', 'Petal length', ...
%!                'Petal width'}, 'DataNormalization', 'zscore');

%!demo
%! ## From a table, with a variable of categories as a coordinate.
%! load carsmall
%! t = table (MPG, Horsepower, Weight, categorical (cellstr (Origin)), ...
%!            'VariableNames', {'MPG', 'Horsepower', 'Weight', 'Origin'});
%! parallelplot (t, 'GroupVariable', 'Origin', 'CoordinateVariables', ...
%!               {'MPG', 'Horsepower', 'Weight'});

%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   parallelplot (magic (4));
%!   tag = 'stats.chart.ParallelCoordinatesPlot';
%!   ax = findall (hf, 'type', 'axes', 'tag', tag);
%!   assert_equal (numel (ax), 1);
%!   assert_equal (gca (), ax);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   plot (1:3);
%!   pax = gca ();
%!   h = parallelplot (magic (4));
%!   assert_equal (ishghandle (pax), false);
%!   plot (1:3);
%!   assert_equal (isempty (h.Parent), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!error<Invalid call to parallelplot> parallelplot ()
%!error<parallelplot: the data must be a table or a numeric matrix.> ...
%! parallelplot ('abc')
%!error<parallelplot: DATA must be a 2-dimensional real numeric matrix.> ...
%! parallelplot (ones (2, 2, 2))
%!error<parallelplot: optional arguments must be in Name, Value pairs.> ...
%! parallelplot (magic (3), 'Title')
%!error<parallelplot: invalid optional paired argument.> ...
%! parallelplot (magic (3), 'Foo', 1)
%!error<parallelplot: 'CoordinateVariables' needs a table.> ...
%! parallelplot (magic (3), 'CoordinateVariables', {'A'})
%!error<parallelplot: 'GroupVariable' needs a table.> ...
%! parallelplot (magic (3), 'GroupVariable', 'G')
%!error<parallelplot: 'CoordinateData' needs a matrix.> ...
%! parallelplot (table ([1; 2]), 'CoordinateData', 1)
%!error<parallelplot: a parallel coordinates plot cannot be placed in axes.> ...
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   parallelplot (axes (hf), magic (3));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!error<stats.chart.ParallelCoordinatesPlot: a parallel coordinates plot cannot be added to axes on which hold is on.> ...
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   plot (1:3);
%!   hold on;
%!   parallelplot (magic (3));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
