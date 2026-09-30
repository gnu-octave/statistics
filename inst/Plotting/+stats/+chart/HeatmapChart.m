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
## @deftp {statistics} stats.chart.HeatmapChart
##
## A heatmap, as @code{heatmap} draws it.
##
## A @code{stats.chart.HeatmapChart} holds a matrix of values, the labels of
## its columns and rows, and the choices it is drawn with, and redraws itself
## whenever one of them is set.  It is what @code{heatmap} returns, and the
## documented way to reach it.
##
## The values come either from a matrix or from a table, never from both.
## From a table, each distinct value of @code{XVariable} is a column and each
## distinct value of @code{YVariable} a row, and each cell holds the number of
## rows of the table falling in it, or the mean, median, sum, minimum or
## maximum of @code{ColorVariable} over them, as @code{ColorMethod} says.
## Missing values of @code{ColorVariable} are left out of the mean, median,
## sum, minimum and maximum, and counted by @qcode{'count'}.
##
## @code{XData} and @code{YData} name the columns and rows of @code{ColorData}.
## @code{XDisplayData} and @code{YDisplayData} say which of them are shown and
## in what order, and @code{XDisplayLabels} and @code{YDisplayLabels} what they
## are called on the chart.  @code{ColorDisplayData} is the matrix as shown.
##
## A heatmap takes the place of the current axes of its parent, as MATLAB's
## does, keeping their outer position, and refuses axes on which
## @code{hold} is on.  A later plot into the figure takes the place of the
## heatmap in turn.
##
## MATLAB's chart is a graphics object of its own.  This one is a handle class
## that draws into axes of its own, which become the current axes and are
## deleted with it; when a later plot clears those axes the chart lets go of
## them, and @code{Parent} and the drawing are gone.  @code{GridVisible},
## @code{ColorbarVisible} and @code{Visible} hold @qcode{'on'} or
## @qcode{'off'} where MATLAB holds an on-off state.
##
## @seealso{heatmap, imagesc}
## @end deftp

classdef HeatmapChart < handle

  properties (Access = public)

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} ColorData
    ##
    ## The matrix of values, one row for each element of @code{YData} and one
    ## column for each element of @code{XData}.  Where the chart was drawn
    ## from a table it is computed from the table and setting it is refused.
    ##
    ## @end deftp
    ColorData = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} XData
    ##
    ## The names of the columns of @code{ColorData}, a column cell array of
    ## distinct character vectors.  By default they are @qcode{'1'},
    ## @qcode{'2'}, @dots{}, and they follow the size of @code{ColorData} until
    ## they are set.
    ##
    ## @end deftp
    XData = cell (0, 1);

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} YData
    ##
    ## The names of the rows of @code{ColorData}, as @code{XData} names its
    ## columns.
    ##
    ## @end deftp
    YData = cell (0, 1);

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} XDisplayData
    ##
    ## The columns shown and their order, as elements of @code{XData}; a value
    ## not in @code{XData} is shown as a column of missing values.  All of
    ## @code{XData} by default.
    ##
    ## @end deftp
    XDisplayData = cell (0, 1);

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} YDisplayData
    ##
    ## The rows shown and their order, as @code{XDisplayData} gives the
    ## columns.
    ##
    ## @end deftp
    YDisplayData = cell (0, 1);

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} XDisplayLabels
    ##
    ## The labels of the columns shown, one for each element of
    ## @code{XDisplayData}, which they are by default.  A label set stays with
    ## its column when @code{XDisplayData} is reordered.
    ##
    ## @end deftp
    XDisplayLabels = cell (0, 1);

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} YDisplayLabels
    ##
    ## The labels of the rows shown, as @code{XDisplayLabels} labels the
    ## columns.
    ##
    ## @end deftp
    YDisplayLabels = cell (0, 1);

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} XLimits
    ##
    ## The first and the last column shown, as a 1-by-2 cell array of elements
    ## of @code{XDisplayData}; the columns between them in
    ## @code{XDisplayData} are shown.  The whole of @code{XDisplayData} by
    ## default.
    ##
    ## @end deftp
    XLimits = cell (1, 0);

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} YLimits
    ##
    ## The first and the last row shown, as @code{XLimits} gives the columns.
    ##
    ## @end deftp
    YLimits = cell (1, 0);

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} SourceTable
    ##
    ## The table the chart was drawn from, or an empty table.
    ##
    ## @end deftp
    SourceTable = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} XVariable
    ##
    ## The variable of @code{SourceTable} whose values are the columns.
    ##
    ## @end deftp
    XVariable = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} YVariable
    ##
    ## The variable of @code{SourceTable} whose values are the rows.
    ##
    ## @end deftp
    YVariable = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} ColorVariable
    ##
    ## The numeric variable of @code{SourceTable} aggregated into each cell,
    ## or empty to count the rows.
    ##
    ## @end deftp
    ColorVariable = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} ColorMethod
    ##
    ## How the rows of @code{SourceTable} falling in a cell make its value:
    ## @qcode{'count'}, @qcode{'mean'}, @qcode{'median'}, @qcode{'sum'},
    ## @qcode{'min'}, @qcode{'max'} or @qcode{'none'}, which takes the one
    ## value there is and is refused where a cell holds more than one.  It
    ## is @qcode{'count'} where no @code{ColorVariable} is named,
    ## @qcode{'mean'} where one is, and @qcode{'none'} for a matrix.
    ##
    ## @end deftp
    ColorMethod = 'none';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} ColorScaling
    ##
    ## How the values are turned into colours: @qcode{'scaled'}, the default,
    ## maps them to the colormap as they are; @qcode{'scaledcolumns'} and
    ## @qcode{'scaledrows'} rescale each column or row to run from 0 to 1;
    ## @qcode{'log'} maps their natural logarithm, a negative value then
    ## being left uncoloured with a warning.  The cell labels show the values
    ## themselves whatever the scaling.
    ##
    ## @end deftp
    ColorScaling = 'scaled';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} ColorLimits
    ##
    ## The values mapped to the first and the last colour of @code{Colormap},
    ## as a 1-by-2 increasing vector.  Until it is set it follows the values
    ## shown: their smallest and largest, @math{v - 1} and @math{v + 1} where
    ## they are all @math{v}, @math{[0, 1]} where none is finite or where the
    ## columns or rows are rescaled.
    ##
    ## @end deftp
    ColorLimits = [0, 1];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} Colormap
    ##
    ## The colours, as an M-by-3 matrix of values between 0 and 1.  By default
    ## 256 colours running from a pale blue to @code{[0, 0.447, 0.741]}.
    ##
    ## @end deftp
    Colormap = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} ColorbarVisible
    ##
    ## Whether the colour bar is drawn, @qcode{'on'} by default or
    ## @qcode{'off'}.
    ##
    ## @end deftp
    ColorbarVisible = 'on';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} MissingDataColor
    ##
    ## The colour of a cell holding a missing value, as an RGB triplet or a
    ## colour name.  @code{[0.15, 0.15, 0.15]} by default.
    ##
    ## @end deftp
    MissingDataColor = [0.15, 0.15, 0.15];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} MissingDataLabel
    ##
    ## The label of a cell holding a missing value, @qcode{'NaN'} by default.
    ##
    ## @end deftp
    MissingDataLabel = 'NaN';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} CellLabelFormat
    ##
    ## The format each value is written into its cell with, as
    ## @code{sprintf} takes it.  @qcode{'%0.4g'} by default.
    ##
    ## @end deftp
    CellLabelFormat = '%0.4g';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} CellLabelColor
    ##
    ## The colour of the cell labels, as an RGB triplet or a colour name;
    ## @qcode{'auto'}, the default, writes black on a light cell and white on
    ## a dark one, and @qcode{'none'} writes no labels.
    ##
    ## @end deftp
    CellLabelColor = 'auto';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} GridVisible
    ##
    ## Whether lines are drawn between the cells, @qcode{'on'} by default or
    ## @qcode{'off'}.
    ##
    ## @end deftp
    GridVisible = 'on';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} Title
    ##
    ## The title, a character vector or a cell array of them, one for each
    ## line.  From a table it says what the cells hold, such as
    ## @qcode{'Count of Y vs. X'} or @qcode{'Mean of V'}, until it is set.
    ##
    ## @end deftp
    Title = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} XLabel
    ##
    ## The label of the horizontal axis; from a table the name of
    ## @code{XVariable}, until it is set.
    ##
    ## @end deftp
    XLabel = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} YLabel
    ##
    ## The label of the vertical axis; from a table the name of
    ## @code{YVariable}, until it is set.
    ##
    ## @end deftp
    YLabel = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} FontName
    ##
    ## The font of every text on the chart, @qcode{'Helvetica'} by default.
    ##
    ## @end deftp
    FontName = 'Helvetica';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} FontSize
    ##
    ## The size of every text on the chart, in points, 10 by default.
    ##
    ## @end deftp
    FontSize = 10;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} FontColor
    ##
    ## The colour of the title, the axis labels and the tick labels.
    ## @code{[0.15, 0.15, 0.15]} by default.
    ##
    ## @end deftp
    FontColor = [0.15, 0.15, 0.15];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} Interpreter
    ##
    ## How the texts are interpreted: @qcode{'tex'}, the default,
    ## @qcode{'latex'} or @qcode{'none'}.
    ##
    ## @end deftp
    Interpreter = 'tex';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} Position
    ##
    ## The position of the grid of cells in its parent, as
    ## @code{[left, bottom, width, height]} in @code{Units}.  While
    ## @code{PositionConstraint} is @qcode{'outerposition'} the chart places
    ## the cells within @code{OuterPosition} itself, leaving room for the
    ## labels and the colour bar, @code{[0.13, 0.11, 0.732143, 0.815]} of it
    ## at least; setting @code{Position} fixes the cells there and turns
    ## @code{PositionConstraint} to @qcode{'innerposition'}.
    ## @code{InnerPosition} is the same.
    ##
    ## @end deftp
    Position = [0.13, 0.11, 0.732143, 0.815];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} OuterPosition
    ##
    ## The position of the whole chart in its parent, labels included,
    ## @code{[0, 0, 1, 1]} by default, or the outer position of the axes it
    ## took the place of.  Setting it turns @code{PositionConstraint} to
    ## @qcode{'outerposition'}.
    ##
    ## @end deftp
    OuterPosition = [0, 0, 1, 1];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} PositionConstraint
    ##
    ## Which position the chart keeps when its labels change:
    ## @qcode{'outerposition'}, the default, or @qcode{'innerposition'}.
    ##
    ## @end deftp
    PositionConstraint = 'outerposition';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} Units
    ##
    ## The units of the positions, @qcode{'normalized'} by default.
    ##
    ## @end deftp
    Units = 'normalized';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} Visible
    ##
    ## Whether the chart is shown, @qcode{'on'} by default or @qcode{'off'}.
    ##
    ## @end deftp
    Visible = 'on';

  endproperties

  properties (Dependent)

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} ColorDisplayData
    ##
    ## @code{ColorData} as it is shown: its rows and columns in the order of
    ## @code{YDisplayData} and @code{XDisplayData}, missing where a displayed
    ## value is not in the data.  This property is read-only.
    ##
    ## @end deftp
    ColorDisplayData;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} InnerPosition
    ##
    ## The same as @code{Position}.
    ##
    ## @end deftp
    InnerPosition;

  endproperties

  properties (GetAccess = public, SetAccess = private)

    ## -*- texinfo -*-
    ## @deftp {stats.chart.HeatmapChart} {property} Parent
    ##
    ## The figure or panel the chart is placed in.  This property is
    ## read-only.
    ##
    ## @end deftp
    Parent = [];

  endproperties

  properties (Access = private, Hidden)
    Axes_ = [];               # the axes the chart draws into
    Colorbar_ = [];           # its colour bar
    Drawn_ = false;           # false while the constructor fills the fields
    Redrawing_ = false;       # true while the chart clears its own drawing
    Aggregating_ = false;     # true while the table is being read
    Internal_ = false;        # true while the chart sets a value itself
    XDataAuto_ = true;        # XData follows the size of ColorData
    YDataAuto_ = true;
    XDisplayAuto_ = true;     # XDisplayData follows XData
    YDisplayAuto_ = true;
    XLabelMap_ = cell (0, 2); # labels set by hand, value against label
    YLabelMap_ = cell (0, 2);
    XLimitsAuto_ = true;
    YLimitsAuto_ = true;
    ColorLimitsAuto_ = true;
    TitleAuto_ = true;
    XLabelAuto_ = true;
    YLabelAuto_ = true;
    MethodAuto_ = true;       # ColorMethod follows ColorVariable
  endproperties

  methods (Hidden)

    function disp (this)
      if (isempty (this.Title))
        printf ('  HeatmapChart with properties:\n\n');
      else
        t = this.Title;
        if (iscell (t))
          t = strjoin (t, ' ');
        endif
        printf ('  HeatmapChart (%s) with properties:\n\n', t);
      endif
      if (isempty (this.SourceTable))
        printf ('%13s: {%dx1 cell}\n', 'XData', numel (this.XData));
        printf ('%13s: {%dx1 cell}\n', 'YData', numel (this.YData));
        printf ('%13s: [%dx%d double]\n', 'ColorData', ...
                rows (this.ColorData), columns (this.ColorData));
      else
        printf ('%13s: [%dx%d table]\n', 'SourceTable', ...
                size (this.SourceTable));
        printf ('%13s: ''%s''\n', 'XVariable', this.XVariable);
        printf ('%13s: ''%s''\n', 'YVariable', this.YVariable);
        printf ('%13s: ''%s''\n', 'ColorVariable', this.ColorVariable);
        printf ('%13s: ''%s''\n', 'ColorMethod', this.ColorMethod);
      endif
      printf ('\n');
    endfunction

    function display (this)
      disp (this);
    endfunction

  endmethods

  methods (Hidden)

    function set.ColorData (this, val)
      if (this.Drawn_ && ! isempty (this.SourceTable))
        error (strcat ("stats.chart.HeatmapChart: setting 'ColorData'", ...
                       " while 'SourceTable' holds a table is not", ...
                       " supported."));
      endif
      this.ColorData = hmCheckMatrix (val);
      fitData (this);
      redraw (this);
    endfunction

    function set.XData (this, val)
      this.XData = hmCheckNames (val, 'XData');
      this.XDataAuto_ = false;
      fitData (this);
      redraw (this);
    endfunction

    function set.YData (this, val)
      this.YData = hmCheckNames (val, 'YData');
      this.YDataAuto_ = false;
      fitData (this);
      redraw (this);
    endfunction

    function set.XDisplayData (this, val)
      this.XDisplayData = hmCheckNames (val, 'XDisplayData');
      this.XDisplayAuto_ = false;
      this.Internal_ = true;
      unwind_protect
        this.XLimits = hmFitLimits (this.XLimits, this.XDisplayData, ...
                                    this.XLimitsAuto_, 'XLimits');
      unwind_protect_cleanup
        this.Internal_ = false;
      end_unwind_protect
      redraw (this);
    endfunction

    function set.YDisplayData (this, val)
      this.YDisplayData = hmCheckNames (val, 'YDisplayData');
      this.YDisplayAuto_ = false;
      this.Internal_ = true;
      unwind_protect
        this.YLimits = hmFitLimits (this.YLimits, this.YDisplayData, ...
                                    this.YLimitsAuto_, 'YLimits');
      unwind_protect_cleanup
        this.Internal_ = false;
      end_unwind_protect
      redraw (this);
    endfunction

    function val = get.XDisplayLabels (this)
      val = hmLabels (this.XDisplayData, this.XLabelMap_);
    endfunction

    function set.XDisplayLabels (this, val)
      val = hmCheckLabels (val, numel (this.XDisplayData), 'XDisplayLabels');
      this.XLabelMap_ = hmSetLabels (this.XLabelMap_, this.XDisplayData, val);
      redraw (this);
    endfunction

    function val = get.YDisplayLabels (this)
      val = hmLabels (this.YDisplayData, this.YLabelMap_);
    endfunction

    function set.YDisplayLabels (this, val)
      val = hmCheckLabels (val, numel (this.YDisplayData), 'YDisplayLabels');
      this.YLabelMap_ = hmSetLabels (this.YLabelMap_, this.YDisplayData, val);
      redraw (this);
    endfunction

    function set.XLimits (this, val)
      this.XLimits = hmCheckLimits (val, this.XDisplayData, 'XLimits');
      if (! this.Internal_)
        if (this.Drawn_)
          this.XLimitsAuto_ = false;
        endif
        redraw (this);
      endif
    endfunction

    function set.YLimits (this, val)
      this.YLimits = hmCheckLimits (val, this.YDisplayData, 'YLimits');
      if (! this.Internal_)
        if (this.Drawn_)
          this.YLimitsAuto_ = false;
        endif
        redraw (this);
      endif
    endfunction

    function set.SourceTable (this, val)
      if (! (isempty (val) || istable (val)))
        error ("stats.chart.HeatmapChart: 'SourceTable' must be a table.");
      endif
      this.SourceTable = val;
      aggregate (this);
    endfunction

    function set.XVariable (this, val)
      this.XVariable = hmCheckVarName (val, 'XVariable');
      aggregate (this);
    endfunction

    function set.YVariable (this, val)
      this.YVariable = hmCheckVarName (val, 'YVariable');
      aggregate (this);
    endfunction

    function set.ColorVariable (this, val)
      this.ColorVariable = hmCheckVarName (val, 'ColorVariable');
      aggregate (this);
    endfunction

    function set.ColorMethod (this, val)
      this.ColorMethod = hmCheckOneOf (val, {'count', 'mean', 'median', ...
                                             'sum', 'min', 'max', 'none'}, ...
                                       'ColorMethod');
      if (! this.Aggregating_)
        this.MethodAuto_ = false;
        aggregate (this);
      endif
    endfunction

    function set.ColorScaling (this, val)
      this.ColorScaling = hmCheckOneOf (val, {'scaled', 'scaledcolumns', ...
                                              'scaledrows', 'log'}, ...
                                        'ColorScaling');
      redraw (this);
    endfunction

    function set.ColorLimits (this, val)
      if (! (isnumeric (val) && isreal (val) && numel (val) == 2
             && all (isfinite (val)) && val(2) > val(1)))
        error (strcat ("stats.chart.HeatmapChart: 'ColorLimits' must be", ...
                       " a 1-by-2 increasing numeric vector."));
      endif
      this.ColorLimits = double (val(:)');
      if (! this.Internal_)
        if (this.Drawn_)
          this.ColorLimitsAuto_ = false;
        endif
        redraw (this);
      endif
    endfunction

    function set.Colormap (this, val)
      if (! (isnumeric (val) && isreal (val) && ndims (val) == 2
             && columns (val) == 3 && rows (val) > 0
             && all (val(:) >= 0) && all (val(:) <= 1)))
        error (strcat ("stats.chart.HeatmapChart: 'Colormap' must be an", ...
                       " M-by-3 matrix of values between 0 and 1."));
      endif
      this.Colormap = double (val);
      redraw (this);
    endfunction

    function set.ColorbarVisible (this, val)
      this.ColorbarVisible = hmCheckOneOf (val, {'on', 'off'}, ...
                                           'ColorbarVisible');
      redraw (this);
    endfunction

    function set.MissingDataColor (this, val)
      this.MissingDataColor = hmCheckColor (val, 'MissingDataColor', false);
      redraw (this);
    endfunction

    function set.MissingDataLabel (this, val)
      this.MissingDataLabel = hmCheckText (val, 'MissingDataLabel');
      redraw (this);
    endfunction

    function set.CellLabelFormat (this, val)
      if (! (ischar (val) && (isrow (val) || isempty (val)))
          && ! (isa (val, 'string') && isscalar (val)))
        error (strcat ("stats.chart.HeatmapChart: 'CellLabelFormat' must", ...
                       " be a character vector."));
      endif
      this.CellLabelFormat = char (val);
      redraw (this);
    endfunction

    function set.CellLabelColor (this, val)
      this.CellLabelColor = hmCheckColor (val, 'CellLabelColor', true);
      redraw (this);
    endfunction

    function set.GridVisible (this, val)
      this.GridVisible = hmCheckOneOf (val, {'on', 'off'}, 'GridVisible');
      redraw (this);
    endfunction

    function set.Title (this, val)
      this.Title = hmCheckText (val, 'Title');
      if (this.Drawn_)
        this.TitleAuto_ = false;
      endif
      redraw (this);
    endfunction

    function set.XLabel (this, val)
      this.XLabel = hmCheckText (val, 'XLabel');
      if (this.Drawn_)
        this.XLabelAuto_ = false;
      endif
      redraw (this);
    endfunction

    function set.YLabel (this, val)
      this.YLabel = hmCheckText (val, 'YLabel');
      if (this.Drawn_)
        this.YLabelAuto_ = false;
      endif
      redraw (this);
    endfunction

    function set.FontName (this, val)
      this.FontName = hmCheckText (val, 'FontName');
      redraw (this);
    endfunction

    function set.FontSize (this, val)
      if (! (isnumeric (val) && isreal (val) && isscalar (val)
             && isfinite (val) && val > 0))
        error (strcat ("stats.chart.HeatmapChart: 'FontSize' must be a", ...
                       " positive number."));
      endif
      this.FontSize = double (val);
      redraw (this);
    endfunction

    function set.FontColor (this, val)
      this.FontColor = hmCheckColor (val, 'FontColor', false);
      redraw (this);
    endfunction

    function set.Interpreter (this, val)
      this.Interpreter = hmCheckOneOf (val, {'tex', 'latex', 'none'}, ...
                                       'Interpreter');
      redraw (this);
    endfunction

    function set.Position (this, val)
      this.Position = hmCheckPosition (val, 'Position');
      if (! this.Internal_)
        this.PositionConstraint = 'innerposition';
      endif
    endfunction

    function val = get.InnerPosition (this)
      val = this.Position;
    endfunction

    function set.InnerPosition (this, val)
      this.Position = val;
    endfunction

    function set.OuterPosition (this, val)
      this.OuterPosition = hmCheckPosition (val, 'OuterPosition');
      if (! this.Internal_)
        this.PositionConstraint = 'outerposition';
      endif
    endfunction

    function set.PositionConstraint (this, val)
      this.PositionConstraint = hmCheckOneOf (val, {'outerposition', ...
                                                    'innerposition'}, ...
                                              'PositionConstraint');
      redraw (this);
    endfunction

    function set.Units (this, val)
      this.Units = hmCheckOneOf (val, {'normalized', 'inches', ...
                                       'centimeters', 'points', 'pixels', ...
                                       'characters'}, 'Units');
      redraw (this);
    endfunction

    function set.Visible (this, val)
      this.Visible = hmCheckOneOf (val, {'on', 'off'}, 'Visible');
      redraw (this);
    endfunction

    function val = get.ColorDisplayData (this)
      val = hmDisplayMatrix (this.ColorData, this.XData, this.YData, ...
                             this.XDisplayData, this.YDisplayData);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {stats.chart.HeatmapChart} {@var{obj} =} stats.chart.HeatmapChart (@var{parent}, @var{spec}, @var{args})
    ##
    ## Create a @code{stats.chart.HeatmapChart} object.
    ##
    ## @var{parent} is the figure or panel to place the chart in, or empty for
    ## the current figure, which is resolved only once every value has been
    ## accepted.  @var{spec} is a structure carrying the data as
    ## @code{heatmap} resolved it, with the fields @qcode{ColorData},
    ## @qcode{XData}, @qcode{YData}, @qcode{SourceTable}, @qcode{XVariable},
    ## @qcode{YVariable} and @qcode{ColorVariable}, and @var{args} the
    ## name-value pairs left to set.  The documented way to reach this
    ## constructor is @code{heatmap}.
    ##
    ## @seealso{heatmap}
    ## @end deftypefn
    function this = HeatmapChart (parent, spec, args)

      if (nargin < 2)
        error ("stats.chart.HeatmapChart: too few input arguments.");
      endif
      if (nargin < 3)
        args = {};
      endif

      this.Colormap = hmDefaultColormap ();
      if (isempty (spec.SourceTable))
        this.ColorData = spec.ColorData;
        if (! isempty (spec.XData))
          this.XData = spec.XData;
        endif
        if (! isempty (spec.YData))
          this.YData = spec.YData;
        endif
      else
        this.SourceTable = spec.SourceTable;
        this.XVariable = spec.XVariable;
        this.YVariable = spec.YVariable;
        this.ColorVariable = spec.ColorVariable;
      endif

      ## Whatever was named is set now, so that one drawing covers them all;
      ## a value set here counts as chosen, as it would once drawn
      placed = false;
      for k = 1:2:numel (args)
        name = args{k};
        this.(name) = args{k+1};
        switch (name)
          case {'XLimits', 'YLimits', 'ColorLimits', 'Title', 'XLabel', ...
                'YLabel', 'ColorMethod'}
            chosen (this, name);
          case {'Position', 'InnerPosition', 'OuterPosition'}
            placed = true;
        endswitch
      endfor
      aggregate (this);

      ## The parent comes last, so a rejected value leaves no figure behind.
      ## The chart takes the place of the current axes of its parent, as
      ## MATLAB's does, and refuses axes that are held
      if (isempty (parent))
        parent = gcf ();
      endif
      fig = ancestor (parent, 'figure');
      ca = get (fig, 'currentaxes');
      if (! isempty (ca) && get (ca, 'parent') == parent)
        if (ishold (ca))
          error (strcat ("stats.chart.HeatmapChart: a heatmap cannot be", ...
                         " added to axes on which hold is on."));
        endif
        if (! placed)
          units = get (ca, 'units');
          set (ca, 'units', this.Units);
          this.OuterPosition = get (ca, 'outerposition');
          set (ca, 'units', units);
        endif
        delete (ca);
      endif
      this.Parent = parent;
      this.Axes_ = axes ('parent', parent, 'tag', 'stats.chart.HeatmapChart');
      set (fig, 'currentaxes', this.Axes_);
      ## The chart goes with its axes, and its axes with the chart
      set (this.Axes_, 'deletefcn', @(~, ~) delete (this));
      this.Drawn_ = true;
      redraw (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {stats.chart.HeatmapChart} {} sortx (@var{h})
    ## @deftypefnx {stats.chart.HeatmapChart} {} sortx (@var{h}, @var{row})
    ## @deftypefnx {stats.chart.HeatmapChart} {} sortx (@var{h}, @var{row}, @var{direction})
    ## @deftypefnx {stats.chart.HeatmapChart} {} sortx (@dots{}, @qcode{'MissingPlacement'}, @var{mp})
    ##
    ## Reorder the columns of a heatmap.
    ##
    ## @code{sortx (@var{h})} puts the columns in the order of their names.
    ## @code{sortx (@var{h}, @var{row})} puts them in the order of their values
    ## in the row named @var{row}, an element of @code{YData}; a cell array of
    ## rows sorts by the first, then by the next where the first ties.
    ## @var{direction} is @qcode{'ascend'}, the default, or @qcode{'descend'}.
    ## Missing values go last when ascending and first when descending,
    ## unless @var{mp}, @qcode{'first'}, @qcode{'last'} or @qcode{'auto'}, says
    ## otherwise.  The new order is set as @code{XDisplayData}.
    ##
    ## @end deftypefn
    function sortx (this, varargin)
      this.XDisplayData = hmSortOrder (this, varargin, true);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {stats.chart.HeatmapChart} {} sorty (@var{h})
    ## @deftypefnx {stats.chart.HeatmapChart} {} sorty (@var{h}, @var{column})
    ## @deftypefnx {stats.chart.HeatmapChart} {} sorty (@var{h}, @var{column}, @var{direction})
    ## @deftypefnx {stats.chart.HeatmapChart} {} sorty (@dots{}, @qcode{'MissingPlacement'}, @var{mp})
    ##
    ## Reorder the rows of a heatmap, as @code{sortx} reorders its columns, by
    ## the values in the column named @var{column}.
    ##
    ## @end deftypefn
    function sorty (this, varargin)
      this.YDisplayData = hmSortOrder (this, varargin, false);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {stats.chart.HeatmapChart} {@var{limits} =} xlim (@var{h})
    ## @deftypefnx {stats.chart.HeatmapChart} {} xlim (@var{h}, @var{limits})
    ##
    ## The first and last column shown, as @code{XLimits} holds them; given
    ## @var{limits}, set them.
    ##
    ## @end deftypefn
    function out = xlim (this, lims)
      if (nargin > 1)
        this.XLimits = lims;
      else
        out = this.XLimits;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {stats.chart.HeatmapChart} {@var{limits} =} ylim (@var{h})
    ## @deftypefnx {stats.chart.HeatmapChart} {} ylim (@var{h}, @var{limits})
    ##
    ## The first and last row shown, as @code{YLimits} holds them; given
    ## @var{limits}, set them.
    ##
    ## @end deftypefn
    function out = ylim (this, lims)
      if (nargin > 1)
        this.YLimits = lims;
      else
        out = this.YLimits;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {stats.chart.HeatmapChart} {} delete (@var{h})
    ##
    ## Delete the heatmap and the axes it is drawn in.
    ##
    ## @end deftypefn
    function delete (this)
      ax = this.Axes_;
      cb = this.Colorbar_;
      this.Axes_ = [];
      this.Colorbar_ = [];
      this.Parent = [];
      if (! isempty (cb) && ishghandle (cb))
        delete (cb);
      endif
      if (! isempty (ax) && ishghandle (ax)
          && strcmp (get (ax, 'beingdeleted'), 'off'))
        set (ax, 'deletefcn', '');
        delete (ax);
      endif
    endfunction

  endmethods

  methods (Hidden)

    ## Let go of the axes when something other than the chart clears them,
    ## as a later plot does, leaving them to that plot.
    function release (this)
      if (this.Redrawing_ || isempty (this.Axes_))
        return;
      endif
      ax = this.Axes_;
      this.Axes_ = [];
      this.Parent = [];
      if (ishghandle (ax))
        set (ax, 'deletefcn', '', 'tag', '');
      endif
      if (! isempty (this.Colorbar_) && ishghandle (this.Colorbar_))
        delete (this.Colorbar_);
      endif
      this.Colorbar_ = [];
    endfunction

  endmethods

  methods (Access = private)

    ## Note a value the caller chose, so that later changes leave it alone.
    function chosen (this, name)
      switch (name)
        case 'XLimits'
          this.XLimitsAuto_ = false;
        case 'YLimits'
          this.YLimitsAuto_ = false;
        case 'ColorLimits'
          this.ColorLimitsAuto_ = false;
        case 'Title'
          this.TitleAuto_ = false;
        case 'XLabel'
          this.XLabelAuto_ = false;
        case 'YLabel'
          this.YLabelAuto_ = false;
        case 'ColorMethod'
          this.MethodAuto_ = false;
      endswitch
    endfunction

    ## Bring XData, YData and what follows them in line with ColorData.
    function fitData (this)
      [r, c] = size (this.ColorData);
      if (this.XDataAuto_)
        this.XData = hmNumbered (c);
        this.XDataAuto_ = true;
      endif
      if (this.YDataAuto_)
        this.YData = hmNumbered (r);
        this.YDataAuto_ = true;
      endif
      if (this.XDisplayAuto_)
        this.XDisplayData = this.XData;
        this.XDisplayAuto_ = true;
      endif
      if (this.YDisplayAuto_)
        this.YDisplayData = this.YData;
        this.YDisplayAuto_ = true;
      endif
    endfunction

    ## Compute ColorData, its names and the automatic labels from the table.
    function aggregate (this)

      t = this.SourceTable;
      if (this.Aggregating_ || isempty (t) || isempty (this.XVariable)
          || isempty (this.YVariable))
        return;
      endif
      this.Aggregating_ = true;
      unwind_protect
        aggregateTable (this, t);
      unwind_protect_cleanup
        this.Aggregating_ = false;
      end_unwind_protect
      redraw (this);

    endfunction

    ## Read the table into ColorData, its names and the automatic labels.
    function aggregateTable (this, t)

      names = t.Properties.VariableNames;
      for v = {'XVariable', 'YVariable', 'ColorVariable'}
        n = this.(v{1});
        if (! isempty (n) && ! any (strcmp (names, n)))
          error (strcat ("stats.chart.HeatmapChart: '%s' does not name a", ...
                         " variable of the source table."), v{1});
        endif
      endfor
      if (this.MethodAuto_)
        if (isempty (this.ColorVariable))
          this.ColorMethod = 'count';
        else
          this.ColorMethod = 'mean';
        endif
        this.MethodAuto_ = true;
      endif
      method = this.ColorMethod;

      [gx, xn] = hmGroups (t.(this.XVariable), 'XVariable');
      [gy, yn] = hmGroups (t.(this.YVariable), 'YVariable');
      keep = gx > 0 & gy > 0;
      if (strcmp (method, 'count') || isempty (this.ColorVariable))
        v = ones (numel (gx), 1);
        method = 'count';
      else
        v = t.(this.ColorVariable);
        if (! (isnumeric (v) || islogical (v)) || ! isreal (v))
          error (strcat ("stats.chart.HeatmapChart: 'ColorVariable' must", ...
                         " name a real numeric variable of the source", ...
                         " table."));
        endif
        v = double (v(:));
      endif
      C = hmAggregate (gy(keep), gx(keep), v(keep), numel (yn), ...
                       numel (xn), method);

      ## Written straight to the fields: the setters would refuse ColorData
      ## on a table, and each would redraw
      drawn = this.Drawn_;
      this.Drawn_ = false;
      unwind_protect
        this.ColorData = C;
        this.XData = xn;
        this.YData = yn;
        this.XDataAuto_ = false;
        this.YDataAuto_ = false;
        this.XDisplayAuto_ = true;
        this.YDisplayAuto_ = true;
        fitData (this);
        if (this.TitleAuto_)
          this.Title = hmTableTitle (method, this.XVariable, ...
                                     this.YVariable, this.ColorVariable);
        endif
        if (this.XLabelAuto_)
          this.XLabel = this.XVariable;
        endif
        if (this.YLabelAuto_)
          this.YLabel = this.YVariable;
        endif
      unwind_protect_cleanup
        this.Drawn_ = drawn;
      end_unwind_protect

    endfunction

    ## Draw the chart from scratch into its axes.
    function redraw (this)

      if (! this.Drawn_ || isempty (this.Axes_) || ! ishghandle (this.Axes_))
        return;
      endif
      ax = this.Axes_;
      this.Redrawing_ = true;
      unwind_protect
        if (! isempty (this.Colorbar_) && ishghandle (this.Colorbar_))
          delete (this.Colorbar_);
        endif
        this.Colorbar_ = [];
        delete (get (ax, 'children'));
      unwind_protect_cleanup
        this.Redrawing_ = false;
      end_unwind_protect
      ## A child whose deletion by anyone else releases the axes
      text (ax, NaN, NaN, '', 'deletefcn', @(~, ~) release (this));

      ## The columns and rows shown, within the limits
      xd = this.XDisplayData;
      yd = this.YDisplayData;
      xi = hmWithin (xd, this.XLimits);
      yi = hmWithin (yd, this.YLimits);
      C = this.ColorDisplayData;
      C = C(yi, xi);
      xl = this.XDisplayLabels;
      yl = this.YDisplayLabels;
      xl = xl(xi);
      yl = yl(yi);
      [nr, nc] = size (C);

      ## The values as colours, and the limits they follow until set
      V = hmScale (C, this.ColorScaling);
      if (this.ColorLimitsAuto_)
        this.Internal_ = true;
        unwind_protect
          this.ColorLimits = hmAutoLimits (V, this.ColorScaling);
        unwind_protect_cleanup
          this.Internal_ = false;
        end_unwind_protect
      endif

      set (ax, 'units', this.Units, 'position', this.Position, ...
           'visible', this.Visible, 'color', this.MissingDataColor, ...
           'ydir', 'reverse', 'xlim', [0.5, nc + 0.5], ...
           'ylim', [0.5, nr + 0.5], 'xtick', 1:nc, 'ytick', 1:nr, ...
           'xticklabel', xl, 'yticklabel', yl, 'ticklength', [0, 0], ...
           'box', 'on', 'layer', 'top', 'clim', this.ColorLimits, ...
           'fontname', this.FontName, 'fontsize', this.FontSize, ...
           'xcolor', this.FontColor, 'ycolor', this.FontColor, ...
           'ticklabelinterpreter', this.Interpreter);
      colormap (ax, this.Colormap);
      if (nr > 0 && nc > 0)
        image (ax, 'cdata', V, 'cdatamapping', 'scaled', ...
               'alphadata', double (isfinite (V)), 'xdata', [1, nc], ...
               'ydata', [1, nr], 'visible', this.Visible);
      endif

      ## Lines between the cells
      if (strcmp (this.GridVisible, 'on'))
        gx = [(1.5:nc)', (1.5:nc)'; zeros(0, 2)];
        gy = [(1.5:nr)', (1.5:nr)'; zeros(0, 2)];
        for k = 1:rows (gx)
          line (ax, gx(k,:), [0.5, nr + 0.5], 'color', [1, 1, 1], ...
                'visible', this.Visible);
        endfor
        for k = 1:rows (gy)
          line (ax, [0.5, nc + 0.5], gy(k,:), 'color', [1, 1, 1], ...
                'visible', this.Visible);
        endfor
      endif

      ## The value of each cell written into it
      if (! (ischar (this.CellLabelColor)
             && strcmp (this.CellLabelColor, 'none')))
        cmap = this.Colormap;
        cl = this.ColorLimits;
        for i = 1:nr
          for j = 1:nc
            if (isnan (C(i,j)))
              s = this.MissingDataLabel;
              bg = this.MissingDataColor;
            else
              s = sprintf (this.CellLabelFormat, C(i,j));
              bg = hmCellColor (V(i,j), cl, cmap, this.MissingDataColor);
            endif
            if (ischar (this.CellLabelColor))
              fg = [0, 0, 0];
              if (bg * [0.299; 0.587; 0.114] < 0.5)
                fg = [1, 1, 1];
              endif
            else
              fg = this.CellLabelColor;
            endif
            text (ax, j, i, s, 'horizontalalignment', 'center', ...
                  'verticalalignment', 'middle', 'color', fg, ...
                  'fontname', this.FontName, 'fontsize', this.FontSize, ...
                  'interpreter', this.Interpreter, 'visible', this.Visible);
          endfor
        endfor
      endif

      title (ax, this.Title, 'color', this.FontColor, 'interpreter', ...
             this.Interpreter, 'fontname', this.FontName);
      xlabel (ax, this.XLabel, 'color', this.FontColor, 'interpreter', ...
              this.Interpreter, 'fontname', this.FontName);
      ylabel (ax, this.YLabel, 'color', this.FontColor, 'interpreter', ...
              this.Interpreter, 'fontname', this.FontName);

      ## The cells placed within the outer position, leaving the labels the
      ## room they take and the colour bar its own, or where they were put
      showBar = strcmp (this.ColorbarVisible, 'on');
      if (strcmp (this.PositionConstraint, 'outerposition'))
        this.Internal_ = true;
        unwind_protect
          this.Position = hmFitInner (this.OuterPosition, ...
                                      get (ax, 'tightinset'), showBar);
        unwind_protect_cleanup
          this.Internal_ = false;
        end_unwind_protect
        set (ax, 'position', this.Position);
      else
        this.Internal_ = true;
        unwind_protect
          this.OuterPosition = hmOuterFromInner (this.Position);
        unwind_protect_cleanup
          this.Internal_ = false;
        end_unwind_protect
      endif

      ## The colour bar, placed beside the cells rather than taking room from
      ## them; it stays visible to handles, which printing needs
      if (showBar)
        p = this.Position;
        this.Colorbar_ = colorbar (ax, 'units', this.Units, 'position', ...
                                   [p(1) + p(3) + 0.02 * p(3), p(2), ...
                                    0.04 * p(3), p(4)]);
        set (this.Colorbar_, 'visible', this.Visible, 'fontname', ...
             this.FontName, 'fontsize', this.FontSize, 'ycolor', ...
             this.FontColor);
      endif

    endfunction

  endmethods

endclassdef

## The position of the cells within the outer position O: MATLAB's default
## margins, widened where the labels, whose extent TI gives, need more
function p = hmFitInner (o, ti, bar)
  pad = 0.02;
  L = max (0.13 * o(3), ti(1) + pad * o(3));
  B = max (0.11 * o(4), ti(2) + pad * o(4));
  T = max (0.075 * o(4), ti(4) + pad * o(4));
  if (bar)
    R = 0.137857 * o(3);
  else
    R = max (0.095 * o(3), ti(3) + pad * o(3));
  endif
  w = max (o(3) - L - R, 0.01 * o(3));
  h = max (o(4) - B - T, 0.01 * o(4));
  p = [o(1) + L, o(2) + B, w, h];
endfunction

## The outer position that holds the cells at P
function o = hmOuterFromInner (p)
  o = [0, 0, p(3) / 0.732143, p(4) / 0.815];
  o(1) = p(1) - 0.13 * o(3);
  o(2) = p(2) - 0.11 * o(4);
endfunction

## Default colormap: 256 colours from a pale blue to the first line colour
function map = hmDefaultColormap ()
  a = [0.9, 0.9447, 0.9741];
  b = [0, 0.447, 0.741];
  map = a + (0:255)' / 255 .* (b - a);
endfunction

## The names '1', '2', ... up to N, as a column cell array
function out = hmNumbered (n)
  out = arrayfun (@(k) sprintf ('%d', k), (1:n)', 'UniformOutput', false);
endfunction

## A 2-D real numeric matrix
function v = hmCheckMatrix (val)
  if (! ((isnumeric (val) || islogical (val)) && isreal (val)
         && ndims (val) == 2))
    error (strcat ("stats.chart.HeatmapChart: 'ColorData' must be a", ...
                   " 2-dimensional real numeric matrix."));
  endif
  v = double (val);
endfunction

## Distinct names as a column cell array of character vectors, from a cell
## array of text, a string array, a categorical array or numbers
function out = hmCheckNames (val, name)
  out = hmText (val);
  if (isempty (out) && ! isempty (val))
    error (strcat ("stats.chart.HeatmapChart: '%s' must be a vector of", ...
                   " text, numbers or categories."), name);
  endif
  if (numel (unique (out)) != numel (out))
    error ("stats.chart.HeatmapChart: '%s' holds duplicate values.", name);
  endif
endfunction

## Values of any supported kind as a column cell array of text; numbers are
## written with %g, as MATLAB writes them
function out = hmText (val)
  out = cell (0, 1);
  if (isempty (val))
    return;
  elseif (iscellstr (val))
    out = val(:);
  elseif (isa (val, 'string'))
    out = cellstr (val(:));
  elseif (isa (val, 'categorical'))
    out = cellstr (val(:));
  elseif ((isnumeric (val) || islogical (val)) && isvector (val))
    out = arrayfun (@(v) sprintf ('%g', v), double (val(:)), ...
                    'UniformOutput', false);
  elseif (ischar (val) && isrow (val))
    out = {val};
  endif
endfunction

## Labels: text, one for each of N values
function out = hmCheckLabels (val, n, name)
  out = hmText (val);
  if (numel (out) != n)
    error (strcat ("stats.chart.HeatmapChart: '%s' must hold one label", ...
                   " for each displayed value."), name);
  endif
endfunction

## The labels of the displayed values: those set by hand, else the values
function out = hmLabels (shown, map)
  out = shown;
  for k = 1:numel (shown)
    i = find (strcmp (map(:,1), shown{k}), 1);
    if (! isempty (i))
      out{k} = map{i,2};
    endif
  endfor
endfunction

## Record LABELS against the displayed values DISP
function map = hmSetLabels (map, shown, labels)
  for k = 1:numel (shown)
    i = find (strcmp (map(:,1), shown{k}), 1);
    if (isempty (i))
      map(end+1,:) = {shown{k}, labels{k}};
    else
      map{i,2} = labels{k};
    endif
  endfor
endfunction

## Limits: a pair of displayed values, the first not after the second
function out = hmCheckLimits (val, shown, name)
  if (isempty (val) || isempty (shown))
    out = hmFullLimits (shown);
    return;
  endif
  val = hmText (val);
  if (numel (val) != 2)
    error ("stats.chart.HeatmapChart: '%s' must hold two values.", name);
  endif
  i = find (strcmp (shown, val{1}), 1);
  j = find (strcmp (shown, val{2}), 1);
  if (isempty (i) || isempty (j) || i > j)
    error (strcat ("stats.chart.HeatmapChart: '%s' must name two", ...
                   " displayed values, the first not after the second."), ...
           name);
  endif
  out = val(:)';
endfunction

## The first and last of the displayed values
function out = hmFullLimits (shown)
  if (isempty (shown))
    out = cell (1, 0);
  else
    out = {shown{1}, shown{end}};
  endif
endfunction

## Limits kept where they still fit the displayed values, reset otherwise,
## with a warning where they had been chosen
function out = hmFitLimits (lims, shown, auto, name)
  full = hmFullLimits (shown);
  if (auto || isempty (lims))
    out = full;
    return;
  endif
  i = find (strcmp (shown, lims{1}), 1);
  j = find (strcmp (shown, lims{2}), 1);
  if (isempty (i) || isempty (j) || i > j)
    warning (strcat ("stats.chart.HeatmapChart: the '%s' value does not", ...
                     " fit the new display order and is reset to its", ...
                     " full range."), name);
    out = full;
  else
    out = lims;
  endif
endfunction

## The indices of the displayed values from the first limit to the second
function idx = hmWithin (shown, lims)
  idx = 1:numel (shown);
  if (numel (lims) == 2)
    i = find (strcmp (shown, lims{1}), 1);
    j = find (strcmp (shown, lims{2}), 1);
    if (! isempty (i) && ! isempty (j))
      idx = i:j;
    endif
  endif
endfunction

## ColorData rearranged to the displayed columns and rows
function M = hmDisplayMatrix (C, xd, yd, xs, ys)
  M = NaN (numel (ys), numel (xs));
  [okx, jx] = ismember (xs, xd);
  [oky, jy] = ismember (ys, yd);
  if (! isempty (C))
    M(oky, okx) = C(jy(oky), jx(okx));
  endif
endfunction

## The values as colours
function V = hmScale (C, scaling)
  switch (scaling)
    case 'scaledcolumns'
      lo = min (C, [], 1);
      span = max (C, [], 1) - lo;
      span(span == 0) = 1;
      V = (C - lo) ./ span;
    case 'scaledrows'
      lo = min (C, [], 2);
      span = max (C, [], 2) - lo;
      span(span == 0) = 1;
      V = (C - lo) ./ span;
    case 'log'
      if (any (C(:) < 0))
        warning ("stats.chart.HeatmapChart: negative color data ignored.");
      endif
      V = C;
      V(V < 0) = NaN;
      V = log (V);
    otherwise
      V = C;
  endswitch
endfunction

## The colour limits the values shown call for
function cl = hmAutoLimits (V, scaling)
  if (any (strcmp (scaling, {'scaledcolumns', 'scaledrows'})))
    cl = [0, 1];
    return;
  endif
  v = V(isfinite (V));
  if (isempty (v))
    cl = [0, 1];
  elseif (min (v) == max (v))
    cl = min (v) + [-1, 1];
  else
    cl = [min(v), max(v)];
  endif
endfunction

## The colour a value takes, for the colour of its label
function c = hmCellColor (v, cl, cmap, missing)
  if (! isfinite (v))
    c = missing;
    return;
  endif
  k = round ((v - cl(1)) / (cl(2) - cl(1)) * (rows (cmap) - 1)) + 1;
  c = cmap(min (max (k, 1), rows (cmap)),:);
endfunction

## Variable names: a character vector or a string scalar
function out = hmCheckVarName (val, name)
  if (isa (val, 'string') && isscalar (val))
    val = char (val);
  endif
  if (! (ischar (val) && (isrow (val) || isempty (val))))
    error (strcat ("stats.chart.HeatmapChart: '%s' must be a character", ...
                   " vector."), name);
  endif
  out = val;
endfunction

## Text: a character vector, a string, or a cell array of them for lines
function out = hmCheckText (val, name)
  if (isa (val, 'string'))
    val = cellstr (val);
    if (numel (val) == 1)
      val = val{1};
    endif
  endif
  if (! ((ischar (val) && (isrow (val) || isempty (val))) || iscellstr (val)))
    error ("stats.chart.HeatmapChart: '%s' must be text.", name);
  endif
  out = val;
  if (iscell (out))
    out = out(:);
  endif
endfunction

## One of a list of names, case free
function out = hmCheckOneOf (val, list, name)
  if (isa (val, 'string') && isscalar (val))
    val = char (val);
  endif
  if (! (ischar (val) && any (strcmpi (val, list))))
    error ("stats.chart.HeatmapChart: '%s' must be one of %s.", name, ...
           strjoin (strcat ("'", list, "'"), ', '));
  endif
  out = lower (val);
endfunction

## A colour, as an RGB triplet or a name; 'auto' and 'none' too if AUTO
function out = hmCheckColor (val, name, auto)
  if (isa (val, 'string') && isscalar (val))
    val = char (val);
  endif
  if (auto && ischar (val) && any (strcmpi (val, {'auto', 'none'})))
    out = lower (val);
    return;
  endif
  if (ischar (val))
    names = {'red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'black', ...
             'white', 'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'};
    rgb = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0 0 0; 1 1 1];
    i = find (strcmpi (val, names), 1);
    if (! isempty (i))
      out = rgb(mod (i - 1, 8) + 1,:);
      return;
    endif
  elseif (isnumeric (val) && isreal (val) && numel (val) == 3
          && all (val(:) >= 0) && all (val(:) <= 1))
    out = double (val(:)');
    return;
  endif
  error (strcat ("stats.chart.HeatmapChart: '%s' must be an RGB", ...
                 " triplet or a colour name."), name);
endfunction

## A position, as a 1-by-4 vector with positive width and height
function out = hmCheckPosition (val, name)
  if (! (isnumeric (val) && isreal (val) && numel (val) == 4
         && all (isfinite (val)) && val(3) > 0 && val(4) > 0))
    error (strcat ("stats.chart.HeatmapChart: '%s' must be a 1-by-4", ...
                   " vector with positive width and height."), name);
  endif
  out = double (val(:)');
endfunction

## The group of each row of a table variable and the names of the groups,
## 0 for a missing value: categories in their order, unused ones included;
## numbers, text and logical values sorted
function [g, names] = hmGroups (v, name)
  v = v(:);
  if (isa (v, 'categorical'))
    names = categories (v);
    [~, g] = ismember (cellstr (v), names);
    g(isundefined (v)) = 0;
  elseif (isnumeric (v) || islogical (v))
    if (islogical (v))
      u = unique (v);
      names = cell (numel (u), 1);
      names(u == 0) = {'false'};
      names(u == 1) = {'true'};
      [~, g] = ismember (v, u);
    else
      ok = ! isnan (v);
      u = unique (v(ok));
      names = arrayfun (@(x) sprintf ('%g', x), u, 'UniformOutput', false);
      g = zeros (numel (v), 1);
      [~, g(ok)] = ismember (v(ok), u);
    endif
  elseif (iscellstr (v) || isa (v, 'string'))
    t = cellstr (v);
    ok = ! cellfun (@isempty, t);
    if (isa (v, 'string'))
      ok &= ! ismissing (v);
    endif
    names = unique (t(ok));
    g = zeros (numel (t), 1);
    [~, g(ok)] = ismember (t(ok), names);
  else
    error (strcat ("stats.chart.HeatmapChart: '%s' must name a variable", ...
                   " of categories, numbers, text or logical values."), name);
  endif
  names = names(:);
endfunction

## The values V of the rows falling in each cell, combined by METHOD
function C = hmAggregate (gy, gx, v, nr, nc, method)
  if (strcmp (method, 'count'))
    C = accumarray ([gy, gx], 1, [nr, nc]);
    return;
  endif
  if (strcmp (method, 'none'))
    n = accumarray ([gy, gx], 1, [nr, nc]);
    if (any (n(:) > 1))
      error (strcat ("stats.chart.HeatmapChart: 'ColorMethod' 'none'", ...
                     " cannot combine the several rows of the table that", ...
                     " fall in one cell."));
    endif
    C = NaN (nr, nc);
    C(sub2ind ([nr, nc], gy, gx)) = v;
    return;
  endif
  ## Missing values left out; a cell with none left is missing, but sums to
  ## zero, as MATLAB sums it
  ok = ! isnan (v);
  fill = NaN;
  switch (method)
    case 'mean'
      f = @(x) mean (x);
    case 'median'
      f = @(x) median (x);
    case 'sum'
      f = @(x) sum (x);
      fill = 0;
    case 'min'
      f = @(x) min (x);
    case 'max'
      f = @(x) max (x);
  endswitch
  C = accumarray ([gy(ok), gx(ok)], v(ok), [nr, nc], f, fill);
endfunction

## The title a table chart takes until one is set
function t = hmTableTitle (method, xv, yv, cv)
  switch (method)
    case 'count'
      t = sprintf ('Count of %s vs. %s', yv, xv);
    case 'mean'
      t = sprintf ('Mean of %s', cv);
    case 'median'
      t = sprintf ('Median of %s', cv);
    case 'sum'
      t = sprintf ('Sum of %s', cv);
    case 'min'
      t = sprintf ('Minimum of %s', cv);
    case 'max'
      t = sprintf ('Maximum of %s', cv);
    otherwise
      t = cv;
  endswitch
endfunction

## The new display order sortx and sorty ask for
function order = hmSortOrder (this, args, byColumns)

  if (byColumns)
    shown = this.XDisplayData;
    other = this.YData;
    otherDisp = this.YDisplayData;
    fname = 'sortx';
  else
    shown = this.YDisplayData;
    other = this.XData;
    otherDisp = this.XDisplayData;
    fname = 'sorty';
  endif
  mp = 'auto';
  isMP = @(a) (ischar (a) || isa (a, 'string')) ...
              && strcmpi (a, 'MissingPlacement');
  k = find (cellfun (isMP, args), 1);
  if (! isempty (k))
    if (k == numel (args))
      error ("%s: 'MissingPlacement' needs a value.", fname);
    endif
    mp = hmCheckOneOf (args{k+1}, {'auto', 'first', 'last'}, ...
                       'MissingPlacement');
    args(k:k+1) = [];
  endif
  if (numel (args) > 2)
    error ("%s: too many input arguments.", fname);
  endif
  direction = 'ascend';
  if (numel (args) > 1)
    direction = hmCheckOneOf (args{2}, {'ascend', 'descend'}, 'direction');
  endif

  if (isempty (args) || isempty (args{1}))
    [~, i] = sort (shown);
    if (strcmp (direction, 'descend'))
      i = flipud (i(:));
    endif
    order = shown(i);
    return;
  endif

  ## The values of the named rows (or columns), one column of keys each
  keys = hmText (args{1});
  if (isempty (keys) || ! all (ismember (keys, [other; otherDisp])))
    error (strcat ("%s: every element to sort by must be in 'XData' or", ...
                   " 'XDisplayData'."), fname);
  endif
  M = hmDisplayMatrix (this.ColorData, this.XData, this.YData, ...
                       this.XDisplayData, this.YDisplayData);
  if (byColumns)
    [~, r] = ismember (keys, this.YDisplayData);
    K = M(r,:)';
  else
    [~, c] = ismember (keys, this.XDisplayData);
    K = M(:,c);
  endif
  if (strcmp (mp, 'auto'))
    mp = 'last';
    if (strcmp (direction, 'descend'))
      mp = 'first';
    endif
  endif
  ## Missing values placed as asked, then sorted by each key in turn
  if (strcmp (direction, 'descend'))
    K = -K;
  endif
  if (strcmp (mp, 'first'))
    K(isnan (K)) = -Inf;
  else
    K(isnan (K)) = Inf;
  endif
  [~, i] = sortrows ([K, (1:rows (K))']);
  order = shown(i);

endfunction

%!shared tbl
%! tbl = table (categorical ({'a';'b';'a';'c';'b';'a';'c';'c'}), ...
%!              categorical ({'u';'u';'v';'v';'u';'u';'v';'u'}), ...
%!              [1;2;3;4;5;6;7;8], 'VariableNames', {'X', 'Y', 'V'});

%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap (magic (4));
%!   assert_equal (class (h), 'stats.chart.HeatmapChart');
%!   assert_equal (h.XData, {'1'; '2'; '3'; '4'});
%!   assert_equal (h.YData, {'1'; '2'; '3'; '4'});
%!   assert_equal (h.ColorData, magic (4));
%!   assert_equal (h.ColorDisplayData, magic (4));
%!   assert_equal (h.ColorLimits, [1, 16]);
%!   assert_equal (h.ColorMethod, 'none');
%!   assert_equal (h.ColorScaling, 'scaled');
%!   assert_equal (h.CellLabelFormat, '%0.4g');
%!   assert_equal (h.CellLabelColor, 'auto');
%!   assert_equal (h.MissingDataColor, [0.15, 0.15, 0.15]);
%!   assert_equal (h.MissingDataLabel, 'NaN');
%!   assert_equal (h.FontSize, 10);
%!   assert_equal (h.Position, [0.13, 0.11, 0.732143, 0.815], -1e-12);
%!   assert_equal (h.OuterPosition, [0, 0, 1, 1]);
%!   assert_equal (h.Parent, hf);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap (magic (3));
%!   assert_equal (size (h.Colormap), [256, 3]);
%!   assert_equal (h.Colormap(1,:), [0.9, 0.9447, 0.9741]);
%!   assert_equal (h.Colormap(2,:), ...
%!                 [0.896470588235294, 0.942748235294118, ...
%!                  0.973185882352941], -1e-14);
%!   assert_equal (h.Colormap(end,:), [0, 0.447, 0.741], -1e-15);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ({'a', 'b', 'c'}, {'r1', 'r2'}, [1, 2, 3; 4, 5, 6]);
%!   assert_equal (h.XData, {'a'; 'b'; 'c'});
%!   assert_equal (h.YData, {'r1'; 'r2'});
%!   assert_equal (h.XDisplayLabels, {'a'; 'b'; 'c'});
%!   assert_equal (h.XLimits, {'a', 'c'});
%!   assert_equal (h.ColorLimits, [1, 6]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ([1, 2, 3], [10, 20], [1, 2, 3; 4, 5, 6]);
%!   assert_equal (h.XData, {'1'; '2'; '3'});
%!   assert_equal (h.YData, {'10'; '20'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap (tbl, 'X', 'Y');
%!   assert_equal (h.ColorData, [2, 2, 1; 1, 0, 2]);
%!   assert_equal (h.XData, {'a'; 'b'; 'c'});
%!   assert_equal (h.YData, {'u'; 'v'});
%!   assert_equal (h.Title, 'Count of Y vs. X');
%!   assert_equal (h.XLabel, 'X');
%!   assert_equal (h.YLabel, 'Y');
%!   assert_equal (h.ColorMethod, 'count');
%!   assert_equal (h.ColorVariable, '');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap (tbl, 'X', 'Y', 'ColorVariable', 'V');
%!   assert_equal (h.ColorData, [3.5, 3.5, 8; 3, NaN, 5.5]);
%!   assert_equal (h.Title, 'Mean of V');
%!   assert_equal (h.ColorMethod, 'mean');
%!   h.ColorMethod = 'sum';
%!   assert_equal (h.ColorData, [7, 7, 8; 3, 0, 11]);
%!   assert_equal (h.Title, 'Sum of V');
%!   h.ColorMethod = 'median';
%!   assert_equal (h.ColorData, [3.5, 3.5, 8; 3, NaN, 5.5]);
%!   assert_equal (h.Title, 'Median of V');
%!   h.ColorMethod = 'max';
%!   assert_equal (h.ColorData, [6, 5, 8; 3, NaN, 7]);
%!   assert_equal (h.Title, 'Maximum of V');
%!   h.ColorMethod = 'min';
%!   assert_equal (h.Title, 'Minimum of V');
%!   h.ColorMethod = 'count';
%!   assert_equal (h.Title, 'Count of Y vs. X');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   t = table (categorical ({'a';'b';'c';'a'}), ...
%!              categorical ({'u';'u';'v';'v'}), [1;2;3;4], ...
%!              'VariableNames', {'X', 'Y', 'V'});
%!   h = heatmap (t, 'X', 'Y', 'ColorVariable', 'V', 'ColorMethod', 'none');
%!   assert_equal (h.ColorData, [1, 2, NaN; 4, NaN, 3]);
%!   assert_equal (h.Title, 'V');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   t = table ({'a';'a';'a';'b';'b'}, {'u';'u';'u';'u';'u'}, ...
%!              [1;NaN;5;2;4], 'VariableNames', {'X', 'Y', 'V'});
%!   h = heatmap (t, 'X', 'Y', 'ColorVariable', 'V');
%!   assert_equal (h.ColorData, [3, 3]);
%!   h.ColorMethod = 'count';
%!   assert_equal (h.ColorData, [3, 2]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   t = table (categorical ({'lo';'hi';'lo'}, {'lo', 'mid', 'hi'}), ...
%!              categorical ({'u';'u';'v'}), [1;2;3], ...
%!              'VariableNames', {'X', 'Y', 'V'});
%!   h = heatmap (t, 'X', 'Y');
%!   assert_equal (h.XData, {'lo'; 'mid'; 'hi'});
%!   assert_equal (h.ColorData, [1, 0, 1; 1, 0, 0]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   t = table ([0.5; -2; 1e6; NaN; 3], {'u';'v';'u';'u';''}, ...
%!              [1;2;3;4;5], 'VariableNames', {'X', 'Y', 'V'});
%!   h = heatmap (t, 'X', 'Y');
%!   assert_equal (h.XData, {'-2'; '0.5'; '3'; '1e+06'});
%!   assert_equal (h.YData, {'u'; 'v'});
%!   assert_equal (h.ColorData, [0, 1, 0, 1; 1, 0, 0, 0]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   t = table ([3;1;2;1;3], {'p';'q';'p';'p';'q'}, [true; false; true; ...
%!              true; false], 'VariableNames', {'N', 'S', 'L'});
%!   h = heatmap (t, 'N', 'S');
%!   assert_equal (h.XData, {'1'; '2'; '3'});
%!   assert_equal (h.ColorData, [1, 1, 1; 1, 0, 1]);
%!   h.XVariable = 'L';
%!   assert_equal (h.XData, {'false'; 'true'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ([1, NaN, 3; 4, 5, NaN]);
%!   assert_equal (h.ColorLimits, [1, 5]);
%!   assert_equal (h.ColorDisplayData, [1, NaN, 3; 4, 5, NaN]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ([1, 2, 3; 4, 5, 16]);
%!   h.ColorScaling = 'scaledcolumns';
%!   assert_equal (h.ColorDisplayData, [1, 2, 3; 4, 5, 16]);
%!   assert_equal (h.ColorLimits, [0, 1]);
%!   h.ColorScaling = 'scaledrows';
%!   assert_equal (h.ColorLimits, [0, 1]);
%!   h.ColorScaling = 'log';
%!   assert_equal (h.ColorLimits, [0, log(16)], -1e-15);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ([2, 2; 2, 2]);
%!   assert_equal (h.ColorLimits, [1, 3]);
%!   h.ColorData = [NaN, NaN; NaN, NaN];
%!   assert_equal (h.ColorLimits, [0, 1]);
%!   h.ColorData = [5, 6; 7, 50];
%!   assert_equal (h.ColorLimits, [5, 50]);
%!   h.ColorLimits = [0, 20];
%!   h.ColorData = [1, 100; 2, 3];
%!   assert_equal (h.ColorLimits, [0, 20]);
%!   h.ColorScaling = 'scaledcolumns';
%!   assert_equal (h.ColorLimits, [0, 20]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!warning<stats.chart.HeatmapChart: negative color data ignored.> ...
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ([1, -2, 3; 4, 5, 6]);
%!   h.ColorScaling = 'log';
%!   assert_equal (h.ColorLimits, [0, log(6)], -1e-15);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ({'c', 'a', 'b'}, {'y', 'x'}, [3, 1, 2; 6, 4, 5]);
%!   sortx (h);
%!   assert_equal (h.XDisplayData, {'a'; 'b'; 'c'});
%!   assert_equal (h.ColorDisplayData, [1, 2, 3; 4, 5, 6]);
%!   sortx (h, 'x', 'descend');
%!   assert_equal (h.XDisplayData, {'c'; 'b'; 'a'});
%!   assert_equal (h.ColorDisplayData, [3, 2, 1; 6, 5, 4]);
%!   sorty (h, 'a');
%!   assert_equal (h.YDisplayData, {'y'; 'x'});
%!   sorty (h);
%!   assert_equal (h.YDisplayData, {'x'; 'y'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ({'c', 'a', 'b'}, {'y', 'x'}, [3, NaN, 2; 6, 4, 5]);
%!   sortx (h, 'y');
%!   assert_equal (h.XDisplayData, {'b'; 'c'; 'a'});
%!   sortx (h, 'y', 'descend');
%!   assert_equal (h.XDisplayData, {'a'; 'c'; 'b'});
%!   sortx (h, 'y', 'ascend', 'MissingPlacement', 'first');
%!   assert_equal (h.XDisplayData, {'a'; 'b'; 'c'});
%!   sortx (h, {'y', 'x'});
%!   assert_equal (h.XDisplayData, {'b'; 'c'; 'a'});
%!   sorty (h, 'b', 'descend');
%!   assert_equal (h.YDisplayData, {'x'; 'y'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ({'a', 'b', 'c'}, {'r1', 'r2'}, [1, 2, 3; 4, 5, 6]);
%!   h.XDisplayData = {'c'; 'a'};
%!   assert_equal (h.ColorDisplayData, [3, 1; 6, 4]);
%!   h.XDisplayLabels = {'C'; 'A'};
%!   assert_equal (h.XDisplayLabels, {'C'; 'A'});
%!   assert_equal (h.XDisplayData, {'c'; 'a'});
%!   h.XDisplayData = {'a'; 'c'; 'b'};
%!   assert_equal (h.XDisplayLabels, {'A'; 'C'; 'b'});
%!   h.XDisplayData = {'9'};
%!   assert_equal (h.ColorDisplayData, [NaN; NaN]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ({'c', 'a', 'b'}, {'y', 'x'}, [3, 1, 2; 6, 4, 5]);
%!   assert_equal (xlim (h), {'c', 'b'});
%!   xlim (h, {'a', 'b'});
%!   assert_equal (h.XLimits, {'a', 'b'});
%!   assert_equal (h.XDisplayData, {'c'; 'a'; 'b'});
%!   assert_equal (ylim (h), {'y', 'x'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!warning<stats.chart.HeatmapChart: the 'XLimits' value does not fit the new display order and is reset to its full range.> ...
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ({'c', 'a', 'b'}, {'y', 'x'}, [3, 1, 2; 6, 4, 5]);
%!   xlim (h, {'a', 'b'});
%!   h.XDisplayData = {'b'; 'c'; 'a'};
%!   assert_equal (h.XLimits, {'b', 'a'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ({'c', 'a', 'b'}, {'y', 'x'}, [3, 1, 2; 6, 4, 5]);
%!   sortx (h);
%!   sorty (h);
%!   lastwarn ('');
%!   sortx (h, 'x', 'descend');
%!   sorty (h, 'a');
%!   assert_equal (isempty (strfind (lastwarn (), 'HeatmapChart')), true);
%!   assert_equal (h.XLimits, {'c', 'a'});
%!   assert_equal (h.YLimits, {'y', 'x'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap ([1, 2; 3, 4]);
%!   h.ColorData = [1, 2, 3; 4, 5, 6];
%!   assert_equal (h.XData, {'1'; '2'; '3'});
%!   h.XData = {'p'; 'q'; 'r'};
%!   assert_equal (h.XDisplayData, {'p'; 'q'; 'r'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap (tbl, 'X', 'Y', 'ColorVariable', 'V');
%!   h.Title = 'Mine';
%!   h.ColorMethod = 'sum';
%!   assert_equal (h.Title, 'Mine');
%!   h.XLabel = 'XX';
%!   h.XVariable = 'V';
%!   assert_equal (h.XLabel, 'XX');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap (magic (3), 'Colormap', gray (5), 'ColorbarVisible', 'off', ...
%!                'CellLabelColor', 'none', 'GridVisible', 'off', ...
%!                'Title', 'T', 'XLabel', 'XL', 'YLabel', 'YL', ...
%!                'FontSize', 14, 'CellLabelFormat', '%.2f');
%!   assert_equal (h.Colormap, gray (5));
%!   assert_equal (h.ColorbarVisible, 'off');
%!   assert_equal (h.CellLabelColor, 'none');
%!   assert_equal (h.GridVisible, 'off');
%!   assert_equal (h.Title, 'T');
%!   assert_equal (h.FontSize, 14);
%!   assert_equal (h.CellLabelFormat, '%.2f');
%!   h.Title = {'one', 'two'};
%!   assert_equal (h.Title, {'one'; 'two'});
%!   h.MissingDataColor = 'r';
%!   assert_equal (h.MissingDataColor, [1, 0, 0]);
%!   h.FontColor = 'g';
%!   assert_equal (h.FontColor, [0, 1, 0]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = heatmap (magic (3));
%!   ax = findall (hf, 'type', 'axes', 'tag', 'stats.chart.HeatmapChart');
%!   assert_equal (numel (ax), 1);
%!   assert_equal (gca (), ax);
%!   im = findall (ax, 'type', 'image');
%!   assert_equal (get (im, 'cdata'), magic (3));
%!   assert_equal (get (ax, 'clim'), [1, 9]);
%!   assert_equal (numel (findall (ax, 'type', 'text', 'string', '5')), 1);
%!   delete (h);
%!   tag = 'stats.chart.HeatmapChart';
%!   assert_equal (isempty (findall (hf, 'tag', tag)), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   hp = uipanel (hf);
%!   h = heatmap (hp, magic (3));
%!   assert_equal (h.Parent, hp);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!shared hmS, hmT
%! hmS = struct ('ColorData', magic (3), 'XData', {{}}, 'YData', {{}}, ...
%!               'SourceTable', [], 'XVariable', '', 'YVariable', '', ...
%!               'ColorVariable', '');
%! hmT = struct ('ColorData', [], 'XData', {{}}, 'YData', {{}}, ...
%!               'SourceTable', table ({'a'; 'a'}, {'u'; 'u'}, [1; 2], ...
%!                                     'VariableNames', {'X', 'Y', 'V'}), ...
%!               'XVariable', 'X', 'YVariable', 'Y', 'ColorVariable', 'V');
%!error<stats.chart.HeatmapChart: too few input arguments.> ...
%! stats.chart.HeatmapChart ([])
%!error<stats.chart.HeatmapChart: 'ColorData' must be a 2-dimensional real numeric matrix.> ...
%! stats.chart.HeatmapChart ([], hmS, {'ColorData', ones(2, 2, 2)})
%!error<stats.chart.HeatmapChart: 'XData' holds duplicate values.> ...
%! stats.chart.HeatmapChart ([], hmS, {'XData', {'a', 'a', 'b'}})
%!error<stats.chart.HeatmapChart: 'XData' must be a vector of text, numbers or categories.> ...
%! stats.chart.HeatmapChart ([], hmS, {'XData', {1, 2}})
%!error<stats.chart.HeatmapChart: 'XDisplayLabels' must hold one label for each displayed value.> ...
%! stats.chart.HeatmapChart ([], hmS, {'XDisplayLabels', {'a'}})
%!error<stats.chart.HeatmapChart: 'XLimits' must hold two values.> ...
%! stats.chart.HeatmapChart ([], hmS, {'XLimits', {'1'}})
%!error<stats.chart.HeatmapChart: 'XLimits' must name two displayed values, the first not after the second.> ...
%! stats.chart.HeatmapChart ([], hmS, {'XLimits', {'3', '1'}})
%!error<stats.chart.HeatmapChart: 'SourceTable' must be a table.> ...
%! stats.chart.HeatmapChart ([], hmS, {'SourceTable', 5})
%!error<stats.chart.HeatmapChart: 'ColorVariable' must be a character vector.> ...
%! stats.chart.HeatmapChart ([], hmS, {'ColorVariable', 5})
%!error<stats.chart.HeatmapChart: 'ColorMethod' must be one of 'count', 'mean', 'median', 'sum', 'min', 'max', 'none'.> ...
%! stats.chart.HeatmapChart ([], hmS, {'ColorMethod', 'mode'})
%!error<stats.chart.HeatmapChart: 'ColorScaling' must be one of 'scaled', 'scaledcolumns', 'scaledrows', 'log'.> ...
%! stats.chart.HeatmapChart ([], hmS, {'ColorScaling', 'foo'})
%!error<stats.chart.HeatmapChart: 'ColorLimits' must be a 1-by-2 increasing numeric vector.> ...
%! stats.chart.HeatmapChart ([], hmS, {'ColorLimits', [5, 1]})
%!error<stats.chart.HeatmapChart: 'Colormap' must be an M-by-3 matrix of values between 0 and 1.> ...
%! stats.chart.HeatmapChart ([], hmS, {'Colormap', 'jet'})
%!error<stats.chart.HeatmapChart: 'ColorbarVisible' must be one of 'on', 'off'.> ...
%! stats.chart.HeatmapChart ([], hmS, {'ColorbarVisible', 'yes'})
%!error<stats.chart.HeatmapChart: 'MissingDataColor' must be an RGB triplet or a colour name.> ...
%! stats.chart.HeatmapChart ([], hmS, {'MissingDataColor', 'auto'})
%!error<stats.chart.HeatmapChart: 'MissingDataLabel' must be text.> ...
%! stats.chart.HeatmapChart ([], hmS, {'MissingDataLabel', 5})
%!error<stats.chart.HeatmapChart: 'CellLabelFormat' must be a character vector.> ...
%! stats.chart.HeatmapChart ([], hmS, {'CellLabelFormat', 5})
%!error<stats.chart.HeatmapChart: 'CellLabelColor' must be an RGB triplet or a colour name.> ...
%! stats.chart.HeatmapChart ([], hmS, {'CellLabelColor', [2, 0, 0]})
%!error<stats.chart.HeatmapChart: 'GridVisible' must be one of 'on', 'off'.> ...
%! stats.chart.HeatmapChart ([], hmS, {'GridVisible', 1})
%!error<stats.chart.HeatmapChart: 'Title' must be text.> ...
%! stats.chart.HeatmapChart ([], hmS, {'Title', 5})
%!error<stats.chart.HeatmapChart: 'FontSize' must be a positive number.> ...
%! stats.chart.HeatmapChart ([], hmS, {'FontSize', 0})
%!error<stats.chart.HeatmapChart: 'FontColor' must be an RGB triplet or a colour name.> ...
%! stats.chart.HeatmapChart ([], hmS, {'FontColor', 'none'})
%!error<stats.chart.HeatmapChart: 'Interpreter' must be one of 'tex', 'latex', 'none'.> ...
%! stats.chart.HeatmapChart ([], hmS, {'Interpreter', 'html'})
%!error<stats.chart.HeatmapChart: 'Position' must be a 1-by-4 vector with positive width and height.> ...
%! stats.chart.HeatmapChart ([], hmS, {'Position', [0, 0, 0, 1]})
%!error<stats.chart.HeatmapChart: 'Units' must be one of 'normalized', 'inches', 'centimeters', 'points', 'pixels', 'characters'.> ...
%! stats.chart.HeatmapChart ([], hmS, {'Units', 'miles'})
%!error<stats.chart.HeatmapChart: 'Visible' must be one of 'on', 'off'.> ...
%! stats.chart.HeatmapChart ([], hmS, {'Visible', 'no'})
%!error<stats.chart.HeatmapChart: 'XVariable' does not name a variable of the source table.> ...
%! stats.chart.HeatmapChart ([], hmT, {'XVariable', 'Q'})
%!error<stats.chart.HeatmapChart: 'ColorMethod' 'none' cannot combine the several rows of the table that fall in one cell.> ...
%! stats.chart.HeatmapChart ([], hmT, {'ColorMethod', 'none'})
%!error<stats.chart.HeatmapChart: 'ColorVariable' must name a real numeric variable of the source table.> ...
%! stats.chart.HeatmapChart ([], hmT, {'ColorVariable', 'X'})
