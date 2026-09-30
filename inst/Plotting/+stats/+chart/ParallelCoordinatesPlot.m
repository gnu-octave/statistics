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
## @deftp {statistics} stats.chart.ParallelCoordinatesPlot
##
## A parallel coordinates plot, as @code{parallelplot} draws it.
##
## A @code{stats.chart.ParallelCoordinatesPlot} holds observations of several
## variables, optionally split into groups, and the choices they are drawn
## with, and redraws itself whenever one of them is set.  It is what
## @code{parallelplot} returns, and the documented way to reach it.  Each
## variable is a vertical ruler, a coordinate, and each observation a line
## joining its values across the rulers.
##
## The observations come either from a numeric matrix, @code{Data}, whose
## columns @code{CoordinateData} chooses, or from a table,
## @code{SourceTable}, whose variables @code{CoordinateVariables} chooses,
## never from both.  A table variable of categories, text or logical values
## places its distinct values evenly along its ruler, whatever
## @code{DataNormalization} says, and @code{Jitter} spreads the lines passing
## through each value.  A missing value breaks its line.  The groups take
## the colours of the colour order, in the order in which each first
## appears in @code{GroupData}.
##
## A chart takes the place of the current axes of its parent, keeping their
## outer position, and refuses axes on which @code{hold} is on; a later plot
## into the figure takes the place of the chart in turn.  MATLAB's chart is
## a graphics object of its own; this one is a handle class drawing into
## axes of its own, made current, and it lets go of them when a later plot
## clears them.  Octave's lines carry no transparency, so @code{LineAlpha}
## blends each line's colour toward the background of the axes instead,
## which matches MATLAB where lines do not cross.  @code{LegendVisible} and
## @code{Visible} hold @qcode{'on'} or @qcode{'off'}, and @code{LineStyle}
## and @code{MarkerStyle} a character vector or a cell array of them, where
## MATLAB holds an on-off state and a string array.
##
## @seealso{parallelplot, parallelcoords}
## @end deftp

classdef ParallelCoordinatesPlot < handle

  properties (Access = public)

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} Data
    ##
    ## The observations, a numeric matrix with one row for each observation
    ## and one column for each variable; empty where the chart was drawn from
    ## a table.
    ##
    ## @end deftp
    Data = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} CoordinateData
    ##
    ## The columns of @code{Data} drawn as coordinates, in the order drawn,
    ## as indices or a logical vector; all of them by default.  Setting it on
    ## a chart drawn from a table is refused.
    ##
    ## @end deftp
    CoordinateData = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} SourceTable
    ##
    ## The table the chart was drawn from, or empty.
    ##
    ## @end deftp
    SourceTable = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} CoordinateVariables
    ##
    ## The variables of @code{SourceTable} drawn as coordinates, in the order
    ## drawn, as names, indices or a logical vector, kept in the form given;
    ## all of them by default.  Setting it on a chart drawn from a matrix is
    ## refused.
    ##
    ## @end deftp
    CoordinateVariables = {};

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} GroupData
    ##
    ## The group of each observation: a vector with one element for each row
    ## of @code{Data} of categories, text, numbers or logical values, or
    ## empty for no groups.
    ##
    ## @end deftp
    GroupData = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} GroupVariable
    ##
    ## The variable of @code{SourceTable} holding the groups, or empty.
    ##
    ## @end deftp
    GroupVariable = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} DataNormalization
    ##
    ## How the values are placed along the rulers: @qcode{'range'}, the
    ## default, as they are on rulers with limits of their own;
    ## @qcode{'none'} as they are on rulers sharing one scale;
    ## @qcode{'zscore'} as z-scores; @qcode{'scale'} divided by the standard
    ## deviation; @qcode{'center'} less the mean; @qcode{'norm'} divided by
    ## the 2-norm of the variable.  All but @qcode{'range'} share one scale.
    ##
    ## @end deftp
    DataNormalization = 'range';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} CoordinateTickLabels
    ##
    ## The label of each coordinate, a column cell array; by default the
    ## numbers of the columns or the names of the variables.
    ##
    ## @end deftp
    CoordinateTickLabels = cell (0, 1);

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} Jitter
    ##
    ## How far the lines through one value of a variable of categories are
    ## spread, from 0 to 1, where 1 fills the room up to the next value.  0.1
    ## by default.
    ##
    ## @end deftp
    Jitter = 0.1;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} Color
    ##
    ## The colour of each group, one row of an RGB matrix for each; by
    ## default the colours of the colour order.
    ##
    ## @end deftp
    Color = [0, 0.447, 0.741];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} LineStyle
    ##
    ## The line style of each group, a character vector or a cell array of
    ## them, @qcode{'-'} by default.
    ##
    ## @end deftp
    LineStyle = '-';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} LineWidth
    ##
    ## The width of the lines, a scalar or one for each group, 1 by default.
    ##
    ## @end deftp
    LineWidth = 1;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} LineAlpha
    ##
    ## How opaque the lines are, from 0 to 1, a scalar or one for each group,
    ## 0.7 by default.
    ##
    ## @end deftp
    LineAlpha = 0.7;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} MarkerStyle
    ##
    ## The marker at each value, for each group, @qcode{'none'} by default.
    ##
    ## @end deftp
    MarkerStyle = 'none';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} MarkerSize
    ##
    ## The size of the markers in points, a scalar or one for each group, 6
    ## by default.
    ##
    ## @end deftp
    MarkerSize = 6;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} LegendVisible
    ##
    ## Whether the legend of the groups is shown, @qcode{'on'} or
    ## @qcode{'off'}; on where there are groups, until set.
    ##
    ## @end deftp
    LegendVisible = 'off';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} LegendTitle
    ##
    ## The title of the legend; from a table the name of
    ## @code{GroupVariable}, until set.
    ##
    ## @end deftp
    LegendTitle = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} CoordinateLabel
    ##
    ## The label under the coordinates.
    ##
    ## @end deftp
    CoordinateLabel = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} DataLabel
    ##
    ## The label beside the rulers.
    ##
    ## @end deftp
    DataLabel = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} Title
    ##
    ## The title of the chart.
    ##
    ## @end deftp
    Title = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} FontName
    ##
    ## The font of every text on the chart, @qcode{'Helvetica'} by default.
    ##
    ## @end deftp
    FontName = 'Helvetica';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} FontSize
    ##
    ## The size of every text on the chart, in points, 10 by default.
    ##
    ## @end deftp
    FontSize = 10;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} Position
    ##
    ## The position of the rulers in their parent, as
    ## @code{[left, bottom, width, height]} in @code{Units}.  While
    ## @code{PositionConstraint} is @qcode{'outerposition'} the chart places
    ## them within @code{OuterPosition} itself; setting it fixes them there
    ## and turns @code{PositionConstraint} to @qcode{'innerposition'}.
    ## @code{InnerPosition} is the same.
    ##
    ## @end deftp
    Position = [0.13, 0.11, 0.775, 0.815];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} OuterPosition
    ##
    ## The position of the whole chart in its parent, labels included,
    ## @code{[0, 0, 1, 1]} by default, or the outer position of the axes it
    ## took the place of.
    ##
    ## @end deftp
    OuterPosition = [0, 0, 1, 1];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} PositionConstraint
    ##
    ## Which position the chart keeps: @qcode{'outerposition'}, the default,
    ## or @qcode{'innerposition'}.
    ##
    ## @end deftp
    PositionConstraint = 'outerposition';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} Units
    ##
    ## The units of the positions, @qcode{'normalized'} by default.
    ##
    ## @end deftp
    Units = 'normalized';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} Visible
    ##
    ## Whether the chart is shown, @qcode{'on'} by default or @qcode{'off'}.
    ##
    ## @end deftp
    Visible = 'on';

  endproperties

  properties (Dependent)

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} InnerPosition
    ##
    ## The same as @code{Position}.
    ##
    ## @end deftp
    InnerPosition;

  endproperties

  properties (GetAccess = public, SetAccess = private)

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ParallelCoordinatesPlot} {property} Parent
    ##
    ## The figure or panel the chart is placed in.  This property is
    ## read-only.
    ##
    ## @end deftp
    Parent = [];

  endproperties

  properties (Access = private, Hidden)
    Axes_ = [];
    Legend_ = [];
    Drawn_ = false;
    Redrawing_ = false;
    Internal_ = false;       # true while the chart sets a value itself
    Filling_ = false;        # true while the constructor sets the data
    StyleAuto_ = true;       # colours and styles follow the groups
    LegendAuto_ = true;
    LabelsAuto_ = true;      # CoordinateTickLabels follow the coordinates
    TitleAuto_ = true;       # LegendTitle follows GroupVariable
    Offsets_ = [];           # the jitter of each observation, drawn once
  endproperties

  methods (Hidden)

    function disp (this)
      printf ('  ParallelCoordinatesPlot with properties:\n\n');
      if (isempty (this.SourceTable))
        printf ('%14s: [%dx%d double]\n', 'Data', size (this.Data));
        printf ('%14s: %s\n', 'CoordinateData', mat2str (this.CoordinateData));
        if (isempty (this.GroupData))
          printf ('%14s: []\n', 'GroupData');
        else
          printf ('%14s: [%dx1 %s]\n', 'GroupData', numel (this.GroupData), ...
                  class (this.GroupData));
        endif
      else
        printf ('%19s: [%dx%d table]\n', 'SourceTable', ...
                size (this.SourceTable));
        cv = this.CoordinateVariables;
        if (iscellstr (cv))
          cv = ['{', strjoin(strcat ("'", cv, "'"), '  '), '}'];
        else
          cv = mat2str (cv);
        endif
        printf ('%19s: %s\n', 'CoordinateVariables', cv);
        printf ('%19s: ''%s''\n', 'GroupVariable', this.GroupVariable);
      endif
      printf ('\n');
    endfunction

    function display (this)
      disp (this);
    endfunction

  endmethods

  methods (Hidden)

    function set.Data (this, val)
      if (! ((isnumeric (val) || islogical (val)) && isreal (val)
             && ndims (val) == 2))
        error (strcat ("stats.chart.ParallelCoordinatesPlot: 'Data' must", ...
                       " be a real numeric matrix."));
      endif
      if (this.Drawn_ && ! isempty (this.SourceTable))
        error (strcat ("stats.chart.ParallelCoordinatesPlot: setting", ...
                       " 'Data' while 'SourceTable' holds a table is not", ...
                       " supported."));
      endif
      this.Data = double (val);
      if (! this.Filling_)
        this.CoordinateData = 1:columns (this.Data);
      endif
    endfunction

    function set.CoordinateData (this, val)
      if (! isempty (this.SourceTable))
        error (strcat ("stats.chart.ParallelCoordinatesPlot: setting", ...
                       " 'CoordinateData' while 'SourceTable' holds a", ...
                       " table is not supported."));
      endif
      this.CoordinateData = pcCheckIndex (val, columns (this.Data), ...
                                          'CoordinateData');
      refit (this);
    endfunction

    function set.SourceTable (this, val)
      if (! (isempty (val) || istable (val)))
        error (strcat ("stats.chart.ParallelCoordinatesPlot:", ...
                       " 'SourceTable' must be a table."));
      endif
      this.SourceTable = val;
      if (! this.Filling_ && ! isempty (val))
        this.CoordinateVariables = val.Properties.VariableNames;
      endif
    endfunction

    function set.CoordinateVariables (this, val)
      t = this.SourceTable;
      if (isempty (t))
        error (strcat ("stats.chart.ParallelCoordinatesPlot: setting", ...
                       " 'CoordinateVariables' after setting 'Data' is", ...
                       " not supported."));
      endif
      names = t.Properties.VariableNames;
      if (isa (val, 'string'))
        val = cellstr (val);
      elseif (ischar (val))
        val = {val};
      endif
      if (iscellstr (val))
        if (! all (ismember (val, names)))
          error (strcat ("stats.chart.ParallelCoordinatesPlot:", ...
                         " 'CoordinateVariables' does not name variables", ...
                         " of the source table."));
        endif
        val = val(:)';
      else
        pcCheckIndex (val, numel (names), 'CoordinateVariables');
      endif
      this.CoordinateVariables = val;
      refit (this);
    endfunction

    function set.GroupData (this, val)
      if (! isempty (val) && ! isvector (val))
        error (strcat ("stats.chart.ParallelCoordinatesPlot: 'GroupData'", ...
                       " must be a vector."));
      endif
      if (isempty (val))
        this.GroupData = [];
      else
        this.GroupData = val(:);
      endif
      refit (this);
    endfunction

    function set.GroupVariable (this, val)
      if (isa (val, 'string') && isscalar (val))
        val = char (val);
      endif
      if (! (ischar (val) && (isrow (val) || isempty (val))))
        error (strcat ("stats.chart.ParallelCoordinatesPlot:", ...
                       " 'GroupVariable' must be a character vector."));
      endif
      t = this.SourceTable;
      if (! isempty (val)
          && (isempty (t) || ! any (strcmp (t.Properties.VariableNames, val))))
        error (strcat ("stats.chart.ParallelCoordinatesPlot:", ...
                       " 'GroupVariable' does not name a variable of the", ...
                       " source table."));
      endif
      this.GroupVariable = val;
      if (this.TitleAuto_)
        this.Internal_ = true;
        this.LegendTitle = val;
        this.Internal_ = false;
      endif
      refit (this);
    endfunction

    function set.DataNormalization (this, val)
      this.DataNormalization = pcCheckOneOf (val, {'range', 'none', ...
                                                   'zscore', 'scale', ...
                                                   'center', 'norm'}, ...
                                             'DataNormalization');
      redraw (this);
    endfunction

    function set.CoordinateTickLabels (this, val)
      if (isa (val, 'string'))
        val = cellstr (val);
      endif
      if (! iscellstr (val))
        error (strcat ("stats.chart.ParallelCoordinatesPlot:", ...
                       " 'CoordinateTickLabels' must be a cell array of", ...
                       " character vectors."));
      endif
      this.CoordinateTickLabels = val(:);
      if (! this.Internal_)
        this.LabelsAuto_ = false;
        redraw (this);
      endif
    endfunction

    function set.Jitter (this, val)
      if (! (isnumeric (val) && isreal (val) && isscalar (val)
             && val >= 0 && val <= 1))
        error (strcat ("stats.chart.ParallelCoordinatesPlot: 'Jitter' must", ...
                       " be a number from 0 to 1."));
      endif
      this.Jitter = double (val);
      redraw (this);
    endfunction

    function set.Color (this, val)
      if (! (isnumeric (val) && isreal (val) && columns (val) == 3
             && rows (val) > 0 && all (val(:) >= 0) && all (val(:) <= 1)))
        error (strcat ("stats.chart.ParallelCoordinatesPlot: 'Color' must", ...
                       " be an RGB triplet or a matrix of them."));
      endif
      this.Color = double (val);
      if (! this.Internal_)
        this.StyleAuto_ = false;
        redraw (this);
      endif
    endfunction

    function set.LineStyle (this, val)
      this.LineStyle = pcCheckStyles (val, {'-', '--', ':', '-.', 'none'}, ...
                                      'LineStyle');
      if (! this.Internal_)
        redraw (this);
      endif
    endfunction

    function set.LineWidth (this, val)
      this.LineWidth = pcCheckPositive (val, 'LineWidth');
      if (! this.Internal_)
        redraw (this);
      endif
    endfunction

    function set.LineAlpha (this, val)
      if (! (isnumeric (val) && isreal (val) && isvector (val)
             && all (val >= 0) && all (val <= 1)))
        error (strcat ("stats.chart.ParallelCoordinatesPlot: 'LineAlpha'", ...
                       " must be from 0 to 1."));
      endif
      this.LineAlpha = double (val(:)');
      if (! this.Internal_)
        redraw (this);
      endif
    endfunction

    function set.MarkerStyle (this, val)
      this.MarkerStyle = pcCheckStyles (val, {'o', '+', '*', '.', 'x', 's', ...
                                              'd', '^', 'v', '>', '<', 'p', ...
                                              'h', 'none'}, 'MarkerStyle');
      if (! this.Internal_)
        redraw (this);
      endif
    endfunction

    function set.MarkerSize (this, val)
      this.MarkerSize = pcCheckPositive (val, 'MarkerSize');
      if (! this.Internal_)
        redraw (this);
      endif
    endfunction

    function set.LegendVisible (this, val)
      this.LegendVisible = pcCheckOneOf (val, {'on', 'off'}, 'LegendVisible');
      if (! this.Internal_)
        this.LegendAuto_ = false;
        redraw (this);
      endif
    endfunction

    function set.LegendTitle (this, val)
      this.LegendTitle = pcCheckText (val, 'LegendTitle');
      if (! this.Internal_)
        this.TitleAuto_ = false;
        redraw (this);
      endif
    endfunction

    function set.CoordinateLabel (this, val)
      this.CoordinateLabel = pcCheckText (val, 'CoordinateLabel');
      redraw (this);
    endfunction

    function set.DataLabel (this, val)
      this.DataLabel = pcCheckText (val, 'DataLabel');
      redraw (this);
    endfunction

    function set.Title (this, val)
      this.Title = pcCheckText (val, 'Title');
      redraw (this);
    endfunction

    function set.FontName (this, val)
      this.FontName = pcCheckText (val, 'FontName');
      redraw (this);
    endfunction

    function set.FontSize (this, val)
      if (! (isnumeric (val) && isreal (val) && isscalar (val)
             && isfinite (val) && val > 0))
        error (strcat ("stats.chart.ParallelCoordinatesPlot: 'FontSize'", ...
                       " must be a positive number."));
      endif
      this.FontSize = double (val);
      redraw (this);
    endfunction

    function set.Position (this, val)
      this.Position = pcCheckPosition (val, 'Position');
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
      this.OuterPosition = pcCheckPosition (val, 'OuterPosition');
      if (! this.Internal_)
        this.PositionConstraint = 'outerposition';
      endif
    endfunction

    function set.PositionConstraint (this, val)
      this.PositionConstraint = pcCheckOneOf (val, {'outerposition', ...
                                                    'innerposition'}, ...
                                              'PositionConstraint');
      redraw (this);
    endfunction

    function set.Units (this, val)
      this.Units = pcCheckOneOf (val, {'normalized', 'inches', ...
                                       'centimeters', 'points', 'pixels', ...
                                       'characters'}, 'Units');
      redraw (this);
    endfunction

    function set.Visible (this, val)
      this.Visible = pcCheckOneOf (val, {'on', 'off'}, 'Visible');
      redraw (this);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {stats.chart.ParallelCoordinatesPlot} {@var{obj} =} stats.chart.ParallelCoordinatesPlot (@var{parent}, @var{spec}, @var{args})
    ##
    ## Create a @code{stats.chart.ParallelCoordinatesPlot} object.
    ##
    ## @var{parent} is the figure or panel to place the chart in, or empty for
    ## the current figure, which is resolved only once every value has been
    ## accepted.  @var{spec} is a structure carrying the data as
    ## @code{parallelplot} resolved it, with the fields @qcode{Data},
    ## @qcode{SourceTable}, @qcode{CoordinateData},
    ## @qcode{CoordinateVariables}, @qcode{GroupData} and
    ## @qcode{GroupVariable}, and @var{args} the name-value pairs left to set.
    ## The documented way to reach this constructor is @code{parallelplot}.
    ##
    ## @seealso{parallelplot}
    ## @end deftypefn
    function this = ParallelCoordinatesPlot (parent, spec, args)

      if (nargin < 2)
        error (strcat ("stats.chart.ParallelCoordinatesPlot: too few input", ...
                       " arguments."));
      endif
      if (nargin < 3)
        args = {};
      endif

      this.Filling_ = true;
      if (isempty (spec.SourceTable))
        this.Data = spec.Data;
        this.Filling_ = false;
        if (isempty (spec.CoordinateData))
          this.CoordinateData = 1:columns (this.Data);
        else
          this.CoordinateData = spec.CoordinateData;
        endif
        this.GroupData = spec.GroupData;
      else
        this.SourceTable = spec.SourceTable;
        this.Filling_ = false;
        if (isempty (spec.CoordinateVariables))
          this.CoordinateVariables = spec.SourceTable.Properties.VariableNames;
        else
          this.CoordinateVariables = spec.CoordinateVariables;
        endif
        this.GroupVariable = spec.GroupVariable;
      endif

      ## Whatever was named is set now, so that one drawing covers them all
      placed = false;
      for k = 1:2:numel (args)
        name = args{k};
        this.(name) = args{k+1};
        if (any (strcmp (name, {'Position', 'InnerPosition', ...
                                'OuterPosition'})))
          placed = true;
        endif
      endfor
      [~, ~, n] = observations (this);
      if (! isempty (this.GroupData) && numel (this.GroupData) != n)
        error (strcat ("stats.chart.ParallelCoordinatesPlot: 'GroupData'", ...
                       " must hold one value for each observation."));
      endif
      if (! this.LabelsAuto_
          && numel (this.CoordinateTickLabels) != numel (coordinates (this)))
        error (strcat ("stats.chart.ParallelCoordinatesPlot:", ...
                       " 'CoordinateTickLabels' must hold one label for", ...
                       " each coordinate."));
      endif
      refit (this);

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
          error (strcat ("stats.chart.ParallelCoordinatesPlot: a parallel", ...
                         " coordinates plot cannot be added to axes on", ...
                         " which hold is on."));
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
      this.Axes_ = axes ('parent', parent, ...
                         'tag', 'stats.chart.ParallelCoordinatesPlot');
      set (fig, 'currentaxes', this.Axes_);
      ## The chart goes with its axes, and its axes with the chart
      set (this.Axes_, 'deletefcn', @(~, ~) delete (this));
      this.Drawn_ = true;
      redraw (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {stats.chart.ParallelCoordinatesPlot} {} delete (@var{h})
    ##
    ## Delete the chart and the axes it is drawn in.
    ##
    ## @end deftypefn
    function delete (this)
      ax = this.Axes_;
      lg = this.Legend_;
      this.Axes_ = [];
      this.Legend_ = [];
      this.Parent = [];
      if (! isempty (lg) && ishghandle (lg))
        delete (lg);
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
      if (! isempty (this.Legend_) && ishghandle (this.Legend_))
        delete (this.Legend_);
      endif
      this.Legend_ = [];
    endfunction

  endmethods

  methods (Access = private)

    ## The columns or variables drawn, as indices into the matrix or table.
    function idx = coordinates (this)
      if (isempty (this.SourceTable))
        idx = this.CoordinateData;
        n = columns (this.Data);
      else
        idx = this.CoordinateVariables;
        names = this.SourceTable.Properties.VariableNames;
        n = numel (names);
        if (iscellstr (idx))
          [~, idx] = ismember (idx, names);
        endif
      endif
      if (islogical (idx))
        idx = find (idx);
      endif
      idx = idx(idx >= 1 & idx <= n);
      idx = idx(:)';
    endfunction

    ## The values of each coordinate: a matrix of positions, one column for
    ## each, the levels of the columns holding categories, and the number of
    ## observations.
    function [V, levels, n] = observations (this)
      idx = coordinates (this);
      levels = cell (1, numel (idx));
      if (isempty (this.SourceTable))
        V = this.Data(:,idx);
        n = rows (this.Data);
        return;
      endif
      t = this.SourceTable;
      n = height (t);
      V = NaN (n, numel (idx));
      for j = 1:numel (idx)
        v = t{:,idx(j)};
        if ((isnumeric (v) || islogical (v)) && ! islogical (v))
          V(:,j) = double (v(:,1));
        else
          [code, levels{j}] = pcLevels (v);
          V(:,j) = code;
        endif
      endfor
    endfunction

    ## Bring the labels, styles and legend that follow the data in line with
    ## it, then draw.
    function refit (this)

      if (this.Filling_)
        return;
      endif
      [names, ~] = pcGroups (currentGroups (this));
      ng = max (numel (names), 1);
      this.Internal_ = true;
      unwind_protect
        if (this.LabelsAuto_)
          if (isempty (this.SourceTable))
            this.CoordinateTickLabels = arrayfun (@(k) sprintf ('%d', k), ...
                                                  coordinates (this)', ...
                                                  'UniformOutput', false);
          else
            names_t = this.SourceTable.Properties.VariableNames;
            this.CoordinateTickLabels = names_t(coordinates (this))';
          endif
        endif
        if (this.StyleAuto_)
          co = get (groot, 'defaultaxescolororder');
          this.Color = co(mod ((1:ng) - 1, rows (co)) + 1,:);
          if (ng > 1)
            this.LineStyle = repmat ({'-'}, 1, ng);
            this.MarkerStyle = repmat ({'none'}, 1, ng);
            this.LineWidth = ones (1, ng);
            this.LineAlpha = 0.7 * ones (1, ng);
            this.MarkerSize = 6 * ones (1, ng);
          endif
        endif
        if (this.LegendAuto_)
          this.LegendVisible = pcOnOff (! isempty (names));
        endif
      unwind_protect_cleanup
        this.Internal_ = false;
      end_unwind_protect
      this.Offsets_ = [];
      redraw (this);

    endfunction

    ## The groups of the observations, from the table or GroupData.
    function g = currentGroups (this)
      if (! isempty (this.SourceTable) && ! isempty (this.GroupVariable))
        g = this.SourceTable.(this.GroupVariable);
      else
        g = this.GroupData;
      endif
    endfunction

    ## Draw the chart from scratch into its axes.
    function redraw (this)

      if (! this.Drawn_ || isempty (this.Axes_) || ! ishghandle (this.Axes_))
        return;
      endif
      ax = this.Axes_;
      this.Redrawing_ = true;
      unwind_protect
        if (! isempty (this.Legend_) && ishghandle (this.Legend_))
          delete (this.Legend_);
        endif
        this.Legend_ = [];
        delete (get (ax, 'children'));
      unwind_protect_cleanup
        this.Redrawing_ = false;
      end_unwind_protect
      ## A child whose deletion by anyone else releases the axes
      text (ax, NaN, NaN, '', 'deletefcn', @(~, ~) release (this));

      [V, levels, n] = observations (this);
      nc = columns (V);
      iscat = ! cellfun (@isempty, levels);
      [names, g] = pcGroups (currentGroups (this));
      ng = max (numel (names), 1);
      if (isempty (names))
        g = ones (n, 1);
      endif

      ## Each value as a position along its ruler: on rulers of their own,
      ## from 0 at the smallest to 1 at the largest, or on one shared scale
      how = this.DataNormalization;
      P = NaN (n, nc);
      ticks = cell (1, nc);
      for j = 1:nc
        v = V(:,j);
        if (iscat(j))
          P(:,j) = v;
        else
          P(:,j) = pcNormalize (v, how);
        endif
      endfor
      if (strcmp (how, 'range'))
        for j = 1:nc
          [P(:,j), ticks{j}] = pcRuler (P(:,j), levels{j});
        endfor
        yl = [0, 1];
      else
        num = P(:,! iscat);
        lo = min (num(:));
        hi = max (num(:));
        if (isempty (lo) || ! isfinite (lo))
          lo = 0;
          hi = 1;
        elseif (lo == hi)
          lo -= 0.5;
          hi += 0.5;
        endif
        yl = [lo, hi];
        for j = find (iscat)
          k = numel (levels{j});
          P(:,j) = lo + (P(:,j) - 1) / max (k - 1, 1) * (hi - lo);
          tv = lo + ((1:k) - 1) / max (k - 1, 1) * (hi - lo);
          ticks{j} = {tv, levels{j}};
        endfor
      endif

      ## Lines through the values of a variable of categories, spread over a
      ## share of the room up to the next value, drawn once for each set of
      ## observations
      if (rows (this.Offsets_) != n || columns (this.Offsets_) != nc)
        this.Offsets_ = rand (n, nc) - 0.5;
      endif
      for j = find (iscat)
        k = numel (levels{j});
        step = (yl(2) - yl(1)) / max (k - 1, 1);
        P(:,j) += this.Jitter * step * this.Offsets_(:,j);
      endfor

      ## The lines of each group, observations joined across the rulers and
      ## parted by a missing value
      bg = get (ax, 'color');
      if (ischar (bg))
        bg = [1, 1, 1];
      endif
      hl = zeros (1, ng);
      xs = [1:nc, NaN];
      for k = 1:ng
        take = find (g == k);
        px = repmat (xs, numel (take), 1)';
        py = [P(take,:), NaN(numel (take), 1)]';
        c = pcPick (this.Color, k);
        a = pcPick (this.LineAlpha(:), k);
        hl(k) = line (ax, px(:), py(:), 'color', a * c + (1 - a) * bg, ...
                      'linestyle', pcPickStyle (this.LineStyle, k), ...
                      'linewidth', pcPick (this.LineWidth(:), k), ...
                      'marker', pcPickStyle (this.MarkerStyle, k), ...
                      'markersize', pcPick (this.MarkerSize(:), k), ...
                      'visible', this.Visible);
        if (! isempty (names))
          set (hl(k), 'displayname', names{k});
        endif
      endfor

      ## The rulers, one for each coordinate, with their own ticks where they
      ## have their own scale or hold categories
      for j = 1:nc
        line (ax, [j, j], yl, 'color', [0.15, 0.15, 0.15], ...
              'visible', this.Visible);
        if (! isempty (ticks{j}))
          tv = ticks{j}{1};
          tl = ticks{j}{2};
          for i = 1:numel (tv)
            text (ax, j, tv(i), tl{i}, 'horizontalalignment', 'right', ...
                  'fontsize', 0.9 * this.FontSize, 'fontname', ...
                  this.FontName, 'backgroundcolor', bg, 'margin', 1, ...
                  'visible', this.Visible);
          endfor
        endif
      endfor
      set (ax, 'xlim', [0.5, nc + 0.5], 'ylim', yl, 'xtick', 1:nc, ...
           'xticklabel', this.CoordinateTickLabels, 'box', 'off', ...
           'xgrid', 'off', 'units', this.Units, 'visible', this.Visible, ...
           'fontname', this.FontName, 'fontsize', this.FontSize);
      if (strcmp (how, 'range'))
        set (ax, 'ytick', [], 'ycolor', 'none');
      else
        set (ax, 'ytickmode', 'auto', 'ycolor', [0.15, 0.15, 0.15]);
      endif
      title (ax, this.Title, 'fontname', this.FontName);
      xlabel (ax, this.CoordinateLabel, 'fontname', this.FontName);
      ylabel (ax, this.DataLabel, 'fontname', this.FontName);

      ## Where the axes goes: within the outer position, leaving the labels
      ## the room they take and the legend its own
      showLegend = strcmp (this.LegendVisible, 'on') && ! isempty (names);
      if (strcmp (this.PositionConstraint, 'outerposition'))
        o = this.OuterPosition;
        ti = get (ax, 'tightinset');
        pad = 0.02;
        L = max (0.13 * o(3), ti(1) + pad * o(3));
        B = max (0.11 * o(4), ti(2) + pad * o(4));
        T = max (0.075 * o(4), ti(4) + pad * o(4));
        R = 0.095 * o(3);
        if (showLegend)
          R = 0.2 * o(3);
        endif
        w = max (o(3) - L - R, 0.01 * o(3));
        h = max (o(4) - B - T, 0.01 * o(4));
        this.Internal_ = true;
        unwind_protect
          this.Position = [o(1) + L, o(2) + B, w, h];
        unwind_protect_cleanup
          this.Internal_ = false;
        end_unwind_protect
      endif
      set (ax, 'position', this.Position);

      ## The legend of the groups, beside the rulers
      if (showLegend)
        this.Legend_ = legend (ax, hl, names(:)', 'location', ...
                               'northeastoutside');
        set (this.Legend_, 'fontname', this.FontName, 'visible', this.Visible);
        if (! isempty (this.LegendTitle))
          title (this.Legend_, this.LegendTitle);
        endif
        set (ax, 'position', this.Position);
      endif

    endfunction

  endmethods

endclassdef

## Values of a variable normalized as asked
function v = pcNormalize (v, how)
  ok = isfinite (v);
  switch (how)
    case 'zscore'
      s = std (v(ok));
      v = (v - mean (v(ok))) / (s + (s == 0));
    case 'scale'
      s = std (v(ok));
      v = v / (s + (s == 0));
    case 'center'
      v = v - mean (v(ok));
    case 'norm'
      s = norm (v(ok));
      v = v / (s + (s == 0));
  endswitch
endfunction

## Positions from 0 to 1 along a ruler of its own, with the ticks it carries:
## the levels of a variable of categories, evenly spaced, or its smallest,
## largest and a few values between
function [p, ticks] = pcRuler (v, levels)
  if (! isempty (levels))
    k = numel (levels);
    p = (v - 1) / max (k - 1, 1);
    tv = ((1:k) - 1) / max (k - 1, 1);
    ticks = {tv, levels};
    return;
  endif
  ok = isfinite (v);
  if (! any (ok))
    p = v;
    ticks = {[], {}};
    return;
  endif
  lo = min (v(ok));
  hi = max (v(ok));
  if (lo == hi)
    p = 0.5 * ones (size (v));
    p(! ok) = NaN;
    ticks = {0.5, {sprintf('%.4g', lo)}};
    return;
  endif
  p = (v - lo) / (hi - lo);
  tv = pcNiceTicks (lo, hi);
  tl = arrayfun (@(t) sprintf ('%.4g', t), tv, 'UniformOutput', false);
  ticks = {(tv - lo) / (hi - lo), tl};
endfunction

## Ticks between LO and HI at a round step of 1, 2 or 5 times a power of
## ten, three to seven of them
function tv = pcNiceTicks (lo, hi)
  raw = (hi - lo) / 5;
  e = 10 ^ floor (log10 (raw));
  m = [1, 2, 5, 10];
  step = e * m(find (m * e >= raw, 1));
  tv = (ceil (lo / step - 1e-9):floor (hi / step + 1e-9)) * step;
  tv(abs (tv) < step * 1e-9) = 0;
endfunction

## The codes and names of the levels of a variable of categories, text or
## logical values: categories in their order, the rest sorted; NaN where
## missing
function [code, names] = pcLevels (v)
  v = v(:,1);
  if (isa (v, 'categorical'))
    names = categories (v)';
    [~, code] = ismember (cellstr (v), names);
    code = double (code);
    code(isundefined (v)) = NaN;
  elseif (islogical (v))
    names = {'false', 'true'};
    code = double (v) + 1;
  else
    t = cellstr (v);
    miss = cellfun (@isempty, t);
    if (isa (v, 'string'))
      miss |= ismissing (v);
    endif
    names = unique (t(! miss))';
    code = NaN (numel (t), 1);
    [~, code(! miss)] = ismember (t(! miss), names);
  endif
endfunction

## The names of the groups in order of first appearance, and the group of
## each observation, 0 where it is missing
function [names, g] = pcGroups (gd)
  names = cell (0, 1);
  g = zeros (numel (gd), 1);
  if (isempty (gd))
    return;
  endif
  if (isa (gd, 'categorical'))
    miss = isundefined (gd);
    t = cellstr (gd);
  elseif (iscellstr (gd) || isa (gd, 'string'))
    t = cellstr (gd);
    miss = cellfun (@isempty, t);
    if (isa (gd, 'string'))
      miss |= ismissing (gd);
    endif
  elseif (isnumeric (gd) || islogical (gd))
    miss = isnan (double (gd));
    t = arrayfun (@(v) sprintf ('%g', v), double (gd), 'UniformOutput', false);
    if (islogical (gd))
      t(gd) = {'true'};
      t(! gd) = {'false'};
    endif
  else
    error (strcat ("stats.chart.ParallelCoordinatesPlot: 'GroupData'", ...
                   " must hold categories, text, numbers or logical values."));
  endif
  t = t(:);
  miss = miss(:);
  [names, i] = unique (t(! miss), 'first');
  [~, o] = sort (i);
  names = names(o);
  [~, g(! miss)] = ismember (t(! miss), names);
endfunction

## Indices into N columns: positive integers not above N, or a logical
## vector of length N
function v = pcCheckIndex (val, n, name)
  if (islogical (val) && isvector (val) && numel (val) == n)
    v = val(:)';
  elseif (isnumeric (val) && isreal (val) && (isvector (val) || isempty (val))
          && all (val == fix (val)) && all (val >= 1) && all (val <= n))
    v = double (val(:)');
  else
    error (strcat ("stats.chart.ParallelCoordinatesPlot: '%s' must be", ...
                   " indices of columns or a logical vector over them."), name);
  endif
endfunction

## Row K of a matrix of values for the groups, the last one repeated
function v = pcPick (vals, k)
  v = vals(min (k, rows (vals)),:);
endfunction

## Style K of a character vector or a cell array of styles
function s = pcPickStyle (styles, k)
  if (iscell (styles))
    s = styles{min (k, numel (styles))};
  else
    s = styles;
  endif
endfunction

## 'on' or 'off'
function s = pcOnOff (tf)
  if (tf)
    s = 'on';
  else
    s = 'off';
  endif
endfunction

## Styles: one of LIST, or a cell array or string array of them
function v = pcCheckStyles (val, list, name)
  if (isa (val, 'string'))
    val = cellstr (val);
    if (numel (val) == 1)
      val = val{1};
    endif
  endif
  if (ischar (val) && any (strcmp (val, list)))
    v = val;
  elseif (iscellstr (val) && ! isempty (val) && all (ismember (val, list)))
    v = val(:)';
  else
    error ("stats.chart.ParallelCoordinatesPlot: '%s' must be one of %s.", ...
           name, strjoin (strcat ("'", list, "'"), ', '));
  endif
endfunction

## Positive values, a scalar or one for each group
function v = pcCheckPositive (val, name)
  if (! (isnumeric (val) && isreal (val) && isvector (val)
         && all (isfinite (val)) && all (val > 0)))
    error (strcat ("stats.chart.ParallelCoordinatesPlot: '%s' must be", ...
                   " positive."), name);
  endif
  v = double (val(:)');
endfunction

## Text: a character vector, a string, or a cell array of them for lines
function out = pcCheckText (val, name)
  if (isa (val, 'string'))
    val = cellstr (val);
    if (numel (val) == 1)
      val = val{1};
    endif
  endif
  if (! ((ischar (val) && (isrow (val) || isempty (val))) || iscellstr (val)))
    error ("stats.chart.ParallelCoordinatesPlot: '%s' must be text.", name);
  endif
  out = val;
endfunction

## One of a list of names, case free
function out = pcCheckOneOf (val, list, name)
  if (isa (val, 'string') && isscalar (val))
    val = char (val);
  endif
  if (! (ischar (val) && any (strcmpi (val, list))))
    error ("stats.chart.ParallelCoordinatesPlot: '%s' must be one of %s.", ...
           name, strjoin (strcat ("'", list, "'"), ', '));
  endif
  out = lower (val);
endfunction

## A position, as a 1-by-4 vector with positive width and height
function out = pcCheckPosition (val, name)
  if (! (isnumeric (val) && isreal (val) && numel (val) == 4
         && all (isfinite (val)) && val(3) > 0 && val(4) > 0))
    error (strcat ("stats.chart.ParallelCoordinatesPlot: '%s' must be a", ...
                   " 1-by-4 vector with positive width and height."), name);
  endif
  out = double (val(:)');
endfunction

%!shared X, g, T
%! X = [2.1, 1.2, 10; 3.4, 2.8, 12; 1.9, 0.9, 11; 5.6, 2.2, 15; ...
%!      4.4, 3.1, 9; 3.8, 1.7, 13];
%! g = categorical ({'a'; 'b'; 'a'; 'b'; 'c'; 'a'});
%! T = table (X(:,1), X(:,2), X(:,3), g, 'VariableNames', {'A', 'B', 'C', 'G'});
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot (X);
%!   assert_equal (class (h), 'stats.chart.ParallelCoordinatesPlot');
%!   assert_equal (h.Data, X);
%!   assert_equal (h.CoordinateData, [1, 2, 3]);
%!   assert_equal (h.GroupData, []);
%!   assert_equal (h.CoordinateTickLabels, {'1'; '2'; '3'});
%!   assert_equal (h.DataNormalization, 'range');
%!   assert_equal (h.Color, [0, 0.447, 0.741]);
%!   assert_equal (h.LineWidth, 1);
%!   assert_equal (h.LineStyle, '-');
%!   assert_equal (h.LineAlpha, 0.7);
%!   assert_equal (h.MarkerStyle, 'none');
%!   assert_equal (h.MarkerSize, 6);
%!   assert_equal (h.Jitter, 0.1);
%!   assert_equal (h.LegendVisible, 'off');
%!   assert_equal (h.FontSize, 10);
%!   assert_equal (h.OuterPosition, [0, 0, 1, 1]);
%!   assert_equal (h.Parent, hf);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot (X, 'GroupData', g);
%!   assert_equal (h.Color, [0, 0.447, 0.741; 0.85, 0.325, 0.098; ...
%!                           0.929, 0.694, 0.125]);
%!   assert_equal (h.LineStyle, {'-', '-', '-'});
%!   assert_equal (h.LineWidth, [1, 1, 1]);
%!   assert_equal (h.LineAlpha, [0.7, 0.7, 0.7]);
%!   assert_equal (h.MarkerStyle, {'none', 'none', 'none'});
%!   assert_equal (h.MarkerSize, [6, 6, 6]);
%!   assert_equal (h.LegendVisible, 'on');
%!   assert_equal (h.LegendTitle, '');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot (X, 'GroupData', {'p'; 'q'; 'p'; 'q'; 'p'; 'q'});
%!   assert_equal (rows (h.Color), 2);
%!   assert_equal (h.LegendVisible, 'on');
%!   lg = findall (hf, 'type', 'axes', 'tag', 'legend');
%!   s = get (lg, 'string');   # a row under gnuplot, a column otherwise
%!   assert_equal (s(:), {'p'; 'q'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot (T);
%!   assert_equal (h.CoordinateVariables, {'A', 'B', 'C', 'G'});
%!   assert_equal (h.CoordinateTickLabels, {'A'; 'B'; 'C'; 'G'});
%!   assert_equal (h.GroupVariable, '');
%!   assert_equal (h.LegendTitle, '');
%!   assert_equal (h.LegendVisible, 'off');
%!   assert_equal (h.Data, []);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot (T, 'GroupVariable', 'G');
%!   assert_equal (h.LegendTitle, 'G');
%!   assert_equal (h.LegendVisible, 'on');
%!   assert_equal (rows (h.Color), 3);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot (T, 'CoordinateVariables', {'C', 'A'});
%!   assert_equal (h.CoordinateVariables, {'C', 'A'});
%!   assert_equal (h.CoordinateTickLabels, {'C'; 'A'});
%!   h.CoordinateVariables = [3, 1];
%!   assert_equal (h.CoordinateVariables, [3, 1]);
%!   h.CoordinateVariables = [true, false, true, false];
%!   assert_equal (h.CoordinateVariables, [true, false, true, false]);
%!   assert_equal (h.CoordinateTickLabels, {'A'; 'C'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot (X, 'CoordinateTickLabels', {'x', 'y', 'z'}, ...
%!                     'Jitter', 0.3, 'DataNormalization', 'zscore');
%!   assert_equal (h.CoordinateTickLabels, {'x'; 'y'; 'z'});
%!   assert_equal (h.Jitter, 0.3);
%!   assert_equal (h.DataNormalization, 'zscore');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot (T, 'GroupVariable', 'G', 'Color', ...
%!                     [1, 0, 0; 0, 1, 0; 0, 0, 1], 'LineStyle', ...
%!                     {'-', '--', ':'}, 'LineWidth', [1, 2, 3], ...
%!                     'MarkerStyle', {'o', 'x', 's'}, 'MarkerSize', ...
%!                     [4, 6, 8], 'LineAlpha', 0.3);
%!   assert_equal (h.Color, [1, 0, 0; 0, 1, 0; 0, 0, 1]);
%!   assert_equal (h.LineStyle, {'-', '--', ':'});
%!   assert_equal (h.LineWidth, [1, 2, 3]);
%!   assert_equal (h.MarkerStyle, {'o', 'x', 's'});
%!   assert_equal (h.MarkerSize, [4, 6, 8]);
%!   assert_equal (h.LineAlpha, 0.3);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot ([X; NaN, 1, 1]);
%!   assert_equal (size (h.Data), [7, 3]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! ## Values on rulers of their own run from 0 at the smallest to 1 at the
%! ## largest, and a missing value breaks its line
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot ([1, 10; 3, 20; 2, NaN]);
%!   tag = 'stats.chart.ParallelCoordinatesPlot';
%!   ax = findall (hf, 'type', 'axes', 'tag', tag);
%!   ln = findobj (ax, 'type', 'line', '-and', '-not', 'linestyle', 'none');
%!   yd = get (ln(end), 'ydata');
%!   assert_equal (yd(:)', [0, 0, NaN, 1, 1, NaN, 0.5, NaN, NaN], -1e-15);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = parallelplot ([1, 10; 3, 20; 2, 30], 'DataNormalization', 'zscore');
%!   tag = 'stats.chart.ParallelCoordinatesPlot';
%!   ax = findall (hf, 'type', 'axes', 'tag', tag);
%!   ln = findobj (ax, 'type', 'line');
%!   yd = get (ln(end), 'ydata');
%!   assert_equal (yd(:)', [-1, -1, NaN, 1, 0, NaN, 0, 1, NaN], -1e-15);
%!   h.DataNormalization = 'center';
%!   ln = findobj (ax, 'type', 'line');
%!   yd = get (ln(end), 'ydata');
%!   assert_equal (yd(:)', [-1, -10, NaN, 1, 0, NaN, 0, 10, NaN], -1e-14);
%!   h.DataNormalization = 'none';
%!   assert_equal (get (ax, 'ylim'), [1, 30]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! ## A variable of categories places its levels evenly along its ruler
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   t = table ([1; 2; 3], categorical ({'lo'; 'hi'; 'lo'}, {'lo', 'mid', ...
%!              'hi'}), 'VariableNames', {'N', 'C'});
%!   h = parallelplot (t, 'Jitter', 0);
%!   tag = 'stats.chart.ParallelCoordinatesPlot';
%!   ax = findall (hf, 'type', 'axes', 'tag', tag);
%!   ln = findobj (ax, 'type', 'line');
%!   yd = get (ln(end), 'ydata');
%!   assert_equal (yd(:)', [0, 0, NaN, 0.5, 1, NaN, 1, 0, NaN], -1e-15);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!shared pcS
%! pcS = struct ('Data', magic (3), 'SourceTable', [], ...
%!               'CoordinateData', [], 'CoordinateVariables', {{}}, ...
%!               'GroupData', [], 'GroupVariable', '');
%!error<stats.chart.ParallelCoordinatesPlot: too few input arguments.> ...
%! stats.chart.ParallelCoordinatesPlot ([])
%!error<stats.chart.ParallelCoordinatesPlot: 'Data' must be a real numeric matrix.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'Data', {1}})
%!error<stats.chart.ParallelCoordinatesPlot: 'CoordinateData' must be indices of columns or a logical vector over them.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'CoordinateData', 4})
%!error<stats.chart.ParallelCoordinatesPlot: setting 'CoordinateVariables' after setting 'Data' is not supported.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'CoordinateVariables', 1})
%!error<stats.chart.ParallelCoordinatesPlot: 'GroupData' must hold one value for each observation.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'GroupData', [1; 2]})
%!error<stats.chart.ParallelCoordinatesPlot: 'GroupData' must be a vector.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'GroupData', ones(3, 3)})
%!error<stats.chart.ParallelCoordinatesPlot: 'DataNormalization' must be one of 'range', 'none', 'zscore', 'scale', 'center', 'norm'.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'DataNormalization', 'norm1'})
%!error<stats.chart.ParallelCoordinatesPlot: 'CoordinateTickLabels' must be a cell array of character vectors.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'CoordinateTickLabels', 5})
%!error<stats.chart.ParallelCoordinatesPlot: 'CoordinateTickLabels' must hold one label for each coordinate.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, ...
%!                                     {'CoordinateTickLabels', {'x'}})
%!error<stats.chart.ParallelCoordinatesPlot: 'Jitter' must be a number from 0 to 1.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'Jitter', 2})
%!error<stats.chart.ParallelCoordinatesPlot: 'Color' must be an RGB triplet or a matrix of them.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'Color', [2, 0, 0]})
%!error<stats.chart.ParallelCoordinatesPlot: 'LineStyle' must be one of '-', '--', ':', '-.', 'none'.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'LineStyle', 'o'})
%!error<stats.chart.ParallelCoordinatesPlot: 'LineWidth' must be positive.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'LineWidth', 0})
%!error<stats.chart.ParallelCoordinatesPlot: 'LineAlpha' must be from 0 to 1.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'LineAlpha', 2})
## test.m ends a pattern at the first '>' of its line, and the list of
## markers holds one, so the message is matched only as far as 'v'
%!error<stats.chart.ParallelCoordinatesPlot: 'MarkerStyle' must be one of 'o', '\+', '\*', '.', 'x', 's', 'd', '\^', 'v', > ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'MarkerStyle', 'q'})
%!error<stats.chart.ParallelCoordinatesPlot: 'MarkerSize' must be positive.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'MarkerSize', -1})
%!error<stats.chart.ParallelCoordinatesPlot: 'LegendVisible' must be one of 'on', 'off'.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'LegendVisible', 1})
%!error<stats.chart.ParallelCoordinatesPlot: 'Title' must be text.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'Title', 5})
%!error<stats.chart.ParallelCoordinatesPlot: 'FontSize' must be a positive number.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'FontSize', 0})
%!error<stats.chart.ParallelCoordinatesPlot: 'Position' must be a 1-by-4 vector with positive width and height.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'Position', [0, 0, 1, 0]})
%!error<stats.chart.ParallelCoordinatesPlot: 'Visible' must be one of 'on', 'off'.> ...
%! stats.chart.ParallelCoordinatesPlot ([], pcS, {'Visible', 'maybe'})
%!error<stats.chart.ParallelCoordinatesPlot: 'CoordinateVariables' does not name variables of the source table.> ...
%! stats.chart.ParallelCoordinatesPlot ([], ...
%!   struct ('Data', [], 'SourceTable', ...
%!           table ([1; 2], 'VariableNames', {'A'}), ...
%!           'CoordinateData', [], 'CoordinateVariables', {{'Q'}}, ...
%!           'GroupData', [], 'GroupVariable', ''))
%!error<stats.chart.ParallelCoordinatesPlot: 'GroupVariable' does not name a variable of the source table.> ...
%! stats.chart.ParallelCoordinatesPlot ([], ...
%!   struct ('Data', [], 'SourceTable', ...
%!           table ([1; 2], 'VariableNames', {'A'}), ...
%!           'CoordinateData', [], 'CoordinateVariables', {{}}, ...
%!           'GroupData', [], 'GroupVariable', 'Q'))
