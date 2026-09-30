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
## @deftp {statistics} stats.chart.ScatterHistogramChart
##
## A scatter plot with a histogram of each variable along its sides, as
## @code{scatterhistogram} draws it.
##
## A @code{stats.chart.ScatterHistogramChart} holds two samples, optionally
## split into groups, and the choices they are drawn with, and redraws itself
## whenever one of them is set.  It is what @code{scatterhistogram} returns,
## and the documented way to reach it.
##
## The samples come either from vectors or from a table, never from both.
## The histograms are normalized as probability densities, one for each group,
## with the bins @code{histcounts} chooses for each group unless
## @code{NumBins} or @code{BinWidths} say otherwise.  The groups take the
## colours of the colour order and the line styles @qcode{'-'}, @qcode{':'},
## @qcode{'-.'} and @qcode{'--'} in turn, in the order in which each first
## appears in @code{GroupData}.
##
## A chart takes the place of the current axes of its parent, keeping their
## outer position, and refuses axes on which @code{hold} is on; a later plot
## into the figure takes the place of the chart in turn.  MATLAB's chart is
## a graphics object of its own; this one is a handle class drawing into
## three axes of its own, the scatter plot's being made current, and it lets
## go of them when a later plot clears them.  @code{LegendVisible},
## @code{MarkerFilled} and @code{Visible} hold @qcode{'on'} or @qcode{'off'},
## and @code{LineStyle} and @code{MarkerStyle} a character vector or a cell
## array of them, where MATLAB holds an on-off state and a string array.
##
## @seealso{scatterhistogram, scatterhist, histcounts, ksdensity}
## @end deftp

classdef ScatterHistogramChart < handle

  properties (Access = public)

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} XData
    ##
    ## The values along the horizontal axis, a numeric column vector.  Where
    ## the chart was drawn from a table it is read from the table and setting
    ## it is refused.
    ##
    ## @end deftp
    XData = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} YData
    ##
    ## The values along the vertical axis, a numeric column vector as long as
    ## @code{XData}.
    ##
    ## @end deftp
    YData = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} GroupData
    ##
    ## The group of each observation: a vector as long as @code{XData} of
    ## categories, text, numbers or logical values, or empty for no groups.
    ##
    ## @end deftp
    GroupData = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} SourceTable
    ##
    ## The table the chart was drawn from, or empty.
    ##
    ## @end deftp
    SourceTable = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} XVariable
    ##
    ## The variable of @code{SourceTable} holding @code{XData}.
    ##
    ## @end deftp
    XVariable = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} YVariable
    ##
    ## The variable of @code{SourceTable} holding @code{YData}.
    ##
    ## @end deftp
    YVariable = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} GroupVariable
    ##
    ## The variable of @code{SourceTable} holding @code{GroupData}, or empty.
    ##
    ## @end deftp
    GroupVariable = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} HistogramDisplayStyle
    ##
    ## How the histograms are drawn: @qcode{'stairs'}, the default, as their
    ## outline; @qcode{'bar'} as bars; @qcode{'smooth'} as a kernel density
    ## estimate, as @code{ksdensity} computes it.
    ##
    ## @end deftp
    HistogramDisplayStyle = 'stairs';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} NumBins
    ##
    ## The number of bins of each histogram, a 2-by-1 vector for the
    ## horizontal and the vertical variable, or a 2-by-N matrix with a column
    ## for each of N groups.  Setting a scalar sets both.  Empty where
    ## @code{HistogramDisplayStyle} is @qcode{'smooth'}.
    ##
    ## @end deftp
    NumBins = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} BinWidths
    ##
    ## The width of the bins, shaped as @code{NumBins}.  Setting it takes
    ## precedence over @code{NumBins}, which then follows from it.
    ##
    ## @end deftp
    BinWidths = [];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} ScatterPlotLocation
    ##
    ## Which corner the scatter plot takes, the histograms lining the two
    ## sides away from it: @qcode{'SouthWest'}, the default,
    ## @qcode{'SouthEast'}, @qcode{'NorthEast'} or @qcode{'NorthWest'}.
    ##
    ## @end deftp
    ScatterPlotLocation = 'SouthWest';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} ScatterPlotProportion
    ##
    ## The share of the width and of the height the scatter plot takes, a
    ## number between 0 and 1, 0.75 by default.
    ##
    ## @end deftp
    ScatterPlotProportion = 0.75;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} XHistogramDirection
    ##
    ## Which way the bars of the histogram of @code{XData} rise:
    ## @qcode{'up'}, the default, or @qcode{'down'}.
    ##
    ## @end deftp
    XHistogramDirection = 'up';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} YHistogramDirection
    ##
    ## Which way the bars of the histogram of @code{YData} rise:
    ## @qcode{'right'}, the default, or @qcode{'left'}.
    ##
    ## @end deftp
    YHistogramDirection = 'right';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} XLimits
    ##
    ## The limits of the horizontal axis, the range of @code{XData} until
    ## set.
    ##
    ## @end deftp
    XLimits = [0, 1];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} YLimits
    ##
    ## The limits of the vertical axis, the range of @code{YData} until set.
    ##
    ## @end deftp
    YLimits = [0, 1];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} Color
    ##
    ## The colour of each group, one row of an RGB matrix for each; by
    ## default the colours of the colour order.
    ##
    ## @end deftp
    Color = [0, 0.447, 0.741];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} LineStyle
    ##
    ## The line style of the histograms of each group, a character vector or
    ## a cell array of them.
    ##
    ## @end deftp
    LineStyle = '-';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} LineWidth
    ##
    ## The width of the histogram lines, a scalar or one for each group, 0.5
    ## by default.
    ##
    ## @end deftp
    LineWidth = 0.5;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} MarkerStyle
    ##
    ## The marker of each group, a character vector or a cell array of them,
    ## @qcode{'o'} by default.
    ##
    ## @end deftp
    MarkerStyle = 'o';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} MarkerSize
    ##
    ## The area of the markers in square points, a scalar or one for each
    ## group, 36 by default.
    ##
    ## @end deftp
    MarkerSize = 36;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} MarkerFilled
    ##
    ## Whether the markers are filled, @qcode{'on'} by default or
    ## @qcode{'off'}.
    ##
    ## @end deftp
    MarkerFilled = 'on';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} MarkerAlpha
    ##
    ## How opaque the markers are, from 0 to 1, 1 by default.
    ##
    ## @end deftp
    MarkerAlpha = 1;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} LegendVisible
    ##
    ## Whether the legend of the groups is shown, @qcode{'on'} or
    ## @qcode{'off'}; on where there are groups, until set.
    ##
    ## @end deftp
    LegendVisible = 'off';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} LegendTitle
    ##
    ## The title of the legend; from a table the name of
    ## @code{GroupVariable}, until set.
    ##
    ## @end deftp
    LegendTitle = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} Title
    ##
    ## The title of the chart.
    ##
    ## @end deftp
    Title = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} XLabel
    ##
    ## The label of the horizontal axis; from a table the name of
    ## @code{XVariable}, until set.
    ##
    ## @end deftp
    XLabel = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} YLabel
    ##
    ## The label of the vertical axis; from a table the name of
    ## @code{YVariable}, until set.
    ##
    ## @end deftp
    YLabel = '';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} FontName
    ##
    ## The font of every text on the chart, @qcode{'Helvetica'} by default.
    ##
    ## @end deftp
    FontName = 'Helvetica';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} FontSize
    ##
    ## The size of every text on the chart, in points, 10 by default.
    ##
    ## @end deftp
    FontSize = 10;

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} Position
    ##
    ## The position of the scatter plot in its parent, as
    ## @code{[left, bottom, width, height]} in @code{Units}.  While
    ## @code{PositionConstraint} is @qcode{'outerposition'} the chart places
    ## it within @code{OuterPosition} itself; setting it fixes the scatter
    ## plot there and turns @code{PositionConstraint} to
    ## @qcode{'innerposition'}.  @code{InnerPosition} is the same.
    ##
    ## @end deftp
    Position = [0.13, 0.11, 0.6, 0.6];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} OuterPosition
    ##
    ## The position of the whole chart in its parent, labels included,
    ## @code{[0, 0, 1, 1]} by default, or the outer position of the axes it
    ## took the place of.
    ##
    ## @end deftp
    OuterPosition = [0, 0, 1, 1];

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} PositionConstraint
    ##
    ## Which position the chart keeps: @qcode{'outerposition'}, the default,
    ## or @qcode{'innerposition'}.
    ##
    ## @end deftp
    PositionConstraint = 'outerposition';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} Units
    ##
    ## The units of the positions, @qcode{'normalized'} by default.
    ##
    ## @end deftp
    Units = 'normalized';

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} Visible
    ##
    ## Whether the chart is shown, @qcode{'on'} by default or @qcode{'off'}.
    ##
    ## @end deftp
    Visible = 'on';

  endproperties

  properties (Dependent)

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} InnerPosition
    ##
    ## The same as @code{Position}.
    ##
    ## @end deftp
    InnerPosition;

  endproperties

  properties (GetAccess = public, SetAccess = private)

    ## -*- texinfo -*-
    ## @deftp {stats.chart.ScatterHistogramChart} {property} Parent
    ##
    ## The figure or panel the chart is placed in.  This property is
    ## read-only.
    ##
    ## @end deftp
    Parent = [];

  endproperties

  properties (Access = private, Hidden)
    Axes_ = [];              # scatter plot, histogram of X, histogram of Y
    Legend_ = [];
    Drawn_ = false;
    Redrawing_ = false;
    Internal_ = false;       # true while the chart sets a value itself
    Reading_ = false;        # true while the table is being read
    BinsAuto_ = true;        # NumBins and BinWidths chosen by histcounts
    WidthsSet_ = false;      # BinWidths set, NumBins following from it
    LimitsAuto_ = [true, true];
    StyleAuto_ = true;       # colours and styles follow the groups
    LegendAuto_ = true;
    LabelsAuto_ = [true, true, true];   # XLabel, YLabel, LegendTitle
  endproperties

  methods (Hidden)

    function disp (this)
      printf ('  ScatterHistogramChart with properties:\n\n');
      if (isempty (this.SourceTable))
        printf ('%13s: [%dx1 double]\n', 'XData', numel (this.XData));
        printf ('%13s: [%dx1 double]\n', 'YData', numel (this.YData));
        if (! isempty (this.GroupData))
          printf ('%13s: [%dx1 %s]\n', 'GroupData', numel (this.GroupData), ...
                  class (this.GroupData));
        endif
      else
        printf ('%13s: [%dx%d table]\n', 'SourceTable', ...
                size (this.SourceTable));
        printf ('%13s: ''%s''\n', 'XVariable', this.XVariable);
        printf ('%13s: ''%s''\n', 'YVariable', this.YVariable);
        printf ('%13s: ''%s''\n', 'GroupVariable', this.GroupVariable);
      endif
      printf ('\n');
    endfunction

    function display (this)
      disp (this);
    endfunction

  endmethods

  methods (Hidden)

    function set.XData (this, val)
      if (this.Drawn_ && ! isempty (this.SourceTable))
        error (strcat ("stats.chart.ScatterHistogramChart: setting", ...
                       " 'XData' while 'SourceTable' holds a table is not", ...
                       " supported."));
      endif
      this.XData = shCheckData (val, 'XData');
      refit (this);
    endfunction

    function set.YData (this, val)
      if (this.Drawn_ && ! isempty (this.SourceTable))
        error (strcat ("stats.chart.ScatterHistogramChart: setting", ...
                       " 'YData' while 'SourceTable' holds a table is not", ...
                       " supported."));
      endif
      this.YData = shCheckData (val, 'YData');
      refit (this);
    endfunction

    function set.GroupData (this, val)
      if (! isempty (val) && ! isvector (val))
        error (strcat ("stats.chart.ScatterHistogramChart: 'GroupData'", ...
                       " must be a vector."));
      endif
      if (isempty (val))
        this.GroupData = [];
      else
        this.GroupData = val(:);
      endif
      refit (this);
    endfunction

    function set.SourceTable (this, val)
      if (! (isempty (val) || istable (val)))
        error (strcat ("stats.chart.ScatterHistogramChart: 'SourceTable'", ...
                       " must be a table."));
      endif
      this.SourceTable = val;
      readTable (this);
    endfunction

    function set.XVariable (this, val)
      this.XVariable = shCheckName (val, 'XVariable');
      readTable (this);
    endfunction

    function set.YVariable (this, val)
      this.YVariable = shCheckName (val, 'YVariable');
      readTable (this);
    endfunction

    function set.GroupVariable (this, val)
      this.GroupVariable = shCheckName (val, 'GroupVariable');
      readTable (this);
    endfunction

    function set.HistogramDisplayStyle (this, val)
      this.HistogramDisplayStyle = shCheckOneOf (val, {'stairs', 'bar', ...
                                                       'smooth'}, ...
                                                 'HistogramDisplayStyle');
      refit (this);
    endfunction

    function set.NumBins (this, val)
      if (this.Internal_)
        this.NumBins = val;
        return;
      endif
      this.NumBins = shCheckBins (val, numel (shGroups (this.GroupData)), ...
                                  'NumBins', true);
      this.BinsAuto_ = false;
      this.WidthsSet_ = false;
      refit (this);
    endfunction

    function set.BinWidths (this, val)
      if (this.Internal_)
        this.BinWidths = val;
        return;
      endif
      this.BinWidths = shCheckBins (val, numel (shGroups (this.GroupData)), ...
                                    'BinWidths', false);
      this.BinsAuto_ = false;
      this.WidthsSet_ = true;
      refit (this);
    endfunction

    function set.ScatterPlotLocation (this, val)
      list = {'SouthWest', 'SouthEast', 'NorthEast', 'NorthWest'};
      v = shCheckOneOf (val, list, 'ScatterPlotLocation');
      this.ScatterPlotLocation = list{strcmpi (v, list)};
      redraw (this);
    endfunction

    function set.ScatterPlotProportion (this, val)
      if (! (isnumeric (val) && isreal (val) && isscalar (val)
             && val > 0 && val < 1))
        error (strcat ("stats.chart.ScatterHistogramChart:", ...
                       " 'ScatterPlotProportion' must be a number between", ...
                       " 0 and 1."));
      endif
      this.ScatterPlotProportion = double (val);
      redraw (this);
    endfunction

    function set.XHistogramDirection (this, val)
      this.XHistogramDirection = shCheckOneOf (val, {'up', 'down'}, ...
                                               'XHistogramDirection');
      redraw (this);
    endfunction

    function set.YHistogramDirection (this, val)
      this.YHistogramDirection = shCheckOneOf (val, {'right', 'left'}, ...
                                               'YHistogramDirection');
      redraw (this);
    endfunction

    function set.XLimits (this, val)
      this.XLimits = shCheckLimits (val, 'XLimits');
      if (! this.Internal_)
        this.LimitsAuto_(1) = false;
        redraw (this);
      endif
    endfunction

    function set.YLimits (this, val)
      this.YLimits = shCheckLimits (val, 'YLimits');
      if (! this.Internal_)
        this.LimitsAuto_(2) = false;
        redraw (this);
      endif
    endfunction

    function set.Color (this, val)
      if (! (isnumeric (val) && isreal (val) && columns (val) == 3
             && rows (val) > 0 && all (val(:) >= 0) && all (val(:) <= 1)))
        error (strcat ("stats.chart.ScatterHistogramChart: 'Color' must be", ...
                       " an RGB triplet or a matrix of them."));
      endif
      this.Color = double (val);
      if (! this.Internal_)
        this.StyleAuto_ = false;
        redraw (this);
      endif
    endfunction

    function set.LineStyle (this, val)
      this.LineStyle = shCheckStyles (val, {'-', '--', ':', '-.', 'none'}, ...
                                      'LineStyle');
      if (! this.Internal_)
        this.StyleAuto_ = false;
        redraw (this);
      endif
    endfunction

    function set.LineWidth (this, val)
      this.LineWidth = shCheckPositive (val, 'LineWidth');
      redraw (this);
    endfunction

    function set.MarkerStyle (this, val)
      this.MarkerStyle = shCheckStyles (val, {'o', '+', '*', '.', 'x', ...
                                              's', 'd', '^', 'v', '>', '<', ...
                                              'p', 'h', 'none'}, ...
                                        'MarkerStyle');
      redraw (this);
    endfunction

    function set.MarkerSize (this, val)
      this.MarkerSize = shCheckPositive (val, 'MarkerSize');
      redraw (this);
    endfunction

    function set.MarkerFilled (this, val)
      this.MarkerFilled = shCheckOneOf (val, {'on', 'off'}, 'MarkerFilled');
      redraw (this);
    endfunction

    function set.MarkerAlpha (this, val)
      if (! (isnumeric (val) && isreal (val) && isscalar (val)
             && val >= 0 && val <= 1))
        error (strcat ("stats.chart.ScatterHistogramChart: 'MarkerAlpha'", ...
                       " must be a number from 0 to 1."));
      endif
      this.MarkerAlpha = double (val);
      redraw (this);
    endfunction

    function set.LegendVisible (this, val)
      this.LegendVisible = shCheckOneOf (val, {'on', 'off'}, 'LegendVisible');
      if (! this.Internal_)
        this.LegendAuto_ = false;
        redraw (this);
      endif
    endfunction

    function set.LegendTitle (this, val)
      this.LegendTitle = shCheckText (val, 'LegendTitle');
      if (! this.Internal_)
        this.LabelsAuto_(3) = false;
        redraw (this);
      endif
    endfunction

    function set.Title (this, val)
      this.Title = shCheckText (val, 'Title');
      redraw (this);
    endfunction

    function set.XLabel (this, val)
      this.XLabel = shCheckText (val, 'XLabel');
      if (! this.Internal_)
        this.LabelsAuto_(1) = false;
        redraw (this);
      endif
    endfunction

    function set.YLabel (this, val)
      this.YLabel = shCheckText (val, 'YLabel');
      if (! this.Internal_)
        this.LabelsAuto_(2) = false;
        redraw (this);
      endif
    endfunction

    function set.FontName (this, val)
      this.FontName = shCheckText (val, 'FontName');
      redraw (this);
    endfunction

    function set.FontSize (this, val)
      if (! (isnumeric (val) && isreal (val) && isscalar (val)
             && isfinite (val) && val > 0))
        error (strcat ("stats.chart.ScatterHistogramChart: 'FontSize'", ...
                       " must be a positive number."));
      endif
      this.FontSize = double (val);
      redraw (this);
    endfunction

    function set.Position (this, val)
      this.Position = shCheckPosition (val, 'Position');
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
      this.OuterPosition = shCheckPosition (val, 'OuterPosition');
      if (! this.Internal_)
        this.PositionConstraint = 'outerposition';
      endif
    endfunction

    function set.PositionConstraint (this, val)
      this.PositionConstraint = shCheckOneOf (val, {'outerposition', ...
                                                    'innerposition'}, ...
                                              'PositionConstraint');
      redraw (this);
    endfunction

    function set.Units (this, val)
      this.Units = shCheckOneOf (val, {'normalized', 'inches', ...
                                       'centimeters', 'points', 'pixels', ...
                                       'characters'}, 'Units');
      redraw (this);
    endfunction

    function set.Visible (this, val)
      this.Visible = shCheckOneOf (val, {'on', 'off'}, 'Visible');
      redraw (this);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {stats.chart.ScatterHistogramChart} {@var{obj} =} stats.chart.ScatterHistogramChart (@var{parent}, @var{spec}, @var{args})
    ##
    ## Create a @code{stats.chart.ScatterHistogramChart} object.
    ##
    ## @var{parent} is the figure or panel to place the chart in, or empty for
    ## the current figure, which is resolved only once every value has been
    ## accepted.  @var{spec} is a structure carrying the data as
    ## @code{scatterhistogram} resolved it, with the fields @qcode{XData},
    ## @qcode{YData}, @qcode{GroupData}, @qcode{SourceTable},
    ## @qcode{XVariable}, @qcode{YVariable} and @qcode{GroupVariable}, and
    ## @var{args} the name-value pairs left to set.  The documented way to
    ## reach this constructor is @code{scatterhistogram}.
    ##
    ## @seealso{scatterhistogram}
    ## @end deftypefn
    function this = ScatterHistogramChart (parent, spec, args)

      if (nargin < 2)
        error ("stats.chart.ScatterHistogramChart: too few input arguments.");
      endif
      if (nargin < 3)
        args = {};
      endif

      if (isempty (spec.SourceTable))
        this.XData = spec.XData;
        this.YData = spec.YData;
        this.GroupData = spec.GroupData;
      else
        this.SourceTable = spec.SourceTable;
        this.XVariable = spec.XVariable;
        this.YVariable = spec.YVariable;
        this.GroupVariable = spec.GroupVariable;
      endif

      ## Whatever was named is set now, so that one drawing covers them all;
      ## a value set here counts as chosen, as it would once drawn
      placed = false;
      for k = 1:2:numel (args)
        name = args{k};
        this.(name) = args{k+1};
        if (any (strcmp (name, {'Position', 'InnerPosition', ...
                                'OuterPosition'})))
          placed = true;
        endif
      endfor
      if (numel (this.XData) != numel (this.YData))
        error (strcat ("stats.chart.ScatterHistogramChart: 'XData' and", ...
                       " 'YData' must hold the same number of values."));
      endif
      if (! isempty (this.GroupData)
          && numel (this.GroupData) != numel (this.XData))
        error (strcat ("stats.chart.ScatterHistogramChart: 'GroupData'", ...
                       " must hold one value for each observation."));
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
          error (strcat ("stats.chart.ScatterHistogramChart: a scatter", ...
                         " histogram cannot be added to axes on which hold", ...
                         " is on."));
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
      tag = 'stats.chart.ScatterHistogramChart';
      this.Axes_ = [axes('parent', parent, 'tag', tag), ...
                    axes('parent', parent, 'tag', tag), ...
                    axes('parent', parent, 'tag', tag)];
      set (fig, 'currentaxes', this.Axes_(1));
      ## The chart goes with its axes, and its axes with the chart
      set (this.Axes_, 'deletefcn', @(~, ~) delete (this));
      this.Drawn_ = true;
      redraw (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {stats.chart.ScatterHistogramChart} {@var{limits} =} xlim (@var{h})
    ## @deftypefnx {stats.chart.ScatterHistogramChart} {} xlim (@var{h}, @var{limits})
    ##
    ## The limits of the horizontal axis, as @code{XLimits} holds them; given
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
    ## @deftypefn  {stats.chart.ScatterHistogramChart} {@var{limits} =} ylim (@var{h})
    ## @deftypefnx {stats.chart.ScatterHistogramChart} {} ylim (@var{h}, @var{limits})
    ##
    ## The limits of the vertical axis, as @code{YLimits} holds them; given
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
    ## @deftypefn {stats.chart.ScatterHistogramChart} {} delete (@var{h})
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
      for k = 1:numel (ax)
        if (ishghandle (ax(k)) && strcmp (get (ax(k), 'beingdeleted'), 'off'))
          set (ax(k), 'deletefcn', '');
          delete (ax(k));
        endif
      endfor
    endfunction

  endmethods

  methods (Hidden)

    ## Let go of the axes when something other than the chart clears the
    ## scatter plot, as a later plot does, leaving it to that plot and
    ## deleting the histograms.
    function release (this)
      if (this.Redrawing_ || isempty (this.Axes_))
        return;
      endif
      ax = this.Axes_;
      this.Axes_ = [];
      this.Parent = [];
      if (ishghandle (ax(1)))
        set (ax(1), 'deletefcn', '', 'tag', '');
      endif
      for k = 2:numel (ax)
        if (ishghandle (ax(k)))
          set (ax(k), 'deletefcn', '');
          delete (ax(k));
        endif
      endfor
      if (! isempty (this.Legend_) && ishghandle (this.Legend_))
        delete (this.Legend_);
      endif
      this.Legend_ = [];
    endfunction

  endmethods

  methods (Access = private)

    ## Read the samples from the table, with the labels they bring.
    function readTable (this)

      t = this.SourceTable;
      if (this.Reading_ || isempty (t) || isempty (this.XVariable)
          || isempty (this.YVariable))
        return;
      endif
      names = t.Properties.VariableNames;
      for v = {'XVariable', 'YVariable', 'GroupVariable'}
        n = this.(v{1});
        if (! isempty (n) && ! any (strcmp (names, n)))
          error (strcat ("stats.chart.ScatterHistogramChart: '%s' does not", ...
                         " name a variable of the source table."), v{1});
        endif
      endfor
      this.Reading_ = true;
      drawn = this.Drawn_;
      this.Drawn_ = false;
      unwind_protect
        this.XData = shCheckData (t.(this.XVariable), 'XVariable');
        this.YData = shCheckData (t.(this.YVariable), 'YVariable');
        if (isempty (this.GroupVariable))
          this.GroupData = [];
        else
          this.GroupData = t.(this.GroupVariable);
        endif
        this.Internal_ = true;
        if (this.LabelsAuto_(1))
          this.XLabel = this.XVariable;
        endif
        if (this.LabelsAuto_(2))
          this.YLabel = this.YVariable;
        endif
        if (this.LabelsAuto_(3))
          this.LegendTitle = this.GroupVariable;
        endif
      unwind_protect_cleanup
        this.Internal_ = false;
        this.Drawn_ = drawn;
        this.Reading_ = false;
      end_unwind_protect
      refit (this);

    endfunction

    ## Bring the bins, limits, styles and legend that follow the data in line
    ## with it, then draw.
    function refit (this)

      if (this.Reading_)
        return;
      endif
      x = this.XData;
      y = this.YData;
      [names, g] = shGroups (this.GroupData);
      ng = max (numel (names), 1);
      this.Internal_ = true;
      unwind_protect
        ## Limits: the range of each sample
        if (this.LimitsAuto_(1))
          this.XLimits = shRange (x);
        endif
        if (this.LimitsAuto_(2))
          this.YLimits = shRange (y);
        endif
        ## Colours and line styles, one for each group
        if (this.StyleAuto_)
          co = get (groot, 'defaultaxescolororder');
          this.Color = co(mod ((1:ng) - 1, rows (co)) + 1,:);
          ls = {'-', ':', '-.', '--'};
          styles = ls(mod ((1:ng) - 1, 4) + 1);
          if (ng == 1)
            this.LineStyle = styles{1};
          else
            this.LineStyle = styles;
          endif
        endif
        if (this.LegendAuto_)
          this.LegendVisible = shOnOff (! isempty (names));
        endif
        ## Bins, for each variable and each group, as histcounts chooses them
        ## or from the numbers or widths set
        if (strcmp (this.HistogramDisplayStyle, 'smooth'))
          this.NumBins = [];
          this.BinWidths = [];
        elseif (numel (x) == numel (y))
          [nb, bw] = shBins (x, y, g, ng, this.NumBins, this.BinWidths, ...
                             this.BinsAuto_, this.WidthsSet_);
          this.NumBins = nb;
          this.BinWidths = bw;
        endif
      unwind_protect_cleanup
        this.Internal_ = false;
      end_unwind_protect
      redraw (this);

    endfunction

    ## Draw the chart from scratch into its axes.
    function redraw (this)

      if (! this.Drawn_ || isempty (this.Axes_)
          || ! all (ishghandle (this.Axes_)))
        return;
      endif
      ax = this.Axes_;
      this.Redrawing_ = true;
      unwind_protect
        if (! isempty (this.Legend_) && ishghandle (this.Legend_))
          delete (this.Legend_);
        endif
        this.Legend_ = [];
        for k = 1:3
          delete (get (ax(k), 'children'));
        endfor
      unwind_protect_cleanup
        this.Redrawing_ = false;
      end_unwind_protect
      ## A child whose deletion by anyone else releases the axes
      text (ax(1), NaN, NaN, '', 'deletefcn', @(~, ~) release (this));

      x = this.XData;
      y = this.YData;
      [names, g] = shGroups (this.GroupData);
      ng = max (numel (names), 1);
      if (isempty (names))
        g = ones (numel (x), 1);
      endif
      vis = this.Visible;
      set (ax, 'units', this.Units, 'visible', vis, 'fontname', ...
           this.FontName, 'fontsize', this.FontSize, 'nextplot', 'add');

      ## The scatter plot, one set of markers for each group
      hs = zeros (1, ng);
      for k = 1:ng
        take = (g == k);
        c = shPick (this.Color, k);
        m = shPickStyle (this.MarkerStyle, k);
        sz = shPick (this.MarkerSize(:), k);
        if (strcmp (this.MarkerFilled, 'on'))
          hs(k) = scatter (ax(1), x(take), y(take), sz, c, m, 'filled');
        else
          hs(k) = scatter (ax(1), x(take), y(take), sz, c, m);
        endif
        set (hs(k), 'markerfacealpha', this.MarkerAlpha, 'markeredgealpha', ...
             this.MarkerAlpha, 'visible', vis);
        if (! isempty (names))
          set (hs(k), 'displayname', names{k});
        endif
      endfor
      set (ax(1), 'xlim', this.XLimits, 'ylim', this.YLimits, 'box', 'on');
      xlabel (ax(1), this.XLabel, 'fontname', this.FontName);
      ylabel (ax(1), this.YLabel, 'fontname', this.FontName);

      ## The histograms, as densities, one outline, set of bars or curve for
      ## each group
      style = this.HistogramDisplayStyle;
      for k = 1:ng
        take = (g == k);
        c = shPick (this.Color, k);
        ls = shPickStyle (this.LineStyle, k);
        lw = shPick (this.LineWidth(:), k);
        for d = 1:2
          if (d == 1)
            v = x(take);
            lims = this.XLimits;
          else
            v = y(take);
            lims = this.YLimits;
          endif
          v = v(isfinite (v));
          if (isempty (v))
            continue;
          endif
          if (strcmp (style, 'smooth'))
            [f, xi] = ksdensity (v);
            px = xi(:)';
            py = f(:)';
          else
            w = this.BinWidths(d, min (k, columns (this.BinWidths)));
            [n, e] = histcounts (v, 'BinWidth', w);
            f = n / (numel (v) * w);
            [px, py] = shStairs (e, f);
          endif
          hax = ax(d + 1);
          if (d == 2)
            [px, py] = deal (py, px);
          endif
          if (strcmp (style, 'bar'))
            if (d == 1)
              [bx, by] = shBars (e, f);
            else
              [by, bx] = shBars (e, f);
            endif
            patch (hax, bx, by, c, 'facealpha', 0.5, 'edgecolor', c, ...
                   'linestyle', ls, 'linewidth', lw, 'visible', vis);
          else
            line (hax, px, py, 'color', c, 'linestyle', ls, 'linewidth', lw, ...
                  'visible', vis);
          endif
        endfor
      endfor

      ## Where each axes goes, the histograms beside the scatter plot and
      ## sharing its limits
      [ps, px, py, pl] = layout (this);
      set (ax(1), 'position', ps);
      set (ax(2), 'position', px, 'xlim', this.XLimits, 'xticklabel', {}, ...
           'box', 'off');
      set (ax(3), 'position', py, 'ylim', this.YLimits, 'yticklabel', {}, ...
           'box', 'off');
      top = any (strcmp (this.ScatterPlotLocation, {'SouthWest', 'SouthEast'}));
      right = any (strcmp (this.ScatterPlotLocation, {'SouthWest', ...
                                                      'NorthWest'}));
      if (top)
        set (ax(2), 'xaxislocation', 'bottom');
      else
        set (ax(2), 'xaxislocation', 'top');
      endif
      if (right)
        set (ax(3), 'yaxislocation', 'left');
      else
        set (ax(3), 'yaxislocation', 'right');
      endif
      if (strcmp (this.XHistogramDirection, 'down'))
        set (ax(2), 'ydir', 'reverse');
      else
        set (ax(2), 'ydir', 'normal');
      endif
      if (strcmp (this.YHistogramDirection, 'left'))
        set (ax(3), 'xdir', 'reverse');
      else
        set (ax(3), 'xdir', 'normal');
      endif
      set (ax(2:3), 'color', 'none');
      set (ax(2), 'ytick', [], 'ycolor', 'none');
      set (ax(3), 'xtick', [], 'xcolor', 'none');

      ## The title over the whole chart
      title (ax(1), '');
      title (ax(2), '');
      if (top)
        title (ax(2), this.Title, 'fontname', this.FontName);
      else
        title (ax(1), this.Title, 'fontname', this.FontName);
      endif

      ## Each group was added to the one before; the axes are left as plain
      ## axes are, so a later plot replaces them rather than adding to them
      set (ax, 'nextplot', 'replace');

      ## The legend of the groups, in the corner the histograms leave empty
      if (strcmp (this.LegendVisible, 'on') && ! isempty (names))
        this.Legend_ = legend (ax(1), hs, names(:)');
        set (this.Legend_, 'units', this.Units, 'position', pl, ...
             'visible', vis, 'fontname', this.FontName);
        if (! isempty (this.LegendTitle))
          title (this.Legend_, this.LegendTitle);
        endif
      endif

    endfunction

    ## The positions of the scatter plot, the two histograms and the legend.
    function [ps, px, py, pl] = layout (this)

      p = this.ScatterPlotProportion;
      if (strcmp (this.PositionConstraint, 'outerposition'))
        o = this.OuterPosition;
        ti = get (this.Axes_(1), 'tightinset');
        pad = 0.02;
        L = max (0.1 * o(3), ti(1) + pad * o(3));
        B = max (0.1 * o(4), ti(2) + pad * o(4));
        R = 0.05 * o(3);
        T = 0.08 * o(4);
        box = [o(1) + L, o(2) + B, o(3) - L - R, o(4) - B - T];
        box(3:4) = max (box(3:4), 0.01);
        ps = [0, 0, p * box(3), p * box(4)];
      else
        ps = this.Position;
        box = [0, 0, ps(3) / p, ps(4) / p];
      endif
      gap = 0.01;
      w = box(3) * (1 - p) - gap;
      h = box(4) * (1 - p) - gap;
      loc = this.ScatterPlotLocation;
      west = any (strcmp (loc, {'SouthWest', 'NorthWest'}));
      south = any (strcmp (loc, {'SouthWest', 'SouthEast'}));
      if (strcmp (this.PositionConstraint, 'outerposition'))
        ps(1) = box(1) + ! west * (box(3) - ps(3));
        ps(2) = box(2) + ! south * (box(4) - ps(4));
      else
        box(1) = ps(1) - ! west * (box(3) - ps(3));
        box(2) = ps(2) - ! south * (box(4) - ps(4));
      endif
      if (south)
        px = [ps(1), ps(2) + ps(4) + gap, ps(3), h];
      else
        px = [ps(1), ps(2) - gap - h, ps(3), h];
      endif
      if (west)
        py = [ps(1) + ps(3) + gap, ps(2), w, ps(4)];
      else
        py = [ps(1) - gap - w, ps(2), w, ps(4)];
      endif
      pl = [py(1), px(2), w, h];
      this.Internal_ = true;
      unwind_protect
        this.Position = ps;
      unwind_protect_cleanup
        this.Internal_ = false;
      end_unwind_protect

    endfunction

  endmethods

endclassdef

## The observations, as a numeric column
function v = shCheckData (val, name)
  if (! ((isnumeric (val) || islogical (val)) && isreal (val)
         && (isvector (val) || isempty (val))))
    error (strcat ("stats.chart.ScatterHistogramChart: '%s' must be a", ...
                   " real numeric vector."), name);
  endif
  v = double (val(:));
endfunction

## The names of the groups in order of first appearance, and the group of
## each observation, 0 where it is missing
function [names, g] = shGroups (gd)
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
    error (strcat ("stats.chart.ScatterHistogramChart: 'GroupData' must", ...
                   " hold categories, text, numbers or logical values."));
  endif
  t = t(:);
  [names, i] = unique (t(! miss), 'first');
  [~, o] = sort (i);
  names = names(o);
  [~, g(! miss)] = ismember (t(! miss), names);
endfunction

## Bins for each variable and group: what histcounts chooses, or from the
## numbers or widths set, one column for each group
function [nb, bw] = shBins (x, y, g, ng, nbSet, bwSet, auto, widths)
  nb = zeros (2, ng);
  bw = zeros (2, ng);
  if (all (g == 0))
    g = ones (numel (x), 1);
  endif
  for k = 1:ng
    for d = 1:2
      if (d == 1)
        v = x(g == k);
      else
        v = y(g == k);
      endif
      v = v(isfinite (v));
      if (isempty (v))
        continue;
      endif
      if (auto)
        [~, e] = histcounts (v);
      elseif (widths)
        [~, e] = histcounts (v, 'BinWidth', bwSet(d, min (k, columns (bwSet))));
      else
        [~, e] = histcounts (v, nbSet(d, min (k, columns (nbSet))));
      endif
      nb(d,k) = numel (e) - 1;
      bw(d,k) = e(2) - e(1);
    endfor
  endfor
endfunction

## The outline of a histogram with edges E and heights F
function [px, py] = shStairs (e, f)
  n = numel (f);
  px = reshape ([e(1:n); e(2:n+1)], 1, []);
  py = reshape ([f; f], 1, []);
  px = [e(1), px, e(end)];
  py = [0, py, 0];
endfunction

## The patches of the bars of a histogram with edges E and heights F
function [bx, by] = shBars (e, f)
  n = numel (f);
  bx = [e(1:n); e(2:n+1); e(2:n+1); e(1:n)];
  by = [zeros(1, n); zeros(1, n); f(:)'; f(:)'];
endfunction

## The range of the finite values, or [0, 1]
function lims = shRange (v)
  v = v(isfinite (v));
  if (isempty (v))
    lims = [0, 1];
  elseif (min (v) == max (v))
    lims = min (v) + [-0.5, 0.5];
  else
    lims = [min(v), max(v)];
  endif
endfunction

## Row K of a matrix of values for the groups, the last one repeated
function v = shPick (vals, k)
  v = vals(min (k, rows (vals)),:);
endfunction

## Style K of a character vector or a cell array of styles
function s = shPickStyle (styles, k)
  if (iscell (styles))
    s = styles{min (k, numel (styles))};
  else
    s = styles;
  endif
endfunction

## 'on' or 'off'
function s = shOnOff (tf)
  if (tf)
    s = 'on';
  else
    s = 'off';
  endif
endfunction

## Numbers or widths of bins: a positive scalar, a 2-by-1 vector or a
## 2-by-N matrix for N groups
function v = shCheckBins (val, ng, name, integer)
  ok = isnumeric (val) && isreal (val) && ! isempty (val) ...
       && all (isfinite (val(:))) && all (val(:) > 0);
  if (ok && integer)
    ok = all (val(:) == fix (val(:)));
  endif
  if (ok && isscalar (val))
    v = [val; val];
  elseif (ok && isequal (size (val), [2, 1]))
    v = val;
  elseif (ok && rows (val) == 2 && columns (val) == max (ng, 1))
    v = val;
  else
    error (strcat ("stats.chart.ScatterHistogramChart: '%s' must be a", ...
                   " positive scalar, a 2-by-1 vector or a 2-by-N matrix", ...
                   " for N groups."), name);
  endif
  v = double (v);
endfunction

## Limits: a 1-by-2 increasing vector
function v = shCheckLimits (val, name)
  if (! (isnumeric (val) && isreal (val) && numel (val) == 2
         && all (isfinite (val)) && val(2) > val(1)))
    error (strcat ("stats.chart.ScatterHistogramChart: '%s' must be a", ...
                   " 1-by-2 increasing numeric vector."), name);
  endif
  v = double (val(:)');
endfunction

## Styles: one of LIST, or a cell array or string array of them
function v = shCheckStyles (val, list, name)
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
    error ("stats.chart.ScatterHistogramChart: '%s' must be one of %s.", ...
           name, strjoin (strcat ("'", list, "'"), ', '));
  endif
endfunction

## Positive values, a scalar or one for each group
function v = shCheckPositive (val, name)
  if (! (isnumeric (val) && isreal (val) && isvector (val)
         && all (isfinite (val)) && all (val > 0)))
    error (strcat ("stats.chart.ScatterHistogramChart: '%s' must be", ...
                   " positive."), name);
  endif
  v = double (val(:)');
endfunction

## Variable names: a character vector or a string scalar
function out = shCheckName (val, name)
  if (isa (val, 'string') && isscalar (val))
    val = char (val);
  endif
  if (! (ischar (val) && (isrow (val) || isempty (val))))
    error (strcat ("stats.chart.ScatterHistogramChart: '%s' must be a", ...
                   " character vector."), name);
  endif
  out = val;
endfunction

## Text: a character vector, a string, or a cell array of them for lines
function out = shCheckText (val, name)
  if (isa (val, 'string'))
    val = cellstr (val);
    if (numel (val) == 1)
      val = val{1};
    endif
  endif
  if (! ((ischar (val) && (isrow (val) || isempty (val))) || iscellstr (val)))
    error ("stats.chart.ScatterHistogramChart: '%s' must be text.", name);
  endif
  out = val;
endfunction

## One of a list of names, case free
function out = shCheckOneOf (val, list, name)
  if (isa (val, 'string') && isscalar (val))
    val = char (val);
  endif
  if (! (ischar (val) && any (strcmpi (val, list))))
    error ("stats.chart.ScatterHistogramChart: '%s' must be one of %s.", ...
           name, strjoin (strcat ("'", list, "'"), ', '));
  endif
  out = lower (val);
endfunction

## A position, as a 1-by-4 vector with positive width and height
function out = shCheckPosition (val, name)
  if (! (isnumeric (val) && isreal (val) && numel (val) == 4
         && all (isfinite (val)) && val(3) > 0 && val(4) > 0))
    error (strcat ("stats.chart.ScatterHistogramChart: '%s' must be a", ...
                   " 1-by-4 vector with positive width and height."), name);
  endif
  out = double (val(:)');
endfunction

%!shared x, y, g
%! x = [2.1; 3.4; 1.9; 5.6; 4.4; 3.8; 2.7; 6.1; 3.3; 4.9; 1.2; 2.8];
%! y = [1.2; 2.8; 0.9; 2.2; 3.1; 1.7; 2.5; 0.4; 1.9; 3.6; 2.0; 1.1];
%! g = categorical ({'a';'b';'a';'b';'a';'b';'a';'b';'a';'c';'c';'c'});
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y);
%!   assert_equal (class (h), 'stats.chart.ScatterHistogramChart');
%!   assert_equal (h.XData, x);
%!   assert_equal (h.YData, y);
%!   assert_equal (h.GroupData, []);
%!   assert_equal (h.NumBins, [4; 4]);
%!   assert_equal (h.BinWidths, [2; 1]);
%!   assert_equal (h.XLimits, [1.2, 6.1]);
%!   assert_equal (h.YLimits, [0.4, 3.6]);
%!   assert_equal (h.Color, [0, 0.447, 0.741]);
%!   assert_equal (h.HistogramDisplayStyle, 'stairs');
%!   assert_equal (h.ScatterPlotLocation, 'SouthWest');
%!   assert_equal (h.ScatterPlotProportion, 0.75);
%!   assert_equal (h.XHistogramDirection, 'up');
%!   assert_equal (h.YHistogramDirection, 'right');
%!   assert_equal (h.LineStyle, '-');
%!   assert_equal (h.LineWidth, 0.5);
%!   assert_equal (h.MarkerStyle, 'o');
%!   assert_equal (h.MarkerSize, 36);
%!   assert_equal (h.MarkerFilled, 'on');
%!   assert_equal (h.MarkerAlpha, 1);
%!   assert_equal (h.LegendVisible, 'off');
%!   assert_equal (h.OuterPosition, [0, 0, 1, 1]);
%!   assert_equal (h.Parent, hf);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y, 'GroupData', g);
%!   assert_equal (h.LegendVisible, 'on');
%!   assert_equal (h.Color, [0, 0.447, 0.741; 0.85, 0.325, 0.098; ...
%!                           0.929, 0.694, 0.125]);
%!   assert_equal (h.LineStyle, {'-', ':', '-.'});
%!   assert_equal (h.NumBins, [3, 2, 1; 2, 2, 2]);
%!   assert_equal (h.BinWidths, [2, 3, 5; 2, 2, 3]);
%!   tag = 'stats.chart.ScatterHistogramChart';
%!   lg = findall (hf, 'type', 'axes', 'tag', 'legend');
%!   assert_equal (get (lg, 'string'), {'a'; 'b'; 'c'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y, 'GroupData', {'q';'p';'q';'p';'q';'p'; ...
%!                                            'q';'p';'q';'p';'q';'p'});
%!   assert_equal (rows (h.Color), 2);
%!   lg = findall (hf, 'type', 'axes', 'tag', 'legend');
%!   assert_equal (get (lg, 'string'), {'q'; 'p'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y, 'NumBins', 4);
%!   assert_equal (h.NumBins, [4; 4]);
%!   assert_equal (h.BinWidths, [1.3; 0.9], -1e-12);
%!   h.NumBins = [3; 5];
%!   assert_equal (h.NumBins, [3; 5]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y, 'NumBins', 4, 'BinWidths', 1);
%!   assert_equal (h.NumBins, [6; 4]);
%!   assert_equal (h.BinWidths, [1; 1]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y, 'HistogramDisplayStyle', 'smooth');
%!   assert_equal (h.NumBins, []);
%!   assert_equal (h.BinWidths, []);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! ## Each histogram is a density: the area under it is one
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y);
%!   tag = 'stats.chart.ScatterHistogramChart';
%!   ax = findall (hf, 'type', 'axes', 'tag', tag);
%!   ln = [findall(ax, 'type', 'line'); findall(ax, 'type', 'patch')];
%!   for k = 1:numel (ln)
%!     px = get (ln(k), 'xdata');
%!     py = get (ln(k), 'ydata');
%!     if (numel (px) > 2)
%!       area = max (abs (trapz (px(:), py(:))), abs (trapz (py(:), px(:))));
%!       assert_equal (area, 1, -1e-12);
%!     endif
%!   endfor
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y);
%!   h.XLimits = [0, 10];
%!   assert_equal (h.XLimits, [0, 10]);
%!   assert_equal (xlim (h), [0, 10]);
%!   ylim (h, [-1, 5]);
%!   assert_equal (h.YLimits, [-1, 5]);
%!   h.XData = x + 100;
%!   assert_equal (h.XLimits, [0, 10]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram ([x; NaN], [y; 2]);
%!   assert_equal (numel (h.XData), 13);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! ## The histograms line the sides away from the scatter plot's corner
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y);
%!   tag = 'stats.chart.ScatterHistogramChart';
%!   ax = findall (hf, 'type', 'axes', 'tag', tag);
%!   ps = h.Position;
%!   p = cell2mat (get (ax, 'position'));
%!   top = p(:,2) > ps(2) + ps(4) - 1e-9;
%!   right = p(:,1) > ps(1) + ps(3) - 1e-9;
%!   assert_equal ([sum(top), sum(right)], [1, 1]);
%!   h.ScatterPlotLocation = 'northeast';
%!   assert_equal (h.ScatterPlotLocation, 'NorthEast');
%!   ps = h.Position;
%!   p = cell2mat (get (ax, 'position'));
%!   below = p(:,2) + p(:,4) < ps(2) + 1e-9;
%!   left = p(:,1) + p(:,3) < ps(1) + 1e-9;
%!   assert_equal ([sum(below), sum(left)], [1, 1]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y, 'XHistogramDirection', 'down', ...
%!                         'YHistogramDirection', 'left', ...
%!                         'MarkerFilled', 'off', 'MarkerAlpha', 0.5, ...
%!                         'Title', 'T', 'XLabel', 'XL', 'YLabel', 'YL');
%!   assert_equal (h.XHistogramDirection, 'down');
%!   assert_equal (h.YHistogramDirection, 'left');
%!   assert_equal (h.MarkerFilled, 'off');
%!   assert_equal (h.MarkerAlpha, 0.5);
%!   assert_equal (h.Title, 'T');
%!   assert_equal (get (get (gca (), 'xlabel'), 'string'), 'XL');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = scatterhistogram (x, y, 'GroupData', g, 'Color', ...
%!                         [1, 0, 0; 0, 1, 0; 0, 0, 1], 'LineStyle', ...
%!                         {'-', '--', ':'}, 'MarkerStyle', {'o', 'x', 's'}, ...
%!                         'MarkerSize', [4, 6, 8], 'LineWidth', [1, 2, 3]);
%!   assert_equal (h.Color, [1, 0, 0; 0, 1, 0; 0, 0, 1]);
%!   assert_equal (h.LineStyle, {'-', '--', ':'});
%!   assert_equal (h.MarkerStyle, {'o', 'x', 's'});
%!   assert_equal (h.MarkerSize, [4, 6, 8]);
%!   assert_equal (h.LineWidth, [1, 2, 3]);
%!   h.GroupData = g(end:-1:1);
%!   assert_equal (h.Color, [1, 0, 0; 0, 1, 0; 0, 0, 1]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!shared shS
%! shS = struct ('XData', [1; 2; 3], 'YData', [4; 5; 6], 'GroupData', [], ...
%!               'SourceTable', [], 'XVariable', '', 'YVariable', '', ...
%!               'GroupVariable', '');
%!error<stats.chart.ScatterHistogramChart: too few input arguments.> ...
%! stats.chart.ScatterHistogramChart ([])
%!error<stats.chart.ScatterHistogramChart: 'XData' must be a real numeric vector.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'XData', {1, 2}})
%!error<stats.chart.ScatterHistogramChart: 'XData' and 'YData' must hold the same number of values.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'YData', [1; 2]})
%!error<stats.chart.ScatterHistogramChart: 'GroupData' must hold one value for each observation.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'GroupData', [1; 2]})
%!error<stats.chart.ScatterHistogramChart: 'GroupData' must be a vector.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'GroupData', ones(2, 2)})
%!error<stats.chart.ScatterHistogramChart: 'SourceTable' must be a table.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'SourceTable', 5})
%!error<stats.chart.ScatterHistogramChart: 'HistogramDisplayStyle' must be one of 'stairs', 'bar', 'smooth'.> ...
%! stats.chart.ScatterHistogramChart ([], shS, ...
%!                                   {'HistogramDisplayStyle', 'area'})
%!error<stats.chart.ScatterHistogramChart: 'NumBins' must be a positive scalar, a 2-by-1 vector or a 2-by-N matrix for N groups.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'NumBins', 0})
%!error<stats.chart.ScatterHistogramChart: 'NumBins' must be a positive scalar, a 2-by-1 vector or a 2-by-N matrix for N groups.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'NumBins', [3, 5]})
%!error<stats.chart.ScatterHistogramChart: 'BinWidths' must be a positive scalar, a 2-by-1 vector or a 2-by-N matrix for N groups.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'BinWidths', [1, 0.5]})
%!error<stats.chart.ScatterHistogramChart: 'ScatterPlotLocation' must be one of 'SouthWest', 'SouthEast', 'NorthEast', 'NorthWest'.> ...
%! stats.chart.ScatterHistogramChart ([], shS, ...
%!                                   {'ScatterPlotLocation', 'Middle'})
%!error<stats.chart.ScatterHistogramChart: 'ScatterPlotProportion' must be a number between 0 and 1.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'ScatterPlotProportion', 1.2})
%!error<stats.chart.ScatterHistogramChart: 'XHistogramDirection' must be one of 'up', 'down'.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'XHistogramDirection', 'left'})
%!error<stats.chart.ScatterHistogramChart: 'YHistogramDirection' must be one of 'right', 'left'.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'YHistogramDirection', 'up'})
%!error<stats.chart.ScatterHistogramChart: 'XLimits' must be a 1-by-2 increasing numeric vector.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'XLimits', [2, 1]})
%!error<stats.chart.ScatterHistogramChart: 'Color' must be an RGB triplet or a matrix of them.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'Color', 'r'})
%!error<stats.chart.ScatterHistogramChart: 'LineStyle' must be one of '-', '--', ':', '-.', 'none'.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'LineStyle', '*'})
%!error<stats.chart.ScatterHistogramChart: 'LineWidth' must be positive.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'LineWidth', 0})
## test.m ends a pattern at the first '>' of its line, and the list of
## markers holds one, so the message is matched only as far as 'v'
%!error<stats.chart.ScatterHistogramChart: 'MarkerStyle' must be one of 'o', '\+', '\*', '.', 'x', 's', 'd', '\^', 'v', > ...
%! stats.chart.ScatterHistogramChart ([], shS, {'MarkerStyle', 'q'})
%!error<stats.chart.ScatterHistogramChart: 'MarkerFilled' must be one of 'on', 'off'.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'MarkerFilled', 1})
%!error<stats.chart.ScatterHistogramChart: 'MarkerAlpha' must be a number from 0 to 1.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'MarkerAlpha', 2})
%!error<stats.chart.ScatterHistogramChart: 'LegendVisible' must be one of 'on', 'off'.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'LegendVisible', 'yes'})
%!error<stats.chart.ScatterHistogramChart: 'Title' must be text.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'Title', 5})
%!error<stats.chart.ScatterHistogramChart: 'FontSize' must be a positive number.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'FontSize', -1})
%!error<stats.chart.ScatterHistogramChart: 'Position' must be a 1-by-4 vector with positive width and height.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'Position', [0, 0, 1]})
%!error<stats.chart.ScatterHistogramChart: 'Visible' must be one of 'on', 'off'.> ...
%! stats.chart.ScatterHistogramChart ([], shS, {'Visible', 'maybe'})
%!error<stats.chart.ScatterHistogramChart: 'YVariable' does not name a variable of the source table.> ...
%! stats.chart.ScatterHistogramChart ([], shS, ...
%!   {'SourceTable', table([1; 2], 'VariableNames', {'X'}), ...
%!    'XVariable', 'X', 'YVariable', 'Q'})
