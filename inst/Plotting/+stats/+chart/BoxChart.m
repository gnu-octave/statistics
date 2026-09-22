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
## @deftp {statistics} stats.chart.BoxChart
##
## A box chart, as @code{boxchart} draws it.
##
## A @code{stats.chart.BoxChart} holds the data a box chart is drawn from and
## the choices it is drawn with, and redraws itself whenever one of them is
## set.  It is what @code{boxchart} returns, and the documented way to reach
## it; one is returned per colour group where @qcode{'GroupByColor'} names
## one.
##
## The data comes either from vectors or from a table, never from both, and
## the object says which through @code{XDataMode} and @code{YDataMode}.
## Where the data came from vectors these read @qcode{'manual'} and
## @code{XVariable} and @code{YVariable} are empty; where it came from a
## table they read @qcode{'auto'} and naming a variable is what chooses the
## column.  Setting @code{YData} on an object holding a table is refused, and
## so is naming a variable on one holding vectors, as MATLAB R2026a refuses
## them.  Assigning @code{SourceTable} again re-reads the columns, so a table
## changed after the chart was drawn is picked up.
##
## The box spans the lower and upper quartiles with the median across it, the
## whiskers reach the furthest values within one and a half interquartile
## ranges of the box, and whatever lies beyond is drawn as a marker.  These
## are the quartiles @code{quantile} computes and the rule
## @code{isoutlier} applies with @qcode{'quartiles'}.
##
## MATLAB's object carries the properties every graphics object carries,
## @code{Annotation} and @code{DataTipTemplate} among them.  This one is a
## handle class drawing with ordinary primitives rather than a graphics
## object of its own, so it carries what belongs to a box chart and leaves
## the rest to the primitives it draws.
##
## @seealso{boxchart, swarmchart, boxplot}
## @end deftp

classdef BoxChart < handle

  properties (Access = public)

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} BoxWidth
    ##
    ## Width of each box, in the units of the horizontal axis.  The default
    ## is 0.5.
    ##
    ## @end deftp
    BoxWidth = 0.5;

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} BoxFaceColor
    ##
    ## Colour the boxes are filled with, as an RGB triplet.  Setting it puts
    ## @code{BoxFaceColorMode} at @qcode{'manual'}.
    ##
    ## @end deftp
    BoxFaceColor = [0, 0.447, 0.741];

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} BoxFaceAlpha
    ##
    ## How opaque the box fill is, from 0 to 1.  The default is 0.2.
    ##
    ## @end deftp
    BoxFaceAlpha = 0.2;

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} BoxEdgeColor
    ##
    ## Colour of the box outline.  Setting it puts
    ## @code{BoxEdgeColorMode} at @qcode{'manual'}.
    ##
    ## @end deftp
    BoxEdgeColor = [0, 0.447, 0.741];

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} BoxMedianLineColor
    ##
    ## Colour of the line marking the median.  Setting it puts
    ## @code{BoxMedianLineColorMode} at @qcode{'manual'}.
    ##
    ## @end deftp
    BoxMedianLineColor = [0, 0.447, 0.741];

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} WhiskerLineColor
    ##
    ## Colour of the whiskers and their caps.
    ##
    ## @end deftp
    WhiskerLineColor = [0, 0.447, 0.741];

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} WhiskerLineStyle
    ##
    ## Line style of the whiskers: @qcode{'-'} by default, or @qcode{'--'},
    ## @qcode{':'}, @qcode{'-.'} or @qcode{'none'}.
    ##
    ## @end deftp
    WhiskerLineStyle = '-';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} CapWidth
    ##
    ## Width of the cap at the end of each whisker, in the units of the
    ## horizontal axis.  The default is 0.25.
    ##
    ## @end deftp
    CapWidth = 0.25;

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} LineWidth
    ##
    ## Width of the box and whisker lines.  The default is 1.
    ##
    ## @end deftp
    LineWidth = 1;

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} MarkerStyle
    ##
    ## Marker an outlier is drawn with, @qcode{'o'} by default, or
    ## @qcode{'none'} to draw none.
    ##
    ## @end deftp
    MarkerStyle = 'o';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} MarkerSize
    ##
    ## Size of an outlier marker.  The default is 6.
    ##
    ## @end deftp
    MarkerSize = 6;

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} MarkerColor
    ##
    ## Colour of an outlier marker.  Setting it puts
    ## @code{MarkerColorMode} at @qcode{'manual'}.
    ##
    ## @end deftp
    MarkerColor = [0, 0.447, 0.741];

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} Notch
    ##
    ## Whether the box is notched about its median, @qcode{'off'} by default.
    ## A notch spans the median plus and minus
    ## @math{1.57 IQR / sqrt (n)}, so two boxes whose notches do not overlap
    ## have medians that differ at roughly the five per cent level.
    ##
    ## @end deftp
    Notch = 'off';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} JitterOutliers
    ##
    ## Whether outlier markers are spread across the width of the box rather
    ## than drawn in a line, @qcode{'off'} by default.
    ##
    ## @end deftp
    JitterOutliers = 'off';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} Orientation
    ##
    ## Whether the boxes stand @qcode{'vertical'}, the default, or lie
    ## @qcode{'horizontal'}.
    ##
    ## @end deftp
    Orientation = 'vertical';

  endproperties

  properties (Access = public)

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} XData
    ##
    ## Where each observation's box stands, one value per observation: a
    ## grouping variable, numeric or categorical.  Where no grouping was
    ## given every observation shares one position and this is a
    ## @code{categorical} of one level.  Setting it is refused while
    ## @code{XVariable} names a column, the data coming from the table
    ## instead.
    ##
    ## @end deftp
    XData = [];

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} YData
    ##
    ## The observations themselves, one value per observation.  Setting it is
    ## refused while @code{YVariable} names a column.
    ##
    ## @end deftp
    YData = [];

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} SourceTable
    ##
    ## The table the data is read from, empty where it came from vectors.
    ## Assigning it again re-reads the columns @code{XVariable} and
    ## @code{YVariable} name, so a table changed after the chart was drawn is
    ## picked up.
    ##
    ## @end deftp
    SourceTable = [];

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} XVariable
    ##
    ## The column of @code{SourceTable} the boxes are grouped by, empty where
    ## none was named.  Naming one is refused while the data came from
    ## vectors.
    ##
    ## @end deftp
    XVariable = [];

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} YVariable
    ##
    ## The column of @code{SourceTable} holding the observations.  Naming one
    ## is refused while the data came from vectors.
    ##
    ## @end deftp
    YVariable = [];

  endproperties

  properties (GetAccess = public, SetAccess = private)

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} XDataMode
    ##
    ## Where @code{XData} came from: @qcode{'manual'} where it was given as a
    ## grouping vector, and @qcode{'auto'} where the chart derived it, from a
    ## table or from there being no grouping at all.  This property is
    ## read-only, where MATLAB allows it to be assigned.
    ##
    ## @end deftp
    XDataMode = 'auto';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} YDataMode
    ##
    ## Where @code{YData} came from: @qcode{'manual'} where it was given as a
    ## vector and @qcode{'auto'} where it was read from a table.  This
    ## property is read-only, where MATLAB allows it to be assigned.
    ##
    ## @end deftp
    YDataMode = 'manual';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} BoxFaceColorMode
    ##
    ## Whether @code{BoxFaceColor} was chosen, @qcode{'manual'}, or left to
    ## the chart, @qcode{'auto'}.  This property is read-only.
    ##
    ## @end deftp
    BoxFaceColorMode = 'auto';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} BoxEdgeColorMode
    ##
    ## Whether @code{BoxEdgeColor} was chosen.  This property is read-only.
    ##
    ## @end deftp
    BoxEdgeColorMode = 'auto';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} BoxMedianLineColorMode
    ##
    ## Whether @code{BoxMedianLineColor} was chosen.  This property is
    ## read-only.
    ##
    ## @end deftp
    BoxMedianLineColorMode = 'auto';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} MarkerColorMode
    ##
    ## Whether @code{MarkerColor} was chosen.  This property is read-only.
    ##
    ## @end deftp
    MarkerColorMode = 'auto';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} CapWidthMode
    ##
    ## Whether @code{CapWidth} was chosen.  This property is read-only.
    ##
    ## @end deftp
    CapWidthMode = 'auto';

    ## -*- texinfo -*-
    ## @deftp {BoxChart} {property} Parent
    ##
    ## The axes the chart is drawn into.  This property is read-only.
    ##
    ## @end deftp
    Parent = [];

  endproperties

  properties (Access = private, Hidden)
    Handles_ = [];      # every primitive the last drawing made
    Drawn_ = false;     # false while the constructor is still filling fields
  endproperties

  methods (Hidden)

    function disp (this)
      printf ('\n  BoxChart with properties:\n\n');
      printf ('%22s: %d\n', 'observations', numel (this.YData));
      printf ('%22s: %d\n', 'boxes', numel (bcPositions (this)));
      printf ('%22s: %s\n', 'Orientation', this.Orientation);
      printf ('%22s: %s\n', 'Notch', this.Notch);
      if (! isempty (this.YVariable))
        printf ('%22s: %s\n', 'YVariable', this.YVariable);
      endif
      printf ('\n');
    endfunction

    function display (this)
      disp (this);
    endfunction

  endmethods

  methods (Hidden)

    ## The data comes from vectors or from a table, never from both, and
    ## R2026a refuses a mixture in both directions.  Measured 2026-09-17.
    ## The rules below bite only once the chart is drawn: the constructor
    ## fills both sides through these same setters and must not trip them.
    function set.YData (this, val)
      if (this.Drawn_ && ! isempty (this.YVariable))
        error (strcat ("stats.chart.BoxChart: setting 'YData' while", ...
                       " 'YVariable' names a column is not supported."));
      endif
      this.YData = bcCheckData (val, 'YData');
      this.YDataMode = 'manual';
      redraw (this);
    endfunction

    function set.XData (this, val)
      if (this.Drawn_ && ! isempty (this.XVariable))
        error (strcat ("stats.chart.BoxChart: setting 'XData' while", ...
                       " 'XVariable' names a column is not supported."));
      endif
      this.XData = val;
      if (! isempty (val))
        this.XDataMode = 'manual';
      endif
      redraw (this);
    endfunction

    function set.YVariable (this, val)
      ## Clearing a name is always allowed; it is naming a column on data
      ## that came from vectors that R2026a refuses
      if (isempty (val))
        this.YVariable = [];
        return;
      endif
      if (this.Drawn_ && strcmp (this.YDataMode, 'manual')
          && ! isempty (this.YData))
        error (strcat ("stats.chart.BoxChart: setting 'YVariable' while", ...
                       " 'YDataMode' is 'manual' is not supported."));
      endif
      this.YVariable = bcCheckName (val, 'YVariable');
      reread (this);
    endfunction

    function set.XVariable (this, val)
      if (isempty (val))
        this.XVariable = [];
        return;
      endif
      if (this.Drawn_ && strcmp (this.XDataMode, 'manual')
          && ! isempty (this.XData))
        error (strcat ("stats.chart.BoxChart: setting 'XVariable' while", ...
                       " 'XDataMode' is 'manual' is not supported."));
      endif
      this.XVariable = bcCheckName (val, 'XVariable');
      reread (this);
    endfunction

    ## Assigning the table again re-reads the columns, so a table changed
    ## after the chart was drawn is picked up, as R2026a picks it up
    function set.SourceTable (this, val)
      if (! (isempty (val) || istable (val)))
        error ("stats.chart.BoxChart: 'SourceTable' must be a table.");
      endif
      this.SourceTable = val;
      reread (this);
    endfunction

    function set.BoxWidth (this, val)
      this.BoxWidth = bcCheckPositive (val, 'BoxWidth');
      redraw (this);
    endfunction

    function set.CapWidth (this, val)
      this.CapWidth = bcCheckPositive (val, 'CapWidth');
      this.CapWidthMode = 'manual';
      redraw (this);
    endfunction

    function set.LineWidth (this, val)
      this.LineWidth = bcCheckPositive (val, 'LineWidth');
      redraw (this);
    endfunction

    function set.MarkerSize (this, val)
      this.MarkerSize = bcCheckPositive (val, 'MarkerSize');
      redraw (this);
    endfunction

    function set.BoxFaceAlpha (this, val)
      this.BoxFaceAlpha = bcCheckUnit (val, 'BoxFaceAlpha');
      redraw (this);
    endfunction

    function set.BoxFaceColor (this, val)
      this.BoxFaceColor = bcCheckColor (val, 'BoxFaceColor');
      this.BoxFaceColorMode = 'manual';
      redraw (this);
    endfunction

    function set.BoxEdgeColor (this, val)
      this.BoxEdgeColor = bcCheckColor (val, 'BoxEdgeColor');
      this.BoxEdgeColorMode = 'manual';
      redraw (this);
    endfunction

    function set.BoxMedianLineColor (this, val)
      this.BoxMedianLineColor = bcCheckColor (val, 'BoxMedianLineColor');
      this.BoxMedianLineColorMode = 'manual';
      redraw (this);
    endfunction

    function set.MarkerColor (this, val)
      this.MarkerColor = bcCheckColor (val, 'MarkerColor');
      this.MarkerColorMode = 'manual';
      redraw (this);
    endfunction

    function set.WhiskerLineColor (this, val)
      this.WhiskerLineColor = bcCheckColor (val, 'WhiskerLineColor');
      redraw (this);
    endfunction

    function set.WhiskerLineStyle (this, val)
      this.WhiskerLineStyle = bcCheckOneOf (val, ...
                                {'-', '--', ':', '-.', 'none'}, ...
                                'WhiskerLineStyle');
      redraw (this);
    endfunction

    function set.MarkerStyle (this, val)
      this.MarkerStyle = bcCheckOneOf (val, ...
                           {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', ...
                            '>', '<', 'p', 'h', 'none'}, 'MarkerStyle');
      redraw (this);
    endfunction

    function set.Notch (this, val)
      this.Notch = bcCheckOneOf (val, {'on', 'off'}, 'Notch');
      redraw (this);
    endfunction

    function set.JitterOutliers (this, val)
      this.JitterOutliers = bcCheckOneOf (val, {'on', 'off'}, ...
                                          'JitterOutliers');
      redraw (this);
    endfunction

    function set.Orientation (this, val)
      this.Orientation = bcCheckOneOf (val, {'vertical', 'horizontal'}, ...
                                       'Orientation');
      redraw (this);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {stats.chart.BoxChart} {@var{obj} =} stats.chart.BoxChart (@var{hax}, @var{spec}, @var{args})
    ##
    ## Create a @code{stats.chart.BoxChart} object.
    ##
    ## @var{hax} is the axes to draw into, or empty for the current axes,
    ## which is resolved only once every value has been accepted.  @var{spec}
    ## is a structure carrying
    ## the data as @code{boxchart} resolved it, with the fields
    ## @qcode{XData}, @qcode{YData}, @qcode{SourceTable}, @qcode{XVariable},
    ## @qcode{YVariable} and @qcode{XDataMode}, and @var{args} the
    ## name-value pairs left to set.  The documented way to reach this
    ## constructor is @code{boxchart}.
    ##
    ## @seealso{boxchart}
    ## @end deftypefn
    function this = BoxChart (hax, spec, args)

      if (nargin < 2)
        error ("stats.chart.BoxChart: too few input arguments.");
      endif
      if (nargin < 3)
        args = {};
      endif

      this.XData = spec.XData;
      this.YData = spec.YData;
      this.SourceTable = spec.SourceTable;
      this.XVariable = spec.XVariable;
      this.YVariable = spec.YVariable;
      this.XDataMode = spec.XDataMode;
      this.YDataMode = spec.YDataMode;

      ## Whatever was named is set now, so that one drawing covers them all
      for k = 1:2:numel (args)
        name = args{k};
        value = args{k+1};
        this.(name) = value;
      endfor

      ## The axes comes last, so a rejected value leaves no figure behind
      if (isempty (hax))
        hax = gca ();
      endif
      this.Parent = hax;

      this.Drawn_ = true;
      redraw (this);

    endfunction

  endmethods

  methods (Access = private)

    ## Draw the chart from scratch, discarding whatever was drawn before.
    function redraw (this)

      if (! this.Drawn_)
        return;
      endif
      h = this.Handles_;
      for k = 1:numel (h)
        if (ishghandle (h(k)))
          delete (h(k));
        endif
      endfor
      this.Handles_ = [];
      if (isempty (this.YData))
        return;
      endif

      ax = this.Parent;
      held = ishold (ax);
      hold (ax, 'on');
      unwind_protect
        pos = bcPositions (this);
        for k = 1:numel (pos)
          drawOneBox (this, ax, pos(k), k);
        endfor
      unwind_protect_cleanup
        if (! held)
          hold (ax, 'off');
        endif
      end_unwind_protect

    endfunction

    ## Draw the box, whiskers, caps, median and outliers of one group.
    function drawOneBox (this, ax, at, k)

      y = bcGroupData (this, k);
      y = y(! isnan (y));
      if (isempty (y))
        return;
      endif

      q1 = quantile (y, 0.25);
      q3 = quantile (y, 0.75);
      md = median (y);
      sp = q3 - q1;
      out = isoutlier (y, 'quartiles');
      inside = y(! out);
      if (isempty (inside))
        lo = md;
        hi = md;
      else
        lo = min (inside);
        hi = max (inside);
      endif

      hw = this.BoxWidth / 2;
      cw = this.CapWidth / 2;

      ## The box, notched about the median where asked for
      if (strcmp (this.Notch, 'on'))
        n = numel (y);
        d = 1.57 * sp / sqrt (n);
        nlo = max (md - d, q1);
        nhi = min (md + d, q3);
        bx = [at-hw, at+hw, at+hw, at+hw*0.5, at+hw, at+hw, ...
              at-hw, at-hw, at-hw*0.5, at-hw];
        by = [q1, q1, nlo, md, nhi, q3, q3, nhi, md, nlo];
      else
        bx = [at-hw, at+hw, at+hw, at-hw];
        by = [q1, q1, q3, q3];
      endif
      this.add (bcPatch (this, ax, bx, by));

      ## The median, the whiskers and their caps
      this.add (bcLine (this, ax, [at-hw, at+hw], [md, md], ...
                        this.BoxMedianLineColor, '-'));
      this.add (bcLine (this, ax, [at, at], [q3, hi], ...
                        this.WhiskerLineColor, this.WhiskerLineStyle));
      this.add (bcLine (this, ax, [at, at], [q1, lo], ...
                        this.WhiskerLineColor, this.WhiskerLineStyle));
      this.add (bcLine (this, ax, [at-cw, at+cw], [hi, hi], ...
                        this.WhiskerLineColor, '-'));
      this.add (bcLine (this, ax, [at-cw, at+cw], [lo, lo], ...
                        this.WhiskerLineColor, '-'));

      ## Whatever the whiskers do not reach
      if (any (out) && ! strcmp (this.MarkerStyle, 'none'))
        ox = at * ones (sum (out), 1);
        if (strcmp (this.JitterOutliers, 'on'))
          ox += (rand (sum (out), 1) - 0.5) * this.BoxWidth;
        endif
        this.add (bcMarks (this, ax, ox, y(out)));
      endif

    endfunction

    ## Keep a handle so the next drawing can take it away again.
    function add (this, h)
      this.Handles_(end+1) = h;
    endfunction

    ## Take the data from the table again.  The assignments below each ask
    ## for a drawing of their own, so the flag holds them off until the last
    ## one is in and the chart is drawn once.
    function reread (this)

      if (isempty (this.SourceTable) || isempty (this.YVariable))
        return;
      endif
      t = this.SourceTable;
      names = t.Properties.VariableNames;
      if (! any (strcmp (names, this.YVariable)))
        error (strcat ("stats.chart.BoxChart: the table holds no", ...
                       " variable '%s'."), this.YVariable);
      endif
      drawn = this.Drawn_;
      this.Drawn_ = false;
      unwind_protect
        this.YData = bcCheckData (t.(this.YVariable), 'YVariable');
        this.YDataMode = 'auto';
        if (! isempty (this.XVariable))
          if (! any (strcmp (names, this.XVariable)))
            error (strcat ("stats.chart.BoxChart: the table holds no", ...
                           " variable '%s'."), this.XVariable);
          endif
          this.XData = t.(this.XVariable);
          this.XDataMode = 'auto';
        endif
      unwind_protect_cleanup
        this.Drawn_ = drawn;
      end_unwind_protect
      redraw (this);

    endfunction

  endmethods

endclassdef

## The observations, as a column of real numbers.
function v = bcCheckData (val, name)

  if (! (isnumeric (val) && isreal (val) && isvector (val)) && ! isempty (val))
    error ("stats.chart.BoxChart: '%s' must be a real numeric vector.", name);
  endif
  v = double (val(:));

endfunction

## The name of a table variable.
function v = bcCheckName (val, name)

  if (isempty (val))
    v = [];
    return;
  endif
  if (ischar (val) && isrow (val))
    v = val;
  elseif (iscellstr (val) && isscalar (val))
    v = val{1};
  elseif (isa (val, 'string') && isscalar (val))
    v = char (val);
  else
    error ("stats.chart.BoxChart: '%s' must name one table variable.", name);
  endif

endfunction

## A width, a size or any other positive scalar.
function v = bcCheckPositive (val, name)

  if (! (isnumeric (val) && isscalar (val) && isreal (val) && val > 0))
    error ("stats.chart.BoxChart: '%s' must be a positive scalar.", name);
  endif
  v = double (val);

endfunction

## A share, from none to all of it.
function v = bcCheckUnit (val, name)

  if (! (isnumeric (val) && isscalar (val) && isreal (val)
         && val >= 0 && val <= 1))
    error (strcat ("stats.chart.BoxChart: '%s' must be a scalar between", ...
                   " 0 and 1."), name);
  endif
  v = double (val);

endfunction

## A colour, as a triplet or as one of the names every toolkit knows.
function v = bcCheckColor (val, name)

  if (isnumeric (val) && isreal (val) && isequal (size (val), [1, 3])
      && all (val >= 0) && all (val <= 1))
    v = double (val);
    return;
  endif
  if (ischar (val) && isrow (val))
    known = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'};
    rgb = [1, 0, 0; 0, 1, 0; 0, 0, 1; 0, 1, 1; ...
           1, 0, 1; 1, 1, 0; 0, 0, 0; 1, 1, 1];
    long = {'red', 'green', 'blue', 'cyan', 'magenta', 'yellow', ...
            'black', 'white'};
    j = find (strcmp (known, val) | strcmp (long, val), 1);
    if (! isempty (j))
      v = rgb(j,:);
      return;
    endif
  endif
  error (strcat ("stats.chart.BoxChart: '%s' must be an RGB triplet", ...
                 " or a colour name."), name);

endfunction

## One of a short list of words.
function v = bcCheckOneOf (val, allowed, name)

  if (! (ischar (val) && isrow (val) && any (strcmpi (allowed, val))))
    error ("stats.chart.BoxChart: '%s' must be one of %s.", name, ...
           strjoin (strcat ("'", allowed, "'"), ', '));
  endif
  v = lower (val);
  j = find (strcmpi (allowed, val), 1);
  v = allowed{j};

endfunction

## Where each box stands along the grouping axis, one position per group.
function pos = bcPositions (this)

  x = this.XData;
  if (isempty (x))
    pos = 1;
    return;
  endif
  if (isa (x, 'categorical'))
    pos = 1:numel (categories (x));
  else
    pos = unique (x(! isnan (x)))';
  endif
  if (isempty (pos))
    pos = 1;
  endif

endfunction

## The observations belonging to the kth box.
function y = bcGroupData (this, k)

  x = this.XData;
  y = this.YData(:);
  if (isempty (x))
    return;
  endif
  if (isa (x, 'categorical'))
    c = categories (x);
    y = y(x == c{k});
  else
    pos = bcPositions (this);
    y = y(x(:) == pos(k));
  endif

endfunction

## The box itself, filled and outlined as the properties ask.
function h = bcPatch (this, ax, bx, by)

  [bx, by] = bcOrient (this, bx, by);
  h = patch (ax, bx, by, this.BoxFaceColor, ...
             'FaceAlpha', this.BoxFaceAlpha, ...
             'EdgeColor', this.BoxEdgeColor, ...
             'LineWidth', this.LineWidth);

endfunction

## One straight line of the chart, in the colour and style asked for.
function h = bcLine (this, ax, bx, by, col, style)

  [bx, by] = bcOrient (this, bx, by);
  h = line (ax, bx, by, 'Color', col, 'LineStyle', style, ...
            'LineWidth', this.LineWidth);

endfunction

## The markers standing for whatever the whiskers do not reach.
function h = bcMarks (this, ax, bx, by)

  [bx, by] = bcOrient (this, bx, by);
  h = line (ax, bx, by, 'LineStyle', 'none', 'Marker', this.MarkerStyle, ...
            'MarkerSize', this.MarkerSize, ...
            'MarkerEdgeColor', this.MarkerColor, ...
            'MarkerFaceColor', 'none');

endfunction

## A horizontal chart is the same drawing with the axes exchanged, which is
## done in one place so that nothing above has to know which way it lies.
function [ox, oy] = bcOrient (this, bx, by)

  if (strcmp (this.Orientation, 'horizontal'))
    ox = by;
    oy = bx;
  else
    ox = bx;
    oy = by;
  endif

endfunction

%!shared bxY, bxG, bxS, bxR
%! bxY = [1; 2; 3; 4; 5; 6; 7; 8; 9; 100];
%! bxG = [2; 2; 2; 2; 2; 5; 5; 5; 5; 5];
%! bxS = struct ('XData', [], 'YData', bxY, 'SourceTable', [], ...
%!               'XVariable', [], 'YVariable', [], ...
%!               'XDataMode', 'auto', 'YDataMode', 'manual');
%! bxR = struct ('XData', [], 'YData', [], ...
%!               'SourceTable', table (bxG, bxY, ...
%!                                     'VariableNames', {'grp', 'val'}), ...
%!               'XVariable', 'grp', 'YVariable', 'val', ...
%!               'XDataMode', 'auto', 'YDataMode', 'auto');

%!test  # the class is built from an axes and a data specification
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS);
%!   assert_equal (class (b), 'stats.chart.BoxChart');
%!   assert_equal (b.Parent, gca ());
%!   assert_equal (b.YData, bxY);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # an empty axes is the current one, and the pairs are set before drawing
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart ([], bxS, {'Notch', 'on', 'BoxWidth', 0.8});
%!   assert_equal (b.Parent, gca ());
%!   assert_equal (b.Notch, 'on');
%!   assert_equal (b.BoxWidth, 0.8);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a rejected value leaves no figure, the axes being resolved last
%! n = numel (findall (0, 'type', 'figure'));
%! msg = '';
%! try
%!   stats.chart.BoxChart ([], bxS, {'BoxWidth', -1});
%! catch err
%!   msg = err.message;
%! end_try_catch
%! assert_equal (isempty (msg), false);
%! assert_equal (numel (findall (0, 'type', 'figure')), n);

%!test  # the colours every part is drawn in default to the same triplet
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS);
%!   assert_equal (b.BoxFaceColor, [0, 0.447, 0.741]);
%!   assert_equal (b.BoxEdgeColor, [0, 0.447, 0.741]);
%!   assert_equal (b.BoxMedianLineColor, [0, 0.447, 0.741]);
%!   assert_equal (b.WhiskerLineColor, [0, 0.447, 0.741]);
%!   assert_equal (b.MarkerColor, [0, 0.447, 0.741]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # nothing was chosen, so every mode reads 'auto'
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS);
%!   assert_equal (b.BoxFaceColorMode, 'auto');
%!   assert_equal (b.BoxEdgeColorMode, 'auto');
%!   assert_equal (b.BoxMedianLineColorMode, 'auto');
%!   assert_equal (b.MarkerColorMode, 'auto');
%!   assert_equal (b.CapWidthMode, 'auto');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # setting one of them records that it was chosen
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS);
%!   b.BoxEdgeColor = [1, 0, 0];
%!   assert_equal (b.BoxEdgeColorMode, 'manual');
%!   b.BoxMedianLineColor = [1, 0, 0];
%!   assert_equal (b.BoxMedianLineColorMode, 'manual');
%!   b.MarkerColor = [1, 0, 0];
%!   assert_equal (b.MarkerColorMode, 'manual');
%!   b.CapWidth = 0.4;
%!   assert_equal (b.CapWidthMode, 'manual');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a colour name is kept as the triplet it stands for
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS);
%!   b.BoxFaceColor = 'r';
%!   assert_equal (b.BoxFaceColor, [1, 0, 0]);
%!   b.WhiskerLineColor = 'black';
%!   assert_equal (b.WhiskerLineColor, [0, 0, 0]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a word is taken whatever its case and kept in the documented spelling
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS);
%!   b.Notch = 'ON';
%!   assert_equal (b.Notch, 'on');
%!   b.Orientation = 'HORIZONTAL';
%!   assert_equal (b.Orientation, 'horizontal');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## The modes and the parent are read-only here, where MATLAB assigns them
%!test
%! m = meta.class.fromName ('stats.chart.BoxChart');
%! nm = cellfun (@(p) p.Name, m.PropertyList, 'UniformOutput', false);
%! i = find (strcmp (nm, 'XDataMode'), 1);
%! assert_equal (m.PropertyList{i}.SetAccess, 'private');
%! assert_equal (m.PropertyList{i}.GetAccess, 'public');
%! i = find (strcmp (nm, 'Parent'), 1);
%! assert_equal (m.PropertyList{i}.SetAccess, 'private');

%!test  # the box spans the quartiles and the median lies across it
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   stats.chart.BoxChart (gca (), bxS);
%!   k = flipud (get (gca (), 'children'));
%!   assert_equal (get (k(1), 'xdata'), [0.75; 1.25; 1.25; 0.75]);
%!   assert_equal (get (k(1), 'ydata'), [3; 3; 8; 8]);
%!   assert_equal (get (k(2), 'ydata'), [5.5, 5.5]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## The whiskers reach the furthest values within one and a half IQR of the box
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   stats.chart.BoxChart (gca (), bxS);
%!   k = flipud (get (gca (), 'children'));
%!   assert_equal (get (k(3), 'ydata'), [8, 9]);
%!   assert_equal (get (k(4), 'ydata'), [3, 1]);
%!   assert_equal (get (k(5), 'xdata'), [0.875, 1.125]);
%!   assert_equal (get (k(6), 'xdata'), [0.875, 1.125]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # whatever the whiskers do not reach is drawn as a marker
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS);
%!   k = flipud (get (gca (), 'children'));
%!   assert_equal (numel (k), 7);
%!   assert_equal (get (k(7), 'ydata'), 100);
%!   assert_equal (get (k(7), 'marker'), 'o');
%!   b.MarkerStyle = 'none';
%!   assert_equal (numel (get (gca (), 'children')), 6);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a jittered marker is spread across the width of the box
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS);
%!   k = flipud (get (gca (), 'children'));
%!   assert_equal (get (k(7), 'xdata'), 1);
%!   b.JitterOutliers = 'on';
%!   k = flipud (get (gca (), 'children'));
%!   assert_equal (abs (get (k(7), 'xdata') - 1) <= 0.25, true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a horizontal chart is the same drawing with the axes exchanged
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS, {'Orientation', 'horizontal'});
%!   k = flipud (get (gca (), 'children'));
%!   assert_equal (get (k(1), 'xdata'), [3; 3; 8; 8]);
%!   assert_equal (get (k(1), 'ydata'), [0.75; 1.25; 1.25; 0.75]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## A notch spans the median plus and minus 1.57 IQR / sqrt (n)
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   stats.chart.BoxChart (gca (), bxS, {'Notch', 'on'});
%!   k = flipud (get (gca (), 'children'));
%!   by = get (k(1), 'ydata');
%!   d = 1.57 * (8 - 3) / sqrt (10);
%!   assert_equal (numel (by), 10);
%!   assert_equal (by(3), 5.5 - d, eps);
%!   assert_equal (by(5), 5.5 + d, eps);
%!   assert_equal (by(4), 5.5);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a numeric grouping puts each box at the value it was given
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = bxS;
%!   s.XData = bxG;
%!   s.XDataMode = 'manual';
%!   stats.chart.BoxChart (gca (), s);
%!   p = findobj (gca (), 'type', 'patch');
%!   x1 = get (p(1), 'xdata');
%!   x2 = get (p(2), 'xdata');
%!   assert_equal (sort ([x1(1), x2(1)]), [1.75, 4.75]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a categorical grouping draws one box per level
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = bxS;
%!   s.XData = categorical ({'a'; 'a'; 'a'; 'b'; 'b'; 'b'; 'c'; 'c'; 'c'; 'c'});
%!   s.XDataMode = 'manual';
%!   stats.chart.BoxChart (gca (), s);
%!   assert_equal (numel (findobj (gca (), 'type', 'patch')), 3);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # no observations, nothing drawn, and the data may be given afterwards
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = bxS;
%!   s.YData = [];
%!   b = stats.chart.BoxChart (gca (), s);
%!   assert_equal (isempty (get (gca (), 'children')), true);
%!   b.YData = bxY;
%!   assert_equal (numel (get (gca (), 'children')), 7);
%!   assert_equal (b.YDataMode, 'manual');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a row of observations is kept as a column
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = bxS;
%!   s.YData = bxY';
%!   b = stats.chart.BoxChart (gca (), s);
%!   assert_equal (b.YData, bxY);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the chart is a handle, so two names are one chart
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS);
%!   c = b;
%!   c.BoxWidth = 0.9;
%!   assert_equal (b.BoxWidth, 0.9);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a table fills both sides and puts both modes at 'auto'
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxR);
%!   assert_equal (b.YData, bxY);
%!   assert_equal (b.XDataMode, 'auto');
%!   assert_equal (b.YDataMode, 'auto');
%!   assert_equal (numel (findobj (gca (), 'type', 'patch')), 2);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # clearing a variable name is allowed where naming one would not be
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxR);
%!   b.XVariable = [];
%!   assert_equal (isempty (b.XVariable), true);
%!   assert_equal (b.YData, bxY);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # disp heads with the short name and says what was drawn
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxS);
%!   s = evalc ('disp (b)');
%!   assert_equal (isempty (strfind (s, 'BoxChart with properties:')), false);
%!   assert_equal (isempty (strfind (s, 'observations: 10')), false);
%!   assert_equal (isempty (strfind (s, 'boxes: 1')), false);
%!   assert_equal (isempty (strfind (s, 'Orientation: vertical')), false);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a chart drawn from a table also names the column it was drawn from
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = stats.chart.BoxChart (gca (), bxR);
%!   s = evalc ('disp (b)');
%!   assert_equal (isempty (strfind (s, 'boxes: 2')), false);
%!   assert_equal (isempty (strfind (s, 'YVariable: val')), false);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Input validation
%!error<stats.chart.BoxChart: too few input arguments.> ...
%! stats.chart.BoxChart ()

%!error<stats.chart.BoxChart: 'YData' must be a real numeric vector.> ...
%! stats.chart.BoxChart ([], setfield (bxS, 'YData', {1, 2}))

%!error<stats.chart.BoxChart: 'SourceTable' must be a table.> ...
%! stats.chart.BoxChart ([], setfield (bxS, 'SourceTable', 5))

%!error<stats.chart.BoxChart: 'XVariable' must name one table variable.> ...
%! stats.chart.BoxChart ([], setfield (bxR, 'XVariable', 5))

%!error<stats.chart.BoxChart: the table holds no variable 'nope'.> ...
%! stats.chart.BoxChart ([], setfield (bxR, 'YVariable', 'nope'))

%!error<stats.chart.BoxChart: the table holds no variable 'nope'.> ...
%! stats.chart.BoxChart ([], setfield (bxR, 'XVariable', 'nope'))

%!error<stats.chart.BoxChart: 'BoxWidth' must be a positive scalar.> ...
%! stats.chart.BoxChart ([], bxS, {'BoxWidth', 0})

%!error<stats.chart.BoxChart: 'CapWidth' must be a positive scalar.> ...
%! stats.chart.BoxChart ([], bxS, {'CapWidth', [1, 2]})

%!error<stats.chart.BoxChart: 'LineWidth' must be a positive scalar.> ...
%! stats.chart.BoxChart ([], bxS, {'LineWidth', -1})

%!error<stats.chart.BoxChart: 'MarkerSize' must be a positive scalar.> ...
%! stats.chart.BoxChart ([], bxS, {'MarkerSize', 'big'})

%!error<stats.chart.BoxChart: 'BoxFaceAlpha' must be a scalar between 0 and 1.> ...
%! stats.chart.BoxChart ([], bxS, {'BoxFaceAlpha', 1.5})

%!error<stats.chart.BoxChart: 'BoxFaceColor' must be an RGB triplet or a colour name.> ...
%! stats.chart.BoxChart ([], bxS, {'BoxFaceColor', [2, 0, 0]})

%!error<stats.chart.BoxChart: 'WhiskerLineColor' must be an RGB triplet or a colour name.> ...
%! stats.chart.BoxChart ([], bxS, {'WhiskerLineColor', 'puce'})

%!error<stats.chart.BoxChart: 'WhiskerLineStyle' must be one of '-', '--', ':', '-.', 'none'.> ...
%! stats.chart.BoxChart ([], bxS, {'WhiskerLineStyle', 'wavy'})

## The pattern of the next block stops at 'v': 'test.m' ends it at the first
## '>' of the line, with 'index (str, ">")', so the '>' among the markers this
## message lists cannot be matched.  Measured on Octave 11.3.0.
%!error<stats.chart.BoxChart: 'MarkerStyle' must be one of 'o', '\+', '\*', '.', 'x', 's', 'd', '\^', 'v'> ...
%! stats.chart.BoxChart ([], bxS, {'MarkerStyle', 'q'})

%!error<stats.chart.BoxChart: 'Notch' must be one of 'on', 'off'.> ...
%! stats.chart.BoxChart ([], bxS, {'Notch', 'maybe'})

%!error<stats.chart.BoxChart: 'JitterOutliers' must be one of 'on', 'off'.> ...
%! stats.chart.BoxChart ([], bxS, {'JitterOutliers', 'yes'})

%!error<stats.chart.BoxChart: 'Orientation' must be one of 'vertical', 'horizontal'.> ...
%! stats.chart.BoxChart ([], bxS, {'Orientation', 'sideways'})

## The data comes from vectors or from a table, never from both
%!error<stats.chart.BoxChart: setting 'YData' while 'YVariable' names a column is not supported.>
%! hf = figure ('visible', 'off');
%! b = stats.chart.BoxChart (gca (), bxR);
%! close (hf);
%! b.YData = bxY;

%!error<stats.chart.BoxChart: setting 'XData' while 'XVariable' names a column is not supported.>
%! hf = figure ('visible', 'off');
%! b = stats.chart.BoxChart (gca (), bxR);
%! close (hf);
%! b.XData = bxG;

%!error<stats.chart.BoxChart: setting 'YVariable' while 'YDataMode' is 'manual' is not supported.>
%! hf = figure ('visible', 'off');
%! b = stats.chart.BoxChart (gca (), bxS);
%! close (hf);
%! b.YVariable = 'val';

%!error<stats.chart.BoxChart: setting 'XVariable' while 'XDataMode' is 'manual' is not supported.>
%! hf = figure ('visible', 'off');
%! s = bxS;
%! s.XData = bxG;
%! s.XDataMode = 'manual';
%! b = stats.chart.BoxChart (gca (), s);
%! close (hf);
%! b.XVariable = 'grp';
