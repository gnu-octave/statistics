## Copyright (C) 2020-2021 Stefano Guidoni <ilguido@users.sf.net>
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

##
## Author: Stefano Guidoni <ilguido@users.sf.net>

classdef ConfusionMatrixChart < handle

  ## -*- texinfo -*-
  ## @deftypefn {} {@var{p} =} ConfusionMatrixChart ()
  ## Create object @var{p}, a Confusion Matrix Chart object.
  ##
  ## @table @asis
  ## @item @qcode{"DiagonalColor"}
  ## The color of the patches on the diagonal, default is [0.0, 0.4471, 0.7412].
  ##
  ## @item @qcode{"OffDiagonalColor"}
  ## The color of the patches off the diagonal, default is [0.851, 0.3255, 0.098].
  ##
  ## @item @qcode{"GridVisible"}
  ## Available values: @qcode{on} (default), @qcode{off}.
  ##
  ## @item @qcode{"Normalization"}
  ## Available values: @qcode{absolute} (default), @qcode{column-normalized},
  ## @qcode{row-normalized}, @qcode{total-normalized}.
  ##
  ## @item @qcode{"ColumnSummary"}
  ## Available values: @qcode{off} (default), @qcode{absolute},
  ## @qcode{column-normalized},@qcode{total-normalized}.
  ##
  ## @item @qcode{"RowSummary"}
  ## Available values: @qcode{off} (default), @qcode{absolute},
  ## @qcode{row-normalized}, @qcode{total-normalized}.
  ## @end table
  ##
  ## MATLAB compatibility -- the not implemented properties are: FontColor,
  ## PositionConstraint, InnerPosition, Layout.
  ##
  ## @seealso{confusionchart}
  ## @end deftypefn

  properties (Access = public)
    ## text properties
    XLabel = "Predicted Class";
    YLabel = "True Class";
    Title  = "";

    FontName  = "";
    FontSize  = 0;

    ## chart colours
    DiagonalColor = [0 0.4471 0.7412];
    OffDiagonalColor = [0.8510 0.3255 0.0980];

    ## data visualization
    Normalization = "absolute";
    ColumnSummary = "off";
    RowSummary = "off";

    GridVisible = "on";

    HandleVisibility = "";
    OuterPosition = [];
    Position = [];
    Units = "";
  endproperties

  properties (GetAccess = public, SetAccess = private)
    ClassLabels = {}; # a string cell array of classes
    NormalizedValues = []; # the normalized confusion matrix
    Parent = 0; # a handle to the parent object
  endproperties

  properties (Access = protected)
    hax = 0.0; # a handle to the axes
    ClassN = 0; # the number of classes
    AbsoluteValues = []; # the original confusion matrix
    ColumnSummaryAbsoluteValues = []; # default values of the column summary
    RowSummaryAbsoluteValues = []; # default values of the row summary
  endproperties

  methods (Access = public)
    ## class constructor
    ## inputs: axis handle, a confusion matrix, a list of class labels,
    ##         an array of optional property-value pairs.
    function this = ConfusionMatrixChart (hax, cm, cl, args)
      ## class initialization
      this.hax = hax;
      this.Parent = get (this.hax, "parent");
      this.ClassLabels = cl;
      this.NormalizedValues = cm;
      this.AbsoluteValues = cm;
      this.ClassN = rows (cm);
      this.FontName = get (this.hax, "fontname");
      this.FontSize = get (this.hax, "fontsize");

      set (this.hax, "xlabel", this.XLabel);
      set (this.hax, "ylabel", this.YLabel);

      ## draw the chart
      draw (this);

      ## apply paired properties
      if (! isempty (args))
        pair_idx = 1;
        while (pair_idx < length (args))
          switch (args{pair_idx})
            case "XLabel"
              this.XLabel = args{pair_idx + 1};
            case "YLabel"
              this.YLabel = args{pair_idx + 1};
            case "Title"
              this.Title = args{pair_idx + 1};
            case "FontName"
              this.FontName = args{pair_idx + 1};
            case "FontSize"
              this.FontSize = args{pair_idx + 1};
            case "DiagonalColor"
              this.DiagonalColor = args{pair_idx + 1};
            case "OffDiagonalColor"
              this.OffDiagonalColor = args{pair_idx + 1};
            case "Normalization"
              this.Normalization = args{pair_idx + 1};
            case "ColumnSummary"
              this.ColumnSummary = args{pair_idx + 1};
            case "RowSummary"
              this.RowSummary = args{pair_idx + 1};
            case "GridVisible"
              this.GridVisible = args{pair_idx + 1};
            case "HandleVisibility"
              this.HandleVisibility = args{pair_idx + 1};
            case "OuterPosition"
              this.OuterPosition = args{pair_idx + 1};
            case "Position"
              this.Position = args{pair_idx + 1};
            case "Units"
              this.Units = args{pair_idx + 1};
          otherwise
              close (this.Parent);
              error ("confusionchart: invalid property %s", args{pair_idx});
          endswitch

          pair_idx += 2;
        endwhile
      endif

      ## init the color map
      updateColorMap (this);
    endfunction

    ## set functions
    function set.XLabel (this, string)
      if (! ischar (string))
        close (this.Parent);
        error ("confusionchart: XLabel must be a string.");
      endif

      this.XLabel = updateAxesProperties (this, "xlabel", string);
    endfunction

    function set.YLabel (this, string)
      if (! ischar (string))
        close (this.Parent);
        error ("confusionchart: YLabel must be a string.");
      endif

      this.YLabel = updateAxesProperties (this, "ylabel", string);
    endfunction

    function set.Title (this, string)
      if (! ischar (string))
        close (this.Parent);
        error ("confusionchart: Title must be a string.");
      endif

      this.Title = updateAxesProperties (this, "title", string);
    endfunction

    function set.FontName (this, string)
      if (! ischar (string))
        close (this.Parent);
        error ("confusionchart: FontName must be a string.");
      endif

      this.FontName = updateTextProperties (this, "fontname", string);
    endfunction

    function set.FontSize (this, value)
      if (! isnumeric (value))
        close (this.Parent);
        error ("confusionchart: FontSize must be numeric.");
      endif

      this.FontSize = updateTextProperties (this, "fontsize", value);
    endfunction

    function set.DiagonalColor (this, color)
      if (ischar (color))
        color = this.convertNamedColor (color);
      endif

      if (! (isvector (color) && length (color) == 3 ))
        close (this.Parent);
        error ("confusionchart: DiagonalColor must be a color.");
      endif

      this.DiagonalColor = color;
      updateColorMap (this);
    endfunction

    function set.OffDiagonalColor (this, color)
      if (ischar (color))
        color = this.convertNamedColor (color);
      endif

      if (! (isvector (color) && length (color) == 3))
        close (this.Parent);
        error ("confusionchart: OffDiagonalColor must be a color.");
      endif

      this.OffDiagonalColor = color;
      updateColorMap (this);
    endfunction

    function set.Normalization (this, string)
      if (! any (strcmp (string, {"absolute", "column-normalized",...
        "row-normalized", "total-normalized"})))
        close (this.Parent);
        error ("confusionchart: invalid value for Normalization.");
      endif

      this.Normalization = string;
      updateChart (this);
    endfunction

    function set.ColumnSummary (this, string)
      if (! any (strcmp (string, {"off", "absolute", "column-normalized",...
        "total-normalized"})))
        close (this.Parent);
        error ("confusionchart: invalid value for ColumnSummary.");
      endif

      this.ColumnSummary = string;
      updateChart (this);
    endfunction

    function set.RowSummary (this, string)
      if (! any (strcmp (string, {"off", "absolute", "row-normalized",...
        "total-normalized"})))
        close (this.Parent);
        error ("confusionchart: invalid value for RowSummary.");
      endif

      this.RowSummary = string;
      updateChart (this);
    endfunction

    function set.GridVisible (this, string)
      if (! any (strcmp (string, {"off", "on"})))
        close (this.Parent);
        error ("confusionchart: invalid value for GridVisible.");
      endif

      this.GridVisible = string;
      setGridVisibility (this);
    endfunction

    function set.HandleVisibility (this, string)
      if (! any (strcmp (string, {"off", "on", "callback"})))
        close (this.Parent);
        error ("confusionchart: invalid value for HandleVisibility");
      endif

      set (this.hax, "handlevisibility", string);
    endfunction

    function set.OuterPosition (this, vector)
      if (! isvector (vector) || ! isnumeric (vector) || length (vector) != 4)
        close (this.Parent);
        error ("confusionchart: invalid value for OuterPosition");
      endif

      set (this.hax, "outerposition", vector);
    endfunction

    function set.Position (this, vector)
      if (! isvector (vector) || ! isnumeric (vector) || length (vector) != 4)
        close (this.Parent);
        error ("confusionchart: invalid value for Position");
      endif

      set (this.hax, "position", vector);
    endfunction

    function set.Units (this, string)
      if (! any (strcmp (string, {"centimeters", "characters", "inches", ...
                                  "normalized", "pixels", "points"})))
        close (this.Parent);
        error ("confusionchart: invalid value for Units");
      endif

      set (this.hax, "units", string);
    endfunction

    ## display method
    ## MATLAB compatibility, this tries to mimic the MATLAB behaviour
    function disp (this)
      nv_sizes = size (this.NormalizedValues);
      cl_sizes = size (this.ClassLabels);

      printf ("%s with properties:\n\n", class (this));
      printf ("\tNormalizedValues: [ %dx%d %s ]\n", nv_sizes(1), nv_sizes(2),...
        class (this.NormalizedValues));
      printf ("\tClassLabels: { %dx%d %s }\n\n", cl_sizes(1), cl_sizes(2),...
        class (this.ClassLabels));
    endfunction

    ## sortClasses
    ## reorder the chart
    function sortClasses (this, order)
      ## -*- texinfo -*-
      ## @deftypefn  {} {} sortClasses (@var{cm},@var{order})
      ## Sort the classes of the @code{ConfusionMatriChart} object @var{cm}
      ## according to @var{order}.
      ##
      ## Valid values for @var{order} can be an array or cell array including
      ## the same class labels as @var{cm}, or a value like @code{'auto'},
      ## @code{'ascending-diagonal'}, @code{'descending-diagonal'} and
      ## @code{'cluster'}.
      ##
      ## @end deftypefn
      ##
      ## @seealso{confusionchart, linkage, pdist}

      ## check the input parameters
      if (nargin != 2)
        print_usage ();
      endif

      cl = this.ClassLabels;
      cm_size = this.ClassN;
      nv = this.NormalizedValues;
      av = this.AbsoluteValues;
      cv = this.ColumnSummaryAbsoluteValues;
      rv = this.RowSummaryAbsoluteValues;

      scl = {};
      Idx = [];

      if (strcmp (order, "auto"))
        [scl, Idx] = sort (cl);
      elseif (strcmp (order, "ascending-diagonal"))
        [s, Idx] = sort (diag (nv));
        scl = cl(Idx);
      elseif (strcmp (order, "descending-diagonal"))
        [s, Idx] = sort (diag (nv));
        Idx = flip (Idx);
        scl = cl(Idx);
      elseif (strcmp (order, "cluster"))
        ## the classes are all grouped together
        ## this way one can visually evaluate which are the most similar classes
        ## according to the learning algorithm
        D = zeros (1, ((cm_size - 1) * cm_size / 2)); # a pdist like vector
        maxD = 2 * max (max (av));
        k = 1; # better than computing the index at every cycle
        for i = 1 : (cm_size - 1)
          for j = (i + 1) : cm_size
            D(k++) = maxD - (av(i, j) + av(j, i)); # distance
          endfor
        endfor
        tree = linkage (D, "average"); # clustering
        ## we could have optimal leaf ordering with
        Idx = optimalleaforder (tree, D); # optimal clustering
        ## [sorted_v Idx] = sort (cluster (tree, ));
        nodes_to_visit = 2 * cm_size - 1;
        nodecount = 0;
        while (! isempty (nodes_to_visit))
          current_node = nodes_to_visit(1);
          nodes_to_visit(1) = [];
          if (current_node > cm_size)
            node = current_node - cm_size;
            nodes_to_visit = [tree(node,[2 1]) nodes_to_visit];
          end

          if (current_node <= cm_size)
            nodecount++;
            Idx(nodecount) = current_node;
          end
        end
        ##
        scl = cl(Idx);
      else
        ## must be an array or cell array of labels
        if (! iscellstr (order))
          if (! ischar (order))
            if (isrow (order))
              order = vec (order);
            endif
            order = num2str (order);
          endif

          scl = cellstr (order);
        endif

        if (length (scl) != length (cl))
          error ("sortClasses: wrong size for order.")
        endif

        Idx = zeros (length (scl), 1);

        for i = 1 : length (scl)
          Idx(i) = find (strcmp (cl, scl{i}));
        endfor
      endif

      ## rearrange the normalized values...
      nv = nv(Idx, :);
      nv = nv(:, Idx);
      this.NormalizedValues = nv;

      ## ...and the absolute values...
      av = av(Idx, :);
      av = av(:, Idx);
      this.AbsoluteValues = av;

      cv = cv([Idx ( Idx + cm_size )]);
      this.ColumnSummaryAbsoluteValues = cv;

      rv = rv([Idx ( Idx + cm_size )]);
      this.RowSummaryAbsoluteValues = rv;

      ## ...and the class labels
      this.ClassLabels = scl;

      ## update the axes
      set (this.hax, "xtick", (0.5 : 1 : (cm_size - 0.5)), "xticklabel", scl,...
          "ytick", (0.5 : 1 : (cm_size - 0.5)), "yticklabel", scl);

      ## get text and patch handles
      kids = get (this.hax, "children");
      t_kids = kids(find (isprop (kids, "fontname"))); # hack to find texts
      m_kid = kids(find (strcmp (get (kids, "userdata"), "MainChart")));
      c_kid = kids(find (strcmp (get (kids, "userdata"), "ColumnSummary")));
      r_kid = kids(find (strcmp (get (kids, "userdata"), "RowSummary")));

      ## re-assign colors to the main chart
      cdata_m = reshape (get (m_kid, "cdata"), cm_size, cm_size);
      cdata_m = cdata_m(Idx, :);
      cdata_m = cdata_m(:, Idx);

      cdata_v = vec (cdata_m);

      set (m_kid, "cdata", cdata_v);

      ## re-assign colors to the column summary
      cdata_m = reshape (transpose (get (c_kid, "cdata")), cm_size, 2);
      cdata_m = cdata_m(Idx, :);

      cdata_v = vec (cdata_m);

      set (c_kid, "cdata", cdata_v);

      ## re-assign colors to the row summary
      cdata_m = reshape (get (r_kid, "cdata"), cm_size, 2);
      cdata_m = cdata_m(Idx, :);

      cdata_v = vec (cdata_m);

      set (r_kid, "cdata", cdata_v);

      ## move the text labels
      for i = 1:length (t_kids)
        t_pos = get (t_kids(i), "userdata");

        if (t_pos(2) > cm_size)
          ## row summary
          t_pos(1) = find (Idx == (t_pos(1) + 1)) - 1;
          set (t_kids(i), "userdata", t_pos);

          t_pos = t_pos([2 1]) + 0.5;
          set (t_kids(i), "position", t_pos);
        elseif (t_pos(1) > cm_size)
          ## column summary
          t_pos(2) = find (Idx == (t_pos(2) + 1)) - 1;
          set (t_kids(i), "userdata", t_pos);

          t_pos = t_pos([2 1]) + 0.5;
          set (t_kids(i), "position", t_pos);
        else
          ## main chart
          t_pos(1) = find (Idx == (t_pos(1) + 1)) - 1;
          t_pos(2) = find (Idx == (t_pos(2) + 1)) - 1;
          set (t_kids(i), "userdata", t_pos);

          t_pos = t_pos([2 1]) + 0.5;
          set (t_kids(i), "position", t_pos);
        endif
      endfor

      updateChart (this);
    endfunction
  endmethods

  methods (Access = private)
    ## convertNamedColor
    ## convert a named colour to a colour triplet
    function ret = convertNamedColor (this, color)
      vColorNames = ["ymcrgbwk"]';
      vColorTriplets = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 1 1 1; 0 0 0];
      if (strcmp (color, "black"))
        color = 'k';
      endif

      index = find (vColorNames == color(1));
      if (! isempty (index))
        ret = vColorTriplets(index, :);
      else
        ret = []; # trigger an error message
      endif
    endfunction

    ## updateAxesProperties
    ## update the properties of the axes
    function ret = updateAxesProperties (this, prop, value)
      set (this.hax, prop, value);

      ret = value;
    endfunction

    ## updateTextProperties
    ## set the properties of the texts
    function ret = updateTextProperties (this, prop, value)
      hax_kids = get (this.hax, "children");
      text_kids = hax_kids(isprop (hax_kids , "fontname")); # hack to find texts
      text_kids(end + 1) = get (this.hax, "xlabel");
      text_kids(end + 1) = get (this.hax, "ylabel");
      text_kids(end + 1) = get (this.hax, "title");

      updateAxesProperties (this, prop, value);
      set (text_kids, prop, value);

      ret = value;
    endfunction

    ## setGridVisibility
    ## toggle the visibility of the grid
    function setGridVisibility (this)
      kids = get (this.hax, "children");
      kids = kids(find (isprop (kids, "linestyle")));

      if (strcmp (this.GridVisible, "on"))
        set (kids, "linestyle", "-");
      else
        set (kids, "linestyle", "none");
      endif
    endfunction

    ## updateColorMap
    ## change the colormap and, accordingly, the text colors
    function updateColorMap (this)
      cm_size = this.ClassN;
      d_color = this.DiagonalColor;
      o_color = this.OffDiagonalColor;

      ## quick hack
      d_color(find (d_color == 1.0)) = 0.999;
      o_color(find (o_color == 1.0)) = 0.999;

      ## 64 shades for each color
      cm_colormap(1:64,:) = [1.0 : (-(1.0 - o_color(1)) / 63) : o_color(1);...
        1.0 : (-(1.0 - o_color(2)) / 63) : o_color(2);...
        1.0 : (-(1.0 - o_color(3)) / 63) : o_color(3)]';
      cm_colormap(65:128,:) = [1.0 : (-(1.0 - d_color(1)) / 63) : d_color(1);...
        1.0 : (-(1.0 - d_color(2)) / 63) : d_color(2);...
        1.0 : (-(1.0 - d_color(3)) / 63) : d_color(3)]';

      colormap (this.hax, cm_colormap);

      ## update text colors
      kids = get (this.hax, "children");
      t_kids = kids(find (isprop (kids, "fontname"))); # hack to find texts
      m_patch = kids(find (strcmp (get (kids, "userdata"), "MainChart")));
      c_patch = kids(find (strcmp (get (kids, "userdata"), "ColumnSummary")));
      r_patch = kids(find (strcmp (get (kids, "userdata"), "RowSummary")));

      m_colors = get (m_patch, "cdata");
      c_colors = get (c_patch, "cdata");
      r_colors = get (r_patch, "cdata");

      ## when a patch is dark, let's use a pale color for the text
      for i = 1 : length (t_kids)
        t_pos = get (t_kids(i), "userdata");
        color_idx = 1;

        if (t_pos(2) > cm_size)
          ## row summary
          idx = (t_pos(2) - cm_size - 1) * cm_size + t_pos(1) + 1;
          color_idx = r_colors(idx) + 1;
        elseif (t_pos(1) > cm_size)
          ## column summary
          idx = (t_pos(1) - cm_size - 1) * cm_size + t_pos(2) + 1;
          color_idx = c_colors(idx) + 1;
        else
          ## main chart
          idx = t_pos(2) * cm_size + t_pos(1) + 1;
          color_idx = m_colors(idx) + 1;
        endif

        if (sum (cm_colormap(color_idx, :)) < 1.8)
          set (t_kids(i), "color", [.97 .97 1.0]);
        else
          set (t_kids(i), "color", [.15 .15 .15]);
        endif
      endfor
    endfunction

    ## updateChart
    ## update the text labels and the NormalizedValues property
    function updateChart (this)
      cm_size = this.ClassN;
      cm = this.AbsoluteValues;
      l_cs = this.ColumnSummaryAbsoluteValues;
      l_rs = this.RowSummaryAbsoluteValues;

      kids = get (this.hax, "children");
      t_kids = kids(find (isprop (kids, "fontname"))); # hack to find texts

      normalization = this.Normalization;
      column_summary = this.ColumnSummary;
      row_summary = this.RowSummary;

      ## normalization for labelling
      row_totals = sum (cm, 2);
      col_totals = sum (cm, 1);
      mat_total = sum (col_totals);
      cm_labels = cm;
      add_percent = true;

      if (strcmp (normalization, "column-normalized"))
        for i = 1 : cm_size
          cm_labels(:,i) = cm_labels(:,i) ./ col_totals(i);
        endfor
      elseif (strcmp (normalization, "row-normalized"))
        for i = 1 : cm_size
          cm_labels(i,:) = cm_labels(i,:) ./ row_totals(i);
        endfor
      elseif (strcmp (normalization, "total-normalized"))
        cm_labels = cm_labels ./ mat_total;
      else
        add_percent = false;
      endif

      ## update NormalizedValues
      this.NormalizedValues = cm_labels;

      ## update axes
      last_row = cm_size;
      last_col = cm_size;
      userdata = cell2mat (get (t_kids, "userdata"));

      cs_kids = t_kids(find (userdata(:,1) > cm_size));
      cs_kids(end + 1) = kids(find (strcmp (get (kids, "userdata"),...
          "ColumnSummary")));

      if (! strcmp ("off", column_summary))
        set (cs_kids, "visible", "on");
        last_row += 3;
      else
        set (cs_kids, "visible", "off");
      endif

      rs_kids = t_kids(find (userdata(:,2) > cm_size));
      rs_kids(end + 1) = kids(find (strcmp (get (kids, "userdata"),...
          "RowSummary")));

      if (! strcmp ("off", row_summary))
        set (rs_kids, "visible", "on");
        last_col += 3;
      else
        set (rs_kids, "visible", "off");
      endif

      axis (this.hax, [0 last_col 0 last_row]);

      ## update column summary data
      cs_add_percent = true;
      if (! strcmp (column_summary, "off"))
        if (strcmp (column_summary, "column-normalized"))
          for i = 1 : cm_size
            if (col_totals(i) == 0)
              ## avoid division by zero
              l_cs([i (cm_size + i)]) = 0;
            else
              l_cs([i, cm_size + i]) = l_cs([i, cm_size + i]) ./ col_totals(i);
            endif
          endfor
        elseif strcmp (column_summary, "total-normalized")
          l_cs = l_cs ./ mat_total;
        else
          cs_add_percent = false;
        endif
      endif

      ## update row summary data
      rs_add_percent = true;
      if (! strcmp (row_summary, "off"))
        if (strcmp (row_summary, "row-normalized"))
          for i = 1 : cm_size
            if (row_totals(i) == 0)
              ## avoid division by zero
              l_rs([i (cm_size + i)]) = 0;
            else
              l_rs([i, cm_size + i]) = l_rs([i, cm_size + i]) ./ row_totals(i);
            endif
          endfor
        elseif (strcmp (row_summary, "total-normalized"))
          l_rs = l_rs ./ mat_total;
        else
          rs_add_percent = false;
        endif
      endif

      ## update text
      label_list = vec (cm_labels);

      for i = 1 : length (t_kids)
        t_pos = get (t_kids(i), "userdata");
        new_string = "";

        if (t_pos(2) > cm_size)
          ## this is the row summary
          idx = (t_pos(2) - cm_size - 1) * cm_size + t_pos(1) + 1;

          if (rs_add_percent)
            new_string = num2str (100.0 * l_rs(idx), "%3.1f");
            new_string = [new_string "%"];
          else
            new_string = num2str (l_rs(idx));
          endif
        elseif (t_pos(1) > cm_size)
          ## this is the column summary
          idx = (t_pos(1) - cm_size - 1) * cm_size + t_pos(2) + 1;

          if (cs_add_percent)
            new_string = num2str (100.0 * l_cs(idx), "%3.1f");
            new_string = [new_string "%"];
          else
            new_string = num2str (l_cs(idx));
          endif
        else
          ## this is the main chart
          idx = t_pos(2) * cm_size + t_pos(1) + 1;

          if (add_percent)
            new_string = num2str (100.0 * label_list(idx), "%3.1f");
            new_string = [new_string "%"];
          else
            new_string = num2str (label_list(idx));
          endif
        endif

        set (t_kids(i), "string", new_string);
      endfor

    endfunction

    ## draw
    ## draw the chart
    function draw (this)
      cm = this.AbsoluteValues;
      cl = this.ClassLabels;
      cm_size = this.ClassN;

      ## set up the axes
      set (this.hax, "xtick", (0.5 : 1 : (cm_size - 0.5)), "xticklabel",  cl,...
          "ytick", (0.5 : 1 : (cm_size - 0.5)), "yticklabel",  cl );
      axis ("ij");
      axis (this.hax, [0 cm_size 0 cm_size]);

      ## prepare the patches
      indices_b = 0 : (cm_size -1);
      indices_v = repmat (indices_b, cm_size, 1);
      indices_vx = transpose (vec (indices_v));
      indices_vy = vec (indices_v', 2);
      indices_ex = vec ((cm_size + 1) * [1; 2] .* ones (2, cm_size), 2);

      ## normalization for colorization
      ## it is used a colormap of 128 shades of two colors, 64 shades for each
      ## color
      normal = max (max (cm));
      cm_norm = round (63 * cm ./ normal);
      cm_norm = cm_norm + 64 * eye (cm_size);

      ## default normalization: absolute
      cm_labels = vec (cm);

      ## the patches of the main chart
      x_patch = [indices_vx;
                ( indices_vx + 1 );
                ( indices_vx + 1 );
                indices_vx];
      y_patch = [indices_vy;
                indices_vy;
                ( indices_vy + 1 );
                ( indices_vy + 1 )];
      c_patch = vec (cm_norm(1 : cm_size, 1 : cm_size));

      ## display the patches
      ph = patch (this.hax, x_patch, y_patch, c_patch);

      set (ph, "userdata", "MainChart");

      ## display the labels
      userdata = [indices_vy; indices_vx]';
      nonzero_idx = find (cm_labels != 0);
      th = text ((x_patch(1, nonzero_idx) + 0.5), (y_patch(1, nonzero_idx) +...
          0.5), num2str (cm_labels(nonzero_idx)), "parent", this.hax );

      set (th, "horizontalalignment", "center");
      for i = 1 : length (nonzero_idx)
        set (th(i), "userdata", userdata(nonzero_idx(i), :));
      endfor

      ## patches for the summaries
      main_values = diag (cm);
      ct_values = sum (cm)';
      rt_values = sum (cm, 2);
      cd_values = ct_values - main_values;
      rd_values = rt_values - main_values;

      ## column summary
      x_cs = [[indices_b indices_b];
              ( [indices_b indices_b] + 1 );
              ( [indices_b indices_b] + 1 );
              [indices_b indices_b]];
      y_cs = [(repmat ([1 1 2 2]', 1, cm_size)) (repmat ([2 2 3 3]', 1, cm_size))] +...
          cm_size;
      c_cs = [(round (63 * (main_values ./ ct_values)) + 64);
              (round (63 * (cd_values ./ ct_values)))];
      c_cs(isnan (c_cs)) = 0;
      l_cs = [main_values; cd_values];

      ph = patch (this.hax, x_cs, y_cs, c_cs);

      set (ph, "userdata", "ColumnSummary");
      set (ph, "visible", "off" );

      userdata = [y_cs(1,:); x_cs(1,:)]';
      nonzero_idx = find (l_cs != 0);
      th = text ((x_cs(1,nonzero_idx) + 0.5), (y_cs(1,nonzero_idx) + 0.5),...
          num2str (l_cs(nonzero_idx)), "parent", this.hax);

      set (th, "horizontalalignment", "center");
      for i = 1 : length (nonzero_idx)
        set (th(i), "userdata", userdata(nonzero_idx(i), :));
      endfor
      set (th, "visible", "off");

      ## row summary
      x_rs = y_cs;
      y_rs = x_cs;
      c_rs = [(round (63 * (main_values ./ rt_values)) + 64);
              (round (63 * (rd_values ./ rt_values)))];
      c_rs(isnan (c_rs)) = 0;
      l_rs = [main_values; rd_values];

      ph = patch (this.hax, x_rs, y_rs, c_rs);

      set (ph, "userdata", "RowSummary");
      set (ph, "visible", "off");

      userdata = [y_rs(1,:); x_rs(1,:)]';
      nonzero_idx = find (l_rs != 0);
      th = text ((x_rs(1,nonzero_idx) + 0.5), (y_rs(1,nonzero_idx) + 0.5),...
          num2str (l_rs(nonzero_idx)), "parent", this.hax);

      set (th, "horizontalalignment", "center");
      for i = 1 : length (nonzero_idx)
        set (th(i), "userdata", userdata(nonzero_idx(i), :));
      endfor
      set (th, "visible", "off");

      this.ColumnSummaryAbsoluteValues = l_cs;
      this.RowSummaryAbsoluteValues = l_rs;
    endfunction
  endmethods

endclassdef

