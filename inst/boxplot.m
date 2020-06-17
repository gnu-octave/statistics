## Copyright (C) 2002 Alberto Terruzzi <t-albert@libero.it>
## Copyright (C) 2006 Alberto Pose <apose@alu.itba.edu.ar>
## Copyright (C) 2011 Pascal Dupuis <Pascal.Dupuis@worldonline.be>
## Copyright (C) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
## Copyright (C) 2016 Pascal Dupuis <cdemills@gmail.com>
## Copyright (C) 2020 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

## -*- texinfo -*-
## @deftypefn {Function File} {@var{s} =} boxplot (@var{data})
## @deftypefnx {Function File} {@var{s} =} boxplot (@var{data}, @var{group})
## @deftypefnx {Function File} {@var{s} =} boxplot (@var{data}, @var{notched}, @var{symbol}, @var{orientation}, @var{whisker}, @dots{})
## @deftypefnx {Function File} {@var{s} =} boxplot (@var{data}, @var{group}, @var{notched}, @var{symbol}, @var{orientation}, @var{whisker}, @dots{})
## @deftypefnx {Function File} {@var{s} =} boxplot (@var{data}, @var{options})
## @deftypefnx {Function File} {@var{s} =} boxplot (@var{data}, @var{group}, @var{options}, @dots{})
## @deftypefnx {Function File} {[@dots{} @var{h}]=} boxplot (@var{data}, @dots{})
## Produce a box plot.
##
## A box plot is a graphical display that simultaneously describes several
## important features of a data set, such as center, spread, departure from
## symmetry, and identification of observations that lie unusually far from
## the bulk of the data.
##
## Input arguments recognized by boxplot are:
##
## @itemize
## @item
## @var{data} is a matrix with one column for each data set, or a cell vector
## with one cell for each data set.  Each cell must contain a numerical row or
## column vector (NaN and NA are ignored) and not a nested vector of cells.
##
## @item
## @var{notched} = 1 produces a notched-box plot.  Notches represent a robust
## estimate of the uncertainty about the median.
##
## @var{notched} = 0 (default) produces a rectangular box plot.
##
## @var{notched} within the interval (0,1) produces a notch of the specified
## depth.  Notched values outside (0,1) are amusing if not exactly impractical.
##
## @item
## @var{symbol} sets the symbol for the outlier values.  The default symbol
## for points that lie outside 3 times the interquartile range is 'o';
## the default symbol for points between 1.5 and 3 times the interquartile
## range is '+'. @*
## Alternative @var{symbol} settings:
##
## @var{symbol} = '.': points between 1.5 and 3 times the IQR are marked with
## '.' and points outside 3 times IQR with 'o'.
##
## @var{symbol} = ['x','*']: points between 1.5 and 3 times the IQR are marked
## with 'x' and points outside 3 times IQR with '*'.
##
## @item
## @var{orientation} = 0 makes the boxes horizontal (default).
## @var{orientation} = 1 plots the boxes vertically.  Alternatively,
## orientation options can be passed as a string, e.g., 'vertical' or
## 'horizontal'.
##
## @item
## @var{whisker} defines the length of the whiskers as a function of the IQR
## (default = 1.5).  If @var{whisker} = 0 then @code{boxplot} displays all data
## values outside the box using the plotting symbol for points that lie
## outside 3 times the IQR.
##
## @item
## @var{group} may be passed as an optional argument only in the second
## position after @var{data}.  @var{group} contains a numerical vector defining
## separate categories, each plotted in a different box, for each set of
## @var{DATA} values that share the same @var{group} value or values.  With
## the formalism (@var{data}, @var{group}), both must be vectors of the same
## length.
##
## @item
## @var{options} are additional paired arguments passed with the formalism
## (Name, Value) that provide extra functionality as listed below.
## @var{options} can be passed at any order after the initial arguments.
##
## @multitable {Name} {Value} {description} @columnfractions .2 .2 .6
## @item 'Notch' @tab  'on' @tab Notched by 0.25 of the boxes width.
## @item @tab 'off' @tab Produces a straight box.
## @item @tab scalar @tab Proportional width of the notch.
##
## @item 'Symbol' @tab '.' @tab Defines only outliers between 1.5 and 3 IQR.
## @item @tab ['x','*'] @tab 2nd character defines outliers > 3 IQR.  Valid
## characters are: '+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '<', '>',
## 'p', and 'h'.  See @function{plot} markers for more info.
##
## @item 'Orientation' @tab 'vertical' @tab Default value, can also be defined
## with numerical 1.
## @item @tab 'horizontal' @tab Can also be defined with numerical 0.
##
## @item 'Whisker' @tab scalar @tab Multiplier of IQR (defualt is 1.5).
##
## @item 'OutlierTags' @tab 'on' @tab Plot the vector index of the outlier
## value next to its point.
## @item @tab 'off' @tab No tags are plotted (default value).
##
## @item 'Sample_IDs' @tab 'cell' @tab A cell vector with one cell for each
## data set containing a nested cell vector with each sample's ID (should be
## a string).  If this option is passed, then all outliers are tagged with
## their respective sample's ID string instead of their vector's index.
##
## @item 'BoxWidth' @tab 'proportional' @tab Create boxes with their width
## proportional to the number of samples in their respective dataset (default
## value).
## @item @tab 'fixed' @tab Make all boxes with equal width.
##
## @item 'Widths' @tab scalar @tab Scaling factor for box widths (default
## value is 0.4).
##
## @item 'BoxStyle' @tab 'outline' @tab Draw boxes as outlines (default value).
## @item @tab 'filled' @tab Fill boxes with a color (outlines are still
## plotted).
##
## @item 'Positions' @tab vector @tab Numerical vector that defines the
## position of each data set.  It must have the same length as the number of
## in a desired manner.  This vector merely defines the points along the group
## axis, which by default is [1:number of groups].
##
## @item 'Labels' @tab cell @tab A cell vector of strings containing the names
## of each group.  By default each group is labeled numerically according to
## its order in the data set
## @end multitable
## @end itemize
##
## Supplemental arguments not described above (@dots{}) are concatenated and
## passed to the plot() function.
##
## The returned matrix @var{s} has one column for each data set as follows:
##
## @multitable @columnfractions .1 .8
## @item 1 @tab Minimum
## @item 2 @tab 1st quartile
## @item 3 @tab 2nd quartile (median)
## @item 4 @tab 3rd quartile
## @item 5 @tab Maximum
## @item 6 @tab Lower confidence limit for median
## @item 7 @tab Upper confidence limit for median
## @end multitable
##
## The returned structure @var{h} contains handles to the plot elements,
## allowing customization of the visualization using set/get functions.
##
## Example
##
## @example
## title ("Grade 3 heights");
## axis ([0,3]);
## set(gca (), "xtick", [1 2], "xticklabel", @{"girls", "boys"@});
## boxplot (@{randn(10,1)*5+140, randn(13,1)*8+135@});
## @end example
##
## @end deftypefn

function [s_o, hs_o] = boxplot (data, varargin)

  ## Assign parameter defaults
  if (nargin < 1)
    print_usage;
  endif

  ## Default values
  maxwhisker = 1.5;
  orientation = 1;
  symbol = ["+", "o"];
  notched = 0;
  plot_opts = {};
  groups = [];
  sample_IDs = {};
  outlier_tags = 0;
  box_width = "proportional";
  widths = 0.4;
  box_style = 0;
  positions = [];
  labels = {};

  ## Optional arguments analysis
  numarg = nargin - 1;
  option_args = {"Notch"; "Symbol"; "Orientation"; "Whisker"; ...
                 "OutlierTags"; "Sample_IDs"; "BoxWidth"; "Widths"; ...
                 "BoxStyle"; "Positions"; "Labels"};
  indopt = 1;
  group_exists = 0;
  while (numarg)
    dummy = varargin{indopt++};
    if (! ischar (dummy) && indopt < 6)
      ## MATLAB allows passing the second argument as a grouping vector
      if (length (dummy) > 1)
        if (2 != indopt)
          error ("boxplot: grouping vector may only be passed as second arg");
        endif
        groups = dummy;
        group_exists = 1;
      elseif (length (dummy) == 1)
        ## Old way: positional argument
        switch indopt - group_exists
          case 2
            notched = dummy;
          case 4
            orientation = dummy;
          case 5
            maxwhisker = dummy;
          otherwise
            error("boxplot: no positional argument allowed at position %d", ...
                  --indopt);
        endswitch
      endif
      numarg--;
      continue;
    else
      if (3 == (indopt - group_exists) && length (dummy) <= 2)
        symbol = dummy;
        numarg--;
        continue;
      else
        ## Check for additional paired arguments
        tt = find (strcmpi (dummy, option_args));
        switch (tt)
          case 1
            notched = varargin{indopt};
            ## Check for string input: "on" or "off"
            if strcmpi (notched, "on")
              notched = 1;
            elseif strcmpi (notched, "off")
              notched = 0;
            else
              msg = ["boxplot: 'Notch' input argument takes only 'on'", ...
                    " or 'off' as input parameter"];
              error (msg);
            endif
          case 2
            symbol = varargin{indopt};
          case 3
            orientation = varargin{indopt};
            ## Check for string input: "vertical" or "horizontal"
            if strcmpi (orientation, "vertical")
              notched = 1;
            elseif strcmpi (orientation, "horizontal")
              notched = 0;
            else
              msg = ["boxplot: 'Orientation' input argument takes only", ...
                    " 'vertical' or 'horizontal' as input parameter"];
              error (msg);
            endif
          case 4
            maxwhisker = varargin{indopt};
            ## Check for scalar value
            if (! isscalar (maxwhisker))
              msg = ["boxplot: 'Whisker' input argument takes only", ...
                    " a real scalar value as input parameter"];
              error (msg);
            endif
          case 5
            outlier_tags = varargin{indopt};
            ## Check for string input: "on" or "off"
            if strcmpi (outlier_tags, "on")
              outlier_tags = 1;
            elseif strcmpi (outlier_tags, "off")
              outlier_tags = 0;
            else
              msg = ["boxplot: 'OutlierTags' input argument takes only", ...
                    " 'on' or 'off' as input parameter"];
              error (msg);
            endif
          case 6
            sample_IDs = varargin{indopt};
            ## Check for cell array
            if (! iscell (sample_IDs))
              msg = ["boxplot: 'Labels' input argument takes only", ...
                    " a cell array as input parameter"];
              error (msg);
            endif
            outlier_tags = 1;
          case 7
            box_width = varargin{indopt};
            ## Check for string input: "fixed" or "proportional" (default
            ## if mis-spelled)
            if strcmpi (box_width, "fixed")
              box_width = "fixed";
            elseif strcmpi (box_width, "proportional")
              box_width = "proportional";
            else
              msg = ["boxplot: 'BoxWidth' input argument takes only", ...
                    " 'fixed' or 'proportional' as input parameter"];
              error (msg);
            endif
          case 8
            widths = varargin{indopt};
            ## Check for scalar value
            if (! isscalar (widths))
              msg = ["boxplot: 'Widths' input argument takes only", ...
                    " a real scalar value as input parameter"];
              error (msg);
            endif
          case 9
            box_style = varargin{indopt};
            ## Check for string input: "outline" or "filled"
            if strcmpi (box_style, "outline")
              box_style = 0;
            elseif strcmpi (box_style, "filled")
              box_style = 1;
            else
              msg = ["boxplot: 'BoxStyle' input argument takes only", ...
                    " 'outline' or 'filled' as input parameter"];
              error (msg);
            endif
          case 10
            positions = varargin{indopt};
            ## Check for numeric vector
            if (! isvector (positions) || ! isnumeric (positions))
              msg = ["boxplot: 'Positions' input argument takes only", ...
                    " a numeric vector as input parameter"];
              error (msg);
            endif
          case 11
            labels = varargin{indopt};
            ## Check for cell array
            if (! iscell (labels))
              msg = ["boxplot: 'Labels' input argument takes only", ...
                    " a cell array as input parameter"];
              error (msg);
            endif
          otherwise
            ## Take two args and append them to plot_opts
            plot_opts(1, end+1:end+2) = {dummy, varargin{indopt}};
        endswitch
      endif
      numarg -= 2;
      indopt++;
    endif
  endwhile

  if (1 == length (symbol))
    symbol(2) = symbol(1);
  endif

  if (1 == notched)
    notched = 0.25;
  endif
  a = 1-notched;

  ## Figure out how many data sets we have
  if (isempty (groups))
    if (iscell (data))
      nc = length (data);
      for ind_c = (1:nc)
        lc(ind_c) = length (data{ind_c});
      endfor
    else
      if (isvector (data))
        data = data(:);
      endif
      nc = columns (data);
      lc = ones (1,nc) * rows (data);
    endif
    groups = (1:nc);
    ## Check if sample_IDs exists that it has same size with data
    if (! isempty (sample_IDs) && length (sample_IDs) == 1)
      for ind_c = (1:nc)
        if (lc(ind_c) != length (sample_IDs))
          error ("boxplot: Sample_IDs must match the data");
        endif
      endfor
    elseif (! isempty (sample_IDs) && length (sample_IDs) == nc)
      for ind_c = (1:nc)
        if (lc(ind_c) != length (sample_IDs{ind_c}))
          error ("boxplot: Sample_IDs must match the data");
        endif
      endfor
    elseif (! isempty (sample_IDs) && length (sample_IDs) != nc)
      error ("boxplot: Sample_IDs must match the data");
    endif
    ## Create labels according to number of datasets as ordered in data
    ## in case they are not provided by the user as optional argument
    if (isempty (labels))
      for i = 1:nc
        column_label = num2str (groups(i));
        labels(i) = {column_label};
      endfor
    endif
  else
    if (! isvector (data))
      error ("boxplot: with the formalism (data, group), both must be vectors");
    endif
    ## Check if sample IDs exist that they have same size with data
    if (! isempty (sample_IDs))
      if (length (sample_IDs) != 1 || length (sample_IDs{1}) != length (data))
        error ("boxplot: Sample_IDs must match the data");
      endif
      nc = unique (groups);
      dummy_data = cell (1, length (nc));
      dummy_sIDs = cell (1, length (nc));
      ## Check if groups are parsed as a numeric vector
      if isnumeric (groups)
        for ind_c = (1:length (nc))
          dummy_data(ind_c) = data(groups == nc(ind_c));
          dummy_sIDs(ind_c) = {sample_IDs{1}(groups == nc(ind_c))};
        endfor
        ## Create labels according to unique numeric groups in case
        ## they are not provided by the user as optional argument
        if (isempty (labels))
          for i = 1:nc
            column_label = num2str (groups(i));
            labels(i) = {column_label};
          endfor
        endif
      ## Check if groups are parsed as a cell string vector
      elseif iscellstr (groups)
        for ind_c = (1:length (nc))
          dummy_data(ind_c) = data(ismember (group, nc(ind_c)));
          dummy_sIDs(ind_c) = {sample_IDs{1}(ismember (group, nc(ind_c)))};
        endfor
        ## Create labels according to unique cell string groups in case
        ## they are not provided by the user as optional argument
        if (isempty (labels))
          labels = nc;
        endif
      else
        error ("boxplot: group argument must be numeric or cell string vector");
      endif
      data = dummy_data;
      sample_IDs = dummy_sIDs;
      nc = length (nc);
      groups = (1:nc);
    else
      nc = unique (groups);
      dummy_data = cell (1, length (nc));
      ## Check if groups are parsed as a numeric vector
      if isnumeric (groups)
        for ind_c = (1:length (nc))
          dummy_data(ind_c) = data(groups == nc(ind_c));
        endfor
        ## Create labels according to unique numeric groups in case
        ## they are not provided by the user as optional argument
        if (isempty (labels))
          for i = 1:nc
            column_label = num2str (groups(i));
            labels(i) = {column_label};
          endfor
        endif
      ## Check if groups are parsed as a cell string vector
      elseif iscellstr (groups)
        for ind_c = (1:length (nc))
          dummy_data(ind_c) = data(ismember (group, nc(ind_c)));
        endfor
        ## Create labels according to unique cell string groups in case
        ## they are not provided by the user as optional argument
        if (isempty (labels))
          labels = nc;
        endif
      else
        error ("boxplot: group argument must be numeric vector or cell string");
      endif
      data = dummy_data;
      nc = length (nc);
      groups = (1:nc);
    endif
  endif

  ## Compute statistics.
  ## s will contain
  ##    1,5    min and max
  ##    2,3,4  1st, 2nd and 3rd quartile
  ##    6,7    lower and upper confidence intervals for median
  s = zeros (7, nc);
  box = zeros (1, nc);
  ## Arrange the boxes into desired positions (if requested, otherwise leave
  ## default 1:nc)
  if (! isempty (positions))
    groups = positions;
  endif
  ## Initialize whisker matrices to correct size and all necessary outlier
  ## variables
  whisker_x = ones (2, 1) * [groups, groups];
  whisker_y = zeros (2, 2 * nc);
  outliers_x = [];
  outliers_y = [];
  outliers_idx = [];
  outliers_IDs = {};
  outliers2_x = [];
  outliers2_y = [];
  outliers2_idx = [];
  outliers2_IDs = {};

  for indi = (1:nc)
    ## Get the next data set from the array or cell array
    if (iscell (data))
      col = data{indi}(:);
      if (!isempty (sample_IDs))
        sIDs = sample_IDs{indi};
      else
        sIDs = num2cell([1:length(col)]);
      endif
    else
      col = data(:, indi);
      sIDs = num2cell([1:length(col)]);
    endif
    ## Skip missing data (NaN, NA) and remove respective sample IDs.
    ## Do this only on nonempty data
    if length(col) > 0
      remove_samples = find (col(isnan (col) | isna (col)));
      if length (remove_samples) > 0
        col(remove_samples) = [];
        sIDs(remove_samples) = [];
      endif
    endif
    ## Remember data length
    nd = length (col);
    box(indi) = nd;
    if (nd > 1)
      ## Min, max and quartiles
      s(1:5, indi) = statistics (col)(1:5);
      ## Confidence interval for the median
      est = 1.57 * (s(4, indi) - s(2, indi)) / sqrt (nd);
      s(6, indi) = max ([s(3, indi) - est, s(2, indi)]);
      s(7, indi) = min ([s(3, indi) + est, s(4, indi)]);
      ## Whiskers out to the last point within the desired inter-quartile range
      IQR = maxwhisker * (s(4, indi) - s(2, indi));
      whisker_y(:, indi) = [min(col(col >= s(2, indi) - IQR)); s(2, indi)];
      whisker_y(:, nc+indi) = [max(col(col <= s(4, indi) + IQR)); s(4, indi)];
      ## Outliers beyond 1 and 2 inter-quartile ranges
      outliers = col((col < s(2, indi) - IQR & col >= s(2, indi) - 2 * IQR) | ...
                     (col > s(4, indi) + IQR & col <= s(4, indi) + 2 * IQR));
      outliers2 = col(col < s(2, indi) - 2 * IQR | col > s(4, indi) + 2 * IQR);
      ## Get outliers indices from this dataset
      if (length (outliers) > 0)
        for out_i = 1:length (outliers)
          outliers_idx = [outliers_idx; (find (col == outliers(out_i)))];
          outliers_IDs = {outliers_IDs{:}, sIDs{(find (col == outliers(out_i)))}};
        endfor
      endif
      if (length (outliers2) > 0)
        for out_i = 1:length (outliers2)
          outliers2_idx = [outliers2_idx; find(col == outliers2(out_i))];
          outliers2_IDs = {outliers2_IDs{:}, sIDs{find(col == outliers2(out_i))}};
        endfor
      endif
      outliers_x = [outliers_x; (groups(indi) * ones (size (outliers)))];
      outliers_y = [outliers_y; outliers];
      outliers2_x = [outliers2_x; (groups(indi) * ones (size (outliers2)))];
      outliers2_y = [outliers2_y; outliers2];
    elseif (1 == nd)
      ## All statistics collapse to the value of the point
      s(:, indi) = col;
      ## Single point data sets are plotted as outliers.
      outliers_x = [outliers_x; groups(indi)];
      outliers_y = [outliers_y; col];
      ## Append the single point's index to keep the outliers' vector aligned
      outliers_idx = [outliers_idx; 1];
      outliers_IDs = {outliers_IDs{:}, sIDs{:}};
    else
      ## No statistics if no points
      s(:, indi) = NaN;
    endif
  endfor

  ## Note which boxes don't have enough stats
  chop = find (box <= 1);

  ## Replicate widths (if scalar or shorter vector) to match the number of boxes
  widths = widths(repmat (1:length (widths), 1, nc));
  ## Truncate just in case :)
  widths([nc+1:end]) = [];
  ## Draw a box around the quartiles, with box width being fixed or proportional
  ## to the number of items in the box.
  if (strcmpi (box_width, "proportional"))
    box = box .* (widths ./ max (box));
  else
    box = box .* (widths ./ box);
  endif
  ## Draw notches if desired.
  quartile_x = ones (11, 1) * groups + ...
               [-a; -1; -1; 1 ; 1; a; 1; 1; -1; -1; -a] * box;
  quartile_y = s([3, 7, 4, 4, 7, 3, 6, 2, 2, 6, 3], :);

  ## Draw a line through the median
  median_x = ones (2, 1) * groups + [-a; +a] * box;
  median_y = s([3, 3], :);

  ## Chop all boxes which don't have enough stats
  quartile_x(:, chop) = [];
  quartile_y(:, chop) = [];
  whisker_x(:, [chop, chop + nc]) = [];
  whisker_y(:, [chop, chop + nc]) = [];
  median_x(:, chop) = [];
  median_y(:, chop) = [];
  box(chop) = [];

  ## Add caps to the remaining whiskers
  cap_x = whisker_x;
  if (strcmpi (box_width, "proportional"))
    cap_x(1, :) -= repmat (((box .* (widths ./ max (box))) / 8), 1, 2);
    cap_x(2, :) += repmat (((box .* (widths ./ max (box))) / 8), 1, 2);
  else
    cap_x(1, :) -= repmat ((widths / 8), 1, 2);
    cap_x(2, :) += repmat ((widths / 8), 1, 2);
  endif
  cap_y = whisker_y([1, 1], :);

  ## Calculate coordinates for outlier tags
  outliers_tags_x = outliers_x + 0.08;
  outliers_tags_y = outliers_y;
  outliers2_tags_x = outliers2_x + 0.08;
  outliers2_tags_y = outliers2_y;

  ## Do the plot
  if (orientation)
    ## Define outlier_tags' vertical alignment
    outlier_tags_alignment = {"horizontalalignment", "left"};
    if (isempty (plot_opts))
      if (box_style)
        f = fill(quartile_x, quartile_y,"y");
        hold on;
      endif
      h = plot (quartile_x, quartile_y, "b;;",
                whisker_x, whisker_y, "b;;",
                cap_x, cap_y, "b;;",
                median_x, median_y, "r;;",
                outliers_x, outliers_y, [symbol(1), "r;;"],
                outliers2_x, outliers2_y, [symbol(2), "r;;"]);
      ## Print outlier tags
      if (outlier_tags == 1 && outliers_x > 0)
        t1 = plot_tags (outliers_tags_x, outliers_tags_y, outliers_idx,
                        outliers_IDs, sample_IDs, outlier_tags_alignment);
      endif
      if (outlier_tags == 1 && outliers2_x > 0)
        t2 = plot_tags (outliers2_tags_x, outliers2_tags_y, outliers2_idx,
                        outliers2_IDs, sample_IDs, outlier_tags_alignment);
      endif
    else
      if (box_style)
        f = fill(quartile_x, quartile_y,"y");
        hold on;
      endif
      h = plot (quartile_x, quartile_y, "b;;",
                whisker_x, whisker_y, "b;;",
                cap_x, cap_y, "b;;",
                median_x, median_y, "r;;",
                outliers_x, outliers_y, [symbol(1), "r;;"],
                outliers2_x, outliers2_y, [symbol(2), "r;;"], plot_opts{:});
      ## Print outlier tags
      if (outlier_tags == 1 && outliers_x > 0)
        t1 = plot_tags (outliers_tags_x, outliers_tags_y, outliers_idx,
                        outliers_IDs, sample_IDs, outlier_tags_alignment);
      endif
      if (outlier_tags == 1 && outliers2_x > 0)
        t2 = plot_tags (outliers2_tags_x, outliers2_tags_y, outliers2_idx,
                        outliers2_IDs, sample_IDs, outlier_tags_alignment);
      endif
    endif
  else
    ## Define outlier_tags' horizontal alignment
    outlier_tags_alignment = {"horizontalalignment", "left", "rotation", 90};
    if (isempty (plot_opts))
      h = plot (quartile_y, quartile_x, "b;;",
               whisker_y, whisker_x, "b;;",
               cap_y, cap_x, "b;;",
               median_y, median_x, "r;;",
               outliers_y, outliers_x, [symbol(1), "r;;"],
               outliers2_y, outliers2_x, [symbol(2), "r;;"]);
      if (box_style)
        f = fill (quartile_y, quartile_x, "y");
        hold on;
      endif
      h = plot (quartile_y, quartile_x, "b;;",
                whisker_y, whisker_x, "b;;",
                cap_y, cap_x, "b;;",
                median_y, median_x, "r;;",
                outliers_y, outliers_x, [symbol(1), "r;;"],
                outliers2_y, outliers2_x, [symbol(2), "r;;"]);
      ## Print outlier tags
      if (outlier_tags == 1 && outliers_x > 0)
        t1 = plot_tags (outliers_tags_y, outliers_tags_x, outliers_idx,
                        outliers_IDs, sample_IDs, outlier_tags_alignment);
      endif
      if (outlier_tags == 1 && outliers2_x > 0)
        t2 = plot_tags (outliers2_tags_y, outliers2_tags_x, outliers2_idx,
                        outliers2_IDs, sample_IDs, outlier_tags_alignment);
      endif
    else
      if (box_style)
        f = fill (quartile_y, quartile_x,"y");
        hold on;
      endif
      h = plot (quartile_y, quartile_x, "b;;",
                whisker_y, whisker_x, "b;;",
                cap_y, cap_x, "b;;",
                median_y, median_x, "r;;",
                outliers_y, outliers_x, [symbol(1), "r;;"],
                outliers2_y, outliers2_x, [symbol(2), "r;;"], plot_opts{:});
      ## Print outlier tags
      if (outlier_tags == 1 && outliers_x > 0)
        t1 = plot_tags (outliers_tags_y, outliers_tags_x, outliers_idx,
                        outliers_IDs, sample_IDs, outlier_tags_alignment);
      endif
      if (outlier_tags == 1 && outliers2_x > 0)
        t2 = plot_tags (outliers2_tags_y, outliers2_tags_x, outliers2_idx,
                        outliers2_IDs, sample_IDs, outlier_tags_alignment);
      endif
    endif
  endif

  ## Distribute handles for box outlines and box fill (if any)
  nq = 1 : size(quartile_x, 2);
  hs.box = h(nq);
  if (box_style)
    nf = 1 : length(groups);
    hs.box_fill = f(nf);
  else
    hs.box_fill = [];
  endif

  ## Distribute handles for whiskers (including caps) and median lines
  nw = nq(end) + [1 : 2 * (size (whisker_x, 2))];
  hs.whisker = h(nw);
  nm = nw(end)+ [1 : (size (median_x, 2))];
  hs.median = h(nm);
  ## Distribute handles for outliers (if any) and their respective tags
  ## (if applicable)
  no = nm;
  if (! isempty (outliers_y))
    no = nm(end) + [1 : size(outliers_y, 2)];
    hs.outliers = h(no);
    if (outlier_tags == 1)
      nt = 1 : length (outliers_tags_y);
      hs.out_tags = t1(nt);
    else
      hs.out_tags = [];
    endif
  else
    hs.outliers = [];
    hs.out_tags = [];
  endif
  ## Distribute handles for extreme outliers (if any) and their respective tags
  ## (if applicable)
  if (! isempty (outliers2_y))
    no2 = no(end) + [1 : size(outliers2_y, 2)];
    hs.outliers2 = h(no2);
    if (outlier_tags == 1)
      nt2 = 1 : length (outliers2_tags_y);
      hs.out_tags2 = t2(nt2);
    else
      hs.out_tags2 = [];
    endif
  else
    hs.outliers2 = [];
    hs.out_tags2 = [];
  end

  ## Redraw the median lines to avoid colour overlapping in case of 'filled'
  ## BoxStyle
  if (box_style)
    set (hs.median, "color", "r");
  endif

  ## Print labels according to orientation and return handle
  if (orientation)
    set (gca(), "xtick", groups, "xticklabel", labels);
    hs.labels = get (gcf, "currentaxes");
  else
    set (gca(), "ytick", groups, "yticklabel", labels);
    hs.labels = get (gcf, "currentaxes");
  endif
  hold off;

  ## Return output arguments if desired
  if (nargout >= 1)
    s_o = s;
  endif
  if (nargout == 2)
    hs_o = hs;
  endif

endfunction


function htags = plot_tags (out_tags_x, out_tags_y, out_idx, out_IDs, ...
                 sample_IDs, opt)

  for i=1 : length (out_tags_x)
    if (! isempty (sample_IDs))
      htags(i) = text (out_tags_x(i), out_tags_y(i), out_IDs{i}, opt{:});
    else
      htags(i) = text (out_tags_x(i), out_tags_y(i), num2str (out_idx(i)), ...
                       opt{:});
    endif
  endfor

endfunction


%!demo
%! axis ([0, 3]);
%! boxplot ({(randn(10, 1) * 5 + 140), (randn (13, 1) * 8 + 135)});
%! set (gca (), "xtick", [1 2], "xticklabel", {"girls", "boys"})
%! title ("Grade 3 heights");

%!demo
%! data = [(randn (10, 1) * 5 + 140); (randn (25, 1) * 8 + 135); ...
%!         (randn (20, 1) * 6 + 165)];
%! groups = [(ones (10, 1)); (ones (25, 1) * 2); (ones (20, 1) * 3)];
%! labels = {"Team A", "Team B", "Team C"};
%! pos = [2, 1, 3];
%! boxplot (data, groups, "Notch", "on", "Labels", labels, "Positions", pos, ...
%!          "OutlierTags", "on", "BoxStyle", "filled");
%! title ("Example of Group splitting with paired vectors");
