## Copyright (C) 2002 Alberto Terruzzi <t-albert@libero.it>
## Copyright (C) 2006 Alberto Pose <apose@alu.itba.edu.ar>
## Copyright (C) 2011 Pascal Dupuis <Pascal.Dupuis@worldonline.be>
## Copyright (C) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
## Copyright (C) 2016 Pascal Dupuis <cdemills@gmail.com>
## Copyright (C) 2020 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2020 Philip Nienhuis (prnienhuis@users.sf.net)
## Copyright (C) 2026 Avanish Salunke <avanishsalunke16@gmail.com>
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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{s} =} boxplot (@var{data})
## @deftypefnx {statistics} {@var{s} =} boxplot (@var{data}, @var{group})
## @deftypefnx {statistics} {@var{s} =} boxplot (@var{data}, @var{notched}, @var{symbol}, @var{orientation}, @var{whisker}, @dots{})
## @deftypefnx {statistics} {@var{s} =} boxplot (@var{data}, @var{group}, @var{notched}, @var{symbol}, @var{orientation}, @var{whisker}, @dots{})
## @deftypefnx {statistics} {@var{s} =} boxplot (@var{data}, @var{options})
## @deftypefnx {statistics} {@var{s} =} boxplot (@var{data}, @var{group}, @var{options}, @dots{})
## @deftypefnx {statistics} {[@dots{}, @var{h}] =} boxplot (@var{data}, @dots{})
##
## Produce a box plot.
##
## A box plot is a graphical display that simultaneously describes several
## important features of a data set, such as center, spread, departure from
## symmetry, and identification of observations that lie unusually far from
## the bulk of the data.
##
## Input arguments (case-insensitive) recognized by boxplot are:
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
## @var{orientation} = 0 makes the boxes horizontally. @*
## @var{orientation} = 1 plots the boxes vertically (default).  Alternatively,
## orientation can be passed as a string, e.g., 'vertical' or 'horizontal'.
##
## @item
## @var{whisker} defines the length of the whiskers as a function of the IQR
## (default = 1.5).  If @var{whisker} = 0 then @code{boxplot} displays all data
## values outside the box using the plotting symbol for points that lie
## outside 3 times the IQR.
##
## @item
## @var{group} may be passed as an optional argument only in the second
## position after @var{data}. @var{group} can be a numeric, character,
## string, or categorical vector defining separate categories. To group by
## multiple variables simultaneously, pass a cell array of grouping vectors
## (e.g., @code{@{group1, group2@}}). A separate box is plotted for each
## unique combination of group values. All grouping variables must have the
## same length as @var{data}.
##
## @item
## @var{options} are additional paired arguments passed with the formalism
## (Name, Value) that provide extra functionality as listed below.
## @var{options} can be passed at any order after the initial arguments and
## are case-insensitive.
##
## @multitable {Name} {Value} {description} @columnfractions .2 .2 .6
## @item 'Notch' @tab  'on' @tab Notched by 0.25 of the boxes width.
## @item @tab 'off' @tab Produces a straight box.
## @item @tab scalar @tab Proportional width of the notch.
##
## @item 'Symbol' @tab '.' @tab Defines only outliers between 1.5 and 3 IQR.
## @item @tab ['x','*'] @tab 2nd character defines outliers > 3 IQR
##
## @item 'Orientation' @tab 'vertical' @tab Default value, can also be defined
## with numerical 1.
## @item @tab 'horizontal' @tab Can also be defined with numerical 0.
##
## @item 'Whisker' @tab scalar @tab Multiplier of IQR (default is 1.5).
##
## @item 'OutlierTags' @tab 'on' or 1 @tab Plot the vector index of the outlier
## value next to its point.
## @item @tab 'off' or 0 @tab No tags are plotted (default value).
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
## @item 'CapWidths' @tab scalar @tab Scaling factor for whisker cap widths
## (default value is 1, which results to 'Widths'/8 halflength)
##
## @item 'BoxStyle' @tab 'outline' @tab Draw boxes as outlines (default value).
## @item @tab 'filled' @tab Fill boxes with a color (outlines are still
## plotted).
##
## @item 'Positions' @tab vector @tab Numerical vector that defines the
## position of each data set.  It must have the same length as the number of
## groups in a desired manner.  This vector merely defines the points along
## the group axis, which by default is [1:number of groups].
##
## @item 'Labels' @tab cell @tab A cell vector of strings containing the names
## of each group. By default each group is labeled numerically. If multiple
## grouping variables are provided, default labels are automatically generated
## by joining the category names and stacked hierarchically.
##
## @item 'Colors' @tab character string or Nx3 numerical matrix @tab If just
## one character or 1x3 vector of RGB values, specify the fill color of all
## boxes when BoxStyle = 'filled'.  If a character string or Nx3 matrix is
## entered, box #1's fill color corresponds to the first character or first
## matrix row, and the next boxes' fill colors corresponds to the next
## characters or rows.  If the char string or Nx3 array is exhausted the color
## selection wraps around.
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
## The quartiles are those of @code{quantile} at its default method, which is
## also what @code{prctile} returns, so the box edges of a data set always
## agree with @code{prctile (@var{data}, [25, 75])}.  They set the
## inter-quartile range, and so the whisker fences and which observations are
## reported as outliers.
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

  ## Check data
  if (! (isnumeric (data) || iscell (data)))
    error ("boxplot: numerical array or cell array containing data expected.");
  elseif (iscell (data))
    ## Check if cell contain numerical data
    if (! all (cellfun ('isnumeric', data)))
      error ("boxplot: data cells must contain numerical data.");
    endif
  endif

  ## Integer observations are perfectly good data, and MATLAB accepts them, but
  ## the quartiles are taken through statistics, whose call to var insists on
  ## floating point.  Convert rather than refuse.  Unlike MATLAB the quartiles
  ## are then computed in double, so they are not rounded back to the integer
  ## grid: the quartiles of int32 (1:7) are 2.25 and 5.75, not 3 and 6.
  if (isinteger (data))
    data = double (data);
  elseif (iscell (data) && any (cellfun (@isinteger, data)))
    data = cellfun (@double, data, 'UniformOutput', false);
  endif

  ## Default values
  maxwhisker = 1.5;
  orientation = 1;
  symbol = ['+', 'o'];
  notched = 0;
  plot_opts = {};
  groups = [];
  sample_IDs = {};
  outlier_tags = 0;
  box_width = 'proportional';
  widths = 0.4;
  capwid = 1;
  box_style = 0;
  positions = [];
  labels = {};
  nug = 0;
  bcolor = 'y';

  ## Optional arguments analysis
  numarg = nargin - 1;
  indopt = 1;
  group_exists = 0;
  while (numarg)
    dummy = varargin{indopt++};
    if ((! ischar (dummy) || iscellstr (dummy)) && indopt < 6)
      ## MATLAB allows passing the second argument as a grouping vector
      if (length (dummy) > 1)
        if (2 != indopt)
          error ("boxplot: grouping vector may only be passed as second arg.");
        endif
        if (isnumeric (dummy) || ischar (dummy) || iscell (dummy) || ...
            iscategorical (dummy) || isa (dummy, 'string'))
          groups = dummy;
          group_exists = 1;
        else
          error ("boxplot: grouping variable must be numeric, character, string, categorical, or cell array.");
        endif
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
            error ("boxplot: no positional argument allowed at position %d", ...
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
        switch lower (dummy)
          case 'notch'
            notched = varargin{indopt};
            ## Check for string input: "on" or "off"
            if (ischar (notched))
              if (strcmpi (notched, 'on'))
                notched = 1;
              elseif (strcmpi (notched, 'off'))
                notched = 0;
              else
                error (strcat ("boxplot: 'Notch' input argument accepts", ...
                               " only 'on', 'off' or a numeric scalar value."));
              endif
            elseif (! (isnumeric (notched) && isreal (notched)))
              error ("boxplot: invalid 'Notch' value.");
            endif

          case 'symbol'
            symbol = varargin{indopt};
            if (! ischar (symbol))
              error ("boxplot; Symbol(s) must be character(s)");
            endif

          case 'orientation'
            orientation = varargin{indopt};
            if (ischar (orientation))
              ## Check for string input: "vertical" or "horizontal"
              if (strcmpi (orientation, 'vertical'))
                orientation = 1;
              elseif (strcmpi (orientation, 'horizontal'))
                orientation = 0;
              else
              error (strcat ("boxplot: 'Orientation' input argument", ...
                             " accepts only 'vertical' (or 1) or", ...
                             " 'horizontal' (or 0) as value."));
              endif
            elseif (! (isnumeric (orientation) && isreal (orientation)))
              error ("boxplot: invalid 'Orientation' value.");
            endif

          case 'whisker'
            maxwhisker = varargin{indopt};
            if (! isscalar (maxwhisker) ||
                ! (isnumeric (maxwhisker) && isreal (maxwhisker)))
              error (strcat ("boxplot: 'Whisker' input argument accepts", ...
                             " only a real scalar value as input parameter."));
            endif

          case 'outliertags'
            outlier_tags = varargin{indopt};
            ## Check for string input: "on" or "off"
            if (ischar (outlier_tags))
              if (strcmpi (outlier_tags, 'on'))
                outlier_tags = 1;
              elseif (strcmpi (outlier_tags, 'off'))
                outlier_tags = 0;
              else
              error (strcat ("boxplot: 'OutlierTags' input argument accepts", ...
                             " only 'on' (or 1) or 'off' (or 0) as value."));
              endif
            elseif (! (isnumeric (outlier_tags) && isreal (outlier_tags)))
              error ("boxplot: invalid 'OutlierTags' value.");
            endif

          case 'sample_ids'
            sample_IDs = varargin{indopt};
            if (! iscell (sample_IDs))
              error (strcat ("boxplot: 'Sample_IDs' input argument", ...
                             " accepts only a cell array as value."));
            endif
            outlier_tags = 1;

          case 'boxwidth'
            box_width = varargin{indopt};
            ## Check for string input: "fixed" or "proportional"
            if (! ischar (box_width) ||
                ! ismember (lower (box_width), {'fixed', 'proportional'}))
              error (strcat ("boxplot: 'BoxWidth' input argument accepts", ...
                             " only 'fixed' or 'proportional' as value."));
            endif
            box_width = lower (box_width);

          case 'widths'
            widths = varargin{indopt};
            if (! isscalar (widths) || ! (isnumeric (widths) && isreal (widths)))
              error (strcat ("boxplot: 'Widths' input argument accepts", ...
                             " only a real scalar value as value."));
            endif

          case 'capwidths'
            capwid = varargin{indopt};
            if (! isscalar (capwid) || ! (isnumeric (capwid) && isreal (capwid)))
              error (strcat ("boxplot: 'CapWidths' input argument accepts", ...
                             " only a real scalar value as value."));
            endif

          case 'boxstyle'
            box_style = varargin{indopt};
            ## Check for string input: "outline" or "filled"
            if (! ischar (box_style) ||
                ! ismember (lower (box_style), {'outline', 'filled'}))
              error (strcat ("boxplot: 'BoxStyle' input argument accepts", ...
                             " only 'outline' or 'filled' as value."));
            endif
            box_style = lower (box_style);

          case 'positions'
            positions = varargin{indopt};
            if (! isvector (positions) || ! isnumeric (positions))
              error (strcat ("boxplot: 'Positions' input argument accepts", ...
                             " only a numeric vector as value."));
            endif

          case 'labels'
            labels = varargin{indopt};
            if (! iscellstr (labels))
              error (strcat ("boxplot: 'Labels' input argument accepts", ...
                             " only a cellstr array as value."));
            endif

          case 'colors'
            bcolor = varargin{indopt};
            if (! (ischar (bcolor) ||
                (isnumeric (bcolor) && size (bcolor, 2) == 3)))
              error (strcat ("boxplot: 'Colors' input argument accepts", ...
                             " only a character vector or Nx3 numeric", ...
                             " array as value."));
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
      nc = nug = length (data);
      for ind_c = (1:nc)
        lc(ind_c) = length (data{ind_c});
      endfor
    else
      if (isvector (data))
        data = data(:);
      endif
      nc = nug = columns (data);
      lc = ones (1, nc) * rows (data);
    endif
    groups = (1:nc);
    ## In case sample_IDs exists. check that it has same size as data
    if (! isempty (sample_IDs) && length (sample_IDs) == 1)
      for ind_c = (1:nc)
        if (lc(ind_c) != length (sample_IDs))
          error ("boxplot: Sample_IDs must match the data.");
        endif
      endfor
    elseif (! isempty (sample_IDs) && length (sample_IDs) == nc)
      for ind_c = (1:nc)
        if (lc(ind_c) != length (sample_IDs{ind_c}))
          error ("boxplot: Sample_IDs must match the data.");
        endif
      endfor
    elseif (! isempty (sample_IDs) && length (sample_IDs) != nc)
      error ("boxplot: Sample_IDs must match the data.");
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
      error ("boxplot: with the formalism (data, group), both must be vectors.");
    endif

    ## Normalize groups into a cell array of grouping variables
    if (! iscell (groups) || iscellstr (groups))
      groups_cell = {groups};
    else
      groups_cell = groups;
    endif

    ## Validate dimensions against data
    n_obs = numel (data);
    for i = 1:numel (groups_cell)
      curr_group = groups_cell{i};
      if (size (curr_group, 1) != n_obs && size (curr_group, 2) != n_obs)
        error ("boxplot: all grouping variables must have the same length as the data.");
      endif
      ## Ensure column vector format for consistent processing later
      if (isrow (curr_group) && ! ischar (curr_group))
        groups_cell{i} = curr_group(:);
      endif
    endfor

    ## If sample IDs given, check that their size matches the data
    if (! isempty (sample_IDs))
      if (length (sample_IDs) != 1 || length (sample_IDs{1}) != length (data))
        error ("boxplot: Sample_IDs must match the data");
      endif
    endif

    ## map native groups to integers
    n_groups = numel (groups_cell);
    idx_matrix = zeros (n_obs, n_groups);
    raw_uniques = cell (1, n_groups);

    for i = 1:n_groups
      [unq_vals, ~, col_idx] = unique (groups_cell{i});
      idx_matrix(:, i) = col_idx;
      raw_uniques{i} = unq_vals;
    endfor

    ## find valid group combinations.
    [valid_combinations, ~, final_numeric_groups] = unique (idx_matrix, 'rows');
    nc = size (valid_combinations, 1);
    nug = 1:nc;

    dummy_data = cell (1, nc);
    dummy_sIDs = cell (1, nc);

    ## apply mask and populate dummy arrays
    for i = 1:nc
      mask = (final_numeric_groups == i);
      dummy_data{i} = data(mask);
      if (! isempty (sample_IDs))
        dummy_sIDs{i} = sample_IDs{1}(mask);
      endif
    endfor

    ## generate labels based on native datatypes
    if (isempty (labels))
      labels = cell (1, nc);
      for i = 1:nc
        comb = valid_combinations(i, :);
        lbl_parts = cell (1, n_groups);
        for j = 1:n_groups
          if (iscell (raw_uniques{j}))
            val = raw_uniques{j}{comb(j)};
          else
            val = raw_uniques{j}(comb(j));
          endif

          if (isnumeric (val) || islogical (val))
            lbl_parts{j} = num2str (val);
          elseif (ischar (val))
            lbl_parts{j} = val;
          elseif (isa (val, 'string') || iscategorical (val))
            lbl_parts{j} = char (val);
          else
            lbl_parts{j} = '';
          endif
        endfor
        ## Join multiple group labels with a newline for hierarchical stacking
        labels{i} = strjoin (lbl_parts, char (10));
      endfor
    endif

    data = dummy_data;
    groups = nug(:).';
    if (! isempty (sample_IDs))
      sample_IDs = dummy_sIDs;
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
      if (! isempty (sample_IDs))
        sIDs = sample_IDs{indi};
      else
        sIDs = num2cell ([1:length(col)]);
      endif
    else
      col = data(:, indi);
      sIDs = num2cell ([1:length(col)]);
    endif
    ## Skip missing data (NaN, NA) and remove respective sample IDs.
    ## Do this only on nonempty data
    if (length (col) > 0)
      remove_samples = find (isnan (col) | isna (col));
      if (length (remove_samples) > 0)
        col(remove_samples) = [];
        sIDs(remove_samples) = [];
      endif
    endif
    ## Remember data length
    nd = length (col);
    box(indi) = nd;
    if (nd > 1)
      ## Min, max and quartiles.  These come from quantile at its own default
      ## method, which is the one prctile and MATLAB both use.  Taking them
      ## from the core statistics function instead put the quartiles on a
      ## different definition (it asks quantile for method 7), so a box drawn
      ## here disagreed with prctile called on the same data, and the shifted
      ## inter-quartile range moved the whisker fences with it.
      s(1:5, indi) = quantile (col, [0, 0.25, 0.5, 0.75, 1])(:);
      ## Confidence interval for the median
      est = 1.57 * (s(4, indi) - s(2, indi)) / sqrt (nd);
      s(6, indi) = max ([s(3, indi) - est, s(2, indi)]);
      s(7, indi) = min ([s(3, indi) + est, s(4, indi)]);
      ## Whiskers out to the last point within the desired inter-quartile range
      IQR = maxwhisker * (s(4, indi) - s(2, indi));
      lo_adj = min (col(col >= s(2, indi) - IQR));
      hi_adj = max (col(col <= s(4, indi) + IQR));
      ## A whisker never reaches back across its own quartile.  When no
      ## observation lies between a quartile and its fence, the whisker
      ## collapses onto the quartile rather than being drawn into the box.
      lo_adj = min (lo_adj, s(2, indi));
      hi_adj = max (hi_adj, s(4, indi));
      whisker_y(:, indi) = [lo_adj; s(2, indi)];
      whisker_y(:, nc+indi) = [hi_adj; s(4, indi)];
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
  if (strcmpi (box_width, 'proportional'))
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
  ## The cap widths below are built from BOX and WIDTHS together, so WIDTHS has
  ## to lose the same entries.  Left at its full length it made every grouped
  ## plot with a one-point group fail on a nonconformant subtraction.
  widths(chop) = [];

  ## Add caps to the remaining whiskers
  cap_x = whisker_x;
  if (strcmpi (box_width, 'proportional'))
    cap_x(1, :) -= repmat (((capwid * box .* (widths ./ max (box))) / 8), 1, 2);
    cap_x(2, :) += repmat (((capwid * box .* (widths ./ max (box))) / 8), 1, 2);
  else
    cap_x(1, :) -= repmat ((capwid * widths / 8), 1, 2);
    cap_x(2, :) += repmat ((capwid * widths / 8), 1, 2);
  endif
  cap_y = whisker_y([1, 1], :);

  ## Calculate coordinates for outlier tags
  outliers_tags_x = outliers_x + 0.08;
  outliers_tags_y = outliers_y;
  outliers2_tags_x = outliers2_x + 0.08;
  outliers2_tags_y = outliers2_y;

  ## Do the plot
  hold_status = ishold ();
  if (orientation)
    ## Define outlier_tags' vertical alignment
    outlier_tags_alignment = {'horizontalalignment', 'left'};
    if (box_style)
      f = fillbox (quartile_x, quartile_y, bcolor);
    endif
    h = plot (quartile_x, quartile_y, 'b;;',
              whisker_x, whisker_y, 'b;;',
              cap_x, cap_y, 'b;;',
              median_x, median_y, 'r;;',
              outliers_x, outliers_y, [symbol(1), 'r;;'],
              outliers2_x, outliers2_y, [symbol(2), 'r;;'], plot_opts{:});
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
   ## Define outlier_tags' horizontal alignment
    outlier_tags_alignment = {'horizontalalignment', 'left', 'rotation', 90};
      if (box_style)
        f = fillbox (quartile_y, quartile_x, bcolor);
      endif
      h = plot (quartile_y, quartile_x, 'b;;',
               whisker_y, whisker_x, 'b;;',
               cap_y, cap_x, 'b;;',
               median_y, median_x, 'r;;',
               outliers_y, outliers_x, [symbol(1), 'r;;'],
               outliers2_y, outliers2_x, [symbol(2), 'r;;'], plot_opts{:});
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

  ## Distribute the handles that plot returned.  They come back in the order
  ## the segments were passed and a segment with no columns contributes none,
  ## so walk them with a running offset.  Chaining each block off the last index
  ## of the one before it broke as soon as a block was empty: a single-point,
  ## an all-NaN, or an empty data set leaves no box at all, and the chain then
  ## indexed the previous block at zero.
  n_box  = columns (quartile_x);
  n_whis = 2 * columns (whisker_x);   # the caps repeat the whisker columns
  n_med  = columns (median_x);
  ## The outliers of every group are plotted as one series, so they account for
  ## a single handle when there are any.  Counting their columns would say one
  ## for the 0x1 empty they collapse to when there are none.
  n_out  = double (! isempty (outliers_y));
  n_out2 = double (! isempty (outliers2_y));
  used   = 0;

  ## Box outlines and box fill (if any).  The fill follows the boxes that were
  ## actually drawn, which is not the number of groups once any were chopped.
  hs.box = h(used + [1 : n_box]);
  used += n_box;
  if (box_style)
    hs.box_fill = f(1 : n_box);
  else
    hs.box_fill = [];
  endif

  ## Whiskers (including caps) and median lines
  hs.whisker = h(used + [1 : n_whis]);
  used += n_whis;
  hs.median = h(used + [1 : n_med]);
  used += n_med;

  ## Outliers (if any) and their respective tags (if applicable)
  if (n_out > 0)
    hs.outliers = h(used + [1 : n_out]);
    used += n_out;
    if (outlier_tags == 1)
      hs.out_tags = t1(1 : length (outliers_tags_y));
    else
      hs.out_tags = [];
    endif
  else
    hs.outliers = [];
    hs.out_tags = [];
  endif

  ## Extreme outliers (if any) and their respective tags (if applicable)
  if (n_out2 > 0)
    hs.outliers2 = h(used + [1 : n_out2]);
    used += n_out2;
    if (outlier_tags == 1)
      hs.out_tags2 = t2(1 : length (outliers2_tags_y));
    else
      hs.out_tags2 = [];
    endif
  else
    hs.outliers2 = [];
    hs.out_tags2 = [];
  endif

  ## Redraw the median lines to avoid colour overlapping in case of 'filled'
  ## BoxStyle
  if (box_style)
    set (hs.median, 'color', 'r');
  endif

  ## Print labels according to orientation and return handle
  if (orientation)
    set (gca (), 'xtick', groups, 'xticklabel', labels);
    hs.labels = get (gcf, 'currentaxes');
  else
    set (gca (), 'ytick', groups, 'yticklabel', labels);
    hs.labels = get (gcf, 'currentaxes');
  endif

  ## retain original ishold status.
  if (! hold_status)
    hold off;
  endif

  ## return output arguments if desired.
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


function f = fillbox (quartile_y, quartile_x, bcolor)

  f = [];
  for icol = 1 : columns (quartile_x)
    if (ischar (bcolor))
      f = [ f; fill(quartile_y(:, icol), quartile_x(:, icol), ...
                    bcolor(mod (icol - 1, numel (bcolor)) + 1)) ];
    else
      f = [ f; fill(quartile_y(:, icol), quartile_x(:, icol), ...
                    bcolor(mod (icol - 1, size (bcolor, 1)) + 1, :)) ];
    endif
    hold on;
  endfor

endfunction

%!demo
%! axis ([0, 3]);
%! randn ('seed', 1);    # for reproducibility
%! girls = randn (10, 1) * 5 + 140;
%! randn ('seed', 2);    # for reproducibility
%! boys = randn (13, 1) * 8 + 135;
%! boxplot ({girls, boys});
%! set (gca (), 'xtick', [1 2], 'xticklabel', {'girls', 'boys'})
%! title ('Grade 3 heights');

%!demo
%! randn ('seed', 7);    # for reproducibility
%! A = randn (10, 1) * 5 + 140;
%! randn ('seed', 8);    # for reproducibility
%! B = randn (25, 1) * 8 + 135;
%! randn ('seed', 9);    # for reproducibility
%! C = randn (20, 1) * 6 + 165;
%! data = [A; B; C];
%! groups = [(ones (10, 1)); (ones (25, 1) * 2); (ones (20, 1) * 3)];
%! labels = {'Team A', 'Team B', 'Team C'};
%! pos = [2, 1, 3];
%! boxplot (data, groups, 'Notch', 'on', 'Labels', labels, 'Positions', pos, ...
%!          'OutlierTags', 'on', 'BoxStyle', 'filled');
%! title ('Example of Group splitting with paired vectors');

%!demo
%! randn ('seed', 1);    # for reproducibility
%! data = randn (100, 9);
%! boxplot (data, 'notch', 'on', 'boxstyle', 'filled', ...
%!          'colors', 'ygcwkmb', 'whisker', 1.2);
%! title ('Example of different colors specified with characters');

%!demo
%! randn ('seed', 5);    # for reproducibility
%! data = randn (100, 13);
%! colors = [0.7 0.7 0.7; ...
%!           0.0 0.4 0.9; ...
%!           0.7 0.4 0.3; ...
%!           0.7 0.1 0.7; ...
%!           0.8 0.7 0.4; ...
%!           0.1 0.8 0.5; ...
%!           0.9 0.9 0.2];
%! boxplot (data, 'notch', 'on', 'boxstyle', 'filled', ...
%!          'colors', colors, 'whisker', 1.3, 'boxwidth', 'proportional');
%! title ('Example of different colors specified as RGB values');

%!demo
%! randn ('seed', 11); # for reproducibility
%! data = randn (30, 1);
%! ## Using modern string arrays
%! str_groups = string (repmat (['Control'; 'TreatmentA'; 'TreatmentB'], 10, 1));
%! boxplot (data, str_groups, 'colors', 'rgb');
%! title ('Example using modern string arrays for grouping');

%!demo
%! randn ('seed', 10); # for reproducibility
%! data = randn (40, 1) * 5 + 50;
%! ## Create two different grouping variables
%! group1 = repmat ({'Alpha'; 'Beta'}, 20, 1);
%! group2 = repmat ([2022; 2022; 2023; 2023], 10, 1);
%! ## Pass them together as a cell array
%! boxplot (data, {group1, group2});
%! title ('Example of Multiple Grouping Variables (Model & Year)');

## Input data validation
%!error <numerical array or cell array containing> boxplot ('a')
%!error <data cells must contain> boxplot ({[1 2 3], 'a'})
%!error <grouping vector may only be passed> boxplot ([1 2 3], 1, {2, 3})
%!error <all grouping variables must have the same length> boxplot ([1 2 3], {'a', 'b'})
%!error <'Notch' input argument accepts> boxplot ([1:10], 'notch', 'any')
%!error <boxplot: invalid 'Notch' value.> boxplot ([1:10], 'notch', i)
%!error <boxplot: invalid 'Notch' value.> boxplot ([1:10], 'notch', {})
%!error <must be character> boxplot (1, 'symbol', 1)
%!error <'Orientation' input argument accepts only> boxplot (1, 'orientation', 'diagonal')
%!error <boxplot: invalid 'Orientation' value.> boxplot (1, 'orientation', {})
%!error <'Whisker' input argument accepts only> boxplot (1, 'whisker', 'a')
%!error <'Whisker' input argument accepts only> boxplot (1, 'whisker', [1 3])
%!error <'OutlierTags' input argument accepts only> boxplot (3, 'OutlierTags', 'maybe')
%!error <boxplot: invalid 'OutlierTags' value.> boxplot (3, 'OutlierTags', {})
%!error <'Sample_IDs' input argument accepts only> boxplot (1, 'sample_IDs', 1)
%!error <'BoxWidth' input argument accepts only> boxplot (1, 'boxwidth', 2)
%!error <'BoxWidth' input argument accepts only> boxplot (1, 'boxwidth', 'anything')
%!error <'Widths' input argument accepts only> boxplot (5, 'widths', 'a')
%!error <'Widths' input argument accepts only> boxplot (5, 'widths', [1:4])
%!error <'Widths' input argument accepts only> boxplot (5, 'widths', [])
%!error <'CapWidths' input argument accepts only> boxplot (5, 'capwidths', 'a')
%!error <'CapWidths' input argument accepts only> boxplot (5, 'capwidths', [1:4])
%!error <'CapWidths' input argument accepts only> boxplot (5, 'capwidths', [])
%!error <'BoxStyle' input argument accepts only> boxplot (1, 'Boxstyle', 1)
%!error <'BoxStyle' input argument accepts only> boxplot (1, 'Boxstyle', 'garbage')
%!error <'Positions' input argument accepts only> boxplot (1, 'positions', 'aa')
%!error <'Labels' input argument accepts only> boxplot (3, 'labels', [1 5])
%!error <'Colors' input argument accepts only> boxplot (1, 'colors', {})
%!error <'Colors' input argument accepts only> boxplot (2, 'colors', [1 2 3 4])
%!error <Sample_IDs must match the data> boxplot (randn (10, 3), 'Sample_IDs', {'a', 'b'})
%!error <with the formalism> boxplot (rand (3, 3), [1 2])

## Test plotting
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   [a, b] = boxplot (rand (10, 3));
%!   assert_equal (size (a), [7, 3]);
%!   assert_equal (numel (b.box), 3);
%!   assert_equal (numel (b.whisker), 12);
%!   assert_equal (numel (b.median), 3);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   [~, b] = boxplot (rand (10, 3), 'BoxStyle', 'filled', 'colors', 'ybc');
%!   assert_equal (numel (b.box_fill), 3);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   hold on
%!   [a, b] = boxplot (rand (10, 3));
%!   assert_equal (ishold, true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! ## Test multi-variable grouping.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   data = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12];
%!   g1 = [1; 1; 2; 2; 3; 3; 1; 1; 2; 2; 3; 3];
%!   g2 = string ({'A'; 'B'; 'A'; 'B'; 'A'; 'B'; 'A'; 'B'; 'A'; 'B'; 'A'; 'B'});
%!   g3 = categorical ({'X'; 'X'; 'Y'; 'Y'; 'Z'; 'Z'; 'X'; 'X'; 'Y'; 'Y'; 'Z'; 'Z'});
%!   [a, b] = boxplot (data, {g1, g2, g3});
%!   assert_equal (size (a, 2), 6);
%! unwind_protect_cleanup
%!     close (hf);
%! end_unwind_protect
%!test
%! ## Test multi-variable grouping with empty intersections dropping correctly.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   data = [1; 2; 3; 4];
%!   g1 = [1; 1; 2; 2];
%!   g2 = string ({'A'; 'A'; 'B'; 'B'});
%!   [a, b] = boxplot (data, {g1, g2});
%!   assert_equal (size (a, 2), 2);
%! unwind_protect_cleanup
%!     close (hf);
%! end_unwind_protect

## Inputs that used to die inside boxplot rather than be plotted or refused.

%!test
%! ## A single observation is plotted, not an internal subscript error.  Every
%! ## statistic collapses to the value, as it does in MATLAB.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   [s, hs] = boxplot (5);
%!   assert_equal (s, 5 * ones (7, 1));
%!   assert_equal (isempty (hs.box), true);
%!   assert_equal (numel (hs.outliers), 1);
%!   assert_equal (get (hs.outliers, 'YData'), 5);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## An empty input draws nothing and returns no statistics, as in MATLAB,
%! ## rather than indexing an empty handle list at zero.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   [s, hs] = boxplot ([]);
%!   assert_equal (size (s, 2), 0);
%!   assert_equal (isempty (hs.box), true);
%!   assert_equal (isempty (hs.whisker), true);
%!   assert_equal (isempty (hs.median), true);
%!   assert_equal (isempty (hs.outliers), true);
%!   clf;
%!   [s, hs] = boxplot ({});
%!   assert_equal (size (s, 2), 0);
%!   assert_equal (isempty (hs.box), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## A variable that is entirely missing leaves no box, which used to take the
%! ## handle bookkeeping with it.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   [s, hs] = boxplot ([NaN; NaN; NaN]);
%!   assert_equal (all (isnan (s)), true);
%!   assert_equal (isempty (hs.box), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## Integer observations are accepted and agree with the same numbers in
%! ## double; the internal call to var used to refuse them outright.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   sd = boxplot ([1; 2; 3; 4; 5; 6; 7]);
%!   clf;
%!   si = boxplot (int32 ([1; 2; 3; 4; 5; 6; 7]));
%!   clf;
%!   su = boxplot (uint8 ([1; 2; 3; 4; 5; 6; 7]));
%!   assert_equal (si, sd);
%!   assert_equal (su, sd);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## The same holds for an integer variable inside a cell.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   sd = boxplot ({[1; 2; 3; 4; 5], [2; 3; 4; 5; 6]});
%!   clf;
%!   si = boxplot ({int32([1; 2; 3; 4; 5]), int32([2; 3; 4; 5; 6])});
%!   assert_equal (si, sd);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## A group holding a single observation leaves fewer boxes than groups.  The
%! ## cap widths were built from the full-length widths vector against the
%! ## chopped box vector, so the subtraction did not conform.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   data = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 42];
%!   grp  = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
%!   [s, hs] = boxplot (data, grp);
%!   assert_equal (size (s, 2), 2);
%!   assert_equal (s(:, 2), 42 * ones (7, 1));
%!   assert_equal (numel (hs.box), 1);
%!   assert_equal (numel (hs.median), 1);
%!   assert_equal (numel (hs.whisker), 4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## The box fill follows the boxes actually drawn, not the number of groups.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   data = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 42];
%!   grp  = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
%!   [s, hs] = boxplot (data, grp, 'BoxStyle', 'filled');
%!   assert_equal (numel (hs.box_fill), numel (hs.box));
%!   assert_equal (get (hs.box_fill(1), 'Type'), 'patch');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## Both kinds of outlier keep their own handles, and the ones after them
%! ## stay in step.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   data = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 42];
%!   grp  = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2];
%!   [s, hs] = boxplot ([data; 60], [grp; 1]);
%!   assert_equal (numel (hs.outliers), 1);
%!   assert_equal (numel (hs.outliers2), 1);
%!   assert_equal (get (hs.median, 'Type'), 'line');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Logical and character data stay refused, as they are in MATLAB.
%!error <numerical array or cell array containing data expected.> ...
%! boxplot ([true; false; true; true])
%!error <numerical array or cell array containing data expected.> ...
%! boxplot ('abcde')

## Quartiles follow quantile's own default method, which is prctile's and
## MATLAB's.  Values below are MATLAB R2024a's, read off the box, whisker and
## outlier objects it draws.

%!test
%! ## The box edges are the 25th and 75th percentiles, for odd and even counts
%! ## alike, and agree with prctile called on the same data.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = boxplot ((1:7)');
%!   assert_equal (s(1:5)', [1, 2.25, 4, 5.75, 7]);
%!   clf;
%!   s = boxplot ((1:8)');
%!   assert_equal (s(1:5)', [1, 2.5, 4.5, 6.5, 8]);
%!   clf;
%!   x = [2; 4; 4; 4; 5; 5; 7; 9];
%!   s = boxplot (x);
%!   assert_equal (s(1:5)', [2, 4, 4.5, 6, 9]);
%!   assert_equal (s(2:4)', [prctile(x, 25), median(x), prctile(x, 75)]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## Small samples, where the two quantile definitions differ most.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = boxplot ([1; 2]);
%!   assert_equal (s(1:5)', [1, 1, 1.5, 2, 2]);
%!   clf;
%!   s = boxplot ([1; 2; 3]);
%!   assert_equal (s(1:5)', [1, 1.25, 2, 2.75, 3]);
%!   clf;
%!   s = boxplot ([1; 2; 3; 4]);
%!   assert_equal (s(1:5)', [1, 1.5, 2.5, 3.5, 4]);
%!   clf;
%!   s = boxplot ([1; 2; 3; 4; 5]);
%!   assert_equal (s(1:5)', [1, 1.75, 3, 4.25, 5]);
%!   clf;
%!   s = boxplot ([1; 2; 3; 4; 5; 6]);
%!   assert_equal (s(1:5)', [1, 2, 3.5, 5, 6]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## Missing values are dropped before the quartiles are taken.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   s = boxplot ([1; 2; NaN; 4; 5; 6; 7]);
%!   assert_equal (s(1:5)', [1, 2, 4.5, 6, 7]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## The inter-quartile range sets the fences, so the quartile definition
%! ## decides which points are outliers.  33 lies beyond the fence and 31.5
%! ## does not; the old quartiles put the fence between them and called both
%! ## outliers.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   [s, hs] = boxplot ([(1:20)'; 31.5]);
%!   assert_equal (s(1:5)', [1, 5.75, 11, 16.25, 31.5]);
%!   assert_equal (isempty (hs.outliers), true);
%!   assert_equal (isempty (hs.outliers2), true);
%!   clf;
%!   [s, hs] = boxplot ([(1:20)'; 33]);
%!   assert_equal (s(1:5)', [1, 5.75, 11, 16.25, 33]);
%!   assert_equal (isempty (hs.outliers) && isempty (hs.outliers2), false);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## When no observation lies between a quartile and its fence the whisker
%! ## collapses onto the quartile instead of being drawn back into the box.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   [s, hs] = boxplot ([1; 1; 2; 2; 3; 3; 4; 40; 50]);
%!   assert_equal (s(1:5)', [1, 1.75, 3, 13, 50]);
%!   yd = get (hs.whisker, 'YData');
%!   assert_equal (max ([yd{3}, yd{4}]), 13);
%!   clf;
%!   [s, hs] = boxplot ([-50; -40; -4; -3; -3; -2; -2; -1; -1]);
%!   assert_equal (s(1:5)', [-50, -13, -3, -1.75, -1]);
%!   yd = get (hs.whisker, 'YData');
%!   assert_equal (min ([yd{1}, yd{2}]), -13);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test
%! ## A whisker that does reach real data is not clamped.
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   [s, hs] = boxplot ([-50; -40; -4; -3; -3; -2; -2; 40; 50]);
%!   assert_equal (s(1:5)', [-50, -13, -3, 8.5, 50]);
%!   yd = get (hs.whisker, 'YData');
%!   assert_equal (min ([yd{1}, yd{2}]), -40);
%!   assert_equal (max ([yd{3}, yd{4}]), 40);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
