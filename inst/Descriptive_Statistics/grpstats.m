## Copyright (C) 2022-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 jayantchauhan <0001jayant@gmail.com>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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
## @deftypefn  {statistics} {@var{stats} =} grpstats (@var{x})
## @deftypefnx {statistics} {@var{stats} =} grpstats (@var{x}, @var{group})
## @deftypefnx {statistics} {[@var{means}, @var{sem}, @var{counts}, @var{gname}] =} @
## grpstats (@var{x}, @var{group})
## @deftypefnx {statistics} {[@var{stats1}, @dots{}, @var{statsN}] =} grpstats @
##(@var{x}, @var{group}, @var{whichstats})
## @deftypefnx {statistics} {[@var{stats1}, @dots{}, @var{statsN}] =} grpstats @
## (@var{x}, @var{group}, @var{whichstats}, @qcode{'Alpha'}, @var{alpha})
## @deftypefnx {statistics} {@var{tblstats} =} grpstats (@var{tbl})
## @deftypefnx {statistics} {@var{tblstats} =} grpstats (@var{tbl}, @var{groupvars})
## @deftypefnx {statistics} {@var{tblstats} =} grpstats (@var{tbl}, @var{groupvars}, @var{whichstats})
## @deftypefnx {statistics} {@var{tblstats} =} grpstats (@var{tbl}, @var{groupvars}, @
## @var{whichstats}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {} grpstats (@var{x}, @var{group}, @var{alpha})
## @deftypefnx {statistics} {@var{h} =} grpstats (@var{x}, @var{group}, @var{alpha})
##
## Summary statistics by group.
##
## @code{grpstats} computes groupwise summary statistics for the data in
## @var{x}, which can be a numeric matrix or a table.  Numeric vectors are
## treated as a single column matrix.  @qcode{NaN}s are treated as missing
## values and removed from calculations.
##
## @subheading Syntax for Numeric Input
##
## @code{@var{stats} = grpstats (@var{x})} calculates the mean statistic for
## each column in @var{x} and returns it as row vector in @var{stats}.
##
## @code{@var{stats} = grpstats (@var{x}, @var{group})} calculates the mean
## statistic for each column in @var{x} grouped by @var{group}.  The returned
## argument, @var{stats}, is also a matrix with equal columns as @var{x} and the
## number of rows is equal to the groups specified by @var{group}.
##
## The grouping variable, @var{group} can be a vector of any data type
## supported by the @code{grp2idx} function.  Alternatively, it can be a cell
## vector specifying multiple grouping variables with each cell element
## containing any of the aforementioned supported grouping vectors.  If
## @var{group} is empty (@code{[]}), then input @var{x} is treated as a single
## group.  A group whose values are all missing is dropped, so a grouping
## variable that is entirely missing leaves no groups at all.
##
## @code{[@var{means}, @var{sem}, @var{counts}, @var{gname}] = grpstats
## (@var{x}, @var{group})} returns, without naming any statistic, the mean, the
## standard error of the mean, the number of non-@qcode{NaN} values and the
## group names, in that order.  Up to four outputs may be requested this way.
##
## An empty @var{x} produces empty output with as many columns as @var{x} and
## no rows, one row per group being the general rule and an empty @var{x}
## having no groups.
##
## @code{[@var{stats1}, @dots{}, @var{statsN}] = grpstats (@var{x}, @var{group},
## @var{whichstats})} calculates the summary statistics specified by the
## @var{whichstats} argument, which can include any of the available statistics
## shown below.  The number of output arguments must match the number of
## requested statistics specified in @var{whichstats}.
## computes summary statistics for the numeric matrix @var{x} grouped by
## @var{group}.
##
## @var{x} must be a numeric vector or a 2-D matrix.  Vectors are treated as
## a single-column matrix.
##
## @var{group} is a grouping variable that defines the groups for the rows of
## @var{x}.  It can be a categorical variable, numeric vector, string array, or
## cell array of strings.  @var{group} can also be a cell array containing
## multiple grouping variables.  If @var{group} is empty (@code{[]}) or omitted,
## all of @var{x} is treated as a single group.
##
## @var{whichstats} specifies the statistics to compute. It can be either a
## string array or a cell array of strings specifying any of the following
## builtin statistics.  If omitted, the default is @qcode{'mean'}.
## @var{whichstats} can also contain function handles for custom statistics.
##
## The available statistics are:
## @multitable @columnfractions 0.2 0.75
## @item @qcode{'mean'} @tab Mean of each group.
## @item @qcode{'median'} @tab Median of each group.
## @item @qcode{'sem'} @tab Standard error of the mean for each group.
## @item @qcode{'std'} @tab Standard deviation of each group.
## @item @qcode{'var'} @tab Variance of each group.
## @item @qcode{'min'} @tab Minimum value in each group.
## @item @qcode{'max'} @tab Maximum value in each group.
## @item @qcode{'range'} @tab Difference between max and min in each
## group.
## @item @qcode{'numel'} @tab Number of elements (count) in each group.
## @item @qcode{'meanci'} @tab Confidence interval for the mean.
## @item @qcode{'predci'} @tab Prediction interval for a new observation.
## @item @qcode{'gname'} @tab Group names.
## @end multitable
##
## @code{[@dots{}] = grpstats (@dots{}, @qcode{'Alpha'}, @var{alpha})} specifies
## the significance level for the confidence intervals (@qcode{'meanci'} and
## @qcode{'predci'}) as @code{100 * (1-@var{alpha})@@%}.  @var{alpha} must be a
## scalar between 0 and 1.  When not specified, it defaults to 0.05.  Note that
## this paired input argument is also valid for table input.
##
## @subheading Syntax for Table Input
##
## @code{@var{tblstats} = grpstats (@var{tbl}, @var{groupvars})} computes the
## summary statistics for the data in table @var{tbl}, grouped by the variables
## specified in @var{groupvars}.  If @var{groupvars} is empty or omitted, then
## all of @var{tbl} is treated as a single group.  @var{groupvars} can be a cell
## array of character vectors or a string array specifying one or more variable
## names in @var{tbl} to be used as grouping variables.  Alternatively, all
## valid methods for indexing table variables are supported (e.g. @code{vartype}
## object, logical vector, function handle).
##
## The output @var{tblstats} is a table with one row for each group. It contains
## the grouping variables, an additional  @qcode{'GroupCount'} variable, and the
## specified summary statistics for the variables in @var{tbl}, expect for those
## specified as grouping variables.  When input is a table, only a single output
## variable, @var{tblstats} can be specified.  The output @var{tblstats} also
## contains @qcode{RowNames}, which are the unique combinations of the specified
## groups, for which data are available in @var{tbl}.  When no groups are
## specified, the row name of the single row output table defaults to
## @qcode{'All'}.
##
## @code{@var{tblstats} = grpstats (@var{tbl}, @var{groupvars},
## @var{whichstats})} specifies which statistics to calculate for the variables
## in @var{tbl}.  Unless specified, the mean is calculated for each variable.
## When specifying more than one statistic, @var{tblstats} contains multiple
## variables for each variable in @var{tbl} and each is named by combining the
## applied statistic with the name of the original variable.  When a function
## handle is applied, its string representation is used instead.
##
## For table input specifically, @code{grpstats} also accepts the following
## paired arguments.
##
## @multitable @columnfractions 0.2 0.75
## @headitem Name @tab Value
## @item @qcode{'DataVars'} @tab A vector specifying the variables in
## @var{tbl}, for which to calculate the specified statistics.  The vector can
## be any of the valid options for indexing table variables.
##
## @item @qcode{'VarNames'} @tab A cell array of character vectors or a
## string array specifying the names of the variables in the output table.  The
## number of specified names must match the number of expected variables in the
## output table.
## @end multitable
##
## @subheading Plotting Syntax
##
## The syntax @code{grpstats (@var{x}, @var{group}, @var{alpha})} generates an
## @code{errorbar} plot with the group means and their respective confidence
## intervals.  @var{x} must be a numeric vector or matrix.  @var{alpha} is a
## scalar between 0 and 1 that determines the confidence level.  This syntax is
## an alternative to calling @code{errorbar} after computing @qcode{'mean'} and
## @qcode{'meanci'} statistics.  The optional output @var{h} is a handle to the
## hggroup object representing the data plot and errorbars.
##
## @seealso{grp2idx}
## @end deftypefn

function [varargout] = grpstats (x, group = [], whichstats = [], varargin)

  ## Check data input
  if (nargin < 1)
    print_usage ();
  endif
  if (ndims (x) != 2)
     error ("grpstats: X must be a matrix or a table.");
  endif
  if (! istable (x) && isvector (x))
    x = x(:);
  endif
  [r, c] = size (x);

  ## For table input no more than 1 output argument is allowed
  if (istable (x) && nargout > 1)
    error ("grpstats: only one output argument in allowed when X is a table.");
  endif

  ## Add default grouping variable
  no_group = true;
  if (isempty (group))
    grp_idx = ones (r, 1);
    if (r == 0)
      r_names = cell (0, 1);
      g_names = cell (0, 1);
      ngroups = 0;
    else
      r_names = {'All'};
      g_names = {'1'};
      ngroups = 1;
    endif
    no_group = false;
  endif

  ## Check for plotting functional form with three input arguments
  if (nargin == 3 && isnumeric (whichstats) && isscalar (whichstats) ...
      && whichstats > 0 && whichstats < 1)
    do_plot = true;
    alpha = whichstats;
    whichstats = [];
    fcn_names = {'mean', 'meanci'};
    if (! (isnumeric (x)))
      error ("grpstats: X must be numeric to plot mean and CI for each group.");
    endif
  else
    do_plot = false;
  endif

  ## Parse statistical functions
  if (isempty (whichstats))
    ## Without WHICHSTATS the outputs are the mean, the standard error, the
    ## group counts and the group names, in that order.
    default_fcn = {'mean', 'sem', 'numel', 'gname'};
    if (nargout > numel (default_fcn))
      error ("grpstats: too many output arguments.");
    endif
    fcn_names = default_fcn(1:max (1, nargout));
  else
    if (ischar (whichstats) || isstring (whichstats))
      fcn_names = cellstr (whichstats);
    elseif (iscellstr (whichstats))
      fcn_names = whichstats;
    elseif (is_function_handle (whichstats))
      fcn_names = {whichstats};
    elseif (iscell (whichstats))
      TF = cellfun (@(x) ischar (x) || is_function_handle (x), whichstats);
      if (all (TF))
        fcn_names = whichstats;
      else
        error ("grpstats: invalid WHICHSTATS specification in cell array.");
      endif
    else #if (! (isscalar (whichstats) && whichstats > 0 && whichstats < 1))
      error ("grpstats: invalid WHICHSTATS data type.");
    endif
    ## At this point we have cell array with either function names or
    ## function handles which need to be tested before going any further
    valid_fcn = {'mean', 'median', 'sem', 'std', 'var', 'min', 'max', ...
                 'range', 'numel', 'meanci', 'predci', 'gname'};
    is_char = cellfun ('ischar', fcn_names);
    if (any (is_char))
      isvalid = cellfun (@(x) ismember (x, valid_fcn), fcn_names(is_char));
      if (! all (isvalid))
        error ("grpstats: unrecognized function names in WHICHSTATS.");
      endif
    endif
    ## Function handles need to be tested on the actual data with try..catch
    ## blocks, because there are numerous combination in I/O size which does
    ## not worth to test preemptively.
  endif

  ## Parse optional Name-Value paired arguments
  optNames = {'Alpha', 'DataVars', 'VarNames'};
  dfValues = {0.05, [], {}};
  [alpha, DataVars, VarNames, args] = parsePairedArguments (optNames, dfValues, ...
                                                  varargin(:));
  if (! isempty (args))
    tmp = args{1};
    if (isscalar (tmp) && tmp > 0 && tmp < 1)
      alpha = tmp;
    else
      error ("grpstats: unrecognized input arguments.");
    endif
  endif
  ## Check alpha value
  if (! (isscalar (alpha) && alpha > 0 && alpha < 1))
    error ("grpstats: 'alpha' must be a real scalar in the range (0,1).");
  endif
  ## Force VarNames to cell array of character vectors
  if (! isempty (VarNames))
    if (ischar (VarNames) || isstring (VarNames))
      VarNames = cellstr (VarNames);
    endif
    if (! iscellstr (VarNames))
      error ("grpstats: invalid data types for 'VarNames'.");
    endif
  endif

  ## Handle group and whichstats for matrices and tables
  if (istable (x))
    ## Get grouping variables for table input data
    if (no_group)
      try
        grp_vars = x(:, group);  # () returns groupvars as a table
      catch
        error ("grpstats: cannot resolve GROUPVARS in input table.");
      end_try_catch
      ## Get unique groups and their indices to table data
      [g_names, ~, grp_idx] = unique (grp_vars);
      ## Get number of groups
      ngroups = rows (g_names);
      ## Convert logical to double
      c_names = convertvars (g_names, @islogical, 'double');
      ## Create unique names for each group
      cstrtbl = cellstr (string (table2cell (c_names)));
      r_names = cstrtbl(:,1);
      [cr, cc] = size (cstrtbl);
      tmp_tmp = repmat ({'_'}, cr, 1);
      for idx = 2:cc
        r_names = strcat (r_names, tmp_tmp, cstrtbl(:,idx));
      endfor
    endif

    ## Compute group count by default in tables
    GroupCount = zeros (ngroups, 1);
    for idx = 1:ngroups
      GroupCount(idx,:) = sum (grp_idx == idx);
    endfor

    ## Add group count and row names
    if (no_group)
      g_names = addvars (g_names, GroupCount);
      g_names.Properties.RowNames = r_names;
    else
      g_names = table (GroupCount, 'RowNames', r_names);
    endif

    ## Check that DataVars exist (if given)
    if (! isempty (DataVars))
      ## Octave extension allows DataVars to be specified as a vartype object
      ## or a function handle to select variables based on their data types.
      ## To avoid multiple checking for all possible referencing schemes,
      ## let's add another try...catch blocks
      try
        work_tbl = x(:, DataVars);
      catch
        error ("grpstats: invalid 'DataVars' reference to table X.");
      end_try_catch
    elseif (isempty (group))
      ## No grouping variables, so every variable of the table gets used
      work_tbl = x;
    else
      ## No specified DataVars, all table except grouping variables gets used
      work_tbl = removevars (x, group);
    endif

    ## From this point we can start applying functions on table variables

    ## We have to apply each function to every working variable
    ## We name the variable accordingly and append the output of the
    ## function at the g_names table.
    ## The function is applied on each group of each variable
    ## Handle functions need try..catch block to ensure they are compatible
    for var_idx = 1:columns (work_tbl)
      vname = work_tbl.Properties.VariableNames(var_idx);
      vdata = work_tbl{:,var_idx};
      ncols = columns (vdata);
      for fcn_idx = 1:numel (fcn_names)
        new_vdata = [];
        fcn_op = fcn_names{fcn_idx};
        if (is_function_handle (fcn_op))
          fname = {func2str(fcn_op)};
          new_vname = strcat (fname, '_', vname);
          try
            new_vdata = applyhandle (vdata, grp_idx, ngroups, ncols, fcn_op);
          catch
            # no-op
          end_try_catch
          if (! isempty (new_vdata))
            g_names = addvars (g_names, new_vdata, 'NewVariableNames', new_vname);
          endif
        else  # it must be a function name
          switch (fcn_op)
            case 'mean'
              new_vname = strcat ('mean_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = mean (gdata, 1, 'omitnan');
              endfor
            case 'median'
              new_vname = strcat ('median_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = median (gdata, 1, 'omitnan');
              endfor
            case 'sem'
              new_vname = strcat ('sem_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = std (gdata, 0, 1, 'omitnan') ./ ...
                                   sqrt (size (gdata, 1) - sum (isnan (gdata), 1));
              endfor
            case 'std'
              new_vname = strcat ('std_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = std (gdata, 0, 1, 'omitnan');
              endfor
            case 'var'
              new_vname = strcat ('var_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = var (gdata, 0, 1, 'omitnan');
              endfor
            case 'min'
              new_vname = strcat ('min_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = min (gdata, [], 1);
              endfor
            case 'max'
              new_vname = strcat ('max_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = max (gdata, [], 1);
              endfor
            case 'range'
              new_vname = strcat ('range_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = range (gdata, 1);
              endfor
            case 'numel'
              new_vname = strcat ('numel_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = size (gdata, 1) - sum (isnan (gdata), 1);
              endfor
            case 'meanci'
              new_vname = strcat ('meanci_', vname);
              ## Preallocate twice the columns in variable (lower, upper)
              new_vdata = NaN (ngroups, ncols * 2);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                m = mean (gdata, 1, 'omitnan');
                n = size (gdata, 1) - sum (isnan (gdata), 1);
                s = std (gdata, 0, 1, 'omitnan') ./ sqrt (max (n,1));
                ## Avoid invalid tinv calls for degenerate df
                df = max (n - 1, 0);
                tval = zeros (1, size (gdata, 2));
                pos = (df > 0);
                if (any (pos))
                  tval(pos) = - tinv (alpha / 2, df(pos));
                endif
                d = s .* tval;
                d(n < 2) = NaN;
                new_vdata(idx,[1:2:end]) = m - d;
                new_vdata(idx,[2:2:end]) = m + d;
              endfor

            case 'predci'
              new_vname = strcat ('predci_', vname);
              ## Preallocate twice the columns in variable (lower, upper)
              new_vdata = NaN (ngroups, ncols * 2);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                m = mean (gdata, 1, 'omitnan');
                n = size (gdata, 1) - sum (isnan (gdata), 1);
                s = std (gdata, 0, 1, 'omitnan') .* sqrt (1 + (1 ./ max (n,1)));
                df = max (n - 1, 0);
                tval = zeros (1, size (gdata, 2));
                pos = (df > 0);
                if (any (pos))
                  tval(pos) = - tinv (alpha / 2, df(pos));
                endif
                d = s .* tval;
                d(n < 2) = NaN;
                new_vdata(idx,[1:2:end]) = m - d;
                new_vdata(idx,[2:2:end]) = m + d;
              endfor
          endswitch
          if (! isempty (new_vdata))
            g_names = addvars (g_names, new_vdata, 'NewVariableNames', new_vname);
          endif
        endif
      endfor
    endfor

    ## Last check. VarNames must equal the number of expected variables
    if (! isempty (VarNames))
      if (numel (VarNames) != columns (g_names))
        error ("grpstats: 'VarNames' do not match expected variables.");
      endif
      g_names = renamevars (g_names, g_names.Properties.VariableNames, ...
                            VarNames);
    endif

    varargout{1} = g_names;
    return;
  else
    ## Get groups for array input data
    if (no_group)
      if (iscell (group) && ! iscellstr (group))
        ## Multiple grouping variables in cell array
        grp_vars = numel (group);
        [grp_idx, g_names1] = grp2idx (group{1});
        if (numel (grp_idx) != r)
          error ("grpstats: samples in X and GROUPS mismatch.");
        endif
        all_g_names = {g_names1};
        if (grp_vars > 1)
          for g_idx = 2:grp_vars
            [tmp_grp_idx, tmp_g_names] = grp2idx (group{g_idx});
            if (numel (tmp_grp_idx) != r)
              error ("grpstats: samples in X and GROUPS mismatch.");
            endif
            grp_idx = [grp_idx, tmp_grp_idx];
            all_g_names{g_idx} = tmp_g_names;
          endfor
        endif
        ## Get combination of unique groups and their common index to X.
        ## Rows missing in any grouping variable belong to no group.
        keep = ! any (isnan (grp_idx), 2);
        [g_names_idx, ~, keep_idx] = unique (grp_idx(keep,:), 'rows');
        grp_idx = NaN (r, 1);
        grp_idx(keep) = keep_idx;
        ngroups = rows (g_names_idx);
        c_names = cell (ngroups, grp_vars);
        if (ngroups > 0)
          for gvar_idx = 1:grp_vars
            gn = all_g_names{gvar_idx};
            c_names(:, gvar_idx) = gn(g_names_idx(:,gvar_idx));
          endfor
        endif
        g_names = c_names;
      else
        [grp_idx, g_names] = grp2idx (group);
        ngroups = numel (g_names);
        if (numel (grp_idx) != r)
          error ("grpstats: samples in X and GROUPS mismatch.");
        endif
      endif
    endif

    ## Check for plot option
    if (do_plot)
      if (ngroups == 0)
        error ("grpstats: no groups to plot.");
      endif
      ## Calculate mean and ci
      mu = NaN (ngroups, c);
      ci = NaN (ngroups, c);
      for idx = 1:ngroups
        group_x = x(find (grp_idx == idx), :);
        mu(idx,:) = mean (group_x, 1, 'omitnan');
        n = size (group_x, 1) - sum (isnan (group_x), 1);
        s = std (group_x, 0, 1, 'omitnan') ./ sqrt (max (n,1));
        ## Avoid invalid tinv calls for degenerate df
        df = max (n - 1, 0);
        tval = zeros (1, size (group_x, 2));
        pos = (df > 0);
        if (any (pos))
          tval(pos) = - tinv (alpha / 2, df(pos));
        endif
        d = s .* tval;
        d(n < 2) = NaN;
        ci(idx,:) = d;
      endfor

      ## Plot the error bars
      h = errorbar ([1:ngroups]', mu, ci);
      xticks ([1:ngroups]);
      xlim ([0, ngroups+1]);

      ## For multiple grouping variables, the ticklabels need some work
      labels = g_names(:,1);
      if (columns (g_names) > 1)
        newlines = repmat ("\n", ngroups, 1);
        for ic = 2:columns (g_names)
          labels = strcat (labels, newlines, g_names(:,ic));
        endfor
      endif
      xticklabels (labels);

      ## MATLAB functional form does not return an output.  In Octave,
      ## if output is requested, we return an axes handle to the plot.
      if (nargout > 0)
        varargout{1} = h;
      endif
      return;
    endif

    ## Check consistent number of output arguments
    fcn_num = numel (fcn_names);
    if (! (nargout == 0 && fcn_num == 1) && nargout != fcn_num)
      error ("grpstats: inconsistent number of output arguments.");
    endif

    ## From this point we can start applying functions on the entire array
    for fcn_idx = 1:fcn_num
      if (is_function_handle (fcn_names{fcn_idx}))
        varargout{fcn_idx} = applyhandle (x, grp_idx, ngroups, c, ...
                                          fcn_names{fcn_idx});
        continue;
      endif
      switch (fcn_names{fcn_idx})
        case 'mean'
          group_mean = NaN (ngroups, c);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_mean(idx,:) = mean (group_x, 1, 'omitnan');
          endfor
          varargout{fcn_idx} = group_mean;
        case 'median'
          group_mean = NaN (ngroups, c);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_mean(idx,:) = median (group_x, 1, 'omitnan');
          endfor
          varargout{fcn_idx} = group_mean;
        case 'sem'
          group_sem = NaN (ngroups, c);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_sem(idx,:) = std (group_x, 0, 1, 'omitnan') ./ ...
                             sqrt (size (group_x, 1) - sum (isnan (group_x), 1));
          endfor
          varargout{fcn_idx} = group_sem;
        case 'std'
          group_std = NaN (ngroups, c);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_std(idx,:) = std (group_x, 0, 1, 'omitnan');
          endfor
          varargout{fcn_idx} = group_std;
        case 'var'
          group_var = NaN (ngroups, c);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_var(idx,:) = var (group_x, 0, 1, 'omitnan');
          endfor
          varargout{fcn_idx} = group_var;
        case 'min'
          group_min = NaN (ngroups, c);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_min(idx,:) = nanmin (group_x);
          endfor
          varargout{fcn_idx} = group_min;
        case 'max'
          group_max = NaN (ngroups, c);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_max(idx,:) = nanmax (group_x);
          endfor
          varargout{fcn_idx} = group_max;
        case 'range'
          func_handle = @(x) range (x, 1);
          group_range = NaN (ngroups, c);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_range(idx,:) = range (group_x, 1);
          endfor
          varargout{fcn_idx} = group_range;
        case 'numel'
          group_numel = NaN (ngroups, c);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_numel(idx,:) = size (group_x, 1) - sum (isnan (group_x), 1);
          endfor
          varargout{fcn_idx} = group_numel;
        case 'meanci'
          ## Allocate as 3-D: [ngroups x c x 2] (lower, upper)
          group_meanci = NaN (ngroups, c, 2);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            m = mean (group_x, 1, 'omitnan');
            n = size (group_x, 1) - sum (isnan (group_x), 1);
            s = std (group_x, 0, 1, 'omitnan') ./ sqrt (max (n,1));
            ## Avoid invalid tinv calls for degenerate df
            df = max (n - 1, 0);
            tval = zeros (1, size (group_x, 2));
            pos = (df > 0);
            if (any (pos))
              tval(pos) = - tinv (alpha / 2, df(pos));
            endif
            d = s .* tval;
            d(n < 2) = NaN;
            group_meanci(idx, :, 1) = m - d;
            group_meanci(idx, :, 2) = m + d;
          endfor
          ## MATLAB returns [ngroups x 2] when nvars == 1
          ## Octave used canonical 3-D.
          if (c == 1)
            ## Reshape to [ngroups x 2]
            varargout{fcn_idx} = reshape (group_meanci, ngroups, 2);
          else
            varargout{fcn_idx} = group_meanci;
          endif

        case 'predci'
          ## Allocate as 3-D: [ngroups x c x 2] (lower, upper)
          group_predci = NaN (ngroups, c, 2);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            m = mean (group_x, 1, 'omitnan');
            n = size (group_x, 1) - sum (isnan (group_x), 1);
            s = std (group_x, 0, 1, 'omitnan') .* sqrt (1 + (1 ./ max (n,1)));
            df = max (n - 1, 0);
            tval = zeros (1, size (group_x, 2));
            pos = (df > 0);
            if (any (pos))
              tval(pos) = - tinv (alpha / 2, df(pos));
            endif
            d = s .* tval;
            d(n < 2) = NaN;
            group_predci(idx, :, 1) = m - d;
            group_predci(idx, :, 2) = m + d;
          endfor
          if (c == 1)
            varargout{fcn_idx} = reshape (group_predci, ngroups, 2);
          else
            varargout{fcn_idx} = group_predci;
          endif
        case 'gname'
          varargout{fcn_idx} = g_names;
      endswitch
    endfor
  endif

endfunction

## A function handle in WHICHSTATS is applied to each column of each group, on
## the values of that column which are not NaN, and returns one value per call.
function out = applyhandle (x, grp_idx, ngroups, c, fcn_op)

  out = NaN (ngroups, c);
  for idx = 1:ngroups
    group_x = x(find (grp_idx == idx), :);
    for col_idx = 1:c
      col_x = group_x(:, col_idx);
      value = fcn_op (col_x(! isnan (col_x)));
      if (! isscalar (value))
        error (strcat ("grpstats: function in WHICHSTATS must return", ...
                       " a scalar for each column of each group."));
      endif
      out(idx, col_idx) = value;
    endfor
  endfor

endfunction

%!demo
%! load carsmall;
%! [m, p, g] = grpstats (Weight, Model_Year, {'mean', 'predci', 'gname'})
%! n = length (m);
%! errorbar ((1:n)',m,p(:,2)-m);
%! set (gca, 'xtick', 1:n, 'xticklabel', g);
%! title ('95% prediction intervals for mean weight by year');

%!demo
%! load carsmall;
%! [m, p, g] = grpstats ([Acceleration,Weight/1000],Cylinders, ...
%!                       {'mean', 'meanci', 'gname'}, 0.05)
%! [c, r] = size (m);
%! errorbar ((1:c)'.*ones (c,r),m,p(:,[(1:r)])-m);
%! set (gca, 'xtick', 1:c, 'xticklabel', g);
%! title ('95% prediction intervals for mean weight by year');

%!demo
%! ## Plot mean and 95% CI for a single grouping variable
%! load carsmall;
%! grpstats (Weight, Model_Year, 0.05);
%! title ('Mean Weight by Model Year');

%!demo
%! ## Plot mean and 95% CI for two grouping variables
%! load carsmall;
%! grpstats (Weight, {Origin, Cylinders}, 0.05);
%! title ('Mean Weight by Origin and Number of Cylinders');

%!test
%! load carsmall
%! means = grpstats (Acceleration, Origin);
%! assert_equal (means, [14.4377; 18.0500; 15.8867; 16.3778; 16.6000; 15.5000], 0.001);
%!test
%! load carsmall
%! [grpMin, grpMax, grp] = grpstats (Acceleration, Origin, {'min', 'max', ...
%!                                                          'gname'});
%! assert_equal (grpMin, [8.0; 15.3; 13.9; 12.2; 15.7; 15.5]);
%! assert_equal (grpMax, [22.2; 21.9; 18.2; 24.6; 17.5; 15.5]);
%!test
%! load carsmall
%! [grpMin, grpMax, grp] = grpstats (Acceleration, Origin, {'min', 'max', ...
%!                                                          'gname'});
%! assert_equal (grp', {'USA', 'France', 'Japan', 'Germany', 'Sweden', 'Italy'});
%!test
%! load carsmall
%! [m, p, g] = grpstats ([Acceleration, Weight/1000], Cylinders, ...
%!                       {'mean', 'meanci', 'gname'}, 0.05);
%! ## check meanci lower bounds (first slice) with tolerance
%! expected_lower = [15.9163; 15.6622; 10.7968];
%! expected_upper = [17.4249; 17.2907; 12.4845];
%! assert_equal (abs (p(:,1,1)), expected_lower, 1e-3);
%! assert_equal (abs (p(:,1,2)), expected_upper, 1e-3);
%!test
%! [mC, g] = grpstats ([], []);
%! assert_equal (isempty (mC), true);
%! assert_equal (isempty (g), true);
%!test
%! ## column vector, no group
%! x = [1; 2; 3; 4; 5];
%! m = grpstats (x);
%! expected = 3;
%! assert_equal (m, expected);
%!test
%! ## row vector, no group
%! x = [1 2 3 4 5];
%! m = grpstats (x);
%! expected = 3;
%! assert_equal (m, expected);
%!test
%! ## matrix, no group
%! x = [1 2; 3 4; 5 6];
%! m = grpstats (x);
%! expected = [3 4];
%! assert_equal (m, expected);
%!test
%! ## vector, numeric groups
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g);
%! expected = [15; 35; 55];
%! assert_equal (m, expected);
%!test
%! ## vector, cellstr groups
%! x = [10; 20; 30; 40; 50; 60];
%! g = {'A'; 'A'; 'B'; 'B'; 'C'; 'C'};
%! m = grpstats (x, g);
%! expected = [15; 35; 55];
%! assert_equal (m, expected);
%!test
%! ## matrix, numeric groups
%! x = [1 10; 2 20; 3 30; 4 40; 5 50; 6 60];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g);
%! expected = [1.5 15; 3.5 35; 5.5 55];
%! assert_equal (m, expected);
%!test
%! ## NaN handling
%! x = [1; NaN; 3; 4; NaN; 6];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g);
%! expected = [1; 3.5; 6];
%! assert_equal (m, expected);
%!test
%! ## single group
%! x = [1; 2; 3; 4; 5];
%! g = ones (5, 1);
%! m = grpstats (x, g);
%! expected = 3;
%! assert_equal (m, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g, 'mean');
%! expected = [15; 35; 55];
%! assert_equal (m, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g, 'median');
%! expected = [15; 35; 55];
%! assert_equal (m, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! s = grpstats (x, g, 'std');
%! expected = [7.07106781186548; 7.07106781186548; 7.07106781186548];
%! assert_equal (s, expected, 1e-14);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! v = grpstats (x, g, 'var');
%! expected = [50; 50; 50];
%! assert_equal (v, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! s = grpstats (x, g, 'sem');
%! expected = [5; 5; 5];
%! assert_equal (s, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! mn = grpstats (x, g, 'min');
%! expected = [10; 30; 50];
%! assert_equal (mn, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! mx = grpstats (x, g, 'max');
%! expected = [20; 40; 60];
%! assert_equal (mx, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! r = grpstats (x, g, 'range');
%! expected = [10; 10; 10];
%! assert_equal (r, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! n = grpstats (x, g, 'numel');
%! expected = [2; 2; 2];
%! assert_equal (n, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = {'A'; 'A'; 'B'; 'B'; 'C'; 'C'};
%! names = grpstats (x, g, 'gname');
%! expected = {'A'; 'B'; 'C'};
%! assert_equal (names, expected);
%!test
%! ## single statistic (default alpha)
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'meanci');
%! expected = [-48.5310236808735 78.5310236808735; -28.5310236808735 ...
%!             98.5310236808735; -8.53102368087348 118.531023680873];
%! assert_equal (ci, expected, 1e-12);
%!test
%! ## single statistic (default alpha)
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'predci');
%! expected = [-95.0389608721344 125.038960872134; -75.0389608721344 ...
%!             145.038960872134; -55.0389608721344 165.038960872134];
%! assert_equal (ci, expected, 1e-12);
%!test
%! ## mean + std
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [m, s] = grpstats (x, g, {'mean', 'std'});
%! expected_m = [15; 35; 55];
%! expected_s = [7.07106781186548; 7.07106781186548; 7.07106781186548];
%! assert_equal (m, expected_m, 1e-14);
%! assert_equal (s, expected_s, 1e-14);
%!test
%! ## min + max + range
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [mn, mx, r] = grpstats (x, g, {'min', 'max', 'range'});
%! expected_mn = [10; 30; 50];
%! expected_mx = [20; 40; 60];
%! expected_r  = [10; 10; 10];
%! assert_equal (mn, expected_mn);
%! assert_equal (mx, expected_mx);
%! assert_equal (r, expected_r);
%!test
%! ## mean + median + numel + gname
%! x = [10; 20; 30; 40; 50; 60];
%! g = {'A'; 'A'; 'B'; 'B'; 'C'; 'C'};
%! [m, med, n, names] = grpstats (x, g, {'mean', 'median', 'numel', 'gname'});
%! expected_m   = [15; 35; 55];
%! expected_med = [15; 35; 55];
%! expected_n   = [2; 2; 2];
%! expected_names = {'A'; 'B'; 'C'};
%! assert_equal (m, expected_m);
%! assert_equal (med, expected_med);
%! assert_equal (n, expected_n);
%! assert_equal (names, expected_names);
%!test
%! ## all basic statistics
%! x = [10; 20; 30; 40; 50; 60; 70; 80];
%! g = [1; 1; 2; 2; 2; 2; 3; 3];
%! [m, med, s, v, se, mn, mx, r, n] = grpstats (x, g, {'mean', 'median', ...
%!                                                     'std', 'var', 'sem', ...
%!                                                     'min', 'max', 'range', ...
%!                                                     'numel'});
%! expected_m   = [15; 45; 75];
%! expected_med = [15; 45; 75];
%! expected_s   = [7.07106781186548; 12.9099444873581; 7.07106781186548];
%! expected_v   = [50; 166.666666666667; 50];
%! expected_se  = [5; 6.45497224367903; 5];
%! expected_mn  = [10; 30; 70];
%! expected_mx  = [20; 60; 80];
%! expected_r   = [10; 30; 10];
%! expected_n   = [2; 4; 2];
%! assert_equal (m,   expected_m);
%! assert_equal (med, expected_med);
%! assert_equal (s,   expected_s, 1e-13);
%! assert_equal (v,   expected_v, 1e-12);
%! assert_equal (se,  expected_se, 1e-14);
%! assert_equal (mn,  expected_mn);
%! assert_equal (mx,  expected_mx);
%! assert_equal (r,   expected_r);
%! assert_equal (n,   expected_n);
%!test
%! ## meanci-alpha-0.1
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'meanci', 0.1);
%! expected = [-16.5687575733752 46.5687575733752; 3.4312424266248 ...
%!             66.5687575733752; 23.4312424266248 86.5687575733752];
%! assert_equal (ci, expected, 1e-13);
%!test
%! ## predci-alpha-0.1
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'predci', 0.1);
%! expected = [-39.6786920489106 69.6786920489106; -19.6786920489106 ...
%!             89.6786920489106; 0.321307951089366 109.678692048911];
%! assert_equal (ci, expected, 1e-12);
%!test
%! ## meanci-alpha-0.01
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'meanci', 0.01);
%! expected = [-303.283705814358 333.283705814358; -283.283705814358 ...
%!             353.283705814358; -263.283705814358 373.283705814358];
%! assert_equal (ci, expected, 3e-8);
%!test
%! ## predci-alpha-0.01
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'predci', 0.01);
%! expected = [-536.283549691775 566.283549691775; -516.283549691775 ...
%!             586.283549691775; -496.283549691775 606.283549691775];
%! assert_equal (ci, expected, 3e-8);
%!test
%! ## meanci-alpha-0.2
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'meanci', 0.2);
%! expected = [-0.388417685876263 30.3884176858763; 19.6115823141237 ...
%!             50.3884176858763; 39.6115823141237 70.3884176858763];
%! assert_equal (ci, expected, 1e-13);
%!test
%! ## predci-alpha-0.2
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'predci', 0.2);
%! expected = [-11.6535212800292 41.6535212800292; 8.34647871997083 ...
%!             61.6535212800292; 28.3464787199708 81.6535212800292];
%! assert_equal (ci, expected, 1e-13);
%!test
%! ## meanci, name-value alpha=0.2
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'meanci', 'alpha', 0.2);
%! expected = [-0.388417685876263 30.3884176858763; 19.6115823141237 ...
%!             50.3884176858763; 39.6115823141237 70.3884176858763];
%! assert_equal (ci, expected, 1e-13);
%!test
%! ## meanci + predci, alpha=0.01
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [ci_m, ci_p] = grpstats (x, g, {'meanci', 'predci'}, 0.01);
%! expected_m = [-303.283705814358 333.283705814358; -283.283705814358 ...
%!               353.283705814358; -263.283705814358 373.283705814358];
%! expected_p = [-536.283549691775 566.283549691775; -516.283549691775 ...
%!               586.283549691775; -496.283549691775 606.283549691775];
%! assert_equal (ci_m, expected_m, 3e-8);
%! assert_equal (ci_p, expected_p, 3e-8);
%!test
%! ## matrix, mean+std+numel
%! x = [1 10; 2 20; 3 30; 4 40; 5 50; 6 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [m, s, n] = grpstats (x, g, {'mean', 'std', 'numel'});
%! expected_m = [1.5 15; 3.5 35; 5.5 55];
%! expected_s = [0.707106781186548 7.07106781186548; 0.707106781186548 ...
%!               7.07106781186548; 0.707106781186548 7.07106781186548];
%! expected_n = [2 2; 2 2; 2 2];
%! assert_equal (m, expected_m);
%! assert_equal (s, expected_s, 1e-14);
%! assert_equal (n, expected_n);
%!test
%! ## matrix with NaN, mean+numel
%! x = [1 10; NaN 20; 3 NaN; 4 40; 5 50; 6 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [m, n] = grpstats (x, g, {'mean', 'numel'});
%! expected_m = [1   15; 3.5 40; 5.5 55];
%! expected_n = [1 2; 2 1; 2 2];
%! assert_equal (m, expected_m);
%! assert_equal (n, expected_n);
%!test
%! ## 3-column matrix, mean+min+max
%! x = [1 100 1000; 2 200 2000; 3 300 3000; 4 400 4000];
%! g = [1; 1; 2; 2];
%! [m, mn, mx] = grpstats (x, g, {'mean', 'min', 'max'});
%! expected_m  = [1.5 150 1500; 3.5 350 3500];
%! expected_mn = [1 100 1000; 3 300 3000];
%! expected_mx = [2 200 2000; 4 400 4000];
%! assert_equal (m,  expected_m);
%! assert_equal (mn, expected_mn);
%! assert_equal (mx, expected_mx);
%!test
%! ## one element per group
%! x = [1; 2; 3];
%! g = [1; 2; 3];
%! [m, s, n] = grpstats (x, g, {'mean', 'std', 'numel'});
%! expected_m = [1; 2; 3];
%! expected_s = [0; 0; 0];
%! expected_n = [1; 1; 1];
%! assert_equal (m, expected_m);
%! assert_equal (s, expected_s);
%! assert_equal (n, expected_n);
%!test
%! ## group with all NaN
%! x = [1; 2; NaN; NaN; 5; 6];
%! g = [1; 1; 2; 2; 3; 3];
%! [m, s, n] = grpstats (x, g, {'mean', 'std', 'numel'});
%! expected_m = [1.5; NaN; 5.5];
%! expected_s = [0.707106781186548; NaN; 0.707106781186548];
%! expected_n = [2; 0; 2];
%! assert_equal (m, expected_m);
%! assert_equal (s, expected_s, 1e-14);
%! assert_equal (n, expected_n);
%!test
%! ## unequal group sizes
%! x = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%! g = [1; 1; 1; 1; 2; 2; 2; 3; 3; 3];
%! [m, v, n] = grpstats (x, g, {'mean', 'var', 'numel'});
%! expected_m = [2.5; 6; 9];
%! expected_v = [1.66666666666667; 1; 1];
%! expected_n = [4; 3; 3];
%! assert_equal (m, expected_m);
%! assert_equal (v, expected_v, 1e-14);
%! assert_equal (n, expected_n);
%!test
%! ## non-consecutive numeric groups
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 5; 5; 10; 10];
%! [m, names] = grpstats (x, g, {'mean', 'gname'});
%! expected_m = [15; 35; 55];
%! expected_names = {'1'; '5'; '10'};
%! assert_equal (m, expected_m);
%! assert_equal (names, expected_names);
%!test
%! ## unsorted string groups
%! x = [30; 10; 40; 20; 60; 50];
%! g = {'C'; 'A'; 'C'; 'A'; 'B'; 'B'};
%! [m, names] = grpstats (x, g, {'mean', 'gname'});
%! expected_m = [35; 15; 55];
%! expected_names = {'C'; 'A'; 'B'};
%! assert_equal (m, expected_m);
%! assert_equal (names, expected_names);
%!test
%! ## 20 groups, one element each
%! x = (1:20)';
%! g = (1:20)';
%! [m, n] = grpstats (x, g, {'mean', 'numel'});
%! expected_m = (1:20)';
%! expected_n = ones (20, 1);
%! assert_equal (m, expected_m);
%! assert_equal (n, expected_n);
%!test
%! ## large sample meanci
%! x = (1:50)';
%! g = [ones(25, 1); 2 * ones(25, 1)];
%! ci = grpstats (x, g, 'meanci');
%! expected = [9.96202357522388 16.0379764247761; 34.9620235752239 ...
%!             41.0379764247761];
%! assert_equal (ci, expected, 1e-13);
%!test
%! ## large sample predci
%! x = (1:50)';
%! g = [ones(25, 1); 2 * ones(25, 1)];
%! ci = grpstats (x, g, 'predci');
%! expected = [-2.49070107176829 28.4907010717683; 22.5092989282317 ...
%!             53.4907010717683];
%! assert_equal (ci, expected, 2e-14);
%!test
%! Y = [5; 6; 7; 4; 9; 8];
%! X = [1; 2; 3; 4; 5; 6];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, X, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y', 'numel_Y', ...
%!                                              'mean_X', 'numel_X'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Y, [5.5; 5.5; 8.5]);
%! assert_equal (stats_tbl.numel_Y, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_X, [1.5; 3.5; 5.5]);
%! assert_equal (stats_tbl.numel_X, [2; 2; 2]);
%!test
%! Y = [5; 6; 7; 4; 9; 8];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Y, [5.5; 5.5; 8.5]);
%!test
%! Y = [10; 20; 30; 40];
%! X = [100; 200; 300; 400];
%! Z = [1000; 2000; 3000; 4000];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'});
%! tbl = table (Y, X, Z, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!         'mean_Y', 'numel_Y', 'mean_X', 'numel_X', 'mean_Z', 'numel_Z'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert_equal (stats_tbl.GroupCount, [2; 2]);
%! assert_equal (stats_tbl.mean_Y, [15; 35]);
%! assert_equal (stats_tbl.numel_Y, [2; 2]);
%! assert_equal (stats_tbl.mean_X, [150; 350]);
%! assert_equal (stats_tbl.numel_X, [2; 2]);
%! assert_equal (stats_tbl.mean_Z, [1500; 3500]);
%! assert_equal (stats_tbl.numel_Z, [2; 2]);
%!test
%! Y = [1; 2; 3; 4; 5; 6; 7; 8];
%! Group = categorical ({'A'; 'A'; 'A'; 'A'; 'B'; 'B'; 'B'; 'B'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert_equal (stats_tbl.GroupCount, [4; 4]);
%! assert_equal (stats_tbl.mean_Y, [2.5; 6.5]);
%!test
%! Y = [1; 2; 3; 4; 5; 6; 7];
%! Group = categorical ({'A'; 'A'; 'A'; 'A'; 'A'; 'B'; 'B'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y', 'numel_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert_equal (stats_tbl.GroupCount, [5; 2]);
%! assert_equal (stats_tbl.mean_Y, [3; 6.5]);
%! assert_equal (stats_tbl.numel_Y, [5; 2]);
%!test
%! Y = [10; 20; 30];
%! Group = categorical ({'A'; 'B'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [1; 1; 1]);
%! assert_equal (stats_tbl.mean_Y, [10; 20; 30]);
%!test
%! Y = [5; 5; 5; 5];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert_equal (stats_tbl.GroupCount, [2; 2]);
%! assert_equal (stats_tbl.mean_Y, [5; 5]);
%!test
%! Y = [1; NaN; 3; 4; NaN; 6];
%! X = [10; 20; NaN; 40; 50; NaN];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, X, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y', 'mean_X'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Y, [1; 3.5; 6]);
%! assert_equal (stats_tbl.mean_X, [15; 40; 50]);
%!test
%! Y = [1; NaN; 3; 4; 5; 6];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y', 'numel_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Y, [1; 3.5; 5.5]);
%! assert_equal (stats_tbl.numel_Y, [1; 2; 2]);
%!test
%! Y = [100; 200; 300; 400; 500; 600];
%! Group = categorical ({'Group1'; 'Group1'; 'Group2'; 'Group2'; ...
%!                       'Group3'; 'Group3'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'Group1'; 'Group2'; 'Group3'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Y, [150; 350; 550]);
%!test
%! Var1 = [1; 2; 3; 4];
%! Var2 = [10; 20; 30; 40];
%! Var3 = [100; 200; 300; 400];
%! Var4 = [1000; 2000; 3000; 4000];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'});
%! tbl = table (Var1, Var2, Var3, Var4, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!         'mean_Var1', 'mean_Var2', 'mean_Var3', 'mean_Var4'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert_equal (stats_tbl.GroupCount, [2; 2]);
%! assert_equal (stats_tbl.mean_Var1, [1.5; 3.5]);
%! assert_equal (stats_tbl.mean_Var2, [15; 35]);
%! assert_equal (stats_tbl.mean_Var3, [150; 350]);
%! assert_equal (stats_tbl.mean_Var4, [1500; 3500]);
%!test
%! Y = [1.5; 2.5; 3.5; 4.5; 5.5; 6.5];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Y, [2; 4; 6]);
%!test
%! Y = [-10; -20; 30; 40; 50; 60];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Y, [-15; 35; 55]);
%!test
%! Y = [0; 0; 0; 0; 0; 0];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Y, [0; 0; 0]);
%!test
%! Y = [1e6; 2e6; 3e6; 4e6; 5e6; 6e6];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Y, [1.5e6; 3.5e6; 5.5e6]);
%!test
%! Y = (1:10)';
%! Group = categorical (repmat ({'A'; 'B'}, 5, 1));
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, ...
%!   {'Group', 'GroupCount', 'mean_Y', 'numel_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert_equal (stats_tbl.GroupCount, [5; 5]);
%! assert_equal (stats_tbl.mean_Y, [5; 6]);
%! assert_equal (stats_tbl.numel_Y, [5; 5]);
%!test
%! Y = (1:20)';
%! Group = categorical (repmat ({'A'; 'B'; 'C'; 'D'}, 5, 1));
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'; 'D'});
%! assert_equal (stats_tbl.GroupCount, [5; 5; 5; 5]);
%! assert_equal (stats_tbl.mean_Y, [9; 10; 11; 12]);
%!test
%! Y = [1; 2; 3; 4; 5];
%! Group = categorical ({'A'; 'B'; 'C'; 'D'; 'E'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'; 'D'; 'E'});
%! assert_equal (stats_tbl.GroupCount, [1; 1; 1; 1; 1]);
%! assert_equal (stats_tbl.mean_Y, [1; 2; 3; 4; 5]);
%!test
%! Score1 = [85; 90; 78; 92; 88; 76];
%! Score2 = [82; 88; 75; 90; 85; 73];
%! Group = categorical ({'High'; 'High'; 'Med'; 'Med'; 'Low'; 'Low'});
%! tbl = table (Score1, Score2, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Score1', 'mean_Score2'});
%! assert_equal (stats_tbl.Properties.RowNames, {'High'; 'Low'; 'Med'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Score1, [87.5; 82; 85]);
%! assert_equal (stats_tbl.mean_Score2, [85; 79; 82.5]);
%!test
%! Height = [170; 175; 165; 180; 160; 185];
%! Weight = [70; 75; 65; 80; 60; 85];
%! Category = categorical ({'M'; 'M'; 'F'; 'F'; 'M'; 'M'});
%! tbl = table (Height, Weight, Category);
%! stats_tbl = grpstats (tbl, 'Category', {'mean', 'numel'});
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Category', 'GroupCount', ...
%!         'mean_Height', 'numel_Height', 'mean_Weight', 'numel_Weight'});
%! assert_equal (stats_tbl.Properties.RowNames, {'F'; 'M'});
%! assert_equal (stats_tbl.GroupCount, [2; 4]);
%! assert_equal (stats_tbl.mean_Height, [172.5; 172.5]);
%! assert_equal (stats_tbl.numel_Height, [2; 4]);
%! assert_equal (stats_tbl.mean_Weight, [72.5; 72.5]);
%! assert_equal (stats_tbl.numel_Weight, [2; 4]);
%!test
%! Value = [10.5; 11.2; 9.8; 10.1; 11.5; 10.8];
%! Type = categorical ({'A'; 'A'; 'A'; 'B'; 'B'; 'B'});
%! tbl = table (Value, Type);
%! stats_tbl = grpstats (tbl, 'Type', 'numel');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Type', 'GroupCount', ...
%!                                              'numel_Value'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert_equal (stats_tbl.GroupCount, [3; 3]);
%! assert_equal (stats_tbl.numel_Value, [3; 3]);
%!test
%! Data = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%! Label = categorical ({'A'; 'A'; 'A'; 'A'; 'A'; 'B'; 'B'; 'B'; 'B'; 'B'});
%! tbl = table (Data, Label);
%! stats_tbl = grpstats (tbl, 'Label', {'mean', 'numel'});
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Label', 'GroupCount', ...
%!                                              'mean_Data', 'numel_Data'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert_equal (stats_tbl.GroupCount, [5; 5]);
%! assert_equal (stats_tbl.mean_Data, [3; 8]);
%! assert_equal (stats_tbl.numel_Data, [5; 5]);
%!test
%! X1 = [1; 2; 3; 4];
%! X2 = [5; 6; 7; 8];
%! G = categorical ({'A'; 'A'; 'B'; 'B'});
%! tbl = table (X1, X2, G);
%! stats_tbl = grpstats (tbl, 'G', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'G', 'GroupCount', ...
%!                                              'mean_X1', 'mean_X2'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert_equal (stats_tbl.GroupCount, [2; 2]);
%! assert_equal (stats_tbl.mean_X1, [1.5; 3.5]);
%! assert_equal (stats_tbl.mean_X2, [5.5; 7.5]);
%!test
%! Measurement = [100; 150; 200; 250; 300; 350];
%! GroupVar = categorical ({'Control'; 'Control'; 'Treatment'; 'Treatment'; ...
%!                          'Placebo'; 'Placebo'});
%! tbl = table (Measurement, GroupVar);
%! stats_tbl = grpstats (tbl, 'GroupVar', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'GroupVar', 'GroupCount', ...
%!                                              'mean_Measurement'});
%! assert_equal (stats_tbl.Properties.RowNames, {'Control'; 'Placebo'; 'Treatment'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_Measurement, [125; 325; 225]);
%!test
%! Y = [NaN; NaN; 3; 4; 5; 6];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (isequaln (stats_tbl.mean_Y, [NaN; 3.5; 5.5]), true);
%!test
%! Y = [1; 2; NaN; NaN; NaN; NaN];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y', 'numel_Y'});
%! assert_equal (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (isequaln (stats_tbl.mean_Y, [1.5; NaN; NaN]), true);
%! assert_equal (stats_tbl.numel_Y, [2; 0; 0]);
%!test
%! Val = [5.5; 6.5; 7.5; 8.5];
%! Cat = categorical ({'X'; 'X'; 'Y'; 'Y'});
%! tbl = table (Val, Cat);
%! stats_tbl = grpstats (tbl, 'Cat', 'numel');
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Cat', 'GroupCount', ...
%!                                              'numel_Val'});
%! assert_equal (stats_tbl.Properties.RowNames, {'X'; 'Y'});
%! assert_equal (stats_tbl.GroupCount, [2; 2]);
%! assert_equal (stats_tbl.numel_Val, [2; 2]);
%!test
%! A = [1; 2; 3; 4; 5; 6];
%! B = [10; 20; 30; 40; 50; 60];
%! C = [100; 200; 300; 400; 500; 600];
%! Grp = categorical ({'G1'; 'G1'; 'G2'; 'G2'; 'G3'; 'G3'});
%! tbl = table (A, B, C, Grp);
%! stats_tbl = grpstats (tbl, 'Grp', {'mean', 'numel'});
%! assert_equal (istable (stats_tbl), true);
%! assert_equal (stats_tbl.Properties.VariableNames, {'Grp', 'GroupCount', ...
%!         'mean_A', 'numel_A', 'mean_B', 'numel_B', 'mean_C', 'numel_C'});
%! assert_equal (stats_tbl.Properties.RowNames, {'G1'; 'G2'; 'G3'});
%! assert_equal (stats_tbl.GroupCount, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_A, [1.5; 3.5; 5.5]);
%! assert_equal (stats_tbl.numel_A, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_B, [15; 35; 55]);
%! assert_equal (stats_tbl.numel_B, [2; 2; 2]);
%! assert_equal (stats_tbl.mean_C, [150; 350; 550]);
%! assert_equal (stats_tbl.numel_C, [2; 2; 2]);
%!test
%! x = [1; NaN; 3; 4];
%! g = [1; 1; 2; 2];
%! muci = grpstats (x, g, 'meanci');
%! assert_equal (muci, [NaN, NaN; -2.8531, 9.8531], 1e-4);
%!test
%! x = [1; NaN; 3; 4; 5; 6];
%! g = [1; 1; 1; 2; 2; 2];
%! predci = grpstats (x, g, 'predci');
%! assert_equal (predci, [-20.0078, 24.0078; 0.0317, 9.9683], 1e-4);

%!test
%! ## An empty X keeps its columns and has no groups
%! assert_equal (grpstats (zeros (0, 3), zeros (0, 1)), zeros (0, 3));
%! assert_equal (grpstats (zeros (0, 1), zeros (0, 1)), zeros (0, 1));
%! assert_equal (grpstats (zeros (0, 3)), zeros (0, 3));
%! assert_equal (grpstats ([], []), []);
%!test
%! ## Every reducing statistic keeps the columns of an empty X
%! g = zeros (0, 1);
%! assert_equal (grpstats (zeros (0, 3), g, 'sem'), zeros (0, 3));
%! assert_equal (grpstats (zeros (0, 3), g, 'std'), zeros (0, 3));
%! assert_equal (grpstats (zeros (0, 3), g, 'var'), zeros (0, 3));
%! assert_equal (grpstats (zeros (0, 3), g, 'min'), zeros (0, 3));
%! assert_equal (grpstats (zeros (0, 3), g, 'max'), zeros (0, 3));
%! assert_equal (grpstats (zeros (0, 3), g, 'range'), zeros (0, 3));
%! assert_equal (grpstats (zeros (0, 3), g, 'numel'), zeros (0, 3));
%!test
%! ## meanci and predci of an empty X keep the columns and the bound pair
%! g = zeros (0, 1);
%! assert_equal (size (grpstats (zeros (0, 3), g, 'meanci')), [0, 3, 2]);
%! assert_equal (size (grpstats (zeros (0, 3), g, 'predci')), [0, 3, 2]);
%! assert_equal (size (grpstats (zeros (0, 1), g, 'meanci')), [0, 2]);
%! assert_equal (size (grpstats ([], [], 'meanci')), [0, 0, 2]);
%!test
%! ## gname of an empty X is an empty cell column
%! assert_equal (grpstats (zeros (0, 3), zeros (0, 1), 'gname'), cell (0, 1));
%! assert_equal (grpstats ([], [], 'gname'), cell (0, 1));
%!test
%! ## Without a grouping variable X is a single group named '1'
%! assert_equal (grpstats ([1; 2; 3], [], 'gname'), {'1'});
%! assert_equal (grpstats ([1; 2; 3], [], 'mean'), 2);
%!test
%! ## A grouping variable that is entirely missing leaves no groups
%! assert_equal (grpstats ([1; 2; 3], [NaN; NaN; NaN]), zeros (0, 1));
%! assert_equal (grpstats ([1; 2; 3], [NaN; NaN; NaN], 'std'), zeros (0, 1));
%! assert_equal (grpstats ([1; 2; 3], [NaN; NaN; NaN], 'gname'), cell (0, 1));
%! assert_equal (grpstats ([1; 2; 3], {''; ''; ''}), zeros (0, 1));
%!test
%! ## A row missing in any of several grouping variables belongs to no group
%! y = [1; 2; 3; 4];
%! assert_equal (grpstats (y, {[NaN; 1; 2; 2], [1; 2; 1; 2]}), [2; 3; 4]);
%! gone = {[NaN; NaN; NaN; NaN], [1; 2; 1; 2]};
%! assert_equal (grpstats (y, gone), zeros (0, 1));
%!test
%! ## gname over several grouping variables has one column per variable
%! y = [1; 2; 3; 4];
%! gn = grpstats (y, {[NaN; NaN; NaN; NaN], [1; 2; 1; 2]}, 'gname');
%! assert_equal (gn, cell (0, 2));
%! gn = grpstats (y, {[1; 1; 2; 2], [1; 2; 1; 2]}, 'gname');
%! assert_equal (gn, {'1', '1'; '1', '2'; '2', '1'; '2', '2'});
%!test
%! ## Several empty grouping variables keep the columns of X
%! g = zeros (0, 1);
%! assert_equal (grpstats (zeros (0, 3), {g, g}), zeros (0, 3));
%! assert_equal (grpstats (zeros (0, 3), {g, g}, 'gname'), cell (0, 2));
%!test
%! ## Without WHICHSTATS the outputs are mean, sem, counts and names
%! [m, s, n, g] = grpstats ([1; 2; 3; 4], [1; 1; 2; 2]);
%! assert_equal (m, [1.5; 3.5]);
%! assert_equal (s, [0.5; 0.5]);
%! assert_equal (n, [2; 2]);
%! assert_equal (g, {'1'; '2'});
%!test
%! ## The default counts and standard error ignore NaNs
%! [m, s, n] = grpstats ([1; 2; NaN; 4], [1; 1; 2; 2]);
%! assert_equal (m, [1.5; 4]);
%! assert_equal (s, [0.5; 0]);
%! assert_equal (n, [2; 1]);
%!test
%! ## The default outputs are computed column by column
%! [m, s, n] = grpstats ([1, 10; 3, 20; 5, 30; 7, 40], [1; 1; 2; 2]);
%! assert_equal (m, [2, 15; 6, 35]);
%! assert_equal (s, [1, 5; 1, 5], 1e-14);
%! assert_equal (n, [2, 2; 2, 2]);
%!test
%! ## The default outputs of an empty X keep its columns
%! [m, s, n, g] = grpstats (zeros (0, 3), zeros (0, 1));
%! assert_equal (m, zeros (0, 3));
%! assert_equal (s, zeros (0, 3));
%! assert_equal (n, zeros (0, 3));
%! assert_equal (g, cell (0, 1));
%!test
%! ## The standard error is taken column by column
%! assert_equal (grpstats ([1, 10; 2, 20; 3, 30], [1; 1; 1], 'sem'), ...
%!               [0.5773502691896258, 5.773502691896258], 1e-14);
%!test
%! ## The standard error divides each column by its own count
%! assert_equal (grpstats ([1, 10; 2, NaN; 3, 30], [1; 1; 1], 'sem'), ...
%!               [0.5773502691896258, 10], 1e-14);
%!test
%! ## 'VarNames' renames every variable of the output table
%! Y = [1; 2; 3; 4];
%! Z = [10; 20; 30; 40];
%! G = {'a'; 'a'; 'b'; 'b'};
%! tbl = grpstats (table (Y, Z, G), 'G', 'mean', ...
%!                 'VarNames', {'grp', 'n', 'avgY', 'avgZ'});
%! assert_equal (tbl.Properties.VariableNames, {'grp', 'n', 'avgY', 'avgZ'});
%! assert_equal (tbl.n, [2; 2]);
%! assert_equal (tbl.avgY, [1.5; 3.5]);
%! assert_equal (tbl.avgZ, [15; 35]);
%!test
%! ## Without 'VarNames' the output table keeps the default names
%! Y = [1; 2; 3; 4];
%! G = {'a'; 'a'; 'b'; 'b'};
%! tbl = grpstats (table (Y, G), 'G', 'mean');
%! assert_equal (tbl.Properties.VariableNames, {'G', 'GroupCount', 'mean_Y'});
%!test
%! ## A function handle is applied to each column of each group
%! assert_equal (grpstats ([1; 2; 3; 4], [1; 1; 2; 2], @mean), [1.5; 3.5]);
%! x = [1, 10; 3, 20; 5, 30; 7, 40];
%! assert_equal (grpstats (x, [1; 1; 2; 2], @mean), [2, 15; 6, 35]);
%!test
%! ## A function handle sees only the values that are not NaN
%! assert_equal (grpstats ([1; 2; NaN; 4], [1; 1; 2; 2], @mean), [1.5; 4]);
%! assert_equal (grpstats ([1; 2; NaN; 4], [1; 1; 2; 2], @numel), [2; 1]);
%! x = [1, 10; 2, NaN; 3, 30; 4, 40];
%! assert_equal (grpstats (x, [1; 1; 1; 1], @mean), [2.5, 80/3], 1e-14);
%!test
%! ## A function handle may be mixed with named statistics
%! [m, s] = grpstats ([1; 2; 3; 4], [1; 1; 2; 2], {@mean, 'sem'});
%! assert_equal (m, [1.5; 3.5]);
%! assert_equal (s, [0.5; 0.5]);
%!test
%! ## An anonymous handle is applied like any other
%! fcn = @(v) max (v) - min (v);
%! assert_equal (grpstats ([1; 2; 3; 4], [1; 1; 2; 2], fcn), [1; 1]);
%! assert_equal (size (grpstats (zeros (0, 3), zeros (0, 1), @mean)), [0, 3]);
%!test
%! ## A table variable takes the handle by its own name, NaNs removed
%! Y = [1; 2; NaN; 4];
%! G = {'a'; 'a'; 'b'; 'b'};
%! tbl = grpstats (table (Y, G), 'G', @mean);
%! assert_equal (tbl.Properties.VariableNames, {'G', 'GroupCount', 'mean_Y'});
%! assert_equal (tbl.mean_Y, [1.5; 4]);
%! tbl = grpstats (table (Y, G), 'G', @numel);
%! assert_equal (tbl.numel_Y, [2; 1]);
%!test
%! ## A table with no grouping variable is a single group named 'All'
%! tbl = grpstats (table ([1; 2; 3], 'VariableNames', {'v'}));
%! assert_equal (tbl.Properties.VariableNames, {'GroupCount', 'mean_v'});
%! assert_equal (tbl.Properties.RowNames, {'All'});
%! assert_equal (tbl.GroupCount, 3);
%! assert_equal (tbl.mean_v, 2);
%!test
%! ## Every variable of an ungrouped table is a data variable
%! Y = [1; 2; 3; 4];
%! Z = [10; 20; 30; 40];
%! tbl = grpstats (table (Y, Z));
%! assert_equal (tbl.Properties.VariableNames, ...
%!               {'GroupCount', 'mean_Y', 'mean_Z'});
%! assert_equal (tbl.mean_Y, 2.5);
%! assert_equal (tbl.mean_Z, 25);
%!test
%! ## An ungrouped table takes WHICHSTATS and 'DataVars' as a grouped one does
%! Y = [1; 2; 3; 4];
%! Z = [10; 20; 30; 40];
%! tbl = grpstats (table (Y, Z), [], {'mean', 'sem'});
%! assert_equal (tbl.Properties.VariableNames, ...
%!               {'GroupCount', 'mean_Y', 'sem_Y', 'mean_Z', 'sem_Z'});
%! assert_equal (tbl.sem_Y, 0.6454972243679028, 1e-14);
%! tbl = grpstats (table (Y, Z), [], 'mean', 'DataVars', 'Y');
%! assert_equal (tbl.Properties.VariableNames, {'GroupCount', 'mean_Y'});
## Test input validation
%!error <grpstats: X must be a matrix or a table.> grpstats (ones (2, 2, 2))
%!error <grpstats: only one output argument in allowed when X is a table.> ...
%!       [a, b] = grpstats (table (1))
%!error <grpstats: invalid WHICHSTATS specification in cell array.> ...
%!       grpstats (ones (6, 2), [1; 1; 1; 2; 2; 2], {'mean', 1.5})
%!error <grpstats: invalid WHICHSTATS data type.> ...
%!       grpstats (ones (6, 2), [1; 1; 1; 2; 2; 2], 1.5)
%!error <grpstats: unrecognized function names in WHICHSTATS.> ...
%!       grpstats (ones (6, 2), [1; 1; 1; 2; 2; 2], 'some_function')
%!error <grpstats: unrecognized input arguments.> ...
%!       grpstats (ones (6, 2), [1; 1; 1; 2; 2; 2], 'mean', 35)
%!error <grpstats: unrecognized input arguments.> ...
%!       grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, 'predci', 'somename', -0.1);
%!error <grpstats: invalid data types for 'VarNames'.> ...
%!       grpstats (ones (6, 2), [1; 1; 1; 2; 2; 2], 'mean', 'VarNames', 3)
%!error <grpstats: X must be numeric to plot mean and CI for each group.> ...
%!       grpstats ({ones(6, 2)}, [], 0.05)
%!error <grpstats: 'alpha' must be a real scalar in the range \(0,1\).> ...
%!       grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, 'predci', 'alpha', -0.1);
%!error <grpstats: cannot resolve GROUPVARS in input table.> ...
%!       grpstats (table ([1:5]'), {'Var_5'})
%!error <grpstats: invalid 'DataVars' reference to table X.> ...
%!       grpstats (table ([1:5]'), {'Var1'}, [], 'DataVars', 'Var5')
%!error <grpstats: 'VarNames' do not match expected variables.> ...
%!       grpstats (table ([1:5]', [1:5]'), {'Var1'}, [], 'VarNames', {'A', 'B'})
%!error <grpstats: samples in X and GROUPS mismatch.> ...
%!       grpstats ([1:5]', {'A'; 'B'; 'A'; 'B'})
%!error <grpstats: inconsistent number of output arguments.> ...
%!       m = grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, {'mean', 'std'})
%!error <grpstats: too many output arguments.> ...
%!       [a, b, c, d, e] = grpstats ([1; 2; 3; 4], [1; 1; 2; 2])
%!error <grpstats: no groups to plot.> ...
%!       grpstats (zeros (0, 3), zeros (0, 1), 0.05)
%!error <grpstats: no groups to plot.> ...
%!       grpstats ([1; 2; 3], [NaN; NaN; NaN], 0.05)
%!error <grpstats: function in WHICHSTATS must return a scalar for each column of each group.> ...
%!       grpstats ([1; 2; 3; 4], [1; 1; 2; 2], @(v) [min(v), max(v)])
