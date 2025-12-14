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
## @deftypefn  {statistics} {@var{mean} =} grpstats (@var{x})
## @deftypefnx {statistics} {@var{mean} =} grpstats (@var{x}, @var{group})
## @deftypefnx {statistics} {[@var{a}, @var{b}, @dots{}] =} grpstats (@var{x}, @var{group}, @var{whichstats})
## @deftypefnx {statistics} {[@var{a}, @var{b}, @dots{}] =} grpstats (@var{x}, @var{group}, @var{whichstats}, @var{alpha})
## @deftypefnx {statistics} {[@var{a}, @var{b}, @dots{}] =} grpstats (@var{x}, @var{group}, @
## @var{whichstats}, @qcode{"alpha"}, @var{a})
##
## Summary statistics by group.
##
## @code{grpstats} computes groupwise summary statistics for the input data in
## @var{x}.  @code{grpstats} treats NaNs as missing values, and removes them.
##
## @code{@var{means} = grpstats (@var{x}, @var{group})}, when X is a matrix of
## observations, returns the means of each column of @var{x} by @var{group}.
## @var{group} is a grouping variable defined as a categorical variable,
## numeric, string array, or cell array of strings.  @var{group} can be [] or
## omitted to compute the mean of the entire sample without grouping.
##
## When the first input @var{x} is a table, @code{grpstats} computes groupwise
## summary statistics for the numeric variables in the table and returns the
## results in a new table.  In this case, the grouping variable @var{group}
## must be given as the name of a table variable, which is typically
## categorical.  Currently table input supports a subset of the statistics
## available for matrix input, such as @qcode{"mean"} and @qcode{"numel"}.
##
## @code{[@var{a}, @var{b}, @dots{}] = grpstats (@var{x}, @var{group},
## @var{whichstats})}, for a numeric matrix X, returns the statistics specified
## by @var{whichstats}, as separate arrays @var{a}, @var{b}, @dots{}.
## @var{whichstats} can be a single function name, or a cell array containing
## multiple function names.  The number of output arguments must match the
## number of function names in @var{whichstats}.
## Names in @var{whichstats} can be chosen from among the following:
##
## @multitable @columnfractions 0.05 0.2 0.75
## @item @tab "mean" @tab mean
## @item @tab "median" @tab median
## @item @tab "sem" @tab standard error of the mean
## @item @tab "std" @tab standard deviation
## @item @tab "var" @tab variance
## @item @tab "min" @tab minimum value
## @item @tab "max" @tab maximum value
## @item @tab "range" @tab maximum - minimum
## @item @tab "numel" @tab count, or number of elements
## @item @tab "meanci" @tab 95% confidence interval for the mean
## @item @tab "predci" @tab 95% prediction interval for a new observation
## @item @tab "gname" @tab group name
## @end multitable
##
## @code{[@dots{}] = grpstats (@var{x}, @var{group}, @var{whichstats},
## @var{alpha})} specifies the confidence level as 100(1-ALPHA)% for the "meanci"
## and "predci" options.  Default value for @var{alpha} is 0.05.  The
## significance can also be specified using the Name-Value pair argument syntax.
##
## Examples:
##
## @example
## load carsmall;
## [m,p,g] = grpstats (Weight, Model_Year, @{"mean", "predci", "gname"@})
## n = length(m);
## errorbar((1:n)',m,p(:,2)-m)
## set (gca, "xtick", 1:n, "xticklabel", g);
## title ("95% prediction intervals for mean weight by year")
## @end example
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
  if (isempty (x))
    [varargout] = repmat ({[]}, nargout, 1);
    return;
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
    r_names = {'All'};
    ngroups = 1;
    no_group = false;
  endif

  ## Check for plotting functional form with three input arguments
  if (nargin == 3 && isscalar (whichstats) && whichstats > 0 && whichstats < 1)
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
    fcn_names = {'mean'};
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
  [alpha, DataVars, VarNames, args] = pairedArgs (optNames, dfValues, ...
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
          fname = cellstr (fcn_op);
          new_vname = strcat (fname, '_', vname);
          try
            for idx = 1:ngroups
              new_vdata(idx,:) = fcn_op (vdata(grp_idx == idx, :));
            endfor
          catch
            # no-op
          end_try_catch
          if (! isempty (new_vdata))
            g_names = addvars (g_names, new_vdata, 'NewVariableNames', new_vname);
          endif
        else  # it must be a function name
          switch (fcn_op)
            case "mean"
              new_vname = strcat ('mean_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = mean (gdata, 1, "omitnan");
              endfor
            case "median"
              new_vname = strcat ('median_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = median (gdata, 1, "omitnan");
              endfor
            case "sem"
              new_vname = strcat ('sem_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = std (gdata, 0, 1, "omitnan") ./ ...
                                   sqrt (size (gdata, 1) - sum (isnan (gdata), 1));
              endfor
            case "std"
              new_vname = strcat ('std_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = std (gdata, 0, 1, "omitnan");
              endfor
            case "var"
              new_vname = strcat ('var_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = var (gdata, 0, 1, "omitnan");
              endfor
            case "min"
              new_vname = strcat ('min_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = min (gdata, [], 1);
              endfor
            case "max"
              new_vname = strcat ('max_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = max (gdata, [], 1);
              endfor
            case "range"
              new_vname = strcat ('range_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = range (gdata, 1);
              endfor
            case "numel"
              new_vname = strcat ('numel_', vname);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                new_vdata(idx,:) = size (gdata, 1) - sum (isnan (gdata), 1);
              endfor
            case "meanci"
              new_vname = strcat ('meanci_', vname);
              ## Preallocate twice the columns in variable (lower, upper)
              new_vdata = NaN (ngroups, ncols * 2);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                m = mean (gdata, 1, "omitnan");
                n = size (gdata, 1) - sum (isnan (gdata), 1);
                s = std (gdata, 0, 1, "omitnan") ./ sqrt (max (n,1));
                ## Avoid invalid tinv calls for degenerate df
                df = max (n - 1, 0);
                tval = zeros (1, size (gdata, 2));
                pos = (df > 0);
                if (any (pos))
                  tval(pos) = - tinv (alpha / 2, df(pos));
                endif
                d = s .* tval;
                new_vdata(idx,[1:2:end]) = m - d;
                new_vdata(idx,[2:2:end]) = m + d;
              endfor

            case "predci"
              new_vname = strcat ('predci_', vname);
              ## Preallocate twice the columns in variable (lower, upper)
              new_vdata = NaN (ngroups, ncols * 2);
              for idx = 1:ngroups
                gdata = vdata(grp_idx == idx, :);
                m = mean (gdata, 1, "omitnan");
                n = size (gdata, 1) - sum (isnan (gdata), 1);
                s = std (gdata, 0, 1, "omitnan") .* sqrt (1 + (1 ./ max (n,1)));
                df = max (n - 1, 0);
                tval = zeros (1, size (gdata, 2));
                pos = (df > 0);
                if (any (pos))
                  tval(pos) = - tinv (alpha / 2, df(pos));
                endif
                d = s .* tval;
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
      g_names - renamevars (g_names, VarNames);
    endif

    varargout{1} = g_names;
    return;
  else
    ## Get groups for array input data
    if (no_group)
      if (iscell (group) && ! iscellstr (group))
        ## Multilple grouping variables in cell array
        [grp_idx, g_names] = grp2idx (group{1});
        if (numel (grp_idx) != r)
          error ("grpstats: samples in X and GROUPS mismatch.");
        endif
        grp_vars = numel (group);
        if (grp_vars > 1)
          for g_idx = 2:grp_vars
            [tmp_grp_idx, tmp_g_names] = grp2idx (group{g_idx});
            if (numel (tmp_grp_idx) != r)
              error ("grpstats: samples in X and GROUPS mismatch.");
            endif
            grp_idx = [grp_idx, tmp_grp_idx];
            g_names = [g_names, tmp_g_names];
          endfor
        endif
        ## Get combination of unique groups and thei common index to X
        [g_names_idx, ~, grp_idx] = unique (grp_idx, 'rows');
        g_names = g_names(g_names_idx);
        ngroups = rows (g_names);
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
      ## MATLAB functional form does not return an output.  In Octave,
      ## if output is requested, we return an axes handle to the plot.

      ## Calculate mean and ci
      for idx = 1:ngroups
        group_x = x(find (grp_idx == idx), :);
        mu(idx,:) = mean (group_x, 1, "omitnan");
        n = size (group_x, 1) - sum (isnan (group_x), 1);
        s = std (group_x, 0, 1, "omitnan") ./ sqrt (max (n,1));
        ## Avoid invalid tinv calls for degenerate df
        df = max (n - 1, 0);
        tval = zeros (1, size (group_x, 2));
        pos = (df > 0);
        if (any (pos))
          tval(pos) = - tinv (alpha / 2, df(pos));
        endif
        ci(idx,:) = s .* tval;
      endfor



      return;
    endif

    ## Check consistent number of output arguments
    fcn_num = numel (fcn_names);
    if (! (nargout == 0 && fcn_num == 1) && nargout != fcn_num)
      error ("grpstats: inconsistent number of output arguments.");
    endif

    ## From this point we can start applying functions on the entire array
    for fcn_idx = 1:fcn_num
      switch (fcn_names{fcn_idx})
        case "mean"
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_mean(idx,:) = mean (group_x, 1, "omitnan");
          endfor
          varargout{fcn_idx} = group_mean;
        case "median"
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_mean(idx,:) = median (group_x, 1, "omitnan");
          endfor
          varargout{fcn_idx} = group_mean;
        case "sem"
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_sem(idx,:) = std (group_x, 0, 1, "omitnan") / ...
                             sqrt (size (group_x, 1) - sum (isnan (group_x), 1));
          endfor
          varargout{fcn_idx} = group_sem;
        case "std"
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_std(idx,:) = std (group_x, 0, 1, "omitnan");
          endfor
          varargout{fcn_idx} = group_std;
        case "var"
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_var(idx,:) = var (group_x, 0, 1, "omitnan");
          endfor
          varargout{fcn_idx} = group_var;
        case "min"
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_min(idx,:) = nanmin (group_x);
          endfor
          varargout{fcn_idx} = group_min;
        case "max"
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_max(idx,:) = nanmax (group_x);
          endfor
          varargout{fcn_idx} = group_max;
        case "range"
          func_handle = @(x) range (x, 1);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_range(idx,:) = range (group_x, 1);
          endfor
          varargout{fcn_idx} = group_range;
        case "numel"
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            group_numel(idx,:) = size (group_x, 1) - sum (isnan (group_x), 1);
          endfor
          varargout{fcn_idx} = group_numel;
        case "meanci"
          ## Allocate as 3-D: [ngroups x c x 2] (lower, upper)
          group_meanci = NaN (ngroups, c, 2);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            m = mean (group_x, 1, "omitnan");
            n = size (group_x, 1) - sum (isnan (group_x), 1);
            s = std (group_x, 0, 1, "omitnan") ./ sqrt (max (n,1));
            ## Avoid invalid tinv calls for degenerate df
            df = max (n - 1, 0);
            tval = zeros (1, size (group_x, 2));
            pos = (df > 0);
            if (any (pos))
              tval(pos) = - tinv (alpha / 2, df(pos));
            endif
            d = s .* tval;
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

        case "predci"
          ## Allocate as 3-D: [ngroups x c x 2] (lower, upper)
          group_predci = NaN (ngroups, c, 2);
          for idx = 1:ngroups
            group_x = x(find (grp_idx == idx), :);
            m = mean (group_x, 1, "omitnan");
            n = size (group_x, 1) - sum (isnan (group_x), 1);
            s = std (group_x, 0, 1, "omitnan") .* sqrt (1 + (1 ./ max (n,1)));
            df = max (n - 1, 0);
            tval = zeros (1, size (group_x, 2));
            pos = (df > 0);
            if (any (pos))
              tval(pos) = - tinv (alpha / 2, df(pos));
            endif
            d = s .* tval;
            group_predci(idx, :, 1) = m - d;
            group_predci(idx, :, 2) = m + d;
          endfor
          if (c == 1)
            varargout{fcn_idx} = reshape (group_predci, ngroups, 2);
          else
            varargout{fcn_idx} = group_predci;
          endif
        case "gname"
          varargout{fcn_idx} = g_names;
      endswitch
    endfor
  endif

endfunction

%!demo
%! load carsmall;
%! [m, p, g] = grpstats (Weight, Model_Year, {'mean', 'predci', 'gname'})
%! n = length (m);
%! errorbar ((1:n)',m,p(:,2)-m);
%! set (gca, "xtick", 1:n, "xticklabel", g);
%! title ("95% prediction intervals for mean weight by year");

%!demo
%! load carsmall;
%! [m, p, g] = grpstats ([Acceleration,Weight/1000],Cylinders, ...
%!                       {'mean', 'meanci', 'gname'}, 0.05)
%! [c, r] = size (m);
%! errorbar ((1:c)'.*ones(c,r),m,p(:,[(1:r)])-m);
%! set (gca, "xtick", 1:c, "xticklabel", g);
%! title ("95% prediction intervals for mean weight by year");

%!test
%! load carsmall
%! means = grpstats (Acceleration, Origin);
%! assert (means, [14.4377; 18.0500; 15.8867; 16.3778; 16.6000; 15.5000], 0.001);
%!test
%! load carsmall
%! [grpMin, grpMax, grp] = grpstats (Acceleration, Origin, {'min', 'max', ...
%!                                                          'gname'});
%! assert (grpMin, [8.0; 15.3; 13.9; 12.2; 15.7; 15.5]);
%! assert (grpMax, [22.2; 21.9; 18.2; 24.6; 17.5; 15.5]);
%!test
%! load carsmall
%! [grpMin, grpMax, grp] = grpstats (Acceleration, Origin, {'min', 'max', ...
%!                                                          'gname'});
%! assert (grp', {'USA', 'France', 'Japan', 'Germany', 'Sweden', 'Italy'});
%!test
%! load carsmall
%! [m, p, g] = grpstats ([Acceleration, Weight/1000], Cylinders, ...
%!                       {'mean', 'meanci', 'gname'}, 0.05);
%! ## check meanci lower bounds (first slice) with tolerance
%! expected_lower = [15.9163; 15.6622; 10.7968];
%! expected_upper = [17.4249; 17.2907; 12.4845];
%! assert (abs (p(:,1,1)), expected_lower, 1e-3);
%! assert (abs (p(:,1,2)), expected_upper, 1e-3);
%!test
%! [mC, g] = grpstats ([], []);
%! assert (isempty (mC), true);
%! assert (isempty (g), true);
%!test
%! ## column vector, no group
%! x = [1; 2; 3; 4; 5];
%! m = grpstats (x);
%! expected = 3;
%! assert (m, expected);
%!test
%! ## row vector, no group
%! x = [1 2 3 4 5];
%! m = grpstats (x);
%! expected = 3;
%! assert (m, expected);
%!test
%! ## matrix, no group
%! x = [1 2; 3 4; 5 6];
%! m = grpstats (x);
%! expected = [3 4];
%! assert (m, expected);
%!test
%! ## vector, numeric groups
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g);
%! expected = [15; 35; 55];
%! assert (m, expected);
%!test
%! ## vector, cellstr groups
%! x = [10; 20; 30; 40; 50; 60];
%! g = {'A'; 'A'; 'B'; 'B'; 'C'; 'C'};
%! m = grpstats (x, g);
%! expected = [15; 35; 55];
%! assert (m, expected);
%!test
%! ## matrix, numeric groups
%! x = [1 10; 2 20; 3 30; 4 40; 5 50; 6 60];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g);
%! expected = [1.5 15; 3.5 35; 5.5 55];
%! assert (m, expected);
%!test
%! ## NaN handling
%! x = [1; NaN; 3; 4; NaN; 6];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g);
%! expected = [1; 3.5; 6];
%! assert (m, expected);
%!test
%! ## single group
%! x = [1; 2; 3; 4; 5];
%! g = ones (5, 1);
%! m = grpstats (x, g);
%! expected = 3;
%! assert (m, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g, 'mean');
%! expected = [15; 35; 55];
%! assert (m, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g, 'median');
%! expected = [15; 35; 55];
%! assert (m, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! s = grpstats (x, g, 'std');
%! expected = [7.07106781186548; 7.07106781186548; 7.07106781186548];
%! assert (s, expected, 1e-14);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! v = grpstats (x, g, 'var');
%! expected = [50; 50; 50];
%! assert (v, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! s = grpstats (x, g, 'sem');
%! expected = [5; 5; 5];
%! assert (s, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! mn = grpstats (x, g, 'min');
%! expected = [10; 30; 50];
%! assert (mn, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! mx = grpstats (x, g, 'max');
%! expected = [20; 40; 60];
%! assert (mx, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! r = grpstats (x, g, 'range');
%! expected = [10; 10; 10];
%! assert (r, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! n = grpstats (x, g, 'numel');
%! expected = [2; 2; 2];
%! assert (n, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = {"A"; "A"; "B"; "B"; "C"; "C"};
%! names = grpstats (x, g, 'gname');
%! expected = {"A"; "B"; "C"};
%! assert (names, expected);
%!test
%! ## single statistic (default alpha)
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'meanci');
%! expected = [-48.5310236808735 78.5310236808735; -28.5310236808735 ...
%!             98.5310236808735; -8.53102368087348 118.531023680873];
%! assert (ci, expected, 1e-12);
%!test
%! ## single statistic (default alpha)
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'predci');
%! expected = [-95.0389608721344 125.038960872134; -75.0389608721344 ...
%!             145.038960872134; -55.0389608721344 165.038960872134];
%! assert (ci, expected, 1e-12);
%!test
%! ## mean + std
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [m, s] = grpstats (x, g, {'mean', 'std'});
%! expected_m = [15; 35; 55];
%! expected_s = [7.07106781186548; 7.07106781186548; 7.07106781186548];
%! assert (m, expected_m, 1e-14);
%! assert (s, expected_s, 1e-14);
%!test
%! ## min + max + range
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [mn, mx, r] = grpstats (x, g, {'min', 'max', 'range'});
%! expected_mn = [10; 30; 50];
%! expected_mx = [20; 40; 60];
%! expected_r  = [10; 10; 10];
%! assert (mn, expected_mn);
%! assert (mx, expected_mx);
%! assert (r, expected_r);
%!test
%! ## mean + median + numel + gname
%! x = [10; 20; 30; 40; 50; 60];
%! g = {'A'; 'A'; 'B'; 'B'; 'C'; 'C'};
%! [m, med, n, names] = grpstats (x, g, {'mean', 'median', 'numel', 'gname'});
%! expected_m   = [15; 35; 55];
%! expected_med = [15; 35; 55];
%! expected_n   = [2; 2; 2];
%! expected_names = {'A'; 'B'; 'C'};
%! assert (m, expected_m);
%! assert (med, expected_med);
%! assert (n, expected_n);
%! assert (names, expected_names);
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
%! assert (m,   expected_m);
%! assert (med, expected_med);
%! assert (s,   expected_s, 1e-13);
%! assert (v,   expected_v, 1e-12);
%! assert (se,  expected_se, 1e-14);
%! assert (mn,  expected_mn);
%! assert (mx,  expected_mx);
%! assert (r,   expected_r);
%! assert (n,   expected_n);
%!test
%! ## meanci-alpha-0.1
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'meanci', 0.1);
%! expected = [-16.5687575733752 46.5687575733752; 3.4312424266248 ...
%!             66.5687575733752; 23.4312424266248 86.5687575733752];
%! assert (ci, expected, 1e-13);
%!test
%! ## predci-alpha-0.1
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'predci', 0.1);
%! expected = [-39.6786920489106 69.6786920489106; -19.6786920489106 ...
%!             89.6786920489106; 0.321307951089366 109.678692048911];
%! assert (ci, expected, 1e-12);
%!test
%! ## meanci-alpha-0.01
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'meanci', 0.01);
%! expected = [-303.283705814358 333.283705814358; -283.283705814358 ...
%!             353.283705814358; -263.283705814358 373.283705814358];
%! assert (ci, expected, 3e-8);
%!test
%! ## predci-alpha-0.01
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'predci', 0.01);
%! expected = [-536.283549691775 566.283549691775; -516.283549691775 ...
%!             586.283549691775; -496.283549691775 606.283549691775];
%! assert (ci, expected, 3e-8);
%!test
%! ## meanci-alpha-0.2
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'meanci', 0.2);
%! expected = [-0.388417685876263 30.3884176858763; 19.6115823141237 ...
%!             50.3884176858763; 39.6115823141237 70.3884176858763];
%! assert (ci, expected, 1e-13);
%!test
%! ## predci-alpha-0.2
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'predci', 0.2);
%! expected = [-11.6535212800292 41.6535212800292; 8.34647871997083 ...
%!             61.6535212800292; 28.3464787199708 81.6535212800292];
%! assert (ci, expected, 1e-13);
%!test
%! ## meanci, name-value alpha=0.2
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, 'meanci', 'alpha', 0.2);
%! expected = [-0.388417685876263 30.3884176858763; 19.6115823141237 ...
%!             50.3884176858763; 39.6115823141237 70.3884176858763];
%! assert (ci, expected, 1e-13);
%!test
%! ## meanci + predci, alpha=0.01
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [ci_m, ci_p] = grpstats (x, g, {'meanci', 'predci'}, 0.01);
%! expected_m = [-303.283705814358 333.283705814358; -283.283705814358 ...
%!               353.283705814358; -263.283705814358 373.283705814358];
%! expected_p = [-536.283549691775 566.283549691775; -516.283549691775 ...
%!               586.283549691775; -496.283549691775 606.283549691775];
%! assert (ci_m, expected_m, 3e-8);
%! assert (ci_p, expected_p, 3e-8);
%!test
%! ## matrix, mean+std+numel
%! x = [1 10; 2 20; 3 30; 4 40; 5 50; 6 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [m, s, n] = grpstats (x, g, {'mean', 'std', 'numel'});
%! expected_m = [1.5 15; 3.5 35; 5.5 55];
%! expected_s = [0.707106781186548 7.07106781186548; 0.707106781186548 ...
%!               7.07106781186548; 0.707106781186548 7.07106781186548];
%! expected_n = [2 2; 2 2; 2 2];
%! assert (m, expected_m);
%! assert (s, expected_s, 1e-14);
%! assert (n, expected_n);
%!test
%! ## matrix with NaN, mean+numel
%! x = [1 10; NaN 20; 3 NaN; 4 40; 5 50; 6 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [m, n] = grpstats (x, g, {'mean', 'numel'});
%! expected_m = [1   15; 3.5 40; 5.5 55];
%! expected_n = [1 2; 2 1; 2 2];
%! assert (m, expected_m);
%! assert (n, expected_n);
%!test
%! ## 3-column matrix, mean+min+max
%! x = [1 100 1000; 2 200 2000; 3 300 3000; 4 400 4000];
%! g = [1; 1; 2; 2];
%! [m, mn, mx] = grpstats (x, g, {'mean', 'min', 'max'});
%! expected_m  = [1.5 150 1500; 3.5 350 3500];
%! expected_mn = [1 100 1000; 3 300 3000];
%! expected_mx = [2 200 2000; 4 400 4000];
%! assert (m,  expected_m);
%! assert (mn, expected_mn);
%! assert (mx, expected_mx);
%!test
%! ## one element per group
%! x = [1; 2; 3];
%! g = [1; 2; 3];
%! [m, s, n] = grpstats (x, g, {'mean', 'std', 'numel'});
%! expected_m = [1; 2; 3];
%! expected_s = [0; 0; 0];
%! expected_n = [1; 1; 1];
%! assert (m, expected_m);
%! assert (s, expected_s);
%! assert (n, expected_n);
%!test
%! ## group with all NaN
%! x = [1; 2; NaN; NaN; 5; 6];
%! g = [1; 1; 2; 2; 3; 3];
%! [m, s, n] = grpstats (x, g, {'mean', 'std', 'numel'});
%! expected_m = [1.5; NaN; 5.5];
%! expected_s = [0.707106781186548; NaN; 0.707106781186548];
%! expected_n = [2; 0; 2];
%! assert (m, expected_m);
%! assert (s, expected_s, 1e-14);
%! assert (n, expected_n);
%!test
%! ## unequal group sizes
%! x = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%! g = [1; 1; 1; 1; 2; 2; 2; 3; 3; 3];
%! [m, v, n] = grpstats (x, g, {'mean', 'var', "numel"});
%! expected_m = [2.5; 6; 9];
%! expected_v = [1.66666666666667; 1; 1];
%! expected_n = [4; 3; 3];
%! assert (m, expected_m);
%! assert (v, expected_v, 1e-14);
%! assert (n, expected_n);
%!test
%! ## non-consecutive numeric groups
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 5; 5; 10; 10];
%! [m, names] = grpstats (x, g, {'mean', 'gname'});
%! expected_m = [15; 35; 55];
%! expected_names = {'1'; '5'; '10'};
%! assert (m, expected_m);
%! assert (names, expected_names);
%!test
%! ## unsorted string groups
%! x = [30; 10; 40; 20; 60; 50];
%! g = {'C'; 'A'; 'C'; 'A'; 'B'; 'B'};
%! [m, names] = grpstats (x, g, {"mean", "gname"});
%! expected_m = [35; 15; 55];
%! expected_names = {'C'; 'A'; 'B'};
%! assert (m, expected_m);
%! assert (names, expected_names);
%!test
%! ## 20 groups, one element each
%! x = (1:20)';
%! g = (1:20)';
%! [m, n] = grpstats (x, g, {'mean', 'numel'});
%! expected_m = (1:20)';
%! expected_n = ones (20, 1);
%! assert (m, expected_m);
%! assert (n, expected_n);
%!test
%! ## large sample meanci
%! x = (1:50)';
%! g = [ones(25, 1); 2 * ones(25, 1)];
%! ci = grpstats (x, g, 'meanci');
%! expected = [9.96202357522388 16.0379764247761; 34.9620235752239 ...
%!             41.0379764247761];
%! assert (ci, expected, 1e-13);
%!test
%! ## large sample predci
%! x = (1:50)';
%! g = [ones(25, 1); 2 * ones(25, 1)];
%! ci = grpstats (x, g, 'predci');
%! expected = [-2.49070107176829 28.4907010717683; 22.5092989282317 ...
%!             53.4907010717683];
%! assert (ci, expected, 1e-14);
%!test
%! Y = [5; 6; 7; 4; 9; 8];
%! X = [1; 2; 3; 4; 5; 6];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, X, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y', 'numel_Y', ...
%!                                              'mean_X', 'numel_X'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Y, [5.5; 5.5; 8.5]);
%! assert (stats_tbl.numel_Y, [2; 2; 2]);
%! assert (stats_tbl.mean_X, [1.5; 3.5; 5.5]);
%! assert (stats_tbl.numel_X, [2; 2; 2]);
%!test
%! Y = [5; 6; 7; 4; 9; 8];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Y, [5.5; 5.5; 8.5]);
%!test
%! Y = [10; 20; 30; 40];
%! X = [100; 200; 300; 400];
%! Z = [1000; 2000; 3000; 4000];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'});
%! tbl = table (Y, X, Z, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!         'mean_Y', 'numel_Y', 'mean_X', 'numel_X', 'mean_Z', 'numel_Z'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert (stats_tbl.GroupCount, [2; 2]);
%! assert (stats_tbl.mean_Y, [15; 35]);
%! assert (stats_tbl.numel_Y, [2; 2]);
%! assert (stats_tbl.mean_X, [150; 350]);
%! assert (stats_tbl.numel_X, [2; 2]);
%! assert (stats_tbl.mean_Z, [1500; 3500]);
%! assert (stats_tbl.numel_Z, [2; 2]);
%!test
%! Y = [1; 2; 3; 4; 5; 6; 7; 8];
%! Group = categorical ({'A'; 'A'; 'A'; 'A'; 'B'; 'B'; 'B'; 'B'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert (stats_tbl.GroupCount, [4; 4]);
%! assert (stats_tbl.mean_Y, [2.5; 6.5]);
%!test
%! Y = [1; 2; 3; 4; 5; 6; 7];
%! Group = categorical ({'A'; 'A'; 'A'; 'A'; 'A'; 'B'; 'B'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y', 'numel_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert (stats_tbl.GroupCount, [5; 2]);
%! assert (stats_tbl.mean_Y, [3; 6.5]);
%! assert (stats_tbl.numel_Y, [5; 2]);
%!test
%! Y = [10; 20; 30];
%! Group = categorical ({'A'; 'B'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [1; 1; 1]);
%! assert (stats_tbl.mean_Y, [10; 20; 30]);
%!test
%! Y = [5; 5; 5; 5];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert (stats_tbl.GroupCount, [2; 2]);
%! assert (stats_tbl.mean_Y, [5; 5]);
%!test
%! Y = [1; NaN; 3; 4; NaN; 6];
%! X = [10; 20; NaN; 40; 50; NaN];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, X, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y', 'mean_X'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Y, [1; 3.5; 6]);
%! assert (stats_tbl.mean_X, [15; 40; 50]);
%!test
%! Y = [1; NaN; 3; 4; 5; 6];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y', 'numel_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Y, [1; 3.5; 5.5]);
%! assert (stats_tbl.numel_Y, [1; 2; 2]);
%!test
%! Y = [100; 200; 300; 400; 500; 600];
%! Group = categorical ({'Group1'; 'Group1'; 'Group2'; 'Group2'; ...
%!                       'Group3'; 'Group3'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'Group1'; 'Group2'; 'Group3'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Y, [150; 350; 550]);
%!test
%! Var1 = [1; 2; 3; 4];
%! Var2 = [10; 20; 30; 40];
%! Var3 = [100; 200; 300; 400];
%! Var4 = [1000; 2000; 3000; 4000];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'});
%! tbl = table (Var1, Var2, Var3, Var4, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!         'mean_Var1', 'mean_Var2', 'mean_Var3', 'mean_Var4'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert (stats_tbl.GroupCount, [2; 2]);
%! assert (stats_tbl.mean_Var1, [1.5; 3.5]);
%! assert (stats_tbl.mean_Var2, [15; 35]);
%! assert (stats_tbl.mean_Var3, [150; 350]);
%! assert (stats_tbl.mean_Var4, [1500; 3500]);
%!test
%! Y = [1.5; 2.5; 3.5; 4.5; 5.5; 6.5];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Y, [2; 4; 6]);
%!test
%! Y = [-10; -20; 30; 40; 50; 60];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Y, [-15; 35; 55]);
%!test
%! Y = [0; 0; 0; 0; 0; 0];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Y, [0; 0; 0]);
%!test
%! Y = [1e6; 2e6; 3e6; 4e6; 5e6; 6e6];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Y, [1.5e6; 3.5e6; 5.5e6]);
%!test
%! Y = (1:10)';
%! Group = categorical (repmat ({'A'; 'B'}, 5, 1));
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, ...
%!   {'Group', 'GroupCount', 'mean_Y', 'numel_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert (stats_tbl.GroupCount, [5; 5]);
%! assert (stats_tbl.mean_Y, [5; 6]);
%! assert (stats_tbl.numel_Y, [5; 5]);
%!test
%! Y = (1:20)';
%! Group = categorical (repmat ({'A'; 'B'; 'C'; 'D'}, 5, 1));
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'; 'D'});
%! assert (stats_tbl.GroupCount, [5; 5; 5; 5]);
%! assert (stats_tbl.mean_Y, [9; 10; 11; 12]);
%!test
%! Y = [1; 2; 3; 4; 5];
%! Group = categorical ({'A'; 'B'; 'C'; 'D'; 'E'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'; 'D'; 'E'});
%! assert (stats_tbl.GroupCount, [1; 1; 1; 1; 1]);
%! assert (stats_tbl.mean_Y, [1; 2; 3; 4; 5]);
%!test
%! Score1 = [85; 90; 78; 92; 88; 76];
%! Score2 = [82; 88; 75; 90; 85; 73];
%! Group = categorical ({'High'; 'High'; 'Med'; 'Med'; 'Low'; 'Low'});
%! tbl = table (Score1, Score2, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Score1', 'mean_Score2'});
%! assert (stats_tbl.Properties.RowNames, {'High'; 'Low'; 'Med'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Score1, [87.5; 82; 85]);
%! assert (stats_tbl.mean_Score2, [85; 79; 82.5]);
%!test
%! Height = [170; 175; 165; 180; 160; 185];
%! Weight = [70; 75; 65; 80; 60; 85];
%! Category = categorical ({'M'; 'M'; 'F'; 'F'; 'M'; 'M'});
%! tbl = table (Height, Weight, Category);
%! stats_tbl = grpstats (tbl, 'Category', {'mean', 'numel'});
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Category', 'GroupCount', ...
%!         'mean_Height', 'numel_Height', 'mean_Weight', 'numel_Weight'});
%! assert (stats_tbl.Properties.RowNames, {'F'; 'M'});
%! assert (stats_tbl.GroupCount, [2; 4]);
%! assert (stats_tbl.mean_Height, [172.5; 172.5]);
%! assert (stats_tbl.numel_Height, [2; 4]);
%! assert (stats_tbl.mean_Weight, [72.5; 72.5]);
%! assert (stats_tbl.numel_Weight, [2; 4]);
%!test
%! Value = [10.5; 11.2; 9.8; 10.1; 11.5; 10.8];
%! Type = categorical ({'A'; 'A'; 'A'; 'B'; 'B'; 'B'});
%! tbl = table (Value, Type);
%! stats_tbl = grpstats (tbl, 'Type', 'numel');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Type', 'GroupCount', ...
%!                                              'numel_Value'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert (stats_tbl.GroupCount, [3; 3]);
%! assert (stats_tbl.numel_Value, [3; 3]);
%!test
%! Data = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10];
%! Label = categorical ({'A'; 'A'; 'A'; 'A'; 'A'; 'B'; 'B'; 'B'; 'B'; 'B'});
%! tbl = table (Data, Label);
%! stats_tbl = grpstats (tbl, 'Label', {'mean', 'numel'});
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Label', 'GroupCount', ...
%!                                              'mean_Data', 'numel_Data'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert (stats_tbl.GroupCount, [5; 5]);
%! assert (stats_tbl.mean_Data, [3; 8]);
%! assert (stats_tbl.numel_Data, [5; 5]);
%!test
%! X1 = [1; 2; 3; 4];
%! X2 = [5; 6; 7; 8];
%! G = categorical ({'A'; 'A'; 'B'; 'B'});
%! tbl = table (X1, X2, G);
%! stats_tbl = grpstats (tbl, 'G', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'G', 'GroupCount', ...
%!                                              'mean_X1', 'mean_X2'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'});
%! assert (stats_tbl.GroupCount, [2; 2]);
%! assert (stats_tbl.mean_X1, [1.5; 3.5]);
%! assert (stats_tbl.mean_X2, [5.5; 7.5]);
%!test
%! Measurement = [100; 150; 200; 250; 300; 350];
%! GroupVar = categorical ({'Control'; 'Control'; 'Treatment'; 'Treatment'; ...
%!                          'Placebo'; 'Placebo'});
%! tbl = table (Measurement, GroupVar);
%! stats_tbl = grpstats (tbl, 'GroupVar', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'GroupVar', 'GroupCount', ...
%!                                              'mean_Measurement'});
%! assert (stats_tbl.Properties.RowNames, {'Control'; 'Placebo'; 'Treatment'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_Measurement, [125; 325; 225]);
%!test
%! Y = [NaN; NaN; 3; 4; 5; 6];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', 'mean');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (isequaln (stats_tbl.mean_Y, [NaN; 3.5; 5.5]));
%!test
%! Y = [1; 2; NaN; NaN; NaN; NaN];
%! Group = categorical ({'A'; 'A'; 'B'; 'B'; 'C'; 'C'});
%! tbl = table (Y, Group);
%! stats_tbl = grpstats (tbl, 'Group', {'mean', 'numel'});
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Group', 'GroupCount', ...
%!                                              'mean_Y', 'numel_Y'});
%! assert (stats_tbl.Properties.RowNames, {'A'; 'B'; 'C'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (isequaln (stats_tbl.mean_Y, [1.5; NaN; NaN]));
%! assert (stats_tbl.numel_Y, [2; 0; 0]);
%!test
%! Val = [5.5; 6.5; 7.5; 8.5];
%! Cat = categorical ({'X'; 'X'; 'Y'; 'Y'});
%! tbl = table (Val, Cat);
%! stats_tbl = grpstats (tbl, 'Cat', 'numel');
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Cat', 'GroupCount', ...
%!                                              'numel_Val'});
%! assert (stats_tbl.Properties.RowNames, {'X'; 'Y'});
%! assert (stats_tbl.GroupCount, [2; 2]);
%! assert (stats_tbl.numel_Val, [2; 2]);
%!test
%! A = [1; 2; 3; 4; 5; 6];
%! B = [10; 20; 30; 40; 50; 60];
%! C = [100; 200; 300; 400; 500; 600];
%! Grp = categorical ({'G1'; 'G1'; 'G2'; 'G2'; 'G3'; 'G3'});
%! tbl = table (A, B, C, Grp);
%! stats_tbl = grpstats (tbl, 'Grp', {'mean', 'numel'});
%! assert (istable (stats_tbl));
%! assert (stats_tbl.Properties.VariableNames, {'Grp', 'GroupCount', ...
%!         'mean_A', 'numel_A', 'mean_B', 'numel_B', 'mean_C', 'numel_C'});
%! assert (stats_tbl.Properties.RowNames, {'G1'; 'G2'; 'G3'});
%! assert (stats_tbl.GroupCount, [2; 2; 2]);
%! assert (stats_tbl.mean_A, [1.5; 3.5; 5.5]);
%! assert (stats_tbl.numel_A, [2; 2; 2]);
%! assert (stats_tbl.mean_B, [15; 35; 55]);
%! assert (stats_tbl.numel_B, [2; 2; 2]);
%! assert (stats_tbl.mean_C, [150; 350; 550]);
%! assert (stats_tbl.numel_C, [2; 2; 2]);

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
%!       grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, "predci", "somename", -0.1);
%!error <grpstats: invalid data types for 'VarNames'.> ...
%!       grpstats (ones (6, 2), [1; 1; 1; 2; 2; 2], 'mean', 'VarNames', 3)
%!error <grpstats: X must be numeric to plot mean and CI for each group.> ...
%!       grpstats ({ones(6, 2)}, [], 0.05)
%!error <grpstats: 'alpha' must be a real scalar in the range \(0,1\).> ...
%!       grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, "predci", "alpha", -0.1);
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
