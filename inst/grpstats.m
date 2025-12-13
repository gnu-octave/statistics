## Copyright (C) 2022-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 jayantchauhan <0001jayant@gmail.com>
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
## Summary statistics by group.  @code{grpstats} computes groupwise summary
## statistics, for data in a matrix @var{x}.  @code{grpstats} treats NaNs as
## missing values, and removes them.
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

function [varargout] = grpstats (x, group, whichstats, varargin)
  ## Check input arguments
  narginchk (1, 5)
  ## Table input (datatypes integration)
  if (istable (x))
    if (nargout > 1)
      error ("grpstats: table input currently supports a single output argument.");
    endif
    varargout{1} = __grpstats_table__ (x, group, whichstats, varargin{:});
    return;
  endif

  ## Numeric matrix input (existing behaviour)
  ## Check X being a vector or 2d matrix of real values
  if (ndims (x) > 2 || ! isnumeric (x) || islogical (x))
    error ("grpstats: X must be a vector or 2d matrix.");
  endif
  ## If X is a vector, make it a column vector
  if (isvector (x))
    x = x(:);
  endif
  ## Check groups and if empty make a single group for all X
  [r, c] = size (x);
  if (nargin < 2 || isempty (group))
    if (isempty (x))
      [varargout] = repmat ({[]}, nargout, 1);
      return
    endif
    group = ones (r, 1);
  endif
  ## Get group names and indices
  [group_idx, group_names] = grp2idx (group);
  ngroups = length (group_names);
  if (length (group_idx) != r)
    error ("grpstats: samples in X and GROUPS mismatch.");
  endif
  ## Add default for whichstats and check for 3rd input argument
  func_names = {};
  if (nargin > 2 && ischar (whichstats))
    func_names = {whichstats};
  elseif (nargin > 2 && iscell (whichstats))
    func_names = whichstats;
  endif
  ## Add default for alpha and check for 4th input argument
  if (nargin > 3)
    if (ischar (varargin{1}))
      if (strcmpi (varargin{1}, "alpha"))
        if (nargin > 4)
          alpha = varargin{2};
        else
          error ("grpstats: missing VALUE for optional 'alpha' parameter.");
        endif
      else
        error ("grpstats: invalid NAME for optional 'alpha' parameter.");
      endif
    elseif (isnumeric (varargin{1}))
      alpha = varargin{1};
    else
      error ("grpstats: invalid fourth input argument.");
    endif
    if (! (isnumeric (alpha) && isscalar (alpha) && alpha > 0 && alpha < 1))
      error ("grpstats: 'alpha' must be a real scalar in the range (0,1).");
    endif
  else
    alpha = 0.05;
  endif

  ## Calculate functions
  if (isempty (func_names))
    ## Check consistent number of output arguments
    if (nargout == 1 || nargout == 0)
      for j = 1:ngroups
        group_x = x(find (group_idx == j), :);
        group_mean(j,:) = mean (group_x, 1, "omitnan") ;
      endfor
      varargout{1} = group_mean;
    else
      error ("grpstats: inconsistent number of output arguments.");
    endif
  else
    func_num = length (func_names);
    ## Check consistent number of output arguments
    if (! (nargout == 0 && func_num == 1) && nargout != func_num)
      error ("grpstats: inconsistent number of output arguments.");
    endif
    for l = 1:func_num
      switch (func_names{l})
        case "mean"
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            group_mean(j,:) = mean (group_x, 1, "omitnan");
          endfor
          varargout{l} = group_mean;
        case "median"
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            group_mean(j,:) = median (group_x, 1, "omitnan");
          endfor
          varargout{l} = group_mean;
        case "sem"
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            group_sem(j,:) = std (group_x, 0, 1, "omitnan") / ...
                            sqrt (size (group_x, 1) - sum (isnan (group_x), 1));
          endfor
          varargout{l} = group_sem;
        case "std"
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            group_std(j,:) = std (group_x, 0, 1, "omitnan");
          endfor
          varargout{l} = group_std;
        case "var"
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            group_var(j,:) = var (group_x, 0, 1, "omitnan");
          endfor
          varargout{l} = group_var;
        case "min"
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            group_min(j,:) = nanmin (group_x);
          endfor
          varargout{l} = group_min;
        case "max"
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            group_max(j,:) = nanmax (group_x);
          endfor
          varargout{l} = group_max;
        case "range"
          func_handle = @(x) range (x, 1);
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            group_range(j,:) = range (group_x, 1);
          endfor
          varargout{l} = group_range;
        case "numel"
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            group_numel(j,:) = size (group_x, 1) - sum (isnan (group_x), 1);
          endfor
          varargout{l} = group_numel;
        case "meanci"
          ## allocate as 3-D: [ngroups x c x 2] (lower, upper)
          group_meanci = NaN (ngroups, c, 2);
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            m = mean (group_x, 1, "omitnan");
            n = size (group_x, 1) - sum (isnan (group_x), 1);
            s = std (group_x, 0, 1, "omitnan") ./ sqrt (max (n,1));
            ## avoid invalid tinv calls for degenerate df
            df = max (n - 1, 0);
            tval = zeros (1, size (group_x, 2));
            pos = (df > 0);
            if (any (pos))
              tval(pos) = - tinv (alpha / 2, df(pos));
            endif
            d = s .* tval;
            group_meanci(j, :, 1) = m - d;
            group_meanci(j, :, 2) = m + d;
          endfor
          ## MATLAB returns [ngroups x 2] when nvars == 1; Octave used canonical 3-D.
          if (c == 1)
            ## reshape to [ngroups x 2]
            varargout{l} = reshape (group_meanci, ngroups, 2);
          else
            varargout{l} = group_meanci;
          endif

        case "predci"
          ## allocate as 3-D: [ngroups x c x 2] (lower, upper)
          group_predci = NaN (ngroups, c, 2);
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            m = mean (group_x, 1, "omitnan");
            n = size (group_x, 1) - sum (isnan (group_x), 1);
            s = std (group_x, 0, 1, "omitnan") ./ sqrt (1 + (1 ./ max (n,1)));
            df = max (n - 1, 0);
            tval = zeros (1, size (group_x, 2));
            pos = (df > 0);
            if (any (pos))
              tval(pos) = - tinv (alpha / 2, df(pos));
            endif
            d = s .* tval;
            group_predci(j, :, 1) = m - d;
            group_predci(j, :, 2) = m + d;
          endfor
          if (c == 1)
            varargout{l} = reshape (group_predci, ngroups, 2);
          else
            varargout{l} = group_predci;
          endif
        case "gname"
          varargout{l} = group_names;
        otherwise
          error ("grpstats: wrong whichstats option.");
      endswitch
    endfor
  endif

endfunction

function stats_tbl = __grpstats_table__ (tbl, group, whichstats, varargin)

  if (! istable (tbl))
    error ("grpstats: internal error, expected table input.");
  endif

  ## GROUP must be a variable name (char or string scalar)
  if (ischar (group))
    group_name = group;
  elseif (isstring (group) && isscalar (group))
    group_name = char (group);
  else
    error ("grpstats: for table input, GROUP must be a variable name.");
  endif

  ## Get grouping column
  try
    gcol = tbl.(group_name);
  catch
    error ("grpstats: grouping variable '%s' not found in table.", group_name);
  end_try_catch

  ## Currently only support categorical grouping variable
  if (! iscategorical (gcol))
    error ("grpstats: for table input, grouping variable must be categorical.");
  endif

  ## Normalise whichstats for table input
  if (nargin < 3 || isempty (whichstats))
    func_names = {"mean"};
  elseif (ischar (whichstats))
    func_names = {whichstats};
  elseif (isstring (whichstats) && isscalar (whichstats))
    func_names = {char (whichstats)};
  elseif (iscell (whichstats))
    func_names = whichstats;
  else
    error ("grpstats: invalid WHICHSTATS for table input.");
  endif

  ## Only support a subset initially for table input
  n_funcs = numel (func_names);
  for k = 1:n_funcs
    fname = func_names{k};
    if (! any (strcmp (fname, {"mean", "numel"})))
      error ("grpstats: table input currently supports only 'mean' and 'numel'.");
    endif
  endfor

  ## Group indices and names from categorical column
  group_names = categories (gcol);
  group_idx   = double (gcol(:));
  ngroups     = numel (group_names);

  ## Collect numeric data columns (excluding the grouping variable)
  vnames = tbl.Properties.VariableNames;
  data_var_names = {};
  data_mat = [];

  for k = 1:numel (vnames)
    vname = vnames{k};
    if (strcmp (vname, group_name))
      continue;
    endif
    col = tbl.(vname);
    if (isnumeric (col))
      data_mat = [data_mat, col(:)];
      data_var_names{end+1} = vname;
    endif
  endfor

  if (isempty (data_mat))
    error ("grpstats: no numeric variables found in table (apart from grouping).");
  endif

  nvars = columns (data_mat);

  do_mean  = any (strcmp ("mean",  func_names));
  do_numel = any (strcmp ("numel", func_names));

  if (do_mean)
    mean_vals = NaN (ngroups, nvars);
  endif
  if (do_numel)
    group_count = accumarray (group_idx(:), 1, [ngroups, 1]);
  endif

  ## Compute statistics per group
  for g = 1:ngroups
    idx = (group_idx == g);
    group_data = data_mat(idx, :);
    if (do_mean)
      mean_vals(g,:) = mean (group_data, 1, "omitnan");
    endif
  endfor

  ## Build output table
  ## Group column as categorical using group names
  gcat = categorical (group_names);

  varnames_out = {"Group"};
  data_out = {gcat};

  if (do_numel)
    varnames_out{end+1} = "GroupCount";
    data_out{end+1} = group_count;
  endif

  if (do_mean)
    for k = 1:nvars
      newname = ["mean_" data_var_names{k}];
      varnames_out{end+1} = newname;
      data_out{end+1} = mean_vals(:, k);
    endfor
  endif

  stats_tbl = table (data_out{:}, "VariableNames", varnames_out);

endfunction

%!demo
%! load carsmall;
%! [m,p,g] = grpstats (Weight, Model_Year, {"mean", "predci", "gname"})
%! n = length(m);
%! errorbar((1:n)',m,p(:,2)-m);
%! set (gca, "xtick", 1:n, "xticklabel", g);
%! title ("95% prediction intervals for mean weight by year");

%!demo
%! load carsmall;
%! [m,p,g] = grpstats ([Acceleration,Weight/1000],Cylinders, ...
%!                     {"mean", "meanci", "gname"}, 0.05)
%! [c,r] = size (m);
%! errorbar((1:c)'.*ones(c,r),m,p(:,[(1:r)])-m);
%! set (gca, "xtick", 1:c, "xticklabel", g);
%! title ("95% prediction intervals for mean weight by year");

%!test
%! load carsmall
%! means = grpstats (Acceleration, Origin);
%! assert (means, [14.4377; 18.0500; 15.8867; 16.3778; 16.6000; 15.5000], 0.001);
%!test
%! load carsmall
%! [grpMin,grpMax,grp] = grpstats (Acceleration, Origin, {"min","max","gname"});
%! assert (grpMin, [8.0; 15.3; 13.9; 12.2; 15.7; 15.5]);
%! assert (grpMax, [22.2; 21.9; 18.2; 24.6; 17.5; 15.5]);
%!test
%! load carsmall
%! [grpMin,grpMax,grp] = grpstats (Acceleration, Origin, {"min","max","gname"});
%! assert (grp', {"USA", "France", "Japan", "Germany", "Sweden", "Italy"});
%!test
%! load carsmall
%! [m,p,g] = grpstats ([Acceleration,Weight/1000], Cylinders, ...
%!                     {"mean", "meanci", "gname"}, 0.05);
%! ## check meanci lower bounds (first slice) with tolerance
%! expected_lower = [15.9163; 15.6622; 10.7968]; 
%! expected_upper = [17.4249; 17.2907; 12.4845]; 
%! assert (abs (p(:,1,1) - expected_lower) < 1e-3);
%! assert (abs (p(:,1,2) - expected_upper) < 1e-3);
%!test
%! [mC, g] = grpstats ([], []);
%! assert (isempty (mC), true);
%! assert (isempty (g), true);

## Table input tests (datatypes integration)
%!test
%! Y = [5; 6; 7; 4; 9; 8];
%! X = [1; 2; 3; 4; 5; 6];
%! Group = categorical ({"A"; "A"; "B"; "B"; "C"; "C"});
%! tbl = table (Y, X, Group, "VariableNames", {"Y","X","Group"});
%! stats_tbl = grpstats (tbl, "Group", {"mean","numel"});
%! assert (istable (stats_tbl));
%! assert (isequal (stats_tbl.Properties.VariableNames, ...
%!                  {"Group", "GroupCount", "mean_Y", "mean_X"}));
%! assert (isequal (stats_tbl.GroupCount, [2; 2; 2]));
%!test
%! Y     = [5; 6; 7; 4; 9; 8];
%! Group = categorical ({"A"; "A"; "B"; "B"; "C"; "C"});
%! tbl = table (Y, Group, "VariableNames", {"Y","Group"});
%! stats_tbl = grpstats (tbl, "Group", "mean");
%! assert (istable (stats_tbl));
%! assert (isequal (stats_tbl.Properties.VariableNames, ...
%!                  {"Group", "mean_Y"}));

%!error<grpstats: for table input, grouping variable must be categorical.> ...
%! Y = [1; 2; 3];
%! G = [1; 1; 2];
%! tbl = table (Y, G, "VariableNames", {"Y","G"});
%! grpstats (tbl, "G", {"mean","numel"});
%!error<grpstats: X must be a vector or 2d matrix.> ...
%! grpstats (ones (3, 3, 3));
%!error<grpstats: samples in X and GROUPS mismatch.> ...
%! grpstats ([], {'A'; 'B'; 'A'; 'B'})
%!error<grpstats: missing VALUE for optional 'alpha' parameter.> ...
%! grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, "predci", "alpha");
%!error<grpstats: invalid NAME for optional 'alpha' parameter.> ...
%! grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, "predci", "somename", -0.1);
%!error<grpstats: invalid fourth input argument.> ...
%! grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, "predci", {2, 3}, -0.1);
%!error<grpstats: 'alpha' must be a real scalar in the range \(0,1\).> ...
%! grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, "predci", "alpha", -0.1);
%!error<grpstats: inconsistent number of output arguments.> ...
%! grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, {'mean', 'sum'});
%!error<grpstats: inconsistent number of output arguments.> ...
%! [q, w] = grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'});
%!error<grpstats: wrong whichstats option.> ...
%! grpstats ([1:4]', {'A'; 'B'; 'A'; 'B'}, "whatever");
