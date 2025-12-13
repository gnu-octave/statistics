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
      error ("grpstats: table input supports a single output argument.");
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
            s = std (group_x, 0, 1, "omitnan") .* sqrt (1 + (1 ./ max (n,1)));
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
%! m = grpstats (x, g, "mean");
%! expected = [15; 35; 55];
%! assert (m, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! m = grpstats (x, g, "median");
%! expected = [15; 35; 55];
%! assert (m, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! s = grpstats (x, g, "std");
%! expected = [7.07106781186548; 7.07106781186548; 7.07106781186548];
%! assert (s, expected, 1e-14);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! v = grpstats (x, g, "var");
%! expected = [50; 50; 50];
%! assert (v, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! s = grpstats (x, g, "sem");
%! expected = [5; 5; 5];
%! assert (s, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! mn = grpstats (x, g, "min");
%! expected = [10; 30; 50];
%! assert (mn, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! mx = grpstats (x, g, "max");
%! expected = [20; 40; 60];
%! assert (mx, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! r = grpstats (x, g, "range");
%! expected = [10; 10; 10];
%! assert (r, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! n = grpstats (x, g, "numel");
%! expected = [2; 2; 2];
%! assert (n, expected);
%!test
%! ## single statistic
%! x = [10; 20; 30; 40; 50; 60];
%! g = {"A"; "A"; "B"; "B"; "C"; "C"};
%! names = grpstats (x, g, "gname");
%! expected = {"A"; "B"; "C"};
%! assert (names, expected);
%!test
%! ## single statistic (default alpha)
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, "meanci");
%! expected = [-48.5310236808735 78.5310236808735; -28.5310236808735 ...
%!             98.5310236808735; -8.53102368087348 118.531023680873];
%! assert (ci, expected, 1e-12);
%!test
%! ## single statistic (default alpha)
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, "predci");
%! expected = [-95.0389608721344 125.038960872134; -75.0389608721344 ...
%!             145.038960872134; -55.0389608721344 165.038960872134];
%! assert (ci, expected, 1e-12);
%!test
%! ## mean + std
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [m, s] = grpstats (x, g, {"mean", "std"});
%! expected_m = [15; 35; 55];
%! expected_s = [7.07106781186548; 7.07106781186548; 7.07106781186548];
%! assert (m, expected_m, 1e-14);
%! assert (s, expected_s, 1e-14);
%!test
%! ## min + max + range
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [mn, mx, r] = grpstats (x, g, {"min", "max", "range"});
%! expected_mn = [10; 30; 50];
%! expected_mx = [20; 40; 60];
%! expected_r  = [10; 10; 10];
%! assert (mn, expected_mn);
%! assert (mx, expected_mx);
%! assert (r, expected_r);
%!test
%! ## mean + median + numel + gname
%! x = [10; 20; 30; 40; 50; 60];
%! g = {"A"; "A"; "B"; "B"; "C"; "C"};
%! [m, med, n, names] = grpstats (x, g, {"mean", "median", "numel", "gname"});
%! expected_m   = [15; 35; 55];
%! expected_med = [15; 35; 55];
%! expected_n   = [2; 2; 2];
%! expected_names = {"A"; "B"; "C"};
%! assert (m, expected_m);
%! assert (med, expected_med);
%! assert (n, expected_n);
%! assert (names, expected_names);
%!test
%! ## all basic statistics
%! x = [10; 20; 30; 40; 50; 60; 70; 80];
%! g = [1; 1; 2; 2; 2; 2; 3; 3];
%! [m, med, s, v, se, mn, mx, r, n] = grpstats (x, g, {"mean", "median", ...
%!                                                     "std", "var", "sem", ...
%!                                                     "min", "max", "range", ...
%!                                                     "numel"});
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
%! ci = grpstats (x, g, "meanci", 0.1);
%! expected = [-16.5687575733752 46.5687575733752; 3.4312424266248 ...
%!             66.5687575733752; 23.4312424266248 86.5687575733752];
%! assert (ci, expected, 1e-13);
%!test
%! ## predci-alpha-0.1
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, "predci", 0.1);
%! expected = [-39.6786920489106 69.6786920489106; -19.6786920489106 ...
%!             89.6786920489106; 0.321307951089366 109.678692048911];
%! assert (ci, expected, 1e-12);
%!test
%! ## meanci-alpha-0.01
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, "meanci", 0.01);
%! expected = [-303.283705814358 333.283705814358; -283.283705814358 ...
%!             353.283705814358; -263.283705814358 373.283705814358];
%! assert (ci, expected, 3e-8);
%!test
%! ## predci-alpha-0.01
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, "predci", 0.01);
%! expected = [-536.283549691775 566.283549691775; -516.283549691775 ...
%!             586.283549691775; -496.283549691775 606.283549691775];
%! assert (ci, expected, 3e-8);
%!test
%! ## meanci-alpha-0.2
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, "meanci", 0.2);
%! expected = [-0.388417685876263 30.3884176858763; 19.6115823141237 ...
%!             50.3884176858763; 39.6115823141237 70.3884176858763];
%! assert (ci, expected, 1e-13);
%!test
%! ## predci-alpha-0.2
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, "predci", 0.2);
%! expected = [-11.6535212800292 41.6535212800292; 8.34647871997083 ...
%!             61.6535212800292; 28.3464787199708 81.6535212800292];
%! assert (ci, expected, 1e-13);
%!test
%! ## meanci, name-value alpha=0.2
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! ci = grpstats (x, g, "meanci", "alpha", 0.2);
%! expected = [-0.388417685876263 30.3884176858763; 19.6115823141237 ...
%!             50.3884176858763; 39.6115823141237 70.3884176858763];
%! assert (ci, expected, 1e-13);
%!test
%! ## meanci + predci, alpha=0.01
%! x = [10; 20; 30; 40; 50; 60];
%! g = [1; 1; 2; 2; 3; 3];
%! [ci_m, ci_p] = grpstats (x, g, {"meanci", "predci"}, 0.01);
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
%! [m, s, n] = grpstats (x, g, {"mean", "std", "numel"});
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
%! [m, n] = grpstats (x, g, {"mean", "numel"});
%! expected_m = [1   15; 3.5 40; 5.5 55];
%! expected_n = [1 2; 2 1; 2 2];
%! assert (m, expected_m);
%! assert (n, expected_n);
%!test
%! ## 3-column matrix, mean+min+max
%! x = [1 100 1000; 2 200 2000; 3 300 3000; 4 400 4000];
%! g = [1; 1; 2; 2];
%! [m, mn, mx] = grpstats (x, g, {"mean", "min", "max"});
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
%! [m, s, n] = grpstats (x, g, {"mean", "std", "numel"});
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
%! [m, s, n] = grpstats (x, g, {"mean", "std", "numel"});
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
%! [m, v, n] = grpstats (x, g, {"mean", "var", "numel"});
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
%! [m, names] = grpstats (x, g, {"mean", "gname"});
%! expected_m = [15; 35; 55];
%! expected_names = {"1"; "5"; "10"};
%! assert (m, expected_m);
%! assert (names, expected_names);
%!test
%! ## unsorted string groups
%! x = [30; 10; 40; 20; 60; 50];
%! g = {"C"; "A"; "C"; "A"; "B"; "B"};
%! [m, names] = grpstats (x, g, {"mean", "gname"});
%! expected_m = [35; 15; 55];
%! expected_names = {"C"; "A"; "B"};
%! assert (m, expected_m);
%! assert (names, expected_names);
%!test
%! ## 20 groups, one element each
%! x = (1:20)';
%! g = (1:20)';
%! [m, n] = grpstats (x, g, {"mean", "numel"});
%! expected_m = (1:20)';
%! expected_n = ones (20, 1);
%! assert (m, expected_m);
%! assert (n, expected_n);
%!test
%! ## large sample meanci
%! x = (1:50)';
%! g = [ones(25, 1); 2 * ones(25, 1)];
%! ci = grpstats (x, g, "meanci");
%! expected = [9.96202357522388 16.0379764247761; 34.9620235752239 ...
%!             41.0379764247761];
%! assert (ci, expected, 1e-13);
%!test
%! ## large sample predci
%! x = (1:50)';
%! g = [ones(25, 1); 2 * ones(25, 1)];
%! ci = grpstats (x, g, "predci");
%! expected = [-2.49070107176829 28.4907010717683; 22.5092989282317 ...
%!             53.4907010717683];
%! assert (ci, expected, 1e-14);

## Test input validation
%!error<grpstats: table input supports a single output argument.> ...
%! [a, b] = grpstats (table (1));
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
