## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## and "predci" options.  Default value for @var{alpha} is 0.05.
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

function [varargout] = grpstats (x ,group, whichstats, alpha)
  ## Check input arguments
  narginchk (1, 4)
  ## Check X being a vector or 2d matrix of real values
  if (ndims (x) > 2 || ! isnumeric (x) || islogical (x))
    error ("grpstats: X must be a vector or 2d matrix of real values.");
  endif
  ## If X is a vector, make it a column vector
  if (isvector (x))
    x = x(:);
  endif
  ## Check groups and if empty make a single group for all X
  [r, c] = size (x);
  if (nargin < 2 || isempty (group))
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
    if (! (isnumeric (alpha) && isscalar (alpha) ...
          && alpha > 0 && alpha < 1))
      error ("grpstats: ALPHA must be a real scalar in the range (0,1).");
    endif
  else
    alpha = 0.05;
  endif

  ## Calculate functions
  if (isempty (func_names))
    ## Check consistent number of output arguments
    if (nargout == 1)
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
    if (nargout != func_num)
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
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            m = mean (group_x, 1, "omitnan") ;
            n = size (x, 1) - sum (isnan (group_x), 1);
            s = std (group_x, 0, 1, "omitnan") ./ sqrt (n);
            d = s .* - tinv (alpha / 2, max (n - 1, [], 1));
            group_meanci(j,:) = [m-d, m+d];
          endfor
          varargout{l} = group_meanci;
        case "predci"
          for j = 1:ngroups
            group_x = x(find (group_idx == j), :);
            m = mean (group_x, 1, "omitnan") ;
            n = size (x, 1) - sum (isnan (group_x), 1);
            s = std (group_x, 0, 1, "omitnan") ./ sqrt (1 + (1 ./ n));
            d = s .* - tinv (alpha / 2, max (n - 1, [], 1));
            group_predci(j,:) = [m-d, m+d];
          endfor
          varargout{l} = group_predci;
        case "gname"
          varargout{l} = group_names;
        otherwise
          error ("grpstats: wrong whichstats option.");
      endswitch
    endfor
  endif
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
%! assert (p(:,1), [11.17621760075134, 16.13845847655224, 16.16222663683362]', ...
%!                 [1e-14, 2e-14, 1e-14]');
