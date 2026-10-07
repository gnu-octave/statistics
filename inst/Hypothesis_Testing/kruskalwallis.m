## Copyright (C) 2021 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{p} =} kruskalwallis (@var{x})
## @deftypefnx {statistics} {@var{p} =} kruskalwallis (@var{x}, @var{group})
## @deftypefnx {statistics} {@var{p} =} kruskalwallis (@var{x}, @var{group}, @var{displayopt})
## @deftypefnx {statistics} {@var{p} =} kruskalwallis (@var{x}, @var{group}, @var{displayopt}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{p}, @var{tbl}] =} kruskalwallis (@var{x}, @dots{})
## @deftypefnx {statistics} {[@var{p}, @var{tbl}, @var{stats}] =} kruskalwallis (@var{x}, @dots{})
##
## Perform a Kruskal-Wallis test, the non-parametric alternative of a one-way
## analysis of variance (ANOVA), for comparing the means of two or more groups
## of data under the null hypothesis that the groups are drawn from the same
## population, i.e. the group means are equal.
##
## kruskalwallis can take up to three input arguments:
##
## @itemize
## @item
## @var{x} contains the data and it can either be a vector or matrix.
## If @var{x} is a matrix, then each column is treated as a separate group.
## If @var{x} is a vector, it is a single group whichever way it lies, unless
## @var{group} is given; with one group there is nothing to compare, so
## @var{p} is @qcode{NaN}.
## @item
## @var{group} contains the names for each group.  If @var{x} is a matrix, then
## @var{group} can either be a cell array of strings of a character array, with
## one row per column of @var{x}.  If you want to omit this argument, enter an
## empty array ([]).  If @var{x} is a vector, then @var{group} must be a vector
## of the same length, or a string array or cell array of strings with one row
## for each element of @var{x}.  @var{x} values corresponding to the same value
## of @var{group} are placed in the same group.
## @item
## @var{displayopt} is an optional parameter for displaying the groups contained
## in the data in a boxplot.  If omitted, it is 'on' by default.  If group names
## are defined in @var{group}, these are used to identify the groups in the
## boxplot. Use 'off' to omit displaying this figure.
## @end itemize
##
## kruskalwallis can return up to three output arguments:
##
## @itemize
## @item
## @var{p} is the p-value of the null hypothesis that all group means are equal.
## @item
## @var{tbl} is a cell array containing the results in a standard ANOVA table.
## @item
## @var{stats} is a structure containing statistics useful for performing
## a multiple comparison of means with the MULTCOMPARE function.
## @end itemize
##
## @var{stats} also holds the effect size of the test, an Octave extension:
## @qcode{EtaSquared}, the rank eta squared
## @math{max (0, (H - k + 1) / (n - k))} for the statistic @math{H} over
## @math{k} groups and @math{n} observations, which is the share of the
## variance of the ranks the groups account for, corrected for bias, and
## @qcode{EtaSquaredCI}, its two-sided confidence interval at the level
## @math{100 (1 - alpha)} percent.  By default the interval inverts the
## noncentral chi-square distribution of @math{H} on @math{k - 1} degrees of
## freedom, taken as @math{H}'s distribution under the alternative, which it
## is for large samples.  Three options after @var{displayopt}, also an
## Octave extension, set it:
##
## @multitable @columnfractions 0.32 0.66
## @headitem @var{Name} @tab @var{Value}
## @item @qcode{'Alpha'} @tab The significance level of the interval, a
## scalar between 0 and 1, 0.05 by default.
## @item @qcode{'ConfidenceIntervalType'} @tab @qcode{'exact'}, the default,
## or @qcode{'bootstrap'} for the bias-corrected and accelerated bootstrap
## interval, resampling within each group.
## @item @qcode{'NumBootstraps'} @tab The number of bootstrap replicates, a
## positive integer, 1000 by default.
## @end multitable
##
## If kruskalwallis is called without any output arguments, then it prints the
## results in a one-way ANOVA table to the standard output.  It is also printed
## when @var{displayopt} is 'on'.
##
## Examples:
##
## @example
## x = meshgrid (1:6);
## x = x + normrnd (0, 1, 6, 6);
## [p, atab] = kruskalwallis(x);
## @end example
##
##
## @example
## x = ones (50, 4) .* [-2, 0, 1, 5];
## x = x + normrnd (0, 2, 50, 4);
## group = @{"A", "B", "C", "D"@};
## kruskalwallis (x, group);
## @end example
##
## @seealso{anova1, friedman, ranksum, multcompare}
## @end deftypefn

function [p, tbl, stats] = kruskalwallis (x, group, displayopt, varargin)

  ## check for valid number of input arguments
  if (nargin < 1)
    print_usage ();
  endif
  ## add defaults
  if (nargin < 2)
    group = [];
  endif
  if (nargin < 3)
    displayopt = 'on';
  endif
  [alpha, citype, nboot, args] = parsePairedArguments ...
               ({'Alpha', 'ConfidenceIntervalType', 'NumBootstraps'}, ...
                {0.05, 'exact', 1000}, varargin(:));
  if (! isempty (args))
    error ("kruskalwallis: invalid optional paired argument.");
  endif
  if (! (isnumeric (alpha) && isreal (alpha) && isscalar (alpha)
         && alpha > 0 && alpha < 1))
    error ("kruskalwallis: 'Alpha' must be a scalar between 0 and 1.");
  endif
  if (! (ischar (citype) && any (strcmpi (citype, {'exact', 'bootstrap'}))))
    error (strcat ("kruskalwallis: 'ConfidenceIntervalType' must be", ...
                   " 'exact' or 'bootstrap'."));
  endif
  if (! (isnumeric (nboot) && isreal (nboot) && isscalar (nboot)
         && isfinite (nboot) && nboot >= 1 && nboot == fix (nboot)))
    error ("kruskalwallis: 'NumBootstraps' must be a positive integer.");
  endif
  plotdata = ! (strcmp (displayopt, 'off'));

  ## Convert group to cell array from character array, make it a column
  if (! isempty (group) && ischar (group))
    group = cellstr (group);
  endif
  if (size (group, 1) == 1)
    group = group';
  endif

  ## A vector is one group whichever way it lies, and each column of a matrix
  ## is a group of its own.  With GROUP given, a vector is grouped by it
  ## instead, one element per observation.
  if (isvector (x))
    x = x(:);
  endif
  if (isempty (group) || ! isvector (x))
    [n, m] = size (x);
    x = x(:);
    gi = reshape (repmat ((1:m), n, 1), n*m, 1);
    if (length (group) == 0)          ## no group names are provided
      group = gi;
    elseif (size (group, 1) == m)     ## group names exist and match columns
      group = group(gi,:);
    else
      error ("X columns and GROUP length do not match.");
    endif
  endif

  ## Identify NaN values (if any) and remove them from X along with
  ## their corresponding values from group vector
  nonan = ! isnan (x);
  x = x(nonan);
  group = group(nonan, :);

  ## Convert group to indices and separate names
  [group_id, group_names] = grp2idx (group);
  group_id = group_id(:);
  named = 1;

  ## Rank data for non-parametric analysis
  [xr, tieadj] = tieranks (x);

  ## Get group size and mean for each group
  groups = size (group_names, 1);
  xs = zeros (1, groups);
  xm = xs;
  for j = 1:groups
    group_size = find (group_id == j);
    xs(j) = length (group_size);
    xm(j) = mean (xr(group_size));
  endfor

  ## Calculate statistics
  lx = length (xr);                       ## Number of samples in groups
  gm = mean (xr);                         ## Grand mean of groups
  dfm = length (xm) - 1;                  ## degrees of freedom for model
  dfe = lx - dfm - 1;                     ## degrees of freedom for error
  ## With no observations there is nothing to partition.  Counting the groups
  ## less one, and the samples less that, puts both below zero, which has no
  ## reading; MATLAB reports the negative numbers.
  if (lx == 0)
    dfm = 0;
    dfe = 0;
  endif
  SSM = xs .* (xm - gm) * (xm - gm)';     ## Sum of Squares for Model
  SST = (xr(:) - gm)' * (xr(:) - gm);     ## Sum of Squares Total
  SSE = SST - SSM;                        ## Sum of Squares Error
  if (dfm > 0)
      MSM = SSM / dfm;                    ## Mean Square for Model
  else
      MSM = NaN;
  endif
  if (dfe > 0)
    MSE = SSE / dfe;                      ## Mean Squared Error
  else
    MSE = NaN;
  endif
  ## Calculate Chi-sq statistic
  ChiSq = (12 * SSM) / (lx * (lx + 1));
  if (tieadj > 0)
    ChiSq = ChiSq / (1 - 2 * tieadj / (lx ^ 3 - lx));
  endif
  if (dfm <= 0)
    ## Fewer than two groups leaves nothing to compare, so the statistic has
    ## no reference distribution and the test is undefined rather than
    ## significant.
    ChiSq = NaN;
    p = NaN;
  else
    p = chi2cdf (ChiSq, dfm, 'upper');
  endif

  ## Create results table (if requested)
  if (nargout > 1)
    tbl = {'Source', 'SS', 'df', 'MS', 'Chi-sq', 'Prob>Chi-sq'; ...
           'Groups', SSM, dfm, MSM, ChiSq, p; ...
           'Error', SSE, dfe, MSE, '', ''; ...
           'Total', SST, dfm + dfe, '', '', ''};
  endif
  ## Create stats structure (if requested) for MULTCOMPARE
  if (nargout > 2)
    if (length (group_names) > 0)
        stats.gnames = group_names;
    else
        stats.gnames = strjust (num2str ((1:length (xm))'), 'left');
    endif
    stats.n = xs;
    stats.source = 'kruskalwallis';
    stats.meanranks = xm;
    stats.sumt = 2 * tieadj;
    ## The rank eta squared, corrected for bias as H - df estimates the
    ## noncentrality, and its interval
    stats.EtaSquared = NaN;
    stats.EtaSquaredCI = [NaN, NaN];
    if (dfm > 0 && dfe > 0)
      stats.EtaSquared = max (0, (ChiSq - dfm) / dfe);
      if (strcmpi (citype, 'exact'))
        stats.EtaSquaredCI = min (__ncx2bounds__ (ChiSq, dfm, alpha) / dfe, 1);
      else
        S = arrayfun (@(g) x(group_id == g), 1:groups, 'UniformOutput', false);
        stats.EtaSquaredCI = __bootci__ (@etaSquared, S, alpha, nboot, ...
                                         'strata');
      endif
    endif
  endif
  ## Print results table on screen if no output argument was requested
  if (nargout == 0 || plotdata)
    printf ("              Kruskal-Wallis ANOVA Table\n");
    printf ("Source        SS      df      MS      Chi-sq  Prob>Chi-sq\n");
    printf ("---------------------------------------------------------\n");
    printf ("Columns %10.2f %5.0f %10.2f %8.2f  %11.5e\n", ...
           SSM, dfm, MSM, ChiSq, p);
    printf ("Error   %10.2f %5.0f %10.2f\n", SSE, dfe, MSE);
    printf ("Total   %10.2f %5.0f\n", SST, dfm + dfe);
  endif
  ## Plot data using BOXPLOT (unless opted out)
  if (plotdata)
    boxplot (x, group_id, 'Notch', 'on', 'Labels', group_names);
  endif
endfunction

## The rank eta squared of groups whose columns are resampled copies of each
## sample, one value per column
function e = etaSquared (varargin)
  k = numel (varargin);
  ns = cellfun (@rows, varargin);
  n = sum (ns);
  g = repelem ((1:k)', ns(:));
  X = vertcat (varargin{:});
  e = NaN (1, columns (X));
  if (n <= k)
    return;
  endif
  for j = 1:columns (X)
    [r, tieadj] = tiedrank (X(:,j));
    m = accumarray (g, r) ./ ns(:);
    H = 12 * sum (ns(:) .* (m - (n + 1) / 2) .^ 2) / (n * (n + 1));
    if (tieadj > 0)
      H /= 1 - 2 * tieadj / (n ^ 3 - n);
    endif
    e(j) = max (0, (H - k + 1) / (n - k));
  endfor
endfunction

## local function for computing tied ranks on column vectors
function [r, tieadj] = tieranks (x)
  ## Sort data
  [value, x_idx] = sort (x);
  epsx = zeros (size (x));
  epsx = epsx(x_idx);
  x_l = numel (x);
  ## Count ranks from start (min value)
  ranks = [1:x_l]';
  ## Initialize tie adjustments
  tieadj = 0;
  ## Adjust for ties.
  ties = value(1:x_l-1) + epsx(1:x_l-1) >= value(2:x_l) - epsx(2:x_l);
  t_idx = find (ties);
  t_idx(end+1) = 0;
  maxTies = numel (t_idx);
  ## Calculate tie adjustments
  tiecount = 1;
  while (tiecount < maxTies)
    tiestart = t_idx(tiecount);
    ntied = 2;
    while (t_idx(tiecount+1) == t_idx(tiecount) + 1)
      tiecount = tiecount + 1;
      ntied = ntied + 1;
    endwhile
    ## Check for tieflag
    tieadj = tieadj + ntied * (ntied - 1) * (ntied + 1) / 2;
    ## Average tied ranks
    ranks(tiestart:tiestart + ntied - 1) = ...
                  sum (ranks(tiestart:tiestart + ntied - 1)) / ntied;
    tiecount = tiecount + 1;
  endwhile
  ## Remap data to original dimensions
  r(x_idx) = ranks;
endfunction


%!demo
%! rng (42);
%! x = meshgrid (1:6);
%! x = x + normrnd (0, 1, 6, 6);
%! kruskalwallis (x, [], 'off');

%!demo
%! rng (42);
%! x = meshgrid (1:6);
%! x = x + normrnd (0, 1, 6, 6);
%! [p, atab] = kruskalwallis (x);

%!demo
%! rng (42);
%! x = ones (30, 4) .* [-2, 0, 1, 5];
%! x = x + normrnd (0, 2, 30, 4);
%! group = {'A', 'B', 'C', 'D'};
%! kruskalwallis (x, group);

## testing results against SPSS and R on the GEAR.DAT data file available from
## https://www.itl.nist.gov/div898/handbook/eda/section3/eda354.htm
%!test
%! data = [1.006, 0.996, 0.998, 1.000, 0.992, 0.993, 1.002, 0.999, 0.994, 1.000, ...
%!         0.998, 1.006, 1.000, 1.002, 0.997, 0.998, 0.996, 1.000, 1.006, 0.988, ...
%!         0.991, 0.987, 0.997, 0.999, 0.995, 0.994, 1.000, 0.999, 0.996, 0.996, ...
%!         1.005, 1.002, 0.994, 1.000, 0.995, 0.994, 0.998, 0.996, 1.002, 0.996, ...
%!         0.998, 0.998, 0.982, 0.990, 1.002, 0.984, 0.996, 0.993, 0.980, 0.996, ...
%!         1.009, 1.013, 1.009, 0.997, 0.988, 1.002, 0.995, 0.998, 0.981, 0.996, ...
%!         0.990, 1.004, 0.996, 1.001, 0.998, 1.000, 1.018, 1.010, 0.996, 1.002, ...
%!         0.998, 1.000, 1.006, 1.000, 1.002, 0.996, 0.998, 0.996, 1.002, 1.006, ...
%!         1.002, 0.998, 0.996, 0.995, 0.996, 1.004, 1.004, 0.998, 0.999, 0.991, ...
%!         0.991, 0.995, 0.984, 0.994, 0.997, 0.997, 0.991, 0.998, 1.004, 0.997];
%! group = [1:10] .* ones (10,10);
%! group = group(:);
%! [p, tbl] = kruskalwallis (data, group, 'off');
%! assert_equal (p, 0.048229, 1e-6);
%! assert_equal (tbl{2,5}, 17.03124, 1e-5);
%! assert_equal (tbl{2,3}, 9, 0);
%! assert_equal (tbl{4,2}, 82655.5, 1e-16);
%! data = reshape (data, 10, 10);
%! [p, tbl, stats] = kruskalwallis (data, [], 'off');
%! assert_equal (p, 0.048229, 1e-6);
%! assert_equal (tbl{2,5}, 17.03124, 1e-5);
%! assert_equal (tbl{2,3}, 9, 0);
%! assert_equal (tbl{4,2}, 82655.5, 1e-16);
%! means = [51.85, 60.45, 37.6, 51.1, 29.5, 54.25, 64.55, 66.7, 53.65, 35.35];
%! N = 10 * ones (1, 10);
%! assert_equal (stats.meanranks, means, 1e-6);
%! assert_equal (length (stats.gnames), 10, 0);
%! assert_equal (stats.n, N, 0);

## A vector without GROUP is one group whichever way it lies, and an empty X
## is answered rather than raised, with no degrees of freedom below zero.
%!test
%! [p, tbl] = kruskalwallis ((1:5)', [], 'off');
%! assert_equal (p, NaN);
%! assert_equal (cell2mat (tbl(2:4,2:3)), [0, 0; 10, 4; 10, 4]);

%!test
%! [p, tbl] = kruskalwallis (1:5, [], 'off');
%! assert_equal (p, NaN);
%! assert_equal (cell2mat (tbl(2:4,2:3)), [0, 0; 10, 4; 10, 4]);

%!test
%! [p, tbl] = kruskalwallis ([], [], 'off');
%! assert_equal (p, NaN);
%! assert_equal (cell2mat (tbl(2:4,2:3)), zeros (3, 2));
%!test
%! ## Below the resolution of 1 - chi2cdf, values from MATLAB R2024a
%! p = kruskalwallis ((1:200)', [ones(100, 1); 2 * ones(100, 1)], 'off');
%! assert_equal (p, 2.5239394239903e-34, -1e-10);
%!shared kw, kwg
%! kw = [2.1, 3.4, 1.9, 5.6, 4.4, 3.8, 2.7, 6.1, 3.3, 4.9, 1.2, 2.8, 0.9, ...
%!       2.2, 3.1, 1.7, 2.5, 0.4, 1.9, 3.6, 2.0, 1.1, 1.8, 2.9, 2.2, 4.1, ...
%!       4.9, 2.6, 3.0, 4.8]';
%! kwg = [ones(10, 1); 2 * ones(12, 1); 3 * ones(8, 1)];
## Effect size against the R package effectsize 1.0.3, rank_eta_squared
%!test
%! [~, tbl, stats] = kruskalwallis (kw, kwg, 'off');
%! assert_equal (stats.EtaSquared, 0.29653567934214958, -1e-13);
%! ## The bounds put H at the upper and lower 2.5% points of the noncentral
%! ## chi-square distribution
%! lambda = stats.EtaSquaredCI * tbl{3,3};
%! assert_equal (ncx2cdf (tbl{2,5}, 2, lambda), [0.975, 0.025], -1e-10);
%!test
%! [~, tbl, stats] = kruskalwallis (kw, kwg, 'off', 'Alpha', 0.1);
%! lambda = stats.EtaSquaredCI * tbl{3,3};
%! assert_equal (ncx2cdf (tbl{2,5}, 2, lambda), [0.95, 0.05], -1e-10);
%!test
%! x = [1, 2, 2, 3, 3, 3, 4, 5, 2, 2, 3, 3, 4, 6, 5, 6, 6, 7]';
%! g = [ones(8, 1); 2 * ones(6, 1); 3 * ones(4, 1)];
%! [~, ~, stats] = kruskalwallis (x, g, 'off');
%! assert_equal (stats.EtaSquared, 0.38833809693938293, -1e-13);
%! assert_equal (stats.EtaSquaredCI(2), 1);
%!test
%! [~, ~, stats] = kruskalwallis ([1, 5, 2, 8, 4, 3, 7, 6, 2, 9, 1, 3]', ...
%!                                [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]', 'off');
%! assert_equal (stats.EtaSquared, 0);
%! assert_equal (stats.EtaSquaredCI(1), 0);
%!test
%! rand ('state', 1);
%! [~, ~, stats] = kruskalwallis (kw, kwg, 'off', 'ConfidenceIntervalType', ...
%!                                'bootstrap', 'NumBootstraps', 200);
%! ci = stats.EtaSquaredCI;
%! assert_equal ([ci(1) <= ci(2), ci(1) >= 0, ci(2) <= 1], true (1, 3));

%!error<kruskalwallis: invalid optional paired argument.> ...
%! kruskalwallis ([1, 2, 3, 4], [1, 1, 2, 2], 'off', 'Foo', 1)
%!error<kruskalwallis: 'Alpha' must be a scalar between 0 and 1.> ...
%! kruskalwallis ([1, 2, 3, 4], [1, 1, 2, 2], 'off', 'Alpha', 1)
%!error<kruskalwallis: 'ConfidenceIntervalType' must be 'exact' or 'bootstrap'.> ...
%! kruskalwallis ([1, 2, 3, 4], [1, 1, 2, 2], 'off', ...
%!                'ConfidenceIntervalType', 'none')
%!error<kruskalwallis: 'NumBootstraps' must be a positive integer.> ...
%! kruskalwallis ([1, 2, 3, 4], [1, 1, 2, 2], 'off', 'NumBootstraps', 0)
