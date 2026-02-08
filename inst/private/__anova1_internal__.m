function [p, tbl, stats] = __anova1_internal__ (x, group)
  % __anova1_internal__  Internal numeric engine for one-way ANOVA
  %
  %   [p, tbl, stats] = __anova1_internal__ (x, group)
  %
  % Inputs:
  %   x     = numeric vector of data
  %   group = grouping variable (vector, cellstr, or categorical)
  %
  % Outputs:
  %   p     = p-value of one-way ANOVA
  %   tbl   = ANOVA table (cell array)
  %   stats = structure for multcompare()

  % ---------------------------------------------------
  % Input validation
  % ---------------------------------------------------
  if (nargin < 1)
    error ("__anova1_internal__: x is required");
  end

  if (isempty (group))
    % Auto-generate a single group if none supplied
    group = ones (numel (x), 1);
  end

  x = x(:);
  group = group(:);

  if (numel (x) ~= numel (group))
    error ("__anova1_internal__: x and group must have same length");
  end

  % Convert grouping variable to numeric indices
  [gnames, ~, idx] = unique (group);
  k = numel (gnames);    % number of groups
  n = numel (x);         % total samples

  % ---------------------------------------------------
  % Group-level statistics
  % ---------------------------------------------------
  group_means = accumarray (idx, x, [], @mean);
  group_counts = accumarray (idx, 1);
  grand_mean = mean (x);

  % ---------------------------------------------------
  % Sum of Squares
  % ---------------------------------------------------
  ss_between = sum (group_counts .* (group_means - grand_mean).^2);
  ss_within  = sum ((x - group_means(idx)).^2);

  % Degrees of freedom
  df_between = k - 1;
  df_within  = n - k;

  % Mean squares
  ms_between = ss_between / df_between;
  ms_within  = ss_within / df_within;

  % ---------------------------------------------------
  % F-statistic and p-value
  % ---------------------------------------------------
  F = ms_between / ms_within;
  p = 1 - fcdf (F, df_between, df_within);

  % ---------------------------------------------------
  % ANOVA table (MATLAB style)
  % ---------------------------------------------------
  tbl = {
    'Source',    'SS',          'df',          'MS',          'F',     'Prob>F';
    'Groups',    ss_between,    df_between,    ms_between,    F,      p;
    'Error',     ss_within,     df_within,     ms_within,     NaN,    NaN;
    'Total',     ss_between + ss_within,   n - 1,   NaN,     NaN,    NaN
  };

  % ---------------------------------------------------
  % Stats struct for multcompare()
  % ---------------------------------------------------
  stats.gnames = cellstr (string (gnames));   % group names
  stats.n = group_counts;                     % sample size per group
  stats.means = group_means;                  % means of each group
  stats.df = df_within;                       % degrees of freedom
  stats.s = sqrt (ms_within);                 % pooled standard deviation
  stats.source = 'anova1';                    % metadata

end
