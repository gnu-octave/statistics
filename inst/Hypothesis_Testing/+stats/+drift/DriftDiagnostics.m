## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftp {statistics} stats.drift.DriftDiagnostics
##
## The drift found between two data sets, as @code{detectdrift} returns it.
##
## A @code{stats.drift.DriftDiagnostics} object holds, for each variable
## compared, the metric measuring how far the target data has moved from the
## baseline data, the p-value of the permutation test of that change with its
## confidence interval, and the drift status the interval gives, together with
## the drift status of the data as a whole.  Every property is read-only.
## @code{detectdrift} is the documented way to create one.
##
## @code{summary} tabulates the results, @code{ecdf} and @code{histcounts}
## return the distributions compared, and @code{plotDriftStatus},
## @code{plotEmpiricalCDF}, @code{plotHistogram} and
## @code{plotPermutationResults} draw them.
##
## @seealso{detectdrift}
## @end deftp

classdef DriftDiagnostics

  properties (GetAccess = public, SetAccess = private)

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} Baseline
    ##
    ## The baseline data, as given to @code{detectdrift}.
    ##
    ## @end deftp
    Baseline = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} ConfidenceIntervals
    ##
    ## The 95% Clopper-Pearson confidence interval of each p-value, a two-row
    ## matrix holding the lower bounds over the upper, one column per
    ## variable, @code{NaN} where p-values were not estimated.
    ##
    ## @end deftp
    ConfidenceIntervals = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} DriftStatus
    ##
    ## The drift status of each variable, a string array holding
    ## @qcode{"Drift"}, @qcode{"Warning"} or @qcode{"Stable"}, missing where
    ## p-values were not estimated.
    ##
    ## @end deftp
    DriftStatus = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} Metrics
    ##
    ## The metric each variable was measured by, a string array.
    ##
    ## @end deftp
    Metrics = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} MetricValues
    ##
    ## The metric of each variable between the baseline and the target data,
    ## a row vector.
    ##
    ## @end deftp
    MetricValues = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} MultipleTestCorrection
    ##
    ## The correction for testing every variable, @qcode{"Bonferroni"} or
    ## @qcode{"FalseDiscoveryRate"}.
    ##
    ## @end deftp
    MultipleTestCorrection = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} MultipleTestDriftStatus
    ##
    ## The drift status of the data as a whole, from the corrected p-values,
    ## missing where p-values were not estimated.
    ##
    ## @end deftp
    MultipleTestDriftStatus = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} NumVariables
    ##
    ## The number of variables compared.
    ##
    ## @end deftp
    NumVariables = 0;

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} PValues
    ##
    ## The p-value of each variable, a row vector, @code{NaN} where p-values
    ## were not estimated.
    ##
    ## @end deftp
    PValues = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} Target
    ##
    ## The target data, as given to @code{detectdrift}.
    ##
    ## @end deftp
    Target = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} VariableNames
    ##
    ## The names of the variables compared, a string array, @qcode{"x1"},
    ## @qcode{"x2"}, @dots{} where the data came as arrays.
    ##
    ## @end deftp
    VariableNames = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} CategoricalVariables
    ##
    ## The indices of the variables holding levels, among those compared.
    ##
    ## @end deftp
    CategoricalVariables = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} NumPermutations
    ##
    ## The number of permutations drawn for each variable, a row vector, 1
    ## where p-values were not estimated.
    ##
    ## @end deftp
    NumPermutations = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} PermutationResults
    ##
    ## A table with one row per variable, named by it, holding in its one
    ## variable @code{PermutationResults} a cell with the metric over every
    ## permutation, the observed arrangement first.
    ##
    ## @end deftp
    PermutationResults = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} DriftThreshold
    ##
    ## The threshold below which a p-value means drift.
    ##
    ## @end deftp
    DriftThreshold = [];

    ## -*- texinfo -*-
    ## @deftp {stats.drift.DriftDiagnostics} {property} WarningThreshold
    ##
    ## The threshold below which a p-value means a warning.
    ##
    ## @end deftp
    WarningThreshold = [];

  endproperties

  properties (GetAccess = public, SetAccess = private, Hidden)
    Data_ = [];         # the pooled data, baseline rows first, levels coded
    NumBaseline_ = 0;   # how many rows of Data_ are the baseline's
    Labels_ = {};       # the name of each code of a variable holding levels
    Estimated_ = false; # whether p-values were estimated
  endproperties

  methods (Hidden)

    ## -*- texinfo -*-
    ## @deftypefn {stats.drift.DriftDiagnostics} {@var{obj} =} stats.drift.DriftDiagnostics (@var{S})
    ##
    ## Create a @code{stats.drift.DriftDiagnostics} object from the structure
    ## @var{S} @code{detectdrift} fills, holding every property and the pooled
    ## data the methods read.  The documented way to reach this constructor
    ## is @code{detectdrift}.
    ##
    ## @seealso{detectdrift}
    ## @end deftypefn
    function this = DriftDiagnostics (S)
      if (nargin == 0)
        return;
      endif
      K = numel (S.VariableNames);
      this.Baseline = S.Baseline;
      this.Target = S.Target;
      this.VariableNames = string (S.VariableNames(:)');
      this.CategoricalVariables = S.CategoricalVariables(:)';
      if (isempty (this.CategoricalVariables))
        this.CategoricalVariables = [];
      endif
      this.Metrics = string (S.Metrics(:)');
      this.MetricValues = S.MetricValues(:)';
      this.PValues = S.PValues(:)';
      this.ConfidenceIntervals = S.ConfidenceIntervals;
      this.NumPermutations = S.NumPermutations(:)';
      this.NumVariables = K;
      this.PermutationResults = table (S.PermutationResults(:), ...
                                       'VariableNames', ...
                                       {'PermutationResults'}, ...
                                       'RowNames', S.VariableNames(:));
      if (S.Estimated)
        this.DriftStatus = string (S.DriftStatus(:)');
        this.MultipleTestDriftStatus = string (S.MultipleTestDriftStatus);
      else
        this.DriftStatus = repmat (string (missing), 1, K);
        this.MultipleTestDriftStatus = string (missing);
      endif
      if (strcmp (S.MultipleTestCorrection, 'bonferroni'))
        this.MultipleTestCorrection = string ('Bonferroni');
      else
        this.MultipleTestCorrection = string ('FalseDiscoveryRate');
      endif
      this.DriftThreshold = S.DriftThreshold;
      this.WarningThreshold = S.WarningThreshold;
      this.Data_ = S.Data;
      this.NumBaseline_ = S.NumBaseline;
      this.Labels_ = S.Labels;
      this.Estimated_ = S.Estimated;
    endfunction

    function disp (this)
      printf ('\n  DriftDiagnostics\n\n');
      printf ('%27s: %s\n', 'VariableNames', ddJoin (this.VariableNames));
      printf ('%27s: %s\n', 'CategoricalVariables', ...
              mat2str (this.CategoricalVariables));
      if (this.Estimated_)
        printf ('%27s: %s\n', 'DriftStatus', ddJoin (this.DriftStatus));
        printf ('%27s: %s\n', 'PValues', mat2str (this.PValues, 4));
        printf ('%27s: [2x%d double]\n', 'ConfidenceIntervals', ...
                this.NumVariables);
        printf ('%27s: "%s"\n', 'MultipleTestDriftStatus', ...
                char (this.MultipleTestDriftStatus));
        printf ('%27s: %g\n', 'DriftThreshold', this.DriftThreshold);
        printf ('%27s: %g\n', 'WarningThreshold', this.WarningThreshold);
      else
        printf ('%27s: %s\n', 'Metrics', ddJoin (this.Metrics));
        printf ('%27s: %s\n', 'MetricValues', mat2str (this.MetricValues, 4));
      endif
      printf ('\n');
    endfunction

    function display (this)
      disp (this);
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {stats.drift.DriftDiagnostics} {} summary (@var{obj})
    ## @deftypefnx {stats.drift.DriftDiagnostics} {@var{tbl} =} summary (@var{obj})
    ##
    ## Tabulate the drift found.
    ##
    ## Where p-values were estimated, @code{summary} prints the drift status
    ## of the data as a whole and a table with one row per variable holding
    ## @code{DriftStatus}, @code{PValue} and @code{ConfidenceInterval}.
    ## @var{tbl} is that table with a last row, @qcode{MultipleTest}, holding
    ## the status of the data as a whole and @code{NaN} for its p-value and
    ## interval.  Where they were not, the table holds @code{MetricValue} and
    ## @code{Metric}.
    ##
    ## @seealso{detectdrift}
    ## @end deftypefn
    function tbl = summary (this)
      names = cellstr (this.VariableNames(:));
      if (! this.Estimated_)
        T = table (this.MetricValues(:), this.Metrics(:), ...
                   'VariableNames', {'MetricValue', 'Metric'}, ...
                   'RowNames', names);
      else
        T = table (this.DriftStatus(:), this.PValues(:), ...
                   this.ConfidenceIntervals', 'VariableNames', ...
                   {'DriftStatus', 'PValue', 'ConfidenceInterval'}, ...
                   'RowNames', names);
      endif
      if (nargout == 0)
        if (this.Estimated_)
          printf ('\n  Multiple Test Correction Drift Status: %s\n\n', ...
                  char (this.MultipleTestDriftStatus));
        endif
        disp (T);
        return;
      endif
      if (this.Estimated_)
        T = table ([this.DriftStatus(:); this.MultipleTestDriftStatus], ...
                   [this.PValues(:); NaN], ...
                   [this.ConfidenceIntervals'; NaN, NaN], 'VariableNames', ...
                   {'DriftStatus', 'PValue', 'ConfidenceInterval'}, ...
                   'RowNames', [names; {'MultipleTest'}]);
      endif
      tbl = T;
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {stats.drift.DriftDiagnostics} {@var{tbl} =} ecdf (@var{obj})
    ## @deftypefnx {stats.drift.DriftDiagnostics} {@var{tbl} =} ecdf (@var{obj}, @qcode{'Variable'}, @var{var})
    ##
    ## The empirical distribution functions of the baseline and the target.
    ##
    ## @var{tbl} holds one row per variable, named by it, and three variables
    ## of cells: @code{x}, the pooled values sorted with the first repeated,
    ## and @code{F_Baseline} and @code{F_Target}, the two distribution
    ## functions there, starting from 0.  Where values tie, only the last copy
    ## takes the full value, so that the points draw the two functions as
    ## steps.  A variable holding levels has @code{NaN} in each cell.
    ## @var{var} names the variables, by name or index, all of them by
    ## default.
    ##
    ## @seealso{stats.drift.DriftDiagnostics.plotEmpiricalCDF}
    ## @end deftypefn
    function tbl = ecdf (this, varargin)
      idx = ddVariables (this, varargin, 'ecdf');
      n = numel (idx);
      xs = cell (n, 1);
      fb = cell (n, 1);
      ft = cell (n, 1);
      for i = 1:n
        j = idx(i);
        if (any (this.CategoricalVariables == j))
          [xs{i}, fb{i}, ft{i}] = deal (NaN);
        else
          [xs{i}, fb{i}, ft{i}] = ddSteps (this, j);
        endif
      endfor
      tbl = table (xs, fb, ft, 'VariableNames', ...
                   {'x', 'F_Baseline', 'F_Target'}, ...
                   'RowNames', cellstr (this.VariableNames(idx)(:)));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {stats.drift.DriftDiagnostics} {@var{tbl} =} histcounts (@var{obj})
    ## @deftypefnx {stats.drift.DriftDiagnostics} {@var{tbl} =} histcounts (@var{obj}, @qcode{'Variable'}, @var{var})
    ##
    ## The histograms of the baseline and the target, in percent.
    ##
    ## @var{tbl} holds one row per variable, named by it, and three variables
    ## of cells: @code{Bins}, and @code{Counts_Baseline} and
    ## @code{Counts_Target}, each the percentage of its sample in each bin.
    ## For a continuous variable @code{Bins} holds the bin edges
    ## @code{histcounts} chooses for the two samples pooled.  For a variable
    ## holding levels it holds the levels as a @code{categorical} array, and
    ## each count is increased by 0.5 before the percentages are taken, as the
    ## metrics take them.  @var{var} names the variables, by name or index,
    ## all of them by default.
    ##
    ## @seealso{stats.drift.DriftDiagnostics.plotHistogram}
    ## @end deftypefn
    function tbl = histcounts (this, varargin)
      idx = ddVariables (this, varargin, 'histcounts');
      n = numel (idx);
      bins = cell (n, 1);
      cb = cell (n, 1);
      ct = cell (n, 1);
      for i = 1:n
        [bins{i}, cb{i}, ct{i}] = ddCounts (this, idx(i));
      endfor
      tbl = table (bins, cb, ct, 'VariableNames', ...
                   {'Bins', 'Counts_Baseline', 'Counts_Target'}, ...
                   'RowNames', cellstr (this.VariableNames(idx)(:)));
    endfunction


    ## -*- texinfo -*-
    ## @deftypefn  {stats.drift.DriftDiagnostics} {} plotDriftStatus (@var{obj})
    ## @deftypefnx {stats.drift.DriftDiagnostics} {@var{h} =} plotDriftStatus (@var{obj})
    ##
    ## Plot the p-value and confidence interval of every variable.
    ##
    ## Each variable is drawn as its p-value with its confidence interval as
    ## horizontal error bars, one set of bars for each drift status, the
    ## variables in alphabetical order up the vertical axis, with a line at
    ## each threshold.  @var{h} holds the three sets of bars, for
    ## @qcode{"Stable"}, @qcode{"Warning"} and @qcode{"Drift"} in turn, each
    ## spanning every variable with @code{NaN} where a variable has another
    ## status.  MATLAB draws on a categorical axis, which Octave does not have,
    ## so the variables sit at 1, 2, @dots{} with their names as tick labels.
    ##
    ## @seealso{stats.drift.DriftDiagnostics.summary}
    ## @end deftypefn
    function h = plotDriftStatus (this)
      if (! this.Estimated_)
        error (strcat ("stats.drift.DriftDiagnostics.plotDriftStatus:", ...
                       " p-values were not estimated."));
      endif
      [names, o] = sort (cellstr (this.VariableNames));
      K = numel (names);
      p = this.PValues(o);
      lo = p - this.ConfidenceIntervals(1,o);
      hi = this.ConfidenceIntervals(2,o) - p;
      st = cellstr (this.DriftStatus(o));
      groups = {'Stable', 'Warning', 'Drift'};
      colors = ddStatusColors ();
      hax = newplot ();
      hold (hax, 'on');
      hh = zeros (3, 1);
      for g = 1:3
        x = p;
        x(! strcmp (st, groups{g})) = NaN;
        hh(g) = errorbar (hax, x, 1:K, lo, hi, '>');
        set (hh(g), 'color', colors(g,:), 'linestyle', 'none', ...
                    'marker', 'o');
      endfor
      xline (hax, this.WarningThreshold, ':', 'color', colors(2,:));
      xline (hax, this.DriftThreshold, ':', 'color', colors(3,:));
      hold (hax, 'off');
      set (hax, 'xlim', [0, 1], 'xtick', 0:0.1:1, 'ylim', [0.5, K + 0.5], ...
                'ytick', 1:K, 'yticklabel', names);
      title (hax, 'Estimated P-Values and Confidence Intervals');
      xlabel (hax, 'P-Values');
      ylabel (hax, 'Variable Names');
      legend (hax, [groups, {'Warning Threshold', 'Drift Threshold'}]);
      if (nargout > 0)
        h = hh;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {stats.drift.DriftDiagnostics} {} plotEmpiricalCDF (@var{obj})
    ## @deftypefnx {stats.drift.DriftDiagnostics} {} plotEmpiricalCDF (@var{obj}, @qcode{'Variable'}, @var{var})
    ## @deftypefnx {stats.drift.DriftDiagnostics} {@var{h} =} plotEmpiricalCDF (@dots{})
    ##
    ## Plot the empirical distribution functions of one continuous variable.
    ##
    ## The baseline and the target are drawn as steps, from the points
    ## @code{ecdf} returns.  @var{var} names the variable, by name or index;
    ## by default it is the one with the smallest p-value, or the first where
    ## p-values were not estimated.  @var{h} holds the two lines, baseline
    ## first.
    ##
    ## @seealso{stats.drift.DriftDiagnostics.ecdf}
    ## @end deftypefn
    function h = plotEmpiricalCDF (this, varargin)
      j = ddOneVariable (this, varargin, 'plotEmpiricalCDF');
      if (any (this.CategoricalVariables == j))
        error (strcat ("stats.drift.DriftDiagnostics.plotEmpiricalCDF:", ...
                       " 'Variable' must be a continuous variable."));
      endif
      [x, fb, ft] = ddSteps (this, j);
      colors = ddSampleColors ();
      name = char (this.VariableNames(j));
      hax = newplot ();
      hh = stairs (hax, [x, x], [fb, ft]);
      set (hh(1), 'color', colors(1,:));
      set (hh(2), 'color', colors(2,:));
      title (hax, ['ECDF for ', name]);
      xlabel (hax, name);
      ylabel (hax, 'Cumulative Probability');
      legend (hax, {'Baseline', 'Target'});
      if (nargout > 0)
        h = hh(:);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {stats.drift.DriftDiagnostics} {} plotHistogram (@var{obj})
    ## @deftypefnx {stats.drift.DriftDiagnostics} {} plotHistogram (@var{obj}, @qcode{'Variable'}, @var{var})
    ## @deftypefnx {stats.drift.DriftDiagnostics} {@var{h} =} plotHistogram (@dots{})
    ##
    ## Plot the histograms of one variable in the baseline and the target.
    ##
    ## The percentages @code{histcounts} returns are drawn as bars side by
    ## side, over the bins of a continuous variable or the levels of one
    ## holding levels.  @var{var} names the variable, by name or index; by
    ## default it is the one with the smallest p-value, or the first where
    ## p-values were not estimated.  @var{h} holds the two sets of bars,
    ## baseline first.
    ##
    ## @seealso{stats.drift.DriftDiagnostics.histcounts}
    ## @end deftypefn
    function h = plotHistogram (this, varargin)
      j = ddOneVariable (this, varargin, 'plotHistogram');
      [bins, cb, ct] = ddCounts (this, j);
      colors = ddSampleColors ();
      name = char (this.VariableNames(j));
      hax = newplot ();
      if (iscategorical (bins))
        x = 1:numel (bins);
      else
        x = (bins(1:end-1) + bins(2:end)) / 2;
      endif
      hh = bar (hax, x, [cb(:), ct(:)]);
      set (hh(1), 'facecolor', colors(1,:));
      set (hh(2), 'facecolor', colors(2,:));
      if (iscategorical (bins))
        set (hax, 'xtick', x, 'xticklabel', cellstr (bins));
      endif
      title (hax, ['Histogram for ', name]);
      xlabel (hax, [name, ' Bins']);
      ylabel (hax, 'Distribution (%)');
      legend (hax, {'Baseline', 'Target'});
      if (nargout > 0)
        h = hh(:);
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {stats.drift.DriftDiagnostics} {} plotPermutationResults (@var{obj})
    ## @deftypefnx {stats.drift.DriftDiagnostics} {} plotPermutationResults (@var{obj}, @qcode{'Variable'}, @var{var})
    ## @deftypefnx {stats.drift.DriftDiagnostics} {@var{h} =} plotPermutationResults (@dots{})
    ##
    ## Plot the metric of one variable over its permutations.
    ##
    ## The percentage of permutations in each bin is drawn as bars, split at
    ## the observed metric, which is marked by a line: the bars below it and
    ## those at or above it, whose share is the p-value.  @var{var} names the
    ## variable, by name or index; by default it is the one with the smallest
    ## p-value.  @var{h} holds the two sets of bars, those below first.
    ## MATLAB draws them as @code{histogram} objects, which Octave does not
    ## have.
    ##
    ## @seealso{detectdrift}
    ## @end deftypefn
    function h = plotPermutationResults (this, varargin)
      if (! this.Estimated_)
        error (strcat ("stats.drift.DriftDiagnostics.", ...
                       "plotPermutationResults: p-values were not", ...
                       " estimated."));
      endif
      j = ddOneVariable (this, varargin, 'plotPermutationResults');
      v = this.PermutationResults.PermutationResults{j};
      obs = this.MetricValues(j);
      N = numel (v);
      ## One bin width for both sides, a whole number of bins below the
      ## observed metric and as many above as reach the largest value
      [~, e] = histcounts (v);
      w = e(2) - e(1);
      lo = min (v);
      hi = max (v);
      if (obs > lo)
        nb = ceil ((obs - lo) / w);
        w = (obs - lo) / nb;
        eb = linspace (lo, obs, nb + 1);
      else
        eb = [];
      endif
      ea = obs:w:hi;
      if (isempty (ea) || ea(end) < hi)
        ea(end+1) = hi;
      endif
      if (numel (ea) < 2)
        ea = [obs, obs + w];
      endif
      below = [];
      xb = [];
      if (! isempty (eb))
        below = 100 * histcounts (v(v < obs), eb) / N;
        xb = (eb(1:end-1) + eb(2:end)) / 2;
      endif
      above = 100 * histcounts (v(v >= obs), ea) / N;
      xa = (ea(1:end-1) + ea(2:end)) / 2;
      colors = ddSampleColors ();
      hax = newplot ();
      hold (hax, 'on');
      hh = zeros (2, 1);
      if (isempty (xb))
        hh(1) = bar (hax, NaN, NaN, 1, 'facecolor', colors(1,:));
      else
        hh(1) = bar (hax, xb, below, 1, 'facecolor', colors(1,:));
      endif
      hh(2) = bar (hax, xa, above, 1, 'facecolor', colors(2,:));
      xline (hax, obs, ':', 'color', [0.15, 0.15, 0.15]);
      hold (hax, 'off');
      metric = char (this.Metrics(j));
      title (hax, ['Permutation Results for ', ...
                   char(this.VariableNames(j))]);
      xlabel (hax, [metric, ' Metric Values']);
      ylabel (hax, 'Distribution (%)');
      legend (hax, {sprintf('< %g', obs), sprintf('\\geq %g', obs)});
      if (nargout > 0)
        h = hh;
      endif
    endfunction

  endmethods

endclassdef

## The variable names joined for display.
function s = ddJoin (names)
  s = ['[', strjoin(strcat ('"', cellstr (names), '"'), '  '), ']'];
endfunction

## The indices of the variables a method is asked about: 'Variable' naming
## them by name or index, all of them where it is absent.
function idx = ddVariables (this, args, caller)
  K = this.NumVariables;
  if (isempty (args))
    idx = 1:K;
    return;
  endif
  if (numel (args) != 2 || ! (ischar (args{1}) || isa (args{1}, 'string'))
      || ! strcmpi (args{1}, 'Variable'))
    error (strcat ("stats.drift.DriftDiagnostics.%s: invalid optional", ...
                   " paired argument."), caller);
  endif
  v = args{2};
  names = cellstr (this.VariableNames);
  if (isa (v, 'string') || ischar (v) || iscellstr (v))
    v = cellstr (v);
    [ok, idx] = ismember (v, names);
    if (! all (ok))
      error (strcat ("stats.drift.DriftDiagnostics.%s: 'Variable' must", ...
                     " name variables compared."), caller);
    endif
  elseif (islogical (v) && numel (v) == K)
    idx = find (v);
  elseif (isnumeric (v) && isreal (v) && all (v(:) == fix (v(:)))
          && all (v(:) >= 1) && all (v(:) <= K))
    idx = v;
  else
    error (strcat ("stats.drift.DriftDiagnostics.%s: 'Variable' must be", ...
                   " names or indices of variables compared."), caller);
  endif
  idx = idx(:)';
endfunction

## The staircase points of one continuous variable.
function [x, fb, ft] = ddSteps (this, j)
  z = this.Data_(:,j);
  mx = this.NumBaseline_;
  base = [true(mx, 1); false(numel (z) - mx, 1)];
  [s, o] = sort (z);
  b = base(o);
  N = numel (s);
  ## Within a run of ties, every copy but the last takes the count before
  ## the run, and the last the count at its end
  last = [s(1:end-1) != s(2:end); true];
  first = [true; s(2:end) != s(1:end-1)];
  start = cummax (first .* (1:N)');
  cb = [0; cumsum(b)];
  ct = [0; cumsum(! b)];
  fb = cb(start);
  ft = ct(start);
  fb(last) = cb(find (last) + 1);
  ft(last) = ct(find (last) + 1);
  fb = fb / sum (b);
  ft = ft / sum (! b);
  x = [s(1); s];
  fb = [0; fb];
  ft = [0; ft];
endfunction

## The bins and the percentage of each sample in each, for one variable.
function [bins, cb, ct] = ddCounts (this, j)
  z = this.Data_(:,j);
  mx = this.NumBaseline_;
  if (any (this.CategoricalVariables == j))
    lab = this.Labels_{j};
    L = numel (lab);
    nb = accumarray (z(1:mx), 1, [L, 1])' + 0.5;
    nt = accumarray (z(mx+1:end), 1, [L, 1])' + 0.5;
    bins = categorical (lab(:), lab(:));
    cb = 100 * nb / sum (nb);
    ct = 100 * nt / sum (nt);
  else
    [~, bins] = histcounts (z);
    cb = 100 * histcounts (z(1:mx), bins) / mx;
    ct = 100 * histcounts (z(mx+1:end), bins) / (numel (z) - mx);
  endif
endfunction

## The variable a plot is of: the one named, or the one with the smallest
## p-value, the first where p-values were not estimated.
function j = ddOneVariable (this, args, caller)
  if (isempty (args))
    if (this.Estimated_)
      [~, j] = min (this.PValues);
    else
      j = 1;
    endif
    return;
  endif
  j = ddVariables (this, args, caller);
  if (numel (j) != 1)
    error (strcat ("stats.drift.DriftDiagnostics.%s: 'Variable' must", ...
                   " name one variable."), caller);
  endif
endfunction

## The colours of the three drift statuses, Stable, Warning and Drift.
function c = ddStatusColors ()
  c = [0, 0.4470, 0.7410; 0.9290, 0.6940, 0.1250; 0.8500, 0.3250, 0.0980];
endfunction

## The colours of the baseline and the target.
function c = ddSampleColors ()
  c = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980];
endfunction

## Expected values from MATLAB R2024a
%!shared ddE, ddD
%! ddE = detectdrift ([1; 2; 2; 4], [3; 5], 'EstimatePValues', false);
%! rand ('seed', 1);
%! x = (1:60)' / 10;
%! ddD = detectdrift (table (x, x, 'VariableNames', {'Same', 'Far'}), ...
%!                    table (x([2:60, 1]), x + 10, ...
%!                           'VariableNames', {'Same', 'Far'}));
%!test
%! e = ecdf (ddE);
%! assert_equal (e.x{1}', [1, 1, 2, 2, 3, 4, 5]);
%! assert_equal (e.F_Baseline{1}', [0, 0.25, 0.25, 0.75, 0.75, 1, 1]);
%! assert_equal (e.F_Target{1}', [0, 0, 0, 0, 0.5, 0.5, 1]);
%!test
%! h = histcounts (ddE);
%! assert_equal (h.Bins{1}, 0.5:5.5);
%! assert_equal (h.Counts_Baseline{1}, [25, 50, 0, 25, 0]);
%! assert_equal (h.Counts_Target{1}, [0, 0, 50, 0, 50]);
%!test
%! ## Levels: each count increased by 0.5
%! C = detectdrift (categorical ([1; 1; 2; 3]), categorical ([1; 2; 2; 2]), ...
%!                  'EstimatePValues', false);
%! h = histcounts (C);
%! assert_equal (cellstr (h.Bins{1}), {'1'; '2'; '3'});
%! assert_equal (h.Counts_Baseline{1}, [2.5, 1.5, 1.5] * 100 / 5.5, -1e-14);
%! assert_equal (C.VariableNames, string ('x1'));
%! assert_equal (C.CategoricalVariables, 1);
%!test
%! T = table ([1; 2; 3], categorical ([1; 2; 1]), 'VariableNames', {'A', 'B'});
%! e = ecdf (detectdrift (T, T, 'EstimatePValues', false), 'Variable', 'B');
%! assert_equal (e.x{1}, NaN);
%! assert_equal (e.Properties.RowNames, {'B'});
%!test
%! e = ecdf (ddD, 'Variable', 2);
%! assert_equal (e.Properties.RowNames, {'Far'});
%!test
%! t = summary (ddE);
%! assert_equal (t.Properties.VariableNames, {'MetricValue', 'Metric'});
%! assert_equal (t.MetricValue, 1.75, -1e-14);
%!test
%! t = summary (ddD);
%! assert_equal (t.Properties.RowNames, {'Same'; 'Far'; 'MultipleTest'});
%! assert_equal (t.DriftStatus, string ({'Stable'; 'Drift'; 'Drift'}));
%! assert_equal (t.PValue, [1; 0.001; NaN]);
%! assert_equal (size (t.ConfidenceInterval), [3, 2]);
%!assert_equal (ddD.NumVariables, 2)
%!assert_equal (ddD.WarningThreshold, 0.1)
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = plotDriftStatus (ddD);
%!   assert_equal (numel (h), 3);
%!   ## Alphabetical, so Far comes before Same
%!   assert_equal (get (gca (), 'yticklabel'), {'Far'; 'Same'});
%!   assert_equal (get (h(3), 'xdata')(:)', [0.001, NaN]);
%!   assert_equal (get (get (gca (), 'title'), 'string'), ...
%!                 'Estimated P-Values and Confidence Intervals');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = plotEmpiricalCDF (ddE);
%!   assert_equal (numel (h), 2);
%!   assert_equal (get (get (gca (), 'title'), 'string'), 'ECDF for x1');
%!   assert_equal (get (get (gca (), 'ylabel'), 'string'), ...
%!                 'Cumulative Probability');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = plotHistogram (ddE);
%!   assert_equal (get (h(1), 'xdata')(:)', 1:5);
%!   assert_equal (get (h(1), 'ydata')(:)', [25, 50, 0, 25, 0]);
%!   assert_equal (get (get (gca (), 'xlabel'), 'string'), 'x1 Bins');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! ## The variable with the smallest p-value by default
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   plotHistogram (ddD);
%!   assert_equal (get (get (gca (), 'title'), 'string'), 'Histogram for Far');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
%!test
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = plotPermutationResults (ddD, 'Variable', 'Same');
%!   assert_equal (numel (h), 2);
%!   assert_equal (sum (get (h(2), 'ydata')(:)), 100, -1e-12);
%!   assert_equal (get (get (gca (), 'xlabel'), 'string'), ...
%!                 'Wasserstein Metric Values');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!error<stats.drift.DriftDiagnostics.ecdf: 'Variable' must name variables compared.> ...
%! ecdf (ddE, 'Variable', 'nope')
%!error<stats.drift.DriftDiagnostics.ecdf: invalid optional paired argument.> ...
%! ecdf (ddE, 'Foo', 1)
%!error<stats.drift.DriftDiagnostics.histcounts: 'Variable' must be names or indices of variables compared.> ...
%! histcounts (ddE, 'Variable', 3)
%!error<stats.drift.DriftDiagnostics.plotDriftStatus: p-values were not estimated.> ...
%! plotDriftStatus (ddE)
%!error<stats.drift.DriftDiagnostics.plotPermutationResults: p-values were not estimated.> ...
%! plotPermutationResults (ddE)
%!error<stats.drift.DriftDiagnostics.plotEmpiricalCDF: 'Variable' must be a continuous variable.> ...
%! plotEmpiricalCDF (detectdrift (categorical ([1; 2]), categorical ([2; 2]), ...
%!                                'EstimatePValues', false))
%!error<stats.drift.DriftDiagnostics.plotHistogram: 'Variable' must name one variable.> ...
%! plotHistogram (ddD, 'Variable', [1, 2])
