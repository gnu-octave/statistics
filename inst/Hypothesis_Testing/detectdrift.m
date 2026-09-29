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
## @deftypefn  {statistics} {@var{DDiagnostics} =} detectdrift (@var{Baseline}, @var{Target})
## @deftypefnx {statistics} {@var{DDiagnostics} =} detectdrift (@var{Baseline}, @var{Target}, @var{Name}, @var{Value})
##
## Detect drift between two data sets, one variable at a time.
##
## @code{@var{DDiagnostics} = detectdrift (@var{Baseline}, @var{Target})}
## measures, for each variable, how far the distribution in @var{Target} has
## moved from the one in @var{Baseline}, tests the change by permutation, and
## returns the results as a @code{stats.drift.DriftDiagnostics} object.
## @var{Baseline} and @var{Target} are both numeric matrices with the same
## number of columns, both @code{categorical} arrays, or both tables, read
## over the variables they share, which must be all the variables of one of
## them.  An observation holding a missing value in a variable used is left
## out.
##
## Each variable is measured by a metric: a continuous one by
## @qcode{'ContinuousMetric'}, one holding levels by
## @qcode{'CategoricalMetric'}.  With @math{F} and @math{G} the empirical
## distribution functions of the two samples and @math{H} that of the two
## pooled, the metrics are:
##
## @multitable @columnfractions 0.24 0.02 0.74
## @headitem Metric @tab @tab Definition
## @item @qcode{'wasserstein'} @tab @tab The area between @math{F} and
## @math{G}, the default for a continuous variable.
## @item @qcode{'ks'} @tab @tab The largest difference between @math{F}
## and @math{G}.
## @item @qcode{'ad'} @tab @tab The Anderson-Darling statistic, the sum over
## the distinct pooled values but the largest of
## @math{(F - G)^2 / (H (1 - H))}, divided by the number of pooled
## observations.
## @item @qcode{'energy'} @tab @tab The energy distance,
## @math{sqrt (2 E|x - y| - E|x - x'| - E|y - y'|)}.
## @item @qcode{'hellinger'} @tab @tab @math{sqrt (1 - sum (sqrt (p q)))},
## the default for a variable holding levels.
## @item @qcode{'bhattacharyya'} @tab @tab @math{-log (sum (sqrt (p q)))}.
## @item @qcode{'tv'} @tab @tab The total variation distance,
## @math{sum (abs (p - q)) / 2}.
## @item @qcode{'psi'} @tab @tab The population stability index,
## @math{sum ((p - q) log (p / q))}.
## @item @qcode{'chi2'} @tab @tab @math{sum ((p - q)^2 / p)}.
## @end multitable
##
## For a variable holding levels, @math{p} and @math{q} are the shares of each
## level in @var{Baseline} and @var{Target}, over the levels either holds,
## each count increased by 0.5 first so that a level one sample lacks leaves
## every metric finite.
##
## The p-value of a variable is the share of permutations whose metric is at
## least the observed one, the observed arrangement counted as the first
## permutation, so it is never below one over their number.  Its 95%
## confidence interval is the Clopper-Pearson interval.  The drift status is
## @qcode{"Drift"} where the interval lies below @qcode{'DriftThreshold'},
## @qcode{"Stable"} where it lies above @qcode{'WarningThreshold'}, and
## @qcode{"Warning"} otherwise.  Permutations are drawn in stages, 1000 or
## @qcode{'MaxNumPermutations'} if fewer, then four times as many at each
## stage, until the interval lies within a single one of the three regions or
## the next stage would pass @qcode{'MaxNumPermutations'}.  They are drawn with
## @code{randperm}, so the state of @code{rand} decides them.
##
## The drift status of the data as a whole comes from the p-values corrected
## for testing every variable: by Bonferroni, the smallest p-value times the
## number of variables, or by the false discovery rate, the smallest
## Benjamini-Hochberg adjusted p-value.  It is @qcode{"Drift"} below
## @qcode{'DriftThreshold'}, @qcode{"Warning"} below
## @qcode{'WarningThreshold'}, and @qcode{"Stable"} otherwise.
##
## @code{@var{DDiagnostics} = detectdrift (@dots{}, @var{Name}, @var{Value})}
## takes the following options.
##
## @multitable @columnfractions 0.3 0.02 0.68
## @headitem @var{Name} @tab @tab @var{Value}
## @item @qcode{'VariableNames'} @tab @tab The variables to use, among those
## the two tables share.  All of them by default.
## @item @qcode{'CategoricalVariables'} @tab @tab The variables holding
## levels: @qcode{'all'}, their indices, a logical vector over the variables,
## or, for tables, their names.  A table variable holding logical values, an
## unordered @code{categorical} array, a @code{string} array or a cell array of
## character vectors, and every column of a @code{categorical} array, holds
## levels whether named here or not.
## @item @qcode{'ContinuousMetric'} @tab @tab @qcode{'wasserstein'}, the
## default, @qcode{'ks'}, @qcode{'ad'} or @qcode{'energy'}.
## @item @qcode{'CategoricalMetric'} @tab @tab @qcode{'hellinger'}, the
## default, @qcode{'bhattacharyya'}, @qcode{'tv'}, @qcode{'psi'} or
## @qcode{'chi2'}.
## @item @qcode{'DriftThreshold'} @tab @tab 0.05 by default.
## @item @qcode{'WarningThreshold'} @tab @tab 0.1 by default, and above
## @qcode{'DriftThreshold'}.
## @item @qcode{'MultipleTestCorrection'} @tab @tab @qcode{'bonferroni'}, the
## default, or @qcode{'fdr'}.
## @item @qcode{'MaxNumPermutations'} @tab @tab The most permutations drawn
## for a variable, 1000 by default.
## @item @qcode{'EstimatePValues'} @tab @tab @code{true} by default; with
## @code{false} only the metrics are computed.
## @item @qcode{'Options'} @tab @tab A structure as @code{statset} returns.
## Parallel computing and random streams are not implemented, and are refused
## where asked for.
## @end multitable
##
## MATLAB draws at least 1000 permutations whatever
## @qcode{'MaxNumPermutations'} says, and finishes a stage that passes it, so
## that a maximum of 3000 can end at 4001; here no variable is given more
## permutations than @qcode{'MaxNumPermutations'}.  MATLAB's stages are
## otherwise its own, so the number of permutations and the p-values differ
## from MATLAB's, as they would between any two random draws.
##
## @seealso{stats.drift.DriftDiagnostics, knntest, mmdtest}
## @end deftypefn

function DDiagnostics = detectdrift (Baseline, Target, varargin)

  ## Input validation
  if (nargin < 2)
    error ("detectdrift: too few input arguments.");
  endif
  optNames = {'VariableNames', 'CategoricalVariables', 'ContinuousMetric', ...
              'CategoricalMetric', 'DriftThreshold', 'WarningThreshold', ...
              'MultipleTestCorrection', 'MaxNumPermutations', ...
              'EstimatePValues', 'Options'};
  dfValues = {[], [], 'wasserstein', 'hellinger', 0.05, 0.1, 'bonferroni', ...
              1000, true, []};
  [vnames, catvars, cmetric, lmetric, dthr, wthr, mtc, maxperm, estp, ...
   opts, args] = parsePairedArguments (optNames, dfValues, varargin(:));
  if (! isempty (args))
    error ("detectdrift: invalid optional paired argument.");
  endif
  cmetric = ddWord (cmetric, {'wasserstein', 'ks', 'ad', 'energy'}, ...
                    'ContinuousMetric');
  lmetric = ddWord (lmetric, {'hellinger', 'bhattacharyya', 'tv', 'psi', ...
                              'chi2'}, 'CategoricalMetric');
  mtc = ddWord (mtc, {'bonferroni', 'fdr'}, 'MultipleTestCorrection');
  if (! (isnumeric (dthr) && isreal (dthr) && isscalar (dthr)
         && dthr > 0 && dthr < 1))
    error ("detectdrift: 'DriftThreshold' must be a scalar between 0 and 1.");
  endif
  if (! (isnumeric (wthr) && isreal (wthr) && isscalar (wthr)
         && wthr > 0 && wthr < 1))
    error (strcat ("detectdrift: 'WarningThreshold' must be a scalar", ...
                   " between 0 and 1."));
  endif
  if (! (dthr < wthr))
    error (strcat ("detectdrift: 'DriftThreshold' must be smaller than", ...
                   " 'WarningThreshold'."));
  endif
  if (! (isnumeric (maxperm) && isreal (maxperm) && isscalar (maxperm)
         && isfinite (maxperm) && maxperm >= 1 && maxperm == fix (maxperm)))
    error ("detectdrift: 'MaxNumPermutations' must be a positive integer.");
  endif
  if (! ((islogical (estp) || isnumeric (estp)) && isscalar (estp)))
    error ("detectdrift: 'EstimatePValues' must be a logical scalar.");
  endif
  if (! isempty (opts))
    if (! isstruct (opts))
      error ("detectdrift: 'Options' must be a structure.");
    endif
    if (isfield (opts, 'UseParallel') && ! isempty (opts.UseParallel)
        && ! strcmpi (opts.UseParallel, 'never')
        && ! isequal (opts.UseParallel, false))
      error ("detectdrift: parallel computing is not implemented.");
    endif
    if (isfield (opts, 'Streams') && ! isempty (opts.Streams))
      error ("detectdrift: random streams are not implemented.");
    endif
  endif

  ## A categorical array is read as a table of its columns, each of levels
  X = Baseline;
  Y = Target;
  if (isa (X, 'categorical') != isa (Y, 'categorical'))
    error (strcat ("detectdrift: Baseline and Target must both be", ...
                   " categorical arrays, or neither."));
  endif
  if (isa (X, 'categorical'))
    X = ddColumns (X);
    Y = ddColumns (Y);
  endif

  ## The pooled sample coded one column per variable, rows holding a missing
  ## value left out, and which variables hold levels
  [Z, iscat, mx, my, errmsg, names, labels] = ...
                __twosample__ (X, Y, vnames, catvars, {'Baseline', 'Target'});
  if (! isempty (errmsg))
    error ("detectdrift: %s", errmsg);
  endif
  K = columns (Z);
  m = mx + my;
  if (! istable (Baseline))
    names = arrayfun (@(j) sprintf ('x%d', j), 1:K, 'UniformOutput', false);
  endif

  ## The metric of each variable, observed and over the permutations
  metrics = cell (1, K);
  for j = 1:K
    if (iscat(j))
      metrics{j} = lmetric;
    else
      metrics{j} = cmetric;
    endif
  endfor
  obs = zeros (1, K);
  for j = 1:K
    obs(j) = ddMetric (Z(:,j), 1:mx, mx+1:m, metrics{j});
  endfor

  P = NaN (1, K);
  CI = NaN (2, K);
  NP = ones (1, K);
  perm = cell (K, 1);
  status = repmat ({''}, 1, K);
  overall = '';
  for j = 1:K
    perm{j} = obs(j);
  endfor
  if (estp)
    for j = 1:K
      [perm{j}, P(j), CI(:,j), status{j}] = ddPermute (Z(:,j), mx, ...
                                  metrics{j}, obs(j), maxperm, dthr, wthr);
      NP(j) = numel (perm{j});
    endfor
    overall = ddOverall (P, mtc, dthr, wthr);
  endif

  S = struct ();
  S.Baseline = Baseline;
  S.Target = Target;
  S.VariableNames = names;
  S.CategoricalVariables = find (iscat);
  S.Metrics = cellfun (@ddLabel, metrics, 'UniformOutput', false);
  S.MetricValues = obs;
  S.PValues = P;
  S.ConfidenceIntervals = CI;
  S.NumPermutations = NP;
  S.PermutationResults = perm;
  S.DriftStatus = status;
  S.MultipleTestCorrection = mtc;
  S.MultipleTestDriftStatus = overall;
  S.DriftThreshold = dthr;
  S.WarningThreshold = wthr;
  S.Data = Z;
  S.NumBaseline = mx;
  S.Labels = labels;
  S.Estimated = logical (estp);
  DDiagnostics = stats.drift.DriftDiagnostics (S);

endfunction

## An option holding one word of a list, matched in any case.
function w = ddWord (w, list, name)
  if (isa (w, 'string') && isscalar (w))
    w = char (w);
  endif
  if (! (ischar (w) && isrow (w) && any (strcmpi (w, list))))
    error ("detectdrift: '%s' must be one of %s.", name, ...
           strjoin (strcat ("'", list, "'"), ', '));
  endif
  w = lower (w);
endfunction

## The name a metric is reported by.
function s = ddLabel (name)
  switch (name)
    case 'wasserstein',   s = 'Wasserstein';
    case 'ks',            s = 'KolmogorovSmirnov';
    case 'ad',            s = 'AndersonDarling';
    case 'energy',        s = 'Energy';
    case 'hellinger',     s = 'Hellinger';
    case 'bhattacharyya', s = 'Bhattacharyya';
    case 'tv',            s = 'TotalVariation';
    case 'psi',           s = 'PopulationStabilityIndex';
    case 'chi2',          s = 'ChiSquare';
  endswitch
endfunction

## A categorical array as a table of its columns.
function T = ddColumns (C)
  K = columns (C);
  cols = cell (1, K);
  for j = 1:K
    cols{j} = C(:,j);
  endfor
  T = table (cols{:});
endfunction

## The metric of one variable with the rows A as the baseline and B as the
## target.
function v = ddMetric (z, a, b, name)

  switch (name)
    case {'wasserstein', 'ks', 'ad'}
      x = sort (z(a));
      y = sort (z(b));
      u = unique ([x; y]);
      F = ddCount (x, u) / numel (x);
      G = ddCount (y, u) / numel (y);
      switch (name)
        case 'wasserstein'
          v = sum (abs (F(1:end-1) - G(1:end-1)) .* diff (u));
        case 'ks'
          v = max (abs (F - G));
        case 'ad'
          H = (ddCount (x, u) + ddCount (y, u)) / (numel (x) + numel (y));
          k = H < 1;
          v = sum ((F(k) - G(k)) .^ 2 ./ (H(k) .* (1 - H(k)))) ...
              / (numel (x) + numel (y));
      endswitch
    case 'energy'
      x = sort (z(a));
      y = sort (z(b));
      exy = ddMeanAbs (x, y);
      exx = ddMeanAbs (x, x);
      eyy = ddMeanAbs (y, y);
      v = sqrt (max (0, 2 * exy - exx - eyy));
    otherwise
      ## Levels 1..L, each count increased by 0.5
      L = max (z);
      cx = accumarray (z(a), 1, [L, 1]) + 0.5;
      cy = accumarray (z(b), 1, [L, 1]) + 0.5;
      p = cx / sum (cx);
      q = cy / sum (cy);
      switch (name)
        case 'hellinger'
          v = sqrt (max (0, 1 - sum (sqrt (p .* q))));
        case 'bhattacharyya'
          v = -log (sum (sqrt (p .* q)));
        case 'tv'
          v = sum (abs (p - q)) / 2;
        case 'psi'
          v = sum ((p - q) .* log (p ./ q));
        case 'chi2'
          v = sum ((p - q) .^ 2 ./ p);
      endswitch
  endswitch

endfunction

## How many of the sorted values S are at most each value of U.
function c = ddCount (s, u)
  c = lookup (s, u);
  c = c(:);
endfunction

## The mean of |x_i - y_j| over every pair, from the sorted X and Y.
function e = ddMeanAbs (x, y)
  ## Each x_i against the y at most it and the y above it
  cy = [0; cumsum(y)];
  k = lookup (y, x);
  k = k(:);
  below = x .* k - cy(k + 1);
  above = (cy(end) - cy(k + 1)) - x .* (numel (y) - k);
  e = sum (below + above) / (numel (x) * numel (y));
endfunction

## The permutations of one variable, in stages, and what they say.
function [vals, p, ci, status] = ddPermute (z, mx, name, obs, maxperm, ...
                                           dthr, wthr)

  m = numel (z);
  vals = obs;
  count = 1;
  tol = 1e-12 * abs (obs);
  target = min (1000, maxperm);
  while (true)
    add = target - numel (vals);
    new = zeros (add, 1);
    for i = 1:add
      r = randperm (m);
      new(i) = ddMetric (z, r(1:mx), r(mx+1:end), name);
    endfor
    vals = [vals; new];
    count += sum (new >= obs - tol);
    N = numel (vals);
    p = count / N;
    ci = ddClopperPearson (count, N);
    done = ci(2) < dthr || ci(1) > wthr || (ci(1) > dthr && ci(2) < wthr);
    if (done || N >= maxperm)
      break;
    endif
    target = min (4 * N, maxperm);
  endwhile
  if (ci(2) < dthr)
    status = 'Drift';
  elseif (ci(1) > wthr)
    status = 'Stable';
  else
    status = 'Warning';
  endif

endfunction

## The 95% Clopper-Pearson interval of X successes in N trials.
function ci = ddClopperPearson (x, N)
  ci = [0; 1];
  if (x > 0)
    ci(1) = betaincinv (0.025, x, N - x + 1);
  endif
  if (x < N)
    ci(2) = betaincinv (0.975, x + 1, N - x);
  endif
endfunction

## The drift status of the data as a whole.
function s = ddOverall (P, mtc, dthr, wthr)
  k = numel (P);
  if (strcmp (mtc, 'bonferroni'))
    pc = min (P) * k;
  else
    ps = sort (P);
    pc = min (ps .* k ./ (1:k));
  endif
  if (pc < dthr)
    s = 'Drift';
  elseif (pc < wthr)
    s = 'Warning';
  else
    s = 'Stable';
  endif
endfunction

## Expected values from MATLAB R2024a
%!shared X, Y, C, T
%! X = [0; 1; 2; 3; 7];
%! Y = [1; 2; 5; 6];
%! C = categorical ([1; 1; 2; 3; 3; 3]);
%! T = categorical ([1; 2; 2; 2; 3]);
%!test
%! D = detectdrift (X, Y, 'EstimatePValues', false);
%! assert_equal (D.MetricValues, 1.3, -1e-14);
%! assert_equal (D.Metrics, string ('Wasserstein'));
%! assert_equal (D.VariableNames, string ('x1'));
%!assert_equal (detectdrift (X, Y, 'EstimatePValues', false, ...
%!                          'ContinuousMetric', 'ks').MetricValues, 0.3, -1e-14)
%!assert_equal (detectdrift (X, Y, 'EstimatePValues', false, ...
%!                          'ContinuousMetric', 'AD').MetricValues, ...
%!              0.152357142857143, -1e-13)
%!assert_equal (detectdrift (X, Y, 'EstimatePValues', false, ...
%!                          'ContinuousMetric', 'energy').MetricValues, ...
%!              0.768114574786861, -1e-13)
%!assert_equal (detectdrift (C, T, 'EstimatePValues', false).MetricValues, ...
%!              0.257526267734947, -1e-13)
%!assert_equal (detectdrift (C, T, 'EstimatePValues', false, ...
%!                          'CategoricalMetric', 'tv').MetricValues, ...
%!              0.338461538461538, -1e-13)
%!assert_equal (detectdrift (C, T, 'EstimatePValues', false, ...
%!                          'CategoricalMetric', 'psi').MetricValues, ...
%!              0.539045501736854, -1e-13)
%!assert_equal (detectdrift (C, T, 'EstimatePValues', false, ...
%!                          'CategoricalMetric', 'chi2').MetricValues, ...
%!              0.723584108199493, -1e-13)
%!assert_equal (detectdrift (C, T, 'EstimatePValues', false, ...
%!                          'CategoricalMetric', 'bhattacharyya').MetricValues, ...
%!              0.068621274723466, -1e-12)
%!test
%! ## A level the target lacks leaves every metric finite
%! T2 = categorical ([1; 2; 2; 2]);
%! D = detectdrift (C, T2, 'EstimatePValues', false, 'CategoricalMetric', 'psi');
%! assert_equal (D.MetricValues, 1.131879584405866, -1e-13);
%! D = detectdrift (C, T2, 'EstimatePValues', false, 'CategoricalMetric', 'chi2');
%! assert_equal (D.MetricValues, 1.265643447461629, -1e-13);
%! D = detectdrift (C, T2, 'EstimatePValues', false);
%! assert_equal (D.MetricValues, 0.368461885678999, -1e-13);
%!test
%! D = detectdrift ([0, 5; 1, 3; 2, 8; 3, 1; 7, 4], [1, 9; 2, 7; 5, 8; 6, 6], ...
%!                  'EstimatePValues', false);
%! assert_equal (D.MetricValues, [1.3, 3.3], -1e-14);
%! assert_equal (D.VariableNames, string ({'x1', 'x2'}));
%!test
%! TB = table ([0; 1; 2; 3; 7; 4], categorical ([1; 1; 2; 3; 3; 3]), ...
%!             'VariableNames', {'A', 'B'});
%! TT = table ([1; 2; 5; 6; 3], categorical ([1; 2; 2; 2; 3]), ...
%!             'VariableNames', {'A', 'B'});
%! D = detectdrift (TB, TT, 'EstimatePValues', false);
%! assert_equal (D.MetricValues, [0.9, 0.257526267734947], -1e-13);
%! assert_equal (D.Metrics, string ({'Wasserstein', 'Hellinger'}));
%! assert_equal (D.CategoricalVariables, 2);
%!test
%! D = detectdrift (X, Y, 'EstimatePValues', false);
%! assert_equal (ismissing (D.DriftStatus), true);
%! assert_equal (ismissing (D.MultipleTestDriftStatus), true);
%! assert_equal (D.PValues, NaN);
%! assert_equal (D.ConfidenceIntervals, [NaN; NaN]);
%! assert_equal (D.NumPermutations, 1);
%! assert_equal (D.PermutationResults.PermutationResults, {1.3});
%!test
%! ## A shift no permutation reaches: the observed arrangement alone counts
%! rand ('seed', 1);
%! x = (1:60)' / 10;
%! D = detectdrift ([x, x], [x + 10, x([2:60, 1])]);
%! assert_equal (D.PValues, [0.001, 1]);
%! assert_equal (D.ConfidenceIntervals, [2.53174874912919e-05, ...
%!               0.996317916103134; 0.00555892427982512, 1], -1e-10);
%! assert_equal (D.DriftStatus, string ({'Drift', 'Stable'}));
%! assert_equal (D.NumPermutations, [1000, 1000]);
%! assert_equal (D.PermutationResults.PermutationResults{1}(1), 10, -1e-14);
%!test
%! ## MaxNumPermutations is never passed
%! rand ('seed', 1);
%! x = (1:60)' / 10;
%! D = detectdrift (x, x + 10, 'MaxNumPermutations', 50);
%! assert_equal (D.NumPermutations, 50);
%! assert_equal (D.PValues, 0.02);
%! assert_equal (D.DriftStatus, string ('Warning'));
%!test
%! ## Bonferroni takes 4 * 0.001 and the false discovery rate 0.001
%! rand ('seed', 1);
%! x = (1:60)' / 10;
%! Z = repmat (x, 1, 4);
%! B = detectdrift (Z, Z + 10, 'DriftThreshold', 0.002, 'WarningThreshold', 0.005);
%! assert_equal (B.MultipleTestDriftStatus, string ('Warning'));
%! assert_equal (B.MultipleTestCorrection, string ('Bonferroni'));
%! F = detectdrift (Z, Z + 10, 'DriftThreshold', 0.002, ...
%!                  'WarningThreshold', 0.005, 'MultipleTestCorrection', 'fdr');
%! assert_equal (F.MultipleTestDriftStatus, string ('Drift'));
%! assert_equal (F.MultipleTestCorrection, string ('FalseDiscoveryRate'));

%!error<detectdrift: too few input arguments.> detectdrift (1)
%!error<detectdrift: invalid optional paired argument.> detectdrift (X, Y, 'Tail', 1)
%!error<detectdrift: 'ContinuousMetric' must be one of 'wasserstein', 'ks', 'ad', 'energy'.> ...
%! detectdrift (X, Y, 'ContinuousMetric', 'kl')
%!error<detectdrift: 'CategoricalMetric' must be one of 'hellinger', 'bhattacharyya', 'tv', 'psi', 'chi2'.> ...
%! detectdrift (X, Y, 'CategoricalMetric', 'kl')
%!error<detectdrift: 'MultipleTestCorrection' must be one of 'bonferroni', 'fdr'.> ...
%! detectdrift (X, Y, 'MultipleTestCorrection', 'holm')
%!error<detectdrift: 'DriftThreshold' must be a scalar between 0 and 1.> ...
%! detectdrift (X, Y, 'DriftThreshold', 1)
%!error<detectdrift: 'WarningThreshold' must be a scalar between 0 and 1.> ...
%! detectdrift (X, Y, 'WarningThreshold', 0)
%!error<detectdrift: 'DriftThreshold' must be smaller than 'WarningThreshold'.> ...
%! detectdrift (X, Y, 'DriftThreshold', 0.1, 'WarningThreshold', 0.1)
%!error<detectdrift: 'MaxNumPermutations' must be a positive integer.> ...
%! detectdrift (X, Y, 'MaxNumPermutations', 0)
%!error<detectdrift: 'EstimatePValues' must be a logical scalar.> ...
%! detectdrift (X, Y, 'EstimatePValues', 'no')
%!error<detectdrift: 'Options' must be a structure.> ...
%! detectdrift (X, Y, 'Options', 1)
%!error<detectdrift: parallel computing is not implemented.> ...
%! detectdrift (X, Y, 'Options', struct ('UseParallel', true))
%!error<detectdrift: random streams are not implemented.> ...
%! detectdrift (X, Y, 'Options', struct ('Streams', 1))
%!error<detectdrift: Baseline and Target must both be categorical arrays, or neither.> ...
%! detectdrift (C, X)
%!error<detectdrift: Baseline and Target must have the same number of columns.> ...
%! detectdrift ([X, X], Y)
%!error<detectdrift: Baseline and Target must each hold an observation with no missing value.> ...
%! detectdrift (NaN, Y)
