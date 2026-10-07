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
## @deftypefn  {statistics} {@var{Effect} =} meanEffectSize (@var{X})
## @deftypefnx {statistics} {@var{Effect} =} meanEffectSize (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{Effect} =} meanEffectSize (@dots{}, @var{Name}, @var{Value})
##
## Effect sizes for the difference between two means, with their confidence
## intervals.
##
## @code{@var{Effect} = meanEffectSize (@var{X})} returns the difference
## between the mean of the sample @var{X} and a known population mean, zero by
## default, with its confidence interval.
##
## @code{@var{Effect} = meanEffectSize (@var{X}, @var{Y})} returns the
## difference between the means of the samples @var{X} and @var{Y}, with its
## confidence interval.
##
## @var{X} and @var{Y} are vectors of type double or single.  A missing value
## (@qcode{NaN}) is left out of its sample, and out of both where the samples
## are paired.  Where either sample is single, so are the results.
##
## @var{Effect} is a table with one row for each effect size asked for, in the
## order asked, named as listed below.  Its variable @qcode{Effect} holds the
## effect size and its variable @qcode{ConfidenceIntervals} the lower and upper
## bounds of its confidence interval, the latter absent where
## @qcode{'ConfidenceIntervalType'} is @qcode{'none'}.
##
## @code{@var{Effect} = meanEffectSize (@dots{}, @var{Name}, @var{Value})}
## takes the following options.
##
## @multitable @columnfractions 0.28 0.02 0.70
## @headitem @var{Name} @tab @tab @var{Value}
## @item @qcode{'Effect'} @tab @tab The effect sizes to compute, as a
## character vector, a string array or a cell array of character vectors,
## chosen among those below.  A name may be shortened to any prefix that
## names one effect alone, and case does not matter.
## @qcode{'meandiff'} by default.
## @item @qcode{'Mean'} @tab @tab The known population mean that one sample
## is compared with, a real scalar, 0 by default.  It is ignored where there
## are two samples.
## @item @qcode{'Paired'} @tab @tab Whether @var{X} and @var{Y} are paired
## observations, @qcode{true} or @qcode{false} (also @qcode{'on'} or
## @qcode{'off'}), @qcode{false} by default.  Paired samples must have the
## same number of elements.
## @item @qcode{'VarianceType'} @tab @tab Whether two unpaired samples are
## taken to come from populations of @qcode{'equal'} or @qcode{'unequal'}
## variance, @qcode{'equal'} by default.  Paired samples take
## @qcode{'equal'} only; one sample ignores it.
## @item @qcode{'Alpha'} @tab @tab The significance level, a scalar between 0
## and 1, 0.05 by default, so that each interval covers
## @math{100(1 - @var{Alpha})} percent.
## @item @qcode{'ConfidenceIntervalType'} @tab @tab @qcode{'exact'},
## @qcode{'bootstrap'} or @qcode{'none'}.  By default each effect takes
## @qcode{'exact'} where it has an exact interval and @qcode{'bootstrap'}
## where it has not.
## @item @qcode{'NumBootstraps'} @tab @tab The number of bootstrap replicates,
## a positive integer, 1000 by default.
## @item @qcode{'BootstrapOptions'} @tab @tab A structure, as returned by
## @code{statset}, accepted for compatibility.  The replicates are always
## computed serially and drawn from the generator of @code{rand}.
## @item @qcode{'Resampling'} @tab @tab How the bootstrap resamples two
## unpaired samples, @qcode{'pooled'} or @qcode{'stratified'}; see below.
## @qcode{'pooled'} by default.  This option is an Octave extension.
## @end multitable
##
## The effect sizes are the following, where @math{J(v)} is Hedges'
## correction for bias, @math{J(v) = Gamma(v/2) / (sqrt(v/2) Gamma((v-1)/2))}.
##
## @multitable @columnfractions 0.14 0.35 0.51
## @headitem @var{Effect} @tab Row name @tab Description
## @item @qcode{'meandiff'} @tab @qcode{MeanDifference} @tab The difference
## between the means, @math{mean(X) - mean(Y)}, or @math{mean(X) - mu} for
## one sample.  Its exact interval is Student's t interval: pooled for equal
## variances, Welch's for unequal variances, and on the differences for
## paired samples.
## @item @qcode{'cohen'} @tab @qcode{CohensD} @tab Cohen's d corrected for
## bias, that is Hedges' g: the difference between the means divided by a
## standard deviation, times @math{J(v)}.  See below.
## @item @qcode{'glass'} @tab @qcode{GlasssDelta} @tab Glass's delta,
## @math{J(n_x - 1) (mean(X) - mean(Y)) / s_x}, the difference between the
## means divided by the standard deviation of the control sample @var{X}.
## Two unpaired samples only; @qcode{'VarianceType'} does not change it.
## @item @qcode{'cliff'} @tab @qcode{CliffsDelta} @tab Cliff's delta, the
## share of pairs @math{(x_i, y_j)} with @math{x_i > y_j} less the share with
## @math{x_i < y_j}.  For paired samples the pairs are those of different
## observations, @math{i != j}.  Two samples only.
## @item @qcode{'mediandiff'} @tab @qcode{MedianDifference} @tab The
## difference between the medians, @math{median(X) - median(Y)}, paired
## samples included.  Two samples only.
## @item @qcode{'robustcohen'} @tab @qcode{RobustCohensD} @tab A robust
## Cohen's d, as MATLAB computes it; see below.
## @item @qcode{'akpcohen'} @tab @qcode{AKPCohensD} @tab The robust Cohen's d
## of Algina, Keselman and Penfield; see below.  This effect is an Octave
## extension.
## @item @qcode{'kstest'} @tab @qcode{KolmogorovSmirnovStatistic} @tab The
## two-sample Kolmogorov-Smirnov statistic, the largest distance between the
## empirical distribution functions of @var{X} and @var{Y}, paired samples
## included.  Two samples only.
## @end multitable
##
## @qcode{'meandiff'}, @qcode{'cohen'}, @qcode{'glass'} and @qcode{'cliff'}
## have exact intervals; the others are bootstrap only.
##
## @strong{Cohen's d.}  For one sample it is
## @math{J(n - 1) (mean(X) - mu) / s_x}, and for two unpaired samples of
## equal variance @math{J(n_x + n_y - 2) (mean(X) - mean(Y)) / s_p}, with
## @math{s_p} the pooled standard deviation.  Its interval inverts the
## noncentral t distribution of the one-sample or the two-sample t statistic.
## For unequal variances it is the @math{g^*} of Delacre and colleagues: the
## difference between the means divided by
## @math{s_a = sqrt((s_x^2 + s_y^2) / 2)}, times @math{J(v)} with
## @math{v = (n_x - 1)(n_y - 1)(s_x^2 + s_y^2)^2 / ((n_y - 1) s_x^4 +
## (n_x - 1) s_y^4)}, and its interval inverts the noncentral t distribution
## of Welch's statistic on @math{v} degrees of freedom.  For paired samples it
## is @math{J(n - 1) mean(X - Y) / s_a}, and its interval is the MAG interval
## of Cousineau and Goulet-Pelletier.  The interval of Glass's delta inverts
## the noncentral t distribution of Welch's statistic on @math{n_x - 1}
## degrees of freedom.
##
## @strong{Cliff's delta.}  Its interval is Cliff's asymmetric interval from
## the normal distribution, with the variance of @math{delta} estimated as
## @math{((n_y - 1) var(d_i) + (n_x - 1) var(d_j) + S / ((n_x - 1)(n_y - 1)))
## / (n_x n_y)}, where @math{d_i} and @math{d_j} are the row and column means
## of the dominance matrix @math{d_ij = sign(x_i - y_j)} and @math{S} is the
## sum of @math{(d_ij - delta)^2}.  For paired samples the interval is
## @math{delta +/- z s} with the variance of the U statistic over the
## @math{n(n - 1)} pairs, which needs four pairs or more.
##
## @strong{Robust Cohen's d.}  @qcode{'robustcohen'} is MATLAB's:
## @math{0.642 J(v) (T(X) - T(Y)) / s_w}, where @math{T} is
## @code{trimmean (@dots{}, 20)}, which trims 10 percent of each tail, and
## @math{s_w} the standard deviation of the values clipped at their 20th and
## 80th percentiles, pooled over the two samples as @math{s_p} is, or averaged
## as @math{s_a} is for unequal variances and paired samples.  @math{v} is the
## degrees of freedom of the matching @qcode{'cohen'}.  For one sample
## @math{T(Y)} is replaced by @math{mu}.
##
## @qcode{'akpcohen'} is the estimator of Algina, Keselman and Penfield, as
## the R package WRS2 computes it: @math{c (T(X) - T(Y)) / s_w}, where
## @math{T} trims 20 percent of each tail, @math{floor(0.2 n)} values, and
## @math{s_w} is the standard deviation of the sample winsorized at the same
## order statistics.  Two samples of equal variance pool @math{s_w} as
## @math{s_p} is pooled; for unequal variances @math{s_w} is that of @var{X}
## alone.  One sample takes @math{c (T(X) - mu) / s_w}, and paired samples
## the same over the differences @math{X - Y}, with @math{mu} zero.  There is
## no correction for bias.  The constant @math{c = 0.64194}, which the
## authors round to 0.642, is the standard deviation of a standard normal
## distribution winsorized at 20 percent, so that both estimators estimate
## Cohen's d where the data are normal; @qcode{'akpcohen'} trims more of each
## tail and is the more robust of the two.
##
## @strong{Bootstrap.}  The bootstrap interval is the bias-corrected and
## accelerated (BCa) percentile interval, with the acceleration estimated by
## the jackknife.  One sample is resampled with replacement, and paired
## samples by pairs.  Two unpaired samples are resampled as
## @qcode{'Resampling'} says.  @qcode{'pooled'} resamples the observations of
## both samples together, each keeping its sample, so that the size of each
## sample varies from one replicate to the next; this is what MATLAB does,
## and it suits data where the sample an observation falls in is itself
## random.  @qcode{'stratified'} resamples each sample on its own, keeping
## their sizes; it suits data where the sizes were fixed by design, as in an
## experiment.  A replicate in which a sample comes out empty is left out.
##
## An interval that cannot be computed, because a sample holds too few
## observations or has zero variance, is @qcode{NaN}, with a warning.
##
## MATLAB returns a paired Cliff's delta interval of zero width,
## @math{[delta, delta]}, for fewer than four pairs, where the variance
## cannot be estimated; here that interval is @qcode{NaN}.
##
## References: J. Algina, H. J. Keselman and R. D. Penfield (2005).  An
## alternative to Cohen's standardized mean difference effect size: a robust
## parameter and confidence interval in the two independent groups case.
## Psychological Methods, 10(3), 317-328.  D. Cousineau and
## J.-C. Goulet-Pelletier (2021).  A study of confidence intervals for Cohen's
## dp in within-subject designs with new proposals.  The Quantitative Methods
## for Psychology, 17(1), 51-75.  M. Delacre, D. Lakens, C. Ley, L. Liu and
## C. Leys (2021).  Why Hedges' g*s based on the non-pooled standard deviation
## should be reported with Welch's t-test.  PsyArXiv.  N. Cliff (1993).
## Dominance statistics: ordinal analyses to answer ordinal questions.
## Psychological Bulletin, 114(3), 494-509.
##
## @seealso{ttest, ttest2, ranksum, kstest2, trimmean}
## @end deftypefn

function Effect = meanEffectSize (X, varargin)

  ## Input validation
  if (nargin < 1)
    error ("meanEffectSize: too few input arguments.");
  endif
  if (! isSample (X))
    error (strcat ("meanEffectSize: X must be a nonempty vector of type", ...
                   " double or single."));
  endif
  Y = [];
  if (numel (varargin) > 0 && ! (ischar (varargin{1})
                                 || isa (varargin{1}, 'string')))
    Y = varargin{1};
    varargin(1) = [];
    if (! isSample (Y))
      error (strcat ("meanEffectSize: Y must be a nonempty vector of type", ...
                     " double or single."));
    endif
  endif
  optNames = {'Effect', 'Mean', 'Paired', 'VarianceType', 'Alpha', ...
              'ConfidenceIntervalType', 'NumBootstraps', ...
              'BootstrapOptions', 'Resampling'};
  ## An empty ConfidenceIntervalType resolves to each effect's default
  dfValues = {'meandiff', 0, false, 'equal', 0.05, [], 1000, [], 'pooled'};
  [effects, mu, paired, vartype, alpha, citype, nboot, bopts, resampling, ...
   args] = parsePairedArguments (optNames, dfValues, varargin(:));
  if (! isempty (args))
    error ("meanEffectSize: invalid optional paired argument.");
  endif
  names = {'meandiff', 'cohen', 'glass', 'cliff', 'mediandiff', ...
           'robustcohen', 'akpcohen', 'kstest'};
  effects = matchEffects (effects, names);
  if (isempty (effects))
    error ("meanEffectSize: 'Effect' must be one or more of %s.", ...
           strjoin (strcat ("'", names, "'"), ', '));
  endif
  if (! (isnumeric (mu) && isreal (mu) && isscalar (mu)))
    error ("meanEffectSize: 'Mean' must be a real scalar.");
  endif
  paired = matchSwitch (paired);
  if (isempty (paired))
    error ("meanEffectSize: 'Paired' must be true, false, 'on' or 'off'.");
  endif
  vartype = matchOption (vartype, {'equal', 'unequal'});
  if (isempty (vartype))
    error ("meanEffectSize: 'VarianceType' must be 'equal' or 'unequal'.");
  endif
  if (! (isnumeric (alpha) && isreal (alpha) && isscalar (alpha)
         && alpha > 0 && alpha < 1))
    error ("meanEffectSize: 'Alpha' must be a scalar between 0 and 1.");
  endif
  if (! isempty (citype))
    citype = matchOption (citype, {'exact', 'bootstrap', 'none'});
    if (isempty (citype))
      error (strcat ("meanEffectSize: 'ConfidenceIntervalType' must be", ...
                     " 'exact', 'bootstrap' or 'none'."));
    endif
  endif
  if (! (isnumeric (nboot) && isreal (nboot) && isscalar (nboot)
         && isfinite (nboot) && nboot >= 1 && nboot == fix (nboot)))
    error ("meanEffectSize: 'NumBootstraps' must be a positive integer.");
  endif
  resampling = matchOption (resampling, {'pooled', 'stratified'});
  if (isempty (resampling))
    error ("meanEffectSize: 'Resampling' must be 'pooled' or 'stratified'.");
  endif
  ## BootstrapOptions is checked and otherwise unused
  if (! (isempty (bopts) || isstruct (bopts)))
    error ("meanEffectSize: 'BootstrapOptions' must be a structure.");
  endif

  ## The design, and which effects it admits
  if (isempty (Y))
    mode = 'one';
    if (paired)
      error ("meanEffectSize: 'Paired' needs two samples.");
    endif
    twoOnly = {'glass', 'cliff', 'mediandiff', 'kstest'};
    bad = effects(ismember (effects, twoOnly));
    if (! isempty (bad))
      error ("meanEffectSize: effect '%s' needs two samples.", bad{1});
    endif
  elseif (paired)
    mode = 'paired';
    if (numel (X) != numel (Y))
      error (strcat ("meanEffectSize: paired samples X and Y must have", ...
                     " the same number of elements."));
    endif
    if (strcmp (vartype, 'unequal'))
      error (strcat ("meanEffectSize: paired samples take 'VarianceType'", ...
                     " 'equal' only."));
    endif
    if (any (strcmp (effects, 'glass')))
      error (strcat ("meanEffectSize: effect 'glass' does not apply to", ...
                     " paired samples."));
    endif
  else
    mode = 'two';
  endif
  isExact = {'meandiff', 'cohen', 'glass', 'cliff'};
  if (strcmp (citype, 'exact'))
    bad = effects(! ismember (effects, isExact));
    if (! isempty (bad))
      error (strcat ("meanEffectSize: effect '%s' has no exact confidence", ...
                     " interval; use 'bootstrap'."), bad{1});
    endif
  endif

  ## Missing values left out, of both samples where they are paired
  outclass = 'double';
  if (isa (X, 'single') || isa (Y, 'single'))
    outclass = 'single';
  endif
  x = double (X(:));
  y = double (Y(:));
  if (strcmp (mode, 'paired'))
    keep = ! (isnan (x) | isnan (y));
    x = x(keep);
    y = y(keep);
  else
    x = x(! isnan (x));
    y = y(! isnan (y));
  endif

  ## Each effect and its interval
  ne = numel (effects);
  E = zeros (ne, 1);
  CI = NaN (ne, 2);
  rowNames = cell (ne, 1);
  for ii = 1:ne
    name = effects{ii};
    rowNames{ii} = rowName (name);
    stat = @(a, b) effectValue (name, a, b, mode, vartype, mu);
    E(ii) = stat (x, y);
    ci = citype;
    if (isempty (ci))
      if (ismember (name, isExact))
        ci = 'exact';
      else
        ci = 'bootstrap';
      endif
    endif
    switch (ci)
      case 'exact'
        [CI(ii,:), why] = exactInterval (name, x, y, mode, vartype, mu, alpha);
        if (! isempty (why))
          warning (strcat ("meanEffectSize: the confidence interval of", ...
                           " effect '%s' is NaN, since %s."), name, why);
        endif
      case 'bootstrap'
        switch (mode)
          case 'one'
            CI(ii,:) = __bootci__ (@(a) stat (a, y), {x}, alpha, nboot, 'rows');
          case 'paired'
            CI(ii,:) = __bootci__ (stat, {x, y}, alpha, nboot, 'rows');
          otherwise
            if (strcmp (resampling, 'pooled'))
              CI(ii,:) = __bootci__ (stat, {x, y}, alpha, nboot, 'pooled');
            else
              CI(ii,:) = __bootci__ (stat, {x, y}, alpha, nboot, 'strata');
            endif
        endswitch
    endswitch
  endfor

  E = cast (E, outclass);
  if (strcmp (citype, 'none'))
    Effect = table (E, 'VariableNames', {'Effect'}, 'RowNames', rowNames);
  else
    CI = cast (CI, outclass);
    Effect = table (E, CI, 'VariableNames', {'Effect', ...
                    'ConfidenceIntervals'}, 'RowNames', rowNames);
  endif

endfunction

## A nonempty vector of type double or single
function tf = isSample (v)
  tf = isfloat (v) && isreal (v) && isvector (v) && ! isempty (v);
endfunction

## The effects asked for, in the order asked, each once; empty if any name
## is not recognized
function out = matchEffects (in, names)
  out = {};
  if (isa (in, 'string'))
    in = cellstr (in);
  elseif (ischar (in) && (isrow (in) || isempty (in)))
    in = {in};
  endif
  if (! iscellstr (in) || isempty (in))
    return;
  endif
  for ii = 1:numel (in)
    name = matchOption (in{ii}, names);
    if (isempty (name))
      out = {};
      return;
    endif
    if (! any (strcmp (name, out)))
      out{end+1} = name;
    endif
  endfor
endfunction

## The option in LIST that VAL names, case free and possibly shortened;
## empty if none or several match
function out = matchOption (val, list)
  out = '';
  if (isa (val, 'string') && isscalar (val))
    val = char (val);
  endif
  if (! (ischar (val) && isrow (val)))
    return;
  endif
  idx = strcmpi (val, list);
  if (! any (idx))
    idx = strncmpi (val, list, numel (val));
  endif
  if (sum (idx) == 1)
    out = list{idx};
  endif
endfunction

## A logical switch given as true, false, 0, 1, 'on' or 'off'; empty if not
function out = matchSwitch (val)
  out = [];
  if ((islogical (val) || isnumeric (val)) && isscalar (val)
      && (val == 0 || val == 1))
    out = logical (val);
  else
    val = matchOption (val, {'on', 'off'});
    if (! isempty (val))
      out = strcmp (val, 'on');
    endif
  endif
endfunction

## The row name of each effect
function out = rowName (name)
  switch (name)
    case 'meandiff'
      out = 'MeanDifference';
    case 'cohen'
      out = 'CohensD';
    case 'glass'
      out = 'GlasssDelta';
    case 'cliff'
      out = 'CliffsDelta';
    case 'mediandiff'
      out = 'MedianDifference';
    case 'robustcohen'
      out = 'RobustCohensD';
    case 'akpcohen'
      out = 'AKPCohensD';
    case 'kstest'
      out = 'KolmogorovSmirnovStatistic';
  endswitch
endfunction

## Hedges' correction for bias on V degrees of freedom, 0 at its limit V = 1
## and undefined below it; elementwise
function J = hedgesJ (v)
  J = NaN (size (v));
  k = v >= 1;
  J(k) = exp (gammaln (v(k) / 2) - log (sqrt (v(k) / 2)) ...
              - gammaln ((v(k) - 1) / 2));
endfunction

## Delacre's degrees of freedom for two samples of sizes N1, N2 and
## variances V1, V2; elementwise
function v = delacreDF (n1, n2, v1, v2)
  v = (n1 - 1) * (n2 - 1) * (v1 + v2) .^ 2 ./ ((n2 - 1) * v1 .^ 2 ...
                                               + (n1 - 1) * v2 .^ 2);
endfunction

## The trimmed means and winsorized variances of the columns of V, MATLAB's
## way (trimmean (V, 20), which trims 10 percent of each tail, and the
## values clipped at the 20th and 80th percentiles) or that of Algina,
## Keselman and Penfield (20 percent of each tail by order statistics)
function [tm, wv] = robustParts (v, akp)
  n = rows (v);
  v = sort (v, 1);
  if (akp)
    g = floor (0.2 * n);
    tm = mean (v(g+1:n-g,:), 1);
    v(1:g,:) = repmat (v(g+1,:), g, 1);
    v(n-g+1:n,:) = repmat (v(n-g,:), g, 1);
  else
    ## As trimmean (V, 20) rounds, a half downwards
    k = round (n * 0.1 - eps (n * 0.1));
    tm = mean (v(k+1:n-k,:), 1);
    b = quantile (v, [0.2; 0.8], 1, 5);
    v = min (max (v, b(1,:)), b(2,:));
  endif
  wv = var (v, 0, 1);
endfunction

## Cliff's delta of two unpaired samples, with the row and column means of
## the dominance matrix and the sum of the squares of its deviations from
## delta, all by counting
function [d, di, dj, S] = cliffParts (x, y)
  n1 = numel (x);
  n2 = numel (y);
  ## For each x_i the number of y at or below it and at or above it
  yle = lookup (sort (y), x);
  yge = lookup (sort (-y), -x);
  d = sum (yle - yge) / (n1 * n2);
  if (nargout > 1)
    di = (yle - yge) / n2;
    dj = (lookup (sort (-x), -y) - lookup (sort (x), y)) / n1;
    ## Squares of d_ij are 1 unless tied, so their sum counts untied pairs
    S = n1 * n2 - sum (yle + yge - n2) - n1 * n2 * d ^ 2;
  endif
endfunction

## Cliff's delta over the pairs of different observations of paired samples,
## with the row and column means of the dominance matrix off its diagonal,
## the sum of the squares of its deviations from delta, S2, and the sum of
## the products of each deviation with its transpose, S3
function [d, di, dj, S2, S3] = cliffPaired (x, y)
  n = numel (x);
  m = n * (n - 1);
  s = sign (x - y);
  yle = lookup (sort (y), x);
  yge = lookup (sort (-y), -x);
  d = (sum (yle - yge) - sum (s)) / m;
  if (nargout > 1)
    di = (yle - yge - s) / (n - 1);
    dj = (lookup (sort (-x), -y) - lookup (sort (x), y) - s) / (n - 1);
    ## Untied pairs off the diagonal
    S2 = n * n - sum (yle + yge - n) - sum (s != 0) - m * d ^ 2;
    ## The sum of d_ij d_ji off the diagonal, over blocks of rows
    P = 0;
    step = max (1, floor (1e6 / n));
    for i1 = 1:step:n
      i = i1:min (i1 + step - 1, n);
      P += sum (sum (sign (x(i) - y') .* sign (x' - y(i))));
    endfor
    P -= sum (s .^ 2);
    S3 = P - 2 * d * d * m + m * d ^ 2;
  endif
endfunction

## The two-sample Kolmogorov-Smirnov statistic
function D = ksStat (x, y)
  z = [x; y];
  D = max (abs (lookup (sort (x), z) / numel (x) ...
                - lookup (sort (y), z) / numel (y)));
endfunction

## The value of one effect over each column of X, and of Y where there are
## two samples; MODE is 'one', 'paired' or 'two'
function e = effectValue (name, x, y, mode, vartype, mu)
  n1 = rows (x);
  n2 = rows (y);
  switch (name)
    case 'meandiff'
      switch (mode)
        case 'one'
          e = mean (x, 1) - mu;
        case 'paired'
          e = mean (x - y, 1);
        otherwise
          e = mean (x, 1) - mean (y, 1);
      endswitch
    case 'cohen'
      switch (mode)
        case 'one'
          e = hedgesJ (n1 - 1) * (mean (x, 1) - mu) ./ std (x, 0, 1);
        case 'paired'
          e = hedgesJ (n1 - 1) * mean (x - y, 1) ...
              ./ sqrt ((var (x, 0, 1) + var (y, 0, 1)) / 2);
        otherwise
          v1 = var (x, 0, 1);
          v2 = var (y, 0, 1);
          if (strcmp (vartype, 'equal'))
            v = n1 + n2 - 2;
            s = sqrt (((n1 - 1) * v1 + (n2 - 1) * v2) / v);
          else
            v = delacreDF (n1, n2, v1, v2);
            s = sqrt ((v1 + v2) / 2);
          endif
          e = hedgesJ (v) .* (mean (x, 1) - mean (y, 1)) ./ s;
      endswitch
    case 'glass'
      e = hedgesJ (n1 - 1) * (mean (x, 1) - mean (y, 1)) ./ std (x, 0, 1);
    case 'cliff'
      e = zeros (1, columns (x));
      for k = 1:columns (x)
        if (strcmp (mode, 'paired'))
          e(k) = cliffPaired (x(:,k), y(:,k));
        else
          e(k) = cliffParts (x(:,k), y(:,k));
        endif
      endfor
    case 'mediandiff'
      e = median (x, 1) - median (y, 1);
    case 'robustcohen'
      [t1, w1] = robustParts (x, false);
      if (strcmp (mode, 'one'))
        v = n1 - 1;
        e = 0.642 * (t1 - mu) ./ sqrt (w1);
      else
        [t2, w2] = robustParts (y, false);
        if (strcmp (mode, 'paired'))
          v = n1 - 1;
          s = sqrt ((w1 + w2) / 2);
        elseif (strcmp (vartype, 'equal'))
          v = n1 + n2 - 2;
          s = sqrt (((n1 - 1) * w1 + (n2 - 1) * w2) / v);
        else
          v = delacreDF (n1, n2, w1, w2);
          s = sqrt ((w1 + w2) / 2);
        endif
        e = 0.642 * (t1 - t2) ./ s;
      endif
      e = e .* hedgesJ (v);
    case 'akpcohen'
      ## As akp.effect and D.akp.effect of the R package WRS2 compute it, with
      ## the constant exact where WRS2 integrates it numerically
      c = norminv (0.8);
      c = sqrt (0.6 - 2 * c * normpdf (c) + 0.4 * c ^ 2);
      switch (mode)
        case 'one'
          [t1, w1] = robustParts (x, true);
          e = c * (t1 - mu) ./ sqrt (w1);
        case 'paired'
          [t1, w1] = robustParts (x - y, true);
          e = c * t1 ./ sqrt (w1);
        otherwise
          [t1, w1] = robustParts (x, true);
          [t2, w2] = robustParts (y, true);
          if (strcmp (vartype, 'equal'))
            s = sqrt (((n1 - 1) * w1 + (n2 - 1) * w2) / (n1 + n2 - 2));
          else
            s = sqrt (w1);
          endif
          e = c * (t1 - t2) ./ s;
      endswitch
    case 'kstest'
      e = zeros (1, columns (x));
      for k = 1:columns (x)
        e(k) = ksStat (x(:,k), y(:,k));
      endfor
  endswitch
endfunction

## The noncentrality at which the noncentral t distribution on DF degrees of
## freedom puts T at its P quantile
function lambda = nctBound (t, df, p)
  f = @(l) nctcdf (t, df, l) - p;
  lo = t - 10;
  hi = t + 10;
  while (f (lo) < 0)
    lo -= 2 * (hi - lo);
  endwhile
  while (f (hi) > 0)
    hi += 2 * (hi - lo);
  endwhile
  lambda = fzero (f, [lo, hi], optimset ('TolX', eps));
endfunction

## Both bounds, on the scale of the noncentrality, of the interval that
## inverts the noncentral t distribution of the statistic T
function b = nctInterval (t, df, alpha)
  if (! isfinite (t))
    b = [NaN, NaN];
    return;
  endif
  b = [nctBound(t, df, 1 - alpha / 2), nctBound(t, df, alpha / 2)];
endfunction

## The exact interval of one effect, or NaN with the reason WHY
function [ci, why] = exactInterval (name, x, y, mode, vartype, mu, alpha)
  ci = [NaN, NaN];
  why = '';
  n1 = numel (x);
  n2 = numel (y);
  q = [alpha / 2, 1 - alpha / 2];
  few = 'a sample holds too few observations';
  flat = 'a sample has zero variance';
  switch (name)
    case 'meandiff'
      switch (mode)
        case {'one', 'paired'}
          if (strcmp (mode, 'one'))
            d = x - mu;
          else
            d = x - y;
          endif
          if (n1 < 2)
            why = few;
            return;
          endif
          ci = mean (d) + tinv (q, n1 - 1) * std (d) / sqrt (n1);
        otherwise
          if (strcmp (vartype, 'equal'))
            v = n1 + n2 - 2;
            if (v < 1)
              why = few;
              return;
            endif
            se = sqrt (((n1 - 1) * var (x) + (n2 - 1) * var (y)) / v ...
                       * (1 / n1 + 1 / n2));
          else
            if (n1 < 2 || n2 < 2)
              why = few;
              return;
            endif
            a = var (x) / n1;
            b = var (y) / n2;
            se = sqrt (a + b);
            v = (a + b) ^ 2 / (a ^ 2 / (n1 - 1) + b ^ 2 / (n2 - 1));
          endif
          ci = mean (x) - mean (y) + tinv (q, v) * se;
      endswitch

    case 'cohen'
      switch (mode)
        case 'one'
          v = n1 - 1;
          if (v < 1)
            why = few;
            return;
          elseif (std (x) == 0)
            why = flat;
            return;
          endif
          t = (mean (x) - mu) / (std (x) / sqrt (n1));
          ci = nctInterval (t, v, alpha) * hedgesJ (v) / sqrt (n1);
        case 'paired'
          ## The MAG interval of Cousineau and Goulet-Pelletier (2021)
          if (n1 < 2)
            why = few;
            return;
          endif
          sx = std (x);
          sy = std (y);
          sa = sqrt ((sx ^ 2 + sy ^ 2) / 2);
          rW = corr (x, y) * sx * sy / sa ^ 2;
          if (! (sx > 0 && sy > 0))
            why = flat;
            return;
          elseif (rW >= 1)
            why = 'the paired differences have zero variance';
            return;
          endif
          sc = sqrt (n1 / (2 * (1 - rW)));
          lambda = mean (x - y) / sa * hedgesJ (n1 - 1) ^ 2 * sc;
          ci = nctinv (q, n1 - 1, lambda) / sc;
        otherwise
          if (strcmp (vartype, 'equal'))
            v = n1 + n2 - 2;
            if (v < 1)
              why = few;
              return;
            endif
            sp = sqrt (((n1 - 1) * var (x) + (n2 - 1) * var (y)) / v);
            if (sp == 0)
              why = flat;
              return;
            endif
            r = sqrt (1 / n1 + 1 / n2);
            t = (mean (x) - mean (y)) / (sp * r);
            ci = nctInterval (t, v, alpha) * hedgesJ (v) * r;
          else
            if (n1 < 2 || n2 < 2)
              why = few;
              return;
            elseif (var (x) == 0 || var (y) == 0)
              why = flat;
              return;
            endif
            v = delacreDF (n1, n2, var (x), var (y));
            se = sqrt (var (x) / n1 + var (y) / n2);
            sa = sqrt ((var (x) + var (y)) / 2);
            t = (mean (x) - mean (y)) / se;
            ci = nctInterval (t, v, alpha) * hedgesJ (v) * se / sa;
          endif
      endswitch

    case 'glass'
      if (n1 < 2 || n2 < 2)
        why = few;
        return;
      elseif (var (x) == 0 || var (y) == 0)
        why = flat;
        return;
      endif
      v = n1 - 1;
      se = sqrt (var (x) / n1 + var (y) / n2);
      t = (mean (x) - mean (y)) / se;
      ci = nctInterval (t, v, alpha) * hedgesJ (v) * se / std (x);

    case 'cliff'
      z = norminv (1 - alpha / 2);
      if (strcmp (mode, 'paired'))
        [d, di, dj, S2, S3] = cliffPaired (x, y);
        if (S2 == 0)
          s2 = 0;
        elseif (n1 < 4)
          ## MATLAB returns [d, d] here, a zero width the data cannot support
          why = 'fewer than four pairs leave its variance undefined';
          return;
        else
          s2 = ((n1 - 1) ^ 2 * sum ((di + dj - 2 * d) .^ 2) - S2 - S3) ...
               / (n1 * (n1 - 1) * (n1 - 2) * (n1 - 3));
          s2 = max (s2, 0);
        endif
        ci = d + [-z, z] * sqrt (s2);
      else
        [d, di, dj, S] = cliffParts (x, y);
        if (S == 0)
          s2 = 0;
        elseif (n1 < 2 || n2 < 2)
          why = few;
          return;
        else
          s2 = ((n2 - 1) * var (di) + (n1 - 1) * var (dj) ...
                + S / ((n1 - 1) * (n2 - 1))) / (n1 * n2);
        endif
        if (s2 == 0)
          ci = [d, d];
        else
          c = z * sqrt (s2) * sqrt ((1 - d ^ 2) ^ 2 + z ^ 2 * s2);
          ci = (d - d ^ 3 + [-c, c]) / (1 - d ^ 2 + z ^ 2 * s2);
        endif
      endif
  endswitch
endfunction

%!shared x, y, yp
%! x = [2.1; 3.4; 1.9; 5.6; 4.4; 3.8; 2.7; 6.1; 3.3; 4.9];
%! y = [1.2; 2.8; 0.9; 2.2; 3.1; 1.7; 2.5; 0.4; 1.9; 3.6; 2.0; 1.1];
%! yp = [1.8; 2.9; 2.2; 4.1; 4.9; 2.6; 3.0; 4.8; 2.1; 4.0];
%!test
%! T = meanEffectSize (x, y);
%! assert_equal (T.Properties.VariableNames, {'Effect', 'ConfidenceIntervals'});
%! assert_equal (T.Properties.RowNames, {'MeanDifference'});
%! assert_equal (T.Effect, 1.87, -1e-14);
%! assert_equal (T.ConfidenceIntervals, ...
%!               [0.8093229690413275, 2.9306770309586709], -1e-13);
%!test
%! T = meanEffectSize (x, 'Mean', 3);
%! assert_equal (T{:,:}, [0.8199999999999994, -0.19771934142469916, ...
%!                        1.8377193414246979], -1e-13);
%!test
%! T = meanEffectSize (x, 'Effect', 'cohen');
%! assert_equal (T{:,:}, [2.4538321615715333, 1.1921466698808931, ...
%!                        3.6909287967065962], -1e-13);
%!test
%! T = meanEffectSize (x, y, 'VarianceType', 'unequal');
%! assert_equal (T{:,:}, [1.87, 0.74758396809829564, 2.992416031901703], ...
%!               -1e-13);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'cohen');
%! assert_equal (T.Properties.RowNames, {'CohensD'});
%! assert_equal (T{:,:}, [1.51473229411906, 0.5682402309164547, ...
%!                        2.4328200098432089], -1e-13);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'cohen', 'VarianceType', 'unequal');
%! assert_equal (T{:,:}, [1.4716714731499507, 0.49927861679878149, ...
%!                        2.411323051341753], -1e-13);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'glass');
%! assert_equal (T.Properties.RowNames, {'GlasssDelta'});
%! assert_equal (T{:,:}, [1.2012215031776874, 0.32205950867366817, ...
%!                        2.0417713231164596], -1e-13);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'glass', 'VarianceType', 'unequal');
%! assert_equal (T{:,:}, [1.2012215031776874, 0.32205950867366817, ...
%!                        2.0417713231164596], -1e-13);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'cliff');
%! assert_equal (T.Properties.RowNames, {'CliffsDelta'});
%! assert_equal (T{:,:}, [0.725, 0.27352437233891558, ...
%!                        0.91469527893564417], -1e-13);
%!test
%! T = meanEffectSize ([1; 2; 2; 3; 3; 3; 4; 5], [2; 2; 3; 3; 4; 6], ...
%!                     'Effect', 'cliff');
%! assert_equal (T{:,:}, [-0.14583333333333334, -0.64640108601125412, ...
%!                        0.44249639444771338], -1e-13);
%!test
%! T = meanEffectSize ([5; 6; 7], [1; 2; 3], 'Effect', 'cliff');
%! assert_equal (T{:,:}, [1, 1, 1]);
%!test
%! T = meanEffectSize (1, [3; 4; 5], 'Effect', 'cliff');
%! assert_equal (T{:,:}, [-1, -1, -1]);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'cliff', 'Alpha', 0.1);
%! assert_equal (T.ConfidenceIntervals, ...
%!               [0.35705295319266317, 0.89817708994266032], -1e-13);
%!test
%! T = meanEffectSize (x, y, 'Effect', {'cliff', 'meandiff', 'glass', ...
%!                                      'cohen'});
%! assert_equal (T.Properties.RowNames, {'CliffsDelta'; 'MeanDifference'; ...
%!                                       'GlasssDelta'; 'CohensD'});
%! assert_equal (T.Effect, [0.725; 1.87; 1.2012215031776874; ...
%!                          1.51473229411906], -1e-13);
%!test
%! T = meanEffectSize (x, y, 'Effect', {'cohen', 'cohen'});
%! assert_equal (T.Properties.RowNames, {'CohensD'});
%!test
%! T = meanEffectSize (x, y, 'Effect', 'COH');
%! assert_equal (T.Properties.RowNames, {'CohensD'});
%!test
%! T = meanEffectSize (x, y, 'Effect', 'cohen', 'Mean', 3);
%! assert_equal (T.Effect, 1.51473229411906, -1e-13);
%!test
%! T = meanEffectSize (x', y', 'Effect', 'cohen');
%! assert_equal (T.Effect, 1.51473229411906, -1e-13);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'cohen', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Properties.VariableNames, {'Effect'});
%! assert_equal (T.Effect, 1.51473229411906, -1e-13);
%!test
%! T = meanEffectSize ([x; NaN], [NaN; y], 'Effect', 'cohen');
%! assert_equal (T{:,:}, [1.51473229411906, 0.5682402309164547, ...
%!                        2.4328200098432089], -1e-13);
%!test
%! T = meanEffectSize (single (x), y, 'Effect', 'cohen');
%! assert_equal (class (T.Effect), 'single');
%! assert_equal (class (T.ConfidenceIntervals), 'single');
%!test
%! T = meanEffectSize (x, yp, 'Paired', true);
%! assert_equal (T{:,:}, [0.57999999999999874, 0.044888382079159239, ...
%!                        1.1151116179208382], -1e-13);
%!test
%! T = meanEffectSize (x, yp, 'Paired', 'on', 'Effect', 'cohen');
%! assert_equal (T{:,:}, [0.41222519124990314, 0.016525506916817929, ...
%!                        0.94928510130451305], -1e-12);
%!test
%! T = meanEffectSize (x, yp, 'Paired', true, 'Effect', 'cohen', ...
%!                     'Alpha', 0.01);
%! assert_equal (T.ConfidenceIntervals, ...
%!               [-0.10395569010229536, 1.2369598804778004], -1e-12);
%!test
%! T = meanEffectSize ((1:8)', [8.2; 6.9; 7.4; 5.1; 5.8; 3.6; 4.4; 2.3], ...
%!                     'Paired', true, 'Effect', 'cohen');
%! assert_equal (T{:,:}, [-0.38188801178758269, -2.124735313427935, ...
%!                        1.1823775856056995], -1e-12);
%!test
%! T = meanEffectSize ([1; 2], [2; 4], 'Paired', true, 'Effect', 'cohen');
%! assert_equal (T{:,:}, [0, -5.6823875052232813, 5.6823875052232813], ...
%!               -1e-12);
%!test
%! T = meanEffectSize ([x; NaN], [yp; 2], 'Paired', true, 'Effect', 'cohen');
%! assert_equal (T.Effect, 0.41222519124990314, -1e-13);
%!test
%! T = meanEffectSize (x, yp, 'Paired', true, 'Effect', 'cliff');
%! assert_equal (T{:,:}, [0.22222222222222221, 0.11569417468908846, ...
%!                        0.32875026975535593], -1e-13);
%!test
%! T = meanEffectSize ([1; 2; 3; 4], [3; 1; 2; 5], 'Paired', true, ...
%!                     'Effect', 'cliff');
%! assert_equal (T{:,:}, [-0.083333333333333329, -0.24666366537833756, ...
%!                        0.079996998711670889], -1e-13);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'mediandiff', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Properties.RowNames, {'MedianDifference'});
%! assert_equal (T.Effect, 1.65, -1e-14);
%!test
%! T = meanEffectSize (x, yp, 'Paired', true, 'Effect', 'mediandiff', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Effect, 0.65, -1e-14);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'kstest', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Properties.RowNames, {'KolmogorovSmirnovStatistic'});
%! assert_equal (T.Effect, 0.6166666666666667, -1e-14);
%!test
%! T = meanEffectSize (x, yp, 'Paired', true, 'Effect', 'kstest', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Effect, 0.3, -1e-14);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'robustcohen', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Properties.RowNames, {'RobustCohensD'});
%! assert_equal (T.Effect, 1.2360631043171186, -1e-13);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'robustcohen', ...
%!                     'VarianceType', 'unequal', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Effect, 1.1952733752449272, -1e-13);
%!test
%! T = meanEffectSize (x, yp, 'Paired', true, 'Effect', 'robustcohen', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Effect, 0.31641941202232715, -1e-13);
%!test
%! T = meanEffectSize (x, 'Effect', 'robustcohen', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Effect, 1.9756476082746626, -1e-13);
## akpcohen against WRS2 1.1.7, whose constant is integrated numerically
## and differs from the exact one by 1.2e-9 relative
%!test
%! T = meanEffectSize (x, y, 'Effect', 'akpcohen', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Properties.RowNames, {'AKPCohensD'});
%! assert_equal (T.Effect, 1.4340672042318023, -2e-9);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'akpcohen', 'VarianceType', ...
%!                     'unequal', 'ConfidenceIntervalType', 'none');
%! assert_equal (T.Effect, 1.2409782758085974, -2e-9);
%!test
%! T = meanEffectSize ([1; 2; 3; 4; 5; 6; 7; 8; 20; 30], 2 * (1:10)', ...
%!                     'Effect', 'akpcohen', 'ConfidenceIntervalType', 'none');
%! assert_equal (T.Effect, -1.0275756121741466, -2e-9);
%!test
%! T = meanEffectSize (x, 'Effect', 'akpcohen', 'Mean', 3, ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Effect, 0.50999107225010853, -2e-9);
%!test
%! T = meanEffectSize ([1; 2; 3; 4; 5; 6; 7; 8; 20; 30], ...
%!                     'Effect', 'akpcohen', 'ConfidenceIntervalType', 'none');
%! assert_equal (T.Effect, 1.6247397012560751, -2e-9);
%!test
%! T = meanEffectSize (x, yp, 'Paired', true, 'Effect', 'akpcohen', ...
%!                     'ConfidenceIntervalType', 'none');
%! assert_equal (T.Effect, 0.60651611028436192, -2e-9);
%!test
%! T = meanEffectSize (x, y, 'Effect', 'cohen', 'NumBootstraps', 10);
%! assert_equal (T.ConfidenceIntervals, ...
%!               [0.5682402309164547, 2.4328200098432089], -1e-13);
%!test
%! rand ('state', 1);
%! T = meanEffectSize ([1; 2; 3], [10; 10; 10; 10], 'NumBootstraps', 5000, ...
%!                     'ConfidenceIntervalType', 'bootstrap');
%! assert_equal (T.ConfidenceIntervals, [-9, -7]);
%!test
%! rand ('state', 1);
%! T = meanEffectSize ([1; 2; 3], [10; 10; 10; 10], 'NumBootstraps', 5000, ...
%!                     'ConfidenceIntervalType', 'bootstrap', 'Alpha', 0.4);
%! assert_equal (T.ConfidenceIntervals, [-8.5, -7.5]);
%!test
%! rand ('state', 1);
%! T = meanEffectSize ([1; 2; 3], [10; 10; 10; 10], 'NumBootstraps', 5000, ...
%!                     'ConfidenceIntervalType', 'bootstrap', 'Alpha', 0.4, ...
%!                     'Resampling', 'stratified');
%! assert_equal (T.ConfidenceIntervals, [-25, -23] / 3, -1e-14);
%!test
%! rand ('state', 1);
%! T = meanEffectSize ([10; 10; 10], [1; 2; 3; 4], 'NumBootstraps', 5000, ...
%!                     'ConfidenceIntervalType', 'bootstrap', 'Alpha', 0.4, ...
%!                     'BootstrapOptions', statset ('UseParallel', true));
%! assert_equal (T.ConfidenceIntervals, [7, 8]);
%!warning<meanEffectSize: the confidence interval of effect 'cohen' is NaN, since a sample has zero variance.> ...
%! T = meanEffectSize (ones (5, 1), 2 * ones (6, 1), 'Effect', 'cohen');
%! assert_equal (T{:,:}, [-Inf, NaN, NaN]);
%!warning<meanEffectSize: the confidence interval of effect 'cohen' is NaN, since a sample holds too few observations.> ...
%! T = meanEffectSize (5, 'Effect', 'cohen');
%! assert_equal (T{:,:}, [NaN, NaN, NaN]);
%!warning<meanEffectSize: the confidence interval of effect 'cohen' is NaN, since the paired differences have zero variance.> ...
%! T = meanEffectSize ((1:8)', (2:9)', 'Paired', true, 'Effect', 'cohen');
%! assert_equal (T.ConfidenceIntervals, [NaN, NaN]);
%!warning<meanEffectSize: the confidence interval of effect 'cliff' is NaN, since fewer than four pairs leave its variance undefined.> ...
%! T = meanEffectSize ([1; 2; 3], [3; 1; 2], 'Paired', true, 'Effect', 'cliff');
%! assert_equal (T{:,:}, [-1 / 6, NaN, NaN], -1e-14);

%!error<meanEffectSize: too few input arguments.> meanEffectSize ()
%!error<meanEffectSize: X must be a nonempty vector of type double or single.> ...
%! meanEffectSize ([1, 2; 3, 4])
%!error<meanEffectSize: X must be a nonempty vector of type double or single.> ...
%! meanEffectSize ([])
%!error<meanEffectSize: X must be a nonempty vector of type double or single.> ...
%! meanEffectSize (int32 ([1; 2; 3]))
%!error<meanEffectSize: Y must be a nonempty vector of type double or single.> ...
%! meanEffectSize ([1; 2; 3], true (3, 1))
%!error<meanEffectSize: invalid optional paired argument.> ...
%! meanEffectSize ([1; 2; 3], 'Foo', 1)
%!error<meanEffectSize: 'Effect' must be one or more of 'meandiff', 'cohen', 'glass', 'cliff', 'mediandiff', 'robustcohen', 'akpcohen', 'kstest'.> ...
%! meanEffectSize ([1; 2; 3], 'Effect', 'foo')
%!error<meanEffectSize: 'Effect' must be one or more of 'meandiff', 'cohen', 'glass', 'cliff', 'mediandiff', 'robustcohen', 'akpcohen', 'kstest'.> ...
%! meanEffectSize ([1; 2; 3], 'Effect', 'c')
%!error<meanEffectSize: 'Effect' must be one or more of 'meandiff', 'cohen', 'glass', 'cliff', 'mediandiff', 'robustcohen', 'akpcohen', 'kstest'.> ...
%! meanEffectSize ([1; 2; 3], 'Effect', {})
%!error<meanEffectSize: 'Mean' must be a real scalar.> ...
%! meanEffectSize ([1; 2; 3], 'Mean', [1, 2])
%!error<meanEffectSize: 'Paired' must be true, false, 'on' or 'off'.> ...
%! meanEffectSize ([1; 2; 3], [1; 2; 3], 'Paired', 2)
%!error<meanEffectSize: 'VarianceType' must be 'equal' or 'unequal'.> ...
%! meanEffectSize ([1; 2; 3], [1; 2; 3], 'VarianceType', 'foo')
%!error<meanEffectSize: 'Alpha' must be a scalar between 0 and 1.> ...
%! meanEffectSize ([1; 2; 3], 'Alpha', 1)
%!error<meanEffectSize: 'ConfidenceIntervalType' must be 'exact', 'bootstrap' or 'none'.> ...
%! meanEffectSize ([1; 2; 3], 'ConfidenceIntervalType', 'foo')
%!error<meanEffectSize: 'NumBootstraps' must be a positive integer.> ...
%! meanEffectSize ([1; 2; 3], 'NumBootstraps', 2.5)
%!error<meanEffectSize: 'Resampling' must be 'pooled' or 'stratified'.> ...
%! meanEffectSize ([1; 2; 3], 'Resampling', 'foo')
%!error<meanEffectSize: 'BootstrapOptions' must be a structure.> ...
%! meanEffectSize ([1; 2; 3], 'BootstrapOptions', 1)
%!error<meanEffectSize: 'Paired' needs two samples.> ...
%! meanEffectSize ([1; 2; 3], 'Paired', true)
%!error<meanEffectSize: effect 'glass' needs two samples.> ...
%! meanEffectSize ([1; 2; 3], 'Effect', 'glass')
%!error<meanEffectSize: paired samples X and Y must have the same number of elements.> ...
%! meanEffectSize ([1; 2; 3], [1; 2], 'Paired', true)
%!error<meanEffectSize: paired samples take 'VarianceType' 'equal' only.> ...
%! meanEffectSize ([1; 2; 3], [3; 2; 1], 'Paired', true, ...
%!                 'VarianceType', 'unequal')
%!error<meanEffectSize: effect 'glass' does not apply to paired samples.> ...
%! meanEffectSize ([1; 2; 3], [3; 2; 1], 'Paired', true, 'Effect', 'glass')
%!error<meanEffectSize: effect 'mediandiff' has no exact confidence interval; use 'bootstrap'.> ...
%! meanEffectSize ([1; 2; 3], [3; 2; 1], 'Effect', 'mediandiff', ...
%!                 'ConfidenceIntervalType', 'exact')
