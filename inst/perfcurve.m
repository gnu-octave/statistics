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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {[@var{X}, @var{Y}] =} perfcurve (@var{labels}, @var{scores}, @var{posclass})
## @deftypefnx {statistics} {[@var{X}, @var{Y}, @var{T}, @var{AUC}, @var{OPTROCPT}] =} perfcurve (@dots{})
## @deftypefnx {statistics} {[@dots{}] =} perfcurve (@dots{}, @var{Name}, @var{Value})
##
## Receiver operating characteristic (ROC) and other classifier performance
## curves.
##
## @code{[@var{X}, @var{Y}] = perfcurve (@var{labels}, @var{scores},
## @var{posclass})} returns the ROC curve for the classifier scores in
## @var{scores} given the true class @var{labels} and the positive class
## @var{posclass}.  @var{labels} is a numeric vector or a cell array of
## character vectors; @var{scores} is a numeric vector of the same length, where
## larger values indicate stronger evidence for the positive class.  By default
## @var{X} is the false positive rate and @var{Y} the true positive rate.
##
## @code{[@var{X}, @var{Y}, @var{T}, @var{AUC}, @var{OPTROCPT}] = perfcurve
## (@dots{})} also returns the thresholds @var{T} on the scores, the area
## @var{AUC} under the (@var{X}, @var{Y}) curve, and the optimal operating point
## @var{OPTROCPT} @code{= [FPR, TPR]} of the ROC curve.
##
## The following @qcode{Name-Value} pairs are supported:
##
## @multitable @columnfractions 0.18 0.82
## @headitem @var{Name} @tab @var{Value}
##
## @item @qcode{'XCrit'} @tab The criterion for @var{X}.  The default is
## @qcode{'FPR'}.
##
## @item @qcode{'YCrit'} @tab The criterion for @var{Y}.  The default is
## @qcode{'TPR'}.  Supported criteria are @qcode{'TPR'} (@qcode{'sens'},
## @qcode{'reca'}), @qcode{'FNR'}, @qcode{'FPR'} (@qcode{'fall'}), @qcode{'TNR'}
## (@qcode{'spec'}), @qcode{'PPV'} (@qcode{'prec'}), @qcode{'NPV'},
## @qcode{'accu'}, the counts @qcode{'TP'}, @qcode{'FN'}, @qcode{'FP'},
## @qcode{'TN'}, and the rates @qcode{'RPP'}, @qcode{'RNP'}.
##
## @item @qcode{'NegClass'} @tab The negative class(es).  The default,
## @qcode{'all'}, treats every label other than @var{posclass} as negative.
##
## @item @qcode{'Weights'} @tab A vector of non-negative observation weights.
##
## @item @qcode{'Cost'} @tab A @math{2*2} misclassification-cost matrix
## @code{[C(P|P) C(N|P); C(P|N) C(N|N)]} used for @var{OPTROCPT}.  The default
## is @code{[0 1; 1 0]}.
##
## @item @qcode{'XVals'} @tab Values of the @var{X} criterion at which to return
## the curve.  @qcode{'TVals'} does the same for the thresholds.
##
## @item @qcode{'ProcessNaN'} @tab How to treat @code{NaN} scores:
## @qcode{'ignore'} (default) or @qcode{'addtofalse'}.
##
## @item @qcode{'NBoot'} @tab Number of bootstrap replicates for confidence
## bounds on @var{Y} and @var{AUC}.  The default @math{0} computes no bounds.
##
## @item @qcode{'BootType'} @tab The bootstrap interval: @qcode{'bca'} (default,
## bias-corrected and accelerated), @qcode{'percentile'}, or @qcode{'normal'}.
##
## @item @qcode{'Alpha'} @tab The significance level for the bounds, so the
## confidence level is @math{1 - @var{Alpha}}.  The default is @math{0.05}.
## @end multitable
##
## With @qcode{'NBoot'} greater than zero, @var{Y} is returned as an
## @math{m*3} array @code{[@var{Y}, @var{Ylow}, @var{Yhigh}]} and @var{AUC} as
## @code{[@var{AUC}, @var{AUClow}, @var{AUChigh}]}.  The bootstrap uses an
## independent random stream, so the bounds do not match MATLAB numerically.
## @code{[@dots{}, @var{SUBY}, @var{SUBYNAMES}] = perfcurve (@dots{})} returns
## the @var{Y} values for each negative subclass and their names.
##
## When called with no output arguments the curve is plotted.
##
## @seealso{fitcsvm, fitcknn, glmfit}
## @end deftypefn

function [X, Y, T, AUC, OPTROCPT, SUBY, SUBYNAMES] = perfcurve (labels, ...
                                             scores, posclass, varargin)

  if (nargin < 3)
    print_usage ();
  endif
  if (! (isnumeric (scores) && isreal (scores) && isvector (scores)))
    error ("perfcurve: SCORES must be a vector of real values.");
  endif
  scores = scores(:);
  n = numel (scores);
  if (! ((isnumeric (labels) || iscellstr (labels) || ischar (labels))
         && numel_labels (labels) == n))
    error ("perfcurve: LABELS must have one element per score.");
  endif

  ## Defaults and Name-Value parsing.
  xcrit = "FPR";
  ycrit = "TPR";
  negclass = "all";
  weights = [];
  cost = [0 1; 1 0];
  xvals = [];
  tvals = [];
  usenearest = true;
  prior = "empirical";
  processnan = "ignore";
  nboot = 0;
  boottype = "bca";
  alpha = 0.05;
  if (mod (numel (varargin), 2) != 0)
    error ("perfcurve: optional arguments must be Name-Value pairs.");
  endif
  for k = 1:2:numel (varargin)
    if (! ischar (varargin{k}))
      error ("perfcurve: parameter names must be character vectors.");
    endif
    switch (lower (varargin{k}))
      case "xcrit"
        xcrit = varargin{k+1};
      case "ycrit"
        ycrit = varargin{k+1};
      case "negclass"
        negclass = varargin{k+1};
      case "weights"
        weights = varargin{k+1};
      case "cost"
        cost = varargin{k+1};
      case "xvals"
        xvals = varargin{k+1};
      case "tvals"
        tvals = varargin{k+1};
      case "usenearest"
        uv = varargin{k+1};
        usenearest = strcmpi (uv, "on") || isequal (uv, true);
      case "prior"
        prior = varargin{k+1};
      case "processnan"
        processnan = lower (varargin{k+1});
      case "nboot"
        nboot = varargin{k+1};
      case "boottype"
        boottype = lower (varargin{k+1});
      case "alpha"
        alpha = varargin{k+1};
      otherwise
        error ("perfcurve: unknown parameter name '%s'.", varargin{k});
    endswitch
  endfor
  if (! isempty (xvals) && ! isempty (tvals))
    error ("perfcurve: 'XVals' and 'TVals' cannot be used together.");
  endif
  if (! (isnumeric (cost) && isequal (size (cost), [2, 2])))
    error ("perfcurve: 'Cost' must be a 2-by-2 matrix.");
  endif
  if (! (isnumeric (nboot) && isscalar (nboot) && nboot >= 0
         && nboot == fix (nboot)))
    error ("perfcurve: 'NBoot' must be a non-negative integer.");
  endif
  if (! (isnumeric (alpha) && isscalar (alpha) && alpha > 0 && alpha < 1))
    error ("perfcurve: 'Alpha' must be a scalar in the range (0,1).");
  endif

  ## Positive / negative membership (one-vs-all by default).
  ispos = label_match (labels, posclass);
  if (ischar (negclass) && strcmpi (negclass, "all"))
    isneg = ! ispos;
  else
    isneg = label_match (labels, negclass);
  endif
  SUBYNAMES = negclass_names (labels, isneg, negclass);

  ## Weights.
  if (isempty (weights))
    w = ones (n, 1);
  else
    if (! (isnumeric (weights) && isreal (weights) && isvector (weights)
           && numel (weights) == n && all (weights >= 0)))
      error ("perfcurve: 'Weights' must be a non-negative vector with one element per score.");
    endif
    w = weights(:);
  endif

  ## NaN score handling.
  isnanscore = isnan (scores);
  if (any (isnanscore))
    switch (processnan)
      case "ignore"
        keep = ! isnanscore;
        scores = scores(keep);  ispos = ispos(keep);
        isneg = isneg(keep);    w = w(keep);
      case "addtofalse"
        scores(isnanscore) = -Inf;    ## always classified negative
      otherwise
        error ("perfcurve: unrecognised 'ProcessNaN' value.");
    endswitch
  endif

  ## Restrict to labelled (positive or negative) observations.
  use = ispos | isneg;
  scores = scores(use);  ispos = ispos(use);  isneg = isneg(use);  w = w(use);
  if (! any (ispos) || ! any (isneg))
    error ("perfcurve: LABELS must contain both positive and negative classes.");
  endif

  ## Optionally reweight the observations to a requested prior.
  [~, ~, w] = apply_prior (prior, sum (w(ispos)), sum (w(isneg)), w, ...
                           ispos, isneg);

  ## Full-sample curve for the requested criteria plus the underlying ROC.
  [X, Y, T, tpr, fpr, P, N] = build_curve (scores, ispos, isneg, w, ...
                                           xcrit, ycrit);

  ## Area under the (X, Y) curve and the cost-optimal ROC operating point.
  AUC = trapz (X, Y);
  ecost = P * (1 - tpr) * cost(1,2) + N * fpr * cost(2,1);
  best = find (abs (ecost - min (ecost)) < 1e-12);
  [~, bk] = min (fpr(best));
  OPTROCPT = [fpr(best(bk)), tpr(best(bk))];

  ## Restrict the returned curve to requested X or threshold values.
  if (! isempty (xvals))
    idx = select_xvals (X, Y, xvals(:), usenearest);
  elseif (! isempty (tvals))
    idx = select_tvals (T, tvals(:));
  else
    idx = (1:numel (X))';
  endif
  X = X(idx);  Y = Y(idx);  T = T(idx);
  SUBY = Y;

  ## Bootstrap confidence bounds on Y and AUC (self-contained resampler).
  if (nboot > 0)
    [Ylo, Yhi, alo, ahi] = perfcurve_boot_ (scores, ispos, isneg, w, xcrit, ...
                             ycrit, X, Y, AUC, nboot, boottype, alpha);
    Y = [Y, Ylo, Yhi];
    SUBY = Y;
    AUC = [AUC, alo, ahi];
  endif

  ## Plot when no output is requested.
  if (nargout == 0)
    plot (X, Y(:,1));
    xlabel (xcrit);  ylabel (ycrit);
    clear X Y T AUC OPTROCPT SUBY SUBYNAMES
  endif

endfunction

## Build the full performance curve for a data subset: returns the X and Y
## criteria, thresholds T, the underlying ROC rates, and the class totals.
function [X, Y, T, tpr, fpr, P, N] = build_curve (scores, ispos, isneg, w, ...
                                                  xcrit, ycrit)
  P = sum (w(ispos));
  N = sum (w(isneg));
  [ss, ord] = sort (scores, "descend");
  pw = cumsum (w(ord) .* ispos(ord));
  nw = cumsum (w(ord) .* isneg(ord));
  chg = find (diff (ss) != 0);
  last = [chg; numel(ss)];               ## last index of each score group
  tp = [0; pw(last)];
  fp = [0; nw(last)];
  T = [ss(1); ss(last)];
  X = criterion (xcrit, tp, fp, P, N);
  Y = criterion (ycrit, tp, fp, P, N);
  tpr = tp / P;
  fpr = fp / N;
endfunction

## Number of labels in a numeric vector or cell array of character vectors.
function m = numel_labels (labels)
  if (ischar (labels))
    m = rows (labels);
  else
    m = numel (labels);
  endif
endfunction

## Logical membership of LABELS in CLS (numeric equality or string match).
function tf = label_match (labels, cls)
  if (isnumeric (labels))
    tf = ismember (labels(:), cls(:));
  else
    if (ischar (labels))
      labels = cellstr (labels);
    endif
    if (ischar (cls))
      cls = {cls};
    endif
    tf = ismember (labels(:), cls(:));
  endif
endfunction

## Rescale weights so the class totals follow a requested prior.
function [P, N, w] = apply_prior (prior, P, N, w, ispos, isneg)
  if (ischar (prior))
    switch (lower (prior))
      case "empirical"
        return;
      case "uniform"
        pn = [0.5, 0.5];
      otherwise
        error ("perfcurve: unrecognised 'Prior' value.");
    endswitch
  elseif (isnumeric (prior) && numel (prior) == 2 && all (prior > 0))
    pn = prior(:)' / sum (prior);
  else
    error ("perfcurve: 'Prior' must be 'empirical', 'uniform', or [P N].");
  endif
  w(ispos) = w(ispos) * (pn(1) / P);
  w(isneg) = w(isneg) * (pn(2) / N);
  P = pn(1);  N = pn(2);
endfunction

## Evaluate a named performance criterion from the count arrays.
function v = criterion (name, tp, fp, P, N)
  tn = N - fp;
  fn = P - tp;
  total = P + N;
  pp = tp + fp;
  switch (lower (name))
    case {"tpr", "sens", "reca", "recall"}
      v = tp / P;
    case {"fnr", "miss"}
      v = fn / P;
    case {"fpr", "fall"}
      v = fp / N;
    case {"tnr", "spec"}
      v = tn / N;
    case {"ppv", "prec", "precision"}
      v = tp ./ pp;
    case "npv"
      v = tn ./ (tn + fn);
    case {"accu", "acc"}
      v = (tp + tn) / total;
    case "tp"
      v = tp;
    case "fn"
      v = fn;
    case "fp"
      v = fp;
    case "tn"
      v = tn;
    case "rpp"
      v = pp / total;
    case "rnp"
      v = (tn + fn) / total;
    otherwise
      error ("perfcurve: unsupported criterion '%s'.", name);
  endswitch
endfunction

## Indices of the curve returned for requested X values: the (0,0) start plus,
## for each value, the upper-envelope point at the largest actual X not above.
function idx = select_xvals (X, Y, xvals, usenearest)
  idx = 1;
  for v = xvals'
    if (usenearest)
      cand = find (X <= v + 1e-12);
      if (isempty (cand))
        cand = 1;
      endif
      xstar = max (X(cand));
      here = find (abs (X - xstar) < 1e-12);
      idx(end+1) = here(end);           ## upper envelope (largest Y at this X)
    else
      idx(end+1) = interp_index (X, v);
    endif
  endfor
  idx = idx(:);
endfunction

## Indices for requested threshold values: nearest actual threshold (Inf -> the
## first, "reject all" point).
function idx = select_tvals (T, tvals)
  idx = zeros (numel (tvals), 1);
  for k = 1:numel (tvals)
    if (tvals(k) >= T(1))
      idx(k) = 1;
    else
      [~, idx(k)] = min (abs (T - tvals(k)));
    endif
  endfor
endfunction

## Nearest index (used for UseNearest 'off' fallback).
function i = interp_index (X, v)
  [~, i] = min (abs (X - v));
endfunction

## Names of the negative class(es), for SUBYNAMES.
function names = negclass_names (labels, isneg, negclass)
  if (ischar (negclass) && strcmpi (negclass, "all"))
    if (isnumeric (labels))
      names = arrayfun (@(v) num2str (v), unique (labels(isneg)), ...
                        "UniformOutput", false);
    else
      if (ischar (labels))
        labels = cellstr (labels);
      endif
      names = unique (labels(isneg));
    endif
  elseif (isnumeric (negclass))
    names = arrayfun (@(v) num2str (v), negclass(:), "UniformOutput", false);
  elseif (ischar (negclass))
    names = {negclass};
  else
    names = negclass(:);
  endif
endfunction

## Nonparametric bootstrap confidence bounds on the Y criterion (at the returned
## reference X values REFX, point estimates REFY) and on the AUC.  Observations
## are resampled with replacement; percentile or bias-corrected-and-accelerated
## (BCa) bounds are formed.  Bounds do not match MATLAB numerically because the
## bootstrap uses an independent random stream.
function [Ylo, Yhi, alo, ahi] = perfcurve_boot_ (scores, ispos, isneg, w, ...
                                  xcrit, ycrit, refX, refY, aucHat, ...
                                  nboot, boottype, alpha)
  n = numel (scores);
  m = numel (refX);
  Yb = zeros (nboot, m);
  aucB = zeros (nboot, 1);
  for b = 1:nboot
    ii = resample_both (n, ispos, isneg);
    [Xr, Yr] = build_curve (scores(ii), ispos(ii), isneg(ii), w(ii), ...
                            xcrit, ycrit);
    Yb(b,:) = interp_curve (Xr, Yr, refX);
    aucB(b) = trapz (Xr, Yr);
  endfor
  switch (boottype)
    case {"percentile", "per"}
      [Ylo, Yhi] = ci_percentile (Yb, alpha);
      [alo, ahi] = ci_percentile (aucB, alpha);
    case "bca"
      Yj = zeros (n, m);
      aucJ = zeros (n, 1);
      valid = true (n, 1);
      for i = 1:n
        keep = true (n, 1);
        keep(i) = false;
        if (! any (ispos(keep)) || ! any (isneg(keep)))
          valid(i) = false;
          continue;
        endif
        [Xj, Yjc] = build_curve (scores(keep), ispos(keep), isneg(keep), ...
                                 w(keep), xcrit, ycrit);
        Yj(i,:) = interp_curve (Xj, Yjc, refX);
        aucJ(i) = trapz (Xj, Yjc);
      endfor
      [Ylo, Yhi] = ci_bca (Yb, Yj(valid,:), refY(:)', alpha);
      [alo, ahi] = ci_bca (aucB, aucJ(valid), aucHat, alpha);
    case {"normal", "norm"}
      z = norminv (1 - alpha / 2);
      Ylo = (mean (Yb) - z * std (Yb))';
      Yhi = (mean (Yb) + z * std (Yb))';
      alo = mean (aucB) - z * std (aucB);
      ahi = mean (aucB) + z * std (aucB);
    otherwise
      error ("perfcurve: unrecognised 'BootType' value.");
  endswitch
  Ylo = Ylo(:);  Yhi = Yhi(:);
endfunction

## Resample N observation indices with replacement, ensuring both classes are
## present in the draw.
function ii = resample_both (n, ispos, isneg)
  ii = randi (n, n, 1);
  tries = 0;
  while ((! any (ispos(ii)) || ! any (isneg(ii))) && tries < 100)
    ii = randi (n, n, 1);
    tries++;
  endwhile
endfunction

## Evaluate a curve (X, Y) at query abscissae XQ using its upper envelope (max Y
## at each distinct X) and linear interpolation, clamped to the endpoints.
function yq = interp_curve (X, Y, xq)
  [xu, ~, ic] = unique (X(:));
  yu = accumarray (ic, Y(:), [], @max);
  yq = interp1 (xu, yu, xq(:), "linear");
  yq(xq(:) <= xu(1)) = yu(1);
  yq(xq(:) >= xu(end)) = yu(end);
  yq = yq';
endfunction

## Percentile confidence bounds (over rows) of the bootstrap replicates B.
function [lo, hi] = ci_percentile (B, alpha)
  lo = quantile (B, alpha / 2, 1)';
  hi = quantile (B, 1 - alpha / 2, 1)';
endfunction

## Bias-corrected and accelerated (BCa) confidence bounds.  B is nboot-by-m,
## J the jackknife replicates (n-by-m), THETAHAT the point estimate (1-by-m).
function [lo, hi] = ci_bca (B, J, thetahat, alpha)
  nb = rows (B);
  m = columns (B);
  frac = sum (B < thetahat, 1) / nb;
  frac = min (max (frac, 1 / (nb + 1)), nb / (nb + 1));
  z0 = norminv (frac);
  d = mean (J, 1) - J;
  a = sum (d .^ 3, 1) ./ (6 * (sum (d .^ 2, 1)) .^ 1.5);
  a(! isfinite (a)) = 0;
  za = norminv (alpha / 2);
  zb = norminv (1 - alpha / 2);
  a1 = normcdf (z0 + (z0 + za) ./ (1 - a .* (z0 + za)));
  a2 = normcdf (z0 + (z0 + zb) ./ (1 - a .* (z0 + zb)));
  a1(! isfinite (a1)) = alpha / 2;
  a2(! isfinite (a2)) = 1 - alpha / 2;
  lo = zeros (m, 1);
  hi = zeros (m, 1);
  for j = 1:m
    lo(j) = quantile (B(:,j), a1(j));
    hi(j) = quantile (B(:,j), a2(j));
  endfor
endfunction

%!demo
%! ## ROC curve for scores with a known positive class
%! scores = [0.9 0.8 0.7 0.6 0.55 0.5 0.4 0.3 0.2 0.1];
%! labels = [1 1 0 1 0 1 0 0 1 0];
%! [X, Y, T, AUC] = perfcurve (labels, scores, 1);
%! plot (X, Y, "b-o");  xlabel ("FPR");  ylabel ("TPR");
%! title (sprintf ("ROC curve (AUC = %.2f)", AUC));

%!shared labels, scores
%! labels = [1 1 0 1 0 1 0 0 1 0];
%! scores = [0.95 0.9 0.85 0.7 0.7 0.6 0.55 0.4 0.35 0.2];

%!test  # MATLAB parity: default ROC curve, thresholds, AUC and OPTROCPT
%! [X, Y, T, AUC, OPT] = perfcurve (labels, scores, 1);
%! assert (X, [0 0 0 0.2 0.4 0.4 0.6 0.8 0.8 1]', 1e-12);
%! assert (Y, [0 0.2 0.4 0.4 0.6 0.8 0.8 0.8 1 1]', 1e-12);
%! assert (T, [0.95 0.95 0.9 0.85 0.7 0.6 0.55 0.4 0.35 0.2]', 1e-12);
%! assert (AUC, 0.7, 1e-12);
%! assert (OPT, [0 0.4], 1e-12);

%!test  # MATLAB parity: precision-recall and its NaN at the origin
%! [Xr, Yr] = perfcurve (labels, scores, 1, "XCrit", "reca", "YCrit", "prec");
%! assert (Xr, [0 0.2 0.4 0.4 0.6 0.8 0.8 0.8 1 1]', 1e-12);
%! assert (Yr(2:end), [1 1 2/3 0.6 2/3 4/7 0.5 5/9 0.5]', 1e-12);
%! assert (isnan (Yr(1)));

%!test  # MATLAB parity: accuracy and specificity/sensitivity criteria
%! [~, Ya] = perfcurve (labels, scores, 1, "YCrit", "accu");
%! assert (Ya, [0.5 0.6 0.7 0.6 0.6 0.7 0.6 0.5 0.6 0.5]', 1e-12);
%! [Xs, Ys] = perfcurve (labels, scores, 1, "XCrit", "spec", "YCrit", "sens");
%! assert (Xs, [1 1 1 0.8 0.6 0.6 0.4 0.2 0.2 0]', 1e-12);
%! assert (Ys, [0 0.2 0.4 0.4 0.6 0.8 0.8 0.8 1 1]', 1e-12);

%!test  # MATLAB parity: weighted curve, AUC and OPTROCPT
%! w = [2 1 1 1 1 1 1 1 1 2];
%! [Xw, Yw, Tw, AUCw, OPTw] = perfcurve (labels, scores, 1, "Weights", w);
%! assert (Xw, [0 0 0 1 2 2 3 4 4 6]'/6, 1e-12);
%! assert (Yw, [0 2 3 3 4 5 5 5 6 6]'/6, 1e-12);
%! assert (AUCw, 0.791666666666667, 1e-12);
%! assert (OPTw, [0 0.5], 1e-12);

%!test  # MATLAB parity: OPTROCPT shifts with a non-default cost matrix
%! [~, ~, ~, ~, OPTc] = perfcurve (labels, scores, 1, "Cost", [0 1; 2 0]);
%! assert (OPTc, [0 0.4], 1e-12);

%!test  # MATLAB parity: TVals selects nearest thresholds
%! [Xt, Yt, Tt] = perfcurve (labels, scores, 1, "TVals", [Inf 0.8 0.6 0.4 0.2]);
%! assert (Xt, [0 0.2 0.4 0.8 1]', 1e-12);
%! assert (Yt, [0 0.4 0.8 0.8 1]', 1e-12);
%! assert (Tt, [0.95 0.85 0.6 0.4 0.2]', 1e-12);

%!test  # MATLAB parity: XVals selects floor points plus the origin
%! [Xv, Yv, Tv] = perfcurve (labels, scores, 1, "XVals", [0 0.25 0.5 0.75 1]);
%! assert (Xv, [0 0 0.2 0.4 0.6 1]', 1e-12);
%! assert (Yv, [0 0.4 0.4 0.8 0.8 1]', 1e-12);
%! assert (Tv, [0.95 0.9 0.85 0.6 0.55 0.2]', 1e-12);

%!test  # MATLAB parity: cell-array string labels
%! labels2 = {"g","g","b","g","b","g","b","b","g","b"};
%! [X, Y, T, AUC] = perfcurve (labels2, scores, "g");
%! assert (AUC, 0.7, 1e-12);
%! assert (X, [0 0 0 0.2 0.4 0.4 0.6 0.8 0.8 1]', 1e-12);

%!test  # AUC of a perfectly separable problem is 1
%! sc = [0.9 0.8 0.7 0.6 0.3 0.2 0.1 0.05];
%! lb = [1 1 1 1 0 0 0 0];
%! [~, ~, ~, auc] = perfcurve (lb, sc, 1);
%! assert (auc, 1, 1e-12);

%!test  # bootstrap: shapes, point estimate preserved, bounds bracket it (BCa)
%! rand ("seed", 42);  randn ("seed", 42);
%! [X, Y, T, AUC, OPT, SUBY, SUBYN] = ...
%!   perfcurve (labels, scores, 1, "NBoot", 200);
%! assert (columns (Y), 3);
%! assert (size (AUC), [1, 3]);
%! assert (Y(:,1), [0 0.2 0.4 0.4 0.6 0.8 0.8 0.8 1 1]', 1e-12);
%! assert (all (Y(:,2) <= Y(:,1) + 1e-9));                       # lower <= est
%! assert (all (Y(:,3) >= Y(:,1) - 1e-9));                       # upper >= est
%! assert (AUC(2) <= AUC(1) && AUC(1) <= AUC(3));
%! assert (isequal (SUBY, Y));

%!test  # bootstrap: percentile bounds also run and bracket the AUC estimate
%! rand ("seed", 7);  randn ("seed", 7);
%! [~, Y, ~, AUC] = perfcurve (labels, scores, 1, "NBoot", 200, ...
%!                             "BootType", "percentile");
%! assert (columns (Y), 3);
%! assert (AUC(2) <= AUC(1) && AUC(1) <= AUC(3));

## Test input validation
%!error <Invalid call to perfcurve> perfcurve (1, 2)
%!error <perfcurve: SCORES must be a vector of real values.> ...
%! perfcurve ([1 0], ones (2, 2), 1)
%!error <perfcurve: LABELS must have one element per score.> ...
%! perfcurve ([1 0 1], [0.5 0.4], 1)
%!error <perfcurve: LABELS must contain both positive and negative classes.> ...
%! perfcurve ([1 1 1], [0.5 0.4 0.3], 1)
%!error <perfcurve: unsupported criterion 'foo'.> ...
%! perfcurve ([1 0], [0.5 0.4], 1, "YCrit", "foo")
%!error <perfcurve: 'Cost' must be a 2-by-2 matrix.> ...
%! perfcurve ([1 0], [0.5 0.4], 1, "Cost", [0 1 0])
%!error <perfcurve: 'NBoot' must be a non-negative integer.> ...
%! perfcurve ([1 0], [0.5 0.4], 1, "NBoot", -5)
%!error <perfcurve: 'Alpha' must be a scalar in the range .0,1..> ...
%! perfcurve ([1 0], [0.5 0.4], 1, "Alpha", 1.5)
