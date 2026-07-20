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

classdef rocmetrics
  ## -*- texinfo -*-
  ## @deftp {statistics} rocmetrics
  ##
  ## Receiver operating characteristic (ROC) metrics for classifier output.
  ##
  ## The @code{rocmetrics} class evaluates a classifier's performance by
  ## computing, for each class, a one-versus-all ROC curve together with a set
  ## of threshold-dependent performance metrics.  It stores the results in the
  ## @code{Metrics} table and the per-class area under the curve in @code{AUC},
  ## and provides the @code{addMetrics}, @code{average}, and @code{plot}
  ## methods for follow-up analysis.
  ##
  ## For a problem with @math{K} classes and scores supplied as an
  ## @math{N}-by-@math{K} matrix, the discriminant score used for class
  ## @var{k} is the one-versus-all margin
  ## @code{Scores(:,k) - max (Scores(:,j))} over @code{j != k}, matching
  ## MATLAB's @code{rocmetrics}.  Every metric is evaluated at each distinct
  ## value of that margin.
  ##
  ## @seealso{perfcurve, confusionmat, confusionchart}
  ## @end deftp

  properties (GetAccess = public, SetAccess = private)

    ## -*- texinfo -*-
    ## @deftp {rocmetrics} {property} Metrics
    ##
    ## Performance metrics table
    ##
    ## Table of performance metrics, vertically concatenated across the classes
    ## in @code{ClassNames} order with one row per distinct threshold.  The
    ## standard variables are @code{ClassName}, @code{Threshold},
    ## @code{FalsePositiveRate}, and @code{TruePositiveRate}, followed by one
    ## variable for each metric requested through @qcode{'AdditionalMetrics'}.
    ## This property is read-only.
    ##
    ## @end deftp
    Metrics = [];

    ## -*- texinfo -*-
    ## @deftp {rocmetrics} {property} AUC
    ##
    ## Area under the ROC curve
    ##
    ## Row vector holding the area under the one-versus-all ROC curve for each
    ## class, in @code{ClassNames} order.  This property is read-only.
    ##
    ## @end deftp
    AUC = [];

    ## -*- texinfo -*-
    ## @deftp {rocmetrics} {property} ClassNames
    ##
    ## Class names
    ##
    ## Class names for which the ROC metrics are computed, in the column order
    ## of @code{Scores}.  This property is read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {rocmetrics} {property} Cost
    ##
    ## Misclassification cost matrix
    ##
    ## Square misclassification-cost matrix, with zero diagonal and unit
    ## off-diagonal entries by default.  It is used by the @code{ExpectedCost}
    ## metric.  This property is read-only.
    ##
    ## @end deftp
    Cost = [];

    ## -*- texinfo -*-
    ## @deftp {rocmetrics} {property} Prior
    ##
    ## Prior class probabilities
    ##
    ## Row vector of prior class probabilities, in @code{ClassNames} order,
    ## summing to one.  This property is read-only.
    ##
    ## @end deftp
    Prior = [];

    ## -*- texinfo -*-
    ## @deftp {rocmetrics} {property} Labels
    ##
    ## Observation labels
    ##
    ## True class labels supplied at construction, one per observation.  This
    ## property is read-only.
    ##
    ## @end deftp
    Labels = [];

    ## -*- texinfo -*-
    ## @deftp {rocmetrics} {property} Scores
    ##
    ## Classification scores
    ##
    ## Classification scores supplied at construction, as an
    ## @math{N}-by-@math{K} matrix.  This property is read-only.
    ##
    ## @end deftp
    Scores = [];

    ## -*- texinfo -*-
    ## @deftp {rocmetrics} {property} Weights
    ##
    ## Observation weights
    ##
    ## Non-negative observation weights, one per observation.  Defaults to a
    ## vector of ones.  This property is read-only.
    ##
    ## @end deftp
    Weights = [];
  endproperties

  properties (GetAccess = public, SetAccess = private, Hidden)
    ## Per-class raw curve data retained so that addMetrics can append columns
    ## without recomputing, and so average/plot can reuse the curves.  Each
    ## element has fields Threshold, tp, fp, tn, fn, P, N, priorPos, cNP, cPN.
    ClassData_ = struct ([]);
    ## Names of the additional metrics currently materialised in Metrics.
    AddMetrics_ = {};
  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {rocmetrics} {@var{obj} =} rocmetrics (@var{labels}, @var{scores}, @var{classnames})
    ## @deftypefnx {rocmetrics} {@var{obj} =} rocmetrics (@dots{}, @var{Name}, @var{Value})
    ##
    ## Create a @code{rocmetrics} object from labels and classification scores.
    ##
    ## @var{labels} is a vector of true class labels with one element per
    ## observation; it may be numeric, logical, a character matrix, or a cell
    ## array of character vectors.  @var{scores} is an @math{N}-by-@math{K}
    ## numeric matrix of classification scores, where @code{scores(i,k)} is the
    ## score of observation @var{i} for the class @code{classnames(k)}.
    ## @var{classnames} lists the @math{K} classes in the column order of
    ## @var{scores}.
    ##
    ## The following @qcode{Name-Value} pairs are supported:
    ##
    ## @multitable @columnfractions 0.25 0.75
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'AdditionalMetrics'} @tab A character vector or cell array
    ## of metric names to append to @code{Metrics}.  Supported names are
    ## @qcode{'TruePositives'}, @qcode{'FalseNegatives'},
    ## @qcode{'FalsePositives'}, @qcode{'TrueNegatives'},
    ## @qcode{'SumOfTrueAndFalsePositives'},
    ## @qcode{'RateOfPositivePredictions'},
    ## @qcode{'RateOfNegativePredictions'}, @qcode{'Accuracy'},
    ## @qcode{'FalseNegativeRate'}, @qcode{'TrueNegativeRate'},
    ## @qcode{'PositivePredictiveValue'}, @qcode{'NegativePredictiveValue'},
    ## @qcode{'ExpectedCost'}, and @qcode{'f1score'}.
    ##
    ## @item @qcode{'Prior'} @tab Prior class probabilities, given as
    ## @qcode{'empirical'} (default), @qcode{'uniform'}, or a numeric vector
    ## with one value per class.
    ##
    ## @item @qcode{'Cost'} @tab A @math{K}-by-@math{K} misclassification-cost
    ## matrix used by the @code{ExpectedCost} metric.  The default has zero
    ## diagonal and unit off-diagonal entries.
    ##
    ## @item @qcode{'Weights'} @tab A vector of non-negative observation
    ## weights.  The default is a vector of ones.
    ##
    ## @item @qcode{'NaNFlag'} @tab How to treat @code{NaN} scores:
    ## @qcode{'omitnan'} (default) drops the affected observations, while
    ## @qcode{'includenan'} treats them as always classified negative.
    ##
    ## @item @qcode{'FixedMetricValues'} @tab @qcode{'all'} (default) to use
    ## every distinct threshold, or a numeric vector of threshold values at
    ## which to report the curve (nearest actual thresholds are returned).
    ## @end multitable
    ##
    ## Construction from a trained model object, the @qcode{'FixedMetric'}
    ## grids other than thresholds, and bootstrap confidence intervals are not
    ## implemented.
    ##
    ## @end deftypefn
    function obj = rocmetrics (labels, scores, classnames, varargin)

      if (nargin < 3)
        error ("rocmetrics: too few input arguments.");
      endif
      if (isa (labels, 'ClassificationSVM') || isobject (labels))
        error (strcat (["rocmetrics: construction from a trained model"], ...
                       [" object is not supported; supply labels, scores,"], ...
                       [" and class names."]));
      endif

      ## Validate scores.
      if (! (isnumeric (scores) && isreal (scores) && ! isempty (scores)))
        error ("rocmetrics: SCORES must be a real numeric matrix.");
      endif
      if (isvector (scores))
        scores = scores(:);
      endif
      [n, K] = size (scores);

      ## Validate labels.
      if (rocmetrics.numelLabels_ (labels) != n)
        error ("rocmetrics: LABELS must have one element per observation.");
      endif

      ## Validate class names against the score columns.
      cn = rocmetrics.tocol_ (classnames);
      if (numel (cn) != K)
        error (strcat (["rocmetrics: CLASSNAMES must have one element per"], ...
                       [" column of SCORES."]));
      endif

      ## Defaults and Name-Value parsing.
      addmetrics = {};
      prior = 'empirical';
      cost = [];
      weights = [];
      nanflag = 'omitnan';
      fixedvals = 'all';
      if (mod (numel (varargin), 2) != 0)
        error ("rocmetrics: optional arguments must be Name-Value pairs.");
      endif
      for k = 1:2:numel (varargin)
        name = varargin{k};
        if (! (ischar (name) && isrow (name)))
          error ("rocmetrics: parameter names must be character vectors.");
        endif
        value = varargin{k+1};
        switch (lower (name))
          case 'additionalmetrics'
            addmetrics = value;
          case 'prior'
            prior = value;
          case 'cost'
            cost = value;
          case 'weights'
            weights = value;
          case 'nanflag'
            nanflag = lower (value);
          case 'fixedmetricvalues'
            fixedvals = value;
          case 'fixedmetric'
            if (! (ischar (value) && strcmpi (value, 'Thresholds')))
              error (strcat (["rocmetrics: only the 'Thresholds' value of"], ...
                             [" 'FixedMetric' is supported."]));
            endif
          case {'numbootstraps', 'bootstraptype', 'bootstrapoptions'}
            error (strcat (["rocmetrics: bootstrap confidence intervals are"], ...
                           [" not supported."]));
          otherwise
            error ("rocmetrics: unknown parameter name '%s'.", name);
        endswitch
      endfor

      addmetrics = rocmetrics.resolveMetricNames_ (addmetrics);
      if (! any (strcmpi (nanflag, {'omitnan', 'includenan'})))
        error ("rocmetrics: 'NaNFlag' must be 'omitnan' or 'includenan'.");
      endif

      ## Weights.
      if (isempty (weights))
        w = ones (n, 1);
      else
        if (! (isnumeric (weights) && isreal (weights) && isvector (weights)
               && numel (weights) == n && all (weights >= 0)))
          error (strcat (["rocmetrics: 'Weights' must be a non-negative"], ...
                         [" vector with one element per observation."]));
        endif
        w = double (weights(:));
      endif

      ## Prior probabilities (in ClassNames order).
      priorvec = rocmetrics.computePrior_ (prior, labels, cn, w, K);

      ## Misclassification-cost matrix.
      cost = rocmetrics.validateCost_ (cost, K);

      ## Build the per-class curves and assemble the Metrics table.
      obj.Labels = labels;
      obj.Scores = scores;
      obj.Weights = w;
      obj.ClassNames = classnames;
      obj.Prior = priorvec;
      obj.Cost = cost;
      obj.AddMetrics_ = addmetrics;
      obj = obj.buildCurves_ (cn, nanflag, fixedvals);
      obj = obj.assembleTable_ ();

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {rocmetrics} {@var{obj} =} addMetrics (@var{obj}, @var{metrics})
    ##
    ## Append additional performance metrics to an existing @code{rocmetrics}
    ## object.
    ##
    ## @var{metrics} is a metric name or a cell array of metric names, chosen
    ## from the list supported by the @qcode{'AdditionalMetrics'} constructor
    ## argument.  The named metrics are appended as new variables of the
    ## @code{Metrics} table; metrics already present are left unchanged.
    ##
    ## @end deftypefn
    function obj = addMetrics (obj, metrics)
      if (nargin != 2)
        error ("rocmetrics.addMetrics: METRICS input is required.");
      endif
      newmetrics = rocmetrics.resolveMetricNames_ (metrics);
      for k = 1:numel (newmetrics)
        if (! any (strcmp (newmetrics{k}, obj.AddMetrics_)))
          obj.AddMetrics_{end+1} = newmetrics{k};
        endif
      endfor
      obj = obj.assembleTable_ ();
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {rocmetrics} {[@var{FPR}, @var{TPR}, @var{Thresholds}, @var{AUC}] =} average (@var{obj}, @var{type})
    ##
    ## Compute an averaged ROC curve across the classes of a @code{rocmetrics}
    ## object.
    ##
    ## @var{type} selects the averaging method: @qcode{'macro'} (unweighted
    ## mean of the per-class curves), @qcode{'micro'} (a single curve pooling
    ## every one-versus-all instance), or @qcode{'weighted'} (mean of the
    ## per-class curves weighted by @code{Prior}).  The function returns the
    ## averaged false and true positive rates @var{FPR} and @var{TPR}, the
    ## corresponding @var{Thresholds}, and the area @var{AUC} under the
    ## averaged curve.
    ##
    ## The averaged curve is evaluated on the union of the per-class
    ## thresholds.  MATLAB inserts additional staircase points when building
    ## the averaged curve, so the exact rows and the averaged @var{AUC} may
    ## differ slightly from MATLAB.
    ##
    ## @end deftypefn
    function [FPR, TPR, Thresholds, AUC] = average (obj, type)
      if (nargin < 2)
        type = 'macro';
      endif
      [FPR, TPR, Thresholds, AUC] = obj.averageCurve_ (lower (type));
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {rocmetrics} {} plot (@var{obj})
    ## @deftypefnx {rocmetrics} {@var{h} =} plot (@var{obj})
    ##
    ## Plot the per-class ROC curves of a @code{rocmetrics} object.
    ##
    ## Each class in @code{ClassNames} contributes one true-positive-rate
    ## versus false-positive-rate curve.  A handle to the line objects is
    ## returned in @var{h} when requested.
    ##
    ## @end deftypefn
    function varargout = plot (obj)
      cn = rocmetrics.tocol_ (obj.ClassNames);
      h = [];
      newplot ();
      hold_state = ishold ();
      hold ("on");
      leg = cell (1, numel (obj.ClassData_));
      for k = 1:numel (obj.ClassData_)
        cd = obj.ClassData_(k);
        h(end+1) = plot (cd.fp / cd.N, cd.tp / cd.P, '-o');
        leg{k} = sprintf ("%s (AUC = %.4f)", ...
                          rocmetrics.nameStr_ (cn, k), obj.AUC(k));
      endfor
      plot ([0, 1], [0, 1], 'k:');
      if (! hold_state)
        hold ("off");
      endif
      xlabel ("False Positive Rate");
      ylabel ("True Positive Rate");
      legend (leg, "location", "southeast");
      if (nargout > 0)
        varargout{1} = h;
      endif
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {rocmetrics} {} disp (@var{obj})
    ##
    ## Display a summary of a @code{rocmetrics} object.
    ##
    ## @end deftypefn
    function disp (obj)
      printf ("  rocmetrics with properties:\n\n");
      printf ("    ClassNames: %d class(es)\n", numel (obj.ClassData_));
      printf ("           AUC: [%s]\n", ...
              strtrim (sprintf ("%.4f ", obj.AUC)));
      printf ("       Metrics: [%dx%d table]\n\n", size (obj.Metrics, 1), ...
              size (obj.Metrics, 2));
    endfunction

  endmethods

  methods (Access = private)

    ## Build the one-versus-all curve data for every class.
    function obj = buildCurves_ (obj, cn, nanflag, fixedvals)
      scores = obj.Scores;
      [n, K] = size (scores);
      w = obj.Weights;
      cd = struct ("Threshold", {}, "tp", {}, "fp", {}, "tn", {}, ...
                   "fn", {}, "P", {}, "N", {}, "priorPos", {}, ...
                   "cNP", {}, "cPN", {});
      auc = zeros (1, K);
      for k = 1:K
        ## One-versus-all margin score.
        if (K == 1)
          margin = scores(:,1);
        else
          other = max (scores(:,[1:k-1, k+1:K]), [], 2);
          margin = scores(:,k) - other;
        endif
        ## Collapse floating-point noise introduced by the margin so that
        ## thresholds that are equal in exact arithmetic group together.
        margin = rocmetrics.roundSig_ (margin, 14);
        ispos = rocmetrics.labelEq_ (obj.Labels, cn, k);
        isneg = ! ispos;
        wk = w;

        ## NaN score handling.
        isnanm = isnan (margin);
        if (any (isnanm))
          if (strcmp (nanflag, "omitnan"))
            keep = ! isnanm;
            margin = margin(keep);  ispos = ispos(keep);
            isneg = isneg(keep);    wk = wk(keep);
          else
            margin(isnanm) = -Inf;
          endif
        endif

        if (! any (ispos) || ! any (isneg))
          error (strcat (["rocmetrics: class '%s' must have both positive"], ...
                         [" and negative observations."]), ...
                 rocmetrics.nameStr_ (cn, k));
        endif

        ## Prior-adjusted effective weights (verified MATLAB normalisation).
        pp = obj.Prior(k);
        pn = 1 - pp;
        Wpos = sum (wk(ispos));
        Wneg = sum (wk(isneg));
        S = 1 / (pp / Wpos + pn / Wneg);
        we = zeros (size (wk));
        we(ispos) = wk(ispos) * (pp / Wpos) * S;
        we(isneg) = wk(isneg) * (pn / Wneg) * S;

        ## Cumulative counts at each distinct threshold, descending.
        [ss, ord] = sort (margin, "descend");
        pw = cumsum (we(ord) .* ispos(ord));
        nw = cumsum (we(ord) .* isneg(ord));
        chg = find (diff (ss) != 0);
        last = [chg; numel(ss)];
        tp = [0; pw(last)];
        fp = [0; nw(last)];
        T = [ss(1); ss(last)];
        P = pp * S;
        N = pn * S;

        ## Optionally restrict to requested threshold values.
        if (! (ischar (fixedvals) && strcmpi (fixedvals, "all")))
          idx = rocmetrics.selectThresholds_ (T, fixedvals(:));
          tp = tp(idx);  fp = fp(idx);  T = T(idx);
        endif

        [cNP, cPN] = rocmetrics.ovaCost_ (obj.Cost, K, k);
        cd(k).Threshold = T;
        cd(k).tp = tp;  cd(k).fp = fp;
        cd(k).tn = N - fp;  cd(k).fn = P - tp;
        cd(k).P = P;  cd(k).N = N;
        cd(k).priorPos = pp;  cd(k).cNP = cNP;  cd(k).cPN = cPN;
        auc(k) = trapz (fp / N, tp / P);
      endfor
      obj.ClassData_ = cd;
      obj.AUC = auc;
    endfunction

    ## Assemble the Metrics table from the per-class curve data.
    function obj = assembleTable_ (obj)
      cn = rocmetrics.tocol_ (obj.ClassNames);
      names = {'ClassName', 'Threshold', 'FalsePositiveRate', ...
               'TruePositiveRate', obj.AddMetrics_{:}};
      ncol = numel (names);
      cols = cell (1, ncol);
      for c = 1:ncol
        cols{c} = [];
      endfor
      classcol = {};
      classnum = [];
      numericnames = isnumeric (cn) || islogical (cn);
      for k = 1:numel (obj.ClassData_)
        cd = obj.ClassData_(k);
        m = numel (cd.Threshold);
        if (numericnames)
          classnum = [classnum; repmat(double (cn(k)), m, 1)];
        else
          classcol = [classcol; repmat({rocmetrics.nameStr_(cn, k)}, m, 1)];
        endif
        cols{2} = [cols{2}; cd.Threshold];
        cols{3} = [cols{3}; cd.fp / cd.N];
        cols{4} = [cols{4}; cd.tp / cd.P];
        for c = 5:ncol
          cols{c} = [cols{c}; rocmetrics.metricValue_(names{c}, cd)];
        endfor
      endfor
      if (numericnames)
        cols{1} = classnum;
      else
        cols{1} = classcol;
      endif
      obj.Metrics = table (cols{:}, 'VariableNames', names);
    endfunction

    ## Averaged ROC curve for a given averaging type.
    function [FPR, TPR, Thresholds, AUC] = averageCurve_ (obj, type)
      cd = obj.ClassData_;
      K = numel (cd);
      if (strcmp (type, "micro"))
        ## Pool every one-versus-all instance into one binary problem.
        allthr = [];  alltp = [];  allfp = [];  P = 0;  N = 0;
        for k = 1:K
          allthr = [allthr; cd(k).Threshold];
          P = P + cd(k).P;  N = N + cd(k).N;
        endfor
        Thresholds = unique (allthr, "sorted");
        Thresholds = flipud (Thresholds(:));
        tp = zeros (size (Thresholds));
        fp = zeros (size (Thresholds));
        for k = 1:K
          tp = tp + rocmetrics.stepAt_ (cd(k).Threshold, cd(k).tp, Thresholds);
          fp = fp + rocmetrics.stepAt_ (cd(k).Threshold, cd(k).fp, Thresholds);
        endfor
        FPR = fp / N;  TPR = tp / P;
      else
        ## Macro / weighted: mean of the per-class rates on the union grid.
        if (strcmp (type, "weighted"))
          wt = obj.Prior(:)';
        elseif (strcmp (type, "macro"))
          wt = ones (1, K) / K;
        else
          error ("rocmetrics.average: unknown averaging type '%s'.", type);
        endif
        wt = wt / sum (wt);
        allthr = [];
        for k = 1:K
          allthr = [allthr; cd(k).Threshold];
        endfor
        Thresholds = unique (allthr, "sorted");
        Thresholds = flipud (Thresholds(:));
        FPR = zeros (size (Thresholds));
        TPR = zeros (size (Thresholds));
        for k = 1:K
          fk = rocmetrics.stepAt_ (cd(k).Threshold, cd(k).fp / cd(k).N, ...
                                   Thresholds);
          tk = rocmetrics.stepAt_ (cd(k).Threshold, cd(k).tp / cd(k).P, ...
                                   Thresholds);
          FPR = FPR + wt(k) * fk;
          TPR = TPR + wt(k) * tk;
        endfor
      endif
      AUC = trapz (FPR, TPR);
    endfunction

  endmethods

  methods (Static, Access = private)

    ## Number of labels in a vector or character matrix.
    function m = numelLabels_ (labels)
      if (ischar (labels))
        m = rows (labels);
      else
        m = numel (labels);
      endif
    endfunction

    ## Column form of a class-name list.
    function c = tocol_ (classnames)
      if (ischar (classnames))
        c = cellstr (classnames);
      elseif (iscell (classnames))
        c = classnames(:);
      else
        c = classnames(:);
      endif
    endfunction

    ## Logical membership of the labels in the k-th class.
    function tf = labelEq_ (labels, cn, k)
      if (isnumeric (cn) || islogical (cn))
        tf = (labels(:) == cn(k));
      else
        if (ischar (labels))
          labels = cellstr (labels);
        endif
        tf = strcmp (labels(:), cn{k});
      endif
    endfunction

    ## Printable name of the k-th class.
    function s = nameStr_ (cn, k)
      if (iscell (cn))
        s = cn{k};
      elseif (ischar (cn))
        s = cn(k,:);
      else
        s = num2str (cn(k));
      endif
    endfunction

    ## Prior probabilities in ClassNames order, summing to one.
    function p = computePrior_ (prior, labels, cn, w, K)
      if (ischar (prior))
        switch (lower (prior))
          case 'empirical'
            p = zeros (1, K);
            for k = 1:K
              p(k) = sum (w(rocmetrics.labelEq_ (labels, cn, k)));
            endfor
            if (sum (p) == 0)
              error ("rocmetrics: empirical prior is undefined (no labels).");
            endif
            p = p / sum (p);
          case 'uniform'
            p = ones (1, K) / K;
          otherwise
            error ("rocmetrics: 'Prior' must be 'empirical', 'uniform', or a vector.");
        endswitch
      elseif (isnumeric (prior) && isvector (prior) && numel (prior) == K
              && all (prior >= 0) && sum (prior) > 0)
        p = prior(:)' / sum (prior);
      else
        error ("rocmetrics: 'Prior' must be 'empirical', 'uniform', or a vector.");
      endif
    endfunction

    ## Validate or default the misclassification-cost matrix.
    function cost = validateCost_ (cost, K)
      if (isempty (cost))
        cost = ones (K) - eye (K);
      elseif (! (isnumeric (cost) && isequal (size (cost), [K, K])))
        error ("rocmetrics: 'Cost' must be a %d-by-%d matrix.", K, K);
      endif
    endfunction

    ## One-versus-all costs C(N|P) and C(P|N) for the k-th class.
    function [cNP, cPN] = ovaCost_ (cost, K, k)
      if (K == 2)
        other = 3 - k;
        cNP = cost(k, other);
        cPN = cost(other, k);
      else
        ## For K != 2 only the default identity cost is validated; the
        ## off-diagonal cost of the k-th class is used symmetrically.
        off = cost(k, [1:k-1, k+1:K]);
        cNP = off(1);
        cPN = off(1);
      endif
    endfunction

    ## Value of a named metric for a class-data record.
    function v = metricValue_ (name, cd)
      tp = cd.tp;  fp = cd.fp;  tn = cd.tn;  fn = cd.fn;
      P = cd.P;  N = cd.N;  total = P + N;
      pp = cd.priorPos;  pn = 1 - pp;
      switch (name)
        case 'TruePositives'
          v = tp;
        case 'FalseNegatives'
          v = fn;
        case 'FalsePositives'
          v = fp;
        case 'TrueNegatives'
          v = tn;
        case 'SumOfTrueAndFalsePositives'
          v = tp + fp;
        case 'RateOfPositivePredictions'
          v = (tp + fp) / total;
        case 'RateOfNegativePredictions'
          v = (tn + fn) / total;
        case 'Accuracy'
          v = (tp + tn) / total;
        case 'FalseNegativeRate'
          v = fn / P;
        case 'TrueNegativeRate'
          v = tn / N;
        case 'PositivePredictiveValue'
          v = tp ./ (tp + fp);
        case 'NegativePredictiveValue'
          v = tn ./ (tn + fn);
        case 'ExpectedCost'
          v = pp * pn * (pp * (fn / P) * cd.cNP + pn * (fp / N) * cd.cPN);
        case 'f1score'
          v = 2 * tp ./ (2 * tp + fp + fn);
        otherwise
          error ("rocmetrics: unsupported metric '%s'.", name);
      endswitch
    endfunction

    ## Canonicalise a metric name or list of metric names.
    function names = resolveMetricNames_ (metrics)
      if (isempty (metrics))
        names = {};
        return;
      endif
      if (ischar (metrics))
        metrics = {metrics};
      elseif (! iscell (metrics))
        error ("rocmetrics: metric names must be a string or cell array.");
      endif
      canon = {'TruePositives', 'FalseNegatives', 'FalsePositives', ...
               'TrueNegatives', 'SumOfTrueAndFalsePositives', ...
               'RateOfPositivePredictions', 'RateOfNegativePredictions', ...
               'Accuracy', 'FalseNegativeRate', 'TrueNegativeRate', ...
               'PositivePredictiveValue', 'NegativePredictiveValue', ...
               'ExpectedCost', 'f1score'};
      names = cell (1, numel (metrics));
      for k = 1:numel (metrics)
        idx = find (strcmpi (metrics{k}, canon));
        if (isempty (idx))
          error ("rocmetrics: unrecognised metric '%s'.", metrics{k});
        endif
        names{k} = canon{idx};
      endfor
    endfunction

    ## Indices of the thresholds nearest the requested values.
    function idx = selectThresholds_ (T, vals)
      idx = zeros (numel (vals), 1);
      for k = 1:numel (vals)
        if (vals(k) >= T(1))
          idx(k) = 1;
        else
          [~, idx(k)] = min (abs (T - vals(k)));
        endif
      endfor
    endfunction

    ## Round to SIG significant digits, collapsing floating-point noise while
    ## preserving genuinely distinct values.
    function y = roundSig_ (x, sig)
      y = x;
      nz = (x != 0) & isfinite (x);
      e = floor (log10 (abs (x(nz))));
      f = 10 .^ (sig - 1 - e);
      y(nz) = round (x(nz) .* f) ./ f;
    endfunction

    ## Step-function value of a curve (thresholds T descending, values Y) at
    ## the query thresholds TQ, using the operating point for score >= TQ.
    function yq = stepAt_ (T, Y, TQ)
      yq = zeros (size (TQ));
      for i = 1:numel (TQ)
        j = find (T >= TQ(i), 1, "last");
        if (isempty (j))
          yq(i) = Y(1);
        else
          yq(i) = Y(j);
        endif
      endfor
    endfunction

  endmethods

endclassdef

%!demo
%! ## One-versus-all ROC curves for a three-class problem
%! labels = [1 1 2 2 3 3]';
%! scores = [0.9 0.05 0.05; 0.6 0.3 0.1; 0.2 0.7 0.1; ...
%!           0.1 0.6 0.3; 0.2 0.2 0.6; 0.1 0.3 0.6];
%! rocObj = rocmetrics (labels, scores, [1 2 3]);
%! disp (rocObj.AUC)
%! plot (rocObj);

%!demo
%! ## Binary ROC metrics with additional performance metrics
%! labels = [1 1 1 1 0 0]';
%! p = [0.9 0.7 0.4 0.3 0.8 0.2]';
%! scores = [1-p p];
%! rocObj = rocmetrics (labels, scores, [0 1], ...
%!                      "AdditionalMetrics", {"Accuracy", "ExpectedCost"});
%! head = rocObj.Metrics(1:4,:);
%! disp (head);

%!shared labels, scores, cn
%! labels = [1 1 2 2 3 3]';
%! scores = [0.9 0.05 0.05; 0.6 0.3 0.1; 0.2 0.7 0.1; ...
%!           0.1 0.6 0.3; 0.2 0.2 0.6; 0.1 0.3 0.6];
%! cn = [1 2 3];

%!test  # MATLAB parity: default Metrics columns, thresholds, FPR/TPR, AUC
%! r = rocmetrics (labels, scores, cn);
%! assert_equal (r.Metrics.Properties.VariableNames, ...
%!               {'ClassName', 'Threshold', 'FalsePositiveRate', 'TruePositiveRate'});
%! assert_equal (double (r.Metrics.ClassName)', ...
%!               [1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 3], 0);
%! assert_equal (r.Metrics.Threshold, ...
%!               [0.85;0.85;0.3;-0.4;-0.5;0.5;0.5;0.3;-0.3;-0.4;-0.85; ...
%!                0.4;0.4;0.3;-0.3;-0.5;-0.6;-0.85], 1e-12);
%! assert_equal (r.Metrics.TruePositiveRate, ...
%!               [0;0.5;1;1;1;0;0.5;1;1;1;1;0;0.5;1;1;1;1;1], 1e-12);
%! assert_equal (r.AUC, [1 1 1], 1e-12);
%! assert_equal (r.Prior, [1 1 1]/3, 1e-12);

%!test  # MATLAB parity: one-versus-all margin score sets the thresholds
%! r = rocmetrics (labels, scores, cn);
%! m1 = r.Metrics(double (r.Metrics.ClassName) == 1, :);
%! assert_equal (m1.Threshold, [0.85;0.85;0.3;-0.4;-0.5], 1e-12);
%! assert_equal (m1.FalsePositiveRate, [0;0;0;0.25;1], 1e-12);

%!test  # MATLAB parity: additional count and rate metrics (imbalanced binary)
%! y = [1 1 1 1 0 0]';
%! p = [0.9 0.7 0.4 0.3 0.8 0.2]';
%! s = [1-p p];
%! mets = {"TruePositives","FalsePositives","TrueNegatives","FalseNegatives"};
%! r = rocmetrics (y, s, [0 1], "AdditionalMetrics", mets);
%! m0 = r.Metrics(double (r.Metrics.ClassName) == 0, :);
%! assert_equal (m0.TruePositives, [0;0.5;0.5;0.5;0.5;1;1], 1e-12);
%! assert_equal (m0.FalsePositives, [0;0;0.5;1;1.5;1.5;2], 1e-12);
%! assert_equal (m0.TrueNegatives, [2;2;1.5;1;0.5;0.5;0], 1e-12);
%! assert_equal (m0.FalseNegatives, [1;0.5;0.5;0.5;0.5;0;0], 1e-12);
%! assert_equal (r.AUC, [0.625 0.625], 1e-12);

%!test  # MATLAB parity: ExpectedCost with default and custom cost
%! y = [1 1 1 1 0 0]';
%! p = [0.9 0.7 0.4 0.3 0.8 0.2]';
%! s = [1-p p];
%! r = rocmetrics (y, s, [0 1], "AdditionalMetrics", "ExpectedCost");
%! ec0 = r.Metrics.ExpectedCost(double (r.Metrics.ClassName) == 0);
%! assert_equal (ec0(1), 2/27, 1e-12);
%! assert_equal (ec0(5), 4/27, 1e-12);
%! rc = rocmetrics (y, s, [0 1], "Cost", [0 2; 1 0], ...
%!                  "AdditionalMetrics", "ExpectedCost");
%! ecc = rc.Metrics.ExpectedCost(double (rc.Metrics.ClassName) == 0);
%! assert_equal (ecc(1), 4/27, 1e-12);
%! assert_equal (ecc(5), 5/27, 1e-12);

%!test  # MATLAB parity: uniform prior reweights the counts
%! y = [1 1 1 1 0 0]';
%! p = [0.9 0.7 0.4 0.3 0.8 0.2]';
%! s = [1-p p];
%! r = rocmetrics (y, s, [0 1], "Prior", "uniform", ...
%!                 "AdditionalMetrics", "TruePositives");
%! assert_equal (r.Prior, [0.5 0.5], 1e-12);
%! tp0 = r.Metrics.TruePositives(double (r.Metrics.ClassName) == 0);
%! assert_equal (max (tp0), 4/3, 1e-12);

%!test  # explicit weights renormalise the effective counts
%! y = [1 1 1 1 0 0]';
%! p = [0.9 0.7 0.4 0.3 0.8 0.2]';
%! s = [1-p p];
%! w = [1 1 1 1 2 2]';
%! r = rocmetrics (y, s, [0 1], "Weights", w, ...
%!                 "AdditionalMetrics", "TruePositives");
%! assert_equal (r.Prior, [0.5 0.5], 1e-12);
%! tp0 = r.Metrics.TruePositives(double (r.Metrics.ClassName) == 0);
%! assert_equal (max (tp0), 2, 1e-12);

%!test  # addMetrics appends without recomputing the curve
%! r = rocmetrics (labels, scores, cn);
%! r = addMetrics (r, "Accuracy");
%! assert_equal (any (strcmp ("Accuracy", ...
%!               r.Metrics.Properties.VariableNames)), true);
%! acc1 = r.Metrics.Accuracy(1);
%! assert_equal (acc1, 2/3, 1e-12);

%!test  # PPV is NaN at the origin, NPV is NaN at the all-positive end
%! r = rocmetrics (labels, scores, cn, "AdditionalMetrics", ...
%!                 {"PositivePredictiveValue","NegativePredictiveValue"});
%! m1 = r.Metrics(double (r.Metrics.ClassName) == 1, :);
%! assert_equal (isnan (m1.PositivePredictiveValue(1)), true);
%! assert_equal (isnan (m1.NegativePredictiveValue(end)), true);

%!test  # cell-array string labels and names
%! y = {"a","a","b","b","c","c"};
%! r = rocmetrics (y, scores, {"a","b","c"});
%! assert_equal (r.AUC, [1 1 1], 1e-12);
%! assert_equal (iscellstr (r.Metrics.ClassName), true);

%!test  # FixedMetricValues selects nearest thresholds
%! r = rocmetrics (labels, scores, cn, "FixedMetricValues", [Inf 0.4 -0.5]);
%! m1 = r.Metrics(double (r.Metrics.ClassName) == 1, :);
%! assert_equal (numel (m1.Threshold), 3);
%! assert_equal (m1.Threshold(1), 0.85, 1e-12);

%!test  # average returns a valid curve for each averaging type
%! r = rocmetrics (labels, scores, cn);
%! [fpr, tpr, thr, auc] = average (r, "macro");
%! assert_equal (fpr(1), 0, 1e-12);
%! assert_equal (tpr(end), 1, 1e-12);
%! assert_equal (auc >= 0 && auc <= 1, true);
%! [~, ~, ~, aucmi] = average (r, "micro");
%! assert_equal (aucmi >= 0 && aucmi <= 1, true);

## Test input validation
%!error <rocmetrics: too few input arguments.> rocmetrics ([1 0], [0.4 0.6])
%!error <rocmetrics: SCORES must be a real numeric matrix.> ...
%! rocmetrics ([1 0], {1, 2}, [0 1])
%!error <rocmetrics: LABELS must have one element per observation.> ...
%! rocmetrics ([1 0 1], [0.4 0.6; 0.5 0.5], [0 1])
%!error <rocmetrics: CLASSNAMES must have one element per column of SCORES.> ...
%! rocmetrics ([1 0], [0.4 0.6; 0.5 0.5], [0 1 2])
%!error <rocmetrics: unrecognised metric 'foo'.> ...
%! rocmetrics ([1 0], [0.4 0.6; 0.5 0.5], [0 1], "AdditionalMetrics", "foo")
%!error <rocmetrics: 'Prior' must be 'empirical', 'uniform', or a vector.> ...
%! rocmetrics ([1 0], [0.4 0.6; 0.5 0.5], [0 1], "Prior", "bogus")
%!error <rocmetrics: bootstrap confidence intervals are not supported.> ...
%! rocmetrics ([1 0], [0.4 0.6; 0.5 0.5], [0 1], "NumBootstraps", 100)
%!error <rocmetrics: unknown parameter name 'Zzz'.> ...
%! rocmetrics ([1 0], [0.4 0.6; 0.5 0.5], [0 1], "Zzz", 1)
