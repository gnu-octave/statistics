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

classdef shapley
  ## -*- texinfo -*-
  ## @deftp {statistics} shapley
  ##
  ## Shapley values for a fitted model at one or more query points.
  ##
  ## A Shapley value says how much one predictor contributed to the deviation
  ## of a prediction from the average prediction.  The values of a query point
  ## sum to that deviation exactly, which is the property the computation is
  ## built to keep.
  ##
  ## @code{@var{explainer} = shapley (@var{Mdl})} builds an explainer for the
  ## fitted model @var{Mdl} over the observations it was fitted on.  A compact
  ## model keeps none, so it must be given them as @var{X}, and so must a
  ## function handle.  Nothing is computed until query points are given, either
  ## to the constructor as @qcode{'QueryPoints'} or afterwards to
  ## @code{fit}.
  ##
  ## @code{@var{explainer} = shapley (@var{Mdl}, @var{X})} takes the
  ## observations to average over as @var{X}, a real numeric matrix of one
  ## column per predictor.
  ##
  ## @var{X} may also be a table, and so may @qcode{'QueryPoints'} and what
  ## @code{fit} is given.  Where the model names its predictors the table is
  ## read by those names and not by the order its columns come in, so a
  ## column the model was not fitted on is passed over and a value holding a
  ## level is coded as that level was coded at fitting.  A function handle
  ## names nothing, so a table given for one names the predictors itself and
  ## its columns are taken in the order they come.
  ##
  ## @code{@var{explainer} = shapley (@var{fun}, @var{X})} takes a function
  ## handle in place of a model.  @var{fun} is called with a matrix of
  ## observations and answers with one real numeric column holding one value
  ## for each, so an explainer built on a handle always has a single column of
  ## values.
  ##
  ## @multitable @columnfractions 0.28 0.02 0.7
  ## @headitem @var{Name} @tab @tab @var{Value}
  ##
  ## @item @qcode{'QueryPoints'} @tab @tab The observations to explain, one per
  ## row, with one column per predictor.  The default is none, which leaves the
  ## values unfitted.
  ##
  ## @item @qcode{'NumObservationsToSample'} @tab @tab How many observations to
  ## draw, without replacement, from those averaged over, or @qcode{'all'} for
  ## every one of them.  The default is 100, and so is any number reaching or
  ## exceeding how many there are.  A drawn sample makes the values differ from
  ## one call to the next; @qcode{'all'} is what makes them reproducible.
  ##
  ## @item @qcode{'CategoricalPredictors'} @tab @tab The predictors whose
  ## values are levels, taken as by every learner of this package.  It applies
  ## only to a function handle, a model being asked for its own.
  ##
  ## @item @qcode{'MaxNumSubsets'} @tab @tab How many predictor subsets at
  ## most to compute over, an integer above 1.  The default is the lesser of
  ## @math{2^M}, which is every subset of the @math{M} predictors, and 1024.
  ## Every subset gives the values exactly; fewer estimates them, and fewer
  ## than @math{2M+2} estimates them poorly enough to warn about.  Giving it
  ## at all asks for the subsets, so a linear model or a decision tree that
  ## would otherwise be answered from its own structure is answered over
  ## them instead.
  ##
  ## @item @qcode{'Method'} @tab @tab The algorithm, @qcode{'interventional'}
  ## by default, which averages over the observations as they stand.
  ## @qcode{'conditional'} averages instead over the tenth of them lying
  ## nearest the query point in the predictors being held, which stands in
  ## for conditioning on those predictors.  It asks more of the data and is
  ## the dearer of the two.
  ## @end multitable
  ##
  ## @code{'UseParallel'} is not implemented and is refused rather than
  ## ignored.
  ##
  ## @seealso{partialDependence, plotPartialDependence, PredictiveModel}
  ## @end deftp

  properties (SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} BlackboxModel
    ##
    ## The model being explained
    ##
    ## The fitted model, or the function handle, the explainer was built on.
    ## This property is read-only.
    ##
    ## @end deftp
    BlackboxModel = [];

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} X
    ##
    ## The observations averaged over
    ##
    ## A real numeric matrix of one row per observation and one column per
    ## predictor, either given outright or taken from the model.  It is the
    ## whole of what was given, before any sampling.  This property is
    ## read-only.
    ##
    ## Where the observations were given as a table, they are the coded
    ## matrix and not the table: a variable holding levels is stored as its
    ## level codes, and the coding is kept with the object.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} QueryPoints
    ##
    ## The observations explained
    ##
    ## A real numeric matrix of one row per query point, empty until query
    ## points are given.  This property is read-only.
    ##
    ## @end deftp
    QueryPoints = [];

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} BlackboxFitted
    ##
    ## What the model answers at the query points
    ##
    ## The response of a regression model or the predicted label of a
    ## classifier, one for each query point, empty until query points are
    ## given.  A label keeps the type of the response the model was fitted
    ## with, as everywhere in this package.  This property is read-only.
    ##
    ## @end deftp
    BlackboxFitted = [];

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} Shapley
    ##
    ## The Shapley values
    ##
    ## A table of one row per predictor, holding the predictor names in
    ## @qcode{Predictor} and the values in @qcode{Value} for a regression
    ## model or a function handle, and in one variable per class, named after
    ## it, for a classifier.  Each such variable holds one column per query
    ## point.  It is empty until query points are given.  This property is
    ## read-only.
    ##
    ## @end deftp
    Shapley = [];

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} MeanAbsoluteShapley
    ##
    ## The mean absolute Shapley value of each predictor
    ##
    ## A table laid out as @qcode{Shapley}, holding the mean over the query
    ## points of the absolute values.  With one query point it is their
    ## absolute value.  This property is read-only.
    ##
    ## @end deftp
    MeanAbsoluteShapley = [];

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} Intercept
    ##
    ## The average prediction
    ##
    ## The mean of what the model answers over the observations averaged over,
    ## a scalar for a regression model or a function handle and one value per
    ## class for a classifier.  The values of a query point sum to the
    ## deviation of its prediction from this.  This property is read-only.
    ##
    ## @end deftp
    Intercept = [];

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} Method
    ##
    ## The algorithm the values were computed with
    ##
    ## A character vector.  @qcode{'interventional-linear'} where the model
    ## predicts a weighted sum of its predictors, which is answered from the
    ## weights alone; @qcode{'interventional-tree'} for a decision tree and
    ## for an ensemble of them, which is answered leaf by leaf;
    ## @qcode{'interventional-kernel'} for
    ## every other model, which enumerates every subset of the predictors
    ## where the budget allows it and estimates the values by weighted least
    ## squares where it does not; and @qcode{'conditional-kernel'} where
    ## @qcode{'Method'} asked for conditioning.  A tree and a linear model
    ## are answered exactly however many predictors they have, where the
    ## budget stops the subsets at 1024.  This property is read-only.
    ##
    ## @end deftp
    Method = '';

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} NumSubsets
    ##
    ## How many predictor subsets the values were computed over
    ##
    ## The lesser of what @qcode{'MaxNumSubsets'} allowed and two raised to
    ## the number of predictors, which is every subset.  This property is
    ## read-only.
    ##
    ## @end deftp
    NumSubsets = [];

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} CategoricalPredictors
    ##
    ## The categorical predictors
    ##
    ## The indices of the predictors whose values are levels, empty where
    ## there are none.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} SampledObservationIndices
    ##
    ## The rows of @qcode{X} that were averaged over
    ##
    ## A sorted column of indices into @qcode{X}, holding every row where
    ## @qcode{'NumObservationsToSample'} did not draw a sample.  This property
    ## is read-only.
    ##
    ## @end deftp
    SampledObservationIndices = [];

  endproperties

  properties (Access = protected, Hidden)
    IsClass = false;        # whether the model answers with class scores
    ClassNames = [];        # the classes, in the model's own order
    PredictorNames = {};    # one name per column of X
  endproperties

  methods (Hidden)

    function disp (this)

      if (isempty (this.Shapley))
        printf ("  shapley explainer with no Shapley values fitted.\n");
        printf ("  Use the fit method to fit Shapley values.\n");
      else
        printf ("  shapley explainer with the following ");
        printf ("local Shapley values:\n\n");
        disp (this.Shapley);
      endif

    endfunction

    function display (this)

      inName = inputname (1);
      if (! isempty (inName))
        printf ("%s =\n\n", inName);
      endif
      disp (this);

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {shapley} {@var{obj} =} shapley (@var{Mdl})
    ## @deftypefnx {shapley} {@var{obj} =} shapley (@var{Mdl}, @var{X})
    ## @deftypefnx {shapley} {@var{obj} =} shapley (@var{fun}, @var{X})
    ## @deftypefnx {shapley} {@var{obj} =} shapley (@dots{}, @var{name}, @var{value})
    ##
    ## Build a Shapley explainer.
    ##
    ## The arguments are those described for the class.  Where
    ## @qcode{'QueryPoints'} is given the values are computed at once,
    ## otherwise they are left to @code{fit}.
    ##
    ## @seealso{shapley, shapley.fit}
    ## @end deftypefn

    function this = shapley (blackbox, varargin)

      ## Input validation
      if (nargin < 1)
        error ("shapley: too few input arguments.");
      endif
      if (! (isa (blackbox, 'PredictiveModel')
             || is_function_handle (blackbox)))
        error (strcat ("shapley: BLACKBOX must be a fitted model that", ...
                       " predicts, or a function handle."));
      endif

      ## The observations come before the options, as MATLAB orders them
      args = varargin;
      Data = [];
      if (! isempty (args) && ! (ischar (args{1}) || isa (args{1}, 'string')))
        Data = args{1};
        args(1) = [];
      endif

      optNames = {'QueryPoints', 'NumObservationsToSample', ...
                  'CategoricalPredictors', 'Method', 'MaxNumSubsets', ...
                  'UseParallel'};
      dfValues = {[], [], [], [], [], []};
      [QP, NumObs, CatPred, Method, MaxSub, Par, rem] = ...
                          parsePairedArguments (optNames, dfValues, args);
      if (! isempty (rem))
        error (strcat ("shapley: unknown optional argument or", ...
                       " misplaced value."));
      endif
      if (! isempty (Par))
        error ("shapley: 'UseParallel' is not implemented.");
      endif

      [base, errmsg] = shapMethod (Method);
      if (! isempty (errmsg))
        error ("shapley: %s", errmsg);
      endif

      [F, errmsg] = shapFrame (blackbox, Data, CatPred, NumObs, MaxSub);
      if (! isempty (errmsg))
        error ("shapley: %s", errmsg);
      endif
      method = shapAlgorithm (base, F, MaxSub);

      this.BlackboxModel = blackbox;
      this.X = F.X;
      this.CategoricalPredictors = F.Cat;
      this.SampledObservationIndices = F.Idx;
      this.Method = method;
      this.NumSubsets = F.NumSubsets;
      this.IsClass = F.IsClass;
      this.ClassNames = F.ClassNames;
      this.PredictorNames = F.PredictorNames;

      if (! isempty (QP))
        this = fit (this, QP);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {shapley} {@var{obj} =} fit (@var{obj}, @var{QueryPoints})
    ##
    ## Compute the Shapley values at the given query points.
    ##
    ## @var{QueryPoints} is a real numeric matrix of one row per query point
    ## and one column per predictor, or a table read by the names the
    ## explainer holds.  The values already held are replaced, not added to.
    ##
    ## @seealso{shapley}
    ## @end deftypefn

    function this = fit (this, QueryPoints)

      ## Input validation
      if (nargin != 2)
        error ("shapley.fit: invalid number of input arguments.");
      endif
      M = numel (this.PredictorNames);

      ## A table is read by the names the explainer holds, as X was
      if (istable (QueryPoints))
        QueryPoints = tableToMatrix (this.PredictorNames, ...
                                     shapLevels (this.BlackboxModel), ...
                                     QueryPoints, 'shapley.fit');
      endif

      if (! (isnumeric (QueryPoints) && isreal (QueryPoints)
             && ismatrix (QueryPoints) && ndims (QueryPoints) == 2
             && ! isempty (QueryPoints)))
        error (strcat ("shapley.fit: QUERYPOINTS must be a real numeric", ...
                       " matrix."));
      endif
      if (columns (QueryPoints) != M)
        error (strcat ("shapley.fit: QUERYPOINTS must have one column per", ...
                       " predictor of the model."));
      endif

      Xs = this.X(this.SampledObservationIndices,:);
      sfcn = shapScoreFcn (this.BlackboxModel, this.IsClass);
      K = columns (sfcn (Xs(1,:)));
      nq = rows (QueryPoints);

      ## A weighted sum of the predictors, and a decision tree, are each
      ## answered from what the model is rather than over its subsets.  A
      ## shortcut is checked against the whole deviation before it is kept,
      ## so a model read wrongly falls back rather than answering wrongly.
      linear = strcmp (this.Method, 'interventional-linear');
      tree = strcmp (this.Method, 'interventional-tree');
      exact = (M <= 20 && this.NumSubsets == 2 ^ M);
      Tree = [];
      leafval = [];
      Z = [];
      treephi = [];
      if (tree)
        ## The compiled kernel answers every query point in one call
        [Tree, leafval] = shapTreesOf (this.BlackboxModel);
        treephi = reshape (__shapleytree__ (shapTreeStructs (Tree, leafval), ...
                                            Xs, QueryPoints), M, K, nq);
        tree = shapAccounts (sum (treephi(:,:,1), 1), sfcn, Xs, ...
                             QueryPoints(1,:));
      elseif (linear)
        Z = shapLinearMasks (M);
        V1 = shapSubsetValues (sfcn, Xs, QueryPoints(1,:), Z, K, []);
        parts = sum (V1(3:end,:) - repmat (V1(1,:), M, 1), 1);
        linear = shapAccounts (parts, sfcn, Xs, QueryPoints(1,:));
      endif

      ## Whatever is left over is answered over the subsets
      if (! (linear || tree))
        if (strcmp (this.Method, 'interventional-linear')
            || strcmp (this.Method, 'interventional-tree'))
          this.Method = 'interventional-kernel';
        endif
        if (exact)
          Z = shapAllMasks (M);
        else
          Z = shapSelectMasks (M, this.NumSubsets);
          if (this.NumSubsets < 2 * M + 2)
            warning (strcat ("shapley.fit: the values may be unreliable", ...
                             " because 'MaxNumSubsets' is too small."));
          endif
        endif
      endif

      ## Conditioning stands the nearest observations in for the rest of the
      ## distribution, so what it needs is measured once over the sample
      cnd = shapConditioning (Xs, this.Method);

      phi = zeros (M, K, nq);
      icept = [];
      if (tree)
        phi = treephi;
        icept = mean (sfcn (Xs), 1);
      endif
      for ii = 1:nq
        if (tree)
          break;
        endif
        if (! isempty (cnd))
          cnd.qz = (QueryPoints(ii,:) - cnd.Mu) ./ cnd.Sigma;
        endif
        V = shapSubsetValues (sfcn, Xs, QueryPoints(ii,:), Z, K, cnd);
        v0 = V(1,:);
        if (linear)
          phi(:,:,ii) = V(3:end,:) - repmat (v0, M, 1);
        elseif (exact)
          phi(:,:,ii) = shapFromValues (V, M, K);
        else
          phi(:,:,ii) = shapKernelSolve (Z, V, M, K);
        endif
        if (ii == 1)
          icept = v0;
        endif
      endfor

      this.QueryPoints = QueryPoints;
      this.Intercept = icept;
      this.BlackboxFitted = shapFitted (this.BlackboxModel, QueryPoints, ...
                                        this.IsClass);
      this.Shapley = shapTable (phi, this.PredictorNames, this.ClassNames, ...
                                this.IsClass, false);
      this.MeanAbsoluteShapley = shapTable (phi, this.PredictorNames, ...
                                            this.ClassNames, this.IsClass, ...
                                            true);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {shapley} {} plot (@var{obj})
    ## @deftypefnx {shapley} {} plot (@var{obj}, @var{name}, @var{value})
    ## @deftypefnx {shapley} {} plot (@var{ax}, @dots{})
    ## @deftypefnx {shapley} {@var{b} =} plot (@dots{})
    ##
    ## Plot the Shapley values as a horizontal bar chart.
    ##
    ## One bar per predictor, the least important at the bottom.  Over one
    ## query point the bars hold the values themselves and the chart is
    ## titled @qcode{'Shapley Explanation'}; over several they hold the mean
    ## of the absolute values and it is titled
    ## @qcode{'Shapley Importance Plot'}.
    ##
    ## @var{ax} is the axes to draw into, the current one where none is
    ## given.  @var{b} holds one bar series per class drawn.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'NumImportantPredictors'} @tab @tab How many predictors
    ## to draw on their own, the ten most important by default.  Over
    ## several query points whatever is left over is drawn as one further
    ## bar holding its sum; over one query point it is left out.
    ##
    ## @item @qcode{'ClassNames'} @tab @tab The classes to draw, for a
    ## classification model.  The default is the predicted class over one
    ## query point and every class over several.
    ##
    ## @item @qcode{'QueryPointIndices'} @tab @tab Which query points to
    ## draw, all of them by default.
    ## @end multitable
    ##
    ## @seealso{shapley, shapley.boxchart, shapley.swarmchart}
    ## @end deftypefn

    function varargout = plot (varargin)

      [ax, this, args] = shapDrawArgs (varargin, 'plot');
      optNames = {'NumImportantPredictors', 'ClassNames', ...
                  'QueryPointIndices'};
      [nip, cn, qpi, rem] = parsePairedArguments (optNames, {[], [], []}, ...
                                                  args);
      if (! isempty (rem))
        error (strcat ("shapley.plot: unknown optional argument or", ...
                       " misplaced value."));
      endif

      cnames = shapClassVars (this);
      qpi = shapCheckPoints (qpi, rows (this.QueryPoints), 'plot');
      idx = shapCheckClasses (cnames, cn, 'plot', 'ClassNames', ...
                    shapDefaultClasses (cnames, this.BlackboxFitted, qpi));
      V = shapDrawValues (this.Shapley, cnames, idx, qpi);
      M = rows (V);
      nip = shapCheckCount (nip, M, 'plot');

      ## Over one query point the values stand as they are; over several it
      ## is the mean of their absolute values that MATLAB draws
      one = (numel (qpi) == 1);
      if (one)
        vals = reshape (V(:,1,:), M, []);
      else
        vals = reshape (mean (abs (V), 2), M, []);
      endif

      ## The least important at the bottom, and over several query points
      ## whatever is left over is summed into one bar of its own
      [~, ord] = sort (sum (abs (vals), 2), 'ascend');
      names = cellstr (this.Shapley.Predictor);
      names = names(ord);
      vals = vals(ord,:);
      if (nip < M)
        kept = (M - nip + 1):M;
        if (one)
          vals = vals(kept,:);
          names = names(kept);
        else
          lump = sum (vals(1:(M - nip),:), 1);
          vals = [lump; vals(kept,:)];
          label = sprintf ('Sum of other %d predictor(s)', M - nip);
          names = [{label}; names(kept)];
        endif
      endif

      if (isempty (ax))
        ax = gca ();
      endif
      n = rows (vals);
      h = barh (ax, 1:n, vals);
      set (ax, 'ytick', 1:n, 'yticklabel', names);
      if (one)
        title (ax, 'Shapley Explanation');
        xlabel (ax, 'Shapley Value');
      else
        title (ax, 'Shapley Importance Plot');
        xlabel (ax, 'Mean of Absolute Shapley Values');
      endif
      ylabel (ax, 'Predictor');
      if (numel (idx) > 1)
        legend (ax, cnames{idx});
      endif

      if (nargout > 0)
        varargout{1} = h;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {shapley} {} boxchart (@var{obj})
    ## @deftypefnx {shapley} {} boxchart (@var{obj}, @var{name}, @var{value})
    ## @deftypefnx {shapley} {} boxchart (@var{ax}, @dots{})
    ## @deftypefnx {shapley} {@var{b} =} boxchart (@dots{})
    ##
    ## Draw a box chart of the Shapley values over the query points.
    ##
    ## One box per predictor, the least important at the bottom, spread over
    ## the query points the values were fitted at.  The chart is titled
    ## @qcode{'Shapley Summary Plot'} and lies horizontally.
    ##
    ## @var{ax} is the axes to draw into, the current one where none is
    ## given.  @var{b} is the @code{stats.chart.BoxChart} drawn.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'NumImportantPredictors'} @tab @tab How many predictors
    ## to draw, the ten most important by default.  Whatever is left over is
    ## left out rather than summed, a box over a sum meaning nothing.
    ##
    ## @item @qcode{'ClassName'} @tab @tab The one class to draw, for a
    ## classification model.  The default is the first class of the model.
    ##
    ## @item @qcode{'JitterOutliers'} @tab @tab Whether outlier markers are
    ## spread across the width of the box, @qcode{'off'} by default.
    ## @end multitable
    ##
    ## @seealso{shapley, shapley.plot, shapley.swarmchart, stats.chart.BoxChart}
    ## @end deftypefn

    function varargout = boxchart (varargin)

      [ax, this, args] = shapDrawArgs (varargin, 'boxchart');
      optNames = {'NumImportantPredictors', 'ClassName', 'JitterOutliers'};
      [nip, cn, jit, rem] = parsePairedArguments (optNames, {[], [], []}, ...
                                                  args);
      if (! isempty (rem))
        error (strcat ("shapley.boxchart: unknown optional argument or", ...
                       " misplaced value."));
      endif

      if (! isempty (jit))
        jit = shapCheckWord (jit, {'on', 'off'}, 'boxchart', ...
                             'JitterOutliers');
      endif
      [vals, names] = shapOneClass (this.Shapley, shapClassVars (this), ...
                                    cn, nip, 'boxchart');
      n = rows (vals);
      nq = columns (vals);

      if (isempty (ax))
        ax = gca ();
      endif
      grp = categorical (repmat (names, 1, nq)(:), names);
      args = {'Orientation', 'horizontal', 'BoxWidth', 0.8};
      if (! isempty (jit))
        args = [args, {'JitterOutliers', jit}];
      endif
      b = boxchart (ax, grp, vals(:), args{:});
      set (ax, 'ytick', 1:n, 'yticklabel', names);
      title (ax, 'Shapley Summary Plot');
      xlabel (ax, 'Shapley Value');
      ylabel (ax, 'Predictor');

      if (nargout > 0)
        varargout{1} = b;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {shapley} {} swarmchart (@var{obj})
    ## @deftypefnx {shapley} {} swarmchart (@var{obj}, @var{name}, @var{value})
    ## @deftypefnx {shapley} {} swarmchart (@var{ax}, @dots{})
    ## @deftypefnx {shapley} {@var{s} =} swarmchart (@dots{})
    ##
    ## Draw a swarm chart of the Shapley values over the query points.
    ##
    ## One row of points per predictor, the least important at the bottom,
    ## one point per query point spread vertically by how crowded its
    ## neighbourhood is.  Each point is coloured by the value the predictor
    ## takes at that query point, the least of them at one end of the
    ## colour map and the greatest at the other.  The chart is titled
    ## @qcode{'Shapley Summary Plot'}.
    ##
    ## @var{ax} is the axes to draw into, the current one where none is
    ## given.  @var{s} holds one scatter object per predictor drawn.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'NumImportantPredictors'} @tab @tab How many predictors
    ## to draw, the ten most important by default.
    ##
    ## @item @qcode{'ClassName'} @tab @tab The one class to draw, for a
    ## classification model.  The default is the first class of the model.
    ##
    ## @item @qcode{'YJitter'} @tab @tab How the points of a row are spread,
    ## @qcode{'density'} by default, or @qcode{'rand'}, @qcode{'randn'} or
    ## @qcode{'none'}.
    ##
    ## @item @qcode{'ColorMap'} @tab @tab The colour map the predictor
    ## values are read through, as a name or as a matrix of one RGB triplet
    ## per row.  The default is the one the axes already carries.
    ## @end multitable
    ##
    ## @seealso{shapley, shapley.plot, shapley.boxchart, swarmchart}
    ## @end deftypefn

    function varargout = swarmchart (varargin)

      [ax, this, args] = shapDrawArgs (varargin, 'swarmchart');
      optNames = {'NumImportantPredictors', 'ClassName', 'YJitter', ...
                  'ColorMap'};
      [nip, cn, yj, cmap, rem] = parsePairedArguments (optNames, ...
                                                  {[], [], [], []}, args);
      if (! isempty (rem))
        error (strcat ("shapley.swarmchart: unknown optional argument or", ...
                       " misplaced value."));
      endif
      if (isempty (yj))
        yj = 'density';
      else
        yj = shapCheckWord (yj, {'none', 'density', 'rand', 'randn'}, ...
                            'swarmchart', 'YJitter');
      endif
      if (! isempty (cmap))
        cmap = shapCheckMap (cmap, 'swarmchart');
      endif

      [vals, names, ord] = shapOneClass (this.Shapley, ...
                              shapClassVars (this), cn, nip, 'swarmchart');
      n = rows (vals);
      nq = columns (vals);

      if (isempty (ax))
        ax = gca ();
      endif
      held = ishold (ax);
      hold (ax, 'on');
      h = zeros (n, 1);
      unwind_protect
        for k = 1:n
          col = shapColorValues (this.QueryPoints(:,ord(k)));
          at = k * ones (nq, 1);
          h(k) = swarmchart (ax, vals(k,:)', at, 36, col, ...
                             'XJitter', 'none', 'YJitter', yj);
        endfor
      unwind_protect_cleanup
        if (! held)
          hold (ax, 'off');
        endif
      end_unwind_protect
      set (ax, 'ytick', 1:n, 'yticklabel', names, 'clim', [0, 1]);
      if (! isempty (cmap))
        colormap (ax, cmap);
      endif
      title (ax, 'Shapley Summary Plot');
      xlabel (ax, 'Shapley Value');
      ylabel (ax, 'Predictor');

      if (nargout > 0)
        varargout{1} = h;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {shapley} {} plotDependence (@var{obj}, @var{predictor})
    ## @deftypefnx {shapley} {} plotDependence (@var{obj}, @var{predictor}, @var{name}, @var{value})
    ## @deftypefnx {shapley} {} plotDependence (@var{ax}, @dots{})
    ## @deftypefnx {shapley} {@var{p} =} plotDependence (@dots{})
    ##
    ## Draw the Shapley values of one predictor against its own values.
    ##
    ## @var{predictor} names or indexes the predictor.  Where it holds
    ## numbers the chart is a scatter of its value at each query point
    ## against the value it was given there; where it holds levels the chart
    ## is a box of the values at each level.
    ##
    ## The chart is titled @qcode{'Shapley Dependence Plot'} and its
    ## vertical axis is labelled after the predictor drawn.
    ##
    ## @var{ax} is the axes to draw into, the current one where none is
    ## given.  @var{p} is the scatter object, or the
    ## @code{stats.chart.BoxChart}, drawn.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'ClassName'} @tab @tab The one class to draw, for a
    ## classification model.  The default is the first class of the model.
    ##
    ## @item @qcode{'ColorPredictor'} @tab @tab A second predictor to colour
    ## the points by, none by default.  Its values are read as they stand,
    ## the range of the axes carrying the scale, and a colour bar is drawn
    ## beside the chart.  It applies only where the predictor drawn holds
    ## numbers.
    ##
    ## @item @qcode{'ColorMap'} @tab @tab The colour map the colouring
    ## predictor is read through, as a name or as a matrix of one RGB
    ## triplet per row.  The default is the one the axes already carries.
    ## @end multitable
    ##
    ## @seealso{shapley, shapley.plot, shapley.swarmchart}
    ## @end deftypefn

    function varargout = plotDependence (varargin)

      [ax, this, args] = shapDrawArgs (varargin, 'plotDependence');
      if (isempty (args))
        error ("shapley.plotDependence: too few input arguments.");
      endif
      pred = args{1};
      optNames = {'ClassName', 'ColorPredictor', 'ColorMap'};
      [cn, cp, cmap, rem] = parsePairedArguments (optNames, {[], [], []}, ...
                                                  args(2:end));
      if (! isempty (rem))
        error (strcat ("shapley.plotDependence: unknown optional argument", ...
                       " or misplaced value."));
      endif

      if (! isempty (cmap))
        cmap = shapCheckMap (cmap, 'plotDependence');
      endif
      names = cellstr (this.Shapley.Predictor);
      j = shapPredictorIndex (names, pred, 'plotDependence', 'PREDICTOR');
      cnames = shapClassVars (this);
      idx = shapCheckClasses (cnames, cn, 'plotDependence', 'ClassName', ...
                              shapFirstClass (cnames));
      if (numel (idx) > 1)
        error ("shapley.plotDependence: 'ClassName' must name one class.");
      endif
      V = shapDrawValues (this.Shapley, cnames, idx, ...
                          1:rows (this.QueryPoints));
      y = V(j,:)';
      x = this.QueryPoints(:,j);
      lvl = any (this.CategoricalPredictors == j);

      if (isempty (ax))
        ax = gca ();
      endif
      if (lvl)
        if (! isempty (cp))
          error (strcat ("shapley.plotDependence: 'ColorPredictor' does", ...
                         " not apply to a predictor holding levels."));
        endif
        p = boxchart (ax, x, y);
      elseif (isempty (cp))
        p = scatter (ax, x, y, 36);
      else
        k = shapPredictorIndex (names, cp, 'plotDependence', ...
                                'ColorPredictor');
        ## The colouring predictor is read as it stands here, the range of
        ## the axes carrying the scale, where a swarm normalizes instead
        col = this.QueryPoints(:,k);
        p = scatter (ax, x, y, 36, col);
        lo = min (col);
        hi = max (col);
        if (hi > lo)
          set (ax, 'clim', [lo, hi]);
        endif
        if (! isempty (cmap))
          colormap (ax, cmap);
        endif
        colorbar (ax);
      endif
      title (ax, 'Shapley Dependence Plot');
      xlabel (ax, names{j});
      ylabel (ax, sprintf ('Shapley Values for %s', names{j}));

      if (nargout > 0)
        varargout{1} = p;
      endif

    endfunction

  endmethods

  methods (Access = private)

    ## The variables the Shapley table holds one per class, empty for a
    ## regression model or a function handle, which hold one value only.
    function cnames = shapClassVars (this)

      cnames = {};
      if (this.IsClass)
        cnames = shapClassText (this.ClassNames);
      endif

    endfunction

  endmethods

endclassdef

## A table as the numbers the explainer works over.  A fitted model names
## the predictors and holds the levels they were coded through; a function
## handle has neither, so the table itself names them.
function [X, tnames, errmsg] = shapReadTable (blackbox, isfh, T, caller)

  X = [];
  tnames = {};
  errmsg = '';
  if (isfh)
    [X, ~, args, ~, errmsg] = tableFrame (T, [], {});
    if (isempty (errmsg))
      tnames = args{2};
    endif
    return;
  endif
  props = properties (blackbox);
  if (! any (strcmp (props, 'PredictorNames')))
    errmsg = strcat ("X may be a table only for a model that names its", ...
                     " predictors.");
    return;
  endif
  ## PredictorLevels is hidden, so properties () does not list it; the
  ## class that defines it is what says whether it is there
  lev = {};
  if (isa (blackbox, 'PredictiveModel'))
    lev = blackbox.PredictorLevels;
  endif
  try
    X = tableToMatrix (blackbox.PredictorNames, lev, T, caller);
  catch err
    errmsg = regexprep (err.message, '^shapley: ', '');
    return;
  end_try_catch

endfunction

## The levels a model's predictors were coded through, where it kept any.
function lev = shapLevels (blackbox)

  lev = {};
  if (isa (blackbox, 'PredictiveModel'))
    lev = blackbox.PredictorLevels;
  endif

endfunction

## The algorithm asked for, of those implemented.
function [base, errmsg] = shapMethod (Method)

  base = 'interventional';
  errmsg = '';
  if (isempty (Method))
    return;
  endif
  if (! (ischar (Method) || isa (Method, 'string')) || ! isrow (char (Method)))
    errmsg = "'Method' must be a character vector or a string scalar.";
    return;
  endif
  switch (lower (char (Method)))
    case 'interventional'
      base = 'interventional';
    case 'conditional'
      base = 'conditional';
    otherwise
      errmsg = strcat ("'Method' must be 'interventional' or", ...
                       " 'conditional'.");
  endswitch

endfunction

## Which algorithm answers for this model.  A model whose prediction is a
## weighted sum of its predictors is answered outright, and everything else
## through the subsets.  Asking for a budget of subsets asks for the subsets.
function method = shapAlgorithm (base, F, MaxSub)

  if (strcmp (base, 'conditional'))
    method = 'conditional-kernel';
    return;
  endif
  method = 'interventional-kernel';
  if (! (F.IsLinear || F.IsTree))
    return;
  endif
  if (! isempty (MaxSub))
    warning (strcat ("shapley: 'MaxNumSubsets' is given, so the values", ...
                     " are taken over subsets rather than from the", ...
                     " structure of the model."));
    return;
  endif
  if (F.IsLinear)
    method = 'interventional-linear';
    return;
  endif

  method = 'interventional-tree';

endfunction

## Resolve the observations, the predictors and the classes.
function [F, errmsg] = shapFrame (blackbox, Data, CatPred, NumObs, MaxSub)

  F = [];
  errmsg = '';
  isfh = is_function_handle (blackbox);
  if (isfh)
    props = {};
  else
    props = properties (blackbox);
  endif
  has = @(n) any (strcmp (props, n));

  ## The observations averaged over, from the model where it kept them
  if (isempty (Data))
    if (isfh)
      errmsg = "X is required when the model is a function handle.";
      return;
    endif
    if (! (has ('X') && ! isempty (blackbox.X)))
      errmsg = strcat ("X is required for a model that does not keep the", ...
                       " observations it was fitted on.");
      return;
    endif
    Data = blackbox.X;
  endif

  ## A table is read by the names the model was fitted on, so a column it
  ## was not fitted on is passed over and a level carries the code it
  ## carried then.  A handle has no model behind it, so its table names its
  ## own predictors and every column is one.
  tnames = {};
  if (istable (Data))
    [Data, tnames, errmsg] = shapReadTable (blackbox, isfh, Data, 'shapley');
    if (! isempty (errmsg))
      return;
    endif
  endif

  if (! (isnumeric (Data) && isreal (Data) && ismatrix (Data)
         && ndims (Data) == 2 && ! isempty (Data)))
    errmsg = "X must be a real numeric matrix.";
    return;
  endif
  p = columns (Data);

  ## The predictor names, from the model where it has them
  if (has ('PredictorNames') && ! isempty (blackbox.PredictorNames))
    pnames = blackbox.PredictorNames(:)';
    if (numel (pnames) != p)
      errmsg = "X must have one column per predictor of the model.";
      return;
    endif
  elseif (! isempty (tnames))
    pnames = tnames;
  else
    pnames = arrayfun (@(k) sprintf ('x%d', k), 1:p, 'UniformOutput', false);
  endif

  ## The categorical predictors: the model's own, or those named for a handle
  if (isfh)
    spec = CatPred;
  else
    if (has ('CategoricalPredictors'))
      spec = blackbox.CategoricalPredictors;
    else
      spec = [];
    endif
    if (! isempty (CatPred) && ! isequal (CatPred, spec))
      errmsg = strcat ("'CategoricalPredictors' must match the", ...
                       " CategoricalPredictors property of the model.");
      return;
    endif
  endif
  cat = [];
  if (! isempty (spec))
    [C, msg] = dummyCoding (Data, spec, pnames);
    if (! isempty (msg))
      errmsg = msg;
      return;
    endif
    cat = C.Index;
  endif

  ## How many subsets to compute over.  Every one of them gives the values
  ## exactly, and MATLAB stops at 1024 of them by default.
  full = 2 ^ p;
  if (isempty (MaxSub))
    nsub = min (full, 1024);
  elseif (isnumeric (MaxSub) && isscalar (MaxSub) && isreal (MaxSub)
          && MaxSub > 1 && MaxSub == fix (MaxSub))
    nsub = min (full, double (MaxSub));
  else
    errmsg = "'MaxNumSubsets' must be an integer greater than 1.";
    return;
  endif

  ## The observations drawn from those given
  n = rows (Data);
  if (isempty (NumObs))
    k = min (100, n);
  elseif ((ischar (NumObs) || isa (NumObs, 'string'))
          && strcmpi (char (NumObs), 'all'))
    k = n;
  elseif (isnumeric (NumObs) && isscalar (NumObs) && isreal (NumObs)
          && NumObs >= 1 && NumObs == fix (NumObs))
    k = min (double (NumObs), n);
  else
    errmsg = strcat ("'NumObservationsToSample' must be a positive", ...
                     " integer or 'all'.");
    return;
  endif
  if (k < n)
    idx = sort (randperm (n, k))(:);
  else
    idx = (1:n)(:);
  endif

  ## Assigned field by field: struct () spreads a cell argument into a struct
  ## array, and a cellstr ClassNames would make one
  F.IsClass = has ('ClassNames') && ! isempty (blackbox.ClassNames);
  if (F.IsClass)
    F.ClassNames = blackbox.ClassNames;
  else
    F.ClassNames = [];
  endif
  F.PredictorNames = pnames;
  F.NumPredictors = p;
  F.NumSubsets = nsub;
  F.IsLinear = shapIsLinear (blackbox, isfh, has);
  F.NumTrees = numel (shapTreesOf (blackbox));
  F.IsTree = (F.NumTrees > 0);
  F.Cat = cat;
  F.X = Data;
  F.Idx = idx;

endfunction

## The decision trees a model answers through, and what each answers at each
## of its nodes, already carrying whatever weight the model gives it.  Empty
## where the model does not come apart into trees, which is what sends it to
## the subsets instead.
##
## A single tree answers with its own leaves.  An ensemble answers with the
## weighted sum of its learners' outputs, so it comes apart too, provided
## every learner is a tree, the learners are combined by a constant weight
## each, and nothing is applied to the total afterwards.
function [trees, leafvals] = shapTreesOf (blackbox)

  trees = {};
  leafvals = {};
  if (is_function_handle (blackbox))
    return;
  endif

  switch (class (blackbox))

    case {'ClassificationTree', 'CompactClassificationTree', ...
          'RegressionTree', 'CompactRegressionTree'}
      trees = {blackbox};
      leafvals = {shapLeafValues(blackbox)};

    case {'ClassificationEnsemble', 'CompactClassificationEnsemble', ...
          'ClassificationBaggedEnsemble', 'RegressionEnsemble', ...
          'CompactRegressionEnsemble', 'RegressionBaggedEnsemble'}
      [trees, leafvals] = shapEnsembleTrees (blackbox);

  endswitch

endfunction

## What one tree answers at each of its nodes, transformed as it answers.
function L = shapLeafValues (T)

  props = properties (T);
  if (any (strcmp (props, 'ClassProbability')))
    L = T.ClassProbability;
    if (! isempty (T.STfun))
      L = T.STfun (L);
    endif
  else
    L = T.NodeMean(:);
    if (! isempty (T.RTfun))
      L = T.RTfun (L);
    endif
  endif

endfunction

## The learners of an ensemble and what each contributes to its answer.
##
## A regression ensemble adds up what its learners predict, whatever method
## grew them.  A classification ensemble adds up their class scores only for
## the methods that read a learner that way; the boosting methods that read a
## learner's label, or its margin, contribute something that is not the
## learner's own scores, and those are left to the subsets.  Anything applied
## to the total, a transform of the score or of the response, would not come
## apart and is left to the subsets too.
function [trees, leafvals] = shapEnsembleTrees (Mdl)

  trees = {};
  leafvals = {};
  props = properties (Mdl);
  has = @(n) any (strcmp (props, n));
  if (! has ('Trained') || isempty (Mdl.Trained))
    return;
  endif
  isclass = has ('ClassNames') && ! isempty (Mdl.ClassNames);
  if (isclass)
    trans = Mdl.ScoreTransform;
    ## The methods whose learner output is the learner's own class scores
    if (! any (strcmp (Mdl.Method, {'Bag', 'AdaBoostM2', 'RUSBoost'})))
      return;
    endif
  else
    trans = Mdl.ResponseTransform;
  endif
  if (! (ischar (trans) && any (strcmpi (trans, {'none', 'identity'}))))
    return;
  endif

  learners = Mdl.Trained(:)';
  T = numel (learners);
  istree = cellfun (@(t) any (strcmp (properties (t), 'Children')), learners);
  if (! all (istree))
    return;
  endif

  w = Mdl.TrainedWeights(:)';
  if (isempty (w))
    w = ones (1, T);
  endif
  if (strcmp (Mdl.CombineWeights, 'WeightedAverage'))
    tot = sum (w);
    if (! (tot > 0))
      return;
    endif
    w = w / tot;
  endif

  if (isclass)
    K = numel (Mdl.ClassNames);
  endif
  leafvals = cell (1, T);
  for t = 1:T
    L = shapLeafValues (learners{t});
    if (isclass)
      ## A learner need not know every class of the ensemble, so its columns
      ## are put where the ensemble keeps them
      col = labelIndices (Mdl.ClassNames, learners{t}.ClassNames);
      if (any (col == 0))
        trees = {};
        leafvals = {};
        return;
      endif
      G = zeros (rows (L), K);
      G(:,col) = L;
      L = G;
    endif
    leafvals{t} = w(t) * L;
  endfor
  trees = learners;

endfunction

## The trees as the compiled kernel takes them: plain node tables, with the
## levels of a categorical cut carried beside the thresholds.
function C = shapTreeStructs (trees, leafvals)

  C = cell (1, numel (trees));
  for ti = 1:numel (trees)
    T = trees{ti};
    st = struct ();
    st.Children = double (T.Children);
    st.CutVar = double (T.CutPredictorIndex(:));
    st.CutPoint = double (T.CutPoint(:));
    st.Leaf = double (leafvals{ti});
    if (! isempty (T.CutCategories))
      st.CatLeft = T.CutCategories(:,1);
      st.CatRight = T.CutCategories(:,2);
    else
      st.CatLeft = {};
      st.CatRight = {};
    endif
    C{ti} = st;
  endfor

endfunction

## Whether the model predicts a weighted sum of its predictors.  Beta holds
## the weights of the linear learners and of a support vector machine with a
## linear kernel, and is empty for every other kernel, so its presence is the
## question.  A transform of the score breaks the sum, and a model carrying
## one is answered through the subsets instead.
function tf = shapIsLinear (blackbox, isfh, has)

  tf = false;
  if (isfh || ! has ('Beta') || isempty (blackbox.Beta))
    return;
  endif
  if (has ('ScoreTransform'))
    trans = blackbox.ScoreTransform;
  elseif (has ('ResponseTransform'))
    trans = blackbox.ResponseTransform;
  else
    trans = 'none';
  endif
  tf = (ischar (trans) && any (strcmpi (trans, {'none', 'identity'})));

endfunction

## What the model answers for a set of observations, as a column per class.
function sfcn = shapScoreFcn (blackbox, isclass)

  if (is_function_handle (blackbox))
    sfcn = @(Z) shapHandleScore (blackbox, Z);
  elseif (isclass)
    sfcn = @(Z) shapClassScore (blackbox, Z);
  else
    sfcn = @(Z) shapResponse (blackbox, Z);
  endif

endfunction

function s = shapHandleScore (fun, Z)

  s = fun (Z);
  if (! (isnumeric (s) && isreal (s) && columns (s) == 1
         && rows (s) == rows (Z)))
    error (strcat ("shapley: the function must answer with one real", ...
                   " numeric column holding one value per observation."));
  endif
  s = double (s);

endfunction

function s = shapClassScore (Mdl, Z)

  [~, s] = predict (Mdl, Z);
  s = double (s);

endfunction

function s = shapResponse (Mdl, Z)

  s = double (predict (Mdl, Z))(:);

endfunction

## Every subset of the predictors, as one logical row each.  Bit II of the row
## index, counted from zero, says whether predictor II is in the subset, which
## is the order shapFromValues reads them in.
function Z = shapAllMasks (M)

  nS = 2 ^ M;
  Z = false (nS, M);
  for ii = 1:M
    Z(:,ii) = logical (bitget ((0:(nS - 1))', ii));
  endfor

endfunction

## The subsets a weighted sum of the predictors is answered from: the empty
## set, the full one, and each predictor on its own.  Holding one predictor
## at the query point moves such a prediction by that predictor's whole
## contribution and by nothing else, so the difference from the average is
## the Shapley value itself and the other subsets say nothing new.
function Z = shapLinearMasks (M)

  none = false (1, M);
  every = true (1, M);
  singles = logical (eye (M));
  Z = [none; every; singles];

endfunction

## Whether a shortcut's values account for the whole deviation of the
## prediction from the average, which is what says the shortcut applies.
function tf = shapAccounts (parts, sfcn, Xs, q)

  dev = sfcn (q) - mean (sfcn (Xs), 1);
  tf = all (abs (parts - dev) <= 1e-8 * max (1, max (abs (dev))));

endfunction

## The subsets to compute over when there is no room for all of them: the
## empty set and the full one, which carry the average prediction and the
## prediction itself, and then as many as the budget allows in decreasing
## order of their kernel weight.  That weight depends only on how many
## predictors a subset holds and falls as it moves away from the two
## extremes, so the cardinalities are taken in the pairs (1, M-1), (2, M-2)
## and so on, the smaller one lexicographically and the larger one in the
## order of the complements.
##
## The single predictors are taken as far as the budget reaches whatever is
## left of it.  A later pair is taken entire while it fits, and the first
## that does not is filled with complementary pairs of subsets instead, so
## that what is left of the budget still spans both ends rather than sitting
## at one of them.  Measured against MATLAB R2024a at seven budgets for four
## predictors and ten for five.
function Z = shapSelectMasks (M, nSub)

  none = false (1, M);
  every = true (1, M);
  Z = [none; every];
  budget = nSub - 2;
  if (budget < 1 || M < 2)
    return;
  endif

  sel = zeros (budget, M);
  nr = 0;
  lo = 1;
  hi = M - 1;
  while (lo <= hi && nr < budget)
    left = budget - nr;
    nlo = shapLevelCount (M, lo, left);
    if (lo == hi)
      ## The middle cardinality is its own complement
      [sel, nr] = shapFillPairs (sel, nr, budget, M, lo, true);
    elseif (lo == 1 || 2 * nlo <= left)
      C = shapCombinations (M, lo, left);
      for ii = 1:rows (C)
        if (nr >= budget)
          break;
        endif
        nr++;
        sel(nr,:) = C(ii,:);
      endfor
      for ii = 1:rows (C)
        if (nr >= budget)
          break;
        endif
        nr++;
        sel(nr,:) = ! C(ii,:);
      endfor
    else
      [sel, nr] = shapFillPairs (sel, nr, budget, M, lo, false);
    endif
    lo++;
    hi--;
  endwhile
  taken = logical (sel(1:nr,:));
  Z = [Z; taken];

endfunction

## Fill what is left of the budget with complementary pairs of subsets of S
## predictors, each subset followed by the one holding the rest.  SELF says
## the complement has S predictors too, so a pair can be reached twice and the
## second reading of it is passed over.
function [sel, nr] = shapFillPairs (sel, nr, budget, M, s, self)

  C = shapCombinations (M, s, budget);
  used = false (rows (C), 1);
  for ii = 1:rows (C)
    if (nr >= budget)
      break;
    endif
    if (self && used(ii))
      continue;
    endif
    used(ii) = true;
    nr++;
    sel(nr,:) = C(ii,:);
    comp = ! C(ii,:);
    if (nr >= budget)
      break;
    endif
    if (self)
      jj = find (all (C == repmat (comp, rows (C), 1), 2), 1);
      if (isempty (jj) || used(jj))
        continue;
      endif
      used(jj) = true;
    endif
    nr++;
    sel(nr,:) = comp;
  endfor

endfunction

## How many subsets of S predictors there are out of M, counted up to CAP so
## that a wide model does not overflow the binomial coefficient it never
## needs the whole of.
function n = shapLevelCount (M, s, cap)

  n = 1;
  for ii = 1:s
    n = n * (M - s + ii) / ii;
    if (n > cap)
      n = cap + 1;
      return;
    endif
  endfor
  n = round (n);

endfunction

## The first N combinations of S indices out of M, in lexicographic order, as
## one logical row each.  They are generated rather than enumerated, so a wide
## model does not build a list it has no room for.
function C = shapCombinations (M, s, n)

  C = zeros (n, M);
  nr = 0;
  c = 1:s;
  while (nr < n)
    nr++;
    C(nr,c) = 1;
    jj = s;
    while (jj >= 1 && c(jj) == M - s + jj)
      jj--;
    endwhile
    if (jj < 1)
      break;
    endif
    c(jj)++;
    c((jj + 1):s) = c(jj) + (1:(s - jj));
  endwhile
  C = C(1:nr,:);

endfunction

## The value function over the given subsets.  A subset's columns are held at
## the query point and the rest stay as the observations have them, which is
## what makes the average an interventional one.
function V = shapSubsetValues (sfcn, Xs, q, Z, K, cnd)

  L = rows (Z);
  V = zeros (L, K);
  for ii = 1:L
    mask = Z(ii,:);
    if (isempty (cnd) || ! any (mask) || all (mask))
      ## Conditioning on nothing, or on everything, leaves the whole sample:
      ## the first averages over it and the second replaces all of it
      W = Xs;
    else
      W = Xs(shapNeighbours (cnd, mask),:);
    endif
    if (any (mask))
      W(:,mask) = repmat (q(mask), rows (W), 1);
    endif
    V(ii,:) = mean (sfcn (W), 1);
  endfor

endfunction

## What conditioning needs, measured once over the observations: the column
## means and deviations the distances are taken in, so that a predictor does
## not count for more merely because it is recorded on a wider scale, and how
## many neighbours stand in for the conditional distribution.
function cnd = shapConditioning (Xs, method)

  cnd = [];
  if (! strcmp (method, 'conditional-kernel'))
    return;
  endif
  n = rows (Xs);
  cnd.Mu = mean (Xs, 1);
  sd = std (Xs, 0, 1);
  sd(sd == 0 | ! isfinite (sd)) = 1;
  cnd.Sigma = sd;
  cnd.Xz = (Xs - repmat (cnd.Mu, n, 1)) ./ repmat (sd, n, 1);
  ## A tenth of the observations, rounded up.  Divided rather than multiplied
  ## by a tenth, which is not exact in binary and would round a whole tenth up
  cnd.NumNeighbors = max (1, ceil (n / 10));
  cnd.qz = [];

endfunction

## The observations whose values of the predictors in MASK lie nearest the
## query point, which are the ones conditioning on those predictors keeps.
function nb = shapNeighbours (cnd, mask)

  n = rows (cnd.Xz);
  D = cnd.Xz(:,mask) - repmat (cnd.qz(mask), n, 1);
  [~, ord] = sort (sum (D .^ 2, 2));
  nb = ord(1:cnd.NumNeighbors);

endfunction

## The Shapley values estimated from a budget of subsets, by the weighted
## least squares kernel SHAP solves.  Row 1 of Z is the empty subset and row 2
## the full one, so the values are tied to sum to the deviation of the
## prediction from the average; that is imposed by substitution rather than by
## a penalty, so it holds exactly however few subsets there are.
function phi = shapKernelSolve (Z, V, M, K)

  v0 = V(1,:);
  c = V(2,:) - v0;
  Zi = Z(3:end,:);
  L = rows (Zi);
  if (L == 0)
    ## Nothing but the two extremes says nothing about any one predictor
    phi = NaN (M, K);
    return;
  endif
  if (M == 1)
    phi = c;
    return;
  endif

  ## The Shapley kernel, in logs so that a wide model does not overflow the
  ## binomial coefficient, and only up to a constant since it is a weight
  s = sum (Zi, 2);
  lchoose = gammaln (M + 1) - gammaln (s + 1) - gammaln (M - s + 1);
  lw = log (M - 1) - lchoose - log (s) - log (M - s);
  w = exp (lw - max (lw));

  ## Substituting the last value out of the sum imposes the constraint
  last = Zi(:,M);
  A = Zi(:,1:(M - 1)) - repmat (last, 1, M - 1);
  b = (V(3:end,:) - repmat (v0, L, 1)) - last * c;
  sw = sqrt (w);

  ## A budget too small to determine every value leaves the system rank
  ## deficient, which is the case the caller has already warned about
  warning ("off", "Octave:singular-matrix", "local");
  warning ("off", "Octave:rank-deficient-matrix", "local");
  psi = (repmat (sw, 1, M - 1) .* A) \ (repmat (sw, 1, K) .* b);
  lastPhi = c - sum (psi, 1);
  phi = [psi; lastPhi];

endfunction

## The Shapley value of each predictor, summed over every subset that leaves
## it out.  The weight of a subset of S predictors is 1 / (M * C(M-1, S)),
## which is the S! (M-S-1)! / M! of the definition written so that it does not
## overflow.
function phi = shapFromValues (V, M, K)

  w = zeros (1, M);
  for s = 0:(M - 1)
    w(s + 1) = 1 / (M * nchoosek (M - 1, s));
  endfor

  phi = zeros (M, K);
  nS = 2 ^ M;
  for ii = 1:M
    bit = 2 ^ (ii - 1);
    for m = 0:(nS - 1)
      if (bitand (m, bit) == 0)
        s = sum (bitget (m, 1:M));
        phi(ii,:) += w(s + 1) * (V(m + bit + 1,:) - V(m + 1,:));
      endif
    endfor
  endfor

endfunction

## What the model answers at the query points.
function fitted = shapFitted (blackbox, Q, isclass)

  if (is_function_handle (blackbox))
    fitted = shapHandleScore (blackbox, Q);
  elseif (isclass)
    fitted = predict (blackbox, Q);
  else
    fitted = double (predict (blackbox, Q))(:);
  endif

endfunction

## The values as a table of one row per predictor, holding one column per
## query point, or their mean absolute value where MEANABS is true.
function T = shapTable (phi, pnames, cn, isclass, meanabs)

  M = rows (phi);
  K = columns (phi);
  if (meanabs)
    phi = mean (abs (phi), 3);
  endif
  P = string (pnames(:));

  if (! isclass)
    T = table (P, reshape (phi(:,1,:), M, []), ...
               'VariableNames', {'Predictor', 'Value'});
    return;
  endif

  names = shapClassText (cn);
  cols = cell (1, K + 1);
  cols{1} = P;
  for k = 1:K
    cols{k + 1} = reshape (phi(:,k,:), M, []);
  endfor
  T = table (cols{:}, 'VariableNames', [{'Predictor'}, names]);

endfunction

## The class names as the text a table variable is named with.
function names = shapClassText (cn)

  if (iscellstr (cn))
    names = cn(:)';
  elseif (ischar (cn))
    names = cellstr (cn)(:)';
  elseif (isa (cn, 'string') || isa (cn, 'categorical'))
    names = cellstr (cn)(:)';
  elseif (islogical (cn))
    names = arrayfun (@(v) mat2str (v), cn(:)', 'UniformOutput', false);
  else
    names = arrayfun (@(v) sprintf ('%g', v), cn(:)', 'UniformOutput', false);
  endif

endfunction

## The axes and the explainer, in whichever order the call put them.  A
## method reached with an axes first has the explainer second, dispatch
## having found it wherever it sits.
function [ax, this, args] = shapDrawArgs (in, caller)

  ax = [];
  if (numel (in) > 1 && isscalar (in{1}) && ishghandle (in{1})
      && isaxes (in{1}))
    ax = in{1};
    in(1) = [];
  endif
  this = in{1};
  args = in(2:end);
  if (isempty (this.Shapley))
    error (strcat ("shapley.%s: the Shapley values are not fitted; use", ...
                   " fit to compute them."), caller);
  endif

endfunction

## The values of the classes asked for at the query points asked for, as one
## page per class.  A regression model and a function handle hold one value
## per predictor and give one page.
function V = shapDrawValues (T, cnames, idx, qpi)

  if (isempty (cnames))
    V = T.Value(:,qpi);
    return;
  endif
  M = size (T, 1);
  V = zeros (M, numel (qpi), numel (idx));
  for k = 1:numel (idx)
    A = T.(cnames{idx(k)});
    V(:,:,k) = A(:,qpi);
  endfor

endfunction

## The classes drawn where none was named: the predicted class over one
## query point and every one of them over several, as MATLAB chooses them.
function idx = shapDefaultClasses (cnames, fitted, qpi)

  idx = [];
  if (isempty (cnames))
    return;
  endif
  if (numel (qpi) != 1)
    idx = 1:numel (cnames);
    return;
  endif
  at = shapClassText (fitted(qpi));
  idx = find (strcmp (cnames, at{1}), 1);
  if (isempty (idx))
    idx = 1;
  endif

endfunction

## The class a chart of one class draws where none was named.
function idx = shapFirstClass (cnames)

  idx = [];
  if (! isempty (cnames))
    idx = 1;
  endif

endfunction

## The classes named, as indices into the model's own order.
function idx = shapCheckClasses (cnames, want, caller, optname, dflt)

  if (isempty (want))
    idx = dflt;
    return;
  endif
  if (isempty (cnames))
    error (strcat ("shapley.%s: '%s' is valid only for a classification", ...
                   " model."), caller, optname);
  endif
  want = shapClassText (want);
  idx = zeros (1, numel (want));
  for k = 1:numel (want)
    j = find (strcmp (cnames, want{k}), 1);
    if (isempty (j))
      error ("shapley.%s: '%s' does not name a class of the model.", ...
             caller, want{k});
    endif
    idx(k) = j;
  endfor

endfunction

## The query points to draw, all of them where none were named.
function qpi = shapCheckPoints (qpi, nq, caller)

  if (isempty (qpi))
    qpi = 1:nq;
    return;
  endif
  if (! (isnumeric (qpi) && isreal (qpi) && isvector (qpi)
         && all (isfinite (qpi)) && all (qpi == fix (qpi)) && all (qpi >= 1)))
    error (strcat ("shapley.%s: 'QueryPointIndices' must be positive", ...
                   " integers."), caller);
  endif
  if (any (qpi > nq))
    error (strcat ("shapley.%s: 'QueryPointIndices' must not exceed the", ...
                   " %d query points fitted."), caller, nq);
  endif
  qpi = double (qpi(:)');

endfunction

## How many predictors to draw on their own, ten of them by default, and
## never more than the model has.
function nip = shapCheckCount (nip, M, caller)

  if (isempty (nip))
    nip = min (M, 10);
    return;
  endif
  if (! (isnumeric (nip) && isscalar (nip) && isreal (nip)
         && isfinite (nip) && nip == fix (nip) && nip > 0))
    error (strcat ("shapley.%s: 'NumImportantPredictors' must be a", ...
                   " positive integer."), caller);
  endif
  nip = min (double (nip), M);

endfunction

## The values of one class at every query point, the predictors ordered as
## the charts of one class order them and cut to the ones drawn.
function [vals, names, ord] = shapOneClass (T, cnames, cn, nip, caller)

  idx = shapCheckClasses (cnames, cn, caller, 'ClassName', ...
                          shapFirstClass (cnames));
  if (numel (idx) > 1)
    error ("shapley.%s: 'ClassName' must name one class.", caller);
  endif
  if (isempty (cnames))
    nq = columns (T.Value);
  else
    nq = columns (T.(cnames{1}));
  endif
  V = shapDrawValues (T, cnames, idx, 1:nq);
  M = size (T, 1);
  nip = shapCheckCount (nip, M, caller);
  [~, ord] = sort (mean (abs (V), 2), 'ascend');
  ord = ord((M - nip + 1):M);
  names = cellstr (T.Predictor);
  names = names(ord);
  vals = V(ord,:);

endfunction

## A predictor's values as the share of its range each one stands at, which
## is what a colour map is read with.  A predictor that never changes sits
## in the middle of the map, having no range to spread over.
function c = shapColorValues (v)

  v = double (v(:));
  lo = min (v);
  hi = max (v);
  if (hi > lo)
    c = (v - lo) / (hi - lo);
  else
    c = 0.5 * ones (size (v));
  endif

endfunction

## One of a short list of words, named by the method that took it rather
## than by whatever the drawing is handed on to.
function v = shapCheckWord (val, allowed, caller, optname)

  if (! (ischar (val) && isrow (val) && any (strcmpi (allowed, val))))
    error ("shapley.%s: '%s' must be one of %s.", caller, optname, ...
           strjoin (strcat ("'", allowed, "'"), ', '));
  endif
  j = find (strcmpi (allowed, val), 1);
  v = allowed{j};

endfunction

## A colour map, as a name or as one RGB triplet per row.
function v = shapCheckMap (val, caller)

  if (ischar (val) && isrow (val))
    v = val;
    return;
  endif
  if (isnumeric (val) && isreal (val) && ismatrix (val)
      && columns (val) == 3 && rows (val) > 0
      && all (val(:) >= 0) && all (val(:) <= 1))
    v = double (val);
    return;
  endif
  error (strcat ("shapley.%s: 'ColorMap' must be a colour map name or a", ...
                 " matrix of RGB triplets."), caller);

endfunction

## The one predictor named or indexed.
function j = shapPredictorIndex (names, pred, caller, argname)

  if (isnumeric (pred) && isscalar (pred) && isreal (pred)
      && isfinite (pred) && pred == fix (pred) && pred >= 1
      && pred <= numel (names))
    j = double (pred);
    return;
  endif
  if (ischar (pred) && isrow (pred))
    pred = {pred};
  elseif (isa (pred, 'string') && isscalar (pred))
    pred = cellstr (pred);
  endif
  if (iscellstr (pred) && isscalar (pred))
    j = find (strcmp (names, pred{1}), 1);
    if (! isempty (j))
      return;
    endif
  endif
  error ("shapley.%s: %s must name or index one predictor of the model.", ...
         caller, argname);

endfunction

%!demo
%! ## Explain a model fitted from a table
%!
%! load fisheriris
%! T = table (meas(:,2), meas(:,3), meas(:,4), meas(:,1), ...
%!            'VariableNames', {'SW', 'PL', 'PW', 'SL'});
%! T.Wide = categorical (meas(:,2) > 3, [false true], {'narrow', 'wide'});
%! Mdl = fitrtree (T, 'SL');
%!
%! ## The observations and the query points may be tables too, read by the
%! ## names the model was fitted on rather than by the order of the columns
%! s = shapley (Mdl, T(:, [1, 2, 3, 5]), 'QueryPoints', T(1, [5, 3, 2, 1]), ...
%!              'NumObservationsToSample', 'all');
%! s.Shapley
%!
%! ## A column the model was not fitted on is passed over, so the whole
%! ## table, response and all, gives the same explanation
%! t = shapley (Mdl, T, 'QueryPoints', T(1,:), ...
%!              'NumObservationsToSample', 'all');
%! isequal (t.Shapley.Value, s.Shapley.Value)

## A linear function's exact Shapley value is the coefficient times the
## deviation of the predictor from its mean, which is an expectation the
## engine cannot influence.
%!test
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! s = shapley (f, X, 'QueryPoints', [3, 20], ...
%!              'NumObservationsToSample', 'all');
%! b = [2, -3];
%! q = [3, 20];
%! assert_equal (s.Shapley.Value, (b .* (q - mean (X)))', 1e-12);

%!test  # the values sum to the deviation of the prediction from the average
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! s = shapley (f, X, 'QueryPoints', [3, 20], ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (sum (s.Shapley.Value), f ([3, 20]) - s.Intercept, 1e-12);

%!test  # the intercept is the average prediction over the observations
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! s = shapley (f, X, 'QueryPoints', [3, 20], ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Intercept, mean (f (X)), 1e-12);

%!test  # every subset is used, so the count is two to the predictors
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! f = @(Z) Z(:,1);
%! s = shapley (f, X);
%! assert_equal (s.NumSubsets, 4);

%!test  # a predictor the function ignores contributes nothing
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! f = @(Z) Z(:,1);
%! s = shapley (f, X, 'QueryPoints', [3, 20], ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.Value(2), 0, 1e-12);

%!test  # the algorithm is reported
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! s = shapley (@(Z) Z(:,1), X);
%! assert_equal (s.Method, 'interventional-kernel');

%!test  # a handle names its predictors as MATLAB does
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! s = shapley (@(Z) Z(:,1), X, 'QueryPoints', [3, 20]);
%! assert_equal (cellstr (s.Shapley.Predictor), {'x1'; 'x2'});

%!test  # the table names its one value variable as MATLAB R2026a does
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! s = shapley (@(Z) Z(:,1), X, 'QueryPoints', [3, 20]);
%! assert_equal (s.Shapley.Properties.VariableNames, {'Predictor', 'Value'});

%!test  # one column of values per query point
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! s = shapley (@(Z) Z(:,1), X, 'QueryPoints', [3, 20; 2, 30], ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (size (s.Shapley.Value), [2, 2]);

%!test  # with one query point the mean absolute value is the absolute value
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! s = shapley (f, X, 'QueryPoints', [3, 20], ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.MeanAbsoluteShapley.Value, abs (s.Shapley.Value), 1e-12);

%!test  # fit replaces the query points rather than adding to them
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! s = shapley (@(Z) Z(:,1), X, 'QueryPoints', [3, 20; 2, 30]);
%! s = fit (s, [1, 10]);
%! assert_equal (s.QueryPoints, [1, 10]);
%! assert_equal (size (s.Shapley.Value), [2, 1]);

%!test  # 'all' averages over every row, and says which ones
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! s = shapley (@(Z) Z(:,1), X, 'NumObservationsToSample', 'all');
%! assert_equal (s.SampledObservationIndices, (1:4)');

%!test  # a sample of the observations is drawn without replacement
%! X = (1:200)';
%! X = [X, 2 * X];
%! s = shapley (@(Z) Z(:,1), X);
%! assert_equal (numel (s.SampledObservationIndices), 100);
%! assert_equal (numel (unique (s.SampledObservationIndices)), 100);

%!test  # a handle answers what it was asked at the query points
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! s = shapley (f, X, 'QueryPoints', [3, 20]);
%! assert_equal (s.BlackboxFitted, f ([3, 20]), 1e-12);

%!test  # a model supplies the observations it was fitted on
%! load fisheriris
%! Mdl = fitrtree (meas(:,2:4), meas(:,1));
%! s = shapley (Mdl, 'NumObservationsToSample', 'all');
%! assert_equal (size (s.X), [150, 3]);

%!test  # and its own predictor names
%! load fisheriris
%! Mdl = fitrtree (meas(:,2:4), meas(:,1), ...
%!                 'PredictorNames', {'a', 'b', 'c'});
%! s = shapley (Mdl, 'QueryPoints', meas(1,2:4), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (cellstr (s.Shapley.Predictor), {'a'; 'b'; 'c'});

%!test  # a classifier names one variable per class
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.Properties.VariableNames, ...
%!               {'Predictor', 'setosa', 'versicolor', 'virginica'});

%!test  # the two classes of a binary classifier are exact negatives
%! load fisheriris
%! Mdl = fitctree (meas(51:150,:), species(51:150));
%! s = shapley (Mdl, 'QueryPoints', meas(51,:), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.versicolor, -s.Shapley.virginica, 1e-12);

%!test  # a classifier's values sum to the deviation of its scores
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), ...
%!              'NumObservationsToSample', 'all');
%! [~, sc] = predict (Mdl, meas(1,:));
%! v = [s.Shapley.setosa, s.Shapley.versicolor, s.Shapley.virginica];
%! assert_equal (sum (v, 1), sc - s.Intercept, 1e-12);

## MATLAB parity: the values of a classification tree, measured on R2024a and
## confirmed on R2026a with every observation averaged over
%!test
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), ...
%!              'NumObservationsToSample', 'all');
%! v = [s.Shapley.setosa, s.Shapley.versicolor, s.Shapley.virginica];
%! assert_equal (v, [0, 0, 0; 0, 0, 0; ...
%!                   0.666666666666667, -0.397777777777778, ...
%!                   -0.268888888888889; ...
%!                   0, 0.0644444444444445, -0.0644444444444444], 1e-12);

## MATLAB parity: a nearest-neighbour classifier, which is the model MATLAB
## answers for with the kernel algorithm this computes with
%!test
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), ...
%!              'NumObservationsToSample', 'all');
%! v = [s.Shapley.setosa, s.Shapley.versicolor, s.Shapley.virginica];
%! assert_equal (v, [0, -0.000555555555555559, 0.000555555555555538; ...
%!                   0.00111111111111099, 0.0283333333333334, ...
%!                   -0.0294444444444444; ...
%!                   0.664444444444445, -0.433888888888889, ...
%!                   -0.230555555555556; ...
%!                   0.00111111111111106, 0.0727777777777778, ...
%!                   -0.0738888888888889], 1e-12);

## The budget is capped at every subset, and a model wider than ten
## predictors takes the default of 1024 rather than being refused
%!test
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! s = shapley (@(Z) Z(:,1), X, 'MaxNumSubsets', 100);
%! assert_equal (s.NumSubsets, 4);

%!test  # a wide model takes MATLAB's default budget
%! X = repmat ((1:20)', 1, 16);
%! s = shapley (@(Z) Z(:,1), X);
%! assert_equal (s.NumSubsets, 1024);

%!test  # the weighted least squares reproduces the exact values on a linear
%!      # function, which it fits with no residual whatever subsets it is given
%! X = [1, 10, 2; 2, 20, 5; 3, 30, 1; 4, 45, 7];
%! b = [2, -3, 5];
%! f = @(Z) Z * b';
%! q = [3, 20, 4];
%! s = shapley (f, X, 'QueryPoints', q, 'MaxNumSubsets', 6, ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.Value, (b .* (q - mean (X)))', 1e-10);

%!test  # and a budget short of every subset still sums to the deviation
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), 'MaxNumSubsets', 10, ...
%!              'NumObservationsToSample', 'all');
%! [~, sc] = predict (Mdl, meas(1,:));
%! v = [s.Shapley.setosa, s.Shapley.versicolor, s.Shapley.virginica];
%! assert_equal (sum (v, 1), sc - s.Intercept, 1e-10);

%!test  # a budget of every subset agrees with the exact enumeration
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.setosa(3), 0.664444444444445, 1e-12);

## MATLAB parity: a budget short of every subset, where the four predictors
## are covered but the pairs of them are not and are taken with their
## complements.  Measured on R2024a
%!test
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), 'MaxNumSubsets', 12, ...
%!              'NumObservationsToSample', 'all');
%! v = [s.Shapley.setosa, s.Shapley.versicolor, s.Shapley.virginica];
%! assert_equal (v, [0.000555555555556, -0.007777777777778, ...
%!                   0.007222222222222; ...
%!                   0.000555555555556, 0.035555555555556, ...
%!                   -0.036111111111111; ...
%!                   0.664444444444444, -0.435555555555555, ...
%!                   -0.228888888888889; ...
%!                   0.001111111111111, 0.074444444444445, ...
%!                   -0.075555555555555], 1e-12);

## MATLAB parity: five predictors, where the budget runs out inside the
## second pair of cardinalities.  Measured on R2024a
%!test
%! load fisheriris
%! X5 = [meas, meas(:,1) .* meas(:,4)];
%! Mdl = fitcknn (X5, species);
%! s = shapley (Mdl, 'QueryPoints', X5(1,:), 'MaxNumSubsets', 24, ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.setosa', ...
%!               [0.000625, 0.002239583333333, 0.003510416666667, ...
%!                0.002843750000000, 0.657447916666667], 1e-12);

## MATLAB parity: conditioning on a predictor keeps the tenth of the
## observations nearest the query point in it.  A hundred and one
## observations is what fixes that tenth as rounded up rather than to
## nearest.  Measured on R2024a
%!test
%! t = (1:101)';
%! X = [t, mod(t * 37, 101), mod(t * 53, 101)];
%! s = shapley (@(Z) Z(:,1), X, 'QueryPoints', [10, 50, 90], ...
%!              'Method', 'conditional', 'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.Value', [-45.560606060606048, ...
%!               1.621212121212129, 2.939393939393916], 1e-10);

## MATLAB parity: a classifier answered for by conditioning.  Measured on
## R2024a
%!test
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), 'Method', 'conditional', ...
%!              'NumObservationsToSample', 'all');
%! v = [s.Shapley.setosa, s.Shapley.versicolor, s.Shapley.virginica];
%! assert_equal (v, [0.15, -0.0666666666666667, -0.0833333333333334; ...
%!                   0.172222222222222, -0.0888888888888889, ...
%!                   -0.0833333333333333; ...
%!                   0.172222222222222, -0.0888888888888889, ...
%!                   -0.0833333333333333; ...
%!                   0.172222222222222, -0.0888888888888889, ...
%!                   -0.0833333333333334], 1e-12);

%!test  # conditioning is reported as its own algorithm
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! s = shapley (@(Z) Z(:,1), X, 'Method', 'conditional');
%! assert_equal (s.Method, 'conditional-kernel');

%!test  # and still sums to the deviation of the prediction from the average
%! t = (1:100)';
%! X = [t, mod(t * 37, 101), mod(t * 53, 101)];
%! f = @(Z) Z(:,1);
%! s = shapley (f, X, 'QueryPoints', [10, 50, 90], ...
%!              'Method', 'conditional', 'NumObservationsToSample', 'all');
%! assert_equal (sum (s.Shapley.Value), f ([10, 50, 90]) - s.Intercept, 1e-10);

%!test  # a predictor recorded on a wider scale does not pull the neighbours
%! t = (1:100)';
%! X = [t, mod(t * 37, 101), mod(t * 53, 101)];
%! W = X;
%! W(:,2) = W(:,2) * 1000;
%! s = shapley (@(Z) Z(:,1), X, 'QueryPoints', [10, 50, 90], ...
%!              'Method', 'conditional', 'NumObservationsToSample', 'all');
%! r = shapley (@(Z) Z(:,1), W, 'QueryPoints', [10, 50000, 90], ...
%!              'Method', 'conditional', 'NumObservationsToSample', 'all');
%! assert_equal (r.Shapley.Value, s.Shapley.Value, 1e-10);

%!test  # nothing but the two extremes says nothing about any one predictor
%! X = [1, 10; 2, 20; 3, 30; 4, 45];
%! s = shapley (@(Z) Z(:,1), X, 'QueryPoints', [3, 20], ...
%!              'MaxNumSubsets', 2, 'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.Value, [NaN; NaN]);

%!test  # a model predicting a weighted sum is answered from the weights
%! load fisheriris
%! Mdl = fitrlinear (meas(:,2:4), meas(:,1));
%! s = shapley (Mdl, meas(:,2:4), 'QueryPoints', meas(1,2:4), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-linear');

%!test  # and gives what the subsets would have given
%! load fisheriris
%! Mdl = fitrlinear (meas(:,2:4), meas(:,1));
%! s = shapley (Mdl, meas(:,2:4), 'QueryPoints', meas(1,2:4), ...
%!              'NumObservationsToSample', 'all');
%! k = shapley (Mdl, meas(:,2:4), 'QueryPoints', meas(1,2:4), ...
%!              'MaxNumSubsets', 8, 'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.Value, k.Shapley.Value, 1e-10);

%!test  # which is each weight times the deviation of its predictor
%! load fisheriris
%! Mdl = fitrlinear (meas(:,2:4), meas(:,1));
%! s = shapley (Mdl, meas(:,2:4), 'QueryPoints', meas(1,2:4), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.Value, ...
%!               Mdl.Beta(:) .* (meas(1,2:4) - mean (meas(:,2:4)))', 1e-10);

%!test  # a support vector machine on any other kernel keeps no such weights
%! load fisheriris
%! Mdl = fitrsvm (meas(:,2:4), meas(:,1), 'KernelFunction', 'gaussian');
%! s = shapley (Mdl, meas(:,2:4), 'QueryPoints', meas(1,2:4), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-kernel');

## A transform of the score breaks the weighted sum, so the subsets answer
## and the values still sum to the deviation.  MATLAB keeps the short answer
## here and reports an intercept on the transformed scale beside values on
## the untransformed one, which do not sum to it
%!test
%! load fisheriris
%! Mdl = fitcsvm (meas(51:150,:), species(51:150), ...
%!                'ScoreTransform', 'logit');
%! s = shapley (Mdl, meas(51:150,:), 'QueryPoints', meas(51,:), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-kernel');
%! [~, sc] = predict (Mdl, meas(51,:));
%! v = [s.Shapley.versicolor, s.Shapley.virginica];
%! assert_equal (sum (v, 1), sc - s.Intercept, 1e-10);

%!test  # asking for a budget of subsets asks for the subsets
%! load fisheriris
%! Mdl = fitrlinear (meas(:,2:4), meas(:,1));
%! s = shapley (Mdl, meas(:,2:4), 'QueryPoints', meas(1,2:4), ...
%!              'MaxNumSubsets', 8, 'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-kernel');

%!test  # a decision tree is answered leaf by leaf
%! load fisheriris
%! Mdl = fitrtree (meas(:,2:4), meas(:,1));
%! s = shapley (Mdl, 'QueryPoints', meas(1,2:4), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-tree');

%!test  # and gives what the subsets would have given
%! load fisheriris
%! Mdl = fitrtree (meas(:,2:4), meas(:,1));
%! s = shapley (Mdl, 'QueryPoints', meas(1,2:4), ...
%!              'NumObservationsToSample', 'all');
%! k = shapley (Mdl, meas(:,2:4), 'QueryPoints', meas(1,2:4), ...
%!              'MaxNumSubsets', 8, 'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.Value, k.Shapley.Value, 1e-10);

%!test  # for a classifier too, one column of values per class
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), ...
%!              'NumObservationsToSample', 'all');
%! k = shapley (Mdl, meas, 'QueryPoints', meas(1,:), ...
%!              'MaxNumSubsets', 16, 'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.versicolor, k.Shapley.versicolor, 1e-10);

%!test  # and for a compact one, which keeps the nodes and not the data
%! load fisheriris
%! Mdl = compact (fitctree (meas, species));
%! s = shapley (Mdl, meas, 'QueryPoints', meas(1,:), ...
%!              'NumObservationsToSample', 'all');
%! k = shapley (Mdl, meas, 'QueryPoints', meas(1,:), ...
%!              'MaxNumSubsets', 16, 'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.setosa, k.Shapley.setosa, 1e-10);

%!test  # a cut on the levels of a categorical predictor is followed too
%! t = (1:120)';
%! X = [mod(t, 4), double(t), mod(t * 7, 11)];
%! y = 3 * (X(:,1) == 1) - 2 * (X(:,1) == 3) + 0.01 * mod (t, 5);
%! Mdl = fitrtree (X, y, 'CategoricalPredictors', 1);
%! s = shapley (Mdl, 'QueryPoints', X(3,:), ...
%!              'NumObservationsToSample', 'all');
%! k = shapley (Mdl, X, 'QueryPoints', X(3,:), 'MaxNumSubsets', 8, ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Shapley.Value, k.Shapley.Value, 1e-10);

%!test  # a tree is exact however many predictors it has, where the subsets
%!      # would have stopped at the budget of 1024
%! rand ('seed', 7);
%! randn ('seed', 7);
%! X = randn (200, 14);
%! y = X(:,1) + 2 * X(:,5) - X(:,9);
%! Mdl = fitrtree (X, y);
%! s = shapley (Mdl, 'QueryPoints', X(1,:), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-tree');
%! assert_equal (sum (s.Shapley.Value), ...
%!               predict (Mdl, X(1,:)) - s.Intercept, 1e-10);

## A transform of the score is applied to each leaf, and exactly one leaf
## answers for a row, so the tree still comes apart and the values still sum
## to the deviation.  MATLAB reports an intercept on the transformed scale
## beside values on the untransformed one, which do not sum to it
%!test
%! load fisheriris
%! Mdl = fitctree (meas(51:150,:), species(51:150), ...
%!                 'ScoreTransform', 'logit');
%! s = shapley (Mdl, 'QueryPoints', meas(51,:), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-tree');
%! [~, sc] = predict (Mdl, meas(51,:));
%! v = [s.Shapley.versicolor, s.Shapley.virginica];
%! assert_equal (sum (v, 1), sc - s.Intercept, 1e-10);

%!test  # asking for a budget of subsets asks for the subsets
%! load fisheriris
%! Mdl = fitrtree (meas(:,2:4), meas(:,1));
%! s = shapley (Mdl, meas(:,2:4), 'QueryPoints', meas(1,2:4), ...
%!              'MaxNumSubsets', 8, 'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-kernel');

%!test  # an ensemble of trees is walked leaf by leaf, as one tree is
%! rand ('seed', 11);
%! randn ('seed', 11);
%! X = randn (60, 11);
%! y = X(:,1) + 2 * X(:,2) - X(:,3);
%! Mdl = fitrensemble (X, y, 'Method', 'Bag', 'NumLearningCycles', 4);
%! s = shapley (Mdl, 'QueryPoints', X(1,:), ...
%!              'NumObservationsToSample', 'all');
%! k = shapley (Mdl, X, 'QueryPoints', X(1,:), 'MaxNumSubsets', 2048, ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-tree');
%! assert_equal (s.Shapley.Value, k.Shapley.Value, 1e-10);

%!test  # and a classifier's, where a learner's classes are put where the
%!      # ensemble keeps them
%! rand ('seed', 11);
%! randn ('seed', 11);
%! X = randn (60, 11);
%! g = repmat ({'a'; 'b'}, 30, 1);
%! Mdl = fitcensemble (X, g, 'Method', 'Bag', 'NumLearningCycles', 4);
%! s = shapley (Mdl, 'QueryPoints', X(1,:), ...
%!              'NumObservationsToSample', 'all');
%! k = shapley (Mdl, X, 'QueryPoints', X(1,:), 'MaxNumSubsets', 2048, ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-tree');
%! assert_equal (s.Shapley.a, k.Shapley.a, 1e-10);

%!test  # a narrow ensemble is walked too, and agrees with every subset
%! load fisheriris
%! Mdl = fitcensemble (meas, species, 'Method', 'Bag', ...
%!                     'NumLearningCycles', 10);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), ...
%!              'NumObservationsToSample', 'all');
%! k = shapley (Mdl, meas, 'QueryPoints', meas(1,:), 'MaxNumSubsets', 16, ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-tree');
%! assert_equal (s.Shapley.setosa, k.Shapley.setosa, 1e-10);

%!test  # a boosting method reading a learner's label, not its scores, too
%! rand ('seed', 3);
%! randn ('seed', 3);
%! X = randn (60, 11);
%! g = repmat ({'a'; 'b'}, 30, 1);
%! Mdl = fitcensemble (X, g, 'Method', 'AdaBoostM1', ...
%!                     'NumLearningCycles', 4);
%! s = shapley (Mdl, 'QueryPoints', X(1,:), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.Method, 'interventional-kernel');

%!test  # a classifier's fitted label keeps the type of the response
%! load fisheriris
%! Mdl = fitctree (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (s.BlackboxFitted, predict (Mdl, meas(1,:)));

%!error<shapley: too few input arguments.> shapley ()
%!error<shapley: BLACKBOX must be a fitted model that predicts, or a function handle.> shapley (42)
%!error<shapley: X is required when the model is a function handle.> shapley (@(Z) Z(:,1))
%!error<shapley: X must be a real numeric matrix.> shapley (@(Z) Z(:,1), {1, 2})
%!error<shapley: 'UseParallel' is not implemented.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'UseParallel', true)
%!error<shapley: 'MaxNumSubsets' must be an integer greater than 1.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'MaxNumSubsets', 1)
%!error<shapley: 'MaxNumSubsets' must be an integer greater than 1.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'MaxNumSubsets', 2.5)
%!error<shapley: 'Method' must be 'interventional' or 'conditional'.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'Method', 'marginal')
%!error<shapley: 'NumObservationsToSample' must be a positive integer or 'all'.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'NumObservationsToSample', 0)
%!error<shapley: unknown optional argument or misplaced value.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'NoSuchThing', 1)
%!error<shapley.fit: QUERYPOINTS must be a real numeric matrix.> fit (shapley (@(Z) Z(:,1), [1, 2; 3, 4]), 'abc')
%!error<shapley.fit: QUERYPOINTS must have one column per predictor of the model.> fit (shapley (@(Z) Z(:,1), [1, 2; 3, 4]), [1, 2, 3])

## The four drawing methods, measured on R2024a 2026-09-22
%!test  # one query point: one bar per predictor, the least important first
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = plot (s);
%!   assert_equal (numel (b), 1);
%!   assert_equal (get (b, 'horizontal'), 'on');
%!   assert_equal (get (b, 'ydata')', [-2.8333333333333333, -5, 80], 1e-12);
%!   assert_equal (get (gca (), 'yticklabel')', {'x3', 'x1', 'x2'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # one query point is titled as an explanation of that point
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   plot (s);
%!   ax = gca ();
%!   assert_equal (get (get (ax, 'title'), 'string'), 'Shapley Explanation');
%!   assert_equal (get (get (ax, 'xlabel'), 'string'), 'Shapley Value');
%!   assert_equal (get (get (ax, 'ylabel'), 'string'), 'Predictor');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Several query points are summarised by the mean of the absolute values
%!test
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = plot (s);
%!   ax = gca ();
%!   assert_equal (get (get (ax, 'title'), 'string'), ...
%!                 'Shapley Importance Plot');
%!   assert_equal (get (get (ax, 'xlabel'), 'string'), ...
%!                 'Mean of Absolute Shapley Values');
%!   assert_equal (get (b, 'ydata')', [2.5, 2.75, 43.75], 1e-12);
%!   assert_equal (get (ax, 'yticklabel')', {'x1', 'x3', 'x2'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Over several query points whatever is left over is summed into one bar
%!test
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = plot (s, 'NumImportantPredictors', 2);
%!   assert_equal (get (b, 'ydata')', [2.5, 2.75, 43.75], 1e-12);
%!   assert_equal (get (gca (), 'yticklabel')', ...
%!                 {'Sum of other 1 predictor(s)', 'x3', 'x2'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # over one query point it is left out rather than summed
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = plot (s, 'NumImportantPredictors', 2);
%!   assert_equal (get (b, 'ydata')', [-5, 80], 1e-12);
%!   assert_equal (get (gca (), 'yticklabel')', {'x1', 'x2'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # naming one query point explains that point rather than summarising
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = plot (s, 'QueryPointIndices', 2);
%!   ax = gca ();
%!   assert_equal (get (get (ax, 'title'), 'string'), 'Shapley Explanation');
%!   assert_equal (get (b, 'ydata')', [2.1666666666666667, -3, 50], 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the axes to draw into may be given first, as for any plot
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   a1 = subplot (1, 2, 1);
%!   a2 = subplot (1, 2, 2);
%!   b = plot (a2, s);
%!   assert_equal (get (b, 'parent'), a2);
%!   assert_equal (isempty (get (a1, 'children')), true);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## MATLAB parity: one query point of a classifier explains the class predicted
%!test
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas(1,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = plot (s);
%!   assert_equal (numel (b), 1);
%!   assert_equal (s.BlackboxFitted, {'setosa'});
%!   assert_equal (sort (get (b, 'ydata')'), sort (s.Shapley.setosa'), 1e-12);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## MATLAB parity: several query points of a classifier draw every class
%!test
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas([1, 60, 120],:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = plot (s);
%!   assert_equal (numel (b), 3);
%!   assert_equal (get (b(1), 'displayname'), 'setosa');
%!   assert_equal (get (b(3), 'displayname'), 'virginica');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the classes to draw may be named
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas([1, 60, 120],:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = plot (s, 'ClassNames', {'setosa', 'virginica'});
%!   assert_equal (numel (b), 2);
%!   assert_equal (get (b(2), 'displayname'), 'virginica');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # boxchart draws one box per predictor, lying horizontally
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (s);
%!   assert_equal (class (b), 'stats.chart.BoxChart');
%!   assert_equal (b.Orientation, 'horizontal');
%!   assert_equal (class (b.XData), 'categorical');
%!   assert_equal (numel (findobj (gca (), 'type', 'patch')), 3);
%!   assert_equal (get (gca (), 'yticklabel')', {'x1', 'x3', 'x2'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # boxchart is titled as a summary over the query points
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   boxchart (s);
%!   ax = gca ();
%!   assert_equal (get (get (ax, 'title'), 'string'), 'Shapley Summary Plot');
%!   assert_equal (get (get (ax, 'xlabel'), 'string'), 'Shapley Value');
%!   assert_equal (get (get (ax, 'ylabel'), 'string'), 'Predictor');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## MATLAB parity: boxchart leaves the rest out rather than summing them
%!test
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   boxchart (s, 'NumImportantPredictors', 2);
%!   assert_equal (get (gca (), 'yticklabel')', {'x3', 'x2'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## MATLAB parity: a chart of one class orders the predictors by that class
%!test
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! s = shapley (Mdl, 'QueryPoints', meas([1, 60, 120],:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   boxchart (s, 'ClassName', 'virginica');
%!   [~, ord] = sort (s.MeanAbsoluteShapley.virginica, 'ascend');
%!   assert_equal (get (gca (), 'yticklabel')', ...
%!                 cellstr (s.Shapley.Predictor(ord))');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the outlier markers may be spread across the width of the box
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   b = boxchart (s, 'JitterOutliers', 'on');
%!   assert_equal (b.JitterOutliers, 'on');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # swarmchart draws one row of points per predictor
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = swarmchart (s);
%!   assert_equal (numel (h), 3);
%!   assert_equal (get (h(1), 'xdata')', [-5, -3, -1, 1], 1e-12);
%!   assert_equal (get (gca (), 'yticklabel')', {'x1', 'x3', 'x2'});
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## MATLAB parity: a point is coloured by where the predictor's value
## stands in its own range over the query points
%!test
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = swarmchart (s);
%!   assert_equal (get (h(1), 'cdata')', [0, 1/3, 2/3, 1], 1e-12);
%!   assert_equal (get (gca (), 'clim'), [0, 1]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # the colour map the values are read through may be given
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   swarmchart (s, 'ColorMap', [1, 0, 0; 0, 1, 0]);
%!   assert_equal (get (gca (), 'colormap'), [1, 0, 0; 0, 1, 0]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # how the points of a row are spread may be chosen
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   h = swarmchart (s, 'YJitter', 'none');
%!   assert_equal (get (h(1), 'ydata')', [1, 1, 1, 1]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # plotDependence draws one predictor against its own values
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   p = plotDependence (s, 'x1');
%!   assert_equal (get (p, 'xdata')', [1, 2, 3, 4]);
%!   assert_equal (get (p, 'ydata')', [-5, -3, -1, 1], 1e-12);
%!   ax = gca ();
%!   assert_equal (get (get (ax, 'title'), 'string'), ...
%!                 'Shapley Dependence Plot');
%!   assert_equal (get (get (ax, 'xlabel'), 'string'), 'x1');
%!   assert_equal (get (get (ax, 'ylabel'), 'string'), ...
%!                 'Shapley Values for x1');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a predictor may be indexed rather than named
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   p = plotDependence (s, 2);
%!   assert_equal (get (p, 'xdata')', [10, 20, 30, 45]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## MATLAB parity: the colouring predictor is read as it stands, the range
## of the axes carrying the scale, where a swarm normalizes instead
%!test
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   p = plotDependence (s, 1, 'ColorPredictor', 'x2');
%!   assert_equal (get (p, 'cdata')', [10, 20, 30, 45]);
%!   assert_equal (get (gca (), 'clim'), [10, 45]);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a colour bar is drawn beside a chart whose points are coloured
%! X = [1, 10, 100; 2, 20, 150; 3, 30, 120; 4, 45, 180; 5, 50, 90; ...
%!      6, 65, 130];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2) + 0.1 * Z(:,3);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   plotDependence (s, 1);
%!   none = numel (findobj (hf, 'tag', 'colorbar'));
%!   plotDependence (s, 1, 'ColorPredictor', 'x2');
%!   assert_equal (none, 0);
%!   assert_equal (numel (findobj (hf, 'tag', 'colorbar')), 1);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

%!test  # a predictor holding levels is drawn as a box over each level
%! X = [1, 0; 2, 0; 3, 1; 4, 1; 5, 0; 6, 1];
%! f = @(Z) 2 * Z(:,1) - 3 * Z(:,2);
%! s = shapley (f, X, 'QueryPoints', X(1:4,:), ...
%!              'CategoricalPredictors', 2, ...
%!              'NumObservationsToSample', 'all');
%! hf = figure ('visible', 'off');
%! unwind_protect
%!   p = plotDependence (s, 2);
%!   assert_equal (class (p), 'stats.chart.BoxChart');
%!   assert_equal (numel (findobj (gca (), 'type', 'patch')), 2);
%!   assert_equal (get (get (gca (), 'title'), 'string'), ...
%!                 'Shapley Dependence Plot');
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Input validation
%!error<shapley.plot: the Shapley values are not fitted; use fit to compute them.>
%! plot (shapley (@(Z) Z(:,1), [1, 2; 3, 4]))

%!error<shapley.boxchart: the Shapley values are not fitted; use fit to compute them.>
%! boxchart (shapley (@(Z) Z(:,1), [1, 2; 3, 4]))

%!error<shapley.swarmchart: the Shapley values are not fitted; use fit to compute them.>
%! swarmchart (shapley (@(Z) Z(:,1), [1, 2; 3, 4]))

%!error<shapley.plotDependence: the Shapley values are not fitted; use fit to compute them.>
%! plotDependence (shapley (@(Z) Z(:,1), [1, 2; 3, 4]), 1)

%!error<shapley.plot: 'NumImportantPredictors' must be a positive integer.>
%! plot (shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'QueryPoints', [1, 2]), ...
%!       'NumImportantPredictors', 0)

%!error<shapley.plot: 'QueryPointIndices' must be positive integers.>
%! plot (shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'QueryPoints', [1, 2]), ...
%!       'QueryPointIndices', 1.5)

%!error<shapley.plot: 'QueryPointIndices' must not exceed the 1 query points fitted.>
%! plot (shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'QueryPoints', [1, 2]), ...
%!       'QueryPointIndices', 2)

%!error<shapley.plot: 'ClassNames' is valid only for a classification model.>
%! plot (shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'QueryPoints', [1, 2]), ...
%!       'ClassNames', 'setosa')

%!error<shapley.boxchart: 'ClassName' is valid only for a classification model.>
%! boxchart (shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'QueryPoints', [1, 2]), ...
%!           'ClassName', 'setosa')

%!error<shapley.boxchart: 'JitterOutliers' must be one of 'on', 'off'.>
%! boxchart (shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'QueryPoints', [1, 2]), ...
%!           'JitterOutliers', 'maybe')

%!error<shapley.swarmchart: 'YJitter' must be one of 'none', 'density', 'rand', 'randn'.>
%! swarmchart (shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'QueryPoints', [1, 2]), ...
%!             'YJitter', 'sideways')

%!error<shapley.swarmchart: 'ColorMap' must be a colour map name or a matrix of RGB triplets.>
%! swarmchart (shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'QueryPoints', [1, 2]), ...
%!             'ColorMap', 5)

%!error<shapley.plotDependence: too few input arguments.>
%! plotDependence (shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'QueryPoints', [1, 2]))

%!error<shapley.plotDependence: PREDICTOR must name or index one predictor of the model.>
%! plotDependence (shapley (@(Z) Z(:,1), [1, 2; 3, 4], ...
%!                          'QueryPoints', [1, 2]), 'nope')

%!error<shapley.plot: unknown optional argument or misplaced value.>
%! plot (shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'QueryPoints', [1, 2]), ...
%!       'NoSuchThing', 1)

## Table input
%!shared stT, stM
%! load fisheriris
%! stT = table (meas(:,2), meas(:,3), meas(:,4), meas(:,1), ...
%!              'VariableNames', {'SW', 'PL', 'PW', 'SL'});
%! stT.Wide = categorical (meas(:,2) > 3, [false true], {'narrow', 'wide'});
%! stM = fitrtree (stT, 'SL');

%!test  # the observations to average over may be given as a table
%! s = shapley (stM, stT(:,[1, 2, 3, 5]), 'QueryPoints', ...
%!              stT(1,[1, 2, 3, 5]), 'NumObservationsToSample', 'all');
%! assert_equal (size (s.X), [150, 4]);
%! assert_equal (cellstr (s.Shapley.Predictor)', ...
%!               {'SW', 'PL', 'PW', 'Wide'});

## MATLAB parity: a table is read by the names the model was fitted on, so
## the order of its columns does not matter
%!test
%! a = shapley (stM, stT(:,[1, 2, 3, 5]), 'QueryPoints', ...
%!              stT(1,[1, 2, 3, 5]), 'NumObservationsToSample', 'all');
%! b = shapley (stM, stT(:,[5, 3, 2, 1]), 'QueryPoints', ...
%!              stT(1,[5, 3, 2, 1]), 'NumObservationsToSample', 'all');
%! assert_equal (b.Shapley.Value, a.Shapley.Value);

## MATLAB parity: a column the model was not fitted on is passed over
%!test
%! a = shapley (stM, stT(:,[1, 2, 3, 5]), 'QueryPoints', ...
%!              stT(1,[1, 2, 3, 5]), 'NumObservationsToSample', 'all');
%! b = shapley (stM, stT, 'QueryPoints', stT(1,:), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (b.Shapley.Value, a.Shapley.Value);

%!test  # fit takes a table too, read the same way
%! a = shapley (stM, stT(:,[1, 2, 3, 5]), 'QueryPoints', ...
%!              stT(1,[1, 2, 3, 5]), 'NumObservationsToSample', 'all');
%! b = shapley (stM, stT(:,[1, 2, 3, 5]), ...
%!              'NumObservationsToSample', 'all');
%! b = fit (b, stT(1,[5, 3, 2, 1]));
%! assert_equal (b.Shapley.Value, a.Shapley.Value);

%!test  # a level carries the code it carried at fitting
%! T = stT(1:60,:);
%! T.Wide = categorical (repmat ({'wide'}, 60, 1), {'narrow', 'wide'});
%! s = shapley (stM, T(:,[1, 2, 3, 5]), 'QueryPoints', ...
%!              T(1,[1, 2, 3, 5]), 'NumObservationsToSample', 'all');
%! assert_equal (rows (s.Shapley), 4);

## A function handle names no predictors, so a table given for one names
## them itself and its columns are taken in the order they come
%!test
%! f = @(Z) Z(:,1);
%! s = shapley (f, stT(:,1:3), 'QueryPoints', stT(1,1:3), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (cellstr (s.Shapley.Predictor)', {'SW', 'PL', 'PW'});
%! r = shapley (f, stT(:,[3, 2, 1]), 'QueryPoints', stT(1,[3, 2, 1]), ...
%!              'NumObservationsToSample', 'all');
%! assert_equal (cellstr (r.Shapley.Predictor)', {'PW', 'PL', 'SW'});

%!error<shapley: the table holds no predictor 'PL'.> ...
%! shapley (stM, stT(:,[1, 3, 5]), 'QueryPoints', stT(1,[1, 3, 5]))

%!error<shapley.fit: the table holds no predictor 'PL'.> ...
%! fit (shapley (stM, stT(:,[1, 2, 3, 5]), ...
%!               'NumObservationsToSample', 'all'), stT(1,[1, 3, 5]))
