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
    ## weights alone; @qcode{'interventional-tree'} for a decision tree,
    ## which is answered leaf by leaf; @qcode{'interventional-kernel'} for
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
    ## @seealso{shapley, fit}
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
    ## and one column per predictor.  The values already held are replaced,
    ## not added to.
    ##
    ## @seealso{shapley}
    ## @end deftypefn

    function this = fit (this, QueryPoints)

      ## Input validation
      if (nargin != 2)
        error ("shapley.fit: invalid number of input arguments.");
      endif
      M = numel (this.PredictorNames);
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
      if (tree)
        [Tree, leafval] = shapTreeOf (this.BlackboxModel);
        phi1 = shapTreeValues (Tree, leafval, Xs, QueryPoints(1,:), M, K);
        tree = shapAccounts (sum (phi1, 1), sfcn, Xs, QueryPoints(1,:));
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
      for ii = 1:nq
        if (tree)
          phi(:,:,ii) = shapTreeValues (Tree, leafval, Xs, ...
                                        QueryPoints(ii,:), M, K);
          if (ii == 1)
            icept = mean (sfcn (Xs), 1);
          endif
          continue;
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

  endmethods

endclassdef

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
  else
    method = 'interventional-tree';
  endif

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
  F.IsTree = ! isempty (shapTreeOf (blackbox));
  F.Cat = cat;
  F.X = Data;
  F.Idx = idx;

endfunction

## The decision tree a model answers through, and what it answers at each of
## its leaves, or empty where the model is not one tree.  A tree is answered
## exactly without enumerating the subsets, so this is worth asking.
function [T, leafval] = shapTreeOf (blackbox)

  T = [];
  leafval = [];
  if (is_function_handle (blackbox))
    return;
  endif
  switch (class (blackbox))
    case {'ClassificationTree', 'CompactClassificationTree', ...
          'RegressionTree', 'CompactRegressionTree'}
      T = blackbox;
    otherwise
      return;
  endswitch

  props = properties (T);
  if (any (strcmp (props, 'ClassProbability')))
    leafval = T.ClassProbability;
    if (! isempty (T.STfun))
      leafval = T.STfun (leafval);
    endif
  else
    leafval = T.NodeMean(:);
    if (! isempty (T.RTfun))
      leafval = T.RTfun (leafval);
    endif
  endif

endfunction

## Every leaf of a tree, with the nodes walked to reach it and the child each
## of them was entered by.  A leaf's answer holds for exactly the rows that
## satisfy every one of those, which is what lets the tree be taken apart.
function [lnode, lpath] = shapLeafPaths (T)

  lnode = [];
  lpath = {};
  root = zeros (0, 2);
  stack = {{1, root}};
  while (! isempty (stack))
    cur = stack{end};
    stack(end) = [];
    n = cur{1};
    pth = cur{2};
    if (T.CutPredictorIndex(n) == 0)
      lnode(end+1) = n;
      lpath{end+1} = pth;
      continue;
    endif
    for kid = 1:2
      stack{end+1} = {T.Children(n,kid), [pth; n, kid]};
    endfor
  endwhile

endfunction

## Which child each value of the cut predictor enters at node N, zero where it
## enters neither, which is how a missing value and an unknown level behave.
function c = shapChildFor (T, n, vals)

  L = T.Children(n,1);
  R = T.Children(n,2);
  c = zeros (numel (vals), 1);
  vals = vals(:);
  if (! isempty (T.CutCategories) && ! isempty (T.CutCategories{n,1}))
    c(ismember (vals, T.CutCategories{n,1})) = L;
    c(ismember (vals, T.CutCategories{n,2})) = R;
  else
    ok = ! isnan (vals);
    c(ok & vals < T.CutPoint(n)) = L;
    c(ok & vals >= T.CutPoint(n)) = R;
  endif

endfunction

## The sum, over every way of adding T of the M-A-B predictors that the leaf
## does not constrain, of the weight the Shapley definition gives a subset of
## that size.  Taken in logs, so that neither the binomial coefficient nor the
## factorials overflow however wide the model is.
function g = shapKernelSum (k, m, M)

  t = (0:m)';
  sz = k + t;
  lt = gammaln (m + 1) - gammaln (t + 1) - gammaln (m - t + 1) ...
       + gammaln (sz + 1) + gammaln (M - sz) - gammaln (M + 1);
  g = sum (exp (lt));

endfunction

## The Shapley values of a decision tree, taken leaf by leaf rather than
## subset by subset.
##
## A leaf answers for the rows meeting every condition on the way to it.  Hold
## the predictors of a subset at the query point and leave the rest as an
## observation has them, and that row reaches the leaf exactly when every
## predictor the query fails is outside the subset and every predictor the
## observation fails is inside it.  So, writing A for the predictors the query
## passes and the observation fails and B for the other way about, the leaf is
## reached exactly for the subsets holding all of A and none of B.  Only the
## predictors in A and B can move that, and by how much depends on nothing but
## how many there are, which collapses the sum over subsets into one term.
function phi = shapTreeValues (T, leafval, Xs, q, M, K)

  n = rows (Xs);
  phi = zeros (M, K);
  [lnode, lpath] = shapLeafPaths (T);

  for li = 1:numel (lnode)
    P = lpath{li};
    if (isempty (P))
      continue;
    endif
    nodes = P(:,1);
    kids = P(:,2);
    feats = T.CutPredictorIndex(nodes);
    uf = unique (feats(:))';
    nf = numel (uf);

    ## Which of those predictors the query passes, and which each observation
    ## passes, taking every condition on the way as one demand
    passX = true (1, nf);
    passZ = true (n, nf);
    for jj = 1:nf
      cc = find (feats == uf(jj));
      for c = cc(:)'
        want = T.Children(nodes(c), kids(c));
        passX(jj) = passX(jj) ...
                    && (shapChildFor (T, nodes(c), q(uf(jj))) == want);
        passZ(:,jj) = passZ(:,jj) ...
                      & (shapChildFor (T, nodes(c), Xs(:,uf(jj))) == want);
      endfor
    endfor

    ## An observation failing a predictor the query fails too puts the leaf
    ## out of reach whatever the subset is
    inB = find (! passX);
    inA = find (passX);
    if (isempty (inB))
      keep = true (n, 1);
    else
      keep = all (passZ(:,inB), 2);
    endif
    if (! any (keep))
      continue;
    endif
    b = numel (inB);
    passA = passZ(keep,inA);
    aCount = numel (inA) - sum (passA, 2);
    vL = leafval(lnode(li),:);
    fB = uf(inB);
    fA = uf(inA);

    ## The weight depends on an observation only through how many predictors
    ## fall in A, so it is worked out once for each count that occurs
    for av = unique (aCount)'
      atRow = (aCount == av);
      mv = M - av - b;
      if (b > 0)
        gB = shapKernelSum (av, mv, M);
        phi(fB,:) -= sum (atRow) * gB * repmat (vL, b, 1);
      endif
      if (av >= 1)
        gA = shapKernelSum (av - 1, mv, M);
        cnt = sum (! passA(atRow,:), 1);
        phi(fA,:) += gA * (cnt(:) * vL);
      endif
    endfor
  endfor

  phi = phi / n;

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
