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
## @deftypefn {Private Function} {[@var{pd}, @var{x}, @var{y}, @var{opts}] =} pdCompute (@var{fname}, @var{Mdl}, @var{Vars}, @dots{})
##
## Compute partial dependence, for whichever entry point was called.
##
## @code{partialDependence} and the method of the same name on
## @code{PredictiveModel} both come here, so that the two cannot drift apart.
## @var{fname} names the caller in every error message, and the rest is what
## that caller was given.  @code{plotPartialDependence} comes here too, for
## the values it draws.
##
## @var{opts} holds what the drawing needs and the values do not: the fields
## @qcode{Conditional} and @qcode{Parent} are the two options only the plot
## takes, @qcode{PredArgs} the pairs given to @code{predict}, and
## @qcode{Frame} the frame the values were averaged over.
##
## @end deftypefn

function [pd, x, y, opts] = pdCompute (fname, Mdl, Vars, varargin)

  if (! (isa (Mdl, 'PredictiveModel') || is_function_handle (Mdl)))
    error (strcat ("%s: MDL must be a fitted model that", ...
                   " predicts, or a function handle."), fname);
  endif

  ## A classifier is told which classes to answer for, and the class names
  ## come before the data, as MATLAB orders them
  isclass = (! is_function_handle (Mdl)
             && any (strcmp (properties (Mdl), 'ClassNames'))
             && ! isempty (Mdl.ClassNames));
  args = varargin;
  Labels = [];
  if (isclass)
    if (isempty (args))
      error (strcat ("%s: LABELS is required for a", ...
                     " classification model."), fname);
    endif
    Labels = args{1};
    args(1) = [];
  endif
  Data = [];
  if (! isempty (args) && isnumeric (args{1}))
    Data = args{1};
    args(1) = [];
  elseif (! isempty (args)
          && ! (ischar (args{1}) || isa (args{1}, 'string')))
    ## Neither the data nor the name of an option: for a regression model the
    ## likeliest mistake is classes it has none of, and for a classifier,
    ## whose classes were taken above, it is the data
    if (isclass)
      error ("%s: DATA must be a real numeric matrix.", fname);
    else
      error (strcat ("%s: LABELS applies only to a classification", ...
                     " model."), fname);
    endif
  endif

  optNames = {'QueryPoints', 'NumObservationsToSample', ...
              'CategoricalPredictors', 'IncludeInteractions', ...
              'IncludeIntercept', 'OutputColumns', 'UseParallel', ...
              'PredictionForMissingValue', 'Conditional', 'Parent'};
  dfValues = {[], [], [], [], true, 'all', [], [], [], []};
  [QP, NumObs, CatPred, Inter, Icept, OutCols, Par, Miss, Cond, Ax, rem] = ...
                              parsePairedArguments (optNames, dfValues, args);
  if (! isempty (rem))
    error (strcat ("%s: unknown optional argument or", ...
                   " misplaced value."), fname);
  endif
  if (! isempty (Par))
    error ("%s: 'UseParallel' is not implemented.", fname);
  endif
  if (! isempty (Miss))
    error (strcat ("%s: 'PredictionForMissingValue' is not", ...
                   " implemented."), fname);
  endif
  ## Both entry points come here, and the method wins the dispatch whenever
  ## the first argument is a model, so this is the only place that sees them
  if (! strcmp (fname, 'plotPartialDependence')
      && (! isempty (Cond) || ! isempty (Ax)))
    error (strcat ("%s: 'Conditional' and 'Parent' belong to", ...
                   " plotPartialDependence."), fname);
  endif

  [predArgs, subIcept, OutCols, errmsg] = pdOptions (Mdl, Inter, Icept, ...
                                                     OutCols, CatPred);
  if (! isempty (errmsg))
    error ("%s: %s", fname, errmsg);
  endif

  [F, errmsg] = pdFrame (Mdl, Vars, Labels, Data, QP, NumObs, CatPred, ...
                         OutCols);
  if (! isempty (errmsg))
    error ("%s: %s", fname, errmsg);
  endif

  pd = pdValues (Mdl, F, predArgs, subIcept);
  x = F.QP{1};
  if (numel (F.Vars) == 2)
    y = F.QP{2};
  else
    y = [];
  endif

  ## What the drawing needs and the values do not: the frame it averaged over,
  ## the pairs predict was given, and the two options only the plot takes
  opts.Conditional = Cond;
  opts.Parent = Ax;
  opts.PredArgs = predArgs;
  opts.Frame = F;

endfunction

## The options that belong to a kind of model rather than to every one.
function [predArgs, subIcept, OutCols, errmsg] = pdOptions (Mdl, Inter, ...
                                                            Icept, OutCols, ...
                                                            CatPred)

  predArgs = {};
  subIcept = false;
  errmsg = '';
  isfh = is_function_handle (Mdl);
  if (isfh)
    props = {};
  else
    props = properties (Mdl);
  endif
  isgam = any (strcmp (props, 'Interactions'));

  if (! isempty (Inter))
    if (! isgam)
      errmsg = strcat ("'IncludeInteractions' applies only to a", ...
                       " generalized additive model.");
      return;
    endif
    if (! (islogical (Inter) && isscalar (Inter)))
      errmsg = "'IncludeInteractions' must be true or false.";
      return;
    endif
    predArgs = {'IncludeInteractions', Inter};
  endif

  if (! (islogical (Icept) && isscalar (Icept)))
    errmsg = "'IncludeIntercept' must be true or false.";
    return;
  endif
  if (! Icept)
    if (! isgam)
      errmsg = strcat ("'IncludeIntercept' applies only to a generalized", ...
                       " additive model.");
      return;
    endif
    subIcept = true;
  endif

  if (! isfh && ! isempty (CatPred))
    errmsg = strcat ("'CategoricalPredictors' applies only to a function", ...
                     " handle, a model carrying its own.");
    return;
  endif

  if (ischar (OutCols) && strcmpi (OutCols, 'all'))
    OutCols = [];
  elseif (isnumeric (OutCols) && isreal (OutCols) && isvector (OutCols)
          && all (OutCols == fix (OutCols)) && all (OutCols >= 1))
    OutCols = double (OutCols(:)');
  else
    errmsg = strcat ("'OutputColumns' must be a vector of positive", ...
                     " integers or 'all'.");
    return;
  endif
  if (! isfh && ! isempty (OutCols))
    errmsg = "'OutputColumns' applies only to a function handle.";
  endif

endfunction

## Resolve everything partial dependence is computed over.
##
## This is the one place the arguments of partialDependence and
## plotPartialDependence are read, so that the two cannot drift apart and so
## that a later change to what a model accepts as Data is made here and nowhere
## else.
##
## Mdl is a fitted model or a function handle.  Vars names one or two
## predictors, by index or by name.  Labels holds the classes to answer for and
## is empty for a regression model.  Data is the data to average over, empty to
## take the model's own.  QP holds the query points, empty to take the default.
## NumObs is how many observations to sample, empty for all.  CatPred names the
## categorical predictors of a function handle, which has no model to carry
## them.  OutCols holds the columns to keep of what a function handle answers,
## empty for all of them.
##
## F is a structure with the fields IsClass, whether scores are averaged rather
## than responses; PredictorNames and NumPredictors; Cat, the indices of the
## categorical predictors; X, the observations averaged over, already sampled;
## Vars, the one or two resolved indices; QP, a cell holding the query points
## of each; ClassNames and LabelIdx, the classes answered for and where they
## sit among the model's own; and OutCols, as given.
##
## errmsg is the body of the message the caller should raise, empty when there
## is nothing wrong; the caller emits it under its own name.
function [F, errmsg] = pdFrame (Mdl, Vars, Labels, Data, QP, NumObs, ...
                                CatPred, OutCols)

  F = [];
  errmsg = '';
  if (nargin < 8)
    OutCols = [];
  endif
  isfh = is_function_handle (Mdl);

  ## What the model can tell us about itself
  if (isfh)
    props = {};
  else
    props = properties (Mdl);
  endif
  has = @(n) any (strcmp (props, n));

  ## A generalized additive model answers over the observations it was fitted
  ## on, so a compact one cannot answer at all, and DATA does not help:
  ## R2024a passes over DATA for such a model, measured 2026-09-17
  if (! isfh && any (strcmp (props, 'Interactions'))
      && ! (any (strcmp (props, 'X')) && ! isempty (Mdl.X)))
    errmsg = strcat ("partial dependence of a generalized additive model", ...
                     " is taken over the observations it was fitted on,", ...
                     " which a compact one does not keep; use the model it", ...
                     " was compacted from.");
    return;
  endif

  ## The data averaged over.  A model that kept its training data supplies it,
  ## a compact one cannot and must be given it.
  if (isempty (Data))
    if (isfh)
      errmsg = "DATA is required when the model is a function handle.";
      return;
    endif
    [Data, errmsg] = pdModelData (Mdl, props);
    if (! isempty (errmsg))
      return;
    endif
  endif
  if (! (isnumeric (Data) && isreal (Data) && ismatrix (Data)
         && ndims (Data) == 2 && ! isempty (Data)))
    errmsg = "DATA must be a real numeric matrix.";
    return;
  endif
  p = columns (Data);

  ## The predictor names, from the model where it has them
  if (has ('PredictorNames') && ! isempty (Mdl.PredictorNames))
    pnames = Mdl.PredictorNames(:)';
    if (numel (pnames) != p)
      errmsg = strcat ("DATA must have one column per predictor of the", ...
                       " model.");
      return;
    endif
  else
    pnames = arrayfun (@(k) sprintf ('x%d', k), 1:p, 'UniformOutput', false);
  endif

  ## The categorical predictors: the model's own, or those named for a handle
  if (isfh)
    spec = CatPred;
  elseif (has ('CategoricalPredictors'))
    spec = Mdl.CategoricalPredictors;
  else
    spec = [];
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

  ## VARS, by index or by name, one predictor or two
  [vidx, errmsg] = pdVars (Vars, pnames, p);
  if (! isempty (errmsg))
    return;
  endif

  ## A classifier answers with scores, and needs to be told for which classes
  isclass = has ('ClassNames') && ! isempty (Mdl.ClassNames);
  lidx = [];
  cn = [];
  if (isclass)
    cn = Mdl.ClassNames;
    if (isempty (Labels))
      errmsg = "LABELS is required for a classification model.";
      return;
    endif
    ## namedClasses reports against a 'ClassNames' option, which is not what
    ## is being named here, so the message is this function's own
    [lidx, msg] = namedClasses (cn, Labels);
    if (! isempty (msg) || isempty (lidx) || any (lidx == 0))
      errmsg = "LABELS names a class the model was not fitted with.";
      return;
    endif
  elseif (! isempty (Labels))
    errmsg = "LABELS applies only to a classification model.";
    return;
  endif

  ## Sample the observations before the query points are taken from them, as
  ## the default query points span what was sampled and not what was given
  if (! isempty (NumObs))
    if (! (isnumeric (NumObs) && isscalar (NumObs) && isreal (NumObs)
           && NumObs >= 1 && NumObs == fix (NumObs)))
      errmsg = "'NumObservationsToSample' must be a positive integer.";
      return;
    endif
    n = rows (Data);
    if (NumObs < n)
      Data = Data(randperm (n, double (NumObs)), :);
    endif
  endif

  ## The query points of each variable
  [qp, errmsg] = pdQueryPoints (QP, Data, vidx, cat);
  if (! isempty (errmsg))
    return;
  endif

  ## Assigned field by field: struct () spreads a cell argument into a struct
  ## array, and a cellstr ClassNames would make one
  F.IsClass = isclass;
  F.PredictorNames = pnames;
  F.NumPredictors = p;
  F.Cat = cat;
  F.X = Data;
  F.Vars = vidx;
  F.QP = qp;
  F.ClassNames = cn;
  F.LabelIdx = lidx;
  F.OutCols = OutCols;

endfunction

## The observations a model kept, which the learners hold as a matrix and the
## regression models fitted from a formula hold as a table of every variable.
function [Data, errmsg] = pdModelData (Mdl, props)

  Data = [];
  errmsg = '';
  has = @(n) any (strcmp (props, n));

  if (has ('X') && ! isempty (Mdl.X))
    Data = Mdl.X;
    return;
  endif

  if (has ('Variables') && has ('PredictorNames') && istable (Mdl.Variables))
    pn = Mdl.PredictorNames(:)';
    V = Mdl.Variables;
    cols = cell (1, numel (pn));
    for k = 1:numel (pn)
      c = V.(pn{k});
      if (! ((isnumeric (c) || islogical (c)) && isreal (c)
             && columns (c) == 1))
        errmsg = strcat ("DATA is required where the predictors the model", ...
                         " kept are not all numeric.");
        return;
      endif
      cols{k} = double (c);
    endfor
    Data = [cols{:}];
    return;
  endif

  errmsg = strcat ("DATA is required for a model that does not keep the", ...
                   " observations it was fitted on.");

endfunction

## The query points a predictor takes from the observations themselves.
function [q, errmsg] = pdDefaultPoints (Data, v, iscat)

  q = [];
  errmsg = '';
  col = Data(:,v);
  col = col(! isnan (col));
  if (isempty (col))
    errmsg = "DATA holds no value of a variable to vary.";
    return;
  endif
  if (iscat)
    q = unique (col)(:);
  elseif (min (col) == max (col))
    q = min (col);
  else
    q = linspace (min (col), max (col), 100)(:);
  endif

endfunction

## VARS names one predictor or two, by index or by name.
function [vidx, errmsg] = pdVars (Vars, pnames, p)

  vidx = [];
  errmsg = '';
  if (isempty (Vars))
    errmsg = "VARS must name one or two predictors.";
    return;
  endif

  if (isnumeric (Vars) && isreal (Vars) && isvector (Vars))
    vidx = double (Vars(:)');
    if (! (all (vidx == fix (vidx)) && all (vidx >= 1) && all (vidx <= p)))
      errmsg = "VARS must index the predictors of the model.";
      vidx = [];
      return;
    endif
  elseif (ischar (Vars) || iscellstr (Vars) || isa (Vars, 'string'))
    if (ischar (Vars))
      wanted = cellstr (Vars)';
    elseif (iscellstr (Vars))
      wanted = Vars(:)';
    else
      wanted = cellstr (Vars(:)');
    endif
    vidx = zeros (1, numel (wanted));
    for k = 1:numel (wanted)
      j = find (strcmp (pnames, wanted{k}));
      if (isempty (j))
        errmsg = sprintf ("VARS does not name a predictor: '%s'", wanted{k});
        vidx = [];
        return;
      endif
      vidx(k) = j(1);
    endfor
  else
    errmsg = strcat ("VARS must be one or two predictor indices or", ...
                     " predictor names.");
    return;
  endif

  if (numel (vidx) < 1 || numel (vidx) > 2)
    errmsg = "VARS must name one or two predictors.";
    vidx = [];
  elseif (numel (vidx) == 2 && vidx(1) == vidx(2))
    errmsg = "VARS must name two different predictors.";
    vidx = [];
  endif

endfunction

## The query points of each variable, given or taken from the observations.
function [qp, errmsg] = pdQueryPoints (QP, Data, vidx, cat)

  qp = {};
  errmsg = '';
  nv = numel (vidx);

  ## The levels of a categorical predictor are its query points, and stay so
  ## whatever was asked for: R2024a takes all of them and passes over what it
  ## was given, measured 2026-09-17
  isc = false (1, nv);
  for k = 1:nv
    isc(k) = any (cat == vidx(k));
  endfor

  if (isempty (QP) || all (isc))
    qp = cell (1, nv);
    for k = 1:nv
      [qp{k}, errmsg] = pdDefaultPoints (Data, vidx(k), isc(k));
      if (! isempty (errmsg))
        qp = {};
        return;
      endif
    endfor
    return;
  endif

  ## Given: a vector for one variable, and for two either a two-column matrix
  ## or a cell holding a vector each, which is how they may differ in length
  if (iscell (QP))
    if (numel (QP) != nv)
      errmsg = strcat ("'QueryPoints' given as a cell must hold one vector", ...
                       " per variable.");
      return;
    endif
    qp = cell (1, nv);
    for k = 1:nv
      if (! (isnumeric (QP{k}) && isreal (QP{k}) && isvector (QP{k})
             && ! isempty (QP{k})))
        errmsg = "each 'QueryPoints' vector must be a real numeric vector.";
        qp = {};
        return;
      endif
      if (isc(k))
        [qp{k}, errmsg] = pdDefaultPoints (Data, vidx(k), true);
        if (! isempty (errmsg))
          qp = {};
          return;
        endif
      else
        qp{k} = double (QP{k}(:));
      endif
    endfor
    return;
  endif

  if (! (isnumeric (QP) && isreal (QP) && ! isempty (QP)))
    errmsg = "'QueryPoints' must be a real numeric vector or matrix.";
    return;
  endif
  if (nv == 1)
    if (! isvector (QP))
      errmsg = "'QueryPoints' must be a vector for a single variable.";
      return;
    endif
    v = double (QP(:));
    qp = {v};
  else
    if (columns (QP) != 2)
      errmsg = strcat ("'QueryPoints' must have one column per variable,", ...
                       " or be a cell holding a vector for each.");
      return;
    endif
    v1 = double (QP(:,1));
    v2 = double (QP(:,2));
    qp = {v1, v2};
    for k = 1:2
      if (isc(k))
        [qp{k}, errmsg] = pdDefaultPoints (Data, vidx(k), true);
        if (! isempty (errmsg))
          qp = {};
          return;
        endif
      endif
    endfor
  endif

endfunction

## Average what a model answers over a set of observations, one query point of
## the varied predictors at a time.
##
## F is the frame pdFrame resolved.  predArgs holds the pairs forwarded to
## predict, which is where a generalized additive model is told whether to
## include its interaction terms.  subIntercept takes the model's intercept off
## the result, which is what excluding it means for such a model.
##
## pd is 1-by-numX for a regression model varying one predictor, numY-by-numX
## for two, and gains a leading dimension of one row per class where the model
## is a classifier, as MATLAB R2024a returns them.
function pd = pdValues (Mdl, F, predArgs, subIntercept)

  if (nargin < 3)
    predArgs = {};
  endif
  if (nargin < 4)
    subIntercept = false;
  endif

  qx = F.QP{1};
  nx = numel (qx);
  v1 = F.Vars(1);
  two = numel (F.Vars) == 2;
  if (two)
    qy = F.QP{2};
    ny = numel (qy);
    v2 = F.Vars(2);
  else
    ny = 1;
  endif

  ## One column of answers per class, or the single column of a response
  m = pdWidth (Mdl, F, predArgs);

  ## A tree, and an ensemble of bagged trees, answer over the distribution
  ## they were fitted on rather than over any set of observations
  [trees, leaf, tw] = pdTrees (Mdl);
  overTrees = ! isempty (trees);

  ## A generalized additive model answers over the distribution it was fitted
  ## on as a tree does, measured on R2024a 2026-09-17, but carries no node
  ## table to walk, so the observations it kept stand in for that
  ## distribution.  The query points still come from the data given.
  own = pdOwnData (Mdl);
  pd = zeros (m, ny, nx);
  if (isempty (own))
    Z = F.X;
  else
    Z = own;
  endif
  for ii = 1:nx
    if (! overTrees)
      Z(:,v1) = qx(ii);
    endif
    for jj = 1:ny
      if (overTrees)
        if (two)
          qv = [qx(ii), qy(jj)];
        else
          qv = qx(ii);
        endif
        s = zeros (1, columns (leaf{1}));
        for t = 1:numel (trees)
          s += tw(t) * pdTraverse (trees{t}, F.Vars, qv, leaf{t});
        endfor
        if (F.IsClass)
          s = s(F.LabelIdx);
        endif
        pd(:,jj,ii) = s(:);
      else
        if (two)
          Z(:,v2) = qy(jj);
        endif
        s = pdScore (Mdl, Z, F, predArgs);
        pd(:,jj,ii) = mean (s, 1)(:);
      endif
    endfor
  endfor

  if (subIntercept)
    pd -= Mdl.Intercept;
  endif

  ## A regression model keeps no leading dimension, and one variable no middle
  ## one, so that the result reads as MATLAB documents it
  if (m == 1)
    pd = reshape (pd, ny, nx);
  elseif (! two)
    pd = reshape (pd, m, nx);
  endif

endfunction

## The observations a model takes its own average over, empty where it takes
## the average over the data it was given.
function D = pdOwnData (Mdl)

  D = [];
  if (is_function_handle (Mdl))
    return;
  endif
  props = properties (Mdl);
  isgam = any (strcmp (props, 'Interactions'));
  if (! isgam)
    return;
  endif
  if (any (strcmp (props, 'X')) && ! isempty (Mdl.X))
    D = Mdl.X;
  endif

endfunction

## How many columns of answers a single call gives back.
function m = pdWidth (Mdl, F, predArgs)

  if (F.IsClass)
    m = numel (F.LabelIdx);
  elseif (is_function_handle (Mdl))
    s = pdScore (Mdl, F.X, F, predArgs);
    m = columns (s);
  else
    m = 1;
  endif

endfunction

## What the model answers for a set of observations, as a column per class.
function s = pdScore (Mdl, Z, F, predArgs)

  if (is_function_handle (Mdl))
    s = Mdl (Z);
    if (! (isnumeric (s) && isreal (s) && rows (s) == rows (Z)))
      error (strcat ("the function must answer with a real numeric column", ...
                     " per observation."));
    endif
    if (! isempty (F.OutCols))
      if (any (F.OutCols > columns (s)))
        error (strcat ("'OutputColumns' must index the columns the", ...
                       " function answers with."));
      endif
      s = s(:, F.OutCols);
    endif
  elseif (F.IsClass)
    [~, s] = predict (Mdl, Z, predArgs{:});
    s = s(:, F.LabelIdx);
  else
    s = predict (Mdl, Z, predArgs{:});
    s = s(:);
  endif

endfunction

## The trees a model answers through, where partial dependence is taken over
## them rather than over a set of observations.
##
## MATLAB R2024a takes partial dependence over the fitted distribution, and not
## over the data it is handed, for a decision tree and for an ensemble of
## bagged trees.  Measured 2026-09-17: the value does not move when the data is
## replaced, and differs from the average over it.  Every other model averages
## over the data, a boosted ensemble and TreeBagger included, so this returns
## empty for them and the caller averages instead.
##
## trees is a cell array of the trees, leaf a cell array holding what each
## answers at each of its nodes, already transformed, and w the weight of each
## tree, summing to one.
function [trees, leaf, w] = pdTrees (Mdl)

  trees = {};
  leaf = {};
  w = [];
  if (is_function_handle (Mdl))
    return;
  endif

  cls = class (Mdl);
  switch (cls)

    case {'ClassificationTree', 'CompactClassificationTree', ...
          'RegressionTree', 'CompactRegressionTree'}
      trees = {Mdl};
      leaf = {pdLeaf(Mdl)};
      w = 1;

    case {'ClassificationBaggedEnsemble', 'RegressionBaggedEnsemble', ...
          'ClassificationEnsemble', 'RegressionEnsemble', ...
          'CompactClassificationEnsemble', 'CompactRegressionEnsemble'}
      ## Only a bagged ensemble is taken over its trees; a boosted one is
      ## averaged over the data, as MATLAB documents and as measured.
      ## 'CombineWeights' is what says which, and unlike 'Method' every one of
      ## the six classes carries it, the compact regression one included.
      props = properties (Mdl);
      if (any (strcmp (props, 'CombineWeights'))
          && strcmp (Mdl.CombineWeights, 'WeightedAverage'))
        [trees, leaf, w] = pdEnsembleTrees (Mdl);
      endif

  endswitch

endfunction

## What a tree answers at each of its nodes, with its transform applied.
function L = pdLeaf (T)

  props = properties (T);
  if (any (strcmp (props, 'ClassProbability')))
    L = T.ClassProbability;
    if (any (strcmp (props, 'ScoreTransform')) && ! isempty (T.STfun))
      L = T.STfun (L);
    endif
  else
    L = T.NodeMean(:);
    if (any (strcmp (props, 'ResponseTransform')) && ! isempty (T.RTfun))
      L = T.RTfun (L);
    endif
  endif

endfunction

## The learners of a bagged ensemble and the weight of each.
function [trees, leaf, w] = pdEnsembleTrees (Mdl)

  trees = {};
  leaf = {};
  w = [];
  props = properties (Mdl);
  if (! any (strcmp (props, 'Trained')) || isempty (Mdl.Trained))
    return;
  endif
  trees = Mdl.Trained(:)';
  n = numel (trees);
  ## Bagging takes trees, but an ensemble is not obliged to hold them, and
  ## anything else is averaged over the data rather than walked
  istree = cellfun (@(t) any (strcmp (properties (t), 'Children')), trees);
  if (! all (istree))
    trees = {};
    return;
  endif
  leaf = cell (1, n);
  for k = 1:n
    leaf{k} = pdLeaf (trees{k});
  endfor
  if (any (strcmp (props, 'TrainedWeights')) && ! isempty (Mdl.TrainedWeights))
    w = Mdl.TrainedWeights(:)';
  else
    w = ones (1, n);
  endif
  w = w / sum (w);

endfunction

## What a decision tree answers on average, with some predictors held fixed.
##
## This is Friedman's weighted traversal, and it is what partial dependence
## means for a tree: the average of the tree's answer over the distribution it
## was fitted on, with the predictors named by vars held at q.  It is not an
## average over any set of observations and takes none, which is why it answers
## for a compact model, and why MATLAB R2024a gives the same value for a tree
## whatever Data it is handed.
##
## A node that cuts one of the held predictors is followed to the child that
## value reaches.  A node that cuts any other is followed both ways, each child
## taking the share of the training weight that reached it.  A held value that
## matches no level of a categorical cut stops there, as it stops an
## observation in predict.
##
## T is the tree, vars the held predictors and q their values in the same
## order.  leaf holds what the tree answers at each node, one row per node: a
## column of means for a regression tree and a column per class for a
## classifier.  v is a row, one element per column of leaf.
function v = pdTraverse (T, vars, q, leaf)

  v = zeros (1, columns (leaf));
  node = 1;
  w = 1;
  while (! isempty (node))
    n = node(1);
    wn = w(1);
    node(1) = [];
    w(1) = [];
    j = T.CutPredictorIndex(n);
    if (j == 0)
      v += wn * leaf(n,:);
      continue;
    endif
    L = T.Children(n,1);
    R = T.Children(n,2);
    k = find (vars == j, 1);
    if (isempty (k))
      ## Not held: both ways, in the proportion the training weight took
      pn = T.NodeProbability(n);
      if (pn > 0)
        pL = T.NodeProbability(L) / pn;
      else
        pL = 0.5;
      endif
      node(end+1:end+2) = [L, R];
      w(end+1:end+2) = [wn * pL, wn * (1 - pL)];
    else
      [child, stop] = pdChild (T, n, q(k), L, R);
      if (stop)
        v += wn * leaf(n,:);
      else
        node(end+1) = child;
        w(end+1) = wn;
      endif
    endif
  endwhile

endfunction

## Which child a held value reaches, and whether it reaches neither.
function [child, stop] = pdChild (T, n, val, L, R)

  stop = false;
  child = 0;
  if (! isempty (T.CutCategories) && ! isempty (T.CutCategories{n,1}))
    if (any (T.CutCategories{n,1} == val))
      child = L;
    elseif (any (T.CutCategories{n,2} == val))
      child = R;
    else
      stop = true;
    endif
  elseif (isnan (val))
    stop = true;
  elseif (val < T.CutPoint(n))
    child = L;
  else
    child = R;
  endif

endfunction
