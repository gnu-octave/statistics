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
  ## @item @qcode{'Method'} @tab @tab The algorithm, @qcode{'interventional'}
  ## by default.
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
    ## A character vector.  Only @qcode{'interventional-kernel'} is
    ## implemented, which computes the values exactly by enumerating every
    ## subset of the predictors.  This property is read-only.
    ##
    ## @end deftp
    Method = '';

    ## -*- texinfo -*-
    ## @deftp {shapley} {property} NumSubsets
    ##
    ## How many predictor subsets the values were computed over
    ##
    ## Every subset is used, so this is two raised to the number of
    ## predictors.  This property is read-only.
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
      if (! (isa (blackbox, 'PredictiveModel') || is_function_handle (blackbox)))
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

      [method, errmsg] = shapMethod (Method, MaxSub);
      if (! isempty (errmsg))
        error ("shapley: %s", errmsg);
      endif

      [F, errmsg] = shapFrame (blackbox, Data, CatPred, NumObs);
      if (! isempty (errmsg))
        error ("shapley: %s", errmsg);
      endif

      this.BlackboxModel = blackbox;
      this.X = F.X;
      this.CategoricalPredictors = F.Cat;
      this.SampledObservationIndices = F.Idx;
      this.Method = method;
      this.NumSubsets = 2 ^ F.NumPredictors;
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

      ## Every subset of the predictors, which is what makes the values exact
      phi = zeros (M, K, nq);
      icept = [];
      for ii = 1:nq
        V = shapSubsetValues (sfcn, Xs, QueryPoints(ii,:), M, K);
        phi(:,:,ii) = shapFromValues (V, M, K);
        if (ii == 1)
          icept = V(1,:);
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
function [method, errmsg] = shapMethod (Method, MaxSub)

  method = 'interventional-kernel';
  errmsg = '';
  if (! isempty (MaxSub))
    errmsg = strcat ("'MaxNumSubsets' is not implemented; every subset", ...
                     " is used.");
    return;
  endif
  if (isempty (Method))
    return;
  endif
  if (! (ischar (Method) || isa (Method, 'string')) || ! isrow (char (Method)))
    errmsg = "'Method' must be a character vector or a string scalar.";
    return;
  endif
  switch (lower (char (Method)))
    case 'interventional'
      ## the only one implemented
    case 'conditional'
      errmsg = "'Method' value 'conditional' is not implemented.";
    otherwise
      errmsg = strcat ("'Method' must be 'interventional' or", ...
                       " 'conditional'.");
  endswitch

endfunction

## Resolve the observations, the predictors and the classes.
function [F, errmsg] = shapFrame (blackbox, Data, CatPred, NumObs)

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

  ## Enumerating every subset is what keeps the values exact, and the cap
  ## MATLAB puts on the count belongs to the sampling this does not do
  if (p > 10)
    errmsg = strcat ("a model of more than 10 predictors needs the", ...
                     " subset sampling that 'MaxNumSubsets' asks for,", ...
                     " which is not implemented.");
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
  F.Cat = cat;
  F.X = Data;
  F.Idx = idx;

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

## The value function over every subset of the predictors.  Bit II of the row
## index, counted from zero, says whether predictor II is held at the query
## point; the rest of the columns stay as the observations have them, which is
## what makes the average an interventional one.
function V = shapSubsetValues (sfcn, Xs, q, M, K)

  nS = 2 ^ M;
  n = rows (Xs);
  V = zeros (nS, K);
  for m = 0:(nS - 1)
    mask = logical (bitget (m, 1:M));
    Z = Xs;
    if (any (mask))
      Z(:,mask) = repmat (q(mask), n, 1);
    endif
    V(m + 1,:) = mean (sfcn (Z), 1);
  endfor

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
%!error<shapley: 'MaxNumSubsets' is not implemented; every subset is used.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'MaxNumSubsets', 2)
%!error<shapley: 'Method' value 'conditional' is not implemented.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'Method', 'conditional')
%!error<shapley: 'Method' must be 'interventional' or 'conditional'.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'Method', 'marginal')
%!error<shapley: 'NumObservationsToSample' must be a positive integer or 'all'.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'NumObservationsToSample', 0)
%!error<shapley: a model of more than 10 predictors needs the subset sampling that 'MaxNumSubsets' asks for, which is not implemented.> shapley (@(Z) Z(:,1), ones (3, 11))
%!error<shapley: unknown optional argument or misplaced value.> shapley (@(Z) Z(:,1), [1, 2; 3, 4], 'NoSuchThing', 1)
%!error<shapley.fit: QUERYPOINTS must be a real numeric matrix.> fit (shapley (@(Z) Z(:,1), [1, 2; 3, 4]), 'abc')
%!error<shapley.fit: QUERYPOINTS must have one column per predictor of the model.> fit (shapley (@(Z) Z(:,1), [1, 2; 3, 4]), [1, 2, 3])
