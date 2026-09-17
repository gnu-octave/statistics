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
## @deftypefn {Private Function} {[@var{F}, @var{errmsg}] =} pdFrame (@var{Mdl}, @var{Vars}, @var{Labels}, @var{Data}, @var{QP}, @var{NumObs}, @var{CatPred})
##
## Resolve everything partial dependence is computed over.
##
## This is the one place the arguments of @code{partialDependence} and
## @code{plotPartialDependence} are read, so that the two cannot drift apart
## and so that a later change to what a model accepts as @var{Data} is made
## here and nowhere else.
##
## @var{Mdl} is a fitted model or a function handle.  @var{Vars} names one or
## two predictors, by index or by name.  @var{Labels} holds the classes to
## answer for and is empty for a regression model.  @var{Data} is the data to
## average over, empty to take the model's own.  @var{QP} holds the query
## points, empty to take the default.  @var{NumObs} is how many observations
## to sample, empty for all.  @var{CatPred} names the categorical predictors
## of a function handle, which has no model to carry them.
##
## @var{F} is a structure with the fields @qcode{IsClass}, whether scores are
## averaged rather than responses; @qcode{PredictorNames} and
## @qcode{NumPredictors}; @qcode{Cat}, the indices of the categorical
## predictors; @qcode{X}, the observations averaged over, already sampled;
## @qcode{Vars}, the one or two resolved indices; @qcode{QP}, a cell holding
## the query points of each; @qcode{ClassNames} and @qcode{LabelIdx}, the
## classes answered for and where they sit among the model's own.
##
## @var{errmsg} is the body of the message the caller should raise, empty when
## there is nothing wrong; the caller emits it under its own name.
##
## @end deftypefn

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
