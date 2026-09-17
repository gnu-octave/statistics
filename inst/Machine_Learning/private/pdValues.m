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
## @deftypefn {Private Function} {@var{pd} =} pdValues (@var{Mdl}, @var{F}, @var{predArgs}, @var{subIntercept})
##
## Average what a model answers over a set of observations, one query point of
## the varied predictors at a time.
##
## @var{F} is the frame @code{pdFrame} resolved.  @var{predArgs} holds the
## pairs forwarded to @code{predict}, which is where a generalized additive
## model is told whether to include its interaction terms.  @var{subIntercept}
## takes the model's intercept off the result, which is what excluding it
## means for such a model.
##
## @var{pd} is @math{1xnumX} for a regression model varying one predictor,
## @math{numYxnumX} for two, and gains a leading dimension of one row per
## class where the model is a classifier, as MATLAB R2024a returns them.
##
## @end deftypefn

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
