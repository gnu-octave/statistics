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
## @deftypefn {Private Function} {[@var{trees}, @var{leaf}, @var{w}] =} pdTrees (@var{Mdl})
##
## The trees a model answers through, where partial dependence is taken over
## them rather than over a set of observations.
##
## MATLAB R2024a takes partial dependence over the fitted distribution, and
## not over the data it is handed, for a decision tree and for an ensemble of
## bagged trees.  Measured 2026-09-17: the value does not move when the data
## is replaced, and differs from the average over it.  Every other model
## averages over the data, a boosted ensemble and @code{TreeBagger} included,
## so this returns empty for them and the caller averages instead.
##
## @var{trees} is a cell array of the trees, @var{leaf} a cell array holding
## what each answers at each of its nodes, already transformed, and @var{w}
## the weight of each tree, summing to one.
##
## @end deftypefn

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
