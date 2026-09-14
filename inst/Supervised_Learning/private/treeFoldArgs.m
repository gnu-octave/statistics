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
## @deftypefn {Private Function} {@var{args} =} treeFoldArgs (@var{Mdl})
##
## The Name-Value arguments a fold of a @code{ClassificationTree} or a
## @code{RegressionTree} is grown with.
##
## A classifier's fold is given the parent's class names, prior and cost
## rather than being left to re-derive them from its own rows; both are given
## the growth parameters the parent actually used.  The observation weights
## are not here, being sliced per fold by the caller.
##
## @qcode{'ScoreTransform'} and @qcode{'ResponseTransform'} are deliberately
## absent: a transform is applied once to the assembled answer, and a fold
## carrying one would apply it twice.
##
## @qcode{'MaxNumSplits'} defaults to one less than the number of observations,
## so a fold works its own out; a budget the caller actually asked for is passed
## on.  @qcode{'NumVariablesToSample'} is passed on too, so a fold samples its
## predictors as the parent did.  @qcode{'MinParentSize'} is given the value the
## fit settled on rather than the one asked for, which reproduces it: the
## constructor takes the larger of it and twice @qcode{'MinLeafSize'}, and the
## larger is already there.  A tree with categorical predictors passes them on,
## with @qcode{'MaxNumCategories'} and, for a classifier that was given one,
## @qcode{'AlgorithmForCategorical'}.
##
## @seealso{ClassificationTree, ClassificationPartitionedModel}
## @end deftypefn

function args = treeFoldArgs (Mdl)

  MP = Mdl.ModelParameters;
  args = {'PredictorNames', Mdl.PredictorNames, ...
          'ResponseName', Mdl.ResponseName};

  if (strcmp (MP.Type, 'classification'))
    args = [args, {'ClassNames', Mdl.ClassNames, ...
                   'Prior', Mdl.Prior, 'Cost', Mdl.Cost}];
  else
    args = [args, {'QuadraticErrorTolerance', MP.QEToler}];
  endif

  args = [args, {'SplitCriterion', MP.SplitCriterion, ...
                 'MinParentSize', MP.MinParent, ...
                 'MinLeafSize', MP.MinLeaf, ...
                 'MergeLeaves', MP.MergeLeaves, ...
                 'NumVariablesToSample', MP.NVarToSample, ...
                 'Prune', MP.Prune, ...
                 'PruneCriterion', MP.PruneCriterion}];

  if (MP.MaxSplits != Mdl.NumObservations - 1)
    args = [args, {'MaxNumSplits', MP.MaxSplits}];
  endif

  if (! isempty (Mdl.CategoricalPredictors))
    args = [args, {'CategoricalPredictors', Mdl.CategoricalPredictors, ...
                   'MaxNumCategories', MP.MaxCat}];
    if (strcmp (MP.Type, 'classification') && ! strcmpi (MP.AlgCat, 'auto'))
      args = [args, {'AlgorithmForCategorical', MP.AlgCat}];
    endif
  endif

endfunction
