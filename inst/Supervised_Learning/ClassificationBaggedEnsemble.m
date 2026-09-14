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

classdef ClassificationBaggedEnsemble < ClassificationEnsemble
  ## -*- texinfo -*-
  ## @deftp {statistics} ClassificationBaggedEnsemble
  ##
  ## Bagged ensemble of decision trees for classification
  ##
  ## A @code{ClassificationBaggedEnsemble} object holds decision trees each
  ## grown on a sample drawn from the training data in proportion to the
  ## observation weights, each split chosen from @code{ceil (sqrt (P))}
  ## predictors drawn afresh at every node.  It predicts by averaging its
  ## trees' class probabilities.
  ##
  ## Create one with @code{fitcensemble} and @qcode{'Method'} set to
  ## @qcode{'Bag'}.  It carries everything a @code{ClassificationEnsemble}
  ## does, and which rows each tree drew.  A boosting method that resamples,
  ## asked for with @qcode{'Resample'}, @qcode{'FResample'} or
  ## @qcode{'Replace'}, returns this class too: its trees are boosted and
  ## summed by their weights, each grown on the rows it drew.
  ##
  ## @seealso{fitcensemble, ClassificationEnsemble,
  ## CompactClassificationEnsemble, TreeBagger}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationBaggedEnsemble} {property} FResample
    ##
    ## Share of the observations each tree draws
    ##
    ## A number greater than 0 and no greater than 1, each tree drawing
    ## @code{ceil (FResample * N)} observations.  This property is read-only.
    ##
    ## @end deftp
    FResample = 1;

    ## -*- texinfo -*-
    ## @deftp {ClassificationBaggedEnsemble} {property} Replace
    ##
    ## Whether the trees draw with replacement
    ##
    ## A logical scalar, true by default.  This property is read-only.
    ##
    ## @end deftp
    Replace = true;

    ## -*- texinfo -*-
    ## @deftp {ClassificationBaggedEnsemble} {property} UseObsForLearner
    ##
    ## Which observations each tree drew
    ##
    ## An @math{NxNumTrained} logical matrix, true where a tree's sample holds
    ## an observation.  This property is read-only.
    ##
    ## @end deftp
    UseObsForLearner = [];

  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{obj} =} ClassificationBaggedEnsemble (@var{X}, @var{Y})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{obj} =} ClassificationBaggedEnsemble (@dots{}, @var{name}, @var{value})
    ##
    ## Fit a bagged ensemble of decision trees.
    ##
    ## @code{fitcensemble} with @qcode{'Method'} set to @qcode{'Bag'} is the
    ## documented way in, and its help lists the options both take.
    ##
    ## @seealso{fitcensemble, ClassificationEnsemble}
    ## @end deftypefn
    function this = ClassificationBaggedEnsemble (X, Y, varargin)

      if (nargin < 2)
        error ("ClassificationBaggedEnsemble: too few input arguments.");
      endif
      if (! any (cellfun (@(a) ischar (a) && strcmpi (a, 'Method'), ...
                          varargin(1:2:end))))
        varargin(end+1:end+2) = {'Method', 'Bag'};
      endif
      this = this@ClassificationEnsemble (X, Y, varargin{:});
      this = bagProperties (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate a bagged ensemble.
    ##
    ## Behaves as @code{ClassificationEnsemble.crossval}, returning a
    ## @code{ClassificationPartitionedEnsemble}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, ClassificationPartitionedEnsemble}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)
      CVMdl = crossval@ClassificationEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationBaggedEnsemble} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Drop the training data from a bagged ensemble.
    ##
    ## Returns a @code{CompactClassificationEnsemble}, as
    ## @code{ClassificationEnsemble.compact} does.
    ##
    ## @seealso{ClassificationBaggedEnsemble, CompactClassificationEnsemble}
    ## @end deftypefn
    function CMdl = compact (this)
      CMdl = compact@ClassificationEnsemble (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{B} =} resume (@var{obj}, @var{NumLearningCycles})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{B} =} resume (@dots{}, 'NPrint', @var{n})
    ##
    ## Grow more trees in a bagged ensemble.
    ##
    ## Behaves as @code{ClassificationEnsemble.resume}, the new trees' samples
    ## added to @code{UseObsForLearner}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, ClassificationEnsemble.resume}
    ## @end deftypefn
    function this = resume (this, varargin)
      this = resume@ClassificationEnsemble (this, varargin{:});
      this = bagProperties (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{label} =} predict (@var{obj}, @var{X})
    ## @deftypefnx {ClassificationBaggedEnsemble} {[@var{label}, @var{scores}] =} predict (@dots{})
    ## @deftypefnx {ClassificationBaggedEnsemble} {[@dots{}] =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Classify new data with a bagged ensemble.
    ##
    ## Behaves as @code{CompactClassificationEnsemble.predict}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, CompactClassificationEnsemble.predict}
    ## @end deftypefn
    function [label, scores] = predict (this, varargin)
      [label, scores] = predict@ClassificationEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of a bagged ensemble.
    ##
    ## Behaves as @code{CompactClassificationEnsemble.loss}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, CompactClassificationEnsemble.loss}
    ## @end deftypefn
    function L = loss (this, varargin)
      L = loss@ClassificationEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{e} =} edge (@dots{}, @var{name}, @var{value})
    ##
    ## Classification edge of a bagged ensemble.
    ##
    ## Behaves as @code{CompactClassificationEnsemble.edge}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, CompactClassificationEnsemble.edge}
    ## @end deftypefn
    function e = edge (this, varargin)
      e = edge@ClassificationEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{m} =} margin (@dots{}, @var{name}, @var{value})
    ##
    ## Classification margins of a bagged ensemble.
    ##
    ## Behaves as @code{CompactClassificationEnsemble.margin}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, CompactClassificationEnsemble.margin}
    ## @end deftypefn
    function m = margin (this, varargin)
      m = margin@ClassificationEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{label} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {[@var{label}, @var{scores}] =} resubPredict (@dots{})
    ##
    ## Classify the training data with a bagged ensemble.
    ##
    ## Behaves as @code{ClassificationEnsemble.resubPredict}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, ClassificationEnsemble.resubPredict}
    ## @end deftypefn
    function [label, scores] = resubPredict (this, varargin)
      [label, scores] = resubPredict@ClassificationEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of a bagged ensemble on the training data.
    ##
    ## Behaves as @code{ClassificationEnsemble.resubLoss}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, ClassificationEnsemble.resubLoss}
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      L = resubLoss@ClassificationEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{e} =} resubEdge (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{e} =} resubEdge (@dots{}, @var{name}, @var{value})
    ##
    ## Classification edge of a bagged ensemble on the training data.
    ##
    ## Behaves as @code{ClassificationEnsemble.resubEdge}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, ClassificationEnsemble.resubEdge}
    ## @end deftypefn
    function e = resubEdge (this, varargin)
      e = resubEdge@ClassificationEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{m} =} resubMargin (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{m} =} resubMargin (@dots{}, @var{name}, @var{value})
    ##
    ## Classification margins of the training data of a bagged ensemble.
    ##
    ## Behaves as @code{ClassificationEnsemble.resubMargin}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, ClassificationEnsemble.resubMargin}
    ## @end deftypefn
    function m = resubMargin (this, varargin)
      m = resubMargin@ClassificationEnsemble (this, varargin{:});
    endfunction


    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{imp} =} predictorImportance (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {[@var{imp}, @var{ma}] =} predictorImportance (@var{obj})
    ##
    ## Estimate the importance of each predictor.
    ##
    ## The mean over the trees of each tree's @code{predictorImportance}, as
    ## @code{CompactClassificationEnsemble.predictorImportance} computes it.
    ##
    ## @seealso{ClassificationBaggedEnsemble,
    ## ClassificationBaggedEnsemble.oobPermutedPredictorImportance}
    ## @end deftypefn
    function [imp, ma] = predictorImportance (this)
      [imp, ma] = predictorImportance@ClassificationEnsemble (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{label} =} oobPredict (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {[@var{label}, @var{scores}] =} oobPredict (@dots{})
    ## @deftypefnx {ClassificationBaggedEnsemble} {[@dots{}] =} oobPredict (@dots{}, 'Learners', @var{idx})
    ##
    ## Out-of-bag predictions for the training data.
    ##
    ## Each training observation is classified by the trees whose samples
    ## left it out, as @code{predict} classifies it with
    ## @qcode{'UseObsForLearner'} set to @code{! UseObsForLearner}.  An
    ## observation in the sample of every tree used has @code{NaN} scores and
    ## is given the class of greatest prior probability.  @qcode{'Learners'}
    ## restricts the trees.
    ##
    ## @seealso{ClassificationBaggedEnsemble, ClassificationBaggedEnsemble.oobLoss}
    ## @end deftypefn
    function [label, scores] = oobPredict (this, varargin)
      caller = 'ClassificationBaggedEnsemble.oobPredict';
      args = oobArgs (varargin, {'Learners'}, caller);
      [label, scores] = ensemblePredict (compact (this), this.X, ...
                                         [{'UseObsForLearner', ...
                                           ! this.UseObsForLearner}, args], ...
                                         caller);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{L} =} oobLoss (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{L} =} oobLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-bag classification loss.
    ##
    ## The loss of the out-of-bag scores against @code{Y}, weighted by
    ## @code{W}, an observation in the sample of every tree used being left
    ## out.  @qcode{'LossFun'} and @qcode{'Mode'} are taken as by
    ## @code{CompactClassificationEnsemble.loss}, and @qcode{'Learners'}
    ## restricts the trees.
    ##
    ## @seealso{ClassificationBaggedEnsemble, ClassificationBaggedEnsemble.oobPredict}
    ## @end deftypefn
    function L = oobLoss (this, varargin)
      caller = 'ClassificationBaggedEnsemble.oobLoss';
      args = oobArgs (varargin, {'Learners', 'LossFun', 'Mode'}, caller);
      L = ensembleLoss (compact (this), this.X, this.Y, ...
                        [{'UseObsForLearner', ! this.UseObsForLearner, ...
                          'Weights', this.W}, args], caller);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{e} =} oobEdge (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{e} =} oobEdge (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-bag classification edge.
    ##
    ## The weighted mean of the out-of-bag margins, weighted by @code{W}.
    ## @qcode{'Mode'} and @qcode{'Learners'} are taken as by @code{oobLoss}.
    ##
    ## @seealso{ClassificationBaggedEnsemble, ClassificationBaggedEnsemble.oobMargin}
    ## @end deftypefn
    function e = oobEdge (this, varargin)
      caller = 'ClassificationBaggedEnsemble.oobEdge';
      args = oobArgs (varargin, {'Learners', 'Mode'}, caller);
      e = ensembleEdge (compact (this), this.X, this.Y, ...
                        [{'UseObsForLearner', ! this.UseObsForLearner, ...
                          'Weights', this.W}, args], caller);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{m} =} oobMargin (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{m} =} oobMargin (@dots{}, 'Learners', @var{idx})
    ##
    ## Out-of-bag classification margins.
    ##
    ## The margin of each training observation under its out-of-bag scores,
    ## @code{NaN} for one in the sample of every tree used.
    ##
    ## @seealso{ClassificationBaggedEnsemble, ClassificationBaggedEnsemble.oobEdge}
    ## @end deftypefn
    function m = oobMargin (this, varargin)
      caller = 'ClassificationBaggedEnsemble.oobMargin';
      args = oobArgs (varargin, {'Learners'}, caller);
      m = ensembleMargin (compact (this), this.X, this.Y, ...
                          [{'UseObsForLearner', ...
                            ! this.UseObsForLearner}, args], caller);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationBaggedEnsemble} {@var{imp} =} oobPermutedPredictorImportance (@var{obj})
    ## @deftypefnx {ClassificationBaggedEnsemble} {@var{imp} =} oobPermutedPredictorImportance (@dots{}, 'Learners', @var{idx})
    ##
    ## Out-of-bag predictor importance by permutation.
    ##
    ## For each tree, the values of each predictor are permuted among the
    ## observations out of its bag, and the tree's misclassification rate on
    ## them, weighted by @code{W}, is taken before and after.  @var{imp} holds,
    ## for each predictor, the mean of the rise over the trees divided by its
    ## standard deviation over the trees, zero where the mean is zero.
    ## @qcode{'Learners'} restricts the trees.
    ##
    ## @seealso{ClassificationBaggedEnsemble,
    ## ClassificationBaggedEnsemble.predictorImportance}
    ## @end deftypefn
    function imp = oobPermutedPredictorImportance (this, varargin)

      caller = 'ClassificationBaggedEnsemble.oobPermutedPredictorImportance';
      args = oobArgs (varargin, {'Learners', 'Options'}, caller);
      learners = 1:this.NumTrained;
      for i = 1:2:numel (args)
        if (strcmpi (args{i}, 'Options'))
          error ("%s: 'Options' is not implemented.", caller);
        endif
        learners = oobLearners (args{i+1}, this.NumTrained, caller);
      endfor
      p = columns (this.X);
      gY = labelIndices (this.ClassNames, this.Y);
      D = zeros (numel (learners), p);
      for j = 1:numel (learners)
        t = learners(j);
        r = find (! this.UseObsForLearner(:,t));
        w = this.W(r);
        if (isempty (r) || ! (sum (w) > 0))
          continue;
        endif
        w /= sum (w);
        tree = this.Trained{t};
        Xo = this.X(r,:);
        miss = @(Z) labelIndices (this.ClassNames, predict (tree, Z)) != gY(r);
        e0 = sum (w .* miss (Xo));
        for v = 1:p
          Xp = Xo;
          Xp(:,v) = Xo(randperm (numel (r)), v);
          D(j,v) = sum (w .* miss (Xp)) - e0;
        endfor
      endfor
      imp = mean (D, 1) ./ std (D, 0, 1);
      imp(mean (D, 1) == 0) = 0;

    endfunction

  endmethods

  methods (Access = private)

    ## The public bag properties, from what the fit keeps in the parent.
    function this = bagProperties (this)
      this.FResample = this.BagFResample;
      this.Replace = this.BagReplace;
      this.UseObsForLearner = this.BagInBag;
    endfunction

  endmethods

endclassdef

## The Name-Value pairs of an out-of-bag method, refused unless in ALLOWED.
function args = oobArgs (args, allowed, caller)

  if (mod (numel (args), 2) != 0)
    error ("%s: name-value arguments must be in pairs.", caller);
  endif
  for i = 1:2:numel (args)
    if (! (ischar (args{i}) && any (strcmpi (args{i}, allowed))))
      error ("%s: invalid parameter name in optional pair arguments.", ...
             caller);
    endif
  endfor

endfunction

## The indices of the learners 'Learners' names.
function learners = oobLearners (val, T, caller)

  if (! (isnumeric (val) && isvector (val) && isreal (val)
         && all (val >= 1) && all (val <= T) && all (val == fix (val))))
    error (strcat ("%s: 'Learners' must be a vector of indices of", ...
                   " trained learners."), caller);
  endif
  learners = double (val(:)');

endfunction

## Test output
%!test  # MATLAB parity: the properties of a bagged ensemble
%! load fisheriris
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 3);
%! assert_equal (numel (properties (Mdl)), 29);
%! assert_equal (Mdl.Method, 'Bag');
%! assert_equal (Mdl.TrainedWeights, [1; 1; 1]);
%! assert_equal (Mdl.FitInfo, []);
%! assert_equal (Mdl.FitInfoDescription, 'None');
%! assert_equal (Mdl.FResample, 1);
%! assert_equal (Mdl.Replace, true);
%! assert_equal (size (Mdl.UseObsForLearner), [150, 3]);
%! assert_equal (isfield (Mdl.ModelParameters, 'LearnRate'), false);

%!test  # MATLAB parity: the scores are the mean of the trees' probabilities
%! load fisheriris
%! rng (5);
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 5);
%! [~, s] = predict (Mdl, meas);
%! P = zeros (150, 3);
%! for t = 1:5
%!   [~, st] = predict (Mdl.Trained{t}, meas);
%!   [~, k] = ismember (Mdl.Trained{t}.ClassNames, Mdl.ClassNames);
%!   P(:,k) += st;
%! endfor
%! assert_equal (s, P / 5, 1e-15);

%!test  # MATLAB parity: a subset of trees is averaged over the subset
%! load fisheriris
%! rng (5);
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 3);
%! [~, s] = predict (Mdl, meas(1:2,:), 'Learners', [1, 3]);
%! [~, a] = predict (Mdl.Trained{1}, meas(1:2,:));
%! [~, b] = predict (Mdl.Trained{3}, meas(1:2,:));
%! assert_equal (s, (a + b) / 2, 1e-15);

%!test  # MATLAB parity: each tree draws FResample of the rows
%! load fisheriris
%! rng (8);
%! Mdl = ClassificationBaggedEnsemble (meas, species, ...
%!                                     'NumLearningCycles', 3, ...
%!                                     'Replace', 'off', 'FResample', 0.5);
%! assert_equal (sum (Mdl.UseObsForLearner), [75, 75, 75]);
%! assert_equal (Mdl.Trained{1}.NodeSize(1), 75);
%! R = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 1, ...
%!                                   'FResample', 0.5);
%! assert_equal (R.Trained{1}.NodeSize(1), 75);

%!test  # MATLAB parity: the bootstrap draws in proportion to the weights
%! load fisheriris
%! rng (9);
%! w = [10 * ones(50, 1); ones(100, 1)];
%! Mdl = ClassificationBaggedEnsemble (meas, species, ...
%!                                     'NumLearningCycles', 20, 'Weights', w);
%! assert_equal (Mdl.Prior, [5, 0.5, 0.5] / 6, 1e-15);
%! share = mean (sum (Mdl.UseObsForLearner(1:50,:)) / 50);
%! assert_equal (share > 0.85, true);

%!test  # resuming adds the new trees' samples
%! load fisheriris
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 2);
%! Mdl = resume (Mdl, 2);
%! assert_equal (class (Mdl), 'ClassificationBaggedEnsemble');
%! assert_equal (Mdl.NumTrained, 4);
%! assert_equal (size (Mdl.UseObsForLearner), [150, 4]);
%! assert_equal (Mdl.ModelParameters.NLearn, 4);

%!test  # MATLAB parity: the cumulative loss of the bagged trees
%! load fisheriris
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 3);
%! L = resubLoss (Mdl, 'Mode', 'cumulative');
%! assert_equal (size (L), [3, 1]);
%! assert_equal (L(end), resubLoss (Mdl), 1e-15);
%! assert_equal (class (compact (Mdl)), 'CompactClassificationEnsemble');

## Test input validation
%!error<ClassificationBaggedEnsemble: too few input arguments.> ...
%! ClassificationBaggedEnsemble (1)
%!error<ClassificationBaggedEnsemble: 'Method' must be 'Bag' unless the ensemble resamples.> ...
%! load fisheriris
%! ClassificationBaggedEnsemble (meas, species, 'Method', 'AdaBoostM2')
%!error<ClassificationBaggedEnsemble: 'LearnRate' cannot be used with the 'Bag' method.> ...
%! load fisheriris
%! ClassificationBaggedEnsemble (meas, species, 'LearnRate', 0.5)
%!error<ClassificationBaggedEnsemble.predict: too few input arguments.> ...
%! load fisheriris
%! predict (ClassificationBaggedEnsemble (meas, species, ...
%!                                        'NumLearningCycles', 1))

%!test  # MATLAB parity: the importance of a bagged ensemble is the mean
%! load fisheriris
%! rng (1);
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 5);
%! I = cell2mat (cellfun (@(t) predictorImportance (t), Mdl.Trained, ...
%!                        'UniformOutput', false));
%! assert_equal (predictorImportance (Mdl), mean (I), 1e-15);

%!test  # MATLAB parity: out-of-bag predictions invert UseObsForLearner
%! load fisheriris
%! rng (2);
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 6);
%! U = ! Mdl.UseObsForLearner;
%! [l1, s1] = oobPredict (Mdl);
%! [l2, s2] = predict (Mdl, meas, 'UseObsForLearner', U);
%! assert_equal (l1, l2);
%! assert_equal (isequaln (s1, s2), true);
%! r = find (all (Mdl.UseObsForLearner, 2));
%! assert_equal (all (isnan (s1(r,:))(:)), true);

%!test  # MATLAB parity: the out-of-bag loss leaves out rows in every bag
%! load fisheriris
%! rng (2);
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 6);
%! [~, s] = oobPredict (Mdl);
%! have = ! any (isnan (s), 2);
%! g = grp2idx (species);
%! st = s(sub2ind (size (s), (1:150)', g));
%! so = s;
%! so(sub2ind (size (s), (1:150)', g)) = -Inf;
%! miss = st <= max (so, [], 2);
%! assert_equal (oobLoss (Mdl), mean (miss(have)), 1e-15);
%! U = ! Mdl.UseObsForLearner;
%! assert_equal (oobLoss (Mdl, 'Mode', 'cumulative'), ...
%!               loss (Mdl, meas, species, 'UseObsForLearner', U, ...
%!                     'Mode', 'cumulative'), 1e-15);
%! assert_equal (oobLoss (Mdl, 'Learners', [2, 4]), ...
%!               loss (Mdl, meas, species, 'UseObsForLearner', U, ...
%!                     'Learners', [2, 4]), 1e-15);

%!test  # MATLAB parity: the out-of-bag edge and margins
%! load fisheriris
%! rng (2);
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 6);
%! U = ! Mdl.UseObsForLearner;
%! assert_equal (oobEdge (Mdl), edge (Mdl, meas, species, ...
%!                                    'UseObsForLearner', U), 1e-15);
%! m = oobMargin (Mdl);
%! assert_equal (isequaln (m, margin (Mdl, meas, species, ...
%!                                    'UseObsForLearner', U)), true);
%! assert_equal (oobEdge (Mdl), mean (m(! isnan (m))), 1e-14);

%!test  # MATLAB parity: permuted importance of one tree is zero or infinite
%! load fisheriris
%! rng (4);
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 1);
%! imp = oobPermutedPredictorImportance (Mdl);
%! assert_equal (size (imp), [1, 4]);
%! assert_equal (all (imp == 0 | isinf (imp)), true);

%!test  # MATLAB parity: a constant predictor has zero permuted importance
%! load fisheriris
%! rng (5);
%! Mdl = ClassificationBaggedEnsemble ([meas, ones(150, 1)], species, ...
%!                                     'NumLearningCycles', 20);
%! imp = oobPermutedPredictorImportance (Mdl);
%! assert_equal (imp(5), 0);
%! assert_equal (size (oobPermutedPredictorImportance (Mdl, ...
%!                                                     'Learners', 1:5)), ...
%!               [1, 5]);

%!test  # MATLAB parity: permuting a petal measurement matters most
%! load fisheriris
%! rng (3);
%! Mdl = ClassificationBaggedEnsemble (meas, species, 'NumLearningCycles', 60);
%! imp = oobPermutedPredictorImportance (Mdl);
%! assert_equal (min (imp(3:4)) > max (imp(1:2)), true);

%!error<ClassificationBaggedEnsemble.oobPredict: invalid parameter name in optional pair arguments.> ...
%! load fisheriris
%! oobPredict (ClassificationBaggedEnsemble (meas, species, ...
%!                                           'NumLearningCycles', 1), ...
%!             'UseObsForLearner', true (150, 1))
%!error<ClassificationBaggedEnsemble.oobLoss: name-value arguments must be in pairs.> ...
%! load fisheriris
%! oobLoss (ClassificationBaggedEnsemble (meas, species, ...
%!                                        'NumLearningCycles', 1), 'Mode')
%!error<ClassificationBaggedEnsemble.oobLoss: invalid parameter name in optional pair arguments.> ...
%! load fisheriris
%! oobLoss (ClassificationBaggedEnsemble (meas, species, ...
%!                                        'NumLearningCycles', 1), ...
%!          'Weights', ones (150, 1))
%!error<ClassificationBaggedEnsemble.oobEdge: invalid parameter name in optional pair arguments.> ...
%! load fisheriris
%! oobEdge (ClassificationBaggedEnsemble (meas, species, ...
%!                                        'NumLearningCycles', 1), ...
%!          'LossFun', 'hinge')
%!error<ClassificationBaggedEnsemble.oobMargin: invalid parameter name in optional pair arguments.> ...
%! load fisheriris
%! oobMargin (ClassificationBaggedEnsemble (meas, species, ...
%!                                          'NumLearningCycles', 1), ...
%!            'Mode', 'cumulative')
%!error<ClassificationBaggedEnsemble.oobPermutedPredictorImportance: 'Learners' must be a vector of indices of trained learners.> ...
%! load fisheriris
%! oobPermutedPredictorImportance (ClassificationBaggedEnsemble (meas, ...
%!                                 species, 'NumLearningCycles', 1), ...
%!                                 'Learners', 2)
%!error<ClassificationBaggedEnsemble.oobPermutedPredictorImportance: 'Options' is not implemented.> ...
%! load fisheriris
%! oobPermutedPredictorImportance (ClassificationBaggedEnsemble (meas, ...
%!                                 species, 'NumLearningCycles', 1), ...
%!                                 'Options', struct ())

%!shared X2, Y2, S
%! load fisheriris
%! X2 = meas(51:150,:);
%! Y2 = species(51:150);
%! S = templateTree ('MaxNumSplits', 1);

%!test  # MATLAB parity: drawing every row without replacement is plain boosting
%! M = ClassificationBaggedEnsemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                                   'NumLearningCycles', 5, 'Learners', S, ...
%!                                   'FResample', 1, 'Replace', 'off');
%! assert_equal (M.Method, 'AdaBoostM1');
%! assert_equal (M.CombineWeights, 'WeightedSum');
%! assert_equal (M.UseObsForLearner, true (100, 5));
%! assert_equal (M.TrainedWeights', [1.37576765652097, 0.99353411077441, ...
%!                                   0.883434373403997, 0.554364192108759, ...
%!                                   0.268194447432049], 1e-12);

%!test  # MATLAB parity: a resampling boosting method draws with replacement
%! M = ClassificationBaggedEnsemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                                   'NumLearningCycles', 5, 'Learners', S, ...
%!                                   'Resample', 'on');
%! assert_equal ([M.FResample, M.Replace], [1, true]);
%! assert_equal (size (M.UseObsForLearner), [100, M.NumTrained]);
%! assert_equal (M.FitInfo * 100, round (M.FitInfo * 100), 1e-9);

%!test  # MATLAB parity: AdaBoostM1 reweights only the rows it drew
%! M = ClassificationBaggedEnsemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                                   'NumLearningCycles', 4, 'Learners', S, ...
%!                                   'FResample', 0.5, 'Replace', 'off');
%! y = 2 * strcmp (Y2, M.ClassNames{1}) - 1;
%! d = M.W / sum (M.W);
%! for t = 1:M.NumTrained
%!   u = M.UseObsForLearner(:,t);
%!   h = 2 * strcmp (predict (M.Trained{t}, X2), M.ClassNames{1}) - 1;
%!   e = sum (d(u) .* (h(u) != y(u))) / sum (d(u));
%!   assert_equal (M.FitInfo(t), e, 1e-12);
%!   s0 = sum (d(u));
%!   d(u) = d(u) .* exp (-M.TrainedWeights(t) * y(u) .* h(u));
%!   d(u) = d(u) / sum (d(u)) * s0;
%! endfor

%!test  # MATLAB parity: GentleBoost reweights only the rows it drew
%! M = ClassificationBaggedEnsemble (X2, Y2, 'Method', 'GentleBoost', ...
%!                                   'NumLearningCycles', 4, 'Learners', S, ...
%!                                   'FResample', 0.5, 'Replace', 'off');
%! y = 2 * strcmp (Y2, M.ClassNames{1}) - 1;
%! d = M.W / sum (M.W);
%! for t = 1:M.NumTrained
%!   u = M.UseObsForLearner(:,t);
%!   h = predict (M.Trained{t}, X2);
%!   assert_equal (M.FitInfo(t), ...
%!                 sum (d(u) .* (y(u) - h(u)) .^ 2) / sum (d(u)), 1e-12);
%!   s0 = sum (d(u));
%!   d(u) = d(u) .* exp (-y(u) .* h(u));
%!   d(u) = d(u) / sum (d(u)) * s0;
%! endfor

%!test  # MATLAB parity: LogitBoost moves the score of the rows it drew
%! M = ClassificationBaggedEnsemble (X2, Y2, 'Method', 'LogitBoost', ...
%!                                   'NumLearningCycles', 4, 'Learners', S, ...
%!                                   'FResample', 0.5, 'Replace', 'off');
%! y01 = double (strcmp (Y2, M.ClassNames{1}));
%! w0 = M.W / sum (M.W);
%! w = w0;
%! F = zeros (100, 1);
%! for t = 1:M.NumTrained
%!   u = M.UseObsForLearner(:,t);
%!   p = 1 ./ (1 + exp (-F));
%!   z = (y01 - p) ./ (p .* (1 - p));
%!   h = predict (M.Trained{t}, X2);
%!   assert_equal (M.FitInfo(t), ...
%!                 sum (w(u) .* (z(u) - h(u)) .^ 2) / sum (w(u)), 1e-10);
%!   F(u) += h(u) / 2;
%!   s0 = sum (w(u));
%!   pu = 1 ./ (1 + exp (-F(u)));
%!   w(u) = w0(u) .* pu .* (1 - pu);
%!   w(u) = w(u) / sum (w(u)) * s0;
%! endfor

%!test  # MATLAB parity: AdaBoostM2 keeps one weight per observation
%! load fisheriris
%! M = ClassificationBaggedEnsemble (meas, species, 'Method', 'AdaBoostM2', ...
%!                                   'NumLearningCycles', 3, 'Learners', S, ...
%!                                   'FResample', 0.5, 'Replace', 'off');
%! g = grp2idx (species);
%! tru = sub2ind ([150, 3], (1:150)', g);
%! w = M.W / sum (M.W);
%! for t = 1:M.NumTrained
%!   u = M.UseObsForLearner(:,t);
%!   [~, P] = predict (M.Trained{t}, meas);
%!   hy = P(tru);
%!   L = 1 - hy + P;
%!   L(tru) = 0;
%!   e = sum (w(u) .* sum (L(u,:), 2)) / sum (w(u)) / 4;
%!   assert_equal (M.FitInfo(t), e, 1e-12);
%!   E = exp (-M.TrainedWeights(t) * (1 + hy - P));
%!   E(tru) = 0;
%!   s0 = sum (w(u));
%!   w(u) = w(u) .* sum (E(u,:), 2);
%!   w(u) = w(u) / sum (w(u)) * s0;
%! endfor

%!test  # MATLAB parity: out-of-bag scores sum the learners that did not draw
%! M = ClassificationBaggedEnsemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                                   'NumLearningCycles', 8, 'Learners', S, ...
%!                                   'Resample', 'on');
%! [~, s] = oobPredict (M);
%! out = ! M.UseObsForLearner;
%! man = zeros (100, 1);
%! for t = 1:M.NumTrained
%!   h = 2 * strcmp (predict (M.Trained{t}, X2), M.ClassNames{1}) - 1;
%!   man += M.TrainedWeights(t) * h .* out(:,t);
%! endfor
%! k = any (out, 2);
%! assert_equal (s(k,1), man(k), 1e-12);

%!test  # resume grows the record of the rows each learner drew
%! M = ClassificationBaggedEnsemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                                   'NumLearningCycles', 2, 'Learners', S, ...
%!                                   'FResample', 0.7, 'Replace', 'off');
%! R = resume (M, 2);
%! assert_equal (size (R.UseObsForLearner), [100, R.NumTrained]);

%!test  # a resampled boosting ensemble cross-validates with its learning rate
%! M = ClassificationBaggedEnsemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                                   'NumLearningCycles', 2, 'Learners', S, ...
%!                                   'FResample', 0.7, 'LearnRate', 0.5);
%! CV = crossval (M, 'KFold', 3);
%! assert_equal (CV.Trainable{1}.LearnRate, 0.5);
%! assert_equal (class (CV.Trainable{1}), 'ClassificationBaggedEnsemble');
