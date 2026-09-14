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
  ## does, and which rows each tree drew.
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
%!error<ClassificationBaggedEnsemble: 'Method' must be 'Bag'.> ...
%! load fisheriris
%! ClassificationBaggedEnsemble (meas, species, 'Method', 'AdaBoostM2')
%!error<ClassificationBaggedEnsemble: 'LearnRate' cannot be used with the 'Bag' method.> ...
%! load fisheriris
%! ClassificationBaggedEnsemble (meas, species, 'LearnRate', 0.5)
%!error<ClassificationBaggedEnsemble.predict: too few input arguments.> ...
%! load fisheriris
%! predict (ClassificationBaggedEnsemble (meas, species, ...
%!                                        'NumLearningCycles', 1))
