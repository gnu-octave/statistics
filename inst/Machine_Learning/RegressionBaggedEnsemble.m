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

classdef RegressionBaggedEnsemble < RegressionEnsemble
  ## -*- texinfo -*-
  ## @deftp {statistics} RegressionBaggedEnsemble
  ##
  ## Bagged ensemble of regression trees
  ##
  ## A @code{RegressionBaggedEnsemble} object holds regression trees each
  ## grown on a sample drawn from the training data in proportion to the
  ## observation weights, each split chosen from @code{ceil (P / 3)}
  ## predictors drawn afresh at every node.  It predicts by averaging its
  ## trees.
  ##
  ## Create one with @code{fitrensemble} and @qcode{'Method'} set to
  ## @qcode{'Bag'}.  It carries everything a @code{RegressionEnsemble} does,
  ## and which rows each tree drew.  LSBoost that resamples, asked for with
  ## @qcode{'Resample'}, @qcode{'FResample'} or @qcode{'Replace'}, returns this
  ## class too: its trees are boosted and summed by their weights, each fitted
  ## on the rows it drew.
  ##
  ## @seealso{fitrensemble, RegressionEnsemble, CompactRegressionEnsemble,
  ## TreeBagger}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RegressionBaggedEnsemble} {property} FResample
    ##
    ## Share of the observations each tree draws
    ##
    ## A number greater than 0 and no greater than 1, each tree drawing
    ## @code{ceil (FResample * N)} observations.  This property is read-only.
    ##
    ## @end deftp
    FResample = 1;

    ## -*- texinfo -*-
    ## @deftp {RegressionBaggedEnsemble} {property} Replace
    ##
    ## Whether the trees draw with replacement
    ##
    ## A logical scalar, true by default.  This property is read-only.
    ##
    ## @end deftp
    Replace = true;

    ## -*- texinfo -*-
    ## @deftp {RegressionBaggedEnsemble} {property} UseObsForLearner
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
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{obj} =} RegressionBaggedEnsemble (@var{X}, @var{Y})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{obj} =} RegressionBaggedEnsemble (@dots{}, @var{name}, @var{value})
    ##
    ## Fit a bagged ensemble of regression trees.
    ##
    ## @code{fitrensemble} with @qcode{'Method'} set to @qcode{'Bag'} is the
    ## documented way in, and its help lists the options both take.
    ##
    ## @seealso{fitrensemble, RegressionEnsemble}
    ## @end deftypefn
    function this = RegressionBaggedEnsemble (X, Y, varargin)

      if (nargin < 2)
        error ("RegressionBaggedEnsemble: too few input arguments.");
      endif

      ## A table names its own predictors and says which hold levels
      [this, X, Y, varargin] = resolveTable (this, 'RegressionBaggedEnsemble', ...
                                             X, Y, varargin);
      if (! any (cellfun (@(a) ischar (a) && strcmpi (a, 'Method'), ...
                          varargin(1:2:end))))
        varargin(end+1:end+2) = {'Method', 'Bag'};
      endif
      this = this@RegressionEnsemble (X, Y, varargin{:});
      this = bagProperties (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate a bagged ensemble.
    ##
    ## Behaves as @code{RegressionEnsemble.crossval}, returning a
    ## @code{RegressionPartitionedEnsemble}.
    ##
    ## @seealso{RegressionBaggedEnsemble, RegressionPartitionedEnsemble}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)
      CVMdl = crossval@RegressionEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{B} =} regularize (@var{obj})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{B} =} regularize (@dots{}, @var{name}, @var{value})
    ##
    ## Find lasso weights for the trees of a bagged ensemble.
    ##
    ## Behaves as @code{RegressionEnsemble.regularize}.
    ##
    ## @seealso{RegressionBaggedEnsemble, RegressionEnsemble.regularize}
    ## @end deftypefn
    function this = regularize (this, varargin)
      this = regularize@RegressionEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{C} =} shrink (@var{obj})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{C} =} shrink (@dots{}, @var{name}, @var{value})
    ##
    ## Keep the trees of a bagged ensemble that a lasso weight retains.
    ##
    ## Behaves as @code{RegressionEnsemble.shrink}, returning a
    ## @code{CompactRegressionEnsemble}.
    ##
    ## @seealso{RegressionBaggedEnsemble, RegressionEnsemble.shrink}
    ## @end deftypefn
    function C = shrink (this, varargin)
      C = shrink@RegressionEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {[@var{vals}, @var{nlearn}] =} cvshrink (@var{obj})
    ## @deftypefnx {RegressionBaggedEnsemble} {[@var{vals}, @var{nlearn}] =} cvshrink (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate the shrinking of a bagged ensemble.
    ##
    ## Behaves as @code{RegressionEnsemble.cvshrink}.
    ##
    ## @seealso{RegressionBaggedEnsemble, RegressionEnsemble.cvshrink}
    ## @end deftypefn
    function [vals, nlearn] = cvshrink (this, varargin)
      [vals, nlearn] = cvshrink@RegressionEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionBaggedEnsemble} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Drop the training data from a bagged regression ensemble.
    ##
    ## Returns a @code{CompactRegressionEnsemble}, as
    ## @code{RegressionEnsemble.compact} does.
    ##
    ## @seealso{RegressionBaggedEnsemble, CompactRegressionEnsemble}
    ## @end deftypefn
    function CMdl = compact (this)
      CMdl = compact@RegressionEnsemble (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{B} =} resume (@var{obj}, @var{NumLearningCycles})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{B} =} resume (@dots{}, 'NPrint', @var{n})
    ##
    ## Grow more trees in a bagged regression ensemble.
    ##
    ## Behaves as @code{RegressionEnsemble.resume}, the new trees' samples
    ## added to @code{UseObsForLearner}.
    ##
    ## @seealso{RegressionBaggedEnsemble, RegressionEnsemble.resume}
    ## @end deftypefn
    function this = resume (this, varargin)
      this = resume@RegressionEnsemble (this, varargin{:});
      this = bagProperties (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{yfit} =} predict (@var{obj}, @var{X})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{yfit} =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict the response with a bagged regression ensemble.
    ##
    ## Behaves as @code{CompactRegressionEnsemble.predict}.
    ##
    ## @seealso{RegressionBaggedEnsemble, CompactRegressionEnsemble.predict}
    ## @end deftypefn
    function yfit = predict (this, varargin)
      yfit = predict@RegressionEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{L} =} loss (@var{obj}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{L} =} loss (@var{obj}, @var{Tbl})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss of a bagged ensemble.
    ##
    ## Behaves as @code{CompactRegressionEnsemble.loss}.
    ##
    ## @var{X} may also be a table @var{Tbl}, whose variables are matched to
    ## the predictors the model was fitted on by name and not by position.
    ## @code{loss (@var{obj}, @var{Tbl}, @var{ResponseVarName})} takes the
    ## response from the variable @var{ResponseVarName} names, and
    ## @code{loss (@var{obj}, @var{Tbl})} from the variable the model was
    ## fitted on.  The response may also be given beside the table as
    ## @var{Y}.
    ##
    ## @seealso{RegressionBaggedEnsemble, CompactRegressionEnsemble.loss}
    ## @end deftypefn
    function L = loss (this, varargin)
      L = loss@RegressionEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{yfit} =} resubPredict (@var{obj})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{yfit} =} resubPredict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict the response of the training data with a bagged ensemble.
    ##
    ## Behaves as @code{RegressionEnsemble.resubPredict}.
    ##
    ## @seealso{RegressionBaggedEnsemble, RegressionEnsemble.resubPredict}
    ## @end deftypefn
    function yfit = resubPredict (this, varargin)
      yfit = resubPredict@RegressionEnsemble (this, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss of a bagged ensemble on the training data.
    ##
    ## Behaves as @code{RegressionEnsemble.resubLoss}.
    ##
    ## @seealso{RegressionBaggedEnsemble, RegressionEnsemble.resubLoss}
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      L = resubLoss@RegressionEnsemble (this, varargin{:});
    endfunction


    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{imp} =} predictorImportance (@var{obj})
    ## @deftypefnx {RegressionBaggedEnsemble} {[@var{imp}, @var{ma}] =} predictorImportance (@var{obj})
    ##
    ## Estimate the importance of each predictor.
    ##
    ## The mean over the trees of each tree's @code{predictorImportance}, as
    ## @code{CompactRegressionEnsemble.predictorImportance} computes it.
    ##
    ## @seealso{RegressionBaggedEnsemble,
    ## RegressionBaggedEnsemble.oobPermutedPredictorImportance}
    ## @end deftypefn
    function [imp, ma] = predictorImportance (this)
      [imp, ma] = predictorImportance@RegressionEnsemble (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{yfit} =} oobPredict (@var{obj})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{yfit} =} oobPredict (@dots{}, 'Learners', @var{idx})
    ##
    ## Out-of-bag predictions for the training data.
    ##
    ## Each training observation is predicted by the trees whose samples left
    ## it out, as @code{predict} predicts it with @qcode{'UseObsForLearner'}
    ## set to @code{! UseObsForLearner}; one in the sample of every tree used
    ## is @code{NaN}.  @qcode{'Learners'} restricts the trees.
    ##
    ## @seealso{RegressionBaggedEnsemble, RegressionBaggedEnsemble.oobLoss}
    ## @end deftypefn
    function yfit = oobPredict (this, varargin)
      caller = 'RegressionBaggedEnsemble.oobPredict';
      args = oobArgs (varargin, {'Learners'}, caller);
      U = ! this.UseObsForLearner;
      yfit = ensemblePredict (compact (this), this.X, ...
                              [{'UseObsForLearner', U}, args], caller);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{L} =} oobLoss (@var{obj})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{L} =} oobLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-bag regression loss.
    ##
    ## The loss of the out-of-bag predictions against @code{Y}, weighted by
    ## @code{W}, an observation in the sample of every tree used being left
    ## out.  @qcode{'LossFun'} and @qcode{'Mode'} are taken as by
    ## @code{CompactRegressionEnsemble.loss}, and @qcode{'Learners'}
    ## restricts the trees.
    ##
    ## @seealso{RegressionBaggedEnsemble, RegressionBaggedEnsemble.oobPredict}
    ## @end deftypefn
    function L = oobLoss (this, varargin)
      caller = 'RegressionBaggedEnsemble.oobLoss';
      args = oobArgs (varargin, {'Learners', 'LossFun', 'Mode'}, caller);
      L = ensembleLoss (compact (this), this.X, this.Y, ...
                        [{'UseObsForLearner', ! this.UseObsForLearner, ...
                          'Weights', this.W}, args], caller);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionBaggedEnsemble} {@var{imp} =} oobPermutedPredictorImportance (@var{obj})
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{imp} =} oobPermutedPredictorImportance (@dots{}, 'Learners', @var{idx})
    ##
    ## Out-of-bag predictor importance by permutation.
    ##
    ## For each tree, the values of each predictor are permuted among the
    ## observations out of its bag, and the tree's mean squared error on them,
    ## weighted by @code{W}, is taken before and after.  @var{imp} holds, for
    ## each predictor, the mean of the rise over the trees divided by its
    ## standard deviation over the trees, zero where the mean is zero.
    ## @qcode{'Learners'} restricts the trees.
    ##
    ## @seealso{RegressionBaggedEnsemble,
    ## RegressionBaggedEnsemble.predictorImportance}
    ## @end deftypefn
    function imp = oobPermutedPredictorImportance (this, varargin)

      caller = 'RegressionBaggedEnsemble.oobPermutedPredictorImportance';
      if (mod (numel (varargin), 2) != 0)
        error ("%s: name-value arguments must be in pairs.", caller);
      endif

      ## Parse optional paired arguments; empty 'Learners' stands for every
      ## trained learner
      [Learners, Options, args] = ...
             parsePairedArguments ({'Learners', 'Options'}, {[], []}, ...
                                   varargin(:));

      ## Validate optional paired arguments
      if (isempty (Learners))
        learners = 1:this.NumTrained;
      else
        learners = oobLearners (Learners, this.NumTrained, caller);
      endif
      if (! isempty (Options))
        error ("%s: 'Options' is not implemented.", caller);
      endif

      if (! isempty (args))
        error ("%s: invalid optional paired argument.", caller);
      endif
      p = columns (this.X);
      D = zeros (numel (learners), p);
      for j = 1:numel (learners)
        t = learners(j);
        r = find (! this.UseObsForLearner(:,t));
        w = double (this.W(r));
        if (isempty (r) || ! (sum (w) > 0))
          continue;
        endif
        w /= sum (w);
        tree = this.Trained{t};
        Xo = this.X(r,:);
        yo = this.Y(r);
        e0 = sum (w .* (predict (tree, Xo) - yo) .^ 2);
        for v = 1:p
          Xp = Xo;
          Xp(:,v) = Xo(randperm (numel (r)), v);
          D(j,v) = sum (w .* (predict (tree, Xp) - yo) .^ 2) - e0;
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
      error ("%s: invalid optional paired argument.", caller);
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
%!test  # MATLAB parity: the properties of a bagged regression ensemble
%! load fisheriris
%! Mdl = RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                 'NumLearningCycles', 3);
%! assert_equal (numel (properties (Mdl)), 27);
%! assert_equal (Mdl.Method, 'Bag');
%! assert_equal (Mdl.TrainedWeights, [1; 1; 1]);
%! assert_equal (Mdl.FitInfo, []);
%! assert_equal (Mdl.FitInfoDescription, 'None');
%! assert_equal (Mdl.FResample, 1);
%! assert_equal (Mdl.Replace, true);
%! assert_equal (size (Mdl.UseObsForLearner), [150, 3]);
%! assert_equal (isfield (Mdl.ModelParameters, 'LearnRate'), false);

%!test  # MATLAB parity: the prediction is the mean of the trees
%! load fisheriris
%! rng (5);
%! Mdl = RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                 'NumLearningCycles', 4);
%! P = zeros (150, 1);
%! for t = 1:4
%!   P += predict (Mdl.Trained{t}, meas(:,2:4));
%! endfor
%! assert_equal (predict (Mdl, meas(:,2:4)), P / 4, 1e-14);
%! assert_equal (min (cellfun (@(t) min (t.NodeSize), Mdl.Trained)) >= 5, true);

%!test  # MATLAB parity: each tree draws FResample of the rows
%! load fisheriris
%! rng (5);
%! Mdl = RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                 'NumLearningCycles', 2, ...
%!                                 'FResample', 0.5, 'Replace', 'off');
%! assert_equal (sum (Mdl.UseObsForLearner), [75, 75]);
%! assert_equal (Mdl.Trained{1}.NodeSize(1), 75);

%!test  # MATLAB parity: the bootstrap draws in proportion to the weights
%! load fisheriris
%! rng (4);
%! w = [10 * ones(50, 1); ones(100, 1)];
%! Mdl = RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                 'NumLearningCycles', 20, 'Weights', w);
%! assert_equal ([sum(Mdl.W(1:50)), sum(Mdl.W(51:150))], [5, 1] / 6, 1e-15);
%! share = mean (sum (Mdl.UseObsForLearner(1:50,:)) / 50);
%! assert_equal (share > 0.85, true);

%!test  # MATLAB parity: a row no tree may predict is NaN
%! load fisheriris
%! rng (5);
%! Mdl = RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                 'NumLearningCycles', 4);
%! U = false (2, 4);
%! U(2,1) = true;
%! yf = predict (Mdl, meas(1:2,2:4), 'UseObsForLearner', U);
%! assert_equal (isnan (yf(1)), true);
%! assert_equal (yf(2), predict (Mdl.Trained{1}, meas(2,2:4)), 1e-15);

%!test  # resuming adds the new trees' samples
%! load fisheriris
%! Mdl = RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                 'NumLearningCycles', 2);
%! Mdl = resume (Mdl, 2);
%! assert_equal (class (Mdl), 'RegressionBaggedEnsemble');
%! assert_equal (Mdl.NumTrained, 4);
%! assert_equal (size (Mdl.UseObsForLearner), [150, 4]);
%! assert_equal (class (compact (Mdl)), 'CompactRegressionEnsemble');

## Test input validation
%!error<RegressionBaggedEnsemble: too few input arguments.> ...
%! RegressionBaggedEnsemble (1)
%!error<RegressionBaggedEnsemble: 'Method' must be 'Bag' unless the ensemble resamples.> ...
%! load fisheriris
%! RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), 'Method', 'LSBoost')
%!error<RegressionBaggedEnsemble: 'LearnRate' cannot be used with the 'Bag' method.> ...
%! load fisheriris
%! RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), 'LearnRate', 0.5)
%!error<RegressionBaggedEnsemble.predict: too few input arguments.> ...
%! load fisheriris
%! predict (RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                    'NumLearningCycles', 1))

%!test  # MATLAB parity: out-of-bag predictions invert UseObsForLearner
%! load fisheriris
%! rng (5);
%! Mdl = RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                 'NumLearningCycles', 6);
%! U = ! Mdl.UseObsForLearner;
%! assert_equal (isequaln (oobPredict (Mdl), ...
%!                         predict (Mdl, meas(:,2:4), ...
%!                                  'UseObsForLearner', U)), true);
%! assert_equal (oobLoss (Mdl), loss (Mdl, meas(:,2:4), meas(:,1), ...
%!                                    'UseObsForLearner', U), 1e-15);
%! assert_equal (oobLoss (Mdl, 'Mode', 'cumulative'), ...
%!               loss (Mdl, meas(:,2:4), meas(:,1), ...
%!                     'UseObsForLearner', U, 'Mode', 'cumulative'), 1e-14);

%!test  # MATLAB parity: the out-of-bag loss leaves out rows in every bag
%! load fisheriris
%! rng (5);
%! Mdl = RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                 'NumLearningCycles', 6);
%! yf = oobPredict (Mdl);
%! have = ! isnan (yf);
%! assert_equal (oobLoss (Mdl), mean ((yf(have) - meas(have,1)) .^ 2), 1e-14);

%!test  # MATLAB parity: the importance of a bagged regression ensemble
%! load fisheriris
%! rng (6);
%! Mdl = RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                 'NumLearningCycles', 5);
%! I = cell2mat (cellfun (@(t) predictorImportance (t), Mdl.Trained, ...
%!                        'UniformOutput', false));
%! assert_equal (predictorImportance (Mdl), mean (I), 1e-15);

%!test  # MATLAB parity: a constant predictor has zero permuted importance
%! load fisheriris
%! rng (7);
%! Mdl = RegressionBaggedEnsemble ([meas(:,2:4), ones(150, 1)], meas(:,1), ...
%!                                 'NumLearningCycles', 20);
%! imp = oobPermutedPredictorImportance (Mdl);
%! assert_equal (size (imp), [1, 4]);
%! assert_equal (imp(4), 0);
%! assert_equal (all (imp(1:3) > 0), true);

%!error<RegressionBaggedEnsemble.oobPredict: invalid optional paired argument.> ...
%! load fisheriris
%! oobPredict (RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                       'NumLearningCycles', 1), ...
%!             'Mode', 'ensemble')
%!error<RegressionBaggedEnsemble.oobLoss: invalid optional paired argument.> ...
%! load fisheriris
%! oobLoss (RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                    'NumLearningCycles', 1), ...
%!          'Weights', ones (150, 1))
%!error<RegressionBaggedEnsemble.oobPermutedPredictorImportance: 'Options' is not implemented.> ...
%! load fisheriris
%! oobPermutedPredictorImportance (RegressionBaggedEnsemble (meas(:,2:4), ...
%!                                 meas(:,1), 'NumLearningCycles', 1), ...
%!                                 'Options', struct ())
%!error<RegressionBaggedEnsemble.oobPermutedPredictorImportance: invalid optional paired argument.> ...
%! load fisheriris
%! oobPermutedPredictorImportance (RegressionBaggedEnsemble (meas(:,2:4), ...
%!                                 meas(:,1), 'NumLearningCycles', 1), ...
%!                                 'Bogus', struct ())

%!shared Xr, yr, tr
%! load fisheriris
%! Xr = meas(:,2:4);
%! yr = meas(:,1);
%! tr = templateTree ('MaxNumSplits', 3);

%!test  # MATLAB parity: LSBoost drawing every row without replacement is plain
%! M = RegressionBaggedEnsemble (Xr, yr, 'Method', 'LSBoost', ...
%!                               'NumLearningCycles', 4, 'Learners', tr, ...
%!                               'FResample', 1, 'Replace', 'off');
%! P = fitrensemble (Xr, yr, 'NumLearningCycles', 4, 'Learners', tr);
%! assert_equal (M.FitInfo, P.FitInfo, 1e-12);
%! assert_equal (predict (M, Xr(1:5,:)), predict (P, Xr(1:5,:)), 1e-12);
%! assert_equal ([M.FResample, M.Replace], [1, false]);

%!test  # MATLAB parity: resampled LSBoost measures its fit over every row
%! M = RegressionBaggedEnsemble (Xr, yr, 'Method', 'LSBoost', ...
%!                               'NumLearningCycles', 4, 'Learners', tr, ...
%!                               'Resample', 'on');
%! assert_equal (size (M.UseObsForLearner), [150, 4]);
%! W = M.W / sum (M.W);
%! F = zeros (150, 1);
%! for t = 1:4
%!   h = predict (M.Trained{t}, Xr);
%!   assert_equal (M.FitInfo(t), sum (W .* (yr - F - h) .^ 2), 1e-12);
%!   F += M.TrainedWeights(t) * h;
%! endfor

%!test  # a resampled LSBoost cross-validates with its learning rate
%! M = RegressionBaggedEnsemble (Xr, yr, 'Method', 'LSBoost', ...
%!                               'NumLearningCycles', 2, 'Learners', tr, ...
%!                               'FResample', 0.7, 'LearnRate', 0.5);
%! CV = crossval (M, 'KFold', 3);
%! assert_equal (CV.Trainable{1}.LearnRate, 0.5);

## A table at loss
%!test  # the response is named, left out, or given beside the table
%! load fisheriris
%! X = meas(:,2:3);
%! y = meas(:,1);
%! T = table (X(:,1), X(:,2), 'VariableNames', {'SW', 'PL'});
%! T.SL = y;
%! Mdl = fitrensemble (T, 'SL', 'Method', 'Bag');
%! a = loss (Mdl, X, y);
%! assert_equal (loss (Mdl, T(:,1:2), y), a);
%! assert_equal (loss (Mdl, T, 'SL'), a);
%! assert_equal (loss (Mdl, T), a);
