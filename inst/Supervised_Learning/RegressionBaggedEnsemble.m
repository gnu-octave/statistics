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
  ## and which rows each tree drew.
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
      if (! any (cellfun (@(a) ischar (a) && strcmpi (a, 'Method'), ...
                          varargin(1:2:end))))
        varargin(end+1:end+2) = {'Method', 'Bag'};
      endif
      this = this@RegressionEnsemble (X, Y, varargin{:});
      this = bagProperties (this);

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
    ## @deftypefnx {RegressionBaggedEnsemble} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss of a bagged ensemble.
    ##
    ## Behaves as @code{CompactRegressionEnsemble.loss}.
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
%!error<RegressionBaggedEnsemble: 'Method' must be 'Bag'.> ...
%! load fisheriris
%! RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), 'Method', 'LSBoost')
%!error<RegressionBaggedEnsemble: 'LearnRate' cannot be used with the 'Bag' method.> ...
%! load fisheriris
%! RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), 'LearnRate', 0.5)
%!error<RegressionBaggedEnsemble.predict: too few input arguments.> ...
%! load fisheriris
%! predict (RegressionBaggedEnsemble (meas(:,2:4), meas(:,1), ...
%!                                    'NumLearningCycles', 1))
