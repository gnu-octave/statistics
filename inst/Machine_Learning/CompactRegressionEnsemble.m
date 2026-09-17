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

classdef CompactRegressionEnsemble < PredictiveModel
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactRegressionEnsemble
  ##
  ## Compact ensemble of regression trees
  ##
  ## A @code{CompactRegressionEnsemble} object carries the trained trees of a
  ## boosted or bagged regression ensemble and what prediction needs, but not
  ## the observations it was fitted on.  It predicts new data identically to
  ## the ensemble it came from, and trees can be removed from it.
  ##
  ## Create one with the @code{compact} method of a
  ## @code{RegressionEnsemble} or @code{RegressionBaggedEnsemble} object.
  ##
  ## @seealso{fitrensemble, RegressionEnsemble, RegressionBaggedEnsemble}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionEnsemble} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionEnsemble} {property} CategoricalPredictors
    ##
    ## Indices of categorical predictors
    ##
    ## The predictors every tree treats as categorical, empty when none
    ## is.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionEnsemble} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName = 'Y';

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionEnsemble} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the learners saw them
    ##
    ## The same as @code{PredictorNames}.  This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionEnsemble} {property} UsePredForLearner
    ##
    ## Which predictors each learner uses
    ##
    ## Always empty, as MATLAB returns it for tree learners.  This property
    ## is read-only.
    ##
    ## @end deftp
    UsePredForLearner = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionEnsemble} {property} NumTrained
    ##
    ## Number of trained trees
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumTrained = 0;

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionEnsemble} {property} Trained
    ##
    ## Trained trees
    ##
    ## A column cell array of @code{CompactRegressionTree} objects.  This
    ## property is read-only.
    ##
    ## @end deftp
    Trained = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionEnsemble} {property} TrainedWeights
    ##
    ## Weights of the trained trees
    ##
    ## A column with one weight per tree.  This property is read-only.
    ##
    ## @end deftp
    TrainedWeights = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionEnsemble} {property} CombineWeights
    ##
    ## How the trees are combined
    ##
    ## @qcode{'WeightedSum'} for LSBoost, whose prediction is the sum of each
    ## tree's prediction times its weight, or @qcode{'WeightedAverage'} for
    ## Bag, whose prediction is the weighted average of its trees'.  This
    ## property is read-only.
    ##
    ## @end deftp
    CombineWeights = 'WeightedSum';

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionEnsemble} {property} ResponseTransform
    ##
    ## Transform applied to the predicted response
    ##
    ## @qcode{'none'} (default), @qcode{'exp'}, @qcode{'log'} or a function
    ## handle.  Predictions and losses use the transformed response.
    ##
    ## @end deftp
    ResponseTransform = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    RTfun = @(y) y;      # the response transform as a function
  endproperties

  methods (Hidden)

    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = ...
        parseResponseTransform (val, 'CompactRegressionEnsemble');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactRegressionEnsemble} {@var{obj} =} CompactRegressionEnsemble (@var{Mdl})
    ##
    ## Create a @code{CompactRegressionEnsemble} object.
    ##
    ## @var{Mdl} is the @code{RegressionEnsemble} or
    ## @code{RegressionBaggedEnsemble} object to compact.  The documented way
    ## to reach this constructor is the @code{compact} method.
    ##
    ## @end deftypefn
    function this = CompactRegressionEnsemble (Mdl = [])

      if (isempty (Mdl))
        return;
      endif
      if (! isa (Mdl, 'RegressionEnsemble'))
        error (strcat ("CompactRegressionEnsemble: MDL must be a", ...
                       " 'RegressionEnsemble' object."));
      endif

      this.PredictorNames         = Mdl.PredictorNames;
      this.CategoricalPredictors  = Mdl.CategoricalPredictors;
      this.ResponseName           = Mdl.ResponseName;
      this.ExpandedPredictorNames = Mdl.ExpandedPredictorNames;
      this.UsePredForLearner      = Mdl.UsePredForLearner;
      this.NumTrained             = Mdl.NumTrained;
      this.Trained                = Mdl.Trained;
      this.TrainedWeights         = Mdl.TrainedWeights;
      this.CombineWeights         = Mdl.CombineWeights;
      ## A transform given as a function handle is held as its text, which
      ## names no built-in transform, so the function itself is copied.
      if (any (strcmp (Mdl.ResponseTransform, {'none', 'exp', 'log'})))
        this.ResponseTransform    = Mdl.ResponseTransform;
      else
        this.ResponseTransform    = Mdl.RTfun;
      endif

    endfunction

    ## Keep the learners IDX, in that order, with the weights WEIGHTS and the
    ## combination rule COMBINE; shrink builds its ensemble with this.
    function this = keepLearners (this, idx, weights, combine)

      this.Trained = this.Trained(idx);
      this.TrainedWeights = weights(:);
      this.NumTrained = numel (idx);
      this.CombineWeights = combine;

    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ("\n  CompactRegressionEnsemble\n\n");
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %s\n", 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ("%+25s: '%s'\n", 'ResponseTransform', this.ResponseTransform);
      fprintf ("%+25s: %d\n", 'NumTrained', this.NumTrained);
      fprintf ("\n");
    endfunction

    ## The shared bodies of predict and loss.  The full and bagged classes
    ## call them with their own names, so an error names the method called.

    function yfit = ensemblePredict (this, X, args, caller)

      o = ensembleArgs (this, X, args, {'Learners', 'UseObsForLearner'}, ...
                        caller);
      yfit = this.RTfun (ensembleResponse (this, X, o, 'ensemble'));

    endfunction

    function L = ensembleLoss (this, X, Y, args, caller)

      o = ensembleArgs (this, X, args, {'LossFun', 'Learners', 'Mode', ...
                                        'UseObsForLearner', 'Weights'}, ...
                        caller);
      if (! (isnumeric (Y) && isreal (Y) && isvector (Y)))
        error ("%s: Y must be a real numeric vector.", caller);
      endif
      Y = double (Y(:));
      if (numel (Y) != rows (X))
        error ("%s: X and Y must have the same number of rows.", caller);
      endif
      w = o.Weights;
      if (isempty (w))
        w = ones (rows (X), 1);
      endif
      ## A row with no response cannot be judged, so it is left out and the
      ## weights are normalized over the rest, as MATLAB does.
      keep = ! isnan (Y);
      if (! (sum (w(keep)) > 0))
        error (strcat ("%s: 'Weights' must not be zero for every", ...
                       " observation with a response."), caller);
      endif
      o.U = o.U(keep,:);
      Yf = this.RTfun (ensembleResponse (this, X(keep,:), o, o.Mode));
      Y = Y(keep);
      w = w(keep) / sum (w(keep));
      L = zeros (columns (Yf), 1);
      for k = 1:columns (Yf)
        ## A row no tree may predict is left out with its weight, the rest
        ## renormalized, as MATLAB does.
        have = ! isnan (Yf(:,k));
        wk = w(have);
        if (! (sum (wk) > 0))
          L(k) = NaN;
          continue;
        endif
        wk /= sum (wk);
        if (is_function_handle (o.LossFun))
          L(k) = o.LossFun (Y(have), Yf(have,k), wk);
        else
          L(k) = sum (wk .* (Yf(have,k) - Y(have)) .^ 2);
        endif
      endfor

    endfunction

    ## The untransformed predictions for the rows of X over the first t
    ## trees, for every t, as an NxNumTrained matrix.
    function Y = ensembleSteps (this, X)

      T = this.NumTrained;
      o = struct ('Learners', 1:T, 'U', true (rows (X), T));
      Y = ensembleResponse (this, X, o, 'cumulative');

    endfunction

    function [imp, ma] = ensembleImportance (this)

      imp = zeros (1, numel (this.PredictorNames));
      for t = 1:this.NumTrained
        imp += this.TrainedWeights(t) * predictorImportance (this.Trained{t});
      endfor
      if (sum (this.TrainedWeights) > 0)
        imp /= sum (this.TrainedWeights);
      endif
      ma = [];

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionEnsemble} {@var{yfit} =} predict (@var{obj}, @var{X})
    ## @deftypefnx {CompactRegressionEnsemble} {@var{yfit} =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict the response with a compact regression ensemble.
    ##
    ## @var{yfit} holds, for each row of @var{X}, the sum over the trees of
    ## each tree's prediction times its weight for LSBoost, or the weighted
    ## average of the trees' predictions for Bag, after
    ## @code{ResponseTransform}.  A row that no tree may predict is
    ## @code{NaN}.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Learners'} @tab @tab A vector of indices of the trees to
    ## use.  The default is all of them.
    ## @item @qcode{'UseObsForLearner'} @tab @tab An @math{NxNumTrained}
    ## logical matrix saying which tree may predict which row.  The default
    ## lets every tree predict every row.
    ## @end multitable
    ##
    ## @seealso{CompactRegressionEnsemble, fitrensemble}
    ## @end deftypefn
    function yfit = predict (this, X, varargin)

      if (nargin < 2)
        error ("CompactRegressionEnsemble.predict: too few input arguments.");
      endif
      yfit = ensemblePredict (this, X, varargin, ...
                              'CompactRegressionEnsemble.predict');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionEnsemble} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactRegressionEnsemble} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss of a compact ensemble.
    ##
    ## @var{L} is the weighted mean squared error of the predictions for the
    ## rows of @var{X} against @var{Y}, the weights normalized to sum to one
    ## over the rows that have a response; a row whose response is missing is
    ## left out.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'LossFun'} @tab @tab @qcode{'mse'} (default) or a
    ## function handle called as @code{lossfun (Y, Yfit, W)}, with column
    ## vectors of the responses, the predictions and the normalized weights,
    ## returning a scalar.
    ## @item @qcode{'Mode'} @tab @tab @qcode{'ensemble'} (default) for one
    ## loss over the trees used, @qcode{'cumulative'} for a column whose
    ## element @math{j} uses the first @math{j} of them, or
    ## @qcode{'individual'} for a column with the loss of each on its own.
    ## @item @qcode{'Weights'} @tab @tab A nonnegative vector with one weight
    ## per row.  The default is uniform.
    ## @end multitable
    ##
    ## @qcode{'Learners'} and @qcode{'UseObsForLearner'} are taken as by
    ## @code{predict}.  A row that no tree may predict is left out and the
    ## weights are renormalized over the rest.
    ##
    ## @seealso{CompactRegressionEnsemble, CompactRegressionEnsemble.predict}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      if (nargin < 3)
        error ("CompactRegressionEnsemble.loss: too few input arguments.");
      endif
      L = ensembleLoss (this, X, Y, varargin, ...
                        'CompactRegressionEnsemble.loss');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionEnsemble} {@var{imp} =} predictorImportance (@var{obj})
    ## @deftypefnx {CompactRegressionEnsemble} {[@var{imp}, @var{ma}] =} predictorImportance (@var{obj})
    ##
    ## Estimate the importance of each predictor.
    ##
    ## @var{imp} is a row vector with one element per predictor, the average
    ## over the trees of each tree's @code{predictorImportance}, weighted by
    ## @code{TrainedWeights}.  @var{ma}, the predictive measure of
    ## association between the predictors, is empty, the trees growing no
    ## surrogate splits.
    ##
    ## @seealso{CompactRegressionEnsemble}
    ## @end deftypefn
    function [imp, ma] = predictorImportance (this)

      [imp, ma] = ensembleImportance (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactRegressionEnsemble} {@var{C} =} removeLearners (@var{obj}, @var{idx})
    ##
    ## Remove trees from a compact regression ensemble.
    ##
    ## @var{C} is the ensemble without the trees whose indices @var{idx}
    ## holds, their weights removed with them.
    ##
    ## @seealso{CompactRegressionEnsemble}
    ## @end deftypefn
    function this = removeLearners (this, idx)

      if (nargin < 2)
        error (strcat ("CompactRegressionEnsemble.removeLearners: too few", ...
                       " input arguments."));
      endif
      if (! (isnumeric (idx) && isvector (idx) && isreal (idx)
             && all (idx >= 1) && all (idx <= this.NumTrained)
             && all (idx == fix (idx))))
        error (strcat ("CompactRegressionEnsemble.removeLearners: IDX must", ...
                       " be a vector of indices of trained learners."));
      endif
      keep = true (1, this.NumTrained);
      keep(idx) = false;
      this.Trained = this.Trained(keep);
      this.TrainedWeights = this.TrainedWeights(keep);
      this.NumTrained = sum (keep);

    endfunction

  endmethods

  methods (Access = private)

    ## The Name-Value arguments of predict and loss, those not in ALLOWED
    ## refused.
    function o = ensembleArgs (this, X, args, allowed, caller)

      if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
        error ("%s: X must be a real numeric matrix.", caller);
      endif
      if (columns (X) != numel (this.PredictorNames))
        error ("%s: X must have one column per predictor.", caller);
      endif
      if (mod (numel (args), 2) != 0)
        error ("%s: name-value arguments must be in pairs.", caller);
      endif
      T = this.NumTrained;
      o = struct ('Learners', 1:T, 'U', true (rows (X), T), ...
                  'Mode', 'ensemble', 'Weights', [], 'LossFun', 'mse');
      for i = 1:2:numel (args)
        name = args{i};
        val = args{i+1};
        if (! (ischar (name) && any (strcmpi (name, allowed))))
          error ("%s: invalid parameter name in optional pair arguments.", ...
                 caller);
        endif
        switch (tolower (name))
          case 'learners'
            if (! (isnumeric (val) && isvector (val) && isreal (val)
                   && all (val >= 1) && all (val <= T)
                   && all (val == fix (val))))
              error (strcat ("%s: 'Learners' must be a vector of indices", ...
                             " of trained learners."), caller);
            endif
            o.Learners = double (val(:)');
          case 'useobsforlearner'
            if (! (islogical (val) && isequal (size (val), [rows(X), T])))
              error (strcat ("%s: 'UseObsForLearner' must be a logical", ...
                             " matrix with one row per observation and", ...
                             " one column per trained learner."), caller);
            endif
            o.U = val;
          case 'mode'
            if (! (ischar (val) && any (strcmpi (val, {'ensemble', ...
                                                     'cumulative', ...
                                                     'individual'}))))
              error (strcat ("%s: 'Mode' must be 'ensemble', 'cumulative'", ...
                             " or 'individual'."), caller);
            endif
            o.Mode = tolower (val);
          case 'weights'
            if (! (isnumeric (val) && isvector (val) && isreal (val)
                   && numel (val) == rows (X) && all (val >= 0)))
              error (strcat ("%s: 'Weights' must be a nonnegative numeric", ...
                             " vector with one element per observation."), ...
                     caller);
            endif
            o.Weights = double (val(:));
          case 'lossfun'
            if (ischar (val) && strcmpi (val, 'mse'))
              o.LossFun = 'mse';
            elseif (is_function_handle (val))
              o.LossFun = val;
            else
              error ("%s: 'LossFun' must be 'mse' or a function handle.", ...
                     caller);
            endif
        endswitch
      endfor

    endfunction

    ## The untransformed predictions over the trees O names: a column for
    ## MODE 'ensemble', one column per tree otherwise.  A row no tree may
    ## predict is NaN.
    function Y = ensembleResponse (this, X, o, mode)

      n = rows (X);
      T = numel (o.Learners);
      A = zeros (n, T);
      V = zeros (n, T);
      for j = 1:T
        t = o.Learners(j);
        u = o.U(:,t);
        A(:,j) = this.TrainedWeights(t) * predict (this.Trained{t}, X) .* u;
        V(:,j) = this.TrainedWeights(t) * u;
      endfor
      U = double (o.U(:,o.Learners));
      switch (mode)
        case 'ensemble'
          num = sum (A, 2);
          den = sum (V, 2);
          used = sum (U, 2) > 0;
        case 'cumulative'
          num = cumsum (A, 2);
          den = cumsum (V, 2);
          used = cumsum (U, 2) > 0;
        otherwise
          num = A;
          den = V;
          used = U > 0;
      endswitch
      if (strcmp (this.CombineWeights, 'WeightedAverage'))
        Y = num ./ den;
      else
        Y = num;
      endif
      Y(! used) = NaN;

    endfunction

  endmethods

endclassdef

## Test output
%!shared X, y, C
%! load fisheriris
%! X = meas(:,2:4);
%! y = meas(:,1);
%! C = compact (fitrensemble (X, y, 'NumLearningCycles', 4, ...
%!                            'Learners', templateTree ('MaxNumSplits', 1)));

%!test  # MATLAB parity: the properties of a compact regression ensemble
%! assert_equal (numel (properties (C)), 10);
%! assert_equal (C.NumTrained, 4);
%! assert_equal (isprop (C, 'X'), false);

%!test  # MATLAB parity: predictions over a subset of the trees
%! assert_equal (predict (C, X(1:2,:), 'Learners', [2, 4]), ...
%!               [0.061599937710460; 0.061599937710460], 1e-13);

%!test  # MATLAB parity: a row no tree may predict is NaN
%! U = true (2, 4);
%! U(1,:) = false;
%! U(2,[1, 3]) = false;
%! yf = predict (C, X(1:2,:), 'UseObsForLearner', U);
%! assert_equal (isnan (yf(1)), true);
%! assert_equal (yf(2), 0.061599937710460, 1e-13);

%!test  # MATLAB parity: the loss in its three modes
%! assert_equal (loss (C, X, y), 0.159714716160712, 1e-13);
%! assert_equal (loss (C, X, y, 'Mode', 'cumulative'), ...
%!               [0.263279369032794; 0.185334478476685; ...
%!                0.176267398277276; 0.159714716160713], 1e-13);
%! assert_equal (loss (C, X, y, 'Mode', 'individual'), ...
%!               [0.263279369032794; 34.658932998586877; ...
%!                34.730819835683143; 35.008549381231070], 1e-11);
%! assert_equal (loss (C, X, y, 'Learners', [1, 3]), 0.245973792973264, 1e-13);

%!test  # MATLAB parity: weights are normalized to sum to one
%! w = [2 * ones(50, 1); ones(100, 1)];
%! assert_equal (loss (C, X, y, 'Weights', w), 0.150822300473683, 1e-13);

%!test  # MATLAB parity: a custom loss gets normalized weights
%! f = @(Y, Yf, W) sum (W .* abs (Y - Yf)) / sum (W);
%! assert_equal (loss (C, X, y, 'LossFun', f), 0.321254856340782, 1e-13);
%! w = [3 * ones(50, 1); ones(100, 1)];
%! assert_equal (loss (C, X, y, 'LossFun', @(Y, Yf, W) sum (W), ...
%!                     'Weights', w), 1, 1e-14);
%! assert_equal (loss (C, X, y, 'LossFun', @(Y, Yf, W) W(1), ...
%!                     'Weights', w), 0.012, 1e-15);

%!test  # MATLAB parity: a missing response is left out of the loss
%! C3 = removeLearners (C, 4);
%! yn = y;
%! yn(3) = NaN;
%! assert_equal (loss (C3, X, yn), 0.177018571063146, 1e-13);
%! assert_equal (loss (C3, X, y, 'Mode', 'cumulative', 'Weights', ...
%!                     [3 * ones(50, 1); ones(100, 1)]), ...
%!               [0.218707467544654; 0.164143043739940; ...
%!                0.149902130360542], 1e-13);

%!test  # MATLAB parity: removing trees
%! D = removeLearners (C, [1, 3]);
%! assert_equal (D.NumTrained, 2);
%! assert_equal (predict (D, X(1:2,:)), [0.061599937710460; ...
%!                                       0.061599937710460], 1e-13);

## Test input validation
%!error<CompactRegressionEnsemble: MDL must be a 'RegressionEnsemble' object.> ...
%! CompactRegressionEnsemble (1)
%!error<CompactRegressionEnsemble.predict: too few input arguments.> ...
%! predict (C)
%!error<CompactRegressionEnsemble.predict: X must be a real numeric matrix.> ...
%! predict (C, {1})
%!error<CompactRegressionEnsemble.predict: X must have one column per predictor.> ...
%! predict (C, ones (2, 4))
%!error<CompactRegressionEnsemble.predict: name-value arguments must be in pairs.> ...
%! predict (C, X, 'Learners')
%!error<CompactRegressionEnsemble.predict: invalid parameter name in optional pair arguments.> ...
%! predict (C, X, 'Mode', 'ensemble')
%!error<CompactRegressionEnsemble.predict: 'Learners' must be a vector of indices of trained learners.> ...
%! predict (C, X, 'Learners', 5)
%!error<CompactRegressionEnsemble.predict: 'UseObsForLearner' must be a logical matrix with one row per observation and one column per trained learner.> ...
%! predict (C, X, 'UseObsForLearner', true (2, 4))
%!error<CompactRegressionEnsemble.loss: too few input arguments.> ...
%! loss (C, X)
%!error<CompactRegressionEnsemble.loss: 'Mode' must be 'ensemble', 'cumulative' or 'individual'.> ...
%! loss (C, X, y, 'Mode', 'all')
%!error<CompactRegressionEnsemble.loss: 'Weights' must be a nonnegative numeric vector with one element per observation.> ...
%! loss (C, X, y, 'Weights', ones (3, 1))
%!error<CompactRegressionEnsemble.loss: 'LossFun' must be 'mse' or a function handle.> ...
%! loss (C, X, y, 'LossFun', 'mae')
%!error<CompactRegressionEnsemble.loss: Y must be a real numeric vector.> ...
%! loss (C, X, {1})
%!error<CompactRegressionEnsemble.loss: X and Y must have the same number of rows.> ...
%! loss (C, X, y(1:3))
%!error<CompactRegressionEnsemble.loss: 'Weights' must not be zero for every observation with a response.> ...
%! loss (C, X, y, 'Weights', zeros (150, 1))
%!error<CompactRegressionEnsemble.removeLearners: too few input arguments.> ...
%! removeLearners (C)
%!error<CompactRegressionEnsemble.removeLearners: IDX must be a vector of indices of trained learners.> ...
%! removeLearners (C, 5)
%!error<CompactRegressionEnsemble: 'ResponseTransform' must be a character vector or a function handle.> ...
%! D = C;
%! D.ResponseTransform = 1;

%!test  # MATLAB parity: importance is the weighted average over the trees
%! M = fitrensemble (X, y, 'NumLearningCycles', 3, 'LearnRate', 0.5, ...
%!                   'Learners', templateTree ('MaxNumSplits', 1));
%! [imp, ma] = predictorImportance (compact (M));
%! assert_equal (imp, [0, 0.223154680811778, 0], 1e-13);
%! assert_equal (ma, []);

%!test  # MATLAB parity: a row no tree may predict is left out of the loss
%! U = true (150, 4);
%! U(1,:) = false;
%! e = mean ((predict (C, X(2:end,:)) - y(2:end)) .^ 2);
%! assert_equal (loss (C, X, y, 'UseObsForLearner', U), e, 1e-14);
