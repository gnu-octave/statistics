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

classdef RegressionPartitionedEnsemble
  ## -*- texinfo -*-
  ## @deftp {statistics} RegressionPartitionedEnsemble
  ##
  ## Cross-validated ensemble of regression trees
  ##
  ## A @code{RegressionPartitionedEnsemble} object holds one ensemble per fold
  ## of a partition, each fitted on the observations the fold keeps for
  ## training, and predicts every observation with the ensemble of the fold
  ## that held it out.
  ##
  ## Create one with @code{fitrensemble} given a cross-validation option, or
  ## with the @code{crossval} method of a @code{RegressionEnsemble} or
  ## @code{RegressionBaggedEnsemble}.
  ##
  ## @seealso{fitrensemble, RegressionEnsemble.crossval,
  ## CompactRegressionEnsemble, cvpartition}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} Trainable
    ##
    ## The full ensemble of each fold
    ##
    ## A column cell array with one @code{RegressionEnsemble} or
    ## @code{RegressionBaggedEnsemble} per fold, which @code{resume} grows.
    ## This property is read-only.
    ##
    ## @end deftp
    Trainable = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} NumTrainedPerFold
    ##
    ## Number of trees in each fold
    ##
    ## A row vector.  This property is read-only.
    ##
    ## @end deftp
    NumTrainedPerFold = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} CrossValidatedModel
    ##
    ## Method of the cross-validated ensemble
    ##
    ## @qcode{'LSBoost'} or @qcode{'Bag'}.  This property is read-only.
    ##
    ## @end deftp
    CrossValidatedModel = '';

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} CategoricalPredictors
    ##
    ## Indices of categorical predictors
    ##
    ## Always empty.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName = 'Y';

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations = 0;

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} X
    ##
    ## Predictor data
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} Y
    ##
    ## Response data
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} W
    ##
    ## Observation weights
    ##
    ## The weights of the ensemble that was cross-validated, by which the
    ## losses weigh the held-out observations.  This property is read-only.
    ##
    ## @end deftp
    W = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} ModelParameters
    ##
    ## Parameters of the cross-validation
    ##
    ## A structure with the fields @qcode{Type}, @qcode{Method},
    ## @qcode{'PartitionedEnsemble'}, @qcode{LearnerTemplates}, the
    ## @code{ModelParameters} of the ensemble the folds were fitted as, and
    ## @qcode{NLearn}, the number of folds.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} Trained
    ##
    ## The compact ensemble of each fold
    ##
    ## A column cell array with one @code{CompactRegressionEnsemble} per fold.
    ## This property is read-only.
    ##
    ## @end deftp
    Trained = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} KFold
    ##
    ## Number of folds
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    KFold = 0;

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} Partition
    ##
    ## The partition of the observations
    ##
    ## A @code{cvpartition} object.  This property is read-only.
    ##
    ## @end deftp
    Partition = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## Always empty.  This property is read-only.
    ##
    ## @end deftp
    BinEdges = {};

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedEnsemble} {property} ResponseTransform
    ##
    ## Transform applied to the out-of-fold predictions
    ##
    ## The folds carry none; this one is applied once to what they predict.
    ## MATLAB R2024a leaves the transform on the folds as well and so applies
    ## it twice.
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
        parseResponseTransform (val, 'RegressionPartitionedEnsemble');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionPartitionedEnsemble} {@var{CVMdl} =} RegressionPartitionedEnsemble (@var{Mdl}, @var{Partition})
    ##
    ## Cross-validate a regression ensemble.
    ##
    ## @var{Mdl} is the ensemble to cross-validate and @var{Partition} a
    ## @code{cvpartition} of its observations.  The documented ways in are
    ## @code{fitrensemble} and @code{crossval}.
    ##
    ## @end deftypefn
    function this = RegressionPartitionedEnsemble (Mdl, Partition)

      if (nargin < 2)
        error (strcat ("RegressionPartitionedEnsemble: too few input", ...
                       " arguments."));
      endif
      if (! isa (Mdl, 'RegressionEnsemble'))
        error (strcat ("RegressionPartitionedEnsemble: MDL must be a", ...
                       " 'RegressionEnsemble' object."));
      endif
      if (! isa (Partition, 'cvpartition'))
        error (strcat ("RegressionPartitionedEnsemble: PARTITION must be a", ...
                       " 'cvpartition' object."));
      endif
      if (numel (test (Partition, 1)) != Mdl.NumObservations)
        error (strcat ("RegressionPartitionedEnsemble: PARTITION must", ...
                       " partition the observations MDL was fitted on."));
      endif

      this.CrossValidatedModel   = Mdl.Method;
      this.PredictorNames        = Mdl.PredictorNames;
      this.CategoricalPredictors = Mdl.CategoricalPredictors;
      this.ResponseName          = Mdl.ResponseName;
      this.NumObservations       = Mdl.NumObservations;
      this.X                     = Mdl.X;
      this.Y                     = Mdl.Y;
      this.W                     = Mdl.W;
      this.Partition             = Partition;
      this.KFold                 = Partition.NumTestSets;
      this.ModelParameters = struct ('Type', 'regression', ...
                                     'Method', 'PartitionedEnsemble', ...
                                     'LearnerTemplates', ...
                                     Mdl.ModelParameters, ...
                                     'NLearn', this.KFold);
      if (any (strcmp (Mdl.ResponseTransform, {'none', 'exp', 'log'})))
        this.ResponseTransform = Mdl.ResponseTransform;
      else
        this.ResponseTransform = Mdl.RTfun;
      endif

      mp = Mdl.ModelParameters;
      fargs = {'Method', Mdl.Method, 'NumLearningCycles', mp.NLearn, ...
               'Learners', mp.LearnerTemplates, ...
               'PredictorNames', Mdl.PredictorNames, ...
               'ResponseName', Mdl.ResponseName};
      if (isa (Mdl, 'RegressionBaggedEnsemble'))
        onoff = {'off', 'on'};
        fargs(end+1:end+4) = {'FResample', Mdl.FResample, ...
                              'Replace', onoff{Mdl.Replace + 1}};
      endif
      if (strcmp (Mdl.Method, 'LSBoost'))
        fargs(end+1:end+2) = {'LearnRate', mp.LearnRate};
      endif
      ## Assigned field by field: struct () given a cell array of labels
      ## would build one structure per label.
      F = struct ();
      F.X = Mdl.X;
      F.Y = Mdl.Y;
      F.Weights = Mdl.W;
      this.Trainable = foldModels (class (Mdl), F, Partition, fargs);
      this = compactFolds (this);

    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ("\n  RegressionPartitionedEnsemble\n\n");
      fprintf ("%+25s: '%s'\n", 'CrossValidatedModel', ...
               this.CrossValidatedModel);
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'KFold', this.KFold);
      fprintf ("%+25s: %s\n", 'NumTrainedPerFold', ...
               mat2str (this.NumTrainedPerFold));
      fprintf ("%+25s: '%s'\n", 'ResponseTransform', this.ResponseTransform);
      fprintf ("\n");
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {RegressionPartitionedEnsemble} {@var{yfit} =} kfoldPredict (@var{obj})
    ##
    ## Predict each observation with the fold that held it out.
    ##
    ## @var{yfit} holds each observation's prediction from the ensemble of the
    ## fold that held it out, after @code{ResponseTransform}, and @code{NaN}
    ## for an observation no fold held out.
    ##
    ## @seealso{RegressionPartitionedEnsemble,
    ## RegressionPartitionedEnsemble.kfoldLoss}
    ## @end deftypefn
    function yfit = kfoldPredict (this)

      yfit = this.RTfun (foldResponse (this, 1:this.KFold));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionPartitionedEnsemble} {@var{L} =} kfoldLoss (@var{obj})
    ## @deftypefnx {RegressionPartitionedEnsemble} {@var{L} =} kfoldLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validated regression loss.
    ##
    ## @var{L} is the weighted mean squared error of the out-of-fold
    ## predictions, the held-out observations weighted by @code{W}.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.2 0.02 0.78
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Folds'} @tab @tab The folds to use, pooled.  The default
    ## is all of them.
    ## @item @qcode{'LossFun'} @tab @tab @qcode{'mse'} (default) or a function
    ## handle called as @code{lossfun (Y, Yfit, W)} with normalized weights.
    ## @item @qcode{'Mode'} @tab @tab @qcode{'average'} (default) for one loss
    ## over the observations of every fold used, @qcode{'individual'} for a
    ## column with the loss of each fold, or @qcode{'cumulative'} for a column
    ## whose element @math{t} uses the first @math{t} trees of every fold.
    ## @end multitable
    ##
    ## @seealso{RegressionPartitionedEnsemble,
    ## RegressionPartitionedEnsemble.kfoldPredict}
    ## @end deftypefn
    function L = kfoldLoss (this, varargin)

      o = kfoldOpts (varargin, {'mse'}, 'RegressionPartitionedEnsemble', ...
                     'kfoldLoss', this.KFold, ...
                     {'average', 'individual', 'cumulative'}, true);
      LossFun = o.LossFun;
      if (isempty (LossFun))
        LossFun = 'mse';
      endif
      Folds = o.Folds;
      Mode = o.Mode;
      n = this.NumObservations;

      if (! strcmp (Mode, 'cumulative'))
        sets = foldSets (this.Partition, Folds, Mode, n);
        yf = foldResponse (this, Folds);
        L = zeros (numel (sets), 1);
        for i = 1:numel (sets)
          L(i) = setLoss (this, yf, sets{i}, LossFun);
        endfor
        return;
      endif
      ## The cumulative mode pools the folds used and grows every fold's
      ## trees together.
      pooled = foldSets (this.Partition, Folds, 'average', n){1};
      T = max ([0, this.NumTrainedPerFold(Folds)]);
      Yf = NaN (n, T);
      for k = Folds
        te = test (this.Partition, k);
        Tk = this.NumTrainedPerFold(k);
        if (any (te) && Tk > 0)
          Yk = ensembleSteps (this.Trained{k}, this.X(te,:));
          Yf(te,:) = [Yk, repmat(Yk(:,end), 1, T - Tk)];
        endif
      endfor
      L = zeros (T, 1);
      for t = 1:T
        L(t) = setLoss (this, Yf(:,t), pooled, LossFun);
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionPartitionedEnsemble} {@var{vals} =} kfoldfun (@var{obj}, @var{fun})
    ##
    ## Apply a function to each fold.
    ##
    ## @var{fun} is called once per fold as
    ## @code{@var{fun} (CMP, Xtrain, Ytrain, Wtrain, Xtest, Ytest, Wtest)},
    ## @var{CMP} being the fold's compact ensemble and the weights those of
    ## @code{W}, and must return a row.  @var{vals} stacks the rows.
    ##
    ## @seealso{RegressionPartitionedEnsemble}
    ## @end deftypefn
    function vals = kfoldfun (this, fun)

      if (nargin < 2)
        error (strcat ("RegressionPartitionedEnsemble.kfoldfun: too few", ...
                       " input arguments."));
      endif
      if (! is_function_handle (fun))
        error (strcat ("RegressionPartitionedEnsemble.kfoldfun: FUN must", ...
                       " be a function handle."));
      endif
      vals = [];
      for k = 1:this.KFold
        tr = training (this.Partition, k);
        te = test (this.Partition, k);
        v = fun (this.Trained{k}, this.X(tr,:), this.Y(tr), this.W(tr), ...
                 this.X(te,:), this.Y(te), this.W(te));
        vals = [vals; v];
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionPartitionedEnsemble} {@var{CVMdl} =} resume (@var{obj}, @var{NumLearningCycles})
    ## @deftypefnx {RegressionPartitionedEnsemble} {@var{CVMdl} =} resume (@dots{}, 'NPrint', @var{n})
    ##
    ## Grow more trees in every fold.
    ##
    ## Each fold's ensemble is resumed, as @code{RegressionEnsemble.resume}
    ## does, by @var{NumLearningCycles} trees.
    ##
    ## @seealso{RegressionPartitionedEnsemble, RegressionEnsemble.resume}
    ## @end deftypefn
    function this = resume (this, varargin)

      if (nargin < 2)
        error (strcat ("RegressionPartitionedEnsemble.resume: too few", ...
                       " input arguments."));
      endif
      for k = 1:this.KFold
        this.Trainable{k} = resume (this.Trainable{k}, varargin{:});
      endfor
      this = compactFolds (this);

    endfunction

  endmethods

  methods (Access = private)

    function this = compactFolds (this)
      this.Trained = cellfun (@(M) compact (M), this.Trainable, ...
                              'UniformOutput', false);
      this.NumTrainedPerFold = cellfun (@(M) M.NumTrained, this.Trainable)';
    endfunction

    ## The untransformed predictions for the rows the FOLDS held out, NaN
    ## elsewhere.
    function yf = foldResponse (this, folds)
      yf = NaN (this.NumObservations, 1);
      for k = folds
        te = test (this.Partition, k);
        if (any (te) && this.Trained{k}.NumTrained > 0)
          C = this.Trained{k};
          C.ResponseTransform = 'none';
          yf(te) = predict (C, this.X(te,:));
        endif
      endfor
    endfunction

    ## The weighted loss of the untransformed predictions YF over ROWS.
    function L = setLoss (this, yf, rows, LossFun)
      f = this.RTfun (yf(rows));
      y = this.Y(rows);
      w = this.W(rows);
      have = ! isnan (f);
      w = w(have);
      if (! (sum (w) > 0))
        L = NaN;
        return;
      endif
      w /= sum (w);
      if (is_function_handle (LossFun))
        L = LossFun (y(have), f(have), w);
      else
        L = sum (w .* (f(have) - y(have)) .^ 2);
      endif
    endfunction

  endmethods

endclassdef

## Test output
%!shared X, y, S, c, CV
%! load fisheriris
%! X = meas(:,2:4);
%! y = meas(:,1);
%! S = templateTree ('MaxNumSplits', 1);
%! c = cvpartition (150, 'KFold', 5);
%! CV = fitrensemble (X, y, 'NumLearningCycles', 3, 'Learners', S, ...
%!                    'CVPartition', c);

%!test  # MATLAB parity: the properties of a cross-validated ensemble
%! assert_equal (class (CV), 'RegressionPartitionedEnsemble');
%! assert_equal (numel (properties (CV)), 16);
%! assert_equal (CV.KFold, 5);
%! assert_equal (CV.NumTrainedPerFold, [3, 3, 3, 3, 3]);
%! assert_equal (CV.CrossValidatedModel, 'LSBoost');
%! assert_equal (class (CV.Trained{1}), 'CompactRegressionEnsemble');
%! assert_equal (class (CV.Trainable{1}), 'RegressionEnsemble');

%!test  # MATLAB parity: each row is predicted by the fold that held it out
%! yf = zeros (150, 1);
%! for k = 1:5
%!   te = test (c, k);
%!   yf(te) = predict (CV.Trained{k}, X(te,:));
%! endfor
%! assert_equal (kfoldPredict (CV), yf);
%! assert_equal (kfoldLoss (CV), mean ((yf - y) .^ 2), 1e-14);
%! Li = kfoldLoss (CV, 'Mode', 'individual');
%! assert_equal (Li(4), mean ((yf(test (c, 4)) - y(test (c, 4))) .^ 2), 1e-14);

%!test  # MATLAB parity: the cumulative loss grows the trees of every fold
%! Lc = kfoldLoss (CV, 'Mode', 'cumulative');
%! assert_equal (size (Lc), [3, 1]);
%! yf = zeros (150, 1);
%! for k = 1:5
%!   te = test (c, k);
%!   yf(te) = predict (CV.Trained{k}, X(te,:), 'Learners', 1:2);
%! endfor
%! assert_equal (Lc(2), mean ((yf - y) .^ 2), 1e-14);
%! assert_equal (Lc(3), kfoldLoss (CV), 1e-14);

%!test  # the response transform is applied once, where MATLAB applies it twice
%! M = fitrensemble (X, y, 'NumLearningCycles', 3, 'Learners', S, ...
%!                   'CVPartition', c, 'ResponseTransform', @(z) 2 * z);
%! assert_equal (M.Trained{1}.ResponseTransform, 'none');
%! assert_equal (kfoldPredict (M), 2 * kfoldPredict (CV), 1e-14);

%!test  # MATLAB parity: kfoldfun and resume
%! f = kfoldfun (CV, @(C, Xtr, Ytr, Wtr, Xte, Yte, Wte) [rows(Xtr), rows(Xte)]);
%! assert_equal (f(1,:), [120, 30]);
%! R = resume (CV, 1);
%! assert_equal (R.NumTrainedPerFold, [4, 4, 4, 4, 4]);

%!test  # MATLAB parity: the losses weigh the held-out rows by W
%! w = [5 * ones(50, 1); ones(100, 1)];
%! M = fitrensemble (X, y, 'NumLearningCycles', 3, 'Learners', S, ...
%!                   'CVPartition', c, 'Weights', w);
%! yf = kfoldPredict (M);
%! assert_equal (kfoldLoss (M), sum (M.W .* (yf - y) .^ 2) / sum (M.W), 1e-14);

## Test input validation
%!error<RegressionPartitionedEnsemble: too few input arguments.> ...
%! RegressionPartitionedEnsemble (1)
%!error<RegressionPartitionedEnsemble: MDL must be a 'RegressionEnsemble' object.> ...
%! RegressionPartitionedEnsemble (1, c)
%!error<RegressionPartitionedEnsemble: PARTITION must be a 'cvpartition' object.> ...
%! RegressionPartitionedEnsemble (CV.Trainable{1}, 1)
%!error<RegressionPartitionedEnsemble.kfoldLoss: optional arguments must be given in Name-Value pairs.> ...
%! kfoldLoss (CV, 'Mode')
%!error<RegressionPartitionedEnsemble.kfoldLoss: invalid parameter name in optional pair arguments.> ...
%! kfoldLoss (CV, 'Learners', 1)
%!error<RegressionPartitionedEnsemble.kfoldLoss: 'Mode' must be 'average', 'individual', or 'cumulative'.> ...
%! kfoldLoss (CV, 'Mode', 'ensemble')
%!error<RegressionPartitionedEnsemble.kfoldLoss: 'Folds' must hold integers between 1 and 5.> ...
%! kfoldLoss (CV, 'Folds', 0)
%!error<RegressionPartitionedEnsemble.kfoldLoss: 'LossFun' must be 'mse' or a function handle.> ...
%! kfoldLoss (CV, 'LossFun', 'mae')
%!error<RegressionPartitionedEnsemble.kfoldfun: FUN must be a function handle.> ...
%! kfoldfun (CV, 1)
%!error<RegressionPartitionedEnsemble.resume: too few input arguments.> ...
%! resume (CV)
