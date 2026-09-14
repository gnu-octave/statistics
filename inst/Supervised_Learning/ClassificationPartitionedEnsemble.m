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

classdef ClassificationPartitionedEnsemble
  ## -*- texinfo -*-
  ## @deftp {statistics} ClassificationPartitionedEnsemble
  ##
  ## Cross-validated ensemble of decision trees for classification
  ##
  ## A @code{ClassificationPartitionedEnsemble} object holds one ensemble per
  ## fold of a partition, each fitted on the observations the fold keeps for
  ## training, and answers for every observation with the ensemble of the fold
  ## that held it out.
  ##
  ## Create one with @code{fitcensemble} given a cross-validation option, or
  ## with the @code{crossval} method of a @code{ClassificationEnsemble} or
  ## @code{ClassificationBaggedEnsemble}.
  ##
  ## @seealso{fitcensemble, ClassificationEnsemble.crossval,
  ## CompactClassificationEnsemble, cvpartition}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} Trainable
    ##
    ## The full ensemble of each fold
    ##
    ## A column cell array with one @code{ClassificationEnsemble} or
    ## @code{ClassificationBaggedEnsemble} per fold, which @code{resume}
    ## grows.  This property is read-only.
    ##
    ## @end deftp
    Trainable = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} NumTrainedPerFold
    ##
    ## Number of learners in each fold
    ##
    ## A row vector.  This property is read-only.
    ##
    ## @end deftp
    NumTrainedPerFold = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} CrossValidatedModel
    ##
    ## Method of the cross-validated ensemble
    ##
    ## The @code{Method} of the ensemble, such as @qcode{'AdaBoostM1'} or
    ## @qcode{'Bag'}.  This property is read-only.
    ##
    ## @end deftp
    CrossValidatedModel = '';

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} CategoricalPredictors
    ##
    ## Indices of categorical predictors
    ##
    ## Always empty.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName = 'Y';

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations = 0;

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} X
    ##
    ## Predictor data
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} Y
    ##
    ## Class labels
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} W
    ##
    ## Observation weights
    ##
    ## The weights of the ensemble that was cross-validated.  The losses and
    ## edges weigh the held-out observations by them.  This property is
    ## read-only.
    ##
    ## @end deftp
    W = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} ModelParameters
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
    ## @deftp {ClassificationPartitionedEnsemble} {property} Trained
    ##
    ## The compact ensemble of each fold
    ##
    ## A column cell array with one @code{CompactClassificationEnsemble} per
    ## fold.  This property is read-only.
    ##
    ## @end deftp
    Trained = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} KFold
    ##
    ## Number of folds
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    KFold = 0;

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} Partition
    ##
    ## The partition of the observations
    ##
    ## A @code{cvpartition} object.  This property is read-only.
    ##
    ## @end deftp
    Partition = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## Always empty.  This property is read-only.
    ##
    ## @end deftp
    BinEdges = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} ClassNames
    ##
    ## Names of the classes
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} Cost
    ##
    ## Misclassification costs
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Cost = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} Prior
    ##
    ## Prior probabilities of the classes
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Prior = [];

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedEnsemble} {property} ScoreTransform
    ##
    ## Transform applied to the out-of-fold scores
    ##
    ## The folds carry none; this one is applied once to the scores they
    ## return.
    ##
    ## @end deftp
    ScoreTransform = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    STfun = @(x) x;      # the score transform as a function
    gY = [];             # class index of each observation
    DefaultIndex = 1;    # class given to an observation no fold held out
  endproperties

  methods (Hidden)

    function this = set.ScoreTransform (this, val)
      try
        [this.STfun, this.ScoreTransform] = parseScoreTransform (val, ...
                                      'ClassificationPartitionedEnsemble');
      catch
        error (strcat ("ClassificationPartitionedEnsemble.subsasgn:", ...
                       " 'ScoreTransform' must be a character vector or a", ...
                       " 'function_handle' object."));
      end_try_catch
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationPartitionedEnsemble} {@var{CVMdl} =} ClassificationPartitionedEnsemble (@var{Mdl}, @var{Partition})
    ##
    ## Cross-validate a classification ensemble.
    ##
    ## @var{Mdl} is the ensemble to cross-validate and @var{Partition} a
    ## @code{cvpartition} of its observations.  The documented ways in are
    ## @code{fitcensemble} and @code{crossval}.
    ##
    ## @end deftypefn
    function this = ClassificationPartitionedEnsemble (Mdl, Partition)

      if (nargin < 2)
        error (strcat ("ClassificationPartitionedEnsemble: too few input", ...
                       " arguments."));
      endif
      if (! isa (Mdl, 'ClassificationEnsemble'))
        error (strcat ("ClassificationPartitionedEnsemble: MDL must be a", ...
                       " 'ClassificationEnsemble' object."));
      endif
      if (! isa (Partition, 'cvpartition'))
        error (strcat ("ClassificationPartitionedEnsemble: PARTITION must", ...
                       " be a 'cvpartition' object."));
      endif
      if (numel (test (Partition, 1)) != Mdl.NumObservations)
        error (strcat ("ClassificationPartitionedEnsemble: PARTITION must", ...
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
      this.ClassNames            = Mdl.ClassNames;
      this.Cost                  = Mdl.Cost;
      this.Prior                 = Mdl.Prior;
      this.Partition             = Partition;
      this.KFold                 = Partition.NumTestSets;
      this.gY                    = Mdl.gY;
      this.DefaultIndex          = Mdl.DefaultIndex;
      this.ModelParameters = struct ('Type', 'classification', ...
                                     'Method', 'PartitionedEnsemble', ...
                                     'LearnerTemplates', ...
                                     Mdl.ModelParameters, ...
                                     'NLearn', this.KFold);
      if (any (strcmp (Mdl.ScoreTransform, {'doublelogit', 'invlogit', ...
                       'ismax', 'logit', 'none', 'identity', 'sign', ...
                       'symmetric', 'symmetricismax', 'symmetriclogit'})))
        this.ScoreTransform = Mdl.ScoreTransform;
      else
        this.ScoreTransform = Mdl.STfun;
      endif

      ## Each fold is fitted with the options the parent settled, its classes,
      ## prior and cost included, and carries no score transform of its own:
      ## this object applies one, once, to what the folds return.
      mp = Mdl.ModelParameters;
      fargs = {'Method', Mdl.Method, 'NumLearningCycles', mp.NLearn, ...
               'Learners', mp.LearnerTemplates, ...
               'ClassNames', Mdl.ClassNames, 'Prior', Mdl.Prior, ...
               'Cost', Mdl.Cost, 'PredictorNames', Mdl.PredictorNames, ...
               'ResponseName', Mdl.ResponseName};
      if (isa (Mdl, 'ClassificationBaggedEnsemble'))
        onoff = {'off', 'on'};
        fargs(end+1:end+4) = {'FResample', Mdl.FResample, ...
                              'Replace', onoff{Mdl.Replace + 1}};
      endif
      if (strcmp (Mdl.Method, 'Subspace'))
        fargs(end+1:end+2) = {'NPredToSample', Mdl.NPredToSample};
        if (Mdl.AllCombinations)
          fargs{4} = 'AllPredictorCombinations';
        endif
      elseif (! strcmp (Mdl.Method, 'Bag'))
        fargs(end+1:end+2) = {'LearnRate', mp.LearnRate};
        if (strcmp (Mdl.Method, 'RUSBoost'))
          fargs(end+1:end+2) = {'RatioToSmallest', Mdl.RatioToSmallest};
        endif
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
      fprintf ("\n  ClassificationPartitionedEnsemble\n\n");
      fprintf ("%+25s: '%s'\n", 'CrossValidatedModel', ...
               this.CrossValidatedModel);
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'KFold', this.KFold);
      fprintf ("%+25s: %s\n", 'NumTrainedPerFold', ...
               mat2str (this.NumTrainedPerFold));
      fprintf ("%+25s: %s\n", 'ClassNames', classNameListing (this.ClassNames));
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("\n");
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedEnsemble} {@var{label} =} kfoldPredict (@var{obj})
    ## @deftypefnx {ClassificationPartitionedEnsemble} {[@var{label}, @var{scores}] =} kfoldPredict (@var{obj})
    ##
    ## Classify each observation with the fold that held it out.
    ##
    ## @var{scores} holds each observation's scores from the ensemble of the
    ## fold that held it out, after @code{ScoreTransform}, and @var{label}
    ## the class of highest score.  An observation no fold held out, as under
    ## a holdout partition, has @code{NaN} scores and the class of greatest
    ## prior probability.
    ##
    ## @seealso{ClassificationPartitionedEnsemble,
    ## ClassificationPartitionedEnsemble.kfoldLoss}
    ## @end deftypefn
    function [label, scores] = kfoldPredict (this)

      scores = this.STfun (foldScores (this, 1:this.KFold));
      [~, k] = max (scores, [], 2);
      k(any (isnan (scores), 2)) = this.DefaultIndex;
      label = labelsFromIndex (this.ClassNames, k);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedEnsemble} {@var{L} =} kfoldLoss (@var{obj})
    ## @deftypefnx {ClassificationPartitionedEnsemble} {@var{L} =} kfoldLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validated classification loss.
    ##
    ## @var{L} is the loss of the out-of-fold scores, the held-out observations
    ## weighted by @code{W}, those without scores left out.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.2 0.02 0.78
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Folds'} @tab @tab The folds to use, pooled.  The default
    ## is all of them.
    ## @item @qcode{'LossFun'} @tab @tab A loss
    ## @code{CompactClassificationEnsemble.loss} accepts.  The default is
    ## @qcode{'classiferror'}.
    ## @item @qcode{'Mode'} @tab @tab @qcode{'average'} (default) for one loss
    ## over the observations of every fold used, @qcode{'individual'} for a
    ## column with the loss of each fold, or @qcode{'cumulative'} for a column
    ## whose element @math{t} uses the first @math{t} learners of every fold.
    ## @end multitable
    ##
    ## @seealso{ClassificationPartitionedEnsemble,
    ## ClassificationPartitionedEnsemble.kfoldEdge}
    ## @end deftypefn
    function L = kfoldLoss (this, varargin)

      o = kfoldOpts (varargin, CLASSIFICATION_LOSSES, ...
                     'ClassificationPartitionedEnsemble', 'kfoldLoss', ...
                     this.KFold, ENSEMBLE_MODES, true);
      if (isempty (o.LossFun))
        o.LossFun = 'classiferror';
      endif
      [sets, S] = foldSets_ (this, o);
      L = zeros (numel (sets), 1);
      for i = 1:numel (sets)
        L(i) = setLoss (this, S(:,:,min (i, size (S, 3))), sets{i}, o.LossFun);
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedEnsemble} {@var{e} =} kfoldEdge (@var{obj})
    ## @deftypefnx {ClassificationPartitionedEnsemble} {@var{e} =} kfoldEdge (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validated classification edge.
    ##
    ## The weighted mean of the out-of-fold margins, weighted by @code{W}.
    ## @qcode{'Folds'} and @qcode{'Mode'} are taken as by @code{kfoldLoss}.
    ##
    ## @seealso{ClassificationPartitionedEnsemble,
    ## ClassificationPartitionedEnsemble.kfoldMargin}
    ## @end deftypefn
    function e = kfoldEdge (this, varargin)

      o = kfoldOpts (varargin, {}, 'ClassificationPartitionedEnsemble', ...
                     'kfoldEdge', this.KFold, ENSEMBLE_MODES);
      [sets, S] = foldSets_ (this, o);
      e = zeros (numel (sets), 1);
      for i = 1:numel (sets)
        rows = sets{i};
        m = marginsOf (this.STfun (S(rows,:,min (i, size (S, 3)))), ...
                       this.gY(rows), 1);
        w = this.W(rows);
        have = ! isnan (m);
        e(i) = sum (w(have) .* m(have)) / sum (w(have));
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationPartitionedEnsemble} {@var{m} =} kfoldMargin (@var{obj})
    ##
    ## Cross-validated classification margins.
    ##
    ## @var{m} holds each observation's margin under its out-of-fold scores,
    ## @code{NaN} for one no fold held out.
    ##
    ## @seealso{ClassificationPartitionedEnsemble,
    ## ClassificationPartitionedEnsemble.kfoldEdge}
    ## @end deftypefn
    function m = kfoldMargin (this)

      [~, scores] = kfoldPredict (this);
      m = marginsOf (scores, this.gY, 1);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationPartitionedEnsemble} {@var{vals} =} kfoldfun (@var{obj}, @var{fun})
    ##
    ## Apply a function to each fold.
    ##
    ## @var{fun} is called once per fold as
    ## @code{@var{fun} (CMP, Xtrain, Ytrain, Wtrain, Xtest, Ytest, Wtest)},
    ## @var{CMP} being the fold's compact ensemble and the weights those of
    ## @code{W}, and must return a row.  @var{vals} stacks the rows.
    ##
    ## @seealso{ClassificationPartitionedEnsemble}
    ## @end deftypefn
    function vals = kfoldfun (this, fun)

      if (nargin < 2)
        error (strcat ("ClassificationPartitionedEnsemble.kfoldfun: too", ...
                       " few input arguments."));
      endif
      if (! is_function_handle (fun))
        error (strcat ("ClassificationPartitionedEnsemble.kfoldfun: FUN", ...
                       " must be a function handle."));
      endif
      vals = [];
      for k = 1:this.KFold
        tr = training (this.Partition, k);
        te = test (this.Partition, k);
        v = fun (this.Trained{k}, this.X(tr,:), this.Y(tr,:), this.W(tr), ...
                 this.X(te,:), this.Y(te,:), this.W(te));
        vals = [vals; v];
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedEnsemble} {@var{CVMdl} =} resume (@var{obj}, @var{NumLearningCycles})
    ## @deftypefnx {ClassificationPartitionedEnsemble} {@var{CVMdl} =} resume (@dots{}, 'NPrint', @var{n})
    ##
    ## Grow more learners in every fold.
    ##
    ## Each fold's ensemble is resumed, as @code{ClassificationEnsemble.resume}
    ## does, by @var{NumLearningCycles} learners.
    ##
    ## @seealso{ClassificationPartitionedEnsemble, ClassificationEnsemble.resume}
    ## @end deftypefn
    function this = resume (this, varargin)

      if (nargin < 2)
        error (strcat ("ClassificationPartitionedEnsemble.resume: too few", ...
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

    ## The untransformed scores of the rows the FOLDS held out, from the
    ## ensemble of the fold that held each out, NaN elsewhere.
    function S = foldScores (this, folds)
      K = classCount (this.ClassNames);
      S = NaN (this.NumObservations, K);
      for k = folds
        te = test (this.Partition, k);
        if (any (te) && this.Trained{k}.NumTrained > 0)
          [~, S(te,:)] = predict (this.Trained{k}, this.X(te,:));
        endif
      endfor
    endfunction

    ## The row sets a kfold method reports over and the untransformed scores
    ## to read them with: one pooled set for 'average', one per fold for
    ## 'individual', and for 'cumulative' one pooled set per learner count,
    ## the scores then NxKxT.
    function [sets, S] = foldSets_ (this, o)
      n = this.NumObservations;
      if (! strcmp (o.Mode, 'cumulative'))
        sets = foldSets (this.Partition, o.Folds, o.Mode, n);
        S = foldScores (this, o.Folds);
        return;
      endif
      pooled = foldSets (this.Partition, o.Folds, 'average', n);
      K = classCount (this.ClassNames);
      T = max ([0, this.NumTrainedPerFold(o.Folds)]);
      S = NaN (n, K, T);
      for k = o.Folds
        te = test (this.Partition, k);
        Tk = this.NumTrainedPerFold(k);
        if (any (te) && Tk > 0)
          Sk = ensembleSteps (this.Trained{k}, this.X(te,:));
          S(te,:,:) = cat (3, Sk, repmat (Sk(:,:,end), 1, 1, T - Tk));
        endif
      endfor
      sets = repmat (pooled, T, 1);
    endfunction

    ## The weighted loss of the untransformed scores S over ROWS.
    function L = setLoss (this, S, rows, LossFun)
      S = this.STfun (S(rows,:));
      g = this.gY(rows);
      w = this.W(rows);
      have = ! any (isnan (S), 2);
      w = w(have);
      if (! (sum (w) > 0))
        L = NaN;
        return;
      endif
      w /= sum (w);
      S = S(have,:);
      g = g(have);
      if (is_function_handle (LossFun))
        C = false (size (S));
        C(sub2ind (size (S), (1:rows (S))', g)) = true;
        L = LossFun (C, S, w, this.Cost);
      else
        L = classificationLoss (LossFun, S, g, w, this.Cost);
      endif
    endfunction

  endmethods

endclassdef

## The modes of an ensemble's kfold methods, the first the default.
function m = ENSEMBLE_MODES ()
  m = {'average', 'individual', 'cumulative'};
endfunction

## The classification losses a kfold method takes by name.
function l = CLASSIFICATION_LOSSES ()
  l = {'binodeviance', 'classifcost', 'classiferror', 'exponential', ...
       'hinge', 'logit', 'mincost', 'quadratic'};
endfunction

## Test output
%!shared X2, Y2, S, c, CV
%! load fisheriris
%! X2 = meas(51:150,:);
%! Y2 = species(51:150);
%! S = templateTree ('MaxNumSplits', 1);
%! c = cvpartition (Y2, 'KFold', 5);
%! CV = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                    'NumLearningCycles', 4, 'Learners', S, ...
%!                    'CVPartition', c);

%!test  # MATLAB parity: the properties of a cross-validated ensemble
%! assert_equal (class (CV), 'ClassificationPartitionedEnsemble');
%! assert_equal (numel (properties (CV)), 19);
%! assert_equal (CV.KFold, 5);
%! assert_equal (CV.NumTrainedPerFold, [4, 4, 4, 4, 4]);
%! assert_equal (CV.CrossValidatedModel, 'AdaBoostM1');
%! assert_equal (size (CV.Trained), [5, 1]);
%! assert_equal (class (CV.Trained{1}), 'CompactClassificationEnsemble');
%! assert_equal (class (CV.Trainable{1}), 'ClassificationEnsemble');
%! assert_equal (CV.ModelParameters.Method, 'PartitionedEnsemble');
%! assert_equal (CV.ModelParameters.NLearn, 5);

%!test  # MATLAB parity: each row is predicted by the fold that held it out
%! S0 = zeros (100, 2);
%! for k = 1:5
%!   te = test (c, k);
%!   [~, S0(te,:)] = predict (CV.Trained{k}, X2(te,:));
%! endfor
%! [label, s] = kfoldPredict (CV);
%! assert_equal (s, S0);
%! [~, j] = max (S0, [], 2);
%! assert_equal (label, CV.ClassNames(j));

%!test  # MATLAB parity: the loss pooled over the folds and fold by fold
%! [~, s] = kfoldPredict (CV);
%! g = 1 + strcmp (Y2, 'virginica');
%! miss = s(sub2ind (size (s), (1:100)', g)) <= s(sub2ind (size (s), ...
%!                                                        (1:100)', 3 - g));
%! assert_equal (kfoldLoss (CV), mean (miss), 1e-15);
%! Li = kfoldLoss (CV, 'Mode', 'individual');
%! assert_equal (size (Li), [5, 1]);
%! assert_equal (Li(2), mean (miss(test (c, 2))), 1e-15);
%! assert_equal (kfoldLoss (CV, 'Folds', [1, 3]), ...
%!               mean (miss(test (c, 1) | test (c, 3))), 1e-15);

%!test  # MATLAB parity: the cumulative loss grows the learners of every fold
%! g = 1 + strcmp (Y2, 'virginica');
%! Lc = kfoldLoss (CV, 'Mode', 'cumulative');
%! assert_equal (size (Lc), [4, 1]);
%! s = zeros (100, 2);
%! for k = 1:5
%!   te = test (c, k);
%!   [~, s(te,:)] = predict (CV.Trained{k}, X2(te,:), 'Learners', 1:2);
%! endfor
%! miss = s(sub2ind (size (s), (1:100)', g)) <= s(sub2ind (size (s), ...
%!                                                        (1:100)', 3 - g));
%! assert_equal (Lc(2), mean (miss), 1e-15);
%! assert_equal (Lc(4), kfoldLoss (CV), 1e-15);

%!test  # MATLAB parity: the losses weigh the held-out rows by W
%! w = [5 * ones(50, 1); ones(50, 1)];
%! M = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', 'NumLearningCycles', 4, ...
%!                   'Learners', S, 'CVPartition', c, 'Weights', w);
%! [~, s] = kfoldPredict (M);
%! g = 1 + strcmp (Y2, 'virginica');
%! miss = s(sub2ind (size (s), (1:100)', g)) <= s(sub2ind (size (s), ...
%!                                                        (1:100)', 3 - g));
%! assert_equal (kfoldLoss (M), sum (M.W .* miss) / sum (M.W), 1e-15);

%!test  # MATLAB parity: the edge and margins of the out-of-fold scores
%! [~, s] = kfoldPredict (CV);
%! g = 1 + strcmp (Y2, 'virginica');
%! m = s(sub2ind (size (s), (1:100)', g)) - s(sub2ind (size (s), ...
%!                                                     (1:100)', 3 - g));
%! assert_equal (kfoldMargin (CV), m, 1e-14);
%! assert_equal (kfoldEdge (CV), mean (m), 1e-13);
%! assert_equal (size (kfoldEdge (CV, 'Mode', 'cumulative')), [4, 1]);
%! assert_equal (size (kfoldEdge (CV, 'Mode', 'individual')), [5, 1]);

%!test  # MATLAB parity: the score transform is applied once, by the parent
%! M = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', 'NumLearningCycles', 4, ...
%!                   'Learners', S, 'CVPartition', c, ...
%!                   'ScoreTransform', 'doublelogit');
%! assert_equal (M.Trained{1}.ScoreTransform, 'none');
%! [~, s0] = kfoldPredict (CV);
%! [~, s] = kfoldPredict (M);
%! assert_equal (s, 1 ./ (1 + exp (-2 * s0)), 1e-15);

%!test  # MATLAB parity: kfoldfun gets each fold's data and weights
%! f = kfoldfun (CV, @(C, Xtr, Ytr, Wtr, Xte, Yte, Wte) ...
%!               [rows(Xtr), sum(Wtr), rows(Xte), sum(Wte), ...
%!                isa(C, 'CompactClassificationEnsemble')]);
%! assert_equal (f(1,:), [80, 0.8, 20, 0.2, 1], 1e-15);
%! assert_equal (rows (f), 5);

%!test  # MATLAB parity: resuming grows every fold
%! R = resume (CV, 2);
%! assert_equal (R.NumTrainedPerFold, [6, 6, 6, 6, 6]);
%! assert_equal (class (R.Trained{3}), 'CompactClassificationEnsemble');

%!test  # MATLAB parity: a holdout leaves the training rows unscored
%! load fisheriris
%! rng (3);
%! M = fitcensemble (meas, species, 'Method', 'Bag', 'NumLearningCycles', 3, ...
%!                   'Holdout', 0.3);
%! assert_equal (class (M), 'ClassificationPartitionedEnsemble');
%! assert_equal (M.KFold, 1);
%! [~, s] = kfoldPredict (M);
%! assert_equal (sum (any (isnan (s), 2)), sum (training (M.Partition, 1)));

## Test input validation
%!error<ClassificationPartitionedEnsemble: too few input arguments.> ...
%! ClassificationPartitionedEnsemble (1)
%!error<ClassificationPartitionedEnsemble: MDL must be a 'ClassificationEnsemble' object.> ...
%! ClassificationPartitionedEnsemble (1, c)
%!error<ClassificationPartitionedEnsemble: PARTITION must be a 'cvpartition' object.> ...
%! ClassificationPartitionedEnsemble (CV.Trainable{1}, 1)
%!error<ClassificationPartitionedEnsemble: PARTITION must partition the observations MDL was fitted on.> ...
%! ClassificationPartitionedEnsemble (CV.Trainable{1}, c)
%!error<ClassificationPartitionedEnsemble.kfoldLoss: optional arguments must be given in Name-Value pairs.> ...
%! kfoldLoss (CV, 'Mode')
%!error<ClassificationPartitionedEnsemble.kfoldLoss: invalid parameter name in optional pair arguments.> ...
%! kfoldLoss (CV, 'Learners', 1)
%!error<ClassificationPartitionedEnsemble.kfoldLoss: 'Mode' must be 'average', 'individual', or 'cumulative'.> ...
%! kfoldLoss (CV, 'Mode', 'ensemble')
%!error<ClassificationPartitionedEnsemble.kfoldLoss: 'Folds' must hold integers between 1 and 5.> ...
%! kfoldLoss (CV, 'Folds', 6)
%!error<ClassificationPartitionedEnsemble.kfoldLoss: 'LossFun' must be 'binodeviance', 'classifcost', 'classiferror', 'exponential', 'hinge', 'logit', 'mincost', 'quadratic', or a function handle.> ...
%! kfoldLoss (CV, 'LossFun', 'mse')
%!error<ClassificationPartitionedEnsemble.kfoldEdge: invalid parameter name in optional pair arguments.> ...
%! kfoldEdge (CV, 'LossFun', 'hinge')
%!error<ClassificationPartitionedEnsemble.kfoldfun: too few input arguments.> ...
%! kfoldfun (CV)
%!error<ClassificationPartitionedEnsemble.kfoldfun: FUN must be a function handle.> ...
%! kfoldfun (CV, 1)
%!error<ClassificationPartitionedEnsemble.resume: too few input arguments.> ...
%! resume (CV)
%!error<ClassificationPartitionedEnsemble.subsasgn: 'ScoreTransform' must be a character vector or a 'function_handle' object.> ...
%! M = CV;
%! M.ScoreTransform = 1;
