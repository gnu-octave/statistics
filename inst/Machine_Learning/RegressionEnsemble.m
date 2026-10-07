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

classdef RegressionEnsemble < PredictiveModel
  ## -*- texinfo -*-
  ## @deftp {statistics} RegressionEnsemble
  ##
  ## Boosted ensemble of regression trees
  ##
  ## A @code{RegressionEnsemble} object holds the regression trees LSBoost
  ## grew one after another, each fitted to the residual the trees before it
  ## left, together with the data it was fitted on.
  ##
  ## Create one with @code{fitrensemble}.  A bagged ensemble is a
  ## @code{RegressionBaggedEnsemble}, and @code{compact} returns a
  ## @code{CompactRegressionEnsemble} without the data.
  ##
  ## @seealso{fitrensemble, CompactRegressionEnsemble,
  ## RegressionBaggedEnsemble}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} Regularization
    ##
    ## Result of regularizing the ensemble
    ##
    ## Empty until @code{regularize} fills it with a structure of lasso
    ## weights for the trees, and emptied again by @code{resume}.  This
    ## property is read-only.
    ##
    ## @end deftp
    Regularization = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} Y
    ##
    ## Response data
    ##
    ## The response the ensemble was fitted on, a row missing a value having
    ## been left out.  This property is read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} X
    ##
    ## Predictor data
    ##
    ## The predictors the ensemble was fitted on, one row per observation.
    ##
    ## Where the model was fitted from a table, the predictors are the coded
    ## matrix and not the table: a variable holding levels is stored as its
    ## level codes, and the coding is kept with the model.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} RowsUsed
    ##
    ## Rows of the data that were used
    ##
    ## A logical column over the rows as supplied.  This property is
    ## read-only.
    ##
    ## @end deftp
    RowsUsed = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} W
    ##
    ## Observation weights
    ##
    ## The weights given, normalized to sum to one.  It has the class of the
    ## @qcode{'Weights'} given, single or double.  This property is read-only.
    ##
    ## @end deftp
    W = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} ModelParameters
    ##
    ## Parameters of the fit
    ##
    ## A structure with the fields @qcode{Type}, @qcode{Method},
    ## @qcode{LearnerTemplates}, @qcode{NLearn}, the number of learning
    ## cycles asked for in all, and for LSBoost @qcode{LearnRate}.  This
    ## property is read-only.
    ##
    ## @end deftp
    ModelParameters = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} NumObservations
    ##
    ## Number of observations used
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations = 0;

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## Always empty, binning not being implemented.  This property is
    ## read-only.
    ##
    ## @end deftp
    BinEdges = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} HyperparameterOptimizationResults
    ##
    ## Results of optimizing the hyperparameters
    ##
    ## Always empty, such optimization not being implemented.  This property
    ## is read-only.
    ##
    ## @end deftp
    HyperparameterOptimizationResults = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} CategoricalPredictors
    ##
    ## Indices of categorical predictors
    ##
    ## The predictors every tree treats as categorical, empty when none
    ## is.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName = 'Y';

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the learners saw them
    ##
    ## The same as @code{PredictorNames}.  This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} Method
    ##
    ## Ensemble method
    ##
    ## @qcode{'LSBoost'} or @qcode{'Bag'}.  This property is read-only.
    ##
    ## @end deftp
    Method = '';

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} LearnerNames
    ##
    ## Names of the weak learners
    ##
    ## Always @code{@{'Tree'@}}.  This property is read-only.
    ##
    ## @end deftp
    LearnerNames = {'Tree'};

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} ReasonForTermination
    ##
    ## Why the fit stopped adding trees
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ReasonForTermination = '';

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} FitInfo
    ##
    ## Fit information
    ##
    ## For LSBoost, a column with the weighted mean squared error of each
    ## tree against the residual it was fitted to.  Empty for Bag.  This
    ## property is read-only.
    ##
    ## @end deftp
    FitInfo = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} FitInfoDescription
    ##
    ## Description of @code{FitInfo}
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    FitInfoDescription = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} UsePredForLearner
    ##
    ## Which predictors each learner uses
    ##
    ## Always empty, as MATLAB returns it for tree learners.  This property
    ## is read-only.
    ##
    ## @end deftp
    UsePredForLearner = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} NumTrained
    ##
    ## Number of trained trees
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumTrained = 0;

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} Trained
    ##
    ## Trained trees
    ##
    ## A column cell array of @code{CompactRegressionTree} objects.  This
    ## property is read-only.
    ##
    ## @end deftp
    Trained = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} TrainedWeights
    ##
    ## Weights of the trained trees
    ##
    ## A column with one weight per tree.  This property is read-only.
    ##
    ## @end deftp
    TrainedWeights = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} CombineWeights
    ##
    ## How the trees are combined
    ##
    ## @qcode{'WeightedSum'} for LSBoost, @qcode{'WeightedAverage'} for Bag.
    ## This property is read-only.
    ##
    ## @end deftp
    CombineWeights = 'WeightedSum';

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {RegressionEnsemble} {property} ResponseTransform
    ##
    ## Transform applied to the predicted response
    ##
    ## See @code{CompactRegressionEnsemble.ResponseTransform}.
    ##
    ## @end deftp
    ResponseTransform = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    RTfun = @(y) y;      # the response transform as a function
    LearnRate = 1;       # shrinkage of LSBoost
    TreeArgs = {};       # the options given in the tree template
    F = [];              # what LSBoost has predicted so far, for resume
    BagFResample = 1;    # share of the observations each bag draws
    BagReplace = true;   # whether the bags draw with replacement
    BagInBag = [];       # NxNumTrained logical, the rows each bag drew
    Resampling = false;  # whether LSBoost resamples its rows
  endproperties

  methods (Hidden)

    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = ...
        parseResponseTransform (val, 'RegressionEnsemble');
    endfunction

    function disp (this)
      fprintf ("\n  %s\n\n", class (this));
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %s\n", 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ("%+25s: '%s'\n", 'ResponseTransform', this.ResponseTransform);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumTrained', this.NumTrained);
      fprintf ("%+25s: '%s'\n", 'Method', this.Method);
      fprintf ("%+25s: '%s'\n", 'ReasonForTermination', ...
               this.ReasonForTermination);
      fprintf ("%+25s: [%dx%d double]\n", 'FitInfo', rows (this.FitInfo), ...
               columns (this.FitInfo));
      if (isa (this, 'RegressionBaggedEnsemble'))
        fprintf ("%+25s: %g\n", 'FResample', this.FResample);
        fprintf ("%+25s: %d\n", 'Replace', this.Replace);
        fprintf ("%+25s: [%dx%d logical]\n", 'UseObsForLearner', ...
                 rows (this.UseObsForLearner), ...
                 columns (this.UseObsForLearner));
      endif
      fprintf ("\n");
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{obj} =} RegressionEnsemble (@var{X}, @var{Y})
    ## @deftypefnx {RegressionEnsemble} {@var{obj} =} RegressionEnsemble (@var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {RegressionEnsemble} {@var{obj} =} RegressionEnsemble (@var{Tbl}, @var{formula})
    ## @deftypefnx {RegressionEnsemble} {@var{obj} =} RegressionEnsemble (@var{Tbl}, @var{Y})
    ## @deftypefnx {RegressionEnsemble} {@var{obj} =} RegressionEnsemble (@dots{}, @var{name}, @var{value})
    ##
    ## Fit an ensemble of regression trees by LSBoost.
    ##
    ## @code{fitrensemble} is the documented way in, and its help lists the
    ## options both take.  A bagged ensemble is fitted by
    ## @code{RegressionBaggedEnsemble}.
    ##
    ## @seealso{fitrensemble, RegressionBaggedEnsemble}
    ## @end deftypefn
    function this = RegressionEnsemble (X, Y, varargin)

      if (nargin < 2)
        error ("RegressionEnsemble: too few input arguments.");
      endif

      ## A table names its own predictors and says which hold levels
      [this, X, Y, varargin] = resolveTable (this, 'RegressionEnsemble', ...
                                             X, Y, varargin);
      if (mod (numel (varargin), 2) != 0)
        error ("RegressionEnsemble: name-value arguments must be in pairs.");
      endif
      bagged = isa (this, 'RegressionBaggedEnsemble');
      if (bagged)
        caller = 'RegressionBaggedEnsemble';
      else
        caller = 'RegressionEnsemble';
      endif

      ## Parse optional paired arguments.  'LearnRate', 'FResample' and
      ## 'Replace' stay empty when not given, which is how the checks below
      ## tell them apart from a value, and 'Resample' is off.
      optNames = {'Method', 'NumLearningCycles', 'Learners', 'LearnRate', ...
                  'NPrint', 'Weights', 'PredictorNames', 'ResponseName', ...
                  'ResponseTransform', 'FResample', 'Replace', 'Resample', ...
                  'CategoricalPredictors'};
      dfValues = {'LSBoost', 100, 'tree', [], 'off', [], {}, 'Y', 'none', ...
                  [], [], [], []};
      [Method, NLearn, Learners, LearnRate, NPrint, Weights, ...
       PredictorNames, ResponseName, ResponseTransform, FResample, ...
       Replace, Resample, CatPreds, args] = ...
                 parsePairedArguments (optNames, dfValues, varargin(:));

      ## Validate optional paired arguments
      if (! (ischar (Method) && isrow (Method)))
        error ("%s: 'Method' must be a character vector.", caller);
      elseif (strcmpi (Method, 'LSBoost'))
        Method = 'LSBoost';
      elseif (strcmpi (Method, 'Bag'))
        Method = 'Bag';
      else
        error ("%s: '%s' is not a valid ensemble method.", caller, Method);
      endif
      if (! (isnumeric (NLearn) && isscalar (NLearn) && isreal (NLearn)
             && NLearn >= 1 && NLearn == fix (NLearn)))
        error ("%s: 'NumLearningCycles' must be a positive integer.", caller);
      endif
      NLearn = double (NLearn);
      if (! isempty (LearnRate)
          && ! (isnumeric (LearnRate) && isscalar (LearnRate)
                && isreal (LearnRate) && LearnRate > 0 && LearnRate <= 1))
        error (strcat ("%s: 'LearnRate' must be a number greater than 0", ...
                       " and no greater than 1."), caller);
      endif
      LearnRate = double (LearnRate);
      if (ischar (NPrint) && strcmpi (NPrint, 'off'))
        NPrint = 0;
      elseif (isnumeric (NPrint) && isscalar (NPrint) && isreal (NPrint)
              && NPrint >= 1 && NPrint == fix (NPrint))
        NPrint = double (NPrint);
      else
        error ("%s: 'NPrint' must be a positive integer or 'off'.", caller);
      endif
      if (! (ischar (ResponseName) && isrow (ResponseName)))
        error ("%s: 'ResponseName' must be a character vector.", caller);
      endif
      if (! isempty (FResample)
          && ! (isnumeric (FResample) && isscalar (FResample)
                && isreal (FResample) && FResample > 0 && FResample <= 1))
        error (strcat ("%s: 'FResample' must be a number greater than 0", ...
                       " and no greater than 1."), caller);
      endif
      FResample = double (FResample);
      if (! isempty (Replace))
        [Replace, ok] = onOff (Replace);
        if (! ok)
          error ("%s: 'Replace' must be 'on' or 'off'.", caller);
        endif
      endif
      if (isempty (Resample))
        Resample = false;
      else
        [Resample, ok] = onOff (Resample);
        if (! ok)
          error ("%s: 'Resample' must be 'on' or 'off'.", caller);
        endif
      endif

      ## Options MATLAB takes that this class does not implement are named
      ## one by one, so that asking for one is refused rather than quietly
      ## doing nothing; anything else left over is unknown.
      notImpl = {'NumBins', 'OptimizeHyperparameters', ...
                 'HyperparameterOptimizationOptions', 'Options'};
      for i = 1:2:numel (args)
        if (ischar (args{i}) && any (strcmpi (args{i}, notImpl)))
          error ("%s: '%s' is not implemented.", caller, args{i});
        endif
      endfor
      if (! isempty (args))
        error ("%s: invalid optional paired argument.", caller);
      endif

      if (ischar (Learners) && strcmpi (Learners, 'tree'))
        tmpl = templateTree ();
      elseif (isstruct (Learners) && isscalar (Learners)
              && isfield (Learners, 'Method')
              && strcmpi (Learners.Method, 'Tree'))
        tmpl = Learners;
      else
        error ("%s: 'Learners' must be 'tree' or a tree template.", caller);
      endif
      tmpl.Type = 'regression';
      TreeArgs = {};
      for [val, key] = tmpl
        if (! any (strcmp (key, {'Method', 'Type'})))
          TreeArgs(end+1:end+2) = {key, val};
        endif
      endfor

      isbag = strcmp (Method, 'Bag');
      resampled = Resample || ! isempty (FResample) || ! isempty (Replace);
      if (isbag && ! bagged)
        error (strcat ("RegressionEnsemble: a bagged ensemble is fitted by", ...
                       " RegressionBaggedEnsemble."));
      elseif (! isbag && bagged && ! resampled)
        error (strcat ("RegressionBaggedEnsemble: 'Method' must be 'Bag'", ...
                       " unless the ensemble resamples."));
      elseif (! isbag && ! bagged && resampled)
        error (strcat ("RegressionEnsemble: a resampled ensemble is fitted", ...
                       " by RegressionBaggedEnsemble."));
      endif
      if (isbag)
        if (! isempty (LearnRate))
          error ("%s: 'LearnRate' cannot be used with the 'Bag' method.", ...
                 caller);
        endif
        if (isempty (FResample))
          FResample = 1;
        endif
        if (isempty (Replace))
          Replace = true;
        endif
      else
        if (resampled && isempty (FResample))
          FResample = 1;
        endif
        if (resampled && isempty (Replace))
          Replace = true;
        endif
        if (isempty (LearnRate))
          LearnRate = 1;
        endif
      endif

      F = regFrame (X, Y, Weights, caller);
      if (any (F.Weights < 0) || ! (sum (F.Weights) > 0))
        error (strcat ("%s: 'Weights' must be nonnegative and not all", ...
                       " zero."), caller);
      endif
      if (isempty (PredictorNames))
        PredictorNames = arrayfun (@(k) sprintf ('x%d', k), 1:F.p, ...
                                   'UniformOutput', false);
      elseif (! (iscellstr (PredictorNames) && numel (PredictorNames) == F.p))
        error (strcat ("%s: 'PredictorNames' must be a cell array of", ...
                       " character vectors with one element per column", ...
                       " of X."), caller);
      endif
      ## Categorical predictors reach every tree
      [Cod, errmsg] = dummyCoding (F.X, CatPreds, PredictorNames);
      if (! isempty (errmsg))
        error ("%s: %s", caller, errmsg);
      endif
      CatIdx = [];
      if (! isempty (Cod.Index))
        CatIdx = Cod.Index;
        TreeArgs = [{'CategoricalPredictors', CatIdx}, TreeArgs];
      endif

      this.X = F.X;
      this.Y = double (F.Y);
      this.RowsUsed = F.RowsUsed;
      ## The weights keep their class in the model; every computation runs on
      ## them as double.
      this.W = cast (F.W, F.WeightsClass);
      this.NumObservations = F.n;
      this.PredictorNames = PredictorNames(:)';
      this.ExpandedPredictorNames = this.PredictorNames;
      this.ResponseName = ResponseName;
      this.ResponseTransform = ResponseTransform;
      this.Method = Method;
      this.TreeArgs = TreeArgs;
      this.CategoricalPredictors = CatIdx;
      this.Trained = cell (0, 1);
      this.TrainedWeights = zeros (0, 1);
      this.ModelParameters = struct ('Type', 'regression', ...
                                     'Method', Method, ...
                                     'LearnerTemplates', tmpl, ...
                                     'NLearn', 0);
      if (isbag)
        this.CombineWeights = 'WeightedAverage';
        this.FitInfo = [];
        this.FitInfoDescription = 'None';
        this.BagFResample = FResample;
        this.BagReplace = Replace;
        this.BagInBag = false (F.n, 0);
      else
        this.LearnRate = LearnRate;
        this.ModelParameters.LearnRate = LearnRate;
        this.Resampling = resampled;
        if (resampled)
          this.BagFResample = FResample;
          this.BagReplace = Replace;
          this.BagInBag = false (F.n, 0);
        endif
        this.FitInfo = zeros (0, 1);
        this.FitInfoDescription = ...
          {strcat("Vector of length NumTrained, where NumTrained is the", ...
                  " number of learned weak hypotheses."); ...
           strcat("Element t of this vector is the weighted residual from", ...
                  " learner t.")};
        this.F = zeros (F.n, 1);
      endif

      this = growLearners (this, NLearn, NPrint);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {RegressionEnsemble} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate an ensemble.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} refits the ensemble on the
    ## training part of each of ten folds, and returns a
    ## @code{RegressionPartitionedEnsemble}.  One of @qcode{'KFold'}, an
    ## integer greater than 1, @qcode{'Holdout'}, a number between 0 and 1,
    ## @qcode{'Leaveout'}, @qcode{'on'} for one fold per observation, or
    ## @qcode{'CVPartition'}, a @code{cvpartition} object, may choose the
    ## partition instead.
    ##
    ## @seealso{RegressionEnsemble, RegressionPartitionedEnsemble, cvpartition}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)

      [P, errmsg] = ensemblePartition (varargin, this.Y, ...
                                       this.NumObservations, ...
                                       false);
      if (! isempty (errmsg))
        error ("%s.crossval: %s", class (this), errmsg);
      endif
      if (isempty (P))
        error (strcat ("%s.crossval: no partition was asked for; set", ...
                       " 'CrossVal' to 'on' or give 'KFold'."), class (this));
      endif
      CVMdl = RegressionPartitionedEnsemble (this, P);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionEnsemble} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Drop the training data from a regression ensemble.
    ##
    ## @code{@var{CMdl} = compact (@var{obj})} returns a
    ## @code{CompactRegressionEnsemble} holding the trees and what prediction
    ## needs.  It predicts new data identically.
    ##
    ## @seealso{RegressionEnsemble, CompactRegressionEnsemble}
    ## @end deftypefn
    function CMdl = compact (this)

      CMdl = CompactRegressionEnsemble (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{B} =} resume (@var{obj}, @var{NumLearningCycles})
    ## @deftypefnx {RegressionEnsemble} {@var{B} =} resume (@dots{}, 'NPrint', @var{n})
    ##
    ## Grow more trees.
    ##
    ## @var{B} is the ensemble with @var{NumLearningCycles} further trees
    ## grown as though the fit had asked for them from the start.
    ## @qcode{'NPrint'} is taken as by @code{fitrensemble}.  A
    ## @code{Regularization} is emptied, its weights no longer covering every
    ## tree; MATLAB keeps it and applies its weights to the first trees.
    ##
    ## @seealso{RegressionEnsemble, fitrensemble}
    ## @end deftypefn
    function this = resume (this, NumLearningCycles, varargin)

      caller = sprintf ('%s.resume', class (this));
      if (nargin < 2)
        error ("%s: too few input arguments.", caller);
      endif
      if (! (isnumeric (NumLearningCycles) && isscalar (NumLearningCycles)
             && isreal (NumLearningCycles) && NumLearningCycles >= 1
             && NumLearningCycles == fix (NumLearningCycles)))
        error ("%s: NUMLEARNINGCYCLES must be a positive integer.", caller);
      endif
      if (mod (numel (varargin), 2) != 0)
        error ("%s: name-value arguments must be in pairs.", caller);
      endif
      ## Parse optional paired arguments
      [NPrint, args] = parsePairedArguments ({'NPrint'}, {'off'}, ...
                                             varargin(:));

      ## Validate optional paired arguments
      if (ischar (NPrint) && strcmpi (NPrint, 'off'))
        NPrint = 0;
      elseif (isnumeric (NPrint) && isscalar (NPrint) && isreal (NPrint)
              && NPrint >= 1 && NPrint == fix (NPrint))
        NPrint = double (NPrint);
      else
        error ("%s: 'NPrint' must be a positive integer or 'off'.", caller);
      endif

      if (! isempty (args))
        error ("%s: invalid optional paired argument.", caller);
      endif
      this = growLearners (this, double (NumLearningCycles), NPrint);
      this.Regularization = [];

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{yfit} =} predict (@var{obj}, @var{X})
    ## @deftypefnx {RegressionEnsemble} {@var{yfit} =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict the response with a regression ensemble.
    ##
    ## Behaves as @code{CompactRegressionEnsemble.predict} and takes the same
    ## Name-Value arguments.
    ##
    ##
    ## The new data may be a table, whose variables are matched to the
    ## predictors the model was fitted on by name and not by position:
    ## one the model was not fitted on is passed over, one it needs and
    ## cannot find is named, and a value holding a level is coded as that
    ## level was coded at fitting.
    ## @seealso{RegressionEnsemble, CompactRegressionEnsemble.predict}
    ## @end deftypefn
    function yfit = predict (this, X, varargin)

      if (nargin < 2)
        error ("%s.predict: too few input arguments.", class (this));
      endif

      ## A table is read by the names the model was fitted on
      X = tableColumns (this, 'RegressionEnsemble.predict', X);
      yfit = ensemblePredict (compact (this), X, varargin, ...
                              [class(this), '.predict']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionEnsemble} {@var{L} =} loss (@var{obj}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {RegressionEnsemble} {@var{L} =} loss (@var{obj}, @var{Tbl})
    ## @deftypefnx {RegressionEnsemble} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss of an ensemble.
    ##
    ## Behaves as @code{CompactRegressionEnsemble.loss} and takes the same
    ## Name-Value arguments.
    ##
    ## @var{X} may also be a table @var{Tbl}, whose variables are matched to
    ## the predictors the model was fitted on by name and not by position.
    ## @code{loss (@var{obj}, @var{Tbl}, @var{ResponseVarName})} takes the
    ## response from the variable @var{ResponseVarName} names, and
    ## @code{loss (@var{obj}, @var{Tbl})} from the variable the model was
    ## fitted on.  The response may also be given beside the table as
    ## @var{Y}.
    ##
    ## @seealso{RegressionEnsemble, CompactRegressionEnsemble.loss}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      if (nargin < 3 && ! (nargin > 1 && istable (X)))
        error ("%s.loss: too few input arguments.", class (this));
      endif

      ## A table carries the response: named in the call, given beside
      ## the table, or the variable the model was fitted on
      if (nargin < 3)
        Y = [];
      endif
      [X, Y, varargin] = tableResponse (this, 'loss', X, Y, varargin, ...
                                        nargin > 2);
      L = ensembleLoss (compact (this), X, Y, varargin, ...
                        [class(this), '.loss']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{imp} =} predictorImportance (@var{obj})
    ## @deftypefnx {RegressionEnsemble} {[@var{imp}, @var{ma}] =} predictorImportance (@var{obj})
    ##
    ## Estimate the importance of each predictor.
    ##
    ## Behaves as @code{CompactRegressionEnsemble.predictorImportance}.
    ##
    ## @seealso{RegressionEnsemble,
    ## CompactRegressionEnsemble.predictorImportance}
    ## @end deftypefn
    function [imp, ma] = predictorImportance (this)

      [imp, ma] = ensembleImportance (compact (this));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{yfit} =} resubPredict (@var{obj})
    ## @deftypefnx {RegressionEnsemble} {@var{yfit} =} resubPredict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict the response of the training data.
    ##
    ## @code{predict} on @code{X}, taking the same Name-Value arguments.
    ##
    ## @seealso{RegressionEnsemble, RegressionEnsemble.predict}
    ## @end deftypefn
    function yfit = resubPredict (this, varargin)

      yfit = ensemblePredict (compact (this), this.X, varargin, ...
                              [class(this), '.resubPredict']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {RegressionEnsemble} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss on the training data.
    ##
    ## @code{loss} on @code{X} and @code{Y}, weighted by @code{W} unless
    ## @qcode{'Weights'} are given.
    ##
    ## @seealso{RegressionEnsemble, RegressionEnsemble.loss}
    ## @end deftypefn
    function L = resubLoss (this, varargin)

      L = ensembleLoss (compact (this), this.X, this.Y, ...
                        [{'Weights', this.W}, varargin], ...
                        [class(this), '.resubLoss']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{B} =} regularize (@var{obj})
    ## @deftypefnx {RegressionEnsemble} {@var{B} =} regularize (@dots{}, @var{name}, @var{value})
    ##
    ## Find lasso weights for the trees of a regression ensemble.
    ##
    ## @code{@var{B} = regularize (@var{obj})} fits the training response by
    ## the trees' predictions with a lasso that has no intercept and whose
    ## weights may not be negative, over a path of penalties, and returns the
    ## ensemble with the result in @code{Regularization}.  For a penalty
    ## @var{lambda} the tree weights @var{a} minimize
    ## @code{sum (W .* (Y - P * a) .^ 2) / 2 + lambda * sum (a)}, where @var{P}
    ## holds one column of training predictions per tree and @code{W} is the
    ## observation weights, which sum to one.  @code{TrainedWeights} is left as
    ## it was.
    ##
    ## @multitable @columnfractions 0.2 0.02 0.78
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Lambda'} @tab @tab A vector of non-negative penalties.  The
    ## default is 0 followed by nine values spaced evenly on a log scale from
    ## @code{lmax / 1000} to @code{lmax}, the smallest penalty that sets every
    ## weight to zero, @code{lmax = max (abs (P' * (W .* Y)))}.
    ## @item @qcode{'MaxIter'} @tab @tab The most passes of coordinate descent
    ## over the trees for each penalty.  The default is 1e5.
    ## @item @qcode{'RelTol'} @tab @tab The descent stops after a pass in which
    ## no weight changed by more than this times the largest weight, or this
    ## when that weight is below one.  The default is 1e-10.
    ## @end multitable
    ##
    ## @code{Regularization} is a structure with the fields @qcode{Method},
    ## @qcode{'Lasso'}; @qcode{TrainedWeights}, one column per penalty;
    ## @qcode{Lambda}; @qcode{ResubstitutionMSE}, the weighted mean squared
    ## error of each column; and @qcode{CombineWeights},
    ## @qcode{'WeightedSum'}.
    ##
    ## MATLAB's solver can stop well short of the minimum at small penalties,
    ## so its weights there differ from these, which are the minimum.  Its
    ## @qcode{CombineWeights} is a function handle, and it also takes
    ## @qcode{'Npass'} and @qcode{'Verbose'}, which are not taken here.
    ##
    ## @seealso{RegressionEnsemble, RegressionEnsemble.shrink,
    ## RegressionEnsemble.cvshrink, lasso}
    ## @end deftypefn
    function this = regularize (this, varargin)

      caller = sprintf ('%s.regularize', class (this));
      [o, errmsg] = shrinkOptions (varargin, {'Lambda', 'MaxIter', 'RelTol'});
      if (! isempty (errmsg))
        error ("%s: %s", caller, errmsg);
      endif

      T = this.NumTrained;
      Y = this.Y;
      W = double (this.W) / sum (double (this.W));
      P = zeros (numel (Y), T);
      for t = 1:T
        P(:,t) = predict (this.Trained{t}, this.X);
      endfor
      Lambda = o.Lambda;
      if (isempty (Lambda))
        lmax = 0;
        if (T > 0)
          lmax = max (abs (P' * (W .* Y)));
        endif
        Lambda = [0, lmax * 10 .^ linspace(-3, 0, 9)];
      endif

      ## From the largest penalty down, each starting from the last answer.
      L = numel (Lambda);
      TW = zeros (T, L);
      mse = zeros (1, L);
      [~, order] = sort (Lambda, 'descend');
      w = zeros (T, 1);
      for k = order
        w = nnLasso (P, Y, W, Lambda(k), w, o.MaxIter, o.RelTol);
        TW(:,k) = w;
        mse(k) = sum (W .* (Y - P * w) .^ 2);
      endfor
      this.Regularization = struct ('Method', 'Lasso', ...
                                    'TrainedWeights', TW, ...
                                    'Lambda', Lambda, ...
                                    'ResubstitutionMSE', mse, ...
                                    'CombineWeights', 'WeightedSum');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{C} =} shrink (@var{obj})
    ## @deftypefnx {RegressionEnsemble} {@var{C} =} shrink (@dots{}, @var{name}, @var{value})
    ##
    ## Keep the trees a lasso weight retains.
    ##
    ## @code{@var{C} = shrink (@var{obj})} returns a
    ## @code{CompactRegressionEnsemble} of the trees whose weight in a column
    ## of @code{Regularization.TrainedWeights} is above a threshold, ordered
    ## from the largest weight down, each carrying that weight and their
    ## predictions summed.  An ensemble that has not been regularized is
    ## thresholded on its @code{TrainedWeights} and keeps its way of combining
    ## its trees.
    ##
    ## @multitable @columnfractions 0.2 0.02 0.78
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'WeightColumn'} @tab @tab The column of weights to use, a
    ## positive integer.  The default is 1.
    ## @item @qcode{'Threshold'} @tab @tab A non-negative number; a tree whose
    ## weight is not above it is dropped.  The default is 0.
    ## @item @qcode{'Lambda'} @tab @tab Penalties to regularize with first, as
    ## @code{regularize} does, which also takes @qcode{'MaxIter'} and
    ## @qcode{'RelTol'} here.
    ## @end multitable
    ##
    ## MATLAB accepts a @qcode{'WeightColumn'} that is not a whole number,
    ## which is refused here.
    ##
    ## @seealso{RegressionEnsemble, RegressionEnsemble.regularize,
    ## RegressionEnsemble.cvshrink}
    ## @end deftypefn
    function C = shrink (this, varargin)

      caller = sprintf ('%s.shrink', class (this));
      [o, errmsg] = shrinkOptions (varargin, {'WeightColumn', 'Threshold', ...
                                              'Lambda', 'MaxIter', 'RelTol'});
      if (! isempty (errmsg))
        error ("%s: %s", caller, errmsg);
      endif
      if (! isscalar (o.Threshold))
        error ("%s: 'Threshold' must be a scalar.", caller);
      endif
      if (! isempty (o.Lambda))
        this = regularize (this, 'Lambda', o.Lambda, 'MaxIter', o.MaxIter, ...
                           'RelTol', o.RelTol);
      endif

      C = compact (this);
      if (isempty (this.Regularization))
        weights = this.TrainedWeights(:);
        combine = this.CombineWeights;
      else
        TW = this.Regularization.TrainedWeights;
        if (o.WeightColumn > columns (TW))
          error (strcat ("%s: 'WeightColumn' must not exceed the number of", ...
                         " values in 'Lambda'."), caller);
        endif
        weights = TW(:,o.WeightColumn);
        combine = 'WeightedSum';
      endif
      keep = find (weights > o.Threshold);
      [~, order] = sort (weights(keep), 'descend');
      idx = keep(order);
      C = keepLearners (C, idx, weights(idx), combine);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {[@var{vals}, @var{nlearn}] =} cvshrink (@var{obj})
    ## @deftypefnx {RegressionEnsemble} {[@var{vals}, @var{nlearn}] =} cvshrink (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate the shrinking of a regression ensemble.
    ##
    ## @code{[@var{vals}, @var{nlearn}] = cvshrink (@var{obj})} grows the
    ## ensemble again on the training part of each fold, as @code{crossval}
    ## does, regularizes it with each penalty, shrinks it at each threshold,
    ## and predicts the fold's held-out observations.  @var{vals} holds one row
    ## per penalty and one column per threshold: the weighted mean squared
    ## error pooled over every held-out observation.  @var{nlearn} holds the
    ## matching mean number of trees kept per fold.
    ##
    ## @multitable @columnfractions 0.2 0.02 0.78
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Lambda'} @tab @tab A vector of non-negative penalties.  The
    ## default is @code{Regularization.Lambda}; an ensemble that has not been
    ## regularized must be given one.
    ## @item @qcode{'Threshold'} @tab @tab A vector of non-negative thresholds,
    ## as @code{shrink} takes them.  The default is 0.
    ## @item @qcode{'KFold'}, @qcode{'Holdout'}, @qcode{'Leaveout'},
    ## @qcode{'CVPartition'} @tab @tab The partition, as @code{crossval} takes
    ## it, only one of them.  The default is ten folds.
    ## @item @qcode{'MaxIter'}, @qcode{'RelTol'} @tab @tab As @code{regularize}
    ## takes them.
    ## @end multitable
    ##
    ## MATLAB warns and returns empty outputs when an ensemble that has not
    ## been regularized is given no penalties; here that is an error.
    ##
    ## @seealso{RegressionEnsemble, RegressionEnsemble.shrink,
    ## RegressionEnsemble.crossval}
    ## @end deftypefn
    function [vals, nlearn] = cvshrink (this, varargin)

      caller = sprintf ('%s.cvshrink', class (this));
      if (mod (numel (varargin), 2) != 0)
        error ("%s: name-value arguments must be in pairs.", caller);
      endif
      cvnames = {'kfold', 'holdout', 'leaveout', 'cvpartition'};
      cvargs = {};
      rest = {};
      for i = 1:2:numel (varargin)
        if (ischar (varargin{i}) && any (strcmpi (varargin{i}, cvnames)))
          cvargs(end+1:end+2) = varargin(i:i+1);
        else
          rest(end+1:end+2) = varargin(i:i+1);
        endif
      endfor
      [o, errmsg] = shrinkOptions (rest, {'Lambda', 'Threshold', ...
                                          'MaxIter', 'RelTol'});
      if (! isempty (errmsg))
        error ("%s: %s", caller, errmsg);
      endif
      Lambda = o.Lambda;
      if (isempty (Lambda))
        if (isempty (this.Regularization))
          error (strcat ("%s: 'Lambda' must be given for an ensemble that", ...
                         " has not been regularized."), caller);
        endif
        Lambda = this.Regularization.Lambda;
      endif
      [P, errmsg] = ensemblePartition (cvargs, this.Y, this.NumObservations, ...
                                       false);
      if (! isempty (errmsg))
        error ("%s: %s", caller, errmsg);
      endif

      CV = RegressionPartitionedEnsemble (this, P);
      L = numel (Lambda);
      M = numel (o.Threshold);
      sse = zeros (L, M);
      counts = zeros (L, M);
      wsum = 0;
      K = numel (CV.Trainable);
      for k = 1:K
        te = test (CV.Partition, k);
        Ek = regularize (CV.Trainable{k}, 'Lambda', Lambda, ...
                         'MaxIter', o.MaxIter, 'RelTol', o.RelTol);
        wsum += sum (double (this.W(te)));
        for a = 1:L
          for b = 1:M
            Ck = shrink (Ek, 'WeightColumn', a, 'Threshold', o.Threshold(b));
            r = this.Y(te) - predict (Ck, this.X(te,:));
            sse(a,b) += sum (double (this.W(te)) .* r .^ 2);
            counts(a,b) += Ck.NumTrained;
          endfor
        endfor
      endfor
      vals = sse / wsum;
      nlearn = counts / K;

    endfunction

  endmethods

  methods (Access = protected)

    ## Grow N trees, LSBoost carrying what it has predicted so far in F so
    ## that a resume continues the fit exactly.
    function this = growLearners (this, N, NPrint)

      n = this.NumObservations;
      p = columns (this.X);
      this.ModelParameters.NLearn += N;
      if (NPrint > 0)
        printf ("Training %s...\n", this.Method);
      endif
      for c = 1:N
        if (strcmp (this.Method, 'LSBoost'))
          r = this.Y - this.F;
          targs = {'PredictorNames', this.PredictorNames, ...
                   'MaxNumSplits', 10, 'MinParentSize', 10, ...
                   'MinLeafSize', 5, 'Prune', 'off', 'MergeLeaves', 'off', ...
                   this.TreeArgs{:}};
          if (this.Resampling)
            ## The tree sees the rows drawn; the prediction and the fit
            ## information run over every row, as in MATLAB R2024a.
            [idx, sw, cnt] = boostSample (double (this.W), ...
                                          ceil (this.BagFResample * n), ...
                                          this.BagReplace);
            T = compact (RegressionTree (this.X(idx,:), r(idx), ...
                                         'Weights', sw, targs{:}));
            this.BagInBag(:,end+1) = cnt > 0;
          else
            T = compact (RegressionTree (this.X, r, ...
                                         'Weights', double (this.W), ...
                                         targs{:}));
          endif
          h = predict (T, this.X);
          this.FitInfo(end+1,1) = sum (double (this.W) .* (r - h) .^ 2);
          this.F += this.LearnRate * h;
          this = addLearner (this, T, this.LearnRate);
        else
          m = ceil (this.BagFResample * n);
          if (this.BagReplace)
            cw = [0; cumsum(double (this.W))];
            cw /= cw(end);
            idx = lookup (cw, rand (m, 1));
          else
            [~, order] = sort (rand (n, 1) .^ (1 ./ double (this.W)), ...
                               'descend');
            idx = order(1:m);
          endif
          T = compact (RegressionTree (this.X(idx,:), this.Y(idx), ...
                                       'PredictorNames', ...
                                       this.PredictorNames, ...
                                       'MinParentSize', 10, ...
                                       'MinLeafSize', 5, 'Prune', 'off', ...
                                       'MergeLeaves', 'off', ...
                                       'NumVariablesToSample', ...
                                       ceil (p / 3), this.TreeArgs{:}));
          inbag = false (n, 1);
          inbag(idx) = true;
          this.BagInBag(:,end+1) = inbag;
          this = addLearner (this, T, 1);
        endif
        if (NPrint > 0 && mod (this.NumTrained, NPrint) == 0)
          printf ("Grown weak learners: %d\n", this.NumTrained);
        endif
      endfor
      this.ReasonForTermination = strcat ("Terminated normally after", ...
                                          " completing the requested", ...
                                          " number of training cycles.");

    endfunction

    function this = addLearner (this, mdl, a)
      this.Trained{end+1,1} = mdl;
      this.TrainedWeights(end+1,1) = a;
      this.NumTrained = numel (this.Trained);
    endfunction

  endmethods

endclassdef

## The options of regularize, shrink and cvshrink that ALLOWED names, checked
## and with their defaults.  ERRMSG is the body of the message the caller
## raises, or empty.
function [o, errmsg] = shrinkOptions (args, allowed)

  o = struct ('Lambda', [], 'MaxIter', 1e5, 'RelTol', 1e-10, ...
              'WeightColumn', 1, 'Threshold', 0);
  errmsg = "";
  if (mod (numel (args), 2) != 0)
    errmsg = "name-value arguments must be in pairs.";
    return;
  endif

  ## Parse optional paired arguments.  Each caller accepts the subset of
  ## them ALLOWED lists; an option left empty keeps the default above.
  optNames = {'Lambda', 'MaxIter', 'RelTol', 'WeightColumn', 'Threshold'};
  [Lambda, MaxIter, RelTol, WeightColumn, Threshold, args] = ...
         parsePairedArguments (optNames, {[], [], [], [], []}, args(:));
  given = optNames(! cellfun (@isempty, {Lambda, MaxIter, RelTol, ...
                                          WeightColumn, Threshold}));
  if (! isempty (args) || ! all (ismember (lower (given), lower (allowed))))
    errmsg = "invalid optional paired argument.";
    return;
  endif

  ## Validate optional paired arguments
  if (! isempty (Lambda))
    if (! (isnumeric (Lambda) && isreal (Lambda) && isvector (Lambda)
           && all (isfinite (Lambda)) && all (Lambda >= 0)))
      errmsg = "'Lambda' must be a vector of non-negative numbers.";
      return;
    endif
    o.Lambda = double (Lambda(:)');
  endif
  if (! isempty (MaxIter))
    if (! (isnumeric (MaxIter) && isscalar (MaxIter) && isreal (MaxIter)
           && MaxIter >= 1 && MaxIter == fix (MaxIter)))
      errmsg = "'MaxIter' must be a positive integer.";
      return;
    endif
    o.MaxIter = double (MaxIter);
  endif
  if (! isempty (RelTol))
    if (! (isnumeric (RelTol) && isscalar (RelTol) && isreal (RelTol)
           && isfinite (RelTol) && RelTol > 0))
      errmsg = "'RelTol' must be a positive number.";
      return;
    endif
    o.RelTol = double (RelTol);
  endif
  if (! isempty (WeightColumn))
    if (! (isnumeric (WeightColumn) && isscalar (WeightColumn)
           && isreal (WeightColumn) && WeightColumn >= 1
           && WeightColumn == fix (WeightColumn)))
      errmsg = "'WeightColumn' must be a positive integer.";
      return;
    endif
    o.WeightColumn = double (WeightColumn);
  endif
  if (! isempty (Threshold))
    if (! (isnumeric (Threshold) && isreal (Threshold) && isvector (Threshold)
           && all (isfinite (Threshold)) && all (Threshold >= 0)))
      errmsg = "'Threshold' must hold non-negative numbers.";
      return;
    endif
    o.Threshold = double (Threshold(:)');
  endif

endfunction

## The non-negative weights W over the columns of P that minimize
## sum (V .* (y - P * w) .^ 2) / 2 + lambda * sum (w), V summing to one, by
## cyclic coordinate descent from the given W.
function w = nnLasso (P, y, V, lambda, w, maxiter, reltol)

  VP = V .* P;
  d = sum (VP .* P, 1)';
  r = y - P * w;
  for it = 1:maxiter
    step = 0;
    for j = 1:numel (w)
      if (d(j) <= 0)
        continue;
      endif
      new = max (0, (VP(:,j)' * r + d(j) * w(j) - lambda) / d(j));
      if (new != w(j))
        r -= P(:,j) * (new - w(j));
        step = max (step, abs (new - w(j)));
        w(j) = new;
      endif
    endfor
    if (step <= reltol * max ([abs(w); 1]))
      break;
    endif
  endfor

endfunction

## 'on' or 'off' as a logical scalar, OK false for anything else.
function [tf, ok] = onOff (val)

  ok = ischar (val) && any (strcmpi (val, {'on', 'off'}));
  tf = ok && strcmpi (val, 'on');

endfunction

## Test output
%!shared X, y, S
%! load fisheriris
%! X = meas(:,2:4);
%! y = meas(:,1);
%! S = templateTree ('MaxNumSplits', 1);

%!test  # MATLAB parity: the properties of a boosted regression ensemble
%! Mdl = RegressionEnsemble (X, y, 'NumLearningCycles', 2, 'Learners', S);
%! assert_equal (numel (properties (Mdl)), 24);
%! assert_equal (Mdl.Regularization, []);
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (Mdl.LearnerNames, {'Tree'});
%! assert_equal (Mdl.UsePredForLearner, []);
%! assert_equal (Mdl.ReasonForTermination, ["Terminated normally after ", ...
%!               "completing the requested number of training cycles."]);
%! assert_equal (Mdl.FitInfoDescription{2}, ["Element t of this vector ", ...
%!               "is the weighted residual from learner t."]);

%!test  # MATLAB parity: the parameters of an LSBoost fit
%! Mdl = RegressionEnsemble (X, y, 'NumLearningCycles', 2, 'Learners', S, ...
%!                           'LearnRate', 0.5);
%! mp = Mdl.ModelParameters;
%! assert_equal (mp.Type, 'regression');
%! assert_equal (mp.Method, 'LSBoost');
%! assert_equal (mp.NLearn, 2);
%! assert_equal (mp.LearnRate, 0.5);

%!test  # MATLAB parity: resuming continues the fit exactly
%! M4 = RegressionEnsemble (X, y, 'NumLearningCycles', 4, 'Learners', S, ...
%!                          'LearnRate', 0.5);
%! M2 = RegressionEnsemble (X, y, 'NumLearningCycles', 2, 'Learners', S, ...
%!                          'LearnRate', 0.5);
%! R = resume (M2, 2);
%! assert_equal (R.NumTrained, 4);
%! assert_equal (R.ModelParameters.NLearn, 4);
%! assert_equal (R.FitInfo, M4.FitInfo, 1e-14);
%! assert_equal (R.TrainedWeights, M4.TrainedWeights);

%!test  # MATLAB parity: resubstitution with the training weights
%! Mdl = RegressionEnsemble (X, y, 'NumLearningCycles', 3, 'Learners', S);
%! assert_equal (resubLoss (Mdl), 0.176267398277276, 1e-13);
%! yf = resubPredict (Mdl);
%! assert_equal (yf(1:2), [4.953658536585365; 4.953658536585365], 1e-13);
%! w = [5 * ones(50, 1); ones(100, 1)];
%! Mw = RegressionEnsemble (X, y, 'NumLearningCycles', 2, 'Learners', S, ...
%!                          'Weights', w);
%! assert_equal (resubLoss (Mw), 0.147589707296207, 1e-13);
%! assert_equal (loss (Mw, X, y), 0.181635131105396, 1e-13);

%!test  # MATLAB parity: a missing response drops its row from the fit
%! yn = y;
%! yn(3) = NaN;
%! Mdl = RegressionEnsemble (X, yn, 'NumLearningCycles', 2, 'Learners', S);
%! assert_equal (Mdl.NumObservations, 149);
%! assert_equal (sum (Mdl.RowsUsed), 149);

%!test  # MATLAB parity: the response transform applies to predictions only
%! M0 = RegressionEnsemble (X, y, 'NumLearningCycles', 3, 'Learners', S);
%! Mt = RegressionEnsemble (X, y, 'NumLearningCycles', 3, 'Learners', S, ...
%!                          'ResponseTransform', @(z) 2 * z);
%! assert_equal (Mt.FitInfo, M0.FitInfo);
%! assert_equal (predict (Mt, X(1,:)), 9.907317073170731, 1e-13);
%! assert_equal (loss (Mt, X, y), 35.174803723282551, 1e-11);

%!test  # the response transform can be set on a fitted ensemble
%! Mdl = RegressionEnsemble (X, y, 'NumLearningCycles', 3, 'Learners', S);
%! Mdl.ResponseTransform = @(z) z + 1;
%! assert_equal (predict (Mdl, X(1,:)), 5.953658536585365, 1e-13);
%! assert_equal (loss (Mdl, X, y), 1.176267398277276, 1e-12);

%!test  # compact keeps the trees and drops the data
%! Mdl = RegressionEnsemble (X, y, 'NumLearningCycles', 3, 'Learners', S);
%! C = compact (Mdl);
%! assert_equal (class (C), 'CompactRegressionEnsemble');
%! assert_equal (predict (C, X), predict (Mdl, X));

## Test input validation
%!error<RegressionEnsemble: too few input arguments.> RegressionEnsemble (X)
%!error<RegressionEnsemble: name-value arguments must be in pairs.> ...
%! RegressionEnsemble (X, y, 'Method')
%!error<RegressionEnsemble: invalid optional paired argument.> ...
%! RegressionEnsemble (X, y, 'Foo', 1)
%!error<RegressionEnsemble: invalid optional paired argument.> ...
%! RegressionEnsemble (X, y, 1, 1)
%!error<RegressionEnsemble: 'Method' must be a character vector.> ...
%! RegressionEnsemble (X, y, 'Method', 1)
%!error<RegressionEnsemble: 'AdaBoostM1' is not a valid ensemble method.> ...
%! RegressionEnsemble (X, y, 'Method', 'AdaBoostM1')
%!error<RegressionEnsemble: 'NumLearningCycles' must be a positive integer.> ...
%! RegressionEnsemble (X, y, 'NumLearningCycles', 0)
%!error<RegressionEnsemble: 'LearnRate' must be a number greater than 0 and no greater than 1.> ...
%! RegressionEnsemble (X, y, 'LearnRate', 0)
%!error<RegressionEnsemble: 'NPrint' must be a positive integer or 'off'.> ...
%! RegressionEnsemble (X, y, 'NPrint', 'on')
%!error<RegressionEnsemble: 'ResponseName' must be a character vector.> ...
%! RegressionEnsemble (X, y, 'ResponseName', 1)
%!error<RegressionEnsemble: 'FResample' must be a number greater than 0 and no greater than 1.> ...
%! RegressionEnsemble (X, y, 'FResample', 2)
%!error<RegressionEnsemble: 'Replace' must be 'on' or 'off'.> ...
%! RegressionEnsemble (X, y, 'Replace', true)
%!error<RegressionEnsemble: 'Resample' must be 'on' or 'off'.> ...
%! RegressionEnsemble (X, y, 'Resample', true)
%!error<RegressionEnsemble: 'NumBins' is not implemented.> ...
%! RegressionEnsemble (X, y, 'NumBins', 10)
%!error<RegressionEnsemble: 'Learners' must be 'tree' or a tree template.> ...
%! RegressionEnsemble (X, y, 'Learners', 'svm')
%!error<RegressionEnsemble: a bagged ensemble is fitted by RegressionBaggedEnsemble.> ...
%! RegressionEnsemble (X, y, 'Method', 'Bag')
%!error<RegressionEnsemble: a resampled ensemble is fitted by RegressionBaggedEnsemble.> ...
%! RegressionEnsemble (X, y, 'Resample', 'on')
%!error<RegressionEnsemble: 'PredictorNames' must be a cell array of character vectors with one element per column of X.> ...
%! RegressionEnsemble (X, y, 'PredictorNames', {'a'})
%!error<RegressionEnsemble: 'Weights' must be nonnegative and not all zero.> ...
%! RegressionEnsemble (X, y, 'Weights', -ones (150, 1))
%!error<RegressionEnsemble.resume: too few input arguments.> ...
%! resume (RegressionEnsemble (X, y, 'NumLearningCycles', 1))
%!error<RegressionEnsemble.resume: NUMLEARNINGCYCLES must be a positive integer.> ...
%! resume (RegressionEnsemble (X, y, 'NumLearningCycles', 1), 0)
%!error<RegressionEnsemble.resume: name-value arguments must be in pairs.> ...
%! resume (RegressionEnsemble (X, y, 'NumLearningCycles', 1), 1, 'NPrint')
%!error<RegressionEnsemble.resume: invalid optional paired argument.> ...
%! resume (RegressionEnsemble (X, y, 'NumLearningCycles', 1), 1, 'Foo', 1)
%!error<RegressionEnsemble.resume: 'NPrint' must be a positive integer or 'off'.> ...
%! resume (RegressionEnsemble (X, y, 'NumLearningCycles', 1), 1, ...
%!         'NPrint', 0)
%!error<RegressionEnsemble.predict: too few input arguments.> ...
%! predict (RegressionEnsemble (X, y, 'NumLearningCycles', 1))
%!error<RegressionEnsemble.loss: too few input arguments.> ...
%! loss (RegressionEnsemble (X, y, 'NumLearningCycles', 1), X)
%!error<RegressionEnsemble.resubLoss: invalid optional paired argument.> ...
%! resubLoss (RegressionEnsemble (X, y, 'NumLearningCycles', 1), 'Foo', 1)

%!test  # MATLAB parity: crossval equals cross-validating at fit time
%! c = cvpartition (150, 'KFold', 4);
%! M = RegressionEnsemble (X, y, 'NumLearningCycles', 2, 'Learners', S);
%! CV = crossval (M, 'CVPartition', c);
%! F = fitrensemble (X, y, 'NumLearningCycles', 2, 'Learners', S, ...
%!                   'CVPartition', c);
%! assert_equal (class (CV), 'RegressionPartitionedEnsemble');
%! assert_equal (kfoldLoss (CV), kfoldLoss (F), 1e-15);

%!error<RegressionEnsemble.crossval: invalid parameter name in optional pair arguments.> ...
%! crossval (RegressionEnsemble (X, y, 'NumLearningCycles', 1), 'Foo', 1)

%!shared X, y, tt, E
%! load fisheriris
%! X = meas(:,2:4);
%! y = meas(:,1);
%! tt = templateTree ('MaxNumSplits', 3);
%! E = fitrensemble (X, y, 'Method', 'LSBoost', 'NumLearningCycles', 20, ...
%!                   'Learners', tt);

%!test  # MATLAB parity: the default penalties of regularize
%! R = regularize (E).Regularization;
%! assert_equal (fieldnames (R), {'Method'; 'TrainedWeights'; 'Lambda'; ...
%!                                'ResubstitutionMSE'; 'CombineWeights'});
%! assert_equal (R.Method, 'Lasso');
%! assert_equal (size (R.TrainedWeights), [20, 10]);
%! assert_equal (R.Lambda(1), 0);
%! assert_equal (R.Lambda([2, 10]), [0.0346843052641099, 34.6843052641099], ...
%!               1e-12);

%!test  # MATLAB parity: lasso weights where R2024a converges
%! R = regularize (E, 'Lambda', [0.01, 0.05, 0.1]).Regularization;
%! assert_equal (R.TrainedWeights(:,[2, 3]), ...
%!               [0.998558425788861, 0.997116851577722; zeros(19, 2)], 1e-9);
%! assert_equal (R.TrainedWeights(:,1), [0.99983699976111; ...
%!               0.688155074556548; 0.205495438447582; 0.0983109617173514; ...
%!               0; 0.314714964487384; zeros(14, 1)], 1e-6);
%! assert_equal (R.ResubstitutionMSE(2:3), ...
%!               [0.141433481267387, 0.141649717399058], 1e-12);
%! assert_equal (R.ResubstitutionMSE(1), 0.0981895827257545, 1e-6);

%!test  # MATLAB parity: observation weights enter the fit and the error
%! Ew = fitrensemble (X, y, 'Method', 'LSBoost', 'NumLearningCycles', 20, ...
%!                    'Learners', tt, 'Weights', (1:150)' / 150);
%! R = regularize (Ew, 'Lambda', 0.1).Regularization;
%! assert_equal (R.TrainedWeights, [0.997413841751991; zeros(19, 1)], 1e-9);
%! assert_equal (R.ResubstitutionMSE, 0.136057521685236, 1e-12);
%! assert_equal (regularize (Ew).Regularization.Lambda(end), ...
%!               38.6673940301219, 1e-12);

%!test  # MATLAB parity: the weights are never negative
%! Elr = fitrensemble (X, y, 'Method', 'LSBoost', 'NumLearningCycles', 20, ...
%!                     'Learners', tt, 'LearnRate', 0.1);
%! R = regularize (Elr, 'Lambda', [0, 0.01]).Regularization;
%! assert_equal (all (R.TrainedWeights(:) >= 0), true);

%!test  # regularize leaves the trained weights as they were
%! Er = regularize (E, 'Lambda', 0.1);
%! assert_equal (Er.TrainedWeights, E.TrainedWeights);

%!test  # MATLAB parity: shrink keeps the weighted trees, largest first
%! C = shrink (regularize (E, 'Lambda', 0.01));
%! assert_equal (class (C), 'CompactRegressionEnsemble');
%! assert_equal (C.CombineWeights, 'WeightedSum');
%! assert_equal (C.NumTrained, 5);
%! assert_equal (C.TrainedWeights, [0.99983699976111; 0.688155074556548; ...
%!               0.314714964487384; 0.205495438447582; ...
%!               0.0983109617173514], 1e-6);

%!test  # MATLAB parity: a weight equal to the threshold is dropped
%! Er = regularize (E, 'Lambda', 0.01);
%! C = shrink (Er);
%! assert_equal (shrink (Er, 'Threshold', C.TrainedWeights(5)).NumTrained, 4);

%!test  # MATLAB parity: shrink regularizes first when given penalties
%! assert_equal (shrink (E, 'Lambda', 0.01).NumTrained, 5);

%!test  # MATLAB parity: an ensemble not regularized keeps its trees
%! C = shrink (E);
%! assert_equal (C.NumTrained, 20);
%! assert_equal (C.TrainedWeights, ones (20, 1));

%!test  # MATLAB parity: a shrunk bagged ensemble sums its weighted trees
%! B = fitrensemble (X, y, 'Method', 'Bag', 'NumLearningCycles', 5);
%! Br = regularize (B, 'Lambda', 0.001);
%! C = shrink (Br);
%! assert_equal (C.CombineWeights, 'WeightedSum');
%! P = zeros (3, 5);
%! for t = 1:5
%!   P(:,t) = predict (B.Trained{t}, X(1:3,:));
%! endfor
%! assert_equal (predict (C, X(1:3,:)), ...
%!               P * Br.Regularization.TrainedWeights, 1e-12);

%!test  # MATLAB parity: cvshrink pools the held-out error over the folds
%! cvp = cvpartition (150, 'KFold', 3);
%! [vals, nlearn] = cvshrink (E, 'CVPartition', cvp, 'Lambda', [0.01, 0.1], ...
%!                            'Threshold', [0, 0.5]);
%! sse = zeros (2, 2);
%! counts = zeros (2, 2);
%! for k = 1:3
%!   tr = training (cvp, k);
%!   te = test (cvp, k);
%!   Ek = regularize (fitrensemble (X(tr,:), y(tr), 'Method', 'LSBoost', ...
%!                                  'NumLearningCycles', 20, ...
%!                                  'Learners', tt), 'Lambda', [0.01, 0.1]);
%!   for a = 1:2
%!     for b = 1:2
%!       thr = [0, 0.5];
%!       Ck = shrink (Ek, 'WeightColumn', a, 'Threshold', thr(b));
%!       sse(a,b) += sum ((y(te) - predict (Ck, X(te,:))) .^ 2) / 150;
%!       counts(a,b) += Ck.NumTrained;
%!     endfor
%!   endfor
%! endfor
%! assert_equal (vals, sse, 1e-12);
%! assert_equal (nlearn, counts / 3);

%!test  # MATLAB parity: cvshrink takes the penalties of a regularized ensemble
%! cvp = cvpartition (150, 'KFold', 3);
%! vals = cvshrink (regularize (E, 'Lambda', [0.01, 0.1]), 'CVPartition', cvp);
%! assert_equal (size (vals), [2, 1]);

%!test  # resume empties the regularization
%! Er = resume (regularize (E, 'Lambda', 0.1), 2);
%! assert_equal (Er.NumTrained, 22);
%! assert_equal (Er.Regularization, []);

%!error<RegressionEnsemble.regularize: 'Lambda' must be a vector of non-negative numbers.> ...
%! regularize (E, 'Lambda', -1)
%!error<RegressionEnsemble.regularize: 'MaxIter' must be a positive integer.> ...
%! regularize (E, 'MaxIter', 0)
%!error<RegressionEnsemble.regularize: 'RelTol' must be a positive number.> ...
%! regularize (E, 'RelTol', 0)
%!error<RegressionEnsemble.regularize: invalid optional paired argument.> ...
%! regularize (E, 'Npass', 3)
%!error<RegressionEnsemble.regularize: invalid optional paired argument.> ...
%! regularize (E, 'Threshold', 0.1)
%!error<RegressionEnsemble.regularize: name-value arguments must be in pairs.> ...
%! regularize (E, 'Lambda')
%!error<RegressionEnsemble.shrink: 'WeightColumn' must be a positive integer.> ...
%! shrink (E, 'WeightColumn', 1.5)
%!error<RegressionEnsemble.shrink: 'WeightColumn' must not exceed the number of values in 'Lambda'.> ...
%! shrink (regularize (E, 'Lambda', 0.1), 'WeightColumn', 2)
%!error<RegressionEnsemble.shrink: 'Threshold' must hold non-negative numbers.> ...
%! shrink (E, 'Threshold', -1)
%!error<RegressionEnsemble.shrink: 'Threshold' must be a scalar.> ...
%! shrink (E, 'Threshold', [0, 1])
%!error<RegressionEnsemble.cvshrink: 'Lambda' must be given for an ensemble that has not been regularized.> ...
%! cvshrink (E)
%!error<RegressionEnsemble.cvshrink: invalid optional paired argument.> ...
%! cvshrink (E, 'Lambda', 0.1, 'Foo', 1)

%!shared X, yr
%! k = (0:119)';
%! c = mod (k, 4) + 1;
%! j = floor (k / 4);
%! x2 = mod (k * 7, 10);
%! X = [c, x2];
%! yb = (c == 1) | (c == 2 & mod (j, 4) != 0) | (c == 3 & mod (j, 4) == 0) ...
%!      | (x2 > 7);
%! yr = [3; 1; 4; 1.5];
%! yr = yr(c) + 0.1 * sin (k) + 0.2 * x2;
%! Xq = [1, 0; 3, 5; 5, 0; NaN, 2; 2.5, 9];

%!test  # bagged regression trees and folds take the categorical predictors
%! Mdl = fitrensemble (X, yr, 'Method', 'Bag', 'NumLearningCycles', 3, ...
%!                     'CategoricalPredictors', 1);
%! assert_equal (Mdl.CategoricalPredictors, 1);
%! assert_equal (Mdl.Trained{1}.CategoricalPredictors, 1);
%! CV = crossval (Mdl, 'KFold', 3);
%! assert_equal (CV.Trained{1}.Trained{1}.CategoricalPredictors, 1);

%!error<RegressionEnsemble: 'CategoricalPredictors' indices must not exceed the number of predictors.> ...
%! RegressionEnsemble (X, yr, 'CategoricalPredictors', 3)

## A table at loss
%!test  # the response is named, left out, or given beside the table
%! load fisheriris
%! X = meas(:,2:3);
%! y = meas(:,1);
%! T = table (X(:,1), X(:,2), 'VariableNames', {'SW', 'PL'});
%! T.SL = y;
%! Mdl = fitrensemble (T, 'SL');
%! a = loss (Mdl, X, y);
%! assert_equal (loss (Mdl, T(:,1:2), y), a);
%! assert_equal (loss (Mdl, T, 'SL'), a);
%! assert_equal (loss (Mdl, T), a);

## Observation weights of class single or double
%!error <RegressionEnsemble: 'Weights' must be a real vector of class single or double.> ...
%! fitrensemble ([1, 2; 3, 4; 5, 6; 7, 8], (1:4)', 'Weights', ...
%!               int8 ([1; 1; 1; 1]))
%!error <RegressionEnsemble: 'Weights' must be a real vector of class single or double.> ...
%! fitrensemble ([1, 2; 3, 4; 5, 6; 7, 8], (1:4)', 'Weights', true (4, 1))
%!test
%! ## Single weights are stored single, summing to one
%! load fisheriris
%! X = meas(:,2:4);
%! y = meas(:,1);
%! w = 1 + (1:150)' / 7;
%! Mdl = fitrensemble (X, y, 'Weights', single (w), 'NumLearningCycles', 10);
%! assert_equal (class (Mdl.W), 'single');
%! assert_equal (sum (double (Mdl.W)), 1, 1e-6);
