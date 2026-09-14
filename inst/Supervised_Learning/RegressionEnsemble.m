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

classdef RegressionEnsemble
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
    ## Always empty, @code{regularize} not being implemented.  This property
    ## is read-only.
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
    ## The weights given, normalized to sum to one.  This property is
    ## read-only.
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
    ## Always empty, categorical predictors not being implemented.  This
    ## property is read-only.
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
  endproperties

  methods (Hidden)

    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = ...
        parseResponseTransform (val, 'RegressionEnsemble');
    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
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
      if (mod (numel (varargin), 2) != 0)
        error ("RegressionEnsemble: name-value arguments must be in pairs.");
      endif
      bagged = isa (this, 'RegressionBaggedEnsemble');
      if (bagged)
        caller = 'RegressionBaggedEnsemble';
      else
        caller = 'RegressionEnsemble';
      endif

      Method = 'LSBoost'; NLearn = 100; Learners = 'tree'; LearnRate = [];
      NPrint = 0; Weights = []; PredictorNames = {}; ResponseName = 'Y';
      ResponseTransform = 'none'; FResample = []; Replace = [];
      Resample = false;

      for i = 1:2:numel (varargin)
        name = varargin{i};
        val = varargin{i+1};
        if (! ischar (name))
          error (strcat ("%s: invalid parameter name in optional pair", ...
                         " arguments."), caller);
        endif
        switch (tolower (name))
          case 'method'
            if (! (ischar (val) && isrow (val)))
              error ("%s: 'Method' must be a character vector.", caller);
            elseif (strcmpi (val, 'LSBoost'))
              Method = 'LSBoost';
            elseif (strcmpi (val, 'Bag'))
              Method = 'Bag';
            else
              error ("%s: '%s' is not a valid ensemble method.", caller, val);
            endif
          case 'numlearningcycles'
            if (! (isnumeric (val) && isscalar (val) && isreal (val)
                   && val >= 1 && val == fix (val)))
              error (strcat ("%s: 'NumLearningCycles' must be a positive", ...
                             " integer."), caller);
            endif
            NLearn = double (val);
          case 'learners'
            Learners = val;
          case 'learnrate'
            if (! (isnumeric (val) && isscalar (val) && isreal (val)
                   && val > 0 && val <= 1))
              error (strcat ("%s: 'LearnRate' must be a number greater", ...
                             " than 0 and no greater than 1."), caller);
            endif
            LearnRate = double (val);
          case 'nprint'
            if (ischar (val) && strcmpi (val, 'off'))
              NPrint = 0;
            elseif (isnumeric (val) && isscalar (val) && isreal (val)
                    && val >= 1 && val == fix (val))
              NPrint = double (val);
            else
              error ("%s: 'NPrint' must be a positive integer or 'off'.", ...
                     caller);
            endif
          case 'weights'
            Weights = val;
          case 'predictornames'
            PredictorNames = val;
          case 'responsename'
            if (! (ischar (val) && isrow (val)))
              error ("%s: 'ResponseName' must be a character vector.", ...
                     caller);
            endif
            ResponseName = val;
          case 'responsetransform'
            ResponseTransform = val;
          case 'fresample'
            if (! (isnumeric (val) && isscalar (val) && isreal (val)
                   && val > 0 && val <= 1))
              error (strcat ("%s: 'FResample' must be a number greater", ...
                             " than 0 and no greater than 1."), caller);
            endif
            FResample = double (val);
          case 'replace'
            [Replace, ok] = onOff (val);
            if (! ok)
              error ("%s: 'Replace' must be 'on' or 'off'.", caller);
            endif
          case 'resample'
            [Resample, ok] = onOff (val);
            if (! ok)
              error ("%s: 'Resample' must be 'on' or 'off'.", caller);
            endif
          case 'categoricalpredictors'
            if (! isempty (val))
              error ("%s: 'CategoricalPredictors' is not implemented.", ...
                     caller);
            endif
          case {'numbins', 'optimizehyperparameters', ...
                'hyperparameteroptimizationoptions', 'options'}
            error ("%s: '%s' is not implemented.", caller, name);
          otherwise
            error (strcat ("%s: invalid parameter name in optional pair", ...
                           " arguments."), caller);
        endswitch
      endfor

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
      if (isbag && ! bagged)
        error (strcat ("RegressionEnsemble: a bagged ensemble is fitted by", ...
                       " RegressionBaggedEnsemble."));
      elseif (! isbag && bagged)
        error ("RegressionBaggedEnsemble: 'Method' must be 'Bag'.");
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
        if (Resample || ! isempty (FResample) || ! isempty (Replace))
          error ("%s: resampling in LSBoost is not implemented.", caller);
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

      this.X = F.X;
      this.Y = double (F.Y);
      this.RowsUsed = F.RowsUsed;
      this.W = F.W;
      this.NumObservations = F.n;
      this.PredictorNames = PredictorNames(:)';
      this.ExpandedPredictorNames = this.PredictorNames;
      this.ResponseName = ResponseName;
      this.ResponseTransform = ResponseTransform;
      this.Method = Method;
      this.TreeArgs = TreeArgs;
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
    ## @qcode{'NPrint'} is taken as by @code{fitrensemble}.
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
      NPrint = 0;
      for i = 1:2:numel (varargin)
        if (! (ischar (varargin{i}) && strcmpi (varargin{i}, 'NPrint')))
          error (strcat ("%s: invalid parameter name in optional pair", ...
                         " arguments."), caller);
        endif
        val = varargin{i+1};
        if (ischar (val) && strcmpi (val, 'off'))
          NPrint = 0;
        elseif (isnumeric (val) && isscalar (val) && isreal (val)
                && val >= 1 && val == fix (val))
          NPrint = double (val);
        else
          error ("%s: 'NPrint' must be a positive integer or 'off'.", caller);
        endif
      endfor
      this = growLearners (this, double (NumLearningCycles), NPrint);

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
    ## @seealso{RegressionEnsemble, CompactRegressionEnsemble.predict}
    ## @end deftypefn
    function yfit = predict (this, X, varargin)

      if (nargin < 2)
        error ("%s.predict: too few input arguments.", class (this));
      endif
      yfit = ensemblePredict (compact (this), X, varargin, ...
                              [class(this), '.predict']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionEnsemble} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionEnsemble} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss of an ensemble.
    ##
    ## Behaves as @code{CompactRegressionEnsemble.loss} and takes the same
    ## Name-Value arguments.
    ##
    ## @seealso{RegressionEnsemble, CompactRegressionEnsemble.loss}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      if (nargin < 3)
        error ("%s.loss: too few input arguments.", class (this));
      endif
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
          T = compact (RegressionTree (this.X, r, 'Weights', this.W, ...
                                       'PredictorNames', ...
                                       this.PredictorNames, ...
                                       'MaxNumSplits', 10, ...
                                       'MinParentSize', 10, ...
                                       'MinLeafSize', 5, 'Prune', 'off', ...
                                       'MergeLeaves', 'off', ...
                                       this.TreeArgs{:}));
          h = predict (T, this.X);
          this.FitInfo(end+1,1) = sum (this.W .* (r - h) .^ 2);
          this.F += this.LearnRate * h;
          this = addLearner (this, T, this.LearnRate);
        else
          m = ceil (this.BagFResample * n);
          if (this.BagReplace)
            cw = [0; cumsum(this.W)];
            cw /= cw(end);
            idx = lookup (cw, rand (m, 1));
          else
            [~, order] = sort (rand (n, 1) .^ (1 ./ this.W), 'descend');
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
%!error<RegressionEnsemble: invalid parameter name in optional pair arguments.> ...
%! RegressionEnsemble (X, y, 'Foo', 1)
%!error<RegressionEnsemble: invalid parameter name in optional pair arguments.> ...
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
%!error<RegressionEnsemble: 'CategoricalPredictors' is not implemented.> ...
%! RegressionEnsemble (X, y, 'CategoricalPredictors', 1)
%!error<RegressionEnsemble: 'NumBins' is not implemented.> ...
%! RegressionEnsemble (X, y, 'NumBins', 10)
%!error<RegressionEnsemble: 'Learners' must be 'tree' or a tree template.> ...
%! RegressionEnsemble (X, y, 'Learners', 'svm')
%!error<RegressionEnsemble: a bagged ensemble is fitted by RegressionBaggedEnsemble.> ...
%! RegressionEnsemble (X, y, 'Method', 'Bag')
%!error<RegressionEnsemble: resampling in LSBoost is not implemented.> ...
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
%!error<RegressionEnsemble.resume: invalid parameter name in optional pair arguments.> ...
%! resume (RegressionEnsemble (X, y, 'NumLearningCycles', 1), 1, 'Foo', 1)
%!error<RegressionEnsemble.resume: 'NPrint' must be a positive integer or 'off'.> ...
%! resume (RegressionEnsemble (X, y, 'NumLearningCycles', 1), 1, ...
%!         'NPrint', 0)
%!error<RegressionEnsemble.predict: too few input arguments.> ...
%! predict (RegressionEnsemble (X, y, 'NumLearningCycles', 1))
%!error<RegressionEnsemble.loss: too few input arguments.> ...
%! loss (RegressionEnsemble (X, y, 'NumLearningCycles', 1), X)
%!error<RegressionEnsemble.resubLoss: invalid parameter name in optional pair arguments.> ...
%! resubLoss (RegressionEnsemble (X, y, 'NumLearningCycles', 1), 'Foo', 1)
