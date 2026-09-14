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

classdef ClassificationEnsemble
  ## -*- texinfo -*-
  ## @deftp {statistics} ClassificationEnsemble
  ##
  ## Boosted ensemble of decision trees for classification
  ##
  ## A @code{ClassificationEnsemble} object holds the weak learners a
  ## boosting method grew one after another, each on the observations
  ## reweighted by the errors of those before it, together with the data it
  ## was fitted on.  AdaBoostM1, GentleBoost and LogitBoost fit two classes,
  ## AdaBoostM2 more than two.
  ##
  ## Create one with @code{fitcensemble}.  A bagged ensemble is a
  ## @code{ClassificationBaggedEnsemble}, and @code{compact} returns a
  ## @code{CompactClassificationEnsemble} without the data.
  ##
  ## @seealso{fitcensemble, CompactClassificationEnsemble,
  ## ClassificationBaggedEnsemble}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} X
    ##
    ## Predictor data
    ##
    ## The predictors the ensemble was fitted on, one row per observation, a
    ## row missing a value or a class having been left out.  This property is
    ## read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} Y
    ##
    ## Class labels
    ##
    ## The labels the ensemble was fitted on, in the type they were given in.
    ## This property is read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} RowsUsed
    ##
    ## Rows of the data that were used
    ##
    ## A logical column over the rows as supplied.  This property is
    ## read-only.
    ##
    ## @end deftp
    RowsUsed = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} W
    ##
    ## Observation weights
    ##
    ## The weights given, scaled so that the observations of each class sum
    ## to its prior probability.  A cost matrix is not folded into them.
    ## This property is read-only.
    ##
    ## @end deftp
    W = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} ModelParameters
    ##
    ## Parameters of the fit
    ##
    ## A structure with the fields @qcode{Type}, @qcode{Method},
    ## @qcode{LearnerTemplates}, the tree template the learners were grown
    ## from, @qcode{NLearn}, the number of learning cycles asked for in all,
    ## and @qcode{LearnRate}.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} NumObservations
    ##
    ## Number of observations used
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations = 0;

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## Always empty, binning not being implemented.  This property is
    ## read-only.
    ##
    ## @end deftp
    BinEdges = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} HyperparameterOptimizationResults
    ##
    ## Results of optimizing the hyperparameters
    ##
    ## Always empty, such optimization not being implemented.  This property
    ## is read-only.
    ##
    ## @end deftp
    HyperparameterOptimizationResults = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} CategoricalPredictors
    ##
    ## Indices of categorical predictors
    ##
    ## Always empty, categorical predictors not being implemented.  This
    ## property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName = 'Y';

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the learners saw them
    ##
    ## The same as @code{PredictorNames}.  This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} ClassNames
    ##
    ## Names of the classes
    ##
    ## The classes, in the type of the response and in the order the columns
    ## of the scores take them.  This property is read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} Prior
    ##
    ## Prior probabilities of the classes
    ##
    ## A row vector with one probability per class.  This property is
    ## read-only.
    ##
    ## @end deftp
    Prior = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} Cost
    ##
    ## Misclassification costs
    ##
    ## A square matrix, @code{Cost(i,j)} being the cost of classifying an
    ## observation of class @math{i} as class @math{j}.  A boosting method
    ## starts from observation weights multiplied by the total cost of
    ## misclassifying each observation's class.  This property is read-only.
    ##
    ## @end deftp
    Cost = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} Method
    ##
    ## Ensemble method
    ##
    ## @qcode{'AdaBoostM1'}, @qcode{'AdaBoostM2'}, @qcode{'GentleBoost'},
    ## @qcode{'LogitBoost'} or @qcode{'Bag'}.  This property is read-only.
    ##
    ## @end deftp
    Method = '';

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} LearnerNames
    ##
    ## Names of the weak learners
    ##
    ## Always @code{@{'Tree'@}}.  This property is read-only.
    ##
    ## @end deftp
    LearnerNames = {'Tree'};

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} ReasonForTermination
    ##
    ## Why the fit stopped adding learners
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    ReasonForTermination = '';

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} FitInfo
    ##
    ## Fit information
    ##
    ## A column with one element per learner: the weighted classification
    ## error for AdaBoostM1, the weighted pseudo-loss for AdaBoostM2, and the
    ## weighted mean squared error of the regression tree for GentleBoost and
    ## LogitBoost.  Empty for Bag.  This property is read-only.
    ##
    ## @end deftp
    FitInfo = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} FitInfoDescription
    ##
    ## Description of @code{FitInfo}
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    FitInfoDescription = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} UsePredForLearner
    ##
    ## Which predictors each learner uses
    ##
    ## Always empty, as MATLAB returns it for tree learners.  This property
    ## is read-only.
    ##
    ## @end deftp
    UsePredForLearner = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} NumTrained
    ##
    ## Number of trained weak learners
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumTrained = 0;

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} Trained
    ##
    ## Trained weak learners
    ##
    ## A column cell array of compact trees, as described under
    ## @code{CompactClassificationEnsemble.Trained}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Trained = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} TrainedWeights
    ##
    ## Weights of the trained weak learners
    ##
    ## A column with one weight per learner.  This property is read-only.
    ##
    ## @end deftp
    TrainedWeights = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} CombineWeights
    ##
    ## How the weak learners are combined
    ##
    ## @qcode{'WeightedSum'} for boosting, @qcode{'WeightedAverage'} for
    ## Bag.  This property is read-only.
    ##
    ## @end deftp
    CombineWeights = 'WeightedSum';

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} ScoreTransform
    ##
    ## Transform applied to the predicted scores
    ##
    ## See @code{CompactClassificationEnsemble.ScoreTransform}.
    ##
    ## @end deftp
    ScoreTransform = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    STfun = @(x) x;      # the score transform as a function
    DefaultIndex = 1;    # class given to a row no learner may score
    gY = [];             # class index of each observation
    LearnRate = 1;       # shrinkage of the boosting methods
    TreeArgs = {};       # the options given in the tree template
    State = [];          # what boosting carries from one learner to the next
    Stopped = false;     # whether a learner ended the fit
    BagFResample = 1;    # share of the observations each bag draws
    BagReplace = true;   # whether the bags draw with replacement
    BagInBag = [];       # NxNumTrained logical, the rows each bag drew
  endproperties

  methods (Hidden)

    function this = set.ScoreTransform (this, val)
      try
        [this.STfun, this.ScoreTransform] = parseScoreTransform (val, ...
                                            'ClassificationEnsemble');
      catch
        error (strcat ("ClassificationEnsemble.subsasgn: 'ScoreTransform'", ...
                       " must be a character vector or a", ...
                       " 'function_handle' object."));
      end_try_catch
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
      fprintf ("%+25s: %s\n", 'ClassNames', classNameListing (this.ClassNames));
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumTrained', this.NumTrained);
      fprintf ("%+25s: '%s'\n", 'Method', this.Method);
      fprintf ("%+25s: '%s'\n", 'ReasonForTermination', ...
               this.ReasonForTermination);
      fprintf ("%+25s: [%dx%d double]\n", 'FitInfo', rows (this.FitInfo), ...
               columns (this.FitInfo));
      if (isa (this, 'ClassificationBaggedEnsemble'))
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
    ## @deftypefn  {ClassificationEnsemble} {@var{obj} =} ClassificationEnsemble (@var{X}, @var{Y})
    ## @deftypefnx {ClassificationEnsemble} {@var{obj} =} ClassificationEnsemble (@dots{}, @var{name}, @var{value})
    ##
    ## Fit a boosted ensemble of decision trees.
    ##
    ## @code{fitcensemble} is the documented way in, and its help lists the
    ## options both take.  A bagged ensemble is fitted by
    ## @code{ClassificationBaggedEnsemble}.
    ##
    ## @seealso{fitcensemble, ClassificationBaggedEnsemble}
    ## @end deftypefn
    function this = ClassificationEnsemble (X, Y, varargin)

      if (nargin < 2)
        error ("ClassificationEnsemble: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationEnsemble: name-value arguments must", ...
                       " be in pairs."));
      endif
      bagged = isa (this, 'ClassificationBaggedEnsemble');
      if (bagged)
        caller = 'ClassificationBaggedEnsemble';
      else
        caller = 'ClassificationEnsemble';
      endif

      Method = ''; NLearn = 100; Learners = 'tree'; LearnRate = [];
      NPrint = 0; ClassNames = []; Cost = []; Prior = []; Weights = [];
      PredictorNames = {}; ResponseName = 'Y'; ScoreTransform = 'none';
      FResample = []; Replace = []; Resample = false;

      for i = 1:2:numel (varargin)
        name = varargin{i};
        val = varargin{i+1};
        if (! ischar (name))
          error (strcat ("%s: invalid parameter name in optional pair", ...
                         " arguments."), caller);
        endif
        switch (tolower (name))
          case 'method'
            Method = val;
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
          case 'classnames'
            ClassNames = val;
          case 'cost'
            Cost = val;
          case 'prior'
            Prior = val;
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
          case 'scoretransform'
            ScoreTransform = val;
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
          case {'ratiotosmallest', 'marginprecision', 'robusterrorgoal', ...
                'robustmaxmargin', 'robustmarginsigma', 'numbins', ...
                'optimizehyperparameters', ...
                'hyperparameteroptimizationoptions', 'options'}
            error ("%s: '%s' is not implemented.", caller, name);
          otherwise
            error (strcat ("%s: invalid parameter name in optional pair", ...
                           " arguments."), caller);
        endswitch
      endfor

      ## The learners, by name or as a tree template.
      if (ischar (Learners) && strcmpi (Learners, 'tree'))
        tmpl = templateTree ();
      elseif (isstruct (Learners) && isscalar (Learners)
              && isfield (Learners, 'Method')
              && strcmpi (Learners.Method, 'Tree'))
        tmpl = Learners;
      else
        error ("%s: 'Learners' must be 'tree' or a tree template.", caller);
      endif
      TreeArgs = {};
      for [val, key] = tmpl
        if (! any (strcmp (key, {'Method', 'Type'})))
          TreeArgs(end+1:end+2) = {key, val};
        endif
      endfor

      F = classFrame (X, Y, ClassNames, Prior, Cost, Weights, caller, false);
      K = classCount (F.ClassNames);

      ## The method, by default the one MATLAB chooses for the class count.
      methods2 = {'AdaBoostM1', 'GentleBoost', 'LogitBoost'};
      known = [methods2, {'AdaBoostM2', 'Bag'}];
      later = {'Subspace', 'LPBoost', 'TotalBoost', 'RobustBoost', 'RUSBoost'};
      if (isempty (Method))
        if (K == 2)
          Method = 'LogitBoost';
        else
          Method = 'AdaBoostM2';
        endif
      elseif (! (ischar (Method) && isrow (Method)))
        error ("%s: 'Method' must be a character vector.", caller);
      elseif (any (strcmpi (Method, known)))
        Method = known{strcmpi (Method, known)};
      elseif (any (strcmpi (Method, later)))
        error ("%s: the '%s' method is not implemented.", caller, ...
               later{strcmpi (Method, later)});
      else
        error ("%s: '%s' is not a valid ensemble method.", caller, Method);
      endif
      isbag = strcmp (Method, 'Bag');
      if (isbag && ! bagged)
        error (strcat ("ClassificationEnsemble: a bagged ensemble is", ...
                       " fitted by ClassificationBaggedEnsemble."));
      elseif (! isbag && bagged)
        error (strcat ("ClassificationBaggedEnsemble: 'Method' must be", ...
                       " 'Bag'."));
      endif
      if (any (strcmp (Method, methods2)) && K != 2)
        error ("%s: the '%s' method fits exactly two classes.", caller, Method);
      elseif (strcmp (Method, 'AdaBoostM2') && K < 3)
        error (strcat ("%s: the 'AdaBoostM2' method fits more than two", ...
                       " classes."), caller);
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
          error ("%s: resampling in a boosting method is not implemented.", ...
                 caller);
        endif
        if (isempty (LearnRate))
          LearnRate = 1;
        endif
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
      this.Y = F.Y;
      this.RowsUsed = F.RowsUsed;
      this.NumObservations = F.n;
      this.gY = F.gY(:);
      this.ClassNames = F.ClassNames;
      this.Prior = F.Prior;
      [this.Cost, errmsg] = costMatrix (F.Cost, F.ClassNames);
      if (! isempty (errmsg))
        error ("%s: %s", caller, errmsg);
      endif
      this.W = priorNormalize (F.Weights, this.gY, this.Prior);
      this.W /= sum (this.W);
      [~, this.DefaultIndex] = max (this.Prior);
      this.PredictorNames = PredictorNames(:)';
      this.ExpandedPredictorNames = this.PredictorNames;
      this.ResponseName = ResponseName;
      this.ScoreTransform = ScoreTransform;
      this.Method = Method;
      this.TreeArgs = TreeArgs;
      this.Trained = cell (0, 1);
      this.TrainedWeights = zeros (0, 1);
      this.ModelParameters = struct ('Type', 'classification', ...
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
        this.FitInfoDescription = ClassificationEnsemble.fitInfoText (Method);
      endif

      this = growLearners (this, NLearn, NPrint, caller);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationEnsemble} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Drop the training data from an ensemble.
    ##
    ## @code{@var{CMdl} = compact (@var{obj})} returns a
    ## @code{CompactClassificationEnsemble} holding the learners and what
    ## prediction needs.  It predicts new data identically.
    ##
    ## @seealso{ClassificationEnsemble, CompactClassificationEnsemble}
    ## @end deftypefn
    function CMdl = compact (this)

      CMdl = CompactClassificationEnsemble (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationEnsemble} {@var{B} =} resume (@var{obj}, @var{NumLearningCycles})
    ## @deftypefnx {ClassificationEnsemble} {@var{B} =} resume (@dots{}, 'NPrint', @var{n})
    ##
    ## Grow more weak learners.
    ##
    ## @var{B} is the ensemble with up to @var{NumLearningCycles} further
    ## learners grown as though the fit had asked for them from the start: the
    ## boosting weights carry on from where it stopped.  An ensemble whose fit
    ## ended on a learner that classified the data perfectly grows no more.
    ## @qcode{'NPrint'} is taken as by @code{fitcensemble}.
    ##
    ## MATLAB R2024a restarts the mislabel weights of AdaBoostM2 on a resume,
    ## so its learner weights then differ from those of one longer fit; here
    ## they are the same.
    ##
    ## @seealso{ClassificationEnsemble, fitcensemble}
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
      this = growLearners (this, double (NumLearningCycles), NPrint, caller);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationEnsemble} {@var{label} =} predict (@var{obj}, @var{X})
    ## @deftypefnx {ClassificationEnsemble} {[@var{label}, @var{scores}] =} predict (@dots{})
    ## @deftypefnx {ClassificationEnsemble} {[@dots{}] =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Classify new data with an ensemble.
    ##
    ## Behaves as @code{CompactClassificationEnsemble.predict} and takes the
    ## same Name-Value arguments.
    ##
    ## @seealso{ClassificationEnsemble, CompactClassificationEnsemble.predict}
    ## @end deftypefn
    function [label, scores] = predict (this, X, varargin)

      if (nargin < 2)
        error ("%s.predict: too few input arguments.", class (this));
      endif
      [label, scores] = ensemblePredict (compact (this), X, varargin, ...
                                         [class(this), '.predict']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationEnsemble} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationEnsemble} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of an ensemble.
    ##
    ## Behaves as @code{CompactClassificationEnsemble.loss} and takes the
    ## same Name-Value arguments.
    ##
    ## @seealso{ClassificationEnsemble, CompactClassificationEnsemble.loss}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      if (nargin < 3)
        error ("%s.loss: too few input arguments.", class (this));
      endif
      L = ensembleLoss (compact (this), X, Y, varargin, ...
                        [class(this), '.loss']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationEnsemble} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationEnsemble} {@var{e} =} edge (@dots{}, @var{name}, @var{value})
    ##
    ## Classification edge of an ensemble.
    ##
    ## Behaves as @code{CompactClassificationEnsemble.edge} and takes the
    ## same Name-Value arguments.
    ##
    ## @seealso{ClassificationEnsemble, CompactClassificationEnsemble.edge}
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      if (nargin < 3)
        error ("%s.edge: too few input arguments.", class (this));
      endif
      e = ensembleEdge (compact (this), X, Y, varargin, ...
                        [class(this), '.edge']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationEnsemble} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationEnsemble} {@var{m} =} margin (@dots{}, @var{name}, @var{value})
    ##
    ## Classification margins of an ensemble.
    ##
    ## Behaves as @code{CompactClassificationEnsemble.margin} and takes the
    ## same Name-Value arguments.
    ##
    ## @seealso{ClassificationEnsemble, CompactClassificationEnsemble.margin}
    ## @end deftypefn
    function m = margin (this, X, Y, varargin)

      if (nargin < 3)
        error ("%s.margin: too few input arguments.", class (this));
      endif
      m = ensembleMargin (compact (this), X, Y, varargin, ...
                          [class(this), '.margin']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationEnsemble} {@var{label} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationEnsemble} {[@var{label}, @var{scores}] =} resubPredict (@dots{})
    ## @deftypefnx {ClassificationEnsemble} {[@dots{}] =} resubPredict (@dots{}, @var{name}, @var{value})
    ##
    ## Classify the training data.
    ##
    ## @code{predict} on @code{X}, taking the same Name-Value arguments.
    ##
    ## @seealso{ClassificationEnsemble, ClassificationEnsemble.predict}
    ## @end deftypefn
    function [label, scores] = resubPredict (this, varargin)

      [label, scores] = ensemblePredict (compact (this), this.X, varargin, ...
                                         [class(this), '.resubPredict']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationEnsemble} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {ClassificationEnsemble} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss on the training data.
    ##
    ## @code{loss} on @code{X} and @code{Y}, weighted by @code{W} unless
    ## @qcode{'Weights'} are given.
    ##
    ## @seealso{ClassificationEnsemble, ClassificationEnsemble.loss}
    ## @end deftypefn
    function L = resubLoss (this, varargin)

      L = ensembleLoss (compact (this), this.X, this.Y, ...
                        [{'Weights', this.W}, varargin], ...
                        [class(this), '.resubLoss']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationEnsemble} {@var{e} =} resubEdge (@var{obj})
    ## @deftypefnx {ClassificationEnsemble} {@var{e} =} resubEdge (@dots{}, @var{name}, @var{value})
    ##
    ## Classification edge on the training data.
    ##
    ## @code{edge} on @code{X} and @code{Y}, weighted by @code{W} unless
    ## @qcode{'Weights'} are given.
    ##
    ## @seealso{ClassificationEnsemble, ClassificationEnsemble.edge}
    ## @end deftypefn
    function e = resubEdge (this, varargin)

      e = ensembleEdge (compact (this), this.X, this.Y, ...
                        [{'Weights', this.W}, varargin], ...
                        [class(this), '.resubEdge']);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationEnsemble} {@var{m} =} resubMargin (@var{obj})
    ## @deftypefnx {ClassificationEnsemble} {@var{m} =} resubMargin (@dots{}, @var{name}, @var{value})
    ##
    ## Classification margins of the training data.
    ##
    ## @code{margin} on @code{X} and @code{Y}, taking the same Name-Value
    ## arguments.
    ##
    ## @seealso{ClassificationEnsemble, ClassificationEnsemble.margin}
    ## @end deftypefn
    function m = resubMargin (this, varargin)

      m = ensembleMargin (compact (this), this.X, this.Y, varargin, ...
                          [class(this), '.resubMargin']);

    endfunction

  endmethods

  methods (Access = protected)

    ## Grow up to N learners, carrying the boosting weights in State from one
    ## call to the next, so that a resume continues the fit exactly.
    function this = growLearners (this, N, NPrint, caller)

      n = this.NumObservations;
      K = classCount (this.ClassNames);
      g = this.gY;
      LR = this.LearnRate;
      this.ModelParameters.NLearn += N;
      if (NPrint > 0)
        printf ("Training %s...\n", this.Method);
      endif
      if (isempty (this.State) && ! strcmp (this.Method, 'Bag'))
        csum = sum (this.Cost, 2);
        d0 = this.W .* csum(g);
        d0 /= sum (d0);
        switch (this.Method)
          case 'AdaBoostM2'
            D = repmat (d0 / (K - 1), 1, K);
            D(sub2ind ([n, K], (1:n)', g)) = 0;
            this.State = struct ('D', D);
          case 'LogitBoost'
            this.State = struct ('d0', d0, 'F', zeros (n, 1));
          otherwise
            this.State = struct ('d', d0);
        endswitch
      endif
      ctree = {'PredictorNames', this.PredictorNames, 'MaxNumSplits', 10, ...
               'MinParentSize', 2, 'MinLeafSize', 1, 'Prune', 'off', ...
               'MergeLeaves', 'off'};
      rtree = {'PredictorNames', this.PredictorNames, 'MaxNumSplits', 10, ...
               'MinParentSize', 10, 'MinLeafSize', 1, 'Prune', 'off', ...
               'MergeLeaves', 'off'};
      y = double (g == 1) - double (g == 2);

      for c = 1:N
        if (this.Stopped)
          break;
        endif
        switch (this.Method)

          case 'AdaBoostM1'
            d = this.State.d;
            T = compact (ClassificationTree (this.X, this.Y, 'Weights', d, ...
                                             'ClassNames', this.ClassNames, ...
                                             ctree{:}, this.TreeArgs{:}));
            gh = labelIndices (this.ClassNames, predict (T, this.X));
            h = double (gh == 1) - double (gh == 2);
            e = sum (d .* (h != y));
            if (e > 0.5)
              this.Stopped = true;
              this.ReasonForTermination = sprintf (strcat ("Classification", ...
                " error from the last weak learner is too high: err=%g"), e);
              break;
            endif
            a = LR * log ((1 - e) / max (e, eps)) / 2;
            this = addLearner (this, T, a, e);
            if (e <= 0)
              this = stopPerfect (this, strcat ("Classification error from", ...
                                                " the last weak learner is", ...
                                                " zero."));
              break;
            endif
            d = d .* exp (-a * y .* h);
            this.State.d = d / sum (d);

          case 'AdaBoostM2'
            D = this.State.D;
            tru = sub2ind ([n, K], (1:n)', g);
            T = compact (ClassificationTree (this.X, this.Y, ...
                                             'Weights', sum (D, 2), ...
                                             'ClassNames', this.ClassNames, ...
                                             ctree{:}, this.TreeArgs{:}));
            [~, s] = predict (T, this.X);
            P = zeros (n, K);
            P(:, labelIndices (this.ClassNames, T.ClassNames)) = s;
            hy = P(tru);
            e = sum (sum (D .* (1 - hy + P))) / 2;
            if (e > 0.5)
              this.Stopped = true;
              this.ReasonForTermination = sprintf (strcat ("Pseudo-loss", ...
                " from the last weak learner is too high: err=%g"), e);
              break;
            endif
            a = LR * log ((1 - e) / max (e, eps)) / 2;
            this = addLearner (this, T, a, e);
            if (e <= 0)
              this = stopPerfect (this, strcat ("Pseudo-loss from the last", ...
                                                " weak learner is zero."));
              break;
            endif
            D = D .* exp (-a * (1 + hy - P));
            D(tru) = 0;
            this.State.D = D / sum (D(:));

          case 'GentleBoost'
            d = this.State.d;
            R = compact (RegressionTree (this.X, y, 'Weights', d, rtree{:}, ...
                                         this.TreeArgs{:}));
            h = predict (R, this.X);
            this = addLearner (this, R, LR, sum (d .* (y - h) .^ 2));
            d = d .* exp (-LR * y .* h);
            this.State.d = d / sum (d);

          case 'LogitBoost'
            p = 1 ./ (1 + exp (-this.State.F));
            p = min (max (p, eps), 1 - eps);
            d = this.State.d0 .* p .* (1 - p);
            d /= sum (d);
            z = (double (y > 0) - p) ./ (p .* (1 - p));
            R = compact (RegressionTree (this.X, z, 'Weights', d, rtree{:}, ...
                                         this.TreeArgs{:}));
            h = predict (R, this.X);
            this = addLearner (this, R, LR / 2, sum (d .* (z - h) .^ 2));
            this.State.F += LR / 2 * h;

          case 'Bag'
            m = ceil (this.BagFResample * n);
            if (this.BagReplace)
              cw = [0; cumsum(this.W)];
              cw /= cw(end);
              idx = lookup (cw, rand (m, 1));
            else
              [~, order] = sort (rand (n, 1) .^ (1 ./ this.W), 'descend');
              idx = order(1:m);
            endif
            present = uniqueLabels (labelsFromIndex (this.ClassNames, ...
                                                     unique (g(idx))));
            p = columns (this.X);
            T = compact (ClassificationTree (this.X(idx,:), this.Y(idx,:), ...
                                             'ClassNames', present, ...
                                             'PredictorNames', ...
                                             this.PredictorNames, ...
                                             'MinParentSize', 2, ...
                                             'MinLeafSize', 1, ...
                                             'Prune', 'off', ...
                                             'MergeLeaves', 'off', ...
                                             'NumVariablesToSample', ...
                                             ceil (sqrt (p)), ...
                                             this.TreeArgs{:}));
            inbag = false (n, 1);
            inbag(idx) = true;
            this.BagInBag(:,end+1) = inbag;
            this = addLearner (this, T, 1, []);

        endswitch
        if (NPrint > 0 && mod (this.NumTrained, NPrint) == 0)
          printf ("Grown weak learners: %d\n", this.NumTrained);
        endif
      endfor
      if (! this.Stopped)
        this.ReasonForTermination = strcat ("Terminated normally after", ...
                                            " completing the requested", ...
                                            " number of training cycles.");
      endif

    endfunction

    function this = addLearner (this, mdl, a, info)
      this.Trained{end+1,1} = mdl;
      this.TrainedWeights(end+1,1) = a;
      if (! isempty (info))
        this.FitInfo(end+1,1) = info;
      endif
      this.NumTrained = numel (this.Trained);
    endfunction

    function this = stopPerfect (this, reason)
      this.Stopped = true;
      this.ReasonForTermination = reason;
    endfunction

  endmethods

  methods (Static, Access = private)

    function txt = fitInfoText (Method)
      head = strcat ("Vector of length NumTrained, where NumTrained is the", ...
                     " number of learned weak hypotheses.");
      switch (Method)
        case 'AdaBoostM1'
          tail = strcat ("Element t of this vector is the weighted", ...
                         " classification error from hypothesis t.");
        case 'AdaBoostM2'
          tail = strcat ("Element t of this vector is the weighted", ...
                         " pseudoloss from hypothesis t.");
        otherwise
          tail = strcat ("Element t of this vector is the weighted mean", ...
                         " squared error from regression hypothesis t.");
      endswitch
      txt = {head; tail};
    endfunction

  endmethods

endclassdef

## 'on' or 'off' as a logical scalar, OK false for anything else.
function [tf, ok] = onOff (val)

  ok = ischar (val) && any (strcmpi (val, {'on', 'off'}));
  tf = ok && strcmpi (val, 'on');

endfunction

## Test output
%!shared X2, Y2, S
%! load fisheriris
%! X2 = meas(51:150,:);
%! Y2 = species(51:150);
%! S = templateTree ('MaxNumSplits', 1);

%!test  # MATLAB parity: the properties of a boosted ensemble
%! Mdl = ClassificationEnsemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                               'NumLearningCycles', 2, 'Learners', S);
%! assert_equal (numel (properties (Mdl)), 26);
%! assert_equal (Mdl.NumObservations, 100);
%! assert_equal (Mdl.LearnerNames, {'Tree'});
%! assert_equal (Mdl.UsePredForLearner, []);
%! assert_equal (Mdl.ReasonForTermination, ["Terminated normally after ", ...
%!               "completing the requested number of training cycles."]);
%! assert_equal (size (Mdl.FitInfoDescription), [2, 1]);

%!test  # MATLAB parity: the parameters of a boosted fit
%! Mdl = ClassificationEnsemble (X2, Y2, 'Method', 'LogitBoost', ...
%!                               'NumLearningCycles', 2, 'Learners', S, ...
%!                               'LearnRate', 0.5);
%! mp = Mdl.ModelParameters;
%! assert_equal (mp.Type, 'classification');
%! assert_equal (mp.Method, 'LogitBoost');
%! assert_equal (mp.NLearn, 2);
%! assert_equal (mp.LearnRate, 0.5);
%! assert_equal (mp.LearnerTemplates.MaxNumSplits, 1);

%!test  # MATLAB parity: resuming continues the fit exactly
%! for m = {'AdaBoostM1', 'GentleBoost', 'LogitBoost'}
%!   M5 = ClassificationEnsemble (X2, Y2, 'Method', m{1}, ...
%!                                'NumLearningCycles', 5, 'Learners', S, ...
%!                                'LearnRate', 0.5);
%!   M3 = ClassificationEnsemble (X2, Y2, 'Method', m{1}, ...
%!                                'NumLearningCycles', 3, 'Learners', S, ...
%!                                'LearnRate', 0.5);
%!   R = resume (M3, 2);
%!   assert_equal (R.NumTrained, 5);
%!   assert_equal (R.ModelParameters.NLearn, 5);
%!   assert_equal (R.TrainedWeights, M5.TrainedWeights, 1e-14);
%!   assert_equal (R.FitInfo, M5.FitInfo, 1e-14);
%! endfor

%!test  # resuming AdaBoostM2 continues exactly, where MATLAB restarts
%! load fisheriris
%! M5 = ClassificationEnsemble (meas, species, 'Method', 'AdaBoostM2', ...
%!                              'NumLearningCycles', 5, 'Learners', S);
%! R = resume (ClassificationEnsemble (meas, species, ...
%!                                     'Method', 'AdaBoostM2', ...
%!                                     'NumLearningCycles', 2, ...
%!                                     'Learners', S), 3);
%! assert_equal (R.TrainedWeights, M5.TrainedWeights, 1e-14);

%!test  # an ensemble ended by a perfect learner grows no more
%! load fisheriris
%! Mdl = ClassificationEnsemble (meas, strcmp (species, 'setosa'), ...
%!                               'Method', 'AdaBoostM1', ...
%!                               'NumLearningCycles', 2, 'Learners', S);
%! assert_equal (resume (Mdl, 3).NumTrained, 1);

%!test  # MATLAB parity: resubstitution with the training weights
%! Mdl = ClassificationEnsemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                               'NumLearningCycles', 5, 'Learners', S);
%! assert_equal (resubLoss (Mdl), 0.04, 1e-15);
%! assert_equal (resubLoss (Mdl, 'LossFun', 'exponential'), ...
%!               0.182079166822146, 1e-13);
%! assert_equal (resubEdge (Mdl), 5.115272254726563, 1e-12);
%! m = resubMargin (Mdl);
%! assert_equal (m([1, 21]), [3.544074277136197; -1.958996348947703], 1e-12);
%! assert_equal (resubPredict (Mdl), predict (Mdl, X2));

%!test  # compact keeps the learners and drops the data
%! Mdl = ClassificationEnsemble (X2, Y2, 'Method', 'GentleBoost', ...
%!                               'NumLearningCycles', 3, 'Learners', S);
%! C = compact (Mdl);
%! assert_equal (class (C), 'CompactClassificationEnsemble');
%! [l1, s1] = predict (Mdl, X2);
%! [l2, s2] = predict (C, X2);
%! assert_equal (l2, l1);
%! assert_equal (s2, s1);

%!test  # the score transform can be set on a fitted ensemble
%! Mdl = ClassificationEnsemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                               'NumLearningCycles', 5, 'Learners', S);
%! Mdl.ScoreTransform = 'doublelogit';
%! [~, s] = predict (Mdl, X2(1,:));
%! assert_equal (s, [0.971916134164721, 0.028083865835279], 1e-13);

## Test input validation
%!error<ClassificationEnsemble: too few input arguments.> ...
%! ClassificationEnsemble (X2)
%!error<ClassificationEnsemble: name-value arguments must be in pairs.> ...
%! ClassificationEnsemble (X2, Y2, 'Method')
%!error<ClassificationEnsemble: invalid parameter name in optional pair arguments.> ...
%! ClassificationEnsemble (X2, Y2, 'Foo', 1)
%!error<ClassificationEnsemble: invalid parameter name in optional pair arguments.> ...
%! ClassificationEnsemble (X2, Y2, 1, 1)
%!error<ClassificationEnsemble: 'NumLearningCycles' must be a positive integer.> ...
%! ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 0)
%!error<ClassificationEnsemble: 'LearnRate' must be a number greater than 0 and no greater than 1.> ...
%! ClassificationEnsemble (X2, Y2, 'LearnRate', 2)
%!error<ClassificationEnsemble: 'NPrint' must be a positive integer or 'off'.> ...
%! ClassificationEnsemble (X2, Y2, 'NPrint', 0)
%!error<ClassificationEnsemble: 'ResponseName' must be a character vector.> ...
%! ClassificationEnsemble (X2, Y2, 'ResponseName', 1)
%!error<ClassificationEnsemble: 'FResample' must be a number greater than 0 and no greater than 1.> ...
%! ClassificationEnsemble (X2, Y2, 'FResample', 0)
%!error<ClassificationEnsemble: 'Replace' must be 'on' or 'off'.> ...
%! ClassificationEnsemble (X2, Y2, 'Replace', 1)
%!error<ClassificationEnsemble: 'Resample' must be 'on' or 'off'.> ...
%! ClassificationEnsemble (X2, Y2, 'Resample', 1)
%!error<ClassificationEnsemble: 'CategoricalPredictors' is not implemented.> ...
%! ClassificationEnsemble (X2, Y2, 'CategoricalPredictors', 1)
%!error<ClassificationEnsemble: 'RatioToSmallest' is not implemented.> ...
%! ClassificationEnsemble (X2, Y2, 'RatioToSmallest', 1)
%!error<ClassificationEnsemble: 'Learners' must be 'tree' or a tree template.> ...
%! ClassificationEnsemble (X2, Y2, 'Learners', 'knn')
%!error<ClassificationEnsemble: 'Method' must be a character vector.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 1)
%!error<ClassificationEnsemble: 'Boost' is not a valid ensemble method.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'Boost')
%!error<ClassificationEnsemble: the 'RUSBoost' method is not implemented.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'rusboost')
%!error<ClassificationEnsemble: a bagged ensemble is fitted by ClassificationBaggedEnsemble.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'Bag')
%!error<ClassificationEnsemble: the 'AdaBoostM1' method fits exactly two classes.> ...
%! load fisheriris
%! ClassificationEnsemble (meas, species, 'Method', 'AdaBoostM1')
%!error<ClassificationEnsemble: the 'AdaBoostM2' method fits more than two classes.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'AdaBoostM2')
%!error<ClassificationEnsemble: resampling in a boosting method is not implemented.> ...
%! ClassificationEnsemble (X2, Y2, 'FResample', 0.5)
%!error<ClassificationEnsemble: 'PredictorNames' must be a cell array of character vectors with one element per column of X.> ...
%! ClassificationEnsemble (X2, Y2, 'PredictorNames', {'a'})
%!error<ClassificationEnsemble.resume: too few input arguments.> ...
%! resume (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1))
%!error<ClassificationEnsemble.resume: NUMLEARNINGCYCLES must be a positive integer.> ...
%! resume (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1), 0)
%!error<ClassificationEnsemble.resume: name-value arguments must be in pairs.> ...
%! resume (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1), 1, 'NPrint')
%!error<ClassificationEnsemble.resume: invalid parameter name in optional pair arguments.> ...
%! resume (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1), 1, 'Foo', 1)
%!error<ClassificationEnsemble.resume: 'NPrint' must be a positive integer or 'off'.> ...
%! resume (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1), 1, ...
%!         'NPrint', 0)
%!error<ClassificationEnsemble.predict: too few input arguments.> ...
%! predict (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1))
%!error<ClassificationEnsemble.loss: too few input arguments.> ...
%! loss (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1), X2)
%!error<ClassificationEnsemble.edge: too few input arguments.> ...
%! edge (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1), X2)
%!error<ClassificationEnsemble.margin: too few input arguments.> ...
%! margin (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1), X2)
%!error<ClassificationEnsemble.resubLoss: invalid parameter name in optional pair arguments.> ...
%! resubLoss (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1), 'Foo', 1)
%!error<ClassificationEnsemble.subsasgn: 'ScoreTransform' must be a character vector or a 'function_handle' object.> ...
%! Mdl = ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1);
%! Mdl.ScoreTransform = 1;
