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

classdef ClassificationEnsemble < PredictiveModel
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
    ## The predictors every tree treats as categorical, empty when none
    ## is.  This property is read-only.
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
    ## The classes, in the type of the response, sorted or in the order given
    ## by the @qcode{'ClassNames'} option; the columns of the scores take them
    ## in that order.  This property is read-only.
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
    ## @qcode{'LogitBoost'}, @qcode{'RUSBoost'}, @qcode{'LPBoost'},
    ## @qcode{'TotalBoost'}, @qcode{'Bag'} or @qcode{'Subspace'}.  This
    ## property is read-only.
    ##
    ## @end deftp
    Method = '';

    ## -*- texinfo -*-
    ## @deftp {ClassificationEnsemble} {property} LearnerNames
    ##
    ## Names of the weak learners
    ##
    ## @code{@{'Tree'@}}, or for Subspace @code{@{'KNN'@}} or
    ## @code{@{'Discriminant'@}}.  This property is read-only.
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
    ## A column with one element per learner: the weighted classification error
    ## for AdaBoostM1, the weighted pseudo-loss for AdaBoostM2 and RUSBoost, and
    ## the weighted mean squared error of the regression tree for GentleBoost
    ## and LogitBoost.  For LPBoost and TotalBoost a matrix with a row per
    ## learner: its margin on each observation, and its edge last.  Empty for
    ## Bag.  This property is read-only.
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
    ## For the @qcode{'Subspace'} method, a logical matrix with one row per
    ## predictor and one column per learner.  Empty for tree learners, as
    ## MATLAB returns it.  This property is read-only.
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
    NPredToSample = [];  # predictors each Subspace learner is fitted on
    AllCombinations = false;  # whether Subspace takes every combination
    RatioToSmallest = [];  # RUSBoost's class sample sizes over the smallest
    Resampling = false;  # whether a boosting method resamples its rows
    MarginPrecision = 0.01;  # LPBoost's and TotalBoost's edge tolerance
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

      Method = ''; NLearn = 100; Learners = []; LearnRate = [];
      NPred = []; AllCombinations = false;
      NPrint = 0; ClassNames = []; Cost = []; Prior = []; Weights = [];
      PredictorNames = {}; ResponseName = 'Y'; ScoreTransform = 'none';
      CatPreds = [];
      FResample = []; Replace = []; Resample = false; Ratio = [];
      MarginPrecision = [];

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
            if (ischar (val) && strcmpi (val, 'AllPredictorCombinations'))
              AllCombinations = true;
            elseif (isnumeric (val) && isscalar (val) && isreal (val)
                    && val >= 1 && val == fix (val))
              NLearn = double (val);
            else
              error (strcat ("%s: 'NumLearningCycles' must be a positive", ...
                             " integer or 'AllPredictorCombinations'."), ...
                     caller);
            endif
          case 'npredtosample'
            if (! (isnumeric (val) && isscalar (val) && isreal (val)
                   && val >= 1 && val == fix (val)))
              error ("%s: 'NPredToSample' must be a positive integer.", ...
                     caller);
            endif
            NPred = double (val);
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
            CatPreds = val;
          case 'ratiotosmallest'
            if (! (isnumeric (val) && isvector (val) && isreal (val)
                   && all (isfinite (val)) && all (val >= 0) && any (val > 0)))
              error (strcat ("%s: 'RatioToSmallest' must be a vector of", ...
                             " nonnegative numbers with a positive", ...
                             " element."), caller);
            endif
            Ratio = double (val(:)');
          case 'marginprecision'
            if (! (isnumeric (val) && isscalar (val) && isreal (val)
                   && val >= 0 && val <= 1))
              error ("%s: 'MarginPrecision' must be a number from 0 to 1.", ...
                     caller);
            endif
            MarginPrecision = double (val);
          case {'robusterrorgoal', ...
                'robustmaxmargin', 'robustmarginsigma', 'numbins', ...
                'optimizehyperparameters', ...
                'hyperparameteroptimizationoptions', 'options'}
            error ("%s: '%s' is not implemented.", caller, name);
          otherwise
            error (strcat ("%s: invalid parameter name in optional pair", ...
                           " arguments."), caller);
        endswitch
      endfor

      F = classFrame (X, Y, ClassNames, Prior, Cost, Weights, caller, false);
      K = classCount (F.ClassNames);

      ## The method, by default the one MATLAB chooses for the class count.
      methods2 = {'AdaBoostM1', 'GentleBoost', 'LogitBoost'};
      known = [methods2, {'AdaBoostM2', 'RUSBoost', 'LPBoost', ...
                          'TotalBoost', 'Bag', 'Subspace'}];
      later = {'RobustBoost'};
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
      issub = strcmp (Method, 'Subspace');
      resampled = Resample || ! isempty (FResample) || ! isempty (Replace);
      ismargin = any (strcmp (Method, {'LPBoost', 'TotalBoost'}));
      if (resampled && strcmp (Method, 'RUSBoost'))
        error (strcat ("%s: the 'RUSBoost' method cannot resample the", ...
                       " observations."), caller);
      elseif (resampled && ismargin)
        error (strcat ("%s: resampling with the '%s' method is not", ...
                       " implemented."), caller, Method);
      elseif (resampled && issub)
        error (strcat ("%s: the 'Subspace' method cannot resample the", ...
                       " observations."), caller);
      endif
      if (isbag && ! bagged)
        error (strcat ("ClassificationEnsemble: a bagged ensemble is", ...
                       " fitted by ClassificationBaggedEnsemble."));
      elseif (! isbag && bagged && ! resampled)
        error (strcat ("ClassificationBaggedEnsemble: 'Method' must be", ...
                       " 'Bag' unless the ensemble resamples."));
      elseif (! isbag && ! bagged && resampled)
        error (strcat ("ClassificationEnsemble: a resampled ensemble is", ...
                       " fitted by ClassificationBaggedEnsemble."));
      endif
      if (any (strcmp (Method, methods2)) && K != 2)
        error ("%s: the '%s' method fits exactly two classes.", caller, Method);
      elseif (strcmp (Method, 'AdaBoostM2') && K < 3)
        error (strcat ("%s: the 'AdaBoostM2' method fits more than two", ...
                       " classes."), caller);
      endif

      ## The learners, by name or as a template: trees for every method but
      ## Subspace, which takes nearest neighbours or discriminants instead.
      if (isempty (Learners))
        Learners = ifelse_learner (issub);
      endif
      if (ischar (Learners) && isrow (Learners))
        names = {'tree', 'knn', 'discriminant'};
        makers = {@templateTree, @templateKNN, @templateDiscriminant};
        k = find (strcmpi (Learners, names));
        if (isempty (k))
          tmpl = [];
        else
          tmpl = makers{k} ();
        endif
      elseif (isstruct (Learners) && isscalar (Learners)
              && isfield (Learners, 'Method'))
        tmpl = Learners;
      else
        tmpl = [];
      endif
      if (issub)
        if (isempty (tmpl) || ! any (strcmpi (tmpl.Method, ...
                                             {'KNN', 'Discriminant', 'Tree'})))
          error (strcat ("%s: 'Learners' must be 'knn', 'discriminant' or", ...
                         " a template of either for the 'Subspace'", ...
                         " method."), caller);
        elseif (strcmpi (tmpl.Method, 'Tree'))
          error (strcat ("%s: trees cannot be the learners of the", ...
                         " 'Subspace' method."), caller);
        endif
      elseif (isempty (tmpl) || ! strcmpi (tmpl.Method, 'Tree'))
        error ("%s: 'Learners' must be 'tree' or a tree template.", caller);
      endif
      TreeArgs = {};
      for [val, key] = tmpl
        if (! any (strcmp (key, {'Method', 'Type'})))
          TreeArgs(end+1:end+2) = {key, val};
        endif
      endfor
      if (! issub && ! isempty (NPred))
        error ("%s: 'NPredToSample' applies only to the 'Subspace' method.", ...
               caller);
      endif
      if (! issub && AllCombinations)
        error (strcat ("%s: 'AllPredictorCombinations' applies only to the", ...
                       " 'Subspace' method."), caller);
      endif
      if (strcmp (Method, 'RUSBoost'))
        if (isempty (Ratio))
          Ratio = ones (1, K);
        elseif (isscalar (Ratio))
          Ratio = repmat (Ratio, 1, K);
        elseif (numel (Ratio) != K)
          error ("%s: 'RatioToSmallest' must have one element per class.", ...
                 caller);
        endif
      elseif (! isempty (Ratio))
        error (strcat ("%s: 'RatioToSmallest' applies only to the", ...
                       " 'RUSBoost' method."), caller);
      endif
      if (ismargin)
        if (! isempty (LearnRate))
          error ("%s: 'LearnRate' cannot be used with the '%s' method.", ...
                 caller, Method);
        endif
        if (exist ('__glpk__') == 0)
          error (strcat ("%s: the '%s' method needs GLPK, which this", ...
                         " Octave was built without."), caller, Method);
        endif
        if (isempty (MarginPrecision))
          MarginPrecision = 0.01;
        endif
      elseif (! isempty (MarginPrecision))
        error (strcat ("%s: 'MarginPrecision' applies only to the", ...
                       " 'LPBoost' and 'TotalBoost' methods."), caller);
      endif

      if (issub)
        if (! isempty (LearnRate))
          error (strcat ("%s: 'LearnRate' cannot be used with the", ...
                         " 'Subspace' method."), caller);
        endif
      elseif (isbag)
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

      if (isempty (PredictorNames))
        PredictorNames = arrayfun (@(k) sprintf ('x%d', k), 1:F.p, ...
                                   'UniformOutput', false);
      elseif (! (iscellstr (PredictorNames) && numel (PredictorNames) == F.p))
        error (strcat ("%s: 'PredictorNames' must be a cell array of", ...
                       " character vectors with one element per column", ...
                       " of X."), caller);
      endif
      ## Categorical predictors reach every tree, and no other learner
      [Cod, errmsg] = dummyCoding (F.X, CatPreds, PredictorNames);
      if (! isempty (errmsg))
        error ("%s: %s", caller, errmsg);
      endif
      CatIdx = [];
      if (! isempty (Cod.Index))
        if (issub)
          error (strcat ("%s: 'CategoricalPredictors' cannot be used with", ...
                         " the 'Subspace' method."), caller);
        endif
        CatIdx = Cod.Index;
        TreeArgs = [{'CategoricalPredictors', CatIdx}, TreeArgs];
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
      if (issub)
        ## The nearest neighbour and discriminant learners take no weights,
        ## so weights that are not uniform cannot reach them.
        if (max (F.Weights) - min (F.Weights) > 1e-12 * max (F.Weights))
          error (strcat ("%s: 'Weights' that are not uniform cannot be", ...
                         " used with the 'Subspace' method."), caller);
        endif
        if (isempty (NPred))
          NPred = 1;
        endif
        if (NPred >= F.p)
          error (strcat ("%s: 'NPredToSample' must be less than the number", ...
                         " of predictors."), caller);
        endif
        if (AllCombinations)
          NLearn = nchoosek (F.p, NPred);
        endif
      endif
      this.PredictorNames = PredictorNames(:)';
      this.ExpandedPredictorNames = this.PredictorNames;
      this.ResponseName = ResponseName;
      this.ScoreTransform = ScoreTransform;
      this.Method = Method;
      this.TreeArgs = TreeArgs;
      this.CategoricalPredictors = CatIdx;
      this.Trained = cell (0, 1);
      this.TrainedWeights = zeros (0, 1);
      this.ModelParameters = struct ('Type', 'classification', ...
                                     'Method', Method, ...
                                     'LearnerTemplates', tmpl, ...
                                     'NLearn', 0);
      if (issub)
        this.CombineWeights = 'WeightedAverage';
        this.FitInfo = [];
        this.FitInfoDescription = 'None';
        this.LearnerNames = {tmpl.Method};
        this.NPredToSample = NPred;
        this.AllCombinations = AllCombinations;
        this.UsePredForLearner = false (F.p, 0);
      elseif (isbag)
        this.CombineWeights = 'WeightedAverage';
        this.FitInfo = [];
        this.FitInfoDescription = 'None';
        this.BagFResample = FResample;
        this.BagReplace = Replace;
        this.BagInBag = false (F.n, 0);
      else
        this.LearnRate = LearnRate;
        this.RatioToSmallest = Ratio;
        this.Resampling = resampled;
        if (resampled)
          this.BagFResample = FResample;
          this.BagReplace = Replace;
          this.BagInBag = false (F.n, 0);
        endif
        if (ismargin)
          this.MarginPrecision = MarginPrecision;
          this.ModelParameters.MarginPrecision = MarginPrecision;
          this.FitInfo = zeros (0, F.n + 1);
        else
          this.ModelParameters.LearnRate = LearnRate;
          this.FitInfo = zeros (0, 1);
        endif
        this.FitInfoDescription = ClassificationEnsemble.fitInfoText (Method);
      endif

      this = growLearners (this, NLearn, NPrint, caller);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationEnsemble} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationEnsemble} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate an ensemble.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} refits the ensemble on the
    ## training part of each of ten folds, stratified by class, and returns a
    ## @code{ClassificationPartitionedEnsemble}.  One of @qcode{'KFold'}, an
    ## integer greater than 1, @qcode{'Holdout'}, a number between 0 and 1,
    ## @qcode{'Leaveout'}, @qcode{'on'} for one fold per observation, or
    ## @qcode{'CVPartition'}, a @code{cvpartition} object, may choose the
    ## partition instead.
    ##
    ## @seealso{ClassificationEnsemble, ClassificationPartitionedEnsemble, cvpartition}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)

      [P, errmsg] = ensemblePartition (varargin, this.Y, ...
                                       this.NumObservations, ...
                                       true);
      if (! isempty (errmsg))
        error ("%s.crossval: %s", class (this), errmsg);
      endif
      if (isempty (P))
        error (strcat ("%s.crossval: no partition was asked for; set", ...
                       " 'CrossVal' to 'on' or give 'KFold'."), class (this));
      endif
      CVMdl = ClassificationPartitionedEnsemble (this, P);

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
      if (this.AllCombinations)
        error (strcat ("%s: an ensemble of all predictor combinations", ...
                       " cannot grow further."), caller);
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
    ## @deftypefn  {ClassificationEnsemble} {@var{imp} =} predictorImportance (@var{obj})
    ## @deftypefnx {ClassificationEnsemble} {[@var{imp}, @var{ma}] =} predictorImportance (@var{obj})
    ##
    ## Estimate the importance of each predictor.
    ##
    ## Behaves as @code{CompactClassificationEnsemble.predictorImportance}.
    ##
    ## @seealso{ClassificationEnsemble,
    ## CompactClassificationEnsemble.predictorImportance}
    ## @end deftypefn
    function [imp, ma] = predictorImportance (this)

      [imp, ma] = ensembleImportance (compact (this), ...
                                      [class(this), '.predictorImportance']);

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
      if (isempty (this.State) && ! any (strcmp (this.Method, ...
                                                 {'Bag', 'Subspace'})))
        csum = sum (this.Cost, 2);
        d0 = this.W .* csum(g);
        d0 /= sum (d0);
        switch (this.Method)
          case 'AdaBoostM2'
            if (this.Resampling)
              this.State = struct ('d', d0);
            else
              D = repmat (d0 / (K - 1), 1, K);
              D(sub2ind ([n, K], (1:n)', g)) = 0;
              this.State = struct ('D', D);
            endif
          case 'LogitBoost'
            this.State = struct ('d0', d0, 'F', zeros (n, 1), 'w', d0);
          case 'LPBoost'
            this.State = struct ('d', d0, 'ghat', Inf, 'U', zeros (0, n));
          case 'TotalBoost'
            this.State = struct ('d0', d0, 'd', d0, 'ghat', Inf, ...
                                 'U', zeros (0, n));
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
      rs = this.Resampling;
      m = ceil (this.BagFResample * n);

      for c = 1:N
        if (this.Stopped)
          break;
        endif
        switch (this.Method)

          case 'AdaBoostM1'
            d = this.State.d;
            if (rs)
              [idx, sw, cnt] = boostSample (d, m, this.BagReplace);
              present = labelsFromIndex (this.ClassNames, unique (g(idx)));
              T = compact (ClassificationTree (this.X(idx,:), ...
                                               this.Y(idx,:), ...
                                               'Weights', sw, ...
                                               'ClassNames', present, ...
                                               ctree{:}, this.TreeArgs{:}));
            else
              T = compact (ClassificationTree (this.X, this.Y, ...
                                               'Weights', d, ...
                                               'ClassNames', ...
                                               this.ClassNames, ...
                                               ctree{:}, this.TreeArgs{:}));
            endif
            gh = labelIndices (this.ClassNames, predict (T, this.X));
            h = double (gh == 1) - double (gh == 2);
            if (rs)
              e = sum (sw .* (h(idx) != y(idx)));
            else
              e = sum (d .* (h != y));
            endif
            if (e > 0.5)
              this.Stopped = true;
              this.ReasonForTermination = sprintf (strcat ("Classification", ...
                " error from the last weak learner is too high: err=%g"), e);
              break;
            endif
            a = LR * log ((1 - e) / max (e, eps)) / 2;
            this = addLearner (this, T, a, e);
            if (rs)
              this.BagInBag(:,end+1) = cnt > 0;
            endif
            if (e <= 0)
              this = stopPerfect (this, strcat ("Classification error from", ...
                                                " the last weak learner is", ...
                                                " zero."));
              break;
            endif
            if (rs)
              this.State.d = boostRescale (d, d .* exp (-a * y .* h), cnt);
            else
              d = d .* exp (-a * y .* h);
              this.State.d = d / sum (d);
            endif

          case 'AdaBoostM2'
            tru = sub2ind ([n, K], (1:n)', g);
            if (rs)
              ## A resampled AdaBoostM2 keeps one weight per observation,
              ## as RUSBoost does, and as MATLAB R2024a does.
              d = this.State.d;
              [idx, sw, cnt] = boostSample (d, m, this.BagReplace);
              present = labelsFromIndex (this.ClassNames, unique (g(idx)));
              T = compact (ClassificationTree (this.X(idx,:), ...
                                               this.Y(idx,:), ...
                                               'Weights', sw, ...
                                               'ClassNames', present, ...
                                               ctree{:}, this.TreeArgs{:}));
            else
              D = this.State.D;
              T = compact (ClassificationTree (this.X, this.Y, ...
                                               'Weights', sum (D, 2), ...
                                               'ClassNames', ...
                                               this.ClassNames, ...
                                               ctree{:}, this.TreeArgs{:}));
            endif
            [~, s] = predict (T, this.X);
            P = zeros (n, K);
            P(:, labelIndices (this.ClassNames, T.ClassNames)) = s;
            hy = P(tru);
            if (rs)
              L = 1 - hy + P;
              L(tru) = 0;
              e = sum (sw .* sum (L(idx,:), 2)) / (2 * (K - 1));
            else
              e = sum (sum (D .* (1 - hy + P))) / 2;
            endif
            if (e > 0.5)
              this.Stopped = true;
              this.ReasonForTermination = sprintf (strcat ("Pseudo-loss", ...
                " from the last weak learner is too high: err=%g"), e);
              break;
            endif
            a = LR * log ((1 - e) / max (e, eps)) / 2;
            this = addLearner (this, T, a, e);
            if (rs)
              this.BagInBag(:,end+1) = cnt > 0;
            endif
            if (e <= 0)
              this = stopPerfect (this, strcat ("Pseudo-loss from the last", ...
                                                " weak learner is zero."));
              break;
            endif
            if (rs)
              E = exp (-a * (1 + hy - P));
              E(tru) = 0;
              this.State.d = boostRescale (d, d .* sum (E, 2) / (K - 1), cnt);
            else
              D = D .* exp (-a * (1 + hy - P));
              D(tru) = 0;
              this.State.D = D / sum (D(:));
            endif

          case 'RUSBoost'
            d = this.State.d;
            idx = rusSample (g, d, this.RatioToSmallest);
            present = labelsFromIndex (this.ClassNames, unique (g(idx)));
            T = compact (ClassificationTree (this.X(idx,:), this.Y(idx,:), ...
                                             'ClassNames', present, ...
                                             ctree{:}, this.TreeArgs{:}));
            [~, s] = predict (T, this.X);
            P = zeros (n, K);
            P(:, labelIndices (this.ClassNames, T.ClassNames)) = s;
            tru = sub2ind ([n, K], (1:n)', g);
            hy = P(tru);
            L = 1 - hy + P;
            L(tru) = 0;
            e = sum (d .* sum (L, 2)) / (2 * (K - 1));
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
            E = exp (-a * (1 + hy - P));
            E(tru) = 0;
            d = d .* sum (E, 2);
            this.State.d = d / sum (d);

          case 'GentleBoost'
            d = this.State.d;
            if (rs)
              [idx, sw, cnt] = boostSample (d, m, this.BagReplace);
              R = compact (RegressionTree (this.X(idx,:), y(idx), ...
                                           'Weights', sw, rtree{:}, ...
                                           this.TreeArgs{:}));
              h = predict (R, this.X);
              this = addLearner (this, R, LR, ...
                                 sum (sw .* (y(idx) - h(idx)) .^ 2));
              this.BagInBag(:,end+1) = cnt > 0;
              this.State.d = boostRescale (d, d .* exp (-LR * y .* h), cnt);
            else
              R = compact (RegressionTree (this.X, y, 'Weights', d, ...
                                           rtree{:}, this.TreeArgs{:}));
              h = predict (R, this.X);
              this = addLearner (this, R, LR, sum (d .* (y - h) .^ 2));
              d = d .* exp (-LR * y .* h);
              this.State.d = d / sum (d);
            endif

          case 'LogitBoost'
            p = 1 ./ (1 + exp (-this.State.F));
            p = min (max (p, eps), 1 - eps);
            z = (double (y > 0) - p) ./ (p .* (1 - p));
            if (rs)
              ## Each row keeps its own score and weight; only the rows drawn
              ## move, as in MATLAB R2024a.
              w = this.State.w;
              [idx, sw, cnt] = boostSample (w, m, this.BagReplace);
              R = compact (RegressionTree (this.X(idx,:), z(idx), ...
                                           'Weights', sw, rtree{:}, ...
                                           this.TreeArgs{:}));
              h = predict (R, this.X);
              this = addLearner (this, R, LR / 2, ...
                                 sum (sw .* (z(idx) - h(idx)) .^ 2));
              this.BagInBag(:,end+1) = cnt > 0;
              u = cnt > 0;
              this.State.F(u) += LR / 2 * h(u);
              pu = 1 ./ (1 + exp (-this.State.F(u)));
              pu = min (max (pu, eps), 1 - eps);
              w1 = w;
              w1(u) = this.State.d0(u) .* pu .* (1 - pu);
              this.State.w = boostRescale (w, w1, cnt);
            else
              d = this.State.d0 .* p .* (1 - p);
              d /= sum (d);
              R = compact (RegressionTree (this.X, z, 'Weights', d, ...
                                           rtree{:}, this.TreeArgs{:}));
              h = predict (R, this.X);
              this = addLearner (this, R, LR / 2, sum (d .* (z - h) .^ 2));
              this.State.F += LR / 2 * h;
            endif

          case 'LPBoost'
            ## As MATLAB R2024a grows LPBoost: the linear program over the
            ## trees with the new one gives the least largest edge; when the
            ## smallest edge is no more than MarginPrecision above it, the new
            ## tree is not kept and the fit stops.  Otherwise the program gives
            ## the learner weights, and as its dual the next observation
            ## weights.
            d = this.State.d;
            T = compact (ClassificationTree (this.X, this.Y, 'Weights', d, ...
                                             'ClassNames', this.ClassNames, ...
                                             ctree{:}, this.TreeArgs{:}));
            mg = treeMargins (T, this.X, this.ClassNames, g);
            e = mg * d;
            ghat = min (this.State.ghat, e);
            U = [this.State.U; mg];
            [gamma, a, dn] = minMaxEdge (U);
            if (ghat - gamma <= this.MarginPrecision)
              this = stopPerfect (this, strcat ("No improvement in the", ...
                                                " last iteration."));
              break;
            endif
            this.Trained{end+1,1} = T;
            this.TrainedWeights = a;
            this.FitInfo(end+1,:) = [mg, e];
            this.NumTrained = numel (this.Trained);
            this.State.d = dn;
            this.State.ghat = ghat;
            this.State.U = U;

          case 'TotalBoost'
            ## As MATLAB R2024a grows TotalBoost: the tree's margins give its
            ## edge; the fit stops when the least largest edge over the trees
            ## kept rises above the smallest edge less MarginPrecision; the
            ## weights on the observations take one quadratic step towards
            ## the least relative entropy to the starting weights, and the
            ## learner weights maximise the smallest margin.
            d = this.State.d;
            T = compact (ClassificationTree (this.X, this.Y, 'Weights', d, ...
                                             'ClassNames', this.ClassNames, ...
                                             ctree{:}, this.TreeArgs{:}));
            mg = treeMargins (T, this.X, this.ClassNames, g);
            e = mg * d;
            ghat = min (this.State.ghat, e);
            U = [this.State.U; mg];
            [gamma, a] = minMaxEdge (U);
            if (gamma > ghat - this.MarginPrecision)
              this = stopPerfect (this, strcat ("No improvement in the", ...
                                                " last iteration."));
              break;
            endif
            this.Trained{end+1,1} = T;
            this.TrainedWeights = a;
            this.FitInfo(end+1,:) = [mg, e];
            this.NumTrained = numel (this.Trained);
            p = max (d, 1e-12);
            H = diag (1 ./ p);
            q = log (p ./ max (this.State.d0, 1e-12)) + 1 - H * p;
            dn = qp (d, H, q, ones (1, n), 1, zeros (n, 1), ones (n, 1), ...
                     [], U, (ghat - this.MarginPrecision) * ones (rows (U), 1));
            ## The solver may leave rounding below zero, which no tree takes.
            dn = max (dn, 0);
            this.State.d = dn / sum (dn);
            this.State.ghat = ghat;
            this.State.U = U;

          case 'Subspace'
            p = columns (this.X);
            if (this.AllCombinations)
              ## nchoosek's order; MATLAB R2024a reverses it for some sizes.
              combos = nchoosek (1:p, this.NPredToSample);
              if (this.NumTrained >= rows (combos))
                break;
              endif
              cols = combos(this.NumTrained + 1,:);
            else
              perm = randperm (p);
              cols = sort (perm(1:this.NPredToSample));
            endif
            largs = [{'PredictorNames', this.PredictorNames(cols), ...
                      'ClassNames', this.ClassNames, 'Prior', this.Prior, ...
                      'Cost', this.Cost}, this.TreeArgs];
            if (strcmp (this.LearnerNames{1}, 'KNN'))
              mdl = ClassificationKNN (this.X(:,cols), this.Y, largs{:});
            else
              mdl = compact (ClassificationDiscriminant (this.X(:,cols), ...
                                                         this.Y, largs{:}));
            endif
            use = false (p, 1);
            use(cols) = true;
            this.UsePredForLearner(:,end+1) = use;
            this = addLearner (this, mdl, 1, []);

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
            present = labelsFromIndex (this.ClassNames, unique (g(idx)));
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
        case 'RUSBoost'
          tail = strcat ("Element t of this vector is the weighted loss", ...
                         " from hypothesis t.");
        case {'LPBoost', 'TotalBoost'}
          l1 = strcat ("Matrix of size NumTrained-by-(N+1), where N is the", ...
                       " number of training observations");
          l2 = strcat ("and NumTrained is the number of learned weak", ...
                       " hypotheses.");
          l3 = strcat ("Row t of this matrix stores N classification", ...
                       " margins from hypothesis t");
          l4 = "and the weighted edge of hypothesis t in the last element.";
          txt = {l1; l2; l3; l4};
          return;
        otherwise
          tail = strcat ("Element t of this vector is the weighted mean", ...
                         " squared error from regression hypothesis t.");
      endswitch
      txt = {head; tail};
    endfunction

  endmethods

endclassdef

## The default learner: nearest neighbours for Subspace, trees otherwise.
function name = ifelse_learner (issub)

  if (issub)
    name = 'knn';
  else
    name = 'tree';
  endif

endfunction

## The rows of one RUSBoost tree: from class k, round (RATIO(k) times the size
## of the smallest class) rows drawn in proportion to the weights D, without
## replacement unless the class holds fewer rows than that.
function idx = rusSample (g, d, ratio)

  K = numel (ratio);
  cnt = accumarray (g, 1, [K, 1]);
  nmin = min (cnt(cnt > 0));
  idx = zeros (0, 1);
  for k = 1:K
    rk = find (g == k);
    m = round (ratio(k) * nmin);
    if (m == 0 || isempty (rk))
      continue;
    endif
    w = d(rk);
    if (! any (w > 0))
      w = ones (size (w));
    endif
    if (m <= numel (rk))
      [~, order] = sort (rand (numel (rk), 1) .^ (1 ./ w), 'descend');
      idx = [idx; rk(order(1:m))];
    else
      cw = [0; cumsum(w)];
      cw /= cw(end);
      idx = [idx; rk(lookup (cw, rand (m, 1)))];
    endif
  endfor
  idx = sort (idx);

endfunction

## The margins of tree T on X, a row: the probability it gives each
## observation's class, of index G, less the largest it gives another.
function mg = treeMargins (T, X, ClassNames, g)

  n = rows (X);
  K = classCount (ClassNames);
  [~, s] = predict (T, X);
  P = zeros (n, K);
  P(:, labelIndices (ClassNames, T.ClassNames)) = s;
  tru = sub2ind ([n, K], (1:n)', g);
  other = P;
  other(tru) = -Inf;
  mg = (P(tru) - max (other, [], 2))';

endfunction

## The least, over distributions D on the observations, of the largest edge
## U * D of the learners, one row of margins each, as GAMMA, with such a D,
## and from its dual the learner weights A that maximise the smallest margin,
## summing to one.
function [gamma, a, d] = minMaxEdge (U)

  [T, n] = size (U);
  c = [zeros(n, 1); 1];
  A = [U, -ones(T, 1); ones(1, n), 0];
  b = [zeros(T, 1); 1];
  [x, gamma, ~, extra] = glpk (c, A, b, [zeros(n, 1); -Inf], [], ...
                               [repmat('U', 1, T), 'S'], ...
                               repmat ('C', 1, n + 1), 1);
  a = max (-extra.lambda(1:T), 0) + 0;
  if (sum (a) > 0)
    a /= sum (a);
  else
    a = ones (T, 1) / T;
  endif
  d = max (x(1:n), 0);
  d /= sum (d);

endfunction

## The weights W after a resampled learner: the rows drawn take their new
## weights W1, rescaled so that over their draws they carry what they carried
## before, and the whole is normalised.
function w = boostRescale (w, w1, cnt)

  u = cnt > 0;
  s1 = sum (cnt(u) .* w1(u));
  if (s1 > 0)
    w(u) = w1(u) * (sum (cnt(u) .* w(u)) / s1);
  endif
  w /= sum (w);

endfunction

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
%!error<ClassificationEnsemble: 'NumLearningCycles' must be a positive integer or 'AllPredictorCombinations'.> ...
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
%!error<ClassificationEnsemble: 'RatioToSmallest' must be a vector of nonnegative numbers with a positive element.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'RUSBoost', ...
%!                         'RatioToSmallest', [0, 0])
%!error<ClassificationEnsemble: 'RatioToSmallest' must be a vector of nonnegative numbers with a positive element.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'RUSBoost', ...
%!                         'RatioToSmallest', [1, NaN])
%!error<ClassificationEnsemble: 'RatioToSmallest' must have one element per class.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'RUSBoost', ...
%!                         'RatioToSmallest', [1, 2, 3])
%!error<ClassificationEnsemble: 'RatioToSmallest' applies only to the 'RUSBoost' method.> ...
%! ClassificationEnsemble (X2, Y2, 'RatioToSmallest', 1)
%!error<ClassificationEnsemble: the 'RUSBoost' method cannot resample the observations.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'RUSBoost', 'Resample', 'on')
%!error<ClassificationEnsemble: 'Learners' must be 'tree' or a tree template.> ...
%! ClassificationEnsemble (X2, Y2, 'Learners', 'knn')
%!error<ClassificationEnsemble: 'Method' must be a character vector.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 1)
%!error<ClassificationEnsemble: 'Boost' is not a valid ensemble method.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'Boost')
%!error<ClassificationEnsemble: the 'RobustBoost' method is not implemented.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'robustboost')
%!error<ClassificationEnsemble: a bagged ensemble is fitted by ClassificationBaggedEnsemble.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'Bag')
%!error<ClassificationEnsemble: the 'AdaBoostM1' method fits exactly two classes.> ...
%! load fisheriris
%! ClassificationEnsemble (meas, species, 'Method', 'AdaBoostM1')
%!error<ClassificationEnsemble: the 'AdaBoostM2' method fits more than two classes.> ...
%! ClassificationEnsemble (X2, Y2, 'Method', 'AdaBoostM2')
%!error<ClassificationEnsemble: a resampled ensemble is fitted by ClassificationBaggedEnsemble.> ...
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

%!test  # MATLAB parity: GentleBoost importance over its regression trees
%! Mdl = ClassificationEnsemble (X2, Y2, 'Method', 'GentleBoost', ...
%!                               'NumLearningCycles', 3, 'Learners', S);
%! assert_equal (predictorImportance (Mdl), ...
%!               [0, 0, 0.222096984079193, 0.388833770583033], 1e-13);

%!test  # MATLAB parity: crossval equals cross-validating at fit time
%! load fisheriris
%! c = cvpartition (Y2, 'KFold', 5);
%! M = ClassificationEnsemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                             'NumLearningCycles', 4, 'Learners', S);
%! CV = crossval (M, 'CVPartition', c);
%! F = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                   'NumLearningCycles', 4, 'Learners', S, 'CVPartition', c);
%! assert_equal (kfoldLoss (CV), kfoldLoss (F), 1e-15);
%! assert_equal (crossval (M).KFold, 10);
%! assert_equal (crossval (M, 'KFold', 3).KFold, 3);

%!error<ClassificationEnsemble.crossval: no partition was asked for; set 'CrossVal' to 'on' or give 'KFold'.> ...
%! crossval (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1), ...
%!           'CrossVal', 'off')
%!error<ClassificationEnsemble.crossval: invalid parameter name in optional pair arguments.> ...
%! crossval (ClassificationEnsemble (X2, Y2, 'NumLearningCycles', 1), 'Foo', 1)

%!error<ClassificationEnsemble: trees cannot be the learners of the 'Subspace' method.> ...
%! load fisheriris
%! ClassificationEnsemble (meas, species, 'Method', 'Subspace', ...
%!                         'Learners', 'tree')
%!error<ClassificationEnsemble: 'Learners' must be 'knn', 'discriminant' or a template of either for the 'Subspace' method.> ...
%! load fisheriris
%! ClassificationEnsemble (meas, species, 'Method', 'Subspace', ...
%!                         'Learners', 'svm')
%!error<ClassificationEnsemble: 'NPredToSample' must be a positive integer.> ...
%! ClassificationEnsemble (X2, Y2, 'NPredToSample', 0)
%!error<ClassificationEnsemble: 'NPredToSample' must be less than the number of predictors.> ...
%! load fisheriris
%! ClassificationEnsemble (meas, species, 'Method', 'Subspace', ...
%!                         'NPredToSample', 4)
%!error<ClassificationEnsemble: 'NPredToSample' applies only to the 'Subspace' method.> ...
%! ClassificationEnsemble (X2, Y2, 'NPredToSample', 2)
%!error<ClassificationEnsemble: 'AllPredictorCombinations' applies only to the 'Subspace' method.> ...
%! ClassificationEnsemble (X2, Y2, 'NumLearningCycles', ...
%!                         'AllPredictorCombinations')
%!error<ClassificationEnsemble: 'LearnRate' cannot be used with the 'Subspace' method.> ...
%! load fisheriris
%! ClassificationEnsemble (meas, species, 'Method', 'Subspace', ...
%!                         'LearnRate', 0.5)
%!error<ClassificationEnsemble: the 'Subspace' method cannot resample the observations.> ...
%! load fisheriris
%! ClassificationEnsemble (meas, species, 'Method', 'Subspace', ...
%!                         'FResample', 0.5)
%!error<ClassificationEnsemble: 'Weights' that are not uniform cannot be used with the 'Subspace' method.> ...
%! load fisheriris
%! ClassificationEnsemble (meas, species, 'Method', 'Subspace', ...
%!                         'Weights', [5 * ones(50, 1); ones(100, 1)])
%!error<ClassificationEnsemble.predictorImportance: predictor importance is defined only for ensembles of trees.> ...
%! load fisheriris
%! predictorImportance (ClassificationEnsemble (meas, species, ...
%!                      'Method', 'Subspace', 'NumLearningCycles', 2))
%!error<ClassificationEnsemble.resume: an ensemble of all predictor combinations cannot grow further.> ...
%! load fisheriris
%! resume (ClassificationEnsemble (meas, species, 'Method', 'Subspace', ...
%!         'NumLearningCycles', 'AllPredictorCombinations'), 1)

%!shared X, yb
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

%!test  # the categorical predictors travel to compact models and folds
%! Mdl = ClassificationEnsemble (X, yb, 'Method', 'AdaBoostM1', ...
%!                               'NumLearningCycles', 3, ...
%!                               'CategoricalPredictors', logical ([1, 0]));
%! assert_equal (compact (Mdl).CategoricalPredictors, 1);
%! CV = crossval (Mdl, 'KFold', 3);
%! assert_equal (CV.CategoricalPredictors, 1);
%! assert_equal (CV.Trained{1}.Trained{1}.CategoricalPredictors, 1);

%!test  # a tree template's categorical options reach every tree
%! t = templateTree ('MaxNumCategories', 3, 'MaxNumSplits', 3);
%! Mdl = ClassificationEnsemble (X, yb, 'Method', 'AdaBoostM1', ...
%!                               'NumLearningCycles', 2, 'Learners', t, ...
%!                               'CategoricalPredictors', 1);
%! assert_equal (Mdl.Trained{1}.CategoricalPredictors, 1);

%!error<ClassificationEnsemble: 'CategoricalPredictors' indices must not exceed the number of predictors.> ...
%! ClassificationEnsemble (X, yb, 'CategoricalPredictors', 3)
