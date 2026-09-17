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

classdef CompactClassificationEnsemble < PredictiveModel
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactClassificationEnsemble
  ##
  ## Compact ensemble of weak learners for classification
  ##
  ## A @code{CompactClassificationEnsemble} object carries the trained weak
  ## learners of a boosted or bagged ensemble and what prediction needs, but
  ## not the observations it was fitted on.  It predicts new data identically
  ## to the ensemble it came from, and weak learners can be removed from it.
  ##
  ## Create one with the @code{compact} method of a
  ## @code{ClassificationEnsemble} or @code{ClassificationBaggedEnsemble}
  ## object.
  ##
  ## @seealso{fitcensemble, ClassificationEnsemble,
  ## ClassificationBaggedEnsemble}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} CategoricalPredictors
    ##
    ## Indices of categorical predictors
    ##
    ## The predictors every tree treats as categorical, empty when none
    ## is.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} ClassNames
    ##
    ## Names of the classes
    ##
    ## The classes, in the type of the response and in the order the columns
    ## of the scores take them.  This property is read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} CombineWeights
    ##
    ## How the weak learners are combined
    ##
    ## @qcode{'WeightedSum'} for a boosted ensemble, whose scores are the sum
    ## of each learner's output times its weight, or @qcode{'WeightedAverage'}
    ## for a bagged one, whose scores are the weighted average of its trees'
    ## class probabilities.  This property is read-only.
    ##
    ## @end deftp
    CombineWeights = 'WeightedSum';

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} Cost
    ##
    ## Misclassification costs
    ##
    ## A square matrix, @code{Cost(i,j)} being the cost of classifying an
    ## observation of class @math{i} as class @math{j}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Cost = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the learners saw them
    ##
    ## The same as @code{PredictorNames}, no predictor being expanded.  This
    ## property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} NumTrained
    ##
    ## Number of trained weak learners
    ##
    ## A nonnegative integer.  This property is read-only.
    ##
    ## @end deftp
    NumTrained = 0;

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} Prior
    ##
    ## Prior probabilities of the classes
    ##
    ## A row vector with one probability per class, in the order of
    ## @code{ClassNames}.  This property is read-only.
    ##
    ## @end deftp
    Prior = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector.  This property is read-only.
    ##
    ## @end deftp
    ResponseName = 'Y';

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} Trained
    ##
    ## Trained weak learners
    ##
    ## A column cell array with one compact model per learner: a
    ## @code{CompactClassificationTree} for AdaBoostM1, AdaBoostM2, RUSBoost and
    ## Bag, and a @code{CompactRegressionTree} for GentleBoost and LogitBoost,
    ## which fit regression trees.  MATLAB wraps those regression trees in a
    ## classifier object of its own; here they are held as they are.  This
    ## property is read-only.
    ##
    ## @end deftp
    Trained = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} TrainedWeights
    ##
    ## Weights of the trained weak learners
    ##
    ## A column with one weight per learner.  This property is read-only.
    ##
    ## @end deftp
    TrainedWeights = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} UsePredForLearner
    ##
    ## Which predictors each learner uses
    ##
    ## For the @qcode{'Subspace'} method, a logical matrix with one row per
    ## predictor and one column per learner.  Empty for tree learners, as
    ## MATLAB returns it.  This property is read-only.
    ##
    ## @end deftp
    UsePredForLearner = [];

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationEnsemble} {property} ScoreTransform
    ##
    ## Transform applied to the predicted scores
    ##
    ## A character vector naming a built-in transform, such as
    ## @qcode{'none'} (default) or @qcode{'doublelogit'}, or a function
    ## handle.  The labels, losses, edges and margins are computed from the
    ## transformed scores.
    ##
    ## @end deftp
    ScoreTransform = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    Method = '';         # the ensemble method, deciding a learner's output
    STfun = @(x) x;      # the score transform as a function
    DefaultIndex = 1;    # class given to a row no learner may score
  endproperties

  methods (Hidden)

    function this = set.ScoreTransform (this, val)
      try
        [this.STfun, this.ScoreTransform] = parseScoreTransform (val, ...
                                            'CompactClassificationEnsemble');
      catch
        error (strcat ("CompactClassificationEnsemble.subsasgn:", ...
                       " 'ScoreTransform' must be a character vector or a", ...
                       " 'function_handle' object."));
      end_try_catch
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationEnsemble} {@var{obj} =} CompactClassificationEnsemble (@var{Mdl})
    ##
    ## Create a @code{CompactClassificationEnsemble} object.
    ##
    ## @var{Mdl} is the @code{ClassificationEnsemble} or
    ## @code{ClassificationBaggedEnsemble} object to compact.  The documented
    ## way to reach this constructor is the @code{compact} method.
    ##
    ## @end deftypefn
    function this = CompactClassificationEnsemble (Mdl = [])

      if (isempty (Mdl))
        return;
      endif
      if (! isa (Mdl, 'ClassificationEnsemble'))
        error (strcat ("CompactClassificationEnsemble: MDL must be a", ...
                       " 'ClassificationEnsemble' object."));
      endif

      this.CategoricalPredictors  = Mdl.CategoricalPredictors;
      this.ClassNames             = Mdl.ClassNames;
      this.CombineWeights         = Mdl.CombineWeights;
      this.Cost                   = Mdl.Cost;
      this.ExpandedPredictorNames = Mdl.ExpandedPredictorNames;
      this.NumTrained             = Mdl.NumTrained;
      this.PredictorNames         = Mdl.PredictorNames;
      this.Prior                  = Mdl.Prior;
      this.ResponseName           = Mdl.ResponseName;
      this.Trained                = Mdl.Trained;
      this.TrainedWeights         = Mdl.TrainedWeights;
      this.UsePredForLearner      = Mdl.UsePredForLearner;
      this.Method                 = Mdl.Method;
      this.DefaultIndex           = Mdl.DefaultIndex;
      ## A transform given as a function handle is held as its text, which
      ## names no built-in transform, so the function itself is copied.
      named = {'doublelogit', 'invlogit', 'ismax', 'logit', 'none', ...
               'identity', 'sign', 'symmetric', 'symmetricismax', ...
               'symmetriclogit'};
      if (any (strcmp (Mdl.ScoreTransform, named)))
        this.ScoreTransform       = Mdl.ScoreTransform;
      else
        this.ScoreTransform       = Mdl.STfun;
      endif

    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ("\n  CompactClassificationEnsemble\n\n");
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %s\n", 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ("%+25s: %s\n", 'ClassNames', classNameListing (this.ClassNames));
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("%+25s: %d\n", 'NumTrained', this.NumTrained);
      fprintf ("\n");
    endfunction

    ## The shared bodies of the scoring methods.  The full and bagged
    ## classes call them with their own names, so an error names the method
    ## that was called.

    function [label, S] = ensemblePredict (this, X, args, caller)

      o = ensembleArgs (this, X, args, {'Learners', 'UseObsForLearner'}, ...
                        caller);
      S = this.STfun (ensembleScores (this, X, o, 'ensemble'));
      label = scoreLabels (this, S);

    endfunction

    function L = ensembleLoss (this, X, Y, args, caller)

      o = ensembleArgs (this, X, args, {'LossFun', 'Learners', 'Mode', ...
                                        'UseObsForLearner', 'Weights'}, ...
                        caller);
      [gY, w] = ensembleResponse (this, X, Y, o, caller);
      S = this.STfun (ensembleScores (this, X, o, o.Mode));
      L = zeros (size (S, 3), 1);
      for k = 1:size (S, 3)
        ## A row no learner may score has no scores to judge, and is left out
        ## with its weight, the rest renormalized, as MATLAB does.
        have = ! any (isnan (S(:,:,k)), 2);
        wk = w(have);
        if (! (sum (wk) > 0))
          L(k) = NaN;
          continue;
        endif
        wk /= sum (wk);
        Sk = S(have,:,k);
        gk = gY(have);
        if (is_function_handle (o.LossFun))
          C = false (size (Sk));
          C(sub2ind (size (Sk), (1:rows (Sk))', gk)) = true;
          L(k) = o.LossFun (C, Sk, wk, this.Cost);
        else
          L(k) = classificationLoss (o.LossFun, Sk, gk, wk, this.Cost);
        endif
      endfor

    endfunction

    function e = ensembleEdge (this, X, Y, args, caller)

      o = ensembleArgs (this, X, args, {'Learners', 'Mode', ...
                                        'UseObsForLearner', 'Weights'}, ...
                        caller);
      [gY, w] = ensembleResponse (this, X, Y, o, caller);
      S = this.STfun (ensembleScores (this, X, o, o.Mode));
      m = marginsOf (S, gY, size (S, 3));
      ## Rows without a margin are left out and the weights renormalized.
      have = ! isnan (m);
      m(! have) = 0;
      W2 = w .* have;
      e = (sum (W2 .* m, 1) ./ sum (W2, 1))(:);

    endfunction

    function m = ensembleMargin (this, X, Y, args, caller)

      o = ensembleArgs (this, X, args, {'Learners', 'UseObsForLearner'}, ...
                        caller);
      gY = ensembleResponse (this, X, Y, o, caller);
      S = this.STfun (ensembleScores (this, X, o, 'ensemble'));
      m = marginsOf (S, gY, 1);

    endfunction

    ## The untransformed scores of the rows of X over the first t learners,
    ## for every t, as an NxKxNumTrained array.
    function S = ensembleSteps (this, X)

      T = this.NumTrained;
      o = struct ('Learners', 1:T, 'U', true (rows (X), T));
      S = ensembleScores (this, X, o, 'cumulative');

    endfunction

    function [imp, ma] = ensembleImportance (this, caller)

      if (strcmp (this.Method, 'Subspace'))
        error (strcat ("%s: predictor importance is defined only for", ...
                       " ensembles of trees."), caller);
      endif
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
    ## @deftypefn  {CompactClassificationEnsemble} {@var{label} =} predict (@var{obj}, @var{X})
    ## @deftypefnx {CompactClassificationEnsemble} {[@var{label}, @var{scores}] =} predict (@dots{})
    ## @deftypefnx {CompactClassificationEnsemble} {[@dots{}] =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Classify new data with a compact ensemble.
    ##
    ## @var{label} holds the class of highest score for each row of @var{X},
    ## in the type of @code{ClassNames}, and @var{scores} the @math{NxK}
    ## scores after @code{ScoreTransform}.
    ##
    ## The scores of a boosted ensemble are the sum over the learners of each
    ## learner's weight times its output.  For AdaBoostM1 the output is +1 for
    ## the class the learner predicts and -1 for the other; for GentleBoost and
    ## LogitBoost it is the regression tree's prediction for the first class and
    ## its negative for the second; for AdaBoostM2 and RUSBoost it is the
    ## learner's class probabilities.  The scores of a bagged ensemble are the
    ## average of its trees' class probabilities.  A row that no learner may
    ## score has @code{NaN} scores and is given the class of greatest prior
    ## probability.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Learners'} @tab @tab A vector of indices of the learners
    ## to use.  The default is all of them.
    ## @item @qcode{'UseObsForLearner'} @tab @tab An @math{NxNumTrained}
    ## logical matrix saying which learner may score which row.  The default
    ## lets every learner score every row.
    ## @end multitable
    ##
    ## @seealso{CompactClassificationEnsemble, fitcensemble}
    ## @end deftypefn
    function [label, scores] = predict (this, X, varargin)

      if (nargin < 2)
        error (strcat ("CompactClassificationEnsemble.predict: too few", ...
                       " input arguments."));
      endif
      [label, scores] = ensemblePredict (this, X, varargin, ...
                          'CompactClassificationEnsemble.predict');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationEnsemble} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationEnsemble} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of a compact ensemble.
    ##
    ## @var{L} is the weighted loss of the scores @code{predict} gives the rows
    ## of @var{X} against the labels @var{Y}.  The weights are normalized so
    ## that each class carries its prior probability.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'LossFun'} @tab @tab @qcode{'classiferror'} (default),
    ## @qcode{'binodeviance'}, @qcode{'classifcost'}, @qcode{'exponential'},
    ## @qcode{'hinge'}, @qcode{'logit'}, @qcode{'mincost'},
    ## @qcode{'quadratic'}, or a function handle called as
    ## @code{lossfun (C, S, W, Cost)}, @var{C} being an @math{NxK} logical
    ## matrix marking each row's class, @var{S} the scores and @var{W} the
    ## normalized weights.
    ## @item @qcode{'Mode'} @tab @tab @qcode{'ensemble'} (default) for one
    ## loss over the learners used, @qcode{'cumulative'} for a column whose
    ## element @math{j} uses the first @math{j} of them, or
    ## @qcode{'individual'} for a column with the loss of each on its own.
    ## @item @qcode{'Weights'} @tab @tab A nonnegative vector with one weight
    ## per row.  The default is uniform.
    ## @end multitable
    ##
    ## @qcode{'Learners'} and @qcode{'UseObsForLearner'} are taken as by
    ## @code{predict}.  A row with @code{NaN} scores, which no learner may
    ## score, is left out and the weights are renormalized over the rest;
    ## @code{edge} does the same.
    ##
    ## @seealso{CompactClassificationEnsemble,
    ## CompactClassificationEnsemble.edge}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      if (nargin < 3)
        error ("CompactClassificationEnsemble.loss: too few input arguments.");
      endif
      L = ensembleLoss (this, X, Y, varargin, ...
                        'CompactClassificationEnsemble.loss');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationEnsemble} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationEnsemble} {@var{e} =} edge (@dots{}, @var{name}, @var{value})
    ##
    ## Classification edge of a compact ensemble.
    ##
    ## @var{e} is the weighted mean of the margins of the rows of @var{X}, the
    ## weights normalized so that each class carries its prior probability.
    ## @qcode{'Mode'} and @qcode{'Weights'} are taken as by @code{loss}, and
    ## @qcode{'Learners'} and @qcode{'UseObsForLearner'} as by @code{predict}.
    ##
    ## @seealso{CompactClassificationEnsemble,
    ## CompactClassificationEnsemble.margin}
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      if (nargin < 3)
        error ("CompactClassificationEnsemble.edge: too few input arguments.");
      endif
      e = ensembleEdge (this, X, Y, varargin, ...
                        'CompactClassificationEnsemble.edge');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationEnsemble} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationEnsemble} {@var{m} =} margin (@dots{}, @var{name}, @var{value})
    ##
    ## Classification margins of a compact ensemble.
    ##
    ## @var{m} holds, for each row of @var{X}, the score of its true class
    ## less the highest score among the other classes.  @qcode{'Learners'}
    ## and @qcode{'UseObsForLearner'} are taken as by @code{predict}.
    ##
    ## @seealso{CompactClassificationEnsemble,
    ## CompactClassificationEnsemble.edge}
    ## @end deftypefn
    function m = margin (this, X, Y, varargin)

      if (nargin < 3)
        error (strcat ("CompactClassificationEnsemble.margin: too few", ...
                       " input arguments."));
      endif
      m = ensembleMargin (this, X, Y, varargin, ...
                          'CompactClassificationEnsemble.margin');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationEnsemble} {@var{imp} =} predictorImportance (@var{obj})
    ## @deftypefnx {CompactClassificationEnsemble} {[@var{imp}, @var{ma}] =} predictorImportance (@var{obj})
    ##
    ## Estimate the importance of each predictor.
    ##
    ## @var{imp} is a row vector with one element per predictor, the average
    ## over the trees of each tree's @code{predictorImportance}, weighted by
    ## @code{TrainedWeights}.  GentleBoost and LogitBoost ensembles take the
    ## importance of their regression trees.  @var{ma}, the predictive
    ## measure of association between the predictors, is empty, the trees
    ## growing no surrogate splits.
    ##
    ## @seealso{CompactClassificationEnsemble}
    ## @end deftypefn
    function [imp, ma] = predictorImportance (this)

      [imp, ma] = ensembleImportance (this, ...
                    'CompactClassificationEnsemble.predictorImportance');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationEnsemble} {@var{C} =} removeLearners (@var{obj}, @var{idx})
    ##
    ## Remove weak learners from a compact ensemble.
    ##
    ## @var{C} is the ensemble without the learners whose indices @var{idx}
    ## holds, their weights and their columns of @code{UsePredForLearner}
    ## removed with them.
    ##
    ## @seealso{CompactClassificationEnsemble}
    ## @end deftypefn
    function this = removeLearners (this, idx)

      if (nargin < 2)
        error (strcat ("CompactClassificationEnsemble.removeLearners: too", ...
                       " few input arguments."));
      endif
      if (! (isnumeric (idx) && isvector (idx) && isreal (idx)
             && all (idx >= 1) && all (idx <= this.NumTrained)
             && all (idx == fix (idx))))
        error (strcat ("CompactClassificationEnsemble.removeLearners: IDX", ...
                       " must be a vector of indices of trained learners."));
      endif
      keep = true (1, this.NumTrained);
      keep(idx) = false;
      this.Trained = this.Trained(keep);
      this.TrainedWeights = this.TrainedWeights(keep);
      if (! isempty (this.UsePredForLearner))
        this.UsePredForLearner = this.UsePredForLearner(:,keep);
      endif
      this.NumTrained = sum (keep);

    endfunction

  endmethods

  methods (Access = private)

    ## The Name-Value arguments of the scoring methods, those not in ALLOWED
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
                  'Mode', 'ensemble', 'Weights', [], ...
                  'LossFun', 'classiferror');
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
                   && numel (val) == rows (X) && all (val >= 0)
                   && any (val > 0)))
              error (strcat ("%s: 'Weights' must be a nonnegative numeric", ...
                             " vector with one element per observation,", ...
                             " not all zero."), caller);
            endif
            o.Weights = double (val(:));
          case 'lossfun'
            losses = {'binodeviance', 'classifcost', 'classiferror', ...
                      'exponential', 'hinge', 'logit', 'mincost', ...
                      'quadratic'};
            if (ischar (val) && any (strcmpi (val, losses)))
              o.LossFun = tolower (val);
            elseif (is_function_handle (val))
              o.LossFun = val;
            else
              error (strcat ("%s: 'LossFun' must be the name of a", ...
                             " classification loss or a function handle."), ...
                     caller);
            endif
        endswitch
      endfor

    endfunction

    ## The class index of each row of Y and the weights, normalized so that
    ## each class carries its prior and the whole sums to one.
    function [gY, w] = ensembleResponse (this, X, Y, o, caller)

      if (rows (Y) != rows (X))
        error ("%s: X and Y must have the same number of rows.", caller);
      endif
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("%s: %s", caller, errmsg);
      endif
      w = o.Weights;
      if (isempty (w))
        w = ones (rows (X), 1);
      endif
      w = priorNormalize (w, gY, this.Prior);
      w /= sum (w);

    endfunction

    ## What learner T says about each row of X, laid out over the classes.
    function G = learnerOutput (this, t, X)

      K = classCount (this.ClassNames);
      mdl = this.Trained{t};
      switch (this.Method)
        case 'AdaBoostM1'
          g = labelIndices (this.ClassNames, predict (mdl, X));
          h = double (g == 1) - double (g == 2);
          G = [h, -h];
        case {'GentleBoost', 'LogitBoost'}
          h = predict (mdl, X);
          G = [h, -h];
        case {'LPBoost', 'TotalBoost'}
          ## An LPBoost or TotalBoost learner scores a class with twice its
          ## probability less one, as in MATLAB R2024a.
          [~, s] = predict (mdl, X);
          G = -ones (rows (X), K);
          G(:, labelIndices (this.ClassNames, mdl.ClassNames)) = 2 * s - 1;
        case 'Subspace'
          [~, s] = predict (mdl, X(:,this.UsePredForLearner(:,t)));
          G = zeros (rows (X), K);
          G(:, labelIndices (this.ClassNames, mdl.ClassNames)) = s;
        otherwise
          [~, s] = predict (mdl, X);
          G = zeros (rows (X), K);
          G(:, labelIndices (this.ClassNames, mdl.ClassNames)) = s;
      endswitch

    endfunction

    ## The untransformed scores over the learners O names: NxK for MODE
    ## 'ensemble', NxKxT otherwise, T being the number of learners.  A row no
    ## learner may score has NaN scores.
    function S = ensembleScores (this, X, o, mode)

      n = rows (X);
      K = classCount (this.ClassNames);
      T = numel (o.Learners);
      A = zeros (n, K, T);
      V = zeros (n, 1, T);
      for j = 1:T
        t = o.Learners(j);
        u = o.U(:,t);
        A(:,:,j) = this.TrainedWeights(t) * learnerOutput (this, t, X) .* u;
        V(:,1,j) = this.TrainedWeights(t) * u;
      endfor
      average = strcmp (this.CombineWeights, 'WeightedAverage');
      switch (mode)
        case 'ensemble'
          num = sum (A, 3);
          den = sum (V, 3);
          used = sum (reshape (o.U(:,o.Learners), n, 1, T), 3) > 0;
        case 'cumulative'
          num = cumsum (A, 3);
          den = cumsum (V, 3);
          used = cumsum (reshape (o.U(:,o.Learners), n, 1, T), 3) > 0;
        otherwise
          num = A;
          den = V;
          used = reshape (o.U(:,o.Learners), n, 1, T);
      endswitch
      if (average)
        S = num ./ den;
      else
        S = num;
      endif
      S(repmat (! used, 1, K)) = NaN;

    endfunction

    ## The class of highest score, and for a row with NaN scores the class of
    ## greatest prior probability.
    function label = scoreLabels (this, S)

      [~, k] = max (S, [], 2);
      k(any (isnan (S), 2)) = this.DefaultIndex;
      label = labelsFromIndex (this.ClassNames, k);

    endfunction

  endmethods

endclassdef

## Test output
%!shared X2, Y2, C
%! load fisheriris
%! X2 = meas(51:150,:);
%! Y2 = species(51:150);
%! C = compact (fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                            'NumLearningCycles', 5, ...
%!                            'Learners', templateTree ('MaxNumSplits', 1)));

%!test  # MATLAB parity: the properties of a compact ensemble
%! assert_equal (numel (properties (C)), 13);
%! assert_equal (C.NumTrained, 5);
%! assert_equal (C.CombineWeights, 'WeightedSum');
%! assert_equal (isprop (C, 'X'), false);

%!test  # MATLAB parity: scores over a subset of the learners
%! [~, s] = predict (C, X2([1, 51],:), 'Learners', [2, 4]);
%! assert_equal (s(:,1), [1.547898302883168; -0.439169918665653], 1e-13);

%!test  # MATLAB parity: a row no learner may score has NaN scores
%! U = true (2, 5);
%! U(1,:) = false;
%! U(2,[1, 3]) = false;
%! [label, s] = predict (C, X2(1:2,:), 'UseObsForLearner', U);
%! assert_equal (isnan (s(1,:)), [true, true]);
%! assert_equal (label{1}, 'versicolor');
%! assert_equal (s(2,1), 1.279703855451121, 1e-13);

%!test  # MATLAB parity: such a row takes the class of greatest prior
%! load fisheriris
%! Y = species(41:150);
%! Y(1:10) = {'virginica'};
%! M = compact (fitcensemble (meas(41:150,:), Y, 'Method', 'AdaBoostM1', ...
%!                            'NumLearningCycles', 3, ...
%!                            'Learners', templateTree ('MaxNumSplits', 1)));
%! label = predict (M, meas(1,:), 'UseObsForLearner', false (1, 3));
%! assert_equal (label, {'virginica'});

%!test  # MATLAB parity: the loss in its three modes
%! assert_equal (loss (C, X2, Y2), 0.04, 1e-15);
%! assert_equal (loss (C, X2, Y2, 'Mode', 'cumulative'), ...
%!               [0.06; 0.06; 0.04; 0.06; 0.04], 1e-15);
%! assert_equal (loss (C, X2, Y2, 'Mode', 'individual'), ...
%!               [0.06; 0.08; 0.21; 0.5; 0.5], 1e-14);
%! assert_equal (loss (C, X2, Y2, 'Learners', [1, 3]), 0.06, 1e-15);

%!test  # MATLAB parity: every built-in loss function
%! f = {'binodeviance', 'exponential', 'hinge', 'logit', 'quadratic', ...
%!      'mincost', 'classifcost'};
%! L = cellfun (@(n) loss (C, X2, Y2, 'LossFun', n), f);
%! assert_equal (L, [0.088288301419851, 0.182079166822146, ...
%!                   0.093152538904937, 0.136695159561291, ...
%!                   3.600631106449406, 0.04, 0.04], 1e-12);

%!test  # a loss function given as a handle
%! f = @(Cl, Sc, W, Cost) sum (W .* (Sc(Cl) < 0));
%! assert_equal (loss (C, X2, Y2, 'LossFun', f), 0.04, 1e-15);

%!test  # MATLAB parity: the edge and its cumulative mode
%! assert_equal (edge (C, X2, Y2), 5.115272254726563, 1e-12);
%! assert_equal (edge (C, X2, Y2, 'Mode', 'cumulative'), ...
%!               [2.421351075476915; 4.090488381577923; ...
%!                5.115272254726563; 5.115272254726564; ...
%!                5.115272254726563], 1e-12);

%!test  # MATLAB parity: margins are twice the score of the true class
%! [~, s] = predict (C, X2([1, 21],:));
%! m = margin (C, X2([1, 21],:), Y2([1, 21]));
%! assert_equal (m, [3.544074277136197; -1.958996348947703], 1e-12);
%! assert_equal (m, 2 * s(:,1), 1e-14);

%!test  # MATLAB parity: losses read the transformed scores
%! D = C;
%! D.ScoreTransform = 'doublelogit';
%! assert_equal (loss (D, X2, Y2, 'LossFun', 'quadratic'), ...
%!               0.028390700687515, 1e-13);

%!test  # a score transform given as a function handle survives compacting
%! load fisheriris
%! Mdl = fitcensemble (X2, Y2, 'Method', 'AdaBoostM1', ...
%!                     'NumLearningCycles', 2, ...
%!                     'Learners', templateTree ('MaxNumSplits', 1), ...
%!                     'ScoreTransform', @(s) 2 * s);
%! [~, s] = predict (Mdl, X2(1,:));
%! [~, s0] = predict (C, X2(1,:), 'Learners', [1, 2]);
%! assert_equal (s, 2 * s0, 1e-14);

%!test  # MATLAB parity: removing learners
%! D = removeLearners (C, [2, 4]);
%! assert_equal (D.NumTrained, 3);
%! assert_equal (D.TrainedWeights, [1.375767656520975; ...
%!               0.883434373403998; 0.268194447432047], 1e-13);

## Test input validation
%!error<CompactClassificationEnsemble: MDL must be a 'ClassificationEnsemble' object.> ...
%! CompactClassificationEnsemble (1)
%!error<CompactClassificationEnsemble.predict: too few input arguments.> ...
%! predict (C)
%!error<CompactClassificationEnsemble.predict: X must be a real numeric matrix.> ...
%! predict (C, {1})
%!error<CompactClassificationEnsemble.predict: X must have one column per predictor.> ...
%! predict (C, ones (2, 3))
%!error<CompactClassificationEnsemble.predict: name-value arguments must be in pairs.> ...
%! predict (C, X2, 'Learners')
%!error<CompactClassificationEnsemble.predict: invalid parameter name in optional pair arguments.> ...
%! predict (C, X2, 'Mode', 'ensemble')
%!error<CompactClassificationEnsemble.predict: 'Learners' must be a vector of indices of trained learners.> ...
%! predict (C, X2, 'Learners', 6)
%!error<CompactClassificationEnsemble.predict: 'UseObsForLearner' must be a logical matrix with one row per observation and one column per trained learner.> ...
%! predict (C, X2, 'UseObsForLearner', true (2, 5))
%!error<CompactClassificationEnsemble.loss: too few input arguments.> ...
%! loss (C, X2)
%!error<CompactClassificationEnsemble.loss: 'Mode' must be 'ensemble', 'cumulative' or 'individual'.> ...
%! loss (C, X2, Y2, 'Mode', 'all')
%!error<CompactClassificationEnsemble.loss: 'Weights' must be a nonnegative numeric vector with one element per observation, not all zero.> ...
%! loss (C, X2, Y2, 'Weights', ones (3, 1))
%!error<CompactClassificationEnsemble.loss: 'LossFun' must be the name of a classification loss or a function handle.> ...
%! loss (C, X2, Y2, 'LossFun', 'mse')
%!error<CompactClassificationEnsemble.loss: X and Y must have the same number of rows.> ...
%! loss (C, X2, Y2(1:3))
%!error<CompactClassificationEnsemble.loss: Y must hold only classes the model was trained on.> ...
%! loss (C, X2(1:2,:), {'rose'; 'versicolor'})
%!error<CompactClassificationEnsemble.edge: too few input arguments.> ...
%! edge (C, X2)
%!error<CompactClassificationEnsemble.edge: invalid parameter name in optional pair arguments.> ...
%! edge (C, X2, Y2, 'LossFun', 'hinge')
%!error<CompactClassificationEnsemble.margin: too few input arguments.> ...
%! margin (C, X2)
%!error<CompactClassificationEnsemble.margin: invalid parameter name in optional pair arguments.> ...
%! margin (C, X2, Y2, 'Mode', 'cumulative')
%!error<CompactClassificationEnsemble.removeLearners: too few input arguments.> ...
%! removeLearners (C)
%!error<CompactClassificationEnsemble.removeLearners: IDX must be a vector of indices of trained learners.> ...
%! removeLearners (C, 0)
%!error<CompactClassificationEnsemble.subsasgn: 'ScoreTransform' must be a character vector or a 'function_handle' object.> ...
%! D = C;
%! D.ScoreTransform = 1;

%!test  # MATLAB parity: importance is the weighted average over the learners
%! [imp, ma] = predictorImportance (C);
%! assert_equal (imp, [0.020544313498471, 0, 0.094343052487178, ...
%!                     0.131555749375574], 1e-13);
%! assert_equal (ma, []);

%!test  # MATLAB parity: the importance of AdaBoostM2 stumps
%! load fisheriris
%! M = fitcensemble (meas, species, 'Method', 'AdaBoostM2', ...
%!                   'NumLearningCycles', 4, ...
%!                   'Learners', templateTree ('MaxNumSplits', 1));
%! assert_equal (predictorImportance (compact (M)), ...
%!               [0, 0, 0.230674233368758, 0.072325623666526], 1e-13);

%!test  # MATLAB parity: boosted regression trees give their importance
%! M = fitcensemble (X2, Y2, 'Method', 'LogitBoost', 'NumLearningCycles', 3, ...
%!                   'Learners', templateTree ('MaxNumSplits', 1), ...
%!                   'LearnRate', 0.5);
%! assert_equal (predictorImportance (M), ...
%!               [0, 0, 0.698823105603178, 1.563129005828868], 1e-13);

%!test  # MATLAB parity: a row no learner may score is left out of the loss
%! U = true (100, 5);
%! U(1,:) = false;
%! miss = ! strcmp (predict (C, X2(2:end,:)), Y2(2:end));
%! assert_equal (loss (C, X2, Y2, 'UseObsForLearner', U), mean (miss), 1e-15);
%! m = margin (C, X2(2:end,:), Y2(2:end));
%! assert_equal (edge (C, X2, Y2, 'UseObsForLearner', U), mean (m), 1e-13);
