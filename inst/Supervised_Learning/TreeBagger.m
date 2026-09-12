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

classdef TreeBagger
  ## -*- texinfo -*-
  ## @deftp {statistics} TreeBagger
  ##
  ## Ensemble of bagged decision trees
  ##
  ## A @code{TreeBagger} object is a random forest: an ensemble of decision
  ## trees, each grown on a bootstrap sample of the training data and each
  ## choosing every split from a random subset of the predictors.  The
  ## ensemble predicts by averaging its trees, their class probabilities for
  ## classification and their responses for regression.
  ##
  ## Every observation a tree's sample leaves out is out of bag for that tree,
  ## and the out-of-bag methods judge the ensemble on those observations
  ## alone, giving an estimate of its error on new data without a separate
  ## test set.
  ##
  ## Create one with the @code{TreeBagger} constructor.  The @code{compact}
  ## method drops the training data and returns a @code{CompactTreeBagger}.
  ##
  ## @seealso{CompactTreeBagger, ClassificationTree, RegressionTree, fitctree,
  ## fitrtree}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} Method
    ##
    ## Type of the ensemble
    ##
    ## @qcode{'classification'} or @qcode{'regression'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Method = 'classification';

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} NumTrees
    ##
    ## Number of trees
    ##
    ## A positive integer, the number of trees in the ensemble.  This property
    ## is read-only.
    ##
    ## @end deftp
    NumTrees = 0;

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} Trees
    ##
    ## The trees of the ensemble
    ##
    ## A column cell array holding one @code{CompactClassificationTree} or
    ## @code{CompactRegressionTree} object per tree.  This property is
    ## read-only.
    ##
    ## @end deftp
    Trees = {};

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} X
    ##
    ## Predictor data
    ##
    ## The predictors the ensemble was fitted on, one row per observation.  A
    ## row whose response is missing is not kept.  This property is
    ## read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} Y
    ##
    ## Response data
    ##
    ## The response the ensemble was fitted on, in the type it was given in,
    ## without the observations whose response is missing.  This property is
    ## read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} W
    ##
    ## Observation weights
    ##
    ## A column of weights summing to one, one per observation.  For
    ## classification each class's weights sum to its prior.  The bootstrap
    ## draws observations in proportion to these weights, and the out-of-bag
    ## error is weighted by them.  This property is read-only.
    ##
    ## @end deftp
    W = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} ClassNames
    ##
    ## Names of the classes
    ##
    ## The classes of a classification ensemble, in the type of the response
    ## and in the order its scores are laid out: sorted, or the order given by
    ## @qcode{'ClassNames'}.  Empty for a regression ensemble.  This property
    ## is read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} Prior
    ##
    ## Prior probabilities of the classes
    ##
    ## A row vector with one probability per class, in the order of
    ## @code{ClassNames}.  Every tree is grown with this prior, restricted to
    ## the classes its sample holds.  Empty for a regression ensemble.  This
    ## property is read-only.
    ##
    ## @end deftp
    Prior = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} Cost
    ##
    ## Misclassification costs
    ##
    ## A square matrix, @code{Cost(i,j)} being the cost of classifying an
    ## observation of class @math{i} as class @math{j}.  The trees are grown
    ## with it; the ensemble's label is the class of highest average score, as
    ## MATLAB documents.  Empty for a regression ensemble.  This property is
    ## read-only.
    ##
    ## @end deftp
    Cost = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} DefaultYfit
    ##
    ## Prediction for an observation no tree may answer for
    ##
    ## For classification, the class of greatest prior probability, in the
    ## type of @code{ClassNames}; for regression, the weighted mean of the
    ## response.  An out-of-bag prediction takes it for an observation that is
    ## in the sample of every tree.  This property is read-only.
    ##
    ## @end deftp
    DefaultYfit = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors naming the columns of @code{X}.  This
    ## property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} NumPredictorsToSample
    ##
    ## Predictors each split is chosen from
    ##
    ## A positive integer, or @qcode{'all'}.  The default is the square root
    ## of the number of predictors for classification and a third of it for
    ## regression, rounded up.  This property is read-only.
    ##
    ## @end deftp
    NumPredictorsToSample = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} MinLeafSize
    ##
    ## Fewest observations a leaf may hold
    ##
    ## A positive integer, 1 by default for classification and 5 for
    ## regression.  A node is split only when it holds at least twice as many.
    ## This property is read-only.
    ##
    ## @end deftp
    MinLeafSize = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} InBagFraction
    ##
    ## Size of each tree's sample
    ##
    ## A number greater than 0 and no greater than 1.  Each tree is grown on
    ## @code{ceil (InBagFraction * N)} observations.  This property is
    ## read-only.
    ##
    ## @end deftp
    InBagFraction = 1;

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} SampleWithReplacement
    ##
    ## Whether the samples are drawn with replacement
    ##
    ## A logical scalar, true by default.  This property is read-only.
    ##
    ## @end deftp
    SampleWithReplacement = true;

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} ComputeOOBPrediction
    ##
    ## Whether out-of-bag information is kept
    ##
    ## A logical scalar, false unless the ensemble was fitted with
    ## @qcode{'OOBPrediction'} set to @qcode{'on'}.  The out-of-bag methods
    ## need it.  This property is read-only.
    ##
    ## @end deftp
    ComputeOOBPrediction = false;

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} ComputeOOBPredictorImportance
    ##
    ## Whether out-of-bag predictor importance is computed
    ##
    ## Always false: out-of-bag predictor importance is not implemented.  This
    ## property is read-only.
    ##
    ## @end deftp
    ComputeOOBPredictorImportance = false;

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} OOBIndices
    ##
    ## Which observations each tree left out
    ##
    ## An @math{NxNumTrees} logical matrix, true where an observation is not
    ## in a tree's sample.  Empty unless @code{ComputeOOBPrediction} is true.
    ## This property is read-only.
    ##
    ## @end deftp
    OOBIndices = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} OOBInstanceWeight
    ##
    ## Number of trees each observation is out of bag for
    ##
    ## A column with one count per observation.  Empty unless
    ## @code{ComputeOOBPrediction} is true.  This property is read-only.
    ##
    ## @end deftp
    OOBInstanceWeight = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} MergeLeaves
    ##
    ## Whether the trees merge leaves
    ##
    ## A logical scalar, false unless @qcode{'MergeLeaves'} was set to
    ## @qcode{'on'}.  This property is read-only.
    ##
    ## @end deftp
    MergeLeaves = false;

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} Prune
    ##
    ## Whether the trees estimate a pruning sequence
    ##
    ## A logical scalar, false unless @qcode{'Prune'} was set to @qcode{'on'}.
    ## The trees are never pruned.  This property is read-only.
    ##
    ## @end deftp
    Prune = false;

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} TreeArguments
    ##
    ## Options given to the trees
    ##
    ## A cell array of the Name-Value pairs given to the constructor that are
    ## passed on to every tree, in the order they were given.  This property is
    ## read-only.
    ##
    ## @end deftp
    TreeArguments = {};

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    TreeClassIdx = {};   # columns of ClassNames each tree's scores fill
    DefaultIndex = 0;    # index of DefaultYfit into ClassNames
    DefaultScore = [];   # scores of an observation no tree may answer for
    GrowArgs = {};       # Name-Value pairs every tree is grown with
    gY = [];             # class index of each observation
  endproperties

  methods (Hidden)

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ('\n  TreeBagger\n\n');
      fprintf ('%30s: %s\n', 'Method', this.Method);
      fprintf ('%30s: %d\n', 'NumTrees', this.NumTrees);
      fprintf ('%30s: [%dx%d]\n', 'X', rows (this.X), columns (this.X));
      fprintf ('%30s: %d\n', 'NumPredictors', columns (this.X));
      fprintf ('%30s: %s\n', 'NumPredictorsToSample', ...
               num2str (this.NumPredictorsToSample));
      fprintf ('%30s: %d\n', 'MinLeafSize', this.MinLeafSize);
      fprintf ('%30s: %g\n', 'InBagFraction', this.InBagFraction);
      fprintf ('%30s: %d\n', 'SampleWithReplacement', ...
               this.SampleWithReplacement);
      fprintf ('%30s: %d\n', 'ComputeOOBPrediction', ...
               this.ComputeOOBPrediction);
      if (strcmp (this.Method, 'classification'))
        fprintf ('%30s: %s\n', 'ClassNames', ...
                 classNameListing (this.ClassNames));
      endif
      fprintf ('\n');
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{B} =} TreeBagger (@var{NumTrees}, @var{X}, @var{Y})
    ## @deftypefnx {TreeBagger} {@var{B} =} TreeBagger (@dots{}, @var{name}, @var{value})
    ##
    ## Grow an ensemble of bagged decision trees.
    ##
    ## @code{@var{B} = TreeBagger (@var{NumTrees}, @var{X}, @var{Y})} grows
    ## @var{NumTrees} classification trees on the @math{NxP} predictor matrix
    ## @var{X} and the response @var{Y}.  @var{Y} holds a class label per row,
    ## as a numeric or logical vector, a categorical, string or character
    ## array, or a cell array of character vectors.  An observation whose
    ## response is missing is left out; one missing a predictor is kept.
    ##
    ## Each tree is grown on a sample of @code{ceil (InBagFraction * N)}
    ## observations, drawn with replacement in proportion to the weights, and
    ## chooses every split from @code{NumPredictorsToSample} predictors drawn
    ## afresh at each node.  The trees are neither pruned nor merged unless
    ## asked.  The random numbers come from Octave's generator, so @code{rng}
    ## reproduces an ensemble.
    ##
    ## Name-Value arguments of the ensemble:
    ##
    ## @multitable @columnfractions 0.3 0.02 0.68
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Method'} @tab @tab @qcode{'classification'} (default) or
    ## @qcode{'regression'}, for a real numeric response.
    ## @item @qcode{'NumPredictorsToSample'} @tab @tab A positive integer or
    ## @qcode{'all'}.  The default is @code{ceil (sqrt (P))} for classification
    ## and @code{ceil (P / 3)} for regression.
    ## @item @qcode{'MinLeafSize'} @tab @tab A positive integer, 1 by default
    ## for classification and 5 for regression.
    ## @item @qcode{'InBagFraction'} @tab @tab The share of the observations in
    ## each sample, greater than 0 and no greater than 1.  The default is 1.
    ## @item @qcode{'SampleWithReplacement'} @tab @tab @qcode{'on'} (default)
    ## or @qcode{'off'}.
    ## @item @qcode{'OOBPrediction'} @tab @tab @qcode{'off'} (default) or
    ## @qcode{'on'}, to keep what the out-of-bag methods need.  Sampling
    ## without replacement at an @qcode{'InBagFraction'} of 1 leaves nothing
    ## out of bag, so the two are refused together.
    ## @item @qcode{'Weights'} @tab @tab A nonnegative vector with one weight
    ## per observation.  The default is uniform.
    ## @item @qcode{'Prior'} @tab @tab @qcode{'empirical'} (default),
    ## @qcode{'uniform'}, a vector with one probability per class, or a
    ## structure with fields @qcode{ClassNames} and @qcode{ClassProbs}.
    ## Classification only.
    ## @item @qcode{'Cost'} @tab @tab A square matrix of misclassification
    ## costs, or a structure with fields @qcode{ClassNames} and
    ## @qcode{ClassificationCosts}.  Classification only.
    ## @item @qcode{'ClassNames'} @tab @tab The classes to fit, in the order
    ## their scores are to be laid out; observations of other classes are left
    ## out.  Classification only.
    ## @item @qcode{'PredictorNames'} @tab @tab A cell array of character
    ## vectors naming the columns of @var{X}.
    ## @item @qcode{'NumPrint'} @tab @tab A nonnegative integer.  After every
    ## that many trees a line saying how many are done is printed.  The
    ## default, 0, prints nothing.
    ## @end multitable
    ##
    ## @qcode{'MaxNumSplits'}, @qcode{'MergeLeaves'}, @qcode{'Prune'},
    ## @qcode{'PruneCriterion'} and @qcode{'SplitCriterion'} are passed on to
    ## every tree, and so is @qcode{'QuadraticErrorTolerance'} for regression;
    ## see @code{fitctree} and @code{fitrtree}.  Merging leaves is allowed but
    ## warned against.
    ##
    ## Categorical predictors, surrogate splits, out-of-bag predictor
    ## importance, parallel growth and tall arrays are not implemented, and an
    ## option asking for one of them is refused.
    ##
    ## MATLAB returns classification labels as a cell array of character
    ## vectors whatever the type of the response; this ensemble returns them,
    ## and its @code{ClassNames} and @code{DefaultYfit}, in the type of the
    ## response, as every other classifier in this package does.  Code that
    ## converts MATLAB's labels with @code{str2double} will get @code{NaN}
    ## from a numeric response here, its labels already being numbers.
    ##
    ## @seealso{TreeBagger, CompactTreeBagger, fitctree, fitrtree}
    ## @end deftypefn
    function this = TreeBagger (NumTrees, X, Y, varargin)

      ## Input validation
      if (nargin < 3)
        error ("TreeBagger: too few input arguments.");
      endif
      if (! (isnumeric (NumTrees) && isscalar (NumTrees)
             && isreal (NumTrees) && NumTrees >= 1
             && NumTrees == fix (NumTrees)))
        error ("TreeBagger: NUMTREES must be a positive integer.");
      endif
      if (! (isnumeric (X) && isreal (X) && ismatrix (X) && ! isempty (X)))
        error ("TreeBagger: X must be a non-empty real numeric matrix.");
      endif
      if (ischar (Y))
        nY = rows (Y);
      else
        nY = numel (Y);
      endif
      if (! (isvector (Y) || ischar (Y)) || nY != rows (X))
        error (strcat ("TreeBagger: Y must be a vector with one element", ...
                       " per row in X."));
      endif
      if (! ischar (Y))
        Y = Y(:);
      endif
      if (mod (numel (varargin), 2) != 0)
        error ("TreeBagger: name-value arguments must be in pairs.");
      endif

      ## Defaults
      Method = 'classification';
      NumVarSample = [];
      MinLeafSize = [];
      InBagFraction = 1;
      WithReplacement = true;
      OOBPrediction = false;
      Cost = [];
      Prior = 'empirical';
      ClassNames = [];
      Weights = [];
      PredictorNames = {};
      NumPrint = 0;
      TreeArgs = {};
      givenClass = {};
      givenReg = {};

      ## Parse optional parameters
      for i = 1:2:numel (varargin)
        name = varargin{i};
        Value = varargin{i+1};
        if (! ischar (name))
          error (strcat ("TreeBagger: invalid parameter name in optional", ...
                         " pair arguments."));
        endif
        switch (tolower (name))

          case 'method'
            if (! (ischar (Value)
                   && any (strcmpi (Value, {'classification', ...
                                            'regression'}))))
              error (strcat ("TreeBagger: 'Method' must be", ...
                             " 'classification' or 'regression'."));
            endif
            Method = tolower (Value);

          case 'numpredictorstosample'
            if (ischar (Value) && strcmpi (Value, 'all'))
              NumVarSample = 'all';
            elseif (isnumeric (Value) && isscalar (Value) && isreal (Value)
                    && Value >= 1 && Value == fix (Value))
              NumVarSample = double (Value);
            else
              error (strcat ("TreeBagger: 'NumPredictorsToSample' must be", ...
                             " a positive integer or 'all'."));
            endif

          case 'minleafsize'
            if (! (isnumeric (Value) && isscalar (Value) && isreal (Value)
                   && Value >= 1 && Value == fix (Value)))
              error (strcat ("TreeBagger: 'MinLeafSize' must be a", ...
                             " positive integer."));
            endif
            MinLeafSize = double (Value);
            TreeArgs(end+1:end+2) = {name, Value};

          case 'inbagfraction'
            if (! (isnumeric (Value) && isscalar (Value) && isreal (Value)
                   && Value > 0 && Value <= 1))
              error (strcat ("TreeBagger: 'InBagFraction' must be a", ...
                             " number greater than 0 and no greater", ...
                             " than 1."));
            endif
            InBagFraction = double (Value);

          case 'samplewithreplacement'
            [WithReplacement, ok] = onOff (Value);
            if (! ok)
              error (strcat ("TreeBagger: 'SampleWithReplacement' must", ...
                             " be 'on' or 'off'."));
            endif

          case 'oobprediction'
            [OOBPrediction, ok] = onOff (Value);
            if (! ok)
              error ("TreeBagger: 'OOBPrediction' must be 'on' or 'off'.");
            endif

          case 'oobpredictorimportance'
            [imp, ok] = onOff (Value);
            if (! ok)
              error (strcat ("TreeBagger: 'OOBPredictorImportance' must", ...
                             " be 'on' or 'off'."));
            endif
            if (imp)
              error (strcat ("TreeBagger: 'OOBPredictorImportance' is not", ...
                             " implemented."));
            endif

          case 'cost'
            Cost = Value;
            givenClass{end+1} = name;

          case 'prior'
            Prior = Value;
            givenClass{end+1} = name;

          case 'classnames'
            ClassNames = Value;
            givenClass{end+1} = name;

          case 'weights'
            Weights = Value;

          case 'predictornames'
            if (! (iscellstr (Value) && numel (Value) == columns (X)))
              error (strcat ("TreeBagger: 'PredictorNames' must be a", ...
                             " cell array of character vectors with one", ...
                             " element per column of X."));
            endif
            PredictorNames = Value(:)';

          case 'categoricalpredictors'
            if (! isempty (Value))
              error (strcat ("TreeBagger: 'CategoricalPredictors' is not", ...
                             " implemented."));
            endif

          case 'numprint'
            if (! (isnumeric (Value) && isscalar (Value) && isreal (Value)
                   && Value >= 0 && Value == fix (Value)))
              error (strcat ("TreeBagger: 'NumPrint' must be a", ...
                             " nonnegative integer."));
            endif
            NumPrint = double (Value);

          case {'maxnumsplits', 'mergeleaves', 'prune', 'prunecriterion', ...
                'splitcriterion'}
            TreeArgs(end+1:end+2) = {name, Value};

          case 'quadraticerrortolerance'
            TreeArgs(end+1:end+2) = {name, Value};
            givenReg{end+1} = name;

          case {'surrogate', 'predictorselection', ...
                'algorithmforcategorical', 'maxnumcategories', 'options', ...
                'chunksize'}
            error ("TreeBagger: '%s' is not implemented.", name);

          otherwise
            error (strcat ("TreeBagger: invalid parameter name in", ...
                           " optional pair arguments."));

        endswitch
      endfor

      isclass = strcmp (Method, 'classification');
      if (isclass && ! isempty (givenReg))
        error ("TreeBagger: '%s' applies only to regression.", givenReg{1});
      endif
      if (! isclass && ! isempty (givenClass))
        error ("TreeBagger: '%s' applies only to classification.", ...
               givenClass{1});
      endif
      if (! WithReplacement && InBagFraction == 1 && OOBPrediction)
        error (strcat ("TreeBagger: sampling without replacement at an", ...
                       " 'InBagFraction' of 1 leaves no observation out", ...
                       " of bag for 'OOBPrediction'."));
      endif

      ## Observation weights, before any row is left out
      N = rows (X);
      if (isempty (Weights))
        RawW = ones (N, 1);
      else
        if (! (isnumeric (Weights) && isvector (Weights) && isreal (Weights)
               && numel (Weights) == N && all (Weights >= 0)))
          error (strcat ("TreeBagger: 'Weights' must be a nonnegative", ...
                         " numeric vector with one element per row in X."));
        endif
        RawW = double (Weights(:));
      endif

      if (isclass)

        ## The classes, in the type of the response, and the class of each
        ## observation.  A named subset keeps its own order and leaves the
        ## other classes' observations out, as MATLAB does.
        g0 = grp2idx (Y);
        U = uniqueLabels (Y(! isnan (g0),:));
        if (isempty (ClassNames))
          C = U;
        else
          if (! (iscellstr (ClassNames) || isnumeric (ClassNames)
                 || islogical (ClassNames) || ischar (ClassNames)
                 || isa (ClassNames, 'categorical')
                 || isa (ClassNames, 'string')))
            error (strcat ("TreeBagger: 'ClassNames' must be a", ...
                           " categorical array, a character array, a", ...
                           " string array, a logical vector, a numeric", ...
                           " vector, or a cell array of character", ...
                           " vectors."));
          endif
          textC = ! (isnumeric (ClassNames) || islogical (ClassNames));
          textU = ! (isnumeric (U) || islogical (U));
          if (textC && textU)
            names = cellstr (ClassNames);
            [tf, loc] = ismember (names(:), cellstr (U));
          elseif (! textC && ! textU)
            [tf, loc] = ismember (ClassNames(:), U(:));
          else
            tf = false;
          endif
          if (! all (tf))
            error ("TreeBagger: not all 'ClassNames' are present in Y.");
          endif
          C = labelsFromIndex (U, loc(:));
        endif
        gY = labelIndices (C, Y);
        keep = gY > 0;
        if (! any (keep))
          error ("TreeBagger: no observations with a known class.");
        endif
        X = X(keep,:);
        Y = Y(keep,:);
        gY = gY(keep);
        RawW = RawW(keep);
        if (! (sum (RawW) > 0))
          error (strcat ("TreeBagger: 'Weights' must not be zero for", ...
                         " every observation used."));
        endif
        K = classCount (C);

        if (isempty (Cost))
          Cost = 1 - eye (K);
        else
          [Cost, errmsg] = costMatrix (Cost, C);
          if (! isempty (errmsg))
            error ("TreeBagger: %s", errmsg);
          endif
        endif

        if (isstruct (Prior))
          P = priorFromStruct (Prior, C, 'TreeBagger');
        elseif (ischar (Prior) && strcmpi (Prior, 'empirical'))
          P = accumarray (gY(:), RawW, [K, 1])';
        elseif (ischar (Prior) && strcmpi (Prior, 'uniform'))
          P = ones (1, K);
        elseif (isnumeric (Prior) && isvector (Prior) && isreal (Prior)
                && numel (Prior) == K && all (Prior >= 0) && sum (Prior) > 0)
          P = double (Prior(:)');
        else
          error (strcat ("TreeBagger: 'Prior' must be 'empirical',", ...
                         " 'uniform', a structure, or a nonnegative", ...
                         " vector with one element per class."));
        endif
        this.Prior = P / sum (P);
        this.Cost = Cost;
        this.ClassNames = C;
        this.W = priorNormalize (RawW, gY, this.Prior);
        this.gY = gY;
        [~, this.DefaultIndex] = max (this.Prior);
        this.DefaultYfit = labelsFromIndex (C, this.DefaultIndex);
        this.DefaultScore = this.Prior;

      else

        if (! (isnumeric (Y) && isreal (Y)))
          error (strcat ("TreeBagger: Y must be a real numeric vector for", ...
                         " a regression ensemble."));
        endif
        Y = double (Y);
        keep = ! isnan (Y);
        if (! any (keep))
          error ("TreeBagger: no observations with a known response.");
        endif
        X = X(keep,:);
        Y = Y(keep);
        RawW = RawW(keep);
        if (! (sum (RawW) > 0))
          error (strcat ("TreeBagger: 'Weights' must not be zero for", ...
                         " every observation used."));
        endif
        this.W = RawW / sum (RawW);
        this.DefaultYfit = sum (this.W .* Y);

      endif

      p = columns (X);
      if (isempty (PredictorNames))
        PredictorNames = arrayfun (@(k) sprintf ('x%d', k), 1:p, ...
                                   'UniformOutput', false);
      endif
      if (isempty (NumVarSample))
        if (isclass)
          NumVarSample = ceil (sqrt (p));
        else
          NumVarSample = ceil (p / 3);
        endif
      endif
      if (isempty (MinLeafSize))
        MinLeafSize = 1 + 4 * ! isclass;
      endif

      ## The pairs every tree is grown with.  The tree options given come
      ## last, so that they override the ensemble's own defaults.
      this.GrowArgs = [{'PredictorNames', PredictorNames, ...
                        'MinLeafSize', MinLeafSize, ...
                        'MinParentSize', 2 * MinLeafSize, ...
                        'MergeLeaves', 'off', 'Prune', 'off', ...
                        'NumVariablesToSample', NumVarSample}, TreeArgs];
      this.MergeLeaves = lastOnOff (TreeArgs, 'mergeleaves');
      this.Prune = lastOnOff (TreeArgs, 'prune');
      if (this.MergeLeaves)
        warning (strcat ("TreeBagger: merging leaves for bagged trees", ...
                         " is not recommended."));
      endif

      this.Method = Method;
      this.X = X;
      this.Y = Y;
      this.PredictorNames = PredictorNames;
      this.NumPredictorsToSample = NumVarSample;
      this.MinLeafSize = MinLeafSize;
      this.InBagFraction = InBagFraction;
      this.SampleWithReplacement = WithReplacement;
      this.ComputeOOBPrediction = OOBPrediction;
      this.TreeArguments = TreeArgs;
      if (OOBPrediction)
        this.OOBIndices = false (rows (X), 0);
      endif

      this = growForest (this, NumTrees, NumPrint);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{label} =} predict (@var{obj}, @var{X})
    ## @deftypefnx {TreeBagger} {[@var{label}, @var{scores}] =} predict (@dots{})
    ## @deftypefnx {TreeBagger} {[@var{label}, @var{scores}, @var{stdevs}] =} predict (@dots{})
    ## @deftypefnx {TreeBagger} {[@var{Yfit}, @var{stdevs}] =} predict (@dots{})
    ## @deftypefnx {TreeBagger} {@dots{} =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict responses with a bagged ensemble.
    ##
    ## Behaves as @code{CompactTreeBagger.predict}, and takes the same
    ## @qcode{'Trees'}, @qcode{'TreeWeights'} and @qcode{'UseInstanceForTree'}
    ## Name-Value arguments.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobPredict, CompactTreeBagger.predict}
    ## @end deftypefn
    function [Yfit, scores, stdevs] = predict (this, X, varargin)

      if (nargin < 2)
        error ("TreeBagger.predict: too few input arguments.");
      endif
      [Yfit, scores, stdevs] = bagPredict (this, X, varargin, ...
                                           'TreeBagger.predict', []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{label} =} oobPredict (@var{obj})
    ## @deftypefnx {TreeBagger} {[@var{label}, @var{scores}, @var{stdevs}] =} oobPredict (@dots{})
    ## @deftypefnx {TreeBagger} {[@var{Yfit}, @var{stdevs}] =} oobPredict (@dots{})
    ## @deftypefnx {TreeBagger} {@dots{} =} oobPredict (@dots{}, 'Trees', @var{trees})
    ##
    ## Out-of-bag predictions for the training data.
    ##
    ## Each observation is predicted by the trees whose samples left it out,
    ## as @code{TreeBagger.predict} would with those trees alone.  An
    ## observation in the sample of every tree used takes @code{DefaultYfit},
    ## with the prior as its scores.  @qcode{'Trees'} restricts the trees.
    ## The ensemble must have been fitted with @qcode{'OOBPrediction'} on.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobError, TreeBagger.predict}
    ## @end deftypefn
    function [Yfit, scores, stdevs] = oobPredict (this, varargin)

      requireOOB (this, 'TreeBagger.oobPredict');
      [Yfit, scores, stdevs] = bagPredict (this, this.X, varargin, ...
                                           'TreeBagger.oobPredict', ...
                                           this.OOBIndices);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{err} =} error (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {TreeBagger} {@var{err} =} error (@dots{}, @var{name}, @var{value})
    ##
    ## Misclassification probability or mean squared error of the ensemble.
    ##
    ## Behaves as @code{CompactTreeBagger.error} and takes the same Name-Value
    ## arguments.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobError, CompactTreeBagger.error}
    ## @end deftypefn
    function err = error (this, X, Y, varargin)

      if (nargin < 3)
        error ("TreeBagger.error: too few input arguments.");
      endif
      err = bagLoss ('error', this, X, Y, varargin, 'TreeBagger.error', ...
                     [], []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{err} =} oobError (@var{obj})
    ## @deftypefnx {TreeBagger} {@var{err} =} oobError (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-bag misclassification probability or mean squared error.
    ##
    ## The error of the out-of-bag predictions on the training data, weighted
    ## by @code{W}.  In @qcode{'individual'} mode each tree is judged on its
    ## own, the observations in its sample taking @code{DefaultYfit}, and
    ## every observation counts in every mode, as MATLAB counts it.
    ## @qcode{'Mode'}, @qcode{'Trees'}, @qcode{'TreeWeights'} and
    ## @qcode{'Weights'} are taken as by @code{CompactTreeBagger.error}.  The
    ## ensemble must have been fitted with @qcode{'OOBPrediction'} on.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobPredict, TreeBagger.oobMeanMargin}
    ## @end deftypefn
    function err = oobError (this, varargin)

      requireOOB (this, 'TreeBagger.oobError');
      err = bagLoss ('error', this, this.X, this.Y, varargin, ...
                     'TreeBagger.oobError', this.OOBIndices, this.W);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {TreeBagger} {@var{m} =} margin (@dots{}, @var{name}, @var{value})
    ##
    ## Classification margin of each observation.
    ##
    ## Behaves as @code{CompactTreeBagger.margin} and takes the same Name-Value
    ## arguments.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobMargin, CompactTreeBagger.margin}
    ## @end deftypefn
    function m = margin (this, X, Y, varargin)

      if (nargin < 3)
        error ("TreeBagger.margin: too few input arguments.");
      endif
      m = bagLoss ('margin', this, X, Y, varargin, 'TreeBagger.margin', ...
                   [], []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{m} =} oobMargin (@var{obj})
    ## @deftypefnx {TreeBagger} {@var{m} =} oobMargin (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-bag classification margin of each training observation.
    ##
    ## The margins of the out-of-bag predictions, with @qcode{'Mode'},
    ## @qcode{'Trees'} and @qcode{'TreeWeights'} taken as by
    ## @code{CompactTreeBagger.margin}.  The ensemble must have been fitted
    ## with @qcode{'OOBPrediction'} on.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobMeanMargin, TreeBagger.margin}
    ## @end deftypefn
    function m = oobMargin (this, varargin)

      requireOOB (this, 'TreeBagger.oobMargin');
      m = bagLoss ('margin', this, this.X, this.Y, varargin, ...
                   'TreeBagger.oobMargin', this.OOBIndices, []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{mm} =} meanMargin (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {TreeBagger} {@var{mm} =} meanMargin (@dots{}, @var{name}, @var{value})
    ##
    ## Weighted mean classification margin.
    ##
    ## Behaves as @code{CompactTreeBagger.meanMargin} and takes the same
    ## Name-Value arguments.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobMeanMargin,
    ## CompactTreeBagger.meanMargin}
    ## @end deftypefn
    function mm = meanMargin (this, X, Y, varargin)

      if (nargin < 3)
        error ("TreeBagger.meanMargin: too few input arguments.");
      endif
      mm = bagLoss ('meanMargin', this, X, Y, varargin, ...
                    'TreeBagger.meanMargin', [], []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{mm} =} oobMeanMargin (@var{obj})
    ## @deftypefnx {TreeBagger} {@var{mm} =} oobMeanMargin (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-bag mean classification margin.
    ##
    ## The mean of @code{TreeBagger.oobMargin}, weighted by @code{W} unless
    ## @qcode{'Weights'} are given.  The ensemble must have been fitted with
    ## @qcode{'OOBPrediction'} on.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobMargin, TreeBagger.oobError}
    ## @end deftypefn
    function mm = oobMeanMargin (this, varargin)

      requireOOB (this, 'TreeBagger.oobMeanMargin');
      mm = bagLoss ('meanMargin', this, this.X, this.Y, varargin, ...
                    'TreeBagger.oobMeanMargin', this.OOBIndices, this.W);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {TreeBagger} {@var{C} =} compact (@var{obj})
    ##
    ## Drop the training data from an ensemble.
    ##
    ## @code{@var{C} = compact (@var{obj})} returns a
    ## @code{CompactTreeBagger} object holding the trees and what prediction
    ## needs.  It predicts new data identically, and has no out-of-bag
    ## methods.
    ##
    ## @seealso{CompactTreeBagger, TreeBagger}
    ## @end deftypefn
    function C = compact (this)

      C = CompactTreeBagger (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{B} =} growTrees (@var{obj}, @var{NumTrees})
    ## @deftypefnx {TreeBagger} {@var{B} =} growTrees (@dots{}, 'NumPrint', @var{n})
    ##
    ## Grow more trees in an ensemble.
    ##
    ## @var{B} is the ensemble with @var{NumTrees} further trees grown exactly
    ## as the first were, their out-of-bag information added.
    ## @qcode{'NumPrint'} is taken as by the @code{TreeBagger} constructor.
    ##
    ## @seealso{TreeBagger, TreeBagger.append}
    ## @end deftypefn
    function this = growTrees (this, NumTrees, varargin)

      if (nargin < 2)
        error ("TreeBagger.growTrees: too few input arguments.");
      endif
      if (! (isnumeric (NumTrees) && isscalar (NumTrees)
             && isreal (NumTrees) && NumTrees >= 1
             && NumTrees == fix (NumTrees)))
        error ("TreeBagger.growTrees: NUMTREES must be a positive integer.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("TreeBagger.growTrees: name-value arguments must", ...
                       " be in pairs."));
      endif
      NumPrint = 0;
      for i = 1:2:numel (varargin)
        name = varargin{i};
        if (ischar (name) && strcmpi (name, 'numprint'))
          NumPrint = varargin{i+1};
          if (! (isnumeric (NumPrint) && isscalar (NumPrint)
                 && isreal (NumPrint) && NumPrint >= 0
                 && NumPrint == fix (NumPrint)))
            error (strcat ("TreeBagger.growTrees: 'NumPrint' must be a", ...
                           " nonnegative integer."));
          endif
        elseif (ischar (name) && strcmpi (name, 'options'))
          error ("TreeBagger.growTrees: 'Options' is not implemented.");
        else
          error (strcat ("TreeBagger.growTrees: invalid parameter name", ...
                         " in optional pair arguments."));
        endif
      endfor
      this = growForest (this, NumTrees, NumPrint);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {TreeBagger} {@var{B} =} append (@var{B1}, @var{B2})
    ##
    ## Add the trees of one ensemble to another.
    ##
    ## @var{B} is @var{B1} with the trees of @var{B2} appended.  The two must
    ## be of the same type, fitted on the same number of observations, with
    ## the same classes and priors, and must agree on whether they keep
    ## out-of-bag information, whose matrices are then joined.
    ##
    ## @seealso{TreeBagger, TreeBagger.growTrees, CompactTreeBagger.combine}
    ## @end deftypefn
    function this = append (this, other)

      if (nargin < 2)
        error ("TreeBagger.append: too few input arguments.");
      endif
      if (! isa (other, 'TreeBagger'))
        error ("TreeBagger.append: B2 must be a TreeBagger object.");
      endif
      if (! strcmp (this.Method, other.Method))
        error (strcat ("TreeBagger.append: the two ensembles must be of", ...
                       " the same type."));
      endif
      if (rows (this.X) != rows (other.X)
          || ! isequal (this.PredictorNames, other.PredictorNames))
        error (strcat ("TreeBagger.append: the two ensembles must be", ...
                       " fitted on the same observations and predictors."));
      endif
      if (! (isequal (this.ClassNames, other.ClassNames)
             && isequal (this.Prior, other.Prior)))
        error (strcat ("TreeBagger.append: the two ensembles must have", ...
                       " the same classes and priors."));
      endif
      if (this.ComputeOOBPrediction != other.ComputeOOBPrediction)
        error (strcat ("TreeBagger.append: the two ensembles have", ...
                       " incompatible out-of-bag flags."));
      endif
      this.Trees = [this.Trees; other.Trees];
      this.TreeClassIdx = [this.TreeClassIdx; other.TreeClassIdx];
      this.NumTrees = numel (this.Trees);
      if (this.ComputeOOBPrediction)
        this.OOBIndices = [this.OOBIndices, other.OOBIndices];
        this.OOBInstanceWeight = sum (this.OOBIndices, 2);
      endif

    endfunction

  endmethods

  methods (Access = private)

    ## Grow n trees and add them to the ensemble.  Each sample is drawn in
    ## proportion to W, with replacement through the cumulative weights and
    ## without it through exponential keys, which leave a zero weight
    ## undrawn in both.  A tree is grown on the classes its sample holds, in
    ## the order the tree sorts them, and TreeClassIdx records where its
    ## scores go among the ensemble's classes.
    function this = growForest (this, n, NumPrint)

      N = rows (this.X);
      m = ceil (this.InBagFraction * N);
      isclass = strcmp (this.Method, 'classification');
      if (this.SampleWithReplacement)
        c = cumsum (this.W);
        c = [0; c / c(end)];
      endif
      for t = 1:n
        if (this.SampleWithReplacement)
          idx = lookup (c, rand (m, 1));
        else
          [~, order] = sort (rand (N, 1) .^ (1 ./ this.W), 'descend');
          idx = order(1:m);
        endif
        if (isclass)
          present = unique (this.gY(idx));
          treeC = uniqueLabels (labelsFromIndex (this.ClassNames, present));
          tidx = labelIndices (this.ClassNames, treeC);
          Mdl = ClassificationTree (this.X(idx,:), this.Y(idx,:), ...
                                    'ClassNames', treeC, ...
                                    'Prior', this.Prior(tidx), ...
                                    'Cost', this.Cost(tidx, tidx), ...
                                    this.GrowArgs{:});
          this.TreeClassIdx{end+1,1} = tidx(:)';
        else
          Mdl = RegressionTree (this.X(idx,:), this.Y(idx), ...
                                this.GrowArgs{:});
          this.TreeClassIdx{end+1,1} = [];
        endif
        this.Trees{end+1,1} = compact (Mdl);
        if (this.ComputeOOBPrediction)
          inbag = false (N, 1);
          inbag(idx) = true;
          this.OOBIndices(:,end+1) = ! inbag;
        endif
        this.NumTrees = numel (this.Trees);
        if (NumPrint > 0 && mod (this.NumTrees, NumPrint) == 0)
          printf ("Tree %d done.\n", this.NumTrees);
        endif
      endfor
      if (this.ComputeOOBPrediction)
        this.OOBInstanceWeight = sum (this.OOBIndices, 2);
      endif

    endfunction

    function requireOOB (this, caller)
      if (! this.ComputeOOBPrediction)
        error (strcat ("%s: out-of-bag information was not kept; fit", ...
                       " with 'OOBPrediction' set to 'on'."), caller);
      endif
    endfunction

  endmethods

endclassdef

## 'on' or 'off' as a logical scalar, OK false for anything else.
function [tf, ok] = onOff (val)

  ok = ischar (val) && any (strcmpi (val, {'on', 'off'}));
  tf = ok && strcmpi (val, 'on');

endfunction

## The last 'on' or 'off' given for a tree option, false when none was.
function tf = lastOnOff (args, name)

  tf = false;
  for i = numel (args) - 1:-2:1
    if (strcmpi (args{i}, name))
      tf = ischar (args{i+1}) && strcmpi (args{i+1}, 'on');
      return;
    endif
  endfor

endfunction

%!demo
%! ## Grow a random forest on the iris data, estimate its error from the
%! ## observations each tree left out, and classify a new flower.
%! load fisheriris
%! rng (42);
%! B = TreeBagger (50, meas, species, 'OOBPrediction', 'on');
%! err = oobError (B);
%! err(end)
%! label = predict (B, [5.8, 2.8, 4.5, 1.4])

%!demo
%! ## A regression forest predicting sepal length from the other three
%! ## measurements, with the spread of its trees around each prediction.
%! load fisheriris
%! rng (42);
%! B = TreeBagger (50, meas(:,2:4), meas(:,1), 'Method', 'regression');
%! [yfit, sd] = predict (B, meas([1, 51, 101], 2:4))

%!test  # MATLAB parity: the defaults of a classification ensemble
%! load fisheriris
%! rng (1);
%! B = TreeBagger (5, meas, species);
%! assert_equal (B.NumTrees, 5);
%! assert_equal (B.Method, 'classification');
%! assert_equal (B.NumPredictorsToSample, 2);
%! assert_equal (B.MinLeafSize, 1);
%! assert_equal (B.InBagFraction, 1);
%! assert_equal (B.SampleWithReplacement, true);
%! assert_equal (B.ComputeOOBPrediction, false);
%! assert_equal (B.MergeLeaves, false);
%! assert_equal (B.Prune, false);
%! assert_equal (B.TreeArguments, {});
%! assert_equal (class (B.Trees{1}), 'CompactClassificationTree');

%!test  # MATLAB parity: the number of predictors sampled by default
%! rng (1);
%! X5 = rand (30, 5);
%! X7 = rand (30, 7);
%! y = [ones(15, 1); 2 * ones(15, 1)];
%! assert_equal (TreeBagger (1, X5, y).NumPredictorsToSample, 3);
%! assert_equal (TreeBagger (1, X7, y).NumPredictorsToSample, 3);
%! r = rand (30, 1);
%! assert_equal (TreeBagger (1, X5, r, 'Method', ...
%!                           'regression').NumPredictorsToSample, 2);
%! assert_equal (TreeBagger (1, X7, r, 'Method', ...
%!                           'regression').NumPredictorsToSample, 3);

%!test  # MATLAB parity: a regression leaf holds five observations by default
%! rng (1);
%! B = TreeBagger (1, rand (30, 3), rand (30, 1), 'Method', 'regression');
%! assert_equal (B.MinLeafSize, 5);
%! assert_equal (class (B.Trees{1}), 'CompactRegressionTree');

%!test  # MATLAB parity: without sampling a tree is the single tree
%! load fisheriris
%! B = TreeBagger (2, meas, species, 'SampleWithReplacement', 'off', ...
%!                 'InBagFraction', 1, 'NumPredictorsToSample', 'all');
%! T = ClassificationTree (meas, species, 'Prune', 'off', ...
%!                         'MergeLeaves', 'off', 'MinParentSize', 2);
%! assert_equal (B.Trees{2}.NumNodes, 17);
%! assert_equal (B.Trees{2}.CutPredictorIndex, T.CutPredictorIndex);
%! assert_equal (isequaln (B.Trees{2}.CutPoint, T.CutPoint), true);

%!test  # MATLAB parity: without sampling a regression tree is the single tree
%! load fisheriris
%! B = TreeBagger (1, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'SampleWithReplacement', 'off', 'InBagFraction', 1, ...
%!                 'NumPredictorsToSample', 'all');
%! T = RegressionTree (meas(:,2:4), meas(:,1), 'Prune', 'off', ...
%!                     'MergeLeaves', 'off', 'MinLeafSize', 5);
%! assert_equal (B.Trees{1}.NumNodes, 47);
%! assert_equal (isequaln (B.Trees{1}.CutPoint, T.CutPoint), true);

%!test  # MATLAB parity: scores average the trees and deviations divide by N
%! load fisheriris
%! rng (1);
%! B = TreeBagger (10, meas, species, 'MinLeafSize', 5);
%! [~, s, sd] = predict (B, meas);
%! S = zeros (150, 3, 10);
%! for t = 1:10
%!   [~, st] = predict (B.Trees{t}, meas);
%!   S(:, B.TreeClassIdx{t}, t) = st;
%! endfor
%! assert_equal (s, mean (S, 3), 1e-14);
%! assert_equal (sd, std (S, 1, 3), 1e-12);

%!test  # MATLAB parity: a regression forest averages its trees
%! load fisheriris
%! rng (1);
%! B = TreeBagger (10, meas(:,2:4), meas(:,1), 'Method', 'regression');
%! [yfit, sd] = predict (B, meas(:,2:4));
%! R = zeros (150, 10);
%! for t = 1:10
%!   R(:,t) = predict (B.Trees{t}, meas(:,2:4));
%! endfor
%! assert_equal (yfit, mean (R, 2), 1e-12);
%! assert_equal (sd, std (R, 1, 2), 1e-12);

%!test  # MATLAB parity: the label is the class of highest score under a cost
%! load fisheriris
%! rng (1);
%! B = TreeBagger (10, meas, species, 'MinLeafSize', 15, ...
%!                 'Cost', [0, 1, 1; 20, 0, 1; 1, 1, 0]);
%! [label, s] = predict (B, meas);
%! [~, k] = max (s, [], 2);
%! assert_equal (label, B.ClassNames(k));

%!test  # labels take the type of the response, where MATLAB gives cellstr
%! load fisheriris
%! rng (1);
%! B = TreeBagger (3, meas, categorical (species));
%! assert_equal (class (predict (B, meas(1,:))), 'categorical');
%! assert_equal (class (B.ClassNames), 'categorical');
%! assert_equal (class (B.DefaultYfit), 'categorical');

%!test  # MATLAB parity: a sample holds ceil (InBagFraction * N) observations
%! load fisheriris
%! rng (1);
%! B = TreeBagger (1, meas, species, 'InBagFraction', 0.334);
%! assert_equal (B.Trees{1}.NodeSize(1), 51);

%!test  # MATLAB parity: the default class is the class of most weight
%! load fisheriris
%! w = [20 * ones(20, 1); ones(100, 1)];
%! B = TreeBagger (1, meas([1:20, 51:150],:), species([1:20, 51:150]), ...
%!                 'Weights', w);
%! assert_equal (B.DefaultYfit, {'setosa'});
%! assert_equal (sum (B.W), 1, 1e-15);

%!test  # MATLAB parity: the regression default is the weighted mean
%! load fisheriris
%! w = [5 * ones(50, 1); ones(100, 1)];
%! B = TreeBagger (1, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'Weights', w);
%! assert_equal (B.DefaultYfit, sum (w .* meas(:,1)) / sum (w), 1e-12);

%!test  # MATLAB parity: a missing response drops its row, a missing X does not
%! load fisheriris
%! y = meas(:,1);
%! y([3, 7]) = NaN;
%! X = meas(:,2:4);
%! X(5, 1) = NaN;
%! B = TreeBagger (1, X, y, 'Method', 'regression', 'OOBPrediction', 'on');
%! assert_equal (size (B.X), [148, 3]);
%! assert_equal (size (B.OOBIndices), [148, 1]);

%!test  # MATLAB parity: 'ClassNames' keeps its order and drops other classes
%! load fisheriris
%! B = TreeBagger (1, meas, species, 'ClassNames', {'virginica'; 'setosa'});
%! assert_equal (B.ClassNames, {'virginica'; 'setosa'});
%! assert_equal (rows (B.X), 100);

%!test  # MATLAB parity: the out-of-bag error in its three modes
%! load fisheriris
%! rng (1);
%! B = TreeBagger (8, meas, species, 'OOBPrediction', 'on', 'MinLeafSize', 5);
%! e = oobError (B);
%! assert_equal (size (e), [8, 1]);
%! assert_equal (size (oobError (B, 'Mode', 'individual')), [8, 1]);
%! assert_equal (e(end), oobError (B, 'Mode', 'ensemble'), 1e-15);
%! assert_equal (B.OOBInstanceWeight, sum (B.OOBIndices, 2));

%!test  # MATLAB parity: a tree judged alone gives its sample the default
%! load fisheriris
%! rng (1);
%! B = TreeBagger (4, meas, species, 'OOBPrediction', 'on', 'MinLeafSize', 5);
%! l = predict (B.Trees{3}, meas);
%! l(! B.OOBIndices(:,3)) = B.DefaultYfit;
%! e = oobError (B, 'Mode', 'individual');
%! assert_equal (e(3), sum (B.W .* ! strcmp (l, species)), 1e-14);

%!test  # MATLAB parity: an observation in every sample takes the default
%! load fisheriris
%! rng (5);
%! B = TreeBagger (1, meas, species, 'OOBPrediction', 'on', 'MinLeafSize', 15);
%! r = find (! B.OOBIndices, 1);
%! [label, s] = oobPredict (B);
%! assert_equal (label(r), B.DefaultYfit);
%! assert_equal (s(r,:), B.Prior, 1e-15);

%!test  # MATLAB parity: margins by tree and their means
%! load fisheriris
%! rng (1);
%! B = TreeBagger (6, meas, species, 'OOBPrediction', 'on');
%! assert_equal (size (margin (B, meas, species)), [150, 6]);
%! assert_equal (size (meanMargin (B, meas, species)), [1, 6]);
%! assert_equal (size (oobMargin (B)), [150, 6]);
%! assert_equal (size (oobMeanMargin (B, 'Mode', 'ensemble')), [1, 1]);

%!test  # MATLAB parity: weights make the error a weighted share
%! load fisheriris
%! rng (1);
%! B = TreeBagger (6, meas, species, 'MinLeafSize', 15);
%! w = [5 * ones(50, 1); ones(100, 1)];
%! miss = ! strcmp (predict (B, meas), species);
%! e = error (B, meas, species, 'Mode', 'ensemble', 'Weights', w);
%! assert_equal (e, sum (w .* miss) / sum (w), 1e-15);

%!test  # MATLAB parity: an observation no tree may use takes the default
%! load fisheriris
%! rng (1);
%! B = TreeBagger (4, meas, species);
%! U = true (2, 4);
%! U(1,:) = false;
%! [label, s, sd] = predict (B, meas(1:2,:), 'UseInstanceForTree', U);
%! assert_equal (label(1), B.DefaultYfit);
%! assert_equal (s(1,:), B.Prior, 1e-15);
%! assert_equal (sd(1,:), NaN (1, 3));

%!test  # tree weights give a weighted average of the trees chosen
%! load fisheriris
%! rng (1);
%! B = TreeBagger (6, meas, species, 'MinLeafSize', 15);
%! [~, s] = predict (B, meas(1:5,:), 'Trees', [2, 5], 'TreeWeights', [3, 1]);
%! [~, a] = predict (B.Trees{2}, meas(1:5,:));
%! [~, b] = predict (B.Trees{5}, meas(1:5,:));
%! A = zeros (5, 3);
%! A(:, B.TreeClassIdx{2}) = a;
%! Bs = zeros (5, 3);
%! Bs(:, B.TreeClassIdx{5}) = b;
%! assert_equal (s, (3 * A + Bs) / 4, 1e-15);

%!test  # the state of the generator reproduces an ensemble
%! load fisheriris
%! rng (3);
%! A = TreeBagger (3, meas, species, 'OOBPrediction', 'on');
%! rng (3);
%! B = TreeBagger (3, meas, species, 'OOBPrediction', 'on');
%! assert_equal (A.OOBIndices, B.OOBIndices);
%! assert_equal (A.Trees{3}.CutPredictorIndex, B.Trees{3}.CutPredictorIndex);

%!test  # MATLAB parity: growing and appending extend the out-of-bag matrix
%! load fisheriris
%! rng (1);
%! B = TreeBagger (3, meas, species, 'OOBPrediction', 'on');
%! B = growTrees (B, 2);
%! assert_equal (B.NumTrees, 5);
%! assert_equal (size (B.OOBIndices), [150, 5]);
%! B = append (B, TreeBagger (2, meas, species, 'OOBPrediction', 'on'));
%! assert_equal (B.NumTrees, 7);
%! assert_equal (size (B.OOBIndices), [150, 7]);

%!test  # MATLAB parity: the tree options given are kept as given
%! load fisheriris
%! B = TreeBagger (1, meas, species, 'MinLeafSize', 3, 'MaxNumSplits', 7);
%! assert_equal (B.TreeArguments, {'MinLeafSize', 3, 'MaxNumSplits', 7});

%!test  # MATLAB parity: progress is printed after every NumPrint trees
%! load fisheriris
%! out = evalc ("TreeBagger (2, meas, species, 'NumPrint', 1);");
%! assert_equal (out, sprintf ("Tree 1 done.\nTree 2 done.\n"));

%!warning<TreeBagger: merging leaves for bagged trees is not recommended.> ...
%! TreeBagger (1, [1; 2; 3; 4], [1; 1; 2; 2], 'MergeLeaves', 'on');

## Test input validation
%!shared x, y, B, C, R
%! load fisheriris
%! x = meas;
%! y = species;
%! B = TreeBagger (3, x, y, 'OOBPrediction', 'on');
%! C = compact (B);
%! R = TreeBagger (2, x(:,2:4), x(:,1), 'Method', 'regression');
%!error<TreeBagger: too few input arguments.> TreeBagger (1, x)
%!error<TreeBagger: NUMTREES must be a positive integer.> TreeBagger (0, x, y)
%!error<TreeBagger: X must be a non-empty real numeric matrix.> ...
%! TreeBagger (1, {1}, y)
%!error<TreeBagger: Y must be a vector with one element per row in X.> ...
%! TreeBagger (1, x, y(1:10))
%!error<TreeBagger: name-value arguments must be in pairs.> ...
%! TreeBagger (1, x, y, 'Method')
%!error<TreeBagger: invalid parameter name in optional pair arguments.> ...
%! TreeBagger (1, x, y, 1, 2)
%!error<TreeBagger: invalid parameter name in optional pair arguments.> ...
%! TreeBagger (1, x, y, 'Foo', 2)
%!error<TreeBagger: 'Method' must be 'classification' or 'regression'.> ...
%! TreeBagger (1, x, y, 'Method', 'cluster')
%!error<TreeBagger: 'NumPredictorsToSample' must be a positive integer or 'all'.> ...
%! TreeBagger (1, x, y, 'NumPredictorsToSample', 0)
%!error<TreeBagger: 'MinLeafSize' must be a positive integer.> ...
%! TreeBagger (1, x, y, 'MinLeafSize', 0)
%!error<TreeBagger: 'InBagFraction' must be a number greater than 0 and no greater than 1.> ...
%! TreeBagger (1, x, y, 'InBagFraction', 1.5)
%!error<TreeBagger: 'SampleWithReplacement' must be 'on' or 'off'.> ...
%! TreeBagger (1, x, y, 'SampleWithReplacement', 1)
%!error<TreeBagger: 'OOBPrediction' must be 'on' or 'off'.> ...
%! TreeBagger (1, x, y, 'OOBPrediction', true)
%!error<TreeBagger: 'OOBPredictorImportance' must be 'on' or 'off'.> ...
%! TreeBagger (1, x, y, 'OOBPredictorImportance', 1)
%!error<TreeBagger: 'OOBPredictorImportance' is not implemented.> ...
%! TreeBagger (1, x, y, 'OOBPredictorImportance', 'on')
%!error<TreeBagger: 'PredictorNames' must be a cell array of character vectors with one element per column of X.> ...
%! TreeBagger (1, x, y, 'PredictorNames', {'a'})
%!error<TreeBagger: 'CategoricalPredictors' is not implemented.> ...
%! TreeBagger (1, x, y, 'CategoricalPredictors', 1)
%!error<TreeBagger: 'NumPrint' must be a nonnegative integer.> ...
%! TreeBagger (1, x, y, 'NumPrint', -1)
%!error<TreeBagger: 'Surrogate' is not implemented.> ...
%! TreeBagger (1, x, y, 'Surrogate', 'on')
%!error<TreeBagger: 'Options' is not implemented.> ...
%! TreeBagger (1, x, y, 'Options', struct ())
%!error<TreeBagger: 'QuadraticErrorTolerance' applies only to regression.> ...
%! TreeBagger (1, x, y, 'QuadraticErrorTolerance', 1e-3)
%!error<TreeBagger: 'Cost' applies only to classification.> ...
%! TreeBagger (1, x(:,2:4), x(:,1), 'Method', 'regression', ...
%!             'Cost', [0, 1; 1, 0])
%!error<TreeBagger: sampling without replacement at an 'InBagFraction' of 1 leaves no observation out of bag for 'OOBPrediction'.> ...
%! TreeBagger (1, x, y, 'SampleWithReplacement', 'off', 'OOBPrediction', 'on')
%!error<TreeBagger: 'Weights' must be a nonnegative numeric vector with one element per row in X.> ...
%! TreeBagger (1, x, y, 'Weights', -ones (150, 1))
%!error<TreeBagger: 'ClassNames' must be a categorical array, a character array, a string array, a logical vector, a numeric vector, or a cell array of character vectors.> ...
%! TreeBagger (1, x, y, 'ClassNames', {1})
%!error<TreeBagger: not all 'ClassNames' are present in Y.> ...
%! TreeBagger (1, x, y, 'ClassNames', {'rose'})
%!error<TreeBagger: no observations with a known class.> ...
%! TreeBagger (1, x(1:3,:), NaN (3, 1))
%!error<TreeBagger: 'Weights' must not be zero for every observation used.> ...
%! TreeBagger (1, x, y, 'Weights', zeros (150, 1))
%!error<TreeBagger: the number of rows and columns in 'Cost' must correspond to selected classes in Y.> ...
%! TreeBagger (1, x, y, 'Cost', [0, 1; 1, 0])
%!error<TreeBagger: 'Prior' must be 'empirical', 'uniform', a structure, or a nonnegative vector with one element per class.> ...
%! TreeBagger (1, x, y, 'Prior', [1, 2])
%!error<TreeBagger: a structure 'Prior' must have 'ClassNames' and 'ClassProbs' fields.> ...
%! TreeBagger (1, x, y, 'Prior', struct ('a', 1))
%!error<TreeBagger: Y must be a real numeric vector for a regression ensemble.> ...
%! TreeBagger (1, x, y, 'Method', 'regression')
%!error<TreeBagger: no observations with a known response.> ...
%! TreeBagger (1, x(1:3,:), NaN (3, 1), 'Method', 'regression')
%!error<TreeBagger.oobError: out-of-bag information was not kept; fit with 'OOBPrediction' set to 'on'.> ...
%! oobError (TreeBagger (1, x, y))
%!error<TreeBagger.oobPredict: out-of-bag information was not kept; fit with 'OOBPrediction' set to 'on'.> ...
%! oobPredict (TreeBagger (1, x, y))
%!error<TreeBagger.oobMargin: out-of-bag information was not kept; fit with 'OOBPrediction' set to 'on'.> ...
%! oobMargin (TreeBagger (1, x, y))
%!error<TreeBagger.oobMeanMargin: out-of-bag information was not kept; fit with 'OOBPrediction' set to 'on'.> ...
%! oobMeanMargin (TreeBagger (1, x, y))
%!error<TreeBagger.predict: too few input arguments.> predict (B)
%!error<TreeBagger.predict: X must be a real numeric matrix.> predict (B, {1})
%!error<TreeBagger.predict: X must have one column per predictor.> ...
%! predict (B, ones (2, 3))
%!error<TreeBagger.predict: name-value arguments must be in pairs.> ...
%! predict (B, x, 'Trees')
%!error<TreeBagger.predict: invalid parameter name in optional pair arguments.> ...
%! predict (B, x, 'Mode', 'ensemble')
%!error<TreeBagger.predict: 'Trees' must be 'all' or a vector of indices of trees in the ensemble.> ...
%! predict (B, x, 'Trees', 4)
%!error<TreeBagger.predict: 'TreeWeights' must be a nonnegative numeric vector with one element per tree used, not all zero.> ...
%! predict (B, x, 'TreeWeights', [1, 1])
%!error<TreeBagger.predict: 'UseInstanceForTree' must be a logical matrix with one row per observation and one column per tree.> ...
%! predict (B, x, 'UseInstanceForTree', true (2, 3))
%!error<TreeBagger.oobPredict: invalid parameter name in optional pair arguments.> ...
%! oobPredict (B, 'UseInstanceForTree', true (150, 3))
%!error<TreeBagger.margin: too few input arguments.> margin (B, x)
%!error<TreeBagger.margin: invalid parameter name in optional pair arguments.> ...
%! margin (B, x, y, 'Weights', ones (150, 1))
%!error<TreeBagger.meanMargin: too few input arguments.> meanMargin (B, x)
%!error<TreeBagger.margin: margins are defined only for classification.> ...
%! margin (R, x(:,2:4), x(:,1))
## Octave's test.m cuts a message through its first 'error:' before matching
## it, so these patterns hold the message after the method's prefix.
%!error<too few input arguments.> error (B, x)
%!error<X and Y must have the same number of rows.> error (B, x, y(1:3))
%!error<'Weights' must be a nonnegative numeric vector with one element per observation, not all zero.> ...
%! error (B, x, y, 'Weights', ones (3, 1))
%!error<'Mode' must be 'cumulative', 'individual' or 'ensemble'.> ...
%! error (B, x, y, 'Mode', 'all')
%!error<'TreeWeights' cannot be used in 'individual' mode.> ...
%! error (B, x, y, 'Mode', 'individual', 'TreeWeights', [1, 1, 1])
%!error<Y must hold only classes the model was trained on.> ...
%! error (B, x(1:2,:), {'rose'; 'setosa'})
%!error<Y must be a real numeric vector.> error (R, x(1:2,2:4), {'a'; 'b'})
%!error<TreeBagger.growTrees: too few input arguments.> growTrees (B)
%!error<TreeBagger.growTrees: NUMTREES must be a positive integer.> ...
%! growTrees (B, 0)
%!error<TreeBagger.growTrees: name-value arguments must be in pairs.> ...
%! growTrees (B, 1, 'NumPrint')
%!error<TreeBagger.growTrees: 'NumPrint' must be a nonnegative integer.> ...
%! growTrees (B, 1, 'NumPrint', -1)
%!error<TreeBagger.growTrees: 'Options' is not implemented.> ...
%! growTrees (B, 1, 'Options', 1)
%!error<TreeBagger.growTrees: invalid parameter name in optional pair arguments.> ...
%! growTrees (B, 1, 'Foo', 1)
%!error<TreeBagger.append: too few input arguments.> append (B)
%!error<TreeBagger.append: B2 must be a TreeBagger object.> append (B, C)
%!error<TreeBagger.append: the two ensembles must be of the same type.> ...
%! append (B, TreeBagger (1, x, x(:,1), 'Method', 'regression'))
%!error<TreeBagger.append: the two ensembles must be fitted on the same observations and predictors.> ...
%! append (B, TreeBagger (1, x(1:100,:), y(1:100), 'OOBPrediction', 'on'))
%!error<TreeBagger.append: the two ensembles must have the same classes and priors.> ...
%! append (B, TreeBagger (1, x, y, 'OOBPrediction', 'on', ...
%!                     'Prior', 'uniform', ...
%!                     'ClassNames', {'virginica'; 'versicolor'; 'setosa'}))
%!error<TreeBagger.append: the two ensembles have incompatible out-of-bag flags.> ...
%! append (B, TreeBagger (1, x, y))
