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

classdef TreeBagger < PredictiveModel
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
    ## Where the model was fitted from a table, the predictors are the coded
    ## matrix and not the table: a variable holding levels is stored as its
    ## level codes, and the coding is kept with the model.
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
    ## @deftp {TreeBagger} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A row vector of column indices into @var{X}, naming the predictors
    ## treated as categorical, empty when none is.  This property is
    ## read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

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
    ## A logical scalar, false unless the ensemble was fitted with
    ## @qcode{'OOBPredictorImportance'} set to @qcode{'on'}, which keeps the
    ## out-of-bag information too.  This property is read-only.
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
    ## @deftp {TreeBagger} {property} DeltaCriterionDecisionSplit
    ##
    ## Split criterion contributions of the predictors
    ##
    ## A row vector with one element per predictor, the mean over the trees
    ## of each tree's @code{predictorImportance}.  This property is read-only.
    ##
    ## @end deftp
    DeltaCriterionDecisionSplit = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} NumPredictorSplit
    ##
    ## Decision splits on each predictor
    ##
    ## A row vector with one element per predictor, the sum over the trees of
    ## the share of each tree's branch nodes that split on the predictor.  A
    ## tree without branch nodes adds nothing.  This property is read-only.
    ##
    ## @end deftp
    NumPredictorSplit = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} SurrogateAssociation
    ##
    ## Predictive association between the predictors
    ##
    ## A square matrix with one row and one column per predictor.  The trees
    ## grow no surrogate splits, so it is the identity matrix.  This property
    ## is read-only.
    ##
    ## @end deftp
    SurrogateAssociation = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} OOBPermutedPredictorDeltaError
    ##
    ## Rise in out-of-bag error when a predictor is permuted
    ##
    ## A row vector with one element per predictor.  For each tree, the
    ## values of the predictor are permuted among the observations out of its
    ## bag, and the tree's error on them, the misclassification share or the
    ## mean squared error weighted by @code{W}, is taken before and after.
    ## The element is the mean of the rise over the trees divided by its
    ## standard deviation over the trees, zero when the mean is zero.  Reading
    ## it is an error unless @code{ComputeOOBPredictorImportance} is true.
    ## This property is read-only.
    ##
    ## @end deftp
    OOBPermutedPredictorDeltaError = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} OOBPermutedPredictorDeltaMeanMargin
    ##
    ## Fall in out-of-bag mean margin when a predictor is permuted
    ##
    ## A row vector with one element per predictor, computed as
    ## @code{OOBPermutedPredictorDeltaError} is from each tree's weighted mean
    ## classification margin on its out-of-bag observations, before the
    ## permutation less after it.  Empty for a regression ensemble.  Reading
    ## it is an error unless @code{ComputeOOBPredictorImportance} is true.
    ## This property is read-only.
    ##
    ## @end deftp
    OOBPermutedPredictorDeltaMeanMargin = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} OOBPermutedPredictorCountRaiseMargin
    ##
    ## Margins lowered less margins raised when a predictor is permuted
    ##
    ## A row vector with one element per predictor, computed as
    ## @code{OOBPermutedPredictorDeltaError} is from the number of each tree's
    ## out-of-bag observations whose margin the permutation lowers less the
    ## number whose margin it raises.  Empty for a regression ensemble.
    ## Reading it is an error unless @code{ComputeOOBPredictorImportance} is
    ## true.  This property is read-only.
    ##
    ## @end deftp
    OOBPermutedPredictorCountRaiseMargin = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} Proximity
    ##
    ## Proximity of the training observations
    ##
    ## A symmetric @math{NxN} matrix whose element @math{(i,j)} is the share
    ## of the trees that bring training observations @math{i} and @math{j} to
    ## the same leaf.  Reading it is an error until @code{fillprox} fills it,
    ## and @code{growTrees} and @code{append} empty it again.  This property
    ## is read-only.
    ##
    ## @end deftp
    Proximity = [];

    ## -*- texinfo -*-
    ## @deftp {TreeBagger} {property} OutlierMeasure
    ##
    ## Outlier measure of the training observations
    ##
    ## A column with one element per observation, computed from
    ## @code{Proximity} as @code{CompactTreeBagger.outlierMeasure} computes
    ## it, within each class for classification and over every observation
    ## for regression.  Reading it is an error until @code{fillprox} fills it,
    ## and @code{growTrees} and @code{append} empty it again.  This property
    ## is read-only.
    ##
    ## @end deftp
    OutlierMeasure = [];

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
    ResponseName = 'Y';  # name the table gave the response, which
                         # a table call looks up; hidden as MATLAB's
                         # TreeBagger carries no ResponseName
    TreeClassIdx = {};   # columns of ClassNames each tree's scores fill
    DefaultIndex = 0;    # index of DefaultYfit into ClassNames
    DefaultScore = [];   # scores of an observation no tree may answer for
    GrowArgs = {};       # Name-Value pairs every tree is grown with
    gY = [];             # class index of each observation
    InBag = [];          # sparse NxNumTrees count of each observation in bag
    PermDelta = [];      # NumTreesxPx3 per-tree permutation deltas
  endproperties

  methods

    function v = get.OOBPermutedPredictorDeltaError (this)
      requireImportance (this);
      v = this.OOBPermutedPredictorDeltaError;
    endfunction

    function v = get.OOBPermutedPredictorDeltaMeanMargin (this)
      requireImportance (this);
      v = this.OOBPermutedPredictorDeltaMeanMargin;
    endfunction

    function v = get.OOBPermutedPredictorCountRaiseMargin (this)
      requireImportance (this);
      v = this.OOBPermutedPredictorCountRaiseMargin;
    endfunction

    function v = get.Proximity (this)
      requireProximity (this, this.Proximity);
      v = this.Proximity;
    endfunction

    function v = get.OutlierMeasure (this)
      requireProximity (this, this.OutlierMeasure);
      v = this.OutlierMeasure;
    endfunction

  endmethods

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
    ## @deftypefnx {TreeBagger} {@var{B} =} TreeBagger (@var{NumTrees}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {TreeBagger} {@var{B} =} TreeBagger (@var{NumTrees}, @var{Tbl}, @var{formula})
    ## @deftypefnx {TreeBagger} {@var{B} =} TreeBagger (@var{NumTrees}, @var{Tbl}, @var{Y})
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
    ## @item @qcode{'OOBPredictorImportance'} @tab @tab @qcode{'off'} (default)
    ## or @qcode{'on'}, to estimate the importance of each predictor by
    ## permuting it among each tree's out-of-bag observations.  It turns
    ## @qcode{'OOBPrediction'} on.
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
    ## @qcode{'CategoricalPredictors'}, @qcode{'MaxNumCategories'},
    ## @qcode{'MaxNumSplits'}, @qcode{'MergeLeaves'}, @qcode{'Prune'},
    ## @qcode{'PruneCriterion'} and @qcode{'SplitCriterion'} are passed on to
    ## every tree, and so are @qcode{'AlgorithmForCategorical'} for
    ## classification and @qcode{'QuadraticErrorTolerance'} for regression;
    ## see @code{fitctree} and @code{fitrtree}.  Merging leaves is allowed but
    ## warned against.
    ##
    ## Surrogate splits, parallel growth and tall arrays are not implemented,
    ## and an option asking for one of them is refused.
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

      ## A table names its own predictors and says which hold levels.  This
      ## class carries no public ResponseName, as MATLAB's does not either,
      ## so the name the table gave the response is kept in the hidden one
      ## rather than passed on to the trees.
      [this, X, Y, varargin] = resolveTable (this, 'TreeBagger', X, Y, ...
                                             varargin);
      for k = numel (varargin) - 1:-1:1
        if (ischar (varargin{k}) && strcmpi (varargin{k}, 'ResponseName'))
          this.ResponseName = varargin{k+1};
          varargin(k:k+1) = [];
        endif
      endfor
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
      Importance = false;
      Cost = [];
      Prior = 'empirical';
      ClassNames = [];
      Weights = [];
      PredictorNames = {};
      CatSpec = [];
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
            [Importance, ok] = onOff (Value);
            if (! ok)
              error (strcat ("TreeBagger: 'OOBPredictorImportance' must", ...
                             " be 'on' or 'off'."));
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
            CatSpec = Value;

          case 'numprint'
            if (! (isnumeric (Value) && isscalar (Value) && isreal (Value)
                   && Value >= 0 && Value == fix (Value)))
              error (strcat ("TreeBagger: 'NumPrint' must be a", ...
                             " nonnegative integer."));
            endif
            NumPrint = double (Value);

          case {'maxnumsplits', 'mergeleaves', 'prune', 'prunecriterion', ...
                'splitcriterion', 'maxnumcategories', ...
                'algorithmforcategorical'}
            TreeArgs(end+1:end+2) = {name, Value};

          case 'quadraticerrortolerance'
            TreeArgs(end+1:end+2) = {name, Value};
            givenReg{end+1} = name;

          case {'surrogate', 'predictorselection', 'options', 'chunksize'}
            error ("TreeBagger: '%s' is not implemented.", name);

          otherwise
            error (strcat ("TreeBagger: invalid parameter name in", ...
                           " optional pair arguments."));

        endswitch
      endfor

      isclass = strcmp (Method, 'classification');
      OOBPrediction = OOBPrediction || Importance;
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

      p = columns (X);
      if (isempty (PredictorNames))
        PredictorNames = arrayfun (@(k) sprintf ('x%d', k), 1:p, ...
                                   'UniformOutput', false);
      endif

      ## The categorical predictors are resolved once the predictor names are
      ## known, so that they may be named as well as indexed, and before any
      ## row is left out, so that a level is not lost with the row carrying it.
      [Cod, errmsg] = dummyCoding (X, CatSpec, PredictorNames);
      if (! isempty (errmsg))
        error ("TreeBagger: %s", errmsg);
      endif
      CatPred = [];
      if (! isempty (Cod.Index))
        CatPred = Cod.Index;
        TreeArgs(end+1:end+2) = {'CategoricalPredictors', CatPred};
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
          [loc, errmsg] = namedClasses (U, ClassNames);
          if (! isempty (errmsg))
            error ("TreeBagger: %s", errmsg);
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
      this.CategoricalPredictors = CatPred;
      this.NumPredictorsToSample = NumVarSample;
      this.MinLeafSize = MinLeafSize;
      this.InBagFraction = InBagFraction;
      this.SampleWithReplacement = WithReplacement;
      this.ComputeOOBPrediction = OOBPrediction;
      this.ComputeOOBPredictorImportance = Importance;
      this.TreeArguments = TreeArgs;
      this.InBag = sparse (rows (X), 0);
      if (OOBPrediction)
        this.OOBIndices = false (rows (X), 0);
      endif
      if (Importance)
        this.PermDelta = zeros (0, p, 3);
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
    ##
    ## The new data may be a table, whose variables are matched to the
    ## predictors the model was fitted on by name and not by position:
    ## one the model was not fitted on is passed over, one it needs and
    ## cannot find is named, and a value holding a level is coded as that
    ## level was coded at fitting.
    ## @seealso{TreeBagger, TreeBagger.oobPredict, CompactTreeBagger.predict}
    ## @end deftypefn
    function [Yfit, scores, stdevs] = predict (this, X, varargin)

      if (nargin < 2)
        error ("TreeBagger.predict: too few input arguments.");
      endif

      ## A table is read by the names the model was fitted on
      X = tableColumns (this, 'TreeBagger.predict', X);
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
    ## @deftypefnx {TreeBagger} {@var{err} =} error (@var{obj}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {TreeBagger} {@var{err} =} error (@var{obj}, @var{Tbl})
    ##
    ## Misclassification probability or mean squared error of the ensemble.
    ##
    ## Behaves as @code{CompactTreeBagger.error} and takes the same Name-Value
    ## arguments.
    ##
    ## @var{X} may also be a table @var{Tbl}, whose variables are matched to
    ## the predictors the model was fitted on by name and not by position.
    ## @code{error (@var{obj}, @var{Tbl}, @var{ResponseVarName})} takes the
    ## response from the variable @var{ResponseVarName} names, and
    ## @code{error (@var{obj}, @var{Tbl})} from the variable the model was
    ## fitted on.  The response may also be given beside the table as
    ## @var{Y}.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobError, CompactTreeBagger.error}
    ## @end deftypefn
    function err = error (this, X, Y, varargin)

      if (nargin < 3 && ! (nargin > 1 && istable (X)))
        error ("TreeBagger.error: too few input arguments.");
      endif

      ## A table carries the response: named in the call, given beside
      ## the table, or the variable the model was fitted on
      if (nargin < 3)
        Y = [];
      endif
      [X, Y, varargin] = tableResponse (this, 'error', X, Y, ...
                                        varargin, nargin > 2);
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
    ## @deftypefnx {TreeBagger} {@var{m} =} margin (@var{obj}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {TreeBagger} {@var{m} =} margin (@var{obj}, @var{Tbl})
    ##
    ## Classification margin of each observation.
    ##
    ## Behaves as @code{CompactTreeBagger.margin} and takes the same Name-Value
    ## arguments.
    ##
    ## @var{X} may also be a table @var{Tbl}, whose variables are matched to
    ## the predictors the model was fitted on by name and not by position.
    ## @code{margin (@var{obj}, @var{Tbl}, @var{ResponseVarName})} takes the
    ## response from the variable @var{ResponseVarName} names, and
    ## @code{margin (@var{obj}, @var{Tbl})} from the variable the model was
    ## fitted on.  The response may also be given beside the table as
    ## @var{Y}.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobMargin, CompactTreeBagger.margin}
    ## @end deftypefn
    function m = margin (this, X, Y, varargin)

      if (nargin < 3 && ! (nargin > 1 && istable (X)))
        error ("TreeBagger.margin: too few input arguments.");
      endif

      ## A table carries the response: named in the call, given beside
      ## the table, or the variable the model was fitted on
      if (nargin < 3)
        Y = [];
      endif
      [X, Y, varargin] = tableResponse (this, 'margin', X, Y, ...
                                        varargin, nargin > 2);
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
    ## @deftypefnx {TreeBagger} {@var{mm} =} meanMargin (@var{obj}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {TreeBagger} {@var{mm} =} meanMargin (@var{obj}, @var{Tbl})
    ##
    ## Weighted mean classification margin.
    ##
    ## Behaves as @code{CompactTreeBagger.meanMargin} and takes the same
    ## Name-Value arguments.
    ##
    ## @var{X} may also be a table @var{Tbl}, whose variables are matched to
    ## the predictors the model was fitted on by name and not by position.
    ## @code{meanMargin (@var{obj}, @var{Tbl}, @var{ResponseVarName})}
    ## takes the response from the variable @var{ResponseVarName} names, and
    ## @code{meanMargin (@var{obj}, @var{Tbl})} from the variable the model
    ## was fitted on.  The response may also be given beside the table as
    ## @var{Y}.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobMeanMargin,
    ## CompactTreeBagger.meanMargin}
    ## @end deftypefn
    function mm = meanMargin (this, X, Y, varargin)

      if (nargin < 3 && ! (nargin > 1 && istable (X)))
        error ("TreeBagger.meanMargin: too few input arguments.");
      endif

      ## A table carries the response: named in the call, given beside
      ## the table, or the variable the model was fitted on
      if (nargin < 3)
        Y = [];
      endif
      [X, Y, varargin] = tableResponse (this, 'meanMargin', X, Y, ...
                                        varargin, nargin > 2);
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
    ## @deftypefn  {TreeBagger} {@var{YFit} =} quantilePredict (@var{obj}, @var{X})
    ## @deftypefnx {TreeBagger} {[@var{YFit}, @var{YW}] =} quantilePredict (@dots{})
    ## @deftypefnx {TreeBagger} {@dots{} =} quantilePredict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict quantiles of the response with a regression ensemble.
    ##
    ## Each tree gives every training observation a weight for each row of
    ## @var{X}: the number of times the tree's sample holds the observation,
    ## divided by the size of the leaf the row comes to rest at, when the
    ## observation is in that leaf, and zero otherwise.  The weights are
    ## averaged over the trees, with their tree weights, into @var{YW}, a
    ## sparse @math{NxM} matrix with one row per training observation and one
    ## column per row of @var{X}, each column summing to one.
    ##
    ## @var{YFit} is an @math{MxQ} matrix holding, for each row of @var{X} and
    ## each quantile probability, the quantile of the training responses under
    ## those weights.  The responses are sorted, each keeping its own weight,
    ## and the quantile is interpolated linearly between them at their
    ## cumulative weights less half their own weight, taking the smallest or
    ## the largest response beyond either end.  A row no tree may answer for
    ## takes @code{quantile (@var{obj}.Y, @var{tau})}, and its column of
    ## @var{YW} weighs every training observation equally.  The observation
    ## weights enter only through the samples they drew.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Quantile'} @tab @tab A vector of probabilities @var{tau}
    ## between 0 and 1.  The default is 0.5, the median.
    ## @item @qcode{'Trees'} @tab @tab @qcode{'all'} (default) or a vector of
    ## indices of the trees to use.
    ## @item @qcode{'TreeWeights'} @tab @tab A nonnegative vector with one
    ## weight per tree used.  The default weighs them equally.
    ## @item @qcode{'UseInstanceForTree'} @tab @tab An @math{MxNumTrees}
    ## logical matrix saying which tree may answer for which row.
    ## @end multitable
    ##
    ## @seealso{TreeBagger, TreeBagger.oobQuantilePredict,
    ## TreeBagger.quantileError, TreeBagger.predict}
    ## @end deftypefn
    function [YFit, YW] = quantilePredict (this, X, varargin)

      if (nargin < 2)
        error ("TreeBagger.quantilePredict: too few input arguments.");
      endif
      [tau, o] = quantileArgs (this, X, varargin, ...
                               {'Trees', 'TreeWeights', ...
                                'UseInstanceForTree'}, ...
                               'TreeBagger.quantilePredict', false);
      [YFit, YW] = quantileSteps (this, X, o.trees, o.tw, o.use, tau, ...
                                  'ensemble');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{YFit} =} oobQuantilePredict (@var{obj})
    ## @deftypefnx {TreeBagger} {[@var{YFit}, @var{YW}] =} oobQuantilePredict (@dots{})
    ## @deftypefnx {TreeBagger} {@dots{} =} oobQuantilePredict (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-bag quantile predictions for the training data.
    ##
    ## Each training observation is predicted as by
    ## @code{TreeBagger.quantilePredict}, by the trees whose samples left it
    ## out.  An observation in the sample of every tree used takes
    ## @code{quantile (@var{obj}.Y, @var{tau})}.  @var{YW} is @math{NxN}.
    ## @qcode{'Quantile'}, @qcode{'Trees'} and @qcode{'TreeWeights'} are taken
    ## as by @code{TreeBagger.quantilePredict}.  The ensemble must have been
    ## fitted with @qcode{'OOBPrediction'} on.
    ##
    ## @seealso{TreeBagger, TreeBagger.quantilePredict,
    ## TreeBagger.oobQuantileError}
    ## @end deftypefn
    function [YFit, YW] = oobQuantilePredict (this, varargin)

      requireOOB (this, 'TreeBagger.oobQuantilePredict');
      [tau, o] = quantileArgs (this, this.X, varargin, ...
                               {'Trees', 'TreeWeights'}, ...
                               'TreeBagger.oobQuantilePredict', false);
      [YFit, YW] = quantileSteps (this, this.X, o.trees, o.tw, ...
                                  this.OOBIndices(:,o.trees), tau, ...
                                  'ensemble');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{err} =} quantileError (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {TreeBagger} {@var{err} =} quantileError (@dots{}, @var{name}, @var{value})
    ## @deftypefnx {TreeBagger} {@var{err} =} quantileError (@var{obj}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {TreeBagger} {@var{err} =} quantileError (@var{obj}, @var{Tbl})
    ##
    ## Quantile loss of a regression ensemble.
    ##
    ## For each quantile probability @math{tau}, @var{err} is the weighted
    ## mean over the observations of the pinball loss, @math{tau (y - q)}
    ## where the response @math{y} is not below the predicted quantile
    ## @math{q} and @math{(1 - tau) (q - y)} where it is.  The quantiles are
    ## predicted as by @code{TreeBagger.quantilePredict}.
    ##
    ## In @qcode{'ensemble'} mode, the default, @var{err} is a row with one
    ## element per quantile.  In @qcode{'cumulative'} mode it has one row per
    ## tree, the loss of the first tree, then of the first two, and so on, and
    ## in @qcode{'individual'} mode one row per tree, each on its own.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Mode'} @tab @tab @qcode{'ensemble'} (default),
    ## @qcode{'cumulative'} or @qcode{'individual'}.
    ## @item @qcode{'Quantile'} @tab @tab A vector of probabilities between 0
    ## and 1.  The default is 0.5.
    ## @item @qcode{'Weights'} @tab @tab A nonnegative vector with one weight
    ## per observation.  The default is uniform.
    ## @end multitable
    ##
    ## @qcode{'Trees'}, @qcode{'TreeWeights'} and @qcode{'UseInstanceForTree'}
    ## are taken as by @code{TreeBagger.quantilePredict}, and
    ## @qcode{'TreeWeights'} may not be given in @qcode{'individual'} mode.
    ##
    ## @var{X} may also be a table @var{Tbl}, whose variables are matched to
    ## the predictors the model was fitted on by name and not by position.
    ## @code{quantileError (@var{obj}, @var{Tbl}, @var{ResponseVarName})}
    ## takes the response from the variable @var{ResponseVarName} names, and
    ## @code{quantileError (@var{obj}, @var{Tbl})} from the variable the model
    ## was fitted on.  The response may also be given beside the table as
    ## @var{Y}.
    ##
    ## @seealso{TreeBagger, TreeBagger.quantilePredict,
    ## TreeBagger.oobQuantileError, TreeBagger.error}
    ## @end deftypefn
    function err = quantileError (this, X, Y, varargin)

      if (nargin < 3 && ! (nargin > 1 && istable (X)))
        error ("TreeBagger.quantileError: too few input arguments.");
      endif

      ## A table carries the response: named in the call, given beside
      ## the table, or the variable the model was fitted on
      if (nargin < 3)
        Y = [];
      endif
      [X, Y, varargin] = tableResponse (this, 'quantileError', X, Y, ...
                                        varargin, nargin > 2);
      [tau, o] = quantileArgs (this, X, varargin, ...
                               {'Mode', 'Trees', 'TreeWeights', ...
                                'UseInstanceForTree', 'Weights'}, ...
                               'TreeBagger.quantileError', true);
      if (! (isnumeric (Y) && isreal (Y) && isvector (Y)
             && numel (Y) == rows (X)))
        error (strcat ("TreeBagger.quantileError: Y must be a real", ...
                       " numeric vector with one element per row in X."));
      endif
      if (isempty (o.w))
        w = ones (rows (X), 1);
      else
        w = o.w;
      endif
      q = quantileSteps (this, X, o.trees, o.tw, o.use, tau, o.mode);
      err = pinballLoss (double (Y(:)), q, tau, w);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{err} =} oobQuantileError (@var{obj})
    ## @deftypefnx {TreeBagger} {@var{err} =} oobQuantileError (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-bag quantile loss of a regression ensemble.
    ##
    ## The loss of @code{TreeBagger.quantileError} on the training data, the
    ## quantiles predicted out of bag as by
    ## @code{TreeBagger.oobQuantilePredict} and the observations weighted by
    ## @code{W}.  @qcode{'Mode'}, @qcode{'Quantile'}, @qcode{'Trees'} and
    ## @qcode{'TreeWeights'} are taken as by @code{TreeBagger.quantileError}.
    ## The ensemble must have been fitted with @qcode{'OOBPrediction'} on.
    ##
    ## In @qcode{'individual'} mode each tree is judged on the observations
    ## out of its own bag alone, and a tree that left nothing out has a
    ## @code{NaN} loss.  MATLAB R2024a fails with an indexing error in that
    ## mode.
    ##
    ## @seealso{TreeBagger, TreeBagger.oobQuantilePredict,
    ## TreeBagger.quantileError, TreeBagger.oobError}
    ## @end deftypefn
    function err = oobQuantileError (this, varargin)

      requireOOB (this, 'TreeBagger.oobQuantileError');
      [tau, o] = quantileArgs (this, this.X, varargin, ...
                               {'Mode', 'Trees', 'TreeWeights'}, ...
                               'TreeBagger.oobQuantileError', true);
      use = this.OOBIndices(:,o.trees);
      q = quantileSteps (this, this.X, o.trees, o.tw, use, tau, o.mode);
      if (strcmp (o.mode, 'individual'))
        err = zeros (numel (o.trees), numel (tau));
        for j = 1:numel (o.trees)
          err(j,:) = pinballLoss (this.Y, q(:,:,j), tau, this.W .* use(:,j));
        endfor
      else
        err = pinballLoss (this.Y, q, tau, this.W);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {@var{B} =} fillprox (@var{obj})
    ## @deftypefnx {TreeBagger} {@var{B} =} fillprox (@dots{}, @var{name}, @var{value})
    ##
    ## Fill the proximity matrix of the training data.
    ##
    ## @var{B} is the ensemble with @code{Proximity} holding the share of the
    ## trees that bring each pair of training observations to the same leaf,
    ## and @code{OutlierMeasure} the outlier measure computed from it.
    ## @qcode{'Trees'} is @qcode{'all'} (default) or a vector of indices of
    ## the trees to use.  @qcode{'NumPrint'} is a nonnegative integer; after
    ## every that many trees a line saying how many are done is printed.
    ##
    ## @seealso{TreeBagger, TreeBagger.mdsprox, CompactTreeBagger.proximity,
    ## CompactTreeBagger.outlierMeasure}
    ## @end deftypefn
    function this = fillprox (this, varargin)

      [NumPrint, args] = takePair (varargin, 'NumPrint', 0);
      if (! (isnumeric (NumPrint) && isscalar (NumPrint) && isreal (NumPrint)
             && NumPrint >= 0 && NumPrint == fix (NumPrint)))
        error (strcat ("TreeBagger.fillprox: 'NumPrint' must be a", ...
                       " nonnegative integer."));
      endif
      [o, errmsg] = bagArgs (args, this.NumTrees, rows (this.X), {'Trees'});
      if (! isempty (errmsg))
        error ("TreeBagger.fillprox: %s", errmsg);
      endif
      this.Proximity = bagProximity (this, this.X, o.trees, NumPrint);
      this.OutlierMeasure = bagOutlier (this.Proximity, this.gY);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {TreeBagger} {[@var{S}, @var{E}] =} mdsprox (@var{obj})
    ## @deftypefnx {TreeBagger} {[@var{S}, @var{E}] =} mdsprox (@dots{}, @var{name}, @var{value})
    ##
    ## Multidimensional scaling of the proximity matrix.
    ##
    ## Applies classical multidimensional scaling, as @code{cmdscale} does, to
    ## the distances @code{1 - Proximity}.  @var{S} holds the scaled
    ## coordinates, one column per positive eigenvalue, and @var{E} the
    ## eigenvalues.  @code{fillprox} must have filled @code{Proximity} first.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.25 0.02 0.73
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Keep'} @tab @tab @qcode{'all'} (default), or a vector of
    ## indices or a logical vector selecting the training observations to
    ## scale.
    ## @item @qcode{'Colors'} @tab @tab A character vector with one color
    ## letter per class.  When given, the scaled coordinates are drawn as
    ## overlaid scatter plots, one per class, a class beyond the number of
    ## letters not drawn; a regression ensemble is drawn in the first color.
    ## @item @qcode{'MDSCoordinates'} @tab @tab Two or three indices of the
    ## columns of @var{S} to draw.  The default is @code{[1, 2]}.  They must
    ## not exceed the number of columns of @var{S} even when nothing is drawn,
    ## as in MATLAB, whose documentation says otherwise.
    ## @end multitable
    ##
    ## @seealso{TreeBagger, TreeBagger.fillprox, CompactTreeBagger.mdsprox,
    ## cmdscale}
    ## @end deftypefn
    function [S, E] = mdsprox (this, varargin)

      P = this.Proximity;
      if (mod (numel (varargin), 2) != 0)
        error ("TreeBagger.mdsprox: name-value arguments must be in pairs.");
      endif
      N = rows (this.X);
      keep = 1:N;
      colors = '';
      coords = [1, 2];
      for i = 1:2:numel (varargin)
        name = varargin{i};
        val = varargin{i+1};
        if (! (ischar (name)
               && any (strcmpi (name, {'Keep', 'Colors', 'MDSCoordinates'}))))
          error (strcat ("TreeBagger.mdsprox: invalid parameter name in", ...
                         " optional pair arguments."));
        endif
        switch (tolower (name))
          case 'keep'
            if (ischar (val) && strcmpi (val, 'all'))
              keep = 1:N;
            elseif (islogical (val) && isvector (val) && numel (val) == N)
              keep = find (val);
            elseif (isnumeric (val) && isvector (val) && isreal (val)
                    && all (val >= 1) && all (val <= N)
                    && all (val == fix (val)))
              keep = double (val(:)');
            else
              error (strcat ("TreeBagger.mdsprox: 'Keep' must be 'all',", ...
                             " a vector of indices of observations, or a", ...
                             " logical vector with one element per", ...
                             " observation."));
            endif
          case 'colors'
            colors = val;
          case 'mdscoordinates'
            coords = val;
        endswitch
      endfor
      g = [];
      if (strcmp (this.Method, 'classification'))
        g = this.gY(keep);
      endif
      [S, E, errmsg] = bagMds (P(keep,keep), g, colors, coords);
      if (! isempty (errmsg))
        error ("TreeBagger.mdsprox: %s", errmsg);
      endif

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
    ## as the first were, their out-of-bag information added.  A proximity
    ## matrix filled by @code{fillprox} is emptied, as it no longer describes
    ## the ensemble; MATLAB keeps it unchanged.
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
    ## out-of-bag information and predictor importance, which are then
    ## joined.  A proximity matrix filled by @code{fillprox} is emptied.
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
          || ! isequal (this.PredictorNames, other.PredictorNames)
          || ! isequal (this.CategoricalPredictors,
                        other.CategoricalPredictors))
        error (strcat ("TreeBagger.append: the two ensembles must be", ...
                       " fitted on the same observations and predictors."));
      endif
      if (! (isequal (this.ClassNames, other.ClassNames)
             && isequal (this.Prior, other.Prior)))
        error (strcat ("TreeBagger.append: the two ensembles must have", ...
                       " the same classes and priors."));
      endif
      if (this.ComputeOOBPrediction != other.ComputeOOBPrediction
          || (this.ComputeOOBPredictorImportance
              != other.ComputeOOBPredictorImportance))
        error (strcat ("TreeBagger.append: the two ensembles have", ...
                       " incompatible out-of-bag flags."));
      endif
      this.Trees = [this.Trees; other.Trees];
      this.TreeClassIdx = [this.TreeClassIdx; other.TreeClassIdx];
      this.NumTrees = numel (this.Trees);
      this.InBag = [this.InBag, other.InBag];
      if (this.ComputeOOBPrediction)
        this.OOBIndices = [this.OOBIndices, other.OOBIndices];
        this.OOBInstanceWeight = sum (this.OOBIndices, 2);
      endif
      if (this.ComputeOOBPredictorImportance)
        this.PermDelta = [this.PermDelta; other.PermDelta];
      endif
      this = updateStats (this);

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
        counts = accumarray (idx(:), 1, [N, 1]);
        this.InBag(:,end+1) = sparse (counts);
        if (this.ComputeOOBPrediction)
          this.OOBIndices(:,end+1) = counts == 0;
        endif
        if (this.ComputeOOBPredictorImportance)
          this.PermDelta(end+1,:,:) = permDelta (this, this.Trees{end}, ...
                                                 this.TreeClassIdx{end}, ...
                                                 counts == 0);
        endif
        this.NumTrees = numel (this.Trees);
        if (NumPrint > 0 && mod (this.NumTrees, NumPrint) == 0)
          printf ("Tree %d done.\n", this.NumTrees);
        endif
      endfor
      if (this.ComputeOOBPrediction)
        this.OOBInstanceWeight = sum (this.OOBIndices, 2);
      endif
      this = updateStats (this);

    endfunction

    ## What permuting each predictor among a tree's out-of-bag observations
    ## does to the tree: a 1xPx3 array holding the rise in weighted error,
    ## the fall in weighted mean margin and the number of margins lowered
    ## less the number raised, all zero when nothing is out of bag.
    function d = permDelta (this, tree, tidx, oob)

      p = columns (this.X);
      d = zeros (1, p, 3);
      r = find (oob);
      w = this.W(r);
      if (isempty (r) || ! (sum (w) > 0))
        return;
      endif
      w /= sum (w);
      X = this.X(r,:);
      n = numel (r);
      if (strcmp (this.Method, 'classification'))
        g = this.gY(r);
        K = classCount (this.ClassNames);
        [e0, m0] = treeMissMargin (tree, X, g, tidx, K);
        for v = 1:p
          Xp = X;
          Xp(:,v) = X(randperm (n), v);
          [e1, m1] = treeMissMargin (tree, Xp, g, tidx, K);
          d(1,v,1) = sum (w .* (e1 - e0));
          d(1,v,2) = sum (w .* (m0 - m1));
          d(1,v,3) = sum (m1 < m0) - sum (m1 > m0);
        endfor
      else
        y = this.Y(r);
        e0 = (predict (tree, X) - y) .^ 2;
        for v = 1:p
          Xp = X;
          Xp(:,v) = X(randperm (n), v);
          d(1,v,1) = sum (w .* ((predict (tree, Xp) - y) .^ 2 - e0));
        endfor
      endif

    endfunction

    ## The statistics derived from all the trees: the split contributions and
    ## counts, and the permuted importances, each the mean over the trees
    ## divided by the standard deviation over the trees.  A filled proximity
    ## matrix describes other trees, and is emptied.
    function this = updateStats (this)

      p = columns (this.X);
      [this.DeltaCriterionDecisionSplit, this.NumPredictorSplit] = ...
        bagSplitStats (this.Trees, p);
      this.SurrogateAssociation = eye (p);
      if (this.ComputeOOBPredictorImportance)
        mu = mean (this.PermDelta, 1);
        sd = std (this.PermDelta, 0, 1);
        imp = mu ./ sd;
        imp(mu == 0) = 0;
        imp = reshape (imp, p, 3)';
        this.OOBPermutedPredictorDeltaError = imp(1,:);
        if (strcmp (this.Method, 'classification'))
          this.OOBPermutedPredictorDeltaMeanMargin = imp(2,:);
          this.OOBPermutedPredictorCountRaiseMargin = imp(3,:);
        else
          this.OOBPermutedPredictorDeltaMeanMargin = zeros (1, 0);
          this.OOBPermutedPredictorCountRaiseMargin = zeros (1, 0);
        endif
      endif
      this.Proximity = [];
      this.OutlierMeasure = [];

    endfunction

    ## Quantiles of the response at the rows of Xq, for the trees given, their
    ## weights and the MxT matrix of which tree may answer for which row.
    ## MODE 'ensemble' gives an MxQ matrix, 'cumulative' and 'individual'
    ## an MxQxT array, the trees taken up to each one or each alone.  YW is
    ## the response weight matrix of the last step.
    function [q, YW] = quantileSteps (this, Xq, trees, tw, use, tau, mode)

      N = rows (this.X);
      n = rows (Xq);
      T = numel (trees);
      Ltr = bagLeaves (this, this.X, trees);
      Lq = bagLeaves (this, Xq, trees);
      sample = quantile (this.Y, tau);
      sample = sample(:)';
      if (strcmp (mode, 'ensemble'))
        q = zeros (n, numel (tau));
      else
        q = zeros (n, numel (tau), T);
      endif
      num = sparse (N, n);
      den = zeros (1, n);
      for j = 1:T
        [i, ~, c] = find (this.InBag(:,trees(j)));
        nn = max ([Ltr(:,j); Lq(:,j)]);
        mass = accumarray (Ltr(i,j), c, [nn, 1]);
        A = sparse (i, Ltr(i,j), c ./ mass(Ltr(i,j)), N, nn);
        k = find (use(:,j) & mass(Lq(:,j)) > 0)';
        Cj = A * sparse (Lq(k,j), k, tw(j), nn, n);
        dj = zeros (1, n);
        dj(k) = tw(j);
        if (strcmp (mode, 'individual'))
          [q(:,:,j), YW] = weighQuantiles (this.Y, Cj, dj, tau, sample);
        else
          num += Cj;
          den += dj;
          if (strcmp (mode, 'cumulative'))
            [q(:,:,j), YW] = weighQuantiles (this.Y, num, den, tau, sample);
          endif
        endif
      endfor
      if (strcmp (mode, 'ensemble'))
        [q, YW] = weighQuantiles (this.Y, num, den, tau, sample);
      endif

    endfunction

    function requireImportance (this)
      if (! this.ComputeOOBPredictorImportance)
        error (strcat ("TreeBagger: out-of-bag permutations were not", ...
                       " kept; fit with 'OOBPredictorImportance' set to", ...
                       " 'on'."));
      endif
    endfunction

    function requireProximity (this, val)
      if (isempty (val))
        error (strcat ("TreeBagger: the proximity matrix was not filled;", ...
                       " call fillprox first."));
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

## The value of a Name-Value pair, removed from ARGS, or DEF when it is not
## given.  Pairs that are not well formed are left for bagArgs to refuse.
function [val, args] = takePair (args, name, def)

  val = def;
  if (mod (numel (args), 2) != 0)
    return;
  endif
  keep = true (size (args));
  for i = 1:2:numel (args)
    if (ischar (args{i}) && strcmpi (args{i}, name))
      val = args{i+1};
      keep(i:i+1) = false;
    endif
  endfor
  args = args(keep);

endfunction

## Validate a quantile method's call: a regression ensemble, the predictors,
## 'Quantile' and the pairs bagArgs takes.  An error method defaults 'Mode' to
## 'ensemble'.
function [tau, o] = quantileArgs (M, X, args, allowed, caller, isError)

  if (! strcmp (M.Method, 'regression'))
    error ("%s: quantiles are defined only for regression.", caller);
  endif
  if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
    error ("%s: X must be a real numeric matrix.", caller);
  endif
  if (columns (X) != numel (M.PredictorNames))
    error ("%s: X must have one column per predictor.", caller);
  endif
  [tau, args] = takePair (args, 'Quantile', 0.5);
  if (! (isnumeric (tau) && isreal (tau) && isvector (tau)
         && all (tau >= 0 & tau <= 1)))
    error (strcat ("%s: 'Quantile' must be a vector of probabilities", ...
                   " between 0 and 1."), caller);
  endif
  tau = double (tau(:)');
  if (isError && mod (numel (args), 2) == 0
      && ! any (cellfun (@(a) ischar (a) && strcmpi (a, 'Mode'), ...
                         args(1:2:end))))
    args(end+1:end+2) = {'Mode', 'ensemble'};
  endif
  [o, errmsg] = bagArgs (args, M.NumTrees, rows (X), allowed);
  if (! isempty (errmsg))
    error ("%s: %s", caller, errmsg);
  endif

endfunction

## Weighted quantiles of Y for each column of the weights NUM ./ DEN, and the
## normalized weights YW.  A column without weight takes SAMPLE and weighs
## every observation equally.
function [q, YW] = weighQuantiles (y, num, den, tau, sample)

  [N, n] = size (num);
  q = repmat (sample, n, 1);
  have = den > 0;
  dinv = zeros (n, 1);
  dinv(have) = 1 ./ den(have);
  YW = num * spdiags (dinv, 0, n, n);
  for k = find (have)
    [i, ~, w] = find (YW(:,k));
    [ys, o] = sort (y(i));
    w = w(o);
    if (numel (ys) == 1)
      q(k,:) = ys;
      continue;
    endif
    F = cumsum (w) - w / 2;
    qk = interp1 (F, ys, tau);
    qk(tau <= F(1)) = ys(1);
    qk(tau >= F(end)) = ys(end);
    q(k,:) = qk;
  endfor
  if (! all (have))
    YW(:,! have) = 1 / N;
  endif

endfunction

## Weighted mean pinball loss of the quantiles Q, MxQxL, against Y, as an LxQ
## matrix.  An observation with a missing response or zero weight is left
## out.
function err = pinballLoss (y, q, tau, w)

  w = w(:);
  w(isnan (y)) = 0;
  d = y - q;
  loss = max (tau .* d, (tau - 1) .* d);
  loss(w == 0,:,:) = 0;
  err = permute (sum (w .* loss, 1), [3, 2, 1]) / sum (w);

endfunction

## Whether a tree misclassifies each observation, of class index G, and the
## margin it gives it, its scores laid out over the K classes of the ensemble.
function [miss, m] = treeMissMargin (tree, X, g, tidx, K)

  [~, s] = predict (tree, X);
  S = zeros (rows (X), K);
  S(:,tidx) = s;
  [~, k] = max (S, [], 2);
  miss = double (k != g);
  m = marginsOf (S, g, 1);

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

%!demo
%! ## A 90% prediction interval for sepal length from the other three
%! ## measurements, read off the quantiles of a regression forest.
%! load fisheriris
%! rng (42);
%! B = TreeBagger (100, meas(:,2:4), meas(:,1), 'Method', 'regression');
%! q = quantilePredict (B, meas([1, 51, 101], 2:4), ...
%!                      'Quantile', [0.05, 0.5, 0.95])

%!demo
%! ## Which measurements a forest relies on: permuting a petal measurement
%! ## among each tree's out-of-bag flowers raises the error the most.
%! load fisheriris
%! rng (42);
%! B = TreeBagger (50, meas, species, 'OOBPredictorImportance', 'on');
%! bar (B.OOBPermutedPredictorDeltaError);
%! set (gca, 'xticklabel', {'SL', 'SW', 'PL', 'PW'});
%! ylabel ('Rise in out-of-bag error');

%!demo
%! ## Scale the proximities of a random forest to two dimensions and draw
%! ## the flowers one color per species.
%! load fisheriris
%! rng (42);
%! B = fillprox (TreeBagger (50, meas, species));
%! mdsprox (B, 'Colors', 'rgb');

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

%!test  # MATLAB parity: the split statistics of a regression ensemble
%! load fisheriris
%! B = TreeBagger (2, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'SampleWithReplacement', 'off', 'InBagFraction', 1, ...
%!                 'NumPredictorsToSample', 'all');
%! assert_equal (B.DeltaCriterionDecisionSplit, ...
%!               [0.001540487516843, 0.025160029558674, ...
%!                0.000497162933685], 1e-14);
%! assert_equal (B.NumPredictorSplit, ...
%!               [0.521739130434783, 1.043478260869565, ...
%!                0.434782608695652], 1e-14);
%! assert_equal (B.SurrogateAssociation, eye (3));

%!test  # MATLAB parity: the split statistics of a classification ensemble
%! load fisheriris
%! B = TreeBagger (2, meas, species, 'SampleWithReplacement', 'off', ...
%!                 'InBagFraction', 1, 'NumPredictorsToSample', 'all', ...
%!                 'MinLeafSize', 5);
%! assert_equal (B.DeltaCriterionDecisionSplit, ...
%!               [0.000477777777778, 0, 0.072985238862050, ...
%!                0.051959205582394], 1e-14);
%! assert_equal (B.NumPredictorSplit, [0.4, 0, 1.2, 0.4], 1e-14);

%!test  # MATLAB parity: split statistics are per-tree shares and means
%! load fisheriris
%! rng (1);
%! B = growTrees (TreeBagger (3, meas, species), 2);
%! n = zeros (1, 4);
%! d = zeros (1, 4);
%! for t = 1:5
%!   tr = B.Trees{t};
%!   br = find (tr.Children(:,1) > 0);
%!   n += accumarray (tr.CutPredictorIndex(br), 1, [4, 1])' / numel (br);
%!   d += predictorImportance (tr) / 5;
%! endfor
%! assert_equal (B.NumPredictorSplit, n, 1e-14);
%! assert_equal (B.DeltaCriterionDecisionSplit, d, 1e-14);

%!test  # MATLAB parity: a tree without splits adds nothing to the statistics
%! load fisheriris
%! B = TreeBagger (3, meas, species, 'MinLeafSize', 100);
%! assert_equal (B.NumPredictorSplit, zeros (1, 4));
%! assert_equal (B.DeltaCriterionDecisionSplit, zeros (1, 4));

%!test  # MATLAB parity: the quantiles of a regression ensemble
%! load fisheriris
%! B = TreeBagger (2, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'SampleWithReplacement', 'off', 'InBagFraction', 1, ...
%!                 'NumPredictorsToSample', 'all');
%! q = quantilePredict (B, meas([1, 51, 101, 150], 2:4), ...
%!                      'Quantile', [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1]);
%! assert_equal (q, [4.6, 4.6, 4.9, 5.0, 5.1, 5.5, 5.5; ...
%!                   6.4, 6.4, 6.625, 6.7, 6.925, 7.0, 7.0; ...
%!                   6.3, 6.3, 6.45, 6.7, 6.825, 6.9, 6.9; ...
%!                   5.9, 5.9, 5.975, 6.1, 6.225, 6.3, 6.3], 1e-14);

%!test  # MATLAB parity: observation weights enter only through the sample
%! load fisheriris
%! w = [10 * ones(50, 1); ones(100, 1)];
%! B = TreeBagger (1, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'SampleWithReplacement', 'off', 'InBagFraction', 1, ...
%!                 'NumPredictorsToSample', 'all', 'Weights', w);
%! q = quantilePredict (B, meas([1, 51, 101], 2:4), ...
%!                      'Quantile', [0.25, 0.5, 0.75]);
%! assert_equal (q, [4.9, 5.0, 5.1; 6.625, 6.7, 6.925; ...
%!                   6.45, 6.7, 6.825], 1e-14);

%!test  # response weights are sample counts over leaf sizes, averaged
%! load fisheriris
%! rng (1);
%! B = TreeBagger (4, meas(:,2:4), meas(:,1), 'Method', 'regression');
%! [~, YW] = quantilePredict (B, meas(1:3,2:4));
%! W = zeros (150, 3);
%! for t = 1:4
%!   [~, ntr] = predict (B.Trees{t}, B.X);
%!   [~, nq] = predict (B.Trees{t}, meas(1:3,2:4));
%!   for k = 1:3
%!     in = (ntr == nq(k)) .* full (B.InBag(:,t));
%!     W(:,k) += in / B.Trees{t}.NodeSize(nq(k)) / 4;
%!   endfor
%! endfor
%! assert_equal (issparse (YW), true);
%! assert_equal (full (YW), W, 1e-15);

%!test  # quantiles interpolate sorted responses at mid cumulative weights
%! load fisheriris
%! rng (1);
%! B = TreeBagger (5, meas(:,2:4), meas(:,1), 'Method', 'regression');
%! tau = [0.001, 0.3, 0.62, 0.999];
%! [q, YW] = quantilePredict (B, meas(7,2:4), 'Quantile', tau);
%! [i, ~, w] = find (YW);
%! [ys, o] = sort (B.Y(i));
%! w = w(o);
%! F = cumsum (w) - w / 2;
%! e = interp1 (F, ys, tau);
%! e(tau <= F(1)) = ys(1);
%! e(tau >= F(end)) = ys(end);
%! assert_equal (q, e, 1e-14);

%!test  # MATLAB parity: a row no tree may answer for takes the sample quantile
%! load fisheriris
%! rng (1);
%! B = TreeBagger (3, meas(:,2:4), meas(:,1), 'Method', 'regression');
%! U = true (2, 3);
%! U(1,:) = false;
%! [q, YW] = quantilePredict (B, meas(1:2,2:4), ...
%!                           'Quantile', [0.1, 0.5, 0.9], ...
%!                           'UseInstanceForTree', U);
%! assert_equal (q(1,:), [4.8, 5.8, 6.9], 1e-14);
%! assert_equal (full (YW(:,1)), ones (150, 1) / 150, 1e-15);

%!test  # MATLAB parity: an observation in every bag takes the sample quantile
%! load fisheriris
%! rng (22);
%! B = TreeBagger (1, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'OOBPrediction', 'on');
%! r = find (! B.OOBIndices, 1);
%! q = oobQuantilePredict (B, 'Quantile', [0.1, 0.5, 0.9]);
%! assert_equal (q(r,:), [4.8, 5.8, 6.9], 1e-14);

%!test  # out-of-bag quantiles use only the trees that left a row out
%! load fisheriris
%! rng (1);
%! B = TreeBagger (6, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'OOBPrediction', 'on');
%! tau = [0.2, 0.8];
%! [qo, Wo] = oobQuantilePredict (B, 'Quantile', tau);
%! [q, W] = quantilePredict (B, B.X, 'Quantile', tau, ...
%!                           'UseInstanceForTree', B.OOBIndices);
%! assert_equal (qo, q);
%! assert_equal (Wo, W);

%!test  # MATLAB parity: the quantile loss is the mean pinball loss
%! load fisheriris
%! B = TreeBagger (2, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'SampleWithReplacement', 'off', 'InBagFraction', 1, ...
%!                 'NumPredictorsToSample', 'all');
%! e = quantileError (B, meas(:,2:4), meas(:,1), 'Quantile', [0.25, 0.5, 0.9]);
%! assert_equal (e, [0.069916666666667, 0.088666666666667, ...
%!                   0.035026666666667], 1e-14);

%!test  # MATLAB parity: the quantile loss weighs the observations
%! load fisheriris
%! B = TreeBagger (2, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'SampleWithReplacement', 'off', 'InBagFraction', 1, ...
%!                 'NumPredictorsToSample', 'all');
%! w = [5 * ones(50, 1); ones(100, 1)];
%! e = quantileError (B, meas(:,2:4), meas(:,1), 'Quantile', [0.25, 0.9], ...
%!                    'Weights', w);
%! assert_equal (e, [0.064821428571429, 0.036268571428571], 1e-14);

%!test  # MATLAB parity: the three modes of the quantile loss
%! load fisheriris
%! rng (1);
%! B = TreeBagger (4, meas(:,2:4), meas(:,1), 'Method', 'regression');
%! X = meas(:,2:4);
%! y = meas(:,1);
%! tau = [0.25, 0.75];
%! e = quantileError (B, X, y, 'Quantile', tau);
%! c = quantileError (B, X, y, 'Quantile', tau, 'Mode', 'cumulative');
%! i = quantileError (B, X, y, 'Quantile', tau, 'Mode', 'individual');
%! assert_equal (size (e), [1, 2]);
%! assert_equal (size (c), [4, 2]);
%! assert_equal (c(end,:), e, 1e-15);
%! assert_equal (i(3,:), quantileError (B, X, y, 'Quantile', tau, ...
%!                                      'Trees', 3), 1e-15);

%!test  # MATLAB parity: the out-of-bag quantile loss in ensemble mode
%! load fisheriris
%! rng (1);
%! B = TreeBagger (6, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'OOBPrediction', 'on');
%! tau = [0.1, 0.5, 0.9];
%! d = B.Y - oobQuantilePredict (B, 'Quantile', tau);
%! L = mean (max (tau .* d, (tau - 1) .* d));
%! assert_equal (oobQuantileError (B, 'Quantile', tau), L, 1e-15);
%! assert_equal (size (oobQuantileError (B, 'Mode', 'cumulative')), [6, 1]);

%!test  # each tree's out-of-bag quantile loss is on its own out-of-bag rows
%! load fisheriris
%! rng (1);
%! B = TreeBagger (3, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'OOBPrediction', 'on');
%! e = oobQuantileError (B, 'Mode', 'individual');
%! r = B.OOBIndices(:,2);
%! d = B.Y(r) - quantilePredict (B, B.X(r,:), 'Trees', 2);
%! assert_equal (size (e), [3, 1]);
%! assert_equal (e(2), mean (abs (d)) / 2, 1e-15);

%!test  # MATLAB parity: out-of-bag importance turns out-of-bag prediction on
%! load fisheriris
%! rng (1);
%! B = TreeBagger (3, meas, species, 'OOBPredictorImportance', 'on');
%! assert_equal (B.ComputeOOBPrediction, true);
%! assert_equal (B.ComputeOOBPredictorImportance, true);
%! assert_equal (size (B.OOBPermutedPredictorDeltaError), [1, 4]);
%! assert_equal (size (B.OOBPermutedPredictorDeltaMeanMargin), [1, 4]);
%! assert_equal (size (B.OOBPermutedPredictorCountRaiseMargin), [1, 4]);

%!test  # MATLAB parity: a regression ensemble has no margin importances
%! load fisheriris
%! rng (1);
%! B = TreeBagger (3, meas(:,2:4), meas(:,1), 'Method', 'regression', ...
%!                 'OOBPredictorImportance', 'on');
%! assert_equal (size (B.OOBPermutedPredictorDeltaError), [1, 3]);
%! assert_equal (B.OOBPermutedPredictorDeltaMeanMargin, zeros (1, 0));
%! assert_equal (B.OOBPermutedPredictorCountRaiseMargin, zeros (1, 0));

%!test  # MATLAB parity: importance is the mean change over its deviation
%! load fisheriris
%! rng (1);
%! B = TreeBagger (5, meas, species, 'OOBPredictorImportance', 'on');
%! D = B.PermDelta;
%! assert_equal (B.OOBPermutedPredictorDeltaError, ...
%!               mean (D(:,:,1)) ./ std (D(:,:,1)), 1e-14);
%! assert_equal (B.OOBPermutedPredictorDeltaMeanMargin, ...
%!               mean (D(:,:,2)) ./ std (D(:,:,2)), 1e-14);
%! assert_equal (B.OOBPermutedPredictorCountRaiseMargin, ...
%!               mean (D(:,:,3)) ./ std (D(:,:,3)), 1e-14);

%!test  # MATLAB parity: a predictor no tree splits on has zero importance
%! load fisheriris
%! rng (1);
%! B = TreeBagger (10, [meas, ones(150, 1)], species, ...
%!                 'OOBPredictorImportance', 'on');
%! assert_equal (B.NumPredictorSplit(5), 0);
%! assert_equal (B.OOBPermutedPredictorDeltaError(5), 0);
%! assert_equal (B.OOBPermutedPredictorDeltaMeanMargin(5), 0);
%! assert_equal (B.OOBPermutedPredictorCountRaiseMargin(5), 0);

%!test  # MATLAB parity: the importance of a single tree is zero or infinite
%! load fisheriris
%! rng (9);
%! B = TreeBagger (1, meas, species, 'OOBPredictorImportance', 'on');
%! imp = B.OOBPermutedPredictorDeltaError;
%! assert_equal (all (imp == 0 | isinf (imp)), true);

%!test  # MATLAB parity: permuting a petal measurement matters most
%! load fisheriris
%! rng (1);
%! B = TreeBagger (30, meas, species, 'OOBPredictorImportance', 'on');
%! imp = B.OOBPermutedPredictorDeltaError;
%! assert_equal (min (imp(3:4)) > max (imp(1:2)), true);
%! assert_equal (all (B.OOBPermutedPredictorCountRaiseMargin(3:4) > 0), true);

%!test  # growing and appending extend the importances
%! load fisheriris
%! rng (1);
%! B = TreeBagger (3, meas, species, 'OOBPredictorImportance', 'on');
%! B = growTrees (B, 2);
%! assert_equal (size (B.PermDelta), [5, 4, 3]);
%! B = append (B, TreeBagger (2, meas, species, ...
%!                            'OOBPredictorImportance', 'on'));
%! assert_equal (size (B.PermDelta), [7, 4, 3]);
%! D = B.PermDelta(:,:,1);
%! assert_equal (B.OOBPermutedPredictorDeltaError, mean (D) ./ std (D), 1e-14);

%!test  # MATLAB parity: the proximity of the training observations
%! load fisheriris
%! B = TreeBagger (2, meas, species, 'SampleWithReplacement', 'off', ...
%!                 'InBagFraction', 1, 'NumPredictorsToSample', 'all', ...
%!                 'MinLeafSize', 5);
%! B = fillprox (B);
%! k = [1, 51, 101, 150];
%! assert_equal (B.Proximity(k,k), [1, 0, 0, 0; 0, 1, 0, 0; ...
%!                                  0, 0, 1, 1; 0, 0, 1, 1]);

%!test  # MATLAB parity: fillprox over chosen trees
%! load fisheriris
%! rng (1);
%! B = fillprox (TreeBagger (6, meas, species, 'MinLeafSize', 5), ...
%!               'Trees', [2, 5]);
%! P = zeros (150);
%! for t = [2, 5]
%!   [~, ~, nd] = predict (B.Trees{t}, meas);
%!   P += nd == nd';
%! endfor
%! assert_equal (B.Proximity, P / 2);

%!test  # the filled matrices agree with the compact ensemble's methods
%! load fisheriris
%! rng (1);
%! B = fillprox (TreeBagger (4, meas, species, 'MinLeafSize', 5));
%! C = compact (B);
%! assert_equal (B.Proximity, proximity (C, meas));
%! assert_equal (B.OutlierMeasure, outlierMeasure (C, meas, 'Labels', species));

%!test  # MATLAB parity: a regression outlier measure is over every row
%! load fisheriris
%! rng (1);
%! B = fillprox (TreeBagger (4, meas(:,2:4), meas(:,1), ...
%!                           'Method', 'regression'));
%! assert_equal (B.OutlierMeasure, outlierMeasure (compact (B), meas(:,2:4)));

%!test  # MATLAB parity: fillprox prints its progress
%! load fisheriris
%! B = TreeBagger (3, meas, species);
%! out = evalc ("fillprox (B, 'NumPrint', 1);");
%! assert_equal (out, sprintf ("Tree 1 done.\nTree 2 done.\nTree 3 done.\n"));

%!test  # MATLAB parity: the scaling of the proximity matrix
%! load fisheriris
%! B = TreeBagger (2, meas, species, 'SampleWithReplacement', 'off', ...
%!                 'InBagFraction', 1, 'NumPredictorsToSample', 'all', ...
%!                 'MinLeafSize', 5);
%! [S, E] = mdsprox (fillprox (B));
%! assert_equal (size (S), [150, 5]);
%! assert_equal (E(1:5)', [23.511509991334723, 20.652625175036377, ...
%!                         5.054560979007729, 3, 2.627970521287837], 1e-12);
%! assert_equal (abs (S([1, 51, 101], 1:2)), ...
%!               [0.539357634260452, 0.083561333480090; ...
%!                0.399117306477463, 0.428703363387418; ...
%!                0.228627699102420, 0.556632583919985], 1e-12);

%!test  # MATLAB parity: 'Keep' scales the chosen observations alone
%! load fisheriris
%! rng (1);
%! B = fillprox (TreeBagger (10, meas, species, 'MinLeafSize', 5));
%! k = 1:3:150;
%! [S, E] = mdsprox (B, 'Keep', k);
%! [S0, E0] = cmdscale (1 - B.Proximity(k,k));
%! assert_equal (S, S0);
%! assert_equal (E, E0);
%! assert_equal (mdsprox (B, 'Keep', mod (1:150, 3) == 1), S0);

%!test  # the scaled coordinates are drawn one class per color
%! load fisheriris
%! rng (1);
%! B = fillprox (TreeBagger (5, meas, species, 'MinLeafSize', 5));
%! h = figure ('visible', 'off');
%! unwind_protect
%!   mdsprox (B, 'Colors', 'rgb');
%!   assert_equal (numel (get (gca, 'children')), 3);
%! unwind_protect_cleanup
%!   close (h);
%! end_unwind_protect

%!test  # a regression ensemble is drawn in the first color
%! load fisheriris
%! rng (1);
%! B = fillprox (TreeBagger (5, meas(:,2:4), meas(:,1), ...
%!                           'Method', 'regression'));
%! h = figure ('visible', 'off');
%! unwind_protect
%!   mdsprox (B, 'Colors', 'rgb');
%!   kids = get (gca, 'children');
%!   assert_equal (numel (kids), 1);
%!   assert_equal (get (kids, 'color'), [1, 0, 0]);
%! unwind_protect_cleanup
%!   close (h);
%! end_unwind_protect

%!error<TreeBagger: the proximity matrix was not filled; call fillprox first.> ...
%! load fisheriris
%! B = growTrees (fillprox (TreeBagger (2, meas, species)), 1);
%! B.Proximity

## Test input validation
%!shared x, y, B, C, R, Q
%! load fisheriris
%! x = meas;
%! y = species;
%! B = TreeBagger (3, x, y, 'OOBPrediction', 'on');
%! C = compact (B);
%! R = TreeBagger (2, x(:,2:4), x(:,1), 'Method', 'regression');
%! Q = TreeBagger (2, x(:,2:4), x(:,1), 'Method', 'regression', ...
%!                 'OOBPrediction', 'on');
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
%!error<TreeBagger: 'PredictorNames' must be a cell array of character vectors with one element per column of X.> ...
%! TreeBagger (1, x, y, 'PredictorNames', {'a'})
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
%!error<TreeBagger.append: the two ensembles have incompatible out-of-bag flags.> ...
%! append (B, TreeBagger (1, x, y, 'OOBPredictorImportance', 'on'))
%!error<TreeBagger: out-of-bag permutations were not kept; fit with 'OOBPredictorImportance' set to 'on'.> ...
%! B.OOBPermutedPredictorDeltaError
%!error<TreeBagger: out-of-bag permutations were not kept; fit with 'OOBPredictorImportance' set to 'on'.> ...
%! B.OOBPermutedPredictorDeltaMeanMargin
%!error<TreeBagger: out-of-bag permutations were not kept; fit with 'OOBPredictorImportance' set to 'on'.> ...
%! B.OOBPermutedPredictorCountRaiseMargin
%!error<TreeBagger: the proximity matrix was not filled; call fillprox first.> ...
%! B.Proximity
%!error<TreeBagger: the proximity matrix was not filled; call fillprox first.> ...
%! B.OutlierMeasure
%!error<TreeBagger.quantilePredict: too few input arguments.> ...
%! quantilePredict (R)
%!error<TreeBagger.quantilePredict: quantiles are defined only for regression.> ...
%! quantilePredict (B, x)
%!error<TreeBagger.quantilePredict: X must be a real numeric matrix.> ...
%! quantilePredict (R, {1})
%!error<TreeBagger.quantilePredict: X must have one column per predictor.> ...
%! quantilePredict (R, x)
%!error<TreeBagger.quantilePredict: 'Quantile' must be a vector of probabilities between 0 and 1.> ...
%! quantilePredict (R, x(:,2:4), 'Quantile', 2)
%!error<TreeBagger.quantilePredict: name-value arguments must be in pairs.> ...
%! quantilePredict (R, x(:,2:4), 'Trees')
%!error<TreeBagger.quantilePredict: invalid parameter name in optional pair arguments.> ...
%! quantilePredict (R, x(:,2:4), 'Mode', 'ensemble')
%!error<TreeBagger.oobQuantilePredict: out-of-bag information was not kept; fit with 'OOBPrediction' set to 'on'.> ...
%! oobQuantilePredict (R)
%!error<TreeBagger.oobQuantilePredict: invalid parameter name in optional pair arguments.> ...
%! oobQuantilePredict (Q, 'UseInstanceForTree', true (150, 2))
%!error<TreeBagger.quantileError: too few input arguments.> ...
%! quantileError (R, x)
%!error<TreeBagger.quantileError: Y must be a real numeric vector with one element per row in X.> ...
%! quantileError (R, x(:,2:4), y)
%!error<TreeBagger.quantileError: 'Mode' must be 'cumulative', 'individual' or 'ensemble'.> ...
%! quantileError (R, x(:,2:4), x(:,1), 'Mode', 'all')
%!error<TreeBagger.quantileError: 'TreeWeights' cannot be used in 'individual' mode.> ...
%! quantileError (R, x(:,2:4), x(:,1), 'Mode', 'individual', ...
%!                'TreeWeights', [1, 1])
%!error<TreeBagger.oobQuantileError: out-of-bag information was not kept; fit with 'OOBPrediction' set to 'on'.> ...
%! oobQuantileError (R)
%!error<TreeBagger.oobQuantileError: quantiles are defined only for regression.> ...
%! oobQuantileError (B)
%!error<TreeBagger.oobQuantileError: invalid parameter name in optional pair arguments.> ...
%! oobQuantileError (Q, 'Weights', ones (150, 1))
%!error<TreeBagger.fillprox: 'NumPrint' must be a nonnegative integer.> ...
%! fillprox (B, 'NumPrint', -1)
%!error<TreeBagger.fillprox: name-value arguments must be in pairs.> ...
%! fillprox (B, 'Trees')
%!error<TreeBagger.fillprox: 'Trees' must be 'all' or a vector of indices of trees in the ensemble.> ...
%! fillprox (B, 'Trees', 9)
%!error<TreeBagger.fillprox: invalid parameter name in optional pair arguments.> ...
%! fillprox (B, 'Keep', 1)
%!error<TreeBagger: the proximity matrix was not filled; call fillprox first.> ...
%! mdsprox (B)
%!error<TreeBagger.mdsprox: name-value arguments must be in pairs.> ...
%! mdsprox (fillprox (B), 'Keep')
%!error<TreeBagger.mdsprox: invalid parameter name in optional pair arguments.> ...
%! mdsprox (fillprox (B), 'Data', 'proximity')
%!error<TreeBagger.mdsprox: 'Keep' must be 'all', a vector of indices of observations, or a logical vector with one element per observation.> ...
%! mdsprox (fillprox (B), 'Keep', 0)
%!error<TreeBagger.mdsprox: 'Colors' must be a character vector or a string scalar.> ...
%! mdsprox (fillprox (B), 'Colors', 1)
%!error<TreeBagger.mdsprox: 'MDSCoordinates' must be a vector of two or three positive integers.> ...
%! mdsprox (fillprox (B), 'MDSCoordinates', [1, 2, 3, 4])
%!error<TreeBagger.mdsprox: 'MDSCoordinates' must not exceed the number of scaled coordinates.> ...
%! mdsprox (fillprox (B), 'MDSCoordinates', [1, 500])

%!error<TreeBagger: not all 'ClassNames' are present in Y.> ...
%! load fisheriris
%! TreeBagger (2, meas, species, 'ClassNames', [3, 1, 2])

%!test  # classes named as text for numeric labels keep the labels' type
%! load fisheriris
%! Y = [ones(50, 1); 2 * ones(50, 1); 3 * ones(50, 1)];
%! B = TreeBagger (2, meas, Y, 'ClassNames', {'3', '1', '2'});
%! assert_equal (B.ClassNames, [3; 1; 2]);

%!shared X, yb, yr
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

%!test  # every bagged tree takes the categorical predictors and options
%! B = TreeBagger (3, X, yb, 'CategoricalPredictors', 1, ...
%!                 'MaxNumCategories', 3, 'AlgorithmForCategorical', 'pca');
%! assert_equal (B.Trees{1}.CategoricalPredictors, 1);
%! assert_equal (any (strcmpi (B.TreeArguments, 'MaxNumCategories')), true);
%! assert_equal (any (strcmpi (B.TreeArguments, 'AlgorithmForCategorical')), ...
%!               true);
%! Br = TreeBagger (3, X, yr, 'Method', 'regression', ...
%!                  'CategoricalPredictors', 1);
%! assert_equal (Br.Trees{1}.CategoricalPredictors, 1);

%!test  # the property records the resolved indices, empty when none is named
%! B = TreeBagger (3, X, yb);
%! assert_equal (B.CategoricalPredictors, []);
%! B = TreeBagger (3, X, yb, 'CategoricalPredictors', 1);
%! assert_equal (B.CategoricalPredictors, 1);

%!test  # a logical vector and 'all' are resolved to indices
%! B = TreeBagger (3, X, yb, 'CategoricalPredictors', [true, false]);
%! assert_equal (B.CategoricalPredictors, 1);
%! B = TreeBagger (3, X, yb, 'CategoricalPredictors', 'all');
%! assert_equal (B.CategoricalPredictors, [1, 2]);

%!test  # a categorical predictor may be named rather than indexed
%! B = TreeBagger (3, X, yb, 'PredictorNames', {'grp', 'val'}, ...
%!                 'CategoricalPredictors', {'grp'});
%! assert_equal (B.CategoricalPredictors, 1);
%! assert_equal (B.Trees{1}.CategoricalPredictors, 1);

%!error<TreeBagger: 'CategoricalPredictors' indices must not exceed the number of predictors.> ...
%! TreeBagger (3, X, yb, 'CategoricalPredictors', 3)

%!error<TreeBagger.append: the two ensembles must be fitted on the same observations and predictors.> ...
%! append (TreeBagger (3, X, yb, 'CategoricalPredictors', 1), ...
%!         TreeBagger (3, X, yb))

## Table input
%!test  # the predictors and the response come from a table
%! load fisheriris
%! T = table (meas(:,1), meas(:,2), 'VariableNames', {'SL', 'SW'});
%! T.Wide = categorical (meas(:,2) > 3, [false true], {'narrow', 'wide'});
%! T.Species = categorical (species);
%! B = TreeBagger (20, T, 'Species');
%! assert_equal (B.PredictorNames, {'SL', 'SW', 'Wide'});
%! assert_equal (B.CategoricalPredictors, 3);

%!test  # a model formula names the response and the predictors together
%! load fisheriris
%! T = table (meas(:,1), meas(:,2), 'VariableNames', {'SL', 'SW'});
%! T.Species = categorical (species);
%! B = TreeBagger (20, T, 'Species ~ SW + SL');
%! assert_equal (B.PredictorNames, {'SW', 'SL'});

%!test  # predict takes a table, matched by name and not by position
%! load fisheriris
%! T = table (meas(:,1), meas(:,2), 'VariableNames', {'SL', 'SW'});
%! T.Species = categorical (species);
%! B = TreeBagger (20, T, 'Species');
%! a = predict (B, T);
%! assert_equal (predict (B, T(:, [3, 2, 1])), a);

%!error<TreeBagger: the table holds no variable 'NoSuch'.> ...
%! TreeBagger (10, table (rand (6, 1), rand (6, 1)), 'NoSuch')

## A table at margin
%!test  # the response is named, left out, or given beside the table
%! load fisheriris
%! X = meas(:,1:2);
%! y = categorical (species);
%! T = table (X(:,1), X(:,2), 'VariableNames', {'SL', 'SW'});
%! T.Species = y;
%! Mdl = TreeBagger (20, T, 'Species');
%! a = margin (Mdl, X, y);
%! assert_equal (margin (Mdl, T(:,1:2), y), a);
%! assert_equal (margin (Mdl, T, 'Species'), a);
%! assert_equal (margin (Mdl, T), a);

## A table at error, meanMargin and quantileError
%!test  # the response is named, left out, or given beside the table
%! load fisheriris
%! X = meas(:,1:2);
%! y = categorical (species);
%! T = table (X(:,1), X(:,2), 'VariableNames', {'SL', 'SW'});
%! T.Species = y;
%! Mdl = TreeBagger (20, T, 'Species');
%! a = error (Mdl, X, y);
%! assert_equal (error (Mdl, T(:,1:2), y), a);
%! assert_equal (error (Mdl, T, 'Species'), a);
%! assert_equal (error (Mdl, T), a);
%! assert_equal (error (Mdl, T, 'Mode', 'ensemble'), a(end));
%!test  # the response is named, left out, or given beside the table
%! load fisheriris
%! X = meas(:,1:2);
%! y = categorical (species);
%! T = table (X(:,1), X(:,2), 'VariableNames', {'SL', 'SW'});
%! T.Species = y;
%! Mdl = TreeBagger (20, T, 'Species');
%! a = meanMargin (Mdl, X, y);
%! assert_equal (meanMargin (Mdl, T(:,1:2), y), a);
%! assert_equal (meanMargin (Mdl, T, 'Species'), a);
%! assert_equal (meanMargin (Mdl, T), a);
%!test  # the response is named, left out, or given beside the table
%! load fisheriris
%! X = meas(:,2:4);
%! y = meas(:,1);
%! T = table (X(:,1), X(:,2), X(:,3), y, ...
%!            'VariableNames', {'SW', 'PL', 'PW', 'SL'});
%! Mdl = TreeBagger (20, T, 'SL', 'Method', 'regression');
%! a = quantileError (Mdl, X, y);
%! assert_equal (quantileError (Mdl, T(:,1:3), y), a);
%! assert_equal (quantileError (Mdl, T, 'SL'), a);
%! assert_equal (quantileError (Mdl, T), a);
%! assert_equal (quantileError (Mdl, T(:,[3, 1, 2, 4])), a);
