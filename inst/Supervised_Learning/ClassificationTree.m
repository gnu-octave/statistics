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

classdef ClassificationTree
  ## -*- texinfo -*-
  ## @deftp {statistics} ClassificationTree
  ##
  ## Binary decision tree for classification
  ##
  ## The @code{ClassificationTree} class implements a CART binary decision
  ## tree.  Growth splits each node on the single predictor and cut point
  ## that lower the impurity of the response the most, and stops when a node
  ## is pure, too small to be a parent, or has no split leaving enough
  ## observations on both sides.  The grown tree is then optionally reduced,
  ## first by merging the leaves that buy no accuracy and then by cost
  ## complexity pruning, which orders the branch nodes by how little risk
  ## their subtrees remove and records that order so a subtree of any size
  ## can be recovered afterwards with @code{prune}.
  ##
  ## Create a @code{ClassificationTree} object by using the @code{fitctree}
  ## function or the class constructor.
  ##
  ## The fit is carried out by the compiled engine @code{treetrain} and
  ## predictions by @code{treepredict}, which the regression tree shares.
  ##
  ## An observation missing the predictor a node cuts on descends to neither
  ## child.  It is counted in that node and in every node above it, and
  ## @code{predict} stops it there and gives it that node's answer, so a row
  ## is never sent down a branch on evidence it does not carry.
  ##
  ## @strong{What this class does not do yet.}  Categorical predictors,
  ## surrogate splits, predictor subsampling and the @qcode{'twoing'} split
  ## criterion are not implemented, and an option asking for one of them is
  ## refused rather than quietly ignored.  @code{CategoricalSplit},
  ## @code{CutCategories} and the six @code{Surrogate} properties are
  ## therefore always empty, as they are in MATLAB on numeric data.
  ##
  ## @seealso{fitctree, treetrain, treepredict}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix holding the predictor data the model was fitted on.
    ## Each column is one predictor and each row one observation.  This
    ## property is read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} Y
    ##
    ## Class labels
    ##
    ## A logical or numeric column vector, a character array, or a cell array
    ## of character vectors with one row per row of @var{X}, holding the
    ## observed class label of each observation.  This property is read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer, the number of observations the model was fitted
    ## on.  It counts the rows kept, so it is smaller than the number of rows
    ## given whenever a response was missing.  This property is read-only.
    ##
    ## @end deftp
    NumObservations = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} RowsUsed
    ##
    ## Rows used for fitting
    ##
    ## A logical column vector with one element per row of the predictor data
    ## as it was given, true for each row used for fitting.  It is empty,
    ## @qcode{[]}, when every row was used, so a non-empty value means that
    ## rows were dropped.  Only a missing response drops a row; a row missing
    ## some of its predictors is kept.  This property is read-only.
    ##
    ## @end deftp
    RowsUsed = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector of the weights the fit used, one per retained
    ## observation.  They are the weights given, scaled so that the
    ## observations of each class sum to that class's prior, and so that all
    ## of them together sum to one.  This property is read-only.
    ##
    ## @end deftp
    W = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} PredictorNames
    ##
    ## Names of the predictor variables
    ##
    ## A cell array of character vectors with one name per column of
    ## @var{X}.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector naming the response.  This property is read-only.
    ##
    ## @end deftp
    ResponseName = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} ClassNames
    ##
    ## Names of the classes
    ##
    ## The distinct class labels, in the type the response was given in and
    ## sorted.  This property is read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A row vector of column indices into @var{X}, naming the predictors
    ## treated as categorical.  Categorical predictors are not implemented,
    ## so this is always empty.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} ExpandedPredictorNames
    ##
    ## Expanded predictor names
    ##
    ## A cell array of character vectors.  It differs from
    ## @code{PredictorNames} only when a categorical predictor has been
    ## expanded into one column per level, which this class does not do, so
    ## the two are always equal.  This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## A cell array with one column vector of bin edges per predictor, empty
    ## unless the predictors were binned before fitting.  Binning is not
    ## implemented, so this is always empty.  This property is read-only.
    ##
    ## @end deftp
    BinEdges = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} ModelParameters
    ##
    ## Parameters the fit was run with
    ##
    ## A structure recording the options the tree was grown under:
    ## @qcode{SplitCriterion}, @qcode{MinParent}, @qcode{MinLeaf},
    ## @qcode{MaxSplits}, @qcode{NVarToSample}, @qcode{MergeLeaves},
    ## @qcode{Prune}, @qcode{PruneCriterion}, @qcode{QEToler},
    ## @qcode{NSurrogate}, @qcode{MaxCat}, @qcode{AlgCat},
    ## @qcode{PredictorSelection}, @qcode{Method} and @qcode{Type}.
    ##
    ## @qcode{MinParent} is the value the fit used, which is
    ## @code{max (MinParentSize, 2 * MinLeafSize)} and so may exceed the
    ## @qcode{'MinParentSize'} asked for.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} HyperparameterOptimizationResults
    ##
    ## Results of a hyperparameter optimization
    ##
    ## Hyperparameter optimization is not implemented, so this is always
    ## empty.  This property is read-only.
    ##
    ## @end deftp
    HyperparameterOptimizationResults = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} NumNodes
    ##
    ## Number of nodes in the tree
    ##
    ## A positive integer, the number of nodes the tree holds, branch nodes
    ## and leaves together.  Nodes are numbered as they are created, so a
    ## parent always carries a lower number than either of its children.
    ## This property is read-only.
    ##
    ## @end deftp
    NumNodes = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} Children
    ##
    ## Child nodes of each node
    ##
    ## A @math{NumNodesx2} matrix naming the left and the right child of each
    ## node.  A leaf carries a zero in both columns.  This property is
    ## read-only.
    ##
    ## @end deftp
    Children = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} Parent
    ##
    ## Parent of each node
    ##
    ## A column vector naming the parent of each node.  The root carries a
    ## zero.  This property is read-only.
    ##
    ## @end deftp
    Parent = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} IsBranchNode
    ##
    ## Which nodes are branch nodes
    ##
    ## A logical column vector, true for each node that carries a split and
    ## false for each leaf.  This property is read-only.
    ##
    ## @end deftp
    IsBranchNode = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} CutPredictor
    ##
    ## Name of the predictor each node cuts on
    ##
    ## A cell array of character vectors with one entry per node, holding the
    ## name of the predictor the node splits on and an empty character vector
    ## at a leaf.  This property is read-only.
    ##
    ## @end deftp
    CutPredictor = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} CutPredictorIndex
    ##
    ## Index of the predictor each node cuts on
    ##
    ## A column vector holding, for each node, the column of @var{X} the node
    ## splits on, and zero at a leaf.  This property is read-only.
    ##
    ## @end deftp
    CutPredictorIndex = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} CutPoint
    ##
    ## Cut point of each node
    ##
    ## A column vector holding, for each node, the value the split compares
    ## the predictor against: an observation goes left when its value is less
    ## than the cut point and right otherwise.  A leaf carries @qcode{NaN}.
    ## This property is read-only.
    ##
    ## @end deftp
    CutPoint = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} CutType
    ##
    ## Type of cut at each node
    ##
    ## A cell array of character vectors holding @qcode{'continuous'} at a
    ## branch node that cuts a numeric predictor at a point,
    ## @qcode{'categorical'} at one that splits a set of levels, and an empty
    ## character vector at a leaf.  Categorical predictors are not
    ## implemented, so every branch node is @qcode{'continuous'}.  This
    ## property is read-only.
    ##
    ## @end deftp
    CutType = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} CutCategories
    ##
    ## Categories used at each branch
    ##
    ## A @math{NumNodesx2} cell array holding, for a node that cuts a
    ## categorical predictor, the levels sent left and the levels sent right.
    ## Categorical predictors are not implemented, so every entry is empty.
    ## This property is read-only.
    ##
    ## @end deftp
    CutCategories = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} CategoricalSplit
    ##
    ## Categorical splits of the tree
    ##
    ## A @math{Nx2} cell array with one row per categorical split.
    ## Categorical predictors are not implemented, so this is always empty.
    ## This property is read-only.
    ##
    ## @end deftp
    CategoricalSplit = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} NodeSize
    ##
    ## Number of observations at each node
    ##
    ## A column vector holding how many training observations reached each
    ## node.  A row missing the predictor its node cuts on is counted at that
    ## node and at none below it, so a parent's size is not in general the
    ## sum of its children's.  This property is read-only.
    ##
    ## @end deftp
    NodeSize = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} NodeClass
    ##
    ## Class assigned to each node
    ##
    ## A cell array of character vectors naming, for each node, the class of
    ## least expected misclassification cost given the node's class
    ## probabilities.  Under the default cost that is simply the most
    ## probable class, with the first of the class names kept on a tie.  This
    ## property is read-only.
    ##
    ## @end deftp
    NodeClass = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} NodeError
    ##
    ## Misclassification cost of each node
    ##
    ## A column vector holding, for each node, the expected misclassification
    ## cost of the class the node is assigned.  Under the default cost that
    ## is the probability that the node's class is wrong, one less the
    ## largest class probability.  This property is read-only.
    ##
    ## @end deftp
    NodeError = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} NodeProbability
    ##
    ## Probability of reaching each node
    ##
    ## A column vector holding, for each node, the total weight of the
    ## observations that reached it, the weights being those in @code{W}.
    ## The root carries one.  This property is read-only.
    ##
    ## @end deftp
    NodeProbability = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} NodeRisk
    ##
    ## Risk of each node
    ##
    ## A column vector holding, for each node, the impurity of the node
    ## weighted by the probability of reaching it, measured by whichever of
    ## @qcode{'gdi'} and @qcode{'deviance'} the tree was grown under.
    ##
    ## A non-default @code{Cost} enters here rather than through the class
    ## probabilities: the weights are scaled class by class by the total cost
    ## of misclassifying that class, and the impurity is measured on the
    ## scaled distribution.  This property is read-only.
    ##
    ## @end deftp
    NodeRisk = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} ClassCount
    ##
    ## Class counts at each node
    ##
    ## A @math{NumNodesxK} matrix holding how many training observations of
    ## each class reached each node.  These are counts and take no notice of
    ## the observation weights.  This property is read-only.
    ##
    ## @end deftp
    ClassCount = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} ClassProbability
    ##
    ## Class probabilities at each node
    ##
    ## A @math{NumNodesxK} matrix holding, for each node, the weight of each
    ## class among the observations that reached it, as a proportion of the
    ## node's total weight.  The root row is the prior.  This property is
    ## read-only.
    ##
    ## @end deftp
    ClassProbability = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} PruneList
    ##
    ## Pruning level of each node
    ##
    ## A column vector holding, for each branch node, the level of the cost
    ## complexity sequence at which it stops being a branch node, and zero at
    ## a leaf.  Pruning the tree to level @var{L} turns every node whose
    ## level is between one and @var{L} into a leaf.  It is empty when
    ## neither @qcode{'Prune'} nor @qcode{'MergeLeaves'} was asked for, since
    ## no sequence was then estimated.  This property is read-only.
    ##
    ## @end deftp
    PruneList = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} PruneAlpha
    ##
    ## Cost complexity parameter of each pruning level
    ##
    ## A column vector with one element per level of the pruning sequence,
    ## the first of which is zero and stands for the unpruned tree.  Level
    ## @var{L} is the smallest subtree that is optimal for every complexity
    ## parameter from @code{PruneAlpha(@var{L}+1)} up to the next one.  This
    ## property is read-only.
    ##
    ## @end deftp
    PruneAlpha = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} SurrogateCutCategories
    ##
    ## Categories of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutCategories = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} SurrogateCutFlip
    ##
    ## Cut assignments of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutFlip = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} SurrogateCutPoint
    ##
    ## Cut points of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutPoint = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} SurrogateCutType
    ##
    ## Types of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutType = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} SurrogateCutPredictor
    ##
    ## Predictors of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutPredictor = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} SurrogatePredictorAssociation
    ##
    ## Predictive measures of association of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogatePredictorAssociation = {};

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} Prior
    ##
    ## Prior probability of each class
    ##
    ## A numeric row vector with one element per class, summing to one.  It
    ## defaults to the weight each class carries in the training data, and
    ## may be reassigned after fitting.
    ##
    ## Reassigning it re-derives @code{W} and every node statistic that
    ## depends on the class weights, so @code{ClassProbability},
    ## @code{NodeProbability}, @code{NodeClass}, @code{NodeError} and
    ## @code{NodeRisk} all follow.  The shape of the tree does not, having
    ## been decided by the prior in force when it was grown.
    ##
    ## @end deftp
    Prior = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} Cost
    ##
    ## Misclassification cost
    ##
    ## A square numeric matrix with one row and column per class, where
    ## @code{Cost(i,j)} is the cost of classifying an observation of class
    ## @var{i} into class @var{j}.  It defaults to @code{1 - eye (K)} and may
    ## be reassigned after fitting.
    ##
    ## Reassigning it re-derives @code{NodeClass}, @code{NodeError} and
    ## @code{NodeRisk}, and changes what @code{predict} answers.  The shape
    ## of the tree does not follow, having been decided by the cost in force
    ## when it was grown.
    ##
    ## @end deftp
    Cost = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationTree} {property} ScoreTransform
    ##
    ## Transform applied to the scores
    ##
    ## A character vector naming the function @code{predict} applies to the
    ## class probabilities before returning them, or a function handle.  The
    ## default is @qcode{'none'}.
    ##
    ## @end deftp
    ScoreTransform = [];

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## The parsed ScoreTransform, applied to the scores by predict.
    STfun = [];

    ## The share of each class's total weight that reached each node, one
    ## row per node and one column per class.  It is what the node
    ## statistics are re-derived from, and it is free of both the prior and
    ## the cost: multiplying a class's column by a constant leaves it
    ## unchanged, so the same table serves whatever Prior and Cost are
    ## assigned afterwards.
    ClassShare = [];

    ## The risk each branch node carries on account of the observations that
    ## stop there, missing the predictor it cuts on.  Those observations are
    ## in neither child, so a subtree's risk is its children's plus this.  It
    ## is zero at a leaf and zero throughout a tree fitted on data with no
    ## missing values.
    HeldRisk = [];

    ## The same, on the scale NodeRisk lives on, which is the impurity of the
    ## cost-adjusted distribution.  The two part company under a non-default
    ## cost, and the pruning sequence wants the first while
    ## predictorImportance wants this one.
    HeldImpurity = [];

    ## The probability of reaching each node on the cost-adjusted weights.
    NodeProbAdj = [];

    ## The observation weights as they were given, before the prior scaled
    ## them.  Reassigning Prior re-derives W from these, which reassigning it
    ## from W could not do: W has already had a prior divided into it.
    RawWeights = [];

  endproperties

  ## Set methods for the properties a user may assign after fitting.
  methods (Hidden)

    function this = set.Cost (this, Cost)
      [C, errmsg] = costMatrix (Cost, this.ClassNames);
      if (! isempty (errmsg))
        error ("ClassificationTree.Cost: %s", errmsg);
      endif
      this.Cost = C;
      this = deriveNodes (this);
    endfunction

    function this = set.Prior (this, Prior)
      P = Prior;
      K = classCount (this.ClassNames);
      if (isstruct (P))
        P = priorFromStruct (P, this.ClassNames, 'ClassificationTree.Prior');
      elseif (ischar (P))
        if (strcmpi (P, 'uniform'))
          P = ones (1, K) / K;
        elseif (! strcmpi (P, 'empirical'))
          error (strcat ("ClassificationTree.Prior: a character vector", ...
                         " must be 'empirical' or 'uniform'."));
        else
          return;   # 'empirical' after fitting is what is already stored
        endif
      endif
      if (! (isnumeric (P) && isvector (P) && isreal (P)))
        error (strcat ("ClassificationTree.Prior: must be a real numeric", ...
                       " vector, a structure, 'empirical', or 'uniform'."));
      endif
      if (numel (P) != K)
        error (strcat ("ClassificationTree.Prior: must have one element", ...
                       " per class."));
      endif
      if (any (P < 0) || ! (sum (P) > 0))
        error (strcat ("ClassificationTree.Prior: must be nonnegative and", ...
                       " must not be all zero."));
      endif
      this.Prior = P(:)' / sum (P);
      ## The weights follow the prior, so reassigning one re-derives the
      ## other, and the node statistics follow both.
      if (! isempty (this.RawWeights))
        gY = labelIndices (this.ClassNames, this.Y);
        this.W = priorNormalize (this.RawWeights, gY, this.Prior);
      endif
      this = deriveNodes (this);
    endfunction

    function this = set.ScoreTransform (this, val)
      [f, st] = parseScoreTransform (val, 'ClassificationTree');
      this.ScoreTransform = st;
      this.STfun = f;
    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ('\n  ClassificationTree\n\n');
      fprintf ('%22s: %s\n', 'ResponseName', this.ResponseName);
      fprintf ('%22s: %s\n', 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ('%22s: %s\n', 'ClassNames', classNameListing (this.ClassNames));
      fprintf ('%22s: %s\n', 'ScoreTransform', this.ScoreTransform);
      fprintf ('%22s: %d\n', 'NumObservations', this.NumObservations);
      fprintf ('%22s: %d\n', 'NumNodes', this.NumNodes);
      fprintf ('\n');
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationTree} {@var{obj} =} ClassificationTree (@var{X}, @var{Y})
    ## @deftypefnx {ClassificationTree} {@var{obj} =} ClassificationTree (@dots{}, @var{name}, @var{value})
    ##
    ## Grow a binary decision tree for classification.
    ##
    ## @code{@var{obj} = ClassificationTree (@var{X}, @var{Y})} grows a tree
    ## on the @math{NxP} numeric matrix @var{X} of predictor data and the
    ## @math{Nx1} response @var{Y}, and returns it as a
    ## @code{ClassificationTree} object.  @var{Y} may be a numeric or logical
    ## vector, a character array, or a cell array of character vectors, and
    ## the class names come back in the type it was given in.
    ##
    ## @code{@var{obj} = ClassificationTree (@dots{}, @var{name},
    ## @var{value})} takes the options below.
    ##
    ## @multitable @columnfractions 0.20 0.78
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'ClassNames'} @tab The classes to fit, of the same type
    ## as @var{Y}.  Observations of any other class are dropped.
    ##
    ## @item @qcode{'Cost'} @tab A square matrix with one row and column per
    ## class, where element @math{(i,j)} is the cost of classifying an
    ## observation of class @math{i} into class @math{j}, or a structure with
    ## fields @qcode{ClassNames} and @qcode{ClassificationCosts}.  The
    ## default is @code{1 - eye (K)}.
    ##
    ## @item @qcode{'MaxNumSplits'} @tab A nonnegative integer, the largest
    ## number of branch nodes the tree may take.  The default is one less
    ## than the number of observations, which is as many as a tree can have.
    ##
    ## @item @qcode{'MergeLeaves'} @tab @qcode{'on'} (default) or
    ## @qcode{'off'}.  When on, a pair of leaves whose parent is no worse
    ## than the two of them together is merged back into that parent.
    ##
    ## @item @qcode{'MinLeafSize'} @tab A positive integer, the fewest
    ## observations a leaf may hold.  The default is 1.  A split leaving
    ## fewer than this on either side is not taken.
    ##
    ## @item @qcode{'MinParentSize'} @tab A positive integer, the fewest
    ## observations a node must hold to be split at all.  The default is 10.
    ## The value the fit uses is @code{max (MinParentSize, 2 * MinLeafSize)},
    ## since a smaller node cannot give both children a legal leaf.
    ##
    ## @item @qcode{'PredictorNames'} @tab A cell array of character vectors
    ## naming the columns of @var{X}.
    ##
    ## @item @qcode{'Prior'} @tab @qcode{'empirical'} (default),
    ## @qcode{'uniform'}, a numeric vector with one element per class, or a
    ## structure with fields @qcode{ClassNames} and
    ## @qcode{ClassProbs}.
    ##
    ## @item @qcode{'Prune'} @tab @qcode{'on'} (default) or @qcode{'off'}.
    ## When on, the cost complexity pruning sequence is estimated and
    ## reported in @code{PruneList} and @code{PruneAlpha}.  The tree returned
    ## is the unpruned one either way; @code{prune} takes a subtree out of
    ## the sequence.
    ##
    ## @item @qcode{'PruneCriterion'} @tab @qcode{'error'}, the only
    ## criterion implemented.
    ##
    ## @item @qcode{'ResponseName'} @tab A character vector naming the
    ## response.  The default is @qcode{'Y'}.
    ##
    ## @item @qcode{'ScoreTransform'} @tab A character vector naming a
    ## transform to apply to the scores, or a function handle.  The default
    ## is @qcode{'none'}.
    ##
    ## @item @qcode{'SplitCriterion'} @tab @qcode{'gdi'} (default), the Gini
    ## diversity index, or @qcode{'deviance'}, the cross entropy.
    ##
    ## @item @qcode{'Weights'} @tab A nonnegative numeric vector with one
    ## element per observation.  The default is uniform.
    ##
    ## @end multitable
    ##
    ## @seealso{fitctree, treetrain, treepredict}
    ## @end deftypefn
    function this = ClassificationTree (X, Y, varargin)

      ## Input validation
      if (nargin < 2)
        error ("ClassificationTree: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationTree: name-value arguments must be", ...
                       " in pairs."));
      endif
      if (! (isnumeric (X) && isreal (X) && ismatrix (X) && ! isempty (X)))
        error (strcat ("ClassificationTree: X must be a non-empty real", ...
                       " numeric matrix."));
      endif
      if (rows (X) != rows (Y))
        error (strcat ("ClassificationTree: number of rows in X and Y", ...
                       " must be equal."));
      endif

      ## Groups in Y, which is how a missing response is found: grp2idx
      ## leaves one NaN, whatever type the response is given in.
      [gY, gnY, glY] = grp2idx (Y);

      ## Defaults.  MaxNumSplits is left empty until the retained rows are
      ## known, its default being one less than their number.
      PredictorNames = {};
      ResponseName   = [];
      ClassNames     = [];
      Prior          = 'empirical';
      Cost           = [];
      Weights        = [];
      MaxNumSplits   = [];
      MergeLeaves    = 'on';
      MinLeafSize    = 1;
      MinParentSize  = 10;
      Prune          = 'on';
      PruneCriterion = 'error';
      SplitCriterion = 'gdi';
      this.ScoreTransform = 'none';

      ## Parse optional parameters
      while (numel (varargin) > 0)
        Value = varargin{2};
        switch (tolower (varargin{1}))

          case 'predictornames'
            PredictorNames = Value;
            if (! iscellstr (PredictorNames))
              error (strcat ("ClassificationTree: 'PredictorNames' must", ...
                             " be supplied as a cellstring array."));
            elseif (numel (PredictorNames) != columns (X))
              error (strcat ("ClassificationTree: 'PredictorNames' must", ...
                             " equal the number of columns in X."));
            endif

          case 'responsename'
            ResponseName = Value;
            if (! (ischar (ResponseName) && isrow (ResponseName)))
              error (strcat ("ClassificationTree: 'ResponseName' must be", ...
                             " a character vector."));
            endif

          case 'classnames'
            ClassNames = Value;
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames)
                   || islogical (ClassNames) || ischar (ClassNames)))
              error (strcat ("ClassificationTree: 'ClassNames' must be a", ...
                             " cell array of character vectors, a logical", ...
                             " vector, a numeric vector, or a character", ...
                             " array."));
            endif
            if (iscellstr (ClassNames) || ischar (ClassNames))
              known = ismember (cellstr (ClassNames), gnY);
            else
              known = ismember (ClassNames(:), glY(:));
            endif
            if (! all (known))
              error (strcat ("ClassificationTree: not all 'ClassNames'", ...
                             " are present in Y."));
            endif

          case 'prior'
            Prior = Value;
            if (! (isstruct (Prior) || isnumeric (Prior) || ischar (Prior)))
              error (strcat ("ClassificationTree: 'Prior' must be a", ...
                             " numeric vector, a structure, or a", ...
                             " character vector."));
            endif
            if (ischar (Prior)
                && ! any (strcmpi (Prior, {'empirical', 'uniform'})))
              error (strcat ("ClassificationTree: 'Prior' must be", ...
                             " 'empirical', 'uniform', a numeric vector,", ...
                             " or a structure."));
            endif

          case 'cost'
            Cost = Value;
            if (! (isstruct (Cost)
                   || (isnumeric (Cost) && issquare (Cost))))
              error (strcat ("ClassificationTree: 'Cost' must be a", ...
                             " numeric square matrix or a structure."));
            endif

          case 'weights'
            Weights = Value;
            if (! (isnumeric (Weights) && isvector (Weights)
                   && isreal (Weights)))
              error (strcat ("ClassificationTree: 'Weights' must be a", ...
                             " real numeric vector."));
            endif
            if (numel (Weights) != rows (X))
              error (strcat ("ClassificationTree: 'Weights' must have one", ...
                             " element per row in X."));
            endif
            if (any (Weights < 0) || ! (sum (Weights) > 0))
              error (strcat ("ClassificationTree: 'Weights' must be", ...
                             " nonnegative and must not be all zero."));
            endif

          case 'scoretransform'
            this.ScoreTransform = Value;

          case 'maxnumsplits'
            MaxNumSplits = Value;
            if (! (isnumeric (MaxNumSplits) && isscalar (MaxNumSplits)
                   && isreal (MaxNumSplits) && MaxNumSplits >= 0
                   && MaxNumSplits == fix (MaxNumSplits)))
              error (strcat ("ClassificationTree: 'MaxNumSplits' must be", ...
                             " a nonnegative integer."));
            endif

          case 'minleafsize'
            MinLeafSize = Value;
            if (! (isnumeric (MinLeafSize) && isscalar (MinLeafSize)
                   && isreal (MinLeafSize) && MinLeafSize >= 1
                   && MinLeafSize == fix (MinLeafSize)))
              error (strcat ("ClassificationTree: 'MinLeafSize' must be a", ...
                             " positive integer."));
            endif

          case 'minparentsize'
            MinParentSize = Value;
            if (! (isnumeric (MinParentSize) && isscalar (MinParentSize)
                   && isreal (MinParentSize) && MinParentSize >= 1
                   && MinParentSize == fix (MinParentSize)))
              error (strcat ("ClassificationTree: 'MinParentSize' must be", ...
                             " a positive integer."));
            endif

          case 'mergeleaves'
            MergeLeaves = Value;
            if (! (ischar (MergeLeaves)
                   && any (strcmpi (MergeLeaves, {'on', 'off'}))))
              error (strcat ("ClassificationTree: 'MergeLeaves' must be", ...
                             " either 'on' or 'off'."));
            endif

          case 'prune'
            Prune = Value;
            if (! (ischar (Prune) && any (strcmpi (Prune, {'on', 'off'}))))
              error (strcat ("ClassificationTree: 'Prune' must be either", ...
                             " 'on' or 'off'."));
            endif

          case 'prunecriterion'
            PruneCriterion = Value;
            if (! (ischar (PruneCriterion)
                   && any (strcmpi (PruneCriterion, {'error', 'impurity'}))))
              error (strcat ("ClassificationTree: 'PruneCriterion' must", ...
                             " be either 'error' or 'impurity'."));
            endif
            if (strcmpi (PruneCriterion, 'impurity'))
              error (strcat ("ClassificationTree: 'PruneCriterion'", ...
                             " 'impurity' is not implemented."));
            endif

          case 'splitcriterion'
            SplitCriterion = Value;
            if (! (ischar (SplitCriterion)
                   && any (strcmpi (SplitCriterion, {'gdi', 'deviance', ...
                                                     'twoing'}))))
              error (strcat ("ClassificationTree: 'SplitCriterion' must", ...
                             " be 'gdi', 'deviance', or 'twoing'."));
            endif
            if (strcmpi (SplitCriterion, 'twoing'))
              error (strcat ("ClassificationTree: 'SplitCriterion'", ...
                             " 'twoing' is not implemented."));
            endif

          case 'categoricalpredictors'
            ## Accepted only when it asks for nothing, so that a caller
            ## passing the MATLAB default is not turned away.
            if (! (isempty (Value)
                   || (ischar (Value) && strcmpi (Value, 'none'))))
              error (strcat ("ClassificationTree: categorical predictors", ...
                             " are not implemented."));
            endif

          ## Options MATLAB takes that this class does not implement.  They
          ## are named one by one so that asking for one is refused rather
          ## than quietly doing nothing.
          case {'surrogate', 'numvariablestosample', 'predictorselection', ...
                'algorithmforcategorical', 'maxnumcategories', 'numbins', ...
                'optimizehyperparameters', ...
                'hyperparameteroptimizationoptions'}
            error ("ClassificationTree: '%s' is not implemented.", ...
                   varargin{1});

          case {'crossval', 'cvpartition', 'holdout', 'kfold', 'leaveout'}
            error (strcat ("ClassificationTree: '%s' is not implemented;", ...
                           " fit the model and cross-validate it", ...
                           " afterwards."), varargin{1});

          otherwise
            error (strcat ("ClassificationTree: invalid parameter name in", ...
                           " optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Default predictor and response names
      if (isempty (PredictorNames))
        PredictorNames = cell (1, columns (X));
        for ii = 1:columns (X)
          PredictorNames{ii} = sprintf ("x%d", ii);
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = 'Y';
      endif
      this.PredictorNames = PredictorNames;
      this.ResponseName = ResponseName;
      this.CategoricalPredictors = [];
      this.ExpandedPredictorNames = PredictorNames;

      ## A class the caller did not ask for takes its observations with it.
      ## Anything textual is matched as whole names, gnY being grp2idx's own
      ## cellstr of them: a character matrix is not a cellstr, and ismember
      ## between two of them compares character by character.
      if (! isempty (ClassNames))
        if (iscellstr (ClassNames) || ischar (ClassNames))
          drop = find (! ismember (gnY, cellstr (ClassNames)));
        else
          drop = find (! ismember (glY, ClassNames));
        endif
        for ii = 1:numel (drop)
          gY(gY == drop(ii)) = NaN;
        endfor
      endif

      ## An observation is dropped only when its response is missing.  A row
      ## whose predictors hold missing values is kept and reported as used;
      ## the fit sends it as far down the tree as the predictors it does
      ## carry allow.
      RowsUsed = ! isnan (gY);
      ## Index the rows and not the elements: a response naming its classes
      ## in the rows of a character matrix has one column per character, and
      ## a linear index would flatten the names into single letters.
      Y = Y(RowsUsed, :);
      X = X(RowsUsed, :);
      if (! isempty (Weights))
        Weights = Weights(RowsUsed);
      endif
      if (isempty (Y))
        error ("ClassificationTree: no observations with a known class.");
      endif

      this.X = X;
      this.Y = Y;
      this.NumObservations = rows (X);
      ## RowsUsed is left empty when every observation was used, as in MATLAB
      if (all (RowsUsed))
        this.RowsUsed = [];
      else
        this.RowsUsed = RowsUsed;
      endif

      ## The classes, sorted and in the type the response was given in
      [this.ClassNames, ~, gY] = uniqueLabels (Y);
      gY = gY(:);
      K = classCount (this.ClassNames);

      ## The raw weights, kept so that reassigning Prior can re-derive W
      if (isempty (Weights))
        RawWeights = ones (this.NumObservations, 1);
      else
        RawWeights = double (Weights(:));
      endif
      this.RawWeights = RawWeights;

      ## Cost first: set.Prior counts the classes and set.Cost does not
      ## depend on the prior, and both are wanted before the weights.
      if (isempty (Cost))
        Cost = 1 - eye (K);
      endif
      this.Cost = Cost;

      ## The prior, empirical over the weights unless it was given
      if (ischar (Prior) && strcmpi (Prior, 'empirical'))
        P = zeros (1, K);
        for k = 1:K
          P(k) = sum (RawWeights(gY == k));
        endfor
        if (! (sum (P) > 0))
          error (strcat ("ClassificationTree: 'Weights' must not be zero", ...
                         " for every observation."));
        endif
        this.Prior = P / sum (P);
      else
        this.Prior = Prior;
      endif

      ## Weights scaled so that each class carries its prior and the whole
      ## sums to one, which is the W MATLAB reports.
      this.W = priorNormalize (RawWeights, gY, this.Prior);

      ## The engine grows the tree on the cost-adjusted weights.  Scaling a
      ## class by the total cost of misclassifying it is the classical way a
      ## cost matrix enters CART, and it is what MATLAB does: the split
      ## criterion then needs no notion of cost, while everything reported
      ## below is measured on the unadjusted weights.
      cAdj = sum (this.Cost, 2)';
      wAdj = this.W .* cAdj(gY)';
      if (sum (wAdj) > 0)
        wAdj = wAdj / sum (wAdj);
      else
        ## Nothing costs anything, which is what a single class means: its
        ## cost matrix is the one by one zero.  There is no misclassification
        ## to weigh, so the fit runs on the weights as they stand.
        wAdj = this.W;
      endif

      if (isempty (MaxNumSplits))
        MaxNumSplits = max (this.NumObservations - 1, 0);
      endif
      ## A node smaller than two leaves cannot be split whatever the parent
      ## size asked for, so the fit uses the larger of the two.
      MinParent = max (MinParentSize, 2 * MinLeafSize);
      mergeOn = strcmpi (MergeLeaves, 'on');
      pruneOn = strcmpi (Prune, 'on');

      ## The engine measures risk as misclassified weight, which is the
      ## expected misclassification cost only under the default cost.  Under
      ## any other, growth still belongs to the engine but the merge and the
      ## pruning sequence are taken over below, on the risk the cost defines.
      plainCost = isequal (this.Cost, 1 - eye (K));

      opts = struct ('NumClasses', K, ...
                     'MinParent', MinParent, ...
                     'MinLeaf', MinLeafSize, ...
                     'MaxSplits', MaxNumSplits, ...
                     'SplitCriterion', tolower (SplitCriterion), ...
                     'MergeLeaves', plainCost && mergeOn, ...
                     'Prune', false);

      T = treetrain (X, gY, wAdj, opts);

      ## The node table the engine returns
      this.NumNodes = T.NumNodes;
      this.Children = T.Children;
      this.Parent = T.Parent;
      this.CutPredictorIndex = T.CutPredictorIndex;
      this.CutPoint = T.CutPoint;
      this.IsBranchNode = logical (T.IsBranchNode);
      this.NodeSize = T.NodeSize;
      this.ClassCount = T.ClassCount;

      ## The share of each class's weight that reached each node.  Dividing
      ## by the root's row cancels both the prior and the cost adjustment,
      ## since each scales a whole column, so what is left is the table the
      ## node statistics are re-derived from whenever either is reassigned.
      this.ClassShare = classShareOf (T.ClassWeight);

      this = fillCuts (this);
      this.ModelParameters = struct ('SplitCriterion', ...
                                     tolower (SplitCriterion), ...
                                     'MinParent', MinParent, ...
                                     'MinLeaf', MinLeafSize, ...
                                     'MaxSplits', MaxNumSplits, ...
                                     'NVarToSample', 'all', ...
                                     'MergeLeaves', tolower (MergeLeaves), ...
                                     'Prune', tolower (Prune), ...
                                     'PruneCriterion', ...
                                     tolower (PruneCriterion), ...
                                     'QEToler', [], ...
                                     'NSurrogate', 0, ...
                                     'MaxCat', 10, ...
                                     'AlgCat', 'auto', ...
                                     'PredictorSelection', 'allsplits', ...
                                     'Method', 'Tree', ...
                                     'Type', 'classification');

      ## deriveNodes reads ModelParameters.SplitCriterion, so it follows the
      ## structure above rather than the node table it works on.
      this = deriveNodes (this);
      this = deriveHeld (this);

      if (! plainCost && mergeOn)
        this = mergeLeaves (this);
      endif

      ## Merging leaves is the first step of the cost complexity sequence, so
      ## a merged tree carries one whether or not pruning was asked for.
      ## Measured on R2024a: only turning both off leaves the two pruning
      ## properties empty.
      if (pruneOn || mergeOn)
        this = pruneSequence (this);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationTree} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationTree} {[@var{label}, @var{score}] =} predict (@dots{})
    ## @deftypefnx {ClassificationTree} {[@var{label}, @var{score}, @var{node}] =} predict (@dots{})
    ## @deftypefnx {ClassificationTree} {[@var{label}, @var{score}, @var{node}, @var{cnum}] =} predict (@dots{})
    ##
    ## Classify new data with a trained @code{ClassificationTree} object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} sends each row of
    ## @var{XC} down the tree and returns the class of the node it comes to
    ## rest at.  @var{XC} must have as many columns as the predictor data the
    ## model was fitted on.
    ##
    ## @code{[@var{label}, @var{score}] = predict (@dots{})} also returns
    ## @var{score}, an @math{NxK} matrix holding the class probabilities of
    ## the node each row landed in, after @code{ScoreTransform}.
    ##
    ## @code{[@var{label}, @var{score}, @var{node}] = predict (@dots{})} also
    ## returns the number of the node each row landed in, and
    ## @code{[@var{label}, @var{score}, @var{node}, @var{cnum}] = predict
    ## (@dots{})} the index of the predicted class into @code{ClassNames}.
    ##
    ## The label is the class of least expected misclassification cost, which
    ## under the default @code{Cost} is the most probable class of the node.
    ##
    ## A row missing the predictor a node cuts on is stopped at that node and
    ## takes its answer, rather than being sent down a branch on evidence the
    ## row does not carry.
    ##
    ## @seealso{ClassificationTree, fitctree}
    ## @end deftypefn
    function [label, score, node, cnum] = predict (this, XC)

      ## Input validation
      if (nargin < 2)
        error ("ClassificationTree.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("ClassificationTree.predict: XC is empty.");
      endif
      if (! (isnumeric (XC) && isreal (XC) && ismatrix (XC)))
        error (strcat ("ClassificationTree.predict: XC must be a real", ...
                       " numeric matrix."));
      endif
      if (columns (this.X) != columns (XC))
        error (strcat ("ClassificationTree.predict: XC must have the same", ...
                       " number of predictors as the trained model."));
      endif

      [score, node] = treepredict (XC, this.Children, ...
                                   this.CutPredictorIndex, this.CutPoint, ...
                                   this.ClassProbability);

      ## The class of least expected cost, ties kept by the first class
      [~, cnum] = min (score * this.Cost, [], 2);
      label = labelsFromIndex (this.ClassNames, cnum);

      ## The transform is applied once, to the assembled probabilities
      score = this.STfun (score);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationTree} {@var{obj2} =} prune (@var{obj})
    ## @deftypefnx {ClassificationTree} {@var{obj2} =} prune (@var{obj}, @qcode{'Level'}, @var{L})
    ## @deftypefnx {ClassificationTree} {@var{obj2} =} prune (@var{obj}, @qcode{'Alpha'}, @var{A})
    ## @deftypefnx {ClassificationTree} {@var{obj2} =} prune (@var{obj}, @qcode{'Nodes'}, @var{N})
    ##
    ## Take a subtree out of the pruning sequence.
    ##
    ## @code{@var{obj2} = prune (@var{obj})} returns the tree unchanged.
    ##
    ## @code{@var{obj2} = prune (@var{obj}, @qcode{'Level'}, @var{L})} turns
    ## every branch node whose @code{PruneList} level is between one and
    ## @var{L} into a leaf and discards everything below it.  Level zero is
    ## the tree itself and the largest level is the root alone.  A level
    ## above the largest prunes to the root and warns.
    ##
    ## @code{@var{obj2} = prune (@var{obj}, @qcode{'Alpha'}, @var{A})} prunes
    ## to the smallest subtree that is optimal for the cost complexity
    ## parameter @var{A}, which is the largest level whose @code{PruneAlpha}
    ## does not exceed it.
    ##
    ## @code{@var{obj2} = prune (@var{obj}, @qcode{'Nodes'}, @var{N})} turns
    ## the branch nodes named in @var{N} into leaves, along with everything
    ## below them, and leaves the rest of the tree alone.
    ##
    ## Pruning renumbers the nodes, so the properties of the returned tree
    ## are those of a tree of that shape and not a subset of the original's.
    ##
    ## @seealso{ClassificationTree, fitctree, ClassificationTree.PruneList,
    ## ClassificationTree.PruneAlpha}
    ## @end deftypefn
    function this = prune (this, varargin)

      ## Input validation
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationTree.prune: name-value arguments", ...
                       " must be in pairs."));
      endif
      if (numel (varargin) > 2)
        error (strcat ("ClassificationTree.prune: specify only one of the", ...
                       " optional name-value paired arguments."));
      endif

      nodes = [];
      if (numel (varargin) == 2)
        Value = varargin{2};
        switch (tolower (varargin{1}))

          case 'level'
            if (! (isnumeric (Value) && isscalar (Value) && isreal (Value)
                   && Value >= 0 && Value == fix (Value)))
              error (strcat ("ClassificationTree.prune: 'Level' must be a", ...
                             " nonnegative integer."));
            endif
            nodes = this.nodesAtLevel (Value, 'prune');

          case 'alpha'
            if (! (isnumeric (Value) && isscalar (Value) && isreal (Value)
                   && Value >= 0))
              error (strcat ("ClassificationTree.prune: 'Alpha' must be a", ...
                             " nonnegative scalar."));
            endif
            if (isempty (this.PruneAlpha))
              error (strcat ("ClassificationTree.prune: the tree carries", ...
                             " no pruning sequence."));
            endif
            ## The largest level whose parameter the value reaches
            L = sum (this.PruneAlpha <= Value) - 1;
            nodes = this.nodesAtLevel (max (L, 0), 'prune');

          case 'nodes'
            if (! (isnumeric (Value) && isreal (Value) && isvector (Value)
                   && all (Value >= 1) && all (Value <= this.NumNodes)
                   && all (Value == fix (Value))))
              error (strcat ("ClassificationTree.prune: 'Nodes' must hold", ...
                             " indices of nodes of the tree."));
            endif
            nodes = Value(:);

          otherwise
            error (strcat ("ClassificationTree.prune: invalid parameter", ...
                           " name in optional pair arguments."));

        endswitch
      endif

      if (isempty (nodes))
        return;
      endif
      this = collapseNodes (this, nodes);
      ## The sequence of the tree that is left, not a slice of the one the
      ## tree came with: pruning renumbers the nodes and a subtree of an
      ## optimal sequence is the optimal sequence of that subtree.
      this = pruneSequence (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationTree} {@var{E} =} cvloss (@var{obj})
    ## @deftypefnx {ClassificationTree} {[@var{E}, @var{SE}, @var{Nleaf}, @var{BestLevel}] =} cvloss (@var{obj})
    ## @deftypefnx {ClassificationTree} {[@dots{}] =} cvloss (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validated loss of a tree and of its subtrees.
    ##
    ## @code{@var{E} = cvloss (@var{obj})} partitions the training data into
    ## ten stratified folds, grows a tree on the training part of each, and
    ## returns the loss of the held-out part.
    ##
    ## @code{[@var{E}, @var{SE}, @var{Nleaf}, @var{BestLevel}] = cvloss
    ## (@dots{})} also returns @var{SE}, the standard error of @var{E} over
    ## the folds, @var{Nleaf}, the number of leaves each subtree holds, and
    ## @var{BestLevel}, the pruning level chosen by @qcode{'TreeSize'}.  Each
    ## has one element per subtree asked for.
    ##
    ## @code{[@dots{}] = cvloss (@dots{}, @var{name}, @var{value})} takes the
    ## options below.
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'SubTrees'} @tab A vector of pruning levels in ascending
    ## order, or @qcode{'all'} for every level of the sequence.  The default
    ## is 0, the unpruned tree.
    ##
    ## @item @qcode{'TreeSize'} @tab @qcode{'se'} (default), the smallest
    ## subtree whose loss is within one standard error of the smallest loss,
    ## or @qcode{'min'}, the smallest subtree of least loss.
    ##
    ## @item @qcode{'KFold'} @tab An integer greater than 1, the number of
    ## folds.  The default is 10.  A value above the number of observations
    ## is reduced to it.
    ##
    ## @end multitable
    ##
    ## A fold's tree is pruned to the level its own sequence gives for the
    ## geometric mean of the parent's two neighbouring complexity parameters,
    ## which is the classical way a fold is matched to a subtree of the whole
    ## tree.  The last level takes every fold's tree back to its root.
    ##
    ## Every fold is grown with the parent's class names, prior, cost and
    ## observation weights, and the loss is weighed by the model's own
    ## weights.
    ##
    ## @strong{The standard error is not MATLAB's.}  This is the standard
    ## error of the loss over the folds, which is what the name means and
    ## which is zero when every fold answers alike, as MATLAB's is.  Its
    ## value is not MATLAB's: on the iris tree MATLAB reports 0.019956 where
    ## the folds give 0.024343, and no formula over the fold losses, the
    ## observation losses or the loss itself reproduces MATLAB's number.
    ## @var{E}, @var{Nleaf} and @var{BestLevel} are measured and match.
    ##
    ## @seealso{ClassificationTree, ClassificationTree.prune,
    ## ClassificationTree.crossval, ClassificationTree.loss}
    ## @end deftypefn
    function [E, SE, Nleaf, BestLevel] = cvloss (this, varargin)

      ## Input validation
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationTree.cvloss: name-value arguments", ...
                       " must be in pairs."));
      endif
      if (isempty (this.PruneAlpha))
        error (strcat ("ClassificationTree.cvloss: the tree carries no", ...
                       " pruning sequence; fit it with 'Prune' or", ...
                       " 'MergeLeaves' on."));
      endif

      maxLevel = numel (this.PruneAlpha) - 1;
      SubTrees = 0;
      TreeSize = 'se';
      KFold = 10;

      while (numel (varargin) > 0)
        Value = varargin{2};
        switch (tolower (varargin{1}))

          case 'subtrees'
            if (ischar (Value) && strcmpi (Value, 'all'))
              SubTrees = 0:maxLevel;
            elseif (isnumeric (Value) && isreal (Value) && isvector (Value)
                    && ! isempty (Value) && all (Value >= 0)
                    && all (Value == fix (Value))
                    && all (diff (Value(:)') > 0))
              SubTrees = Value(:)';
            else
              error (strcat ("ClassificationTree.cvloss: 'SubTrees' must", ...
                             " be 'all' or a vector of nonnegative", ...
                             " integers in ascending order."));
            endif

          case 'treesize'
            if (! (ischar (Value) && any (strcmpi (Value, {'se', 'min'}))))
              error (strcat ("ClassificationTree.cvloss: 'TreeSize' must", ...
                             " be either 'se' or 'min'."));
            endif
            TreeSize = tolower (Value);

          case 'kfold'
            KFold = Value;
            if (! (isnumeric (KFold) && isscalar (KFold) && isreal (KFold)
                   && KFold == fix (KFold) && KFold > 1))
              error (strcat ("ClassificationTree.cvloss: 'KFold' must be", ...
                             " an integer value greater than 1."));
            endif

          otherwise
            error (strcat ("ClassificationTree.cvloss: invalid parameter", ...
                           " name in optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      if (any (SubTrees > maxLevel))
        error (strcat ("ClassificationTree.cvloss: 'SubTrees' must not", ...
                       " exceed the largest pruning level, %d."), maxLevel);
      endif
      if (KFold > this.NumObservations)
        warning (strcat ("ClassificationTree.cvloss: 'KFold' is greater", ...
                         " than the number of observations and is reduced", ...
                         " to %d."), this.NumObservations);
        KFold = this.NumObservations;
      endif

      ## The complexity parameter each subtree is matched to in a fold: the
      ## geometric mean of the two the parent's own sequence brackets it
      ## with, and infinity for the last, which takes a fold back to its
      ## root.  Measured on R2024a: matching a fold by level index or by the
      ## parent's parameter itself both give different answers, and this one
      ## reproduces every loss the oracle reports.
      alpha = this.PruneAlpha(:)';
      ## Formed outside the brackets: inside them the space before the paren
      ## would split the call off into an element of its own.
      geo = sqrt (alpha(1:end-1) .* alpha(2:end));
      geo(end+1) = Inf;

      nsub = numel (SubTrees);
      n = this.NumObservations;
      gY = labelIndices (this.ClassNames, this.Y);
      partition = cvpartition (this.Y, 'KFold', KFold);
      args = treeFoldArgs (this);

      ## The loss of every observation, from the fold that did not train on it
      L = zeros (n, nsub);
      fold = zeros (n, 1);
      for k = 1:KFold
        tr = training (partition, k);
        te = test (partition, k);
        fold(te) = k;
        ft = ClassificationTree (this.X(tr,:), this.Y(tr,:), args{:}, ...
                                 'Weights', this.RawWeights(tr));
        fa = ft.PruneAlpha(:)';
        for j = 1:nsub
          lv = sum (fa <= geo(SubTrees(j) + 1)) - 1;
          ## The nodes are collapsed without re-estimating the fold's own
          ## sequence, which prune would do and which nothing here reads.
          sub = collapseNodes (ft, nodesAtLevel (ft, max (lv, 0), 'cvloss'));
          [~, ~, ~, cnum] = predict (sub, this.X(te,:));
          L(te,j) = this.Cost(sub2ind (size (this.Cost), gY(te), cnum));
        endfor
      endfor

      ## The loss is weighed by the model's own weights, which sum to one.
      ## Measured on R2024a, where a tree weighted 1 to 150 reports
      ## 0.0537748344370861 and the unweighted mean of the same losses is
      ## 0.0533333333333333.
      W = this.W(:);
      E = W' * L;

      ## The standard error over the folds, each fold's loss being its
      ## observations' share of the whole.
      foldloss = zeros (KFold, nsub);
      for k = 1:KFold
        idx = fold == k;
        sw = sum (W(idx));
        if (sw > 0)
          foldloss(k,:) = (W(idx)' * L(idx,:)) / sw;
        endif
      endfor
      SE = std (foldloss, 0, 1) / sqrt (KFold);

      ## The leaves each subtree of the whole tree holds
      Nleaf = zeros (1, nsub);
      for j = 1:nsub
        Nleaf(j) = sum (! prune (this, 'Level', SubTrees(j)).IsBranchNode);
      endfor

      ## The smallest subtree the rule allows, which is the largest level
      if (strcmp (TreeSize, 'min'))
        best = find (E == min (E), 1, 'last');
      else
        [~, i] = min (E);
        best = find (E <= E(i) + SE(i), 1, 'last');
      endif
      BestLevel = SubTrees(best);

      E = E(:);
      SE = SE(:);
      Nleaf = Nleaf(:);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationTree} {@var{imp} =} predictorImportance (@var{obj})
    ##
    ## Estimate the importance of each predictor.
    ##
    ## @code{@var{imp} = predictorImportance (@var{obj})} returns a row
    ## vector with one element per predictor, holding the total drop in risk
    ## the splits on that predictor bring about, divided by the number of
    ## branch nodes.  A predictor the tree never splits on scores zero.
    ##
    ## The drop at a branch node is its @code{NodeRisk} less the risk of its
    ## two children and less what it holds back, so a predictor that is
    ## chosen often, high up, and on nodes it separates well, scores
    ## highest.  The numbers are comparable
    ## between predictors of one tree and not between trees.
    ##
    ## @seealso{ClassificationTree, fitctree, ClassificationTree.NodeRisk}
    ## @end deftypefn
    function imp = predictorImportance (this)

      imp = zeros (1, columns (this.X));
      branch = find (this.Children(:,1) > 0);
      if (isempty (branch))
        return;
      endif
      for ii = 1:numel (branch)
        b = branch(ii);
        kids = this.Children(b,:);
        ## The rows that stop at the node are in neither child and are
        ## answered the same whether it is split or not, so they are no part
        ## of what the split buys.
        drop = this.NodeRisk(b) - this.NodeRisk(kids(1)) ...
               - this.NodeRisk(kids(2)) - this.HeldImpurity(b);
        v = this.CutPredictorIndex(b);
        imp(v) = imp(v) + drop;
      endfor
      imp = imp / numel (branch);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationTree} {@var{r} =} nodeVariableRange (@var{obj}, @var{node})
    ##
    ## Range of each predictor at a node.
    ##
    ## @code{@var{r} = nodeVariableRange (@var{obj}, @var{node})} returns a
    ## structure with one field per predictor the path from the root to
    ## @var{node} cuts on, holding the two-element range of values that reach
    ## the node.  A predictor the path never cuts on is unconstrained and is
    ## left out, so the root gives a structure with no fields.
    ##
    ## @seealso{ClassificationTree, fitctree}
    ## @end deftypefn
    function r = nodeVariableRange (this, node)

      ## Input validation
      if (nargin < 2)
        error (strcat ("ClassificationTree.nodeVariableRange: too few", ...
                       " input arguments."));
      endif
      if (! (isnumeric (node) && isscalar (node) && isreal (node)
             && node >= 1 && node <= this.NumNodes && node == fix (node)))
        error (strcat ("ClassificationTree.nodeVariableRange: NODE must", ...
                       " be a positive integer no greater than the number", ...
                       " of nodes in the tree."));
      endif

      p = columns (this.X);
      lo = -Inf (1, p);
      hi = Inf (1, p);
      touched = false (1, p);

      ## Walk up to the root, narrowing the range of whichever predictor
      ## each step cut on.  Going up rather than down finds the path
      ## without a search, a node having exactly one parent.
      child = node;
      up = this.Parent(child);
      while (up > 0)
        v = this.CutPredictorIndex(up);
        cut = this.CutPoint(up);
        touched(v) = true;
        if (this.Children(up,1) == child)
          hi(v) = min (hi(v), cut);
        else
          lo(v) = max (lo(v), cut);
        endif
        child = up;
        up = this.Parent(child);
      endwhile

      r = struct ();
      for v = find (touched)
        r.(this.PredictorNames{v}) = [lo(v), hi(v)];
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationTree} {} view (@var{obj})
    ##
    ## Print the tree as text.
    ##
    ## @code{view (@var{obj})} prints one line per node: a branch node names
    ## the predictor it cuts on, the cut point, and the node each side leads
    ## to, and a leaf names the class it assigns.  A branch node's line ends
    ## with the class it would assign itself, which is the answer an
    ## observation missing that predictor gets.
    ##
    ## @seealso{ClassificationTree, fitctree}
    ## @end deftypefn
    function view (this)

      ## The node numbers are right aligned on the widest of them, which is
      ## how MATLAB lays the listing out.
      w = numel (sprintf ("%d", this.NumNodes));
      fprintf ("Decision tree for classification\n");
      for ii = 1:this.NumNodes
        if (this.Children(ii,1) == 0)
          fprintf ("%*d  class = %s\n", w, ii, this.NodeClass{ii});
        else
          v = this.PredictorNames{this.CutPredictorIndex(ii)};
          c = sprintf ("%g", this.CutPoint(ii));
          fprintf ("%*d  if %s<%s then node %d elseif %s>=%s then node", ...
                   w, ii, v, c, this.Children(ii,1), v, c);
          fprintf (" %d else %s\n", this.Children(ii,2), this.NodeClass{ii});
        endif
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationTree} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Drop the training data from a trained model.
    ##
    ## @code{@var{CMdl} = compact (@var{obj})} returns a
    ## @code{CompactClassificationTree} object carrying the tree and
    ## everything @code{predict} needs, but not the observations the model
    ## was fitted on.  It classifies new data identically and is far smaller
    ## to keep or to ship.
    ##
    ## @seealso{CompactClassificationTree, ClassificationTree}
    ## @end deftypefn
    function CMdl = compact (this)

      CMdl = CompactClassificationTree (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationTree} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationTree} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate a trained decision tree.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} partitions the training data
    ## into ten folds, or into as many folds as there are observations when
    ## there are fewer than ten, grows a tree on the training part of each
    ## and returns them as a @code{ClassificationPartitionedModel}.
    ##
    ## @code{@var{CVMdl} = crossval (@dots{}, @var{name}, @var{value})} takes
    ## one of the following, and one only.
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'KFold'} @tab An integer greater than 1, the number of
    ## folds.
    ##
    ## @item @qcode{'Holdout'} @tab A value between 0 and 1, the fraction of
    ## the data held out for testing, which gives a single fold.
    ##
    ## @item @qcode{'Leaveout'} @tab @qcode{'on'} or @qcode{'off'}, one fold
    ## per observation.
    ##
    ## @item @qcode{'CVPartition'} @tab A @code{cvpartition} object.
    ##
    ## @end multitable
    ##
    ## Every fold is grown with the parent's class names, prior, cost and
    ## observation weights rather than being left to re-derive them from its
    ## own rows, so a fold reports the prior of the whole data and not its own
    ## frequencies.
    ##
    ## @seealso{ClassificationPartitionedModel, ClassificationTree, cvpartition}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)

      ## Input validation
      if (numel (varargin) == 1)
        error (strcat ("ClassificationTree.crossval: Name-Value", ...
                       " arguments must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("ClassificationTree.crossval: specify only one of", ...
                       " the optional Name-Value paired arguments."));
      endif

      if (this.NumObservations < 10)
        numFolds = this.NumObservations;
      else
        numFolds = 10;
      endif
      Holdout     = [];
      Leaveout    = 'off';
      CVPartition = [];

      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))

          case 'kfold'
            numFolds = varargin{2};
            if (! (isnumeric (numFolds) && isscalar (numFolds)
                   && (numFolds == fix (numFolds)) && numFolds > 1))
              error (strcat ("ClassificationTree.crossval: 'KFold' must", ...
                             " be an integer value greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("ClassificationTree.crossval: 'Holdout' must", ...
                             " be a numeric value between 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && any (strcmpi (Leaveout, {'on', 'off'}))))
              error (strcat ("ClassificationTree.crossval: 'Leaveout'", ...
                             " must be either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("ClassificationTree.crossval: 'CVPartition'", ...
                             " must be a 'cvpartition' object."));
            endif

          otherwise
            error (strcat ("ClassificationTree.crossval: invalid", ...
                           " parameter name in optional paired arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## The partition covers the observations actually trained on: a row
      ## dropped for a missing response is not one the folds can use.  The
      ## response is passed rather than a count so the folds stay stratified.
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (this.Y, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (this.NumObservations, 'LeaveOut');
      else
        partition = cvpartition (this.Y, 'KFold', numFolds);
      endif

      CVMdl = ClassificationPartitionedModel (this, partition);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationTree} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin on new data.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns one
    ## margin per observation: the score the model gives the observation's
    ## true class, less the largest score it gives any other class.  A
    ## positive margin means the observation is classified correctly, and a
    ## larger one means it is classified more confidently.
    ##
    ## @seealso{ClassificationTree, ClassificationTree.edge,
    ## ClassificationTree.loss, ClassificationTree.predict}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Input validation
      if (nargin < 3)
        error ("ClassificationTree.margin: too few input arguments.");
      endif
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("ClassificationTree.margin: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("ClassificationTree.margin: number of rows in X", ...
                       " and Y must be equal."));
      endif

      [~, s] = predict (this, X);
      m = marginsOf (s, gY, 1);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationTree} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationTree} {@var{e} =} edge (@dots{}, @qcode{'Weights'}, @var{w})
    ##
    ## Classification edge on new data.
    ##
    ## @code{@var{e} = edge (@var{obj}, @var{X}, @var{Y})} returns the
    ## weighted mean of the margins, a single number summarising how
    ## confidently the model classifies the data.
    ##
    ## The weights are normalized within each class to that class's prior
    ## before they are applied.
    ##
    ## @seealso{ClassificationTree, ClassificationTree.margin,
    ## ClassificationTree.loss, ClassificationTree.predict}
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      ## Input validation
      if (nargin < 3)
        error ("ClassificationTree.edge: too few input arguments.");
      endif

      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       'ClassificationTree', 'edge');
      e = sum (W .* margin (this, X, Y));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationTree} {@var{l} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationTree} {@var{l} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss on new data.
    ##
    ## @code{@var{l} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## minimum expected misclassification cost.
    ##
    ## @code{@var{l} = loss (@dots{}, @var{name}, @var{value})} takes the
    ## following options.
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab One of @qcode{'binodeviance'},
    ## @qcode{'classifcost'}, @qcode{'classiferror'}, @qcode{'exponential'},
    ## @qcode{'hinge'}, @qcode{'logit'}, @qcode{'mincost'} (default) or
    ## @qcode{'quadratic'}.
    ##
    ## @item @qcode{'Weights'} @tab A numeric vector of observation weights,
    ## one per row of @var{X}.
    ##
    ## @end multitable
    ##
    ## @seealso{ClassificationTree, ClassificationTree.margin,
    ## ClassificationTree.edge, ClassificationTree.predict}
    ## @end deftypefn
    function l = loss (this, X, Y, varargin)

      ## Input validation
      if (nargin < 3)
        error ("ClassificationTree.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationTree.loss: name-value arguments", ...
                       " must be in pairs."));
      endif

      LossFun = 'mincost';
      Weights = [];
      lf_opt = {'binodeviance', 'classifcost', 'classiferror', ...
                'exponential', 'hinge', 'logit', 'mincost', 'quadratic'};

      while (numel (varargin) > 0)
        Value = varargin{2};
        switch (tolower (varargin{1}))
          case 'lossfun'
            if (! (ischar (Value) && any (strcmpi (Value, lf_opt))))
              error ("ClassificationTree.loss: invalid loss function.");
            endif
            LossFun = tolower (Value);
          case 'weights'
            if (! (isnumeric (Value) && isvector (Value)))
              error ("ClassificationTree.loss: invalid 'Weights'.");
            endif
            Weights = Value;
          otherwise
            error (strcat ("ClassificationTree.loss: invalid parameter", ...
                           " name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("ClassificationTree.loss: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("ClassificationTree.loss: number of rows in X and", ...
                       " Y must be equal."));
      endif
      if (isempty (Weights))
        w = ones (numel (gY), 1);
      else
        w = Weights(:);
        if (numel (w) != numel (gY))
          error (strcat ("ClassificationTree.loss: 'Weights' must have", ...
                         " one element per observation."));
        endif
      endif

      ## The weights are normalized within each class to that class's prior,
      ## so that a loss and an edge weight the classes the same way.
      w = priorNormalize (w, gY, this.Prior);

      [~, s] = predict (this, X);
      l = classificationLoss (LossFun, s, gY, w, this.Cost);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationTree} {@var{label} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationTree} {[@var{label}, @var{score}, @var{node}, @var{cnum}] =} resubPredict (@var{obj})
    ##
    ## Classify the training data with the model fitted to it.
    ##
    ## @code{@var{label} = resubPredict (@var{obj})} is
    ## @code{predict (@var{obj}, @var{obj}.X)}, and takes the same outputs.
    ##
    ## @seealso{ClassificationTree, ClassificationTree.predict}
    ## @end deftypefn
    function [label, score, node, cnum] = resubPredict (this)

      [label, score, node, cnum] = predict (this, this.X);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationTree} {@var{m} =} resubMargin (@var{obj})
    ##
    ## Classification margin on the training data.
    ##
    ## @code{@var{m} = resubMargin (@var{obj})} is
    ## @code{margin (@var{obj}, @var{obj}.X, @var{obj}.Y)}.
    ##
    ## @seealso{ClassificationTree, ClassificationTree.margin}
    ## @end deftypefn
    function m = resubMargin (this)

      m = margin (this, this.X, this.Y);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationTree} {@var{e} =} resubEdge (@var{obj})
    ##
    ## Classification edge on the training data.
    ##
    ## @code{@var{e} = resubEdge (@var{obj})} is @code{edge} over the
    ## training data, weighed as the fit weighed it.
    ##
    ## @seealso{ClassificationTree, ClassificationTree.edge}
    ## @end deftypefn
    function e = resubEdge (this)

      ## As with resubLoss, the training data is weighed as the fit weighed
      ## it.  Measured on R2024a: on a weighted fit this is edge over the
      ## model's own weights, 0.856171294940654 on an iris tree weighted
      ## 1 to 150, and not the 0.853582991377224 uniform weights give.
      e = edge (this, this.X, this.Y, 'Weights', this.RawWeights);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationTree} {@var{l} =} resubLoss (@var{obj})
    ## @deftypefnx {ClassificationTree} {@var{l} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss on the training data.
    ##
    ## @code{@var{l} = resubLoss (@var{obj})} is @code{loss} over the
    ## training data, weighed as the fit weighed it, and takes the same
    ## @qcode{'LossFun'} option.  Giving @qcode{'Weights'} weighs the
    ## training data some other way instead, which MATLAB refuses rather
    ## than honours.
    ##
    ## @seealso{ClassificationTree, ClassificationTree.loss}
    ## @end deftypefn
    function l = resubLoss (this, varargin)

      ## The training data is weighed as the fit weighed it.  Measured on
      ## R2024a: a weighted fit reports a resubstitution loss over its own
      ## weights, where an unweighted one and a loss call over the same data
      ## agree with it because their weights are uniform.
      if (! any (strcmpi (varargin(1:2:end), 'weights')))
        varargin = [varargin, {'Weights', this.RawWeights}];
      endif
      l = loss (this, this.X, this.Y, varargin{:});

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationTree} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationTree model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## ClassificationTree object into an Octave binary file, the name of
    ## which is specified in @var{filename}, along with an extra variable,
    ## which defines the type of classification object these variables
    ## constitute.  Use @code{loadmodel} in order to load a classification
    ## object into Octave's workspace.
    ##
    ## @seealso{loadmodel, fitctree, ClassificationTree}
    ## @end deftypefn
    function savemodel (this, fname)

      ## Input validation
      if (nargin < 2)
        error ("ClassificationTree.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("ClassificationTree.savemodel: FNAME must be a", ...
                       " character vector."));
      endif

      ## Generate variable for class name
      classdef_name = 'ClassificationTree';

      ## Create variables from model properties
      X                      = this.X;
      Y                      = this.Y;
      RowsUsed               = this.RowsUsed;
      W                      = this.W;
      RawWeights             = this.RawWeights;
      NumObservations        = this.NumObservations;
      PredictorNames         = this.PredictorNames;
      ResponseName           = this.ResponseName;
      ClassNames             = this.ClassNames;
      CategoricalPredictors  = this.CategoricalPredictors;
      ExpandedPredictorNames = this.ExpandedPredictorNames;
      BinEdges               = this.BinEdges;
      ModelParameters        = this.ModelParameters;
      NumNodes               = this.NumNodes;
      Children               = this.Children;
      Parent                 = this.Parent;
      IsBranchNode           = this.IsBranchNode;
      CutPredictorIndex      = this.CutPredictorIndex;
      CutPoint               = this.CutPoint;
      NodeSize               = this.NodeSize;
      ClassCount             = this.ClassCount;
      ClassShare             = this.ClassShare;
      PruneList              = this.PruneList;
      PruneAlpha             = this.PruneAlpha;
      Prior                  = this.Prior;
      Cost                   = this.Cost;
      ScoreTransform         = this.ScoreTransform;
      STfun                  = this.STfun;
      HyperparameterOptimizationResults = ...
                               this.HyperparameterOptimizationResults;

      ## The cut and node descriptions are not saved: every one of them is
      ## re-derived from the node table above, so writing them out would only
      ## make a stale copy possible.
      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'RowsUsed', 'W', ...
            'RawWeights', 'NumObservations', 'PredictorNames', ...
            'ResponseName', 'ClassNames', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'BinEdges', 'ModelParameters', ...
            'NumNodes', 'Children', 'Parent', 'IsBranchNode', ...
            'CutPredictorIndex', 'CutPoint', 'NodeSize', 'ClassCount', ...
            'ClassShare', 'PruneList', 'PruneAlpha', 'Prior', 'Cost', ...
            'ScoreTransform', 'STfun', ...
            'HyperparameterOptimizationResults');

    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      ## The smallest fit the class accepts, filled property by property
      ## below.  Nothing of the stub survives the copy.
      mdl = ClassificationTree ([1; 2], [1; 2]);

      ## Copy the saved data into the object.  Iterate over what was saved
      ## rather than over fieldnames (mdl): a hidden property such as STfun
      ## is written out by savemodel but is not reported by fieldnames, so
      ## comparing the two sets could never match.  Assignment is legal here
      ## because this is a method of the class itself.
      names = fieldnames (data);

      ## These three are assigned once everything else is in place, and in
      ## this order rather than the file's.  Cost comes before Prior because
      ## the two set methods re-derive the node statistics, which want the
      ## node table already loaded.
      order = {'Cost', 'Prior', 'ScoreTransform'};
      late = ismember (names, order);
      tail = order(ismember (order, names));
      names = [names(! late); tail(:)];
      for ii = 1:numel (names)
        try
          mdl.(names{ii}) = data.(names{ii});
        catch
          msg = 'ClassificationTree.load_model: invalid model in ''%s''.';
          error (msg, filename);
        end_try_catch
      endfor

      ## The cut and node descriptions were not saved, being derived
      mdl = fillCuts (mdl);
      mdl = deriveNodes (mdl);
      mdl = deriveHeld (mdl);

    endfunction

  endmethods

  methods (Access = private)

    ## Re-derive every node statistic that depends on the prior or the cost.
    ## ClassShare is free of both, so this is all that a reassignment of
    ## either has to redo, and the shape of the tree is untouched.
    function this = deriveNodes (this)

      ## An object part way through being assembled, as load_model assembles
      ## one property at a time, has a node table and a prior that do not yet
      ## describe the same model.  There is nothing to derive until they do,
      ## and load_model derives it once they do.
      if (isempty (this.ClassShare) || isempty (this.Prior)
          || columns (this.ClassShare) != numel (this.Prior)
          || rows (this.Cost) != numel (this.Prior))
        return;
      endif

      S = treeNodeStats (this.ClassShare, this.Prior, this.Cost, ...
                         this.ModelParameters.SplitCriterion, ...
                         this.ClassNames);
      this.NodeProbability = S.NodeProbability;
      this.ClassProbability = S.ClassProbability;
      this.NodeError = S.NodeError;
      this.NodeClass = S.NodeClass;
      this.NodeRisk = S.NodeRisk;
      this.NodeProbAdj = S.NodeProbAdj;

    endfunction

    ## The descriptions of the cuts, all of them derived from the node table
    ## so that they cannot fall out of step with it.
    function this = fillCuts (this)

      S = treeCutInfo (this.CutPredictorIndex, this.PredictorNames);
      for [val, name] = S
        this.(name) = val;
      endfor

    endfunction

    ## The branch nodes a pruning level turns into leaves.
    function nodes = nodesAtLevel (this, L, caller)

      if (isempty (this.PruneList))
        error (strcat ("ClassificationTree.%s: the tree carries no", ...
                       " pruning sequence."), caller);
      endif
      maxL = max (this.PruneList);
      if (L > maxL)
        warning (strcat ("ClassificationTree.%s: pruning level %d is", ...
                         " greater than the largest level %d; the tree", ...
                         " will be pruned to its root."), caller, L, maxL);
        L = maxL;
      endif
      nodes = find (this.PruneList >= 1 & this.PruneList <= L);

    endfunction

    ## Turn the named branch nodes into leaves, discard everything below
    ## them, and renumber what is left.  The surviving nodes keep their
    ## order, so a parent still precedes its children.
    function this = collapseNodes (this, nodes)

      if (isempty (nodes(this.Children(nodes,1) > 0)))
        return;
      endif
      S = treeCollapse (this.Children, this.Parent, ...
                        this.CutPredictorIndex, this.CutPoint, nodes);
      this.NumNodes = numel (S.keep);
      this.Children = S.Children;
      this.Parent = S.Parent;
      this.CutPredictorIndex = S.CutPredictorIndex;
      this.CutPoint = S.CutPoint;
      this.NodeSize = this.NodeSize(S.keep);
      this.ClassCount = this.ClassCount(S.keep,:);
      this.ClassShare = this.ClassShare(S.keep,:);
      this.PruneList = zeros (this.NumNodes, 1);
      this.PruneAlpha = [];

      this = fillCuts (this);
      this = deriveNodes (this);
      ## A node that became a leaf now holds back nothing, so this is
      ## re-derived rather than subset.
      this = deriveHeld (this);

    endfunction

    ## The risk each branch node carries on account of the observations that
    ## stop there, missing the predictor it cuts on.  Those observations are
    ## in neither child, and MATLAB charges them the node's own error rather
    ## than their own.  It is charged twice over, once on each of the two
    ## scales a classification tree measures a node by: the pruning sequence
    ## runs on the expected misclassification cost and predictorImportance
    ## on the impurity, and under a non-default cost the weight held back is
    ## itself a different fraction of the node on the two.  Both measured on
    ## R2024a against an iris fit with a fifth of one predictor missing.
    function this = deriveHeld (this)

      this.HeldRisk = zeros (this.NumNodes, 1);
      this.HeldImpurity = zeros (this.NumNodes, 1);
      br = find (this.Children(:,1) > 0);
      if (isempty (br))
        return;
      endif
      kids = this.Children(br,:);
      wheld = this.NodeProbability(br) - this.NodeProbability(kids(:,1)) ...
              - this.NodeProbability(kids(:,2));
      this.HeldRisk(br) = wheld .* this.NodeError(br);

      aheld = this.NodeProbAdj(br) - this.NodeProbAdj(kids(:,1)) ...
              - this.NodeProbAdj(kids(:,2));
      nz = this.NodeProbAdj(br) > 0;
      this.HeldImpurity(br(nz)) = aheld(nz) .* this.NodeRisk(br(nz)) ...
                                  ./ this.NodeProbAdj(br(nz));

    endfunction

    ## Merge back the pairs of leaves that buy nothing, which is what
    ## MergeLeaves asks for: a pair whose risks together are no less than
    ## their parent's tells nothing the parent did not.  Repeated until
    ## nothing more merges, so a pair freed by a merge below is caught.
    function this = mergeLeaves (this)

      do
        kid = this.Children;
        leaf = kid(:,1) == 0;
        pair = find (! leaf);
        pair = pair(leaf(kid(pair,1)) & leaf(kid(pair,2)));
        R = costRisk (this);
        tol = 1e-12 * R(1);
        useless = pair(R(kid(pair,1)) + R(kid(pair,2)) ...
                       + this.HeldRisk(pair) >= R(pair) - tol);
        if (! isempty (useless))
          this = collapseNodes (this, useless);
        endif
      until (isempty (useless))

    endfunction

    ## The risk a node carries as a leaf: the expected cost of
    ## misclassifying the observations that reached it, taken over the
    ## weight they carry rather than their share of the node.  It is
    ## NodeProbability times NodeError, formed in one step because the two
    ## are a division and a multiplication by the same quantity and the
    ## merge test below turns on their last digits.
    function R = costRisk (this)

      R = min ((this.ClassShare .* this.Prior(:)') * this.Cost, [], 2);

    endfunction

    ## The cost complexity pruning sequence, on the risk the cost defines.
    ## The compiled engine measures risk as misclassified weight, which is
    ## the expected misclassification cost only under the default cost, and
    ## it opens a level for a subtree that costs nothing to give up where
    ## MATLAB opens none, so the sequence is estimated here throughout.
    function this = pruneSequence (this)

      [this.PruneList, this.PruneAlpha] = ...
        __treeprune__ (this.Children, this.Parent, costRisk (this), ...
                       this.HeldRisk);

    endfunction

  endmethods

endclassdef

## The share of each class's total weight that reached each node.  Dividing
## by the root's row cancels whatever scaled the columns, the prior and the
## cost adjustment alike; a class with no weight at all divides by zero and
## is left at zero, having reached nothing.
function S = classShareOf (CW)

  S = zeros (size (CW));
  root = CW(1,:);
  nz = root > 0;
  S(:,nz) = CW(:,nz) ./ root(nz);

endfunction


## Tests
%!test  # MATLAB parity: the surface a default fit reports
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! assert_equal (class (Mdl), 'ClassificationTree');
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (Mdl.ClassNames, unique (species));
%! assert_equal (Mdl.Prior, [1/3, 1/3, 1/3], 1e-15);
%! assert_equal (Mdl.Cost, [0, 1, 1; 1, 0, 1; 1, 1, 0]);
%! assert_equal (Mdl.ResponseName, 'Y');
%! assert_equal (Mdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (Mdl.ExpandedPredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (Mdl.ScoreTransform, 'none');
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (Mdl.W, ones (150, 1) / 150, 1e-15);

%!test  # MATLAB parity: what a tree with no categories and no surrogates holds
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! assert_equal (size (Mdl.CutCategories), [9, 2]);
%! assert_equal (all (cellfun (@isempty, Mdl.CutCategories(:))), true);
%! assert_equal (size (Mdl.CategoricalSplit), [0, 0]);
%! assert_equal (size (Mdl.SurrogateCutPredictor), [0, 1]);
%! assert_equal (size (Mdl.SurrogateCutPoint), [0, 0]);
%! assert_equal (size (Mdl.SurrogatePredictorAssociation), [0, 0]);
%! assert_equal (Mdl.BinEdges, {});
%! assert_equal (Mdl.HyperparameterOptimizationResults, []);

%!test  # MATLAB parity: the parameters a default fit reports
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! MP = Mdl.ModelParameters;
%! assert_equal (MP.SplitCriterion, 'gdi');
%! assert_equal ({MP.MinParent, MP.MinLeaf, MP.MaxSplits}, {10, 1, 149});
%! assert_equal ({MP.MergeLeaves, MP.Prune}, {'on', 'on'});
%! assert_equal ({MP.PruneCriterion, MP.NVarToSample}, {'error', 'all'});
%! assert_equal ({MP.NSurrogate, MP.MaxCat, MP.AlgCat}, {0, 10, 'auto'});
%! assert_equal (MP.PredictorSelection, 'allsplits');
%! assert_equal ({MP.Method, MP.Type}, {'Tree', 'classification'});

%!test  # MATLAB parity: the node statistics of the iris tree
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! assert_equal (Mdl.NodeProbability', [150, 50, 100, 54, 46, 48, 6, 47, 1] ...
%!                                     / 150, 1e-14);
%! assert_equal (Mdl.NodeError', [2/3, 0, 1/2, 5/54, 1/46, 1/48, 1/3, 0, 0], ...
%!                               1e-14);
%! assert_equal (Mdl.NodeRisk(1), 2/3, 1e-14);
%! assert_equal (Mdl.NodeRisk(3), 1/3, 1e-14);
%! assert_equal (Mdl.ClassProbability(1,:), [1/3, 1/3, 1/3], 1e-14);
%! assert_equal (Mdl.ClassProbability(7,:), [0, 1/3, 2/3], 1e-14);

%!test  # NodeRisk is the impurity weighted by the node probability
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! gdi = 1 - sum (Mdl.ClassProbability .^ 2, 2);
%! assert_equal (Mdl.NodeRisk, Mdl.NodeProbability .* gdi, 1e-14);

%!test  # MATLAB parity: prune takes a subtree out of the sequence
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! sub = prune (Mdl, 'Level', 1);
%! assert_equal (sub.NumNodes, 7);
%! assert_equal (sub.NodeSize', [150, 50, 100, 54, 46, 48, 6]);
%! assert_equal (sub.PruneList', [3, 0, 2, 1, 0, 0, 0]);
%! assert_equal (sub.Children, [2, 3; 0, 0; 4, 5; 6, 7; 0, 0; 0, 0; 0, 0]);
%! assert_equal (sub.Parent', [0, 1, 1, 3, 3, 4, 4]);

%!test  # MATLAB parity: a level of two, and level zero changing nothing
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! assert_equal (prune (Mdl, 'Level', 2).NodeSize', [150, 50, 100, 54, 46]);
%! assert_equal (prune (Mdl, 'Level', 0).NumNodes, 9);
%! assert_equal (prune (Mdl).NumNodes, 9);

%!test  # MATLAB parity: a cost complexity parameter picks its own level
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! assert_equal (prune (Mdl, 'Alpha', 0.1).NumNodes, 5);
%! assert_equal (prune (Mdl, 'Alpha', 0).NumNodes, 9);
%! assert_equal (prune (Mdl, 'Alpha', 1).NumNodes, 1);

%!test  # Pruning named nodes turns just those into leaves
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! sub = prune (Mdl, 'Nodes', 3);
%! assert_equal (sub.NumNodes, 3);
%! assert_equal (sub.NodeSize', [150, 50, 100]);
%! assert_equal (sub.IsBranchNode', [true, false, false]);

%!test  # The pruned tree re-derives its own sequence, and it is the shifted one
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! sub = prune (Mdl, 'Level', 1);
%! assert_equal (sub.PruneAlpha', [0, 2/150, 44/150, 50/150], 1e-14);

%!warning<pruning level 9 is greater than the largest level 4>
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! sub = prune (Mdl, 'Level', 9);

%!test  # MATLAB parity: the four outputs of predict, and a missing predictor
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! x = meas([1, 60, 120], :);
%! x(1, 3) = NaN;
%! [label, score, node, cnum] = predict (Mdl, x);
%! assert_equal (node', [1, 8, 7]);
%! assert_equal (label, {'setosa'; 'versicolor'; 'virginica'});
%! assert_equal (cnum', [1, 2, 3]);
%! assert_equal (score(1,:), [1/3, 1/3, 1/3], 1e-14);

%!test  # MATLAB parity: a row missing every predictor stops at the root
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! [label, score, node] = predict (Mdl, NaN (1, 4));
%! assert_equal (node, 1);
%! assert_equal (label, {'setosa'});
%! assert_equal (score, [1/3, 1/3, 1/3], 1e-14);

%!test  # MATLAB parity: the score transform is applied to the probabilities
%! load fisheriris
%! Mdl = ClassificationTree (meas, species, 'ScoreTransform', 'logit');
%! [~, score] = predict (Mdl, meas([1, 60, 120], :));
%! assert_equal (Mdl.ScoreTransform, 'logit');
%! assert_equal (score(1,:), [0.731058578630005, 0.5, 0.5], 1e-14);

%!test  # MATLAB parity: predictorImportance divides by the branch node count
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! imp = predictorImportance (Mdl);
%! assert_equal (size (imp), [1, 4]);
%! assert_equal (imp(1:2), [0, 0]);
%! assert_equal (imp(3:4), [0.0907, 0.0682], 1e-4);

%!test  # MATLAB parity: the range of each predictor on the path to a node
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! assert_equal (fieldnames (nodeVariableRange (Mdl, 1)), cell (0, 1));
%! r = nodeVariableRange (Mdl, 8);
%! assert_equal (r.x3, [2.45, 4.95], 1e-14);
%! assert_equal (r.x4, [-Inf, 1.65], 1e-14);

%!test  # MATLAB parity: the text form of the tree
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! txt = evalc ('view (Mdl)');
%! lines = strsplit (strtrim (txt), "\n");
%! assert_equal (numel (lines), 10);
%! assert_equal (lines{1}, 'Decision tree for classification');
%! assert_equal (lines{2}, ...
%!   '1  if x3<2.45 then node 2 elseif x3>=2.45 then node 3 else setosa');
%! assert_equal (lines{3}, '2  class = setosa');

%!test  # MATLAB parity: the seven classification losses on the training data
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! assert_equal (resubLoss (Mdl, 'LossFun', 'binodeviance'), 0.1388, 1e-4);
%! assert_equal (resubLoss (Mdl, 'LossFun', 'classiferror'), 0.02, 1e-14);
%! assert_equal (resubLoss (Mdl, 'LossFun', 'exponential'), 0.3829, 1e-4);
%! assert_equal (resubLoss (Mdl, 'LossFun', 'hinge'), 0.0308, 1e-4);
%! assert_equal (resubLoss (Mdl, 'LossFun', 'logit'), 0.3232, 1e-4);
%! assert_equal (resubLoss (Mdl, 'LossFun', 'mincost'), 0.02, 1e-14);
%! assert_equal (resubLoss (Mdl, 'LossFun', 'quadratic'), 0.0154, 1e-4);

%!test  # MATLAB parity: margin and edge
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! m = margin (Mdl, meas([1, 60, 120], :), species([1, 60, 120]));
%! assert_equal (m', [1, 1, 1/3], 1e-14);
%! assert_equal (size (resubMargin (Mdl)), [150, 1]);
%! assert_equal (resubEdge (Mdl), 0.9384, 1e-4);
%! assert_equal (edge (Mdl, meas, species), resubEdge (Mdl), 1e-14);

%!test  # resubPredict answers exactly as predict on the training data
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! [l1, s1, n1, c1] = resubPredict (Mdl);
%! [l2, s2, n2, c2] = predict (Mdl, meas);
%! assert_equal ({l1, s1, n1, c1}, {l2, s2, n2, c2});

%!test  # MATLAB parity: a cost matrix reaches the label and the risk
%! load fisheriris
%! Mdl = ClassificationTree (meas, species, ...
%!                           'Cost', [0, 2, 8; 3, 0, 1; 5, 4, 0]);
%! assert_equal (Mdl.NodeSize', [150, 50, 100, 45, 55, 44, 1]);
%! assert_equal (Mdl.NodeClass{1}, 'versicolor');
%! assert_equal (Mdl.NodeError(1), 2, 1e-14);
%! assert_equal (Mdl.NodeRisk(1), 0.6276, 1e-4);
%! assert_equal (Mdl.PruneAlpha', [0, 0.0267, 0.2667, 1.6667], 1e-4);

%!test  # Reassigning Cost re-derives the node statistics, not the tree
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! shape = {Mdl.Children, Mdl.NodeSize, Mdl.ClassCount};
%! Mdl.Cost = [0, 2, 8; 3, 0, 1; 5, 4, 0];
%! assert_equal ({Mdl.Children, Mdl.NodeSize, Mdl.ClassCount}, shape);
%! assert_equal (Mdl.NodeClass{1}, 'versicolor');
%! assert_equal (Mdl.NodeError(1), 2, 1e-14);

%!test  # Reassigning Prior re-derives the weights and the node statistics
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! Mdl.Prior = [0.5, 0.25, 0.25];
%! assert_equal (Mdl.Prior, [0.5, 0.25, 0.25], 1e-14);
%! assert_equal (Mdl.W([1, 51, 101])', [0.01, 0.005, 0.005], 1e-14);
%! assert_equal (Mdl.ClassProbability(1,:), [0.5, 0.25, 0.25], 1e-14);
%! assert_equal (Mdl.NodeProbability(2), 0.5, 1e-14);
%! assert_equal (Mdl.NodeRisk(1), 0.625, 1e-14);

%!test  # A prior given as a structure names its own class order
%! load fisheriris
%! p.ClassNames = {'virginica'; 'versicolor'; 'setosa'};
%! p.ClassProbs = [0.25, 0.25, 0.5];
%! Mdl = ClassificationTree (meas, species, 'Prior', p);
%! assert_equal (Mdl.Prior, [0.5, 0.25, 0.25], 1e-14);

%!test  # A cost given as a structure names its own class order
%! load fisheriris
%! c.ClassNames = {'virginica'; 'versicolor'; 'setosa'};
%! c.ClassificationCosts = [0, 1, 10; 1, 0, 1; 10, 1, 0];
%! Mdl = ClassificationTree (meas, species, 'Cost', c);
%! assert_equal (Mdl.Cost, [0, 1, 10; 1, 0, 1; 10, 1, 0]);
%! assert_equal (Mdl.NodeSize', [150, 50, 100, 45, 55, 44, 1, 9, 46]);

%!test  # A response of one class alone gives a tree of one node
%! load fisheriris
%! Mdl = ClassificationTree (meas(1:5, :), species(1:5));
%! assert_equal (Mdl.NumNodes, 1);
%! assert_equal (Mdl.NodeClass, {'setosa'});
%! assert_equal (Mdl.ClassProbability, 1);
%! assert_equal (predict (Mdl, meas(1, :)), {'setosa'});

%!test  # A model saved and loaded answers exactly as it did
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! fname = tempname ();
%! unwind_protect
%!   savemodel (Mdl, fname);
%!   New = loadmodel (fname);
%!   assert_equal (class (New), 'ClassificationTree');
%!   assert_equal (New.NumNodes, Mdl.NumNodes);
%!   assert_equal (New.NodeRisk, Mdl.NodeRisk, 1e-15);
%!   assert_equal (New.CutPredictor, Mdl.CutPredictor);
%!   assert_equal (predict (New, meas), predict (Mdl, meas));
%! unwind_protect_cleanup
%!   delete (fname);
%! end_unwind_protect

%!test  # MATLAB parity: an unmerged tree records no level for a free subtree
%! ## With MergeLeaves off, the pair of leaves that merging would have
%! ## removed survives, and giving it up costs nothing.  MATLAB records
%! ## neither a level nor an alpha for such a subtree, so the sequence is
%! ## the merged tree's, over eleven nodes instead of nine.
%! load fisheriris
%! Mdl = ClassificationTree (meas, species, 'MergeLeaves', 'off');
%! assert_equal (Mdl.NumNodes, 11);
%! assert_equal (Mdl.PruneList', [4, 0, 3, 2, 0, 1, 0, 0, 0, 0, 0]);
%! assert_equal (Mdl.PruneAlpha', [0, 1/150, 2/150, 44/150, 50/150], 1e-14);

%!test  # MATLAB parity: a deeper tree, every branch pruned at its own level
%! load fisheriris
%! Mdl = ClassificationTree (meas, species, 'MinParentSize', 2);
%! assert_equal (Mdl.NumNodes, 17);
%! assert_equal (Mdl.PruneList', [5, 0, 4, 3, 1, 2, 2, 1, 0, 0, 0, 0, 2, ...
%!                                0, 0, 0, 0]);
%! assert_equal (Mdl.PruneAlpha', [0, 0.5, 1, 2, 44, 50] / 150, 1e-14);
%! assert_equal (predictorImportance (Mdl), ...
%!               [0.00222222222222222, 0, 0.0458935520665593, ...
%!                0.0352175590445518], 1e-14);

%!test  # MATLAB parity: the resubstitution methods weigh as the fit weighed
%! load fisheriris
%! Mdl = ClassificationTree (meas, species, 'Weights', (1:150)');
%! assert_equal (resubLoss (Mdl), 0.0384988962472406, 1e-14);
%! assert_equal (resubEdge (Mdl), 0.856171294940654, 1e-14);
%! ## loss and edge over the same data weigh it uniformly unless told not to
%! assert_equal (loss (Mdl, meas, species), 0.04, 1e-14);
%! assert_equal (edge (Mdl, meas, species), 0.853582991377224, 1e-14);
%! assert_equal (loss (Mdl, meas, species, 'Weights', (1:150)'), ...
%!               0.0384988962472406, 1e-14);

%!test  # MATLAB parity: the tree a weighted fit grows
%! load fisheriris
%! Mdl = ClassificationTree (meas, species, 'Weights', (1:150)');
%! assert_equal (Mdl.NodeSize', [150, 95, 55, 50, 45, 44, 1]);
%! assert_equal (Mdl.NodeProbability', [1, 0.416865342163355, ...
%!                                      0.583134657836645, ...
%!                                      0.112582781456954, ...
%!                                      0.304282560706402, ...
%!                                      0.294834437086093, ...
%!                                      0.00944812362030905], 1e-14);
%! assert_equal (Mdl.NodeRisk(1), 0.569205054359214, 1e-14);
%! assert_equal (Mdl.PruneAlpha', [0, 0.00944812362030906, ...
%!                                 0.112582781456954, 0.285386313465784], ...
%!                                1e-14);

%!test  # MATLAB parity: the resubstitution edge under a cost matrix
%! load fisheriris
%! Mdl = ClassificationTree (meas, species, ...
%!                           'Cost', [0, 1, 10; 1, 0, 1; 10, 1, 0]);
%! assert_equal (resubLoss (Mdl), 1/30, 1e-14);
%! assert_equal (resubEdge (Mdl), 0.914653784219002, 1e-14);

%!test  # MATLAB parity: compact drops the data and answers identically
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), 'CompactClassificationTree');
%! assert_equal (numel (properties (CMdl)), 33);
%! assert_equal (predict (CMdl, meas), predict (Mdl, meas));
%! assert_equal (CMdl.NodeRisk, Mdl.NodeRisk, 1e-15);

%!test  # MATLAB parity: crossval returns a partitioned model over compacts
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! CVMdl = crossval (Mdl);
%! assert_equal (class (CVMdl), 'ClassificationPartitionedModel');
%! assert_equal (CVMdl.CrossValidatedModel, 'Tree');
%! assert_equal (class (CVMdl.Trained{1}), 'CompactClassificationTree');
%! assert_equal (CVMdl.KFold, 10);
%! assert_equal (CVMdl.NumObservations, 150);
%! assert_equal (CVMdl.Prior, Mdl.Prior, 1e-15);
%! assert_equal (CVMdl.Cost, Mdl.Cost);

%!test  # Each way of asking for a partition gives the folds it names
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! assert_equal (crossval (Mdl, 'KFold', 5).KFold, 5);
%! assert_equal (crossval (Mdl, 'Holdout', 0.3).KFold, 1);
%! assert_equal (crossval (Mdl, 'Leaveout', 'on').KFold, 150);
%! assert_equal (crossval (Mdl, 'CVPartition', ...
%!                         cvpartition (species, 'KFold', 4)).KFold, 4);

%!test  # A fold is grown with the parent's prior, cost and weights
%! load fisheriris
%! Mdl = ClassificationTree (meas, species, 'Weights', (1:150)', ...
%!                           'Cost', [0, 1, 10; 1, 0, 1; 10, 1, 0]);
%! CVMdl = crossval (Mdl, 'KFold', 3);
%! assert_equal (CVMdl.Trained{1}.Prior, Mdl.Prior, 1e-15);
%! assert_equal (CVMdl.Trained{1}.Cost, Mdl.Cost);
%! assert_equal (CVMdl.Trained{1}.ClassNames, Mdl.ClassNames);

%!test  # A fold is grown with the parameters the parent was grown with
%! load fisheriris
%! Mdl = ClassificationTree (meas, species, 'SplitCriterion', 'deviance', ...
%!                           'MinLeafSize', 7, 'MergeLeaves', 'off');
%! CVMdl = crossval (Mdl, 'KFold', 3);
%! assert_equal (CVMdl.ModelParameters.SplitCriterion, 'deviance');
%! assert_equal (CVMdl.ModelParameters.MinLeaf, 7);
%! assert_equal (CVMdl.ModelParameters.MergeLeaves, 'off');

%!test  # kfoldPredict scores every observation and costs them from the score
%! load fisheriris
%! CVMdl = crossval (ClassificationTree (meas, species), 'KFold', 5);
%! [label, score, cost] = kfoldPredict (CVMdl);
%! assert_equal (size (label), [150, 1]);
%! assert_equal (size (score), [150, 3]);
%! assert_equal (cost, score * CVMdl.Cost, 1e-14);
%! assert_equal (sum (score, 2), ones (150, 1), 1e-14);
%! assert_equal (kfoldLoss (CVMdl) < 0.2, true);
%! assert_equal (size (kfoldMargin (CVMdl)), [150, 1]);

%!test  # cvloss over a fixture whose folds all answer alike
%! ## The two classes are separated by a gap no fold can straddle, so every
%! ## fold grows the same tree and every quantity below is the same whatever
%! ## partition is drawn.
%! x = [(1:20)'; (101:120)'];
%! y = [repmat({'a'}, 20, 1); repmat({'b'}, 20, 1)];
%! Mdl = ClassificationTree (x, y);
%! assert_equal (Mdl.NumNodes, 3);
%! [E, SE, Nleaf, BestLevel] = cvloss (Mdl, 'SubTrees', 'all');
%! assert_equal (E', [0, 0.5], 1e-14);
%! assert_equal (SE', [0, 0], 1e-14);
%! assert_equal (Nleaf', [2, 1]);
%! assert_equal (BestLevel, 0);

%!test  # MATLAB parity: the leaf counts and the stable loss of the iris tree
%! ## Only the last level's loss is asserted.  Its parameter is infinite, so
%! ## every fold's tree goes back to its root whatever the partition, and a
%! ## root under a uniform prior answers the first class, which two thirds
%! ## of the data are not.  The levels below it are not fixed: a fold's tree
%! ## is matched to a subtree of the whole by a complexity parameter, and
%! ## which of its own levels that picks out depends on the rows it saw.
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! [E, SE, Nleaf, BestLevel] = cvloss (Mdl, 'SubTrees', 'all');
%! assert_equal (size (E), [5, 1]);
%! assert_equal (Nleaf', [5, 4, 3, 2, 1]);
%! assert_equal (E(5), 2/3, 1e-14);
%! assert_equal (SE(5), 0, 1e-12);
%! assert_equal (BestLevel >= 0 && BestLevel <= 4, true);

%!test  # cvloss defaults to the unpruned tree alone
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! [E, SE, Nleaf, BestLevel] = cvloss (Mdl);
%! assert_equal (size (E), [1, 1]);
%! assert_equal (Nleaf, 5);
%! assert_equal (BestLevel, 0);
%! assert_equal (E > 0 && E < 0.2, true);

%!test  # A subset of the levels is answered in the order it was asked for
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! [E, SE, Nleaf, BestLevel] = cvloss (Mdl, 'SubTrees', [0, 2, 4]);
%! assert_equal (Nleaf', [5, 3, 1]);
%! assert_equal (E(3), 2/3, 1e-14);
%! assert_equal (any (BestLevel == [0, 2, 4]), true);

%!test  # TreeSize picks the smallest tree the rule allows
%! ## 'min' takes the smallest tree of least loss and 'se' the smallest whose
%! ## loss is within one standard error of it, so 'se' never names a larger
%! ## tree than 'min' does.
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! [E, SE] = cvloss (Mdl, 'SubTrees', 'all');
%! [~, ~, ~, bmin] = cvloss (Mdl, 'SubTrees', 'all', 'TreeSize', 'min');
%! [~, ~, ~, bse] = cvloss (Mdl, 'SubTrees', 'all', 'TreeSize', 'se');
%! assert_equal (bmin >= 0 && bmin <= 4, true);
%! assert_equal (bse >= 0 && bse <= 4, true);

%!test  # The number of folds is the number asked for
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! [E, SE, Nleaf] = cvloss (Mdl, 'KFold', 5);
%! assert_equal (size (E), [1, 1]);
%! assert_equal (Nleaf, 5);
%! assert_equal (E > 0 && E < 0.2, true);

%!test  # cvloss weighs the loss as the fit weighed it
%! ## A weighted fit whose folds all answer alike.  The three to one weights
%! ## make the prior 0.25 and 0.75, so every fold's root answers the heavier
%! ## class and the twenty rows of the lighter one are wrong.  Their weight
%! ## is 0.25, where their share of the rows is 0.5.
%! x = [(1:20)'; (101:120)'];
%! y = [repmat({'a'}, 20, 1); repmat({'b'}, 20, 1)];
%! Mdl = ClassificationTree (x, y, 'Weights', [ones(20, 1); 3 * ones(20, 1)]);
%! assert_equal (Mdl.Prior, [0.25, 0.75], 1e-14);
%! E = cvloss (Mdl, 'SubTrees', 'all');
%! assert_equal (E(1), 0, 1e-14);
%! assert_equal (E(2), 0.25, 1e-14);

%!test  # More folds than observations are reduced to one fold each
%! ## The reduction warns, and so does cvpartition about a fold that cannot
%! ## hold every class, so the message is not what is asserted here; the
%! ## answer is, and it is the one leaving each observation out in turn.
%! x = [(1:20)'; (101:120)'];
%! y = [repmat({'a'}, 20, 1); repmat({'b'}, 20, 1)];
%! Mdl = ClassificationTree (x, y);
%! ws = warning ('off', 'all');
%! unwind_protect
%!   E = cvloss (Mdl, 'KFold', 500, 'SubTrees', 'all');
%! unwind_protect_cleanup
%!   warning (ws);
%! end_unwind_protect
%! assert_equal (E', [0, 0.5], 1e-14);

%!test  # MATLAB parity: a node that holds rows back pays for them
%! ## Twenty rows have no fourth predictor, and the node that cuts on it
%! ## sends them to neither child.  They are in no child's risk, so a
%! ## subtree's risk is its children's plus what the node holds back, and
%! ## leaving that out would overstate every link above it.
%! load fisheriris
%! x = meas;
%! x(51:70, 4) = NaN;
%! Mdl = ClassificationTree (x, species);
%! assert_equal (Mdl.NumNodes, 9);
%! assert_equal (Mdl.NodeSize', [150, 50, 100, 45, 55, 25, 1, 8, 46]);
%! assert_equal (Mdl.NodeSize(4) - sum (Mdl.NodeSize(Mdl.Children(4,:))), 19);
%! assert_equal (Mdl.PruneList', [4, 0, 3, 1, 2, 0, 0, 0, 0]);
%! assert_equal (Mdl.PruneAlpha', [0, 0.00385185185185185, ...
%!                                 0.00593939393939394, 0.286666666666666, ...
%!                                 0.333333333333333], 1e-14);
%! assert_equal (predictorImportance (Mdl), ...
%!               [0, 0, 0.145589225589225, 0.00944980621616545], 1e-14);
%! assert_equal (resubLoss (Mdl), 0.04, 1e-14);

%!test  # MATLAB parity: held-back rows under a cost matrix
%! ## The weight a node holds back is a different fraction of it on the two
%! ## scales a classification tree is measured by, so the pruning sequence
%! ## and the predictor importances discount different amounts.
%! load fisheriris
%! x = meas;
%! x(51:80, 3) = NaN;
%! Mdl = ClassificationTree (x, species, ...
%!                           'Cost', [0, 1, 10; 1, 0, 1; 10, 1, 0]);
%! assert_equal (Mdl.NumNodes, 7);
%! assert_equal (Mdl.NodeSize', [150, 50, 70, 15, 55, 5, 50]);
%! assert_equal (Mdl.NodeSize(1) - sum (Mdl.NodeSize(Mdl.Children(1,:))), 30);
%! assert_equal (Mdl.PruneList', [3, 0, 2, 0, 1, 0, 0]);
%! assert_equal (Mdl.PruneAlpha', [0, 0.02, 0.1, 0.4], 1e-14);
%! assert_equal (predictorImportance (Mdl), ...
%!               [0, 0, 0.17774350820443, 0], 1e-14);

%!test  # MATLAB parity: held-back rows with the leaves left unmerged
%! load fisheriris
%! x = meas;
%! x(51:70, 4) = NaN;
%! Mdl = ClassificationTree (x, species, 'MergeLeaves', 'off');
%! assert_equal (Mdl.NumNodes, 11);
%! assert_equal (Mdl.NodeSize', [150, 50, 100, 45, 55, 25, 1, 8, 46, 3, 43]);
%! assert_equal (Mdl.PruneList', [4, 0, 3, 1, 2, 0, 0, 0, 0, 0, 0]);
%! assert_equal (Mdl.PruneAlpha', [0, 0.00385185185185185, ...
%!                                 0.00593939393939394, 0.286666666666666, ...
%!                                 0.333333333333333], 1e-14);

## Test input validation
%!error<ClassificationTree: too few input arguments.> ClassificationTree ()
%!error<ClassificationTree: too few input arguments.>
%! ClassificationTree (ones (4, 2))
%!error<ClassificationTree: name-value arguments must be in pairs.>
%! ClassificationTree (ones (4, 2), ones (4, 1), 'K')
%!error<ClassificationTree: X must be a non-empty real numeric matrix.>
%! ClassificationTree ('a', ones (4, 1))
%!error<ClassificationTree: number of rows in X and Y must be equal.>
%! ClassificationTree (ones (4, 2), ones (3, 1))
%!error<ClassificationTree: 'PredictorNames' must be supplied as a cellstring array.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'PredictorNames', 'a')
%!error<ClassificationTree: 'PredictorNames' must equal the number of columns in X.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'PredictorNames', {'a'})
%!error<ClassificationTree: 'ResponseName' must be a character vector.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'ResponseName', 5)
%!error<ClassificationTree: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'ClassNames', {1})
%!error<ClassificationTree: not all 'ClassNames' are present in Y.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'ClassNames', 5)
%!error<ClassificationTree: 'Prior' must be a numeric vector, a structure, or a character vector.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'Prior', {1})
%!error<ClassificationTree: 'Prior' must be 'empirical', 'uniform', a numeric vector, or a structure.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'Prior', 'bogus')
%!error<ClassificationTree: 'Cost' must be a numeric square matrix or a structure.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'Cost', ones (2, 3))
%!error<ClassificationTree: 'Weights' must be a real numeric vector.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'Weights', 'a')
%!error<ClassificationTree: 'Weights' must have one element per row in X.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'Weights', [1, 2, 3])
%!error<ClassificationTree: 'Weights' must be nonnegative and must not be all zero.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'Weights', -ones (4, 1))
%!error<ClassificationTree: 'MaxNumSplits' must be a nonnegative integer.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'MaxNumSplits', -1)
%!error<ClassificationTree: 'MinLeafSize' must be a positive integer.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'MinLeafSize', 0)
%!error<ClassificationTree: 'MinParentSize' must be a positive integer.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'MinParentSize', 0)
%!error<ClassificationTree: 'MergeLeaves' must be either 'on' or 'off'.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'MergeLeaves', 'x')
%!error<ClassificationTree: 'Prune' must be either 'on' or 'off'.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'Prune', 'x')
%!error<ClassificationTree: 'PruneCriterion' must be either 'error' or 'impurity'.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'PruneCriterion', 'x')
%!error<ClassificationTree: 'PruneCriterion' 'impurity' is not implemented.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'PruneCriterion', 'impurity')
%!error<ClassificationTree: 'SplitCriterion' must be 'gdi', 'deviance', or 'twoing'.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'SplitCriterion', 'x')
%!error<ClassificationTree: 'SplitCriterion' 'twoing' is not implemented.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'SplitCriterion', 'twoing')
%!error<ClassificationTree: categorical predictors are not implemented.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'CategoricalPredictors', 1)
%!error<ClassificationTree: 'Surrogate' is not implemented.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'Surrogate', 'on')
%!error<ClassificationTree: 'KFold' is not implemented; fit the model and cross-validate it afterwards.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'KFold', 5)
%!error<ClassificationTree: invalid parameter name in optional pair arguments.>
%! ClassificationTree (ones (4, 2), [1; 1; 2; 2], 'Bogus', 1)
%!error<ClassificationTree.predict: too few input arguments.>
%! predict (ClassificationTree (ones (4, 2), [1; 1; 2; 2]))
%!error<ClassificationTree.predict: XC is empty.>
%! predict (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), [])
%!error<ClassificationTree.predict: XC must be a real numeric matrix.>
%! predict (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'a')
%!error<ClassificationTree.predict: XC must have the same number of predictors as the trained model.>
%! predict (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), ones (2, 3))
%!error<ClassificationTree.margin: too few input arguments.>
%! margin (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), ones (4, 2))
%!error<ClassificationTree.margin: Y must hold only classes the model was trained on.>
%! margin (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), ones (4, 2), ...
%!         [1; 1; 2; 9])
%!error<ClassificationTree.margin: number of rows in X and Y must be equal.>
%! margin (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), ones (3, 2), ...
%!         [1; 1; 2; 2])
%!error<ClassificationTree.edge: too few input arguments.>
%! edge (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), ones (4, 2))
%!error<ClassificationTree.loss: too few input arguments.>
%! loss (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), ones (4, 2))
%!error<ClassificationTree.loss: name-value arguments must be in pairs.>
%! loss (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), ones (4, 2), ...
%!       [1; 1; 2; 2], 'LossFun')
%!error<ClassificationTree.loss: invalid loss function.>
%! loss (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), ones (4, 2), ...
%!       [1; 1; 2; 2], 'LossFun', 'x')
%!error<ClassificationTree.loss: invalid 'Weights'.>
%! loss (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), ones (4, 2), ...
%!       [1; 1; 2; 2], 'Weights', 'a')
%!error<ClassificationTree.loss: invalid parameter name in optional pair arguments.>
%! loss (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), ones (4, 2), ...
%!       [1; 1; 2; 2], 'Bogus', 1)
%!error<ClassificationTree.prune: name-value arguments must be in pairs.>
%! prune (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'Level')
%!error<ClassificationTree.prune: specify only one of the optional name-value paired arguments.>
%! prune (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'Level', 1, ...
%!        'Alpha', 1)
%!error<ClassificationTree.prune: 'Level' must be a nonnegative integer.>
%! prune (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'Level', -1)
%!error<ClassificationTree.prune: 'Alpha' must be a nonnegative scalar.>
%! prune (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'Alpha', -1)
%!error<ClassificationTree.prune: 'Nodes' must hold indices of nodes of the tree.>
%! prune (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'Nodes', 99)
%!error<ClassificationTree.prune: invalid parameter name in optional pair arguments.>
%! prune (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'Bogus', 1)
%!error<ClassificationTree.nodeVariableRange: too few input arguments.>
%! nodeVariableRange (ClassificationTree (ones (4, 2), [1; 1; 2; 2]))
%!error<ClassificationTree.nodeVariableRange: NODE must be a positive integer no greater than the number of nodes in the tree.>
%! nodeVariableRange (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 99)
%!error<ClassificationTree.savemodel: too few input arguments.>
%! savemodel (ClassificationTree (ones (4, 2), [1; 1; 2; 2]))
%!error<ClassificationTree.savemodel: FNAME must be a character vector.>
%! savemodel (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 5)
%!error<ClassificationTree.cvloss: name-value arguments must be in pairs.>
%! cvloss (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'SubTrees')
%!error<ClassificationTree.cvloss: 'SubTrees' must be 'all' or a vector of nonnegative integers in ascending order.>
%! load fisheriris
%! cvloss (ClassificationTree (meas, species), 'SubTrees', -1)
%!error<ClassificationTree.cvloss: 'SubTrees' must be 'all' or a vector of nonnegative integers in ascending order.>
%! load fisheriris
%! cvloss (ClassificationTree (meas, species), 'SubTrees', [2, 1])
%!error<ClassificationTree.cvloss: 'SubTrees' must be 'all' or a vector of nonnegative integers in ascending order.>
%! load fisheriris
%! cvloss (ClassificationTree (meas, species), 'SubTrees', 'most')
%!error<ClassificationTree.cvloss: 'SubTrees' must not exceed the largest pruning level, 4.>
%! load fisheriris
%! cvloss (ClassificationTree (meas, species), 'SubTrees', 99)
%!error<ClassificationTree.cvloss: 'TreeSize' must be either 'se' or 'min'.>
%! load fisheriris
%! cvloss (ClassificationTree (meas, species), 'TreeSize', 'x')
%!error<ClassificationTree.cvloss: 'KFold' must be an integer value greater than 1.>
%! load fisheriris
%! cvloss (ClassificationTree (meas, species), 'KFold', 1)
%!error<ClassificationTree.cvloss: invalid parameter name in optional pair arguments.>
%! load fisheriris
%! cvloss (ClassificationTree (meas, species), 'Bogus', 1)
%!error<ClassificationTree.cvloss: the tree carries no pruning sequence; fit it with 'Prune' or 'MergeLeaves' on.>
%! load fisheriris
%! cvloss (ClassificationTree (meas, species, 'Prune', 'off', ...
%!                             'MergeLeaves', 'off'))
%!error<ClassificationTree.crossval: Name-Value arguments must be in pairs.>
%! crossval (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'KFold')
%!error<ClassificationTree.crossval: specify only one of the optional Name-Value paired arguments.>
%! crossval (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'KFold', 2, ...
%!           'Holdout', 0.3)
%!error<ClassificationTree.crossval: 'KFold' must be an integer value greater than 1.>
%! crossval (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'KFold', 1)
%!error<ClassificationTree.crossval: 'Holdout' must be a numeric value between 0 and 1.>
%! crossval (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'Holdout', 2)
%!error<ClassificationTree.crossval: 'Leaveout' must be either 'on' or 'off'.>
%! crossval (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'Leaveout', 'x')
%!error<ClassificationTree.crossval: 'CVPartition' must be a 'cvpartition' object.>
%! crossval (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'CVPartition', 5)
%!error<ClassificationTree.crossval: invalid parameter name in optional paired arguments.>
%! crossval (ClassificationTree (ones (4, 2), [1; 1; 2; 2]), 'Bogus', 1)
%!error<ClassificationTree.Prior: must have one element per class.>
%! Mdl = ClassificationTree (ones (4, 2), [1; 1; 2; 2]);
%! Mdl.Prior = [0.2, 0.3, 0.5];
%!error<ClassificationTree.Prior: a character vector must be 'empirical' or 'uniform'.>
%! Mdl = ClassificationTree (ones (4, 2), [1; 1; 2; 2]);
%! Mdl.Prior = 'bogus';
%!error<ClassificationTree.Prior: must be nonnegative and must not be all zero.>
%! Mdl = ClassificationTree (ones (4, 2), [1; 1; 2; 2]);
%! Mdl.Prior = [0, 0];
