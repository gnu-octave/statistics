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

classdef CompactClassificationTree
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactClassificationTree
  ##
  ## Compact binary decision tree for classification
  ##
  ## A @code{CompactClassificationTree} object carries the tree a
  ## @code{ClassificationTree} model grew and everything @code{predict} needs,
  ## but not the observations it was fitted on.  It classifies new data
  ## identically to the model it came from, and is far smaller to keep or to
  ## ship.
  ##
  ## Create one with the @code{compact} method of a @code{ClassificationTree}
  ## object.  Because it holds no training data, it has no @code{resub}
  ## methods and cannot be cross-validated, and it cannot be pruned: the
  ## pruning sequence is reported but taking a subtree out of it rewrites the
  ## node table, which is work for the model that still has its data.
  ##
  ## @seealso{ClassificationTree, fitctree}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} NumNodes
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
    ## @deftp {CompactClassificationTree} {property} Children
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
    ## @deftp {CompactClassificationTree} {property} Parent
    ##
    ## Parent of each node
    ##
    ## A column vector naming the parent of each node.  The root carries a
    ## zero.  This property is read-only.
    ##
    ## @end deftp
    Parent = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} IsBranchNode
    ##
    ## Which nodes are branch nodes
    ##
    ## A logical column vector, true for each node that carries a split and
    ## false for each leaf.  This property is read-only.
    ##
    ## @end deftp
    IsBranchNode = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} CutPredictor
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
    ## @deftp {CompactClassificationTree} {property} CutPredictorIndex
    ##
    ## Index of the predictor each node cuts on
    ##
    ## A column vector holding, for each node, the column of @var{X} the node
    ## splits on, and zero at a leaf.  This property is read-only.
    ##
    ## @end deftp
    CutPredictorIndex = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} CutPoint
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
    ## @deftp {CompactClassificationTree} {property} CutType
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
    ## @deftp {CompactClassificationTree} {property} CutCategories
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
    ## @deftp {CompactClassificationTree} {property} CategoricalSplit
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
    ## @deftp {CompactClassificationTree} {property} NodeSize
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
    ## @deftp {CompactClassificationTree} {property} NodeClass
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
    ## @deftp {CompactClassificationTree} {property} NodeError
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
    ## @deftp {CompactClassificationTree} {property} NodeProbability
    ##
    ## Probability of reaching each node
    ##
    ## A column vector holding, for each node, the total weight of the
    ## observations that reached it, the weights being those the model it
    ## came from was fitted with.  The root carries one.  This property is
    ## read-only.
    ##
    ## @end deftp
    NodeProbability = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} NodeRisk
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
    ## @deftp {CompactClassificationTree} {property} ClassCount
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
    ## @deftp {CompactClassificationTree} {property} ClassProbability
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
    ## @deftp {CompactClassificationTree} {property} PruneList
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
    ## @deftp {CompactClassificationTree} {property} PruneAlpha
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
    ## @deftp {CompactClassificationTree} {property} SurrogateCutCategories
    ##
    ## Categories of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutCategories = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} SurrogateCutFlip
    ##
    ## Cut assignments of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutFlip = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} SurrogateCutPoint
    ##
    ## Cut points of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutPoint = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} SurrogateCutType
    ##
    ## Types of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutType = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} SurrogateCutPredictor
    ##
    ## Predictors of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutPredictor = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} SurrogatePredictorAssociation
    ##
    ## Predictive measures of association of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogatePredictorAssociation = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} PredictorNames
    ##
    ## Names of the predictor variables
    ##
    ## A cell array of character vectors with one name per column of
    ## @var{X}.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector naming the response.  This property is read-only.
    ##
    ## @end deftp
    ResponseName = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} ClassNames
    ##
    ## Names of the classes
    ##
    ## The distinct class labels, in the type the response was given in and
    ## sorted.  This property is read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} CategoricalPredictors
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
    ## @deftp {CompactClassificationTree} {property} ExpandedPredictorNames
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

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} Prior
    ##
    ## Prior probability of each class
    ##
    ## A numeric row vector with one element per class, summing to one.  It
    ## defaults to the weight each class carries in the training data, and
    ## may be reassigned after fitting.
    ##
    ## Reassigning it re-derives every node statistic that depends on the
    ## class weights, so @code{ClassProbability}, @code{NodeProbability},
    ## @code{NodeClass}, @code{NodeError} and @code{NodeRisk} all follow.
    ## The shape of the tree does not, having been decided by the prior in
    ## force when it was grown.  MATLAB reports this property read-only on a
    ## compact tree and refuses the assignment; it is settable here, as it is
    ## on the package's other compact classifiers.
    ##
    ## @end deftp
    Prior = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} Cost
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
    ## when it was grown.  MATLAB reports this property read-only on a
    ## compact tree and refuses the assignment; it is settable here, as it is
    ## on the package's other compact classifiers.
    ##
    ## @end deftp
    Cost = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationTree} {property} ScoreTransform
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

    ## The split criterion the tree was grown under, which the risk is
    ## measured by.  The full class reads it from ModelParameters, which a
    ## compact model does not carry.
    SplitCriterion = [];

    ## The share of each class's total weight that reached each node, one row
    ## per node and one column per class.  It is what the node statistics are
    ## re-derived from, and it is free of both the prior and the cost.
    ClassShare = [];

  endproperties

  ## Set methods for the properties a user may assign after compacting.
  methods (Hidden)

    function this = set.Cost (this, Cost)
      [C, errmsg] = costMatrix (Cost, this.ClassNames);
      if (! isempty (errmsg))
        error ("CompactClassificationTree.Cost: %s", errmsg);
      endif
      this.Cost = C;
      this = deriveNodes (this);
    endfunction

    function this = set.Prior (this, Prior)
      P = Prior;
      K = classCount (this.ClassNames);
      if (isstruct (P))
        P = priorFromStruct (P, this.ClassNames, ...
                             'CompactClassificationTree.Prior');
      elseif (ischar (P))
        if (strcmpi (P, 'uniform'))
          P = ones (1, K) / K;
        elseif (! strcmpi (P, 'empirical'))
          error (strcat ("CompactClassificationTree.Prior: a character", ...
                         " vector must be 'empirical' or 'uniform'."));
        else
          return;   # 'empirical' after fitting is what is already stored
        endif
      endif
      if (! (isnumeric (P) && isvector (P) && isreal (P)))
        error (strcat ("CompactClassificationTree.Prior: must be a real", ...
                       " numeric vector, a structure, 'empirical', or", ...
                       " 'uniform'."));
      endif
      if (numel (P) != K)
        error (strcat ("CompactClassificationTree.Prior: must have one", ...
                       " element per class."));
      endif
      if (any (P < 0) || ! (sum (P) > 0))
        error (strcat ("CompactClassificationTree.Prior: must be", ...
                       " nonnegative and must not be all zero."));
      endif
      this.Prior = P(:)' / sum (P);
      this = deriveNodes (this);
    endfunction

    function this = set.ScoreTransform (this, val)
      [f, st] = parseScoreTransform (val, 'CompactClassificationTree');
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
      fprintf ('\n  CompactClassificationTree\n\n');
      fprintf ('%22s: %s\n', 'ResponseName', this.ResponseName);
      fprintf ('%22s: %s\n', 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ('%22s: %s\n', 'ClassNames', classNameListing (this.ClassNames));
      fprintf ('%22s: %s\n', 'ScoreTransform', this.ScoreTransform);
      fprintf ('%22s: %d\n', 'NumNodes', this.NumNodes);
      fprintf ('\n');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationTree} {@var{obj} =} CompactClassificationTree (@var{Mdl})
    ##
    ## Create a @code{CompactClassificationTree} object.
    ##
    ## @var{Mdl} is the @code{ClassificationTree} object to compact.  The
    ## documented way to reach this constructor is the @code{compact} method.
    ##
    ## @seealso{ClassificationTree, fitctree}
    ## @end deftypefn
    function this = CompactClassificationTree (Mdl)

      ## Input validation
      if (nargin < 1)
        error (strcat ("CompactClassificationTree: too few input", ...
                       " arguments."));
      endif
      if (! isa (Mdl, 'ClassificationTree'))
        error (strcat ("CompactClassificationTree: MDL must be a", ...
                       " ClassificationTree object."));
      endif

      ## The node table and everything derived from it, then the three
      ## properties whose set methods re-derive the rest.  The criterion the
      ## risk is measured by travels with the tree, there being no
      ## ModelParameters here to look it up in.
      this.NumNodes = Mdl.NumNodes;
      this.Children = Mdl.Children;
      this.Parent = Mdl.Parent;
      this.CutPredictorIndex = Mdl.CutPredictorIndex;
      this.CutPoint = Mdl.CutPoint;
      this.NodeSize = Mdl.NodeSize;
      this.ClassCount = Mdl.ClassCount;
      this.PruneList = Mdl.PruneList;
      this.PruneAlpha = Mdl.PruneAlpha;
      this.ClassShare = Mdl.ClassShare;
      this.SplitCriterion = Mdl.ModelParameters.SplitCriterion;

      this.PredictorNames = Mdl.PredictorNames;
      this.ResponseName = Mdl.ResponseName;
      this.ClassNames = Mdl.ClassNames;
      this.CategoricalPredictors = Mdl.CategoricalPredictors;
      this.ExpandedPredictorNames = Mdl.ExpandedPredictorNames;

      this = fillCuts (this);
      this.Cost = Mdl.Cost;
      this.Prior = Mdl.Prior;
      this.ScoreTransform = Mdl.ScoreTransform;

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationTree} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationTree} {[@var{label}, @var{score}] =} predict (@dots{})
    ## @deftypefnx {CompactClassificationTree} {[@var{label}, @var{score}, @var{node}] =} predict (@dots{})
    ## @deftypefnx {CompactClassificationTree} {[@var{label}, @var{score}, @var{node}, @var{cnum}] =} predict (@dots{})
    ##
    ## Classify new data with a trained @code{CompactClassificationTree} object.
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
    ## @seealso{CompactClassificationTree, fitctree}
    ## @end deftypefn
    function [label, score, node, cnum] = predict (this, XC)

      ## Input validation
      if (nargin < 2)
        error ("CompactClassificationTree.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("CompactClassificationTree.predict: XC is empty.");
      endif
      if (! (isnumeric (XC) && isreal (XC) && ismatrix (XC)))
        error (strcat ("CompactClassificationTree.predict: XC must be", ...
                       " a real numeric matrix."));
      endif
      if (numel (this.PredictorNames) != columns (XC))
        error (strcat ("CompactClassificationTree.predict: XC must have", ...
                       " the same number of predictors as the trained", ...
                       " model."));
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
    ## @deftypefn {CompactClassificationTree} {@var{imp} =} predictorImportance (@var{obj})
    ##
    ## Estimate the importance of each predictor.
    ##
    ## @code{@var{imp} = predictorImportance (@var{obj})} returns a row
    ## vector with one element per predictor, holding the total drop in risk
    ## the splits on that predictor bring about, divided by the number of
    ## branch nodes.  A predictor the tree never splits on scores zero.
    ##
    ## The drop at a branch node is its @code{NodeRisk} less the risk of its
    ## two children, so a predictor that is chosen often, high up, and on
    ## nodes it separates well, scores highest.  The numbers are comparable
    ## between predictors of one tree and not between trees.
    ##
    ## @seealso{CompactClassificationTree, fitctree, NodeRisk}
    ## @end deftypefn
    function imp = predictorImportance (this)

      imp = zeros (1, numel (this.PredictorNames));
      branch = find (this.Children(:,1) > 0);
      if (isempty (branch))
        return;
      endif
      for ii = 1:numel (branch)
        b = branch(ii);
        kids = this.Children(b,:);
        drop = this.NodeRisk(b) - this.NodeRisk(kids(1)) ...
               - this.NodeRisk(kids(2));
        v = this.CutPredictorIndex(b);
        imp(v) = imp(v) + drop;
      endfor
      imp = imp / numel (branch);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationTree} {@var{r} =} nodeVariableRange (@var{obj}, @var{node})
    ##
    ## Range of each predictor at a node.
    ##
    ## @code{@var{r} = nodeVariableRange (@var{obj}, @var{node})} returns a
    ## structure with one field per predictor the path from the root to
    ## @var{node} cuts on, holding the two-element range of values that reach
    ## the node.  A predictor the path never cuts on is unconstrained and is
    ## left out, so the root gives a structure with no fields.
    ##
    ## @seealso{CompactClassificationTree, fitctree}
    ## @end deftypefn
    function r = nodeVariableRange (this, node)

      ## Input validation
      if (nargin < 2)
        error (strcat ("CompactClassificationTree.nodeVariableRange:", ...
                       " too few input arguments."));
      endif
      if (! (isnumeric (node) && isscalar (node) && isreal (node)
             && node >= 1 && node <= this.NumNodes && node == fix (node)))
        error (strcat ("CompactClassificationTree.nodeVariableRange:", ...
                       " NODE must be a positive integer no greater than", ...
                       " the number of nodes in the tree."));
      endif

      p = numel (this.PredictorNames);
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
    ## @deftypefn {CompactClassificationTree} {} view (@var{obj})
    ##
    ## Print the tree as text.
    ##
    ## @code{view (@var{obj})} prints one line per node: a branch node names
    ## the predictor it cuts on, the cut point, and the node each side leads
    ## to, and a leaf names the class it assigns.  A branch node's line ends
    ## with the class it would assign itself, which is the answer an
    ## observation missing that predictor gets.
    ##
    ## @seealso{CompactClassificationTree, fitctree}
    ## @end deftypefn
    function view (this)

      fprintf ("Decision tree for classification\n");
      for ii = 1:this.NumNodes
        if (this.Children(ii,1) == 0)
          fprintf ("%d  class = %s\n", ii, this.NodeClass{ii});
        else
          v = this.PredictorNames{this.CutPredictorIndex(ii)};
          c = num2str (this.CutPoint(ii));
          fprintf ("%d  if %s<%s then node %d elseif %s>=%s then node", ...
                   ii, v, c, this.Children(ii,1), v, c);
          fprintf (" %d else %s\n", this.Children(ii,2), this.NodeClass{ii});
        endif
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationTree} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin on new data.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns one
    ## margin per observation: the score the model gives the observation's
    ## true class, less the largest score it gives any other class.  A
    ## positive margin means the observation is classified correctly, and a
    ## larger one means it is classified more confidently.
    ##
    ## @seealso{CompactClassificationTree, edge, loss, predict}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Input validation
      if (nargin < 3)
        error ("CompactClassificationTree.margin: too few input arguments.");
      endif
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("CompactClassificationTree.margin: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("CompactClassificationTree.margin: number of rows", ...
                       " in X and Y must be equal."));
      endif

      [~, s] = predict (this, X);
      m = marginsOf (s, gY, 1);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationTree} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationTree} {@var{e} =} edge (@dots{}, @qcode{'Weights'}, @var{w})
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
    ## @seealso{CompactClassificationTree, margin, loss, predict}
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      ## Input validation
      if (nargin < 3)
        error ("CompactClassificationTree.edge: too few input arguments.");
      endif

      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       'CompactClassificationTree', 'edge');
      e = sum (W .* margin (this, X, Y));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationTree} {@var{l} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationTree} {@var{l} =} loss (@dots{}, @var{name}, @var{value})
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
    ## @seealso{CompactClassificationTree, margin, edge, predict}
    ## @end deftypefn
    function l = loss (this, X, Y, varargin)

      ## Input validation
      if (nargin < 3)
        error ("CompactClassificationTree.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactClassificationTree.loss: name-value", ...
                       " arguments must be in pairs."));
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
              error ("CompactClassificationTree.loss: invalid loss function.");
            endif
            LossFun = tolower (Value);
          case 'weights'
            if (! (isnumeric (Value) && isvector (Value)))
              error ("CompactClassificationTree.loss: invalid 'Weights'.");
            endif
            Weights = Value;
          otherwise
            error (strcat ("CompactClassificationTree.loss: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("CompactClassificationTree.loss: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("CompactClassificationTree.loss: number of rows in", ...
                       " X and Y must be equal."));
      endif
      if (isempty (Weights))
        w = ones (numel (gY), 1);
      else
        w = Weights(:);
        if (numel (w) != numel (gY))
          error (strcat ("CompactClassificationTree.loss: 'Weights' must", ...
                         " have one element per observation."));
        endif
      endif

      ## The weights are normalized within each class to that class's prior,
      ## so that a loss and an edge weight the classes the same way.
      w = priorNormalize (w, gY, this.Prior);

      [~, s] = predict (this, X);
      l = classificationLoss (LossFun, s, gY, w, this.Cost);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationTree} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a CompactClassificationTree model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## CompactClassificationTree object into an Octave binary file, the name
    ## of which is specified in @var{filename}, along with an extra variable,
    ## which defines the type of classification object these variables
    ## constitute.  Use @code{loadmodel} in order to load a classification
    ## object into Octave's workspace.
    ##
    ## @seealso{loadmodel, fitctree, ClassificationTree}
    ## @end deftypefn
    function savemodel (this, fname)

      ## Input validation
      if (nargin < 2)
        error (strcat ("CompactClassificationTree.savemodel: too few", ...
                       " input arguments."));
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("CompactClassificationTree.savemodel: FNAME must", ...
                       " be a character vector."));
      endif

      ## Generate variable for class name
      classdef_name = 'CompactClassificationTree';

      ## Create variables from model properties
      PredictorNames         = this.PredictorNames;
      ResponseName           = this.ResponseName;
      ClassNames             = this.ClassNames;
      CategoricalPredictors  = this.CategoricalPredictors;
      ExpandedPredictorNames = this.ExpandedPredictorNames;
      NumNodes               = this.NumNodes;
      Children               = this.Children;
      Parent                 = this.Parent;
      CutPredictorIndex      = this.CutPredictorIndex;
      CutPoint               = this.CutPoint;
      NodeSize               = this.NodeSize;
      ClassCount             = this.ClassCount;
      ClassShare             = this.ClassShare;
      SplitCriterion         = this.SplitCriterion;
      PruneList              = this.PruneList;
      PruneAlpha             = this.PruneAlpha;
      Prior                  = this.Prior;
      Cost                   = this.Cost;
      ScoreTransform         = this.ScoreTransform;
      STfun                  = this.STfun;

      ## The cut and node descriptions are not saved: every one of them is
      ## re-derived from the node table above, so writing them out would only
      ## make a stale copy possible.
      save ('-binary', fname, 'classdef_name', 'PredictorNames', ...
            'ResponseName', 'ClassNames', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'NumNodes', 'Children', 'Parent', ...
            'CutPredictorIndex', 'CutPoint', 'NodeSize', 'ClassCount', ...
            'ClassShare', 'SplitCriterion', 'PruneList', 'PruneAlpha', ...
            'Prior', 'Cost', 'ScoreTransform', 'STfun');

    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      ## The compact model is built from a full one and has no training data
      ## of its own, so the smallest fit the full class accepts is compacted
      ## and then filled property by property.
      mdl = CompactClassificationTree (ClassificationTree ([1; 2], [1; 2]));

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
          msg = strcat ("CompactClassificationTree.load_model: invalid", ...
                        " model in '%s'.");
          error (msg, filename);
        end_try_catch
      endfor

      ## The cut and node descriptions were not saved, being derived
      mdl = fillCuts (mdl);
      mdl = deriveNodes (mdl);

    endfunction

  endmethods

  methods (Access = private)

    ## Re-derive every node statistic that depends on the prior or the cost.
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
                         this.SplitCriterion, this.ClassNames);
      this.NodeProbability = S.NodeProbability;
      this.ClassProbability = S.ClassProbability;
      this.NodeError = S.NodeError;
      this.NodeClass = S.NodeClass;
      this.NodeRisk = S.NodeRisk;

    endfunction

    ## The descriptions of the cuts, derived from the node table.
    function this = fillCuts (this)

      S = treeCutInfo (this.CutPredictorIndex, this.PredictorNames);
      for [val, name] = S
        this.(name) = val;
      endfor

    endfunction

  endmethods

endclassdef

## Tests
%!test  # MATLAB parity: the surface a compact tree reports
%! load fisheriris
%! CMdl = compact (ClassificationTree (meas, species));
%! assert_equal (class (CMdl), 'CompactClassificationTree');
%! assert_equal (numel (properties (CMdl)), 33);
%! assert_equal (CMdl.NumNodes, 9);
%! assert_equal (CMdl.ClassNames, unique (species));
%! assert_equal (CMdl.Prior, [1/3, 1/3, 1/3], 1e-15);
%! assert_equal (CMdl.Cost, [0, 1, 1; 1, 0, 1; 1, 1, 0]);
%! assert_equal (CMdl.ResponseName, 'Y');
%! assert_equal (CMdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (CMdl.ScoreTransform, 'none');

%!test  # The compact tree carries the node table the model it came from has
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.Children, Mdl.Children);
%! assert_equal (CMdl.Parent, Mdl.Parent);
%! assert_equal (CMdl.NodeSize, Mdl.NodeSize);
%! assert_equal (CMdl.ClassCount, Mdl.ClassCount);
%! assert_equal (CMdl.NodeRisk, Mdl.NodeRisk, 1e-15);
%! assert_equal (CMdl.NodeError, Mdl.NodeError, 1e-15);
%! assert_equal (CMdl.NodeProbability, Mdl.NodeProbability, 1e-15);
%! assert_equal (CMdl.ClassProbability, Mdl.ClassProbability, 1e-15);
%! assert_equal (CMdl.NodeClass, Mdl.NodeClass);
%! assert_equal (CMdl.PruneList, Mdl.PruneList);
%! assert_equal (CMdl.PruneAlpha, Mdl.PruneAlpha, 1e-15);
%! assert_equal (CMdl.CutPredictor, Mdl.CutPredictor);
%! assert_equal (CMdl.CutType, Mdl.CutType);

%!test  # MATLAB parity: what a tree with no categories and no surrogates holds
%! load fisheriris
%! CMdl = compact (ClassificationTree (meas, species));
%! assert_equal (size (CMdl.CutCategories), [9, 2]);
%! assert_equal (size (CMdl.CategoricalSplit), [0, 0]);
%! assert_equal (size (CMdl.SurrogateCutPredictor), [0, 1]);
%! assert_equal (size (CMdl.SurrogateCutPoint), [0, 0]);
%! assert_equal (size (CMdl.SurrogatePredictorAssociation), [0, 0]);

%!test  # MATLAB parity: predict answers exactly as the full model does
%! load fisheriris
%! Mdl = ClassificationTree (meas, species);
%! CMdl = compact (Mdl);
%! [l1, s1, n1, c1] = predict (Mdl, meas);
%! [l2, s2, n2, c2] = predict (CMdl, meas);
%! assert_equal ({l1, s1, n1, c1}, {l2, s2, n2, c2});
%! assert_equal (predict (CMdl, meas([1, 60, 120], :)), ...
%!               {'setosa'; 'versicolor'; 'virginica'});

%!test  # MATLAB parity: loss, margin and edge on a compact tree
%! load fisheriris
%! CMdl = compact (ClassificationTree (meas, species));
%! assert_equal (loss (CMdl, meas, species), 0.02, 1e-14);
%! assert_equal (edge (CMdl, meas, species), 0.938357487922706, 1e-14);
%! m = margin (CMdl, meas([1, 60, 120], :), species([1, 60, 120]));
%! assert_equal (m', [1, 1, 1/3], 1e-14);
%! assert_equal (loss (CMdl, meas, species, 'LossFun', 'hinge'), ...
%!               0.0308212560386473, 1e-14);

%!test  # MATLAB parity: predictorImportance and nodeVariableRange
%! load fisheriris
%! CMdl = compact (ClassificationTree (meas, species));
%! assert_equal (predictorImportance (CMdl), ...
%!               [0, 0, 0.0907484567901233, 0.0682128958668813], 1e-14);
%! r = nodeVariableRange (CMdl, 8);
%! assert_equal (r.x3, [2.45, 4.95], 1e-14);
%! assert_equal (r.x4, [-Inf, 1.65], 1e-14);
%! assert_equal (fieldnames (nodeVariableRange (CMdl, 1)), cell (0, 1));

%!test  # MATLAB parity: the text form of a compact tree
%! load fisheriris
%! CMdl = compact (ClassificationTree (meas, species));
%! lines = strsplit (strtrim (evalc ('view (CMdl)')), "\n");
%! assert_equal (numel (lines), 10);
%! assert_equal (lines{2}, ...
%!   '1  if x3<2.45 then node 2 elseif x3>=2.45 then node 3 else setosa');
%! assert_equal (lines{10}, '9  class = virginica');

%!test  # The risk follows the criterion the tree it came from was grown under
%! load fisheriris
%! Mdl = ClassificationTree (meas, species, 'SplitCriterion', 'deviance');
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.NodeRisk, Mdl.NodeRisk, 1e-15);
%! assert_equal (CMdl.NodeRisk(1), 0.792481250360577, 1e-14);

%!test  # Reassigning Cost re-derives the node statistics, not the tree
%! load fisheriris
%! CMdl = compact (ClassificationTree (meas, species));
%! shape = {CMdl.Children, CMdl.NodeSize, CMdl.ClassCount};
%! CMdl.Cost = [0, 2, 8; 3, 0, 1; 5, 4, 0];
%! assert_equal ({CMdl.Children, CMdl.NodeSize, CMdl.ClassCount}, shape);
%! assert_equal (CMdl.NodeClass{1}, 'versicolor');
%! assert_equal (CMdl.NodeError(1), 2, 1e-14);

%!test  # Reassigning Prior re-derives the node statistics
%! load fisheriris
%! CMdl = compact (ClassificationTree (meas, species));
%! CMdl.Prior = [0.5, 0.25, 0.25];
%! assert_equal (CMdl.Prior, [0.5, 0.25, 0.25], 1e-14);
%! assert_equal (CMdl.ClassProbability(1,:), [0.5, 0.25, 0.25], 1e-14);
%! assert_equal (CMdl.NodeProbability(2), 0.5, 1e-14);
%! assert_equal (CMdl.NodeRisk(1), 0.625, 1e-14);

%!test  # The score transform travels with the compact model
%! load fisheriris
%! CMdl = compact (ClassificationTree (meas, species, ...
%!                                     'ScoreTransform', 'logit'));
%! assert_equal (CMdl.ScoreTransform, 'logit');
%! [~, score] = predict (CMdl, meas(1, :));
%! assert_equal (score, [0.731058578630005, 0.5, 0.5], 1e-14);

%!test  # A compact model saved and loaded answers exactly as it did
%! load fisheriris
%! CMdl = compact (ClassificationTree (meas, species));
%! fname = tempname ();
%! unwind_protect
%!   savemodel (CMdl, fname);
%!   New = loadmodel (fname);
%!   assert_equal (class (New), 'CompactClassificationTree');
%!   assert_equal (New.NumNodes, CMdl.NumNodes);
%!   assert_equal (New.NodeRisk, CMdl.NodeRisk, 1e-15);
%!   assert_equal (New.CutPredictor, CMdl.CutPredictor);
%!   assert_equal (predict (New, meas), predict (CMdl, meas));
%! unwind_protect_cleanup
%!   delete (fname);
%! end_unwind_protect

## Test input validation
%!error<CompactClassificationTree: too few input arguments.>
%! CompactClassificationTree ()
%!error<CompactClassificationTree: MDL must be a ClassificationTree object.>
%! CompactClassificationTree (5)
%!error<CompactClassificationTree.predict: too few input arguments.>
%! predict (compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2])))
%!error<CompactClassificationTree.predict: XC is empty.>
%! predict (compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2])), [])
%!error<CompactClassificationTree.predict: XC must be a real numeric matrix.>
%! predict (compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2])), 'a')
%!error<CompactClassificationTree.predict: XC must have the same number of predictors as the trained model.>
%! predict (compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2])), ...
%!          ones (2, 3))
%!error<CompactClassificationTree.margin: too few input arguments.>
%! margin (compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2])), ...
%!         ones (4, 2))
%!error<CompactClassificationTree.edge: too few input arguments.>
%! edge (compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2])), ones (4, 2))
%!error<CompactClassificationTree.loss: too few input arguments.>
%! loss (compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2])), ones (4, 2))
%!error<CompactClassificationTree.loss: invalid loss function.>
%! loss (compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2])), ...
%!       ones (4, 2), [1; 1; 2; 2], 'LossFun', 'x')
%!error<CompactClassificationTree.nodeVariableRange: NODE must be a positive integer no greater than the number of nodes in the tree.>
%! nodeVariableRange (compact (ClassificationTree (ones (4, 2), ...
%!                             [1; 1; 2; 2])), 99)
%!error<CompactClassificationTree.savemodel: too few input arguments.>
%! savemodel (compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2])))
%!error<CompactClassificationTree.savemodel: FNAME must be a character vector.>
%! savemodel (compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2])), 5)
%!error<CompactClassificationTree.Prior: must have one element per class.>
%! CMdl = compact (ClassificationTree (ones (4, 2), [1; 1; 2; 2]));
%! CMdl.Prior = [0.2, 0.3, 0.5];
