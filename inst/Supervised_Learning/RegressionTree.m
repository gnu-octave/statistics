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

classdef RegressionTree
  ## -*- texinfo -*-
  ## @deftp {statistics} RegressionTree
  ##
  ## Binary decision tree for regression
  ##
  ## The @code{RegressionTree} class implements a CART binary decision tree.
  ## Growth splits each node on the single predictor and cut point that lower
  ## the squared error of the response the most, and stops when a node is too
  ## small to be a parent, has no split leaving enough observations on both
  ## sides, or already accounts for all but a tolerance of the error the root
  ## carried.  The grown tree is then optionally reduced, first by merging
  ## the leaves that buy no accuracy and then by cost complexity pruning,
  ## which orders the branch nodes by how little error their subtrees remove
  ## and records that order so a subtree of any size can be recovered
  ## afterwards with @code{prune}.
  ##
  ## Create a @code{RegressionTree} object by using the @code{fitrtree}
  ## function or the class constructor.
  ##
  ## The fit is carried out by the compiled engine @code{treetrain} and
  ## predictions by @code{treepredict}, which the classification tree shares.
  ##
  ## An observation missing the predictor a node cuts on descends to neither
  ## child.  It is counted in that node and in every node above it, and
  ## @code{predict} stops it there and gives it that node's answer, so a row
  ## is never sent down a branch on evidence it does not carry.
  ##
  ## @strong{What this class does not do yet.}  Categorical predictors,
  ## surrogate splits and predictor subsampling are not implemented, and an
  ## option asking for one of them is refused rather than quietly ignored.
  ## @code{CategoricalSplit}, @code{CutCategories} and the six
  ## @code{Surrogate} properties are therefore always empty, as they are in
  ## MATLAB on numeric data.
  ##
  ## @seealso{fitrtree, ClassificationTree, treetrain, treepredict}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} X
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
    ## @deftp {RegressionTree} {property} Y
    ##
    ## Response data
    ##
    ## A numeric column vector with one element per row of @var{X}, holding
    ## the observed response of each observation.  This property is
    ## read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} NumObservations
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
    ## @deftp {RegressionTree} {property} RowsUsed
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
    ## @deftp {RegressionTree} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector of the weights the fit used, one per retained
    ## observation.  They are the weights given, scaled to sum to one.  This
    ## property is read-only.
    ##
    ## @end deftp
    W = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} PredictorNames
    ##
    ## Names of the predictor variables
    ##
    ## A cell array of character vectors with one name per column of
    ## @var{X}.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector naming the response.  This property is read-only.
    ##
    ## @end deftp
    ResponseName = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} CategoricalPredictors
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
    ## @deftp {RegressionTree} {property} ExpandedPredictorNames
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
    ## @deftp {RegressionTree} {property} BinEdges
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
    ## @deftp {RegressionTree} {property} ModelParameters
    ##
    ## Parameters the fit was run with
    ##
    ## A structure recording the options the tree was grown under:
    ## @qcode{SplitCriterion}, @qcode{MinParent}, @qcode{MinLeaf},
    ## @qcode{MaxSplits}, @qcode{NVarToSample}, @qcode{MergeLeaves},
    ## @qcode{Prune}, @qcode{PruneCriterion}, @qcode{QEToler},
    ## @qcode{NSurrogate}, @qcode{MaxCat}, @qcode{AlgCat},
    ## @qcode{PredictorSelection}, @qcode{Method} and @qcode{Type}.
    ## @qcode{SplitCriterion} and @qcode{PruneCriterion} are both
    ## @qcode{'mse'}, the only criterion a regression tree has, and
    ## @qcode{QEToler} is the tolerance growth stops at.
    ##
    ## @qcode{MinParent} is the value the fit used, which is
    ## @code{max (MinParentSize, 2 * MinLeafSize)} and so may exceed the
    ## @qcode{'MinParentSize'} asked for.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} HyperparameterOptimizationResults
    ##
    ## Results of a hyperparameter optimization
    ##
    ## Hyperparameter optimization is not implemented, so this is always
    ## empty.  This property is read-only.
    ##
    ## @end deftp
    HyperparameterOptimizationResults = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} NumNodes
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
    ## @deftp {RegressionTree} {property} Children
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
    ## @deftp {RegressionTree} {property} Parent
    ##
    ## Parent of each node
    ##
    ## A column vector naming the parent of each node.  The root carries a
    ## zero.  This property is read-only.
    ##
    ## @end deftp
    Parent = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} IsBranchNode
    ##
    ## Which nodes are branch nodes
    ##
    ## A logical column vector, true for each node that carries a split and
    ## false for each leaf.  This property is read-only.
    ##
    ## @end deftp
    IsBranchNode = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} CutPredictor
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
    ## @deftp {RegressionTree} {property} CutPredictorIndex
    ##
    ## Index of the predictor each node cuts on
    ##
    ## A column vector holding, for each node, the column of @var{X} the node
    ## splits on, and zero at a leaf.  This property is read-only.
    ##
    ## @end deftp
    CutPredictorIndex = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} CutPoint
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
    ## @deftp {RegressionTree} {property} CutType
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
    ## @deftp {RegressionTree} {property} CutCategories
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
    ## @deftp {RegressionTree} {property} CategoricalSplit
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
    ## @deftp {RegressionTree} {property} NodeSize
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
    ## @deftp {RegressionTree} {property} NodeMean
    ##
    ## Mean response at each node
    ##
    ## A column vector holding, for each node, the weighted mean of the
    ## response over the observations that reached it.  It is what
    ## @code{predict} answers for a row that comes to rest there.  This
    ## property is read-only.
    ##
    ## @end deftp
    NodeMean = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} NodeError
    ##
    ## Mean squared error of each node
    ##
    ## A column vector holding, for each node, the weighted mean squared
    ## error of the response about the node's mean.  This property is
    ## read-only.
    ##
    ## @end deftp
    NodeError = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} NodeProbability
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
    ## @deftp {RegressionTree} {property} NodeRisk
    ##
    ## Risk of each node
    ##
    ## A column vector holding, for each node, its mean squared error
    ## weighted by the probability of reaching it, which is the squared error
    ## the node contributes to the whole tree.  This property is read-only.
    ##
    ## @end deftp
    NodeRisk = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} PruneList
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
    ## @deftp {RegressionTree} {property} PruneAlpha
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
    ## @deftp {RegressionTree} {property} SurrogateCutCategories
    ##
    ## Categories of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutCategories = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} SurrogateCutFlip
    ##
    ## Cut assignments of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutFlip = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} SurrogateCutPoint
    ##
    ## Cut points of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutPoint = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} SurrogateCutType
    ##
    ## Types of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutType = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} SurrogateCutPredictor
    ##
    ## Predictors of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutPredictor = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionTree} {property} SurrogatePredictorAssociation
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
    ## @deftp {RegressionTree} {property} ResponseTransform
    ##
    ## Transform applied to the predicted response
    ##
    ## A character vector naming the function @code{predict} applies to the
    ## response it predicts before returning it, or a function handle taking
    ## and returning an array of the same size.  The default is
    ## @qcode{'none'}.
    ##
    ## @end deftp
    ResponseTransform = [];

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## The parsed ResponseTransform, applied to the response by predict.
    RTfun = [];

    ## The total weight of the observations that reached each node, before it
    ## was divided by the root's to give NodeProbability.
    NodeWeight = [];

    ## The risk each branch node carries on account of the observations that
    ## stop there, missing the predictor it cuts on.  Those observations are
    ## in neither child, so a subtree's risk is its children's plus this.  It
    ## is zero at a leaf and zero throughout a tree fitted on data with no
    ## missing values.
    HeldRisk = [];

    ## The observation weights as they were given, before they were scaled to
    ## sum to one.  A fold of a cross-validated model is grown with a slice
    ## of these.
    RawWeights = [];

  endproperties

  ## Set methods for the properties a user may assign after fitting.
  methods (Hidden)

    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = ...
                            parseResponseTransform (val, 'RegressionTree');
    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ('\n  RegressionTree\n\n');
      fprintf ('%22s: %s\n', 'ResponseName', this.ResponseName);
      fprintf ('%22s: %s\n', 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ('%22s: %s\n', 'ResponseTransform', this.ResponseTransform);
      fprintf ('%22s: %d\n', 'NumObservations', this.NumObservations);
      fprintf ('%22s: %d\n', 'NumNodes', this.NumNodes);
      fprintf ('\n');
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionTree} {@var{obj} =} RegressionTree (@var{X}, @var{Y})
    ## @deftypefnx {RegressionTree} {@var{obj} =} RegressionTree (@dots{}, @var{name}, @var{value})
    ##
    ## Grow a binary decision tree for regression.
    ##
    ## @code{@var{obj} = RegressionTree (@var{X}, @var{Y})} grows a tree on
    ## the @math{NxP} numeric matrix @var{X} of predictor data and the
    ## @math{Nx1} numeric response @var{Y}, and returns it as a
    ## @code{RegressionTree} object.
    ##
    ## @code{@var{obj} = RegressionTree (@dots{}, @var{name}, @var{value})}
    ## takes the options below.
    ##
    ## @multitable @columnfractions 0.24 0.74
    ## @headitem @var{Name} @tab @var{Value}
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
    ## @item @qcode{'Prune'} @tab @qcode{'on'} (default) or @qcode{'off'}.
    ## When on, the cost complexity pruning sequence is estimated and
    ## reported in @code{PruneList} and @code{PruneAlpha}.  The tree returned
    ## is the unpruned one either way; @code{prune} takes a subtree out of
    ## the sequence.
    ##
    ## @item @qcode{'PruneCriterion'} @tab @qcode{'mse'}, the only criterion
    ## a regression tree has.
    ##
    ## @item @qcode{'QuadraticErrorTolerance'} @tab A positive scalar.  A
    ## node whose squared error has fallen to this fraction of the root's is
    ## not split further.  The default is 1e-6.
    ##
    ## @item @qcode{'ResponseName'} @tab A character vector naming the
    ## response.  The default is @qcode{'Y'}.
    ##
    ## @item @qcode{'ResponseTransform'} @tab A character vector naming a
    ## transform to apply to the predicted response, or a function handle.
    ## The default is @qcode{'none'}.
    ##
    ## @item @qcode{'SplitCriterion'} @tab @qcode{'mse'}, the only criterion
    ## a regression tree has.
    ##
    ## @item @qcode{'Weights'} @tab A nonnegative numeric vector with one
    ## element per observation.  The default is uniform.
    ##
    ## @end multitable
    ##
    ## @seealso{fitrtree, ClassificationTree, treetrain, treepredict}
    ## @end deftypefn
    function this = RegressionTree (X, Y, varargin)

      ## Input validation
      if (nargin < 2)
        error ("RegressionTree: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionTree: name-value arguments must be in", ...
                       " pairs."));
      endif
      if (! (isnumeric (X) && isreal (X) && ismatrix (X) && ! isempty (X)))
        error (strcat ("RegressionTree: X must be a non-empty real numeric", ...
                       " matrix."));
      endif
      if (! (isnumeric (Y) && isreal (Y) && isvector (Y)))
        error ("RegressionTree: Y must be a real numeric vector.");
      endif
      if (rows (X) != numel (Y))
        error (strcat ("RegressionTree: number of rows in X and Y must be", ...
                       " equal."));
      endif

      Y = double (Y(:));

      ## Defaults.  MaxNumSplits is left empty until the retained rows are
      ## known, its default being one less than their number.
      PredictorNames = {};
      ResponseName   = [];
      Weights        = [];
      MaxNumSplits   = [];
      MergeLeaves    = 'on';
      MinLeafSize    = 1;
      MinParentSize  = 10;
      Prune          = 'on';
      QEToler        = 1e-6;
      this.ResponseTransform = 'none';

      ## Parse optional parameters
      while (numel (varargin) > 0)
        Value = varargin{2};
        switch (tolower (varargin{1}))

          case 'predictornames'
            PredictorNames = Value;
            if (! iscellstr (PredictorNames))
              error (strcat ("RegressionTree: 'PredictorNames' must be", ...
                             " supplied as a cellstring array."));
            elseif (numel (PredictorNames) != columns (X))
              error (strcat ("RegressionTree: 'PredictorNames' must equal", ...
                             " the number of columns in X."));
            endif

          case 'responsename'
            ResponseName = Value;
            if (! (ischar (ResponseName) && isrow (ResponseName)))
              error (strcat ("RegressionTree: 'ResponseName' must be a", ...
                             " character vector."));
            endif

          case 'responsetransform'
            this.ResponseTransform = Value;

          case 'weights'
            Weights = Value;
            if (! (isnumeric (Weights) && isvector (Weights)
                   && isreal (Weights)))
              error (strcat ("RegressionTree: 'Weights' must be a real", ...
                             " numeric vector."));
            endif
            if (numel (Weights) != rows (X))
              error (strcat ("RegressionTree: 'Weights' must have one", ...
                             " element per row in X."));
            endif
            if (any (Weights < 0) || ! (sum (Weights) > 0))
              error (strcat ("RegressionTree: 'Weights' must be", ...
                             " nonnegative and must not be all zero."));
            endif

          case 'maxnumsplits'
            MaxNumSplits = Value;
            if (! (isnumeric (MaxNumSplits) && isscalar (MaxNumSplits)
                   && isreal (MaxNumSplits) && MaxNumSplits >= 0
                   && MaxNumSplits == fix (MaxNumSplits)))
              error (strcat ("RegressionTree: 'MaxNumSplits' must be a", ...
                             " nonnegative integer."));
            endif

          case 'minleafsize'
            MinLeafSize = Value;
            if (! (isnumeric (MinLeafSize) && isscalar (MinLeafSize)
                   && isreal (MinLeafSize) && MinLeafSize >= 1
                   && MinLeafSize == fix (MinLeafSize)))
              error (strcat ("RegressionTree: 'MinLeafSize' must be a", ...
                             " positive integer."));
            endif

          case 'minparentsize'
            MinParentSize = Value;
            if (! (isnumeric (MinParentSize) && isscalar (MinParentSize)
                   && isreal (MinParentSize) && MinParentSize >= 1
                   && MinParentSize == fix (MinParentSize)))
              error (strcat ("RegressionTree: 'MinParentSize' must be a", ...
                             " positive integer."));
            endif

          case 'mergeleaves'
            MergeLeaves = Value;
            if (! (ischar (MergeLeaves)
                   && any (strcmpi (MergeLeaves, {'on', 'off'}))))
              error (strcat ("RegressionTree: 'MergeLeaves' must be either", ...
                             " 'on' or 'off'."));
            endif

          case 'prune'
            Prune = Value;
            if (! (ischar (Prune) && any (strcmpi (Prune, {'on', 'off'}))))
              error (strcat ("RegressionTree: 'Prune' must be either 'on'", ...
                             " or 'off'."));
            endif

          case 'quadraticerrortolerance'
            QEToler = Value;
            if (! (isnumeric (QEToler) && isscalar (QEToler)
                   && isreal (QEToler) && QEToler > 0))
              error (strcat ("RegressionTree: 'QuadraticErrorTolerance'", ...
                             " must be a positive scalar."));
            endif

          ## A regression tree has one criterion either way, so the two
          ## names are taken and checked rather than refused outright.
          case {'splitcriterion', 'prunecriterion'}
            if (! (ischar (Value) && strcmpi (Value, 'mse')))
              error (strcat ("RegressionTree: '%s' must be 'mse' for a", ...
                             " regression tree."), varargin{1});
            endif

          case 'categoricalpredictors'
            ## Accepted only when it asks for nothing, so that a caller
            ## passing the MATLAB default is not turned away.
            if (! (isempty (Value)
                   || (ischar (Value) && strcmpi (Value, 'none'))))
              error (strcat ("RegressionTree: categorical predictors are", ...
                             " not implemented."));
            endif

          ## Options MATLAB takes that this class does not implement.  They
          ## are named one by one so that asking for one is refused rather
          ## than quietly doing nothing.
          case {'surrogate', 'numvariablestosample', 'predictorselection', ...
                'algorithmforcategorical', 'maxnumcategories', 'numbins', ...
                'optimizehyperparameters', ...
                'hyperparameteroptimizationoptions'}
            error ("RegressionTree: '%s' is not implemented.", varargin{1});

          case {'crossval', 'cvpartition', 'holdout', 'kfold', 'leaveout'}
            error (strcat ("RegressionTree: '%s' is not implemented; fit", ...
                           " the model and cross-validate it afterwards."), ...
                   varargin{1});

          otherwise
            error (strcat ("RegressionTree: invalid parameter name in", ...
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

      ## An observation is dropped only when its response is missing.  A row
      ## whose predictors hold missing values is kept and reported as used;
      ## the fit sends it as far down the tree as the predictors it does
      ## carry allow.
      RowsUsed = ! isnan (Y);
      Y = Y(RowsUsed);
      X = X(RowsUsed, :);
      if (! isempty (Weights))
        Weights = Weights(RowsUsed);
      endif
      if (isempty (Y))
        error ("RegressionTree: no observations with a known response.");
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

      ## The raw weights, kept so that a fold can be given a slice of them
      if (isempty (Weights))
        RawWeights = ones (this.NumObservations, 1);
      else
        RawWeights = double (Weights(:));
      endif
      this.RawWeights = RawWeights;

      ## A regression carries no prior, so the weights are simply scaled to
      ## sum to one, which is the W MATLAB reports.
      this.W = RawWeights / sum (RawWeights);

      if (isempty (MaxNumSplits))
        MaxNumSplits = max (this.NumObservations - 1, 0);
      endif
      ## A node smaller than two leaves cannot be split whatever the parent
      ## size asked for, so the fit uses the larger of the two.
      MinParent = max (MinParentSize, 2 * MinLeafSize);
      mergeOn = strcmpi (MergeLeaves, 'on');
      pruneOn = strcmpi (Prune, 'on');

      opts = struct ('NumClasses', 1, ...
                     'MinParent', MinParent, ...
                     'MinLeaf', MinLeafSize, ...
                     'MaxSplits', MaxNumSplits, ...
                     'SplitCriterion', 'mse', ...
                     'MergeLeaves', mergeOn, ...
                     'QEToler', QEToler, ...
                     ## The sequence is estimated below instead, MATLAB
                     ## carrying one whenever leaves were merged whether or
                     ## not pruning was asked for, where the engine builds
                     ## one only when it prunes.
                     'Prune', false);

      T = treetrain (X, this.Y, this.W, opts);

      ## The node table the engine returns
      this.NumNodes = T.NumNodes;
      this.Children = T.Children;
      this.Parent = T.Parent;
      this.CutPredictorIndex = T.CutPredictorIndex;
      this.CutPoint = T.CutPoint;
      this.NodeSize = T.NodeSize;
      this.NodeWeight = T.NodeWeight;
      this.NodeMean = T.NodeMean;
      this.NodeError = T.NodeError;

      this = fillCuts (this);
      this = deriveNodes (this);
      this = deriveHeld (this);

      ## Merging leaves is the first step of the cost complexity sequence, so
      ## a merged tree carries one whether or not pruning was asked for.
      ## Measured on R2024a: only turning both off leaves the two pruning
      ## properties empty.
      if (pruneOn || mergeOn)
        this = pruneSequence (this);
      else
        this.PruneList = [];
        this.PruneAlpha = [];
      endif

      this.ModelParameters = struct ('SplitCriterion', 'mse', ...
                                     'MinParent', MinParent, ...
                                     'MinLeaf', MinLeafSize, ...
                                     'MaxSplits', MaxNumSplits, ...
                                     'NVarToSample', 'all', ...
                                     'MergeLeaves', tolower (MergeLeaves), ...
                                     'Prune', tolower (Prune), ...
                                     'PruneCriterion', 'mse', ...
                                     'QEToler', QEToler, ...
                                     'NSurrogate', 0, ...
                                     'MaxCat', 10, ...
                                     'AlgCat', 'auto', ...
                                     'PredictorSelection', 'allsplits', ...
                                     'Method', 'Tree', ...
                                     'Type', 'regression');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionTree} {@var{yFit} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {RegressionTree} {[@var{yFit}, @var{node}] =} predict (@dots{})
    ##
    ## Predict the response with a trained @code{RegressionTree} object.
    ##
    ## @code{@var{yFit} = predict (@var{obj}, @var{XC})} sends each row of
    ## @var{XC} down the tree and returns the mean response of the node it
    ## comes to rest at, after @code{ResponseTransform}.  @var{XC} must have
    ## as many columns as the predictor data the model was fitted on.
    ##
    ## @code{[@var{yFit}, @var{node}] = predict (@dots{})} also returns the
    ## number of the node each row landed in.
    ##
    ## A row missing the predictor a node cuts on is stopped at that node and
    ## takes its answer, rather than being sent down a branch on evidence the
    ## row does not carry.
    ##
    ## @seealso{RegressionTree, fitrtree}
    ## @end deftypefn
    function [yFit, node] = predict (this, XC)

      ## Input validation
      if (nargin < 2)
        error ("RegressionTree.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("RegressionTree.predict: XC is empty.");
      endif
      if (! (isnumeric (XC) && isreal (XC) && ismatrix (XC)))
        error (strcat ("RegressionTree.predict: XC must be a real numeric", ...
                       " matrix."));
      endif
      if (numel (this.PredictorNames) != columns (XC))
        error (strcat ("RegressionTree.predict: XC must have the same", ...
                       " number of predictors as the trained model."));
      endif

      [yFit, node] = treepredict (XC, this.Children, ...
                                  this.CutPredictorIndex, this.CutPoint, ...
                                  this.NodeMean);
      yFit = this.RTfun (yFit);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionTree} {@var{obj2} =} prune (@var{obj})
    ## @deftypefnx {RegressionTree} {@var{obj2} =} prune (@var{obj}, @qcode{'Level'}, @var{L})
    ## @deftypefnx {RegressionTree} {@var{obj2} =} prune (@var{obj}, @qcode{'Alpha'}, @var{A})
    ## @deftypefnx {RegressionTree} {@var{obj2} =} prune (@var{obj}, @qcode{'Nodes'}, @var{N})
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
    ## @seealso{RegressionTree, fitrtree, RegressionTree.PruneList,
    ## RegressionTree.PruneAlpha}
    ## @end deftypefn
    function this = prune (this, varargin)

      ## Input validation
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionTree.prune: name-value arguments must", ...
                       " be in pairs."));
      endif
      if (numel (varargin) > 2)
        error (strcat ("RegressionTree.prune: specify only one of the", ...
                       " optional name-value paired arguments."));
      endif

      nodes = [];
      if (numel (varargin) == 2)
        Value = varargin{2};
        switch (tolower (varargin{1}))

          case 'level'
            if (! (isnumeric (Value) && isscalar (Value) && isreal (Value)
                   && Value >= 0 && Value == fix (Value)))
              error (strcat ("RegressionTree.prune: 'Level' must be a", ...
                             " nonnegative integer."));
            endif
            nodes = this.nodesAtLevel (Value, 'prune');

          case 'alpha'
            if (! (isnumeric (Value) && isscalar (Value) && isreal (Value)
                   && Value >= 0))
              error (strcat ("RegressionTree.prune: 'Alpha' must be a", ...
                             " nonnegative scalar."));
            endif
            if (isempty (this.PruneAlpha))
              error (strcat ("RegressionTree.prune: the tree carries no", ...
                             " pruning sequence."));
            endif
            ## The largest level whose parameter the value reaches
            L = sum (this.PruneAlpha <= Value) - 1;
            nodes = this.nodesAtLevel (max (L, 0), 'prune');

          case 'nodes'
            if (! (isnumeric (Value) && isreal (Value) && isvector (Value)
                   && all (Value >= 1) && all (Value <= this.NumNodes)
                   && all (Value == fix (Value))))
              error (strcat ("RegressionTree.prune: 'Nodes' must hold", ...
                             " indices of nodes of the tree."));
            endif
            nodes = Value(:);

          otherwise
            error (strcat ("RegressionTree.prune: invalid parameter name", ...
                           " in optional pair arguments."));

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
    ## @deftypefn {RegressionTree} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Drop the training data from a trained model.
    ##
    ## @code{@var{CMdl} = compact (@var{obj})} returns a
    ## @code{CompactRegressionTree} object carrying the tree and everything
    ## @code{predict} needs, but not the observations the model was fitted
    ## on.  It answers new data identically and is far smaller to keep or to
    ## ship.
    ##
    ## @seealso{CompactRegressionTree, RegressionTree}
    ## @end deftypefn
    function CMdl = compact (this)

      CMdl = CompactRegressionTree (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionTree} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {RegressionTree} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate a trained decision tree.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} partitions the training data
    ## into ten folds, or into as many folds as there are observations when
    ## there are fewer than ten, grows a tree on the training part of each
    ## and returns them as a @code{RegressionPartitionedModel}.
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
    ## Every fold is grown with the growth parameters the parent was grown
    ## with and a slice of its observation weights.
    ##
    ## @seealso{RegressionPartitionedModel, RegressionTree, cvpartition}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)

      ## Input validation
      if (numel (varargin) == 1)
        error (strcat ("RegressionTree.crossval: Name-Value arguments", ...
                       " must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("RegressionTree.crossval: specify only one of the", ...
                       " optional Name-Value paired arguments."));
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
              error (strcat ("RegressionTree.crossval: 'KFold' must be an", ...
                             " integer value greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("RegressionTree.crossval: 'Holdout' must be", ...
                             " a numeric value between 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && any (strcmpi (Leaveout, {'on', 'off'}))))
              error (strcat ("RegressionTree.crossval: 'Leaveout' must be", ...
                             " either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("RegressionTree.crossval: 'CVPartition' must", ...
                             " be a 'cvpartition' object."));
            endif

          otherwise
            error (strcat ("RegressionTree.crossval: invalid parameter", ...
                           " name in optional paired arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## The partition covers the observations actually trained on: a row
      ## dropped for a missing response is not one the folds can use.  A
      ## regression has no classes to stratify over, so the partition is
      ## drawn over their number.
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (this.NumObservations, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (this.NumObservations, 'LeaveOut');
      else
        partition = cvpartition (this.NumObservations, 'KFold', numFolds);
      endif

      CVMdl = RegressionPartitionedModel (this, partition);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionTree} {@var{E} =} cvloss (@var{obj})
    ## @deftypefnx {RegressionTree} {[@var{E}, @var{SE}, @var{Nleaf}, @var{BestLevel}] =} cvloss (@var{obj})
    ## @deftypefnx {RegressionTree} {[@dots{}] =} cvloss (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validated loss of a tree and of its subtrees.
    ##
    ## @code{@var{E} = cvloss (@var{obj})} partitions the training data into
    ## ten folds, grows a tree on the training part of each, and returns the
    ## mean squared error of the held-out part.
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
    ## tree.  The last level takes every fold's tree back to its root.  The
    ## partition is drawn over the observations rather than over a response
    ## there is nothing to stratify, and the loss is weighed by the model's
    ## own weights.
    ##
    ## @strong{The standard error is not MATLAB's.}  This is the standard
    ## error of the loss over the folds, which is what the name means.  Its
    ## value is not MATLAB's, whose formula is not recoverable from what it
    ## reports; @var{E}, @var{Nleaf} and @var{BestLevel} are measured and
    ## match.
    ##
    ## @seealso{RegressionTree, RegressionTree.prune, RegressionTree.crossval,
    ## RegressionTree.loss}
    ## @end deftypefn
    function [E, SE, Nleaf, BestLevel] = cvloss (this, varargin)

      ## Input validation
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionTree.cvloss: name-value arguments must", ...
                       " be in pairs."));
      endif
      if (isempty (this.PruneAlpha))
        error (strcat ("RegressionTree.cvloss: the tree carries no pruning", ...
                       " sequence; fit it with 'Prune' or 'MergeLeaves'", ...
                       " on."));
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
              error (strcat ("RegressionTree.cvloss: 'SubTrees' must be", ...
                             " 'all' or a vector of nonnegative integers", ...
                             " in ascending order."));
            endif

          case 'treesize'
            if (! (ischar (Value) && any (strcmpi (Value, {'se', 'min'}))))
              error (strcat ("RegressionTree.cvloss: 'TreeSize' must be", ...
                             " either 'se' or 'min'."));
            endif
            TreeSize = tolower (Value);

          case 'kfold'
            KFold = Value;
            if (! (isnumeric (KFold) && isscalar (KFold) && isreal (KFold)
                   && KFold == fix (KFold) && KFold > 1))
              error (strcat ("RegressionTree.cvloss: 'KFold' must be an", ...
                             " integer value greater than 1."));
            endif

          otherwise
            error (strcat ("RegressionTree.cvloss: invalid parameter name", ...
                           " in optional pair arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      if (any (SubTrees > maxLevel))
        error (strcat ("RegressionTree.cvloss: 'SubTrees' must not exceed", ...
                       " the largest pruning level, %d."), maxLevel);
      endif
      if (KFold > this.NumObservations)
        warning (strcat ("RegressionTree.cvloss: 'KFold' is greater than", ...
                         " the number of observations and is reduced to", ...
                         " %d."), this.NumObservations);
        KFold = this.NumObservations;
      endif

      ## The complexity parameter each subtree is matched to in a fold: the
      ## geometric mean of the two the parent's own sequence brackets it
      ## with, and infinity for the last, which takes a fold back to its
      ## root.  Measured on R2024a, which the classification tree measured
      ## first and which holds here too: matching a fold by level index or
      ## by the parent's parameter itself both give different answers.
      alpha = this.PruneAlpha(:)';
      ## Formed outside the brackets: inside them the space before the paren
      ## would split the call off into an element of its own.
      geo = sqrt (alpha(1:end-1) .* alpha(2:end));
      geo(end+1) = Inf;

      nsub = numel (SubTrees);
      n = this.NumObservations;
      partition = cvpartition (n, 'KFold', KFold);
      args = treeFoldArgs (this);

      ## The squared error of every observation, from the fold that did not
      ## train on it
      L = zeros (n, nsub);
      fold = zeros (n, 1);
      for k = 1:KFold
        tr = training (partition, k);
        te = test (partition, k);
        fold(te) = k;
        ft = RegressionTree (this.X(tr,:), this.Y(tr), args{:}, ...
                             'Weights', this.RawWeights(tr));
        fa = ft.PruneAlpha(:)';
        for j = 1:nsub
          lv = sum (fa <= geo(SubTrees(j) + 1)) - 1;
          ## The nodes are collapsed without re-estimating the fold's own
          ## sequence, which prune would do and which nothing here reads.
          sub = collapseNodes (ft, nodesAtLevel (ft, max (lv, 0), 'cvloss'));
          ## The fold carries no transform, so the parent's is applied here
          ## rather than twice over.
          yf = this.RTfun (predict (sub, this.X(te,:)));
          L(te,j) = (this.Y(te) - yf) .^ 2;
        endfor
      endfor

      ## The loss is weighed by the model's own weights, which sum to one.
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
    ## @deftypefn {RegressionTree} {@var{imp} =} predictorImportance (@var{obj})
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
    ## chosen often, high up, and on nodes it separates well, scores highest.
    ## The numbers are comparable between predictors of one tree and not
    ## between trees.
    ##
    ## @seealso{RegressionTree, fitrtree, RegressionTree.NodeRisk}
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
        ## The rows that stop at the node are in neither child and are
        ## answered the same whether it is split or not, so they are no part
        ## of what the split buys.
        drop = this.NodeRisk(b) - this.NodeRisk(kids(1)) ...
               - this.NodeRisk(kids(2)) - this.HeldRisk(b);
        v = this.CutPredictorIndex(b);
        imp(v) = imp(v) + drop;
      endfor
      imp = imp / numel (branch);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionTree} {@var{r} =} nodeVariableRange (@var{obj}, @var{node})
    ##
    ## Range of each predictor at a node.
    ##
    ## @code{@var{r} = nodeVariableRange (@var{obj}, @var{node})} returns a
    ## structure with one field per predictor the path from the root to
    ## @var{node} cuts on, holding the two-element range of values that reach
    ## the node.  A predictor the path never cuts on is unconstrained and is
    ## left out, so the root gives a structure with no fields.
    ##
    ## @seealso{RegressionTree, fitrtree}
    ## @end deftypefn
    function r = nodeVariableRange (this, node)

      ## Input validation
      if (nargin < 2)
        error (strcat ("RegressionTree.nodeVariableRange: too few input", ...
                       " arguments."));
      endif
      if (! (isnumeric (node) && isscalar (node) && isreal (node)
             && node >= 1 && node <= this.NumNodes && node == fix (node)))
        error (strcat ("RegressionTree.nodeVariableRange: NODE must be a", ...
                       " positive integer no greater than the number of", ...
                       " nodes in the tree."));
      endif

      p = numel (this.PredictorNames);
      lo = -Inf (1, p);
      hi = Inf (1, p);
      touched = false (1, p);

      ## Walk up to the root, narrowing the range of whichever predictor
      ## each step cut on.  Going up rather than down finds the path without
      ## a search, a node having exactly one parent.
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
    ## @deftypefn {RegressionTree} {} view (@var{obj})
    ##
    ## Print the tree as text.
    ##
    ## @code{view (@var{obj})} prints one line per node: a branch node names
    ## the predictor it cuts on, the cut point, and the node each side leads
    ## to, and a leaf names the response it fits.  A branch node's line ends
    ## with the response it would fit itself, which is the answer an
    ## observation missing that predictor gets.
    ##
    ## @seealso{RegressionTree, fitrtree}
    ## @end deftypefn
    function view (this)

      ## The node numbers are right aligned on the widest of them, which is
      ## how MATLAB lays the listing out.
      w = numel (sprintf ("%d", this.NumNodes));
      fprintf ("Decision tree for regression\n");
      for ii = 1:this.NumNodes
        m = sprintf ("%g", this.NodeMean(ii));
        if (this.Children(ii,1) == 0)
          fprintf ("%*d  fit = %s\n", w, ii, m);
        else
          v = this.PredictorNames{this.CutPredictorIndex(ii)};
          c = sprintf ("%g", this.CutPoint(ii));
          fprintf ("%*d  if %s<%s then node %d elseif %s>=%s then node", ...
                   w, ii, v, c, this.Children(ii,1), v, c);
          fprintf (" %d else %s\n", this.Children(ii,2), m);
        endif
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionTree} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {RegressionTree} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss on new data.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## weighted mean squared error of the response the model predicts for
    ## @var{X} against the observed response @var{Y}.  A row whose response
    ## is missing is dropped, as it is when fitting.
    ##
    ## @code{@var{L} = loss (@dots{}, @var{name}, @var{value})} takes the
    ## following options.
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab @qcode{'mse'}, the default, or a function
    ## handle taking the true response, the predicted response and the
    ## weights, and returning a numeric scalar.
    ##
    ## @item @qcode{'Weights'} @tab A numeric vector of observation weights,
    ## one per row of @var{X}, normalized to sum to one before it is applied.
    ##
    ## @end multitable
    ##
    ## @seealso{RegressionTree, fitrtree, RegressionTree.predict}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Input validation
      if (nargin < 3)
        error ("RegressionTree.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionTree.loss: name-value arguments must be", ...
                       " in pairs."));
      endif
      if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
        error ("RegressionTree.loss: X must be a real numeric matrix.");
      endif
      if (! (isnumeric (Y) && isreal (Y) && isvector (Y)))
        error ("RegressionTree.loss: Y must be a real numeric vector.");
      endif
      if (rows (X) != numel (Y))
        error (strcat ("RegressionTree.loss: number of rows in X and Y", ...
                       " must be equal."));
      endif

      LossFun = 'mse';
      Weights = [];

      while (numel (varargin) > 0)
        Value = varargin{2};
        switch (tolower (varargin{1}))
          case 'lossfun'
            if (! (is_function_handle (Value)
                   || (ischar (Value) && isrow (Value))))
              error (strcat ("RegressionTree.loss: 'LossFun' must be a", ...
                             " character vector or a function handle."));
            endif
            if (ischar (Value) && ! strcmpi (Value, 'mse'))
              error ("RegressionTree.loss: unsupported 'LossFun' value.");
            endif
            LossFun = Value;
          case 'weights'
            if (! (isnumeric (Value) && isvector (Value) && isreal (Value)))
              error (strcat ("RegressionTree.loss: 'Weights' must be a", ...
                             " real numeric vector."));
            endif
            if (numel (Value) != rows (X))
              error (strcat ("RegressionTree.loss: 'Weights' must have one", ...
                             " element per observation."));
            endif
            Weights = Value;
          otherwise
            error (strcat ("RegressionTree.loss: invalid parameter name in", ...
                           " optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      if (isempty (Weights))
        W = ones (rows (X), 1);
      else
        W = double (Weights(:));
      endif

      ## An observation with no response is one the loss cannot be measured
      ## on, so it is dropped and the weights are normalized over what is
      ## left, which is what the fit itself does with such a row.  Measured
      ## on R2024a, where the loss of a carsmall tree over the whole hundred
      ## rows is the loss over the ninety-four with a response.
      Y = double (Y(:));
      keep = ! isnan (Y);
      Y = Y(keep);
      X = X(keep, :);
      W = W(keep);
      if (! (sum (W) > 0))
        error (strcat ("RegressionTree.loss: 'Weights' must not be zero", ...
                       " for every observation with a response."));
      endif
      ## Weights are normalized to sum to one, as MATLAB does, so a loss is
      ## a weighted average rather than a weighted sum.
      W = W / sum (W);

      yFit = predict (this, X);
      if (is_function_handle (LossFun))
        L = LossFun (Y, yFit, W);
        if (! (isnumeric (L) && isscalar (L)))
          error (strcat ("RegressionTree.loss: 'LossFun' must return a", ...
                         " numeric scalar."));
        endif
      else
        L = sum (W .* (Y - yFit) .^ 2);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionTree} {@var{yFit} =} resubPredict (@var{obj})
    ##
    ## Predict the training response with the model fitted to it.
    ##
    ## @code{@var{yFit} = resubPredict (@var{obj})} is
    ## @code{predict (@var{obj}, @var{obj}.X)}.
    ##
    ## @seealso{RegressionTree, RegressionTree.predict}
    ## @end deftypefn
    function yFit = resubPredict (this)

      yFit = predict (this, this.X);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionTree} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {RegressionTree} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Regression loss on the training data.
    ##
    ## @code{@var{L} = resubLoss (@var{obj})} is @code{loss} over the
    ## training data, weighed as the fit weighed it, and takes the same
    ## @qcode{'LossFun'} option.  Giving @qcode{'Weights'} weighs the
    ## training data some other way instead.
    ##
    ## @seealso{RegressionTree, RegressionTree.loss}
    ## @end deftypefn
    function L = resubLoss (this, varargin)

      ## The training data is weighed as the fit weighed it.  Measured on
      ## R2024a: a weighted fit reports a resubstitution loss over its own
      ## weights, 6.51484752783409 on a carsmall tree weighted 1 to 100,
      ## where uniform weights give 7.44189886464653.
      if (! any (strcmpi (varargin(1:2:end), 'weights')))
        varargin = [varargin, {'Weights', this.RawWeights}];
      endif
      L = loss (this, this.X, this.Y, varargin{:});

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {RegressionTree} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a RegressionTree model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## RegressionTree object into an Octave binary file, the name of which is
    ## specified in @var{filename}, along with an extra variable, which
    ## defines the type of regression object these variables constitute.  Use
    ## @code{loadmodel} in order to load a regression object into Octave's
    ## workspace.
    ##
    ## @seealso{loadmodel, fitrtree, RegressionTree}
    ## @end deftypefn
    function savemodel (this, fname)

      ## Input validation
      if (nargin < 2)
        error ("RegressionTree.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("RegressionTree.savemodel: FNAME must be a", ...
                       " character vector."));
      endif

      ## Generate variable for class name
      classdef_name = 'RegressionTree';

      ## Create variables from model properties
      X                      = this.X;
      Y                      = this.Y;
      RowsUsed               = this.RowsUsed;
      W                      = this.W;
      RawWeights             = this.RawWeights;
      NumObservations        = this.NumObservations;
      PredictorNames         = this.PredictorNames;
      ResponseName           = this.ResponseName;
      CategoricalPredictors  = this.CategoricalPredictors;
      ExpandedPredictorNames = this.ExpandedPredictorNames;
      BinEdges               = this.BinEdges;
      ModelParameters        = this.ModelParameters;
      NumNodes               = this.NumNodes;
      Children               = this.Children;
      Parent                 = this.Parent;
      CutPredictorIndex      = this.CutPredictorIndex;
      CutPoint               = this.CutPoint;
      NodeSize               = this.NodeSize;
      NodeWeight             = this.NodeWeight;
      NodeMean               = this.NodeMean;
      NodeError              = this.NodeError;
      PruneList              = this.PruneList;
      PruneAlpha             = this.PruneAlpha;
      ResponseTransform      = this.ResponseTransform;
      RTfun                  = this.RTfun;
      HyperparameterOptimizationResults = ...
                               this.HyperparameterOptimizationResults;

      ## The cut descriptions, the two node probabilities and the held risk
      ## are not saved: every one of them is re-derived from the node table
      ## above, so writing them out would only make a stale copy possible.
      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'RowsUsed', 'W', ...
            'RawWeights', 'NumObservations', 'PredictorNames', ...
            'ResponseName', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'BinEdges', 'ModelParameters', ...
            'NumNodes', 'Children', 'Parent', 'CutPredictorIndex', ...
            'CutPoint', 'NodeSize', 'NodeWeight', 'NodeMean', 'NodeError', ...
            'PruneList', 'PruneAlpha', 'ResponseTransform', 'RTfun', ...
            'HyperparameterOptimizationResults');

    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      ## The smallest fit the class accepts, filled property by property
      ## below.  Nothing of the stub survives the copy.
      mdl = RegressionTree ([1; 2], [1; 2]);

      ## Copy the saved data into the object.  Iterate over what was saved
      ## rather than over fieldnames (mdl): a hidden property such as RTfun
      ## is written out by savemodel but is not reported by fieldnames, so
      ## comparing the two sets could never match.  Assignment is legal here
      ## because this is a method of the class itself.
      names = fieldnames (data);
      order = {'ResponseTransform'};
      late = ismember (names, order);
      tail = order(ismember (order, names));
      names = [names(! late); tail(:)];
      for ii = 1:numel (names)
        try
          mdl.(names{ii}) = data.(names{ii});
        catch
          msg = 'RegressionTree.load_model: invalid model in ''%s''.';
          error (msg, filename);
        end_try_catch
      endfor

      ## The cut descriptions and the derived node quantities were not saved
      mdl = fillCuts (mdl);
      mdl = deriveNodes (mdl);
      mdl = deriveHeld (mdl);

    endfunction

  endmethods

  methods (Access = private)

    ## The two node quantities that are derived from the weights the engine
    ## accumulated rather than reported by it.
    function this = deriveNodes (this)

      if (isempty (this.NodeWeight))
        return;
      endif
      if (this.NodeWeight(1) > 0)
        this.NodeProbability = this.NodeWeight / this.NodeWeight(1);
      else
        this.NodeProbability = zeros (size (this.NodeWeight));
      endif
      ## Measured on R2024a, where NodeRisk and NodeProbability times
      ## NodeError agree to every digit on every node of a carsmall tree.
      this.NodeRisk = this.NodeProbability .* this.NodeError;

    endfunction

    ## The descriptions of the cuts, derived from the node table.
    function this = fillCuts (this)

      S = treeCutInfo (this.CutPredictorIndex, this.PredictorNames);
      for [val, name] = S
        this.(name) = val;
      endfor

    endfunction

    ## The risk each branch node carries on account of the observations that
    ## stop there, missing the predictor it cuts on.  Those observations are
    ## in neither child, and MATLAB charges them the node's own mean squared
    ## error rather than their own: measured on R2024a, where one row of
    ## carsmall missing its horsepower gives a node of 58 observations two
    ## children of 40 and 17, and that node's alpha comes to
    ## 5.99325416896717, which is its risk less its children's less one
    ## observation's weight times its own NodeError.
    function this = deriveHeld (this)

      this.HeldRisk = zeros (this.NumNodes, 1);
      br = find (this.Children(:,1) > 0);
      if (isempty (br))
        return;
      endif
      kids = this.Children(br,:);
      wheld = this.NodeProbability(br) - this.NodeProbability(kids(:,1)) ...
              - this.NodeProbability(kids(:,2));
      this.HeldRisk(br) = wheld .* this.NodeError(br);

    endfunction

    ## The branch nodes a pruning level turns into leaves.
    function nodes = nodesAtLevel (this, L, caller)

      if (isempty (this.PruneList))
        error (strcat ("RegressionTree.%s: the tree carries no pruning", ...
                       " sequence."), caller);
      endif
      maxL = max (this.PruneList);
      if (L > maxL)
        warning (strcat ("RegressionTree.%s: pruning level %d is greater", ...
                         " than the largest level %d; the tree will be", ...
                         " pruned to its root."), caller, L, maxL);
        L = maxL;
      endif
      nodes = find (this.PruneList >= 1 & this.PruneList <= L);

    endfunction

    ## Turn the named branch nodes into leaves, discard everything below
    ## them, and renumber what is left.
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
      this.NodeWeight = this.NodeWeight(S.keep);
      this.NodeMean = this.NodeMean(S.keep);
      this.NodeError = this.NodeError(S.keep);
      this.PruneList = zeros (this.NumNodes, 1);
      this.PruneAlpha = [];

      this = fillCuts (this);
      this = deriveNodes (this);
      ## A node that became a leaf now holds back nothing, so this is
      ## re-derived rather than subset.
      this = deriveHeld (this);

    endfunction

    ## The cost complexity pruning sequence, on the squared error each node
    ## would carry as a leaf.
    function this = pruneSequence (this)

      risk = this.NodeProbability .* this.NodeError;
      [this.PruneList, this.PruneAlpha] = ...
        treePruneSequence (this.Children, this.Parent, risk, this.HeldRisk);

    endfunction

  endmethods

endclassdef

## Tests
%!test  # MATLAB parity: the surface a default fit reports
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! assert_equal (class (Mdl), 'RegressionTree');
%! assert_equal (numel (properties (Mdl)), 36);
%! assert_equal (Mdl.NumObservations, 94);
%! assert_equal (Mdl.ResponseName, 'Y');
%! assert_equal (Mdl.PredictorNames, {'x1', 'x2', 'x3'});
%! assert_equal (Mdl.ExpandedPredictorNames, {'x1', 'x2', 'x3'});
%! assert_equal (Mdl.ResponseTransform, 'none');
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (sum (Mdl.W), 1, 1e-14);
%! assert_equal (Mdl.W(1), 1 / 94, 1e-15);

%!test  # MATLAB parity: the parameters a default fit reports
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! MP = RegressionTree (X, MPG).ModelParameters;
%! assert_equal (MP.SplitCriterion, 'mse');
%! assert_equal (MP.PruneCriterion, 'mse');
%! assert_equal ({MP.MinParent, MP.MinLeaf, MP.MaxSplits}, {10, 1, 93});
%! assert_equal ({MP.MergeLeaves, MP.Prune}, {'on', 'on'});
%! assert_equal (MP.QEToler, 1e-6);
%! assert_equal ({MP.NSurrogate, MP.MaxCat, MP.AlgCat}, {0, 10, 'auto'});
%! assert_equal ({MP.Method, MP.Type}, {'Tree', 'regression'});

%!test  # MATLAB parity: what a tree with no categories and no surrogates holds
%! load carsmall
%! Mdl = RegressionTree ([Weight, Cylinders, Horsepower], MPG);
%! assert_equal (size (Mdl.CutCategories), [37, 2]);
%! assert_equal (size (Mdl.CategoricalSplit), [0, 0]);
%! assert_equal (size (Mdl.SurrogateCutPredictor), [0, 1]);
%! assert_equal (size (Mdl.SurrogateCutPoint), [0, 0]);
%! assert_equal (Mdl.BinEdges, {});
%! assert_equal (Mdl.HyperparameterOptimizationResults, []);

%!test  # A missing response drops its row, a missing predictor does not
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! assert_equal (Mdl.NumObservations, 94);
%! assert_equal (numel (Mdl.RowsUsed), 100);
%! assert_equal (sum (Mdl.RowsUsed), 94);
%! assert_equal (Mdl.NodeSize(1), 94);
%! ## the row with no horsepower is kept, and node 2 cuts on horsepower
%! assert_equal (Mdl.NodeSize(2) - sum (Mdl.NodeSize(Mdl.Children(2,:))), 1);

%!test  # MATLAB parity: a row missing the split predictor stops at that node
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! [yFit, node] = predict (Mdl, [NaN, NaN, NaN]);
%! assert_equal (node, 1);
%! assert_equal (yFit, Mdl.NodeMean(1), 1e-15);
%! [yFit, node] = predict (Mdl, [2000, 4, NaN]);
%! assert_equal (node, 2);
%! assert_equal (yFit, Mdl.NodeMean(2), 1e-15);

%!test  # MATLAB parity: prune takes a subtree out of the sequence
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG, 'MinLeafSize', 15);
%! sub = prune (Mdl, 'Level', 1);
%! assert_equal (sub.NumNodes, 7);
%! assert_equal (sub.NodeSize', [94, 58, 36, 40, 17, 18, 22]);
%! assert_equal (sub.PruneList', [3, 2, 0, 1, 0, 0, 0]);
%! assert_equal (sub.PruneAlpha', [0, 1.95238622931442, ...
%!                                 5.99325416896717, 41.4954735525515], 1e-11);

%!test  # MATLAB parity: a level of two, and level zero changing nothing
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG, 'MinLeafSize', 15);
%! assert_equal (prune (Mdl, 'Level', 2).NumNodes, 5);
%! assert_equal (prune (Mdl, 'Level', 0).NumNodes, 9);
%! assert_equal (prune (Mdl).NumNodes, 9);
%! assert_equal (prune (Mdl, 'Alpha', 100).NumNodes, 1);
%! assert_equal (prune (Mdl, 'Alpha', 0).NumNodes, 9);
%! assert_equal (prune (Mdl, 'Nodes', 2).NumNodes, 5);

%!warning<pruning level 9 is greater than the largest level 4>
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! sub = prune (RegressionTree (X, MPG, 'MinLeafSize', 15), 'Level', 9);

%!test  # MATLAB parity: the range of each predictor on the path to a node
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! assert_equal (fieldnames (nodeVariableRange (Mdl, 1)), cell (0, 1));
%! r = nodeVariableRange (Mdl, 4);
%! assert_equal (r.x1, [-Inf, 3085.5], 1e-12);
%! assert_equal (r.x3, [-Inf, 89], 1e-12);

%!test  # MATLAB parity: the text form of the tree
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG, 'MinLeafSize', 15);
%! lines = strsplit (strtrim (evalc ('view (Mdl)')), "\n");
%! assert_equal (numel (lines), 10);
%! assert_equal (lines{1}, 'Decision tree for regression');
%! assert_equal (lines{2}, ...
%!   '1  if x1<3085.5 then node 2 elseif x1>=3085.5 then node 3 else 23.7181');
%! assert_equal (lines{6}, '5  fit = 24.0882');

%!test  # The node numbers are right aligned on the widest of them
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! lines = strsplit (strtrim (evalc ('view (RegressionTree (X, MPG))')), "\n");
%! assert_equal (numel (lines), 38);
%! assert_equal (lines{2}(1:2), ' 1');
%! assert_equal (lines{11}(1:2), '10');

%!test  # The loss takes a function handle over the true, fitted and weights
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! r = Mdl.Y - resubPredict (Mdl);
%! mae = @(y, yf, w) sum (w .* abs (y - yf));
%! assert_equal (loss (Mdl, Mdl.X, Mdl.Y, 'LossFun', mae), ...
%!               mean (abs (r)), 1e-12);
%! assert_equal (loss (Mdl, Mdl.X, Mdl.Y, 'LossFun', 'mse'), ...
%!               resubLoss (Mdl), 1e-14);

%!test  # resubPredict answers exactly as predict on the training data
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! assert_equal (resubPredict (Mdl), predict (Mdl, Mdl.X));

%!test  # The response transform is parsed by name as well as by handle
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG, 'ResponseTransform', 'exp');
%! assert_equal (Mdl.ResponseTransform, 'exp');
%! assert_equal (predict (Mdl, X(1, :)), exp (17.25), 1e-10);
%! Mdl.ResponseTransform = 'none';
%! assert_equal (predict (Mdl, X(1, :)), 17.25, 1e-12);

%!test  # A tree of one node answers its mean everywhere
%! Mdl = RegressionTree ((1:5)', [2; 2; 2; 2; 2]);
%! assert_equal (Mdl.NumNodes, 1);
%! assert_equal (Mdl.NodeMean, 2);
%! assert_equal (Mdl.NodeError, 0);
%! assert_equal (Mdl.NodeRisk, 0);
%! assert_equal (predict (Mdl, 99), 2);
%! assert_equal (predictorImportance (Mdl), 0);

%!test  # A model saved and loaded answers exactly as it did
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! fname = tempname ();
%! unwind_protect
%!   savemodel (Mdl, fname);
%!   New = loadmodel (fname);
%!   assert_equal (class (New), 'RegressionTree');
%!   assert_equal (New.NumNodes, Mdl.NumNodes);
%!   assert_equal (New.NodeRisk, Mdl.NodeRisk, 1e-15);
%!   assert_equal (New.PruneAlpha, Mdl.PruneAlpha, 1e-12);
%!   assert_equal (New.CutPredictor, Mdl.CutPredictor);
%!   assert_equal (predict (New, X), predict (Mdl, X));
%! unwind_protect_cleanup
%!   delete (fname);
%! end_unwind_protect

%!test  # MATLAB parity: a small fixture, every node and every alpha
%! x = (1:12)';
%! y = [1; 1; 1; 1; 5; 5; 5; 5; 20; 20; 20; 20];
%! Mdl = RegressionTree (x, y, 'MinParentSize', 2);
%! assert_equal (Mdl.NumNodes, 5);
%! assert_equal (Mdl.NodeSize', [12, 8, 4, 4, 4]);
%! assert_equal (Mdl.NodeMean', [8.66666666666666, 3, 20, 1, 5], 1e-13);
%! assert_equal (Mdl.NodeError', [66.8888888888889, 4, 0, 0, 0], 1e-12);
%! assert_equal (Mdl.NodeRisk', [66.8888888888889, 2.66666666666667, ...
%!                               0, 0, 0], 1e-12);
%! assert_equal (Mdl.PruneList', [2, 1, 0, 0, 0]);
%! assert_equal (Mdl.PruneAlpha', [0, 2.66666666666667, ...
%!                                 64.2222222222222], 1e-12);
%! assert_equal (predictorImportance (Mdl), 33.4444444444444, 1e-12);

%!test  # MATLAB parity: a fixture grown to one observation per leaf
%! x = (1:8)';
%! y = [1; 2; 10; 11; 30; 31; 60; 61];
%! Mdl = RegressionTree (x, y, 'MinParentSize', 2);
%! assert_equal (Mdl.NumNodes, 15);
%! assert_equal (Mdl.NodeSize', [8, 6, 2, 4, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1]);
%! assert_equal (Mdl.NodeMean(1:5)', [25.75, 14.1666666666667, 60.5, ...
%!                                    6, 30.5], 1e-13);
%! assert_equal (Mdl.NodeRisk(1:5)', [512.9375, 110.354166666667, 0.0625, ...
%!                                    10.25, 0.0625], 1e-12);
%! assert_equal (Mdl.PruneList', [4, 3, 1, 2, 1, 0, 0, 1, 1, 0, 0, 0, 0, ...
%!                                0, 0]);
%! assert_equal (Mdl.PruneAlpha', [0, 0.0625, 10.125, 100.041666666667, ...
%!                                 402.520833333333], 1e-11);
%! assert_equal (predictorImportance (Mdl), 73.2767857142857, 1e-12);

%!test  # MATLAB parity: compact drops the data and answers identically
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), 'CompactRegressionTree');
%! assert_equal (numel (properties (CMdl)), 28);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));
%! assert_equal (CMdl.NodeRisk, Mdl.NodeRisk, 1e-15);
%! assert_equal (predictorImportance (CMdl), predictorImportance (Mdl), 1e-15);

%!test  # MATLAB parity: crossval returns a partitioned model over compacts
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! CVMdl = crossval (Mdl);
%! assert_equal (class (CVMdl), 'RegressionPartitionedModel');
%! assert_equal (CVMdl.CrossValidatedModel, 'Tree');
%! assert_equal (class (CVMdl.Trained{1}), 'CompactRegressionTree');
%! assert_equal (CVMdl.KFold, 10);
%! assert_equal (CVMdl.NumObservations, 94);
%! assert_equal (CVMdl.ResponseName, 'Y');
%! assert_equal (CVMdl.W, Mdl.W, 1e-15);

%!test  # Each way of asking for a partition gives the folds it names
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! assert_equal (crossval (Mdl, 'KFold', 5).KFold, 5);
%! assert_equal (crossval (Mdl, 'Holdout', 0.3).KFold, 1);
%! assert_equal (crossval (Mdl, 'Leaveout', 'on').KFold, 94);
%! assert_equal (crossval (Mdl, 'CVPartition', ...
%!                         cvpartition (94, 'KFold', 4)).KFold, 4);

%!test  # A fold is grown with the parameters the parent was grown with
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG, 'MinLeafSize', 7, 'MergeLeaves', 'off', ...
%!                       'QuadraticErrorTolerance', 0.01);
%! CVMdl = crossval (Mdl, 'KFold', 3);
%! assert_equal (CVMdl.ModelParameters.MinLeaf, 7);
%! assert_equal (CVMdl.ModelParameters.MergeLeaves, 'off');
%! assert_equal (CVMdl.ModelParameters.QEToler, 0.01);
%! assert_equal (CVMdl.Trained{1}.PredictorNames, Mdl.PredictorNames);

%!test  # kfoldPredict answers every observation out of fold
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! CVMdl = crossval (RegressionTree (X, MPG), 'KFold', 5);
%! yFit = kfoldPredict (CVMdl);
%! assert_equal (size (yFit), [94, 1]);
%! assert_equal (any (isnan (yFit)), false);
%! assert_equal (kfoldLoss (CVMdl) > 0, true);

%!test  # cvloss over a fixture whose folds all answer alike
%! ## The two groups are separated by a gap no fold can straddle, so every
%! ## fold grows the same tree and answers the held-out row exactly.  Left
%! ## with one leaf, a fold's tree answers the mean of the thirty-nine rows
%! ## it saw, which is 200/39 away from whichever value was held out,
%! ## whichever group that row came from.
%! x = [(1:20)'; (101:120)'];
%! y = [zeros(20, 1); 10 * ones(20, 1)];
%! Mdl = RegressionTree (x, y);
%! assert_equal (Mdl.NumNodes, 3);
%! [E, SE, Nleaf, BestLevel] = cvloss (Mdl, 'SubTrees', 'all', 'KFold', 40);
%! assert_equal (E', [0, (200 / 39) ^ 2], 1e-12);
%! assert_equal (SE', [0, 0], 1e-12);
%! assert_equal (Nleaf', [2, 1]);
%! assert_equal (BestLevel, 0);

%!test  # cvloss reports one row per subtree asked for
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG, 'MinLeafSize', 15);
%! [E, SE, Nleaf, BestLevel] = cvloss (Mdl, 'SubTrees', 'all');
%! assert_equal (size (E), [5, 1]);
%! assert_equal (size (SE), [5, 1]);
%! assert_equal (Nleaf', [5, 4, 3, 2, 1]);
%! assert_equal (all (E > 0), true);
%! assert_equal (BestLevel >= 0 && BestLevel <= 4, true);
%! ## the root alone answers the mean, so it can do no better than the whole
%! assert_equal (E(5) > E(1), true);

%!test  # cvloss defaults to the unpruned tree alone
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG, 'MinLeafSize', 15);
%! [E, SE, Nleaf, BestLevel] = cvloss (Mdl);
%! assert_equal (size (E), [1, 1]);
%! assert_equal (Nleaf, 5);
%! assert_equal (BestLevel, 0);
%! assert_equal (E > 0, true);

%!test  # A subset of the levels is answered in the order it was asked for
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG, 'MinLeafSize', 15);
%! [E, SE, Nleaf, BestLevel] = cvloss (Mdl, 'SubTrees', [0, 2, 4]);
%! assert_equal (Nleaf', [5, 3, 1]);
%! assert_equal (any (BestLevel == [0, 2, 4]), true);

%!test  # The fold count and the tree size rule are the ones asked for
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG, 'MinLeafSize', 15);
%! assert_equal (size (cvloss (Mdl, 'KFold', 5)), [1, 1]);
%! [~, ~, ~, bmin] = cvloss (Mdl, 'SubTrees', 'all', 'TreeSize', 'min');
%! [~, ~, ~, bse] = cvloss (Mdl, 'SubTrees', 'all', 'TreeSize', 'se');
%! assert_equal (bmin >= 0 && bmin <= 4, true);
%! assert_equal (bse >= 0 && bse <= 4, true);

%!test  # cvloss weighs the loss as the fit weighed it
%! ## A weighted fit whose folds all answer alike.  The heavier group pulls
%! ## the one-leaf answer towards itself, and the loss counts each row by
%! ## its weight rather than by its share of the rows.
%! x = [(1:20)'; (101:120)'];
%! y = [zeros(20, 1); 10 * ones(20, 1)];
%! Mdl = RegressionTree (x, y, 'Weights', [ones(20, 1); 3 * ones(20, 1)]);
%! E = cvloss (Mdl, 'SubTrees', 'all', 'KFold', 40);
%! assert_equal (E(1), 0, 1e-12);
%! assert_equal (E(2) > 0 && E(2) < 100, true);

%!test  # More folds than observations are reduced to one fold each
%! x = [(1:20)'; (101:120)'];
%! y = [zeros(20, 1); 10 * ones(20, 1)];
%! Mdl = RegressionTree (x, y);
%! ws = warning ('off', 'all');
%! unwind_protect
%!   E = cvloss (Mdl, 'KFold', 500, 'SubTrees', 'all');
%! unwind_protect_cleanup
%!   warning (ws);
%! end_unwind_protect
%! assert_equal (E', [0, (200 / 39) ^ 2], 1e-12);

## Test input validation
%!error<RegressionTree: too few input arguments.> RegressionTree ()
%!error<RegressionTree: too few input arguments.> RegressionTree (ones (4, 2))
%!error<RegressionTree: name-value arguments must be in pairs.>
%! RegressionTree (ones (4, 2), (1:4)', 'K')
%!error<RegressionTree: X must be a non-empty real numeric matrix.>
%! RegressionTree ('a', (1:4)')
%!error<RegressionTree: Y must be a real numeric vector.>
%! RegressionTree (ones (4, 2), {'a'; 'b'; 'c'; 'd'})
%!error<RegressionTree: number of rows in X and Y must be equal.>
%! RegressionTree (ones (4, 2), (1:3)')
%!error<RegressionTree: 'PredictorNames' must be supplied as a cellstring array.>
%! RegressionTree (ones (4, 2), (1:4)', 'PredictorNames', 'a')
%!error<RegressionTree: 'PredictorNames' must equal the number of columns in X.>
%! RegressionTree (ones (4, 2), (1:4)', 'PredictorNames', {'a'})
%!error<RegressionTree: 'ResponseName' must be a character vector.>
%! RegressionTree (ones (4, 2), (1:4)', 'ResponseName', 5)
%!error<RegressionTree: 'ResponseTransform' must be a character vector or a function handle.>
%! RegressionTree (ones (4, 2), (1:4)', 'ResponseTransform', 5)
%!error<RegressionTree: unrecognized 'ResponseTransform' function.>
%! RegressionTree (ones (4, 2), (1:4)', 'ResponseTransform', 'bogus')
%!error<RegressionTree: 'Weights' must be a real numeric vector.>
%! RegressionTree (ones (4, 2), (1:4)', 'Weights', 'a')
%!error<RegressionTree: 'Weights' must have one element per row in X.>
%! RegressionTree (ones (4, 2), (1:4)', 'Weights', [1, 2, 3])
%!error<RegressionTree: 'Weights' must be nonnegative and must not be all zero.>
%! RegressionTree (ones (4, 2), (1:4)', 'Weights', -ones (4, 1))
%!error<RegressionTree: 'MaxNumSplits' must be a nonnegative integer.>
%! RegressionTree (ones (4, 2), (1:4)', 'MaxNumSplits', -1)
%!error<RegressionTree: 'MinLeafSize' must be a positive integer.>
%! RegressionTree (ones (4, 2), (1:4)', 'MinLeafSize', 0)
%!error<RegressionTree: 'MinParentSize' must be a positive integer.>
%! RegressionTree (ones (4, 2), (1:4)', 'MinParentSize', 0)
%!error<RegressionTree: 'MergeLeaves' must be either 'on' or 'off'.>
%! RegressionTree (ones (4, 2), (1:4)', 'MergeLeaves', 'x')
%!error<RegressionTree: 'Prune' must be either 'on' or 'off'.>
%! RegressionTree (ones (4, 2), (1:4)', 'Prune', 'x')
%!error<RegressionTree: 'QuadraticErrorTolerance' must be a positive scalar.>
%! RegressionTree (ones (4, 2), (1:4)', 'QuadraticErrorTolerance', -1)
%!error<RegressionTree: 'SplitCriterion' must be 'mse' for a regression tree.>
%! RegressionTree (ones (4, 2), (1:4)', 'SplitCriterion', 'gdi')
%!error<RegressionTree: 'PruneCriterion' must be 'mse' for a regression tree.>
%! RegressionTree (ones (4, 2), (1:4)', 'PruneCriterion', 'error')
%!error<RegressionTree: categorical predictors are not implemented.>
%! RegressionTree (ones (4, 2), (1:4)', 'CategoricalPredictors', 1)
%!error<RegressionTree: 'Surrogate' is not implemented.>
%! RegressionTree (ones (4, 2), (1:4)', 'Surrogate', 'on')
%!error<RegressionTree: 'KFold' is not implemented; fit the model and cross-validate it afterwards.>
%! RegressionTree (ones (4, 2), (1:4)', 'KFold', 5)
%!error<RegressionTree: invalid parameter name in optional pair arguments.>
%! RegressionTree (ones (4, 2), (1:4)', 'Bogus', 1)
%!error<RegressionTree: no observations with a known response.>
%! RegressionTree (ones (4, 2), nan (4, 1))
%!error<RegressionTree.predict: too few input arguments.>
%! predict (RegressionTree (ones (4, 2), (1:4)'))
%!error<RegressionTree.predict: XC is empty.>
%! predict (RegressionTree (ones (4, 2), (1:4)'), [])
%!error<RegressionTree.predict: XC must be a real numeric matrix.>
%! predict (RegressionTree (ones (4, 2), (1:4)'), 'a')
%!error<RegressionTree.predict: XC must have the same number of predictors as the trained model.>
%! predict (RegressionTree (ones (4, 2), (1:4)'), ones (2, 5))
%!error<RegressionTree.loss: too few input arguments.>
%! loss (RegressionTree (ones (4, 2), (1:4)'), ones (4, 2))
%!error<RegressionTree.loss: name-value arguments must be in pairs.>
%! loss (RegressionTree (ones (4, 2), (1:4)'), ones (4, 2), (1:4)', 'LossFun')
%!error<RegressionTree.loss: 'LossFun' must be a character vector or a function handle.>
%! loss (RegressionTree (ones (4, 2), (1:4)'), ones (4, 2), (1:4)', ...
%!       'LossFun', 5)
%!error<RegressionTree.loss: unsupported 'LossFun' value.>
%! loss (RegressionTree (ones (4, 2), (1:4)'), ones (4, 2), (1:4)', ...
%!       'LossFun', 'mad')
%!error<RegressionTree.loss: 'LossFun' must return a numeric scalar.>
%! loss (RegressionTree (ones (4, 2), (1:4)'), ones (4, 2), (1:4)', ...
%!       'LossFun', @(y, f, w) [1, 2])
%!error<RegressionTree.loss: 'Weights' must be a real numeric vector.>
%! loss (RegressionTree (ones (4, 2), (1:4)'), ones (4, 2), (1:4)', ...
%!       'Weights', 'a')
%!error<RegressionTree.loss: 'Weights' must have one element per observation.>
%! loss (RegressionTree (ones (4, 2), (1:4)'), ones (4, 2), (1:4)', ...
%!       'Weights', [1, 2])
%!error<RegressionTree.loss: invalid parameter name in optional pair arguments.>
%! loss (RegressionTree (ones (4, 2), (1:4)'), ones (4, 2), (1:4)', 'Bogus', 1)
%!error<RegressionTree.loss: number of rows in X and Y must be equal.>
%! loss (RegressionTree (ones (4, 2), (1:4)'), ones (4, 2), (1:3)')
%!error<RegressionTree.prune: name-value arguments must be in pairs.>
%! prune (RegressionTree (ones (4, 2), (1:4)'), 'Level')
%!error<RegressionTree.prune: specify only one of the optional name-value paired arguments.>
%! prune (RegressionTree (ones (4, 2), (1:4)'), 'Level', 1, 'Alpha', 1)
%!error<RegressionTree.prune: 'Level' must be a nonnegative integer.>
%! prune (RegressionTree (ones (4, 2), (1:4)'), 'Level', -1)
%!error<RegressionTree.prune: 'Alpha' must be a nonnegative scalar.>
%! prune (RegressionTree (ones (4, 2), (1:4)'), 'Alpha', -1)
%!error<RegressionTree.prune: 'Nodes' must hold indices of nodes of the tree.>
%! prune (RegressionTree (ones (4, 2), (1:4)'), 'Nodes', 999)
%!error<RegressionTree.prune: invalid parameter name in optional pair arguments.>
%! prune (RegressionTree (ones (4, 2), (1:4)'), 'Bogus', 1)
%!error<RegressionTree.nodeVariableRange: too few input arguments.>
%! nodeVariableRange (RegressionTree (ones (4, 2), (1:4)'))
%!error<RegressionTree.nodeVariableRange: NODE must be a positive integer no greater than the number of nodes in the tree.>
%! nodeVariableRange (RegressionTree (ones (4, 2), (1:4)'), 999)
%!error<RegressionTree.cvloss: name-value arguments must be in pairs.>
%! cvloss (RegressionTree (ones (4, 2), (1:4)'), 'SubTrees')
%!error<RegressionTree.cvloss: 'SubTrees' must be 'all' or a vector of nonnegative integers in ascending order.>
%! load carsmall
%! cvloss (RegressionTree ([Weight, Cylinders], MPG), 'SubTrees', -1)
%!error<RegressionTree.cvloss: 'SubTrees' must be 'all' or a vector of nonnegative integers in ascending order.>
%! load carsmall
%! cvloss (RegressionTree ([Weight, Cylinders], MPG), 'SubTrees', [2, 1])
%!error<RegressionTree.cvloss: 'SubTrees' must not exceed the largest pruning level, 4.>
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! cvloss (RegressionTree (X, MPG, 'MinLeafSize', 15), 'SubTrees', 99)
%!error<RegressionTree.cvloss: 'TreeSize' must be either 'se' or 'min'.>
%! load carsmall
%! cvloss (RegressionTree ([Weight, Cylinders], MPG), 'TreeSize', 'x')
%!error<RegressionTree.cvloss: 'KFold' must be an integer value greater than 1.>
%! load carsmall
%! cvloss (RegressionTree ([Weight, Cylinders], MPG), 'KFold', 1)
%!error<RegressionTree.cvloss: invalid parameter name in optional pair arguments.>
%! load carsmall
%! cvloss (RegressionTree ([Weight, Cylinders], MPG), 'Bogus', 1)
%!error<RegressionTree.cvloss: the tree carries no pruning sequence; fit it with 'Prune' or 'MergeLeaves' on.>
%! load carsmall
%! cvloss (RegressionTree ([Weight, Cylinders], MPG, 'Prune', 'off', ...
%!                         'MergeLeaves', 'off'))
%!error<RegressionTree.crossval: Name-Value arguments must be in pairs.>
%! crossval (RegressionTree (ones (4, 2), (1:4)'), 'KFold')
%!error<RegressionTree.crossval: specify only one of the optional Name-Value paired arguments.>
%! crossval (RegressionTree (ones (4, 2), (1:4)'), 'KFold', 2, 'Holdout', 0.3)
%!error<RegressionTree.crossval: 'KFold' must be an integer value greater than 1.>
%! crossval (RegressionTree (ones (4, 2), (1:4)'), 'KFold', 1)
%!error<RegressionTree.crossval: 'Holdout' must be a numeric value between 0 and 1.>
%! crossval (RegressionTree (ones (4, 2), (1:4)'), 'Holdout', 2)
%!error<RegressionTree.crossval: 'Leaveout' must be either 'on' or 'off'.>
%! crossval (RegressionTree (ones (4, 2), (1:4)'), 'Leaveout', 'x')
%!error<RegressionTree.crossval: 'CVPartition' must be a 'cvpartition' object.>
%! crossval (RegressionTree (ones (4, 2), (1:4)'), 'CVPartition', 5)
%!error<RegressionTree.crossval: invalid parameter name in optional paired arguments.>
%! crossval (RegressionTree (ones (4, 2), (1:4)'), 'Bogus', 1)
%!error<RegressionTree.savemodel: too few input arguments.>
%! savemodel (RegressionTree (ones (4, 2), (1:4)'))
%!error<RegressionTree.savemodel: FNAME must be a character vector.>
%! savemodel (RegressionTree (ones (4, 2), (1:4)'), 5)
