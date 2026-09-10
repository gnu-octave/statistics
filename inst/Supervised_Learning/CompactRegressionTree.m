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

classdef CompactRegressionTree
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactRegressionTree
  ##
  ## Compact binary decision tree for regression
  ##
  ## A @code{CompactRegressionTree} object carries the tree a
  ## @code{RegressionTree} model grew and everything @code{predict} needs, but
  ## not the observations it was fitted on.  It answers new data identically
  ## to the model it came from, and is far smaller to keep or to ship.
  ##
  ## Create one with the @code{compact} method of a @code{RegressionTree}
  ## object.  Because it holds no training data, it has no @code{resub}
  ## methods and cannot be cross-validated, and it cannot be pruned: the
  ## pruning sequence is reported but taking a subtree out of it rewrites the
  ## node table, which is work for the model that still has its data.
  ##
  ## @seealso{RegressionTree, fitrtree}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} NumNodes
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
    ## @deftp {CompactRegressionTree} {property} Children
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
    ## @deftp {CompactRegressionTree} {property} Parent
    ##
    ## Parent of each node
    ##
    ## A column vector naming the parent of each node.  The root carries a
    ## zero.  This property is read-only.
    ##
    ## @end deftp
    Parent = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} IsBranchNode
    ##
    ## Which nodes are branch nodes
    ##
    ## A logical column vector, true for each node that carries a split and
    ## false for each leaf.  This property is read-only.
    ##
    ## @end deftp
    IsBranchNode = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} CutPredictor
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
    ## @deftp {CompactRegressionTree} {property} CutPredictorIndex
    ##
    ## Index of the predictor each node cuts on
    ##
    ## A column vector holding, for each node, the column of @var{X} the node
    ## splits on, and zero at a leaf.  This property is read-only.
    ##
    ## @end deftp
    CutPredictorIndex = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} CutPoint
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
    ## @deftp {CompactRegressionTree} {property} CutType
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
    ## @deftp {CompactRegressionTree} {property} CutCategories
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
    ## @deftp {CompactRegressionTree} {property} CategoricalSplit
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
    ## @deftp {CompactRegressionTree} {property} NodeSize
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
    ## @deftp {CompactRegressionTree} {property} NodeMean
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
    ## @deftp {CompactRegressionTree} {property} NodeError
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
    ## @deftp {CompactRegressionTree} {property} NodeProbability
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
    ## @deftp {CompactRegressionTree} {property} NodeRisk
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
    ## @deftp {CompactRegressionTree} {property} PruneList
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
    ## @deftp {CompactRegressionTree} {property} PruneAlpha
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
    ## @deftp {CompactRegressionTree} {property} SurrogateCutCategories
    ##
    ## Categories of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutCategories = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} SurrogateCutFlip
    ##
    ## Cut assignments of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutFlip = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} SurrogateCutPoint
    ##
    ## Cut points of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutPoint = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} SurrogateCutType
    ##
    ## Types of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutType = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} SurrogateCutPredictor
    ##
    ## Predictors of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogateCutPredictor = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} SurrogatePredictorAssociation
    ##
    ## Predictive measures of association of the surrogate splits
    ##
    ## Surrogate splits are not implemented, so this is always empty.  This
    ## property is read-only.
    ##
    ## @end deftp
    SurrogatePredictorAssociation = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} PredictorNames
    ##
    ## Names of the predictor variables
    ##
    ## A cell array of character vectors with one name per column of
    ## @var{X}.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} ResponseName
    ##
    ## Name of the response variable
    ##
    ## A character vector naming the response.  This property is read-only.
    ##
    ## @end deftp
    ResponseName = [];

    ## -*- texinfo -*-
    ## @deftp {CompactRegressionTree} {property} CategoricalPredictors
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
    ## @deftp {CompactRegressionTree} {property} ExpandedPredictorNames
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
    ## @deftp {CompactRegressionTree} {property} ResponseTransform
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

    ## The risk each branch node carries on account of the observations that
    ## stop there, missing the predictor it cuts on.  Those observations are
    ## in neither child, so what a split buys is the drop in risk less this.
    HeldRisk = [];

  endproperties

  ## Set methods for the properties a user may assign after compacting.
  methods (Hidden)

    function this = set.ResponseTransform (this, val)
      [this.RTfun, this.ResponseTransform] = ...
                     parseResponseTransform (val, 'CompactRegressionTree');
    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ('\n  CompactRegressionTree\n\n');
      fprintf ('%22s: %s\n', 'ResponseName', this.ResponseName);
      fprintf ('%22s: %s\n', 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ('%22s: %s\n', 'ResponseTransform', this.ResponseTransform);
      fprintf ('%22s: %d\n', 'NumNodes', this.NumNodes);
      fprintf ('\n');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactRegressionTree} {@var{obj} =} CompactRegressionTree (@var{Mdl})
    ##
    ## Create a @code{CompactRegressionTree} object.
    ##
    ## @var{Mdl} is the @code{RegressionTree} object to compact.  The
    ## documented way to reach this constructor is the @code{compact} method.
    ##
    ## @seealso{RegressionTree, fitrtree}
    ## @end deftypefn
    function this = CompactRegressionTree (Mdl)

      ## Input validation
      if (nargin < 1)
        error ("CompactRegressionTree: too few input arguments.");
      endif
      if (! isa (Mdl, 'RegressionTree'))
        error (strcat ("CompactRegressionTree: MDL must be a", ...
                       " RegressionTree object."));
      endif

      ## The node table and everything derived from it, then the transform.
      this.NumNodes = Mdl.NumNodes;
      this.Children = Mdl.Children;
      this.Parent = Mdl.Parent;
      this.CutPredictorIndex = Mdl.CutPredictorIndex;
      this.CutPoint = Mdl.CutPoint;
      this.NodeSize = Mdl.NodeSize;
      this.NodeMean = Mdl.NodeMean;
      this.NodeError = Mdl.NodeError;
      this.NodeProbability = Mdl.NodeProbability;
      this.NodeRisk = Mdl.NodeRisk;
      this.HeldRisk = Mdl.HeldRisk;
      this.PruneList = Mdl.PruneList;
      this.PruneAlpha = Mdl.PruneAlpha;

      this.PredictorNames = Mdl.PredictorNames;
      this.ResponseName = Mdl.ResponseName;
      this.CategoricalPredictors = Mdl.CategoricalPredictors;
      this.ExpandedPredictorNames = Mdl.ExpandedPredictorNames;

      this = fillCuts (this);
      this.ResponseTransform = Mdl.ResponseTransform;

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactRegressionTree} {@var{yFit} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactRegressionTree} {[@var{yFit}, @var{node}] =} predict (@dots{})
    ##
    ## Predict the response with a trained @code{CompactRegressionTree} object.
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
    ## @seealso{CompactRegressionTree, fitrtree}
    ## @end deftypefn
    function [yFit, node] = predict (this, XC)

      ## Input validation
      if (nargin < 2)
        error ("CompactRegressionTree.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("CompactRegressionTree.predict: XC is empty.");
      endif
      if (! (isnumeric (XC) && isreal (XC) && ismatrix (XC)))
        error (strcat ("CompactRegressionTree.predict: XC must be a", ...
                       " real numeric matrix."));
      endif
      if (numel (this.PredictorNames) != columns (XC))
        error (strcat ("CompactRegressionTree.predict: XC must have the", ...
                       " same number of predictors as the trained", ...
                       " model."));
      endif

      [yFit, node] = treepredict (XC, this.Children, ...
                                  this.CutPredictorIndex, this.CutPoint, ...
                                  this.NodeMean);
      yFit = this.RTfun (yFit);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactRegressionTree} {@var{imp} =} predictorImportance (@var{obj})
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
    ## @seealso{CompactRegressionTree, fitrtree, CompactRegressionTree.NodeRisk}
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
    ## @deftypefn {CompactRegressionTree} {@var{r} =} nodeVariableRange (@var{obj}, @var{node})
    ##
    ## Range of each predictor at a node.
    ##
    ## @code{@var{r} = nodeVariableRange (@var{obj}, @var{node})} returns a
    ## structure with one field per predictor the path from the root to
    ## @var{node} cuts on, holding the two-element range of values that reach
    ## the node.  A predictor the path never cuts on is unconstrained and is
    ## left out, so the root gives a structure with no fields.
    ##
    ## @seealso{CompactRegressionTree, fitrtree}
    ## @end deftypefn
    function r = nodeVariableRange (this, node)

      ## Input validation
      if (nargin < 2)
        error (strcat ("CompactRegressionTree.nodeVariableRange: too few", ...
                       " input arguments."));
      endif
      if (! (isnumeric (node) && isscalar (node) && isreal (node)
             && node >= 1 && node <= this.NumNodes && node == fix (node)))
        error (strcat ("CompactRegressionTree.nodeVariableRange: NODE must", ...
                       " be a positive integer no greater than the", ...
                       " number of nodes in the tree."));
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
    ## @deftypefn {CompactRegressionTree} {} view (@var{obj})
    ##
    ## Print the tree as text.
    ##
    ## @code{view (@var{obj})} prints one line per node: a branch node names
    ## the predictor it cuts on, the cut point, and the node each side leads
    ## to, and a leaf names the response it fits.  A branch node's line ends
    ## with the response it would fit itself, which is the answer an
    ## observation missing that predictor gets.
    ##
    ## @seealso{CompactRegressionTree, fitrtree}
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
    ## @deftypefn  {CompactRegressionTree} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactRegressionTree} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
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
    ## @seealso{CompactRegressionTree, fitrtree, CompactRegressionTree.predict}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Input validation
      if (nargin < 3)
        error ("CompactRegressionTree.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactRegressionTree.loss: name-value arguments", ...
                       " must be in pairs."));
      endif
      if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
        error ("CompactRegressionTree.loss: X must be a real numeric matrix.");
      endif
      if (! (isnumeric (Y) && isreal (Y) && isvector (Y)))
        error ("CompactRegressionTree.loss: Y must be a real numeric vector.");
      endif
      if (rows (X) != numel (Y))
        error (strcat ("CompactRegressionTree.loss: number of rows in X", ...
                       " and Y must be equal."));
      endif

      LossFun = 'mse';
      Weights = [];

      while (numel (varargin) > 0)
        Value = varargin{2};
        switch (tolower (varargin{1}))
          case 'lossfun'
            if (! (is_function_handle (Value)
                   || (ischar (Value) && isrow (Value))))
              error (strcat ("CompactRegressionTree.loss:", ...
                             " 'LossFun' must be a character vector or a", ...
                             " function handle."));
            endif
            if (ischar (Value) && ! strcmpi (Value, 'mse'))
              error (strcat ("CompactRegressionTree.loss:", ...
                             " unsupported 'LossFun' value."));
            endif
            LossFun = Value;
          case 'weights'
            if (! (isnumeric (Value) && isvector (Value) && isreal (Value)))
              error (strcat ("CompactRegressionTree.loss:", ...
                             " 'Weights' must be a real numeric vector."));
            endif
            if (numel (Value) != rows (X))
              error (strcat ("CompactRegressionTree.loss:", ...
                             " 'Weights' must have one element per", ...
                             " observation."));
            endif
            Weights = Value;
          otherwise
            error (strcat ("CompactRegressionTree.loss: invalid parameter", ...
                           " name in optional pair arguments."));
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
        error (strcat ("CompactRegressionTree.loss: 'Weights' must not be", ...
                       " zero for every observation with a response."));
      endif
      ## Weights are normalized to sum to one, as MATLAB does, so a loss is
      ## a weighted average rather than a weighted sum.
      W = W / sum (W);

      yFit = predict (this, X);
      if (is_function_handle (LossFun))
        L = LossFun (Y, yFit, W);
        if (! (isnumeric (L) && isscalar (L)))
          error (strcat ("CompactRegressionTree.loss:", ...
                         " 'LossFun' must return a numeric scalar."));
        endif
      else
        L = sum (W .* (Y - yFit) .^ 2);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactRegressionTree} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a CompactRegressionTree model to a file.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## CompactRegressionTree object into an Octave binary file, the name of
    ## which is specified in @var{filename}, along with an extra variable,
    ## which defines the type of regression object these variables
    ## constitute.  Use @code{loadmodel} in order to load a regression object
    ## into Octave's workspace.
    ##
    ## @seealso{loadmodel, fitrtree, RegressionTree}
    ## @end deftypefn
    function savemodel (this, fname)

      ## Input validation
      if (nargin < 2)
        error (strcat ("CompactRegressionTree.savemodel: too few input", ...
                       " arguments."));
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("CompactRegressionTree.savemodel: FNAME must be a", ...
                       " character vector."));
      endif

      ## Generate variable for class name
      classdef_name = 'CompactRegressionTree';

      ## Create variables from model properties
      PredictorNames         = this.PredictorNames;
      ResponseName           = this.ResponseName;
      CategoricalPredictors  = this.CategoricalPredictors;
      ExpandedPredictorNames = this.ExpandedPredictorNames;
      NumNodes               = this.NumNodes;
      Children               = this.Children;
      Parent                 = this.Parent;
      CutPredictorIndex      = this.CutPredictorIndex;
      CutPoint               = this.CutPoint;
      NodeSize               = this.NodeSize;
      NodeMean               = this.NodeMean;
      NodeError              = this.NodeError;
      NodeProbability        = this.NodeProbability;
      NodeRisk               = this.NodeRisk;
      HeldRisk               = this.HeldRisk;
      PruneList              = this.PruneList;
      PruneAlpha             = this.PruneAlpha;
      ResponseTransform      = this.ResponseTransform;
      RTfun                  = this.RTfun;

      ## The cut descriptions are not saved: every one of them is re-derived
      ## from the node table above, so writing them out would only make a
      ## stale copy possible.
      save ('-binary', fname, 'classdef_name', 'PredictorNames', ...
            'ResponseName', 'CategoricalPredictors', ...
            'ExpandedPredictorNames', 'NumNodes', 'Children', 'Parent', ...
            'CutPredictorIndex', 'CutPoint', 'NodeSize', 'NodeMean', ...
            'NodeError', 'NodeProbability', 'NodeRisk', 'HeldRisk', ...
            'PruneList', 'PruneAlpha', 'ResponseTransform', 'RTfun');

    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      ## The compact model is built from a full one and has no training data
      ## of its own, so the smallest fit the full class accepts is compacted
      ## and then filled property by property.
      mdl = CompactRegressionTree (RegressionTree ([1; 2], [1; 2]));

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
          msg = strcat ("CompactRegressionTree.load_model: invalid model", ...
                        " in '%s'.");
          error (msg, filename);
        end_try_catch
      endfor

      ## The cut descriptions were not saved, being derived
      mdl = fillCuts (mdl);

    endfunction

  endmethods

  methods (Access = private)

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
%!test  # MATLAB parity: the surface a compact regression tree reports
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! CMdl = compact (RegressionTree (X, MPG));
%! assert_equal (class (CMdl), 'CompactRegressionTree');
%! assert_equal (numel (properties (CMdl)), 28);
%! assert_equal (CMdl.NumNodes, 37);
%! assert_equal (CMdl.ResponseName, 'Y');
%! assert_equal (CMdl.PredictorNames, {'x1', 'x2', 'x3'});
%! assert_equal (CMdl.ResponseTransform, 'none');
%! assert_equal (CMdl.CategoricalPredictors, []);

%!test  # The compact tree carries the node table the model it came from has
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.Children, Mdl.Children);
%! assert_equal (CMdl.Parent, Mdl.Parent);
%! assert_equal (CMdl.NodeSize, Mdl.NodeSize);
%! assert_equal (CMdl.NodeMean, Mdl.NodeMean, 1e-15);
%! assert_equal (CMdl.NodeError, Mdl.NodeError, 1e-15);
%! assert_equal (CMdl.NodeProbability, Mdl.NodeProbability, 1e-15);
%! assert_equal (CMdl.NodeRisk, Mdl.NodeRisk, 1e-15);
%! assert_equal (CMdl.PruneList, Mdl.PruneList);
%! assert_equal (CMdl.PruneAlpha, Mdl.PruneAlpha, 1e-12);
%! assert_equal (CMdl.CutPredictor, Mdl.CutPredictor);
%! assert_equal (CMdl.CutType, Mdl.CutType);

%!test  # MATLAB parity: what a tree with no categories and no surrogates holds
%! load carsmall
%! CMdl = compact (RegressionTree ([Weight, Cylinders, Horsepower], MPG));
%! assert_equal (size (CMdl.CutCategories), [37, 2]);
%! assert_equal (size (CMdl.CategoricalSplit), [0, 0]);
%! assert_equal (size (CMdl.SurrogateCutPredictor), [0, 1]);
%! assert_equal (size (CMdl.SurrogateCutPoint), [0, 0]);

%!test  # MATLAB parity: predict answers exactly as the full model does
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG);
%! CMdl = compact (Mdl);
%! [y1, n1] = predict (Mdl, X);
%! [y2, n2] = predict (CMdl, X);
%! assert_equal ({y1, n1}, {y2, n2});
%! assert_equal (predict (CMdl, X([1, 20, 60], :))', ...
%!               [17.25, 12.3333333333333, 29.1], 1e-12);

%!test  # MATLAB parity: loss, importance and the node ranges of a compact tree
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! CMdl = compact (RegressionTree (X, MPG));
%! assert_equal (loss (CMdl, X, MPG), 5.5828069902791, 1e-12);
%! assert_equal (predictorImportance (CMdl), [2.5904, 0.1006, 0.5499], 1e-4);
%! r = nodeVariableRange (CMdl, 4);
%! assert_equal (r.x1, [-Inf, 3085.5], 1e-12);
%! assert_equal (r.x3, [-Inf, 89], 1e-12);
%! assert_equal (fieldnames (nodeVariableRange (CMdl, 1)), cell (0, 1));

%!test  # MATLAB parity: the text form of a compact tree
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! CMdl = compact (RegressionTree (X, MPG, 'MinLeafSize', 15));
%! lines = strsplit (strtrim (evalc ('view (CMdl)')), "\n");
%! assert_equal (numel (lines), 10);
%! assert_equal (lines{1}, 'Decision tree for regression');
%! assert_equal (lines{2}, ...
%!   '1  if x1<3085.5 then node 2 elseif x1>=3085.5 then node 3 else 23.7181');
%! assert_equal (lines{6}, '5  fit = 24.0882');

%!test  # A node that holds rows back is discounted on the compact tree too
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! Mdl = RegressionTree (X, MPG, 'MinLeafSize', 15);
%! assert_equal (predictorImportance (compact (Mdl)), ...
%!               predictorImportance (Mdl), 1e-15);
%! assert_equal (predictorImportance (compact (Mdl)), ...
%!               [11.253155105041, 0, 1.49831354224175], 1e-12);

%!test  # The response transform travels with the compact model
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! CMdl = compact (RegressionTree (X, MPG, 'ResponseTransform', 'exp'));
%! assert_equal (CMdl.ResponseTransform, 'exp');
%! assert_equal (predict (CMdl, X(1, :)), exp (17.25), 1e-10);
%! CMdl.ResponseTransform = @(y) 2 * y;
%! assert_equal (predict (CMdl, X(1, :)), 34.5, 1e-12);

%!test  # A compact model saved and loaded answers exactly as it did
%! load carsmall
%! X = [Weight, Cylinders, Horsepower];
%! CMdl = compact (RegressionTree (X, MPG));
%! fname = tempname ();
%! unwind_protect
%!   savemodel (CMdl, fname);
%!   New = loadmodel (fname);
%!   assert_equal (class (New), 'CompactRegressionTree');
%!   assert_equal (New.NumNodes, CMdl.NumNodes);
%!   assert_equal (New.NodeRisk, CMdl.NodeRisk, 1e-15);
%!   assert_equal (New.CutPredictor, CMdl.CutPredictor);
%!   assert_equal (predict (New, X), predict (CMdl, X));
%!   assert_equal (predictorImportance (New), ...
%!                 predictorImportance (CMdl), 1e-15);
%! unwind_protect_cleanup
%!   delete (fname);
%! end_unwind_protect

## Test input validation
%!error<CompactRegressionTree: too few input arguments.>
%! CompactRegressionTree ()
%!error<CompactRegressionTree: MDL must be a RegressionTree object.>
%! CompactRegressionTree (5)
%!error<CompactRegressionTree.predict: too few input arguments.>
%! predict (compact (RegressionTree (ones (4, 2), (1:4)')))
%!error<CompactRegressionTree.predict: XC is empty.>
%! predict (compact (RegressionTree (ones (4, 2), (1:4)')), [])
%!error<CompactRegressionTree.predict: XC must be a real numeric matrix.>
%! predict (compact (RegressionTree (ones (4, 2), (1:4)')), 'a')
%!error<CompactRegressionTree.predict: XC must have the same number of predictors as the trained model.>
%! predict (compact (RegressionTree (ones (4, 2), (1:4)')), ones (2, 5))
%!error<CompactRegressionTree.loss: too few input arguments.>
%! loss (compact (RegressionTree (ones (4, 2), (1:4)')), ones (4, 2))
%!error<CompactRegressionTree.loss: unsupported 'LossFun' value.>
%! loss (compact (RegressionTree (ones (4, 2), (1:4)')), ones (4, 2), ...
%!       (1:4)', 'LossFun', 'mad')
%!error<CompactRegressionTree.loss: 'Weights' must have one element per observation.>
%! loss (compact (RegressionTree (ones (4, 2), (1:4)')), ones (4, 2), ...
%!       (1:4)', 'Weights', [1, 2])
%!error<CompactRegressionTree.nodeVariableRange: NODE must be a positive integer no greater than the number of nodes in the tree.>
%! nodeVariableRange (compact (RegressionTree (ones (4, 2), (1:4)')), 999)
%!error<CompactRegressionTree.savemodel: too few input arguments.>
%! savemodel (compact (RegressionTree (ones (4, 2), (1:4)')))
%!error<CompactRegressionTree.savemodel: FNAME must be a character vector.>
%! savemodel (compact (RegressionTree (ones (4, 2), (1:4)')), 5)
%!error<CompactRegressionTree: unrecognized 'ResponseTransform' function.>
%! CMdl = compact (RegressionTree (ones (4, 2), (1:4)'));
%! CMdl.ResponseTransform = 'bogus';
