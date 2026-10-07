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

classdef CompactTreeBagger < PredictiveModel
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactTreeBagger
  ##
  ## Compact ensemble of bagged decision trees
  ##
  ## A @code{CompactTreeBagger} object carries the trees a @code{TreeBagger}
  ## ensemble grew and what prediction needs, but not the observations it was
  ## fitted on, and so none of the out-of-bag information.  It predicts new
  ## data identically to the ensemble it came from.
  ##
  ## Create one with the @code{compact} method of a @code{TreeBagger} object.
  ## Two compact ensembles fitted on the same classes are joined with
  ## @code{combine}.
  ##
  ## @seealso{TreeBagger, TreeBagger.compact, CompactClassificationTree,
  ## CompactRegressionTree}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} Method
    ##
    ## Type of the ensemble
    ##
    ## @qcode{'classification'} or @qcode{'regression'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    Method = 'classification';

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} NumTrees
    ##
    ## Number of trees
    ##
    ## A nonnegative integer, the number of trees in the ensemble.  This
    ## property is read-only.
    ##
    ## @end deftp
    NumTrees = 0;

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} Trees
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
    ## @deftp {CompactTreeBagger} {property} ClassNames
    ##
    ## Names of the classes
    ##
    ## The classes of a classification ensemble, in the type of the response
    ## it was fitted on and in the order its scores are laid out.  Empty for a
    ## regression ensemble.  This property is read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} DefaultYfit
    ##
    ## Prediction for an observation no tree may answer for
    ##
    ## For classification, the class of greatest prior probability, in the
    ## type of the class names, or the missing label of that type after
    ## @code{setDefaultYfit} with @qcode{''}.  For regression, the weighted
    ## mean of the training response unless set otherwise.  This property is
    ## read-only; change it with @code{setDefaultYfit}.
    ##
    ## @end deftp
    DefaultYfit = [];

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors naming the columns of the predictor
    ## data.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} CategoricalPredictors
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
    ## @deftp {CompactTreeBagger} {property} DeltaCriterionDecisionSplit
    ##
    ## Split criterion contributions of the predictors
    ##
    ## A row vector with one element per predictor, the mean over the trees
    ## of each tree's @code{predictorImportance}.  This property is read-only.
    ##
    ## @end deftp
    DeltaCriterionDecisionSplit = [];

    ## -*- texinfo -*-
    ## @deftp {CompactTreeBagger} {property} NumPredictorSplit
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
    ## @deftp {CompactTreeBagger} {property} SurrogateAssociation
    ##
    ## Predictive association between the predictors
    ##
    ## A square matrix with one row and one column per predictor.  The trees
    ## grow no surrogate splits, so it is the identity matrix.  This property
    ## is read-only.
    ##
    ## @end deftp
    SurrogateAssociation = [];

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    ResponseName = 'Y';  # name the table gave the response, which
                         # a table call looks up; hidden as MATLAB's
                         # TreeBagger carries no ResponseName
    TreeClassIdx = {};   # columns of ClassNames each tree's scores fill
    DefaultIndex = 0;    # index of DefaultYfit into ClassNames, 0 if missing
    DefaultScore = [];   # scores of an observation no tree may answer for
  endproperties

  methods (Hidden)

    function disp (this)
      fprintf ('\n  CompactTreeBagger\n\n');
      fprintf ('%16s: %s\n', 'Method', this.Method);
      fprintf ('%16s: %d\n', 'NumTrees', this.NumTrees);
      fprintf ('%16s: %d\n', 'NumPredictors', numel (this.PredictorNames));
      if (strcmp (this.Method, 'classification'))
        fprintf ('%16s: %s\n', 'ClassNames', ...
                 classNameListing (this.ClassNames));
      endif
      fprintf ('\n');
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactTreeBagger} {@var{obj} =} CompactTreeBagger (@var{B})
    ##
    ## Create a @code{CompactTreeBagger} object.
    ##
    ## @var{B} is the @code{TreeBagger} object to compact.  The documented way
    ## to reach this constructor is the @code{compact} method.
    ##
    ## @seealso{TreeBagger, TreeBagger.compact}
    ## @end deftypefn
    function this = CompactTreeBagger (B)

      ## Input validation
      if (nargin < 1)
        error ("CompactTreeBagger: too few input arguments.");
      endif
      if (! isa (B, 'TreeBagger'))
        error ("CompactTreeBagger: B must be a TreeBagger object.");
      endif

      ## The levels a predictor read from a table was coded through
      ## travel with the model, so a compact one still reads a table, and
      ## so does the name the table gave the response
      this.PredictorLevels = B.PredictorLevels;
      this.ResponseName = B.ResponseName;

      this.Method = B.Method;
      this.NumTrees = B.NumTrees;
      this.Trees = B.Trees;
      this.ClassNames = B.ClassNames;
      this.DefaultYfit = B.DefaultYfit;
      this.PredictorNames = B.PredictorNames;
      this.CategoricalPredictors = B.CategoricalPredictors;
      this.DeltaCriterionDecisionSplit = B.DeltaCriterionDecisionSplit;
      this.NumPredictorSplit = B.NumPredictorSplit;
      this.SurrogateAssociation = B.SurrogateAssociation;
      this.TreeClassIdx = B.TreeClassIdx;
      this.DefaultIndex = B.DefaultIndex;
      this.DefaultScore = B.DefaultScore;

    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactTreeBagger} {@var{label} =} predict (@var{obj}, @var{X})
    ## @deftypefnx {CompactTreeBagger} {[@var{label}, @var{scores}] =} predict (@dots{})
    ## @deftypefnx {CompactTreeBagger} {[@var{label}, @var{scores}, @var{stdevs}] =} predict (@dots{})
    ## @deftypefnx {CompactTreeBagger} {[@var{Yfit}, @var{stdevs}] =} predict (@dots{})
    ## @deftypefnx {CompactTreeBagger} {@dots{} =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict responses with a compact bagged ensemble.
    ##
    ## For a classification ensemble, @var{label} holds the predicted class of
    ## each row of @var{X}, in the type of @code{ClassNames}.  @var{scores} is
    ## an @math{NxK} matrix, the weighted average over the trees of the class
    ## probability each tree gives, and @var{stdevs} holds their standard
    ## deviations over the trees.  The label is the class of highest score,
    ## whatever cost matrix the trees were grown with, as MATLAB documents.
    ##
    ## For a regression ensemble, @var{Yfit} is the weighted average of the
    ## trees' predictions and @var{stdevs} their standard deviation.
    ##
    ## The standard deviations are population deviations, taken over the
    ## trees that answer for the observation with their weights.  An
    ## observation no tree may answer for is given @code{DefaultYfit}, with
    ## the prior as its scores and @code{NaN} as its deviations.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Trees'} @tab @tab @qcode{'all'} (default) or a vector of
    ## indices of the trees to use.
    ## @item @qcode{'TreeWeights'} @tab @tab A nonnegative vector with one
    ## weight per tree used.  The default weighs them equally.
    ## @item @qcode{'UseInstanceForTree'} @tab @tab An @math{NxNumTrees}
    ## logical matrix saying which tree may answer for which observation.  The
    ## default lets every tree answer for every observation.
    ## @end multitable
    ##
    ## MATLAB returns the labels as a cell array of character vectors whatever
    ## the type of the response; they are returned here in the type of
    ## @code{ClassNames}, as by every other classifier in this package.
    ##
    ##
    ## The new data may be a table, whose variables are matched to the
    ## predictors the model was fitted on by name and not by position:
    ## one the model was not fitted on is passed over, one it needs and
    ## cannot find is named, and a value holding a level is coded as that
    ## level was coded at fitting.
    ## @seealso{CompactTreeBagger, CompactTreeBagger.error, TreeBagger.predict}
    ## @end deftypefn
    function [Yfit, scores, stdevs] = predict (this, X, varargin)

      if (nargin < 2)
        error ("CompactTreeBagger.predict: too few input arguments.");
      endif

      ## A table is read by the names the model was fitted on
      X = tableColumns (this, 'CompactTreeBagger.predict', X);
      [Yfit, scores, stdevs] = bagPredict (this, X, varargin, ...
                                           'CompactTreeBagger.predict', []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactTreeBagger} {@var{err} =} error (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactTreeBagger} {@var{err} =} error (@dots{}, @var{name}, @var{value})
    ## @deftypefnx {CompactTreeBagger} {@var{err} =} error (@var{obj}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {CompactTreeBagger} {@var{err} =} error (@var{obj}, @var{Tbl})
    ##
    ## Misclassification probability or mean squared error of the ensemble.
    ##
    ## For a classification ensemble @var{err} is the weighted share of the
    ## rows of @var{X} whose predicted class differs from @var{Y}; for a
    ## regression ensemble it is the weighted mean squared error.  By default
    ## @var{err} is a column with one element per tree, the error of the
    ## first tree, then of the first two, and so on.
    ##
    ## An observation no tree may answer for is predicted as
    ## @code{DefaultYfit} and counted; when @code{DefaultYfit} is the missing
    ## label it has no prediction and is left out.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Mode'} @tab @tab @qcode{'cumulative'} (default),
    ## @qcode{'individual'} for the error of each tree on its own, or
    ## @qcode{'ensemble'} for a single error over every tree used.
    ## @item @qcode{'Trees'} @tab @tab @qcode{'all'} (default) or a vector of
    ## indices of the trees to use, in the order they are accumulated.
    ## @item @qcode{'TreeWeights'} @tab @tab A nonnegative vector with one
    ## weight per tree used.  It may not be given in @qcode{'individual'} mode.
    ## @item @qcode{'UseInstanceForTree'} @tab @tab An @math{NxNumTrees}
    ## logical matrix saying which tree may answer for which observation.
    ## @item @qcode{'Weights'} @tab @tab A nonnegative vector with one weight
    ## per observation.  The default is uniform.
    ## @end multitable
    ##
    ## @var{X} may also be a table @var{Tbl}, whose variables are matched to
    ## the predictors the model was fitted on by name and not by position.
    ## @code{error (@var{obj}, @var{Tbl}, @var{ResponseVarName})} takes the
    ## response from the variable @var{ResponseVarName} names, and
    ## @code{error (@var{obj}, @var{Tbl})} from the variable the model was
    ## fitted on.  The response may also be given beside the table as
    ## @var{Y}.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.predict,
    ## CompactTreeBagger.meanMargin, TreeBagger.oobError}
    ## @end deftypefn
    function err = error (this, X, Y, varargin)

      if (nargin < 3 && ! (nargin > 1 && istable (X)))
        error ("CompactTreeBagger.error: too few input arguments.");
      endif

      ## A table carries the response: named in the call, given beside
      ## the table, or the variable the model was fitted on
      if (nargin < 3)
        Y = [];
      endif
      [X, Y, varargin] = tableResponse (this, 'error', X, Y, ...
                                        varargin, nargin > 2);
      err = bagLoss ('error', this, X, Y, varargin, ...
                     'CompactTreeBagger.error', [], []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactTreeBagger} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactTreeBagger} {@var{m} =} margin (@dots{}, @var{name}, @var{value})
    ## @deftypefnx {CompactTreeBagger} {@var{m} =} margin (@var{obj}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {CompactTreeBagger} {@var{m} =} margin (@var{obj}, @var{Tbl})
    ##
    ## Classification margin of each observation.
    ##
    ## The margin is the score of the true class less the highest score among
    ## the other classes.  By default @var{m} has one row per observation and
    ## one column per tree, the margin of the first tree, then of the first
    ## two, and so on.  An observation left without a prediction, as described
    ## under @code{CompactTreeBagger.error}, has a @code{NaN} margin.
    ##
    ## @qcode{'Mode'}, @qcode{'Trees'}, @qcode{'TreeWeights'} and
    ## @qcode{'UseInstanceForTree'} are taken as by
    ## @code{CompactTreeBagger.error}.  A regression ensemble has no margins.
    ##
    ## @var{X} may also be a table @var{Tbl}, whose variables are matched to
    ## the predictors the model was fitted on by name and not by position.
    ## @code{margin (@var{obj}, @var{Tbl}, @var{ResponseVarName})} takes the
    ## response from the variable @var{ResponseVarName} names, and
    ## @code{margin (@var{obj}, @var{Tbl})} from the variable the model was
    ## fitted on.  The response may also be given beside the table as
    ## @var{Y}.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.meanMargin,
    ## CompactTreeBagger.error}
    ## @end deftypefn
    function m = margin (this, X, Y, varargin)

      if (nargin < 3 && ! (nargin > 1 && istable (X)))
        error ("CompactTreeBagger.margin: too few input arguments.");
      endif

      ## A table carries the response: named in the call, given beside
      ## the table, or the variable the model was fitted on
      if (nargin < 3)
        Y = [];
      endif
      [X, Y, varargin] = tableResponse (this, 'margin', X, Y, ...
                                        varargin, nargin > 2);
      m = bagLoss ('margin', this, X, Y, varargin, ...
                   'CompactTreeBagger.margin', [], []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactTreeBagger} {@var{mm} =} meanMargin (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactTreeBagger} {@var{mm} =} meanMargin (@dots{}, @var{name}, @var{value})
    ## @deftypefnx {CompactTreeBagger} {@var{mm} =} meanMargin (@var{obj}, @var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {CompactTreeBagger} {@var{mm} =} meanMargin (@var{obj}, @var{Tbl})
    ##
    ## Weighted mean classification margin.
    ##
    ## @var{mm} is the weighted mean over the observations of the margins
    ## @code{CompactTreeBagger.margin} returns, a row with one element per
    ## tree by default.  The Name-Value arguments are those of
    ## @code{CompactTreeBagger.error}, @qcode{'Weights'} included.
    ##
    ## @var{X} may also be a table @var{Tbl}, whose variables are matched to
    ## the predictors the model was fitted on by name and not by position.
    ## @code{meanMargin (@var{obj}, @var{Tbl}, @var{ResponseVarName})}
    ## takes the response from the variable @var{ResponseVarName} names, and
    ## @code{meanMargin (@var{obj}, @var{Tbl})} from the variable the model
    ## was fitted on.  The response may also be given beside the table as
    ## @var{Y}.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.margin}
    ## @end deftypefn
    function mm = meanMargin (this, X, Y, varargin)

      if (nargin < 3 && ! (nargin > 1 && istable (X)))
        error ("CompactTreeBagger.meanMargin: too few input arguments.");
      endif

      ## A table carries the response: named in the call, given beside
      ## the table, or the variable the model was fitted on
      if (nargin < 3)
        Y = [];
      endif
      [X, Y, varargin] = tableResponse (this, 'meanMargin', X, Y, ...
                                        varargin, nargin > 2);
      mm = bagLoss ('meanMargin', this, X, Y, varargin, ...
                    'CompactTreeBagger.meanMargin', [], []);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactTreeBagger} {@var{C} =} combine (@var{C1}, @var{C2})
    ##
    ## Join two compact ensembles.
    ##
    ## @var{C} holds the trees of @var{C1} followed by those of @var{C2}.  The
    ## two must be of the same type and, for classification, have the same
    ## class names and the same default scores.  @var{C} keeps the
    ## @code{DefaultYfit} of @var{C1}.
    ##
    ## @seealso{CompactTreeBagger, TreeBagger.append}
    ## @end deftypefn
    function this = combine (this, other)

      if (nargin < 2)
        error ("CompactTreeBagger.combine: too few input arguments.");
      endif
      if (! isa (other, 'CompactTreeBagger'))
        error (strcat ("CompactTreeBagger.combine: C2 must be a", ...
                       " CompactTreeBagger object."));
      endif
      if (! strcmp (this.Method, other.Method))
        error (strcat ("CompactTreeBagger.combine: the two ensembles must", ...
                       " be of the same type."));
      endif
      if (! (isequal (this.ClassNames, other.ClassNames)
             && isequal (this.DefaultScore, other.DefaultScore)))
        error (strcat ("CompactTreeBagger.combine: the two ensembles must", ...
                       " have the same classes and priors."));
      endif
      if (! isequal (this.PredictorNames, other.PredictorNames))
        error (strcat ("CompactTreeBagger.combine: the two ensembles must", ...
                       " have the same predictors."));
      endif
      this.Trees = [this.Trees; other.Trees];
      this.TreeClassIdx = [this.TreeClassIdx; other.TreeClassIdx];
      this.NumTrees = numel (this.Trees);
      [this.DeltaCriterionDecisionSplit, this.NumPredictorSplit] = ...
        bagSplitStats (this.Trees, numel (this.PredictorNames));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactTreeBagger} {@var{C} =} setDefaultYfit (@var{obj}, @var{Yfit})
    ##
    ## Set the prediction for an observation no tree may answer for.
    ##
    ## For a classification ensemble @var{Yfit} is @qcode{'MostPopular'}, the
    ## class of greatest prior probability, or @qcode{''}, the missing label,
    ## which leaves such an observation without a prediction and out of the
    ## error and the mean margin.  A logical response has no missing label.
    ## For a regression ensemble @var{Yfit} is a numeric scalar.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.predict}
    ## @end deftypefn
    function this = setDefaultYfit (this, Yfit)

      if (nargin < 2)
        error ("CompactTreeBagger.setDefaultYfit: too few input arguments.");
      endif
      if (strcmp (this.Method, 'classification'))
        if (ischar (Yfit) && isempty (Yfit))
          if (islogical (this.ClassNames))
            error (strcat ("CompactTreeBagger.setDefaultYfit: a logical", ...
                           " response has no missing label, so YFIT must", ...
                           " be 'MostPopular'."));
          endif
          this.DefaultIndex = 0;
          this.DefaultYfit = missingLabels (this.ClassNames, 1);
        elseif (ischar (Yfit) && strcmpi (Yfit, 'MostPopular'))
          [~, this.DefaultIndex] = max (this.DefaultScore);
          this.DefaultYfit = labelsFromIndex (this.ClassNames, ...
                                              this.DefaultIndex);
        else
          error (strcat ("CompactTreeBagger.setDefaultYfit: YFIT must be", ...
                         " '' or 'MostPopular' for a classification", ...
                         " ensemble."));
        endif
      else
        if (! (isnumeric (Yfit) && isscalar (Yfit) && isreal (Yfit)))
          error (strcat ("CompactTreeBagger.setDefaultYfit: YFIT must be", ...
                         " a real numeric scalar for a regression", ...
                         " ensemble."));
        endif
        this.DefaultYfit = double (Yfit);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactTreeBagger} {@var{prox} =} proximity (@var{obj}, @var{X})
    ##
    ## Proximity matrix of the observations.
    ##
    ## @var{prox} is a symmetric @math{NxN} matrix, @math{N} being the number
    ## of rows of @var{X}, whose element @math{(i,j)} is the share of the
    ## trees that bring observations @math{i} and @math{j} to the same leaf.
    ## Its diagonal holds ones.
    ##
    ## @var{X} may also be a table, whose variables are matched to the
    ## predictors the model was fitted on by name and not by position: one
    ## the model was not fitted on is passed over, one it needs and cannot
    ## find is named, and a value holding a level is coded as that level was
    ## coded at fitting.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.outlierMeasure,
    ## CompactTreeBagger.mdsprox, TreeBagger.fillprox}
    ## @end deftypefn
    function prox = proximity (this, X)

      if (nargin < 2)
        error ("CompactTreeBagger.proximity: too few input arguments.");
      endif
      o = proxArgs (this, X, {}, {}, 'CompactTreeBagger.proximity');
      prox = o.P;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactTreeBagger} {@var{out} =} outlierMeasure (@var{obj}, @var{X})
    ## @deftypefnx {CompactTreeBagger} {@var{out} =} outlierMeasure (@dots{}, @var{name}, @var{value})
    ##
    ## Outlier measure of each observation.
    ##
    ## The raw measure of an observation is the size of its class divided by
    ## the sum of its squared proximities to the observations of the class,
    ## itself included, all the observations forming one class when no labels
    ## are given.  @var{out} is a column holding, for each observation, the
    ## absolute deviation of its raw measure from the median of its class,
    ## divided by the median absolute deviation of the class.  A large value
    ## marks an observation that the trees seldom group with the rest of its
    ## class.  As in MATLAB, a class whose median absolute deviation is zero
    ## gives the raw measures themselves, and a class of one or two
    ## observations gives zeros.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.2 0.02 0.78
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Data'} @tab @tab @qcode{'predictors'} (default), for
    ## @var{X} holding predictor data, or @qcode{'proximity'}, for @var{X}
    ## holding a proximity matrix such as @code{proximity} returns.
    ## @item @qcode{'Labels'} @tab @tab The class label of each observation,
    ## each one of @code{ClassNames}.  Classification only.
    ## @end multitable
    ##
    ## @var{X} holding predictor data may also be a table, whose variables
    ## are matched to the predictors the model was fitted on by name and not
    ## by position: one the model was not fitted on is passed over, one it
    ## needs and cannot find is named, and a value holding a level is coded as
    ## that level was coded at fitting.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.proximity,
    ## TreeBagger.OutlierMeasure}
    ## @end deftypefn
    function out = outlierMeasure (this, X, varargin)

      if (nargin < 2)
        error ("CompactTreeBagger.outlierMeasure: too few input arguments.");
      endif
      o = proxArgs (this, X, varargin, {'Data', 'Labels'}, ...
                    'CompactTreeBagger.outlierMeasure');
      out = bagOutlier (o.P, o.g);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactTreeBagger} {[@var{S}, @var{E}] =} mdsprox (@var{obj}, @var{X})
    ## @deftypefnx {CompactTreeBagger} {[@var{S}, @var{E}] =} mdsprox (@dots{}, @var{name}, @var{value})
    ##
    ## Multidimensional scaling of the proximity matrix.
    ##
    ## Applies classical multidimensional scaling, as @code{cmdscale} does, to
    ## the distances @code{1 - @var{prox}}, @var{prox} being the proximity
    ## matrix of the rows of @var{X}.  @var{S} holds the scaled coordinates,
    ## one column per positive eigenvalue, and @var{E} the eigenvalues.
    ##
    ## Name-Value arguments:
    ##
    ## @multitable @columnfractions 0.25 0.02 0.73
    ## @headitem @var{Name} @tab @tab @var{Value}
    ## @item @qcode{'Data'} @tab @tab @qcode{'predictors'} (default) or
    ## @qcode{'proximity'}, as for @code{outlierMeasure}.
    ## @item @qcode{'Colors'} @tab @tab A character vector with one color
    ## letter per class.  When given, the scaled coordinates are drawn as
    ## overlaid scatter plots, one per class, and a class beyond the number of
    ## letters is not drawn.
    ## @item @qcode{'Labels'} @tab @tab The class label of each observation,
    ## each one of @code{ClassNames}.  Classification only.  Without labels
    ## every observation is drawn in the first color.
    ## @item @qcode{'MDSCoordinates'} @tab @tab Two or three indices of the
    ## columns of @var{S} to draw.  The default is @code{[1, 2]}.  They must
    ## not exceed the number of columns of @var{S} even when nothing is drawn,
    ## as in MATLAB, whose documentation says otherwise.
    ## @end multitable
    ##
    ## @var{X} holding predictor data may also be a table, whose variables
    ## are matched to the predictors the model was fitted on by name and not
    ## by position: one the model was not fitted on is passed over, one it
    ## needs and cannot find is named, and a value holding a level is coded as
    ## that level was coded at fitting.
    ##
    ## @seealso{CompactTreeBagger, CompactTreeBagger.proximity, cmdscale,
    ## TreeBagger.mdsprox}
    ## @end deftypefn
    function [S, E] = mdsprox (this, X, varargin)

      if (nargin < 2)
        error ("CompactTreeBagger.mdsprox: too few input arguments.");
      endif
      o = proxArgs (this, X, varargin, ...
                    {'Data', 'Colors', 'Labels', 'MDSCoordinates'}, ...
                    'CompactTreeBagger.mdsprox');
      [S, E, errmsg] = bagMds (o.P, o.g, o.colors, o.coords);
      if (! isempty (errmsg))
        error ("CompactTreeBagger.mdsprox: %s", errmsg);
      endif

    endfunction

  endmethods

endclassdef

## The proximity matrix and the Name-Value arguments of the proximity methods.
function o = proxArgs (M, X, args, allowed, caller)

  if (mod (numel (args), 2) != 0)
    error ("%s: name-value arguments must be in pairs.", caller);
  endif
  o = struct ('P', [], 'g', [], 'colors', '', 'coords', [1, 2]);
  isprox = false;
  labels = [];
  for i = 1:2:numel (args)
    name = args{i};
    val = args{i+1};
    if (! (ischar (name) && any (strcmpi (name, allowed))))
      error ("%s: invalid parameter name in optional pair arguments.", ...
             caller);
    endif
    switch (tolower (name))
      case 'data'
        if (! (ischar (val)
               && any (strcmpi (val, {'predictors', 'proximity'}))))
          error ("%s: 'Data' must be 'predictors' or 'proximity'.", caller);
        endif
        isprox = strcmpi (val, 'proximity');
      case 'labels'
        labels = {val};
      case 'colors'
        o.colors = val;
      case 'mdscoordinates'
        o.coords = val;
    endswitch
  endfor

  if (isprox)
    if (! (isnumeric (X) && isreal (X) && issquare (X) && ! isempty (X)))
      error (strcat ("%s: X must be a square real numeric matrix when", ...
                     " 'Data' is 'proximity'."), caller);
    endif
    o.P = double (X);
  else
    ## A table is read by the names the model was fitted on
    X = tableToMatrix (M.PredictorNames, M.PredictorLevels, X, caller);
    if (! (isnumeric (X) && isreal (X) && ismatrix (X)))
      error ("%s: X must be a real numeric matrix.", caller);
    endif
    if (columns (X) != numel (M.PredictorNames))
      error ("%s: X must have one column per predictor.", caller);
    endif
    o.P = bagProximity (M, X, 1:M.NumTrees);
  endif

  if (! isempty (labels))
    lab = labels{1};
    if (ischar (lab))
      nL = rows (lab);
    else
      nL = numel (lab);
    endif
    if (! (isnumeric (lab) || islogical (lab) || ischar (lab)
           || iscellstr (lab) || isa (lab, 'categorical')
           || isa (lab, 'string')) || nL != rows (o.P))
      error ("%s: 'Labels' must hold one class label per observation.", ...
             caller);
    endif
    if (! strcmp (M.Method, 'classification'))
      error ("%s: 'Labels' cannot be used with a regression ensemble.", ...
             caller);
    endif
    [o.g, errmsg] = labelIndices (M.ClassNames, lab);
    if (! isempty (errmsg))
      error (strcat ("%s: 'Labels' must hold only classes the ensemble", ...
                     " was trained on."), caller);
    endif
  endif

endfunction

%!demo
%! ## Compact a random forest and classify new flowers with it.  The compact
%! ## ensemble keeps the trees but not the training data.
%! load fisheriris
%! rng (42);
%! C = compact (TreeBagger (30, meas, species));
%! label = predict (C, [5.0, 3.4, 1.5, 0.2; 6.7, 3.0, 5.2, 2.3])

%!test  # a compact ensemble predicts as the ensemble it came from
%! load fisheriris
%! rng (1);
%! B = TreeBagger (5, meas, species);
%! C = compact (B);
%! assert_equal (class (C), 'CompactTreeBagger');
%! [la, sa] = predict (B, meas);
%! [lc, sc] = predict (C, meas);
%! assert_equal (lc, la);
%! assert_equal (sc, sa);
%! assert_equal (C.ClassNames, B.ClassNames);
%! assert_equal (C.DefaultYfit, B.DefaultYfit);

%!test  # a compact ensemble keeps the ensemble's categorical predictors
%! Xc = [1, 2; 2, 3; 3, 4; 1, 5; 2, 6; 3, 7; 1, 8; 2, 9];
%! yc = [1; 1; 2; 2; 1; 2; 1; 2];
%! C = compact (TreeBagger (3, Xc, yc, 'CategoricalPredictors', 1));
%! assert_equal (C.CategoricalPredictors, 1);
%! C = compact (TreeBagger (3, Xc, yc));
%! assert_equal (C.CategoricalPredictors, []);

%!test  # MATLAB parity: combining keeps the first ensemble's default
%! load fisheriris
%! rng (1);
%! C1 = setDefaultYfit (compact (TreeBagger (3, meas, species)), '');
%! C2 = compact (TreeBagger (4, meas, species));
%! C = combine (C1, C2);
%! assert_equal (C.NumTrees, 7);
%! assert_equal (C.DefaultYfit, {''});

%!test  # a missing default leaves an observation without a prediction
%! load fisheriris
%! rng (1);
%! C = setDefaultYfit (compact (TreeBagger (4, meas, species)), '');
%! U = true (150, 4);
%! U(1,:) = false;
%! label = predict (C, meas(1:2,:), 'UseInstanceForTree', U(1:2,:));
%! assert_equal (label(1), {''});
%! m = margin (C, meas(1:2,:), species(1:2), 'Mode', 'ensemble', ...
%!             'UseInstanceForTree', U(1:2,:));
%! assert_equal (isnan (m), [true; false]);
%! miss = ! strcmp (predict (C, meas(2:end,:)), species(2:end));
%! e = error (C, meas, species, 'Mode', 'ensemble', 'UseInstanceForTree', U);
%! assert_equal (e, sum (miss) / 149, 1e-15);

%!test  # 'MostPopular' restores the class of greatest prior probability
%! load fisheriris
%! rng (1);
%! C = setDefaultYfit (compact (TreeBagger (2, meas, species)), '');
%! C = setDefaultYfit (C, 'MostPopular');
%! assert_equal (C.DefaultYfit, {'setosa'});

%!test  # a regression default is any numeric scalar
%! load fisheriris
%! rng (1);
%! C = compact (TreeBagger (2, meas(:,2:4), meas(:,1), 'Method', 'regression'));
%! C = setDefaultYfit (C, 3);
%! assert_equal (C.DefaultYfit, 3);

%!test  # MATLAB parity: combining adds up the split counts
%! load fisheriris
%! rng (1);
%! C1 = compact (TreeBagger (3, meas, species));
%! C2 = compact (TreeBagger (4, meas, species));
%! C = combine (C1, C2);
%! assert_equal (C.NumPredictorSplit, ...
%!               C1.NumPredictorSplit + C2.NumPredictorSplit, 1e-14);
%! assert_equal (C.DeltaCriterionDecisionSplit, ...
%!               (3 * C1.DeltaCriterionDecisionSplit ...
%!                + 4 * C2.DeltaCriterionDecisionSplit) / 7, 1e-14);
%! assert_equal (C.SurrogateAssociation, eye (4));

%!test  # MATLAB parity: proximity is the share of trees sharing a leaf
%! load fisheriris
%! rng (1);
%! C = compact (TreeBagger (6, meas, species, 'MinLeafSize', 5));
%! X = meas(1:3:end,:);
%! P = zeros (rows (X));
%! for t = 1:6
%!   [~, ~, nd] = predict (C.Trees{t}, X);
%!   P += nd == nd';
%! endfor
%! prox = proximity (C, X);
%! assert_equal (prox, P / 6, 1e-15);
%! assert_equal (diag (prox), ones (50, 1));

%!test  # a regression ensemble has proximities too
%! load fisheriris
%! rng (1);
%! C = compact (TreeBagger (3, meas(:,2:4), meas(:,1), 'Method', 'regression'));
%! [~, n1] = predict (C.Trees{1}, meas(1:2,2:4));
%! [~, n2] = predict (C.Trees{2}, meas(1:2,2:4));
%! [~, n3] = predict (C.Trees{3}, meas(1:2,2:4));
%! p12 = mean ([n1(1) == n1(2), n2(1) == n2(2), n3(1) == n3(2)]);
%! assert_equal (proximity (C, meas(1:2,2:4)), [1, p12; p12, 1], 1e-15);

%!test  # MATLAB parity: the outlier measure within each class
%! load fisheriris
%! rng (1);
%! C = compact (TreeBagger (8, meas, species, 'MinLeafSize', 5));
%! P = proximity (C, meas);
%! g = grp2idx (species);
%! om = zeros (150, 1);
%! for k = 2:3
%!   raw = 50 ./ sum (P(g == k, g == k) .^ 2, 2);
%!   om(g == k) = abs (raw - median (raw)) / median (abs (raw - median (raw)));
%! endfor
%! out = outlierMeasure (C, meas, 'Labels', species);
%! assert_equal (out(g > 1), om(g > 1), 1e-12);

%!test  # MATLAB parity: a proximity matrix may be given in place of the data
%! load fisheriris
%! rng (1);
%! C = compact (TreeBagger (8, meas, species, 'MinLeafSize', 5));
%! P = proximity (C, meas);
%! assert_equal (outlierMeasure (C, P, 'Data', 'proximity', ...
%!                               'Labels', species), ...
%!               outlierMeasure (C, meas, 'Labels', species));

%!test  # MATLAB parity: without labels every observation is one class
%! C = compact (TreeBagger (1, [1; 2; 3; 4], [1; 1; 2; 2]));
%! P = [1, 0.2, 0.6, 0.1, 0.3; 0.2, 1, 0.4, 0.7, 0.5; ...
%!      0.6, 0.4, 1, 0.2, 0.8; 0.1, 0.7, 0.2, 1, 0.3; ...
%!      0.3, 0.5, 0.8, 0.3, 1];
%! assert_equal (outlierMeasure (C, P, 'Data', 'proximity'), ...
%!               [2.482051282051284; 0; 1; 1.609249646059462; ...
%!                0.531400966183575], 1e-12);

%!test  # MATLAB parity: a zero median absolute deviation gives the raw measure
%! C = compact (TreeBagger (1, [1; 2; 3; 4], [1; 1; 2; 2]));
%! P = [1, 1, 1, 0.5; 1, 1, 1, 0.5; 1, 1, 1, 0.5; 0.5, 0.5, 0.5, 1];
%! assert_equal (outlierMeasure (C, P, 'Data', 'proximity'), ...
%!               [1.230769230769231; 1.230769230769231; ...
%!                1.230769230769231; 2.285714285714286], 1e-12);
%! assert_equal (outlierMeasure (C, [1, 1, 0; 1, 1, 0; 0, 0, 1], ...
%!                               'Data', 'proximity'), [1.5; 1.5; 3]);

%!test  # MATLAB parity: a class of one or two observations gives zeros
%! load fisheriris
%! C = compact (TreeBagger (1, meas, species));
%! P = [1, 1, 1, 0, 0, 0; 1, 1, 1, 0, 0, 0; 1, 1, 1, 0, 0, 0; ...
%!      0, 0, 0, 1, 1, 1; 0, 0, 0, 1, 1, 1; 0, 0, 0, 1, 1, 1];
%! lab = {'setosa'; 'setosa'; 'virginica'; 'virginica'; 'virginica'; ...
%!        'virginica'};
%! assert_equal (outlierMeasure (C, P, 'Data', 'proximity', 'Labels', lab), ...
%!               [0; 0; 4; 4/3; 4/3; 4/3], 1e-15);
%! assert_equal (outlierMeasure (C, [1, 0.5; 0.5, 0.8], 'Data', ...
%!                               'proximity'), [0; 0]);

%!test  # MATLAB parity: scaling applies cmdscale to one less the proximity
%! load fisheriris
%! rng (1);
%! C = compact (TreeBagger (5, meas, species, 'MinLeafSize', 5));
%! X = meas(1:5:end,:);
%! [S, E] = mdsprox (C, X);
%! [S0, E0] = cmdscale (1 - proximity (C, X));
%! assert_equal (S, S0);
%! assert_equal (E, E0);
%! [S1, E1] = mdsprox (C, proximity (C, X), 'Data', 'proximity');
%! assert_equal (S1, S0);

%!test  # the scaled coordinates are drawn one class per color
%! load fisheriris
%! rng (1);
%! C = compact (TreeBagger (5, meas, species, 'MinLeafSize', 5));
%! h = figure ('visible', 'off');
%! unwind_protect
%!   mdsprox (C, meas, 'Colors', 'rb', 'Labels', species);
%!   kids = get (gca, 'children');
%!   assert_equal (numel (kids), 2);
%!   assert_equal (sort (cellfun (@numel, get (kids, 'xdata'))), [50; 50]);
%!   assert_equal (ishold (), false);
%! unwind_protect_cleanup
%!   close (h);
%! end_unwind_protect

%!test  # three coordinates are drawn in three dimensions
%! load fisheriris
%! rng (1);
%! C = compact (TreeBagger (5, meas, species, 'MinLeafSize', 5));
%! h = figure ('visible', 'off');
%! unwind_protect
%!   S = mdsprox (C, meas, 'Colors', 'k', 'MDSCoordinates', [1, 2, 3]);
%!   kids = get (gca, 'children');
%!   assert_equal (numel (kids), 1);
%!   assert_equal (get (kids, 'zdata')(:), S(:,3));
%! unwind_protect_cleanup
%!   close (h);
%! end_unwind_protect

## Test input validation
%!shared x, y, C, R
%! load fisheriris
%! x = meas;
%! y = species;
%! C = compact (TreeBagger (2, x, y));
%! R = compact (TreeBagger (2, x(:,2:4), x(:,1), 'Method', 'regression'));
%!error<CompactTreeBagger: too few input arguments.> CompactTreeBagger ()
%!error<CompactTreeBagger: B must be a TreeBagger object.> CompactTreeBagger (1)
%!error<CompactTreeBagger.predict: too few input arguments.> predict (C)
%!error<CompactTreeBagger.margin: too few input arguments.> margin (C, x)
%!error<CompactTreeBagger.meanMargin: too few input arguments.> ...
%! meanMargin (C, x)
## Octave's test.m cuts a message through its first 'error:' before matching
## it, so these patterns hold the message after the method's prefix.
%!error<too few input arguments.> error (C, x)
%!error<CompactTreeBagger.combine: too few input arguments.> combine (C)
%!error<CompactTreeBagger.combine: C2 must be a CompactTreeBagger object.> ...
%! combine (C, 1)
%!error<CompactTreeBagger.combine: the two ensembles must be of the same type.> ...
%! combine (C, R)
%!error<CompactTreeBagger.combine: the two ensembles must have the same classes and priors.> ...
%! combine (C, compact (TreeBagger (1, x, y, 'Prior', 'uniform', ...
%!          'ClassNames', {'virginica'; 'versicolor'; 'setosa'})))
%!error<CompactTreeBagger.combine: the two ensembles must have the same predictors.> ...
%! combine (C, compact (TreeBagger (1, x, y, 'PredictorNames', ...
%!                                  {'a', 'b', 'c', 'd'})))
%!error<CompactTreeBagger.setDefaultYfit: too few input arguments.> ...
%! setDefaultYfit (C)
%!error<CompactTreeBagger.setDefaultYfit: a logical response has no missing label, so YFIT must be 'MostPopular'.> ...
%! setDefaultYfit (compact (TreeBagger (1, x, strcmp (y, 'setosa'))), '')
%!error<CompactTreeBagger.setDefaultYfit: YFIT must be '' or 'MostPopular' for a classification ensemble.> ...
%! setDefaultYfit (C, 'setosa')
%!error<CompactTreeBagger.setDefaultYfit: YFIT must be a real numeric scalar for a regression ensemble.> ...
%! setDefaultYfit (R, 'mean')
%!error<CompactTreeBagger.proximity: too few input arguments.> proximity (C)
%!error<CompactTreeBagger.proximity: X must be a real numeric matrix.> ...
%! proximity (C, {1})
%!error<CompactTreeBagger.proximity: X must have one column per predictor.> ...
%! proximity (C, ones (2, 3))
%!error<CompactTreeBagger.outlierMeasure: too few input arguments.> ...
%! outlierMeasure (C)
%!error<CompactTreeBagger.outlierMeasure: name-value arguments must be in pairs.> ...
%! outlierMeasure (C, x, 'Data')
%!error<CompactTreeBagger.outlierMeasure: invalid parameter name in optional pair arguments.> ...
%! outlierMeasure (C, x, 'Colors', 'r')
%!error<CompactTreeBagger.outlierMeasure: 'Data' must be 'predictors' or 'proximity'.> ...
%! outlierMeasure (C, x, 'Data', 'distance')
%!error<CompactTreeBagger.outlierMeasure: X must be a square real numeric matrix when 'Data' is 'proximity'.> ...
%! outlierMeasure (C, x, 'Data', 'proximity')
%!error<CompactTreeBagger.outlierMeasure: 'Labels' must hold one class label per observation.> ...
%! outlierMeasure (C, x, 'Labels', y(1:10))
%!error<CompactTreeBagger.mdsprox: too few input arguments.> mdsprox (C)
%!error<CompactTreeBagger.mdsprox: 'Colors' must be a character vector or a string scalar.> ...
%! mdsprox (C, x(1:10,:), 'Colors', 1)
%!error<CompactTreeBagger.mdsprox: 'MDSCoordinates' must be a vector of two or three positive integers.> ...
%! mdsprox (C, x(1:10,:), 'MDSCoordinates', 1)
%!error<CompactTreeBagger.mdsprox: 'MDSCoordinates' must not exceed the number of scaled coordinates.> ...
%! mdsprox (C, x(1:10,:), 'Colors', 'r', 'MDSCoordinates', [1, 200])
%!error<CompactTreeBagger.mdsprox: 'MDSCoordinates' must not exceed the number of scaled coordinates.> ...
%! mdsprox (C, x(1:10,:), 'MDSCoordinates', [1, 200])
%!error<CompactTreeBagger.outlierMeasure: 'Labels' must hold only classes the ensemble was trained on.> ...
%! outlierMeasure (C, x(1:2,:), 'Labels', {'rose'; 'setosa'})
%!error<CompactTreeBagger.outlierMeasure: 'Labels' must hold only classes the ensemble was trained on.> ...
%! outlierMeasure (C, x(1:2,:), 'Labels', [1; 2])
%!error<CompactTreeBagger.outlierMeasure: 'Labels' cannot be used with a regression ensemble.> ...
%! outlierMeasure (R, x(1:2,2:4), 'Labels', [1; 2])
%!error<CompactTreeBagger.mdsprox: 'Labels' must hold only classes the ensemble was trained on.> ...
%! mdsprox (C, x(1:2,:), 'Labels', {''; 'setosa'})

## A table at prediction
%!test  # the levels travel with the model, and predict reads a table by name
%! load fisheriris
%! T = table (meas(:,1), meas(:,2), 'VariableNames', {'SL', 'SW'});
%! T.Wide = categorical (meas(:,2) > 3, [false true], {'narrow', 'wide'});
%! T.Species = categorical (species);
%! B = TreeBagger (20, T, 'Species');
%! CB = compact (B);
%! assert_equal (CB.PredictorLevels, B.PredictorLevels);
%! assert_equal (predict (CB, T(:, [4, 3, 2, 1])), predict (CB, T));

## A table at margin
%!test  # the response is named, left out, or given beside the table
%! load fisheriris
%! X = meas(:,1:2);
%! y = categorical (species);
%! T = table (X(:,1), X(:,2), 'VariableNames', {'SL', 'SW'});
%! T.Species = y;
%! Mdl = compact (TreeBagger (20, T, 'Species'));
%! a = margin (Mdl, X, y);
%! assert_equal (margin (Mdl, T(:,1:2), y), a);
%! assert_equal (margin (Mdl, T, 'Species'), a);
%! assert_equal (margin (Mdl, T), a);

## A table at error and meanMargin
%!test  # the response is named, left out, or given beside the table
%! load fisheriris
%! X = meas(:,1:2);
%! y = categorical (species);
%! T = table (X(:,1), X(:,2), 'VariableNames', {'SL', 'SW'});
%! T.Species = y;
%! Mdl = compact (TreeBagger (20, T, 'Species'));
%! a = error (Mdl, X, y);
%! assert_equal (error (Mdl, T(:,1:2), y), a);
%! assert_equal (error (Mdl, T, 'Species'), a);
%! assert_equal (error (Mdl, T), a);
%!test  # the response is named, left out, or given beside the table
%! load fisheriris
%! X = meas(:,1:2);
%! y = categorical (species);
%! T = table (X(:,1), X(:,2), 'VariableNames', {'SL', 'SW'});
%! T.Species = y;
%! Mdl = compact (TreeBagger (20, T, 'Species'));
%! a = meanMargin (Mdl, X, y);
%! assert_equal (meanMargin (Mdl, T(:,1:2), y), a);
%! assert_equal (meanMargin (Mdl, T, 'Species'), a);
%! assert_equal (meanMargin (Mdl, T), a);

## A table at proximity, outlierMeasure and mdsprox
%!test  # a table is matched to the predictors by name
%! load fisheriris
%! T = table (meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!           'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! T.Species = categorical (species);
%! Mdl = compact (TreeBagger (20, T, 'Species'));
%! a = proximity (Mdl, meas);
%! assert_equal (proximity (Mdl, T(:,1:4)), a);
%! assert_equal (proximity (Mdl, T(:,[5, 4, 2, 3, 1])), a);
%!test  # a table is matched to the predictors by name
%! load fisheriris
%! T = table (meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!           'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! T.Species = categorical (species);
%! Mdl = compact (TreeBagger (20, T, 'Species'));
%! a = outlierMeasure (Mdl, meas, 'Labels', species);
%! assert_equal (outlierMeasure (Mdl, T(:,1:4), 'Labels', species), a);
%! assert_equal (outlierMeasure (Mdl, T(:,[5, 4, 2, 3, 1]), ...
%!                               'Labels', species), a);
%!test  # a table is matched to the predictors by name
%! load fisheriris
%! T = table (meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!           'VariableNames', {'SL', 'SW', 'PL', 'PW'});
%! T.Species = categorical (species);
%! Mdl = compact (TreeBagger (20, T, 'Species'));
%! a = mdsprox (Mdl, meas);
%! assert_equal (mdsprox (Mdl, T(:,1:4)), a);
%! assert_equal (mdsprox (Mdl, T(:,[5, 4, 2, 3, 1])), a);
