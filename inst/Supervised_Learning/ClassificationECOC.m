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

classdef ClassificationECOC
## -*- texinfo -*-
## @deftypefn {statistics} ClassificationECOC
##
## A multiclass model built from binary learners.
##
## An error correcting output codes model turns a problem of @math{K} classes
## into a set of two class problems.  A coding matrix gives one column per
## binary learner saying which classes that learner calls +1, which it calls
## -1, and which sit it out; a new observation is sent to every learner and
## given the class whose column of the matrix its scores match most closely.
##
## The fit is carried out by the learners themselves, whichever
## @code{fitcecoc} was asked for, and the decoding by
## @code{CompactClassificationECOC}, which this class holds the data of a fit
## on top of.
##
## @seealso{fitcecoc, CompactClassificationECOC, designecoc}
## @end deftypefn

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} X
    ##
    ## The predictor data the model was fitted on, one row per observation.
    ## This property is read-only.
    ##
    ## @end deftp
    X                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} Y
    ##
    ## The class labels the model was fitted on.  This property is read-only.
    ##
    ## @end deftp
    Y                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} W
    ##
    ## The observation weights, scaled so that each class carries its prior
    ## and the whole sums to one.  This property is read-only.
    ##
    ## @end deftp
    W                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} RowsUsed
    ##
    ## A logical column marking the rows of the original data that were
    ## used, the rest having been dropped as missing.  This property is
    ## read-only.
    ##
    ## @end deftp
    RowsUsed              = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} NumObservations
    ##
    ## The number of observations used.  This property is read-only.
    ##
    ## @end deftp
    NumObservations       = 0;

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} BinaryY
    ##
    ## What each observation was to each binary learner
    ##
    ## An @math{NxL} matrix of -1, 0 and +1, row @math{n} being the row of
    ## @code{CodingMatrix} belonging to that observation's class.  This
    ## property is read-only.
    ##
    ## @end deftp
    BinaryY               = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} BinaryLearners
    ##
    ## The trained binary learners, one per column of @code{CodingMatrix}.
    ## This property is read-only.
    ##
    ## @end deftp
    BinaryLearners        = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} CodingMatrix
    ##
    ## The coding design, a @math{KxL} matrix of -1, 0 and +1.  This
    ## property is read-only.
    ##
    ## @end deftp
    CodingMatrix          = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} CodingName
    ##
    ## The name of the coding design, or @qcode{'custom'} when the matrix was
    ## given outright.  This property is read-only.
    ##
    ## @end deftp
    CodingName            = 'onevsone';

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} LearnerWeights
    ##
    ## The total observation weight each binary learner was trained on.  This
    ## property is read-only.
    ##
    ## @end deftp
    LearnerWeights        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} ClassNames
    ##
    ## The distinct class labels, sorted or in the order given by the
    ## @qcode{'ClassNames'} option; the rows of @code{CodingMatrix},
    ## @code{Prior} and @code{Cost} take them in that order.  This property is
    ## read-only.
    ##
    ## @end deftp
    ClassNames            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} PredictorNames
    ##
    ## The names of the predictors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} ExpandedPredictorNames
    ##
    ## The names of the predictors as the learners saw them.  This property
    ## is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} CategoricalPredictors
    ##
    ## The columns holding categorical predictors, empty when none is.  Every
    ## binary learner was given them.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} ResponseName
    ##
    ## The name of the response variable.  This property is read-only.
    ##
    ## @end deftp
    ResponseName          = 'Y';

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} BinEdges
    ##
    ## The bin edges of the predictors, empty unless the learners binned
    ## them, which none here does.  This property is read-only.
    ##
    ## @end deftp
    BinEdges              = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} HyperparameterOptimizationResults
    ##
    ## The result of optimizing the hyperparameters, always empty here, no
    ## such optimization being implemented.  This property is read-only.
    ##
    ## @end deftp
    HyperparameterOptimizationResults = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} ModelParameters
    ##
    ## A structure of the options the model was fitted with.  This property
    ## is read-only.
    ##
    ## @end deftp
    ModelParameters       = [];

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} BinaryLoss
    ##
    ## The loss that turns a binary learner's score into a cost.  See
    ## @code{CompactClassificationECOC.BinaryLoss}.
    ##
    ## @end deftp
    BinaryLoss            = 'hinge';

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} Prior
    ##
    ## The prior probability of each class.
    ##
    ## @end deftp
    Prior                 = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} Cost
    ##
    ## The cost of misclassification, a @math{KxK} matrix.
    ##
    ## @end deftp
    Cost                  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationECOC} {property} ScoreTransform
    ##
    ## The transform applied to the predicted scores.
    ##
    ## @end deftp
    ScoreTransform        = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    ScoreRange            = [-Inf, Inf];
    STfun                 = @(x) x;
  endproperties

  methods (Hidden)

    function this = set.ScoreTransform (this, val)
      name = 'ClassificationECOC';
      try
        [this.STfun, this.ScoreTransform] = parseScoreTransform (val, name);
      catch
        error (strcat ("ClassificationECOC.subsasgn: 'ScoreTransform' must", ...
                       " be a character vector or a 'function_handle'", ...
                       " object."));
      end_try_catch
    endfunction

    function this = set.BinaryLoss (this, val)
      if (! (ischar (val) && isrow (val)))
        error (strcat ("ClassificationECOC.subsasgn: 'BinaryLoss' must be", ...
                       " a character vector."));
      endif
      [~, errmsg] = ecocDecode (zeros (1, columns (this.CodingMatrix)), ...
                                this.CodingMatrix, tolower (val), ...
                                'lossweighted', this.ScoreRange);
      if (! isempty (errmsg))
        error ("ClassificationECOC.subsasgn: %s", errmsg);
      endif
      this.BinaryLoss = tolower (val);
    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ("\n  ClassificationECOC\n\n");
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %s\n", 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ("%+25s: %s\n", 'ClassNames', classNameListing (this.ClassNames));
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("%+25s: {%dx%d cell}\n", 'BinaryLearners', ...
               numel (this.BinaryLearners), 1);
      fprintf ("%+25s: '%s'\n", 'CodingName', this.CodingName);
      fprintf ("\n");
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationECOC} {@var{obj} =} ClassificationECOC (@var{X}, @var{Y})
    ## @deftypefnx {ClassificationECOC} {@var{obj} =} ClassificationECOC (@dots{}, @var{name}, @var{value})
    ##
    ## Fit a multiclass model from binary learners.
    ##
    ## @code{@var{obj} = ClassificationECOC (@var{X}, @var{Y})} fits one
    ## binary learner per column of a one against one coding design and
    ## returns them as a @code{ClassificationECOC} object.  @code{fitcecoc} is
    ## the documented way in, and its help lists the options both take.
    ##
    ## @seealso{fitcecoc, CompactClassificationECOC, designecoc}
    ## @end deftypefn
    function this = ClassificationECOC (X, Y, varargin)

      if (nargin < 2)
        error ("ClassificationECOC: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error ("ClassificationECOC: name-value arguments must be in pairs.");
      endif
      if (rows (X) != rows (Y))
        error ("ClassificationECOC: number of rows in X and Y must be equal.");
      endif

      ClassNames = []; Cost = []; Prior = []; Weights = [];
      PredictorNames = {}; ResponseName = 'Y'; ScoreTransform = 'none';
      Coding = 'onevsone'; Learners = 'svm'; BinaryLoss = [];
      CatPreds = [];

      for i = 1:2:numel (varargin)
        switch (tolower (varargin{i}))
          case 'classnames'
            ClassNames = varargin{i+1};
          case 'cost'
            Cost = varargin{i+1};
          case 'prior'
            Prior = varargin{i+1};
          case 'weights'
            Weights = varargin{i+1};
          case 'predictornames'
            PredictorNames = varargin{i+1};
          case 'responsename'
            ResponseName = varargin{i+1};
          case 'scoretransform'
            ScoreTransform = varargin{i+1};
          case 'coding'
            Coding = varargin{i+1};
          case 'learners'
            Learners = varargin{i+1};
          case 'binaryloss'
            BinaryLoss = varargin{i+1};
          case 'categoricalpredictors'
            CatPreds = varargin{i+1};
          case 'fitposterior'
            error (strcat ("ClassificationECOC: 'FitPosterior' is not", ...
                           " implemented, the binary learners having no", ...
                           " fitted score transform to install."));
          otherwise
            error (strcat ("ClassificationECOC: invalid parameter name in", ...
                           " optional pair arguments."));
        endswitch
      endfor

      ## The learner, as a template or as a name.
      [tmpl, errmsg] = ClassificationECOC.ecocLearnerTemplate (Learners);
      if (! isempty (errmsg))
        error ("ClassificationECOC: %s", errmsg);
      endif

      F = classFrame (X, Y, ClassNames, Prior, Cost, Weights, ...
                      'ClassificationECOC', false);
      K = classCount (F.ClassNames);

      ## The coding design, named or given outright.
      if (ischar (Coding) && isrow (Coding))
        M = designecoc (K, Coding);
        this.CodingName = tolower (Coding);
      elseif (isnumeric (Coding) && ismatrix (Coding))
        [M, errmsg] = ClassificationECOC.checkCodingMatrix (Coding, K);
        if (! isempty (errmsg))
          error ("ClassificationECOC: %s", errmsg);
        endif
        this.CodingName = 'custom';
      else
        error (strcat ("ClassificationECOC: 'Coding' must be a character", ...
                       " vector or a coding matrix."));
      endif

      this.X                     = F.X;
      this.Y                     = F.Y;
      this.W                     = F.W;
      this.RowsUsed              = F.RowsUsed;
      this.NumObservations       = F.n;
      this.ClassNames            = F.ClassNames;
      this.Prior                 = F.Prior;
      this.Cost                  = F.Cost;
      this.ResponseName          = ResponseName;
      this.CodingMatrix          = M;
      this.BinaryY               = M(F.gY, :);
      this.ScoreRange            = ...
                            ClassificationECOC.ecocScoreRange (tmpl.Method);

      if (isempty (PredictorNames))
        PredictorNames = arrayfun (@(k) sprintf ('x%d', k), ...
                                   1:columns (F.X), 'UniformOutput', false);
      endif
      this.PredictorNames         = PredictorNames;
      this.ExpandedPredictorNames = PredictorNames;

      ## Categorical predictors are checked here and handed to every binary
      ## learner as given, which codes them its own way: the model keeps the
      ## predictor names, the learners report the expanded ones, as in R2024a.
      [Cod, errmsg] = dummyCoding (F.X, CatPreds, PredictorNames);
      if (! isempty (errmsg))
        error ("ClassificationECOC: %s", errmsg);
      endif
      this.CategoricalPredictors = [];
      if (! isempty (Cod.Index))
        this.CategoricalPredictors = Cod.Index;
      endif

      ## One learner per column: the classes that column marks +1 against
      ## those it marks -1, the rest of the rows left out of the fit.
      ## Whether the weights given vary within a class.  Weights that only
      ## carry the prior are constant within each class, and a learner taking
      ## no observation weights can be given that prior instead.
      evenWithin = true;
      for k = 1:classCount (F.ClassNames)
        wk = F.Weights(F.gY == k);
        if (! isempty (wk) && max (wk) - min (wk) > 1e-12 * max (wk))
          evenWithin = false;
        endif
      endfor

      L = columns (M);
      this.BinaryLearners = cell (L, 1);
      this.LearnerWeights = zeros (1, L);
      for j = 1:L
        take = M(F.gY, j) != 0;
        by = M(F.gY(take), j);
        this.BinaryLearners{j} = ...
          ClassificationECOC.ecocFitBinary (tmpl, F.X(take,:), by, ...
                                            F.W(take), PredictorNames, ...
                                            CatPreds, evenWithin);
        this.LearnerWeights(j) = sum (F.W(take));
      endfor

      if (isempty (BinaryLoss))
        BinaryLoss = ClassificationECOC.ecocDefaultLoss (tmpl.Method, ...
                                                         this.ScoreRange);
      endif
      this.BinaryLoss     = BinaryLoss;
      this.ScoreTransform = ScoreTransform;

      this.ModelParameters = struct ('BinaryLearners', tmpl, ...
                                     'Coding', this.CodingName, ...
                                     'FitPosterior', false, ...
                                     'Method', 'ECOC', ...
                                     'Type', 'classification');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationECOC} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationECOC} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross validate a @code{ClassificationECOC}.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} partitions the training data
    ## into ten folds, fits the same model to each fold's training part, and
    ## returns them as a @code{ClassificationPartitionedECOC}.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'KFold'} @tab @tab The number of folds, an integer greater
    ## than one.  The default is ten, or the number of observations when there
    ## are fewer than ten.
    ##
    ## @item @qcode{'Holdout'} @tab @tab The share of the data held out, a
    ## number strictly between 0 and 1.
    ##
    ## @item @qcode{'Leaveout'} @tab @tab @qcode{'on'} for one fold per
    ## observation, @qcode{'off'} otherwise.
    ##
    ## @item @qcode{'CVPartition'} @tab @tab A @code{cvpartition} object,
    ## which names the folds outright.
    ## @end multitable
    ##
    ## Only one of the four may be given.
    ##
    ## @seealso{ClassificationECOC, ClassificationPartitionedECOC, cvpartition}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)

      if (numel (varargin) == 1)
        error (strcat ("ClassificationECOC.crossval: name-value arguments", ...
                       " must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("ClassificationECOC.crossval: specify only one of", ...
                       " the optional name-value paired arguments."));
      endif

      if (this.NumObservations < 10)
        numFolds = this.NumObservations;
      else
        numFolds = 10;
      endif
      Holdout = [];
      Leaveout = 'off';
      CVPartition = [];

      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))
          case 'kfold'
            numFolds = varargin{2};
            if (! (isnumeric (numFolds) && isscalar (numFolds)
                   && numFolds == fix (numFolds) && numFolds > 1))
              error (strcat ("ClassificationECOC.crossval: 'KFold' must", ...
                             " be an integer value greater than 1."));
            endif
          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("ClassificationECOC.crossval: 'Holdout' must", ...
                             " be a numeric value between 0 and 1."));
            endif
          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && any (strcmpi (Leaveout, {'on', 'off'}))))
              error (strcat ("ClassificationECOC.crossval: 'Leaveout' must", ...
                             " be either 'on' or 'off'."));
            endif
          case 'cvpartition'
            CVPartition = varargin{2};
            if (! isa (CVPartition, 'cvpartition'))
              error (strcat ("ClassificationECOC.crossval: 'CVPartition'", ...
                             " must be a 'cvpartition' object."));
            endif
          otherwise
            error (strcat ("ClassificationECOC.crossval: invalid parameter", ...
                           " name in optional paired arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      ## The response is passed rather than a count so that the folds stay
      ## stratified, as every other learner here does.
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (this.Y, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (this.NumObservations, 'LeaveOut');
      else
        partition = cvpartition (this.Y, 'KFold', numFolds);
      endif

      CVMdl = ClassificationPartitionedECOC (this, partition);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationECOC} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Drop the training data from a @code{ClassificationECOC}.
    ##
    ## @code{@var{CMdl} = compact (@var{obj})} returns a
    ## @code{CompactClassificationECOC} carrying the binary learners and the
    ## coding matrix but not @code{X}, @code{Y} or @code{W}, so it predicts
    ## and scores new data but cannot be refitted or cross validated.
    ##
    ## @seealso{ClassificationECOC, CompactClassificationECOC}
    ## @end deftypefn
    function CMdl = compact (this)
      CMdl = CompactClassificationECOC (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationECOC} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationECOC} {[@var{label}, @var{NegLoss}, @var{PBScore}] =} predict (@dots{})
    ## @deftypefnx {ClassificationECOC} {[@dots{}] =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Classify new data with a trained @code{ClassificationECOC}.
    ##
    ## It takes and returns exactly what
    ## @code{CompactClassificationECOC.predict} does, the training data
    ## playing no part in a prediction.
    ##
    ## @seealso{ClassificationECOC, CompactClassificationECOC.predict}
    ## @end deftypefn
    function [label, NegLoss, PBScore] = predict (this, XC, varargin)
      if (nargin < 2)
        error ("ClassificationECOC.predict: too few input arguments.");
      endif
      [label, NegLoss, PBScore] = predict (compact (this), XC, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationECOC} {@var{obj} =} discardSupportVectors (@var{obj})
    ##
    ## Give up the support vectors of the binary learners.
    ##
    ## @code{@var{obj} = discardSupportVectors (@var{obj})} empties
    ## @code{Alpha}, @code{SupportVectors} and @code{SupportVectorLabels} on
    ## every binary learner that is a support vector machine on a linear
    ## kernel, whose linear model stands in for them exactly, so nothing the
    ## model answers changes.  Any other learner is left as it is, a code
    ## being free to mix them, and a model with no linear support vector
    ## machine among its learners warns and is returned unchanged.
    ##
    ## @seealso{ClassificationECOC, ClassificationSVM.discardSupportVectors}
    ## @end deftypefn
    function this = discardSupportVectors (this)

      [learners, n] = ecocDiscardSVs (this.BinaryLearners);
      if (n == 0)
        warning (strcat ("ClassificationECOC.discardSupportVectors: no", ...
                         " binary", ...
                         " learner is a support vector machine on a", ...
                         " linear kernel; nothing was discarded."));
        return;
      endif
      this.BinaryLearners = learners;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationECOC} {@var{sub} =} selectModels (@var{obj}, @var{idx})
    ##
    ## Keep a subset of the fitted regularization strengths.
    ##
    ## @code{@var{sub} = selectModels (@var{obj}, @var{idx})} narrows every
    ## binary learner to the strengths @var{idx} names, which may be indices
    ## into the learner's @qcode{Lambda} or a logical vector over it.  Only a
    ## linear learner is fitted over several strengths, so any other raises.
    ##
    ## @seealso{ClassificationECOC, ClassificationLinear.selectModels}
    ## @end deftypefn
    function sub = selectModels (this, idx)

      if (nargin < 2)
        error ("ClassificationECOC.selectModels: too few input arguments.");
      endif
      [learners, errmsg] = ecocSelectModels (this.BinaryLearners, idx);
      if (! isempty (errmsg))
        error ("ClassificationECOC.selectModels: %s", errmsg);
      endif
      sub = this;
      sub.BinaryLearners = learners;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationECOC} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin of a @code{ClassificationECOC}.  See
    ## @code{CompactClassificationECOC.margin}.
    ##
    ## @seealso{CompactClassificationECOC.margin}
    ## @end deftypefn
    function m = margin (this, X, Y, varargin)
      if (nargin < 3)
        error ("ClassificationECOC.margin: too few input arguments.");
      endif
      m = margin (compact (this), X, Y, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationECOC} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification edge of a @code{ClassificationECOC}.  See
    ## @code{CompactClassificationECOC.edge}.
    ##
    ## @seealso{CompactClassificationECOC.edge}
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)
      if (nargin < 3)
        error ("ClassificationECOC.edge: too few input arguments.");
      endif
      e = edge (compact (this), X, Y, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationECOC} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification loss of a @code{ClassificationECOC}.  See
    ## @code{CompactClassificationECOC.loss}.
    ##
    ## @seealso{CompactClassificationECOC.loss}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)
      if (nargin < 3)
        error ("ClassificationECOC.loss: too few input arguments.");
      endif
      L = loss (compact (this), X, Y, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationECOC} {@var{label} =} resubPredict (@var{obj})
    ##
    ## Classify the training data with a trained @code{ClassificationECOC}.
    ##
    ## @seealso{ClassificationECOC.predict}
    ## @end deftypefn
    function [label, NegLoss, PBScore] = resubPredict (this, varargin)
      [label, NegLoss, PBScore] = predict (this, this.X, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationECOC} {@var{m} =} resubMargin (@var{obj})
    ##
    ## Classification margin on the training data.
    ##
    ## @seealso{ClassificationECOC.margin}
    ## @end deftypefn
    function m = resubMargin (this, varargin)
      m = margin (this, this.X, this.Y, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationECOC} {@var{e} =} resubEdge (@var{obj})
    ##
    ## Classification edge on the training data.
    ##
    ## @seealso{ClassificationECOC.edge}
    ## @end deftypefn
    function e = resubEdge (this, varargin)
      e = edge (this, this.X, this.Y, 'Weights', this.W, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationECOC} {@var{L} =} resubLoss (@var{obj})
    ##
    ## Classification loss on the training data.
    ##
    ## @seealso{ClassificationECOC.loss}
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      L = loss (this, this.X, this.Y, 'Weights', this.W, varargin{:});
    endfunction

  endmethods

  ## Helpers of the fit.  Each has one consumer, so it sits here rather than
  ## in private/, and none of them needs the object.
  methods (Static, Access = private)

    ## The learner, given as a template or by name.
    function [tmpl, errmsg] = ecocLearnerTemplate (Learners)

      tmpl = [];
      errmsg = '';
      known = {'svm', 'tree', 'knn', 'naivebayes', 'discriminant', ...
               'linear', 'kernel'};
      makers = {@templateSVM, @templateTree, @templateKNN, ...
                @templateNaiveBayes, @templateDiscriminant, ...
                @templateLinear, @templateKernel};
      if (ischar (Learners) && isrow (Learners))
        k = find (strcmp (tolower (Learners), known));
        if (strcmpi (Learners, 'ensemble'))
          ## MATLAB's default ensemble for a binary learner.
          tmpl = templateEnsemble ('LogitBoost', 100, 'tree');
        elseif (isempty (k))
          errmsg = strcat ("'", Learners, "' is not a binary learner.");
          return;
        else
          tmpl = makers{k} ();
        endif
      elseif (isstruct (Learners) && isscalar (Learners)
              && isfield (Learners, 'Method'))
        tmpl = Learners;
        if (any (strcmp (tmpl.Method, ClassificationECOC.ensembleMethods ())))
          if (! strcmp (tmpl.Type, 'classification'))
            errmsg = "templates of regression type are not supported.";
            return;
          endif
        elseif (! any (strcmp (tolower (tmpl.Method), known)))
          errmsg = strcat ("'", tmpl.Method, "' is not a binary learner.");
          return;
        endif
      else
        errmsg = strcat ("'Learners' must be a learner name or a template.");
      endif

    endfunction

    ## The interval a learner scores on, which is what decides the losses its
    ## scores can be read with.  Measured on R2024a across all seven learners
    ## and the ensembles: a boosted ensemble scores on the whole line, and a
    ## bagged or random subspace one with class probabilities.  A RUSBoost
    ## one scores in [0,1] once CompactClassificationECOC scales it.
    function range = ecocScoreRange (method)

      if (any (strcmpi (method, {'SVM', 'Linear', 'Kernel', 'AdaBoostM1', ...
                                 'AdaBoostM2', 'GentleBoost', ...
                                 'LogitBoost', 'LPBoost', 'TotalBoost', ...
                                 'RobustBoost'})))
        range = [-Inf, Inf];
      else
        range = [0, 1];
      endif

    endfunction

    ## The loss a learner's scores are read with when none is asked for, as
    ## R2024a chooses it.
    function loss = ecocDefaultLoss (method, range)

      if (any (strcmp (method, {'AdaBoostM1', 'GentleBoost'})))
        loss = 'exponential';
      elseif (strcmp (method, 'LogitBoost'))
        loss = 'binodeviance';
      elseif (isequal (range, [0, 1]))
        loss = 'quadratic';
      else
        loss = 'hinge';
      endif

    endfunction

    ## The methods an ensemble template may name.
    function names = ensembleMethods ()

      names = {'AdaBoostM1', 'AdaBoostM2', 'GentleBoost', 'LogitBoost', ...
               'LPBoost', 'TotalBoost', 'RobustBoost', 'RUSBoost', ...
               'Subspace', 'Bag', 'LSBoost'};

    endfunction

    ## Fit one binary learner.  Its two classes are given outright as -1 and
    ## +1 so that the second is always the one the column calls +1, which is
    ## the score the decoding reads.
    function Mdl = ecocFitBinary (tmpl, X, y, w, pnames, cats = [], ...
                                  evenWithin = true)

      ## An ensemble template carries the method, cycles and learners under
      ## names of its own, which the ensemble takes under fitcensemble's.
      ensemble = any (strcmp (tmpl.Method, ...
                              ClassificationECOC.ensembleMethods ()));
      skip = {'Method', 'Type'};
      args = {};
      if (ensemble)
        skip(end+1:end+2) = {'NLearn', 'LearnerTemplates'};
        args = {'Method', tmpl.Method, 'NumLearningCycles', tmpl.NLearn, ...
                'Learners', tmpl.LearnerTemplates};
      endif
      for [val, name] = tmpl
        if (any (strcmp (name, skip)))
          continue;
        endif
        args(end+1:end+2) = {name, val};
      endfor
      ## Weights are passed only when they carry information.  Uniform ones
      ## say nothing a learner does not assume, and three of the seven
      ## learners take no 'Weights' at all, so passing them regardless would
      ## refuse the commonest fit there is.  A learner that cannot take them
      ## refuses under its own name, which is the right place for it.
      args(end+1:end+4) = {'PredictorNames', pnames, 'ClassNames', [-1; 1]};
      if (! isempty (cats))
        args(end+1:end+2) = {'CategoricalPredictors', cats};
      endif
      ## Weights spread a prior over its classes one class at a time, so even
      ## equal ones differ in the last bits; only a real spread counts.  Four
      ## of the learners take no observation weights.  They are given the
      ## share of the weight each side holds as their prior, which is the
      ## prior MATLAB's learners report, and weights that vary within a class
      ## have no such form and are refused.
      if (max (w) - min (w) > 1e-12 * max (w))
        if (any (strcmpi (tmpl.Method, {'svm', 'knn', 'naivebayes', ...
                                        'discriminant'})))
          if (! evenWithin)
            error (strcat ("ClassificationECOC: the '%s' learners take no", ...
                           " observation weights, so 'Weights' that vary", ...
                           " within a class cannot be used with them."), ...
                   tolower (tmpl.Method));
          endif
          args(end+1:end+2) = {'Prior', [sum(w(y == -1)), sum(w(y == 1))]};
        else
          args(end+1:end+2) = {'Weights', w};
        endif
      endif

      if (ensemble && strcmp (tmpl.Method, 'Bag'))
        Mdl = ClassificationBaggedEnsemble (X, y, args{:});
        return;
      elseif (ensemble)
        Mdl = ClassificationEnsemble (X, y, args{:});
        return;
      endif
      switch (tolower (tmpl.Method))
        case 'svm'
          Mdl = ClassificationSVM (X, y, args{:});
        case 'tree'
          Mdl = ClassificationTree (X, y, args{:});
        case 'knn'
          Mdl = ClassificationKNN (X, y, args{:});
        case 'naivebayes'
          Mdl = ClassificationNaiveBayes (X, y, args{:});
        case 'discriminant'
          Mdl = ClassificationDiscriminant (X, y, args{:});
        case 'linear'
          Mdl = ClassificationLinear (X, y, args{:});
        case 'kernel'
          Mdl = ClassificationKernel (X, y, args{:});
      endswitch

    endfunction

    ## A coding matrix given outright.  A column that marks nothing +1 or
    ## nothing -1 trains no learner, and one that repeats another up to sign
    ## trains the same learner twice; both are refused, as R2024a refuses
    ## them.
    function [M, errmsg] = checkCodingMatrix (M, K)

      errmsg = '';
      if (! (isnumeric (M) && ismatrix (M) && ! isempty (M)
             && rows (M) == K))
        errmsg = strcat ("'Coding' must be a matrix with one row per", ...
                         " class.");
        return;
      endif
      if (! all (ismember (M(:), [-1, 0, 1])))
        errmsg = "'Coding' must hold only -1, 0 and 1.";
        return;
      endif
      if (! all (any (M > 0, 1) & any (M < 0, 1)))
        errmsg = strcat ("every column of 'Coding' must mark at least one", ...
                         " class 1 and at least one -1.");
        return;
      endif
      if (any (all (M == 0, 2)))
        errmsg = "every class must take part in at least one column.";
        return;
      endif
      L = columns (M);
      for a = 1:L-1
        for b = a+1:L
          if (isequal (M(:,a), M(:,b)) || isequal (M(:,a), -M(:,b)))
            errmsg = sprintf (strcat ("columns %d and %d in the coding", ...
                                      " matrix are identical or differ", ...
                                      " only by sign."), a, b);
            return;
          endif
        endfor
      endfor

    endfunction

  endmethods

endclassdef

## Tests
%!test  # MATLAB parity: the property surface a fit reports
%! load fisheriris
%! Mdl = ClassificationECOC (meas, species);
%! assert_equal (class (Mdl), 'ClassificationECOC');
%! assert_equal (numel (properties (Mdl)), 22);
%! assert_equal (Mdl.ResponseName, 'Y');
%! assert_equal (Mdl.ClassNames, unique (species));
%! assert_equal (Mdl.Prior, [1/3, 1/3, 1/3], 1e-14);
%! assert_equal (Mdl.Cost, ones (3) - eye (3));
%! assert_equal (Mdl.ScoreTransform, 'none');
%! assert_equal (Mdl.CategoricalPredictors, []);

%!test  # the binary learners are one per column of the coding matrix
%! load fisheriris
%! Mdl = ClassificationECOC (meas, species, 'Coding', 'ternarycomplete');
%! assert_equal (columns (Mdl.CodingMatrix), 6);
%! assert_equal (numel (Mdl.BinaryLearners), 6);
%! assert_equal (numel (Mdl.LearnerWeights), 6);

%!test  # MATLAB parity: resubstitution answers the training data
%! load fisheriris
%! Mdl = ClassificationECOC (meas, species, 'Learners', 'tree');
%! assert_equal (resubLoss (Mdl), 0.02, 1e-12);
%! label = resubPredict (Mdl);
%! assert_equal (sum (! strcmp (label, species)), 3);

%!test  # a margin is positive exactly where the label was right
%! load fisheriris
%! Mdl = ClassificationECOC (meas, species, 'Learners', 'tree');
%! m = resubMargin (Mdl);
%! right = strcmp (resubPredict (Mdl), species);
%! assert_equal (m > 0, right);
%! assert_equal (sum (! right), 3);

%!test  # compact keeps the learners and drops the data
%! load fisheriris
%! Mdl = ClassificationECOC (meas, species);
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), 'CompactClassificationECOC');
%! assert_equal (CMdl.CodingMatrix, Mdl.CodingMatrix);
%! assert_equal (predict (CMdl, meas(1:5,:)), predict (Mdl, meas(1:5,:)));
%! assert_equal (isprop (CMdl, 'X'), false);

%!test  # a class the coding matrix leaves out of every column is refused
%! load fisheriris
%! M = [1, 1; -1, 0; 0, -1];
%! assert_equal (columns (ClassificationECOC (meas, species, ...
%!                                            'Coding', M).CodingMatrix), 2);

%!test  # crossval returns the ECOC partitioned class, not the general one
%! load fisheriris
%! CV = crossval (ClassificationECOC (meas, species, 'Learners', 'tree'));
%! assert_equal (class (CV), 'ClassificationPartitionedECOC');
%! assert_equal (CV.KFold, 10);
%! assert_equal (CV.CodingMatrix, [1, 1, 0; -1, 0, 1; 0, -1, -1]);

%!test  # discarding support vectors changes nothing the model answers
%! ## The linear model stands in for the vectors exactly, so the labels and
%! ## the loss are what they were and only the memory is gone.
%! load fisheriris
%! Mdl = fitcecoc (meas, species);
%! assert_equal (size (Mdl.BinaryLearners{1}.SupportVectors), [3, 4]);
%! D = discardSupportVectors (Mdl);
%! assert_equal (class (D), 'ClassificationECOC');
%! assert_equal (size (D.BinaryLearners{1}.SupportVectors), [0, 0]);
%! assert_equal (isempty (D.BinaryLearners{1}.Alpha), true);
%! assert_equal (predict (D, meas), predict (Mdl, meas));
%! assert_equal (resubLoss (D), 0.0066666666666667, 1e-12);

%!warning<ClassificationECOC.discardSupportVectors: no binary learner is a support vector machine on a linear kernel; nothing was discarded.> ...
%! load fisheriris; ...
%! discardSupportVectors (fitcecoc (meas, species, 'Learners', 'tree'));

%!test  # a code of trees is returned unchanged by discardSupportVectors
%! load fisheriris
%! Mdl = fitcecoc (meas, species, 'Learners', 'tree');
%! ## The state is saved and put back rather than switched on: 'on' would
%! ## enable warning classes Octave disables by default and leak them into
%! ## every test that runs after this one.
%! w = warning ('off', 'all');
%! D = discardSupportVectors (Mdl);
%! warning (w);
%! assert_equal (predict (D, meas), predict (Mdl, meas));

%!test  # selectModels narrows every binary learner to the same strengths
%! load fisheriris
%! ## The constructor is used and not fitcecoc: a linear learner keeps no
%! ## training data, so the fit route gives back the compact model.
%! LC = ClassificationECOC (meas, species, ...
%!        'Learners', templateLinear ('Lambda', [1e-4, 1e-3, 1e-2]));
%! assert_equal (LC.BinaryLearners{1}.Lambda, [1e-4, 1e-3, 1e-2], 1e-12);
%! S = selectModels (LC, 2);
%! assert_equal (class (S), 'ClassificationECOC');
%! for j = 1:numel (S.BinaryLearners)
%!   assert_equal (S.BinaryLearners{j}.Lambda, 1e-3, 1e-12);
%! endfor

%!error<ClassificationECOC.selectModels: too few input arguments.> ...
%! load fisheriris; ...
%! selectModels (fitcecoc (meas, species))
%!error<ClassificationECOC.selectModels: the binary learners are 'ClassificationTree' models, which are fitted over one regularization strength and have none to select between.> ...
%! load fisheriris; ...
%! selectModels (fitcecoc (meas, species, 'Learners', 'tree'), 1)

%!test  # An unused category of a categorical response is not a class
%! load fisheriris
%! y = categorical (species);
%! Mdl = ClassificationECOC (meas(51:150,:), y(51:150));
%! assert_equal (cellstr (Mdl.ClassNames), {'versicolor'; 'virginica'});
%! assert_equal (size (Mdl.CodingMatrix), [2, 1]);

## Test input validation
%!error<ClassificationECOC.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (ClassificationECOC (ones (4, 2), [1; 2; 1; 2]), 'KFold', 1)
%!error<ClassificationECOC.crossval: specify only one of the optional name-value paired arguments.> ...
%! Mdl = ClassificationECOC (ones (4, 2), [1; 2; 1; 2]); ...
%! crossval (Mdl, 'KFold', 2, 'Holdout', 0.3)
%!error<ClassificationECOC: too few input arguments.> ...
%! ClassificationECOC (ones (4, 2))
%!error<ClassificationECOC: name-value arguments must be in pairs.> ...
%! ClassificationECOC (ones (4, 2), [1; 2; 1; 2], 'Coding')
%!error<ClassificationECOC: number of rows in X and Y must be equal.> ...
%! ClassificationECOC (ones (4, 2), [1; 2; 1])
%!error<ClassificationECOC: 'Coding' must hold only -1, 0 and 1.> ...
%! ClassificationECOC (ones (4, 2), [1; 2; 1; 2], 'Coding', [2; -2])
%!error<ClassificationECOC: every column of 'Coding' must mark at least one class 1 and at least one -1.> ...
%! ClassificationECOC (ones (4, 2), [1; 2; 1; 2], 'Coding', [1; 1])
%!error<ClassificationECOC: columns 1 and 2 in the coding matrix are identical or differ only by sign.> ...
%! ClassificationECOC (ones (4, 2), [1; 2; 1; 2], 'Coding', [1, -1; -1, 1])
%!error<ClassificationECOC: 'Coding' must be a character vector or a coding matrix.> ...
%! ClassificationECOC (ones (4, 2), [1; 2; 1; 2], 'Coding', {1})

%!shared Xc, y3
%! k = (0:119)';
%! c1 = mod (k, 3) + 1;
%! c3 = 10 * (mod (floor (k / 2), 2) + 1);
%! Xc = [c1, sin(k), c3];
%! y3 = mod (c1 + floor (k / 7), 3) + 1;

%!test  # MATLAB parity: every binary learner gets the categorical predictors
%! Mdl = ClassificationECOC (Xc, y3, 'CategoricalPredictors', [1, 3]);
%! assert_equal (Mdl.CategoricalPredictors, [1, 3]);
%! assert_equal (Mdl.ExpandedPredictorNames, {'x1', 'x2', 'x3'});
%! L = Mdl.BinaryLearners{1};
%! assert_equal (L.CategoricalPredictors, [1, 3]);
%! assert_equal (L.ExpandedPredictorNames, {'x1 == 1', 'x1 == 2', ...
%!               'x1 == 3', 'x2', 'x3 == 10', 'x3 == 20'});
%! assert_equal (compact (Mdl).CategoricalPredictors, [1, 3]);

%!test  # tree learners take them too, and so do the folds
%! Mdl = ClassificationECOC (Xc, y3, 'CategoricalPredictors', [1, 3], ...
%!                           'Learners', 'tree');
%! assert_equal (Mdl.BinaryLearners{1}.CategoricalPredictors, [1, 3]);
%! CV = crossval (Mdl, 'KFold', 3);
%! assert_equal (CV.Trained{1}.CategoricalPredictors, [1, 3]);

%!test  # MATLAB parity: nearest neighbour learners compare levels
%! Mdl = ClassificationECOC (Xc(:,[1, 3]), y3, 'CategoricalPredictors', ...
%!                           'all', 'Learners', 'knn');
%! assert_equal (Mdl.BinaryLearners{1}.Distance, 'hamming');

%!error<ClassificationECOC: 'CategoricalPredictors' indices must not exceed the number of predictors.> ...
%! ClassificationECOC (Xc, y3, 'CategoricalPredictors', 4)
%!error<ClassificationDiscriminant: categorical predictors cannot be used for discriminant analysis.> ...
%! ClassificationECOC (Xc, y3, 'CategoricalPredictors', 1, ...
%!                     'Learners', 'discriminant')

%!test  # equal weights that differ only by rounding reach any learner
%! load fisheriris
%! y = [ones(60, 1); 2 * ones(50, 1); 3 * ones(40, 1)];
%! Mdl = ClassificationECOC (meas, y, 'Learners', 'knn');
%! assert_equal (numel (Mdl.BinaryLearners), 3);
%! Mdl = ClassificationECOC (meas, y);
%! assert_equal (class (Mdl.BinaryLearners{1}), 'ClassificationSVM');

%!test  # MATLAB parity: a learner without weights reports its side's prior
%! load fisheriris
%! y = [ones(60, 1); 2 * ones(50, 1); 3 * ones(40, 1)];
%! Mdl = ClassificationECOC (meas, y, 'Learners', 'knn', 'Prior', 'uniform');
%! assert_equal (Mdl.BinaryLearners{1}.Prior, [0.5, 0.5], 1e-12);
%! Mdl = ClassificationECOC (meas(11:150,:), y(11:150));
%! assert_equal (Mdl.BinaryLearners{1}.Prior, [50, 50] / 100, 1e-12);

%!error<ClassificationECOC: the 'svm' learners take no observation weights, so 'Weights' that vary within a class cannot be used with them.> ...
%! load fisheriris
%! ClassificationECOC (meas, species, 'Weights', (1:150)')
