## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
## Copyright (C) 2024 Pallav Purbia <pallavpurbia@gmail.com>
## Copyright (C) 2024-2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

classdef ClassificationPartitionedModel
  ## -*- texinfo -*-
  ## @deftp {statistics} ClassificationPartitionedModel
  ##
  ## Cross-validated classification model
  ##
  ## The @code{ClassificationPartitionedModel} class stores cross-validated
  ## classification models trained on different partitions of the data.
  ## It can predict responses for observations not used for training using
  ## the @code{kfoldPredict} method.
  ##
  ## Create a @code{ClassificationPartitionedModel} object by using the
  ## @code{crossval} function.
  ##
  ## @seealso{crossval}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector with one entry per observation, carried over
    ## from the model that was cross validated.  This property is read-only.
    ##
    ## @end deftp
    W                            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} BinEdges
    ##
    ## Bin edges
    ##
    ## A cell array with one entry per predictor, holding that predictor's
    ## bin edges where the learner discretized it before fitting.  It is
    ## carried over from the model that was cross validated, and is empty
    ## whenever that model did no binning, which is every learner this package
    ## implements: MATLAB fills it only for its GAM, which bins because it is
    ## built from boosted trees where ours is built from splines.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    BinEdges                     = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} CategoricalPredictors
    ##
    ## Indices of categorical predictors
    ##
    ## A vector of positive integers specifying the indices of categorical
    ## predictors.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix containing the unstandardized predictor data.  Each
    ## column of @var{X} represents one predictor (variable), and each row
    ## represents one observation.  This property is read-only.
    ##
    ## @end deftp
    X                            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Y
    ##
    ## Class labels
    ##
    ## Specified as a logical or numeric column vector, or as a character array
    ## or a cell array of character vectors with the same number of rows as the
    ## predictor data.  Each row in @var{Y} is the observed class label for
    ## the corresponding row in @var{X}.  This property is read-only.
    ##
    ## @end deftp
    Y                            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} ClassNames
    ##
    ## Names of classes in the response variable
    ##
    ## An array of unique values of the response variable @var{Y}, which has the
    ## same data types as the data in @var{Y}.  This property is read-only.
    ## @qcode{ClassNames} can have any of the following datatypes:
    ##
    ## @itemize
    ## @item Cell array of character vectors
    ## @item Character array
    ## @item Logical vector
    ## @item Numeric vector
    ## @end itemize
    ##
    ## @end deftp
    ClassNames                   = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} CrossValidatedModel
    ##
    ## Cross-validated model class
    ##
    ## A character vector holding the short name of the learner that was
    ## cross validated, as MATLAB reports it: @qcode{'Discriminant'},
    ## @qcode{'GAM'}, @qcode{'KNN'}, @qcode{'NeuralNetwork'} or
    ## @qcode{'SVM'}.  It is not the class name of that learner, and the
    ## regression side uses the same names.  This property is read-only.
    ##
    ## @end deftp
    CrossValidatedModel          = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} KFold
    ##
    ## Number of cross-validated folds
    ##
    ## A positive integer value specifying the number of cross-validated folds.
    ## This property is read-only.
    ##
    ## @end deftp
    KFold                        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} ModelParameters
    ##
    ## Model parameters
    ##
    ## A structure containing the model parameters used during training.
    ## This includes any model-specific parameters that were configured prior
    ## to training.  This property is read-only.
    ##
    ## @end deftp
    ModelParameters              = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer value specifying the number of observations in the
    ## training dataset used for training the cross-validated model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations              = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Partition
    ##
    ## Partition configuration
    ##
    ## A @code{cvpartition} object specifying the partition configuration used
    ## for cross-validation.  This field stores the cvpartition instance that
    ## describes how the data was split into training and validation sets.
    ## This property is read-only.
    ##
    ## @end deftp
    Partition                    = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which they appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames               = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName                 = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Trained
    ##
    ## Models trained on each fold
    ##
    ## A cell array of models trained on each fold.  Each cell contains a model
    ## trained on the minus-one fold of the data (all but one fold used for
    ## training and the remaining fold used for validation).  This property is
    ## read-only.
    ##
    ## @end deftp
    Trained                      = [];
  endproperties

  ## Properties a user may set.  Each one is validated by its set method.
  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Cost
    ##
    ## Cost of Misclassification
    ##
    ## A square matrix specifying the cost of misclassification of a point.
    ## @qcode{Cost(i,j)} is the cost of classifying a point into class @qcode{j}
    ## if its true class is @qcode{i} (that is, the rows correspond to the true
    ## class and the columns correspond to the predicted class).  The order of
    ## the rows and columns in @qcode{Cost} corresponds to the order of the
    ## classes in @qcode{ClassNames}.  The number of rows and columns in
    ## @qcode{Cost} is the number of unique classes in the response.  By
    ## default, @qcode{Cost(i,j) = 1} if @qcode{i != j}, and
    ## @qcode{Cost(i,j) = 0} if @qcode{i = j}.  In other words, the cost is 0
    ## for correct classification and 1 for incorrect classification.
    ##
    ## Assigning @qcode{Cost} rebuilds it on every fold in @qcode{Trained}, so
    ## @code{kfoldPredict} and @code{kfoldLoss} answer under the new costs.  It
    ## is refused on a cross-validated @code{ClassificationSVM}, whose costs
    ## enter the box constraint while it is being fitted: a model already fitted
    ## under one cost matrix cannot be made to describe another.
    ##
    ##
    ## A cost may also be given as a struct with the fields
    ## @qcode{ClassNames} and @qcode{ClassificationCosts}, which names the
    ## order its own matrix is written in.  That matrix is permuted into the
    ## order of @qcode{ClassNames} above, so a caller need not know which
    ## order the classes were sorted into.  It must name every class.
    ##
    ## A cost must be floating point, not sparse, not complex, non-negative
    ## and zero down its diagonal, and must hold no @qcode{NaN} or
    ## @qcode{Inf}.  A @code{single} is widened to @code{double}.
    ## @end deftp
    Cost                         = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} Prior
    ##
    ## Prior probability for each class
    ##
    ## A numeric vector specifying the prior probabilities for each class.  The
    ## order of the elements in @qcode{Prior} corresponds to the order of the
    ## classes in @qcode{ClassNames}.
    ##
    ## It may be assigned only on a cross-validated
    ## @code{ClassificationDiscriminant}, the one learner that re-derives its
    ## coefficients from the priors it is given; every other learner consumes
    ## them while it fits and cannot revisit them afterwards.  Assigning it
    ## rebuilds the priors on every fold in @qcode{Trained}.
    ##
    ## @end deftp
    Prior                        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedModel} {property} ScoreTransform
    ##
    ## Transformation function for classification scores
    ##
    ## Specified as a function handle for transforming the classification
    ## scores.
    ##
    ## @end deftp
    ScoreTransform               = [];

  endproperties

  ## Copied from the parent model and kept out of the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    STfun = @(x) x;

    ## Raised once the constructor is done.  Cost and Prior are refused for
    ## the learners MATLAB refuses them for, but the constructor has to carry
    ## whatever the cross validated model held, so the guards only apply to a
    ## user's assignment.
    Fitted = false;
  endproperties

  ## Set methods for the properties a user may assign.
  methods (Hidden)

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    ## Shaped after RegressionPartitionedModel's, with the classifier's own
    ## fields: the classes it separates and the transform its scores carry.
    function disp (this)
      fprintf ("\n  ClassificationPartitionedModel\n\n");
      fprintf ("%+25s: '%s'\n", 'CrossValidatedModel', ...
               this.CrossValidatedModel);
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      if (iscellstr (this.ClassNames))
        str = repmat ({'''%s'''}, 1, numel (this.ClassNames));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.ClassNames{:});
      elseif (ischar (this.ClassNames))
        str = repmat ({'''%s'''}, 1, rows (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, cellstr (this.ClassNames){:});
      else   # single, double, logical
        str = repmat ({'%d'}, 1, numel (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, this.ClassNames);
      endif
      fprintf ("%+25s: %s\n", 'ClassNames', str);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'KFold', this.KFold);
      fprintf ("%+25s: '%s'\n\n", 'ScoreTransform', this.ScoreTransform);
    endfunction

    ## MATLAB refuses this one on a cross-validated SVM alone, and so do we:
    ## the SVM's costs enter its box constraint while it is being fitted, so a
    ## model fitted under one cost matrix cannot be made to report another.
    function this = set.Cost (this, val)
      if (this.Fitted && strcmp (this.CrossValidatedModel, 'SVM'))
        error (strcat ("ClassificationPartitionedModel: cannot assign", ...
                       " 'Cost' on a cross-validated ClassificationSVM,", ...
                       " whose costs are consumed while it is fitted."));
      endif
      gnY = this.ClassNames;
      if (isempty (val))
        this.Cost = cast (! eye (classCount (gnY)), 'double');
      else
        ## Everything a cost must be, and the struct form, which
        ## is permuted into this model's class order.
        [val, errmsg] = costMatrix (val, gnY);
        if (! isempty (errmsg))
          error ("ClassificationPartitionedModel: %s", errmsg);
        endif
        this.Cost = val;
      endif
      this = pushToFolds (this, 'Cost', this.Cost);
    endfunction

    ## The discriminant is the only learner that re-derives its coefficients
    ## from the priors it holds, so it is the only one whose folds can honour
    ## a prior assigned after the fit.  MATLAB draws the same line.
    function this = set.Prior (this, val)
      if (this.Fitted && ! strcmp (this.CrossValidatedModel, ...
                                   'Discriminant'))
        error (strcat ("ClassificationPartitionedModel: 'Prior' can only", ...
                       " be assigned on a cross-validated", ...
                       " ClassificationDiscriminant."));
      endif
      if (! this.Fitted)
        this.Prior = val;
        return;
      endif
      if (isstruct (val))
        val = priorFromStruct (val, this.ClassNames, ...
                               'ClassificationPartitionedModel');
      endif
      if (ischar (val) && strcmpi ('uniform', val))
        n = classCount (this.ClassNames);
        this.Prior = ones (1, n) ./ n;
      elseif (isempty (val) || (ischar (val) && strcmpi ('empirical', val)))
        [~, gnY, gY] = uniqueLabels (this.Y);
        pr = accumarray (gY(:), 1, [numel(gnY), 1]);
        this.Prior = pr(:)' ./ sum (pr);
      elseif (isnumeric (val))
        if (classCount (this.ClassNames) != numel (val))
          error (strcat ("ClassificationPartitionedModel: the elements", ...
                         " in 'Prior' do not correspond to the selected", ...
                         " classes in Y."));
        endif
        this.Prior = val(:)' ./ sum (val);
      else
        error (strcat ("ClassificationPartitionedModel: invalid value", ...
                       " for 'Prior'."));
      endif
      this = pushToFolds (this, 'Prior', this.Prior);
    endfunction

    function this = set.ScoreTransform (this, val)
      [f, nm] = parseScoreTransform (val, 'ClassificationPartitionedModel');
      this.ScoreTransform = nm;
      this.STfun = f;
    endfunction

  endmethods

  methods (Access = public)
    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedModel} {@var{this} =} ClassificationPartitionedModel (@var{Mdl}, @var{Partition})
    ##
    ## Create a @code{ClassificationPartitionedModel} class object for
    ## cross-validation of classification models.
    ##
    ## @code{@var{this} = ClassificationPartitionedModel (@var{Mdl},
    ## @var{Partition})} returns a ClassificationPartitionedModel object, with
    ## @var{Mdl} as the trained classification model object and
    ## @var{Partition} as the partitioning object obtained using
    ## @code{cvpartition}
    ## function.
    ##
    ## @seealso{cvpartition}
    ## @end deftypefn
    function this = ClassificationPartitionedModel (Mdl, Partition)

      ## Check input arguments
      if (nargin < 2)
        error ("ClassificationPartitionedModel: too few input arguments.");
      endif

      ## Check for valid Classification object
      validTypes = {'ClassificationDiscriminant', 'ClassificationGAM', ...
                    'ClassificationKNN', 'ClassificationNeuralNetwork', ...
                    'ClassificationSVM'};
      if (! any (strcmp (class (Mdl), validTypes)))
        error ("ClassificationPartitionedModel: unsupported model type.");
      endif

      ## Check for valid cvpartition object
      if (! strcmp (class (Partition), 'cvpartition'))
        error ("ClassificationPartitionedModel: invalid 'cvpartition' object.");
      endif

      ## Set properties.  The rows dropped for missing values are outside the
      ## partition, so they are dropped here too and every index below, the
      ## partition's included, refers to the same set of observations.
      this.X = Mdl.X;
      this.Y = Mdl.Y;
      this.KFold = Partition.NumTestSets;
      this.Trained = cell (this.KFold, 1);
      this.ClassNames = Mdl.ClassNames;
      this.ResponseName = Mdl.ResponseName;
      this.NumObservations = rows (this.X);
      this.PredictorNames = Mdl.PredictorNames;
      this.Partition = Partition;
      ## MATLAB stores a short name here, shared with the regression
      ## side: 'Discriminant', 'GAM', 'KNN', 'NeuralNetwork', 'SVM'.
      this.CrossValidatedModel = strrep (class (Mdl), 'Classification', '');
      this.ScoreTransform = Mdl.ScoreTransform;
      this.STfun = Mdl.STfun;
      ## Every classifier reports a prior and a cost now, so they are carried
      ## whatever was cross validated; this used to name three of the five.
      this.Prior = Mdl.Prior;
      this.Cost = Mdl.Cost;
      this.W = Mdl.W;
      this.BinEdges = Mdl.BinEdges;
      ## Switch Classification object types
      switch (this.CrossValidatedModel)

        case 'Discriminant'
          ## Arguments to pass in fitcdiscr
          args = {};
          ## List of acceptable parameters for fitcdiscr
          DiscrParams = {'PredictorNames', 'ResponseName', 'ClassNames', ...
                         'Cost', 'DiscrimType', 'Gamma'};
          ## Set parameters
          for i = 1:numel (DiscrParams)
            paramName = DiscrParams{i};
            paramValue = Mdl.(paramName);
            if (! isempty (paramValue))
              args = [args, {paramName, paramValue}];
            endif
          endfor
          ## Add 'FillCoeffs' parameter
          if (isempty (Mdl.Coeffs))
            args = [args, {'FillCoeffs', 'off'}];
          endif

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            tmp = fitcdiscr (this.X(idx, :), this.Y(idx,:), args{:});
            this.Trained{k} = compact (tmp);
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct ();
          paramList = {'DiscrimType', 'FillCoeffs', 'Gamma'};
          for i = 1:numel (paramList)
            paramName = paramList{i};
            if (isprop (Mdl, paramName))
              params.(paramName) = Mdl.(paramName);
            endif
          endfor
          this.ModelParameters = params;

        case 'GAM'
          ## Arguments to pass in fitcgam
          args = {};
          ## List of acceptable parameters for fitcdiscr
          GAMparams = {'PredictorNames', 'ResponseName', 'ClassNames', ...
                       'Cost', 'Formula', 'Knots', 'Order', ...
                       'LearningRate', 'NumIterations'};
          ## Set parameters
          for i = 1:numel (GAMparams)
            paramName = GAMparams{i};
            paramValue = Mdl.(paramName);
            if (! isempty (paramValue))
              args = [args, {paramName, paramValue}];
            endif
          endfor

          ## Interactions now holds the fitted pairs, which the constructor
          ## does not take as a specification.  The term matrix does, and it
          ## reproduces the parent's terms exactly rather than re-selecting
          ## them.  A formula names its own terms and is passed instead, so
          ## this must not be passed alongside one.
          if (isempty (Mdl.Formula) && ! isempty (Mdl.IntMatrix))
            args = [args, {'Interactions', Mdl.IntMatrix}];
          endif

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            tmp = fitcgam (this.X(idx, :), this.Y(idx,:), args{:});
            this.Trained{k} = compact (tmp);
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct ();
          paramList = {'Formula', 'Interactions', 'Knots', 'Order', 'DoF', ...
                       'LearningRate', 'NumIterations'};
          for i = 1:numel (paramList)
            paramName = paramList{i};
            if (isprop (Mdl, paramName))
              params.(paramName) = Mdl.(paramName);
            endif
          endfor
          this.ModelParameters = params;

        case 'KNN'
          ## Arguments to pass in fitcknn
          args = {};
          ## List of acceptable parameters for fitcknn.  ScoreTransform is
          ## deliberately absent: the parent applies it to the assembled
          ## scores, so a fold carrying it too would apply it twice.
          KNNparams = {'PredictorNames', 'ResponseName', 'ClassNames', ...
                       'Prior', 'Cost', 'BreakTies', ...
                       'NSMethod', 'BucketSize', 'NumNeighbors', 'Exponent', ...
                       'Scale', 'Cov', 'Distance', 'DistanceWeight', ...
                       'IncludeTies'};
          ## Set parameters
          for i = 1:numel (KNNparams)
            paramName = KNNparams{i};
            if (isprop (Mdl, paramName))
              paramValue = Mdl.(paramName);
              if (! isempty (paramValue))
                args = [args, {paramName, paramValue}];
              endif
            else
              switch (paramName)
                case 'Cov'
                  if (strcmpi (Mdl.Distance, 'mahalanobis') && ...
                      (! isempty (Mdl.DistParameter)))
                    args = [args, {'Cov', Mdl.DistParameter}];
                  endif
                case 'Exponent'
                  if (strcmpi (Mdl.Distance,'minkowski') && ...
                      (! isempty (Mdl.DistParameter)))
                    args = [args, {'Exponent', Mdl.DistParameter}];
                  endif
                case 'Scale'
                  if (strcmpi (Mdl.Distance,'seuclidean') && ...
                      (! isempty (Mdl.DistParameter)))
                    args = [args, {'Scale', Mdl.DistParameter}];
                  endif
              endswitch
            endif
          endfor

          ## Standardization is told by Mu, and is passed on: leaving it out
          ## refitted every fold on the raw scale while the parent reported
          ## itself standardized.  It is passed only when it is on, the
          ## constructor taking at most one of Standardize, Scale and Cov, so
          ## a model carrying either of those was never standardized.
          if (! isempty (Mdl.Mu))
            args = [args, {'Standardize', true}];
          endif

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            this.Trained{k} = fitcknn (this.X(idx, :), this.Y(idx,:), args{:});
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct ();
          paramList = {'NumNeighbors', 'Distance', 'DistParameter', ...
                       'NSMethod', 'DistanceWeight'};
          for i = 1:numel (paramList)
            paramName = paramList{i};
            if (isprop (Mdl, paramName))
              params.(paramName) = Mdl.(paramName);
            endif
          endfor
          this.ModelParameters = params;

        case 'NeuralNetwork'
          ## Arguments to pass in fitcnet
          args = {};
          ## List of acceptable parameters for fitcnet.  ScoreTransform is
          ## deliberately absent, as it is for the KNN above.
          NNparams = {'PredictorNames', 'ResponseName', 'ClassNames', ...
                      'LayerSizes', ...
                      'Activations', 'OutputLayerActivation', ...
                      'LearningRate', 'IterationLimit', 'DisplayInfo'};
          ## Set parameters
          for i = 1:numel (NNparams)
            paramName = NNparams{i};
            paramValue = Mdl.(paramName);
            if (! isempty (paramValue))
              args = [args, {paramName, paramValue}];
            endif
          endfor
          stdz = ! isempty (Mdl.Mu);
          args = [args, {'Standardize', stdz}];

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            tmp = fitcnet (this.X(idx, :), this.Y(idx,:), args{:});
            this.Trained{k} = compact (tmp);
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          params = struct ();
          paramList = {'LayerSizes', 'Activations', 'OutputLayerActivation', ...
                       'LearningRate', 'IterationLimit', 'Solver'};
          for i = 1:numel (paramList)
            paramName = paramList{i};
            if (isprop (Mdl, paramName))
              params.(paramName) = Mdl.(paramName);
            endif
          endfor
          this.ModelParameters = params;

        case 'SVM'
          ## Get ModelParameters structure from ClassificationSVM object
          params = Mdl.ModelParameters;

          ## Train model according to partition object
          for k = 1:this.KFold
            idx = training (this.Partition, k);
            ## Pass all arguments directly to fitcsvm
            tmp = fitcsvm (this.X(idx, :), this.Y(idx,:), ...
                           'Standardize', ! isempty (Mdl.Mu), ...
                           'PredictorNames', Mdl.PredictorNames, ...
                           'ResponseName', Mdl.ResponseName, ...
                           'ClassNames', Mdl.ClassNames, ...
                           'SVMtype', params.SVMtype, ...
                           'KernelFunction', params.KernelFunction, ...
                           'PolynomialOrder', params.PolynomialOrder, ...
                           'KernelScale', params.KernelScale, ...
                           'KernelOffset', params.KernelOffset, ...
                           'BoxConstraint', params.BoxConstraint, ...
                           'Nu', params.Nu, ...
                           'CacheSize', params.CacheSize, ...
                           'Tolerance', params.Tolerance, ...
                           'Shrinking', params.Shrinking);
            this.Trained{k} = compact (tmp);
          endfor

          ## Store ModelParameters to ClassificationPartitionedModel object
          this.ModelParameters = params;

      endswitch

      ## No fold carries the transform: the parent applies it once to the
      ## assembled score.  Cleared here rather than at each fit call so that
      ## a backing added later cannot reintroduce the double application by
      ## inheriting a class default of its own, which is how a cross
      ## validated GAM came to transform twice.  The KNN case stores a full
      ## ClassificationKNN where the other four store compacts, and both
      ## accept the assignment.
      ##
      ## Safe because the transform here is a caller's preference and not
      ## fitted content.  A fitted transform would not be: fitPosterior
      ## installs a sigmoid it has estimated, so if that is ever implemented
      ## through crossval it must not be cleared here.
      this = pushToFolds (this, 'ScoreTransform', 'none');

      this.Fitted = true;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedModel} {@var{label} =} kfoldPredict (@var{this})
    ## @deftypefnx {ClassificationPartitionedModel} {[@var{label}, @var{score}, @var{cost}] =} kfoldPredict (@var{this})
    ##
    ## Predict responses for observations not used for training in a
    ## cross-validated classification model.
    ##
    ## @code{@var{[label, Score, Cost]} = kfoldPredict (@var{this})}
    ## returns the predicted class labels, classification scores, and
    ## classification costs for the data used
    ## to train the cross-validated model @var{this}.
    ##
    ## @var{this} is a @code{ClassificationPartitionedModel} object.
    ## The function predicts the response for each observation that was
    ## held out during training in the cross-validation process.
    ##
    ## An observation that no fold held out is not predicted at all: its
    ## scores and costs are @qcode{NaN} and its label is missing, an empty
    ## character vector for a cell array of strings and @qcode{NaN} for a
    ## numeric response.  Under a @qcode{'Holdout'} partition that is every
    ## observation outside the holdout set.  @strong{This differs from
    ## MATLAB}, which reports @qcode{NaN} scores for those rows as we do but
    ## labels every one of them with the @emph{first} class, whatever their
    ## response: that label is the least-cost class of a row of @qcode{NaN}
    ## costs rather than a prediction any model made, and naming a class for
    ## an observation nothing scored would be wrong.  A logical response has
    ## no missing value to give, so those rows stay @qcode{false}.
    ##
    ## @multitable @columnfractions 0.28 0.7
    ## @headitem @var{Output} @tab @var{Description}
    ##
    ## @item @qcode{label} @tab Predicted class labels, returned as a
    ## vector or cell array. The type of @var{label} matches the type of
    ## @var{Y} in the original training data. Each element of @var{label}
    ## corresponds to the predicted class
    ## label for the corresponding row in @var{X}.
    ##
    ## @item @qcode{Score} @tab Classification scores, returned as a
    ## numeric matrix. Each row of @var{Score} corresponds to an observation,
    ## and each column corresponds to a class. The value in row @var{i} and
    ## column @var{j} is the
    ## classification score for class @var{j} for observation @var{i}.
    ##
    ## @item @qcode{Cost} @tab Classification costs, returned as a
    ## numeric matrix. Each row of @var{Cost} corresponds to an observation,
    ## and each column corresponds to a class. The value in row @var{i}
    ## and column @var{j} is the classification cost for class @var{j} for
    ## observation @var{i}. This output is optional and only returned if
    ## requested.
    ## @end multitable
    ##
    ## @seealso{ClassificationKNN, ClassificationSVM,
    ## ClassificationPartitionedModel}
    ## @end deftypefn
    function [label, Score, Cost] = kfoldPredict (this)

      ## Input validation
      ## Models whose compact form reports no cost.  CompactClassificationGAM
      ## returns a label and a score only, so asking it for a third output
      ## raised and kfoldPredict could not run on a cross-validated GAM.
      no_cost_models = {'GAM', 'NeuralNetwork', 'SVM'};
      no_cost = any (strcmp (this.CrossValidatedModel, no_cost_models));
      if (no_cost && nargout > 2)
        error (strcat ("ClassificationPartitionedModel.kfoldPredict:", ...
                       " 'Cost' output is not supported for %s cross", ...
                       " validated models."), this.CrossValidatedModel);
      endif

      ## The label vector starts missing rather than empty, so an observation
      ## no fold tests is reported as missing instead of carrying a class.
      ## Under a holdout partition that is most of them.  A logical response
      ## has no missing value to give, so those rows stay false.
      if (iscellstr (this.Y))
        label = repmat ({''}, this.NumObservations, 1);
      elseif (islogical (this.Y))
        label = false (this.NumObservations, 1);
      elseif (isnumeric (this.Y))
        label = nan (this.NumObservations, 1);
      elseif (ischar (this.Y))
        label = repmat (' ', this.NumObservations, size (this.Y, 2));
      endif

      ## Initialize the score and cost matrices
      Score = nan (this.NumObservations, classCount (this.ClassNames));
      Cost = nan (this.NumObservations, classCount (this.ClassNames));

      ## Predict label, score, and cost (if applicable) for each KFold partition
      for k = 1:this.KFold

        ## Get data and trained model for this fold
        testIdx = test (this.Partition, k);
        model = this.Trained{k};

        ## Train
        if (no_cost)
          [predictedLabel, score] = predict (model, this.X(testIdx, :));
        else
          [predictedLabel, score, cost] = predict (model, this.X(testIdx, :));
        endif

        ## Convert cell array of labels to appropriate type (if applicable)
        if (iscell (predictedLabel))
          if (isnumeric (this.Y))
            predictedLabel = cellfun (@str2num, predictedLabel);
          elseif (islogical (this.Y))
            predictedLabel = cellfun (@logical, predictedLabel);
          elseif (iscellstr (this.Y))
            predictedLabel = predictedLabel;
          endif
        endif

        ## Get labels, score, and cost (if applicable).  The rows are
        ## assigned and not the elements: a character matrix of labels was
        ## allocated one name per row above, and a linear index would write
        ## down its first column.
        label(testIdx, :) = predictedLabel;
        Score(testIdx, :) = score;
        if (nargout > 2)
          Cost(testIdx, :) = cost;
        endif

      endfor

      ## The folds never carry the transform: assigning ScoreTransform on a
      ## cross-validated model leaves every Trained{k} at 'none', in MATLAB as
      ## here, and it is applied once to the assembled scores instead.  It used
      ## to be accepted and never read, so the property was inert on every
      ## cross-validated classifier.  kfoldLoss reads these scores, so the
      ## margin-based losses follow the transform as they should.
      Score = this.STfun (Score);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedModel} {@var{L} =} kfoldLoss (@var{obj})
    ## @deftypefnx {ClassificationPartitionedModel} {@var{L} =} kfoldLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute the cross-validated classification loss.
    ##
    ## @code{@var{L} = kfoldLoss (@var{obj})} returns the fraction of
    ## observations the folds misclassify, each answered for by the fold's
    ## model that did not see it, which is what @code{kfoldPredict} returns.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationPartitionedModel} object.
    ## @end itemize
    ##
    ## @code{@var{L} = kfoldLoss (@dots{}, @var{name}, @var{value})} accepts
    ## the following @qcode{Name-Value} pairs.
    ##
    ## @multitable @columnfractions 0.24 0.76
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab @qcode{'classiferror'}, the default,
    ## @qcode{'classifcost'}, @qcode{'mincost'}, or a function handle called
    ## as
    ## @code{@var{lossfun} (@var{C}, @var{S}, @var{W}, @var{Cost})}, where
    ## @var{C} is a logical matrix with one true per row marking the true
    ## class, @var{S} the scores, @var{W} the weights and @var{Cost} the
    ## misclassification cost.
    ##
    ## @item @qcode{'Mode'} @tab @qcode{'average'}, the default, which returns
    ## one number over the observations of every fold asked for, or
    ## @qcode{'individual'}, which returns one number per fold.
    ##
    ## @item @qcode{'Folds'} @tab A vector of fold indices to restrict the
    ## loss to.  It defaults to every fold.
    ## @end multitable
    ##
    ## @seealso{ClassificationPartitionedModel, kfoldPredict}
    ## @end deftypefn
    function L = kfoldLoss (this, varargin)

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      ## Defaults, then the optional pairs
      LossFun = 'classiferror';
      Mode    = 'average';
      Folds   = 1:this.KFold;
      ## Only the losses whose value has been measured against R2024a are
      ## offered.  The margin-based ones are not: ours came out 2x MATLAB's
      ## for hinge and 4x for quadratic, the square of the same factor, which
      ## says its margin is not the score difference this would use.  Shipping
      ## a number that close to right and not right is worse than not
      ## shipping it.
      names   = {'classifcost', 'classiferror', 'mincost'};
      while (numel (varargin) > 0)
        if (! (ischar (varargin{1}) && isrow (varargin{1})))
          error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                         " parameter name must be a character vector."));
        endif
        switch (tolower (varargin{1}))

          case 'lossfun'
            LossFun = varargin{2};
            if (! (is_function_handle (LossFun) ||
                   (ischar (LossFun) && isrow (LossFun))))
              error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                             " 'LossFun' must be a character vector or a", ...
                             " function handle."));
            endif
            if (ischar (LossFun))
              LossFun = tolower (LossFun);
              if (! any (strcmp (LossFun, names)))
                error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                               " unsupported 'LossFun' value."));
              endif
            endif

          case 'mode'
            Mode = varargin{2};
            if (! (ischar (Mode) && isrow (Mode) &&
                   any (strcmpi (Mode, {'average', 'individual'}))))
              error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                             " 'Mode' must be either 'average' or", ...
                             " 'individual'."));
            endif

          case 'folds'
            Folds = varargin{2};
            if (! (isnumeric (Folds) && isvector (Folds)
                   && all (Folds == fix (Folds))
                   && all (Folds >= 1) && all (Folds <= this.KFold)))
              error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                             " 'Folds' must be a vector of fold indices", ...
                             " between 1 and KFold."));
            endif

          otherwise
            error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                           " invalid parameter name in optional paired", ...
                           " arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Every observation answered for by the fold that held it out
      [label, Score] = kfoldPredict (this);
      classes = this.ClassNames;
      K = classCount (classes);
      n = rows (this.X);

      ## Index of the true and the predicted class, so the cost matrix and
      ## the score matrix can be addressed by the same column.
      true_idx = zeros (n, 1);
      pred_idx = zeros (n, 1);
      for i = 1:n
        t = find (ismember (classes, this.Y(i)));
        p = find (ismember (classes, label(i)));
        if (! isempty (t))
          true_idx(i) = t(1);
        endif
        if (! isempty (p))
          pred_idx(i) = p(1);
        endif
      endfor

      if (strcmpi (Mode, 'individual'))
        L = nan (numel (Folds), 1);
        for i = 1:numel (Folds)
          idx = test (this.Partition, Folds(i));
          L(i) = foldLoss_ (this, idx, Score, true_idx, pred_idx, LossFun, K);
        endfor
      else
        idx = false (n, 1);
        for i = 1:numel (Folds)
          idx = idx | test (this.Partition, Folds(i));
        endfor
        L = foldLoss_ (this, idx, Score, true_idx, pred_idx, LossFun, K);
      endif

    endfunction


    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedModel} {@var{m} =} kfoldMargin (@var{obj})
    ##
    ## Classification margins of the cross-validated observations.
    ##
    ## @code{@var{m} = kfoldMargin (@var{obj})} returns an @math{Nx1} vector
    ## holding, for every observation, the score its own fold's model gave the
    ## true class less the largest score that model gave any other class.  A
    ## larger margin is a more confident correct answer and a negative one is
    ## a misclassification.  Every observation is scored by the fold that held
    ## it out, so no model answers for a row it was trained on.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationPartitionedModel} object.
    ## @end itemize
    ##
    ## Where the fold that held an observation out produced no score for it,
    ## the margin is @qcode{NaN}.  This method takes no optional arguments,
    ## as MATLAB's does not.
    ##
    ## @end deftypefn
    function m = kfoldMargin (this)

      m = foldMargin_ (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedModel} {@var{e} =} kfoldEdge (@var{obj})
    ## @deftypefnx {ClassificationPartitionedModel} {@var{e} =} kfoldEdge (@dots{}, @var{name}, @var{value})
    ##
    ## Classification edge of the cross-validated observations.
    ##
    ## @code{@var{e} = kfoldEdge (@var{obj})} returns the mean of the
    ## classification margins over every cross-validated observation, which is
    ## the mean of @code{kfoldMargin (@var{obj})}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationPartitionedModel} object.
    ## @end itemize
    ##
    ## @code{@var{e} = kfoldEdge (@dots{}, @var{name}, @var{value})} accepts
    ## the following @qcode{Name-Value} pairs.
    ##
    ## @multitable @columnfractions 0.24 0.76
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'Mode'} @tab @qcode{'average'}, the default, which returns
    ## one number over the observations of every fold asked for, or
    ## @qcode{'individual'}, which returns one number per fold.
    ##
    ## @item @qcode{'Folds'} @tab A vector of fold indices to restrict the
    ## edge to.  It defaults to every fold.
    ## @end multitable
    ##
    ## The observations of a selection are weighted uniformly and normalized
    ## over that selection, so a subset of folds is an average rather than a
    ## sum, exactly as @code{kfoldLoss} does.
    ##
    ## @end deftypefn
    function e = kfoldEdge (this, varargin)

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationPartitionedModel.kfoldEdge:", ...
                       " Name-Value arguments must be in pairs."));
      endif

      ## Defaults, then the optional pairs
      Mode  = 'average';
      Folds = 1:this.KFold;
      while (numel (varargin) > 0)
        if (! (ischar (varargin{1}) && isrow (varargin{1})))
          error (strcat ("ClassificationPartitionedModel.kfoldEdge:", ...
                         " parameter name must be a character vector."));
        endif
        switch (tolower (varargin{1}))

          case 'mode'
            Mode = varargin{2};
            if (! (ischar (Mode) && isrow (Mode) &&
                   any (strcmpi (Mode, {'average', 'individual'}))))
              error (strcat ("ClassificationPartitionedModel.kfoldEdge:", ...
                             " 'Mode' must be either 'average' or", ...
                             " 'individual'."));
            endif

          case 'folds'
            Folds = varargin{2};
            if (! (isnumeric (Folds) && isvector (Folds)
                   && all (Folds == fix (Folds))
                   && all (Folds >= 1) && all (Folds <= this.KFold)))
              error (strcat ("ClassificationPartitionedModel.kfoldEdge:", ...
                             " 'Folds' must be a vector of fold indices", ...
                             " between 1 and KFold."));
            endif

          otherwise
            error (strcat ("ClassificationPartitionedModel.kfoldEdge:", ...
                           " invalid parameter name in optional paired", ...
                           " arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      m = foldMargin_ (this);
      n = rows (this.X);

      if (strcmpi (Mode, 'individual'))
        e = nan (numel (Folds), 1);
        for i = 1:numel (Folds)
          idx = test (this.Partition, Folds(i));
          e(i) = foldEdge_ (this, m, idx);
        endfor
      else
        idx = false (n, 1);
        for i = 1:numel (Folds)
          idx = idx | test (this.Partition, Folds(i));
        endfor
        e = foldEdge_ (this, m, idx);
      endif

    endfunction

  endmethods

  methods (Access = private)

    ## Margin of every observation from the fold that held it out: the score
    ## of the true class less the largest score given any other class.  NaN
    ## where the response is not one of the trained classes, since a margin
    ## needs a true class to measure from.
    function m = foldMargin_ (this)
      [~, Score] = kfoldPredict (this);
      classes = this.ClassNames;
      K = classCount (classes);
      n = rows (this.X);
      m = nan (n, 1);
      for i = 1:n
        t = find (ismember (classes, this.Y(i)));
        if (isempty (t))
          continue;
        endif
        t = t(1);
        if (K < 2)
          ## Degenerate single-class model: nothing competes, so the score
          ## itself is the whole margin.
          m(i) = Score(i, t);
        else
          other = Score(i, [1:t-1, t+1:K]);
          m(i) = Score(i, t) - max (other);
        endif
      endfor
    endfunction

    ## Edge over the observations selected by IDX: the margins weighted
    ## uniformly and normalized over that selection, so a subset of folds is
    ## an average rather than a sum, as foldLoss_ does for the losses.
    function e = foldEdge_ (this, m, idx)
      if (! any (idx))
        e = NaN;
        return;
      endif
      mi = m(idx);
      e = sum (mi) / numel (mi);
    endfunction

    ## Loss over the observations selected by IDX, weighted uniformly and
    ## normalized over that selection, so a subset of folds is an average
    ## rather than a sum.
    function L = foldLoss_ (this, idx, Score, true_idx, pred_idx, LossFun, K)
      if (! any (idx) || ! any (true_idx(idx)))
        L = NaN;
        return;
      endif
      S = Score(idx, :);
      ti = true_idx(idx);
      pi_ = pred_idx(idx);
      m = sum (idx);
      W = ones (m, 1) / m;

      if (is_function_handle (LossFun))
        C = false (m, K);
        C(sub2ind ([m, K], (1:m)', ti)) = true;
        L = LossFun (C, S, W, this.Cost);
        if (! (isnumeric (L) && isscalar (L)))
          error (strcat ("ClassificationPartitionedModel.kfoldLoss:", ...
                         " 'LossFun' must return a numeric scalar."));
        endif
        return;
      endif

      switch (LossFun)
        case 'classiferror'
          L = sum (W .* (pi_ != ti));

        case 'classifcost'
          L = 0;
          for i = 1:m
            L = L + W(i) * this.Cost(ti(i), pi_(i));
          endfor

        case 'mincost'
          L = 0;
          for i = 1:m
            [~, k] = min (S(i,:) * this.Cost);
            L = L + W(i) * this.Cost(ti(i), k);
          endfor

      endswitch
    endfunction

    ## Carry an assigned property down into every fold.  MATLAB's Cost and
    ## Prior reach the fold models, so kfoldPredict and kfoldLoss answer under
    ## what was assigned rather than under what was fitted.  The folds are
    ## absent while the constructor is still running, and the two properties
    ## are only ever pushed for the learners that accept the assignment.
    function this = pushToFolds (this, name, val)
      for k = 1:numel (this.Trained)
        if (! isempty (this.Trained{k}))
          this.Trained{k}.(name) = val;
        endif
      endfor
    endfunction

  endmethods

endclassdef


%!demo
%!
%! load fisheriris
%! x = meas;
%! y = species;
%!
%! ## Create a KNN classifier model
%! obj = fitcknn (x, y, 'NumNeighbors', 5, 'Standardize', 1);
%!
%! ## Create a partition for 5-fold cross-validation
%! partition = cvpartition (y, 'KFold', 5);
%!
%! ## Create the ClassificationPartitionedModel object
%! cvModel = crossval (obj, 'cvPartition', partition)

%!demo
%!
%! load fisheriris
%! x = meas;
%! y = species;
%!
%! ## Create a KNN classifier model
%! obj = fitcknn (x, y, 'NumNeighbors', 5, 'Standardize', 1);
%!
%! ## Create the ClassificationPartitionedModel object
%! cvModel = crossval (obj);
%!
%! ## Predict the class labels for the observations not used for training
%! [label, score, cost] = kfoldPredict (cvModel);
%! fprintf ("Cross-validated accuracy = %1.2f%% (%d/%d)\n", ...
%!          sum (strcmp (label, y)) / numel (y) *100, ...
%!          sum (strcmp (label, y)), numel (y))

## Tests
%!test
%! load fisheriris
%! a = fitcdiscr (meas, species, 'gamma', 0.3);
%! cvModel = crossval (a, 'KFold', 5);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (cvModel.NumObservations, 150);
%! assert_equal (numel (cvModel.Trained), 5);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationDiscriminant");
%! assert_equal (cvModel.CrossValidatedModel, "Discriminant");
%! assert_equal (cvModel.KFold, 5);
%!test
%! load fisheriris
%! a = fitcdiscr (meas, species, 'gamma', 0.5, 'fillcoeffs', 'off');
%! cvModel = crossval (a, 'HoldOut', 0.3);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal ({cvModel.X, cvModel.Y}, {meas, species});
%! assert_equal (cvModel.NumObservations, 150);
%! assert_equal (numel (cvModel.Trained), 1);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationDiscriminant");
%! assert_equal (cvModel.CrossValidatedModel, "Discriminant");
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ['a'; 'a'; 'b'; 'b'];
%! a = fitcgam (x, y, 'Interactions', 'all');
%! cvModel = crossval (a, 'KFold', 2);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (numel (cvModel.Trained), 2);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationGAM");
%! assert_equal (cvModel.CrossValidatedModel, "GAM");
%! assert_equal (cvModel.KFold, 2);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ['a'; 'a'; 'b'; 'b'];
%! a = fitcgam (x, y);
%! cvModel = crossval (a, 'LeaveOut', 'on');
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal ({cvModel.X, cvModel.Y}, {x, y});
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (numel (cvModel.Trained), 4);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationGAM");
%! assert_equal (cvModel.CrossValidatedModel, "GAM");
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ['a'; 'a'; 'b'; 'b'];
%! a = fitcknn (x, y);
%! partition = cvpartition (y, 'KFold', 2);
%! cvModel = ClassificationPartitionedModel (a, partition);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (class (cvModel.Trained{1}), "ClassificationKNN");
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (cvModel.ModelParameters.NumNeighbors, 1);
%! assert_equal (cvModel.ModelParameters.NSMethod, "kdtree");
%! assert_equal (cvModel.ModelParameters.Distance, "euclidean");
%! assert_equal (isempty (cvModel.Trained{1}.Mu), true);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ['a'; 'a'; 'b'; 'b'];
%! a = fitcknn (x, y, 'NSMethod', 'exhaustive');
%! partition = cvpartition (y, 'HoldOut', 0.2);
%! cvModel = ClassificationPartitionedModel (a, partition);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (class (cvModel.Trained{1}), "ClassificationKNN");
%! assert_equal ({cvModel.X, cvModel.Y}, {x, y});
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (cvModel.ModelParameters.NumNeighbors, 1);
%! assert_equal (cvModel.ModelParameters.NSMethod, "exhaustive");
%! assert_equal (cvModel.ModelParameters.Distance, "euclidean");
%! assert_equal (isempty (cvModel.Trained{1}.Mu), true);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = ['a'; 'a'; 'b'; 'b'];
%! k = 2;
%! a = fitcknn (x, y, 'NumNeighbors' ,k);
%! partition = cvpartition (numel (y), 'LeaveOut');
%! cvModel = ClassificationPartitionedModel (a, partition);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (class (cvModel.Trained{1}), "ClassificationKNN");
%! assert_equal ({cvModel.X, cvModel.Y}, {x, y});
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (cvModel.ModelParameters.NumNeighbors, k);
%! assert_equal (cvModel.ModelParameters.NSMethod, "kdtree");
%! assert_equal (cvModel.ModelParameters.Distance, "euclidean");
%! assert_equal (isempty (cvModel.Trained{1}.Mu), true);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = {'a'; 'a'; 'b'; 'b'};
%! a = fitcnet (x, y, 'IterationLimit', 50);
%! cvModel = crossval (a, 'KFold', 2);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (numel (cvModel.Trained), 2);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationNeuralNetwork");
%! assert_equal (cvModel.CrossValidatedModel, "NeuralNetwork");
%! assert_equal (cvModel.KFold, 2);
%!test
%! x = [1, 2, 3; 4, 5, 6; 7, 8, 9; 3, 2, 1];
%! y = {'a'; 'a'; 'b'; 'b'};
%! a = fitcnet (x, y, 'LayerSizes', [5, 3]);
%! cvModel = crossval (a, 'LeaveOut', 'on');
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal ({cvModel.X, cvModel.Y}, {x, y});
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (numel (cvModel.Trained), 4);
%! assert_equal (class (cvModel.Trained{1}), "CompactClassificationNeuralNetwork");
%! assert_equal (cvModel.CrossValidatedModel, "NeuralNetwork");
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%! SVMModel = fitcsvm (x,y);
%! CVMdl = crossval (SVMModel, 'KFold', 5);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (CVMdl.KFold == 5, true)
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationSVM")
%! assert_equal (CVMdl.CrossValidatedModel, "SVM");
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%! obj = fitcsvm (x, y);
%! CVMdl = crossval (obj, 'HoldOut', 0.2);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationSVM")
%! assert_equal (CVMdl.CrossValidatedModel, "SVM");
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'setosa');
%! x = meas(inds, 3:4);
%! y = grp2idx (species(inds));
%! obj = fitcsvm (x, y);
%! CVMdl = crossval (obj, 'LeaveOut', 'on');
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationSVM")
%! assert_equal (CVMdl.CrossValidatedModel, "SVM");

## The KNN short name, the one of the five that no test pinned.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcknn (meas, species), 'KFold', 5);
%! assert_equal (CVMdl.CrossValidatedModel, "KNN");

## Test input validation for ClassificationPartitionedModel
## Cross-validating a GAM rebuilds each fold from the term matrix, the
## interaction pairs no longer being a form the constructor accepts.  Every
## fold carries the parent's terms.
%!test
%! k = (1:60)';
%! X = [mod(k*7,11)-5, mod(k*3,11)-5, mod(k*5,11)-5];
%! y = double (X(:,1).*X(:,2) > 0) + 1;
%! Mdl = fitcgam (X, y, "Interactions", "all");
%! cv = crossval (Mdl, "KFold", 3);
%! assert_equal (numel (cv.Trained), 3);
%! assert_equal (cv.Trained{1}.Interactions, [1, 2; 1, 3; 2, 3]);
%! assert_equal (cv.Trained{3}.Interactions, Mdl.Interactions);

## kfoldPredict writes whole names into its result: the labels are laid out
## one name per row, and assigning by element would fill the first column.
%!test
%! load fisheriris
%! b = ! strcmp (species, "setosa");
%! X = meas(b,:); Ys = species(b); Yc = char (Ys);
%! rand ("state", 1); Mc = fitcdiscr (X, Yc);
%! rand ("state", 1); Ms = fitcdiscr (X, Ys);
%! rand ("state", 2); cvc = crossval (Mc, "KFold", 3);
%! rand ("state", 2); cvs = crossval (Ms, "KFold", 3);
%! p = kfoldPredict (cvc);
%! assert_equal (columns (p), 10);
%! assert_equal (cellstr (p), kfoldPredict (cvs));

%!error<ClassificationPartitionedModel: too few input arguments.> ...
%! ClassificationPartitionedModel ()
%!error<ClassificationPartitionedModel: too few input arguments.> ...
%! ClassificationPartitionedModel (ClassificationKNN (ones (4,2), ones (4,1)))
%!error<ClassificationPartitionedModel: unsupported model type.> ...
%! ClassificationPartitionedModel (RegressionGAM (ones (40,2), ...
%! randi ([1, 2], 40, 1)), cvpartition (randi ([1, 2], 40, 1), 'Holdout', 0.3))
%!error<ClassificationPartitionedModel: invalid 'cvpartition' object.> ...
%! ClassificationPartitionedModel (ClassificationKNN (ones (4,2), ...
%! ones (4,1)), 'Holdout')

## Test for kfoldPredict
%!test
%! load fisheriris
%! a = fitcdiscr (meas, species, 'gamma', 0.5, 'fillcoeffs', 'off');
%! cvModel = crossval (a, 'Kfold', 4);
%! [label, score, cost] = kfoldPredict (cvModel);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal ({cvModel.X, cvModel.Y}, {meas, species});
%! assert_equal (cvModel.NumObservations, 150);
%!# assert_equal (label, {"b"; "b"; "a"; "a"});
%!# assert_equal (score, [4.5380e-01, 5.4620e-01; 2.4404e-01, 7.5596e-01; ...
%!#         9.9392e-01, 6.0844e-03; 9.9820e-01, 1.8000e-03], 1e-4);
%!# assert_equal (cost, [5.4620e-01, 4.5380e-01; 7.5596e-01, 2.4404e-01; ...
%!#         6.0844e-03, 9.9392e-01; 1.8000e-03, 9.9820e-01], 1e-4);
%!test
%! x = ones (4, 11);
%! y = {'a'; 'a'; 'b'; 'b'};
%! k = 3;
%! a = fitcknn (x, y, 'NumNeighbors', k);
%! partition = cvpartition (numel (y), 'LeaveOut');
%! cvModel = ClassificationPartitionedModel (a, partition);
%! [label, score, cost] = kfoldPredict (cvModel);
%! assert_equal (class (cvModel), "ClassificationPartitionedModel");
%! assert_equal ({cvModel.X, cvModel.Y}, {x, y});
%! assert_equal (cvModel.NumObservations, 4);
%! assert_equal (cvModel.ModelParameters.NumNeighbors, k);
%! assert_equal (cvModel.ModelParameters.NSMethod, "exhaustive");
%! assert_equal (cvModel.ModelParameters.Distance, "euclidean");
%! assert_equal (isempty (cvModel.Trained{1}.Mu), true);
%! assert_equal (label, {'b'; 'b'; 'a'; 'a'});
%! assert_equal (score, [0.3333, 0.6667; 0.3333, 0.6667; 0.6667, 0.3333; ...
%!          0.6667, 0.3333], 1e-4);
%! assert_equal (cost, [0.6667, 0.3333; 0.6667, 0.3333; 0.3333, 0.6667; ...
%!          0.3333, 0.6667], 1e-4);

## Test input validation for kfoldPredict
%!error<ClassificationPartitionedModel.kfoldPredict: 'Cost' output is not supported for SVM cross validated models.> ...
%! [label, score, cost] = kfoldPredict (crossval (ClassificationSVM (ones (40,2), randi ([1, 2], 40, 1))))
%!error<ClassificationPartitionedModel.kfoldPredict: 'Cost' output is not supported for NeuralNetwork cross validated models.> ...
%! [label, score, cost] = kfoldPredict (crossval (ClassificationNeuralNetwork (ones (40,2), randi ([1, 2], 40, 1))))

## The partition, the stored data and NumObservations all describe the same
## set of observations.  A row dropped for a missing value is outside all
## three: before this the partition covered it, NumObservations did not, and
## the folds trained on rows their own fit then discarded.
%!test
%! randn ('seed', 42);
%! lab = double (randn (40, 1) > 0) + 1;
%! X = [randn(40, 2); NaN, 1; 2, NaN];
%! Y = [lab; 1; 2];
%! Mdl = fitcsvm (X, Y);
%! assert_equal (Mdl.NumObservations, 42);
%! assert_equal (Mdl.RowsUsed, []);
%! CVMdl = crossval (Mdl, 'KFold', 4);
%! assert_equal (CVMdl.NumObservations, 42);
%! assert_equal (rows (CVMdl.X), 42);
%! assert_equal (rows (CVMdl.Y), 42);
%! assert_equal (CVMdl.Partition.NumObservations, 42);
%! assert_equal (numel (kfoldPredict (CVMdl)), 42);

## The folds stay stratified, the response being passed to cvpartition
## rather than a bare count.
%!test
%! randn ('seed', 7);
%! X = randn (60, 2);
%! Y = [ones(20, 1); 2 * ones(40, 1)];
%! CVMdl = crossval (fitcsvm (X, Y), 'KFold', 4);
%! for k = 1:4
%!   idx = test (CVMdl.Partition, k);
%!   assert_equal (sum (CVMdl.Y(idx) == 1) >= 4, true);
%!   assert_equal (sum (CVMdl.Y(idx) == 2) >= 8, true);
%! endfor

## kfoldLoss defaults to the misclassification rate of the out-of-fold
## answers, which is the identity R2024a shows, and holds for every model
## type the class accepts.
%!test
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! rand ('seed', 42);
%! cvp = cvpartition (Y, 'KFold', 5);
%! for f = {'fitcsvm', 'fitcnet', 'fitcknn', 'fitcdiscr', 'fitcgam'}
%!   CVMdl = crossval (feval (f{1}, X, Y), 'CVPartition', cvp);
%!   assert_equal (kfoldLoss (CVMdl), ...
%!                 mean (! strcmp (kfoldPredict (CVMdl), CVMdl.Y)), 1e-12);
%! endfor

## Leave-one-out on a two-class fixture agrees with R2024a to the digit for
## every loss offered: it reports 0.02 for all three.
%!test
%! load fisheriris
%! idx = [51:75, 101:125];
%! CVMdl = crossval (fitcdiscr (meas(idx,:), species(idx)), 'Leaveout', 'on');
%! assert_equal (CVMdl.KFold, 50);
%! assert_equal (kfoldLoss (CVMdl), 0.02, 1e-12);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', 'classifcost'), 0.02, 1e-12);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', 'mincost'), 0.02, 1e-12);

## 'individual' is one number per fold, 'Folds' restricts to those named, and
## a handle is called with the class indicator, the scores, the weights and
## the cost, its weights summing to one as MATLAB's do.
%!test
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! rand ('seed', 42);
%! cvp = cvpartition (Y, 'KFold', 5);
%! CVMdl = crossval (fitcsvm (X, Y), 'CVPartition', cvp);
%! L = kfoldLoss (CVMdl, 'Mode', 'individual');
%! assert_equal (size (L), [5, 1]);
%! lab = kfoldPredict (CVMdl);
%! k2 = test (cvp, 2);
%! assert_equal (L(2), mean (! strcmp (lab(k2), CVMdl.Y(k2))), 1e-12);
%! sel = test (cvp, 1) | test (cvp, 3);
%! assert_equal (kfoldLoss (CVMdl, 'Folds', [1, 3]), ...
%!               mean (! strcmp (lab(sel), CVMdl.Y(sel))), 1e-12);
%! assert_equal (numel (kfoldLoss (CVMdl, 'Folds', [2, 4], ...
%!                                 'Mode', 'individual')), 2);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', @(C, S, W, Cost) sum (W)), ...
%!               1, 1e-12);

## Cost reaches the cost-aware losses.
%!test
%! load fisheriris
%! X = meas(51:150,:);
%! Y = species(51:150);
%! rand ('seed', 42);
%! CVMdl = crossval (fitcsvm (X, Y, 'Cost', [0, 4; 1, 0]), 'KFold', 5);
%! assert_equal (CVMdl.Cost, [0, 4; 1, 0]);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', 'classifcost') >= ...
%!               kfoldLoss (CVMdl, 'LossFun', 'classiferror'), true);

## Test input validation for kfoldLoss
%!shared CVK
%! load fisheriris
%! rand ('seed', 42);
%! CVK = crossval (fitcsvm (meas(51:150,:), species(51:150)), 'KFold', 4);
%!error<ClassificationPartitionedModel.kfoldLoss: Name-Value arguments must be in pairs.> ...
%! kfoldLoss (CVK, 'Mode')
%!error<ClassificationPartitionedModel.kfoldLoss: parameter name must be a character vector.> ...
%! kfoldLoss (CVK, 5, 1)
%!error<ClassificationPartitionedModel.kfoldLoss: 'LossFun' must be a character vector or a function handle.> ...
%! kfoldLoss (CVK, 'LossFun', 5)
%!error<ClassificationPartitionedModel.kfoldLoss: unsupported 'LossFun' value.> ...
%! kfoldLoss (CVK, 'LossFun', 'hinge')
%!error<ClassificationPartitionedModel.kfoldLoss: 'LossFun' must return a numeric scalar.> ...
%! kfoldLoss (CVK, 'LossFun', @(C, S, W, Cost) [1, 2])
%!error<ClassificationPartitionedModel.kfoldLoss: 'Mode' must be either 'average' or 'individual'.> ...
%! kfoldLoss (CVK, 'Mode', 'nope')
%!error<ClassificationPartitionedModel.kfoldLoss: 'Folds' must be a vector of fold indices between 1 and KFold.> ...
%! kfoldLoss (CVK, 'Folds', 0)
%!error<ClassificationPartitionedModel.kfoldLoss: invalid parameter name in optional paired arguments.> ...
%! kfoldLoss (CVK, 'Nope', 1)

## Standardization reaches the folds.  It did not for the KNN, whose refit
## arguments omitted it, so every fold trained on the raw scale while the
## parent reported itself standardized and nothing said so.
%!test
%! load fisheriris
%! X = meas(:,1:3);
%! X(:,1) = X(:,1) * 1000;
%! CVMdl = crossval (fitcknn (X, species, 'Standardize', true), 'KFold', 5);
%! assert_equal (isempty (CVMdl.Trained{1}.Mu), false);

%!test
%! load fisheriris
%! X = meas(:,1:3);
%! X(:,1) = X(:,1) * 1000;
%! cvp = cvpartition (species, 'KFold', 5);
%! CVMdl = crossval (fitcknn (X, species, 'Standardize', true), ...
%!                   'CVPartition', cvp);
%! hit = 0;
%! for k = 1:5
%!   tr = training (cvp, k);
%!   te = test (cvp, k);
%!   Mdl = fitcknn (X(tr,:), species(tr), 'Standardize', true);
%!   hit += sum (strcmp (predict (Mdl, X(te,:)), species(te)));
%! endfor
%! assert_equal (sum (strcmp (kfoldPredict (CVMdl), species)), hit);

## The weights of the model that was cross validated are carried over, as
## the regression counterpart already did.
%!test
%! load fisheriris
%! Mdl = fitcknn (meas, species);
%! CVMdl = crossval (Mdl, 'KFold', 3);
%! assert_equal (CVMdl.W, Mdl.W);
%! assert_equal (size (CVMdl.W), [150, 1]);

## BinEdges is an empty cell, and a cell rather than an empty matrix: code
## that reaches into it with cellfun works against MATLAB and used to fail
## here on the type alone.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcknn (meas, species), 'KFold', 3);
%! assert_equal (class (CVMdl.BinEdges), 'cell');
%! assert_equal (CVMdl.BinEdges, {});

## An assigned Cost reaches every fold, so the folds predict under what was
## assigned rather than under what was fitted.  Values measured on R2024a.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 3);
%! C = double (! eye (3)); C(1,2) = 4;
%! CVMdl.Cost = C;
%! assert_equal (CVMdl.Trained{1}.Cost, C);
%! assert_equal (CVMdl.Trained{3}.Cost, C);

%!test
%! load fisheriris
%! CVMdl = crossval (fitcknn (meas, species), 'KFold', 3);
%! C = double (! eye (3)); C(1,2) = 4;
%! CVMdl.Cost = C;
%! assert_equal (CVMdl.Trained{2}.Cost, C);

## The discriminant is the only learner whose priors may be assigned after
## the fit, and they reach the folds too.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 3);
%! CVMdl.Prior = [0.6, 0.2, 0.2];
%! assert_equal (CVMdl.Prior, [0.6, 0.2, 0.2], 1e-15);
%! assert_equal (CVMdl.Trained{1}.Prior, [0.6, 0.2, 0.2], 1e-15);

## Priors are normalized, as they are on the learner itself.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 3);
%! CVMdl.Prior = [3, 1, 1];
%! assert_equal (CVMdl.Prior, [0.6, 0.2, 0.2], 1e-15);

%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 3);
%! CVMdl.Prior = 'uniform';
%! assert_equal (CVMdl.Prior, [1, 1, 1] / 3, 1e-15);

## An SVM's costs enter its box constraint while it is fitted, so MATLAB
## refuses them afterwards on the cross-validated model and so do we.
%!error<ClassificationPartitionedModel: cannot assign 'Cost' on a cross-validated ClassificationSVM, whose costs are consumed while it is fitted.> ...
%! load fisheriris; ...
%! b = ismember (species, {'setosa', 'versicolor'}); ...
%! CVMdl = crossval (fitcsvm (meas(b,:), species(b)), 'KFold', 3); ...
%! CVMdl.Cost = [0, 4; 1, 0];
%!error<ClassificationPartitionedModel: 'Prior' can only be assigned on a cross-validated ClassificationDiscriminant.> ...
%! load fisheriris; ...
%! CVMdl = crossval (fitcknn (meas, species), 'KFold', 3); ...
%! CVMdl.Prior = [0.6, 0.2, 0.2];
%!error<ClassificationPartitionedModel: 'Prior' can only be assigned on a cross-validated ClassificationDiscriminant.> ...
%! load fisheriris; ...
%! CVMdl = crossval (fitcnet (meas, species), 'KFold', 3); ...
%! CVMdl.Prior = [0.6, 0.2, 0.2];
%!error<ClassificationPartitionedModel: the elements in 'Prior' do not correspond to the selected classes in Y.> ...
%! load fisheriris; ...
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 3); ...
%! CVMdl.Prior = [0.5, 0.5];

## A GAM backing transforms once, like every other.  Its class default is
## 'logit', so the fold has to be pinned at 'none' or the parent's transform
## lands on scores the fold has already transformed and the rows sum to about
## 1.23 instead of one.
%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! CVMdl = crossval (fitcgam (meas(inds,:), species(inds)), 'KFold', 3);
%! assert_equal (CVMdl.ScoreTransform, 'logit');
%! assert_equal (CVMdl.Trained{1}.ScoreTransform, 'none');
%! [~, s] = kfoldPredict (CVMdl);
%! assert_equal (sum (s, 2), ones (rows (s), 1), 1e-12);

## ScoreTransform is applied to the assembled scores, not carried into the
## folds.  MATLAB leaves every Trained{k} at 'none' and transforms at the
## wrapper; the property used to be stored here and never read.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 3);
%! [~, s0] = kfoldPredict (CVMdl);
%! CVMdl.ScoreTransform = 'doublelogit';
%! [~, s1] = kfoldPredict (CVMdl);
%! assert_equal (CVMdl.Trained{1}.ScoreTransform, 'none');
%! assert_equal (s1, 1 ./ (1 + exp (-2 * s0)), 1e-12);

## The transform reaches kfoldLoss, which reads those scores.  A loss that
## sums them shows it; mincost and classiferror do not move on this fixture,
## since a monotone transform leaves the chosen labels where they were.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 3);
%! f = @(C, S, W, Cost) sum (S(:));
%! [~, s0] = kfoldPredict (CVMdl);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', f), 150, 1e-12);
%! CVMdl.ScoreTransform = 'doublelogit';
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', f), ...
%!               sum (sum (1 ./ (1 + exp (-2 * s0)))), 1e-12);

## 'none' is the identity, so the default transforms nothing.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcknn (meas, species), 'KFold', 3);
%! [~, s0] = kfoldPredict (CVMdl);
%! CVMdl.ScoreTransform = 'none';
%! [~, s1] = kfoldPredict (CVMdl);
%! assert_equal (s1, s0);

## A transform carried by the model being cross-validated moves to the parent
## and is not handed to the folds.  R2024a leaves every Trained@{k@} at 'none'
## for all five backings; the KNN and the network used to be given it as well
## and so applied it a second time.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcknn (meas, species, ...
%!                            'ScoreTransform', 'doublelogit'), 'KFold', 3);
%! assert_equal (CVMdl.ScoreTransform, 'doublelogit');
%! assert_equal (CVMdl.Trained{1}.ScoreTransform, 'none');

%!test
%! load fisheriris
%! CVMdl = crossval (fitcnet (meas, species, ...
%!                            'ScoreTransform', 'doublelogit'), 'KFold', 3);
%! assert_equal (CVMdl.ScoreTransform, 'doublelogit');
%! assert_equal (CVMdl.Trained{1}.ScoreTransform, 'none');

%!test
%! load fisheriris
%! inds = ! strcmp (species, 'virginica');
%! CVMdl = crossval (fitcsvm (meas(inds,:), species(inds), ...
%!                            'ScoreTransform', 'doublelogit'), 'KFold', 3);
%! assert_equal (CVMdl.ScoreTransform, 'doublelogit');
%! assert_equal (CVMdl.Trained{1}.ScoreTransform, 'none');

## And it is applied once, not twice.  A KNN is deterministic given the
## partition, so the two fits share their folds and the untransformed one is
## an exact baseline: doublelogit applied twice would not match it.
%!test
%! load fisheriris
%! c = cvpartition (species, 'KFold', 3);
%! CV0 = crossval (fitcknn (meas, species), 'CVPartition', c);
%! CV1 = crossval (fitcknn (meas, species, ...
%!                          'ScoreTransform', 'doublelogit'), 'CVPartition', c);
%! [~, s0] = kfoldPredict (CV0);
%! [~, s1] = kfoldPredict (CV1);
%! assert_equal (s1, 1 ./ (1 + exp (-2 * s0)), 1e-12);

%!test
%! load fisheriris
%! CVMdl = crossval (fitcnet (meas, species, ...
%!                            'ScoreTransform', 'doublelogit'), 'KFold', 3);
%! [~, s1] = kfoldPredict (CVMdl);
%! CVMdl.ScoreTransform = 'none';
%! [~, s0] = kfoldPredict (CVMdl);
%! assert_equal (s1, 1 ./ (1 + exp (-2 * s0)), 1e-12);

%!error<ClassificationPartitionedModel: the number of rows and columns in 'Cost' must correspond to selected classes in Y.> ...
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 3);
%! CVMdl.Cost = 1:9;

## foldLoss_ is a private helper and stays out of the method list.
%!test
%! assert_equal (any (strcmp (methods ("ClassificationPartitionedModel"), ...
%!                           "foldLoss_")), false);

## kfoldMargin answers for every observation, from the fold that held it out.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcknn (meas, species), 'Leaveout', 'on');
%! m = kfoldMargin (CVMdl);
%! assert_equal (size (m), [150, 1]);
%! assert_equal (all (m >= -1 & m <= 1), true);

## kfoldEdge is the mean of the margins.  R2024a returns 0.92 on this
## leave-one-out KNN and ours reproduces it exactly.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcknn (meas, species), 'Leaveout', 'on');
%! assert_equal (kfoldEdge (CVMdl), 0.92, 1e-14);
%! assert_equal (kfoldEdge (CVMdl), mean (kfoldMargin (CVMdl)), 1e-14);

## 'Mode' selects one number per fold or one over them all.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 5);
%! assert_equal (size (kfoldEdge (CVMdl, 'Mode', 'individual')), [5, 1]);
%! assert_equal (isscalar (kfoldEdge (CVMdl, 'Mode', 'average')), true);

## 'Folds' restricts the selection, and one fold's edge is that fold's mean.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 5);
%! m = kfoldMargin (CVMdl);
%! idx = test (CVMdl.Partition, 2);
%! assert_equal (kfoldEdge (CVMdl, 'Folds', 2), mean (m(idx)), 1e-12);
%! assert_equal (kfoldEdge (CVMdl, 'Folds', 1:5), kfoldEdge (CVMdl), 1e-12);

%!error<ClassificationPartitionedModel.kfoldEdge: Name-Value arguments must be in pairs.> ...
%! load fisheriris
%! kfoldEdge (crossval (fitcdiscr (meas, species), 'KFold', 3), 'Mode')
%!error<ClassificationPartitionedModel.kfoldEdge: parameter name must be a character vector.> ...
%! load fisheriris
%! kfoldEdge (crossval (fitcdiscr (meas, species), 'KFold', 3), 5, 'average')
%!error<ClassificationPartitionedModel.kfoldEdge: 'Mode' must be either 'average' or 'individual'.> ...
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'KFold', 3);
%! kfoldEdge (CVMdl, 'Mode', 'cumulative')
%!error<ClassificationPartitionedModel.kfoldEdge: 'Folds' must be a vector of fold indices between 1 and KFold.> ...
%! load fisheriris
%! kfoldEdge (crossval (fitcdiscr (meas, species), 'KFold', 3), 'Folds', 7)
%!error<ClassificationPartitionedModel.kfoldEdge: invalid parameter name in optional paired arguments.> ...
%! load fisheriris
%! kfoldEdge (crossval (fitcdiscr (meas, species), 'KFold', 3), 'bogus', 1)

## A holdout partition predicts the holdout set and leaves the rest missing.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'Holdout', 0.2);
%! idx = test (CVMdl.Partition, 1);
%! [label, Score] = kfoldPredict (CVMdl);
%! assert_equal (iscellstr (label), true);
%! assert_equal (all (strcmp (label(! idx), '')), true);
%! assert_equal (any (strcmp (label(idx), '')), false);
%! assert_equal (all (all (isnan (Score(! idx, :)))), true);
%! assert_equal (any (any (isnan (Score(idx, :)))), false);

## The holdout predictions are the fold's own, not a stand-in class.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'Holdout', 0.2);
%! idx = test (CVMdl.Partition, 1);
%! label = kfoldPredict (CVMdl);
%! assert_equal (label(idx), predict (CVMdl.Trained{1}, meas(idx,:)));
%! assert_equal (mean (strcmp (label(idx), species(idx))) > 0.9, true);

## A numeric response reports a missing label as NaN.
%!test
%! load fisheriris
%! y = grp2idx (species);
%! CVMdl = crossval (fitcdiscr (meas, y), 'Holdout', 0.2);
%! idx = test (CVMdl.Partition, 1);
%! label = kfoldPredict (CVMdl);
%! assert_equal (all (isnan (label(! idx))), true);
%! assert_equal (any (isnan (label(idx))), false);

## kfoldLoss and kfoldEdge answer over the holdout set alone.
%!test
%! load fisheriris
%! CVMdl = crossval (fitcdiscr (meas, species), 'Holdout', 0.2);
%! idx = test (CVMdl.Partition, 1);
%! label = kfoldPredict (CVMdl);
%! assert_equal (kfoldLoss (CVMdl), ...
%!               mean (! strcmp (label(idx), species(idx))), 1e-12);
%! m = kfoldMargin (CVMdl);
%! assert_equal (all (isnan (m(! idx))), true);
%! assert_equal (kfoldEdge (CVMdl), mean (m(idx)), 1e-12);

## The shared cost guard is in force here too, and the struct form is
## permuted into this model's class order.  The battery is on
## ClassificationDiscriminant.
%!test
%! load fisheriris
%! Mdl = crossval (fitcdiscr (meas, species), 'KFold', 3);
%! S = struct ('ClassNames', {{'virginica'; 'setosa'; 'versicolor'}}, ...
%!             'ClassificationCosts', [0, 1, 2; 3, 0, 4; 5, 6, 0]);
%! Mdl.Cost = S;
%! assert_equal (Mdl.Cost, [0, 4, 3; 6, 0, 5; 1, 2, 0]);

%!error<ClassificationPartitionedModel: 'Cost' must have zeros on its diagonal.> ...
%! load fisheriris
%! Mdl = crossval (fitcdiscr (meas, species), 'KFold', 3);
%! Mdl.Cost = ones (3);
