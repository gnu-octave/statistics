## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftp {statistics} ClassificationPartitionedLinear
##
## Cross-validated linear binary classifier.
##
## A @qcode{ClassificationPartitionedLinear} object holds one
## @qcode{ClassificationLinear} per fold of a partition, each fitted to the
## observations the fold trains on.  Every @code{kfold} method predicts each
## observation with the fold that held it @emph{out}, so the estimate it
## returns is an out-of-sample one.
##
## A @qcode{ClassificationLinear} stores no copy of its training data and so
## has no resubstitution methods and no @code{compact} form.  This class is
## what takes their place: cross-validation is the way a linear model is
## asked how it would do on data it has not seen.
##
## When the fold models carry a whole regularization path, every method
## returns one column per strength, in the order of the @qcode{'Lambda'}
## that was asked for.
##
## Create one with @code{fitclinear} and a cross-validation option, or
## directly.
##
## @seealso{fitclinear, ClassificationLinear, ClassificationPartitionedKernel}
## @end deftp

classdef ClassificationPartitionedLinear

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} ClassNames
    ##
    ## Names of the two classes
    ##
    ## A column of the same type as the response, shared by every fold.
    ## This property is read-only.
    ##
    ## @end deftp
    ClassNames             = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} Cost
    ##
    ## Cost of misclassifying an observation
    ##
    ## A square numeric matrix with one row and one column per class.  It is
    ## handed to every fold rather than re-derived by each.  This property is
    ## read-only.
    ##
    ## @end deftp
    Cost                   = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} Prior
    ##
    ## Prior probability of each class
    ##
    ## A numeric row vector summing to one, in the order of
    ## @qcode{ClassNames}.  Like the cost it is the parent's and is handed to
    ## every fold, so a fold of an unbalanced problem does not quietly adopt
    ## a prior of its own.  This property is read-only.
    ##
    ## @end deftp
    Prior                  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} ScoreTransform
    ##
    ## Transformation applied to the predicted scores
    ##
    ## A character vector naming a transformation, or the text of the
    ## function handle that was supplied.  The fold models carry no
    ## transform of their own; this one is applied once to the assembled
    ## scores.  This property is read-only.
    ##
    ## @end deftp
    ScoreTransform         = 'none';

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} CrossValidatedModel
    ##
    ## Name of the model that was cross-validated
    ##
    ## Always @qcode{'Linear'}, the short name MATLAB uses.  This property is
    ## read-only.
    ##
    ## @end deftp
    CrossValidatedModel    = 'Linear';

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} NumObservations
    ##
    ## Number of observations the partition covers
    ##
    ## A positive integer scalar, counting the rows that survived the removal
    ## of missing values.  This property is read-only.
    ##
    ## @end deftp
    NumObservations        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} Y
    ##
    ## Response of the retained observations
    ##
    ## In the type it was supplied in.  This property is read-only.
    ##
    ## @end deftp
    Y                      = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} W
    ##
    ## Observation weights
    ##
    ## An @math{Nx1} numeric vector summing to one, normalized within each
    ## class to that class's cost-adjusted prior.  This property is
    ## read-only.
    ##
    ## @end deftp
    W                      = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames         = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A row vector of column indices, empty when every predictor is
    ## numeric.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} ResponseName
    ##
    ## Name of the response
    ##
    ## A character vector, defaulting to @qcode{'Y'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    ResponseName           = 'Y';

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} Trained
    ##
    ## The models fitted to the folds
    ##
    ## A cell column with one @qcode{ClassificationLinear} per fold, each
    ## fitted to the observations its fold trains on.  This property is
    ## read-only.
    ##
    ## @end deftp
    Trained                = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} KFold
    ##
    ## Number of folds
    ##
    ## A positive integer scalar.  A holdout partition has one fold and a
    ## leave-one-out partition has as many as there are observations.  This
    ## property is read-only.
    ##
    ## @end deftp
    KFold                  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} Partition
    ##
    ## The partition itself
    ##
    ## A @code{cvpartition} object over the retained observations.  This
    ## property is read-only.
    ##
    ## @end deftp
    Partition              = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedLinear} {property} ModelParameters
    ##
    ## What was cross-validated, and how
    ##
    ## A structure carrying @qcode{Type}, @qcode{Method},
    ## @qcode{LearnerTemplates} and @qcode{NLearn}.  This property is
    ## read-only.
    ##
    ## @end deftp
    ModelParameters        = [];

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## The callable behind ScoreTransform.
    STfun                  = @(s) s;

    ## The predictors of the retained observations.  MATLAB does not report
    ## them and neither do we, but the kfold methods have to predict from
    ## something and the fold models hold no data of their own.
    X_                     = [];

    ## Number of regularization strengths the fold models carry, so that
    ## every method knows how many columns to return without asking one.
    NumLambda_             = 1;

  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedLinear} {@var{obj} =} ClassificationPartitionedLinear (@var{X}, @var{Y})
    ## @deftypefnx {ClassificationPartitionedLinear} {@var{obj} =} ClassificationPartitionedLinear (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate a linear binary classifier.
    ##
    ## @code{@var{obj} = ClassificationPartitionedLinear (@var{X}, @var{Y})}
    ## partitions the data into ten stratified folds and fits a
    ## @qcode{ClassificationLinear} to each.
    ##
    ## @code{@var{obj} = ClassificationPartitionedLinear (@dots{},
    ## @var{name}, @var{value})} takes one of @qcode{'KFold'},
    ## @qcode{'Holdout'}, @qcode{'Leaveout'} and @qcode{'CVPartition'} to say
    ## how to partition, and any option @code{ClassificationLinear} takes to
    ## say how to fit.  @qcode{'CrossVal'} is accepted and has no effect
    ## here, this class being cross-validated by construction.
    ##
    ## The classes, the prior and the cost are resolved once over the whole
    ## data and handed to every fold.  Anything left as @qcode{'auto'} is
    ## not: each fold resolves @qcode{'Lambda'} against its own training
    ## rows, so ten folds of a hundred observations each get one ninetieth
    ## rather than one hundredth.  Both are MATLAB's behaviour, measured.
    ##
    ## @seealso{fitclinear, ClassificationLinear}
    ## @end deftypefn
    function this = ClassificationPartitionedLinear (X, Y, varargin)

      if (nargin < 2)
        error (strcat ("ClassificationPartitionedLinear: too few input", ...
                       " arguments."));
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationPartitionedLinear: optional", ...
                       " arguments must be given in Name-Value pairs."));
      endif

      ## Split the argument list three ways: what says how to partition,
      ## what the parent owns, and what each fold is fitted with.
      [P, args] = partitionedArgs (varargin, ...
                                   'ClassificationPartitionedLinear');

      F = classFrame (X, Y, P.ClassNames, P.Prior, P.Cost, P.Weights, ...
                      'ClassificationPartitionedLinear');

      [part, args] = cvPartitionOf (args, F.Y, F.n, ...
                                    'ClassificationPartitionedLinear');

      ## The score transform stays with the fold models: a logistic fold
      ## reports posteriors and the parent reports no transform of its own.
      ## Measured on R2024a, where a cross-validated logistic model has
      ## ScoreTransform 'none' while every Trained@{k@} has 'logit', and the
      ## assembled scores nonetheless sum to one.
      fargs = [args, {'ClassNames', F.ClassNames, 'Prior', F.Prior, ...
                      'Cost', F.Cost}];
      if (! isempty (P.ScoreTransform))
        fargs = [fargs, {'ScoreTransform', P.ScoreTransform}];
      endif
      ## Built field by field rather than with one struct call: a cell
      ## array response passed to struct () makes a struct array of that
      ## size instead of a scalar struct holding the cell.
      G = struct ();
      G.X = F.X;
      G.Y = F.Y;
      G.Weights = F.Weights;
      this.Trained = foldModels ('ClassificationLinear', G, part, fargs);

      this.ClassNames = F.ClassNames;
      this.Cost = F.Cost;
      this.Prior = F.Prior;
      this.NumObservations = F.n;
      this.Y = F.Y;
      this.W = F.W;
      this.X_ = F.X;
      this.KFold = part.NumTestSets;
      this.Partition = part;
      this.NumLambda_ = numel (this.Trained{1}.Lambda);

      this.PredictorNames = this.Trained{1}.PredictorNames;
      this.CategoricalPredictors = this.Trained{1}.CategoricalPredictors;
      this.ResponseName = this.Trained{1}.ResponseName;

      ## The parent applies nothing further, the folds having applied
      ## whatever transform they carry.  Its property reports that plainly
      ## rather than restating the folds'.
      [this.STfun, this.ScoreTransform] = ...
              parseScoreTransform ('none', ...
                                   'ClassificationPartitionedLinear');

      this.ModelParameters = struct ('Type', 'classification', ...
                                     'Method', 'PartitionedLinear', ...
                                     'LearnerTemplates', 'Linear', ...
                                     'NLearn', this.KFold);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedLinear} {@var{labels} =} kfoldPredict (@var{obj})
    ## @deftypefnx {ClassificationPartitionedLinear} {[@var{labels}, @var{scores}] =} kfoldPredict (@var{obj})
    ##
    ## Out-of-fold class of every observation.
    ##
    ## Each observation is classified by the fold that held it out, so the
    ## labels are out-of-sample.  An observation that no fold held out, which
    ## under a holdout partition is most of them, comes back missing rather
    ## than classified, and its scores come back @qcode{NaN}.
    ##
    ## With @math{L} regularization strengths @var{labels} has one column per
    ## strength and @var{scores} is @math{Nx2xL}.
    ##
    ## @end deftypefn
    function [labels, scores] = kfoldPredict (this)

      [labels, raw] = kfoldScores (this.Trained, this.Partition, this.X_, ...
                                   this.Y, this.ClassNames, ...
                                   this.NumLambda_);
      if (nargout > 1)
        scores = transformScores (this, raw);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedLinear} {@var{m} =} kfoldMargin (@var{obj})
    ##
    ## Out-of-fold classification margin of every observation.
    ##
    ## The score the out-of-fold model gives the true class, less the score
    ## it gives the other one.  An observation no fold held out comes back
    ## @qcode{NaN}.  With @math{L} regularization strengths @var{m} has one
    ## column per strength.
    ##
    ## @end deftypefn
    function m = kfoldMargin (this)

      [~, raw] = kfoldScores (this.Trained, this.Partition, this.X_, ...
                              this.Y, this.ClassNames, this.NumLambda_);
      s = transformScores (this, raw);
      gY = binaryIndex (this.ClassNames, this.Y, ...
                        'ClassificationPartitionedLinear', 'kfoldMargin');
      m = marginsOf (s, gY, this.NumLambda_);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedLinear} {@var{e} =} kfoldEdge (@var{obj})
    ## @deftypefnx {ClassificationPartitionedLinear} {@var{e} =} kfoldEdge (@dots{}, @var{name}, @var{value})
    ##
    ## Weighted mean of the out-of-fold classification margins.
    ##
    ## @code{@var{e} = kfoldEdge (@dots{}, @var{name}, @var{value})} takes
    ## @qcode{'Folds'}, a subset of the folds to average over, and
    ## @qcode{'Mode'}, either @qcode{'average'}, the default, or
    ## @qcode{'individual'}, which returns one row per fold instead.
    ##
    ## @end deftypefn
    function e = kfoldEdge (this, varargin)

      O = kfoldOpts (varargin, {}, 'ClassificationPartitionedLinear', ...
                     'kfoldEdge', this.KFold);
      [~, raw] = kfoldScores (this.Trained, this.Partition, this.X_, ...
                              this.Y, this.ClassNames, this.NumLambda_);
      s = transformScores (this, raw);
      gY = binaryIndex (this.ClassNames, this.Y, ...
                        'ClassificationPartitionedLinear', 'kfoldEdge');
      m = marginsOf (s, gY, this.NumLambda_);

      sets = foldSets (this.Partition, O.Folds, O.Mode, ...
                       this.NumObservations);
      e = zeros (numel (sets), this.NumLambda_);
      for i = 1:numel (sets)
        idx = sets{i};
        w = this.W(idx);
        w = w / sum (w);
        e(i,:) = sum (w .* m(idx,:), 1);
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedLinear} {@var{l} =} kfoldLoss (@var{obj})
    ## @deftypefnx {ClassificationPartitionedLinear} {@var{l} =} kfoldLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-fold classification loss.
    ##
    ## @code{@var{l} = kfoldLoss (@var{obj})} returns the out-of-fold
    ## misclassification rate.
    ##
    ## @code{@var{l} = kfoldLoss (@dots{}, @var{name}, @var{value})} takes
    ## @qcode{'LossFun'}, one of @qcode{'binodeviance'},
    ## @qcode{'classifcost'}, @qcode{'classiferror'}, @qcode{'exponential'},
    ## @qcode{'hinge'}, @qcode{'logit'}, @qcode{'mincost'} and
    ## @qcode{'quadratic'}; @qcode{'Folds'}; and @qcode{'Mode'}.
    ##
    ## @end deftypefn
    function l = kfoldLoss (this, varargin)

      valid = {'binodeviance', 'classifcost', 'classiferror', ...
               'exponential', 'hinge', 'logit', 'mincost', 'quadratic'};
      O = kfoldOpts (varargin, valid, 'ClassificationPartitionedLinear', ...
                     'kfoldLoss', this.KFold);
      if (isempty (O.LossFun))
        O.LossFun = 'classiferror';
      endif

      [~, raw] = kfoldScores (this.Trained, this.Partition, this.X_, ...
                              this.Y, this.ClassNames, this.NumLambda_);
      s = transformScores (this, raw);
      gY = binaryIndex (this.ClassNames, this.Y, ...
                        'ClassificationPartitionedLinear', 'kfoldLoss');

      sets = foldSets (this.Partition, O.Folds, O.Mode, ...
                       this.NumObservations);
      l = zeros (numel (sets), this.NumLambda_);
      for i = 1:numel (sets)
        idx = sets{i};
        w = this.W(idx);
        w = w / sum (w);
        for k = 1:this.NumLambda_
          if (this.NumLambda_ == 1)
            sk = s(idx,:);
          else
            sk = s(idx,:,k);
          endif
          l(i,k) = classificationLoss (O.LossFun, sk, gY(idx), w, this.Cost);
        endfor
      endfor

    endfunction

  endmethods

  methods (Access = public, Hidden)

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        printf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      printf ("\n  ClassificationPartitionedLinear\n\n");
      printf ("%+25s: '%s'\n", 'CrossValidatedModel', ...
              this.CrossValidatedModel);
      printf ("%+25s: %s\n", 'ClassNames', ...
              classNameListing (this.ClassNames));
      printf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      printf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      printf ("%+25s: %d\n", 'KFold', this.KFold);
      printf ("\n");
    endfunction

  endmethods

  methods (Access = private)

    ## Apply the parent's score transform to the assembled scores, whatever
    ## shape they came back in.
    function s = transformScores (this, raw)

      if (this.NumLambda_ == 1)
        s = this.STfun (raw);
      else
        s = raw;
        for k = 1:this.NumLambda_
          s(:,:,k) = this.STfun (raw(:,:,k));
        endfor
      endif

    endfunction

  endmethods

endclassdef

%!demo
%! ## Cross-validate a linear classifier on the two overlapping iris
%! ## species and read the out-of-sample error rate.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'KFold', 5)
%! outOfSample = kfoldLoss (CVMdl)

%!test
%! ## The model reports the surface MATLAB reports
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'KFold', 5);
%! assert_equal (class (CVMdl), 'ClassificationPartitionedLinear');
%! assert_equal (CVMdl.CrossValidatedModel, 'Linear');
%! assert_equal (CVMdl.KFold, 5);
%! assert_equal (CVMdl.NumObservations, 100);
%! assert_equal (CVMdl.ClassNames, {'versicolor'; 'virginica'});
%! assert_equal (CVMdl.Prior, [0.5, 0.5]);
%! assert_equal (CVMdl.Cost, [0, 1; 1, 0]);
%! assert_equal (CVMdl.ScoreTransform, 'none');
%! assert_equal (CVMdl.ResponseName, 'Y');
%! assert_equal (CVMdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (size (CVMdl.Trained), [5, 1]);
%! assert_equal (class (CVMdl.Trained{1}), 'ClassificationLinear');
%! assert_equal (sum (CVMdl.W), 1, 1e-12);

%!test
%! ## The properties are the ones MATLAB lists, in its order
%! load fisheriris
%! CVMdl = ClassificationPartitionedLinear (meas(51:end,:), species(51:end));
%! assert_equal (properties (CVMdl), ...
%!               {'ClassNames'; 'Cost'; 'Prior'; 'ScoreTransform'; ...
%!                'CrossValidatedModel'; 'NumObservations'; 'Y'; 'W'; ...
%!                'PredictorNames'; 'CategoricalPredictors'; ...
%!                'ResponseName'; 'Trained'; 'KFold'; 'Partition'; ...
%!                'ModelParameters'});

%!test
%! ## The default is ten stratified folds
%! load fisheriris
%! CVMdl = ClassificationPartitionedLinear (meas(51:end,:), species(51:end));
%! assert_equal (CVMdl.KFold, 10);
%! assert_equal (numel (CVMdl.Trained), 10);

%!test
%! ## Each fold resolves 'Lambda' against its own training rows, not the
%! ## parent's count: eighty of a hundred observations train each fold of
%! ## five, so the strength is one eightieth.  This is R2024a's number.
%! load fisheriris
%! CVMdl = ClassificationPartitionedLinear (meas(51:end,:), ...
%!                                          species(51:end), 'KFold', 5);
%! assert_equal (CVMdl.Trained{1}.Lambda, 1 / 80, 1e-15);
%! assert_equal (CVMdl.Trained{3}.Lambda, 1 / 80, 1e-15);

%!test
%! ## The class names, the prior and the cost are the parent's and are
%! ## handed to every fold, so an unbalanced problem does not give each fold
%! ## a prior of its own.  Also R2024a's behaviour.
%! load fisheriris
%! X = meas([51:100, 101:120],:);
%! Y = species([51:100, 101:120]);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'KFold', 5);
%! assert_equal (CVMdl.Prior, [50/70, 20/70], 1e-12);
%! assert_equal (CVMdl.Trained{1}.Prior, CVMdl.Prior);
%! assert_equal (CVMdl.Trained{4}.Prior, CVMdl.Prior);
%! assert_equal (CVMdl.Trained{1}.ClassNames, CVMdl.ClassNames);

%!test
%! ## kfoldPredict predicts each observation with the fold that held it out.
%! ## Doing the same partition by hand must give the same answers, which is
%! ## what pins the assembly rather than the fit.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! part = cvpartition (Y, 'KFold', 4);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'CVPartition', part);
%! [label, score] = kfoldPredict (CVMdl);
%! byhand = cell (100, 1);
%! for k = 1:4
%!   tr = training (part, k);
%!   te = test (part, k);
%!   m = ClassificationLinear (X(tr,:), Y(tr), 'ClassNames', ...
%!                             CVMdl.ClassNames, 'Prior', CVMdl.Prior, ...
%!                             'Cost', CVMdl.Cost);
%!   byhand(te) = predict (m, X(te,:));
%! endfor
%! assert_equal (label, byhand);
%! assert_equal (size (score), [100, 2]);

%!test
%! ## The margin is the out-of-fold score of the true class less the other,
%! ## and the edge its weighted mean
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'KFold', 5);
%! [~, score] = kfoldPredict (CVMdl);
%! m = kfoldMargin (CVMdl);
%! virg = strcmp (Y, 'virginica');
%! expect = score(:,1) - score(:,2);
%! expect(virg) = score(virg,2) - score(virg,1);
%! assert_equal (m, expect, 1e-12);
%! assert_equal (mean (m > 0), 1 - kfoldLoss (CVMdl), 1e-12);

%!test
%! ## 'Mode', 'individual' gives one row per fold
%! load fisheriris
%! CVMdl = ClassificationPartitionedLinear (meas(51:end,:), ...
%!                                          species(51:end), 'KFold', 5);
%! assert_equal (size (kfoldLoss (CVMdl, 'Mode', 'individual')), [5, 1]);
%! assert_equal (size (kfoldEdge (CVMdl, 'Mode', 'individual')), [5, 1]);

%!test
%! ## Averaging pools the observations rather than averaging the per-fold
%! ## values, which is not the same thing once the folds differ in size.
%! ## The pooled reading is the one R2024a returns.
%! load fisheriris
%! X = meas([51:100, 101:143],:);
%! Y = species([51:100, 101:143]);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'KFold', 4);
%! [~, score] = kfoldPredict (CVMdl);
%! gY = 1 + strcmp (Y, 'virginica');
%! wrong = mean (score(sub2ind (size (score), (1:93)', gY)) ...
%!               <= score(sub2ind (size (score), (1:93)', 3 - gY)));
%! assert_equal (kfoldLoss (CVMdl), wrong, 1e-12);

%!test
%! ## 'Folds' reports over the observations of the folds it names
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'KFold', 5);
%! two = kfoldLoss (CVMdl, 'Folds', [1, 2]);
%! idx = test (CVMdl.Partition, 1) | test (CVMdl.Partition, 2);
%! [label, ~] = kfoldPredict (CVMdl);
%! assert_equal (two, mean (! strcmp (label(idx), Y(idx))), 1e-12);

%!test
%! ## Every classification loss is offered, and each reads the out-of-fold
%! ## score of the true class
%! load fisheriris
%! CVMdl = ClassificationPartitionedLinear (meas(51:end,:), ...
%!                                          species(51:end), 'KFold', 5);
%! for f = {'binodeviance', 'classifcost', 'classiferror', 'exponential', ...
%!          'hinge', 'logit', 'mincost', 'quadratic'}
%!   assert_equal (isfinite (kfoldLoss (CVMdl, 'LossFun', f{1})), true);
%! endfor

%!test
%! ## A logistic fit reports posteriors, and the transform that produces
%! ## them stays with the fold models: the cross-validated model reports
%! ## none of its own, which is R2024a's arrangement and not the reverse
%! load fisheriris
%! CVMdl = ClassificationPartitionedLinear (meas(51:end,:), ...
%!                                          species(51:end), 'KFold', 5, ...
%!                                          'Learner', 'logistic');
%! assert_equal (CVMdl.ScoreTransform, 'none');
%! assert_equal (CVMdl.Trained{1}.ScoreTransform, 'logit');
%! [~, score] = kfoldPredict (CVMdl);
%! assert_equal (sum (score, 2), ones (100, 1), 1e-12);

%!test
%! ## A whole regularization path gives one column per strength everywhere
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'KFold', 5, ...
%!                                          'Lambda', [0.001, 0.01, 0.1]);
%! [label, score] = kfoldPredict (CVMdl);
%! assert_equal (size (label), [100, 3]);
%! assert_equal (size (score), [100, 2, 3]);
%! assert_equal (size (kfoldMargin (CVMdl)), [100, 3]);
%! assert_equal (size (kfoldLoss (CVMdl)), [1, 3]);
%! assert_equal (size (kfoldEdge (CVMdl)), [1, 3]);
%! assert_equal (size (kfoldLoss (CVMdl, 'Mode', 'individual')), [5, 3]);

%!test
%! ## An observation that no fold held out is not classified: under a
%! ## holdout partition that is the training rows, and they come back
%! ## missing rather than carrying a class
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'Holdout', 0.3);
%! assert_equal (CVMdl.KFold, 1);
%! [label, score] = kfoldPredict (CVMdl);
%! assert_equal (sum (cellfun (@isempty, label)), 70);
%! assert_equal (sum (isnan (score(:,1))), 70);
%! assert_equal (isfinite (kfoldLoss (CVMdl)), true);

%!test
%! ## Leave-one-out gives as many folds as there are observations
%! load fisheriris
%! X = meas([51:70, 101:120],:);
%! Y = species([51:70, 101:120]);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'Leaveout', 'on');
%! assert_equal (CVMdl.KFold, 40);
%! assert_equal (numel (kfoldPredict (CVMdl)), 40);

%!test
%! ## A row with a missing predictor is dropped before the partition, so
%! ## every index below refers to the same set of observations
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! X(3,2) = NaN;
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'KFold', 5);
%! assert_equal (CVMdl.NumObservations, 99);
%! assert_equal (numel (CVMdl.Y), 99);
%! assert_equal (CVMdl.Partition.NumObservations, 99);

%!test
%! ## Observation weights reach both the folds and the reported prior
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedLinear (X, Y, 'KFold', 5, ...
%!                                          'Weights', (1:100)');
%! assert_equal (CVMdl.Prior, [0.252475247524752, 0.747524752475248], 1e-12);
%! assert_equal (CVMdl.Trained{1}.Prior, CVMdl.Prior);
%! assert_equal (sum (CVMdl.W), 1, 1e-12);

%!test
%! ## A character matrix response carries through cross-validation: the
%! ## class names, the fold models, the assembled labels and every kfold
%! ## method answer as they do for the equivalent cell array
%! load fisheriris
%! X = meas(51:end,:);
%! Yc = species(51:end);
%! Ym = char (Yc);
%! part = cvpartition (Yc, 'KFold', 4);
%! CVc = ClassificationPartitionedLinear (X, Yc, 'CVPartition', part);
%! CVm = ClassificationPartitionedLinear (X, Ym, 'CVPartition', part);
%! assert_equal (cellstr (CVm.ClassNames), CVc.ClassNames);
%! assert_equal (size (CVm.ClassNames), [2, 10]);
%! assert_equal (cellstr (CVm.Trained{1}.ClassNames), CVc.ClassNames);
%! assert_equal (cellstr (kfoldPredict (CVm)), kfoldPredict (CVc));
%! assert_equal (kfoldMargin (CVm), kfoldMargin (CVc));
%! assert_equal (kfoldEdge (CVm), kfoldEdge (CVc));
%! assert_equal (kfoldLoss (CVm), kfoldLoss (CVc));
%! assert_equal (kfoldLoss (CVm, 'LossFun', 'hinge'), ...
%!               kfoldLoss (CVc, 'LossFun', 'hinge'));

%!test
%! ## Names of unequal length are padded by the character matrix and the
%! ## padding is not part of the name, through the cross-validated path too
%! load fisheriris
%! X = meas(51:end,:);
%! Y = [repmat({'ab'}, 50, 1); repmat({'abcd'}, 50, 1)];
%! part = cvpartition (Y, 'KFold', 4);
%! CVc = ClassificationPartitionedLinear (X, Y, 'CVPartition', part);
%! CVm = ClassificationPartitionedLinear (X, char (Y), 'CVPartition', part);
%! assert_equal (cellstr (CVm.ClassNames), CVc.ClassNames);
%! assert_equal (cellstr (kfoldPredict (CVm)), kfoldPredict (CVc));

%!test
%! ## A row dropped for a missing predictor is dropped before the partition,
%! ## through the character path, which is the case that exercises the row
%! ## indexing rather than the label comparison
%! load fisheriris
%! X = meas(51:end,:);
%! X(7,3) = NaN;
%! CVMdl = ClassificationPartitionedLinear (X, char (species(51:end)), ...
%!                                          'KFold', 4);
%! assert_equal (CVMdl.NumObservations, 99);
%! assert_equal (size (CVMdl.Y), [99, 10]);
%! assert_equal (numel (kfoldPredict (CVMdl)) / 10, 99);

## Test input validation
%!error<ClassificationPartitionedLinear: too few input arguments.> ...
%! ClassificationPartitionedLinear (ones (10, 2))
%!error<ClassificationPartitionedLinear: optional arguments must be given in Name-Value pairs.> ...
%! ClassificationPartitionedLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'KFold')
%!error<ClassificationPartitionedLinear: 'KFold' must be an integer value greater than 1.> ...
%! ClassificationPartitionedLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'KFold', 1)
%!error<ClassificationPartitionedLinear: 'Holdout' must be a numeric value between 0 and 1.> ...
%! ClassificationPartitionedLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'Holdout', 1.5)
%!error<ClassificationPartitionedLinear: 'Leaveout' must be either 'on' or 'off'.> ...
%! ClassificationPartitionedLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'Leaveout', 1)
%!error<ClassificationPartitionedLinear: 'CVPartition' must be a cvpartition object.> ...
%! ClassificationPartitionedLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'CVPartition', 3)
%!error<ClassificationPartitionedLinear: you can use only one of 'KFold', 'Holdout', 'Leaveout', or 'CVPartition' options.> ...
%! ClassificationPartitionedLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'KFold', 2, 'Holdout', 0.2)
%!error<ClassificationPartitionedLinear: 'CVPartition' must be built over the same number of observations as the data.> ...
%! ClassificationPartitionedLinear (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'CVPartition', cvpartition (20, 'KFold', 2))
%!error<ClassificationPartitionedLinear: Y must name exactly two classes, this being a binary model; use fitcecoc for more than two.> ...
%! ClassificationPartitionedLinear (ones (9, 2), ...
%!                     [1; 1; 1; 2; 2; 2; 3; 3; 3], 'KFold', 3)
%!error<ClassificationPartitionedLinear.kfoldLoss: 'Mode' must be either 'average' or 'individual'.> ...
%! kfoldLoss (ClassificationPartitionedLinear (ones (10, 2), ...
%!                     [ones(5,1); 2*ones(5,1)], 'KFold', 2), 'Mode', 'each')
%!error<ClassificationPartitionedLinear.kfoldLoss: 'Folds' must hold integers between 1 and 2.> ...
%! kfoldLoss (ClassificationPartitionedLinear (ones (10, 2), ...
%!                     [ones(5,1); 2*ones(5,1)], 'KFold', 2), 'Folds', 7)
%!error<ClassificationPartitionedLinear.kfoldLoss: 'LossFun' must be 'binodeviance', 'classifcost', 'classiferror', 'exponential', 'hinge', 'logit', 'mincost', or 'quadratic'.> ...
%! kfoldLoss (ClassificationPartitionedLinear (ones (10, 2), ...
%!                     [ones(5,1); 2*ones(5,1)], 'KFold', 2), 'LossFun', 'mse')
%!error<ClassificationPartitionedLinear.kfoldEdge: invalid parameter name in optional pair arguments.> ...
%! kfoldEdge (ClassificationPartitionedLinear (ones (10, 2), ...
%!                     [ones(5,1); 2*ones(5,1)], 'KFold', 2), 'LossFun', ...
%!                     'hinge')
