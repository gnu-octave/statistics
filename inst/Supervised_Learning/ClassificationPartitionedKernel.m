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
## @deftp {statistics} ClassificationPartitionedKernel
##
## Cross-validated Gaussian kernel binary classifier.
##
## A @qcode{ClassificationPartitionedKernel} object holds one
## @qcode{ClassificationKernel} per fold of a partition, each fitted to the
## observations the fold trains on.  Every @code{kfold} method predicts each
## observation with the fold that held it @emph{out}, so the estimate it
## returns is an out-of-sample one.
##
## A @qcode{ClassificationKernel} stores no copy of its training data and so
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
## @seealso{fitclinear, ClassificationKernel, ClassificationPartitionedKernel}
## @end deftp

classdef ClassificationPartitionedKernel

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} ClassNames
    ##
    ## Names of the two classes
    ##
    ## A column of the same type as the response, shared by every fold.
    ## This property is read-only.
    ##
    ## @end deftp
    ClassNames             = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} Cost
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
    ## @deftp {ClassificationPartitionedKernel} {property} Prior
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

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} ScoreTransform
    ##
    ## Transformation applied to the predicted scores
    ##
    ## A character vector naming a transformation, or the text of the
    ## function handle that was supplied, which may be assigned after the
    ## model is built.  It is applied once to the assembled scores and is
    ## not handed to the folds.  A transform the learner implies, as
    ## @qcode{'logistic'} implies @qcode{'logit'}, does stay with the folds,
    ## and this one is then applied on top of it.
    ##
    ## @end deftp
    ScoreTransform         = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} CrossValidatedModel
    ##
    ## Name of the model that was cross-validated
    ##
    ## Always @qcode{'Linear'}, the short name MATLAB uses.  This property is
    ## read-only.
    ##
    ## @end deftp
    CrossValidatedModel    = 'Kernel';

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} NumObservations
    ##
    ## Number of observations the partition covers
    ##
    ## A positive integer scalar, counting the rows that survived the removal
    ## of missing values.  This property is read-only.
    ##
    ## @end deftp
    NumObservations        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} Y
    ##
    ## Response of the retained observations
    ##
    ## In the type it was supplied in.  This property is read-only.
    ##
    ## @end deftp
    Y                      = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} W
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
    ## @deftp {ClassificationPartitionedKernel} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames         = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A row vector of column indices, empty when every predictor is
    ## numeric.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} ResponseName
    ##
    ## Name of the response
    ##
    ## A character vector, defaulting to @qcode{'Y'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    ResponseName           = 'Y';

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} Trained
    ##
    ## The models fitted to the folds
    ##
    ## A cell column with one @qcode{ClassificationKernel} per fold, each
    ## fitted to the observations its fold trains on.  This property is
    ## read-only.
    ##
    ## @end deftp
    Trained                = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} KFold
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
    ## @deftp {ClassificationPartitionedKernel} {property} Partition
    ##
    ## The partition itself
    ##
    ## A @code{cvpartition} object over the retained observations.  This
    ## property is read-only.
    ##
    ## @end deftp
    Partition              = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedKernel} {property} ModelParameters
    ##
    ## What was cross-validated, and how
    ##
    ## A structure holding the parameters the folds were fitted with, carried
    ## through from the learner that was cross validated, beside
    ## @qcode{NLearn}, the number of folds, and the @qcode{Version},
    ## @qcode{Method} and @qcode{Type} tags of this class, with
    ## @qcode{LearnerTemplates} naming the backing.  The
    ## learner's own tags are replaced rather than kept, so a cross-validated
    ## SVM reports @qcode{Method} as @qcode{'PartitionedKernel'} and not
    ## @qcode{'SVM'}.
    ##
    ## @strong{Deviation from MATLAB.}  MATLAB reports the parameter record of
    ## the cross-validation @emph{ensemble} here rather than of the learner,
    ## so it says nothing at all about how the folds were fitted: of its
    ## eighteen fields only the fold count, its partitioner and a fit template
    ## carry anything, and the rest are boosting settings left inert.  Nor can
    ## the parameters be reached through the folds, a compact model carrying
    ## none in MATLAB.  This class reports the fit instead, which is strictly
    ## more than MATLAB offers, and everything MATLAB's record does carry is
    ## published here as the @qcode{KFold}, @qcode{Partition}, @qcode{X},
    ## @qcode{Y}, @qcode{W} and @qcode{CrossValidatedModel} properties.
    ##
    ## This property is read-only.
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
    ## @deftypefn  {ClassificationPartitionedKernel} {@var{obj} =} ClassificationPartitionedKernel (@var{X}, @var{Y})
    ## @deftypefnx {ClassificationPartitionedKernel} {@var{obj} =} ClassificationPartitionedKernel (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate a linear binary classifier.
    ##
    ## @code{@var{obj} = ClassificationPartitionedKernel (@var{X}, @var{Y})}
    ## partitions the data into ten stratified folds and fits a
    ## @qcode{ClassificationKernel} to each.
    ##
    ## @code{@var{obj} = ClassificationPartitionedKernel (@dots{},
    ## @var{name}, @var{value})} takes one of @qcode{'KFold'},
    ## @qcode{'Holdout'}, @qcode{'Leaveout'} and @qcode{'CVPartition'} to say
    ## how to partition, and any option @code{ClassificationKernel} takes to
    ## say how to fit.  @qcode{'CrossVal'} is accepted and has no effect
    ## here, this class being cross-validated by construction.
    ##
    ## The classes, the prior and the cost are resolved once over the whole
    ## data and handed to every fold.  Anything left as @qcode{'auto'} is
    ## not: each fold resolves @qcode{'Lambda'} and @qcode{'KernelScale'}
    ## against its own training rows, so ten folds of a hundred
    ## observations each get a @qcode{Lambda} of one ninetieth rather than
    ## one hundredth.  Both are MATLAB's behaviour, measured.
    ##
    ## @seealso{fitclinear, ClassificationKernel}
    ## @end deftypefn
    function this = ClassificationPartitionedKernel (X, Y, varargin)

      if (nargin < 2)
        error (strcat ("ClassificationPartitionedKernel: too few input", ...
                       " arguments."));
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationPartitionedKernel: optional", ...
                       " arguments must be given in Name-Value pairs."));
      endif

      ## Split the argument list three ways: what says how to partition,
      ## what the parent owns, and what each fold is fitted with.
      [P, args] = partitionedArgs (varargin, ...
                                   'ClassificationPartitionedKernel');

      F = classFrame (X, Y, P.ClassNames, P.Prior, P.Cost, P.Weights, ...
                      'ClassificationPartitionedKernel');

      [part, args] = cvPartitionOf (args, F.Y, F.n, ...
                                    'ClassificationPartitionedKernel');

      ## A transform asked for by name belongs to the parent, which applies
      ## it once to the assembled scores, and the folds are not given it.
      ## Measured on R2024a, where a cross-validated model fitted with a
      ## ScoreTransform reports it and leaves every Trained@{k@} at 'none'.
      ## A transform the learner implies is a different thing: a logistic
      ## fold reports posteriors and carries 'logit' of its own, which
      ## R2024a also does, and the parent's transform composes on top of it.
      fargs = [args, {'ClassNames', F.ClassNames, 'Prior', F.Prior, ...
                      'Cost', F.Cost}];
      ## Built field by field rather than with one struct call: a cell
      ## array response passed to struct () makes a struct array of that
      ## size instead of a scalar struct holding the cell.
      G = struct ();
      G.X = F.X;
      G.Y = F.Y;
      G.Weights = F.Weights;
      this.Trained = foldModels ('ClassificationKernel', G, part, fargs);

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

      if (isempty (P.ScoreTransform))
        P.ScoreTransform = 'none';
      endif
      [this.STfun, this.ScoreTransform] = ...
              parseScoreTransform (P.ScoreTransform, ...
                                   'ClassificationPartitionedKernel');

      ## The learner's parameters, under this class's own tags.  MATLAB
      ## reports an EnsembleParams here instead and so says nothing about
      ## the fit; see the ModelParameters property for the deviation.
      this.ModelParameters = partitionedModelParams (this.Trained{1}, ...
                               this.KFold, 'PartitionedKernel', ...
                               'classification', 'Kernel');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedKernel} {@var{labels} =} kfoldPredict (@var{obj})
    ## @deftypefnx {ClassificationPartitionedKernel} {[@var{labels}, @var{scores}] =} kfoldPredict (@var{obj})
    ##
    ## Out-of-fold class of every observation.
    ##
    ## Each observation is classified by the fold that held it out, so the
    ## labels are out-of-sample.  An observation that no fold held out, which
    ## under a holdout partition is most of them, comes back missing rather
    ## than classified, and its scores come back @qcode{NaN}.
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
    ## @deftypefn  {ClassificationPartitionedKernel} {@var{m} =} kfoldMargin (@var{obj})
    ##
    ## Out-of-fold classification margin of every observation.
    ##
    ## The score the out-of-fold model gives the true class, less the score
    ## it gives the other one.  An observation no fold held out comes back
    ## @qcode{NaN}.
    ##
    ## @end deftypefn
    function m = kfoldMargin (this)

      [~, raw] = kfoldScores (this.Trained, this.Partition, this.X_, ...
                              this.Y, this.ClassNames, this.NumLambda_);
      s = transformScores (this, raw);
      [gY, errmsg] = labelIndices (this.ClassNames, this.Y);
      if (! isempty (errmsg))
        error ("ClassificationPartitionedKernel.kfoldMargin: %s", errmsg);
      endif
      m = marginsOf (s, gY, this.NumLambda_);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedKernel} {@var{e} =} kfoldEdge (@var{obj})
    ## @deftypefnx {ClassificationPartitionedKernel} {@var{e} =} kfoldEdge (@dots{}, @var{name}, @var{value})
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

      O = kfoldOpts (varargin, {}, 'ClassificationPartitionedKernel', ...
                     'kfoldEdge', this.KFold);
      [~, raw] = kfoldScores (this.Trained, this.Partition, this.X_, ...
                              this.Y, this.ClassNames, this.NumLambda_);
      s = transformScores (this, raw);
      [gY, errmsg] = labelIndices (this.ClassNames, this.Y);
      if (! isempty (errmsg))
        error ("ClassificationPartitionedKernel.kfoldEdge: %s", errmsg);
      endif
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
    ## @deftypefn  {ClassificationPartitionedKernel} {@var{l} =} kfoldLoss (@var{obj})
    ## @deftypefnx {ClassificationPartitionedKernel} {@var{l} =} kfoldLoss (@dots{}, @var{name}, @var{value})
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
      O = kfoldOpts (varargin, valid, 'ClassificationPartitionedKernel', ...
                     'kfoldLoss', this.KFold);
      if (isempty (O.LossFun))
        O.LossFun = 'classiferror';
      endif

      [~, raw] = kfoldScores (this.Trained, this.Partition, this.X_, ...
                              this.Y, this.ClassNames, this.NumLambda_);
      s = transformScores (this, raw);
      [gY, errmsg] = labelIndices (this.ClassNames, this.Y);
      if (! isempty (errmsg))
        error ("ClassificationPartitionedKernel.kfoldLoss: %s", errmsg);
      endif

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

    function this = set.ScoreTransform (this, val)
      [f, nm] = parseScoreTransform (val, 'ClassificationPartitionedKernel');
      this.ScoreTransform = nm;
      this.STfun = f;
    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        printf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      printf ("\n  ClassificationPartitionedKernel\n\n");
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
%! ## Cross-validate a Gaussian kernel classifier on the two overlapping
%! ## iris species and read the out-of-sample error rate.
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedKernel (X, Y, 'KFold', 5)
%! outOfSample = kfoldLoss (CVMdl)

%!test
%! ## The model reports the surface MATLAB reports
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedKernel (X, Y, 'KFold', 5);
%! assert_equal (class (CVMdl), 'ClassificationPartitionedKernel');
%! assert_equal (CVMdl.CrossValidatedModel, 'Kernel');
%! assert_equal (CVMdl.KFold, 5);
%! assert_equal (CVMdl.NumObservations, 100);
%! assert_equal (CVMdl.ClassNames, {'versicolor'; 'virginica'});
%! assert_equal (CVMdl.Prior, [0.5, 0.5]);
%! assert_equal (class (CVMdl.Trained{1}), 'ClassificationKernel');
%! assert_equal (CVMdl.ModelParameters.Method, 'PartitionedKernel');
%! assert_equal (CVMdl.ModelParameters.LearnerTemplates, 'Kernel');
%! assert_equal (CVMdl.ModelParameters.NLearn, 5);

%!test
%! ## The properties are the ones MATLAB lists, in its order
%! load fisheriris
%! CVMdl = ClassificationPartitionedKernel (meas(51:end,:), species(51:end));
%! assert_equal (sort (properties (CVMdl)), ...
%!               sort ({'ClassNames'; 'Cost'; 'Prior'; 'ScoreTransform'; ...
%!                      'CrossValidatedModel'; 'NumObservations'; 'Y'; 'W'; ...
%!                      'PredictorNames'; 'CategoricalPredictors'; ...
%!                      'ResponseName'; 'Trained'; 'KFold'; 'Partition'; ...
%!                      'ModelParameters'}));

%!test
%! ## Each fold resolves 'Lambda' against its own training rows, which is
%! ## one eightieth of five folds over a hundred observations.  R2024a's
%! ## number.
%! load fisheriris
%! CVMdl = ClassificationPartitionedKernel (meas(51:end,:), ...
%!                                          species(51:end), 'KFold', 5);
%! assert_equal (CVMdl.Trained{1}.Lambda, 1 / 80, 1e-15);
%! assert_equal (CVMdl.Trained{1}.NumExpansionDimensions, 128);
%! assert_equal (CVMdl.Trained{1}.KernelScale, 1);

%!test
%! ## Every fold draws its own basis, so two folds hold different
%! ## expansions of the same kernel
%! load fisheriris
%! CVMdl = ClassificationPartitionedKernel (meas(51:end,:), ...
%!                                          species(51:end), 'KFold', 5);
%! X = meas(51:53,:);
%! [~, s1] = predict (CVMdl.Trained{1}, X);
%! [~, s2] = predict (CVMdl.Trained{2}, X);
%! assert_equal (isequal (s1, s2), false);

%!test
%! ## kfoldPredict classifies each observation with the fold that held it
%! ## out, and doing the same partition by hand agrees
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! part = cvpartition (Y, 'KFold', 4);
%! CVMdl = ClassificationPartitionedKernel (X, Y, 'CVPartition', part);
%! label = kfoldPredict (CVMdl);
%! byhand = cell (100, 1);
%! for k = 1:4
%!   te = test (part, k);
%!   byhand(te) = predict (CVMdl.Trained{k}, X(te,:));
%! endfor
%! assert_equal (label, byhand);

%!test
%! ## The fit separates the two species out of sample
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedKernel (X, Y, 'KFold', 5);
%! assert_equal (kfoldLoss (CVMdl) < 0.2, true);
%! assert_equal (kfoldEdge (CVMdl) > 0, true);
%! assert_equal (numel (kfoldMargin (CVMdl)), 100);

%!test
%! ## 'Mode', 'individual' gives one row per fold, and averaging pools the
%! ## observations rather than the fold values
%! load fisheriris
%! X = meas(51:end,:);
%! Y = species(51:end);
%! CVMdl = ClassificationPartitionedKernel (X, Y, 'KFold', 5);
%! assert_equal (size (kfoldLoss (CVMdl, 'Mode', 'individual')), [5, 1]);
%! label = kfoldPredict (CVMdl);
%! assert_equal (kfoldLoss (CVMdl), mean (! strcmp (label, Y)), 1e-12);

%!test
%! ## A logistic fit reports posteriors, and the transform stays with the
%! ## fold models as it does in MATLAB
%! load fisheriris
%! CVMdl = ClassificationPartitionedKernel (meas(51:end,:), ...
%!                                          species(51:end), 'KFold', 5, ...
%!                                          'Learner', 'logistic');
%! assert_equal (CVMdl.ScoreTransform, 'none');
%! assert_equal (CVMdl.Trained{1}.ScoreTransform, 'logit');
%! [~, score] = kfoldPredict (CVMdl);
%! assert_equal (sum (score, 2), ones (100, 1), 1e-12);

%!test
%! ## An observation that no fold held out is not classified
%! load fisheriris
%! CVMdl = ClassificationPartitionedKernel (meas(51:end,:), ...
%!                                          species(51:end), 'Holdout', 0.3);
%! assert_equal (CVMdl.KFold, 1);
%! [label, score] = kfoldPredict (CVMdl);
%! assert_equal (sum (cellfun (@isempty, label)), 70);
%! assert_equal (sum (isnan (score(:,1))), 70);

%!test
%! ## A character matrix response carries through cross-validation
%! load fisheriris
%! X = meas(51:end,:);
%! Yc = species(51:end);
%! part = cvpartition (Yc, 'KFold', 4);
%! CVm = ClassificationPartitionedKernel (X, char (Yc), 'CVPartition', part);
%! assert_equal (size (CVm.ClassNames), [2, 10]);
%! assert_equal (cellstr (CVm.ClassNames), {'versicolor'; 'virginica'});
%! assert_equal (size (kfoldPredict (CVm)), [100, 10]);
%! assert_equal (isfinite (kfoldLoss (CVm)), true);
%! assert_equal (isfinite (kfoldEdge (CVm)), true);

%!test
%! ## A transform asked for by name goes to the parent and not to the folds,
%! ## and is applied once to the assembled scores.  R2024a's arrangement.
%! ## The baseline is read from the same object with the transform switched
%! ## off, a second fit being no baseline at all here: the random feature
%! ## expansion differs between two fits of the same data.
%! load fisheriris
%! CVMdl = ClassificationPartitionedKernel (meas(51:end,:), species(51:end), ...
%!                                          'KFold', 5, ...
%!                                          'ScoreTransform', 'doublelogit');
%! assert_equal (CVMdl.ScoreTransform, 'doublelogit');
%! assert_equal (CVMdl.Trained{1}.ScoreTransform, 'none');
%! [~, s1] = kfoldPredict (CVMdl);
%! CVMdl.ScoreTransform = 'none';
%! [~, s0] = kfoldPredict (CVMdl);
%! assert_equal (s1, 1 ./ (1 + exp (-2 * s0)), 1e-12);

%!test
%! ## It can be assigned after the model is built, and reaches kfoldPredict
%! ## without being carried into the folds
%! load fisheriris
%! CVMdl = ClassificationPartitionedKernel (meas(51:end,:), species(51:end), ...
%!                                          'KFold', 5);
%! [~, s0] = kfoldPredict (CVMdl);
%! CVMdl.ScoreTransform = 'doublelogit';
%! [~, s1] = kfoldPredict (CVMdl);
%! assert_equal (CVMdl.Trained{1}.ScoreTransform, 'none');
%! assert_equal (s1, 1 ./ (1 + exp (-2 * s0)), 1e-12);

%!test
%! ## A transform the learner implies stays with the folds, and an assigned
%! ## one is applied on top of it rather than replacing it.  Measured on
%! ## R2024a, where the folds keep 'logit' and the parent's transform
%! ## composes.
%! load fisheriris
%! CVMdl = ClassificationPartitionedKernel (meas(51:end,:), species(51:end), ...
%!                                          'KFold', 5, ...
%!                                          'Learner', 'logistic');
%! [~, s0] = kfoldPredict (CVMdl);
%! CVMdl.ScoreTransform = 'doublelogit';
%! [~, s1] = kfoldPredict (CVMdl);
%! assert_equal (CVMdl.Trained{1}.ScoreTransform, 'logit');
%! assert_equal (s1, 1 ./ (1 + exp (-2 * s0)), 1e-12);

%!test
%! ## 'none' is the identity, so assigning it transforms nothing
%! load fisheriris
%! CVMdl = ClassificationPartitionedKernel (meas(51:end,:), species(51:end), ...
%!                                          'KFold', 5);
%! [~, s0] = kfoldPredict (CVMdl);
%! CVMdl.ScoreTransform = 'none';
%! [~, s1] = kfoldPredict (CVMdl);
%! assert_equal (s1, s0);

%!error<ClassificationPartitionedKernel: unrecognized 'ScoreTransform' function.> ...
%! load fisheriris
%! CVMdl = ClassificationPartitionedKernel (meas(51:end,:), species(51:end), ...
%!                                          'KFold', 5);
%! CVMdl.ScoreTransform = 'nosuchtransform';

## Test input validation
%!error<ClassificationPartitionedKernel: too few input arguments.> ...
%! ClassificationPartitionedKernel (ones (10, 2))
%!error<ClassificationPartitionedKernel: optional arguments must be given in Name-Value pairs.> ...
%! ClassificationPartitionedKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'KFold')
%!error<ClassificationPartitionedKernel: 'KFold' must be an integer value greater than 1.> ...
%! ClassificationPartitionedKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'KFold', 0)
%!error<ClassificationPartitionedKernel: you can use only one of 'KFold', 'Holdout', 'Leaveout', or 'CVPartition' options.> ...
%! ClassificationPartitionedKernel (ones (10, 2), [ones(5,1); 2*ones(5,1)], ...
%!                     'KFold', 2, 'Leaveout', 'on')
%!error<ClassificationPartitionedKernel.kfoldLoss: 'LossFun' must be 'binodeviance', 'classifcost', 'classiferror', 'exponential', 'hinge', 'logit', 'mincost', or 'quadratic'.> ...
%! kfoldLoss (ClassificationPartitionedKernel (ones (10, 2), ...
%!                     [ones(5,1); 2*ones(5,1)], 'KFold', 2), 'LossFun', 'mse')
%!error<ClassificationPartitionedKernel.kfoldEdge: 'Folds' must hold integers between 1 and 2.> ...
%! kfoldEdge (ClassificationPartitionedKernel (ones (10, 2), ...
%!                     [ones(5,1); 2*ones(5,1)], 'KFold', 2), 'Folds', 0)

## ModelParameters carries the learner's parameters beside the tags this
## class reports for itself.
%!test
%! load fisheriris
%! CVMdl = fitckernel (meas, strcmp (species, 'setosa'), 'KFold', 3);
%! MP = CVMdl.ModelParameters;
%! assert_equal (MP.Method, 'PartitionedKernel');
%! assert_equal (MP.LearnerTemplates, 'Kernel');
%! assert_equal (MP.NLearn, 3);
%! assert_equal (MP.Learner, 'svm');
%! assert_equal (MP.BlockSize, 4000);
