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
## @deftp {statistics} RegressionPartitionedLinear
##
## Cross-validated linear regression model.
##
## A @qcode{RegressionPartitionedLinear} object holds one
## @qcode{RegressionLinear} per fold of a partition, each fitted to the
## observations the fold trains on.  @code{kfoldPredict} predicts each
## observation with the fold that held it @emph{out}, so what it returns is
## an out-of-sample prediction.
##
## A @qcode{RegressionLinear} stores no copy of its training data and so has
## no resubstitution methods and no @code{compact} form.  This class is what
## takes their place: cross-validation is the way a linear model is asked
## how it would do on data it has not seen.
##
## When the fold models carry a whole regularization path, both methods
## return one column per strength, in the order of the @qcode{'Lambda'} that
## was asked for.
##
## Create one with @code{fitrlinear} and a cross-validation option, or
## directly.
##
## @seealso{fitrlinear, RegressionLinear, RegressionPartitionedKernel}
## @end deftp

classdef RegressionPartitionedLinear

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} ResponseTransform
    ##
    ## Transformation applied to the predicted response
    ##
    ## A character vector, or the text of the function handle that was
    ## supplied, which may be assigned after the model is built.  The fold
    ## models carry no transform of their own; this one is applied once to
    ## the assembled predictions.
    ##
    ## @end deftp
    ResponseTransform      = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} CrossValidatedModel
    ##
    ## Name of the model that was cross-validated
    ##
    ## Always @qcode{'Linear'}, the short name MATLAB uses.  This property is
    ## read-only.
    ##
    ## @end deftp
    CrossValidatedModel    = 'Linear';

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} NumObservations
    ##
    ## Number of observations the partition covers
    ##
    ## A positive integer scalar, counting the rows that survived the removal
    ## of missing values.  This property is read-only.
    ##
    ## @end deftp
    NumObservations        = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} Y
    ##
    ## Response of the retained observations
    ##
    ## An @math{Nx1} numeric vector.  This property is read-only.
    ##
    ## @end deftp
    Y                      = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} W
    ##
    ## Observation weights
    ##
    ## An @math{Nx1} numeric vector summing to one.  This property is
    ## read-only.
    ##
    ## @end deftp
    W                      = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames         = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A row vector of column indices, empty when every predictor is
    ## numeric.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors  = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} ResponseName
    ##
    ## Name of the response
    ##
    ## A character vector, defaulting to @qcode{'Y'}.  This property is
    ## read-only.
    ##
    ## @end deftp
    ResponseName           = 'Y';

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} Trained
    ##
    ## The models fitted to the folds
    ##
    ## A cell column with one @qcode{RegressionLinear} per fold, each fitted
    ## to the observations its fold trains on.  This property is read-only.
    ##
    ## @end deftp
    Trained                = {};

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} KFold
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
    ## @deftp {RegressionPartitionedLinear} {property} Partition
    ##
    ## The partition itself
    ##
    ## A @code{cvpartition} object over the retained observations.  This
    ## property is read-only.
    ##
    ## @end deftp
    Partition              = [];

    ## -*- texinfo -*-
    ## @deftp {RegressionPartitionedLinear} {property} ModelParameters
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

    ## The callable behind ResponseTransform.
    RTfun                  = @(y) y;

    ## The predictors of the retained observations.  MATLAB does not report
    ## them and neither do we, but the kfold methods have to predict from
    ## something and the fold models hold no data of their own.
    X_                     = [];

    ## Number of regularization strengths the fold models carry.
    NumLambda_             = 1;

    ## Whether the fold models were fitted with an insensitive band, which
    ## decides whether kfoldLoss will offer the loss that reads it.
    HasEpsilon_            = false;

  endproperties

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionPartitionedLinear} {@var{obj} =} RegressionPartitionedLinear (@var{X}, @var{Y})
    ## @deftypefnx {RegressionPartitionedLinear} {@var{obj} =} RegressionPartitionedLinear (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate a linear regression model.
    ##
    ## @code{@var{obj} = RegressionPartitionedLinear (@var{X}, @var{Y})}
    ## partitions the data into ten folds and fits a
    ## @qcode{RegressionLinear} to each.
    ##
    ## @code{@var{obj} = RegressionPartitionedLinear (@dots{}, @var{name},
    ## @var{value})} takes one of @qcode{'KFold'}, @qcode{'Holdout'},
    ## @qcode{'Leaveout'} and @qcode{'CVPartition'} to say how to partition,
    ## and any option @code{RegressionLinear} takes to say how to fit.
    ## @qcode{'CrossVal'} is accepted and has no effect here, this class
    ## being cross-validated by construction.
    ##
    ## Anything left as @qcode{'auto'} is resolved by each fold against its
    ## own training rows rather than once over the whole data, so ten folds
    ## of a hundred observations each get a @qcode{Lambda} of one ninetieth
    ## rather than one hundredth, and each its own @qcode{Epsilon}.  Both are
    ## MATLAB's behaviour, measured.
    ##
    ## @seealso{fitrlinear, RegressionLinear}
    ## @end deftypefn
    function this = RegressionPartitionedLinear (X, Y, varargin)

      if (nargin < 2)
        error (strcat ("RegressionPartitionedLinear: too few input", ...
                       " arguments."));
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("RegressionPartitionedLinear: optional arguments", ...
                       " must be given in Name-Value pairs."));
      endif

      [P, args] = partitionedArgs (varargin, 'RegressionPartitionedLinear');
      F = regFrame (X, Y, P.Weights, 'RegressionPartitionedLinear');
      [part, args] = cvPartitionOf (args, [], F.n, ...
                                    'RegressionPartitionedLinear');

      ## Every fold is fitted without a response transform, so that the one
      ## the parent carries is applied exactly once to the assembled
      ## predictions.
      fargs = [args, {'ResponseTransform', 'none'}];
      G = struct ();
      G.X = F.X;
      G.Y = F.Y;
      G.Weights = F.Weights;
      this.Trained = foldModels ('RegressionLinear', G, part, fargs);

      this.NumObservations = F.n;
      this.Y = F.Y;
      this.W = F.W;
      this.X_ = F.X;
      this.KFold = part.NumTestSets;
      this.Partition = part;
      this.NumLambda_ = numel (this.Trained{1}.Lambda);
      this.HasEpsilon_ = ! isempty (this.Trained{1}.Epsilon);

      this.PredictorNames = this.Trained{1}.PredictorNames;
      this.CategoricalPredictors = this.Trained{1}.CategoricalPredictors;
      this.ResponseName = this.Trained{1}.ResponseName;

      if (isempty (P.ResponseTransform))
        P.ResponseTransform = 'none';
      endif
      [this.RTfun, this.ResponseTransform] = ...
              parseResponseTransform (P.ResponseTransform, ...
                                      'RegressionPartitionedLinear');

      this.ModelParameters = struct ('Type', 'regression', ...
                                     'Method', 'PartitionedLinear', ...
                                     'LearnerTemplates', 'Linear', ...
                                     'NLearn', this.KFold);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionPartitionedLinear} {@var{yFit} =} kfoldPredict (@var{obj})
    ##
    ## Out-of-fold prediction for every observation.
    ##
    ## Each observation is predicted by the fold that held it out, so the
    ## predictions are out-of-sample.  An observation that no fold held out,
    ## which under a holdout partition is most of them, comes back
    ## @qcode{NaN}.
    ##
    ## With @math{L} regularization strengths @var{yFit} has one column per
    ## strength.
    ##
    ## @end deftypefn
    function yFit = kfoldPredict (this)

      yFit = this.RTfun (kfoldResponse (this.Trained, this.Partition, ...
                                        this.X_, this.NumLambda_));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {RegressionPartitionedLinear} {@var{l} =} kfoldLoss (@var{obj})
    ## @deftypefnx {RegressionPartitionedLinear} {@var{l} =} kfoldLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-fold regression loss.
    ##
    ## @code{@var{l} = kfoldLoss (@var{obj})} returns the out-of-fold mean
    ## squared error.
    ##
    ## @code{@var{l} = kfoldLoss (@dots{}, @var{name}, @var{value})} takes
    ## @qcode{'LossFun'}, either @qcode{'mse'} or
    ## @qcode{'epsiloninsensitive'}; @qcode{'Folds'}, a subset of the folds
    ## to average over; and @qcode{'Mode'}, either @qcode{'average'}, the
    ## default, or @qcode{'individual'}, which returns one row per fold.
    ##
    ## @end deftypefn
    function l = kfoldLoss (this, varargin)

      O = kfoldOpts (varargin, {'mse', 'epsiloninsensitive'}, ...
                     'RegressionPartitionedLinear', 'kfoldLoss', this.KFold);
      if (isempty (O.LossFun))
        O.LossFun = 'mse';
      endif
      if (strcmp (O.LossFun, 'epsiloninsensitive') && ! this.HasEpsilon_)
        error (strcat ("RegressionPartitionedLinear.kfoldLoss: the", ...
                       " 'epsiloninsensitive' loss applies to a support", ...
                       " vector machine only."));
      endif

      yFit = kfoldPredict (this);

      ## Each fold has its own insensitive band, so a residual is judged
      ## against the band of the fold that produced it.  Assembled per
      ## observation, the band follows the prediction.
      band = nan (this.NumObservations, 1);
      for i = 1:this.KFold
        idx = test (this.Partition, i);
        if (any (idx) && ! isempty (this.Trained{i}.Epsilon))
          band(idx) = this.Trained{i}.Epsilon;
        endif
      endfor

      sets = foldSets (this.Partition, O.Folds, O.Mode, ...
                       this.NumObservations);
      l = zeros (numel (sets), this.NumLambda_);
      for i = 1:numel (sets)
        idx = sets{i};
        w = this.W(idx);
        w = w / sum (w);
        for k = 1:this.NumLambda_
          r = this.Y(idx) - yFit(idx,k);
          if (strcmp (O.LossFun, 'mse'))
            l(i,k) = sum (w .* (r .^ 2));
          else
            l(i,k) = sum (w .* max (0, abs (r) - band(idx)));
          endif
        endfor
      endfor

    endfunction

  endmethods

  methods (Access = public, Hidden)

    function this = set.ResponseTransform (this, val)
      [f, nm] = parseResponseTransform (val, 'RegressionPartitionedLinear');
      this.ResponseTransform = nm;
      this.RTfun = f;
    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        printf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      printf ("\n  RegressionPartitionedLinear\n\n");
      printf ("%+25s: '%s'\n", 'CrossValidatedModel', ...
              this.CrossValidatedModel);
      printf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      printf ("%+25s: '%s'\n", 'ResponseTransform', this.ResponseTransform);
      printf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      printf ("%+25s: %d\n", 'KFold', this.KFold);
      printf ("\n");
    endfunction

  endmethods

endclassdef

%!demo
%! ## Cross-validate a linear regression of fuel consumption and read the
%! ## out-of-sample mean squared error.
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! CVMdl = RegressionPartitionedLinear (X(ok,:), MPG(ok), 'KFold', 5)
%! outOfSample = kfoldLoss (CVMdl)

%!shared X, Y
%! load carsmall
%! X = [Acceleration, Displacement, Horsepower, Weight];
%! ok = ! any (isnan ([X, MPG]), 2);
%! X = X(ok,:);
%! Y = MPG(ok);

%!test
%! ## The model reports the surface MATLAB reports
%! CVMdl = RegressionPartitionedLinear (X, Y, 'KFold', 4);
%! assert_equal (class (CVMdl), 'RegressionPartitionedLinear');
%! assert_equal (CVMdl.CrossValidatedModel, 'Linear');
%! assert_equal (CVMdl.KFold, 4);
%! assert_equal (CVMdl.NumObservations, 93);
%! assert_equal (CVMdl.ResponseTransform, 'none');
%! assert_equal (CVMdl.ResponseName, 'Y');
%! assert_equal (CVMdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (size (CVMdl.Trained), [4, 1]);
%! assert_equal (class (CVMdl.Trained{1}), 'RegressionLinear');
%! assert_equal (sum (CVMdl.W), 1, 1e-12);
%! assert_equal (CVMdl.ModelParameters.Method, 'PartitionedLinear');

%!test
%! ## The properties are the ones MATLAB lists, in its order
%! CVMdl = RegressionPartitionedLinear (X, Y);
%! assert_equal (sort (properties (CVMdl)), ...
%!               sort ({'ResponseTransform'; 'CrossValidatedModel'; ...
%!                      'NumObservations'; 'Y'; 'W'; 'PredictorNames'; ...
%!                      'CategoricalPredictors'; 'ResponseName'; 'Trained'; ...
%!                      'KFold'; 'Partition'; 'ModelParameters'}));

%!test
%! ## The default is ten folds
%! CVMdl = RegressionPartitionedLinear (X, Y);
%! assert_equal (CVMdl.KFold, 10);
%! assert_equal (numel (CVMdl.Trained), 10);

%!test
%! ## Each fold resolves its own 'Lambda' and its own 'Epsilon' from its own
%! ## training rows, which is why two folds report different bands
%! CVMdl = RegressionPartitionedLinear (X, Y, 'KFold', 4);
%! n1 = sum (training (CVMdl.Partition, 1));
%! assert_equal (CVMdl.Trained{1}.Lambda, 1 / n1, 1e-15);
%! assert_equal (isequal (CVMdl.Trained{1}.Epsilon, ...
%!                        CVMdl.Trained{4}.Epsilon), false);

%!test
%! ## kfoldPredict predicts each observation with the fold that held it out,
%! ## and doing the same partition by hand gives the same numbers
%! part = cvpartition (93, 'KFold', 4);
%! CVMdl = RegressionPartitionedLinear (X, Y, 'CVPartition', part, ...
%!                                      'Learner', 'leastsquares');
%! yFit = kfoldPredict (CVMdl);
%! byhand = nan (93, 1);
%! for k = 1:4
%!   tr = training (part, k);
%!   te = test (part, k);
%!   m = RegressionLinear (X(tr,:), Y(tr), 'Learner', 'leastsquares');
%!   byhand(te) = predict (m, X(te,:));
%! endfor
%! assert_equal (yFit, byhand, 1e-10);

%!test
%! ## Averaging pools the observations rather than the per-fold values, and
%! ## the two differ once the folds differ in size.  The pooled reading is
%! ## the one R2024a returns.
%! part = cvpartition (93, 'KFold', 4);
%! CVMdl = RegressionPartitionedLinear (X, Y, 'CVPartition', part, ...
%!                                      'Learner', 'leastsquares');
%! yFit = kfoldPredict (CVMdl);
%! assert_equal (kfoldLoss (CVMdl), mean ((Y - yFit) .^ 2), 1e-10);
%! each = kfoldLoss (CVMdl, 'Mode', 'individual');
%! assert_equal (size (each), [4, 1]);
%! assert_equal (isequal (kfoldLoss (CVMdl), mean (each)), false);

%!test
%! ## 'Folds' reports over the observations of the folds it names
%! CVMdl = RegressionPartitionedLinear (X, Y, 'KFold', 4, ...
%!                                      'Learner', 'leastsquares');
%! yFit = kfoldPredict (CVMdl);
%! idx = test (CVMdl.Partition, 1) | test (CVMdl.Partition, 3);
%! assert_equal (kfoldLoss (CVMdl, 'Folds', [1, 3]), ...
%!               mean ((Y(idx) - yFit(idx)) .^ 2), 1e-10);

%!test
%! ## The epsilon-insensitive loss judges each residual against the band of
%! ## the fold that produced it, and is offered by a support vector machine
%! ## alone
%! CVMdl = RegressionPartitionedLinear (X, Y, 'KFold', 4);
%! assert_equal (isfinite (kfoldLoss (CVMdl, 'LossFun', ...
%!                                    'epsiloninsensitive')), true);
%! assert_equal (kfoldLoss (CVMdl, 'LossFun', 'epsiloninsensitive') ...
%!               < kfoldLoss (CVMdl), true);

%!test
%! ## A whole regularization path gives one column per strength
%! CVMdl = RegressionPartitionedLinear (X, Y, 'KFold', 4, ...
%!                                      'Lambda', [0.001, 0.01, 0.1]);
%! assert_equal (size (kfoldPredict (CVMdl)), [93, 3]);
%! assert_equal (size (kfoldLoss (CVMdl)), [1, 3]);
%! assert_equal (size (kfoldLoss (CVMdl, 'Mode', 'individual')), [4, 3]);

%!test
%! ## An observation that no fold held out comes back NaN rather than
%! ## predicted
%! CVMdl = RegressionPartitionedLinear (X, Y, 'Holdout', 0.3);
%! assert_equal (CVMdl.KFold, 1);
%! yFit = kfoldPredict (CVMdl);
%! assert_equal (sum (isnan (yFit)), 65);
%! assert_equal (isfinite (kfoldLoss (CVMdl)), true);

%!test
%! ## A response transform reaches the assembled predictions once, the fold
%! ## models carrying none
%! part = cvpartition (93, 'KFold', 4);
%! plain = RegressionPartitionedLinear (X, Y, 'CVPartition', part, ...
%!                                      'Learner', 'leastsquares');
%! CVexp = RegressionPartitionedLinear (X, Y, 'CVPartition', part, ...
%!                                      'Learner', 'leastsquares', ...
%!                                      'ResponseTransform', 'exp');
%! assert_equal (CVexp.ResponseTransform, 'exp');
%! assert_equal (CVexp.Trained{1}.ResponseTransform, 'none');
%! assert_equal (kfoldPredict (CVexp), exp (kfoldPredict (plain)), 1e-10);

%!test
%! ## A row with a missing value is dropped before the partition
%! Xn = X;
%! Xn(3,2) = NaN;
%! CVMdl = RegressionPartitionedLinear (Xn, Y, 'KFold', 4);
%! assert_equal (CVMdl.NumObservations, 92);
%! assert_equal (CVMdl.Partition.NumObservations, 92);

%!test
%! ## It can be assigned after the model is built, and reaches kfoldPredict
%! ## without being carried into the folds
%! CVMdl = RegressionPartitionedLinear (X, Y, 'KFold', 4);
%! y0 = kfoldPredict (CVMdl);
%! CVMdl.ResponseTransform = @(y) y + 100;
%! y1 = kfoldPredict (CVMdl);
%! assert_equal (CVMdl.Trained{1}.ResponseTransform, 'none');
%! assert_equal (y1, y0 + 100, 1e-10);

%!test
%! ## And it reaches kfoldLoss, which is computed from those predictions
%! CVMdl = RegressionPartitionedLinear (X, Y, 'KFold', 4);
%! before = kfoldLoss (CVMdl);
%! CVMdl.ResponseTransform = @(y) y + 100;
%! assert (kfoldLoss (CVMdl) > before);

%!test
%! ## 'none' is the identity, so assigning it transforms nothing
%! CVMdl = RegressionPartitionedLinear (X, Y, 'KFold', 4);
%! y0 = kfoldPredict (CVMdl);
%! CVMdl.ResponseTransform = 'none';
%! assert_equal (kfoldPredict (CVMdl), y0);

%!error<RegressionPartitionedLinear: unrecognized 'ResponseTransform' function.> ...
%! CVMdl = RegressionPartitionedLinear (X, Y, 'KFold', 4);
%! CVMdl.ResponseTransform = 'nosuchtransform';

## Test input validation
%!error<RegressionPartitionedLinear: too few input arguments.> ...
%! RegressionPartitionedLinear (ones (10, 2))
%!error<RegressionPartitionedLinear: optional arguments must be given in Name-Value pairs.> ...
%! RegressionPartitionedLinear (ones (10, 2), ones (10, 1), 'KFold')
%!error<RegressionPartitionedLinear: 'KFold' must be an integer value greater than 1.> ...
%! RegressionPartitionedLinear (ones (10, 2), ones (10, 1), 'KFold', 1)
%!error<RegressionPartitionedLinear: you can use only one of 'KFold', 'Holdout', 'Leaveout', or 'CVPartition' options.> ...
%! RegressionPartitionedLinear (ones (10, 2), ones (10, 1), 'KFold', 2, ...
%!                     'Holdout', 0.2)
%!error<RegressionPartitionedLinear: invalid values in Y.> ...
%! RegressionPartitionedLinear (ones (10, 2), {1, 2})
%!error<RegressionPartitionedLinear.kfoldLoss: 'LossFun' must be 'mse' or 'epsiloninsensitive'.> ...
%! kfoldLoss (RegressionPartitionedLinear (ones (10, 2), ones (10, 1), ...
%!                     'KFold', 2), 'LossFun', 'hinge')
%!error<RegressionPartitionedLinear.kfoldLoss: the 'epsiloninsensitive' loss applies to a support vector machine only.> ...
%! kfoldLoss (RegressionPartitionedLinear (ones (10, 2), ones (10, 1), ...
%!                     'KFold', 2, 'Learner', 'leastsquares'), 'LossFun', ...
%!                     'epsiloninsensitive')
%!error<RegressionPartitionedLinear.kfoldLoss: 'Mode' must be either 'average' or 'individual'.> ...
%! kfoldLoss (RegressionPartitionedLinear (ones (10, 2), ones (10, 1), ...
%!                     'KFold', 2), 'Mode', 'each')
