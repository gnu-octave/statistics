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

classdef ClassificationPartitionedECOC
## -*- texinfo -*-
## @deftypefn {statistics} ClassificationPartitionedECOC
##
## A cross-validated multiclass model built from binary learners.
##
## Each fold holds out part of the data, fits an error correcting output codes
## model on the rest, and answers the part it held out, so every observation
## is classified by a model that never saw it.
##
## It comes from @code{crossval} on a @code{ClassificationECOC}, and from
## @code{fitcecoc} given any of @qcode{'KFold'}, @qcode{'Holdout'},
## @qcode{'Leaveout'} or @qcode{'CVPartition'}.
##
## This is the one cross-validated class in the package that is not the
## general @code{ClassificationPartitionedModel}.  It carries
## @code{CodingMatrix}, @code{BinaryLoss} and @code{BinaryY}, three things
## the general class has nowhere to put and without which a fold's scores
## cannot be decoded at all.
##
## @seealso{fitcecoc, ClassificationECOC, CompactClassificationECOC}
## @end deftypefn

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} BinaryY
    ##
    ## What each observation was to each binary learner, an @math{NxL} matrix
    ## of -1, 0 and +1.  This property is read-only.
    ##
    ## @end deftp
    BinaryY               = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} CodingMatrix
    ##
    ## The coding design every fold was fitted on, a @math{KxL} matrix.  This
    ## property is read-only.
    ##
    ## @end deftp
    CodingMatrix          = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} BinaryLoss
    ##
    ## The loss the binary scores of a fold are read with.  This property is
    ## read-only.
    ##
    ## @end deftp
    BinaryLoss            = 'hinge';

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} ClassNames
    ##
    ## The distinct class labels.  This property is read-only.
    ##
    ## @end deftp
    ClassNames            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} Cost
    ##
    ## The cost of misclassification.  This property is read-only.
    ##
    ## @end deftp
    Cost                  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} Prior
    ##
    ## The prior probability of each class.  This property is read-only.
    ##
    ## @end deftp
    Prior                 = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} CrossValidatedModel
    ##
    ## The kind of model that was cross-validated, @qcode{'ECOC'} here.  This
    ## property is read-only.
    ##
    ## @end deftp
    CrossValidatedModel   = 'ECOC';

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} PredictorNames
    ##
    ## The names of the predictors.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} CategoricalPredictors
    ##
    ## The columns holding categorical predictors, always empty here.  This
    ## property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} ResponseName
    ##
    ## The name of the response variable.  This property is read-only.
    ##
    ## @end deftp
    ResponseName          = 'Y';

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} NumObservations
    ##
    ## The number of observations used.  This property is read-only.
    ##
    ## @end deftp
    NumObservations       = 0;

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} X
    ##
    ## The predictor data the folds were fitted on.  This property is
    ## read-only.
    ##
    ## @end deftp
    X                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} Y
    ##
    ## The class labels the folds were fitted on.  This property is
    ## read-only.
    ##
    ## @end deftp
    Y                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} W
    ##
    ## The observation weights.  This property is read-only.
    ##
    ## @end deftp
    W                     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} ModelParameters
    ##
    ## A structure of the options the folds were fitted with.  This property
    ## is read-only.
    ##
    ## @end deftp
    ModelParameters       = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} Trained
    ##
    ## The models the folds were fitted to
    ##
    ## A cell column with one @code{CompactClassificationECOC} per fold, the
    ## training data having no place in a model that only ever answers the
    ## rows it did not see.  This property is read-only.
    ##
    ## @end deftp
    Trained               = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} KFold
    ##
    ## The number of folds.  This property is read-only.
    ##
    ## @end deftp
    KFold                 = 0;

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} Partition
    ##
    ## The @code{cvpartition} object that says which rows each fold held out.
    ## This property is read-only.
    ##
    ## @end deftp
    Partition             = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} BinEdges
    ##
    ## The bin edges of the predictors, always empty here.  This property is
    ## read-only.
    ##
    ## @end deftp
    BinEdges              = {};

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationPartitionedECOC} {property} ScoreTransform
    ##
    ## The transform applied to the assembled out-of-fold scores.  The folds
    ## never carry it: it is applied once to what they return, which is what
    ## R2024a does.
    ##
    ## @end deftp
    ScoreTransform        = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    STfun = @(x) x;
  endproperties

  methods (Hidden)

    function this = set.ScoreTransform (this, val)
      name = 'ClassificationPartitionedECOC';
      try
        [this.STfun, this.ScoreTransform] = parseScoreTransform (val, name);
      catch
        error (strcat ("ClassificationPartitionedECOC.subsasgn:", ...
                       " 'ScoreTransform' must be a character vector or a", ...
                       " 'function_handle' object."));
      end_try_catch
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationPartitionedECOC} {@var{obj} =} ClassificationPartitionedECOC (@var{Mdl}, @var{Partition})
    ##
    ## Create a @code{ClassificationPartitionedECOC} object.
    ##
    ## @var{Mdl} is the @code{ClassificationECOC} to cross validate and
    ## @var{Partition} the @code{cvpartition} saying which rows each fold
    ## holds out.  The documented way to reach this constructor is
    ## @code{crossval}.
    ##
    ## @end deftypefn
    function this = ClassificationPartitionedECOC (Mdl, Partition)

      if (nargin < 2)
        error (strcat ("ClassificationPartitionedECOC: too few input", ...
                       " arguments."));
      endif
      if (! strcmp (class (Mdl), 'ClassificationECOC'))
        error (strcat ("ClassificationPartitionedECOC: MDL must be a", ...
                       " 'ClassificationECOC' object."));
      endif
      if (! isa (Partition, 'cvpartition'))
        error (strcat ("ClassificationPartitionedECOC: PARTITION must be a", ...
                       " 'cvpartition' object."));
      endif

      this.ClassNames            = Mdl.ClassNames;
      this.Cost                  = Mdl.Cost;
      this.Prior                 = Mdl.Prior;
      this.PredictorNames        = Mdl.PredictorNames;
      this.CategoricalPredictors = Mdl.CategoricalPredictors;
      this.ResponseName          = Mdl.ResponseName;
      this.NumObservations       = Mdl.NumObservations;
      this.X                     = Mdl.X;
      this.Y                     = Mdl.Y;
      this.W                     = Mdl.W;
      this.BinaryY               = Mdl.BinaryY;
      this.CodingMatrix          = Mdl.CodingMatrix;
      this.BinaryLoss            = Mdl.BinaryLoss;
      this.Partition             = Partition;
      this.KFold                 = Partition.NumTestSets;
      this.ModelParameters       = Mdl.ModelParameters;

      ## Each fold is fitted on the coding matrix the parent settled, not on
      ## one of its own: a random design would differ fold to fold and the
      ## out-of-fold scores could not be put side by side.
      fargs = {'Coding', Mdl.CodingMatrix, ...
               'Learners', Mdl.ModelParameters.BinaryLearners, ...
               'ClassNames', Mdl.ClassNames, 'Prior', Mdl.Prior, ...
               'Cost', Mdl.Cost, 'BinaryLoss', Mdl.BinaryLoss, ...
               'ResponseName', Mdl.ResponseName, ...
               'PredictorNames', Mdl.PredictorNames};

      G = struct ();
      G.X = Mdl.X;
      G.Y = Mdl.Y;
      G.Weights = Mdl.W;
      Trained = foldModels ('ClassificationECOC', G, Partition, fargs);

      ## A fold answers rows it never saw, so it keeps no training data.
      ## R2024a stores a CompactClassificationECOC here too.
      for k = 1:numel (Trained)
        Trained{k} = compact (Trained{k});
      endfor
      this.Trained = Trained;

    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ("\n  ClassificationPartitionedECOC\n\n");
      fprintf ("%+25s: '%s'\n", 'CrossValidatedModel', ...
               this.CrossValidatedModel);
      fprintf ("%+25s: %s\n", 'PredictorNames', ...
               classNameListing (this.PredictorNames));
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'KFold', this.KFold);
      fprintf ("%+25s: %s\n", 'ClassNames', classNameListing (this.ClassNames));
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("\n");
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedECOC} {@var{label} =} kfoldPredict (@var{obj})
    ## @deftypefnx {ClassificationPartitionedECOC} {[@var{label}, @var{NegLoss}] =} kfoldPredict (@var{obj})
    ##
    ## Out-of-fold class of every observation.
    ##
    ## Each observation is classified by the fold that held it out, so the
    ## labels are out-of-sample.  An observation no fold held out, which
    ## under a holdout partition is most of them, comes back missing rather
    ## than classified, and its @var{NegLoss} row comes back @qcode{NaN}.
    ##
    ## @seealso{ClassificationPartitionedECOC,
    ## CompactClassificationECOC.predict}
    ## @end deftypefn
    function [label, NegLoss] = kfoldPredict (this)

      n = this.NumObservations;
      K = numel (this.ClassNames);
      if (iscellstr (this.Y))
        label = repmat ({''}, n, 1);
      elseif (islogical (this.Y))
        label = false (n, 1);
      elseif (ischar (this.Y))
        label = repmat (' ', n, columns (this.Y));
      else
        label = nan (n, 1);
      endif
      NegLoss = nan (n, K);

      for k = 1:this.KFold
        idx = test (this.Partition, k);
        if (! any (idx))
          continue;
        endif
        [lab, NL] = predict (this.Trained{k}, this.X(idx,:));
        label(idx,:) = lab;
        NegLoss(idx,:) = NL;
      endfor

      ## The folds never carry the transform, so it is applied once here.
      NegLoss = this.STfun (NegLoss);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationPartitionedECOC} {@var{m} =} kfoldMargin (@var{obj})
    ##
    ## Out-of-fold classification margin of every observation.
    ##
    ## @seealso{ClassificationPartitionedECOC.kfoldPredict,
    ## ClassificationPartitionedECOC.kfoldEdge}
    ## @end deftypefn
    function m = kfoldMargin (this)

      [~, NegLoss] = kfoldPredict (this);
      [gY, errmsg] = labelIndices (this.ClassNames, this.Y);
      if (! isempty (errmsg))
        error ("ClassificationPartitionedECOC.kfoldMargin: %s", errmsg);
      endif
      m = marginsOf (NegLoss, gY, 1);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationPartitionedECOC} {@var{e} =} kfoldEdge (@var{obj})
    ##
    ## Out-of-fold classification edge, the weighted mean of the margins.
    ##
    ## @seealso{ClassificationPartitionedECOC.kfoldMargin}
    ## @end deftypefn
    function e = kfoldEdge (this)

      m = kfoldMargin (this);
      keep = ! isnan (m);
      w = this.W(:)(keep);
      e = sum (w .* m(keep)) / sum (w);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationPartitionedECOC} {@var{L} =} kfoldLoss (@var{obj})
    ## @deftypefnx {ClassificationPartitionedECOC} {@var{L} =} kfoldLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Out-of-fold classification loss.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab @tab @qcode{'classiferror'} (default),
    ## @qcode{'classifcost'}, @qcode{'mincost'}, @qcode{'binodeviance'},
    ## @qcode{'exponential'}, @qcode{'hinge'}, @qcode{'logit'} or
    ## @qcode{'quadratic'}.
    ## @end multitable
    ##
    ## @seealso{ClassificationPartitionedECOC.kfoldPredict}
    ## @end deftypefn
    function L = kfoldLoss (this, varargin)

      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationPartitionedECOC.kfoldLoss:", ...
                       " name-value arguments must be in pairs."));
      endif
      LossFun = 'classiferror';
      for i = 1:2:numel (varargin)
        switch (tolower (varargin{i}))
          case 'lossfun'
            LossFun = varargin{i+1};
            if (! (ischar (LossFun) && isrow (LossFun)))
              error (strcat ("ClassificationPartitionedECOC.kfoldLoss:", ...
                             " 'LossFun' must be a character vector."));
            endif
            LossFun = tolower (LossFun);
          otherwise
            error (strcat ("ClassificationPartitionedECOC.kfoldLoss:", ...
                           " invalid parameter name in optional pair", ...
                           " arguments."));
        endswitch
      endfor

      [~, NegLoss] = kfoldPredict (this);
      [gY, errmsg] = labelIndices (this.ClassNames, this.Y);
      if (! isempty (errmsg))
        error ("ClassificationPartitionedECOC.kfoldLoss: %s", errmsg);
      endif

      ## An observation no fold held out has no out-of-sample answer and is
      ## left out of the loss rather than counted as an error.
      keep = ! any (isnan (NegLoss), 2);
      L = classificationLoss (LossFun, NegLoss(keep,:), gY(keep), ...
                              this.W(keep), this.Cost);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationPartitionedECOC} {@var{vals} =} kfoldfun (@var{obj}, @var{fun})
    ##
    ## Apply a function to every fold.
    ##
    ## @var{fun} is called once per fold as
    ## @code{@var{fun} (@var{CMdl}, @var{Xtrain}, @var{Ytrain}, @var{Wtrain},
    ## @var{Xtest}, @var{Ytest}, @var{Wtest})} and must return a numeric row
    ## of the same length every time.  @var{vals} has one row per fold.
    ##
    ## @seealso{ClassificationPartitionedECOC.kfoldPredict}
    ## @end deftypefn
    function vals = kfoldfun (this, fun)

      if (nargin < 2)
        error (strcat ("ClassificationPartitionedECOC.kfoldfun: too few", ...
                       " input arguments."));
      endif
      if (! isa (fun, 'function_handle'))
        error (strcat ("ClassificationPartitionedECOC.kfoldfun: FUN must", ...
                       " be a function handle."));
      endif

      vals = [];
      for k = 1:this.KFold
        tr = training (this.Partition, k);
        te = test (this.Partition, k);
        v = fun (this.Trained{k}, this.X(tr,:), this.Y(tr,:), this.W(tr), ...
                 this.X(te,:), this.Y(te,:), this.W(te));
        if (! (isnumeric (v) && isrow (v)))
          error (strcat ("ClassificationPartitionedECOC.kfoldfun: FUN must", ...
                         " return a numeric row vector."));
        endif
        if (! isempty (vals) && numel (v) != columns (vals))
          error (strcat ("ClassificationPartitionedECOC.kfoldfun: FUN must", ...
                         " return the same number of values for every", ...
                         " fold."));
        endif
        vals = [vals; v];
      endfor

    endfunction

  endmethods

endclassdef

## Tests
%!test  # MATLAB parity: the property surface of a cross-validated model
%! load fisheriris
%! c = cvpartition ('CustomPartition', repmat ((1:5)', 30, 1));
%! CV = fitcecoc (meas, species, 'Learners', 'tree', 'CVPartition', c);
%! assert_equal (class (CV), 'ClassificationPartitionedECOC');
%! assert_equal (numel (properties (CV)), 20);
%! assert_equal (CV.CrossValidatedModel, 'ECOC');
%! assert_equal (CV.KFold, 5);
%! assert_equal (numel (CV.Trained), 5);
%! assert_equal (class (CV.Trained{1}), 'CompactClassificationECOC');

%!test  # MATLAB parity: the out-of-fold answers of a discriminant code
%! ## Measured on R2024a over a partition both engines are given outright.
%! ## The learner is a discriminant and not a tree on purpose: a tree fitted
%! ## on four fifths of this fixture has ties among its best splits, which
%! ## the two engines break differently and which no test can pin.
%! load fisheriris
%! c = cvpartition ('CustomPartition', repmat ((1:5)', 30, 1));
%! CV = fitcecoc (meas, species, 'Learners', 'discriminant', 'CVPartition', c);
%! assert_equal (kfoldLoss (CV), 0.0333333333333334, 1e-12);
%! assert_equal (kfoldEdge (CV), 0.923697865527129, 1e-12);
%! m = kfoldMargin (CV);
%! assert_equal (m([1, 20, 51, 70, 101, 130])', ...
%!               [1, 1, 0.999820708981217, 0.999920781404916, ...
%!                0.999999631930005, -0.0807626827878425], 1e-12);

%!test  # MATLAB parity: the labels a tree code gives out of fold
%! load fisheriris
%! c = cvpartition ('CustomPartition', repmat ((1:5)', 30, 1));
%! CV = fitcecoc (meas, species, 'Learners', 'tree', 'CVPartition', c);
%! [label, NegLoss] = kfoldPredict (CV);
%! assert_equal (size (NegLoss), [150, 3]);
%! assert_equal (sum (! strcmp (label, species)), 11);
%! assert_equal (kfoldLoss (CV), 11/150, 1e-12);

%!test  # every observation is answered by the fold that held it out
%! load fisheriris
%! c = cvpartition ('CustomPartition', repmat ((1:5)', 30, 1));
%! CV = fitcecoc (meas, species, 'Learners', 'discriminant', 'CVPartition', c);
%! [label, NegLoss] = kfoldPredict (CV);
%! assert_equal (any (cellfun (@isempty, label)), false);
%! assert_equal (any (any (isnan (NegLoss))), false);

%!test  # a holdout partition leaves the rows no fold tested unanswered
%! load fisheriris
%! CV = crossval (ClassificationECOC (meas, species, 'Learners', ...
%!                                    'discriminant'), 'Holdout', 0.2);
%! [label, NegLoss] = kfoldPredict (CV);
%! assert_equal (sum (cellfun (@isempty, label)), 120);
%! assert_equal (sum (all (isnan (NegLoss), 2)), 120);

%!test  # the margin is positive exactly where the out-of-fold label was right
%! load fisheriris
%! c = cvpartition ('CustomPartition', repmat ((1:5)', 30, 1));
%! CV = fitcecoc (meas, species, 'Learners', 'discriminant', 'CVPartition', c);
%! assert_equal (kfoldMargin (CV) > 0, strcmp (kfoldPredict (CV), species));

%!test  # kfoldfun is called once per fold and its rows are stacked
%! load fisheriris
%! c = cvpartition ('CustomPartition', repmat ((1:5)', 30, 1));
%! CV = fitcecoc (meas, species, 'Learners', 'discriminant', 'CVPartition', c);
%! v = kfoldfun (CV, @(m, xt, yt, wt, xs, ys, ws) [rows(xt), rows(xs)]);
%! assert_equal (v, [120, 30; 120, 30; 120, 30; 120, 30; 120, 30]);

%!test  # crossval and the fit route give the same partitioned class
%! load fisheriris
%! c = cvpartition ('CustomPartition', repmat ((1:5)', 30, 1));
%! a = fitcecoc (meas, species, 'Learners', 'discriminant', 'CVPartition', c);
%! b = crossval (ClassificationECOC (meas, species, 'Learners', ...
%!                                   'discriminant'), 'CVPartition', c);
%! assert_equal (kfoldLoss (a), kfoldLoss (b), 1e-12);

## Test input validation
%!error<ClassificationPartitionedECOC: too few input arguments.> ...
%! ClassificationPartitionedECOC (1)
%!error<ClassificationPartitionedECOC: MDL must be a 'ClassificationECOC' object.> ...
%! ClassificationPartitionedECOC (1, cvpartition (10, 'KFold', 2))
%!error<ClassificationPartitionedECOC: PARTITION must be a 'cvpartition' object.> ...
%! Mdl = ClassificationECOC (ones (4, 2), [1; 2; 1; 2]); ...
%! ClassificationPartitionedECOC (Mdl, 1)
%!error<ClassificationPartitionedECOC.kfoldLoss: name-value arguments must be in pairs.> ...
%! y = [1; 2; 1; 2; 1; 2; 1; 2]; ...
%! CV = crossval (ClassificationECOC (ones (8, 2), y), 'KFold', 2); ...
%! kfoldLoss (CV, 'LossFun')
%!error<ClassificationPartitionedECOC.kfoldfun: FUN must be a function handle.> ...
%! y = [1; 2; 1; 2; 1; 2; 1; 2]; ...
%! CV = crossval (ClassificationECOC (ones (8, 2), y), 'KFold', 2); ...
%! kfoldfun (CV, 1)
