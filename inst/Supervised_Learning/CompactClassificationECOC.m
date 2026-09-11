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

classdef CompactClassificationECOC
## -*- texinfo -*-
## @deftypefn {statistics} CompactClassificationECOC
##
## A multiclass model built from binary learners, without its training data.
##
## A @code{CompactClassificationECOC} carries the binary learners of an error
## correcting output codes model and the coding matrix that says what each of
## them was trained to tell apart, and nothing else: the predictor data, the
## labels and the weights are gone, so it predicts and scores new data but
## cannot be refitted or cross validated.
##
## It comes from @code{compact} on a @code{ClassificationECOC}, and from
## @code{fitcecoc} itself when the binary learners are linear or kernel
## classifiers, which carry no training data of their own.
##
## @seealso{fitcecoc, ClassificationECOC, designecoc}
## @end deftypefn

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} BinaryLearners
    ##
    ## The trained binary learners
    ##
    ## A cell column with one trained model per column of
    ## @code{CodingMatrix}, each telling the classes that column marks +1
    ## from those it marks -1.  This property is read-only.
    ##
    ## @end deftp
    BinaryLearners        = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} CodingMatrix
    ##
    ## The coding design
    ##
    ## A @math{KxL} matrix of -1, 0 and +1 with one row per class and one
    ## column per binary learner.  A class marked 0 took no part in that
    ## learner.  This property is read-only.
    ##
    ## @end deftp
    CodingMatrix          = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} LearnerWeights
    ##
    ## The weight each binary learner was trained on
    ##
    ## A row with one element per binary learner, the total observation
    ## weight of the classes that learner took part in.  This property is
    ## read-only.
    ##
    ## @end deftp
    LearnerWeights        = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} ClassNames
    ##
    ## The class labels
    ##
    ## The distinct labels seen in the training data, in the order the rows
    ## of @code{CodingMatrix}, @code{Prior} and @code{Cost} take them.  This
    ## property is read-only.
    ##
    ## @end deftp
    ClassNames            = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} PredictorNames
    ##
    ## The names of the predictors, one per column of the training data.
    ## This property is read-only.
    ##
    ## @end deftp
    PredictorNames        = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} ExpandedPredictorNames
    ##
    ## The names of the predictors as the learners saw them.  It differs from
    ## @code{PredictorNames} only where a categorical predictor was expanded,
    ## which this package does not do.  This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} CategoricalPredictors
    ##
    ## The columns of the training data that held categorical predictors,
    ## always empty here.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} ResponseName
    ##
    ## The name of the response variable.  This property is read-only.
    ##
    ## @end deftp
    ResponseName          = 'Y';

  endproperties

  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} BinaryLoss
    ##
    ## The loss that turns a binary learner's score into a cost
    ##
    ## One of @qcode{'binodeviance'}, @qcode{'exponential'},
    ## @qcode{'hamming'}, @qcode{'hinge'}, @qcode{'linear'},
    ## @qcode{'logit'} or @qcode{'quadratic'}.
    ##
    ## A learner scoring on @math{(-Inf,+Inf)} takes every one of them but
    ## @qcode{'quadratic'}, and one scoring on @math{[0,1]} takes only
    ## @qcode{'hamming'} and @qcode{'quadratic'}: a loss reads a score
    ## against the interval it was written for, and the other way round it
    ## would read a posterior as a signed score.  Assigning one the learners
    ## cannot take raises.
    ##
    ## @end deftp
    BinaryLoss            = 'hinge';

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} Prior
    ##
    ## The prior probability of each class, in the order of
    ## @code{ClassNames} and summing to one.
    ##
    ## @end deftp
    Prior                 = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} Cost
    ##
    ## The cost of misclassification
    ##
    ## A @math{KxK} matrix whose @math{(i,j)} element is the cost of calling
    ## a member of class @math{i} a member of class @math{j}.
    ##
    ## @end deftp
    Cost                  = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationECOC} {property} ScoreTransform
    ##
    ## The transform applied to the predicted scores
    ##
    ## A character vector naming one of the transforms
    ## @code{parseScoreTransform} accepts, or a function handle.
    ##
    ## @end deftp
    ScoreTransform        = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)
    ## The interval the binary learners score on, which decides what losses
    ## they can be read with.  MATLAB carries it inside each learner; here it
    ## is one property, every learner of a model being of one kind.
    ScoreRange            = [-Inf, Inf];
    STfun                 = @(x) x;
  endproperties

  ## Set methods for the properties a user may assign.
  methods (Hidden)

    function this = set.ScoreTransform (this, val)
      name = 'CompactClassificationECOC';
      try
        [this.STfun, this.ScoreTransform] = parseScoreTransform (val, name);
      catch
        error (strcat ("CompactClassificationECOC.subsasgn:", ...
                       " 'ScoreTransform' must be a character vector or a", ...
                       " 'function_handle' object."));
      end_try_catch
    endfunction

    function this = set.BinaryLoss (this, val)
      if (! (ischar (val) && isrow (val)))
        error (strcat ("CompactClassificationECOC.subsasgn: 'BinaryLoss'", ...
                       " must be a character vector."));
      endif
      [~, errmsg] = ecocDecode (zeros (1, columns (this.CodingMatrix)), ...
                                this.CodingMatrix, tolower (val), ...
                                'lossweighted', this.ScoreRange);
      if (! isempty (errmsg))
        error ("CompactClassificationECOC.subsasgn: %s", errmsg);
      endif
      this.BinaryLoss = tolower (val);
    endfunction

  endmethods

  methods (Hidden)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationECOC} {@var{obj} =} CompactClassificationECOC (@var{Mdl})
    ## @deftypefnx {CompactClassificationECOC} {@var{obj} =} CompactClassificationECOC ()
    ##
    ## Create a @code{CompactClassificationECOC} object.
    ##
    ## @var{Mdl} is the @code{ClassificationECOC} object to compact.  The
    ## documented way to reach this constructor is the @code{compact} method.
    ##
    ## Called with no arguments it returns an object with its properties at
    ## their defaults, which is what @code{loadmodel} fills in.
    ##
    ## @end deftypefn
    function this = CompactClassificationECOC (Mdl = [])

      if (isempty (Mdl))
        return;
      endif
      if (! strcmp (class (Mdl), 'ClassificationECOC'))
        error (strcat ("CompactClassificationECOC: MDL must be a", ...
                       " 'ClassificationECOC' object."));
      endif

      this.BinaryLearners         = Mdl.BinaryLearners;
      this.CodingMatrix           = Mdl.CodingMatrix;
      this.LearnerWeights         = Mdl.LearnerWeights;
      this.ClassNames             = Mdl.ClassNames;
      this.PredictorNames         = Mdl.PredictorNames;
      this.ExpandedPredictorNames = Mdl.ExpandedPredictorNames;
      this.CategoricalPredictors  = Mdl.CategoricalPredictors;
      this.ResponseName           = Mdl.ResponseName;
      this.Prior                  = Mdl.Prior;
      this.Cost                   = Mdl.Cost;
      this.ScoreRange             = Mdl.ScoreRange;
      this.BinaryLoss             = Mdl.BinaryLoss;
      this.ScoreTransform         = Mdl.ScoreTransform;

    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    ## Custom display
    function disp (this)
      fprintf ("\n  CompactClassificationECOC\n\n");
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      fprintf ("%+25s: %s\n", 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ("%+25s: %s\n", 'ClassNames', classNameListing (this.ClassNames));
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("%+25s: {%dx%d cell}\n", 'BinaryLearners', ...
               numel (this.BinaryLearners), 1);
      fprintf ("%+25s: [%dx%d double]\n", 'CodingMatrix', ...
               rows (this.CodingMatrix), columns (this.CodingMatrix));
      fprintf ("\n");
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationECOC} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationECOC} {[@var{label}, @var{NegLoss}] =} predict (@dots{})
    ## @deftypefnx {CompactClassificationECOC} {[@var{label}, @var{NegLoss}, @var{PBScore}] =} predict (@dots{})
    ## @deftypefnx {CompactClassificationECOC} {[@dots{}] =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Classify new data with a trained @code{CompactClassificationECOC}.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} sends each row of
    ## @var{XC} to every binary learner, turns the scores they return into a
    ## cost per class, and returns the class of least cost.  @var{XC} must
    ## have as many columns as the data the model was fitted on.
    ##
    ## @code{[@var{label}, @var{NegLoss}] = predict (@dots{})} also returns
    ## the @math{NxK} negated average loss, the largest entry of a row naming
    ## the class that row was given.
    ##
    ## @code{[@var{label}, @var{NegLoss}, @var{PBScore}] = predict (@dots{})}
    ## also returns the @math{NxL} scores the binary learners gave the class
    ## each was trained to call +1.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'BinaryLoss'} @tab @tab The loss to read the binary
    ## scores with, overriding the @code{BinaryLoss} property for this call.
    ##
    ## @item @qcode{'Decoding'} @tab @tab @qcode{'lossweighted'} (default) or
    ## @qcode{'lossbased'}.  The first averages the loss over the learners a
    ## class took part in, the second over every learner, a class that sat a
    ## column out costing the same there as a class at the decision boundary.
    ## @end multitable
    ##
    ## @seealso{CompactClassificationECOC, ClassificationECOC, fitcecoc}
    ## @end deftypefn
    function [label, NegLoss, PBScore] = predict (this, XC, varargin)

      if (nargin < 2)
        error ("CompactClassificationECOC.predict: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactClassificationECOC.predict: name-value", ...
                       " arguments must be in pairs."));
      endif
      if (isempty (XC))
        error ("CompactClassificationECOC.predict: XC is empty.");
      endif
      if (columns (XC) != numel (this.PredictorNames))
        error (strcat ("CompactClassificationECOC.predict: XC must have", ...
                       " the same number of predictors as the trained", ...
                       " model."));
      endif

      lossname = this.BinaryLoss;
      decoding = 'lossweighted';
      for i = 1:2:numel (varargin)
        switch (tolower (varargin{i}))
          case 'binaryloss'
            lossname = varargin{i+1};
            if (! (ischar (lossname) && isrow (lossname)))
              error (strcat ("CompactClassificationECOC.predict:", ...
                             " 'BinaryLoss' must be a character vector."));
            endif
            lossname = tolower (lossname);
          case 'decoding'
            decoding = varargin{i+1};
            if (! (ischar (decoding) && isrow (decoding)
                   && any (strcmpi (decoding, {'lossweighted', 'lossbased'}))))
              error (strcat ("CompactClassificationECOC.predict:", ...
                             " 'Decoding' must be 'lossweighted' or", ...
                             " 'lossbased'."));
            endif
            decoding = tolower (decoding);
          otherwise
            error (strcat ("CompactClassificationECOC.predict: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
      endfor

      PBScore = binaryScores (this, XC);
      [NegLoss, errmsg] = ecocDecode (PBScore, this.CodingMatrix, lossname, ...
                                      decoding, this.ScoreRange);
      if (! isempty (errmsg))
        error ("CompactClassificationECOC.predict: %s", errmsg);
      endif

      NegLoss = this.STfun (NegLoss);
      [~, idx] = max (NegLoss, [], 2);
      label = labelsFromIndex (this.ClassNames, idx);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationECOC} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationECOC} {@var{m} =} margin (@dots{}, @var{name}, @var{value})
    ##
    ## Classification margin of a @code{CompactClassificationECOC}.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns one
    ## margin per row of @var{X}: the negated loss of the true class less the
    ## largest negated loss among the others.  A positive margin means the
    ## row was classified correctly, and the larger it is the further the
    ## decision was from going the other way.
    ##
    ## It takes the same @qcode{'BinaryLoss'} and @qcode{'Decoding'}
    ## arguments @code{predict} does.
    ##
    ## @seealso{CompactClassificationECOC.predict,
    ## CompactClassificationECOC.edge}
    ## @end deftypefn
    function m = margin (this, X, Y, varargin)

      if (nargin < 3)
        error ("CompactClassificationECOC.margin: too few input arguments.");
      endif
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("CompactClassificationECOC.margin: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("CompactClassificationECOC.margin: X and Y must", ...
                       " have the same number of rows."));
      endif

      [~, NegLoss] = predict (this, X, varargin{:});
      ## marginsOf takes the number of score matrices, which is one here.
      m = marginsOf (NegLoss, gY, 1);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationECOC} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationECOC} {@var{e} =} edge (@dots{}, @var{name}, @var{value})
    ##
    ## Classification edge of a @code{CompactClassificationECOC}.
    ##
    ## @code{@var{e} = edge (@var{obj}, @var{X}, @var{Y})} returns the
    ## weighted mean of the margins, one number for the whole of @var{X}.
    ##
    ## It takes @qcode{'Weights'} beside the arguments @code{predict} takes.
    ##
    ## @seealso{CompactClassificationECOC.margin,
    ## CompactClassificationECOC.loss}
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      if (nargin < 3)
        error ("CompactClassificationECOC.edge: too few input arguments.");
      endif
      [W, rest] = splitWeights (this, varargin, Y, 'edge');
      m = margin (this, X, Y, rest{:});
      e = sum (W(:) .* m(:)) / sum (W);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationECOC} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationECOC} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of a @code{CompactClassificationECOC}.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} returns the
    ## weighted classification loss of @var{X} against the true labels
    ## @var{Y}.
    ##
    ## @multitable @columnfractions 0.28 0.02 0.7
    ## @headitem @var{Name} @tab @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab @tab @qcode{'classiferror'} (default),
    ## @qcode{'classifcost'}, @qcode{'mincost'}, @qcode{'binodeviance'},
    ## @qcode{'exponential'}, @qcode{'hinge'}, @qcode{'logit'} or
    ## @qcode{'quadratic'}.
    ##
    ## @item @qcode{'Weights'} @tab @tab One nonnegative weight per row of
    ## @var{X}.  The default is uniform.
    ## @end multitable
    ##
    ## It also takes the @qcode{'BinaryLoss'} and @qcode{'Decoding'}
    ## arguments @code{predict} takes.
    ##
    ## @seealso{CompactClassificationECOC.predict,
    ## CompactClassificationECOC.edge}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      if (nargin < 3)
        error ("CompactClassificationECOC.loss: too few input arguments.");
      endif
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("CompactClassificationECOC.loss: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("CompactClassificationECOC.loss: X and Y must have", ...
                       " the same number of rows."));
      endif

      LossFun = 'classiferror';
      keep = {};
      W = [];
      for i = 1:2:numel (varargin)
        switch (tolower (varargin{i}))
          case 'lossfun'
            LossFun = varargin{i+1};
            if (! (ischar (LossFun) && isrow (LossFun)))
              error (strcat ("CompactClassificationECOC.loss: 'LossFun'", ...
                             " must be a character vector."));
            endif
            LossFun = tolower (LossFun);
          case 'weights'
            W = varargin{i+1};
          otherwise
            keep(end+1:end+2) = varargin(i:i+1);
        endswitch
      endfor
      W = observationWeights (this, W, gY, 'loss');

      [~, NegLoss] = predict (this, X, keep{:});
      L = classificationLoss (LossFun, NegLoss, gY, W, this.Cost);

    endfunction

  endmethods

  methods (Access = private)

    ## The score each binary learner gives the class it calls +1, which is
    ## always the second of its two classes: the learners are fitted on -1
    ## and +1 and a learner lists its classes in order.
    function S = binaryScores (this, XC)

      L = numel (this.BinaryLearners);
      S = zeros (rows (XC), L);
      for j = 1:L
        [~, sc] = predict (this.BinaryLearners{j}, XC);
        S(:,j) = sc(:,2);
      endfor

    endfunction

    ## Pull 'Weights' out of a name-value list, leaving the rest for predict.
    function [W, rest] = splitWeights (this, args, Y, caller)

      rest = {};
      W = [];
      for i = 1:2:numel (args)
        if (strcmpi (args{i}, 'weights'))
          W = args{i+1};
        else
          rest(end+1:end+2) = args(i:i+1);
        endif
      endfor
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("CompactClassificationECOC.%s: %s", caller, errmsg);
      endif
      W = observationWeights (this, W, gY, caller);

    endfunction

    ## The weights an observation carries, defaulting to the class prior
    ## spread evenly over the rows of that class, as every learner here does.
    function W = observationWeights (this, W, gY, caller)

      n = numel (gY);
      if (isempty (W))
        W = ones (n, 1);
      endif
      if (! (isnumeric (W) && isvector (W) && numel (W) == n
             && all (W >= 0) && any (W > 0)))
        error (strcat ("CompactClassificationECOC.%s: 'Weights' must be a", ...
                       " nonnegative numeric vector with one element per", ...
                       " observation."), caller);
      endif
      W = priorNormalize (W(:), gY, this.Prior);

    endfunction

  endmethods

endclassdef

## Tests
%!test  # MATLAB parity: the property surface of a compact model
%! load fisheriris
%! CMdl = compact (fitcecoc (meas, species));
%! assert_equal (class (CMdl), 'CompactClassificationECOC');
%! assert_equal (numel (properties (CMdl)), 12);
%! assert_equal (CMdl.ClassNames, unique (species));
%! assert_equal (CMdl.CodingMatrix, [1, 1, 0; -1, 0, 1; 0, -1, -1]);

%!test  # MATLAB parity: the three outputs of predict
%! load fisheriris
%! CMdl = compact (fitcecoc (meas, species, 'Learners', 'tree'));
%! [label, NegLoss, PBScore] = predict (CMdl, meas([1, 51, 101], :));
%! assert_equal (label, {'setosa'; 'versicolor'; 'virginica'});
%! assert_equal (size (NegLoss), [3, 3]);
%! assert_equal (size (PBScore), [3, 3]);

%!test  # MATLAB parity: the decoding of a tree code, both schemes
%! ## Measured on R2024a, whose trees this package reproduces on this
%! ## fixture, so the whole path is compared and not only the labels.
%! load fisheriris
%! CMdl = compact (fitcecoc (meas, species, 'Learners', 'tree'));
%! r = [1, 20, 51, 70, 101, 130];
%! [~, W] = predict (CMdl, meas(r,:), 'Decoding', 'lossweighted');
%! assert_equal (W, ...
%!               [0, -1, -2; ...
%!                0, -1, -2; ...
%!                -2, 0, -1; ...
%!                -2, 0, -1; ...
%!                -2, -0.956994328922495, -0.000472589792060491; ...
%!                -2, -0.444444444444445, -0.111111111111111], 1e-12);
%! [~, B] = predict (CMdl, meas(r,:), 'Decoding', 'lossbased');
%! assert_equal (B(1,:), [-0.166666666666667, -0.833333333333333, ...
%!                        -1.5], 1e-12);

%!test  # lossweighted is the default decoding
%! load fisheriris
%! CMdl = compact (fitcecoc (meas, species));
%! [~, a] = predict (CMdl, meas(1:5,:));
%! [~, b] = predict (CMdl, meas(1:5,:), 'Decoding', 'lossweighted');
%! assert_equal (a, b);

%!test  # MATLAB parity: the edge is the weighted mean of the margins
%! load fisheriris
%! CMdl = compact (fitcecoc (meas, species, 'Learners', 'tree'));
%! m = margin (CMdl, meas, species);
%! assert_equal (edge (CMdl, meas, species), mean (m), 1e-12);

%!test  # loss counts the labels it got wrong
%! load fisheriris
%! CMdl = compact (fitcecoc (meas, species, 'Learners', 'tree'));
%! assert_equal (loss (CMdl, meas, species), 0.02, 1e-12);
%! assert_equal (loss (CMdl, meas, species, 'LossFun', 'classiferror'), ...
%!               0.02, 1e-12);

%!test  # a binary loss the learners cannot be read with is refused
%! load fisheriris
%! CMdl = compact (fitcecoc (meas, species));
%! assert_equal (CMdl.BinaryLoss, 'hinge');
%! CMdl.BinaryLoss = 'logit';
%! assert_equal (CMdl.BinaryLoss, 'logit');

## Test input validation
%!shared CMdl
%! CMdl = compact (fitcecoc (ones (4, 2), [1; 2; 1; 2]));
%!error<CompactClassificationECOC.predict: too few input arguments.> ...
%! predict (CMdl)
%!error<CompactClassificationECOC.predict: name-value arguments must be in pairs.> ...
%! predict (CMdl, ones (1, 2), 'Decoding')
%!error<CompactClassificationECOC.predict: 'Decoding' must be 'lossweighted' or 'lossbased'.> ...
%! predict (CMdl, ones (1, 2), 'Decoding', 'nosuch')
%!error<CompactClassificationECOC.predict: you cannot use 'quadratic' loss for binary learners with response in the range \(-Inf,\+Inf\).> ...
%! predict (CMdl, ones (1, 2), 'BinaryLoss', 'quadratic')
%!error<CompactClassificationECOC.subsasgn: you cannot use 'quadratic' loss for binary learners with response in the range \(-Inf,\+Inf\).> ...
%! CMdl.BinaryLoss = 'quadratic';
%!error<CompactClassificationECOC.margin: too few input arguments.> ...
%! margin (CMdl, ones (1, 2))
%!error<CompactClassificationECOC.loss: too few input arguments.> ...
%! loss (CMdl, ones (1, 2))
