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
##

classdef CompactClassificationNaiveBayes

  ## -*- texinfo -*-
  ## @deftp {statistics} CompactClassificationNaiveBayes
  ##
  ## Compact naive Bayes classification
  ##
  ## A @code{CompactClassificationNaiveBayes} object carries the fitted
  ## densities of a @code{ClassificationNaiveBayes} model and everything
  ## @code{predict} needs, but not the observations the model was fitted on.
  ## It classifies new data identically to the model it came from, and is far
  ## smaller to keep or to ship.
  ##
  ## Create one with the @code{compact} method of a
  ## @code{ClassificationNaiveBayes} object.  Because it holds no training
  ## data, it has no @code{resub} methods and cannot be cross-validated.
  ##
  ## @seealso{ClassificationNaiveBayes, fitcnb}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} DistributionNames
    ##
    ## Predictor distributions
    ##
    ## A cell array of character vectors with one entry per predictor, naming
    ## the distribution fitted to it.  This property is read-only.
    ##
    ## @end deftp
    DistributionNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} Mu
    ##
    ## Predictor means
    ##
    ## The means used to center the predictors, when the model standardizes
    ## them, and empty otherwise.  This property is read-only.
    ##
    ## @end deftp
    Mu = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} Sigma
    ##
    ## Predictor standard deviations
    ##
    ## The standard deviations used to scale the predictors, when the model
    ## standardizes them, and empty otherwise.  This property is read-only.
    ##
    ## @end deftp
    Sigma = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} DistributionParameters
    ##
    ## Fitted distribution parameters
    ##
    ## A cell array with one row per class and one column per predictor,
    ## holding the parameters fitted to each.  This property is read-only.
    ##
    ## @end deftp
    DistributionParameters = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} CategoricalLevels
    ##
    ## Levels of the categorical predictors
    ##
    ## A cell array with one entry per predictor, holding the distinct levels
    ## of each categorical predictor and empty for every other.  This property
    ## is read-only.
    ##
    ## @end deftp
    CategoricalLevels = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} Kernel
    ##
    ## Kernel smoothing functions
    ##
    ## A cell array naming the smoothing kernel of each kernel predictor, and
    ## empty for every other.  This property is read-only.
    ##
    ## @end deftp
    Kernel = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} Support
    ##
    ## Kernel smoothing supports
    ##
    ## A cell array giving the support of each kernel predictor's density, and
    ## empty for every other.  This property is read-only.
    ##
    ## @end deftp
    Support = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} Width
    ##
    ## Kernel smoothing bandwidths
    ##
    ## A numeric matrix with one row per class and one column per predictor,
    ## and empty when no predictor uses a kernel density.  This property is
    ## read-only.
    ##
    ## @end deftp
    Width = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} ClassNames
    ##
    ## Class labels of the fitted model
    ##
    ## The distinct classes the model was fitted on, in the order the other
    ## per-class properties use.  This property is read-only.
    ##
    ## @end deftp
    ClassNames = [];

  endproperties

  ## Properties a user may set after the model is fitted.  They sit between
  ## the read-only blocks so that 'properties' reports MATLAB's own order.
  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} Prior
    ##
    ## Class prior probabilities
    ##
    ## A numeric row vector with one entry per class, in the order of
    ## @qcode{ClassNames}, summing to one.
    ##
    ## @end deftp
    Prior = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} Cost
    ##
    ## Misclassification cost
    ##
    ## A square numeric matrix where @code{Cost(i,j)} is the cost of
    ## classifying an observation of class @math{i} into class @math{j}.
    ##
    ## @end deftp
    Cost = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} ScoreTransform
    ##
    ## Score transformation
    ##
    ## A character vector naming the function applied to the posterior
    ## returned by @code{predict}, or a function handle.
    ##
    ## @end deftp
    ScoreTransform = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} PredictorNames
    ##
    ## Predictor variable names
    ##
    ## A cell array of character vectors naming the predictors.  This property
    ## is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} CategoricalPredictors
    ##
    ## Categorical predictor indices
    ##
    ## The column indices treated as categorical, or empty when none is.  This
    ## property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector naming the response variable.  This property is
    ## read-only.
    ##
    ## @end deftp
    ResponseName = 'Y';

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationNaiveBayes} {property} ExpandedPredictorNames
    ##
    ## Expanded predictor variable names
    ##
    ## A cell array of character vectors naming the predictors as the model
    ## sees them.  This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

  endproperties

  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## The parsed ScoreTransform, applied to the posterior by predict.
    STfun = [];

  endproperties

  methods (Hidden)

    function this = set.Cost (this, Cost)
      [C, errmsg] = costMatrix (Cost, this.ClassNames);
      if (! isempty (errmsg))
        error ('CompactClassificationNaiveBayes.Cost: %s', errmsg);
      endif
      this.Cost = C;
    endfunction

    function this = set.Prior (this, Prior)
      P = Prior;
      if (isstruct (P))
        P = priorFromStruct (P, this.ClassNames, ...
                             'CompactClassificationNaiveBayes.Prior');
      elseif (ischar (P))
        if (strcmpi (P, 'uniform'))
          P = ones (1, rows (this.Cost)) / rows (this.Cost);
        elseif (strcmpi (P, 'empirical'))
          return;
        else
          error (strcat ("CompactClassificationNaiveBayes.Prior: a", ...
                         " character vector must be 'empirical' or", ...
                         " 'uniform'."));
        endif
      endif
      if (! (isnumeric (P) && isvector (P) && isreal (P)))
        error (strcat ("CompactClassificationNaiveBayes.Prior: must be a", ...
                       " real numeric vector, a structure, 'empirical', or", ...
                       " 'uniform'."));
      endif
      if (any (P < 0) || ! (sum (P) > 0))
        error (strcat ("CompactClassificationNaiveBayes.Prior: must be", ...
                       " nonnegative and must not be all zero."));
      endif
      this.Prior = P(:)' / sum (P);
    endfunction

    function this = set.ScoreTransform (this, val)
      [f, st] = parseScoreTransform (val, 'CompactClassificationNaiveBayes');
      this.ScoreTransform = st;
      this.STfun = f;
    endfunction

    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    function disp (this)
      fprintf ('\n  CompactClassificationNaiveBayes\n\n');
      fprintf ('%22s: %s\n', 'ResponseName', this.ResponseName);
      fprintf ('%22s: %s\n', 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ('%22s: %s\n', 'ClassNames', classNameListing (this.ClassNames));
      fprintf ('%22s: %s\n', 'ScoreTransform', this.ScoreTransform);
      fprintf ('%22s: %s\n', 'DistributionNames', ...
               classNameListing (this.DistributionNames));
      fprintf ('%22s: {%dx%d cell}\n', 'DistributionParameters', ...
               size (this.DistributionParameters));
      fprintf ('\n');
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationNaiveBayes} {@var{obj} =} CompactClassificationNaiveBayes (@var{Mdl})
    ##
    ## Create a @code{CompactClassificationNaiveBayes} object.
    ##
    ## @var{Mdl} is the @code{ClassificationNaiveBayes} object to compact.  The
    ## documented way to reach this constructor is the @code{compact} method.
    ##
    ## @seealso{ClassificationNaiveBayes}
    ## @end deftypefn
    function this = CompactClassificationNaiveBayes (Mdl)

      if (nargin < 1)
        error (strcat ("CompactClassificationNaiveBayes: too few", ...
                       " input arguments."));
      endif
      if (! isa (Mdl, 'ClassificationNaiveBayes'))
        error (strcat ("CompactClassificationNaiveBayes: MDL must be a", ...
                       " ClassificationNaiveBayes object."));
      endif

      this.DistributionNames      = Mdl.DistributionNames;
      this.Mu                     = Mdl.Mu;
      this.Sigma                  = Mdl.Sigma;
      this.DistributionParameters = Mdl.DistributionParameters;
      this.CategoricalLevels      = Mdl.CategoricalLevels;
      this.Kernel                 = Mdl.Kernel;
      this.Support                = Mdl.Support;
      this.Width                  = Mdl.Width;
      this.ClassNames             = Mdl.ClassNames;
      this.Cost                   = Mdl.Cost;
      this.Prior                  = Mdl.Prior;
      this.ScoreTransform         = Mdl.ScoreTransform;
      this.PredictorNames         = Mdl.PredictorNames;
      this.CategoricalPredictors  = Mdl.CategoricalPredictors;
      this.ResponseName           = Mdl.ResponseName;
      this.ExpandedPredictorNames = Mdl.ExpandedPredictorNames;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationNaiveBayes} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationNaiveBayes} {[@var{label}, @var{score}, @var{cost}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data with a compact naive Bayes model.
    ##
    ## The same classification the model it came from would give: the label of
    ## least expected cost, the posterior of each class, and the expected
    ## misclassification cost of each.
    ##
    ## @end deftypefn
    function [label, score, cost] = predict (this, XC)

      if (nargin < 2)
        error (strcat ("CompactClassificationNaiveBayes.predict: too few", ...
                       " input arguments."));
      endif
      if (isempty (XC))
        error ("CompactClassificationNaiveBayes.predict: XC is empty.");
      endif
      if (! (isnumeric (XC) && isreal (XC) && ismatrix (XC)))
        error (strcat ("CompactClassificationNaiveBayes.predict: XC must", ...
                       " be a real numeric matrix."));
      endif
      if (numel (this.PredictorNames) != columns (XC))
        error (strcat ("CompactClassificationNaiveBayes.predict: XC must", ...
                       " have the same number of predictors as the", ...
                       " trained model."));
      endif

      nCls = numel (this.Prior);
      logscore = nbLogLik (XC, this.DistributionNames, ...
                           this.DistributionParameters, nCls, ...
                           this.CategoricalLevels);
      logscore = logscore + log (this.Prior);
      ## An observation no class can account for, a categorical level the
      ## model never saw among them, leaves the likelihood at -Inf for every
      ## class alike.  Its posterior is the prior: what remains when the data
      ## says nothing.  Measured against R2024a, which returns the prior and
      ## not a uniform distribution, and does so for the whole observation
      ## even when its other predictors are perfectly informative.
      noinfo = all (logscore == -Inf, 2);
      logscore = logscore - max (logscore, [], 2);
      score = exp (logscore);
      score = score ./ sum (score, 2);
      score(isnan (score)) = 0;
      if (any (noinfo))
        score(noinfo,:) = repmat (this.Prior, sum (noinfo), 1);
      endif

      cost = score * this.Cost;
      [~, minIdx] = min (cost, [], 2);
      label = labelsFromIndex (this.ClassNames, minIdx);

      score = this.STfun (score);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationNaiveBayes} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin on new data.
    ##
    ## @end deftypefn
    function m = margin (this, X, Y)

      if (nargin < 3)
        error (strcat ("CompactClassificationNaiveBayes.margin: too few", ...
                       " input arguments."));
      endif
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("CompactClassificationNaiveBayes.margin: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("CompactClassificationNaiveBayes.margin: number of", ...
                       " rows in X and Y must be equal."));
      endif

      [~, s] = predict (this, X);
      m = marginsOf (s, gY, 1);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationNaiveBayes} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationNaiveBayes} {@var{e} =} edge (@dots{}, @qcode{'Weights'}, @var{w})
    ##
    ## Classification edge on new data.
    ##
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      if (nargin < 3)
        error (strcat ("CompactClassificationNaiveBayes.edge: too few", ...
                       " input arguments."));
      endif
      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       'CompactClassificationNaiveBayes', 'edge');
      e = sum (W .* margin (this, X, Y));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationNaiveBayes} {@var{l} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationNaiveBayes} {@var{l} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss on new data.
    ##
    ## Takes the @qcode{'LossFun'} and @qcode{'Weights'} options that
    ## @code{ClassificationNaiveBayes.loss} takes.
    ##
    ## @end deftypefn
    function l = loss (this, X, Y, varargin)

      if (nargin < 3)
        error (strcat ("CompactClassificationNaiveBayes.loss: too few", ...
                       " input arguments."));
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactClassificationNaiveBayes.loss: name-value", ...
                       " arguments must be in pairs."));
      endif

      LossFun = 'mincost';
      Weights = [];
      lf_opt = {'binodeviance', 'classifcost', 'classiferror', ...
                'exponential', 'hinge', 'logit', 'mincost', 'quadratic'};

      while (numel (varargin) > 0)
        Value = varargin{2};
        switch (tolower (varargin{1}))
          case 'lossfun'
            if (! (ischar (Value) && any (strcmpi (Value, lf_opt))))
              error (strcat ("CompactClassificationNaiveBayes.loss:", ...
                             " invalid loss function."));
            endif
            LossFun = tolower (Value);
          case 'weights'
            if (! (isnumeric (Value) && isvector (Value)))
              error (strcat ("CompactClassificationNaiveBayes.loss:", ...
                             " invalid 'Weights'."));
            endif
            Weights = Value;
          otherwise
            error (strcat ("CompactClassificationNaiveBayes.loss: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("CompactClassificationNaiveBayes.loss: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("CompactClassificationNaiveBayes.loss: number of", ...
                       " rows in X and Y must be equal."));
      endif
      if (isempty (Weights))
        w = ones (numel (gY), 1);
      else
        w = Weights(:);
        if (numel (w) != numel (gY))
          error (strcat ("CompactClassificationNaiveBayes.loss:", ...
                         " 'Weights' must have one element per", ...
                         " observation."));
        endif
      endif
      w = priorNormalize (w, gY, this.Prior);

      [~, s] = predict (this, X);
      l = classificationLoss (LossFun, s, gY, w, this.Cost);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationNaiveBayes} {@var{lp} =} logp (@var{obj}, @var{X})
    ##
    ## Log unconditional probability density of new data.
    ##
    ## @end deftypefn
    function lp = logp (this, X)

      if (nargin < 2)
        error (strcat ("CompactClassificationNaiveBayes.logp: too few", ...
                       " input arguments."));
      endif
      if (isempty (X))
        error ("CompactClassificationNaiveBayes.logp: X is empty.");
      endif
      if (numel (this.PredictorNames) != columns (X))
        error (strcat ("CompactClassificationNaiveBayes.logp: X must have", ...
                       " the same number of predictors as the trained", ...
                       " model."));
      endif

      nCls = numel (this.Prior);
      L = nbLogLik (X, this.DistributionNames, ...
                    this.DistributionParameters, nCls, ...
                    this.CategoricalLevels);
      L = L + log (this.Prior);
      ## Factoring out the largest term keeps a density small enough to
      ## underflow contributing its logarithm.  Where every class is
      ## impossible the largest term is -Inf and the factoring is 0/0, so the
      ## answer is written directly: the density really is zero there.
      Lmax = max (L, [], 2);
      lp = Lmax + log (sum (exp (L - Lmax), 2));
      lp(Lmax == -Inf) = -Inf;
      ## predict skips a missing predictor and classifies on the rest, but a
      ## density is not defined for an observation that is not fully observed,
      ## so this reports NaN where predict reports a class.  R2024a does the
      ## same.
      lp(any (isnan (X), 2)) = NaN;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationNaiveBayes} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a CompactClassificationNaiveBayes object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## CompactClassificationNaiveBayes object into an Octave binary file, the
    ## name of which is specified in @var{filename}, along with an extra
    ## variable, which defines the type classification object these variables
    ## constitute.  Use @code{loadmodel} in order to load a classification
    ## object into Octave's workspace.
    ##
    ## @seealso{loadmodel, fitcnb, CompactClassificationNaiveBayes}
    ## @end deftypefn
    function savemodel (this, fname)

      if (nargin < 2)
        error (strcat ("CompactClassificationNaiveBayes.savemodel:", ...
                       " too few input arguments."));
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("CompactClassificationNaiveBayes.savemodel:", ...
                       " FNAME must be a character vector."));
      endif

      ## Generate variable for class name
      classdef_name = 'CompactClassificationNaiveBayes';

      ## Create variables from model properties
      DistributionNames      = this.DistributionNames;
      Mu                     = this.Mu;
      Sigma                  = this.Sigma;
      CategoricalLevels      = this.CategoricalLevels;

      ## A kernel predictor's density is a classdef object, which Octave's
      ## save cannot serialize, so it goes out as the sample it was fitted to
      ## and load_model rebuilds it from that and the recorded bandwidth.
      DistributionParameters = nbKernelPack (this.DistributionParameters, ...
                                             this.DistributionNames);
      Kernel                 = this.Kernel;
      Support                = this.Support;
      Width                  = this.Width;
      ClassNames             = this.ClassNames;
      Prior                  = this.Prior;
      Cost                   = this.Cost;
      ScoreTransform         = this.ScoreTransform;
      PredictorNames         = this.PredictorNames;
      CategoricalPredictors  = this.CategoricalPredictors;
      ResponseName           = this.ResponseName;
      ExpandedPredictorNames = this.ExpandedPredictorNames;
      STfun                  = this.STfun;

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'DistributionNames', 'Mu', ...
            'Sigma', 'DistributionParameters', 'CategoricalLevels', ...
            'Kernel', 'Support', 'Width', 'ClassNames', 'Prior', 'Cost', ...
            'ScoreTransform', 'PredictorNames', 'CategoricalPredictors', ...
            'ResponseName', 'ExpandedPredictorNames', 'STfun');

    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      ## The compact model is built from a full one and has no training data
      ## of its own, so the smallest fit the full class accepts is compacted
      ## and then filled property by property.
      mdl = CompactClassificationNaiveBayes ( ...
                    ClassificationNaiveBayes ([0; 1; 2; 3], [1; 1; 2; 2]));

      ## Copy the saved data into the object.  Iterate over what was saved
      ## rather than over fieldnames (mdl): a private property such as STfun
      ## is written out by savemodel but is not reported by fieldnames, so
      ## comparing the two sets could never match and every load failed.
      ## Assignment is legal here because this is a method of the class
      ## itself.
      names = fieldnames (data);

      ## These three are assigned once everything else is in place, and in
      ## this order rather than the file's.  Cost comes before Prior because
      ## set.Prior counts the classes by the rows of Cost and not by
      ## ClassNames, so a prior arriving first is measured against the stub's.
      order = {'Cost', 'Prior', 'ScoreTransform'};
      late = ismember (names, order);
      tail = order(ismember (order, names));
      names = [names(! late); tail(:)];
      for i = 1:numel (names)
        try
          mdl.(names{i}) = data.(names{i});
        catch
          msg = strcat ("CompactClassificationNaiveBayes.load_model:", ...
                        " invalid model in '%s'.");
          error (msg, filename);
        end_try_catch
      endfor

      mdl.DistributionParameters = nbKernelUnpack ( ...
                    mdl.DistributionParameters, mdl.DistributionNames, ...
                    mdl.Kernel, mdl.Support, mdl.Width);

    endfunction

  endmethods

endclassdef

## Tests
%!test  # MATLAB parity: the surface a compact model reports
%! load fisheriris
%! CMdl = compact (fitcnb (meas, species));
%! assert_equal (class (CMdl), 'CompactClassificationNaiveBayes');
%! assert_equal (CMdl.ClassNames, unique (species));
%! assert_equal (CMdl.Prior, [1/3, 1/3, 1/3], 1e-15);
%! assert_equal (CMdl.Cost, [0, 1, 1; 1, 0, 1; 1, 1, 0]);
%! assert_equal (CMdl.ResponseName, 'Y');
%! assert_equal (CMdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (CMdl.DistributionNames, ...
%!               {'normal', 'normal', 'normal', 'normal'});
%! assert_equal (CMdl.ScoreTransform, 'none');

## A compact model keeps no observations, so it has neither the training data
## nor anything derived from it.
%!test  # the properties a compact model does not carry
%! load fisheriris
%! CMdl = compact (fitcnb (meas, species));
%! assert_equal (isprop (CMdl, 'X'), false);
%! assert_equal (isprop (CMdl, 'Y'), false);
%! assert_equal (isprop (CMdl, 'W'), false);
%! assert_equal (isprop (CMdl, 'NumObservations'), false);
%! assert_equal (isprop (CMdl, 'RowsUsed'), false);
%! assert_equal (isprop (CMdl, 'ModelParameters'), false);

%!test  # MATLAB parity: it classifies exactly as the model it came from
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! CMdl = compact (Mdl);
%! [ml, ms, mc] = predict (Mdl, meas);
%! [cl, cs, cc] = predict (CMdl, meas);
%! assert_equal (cl, ml);
%! assert_equal (cs, ms);
%! assert_equal (cc, mc);

%!test  # MATLAB parity: loss, edge, margin and logp
%! load fisheriris
%! CMdl = compact (fitcnb (meas, species));
%! assert_equal (loss (CMdl, meas, species), 0.04, 1e-14);
%! assert_equal (loss (CMdl, meas, species, 'LossFun', 'hinge'), ...
%!               0.052784701267562, 1e-12);
%! assert_equal (edge (CMdl, meas, species), 0.894430597464877, 1e-12);
%! assert_equal (sum (margin (CMdl, meas, species)), ...
%!               134.164589619731402, 1e-10);
%! assert_equal (logp (CMdl, meas)(1), 1.026591235856343, 1e-12);

%!test  # a kernel model compacts, densities and all
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'DistributionNames', 'kernel');
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.Kernel, {'normal', 'normal', 'normal', 'normal'});
%! assert_equal (CMdl.Width, Mdl.Width);
%! assert_equal (class (CMdl.DistributionParameters{1,1}), ...
%!               'prob.KernelDistribution');
%! assert_equal (predict (CMdl, meas), predict (Mdl, meas));

%!test  # Cost and Prior may be assigned on a compact model too
%! load fisheriris
%! CMdl = compact (fitcnb (meas, species));
%! CMdl.Cost = [0, 2, 2; 2, 0, 2; 2, 2, 0];
%! assert_equal (CMdl.Cost, [0, 2, 2; 2, 0, 2; 2, 2, 0]);
%! CMdl.Prior = [0.2, 0.3, 0.5];
%! assert_equal (CMdl.Prior, [0.2, 0.3, 0.5], 1e-15);
%! CMdl.ScoreTransform = 'logit';
%! assert_equal (CMdl.ScoreTransform, 'logit');

%!test  # a categorical model compacts, levels and all
%! X = [1, 2; 1, 3; 2, 2; 2, 3; 1, 2; 3, 1; 3, 3; 2, 1; 1, 1; 3, 2; ...
%!      2, 2; 1, 3; 3, 1; 2, 3; 1, 1; 3, 2; 2, 1; 1, 2; 3, 3; 2, 2];
%! Y = [repmat({'a'}, 10, 1); repmat({'b'}, 10, 1)];
%! Mdl = fitcnb (X, Y, 'DistributionNames', 'mvmn', ...
%!               'CategoricalPredictors', [1, 2]);
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.CategoricalLevels{1}, [1; 2; 3]);
%! assert_equal (CMdl.CategoricalPredictors, [1, 2]);
%! assert_equal (predict (CMdl, X), predict (Mdl, X));
%! assert_equal (loss (CMdl, X, Y), 0.45, 1e-14);

%!test  # a multinomial model compacts, and keeps its character DistributionNames
%! C = [2, 0, 1; 1, 3, 0; 0, 1, 4; 3, 1, 0; ...
%!      0, 2, 2; 1, 0, 3; 4, 1, 1; 0, 3, 2];
%! L = [repmat({'x'}, 4, 1); repmat({'y'}, 4, 1)];
%! CMdl = compact (fitcnb (C, L, 'DistributionNames', 'mn'));
%! assert_equal (CMdl.DistributionNames, 'mn');
%! [label, score] = predict (CMdl, [1, 1, 1; 4, 0, 0]);
%! assert_equal (label, {'x'; 'x'});
%! assert_equal (score(2,1), 0.769060987976952, 1e-12);

%!test  # an unseen level falls back to the prior on a compact model too
%! Z = [1, 1; 1, 1; 1, 2; 1, 2; 1, 1; 1, 2; 2, 2; 2, 1];
%! G = [repmat({'p'}, 6, 1); repmat({'q'}, 2, 1)];
%! CMdl = compact (fitcnb (Z, G, 'DistributionNames', 'mvmn', ...
%!                         'CategoricalPredictors', [1, 2]));
%! [~, score] = predict (CMdl, [3, 1]);
%! assert_equal (score, [0.75, 0.25], 1e-12);
%! assert_equal (logp (CMdl, [3, 1]), -Inf);

## A compact model survives savemodel and loadmodel, kernel densities and all,
## although it has no training data of its own to be rebuilt from.
%!test
%! load fisheriris
%! CMdl = compact (fitcnb (meas, species));
%! fname = tempname ();
%! savemodel (CMdl, fname);
%! C2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (C2), 'CompactClassificationNaiveBayes');
%! p = properties (CMdl);
%! for i = 1:numel (p)
%!   assert_equal (C2.(p{i}), CMdl.(p{i}));
%! endfor
%! assert_equal (predict (C2, meas(1:10,:)), predict (CMdl, meas(1:10,:)));

%!test
%! load fisheriris
%! CMdl = compact (fitcnb (meas, species, 'DistributionNames', 'kernel'));
%! fname = tempname ();
%! savemodel (CMdl, fname);
%! C2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (C2.DistributionParameters{1,1}), ...
%!               'prob.KernelDistribution');
%! assert_equal (C2.Width, CMdl.Width);
%! assert_equal (predict (C2, meas(1:10,:)), predict (CMdl, meas(1:10,:)));

%!error<CompactClassificationNaiveBayes.savemodel: too few input arguments.> ...
%! savemodel (compact (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2])))
%!error<CompactClassificationNaiveBayes.savemodel: FNAME must be a character vector.> ...
%! savemodel (compact (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2])), 1)
%!error<CompactClassificationNaiveBayes.savemodel: FNAME must be a character vector.> ...
%! savemodel (compact (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2])), ...
%!            ['ab'; 'cd'])

## Test input validation
%!error<CompactClassificationNaiveBayes: too few input arguments.> ...
%! CompactClassificationNaiveBayes ()
%!error<CompactClassificationNaiveBayes: MDL must be a ClassificationNaiveBayes object.> ...
%! CompactClassificationNaiveBayes (5)
%!error<CompactClassificationNaiveBayes.predict: XC is empty.> ...
%! predict (compact (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2])), [])
%!error<CompactClassificationNaiveBayes.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (compact (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2])), ones (2, 3))
%!error<CompactClassificationNaiveBayes.loss: too few input arguments.> ...
%! loss (compact (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2])), [1, 2; 2, 3; 3, 4; 4, 5])
%!error<CompactClassificationNaiveBayes.loss: invalid loss function.> ...
%! loss (compact (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2])), [1, 2; 2, 3; 3, 4; 4, 5], ...
%!       [1; 1; 2; 2], 'LossFun', 'nope')
%!error<CompactClassificationNaiveBayes.margin: too few input arguments.> ...
%! margin (compact (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2])), [1, 2; 2, 3; 3, 4; 4, 5])
%!error<CompactClassificationNaiveBayes.edge: too few input arguments.> ...
%! edge (compact (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2])), [1, 2; 2, 3; 3, 4; 4, 5])
%!error<CompactClassificationNaiveBayes.logp: X is empty.> ...
%! logp (compact (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2])), [])

## Every documented score transform reaches the scores that are reported, and
## none of them moves the label: a transform reshapes what is reported, not
## what is decided.
%!test
%! load fisheriris
%! Mdl = compact (fitcnb (meas, species));
%! Mdl.ScoreTransform = 'none';
%! [label, raw] = predict (Mdl, meas([1, 60, 120],:));
%! T = {'identity', @(x) x; 'doublelogit', @(x) 1 ./ (1 + exp (-2 * x)); ...
%!      'invlogit', @(x) log (x ./ (1 - x)); ...
%!      'logit', @(x) 1 ./ (1 + exp (-x)); ...
%!      'sign', @(x) sign (x); 'symmetric', @(x) 2 * x - 1; ...
%!      'symmetriclogit', @(x) 2 ./ (1 + exp (-x)) - 1};
%! for i = 1:rows (T)
%!   Mdl.ScoreTransform = T{i,1};
%!   [l, s] = predict (Mdl, meas([1, 60, 120],:));
%!   assert_equal (s, T{i,2}(raw), 1e-12);
%!   assert_equal (l, label);
%! endfor
%! ## ismax marks the largest score of each observation, ties to the first.
%! [~, k] = max (raw, [], 2);
%! e = zeros (size (raw));
%! e(sub2ind (size (raw), (1:rows (raw))', k)) = 1;
%! Mdl.ScoreTransform = 'ismax';
%! [~, s] = predict (Mdl, meas([1, 60, 120],:));
%! assert_equal (s, e);
%! Mdl.ScoreTransform = 'symmetricismax';
%! [~, s] = predict (Mdl, meas([1, 60, 120],:));
%! assert_equal (s, 2 * e - 1);

## A function handle is taken as given and applied to the scores.
%!test
%! load fisheriris
%! Mdl = compact (fitcnb (meas, species));
%! Mdl.ScoreTransform = 'none';
%! [label, raw] = predict (Mdl, meas([1, 60, 120],:));
%! Mdl.ScoreTransform = @(x) x .^ 2;
%! [l, s] = predict (Mdl, meas([1, 60, 120],:));
%! assert_equal (s, raw .^ 2, 1e-12);
%! assert_equal (l, label);
