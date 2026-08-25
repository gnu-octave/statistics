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

classdef ClassificationNaiveBayes

  ## -*- texinfo -*-
  ## @deftp {statistics} ClassificationNaiveBayes
  ##
  ## Naive Bayes classification
  ##
  ## The @code{ClassificationNaiveBayes} class implements a naive Bayes
  ## classifier object, which can predict responses for new data using the
  ## @code{predict} method.
  ##
  ## A naive Bayes classifier estimates one univariate density per class and
  ## per predictor, and treats the predictors as conditionally independent
  ## given the class.  The joint likelihood of an observation is therefore the
  ## product of its per-predictor densities, and the posterior follows from
  ## the class prior by Bayes' rule.  The independence assumption is rarely
  ## true, but it costs only one density per predictor rather than one joint
  ## density over all of them, which is what makes the model usable when the
  ## predictors are many and the observations few.
  ##
  ## Create a @code{ClassificationNaiveBayes} object by using the
  ## @code{fitcnb} function or the class constructor.
  ##
  ## Each predictor carries its own distribution, named in
  ## @qcode{DistributionNames}, and the fitted parameters of class @math{k} and
  ## predictor @math{j} are held in @code{DistributionParameters@{k,j@}}.  A
  ## @qcode{'normal'} predictor stores a two element column vector, the class
  ## conditional mean and standard deviation; a @qcode{'kernel'} predictor
  ## stores a @code{prob.KernelDistribution} object.
  ##
  ## @seealso{fitcnb}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} Y
    ##
    ## Class labels
    ##
    ## Specified as a logical or numeric column vector, or as a character array
    ## or a cell array of character vectors with the same number of rows as the
    ## predictor data.  Each row in @var{Y} is the observed class label for the
    ## corresponding row in @var{X}.  This property is read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix containing the predictor data.  Each column of @var{X}
    ## represents one predictor (variable), and each row represents one
    ## observation.  This property is read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} RowsUsed
    ##
    ## Rows used for fitting
    ##
    ## A logical column vector with the same length as the observations in the
    ## original predictor data @var{X}, true for each row that was used for
    ## fitting the model.  It is empty, @qcode{[]}, when every observation was
    ## used, so a non-empty value means that rows holding missing values were
    ## dropped.  This property is read-only.
    ##
    ## @end deftp
    RowsUsed = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector with one entry per observation used for fitting,
    ## summing to one.  Each class contributes its prior, spread evenly over
    ## the observations belonging to it.  This property is read-only.
    ##
    ## @end deftp
    W = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} ModelParameters
    ##
    ## What was fitted, and how
    ##
    ## A structure carrying @qcode{DistributionNames}, @qcode{Kernel},
    ## @qcode{Support}, @qcode{Width}, @qcode{StandardizeData},
    ## @qcode{Version}, @qcode{Method} and @qcode{Type}.
    ##
    ## It records the arguments as they were @emph{given}, where the
    ## properties of the same name record what they were resolved to: a model
    ## fitted with no @qcode{'DistributionNames'} argument reports the single
    ## name @qcode{'normal'} here and one name per predictor there.  The
    ## kernel settings are filled in with their defaults when a kernel density
    ## was asked for, and left empty when none was.  This property is
    ## read-only.
    ##
    ## @end deftp
    ModelParameters = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer specifying the number of observations used to train
    ## the model, after any row holding a missing value has been dropped.  This
    ## property is read-only.
    ##
    ## @end deftp
    NumObservations = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} BinEdges
    ##
    ## Bin edges
    ##
    ## A cell array with one entry per predictor, holding that predictor's bin
    ## edges where the learner discretized it before fitting.  A naive Bayes
    ## model fits a density to each predictor as it stands and bins nothing, so
    ## this is always an empty cell.  It is kept because the cross-validated
    ## model carries it across, and because code that reaches into it with
    ## @code{cellfun} must find a cell rather than an empty matrix.  This
    ## property is read-only.
    ##
    ## @end deftp
    BinEdges = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} PredictorNames
    ##
    ## Predictor variable names
    ##
    ## A cell array of character vectors naming the predictors, in the order in
    ## which they appear in @var{X}.  The default names are @qcode{'x1'},
    ## @qcode{'x2'}, and so on.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} CategoricalPredictors
    ##
    ## Categorical predictor indices
    ##
    ## A numeric row vector of the column indices of @var{X} treated as
    ## categorical, or empty when none is.  This property is read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector naming the response variable, @qcode{'Y'} by default.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName = 'Y';

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} ExpandedPredictorNames
    ##
    ## Expanded predictor variable names
    ##
    ## A cell array of character vectors naming the predictors as the model
    ## sees them.  It equals @qcode{PredictorNames} unless a categorical
    ## predictor has been expanded into indicator variables.  This property is
    ## read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} ClassNames
    ##
    ## Class labels of the fitted model
    ##
    ## A cell array of character vectors, a logical or numeric column vector,
    ## or a character array, holding the distinct classes the model was fitted
    ## on, in the order the other per-class properties use.  This property is
    ## read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} HyperparameterOptimizationResults
    ##
    ## Results of the hyperparameter optimization
    ##
    ## @strong{Always empty.}  It is declared for MATLAB compatibility, where
    ## it holds what an automatic search over the hyperparameters found.  This
    ## class fits the parameters it is given and runs no such search, so there
    ## is nothing to report.  This property is read-only.
    ##
    ## @end deftp
    HyperparameterOptimizationResults = [];

  endproperties

  ## Properties a user may set after the model is fitted.  Each one is
  ## validated by its set method below.  They sit between the two read-only
  ## blocks so that 'properties' reports MATLAB's own order.
  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} Prior
    ##
    ## Class prior probabilities
    ##
    ## A numeric row vector with one entry per class, in the order of
    ## @qcode{ClassNames}, summing to one.  It may be assigned after fitting,
    ## as a numeric vector, as a structure carrying @qcode{ClassNames} and
    ## @qcode{ClassProbs}, or as @qcode{'empirical'} or @qcode{'uniform'}.
    ## Assigning it re-derives @qcode{W}.
    ##
    ## @end deftp
    Prior = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} Cost
    ##
    ## Misclassification cost
    ##
    ## A square numeric matrix with one row and column per class, in the order
    ## of @qcode{ClassNames}.  @code{Cost(i,j)} is the cost of classifying an
    ## observation of class @math{i} into class @math{j}, and the default is
    ## one off the diagonal and zero on it.  It may be assigned after fitting,
    ## as a matrix or as a structure carrying @qcode{ClassNames} and
    ## @qcode{ClassificationCosts}.
    ##
    ## @end deftp
    Cost = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} ScoreTransform
    ##
    ## Score transformation
    ##
    ## A character vector naming the function applied to the posterior returned
    ## by @code{predict}, or a function handle taking and returning a matrix of
    ## the same size.  The default, @qcode{'none'}, leaves the posterior
    ## untouched.
    ##
    ## @end deftp
    ScoreTransform = 'none';

  endproperties

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} DistributionNames
    ##
    ## Predictor distributions
    ##
    ## A cell array of character vectors with one entry per predictor, naming
    ## the distribution fitted to it: @qcode{'normal'} or @qcode{'kernel'}.
    ## This property is read-only.
    ##
    ## @end deftp
    DistributionNames = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} Mu
    ##
    ## Predictor means
    ##
    ## The means used to center the predictors, when the model standardizes
    ## them, and empty otherwise.  These are @emph{not} the class conditional
    ## means, which are held in @qcode{DistributionParameters}.  This property
    ## is read-only.
    ##
    ## @end deftp
    Mu = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} Sigma
    ##
    ## Predictor standard deviations
    ##
    ## The standard deviations used to scale the predictors, when the model
    ## standardizes them, and empty otherwise.  These are @emph{not} the class
    ## conditional standard deviations, which are held in
    ## @qcode{DistributionParameters}.  This property is read-only.
    ##
    ## @end deftp
    Sigma = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} DistributionParameters
    ##
    ## Fitted distribution parameters
    ##
    ## A cell array with one row per class and one column per predictor.
    ## @code{DistributionParameters@{k,j@}} holds the parameters fitted to
    ## predictor @math{j} within class @math{k}: a two element column vector,
    ## the mean and the standard deviation, for a @qcode{'normal'} predictor,
    ## and a @code{prob.KernelDistribution} object for a @qcode{'kernel'} one.
    ## This property is read-only.
    ##
    ## @end deftp
    DistributionParameters = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} CategoricalLevels
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
    ## @deftp {ClassificationNaiveBayes} {property} Kernel
    ##
    ## Kernel smoothing functions
    ##
    ## A cell array with one entry per predictor naming the smoothing kernel
    ## used by a @qcode{'kernel'} predictor, and empty for every other.  This
    ## property is read-only.
    ##
    ## @end deftp
    Kernel = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} Support
    ##
    ## Kernel smoothing supports
    ##
    ## A cell array with one entry per predictor giving the support of a
    ## @qcode{'kernel'} predictor's density, and empty for every other.  This
    ## property is read-only.
    ##
    ## @end deftp
    Support = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationNaiveBayes} {property} Width
    ##
    ## Kernel smoothing bandwidths
    ##
    ## A numeric matrix with one row per class and one column per predictor,
    ## giving the bandwidth of each @qcode{'kernel'} predictor's density, and
    ## empty when no predictor uses one.  This property is read-only.
    ##
    ## @end deftp
    Width = [];

  endproperties

  ## Readable by the compact counterpart, which copies it, and kept out of the
  ## documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## The parsed ScoreTransform, applied to the posterior by predict.
    STfun = [];

  endproperties

  ## Set methods for the properties a user may assign after fitting.
  methods (Hidden)

    function this = set.Cost (this, Cost)
      [C, errmsg] = costMatrix (Cost, this.ClassNames);
      if (! isempty (errmsg))
        error ('ClassificationNaiveBayes.Cost: %s', errmsg);
      endif
      this.Cost = C;
    endfunction

    function this = set.Prior (this, Prior)
      P = Prior;
      if (isstruct (P))
        P = priorFromStruct (P, this.ClassNames, ...
                             'ClassificationNaiveBayes.Prior');
      elseif (ischar (P))
        if (strcmpi (P, 'uniform'))
          P = ones (1, rows (this.Cost)) / rows (this.Cost);
        elseif (! strcmpi (P, 'empirical'))
          error (strcat ("ClassificationNaiveBayes.Prior: a character", ...
                         " vector must be 'empirical' or 'uniform'."));
        else
          return;   # 'empirical' after fitting is what is already stored
        endif
      endif
      if (! (isnumeric (P) && isvector (P) && isreal (P)))
        error (strcat ("ClassificationNaiveBayes.Prior: must be a real", ...
                       " numeric vector, a structure, 'empirical', or", ...
                       " 'uniform'."));
      endif
      if (numel (P) != rows (this.Cost))
        error (strcat ("ClassificationNaiveBayes.Prior: must have one", ...
                       " element per class."));
      endif
      if (any (P < 0) || ! (sum (P) > 0))
        error (strcat ("ClassificationNaiveBayes.Prior: must be", ...
                       " nonnegative and must not be all zero."));
      endif
      this.Prior = P(:)' / sum (P);
      ## The weights follow the prior, so reassigning one re-derives the other.
      if (! isempty (this.Y))
        gY = labelIndices (this.ClassNames, this.Y);
        keep = gY > 0;
        this.W = priorWeights (this.Prior, gY(keep), sum (keep));
      endif
    endfunction

    function this = set.ScoreTransform (this, val)
      [f, st] = parseScoreTransform (val, 'ClassificationNaiveBayes');
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
      fprintf ('\n  ClassificationNaiveBayes\n\n');
      fprintf ('%22s: %s\n', 'ResponseName', this.ResponseName);
      fprintf ('%22s: %s\n', 'CategoricalPredictors', ...
               mat2str (this.CategoricalPredictors));
      fprintf ('%22s: %s\n', 'ClassNames', classNameListing (this.ClassNames));
      fprintf ('%22s: %s\n', 'ScoreTransform', this.ScoreTransform);
      fprintf ('%22s: %d\n', 'NumObservations', this.NumObservations);
      fprintf ('%22s: %s\n', 'DistributionNames', ...
               classNameListing (this.DistributionNames));
      fprintf ('%22s: {%dx%d cell}\n', 'DistributionParameters', ...
               size (this.DistributionParameters));
      fprintf ('\n');
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNaiveBayes} {@var{obj} =} ClassificationNaiveBayes (@var{X}, @var{Y})
    ## @deftypefnx {ClassificationNaiveBayes} {@var{obj} =} ClassificationNaiveBayes (@dots{}, @var{name}, @var{value})
    ##
    ## Create a @code{ClassificationNaiveBayes} object.
    ##
    ## @code{@var{obj} = ClassificationNaiveBayes (@var{X}, @var{Y})} fits a
    ## naive Bayes classifier to the predictor data @var{X} and the class
    ## labels @var{Y}.  The supported @qcode{Name}/@qcode{Value} pairs are
    ## those of @code{fitcnb}, which is the documented way to reach this
    ## constructor.
    ##
    ## @seealso{fitcnb}
    ## @end deftypefn
    function this = ClassificationNaiveBayes (X, Y, varargin)

      ## Check for appropriate number of input arguments
      if (nargin < 2)
        error ("ClassificationNaiveBayes: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationNaiveBayes: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      ## Validate X
      if (! (isnumeric (X) && isreal (X) && ismatrix (X) && ! isempty (X)))
        error (strcat ("ClassificationNaiveBayes: X must be a", ...
                       " non-empty real numeric matrix."));
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error (strcat ("ClassificationNaiveBayes: number of rows in", ...
                       " X and Y must be equal."));
      endif

      nPred = columns (X);

      ## Set default values before parsing optional parameters
      CatPreds       = [];
      ClassNames     = [];
      Cost           = [];
      DistNames      = [];
      Kernel         = [];
      PredictorNames = {};
      Prior          = 'empirical';
      ResponseName   = 'Y';
      ScoreTransform = 'none';
      Support        = [];
      Width          = [];

      ## Parse optional parameters
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("ClassificationNaiveBayes: 'PredictorNames'", ...
                             " must be supplied as a cellstring array."));
            elseif (numel (PredictorNames) != nPred)
              error (strcat ("ClassificationNaiveBayes: 'PredictorNames'", ...
                             " must equal the number of columns in X."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! (ischar (ResponseName) && isrow (ResponseName)))
              error (strcat ("ClassificationNaiveBayes: 'ResponseName'", ...
                             " must be a character vector."));
            endif

          case 'classnames'
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames)
                   || islogical (ClassNames) || ischar (ClassNames)))
              error (strcat ("ClassificationNaiveBayes: 'ClassNames'", ...
                             " must be a cell array of character vectors,", ...
                             " a logical vector, a numeric vector, or a", ...
                             " character array."));
            endif

          case 'prior'
            Prior = varargin{2};

          case 'cost'
            Cost = varargin{2};

          case 'scoretransform'
            ScoreTransform = varargin{2};

          case 'categoricalpredictors'
            CatPreds = varargin{2};

          case 'distributionnames'
            DistNames = varargin{2};

          case 'kernel'
            Kernel = varargin{2};

          case 'support'
            Support = varargin{2};

          case 'width'
            Width = varargin{2};

          otherwise
            error (strcat ("ClassificationNaiveBayes: invalid parameter", ...
                           sprintf (" name '%s'.", varargin{1})));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## Resolve the categorical predictors and the distribution of every
      ## predictor.  The two are tied: naming a predictor categorical makes
      ## it multivariate multinomial unless a distribution is named for it,
      ## and naming that distribution makes the predictor categorical, which
      ## MATLAB warns about rather than silently accepting.
      catidx = nbCategorical (CatPreds, nPred, 'ClassificationNaiveBayes');
      if (isempty (DistNames) && ! isempty (catidx))
        DistNames = repmat ({'normal'}, 1, nPred);
        DistNames(catidx) = {'mvmn'};
      endif
      D = nbDistNames (DistNames, nPred, 'ClassificationNaiveBayes');
      if (iscell (D))
        mvidx = find (strcmp (D, 'mvmn'));
        if (! all (ismember (mvidx, catidx)))
          warning (strcat ("ClassificationNaiveBayes: the 'mvmn'", ...
                           " distribution was named for a predictor that is", ...
                           " not in 'CategoricalPredictors', which is", ...
                           " updated to include every 'mvmn' predictor."));
          catidx = sort (unique ([catidx, mvidx]));
        endif
      endif
      this.DistributionNames = D;
      this.CategoricalPredictors = catidx;

      ## Store the data as supplied, then drop the rows that are not complete
      this.X = X;
      this.Y = Y;
      ok = ! any (isnan (X), 2);
      if (! all (ok))
        this.RowsUsed = ok;
      endif
      Xf = X(ok, :);
      Yf = Y(ok, :);
      if (isempty (Xf))
        error (strcat ("ClassificationNaiveBayes: no observation is", ...
                       " free of missing values."));
      endif

      ## Resolve the classes.  Naming them keeps only those observations, so
      ## the rows of the others are dropped exactly as a missing value is.
      if (isempty (ClassNames))
        C = uniqueLabels (Yf);
      else
        C = ClassNames;
        if (ischar (C) && isrow (C))
          C = cellstr (C);
        endif
        ## The classes are held one per row, whichever orientation they were
        ## given in: every per-class property is counted by rows.
        if (! ischar (C))
          C = C(:);
        endif
        if (! labelsKnown (C, uniqueLabels (Yf)))
          error (strcat ("ClassificationNaiveBayes: not all 'ClassNames'", ...
                         " are present in Y."));
        endif
      endif

      ## A label outside the classes indexes as zero.  That is a fault only
      ## when the classes came from the response itself; naming a subset of
      ## them is a request to keep those observations and drop the rest.
      [gY, errmsg] = labelIndices (C, Yf);
      if (! isempty (errmsg) && isempty (ClassNames))
        error ('ClassificationNaiveBayes: %s', errmsg);
      endif
      inC = gY > 0;
      if (! all (inC))
        Xf = Xf(inC, :);
        Yf = Yf(inC, :);
        gY = gY(inC);
        ok(ok) = inC;
        this.RowsUsed = ok;
        if (isempty (Xf))
          error (strcat ("ClassificationNaiveBayes: no observation", ...
                         " belongs to the named classes."));
        endif
      endif

      this.ClassNames = C;
      nObs = rows (Xf);
      nCls = rows (C);

      this.NumObservations = nObs;
      this.ResponseName = ResponseName;
      if (isempty (PredictorNames))
        PredictorNames = arrayfun (@(k) sprintf ('x%d', k), 1:nPred, ...
                                   'UniformOutput', false);
      endif
      this.PredictorNames = PredictorNames;
      this.ExpandedPredictorNames = PredictorNames;

      ## Cost first: the Prior set method reads its size to count the classes
      if (isempty (Cost))
        this.Cost = ones (nCls) - eye (nCls);
      else
        this.Cost = Cost;
      endif

      ## Prior, and the observation weights that follow from it
      if (ischar (Prior) && strcmpi (Prior, 'empirical'))
        this.Prior = accumarray (gY, 1, [nCls, 1])' / nObs;
      elseif (ischar (Prior) && strcmpi (Prior, 'uniform'))
        this.Prior = ones (1, nCls) / nCls;
      else
        this.Prior = Prior;
      endif

      this.ScoreTransform = ScoreTransform;

      ## Fit one density per class and per predictor
      [this.DistributionParameters, this.Kernel, this.Support, ...
       this.Width, this.CategoricalLevels] = ...
              nbFit (Xf, gY, nCls, this.DistributionNames, this.W, ...
                     Kernel, Support, Width, 'ClassificationNaiveBayes', ...
                     this.ClassNames, this.PredictorNames);

      ## The arguments as they were given.  The fields are assigned one at a
      ## time: a cell handed to struct() would make a structure array of its
      ## elements rather than one structure holding the cell.
      mp = struct ();
      if (isempty (DistNames))
        mp.DistributionNames = 'normal';
      else
        mp.DistributionNames = DistNames;
      endif
      ## The kernel settings carry their defaults only where a kernel density
      ## was actually asked for, and stay empty otherwise.
      if (iscell (this.DistributionNames)
          && any (strcmp (this.DistributionNames, 'kernel')))
        mp.Kernel = nbGivenOr (Kernel, 'normal');
        mp.Support = nbGivenOr (Support, 'unbounded');
        mp.Width = nbGivenOr (Width, NaN);
        mp.StandardizeData = 0;
      else
        mp.Kernel = [];
        mp.Support = [];
        mp.Width = [];
        mp.StandardizeData = [];
      endif
      mp.Version = 1;
      mp.Method = 'NaiveBayes';
      mp.Type = 'classification';
      this.ModelParameters = mp;

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNaiveBayes} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationNaiveBayes} {[@var{label}, @var{score}, @var{cost}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data with a trained @code{ClassificationNaiveBayes} object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the predicted
    ## class label for each row of @var{XC}, which must have as many columns as
    ## the predictor data the model was fitted on.
    ##
    ## @code{[@var{label}, @var{score}, @var{cost}] = predict (@var{obj},
    ## @var{XC})} also returns @var{score}, the posterior probability of each
    ## class, and @var{cost}, the expected misclassification cost of assigning
    ## each observation to each class.  The label of an observation is the
    ## class of least expected cost.
    ##
    ## @end deftypefn
    function [label, score, cost] = predict (this, XC)

      if (nargin < 2)
        error ("ClassificationNaiveBayes.predict: too few input arguments.");
      endif
      if (isempty (XC))
        error ("ClassificationNaiveBayes.predict: XC is empty.");
      endif
      if (! (isnumeric (XC) && isreal (XC) && ismatrix (XC)))
        error (strcat ("ClassificationNaiveBayes.predict: XC must be a", ...
                       " real numeric matrix."));
      endif
      if (columns (this.X) != columns (XC))
        error (strcat ("ClassificationNaiveBayes.predict: XC must have", ...
                       " the same number of predictors as the trained", ...
                       " model."));
      endif

      nCls = numel (this.Prior);

      ## Work in logs: the product over predictors underflows for even a
      ## moderate number of them, and the posterior only needs the differences.
      logscore = nbLogLik (XC, this.DistributionNames, ...
                           this.DistributionParameters, nCls, ...
                           this.CategoricalLevels);
      logscore = logscore + log (this.Prior);

      ## Subtracting the row maximum keeps a well separated observation a
      ## posterior rather than a ratio of two underflowed zeros.
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

      ## Expected misclassification cost, and the label of least cost
      cost = score * this.Cost;
      [~, minIdx] = min (cost, [], 2);
      label = labelsFromIndex (this.ClassNames, minIdx);

      ## The transform is applied once, to the assembled posterior
      score = this.STfun (score);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNaiveBayes} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationNaiveBayes} {@var{CVMdl} =} crossval (@dots{}, @var{name}, @var{value})
    ##
    ## Cross-validate a trained naive Bayes model.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} partitions the training data
    ## into ten folds, or into as many folds as there are observations when
    ## there are fewer than ten, refits the model on each fold's training part
    ## and returns a @code{ClassificationPartitionedModel}.
    ##
    ## @code{@var{CVMdl} = crossval (@dots{}, @var{name}, @var{value})} takes
    ## exactly one of @qcode{'KFold'}, @qcode{'Holdout'}, @qcode{'Leaveout'}
    ## or @qcode{'CVPartition'}.
    ##
    ## @seealso{ClassificationPartitionedModel, cvpartition}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)

      if (nargin < 1)
        error (strcat ("ClassificationNaiveBayes.crossval: too few", ...
                       " input arguments."));
      endif
      if (numel (varargin) == 1)
        error (strcat ("ClassificationNaiveBayes.crossval: Name-Value", ...
                       " arguments must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("ClassificationNaiveBayes.crossval: specify only", ...
                       " one of the optional Name-Value paired arguments."));
      endif

      if (this.NumObservations < 10)
        numFolds = this.NumObservations;
      else
        numFolds = 10;
      endif
      Holdout     = [];
      Leaveout    = 'off';
      CVPartition = [];

      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))

          case 'kfold'
            numFolds = varargin{2};
            if (! (isnumeric (numFolds) && isscalar (numFolds)
                   && (numFolds == fix (numFolds)) && numFolds > 1))
              error (strcat ("ClassificationNaiveBayes.crossval: 'KFold'", ...
                             " must be an integer value greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("ClassificationNaiveBayes.crossval:", ...
                             " 'Holdout' must be a numeric value between", ...
                             " 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && any (strcmpi (Leaveout, {'on', 'off'}))))
              error (strcat ("ClassificationNaiveBayes.crossval:", ...
                             " 'Leaveout' must be either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("ClassificationNaiveBayes.crossval:", ...
                             " 'CVPartition' must be a 'cvpartition'", ...
                             " object."));
            endif

          otherwise
            error (strcat ("ClassificationNaiveBayes.crossval: invalid", ...
                           " parameter name in optional paired arguments."));

        endswitch
        varargin(1:2) = [];
      endwhile

      ## The partition covers the observations actually trained on: a row
      ## dropped for a missing value is not one the folds can use.  The
      ## response is passed rather than a count so the folds stay stratified.
      [~, Yused] = nbTrainX (this);
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (Yused, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (numel (Yused), 'LeaveOut');
      else
        partition = cvpartition (Yused, 'KFold', numFolds);
      endif

      CVMdl = ClassificationPartitionedModel (this, partition);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationNaiveBayes} {@var{CMdl} =} compact (@var{obj})
    ##
    ## Drop the training data from a trained model.
    ##
    ## @code{@var{CMdl} = compact (@var{obj})} returns a
    ## @code{CompactClassificationNaiveBayes} object carrying the fitted
    ## densities and everything @code{predict} needs, but not the observations
    ## the model was fitted on.  It classifies new data identically and is far
    ## smaller to keep or to ship.
    ##
    ## @seealso{CompactClassificationNaiveBayes}
    ## @end deftypefn
    function CMdl = compact (this)

      CMdl = CompactClassificationNaiveBayes (this);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationNaiveBayes} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margin on new data.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns one margin
    ## per observation: the posterior the model gives the observation's true
    ## class, less the largest posterior it gives any other class.  A positive
    ## margin means the observation is classified correctly, and a larger one
    ## means it is classified more confidently.
    ##
    ## @end deftypefn
    function m = margin (this, X, Y)

      if (nargin < 3)
        error ("ClassificationNaiveBayes.margin: too few input arguments.");
      endif
      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("ClassificationNaiveBayes.margin: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("ClassificationNaiveBayes.margin: number of rows", ...
                       " in X and Y must be equal."));
      endif

      [~, s] = predict (this, X);
      m = marginsOf (s, gY, 1);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNaiveBayes} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationNaiveBayes} {@var{e} =} edge (@dots{}, @qcode{'Weights'}, @var{w})
    ##
    ## Classification edge on new data.
    ##
    ## @code{@var{e} = edge (@var{obj}, @var{X}, @var{Y})} returns the weighted
    ## mean of the margins, a single number summarising how confidently the
    ## model classifies the data.
    ##
    ## The weights are normalized within each class to that class's prior
    ## before they are applied.
    ##
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      if (nargin < 3)
        error ("ClassificationNaiveBayes.edge: too few input arguments.");
      endif
      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       'ClassificationNaiveBayes', 'edge');
      e = sum (W .* margin (this, X, Y));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNaiveBayes} {@var{l} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationNaiveBayes} {@var{l} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss on new data.
    ##
    ## @code{@var{l} = loss (@var{obj}, @var{X}, @var{Y})} returns the minimum
    ## expected misclassification cost.
    ##
    ## @code{@var{l} = loss (@dots{}, @var{name}, @var{value})} takes the
    ## following options.
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab One of @qcode{'binodeviance'},
    ## @qcode{'classifcost'}, @qcode{'classiferror'}, @qcode{'exponential'},
    ## @qcode{'hinge'}, @qcode{'logit'}, @qcode{'mincost'} (default) or
    ## @qcode{'quadratic'}.
    ##
    ## @item @qcode{'Weights'} @tab A numeric vector of observation weights,
    ## one per row of @var{X}.
    ##
    ## @end multitable
    ##
    ## @end deftypefn
    function l = loss (this, X, Y, varargin)

      if (nargin < 3)
        error ("ClassificationNaiveBayes.loss: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationNaiveBayes.loss: name-value", ...
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
              error (strcat ("ClassificationNaiveBayes.loss: invalid", ...
                             " loss function."));
            endif
            LossFun = tolower (Value);
          case 'weights'
            if (! (isnumeric (Value) && isvector (Value)))
              error ("ClassificationNaiveBayes.loss: invalid 'Weights'.");
            endif
            Weights = Value;
          otherwise
            error (strcat ("ClassificationNaiveBayes.loss: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      [gY, errmsg] = labelIndices (this.ClassNames, Y);
      if (! isempty (errmsg))
        error ("ClassificationNaiveBayes.loss: %s", errmsg);
      endif
      if (rows (X) != numel (gY))
        error (strcat ("ClassificationNaiveBayes.loss: number of rows in", ...
                       " X and Y must be equal."));
      endif
      if (isempty (Weights))
        w = ones (numel (gY), 1);
      else
        w = Weights(:);
        if (numel (w) != numel (gY))
          error (strcat ("ClassificationNaiveBayes.loss: 'Weights' must", ...
                         " have one element per observation."));
        endif
      endif

      ## The weights are normalized within each class to that class's prior,
      ## so that a loss and an edge weight the classes the same way.
      w = priorNormalize (w, gY, this.Prior);

      [~, s] = predict (this, X);
      l = classificationLoss (LossFun, s, gY, w, this.Cost);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationNaiveBayes} {@var{lp} =} logp (@var{obj}, @var{X})
    ##
    ## Log unconditional probability density of new data.
    ##
    ## @code{@var{lp} = logp (@var{obj}, @var{X})} returns one value per
    ## observation, the logarithm of its density under the fitted model taken
    ## over all the classes, each weighted by its prior.  A markedly low value
    ## marks an observation the model finds unlike anything it was trained on,
    ## whatever class it would be assigned to.
    ##
    ## @end deftypefn
    function lp = logp (this, X)

      if (nargin < 2)
        error ("ClassificationNaiveBayes.logp: too few input arguments.");
      endif
      if (isempty (X))
        error ("ClassificationNaiveBayes.logp: X is empty.");
      endif
      if (columns (this.X) != columns (X))
        error (strcat ("ClassificationNaiveBayes.logp: X must have the", ...
                       " same number of predictors as the trained model."));
      endif

      nCls = numel (this.Prior);
      L = nbLogLik (X, this.DistributionNames, ...
                    this.DistributionParameters, nCls, ...
                    this.CategoricalLevels);
      L = L + log (this.Prior);

      ## Sum the classes in logs: the largest term is factored out so that a
      ## density small enough to underflow still contributes its logarithm.
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
    ## @deftypefn  {ClassificationNaiveBayes} {@var{label} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationNaiveBayes} {[@var{label}, @var{score}, @var{cost}] =} resubPredict (@var{obj})
    ##
    ## Classify the training data with the trained model.
    ##
    ## The same as calling @code{predict} on the data the model was fitted on,
    ## with the rows that were dropped for missing values left out.
    ##
    ## @end deftypefn
    function [label, score, cost] = resubPredict (this)

      if (nargin < 1)
        error (strcat ("ClassificationNaiveBayes.resubPredict:", ...
                       " too few input arguments."));
      endif
      [label, score, cost] = predict (this, nbTrainX (this));

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationNaiveBayes} {@var{m} =} resubMargin (@var{obj})
    ##
    ## Classification margin on the training data.
    ##
    ## @end deftypefn
    function m = resubMargin (this)

      if (nargin < 1)
        error (strcat ("ClassificationNaiveBayes.resubMargin:", ...
                       " too few input arguments."));
      endif
      [X, Y] = nbTrainX (this);
      m = margin (this, X, Y);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationNaiveBayes} {@var{e} =} resubEdge (@var{obj})
    ##
    ## Classification edge on the training data.
    ##
    ## @end deftypefn
    function e = resubEdge (this)

      if (nargin < 1)
        error (strcat ("ClassificationNaiveBayes.resubEdge:", ...
                       " too few input arguments."));
      endif
      [X, Y] = nbTrainX (this);
      e = edge (this, X, Y);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNaiveBayes} {@var{l} =} resubLoss (@var{obj})
    ## @deftypefnx {ClassificationNaiveBayes} {@var{l} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss on the training data.
    ##
    ## Takes the same options as @code{loss}.
    ##
    ## @end deftypefn
    function l = resubLoss (this, varargin)

      if (nargin < 1)
        error (strcat ("ClassificationNaiveBayes.resubLoss:", ...
                       " too few input arguments."));
      endif
      [X, Y] = nbTrainX (this);
      l = loss (this, X, Y, varargin{:});

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationNaiveBayes} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationNaiveBayes object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## ClassificationNaiveBayes object into an Octave binary file, the name of
    ## which is specified in @var{filename}, along with an extra variable, which
    ## defines the type classification object these variables constitute.  Use
    ## @code{loadmodel} in order to load a classification object into Octave's
    ## workspace.
    ##
    ## @seealso{loadmodel, fitcnb, ClassificationNaiveBayes}
    ## @end deftypefn
    function savemodel (this, fname)

      if (nargin < 2)
        error (strcat ("ClassificationNaiveBayes.savemodel:", ...
                       " too few input arguments."));
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error (strcat ("ClassificationNaiveBayes.savemodel:", ...
                       " FNAME must be a character vector."));
      endif

      ## Generate variable for class name
      classdef_name = 'ClassificationNaiveBayes';

      ## Create variables from model properties
      X                      = this.X;
      Y                      = this.Y;
      RowsUsed               = this.RowsUsed;
      W                      = this.W;
      ModelParameters        = this.ModelParameters;
      NumObservations        = this.NumObservations;
      BinEdges               = this.BinEdges;
      PredictorNames         = this.PredictorNames;
      CategoricalPredictors  = this.CategoricalPredictors;
      ResponseName           = this.ResponseName;
      ExpandedPredictorNames = this.ExpandedPredictorNames;
      ClassNames             = this.ClassNames;
      Prior                  = this.Prior;
      Cost                   = this.Cost;
      ScoreTransform         = this.ScoreTransform;
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
      STfun                  = this.STfun;

      ## Save classdef name and all model properties as individual variables
      HyperparameterOptimizationResults = this.HyperparameterOptimizationResults;
      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'RowsUsed', 'W', ...
            'ModelParameters', 'NumObservations', 'BinEdges', ...
            'PredictorNames', 'CategoricalPredictors', 'ResponseName', ...
            'ExpandedPredictorNames', 'ClassNames', 'Prior', 'Cost', ...
            'ScoreTransform', 'DistributionNames', 'Mu', 'Sigma', ...
            'DistributionParameters', 'CategoricalLevels', 'Kernel', ...
            'Support', 'Width', 'STfun', ...
            'HyperparameterOptimizationResults');

    endfunction

  endmethods

  methods (Static, Hidden)

    function mdl = load_model (filename, data)

      ## The smallest fit the class accepts, filled property by property
      ## below.  Two observations per class give every distribution something
      ## to estimate from, and nothing of the stub survives the copy.
      mdl = ClassificationNaiveBayes ([0; 1; 2; 3], [1; 1; 2; 2]);

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
          msg = 'ClassificationNaiveBayes.load_model: invalid model in ''%s''.';
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
%!test  # MATLAB parity: the surface a default fit reports
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! assert_equal (class (Mdl), 'ClassificationNaiveBayes');
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (Mdl.ClassNames, unique (species));
%! assert_equal (Mdl.Prior, [1/3, 1/3, 1/3], 1e-15);
%! assert_equal (Mdl.Cost, [0, 1, 1; 1, 0, 1; 1, 1, 0]);
%! assert_equal (Mdl.ResponseName, 'Y');
%! assert_equal (Mdl.PredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (Mdl.ExpandedPredictorNames, {'x1', 'x2', 'x3', 'x4'});
%! assert_equal (Mdl.DistributionNames, {'normal', 'normal', 'normal', 'normal'});
%! assert_equal (Mdl.ScoreTransform, 'none');
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (Mdl.RowsUsed, []);

## Mu and Sigma are the standardization parameters and stay empty, where the
## class conditional ones live in DistributionParameters.  BinEdges is a cell
## rather than an empty matrix, as code reaching into it with cellfun needs.
%!test  # MATLAB parity: the properties that are deliberately empty
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! assert_equal (Mdl.Mu, []);
%! assert_equal (Mdl.Sigma, []);
%! assert_equal (class (Mdl.BinEdges), 'cell');
%! assert_equal (Mdl.BinEdges, {});
%! assert_equal (Mdl.Width, []);
%! assert_equal (Mdl.CategoricalLevels, cell (1, 4));

%!test  # MATLAB parity: the fitted normal parameters and the weights
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! assert_equal (size (Mdl.DistributionParameters), [3, 4]);
%! assert_equal (Mdl.DistributionParameters{1,1}, ...
%!               [5.005999999999998; 0.352489687213451], 1e-13);
%! assert_equal (Mdl.DistributionParameters{2,3}, ...
%!               [4.260000000000001; 0.469910977239958], 1e-13);
%! assert_equal (Mdl.W, repmat (1/150, 150, 1), 1e-15);
%! assert_equal (sum (Mdl.W), 1, 1e-14);

%!test  # MATLAB parity: predict returns the label, posterior and cost
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! [label, score, cost] = predict (Mdl, meas(1:3,:));
%! assert_equal (label, {'setosa'; 'setosa'; 'setosa'});
%! assert_equal (score, repmat ([1, 0, 0], 3, 1), 1e-12);
%! assert_equal (cost, repmat ([0, 1, 1], 3, 1), 1e-12);

%!test  # MATLAB parity: resubstitution loss, edge and margin
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! assert_equal (resubLoss (Mdl), 0.04, 1e-14);
%! assert_equal (resubEdge (Mdl), 0.894430597464877, 1e-12);
%! assert_equal (sum (resubMargin (Mdl)), 134.164589619731402, 1e-10);
%! assert_equal (resubMargin (Mdl)(1:5), ones (5, 1), 1e-12);

%!test  # MATLAB parity: resubPredict agrees with predict on the training data
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! [rl, rs, rc] = resubPredict (Mdl);
%! [pl, ps, pc] = predict (Mdl, meas);
%! assert_equal (rl, pl);
%! assert_equal (rs, ps);
%! assert_equal (rc, pc);

%!test  # MATLAB parity: every one of the eight loss functions
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! assert_equal (loss (Mdl, meas, species), 0.04, 1e-14);
%! assert_equal (loss (Mdl, meas, species, 'LossFun', 'mincost'), 0.04, 1e-14);
%! assert_equal (loss (Mdl, meas, species, 'LossFun', 'classiferror'), ...
%!               0.04, 1e-14);
%! assert_equal (loss (Mdl, meas, species, 'LossFun', 'classifcost'), ...
%!               0.04, 1e-14);
%! assert_equal (loss (Mdl, meas, species, 'LossFun', 'binodeviance'), ...
%!               0.149545751140664, 1e-12);
%! assert_equal (loss (Mdl, meas, species, 'LossFun', 'exponential'), ...
%!               0.395395632839880, 1e-12);
%! assert_equal (loss (Mdl, meas, species, 'LossFun', 'hinge'), ...
%!               0.052784701267562, 1e-12);
%! assert_equal (loss (Mdl, meas, species, 'LossFun', 'logit'), ...
%!               0.331045591275855, 1e-12);
%! assert_equal (loss (Mdl, meas, species, 'LossFun', 'quadratic'), ...
%!               0.033005648597277, 1e-12);

%!test  # MATLAB parity: weighted loss and edge
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! w = (1:150)' / sum (1:150);
%! assert_equal (loss (Mdl, meas, species, 'Weights', w), ...
%!               0.037013271417641, 1e-12);
%! assert_equal (edge (Mdl, meas, species, 'Weights', w), ...
%!               0.898902916457462, 1e-12);
%! assert_equal (edge (Mdl, meas, species), 0.894430597464877, 1e-12);

%!test  # MATLAB parity: logp over all the classes
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! lp = logp (Mdl, meas);
%! assert_equal (numel (lp), 150);
%! assert_equal (lp(1), 1.026591235856343, 1e-12);
%! assert_equal (lp(5), 0.977098877116533, 1e-12);
%! assert_equal (sum (lp), -309.559846128495394, 1e-10);

%!test  # MATLAB parity: a given prior reweights the observations and the loss
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'Prior', [0.2, 0.3, 0.5]);
%! assert_equal (Mdl.Prior, [0.2, 0.3, 0.5], 1e-15);
%! assert_equal (Mdl.W(1), 0.004, 1e-15);
%! assert_equal (Mdl.W(51), 0.006, 1e-15);
%! assert_equal (Mdl.W(101), 0.010, 1e-15);
%! assert_equal (resubLoss (Mdl), 0.054, 1e-14);
%! assert_equal (edge (Mdl, meas, species), 0.873545897578786, 1e-12);

%!test  # MATLAB parity: a uniform prior, and a given cost
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'Prior', 'uniform');
%! assert_equal (Mdl.Prior, [1/3, 1/3, 1/3], 1e-15);
%! Mdl = fitcnb (meas, species, 'Cost', [0, 1, 2; 1, 0, 1; 2, 1, 0]);
%! assert_equal (Mdl.Cost, [0, 1, 2; 1, 0, 1; 2, 1, 0]);

%!test  # MATLAB parity: ScoreTransform is applied to the posterior
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'ScoreTransform', 'logit');
%! assert_equal (Mdl.ScoreTransform, 'logit');
%! [~, score] = predict (Mdl, meas(1:2,:));
%! assert_equal (score, repmat ([0.731058578630005, 0.5, 0.5], 2, 1), 1e-12);

%!test  # MATLAB parity: the kernel densities and their default bandwidths
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'DistributionNames', 'kernel');
%! assert_equal (Mdl.DistributionNames, ...
%!               {'kernel', 'kernel', 'kernel', 'kernel'});
%! assert_equal (Mdl.Kernel, {'normal', 'normal', 'normal', 'normal'});
%! assert_equal (Mdl.Support, ...
%!               {'unbounded', 'unbounded', 'unbounded', 'unbounded'});
%! assert_equal (class (Mdl.DistributionParameters{1,1}), ...
%!               'prob.KernelDistribution');
%! width = [0.143628884694882, 0.179536105868602, 0.071814442347441, ...
%!          0.242194206816745; ...
%!          0.251350548216043, 0.143628884694882, 0.251350548216043, ...
%!          0.107721663521161; ...
%!          0.287257769389764, 0.143628884694882, 0.323164990563484, ...
%!          0.143628884694882];
%! assert_equal (Mdl.Width, width, 1e-12);

## The bandwidth belongs to the space the density is smoothed in, so a bounded
## support takes it in the transformed one.
%!test  # MATLAB parity: an explicit bandwidth, and a positive support
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'DistributionNames', 'kernel', ...
%!               'Kernel', 'triangle', 'Width', 0.5);
%! assert_equal (Mdl.Kernel, {'triangle', 'triangle', 'triangle', 'triangle'});
%! assert_equal (Mdl.Width, repmat (0.5, 3, 4), 1e-15);
%! x = [0.13; 0.41; 0.22; 1.87; 0.55; 0.09; 2.94; 0.31; 0.68; 0.17; ...
%!      0.44; 3.61; 0.26; 0.72; 0.05; 1.13; 0.38; 0.91; 0.19; 4.52];
%! y = [1.02; 1.44; 0.87; 1.19; 2.63; 1.07; 0.95; 1.31; 1.76; 1.12; ...
%!      0.99; 1.28; 3.41; 1.05; 1.21; 0.93; 1.38; 1.14; 1.09; 2.02];
%! G = [repmat({'a'}, 20, 1); repmat({'b'}, 20, 1)];
%! Mdl = fitcnb ([x; y], G, 'DistributionNames', 'kernel');
%! assert_equal (Mdl.Width, ...
%!               [0.237209723894721; 0.138012930266019], 1e-12);
%! Mdl = fitcnb ([x; y], G, 'DistributionNames', 'kernel', ...
%!               'Support', 'positive');
%! assert_equal (Mdl.Width, ...
%!               [0.675582146901479; 0.127329555528534], 1e-12);

%!test  # MATLAB parity: a kernel model classifies as MATLAB does
%! x = [0.13; 0.41; 0.22; 1.87; 0.55; 0.09; 2.94; 0.31; 0.68; 0.17; ...
%!      0.44; 3.61; 0.26; 0.72; 0.05; 1.13; 0.38; 0.91; 0.19; 4.52];
%! y = [1.02; 1.44; 0.87; 1.19; 2.63; 1.07; 0.95; 1.31; 1.76; 1.12; ...
%!      0.99; 1.28; 3.41; 1.05; 1.21; 0.93; 1.38; 1.14; 1.09; 2.02];
%! G = [repmat({'a'}, 20, 1); repmat({'b'}, 20, 1)];
%! Mdl = fitcnb ([x; y], G, 'DistributionNames', 'kernel');
%! [label, score] = predict (Mdl, [0.5; 1.5; 3.0]);
%! assert_equal (label, {'a'; 'b'; 'a'});
%! assert_equal (score, [0.991485542298633, 0.008514457701367; ...
%!                       0.121739189175926, 0.878260810824074; ...
%!                       0.936549935287367, 0.063450064712633], 1e-12);
%! assert_equal (resubLoss (Mdl), 0.075, 1e-14);

%!test  # MATLAB parity: one distribution per predictor
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'DistributionNames', ...
%!               {'normal', 'kernel', 'normal', 'kernel'});
%! assert_equal (Mdl.DistributionNames, ...
%!               {'normal', 'kernel', 'normal', 'kernel'});
%! assert_equal (class (Mdl.DistributionParameters{1,1}), 'double');
%! assert_equal (class (Mdl.DistributionParameters{1,2}), ...
%!               'prob.KernelDistribution');
%! assert_equal (Mdl.DistributionParameters{1,2}.Bandwidth, ...
%!               0.179536105868602, 1e-12);

## ModelParameters records the arguments as they were given, where the
## properties of the same name record what they were resolved to.
%!test  # MATLAB parity: ModelParameters
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! assert_equal (fieldnames (Mdl.ModelParameters), ...
%!               {'DistributionNames'; 'Kernel'; 'Support'; 'Width'; ...
%!                'StandardizeData'; 'Version'; 'Method'; 'Type'});
%! assert_equal (Mdl.ModelParameters.DistributionNames, 'normal');
%! assert_equal (Mdl.ModelParameters.Kernel, []);
%! assert_equal (Mdl.ModelParameters.Method, 'NaiveBayes');
%! assert_equal (Mdl.ModelParameters.Type, 'classification');
%! assert_equal (Mdl.ModelParameters.Version, 1);
%! Mdl = fitcnb (meas, species, 'DistributionNames', 'kernel');
%! assert_equal (Mdl.ModelParameters.DistributionNames, 'kernel');
%! assert_equal (Mdl.ModelParameters.Kernel, 'normal');
%! assert_equal (Mdl.ModelParameters.Support, 'unbounded');
%! assert_equal (isnan (Mdl.ModelParameters.Width), true);
%! assert_equal (Mdl.ModelParameters.StandardizeData, 0);

%!test  # naming a subset of the classes keeps only those observations
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'ClassNames', {'setosa', 'virginica'});
%! assert_equal (Mdl.ClassNames, {'setosa'; 'virginica'});
%! assert_equal (Mdl.NumObservations, 100);
%! assert_equal (sum (Mdl.RowsUsed), 100);
%! assert_equal (size (Mdl.DistributionParameters), [2, 4]);
%! assert_equal (Mdl.Prior, [0.5, 0.5], 1e-15);

%!test  # a row holding a missing value is dropped, and RowsUsed says so
%! load fisheriris
%! X = meas;
%! X(3,2) = NaN;
%! X(77,4) = NaN;
%! Mdl = fitcnb (X, species);
%! assert_equal (Mdl.NumObservations, 148);
%! assert_equal (sum (Mdl.RowsUsed), 148);
%! assert_equal (Mdl.RowsUsed([3, 77]), [false; false]);
%! assert_equal (size (Mdl.X), [150, 4]);

%!test  # PredictorNames and ResponseName are carried through
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'PredictorNames', {'a', 'b', 'c', 'd'}, ...
%!               'ResponseName', 'flower');
%! assert_equal (Mdl.PredictorNames, {'a', 'b', 'c', 'd'});
%! assert_equal (Mdl.ExpandedPredictorNames, {'a', 'b', 'c', 'd'});
%! assert_equal (Mdl.ResponseName, 'flower');

%!test  # Cost and Prior may be assigned after fitting, and W follows the prior
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! Mdl.Cost = [0, 2, 2; 2, 0, 2; 2, 2, 0];
%! assert_equal (Mdl.Cost, [0, 2, 2; 2, 0, 2; 2, 2, 0]);
%! Mdl.Prior = [0.2, 0.3, 0.5];
%! assert_equal (Mdl.Prior, [0.2, 0.3, 0.5], 1e-15);
%! assert_equal (Mdl.W(1), 0.004, 1e-15);
%! assert_equal (Mdl.W(101), 0.010, 1e-15);

%!test  # a numeric response is classified as a numeric response
%! X = [1, 2; 1.1, 2.1; 5, 6; 5.2, 6.1; 0.9, 1.8; 5.1, 6.2];
%! Y = [1; 1; 2; 2; 1; 2];
%! Mdl = fitcnb (X, Y);
%! assert_equal (Mdl.ClassNames, [1; 2]);
%! assert_equal (predict (Mdl, [1, 2; 5, 6]), [1; 2]);

%!test  # MATLAB parity: compact classifies identically to the model it came from
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! CMdl = compact (Mdl);
%! assert_equal (class (CMdl), 'CompactClassificationNaiveBayes');
%! [ml, ms] = predict (Mdl, meas);
%! [cl, cs] = predict (CMdl, meas);
%! assert_equal (cl, ml);
%! assert_equal (cs, ms);
%! assert_equal (loss (CMdl, meas, species), 0.04, 1e-14);
%! assert_equal (edge (CMdl, meas, species), 0.894430597464877, 1e-12);

%!test  # MATLAB parity: leave-one-out cross-validation
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! CVMdl = crossval (Mdl, 'Leaveout', 'on');
%! assert_equal (class (CVMdl), 'ClassificationPartitionedModel');
%! assert_equal (CVMdl.KFold, 150);
%! assert_equal (class (CVMdl.Trained{1}), 'CompactClassificationNaiveBayes');
%! assert_equal (kfoldLoss (CVMdl), 0.046666666666667, 1e-12);
%! assert_equal (sum (strcmp (kfoldPredict (CVMdl), species)), 143);

%!test  # a k-fold partition holds k folds, each a model of its own
%! load fisheriris
%! CVMdl = crossval (fitcnb (meas, species), 'KFold', 3);
%! assert_equal (CVMdl.KFold, 3);
%! assert_equal (numel (CVMdl.Trained), 3);
%! assert_equal (class (CVMdl.Trained{3}), 'CompactClassificationNaiveBayes');
%! assert_equal (CVMdl.CrossValidatedModel, 'NaiveBayes');

## The categorical distributions.  Every level's count is raised by one before
## it is normalized, so a level a class never took stays possible.
%!test  # MATLAB parity: mvmn fits a smoothed distribution over the levels
%! X = [1, 2; 1, 3; 2, 2; 2, 3; 1, 2; 3, 1; 3, 3; 2, 1; 1, 1; 3, 2; ...
%!      2, 2; 1, 3; 3, 1; 2, 3; 1, 1; 3, 2; 2, 1; 1, 2; 3, 3; 2, 2];
%! Y = [repmat({'a'}, 10, 1); repmat({'b'}, 10, 1)];
%! Mdl = fitcnb (X, Y, 'DistributionNames', 'mvmn', ...
%!               'CategoricalPredictors', [1, 2]);
%! assert_equal (Mdl.DistributionNames, {'mvmn', 'mvmn'});
%! assert_equal (Mdl.CategoricalPredictors, [1, 2]);
%! assert_equal (Mdl.CategoricalLevels{1}, [1; 2; 3]);
%! assert_equal (Mdl.CategoricalLevels{2}, [1; 2; 3]);
%! assert_equal (Mdl.DistributionParameters{1,1}, ...
%!               [0.384615384615385; 0.307692307692308; 0.307692307692308], ...
%!               1e-12);
%! assert_equal (Mdl.DistributionParameters{2,2}(2), ...
%!               0.384615384615385, 1e-12);
%! assert_equal (all (isnan (Mdl.Width(:))), true);
%! assert_equal (size (Mdl.Width), [2, 2]);

%!test  # MATLAB parity: an mvmn model classifies as MATLAB does
%! X = [1, 2; 1, 3; 2, 2; 2, 3; 1, 2; 3, 1; 3, 3; 2, 1; 1, 1; 3, 2; ...
%!      2, 2; 1, 3; 3, 1; 2, 3; 1, 1; 3, 2; 2, 1; 1, 2; 3, 3; 2, 2];
%! Y = [repmat({'a'}, 10, 1); repmat({'b'}, 10, 1)];
%! Mdl = fitcnb (X, Y, 'DistributionNames', 'mvmn', ...
%!               'CategoricalPredictors', [1, 2]);
%! [label, score] = predict (Mdl, [1, 1; 3, 3; 2, 2]);
%! assert_equal (label, {'a'; 'a'; 'b'});
%! assert_equal (score(1,1), 0.555555555555556, 1e-12);
%! assert_equal (score(3,2), 0.555555555555556, 1e-12);
%! assert_equal (resubLoss (Mdl), 0.45, 1e-14);
%! assert_equal (logp (Mdl, [1, 1; 3, 3]), ...
%!               [-2.239526957026909; -2.357309992683292], 1e-12);

%!test  # MATLAB parity: a level a class never took keeps a probability
%! Z = [1, 2; 1, 3; 1, 2; 1, 3; 1, 2; 2, 1; 2, 3; 2, 1; 2, 3; 2, 2];
%! G = [repmat({'p'}, 5, 1); repmat({'q'}, 5, 1)];
%! Mdl = fitcnb (Z, G, 'DistributionNames', 'mvmn', ...
%!               'CategoricalPredictors', [1, 2]);
%! assert_equal (Mdl.CategoricalLevels{1}, [1; 2]);
%! assert_equal (Mdl.DistributionParameters{1,1}, ...
%!               [0.857142857142857; 0.142857142857143], 1e-12);
%! assert_equal (Mdl.DistributionParameters{2,1}, ...
%!               [0.142857142857143; 0.857142857142857], 1e-12);
%! [label, score] = predict (Mdl, [2, 2; 1, 1]);
%! assert_equal (label, {'q'; 'p'});
%! assert_equal (score(1,:), [0.25, 0.75], 1e-12);

%!test  # MATLAB parity: naming a predictor categorical makes it mvmn
%! X = [1, 2; 1, 3; 2, 2; 2, 3; 1, 2; 3, 1; 3, 3; 2, 1; 1, 1; 3, 2; ...
%!      2, 2; 1, 3; 3, 1; 2, 3; 1, 1; 3, 2; 2, 1; 1, 2; 3, 3; 2, 2];
%! Y = [repmat({'a'}, 10, 1); repmat({'b'}, 10, 1)];
%! Mdl = fitcnb (X, Y, 'CategoricalPredictors', [1, 2]);
%! assert_equal (Mdl.DistributionNames, {'mvmn', 'mvmn'});
%! assert_equal (Mdl.CategoricalPredictors, [1, 2]);
%! assert_equal (Mdl.DistributionParameters{1,1}(1), ...
%!               0.384615384615385, 1e-12);

%!test  # MATLAB parity: a normal and a categorical predictor side by side
%! X = [1, 2; 1, 3; 2, 2; 2, 3; 1, 2; 3, 1; 3, 3; 2, 1; 1, 1; 3, 2; ...
%!      2, 2; 1, 3; 3, 1; 2, 3; 1, 1; 3, 2; 2, 1; 1, 2; 3, 3; 2, 2];
%! Y = [repmat({'a'}, 10, 1); repmat({'b'}, 10, 1)];
%! Mdl = fitcnb (X, Y, 'DistributionNames', {'normal', 'mvmn'}, ...
%!               'CategoricalPredictors', 2);
%! assert_equal (Mdl.DistributionNames, {'normal', 'mvmn'});
%! assert_equal (isempty (Mdl.CategoricalLevels{1}), true);
%! assert_equal (Mdl.CategoricalLevels{2}, [1; 2; 3]);
%! assert_equal (Mdl.DistributionParameters{1,1}, ...
%!               [1.900000000000001; 0.875595035770913], 1e-12);
%! assert_equal (Mdl.DistributionParameters{1,2}(2), ...
%!               0.384615384615385, 1e-12);

## The multinomial reads a row as token counts and is one distribution over
## the whole predictor vector, so its DistributionNames stays a character
## vector where every other distribution reports one name per predictor.
%!test  # MATLAB parity: mn over token counts
%! C = [2, 0, 1; 1, 3, 0; 0, 1, 4; 3, 1, 0; ...
%!      0, 2, 2; 1, 0, 3; 4, 1, 1; 0, 3, 2];
%! L = [repmat({'x'}, 4, 1); repmat({'y'}, 4, 1)];
%! Mdl = fitcnb (C, L, 'DistributionNames', 'mn');
%! assert_equal (Mdl.DistributionNames, 'mn');
%! assert_equal (class (Mdl.DistributionNames), 'char');
%! assert_equal (size (Mdl.DistributionParameters), [2, 3]);
%! assert_equal (Mdl.DistributionParameters{1,1}, 0.368421052631579, 1e-12);
%! assert_equal (Mdl.DistributionParameters{1,2}, 0.315789473684211, 1e-12);
%! assert_equal (Mdl.DistributionParameters{2,3}, 0.409090909090909, 1e-12);
%! assert_equal (Mdl.CategoricalPredictors, []);

%!test  # MATLAB parity: an mn model classifies as MATLAB does
%! C = [2, 0, 1; 1, 3, 0; 0, 1, 4; 3, 1, 0; ...
%!      0, 2, 2; 1, 0, 3; 4, 1, 1; 0, 3, 2];
%! L = [repmat({'x'}, 4, 1); repmat({'y'}, 4, 1)];
%! Mdl = fitcnb (C, L, 'DistributionNames', 'mn');
%! [label, score] = predict (Mdl, [1, 1, 1; 4, 0, 0]);
%! assert_equal (label, {'x'; 'x'});
%! assert_equal (score(1,1), 0.508585484679865, 1e-12);
%! assert_equal (score(2,1), 0.769060987976952, 1e-12);
%! assert_equal (resubLoss (Mdl), 0.25, 1e-14);

## An observation carrying a level the model never saw tells it nothing, so
## its posterior is the prior.  One such predictor decides the whole
## observation, whatever the others say.
%!test  # MATLAB parity: an unseen level falls back to the prior
%! Z = [1, 1; 1, 1; 1, 2; 2, 2; 2, 2; 2, 1; 1, 1; 2, 2];
%! G = [repmat({'p'}, 4, 1); repmat({'q'}, 4, 1)];
%! Mdl = fitcnb (Z, G, 'DistributionNames', 'mvmn', ...
%!               'CategoricalPredictors', [1, 2]);
%! [~, score] = predict (Mdl, [1, 1]);
%! assert_equal (score, [0.666666666666667, 0.333333333333333], 1e-12);
%! [~, score, cost] = predict (Mdl, [3, 1]);
%! assert_equal (score, [0.5, 0.5], 1e-12);
%! assert_equal (cost, [0.5, 0.5], 1e-12);
%! [~, score] = predict (Mdl, [3, 3]);
%! assert_equal (score, [0.5, 0.5], 1e-12);
%! assert_equal (logp (Mdl, [1, 1]), -1.386294361119891, 1e-12);
%! assert_equal (logp (Mdl, [3, 1]), -Inf);

%!test  # MATLAB parity: the fallback is the prior, not a uniform distribution
%! Z = [1, 1; 1, 1; 1, 2; 1, 2; 1, 1; 1, 2; 2, 2; 2, 1];
%! G = [repmat({'p'}, 6, 1); repmat({'q'}, 2, 1)];
%! Mdl = fitcnb (Z, G, 'DistributionNames', 'mvmn', ...
%!               'CategoricalPredictors', [1, 2]);
%! assert_equal (Mdl.Prior, [0.75, 0.25], 1e-15);
%! [~, score] = predict (Mdl, [3, 1]);
%! assert_equal (score, [0.75, 0.25], 1e-12);
%! Mdl = fitcnb (Z, G, 'DistributionNames', 'mvmn', ...
%!               'CategoricalPredictors', [1, 2], 'Prior', [0.2, 0.8]);
%! [label, score] = predict (Mdl, [3, 1]);
%! assert_equal (score, [0.2, 0.8], 1e-12);
%! assert_equal (label, {'q'});

%!test  # MATLAB parity: an informative predictor cannot rescue the row
%! Z = [1, 1; 1, 1; 1, 2; 2, 2; 2, 2; 2, 1; 1, 1; 2, 2];
%! G = [repmat({'p'}, 4, 1); repmat({'q'}, 4, 1)];
%! Mdl = fitcnb (Z, G, 'DistributionNames', {'normal', 'mvmn'}, ...
%!               'CategoricalPredictors', 2);
%! [~, score] = predict (Mdl, [1, 3]);
%! assert_equal (score, [0.5, 0.5], 1e-12);
%! assert_equal (logp (Mdl, [1, 3]), -Inf);

## Naming 'mvmn' for a predictor makes it categorical, which MATLAB reports
## rather than doing silently.
%!warning<the 'mvmn' distribution was named for a predictor that is not in 'CategoricalPredictors'> ...
%! fitcnb ([1, 2; 2, 1; 1, 1; 2, 2], [1; 1; 2; 2], ...
%!         'DistributionNames', 'mvmn');

%!test  # the warning is not raised when the predictor is already categorical
%! Mdl = fitcnb ([1, 2; 2, 1; 1, 1; 2, 2], [1; 1; 2; 2], ...
%!               'DistributionNames', 'mvmn', 'CategoricalPredictors', 'all');
%! assert_equal (Mdl.CategoricalPredictors, [1, 2]);
%! assert_equal (Mdl.DistributionNames, {'mvmn', 'mvmn'});

## A fitted model survives savemodel and loadmodel: every property comes back
## as it was and it predicts the same.
%!test
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'ClassificationNaiveBayes');
%! p = properties (Mdl);
%! for i = 1:numel (p)
%!   assert_equal (M2.(p{i}), Mdl.(p{i}));
%! endfor
%! assert_equal (predict (M2, meas(1:10,:)), predict (Mdl, meas(1:10,:)));

## A kernel predictor stores a classdef object, which Octave's save cannot
## serialize.  It is written as the sample and refitted at the recorded
## bandwidth, so the density comes back identical and not merely close.
%!test
%! load fisheriris
%! Mdl = fitcnb (meas, species, 'DistributionNames', 'kernel');
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2.DistributionParameters{1,1}), ...
%!               'prob.KernelDistribution');
%! assert_equal (M2.Width, Mdl.Width);
%! assert_equal (predict (M2, meas(1:10,:)), predict (Mdl, meas(1:10,:)));
%! assert_equal (logp (M2, meas(1:10,:)), logp (Mdl, meas(1:10,:)), 1e-12);

## Mixed distributions, a categorical predictor among them, and a model
## carrying a transform and a non-default cost.
%!test
%! load fisheriris
%! X = [meas, round(meas(:,1))];
%! Mdl = fitcnb (X, species, 'DistributionNames', ...
%!               {'normal', 'kernel', 'normal', 'kernel', 'mvmn'}, ...
%!               'CategoricalPredictors', 5, ...
%!               'Cost', [0, 2, 1; 1, 0, 1; 1, 1, 0]);
%! Mdl.ScoreTransform = 'logit';
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (M2.Cost, Mdl.Cost);
%! assert_equal (M2.ScoreTransform, 'logit');
%! assert_equal (M2.CategoricalLevels, Mdl.CategoricalLevels);
%! [l1, s1] = predict (Mdl, X(1:10,:));
%! [l2, s2] = predict (M2, X(1:10,:));
%! assert_equal (l2, l1);
%! assert_equal (s2, s1);

%!error<ClassificationNaiveBayes.savemodel: too few input arguments.> ...
%! savemodel (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]))
%!error<ClassificationNaiveBayes.savemodel: FNAME must be a character vector.> ...
%! savemodel (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), 1)
%!error<ClassificationNaiveBayes.savemodel: FNAME must be a character vector.> ...
%! savemodel (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), ['ab'; 'cd'])

%!error<ClassificationNaiveBayes: the 'mn' distribution applies to every predictor at once and cannot be named for some of them.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'DistributionNames', {'mn', 'normal'})
%!error<ClassificationNaiveBayes: a character vector 'CategoricalPredictors' must be 'all'.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'CategoricalPredictors', 'some')
%!error<ClassificationNaiveBayes: a logical 'CategoricalPredictors' must have one element per predictor.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'CategoricalPredictors', [true, true, true])
%!error<ClassificationNaiveBayes: 'CategoricalPredictors' must be column indices into X.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'CategoricalPredictors', 5)

## A missing predictor is skipped, not propagated: the predictors are
## conditionally independent given the class, so the others still describe the
## observation.  An observation missing every predictor falls back to the
## prior, by the same arithmetic that gives an unseen level the prior.
%!test  # MATLAB parity: predict skips a NaN predictor
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! [label, score, cost] = predict (Mdl, [NaN, 3, 1, 0.2]);
%! assert_equal (label, {'setosa'});
%! assert_equal (score, [1, 0, 0], 1e-12);
%! assert_equal (cost, [0, 1, 1], 1e-12);
%! [~, s1] = predict (Mdl, [5.1, 3.5, 1.4, 0.2]);
%! [~, s2] = predict (Mdl, [5.1, NaN, 1.4, 0.2]);
%! assert_equal (s2, s1, 1e-12);
%! [~, s3] = predict (Mdl, [NaN, NaN, NaN, NaN]);
%! assert_equal (s3, Mdl.Prior, 1e-12);

## logp reports NaN where predict reports a class: predict can classify on the
## predictors it has, but a density is not defined for a partial observation.
%!test  # MATLAB parity: logp of an incomplete observation is NaN
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! assert_equal (isnan (logp (Mdl, [NaN, 3, 1, 0.2])), true);
%! assert_equal (isnan (logp (Mdl, [5.1, NaN, 1.4, 0.2])), true);
%! assert_equal (logp (Mdl, [5.1, 3.5, 1.4, 0.2]), 1.026591235856343, 1e-12);

%!test  # MATLAB parity: a row with a missing predictor still counts in a loss
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! X = [meas; NaN, 3, 1, 0.2];
%! Y = [species; {'setosa'}];
%! assert_equal (loss (Mdl, X, Y), 0.04, 1e-14);
%! assert_equal (loss (Mdl, X, Y, 'LossFun', 'classiferror'), 0.04, 1e-14);

## A predictor that does not vary within a class has no normal density to fit.
## MATLAB refuses the combination rather than answering from a distribution it
## never fitted, and so do we.  The refusal is per class and predictor, not per
## model: a kernel on the same column fits, which is the caller's way out.
%!test  # MATLAB parity: a kernel fits data a normal cannot
%! X = [1, 2; 2, 3; 3, 4; 10, 20; 10, 25];
%! Y = [1; 1; 1; 2; 2];
%! Mdl = fitcnb (X, Y, 'DistributionNames', 'kernel');
%! assert_equal (class (Mdl), 'ClassificationNaiveBayes');
%! assert_equal (Mdl.Width(2,1), 1, 1e-12);
%! Mdl = fitcnb (X, Y, 'DistributionNames', 'mvmn', ...
%!               'CategoricalPredictors', [1, 2]);
%! assert_equal (Mdl.CategoricalLevels{1}, [1; 2; 3; 10]);

## The row that proves the guard is per combination rather than per model:
## give the degenerate column a kernel and leave the other normal, and the fit
## goes through.  Widening the guard to the whole model would fail this.
%!test  # MATLAB parity: only the degenerate combination is refused
%! X = [1, 2; 2, 3; 3, 4; 10, 20; 10, 25];
%! Y = [1; 1; 1; 2; 2];
%! Mdl = fitcnb (X, Y, 'DistributionNames', {'kernel', 'normal'});
%! assert_equal (Mdl.DistributionNames, {'kernel', 'normal'});
%! assert_equal (Mdl.DistributionParameters{2,2}, [22.5; 3.535533905932738], ...
%!               1e-12);

## Test input validation
%!error<ClassificationNaiveBayes: too few input arguments.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5])
%!error<ClassificationNaiveBayes: Name-Value arguments must be in pairs.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], ones (4, 1), 'Prior')
%!error<ClassificationNaiveBayes: X must be a non-empty real numeric matrix.> ...
%! ClassificationNaiveBayes ('a', ones (4, 1))
%!error<ClassificationNaiveBayes: number of rows in X and Y must be equal.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], ones (3, 1))
%!error<ClassificationNaiveBayes: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'PredictorNames', 5)
%!error<ClassificationNaiveBayes: 'PredictorNames' must equal the number of columns in X.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'PredictorNames', {'a'})
%!error<ClassificationNaiveBayes: 'ResponseName' must be a character vector.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'ResponseName', 5)
%!error<ClassificationNaiveBayes: not all 'ClassNames' are present in Y.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'ClassNames', [1, 3])
%!error<ClassificationNaiveBayes: invalid parameter name 'nope'.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], 'nope', 1)
%!error<ClassificationNaiveBayes: unsupported distribution 'poisson'.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], ...
%!                           'DistributionNames', 'poisson')
%!error<ClassificationNaiveBayes: 'DistributionNames' must name one distribution per predictor.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], ...
%!                           'DistributionNames', {'normal'})
%!error<ClassificationNaiveBayes: unsupported kernel 'cosine'.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], ...
%!                           'DistributionNames', 'kernel', 'Kernel', 'cosine')
%!error<ClassificationNaiveBayes: 'Width' must be positive and real.> ...
%! ClassificationNaiveBayes ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], ...
%!                           'DistributionNames', 'kernel', 'Width', -1)
%!error<ClassificationNaiveBayes.predict: too few input arguments.> ...
%! predict (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]))
%!error<ClassificationNaiveBayes.predict: XC is empty.> ...
%! predict (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [])
%!error<ClassificationNaiveBayes.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), ones (2, 3))
%!error<ClassificationNaiveBayes.loss: too few input arguments.> ...
%! loss (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [1, 2; 2, 3; 3, 4; 4, 5])
%!error<ClassificationNaiveBayes.loss: invalid loss function.> ...
%! loss (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], ...
%!       'LossFun', 'nope')
%!error<ClassificationNaiveBayes.loss: 'Weights' must have one element per observation.> ...
%! loss (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2], ...
%!       'Weights', [1, 2])
%!error<ClassificationNaiveBayes.margin: too few input arguments.> ...
%! margin (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [1, 2; 2, 3; 3, 4; 4, 5])
%!error<ClassificationNaiveBayes.edge: too few input arguments.> ...
%! edge (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [1, 2; 2, 3; 3, 4; 4, 5])
%!error<ClassificationNaiveBayes.logp: X is empty.> ...
%! logp (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), [])
%!error<ClassificationNaiveBayes.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), 'KFold', 1)
%!error<ClassificationNaiveBayes.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), 'Leaveout', 5)
%!error<ClassificationNaiveBayes.crossval: specify only one of the optional Name-Value paired arguments.> ...
%! crossval (fitcnb ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), 'KFold', 2, 'Leaveout', 'on')
%!error<ClassificationNaiveBayes: a normal distribution cannot be fit for the combination of class 2 and predictor x1. The data has zero variance.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 10, 20], [1; 1; 1; 2])
%!error<ClassificationNaiveBayes: a normal distribution cannot be fit for the combination of class 2 and predictor x1. The data has zero variance.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 10, 20; 10, 25], [1; 1; 1; 2; 2])
%!error<ClassificationNaiveBayes: a normal distribution cannot be fit for the combination of class beta and predictor x1. The data has zero variance.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 10, 20; 10, 25], ...
%!         [{'alpha'}; {'alpha'}; {'alpha'}; {'beta'}; {'beta'}])
%!error<ClassificationNaiveBayes: a normal distribution cannot be fit for the combination of class 2 and predictor height. The data has zero variance.> ...
%! fitcnb ([1, 2; 2, 3; 3, 4; 10, 20; 10, 25], [1; 1; 1; 2; 2], ...
%!         'PredictorNames', {'height', 'weight'})
%!error<ClassificationNaiveBayes: a normal distribution cannot be fit for the combination of class 1 and predictor x1. The data has zero variance.> ...
%! fitcnb ([5, 2; 5, 3; 5, 4; 10, 20; 11, 25], [1; 1; 1; 2; 2])

## HyperparameterOptimizationResults is declared for MATLAB compatibility and
## stays empty, this class running no search over its hyperparameters.
%!test
%! load fisheriris
%! Mdl = fitcnb (meas, species);
%! assert_equal (isempty (Mdl.HyperparameterOptimizationResults), true);
