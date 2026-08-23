## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
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

classdef ClassificationDiscriminant
  ## -*- texinfo -*-
  ## @deftp {statistics} ClassificationDiscriminant
  ##
  ## Discriminant analysis classification
  ##
  ## The @code{ClassificationDiscriminant} class implements a
  ## discriminant analysis classifier object, which can predict responses for
  ## new data using the @code{predict} method.
  ##
  ## Discriminant analysis classification is a statistical method used to
  ## classify observations into predefined groups based on their
  ## characteristics.  It estimates the parameters of different distributions
  ## for each class and predicts the class of new observations by finding the
  ## one with the smallest misclassification cost.
  ##
  ## Create a @code{ClassificationDiscriminant} object by using the
  ## @code{fitcdiscr} function or the class constructor.
  ##
  ##
  ## Six discriminant types are available, in two families.  The linear family,
  ## @qcode{'linear'}, @qcode{'diagLinear'} and @qcode{'pseudoLinear'}, pools
  ## one covariance across the classes and separates them with a hyperplane.
  ## The quadratic family, @qcode{'quadratic'}, @qcode{'diagQuadratic'} and
  ## @qcode{'pseudoQuadratic'}, estimates a covariance per class and separates
  ## them with a quadric.  A @qcode{'diag'} type keeps only the variances,
  ## which is the same model as a @qcode{Gamma} of 1, and a @qcode{'pseudo'}
  ## type inverts a singular covariance rather than refusing it.
  ##
  ## @qcode{DiscrimType} may be assigned after fitting, but @emph{only within
  ## its own family}: the family is fixed when the model is fitted, because it
  ## decides which covariances the fit has to estimate.  Assigning it, or
  ## @qcode{Gamma}, re-derives @qcode{Sigma}, @qcode{LogDetSigma} and
  ## @qcode{Coeffs} without refitting.
  ##
  ## @seealso{fitcdiscr}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} W
    ##
    ## Observation weights
    ##
    ## A numeric column vector with one entry per observation used for fitting.
    ## Every observation carries the same weight, so the vector sums
    ## to one.  This property is read-only.
    ##
    ## @end deftp
    W               = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} X
    ##
    ## Predictor data
    ##
    ## A numeric matrix containing the unstandardized predictor data.  Each
    ## column of @var{X} represents one predictor (variable), and each row
    ## represents one observation.  This property is read-only.
    ##
    ## @end deftp
    X = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Y
    ##
    ## Class labels
    ##
    ## Specified as a logical or numeric column vector, or as a character array
    ## or a cell array of character vectors with the same number of rows as the
    ## predictor data.  Each row in @var{Y} is the observed class label for
    ## the corresponding row in @var{X}.  This property is read-only.
    ##
    ## @end deftp
    Y = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## A positive integer value specifying the number of observations in the
    ## training dataset used for training the ClassificationDiscriminant model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} RowsUsed
    ##
    ## Rows used for fitting
    ##
    ## A logical column vector with the same length as the observations in the
    ## original predictor data @var{X}, true for each row that was used for
    ## fitting the ClassificationDiscriminant model.  It is empty, @qcode{[]},
    ## when every observation was used, so a non-empty value means that rows
    ## holding missing values were dropped.  This property is read-only.
    ##
    ## @end deftp
    RowsUsed        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer value specifying the number of predictors in the
    ## training dataset used for training the ClassificationDiscriminant model.
    ## This property is read-only.
    ##
    ## @end deftp
    NumPredictors   = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which they appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames  = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} BetweenSigma
    ##
    ## Between-class covariance matrix
    ##
    ## A @math{P}-by-@math{P} matrix holding the covariance of the class means
    ## about the overall mean, weighted by how many observations each class
    ## contributes.  With @math{n_k} observations in class @math{k},
    ## @math{p_k = n_k / n} and @math{\bar{\mu} = \sum_k p_k \mu_k}, it is
    ##
    ## @example
    ## @group
    ## BetweenSigma = sum_k n_k (Mu(k,:) - mubar)' * (Mu(k,:) - mubar)
    ##                / (n * (1 - sum_k p_k^2))
    ## @end group
    ## @end example
    ##
    ## The denominator is the unbiased one for a weighted covariance, so a
    ## balanced fit divides by @math{n (K-1) / K}.  It reads the @strong{class
    ## sizes}, not @qcode{Prior}: assigning a prior leaves it where it was.  It
    ## is estimated for every discriminant type, the quadratic family included,
    ## since it describes the classes rather than the fit.  This property is
    ## read-only.
    ##
    ## @end deftp
    BetweenSigma    = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} CategoricalPredictors
    ##
    ## Indices of the categorical predictors
    ##
    ## A numeric vector of column indices into @code{X} naming the predictors
    ## treated as categorical, and empty when none is.  This property is
    ## read-only.
    ##
    ## @end deftp
    CategoricalPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} ExpandedPredictorNames
    ##
    ## Names of the predictors as the model expanded them
    ##
    ## A cell array of character vectors.  It matches @code{PredictorNames}
    ## unless a categorical predictor was expanded into indicator variables.
    ## This property is read-only.
    ##
    ## @end deftp
    ExpandedPredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName    = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} ClassNames
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
    ClassNames      = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Sigma
    ##
    ## Within-class covariance
    ##
    ## A numeric array whose shape follows @qcode{DiscrimType}, with @math{P}
    ## predictors and @math{K} classes:
    ##
    ## @multitable @columnfractions 0.4 0.25 0.35
    ## @headitem @var{DiscrimType} @tab @var{Sigma} @tab @var{LogDetSigma}
    ## @item @qcode{'linear'}, @qcode{'pseudoLinear'} @tab @math{PxP}
    ## @tab scalar
    ## @item @qcode{'quadratic'}, @qcode{'pseudoQuadratic'} @tab @math{PxPxK}
    ## @tab @math{Kx1}
    ## @item @qcode{'diagLinear'} @tab @math{1xP} @tab scalar
    ## @item @qcode{'diagQuadratic'} @tab @math{1xPxK} @tab @math{Kx1}
    ## @end multitable
    ##
    ## The linear family pools one covariance across the classes and the
    ## quadratic family estimates one per class.  This property is read-only,
    ## but it is re-derived whenever @qcode{DiscrimType} or @qcode{Gamma} is
    ## assigned.
    ##
    ## @end deftp
    Sigma           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Mu
    ##
    ## Class means
    ##
    ## A @math{K*P} numeric matrix specifying the mean of the multivariate
    ## normal distribution of each corresponding class, where @math{K} is the
    ## number of classes and @math{P} is the number of predictors in @var{X}.
    ## This property is read-only.
    ##
    ## @end deftp
    Mu              = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Coeffs
    ##
    ## Coefficient matrices
    ##
    ## A @math{K*K} structure containing the coefficient matrices, where
    ## @math{K} is the number of classes.  If the @qcode{'FillCoeffs'} parameter
    ## was set to @qcode{'off'} in either the @code{fitcdiscr} function or the
    ## @code{ClassificationDiscriminant} constructor, then @qcode{Coeffs} is
    ## empty @qcode{([])}.  This property is read-only.
    ##
    ## @qcode{Coeffs(i,j)} contains the coefficients of the boundary between
    ## the classes @code{i} and @code{j} in the following fields:
    ##
    ## @itemize
    ## @item @qcode{DiscrimType} - A character vector
    ## @item @qcode{Class1} - @qcode{@var{ClassNames}(i)}
    ## @item @qcode{Class2} - @qcode{@var{ClassNames}(j)}
    ## @item @qcode{Const} - A scalar
    ## @item @qcode{Linear} - A vector with length as the number of predictors.
    ## @item @qcode{Quadratic} - The quadratic family only.  A @math{PxP}
    ## matrix, or a @math{1xP} vector for @qcode{'diagQuadratic'}, following
    ## the shape of @qcode{Sigma}.
    ## @end itemize
    ##
    ## The diagonal entries carry the two class names and nothing else.  The
    ## structure is rebuilt whenever @qcode{DiscrimType}, @qcode{Gamma} or
    ## @qcode{Prior} is assigned.
    ##
    ## @end deftp
    Coeffs          = [];

    ## -*- texinfo -*-
    ## @deftp {%s} {property} DeltaPredictor
    ##
    ## Minimum Delta at which each predictor drops out
    ##
    ## A row vector with one entry per predictor, the value of @qcode{Delta} at
    ## which that predictor's coefficient is zero for every class and the
    ## predictor leaves the model altogether.  It is all zeros for the
    ## quadratic family, which has no linear coefficients to eliminate.
    ##
    ## This property is read-only, and it describes the fit rather than the
    ## threshold: assigning @qcode{Delta} does not move it.
    ##
    ## @end deftp
    DeltaPredictor  = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} MinGamma
    ##
    ## Minimum value for the Gamma regularization parameter
    ##
    ## A scalar from 0 to 1, the least regularization that leaves the
    ## correlation matrix invertible.  It is 0 when the matrix is already
    ## invertible, and positive when the predictors are collinear, in which
    ## case a plain @qcode{'linear'} or @qcode{'quadratic'} fit is raised to it
    ## rather than failing.  Assigning a @qcode{Gamma} below it is refused.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    MinGamma        = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} LogDetSigma
    ##
    ## Logarithm of the determinant of the within-class covariance matrix
    ##
    ## A scalar for the linear family and a @math{Kx1} vector for the quadratic
    ## one, one entry per class.  It is computed in correlation space, as the
    ## sum of the logarithms of the predictor variances plus the log
    ## determinant of the correlation matrix, which is far better conditioned
    ## than the covariance when the data are nearly collinear.  A predictor
    ## with no variance contributes nothing rather than an infinity, and the
    ## @qcode{'pseudo'} types sum only over the directions that carry variance.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    LogDetSigma     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} XCentered
    ##
    ## Predictor data with class means subtracted
    ##
    ## A matrix of the same size as @var{X} and the values in @var{X} with the
    ## corresponding class means subtracted.  This property is read-only.
    ##
    ## @end deftp
    XCentered       = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} BinEdges
    ##
    ## Bin edges of the predictors
    ##
    ## A cell array with one entry per predictor, holding that predictor's bin
    ## edges where the learner discretized it before fitting.  It is empty here
    ## and stays empty: this learner fits the predictors as they are, and
    ## MATLAB's reports an empty cell for it as well.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    BinEdges        = {};

  endproperties

  ## Properties a user may set after the model is fitted.  Each one is
  ## validated by its set method below.
  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} DiscrimType
    ##
    ## Discriminant type
    ##
    ## A character vector naming the discriminant model, one of
    ## @qcode{'linear'}, @qcode{'quadratic'}, @qcode{'diagLinear'},
    ## @qcode{'diagQuadratic'}, @qcode{'pseudoLinear'} or
    ## @qcode{'pseudoQuadratic'}.  A linear type pools one covariance across
    ## the classes; a quadratic type estimates one per class.  A
    ## @qcode{'diag'} type keeps only the variances, and a @qcode{'pseudo'}
    ## type inverts a singular covariance instead of refusing it.
    ##
    ## This property may be assigned, but @emph{only within its own family}:
    ## the three linear types interchange freely and so do the three quadratic
    ## ones, while no assignment moves a model between the two.  The family is
    ## fixed when the model is fitted, because it decides which covariances the
    ## fit has to estimate.  Assigning re-derives @qcode{Sigma},
    ## @qcode{LogDetSigma}, @qcode{Gamma} and @qcode{Coeffs}.
    ##
    ## @end deftp
    DiscrimType     = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Gamma
    ##
    ## Gamma regularization parameter
    ##
    ## A scalar from 0 to 1 shrinking the covariance towards its diagonal.
    ## @qcode{Gamma} and @qcode{DiscrimType} are one state: a value of 1 is the
    ## diagonal type, so assigning it renames @qcode{DiscrimType} to
    ## @qcode{'diagLinear'} or @qcode{'diagQuadratic'}, and assigning a
    ## diagonal type sets @qcode{Gamma} to 1.
    ##
    ## The quadratic family admits 0 and 1 only.  A value below
    ## @qcode{MinGamma} is refused, since it would leave the covariance
    ## singular.  Assigning re-derives @qcode{Sigma}, @qcode{LogDetSigma} and
    ## @qcode{Coeffs}.
    ##
    ## @end deftp
    Gamma           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Delta
    ##
    ## Delta threshold for the linear coefficients
    ##
    ## A nonnegative scalar that eliminates predictors.  A per-class linear
    ## coefficient is set to zero when it falls below @qcode{Delta}, and the
    ## comparison is made on the @strong{standardized} coefficient, the
    ## coefficient times the within-class standard deviation of its predictor.
    ## Scaling matters here: a threshold on the raw coefficients would depend
    ## on the units each predictor is measured in, so the same model in
    ## centimetres and in metres would drop different predictors.
    ##
    ## @qcode{DeltaPredictor} reports, per predictor, the value at which it
    ## drops out of every class at once.
    ##
    ## It applies to the linear family only, a quadratic discriminant having no
    ## linear coefficients to eliminate.  Assigning it rebuilds @qcode{Coeffs}
    ## and changes what @code{predict} answers.
    ##
    ## @end deftp
    Delta           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Cost
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
    ## Add or change the @qcode{Cost} property using dot notation as in:
    ## @itemize
    ## @item @qcode{@var{obj}.Cost = @var{costMatrix}}
    ## @end itemize
    ##
    ## @end deftp
    Cost            = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} Prior
    ##
    ## Prior probability for each class
    ##
    ## A numeric vector specifying the prior probabilities for each class.  The
    ## order of the elements in @qcode{Prior} corresponds to the order of the
    ## classes in @qcode{ClassNames}.
    ##
    ## Add or change the @qcode{Prior} property using dot notation as in:
    ## @itemize
    ## @item @qcode{@var{obj}.Prior = @var{priorVector}}
    ## @end itemize
    ##
    ## Specified as a row vector with one entry per class, in the order of
    ## @qcode{ClassNames}, and rescaled to sum to one.  It may be given as
    ## @qcode{'empirical'}, @qcode{'uniform'}, a numeric vector, or a
    ## structure with @qcode{ClassNames} and @qcode{ClassProbs} fields, which
    ## assigns each probability by class name rather than by position.
    ##
    ## @end deftp
    Prior           = [];

    ## -*- texinfo -*-
    ## @deftp {ClassificationDiscriminant} {property} ScoreTransform
    ##
    ## Transformation function for classification scores
    ##
    ## Specified as a function handle for transforming the classification
    ## scores.  Add or change the @qcode{ScoreTransform} property using dot
    ## notation as in:
    ##
    ## @itemize
    ## @item @qcode{@var{obj}.ScoreTransform = 'function_name'}
    ## @item @qcode{@var{obj}.ScoreTransform = @@function_handle}
    ## @end itemize
    ##
    ## When specified as a character vector, it can be any of the following
    ## built-in functions.  Nevertheless, the @qcode{ScoreTransform} property
    ## always stores their function handle equivalent.
    ##
    ## @multitable @columnfractions 0.2 0.75
    ## @headitem @var{Value} @tab @var{Description}
    ## @item @qcode{'doublelogit'} @tab @math{1 ./ (1 + exp (-2 * x))}
    ## @item @qcode{'invlogit'} @tab @math{log (x ./ (1 - x))}
    ## @item @qcode{'ismax'} @tab Sets the score for the class with the
    ## largest score to 1, and for all other classes to 0
    ## @item @qcode{'logit'} @tab @math{1 ./ (1 + exp (-x))}
    ## @item @qcode{'none'} @tab @math{x} (no transformation)
    ## @item @qcode{'identity'} @tab @math{x} (no transformation)
    ## @item @qcode{'sign'} @tab
    ## @math{-1 for x < 0, 0 for x = 0, 1 for x >
    ## 0}
    ## @item @qcode{'symmetric'} @tab @math{2 * x - 1}
    ## @item @qcode{'symmetricismax'} @tab Sets the score for the class
    ## with the largest score to 1, and for all other classes to -1
    ## @item @qcode{'symmetriclogit'} @tab @math{2 ./ (1 + exp (-x)) - 1}
    ## @end multitable
    ##
    ## @end deftp
    ScoreTransform  = 'none';

  endproperties

  ## Readable by the compact counterpart, which copies it, and kept out of
  ## the documented surface.  MATLAB hides its own equivalent the same way.
  properties (GetAccess = public, SetAccess = protected, Hidden)
    STfun = @(x) x;

    ## The unregularized within-class covariance the fit estimated, PxP for the
    ## linear family and PxPxK for the quadratic one.  Sigma, LogDetSigma and
    ## Coeffs are all derived from it, so assigning DiscrimType within the
    ## family re-derives rather than refits.
    BaseSigma = [];
  endproperties

  ## Set methods for the properties a user may assign after fitting.
  methods (Hidden)

    function this = set.Cost (this, Cost)
      [~, gnY] = unique (this.Y);
      if (isempty (Cost))
        this.Cost = cast (! eye (numel (gnY)), 'double');
      else
        if (numel (gnY) != sqrt (numel (Cost)))
          error (strcat ("ClassificationDiscriminant: the number", ...
                         " of rows and columns in 'Cost' must", ...
                         " correspond to selected classes in Y."));
        endif
        this.Cost = Cost;
      endif
    endfunction

    function this = set.Prior (this, Prior)
      [~, gnY, gY] = unique (this.Y);
      if (isstruct (Prior))
        Prior = priorFromStruct (Prior, this.ClassNames, ...
                                 'ClassificationDiscriminant');
      endif
      ## Set prior
      if (strcmpi ('uniform', Prior))
        this.Prior = ones (1, numel (gnY)) ./ numel (gnY);
      elseif (isempty (Prior) || strcmpi ('empirical', Prior))
        pr = [];
        for i = 1:numel (gnY)
          pr = [pr; sum(gY==i)];
        endfor
        this.Prior = pr(:)' ./ sum (pr);
      elseif (isnumeric (Prior))
        if (numel (gnY) != numel (Prior))
          error (strcat ("ClassificationDiscriminant: the elements", ...
                         " in 'Prior' do not correspond to the", ...
                         " selected classes in Y."));
        endif
        this.Prior = Prior(:)' ./ sum (Prior);
      endif
      ## Rebuild the coefficients, whose constant term reads the priors.
      ## Delegating rather than restating the linear formula is what keeps the
      ## six discriminant types in step: the formula here was linear only, and
      ## a quadratic model's Sigma is not even the same shape.
      if (! isempty (this.Coeffs) && ! isempty (this.BaseSigma))
        ## Derived here rather than read off Sigma, so the covariance and the
        ## type always agree even when the object is part way through a load.
        [Sg, Ld] = discrimderive (this.BaseSigma, this.DiscrimType, ...
                                  this.Gamma, this.MinGamma);
        this.Coeffs = discrimcoeffs (this.Mu, Sg, Ld, this.Prior, ...
                                     this.DiscrimType, this.Delta, ...
                                     this.ClassNames);
      endif
    endfunction

    ## DiscrimType, Gamma and Delta are one state, and these three methods are
    ## the only place it is changed.  Assigning the type or the regularization
    ## re-derives Sigma, LogDetSigma and Coeffs from BaseSigma, the covariance
    ## the fit estimated, so no assignment ever refits.  Each writes its own
    ## property directly, which does not re-enter its own set method, and
    ## reaches the others through theirs, which is what keeps the three
    ## consistent whichever one the user assigns.
    function this = set.DiscrimType (this, val)
      [t, fam] = discrimcanon (val);
      ## The family a fit belongs to is written in BaseSigma's shape: one
      ## covariance for the linear family, one per class for the quadratic.
      ## Reading it there rather than from the stored type is what lets a
      ## loaded model take its own type back, the stub it is loaded into
      ## carrying the default type until the file has been read.
      if (! isempty (this.BaseSigma) && size (this.BaseSigma, 3) > 1)
        cur = 'quadratic';
      elseif (! isempty (this.BaseSigma))
        cur = 'linear';
      else
        [~, cur] = discrimcanon (this.DiscrimType);
      endif
      if (isempty (t) || (! isempty (cur) && ! strcmp (fam, cur)))
        ## The family was fixed by the fit, so only three types are on offer,
        ## and an unrecognized name is refused by the same message.
        if (strcmp (cur, 'quadratic'))
          ok = "quadratic, diagQuadratic, or pseudoQuadratic";
        else
          ok = "linear, diagLinear, or pseudoLinear";
        endif
        error (strcat ("ClassificationDiscriminant: 'DiscrimType' can only", ...
                       " be set to one of: %s."), ok);
      endif
      this.DiscrimType = t;
      if (isempty (this.BaseSigma))
        return;   # not fitted yet, nothing to derive
      endif
      ## A diagonal type is Gamma of 1 and every other type starts from 0,
      ## which is what the oracle reports for each transition.  The derivation
      ## raises it to MinGamma where the covariance needs it.
      g = double (strncmp (t, 'diag', 4));
      [~, ~, g] = discrimderive (this.BaseSigma, t, g, this.MinGamma);
      this.Gamma = g;
    endfunction

    function this = set.Gamma (this, val)
      if (! (isnumeric (val) && isscalar (val) && val >= 0 && val <= 1))
        error (strcat ("ClassificationDiscriminant: 'Gamma' must be a", ...
                       " scalar between 0 and 1."));
      endif
      if (isempty (this.BaseSigma))
        this.Gamma = val;
        return;
      endif
      [~, fam] = discrimcanon (this.DiscrimType);
      if (strcmp (fam, 'quadratic') && val > 0 && val < 1)
        error (strcat ("ClassificationDiscriminant: cannot set 'Gamma' to", ...
                       " any value but 0 or 1 for a quadratic discriminant."));
      endif

      ## Gamma of 1 IS the diagonal type.  Naming the type does the rest, and
      ## comes back here with a value that no longer needs the rename.
      if (val == 1 && ! strncmp (this.DiscrimType, 'diag', 4))
        if (strcmp (fam, 'quadratic'))
          this.DiscrimType = 'diagQuadratic';
        else
          this.DiscrimType = 'diagLinear';
        endif
        return;
      elseif (val < 1 && strncmp (this.DiscrimType, 'diag', 4))
        if (strcmp (fam, 'quadratic'))
          this.DiscrimType = 'quadratic';
        else
          this.DiscrimType = 'linear';
        endif
      endif

      ## The MinGamma floor belongs to the one type that regularizes its way
      ## out of a singular covariance.  A pseudo type is Gamma of 0 by
      ## definition, inverting the rank it has instead; a diagonal type is
      ## Gamma of 1 and never near the floor; and the quadratic family refuses
      ## a singular class covariance outright rather than ridging it.
      Mg = this.MinGamma;
      [Sg, Ld, Ge] = discrimderive (this.BaseSigma, this.DiscrimType, val, Mg);
      if (val < Mg && strcmp (this.DiscrimType, 'linear'))
        error (strcat ("ClassificationDiscriminant: 'Gamma' must be", ...
                       " between %g and 1."), Mg);
      endif
      this.Gamma = Ge;
      this.Sigma = Sg;
      this.LogDetSigma = Ld;
      if (! isempty (this.Coeffs))
        this.Coeffs = discrimcoeffs (this.Mu, Sg, Ld, this.Prior, ...
                                     this.DiscrimType, this.Delta, ...
                                     this.ClassNames);
      endif
    endfunction

    function this = set.Delta (this, val)
      if (! (isnumeric (val) && isscalar (val) && val >= 0))
        error (strcat ("ClassificationDiscriminant: 'Delta' must be a", ...
                       " nonnegative scalar."));
      endif
      [~, fam] = discrimcanon (this.DiscrimType);
      if (val > 0 && strcmp (fam, 'quadratic'))
        error (strcat ("ClassificationDiscriminant: cannot eliminate", ...
                       " linear predictors in a quadratic discriminant."));
      endif
      this.Delta = val;
      if (! isempty (this.BaseSigma) && ! isempty (this.Coeffs))
        this.Coeffs = discrimcoeffs (this.Mu, this.Sigma, this.LogDetSigma, ...
                                     this.Prior, this.DiscrimType, val, ...
                                     this.ClassNames);
      endif
    endfunction

    function this = set.ScoreTransform (this, val)
      [f, nm] = parseScoreTransform (val, 'ClassificationDiscriminant');
      this.ScoreTransform = nm;
      this.STfun = f;
    endfunction

    ## Custom display
    function display (this)
      in_name = inputname (1);
      if (! isempty (in_name))
        fprintf ('%s =\n', in_name);
      endif
      disp (this);
    endfunction

    ## Custom display
    function disp (this)
      fprintf ("\n  ClassificationDiscriminant\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      if (iscellstr (this.ClassNames))
        str = repmat ({'''%s'''}, 1, numel (this.ClassNames));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.ClassNames{:});
      elseif (ischar (this.ClassNames))
        str = repmat ({'''%s'''}, 1, rows (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, cellstr (this.ClassNames){:});
      else # single, double, logical
        str = repmat ({'%d'}, 1, numel (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, this.ClassNames);
      endif
      fprintf ("%+25s: %s\n", 'ClassNames', str);
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("%+25s: %d\n", 'NumObservations', this.NumObservations);
      fprintf ("%+25s: %d\n", 'NumPredictors', this.NumPredictors);
      fprintf ("%+25s: '%s'\n", 'DiscrimType', this.DiscrimType);
      fprintf ("%+25s: [%dx%d double]\n", 'Mu', size (this.Mu));
      ## Coeffs is KxK, which is Sigma's shape only for a linear discriminant.
      fprintf ("%+25s: [%dx%d struct]\n\n", 'Coeffs', ...
               rows (this.ClassNames), rows (this.ClassNames));
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {statistics} {@var{obj} =} ClassificationDiscriminant (@var{X}, @var{Y})
    ## @deftypefnx {statistics} {@var{obj} =} ClassificationDiscriminant (@dots{}, @var{name}, @var{value})
    ##
    ## Create a @qcode{ClassificationDiscriminant} class object containing a
    ## discriminant analysis model.
    ##
    ## @code{@var{obj} = ClassificationDiscriminant (@var{X}, @var{Y})} returns
    ## a ClassificationDiscriminant object, with @var{X} as the predictor data
    ## and @var{Y} containing the class labels of observations in @var{X}.
    ##
    ## @itemize
    ## @item
    ## @code{X} must be a @math{N*P} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.  @var{X} will be used to train the discriminant model.
    ## @item
    ## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}.  @var{Y} can contain any type
    ## of categorical data. @var{Y} must have the same number of rows as
    ## @var{X}.
    ## @end itemize
    ##
    ## @code{@var{obj} = ClassificationDiscriminant (@dots{}, @var{name},
    ## @var{value})} returns a ClassificationDiscriminant object with parameters
    ## specified by the following @qcode{@var{name}, @var{value}} paired input
    ## arguments:
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'PredictorNames'} @tab A cell array of character
    ## vectors specifying the names of the predictors. The length of this array
    ## must match the number of columns in @var{X}.
    ##
    ## @item @qcode{'ResponseName'} @tab A character vector specifying the
    ## name of the response variable.
    ##
    ## @item @qcode{'ClassNames'} @tab Names of the classes in the class
    ## labels, @var{Y}, used for fitting the Discriminant model.
    ## @qcode{ClassNames} are of the same type as the class labels in @var{Y}.
    ##
    ## @item @qcode{'Cost'} @tab An @math{N*R} numeric matrix containing
    ## misclassification cost for the corresponding instances in @var{X}, where
    ## @math{R} is the number of unique categories in @var{Y}.  If an instance
    ## is correctly classified into its category the cost is calculated to be 1,
    ## otherwise 0. The cost matrix can be altered by using
    ## @code{@var{Mdl}.cost = somecost}.  By default, its value is
    ## @qcode{@var{cost} = ones (rows (X), numel (unique (Y)))}.
    ##
    ## @item @qcode{'Prior'} @tab A numeric vector specifying the prior
    ## probabilities for each class.  The order of the elements in @qcode{Prior}
    ## corresponds to the order of the classes in @qcode{ClassNames}.
    ## Alternatively, you can specify @qcode{'empirical'} to use the empirical
    ## class probabilities or @qcode{'uniform'} to assume equal class
    ## probabilities.
    ##
    ## @item @qcode{'ScoreTransform'} @tab A user-defined function handle
    ## or a character vector specifying one of the following builtin functions
    ## specifying the transformation applied to predicted classification scores.
    ## Supported values include @qcode{'doublelogit'}, @qcode{'invlogit'},
    ## @qcode{'ismax'}, @qcode{'logit'}, @qcode{'none'}, @qcode{'identity'},
    ## @qcode{'sign'}, @qcode{'symmetric'}, @qcode{'symmetricismax'}, and
    ## @qcode{'symmetriclogit'}.
    ##
    ## @item @qcode{'DiscrimType'} @tab A character vector or string scalar
    ## specifying the type of discriminant analysis to perform. The only
    ## supported value is @qcode{'linear'}.
    ##
    ## @item @qcode{'FillCoeffs'} @tab A character vector or string scalar
    ## with values @qcode{'on'} or @qcode{'off'} specifying whether to fill the
    ## coefficients after fitting. If set to @qcode{'on'}, the coefficients are
    ## computed during model fitting, which can be useful for prediction.
    ##
    ## @item @qcode{'Gamma'} @tab A numeric scalar specifying the
    ## regularization parameter for the covariance matrix. It adjusts the linear
    ## discriminant analysis to make the model more stable in the presence of
    ## multicollinearity or small sample sizes. A value of 0 corresponds to no
    ## regularization, while a value of 1 corresponds to
    ## a completely regularized model.
    ## @end multitable
    ##
    ## @seealso{fitcdiscr}
    ## @end deftypefn
    function this = ClassificationDiscriminant (X, Y, varargin)

      ## Check for appropriate number of input arguments
      if (nargin < 2)
        error ("ClassificationDiscriminant: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationDiscriminant: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      ## Validate X
      if (! isnumeric (X))
        error ("ClassificationDiscriminant: X must be a numeric matrix.");
      endif

      ## Check X and Y have the same number of observations
      if (rows (X) != rows (Y))
        error (strcat ("ClassificationDiscriminant: number", ...
                       " of rows in X and Y must be equal."));
      endif

      ## Assign original X and Y data
      this.X = X;
      this.Y = Y;

      ## Get groups in Y
      [gY, gnY, glY] = grp2idx (Y);

      ## Set default values before parsing optional parameters
      ClassNames     = [];
      Cost           = [];
      DiscrimType    = 'linear';
      Gamma          = 0;
      Delta          = 0;
      NumPredictors  = [];
      PredictorNames = {};
      ResponseName   = 'Y';
      Prior          = 'empirical';
      FillCoeffs     = 'on';

      ## Parse optional parameters
      while (numel (varargin) > 0)
        switch (lower (varargin{1}))

          case 'predictornames'
            PredictorNames = varargin{2};
            if (! iscellstr (PredictorNames))
              error (strcat ("ClassificationDiscriminant: 'PredictorNames'", ...
                             " must be supplied as a cellstring array."));
            elseif (numel (PredictorNames) != columns (X))
              error (strcat ("ClassificationDiscriminant: 'PredictorNames'", ...
                             " must equal the number of columns in X."));
            endif

          case 'responsename'
            ResponseName = varargin{2};
            if (! ischar (ResponseName))
              error (strcat ("ClassificationDiscriminant: 'ResponseName'", ...
                             " must be a character vector."));
            endif

          case 'classnames'
            ClassNames = varargin{2};
            if (! (iscellstr (ClassNames) || isnumeric (ClassNames) ||
                   islogical (ClassNames) || ischar (ClassNames)))
              error (strcat ("ClassificationDiscriminant: 'ClassNames'", ...
                             " must be a cell array of character vectors,", ...
                             " a logical vector, a numeric vector,", ...
                             " or a character array."));
            endif
            ## Check that all class names are available in gnY
            if (iscellstr (ClassNames) || ischar (ClassNames))
              ClassNames = cellstr (ClassNames);
              if (! all (cell2mat (cellfun (@(x) any (strcmp (x, gnY)),
                                   ClassNames, 'UniformOutput', false))))
                error (strcat ("ClassificationDiscriminant: not all", ...
                               " 'ClassNames' are present in Y."));
              endif
            else
              if (! all (cell2mat (arrayfun (@(x) any (x == glY),
                                   ClassNames, 'UniformOutput', false))))
                error (strcat ("ClassificationDiscriminant: not all", ...
                               " 'ClassNames' are present in Y."));
              endif
            endif

          case 'prior'
            Prior = varargin{2};
            if (! (isstruct (Prior) || (isnumeric (Prior) && isvector (Prior))
                   || (ischar (Prior) && (strcmpi (Prior, 'empirical')
                       || strcmpi (Prior, 'uniform')))))
              error (strcat ("ClassificationDiscriminant: 'Prior' must", ...
                             " be either a numeric or a character vector."));
            endif

          case 'cost'
            Cost = varargin{2};
            if (! (isnumeric (Cost) && issquare (Cost)))
              error (strcat ("ClassificationDiscriminant: 'Cost'", ...
                             " must be a numeric square matrix."));
            endif

          case 'scoretransform'
            this.ScoreTransform = varargin{2};

          case 'discrimtype'
            DiscrimType = discrimcanon (varargin{2});
            if (isempty (DiscrimType))
              error (strcat ("ClassificationDiscriminant: 'DiscrimType'", ...
                             " must be one of the following: linear,", ...
                             " quadratic, diagLinear, diagQuadratic,", ...
                             " pseudoLinear, or pseudoQuadratic."));
            endif

          case 'fillcoeffs'
            FillCoeffs = tolower (varargin{2});
            if (! any (strcmpi (FillCoeffs, {'on', 'off'})))
              error (strcat ("ClassificationDiscriminant: 'FillCoeffs'", ...
                             " must be 'on' or 'off'."));
            endif

          case 'gamma'
            Gamma = varargin{2};
            if (! (isnumeric (Gamma) && isscalar (Gamma)
                   && Gamma >= 0 && Gamma <= 1))
              error (strcat ("ClassificationDiscriminant: 'Gamma'", ...
                             " must be a scalar between 0 and 1."));
            endif

          case 'delta'
            Delta = varargin{2};
            if (! (isnumeric (Delta) && isscalar (Delta) && Delta >= 0))
              error (strcat ("ClassificationDiscriminant: 'Delta'", ...
                             " must be a nonnegative scalar."));
            endif


          otherwise
            error (strcat ("ClassificationDiscriminant: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      ## Generate default predictors and response variable names (if necessary)
      NumPredictors = columns (X);
      if (isempty (PredictorNames))
        for i = 1:NumPredictors
          PredictorNames {i} = strcat ("x", num2str (i));
        endfor
      endif
      if (isempty (ResponseName))
        ResponseName = 'Y';
      endif

      ## Assign predictors and response variable names
      this.NumPredictors  = NumPredictors;
      this.PredictorNames = PredictorNames;
      this.CategoricalPredictors = [];
      this.ExpandedPredictorNames = PredictorNames;
      this.ResponseName   = ResponseName;

      ## Handle class names
      if (! isempty (ClassNames))
        if (iscellstr (ClassNames))
          ru = find (! ismember (gnY, ClassNames));
        else
          ru = find (! ismember (glY, ClassNames));
        endif
        for i = 1:numel (ru)
          gY(gY == ru(i)) = NaN;
        endfor
      endif

      ## An observation is dropped only when its response is missing.  A row
      ## whose predictors hold missing values is kept, and each estimate below
      ## uses whatever part of it is present.
      RowsUsed  = ! isnan (gY);
      Y         = Y(RowsUsed, :);
      X         = X(RowsUsed, :);

      ## Store the retained observations
      this.X = X;
      this.Y = Y;

      ## Renew groups in Y, get classes ordered, keep the same type
      [this.ClassNames, gnY, gY] = unique (Y);

      ## Check X contains valid data
      if (! (isnumeric (X) && ! any (isinf (X(:)))))
        error ("ClassificationDiscriminant: invalid values in X.");
      endif

      ## Assign the number of observations and their corresponding indices
      ## on the original data, which will be used for training the model,
      ## to the ClassificationDiscriminant object
      this.NumObservations = rows (X);
      ## RowsUsed is left empty when every observation was used, as in MATLAB
      if (all (RowsUsed))
        this.RowsUsed = [];
      else
        this.RowsUsed = RowsUsed;
      endif

      ## Handle Cost and Prior
      this.Cost  = Cost;
      this.Prior = Prior;
      ## A discriminant weighs every observation alike; the prior enters
      ## prediction rather than the fit, which is what MATLAB reports.
      this.W = ones (this.NumObservations, 1) / this.NumObservations;

      ## Assign DiscrimType
      ## Reconcile the type with the regularization before anything is
      ## estimated.  Gamma of 1 IS the diagonal type, in both directions, and
      ## the quadratic family admits neither an intermediate Gamma nor a Delta.
      [~, fam] = discrimcanon (DiscrimType);
      if (Gamma == 1 && ! strncmp (DiscrimType, 'diag', 4))
        if (strcmp (fam, 'linear'))
          DiscrimType = 'diagLinear';
        else
          DiscrimType = 'diagQuadratic';
        endif
      elseif (Gamma > 0 && Gamma < 1 && strcmp (fam, 'quadratic'))
        error (strcat ("ClassificationDiscriminant: cannot set 'Gamma' to", ...
                       " any value but 0 or 1 for a quadratic discriminant."));
      endif
      if (Delta > 0 && strcmp (fam, 'quadratic'))
        error (strcat ("ClassificationDiscriminant: cannot set linear", ...
                       " coefficients to zero for a quadratic discriminant."));
      endif

      this.DiscrimType = DiscrimType;
      this.Delta = Delta;
      this.Gamma = Gamma;

      num_classes = rows (this.ClassNames);
      num_features = columns (X);
      ## Each class mean uses every observation where that predictor is
      ## present, so a row missing one predictor still counts towards the rest
      this.Mu = zeros (num_classes, num_features);
      for i = 1:num_classes
        Xi = X(gY == i, :);
        for j = 1:num_features
          xj = Xi(:, j);
          this.Mu(i, j) = mean (xj(! isnan (xj)));
        endfor
      endfor

      ## The between-class covariance of the means about the overall mean,
      ## weighted by class size and divided by the unbiased weighted-covariance
      ## denominator.  It reads the class sizes rather than Prior, which is why
      ## assigning a prior does not move it.
      nk = accumarray (gY(:), 1, [num_classes, 1]);
      pk = nk ./ sum (nk);
      D = this.Mu - pk' * this.Mu;
      this.BetweenSigma = (D' * (D .* nk)) ./ (sum (nk) * (1 - sum (pk .^ 2)));

      ## Center the predictors (XCentered), keeping the missing entries
      this.XCentered = X - this.Mu(gY, :);

      ## Estimate the covariance the family calls for.  A linear family pools
      ## one covariance across the classes and a quadratic family estimates one
      ## per class; every type in a family is derived from what is estimated
      ## here, which is what lets DiscrimType be assigned after the fit without
      ## the model having to keep both.
      [~, fam] = discrimcanon (this.DiscrimType);
      cobs = ! any (isnan (X), 2);

      ## The pooled covariance is estimated whatever the family, because
      ## MinGamma is taken from it.  Only complete observations enter it,
      ## reweighted so that each class keeps the total weight it carried
      ## before any was dropped.
      cw   = zeros (rows (X), 1);
      Wk   = zeros (num_classes, 1);
      for i = 1:num_classes
        Wk(i) = sum (gY == i) / rows (X);
        ci = (gY == i) & cobs;
        cw(ci) = Wk(i) / sum (ci);
      endfor
      den = 1;
      for i = 1:num_classes
        ci = (gY == i) & cobs;
        den -= sum (cw(ci) .^ 2) / Wk(i);
      endfor
      Zc = this.XCentered(cobs, :);
      pooled = (Zc .* cw(cobs))' * Zc / den;
      this.MinGamma = discrimmingamma (pooled);

      if (strcmp (fam, 'linear'))
        this.BaseSigma = pooled;

        ## A predictor with no within-class variance cannot be inverted, so the
        ## plain type refuses it and names the two types that can take it.
        zwcv = find (diag (this.BaseSigma) == 0);
        if (! isempty (zwcv) && strcmp (this.DiscrimType, 'linear'))
          error (strcat ("ClassificationDiscriminant: predictor '%s'", ...
                         " has zero within-class variance.  Either exclude", ...
                         " this predictor or set 'DiscrimType' to", ...
                         " 'pseudoLinear' or 'diagLinear'."), ...
                 PredictorNames{zwcv(1)});
        endif
      else
        this.BaseSigma = zeros (num_features, num_features, num_classes);
        for i = 1:num_classes
          ci = (gY == i) & cobs;
          Zi = this.XCentered(ci, :);
          this.BaseSigma(:,:,i) = (Zi' * Zi) / (sum (ci) - 1);

          zwcv = find (diag (this.BaseSigma(:,:,i)) == 0);
          if (! isempty (zwcv) && strcmp (this.DiscrimType, 'quadratic'))
            ## ClassNames keeps the type of Y, so the name is rendered rather
            ## than indexed as though it were always a cell.
            if (iscellstr (this.ClassNames))
              cname = this.ClassNames{i};
            elseif (ischar (this.ClassNames))
              cname = this.ClassNames(i,:);
            else
              cname = num2str (this.ClassNames(i));
            endif
            error (strcat ("ClassificationDiscriminant: predictor '%s' has", ...
                           " zero variance for class '%s'.  Either exclude", ...
                           " this predictor or set 'DiscrimType' to", ...
                           " 'pseudoQuadratic' or 'diagQuadratic'."), ...
                   PredictorNames{zwcv(1)}, cname);
          endif
        endfor
      endif

      ## Derive Sigma, its log determinant and the regularization from the
      ## estimate above.  The type and Gamma are one state, so both come back.
      [this.Sigma, this.LogDetSigma, this.Gamma] = ...
              discrimderive (this.BaseSigma, this.DiscrimType, this.Gamma, ...
                             this.MinGamma);

      ## A singular class covariance is what the plain quadratic type cannot
      ## take; the pseudo and diagonal ones answer it.
      if (strcmp (this.DiscrimType, 'quadratic'))
        for i = 1:num_classes
          if (rcond (this.Sigma(:,:,i)) < num_features * eps)
            error (strcat ("ClassificationDiscriminant: cannot use", ...
                           " 'quadratic' type because one or more classes", ...
                           " have singular covariance matrices."));
          endif
        endfor
      endif

      ## Delta is a threshold on the linear coefficients, so the quadratic
      ## family reports it as zeros and has nothing to eliminate.
      if (strcmp (fam, 'linear'))
        [~, ~, this.DeltaPredictor] = discrimlinear (this.Mu, this.Sigma, ...
                                       this.Prior, this.DiscrimType, 0);
      else
        this.DeltaPredictor = zeros (1, num_features);
      endif

      if (strcmpi (FillCoeffs, 'on'))
        this.Coeffs = discrimcoeffs (this.Mu, this.Sigma, this.LogDetSigma, ...
                                     this.Prior, this.DiscrimType, ...
                                     this.Delta, this.ClassNames);
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {ClassificationDiscriminant} {[@var{label}, @var{score}, @var{cost}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the discriminant
    ## analysis model from a ClassificationDiscriminant object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## predictor data in @code{obj.X} and corresponding labels, @code{obj.Y},
    ## stored in the ClassificationDiscriminant model, @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{ClassificationDiscriminant} class object.
    ## @item
    ## @var{XC} must be an @math{M*P} numeric matrix with the same number of
    ## features @math{P} as the corresponding predictors of the discriminant
    ## model in @var{obj}.
    ## @end itemize
    ##
    ## @code{[@var{label}, @var{score}, @var{cost}] = predict (@var{obj},
    ## @var{XC})} also returns @var{score}, which contains the predicted class
    ## scores or posterior probabilities for each instance of the corresponding
    ## unique classes, and @var{cost}, which is a matrix containing the expected
    ## cost of the classifications.
    ##
    ## The @var{score} matrix contains the posterior probabilities for each
    ## class, calculated using the multivariate normal probability density
    ## function and the prior probabilities of each class.  These scores are
    ## normalized to ensure they sum to 1 for each observation.
    ##
    ## The @var{cost} matrix contains the expected classification cost for each
    ## class, computed based on the posterior probabilities and the specified
    ## misclassification costs.
    ##
    ## @seealso{ClassificationDiscriminant, fitcdiscr}
    ## @end deftypefn
    function [label, score, cost] = predict (this, XC)
      ## Check for sufficient input arguments
      if (nargin < 2)
        error ("ClassificationDiscriminant.predict: too few input arguments.");
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("ClassificationDiscriminant.predict: XC is empty.");
      elseif (columns (this.X) != columns (XC))
        error (strcat ("ClassificationDiscriminant.predict: XC must have ", ...
                       " the same number of predictors as the trained model."));
      endif

      ## Get training data and labels
      X = this.X;
      Y = this.Y;

      numObservations = rows (XC);
      numClasses = rows (this.ClassNames);
      score = zeros (numObservations, numClasses);
      cost = zeros (numObservations, numClasses);

      ## Score from the inverse covariance and its log determinant rather than
      ## through mvnpdf.  A pseudo type's covariance is deliberately singular
      ## and mvnpdf refuses it, and LogDetSigma has to be the value the
      ## property reports, so the score and the property cannot drift apart.
      [~, fam] = discrimcanon (this.DiscrimType);
      logscore = zeros (numObservations, numClasses);
      if (strcmp (fam, 'linear'))
        ## The linear family scores from its per-class coefficients, which is
        ## the same function with the terms common to every class dropped, and
        ## is the only form in which Delta can eliminate a predictor.
        [Z, b] = discrimlinear (this.Mu, this.Sigma, this.Prior, ...
                                this.DiscrimType, this.Delta);
        logscore = XC * Z' + b;
      else
        SigmaInv = discriminv (this.Sigma, this.DiscrimType);
        nInv = size (SigmaInv, 3);
        logdet = this.LogDetSigma;
        if (isscalar (logdet))
          logdet = repmat (logdet, numClasses, 1);
        endif
        for i = 1:numClasses
          Si = SigmaInv(:,:,min (i, nInv));
          Zc = XC - this.Mu(i, :);
          logscore(:, i) = -0.5 * sum ((Zc * Si) .* Zc, 2) ...
                           - 0.5 * logdet(i) + log (this.Prior(i));
        endfor
      endif

      ## The shared 2*pi factor cancels in the normalization, and subtracting
      ## the row maximum first keeps a well separated observation a posterior
      ## rather than a ratio of two underflowed zeros.
      logscore = logscore - max (logscore, [], 2);
      score = exp (logscore);
      score = score ./ sum (score, 2);
      score(isnan (score)) = 0;

      ## Calculate expected classification cost
      for i = 1:numClasses
        cost(:, i) = sum (bsxfun (@times, score, this.Cost(:, i)'), 2);
      endfor

      ## Predict the class labels based on the minimum cost
      [~, minIdx] = min (cost, [], 2);
      label = this.ClassNames(minIdx,:);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationDiscriminant} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute loss for a trained ClassificationDiscriminant object.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} computes the loss,
    ## @var{L}, using the default loss function @qcode{'mincost'}.
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{ClassificationDiscriminant} object trained on
    ## @code{X} and @code{Y}.
    ## @item
    ## @code{X} must be a @math{N*P} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of Rows as @var{X}.
    ## @end itemize
    ##
    ## @code{@var{L} = loss (@dots{}, @var{name}, @var{value})} allows
    ## additional options specified by @var{name}-@var{value} pairs:
    ##
    ## @multitable @columnfractions 0.18 0.8
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'LossFun'} @tab Specifies the loss function to use.
    ## Can be a function handle with four input arguments (C, S, W, Cost)
    ## which returns a scalar value or one of:
    ## 'binodeviance', 'classifcost', 'classiferror', 'exponential',
    ## 'hinge', 'logit','mincost', 'quadratic'.
    ## @itemize
    ## @item
    ## @code{C} is a logical matrix of size @math{N*K}, where @math{N} is the
    ## number of observations and @math{K} is the number of classes.
    ## The element @code{C(i,j)} is true if the class label of the i-th
    ## observation is equal to the j-th class.
    ## @item
    ## @code{S} is a numeric matrix of size @math{N*K}, where each element
    ## represents the classification score for the corresponding class.
    ## @item
    ## @code{W} is a numeric vector of length @math{N}, representing
    ## the observation weights.
    ## @item
    ## @code{Cost} is a @math{K*K} matrix representing the misclassification
    ## costs.
    ## @end itemize
    ##
    ## @item @qcode{'Weights'} @tab Specifies observation weights, must be
    ## a numeric vector of length equal to the number of rows in X.
    ## Default is @code{ones (size (X, 1))}. loss normalizes the weights so that
    ## observation weights in each class sum to the prior probability of that
    ## class. When you supply Weights, loss computes the weighted
    ## classification loss.
    ##
    ## @end multitable
    ##
    ## @seealso{ClassificationDiscriminant}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationDiscriminant.loss: too few input arguments.");
      elseif (mod (nargin - 3, 2) != 0)
        error (strcat ("ClassificationDiscriminant.loss: name-value", ...
                       " arguments must be in pairs."));
      elseif (nargin > 7)
        error ("ClassificationDiscriminant.loss: too many input arguments.");
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationDiscriminant.loss: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat ("ClassificationDiscriminant.loss: X must have the", ...
                       " same number of predictors as the trained model."));
      endif

      ## Default values
      LossFun = 'mincost';
      Weights = [];

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error ("ClassificationDiscriminant.loss: Y must be of a valid type.");
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat ("ClassificationDiscriminant.loss: Y must", ...
                       " have the same number of rows as X."));
      endif

      ## Parse name-value arguments
      while (numel (varargin) > 0)
        Value = varargin{2};
        switch (tolower (varargin{1}))
          case 'lossfun'
            lf_opt = {'binodeviance', 'classifcost', 'classiferror', ...
                      'exponential', 'hinge','logit', 'mincost', 'quadratic'};
            if (isa (Value, 'function_handle'))
              ## Check if the loss function is valid
              if (nargin (Value) != 4)
                error (strcat ("ClassificationDiscriminant.loss: custom", ...
                               " loss function must accept exactly four", ...
                               " input arguments."));
              endif
              try
                n = 1;
                K = 2;
                C_test = false (n, K);
                S_test = zeros (n, K);
                W_test = ones (n, 1);
                Cost_test = ones (K) - eye (K);
                test_output = Value(C_test, S_test, W_test, Cost_test);
                if (! isscalar (test_output))
                  error (strcat ("ClassificationDiscriminant.loss:", ...
                                 " custom loss function must return", ...
                                 " a scalar value."));
                endif
              catch
                error (strcat ("ClassificationDiscriminant.loss: custom", ...
                               " loss function is not valid or does not", ...
                               " produce correct output."));
              end_try_catch
              LossFun = Value;
            elseif (ischar (Value) && any (strcmpi (Value, lf_opt)))
              LossFun = Value;
            else
              error ("ClassificationDiscriminant.loss: invalid loss function.");
            endif

          case 'weights'
            if (isnumeric (Value) && isvector (Value))
              if (numel (Value) != size (X ,1))
                error (strcat ("ClassificationDiscriminant.loss: number", ...
                               " of 'Weights' must be equal to the", ...
                               " number of rows in X."));
              elseif (numel (Value) == size (X, 1))
                Weights = Value;
              endif
            else
              error ("ClassificationDiscriminant.loss: invalid 'Weights'.");
            endif

          otherwise
            error (strcat ("ClassificationDiscriminant.loss: invalid", ...
                           " parameter name in optional pair arguments."));
        endswitch
        varargin(1:2) = [];
      endwhile

      ## Check for missing values in X
      if (! isa (LossFun, 'function_handle'))
        lossfun = tolower (LossFun);
        if (! strcmp (lossfun, 'mincost') && ! strcmp (lossfun, 'classiferror')
            && ! strcmp (lossfun, 'classifcost') && any (isnan (X(:))))
          L = NaN;
          return;
        endif
      endif

      ## If Y is a char array convert it to a cell array of character vectors
      classes = this.ClassNames;
      if (ischar (Y) && ischar (classes))
        Y = cellstr (Y);
        classes = cellstr (classes);
      endif

      ## Check that Y is of the same type as ClassNames
      if (! strcmp (class (Y), class (classes)))
        error (strcat ("ClassificationDiscriminant.loss: Y must be", ...
                       " the same data type as the model's ClassNames."));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), classes)))
        error (strcat ("ClassificationDiscriminant.loss: Y must contain", ...
                       " only the classes in model's ClassNames."));
      endif

      ## Set default weights if not specified
      if (isempty (Weights))
        Weights = ones (size (X, 1), 1);
      endif

      ## Normalize Weights
      K = numel (classes);
      class_prior_probs = this.Prior;
      norm_weights = zeros (size (Weights));
      for i = 1:K
        class_idx = ismember (Y, classes(i));
        if (sum (Weights(class_idx)) > 0)
          norm_weights(class_idx) = ...
          Weights(class_idx) * class_prior_probs(i) / sum (Weights(class_idx));
        endif
      endfor
      Weights = norm_weights / sum (norm_weights);

      ## Number of observations
      n = size (X, 1);

      ## Predict classification scores
      [label, scores] = predict (this, X);

      ## C is vector of K-1 zeros, with 1 in the
      ## position corresponding to the true class
      C = false (n, K);
      for i = 1:n
        class_idx = find (ismember (classes, Y(i)));
        C(i, class_idx) = true;
      endfor
      Y_new = C';

      ## Compute the loss using custom loss function
      if (isa (LossFun, 'function_handle'))
        L = LossFun(C, scores, Weights, this.Cost);
        return;
      endif

      ## Compute the scalar classification score for each observation
      m_j = zeros (n, 1);
      for i = 1:n
        m_j(i) = scores(i,:) * Y_new(:,i);
      endfor

      ## Compute the loss
      switch (tolower (LossFun))
        case 'binodeviance'
          b = log (1 + exp (-2 * m_j));
          L = (Weights') * b;
        case 'hinge'
          h = max (0, 1 - m_j);
          L = (Weights') * h;
        case 'exponential'
          e = exp (-m_j);
          L = (Weights') * e;
        case 'logit'
          l = log (1 + exp (-m_j));
          L = (Weights') * l;
        case 'quadratic'
          q = (1 - m_j) .^ 2;
          L = (Weights') * q;
        case 'classiferror'
          L = 0;
          for i = 1:n
            L = L + Weights(i) * (! isequal (Y(i), label(i)));
          endfor
        case 'mincost'
          Cost = this.Cost;
          L = 0;
          for i = 1:n
            f_Xj = scores(i, :);
            gamma_jk = f_Xj * Cost;
            [~, min_cost_class] = min (gamma_jk);
            cj = Cost(find (ismember (classes, Y(i))), min_cost_class);
            L = L + Weights(i) * cj;
          endfor
        case 'classifcost'
          Cost = this.Cost;
          L = 0;
          for i = 1:n
            y_idx = find (ismember (classes, Y(i)));
            y_hat_idx = find (ismember (classes, label(i)));
            L = L + Weights(i) * Cost(y_idx, y_hat_idx);
          endfor
        otherwise
          error ("ClassificationDiscriminant.loss: invalid loss function.");
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationDiscriminant} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margins for discriminant analysis classifier.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns
    ## the classification margins for @var{obj} with data @var{X} and
    ## classification @var{Y}. @var{m} is a numeric vector of length size (X,1).
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{ClassificationDiscriminant} object trained on
    ## @code{X}
    ## and @code{Y}.
    ## @item
    ## @code{X} must be a @math{N*P} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of Rows as @var{X}.
    ## @end itemize
    ##
    ## The classification margin for each observation is the difference between
    ## the classification score for the true class and the maximal
    ## classification score for the false classes.
    ##
    ## @seealso{fitcdiscr, ClassificationDiscriminant}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error ("ClassificationDiscriminant.margin: too few input arguments.");
      endif

      ## Check for valid X
      if (isempty (X))
        error ("ClassificationDiscriminant.margin: X is empty.");
      elseif (columns (this.X) != columns (X))
        error (strcat ("ClassificationDiscriminant.margin: X must have the", ...
                       " same number of predictors as the trained model."));
      endif

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error ("ClassificationDiscriminant.margin: Y must be of a valid type.");
      endif

      ## Validate X
      valid_types = {'single', 'double'};
      if (! (any (strcmp (class (X), valid_types))))
        error ("ClassificationDiscriminant.margin: X must be of a valid type.");
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat ("ClassificationDiscriminant.margin: Y must", ...
                       " have the same number of rows as X."));
      endif

      ## If Y is a char array convert it to a cell array of character vectors
      classes = this.ClassNames;
      if (ischar (Y) && ischar (classes))
        Y = cellstr (Y);
        classes = cellstr (classes);
      endif

      ## Check that Y is of the same type as ClassNames
      if (! strcmp (class (Y), class (classes)))
        error (strcat ("ClassificationDiscriminant.margin: Y must be", ...
                       " the same data type as the model's ClassNames."));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), classes)))
        error (strcat ("ClassificationDiscriminant.margin: Y must", ...
                       " contain only the classes in model's ClassNames."));
      endif

      ## Number of Observations
      n = size (X, 1);

      ## Initialize the margin vector
      m = zeros (n, 1);

      ## Calculate the classification scores
      [~, scores] = predict (this, X);

      ## Loop over each observation to compute the margin
      for i = 1:n
        ## True class index
        true_class_idx = find (ismember (classes, Y(i)));

        ## Score for the true class
        true_class_score = scores(i, true_class_idx);

        ## Get the maximal score for the false classes
        scores(i, true_class_idx) = -Inf;              # Temporarily
        max_false_class_score = max (scores(i, :));
        if (max_false_class_score == -Inf)
          m = NaN;
          return;
        endif
        scores(i, true_class_idx) = true_class_score;  # Restore

        ## Calculate the margin
        m(i) = true_class_score - max_false_class_score;
      endfor

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{CVMdl} =} crossval (@var{obj})
    ## @deftypefnx {ClassificationDiscriminant} {@var{CVMdl} =} crossval (@dots{}, @var{Name}, @var{Value})
    ##
    ## Cross Validate a Discriminant classification object.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj})} returns a cross-validated model
    ## object, @var{CVMdl}, from a trained model, @var{obj}, using 10-fold
    ## cross-validation by default.
    ##
    ## @code{@var{CVMdl} = crossval (@var{obj}, @var{name}, @var{value})}
    ## specifies additional name-value pair arguments to customize the
    ## cross-validation process.
    ##
    ## @multitable @columnfractions 0.28 0.7
    ## @headitem @var{Name} @tab @var{Value}
    ##
    ## @item @qcode{'KFold'} @tab Specify the number of folds to use in
    ## k-fold cross-validation.  @code{"KFold", @var{k}}, where @var{k} is an
    ## integer greater than 1.
    ##
    ## @item @qcode{'Holdout'} @tab Specify the fraction of the data to
    ## hold out for testing.  @code{"Holdout", @var{p}}, where @var{p} is a
    ## scalar in the range @math{(0,1)}.
    ##
    ## @item @qcode{'Leaveout'} @tab Specify whether to perform
    ## leave-one-out cross-validation.  @code{"Leaveout", @var{Value}}, where
    ## @var{Value} is 'on' or 'off'.
    ##
    ## @item @qcode{'CVPartition'} @tab Specify a @qcode{cvpartition}
    ## object used for cross-validation.  @code{"CVPartition", @var{cv}}, where
    ## @code{isa (@var{cv}, "cvpartition")} = 1.
    ##
    ## @end multitable
    ##
    ## @seealso{fitcdiscr, ClassificationDiscriminant, cvpartition,
    ## ClassificationPartitionedModel}
    ## @end deftypefn
    function CVMdl = crossval (this, varargin)
      ## Check input
      if (nargin < 1)
        error ("ClassificationDiscriminant.crossval: too few input arguments.");
      endif

      if (numel (varargin) == 1)
        error (strcat ("ClassificationDiscriminant.crossval: Name-Value", ...
                       " arguments must be in pairs."));
      elseif (numel (varargin) > 2)
        error (strcat ("ClassificationDiscriminant.crossval: specify only", ...
                       " one of the optional Name-Value paired arguments."));
      endif

      ## Add default values
      if (this.NumObservations < 10)
        numFolds  = this.NumObservations;
      else
        numFolds  = 10;
      endif
      Holdout     = [];
      Leaveout    = 'off';
      CVPartition = [];

      ## Parse extra parameters
      while (numel (varargin) > 0)
        switch (tolower (varargin {1}))

          case 'kfold'
            numFolds = varargin{2};
            if (! (isnumeric (numFolds) && isscalar (numFolds)
                   && (numFolds == fix (numFolds)) && numFolds > 1))
              error (strcat ("ClassificationDiscriminant.crossval: 'KFold'", ...
                             " must be an integer value greater than 1."));
            endif

          case 'holdout'
            Holdout = varargin{2};
            if (! (isnumeric (Holdout) && isscalar (Holdout) && Holdout > 0
                   && Holdout < 1))
              error (strcat ("ClassificationDiscriminant.crossval: 'Holdout'", ...
                             " must be a numeric value between 0 and 1."));
            endif

          case 'leaveout'
            Leaveout = varargin{2};
            if (! (ischar (Leaveout)
                   && (strcmpi (Leaveout, 'on') || strcmpi (Leaveout, 'off'))))
              error (strcat ("ClassificationDiscriminant.crossval:", ...
                             " 'Leaveout' must be either 'on' or 'off'."));
            endif

          case 'cvpartition'
            CVPartition = varargin{2};
            if (! (isa (CVPartition, 'cvpartition')))
              error (strcat ("ClassificationDiscriminant.crossval:",...
                             " 'CVPartition' must be a 'cvpartition' object."));
            endif

          otherwise
            error (strcat ("ClassificationDiscriminant.crossval: invalid",...
                           " parameter name in optional paired arguments."));
          endswitch
        varargin(1:2) = [];
      endwhile

      ## Determine the cross-validation method to use.  The partition covers
      ## the observations actually trained on: a row dropped for a missing
      ## value is not one the folds can use, and including it would leave the
      ## partition, the stored data and NumObservations disagreeing.  The
      ## response is passed rather than a count so the folds stay stratified.
      Yused = this.Y;
      if (! isempty (CVPartition))
        partition = CVPartition;
      elseif (! isempty (Holdout))
        partition = cvpartition (Yused, 'Holdout', Holdout);
      elseif (strcmpi (Leaveout, 'on'))
        partition = cvpartition (numel (this.Y), 'LeaveOut');
      else
        partition = cvpartition (Yused, 'KFold', numFolds);
      endif

      ## Create a cross-validated model object
      CVMdl = ClassificationPartitionedModel (this, partition);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{CVMdl} =} compact (@var{obj})
    ##
    ## Create a CompactClassificationDiscriminant object.
    ##
    ## @code{@var{CVMdl} = compact (@var{obj})} creates a compact version of the
    ## ClassificationDiscriminant object, @var{obj}.
    ##
    ## @seealso{fitcdiscr, ClassificationDiscriminant,
    ## CompactClassificationDiscriminant}
    ## @end deftypefn
    function CVMdl = compact (this)
      ## Create a compact model
      CVMdl = CompactClassificationDiscriminant (this);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {ClassificationDiscriminant} {@var{e} =} edge (@dots{}, @qcode{"Weights"}, @var{w})
    ##
    ## Classification edge, the mean of the classification margins.
    ##
    ## @code{@var{e} = edge (@var{obj}, @var{X}, @var{Y})} reduces the vector
    ## that @code{margin} returns to a single number, the mean margin over the
    ## rows of @var{X}.  It says how far the model puts the true class ahead of
    ## its nearest rival on average, so a larger edge is a better model, and
    ## unlike a loss it is not bounded above and rewards confidence rather than
    ## bare correctness.
    ##
    ## @code{@var{e} = edge (@dots{}, @qcode{"Weights"}, @var{w})} takes the
    ## weighted mean instead, with one weight per row of @var{X}.
    ##
    ## @end deftypefn
    function e = edge (this, X, Y, varargin)

      if (nargin < 3)
        error ("ClassificationDiscriminant.edge: too few input arguments.");
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("ClassificationDiscriminant.edge: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      ## The weights are parsed before anything is computed, so a bad
      ## Name-Value pair is reported as such rather than after a margin.
      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       "ClassificationDiscriminant", "edge");
      m = margin (this, X, Y);
      e = sum (W .* m(:)) / sum (W);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{label} =} resubPredict (@var{obj})
    ## @deftypefnx {ClassificationDiscriminant} {[@var{label}, @var{score}, @var{cost}] =} resubPredict (@var{obj})
    ##
    ## Classify the training data with the model fitted to it.
    ##
    ## @code{@var{label} = resubPredict (@var{obj})} is @code{predict} applied
    ## to the observations the model was fitted on, which it holds in
    ## @qcode{X}.  Handing them over yourself is not the same thing: a row
    ## dropped for a missing response is not in @qcode{X}, so the original
    ## matrix and the model's own are different data.
    ##
    ## The result measures fit and not generalization, and is optimistic by
    ## construction.  @code{crossval} is what estimates performance on data the
    ## model has not seen.
    ##
    ## @end deftypefn
    function [label, score, cost] = resubPredict (this)
      [label, score, cost] = predict (this, this.X);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationDiscriminant} {@var{m} =} resubMargin (@var{obj})
    ##
    ## Classification margins of the model on its own training data.
    ##
    ## @code{@var{m} = resubMargin (@var{obj})} is @code{margin} applied to the
    ## observations the model was fitted on, one number per observation.  Being
    ## a resubstitution quantity it is optimistic by construction.
    ##
    ## @end deftypefn
    function m = resubMargin (this)
      m = margin (this, this.X, this.Y);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {ClassificationDiscriminant} {@var{e} =} resubEdge (@var{obj})
    ##
    ## Classification edge of the model on its own training data.
    ##
    ## @code{@var{e} = resubEdge (@var{obj})} is @code{edge} applied to the
    ## observations the model was fitted on, the mean of @code{resubMargin}.
    ##
    ## @end deftypefn
    function e = resubEdge (this)
      e = edge (this, this.X, this.Y);
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {@var{L} =} resubLoss (@var{obj})
    ## @deftypefnx {ClassificationDiscriminant} {@var{L} =} resubLoss (@dots{}, @var{name}, @var{value})
    ##
    ## Classification loss of the model on its own training data.
    ##
    ## @code{@var{L} = resubLoss (@var{obj})} is @code{loss} applied to the
    ## observations the model was fitted on, defaulting to
    ## @qcode{'mincost'}, and it accepts the same @qcode{Name-Value} pairs.
    ##
    ## Being a resubstitution quantity it is a lower bound on the error rather
    ## than an estimate of it.  It is worth least on a lazy learner: a
    ## one-neighbour @code{ClassificationKNN} has a resubstitution loss of
    ## exactly zero, every training point being its own nearest neighbour.
    ##
    ## @end deftypefn
    function L = resubLoss (this, varargin)
      L = loss (this, this.X, this.Y, varargin{:});
    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {ClassificationDiscriminant} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a ClassificationDiscriminant object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## ClassificationDiscriminant object into an Octave binary file, the name of
    ## which is specified in @var{filename}, along with an extra variable, which
    ## defines the type classification object these variables constitute.  Use
    ## @code{loadmodel} in order to load a classification object into Octave's
    ## workspace.
    ##
    ## @seealso{loadmodel, fitcdiscr, ClassificationDiscriminant}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error ("ClassificationDiscriminant.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error ("ClassificationDiscriminant.savemodel: FNAME must be a character vector.");
      endif
      ## Generate variable for class name
      classdef_name = 'ClassificationDiscriminant';

      ## Create variables from model properties
      X = this.X;
      Y = this.Y;
      NumObservations = this.NumObservations;
      W               = this.W;
      RowsUsed        = this.RowsUsed;
      NumPredictors   = this.NumPredictors;
      PredictorNames  = this.PredictorNames;
      ResponseName    = this.ResponseName;
      ClassNames      = this.ClassNames;
      Prior           = this.Prior;
      Cost            = this.Cost;
      ScoreTransform  = this.ScoreTransform;
      Sigma           = this.Sigma;
      BinEdges        = this.BinEdges;
      BaseSigma       = this.BaseSigma;
      Mu              = this.Mu;
      Coeffs          = this.Coeffs;
      Delta           = this.Delta;
      DiscrimType     = this.DiscrimType;
      Gamma           = this.Gamma;
      MinGamma        = this.MinGamma;
      LogDetSigma     = this.LogDetSigma;
      XCentered       = this.XCentered;
      BetweenSigma    = this.BetweenSigma;
      CategoricalPredictors  = this.CategoricalPredictors;
      ExpandedPredictorNames = this.ExpandedPredictorNames;
      STfun          = this.STfun;

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'X', 'Y', 'NumObservations', ...
            'W', 'RowsUsed', 'NumPredictors', 'PredictorNames', 'BinEdges', ...
            'ResponseName', ...
            'ClassNames', 'ScoreTransform', 'Prior', 'Cost', 'Sigma', ...
            'BaseSigma', 'Mu', ...
            'Coeffs', 'Delta', 'DiscrimType', 'Gamma', 'MinGamma', ...
            'LogDetSigma', 'XCentered', 'BetweenSigma', ...
            'CategoricalPredictors', 'ExpandedPredictorNames', 'STfun');
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a ClassificationDiscriminant object
      ## Built without coefficients: every set method rebuilds Coeffs when it
      ## finds one, and a stub's would be rebuilt against half a loaded model.
      mdl = ClassificationDiscriminant (1, 1, 'FillCoeffs', 'off');

      ## Copy the saved data into the object.  Iterate over what was
      ## saved rather than over fieldnames (mdl): a private property such
      ## as STfun is written out by savemodel but is not reported by
      ## fieldnames, so comparing the two sets could never match and every
      ## load failed.  Assignment is legal here because this is a method of
      ## the class itself.
      names = fieldnames (data);
      ## The set methods for these read other properties, and some of them
      ## rebuild Sigma or Coeffs, so they are assigned once everything else is
      ## in place rather than in the order the file happens to list them.
      ## Order matters within the late group as well: Prior and Cost before
      ## DiscrimType, because assigning the type rebuilds Coeffs and Coeffs
      ## reads the priors; and DiscrimType before Gamma, because assigning the
      ## type implies a Gamma and would overwrite the saved one if it came
      ## second.  Coeffs is installed last of all: the set methods rebuild it
      ## whenever it is already there, and rebuilding it from a half-loaded
      ## model is how a diagonal Sigma met a linear formula.
      order = {'Cost', 'Prior', 'ScoreTransform', 'ResponseTransform', ...
               'DiscrimType', 'Gamma', 'Delta', 'Coeffs'};
      late = ismember (names, order);
      tail = {};
      for k = 1:numel (order)
        if (any (strcmp (names, order{k})))
          tail{end+1} = order{k};
        endif
      endfor
      names = [names(! late); tail(:)];
      for i = 1:numel (names)
        try
          mdl.(names{i}) = data.(names{i});
        catch
          msg = 'ClassificationDiscriminant.load_model: invalid model in ''%s''.';
          error (msg, filename);
        end_try_catch
      endfor
    endfunction

  endmethods

  methods(Access = private)



  endmethods

endclassdef

%!demo
%! ## Create discriminant classifier
%! ## Evaluate some model predictions on new data.
%!
%! load fisheriris
%! x = meas;
%! y = species;
%! xc = [min(x); mean(x); max(x)];
%! obj = fitcdiscr (x, y);
%! [label, score, cost] = predict (obj, xc);

%!demo
%! load fisheriris
%! model = fitcdiscr (meas, species);
%! X = mean (meas);
%! Y = {'versicolor'};
%! ## Compute loss for discriminant model
%! L = loss (model, X, Y)

%!demo
%! load fisheriris
%! mdl = fitcdiscr (meas, species);
%! X = mean (meas);
%! Y = {'versicolor'};
%! ## Margin for discriminant model
%! m = margin (mdl, X, Y)

%!demo
%! load fisheriris
%! x = meas;
%! y = species;
%! obj = fitcdiscr (x, y, 'gamma', 0.4);
%! ## Cross-validation for discriminant model
%! CVMdl = crossval (obj)

## Test constructor
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! PredictorNames = {'Sepal Length', 'Sepal Width', 'Petal Length', 'Petal Width'};
%! Mdl = ClassificationDiscriminant (x, y, 'PredictorNames', PredictorNames);
%! sigma = [0.265008, 0.092721, 0.167514, 0.038401; ...
%!          0.092721, 0.115388, 0.055244, 0.032710; ...
%!          0.167514, 0.055244, 0.185188, 0.042665; ...
%!          0.038401, 0.032710, 0.042665, 0.041882];
%! mu = [5.0060, 3.4280, 1.4620, 0.2460; ...
%!       5.9360, 2.7700, 4.2600, 1.3260; ...
%!       6.5880, 2.9740, 5.5520, 2.0260];
%! xCentered = [ 9.4000e-02,  7.2000e-02, -6.2000e-02, -4.6000e-02; ...
%!              -1.0600e-01, -4.2800e-01, -6.2000e-02, -4.6000e-02; ...
%!              -3.0600e-01, -2.2800e-01, -1.6200e-01, -4.6000e-02];
%! assert_equal (class (Mdl), "ClassificationDiscriminant");
%! assert_equal ({Mdl.X, Mdl.Y, Mdl.NumObservations}, {x, y, 150})
%! assert_equal ({Mdl.DiscrimType, Mdl.ResponseName}, {'linear', 'Y'})
%! assert_equal ({Mdl.Gamma, Mdl.MinGamma}, {0, 0}, 1e-15)
%! assert_equal (Mdl.ClassNames, unique (species))
%! assert_equal (Mdl.Sigma, sigma, 1e-6)
%! assert_equal (Mdl.Mu, mu, 1e-14)
%! assert_equal (Mdl.XCentered([1:3],:), xCentered, 1e-14)
%! assert_equal (Mdl.LogDetSigma, -9.9585, 1e-4)
%! assert_equal (Mdl.PredictorNames, PredictorNames)
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! Mdl = ClassificationDiscriminant (x, y, 'Gamma', 0.5);
%! sigma = [0.265008, 0.046361, 0.083757, 0.019201; ...
%!          0.046361, 0.115388, 0.027622, 0.016355; ...
%!          0.083757, 0.027622, 0.185188, 0.021333; ...
%!          0.019201, 0.016355, 0.021333, 0.041882];
%! mu = [5.0060, 3.4280, 1.4620, 0.2460; ...
%!       5.9360, 2.7700, 4.2600, 1.3260; ...
%!       6.5880, 2.9740, 5.5520, 2.0260];
%! xCentered = [ 9.4000e-02,  7.2000e-02, -6.2000e-02, -4.6000e-02; ...
%!              -1.0600e-01, -4.2800e-01, -6.2000e-02, -4.6000e-02; ...
%!              -3.0600e-01, -2.2800e-01, -1.6200e-01, -4.6000e-02];
%! assert_equal (class (Mdl), "ClassificationDiscriminant");
%! assert_equal ({Mdl.X, Mdl.Y, Mdl.NumObservations}, {x, y, 150})
%! assert_equal ({Mdl.DiscrimType, Mdl.ResponseName}, {'linear', 'Y'})
%! assert_equal ({Mdl.Gamma, Mdl.MinGamma}, {0.5, 0})
%! assert_equal (Mdl.ClassNames, unique (species))
%! assert_equal (Mdl.Sigma, sigma, 1e-6)
%! assert_equal (Mdl.Mu, mu, 1e-14)
%! assert_equal (Mdl.XCentered([1:3],:), xCentered, 1e-14)
%! assert_equal (Mdl.LogDetSigma, -8.6884, 1e-4)

## Test input validation for constructor
%!shared X, Y, MODEL
%! X = rand (10,2);
%! Y = [ones(5,1);2*ones(5,1)];
%! MODEL = ClassificationDiscriminant (X, Y);
%!error<ClassificationDiscriminant: too few input arguments.> ClassificationDiscriminant ()
%!error<ClassificationDiscriminant: too few input arguments.> ...
%! ClassificationDiscriminant (ones (4, 1))
%!error<ClassificationDiscriminant: Name-Value arguments must be in pairs.> ...
%! ClassificationDiscriminant (X, Y, 'prior')
%!error<ClassificationDiscriminant: number of rows in X and Y must be equal.> ...
%! ClassificationDiscriminant (ones (4,2), ones (1,4))
%!error<ClassificationDiscriminant: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationDiscriminant (X, Y, 'PredictorNames', ['A'])
%!error<ClassificationDiscriminant: 'PredictorNames' must be supplied as a cellstring array.> ...
%! ClassificationDiscriminant (X, Y, 'PredictorNames', 'A')
%!error<ClassificationDiscriminant: 'PredictorNames' must equal the number of columns in X.> ...
%! ClassificationDiscriminant (X, Y, 'PredictorNames', {'A', 'B', 'C'})
%!error<ClassificationDiscriminant: 'ResponseName' must be a character vector.> ...
%! ClassificationDiscriminant (X, Y, 'ResponseName', {'Y'})
%!error<ClassificationDiscriminant: 'ResponseName' must be a character vector.> ...
%! ClassificationDiscriminant (X, Y, 'ResponseName', 1)
%!error<ClassificationDiscriminant: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationDiscriminant (X, Y, 'ClassNames', @(x)x)
%!error<ClassificationDiscriminant: 'ClassNames' must be a cell array of character vectors, a logical vector, a numeric vector, or a character array.> ...
%! ClassificationDiscriminant (X, Y, 'ClassNames', {1})
%!error<ClassificationDiscriminant: not all 'ClassNames' are present in Y.> ...
%! ClassificationDiscriminant (X, ones (10,1), 'ClassNames', [1, 2])
%!error<ClassificationDiscriminant: not all 'ClassNames' are present in Y.> ...
%! ClassificationDiscriminant ([1;2;3;4;5], ['a';'b';'a';'a';'b'], 'ClassNames', ['a';'c'])
%!error<ClassificationDiscriminant: not all 'ClassNames' are present in Y.> ...
%! ClassificationDiscriminant ([1;2;3;4;5], {'a';'b';'a';'a';'b'}, 'ClassNames', {'a','c'})
%!error<ClassificationDiscriminant: not all 'ClassNames' are present in Y.> ...
%! ClassificationDiscriminant (X, logical (ones (10,1)), 'ClassNames', [true, false])
%!error<ClassificationDiscriminant: 'Prior' must be either a numeric or a character vector.> ...
%! ClassificationDiscriminant (X, Y, 'Prior', {'1', '2'})
%!error<ClassificationDiscriminant: the elements in 'Prior' do not correspond to the selected classes in Y.> ...
%! ClassificationDiscriminant (X, ones (10,1), 'Prior', [1 2])
%!error<ClassificationDiscriminant: 'Cost' must be a numeric square matrix.> ...
%! ClassificationDiscriminant (X, Y, 'Cost', [1, 2])
%!error<ClassificationDiscriminant: 'Cost' must be a numeric square matrix.> ...
%! ClassificationDiscriminant (X, Y, 'Cost', 'string')
%!error<ClassificationDiscriminant: 'Cost' must be a numeric square matrix.> ...
%! ClassificationDiscriminant (X, Y, 'Cost', {eye(2)})
%!error<ClassificationDiscriminant: the number of rows and columns in 'Cost' must correspond to selected classes in Y.> ...
%! ClassificationDiscriminant (X, Y, 'Cost', ones (3))
%!error<ClassificationDiscriminant: predictor 'x1' has zero within-class variance.  Either exclude this predictor or set 'DiscrimType' to 'pseudoLinear' or 'diagLinear'.> ...
%! ClassificationDiscriminant (ones (5,2), [1; 1; 2; 2; 2])
%!error<ClassificationDiscriminant: predictor 'A' has zero within-class variance.  Either exclude this predictor or set 'DiscrimType' to 'pseudoLinear' or 'diagLinear'.> ...
%! ClassificationDiscriminant (ones (5,2), [1; 1; 2; 2; 2], 'PredictorNames', {'A', 'B'})
%!error<ClassificationDiscriminant: predictor 'x2' has zero within-class variance.  Either exclude this predictor or set 'DiscrimType' to 'pseudoLinear' or 'diagLinear'.> ...
%! ClassificationDiscriminant ([1,2;2,2;3,2;4,2;5,2], ones (5, 1))
%!error<ClassificationDiscriminant: predictor 'B' has zero within-class variance.  Either exclude this predictor or set 'DiscrimType' to 'pseudoLinear' or 'diagLinear'.> ...
%! ClassificationDiscriminant ([1,2;2,2;3,2;4,2;5,2], ones (5, 1), 'PredictorNames', {'A', 'B'})

## Test predict method
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! Mdl = fitcdiscr (meas, species, 'Gamma', 0.5);
%! [label, score, cost] = predict (Mdl, [2, 2, 2, 2]);
%! assert_equal (label, {'versicolor'})
%! assert_equal (score, [0, 0.9999, 0.0001], 1e-4)
%! assert_equal (cost, [1, 0.0001, 0.9999], 1e-4)
%! [label, score, cost] = predict (Mdl, [2.5, 2.5, 2.5, 2.5]);
%! assert_equal (label, {'versicolor'})
%! assert_equal (score, [0, 0.6368, 0.3632], 1e-4)
%! assert_equal (cost, [1, 0.3632, 0.6368], 1e-4)
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! xc = [min(x); mean(x); max(x)];
%! Mdl = fitcdiscr (x, y);
%! [label, score, cost] = predict (Mdl, xc);
%! l = {'setosa'; 'versicolor'; 'virginica'};
%! s = [1, 0, 0; 0, 1, 0; 0, 0, 1];
%! c = [0, 1, 1; 1, 0, 1; 1, 1, 0];
%! assert_equal (label, l)
%! assert_equal (score, s, 1e-4)
%! assert_equal (cost, c, 1e-4)

## Test input validation for predict method
%!error<ClassificationDiscriminant.predict: too few input arguments.> ...
%! predict (MODEL)
%!error<ClassificationDiscriminant.predict: XC is empty.> ...
%! predict (MODEL, [])
%!error<ClassificationDiscriminant.predict: XC must have the same number of predictors as the trained model.> ...
%! predict (MODEL, 1)

## Test loss method
%!test
%! load fisheriris
%! model = fitcdiscr (meas, species);
%! x = mean (meas);
%! y = {'versicolor'};
%! L = loss (model, x, y);
%! assert_equal (L, 0)
%!test
%! x = [1, 2; 3, 4; 5, 6];
%! y = {'A'; 'B'; 'A'};
%! model = fitcdiscr (x, y, 'Gamma', 0.4);
%! x_test = [1, 6; 3, 3];
%! y_test = {'A'; 'B'};
%! L = loss (model, x_test, y_test);
%! assert_equal (L, 0.3333, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = ['1'; '2'; '3'; '1'];
%! model = fitcdiscr (x, y, 'gamma' , 0.5);
%! x_test = [3, 3];
%! y_test = ['1'];
%! L = loss (model, x_test, y_test, 'LossFun', 'quadratic');
%! assert_equal (L, 0.2423, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = ['1'; '2'; '3'; '1'];
%! model = fitcdiscr (x, y, 'gamma' , 0.5);
%! x_test = [3, 3; 5, 7];
%! y_test = ['1'; '2'];
%! L = loss (model, x_test, y_test, 'LossFun', 'classifcost');
%! assert_equal (L, 0.3333, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = ['1'; '2'; '3'; '1'];
%! model = fitcdiscr (x, y, 'gamma' , 0.5);
%! x_test = [3, 3; 5, 7];
%! y_test = ['1'; '2'];
%! L = loss (model, x_test, y_test, 'LossFun', 'hinge');
%! assert_equal (L, 0.5886, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6; 7, 8];
%! y = ['1'; '2'; '3'; '1'];
%! model = fitcdiscr (x, y, 'gamma' , 0.5);
%! x_test = [3, 3; 5, 7];
%! y_test = ['1'; '2'];
%! W = [1; 2];
%! L = loss (model, x_test, y_test, 'LossFun', 'logit', 'Weights', W);
%! assert_equal (L, 0.5107, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6];
%! y = {'A'; 'B'; 'A'};
%! model = fitcdiscr (x, y, 'gamma' , 0.5);
%! x_with_nan = [1, 2; NaN, 4];
%! y_test = {'A'; 'B'};
%! L = loss (model, x_with_nan, y_test);
%! assert_equal (L, 0.3333, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6];
%! y = {'A'; 'B'; 'A'};
%! model = fitcdiscr (x, y);
%! x_with_nan = [1, 2; NaN, 4];
%! y_test = {'A'; 'B'};
%! L = loss (model, x_with_nan, y_test, 'LossFun', 'logit');
%! assert_equal (isnan (L), true)
%!test
%! x = [1, 2; 3, 4; 5, 6];
%! y = {'A'; 'B'; 'A'};
%! model = fitcdiscr (x, y);
%! customLossFun = @(C, S, W, Cost) sum (W .* sum (abs (C - S), 2));
%! L = loss (model, x, y, 'LossFun', customLossFun);
%! assert_equal (L, 0.8889, 1e-4)
%!test
%! x = [1, 2; 3, 4; 5, 6];
%! y = [1; 2; 1];
%! model = fitcdiscr (x, y);
%! L = loss (model, x, y, 'LossFun', 'classiferror');
%! assert_equal (L, 0.3333, 1e-4)

## Test input validation for loss method
%!error<ClassificationDiscriminant.loss: too few input arguments.> ...
%! loss (MODEL)
%!error<ClassificationDiscriminant.loss: too few input arguments.> ...
%! loss (MODEL, ones (4,2))
%!error<ClassificationDiscriminant.loss: X is empty.> ...
%! loss (MODEL, [], zeros (2))
%!error<ClassificationDiscriminant.loss: X must have the same number of predictors as the trained model.> ...
%! loss (MODEL, 1, zeros (2))
%!error<ClassificationDiscriminant.loss: name-value arguments must be in pairs.> ...
%! loss (MODEL, ones (4,2), ones (4,1), 'LossFun')
%!error<ClassificationDiscriminant.loss: Y must have the same number of rows as X.> ...
%! loss (MODEL, ones (4,2), ones (3,1))
%!error<ClassificationDiscriminant.loss: invalid loss function.> ...
%! loss (MODEL, ones (4,2), ones (4,1), 'LossFun', 'a')
%!error<ClassificationDiscriminant.loss: invalid 'Weights'.> ...
%! loss (MODEL, ones (4,2), ones (4,1), 'Weights', 'w')

## Test margin method
%! load fisheriris
%! mdl = fitcdiscr (meas, species);
%! X = mean (meas);
%! Y = {'versicolor'};
%! m = margin (mdl, X, Y);
%! assert_equal (m, 1, 1e-6)
%!test
%! X = [1, 2; 3, 4; 5, 6];
%! Y = [1; 2; 1];
%! mdl = fitcdiscr (X, Y, 'gamma', 0.5);
%! m = margin (mdl, X, Y);
%! assert_equal (m, [0.3333; -0.3333; 0.3333], 1e-4)

## Test input validation for margin method
%!error<ClassificationDiscriminant.margin: too few input arguments.> ...
%! margin (MODEL)
%!error<ClassificationDiscriminant.margin: too few input arguments.> ...
%! margin (MODEL, ones (4,2))
%!error<ClassificationDiscriminant.margin: X is empty.> ...
%! margin (MODEL, [], zeros (2))
%!error<ClassificationDiscriminant.margin: X must have the same number of predictors as the trained model.> ...
%! margin (MODEL, 1, zeros (2))
%!error<ClassificationDiscriminant.margin: Y must have the same number of rows as X.> ...
%! margin (MODEL, ones (4,2), ones (3,1))

## Test crossval method
%!shared x, y, obj
%! load fisheriris
%! x = meas;
%! y = species;
%! obj = fitcdiscr (x, y, 'gamma', 0.4);
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (obj);
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (CVMdl.KFold == 10, true)
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationDiscriminant")
%! assert_equal (CVMdl.CrossValidatedModel, "ClassificationDiscriminant")
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (obj, 'KFold', 3);
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (CVMdl.KFold == 3, true)
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationDiscriminant")
%! assert_equal (CVMdl.CrossValidatedModel, "ClassificationDiscriminant")
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (obj, 'HoldOut', 0.2);
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationDiscriminant")
%! assert_equal (CVMdl.CrossValidatedModel, "ClassificationDiscriminant")
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! CVMdl = crossval (obj, 'LeaveOut', 'on');
%! warning (status);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal ({CVMdl.X, CVMdl.Y}, {x, y})
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationDiscriminant")
%! assert_equal (CVMdl.CrossValidatedModel, "ClassificationDiscriminant")
%!test
%! status = warning;
%! warning ('off');
%! rand ('seed', 23);
%! partition = cvpartition (y, 'KFold', 3);
%! warning (status);
%! CVMdl = crossval (obj, 'cvPartition', partition);
%! assert_equal (class (CVMdl), "ClassificationPartitionedModel")
%! assert_equal (CVMdl.KFold == 3, true)
%! assert_equal (class (CVMdl.Trained{1}), "CompactClassificationDiscriminant")
%! assert_equal (CVMdl.CrossValidatedModel, "ClassificationDiscriminant")

## Test input validation for crossval method
%!error<ClassificationDiscriminant.crossval: Name-Value arguments must be in pairs.> ...
%! crossval (obj, 'kfold')
%!error<ClassificationDiscriminant.crossval: specify only one of the optional Name-Value paired arguments.>...
%! crossval (obj, 'kfold', 12, 'holdout', 0.2)
%!error<ClassificationDiscriminant.crossval: 'KFold' must be an integer value greater than 1.> ...
%! crossval (obj, 'kfold', 'a')
%!error<ClassificationDiscriminant.crossval: 'Holdout' must be a numeric value between 0 and 1.> ...
%! crossval (obj, 'holdout', 2)
%!error<ClassificationDiscriminant.crossval: 'Leaveout' must be either 'on' or 'off'.> ...
%! crossval (obj, 'leaveout', 1)
%!error<ClassificationDiscriminant.crossval: 'CVPartition' must be a 'cvpartition' object.> ...
%! crossval (obj, 'cvpartition', 1)
%!error <ClassificationDiscriminant.savemodel: too few input arguments.> ...
%! savemodel (ClassificationDiscriminant ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]))
%!error <ClassificationDiscriminant.savemodel: FNAME must be a character vector.> ...
%! savemodel (ClassificationDiscriminant ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), 1)
%!error <ClassificationDiscriminant.savemodel: FNAME must be a character vector.> ...
%! savemodel (ClassificationDiscriminant ([1, 2; 2, 3; 3, 4; 4, 5], [1; 1; 2; 2]), ['ab'; 'cd'])

## A ScoreTransform can be assigned, and is stored as a function handle.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! Mdl.ScoreTransform = 'symmetric';
%! assert_equal (class (Mdl.ScoreTransform), 'char');
%! assert_equal (Mdl.ScoreTransform, 'symmetric');
%! assert_equal (Mdl.STfun (0.25), -0.5);

## RowsUsed is empty when every observation was used.
%!test
%! load fisheriris
%! X = meas;
%! Y = grp2idx (species);
%! Mdl = fitcdiscr (X, Y);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (class (Mdl.RowsUsed), 'double');
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (rows (Mdl.X), 150);

## A missing response drops its observation and RowsUsed marks it.
%!test
%! load fisheriris
%! X = meas;
%! Y = grp2idx (species);
%! Y(5) = NaN;
%! Mdl = fitcdiscr (X, Y);
%! assert_equal (class (Mdl.RowsUsed), 'logical');
%! assert_equal (size (Mdl.RowsUsed), [150, 1]);
%! assert_equal (sum (Mdl.RowsUsed), 149);
%! assert_equal (Mdl.RowsUsed(5), false);
%! assert_equal (Mdl.NumObservations, 149);
%! assert_equal (rows (Mdl.X), 149);

## A missing predictor keeps its observation, so RowsUsed stays empty.
%!test
%! load fisheriris
%! X = meas;
%! X(3,2) = NaN;
%! Y = grp2idx (species);
%! Mdl = fitcdiscr (X, Y);
%! assert_equal (Mdl.RowsUsed, []);
%! assert_equal (Mdl.NumObservations, 150);
%! assert_equal (rows (Mdl.X), 150);
%! assert_equal (sum (isnan (Mdl.X(:))), 1);

## Prior is a row of class frequencies and W stays uniform: a discriminant
## applies the prior when it predicts, not to the fit.  Values from R2024a.
%!test
%! load fisheriris
%! i3 = [1:50, 51:80, 101:120];
%! Mdl = fitcdiscr (meas(i3,:), species(i3));
%! assert_equal (size (Mdl.Prior), [1, 3]);
%! assert_equal (Mdl.Prior, [0.5, 0.3, 0.2], 1e-14);
%! assert_equal (Mdl.W, ones (100, 1) / 100, 1e-14);

## A uniform prior leaves a discriminant's weights alone.
%!test
%! load fisheriris
%! i3 = [1:50, 51:80, 101:120];
%! Mdl = fitcdiscr (meas(i3,:), species(i3), 'Prior', 'uniform');
%! assert_equal (Mdl.Prior, [1, 1, 1] / 3, 1e-14);
%! assert_equal (Mdl.W, ones (100, 1) / 100, 1e-14);

## A structure Prior assigns by class name, whatever order it names them in.
%!test
%! load fisheriris
%! i3 = [1:50, 51:80, 101:120];
%! p = struct ('ClassNames', {{'setosa'; 'virginica'; 'versicolor'}}, ...
%!             'ClassProbs', [0.2, 0.3, 0.5]);
%! Mdl = fitcdiscr (meas(i3,:), species(i3), 'Prior', p);
%! assert_equal (Mdl.Prior, [0.2, 0.5, 0.3], 1e-14);

## An unnormalized Prior is rescaled on assignment.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! Mdl.Prior = [2, 3, 5];
%! assert_equal (Mdl.Prior, [0.2, 0.3, 0.5], 1e-14);
%! Mdl.Cost = [0, 2, 3; 1, 0, 1; 1, 1, 0];
%! assert_equal (Mdl.Cost, [0, 2, 3; 1, 0, 1; 1, 1, 0]);

## A fitted model survives savemodel and loadmodel: the properties come
## back as they were and it predicts the same.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'ClassificationDiscriminant');
%! assert_equal (M2.NumObservations, Mdl.NumObservations);
%! assert_equal (M2.PredictorNames, Mdl.PredictorNames);
%! assert_equal (class (M2.ScoreTransform), class (Mdl.ScoreTransform));
%! assert_equal (predict (M2, meas(1:5,:)), predict (Mdl, meas(1:5,:)));

## The six discriminant types.  Every value below was measured on MATLAB
## R2024a; see DISCRIMINANT_LEDGER.md for the probes that produced them.

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'linear');
%! assert_equal (Mdl.DiscrimType, 'linear');
%! assert_equal (size (Mdl.Sigma), [4, 4]);
%! assert_equal (Mdl.LogDetSigma, -9.9585, 1e-4);
%! assert_equal (Mdl.Gamma, 0);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'quadratic');
%! assert_equal (Mdl.DiscrimType, 'quadratic');
%! assert_equal (size (Mdl.Sigma), [4, 4, 3]);
%! assert_equal (Mdl.LogDetSigma', [-13.0674, -10.8743, -8.9271], 1e-4);
%! assert_equal (Mdl.Gamma, 0);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'diagLinear');
%! assert_equal (Mdl.DiscrimType, 'diagLinear');
%! assert_equal (size (Mdl.Sigma), [1, 4]);
%! assert_equal (Mdl.LogDetSigma, -8.3467, 1e-4);
%! assert_equal (Mdl.Gamma, 1);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'diagQuadratic');
%! assert_equal (Mdl.DiscrimType, 'diagQuadratic');
%! assert_equal (size (Mdl.Sigma), [1, 4, 3]);
%! assert_equal (Mdl.LogDetSigma', [-12.0271, -8.3925, -6.9421], 1e-4);
%! assert_equal (Mdl.Gamma, 1);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'pseudoLinear');
%! assert_equal (Mdl.DiscrimType, 'pseudoLinear');
%! assert_equal (size (Mdl.Sigma), [4, 4]);
%! assert_equal (Mdl.LogDetSigma, -9.9585, 1e-4);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'pseudoQuadratic');
%! assert_equal (Mdl.DiscrimType, 'pseudoQuadratic');
%! assert_equal (size (Mdl.Sigma), [4, 4, 3]);
%! assert_equal (Mdl.LogDetSigma', [-13.0674, -10.8743, -8.9271], 1e-4);

## The posterior of the one iris the model is unsure about separates the six
## types from one another, which a correctly classified row does not.
%!test
%! load fisheriris
%! [~, s] = predict (fitcdiscr (meas, species), meas(71,:));
%! assert_equal (s, [0, 0.2532, 0.7468], 1e-4);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'quadratic');
%! [~, s] = predict (Mdl, meas(71,:));
%! assert_equal (s, [0, 0.3359, 0.6641], 1e-4);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'diagLinear');
%! [~, s] = predict (Mdl, meas(71,:));
%! assert_equal (s, [0, 0.2646, 0.7354], 1e-4);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'diagQuadratic');
%! [~, s] = predict (Mdl, meas(71,:));
%! assert_equal (s, [0, 0.1609, 0.8391], 1e-4);

## Exact coefficients on a fixture small enough that every figure is closed
## form: the pooled covariance is 13.75/6 and each class covariance n_k-1.
%!test
%! x = [1, 2; 2, 1; 3, 4; 4, 3; 5, 7; 7, 5; 8, 9; 9, 8];
%! y = [1; 1; 1; 1; 2; 2; 2; 2];
%! Mdl = fitcdiscr (x, y, 'DiscrimType', 'linear');
%! assert_equal (Mdl.Sigma, [2.2917, 1.125; 1.125, 2.2917], 1e-4);
%! assert_equal (Mdl.LogDetSigma, 1.3828, 1e-4);
%! assert_equal (Mdl.Coeffs(1,2).Const, 13.5549, 1e-4);
%! assert_equal (Mdl.Coeffs(1,2).Linear', [-1.3902, -1.3902], 1e-4);

%!test
%! x = [1, 2; 2, 1; 3, 4; 4, 3; 5, 7; 7, 5; 8, 9; 9, 8];
%! y = [1; 1; 1; 1; 2; 2; 2; 2];
%! Mdl = fitcdiscr (x, y, 'DiscrimType', 'quadratic');
%! assert_equal (Mdl.Sigma(:,:,1), [1.6667, 1; 1, 1.6667], 1e-4);
%! assert_equal (Mdl.Coeffs(1,2).Const, 10.9525, 1e-4);
%! assert_equal (Mdl.Coeffs(1,2).Linear', [-0.8025, -0.8025], 1e-4);
%! assert_equal (Mdl.Coeffs(1,2).Quadratic, ...
%!               [-0.2587, 0.1912; 0.1912, -0.2587], 1e-4);

%!test
%! x = [1, 2; 2, 1; 3, 4; 4, 3; 5, 7; 7, 5; 8, 9; 9, 8];
%! y = [1; 1; 1; 1; 2; 2; 2; 2];
%! Mdl = fitcdiscr (x, y, 'DiscrimType', 'diagQuadratic');
%! assert_equal (Mdl.Sigma(:,:,1), [1.6667, 1.6667], 1e-4);
%! assert_equal (Mdl.Coeffs(1,2).Const, 14.8310, 1e-4);
%! assert_equal (Mdl.Coeffs(1,2).Quadratic, [-0.1286, -0.1286], 1e-4);

## The quadratic family reports a Quadratic field and the linear family does
## not, which is how a user tells the two apart from the structure alone.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! assert_equal (fieldnames (Mdl.Coeffs)', ...
%!               {'DiscrimType', 'Const', 'Linear', 'Class1', 'Class2'});

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'quadratic');
%! assert_equal (fieldnames (Mdl.Coeffs)', {'DiscrimType', 'Const', ...
%!               'Linear', 'Quadratic', 'Class1', 'Class2'});

## DiscrimType is assignable within its family and re-derives the covariance.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! Mdl.DiscrimType = 'diagLinear';
%! assert_equal (size (Mdl.Sigma), [1, 4]);
%! assert_equal (Mdl.Gamma, 1);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'diagLinear');
%! Mdl.DiscrimType = 'linear';
%! assert_equal (size (Mdl.Sigma), [4, 4]);
%! assert_equal (Mdl.Gamma, 0);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'quadratic');
%! Mdl.DiscrimType = 'diagQuadratic';
%! assert_equal (size (Mdl.Sigma), [1, 4, 3]);

## Assigning the type rebuilds the coefficients rather than leaving the ones
## the previous type produced.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! before = Mdl.Coeffs(1,2).Const;
%! Mdl.DiscrimType = 'diagLinear';
%! assert (abs (Mdl.Coeffs(1,2).Const - before) > 1);

## Gamma and the diagonal type are one state, in both directions.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! Mdl.Gamma = 1;
%! assert_equal (Mdl.DiscrimType, 'diagLinear');
%! assert_equal (size (Mdl.Sigma), [1, 4]);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'quadratic');
%! Mdl.Gamma = 1;
%! assert_equal (Mdl.DiscrimType, 'diagQuadratic');

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'Gamma', 1);
%! assert_equal (Mdl.DiscrimType, 'diagLinear');

## An intermediate Gamma shrinks the covariance towards its diagonal, which
## the resubstitution error follows.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! Mdl.Gamma = 0.25;
%! assert_equal (loss (Mdl, meas, species), 0.0267, 1e-4);

## A collinear fit is raised to MinGamma rather than failing, and cannot be
## brought back below it.
%!test
%! load fisheriris
%! Mdl = fitcdiscr ([meas(:,1:3), meas(:,3)], species);
%! assert (Mdl.MinGamma > 0);
%! assert_equal (Mdl.Gamma, Mdl.MinGamma);
%! assert_equal (Mdl.LogDetSigma, -41.3112, 1e-4);

## The pseudo types answer a singular covariance where the plain one is
## regularized and the diagonal one drops the correlations.
%!test
%! load fisheriris
%! Mdl = fitcdiscr ([meas(:,1:3), meas(:,3)], species, ...
%!                  'DiscrimType', 'pseudoLinear');
%! assert_equal (Mdl.LogDetSigma, -7.3470, 1e-4);
%! assert_equal (Mdl.Gamma, 0);

%!test
%! load fisheriris
%! Mdl = fitcdiscr ([meas(:,1:3), meas(:,3)], species, ...
%!                  'DiscrimType', 'pseudoQuadratic');
%! assert_equal (Mdl.LogDetSigma', [-11.2116, -7.2229, -6.4619], 1e-4);

## A predictor with no variance is dropped rather than inverted, so its
## coefficient is zero.
%!test
%! load fisheriris
%! Mdl = fitcdiscr ([meas(:,1:3), ones(150,1)], species, ...
%!                  'DiscrimType', 'pseudoLinear');
%! assert_equal (Mdl.Coeffs(1,2).Linear(4), 0);
%! assert_equal (Mdl.LogDetSigma, -6.3538, 1e-4);

## The compact model answers identically for every type.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'quadratic');
%! CMdl = compact (Mdl);
%! assert_equal (predict (CMdl, meas), predict (Mdl, meas));
%! assert_equal (CMdl.Sigma, Mdl.Sigma);

## A saved quadratic model comes back able to change its type, which it can
## only do if the covariance the fit estimated was saved with it.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'quadratic');
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (predict (M2, meas), predict (Mdl, meas));
%! M2.DiscrimType = 'diagQuadratic';
%! assert_equal (size (M2.Sigma), [1, 4, 3]);

## Five observations per class over two predictors, enough for a quadratic
## discriminant to have a covariance it can invert.
%!function x = dfix ()
%!  x = [1, 2; 2, 1; 3, 4; 4, 3; 2, 3; 5, 7; 7, 5; 8, 9; 9, 8; 7, 8];
%!endfunction

## Delta eliminates predictors, and the threshold is on the standardized
## coefficient: the raw one would depend on the units each predictor is
## measured in.  Values measured on MATLAB R2024a.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! assert_equal (Mdl.DeltaPredictor, ...
%!               [3.2508, 4.1236, 7.2926, 4.2506], 1e-4);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'diagLinear');
%! assert_equal (Mdl.DeltaPredictor, ...
%!               [1.6266, 1.0912, 5.3354, 4.6584], 1e-4);

## The quadratic family has no linear coefficients to eliminate.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'quadratic');
%! assert_equal (Mdl.DeltaPredictor, [0, 0, 0, 0]);

## A coefficient drops out one class at a time, and the boundary follows.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! Mdl.Delta = 2.25;
%! assert_equal (Mdl.Coeffs(1,2).Linear', ...
%!               [6.3148, 12.1393, -16.9464, -20.7701], 1e-4);
%! assert_equal (Mdl.Coeffs(1,2).Const, -14.3792, 1e-4);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! Mdl.Delta = 5;
%! assert_equal (Mdl.Coeffs(1,2).Linear', [0, 0, -16.9464, 0], 1e-4);
%! assert_equal (loss (Mdl, meas, species), 0.08, 1e-10);

## Past the largest DeltaPredictor every predictor is gone and the model
## answers from the priors alone.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! Mdl.Delta = 8;
%! assert_equal (Mdl.Coeffs(1,2).Linear', [0, 0, 0, 0]);
%! assert_equal (loss (Mdl, meas, species), 2/3, 1e-10);

## Delta changes what predict answers, which is the whole point of it.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! before = loss (Mdl, meas, species);
%! Mdl.Delta = 5;
%! assert_equal (before, 0.02, 1e-10);
%! assert (loss (Mdl, meas, species) > before);

## DeltaPredictor describes the fit, not the threshold, so assigning Delta
## leaves it where it was.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! before = Mdl.DeltaPredictor;
%! Mdl.Delta = 5;
%! assert_equal (Mdl.DeltaPredictor, before);

%!error<ClassificationDiscriminant: 'DiscrimType' must be one of the following: linear, quadratic, diagLinear, diagQuadratic, pseudoLinear, or pseudoQuadratic.> ...
%! fitcdiscr (ones (10, 2), [1;1;1;1;1;2;2;2;2;2], 'DiscrimType', 'bogus')
%!error<ClassificationDiscriminant: 'DiscrimType' can only be set to one of: linear, diagLinear, or pseudoLinear.> ...
%! Mdl = fitcdiscr (dfix (), [1;1;1;1;1;2;2;2;2;2]); ...
%! Mdl.DiscrimType = 'quadratic';
%!error<ClassificationDiscriminant: 'DiscrimType' can only be set to one of: quadratic, diagQuadratic, or pseudoQuadratic.> ...
%! Mdl = fitcdiscr (dfix (), [1;1;1;1;1;2;2;2;2;2], ...
%!                  'DiscrimType', 'quadratic'); ...
%! Mdl.DiscrimType = 'linear';
%!error<ClassificationDiscriminant: 'DiscrimType' can only be set to one of: linear, diagLinear, or pseudoLinear.> ...
%! Mdl = fitcdiscr (dfix (), [1;1;1;1;1;2;2;2;2;2]); Mdl.DiscrimType = 'bogus';
%!error<ClassificationDiscriminant: cannot set 'Gamma' to any value but 0 or 1 for a quadratic discriminant.> ...
%! Mdl = fitcdiscr (dfix (), [1;1;1;1;1;2;2;2;2;2], ...
%!                  'DiscrimType', 'quadratic'); ...
%! Mdl.Gamma = 0.5;
%!error<ClassificationDiscriminant: cannot set 'Gamma' to any value but 0 or 1 for a quadratic discriminant.> ...
%! fitcdiscr (dfix (), [1;1;1;1;1;2;2;2;2;2], ...
%!            'DiscrimType', 'quadratic', 'Gamma', 0.5)
%!error<ClassificationDiscriminant: 'Gamma' must be a scalar between 0 and 1.> ...
%! Mdl = fitcdiscr (dfix (), [1;1;1;1;1;2;2;2;2;2]); Mdl.Gamma = 1.5;
%!error<ClassificationDiscriminant: cannot eliminate linear predictors in a quadratic discriminant.> ...
%! Mdl = fitcdiscr (dfix (), [1;1;1;1;1;2;2;2;2;2], ...
%!                  'DiscrimType', 'quadratic'); ...
%! Mdl.Delta = 0.5;
%!error<ClassificationDiscriminant: cannot set linear coefficients to zero for a quadratic discriminant.> ...
%! fitcdiscr (dfix (), [1;1;1;1;1;2;2;2;2;2], ...
%!            'DiscrimType', 'quadratic', 'Delta', 0.5)
%!error<ClassificationDiscriminant: 'Delta' must be a nonnegative scalar.> ...
%! fitcdiscr (dfix (), [1;1;1;1;1;2;2;2;2;2], 'Delta', -1)
%!error<ClassificationDiscriminant: 'Delta' must be a nonnegative scalar.> ...
%! Mdl = fitcdiscr (dfix (), [1;1;1;1;1;2;2;2;2;2]); Mdl.Delta = [1, 2];
%!error<ClassificationDiscriminant: cannot use 'quadratic' type because one or more classes have singular covariance matrices.> ...
%! load fisheriris; ...
%! fitcdiscr ([meas(:,1:3), meas(:,3)], species, 'DiscrimType', 'quadratic')
%!error<ClassificationDiscriminant: predictor 'x4' has zero variance for class 'setosa'.  Either exclude this predictor or set 'DiscrimType' to 'pseudoQuadratic' or 'diagQuadratic'.> ...
%! load fisheriris; ...
%! fitcdiscr ([meas(:,1:3), ones(150,1)], species, 'DiscrimType', 'quadratic')

## MinGamma describes the data rather than the type, so a pseudo fit reports
## the same value as a plain one on the same predictors while regularizing by
## nothing at all.
%!test
%! load fisheriris
%! Mdl = fitcdiscr ([meas(:,1:3), meas(:,3)], species, ...
%!                  'DiscrimType', 'pseudoLinear');
%! assert (Mdl.MinGamma > 0);
%! assert_equal (Mdl.Gamma, 0);

%!test
%! load fisheriris
%! Mdl = fitcdiscr ([meas(:,1:3), meas(:,3)], species, ...
%!                  'DiscrimType', 'diagLinear');
%! assert (Mdl.MinGamma > 0);
%! assert_equal (Mdl.Gamma, 1);

## edge, resubPredict, resubMargin, resubEdge and resubLoss.  Every value
## below was measured on MATLAB R2024a; the discriminant's scores agree with
## it exactly, so these are the oracle's own numbers and not ours.

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! assert_equal (edge (Mdl, meas, species), 0.9454289377, 1e-9);

## The edge is the mean of the margins, by definition.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! m = margin (Mdl, meas, species);
%! assert_equal (edge (Mdl, meas, species), mean (m), 1e-12);

## Weights are normalized within each class to that class's prior, which is
## not the same as dividing by their total: that would give 0.9269539697.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! assert_equal (edge (Mdl, meas, species, 'Weights', (1:150)'), ...
%!               0.9438468986, 1e-9);

## A weight constant within a class therefore changes nothing.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! w = [ones(50,1); 7 * ones(50,1); 0.5 * ones(50,1)];
%! assert_equal (edge (Mdl, meas, species, 'Weights', w), ...
%!               edge (Mdl, meas, species), 1e-12);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! assert_equal (resubEdge (Mdl), 0.9454289377, 1e-9);
%! assert_equal (resubLoss (Mdl), 0.02, 1e-12);
%! assert_equal (sum (resubMargin (Mdl)), 141.8143406564, 1e-8);

## resubPredict is predict on the data the model kept.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! [label, score, cost] = resubPredict (Mdl);
%! [l2, s2, c2] = predict (Mdl, meas);
%! assert_equal (label, l2);
%! assert_equal (score, s2);
%! assert_equal (cost, c2);
%! assert_equal (numel (label), 150);

## The compact model answers the same edge as the full one.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! assert_equal (edge (compact (Mdl), meas, species), 0.9454289377, 1e-9);

%!error<ClassificationDiscriminant.edge: too few input arguments.> ...
%! load fisheriris; edge (fitcdiscr (meas, species), meas)
%!error<ClassificationDiscriminant.edge: Name-Value arguments must be in pairs.> ...
%! load fisheriris; ...
%! edge (fitcdiscr (meas, species), meas, species, 'Weights')
%!error<ClassificationDiscriminant.edge: size of 'Weights' must equal the number of rows in X.> ...
%! load fisheriris; ...
%! edge (fitcdiscr (meas, species), meas, species, 'Weights', ones (3, 1))
%!error<ClassificationDiscriminant.edge: invalid parameter name in optional paired arguments.> ...
%! load fisheriris; ...
%! edge (fitcdiscr (meas, species), meas, species, 'Nope', 1)

## BinEdges is an empty cell, which is what MATLAB reports for this
## learner as well: it fits the predictors as they are.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! assert_equal (class (Mdl.BinEdges), 'cell');
%! assert_equal (Mdl.BinEdges, {});

## CategoricalPredictors and ExpandedPredictorNames, shapes measured on
## R2024a: an empty double and one name per predictor.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! assert_equal (Mdl.CategoricalPredictors, []);
%! assert_equal (size (Mdl.CategoricalPredictors), [0, 0]);
%! assert_equal (Mdl.ExpandedPredictorNames, Mdl.PredictorNames);
%! assert_equal (size (Mdl.ExpandedPredictorNames), [1, 4]);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (M2.CategoricalPredictors, Mdl.CategoricalPredictors);
%! assert_equal (M2.ExpandedPredictorNames, Mdl.ExpandedPredictorNames);

## BetweenSigma, measured on MATLAB R2024a.  It weights the class means by
## class size, not by Prior, and the denominator is the unbiased one for a
## weighted covariance.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! assert_equal (Mdl.BetweenSigma, ...
%!               [0.632121333333333, -0.199526666666667, ...
%!                1.652483999999999, 0.712793333333333;
%!                -0.199526666666667, 0.113449333333333, ...
%!                -0.572395999999999, -0.229326666666666;
%!                1.652483999999999, -0.572395999999999, ...
%!                4.371027999999996, 1.867739999999998;
%!                0.712793333333333, -0.229326666666666, ...
%!                1.867739999999998, 0.804133333333333], 1e-12);

## Unequal class sizes: 50, 30 and 10 of the three species.
%!test
%! load fisheriris
%! k = [1:50, 51:80, 101:110];
%! Mdl = fitcdiscr (meas(k,:), species(k));
%! assert_equal (Mdl.BetweenSigma, ...
%!               [0.651346086956520, -0.299426956521739, ...
%!                1.775435652173911, 0.711567826086955;
%!                -0.299426956521739, 0.160084347826087, ...
%!                -0.811819130434784, -0.318816086956522;
%!                1.775435652173911, -0.811819130434784, ...
%!                4.840319130434781, 1.941198695652173;
%!                0.711567826086955, -0.318816086956522, ...
%!                1.941198695652173, 0.780424347826086], 1e-12);

## A row missing a predictor still counts towards the others, as it does for
## Mu; only a missing response drops one.
%!test
%! load fisheriris
%! X = meas; X(3,2) = NaN; X(77,4) = NaN;
%! Mdl = fitcdiscr (X, species);
%! assert_equal (Mdl.BetweenSigma(1,2), -0.201474748299320, 1e-12);
%! assert_equal (Mdl.BetweenSigma(4,4), 0.803942801055115, 1e-12);

## It describes the classes, not the fit, so the quadratic family reports it
## and assigning a prior does not move it.
%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'quadratic');
%! assert_equal (Mdl.BetweenSigma(1,1), 0.632121333333333, 1e-12);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! before = Mdl.BetweenSigma;
%! Mdl.Prior = [0.6, 0.2, 0.2];
%! assert_equal (Mdl.BetweenSigma, before);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species);
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (M2.BetweenSigma, Mdl.BetweenSigma);
