## Copyright (C) 2024-2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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

classdef CompactClassificationDiscriminant
  ## -*- texinfo -*-
  ## @deftp {statistics} CompactClassificationDiscriminant
  ##
  ## Compact discriminant analysis classification
  ##
  ## The @code{CompactClassificationDiscriminant} class implements a compact
  ## version of a linear discriminant analysis classifier object, which can
  ## predict responses for new data using the @code{predict} method but does not
  ## store the training data.
  ##
  ## A @code{CompactClassificationDiscriminant} object is a compact version of a
  ## discriminant analysis model, @code{ClassificationDiscriminant}.  It does
  ## not include the training data resulting in a smaller classifier size, which
  ## can be used for making predictions from new data, but not for tasks such as
  ## cross validation.  It can only be created from a
  ## @code{ClassificationDiscriminant} model by using the @code{compact} object
  ## method.
  ##
  ## Create a @code{CompactClassificationDiscriminant} object by using the
  ## @code{compact} method of a @code{ClassificationDiscriminant} object.
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
  ## @seealso{fitcdiscr, ClassificationDiscriminant}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## A positive integer value specifying the number of predictors in the
    ## training dataset used for training the CompactClassificationDiscriminant
    ## model.  This property is read-only.
    ##
    ## @end deftp
    NumPredictors   = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property} PredictorNames
    ##
    ## Names of predictor variables
    ##
    ## A cell array of character vectors specifying the names of the predictor
    ## variables.  The names are in the order in which they appear in the
    ## training dataset.  This property is read-only.
    ##
    ## @end deftp
    PredictorNames  = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property} ResponseName
    ##
    ## Response variable name
    ##
    ## A character vector specifying the name of the response variable @var{Y}.
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName    = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property} ClassNames
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
    ## @deftp {CompactClassificationDiscriminant} {property} Sigma
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
    ## @deftp {CompactClassificationDiscriminant} {property} Mu
    ##
    ## Class means
    ##
    ## A @math{K*P} numeric matrix specifying the mean of the multivariate
    ## normal distribution of each corresponding class, where @math{K} is the
    ## number of classes and @math{P} is the number of predictors. This property
    ## is read-only.
    ##
    ## @end deftp
    Mu              = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property} Coeffs
    ##
    ## Coefficient matrices
    ##
    ## A @math{K*K} structure containing the coefficient matrices, where
    ## @math{K} is the number of classes.  If the @qcode{'FillCoeffs'} parameter
    ## was set to @qcode{'off'} in the original
    ## @code{ClassificationDiscriminant} model, then @qcode{Coeffs} is empty
    ## @qcode{([])}. This property is read-only.
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
    ## @deftp {CompactClassificationDiscriminant} {property} MinGamma
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
    ## @deftp {CompactClassificationDiscriminant} {property} LogDetSigma
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
  endproperties

  ## Properties a user may set after the model is built.  Each one is
  ## validated by its set method below.
  properties (GetAccess = public, SetAccess = public)

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property} DiscrimType
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
    ## @deftp {CompactClassificationDiscriminant} {property} Gamma
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
    ## @deftp {CompactClassificationDiscriminant} {property} Delta
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
    ## @deftp {CompactClassificationDiscriminant} {property} Cost
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
    ## This property is read-only.
    ##
    ## @end deftp
    Cost            = [];

    ## -*- texinfo -*-
    ## @deftp {CompactClassificationDiscriminant} {property} Prior
    ##
    ## Prior probability for each class
    ##
    ## A numeric vector specifying the prior probabilities for each class.  The
    ## order of the elements in @qcode{Prior} corresponds to the order of the
    ## classes in @qcode{ClassNames}.
    ##
    ## This property is read-only.
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
    ## @deftp {CompactClassificationDiscriminant} {property} ScoreTransform
    ##
    ## Transformation function for classification scores
    ##
    ## Specified as a function handle for transforming the classification
    ## scores.  This property is read-only.
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

  ## Readable by the counterpart class, which copies it, and kept out of
  ## the documented surface.
  properties (GetAccess = public, SetAccess = protected, Hidden)

    ## The unregularized within-class covariance the fit estimated, carried
    ## over from the full model so that assigning DiscrimType re-derives here
    ## exactly as it does there.
    BaseSigma = [];
    STfun = @(x) x;
  endproperties

  ## Set methods for the properties a user may assign.
  methods (Hidden)

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
        error (strcat ("CompactClassificationDiscriminant:", ...
                       " 'DiscrimType' can only be set to one of: %s."), ok);
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
        error (strcat ("CompactClassificationDiscriminant: 'Gamma' must", ...
                       " be a scalar between 0 and 1."));
      endif
      if (isempty (this.BaseSigma))
        this.Gamma = val;
        return;
      endif
      [~, fam] = discrimcanon (this.DiscrimType);
      if (strcmp (fam, 'quadratic') && val > 0 && val < 1)
        error (strcat ("CompactClassificationDiscriminant: cannot set", ...
                       " 'Gamma' to any value but 0 or 1 for a quadratic", ...
                       " discriminant."));
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
        error (strcat ("CompactClassificationDiscriminant: 'Gamma' must be", ...
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
        error (strcat ("CompactClassificationDiscriminant: 'Delta' must", ...
                       " be a nonnegative scalar."));
      endif
      [~, fam] = discrimcanon (this.DiscrimType);
      if (val > 0 && strcmp (fam, 'quadratic'))
        error (strcat ("CompactClassificationDiscriminant: cannot", ...
                       " eliminate linear predictors in a quadratic", ...
                       " discriminant."));
      endif
      this.Delta = val;
      if (! isempty (this.BaseSigma) && ! isempty (this.Coeffs))
        this.Coeffs = discrimcoeffs (this.Mu, this.Sigma, this.LogDetSigma, ...
                                     this.Prior, this.DiscrimType, val, ...
                                     this.ClassNames);
      endif
    endfunction

    function this = set.ScoreTransform (this, val)
      name = 'CompactClassificationDiscriminant';
      [this.STfun, this.ScoreTransform] = parseScoreTransform (val, name);
    endfunction

    function this = set.Cost (this, val)
      gnY = this.ClassNames;
      if (isempty (val))
        this.Cost = cast (! eye (numel (gnY)), 'double');
      else
        if (numel (gnY) != sqrt (numel (val)))
          error (strcat ("CompactClassificationDiscriminant: the number", ...
                         " of rows and columns in 'Cost' must correspond", ...
                         " to the selected classes in Y."));
        endif
        this.Cost = val;
      endif
    endfunction

    function this = set.Prior (this, val)
      K = numel (this.ClassNames);
      if (isstruct (val))
        val = priorFromStruct (val, this.ClassNames, ...
                               'CompactClassificationDiscriminant');
      endif
      if (ischar (val) && strcmpi (val, 'uniform'))
        this.Prior = ones (1, K) / K;
      elseif (isnumeric (val) && isreal (val) && isvector (val)
              && numel (val) == K && all (val >= 0) && sum (val) > 0)
        this.Prior = val(:)' / sum (val);
      else
        error (strcat ("CompactClassificationDiscriminant: 'Prior' must be", ...
                       " 'uniform' or a non-negative numeric vector with", ...
                       " one entry per class."));
      endif
    endfunction

    ## Constructor
    function this = CompactClassificationDiscriminant (Mdl = [])

      ## Check for appropriate class
      if (isempty (Mdl))
        return;
      elseif (! strcmpi (class (Mdl), 'ClassificationDiscriminant'))
        error (strcat ("CompactClassificationDiscriminant: invalid", ...
                       " classification object."));
      endif

      ## Save properties to compact model
      this.NumPredictors   = Mdl.NumPredictors;
      this.PredictorNames  = Mdl.PredictorNames;
      this.ResponseName    = Mdl.ResponseName;
      this.ClassNames      = Mdl.ClassNames;

      this.Cost            = Mdl.Cost;
      this.Prior           = Mdl.Prior;
      this.ScoreTransform  = Mdl.ScoreTransform;
      this.STfun          = Mdl.STfun;

      this.Sigma           = Mdl.Sigma;
      this.Mu              = Mdl.Mu;
      this.Coeffs          = Mdl.Coeffs;
      this.Delta           = Mdl.Delta;
      this.DiscrimType     = Mdl.DiscrimType;
      this.Gamma           = Mdl.Gamma;
      this.MinGamma        = Mdl.MinGamma;
      this.LogDetSigma     = Mdl.LogDetSigma;
      this.DeltaPredictor  = Mdl.DeltaPredictor;
      this.BaseSigma       = Mdl.BaseSigma;

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
      fprintf ("\n  CompactClassificationDiscriminant\n\n");
      ## Print selected properties
      fprintf ("%+25s: '%s'\n", 'ResponseName', this.ResponseName);
      if (iscellstr (this.ClassNames))
        str = repmat ({'''%s'''}, 1, numel (this.ClassNames));
        str = strcat ('{', strjoin (str, ' '), '}');
        str = sprintf (str, this.ClassNames{:});
      else # numeric
        str = repmat ({'%d'}, 1, numel (this.ClassNames));
        str = strcat ('[', strjoin (str, ' '), ']');
        str = sprintf (str, this.ClassNames);
      endif
      fprintf ("%+25s: %s\n", 'ClassNames', str);
      fprintf ("%+25s: '%s'\n", 'ScoreTransform', this.ScoreTransform);
      fprintf ("%+25s: '%d'\n", 'NumPredictors', this.NumPredictors);
      fprintf ("%+25s: '%s'\n", 'DiscrimType', this.DiscrimType);
      fprintf ("%+25s: [%dx%d double]\n", 'Mu', size (this.Mu));
      ## Coeffs is KxK, which is Sigma's shape only for a linear discriminant.
      fprintf ("%+25s: [%dx%d struct]\n\n", 'Coeffs', ...
               rows (this.ClassNames), rows (this.ClassNames));
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationDiscriminant} {@var{label} =} predict (@var{obj}, @var{XC})
    ## @deftypefnx {CompactClassificationDiscriminant} {[@var{label}, @var{score}, @var{cost}] =} predict (@var{obj}, @var{XC})
    ##
    ## Classify new data points into categories using the discriminant
    ## analysis model from a CompactClassificationDiscriminant object.
    ##
    ## @code{@var{label} = predict (@var{obj}, @var{XC})} returns the vector of
    ## labels predicted for the corresponding instances in @var{XC}, using the
    ## corresponding labels from the trained @qcode{ClassificationDiscriminant},
    ## model, @var{obj}.
    ##
    ## @itemize
    ## @item
    ## @var{obj} must be a @qcode{CompactClassificationDiscriminant} class
    ## object.
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
    ## @seealso{CompactClassificationDiscriminant, fitcdiscr}
    ## @end deftypefn
    function [label, score, cost] = predict (this, XC)

      ## Check for sufficient input arguments
      if (nargin < 2)
        error (strcat ("CompactClassificationDiscriminant.predict:", ...
                       " too few input arguments."));
      endif

      ## Check for valid XC
      if (isempty (XC))
        error ("CompactClassificationDiscriminant.predict: XC is empty.");
      elseif (this.NumPredictors != columns (XC))
        error (strcat ("CompactClassificationDiscriminant.predict: XC", ...
                       " must have the same number of features as the", ...
                       " trained model."));
      endif

      ## Initialize matrices
      numObservations = rows (XC);
      numClasses = numel (this.ClassNames);
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
      label = this.ClassNames(minIdx);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationDiscriminant} {@var{L} =} loss (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationDiscriminant} {@var{L} =} loss (@dots{}, @var{name}, @var{value})
    ##
    ## Compute loss for a trained CompactClassificationDiscriminant object.
    ##
    ## @code{@var{L} = loss (@var{obj}, @var{X}, @var{Y})} computes the loss,
    ## @var{L}, using the default loss function @qcode{'mincost'}.
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{CompactClassificationDiscriminant} object.
    ## @item
    ## @code{X} must be a @math{N*P} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of rows as @var{X}.
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
    ## @seealso{CompactClassificationDiscriminant}
    ## @end deftypefn
    function L = loss (this, X, Y, varargin)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error (strcat ("CompactClassificationDiscriminant.loss:", ...
                       " too few input arguments."));
      elseif (mod (nargin - 3, 2) != 0)
        error (strcat ("CompactClassificationDiscriminant.loss:", ...
                       " name-value arguments must be in pairs."));
      elseif (nargin > 7)
        error (strcat ("CompactClassificationDiscriminant.loss:", ...
                       " too many input arguments."));
      endif

      ## Default values
      LossFun = 'mincost';
      Weights = [];

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error (strcat ("CompactClassificationDiscriminant.loss:", ...
                       " Y must be of a valid type."));
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat ("CompactClassificationDiscriminant.loss: Y must", ...
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
                error (strcat ("CompactClassificationDiscriminant.loss:", ...
                               " custom loss function must accept", ...
                               " exactly four input arguments."));
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
                  error (strcat ("CompactClassificationDiscriminant.loss:", ...
                                 " custom loss function must return", ...
                                 " a scalar value."));
                endif
              catch
                error (strcat ("CompactClassificationDiscriminant.loss:", ...
                               " custom loss function is not valid or", ...
                               " does not produce correct output."));
              end_try_catch
              LossFun = Value;
            elseif (ischar (Value) && any (strcmpi (Value, lf_opt)))
              LossFun = Value;
            else
              error (strcat ("CompactClassificationDiscriminant.loss:", ...
                             " invalid loss function."));
            endif

          case 'weights'
            if (isnumeric (Value) && isvector (Value))
              if (numel (Value) != size (X ,1))
                error (strcat ("CompactClassificationDiscriminant.loss:", ...
                               " number of 'Weights' must be equal to", ...
                               " the number of rows in X."));
              elseif (numel (Value) == size (X, 1))
                Weights = Value;
              endif
            else
              error (strcat ("CompactClassificationDiscriminant.loss:", ...
                             " invalid 'Weights'."));
            endif

          otherwise
            error (strcat ("CompactClassificationDiscriminant.loss:", ...
                           " invalid parameter name in optional pair", ...
                           " arguments."));
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

      ## Convert Y to a cell array of strings
      if (ischar (Y))
        Y = cellstr (Y);
      elseif (isnumeric (Y))
        Y = cellstr (num2str (Y));
      elseif (islogical (Y))
        Y = cellstr (num2str (double (Y)));
      elseif (iscell (Y))
        Y = cellfun (@num2str, Y, 'UniformOutput', false);
      else
        error (strcat ("CompactClassificationDiscriminant.loss: Y must be", ...
                       " a numeric, logical, char, string, or cell array."));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), this.ClassNames)))
        error (strcat ("CompactClassificationDiscriminant.loss: Y must", ...
                       " contain only the classes in ClassNames."));
      endif

      ## Set default weights if not specified
      if (isempty (Weights))
        Weights = ones (size (X, 1), 1);
      endif

      ## Normalize Weights
      unique_classes = this.ClassNames;
      class_prior_probs = this.Prior;
      norm_weights = zeros (size (Weights));
      for i = 1:numel (unique_classes)
        class_idx = ismember (Y, unique_classes{i});
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
      K = numel (this.ClassNames);
      C = false (n, K);
      for i = 1:n
        class_idx = find (ismember (this.ClassNames, Y{i}));
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
            cj = Cost(find (ismember (this.ClassNames, Y(i))), min_cost_class);
            L = L + Weights(i) * cj;
          endfor
        case 'classifcost'
          Cost = this.Cost;
          L = 0;
          for i = 1:n
            y_idx = find (ismember (this.ClassNames, Y(i)));
            y_hat_idx = find (ismember (this.ClassNames, label(i)));
            L = L + Weights(i) * Cost(y_idx, y_hat_idx);
          endfor
        otherwise
          error (strcat ("CompactClassificationDiscriminant.loss:", ...
                         " invalid loss function."));
      endswitch

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn {CompactClassificationDiscriminant} {@var{m} =} margin (@var{obj}, @var{X}, @var{Y})
    ##
    ## Classification margins for discriminant analysis classifier.
    ##
    ## @code{@var{m} = margin (@var{obj}, @var{X}, @var{Y})} returns
    ## the classification margins for @var{obj} with data @var{X} and
    ## classification @var{Y}. @var{m} is a numeric vector of length size (X,1).
    ##
    ## @itemize
    ## @item
    ## @code{obj} is a @var{CompactClassificationDiscriminant} object.
    ## @item
    ## @code{X} must be a @math{N*P} numeric matrix of input data where rows
    ## correspond to observations and columns correspond to features or
    ## variables.
    ## @item
    ## @code{Y} is @math{N*1} matrix or cell matrix containing the class labels
    ## of corresponding predictor data in @var{X}. @var{Y} must have same
    ## numbers of rows as @var{X}.
    ## @end itemize
    ##
    ## The classification margin for each observation is the difference between
    ## the classification score for the true class and the maximal
    ## classification score for the false classes.
    ##
    ## @seealso{fitcdiscr, CompactClassificationDiscriminant}
    ## @end deftypefn
    function m = margin (this, X, Y)

      ## Check for sufficient input arguments
      if (nargin < 3)
        error (strcat ("CompactClassificationDiscriminant.margin:", ...
                       " too few input arguments."));
      endif

      ## Validate Y
      valid_types = {'char', 'string', 'logical', 'single', 'double', 'cell'};
      if (! (any (strcmp (class (Y), valid_types))))
        error (strcat ("CompactClassificationDiscriminant.margin:", ...
                       " Y must be of a valid type."));
      endif

      ## Validate X
      valid_types = {'single', 'double'};
      if (! (any (strcmp (class (X), valid_types))))
        error (strcat ("CompactClassificationDiscriminant.margin:", ...
                       " X must be of a valid type."));
      endif

      ## Validate size of Y
      if (size (Y, 1) != size (X, 1))
        error (strcat ("CompactClassificationDiscriminant.margin: Y must", ...
                       " have the same number of rows as X."));
      endif

      ## Convert Y to a cell array of strings
      if (ischar (Y))
        Y = cellstr (Y);
      elseif (isnumeric (Y))
        Y = cellstr (num2str (Y));
      elseif (islogical (Y))
        Y = cellstr (num2str (double (Y)));
      elseif (iscell (Y))
        Y = cellfun (@num2str, Y, 'UniformOutput', false);
      else
        error (strcat ("CompactClassificationDiscriminant.margin: Y must", ...
                       " be a numeric, logical, char, string, or cell array."));
      endif

      ## Check if Y contains correct classes
      if (! all (ismember (unique (Y), this.ClassNames)))
        error (strcat ("CompactClassificationDiscriminant.margin: Y must", ...
                       " contain only the classes in ClassNames."));
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
        true_class_idx = find (ismember (this.ClassNames, Y{i}));

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
    ## @deftypefn  {CompactClassificationDiscriminant} {@var{e} =} edge (@var{obj}, @var{X}, @var{Y})
    ## @deftypefnx {CompactClassificationDiscriminant} {@var{e} =} edge (@dots{}, @qcode{"Weights"}, @var{w})
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
        error (strcat ("CompactClassificationDiscriminant.edge: too few", ...
                       " input arguments."));
      endif
      if (mod (numel (varargin), 2) != 0)
        error (strcat ("CompactClassificationDiscriminant.edge: Name-Value", ...
                       " arguments must be in pairs."));
      endif

      ## The weights are parsed before anything is computed, so a bad
      ## Name-Value pair is reported as such rather than after a margin.
      W = edgeWeights (varargin, Y, this.ClassNames, this.Prior, ...
                       "CompactClassificationDiscriminant", "edge");
      m = margin (this, X, Y);
      e = sum (W .* m(:)) / sum (W);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {CompactClassificationDiscriminant} {} savemodel (@var{obj}, @var{filename})
    ##
    ## Save a CompactClassificationDiscriminant object.
    ##
    ## @code{savemodel (@var{obj}, @var{filename})} saves each property of a
    ## CompactClassificationDiscriminant object into an Octave binary file, the
    ## name of which is specified in @var{filename}, along with an extra
    ## variable, which defines the type classification object these variables
    ## constitute. Use @code{loadmodel} in order to load a classification object
    ## into Octave's workspace.
    ##
    ## @seealso{loadmodel, fitcdiscr, ClassificationDiscriminant}
    ## @end deftypefn
    function savemodel (this, fname)
      if (nargin < 2)
        error ("CompactClassificationDiscriminant.savemodel: too few input arguments.");
      endif
      if (! (ischar (fname) && isrow (fname) && ! isempty (fname)))
        error ("CompactClassificationDiscriminant.savemodel: FNAME must be a character vector.");
      endif
      ## Generate variable for class name
      classdef_name = 'CompactClassificationDiscriminant';

      ## Create variables from model properties
      NumPredictors   = this.NumPredictors;
      PredictorNames  = this.PredictorNames;
      ResponseName    = this.ResponseName;
      ClassNames      = this.ClassNames;
      Prior           = this.Prior;
      Cost            = this.Cost;
      ScoreTransform  = this.ScoreTransform;
      STfun          = this.STfun;
      Sigma           = this.Sigma;
      BaseSigma       = this.BaseSigma;
      Mu              = this.Mu;
      Coeffs          = this.Coeffs;
      Delta           = this.Delta;
      DiscrimType     = this.DiscrimType;
      Gamma           = this.Gamma;
      MinGamma        = this.MinGamma;
      LogDetSigma     = this.LogDetSigma;

      ## Save classdef name and all model properties as individual variables
      save ('-binary', fname, 'classdef_name', 'NumPredictors', ...
            'PredictorNames', 'ResponseName', 'ClassNames', 'Prior', ...
            'Cost', 'ScoreTransform', 'STfun', 'Sigma', 'BaseSigma', ...
            'Mu', 'Coeffs', ...
            'Delta', 'DiscrimType', 'Gamma', 'MinGamma', 'LogDetSigma');
    endfunction

  endmethods

  methods(Static, Hidden)

    function mdl = load_model (filename, data)
      ## Create a CompactClassificationDiscriminant object
      mdl = CompactClassificationDiscriminant ();

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
          msg = strcat ("CompactClassificationDiscriminant.load_model:", ...
                        " invalid model in '%s'.");
          error (msg, filename);
        end_try_catch
      endfor
    endfunction

  endmethods

endclassdef

%!demo
%! ## Create a discriminant analysis classifier and its compact version
%! # and compare their size
%!
%! load fisheriris
%! X = meas;
%! Y = species;
%!
%! Mdl = fitcdiscr (X, Y, 'ClassNames', unique (species))
%! CMdl = crossval (Mdl)

## Test constructor
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! PredictorNames = {'Sepal Length', 'Sepal Width', 'Petal Length', 'Petal Width'};
%! Mdl = fitcdiscr (x, y, 'PredictorNames', PredictorNames);
%! CMdl = compact (Mdl);
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
%! assert_equal (class (CMdl), "CompactClassificationDiscriminant");
%! assert_equal ({CMdl.DiscrimType, CMdl.ResponseName}, {'linear', 'Y'})
%! assert_equal ({CMdl.Gamma, CMdl.MinGamma}, {0, 0}, 1e-15)
%! assert_equal (CMdl.ClassNames, unique (species))
%! assert_equal (CMdl.Sigma, sigma, 1e-6)
%! assert_equal (CMdl.Mu, mu, 1e-14)
%! assert_equal (CMdl.LogDetSigma, -9.9585, 1e-4)
%! assert_equal (CMdl.PredictorNames, PredictorNames)
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! Mdl = fitcdiscr (x, y, 'Gamma', 0.5);
%! CMdl = compact (Mdl);
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
%! assert_equal (class (CMdl), "CompactClassificationDiscriminant");
%! assert_equal ({CMdl.DiscrimType, CMdl.ResponseName}, {'linear', 'Y'})
%! assert_equal ({CMdl.Gamma, CMdl.MinGamma}, {0.5, 0})
%! assert_equal (CMdl.ClassNames, unique (species))
%! assert_equal (CMdl.Sigma, sigma, 1e-6)
%! assert_equal (CMdl.Mu, mu, 1e-14)
%! assert_equal (CMdl.LogDetSigma, -8.6884, 1e-4)

## Test input validation for constructor
%!error<CompactClassificationDiscriminant: invalid classification object.> ...
%! CompactClassificationDiscriminant (1)

## Test predict method
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! Mdl = fitcdiscr (meas, species, 'Gamma', 0.5);
%! CMdl = compact (Mdl);
%! [label, score, cost] = predict (CMdl, [2, 2, 2, 2]);
%! assert_equal (label, {'versicolor'})
%! assert_equal (score, [0, 0.9999, 0.0001], 1e-4)
%! assert_equal (cost, [1, 0.0001, 0.9999], 1e-4)
%! [label, score, cost] = predict (CMdl, [2.5, 2.5, 2.5, 2.5]);
%! assert_equal (label, {'versicolor'})
%! assert_equal (score, [0, 0.6368, 0.3632], 1e-4)
%! assert_equal (cost, [1, 0.3632, 0.6368], 1e-4)
%!test
%! load fisheriris
%! x = meas;
%! y = species;
%! xc = [min(x); mean(x); max(x)];
%! Mdl = fitcdiscr (x, y);
%! CMdl = compact (Mdl);
%! [label, score, cost] = predict (CMdl, xc);
%! l = {'setosa'; 'versicolor'; 'virginica'};
%! s = [1, 0, 0; 0, 1, 0; 0, 0, 1];
%! c = [0, 1, 1; 1, 0, 1; 1, 1, 0];
%! assert_equal (label, l)
%! assert_equal (score, s, 1e-4)
%! assert_equal (cost, c, 1e-4)

%!shared MODEL
%! X = rand (10,2);
%! Y = [ones(5,1);2*ones(5,1)];
%! MODEL = compact (ClassificationDiscriminant (X, Y));

## Test input validation for predict method
%!error<CompactClassificationDiscriminant.predict: too few input arguments.> ...
%! predict (MODEL)
%!error<CompactClassificationDiscriminant.predict: XC is empty.> ...
%! predict (MODEL, [])
%!error<CompactClassificationDiscriminant.predict: XC must have the same number of features as the trained model.> ...
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
%!error<CompactClassificationDiscriminant.loss: too few input arguments.> ...
%! loss (MODEL)
%!error<CompactClassificationDiscriminant.loss: too few input arguments.> ...
%! loss (MODEL, ones (4,2))
%!error<CompactClassificationDiscriminant.loss: name-value arguments must be in pairs.> ...
%! loss (MODEL, ones (4,2), ones (4,1), 'LossFun')
%!error<CompactClassificationDiscriminant.loss: Y must have the same number of rows as X.> ...
%! loss (MODEL, ones (4,2), ones (3,1))
%!error<CompactClassificationDiscriminant.loss: invalid loss function.> ...
%! loss (MODEL, ones (4,2), ones (4,1), 'LossFun', 'a')
%!error<CompactClassificationDiscriminant.loss: invalid 'Weights'.> ...
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
%!error<CompactClassificationDiscriminant.margin: too few input arguments.> ...
%! margin (MODEL)
%!error<CompactClassificationDiscriminant.margin: too few input arguments.> ...
%! margin (MODEL, ones (4,2))
%!error<CompactClassificationDiscriminant.margin: Y must have the same number of rows as X.> ...
%! margin (MODEL, ones (4,2), ones (3,1))
%!error <CompactClassificationDiscriminant.savemodel: too few input arguments.> ...
%! savemodel (CompactClassificationDiscriminant ())
%!error <CompactClassificationDiscriminant.savemodel: FNAME must be a character vector.> ...
%! savemodel (CompactClassificationDiscriminant (), 1)
%!error <CompactClassificationDiscriminant.savemodel: FNAME must be a character vector.> ...
%! savemodel (CompactClassificationDiscriminant (), ['ab'; 'cd'])

## A ScoreTransform can be assigned, and is stored as a function handle.
%!test
%! load fisheriris
%! CMdl = compact (fitcdiscr (meas, species));
%! CMdl.ScoreTransform = 'symmetric';
%! assert_equal (class (CMdl.ScoreTransform), 'char');
%! assert_equal (CMdl.ScoreTransform, 'symmetric');
%! assert_equal (CMdl.STfun (0.25), -0.5);

## Prior and Cost are settable on the compact model, as they are on the full
## one, and an unnormalized prior is rescaled.
%!test
%! load fisheriris
%! CMdl = compact (fitcdiscr (meas, species));
%! CMdl.Prior = [2, 3, 5];
%! assert_equal (CMdl.Prior, [0.2, 0.3, 0.5], 1e-14);
%! CMdl.Cost = [0, 2, 3; 1, 0, 1; 1, 1, 0];
%! assert_equal (CMdl.Cost, [0, 2, 3; 1, 0, 1; 1, 1, 0]);

## A fitted model survives savemodel and loadmodel: the properties come
## back as they were and it predicts the same.
%!test
%! load fisheriris
%! Mdl = compact (fitcdiscr (meas, species));
%! fname = tempname ();
%! savemodel (Mdl, fname);
%! M2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (class (M2), 'CompactClassificationDiscriminant');
%! assert_equal (M2.PredictorNames, Mdl.PredictorNames);
%! assert_equal (class (M2.ScoreTransform), class (Mdl.ScoreTransform));
%! assert_equal (predict (M2, meas(1:5,:)), predict (Mdl, meas(1:5,:)));

## The compact model carries the discriminant type and answers identically to
## the full one, for every type.  Values measured on MATLAB R2024a; see
## DISCRIMINANT_LEDGER.md.

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'quadratic');
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.DiscrimType, 'quadratic');
%! assert_equal (size (CMdl.Sigma), [4, 4, 3]);
%! assert_equal (predict (CMdl, meas), predict (Mdl, meas));

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'diagQuadratic');
%! CMdl = compact (Mdl);
%! assert_equal (size (CMdl.Sigma), [1, 4, 3]);
%! assert_equal (predict (CMdl, meas), predict (Mdl, meas));

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'pseudoLinear');
%! CMdl = compact (Mdl);
%! assert_equal (CMdl.LogDetSigma, -9.9585, 1e-4);
%! assert_equal (predict (CMdl, meas), predict (Mdl, meas));

## The type is assignable on a compact model too, under the same family rule,
## which it can only do because the fit's covariance came across with it.
%!test
%! load fisheriris
%! CMdl = compact (fitcdiscr (meas, species, 'DiscrimType', 'quadratic'));
%! CMdl.DiscrimType = 'diagQuadratic';
%! assert_equal (size (CMdl.Sigma), [1, 4, 3]);
%! assert_equal (CMdl.Gamma, 1);

%!test
%! load fisheriris
%! CMdl = compact (fitcdiscr (meas, species));
%! CMdl.Gamma = 1;
%! assert_equal (CMdl.DiscrimType, 'diagLinear');
%! assert_equal (size (CMdl.Sigma), [1, 4]);

%!test
%! load fisheriris
%! Mdl = fitcdiscr (meas, species, 'DiscrimType', 'diagQuadratic');
%! CMdl = compact (Mdl);
%! fname = tempname ();
%! savemodel (CMdl, fname);
%! C2 = loadmodel (fname);
%! delete (fname);
%! assert_equal (predict (C2, meas), predict (CMdl, meas));
%! C2.DiscrimType = 'quadratic';
%! assert_equal (size (C2.Sigma), [4, 4, 3]);

%!error<CompactClassificationDiscriminant: 'DiscrimType' can only be set to one of: quadratic, diagQuadratic, or pseudoQuadratic.> ...
%! load fisheriris; ...
%! CMdl = compact (fitcdiscr (meas, species, 'DiscrimType', 'quadratic')); ...
%! CMdl.DiscrimType = 'linear';
%!error<CompactClassificationDiscriminant: cannot set 'Gamma' to any value but 0 or 1 for a quadratic discriminant.> ...
%! load fisheriris; ...
%! CMdl = compact (fitcdiscr (meas, species, 'DiscrimType', 'quadratic')); ...
%! CMdl.Gamma = 0.5;

## edge, measured on MATLAB R2024a, whose discriminant scores this class
## reproduces exactly.
%!test
%! load fisheriris
%! CMdl = compact (fitcdiscr (meas, species));
%! assert_equal (edge (CMdl, meas, species), 0.9454289377, 1e-9);
%! assert_equal (edge (CMdl, meas, species), ...
%!               mean (margin (CMdl, meas, species)), 1e-12);

## Weights are normalized within each class to that class's prior.
%!test
%! load fisheriris
%! CMdl = compact (fitcdiscr (meas, species));
%! assert_equal (edge (CMdl, meas, species, 'Weights', (1:150)'), ...
%!               0.9438468986, 1e-9);

%!error<CompactClassificationDiscriminant.edge: too few input arguments.> ...
%! load fisheriris; edge (compact (fitcdiscr (meas, species)), meas)
