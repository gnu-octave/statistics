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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

classdef MultinomialRegression < PredictiveModel
  ## -*- texinfo -*-
  ## @deftp {statistics} MultinomialRegression
  ##
  ## Multinomial regression model class.
  ##
  ## A @code{MultinomialRegression} object holds a multinomial regression of a
  ## categorical response on one or more predictors, fitted by @code{mnrfit}.
  ## It is normally created with the @code{fitmnr} function.  The response has
  ## @math{K} categories, listed in @code{ClassNames}, and the model is one of
  ## three kinds, named by @code{ModelType}: a nominal model contrasts each
  ## category with the last, an ordinal model fits the cumulative probabilities
  ## of ordered categories, and a hierarchical model fits each category against
  ## the ones after it.
  ##
  ## The most useful properties are @code{Coefficients} (a table of estimates,
  ## standard errors, @math{t}-statistics and p-values), @code{Deviance},
  ## @code{Fitted}, @code{Residuals}, @code{LogLikelihood},
  ## @code{ModelCriterion} and @code{Rsquared}.  @code{ObservationInfo} records
  ## which rows were weighted or missing, and @code{Variables} holds the data
  ## the model was built from.  A fitted model supports @code{predict},
  ## @code{feval} and @code{random}, inference on its coefficients through
  ## @code{coefCI}, @code{coefTest} and @code{testDeviance}, and partial
  ## dependence through @code{partialDependence} and
  ## @code{plotPartialDependence}.
  ##
  ## @code{Fitted}, @code{Residuals} and @code{ObservationInfo} have one row per
  ## @emph{input} row.  A row left out of the fit because a value it needs is
  ## missing has @code{NaN} residuals and a missing label in @code{Fitted}, the
  ## empty character vector, @qcode{<undefined>}, @qcode{<missing>} or
  ## @qcode{NaN} as the type of @code{ClassNames} has it; MATLAB labels such a
  ## row with the first class.
  ##
  ## @seealso{fitmnr, mnrfit, mnrval, GeneralizedLinearModel}
  ## @end deftp

  properties (GetAccess = public, SetAccess = protected)

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} Coefficients
    ##
    ## Coefficient values
    ##
    ## A table with one row per coefficient and four columns:
    ## @itemize
    ## @item @code{Value} - estimated coefficient value
    ## @item @code{SE} - standard error of the estimate
    ## @item @code{tStat} - the estimate divided by its standard error
    ## @item @code{pValue} - p-value of that statistic
    ## @end itemize
    ##
    ## The statistic is referred to the normal distribution, or to a
    ## @math{t}-distribution on @code{DFE} degrees of freedom where the
    ## dispersion is estimated.  The row names are the coefficient names.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Coefficients = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} CoefficientCovariance
    ##
    ## Covariance matrix of the coefficient estimates
    ##
    ## A square matrix with one row and column per coefficient, in the order of
    ## @code{CoefficientNames}.  It is scaled by @code{Dispersion}, so it is the
    ## covariance under the estimated dispersion wherever one was estimated.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    CoefficientCovariance = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} CoefficientNames
    ##
    ## Names of the coefficients
    ##
    ## A string column, one name per coefficient.  A coefficient of one category
    ## joins the term to the category, @qcode{(Intercept_@var{class})} and
    ## @qcode{@var{term}_@var{class}}; a slope shared by every category is
    ## named by its term alone.  A categorical predictor contributes one term
    ## per level bar the first, spelled @qcode{@var{name}_@var{level}}, and an
    ## interaction joins its factors with a colon.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    CoefficientNames = string ([]);

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} NumCoefficients
    ##
    ## Number of coefficients
    ##
    ## A positive integer.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumCoefficients = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} ClassNames
    ##
    ## Names of the response categories
    ##
    ## A @math{K}-by-1 array in the type of the response: numeric, logical, a
    ## cell array of character vectors, a string or a categorical array.  A
    ## character matrix response gives a cell array, and a matrix of counts the
    ## numbers of its columns.  The categories of a numeric, logical or text
    ## response are sorted, and those of a categorical array are in the order of
    ## its categories.  A category without observations is no class.  The last
    ## class is the reference category of a nominal model.  The model is fitted
    ## in this order, so @code{Fitted} and @code{predict} name the class they
    ## estimate; R2024a sorts the classes of text but fits them in the order the
    ## data first shows them, so its labels name the wrong classes unless the
    ## data show them sorted.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ClassNames = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} IncludeClassInteractions
    ##
    ## Whether each category has its own slopes
    ##
    ## A logical scalar: @qcode{true} when each category has its own slopes,
    ## @qcode{false} when every category shares them.  The intercepts always
    ## differ by category.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    IncludeClassInteractions = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} Link
    ##
    ## Link function
    ##
    ## A character vector, @qcode{'logit'}, @qcode{'probit'},
    ## @qcode{'comploglog'} or @qcode{'loglog'}.  A nominal model always uses
    ## the logit link.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Link = '';

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} ModelType
    ##
    ## Type of the model
    ##
    ## A character vector, @qcode{'nominal'}, @qcode{'ordinal'} or
    ## @qcode{'hierarchical'}.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ModelType = '';

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} Deviance
    ##
    ## Deviance of the fit
    ##
    ## A nonnegative scalar, twice the difference between the log-likelihood of
    ## the saturated model and that of the fitted one.  The difference between
    ## the deviances of two nested models follows a chi-square distribution.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Deviance = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} DFE
    ##
    ## Error degrees of freedom
    ##
    ## The number of rows used in the fit times @math{K-1}, less the number of
    ## coefficients; for a hierarchical model each stage counts only the rows
    ## with observations that reach it.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    DFE = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} Dispersion
    ##
    ## Dispersion parameter
    ##
    ## A positive scalar, 1 unless @code{EstimateDispersion} is true.  The
    ## estimate is the square of the @qcode{'sfit'} that @code{mnrfit}
    ## returns, the Pearson statistic over the categories divided by
    ## @code{DFE}, so it differs from MATLAB's, as @code{mnrfit} documents.
    ## The standard errors are scaled by its square root.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Dispersion = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} EstimateDispersion
    ##
    ## Whether the dispersion was estimated
    ##
    ## A logical scalar.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    EstimateDispersion = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} Fitted
    ##
    ## Fitted response values
    ##
    ## An @math{N}-by-1 array of the most probable class of each input row, in
    ## the type of @code{ClassNames}; a row left out of the fit is a missing
    ## label.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Fitted = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} LogLikelihood
    ##
    ## Log-likelihood of the fit
    ##
    ## A scalar, @code{-@var{Deviance} / 2}: the log-likelihood measured from
    ## that of the saturated model, as MATLAB reports it.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    LogLikelihood = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} ModelCriterion
    ##
    ## Information criteria
    ##
    ## A structure with fields @code{AIC}, @code{AICc}, @code{BIC} and
    ## @code{CAIC}, computed from @code{LogLikelihood}, the number of
    ## coefficients @math{m} and the number of rows @math{n} used in the fit:
    ## @code{-2 * LogLikelihood + 2 * m}, its small-sample correction
    ## @code{+ 2 * m * (m + 1) / (n - m - 1)}, @code{-2 * LogLikelihood + m *
    ## log (n)} and @code{-2 * LogLikelihood + m * (log (n) + 1)}.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ModelCriterion = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} Residuals
    ##
    ## Residuals of the fit
    ##
    ## A table with one row per input row and three variables, @code{Raw},
    ## @code{Pearson} and @code{Deviance}, the @qcode{'resid'},
    ## @qcode{'residp'} and @qcode{'residd'} that @code{mnrfit} returns and
    ## documents.  @code{Raw} and @code{Pearson} have a column per category of a
    ## nominal model, and one fewer for an ordinal or hierarchical model;
    ## @code{Deviance} has one.  A row left out of the fit is @code{NaN}.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Residuals = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} Rsquared
    ##
    ## Pseudo R-squared values
    ##
    ## A structure with fields @code{Ordinary}, @code{1 - L / L0}, and
    ## @code{Adjusted}, @code{1 - (L - m) / L0}, where @math{L} is
    ## @code{LogLikelihood}, @math{L0} that of the model with intercepts alone
    ## and @math{m} the number of coefficients.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Rsquared = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} Formula
    ##
    ## Model formula
    ##
    ## A @code{LinearFormula} object over the variables of the model, with the
    ## link on the response, such as @qcode{logit(y) ~ 1 + x1 + x2}.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Formula = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} NumObservations
    ##
    ## Number of observations
    ##
    ## The number of input rows used in the fit, those left out because a value
    ## was missing not counted.  A row of a matrix of counts counts once,
    ## whatever number of observations it holds.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumObservations = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} NumPredictors
    ##
    ## Number of predictors
    ##
    ## The number of predictors used in the fit.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumPredictors = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} NumVariables
    ##
    ## Number of variables
    ##
    ## The number of variables in the input data, the response and any variable
    ## the model does not use included.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    NumVariables = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} ObservationInfo
    ##
    ## Information about the observations
    ##
    ## A table with one row per input row and three variables: @code{Weights},
    ## the observation weight, @code{Missing}, true where a value the fit needs
    ## is missing, and @code{Subset}, true where the row was used in the fit.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ObservationInfo = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} ObservationNames
    ##
    ## Names of the observations
    ##
    ## The row names of an input table, as a cell array of character vectors,
    ## or an empty cell array where there were none.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ObservationNames = {};

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} PredictorNames
    ##
    ## Names of the predictors
    ##
    ## A cell array of character vectors, one per predictor used in the fit.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    PredictorNames = {};

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} ResponseName
    ##
    ## Name of the response
    ##
    ## A character vector, @qcode{'y'} by default.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    ResponseName = '';

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} VariableInfo
    ##
    ## Information about the variables
    ##
    ## A table with one row per variable and four variables: @code{Class}, the
    ## class of the variable, @code{Range}, @code{[min, max]} of a continuous
    ## variable or the levels of a categorical one, @code{InModel}, true for a
    ## predictor used in the fit, and @code{IsCategorical}.  The response is not
    ## in the model and is categorical.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    VariableInfo = [];

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} VariableNames
    ##
    ## Names of the variables
    ##
    ## A cell array of character vectors naming every variable of the input
    ## data, the response included.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    VariableNames = {};

    ## -*- texinfo -*-
    ## @deftp {MultinomialRegression} {property} Variables
    ##
    ## Input data
    ##
    ## A table of the predictors and the response, with any other variable of an
    ## input table.
    ##
    ## This property is read-only.
    ##
    ## @end deftp
    Variables = [];

  endproperties

  properties (Access = private, Hidden)
    ## Coefficients in the layout mnrfit returns them
    B_ = [];
    ## mnrfit's stats structure, for the confidence bounds of predict
    stats_ = [];
    ## Terms over the encoded predictor columns, the intercept row first
    terms_ = [];
    ## Categorical predictors and their levels, to code new data alike
    catinfo_ = [];
    ## Names of the encoded predictor columns that terms_ is written over
    encnames_ = {};
    ## Which predictors are categorical
    iscat_ = [];
    ## Deviance of the model with intercepts alone
    nulldev_ = [];
  endproperties

  methods (Hidden)

    ## The model summary
    function disp (this)
      fprintf ("\n  Multinomial regression with %s responses\n\n", ...
               this.ModelType);
      if (! isempty (this.Coefficients))
        disp (this.Coefficients);
      endif
      fprintf ("\n%d observations, %d error degrees of freedom\n", ...
               this.NumObservations, this.DFE);
      fprintf ("Dispersion: %g\n", this.Dispersion);
      [p, stat] = testDeviance (this);
      if (! isnan (stat))
        if (this.EstimateDispersion)
          fprintf ("F-statistic vs. constant model: %g, p-value = %g\n", ...
                   stat, p);
        else
          fprintf (strcat ("Chi^2-statistic vs. constant model: %g,", ...
                           " p-value = %g\n"), stat, p);
        endif
      endif
    endfunction

    ## Not implemented yet
    function h = plotResiduals (this, varargin)
      error ("MultinomialRegression.plotResiduals: not implemented yet.");
    endfunction

    ## Not implemented yet
    function h = plotSlice (this, varargin)
      error ("MultinomialRegression.plotSlice: not implemented yet.");
    endfunction

  endmethods

  methods (Access = public)

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialRegression} {@var{mdl} =} MultinomialRegression (@var{X}, @var{Y})
    ## @deftypefnx {MultinomialRegression} {@var{mdl} =} MultinomialRegression (@var{Tbl}, @var{Y})
    ## @deftypefnx {MultinomialRegression} {@var{mdl} =} MultinomialRegression (@var{Tbl}, @var{ResponseVarName})
    ## @deftypefnx {MultinomialRegression} {@var{mdl} =} MultinomialRegression (@var{Tbl}, @var{Formula})
    ## @deftypefnx {MultinomialRegression} {@var{mdl} =} MultinomialRegression (@dots{}, @var{name}, @var{value})
    ##
    ## Fit a multinomial regression model.
    ##
    ## This is the constructor that @code{fitmnr} calls, and it takes exactly
    ## what @code{fitmnr} takes.  See @code{fitmnr} for the arguments and the
    ## @var{name}, @var{value} pairs.
    ##
    ## @seealso{fitmnr}
    ## @end deftypefn
    function this = MultinomialRegression (data, resp, varargin)

      if (nargin == 0)
        return;
      endif
      if (nargin < 2)
        error ("MultinomialRegression: too few input arguments.");
      endif

      ## Parse the optional paired arguments
      if (mod (numel (varargin), 2) != 0)
        error ("MultinomialRegression: optional arguments must be in pairs.");
      endif
      optNames = {'CategoricalPredictors', 'EstimateDispersion', ...
                  'IncludeClassInteractions', 'IterationLimit', 'Link', ...
                  'ModelType', 'PredictorNames', 'ResponseName', ...
                  'Tolerance', 'Weights'};
      ## IncludeClassInteractions resolves to a default per model type, and
      ## PredictorNames and ResponseName to the names the data gives
      dfValues = {[], false, [], 100, 'logit', 'nominal', [], [], 1e-6, []};
      [CatPred, EstDisp, Interact, IterLim, LinkName, ModType, PredNames, ...
       RespName, Tol, Wts, args] = parsePairedArguments (optNames, dfValues, ...
                                                         varargin);
      if (! isempty (args))
        error ("MultinomialRegression: invalid optional paired argument.");
      endif
      if (! ((ischar (ModType) && isrow (ModType))
             || (isa (ModType, 'string') && isscalar (ModType))))
        error (strcat ("MultinomialRegression: 'ModelType' must be a", ...
                       " character vector."));
      endif
      ModType = lower (char (ModType));
      if (! any (strcmp (ModType, {'nominal', 'ordinal', 'hierarchical'})))
        error ("MultinomialRegression: unknown 'ModelType' '%s'.", ModType);
      endif
      if (! ((ischar (LinkName) && isrow (LinkName))
             || (isa (LinkName, 'string') && isscalar (LinkName))))
        error ("MultinomialRegression: 'Link' must be a character vector.");
      endif
      LinkName = lower (char (LinkName));
      if (! ((islogical (EstDisp) || isnumeric (EstDisp)) && isscalar (EstDisp)
             && any (EstDisp == [0, 1])))
        error (strcat ("MultinomialRegression: 'EstimateDispersion' must", ...
                       " be true or false."));
      endif
      EstDisp = logical (EstDisp);
      if (isempty (Interact))
        Interact = ! strcmp (ModType, 'ordinal');
      elseif (! ((islogical (Interact) || isnumeric (Interact))
                 && isscalar (Interact) && any (Interact == [0, 1])))
        error (strcat ("MultinomialRegression: 'IncludeClassInteractions'", ...
                       " must be true or false."));
      endif
      Interact = logical (Interact);

      ## Resolve the predictors and the response
      [X_raw, tbl, Y, pred_names, resp_name, var_names, formula, ...
       Variables, obs_names] = MultinomialRegression.intake (data, resp, ...
                                                             PredNames, ...
                                                             RespName);
      n_total = rows (Y);
      p_raw = numel (pred_names);

      ## Which predictors are categorical: any but a numeric one in a table,
      ## and those 'CategoricalPredictors' names
      iscat = false (1, p_raw);
      if (! isempty (tbl))
        for j = 1:p_raw
          col = tbl.(pred_names{j});
          iscat(j) = ! (isnumeric (col) && isreal (col));
        endfor
      endif
      if (! isempty (CatPred))
        iscat |= MultinomialRegression.catflags (CatPred, pred_names);
      endif

      ## Rows missing a predictor, the response or a weight
      if (isempty (tbl))
        missing = any (isnan (X_raw), 2);
      else
        missing = any (ismissing (tbl(:,pred_names)), 2);
      endif
      missing |= MultinomialRegression.missingresponse (Y);
      if (isempty (Wts))
        w_full = ones (n_total, 1);
      else
        if (! (isnumeric (Wts) && isreal (Wts) && isvector (Wts)
                               && numel (Wts) == n_total))
          error (strcat ("MultinomialRegression: 'Weights' must be a", ...
                         " numeric vector with one element per", ...
                         " observation."));
        endif
        w_full = double (Wts(:));
        missing |= isnan (w_full);
      endif
      subset = ! missing;
      if (! any (subset))
        error (strcat ("MultinomialRegression: no observations remain", ...
                       " after removing missing values."));
      endif

      ## The design matrix over every input row, without the intercept, which
      ## mnrfit adds itself
      if (isempty (formula))
        [X_num, cat_levels] = raw_to_codes (data, X_raw, tbl, pred_names, ...
                                            iscat, n_total);
        [X_enc, enc_names, cat_info] = encode_categorical (X_num, iscat, ...
                                                           pred_names, ...
                                                           cat_levels, true);
        q = columns (X_enc);
        terms = [zeros(1, q); eye(q)];
        coef_terms = enc_names(:)';
        X_design = X_enc;
        X_design(missing,:) = NaN;
      else
        tbl_sub = tbl(subset,:);
        [X_sub, ~, coef_names] = parseWilkinsonFormula (formula, ...
                                   'model_matrix', tbl_sub, pred_names(iscat));
        coef_names = coef_names(:)';
        [terms, cat_info, term_cols] = terms_from_coefnames (coef_names, ...
                                         pred_names, iscat, data, tbl_sub);
        enc_names = term_cols(1:end-1);
        terms = terms(:,1:end-1);
        icpt = strcmp (coef_names, '(Intercept)');
        coef_terms = coef_names(! icpt);
        terms = [zeros(1, columns (terms)); terms(! all (terms == 0, 2),:)];
        X_design = NaN (n_total, numel (coef_terms));
        X_design(subset,:) = X_sub(:,! icpt);
      endif

      ## Fit through mnrfit
      if (Interact)
        onoff = 'on';
      else
        onoff = 'off';
      endif
      if (EstDisp)
        estdisp = 'on';
      else
        estdisp = 'off';
      endif
      ## The classes, and the response as the number of its class in each row,
      ## or as the matrix of counts it is, which is what mnrfit is given: its
      ## classes are then numbered as ClassNames lists them
      [classes, C, Yfit] = MultinomialRegression.classesof (Y(subset,:));
      K = numel (classes);
      if (K < 2)
        error (strcat ("MultinomialRegression: the response must hold at", ...
                       " least two categories."));
      endif
      [B, dev, st] = mnrfit (X_design(subset,:), Yfit, 'model', ModType, ...
                             'link', LinkName, 'interactions', onoff, ...
                             'estdisp', estdisp, 'weights', w_full(subset), ...
                             'iterationlimit', IterLim, 'tolerance', Tol);

      ## Coefficient names, in the order of B(:)
      cnames = cellstr (MultinomialRegression.classlabels (classes));

      if (strcmp (ModType, 'nominal') && ! Interact)
        names = strcat ('(Intercept_', cnames(1:K-1), ')');
      elseif (Interact)
        names = cell (numel (coef_terms) + 1, K - 1);
        for j = 1:K-1
          names(:,j) = [{['(Intercept_', cnames{j}, ')']}; ...
                        strcat(coef_terms(:), '_', cnames{j})];
        endfor
      else
        names = [strcat('(Intercept_', cnames(1:K-1), ')'); coef_terms(:)];
      endif
      names = names(:);

      this.B_ = B;
      this.stats_ = st;
      this.terms_ = terms;
      this.catinfo_ = cat_info;
      this.encnames_ = enc_names;
      this.iscat_ = iscat;
      this.ModelType = ModType;
      this.Link = LinkName;
      this.IncludeClassInteractions = Interact;
      this.EstimateDispersion = EstDisp;
      this.ClassNames = classes;
      this.Coefficients = table (B(:), st.se(:), st.t(:), st.p(:), ...
                                 'VariableNames', ...
                                 {'Value', 'SE', 'tStat', 'pValue'}, ...
                                 'RowNames', names);
      this.CoefficientNames = string (names);
      this.CoefficientCovariance = st.covb;
      this.NumCoefficients = numel (B);
      this.Deviance = dev;
      this.DFE = st.dfe;
      this.Dispersion = st.s ^ 2;

      ## The fit of the model with intercepts alone: every row has the same
      ## probabilities, the weighted shares of the categories
      wC = w_full(subset) .* C;
      share = sum (wC, 1) / sum (wC(:));
      d0 = wC .* log (wC ./ max (sum (wC, 2) .* share, realmin));
      d0(wC == 0) = 0;
      this.nulldev_ = 2 * sum (d0(:));

      n_obs = sum (subset);
      m = this.NumCoefficients;
      LL = - dev / 2;
      LL0 = - this.nulldev_ / 2;
      this.LogLikelihood = LL;
      aic = - 2 * LL + 2 * m;
      if (n_obs - m - 1 > 0)
        aicc = aic + 2 * m * (m + 1) / (n_obs - m - 1);
      else
        aicc = Inf;
      endif
      this.ModelCriterion = struct ('AIC', aic, 'AICc', aicc, ...
                                    'BIC', - 2 * LL + m * log (n_obs), ...
                                    'CAIC', - 2 * LL + m * (log (n_obs) + 1));
      this.Rsquared = struct ('Ordinary', 1 - LL / LL0, ...
                              'Adjusted', 1 - (LL - m) / LL0);

      ## Fitted labels and residuals over every input row
      P = this.probabilities (X_design(subset,:), 'category');
      [~, best] = max (P, [], 2);
      this.Fitted = MultinomialRegression.labelrows (classes, best, subset);
      R = cell (1, 3);
      fields = {'resid', 'residp', 'residd'};
      for k = 1:3
        R{k} = NaN (n_total, columns (st.(fields{k})));
        R{k}(subset,:) = st.(fields{k});
      endfor
      this.Residuals = table (R{:}, 'VariableNames', ...
                              {'Raw', 'Pearson', 'Deviance'});

      ## The model's variables
      if (isempty (formula))
        in_model = true (1, p_raw);
      else
        in_model = any (terms(2:end,:) != 0, 1);
        enc2raw = encodednames_to_row (enc_names, pred_names, cat_info);
        used = false (1, p_raw);
        used(enc2raw(enc2raw > 0 & in_model)) = true;
        in_model = used;
      endif
      this.PredictorNames = pred_names(:);
      this.ResponseName = resp_name;
      this.NumPredictors = p_raw;
      this.NumObservations = n_obs;
      this.Variables = Variables;
      this.VariableNames = var_names(:);
      this.NumVariables = numel (var_names);
      this.ObservationNames = obs_names;
      this.ObservationInfo = table (w_full, missing, subset, ...
                                    'VariableNames', ...
                                    {'Weights', 'Missing', 'Subset'});

      n_vars = numel (var_names);
      vi_class = cell (n_vars, 1);
      vi_range = cell (n_vars, 1);
      vi_inmodel = false (n_vars, 1);
      vi_iscat = false (n_vars, 1);
      for j = 1:n_vars
        col = Variables.(var_names{j});
        col = col(subset,:);
        if (isa (col, 'categorical'))
          col = removecats (col);
        endif
        k = find (strcmp (pred_names, var_names{j}), 1);
        if (strcmp (var_names{j}, resp_name))
          vi_iscat(j) = true;
          [vi_class{j}, ~] = variable_class_and_range (col);
          vi_range{j} = MultinomialRegression.classlabels (classes)(:)';
          if (isnumeric (classes) || islogical (classes))
            vi_range{j} = classes(:)';
          endif
        else
          [vi_class{j}, vi_range{j}] = variable_class_and_range (col);
          if (! isempty (k))
            vi_inmodel(j) = in_model(k);
            vi_iscat(j) = iscat(k);
          else
            vi_iscat(j) = ! (isnumeric (col) && isreal (col));
          endif
        endif
      endfor
      this.VariableInfo = table (vi_class, vi_range, vi_inmodel, vi_iscat, ...
                                 'VariableNames', {'Class', 'Range', ...
                                                   'InModel', ...
                                                   'IsCategorical'}, ...
                                 'RowNames', var_names(:));

      this.Formula = MultinomialRegression.buildformula (terms, enc_names, ...
                                                         pred_names, ...
                                                         cat_info, ...
                                                         var_names, ...
                                                         resp_name, LinkName);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialRegression} {@var{label} =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {MultinomialRegression} {[@var{label}, @var{probs}] =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {MultinomialRegression} {[@var{label}, @var{probs}, @var{lower}, @var{upper}] =} predict (@var{mdl}, @var{Xnew})
    ## @deftypefnx {MultinomialRegression} {[@dots{}] =} predict (@dots{}, @var{name}, @var{value})
    ##
    ## Predict the response of a multinomial regression model.
    ##
    ## @code{@var{label} = predict (@var{mdl}, @var{Xnew})} returns the most
    ## probable class at each row of the new predictor data @var{Xnew}, in the
    ## type of @code{ClassNames}.  @var{Xnew} is a table holding every predictor
    ## named in @code{PredictorNames}, or a numeric matrix with one column per
    ## predictor where they are all numeric.  A row holding a missing value is a
    ## missing label.
    ##
    ## @var{probs} holds the probabilities, one row per row of @var{Xnew}, and
    ## @var{lower} and @var{upper} their confidence bounds, from the delta
    ## method.  The following @var{name}, @var{value} pairs are accepted:
    ##
    ## @multitable @columnfractions 0.25 0.7
    ## @headitem Name @tab Value
    ## @item @qcode{'Alpha'} @tab The significance level of the bounds, a
    ## scalar in @math{(0, 1)}; the default 0.05 gives 95% bounds.
    ## @item @qcode{'ProbabilityType'} @tab @qcode{'category'} (default), the
    ## probability of each of the @math{K} categories; @qcode{'cumulative'},
    ## the probability of each of the first @math{K-1} categories or any before
    ## it; or @qcode{'conditional'}, the probability of each of the first
    ## @math{K-1} categories given that the response is not in an earlier one.
    ## @end multitable
    ##
    ## The bounds are at the confidence level @math{1 - @var{Alpha}}.  R2024a
    ## takes @var{Alpha} as the confidence level itself, so its default bounds
    ## are at 5% and far too narrow.
    ##
    ## @seealso{MultinomialRegression.feval, MultinomialRegression.random,
    ## mnrval}
    ## @end deftypefn
    function [label, probs, lower, upper] = predict (this, Xnew, varargin)

      if (nargin < 2)
        error ("MultinomialRegression.predict: too few input arguments.");
      endif
      [Alpha, Type, args] = parsePairedArguments ({'Alpha', ...
                                                   'ProbabilityType'}, ...
                                                  {0.05, 'category'}, varargin);
      if (! isempty (args))
        error (strcat ("MultinomialRegression.predict: invalid optional", ...
                       " paired argument."));
      endif
      if (! (isnumeric (Alpha) && isscalar (Alpha) && Alpha > 0 && Alpha < 1))
        error (strcat ("MultinomialRegression.predict: 'Alpha' must be a", ...
                       " scalar in the range (0, 1)."));
      endif
      if (! ((ischar (Type) && isrow (Type))
             || (isa (Type, 'string') && isscalar (Type))))
        error (strcat ("MultinomialRegression.predict: 'ProbabilityType'", ...
                       " must be a character vector."));
      endif
      Type = lower (char (Type));
      if (! any (strcmp (Type, {'category', 'cumulative', 'conditional'})))
        error (strcat ("MultinomialRegression.predict: unknown", ...
                       " 'ProbabilityType' '%s'."), Type);
      endif

      Xd = this.designof (Xnew, 'MultinomialRegression.predict');
      ok = ! any (isnan (Xd), 2);
      Pc = NaN (rows (Xd), numel (this.ClassNames));
      Pc(ok,:) = this.probabilities (Xd(ok,:), 'category');
      [~, best] = max (Pc(ok,:), [], 2);
      label = MultinomialRegression.labelrows (this.ClassNames, best, ok);
      if (nargout > 1)
        if (strcmp (Type, 'category'))
          probs = Pc;
        else
          probs = NaN (rows (Xd), numel (this.ClassNames) - 1);
          probs(ok,:) = this.probabilities (Xd(ok,:), Type);
        endif
      endif
      if (nargout > 2)
        [~, dlo, dhi] = this.probabilities (Xd(ok,:), Type, 1 - Alpha);
        lower = NaN (size (probs));
        upper = NaN (size (probs));
        lower(ok,:) = probs(ok,:) - dlo;
        upper(ok,:) = probs(ok,:) + dhi;
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialRegression} {@var{label} =} feval (@var{mdl}, @var{X1}, @dots{}, @var{Xp})
    ## @deftypefnx {MultinomialRegression} {@var{label} =} feval (@var{mdl}, @var{Xnew})
    ##
    ## Predict the response of a multinomial regression model.
    ##
    ## @code{@var{label} = feval (@var{mdl}, @var{X1}, @dots{}, @var{Xp})}
    ## takes one argument per predictor, each a column vector or a scalar, and a
    ## scalar stands for every row.  @code{@var{label} = feval (@var{mdl},
    ## @var{Xnew})} takes a table or a matrix, as @code{predict} does.  Either
    ## returns what @code{predict} returns first, the most probable class.
    ##
    ## @seealso{MultinomialRegression.predict}
    ## @end deftypefn
    function label = feval (this, varargin)

      if (numel (varargin) == 1
          && (istable (varargin{1})
              || columns (varargin{1}) == this.NumPredictors))
        label = predict (this, varargin{1});
        return;
      endif
      if (numel (varargin) != this.NumPredictors)
        error ("MultinomialRegression.feval: give one argument per predictor.");
      endif
      n = max (cellfun (@numel, varargin));
      X = zeros (n, this.NumPredictors);
      for j = 1:this.NumPredictors
        v = varargin{j};
        if (! (isnumeric (v) && (isscalar (v) || numel (v) == n)))
          error (strcat ("MultinomialRegression.feval: each predictor", ...
                         " must be a scalar or a vector of one common", ...
                         " length."));
        endif
        X(:,j) = v(:);
      endfor
      label = predict (this, X);

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialRegression} {@var{ysim} =} random (@var{mdl}, @var{Xnew})
    ## @deftypefnx {MultinomialRegression} {@var{ysim} =} random (@var{mdl}, @var{Xnew}, @var{NumTrials})
    ## @deftypefnx {MultinomialRegression} {@var{ysim} =} random (@dots{}, @var{type})
    ##
    ## Simulate responses of a multinomial regression model.
    ##
    ## @code{@var{ysim} = random (@var{mdl}, @var{Xnew})} draws one response at
    ## each row of @var{Xnew} from the probabilities the model gives there.
    ## @var{Xnew} is a table or a matrix, as @code{predict} takes it.
    ## @var{NumTrials}, a positive integer, draws that many at each row.
    ##
    ## @var{type} is @qcode{'samples'} (default), for an
    ## @math{r}-by-@var{NumTrials} array of class labels, in the type of
    ## @code{ClassNames}, or @qcode{'counts'}, for an @math{r}-by-@math{K}
    ## table counting how many draws fell in each class, one variable per
    ## class.
    ##
    ## @seealso{MultinomialRegression.predict}
    ## @end deftypefn
    function ysim = random (this, Xnew, varargin)

      if (nargin < 2)
        error ("MultinomialRegression.random: too few input arguments.");
      endif
      NumTrials = 1;
      type = 'samples';
      for k = 1:numel (varargin)
        a = varargin{k};
        if (isnumeric (a))
          if (! (isscalar (a) && a >= 1 && fix (a) == a))
            error (strcat ("MultinomialRegression.random: NUMTRIALS must", ...
                           " be a positive integer."));
          endif
          NumTrials = a;
        elseif ((ischar (a) && isrow (a))
                || (isa (a, 'string') && isscalar (a)))
          type = lower (char (a));
          if (! any (strcmp (type, {'samples', 'counts'})))
            error ("MultinomialRegression.random: unknown TYPE '%s'.", type);
          endif
        else
          error ("MultinomialRegression.random: invalid input argument.");
        endif
      endfor

      [~, P] = predict (this, Xnew);
      [r, K] = size (P);
      ## Each draw is the first class whose cumulative probability passes a
      ## uniform number
      cumP = cumsum (P, 2);
      idx = zeros (r, NumTrials);
      for t = 1:NumTrials
        u = rand (r, 1);
        idx(:,t) = sum (u > cumP(:,1:K-1), 2) + 1;
      endfor
      if (strcmp (type, 'samples'))
        lab = MultinomialRegression.labelrows (this.ClassNames, idx(:), ...
                                               true (r * NumTrials, 1));
        ysim = reshape (lab, r, NumTrials);
      else
        counts = zeros (r, K);
        for j = 1:K
          counts(:,j) = sum (idx == j, 2);
        endfor
        vn = cellstr (MultinomialRegression.classlabels (this.ClassNames));
        ysim = array2table (counts, 'VariableNames', vn(:)');
      endif

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialRegression} {@var{ci} =} coefCI (@var{mdl})
    ## @deftypefnx {MultinomialRegression} {@var{ci} =} coefCI (@var{mdl}, @var{alpha})
    ##
    ## Confidence intervals for the coefficients.
    ##
    ## @code{@var{ci} = coefCI (@var{mdl})} returns a @math{p}-by-2 matrix, one
    ## row per coefficient in the order of @code{CoefficientNames}, of the lower
    ## and upper limits of its 95% confidence interval,
    ## @code{Value +/- tinv (1 - @var{alpha} / 2, DFE) * SE}.  @var{alpha}, in
    ## @math{(0, 1)}, sets the level to @math{100 (1 - @var{alpha})%}.
    ##
    ## @seealso{MultinomialRegression.coefTest}
    ## @end deftypefn
    function ci = coefCI (this, alpha)

      if (nargin < 2)
        alpha = 0.05;
      endif
      if (! (isnumeric (alpha) && isscalar (alpha) && alpha > 0 && alpha < 1))
        error (strcat ("MultinomialRegression.coefCI: ALPHA must be a", ...
                       " scalar in the range (0, 1)."));
      endif
      b = this.Coefficients.Value;
      se = this.Coefficients.SE;
      t = tinv (1 - alpha / 2, this.DFE);
      ci = [b - t * se, b + t * se];

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialRegression} {@var{p} =} coefTest (@var{mdl})
    ## @deftypefnx {MultinomialRegression} {@var{p} =} coefTest (@var{mdl}, @var{H})
    ## @deftypefnx {MultinomialRegression} {@var{p} =} coefTest (@var{mdl}, @var{H}, @var{C})
    ## @deftypefnx {MultinomialRegression} {[@var{p}, @var{F}, @var{r}] =} coefTest (@dots{})
    ##
    ## Linear hypothesis test on the coefficients.
    ##
    ## @code{@var{p} = coefTest (@var{mdl}, @var{H}, @var{C})} tests the
    ## hypothesis that @code{@var{H} * @var{b} = @var{C}}, where @var{b} is the
    ## column of coefficients in the order of @code{CoefficientNames}.  @var{H}
    ## is a full-rank @math{r}-by-@math{p} matrix and @var{C} a vector of
    ## @math{r} elements, zero by default; without @var{H} the hypothesis is
    ## that every coefficient is zero.  @var{F} is the Wald statistic divided by
    ## @var{r}, and @var{p} its p-value from the @math{F} distribution on
    ## @var{r} and @code{DFE} degrees of freedom.
    ##
    ## @seealso{MultinomialRegression.coefCI,
    ## MultinomialRegression.testDeviance}
    ## @end deftypefn
    function [p, F, r] = coefTest (this, H, C)

      b = this.Coefficients.Value;
      nc = numel (b);
      if (nargin < 2)
        H = eye (nc);
      endif
      if (! ((isnumeric (H) || islogical (H)) && ismatrix (H)
                                              && columns (H) == nc))
        error (strcat ("MultinomialRegression.coefTest: H must have one", ...
                       " column per coefficient."));
      endif
      H = double (H);
      r = rows (H);
      if (nargin < 3)
        C = zeros (r, 1);
      endif
      if (! (isnumeric (C) && isvector (C) && numel (C) == r))
        error (strcat ("MultinomialRegression.coefTest: C must have one", ...
                       " element per row of H."));
      endif
      d = H * b - C(:);
      F = (d' * ((H * this.CoefficientCovariance * H') \ d)) / r;
      p = fcdf (F, r, this.DFE, 'upper');

    endfunction

    ## -*- texinfo -*-
    ## @deftypefn  {MultinomialRegression} {@var{p} =} testDeviance (@var{mdl})
    ## @deftypefnx {MultinomialRegression} {[@var{p}, @var{stat}] =} testDeviance (@var{mdl})
    ##
    ## Deviance test against the model with intercepts alone.
    ##
    ## @code{@var{p} = testDeviance (@var{mdl})} tests whether the model fits
    ## better than one with an intercept for each category and no predictors,
    ## comparing their deviances.  Where the dispersion is fixed, @var{stat} is
    ## the drop in deviance and @var{p} its chi-square p-value on as many
    ## degrees of freedom as the model has coefficients beyond the intercepts;
    ## where it is estimated, @var{stat} is that drop per degree of freedom
    ## divided by @code{Dispersion}, and @var{p} its @math{F} p-value on those
    ## and @code{DFE} degrees of freedom.
    ##
    ## @seealso{MultinomialRegression.coefTest}
    ## @end deftypefn
    function [p, stat] = testDeviance (this)

      df1 = this.NumCoefficients - (numel (this.ClassNames) - 1);
      drop = this.nulldev_ - this.Deviance;
      if (df1 < 1)
        p = NaN;
        stat = NaN;
      elseif (this.EstimateDispersion)
        stat = (drop / df1) / this.Dispersion;
        p = fcdf (stat, df1, this.DFE, 'upper');
      else
        stat = drop;
        p = chi2cdf (stat, df1, 'upper');
      endif

    endfunction

  endmethods

  methods (Access = private)

    ## Probabilities of TYPE at the design rows XD, and with a confidence level
    ## their delta-method bounds
    function [P, dlo, dhi] = probabilities (this, Xd, type, level)
      if (this.IncludeClassInteractions)
        onoff = 'on';
      else
        onoff = 'off';
      endif
      B = this.B_;
      if (strcmp (this.ModelType, 'nominal') && ! this.IncludeClassInteractions)
        ## The intercepts alone
        B = B(:)';
        Xd = zeros (rows (Xd), 0);
        onoff = 'on';
      endif
      args = {'model', this.ModelType, 'link', this.Link, ...
              'interactions', onoff, 'type', type};
      if (nargout > 1)
        [P, dlo, dhi] = mnrval (B, Xd, this.stats_, args{:}, ...
                                'confidence', level);
      else
        P = mnrval (B, Xd, args{:});
      endif
    endfunction

    ## The design matrix, without the intercept, of new predictor data
    function Xd = designof (this, Xnew, caller)
      p = this.NumPredictors;
      if (istable (Xnew))
        X_raw = zeros (height (Xnew), p);
        for j = 1:p
          name = this.PredictorNames{j};
          if (! any (strcmp (Xnew.Properties.VariableNames, name)))
            error ("%s: the table holds no predictor '%s'.", caller, name);
          endif
          col = Xnew.(name);
          k = [];
          if (! isempty (this.catinfo_.names))
            k = find (strcmp (this.catinfo_.names, name));
          endif
          if (! isempty (k))
            levels = this.catinfo_.levels{k};
            [tf, code] = ismember (MultinomialRegression.levelstrings (col), ...
                                   cellstr (levels));
            code(! tf) = NaN;
            X_raw(:,j) = code;
          else
            X_raw(:,j) = double (col);
          endif
        endfor
      else
        if (! (isnumeric (Xnew) && isreal (Xnew) && ismatrix (Xnew)
                                && columns (Xnew) == p))
          error (strcat ("%s: XNEW must be a table or a numeric matrix", ...
                         " with %d columns."), caller, p);
        endif
        if (any (this.iscat_))
          ## A categorical predictor given as numbers is coded by its levels
          X_raw = double (Xnew);
          for j = find (this.iscat_)
            k = find (strcmp (this.catinfo_.names, this.PredictorNames{j}));
            levels = str2double (this.catinfo_.levels{k});
            [tf, code] = ismember (X_raw(:,j), levels);
            code(! tf) = NaN;
            X_raw(:,j) = code;
          endfor
        else
          X_raw = double (Xnew);
        endif
      endif
      X_enc = reencode_predictors (X_raw, this.PredictorNames(:)', ...
                                   this.catinfo_, this.encnames_);
      ## A level the model never saw codes as NaN, and so does its row
      bad = any (isnan (X_raw), 2);
      Xd = build_design ([this.terms_, zeros(rows (this.terms_), 1)], X_enc);
      Xd = Xd(:,2:end);
      Xd(bad,:) = NaN;
    endfunction

  endmethods

  methods (Static, Access = private)

    ## Resolve the predictors and the response from the forms fitmnr takes
    function [X_raw, tbl, Y, pred_names, resp_name, var_names, formula, ...
              Variables, obs_names] = intake (data, resp, PredNames, RespName)
      formula = '';
      obs_names = {};
      if (! istable (data))
        if (! (isnumeric (data) && isreal (data) && ismatrix (data)))
          error (strcat ("MultinomialRegression: X must be a real", ...
                         " numeric matrix or a table."));
        endif
        X_raw = double (data);
        tbl = [];
        p = columns (X_raw);
        Y = MultinomialRegression.responsecolumn (resp);
        if (rows (Y) != rows (X_raw))
          error (strcat ("MultinomialRegression: X and Y must have the", ...
                         " same number of rows."));
        endif
        if (isempty (PredNames))
          pred_names = arrayfun (@(k) sprintf ('x%d', k), 1:p, ...
                                 'UniformOutput', false);
        else
          pred_names = cellstr (PredNames)(:)';
          if (numel (pred_names) != p)
            error (strcat ("MultinomialRegression: 'PredictorNames' must", ...
                           " have %d names."), p);
          endif
        endif
        if (isempty (RespName))
          resp_name = 'y';
        else
          resp_name = char (RespName);
        endif
        cols = num2cell (X_raw, 1);
        Variables = table (cols{:}, Y, ...
                           'VariableNames', [pred_names, {resp_name}]);
        var_names = [pred_names, {resp_name}];
        return;
      endif

      tbl = data;
      X_raw = [];
      vn = tbl.Properties.VariableNames;
      if (! isempty (tbl.Properties.RowNames))
        obs_names = tbl.Properties.RowNames(:);
      endif
      istext = (ischar (resp) && isrow (resp)) ...
               || (isa (resp, 'string') && isscalar (resp));
      if (istext && any (char (resp) == '~'))
        ## A formula names the response and the predictors
        formula = char (resp);
        parts = strsplit (formula, '~');
        resp_name = strtrim (parts{1});
        if (! any (strcmp (vn, resp_name)))
          error ("MultinomialRegression: the table holds no variable '%s'.", ...
                 resp_name);
        endif
        schema = parseWilkinsonFormula (formula, 'matrix');
        used = regexprep (schema.VariableNames, '\^\d+$', '');
        used = unique (used, 'stable');
        pred_names = vn(ismember (vn, used) & ! strcmp (vn, resp_name));
        Y = MultinomialRegression.responsecolumn (tbl.(resp_name));
        Variables = tbl;
      elseif (istext)
        resp_name = char (resp);
        if (! any (strcmp (vn, resp_name)))
          error ("MultinomialRegression: the table holds no variable '%s'.", ...
                 resp_name);
        endif
        Y = MultinomialRegression.responsecolumn (tbl.(resp_name));
        if (isempty (PredNames))
          pred_names = vn(! strcmp (vn, resp_name));
        else
          pred_names = cellstr (PredNames)(:)';
        endif
        Variables = tbl;
      else
        ## The response is given beside the table, all of whose variables are
        ## predictors unless 'PredictorNames' names some
        Y = MultinomialRegression.responsecolumn (resp);
        if (rows (Y) != height (tbl))
          error (strcat ("MultinomialRegression: TBL and Y must have the", ...
                         " same number of rows."));
        endif
        if (isempty (RespName))
          resp_name = 'y';
        else
          resp_name = char (RespName);
        endif
        if (isempty (PredNames))
          pred_names = vn;
        else
          pred_names = cellstr (PredNames)(:)';
        endif
        Variables = tbl;
        Variables.(resp_name) = Y;
      endif
      missing = pred_names(! ismember (pred_names, vn));
      if (! isempty (missing))
        error ("MultinomialRegression: the table holds no variable '%s'.", ...
               missing{1});
      endif
      ## Code a string or logical predictor by its text, as a cell array is
      for j = 1:numel (pred_names)
        col = tbl.(pred_names{j});
        if (isa (col, 'string'))
          tbl.(pred_names{j}) = cellstr (col);
        elseif (islogical (col))
          c = repmat ({'false'}, size (col));
          c(col) = {'true'};
          tbl.(pred_names{j}) = c;
        endif
      endfor
      var_names = Variables.Properties.VariableNames;
    endfunction

    ## The response as one row per observation: a character matrix becomes a
    ## cell array, and a row vector of labels a column
    function Y = responsecolumn (Y)
      if (ischar (Y))
        Y = cellstr (Y);
      elseif ((iscell (Y) || isa (Y, 'string') || isa (Y, 'categorical')
               || isnumeric (Y) || islogical (Y)) && isrow (Y))
        Y = Y(:);
      endif
      if (! (isnumeric (Y) || islogical (Y) || iscellstr (Y)
             || isa (Y, 'string') || isa (Y, 'categorical')))
        error (strcat ("MultinomialRegression: Y must be a numeric,", ...
                       " logical, string or categorical array, a cell", ...
                       " array of character vectors, or a character", ...
                       " matrix."));
      endif
    endfunction

    ## Which rows of the response are missing
    function tf = missingresponse (Y)
      if (isnumeric (Y))
        tf = any (isnan (Y), 2);
      elseif (islogical (Y))
        tf = false (rows (Y), 1);
      else
        tf = ismissing (Y);
        if (iscellstr (Y))
          tf = cellfun (@isempty, Y);
        endif
      endif
    endfunction

    ## The classes of the response, the counts of each class in each row, and
    ## the response as mnrfit is to be given it: the number of its class in
    ## each row, or the matrix of counts.  Text is sorted, as MATLAB sorts it,
    ## and a categorical array keeps the order of the categories it uses.
    function [classes, C, Yfit] = classesof (Y)
      if ((isnumeric (Y) || islogical (Y)) && columns (Y) > 1)
        keep = any (Y > 0, 1);
        classes = find (keep)(:);
        C = double (Y(:,keep));
        Yfit = Y;
        return;
      endif
      if (isa (Y, 'categorical'))
        [~, code] = ismember (cellstr (Y), categories (Y));
        [present, first, num] = unique (code);
        classes = Y(first);
      elseif (iscellstr (Y) || isa (Y, 'string'))
        [~, first, num] = unique (cellstr (Y));
        classes = Y(first);
      else
        [~, first, num] = unique (Y);
        classes = Y(first);
      endif
      classes = classes(:);
      num = num(:);
      C = double (num == (1:numel (classes)));
      Yfit = num;
    endfunction

    ## The labels of the classes as a string column
    function s = classlabels (classes)
      if (islogical (classes))
        s = repmat (string ('false'), numel (classes), 1);
        s(classes(:)) = 'true';
      elseif (isnumeric (classes))
        s = string (arrayfun (@(v) num2str (v), double (classes(:)), ...
                              'UniformOutput', false));
      else
        s = string (cellstr (classes));
        s = s(:);
      endif
    endfunction

    ## The labels at the rows TF, the class numbered IDX each, and a missing
    ## label at every other row
    function lab = labelrows (classes, idx, tf)
      n = numel (tf);
      if (iscellstr (classes))
        lab = repmat ({''}, n, 1);
      elseif (isa (classes, 'string'))
        lab = repmat (string (missing), n, 1);
      elseif (isa (classes, 'categorical'))
        lab = classes(ones (n, 1));
        lab(:) = missing;
      elseif (islogical (classes))
        lab = false (n, 1);
      else
        lab = NaN (n, 1);
      endif
      lab(tf) = classes(idx);
    endfunction

    ## The text of each element of a categorical predictor
    function s = levelstrings (col)
      if (isnumeric (col) || islogical (col))
        s = arrayfun (@(v) num2str (v), double (col(:)), ...
                      'UniformOutput', false);
      else
        s = cellstr (col);
        s = s(:);
      endif
    endfunction

    ## Logical flags of the predictors 'CategoricalPredictors' names
    function tf = catflags (spec, pred_names)
      p = numel (pred_names);
      tf = false (1, p);
      if ((ischar (spec) && strcmpi (spec, 'all'))
          || (isa (spec, 'string') && isscalar (spec) && strcmpi (spec, 'all')))
        tf(:) = true;
      elseif (islogical (spec))
        if (numel (spec) != p)
          error (strcat ("MultinomialRegression: 'CategoricalPredictors'", ...
                         " must have %d elements."), p);
        endif
        tf = spec(:)';
      elseif (isnumeric (spec))
        if (! all (spec >= 1 & spec <= p & fix (spec) == spec))
          error (strcat ("MultinomialRegression: 'CategoricalPredictors'", ...
                         " indices must be between 1 and %d."), p);
        endif
        tf(spec) = true;
      elseif (iscellstr (spec) || isa (spec, 'string') || ischar (spec))
        names = cellstr (spec);
        [ok, k] = ismember (names, pred_names);
        if (! all (ok))
          error (strcat ("MultinomialRegression: 'CategoricalPredictors'", ...
                         " names no predictor '%s'."), names{find (! ok, 1)});
        endif
        tf(k) = true;
      else
        error ("MultinomialRegression: invalid 'CategoricalPredictors'.");
      endif
    endfunction

    ## The model formula over the variables, with the link on the response
    function F = buildformula (terms, enc_names, pred_names, cat_info, ...
                               var_names, resp_name, link)
      n_vars = numel (var_names);
      var_idx = zeros (1, numel (pred_names));
      for j = 1:numel (pred_names)
        var_idx(j) = find (strcmp (var_names, pred_names{j}), 1);
      endfor
      [enc2raw, col_pow] = encodednames_to_row (enc_names, pred_names, ...
                                                cat_info);
      terms_var = variable_level_terms (terms, enc2raw, var_idx, n_vars, ...
                                        col_pow);
      F = LinearFormula (terms_var, var_names(:)', ...
                         'ResponseName', resp_name, 'Link', link);
    endfunction

  endmethods

endclassdef

## Test the fit against MATLAB R2024a
%!shared x, z, y, Xm, m
%! x = [0.2, 1.1, 1.9, 3.2, 3.8, 5.1, 5.9, 7.2, 8.1, 8.8, 1.5, 4.4, 6.6, ...
%!      2.7, 9.3]';
%! z = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9]' / 3;
%! y = [1, 1, 2, 1, 2, 2, 3, 2, 3, 3, 2, 3, 2, 1, 3]';
%! Xm = [x, z];
%! m = fitmnr (Xm, y);
%!test
%! assert_equal (m.Coefficients.Value, [5.9358487983056; -1.30171796766448; ...
%!               -0.474611351178345; 2.68369517887323; -0.725396483366165; ...
%!               0.882324211575015], -1e-10);
%!test
%! assert_equal (m.Coefficients.SE, [3.20733455995645; 0.681543706775439; ...
%!               1.36563950586161; 2.73281100769504; 0.404575904177872; ...
%!               0.92708159102731], -1e-8);
%!test
%! assert_equal (m.Coefficients.pValue(1:2), [0.0642111655098015; ...
%!               0.0561389778598589], -1e-8);
%!test
%! assert_equal (m.CoefficientNames, string ({'(Intercept_1)'; 'x1_1'; ...
%!               'x2_1'; '(Intercept_2)'; 'x1_2'; 'x2_2'}));
%!test
%! assert_equal (m.ClassNames, [1; 2; 3]);
%!test
%! assert_equal ([m.NumCoefficients, m.DFE, m.NumObservations], [6, 24, 15]);
%!test
%! assert_equal (m.LogLikelihood, -9.30901781932902, 1e-10);
%!test
%! assert_equal ([m.Rsquared.Ordinary, m.Rsquared.Adjusted], ...
%!               [0.428116752664665, 0.0595172343688309], 1e-10);
%!test
%! mc = m.ModelCriterion;
%! assert_equal ([mc.AIC, mc.AICc, mc.BIC, mc.CAIC], [30.618035638658, ...
%!               41.118035638658, 34.8663368452713, 40.8663368452713], 1e-10);
%!test
%! assert_equal (char (m.Formula), 'logit(y) ~ 1 + x1 + x2');
%!test
%! assert_equal (m.VariableInfo.Range{3}, [1, 2, 3]);
%!test
%! assert_equal (m.VariableInfo.InModel, [true; true; false]);
%!test
%! assert_equal (m.VariableInfo.IsCategorical, [false; false; true]);
%!test
%! assert_equal (m.Fitted(1:5), [1; 1; 1; 1; 2]);
%!test
%! assert_equal (m.Residuals.Properties.VariableNames, ...
%!               {'Raw', 'Pearson', 'Deviance'});
%!test
%! [~, ~, s] = mnrfit (Xm, y);
%! assert_equal (m.Residuals.Pearson, s.residp);
%!test
%! assert_equal (m.ObservationInfo.Properties.VariableNames, ...
%!               {'Weights', 'Missing', 'Subset'});
%!test
%! assert_equal (m.Dispersion, 1);
%!test
%! [~, ~, s] = mnrfit (Xm, y, 'estdisp', 'on');
%! md = fitmnr (Xm, y, 'EstimateDispersion', true);
%! assert_equal (md.Dispersion, s.s ^ 2, 1e-12);
%!test
%! md = fitmnr (Xm, y, 'EstimateDispersion', true);
%! se = m.Coefficients.SE * sqrt (md.Dispersion);
%! assert_equal (md.Coefficients.SE, se, -1e-12);

## Test the model types against MATLAB R2024a
%!test
%! mo = fitmnr (Xm, y, 'ModelType', 'ordinal');
%! assert_equal (mo.Coefficients.Value, [1.73183634703477; ...
%!               5.07302821549392; -0.795693713914391; ...
%!               0.00358605506018692], -1e-6);
%!test
%! mo = fitmnr (Xm, y, 'ModelType', 'ordinal');
%! assert_equal (mo.CoefficientNames, string ({'(Intercept_1)'; ...
%!               '(Intercept_2)'; 'x1'; 'x2'}));
%!test
%! mo = fitmnr (Xm, y, 'ModelType', 'ordinal');
%! assert_equal ([mo.Rsquared.Ordinary, mo.Rsquared.Adjusted], ...
%!               [0.354457649951417, 0.108724637754195], 1e-8);
%!test
%! mh = fitmnr (Xm, y, 'ModelType', 'hierarchical', 'Link', 'probit');
%! assert_equal (mh.Coefficients.Value, [1.94768532702581; ...
%!               -0.474517402244592; -0.637169549255952; ...
%!               1.43714192930196; -0.419757944547173; ...
%!               0.547210016922966], -1e-10);
%!test
%! mh = fitmnr (Xm, y, 'ModelType', 'hierarchical', 'Link', 'probit');
%! assert_equal (char (mh.Formula), 'probit(y) ~ 1 + x1 + x2');
%!test
%! mo = fitmnr (Xm, y, 'ModelType', 'ordinal');
%! assert_equal (mo.IncludeClassInteractions, false);
%!test
%! mo = fitmnr (Xm, y, 'ModelType', 'ordinal', ...
%!              'IncludeClassInteractions', true);
%! assert_equal (mo.CoefficientNames(1:3), string ({'(Intercept_1)'; ...
%!               'x1_1'; 'x2_1'}));

## Test the classes of each type of response
%!test
%! names = {'mid', 'lo', 'hi'};
%! mt = fitmnr (x, names(y)');
%! assert_equal (mt.ClassNames, {'hi'; 'lo'; 'mid'});
%!test
%! names = {'mid', 'lo', 'hi'};
%! mt = fitmnr (x, names(y)');
%! assert_equal (mt.Fitted(1:3), {'mid'; 'mid'; 'mid'});
%!test
%! names = {'mid', 'lo', 'hi'};
%! mt = fitmnr (x, names(y)');
%! assert_equal (mt.CoefficientNames, string ({'(Intercept_hi)'; 'x1_hi'; ...
%!               '(Intercept_lo)'; 'x1_lo'}));
%!test
%! names = {'mid', 'lo', 'hi'};
%! ms = fitmnr (x, string (names(y)'));
%! assert_equal (ms.ClassNames, string ({'hi'; 'lo'; 'mid'}));
%!test
%! g = [2, 5, 9];
%! mg = fitmnr (x, g(y)');
%! assert_equal (mg.ClassNames, [2; 5; 9]);
%!test
%! names = {'mid', 'lo', 'hi'};
%! yk = categorical (names(y)', {'lo', 'mid', 'unused', 'hi'});
%! mk = fitmnr (x, yk);
%! assert_equal (cellstr (mk.ClassNames), {'lo'; 'mid'; 'hi'});
%!test
%! ml = fitmnr (x, y == 1);
%! assert_equal (ml.CoefficientNames, string ({'(Intercept_false)'; ...
%!               'x1_false'}));
%!test
%! Yc = [3, 1, 0; 2, 2, 1; 1, 2, 2; 0, 2, 3; 1, 0, 4; 0, 1, 1];
%! mc = fitmnr ((1:6)', Yc);
%! assert_equal (mc.Fitted, [1; 1; 2; 3; 3; 3]);
%!test
%! Yc = [3, 1, 0; 2, 2, 1; 1, 2, 2; 0, 2, 3; 1, 0, 4; 0, 1, 1];
%! mc = fitmnr ((1:6)', [Yc(:,1), zeros(6,1), Yc(:,2:3)]);
%! assert_equal (mc.ClassNames, [1; 3; 4]);

## Test tables and formulas against MATLAB R2024a
%!test
%! tbl = table (x, z, y);
%! mtb = fitmnr (tbl, 'y');
%! [~, pr] = predict (mtb, tbl(1:3,:));
%! [~, pm] = predict (m, Xm(1:3,:));
%! assert_equal (pr, pm, 1e-12);
%!test
%! tbl = table (x, z, y);
%! tbl.g = categorical (mod ((1:15)', 3));
%! S = warning ('off', 'all');
%! mtg = fitmnr (tbl, 'y ~ x + g');
%! warning (S);
%! assert_equal (mtg.CoefficientNames, string ({'(Intercept_1)'; 'x_1'; ...
%!               'g_1_1'; 'g_2_1'; '(Intercept_2)'; 'x_2'; 'g_1_2'; ...
%!               'g_2_2'}));
%!test
%! tbl = table (x, z, y);
%! tbl.g = categorical (mod ((1:15)', 3));
%! S = warning ('off', 'all');
%! mtg = fitmnr (tbl, 'y ~ x + g');
%! warning (S);
%! [~, pr] = predict (mtg, tbl(1:3,:));
%! assert_equal (pr(:,2), [0.00316347560363747; 0.193496952442625; ...
%!               0.913115541220257], 1e-7);
%!test
%! tbl = table (x, z, y);
%! mt2 = fitmnr (tbl, 'y ~ x*z');
%! [~, pr] = predict (mt2, tbl(1:4,:));
%! B = mnrfit ([x, z, x .* z], y);
%! assert_equal (pr, mnrval (B, [x(1:4), z(1:4), x(1:4) .* z(1:4)]), 1e-12);
%!test
%! tbl = table (x, z, y);
%! tbl.g = categorical (mod ((1:15)', 3));
%! mt2 = fitmnr (tbl, 'y ~ x*z');
%! assert_equal (mt2.VariableInfo.InModel, [true; true; false; false]);

## Test missing values and weights
%!test
%! Xn = Xm;
%! Xn(2,1) = NaN;
%! mm = fitmnr (Xn, y);
%! assert_equal (mm.Fitted(1:3), [1; NaN; 1]);
%!test
%! Xn = Xm;
%! Xn(2,1) = NaN;
%! mm = fitmnr (Xn, y);
%! assert_equal (mm.ObservationInfo.Subset(1:3), [true; false; true]);
%!test
%! Xn = Xm;
%! Xn(2,1) = NaN;
%! mm = fitmnr (Xn, y);
%! assert_equal (mm.Residuals.Deviance(2), NaN);
%!test
%! w = 1 + mod ((1:15)', 3);
%! mw = fitmnr (Xm, y, 'Weights', w);
%! [~, dev] = mnrfit (Xm, y, 'weights', w);
%! assert_equal (mw.Deviance, dev);

## Test the methods against MATLAB R2024a
%!test
%! [lab, pr] = predict (m, Xm(1:3,:));
%! assert_equal (pr, [0.851674439820704, 0.143631440997612, ...
%!               0.00469411918168434; 0.886831611467944, ...
%!               0.101673529742501, 0.0114948587895541; ...
%!               0.566508298983404, 0.400052723700518, ...
%!               0.0334389773160776], 1e-10);
%!test
%! assert_equal (predict (m, Xm(1:3,:)), [1; 1; 1]);
%!test
%! [~, pr, lo, up] = predict (m, Xm(1:3,:));
%! [B, ~, s] = mnrfit (Xm, y);
%! [~, dlo, dhi] = mnrval (B, Xm(1:3,:), s);
%! assert_equal ([lo, up], [pr - dlo, pr + dhi], 1e-12);
%!test
%! [~, pr, lo] = predict (m, Xm(1:3,:), 'Alpha', 0.1);
%! [B, ~, s] = mnrfit (Xm, y);
%! [~, dlo] = mnrval (B, Xm(1:3,:), s, 'confidence', 0.9);
%! assert_equal (lo, pr - dlo, 1e-12);
%!test
%! [~, pc] = predict (m, Xm(1:3,:), 'ProbabilityType', 'cumulative');
%! assert_equal (pc(:,2), [0.995305880818316; 0.988505141210446; ...
%!               0.966561022683922], 1e-10);
%!test
%! [~, pk] = predict (m, Xm(1:3,:), 'ProbabilityType', 'conditional');
%! assert_equal (pk(:,2), [0.968352594279704; 0.898426946440983; ...
%!               0.922861320671979], 1e-10);
%!test
%! assert_equal (feval (m, x(1:3), z(1:3)), predict (m, Xm(1:3,:)));
%!test
%! assert_equal (feval (m, x(1:3), 1), predict (m, [x(1:3), ones(3,1)]));
%!test
%! [p, F, r] = coefTest (m);
%! assert_equal ([p, F, r], [0.407376657888557, 1.07023588159414, 6], 1e-8);
%!test
%! [p, F, r] = coefTest (m, [0, 1, 0, 0, 0, 0]);
%! assert_equal ([p, F, r], [0.0681608433167721, 3.64792898304325, 1], 1e-8);
%!test
%! [p, F, r] = coefTest (m, [0, 1, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0], [0.5; 0]);
%! assert_equal ([p, F, r], [0.0377095939218408, 3.76922166716901, 2], 1e-8);
%!test
%! [p, stat] = testDeviance (m);
%! assert_equal ([p, stat], [0.00749682238009594, 13.9376227503715], 1e-10);
%!test
%! mo = fitmnr (Xm, y, 'ModelType', 'ordinal');
%! [p, stat] = testDeviance (mo);
%! assert_equal ([p, stat], [0.00312037814994931, 11.5396021651965], 1e-8);
%!test
%! md = fitmnr (Xm, y, 'EstimateDispersion', true);
%! [p, stat] = testDeviance (md);
%! assert_equal (stat, 13.9376227503715 / 4 / md.Dispersion, 1e-10);
%!test
%! assert_equal (coefCI (m), [-0.683764386648367, 12.5554619832596; ...
%!               -2.70835504376494, 0.104919108435981; ...
%!               -3.29315276302852, 2.34393006067183; ...
%!               -2.95654952910981, 8.32393988685627; ...
%!               -1.56040011006823, 0.109607143335903; ...
%!               -1.03107815065807, 2.7957265738081], 1e-8);
%!test
%! ci = coefCI (m, 0.1);
%! t = tinv (0.95, 24);
%! assert_equal (ci(:,2), m.Coefficients.Value + t * m.Coefficients.SE, 1e-12);
%!test
%! rs = random (m, Xm(1:3,:), 5);
%! assert_equal (size (rs), [3, 5]);
%!test
%! rc = random (m, Xm(1:3,:), 4, 'counts');
%! assert_equal (sum (rc{:,:}, 2), [4; 4; 4]);
%!test
%! names = {'mid', 'lo', 'hi'};
%! mt = fitmnr (x, names(y)');
%! rc = random (mt, x(1:3), 2, 'counts');
%! assert_equal (rc.Properties.VariableNames, {'hi', 'lo', 'mid'});
%!test
%! pd = partialDependence (m, 1, 1);
%! assert_equal (size (pd), [1, 100]);

## Test input validation
%!error<MultinomialRegression: too few input arguments.> MultinomialRegression (1)
%!error<MultinomialRegression: invalid optional paired argument.> ...
%! fitmnr (Xm, y, 'Foo', 1)
%!error<MultinomialRegression: unknown 'ModelType' 'linear'.> ...
%! fitmnr (Xm, y, 'ModelType', 'linear')
%!error<MultinomialRegression: 'ModelType' must be a character vector.> ...
%! fitmnr (Xm, y, 'ModelType', 1)
%!error<MultinomialRegression: 'Link' must be a character vector.> ...
%! fitmnr (Xm, y, 'Link', 1)
%!error<MultinomialRegression: 'EstimateDispersion' must be true or false.> ...
%! fitmnr (Xm, y, 'EstimateDispersion', 2)
%!error<MultinomialRegression: 'IncludeClassInteractions' must be true or false.> ...
%! fitmnr (Xm, y, 'IncludeClassInteractions', 'on')
%!error<MultinomialRegression: 'Weights' must be a numeric vector with one element per observation.> ...
%! fitmnr (Xm, y, 'Weights', ones (3, 1))
%!error<MultinomialRegression: X must be a real numeric matrix or a table.> ...
%! fitmnr ({1}, y)
%!error<MultinomialRegression: X and Y must have the same number of rows.> ...
%! fitmnr (Xm, y(1:5))
%!error<MultinomialRegression: 'PredictorNames' must have 2 names.> ...
%! fitmnr (Xm, y, 'PredictorNames', {'a'})
%!error<MultinomialRegression: the table holds no variable 'w'.> ...
%! fitmnr (table (x, y), 'w')
%!error<MultinomialRegression: the table holds no variable 'w'.> ...
%! fitmnr (table (x, y), 'w ~ x')
%!error<MultinomialRegression: TBL and Y must have the same number of rows.> ...
%! fitmnr (table (x), y(1:5))
%!error<MultinomialRegression: the response must hold at least two categories.> ...
%! fitmnr (x, ones (15, 1))
%!error<MultinomialRegression: no observations remain after removing missing values.> ...
%! fitmnr (NaN (15, 1), y)
%!error<MultinomialRegression: 'CategoricalPredictors' indices must be between 1 and 2.> ...
%! fitmnr (Xm, y, 'CategoricalPredictors', 3)
%!error<MultinomialRegression: 'CategoricalPredictors' names no predictor 'w'.> ...
%! fitmnr (Xm, y, 'CategoricalPredictors', {'w'})
%!error<MultinomialRegression: 'CategoricalPredictors' must have 2 elements.> ...
%! fitmnr (Xm, y, 'CategoricalPredictors', true)
%!error<MultinomialRegression.predict: too few input arguments.> predict (m)
%!error<MultinomialRegression.predict: invalid optional paired argument.> ...
%! predict (m, Xm, 'Foo', 1)
%!error<MultinomialRegression.predict: 'Alpha' must be a scalar in the range \(0, 1\)\.> ...
%! predict (m, Xm, 'Alpha', 2)
%!error<MultinomialRegression.predict: unknown 'ProbabilityType' 'other'.> ...
%! predict (m, Xm, 'ProbabilityType', 'other')
%!error<MultinomialRegression.predict: 'ProbabilityType' must be a character vector.> ...
%! predict (m, Xm, 'ProbabilityType', 1)
%!error<MultinomialRegression.predict: XNEW must be a table or a numeric matrix with 2 columns.> ...
%! predict (m, x)
%!error<MultinomialRegression.predict: the table holds no predictor 'x2'.> ...
%! predict (m, table (x, 'VariableNames', {'x1'}))
%!error<MultinomialRegression.feval: give one argument per predictor.> ...
%! feval (m, 1, 2, 3)
%!error<MultinomialRegression.feval: each predictor must be a scalar or a vector of one common length.> ...
%! feval (m, [1; 2], [1; 2; 3])
%!error<MultinomialRegression.random: too few input arguments.> random (m)
%!error<MultinomialRegression.random: NUMTRIALS must be a positive integer.> ...
%! random (m, Xm, 0)
%!error<MultinomialRegression.random: unknown TYPE 'draws'.> ...
%! random (m, Xm, 'draws')
%!error<MultinomialRegression.random: invalid input argument.> ...
%! random (m, Xm, {1})
%!error<MultinomialRegression.coefCI: ALPHA must be a scalar in the range \(0, 1\)\.> ...
%! coefCI (m, 1)
%!error<MultinomialRegression.coefTest: H must have one column per coefficient.> ...
%! coefTest (m, [1, 0])
%!error<MultinomialRegression.coefTest: C must have one element per row of H.> ...
%! coefTest (m, [1, 0, 0, 0, 0, 0], [1; 2])
